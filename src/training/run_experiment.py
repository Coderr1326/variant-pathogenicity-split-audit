#!/usr/bin/env python
from __future__ import annotations
import argparse,json,random,time,copy,sys
from pathlib import Path
import numpy as np,pandas as pd,torch
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,roc_auc_score,precision_recall_fscore_support,confusion_matrix
from torch import nn
from torch.utils.data import Dataset,DataLoader
from transformers import AutoTokenizer,AutoModel
from src.training.transformer_utils import load_transformer,compute_max_length
from src.models.classical import CNN,BiLSTM,CNNBiLSTM,EnsembleCNNRNN
from src.data.splits import load_manifest,apply_split

ROOT=Path(__file__).resolve().parents[2]
def seed(s=42): random.seed(s); np.random.seed(s); torch.manual_seed(s); torch.cuda.manual_seed_all(s)
class PairSet(Dataset):
 def __init__(self,df,kind,tokenizer=None,maxlen=200):
  self.df=df.reset_index(drop=True); self.kind=kind; self.tok=tokenizer; self.maxlen=maxlen
 def __len__(self): return len(self.df)
 def __getitem__(self,i):
  r=self.df.iloc[i]; text=r.ref_sequence+r.mut_sequence
  if self.kind=='classical':
   ids=[{'A':0,'C':1,'G':2,'T':3}.get(x,4) for x in text][:self.maxlen]
   ids += [4]*(self.maxlen-len(ids))
   return torch.tensor(ids,dtype=torch.long),torch.tensor(int(r.label))
  z=self.tok(text,padding='max_length',truncation=True,max_length=self.maxlen,return_tensors='pt'); return {k:v.squeeze(0) for k,v in z.items()},torch.tensor(int(r.label))
def metrics(y,p,prob):
 cm=confusion_matrix(y,p,labels=[0,1]); pr,re,f1,_=precision_recall_fscore_support(y,p,labels=[0,1],zero_division=0); out={'accuracy':accuracy_score(y,p),'roc_auc':roc_auc_score(y,prob) if len(set(y))==2 else float('nan'),'weighted_precision':precision_recall_fscore_support(y,p,average='weighted',zero_division=0)[0],'weighted_recall':precision_recall_fscore_support(y,p,average='weighted',zero_division=0)[1],'weighted_f1':precision_recall_fscore_support(y,p,average='weighted',zero_division=0)[2]}
 for i,n in enumerate(['benign','pathogenic']): out.update({f'{n}_precision':pr[i],f'{n}_recall':re[i],f'{n}_f1':f1[i]})
 out['tn'],out['fp'],out['fn'],out['tp']=map(int,cm.ravel()); assert sum(cm.ravel())==len(y); return out
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--data',required=True); ap.add_argument('--model',choices=['cnn','bilstm','cnn_bilstm','ensemble','dnabert2','nt'],required=True); ap.add_argument('--epochs',type=int,default=10); ap.add_argument('--batch-size',type=int,default=16); ap.add_argument('--tiny',action='store_true'); ap.add_argument('--fraction',type=float,default=1.0); ap.add_argument('--window',type=int,default=100); ap.add_argument('--split-manifest',default=None,help='path to a frozen split manifest (e.g. data/splits/split_manifest_b.json); default omits this and uses the original random 80/20 split (SPLIT-A / S1), unchanged'); args=ap.parse_args(); seed()
 df=pd.read_parquet(args.data)
 # captured before --tiny/--fraction subsampling so a quick gate run still sees the
 # true longest ref+mut strings in the full dataset, not just whatever survived sampling
 full_pairs=df.ref_sequence+df.mut_sequence; full_longest=full_pairs.loc[full_pairs.str.len().nlargest(20).index].tolist()
 if args.tiny: df=df.groupby('label',group_keys=False).head(128)
 if args.fraction < 1.0: df=df.groupby('label',group_keys=False).sample(frac=args.fraction,random_state=42)
 split_tag=''
 if args.split_manifest:
  manifest=load_manifest(args.split_manifest); split_tag='_split'+manifest['split_id'].upper()
  tr,va=apply_split(df,manifest)
 else:
  tr,va=train_test_split(df,test_size=.2,random_state=42,stratify=df.label)
 device=torch.device('cuda' if torch.cuda.is_available() else 'cpu'); use_amp=device.type=='cuda'
 tok=None; kind='classical'; maxlen=2*args.window
 if args.model in ['dnabert2','nt']:
  name='zhihan1996/DNABERT-2-117M' if args.model=='dnabert2' else 'InstaDeepAI/nucleotide-transformer-v2-100m-multi-species'
  tok,model=load_transformer(name)
  kind='transformer'
  if args.model=='nt':
   maxlen=compute_max_length(tok,full_longest)
 else: model={'cnn':CNN,'bilstm':BiLSTM,'cnn_bilstm':CNNBiLSTM,'ensemble':EnsembleCNNRNN}[args.model](vocab=5)
 print('max_length used:',maxlen,file=sys.stderr)
 model.to(device)
 ds1=PairSet(tr,kind,tok,maxlen); ds2=PairSet(va,kind,tok,maxlen); dl1=DataLoader(ds1,args.batch_size,shuffle=True); dl2=DataLoader(ds2,args.batch_size)
 counts=np.bincount(tr.label,minlength=2); weights=torch.tensor(len(tr)/(2*np.maximum(counts,1)),dtype=torch.float,device=device); lossfn=nn.CrossEntropyLoss(weight=weights); opt=torch.optim.AdamW(model.parameters(),lr=2e-5,weight_decay=.01); scaler=torch.amp.GradScaler('cuda',enabled=use_amp); best=-1; beststate=None; history=[]
 for ep in range(args.epochs):
  model.train(); t=time.time()
  for x,y in dl1:
   y=y.to(device); opt.zero_grad(set_to_none=True)
   with torch.autocast(device_type='cuda',dtype=torch.float16,enabled=use_amp): logits=model(**{k:v.to(device) for k,v in x.items()}) if kind=='transformer' else model(x.to(device)); loss=lossfn(logits,y)
   scaler.scale(loss).backward(); scaler.unscale_(opt); torch.nn.utils.clip_grad_norm_(model.parameters(),1.0); scaler.step(opt); scaler.update()
  model.eval(); ys=[];ps=[];probs=[]
  with torch.no_grad():
   for x,y in dl2:
    logits=model(**{k:v.to(device) for k,v in x.items()}) if kind=='transformer' else model(x.to(device)); q=torch.softmax(logits,1)[:,1].cpu().numpy(); probs.extend(q); ps.extend((q>=.5).astype(int)); ys.extend(y.numpy())
  m=metrics(ys,ps,probs); history.append({'epoch':ep+1,**m,'seconds':time.time()-t}); print(json.dumps(history[-1]))
  if m['weighted_f1']>best: best=m['weighted_f1']; beststate=copy.deepcopy(model.state_dict())
 out=ROOT/'results/metrics'; out.mkdir(parents=True,exist_ok=True); (out/f'{args.model}_{args.window}bp{split_tag}.json').write_text(json.dumps({'model':args.model,'window':args.window,'n':len(df),'history':history},indent=2));
 if beststate: torch.save(beststate,ROOT/f'results/{args.model}_{args.window}bp{split_tag}.pt')
if __name__=='__main__': main()
