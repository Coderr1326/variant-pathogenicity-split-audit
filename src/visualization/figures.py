#!/usr/bin/env python
import json
from pathlib import Path
import matplotlib.pyplot as plt,seaborn as sns
ROOT=Path(__file__).resolve().parents[2]; out=ROOT/'results/figures'; out.mkdir(parents=True,exist_ok=True)
for model in ['dnabert2','nt']:
 f=ROOT/'results/metrics'/f'{model}_100bp.json'
 if not f.exists(): continue
 d=json.loads(f.read_text()); h=max(d['history'],key=lambda x:x.get('weighted_f1',-1)); cm=[[h['tn'],h['fp']],[h['fn'],h['tp']]]
 sns.heatmap(cm,annot=True,fmt='d',cmap='Blues',xticklabels=['Benign','Pathogenic'],yticklabels=['Benign','Pathogenic']); plt.xlabel('Predicted'); plt.ylabel('Actual'); plt.title(model+' at 100 bp'); plt.tight_layout(); plt.savefig(out/f'{model}_100bp_confusion_matrix.png',dpi=180); plt.close()
