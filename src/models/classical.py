import torch
from torch import nn

class CNN(nn.Module):
    def __init__(self, vocab=5, emb=32, hidden=64, dropout=.1):
        super().__init__(); self.embed=nn.Embedding(vocab,emb); self.net=nn.Sequential(nn.Conv1d(emb,hidden,7,padding=3),nn.ReLU(),nn.MaxPool1d(2),nn.Conv1d(hidden,hidden,5,padding=2),nn.ReLU(),nn.AdaptiveMaxPool1d(1)); self.head=nn.Sequential(nn.Flatten(),nn.Dropout(dropout),nn.Linear(hidden,2))
    def forward(self,x): return self.head(self.net(self.embed(x).transpose(1,2)))

class BiLSTM(nn.Module):
    def __init__(self,vocab=5,emb=32,hidden=64,dropout=.1):
        super().__init__(); self.embed=nn.Embedding(vocab,emb); self.rnn=nn.LSTM(emb,hidden,batch_first=True,bidirectional=True); self.head=nn.Sequential(nn.Dropout(dropout),nn.Linear(hidden*2,2))
    def forward(self,x): _,(h,_)=self.rnn(self.embed(x)); return self.head(torch.cat([h[-2],h[-1]],1))

class CNNBiLSTM(nn.Module):
    def __init__(self,vocab=5,emb=32,hidden=64,dropout=.1):
        super().__init__(); self.embed=nn.Embedding(vocab,emb); self.conv=nn.Sequential(nn.Conv1d(emb,hidden,7,padding=3),nn.ReLU()); self.rnn=nn.LSTM(hidden,hidden,batch_first=True,bidirectional=True); self.head=nn.Sequential(nn.Dropout(dropout),nn.Linear(hidden*2,2))
    def forward(self,x): x=self.conv(self.embed(x).transpose(1,2)).transpose(1,2); _,(h,_)=self.rnn(x); return self.head(torch.cat([h[-2],h[-1]],1))

class EnsembleCNNRNN(nn.Module):
    """Best-effort interpretation of the paper's ambiguous LSTM/BiLSTM/GRU ensemble."""
    def __init__(self,vocab=5,emb=32,hidden=48,dropout=.1):
        super().__init__(); self.embed=nn.Embedding(vocab,emb); self.conv=nn.Sequential(nn.Conv1d(emb,hidden,7,padding=3),nn.ReLU()); self.lstm=nn.LSTM(hidden,hidden,batch_first=True); self.bilstm=nn.LSTM(hidden,hidden,batch_first=True,bidirectional=True); self.gru=nn.GRU(hidden,hidden,batch_first=True); self.head=nn.Sequential(nn.Linear(hidden*4,hidden),nn.ReLU(),nn.Dropout(dropout),nn.Linear(hidden,2))
    def forward(self,x): x=self.conv(self.embed(x).transpose(1,2)).transpose(1,2); _,(a,_)=self.lstm(x); _,(b,_)=self.bilstm(x); _,c=self.gru(x); return self.head(torch.cat([a[-1],b[-2],b[-1],c[-1]],1))
