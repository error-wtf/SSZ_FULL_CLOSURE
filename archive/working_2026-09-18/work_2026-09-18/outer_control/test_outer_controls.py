import importlib.util, numpy as np, pandas as pd
from pathlib import Path
R=Path('/mnt/data/ssz_work/repo/repo')
spec=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
d0=pd.read_csv(R/'archive/full_working_snapshot/ssz_p5_OUTER_HYBRID_SELECTED_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
# locks
for c in ['v5','c3','e3']: d0[c]=0.0
d0['v12']=-d0.v6/(2*d0.h)
S=d0.S_SVT.to_numpy(float)
weak_c2=float(d0.iloc[-1].c2) if d0.iloc[-1].u<d0.iloc[0].u else float(d0.iloc[0].c2)
# sorted x => u descending? inspect
print('u ends',d0.u.iloc[0],d0.u.iloc[-1])
# map weak endpoint is min u
weak_idx=int(np.argmin(d0.u)); central_idx=int(np.argmax(d0.u))
weak_c2=float(d0.c2.iloc[weak_idx]); central_c2=54.77125224555108
weak_v1=float(d0.v1.iloc[weak_idx]); central_v1=20.09857175106137

def audit(label,d):
 print('\n',label)
 for L in [6,12,20,42,110,420,1000]:
  a=red.canonical_audit(d,L); K=.5*(a['K']+a['K'].transpose(0,2,1)); G=.5*(a['G']+a['G'].transpose(0,2,1))
  mine=np.linalg.eigvalsh(K)[:,0]; rad=[]
  for k,g in zip(K,G):
   w,V=np.linalg.eigh(k)
   if w[0]>0:
    W=V@np.diag(1/np.sqrt(w))@V.T; rad.append(np.linalg.eigvalsh(W@g@W)[0])
   else: rad.append(np.nan)
  print(L,'K',np.nanmin(mine),'neg',np.sum(mine<=0),'rad',np.nanmin(rad))

audit('base',d0)
for variant in range(1,7):
 d=d0.copy()
 if variant==1:
  d['c2']=(1-S)*weak_c2+S*central_c2
 if variant==2:
  d['c2']=(1-S)*weak_c2+S*central_c2; d['v4']=0
 if variant==3:
  d['c2']=(1-S)*weak_c2+S*central_c2; d['v4']=0; d['v1']=(1-S)*weak_v1+S*central_v1
 if variant==4:
  d['c2']=central_c2
 if variant==5:
  d['c2']=np.maximum(d.c2, (1-S)*weak_c2+S*central_c2)
 if variant==6:
  d['c2']=(1-S)*weak_c2+S*central_c2; d['v4']*=0.1
 audit('var'+str(variant),d)
