import importlib.util,numpy as np,pandas as pd
from pathlib import Path
from itertools import product
R=Path('/mnt/data/ssz_work/repo/repo'); spec=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
d0=pd.read_csv(R/'archive/full_working_snapshot/ssz_p5_OUTER_HYBRID_SELECTED_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
for c in ['v5','c3','e3']: d0[c]=0.; d0['v12']=-d0.v6/(2*d0.h)
S=d0.S_SVT.to_numpy(float)
Ls=[6,42,1000]
def score(d):
 vals=[]
 for L in Ls:
  a=red.canonical_audit(d,L); K=.5*(a['K']+a['K'].transpose(0,2,1)); vals.append(np.linalg.eigvalsh(K)[:,0].min())
 return vals
rows=[]
for a,b in product([.1,.25,.5,1,2,4,8,16,32,64],[-10,-5,-2,-1,-.5,0,.5,1,2,5,10]):
 d=d0.copy(); d['v1']=np.maximum(1e-6,a*d0.v1); d['v4']=b*d0.v4
 try: vals=score(d); rows.append([a,b,*vals,min(vals)])
 except Exception: pass
rows=sorted(rows,key=lambda x:x[-1],reverse=True)
for r in rows[:30]: print(r)
