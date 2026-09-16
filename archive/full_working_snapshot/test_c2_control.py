import sys, importlib.util, pandas as pd, numpy as np
from pathlib import Path
P=Path('/mnt/data')
spec=importlib.util.spec_from_file_location('red',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')
red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
d=pd.read_csv(P/'ssz_p5_OUTER_FINAL_STRONGH_HYBRID_41of41_2026-09-16.csv')
c=pd.read_csv(P/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
# central endpoint values
ce=c.iloc[0]
print('endpoint deltas', {k:float(ce[k]-d.iloc[-1][k]) for k in ['v1','v4','c2']})
S=d.S_SVT.to_numpy(float)
# c-infinity-ish use stored S which is already flat endpoints
D=float(ce.c2-d.iloc[-1].c2)
for p in [0.25,0.5,1,1.5,2,3,4,6,8]:
  x=d.copy(); x['c2']=x.c2 + D*S**p
  vals=[]
  for L in [6,12,20,42,110,420,1000]:
    a=red.canonical_audit(x,L)
    K=np.asarray(a['K']); ev=np.linalg.eigvalsh((K+np.swapaxes(K,1,2))/2)[:,0]
    vals.append((L,float(ev.min()),int((ev<=0).sum()),float(x.u.iloc[ev.argmin()])))
  print('p',p,vals)
