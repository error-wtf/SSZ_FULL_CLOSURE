import importlib.util, pandas as pd, numpy as np
from pathlib import Path
P=Path('/mnt/data'); spec=importlib.util.spec_from_file_location('red',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
d=pd.read_csv(P/'ssz_p5_OUTER_FINAL_STRONGH_HYBRID_41of41_2026-09-16.csv'); S=d.S_SVT.to_numpy(float)

def test(x,L=6):
 a=red.canonical_audit(x,L); K=np.asarray(a['K']); ev=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
 return float(ev.min()),int((ev<=0).sum()),float(x.u.iloc[ev.argmin()])
print('baseline',test(d,6),test(d,42))
for A in [-100,-30,-10,-3,-1,-.3,.3,1,3,10,30,100]:
 x=d.copy(); x.v4=x.v4 + A*S
 try: print('v4',A,test(x,6),test(x,42))
 except Exception as e: print('ERR v4',A,e)
for A in [-100,-30,-10,-3,-1,-.3,.3,1,3,10,30,100]:
 x=d.copy(); x.v1=x.v1 + A*S
 if np.min(x.v1)<=0: continue
 try: print('v1',A,test(x,6),test(x,42))
 except Exception as e: print('ERR v1',A,e)
