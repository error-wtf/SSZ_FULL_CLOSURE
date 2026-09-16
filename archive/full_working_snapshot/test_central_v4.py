import importlib.util,pandas as pd,numpy as np
from pathlib import Path
P=Path('/mnt/data'); sp=importlib.util.spec_from_file_location('r',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); r=importlib.util.module_from_spec(sp); sp.loader.exec_module(r)
d=pd.read_csv(P/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
# selected c3,e3 merge? baseline file c3/e3 maybe zero? Need use selected rep where possible
sel=pd.read_csv(P/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
for c in ['c3_selected','e3_selected']:
 if c in sel: d[c.split('_')[0]]=np.interp(d.u,sel.u,sel[c])
for A in [0,1,3,5,10,15,-1,-3]:
 x=d.copy(); x.v4=x.v4+A
 for L in [6,42,1000]:
  a=r.canonical_audit(x,L); K=np.asarray(a['K']); ev=np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[:,0]
  print(A,L,ev.min(),x.u.iloc[ev.argmin()],(ev<=0).sum())
 print()
