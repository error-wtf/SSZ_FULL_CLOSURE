import pandas as pd,numpy as np,importlib.util,itertools
from pathlib import Path
P=Path('/mnt/data'); sp=importlib.util.spec_from_file_location('r',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); r=importlib.util.module_from_spec(sp); sp.loader.exec_module(r)
d=pd.read_csv(P/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv').tail(80).copy().reset_index(drop=True); sel=pd.read_csv(P/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
for name,col in [('c3','c3_selected'),('e3','e3_selected')]: d[name]=np.interp(d.u,sel.u,sel[col])
rows=[]
for q in [0.03,.05,.1,.2,.3,.5,.75,1,1.5,2,3,5,10,20]:
 for A in [0,100,200,300,400,500,600,800,1000,1500,2000]:
  x=d.copy(); x.v1=x.v1*q; x.v4=x.v4+A
  vals=[]
  for L in [6,42,1000]:
   z=r.canonical_audit(x,L); K=np.asarray(z['K']); vals.append(float(np.linalg.eigvalsh((K+K.transpose(0,2,1))/2)[-1,0]))
  rows.append((q,A,*vals,min(vals)))
res=pd.DataFrame(rows,columns=['q','A','e6','e42','e1000','score']).sort_values('score',ascending=False)
print(res.head(25).to_string(index=False)); res.to_csv(P/'ssz_p5_CENTRAL_ENDPOINT_V1V4_SCAN_2026-09-16.csv',index=False)
