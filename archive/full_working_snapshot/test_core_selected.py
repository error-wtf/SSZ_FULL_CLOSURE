from pathlib import Path
import importlib.util,numpy as np,pandas as pd
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('red',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(sp); sp.loader.exec_module(red)
d=pd.read_csv(B/'ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv').sort_values('x').reset_index(drop=True)
# epsY
ph=d.phi_r.to_numpy(float); kap=d.h.to_numpy(float)*ph**2; za=1-.02*kap; d['v1']*=za; d['v10']*=za; d['v7']=d['v2']**2/(4*d['v1']); d['v12']=-d['v6']/(2*d['h'])
for L in [6,42,1000]:
 a=red.canonical_audit(d,float(L)); K=a['K'];G=a['G']; ke=[];ge=[]
 for i in range(len(d)):
  ks=(K[i]+K[i].T)/2; gs=(G[i]+G[i].T)/2; w,V=np.linalg.eigh(ks); ke.append(w[0]);
  if w[0]>0:
   W=V@np.diag(1/np.sqrt(w))@V.T; ge.append(np.linalg.eigvalsh(W@gs@W)[0])
  else: ge.append(np.nan)
 print('L',L,'Kmin',np.nanmin(ke),'u',d.u.iloc[int(np.nanargmin(ke))],'fracpos',np.mean(np.array(ke)>0),'Gmin',np.nanmin(ge),'diag',a['diagnostics'])
