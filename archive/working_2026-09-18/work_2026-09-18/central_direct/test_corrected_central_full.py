import pandas as pd,numpy as np,importlib.util
B='/mnt/data/'
sp=importlib.util.spec_from_file_location('rr',B+'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py');rr=importlib.util.module_from_spec(sp);sp.loader.exec_module(rr)
d=pd.read_csv(B+'ssz_p5_CENTRAL_SELECTED_41of41_CORRECTED_A5_V12_2026-09-16.csv')
p=d[(d.u>=.61)&(d.u<.71)].reset_index(drop=True)
for L in [6,12,20,42,110,420,1000]:
 a=rr.canonical_audit(p,L); K=.5*(a['K']+np.swapaxes(a['K'],1,2)); G=.5*(a['G']+np.swapaxes(a['G'],1,2)); ev=np.linalg.eigvalsh(K); mine=ev[:,0]; rad=[]
 for i in range(len(p)):
  w,V=np.linalg.eigh(K[i])
  if w[0]>0:
   W=V@np.diag(1/np.sqrt(w))@V.T; rad.append(np.linalg.eigvalsh(W@G[i]@W)[0])
 print(L,'Kmin',mine.min(),'neg',int((mine<=0).sum()),'u',p.u.iloc[mine.argmin()],'radmin',min(rad) if rad else np.nan,'diag',a['diagnostics'])
