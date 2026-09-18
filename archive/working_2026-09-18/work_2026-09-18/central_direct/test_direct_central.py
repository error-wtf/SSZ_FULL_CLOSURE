from pathlib import Path
import pandas as pd, numpy as np, importlib.util
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(sp); sp.loader.exec_module(zk)
sp=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); rr=importlib.util.module_from_spec(sp); sp.loader.exec_module(rr)
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
src=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
assert len(raw)==len(src) and np.max(abs(raw.u-src.u))<1e-12
inp=pd.DataFrame({
'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,
'f2X':src.f2X,'f2F':np.ones(len(src)),'f2Y':0.,
'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,
'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':0.,'tf3':src.tilde_f3,
'f4':src.f4,'f4X':src.N4,'f4XX':0.,'f4XXX':0.,'tf4':0.,
})
d=zk.emit(inp,selected_v5=0.,selected_c3=0.,selected_e3=0.)
# production central only, but retain guard rows for derivative? scan trimmed to exact region excludes inner u>=.71
p=d[(d.u>=.61)&(d.u<.71)].reset_index(drop=True)
print('rows',len(p),'v12res',np.max(abs(p.v12 + p.v6/(2*p.h))))
for L in [6,12,20,42,110,420,1000]:
 a=rr.canonical_audit(p,float(L)); K=.5*(a['K']+np.swapaxes(a['K'],1,2)); G=.5*(a['G']+np.swapaxes(a['G'],1,2));
 ev=np.linalg.eigvalsh(K); mine=ev[:,0]; rad=[]
 for i in range(len(p)):
  w,V=np.linalg.eigh(K[i])
  if w[0]>0:
   W=V@np.diag(1/np.sqrt(w))@V.T; rad.append(np.linalg.eigvalsh(W@G[i]@W)[0])
 print(L,'Kmin',mine.min(),'neg',int((mine<=0).sum()),'u',p.u.iloc[mine.argmin()],'radmin',min(rad) if rad else np.nan,'diag',a['diagnostics'])
# compare direct emitter selected key slots to raw/prestaged
ref=pd.read_csv(B/'ssz_p5_CENTRAL_SELECTED_41of41_CORRECTED_A5_V12_2026-09-16.csv')
for k in ['v1','v4','c2','a2','a5','d3','e4','v6','v13']:
 rel=np.abs(d[k]-ref[k])/np.maximum(1,np.abs(ref[k])); print(k,'median',np.median(rel),'p95',np.quantile(rel,.95),'max',np.max(rel))
d.to_csv('/mnt/data/direct_central_emitted.csv',index=False)
