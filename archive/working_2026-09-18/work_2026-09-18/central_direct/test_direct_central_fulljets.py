from pathlib import Path
import pandas as pd, numpy as np, importlib.util
B=Path('/mnt/data')
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(sp); sp.loader.exec_module(zk)
sp=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); rr=importlib.util.module_from_spec(sp); sp.loader.exec_module(rr)
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
src=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
jet=pd.read_csv(B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
ref=pd.read_csv(B/'ssz_p5_CENTRAL_SELECTED_41of41_CORRECTED_A5_V12_2026-09-16.csv')
assert len(raw)==len(src)==len(jet)==len(ref)
inp=pd.DataFrame({
'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,
'f2X':src.f2X,'f2F':np.ones(len(src)),'f2Y':0.,
'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,
'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':jet.f3XX_selected,'tf3':src.tilde_f3,
'f4':src.f4,'f4X':src.N4,'f4XX':jet.f4XX_recovered,'f4XXX':jet.f4XXX_selected,'tf4':0.,
'f3phi':jet.f3phi,'f3phiX':jet.f3phiX,'f4phi':jet.f4phi,'f4phiX':jet.f4phiX,
})
d=zk.emit(inp,selected_v5=ref.v5.to_numpy(float),selected_c3=jet.c3_selected.to_numpy(float),selected_e3=jet.e3_selected.to_numpy(float))
d.to_csv(B/'direct_central_fulljets.csv',index=False)
# compare key slots to corrected ref
print('slot comparison fulljets vs corrected ref')
for k in [*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],*[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],*[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]:
 rel=np.abs(d[k]-ref[k])/np.maximum(1,np.abs(ref[k]));
 if np.nanmedian(rel)>1e-5 or np.nanquantile(rel,.95)>1e-4:
  print(k,'med',np.nanmedian(rel),'p95',np.nanquantile(rel,.95),'max',np.nanmax(rel),'abs',np.nanmax(np.abs(d[k]-ref[k])))
# K2 direct-vs-Eq433 derivative residual using emitted slots and alpha7 from emitted
hj=zk.hj; r=d.x.to_numpy(float); der=lambda y,o=1:hj.local_poly_deriv(r,np.asarray(y,float),order=o,window=9,degree=8)
f=d.f.to_numpy(float);h=d.h.to_numpy(float);Ap=d.A0prime.to_numpy(float);fp=der(f)
for L in [6,42,1000]:
 a4=d.a4.to_numpy(float); v6=d.v6.to_numpy(float); v10=d.v10.to_numpy(float); alpha7=d.alpha7.to_numpy(float)
 K2=(4*r*r*a4/f)/(fp/f+L/(r*h)-2/r+Ap*v6/(2*a4)); K2p_direct=der(K2)
 A=r*r*f*h; B=L*der(f**3/(r*r*h))+2*f*(r*fp-2*f)**2/r**3-4*f*f*(2*f*f*alpha7+Ap*Ap*v10)/(r*a4)
 K2p_eq=der(A)/A*K2-B*K2*K2/(8*r*f*f*a4)
 rel=np.abs(K2p_direct-K2p_eq)/np.maximum(1,np.abs(K2p_direct))
 print('K2p',L,'median',np.median(rel),'p95',np.quantile(rel,.95),'max',np.max(rel))
# reducer production region
p=d[(d.u>=.61)&(d.u<.71)].reset_index(drop=True)
for L in [6,12,20,42,110,420,1000]:
 a=rr.canonical_audit(p,float(L)); K=.5*(a['K']+np.swapaxes(a['K'],1,2)); G=.5*(a['G']+np.swapaxes(a['G'],1,2)); ev=np.linalg.eigvalsh(K); mine=ev[:,0]
 rad=[]
 for i in range(len(p)):
  w,V=np.linalg.eigh(K[i])
  if w[0]>0:
   W=V@np.diag(1/np.sqrt(w))@V.T; rad.append(np.linalg.eigvalsh(W@G[i]@W)[0])
 print('RED',L,'Kmin',mine.min(),'neg',int((mine<=0).sum()),'u',p.u.iloc[mine.argmin()],'radmin',min(rad) if rad else np.nan)
