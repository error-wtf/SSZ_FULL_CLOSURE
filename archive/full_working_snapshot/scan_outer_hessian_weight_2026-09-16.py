from pathlib import Path
import importlib.util, numpy as np, pandas as pd
B=Path('/mnt/data')
SLOTS=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],*[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],*[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]
# modules
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py'); zk=importlib.util.module_from_spec(sp); sp.loader.exec_module(zk)
sr=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); rr=importlib.util.module_from_spec(sr); sr.loader.exec_module(rr)
res=pd.read_csv(B/'ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv').sort_values('u').reset_index(drop=True)
src=pd.read_csv(B/'ssz_p5_outer_svt_repaired_integrable_transition_profile_2026-09-12.csv').sort_values('u').reset_index(drop=True)
mh=pd.read_csv(B/'ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv').sort_values('u').reset_index(drop=True)
u=res.u.to_numpy(float); S=res.S_SVT.to_numpy(float); r=res.x.to_numpy(float); h=res.h.to_numpy(float); X=res.X.to_numpy(float); ph=-np.sqrt(np.maximum(0,-2*X/h))
def ip(c): return np.interp(u,src.u.to_numpy(float),src[c].to_numpy(float))
basecols=dict(u=u,x=r,phi=res.phi,f=res.f,h=res.h,phiprime=ph,A0prime=res.A0prime,X=X,f2X=res.f2X,f2F=res.f2F,f2Y=0.,f3=res.f3,f3X=res.f3X,f3XX=0.,tf3=0.,f4=res.f4,f4X=res.N4,f4XX=0.,f4XXX=0.,tf4=0.)
# shared maps independent of p
d0=pd.DataFrame(basecols)
for z in ['f2XX','f2XF','f2XY','f2FF','f2FY','f2YY']: d0[z]=0.
ehdf=d0.copy()
for c in ['f2X','f2F','f2Y','f2XX','f2XF','f2XY','f2FF','f2FY','f2YY','f3','f3X','f3XX','tf3','f4','f4X','f4XX','f4XXX','tf4']: ehdf[c]=0.
eh=zk.emit(ehdf,selected_v5=0.,selected_c3=0.,selected_e3=0.)
emdf=ehdf.copy(); emdf['f2F']=1.; em=zk.emit(emdf,selected_v5=0.,selected_c3=0.,selected_e3=0.)
rows=[]
for p in [1.0,.75,.5,.35,.25,.15,.10,.05,.02]:
 W=np.where(S>0,S**p,0.)
 d=pd.DataFrame(basecols)
 for c,z in [('HXX','f2XX'),('HXF','f2XF'),('HXY','f2XY'),('HFF','f2FF'),('HFY','f2FY'),('HYY','f2YY')]: d[z]=W*ip(c)
 raw=zk.emit(d,selected_v5=0.,selected_c3=0.,selected_e3=0.)
 tot=pd.DataFrame({'u':u,'x':r,'phi':res.phi,'f':res.f,'h':res.h,'phiprime':ph,'A0prime':res.A0prime,'T_H':res.T_H,'S_SVT':res.S_SVT})
 for k in SLOTS:
   shared=eh[k].to_numpy(float)+S*(em[k].to_numpy(float)-eh[k].to_numpy(float))
   delta=raw[k].to_numpy(float)-shared
   tot[k]=mh[k].to_numpy(float)+delta
 for L in [42,1000]:
   a=rr.canonical_audit(tot,L); K=.5*(a['K']+np.swapaxes(a['K'],1,2)); ev=np.linalg.eigvalsh(K)[:,0]; neg=ev<0
   rows.append(dict(p=p,L=L,minK=float(ev.min()),neg_rows=int(neg.sum()),neg_u_min=float(tot.u[neg].min()) if neg.any() else np.nan,neg_u_max=float(tot.u[neg].max()) if neg.any() else np.nan,max_K_asym=float(np.max(np.abs(a['K']-np.swapaxes(a['K'],1,2))))))
 print('p',p,[r for r in rows if r['p']==p])
pd.DataFrame(rows).to_csv(B/'ssz_p5_OUTER_HESSIAN_WEIGHT_SCAN_2026-09-16.csv',index=False)
