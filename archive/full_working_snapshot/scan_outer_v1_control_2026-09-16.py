from pathlib import Path
import importlib.util,numpy as np,pandas as pd
B=Path('/mnt/data');SLOTS=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],*[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],*[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]
sp=importlib.util.spec_from_file_location('zk',B/'ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py');zk=importlib.util.module_from_spec(sp);sp.loader.exec_module(zk)
sr=importlib.util.spec_from_file_location('rr',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py');rr=importlib.util.module_from_spec(sr);sr.loader.exec_module(rr)
res=pd.read_csv(B/'ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv').sort_values('u').reset_index(drop=True); src=pd.read_csv(B/'ssz_p5_outer_svt_repaired_integrable_transition_profile_2026-09-12.csv').sort_values('u').reset_index(drop=True);mh=pd.read_csv(B/'ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv').sort_values('u').reset_index(drop=True)
u=res.u.values;S=res.S_SVT.values;T=res.T_H.values;r=res.x.values;h=res.h.values;f=res.f.values;X=res.X.values;Ap=res.A0prime.values;ph=-np.sqrt(-2*X/h)
def ip(c):return np.interp(u,src.u,src[c])
base=dict(u=u,x=r,phi=res.phi,f=f,h=h,phiprime=ph,A0prime=Ap,X=X,f2X=res.f2X,f2F=res.f2F,f2Y=0.,f3=res.f3,f3X=res.f3X,f3XX=0.,tf3=0.,f4=res.f4,f4X=res.N4,f4XX=0.,f4XXX=0.,tf4=0.)
# baseline shared maps
ehdf=pd.DataFrame(base)
for c in ['f2X','f2F','f2Y','f3','f3X','f3XX','tf3','f4','f4X','f4XX','f4XXX','tf4']:ehdf[c]=0.
for z in ['f2XX','f2XF','f2XY','f2FF','f2FY','f2YY']:ehdf[z]=0.
eh=zk.emit(ehdf,selected_v5=0,selected_c3=0,selected_e3=0);emdf=ehdf.copy();emdf['f2F']=1.;em=zk.emit(emdf,selected_v5=0,selected_c3=0,selected_e3=0)
pref=r*r*h**1.5*Ap**2/(2*f**1.5)
Bump=S*T
rows=[]
for q in [0,1,2,5,10,20,50,100,200,500]:
 d=pd.DataFrame(base)
 for c,z in [('HXX','f2XX'),('HXF','f2XF'),('HXY','f2XY'),('HFF','f2FF'),('HFY','f2FY'),('HYY','f2YY')]:d[z]=S*ip(c)
 # desired extra delta v1 = q*S*T through f2FF only
 mask=(pref>1e-12)&(Bump>0)
 extra=np.zeros(len(u)); extra[mask]=q*Bump[mask]/pref[mask]
 d['f2FF'] += extra
 raw=zk.emit(d,selected_v5=0,selected_c3=0,selected_e3=0)
 tot=pd.DataFrame({'u':u,'x':r,'phi':res.phi,'f':f,'h':h,'phiprime':ph,'A0prime':Ap,'T_H':T,'S_SVT':S})
 for k in SLOTS:
  shared=eh[k].values+S*(em[k].values-eh[k].values)
  tot[k]=mh[k].values+raw[k].values-shared
 for L in [42,1000]:
  a=rr.canonical_audit(tot,L);K=.5*(a['K']+np.swapaxes(a['K'],1,2)); ev=np.linalg.eigvalsh(K)[:,0];neg=ev<0
  rows.append(dict(q=q,L=L,minK=ev.min(),neg_rows=neg.sum(),umin=tot.u[neg].min() if neg.any() else np.nan,umax=tot.u[neg].max() if neg.any() else np.nan))
 print('q',q,[x for x in rows if x['q']==q])
pd.DataFrame(rows).to_csv(B/'ssz_p5_OUTER_V1_CONTROL_SCAN_2026-09-16.csv',index=False)
