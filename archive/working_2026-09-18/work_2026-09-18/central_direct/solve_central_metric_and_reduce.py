from pathlib import Path
import importlib.util, sys, numpy as np, pandas as pd
ROOT=Path('/mnt/data/ssz_work/repo/repo'); B=ROOT/'archive/full_working_snapshot'
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
zk=load('zk',ROOT/'src/ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py')
red=load('red',ROOT/'src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')
prof=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
hj=pd.read_csv(B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
r=raw.x.to_numpy(float); f=raw.f.to_numpy(float); h=raw.h.to_numpy(float); ph=raw.phiprime.to_numpy(float); Ap=raw.A0prime.to_numpy(float)
fp=zk.dr(r,f,1,9,8); hp=zk.dr(r,h,1,9,8)
f2=prof.f2.to_numpy(float).copy(); f2X=prof.f2X.to_numpy(float).copy(); f2F=np.ones(len(raw)); f2Y=np.zeros(len(raw))
f3=prof.f3.to_numpy(float); f3X=prof.f3X_integrated.to_numpy(float); f4=prof.f4.to_numpy(float); f4X=prof.N4.to_numpy(float); f4XX=hj.f4XX_recovered.to_numpy(float); tf4=np.zeros(len(raw))
A=r*r*f; kap=h*ph*ph
# initial residuals
E00 = r*f*hp - (f*(1-h) + r**2*(f*f2-h*Ap**2*f2F)-2*r*h**2*ph*Ap**2*f3+h*Ap**2*(4*(h-1)*f4-h**2*ph**2*f4X))
E11 = r*h*fp - (f*(1-h) + r**2*(f*f2+f*h*ph**2*f2X-h*Ap**2*f2F)-2*r*h**2*ph*Ap**2*(3*f3-h*ph**2*f3X)+h*Ap**2*(4*(3*h-1)*f4-h*(9*h-4)*ph**2*f4X+h**3*ph**4*f4XX))
f2_new=f2+E00/A
f2X_new=f2X+(E11-E00)/(A*kap)
# check
E00n = r*f*hp - (f*(1-h) + r**2*(f*f2_new-h*Ap**2*f2F)-2*r*h**2*ph*Ap**2*f3+h*Ap**2*(4*(h-1)*f4-h**2*ph**2*f4X))
E11n = r*h*fp - (f*(1-h) + r**2*(f*f2_new+f*h*ph**2*f2X_new-h*Ap**2*f2F)-2*r*h**2*ph*Ap**2*(3*f3-h*ph**2*f3X)+h*Ap**2*(4*(3*h-1)*f4-h*(9*h-4)*ph**2*f4X+h**3*ph**4*f4XX))
# action input with same transverse Hessians for first iteration
D=pd.DataFrame({'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,
 'f2X':f2X_new,'f2F':1.0,'f2Y':0.0,'f2XX':prof.HXX,'f2XF':prof.HXF,'f2XY':prof.HXY,'f2FF':prof.HFF,'f2FY':prof.HFY,'f2YY':prof.HYY,
 'f3':prof.f3,'f3X':prof.f3X_integrated,'f3XX':hj.f3XX_selected,'tf3':prof.tilde_f3,
 'f4':prof.f4,'f4X':prof.N4,'f4XX':hj.f4XX_recovered,'f4XXX':hj.f4XXX_selected,'tf4':0.0,
 'f3phi':hj.f3phi,'f3phiX':hj.f3phiX,'f4phi':hj.f4phi,'f4phiX':hj.f4phiX})
out=zk.emit(D)
# use action-derived selected lower order where available
out['v5']=raw.v5
out['c3']=hj.c3_selected
out['e3']=hj.e3_selected
# a5 reholonomize with v5
out['a5']=zk.dr(r,out.a2.to_numpy(float),1,9,8)-zk.dr(r,out.a1.to_numpy(float),2,9,8)-zk.dr(r,Ap*out.v4.to_numpy(float)/2,1,9,8)+Ap*out.v5.to_numpy(float)/2
out.to_csv('/mnt/data/central_metric_resolved_41.csv',index=False)
pd.DataFrame({'u':raw.u,'x':r,'f2_old':f2,'f2_new':f2_new,'f2X_old':f2X,'f2X_new':f2X_new,'E00_old':E00,'E11_old':E11,'E00_new':E00n,'E11_new':E11n}).to_csv('/mnt/data/central_metric_resolve_audit.csv',index=False)
m=(raw.u>=.61)&(raw.u<.71)
print('resolved max E00,E11',np.max(np.abs(E00n[m])),np.max(np.abs(E11n[m])))
print('f2 delta range',np.min((f2_new-f2)[m]),np.max((f2_new-f2)[m]))
print('f2X delta range',np.min((f2X_new-f2X)[m]),np.max((f2X_new-f2X)[m]))
for L in [6,12,20,42,110,420,1000]:
 a=red.canonical_audit(out.loc[m].reset_index(drop=True),L)
 K=np.asarray(a['K']); G=np.asarray(a['G']); Ks=(K+np.swapaxes(K,1,2))/2; Gs=(G+np.swapaxes(G,1,2))/2
 ke=np.linalg.eigvalsh(Ks)
 # generalized eig where K positive: eig invsqrt K G invsqrt K
 cmins=[]
 for i in range(len(K)):
  w,U=np.linalg.eigh(Ks[i])
  if w.min()<=0: cmins.append(np.nan); continue
  Ki=U@np.diag(1/np.sqrt(w))@U.T
  vals=np.linalg.eigvalsh((Ki@Gs[i]@Ki + (Ki@Gs[i]@Ki).T)/2)
  cmins.append(vals.min())
 print('L',L,'Kmin',float(ke.min()),'negK',int(np.sum(ke[:,0]<=0)),'c2min',float(np.nanmin(cmins)),'negc2',int(np.sum(np.array(cmins)<0)))
