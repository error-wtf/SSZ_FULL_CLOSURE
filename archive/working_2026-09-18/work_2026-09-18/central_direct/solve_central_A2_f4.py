from pathlib import Path
import importlib.util, numpy as np, pandas as pd
ROOT=Path('/mnt/data/ssz_work/repo/repo'); B=ROOT/'archive/full_working_snapshot'
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
zk=load('zk',ROOT/'src/ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py')
red=load('red',ROOT/'src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')
prof=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
hj=pd.read_csv(B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
r=raw.x.to_numpy(float);u=raw.u.to_numpy(float);f=raw.f.to_numpy(float);h=raw.h.to_numpy(float);ph=raw.phiprime.to_numpy(float);Ap=raw.A0prime.to_numpy(float);X=raw.X.to_numpy(float)
fp=zk.dr(r,f,1,9,8);fpp=zk.dr(r,f,2,9,8);hp=zk.dr(r,h,1,9,8);App=zk.dr(r,Ap,1,9,8);Xp=zk.dr(r,X,1,9,8)
a4=raw.a4.to_numpy(float);v6=raw.v6.to_numpy(float);v6p=zk.dr(r,v6,1,9,8);f4X=prof.N4.to_numpy(float);f3X=prof.f3X_integrated.to_numpy(float);f4XX=hj.f4XX_recovered.to_numpy(float)
mask=(u>=.61)&(u<.71); ids=np.where(mask)[0]
# residual of A2 for arbitrary f4 value, fixed f4X and target v6
def vals(i,F4):
 rr=r[i];ff=f[i];hh=h[i];pp=ph[i];aa=Ap[i];FP=fp[i];FPP=fpp[i];HP=hp[i];APP=App[i];A4=a4[i];V6=v6[i];F4X=f4X[i]
 C=2*hh**1.5*aa/(rr*np.sqrt(ff)); F3=(V6/C+4*F4-hh*pp*pp*F4X)/(rr*pp)
 F2F=-(4*rr*hh*pp*F3+8*(1-hh)*F4+2*hh*hh*pp*pp*F4X)/(rr*rr)
 V10=-(np.sqrt(ff*hh)/(2*rr))*(rr*F2F+2*hh*pp*F3+(hh*FP/ff)*(rr*pp*F3-4*F4+hh*pp*pp*F4X))
 A7=(1-(4*hh*aa*aa/ff)*F4)/(4*rr*rr*np.sqrt(ff*hh))
 num=A4*(2*ff*ff*(rr*HP+2*hh)+hh*rr*rr*FP*FP-ff*rr*(rr*FP*HP+2*hh*(rr*FPP+FP)))-ff*hh*rr*rr*(aa*(4*V10*aa+V6*FP)+ff*V6*APP+8*A7*ff*ff)
 rhs=num/(ff*ff*hh*rr*rr*aa)
 return F3,F2F,V10,A7,v6p[i]-rhs
f4new=np.full(len(raw),np.nan); f3new=np.full(len(raw),np.nan); f2Fnew=np.full(len(raw),np.nan); v10new=np.full(len(raw),np.nan)
for i in ids:
 q0=vals(i,0.0); q1=vals(i,1.0); slope=q1[-1]-q0[-1]
 if abs(slope)<1e-12: raise RuntimeError(('tiny slope',i,slope))
 F4=-q0[-1]/slope
 q=vals(i,F4); f4new[i]=F4;f3new[i]=q[0];f2Fnew[i]=q[1];v10new[i]=q[2]
# metric solve f2,f2X
f2old=prof.f2.to_numpy(float);f2Xold=prof.f2X.to_numpy(float);A=r*r*f;kap=h*ph*ph
E00=r*f*hp-(f*(1-h)+r*r*(f*f2old-h*Ap*Ap*f2Fnew)-2*r*h*h*ph*Ap*Ap*f3new+h*Ap*Ap*(4*(h-1)*f4new-h*h*ph*ph*f4X))
E11=r*h*fp-(f*(1-h)+r*r*(f*f2old+f*h*ph*ph*f2Xold-h*Ap*Ap*f2Fnew)-2*r*h*h*ph*Ap*Ap*(3*f3new-h*ph*ph*f3X)+h*Ap*Ap*(4*(3*h-1)*f4new-h*(9*h-4)*ph*ph*f4X+h**3*ph**4*f4XX))
f2new=f2old+E00/A;f2Xnew=f2Xold+(E11-E00)/(A*kap)
# use holonomic f4phi/f3phi on the central subset; f4X is independent jet
rr=r[ids]; f4p=zk.dr(rr,f4new[ids],1,9,8); Xpc=zk.dr(rr,X[ids],1,9,8); php=ph[ids]
f4phi=(f4p-f4X[ids]*Xpc)/php
f3p=zk.dr(rr,f3new[ids],1,9,8); f3phi=(f3p-f3X[ids]*Xpc)/php
# Build initial action and emit with OLD Hessians
base=pd.DataFrame({'u':u[ids],'x':r[ids],'phi':raw.phi.iloc[ids].to_numpy(),'f':f[ids],'h':h[ids],'phiprime':ph[ids],'A0prime':Ap[ids],'X':X[ids],
'f2X':f2Xnew[ids],'f2F':f2Fnew[ids],'f2Y':0.0,'f2XX':prof.HXX.iloc[ids].to_numpy(),'f2XF':prof.HXF.iloc[ids].to_numpy(),'f2XY':prof.HXY.iloc[ids].to_numpy(),'f2FF':prof.HFF.iloc[ids].to_numpy(),'f2FY':prof.HFY.iloc[ids].to_numpy(),'f2YY':prof.HYY.iloc[ids].to_numpy(),
'f3':f3new[ids],'f3X':f3X[ids],'f3XX':hj.f3XX_selected.iloc[ids].to_numpy(),'tf3':prof.tilde_f3.iloc[ids].to_numpy(),'f4':f4new[ids],'f4X':f4X[ids],'f4XX':hj.f4XX_recovered.iloc[ids].to_numpy(),'f4XXX':hj.f4XXX_selected.iloc[ids].to_numpy(),'tf4':0.0,'f3phi':f3phi,'f3phiX':hj.f3phiX.iloc[ids].to_numpy(),'f4phi':f4phi,'f4phiX':hj.f4phiX.iloc[ids].to_numpy()})
out0=zk.emit(base)
# Solve HXX,HXF,HFF to exactly restore old targets v1,v4,c2 while keep HXY/HFY/HYY fixed.
tv1=raw.v1.iloc[ids].to_numpy(); tv4=raw.v4.iloc[ids].to_numpy(); tc2=raw.c2.iloc[ids].to_numpy()
# finite-difference exact linear map columns for HXX,HXF,HFF around current base. Use analytic linearity via unit emissions in batches is costly but okay 3 emits.
def emitted_with(col,delta):
 d=base.copy(); d[col]=d[col].to_numpy()+delta; return zk.emit(d)
# compute unit responses globally
responses=[]
for col in ['f2XX','f2XF','f2FF']:
 g=emitted_with(col,1.0); responses.append(np.stack([g.v1.to_numpy()-out0.v1.to_numpy(),g.v4.to_numpy()-out0.v4.to_numpy(),g.c2.to_numpy()-out0.c2.to_numpy()],axis=1))
# solve per point 3x3
Deltas=np.zeros((len(ids),3)); bad=0
for j in range(len(ids)):
 M=np.column_stack([responses[k][j] for k in range(3)])
 tgt=np.array([tv1[j]-out0.v1.iloc[j],tv4[j]-out0.v4.iloc[j],tc2[j]-out0.c2.iloc[j]])
 try:Deltas[j]=np.linalg.solve(M,tgt)
 except np.linalg.LinAlgError: Deltas[j]=np.linalg.lstsq(M,tgt,rcond=None)[0];bad+=1
for k,col in enumerate(['f2XX','f2XF','f2FF']): base[col]=base[col].to_numpy()+Deltas[:,k]
out=zk.emit(base)
# lower order use action-derived existing selected (doesn't affect principal)
out['v5']=raw.v5.iloc[ids].to_numpy();out['c3']=hj.c3_selected.iloc[ids].to_numpy();out['e3']=hj.e3_selected.iloc[ids].to_numpy();out['a5']=zk.dr(rr,out.a2.to_numpy(float),1,9,8)-zk.dr(rr,out.a1.to_numpy(float),2,9,8)-zk.dr(rr,Ap[ids]*out.v4.to_numpy(float)/2,1,9,8)+Ap[ids]*out.v5.to_numpy(float)/2
# background checks
F2=f2new[ids];F2X=f2Xnew[ids];F2F=f2Fnew[ids];F3=f3new[ids];F4=f4new[ids];F4X=f4X[ids];F3X=f3X[ids];F4XX=f4XX[ids];fc=f[ids];hc=h[ids];pc=ph[ids];Ac=Ap[ids];fpc=fp[ids];hpc=hp[ids]
E00n=rr*fc*hpc-(fc*(1-hc)+rr**2*(fc*F2-hc*Ac**2*F2F)-2*rr*hc**2*pc*Ac**2*F3+hc*Ac**2*(4*(hc-1)*F4-hc**2*pc**2*F4X))
E11n=rr*hc*fpc-(fc*(1-hc)+rr**2*(fc*F2+fc*hc*pc**2*F2X-hc*Ac**2*F2F)-2*rr*hc**2*pc*Ac**2*(3*F3-hc*pc**2*F3X)+hc*Ac**2*(4*(3*hc-1)*F4-hc*(9*hc-4)*pc**2*F4X+hc**3*pc**4*F4XX))
JA=np.sqrt(hc/fc)*Ac*(rr**2*F2F+4*rr*hc*pc*F3+8*(1-hc)*F4+2*hc**2*pc**2*F4X)
# A2 residual recomputed from solved f4
res=[]
for i in ids: res.append(vals(i,f4new[i])[-1])
print('f4 old/new ranges',prof.f4.iloc[ids].min(),prof.f4.iloc[ids].max(),np.nanmin(F4),np.nanmax(F4),'delta median/max',np.median(abs(F4-prof.f4.iloc[ids].to_numpy())),np.max(abs(F4-prof.f4.iloc[ids].to_numpy())))
print('f2F range',np.nanmin(F2F),np.nanmax(F2F),'f3 range',np.nanmin(F3),np.nanmax(F3))
print('bg max',np.max(abs(E00n)),np.max(abs(E11n)),np.max(abs(JA)),'A2 point residual max',np.max(np.abs(np.asarray(res))),'control solve bad',bad)
print('target restoration',np.max(abs(out.v1-tv1)),np.max(abs(out.v4-tv4)),np.max(abs(out.c2-tc2)))
print('H deltas max',np.max(abs(Deltas),axis=0))
out.to_csv('/mnt/data/central_A2_f4_resolved_41.csv',index=False); base.to_csv('/mnt/data/central_A2_f4_action_jets.csv',index=False)
pd.DataFrame({'u':u[ids],'x':rr,'f4_old':prof.f4.iloc[ids].to_numpy(),'f4_new':F4,'f3_new':F3,'f2F_new':F2F,'f2_new':F2,'f2X_new':F2X,'E00':E00n,'E11':E11n,'JA':JA,'A2_res':res}).to_csv('/mnt/data/central_A2_f4_background_audit.csv',index=False)
# reducer gates
rows=[]
for L in [6,12,20,42,110,420,1000]:
 a=red.canonical_audit(out,L); K=np.asarray(a['K']);G=np.asarray(a['G']);Ks=(K+K.transpose(0,2,1))/2;Gs=(G+G.transpose(0,2,1))/2;ke=np.linalg.eigvalsh(Ks)
 cvals=[]
 for j in range(len(K)):
  w,U=np.linalg.eigh(Ks[j]);
  if w.min()<=0: cvals.append(np.nan);continue
  Ki=U@np.diag(1/np.sqrt(w))@U.T; C=Ki@Gs[j]@Ki;cvals.append(np.linalg.eigvalsh((C+C.T)/2).min())
 rows.append((L,float(ke.min()),int(np.sum(ke[:,0]<=0)),float(np.nanmin(cvals)),int(np.sum(np.array(cvals)<0))))
 print('L',rows[-1])
pd.DataFrame(rows,columns=['L','Kmin','negK','c2min','negc2']).to_csv('/mnt/data/central_A2_f4_operator_audit.csv',index=False)
