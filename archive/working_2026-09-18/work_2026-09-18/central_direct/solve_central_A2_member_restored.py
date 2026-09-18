from pathlib import Path
import importlib.util, numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
ROOT=Path('/mnt/data/ssz_work/repo/repo'); B=ROOT/'archive/full_working_snapshot'
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
zk=load('zk',ROOT/'src/ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py')
red=load('red',ROOT/'src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py')
prof=pd.read_csv(B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
raw=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
hj=pd.read_csv(B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
# full arrays, JET geometry derivatives
r=raw.x.to_numpy(float); u=raw.u.to_numpy(float); f=raw.f.to_numpy(float); h=raw.h.to_numpy(float); ph=raw.phiprime.to_numpy(float); Ap=raw.A0prime.to_numpy(float)
fp=zk.dr(r,f,1,9,8); fpp=zk.dr(r,f,2,9,8); hp=zk.dr(r,h,1,9,8); App=zk.dr(r,Ap,1,9,8)
f4=prof.f4.to_numpy(float); f4X=prof.N4.to_numpy(float); a4=raw.a4.to_numpy(float); f3X=prof.f3X_integrated.to_numpy(float); f4XX=hj.f4XX_recovered.to_numpy(float)
# central production subset but integrate on full .61-.715 for smoothness; sorted r increasing
mask=(u>=.61)&(u<=.715)
idx=np.where(mask)[0]; order=idx[np.argsort(r[idx])]
rs=r[order]
arrs={name:PchipInterpolator(rs, arr[order], extrapolate=False) for name,arr in [('f',f),('h',h),('ph',ph),('Ap',Ap),('fp',fp),('fpp',fpp),('hp',hp),('App',App),('f4',f4),('f4X',f4X),('a4',a4)]}
# formulas from chosen f4/f4X, with tf4=f2Y=0, JA=0
def aux(rr, vv):
 ff=float(arrs['f'](rr)); hh=float(arrs['h'](rr)); pp=float(arrs['ph'](rr)); aa=float(arrs['Ap'](rr)); F4=float(arrs['f4'](rr)); F4X=float(arrs['f4X'](rr));
 C=2*hh**1.5*aa/(rr*np.sqrt(ff))
 F3=(vv/C + 4*F4-hh*pp*pp*F4X)/(rr*pp)
 F2F=-(4*rr*hh*pp*F3 + 8*(1-hh)*F4 + 2*hh*hh*pp*pp*F4X)/(rr*rr)
 FP=float(arrs['fp'](rr))
 V10=-(np.sqrt(ff*hh)/(2*rr))*(rr*F2F+2*hh*pp*F3+(hh*FP/ff)*(rr*pp*F3-4*F4+hh*pp*pp*F4X))
 A7=(1-(4*hh*aa*aa/ff)*F4)/(4*rr*rr*np.sqrt(ff*hh))
 return F3,F2F,V10,A7

def ode(rr,y):
 vv=float(y[0]); ff=float(arrs['f'](rr)); hh=float(arrs['h'](rr)); aa=float(arrs['Ap'](rr)); FP=float(arrs['fp'](rr)); FPP=float(arrs['fpp'](rr)); HP=float(arrs['hp'](rr)); APP=float(arrs['App'](rr)); A4=float(arrs['a4'](rr));
 _,_,V10,A7=aux(rr,vv)
 num=A4*(2*ff*ff*(rr*HP+2*hh)+hh*rr*rr*FP*FP-ff*rr*(rr*FP*HP+2*hh*(rr*FPP+FP)))-ff*hh*rr*rr*(aa*(4*V10*aa+vv*FP)+ff*vv*APP+8*A7*ff*ff)
 den=ff*ff*hh*rr*rr*aa
 return [num/den]
# outer boundary at u=.61 = max r
io=np.argmin(abs(u-.61)); r0=r[io]; v0=raw.v6.iloc[io]
r1=rs.min()
sol=solve_ivp(ode,(r0,r1),[v0],rtol=2e-11,atol=2e-12,dense_output=True,max_step=2e-4)
print('solve success',sol.success,sol.message,'nfev',sol.nfev,'r0,r1',r0,r1,'v0,vend',v0,sol.y[0,-1])
v6new=np.full(len(raw),np.nan); v6new[idx]=sol.sol(r[idx])[0]
# action profiles derived pointwise
f3new=np.full(len(raw),np.nan); f2Fnew=np.full(len(raw),np.nan); v10new=np.full(len(raw),np.nan); alpha7new=np.full(len(raw),np.nan)
for i in idx:
 f3new[i],f2Fnew[i],v10new[i],alpha7new[i]=aux(r[i],v6new[i])
# solve metric f2/f2X exactly with fixed f3X/f4XX
f2old=prof.f2.to_numpy(float); f2Xold=prof.f2X.to_numpy(float)
A=r*r*f; kap=h*ph*ph
# residual based on provisional old f2/f2X but new first-order action
E00 = r*f*hp - (f*(1-h)+r*r*(f*f2old-h*Ap*Ap*f2Fnew)-2*r*h*h*ph*Ap*Ap*f3new+h*Ap*Ap*(4*(h-1)*f4-h*h*ph*ph*f4X))
E11 = r*h*fp - (f*(1-h)+r*r*(f*f2old+f*h*ph*ph*f2Xold-h*Ap*Ap*f2Fnew)-2*r*h*h*ph*Ap*Ap*(3*f3new-h*ph*ph*f3X)+h*Ap*Ap*(4*(3*h-1)*f4-h*(9*h-4)*ph*ph*f4X+h**3*ph**4*f4XX))
f2new=f2old+E00/A; f2Xnew=f2Xold+(E11-E00)/(A*kap)
# scalar/tangent higher jets: keep background-null Hessians initially; f3phi holonomic from new f3
Xp=zk.dr(r,raw.X.to_numpy(float),1,9,8)
f3phi_new=(zk.dr(r,np.nan_to_num(f3new,nan=0.0),1,9,8)-f3X*Xp)/ph
# only central values relevant; edges derivative polluted from nan fill -> replace via local derivative on idx in r order separately
# recompute f3phi using full central interpolation values
f3tmp=f3new[idx]; rtmp=r[idx]; Xtmp=raw.X.to_numpy(float)[idx]; phtmp=ph[idx]
f3p=zk.dr(rtmp,f3tmp,1,9,8); Xp2=zk.dr(rtmp,Xtmp,1,9,8); f3phi_new[idx]=(f3p-f3X[idx]*Xp2)/phtmp
D=pd.DataFrame({'u':raw.u,'x':raw.x,'phi':raw.phi,'f':raw.f,'h':raw.h,'phiprime':raw.phiprime,'A0prime':raw.A0prime,'X':raw.X,
'f2X':f2Xnew,'f2F':f2Fnew,'f2Y':0.0,'f2XX':prof.HXX,'f2XF':prof.HXF,'f2XY':prof.HXY,'f2FF':prof.HFF,'f2FY':prof.HFY,'f2YY':prof.HYY,
'f3':f3new,'f3X':f3X,'f3XX':hj.f3XX_selected,'tf3':prof.tilde_f3,'f4':prof.f4,'f4X':prof.N4,'f4XX':hj.f4XX_recovered,'f4XXX':hj.f4XXX_selected,'tf4':0.0,
'f3phi':f3phi_new,'f3phiX':hj.f3phiX,'f4phi':hj.f4phi,'f4phiX':hj.f4phiX})
# emit only central subset separately to avoid nan contamination
central_idx=np.where((u>=.61)&(u<.71))[0]
dc=D.iloc[central_idx].reset_index(drop=True)
out0=zk.emit(dc)
# background-null Hessian right inverse restores trusted v1,v4,c2 target profiles
tv1=raw.v1.iloc[central_idx].to_numpy(); tv4=raw.v4.iloc[central_idx].to_numpy(); tc2=raw.c2.iloc[central_idx].to_numpy()
responses=[]
for col in ['f2XX','f2XF','f2FF']:
    dd=dc.copy(); dd[col]=dd[col].to_numpy()+1.0; gg=zk.emit(dd)
    responses.append(np.stack([gg.v1.to_numpy()-out0.v1.to_numpy(),gg.v4.to_numpy()-out0.v4.to_numpy(),gg.c2.to_numpy()-out0.c2.to_numpy()],axis=1))
deltas=np.zeros((len(dc),3))
for jj in range(len(dc)):
    M=np.column_stack([responses[k][jj] for k in range(3)])
    tgt=np.array([tv1[jj]-out0.v1.iloc[jj],tv4[jj]-out0.v4.iloc[jj],tc2[jj]-out0.c2.iloc[jj]])
    deltas[jj]=np.linalg.solve(M,tgt)
for kk,col in enumerate(['f2XX','f2XF','f2FF']): dc[col]=dc[col].to_numpy()+deltas[:,kk]
out=zk.emit(dc)
print('principal target residuals',np.max(np.abs(out.v1-tv1)),np.max(np.abs(out.v4-tv4)),np.max(np.abs(out.c2-tc2)),'max Hessian deltas',np.max(np.abs(deltas),axis=0))
# lower order use action historical selected for now
out['v5']=raw.v5.iloc[central_idx].to_numpy(); out['c3']=hj.c3_selected.iloc[central_idx].to_numpy(); out['e3']=hj.e3_selected.iloc[central_idx].to_numpy()
rr=out.x.to_numpy(float); Apc=out.A0prime.to_numpy(float)
out['a5']=zk.dr(rr,out.a2.to_numpy(float),1,9,8)-zk.dr(rr,out.a1.to_numpy(float),2,9,8)-zk.dr(rr,Apc*out.v4.to_numpy(float)/2,1,9,8)+Apc*out.v5.to_numpy(float)/2
out.to_csv('/mnt/data/central_A2_resolved_restored_41.csv',index=False)
# audits: background residuals / JA / A2 direct
f_c=f[central_idx];h_c=h[central_idx];ph_c=ph[central_idx];Ap_c=Ap[central_idx];r_c=r[central_idx];fp_c=fp[central_idx];hp_c=hp[central_idx];
f2_c=f2new[central_idx]; f2X_c=f2Xnew[central_idx]; f2F_c=f2Fnew[central_idx]; f3_c=f3new[central_idx]; f3X_c=f3X[central_idx]; f4_c=f4[central_idx]; f4X_c=f4X[central_idx]; f4XX_c=f4XX[central_idx]
E00n=r_c*f_c*hp_c-(f_c*(1-h_c)+r_c**2*(f_c*f2_c-h_c*Ap_c**2*f2F_c)-2*r_c*h_c**2*ph_c*Ap_c**2*f3_c+h_c*Ap_c**2*(4*(h_c-1)*f4_c-h_c**2*ph_c**2*f4X_c))
E11n=r_c*h_c*fp_c-(f_c*(1-h_c)+r_c**2*(f_c*f2_c+f_c*h_c*ph_c**2*f2X_c-h_c*Ap_c**2*f2F_c)-2*r_c*h_c**2*ph_c*Ap_c**2*(3*f3_c-h_c*ph_c**2*f3X_c)+h_c*Ap_c**2*(4*(3*h_c-1)*f4_c-h_c*(9*h_c-4)*ph_c**2*f4X_c+h_c**3*ph_c**4*f4XX_c))
JA=np.sqrt(h_c/f_c)*Ap_c*(r_c*r_c*f2F_c+4*r_c*h_c*ph_c*f3_c+8*(1-h_c)*f4_c+2*h_c*h_c*ph_c*ph_c*f4X_c)
v6e=out.v6.to_numpy(float); v6p=zk.dr(r_c,v6e,1,9,8)
# A2 rhs on emitted values / alpha7 emitter
A4=out.a4.to_numpy(float); V10=out.v10.to_numpy(float); A7=out.alpha7.to_numpy(float); App_c=App[central_idx]; fpp_c=fpp[central_idx]
rhs=(A4*(2*f_c*f_c*(r_c*hp_c+2*h_c)+h_c*r_c*r_c*fp_c*fp_c-f_c*r_c*(r_c*fp_c*hp_c+2*h_c*(r_c*fpp_c+fp_c)))-f_c*h_c*r_c*r_c*(Ap_c*(4*V10*Ap_c+v6e*fp_c)+f_c*v6e*App_c+8*A7*f_c*f_c))/(f_c*f_c*h_c*r_c*r_c*Ap_c)
print('background max',np.max(abs(E00n)),np.max(abs(E11n)),np.max(abs(JA)))
print('A2 residual median/p95/max',np.median(abs(v6p-rhs)),np.quantile(abs(v6p-rhs),.95),np.max(abs(v6p-rhs)))
print('v6 range',v6e.min(),v6e.max(),'f2F range',np.nanmin(f2F_c),np.nanmax(f2F_c),'f3 range',np.nanmin(f3_c),np.nanmax(f3_c))
# reducer gates
rows=[]
for L in [6,12,20,42,110,420,1000]:
 a=red.canonical_audit(out,L); K=np.asarray(a['K']); G=np.asarray(a['G']); Ks=(K+K.transpose(0,2,1))/2; Gs=(G+G.transpose(0,2,1))/2; ke=np.linalg.eigvalsh(Ks)
 cmin=1e99; negc=0
 for i in range(len(K)):
  w,U=np.linalg.eigh(Ks[i])
  if w.min()<=0: continue
  Ki=U@np.diag(1/np.sqrt(w))@U.T; C=(Ki@Gs[i]@Ki); vals=np.linalg.eigvalsh((C+C.T)/2); cmin=min(cmin,vals.min()); negc += int(vals.min()<=0)
 rows.append((L,float(ke.min()),int((ke[:,0]<=0).sum()),float(cmin),negc,float(a['max_h0_quadratic'] if 'max_h0_quadratic' in a else np.nan)))
 print('L',rows[-1])
pd.DataFrame(rows,columns=['L','Kmin','negK','c2min','negc2','h0']).to_csv('/mnt/data/central_A2_restored_operator_audit.csv',index=False)
# source profiles audit
pd.DataFrame({'u':u[central_idx],'x':r_c,'v6':v6new[central_idx],'f3':f3_c,'f2F':f2F_c,'f2':f2_c,'f2X':f2X_c,'E00':E00n,'E11':E11n,'JA':JA,'A2res':v6p-rhs}).to_csv('/mnt/data/central_A2_restored_background_audit.csv',index=False)
