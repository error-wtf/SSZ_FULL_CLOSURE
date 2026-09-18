import importlib.util, numpy as np, pandas as pd, math, time
from pathlib import Path
from scipy.optimize import differential_evolution, minimize
ROOT=Path('/mnt/data/ssz_work/repo/repo'); A=ROOT/'archive/full_working_snapshot'
sp=importlib.util.spec_from_file_location('mhg',ROOT/'src/ssz_p5/coefficients/mh_general_primitives.py'); mhg=importlib.util.module_from_spec(sp); sp.loader.exec_module(mhg)
sp=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(sp); sp.loader.exec_module(red)
SLOTS=red.SLOTS
mh0=pd.read_csv(A/'ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
delta=pd.read_csv(A/'ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
assert np.max(abs(mh0.x-delta.x))<1e-12
u=mh0.u.to_numpy(float); geom=mh0[['u','x','phi','f','h','phiprime','A0prime']].copy(); n=len(u)
base_a1=mh0.a1.to_numpy(float); base_c2=mh0.c2.to_numpy(float); base_c4=mh0.c4.to_numpy(float)
# compact Cinf bumps
def bump(a,b):
 t=(u-a)/(b-a); y=np.zeros(n); q=(t>0)&(t<1); s=2*t[q]-1; y[q]=np.exp(-1/(1-s*s));
 if y.max()>0: y/=y.max()
 return y
intervals=[(.576,.592),(.582,.600),(.588,.606),(.596,.6095)]
B=np.column_stack([bump(a,b) for a,b in intervals]); p=B.shape[1]
r=geom.x.to_numpy(float); f=geom.f.to_numpy(float); h=geom.h.to_numpy(float); ph=geom.phiprime.to_numpy(float); Ap=geom.A0prime.to_numpy(float)
der=lambda y,o=1:red.deriv(r,np.asarray(y,float),o)
fp=der(f); hp=der(h)
Ls=[6,12,20,42,110,420,1000]

def assemble(theta):
 # theta: hlog[p], flog[p], a1add[p]
 hlog=theta[:p]; flog=theta[p:2*p]; aa=theta[2*p:3*p]
 H=np.exp(B@hlog); F=np.exp(B@flog); G=np.ones(n)
 inp=geom.copy(); inp['H_tensor']=H; inp['F_tensor']=F; inp['G_tensor']=G
 inp['a1']=base_a1+B@aa; inp['c4']=base_c4; inp['c2']=base_c2
 mh=mhg.emit_from_primitives(inp,regularize_photon_root=True)
 out=geom.copy();
 for k in SLOTS: out[k]=mh[k].to_numpy(float)+delta[k].to_numpy(float)
 out['v5']=0.;out['c3']=0.;out['e3']=0.;out['v12']=-out.v6/(2*out.h);out['v7']=out.v2**2/(4*out.v1)
 # corrected holonomic a5
 out['a5']=der(out.a2.to_numpy(float),1)-der(out.a1.to_numpy(float),2)-der(.5*Ap*out.v4.to_numpy(float),1)
 return out,H,F

def exactK(d,L):
 a4=d.a4.to_numpy(float);a5=d.a5.to_numpy(float);a6=d.a6.to_numpy(float);b3=d.b3.to_numpy(float);e1=d.e1.to_numpy(float);v1=d.v1.to_numpy(float);v5=d.v5.to_numpy(float);v6=d.v6.to_numpy(float);v10=d.v10.to_numpy(float)
 a4p=der(a4); geom7=(r*(fp*h+f*hp)-2*f*h)/(2*r*f*h*h*ph)*a4-a4p/(h*ph)-Ap*v6/(4*h*ph); alpha7=(a6-geom7)*ph/(2*r*f)
 den1=1-f*v6*v6/(8*a4*v10); K1=-(2*r*r*a4/f)/den1
 den2=fp/f+L/(r*h)-2/r+Ap*v6/(2*a4); K2=(4*r*r*a4/f)/den2
 AA=r*r*f*h; BB=L*der(f**3/(r*r*h))+2*f*(r*fp-2*f)**2/r**3-4*f*f*(2*f*f*alpha7+Ap*Ap*v10)/(r*a4)
 K2p=der(AA)/AA*K2-BB*K2*K2/(8*r*f*f*a4)
 Pp=der(a4)/a4+1/r+L/(2*r*h); C=a5+L*a6-Ap*v5/2; Y=f*b3/(r*a4);Yp=der(Y); YK2p=Yp*K2+Y*K2p;Y2K2p=2*Y*Yp*K2+Y*Y*K2p
 K11=(K1+Pp*K2-.5*K2p)/L; K12=(Y*(K1+.5*Pp*K2)-C*K2/(r*a4)-.5*YK2p)/(2*L); K22=e1+((f*b3/(r*r*a4*a4))*(f*b3*K1-2*C*K2)-.5*Y2K2p)/(4*L); K13=v1/(4*L*r*a4)*(f*v6/v10*K1-2*Ap*K2);K23=f*b3/(2*r*a4)*K13;K33=f*v1*v1/(2*L*r*r*a4*v10)*K1
 K=np.stack([np.stack([K11,K12,K13],-1),np.stack([K12,K22,K23],-1),np.stack([K13,K23,K33],-1)],-2)
 return K,alpha7,den1,den2

def mins(theta, useLs=Ls):
 d,H,F=assemble(theta); vals=[]; meta=[]
 for L in useLs:
  K,a7,d1,d2=exactK(d,L); ev=np.linalg.eigvalsh(K)[:,0]; vals.append(float(np.nanmin(ev))); meta.append((L,float(np.nanmin(ev)),float(u[np.nanargmin(ev)])))
 return vals,meta,d,H,F

def obj(theta):
 try:
  vals,_,d,H,F=mins(theta,[6,20,42,110,1000])
  # log barrier-like quadratic deficits; target small positive margin 1e-5 scaled by L rough
  arr=np.array(vals)
  penalty=np.sum(np.square(np.maximum(0,-arr))) + 0.1*np.sum(np.square(np.maximum(0,1e-6-arr)))
  regu=1e-4*np.dot(theta,theta)
  # keep H,F sensible >0 automatically; penalize extremes
  regu+=1e-3*(np.mean((np.log(H))**2)+np.mean((np.log(F))**2))
  return float(penalty+regu)
 except Exception as e:
  return 1e12

x0=np.zeros(3*p)
print('BASE',mins(x0,[6,20,42,110,1000])[0],flush=True)
# sensitivity probes
for j in range(3*p):
 for amp in [0.1,-0.1] if j<2*p else [0.05,-0.05]:
  x=x0.copy();x[j]=amp
  vals=mins(x,[6,42,1000])[0]
  print('SENS',j,amp,vals,flush=True)
# differential evolution bounded
bounds=[(-1.2,1.2)]*(2*p)+[(-2.0,2.0)]*p
res=differential_evolution(obj,bounds,maxiter=20,popsize=6,tol=1e-4,polish=False,workers=1,updating='immediate',seed=3,disp=True)
print('DE',res.fun,res.x,flush=True); print('MINS',mins(res.x)[1],flush=True)
# polish
r2=minimize(obj,res.x,method='Nelder-Mead',options={'maxiter':300,'xatol':1e-5,'fatol':1e-8,'disp':True})
print('NM',r2.fun,r2.x,flush=True); vals,meta,d,H,F=mins(r2.x); print('FINAL',meta,flush=True)
d['T_H']=delta.T_H;d['S_SVT']=delta.S_SVT;d['control_H']=H;d['control_F']=F;d['construction']='direct outer controlled Horndeski primitives + fixed genuine-SVT delta; Cinf compact interior controls; endpoints unchanged'
d.to_csv('/mnt/data/ssz_work/outer_direct_controlled_41.csv',index=False)
pd.DataFrame({'param':np.arange(len(r2.x)),'value':r2.x}).to_csv('/mnt/data/ssz_work/outer_control_params.csv',index=False)
