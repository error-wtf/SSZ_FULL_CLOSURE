import importlib.util, numpy as np, pandas as pd, math, sys, time
from pathlib import Path
from scipy.optimize import linprog
ROOT=Path('/mnt/data/ssz_work/repo/repo'); A=ROOT/'archive/full_working_snapshot'
spec=importlib.util.spec_from_file_location('mhg',ROOT/'src/ssz_p5/coefficients/mh_general_primitives.py'); mhg=importlib.util.module_from_spec(spec); spec.loader.exec_module(mhg)
spec=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
SLOTS=red.SLOTS
mh0=pd.read_csv(A/'ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
delta=pd.read_csv(A/'ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
u=mh0.u.to_numpy(float); geom=mh0[['u','x','phi','f','h','phiprime','A0prime']].copy(); n=len(u)
r=geom.x.to_numpy(float); f=geom.f.to_numpy(float); h=geom.h.to_numpy(float); ph=geom.phiprime.to_numpy(float); Ap=geom.A0prime.to_numpy(float)
der=lambda y,o=1:red.deriv(r,np.asarray(y,float),o); fp=der(f); hp=der(h)
base_a1=mh0.a1.to_numpy(float); base_c2=mh0.c2.to_numpy(float); base_c4=mh0.c4.to_numpy(float)

def bump(c,w):
 a=c-w; b=c+w; t=(u-a)/(b-a); y=np.zeros(n); q=(t>0)&(t<1); s=2*t[q]-1; y[q]=np.exp(-1/(1-s*s));
 if y.max()>0:y/=y.max()
 return y
centers=np.array([.574,.578,.582,.586,.590,.592,.594,.596,.598,.602,.606])
widths=np.array([.0045,.0045,.0045,.0045,.0045,.004,.004,.004,.0045,.0045,.0035])
B=np.column_stack([bump(c,w) for c,w in zip(centers,widths)]); p=B.shape[1]
# parameters Hlog[p], a1add[p]
def assemble(theta):
 H=np.exp(B@theta[:p]); F=np.ones(n); G=np.ones(n); a1=base_a1+B@theta[p:]
 inp=geom.copy(); inp['H_tensor']=H;inp['F_tensor']=F;inp['G_tensor']=G;inp['a1']=a1;inp['c4']=base_c4;inp['c2']=base_c2
 mh=mhg.emit_from_primitives(inp,regularize_photon_root=True)
 out=geom.copy()
 for k in SLOTS: out[k]=mh[k].to_numpy(float)+delta[k].to_numpy(float)
 out['v5']=0.;out['c3']=0.;out['e3']=0.;out['v12']=-out.v6/(2*out.h);out['v7']=out.v2**2/(4*out.v1)
 out['a5']=der(out.a2.to_numpy(float),1)-der(out.a1.to_numpy(float),2)-der(.5*Ap*out.v4.to_numpy(float),1)
 return out,H,a1

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
 return K
Ls=[6,12,20,42,110,420,1000]
def eval_all(theta):
 d,H,a1=assemble(theta); out={}
 for L in Ls:
  K=exactK(d,L); out[L]=np.linalg.eigvalsh(K)[:,0]
 return d,H,a1,out

theta=np.zeros(2*p)
for it in range(10):
 d,H,a1,E=eval_all(theta)
 print('ITER',it,'mins',[(L,float(E[L].min()),float(u[E[L].argmin()])) for L in Ls],flush=True)
 if min(float(x.min()) for x in E.values())>1e-5: break
 # constraints: 12 worst per L plus all rows below -1 or 0 around bad region
 pairs=[]
 for L in [6,12,42,110,1000]:
  e=E[L]; idx=np.argsort(e)[:20]
  for i in idx: pairs.append((L,int(i)))
 # unique
 pairs=list(dict.fromkeys(pairs)); y0=np.array([E[L][i] for L,i in pairs])
 m=len(theta); J=np.empty((len(pairs),m)); eps=np.r_[np.full(p,0.01),np.full(p,0.01)]
 for j in range(m):
  t=theta.copy();t[j]+=eps[j]
  _,_,_,Ep=eval_all(t)
  J[:,j]=(np.array([Ep[L][i] for L,i in pairs])-y0)/eps[j]
 # LP vars step[m], z[m] abs. require y0+J step >= target, minimize sum z + slight max norm
 target=0.02
 c=np.r_[np.zeros(m),np.ones(m)]
 A_ub=[];b_ub=[]
 # -J step <= y0-target
 for row,yy in zip(J,y0): A_ub.append(np.r_[-row,np.zeros(m)]); b_ub.append(yy-target)
 # abs constraints step-z<=0, -step-z<=0
 for j in range(m):
  row=np.zeros(2*m);row[j]=1;row[m+j]=-1;A_ub.append(row);b_ub.append(0)
  row=np.zeros(2*m);row[j]=-1;row[m+j]=-1;A_ub.append(row);b_ub.append(0)
 bounds=[(-.15,.15)]*p+[(-.15,.15)]*p+[(0,None)]*m
 lp=linprog(c,A_ub=np.array(A_ub),b_ub=np.array(b_ub),bounds=bounds,method='highs')
 if not lp.success:
  print('LP fail',lp.message,flush=True)
  # minimize deficits in linear model with slack variable per constraint
  # vars step,z,slack; target inequalities -Jstep-slack<=y0-target
  nc=len(pairs); cc=np.r_[np.zeros(m),np.ones(m)*.01,np.ones(nc)]
  Aub=[];bub=[]
  for k,(row,yy) in enumerate(zip(J,y0)):
   rr=np.zeros(2*m+nc);rr[:m]=-row;rr[2*m+k]=-1;Aub.append(rr);bub.append(yy-target)
  for j in range(m):
   rr=np.zeros(2*m+nc);rr[j]=1;rr[m+j]=-1;Aub.append(rr);bub.append(0)
   rr=np.zeros(2*m+nc);rr[j]=-1;rr[m+j]=-1;Aub.append(rr);bub.append(0)
  bd=[(-.15,.15)]*m+[(0,None)]*m+[(0,None)]*nc
  lp=linprog(cc,A_ub=np.array(Aub),b_ub=np.array(bub),bounds=bd,method='highs')
  if not lp.success: raise RuntimeError(lp.message)
 step=lp.x[:m]
 # line search actual objective = worst min
 best=None
 for scale in [1,.5,.25,.125]:
  tt=theta+scale*step
  _,_,_,Et=eval_all(tt); worst=min(float(x.min()) for x in Et.values())
  if best is None or worst>best[0]: best=(worst,tt,scale)
 print(' step norm',np.max(abs(step)),'best',best[0],'scale',best[2],flush=True)
 if best[0] <= min(float(x.min()) for x in E.values())+1e-5:
  print('no improvement stop',flush=True);break
 theta=best[1]
d,H,a1,E=eval_all(theta)
print('FINAL theta',theta.tolist(),flush=True)
print('FINAL mins',[(L,float(E[L].min()),float(u[E[L].argmin()])) for L in Ls],flush=True)
d['control_H']=H;d['control_a1']=a1;d['T_H']=delta.T_H;d['S_SVT']=delta.S_SVT;d.to_csv('/mnt/data/ssz_work/outer_sqp_41.csv',index=False)
pd.DataFrame({'parameter':list(range(len(theta))),'value':theta}).to_csv('/mnt/data/ssz_work/outer_sqp_params.csv',index=False)
