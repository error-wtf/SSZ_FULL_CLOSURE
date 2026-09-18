import importlib.util, numpy as np, pandas as pd
from pathlib import Path
ROOT=Path('/mnt/data/ssz_work/repo/repo');A=ROOT/'archive/full_working_snapshot'
sp=importlib.util.spec_from_file_location('mhg',ROOT/'src/ssz_p5/coefficients/mh_general_primitives.py');mhg=importlib.util.module_from_spec(sp);sp.loader.exec_module(mhg)
sp=importlib.util.spec_from_file_location('red','/mnt/data/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py');red=importlib.util.module_from_spec(sp);sp.loader.exec_module(red)
mh0=pd.read_csv(A/'ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True); delta=pd.read_csv(A/'ssz_p5_OUTER_GENUINE_SVT_DELTA_41of41_2026-09-16.csv').sort_values('x').reset_index(drop=True)
u=mh0.u.to_numpy(float);geom=mh0[['u','x','phi','f','h','phiprime','A0prime']].copy();r=geom.x.to_numpy(float);f=geom.f.to_numpy(float);h=geom.h.to_numpy(float);ph=geom.phiprime.to_numpy(float);Ap=geom.A0prime.to_numpy(float);der=lambda y,o=1:red.deriv(r,np.asarray(y,float),o);fp=der(f);hp=der(h);SLOTS=red.SLOTS

def bump(a,b):
 t=(u-a)/(b-a);y=np.zeros(len(u));q=(t>0)&(t<1);s=2*t[q]-1;y[q]=np.exp(-1/(1-s*s));y/=max(y.max(),1e-300);return y
B=bump(.572,.608)

def run(a1amp=0,hlog=0,flog=0):
 inp=geom.copy();inp['H_tensor']=np.exp(hlog*B);inp['F_tensor']=np.exp(flog*B);inp['G_tensor']=1.;inp['a1']=mh0.a1+B*a1amp;inp['c4']=mh0.c4;inp['c2']=mh0.c2
 mh=mhg.emit_from_primitives(inp,regularize_photon_root=True);d=geom.copy()
 for k in SLOTS:d[k]=mh[k].to_numpy(float)+delta[k].to_numpy(float)
 d['v5']=0.;d['c3']=0.;d['e3']=0.;d['v12']=-d.v6/(2*d.h);d['v7']=d.v2**2/(4*d.v1);d['a5']=der(d.a2,1)-der(d.a1,2)-der(.5*Ap*d.v4,1)
 def Kmin(L):
  a4=d.a4.to_numpy(float);a5=d.a5.to_numpy(float);a6=d.a6.to_numpy(float);b3=d.b3.to_numpy(float);e1=d.e1.to_numpy(float);v1=d.v1.to_numpy(float);v5=d.v5.to_numpy(float);v6=d.v6.to_numpy(float);v10=d.v10.to_numpy(float);a4p=der(a4);geom7=(r*(fp*h+f*hp)-2*f*h)/(2*r*f*h*h*ph)*a4-a4p/(h*ph)-Ap*v6/(4*h*ph);alpha7=(a6-geom7)*ph/(2*r*f);den1=1-f*v6*v6/(8*a4*v10);K1=-(2*r*r*a4/f)/den1;den2=fp/f+L/(r*h)-2/r+Ap*v6/(2*a4);K2=(4*r*r*a4/f)/den2;AA=r*r*f*h;BB=L*der(f**3/(r*r*h))+2*f*(r*fp-2*f)**2/r**3-4*f*f*(2*f*f*alpha7+Ap*Ap*v10)/(r*a4);K2p=der(AA)/AA*K2-BB*K2*K2/(8*r*f*f*a4);Pp=der(a4)/a4+1/r+L/(2*r*h);C=a5+L*a6-Ap*v5/2;Y=f*b3/(r*a4);Yp=der(Y);YK2p=Yp*K2+Y*K2p;Y2K2p=2*Y*Yp*K2+Y*Y*K2p;K11=(K1+Pp*K2-.5*K2p)/L;K12=(Y*(K1+.5*Pp*K2)-C*K2/(r*a4)-.5*YK2p)/(2*L);K22=e1+((f*b3/(r*r*a4*a4))*(f*b3*K1-2*C*K2)-.5*Y2K2p)/(4*L);K13=v1/(4*L*r*a4)*(f*v6/v10*K1-2*Ap*K2);K23=f*b3/(2*r*a4)*K13;K33=f*v1*v1/(2*L*r*r*a4*v10)*K1;K=np.stack([np.stack([K11,K12,K13],-1),np.stack([K12,K22,K23],-1),np.stack([K13,K23,K33],-1)],-2);e=np.linalg.eigvalsh(K)[:,0];i=e.argmin();return float(e[i]),float(u[i])
 return [Kmin(L) for L in [6,42,1000]]
for typ in ['a1','H','F']:
 print('\n',typ,flush=True)
 vals=[-2,-1,-.5,-.2,-.1,-.05,0,.05,.1,.2,.5,1,2] if typ!='a1' else [-5,-2,-1,-.5,-.2,-.1,-.05,0,.05,.1,.2,.5,1,2,5]
 for a in vals:
  kw={typ.lower()+'log':a} if typ in ['H','F'] else {'a1amp':a}
  if typ=='H':kw={'hlog':a}
  if typ=='F':kw={'flog':a}
  try: print(a,run(**kw),flush=True)
  except Exception as e: print('ERR',a,e,flush=True)
