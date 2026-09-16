from pathlib import Path
import sys, importlib.util, numpy as np, pandas as pd
from scipy.integrate import cumulative_trapezoid
B=Path('/mnt/data')
# load modules
sp=importlib.util.spec_from_file_location('mh',B/'ssz_p5_mh_luminal_g4phi_emitter_JET9D8_2026-09-16.py'); mh=importlib.util.module_from_spec(sp); sp.loader.exec_module(mh)
sp=importlib.util.spec_from_file_location('red',B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(sp); sp.loader.exec_module(red)
car=pd.read_csv(B/'ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv').sort_values('r_over_rs').reset_index(drop=True)
old=pd.read_csv(B/'ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv').sort_values('x',ascending=False).reset_index(drop=True)
# r descending from .71 endpoint inward
r=old.x.to_numpy(float); f=old.f.to_numpy(float); h=old.h.to_numpy(float); ph=old.phi_r.to_numpy(float)
K0=0.009742; P1=(1+K0)/2
# match strong carrier endpoint
j=np.argmin(abs(car.u-0.71)); c0=car.iloc[j]
mu0=float(c0.mu); Z0=r[0]**2*np.sqrt(f[0]/h[0])/mu0
Z=Z0+np.r_[0,cumulative_trapezoid(P1*np.sqrt(f/h),r)]
mu=r**2*np.sqrt(f/h)/Z
sq=np.sqrt(f*h); a4=.5*sq
a1=sq*(mu/2-r)/ph
G3X=2*a1/(sq*h*ph**2*r**2)
# choose several c2 continuations; start carrier endpoint constant / old core / blend
c2start=float(pd.read_csv(B/'ssz_p5_F3_horndeski_carrier_SELECTED_41of41_2026-09-15.csv').sort_values('u').iloc[-1].c2)
print('Zmin',Z.min(),'mu_end',mu[-1],'G3X endpoints',G3X[0],G3X[-1],'c2start',c2start)
inp=pd.DataFrame({'u':old.u,'x':r,'phi':old.phi,'f':f,'h':h,'phi_r':ph,'A0prime':0.0})
profiles={
 'const':np.full(len(old),c2start),
 'old':old.c2.to_numpy(float),
 'zero':np.zeros(len(old)),
}
for name,c2 in profiles.items():
 d=mh.emit(inp,H=np.ones(len(old)),G4phi=np.zeros(len(old)),G3X=G3X,c2_profile=c2)
 # normalize naming expected is x + slots, okay
 # epsilon Y production deformation
 kap=h*ph**2; ZA=1-.02*kap
 d['v1']*=ZA; d['v10']*=ZA
 # ensure v7 canonical identity (zero vector => v2=0, already zero)
 d['v7']=d['v2']**2/(4*d['v1'])
 out=B/f'test_cubic_core_{name}.csv'; d.to_csv(out,index=False)
 for L in [6,42,1000]:
  a=red.canonical_audit(d,float(L)); K=a['K'];G=a['G']
  ke=[];ge=[]
  for i in range(len(d)):
   ks=(K[i]+K[i].T)/2; gs=(G[i]+G[i].T)/2
   w,V=np.linalg.eigh(ks); ke.append(w[0])
   if w[0]>0:
    W=V@np.diag(1/np.sqrt(w))@V.T
    ge.append(np.linalg.eigvalsh(W@gs@W)[0])
   else: ge.append(np.nan)
  print(name,'L',L,'Kmin',np.nanmin(ke),'at',d.u.iloc[int(np.nanargmin(ke))],'Ggenmin',np.nanmin(ge),'diag',a['diagnostics'])
