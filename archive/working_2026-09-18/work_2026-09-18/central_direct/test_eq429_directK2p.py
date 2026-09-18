import numpy as np,pandas as pd,importlib.util
from pathlib import Path
P=Path('/mnt/data')
d=pd.read_csv(P/'ssz_p5_CENTRAL_SELECTED_41of41_CORRECTED_A5_V12_2026-09-16.csv').copy()
hj_spec=importlib.util.spec_from_file_location('hj',P/'ssz_p5_higher_jet_closure_2026-09-16.py');hj=importlib.util.module_from_spec(hj_spec);hj_spec.loader.exec_module(hj)
red_spec=importlib.util.spec_from_file_location('red',P/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py');red=importlib.util.module_from_spec(red_spec);red_spec.loader.exec_module(red)
r=d.x.to_numpy(float);f=d.f.to_numpy(float);h=d.h.to_numpy(float);Ap=d.A0prime.to_numpy(float)
der=lambda y,o=1:hj.local_poly_deriv(r,np.asarray(y,float),order=o,window=9,degree=8)
fp=der(f)
def Kpub(L):
 a4=d.a4.to_numpy(float); a5=d.a5.to_numpy(float); a6=d.a6.to_numpy(float); b3=d.b3.to_numpy(float); e1=d.e1.to_numpy(float)
 v1=d.v1.to_numpy(float);v5=d.v5.to_numpy(float);v6=d.v6.to_numpy(float);v10=d.v10.to_numpy(float)
 K1=-(2*r*r*a4/f)/(1-f*v6*v6/(8*a4*v10))
 K2=(4*r*r*a4/f)/(fp/f+L/(r*h)-2/r+Ap*v6/(2*a4))
 K2p=der(K2)
 Pp=der(a4)/a4+1/r+L/(2*r*h);C=a5+L*a6-Ap*v5/2;Y=f*b3/(r*a4);Yp=der(Y)
 YK2p=Yp*K2+Y*K2p;Y2K2p=2*Y*Yp*K2+Y*Y*K2p
 K11=(K1+Pp*K2-.5*K2p)/L
 K12=(Y*(K1+.5*Pp*K2)-C*K2/(r*a4)-.5*YK2p)/(2*L)
 K22=e1+((f*b3/(r*r*a4*a4))*(f*b3*K1-2*C*K2)-.5*Y2K2p)/(4*L)
 K13=v1/(4*L*r*a4)*(f*v6/v10*K1-2*Ap*K2)
 K23=f*b3/(2*r*a4)*K13
 K33=f*v1*v1/(2*L*r*r*a4*v10)*K1
 return np.stack([np.stack([K11,K12,K13],-1),np.stack([K12,K22,K23],-1),np.stack([K13,K23,K33],-1)],-2)
for L in [6,12,20,42,110,420,1000]:
 K=Kpub(L); ev=np.linalg.eigvalsh(K)
 A=red.canonical_audit(d,float(L)); Kr=.5*(A['K']+np.swapaxes(A['K'],1,2))
 diff=np.max(np.abs(K-Kr),axis=(1,2))
 print(L,'pub min',ev[:,0].min(),'reducer min',np.linalg.eigvalsh(Kr)[:,0].min(),'maxdiff',diff.max(),'meddiff',np.median(diff))
