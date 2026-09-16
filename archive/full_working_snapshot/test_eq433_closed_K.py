import numpy as np,pandas as pd,importlib.util
from pathlib import Path
P=Path('/mnt/data')
d=pd.read_csv(P/'ssz_p5_CENTRAL_SELECTED_41of41_CORRECTED_A5_V12_2026-09-16.csv').copy()
z=pd.read_csv(P/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv')
sp=importlib.util.spec_from_file_location('hj',P/'ssz_p5_higher_jet_closure_2026-09-16.py');hj=importlib.util.module_from_spec(sp);sp.loader.exec_module(hj)
r=d.x.to_numpy(float);f=d.f.to_numpy(float);h=d.h.to_numpy(float);Ap=d.A0prime.to_numpy(float)
der=lambda y,o=1:hj.local_poly_deriv(r,np.asarray(y,float),order=o,window=9,degree=8)
fp=der(f);alpha7=np.interp(d.u,z.u,z.alpha7)

def K433(L):
 a4=d.a4.to_numpy(float); a5=d.a5.to_numpy(float); a6=d.a6.to_numpy(float); b3=d.b3.to_numpy(float); e1=d.e1.to_numpy(float)
 v1=d.v1.to_numpy(float);v5=d.v5.to_numpy(float);v6=d.v6.to_numpy(float);v10=d.v10.to_numpy(float)
 K1=-(2*r*r*a4/f)/(1-f*v6*v6/(8*a4*v10))
 K2=(4*r*r*a4/f)/(fp/f+L/(r*h)-2/r+Ap*v6/(2*a4))
 A=r*r*f*h
 B=L*der(f**3/(r*r*h))+2*f*(r*fp-2*f)**2/r**3-4*f*f*(2*f*f*alpha7+Ap*Ap*v10)/(r*a4)
 K2p=der(A)/A*K2-B*K2*K2/(8*r*f*f*a4)
 Pp=der(a4)/a4+1/r+L/(2*r*h);C=a5+L*a6-Ap*v5/2;Y=f*b3/(r*a4);Yp=der(Y)
 YK2p=Yp*K2+Y*K2p;Y2K2p=2*Y*Yp*K2+Y*Y*K2p
 K11=(K1+Pp*K2-.5*K2p)/L
 K12=(Y*(K1+.5*Pp*K2)-C*K2/(r*a4)-.5*YK2p)/(2*L)
 K22=e1+((f*b3/(r*r*a4*a4))*(f*b3*K1-2*C*K2)-.5*Y2K2p)/(4*L)
 K13=v1/(4*L*r*a4)*(f*v6/v10*K1-2*Ap*K2);K23=f*b3/(2*r*a4)*K13;K33=f*v1*v1/(2*L*r*r*a4*v10)*K1
 return np.stack([np.stack([K11,K12,K13],-1),np.stack([K12,K22,K23],-1),np.stack([K13,K23,K33],-1)],-2)
rows=[]
for L in [6,12,20,42,110,420,1000]:
 K=K433(L); ev=np.linalg.eigvalsh(K); rows.append([L,ev[:,0].min(),ev[:,1].min(),ev[:,2].min(),np.mean(ev[:,0]<=0)])
 print(L,'mins',ev.min(axis=0),'badfrac',np.mean(ev[:,0]<=0))
for U in [.61,.6666667,.706135,.715]:
 j=np.argmin(abs(d.u-U));print('\nU',d.u.iloc[j]);
 for L in [6,42]: print('L',L,'eig',np.linalg.eigvalsh(K433(L)[j]))
pd.DataFrame(rows,columns=['L','min_eig1','min_eig2','min_eig3','bad_fraction']).to_csv(P/'ssz_p5_ZK_EQ433_CLOSED_FINITE_L_K_2026-09-16.csv',index=False)
