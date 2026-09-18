import importlib.util, numpy as np, pandas as pd
from pathlib import Path
from scipy.sparse import lil_matrix,diags, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import linprog
ROOT=Path('/mnt/data/ssz_work/repo/repo')
spec=importlib.util.spec_from_file_location('red',ROOT/'src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'); red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)
raw=pd.read_csv(ROOT/'data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv'); z=pd.read_csv(ROOT/'data/production/ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv'); outer=pd.read_csv(ROOT/'data/prestaged/direct41/outer_resolved_raw_41of41.csv'); inner=pd.read_csv(ROOT/'data/prestaged/direct41/inner_selected_candidate_41of41.csv')
m=(raw.u>=.61)&(raw.u<.71); d=raw.loc[m].reset_index(drop=True); zz=z.loc[m].reset_index(drop=True)
r=d.x.to_numpy(float); u=d.u.to_numpy(float); f=d.f.to_numpy(float); h=d.h.to_numpy(float); ph=d.phiprime.to_numpy(float); Ap=d.A0prime.to_numpy(float); der=lambda y,o=1:red.deriv(r,np.asarray(y,float),o)
fp=der(f); hp=der(h); fpp=der(fp); App=der(Ap); a4=np.sqrt(f*h)/2; N4=zz.N4.to_numpy(float); oldf4=zz.f4.to_numpy(float); alpha1=np.sqrt(h)/(4*np.sqrt(f)); vstart=float(outer.iloc[-1].v6); target=float(inner.iloc[0].v6)
inds,w=red._jet_weights(r,1,9,8); n=len(r); D=lil_matrix((n,n),dtype=float)
for i in range(n):
 for j,ww in zip(inds[i],w[i]):
  if ww: D[i,int(j)]+=float(ww)
D=D.tocsr()
def bump(a,b):
 x=(u-a)/(b-a); y=np.zeros_like(u); q=(x>0)&(x<1); s=2*x[q]-1; y[q]=np.exp(-1/(1-s*s));
 if y.max()>0:y/=y.max()
 return y
# dense overlapping compact supports to give LP flexibility
intervals=[]
for c in np.linspace(.622,.702,17):
 for width in [.016,.024,.032]:
  a=max(.6105,c-width/2); b=min(.7095,c+width/2)
  if b-a>.008: intervals.append((a,b))
bases=[bump(a,b) for a,b in intervals]; B=np.column_stack(bases); p=B.shape[1]
# geometry A matrix independent of f4
Bcoef=r*np.sqrt(f)/(2*h**1.5*Ap); f3_coef=Bcoef/(r*ph); Q_coef=-(4*r*h*ph*f3_coef)/r**2; pref=-np.sqrt(f*h)/(2*r); v10_coef=pref*(r*Q_coef+2*h*ph*f3_coef+(h*fp/f)*Bcoef)
Geom=2*f*f*(r*hp+2*h)+h*r*r*fp*fp-f*r*(r*fp*hp+2*h*(r*fpp+fp)); den=f*f*h*r*r*Ap
Avec=-(f*h*r*r*(Ap*(4*Ap*v10_coef+fp)+f*App))/den
M=(D-diags(Avec)).tolil(); M[0,:]=0;M[0,0]=1;M=M.tocsr()
P1=(1/f)*der(r*np.sqrt(f)/np.sqrt(h))*a4; coefK=Ap*r/np.sqrt(f*h)
def calc(f4):
 a7=(1-(4*h*Ap**2/f)*f4)/(4*r**2*np.sqrt(f*h)); a6=-(np.sqrt(f*h)/(4*r**2))*(1-(h*Ap**2/f)*(4*f4-h*ph**2*N4))
 f3c=(4*f4-h*ph**2*N4)/(r*ph); Qc=-(4*r*h*ph*f3c+8*(1-h)*f4+2*h*h*ph*ph*N4)/r**2; v10c=pref*(r*Qc+2*h*ph*f3c)
 rhs=(a4*Geom-f*h*r*r*(Ap*(4*Ap*v10c)+8*a7*f*f))/den; rhs[0]=vstart; v=spsolve(M,rhs)
 v10=v10c+v10_coef*v; K=2*P1+coefK*v-4*np.sqrt(f*h)*a7*r*r; odd=v*v/16-alpha1*v10
 return v,a7,a6,K,odd,v10
base=calc(oldf4); v0,a70,a60,K0,odd0,v100=base
# affine responses from each bump; solve matrix RHS all at once by repeated direct sparse factorization? scipy factorized
from scipy.sparse.linalg import factorized
solve=factorized(M.tocsc())
# responses derive by calc full; 51 solves fine using factorization
RV=[];RA7=[];RA6=[];RK=[];RV10=[]
for j in range(p):
 v,a7,a6,K,odd,v10=calc(oldf4+B[:,j]); RV.append(v-v0);RA7.append(a7-a70);RA6.append(a6-a60);RK.append(K-K0);RV10.append(v10-v100)
RV=np.column_stack(RV); RA7=np.column_stack(RA7); RA6=np.column_stack(RA6); RK=np.column_stack(RK); RV10=np.column_stack(RV10)

print('baseline K min',K0.min(),'at',u[K0.argmin()]); q=np.where(K0<0)[0]; print('neg count',len(q),'range',u[q[0]] if len(q) else None,u[q[-1]] if len(q) else None); print('a7 min',a70.min(),'a6 max',a60.max(),'odd min',odd0.min());
