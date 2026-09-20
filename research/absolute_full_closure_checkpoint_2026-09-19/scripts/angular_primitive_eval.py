import sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd
from scipy.interpolate import CubicSpline

ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE')
P=Path('/mnt/data/angular_primitive_candidate_41.csv')
d=pd.read_csv(P).copy()
# trusted central bulk, retain a little margin for derivatives
# sorted x ascending for derivative splines
r=d.x.to_numpy(float); u=d.u.to_numpy(float)
order=np.argsort(r); rs=r[order]
def dr(y):
    y=np.asarray(y,float); sp=CubicSpline(rs,y[order]); z=np.empty_like(y); z[order]=sp(rs,1); return z
spec=importlib.util.spec_from_file_location('red',ROOT/'src/ssz_hybrid_jet_aware_constraint_reducer.py')
red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)

def direct_mass(L):
    h0=red.generalized_psi_h0_reduction(d,L,r_col='x')
    D=h0.D_h1_pivot.to_numpy(float); p=h0.p_a1_over_a3.to_numpy(float); q=h0.q_La4_over_a3.to_numpy(float)
    h0v=np.c_[-h0.C_H2/D,-h0.C_phi/D,-d.v2.to_numpy(float)/D]
    h1v=np.c_[-d.a3.to_numpy(float)/D,-h0.B_eff/D,np.zeros(len(d))]
    ep=np.array([1.,0,0]); ef=np.array([0,1.,0]); ev=np.array([0,0,1.])
    H0=np.empty_like(h0v); H1=np.empty_like(h1v)
    for i in range(len(d)):
        H0[i]=ep-q[i]*h0v[i]
        H1[i]=-q[i]*h1v[i]-p[i]*ef
    Hv=np.concatenate([H0,H1],axis=1); hv=np.concatenate([h0v,h1v],axis=1)
    v1=d.v1.to_numpy(float); v1p=dr(v1)
    Cs=np.zeros((len(d),6,6))
    def add(C,coef,a,b=None):
        if b is None:b=a
        C += .5*coef*(np.outer(a,b)+np.outer(b,a))
    for i,c in d.iterrows():
        C=np.zeros((6,6)); H=Hv[i]; hh=hv[i]
        ephi=np.r_[ef,np.zeros(3)]; ephip=np.r_[np.zeros(3),ef]
        eV=np.r_[ev,np.zeros(3)]; eVp=np.r_[np.zeros(3),ev]
        add(C,float(c.c2),H,ephip)
        add(C,float(c.c3+L*c.c4),H,ephi)
        add(C,L*float(c.c5),H,hh)
        add(C,float(c.c6),H,H)
        add(C,L*float(c.d2),hh,ephip)
        add(C,L*float(c.d3),hh,ephi)
        add(C,L*float(c.d4),hh,hh)
        add(C,float(c.e2),ephip,ephip)
        add(C,float(c.e3+L*c.e4),ephi,ephi)
        S0=float(c.v3)*H+float(c.v4)*ephip+float(c.v5)*ephi+L*float(c.v6)*hh
        add(C,-float(c.v1),eV,eV); add(C,1.,eV,S0); add(C,-1/(4*float(c.v1)),S0,S0)
        J0=-2*v1p[i]*eV-2*float(c.v1)*eVp+L*(float(c.v8)*hh+float(c.v12)*H+float(c.v13)*ephi)
        add(C,-1/(4*L*float(c.v9)),J0,J0)
        Cs[i]=C
    A=Cs[:,:3,:3]; B=Cs[:,:3,3:]; Bs=.5*(B+np.swapaxes(B,1,2))
    Bp=np.empty_like(Bs)
    for ia in range(3):
        for ib in range(3): Bp[:,ia,ib]=dr(Bs[:,ia,ib])
    return A-Bp

# exact leading pieces from action data
f=d.f.to_numpy(float); h=d.h.to_numpy(float); ph=d.phiprime.to_numpy(float)
a4=d.a4.to_numpy(float); a6=d.a6.to_numpy(float); v1=d.v1.to_numpy(float); v6=d.v6.to_numpy(float); v9=d.v9.to_numpy(float); v10=d.v10.to_numpy(float); alpha7=d.alpha7.to_numpy(float)
M11_0=-r*r*v6*v6/(4*v1); M13_0=-r*v6/2; M33_0=-v1
# use two L values and Richardson for 1/L expansions
L1=1e7; L2=2e7
print('computing direct masses...', flush=True)
M1=direct_mass(L1); M2=direct_mass(L2)
# for x(L)=x0+x1/L+x2/L^2, x1 approx 2L1*(x(L1)-x(L2)) ; x0 approx 2*x(L2)-x(L1)
M11_1=2*L1*(M1[:,0,0]-M2[:,0,0])
M13_1=2*L1*(M1[:,0,2]-M2[:,0,2])
M33_1=2*L1*(M1[:,2,2]-M2[:,2,2])
# leading finite terms extrapolated
M12_0=2*M2[:,0,1]-M1[:,0,1]
M23_0=2*M2[:,1,2]-M1[:,1,2]
# M22/L = M22_0 + b/L; extrapolate
Y1=M1[:,1,1]/L1; Y2=M2[:,1,1]/L2
M22_0=2*Y2-Y1
# equations
Mcal1=(M11_0*M33_1 - 2*M13_0*M13_1 + M33_0*M11_1)/(4*r*r*f*h*h*alpha7*M33_0)
Mcal2=ph*(r*h*ph*M22_0*M33_0 + 2*M13_0*M23_0 - 2*M12_0*M33_0)/(4*r*f*h*v1*alpha7*M33_0)
K22_0=d.e1.to_numpy(float)
K1=-(2*r*r*a4/f)/(1-f*v6*v6/(8*a4*v10))
K33_0=f*v1*v1*K1/(2*r*r*a4*v10)
cV=-r*r*M33_0/(f*K33_0)
cond55=-(Mcal1+Mcal2+r*r*M22_0/(f*K22_0))
cross=r*r*(4*f*alpha7*Mcal2+ph*ph*M22_0)**2/(4*f*f*ph*ph*alpha7*K22_0)
cond56=r*r*(Mcal1+Mcal2)*M22_0/(f*K22_0)-cross
disc=cond55*cond55-4*cond56
sq=np.sqrt(np.maximum(disc,0)); cminus=.5*(cond55-sq); cplus=.5*(cond55+sq)
mask=(u>=0.62)&(u<0.70)&np.isfinite(cminus)&np.isfinite(cond56)&(disc>=0)
for name,a in [('cV',cV),('cond55',cond55),('cond56',cond56),('disc',disc),('cminus',cminus),('cplus',cplus)]:
    aa=np.asarray(a); j=np.flatnonzero(mask)[np.nanargmin(aa[mask])]; print(name, float(aa[j]), 'u', float(u[j]))
out=pd.DataFrame({'u':u,'x':r,'M11_1':M11_1,'M13_1':M13_1,'M33_1':M33_1,'M12_0':M12_0,'M23_0':M23_0,'M22_0':M22_0,'Mcal1':Mcal1,'Mcal2':Mcal2,'cond55':cond55,'cond56':cond56,'disc':disc,'cV':cV,'cminus':cminus,'cplus':cplus})
out.to_csv('/mnt/data/angular_primitive_candidate_angular.csv',index=False)
