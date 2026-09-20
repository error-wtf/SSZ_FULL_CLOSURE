from pathlib import Path
import importlib.util
import numpy as np, pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

ROOT=Path('/mnt/data/ssz_final_work/SSZ_FULL_CLOSURE')
P=ROOT/'data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_41.csv'
d=pd.read_csv(P).copy().reset_index(drop=True)
r=d.x.to_numpy(float); u=d.u.to_numpy(float); n=len(d); order=np.argsort(r); rs=r[order]

def dr(y):
    y=np.asarray(y,float); sp=CubicSpline(rs,y[order]); z=np.empty_like(y); z[order]=sp(rs,1); return z

spec=importlib.util.spec_from_file_location('red',ROOT/'src/ssz_hybrid_jet_aware_constraint_reducer.py')
red=importlib.util.module_from_spec(spec); spec.loader.exec_module(red)

ep=np.array([1.,0.,0.]); ef=np.array([0.,1.,0.]); ev=np.array([0.,0.,1.])
ephi=np.r_[ef,np.zeros(3)]; ephip=np.r_[np.zeros(3),ef]
eV=np.r_[ev,np.zeros(3)]; eVp=np.r_[np.zeros(3),ev]

def sym_outer_batch(a,b):
    return .5*(np.einsum('ni,nj->nij',a,b)+np.einsum('ni,nj->nij',b,a))

def prepare_L(L):
    h0=red.generalized_psi_h0_reduction(d,L,r_col='x')
    D=h0.D_h1_pivot.to_numpy(float); pp=h0.p_a1_over_a3.to_numpy(float); qq=h0.q_La4_over_a3.to_numpy(float)
    h0v=np.c_[-h0.C_H2/D,-h0.C_phi/D,-d.v2.to_numpy(float)/D]
    h1v=np.c_[-d.a3.to_numpy(float)/D,-h0.B_eff/D,np.zeros(n)]
    H0=ep[None,:]-qq[:,None]*h0v
    H1=-qq[:,None]*h1v-pp[:,None]*ef[None,:]
    H=np.concatenate([H0,H1],axis=1); hh=np.concatenate([h0v,h1v],axis=1)
    v1=d.v1.to_numpy(float); v1p=dr(v1)
    J0=(-2*v1p[:,None]*eV[None,:]-2*v1[:,None]*eVp[None,:]
        +L*(d.v8.to_numpy(float)[:,None]*hh+d.v12.to_numpy(float)[:,None]*H+d.v13.to_numpy(float)[:,None]*ephi[None,:]))
    # baseline direct mass, vectorized all terms
    C=np.zeros((n,6,6))
    def add_var(coef,a,b):
        nonlocal C
        C += coef[:,None,None]*sym_outer_batch(a,b)
    add_var(d.c2.to_numpy(float),H,np.tile(ephip,(n,1)))
    add_var(d.c3.to_numpy(float)+L*d.c4.to_numpy(float),H,np.tile(ephi,(n,1)))
    add_var(L*d.c5.to_numpy(float),H,hh)
    add_var(d.c6.to_numpy(float),H,H)
    add_var(L*d.d2.to_numpy(float),hh,np.tile(ephip,(n,1)))
    add_var(L*d.d3.to_numpy(float),hh,np.tile(ephi,(n,1)))
    add_var(L*d.d4.to_numpy(float),hh,hh)
    add_var(d.e2.to_numpy(float),np.tile(ephip,(n,1)),np.tile(ephip,(n,1)))
    add_var(d.e3.to_numpy(float)+L*d.e4.to_numpy(float),np.tile(ephi,(n,1)),np.tile(ephi,(n,1)))
    S0=(d.v3.to_numpy(float)[:,None]*H+d.v4.to_numpy(float)[:,None]*np.tile(ephip,(n,1))
        +d.v5.to_numpy(float)[:,None]*np.tile(ephi,(n,1))+L*d.v6.to_numpy(float)[:,None]*hh)
    add_var(-d.v1.to_numpy(float),np.tile(eV,(n,1)),np.tile(eV,(n,1)))
    add_var(np.ones(n),np.tile(eV,(n,1)),S0)
    add_var(-1/(4*d.v1.to_numpy(float)),S0,S0)
    add_var(-1/(4*L*d.v9.to_numpy(float)),J0,J0)
    A=C[:,:3,:3]; B=C[:,:3,3:]; Bs=.5*(B+np.swapaxes(B,1,2)); Bp=np.empty_like(Bs)
    for i in range(3):
        for j in range(3): Bp[:,i,j]=dr(Bs[:,i,j])
    Mbase=A-Bp
    return dict(L=L,H=H,hh=hh,J0=J0,S0=S0,Cbase=C,Mbase=Mbase,Dbase=D,p=pp,q=qq,Cphi=h0.C_phi.to_numpy(float),CH2=h0.C_H2.to_numpy(float),Beff=h0.B_eff.to_numpy(float))

L1=1e7; L2=2e7
PRE1=prepare_L(L1); PRE2=prepare_L(L2)

# fixed background/principal arrays
f=d.f.to_numpy(float); h=d.h.to_numpy(float); ph=d.phiprime.to_numpy(float)
a4=d.a4.to_numpy(float); v1=d.v1.to_numpy(float); v6=d.v6.to_numpy(float); v10=d.v10.to_numpy(float); alpha7=d.alpha7.to_numpy(float)
M11_0=-r*r*v6*v6/(4*v1); M13_0=-r*v6/2; M33_0=-v1
K22_0=d.e1.to_numpy(float)
K1=-(2*r*r*a4/f)/(1-f*v6*v6/(8*a4*v10)); K33_0=f*v1*v1*K1/(2*r*r*a4*v10)
cV=-r*r*M33_0/(f*K33_0)

# C-infinity compact bumps
centers=np.array([0.657,0.667,0.677,0.687,0.696,0.703])
widths=np.array([0.015,0.014,0.014,0.013,0.011,0.007])
def cinf_bump(center,width):
    z=(u-center)/width; out=np.zeros_like(u); m=np.abs(z)<1
    q=1-z[m]**2; raw=np.exp(-1/q)
    out[m]=raw/np.exp(-1.0)
    return out
BAS=np.column_stack([cinf_bump(c,w) for c,w in zip(centers,widths)])
# enforce exactly zero outside central+inner-handover domain
BAS[(u<=0.62)|(u>=0.709985),:]=0
slots=['c4','c5','d3','e4','v13','v5','a9']; ns=len(slots); nb=BAS.shape[1]

# scale each slot so coefficients O(1) have comparable response; empirical native magnitudes
SCALE={'c4':10.0,'c5':10.0,'d3':50.0,'e4':100.0,'v13':200.0,'v5':100.0,'a9':10.0}

def deltas_from_p(p):
    p=np.asarray(p,float).reshape(ns,nb)
    return {slot:SCALE[slot]*(BAS@p[i]) for i,slot in enumerate(slots)}

def mass_with_delta(pre,dd):
    L=pre['L']
    # a9 changes only the algebraic H0/h1 pivot at this stage.
    D=pre['Dbase']+L*dd['a9']
    h0v=np.c_[-pre['CH2']/D,-pre['Cphi']/D,-d.v2.to_numpy(float)/D]
    h1v=np.c_[-d.a3.to_numpy(float)/D,-pre['Beff']/D,np.zeros(n)]
    H0=ep[None,:]-pre['q'][:,None]*h0v
    H1=-pre['q'][:,None]*h1v-pre['p'][:,None]*ef[None,:]
    H=np.concatenate([H0,H1],axis=1); hh=np.concatenate([h0v,h1v],axis=1)
    EPHI=np.tile(ephi,(n,1)); EPHIP=np.tile(ephip,(n,1)); EV=np.tile(eV,(n,1)); EVP=np.tile(eVp,(n,1))
    C=np.zeros((n,6,6))
    def addv(coef,a,b):
        nonlocal C
        C += np.asarray(coef,float)[:,None,None]*sym_outer_batch(a,b)
    addv(d.c2.to_numpy(float),H,EPHIP)
    addv(d.c3.to_numpy(float)+L*(d.c4.to_numpy(float)+dd['c4']),H,EPHI)
    addv(L*(d.c5.to_numpy(float)+dd['c5']),H,hh)
    addv(d.c6.to_numpy(float),H,H)
    addv(L*d.d2.to_numpy(float),hh,EPHIP)
    addv(L*(d.d3.to_numpy(float)+dd['d3']),hh,EPHI)
    addv(L*d.d4.to_numpy(float),hh,hh)
    addv(d.e2.to_numpy(float),EPHIP,EPHIP)
    addv(d.e3.to_numpy(float)+L*(d.e4.to_numpy(float)+dd['e4']),EPHI,EPHI)
    S0=(d.v3.to_numpy(float)[:,None]*H+d.v4.to_numpy(float)[:,None]*EPHIP
        +(d.v5.to_numpy(float)+dd['v5'])[:,None]*EPHI+L*d.v6.to_numpy(float)[:,None]*hh)
    addv(-d.v1.to_numpy(float),EV,EV); addv(np.ones(n),EV,S0); addv(-1/(4*d.v1.to_numpy(float)),S0,S0)
    v1p=dr(d.v1.to_numpy(float))
    J0=(-2*v1p[:,None]*EV-2*d.v1.to_numpy(float)[:,None]*EVP
        +L*(d.v8.to_numpy(float)[:,None]*hh+d.v12.to_numpy(float)[:,None]*H+(d.v13.to_numpy(float)+dd['v13'])[:,None]*EPHI))
    addv(-1/(4*L*d.v9.to_numpy(float)),J0,J0)
    A=C[:,:3,:3]; B=C[:,:3,3:]; Bs=.5*(B+np.swapaxes(B,1,2)); Bp=np.empty_like(Bs)
    for i in range(3):
        for j in range(3): Bp[:,i,j]=dr(Bs[:,i,j])
    return A-Bp

def angular(p):
    dd=deltas_from_p(p)
    A=mass_with_delta(PRE1,dd); B=mass_with_delta(PRE2,dd)
    M11_1=2*L1*(A[:,0,0]-B[:,0,0]); M13_1=2*L1*(A[:,0,2]-B[:,0,2]); M33_1=2*L1*(A[:,2,2]-B[:,2,2])
    M12_0=2*B[:,0,1]-A[:,0,1]; M23_0=2*B[:,1,2]-A[:,1,2]
    M22_0=2*(B[:,1,1]/L2)-(A[:,1,1]/L1)
    Mcal1=(M11_0*M33_1-2*M13_0*M13_1+M33_0*M11_1)/(4*r*r*f*h*h*alpha7*M33_0)
    Mcal2=ph*(r*h*ph*M22_0*M33_0+2*M13_0*M23_0-2*M12_0*M33_0)/(4*r*f*h*v1*alpha7*M33_0)
    c55=-(Mcal1+Mcal2+r*r*M22_0/(f*K22_0))
    cross=r*r*(4*f*alpha7*Mcal2+ph*ph*M22_0)**2/(4*f*f*ph*ph*alpha7*K22_0)
    c56=r*r*(Mcal1+Mcal2)*M22_0/(f*K22_0)-cross
    disc=c55*c55-4*c56; sq=np.sqrt(np.maximum(disc,0)); cm=.5*(c55-sq); cp=.5*(c55+sq)
    return dict(c55=c55,c56=c56,disc=disc,cminus=cm,cplus=cp,cV=cV,delta=dd)

TRUST=(u>=0.622)&(u<0.70)&np.isfinite(cV)

if __name__ == '__main__':
    # compare zero to separately computed baseline
    z=np.zeros(ns*nb); a0=angular(z)
    print('baseline minima')
    for k in ['cminus','cplus','c55','c56','disc','cV']:
        arr=a0[k]; j=np.flatnonzero(TRUST)[np.nanargmin(arr[TRUST])]; print(k,float(arr[j]),'u',float(u[j]))
    
    # finite sensitivity of each slot/bump at current worst cminus point; write table
    j0=np.flatnonzero(TRUST)[np.nanargmin(a0['cminus'][TRUST])]
    rows=[]
    for si,s in enumerate(slots):
        for bi,c in enumerate(centers):
            p=np.zeros(ns*nb); p[si*nb+bi]=1e-3
            aa=angular(p)
            rows.append([s,float(c), (aa['cminus'][j0]-a0['cminus'][j0])/1e-3,
                         (aa['c55'][j0]-a0['c55'][j0])/1e-3,(aa['c56'][j0]-a0['c56'][j0])/1e-3])
    pd.DataFrame(rows,columns=['slot','center','dcminus_dp','dc55_dp','dc56_dp']).to_csv('/mnt/data/angular_mass_control_jacobian.csv',index=False)
    print('wrote jacobian')
