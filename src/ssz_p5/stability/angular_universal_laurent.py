"""Formal universal large-L angular reducer from the common 41-slot action.

This module mirrors the accepted finite-L profile-aware reducer but carries a
truncated Laurent series in eps=1/L through:

    common 41-slot action -> constraint maps -> radial product rules
    -> Euler reduction -> K(eps), M(eps) -> characteristic Laurent polynomial.

It is deliberately theory-agnostic.  Published genuine-SVT and
Maxwell-Horndeski formulas are *oracles*, not hard-coded repair targets.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import comb
from pathlib import Path
from typing import Dict, Mapping
import importlib.util

import numpy as np
import pandas as pd

from ssz_p5.paths import paths as B
from ssz_p5.jets.jet9d8 import profile_derivative
from .laurent import LaurentSeries, outer_product, derivative as lderivative

SLOTS=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],
       *[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],
       *[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]

PMIN=-6
PMAX=14


def C(a): return LaurentSeries.const(np.asarray(a,float),pmin=PMIN,pmax=PMAX)
def E(p,a): return LaurentSeries.monomial(p,np.asarray(a,float),pmin=PMIN,pmax=PMAX)
def L(a): return E(-1,a)
def invL(a): return E(1,a)


def dr(r,a,order=1,window=9,degree=8):
    return profile_derivative(np.asarray(r,float),np.asarray(a,float),order,window,degree)


def ds(r,s:LaurentSeries,order=1,window=9,degree=8):
    return lderivative(s,lambda a: dr(r,a,order,window,degree))


def _series_arrays(df):
    return {k:C(df[k].to_numpy(float)) for k in SLOTS}


def constraint_maps_laurent(df,window=9,degree=8):
    d=df.copy().reset_index(drop=True)
    d['v12']=-d['v6']/(2*d['h'])
    r=d.x.to_numpy(float); n=len(d)
    q=_series_arrays(d)
    one=np.ones(n)

    p=q['a1']/q['a3']
    qq=L(one)*(q['a4']/q['a3'])
    pp=ds(r,p,1,window,degree); qrp=ds(r,qq,1,window,degree)

    Braw=q['a2']-q['v2']*q['v4']/(C(2*one)*q['v1'])
    Cphi=q['a5']+L(d['a6'].to_numpy(float))-q['v2']*q['v5']/(C(2*one)*q['v1'])
    CH2=q['a7']+L(d['a8'].to_numpy(float))-q['v2']*q['v3']/(C(2*one)*q['v1'])
    Ch1=L(d['a9'].to_numpy(float)-d['v2'].to_numpy(float)*d['v6'].to_numpy(float)/(2*d['v1'].to_numpy(float)))
    Beff=Braw-q['a3']*pp-CH2*p
    Dh1=Ch1-q['a3']*qrp-CH2*qq

    # vector-valued helpers: multiplication by scalar series is explicit
    def smul(s:LaurentSeries,v:LaurentSeries): return s*v
    z=np.zeros((n,3)); epsi=z.copy(); epsi[:,0]=1.; ephi=z.copy(); ephi[:,1]=1.; eV=z.copy(); eV[:,2]=1.
    H10 = -(CH2/Dh1) * C(epsi) -(Cphi/Dh1)*C(ephi) -(q['v2']/Dh1)*C(eV)
    H11 = -(q['a3']/Dh1)*C(epsi) -(Beff/Dh1)*C(ephi)
    H20 = C(epsi)-qq*H10
    H21 = -qq*H11-p*C(ephi)

    J10=q['b3']*C(ephi)+q['b4']*H20+L(d['b5'].to_numpy(float))*H10
    J11=q['b2']*C(ephi)+q['b4']*H21+L(d['b5'].to_numpy(float))*H11
    R10=invL(-one)*J10; R11=invL(-one)*J11
    R20=C(d['v6'].to_numpy(float)/2)*H10 + invL(-2*d['v1'].to_numpy(float))*C(eV)
    R21=C(d['v6'].to_numpy(float)/2)*H11

    Delta=C(4*d['b1'].to_numpy(float)*d['v10'].to_numpy(float)-d['v11'].to_numpy(float)**2)
    T10=(C(2*d['v10'].to_numpy(float))*R10-C(d['v11'].to_numpy(float))*R20)/Delta
    T11=(C(2*d['v10'].to_numpy(float))*R11-C(d['v11'].to_numpy(float))*R21)/Delta
    A10=(-C(d['v11'].to_numpy(float))*R10+C(2*d['b1'].to_numpy(float))*R20)/Delta
    A11=(-C(d['v11'].to_numpy(float))*R11+C(2*d['b1'].to_numpy(float))*R21)/Delta

    v1p=dr(r,d['v1'].to_numpy(float),1,window,degree)
    A00=-(q['v8']*H10+q['v12']*H20)/(C(2*d['v9'].to_numpy(float)))
    A01=-(q['v8']*H11+q['v12']*H21)/(C(2*d['v9'].to_numpy(float)))
    A00=A00 + C(-d['v13'].to_numpy(float)/(2*d['v9'].to_numpy(float)))*C(ephi)
    A00=A00 + invL(v1p/d['v9'].to_numpy(float))*C(eV)
    A01=A01 + invL(d['v1'].to_numpy(float)/d['v9'].to_numpy(float))*C(eV)

    return dict(p=p,q=qq,Dh1=Dh1,DeltaV=Delta,
        H10=H10,H11=H11,H20=H20,H21=H21,
        H1dot0=T10,H1dot1=T11,dA1dot0=A10,dA1dot1=A11,
        dA00=A00,dA01=A01)


def field_maps_laurent(df,window=9,degree=8):
    m=constraint_maps_laurent(df,window,degree)
    n=len(df); z=np.zeros((n,3)); ephi=z.copy(); ephi[:,1]=1.; eV=z.copy(); eV[:,2]=1.
    return {
      'dphi': {(0,0):C(ephi)}, 'V':{(0,0):C(eV)},
      'h1': {(0,0):m['H10'],(0,1):m['H11']},
      'H2': {(0,0):m['H20'],(0,1):m['H21']},
      'H1': {(1,0):m['H1dot0'],(1,1):m['H1dot1']},
      'dA1':{(1,0):m['dA1dot0'],(1,1):m['dA1dot1']},
      'dA0':{(0,0):m['dA00'],(0,1):m['dA01']},
    },m


def action_terms_laurent(df):
    d=df.copy().reset_index(drop=True); d['v12']=-d['v6']/(2*d['h'])
    q={k:d[k].to_numpy(float) for k in SLOTS}; n=len(d)
    T=[]
    def add(c,A,da,Bn,db):
        if not isinstance(c,LaurentSeries): c=C(c)
        T.append((c,A,da,Bn,db))
    # metric/scalar non-H0 block
    add(L(q['b1']),'H1',(0,0),'H1',(0,0)); add(q['b2'],'H1',(0,0),'dphi',(1,1))
    add(q['b3'],'H1',(0,0),'dphi',(1,0)); add(q['b4'],'H1',(0,0),'H2',(1,0))
    add(L(q['b5']),'H1',(0,0),'h1',(1,0)); add(q['c1'],'dphi',(1,0),'H2',(1,0))
    add(q['c2'],'H2',(0,0),'dphi',(0,1)); add(C(q['c3'])+L(q['c4']),'H2',(0,0),'dphi',(0,0))
    add(L(q['c5']),'H2',(0,0),'h1',(0,0)); add(q['c6'],'H2',(0,0),'H2',(0,0))
    add(L(q['d1']),'h1',(1,0),'h1',(1,0)); add(L(q['d2']),'h1',(0,0),'dphi',(0,1))
    add(L(q['d3']),'h1',(0,0),'dphi',(0,0)); add(L(q['d4']),'h1',(0,0),'h1',(0,0))
    add(q['e1'],'dphi',(1,0),'dphi',(1,0)); add(q['e2'],'dphi',(0,1),'dphi',(0,1))
    add(C(q['e3'])+L(q['e4']),'dphi',(0,0),'dphi',(0,0))
    # auxiliary vector form
    add(2*q['v1'],'V',(0,0),'dA0',(0,1)); add(-2*q['v1'],'V',(0,0),'dA1',(1,0)); add(-q['v1'],'V',(0,0),'V',(0,0))
    parts=[('H2',C(q['v3'])),('dphi_r',C(q['v4'])),('dphi',C(q['v5'])),('h1',L(q['v6']))]
    pmap={'H2':('H2',(0,0)),'dphi_r':('dphi',(0,1)),'dphi':('dphi',(0,0)),'h1':('h1',(0,0))}
    for name,c in parts:
        f,dd=pmap[name]; add(c,'V',(0,0),f,dd)
    for ia,(na,ca) in enumerate(parts):
        fa,da=pmap[na]
        for ib in range(ia,len(parts)):
            nb,cb=parts[ib]; fb,db=pmap[nb]; fac=2.0 if ib!=ia else 1.0
            add((ca*cb).scale(-fac/(4*q['v1'])) ,fa,da,fb,db)
    add(L(0.5*q['v6']),'h1',(0,0),'dA1',(1,0)); add(L(q['v8']),'h1',(0,0),'dA0',(0,0))
    add(L(q['v9']),'dA0',(0,0),'dA0',(0,0)); add(L(q['v10']),'dA1',(0,0),'dA1',(0,0))
    add(L(q['v11']),'H1',(0,0),'dA1',(0,0)); add(L(q['v12']),'H2',(0,0),'dA0',(0,0)); add(L(q['v13']),'dphi',(0,0),'dA0',(0,0))
    return T


def apply_derivative_laurent(op,dt,drn,r,window=9,degree=8):
    out={}
    for (pt,pr),a in op.items():
        for k in range(drn+1):
            ak=ds(r,a,k,window,degree) if k else a
            key=(pt+dt,pr+drn-k); v=ak.scale(comb(drn,k))
            out[key]=out[key]+v if key in out else v
    return out


def adjoint_product_laurent(A,c,Bop,r,window=9,degree=8):
    out={}
    for (pt,pr),av in A.items():
      for (qt,qr),bv in Bop.items():
        base=outer_product(av,bv,c)
        sgn=(-1.0)**(pt+pr)
        for k in range(pr+1):
            bk=ds(r,base,k,window,degree) if k else base
            key=(pt+qt,pr-k+qr); v=bk.scale(sgn*comb(pr,k))
            out[key]=out[key]+v if key in out else v
    return out


def reduced_operator_laurent(df,window=9,degree=8):
    d=df.copy().reset_index(drop=True); d['v12']=-d['v6']/(2*d['h']); r=d.x.to_numpy(float)
    F,m=field_maps_laurent(d,window,degree); P={}
    for c,A,da,Bn,db in action_terms_laurent(d):
        Aop=apply_derivative_laurent(F[A],da[0],da[1],r,window,degree)
        Bop=apply_derivative_laurent(F[Bn],db[0],db[1],r,window,degree)
        for X in (adjoint_product_laurent(Aop,c,Bop,r,window,degree),adjoint_product_laurent(Bop,c,Aop,r,window,degree)):
            for k,v in X.items(): P[k]=P[k]+v if k in P else v
    return P,m


def canonical_laurent(df,window=9,degree=8):
    P,m=reduced_operator_laurent(df,window,degree); r=df.x.to_numpy(float); n=len(df)
    zero=LaurentSeries.const(np.zeros((n,3,3)),pmin=PMIN,pmax=PMAX)
    get=lambda k:P.get(k,zero)
    P20=get((2,0)); P02=get((0,2)); P01=get((0,1)); P00=get((0,0))
    K=P20.scale(-0.5); G=P02.scale(0.5)
    P02p=ds(r,P02,1,window,degree); S=(P02p-P01).scale(0.5); Sp=ds(r,S,1,window,degree)
    M=(P00+Sp).scale(-0.5)
    return dict(P=P,K=K,G=G,S=S,M=M,maps=m)


def evaluate_blocks(a,eps):
    return {k:a[k].evaluate(eps) for k in ('K','G','S','M')}


def compare_finite_reducer(df,Lvals=(1e3,1e4),window=9,degree=8,mask=None):
    """Regression of formal Laurent blocks against the accepted finite-L reducer.

    ``mask`` restricts the reported norm while retaining the full radial profile
    for all product-rule derivatives.  This matters near handover boundaries,
    where the asymptotic witness itself is not the certified oracle region.
    """
    pth=B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'
    spec=importlib.util.spec_from_file_location('finite_profile',pth); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    la=canonical_laurent(df,window,degree); rows=[]
    use=np.ones(len(df),dtype=bool) if mask is None else np.asarray(mask,dtype=bool)
    for Lv in Lvals:
        fin=mod.canonical_audit(df,float(Lv),window,degree)
        ev=evaluate_blocks(la,1.0/float(Lv))
        row={'L':float(Lv)}
        for nm in ('K','M'):
            f=fin[nm][use]; e=ev[nm][use]
            den=np.maximum(1.0,np.abs(f)); row[f'max_scaled_{nm}']=float(np.max(np.abs(e-f)/den))
            row[f'max_abs_{nm}']=float(np.max(np.abs(e-f)))
        rows.append(row)
    return pd.DataFrame(rows),la


def _poly_add(a,b):
    n=max(len(a),len(b)); out=np.zeros(n,float); out[:len(a)]+=a; out[:len(b)]+=b; return out

def _poly_mul(a,b): return np.convolve(a,b)

def _pl_mul(A:Dict[int,np.ndarray],B:Dict[int,np.ndarray],pmin=-12,pmax=16):
    out={}
    for p,a in A.items():
      for q,b in B.items():
        k=p+q
        if pmin<=k<=pmax:
          v=np.zeros((a.shape[0],a.shape[1]+b.shape[1]-1))
          for i in range(a.shape[1]):
            for j in range(b.shape[1]): v[:,i+j]+=a[:,i]*b[:,j]
          out[k]=out[k]+v if k in out else v
    return out

def _pl_add(A,B,sgn=1.0):
    out={k:np.array(v,copy=True) for k,v in A.items()}
    for k,v in B.items(): out[k]=out[k]+sgn*v if k in out else sgn*np.array(v,copy=True)
    return out


def characteristic_laurent_polynomial(df,la=None,return_scale=False):
    """Return eps-power -> polynomial coefficients [z^0..z^3] for det[(r^2/f) eps M-zK]."""
    if la is None: la=canonical_laurent(df)
    n=len(df); r=df.x.to_numpy(float); f=df.f.to_numpy(float); fac=r*r/f
    entries=[[{} for _ in range(3)] for __ in range(3)]
    for i in range(3):
      for j in range(3):
        D={}
        # eps * (r^2/f) M term, z degree 0
        for p,a in la['M'].coeffs.items():
            poly=np.zeros((n,2)); poly[:,0]=fac*(a[:,i,j]+a[:,j,i])/2
            D[p+1]=D.get(p+1,0)+poly
        # -z K term, z degree 1
        for p,a in la['K'].coeffs.items():
            poly=np.zeros((n,2)); poly[:,1]=-a[:,i,j]
            D[p]=D.get(p,0)+poly
        entries[i][j]=D
    total={}; scale={}
    perms=list(permutations(range(3)))
    for perm in perms:
        inv=sum(perm[a]>perm[b] for a in range(3) for b in range(a+1,3)); s=-1.0 if inv%2 else 1.0
        term={0:np.ones((n,1))}; bound={0:np.ones((n,1))}
        for i,j in enumerate(perm):
            term=_pl_mul(term,entries[i][j])
            bound=_pl_mul(bound,{p:np.abs(v) for p,v in entries[i][j].items()})
        total=_pl_add(total,term,sgn=s)
        scale=_pl_add(scale,bound)
    return (total,scale) if return_scale else total


def leading_characteristic(df,la=None,tol=1e-10):
    P,scale=characteristic_laurent_polynomial(df,la,return_scale=True); n=len(df)
    # choose first eps power with non-negligible polynomial at each row.  The witness uses a uniform power;
    # export per-row powers so limit-sector changes are visible rather than hidden.
    powers=np.full(n,999,dtype=int); coeff=np.zeros((n,4))
    for p in sorted(P):
        arr=P[p]
        if arr.shape[1]<4: arr=np.pad(arr,((0,0),(0,4-arr.shape[1])))
        bound=np.max(scale[p],axis=1)
        # Relative to the absolute determinant expansion before cancellation.
        # This is invariant under a common rescaling of the quadratic action.
        mask=(powers==999)&(np.max(np.abs(arr),axis=1)>tol*np.maximum(bound,np.finfo(float).tiny))
        # Three propagating variables require a resolved cubic kinetic term.
        # Lower-degree roundoff remnants cannot define a three-mode spectrum.
        mask &= np.abs(arr[:,3]) > tol*np.maximum(scale[p][:,3],np.finfo(float).tiny)
        powers[mask]=p; coeff[mask]=arr[mask,:4]
    return powers,coeff,P


def roots_from_coeffs(coeff):
    out=[]
    for c in np.asarray(coeff,float):
        # numpy roots wants descending powers
        if np.max(np.abs(c))==0: out.append([np.nan]*3); continue
        # Preserve small cubic terms: they can encode large physical roots.
        idx=np.flatnonzero(c!=0)
        deg=int(idx.max()) if len(idx) else 0
        rr=np.roots(c[:deg+1][::-1]) if deg else np.array([])
        rreal=[float(x.real) if abs(x.imag)<1e-8 else np.nan for x in rr]
        rreal=sorted(rreal)+[np.nan]*(3-len(rreal)); out.append(rreal[:3])
    return np.asarray(out,float)


def extract_eq83_style_coefficients(df,la):
    """Extract the high-L coefficient set used by the published genuine-SVT angular parametrization.

    Naming follows the project/paper cross-check convention:
      M11_0,M13_0,M33_0 = eps^0 coefficients;
      M11_1,M13_1,M33_1 = eps^1 coefficients;
      M12_0,M23_0 = eps^0 coefficients;
      M22_0 = eps^-1 coefficient; K22_0 = eps^0 coefficient.
    """
    K=la['K']; M=la['M'];
    # finite profile reducer stores project M with the opposite sign to the
    # Zhang-Kase action convention (see canonical_audit comment: '-Y^T M Y').
    # Eq83-style paper coefficients therefore use M_paper = -M_project.
    def mc(p,i,j): return -(M.coeff(p)[:,i,j]+M.coeff(p)[:,j,i])/2
    def kc(p,i,j): return K.coeff(p)[:,i,j]
    return {
      'K22_0':kc(0,1,1),
      'M11_0':mc(0,0,0),'M11_1':mc(1,0,0),
      'M12_0':mc(0,0,1),'M13_0':mc(0,0,2),'M13_1':mc(1,0,2),
      'M22_0':mc(-1,1,1),'M23_0':mc(0,1,2),
      'M33_0':mc(0,2,2),'M33_1':mc(1,2,2),
    }


def published_svt_angular_from_eq83(df,c):
    """Zhang-Kase published angular shortcut evaluated from Eq83-style coefficients.

    This intentionally consumes action-derived coefficients.  A disagreement with the raw determinant
    therefore localizes a convention/mapping issue *before* root sorting.
    """
    r=df.x.to_numpy(float); f=df.f.to_numpy(float); h=df.h.to_numpy(float); ph=df.phiprime.to_numpy(float)
    q={k:df[k].to_numpy(float) for k in SLOTS}
    # alpha7 inversion used consistently in the historical project oracle.
    fp=dr(r,f); hp=dr(r,h); a4=q['a4']; a4p=dr(r,a4); Ap=df.A0prime.to_numpy(float)
    geom=(r*(fp*h+f*hp)-2*f*h)/(2*r*f*h*h*ph)
    alpha7=ph/(2*r*f)*(q['a6']-geom*a4+a4p/(h*ph)+Ap*q['v6']/(4*h*ph))
    M11_0,M11_1=c['M11_0'],c['M11_1']; M13_0,M13_1=c['M13_0'],c['M13_1']; M33_0,M33_1=c['M33_0'],c['M33_1']
    M12_0,M22_0,M23_0=c['M12_0'],c['M22_0'],c['M23_0']; K22=c['K22_0']
    Mc1=(M11_0*M33_1-2*M13_0*M13_1+M33_0*M11_1)/(4*r*r*f*h*h*alpha7*M33_0)
    Mc2=ph*(r*h*ph*M22_0*M33_0+2*M13_0*M23_0-2*M12_0*M33_0)/(4*r*f*h*q['v1']*alpha7*M33_0)
    K1=-(2*r*r*a4/f)/(1-f*q['v6']**2/(8*a4*q['v10']))
    K33=f*q['v1']**2*K1/(2*r*r*a4*q['v10'])
    cV=-r*r*M33_0/(f*K33)
    B1=-(Mc1+Mc2+r*r*M22_0/(f*K22))
    cross=r*r*(4*f*alpha7*Mc2+ph*ph*M22_0)**2/(4*f*f*ph*ph*alpha7*K22)
    B2=r*r*(Mc1+Mc2)*M22_0/(f*K22)-cross
    disc=B1*B1-4*B2
    sq=np.sqrt(np.maximum(disc,0)); cm=.5*(B1-sq); cp=.5*(B1+sq)
    return dict(B1=B1,B2=B2,disc=disc,cV=cV,cminus=cm,cplus=cp,alpha7=alpha7,Mc1=Mc1,Mc2=Mc2)


def action_m5_mass_shortcuts(df,branch='minus'):
    """Action-derived high-L mass shortcuts; plus is retained only as a published-conflict comparator."""
    d=df.copy(); r=d.x.to_numpy(float); h=d.h.to_numpy(float)
    q={k:d[k].to_numpy(float) for k in SLOTS}; Ap=d.A0prime.to_numpy(float); ph=d.phiprime.to_numpy(float)
    a4,a6,v6,v13,v9,v1=q['a4'],q['a6'],q['v6'],q['v13'],q['v9'],q['v1']
    m5=a4*v13 + (a6*v6 if branch=='plus' else -a6*v6)
    m2=2*q['c5']*v1+(Ap*v1+.5*ph*q['v4'])*v6
    m3=a4*dr(r,v1)+.5*Ap*v1*v6; m4=2*Ap*v1+ph*q['v4']
    z1=(a4+r*q['a9']-r*Ap*v6/2)*v6
    m1m=-r*a4*q['v8']+z1
    M11=-r*r*v6*v6/(4*v1); M13=-r*v6/2; M33=-v1
    M22=q['e4']+2*h*q['c4']*a6/a4-m5*m5/(4*a4*a4*v9)
    inner12=2*r*h*q['c4'] + r*(2*q['d2']*v1-q['v4']*v6)/(2*v1) + r*v6*m5/(2*a4*v9)
    M12=(-r*q['d3']/2 - h*q['c4']*z1/(a4*v6) - r*h*a6*m2/(2*a4*v1) + r*q['v5']*v6/(4*v1)
         -m1m*m5/(4*a4*a4*v9)+.25*dr(r,inner12))
    M23=(h*Ap*q['c4']*v1/a4+q['v5']/2+m3*m5/(2*a4*a4*v9)-h*a6*m4/(2*a4)
         -.25*dr(r,q['v4']+v1*m5/(a4*v9)))
    return dict(m5=m5,M11_0=M11,M12_0=M12,M13_0=M13,M22_0=M22,M23_0=M23,M33_0=M33)


def load_genuine_svt_witness():
    # same selected representative used by the accepted finite-L reducer
    pth=B/'ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py'
    spec=importlib.util.spec_from_file_location('finite_profile_load',pth); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod.prepare_central().reset_index(drop=True)
