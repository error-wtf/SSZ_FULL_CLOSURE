#!/usr/bin/env python3
"""Profile-aware SSZ P5 even-parity differential-operator reducer.

Builds the reduced Euler-Lagrange differential operator directly from the
unreduced 41-slot Zhang-Kase action after the accepted generalized-psi
constraint maps.  All radial product rules are retained with JET9D8.

This is an audit/reconstruction layer, not by itself a global QNM claim.
"""
from __future__ import annotations
from pathlib import Path
from math import comb
import importlib.util
import numpy as np
import pandas as pd

from ssz_p5.paths import paths as B

# accepted derivative service
spec=importlib.util.spec_from_file_location('hj',B/'ssz_p5_higher_jet_closure_2026-09-16.py')
hj=importlib.util.module_from_spec(spec); spec.loader.exec_module(hj)
# accepted full algebraic maps (uploaded name may carry (1))
FM=B/'ssz_hybrid_full_constraint_maps_JET9D8(1).py'
if not FM.exists(): FM=B/'ssz_hybrid_full_constraint_maps_JET9D8.py'
spec=importlib.util.spec_from_file_location('fm',FM)
fm=importlib.util.module_from_spec(spec); spec.loader.exec_module(fm)

SLOTS=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],
       *[f'c{i}' for i in range(1,7)],*[f'd{i}' for i in range(1,5)],
       *[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]


def deriv(r, a, order, window=9, degree=8):
    from ssz_p5.jets.jet9d8 import profile_derivative
    return profile_derivative(r, a, order, window, degree)


def addop(dst,key,val):
    if key in dst: dst[key] += val
    else: dst[key] = np.array(val,float,copy=True)


def field_maps(df,L,window=9,degree=8):
    """Return field -> scalar row differential operator on y=(psi,dphi,V)."""
    m=fm.maps(df,L,window,degree)
    n=len(df); z=np.zeros((n,3)); ephi=z.copy(); ephi[:,1]=1.; eV=z.copy(); eV[:,2]=1.
    return {
      'dphi': {(0,0):ephi},
      'V': {(0,0):eV},
      'h1': {(0,0):m['H10'],(0,1):m['H11']},
      'H2': {(0,0):m['H20'],(0,1):m['H21']},
      'H1': {(1,0):m['H1dot0'],(1,1):m['H1dot1']},
      'dA1':{(1,0):m['dA1dot0'],(1,1):m['dA1dot1']},
      'dA0':{(0,0):m['dA00'],(0,1):m['dA01']},
    },m


def apply_derivative(op,dt,dr,r,window=9,degree=8):
    """Apply D_t^dt D_r^dr to a scalar-row differential operator."""
    out={}
    for (pt,pr),a in op.items():
        for k in range(dr+1):
            ak=deriv(r,a,k,window,degree) if k else a
            key=(pt+dt, pr+dr-k)
            addop(out,key,comb(dr,k)*ak)
    return out


def action_terms(df,L):
    """Exact non-H0 terms of the common unreduced even action.

    Every tuple is (coefficient array, field A, deriv A, field B, deriv B).
    H0-proportional terms are absent because the generalized-psi map enforces
    the H0 multiplier constraint identically; the H0^2 coefficient is audited
    separately before this reduction is accepted.
    """
    q={k:df[k].to_numpy(float) for k in SLOTS}; n=len(df)
    T=[]
    def add(c,A,da,B,db): T.append((np.asarray(c,float),A,da,B,db))
    # metric/scalar non-H0 block
    add(L*q['b1'],'H1',(0,0),'H1',(0,0))
    add(q['b2'],'H1',(0,0),'dphi',(1,1))
    add(q['b3'],'H1',(0,0),'dphi',(1,0))
    add(q['b4'],'H1',(0,0),'H2',(1,0))
    add(L*q['b5'],'H1',(0,0),'h1',(1,0))
    add(q['c1'],'dphi',(1,0),'H2',(1,0))
    add(q['c2'],'H2',(0,0),'dphi',(0,1))
    add(q['c3']+L*q['c4'],'H2',(0,0),'dphi',(0,0))
    add(L*q['c5'],'H2',(0,0),'h1',(0,0))
    add(q['c6'],'H2',(0,0),'H2',(0,0))
    add(L*q['d1'],'h1',(1,0),'h1',(1,0))
    add(L*q['d2'],'h1',(0,0),'dphi',(0,1))
    add(L*q['d3'],'h1',(0,0),'dphi',(0,0))
    add(L*q['d4'],'h1',(0,0),'h1',(0,0))
    add(q['e1'],'dphi',(1,0),'dphi',(1,0))
    add(q['e2'],'dphi',(0,1),'dphi',(0,1))
    add(q['e3']+L*q['e4'],'dphi',(0,0),'dphi',(0,0))

    # auxiliary-vector form, H0 set to zero after its exact constraint
    add(2*q['v1'],'V',(0,0),'dA0',(0,1))
    add(-2*q['v1'],'V',(0,0),'dA1',(1,0))
    add(-q['v1'],'V',(0,0),'V',(0,0))
    # V*S0 and -S0^2/(4 v1), S0 excludes v2 H0
    parts=[('H2',q['v3']),('dphi_r',q['v4']),('dphi',q['v5']),('h1',L*q['v6'])]
    pmap={'H2':('H2',(0,0)),'dphi_r':('dphi',(0,1)),'dphi':('dphi',(0,0)),'h1':('h1',(0,0))}
    for name,c in parts:
        f,dd=pmap[name]; add(c,'V',(0,0),f,dd)
    for ia,(na,ca) in enumerate(parts):
        fa,da=pmap[na]
        for ib in range(ia,len(parts)):
            nb,cb=parts[ib]; fb,db=pmap[nb]
            fac=2.0 if ib!=ia else 1.0
            add(-fac*ca*cb/(4*q['v1']),fa,da,fb,db)
    add(0.5*L*q['v6'],'h1',(0,0),'dA1',(1,0))
    add(L*q['v8'],'h1',(0,0),'dA0',(0,0))
    add(L*q['v9'],'dA0',(0,0),'dA0',(0,0))
    add(L*q['v10'],'dA1',(0,0),'dA1',(0,0))
    add(L*q['v11'],'H1',(0,0),'dA1',(0,0))
    add(L*q['v12'],'H2',(0,0),'dA0',(0,0))
    add(L*q['v13'],'dphi',(0,0),'dA0',(0,0))
    return T


def adjoint_product(A,c,Bop,r,window=9,degree=8):
    """A^dagger * c * B, for scalar-row operators A,B; returns 3x3 operator."""
    n=len(r); out={}
    for (pt,pr),av in A.items():
      for (qt,qr),bv in Bop.items():
        base=av[:,:,None]*bv[:,None,:]*c[:,None,None]  # output i, input j
        sgn=(-1.0)**(pt+pr)
        for k in range(pr+1):
            bk=deriv(r,base,k,window,degree) if k else base
            key=(pt+qt, pr-k+qr)
            addop(out,key,sgn*comb(pr,k)*bk)
    return out


def reduced_operator(df,L,window=9,degree=8):
    d=df.copy().reset_index(drop=True)
    # accepted central Appendix-A branch
    d['v12']=-d['v6']/(2*d['h'])
    r=d.x.to_numpy(float)
    F,m=field_maps(d,L,window,degree)
    P={}
    for c,A,da,Bn,db in action_terms(d,L):
        Aop=apply_derivative(F[A],da[0],da[1],r,window,degree)
        Bop=apply_derivative(F[Bn],db[0],db[1],r,window,degree)
        X=adjoint_product(Aop,c,Bop,r,window,degree)
        Y=adjoint_product(Bop,c,Aop,r,window,degree)
        for k,v in X.items(): addop(P,k,v)
        for k,v in Y.items(): addop(P,k,v)
    # h0 quadratic identity (should vanish before setting H0=0)
    h0_quad=d['v7'].to_numpy(float)-d['v2'].to_numpy(float)**2/(4*d['v1'].to_numpy(float))
    return P,m,h0_quad


def canonical_audit(df,L,window=9,degree=8):
    P,m,h0q=reduced_operator(df,L,window,degree)
    r=df.x.to_numpy(float); n=len(df)
    Z=lambda: np.zeros((n,3,3))
    get=lambda key: P.get(key,Z())
    # Euler-operator principal blocks. With target action normalization,
    # K=-P20/2, G=P02/2 if source action is written without overall 1/2.
    # Both raw P and half-normalized blocks are exported for regression.
    P20=get((2,0)); P11=get((1,1)); P10=get((1,0)); P02=get((0,2)); P01=get((0,1)); P00=get((0,0))
    K=-0.5*P20; G=0.5*P02; R=-0.5*P11
    # Zhang--Kase Eq. (4.25) contains an additional antisymmetric radial
    # connection block S:  Y'^T S Y, S^T=-S.  In the project sign convention
    # we keep -Y'^T G Y' but retain +Y'^T S Y.  For the Euler operator,
    #   P01 = P02' - 2 S,   P00 = -2 M - S'
    # where M below is the project convention (minus sign in the action).
    P11p=np.empty_like(P11); P02p=np.empty_like(P02)
    for i in range(3):
      for j in range(3):
        P11p[:,i,j]=deriv(r,P11[:,i,j],1,window,degree)
        P02p[:,i,j]=deriv(r,P02[:,i,j],1,window,degree)
    S=0.5*(P02p-P01)
    Sp=np.empty_like(S)
    for i in range(3):
      for j in range(3):
        Sp[:,i,j]=deriv(r,S[:,i,j],1,window,degree)
    M=-0.5*(P00+Sp)
    diagnostics={
      'max_h0_quadratic':float(np.max(np.abs(h0q))),
      'max_P20_asym':float(np.max(np.abs(P20-np.swapaxes(P20,1,2)))),
      'max_P02_asym':float(np.max(np.abs(P02-np.swapaxes(P02,1,2)))),
      'max_P11_asym':float(np.max(np.abs(P11-np.swapaxes(P11,1,2)))),
      'max_time_first_sym_res':float(np.max(np.abs(0.5*(P10+np.swapaxes(P10,1,2))-0.5*P11p))),
      'max_radial_first_sym_res':float(np.max(np.abs(0.5*(P01+np.swapaxes(P01,1,2))-P02p))),
    }
    high=[k for k in P if k[0]>2 or k[1]>2]
    diagnostics['higher_operator_keys']=str(sorted(high))
    diagnostics['max_high_order']=float(max([np.max(np.abs(P[k])) for k in high],default=0.0))
    diagnostics['max_S_sym']=float(np.max(np.abs(S+np.swapaxes(S,1,2))))
    diagnostics['max_R_abs']=float(np.max(np.abs(R)))
    return dict(P=P,K=K,R=R,G=G,S=S,M=M,maps=m,diagnostics=diagnostics)


def prepare_central():
    d=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv').copy()
    ce=pd.read_csv(B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv')
    d['c3']=ce.c3_selected.to_numpy(float); d['e3']=ce.e3_selected.to_numpy(float)
    d['v12']=-d['v6']/(2*d['h'])
    return d

if __name__=='__main__':
    d=prepare_central(); rows=[]
    for L in (6,12,20,42,110,420,1000):
        a=canonical_audit(d,float(L)); row={'L':L,**a['diagnostics']}
        for nm in ('K','G','S','M'):
            A=a[nm]; row[f'{nm}_maxabs']=float(np.max(np.abs(A)))
        rows.append(row)
    out=pd.DataFrame(rows)
    out.to_csv(B/'ssz_p5_profile_operator_structural_audit_2026-09-16.csv',index=False)
    print(out.to_string(index=False))
