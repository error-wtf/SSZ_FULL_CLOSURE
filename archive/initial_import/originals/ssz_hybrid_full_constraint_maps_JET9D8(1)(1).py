#!/usr/bin/env python3
"""
SSZ P5 full algebraic constraint-map audit (JET9D8), 2026-09-16.

Implements the exact manuscript Sec. 20.4 constraint architecture after the
generalized-psi H0 reduction, using the accepted local polynomial radial jet
service rather than CubicSpline differentiation.

Physical basis: y = (psi, dphi, V), with y_r = (psi', dphi', V').
Maps exported:
  h1  = H10*y + H11*y_r
  H2  = H20*y + H21*y_r
  H1  = T10*dot(y) + T11*dot(y_r)
  dA1 = A10*dot(y) + A11*dot(y_r)
  dA0 = A00*y + A01*y_r
"""
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd

B = Path('/mnt/data')
SLOTS = [*[f'a{i}' for i in range(1,10)],
         *[f'b{i}' for i in range(1,6)],
         *[f'c{i}' for i in range(1,7)],
         *[f'd{i}' for i in range(1,5)],
         *[f'e{i}' for i in range(1,5)],
         *[f'v{i}' for i in range(1,14)]]

spec=importlib.util.spec_from_file_location('hj',B/'ssz_p5_higher_jet_closure_2026-09-16.py')
hj=importlib.util.module_from_spec(spec); spec.loader.exec_module(hj)

def dr(r,y,window=9,degree=8):
    return hj.local_poly_deriv(np.asarray(r,float),np.asarray(y,float),1,window,degree)

def maps(df,L,window=9,degree=8):
    d=df.copy().reset_index(drop=True)
    miss=[k for k in ['x']+SLOTS if k not in d.columns]
    if miss: raise KeyError(miss)
    r=d.x.to_numpy(float); n=len(d)
    q={k:d[k].to_numpy(float) for k in SLOTS}

    p=q['a1']/q['a3']
    qq=L*q['a4']/q['a3']
    pp=dr(r,p,window,degree); qrp=dr(r,qq,window,degree)

    Braw=q['a2']-q['v2']*q['v4']/(2*q['v1'])
    Cphi=q['a5']+L*q['a6']-q['v2']*q['v5']/(2*q['v1'])
    CH2=q['a7']+L*q['a8']-q['v2']*q['v3']/(2*q['v1'])
    Ch1=L*(q['a9']-q['v2']*q['v6']/(2*q['v1']))
    Beff=Braw-q['a3']*pp-CH2*p
    Dh1=Ch1-q['a3']*qrp-CH2*qq

    # h1 = H10 y + H11 y'
    H10=np.zeros((n,3)); H11=np.zeros((n,3))
    H10[:,0]=-CH2/Dh1
    H10[:,1]=-Cphi/Dh1
    H10[:,2]=-q['v2']/Dh1
    H11[:,0]=-q['a3']/Dh1
    H11[:,1]=-Beff/Dh1

    # H2 = psi - q h1 - p dphi'
    H20=np.zeros((n,3)); H21=np.zeros((n,3))
    H20[:,0]=1.0
    H20-=qq[:,None]*H10
    H21-=qq[:,None]*H11
    H21[:,1]-=p

    # J1 = b3 dot dphi + b4 dot H2 + L b5 dot h1.
    J10=np.zeros((n,3)); J11=np.zeros((n,3))
    J10[:,1]+=q['b3']
    J10+=q['b4'][:,None]*H20 + (L*q['b5'])[:,None]*H10
    J11+=q['b4'][:,None]*H21 + (L*q['b5'])[:,None]*H11

    # R1=-J1/L ; R2=-(2v1/L)dot V +(v6/2)dot h1.
    R10=-J10/L; R11=-J11/L
    R20=(q['v6']/2)[:,None]*H10
    R21=(q['v6']/2)[:,None]*H11
    R20[:,2]+=-2*q['v1']/L

    Delta=4*q['b1']*q['v10']-q['v11']**2

    # H1=(2v10 R1-v11 R2)/Delta
    T10=((2*q['v10'])[:,None]*R10-q['v11'][:,None]*R20)/Delta[:,None]
    T11=((2*q['v10'])[:,None]*R11-q['v11'][:,None]*R21)/Delta[:,None]

    # dA1=(-v11 R1+2b1 R2)/Delta
    A10=(-q['v11'][:,None]*R10+(2*q['b1'])[:,None]*R20)/Delta[:,None]
    A11=(-q['v11'][:,None]*R11+(2*q['b1'])[:,None]*R21)/Delta[:,None]

    # dA0=(v1 V)'/(L v9) -(v8 h1+v12 H2+v13 dphi)/(2v9)
    v1p=dr(r,q['v1'],window,degree)
    A00=-(q['v8'][:,None]*H10+q['v12'][:,None]*H20)/(2*q['v9'])[:,None]
    A01=-(q['v8'][:,None]*H11+q['v12'][:,None]*H21)/(2*q['v9'])[:,None]
    A00[:,1]+=-q['v13']/(2*q['v9'])
    A00[:,2]+=v1p/(L*q['v9'])
    A01[:,2]+=q['v1']/(L*q['v9'])

    return dict(
        p=p,q=qq,Dh1=Dh1,DeltaV=Delta,pivotA0=2*L*q['v9'],
        H10=H10,H11=H11,H20=H20,H21=H21,
        H1dot0=T10,H1dot1=T11,dA1dot0=A10,dA1dot1=A11,
        dA00=A00,dA01=A01,
        cancel_dphi2=q['a1']-q['a3']*p,
        cancel_h1prime=L*q['a4']-q['a3']*qq
    )

def audit(df,label,Ls=(6,12,20,30,42,72,110,210,420,1000)):
    rows=[]
    for L in Ls:
        m=maps(df,float(L))
        row={
          'sector':label,'L':L,
          'min_abs_Dh1':float(np.min(np.abs(m['Dh1']))),
          'min_abs_DeltaV':float(np.min(np.abs(m['DeltaV']))),
          'min_abs_2Lv9':float(np.min(np.abs(m['pivotA0']))),
          'max_abs_dphi2_cancel':float(np.max(np.abs(m['cancel_dphi2']))),
          'max_abs_h1prime_cancel':float(np.max(np.abs(m['cancel_h1prime']))),
          'max_abs_h1_map':float(max(np.max(np.abs(m['H10'])),np.max(np.abs(m['H11'])))),
          'max_abs_H1_map':float(max(np.max(np.abs(m['H1dot0'])),np.max(np.abs(m['H1dot1'])))),
          'max_abs_dA1_map':float(max(np.max(np.abs(m['dA1dot0'])),np.max(np.abs(m['dA1dot1'])))),
          'max_abs_dA0_map':float(max(np.max(np.abs(m['dA00'])),np.max(np.abs(m['dA01'])))),
        }
        rows.append(row)
    return pd.DataFrame(rows)

if __name__=='__main__':
    central=pd.read_csv(B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv')
    # Force the Appendix-A branch that is independently validated by the higher-jet audit.
    central['v12']=-central['v6']/(2*central['h'])
    out=audit(central,'central_genuine_SVT')
    out.to_csv(B/'ssz_p5_full_constraint_maps_JET9D8_audit_2026-09-16.csv',index=False)
    print(out.to_string(index=False))
