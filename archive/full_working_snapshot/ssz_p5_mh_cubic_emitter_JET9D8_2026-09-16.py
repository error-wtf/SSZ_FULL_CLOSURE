#!/usr/bin/env python3
"""Simplified Maxwell--Horndeski Appendix-A emitter for the selected P5 outer cubic member.

Selected member:
  G4 = 1/2, G4_X = 0, G5 = 0,
  G2_F = 1, G2_FF = G2_XF = G2_phiF = 0,
  c_r,S^2 = 1 fixes c2 (equivalently one G2_XX transverse jet),
  c3=e3=0 lower-order representative.

Canonical 41-slot mapping embeds Kase--Tsujikawa MH v1..v9 into the common
Zhang--Kase 13-vector schema:
 old v1..v5 -> v1..v5; old v6->v7; old v7->v8; old v8->v9; old v9->v10;
 common v6=v11=v12=v13=0.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd

B=Path('/mnt/data')

def deriv(x,y,order=1,window=9,degree=8):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(x); out=np.empty(n)
    half=window//2
    for i in range(n):
        lo=max(0,i-half); hi=min(n,lo+window); lo=max(0,hi-window)
        xx=x[lo:hi]; yy=y[lo:hi]; x0=x[i]
        sc=float(np.max(np.abs(xx-x0))) or 1.0
        z=(xx-x0)/sc
        deg=min(degree,len(xx)-1)
        co=np.polynomial.polynomial.polyfit(z,yy,deg)
        out[i]=math.factorial(order)*co[order]/sc**order
    return out

def emit(df, g3x=None, target_cr2=1.0, c2_profile=None, window=9, degree=8):
    d=df.copy().reset_index(drop=True)
    r=d['x'].to_numpy(float); f=d['f'].to_numpy(float); h=d['h'].to_numpy(float)
    ph=(d['phi_r'] if 'phi_r' in d else (d['phiprime'] if 'phiprime' in d else pd.Series(-np.sqrt(np.maximum(0,-2*d['X'].to_numpy(float)/d['h'].to_numpy(float)))))).to_numpy(float)
    Ap=(d['A0prime'].to_numpy(float) if 'A0prime' in d else np.zeros(len(d)))
    if g3x is None:
        if 'G3X' not in d: raise ValueError('Need G3X or g3x')
        g3x=d['G3X'].to_numpy(float)
    else: g3x=np.asarray(g3x,float)

    fp=deriv(r,f,1,window,degree); hp=deriv(r,h,1,window,degree)
    fpp=deriv(r,f,2,window,degree); phipp=deriv(r,ph,1,window,degree)
    sqfh=np.sqrt(f*h)
    # Selected constant-G4 cubic Horndeski member: F=G=H=1.
    Ft=np.ones_like(r); Gt=np.ones_like(r); Ht=np.ones_like(r)
    a1=0.5*sqfh*h*g3x*ph**2*r**2
    a4=0.5*sqfh
    a1p=deriv(r,a1,1,window,degree); a4p=deriv(r,a4,1,window,degree)
    a2=sqfh*deriv(r,a1/sqfh,1,window,degree) - (phipp/ph-0.5*fp/f)*a1 + r/ph*(fp/f-hp/h)*a4
    a3=-0.5*ph*a1-r*a4
    a3p=deriv(r,a3,1,window,degree)
    # Maxwell baseline is kept unpartitioned: G2_F=1, no mixed F derivatives.
    old_v1=0.5*r**2*np.sqrt(h/f)
    old_v2=Ap*old_v1
    old_v4=np.zeros_like(r); old_v5=np.zeros_like(r)
    old_v3=-Ap*old_v1
    old_v6=0.25*Ap**2*old_v1
    old_v8=1/(2*sqfh)
    old_v7=-2*h*Ap*old_v8
    old_v9=-f*h*old_v8

    P=phipp/ph+0.5*hp/h
    R=r/ph*(fp/f-hp/h)*a4
    a5=-deriv(r,P,1,window,degree)*a1-P*a1p+deriv(r,R,1,window,degree) # exact holonomic path, v4=v5=0
    a6=np.zeros_like(r)
    a7=a3p-0.5*Ap**2*old_v1
    a8=-0.5*a4/h
    a9=a4p+(1/r-0.5*fp/f)*a4
    b1=a4/(2*f); b2=-2*a1/f; b3=-2*(a2-a1p)/f; b4=-2*a3/f; b5=-2*b1
    c1=-a1/(f*h)

    # Eq. (4.30): P1 and mu, H=1.
    mu=2*(ph*a1+2*r*a4)/sqfh
    Y=f*r**4/(mu**2*h)
    P1=h*mu/(2*f*r**2)*deriv(r,Y,1,window,degree)
    denom=(f*h)**1.5*(2*P1-1.0)*mu**2
    pref=4*ph/denom
    c4=a1/(h*r**2) # exact identity for G4=1/2,G5=0 cubic member
    B0=(8*r**2*h*a4*c4*(ph*a1+r*a4)
        -sqfh*ph*a1**2
        +2*r**2*a4**2*(fp/f*a1 + ph*old_v2**2/(2*old_v1)))
    c2_unit=(target_cr2/pref-B0)/(4*r**2*a4**2)
    if c2_profile is None:
        c2=c2_unit
    elif np.isscalar(c2_profile):
        c2=np.full_like(r,float(c2_profile))
    else:
        c2=np.asarray(c2_profile,float)

    c3=np.zeros_like(r)
    c5=-h*ph*c4-0.5*sqfh/r-0.5*fp/f*a4
    c6=(fp*ph/(8*f))*a1+(fp/(2*r*f))*a4-0.25*ph*c2+0.5*h*ph*r*c4+0.25*sqfh+0.25*Ap**2*old_v1
    d1=a4/(2*f); d2=2*h*c4
    # Appendix A d3. ∂a4/∂phi=0 for constant G4; F=G=1.
    photon=fp*r-2*f
    Q=(2*phipp/(h*ph*r)+fp**2/f**2-fp*hp/(f*h)-2*fp/(f*r)+2*hp/(h*r)+hp/(h**2*r))
    with np.errstate(divide='ignore',invalid='ignore'):
        d3=(-((2*phipp/ph+hp/h)*a1)/r**2
            +2*f/(photon*ph)*Q*a4
            +np.sqrt(f)/(ph*np.sqrt(h)*r**2)
            -f**1.5/(np.sqrt(h)*photon*ph)*(fp/(f*r)+2*phipp/(ph*r)+hp/(h*r)-2/r**2))
    d4=0.5*sqfh/r**2
    e1=((fp/f+0.5*hp/h)*a1-2*a1p+a2)/(ph*f*h)
    e2=-(fp/f*a1+2*c2+4*h*r*c4)/(2*ph)
    e3=np.zeros_like(r)
    c4p=deriv(r,c4,1,window,degree)
    # Appendix A e4 with F=G=H=1, G'=0 and ∂a4/∂phi=0.
    with np.errstate(divide='ignore',invalid='ignore'):
        e4=(c4p/ph
            -0.5*fp*a4p/(f*ph**2*h)
            +(phipp/ph+0.5*hp/h)*a1/(h*ph*r**2)
            +a4/(4*h*ph**2)*((photon-4*f)*fp/(f**2*r) + hp*(photon+6*f)/(h*r*f)
                -4*f*(2*phipp*h+hp*ph)/(ph*h**2*r*photon))
            +0.5*hp*c4/(h*ph)
            +0.5*(fp*h*r-f)/(r**2*np.sqrt(f)*ph**2*h**1.5)
            +0.5*np.sqrt(f)/(r*ph**2*h**1.5)*(f*(2*phipp*h+hp*ph)/(h*ph*photon)
                +0.5*(2*f-fp*h*r)/(f*r)))

    out=pd.DataFrame({
      'u':d['u'].to_numpy(float) if 'u' in d else 1/r,'x':r,
      'phi':d['phi'].to_numpy(float) if 'phi' in d else np.nan,
      'f':f,'h':h,'phiprime':ph,'A0prime':Ap,
      'a1':a1,'a2':a2,'a3':a3,'a4':a4,'a5':a5,'a6':a6,'a7':a7,'a8':a8,'a9':a9,
      'b1':b1,'b2':b2,'b3':b3,'b4':b4,'b5':b5,
      'c1':c1,'c2':c2,'c3':c3,'c4':c4,'c5':c5,'c6':c6,
      'd1':d1,'d2':d2,'d3':d3,'d4':d4,
      'e1':e1,'e2':e2,'e3':e3,'e4':e4,
      'v1':old_v1,'v2':old_v2,'v3':old_v3,'v4':old_v4,'v5':old_v5,
      'v6':np.zeros_like(r),'v7':old_v6,'v8':old_v7,'v9':old_v8,'v10':old_v9,
      'v11':np.zeros_like(r),'v12':np.zeros_like(r),'v13':np.zeros_like(r),
      'mu':mu,'P1':P1,'c2_unit_reference':c2_unit,'selected_crS2':np.full_like(r,target_cr2),'photon_factor':photon,
      'G3X_selected':g3x
    })
    return out

SLOTS=[*(f'a{i}' for i in range(1,10)),*(f'b{i}' for i in range(1,6)),*(f'c{i}' for i in range(1,7)),*(f'd{i}' for i in range(1,5)),*(f'e{i}' for i in range(1,5)),*(f'v{i}' for i in range(1,14))]

def scaled_rel(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return np.abs(a-b)/np.maximum(1.0,np.abs(b))

if __name__=='__main__':
    strong=pd.read_csv(B/'ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv').sort_values('x').reset_index(drop=True)
    r=strong.x.to_numpy(float); f=strong.f.to_numpy(float); h=strong.h.to_numpy(float); ph=strong.phi_r.to_numpy(float)
    g3x=2*strong.a1.to_numpy(float)/(np.sqrt(f*h)*h*ph**2*r**2)
    reb=emit(strong,g3x=g3x)
    rows=[]
    mask=np.abs(reb.photon_factor.to_numpy(float))>2e-3
    for c in SLOTS:
        if c in ['c3','e3'] or c not in strong or strong[c].isna().all(): continue
        sr=scaled_rel(reb[c],strong[c]); use=sr[np.isfinite(sr)&mask]
        rows.append(dict(slot=c,median=float(np.median(use)),p95=float(np.quantile(use,.95)),max=float(np.max(use))))
    reg=pd.DataFrame(rows)
    reb.to_csv(B/'ssz_p5_MH_CUBIC_EMITTER_STRONGH_REBUILT_2026-09-16.csv',index=False)
    reg.to_csv(B/'ssz_p5_MH_CUBIC_EMITTER_STRONGH_REGRESSION_2026-09-16.csv',index=False)
    print(reg.sort_values('max',ascending=False).head(20).to_string(index=False))
