#!/usr/bin/env python3
"""Maxwell--Horndeski Appendix-A emitter for the P5 luminal G4(phi), G5=0 carrier.

Specialization of Kase & Tsujikawa (2023), Appendix A, to
  G4 = G4(phi), G4_X = 0, G5 = 0,
  G2_F = 1, G2_FF = G2_XF = G2_phiF = 0,
with the P5 carrier principal data H=F=G=2 G4 and finite G3_X.

The scalar principal slot c2 is accepted as an explicit target; in the covariant
action it is realized by the remaining background-null G2_XX/G3 normal-jet freedom.
Lower-order c3,e3 are set to the selected representative zero.  a5 uses the
correct Appendix-A identity including -a1'' via a holonomic first-derivative form.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd

from ssz_p5.paths import paths as B

def deriv(x, y, order=1, window=9, degree=8):
    from ssz_p5.jets.jet9d8 import derivative
    return derivative(x, y, order, window=window, degree=degree)


def emit(df, H=None, G4phi=None, G3X=None, c2_profile=None, window=9, degree=8):
    d=df.copy().reset_index(drop=True)
    # This emitter is explicitly G4(phi), G5=0. Never silently discard
    # nonzero quartic/quintic action jets supplied by a core action profile.
    from ssz_p5.policy import numerical_policy
    tolerance = numerical_policy()["exact_zero_abs"]
    for name in ("G4X", "G4XX", "G4phiX", "G5X", "G5phi", "G5phiX", "G5phiphi"):
        if name in d:
            values = d[name].to_numpy(float)
            if not np.isfinite(values).all() or np.max(np.abs(values)) > tolerance:
                i = int(np.argmax(np.abs(values)))
                raise ValueError(
                    f"luminal G4(phi), G5=0 emitter cannot emit nonzero {name}; "
                    f"row={i}, value={values[i]:.17g}"
                )
    r=(d['x'] if 'x' in d else d['r_over_rs']).to_numpy(float)
    f=d['f'].to_numpy(float); h=d['h'].to_numpy(float)
    if 'phi_r' in d: ph=d['phi_r'].to_numpy(float)
    elif 'phiprime' in d: ph=d['phiprime'].to_numpy(float)
    elif 'X' in d: ph=-np.sqrt(np.maximum(0,-2*d['X'].to_numpy(float)/h))
    else: raise ValueError('need phi_r/phiprime/X')
    Ap=d['A0prime'].to_numpy(float) if 'A0prime' in d else np.zeros(len(d))

    if H is None:
        if 'H_equals_F_equals_G' in d: H=d['H_equals_F_equals_G'].to_numpy(float)
        elif 'H_tensor' in d: H=d['H_tensor'].to_numpy(float)
        else: raise ValueError('need H')
    else: H=np.asarray(H,float)
    F=H.copy(); G=H.copy()

    if G4phi is None:
        if 'G4_phi' in d: G4phi=d['G4_phi'].to_numpy(float)
        else: G4phi=0.5*deriv(r,H,1,window,degree)/ph
    else: G4phi=np.asarray(G4phi,float)
    if G3X is None:
        if 'G3X' not in d: raise ValueError('need G3X')
        G3X=d['G3X'].to_numpy(float)
    else: G3X=np.asarray(G3X,float)

    fp=deriv(r,f,1,window,degree); hp=deriv(r,h,1,window,degree)
    phipp=deriv(r,ph,1,window,degree); Hp=deriv(r,H,1,window,degree)
    sqfh=np.sqrt(f*h)

    # Appendix A a1..a4 for G4(phi), G5=0.
    a1=sqfh*(G4phi+0.5*h*G3X*ph**2)*r**2
    a4=0.5*sqfh*H
    a1p=deriv(r,a1,1,window,degree); a4p=deriv(r,a4,1,window,degree)

    # Maxwell baseline, canonical old MH vector coefficients.
    ov1=0.5*r**2*np.sqrt(h/f)
    ov2=Ap*ov1; ov3=-Ap*ov1
    ov4=np.zeros_like(r); ov5=np.zeros_like(r)
    ov6=0.25*Ap**2*ov1
    ov8=1/(2*sqfh); ov7=-2*h*Ap*ov8; ov9=-f*h*ov8

    a2=sqfh*deriv(r,a1/sqfh,1,window,degree) -(phipp/ph-0.5*fp/f)*a1 + r/ph*(fp/f-hp/h)*a4
    a3=-0.5*ph*a1-r*a4
    a3p=deriv(r,a3,1,window,degree)
    a6=-np.sqrt(f)/(2*np.sqrt(h)*ph)*(Hp + H/r - F/r)  # = -sqrt(f/h) H'/(2 phi') here

    # Correct a5 without differentiating rounded a1 twice.
    # From a2 = a1' - P a1 + R (v4=0), a5=a2'-a1''=-P'a1-Pa1'+R'.
    P=phipp/ph+0.5*hp/h
    R=r/ph*(fp/f-hp/h)*a4
    a5=-deriv(r,P,1,window,degree)*a1-P*a1p+deriv(r,R,1,window,degree)

    a7=a3p-0.5*Ap**2*ov1
    a8=-0.5*a4/h
    a9=a4p+(1/r-0.5*fp/f)*a4
    b1=a4/(2*f); b2=-2*a1/f; b3=-2*(a2-a1p)/f; b4=-2*a3/f; b5=-2*b1
    c1=-a1/(f*h)

    if c2_profile is None:
        if 'c2' in d: c2=d['c2'].to_numpy(float)
        else: raise ValueError('need c2_profile/column c2')
    elif np.isscalar(c2_profile): c2=np.full_like(r,float(c2_profile))
    else: c2=np.asarray(c2_profile,float)
    c3=np.zeros_like(r)

    # Appendix A c4; algebraically equals a1/(h r^2) in this subclass.
    c4=np.sqrt(f)/(4*np.sqrt(h))*(4*G4phi+2*h*G3X*ph**2)
    c5=-h*ph*c4-0.5*sqfh/r*G-0.5*fp/f*a4
    c6=(fp*ph/(8*f))*a1+(fp/(2*r*f))*a4-0.25*ph*c2+0.5*h*ph*r*c4+0.25*sqfh*G+0.25*Ap**2*ov1
    d1=a4/(2*f); d2=2*h*c4

    photon=fp*r-2*f
    da4_dphi=sqfh*G4phi  # partial_phi a4 = sqrt(fh)/2 * partial_phi H; H=2G4
    Q=(2*phipp/(h*ph*r)+fp**2/f**2-fp*hp/(f*h)-2*fp/(f*r)+2*hp/(h*r)+hp/(h**2*r))
    with np.errstate(divide='ignore',invalid='ignore'):
        d3=(-((2*phipp/ph+hp/h)*a1)/r**2
            +2*f/(photon*ph)*Q*a4
            +(photon/(f*r))*da4_dphi
            +np.sqrt(f)/(ph*np.sqrt(h)*r**2)*F
            -f**1.5/(np.sqrt(h)*photon*ph)*(fp/(f*r)+2*phipp/(ph*r)+hp/(h*r)-2/r**2)*G)
    d4=0.5*sqfh/r**2*G
    e1=((fp/f+0.5*hp/h)*a1-2*a1p+a2-2*r*h*a6)/(ph*f*h)
    e2=-(fp/f*a1+2*c2+4*h*r*c4)/(2*ph)
    e3=np.zeros_like(r)
    c4p=deriv(r,c4,1,window,degree); Gp=deriv(r,G,1,window,degree)
    with np.errstate(divide='ignore',invalid='ignore'):
        e4=(c4p/ph
            -0.5*fp*a4p/(f*ph**2*h)
            -0.5*np.sqrt(f)*Gp/(ph**2*np.sqrt(h)*r)
            +(phipp/ph+0.5*hp/h)*a1/(h*ph*r**2)
            +a4/(4*h*ph**2)*((photon-4*f)*fp/(f**2*r) + hp*(photon+6*f)/(h*r*f)
                -4*f*(2*phipp*h+hp*ph)/(ph*h**2*r*photon))
            +0.5*hp*c4/(h*ph)
            -0.5*photon*da4_dphi/(f*h*r*ph)
            +0.5*(fp*h*r-f)*F/(r**2*np.sqrt(f)*ph**2*h**1.5)
            +0.5*np.sqrt(f)*G/(r*ph**2*h**1.5)*(f*(2*phipp*h+hp*ph)/(h*ph*photon)
                +0.5*(2*f-fp*h*r)/(f*r)))

    # Correct P1 Eq. (4.30), retaining H powers.
    mu=2*(ph*a1+2*r*a4)/sqfh
    Y=f*r**4*H**4/(mu**2*h)
    P1=h*mu/(2*f*r**2*H**2)*deriv(r,Y,1,window,degree)
    K=2*P1-F

    out=pd.DataFrame({
      'u':d['u'].to_numpy(float) if 'u' in d else 1/r,'x':r,
      'phi':d['phi'].to_numpy(float) if 'phi' in d else np.nan,
      'f':f,'h':h,'phiprime':ph,'A0prime':Ap,
      'a1':a1,'a2':a2,'a3':a3,'a4':a4,'a5':a5,'a6':a6,'a7':a7,'a8':a8,'a9':a9,
      'b1':b1,'b2':b2,'b3':b3,'b4':b4,'b5':b5,
      'c1':c1,'c2':c2,'c3':c3,'c4':c4,'c5':c5,'c6':c6,
      'd1':d1,'d2':d2,'d3':d3,'d4':d4,
      'e1':e1,'e2':e2,'e3':e3,'e4':e4,
      'v1':ov1,'v2':ov2,'v3':ov3,'v4':ov4,'v5':ov5,
      'v6':np.zeros_like(r),'v7':ov6,'v8':ov7,'v9':ov8,'v10':ov9,
      'v11':np.zeros_like(r),'v12':np.zeros_like(r),'v13':np.zeros_like(r),
      'F_tensor':F,'G_tensor':G,'H_tensor':H,'G4':H/2,'G4phi':G4phi,'G3X':G3X,
      'mu':mu,'P1':P1,'K_scalar':K,'photon_factor':photon,
    })
    return out

SLOTS=[*(f'a{i}' for i in range(1,10)),*(f'b{i}' for i in range(1,6)),*(f'c{i}' for i in range(1,7)),*(f'd{i}' for i in range(1,5)),*(f'e{i}' for i in range(1,5)),*(f'v{i}' for i in range(1,14))]

def scaled_rel(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return np.abs(a-b)/np.maximum(1.0,np.abs(b))

if __name__=='__main__':
    carrier=pd.read_csv(B/'ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv').sort_values('r_over_rs').reset_index(drop=True)
    strong=pd.read_csv(B/'ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv').sort_values('x').reset_index(drop=True)
    # Exact same 5000-point carrier grid; use archived principal action data and selected c2.
    inp=pd.DataFrame({'u':strong.u,'x':strong.x,'phi':strong.phi,'f':strong.f,'h':strong.h,'phi_r':strong.phi_r,'A0prime':np.zeros(len(strong))})
    # Align carrier to ascending x.
    car=carrier.sort_values('r_over_rs').reset_index(drop=True)
    if np.max(np.abs(inp.x.to_numpy()-car.r_over_rs.to_numpy()))>1e-10:
        raise RuntimeError('carrier/strong grids do not align')
    reb=emit(inp,H=car.H_equals_F_equals_G,G4phi=car.G4_phi,G3X=car.G3X,c2_profile=strong.c2)
    rows=[]
    mask=np.abs(reb.photon_factor.to_numpy(float))>2e-3
    for c in SLOTS:
        # a5 is intentionally superseded by the Sep-16 corrected Appendix identity; c3/e3 were absent.
        if c in ['a5','c3','e3'] or c not in strong or strong[c].isna().all(): continue
        sr=scaled_rel(reb[c],strong[c]); use=sr[np.isfinite(sr)&mask]
        rows.append(dict(slot=c,median=float(np.median(use)),p95=float(np.quantile(use,.95)),max=float(np.max(use))))
    for c,src in [('F_tensor','F_tensor'),('G_tensor','G_tensor'),('H_tensor','H_tensor'),('mu','mu'),('K_scalar','K_scalar')]:
        sr=scaled_rel(reb[c],strong[src]); use=sr[np.isfinite(sr)&mask]
        rows.append(dict(slot=c,median=float(np.median(use)),p95=float(np.quantile(use,.95)),max=float(np.max(use))))
    reg=pd.DataFrame(rows)
    reb.to_csv(B/'ssz_p5_MH_G4PHI_EMITTER_STRONGH_REBUILT_2026-09-16.csv',index=False)
    reg.to_csv(B/'ssz_p5_MH_G4PHI_EMITTER_STRONGH_REGRESSION_2026-09-16.csv',index=False)
    print(reg.sort_values('max',ascending=False).to_string(index=False))
