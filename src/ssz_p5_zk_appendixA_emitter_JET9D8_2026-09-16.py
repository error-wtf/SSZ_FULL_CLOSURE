#!/usr/bin/env python3
"""Zhang–Kase Appendix-A coefficient emitter for the selected P5 U(1)-SVT branch.

Date: 2026-09-16

Purpose
-------
Reconstruct the canonical 41-slot unreduced even-parity coefficient schema
  a1..a9, b1..b5, c1..c6, d1..d4, e1..e4, v1..v13
from one radial action-jet profile, using the accepted JET9D8 radial derivative
service.  The lower-order transverse controls v5,c3,e3 can be supplied by the
selected member; by default they are set to NaN so that this module can also
regress the historical 37/41 pure-SVT table without silently inventing them.

The implementation follows Zhang & Kase, PRD 110, 044047 (2024), Appendix A,
with the P5 branch conventions M_pl^2=1, f2_Y=0, f2_FtildeFtilde=0,
tilde f4=0 unless explicitly supplied.

Important project-specific holonomic convention
------------------------------------------------
The d3 term contains a phi-derivative of v6.  For the selected on-curve P5
representative this is evaluated as (dv6/dr)/phi', which independently
regresses the archived 37/41 outer table.  This is an on-curve holonomic
selector, not a claim about an arbitrary off-curve extension.
"""
from __future__ import annotations
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd

from ssz_p5.paths import paths as B

# Accepted derivative service.
spec = importlib.util.spec_from_file_location('hj', B/'ssz_p5_higher_jet_closure_2026-09-16.py')
hj = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hj)


def dr(x, y, order=1, window=9, degree=8):
    return hj.local_poly_deriv(np.asarray(x,float), np.asarray(y,float), order, window, degree)


def _arr(df, name, default=0.0):
    if name in df.columns:
        return df[name].to_numpy(float)
    return np.full(len(df), float(default))


def emit(df: pd.DataFrame, *, window=9, degree=8,
         selected_v5=None, selected_c3=None, selected_e3=None) -> pd.DataFrame:
    """Emit Zhang–Kase Appendix-A coefficients for one aligned radial action-jet table.

    Required columns:
      x,f,h,phiprime,A0prime,X,f2X,f2F,f3,f3X,f4,f4X
    Optional higher jets use zero defaults only where the selected branch declares
    them absent.  f2 Hessians use names f2XX,f2XF,f2XY,f2FF,f2FY,f2YY.
    """
    d=df.copy().reset_index(drop=True)
    req=['x','f','h','phiprime','A0prime','X','f2X','f2F','f3','f3X','f4','f4X']
    miss=[k for k in req if k not in d.columns]
    if miss: raise KeyError(f'missing required action/background columns: {miss}')

    r=_arr(d,'x'); f=_arr(d,'f'); h=_arr(d,'h'); ph=_arr(d,'phiprime'); Ap=_arr(d,'A0prime'); X=_arr(d,'X')
    if np.any(f<=0) or np.any(h<=0) or np.any(np.abs(ph)<1e-14):
        raise ValueError('requires f>0, h>0, phiprime != 0')
    fp=dr(r,f,1,window,degree); hp=dr(r,h,1,window,degree)
    phipp=dr(r,ph,1,window,degree); App=dr(r,Ap,1,window,degree); Xp=dr(r,X,1,window,degree)

    f2X=_arr(d,'f2X'); f2F=_arr(d,'f2F'); f2Y=_arr(d,'f2Y')
    f2XX=_arr(d,'f2XX'); f2XF=_arr(d,'f2XF'); f2XY=_arr(d,'f2XY')
    f2FF=_arr(d,'f2FF'); f2FY=_arr(d,'f2FY'); f2YY=_arr(d,'f2YY')
    f2FtFt=_arr(d,'f2FtFt')

    f3=_arr(d,'f3'); f3X=_arr(d,'f3X'); f3XX=_arr(d,'f3XX')
    tf3=_arr(d,'tf3')
    f4=_arr(d,'f4'); f4X=_arr(d,'f4X'); f4XX=_arr(d,'f4XX'); f4XXX=_arr(d,'f4XXX')
    tf4=_arr(d,'tf4')

    # Holonomic mixed phi jets: use supplied values when present, otherwise infer
    # from the on-curve chain rule.
    f3phi=_arr(d,'f3phi', np.nan)
    if np.all(np.isnan(f3phi)): f3phi=(dr(r,f3,1,window,degree)-f3X*Xp)/ph
    f3phiX=_arr(d,'f3phiX', np.nan)
    if np.all(np.isnan(f3phiX)): f3phiX=(dr(r,f3X,1,window,degree)-f3XX*Xp)/ph
    f4phi=_arr(d,'f4phi', np.nan)
    if np.all(np.isnan(f4phi)): f4phi=(dr(r,f4,1,window,degree)-f4X*Xp)/ph
    f4phiX=_arr(d,'f4phiX', np.nan)
    if np.all(np.isnan(f4phiX)): f4phiX=(dr(r,f4X,1,window,degree)-f4XX*Xp)/ph
    tf4phi=_arr(d,'tf4phi')

    # Odd-sector coefficients Eq. (3.6), selected P5 branch.
    alpha1=np.sqrt(h)/(4*np.sqrt(f))
    alpha2=h**1.5*Ap/(2*r*np.sqrt(f))*(r*ph*f3-4*f4+h*ph**2*(f4X+2*tf4))
    alpha4=(r*f2F + 2*h*ph*(f3+h*ph**2*tf3) -4*hp*f4
            +(2*h*phipp+ph*hp)*(r*f3+h*ph*(f4X+2*tf4)))/(2*r*np.sqrt(f*h))
    alpha5=-(np.sqrt(f*h)/(2*r))*(r*(f2F-2*h*ph**2*f2Y)+2*h*ph*f3
             +(h*fp/f)*(r*ph*f3-4*f4+h*ph**2*(f4X+2*tf4)))
    alpha6=-(np.sqrt(f*h)/(4*r**2))*(1-(h*Ap**2/f)*(4*f4-h*ph**2*(f4X+2*tf4)))
    alpha7=(1-(4*h*Ap**2/f)*f4)/(4*r**2*np.sqrt(f*h))

    # Main Hessian-dependent channels in Appendix A.
    c2=(0.5*r**2*np.sqrt(f*h)*ph*(f2X-h*ph**2*f2XX)
        -r**2*h**2.5*ph*Ap**4*(4*h*ph**2*f2YY-f2FY)/(f**1.5)
        -h**1.5*Ap**2/(2*np.sqrt(f))*(
            h**3*ph**5*f4XXX + 2*r*h**2*ph**4*f3XX
            +h*ph**3*(6*r**2*f2XY+(4-13*h)*f4XX)
            -12*r*h*ph**2*f3X
            -ph*(r**2*(6*f2Y+f2XF)+6*(2-5*h)*f4X-20*h*tf4)
            +6*r*f3))

    v1=(r**2*h**1.5*Ap**2/(2*f**1.5)*(f2FF-4*h*ph**2*(f2FY-h*ph**2*f2YY))
        -np.sqrt(h)/(2*np.sqrt(f))*(2*h*ph**2*(r**2*f2Y-h*(f4X+2*tf4))
          -4*r*h*ph*f3-8*(1-h)*f4-r**2*f2F))
    v4=(2*r**2*h**2.5*ph*Ap**3*(2*h*ph**2*f2YY-f2FY)/(f**1.5)
        +h**1.5*Ap/np.sqrt(f)*(2*h*ph**3*(r**2*f2XY-h*f4XX)
          -4*r*h*ph**2*f3X-r**2*ph*(4*f2Y+f2XF)
          +4*ph*((3*h-2)*f4X+2*h*tf4)+4*r*f3))

    v6=4*alpha2; v9=alpha4; v10=alpha5
    v2=Ap*v1
    v3=-Ap*v1-ph*v4/2-r*v6
    v7=Ap**2*v1/4
    v8=((r*fp-2*f)*v6+4*r*Ap*v10)/(2*r*f)
    v11=-v6/2; v12=-v6/(2*h)

    c4=-(np.sqrt(h)*Ap**2/(2*r*np.sqrt(f)))*(h**2*ph**3*f4XX+r*h*ph**2*f3X-3*h*ph*(f4X+2*tf4)-r*f3)
    d2=-(h**1.5*Ap**2/(r*np.sqrt(f)))*(h**2*ph**3*f4XX+r*h*ph**2*f3X-2*h*ph*(3*f4X+2*tf4)-r*f3)

    # v13 Appendix A.
    B1=h*ph*(f4X+2*tf4)+r*f3
    B2=h*h*ph**3*f4XX+r*h*ph**2*f3X-3*h*ph*(f4X+2*tf4)-r*f3
    B3=h*ph**2*(tf3-2*tf4phi-f4phiX)+r*ph*(f2Y-f3phi)-f3
    v13=(-2*h**1.5*phipp*(h*ph**2*f4XX+r*ph*f3X-f4X-2*tf4)/(r*np.sqrt(f))
         +2*np.sqrt(h)*App*B1/(r*np.sqrt(f))
         -np.sqrt(h)*fp*Ap*B1/(r*f**1.5)
         -hp*Ap*B2/(r*np.sqrt(f*h))
         -2*np.sqrt(h)*Ap*B3/(r*np.sqrt(f)))

    a4=np.sqrt(f*h)/2; a4p=dr(r,a4,1,window,degree)
    a1=np.zeros_like(r)
    a2=r*(fp*h-f*hp)/(f*h*ph)*a4+Ap*v4/2+r*Ap*v6/(2*ph)
    a3=-r*a4
    a6=(r*(fp*h+f*hp)-2*f*h)/(2*r*f*h**2*ph)*a4-a4p/(h*ph)-Ap*v6/(4*h*ph)+2*r*f*alpha7/ph
    a7=-dr(r,r*a4,1,window,degree)-Ap*(2*Ap*v1+ph*v4)/4-r*Ap*v6/2
    a8=-a4/(2*h)
    a9=a4p-(r*fp-2*f)*a4/(2*r*f)+Ap*v6/4

    b1=a4/(2*f); b2=np.zeros_like(r); b3=-2*a2/f+Ap*v4/f; b4=2*r*a4/f; b5=-a4/f
    c1=np.zeros_like(r)
    c5=-fp*a4/(2*f)-ph*d2/2-Ap*v6/2+2*r*alpha6
    c6=r*fp*a4/(2*f)-ph*c2/4+r*ph*d2/4+Ap**2*v1/4+ph*Ap*v4/8+r*Ap*v6/2-r**2*alpha6
    d1=a4/(2*f); d4=-2*alpha6

    # The archived P5 outer table fixes this on-curve phi derivative through
    # holonomicity: partial_phi v6 -> (d v6/dr)/phi'.
    v6phi=dr(r,v6,1,window,degree)/ph
    d3=2*(fp*h-f*hp)/(r*f*h*ph)*a4+Ap*v6/(r*ph)+Ap*v6phi/2+4*alpha6/(h*ph)+4*f*alpha7/ph

    e1=(a2-2*r*h*a6-Ap*v4/2)/(f*h*ph)
    e2=-(c2+r*d2+Ap*v4/2)/ph

    # Lower-order selected controls.
    def chosen(v):
        if v is None: return np.full(len(d),np.nan)
        if np.isscalar(v): return np.full(len(d),float(v))
        return np.asarray(v,float)
    v5=chosen(selected_v5); c3=chosen(selected_c3); e3=chosen(selected_e3)

    # a5: selected common holonomic completion including generalized a1 term.
    # This equals Zhang–Kase Appendix A in the pure ZK limit a1=0.
    v4p=dr(r,v4,1,window,degree); a2p=dr(r,a2,1,window,degree); a1pp=dr(r,dr(r,a1,1,window,degree),1,window,degree)
    a5=a2p-a1pp-App*v4/2+Ap*(v5-v4p)/2

    # e4 Appendix A with one common radial-jet service.
    c4p=dr(r,c4,1,window,degree); v6p=dr(r,v6,1,window,degree); alpha6p=dr(r,alpha6,1,window,degree)
    e4=(a4/(r*r*h*ph*ph)*(1-3*r*fp/(2*f)+r*hp/(2*h))
        +(h*phipp+hp*ph)/(h*ph*ph)*c4+c4p/ph
        +(f*App-2*fp*Ap)/(4*f*h*ph*ph)*v6+Ap*v6p/(4*h*ph*ph)
        -Ap*Ap*v9/(ph*ph)-Ap*Ap*v10/(f*h*ph*ph)-Ap*v13/ph
        -(2*r*hp+2*(1-3*h))*alpha6/(h*h*ph*ph)+2*r*alpha6p/(h*ph*ph)
        +2*(r*fp*h-f)*alpha7/(h*ph*ph))

    out=pd.DataFrame(index=d.index)
    for k in ['u','x','phi','f','h','phiprime','A0prime','X']:
        if k in d.columns: out[k]=d[k]
    vals={
      'a1':a1,'a2':a2,'a3':a3,'a4':a4,'a5':a5,'a6':a6,'a7':a7,'a8':a8,'a9':a9,
      'b1':b1,'b2':b2,'b3':b3,'b4':b4,'b5':b5,
      'c1':c1,'c2':c2,'c3':c3,'c4':c4,'c5':c5,'c6':c6,
      'd1':d1,'d2':d2,'d3':d3,'d4':d4,
      'e1':e1,'e2':e2,'e3':e3,'e4':e4,
      'v1':v1,'v2':v2,'v3':v3,'v4':v4,'v5':v5,'v6':v6,'v7':v7,'v8':v8,'v9':v9,'v10':v10,'v11':v11,'v12':v12,'v13':v13,
      'alpha1':alpha1,'alpha2':alpha2,'alpha4':alpha4,'alpha5':alpha5,'alpha6':alpha6,'alpha7':alpha7,
      'f3phi_holonomic':f3phi,'f3phiX_holonomic':f3phiX,'f4phi_holonomic':f4phi,'f4phiX_holonomic':f4phiX,
    }
    for k,v in vals.items(): out[k]=v
    return out


def _scaled_rel(a,b):
    return np.abs(np.asarray(a,float)-np.asarray(b,float))/np.maximum(1.0,np.abs(np.asarray(b,float)))


def regress_outer_archive():
    src=pd.read_csv(B/'ssz_p5_outer_svt_repaired_integrable_transition_profile_2026-09-12.csv')
    ref=pd.read_csv(B/'ssz_p5_F2_outer_transition_unreduced_PRINCIPAL_37of41_2026-09-15.csv')
    if len(src)!=len(ref) or np.max(np.abs(src.u-ref.u))>1e-12:
        raise RuntimeError('outer source/reference grids are not aligned')
    d=pd.DataFrame({
      'u':ref.u,'x':ref.x,'phi':ref.phi,'f':ref.f,'h':ref.h,'phiprime':ref.phiprime,
      'A0prime':ref.A0prime,'X':ref.X,
      'f2X':src.f2X,'f2F':np.ones(len(src)),'f2Y':0.0,
      'f2XX':src.HXX,'f2XF':src.HXF,'f2XY':src.HXY,'f2FF':src.HFF,'f2FY':src.HFY,'f2YY':src.HYY,
      'f3':src.f3,'f3X':src.f3X_integrated,'f3XX':0.0,'tf3':src.tilde_f3,
      'f4':src.f4,'f4X':src.N4,'f4XX':0.0,'f4XXX':0.0,'tf4':0.0,
    })
    got=emit(d)
    slots=[*[f'a{i}' for i in range(1,10)],*[f'b{i}' for i in range(1,6)],*[f'c{i}' for i in range(1,7)],
           *[f'd{i}' for i in range(1,5)],*[f'e{i}' for i in range(1,5)],*[f'v{i}' for i in range(1,14)]]
    rows=[]
    for k in slots:
        if k in ['a5','v5','c3','e3']: continue
        rel=_scaled_rel(got[k],ref[k])
        rows.append([k,float(np.nanmedian(rel)),float(np.nanquantile(rel,.95)),float(np.nanmax(rel)),float(np.nanmax(np.abs(got[k]-ref[k])) )])
    audit=pd.DataFrame(rows,columns=['slot','median_scaled_rel','p95_scaled_rel','max_scaled_rel','max_abs'])
    audit.to_csv(B/'ssz_p5_ZK_APPENDIX_A_EMITTER_OUTER_REGRESSION_2026-09-16.csv',index=False)
    got.to_csv(B/'ssz_p5_ZK_APPENDIX_A_EMITTER_OUTER_REBUILT_2026-09-16.csv',index=False)
    return audit


if __name__=='__main__':
    a=regress_outer_archive()
    print(a.to_string(index=False))
    print('\nworst median:',a.loc[a.median_scaled_rel.idxmax()].to_dict())
    print('worst p95:',a.loc[a.p95_scaled_rel.idxmax()].to_dict())
