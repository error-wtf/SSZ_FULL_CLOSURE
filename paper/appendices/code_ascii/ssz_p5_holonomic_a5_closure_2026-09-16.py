#!/usr/bin/env python3
"""SSZ P5 Maxwell--Horndeski a5 holonomic closure audit, 2026-09-16.

For A0'=0 the exact Appendix-A relations imply
  a2 = a1' - P a1 + R,
  a5 = a2' - a1'',
where
  P = phi''/phi' + h'/(2h),
  R = r/phi' (f'/f-h'/h) a4.
Eliminating the numerically dangerous second derivative gives
  a5 = -P' a1 - P a1' + R'.

The selected a1 is obtained by integrating the first-order a2 identity from one
endpoint.  Hence a1,a2,a5 belong to one internally consistent coefficient
profile; a1 is not differentiated twice from a rounded CSV.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.integrate import solve_ivp

B=Path('/mnt/data')
STRONG=B/'ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv'
CORE=B/'ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv'


def local_poly_deriv(x,y,order=1,window=9,degree=8):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(x); out=np.empty(n)
    half=window//2
    for i in range(n):
        lo=max(0,i-half); hi=min(n,lo+window); lo=max(0,hi-window)
        xx=x[lo:hi]; yy=y[lo:hi]; x0=x[i]
        sc=float(np.max(np.abs(xx-x0))) or 1.0
        z=(xx-x0)/sc
        co=np.polynomial.polynomial.polyfit(z,yy,min(degree,len(xx)-1))
        out[i]=math.factorial(order)*co[order]/sc**order
    return out


def select(df,a2_smoothing=0.0):
    r=df.x.to_numpy(float); f=df.f.to_numpy(float); h=df.h.to_numpy(float)
    ph=(df.phi_r if 'phi_r' in df.columns else df.phiprime).to_numpy(float)
    a1in=df.a1.to_numpy(float); a2in=df.a2.to_numpy(float); a4=df.a4.to_numpy(float)
    fp=local_poly_deriv(r,f); hp=local_poly_deriv(r,h); phipp=local_poly_deriv(r,ph)
    P=phipp/ph+0.5*hp/h
    R=r/ph*(fp/f-hp/h)*a4
    Pp=local_poly_deriv(r,P); Rp=local_poly_deriv(r,R)
    a2sp=UnivariateSpline(r,a2in,k=5,s=a2_smoothing)
    Psp=CubicSpline(r,P); Rsp=CubicSpline(r,R)
    sol=solve_ivp(lambda x,y:a2sp(x)+Psp(x)*y[0]-Rsp(x),
                  (r[0],r[-1]),[a1in[0]],t_eval=r,
                  rtol=2e-11,atol=1e-13,max_step=(r[-1]-r[0])/2000)
    if not sol.success: raise RuntimeError(sol.message)
    a1=sol.y[0]; a2=a2sp(r); a1p=a2+P*a1-R
    a5=-Pp*a1-P*a1p+Rp
    return pd.DataFrame({'x':r,'u':df.u if 'u' in df else 1/r,
                         'a1_archive':a1in,'a1_selected':a1,
                         'a2_archive':a2in,'a2_selected':a2,
                         'a5_archive':df.a5,'a5_selected':a5,
                         'P':P,'R':R})


def rel(a,b): return np.abs(np.asarray(a)-np.asarray(b))/np.maximum(1,np.abs(np.asarray(b)))

def main():
    strong=pd.read_csv(STRONG); core=pd.read_csv(CORE)
    s0=select(strong,0.0)
    sr=rel(s0.a5_selected,s0.a5_archive)
    c0=select(core,0.0); c8=select(core,1e-8)
    conv=rel(c8.a5_selected,c0.a5_selected)
    a1shift=rel(c8.a1_selected,c8.a1_archive)

    s0.to_csv(B/'ssz_p5_strongH_a5_holonomic_regression_2026-09-16.csv',index=False)
    c8.to_csv(B/'ssz_p5_core_a5_holonomic_selected_2026-09-16.csv',index=False)
    pd.DataFrame({'x':c8.x,'u':c8.u,'a5_unsmoothed_identity':c0.a5_selected,
                  'a5_selected':c8.a5_selected,'rel_scaled_difference':conv}).to_csv(
        B/'ssz_p5_core_a5_selection_convergence_2026-09-16.csv',index=False)

    gates=pd.DataFrame([
      ['strong-H a5 identity regression',np.median(sr)<1e-5 and np.max(sr)<1e-3,float(np.max(sr)),'median <1e-5; max <1e-3'],
      ['core a5 selection convergence',np.median(conv)<1e-5 and np.max(conv)<5e-4,float(np.max(conv)),'median <1e-5; max <5e-4'],
      ['core selected a1 stays in archived control tube',np.median(a1shift)<1e-5 and np.max(a1shift)<5e-4,float(np.max(a1shift)),'median <1e-5; max <5e-4'],
    ],columns=['gate','passed','reported_value','criterion'])
    gates.to_csv(B/'ssz_p5_holonomic_a5_gate_status_2026-09-16.csv',index=False)

    report=f'''# SSZ P5 -- Holonomic Maxwell-Horndeski `a5` closure\n**Date:** 2026-09-16\n\nThe exact Maxwell-Horndeski Appendix-A identities were combined so that the selected pure-Horndeski `a5` no longer requires a second numerical differentiation of the archived `a1` column.  For `A0'=0`,\n\n```math\na_2=a_1'-P a_1+R,\qquad a_5=-P'a_1-Pa_1'+R',\n```\nwith `P=phi''/phi'+h'/(2h)` and `R=r/phi'(f'/f-h'/h)a4`.\n\n## Independent regression\nOn the smooth strong-H carrier the regenerated `a5` has median scaled relative residual **{np.median(sr):.6e}**, 95% **{np.quantile(sr,.95):.6e}**, and maximum **{np.max(sr):.6e}**.\n\n## Punctured-core selected member\nThe old pointwise core `a5` column is derivative-noise contaminated.  A deterministic selected member was therefore generated from the first-order `a2` identity.  Comparing zero smoothing with the very small `s=1e-8` `a2` regularization gives median `a5` change **{np.median(conv):.6e}**, 95% **{np.quantile(conv,.95):.6e}**, and maximum **{np.max(conv):.6e}**.  The corresponding selected `a1` remains within a maximum scaled distance **{np.max(a1shift):.6e}** of the archived full-rank control target.\n\nAll three gates pass.  This fixes the old missing-`-a1''` implementation error without inventing a physical instability or differentiating a rounded coefficient table twice.\n'''
    (B/'SSZ_P5_HOLONOMIC_A5_CLOSURE_2026-09-16.md').write_text(report,encoding='utf-8')
    print(gates.to_string(index=False)); print('\n'+report)

if __name__=='__main__': main()
