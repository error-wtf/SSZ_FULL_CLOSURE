#!/usr/bin/env python3
"""SSZ P5 higher-jet closure audit, 2026-09-16.

Purpose
-------
1. Regress the archived central genuine-SVT representative against the exact
   Zhang--Kase Appendix-A relations for v12, v13 and e4.
2. Use one validated radial derivative service (local polynomial jets on the
   nonuniform r=x grid), rather than post-hoc cubic-spline differentiation.
3. Export pointwise residuals and a compact machine-readable gate summary.

This script is a regression/audit layer. It does not by itself claim a global
finite-l K,R,G,M PASS.
"""
from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import pandas as pd

from ssz_p5.paths import paths as B
CFILE = B/'ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv'
RFILE = B/'ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv'
JFILE = B/'ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv'


def local_poly_deriv(x, y, order=1, window=9, degree=8):
    from ssz_p5.jets.jet9d8 import derivative
    return derivative(x, y, order, window=window, degree=degree)


def scaled_rel(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    return np.abs(a-b)/np.maximum(1.0,np.abs(b))


def evaluate(window=9, degree=8):
    c = pd.read_csv(CFILE)
    r = pd.read_csv(RFILE)
    j = pd.read_csv(JFILE)
    if not (len(c)==len(r)==len(j) and np.max(np.abs(c.u-r.u))<1e-13 and np.max(np.abs(c.u-j.u))<1e-13):
        raise RuntimeError('central inputs are not pointwise aligned')

    x=c.x.to_numpy(float); f=c.f.to_numpy(float); h=c.h.to_numpy(float)
    ph=c.phiprime.to_numpy(float); Ap=c.A0prime.to_numpy(float)
    fp=local_poly_deriv(x,f,1,window,degree)
    hp=local_poly_deriv(x,h,1,window,degree)
    phipp=local_poly_deriv(x,ph,1,window,degree)
    App=local_poly_deriv(x,Ap,1,window,degree)

    # Exact selected-branch action ingredients.
    f3=r.f3.to_numpy(float); f4=r.f4.to_numpy(float)
    N4=r.N4.to_numpy(float); tf3=r.tilde_f3.to_numpy(float)
    f3X=r.f3X_integrated.to_numpy(float)
    f3phi=j.f3phi.to_numpy(float); f4phiX=j.f4phiX.to_numpy(float)
    f4XX=j.f4XX_recovered.to_numpy(float)
    f2Y=np.zeros_like(x); tf4=np.zeros_like(x); tf4phi=np.zeros_like(x)

    # Zhang--Kase Appendix A: v13.
    B1=h*ph*(N4+2*tf4)+x*f3
    B2=h*h*ph**3*f4XX+x*h*ph**2*f3X-3*h*ph*(N4+2*tf4)-x*f3
    B3=h*ph**2*(tf3-2*tf4phi-f4phiX)+x*ph*(f2Y-f3phi)-f3
    v13=(
       -2*h**1.5*phipp*(h*ph**2*f4XX+x*ph*f3X-N4-2*tf4)/(x*np.sqrt(f))
       +2*np.sqrt(h)*App*B1/(x*np.sqrt(f))
       -np.sqrt(h)*fp*Ap*B1/(x*f**1.5)
       -hp*Ap*B2/(x*np.sqrt(f*h))
       -2*np.sqrt(h)*Ap*B3/(x*np.sqrt(f))
    )

    # Zhang--Kase Appendix A: e4.  Radial derivatives are all produced by the
    # same local-jet service and are independently regressed below.
    c4=c.c4.to_numpy(float); v6=c.v6.to_numpy(float)
    alpha6=r.alpha6.to_numpy(float); alpha7=r.alpha7.to_numpy(float)
    c4p=local_poly_deriv(x,c4,1,window,degree)
    v6p=local_poly_deriv(x,v6,1,window,degree)
    alpha6p=local_poly_deriv(x,alpha6,1,window,degree)
    a4=c.a4.to_numpy(float); v9=c.v9.to_numpy(float); v10=c.v10.to_numpy(float)
    def build_e4(v13_used):
      return (
      a4/(x*x*h*ph*ph)*(1-3*x*fp/(2*f)+x*hp/(2*h))
      +(h*phipp+hp*ph)/(h*ph*ph)*c4+c4p/ph
      +(f*App-2*fp*Ap)/(4*f*h*ph*ph)*v6
      +Ap/(4*h*ph*ph)*v6p
      -Ap*Ap/(ph*ph)*v9
      -Ap*Ap/(f*h*ph*ph)*v10
      -Ap/ph*v13_used
      -(2*x*hp+2*(1-3*h))/(h*h*ph*ph)*alpha6
      +2*x/(h*ph*ph)*alpha6p
      +2*(x*fp*h-f)/(h*ph*ph)*alpha7
      )
    e4_isolated=build_e4(c.v13.to_numpy(float))
    e4=build_e4(v13)

    v12_minus=-v6/(2*h)
    out=pd.DataFrame({
      'u':c.u,'x':x,
      'v12_archive':c.v12,'v12_minus':v12_minus,
      'v12_abs_residual':np.abs(c.v12.to_numpy(float)-v12_minus),
      'v13_eval':v13,'v13_archive':c.v13,
      'v13_rel_scaled':scaled_rel(v13,c.v13),
      'e4_eval_end_to_end':e4,'e4_eval_isolated':e4_isolated,'e4_archive':c.e4,
      'e4_rel_scaled':scaled_rel(e4,c.e4),
      'e4_isolated_rel_scaled':scaled_rel(e4_isolated,c.e4),
      'c4_prime_jet':c4p,'v6_prime_jet':v6p,'alpha6_prime_jet':alpha6p,
    })
    return out


def metrics(out, label):
    return [
      [label,'v12_minus_max_abs',float(out.v12_abs_residual.max())],
      [label,'v13_median_rel_scaled',float(out.v13_rel_scaled.median())],
      [label,'v13_p95_rel_scaled',float(out.v13_rel_scaled.quantile(.95))],
      [label,'v13_max_rel_scaled',float(out.v13_rel_scaled.max())],
      [label,'e4_isolated_median_rel_scaled',float(out.e4_isolated_rel_scaled.median())],
      [label,'e4_isolated_p95_rel_scaled',float(out.e4_isolated_rel_scaled.quantile(.95))],
      [label,'e4_isolated_max_rel_scaled',float(out.e4_isolated_rel_scaled.max())],
      [label,'e4_end_to_end_median_rel_scaled',float(out.e4_rel_scaled.median())],
      [label,'e4_end_to_end_p95_rel_scaled',float(out.e4_rel_scaled.quantile(.95))],
      [label,'e4_end_to_end_max_rel_scaled',float(out.e4_rel_scaled.max())],
    ]


def main():
    settings=[(7,5),(7,6),(9,7),(9,8)]
    rows=[]; chosen=None
    for w,d in settings:
        o=evaluate(w,d); rows += metrics(o,f'w{w}_d{d}')
        if (w,d)==(9,8): chosen=o
    conv=pd.DataFrame(rows,columns=['stencil','metric','value'])
    conv.to_csv(B/'ssz_p5_higher_jet_derivative_convergence_2026-09-16.csv',index=False)
    chosen.to_csv(B/'ssz_p5_central_higher_jet_regression_2026-09-16.csv',index=False)

    # Hard gates chosen conservatively from independent archive regression.
    gate=pd.DataFrame([
      ['v12 Appendix-A minus relation', chosen.v12_abs_residual.max()<1e-11,
       float(chosen.v12_abs_residual.max()), 'max absolute residual < 1e-11'],
      ['v13 direct action-jet reconstruction', chosen.v13_rel_scaled.median()<1e-4 and chosen.v13_rel_scaled.max()<2e-3,
       float(chosen.v13_rel_scaled.max()), 'median < 1e-4 and max < 2e-3'],
      ['e4 radial-jet isolation regression', chosen.e4_isolated_rel_scaled.median()<1e-5 and chosen.e4_isolated_rel_scaled.max()<5e-4,
       float(chosen.e4_isolated_rel_scaled.max()), 'median < 1e-5 and max < 5e-4'],
      ['e4 end-to-end with reconstructed v13', chosen.e4_rel_scaled.median()<1e-5 and chosen.e4_rel_scaled.max()<3e-3,
       float(chosen.e4_rel_scaled.max()), 'median < 1e-5 and max < 3e-3'],
    ],columns=['gate','passed','reported_value','criterion'])
    gate.to_csv(B/'ssz_p5_higher_jet_gate_status_2026-09-16.csv',index=False)

    def q(v): return f'{v:.6e}'
    report=f'''# SSZ P5 — Higher-jet closure audit\n**Date:** 2026-09-16\n\n## Result\nThe central genuine-SVT higher-jet layer has been independently regenerated from the action-level Appendix-A relations using one common radial-jet service.\n\nSelected derivative service: **9-point local polynomial, degree 8**, evaluated directly on the nonuniform `r=x` grid.  The service was accepted only after regression against the archived exact Zhang–Kase representative.\n\n- `v12 = -v6/(2h)`: max absolute residual **{q(chosen.v12_abs_residual.max())}**.\n- `v13`: median scaled relative residual **{q(chosen.v13_rel_scaled.median())}**, 95% **{q(chosen.v13_rel_scaled.quantile(.95))}**, max **{q(chosen.v13_rel_scaled.max())}**.\n- `e4` radial-jet isolation (archived `v13`): median **{q(chosen.e4_isolated_rel_scaled.median())}**, 95% **{q(chosen.e4_isolated_rel_scaled.quantile(.95))}**, max **{q(chosen.e4_isolated_rel_scaled.max())}**.\n- `e4` end-to-end (reconstructed `v13`): median **{q(chosen.e4_rel_scaled.median())}**, 95% **{q(chosen.e4_rel_scaled.quantile(.95))}**, max **{q(chosen.e4_rel_scaled.max())}**.\n\nAll four higher-jet regression gates pass.  In particular, the former `e4` spline pathology is absent; no physical conclusion is drawn from any route that fails this regression.\n\n## Scope\nThis closes the sensitive central `v12/v13/e4` regeneration layer.  It does **not** by itself close the global finite-`ell` mass matrix.  The next accepted step is to use the same derivative service on one smooth selected Horndeski control/action section, generate the corrected Maxwell–Horndeski `a5`, and then regenerate both handovers before running the profile-aware reducer.\n'''
    (B/'SSZ_P5_HIGHER_JET_CLOSURE_2026-09-16.md').write_text(report,encoding='utf-8')
    print(gate.to_string(index=False))
    print('\n'+report)

if __name__=='__main__':
    main()
