"""G13 -> G14 -> G15 -> G16: electric-branch gate chain on the CURRENT member.

G13: production electric specialization (f2Y = 0 branch).
G14: full light-ring limit (undivided; u_lr / A0_u symbol separation).
G15: Sigma_SVT operator decomposition with the DERIVED Delta22_SVT.
G16: full on-shell light-ring balance (C_bg with genuine-SVT theta-theta).
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp

from ssz_p5.action import light_ring_identity as lri
from ssz_p5.action import svt_direct_variation as sdv
from ssz_p5.action.symbolic_canon import canonical_simplify
from ssz_p5.jets.jet9d8 import profile_derivative
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.svt_background_eom import evaluate_svt_background

ROOT = Path(__file__).resolve().parents[2]


def _member_numeric():
    """Current member frame with the canonical derivative service applied."""
    build = build_onshell_central(ROOT)
    d = build.action
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    h = d.h.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    ap = d.A0prime.to_numpy(float)
    frame = {
        "d": d, "r": r, "f": f, "h": h, "ph": ph, "ap": ap,
        "fp": profile_derivative(r, f, 1, 9, 8),
        "fpp": profile_derivative(r, f, 2, 9, 8),
        "hp": profile_derivative(r, h, 1, 9, 8),
        "phpp": profile_derivative(r, ph, 1, 9, 8),
        "ap2": profile_derivative(r, ap, 1, 9, 8),
        "mask": (d.u.to_numpy(float) > 0.62) & (d.u.to_numpy(float) < 0.70),
    }
    return frame


def _delta22_numeric(fr):
    d22 = sdv.build_delta22_svt(verify=False)
    args = [sdv.r, sdv.f, sdv.h, sdv.ph, sdv.ap, sdv.fp_r, sdv.fpp_r, sdv.hp_r,
            sdv.phpps, sdv.ap2s, sdv.f3, sdv.f3X, sdv.f4, sdv.f4X, sdv.f4XX,
            sdv.tf4, sdv.f3phis, sdv.f4phis, sdv.f4Xphis, sdv.tf4phis]
    d = fr["d"]
    vals = [fr["r"], fr["f"], fr["h"], fr["ph"], fr["ap"], fr["fp"], fr["fpp"],
            fr["hp"], fr["phpp"], fr["ap2"],
            d.f3.to_numpy(float), d.f3X.to_numpy(float), d.f4.to_numpy(float),
            d.f4X.to_numpy(float), d.f4XX.to_numpy(float), d.tf4.to_numpy(float),
            d.f3phi.to_numpy(float), d.f4phi.to_numpy(float),
            d.f4phiX.to_numpy(float), d.tf4phi.to_numpy(float)]
    return np.asarray(sp.lambdify(args, d22)(*vals), float)


def test_G13_production_f2Y_zero_branch():
    """G13: the production member IS the f2Y = 0 branch (column exactly zero),
    the locked HT evaluator accepts it, and the symbolic C_bg machinery keeps
    f2Y symbolic through Test B (specialization only after G11/G12)."""
    fr = _member_numeric()
    d = fr["d"]
    assert float(np.max(np.abs(d.f2Y.to_numpy(float)))) == 0.0
    res = evaluate_svt_background(d, require_f2y_zero=True)  # raises if f2Y != 0
    mask = fr["mask"]
    # NUMERICAL_POLICY.background_residual_abs = 1e-10
    assert np.max(np.abs(res.JA[mask])) < 1e-10
    assert np.max(np.abs(res.E00[mask])) < 1e-10
    assert np.max(np.abs(res.E11[mask])) < 1e-10
    # symbolic protocol intact: Test B exact with f2Y SYMBOLIC
    ok, residual = lri.test_B_check()
    assert ok
    # and the derived Delta22_SVT carries no f2Y channel at all
    assert not sdv.build_delta22_svt(verify=False).has(lri.f2Y)


def test_G14_full_light_ring_limit():
    """G14: undivided light-ring structure on the current member - the ring
    denominator 2f - r f' vanishes at the inner ring while no division is
    performed; u_lr and A0_u remain distinct symbols."""
    fr = _member_numeric()
    d = fr["d"]
    u = d.u.to_numpy(float)
    ring_expr = 2*fr["f"] - fr["r"]*fr["fp"]
    i_lr = int(np.argmin(np.abs(u - 0.706135)))
    # ring denominator vanishes at the inner stable ring (undivided handling)
    assert abs(ring_expr[i_lr]) < 5e-3
    # symbol separation regression stays exact
    dd = lri.lr_ring_limit_expressions()
    S = dd["symbols"]
    assert S["u_lr"] != S["A0_u"]
    assert dd["geometric_ring_term_raw"].has(S["u_lr"])
    assert dd["electric_normalized"].has(S["A0_u"])
    # undivided C85 carries the ring limit by construction (no 1/(2f - r f'))
    assert not sp.sympify(lri.C85_undivided()).has(1/(2*lri.f - lri.r*lri.fp_r))


def test_G15_sigma_svt_decomposition_derived():
    """G15: Sigma_SVT = C_bg - P_MH[C_bg] decomposes exactly into the
    documented operator list with the DERIVED Delta22_SVT; numerically the
    Sigma decomposition on the member equals BETA * Delta22_SVT."""
    d = lri.sigma_svt_decomposition()
    from ssz_p5.action.symbolic_canon import canonical_simplify
    total = d["total"]
    s = sum(sp.sympify(op["symbolic"]) for op in d["operators"])
    assert canonical_simplify(total - s) == 0
    d22 = lri.get_delta22_svt_expr()
    # the theta-theta operator entry is EXACTLY BETA * Delta22_SVT
    assert canonical_simplify(sp.sympify(d["operators"][2]["symbolic"])
                              - lri.BETA*d22) == 0
    # numeric Sigma content on the current member is nontrivial
    fr = _member_numeric()
    Delta = _delta22_numeric(fr)
    BETA = fr["r"]*fr["f"].astype(float)**1.5/np.sqrt(fr["h"])
    sigma_direct = BETA*Delta
    assert np.max(np.abs(sigma_direct[fr["mask"]])) > 0.0


def test_G16_full_onshell_LR_balance():
    """G16: full on-shell light-ring balance on the CURRENT member.

    C_bg = ALPHA*E00_full + BETA*E22_full with the DERIVED genuine-SVT
    theta-theta correction satisfies the on-shell balance at the resolution
    level of the derivative service, while the MH slice ALONE is far off
    (|E22_MH| = O(1) vs |E22_full| ~ 1e-6) - the negative control is built in.
    """
    fr = _member_numeric()
    d = fr["d"]
    mask = fr["mask"]
    res = evaluate_svt_background(d, require_f2y_zero=True)
    Delta = _delta22_numeric(fr)
    E22mh_num = _e22_mh_numeric(fr)
    E22_full = E22mh_num + Delta
    ALPHA = np.sqrt(fr["f"])/(fr["r"]*np.sqrt(fr["h"]))
    BETA = fr["r"]*fr["f"]**sp.Rational(3, 2)/np.sqrt(fr["h"])
    # strict machine-level channel (E00_full), policy background_residual_abs
    assert np.max(np.abs(res.E00[mask])) < 1e-10
    # full balance: |C_bg| = |ALPHA E00 + BETA E22_full| at resolution level
    C_bg_full = ALPHA*res.E00 + BETA*E22_full
    assert np.max(np.abs(C_bg_full[mask])) < 1e-5
    # negative control: WITHOUT Delta22 the balance fails by O(1)
    C_bg_MH = ALPHA*res.E00 + BETA*E22mh_num
    assert np.max(np.abs(C_bg_MH[mask])) > 1.0
    # at the inner light ring specifically
    i_lr = int(np.argmin(np.abs(d.u.to_numpy(float) - 0.706135)))
    assert abs(float(C_bg_full[i_lr])) < 1e-5
    assert abs(float(C_bg_MH[i_lr])) > 1.0


def _e22_mh_numeric(fr):
    d = fr["d"]
    f2 = d.f2.to_numpy(float)
    return np.asarray(sp.lambdify([lri.r, lri.f, lri.h, lri.ph, lri.fp_r,
                                   lri.fpp_r, lri.hp_r, lri.f2],
                                  lri.E22_MH_slice())(
        fr["r"], fr["f"], fr["h"], fr["ph"], fr["fp"], fr["fpp"], fr["hp"], f2),
        float)
