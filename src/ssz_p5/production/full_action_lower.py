"""Direct Appendix-A lower slots from one holonomic f2+f3+f4 action.

This module removes the historical selected-coefficient shortcut for v5/c3/e3.
It works on the static electric Zhang--Kase background and differentiates the
on-curve action jets before any regional trimming.

Conventions follow arXiv:2404.11910v3, Eqs. (2.11)--(2.18) and Appendix A.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..jets.jet9d8 import profile_derivative


def _arr(d: pd.DataFrame, name: str, default=None) -> np.ndarray:
    if name in d:
        return d[name].to_numpy(float)
    if default is None:
        raise KeyError(name)
    if np.isscalar(default):
        return np.full(len(d), float(default))
    return np.asarray(default, float)


def complete_total_action_jets(df: pd.DataFrame, *, window: int = 9, degree: int = 8) -> pd.DataFrame:
    """Complete mixed phi jets by the total-action on-curve chain rules.

    Unlike a background-null *additive* Hessian, a total action obeys H t = g',
    not H t = 0.  The six transverse f2 Hessian entries plus f3XX/f4XX/f4XXX
    determine the mixed jets once the first action jets are fixed.
    """
    d = df.copy().reset_index(drop=True)
    r = _arr(d, "x")
    f, h, X, A = (_arr(d, k) for k in ("f", "h", "X", "A0prime"))
    ph = _arr(d, "phiprime", _arr(d, "phi_r") if "phi_r" in d else None)
    if np.any(np.diff(r) <= 0):
        raise ValueError("x must be strictly increasing")
    if np.any(np.abs(ph) <= 1e-14):
        raise ValueError("nonzero phiprime required")

    def dr(y):
        return profile_derivative(r, np.asarray(y, float), 1, window, degree)

    F = h * A**2 / (2.0 * f)
    Y = 4.0 * X * F
    Xp, Fp, Yp = dr(X), dr(F), dr(Y)

    # f2(phi,X,F,Y): total-action Hessian completion H t = g'.
    for name in ("f2phi", "f2X", "f2F", "f2Y", "f2XX", "f2XF", "f2XY", "f2FF", "f2FY", "f2YY"):
        if name not in d:
            if name == "f2Y":
                d[name] = 0.0
            else:
                raise KeyError(name)
    d["f2phiX"] = (dr(d.f2X) - d.f2XX * Xp - d.f2XF * Fp - d.f2XY * Yp) / ph
    d["f2phiF"] = (dr(d.f2F) - d.f2XF * Xp - d.f2FF * Fp - d.f2FY * Yp) / ph
    d["f2phiY"] = (dr(d.f2Y) - d.f2XY * Xp - d.f2FY * Fp - d.f2YY * Yp) / ph
    d["f2phiphi"] = (
        dr(d.f2phi) - d.f2phiX * Xp - d.f2phiF * Fp - d.f2phiY * Yp
    ) / ph

    # f3(phi,X).
    for name in ("f3", "f3X", "f3XX"):
        if name not in d:
            raise KeyError(name)
    d["f3phi"] = (dr(d.f3) - d.f3X * Xp) / ph
    d["f3phiX"] = (dr(d.f3X) - d.f3XX * Xp) / ph
    d["f3phiphi"] = (dr(d.f3phi) - d.f3phiX * Xp) / ph

    # f4(phi,X); tilde-f4 is zero on the selected P5 branch unless supplied.
    for name in ("f4", "f4X", "f4XX", "f4XXX"):
        if name not in d:
            raise KeyError(name)
    d["f4phi"] = (dr(d.f4) - d.f4X * Xp) / ph
    d["f4phiX"] = (dr(d.f4X) - d.f4XX * Xp) / ph
    d["f4phiXX"] = (dr(d.f4XX) - d.f4XXX * Xp) / ph
    d["f4phiphi"] = (dr(d.f4phi) - d.f4phiX * Xp) / ph
    d["f4phiphiX"] = (dr(d.f4phiX) - d.f4phiXX * Xp) / ph
    if "tf4phi" not in d:
        d["tf4phi"] = 0.0
    if "tf4phiphi" not in d:
        d["tf4phiphi"] = 0.0
    d["Fbg_action"] = F
    d["Ybg_action"] = Y
    return d


def emit_lower_slots(df: pd.DataFrame, *, window: int = 9, degree: int = 8) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return direct Appendix-A (v5,c3,e3) and the completed total-action jets."""
    d = complete_total_action_jets(df, window=window, degree=degree)
    r, f, h, ph, A = (_arr(d, c) for c in ("x", "f", "h", "phiprime", "A0prime"))

    def dr(y):
        return profile_derivative(r, np.asarray(y, float), 1, window, degree)

    # Appendix A v5.
    v5 = np.sqrt(h / f) * A * (
        2.0 * h * ph**2 * (h * (d.f4phiX + 2.0 * d.tf4phi) - r**2 * d.f2phiY)
        + 4.0 * r * h * ph * d.f3phi
        + 8.0 * (1.0 - h) * d.f4phi
        + r**2 * d.f2phiF
    )

    # c3 = -(2 sqrt(fh))^-1 partial_phi E11.  E11 = metric - bracket.
    bphi = (
        r**2
        * (
            f * d.f2phi
            + f * h * ph**2 * d.f2phiX
            - h * A**2 * (d.f2phiF - 4.0 * h * ph**2 * d.f2phiY)
        )
        - 2.0 * r * h**2 * ph * A**2 * (3.0 * d.f3phi - h * ph**2 * d.f3phiX)
        + h
        * A**2
        * (
            4.0 * (3.0 * h - 1.0) * d.f4phi
            - h * (9.0 * h - 4.0) * ph**2 * d.f4phiX
            + h**3 * ph**4 * d.f4phiXX
            - 10.0 * h**2 * ph**2 * d.tf4phi
        )
    )
    c3 = bphi / (2.0 * np.sqrt(f * h))

    # e3 = 1/2 partial_phi E_phi = 1/2[(partial_phi J_phi)' - partial_phi P_phi].
    jphi_bracket = (
        r**2 * (f * d.f2phiX + 2.0 * h * A**2 * d.f2phiY) * ph
        - 2.0 * h * A**2 * (2.0 * h * d.tf4phi + (3.0 * h - 2.0) * d.f4phiX) * ph
        + 2.0 * r * h**2 * A**2 * d.f3phiX * ph**2
        + h**3 * A**2 * d.f4phiXX * ph**3
        - 2.0 * r * h * A**2 * d.f3phi
    )
    partial_j = -np.sqrt(h / f) * jphi_bracket
    partial_p = 1.0 / np.sqrt(f * h) * (
        r**2 * f * d.f2phiphi
        + h
        * A**2
        * (
            4.0 * d.f4phiphi
            + 2.0 * h * (r * ph * d.f3phiphi - 2.0 * d.f4phiphi)
            + h**2 * (d.f4phiphiX + 2.0 * d.tf4phiphi) * ph**2
        )
    )
    e3 = 0.5 * (dr(partial_j) - partial_p)

    out = d[[c for c in ("u", "x") if c in d]].copy()
    out["v5"], out["c3"], out["e3"] = v5, c3, e3
    if not np.isfinite(out[["v5", "c3", "e3"]].to_numpy()).all():
        raise ValueError("nonfinite direct lower slots")
    return out, d


def audit_central_lower(root) -> dict:
    """Audit the historical selected Central lower slots against one total action."""
    from .central_action import central_action_inputs
    from .regional_coefficients import central_selected

    d = central_action_inputs(root).sort_values("x").reset_index(drop=True)
    direct, _ = emit_lower_slots(d)
    ref = central_selected(root).sort_values("x").reset_index(drop=True)
    selected = {"v5": ref.v5.to_numpy(), "c3": d.selected_c3.to_numpy(), "e3": d.selected_e3.to_numpy()}
    mask = (d.u.to_numpy() > 0.62) & (d.u.to_numpy() < 0.70)
    errors = {}
    for name in ("v5", "c3", "e3"):
        a, b = direct[name].to_numpy()[mask], selected[name][mask]
        errors[name] = {
            "max_scaled": float(np.max(np.abs(a - b) / np.maximum(1.0, np.abs(b)))),
            "median_scaled": float(np.median(np.abs(a - b) / np.maximum(1.0, np.abs(b)))),
        }
    compatible = all(errors[n]["max_scaled"] < 1e-5 for n in ("v5", "c3", "e3"))
    return {
        "status": "PASS" if compatible else "CENTRAL_SELECTED_LOWER_NOT_SINGLE_ACTION_REPLAY",
        "scope": "direct Appendix-A lower replay from one holonomic f2+f3+f4 action",
        "errors": errors,
        "selected_member_must_be_replaced_for_absolute_closure": not compatible,
    }
