"""Holonomic f2 Hessian compatibility for the Inner action controls.

The six second derivatives of one smooth f2(phi, X, F) are not six independent
radial controls when the background value and first derivatives of f2 are held
fixed along the selected radial curve.  Writing

    H = d^2 f2 / d(phi, X, F)^2,
    t = (phi', X', F'),

background-nullity of the first on-curve derivatives requires

    H @ t = 0.

Thus the symmetric 3x3 Hessian has only three transverse degrees of freedom.
This module exposes that exact compatibility condition and the unique lower
mixed jets implied by a chosen transverse (XX, XF, FF) block whenever phi' is
nonzero.

This is deliberately an audit/inversion layer.  It does not promote any Inner
candidate or declare Direct-41 closure by itself.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..jets.jet9d8 import profile_derivative

HESSIAN_COLUMNS = (
    "f2phiphi_control",
    "f2Xphi_control",
    "f2phiF_control",
    "f2XX",
    "f2XF",
    "f2FF",
)
CHAIN_NAMES = ("phi_gradient", "X_gradient", "F_gradient")


def _phiprime(background: pd.DataFrame) -> np.ndarray:
    if "phiprime" in background:
        ph = background.phiprime.to_numpy(float)
    elif "X" in background and "h" in background:
        X = background.X.to_numpy(float)
        h = background.h.to_numpy(float)
        radicand = -2.0 * X / h
        if np.any(radicand < -1e-13):
            raise ValueError("cannot infer real phiprime from X,h")
        # Project branch: phi' < 0 throughout the Inner handover.
        ph = -np.sqrt(np.maximum(0.0, radicand))
    else:
        raise KeyError("phiprime or X,h required")
    if not np.isfinite(ph).all() or np.any(np.abs(ph) <= 1e-14):
        raise ValueError("holonomic Hessian audit requires finite nonzero phiprime")
    return ph


def electromagnetic_invariant(background: pd.DataFrame) -> np.ndarray:
    """Return F=-F_{mu nu}F^{mu nu}/4 for the static electric ansatz."""
    f = background.f.to_numpy(float)
    h = background.h.to_numpy(float)
    A = background.A0prime.to_numpy(float)
    if not all(np.isfinite(v).all() for v in (f, h, A)):
        raise ValueError("nonfinite electromagnetic background")
    if np.any(f <= 0) or np.any(h <= 0):
        raise ValueError("requires f>0 and h>0")
    return h * A**2 / (2.0 * f)


def tangent(background: pd.DataFrame, *, window: int = 9, degree: int = 8) -> np.ndarray:
    """Return the on-curve tangent (phi', X', F') on the radial x grid."""
    x = background.x.to_numpy(float)
    X = background.X.to_numpy(float)
    ph = _phiprime(background)
    F = electromagnetic_invariant(background)
    Xp = profile_derivative(x, X, 1, window, degree)
    Fp = profile_derivative(x, F, 1, window, degree)
    t = np.column_stack([ph, Xp, Fp])
    if not np.isfinite(t).all():
        raise ValueError("nonfinite Hessian tangent")
    return t


def assemble_hessian(lower: pd.DataFrame, principal: pd.DataFrame) -> np.ndarray:
    """Assemble H in coordinate order (phi, X, F) from existing control frames."""
    if len(lower) != len(principal):
        raise ValueError("lower/principal controls have different lengths")
    if (
        "x" in lower
        and "x" in principal
        and not np.array_equal(lower.x.to_numpy(float), principal.x.to_numpy(float))
    ):
        # Exact identity is intentional: these are controls on one radial section.
        raise ValueError("lower/principal controls are not on the identical x grid")
    for name in ("f2phiphi_control", "f2Xphi_control", "f2phiF_control"):
        if name not in lower:
            raise KeyError(name)
    for name in ("f2XX", "f2XF", "f2FF"):
        if name not in principal:
            raise KeyError(name)

    H = np.zeros((len(lower), 3, 3), dtype=float)
    H[:, 0, 0] = lower.f2phiphi_control.to_numpy(float)
    H[:, 0, 1] = H[:, 1, 0] = lower.f2Xphi_control.to_numpy(float)
    H[:, 0, 2] = H[:, 2, 0] = lower.f2phiF_control.to_numpy(float)
    H[:, 1, 1] = principal.f2XX.to_numpy(float)
    H[:, 1, 2] = H[:, 2, 1] = principal.f2XF.to_numpy(float)
    H[:, 2, 2] = principal.f2FF.to_numpy(float)
    if not np.isfinite(H).all():
        raise ValueError("nonfinite Hessian controls")
    return H


def chain_residuals(
    background: pd.DataFrame,
    lower: pd.DataFrame,
    principal: pd.DataFrame,
    *,
    window: int = 9,
    degree: int = 8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Ht, a termwise-normalized residual, and the curve tangent.

    The normalized residual is |sum terms|/sum|terms| per chain-rule row and is
    therefore in [0,1] up to floating-point roundoff.  A value near one means
    essentially no cancellation toward the required zero.
    """
    H = assemble_hessian(lower, principal)
    t = tangent(background, window=window, degree=degree)
    residual = np.einsum("nij,nj->ni", H, t)
    terms = H * t[:, None, :]
    denom = np.sum(np.abs(terms), axis=2)
    normalized = np.divide(
        np.abs(residual),
        denom,
        out=np.zeros_like(residual),
        where=denom > 0,
    )
    return residual, normalized, t


def implied_lower_from_transverse(
    background: pd.DataFrame,
    principal: pd.DataFrame,
    *,
    window: int = 9,
    degree: int = 8,
) -> pd.DataFrame:
    """Recover the unique mixed/phi Hessian jets implied by XX,XF,FF and Ht=0.

    With phi' != 0 the last two rows of Ht=0 determine phiX and phiF, and the
    first row then determines phiphi.  This exposes the actual three-degree-of-
    freedom background-null Hessian section.
    """
    t = tangent(background, window=window, degree=degree)
    ph, Xp, Fp = t.T
    XX = principal.f2XX.to_numpy(float)
    XF = principal.f2XF.to_numpy(float)
    FF = principal.f2FF.to_numpy(float)

    phiX = -(Xp * XX + Fp * XF) / ph
    phiF = -(Xp * XF + Fp * FF) / ph
    phiphi = -(Xp * phiX + Fp * phiF) / ph

    out = background[[c for c in ("u", "x") if c in background]].copy()
    out["f2phiphi_control"] = phiphi
    out["f2Xphi_control"] = phiX
    out["f2phiF_control"] = phiF
    out["f2XX"] = XX
    out["f2XF"] = XF
    out["f2FF"] = FF
    return out


def split_local_response(background: pd.DataFrame, H: np.ndarray) -> np.ndarray:
    """Replay the historical *split* local maps [v5,c3,e3,v1,v4,c2].

    This function exists only to state exactly which six coefficient deltas the
    previous independent inverses were trying to realize.  It is not a valid
    common-action emitter: in particular it omits the Appendix-A c3 cross term
    and the derivative contribution to e3.
    """
    r = background.x.to_numpy(float)
    f = background.f.to_numpy(float)
    h = background.h.to_numpy(float)
    A = background.A0prime.to_numpy(float)
    ph = _phiprime(background)
    if H.shape != (len(background), 3, 3):
        raise ValueError("H shape mismatch")

    phiphi = H[:, 0, 0]
    phiX = H[:, 0, 1]
    phiF = H[:, 0, 2]
    XX = H[:, 1, 1]
    XF = H[:, 1, 2]
    FF = H[:, 2, 2]

    v5 = r**2 * np.sqrt(h / f) * A * phiF
    c3 = 0.5 * r**2 * np.sqrt(f * h) * ph**2 * phiX
    e3 = -0.5 * r**2 * np.sqrt(f / h) * phiphi
    v1 = r**2 * h**1.5 * A**2 / (2 * f**1.5) * FF
    v4 = -(h**1.5) * A * r**2 * ph / np.sqrt(f) * XF
    c2 = (
        h**1.5 * A**2 * r**2 * ph / (2 * np.sqrt(f)) * XF
        - 0.5 * r**2 * np.sqrt(f * h) * h * ph**3 * XX
    )
    return np.column_stack([v5, c3, e3, v1, v4, c2])


def response_from_hessian(
    background: pd.DataFrame,
    H: np.ndarray,
    *,
    window: int = 9,
    degree: int = 8,
) -> np.ndarray:
    """Appendix-A f2-only response [v5,c3,e3,v1,v4,c2] of one Hessian.

    For a background-null deformation delta f2(phi,X,F), value and first action
    jets vanish on the selected curve and the Hessian obeys H t = 0.  Zhang--
    Kase Appendix A gives

      v5 = sqrt(h/f) A0' r^2 f2_{phi F},
      c3 = r^2/(2 sqrt(fh)) [f h phi'^2 f2_{phi X}
                              - h A0'^2 f2_{phi F}],

    while e3 = (1/2) partial_phi E_phi.  Using E_phi=J_phi'-P_phi and Ht=0,
    the f2-only contribution reduces to

      e3 = -1/2 d_r[r^2 sqrt(fh) phi' f2_{phi X}]
           -1/2 r^2 sqrt(f/h) f2_{phi phi}.

    The principal responses are the exact Appendix-A f2 Hessian terms already
    used by ``principal_controls.raw_response``.
    """
    r = background.x.to_numpy(float)
    f = background.f.to_numpy(float)
    h = background.h.to_numpy(float)
    A = background.A0prime.to_numpy(float)
    ph = _phiprime(background)
    if H.shape != (len(background), 3, 3):
        raise ValueError("H shape mismatch")

    phiphi = H[:, 0, 0]
    phiX = H[:, 0, 1]
    phiF = H[:, 0, 2]
    XX = H[:, 1, 1]
    XF = H[:, 1, 2]
    FF = H[:, 2, 2]

    v5 = r**2 * np.sqrt(h / f) * A * phiF
    c3 = 0.5 * r**2 * np.sqrt(f * h) * ph**2 * phiX - 0.5 * r**2 * A**2 * np.sqrt(h / f) * phiF
    current_phi = r**2 * np.sqrt(f * h) * ph * phiX
    e3 = (
        -0.5 * profile_derivative(r, current_phi, 1, window, degree)
        - 0.5 * r**2 * np.sqrt(f / h) * phiphi
    )
    v1 = r**2 * h**1.5 * A**2 / (2 * f**1.5) * FF
    v4 = -(h**1.5) * A * r**2 * ph / np.sqrt(f) * XF
    c2 = (
        h**1.5 * A**2 * r**2 * ph / (2 * np.sqrt(f)) * XF
        - 0.5 * r**2 * np.sqrt(f * h) * h * ph**3 * XX
    )
    return np.column_stack([v5, c3, e3, v1, v4, c2])


def audit_existing_split_controls(
    background: pd.DataFrame,
    lower: pd.DataFrame,
    principal: pd.DataFrame,
    *,
    window: int = 9,
    degree: int = 8,
    tolerance: float = 1e-9,
) -> tuple[pd.DataFrame, dict]:
    """Audit whether the existing split controls can be one background-null H."""
    residual, normalized, t = chain_residuals(
        background, lower, principal, window=window, degree=degree
    )
    H = assemble_hessian(lower, principal)
    implied = implied_lower_from_transverse(background, principal, window=window, degree=degree)
    H_implied = assemble_hessian(implied, principal)
    implied_residual = np.einsum("nij,nj->ni", H_implied, t)

    # ``current_response`` is the six-delta target actually requested by the
    # historical split inverses.  ``implied_response`` is what the same
    # transverse principal block produces after enforcing one smooth f2.
    current_response = split_local_response(background, H)
    implied_response = response_from_hessian(background, H_implied, window=window, degree=degree)
    response_error = np.abs(implied_response - current_response) / np.maximum(
        1.0, np.abs(current_response)
    )

    out = background[[c for c in ("u", "x") if c in background]].copy()
    out[["phi_prime", "X_prime", "F_prime"]] = t
    for j, name in enumerate(CHAIN_NAMES):
        out[f"chain_{name}"] = residual[:, j]
        out[f"chain_{name}_normalized"] = normalized[:, j]
    for name in ("f2phiphi_control", "f2Xphi_control", "f2phiF_control"):
        out[f"current_{name}"] = lower[name].to_numpy(float)
        out[f"implied_{name}"] = implied[name].to_numpy(float)
    response_names = ("v5", "c3", "e3", "v1", "v4", "c2")
    for j, name in enumerate(response_names):
        out[f"current_response_{name}"] = current_response[:, j]
        out[f"holonomic_response_{name}"] = implied_response[:, j]
        out[f"response_{name}_scaled_change"] = response_error[:, j]

    first_row_max = float(np.max(normalized[:, 0]))
    max_norm = float(np.max(normalized))
    report = {
        "status": "PASS" if max_norm <= tolerance else "FAIL",
        "criterion": "symmetric f2(phi,X,F) Hessian must satisfy H @ (phi',X',F') = 0",
        "tolerance": tolerance,
        "rows": len(background),
        "max_abs_chain_residual": {
            name: float(np.max(np.abs(residual[:, j]))) for j, name in enumerate(CHAIN_NAMES)
        },
        "median_abs_chain_residual": {
            name: float(np.median(np.abs(residual[:, j]))) for j, name in enumerate(CHAIN_NAMES)
        },
        "max_normalized_chain_residual": {
            name: float(np.max(normalized[:, j])) for j, name in enumerate(CHAIN_NAMES)
        },
        "median_normalized_chain_residual": {
            name: float(np.median(normalized[:, j])) for j, name in enumerate(CHAIN_NAMES)
        },
        "first_chain_rule_lower_only": True,
        "first_chain_rule_max_normalized_residual": first_row_max,
        "principal_controls_can_repair_first_chain_rule": False,
        "max_implied_holonomic_chain_abs": float(np.max(np.abs(implied_residual))),
        "max_scaled_response_change_if_principal_block_is_kept": {
            name: float(np.max(response_error[:, j])) for j, name in enumerate(response_names)
        },
        "conclusion": (
            "existing split lower/principal controls are not one background-null smooth f2 Hessian"
            if max_norm > tolerance
            else "existing split controls satisfy the common Hessian chain rules"
        ),
        "closure_effect": (
            "INNER_DIRECT_41 remains NOT_CERTIFIED until targets are redesigned "
            "on the 3-DOF holonomic Hessian manifold"
        ),
    }
    return out, report
