#!/usr/bin/env python3
"""Diagnostic only: isolate the Electric-Hybrid scalar/metric coupling.

This deliberately does *not* modify the production member.  It solves the
linear scalar identity A q' + B q + C = 0 for q=f3X while freezing the other
background controls, reimposes the holonomic f3_phi chain rule, and then checks
all background equations.  A scalar-only solve is promotable only if E00, E11,
JA and Ephi all stay within the release tolerance; the current checkpoint is
expected to demonstrate why a joint scalar+metric action solve is still needed.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.jets.jet9d8 import profile_derivative
from ssz_p5.policy import numerical_policy
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.svt_background_eom import (
    evaluate_svt_background,
    residual_metrics,
    scalar_ode_identity,
)


def _metrics(action):
    mask = (action.u.to_numpy(float) > 0.62) & (action.u.to_numpy(float) < 0.70)
    vals = evaluate_svt_background(action, window=9, degree=8)
    return residual_metrics(vals, mask), mask


def diagnostic(root: Path = ROOT) -> dict:
    policy = numerical_policy()
    base = build_onshell_central(root).action.copy()
    ident = scalar_ode_identity(base, window=9, degree=8)

    x = base.x.to_numpy(float)
    u = base.u.to_numpy(float)
    q0 = base.f3X.to_numpy(float)
    band = (u >= 0.61) & (u <= 0.71)
    ib = np.flatnonzero(band)
    if len(ib) < 20:
        raise RuntimeError("Electric-Hybrid scalar diagnostic band is empty")

    xb = x[ib]
    A = PchipInterpolator(xb, ident.A[ib], extrapolate=False)
    B = PchipInterpolator(xb, ident.B[ib], extrapolate=False)
    C = PchipInterpolator(xb, ident.C[ib], extrapolate=False)

    i61 = int(np.argmin(np.abs(u - 0.61)))
    i71 = int(np.argmin(np.abs(u - 0.71)))

    def rhs(xx, yy):
        aa = float(A(xx))
        if abs(aa) < 1e-14:
            raise RuntimeError(f"scalar identity A coefficient too small at x={xx}")
        return [-(float(B(xx)) * float(yy[0]) + float(C(xx))) / aa]

    sol = solve_ivp(
        rhs,
        (float(x[i61]), float(x[i71])),
        [float(q0[i61])],
        rtol=2e-11,
        atol=2e-12,
        dense_output=True,
        max_step=2e-4,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    trial = base.copy()
    q = q0.copy()
    q[ib] = sol.sol(x[ib])[0]
    trial["f3X"] = q

    # Keep f3 itself fixed but restore the on-curve chain rule after changing f3X.
    f3r = profile_derivative(x, trial.f3.to_numpy(float), 1, 9, 8)
    Xr = profile_derivative(x, trial.X.to_numpy(float), 1, 9, 8)
    ph = trial.phiprime.to_numpy(float)
    trial["f3phi"] = (f3r - q * Xr) / ph

    baseline_metrics, prodmask = _metrics(base)
    trial_metrics, _ = _metrics(trial)

    qprime = q0.copy()
    qprime[ib] = np.array([rhs(float(xx), [float(qq)])[0] for xx, qq in zip(x[ib], q[ib])])
    ident_trial = scalar_ode_identity(trial, q_prime=qprime, window=9, degree=8)
    scalar_identity_max = float(np.max(np.abs(ident_trial.residual[prodmask])))

    tol = float(policy["background_residual_abs"])
    simultaneous = {
        "E00": trial_metrics["E00"]["max_abs"] < tol,
        "E11": trial_metrics["E11"]["max_abs"] < tol,
        "JA": trial_metrics["JA"]["max_abs"] < tol,
        "Ephi_identity": scalar_identity_max < tol,
    }
    eligible = bool(all(simultaneous.values()))

    report = {
        "status": "PASS_JOINT_BACKGROUND" if eligible else "DIAGNOSTIC_COUPLING_OPEN",
        "scope": "non-production scalar-only f3X identity solve with all other action controls frozen",
        "production_member_modified": False,
        "promotion_eligible": eligible,
        "release_tolerance": tol,
        "solve": {
            "unknown": "q=f3X",
            "identity": "A(r) q'(r) + B(r) q(r) + C(r) = 0",
            "initial_u": float(u[i61]),
            "terminal_u": float(u[i71]),
            "success": bool(sol.success),
            "nfev": int(sol.nfev),
            "max_abs_q_shift_production": float(np.max(np.abs((q - q0)[prodmask]))),
            "scalar_identity_max_abs_production": scalar_identity_max,
        },
        "baseline_background": baseline_metrics,
        "scalar_only_trial_background": trial_metrics,
        "simultaneous_strict_gates": simultaneous,
        "interpretation": (
            "Scalar identity can be solved at strict precision, but it may not be promoted unless "
            "the same action member simultaneously preserves E00, E11 and JA. A failure of E11 here "
            "is evidence for a joint scalar+metric action solve, not permission to loosen tolerances."
        ),
    }
    return report


def main() -> int:
    out = ROOT / "data/generated/absolute_attempt_2026-09-20/ELECTRIC_HYBRID_SCALAR_COUPLING_DIAGNOSTIC.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    report = diagnostic(ROOT)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"WROTE {out}")
    # Diagnostic succeeds when it has made a definite, finite determination.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
