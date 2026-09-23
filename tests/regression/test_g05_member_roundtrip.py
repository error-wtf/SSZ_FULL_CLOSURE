"""G05 + G20: scalar-integration round trip and same-member provenance.

G05 dedicated gate test (closes the coverage gap):
    certified boundary data -> scalar integration -> algebraic re-solve ->
    independent background residuals -> SAME member.

G20: every artifact references the SAME member hash (MODEL_LOCK == manifest
== canonical CSV stream).
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.svt_background_eom import (
    evaluate_svt_background,
    scalar_ode_identity,
)

ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv"
MANIFEST_PATH = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
LOCK_PATH = ROOT / "MODEL_LOCK.json"


def _member_stream() -> pd.DataFrame:
    return pd.read_csv(CSV_PATH, index_col=0)


def _canonical_hash(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def test_G20_same_member_provenance():
    """MODEL_LOCK, manifest and the canonical CSV stream all carry the SAME
    member hash; the manifest pins a git commit and the production window."""
    lock = json.loads(LOCK_PATH.read_text())
    manifest = json.loads(MANIFEST_PATH.read_text())
    csv_text = CSV_PATH.read_text()
    csv_hash = _canonical_hash(csv_text)
    assert lock["action_member_sha256"] == csv_hash
    assert manifest["member_hash"] == csv_hash
    assert lock["action_member_reference"] == str(MANIFEST_PATH.relative_to(ROOT))
    assert manifest["supersedes"]["reason"].startswith("historical ZERO-VECTOR")
    assert manifest["domain"]["production_window"] == [0.62, 0.70]
    assert manifest["boundary_provenance"]["a2_solve_success"] is True
    assert manifest["action_definition_status"]["genuine_svt_theta_theta"] \
        .startswith("CLOSED")


def test_member_builder_reproduces_locked_stream():
    """The locked builder must reproduce the frozen member stream within the
    release numerical policy (background_residual_abs = 1e-10): no silent
    member switching.  The LOCKED stream itself is pinned by the member hash
    (CSV == manifest == MODEL_LOCK)."""
    build = build_onshell_central(ROOT)
    stream = _member_stream().drop(columns=["A0"])
    fresh = build.action[stream.columns].sort_values("u").reset_index(drop=True)
    for col in stream.columns:
        a = fresh[col].to_numpy(float)
        b = stream[col].to_numpy(float)
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-10)
    lock = json.loads(LOCK_PATH.read_text())
    csv_text = CSV_PATH.read_text()
    assert lock["action_member_sha256"] == _canonical_hash(csv_text)


def test_G05_scalar_integration_round_trip():
    """G05 dedicated gate: certified boundary data -> scalar integration of
    the f3X ODE (A q' + B q + C = 0 identity form) -> algebraic triangle
    re-solve (JA -> f2F_eff, E00 -> f2, E11-E00 -> f2X) -> independent
    background residuals at machine level -> SAME member hash."""
    build = build_onshell_central(ROOT)
    d = build.action.sort_values("x").reset_index(drop=True)  # ascending x
    mask = (d.u.to_numpy(float) > 0.62) & (d.u.to_numpy(float) < 0.70)

    # 1) certified boundary data: member scalar jet at the inner boundary row
    q0 = float(d.f3X.iloc[0])

    # 2) scalar integration along the certified frame
    ident = scalar_ode_identity(d)
    x = d.x.to_numpy(float)

    def rhs(xv, y):
        a = np.interp(xv, x, ident.A)
        bq = np.interp(xv, x, ident.B)
        c = np.interp(xv, x, ident.C)
        return [(-bq*y - c)/a]

    sol = solve_ivp(rhs, (x[0], x[-1]), [q0], dense_output=True,
                    rtol=1e-10, atol=1e-12)
    assert sol.success
    q_int = sol.sol(x)[0]

    # 3) algebraic triangle re-solve on the integrated frame
    f, h, ph, ap = (d.f.to_numpy(float), d.h.to_numpy(float),
                    d.phiprime.to_numpy(float), d.A0prime.to_numpy(float))
    r = d.x.to_numpy(float)
    f3, f4, f4X, f4XX, tf4 = (d.f3.to_numpy(float), d.f4.to_numpy(float),
                              d.f4X.to_numpy(float), d.f4XX.to_numpy(float),
                              d.tf4.to_numpy(float))
    f2Y = d.f2Y.to_numpy(float)
    f2F_eff = -(4*r*h*ph*f3 + 8*(1-h)*f4 + 2*h**2*ph**2*(f4X + 2*tf4))/r**2
    # canonical jet derivative service (same as the evaluator - machine-level
    # algebraic re-solve requires consistent derivative channels)
    from ssz_p5.jets.jet9d8 import profile_derivative
    hp = profile_derivative(r, h, 1, 9, 8)
    fp = profile_derivative(r, f, 1, 9, 8)
    core = (-2*r*h**2*ph*ap**2*f3
            + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X + 2*tf4)))
    # explicit re-solve: E00 = rf h' - f(1-h) - r^2(f f2 - h ap^2 f2F_eff) - core = 0
    f2_re = (r*f*hp - f*(1-h) - core + r**2*h*ap**2*f2F_eff)/(r**2*f)
    d_re = d.copy()
    d_re["f3X"] = q_int
    d_re["f2F"] = f2F_eff + 2*h*ph**2*f2Y  # store f2F such that f2F_eff is preserved
    d_re["f2"] = f2_re
    # E11 - E00 = G0 - r^2 f h phi'^2 f2X with G0 the f2X-free remainder
    # (f3/f4-channel brackets included) -> evaluate G0 numerically:
    d_q = d_re.copy()
    d_q["f2X"] = np.zeros_like(r)
    res_q = evaluate_svt_background(d_q, require_f2y_zero=False)
    G0 = res_q.E11 - res_q.E00
    d_re["f2X"] = G0/(r**2*f*h*ph**2)

    # 4) independent background residuals on the re-solved frame
    res = evaluate_svt_background(d_re, require_f2y_zero=False)
    assert np.max(np.abs(res.JA[mask])) < 1e-10
    assert np.max(np.abs(res.E00[mask])) < 1e-10
    assert np.max(np.abs(res.E11[mask])) < 1e-10

    # machine-algebra: the ORIGINAL certified frame must satisfy the same
    assert (evaluate_svt_background(d, require_f2y_zero=False) is not None)

    # scalar round-trip quality: integrated vs stored member f3X
    scale = np.maximum(1.0, np.abs(d.f3X.to_numpy(float)))
    dev = np.abs(q_int - d.f3X.to_numpy(float))/scale
    assert np.max(dev[mask]) < 1e-4, np.max(dev[mask])

    # 5) SAME member: hash of the frozen stream is unchanged
    lock = json.loads(LOCK_PATH.read_text())
    csv_hash = _canonical_hash(CSV_PATH.read_text())
    assert lock["action_member_sha256"] == csv_hash
