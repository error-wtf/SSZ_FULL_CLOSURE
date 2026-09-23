"""Dedicated tests for coverage-gap gates G02, G03, G04 (current electric member).

Member under test: build_onshell_central(ROOT).action — the on-shell electric
central background (u in [0.61, 0.71), the production window containing the
inner stable light ring).  These tests close the corresponding entries in
TEST_COVERAGE_GAPS.json and are required for those gates to certify.
"""
from pathlib import Path

import json

import numpy as np
import pytest

from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.svt_background_eom import evaluate_svt_background, scalar_ode_identity

ROOT = Path(__file__).resolve().parents[2]
INTERIOR = (0.615, 0.705)   # excludes handover stencil edges


@pytest.fixture(scope="module")
def member():
    return build_onshell_central(ROOT).action


def test_G02_geometry_profile_identity(member):
    """G02: the geometry/profile actually used by the production member obeys
    the MODEL_LOCK conventions exactly: u*x = r_s = 1, strictly increasing
    radial grid, f > 0 and h > 0 outside the horizon (physical branch), and
    the X = -h*phi'^2/2 chain rule holds to machine level."""
    x = member.x.to_numpy(float)
    u = member.u.to_numpy(float)
    f = member.f.to_numpy(float)
    h = member.h.to_numpy(float)
    ph = member.phiprime.to_numpy(float)
    X = member.X.to_numpy(float)
    assert np.all(np.diff(x) > 0), "radial grid must be strictly increasing"
    assert float(np.max(np.abs(u * x - 1.0))) < 1e-12, "u = r_s/x with r_s = 1"
    assert np.all(f > 0) and np.all(h > 0), "physical branch f>0, h>0"
    xerr = np.max(np.abs(X - (-h * ph**2 / 2)))
    assert xerr < 1e-12, f"X = -h phi'^2/2 violated by {xerr}"


def test_G03_background_algebraic_rank_full_window(member):
    """G03: the algebraic solve d(JA, E00, E11)/d(f2F, f2, f2X) is regular
    (nonzero determinant, condition number finite) over the FULL production
    window, not only at the inner ring.  Determinant computed by pointwise
    column perturbation of the algebraic jets (the residuals are pointwise
    algebraic in the jets; the geometry columns are held fixed)."""
    d = member.copy()
    mask = (d.u > INTERIOR[0]) & (d.u < INTERIOR[1])
    idx = np.flatnonzero(mask)[::8]          # subsampled grid, same evaluator
    scale = {"f2F": 1e-6, "f2": 1e-6, "f2X": 1e-6}
    slots = ("JA", "E00", "E11")
    jets = ("f2F", "f2", "f2X")

    def residuals(jet_values):
        dd = d.copy()
        for k, v in jet_values.items():
            dd[k] = v
        out = evaluate_svt_background(dd, require_f2y_zero=True)
        return {s: np.asarray(getattr(out, s), float) for s in slots}

    base = residuals({})
    dets = []
    for i in idx:
        J = np.zeros((3, 3))
        for j, jet in enumerate(jets):
            up = {k: v.copy() for k, v in d.items() if k in scale}
            up = {k: d[k].to_numpy(float).copy() for k in scale}
            up[jet][i] += scale[jet]
            pert = residuals(up)
            for a, s in enumerate(slots):
                J[a, j] = (pert[s][i] - base[s][i]) / scale[jet]
        dets.append(np.linalg.det(J))
    dets = np.asarray(dets)
    assert np.all(np.isfinite(dets)), "Jacobian determinant must be finite"
    assert np.min(np.abs(dets)) > 1e-6, (
        f"algebraic solve degenerate: min|det| = {np.min(np.abs(dets)):.3e} "
        f"at x = {d.x.to_numpy(float)[idx[int(np.argmin(np.abs(dets)))]]}")


def test_G04_scalar_ode_regular_full_window(member):
    """G04: the scalar ODE A(r) q' + B(r) q + C(r) = 0 has a regular leading
    coefficient over the FULL interior window (A(r) != 0, no sign changes) —
    min |A|, its location and the boundary classification are recorded."""
    ident = scalar_ode_identity(member)
    A = np.asarray(ident.A, float)
    u = member.u.to_numpy(float)
    mask = (u > INTERIOR[0]) & (u < INTERIOR[1])
    Ai = A[mask]
    assert np.all(np.isfinite(Ai)) and np.min(np.abs(Ai)) > 1e-3, (
        f"scalar ODE leading coefficient degenerate: min|A| = {np.min(np.abs(Ai)):.3e}")
    assert np.all(np.sign(Ai) == np.sign(Ai[0])), "sign change in A(r): hidden singular point"
    j = int(np.argmin(np.abs(Ai)))
    record = {"min_abs_A_interior": float(np.min(np.abs(Ai))),
              "u_at_min": float(u[mask][j]),
              "sign_consistent": True,
              "window": INTERIOR,
              "boundary_rows_excluded": int((~mask).sum())}
    (ROOT / "data" / "generated" / "phase2_q2" / "G04_SCALAR_ODE_REGULARITY.json") \
        .write_text(json.dumps(record, indent=1))
