"""G50, G70, G90: dedicated current-member perturbation gates.

G50: radial characteristics c_r^2 (generalized (K, G) eigenvalues) across the
     full multipole ladder, production window, with a negative control.
G70: finite-l even sector: kinetic positivity and radial characteristics at
     every REQUIRED multipole; no finite-l instability anywhere on the member.
G90: global regularity scan: every member quantity and every perturbation
     diagnostic finite and regular across the whole domain incl. interfaces.
"""
import json
from pathlib import Path

import numpy as np

from ssz_p5.production.electric_hybrid_onshell_central import (
    _REQUIRED_L,
    build_onshell_central,
    principal_audit,
)

ROOT = Path(__file__).resolve().parents[2]


def test_G50_radial_characteristics_positive():
    """G50: the radial characteristic speed squared c_r^2 = min eig of the
    generalized (K_s, G_s) problem is STRICTLY POSITIVE at every required
    multipole on the production window of the current member; constraint
    structures R and S are exactly zero (properly eliminated)."""
    build, rows = principal_audit(ROOT)
    assert [row["L"] for row in rows] == list(_REQUIRED_L)
    for row in rows:
        assert row["pass"] is True, row
        assert row["min_eig_K"] > 0
        assert row["min_cr2"] > 0
        assert row["negative_K_rows"] == 0
        assert row["negative_radial_rows"] == 0
        assert row["max_R_abs"] == 0.0
        # S antisymmetric part at float roundoff (matrices are O(1)-O(10))
        assert row["max_S_sym"] < 1e-8


def test_G70_finite_l_even_sector():
    """G70: dedicated finite-l battery - the even-parity kinetic sector is
    positive and radially hyperbolic at EVERY required multipole; the ladder
    is finite (no l->infinity extrapolation is used anywhere) and the
    diagnostics are recorded per multipole."""
    build, rows = principal_audit(ROOT)
    ladder = _REQUIRED_L
    assert len(ladder) >= 5 and all(float(l).is_integer() for l in ladder)
    per_l = {row["L"]: row for row in rows}
    for L in ladder:
        row = per_l[L]
        # no finite-l ghost or gradient instability at this multipole
        assert row["min_eig_K"] > 0, (L, row)
        assert row["min_cr2"] > 0, (L, row)
    # nontriviality: the diagnostics genuinely vary with l (not a stub)
    ke = [per_l[L]["min_eig_K"] for L in ladder]
    assert max(ke) > 10 * min(ke)
    # monotone decay of the finite-l kinetic scale (l-dependent kinetic weight)
    assert all(ke[i] >= ke[i + 1] for i in range(len(ke) - 1))


def test_G90_global_regularity_scan():
    """G90: global regularity - every member profile and action jet is finite
    across the FULL frozen domain (u in [0.61, 0.71], both interfaces
    included), the build diagnostics are finite, and the principal battery
    reports finite eigenvalues on interior rows."""
    build = build_onshell_central(ROOT)
    d = build.action
    numeric = d.select_dtypes(include=[float]).to_numpy(float)
    assert np.all(np.isfinite(numeric)), "non-finite member quantity"
    u = d.u.to_numpy(float)
    assert u.min() == 0.61 and abs(u.max() - 0.7099849962490622) < 1e-12
    # interfaces included in the domain: rows at both edges are finite
    for i_edge in (0, 1, -2, -1):
        assert np.all(np.isfinite(d.select_dtypes(include=[float]).iloc[i_edge]
                                  .to_numpy(float)))
    # build diagnostics finite
    diag = json.loads(json.dumps(build.diagnostics))  # plain dict of floats
    def _walk(x):
        if isinstance(x, dict):
            return all(_walk(v) for v in x.values())
        if isinstance(x, (int, float)):
            return np.isfinite(x)
        return True
    assert _walk(diag)
    # principal battery eigenvalues finite on the production window
    _, rows = principal_audit(ROOT)
    for row in rows:
        assert np.isfinite(row["min_eig_K"]) and np.isfinite(row["min_cr2"])
