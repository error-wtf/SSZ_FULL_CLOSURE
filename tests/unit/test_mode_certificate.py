"""MODE_CERTIFICATE_V1 axis evaluation: certified / rejected paths."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.qnm.mode_certificate import (  # noqa: E402
    THRESHOLDS,
    certify_mode,
)


def _all_good(**overrides):
    n = 20
    base = dict(
        L=6,
        mode_index=0,
        kinetic_norms=np.ones(n),
        tracking_overlaps=np.full(n, 0.99),
        relative_gaps=np.full(n, 0.1),
        omega_values=np.full(n, 1.0),
        observable_shares=np.full(n, 0.3),
        d2_value=0.95,
        scale_control_ok=True,
        basis_control_ok=True,
    )
    base.update(overrides)
    return certify_mode(**base)


def test_perfect_mode_is_certified():
    c = _all_good()
    assert c.verdict == "CERTIFIED_V1"
    assert c.failed_axes == []


def test_negative_kinetic_norm_rejects_A1():
    kn = np.ones(20)
    kn[5] = -0.1
    c = _all_good(kinetic_norms=kn)
    assert "A1" in c.failed_axes
    assert c.verdict == "REJECTED_V1"


def test_low_tracking_overlap_rejects_A2():
    ov = np.full(20, 0.99)
    ov[3] = 0.5
    c = _all_good(tracking_overlaps=ov)
    assert "A2" in c.failed_axes


def test_degenerate_gap_rejects_A3():
    gaps = np.full(20, 0.1)
    gaps[7] = 1e-5
    c = _all_good(relative_gaps=gaps)
    assert "A3" in c.failed_axes


def test_frequency_branch_jump_rejects_A4():
    # A4 now measures smoothness (a radial omega profile legitimately
    # varies; branch crossings jump).  A 31% single-step jump must fail.
    om = np.array([1.0] * 10 + [1.31] * 10)
    c = _all_good(omega_values=om)
    assert "A4" in c.failed_axes


def test_frequency_smooth_radial_profile_passes_A4():
    om = np.linspace(1.0, 1.2, 20)  # smooth 20% drift, tiny steps
    c = _all_good(omega_values=om)
    assert "A4" not in c.failed_axes


def test_invisible_mode_rejects_A5():
    c = _all_good(observable_shares=np.full(20, 1e-9))
    assert "A5" in c.failed_axes


def test_nan_d2_rejects_A6():
    c = _all_good(d2_value=float("nan"))
    assert "A6" in c.failed_axes


def test_failed_controls_reject_A7_A8():
    c = _all_good(scale_control_ok=False, basis_control_ok=False)
    assert "A7" in c.failed_axes and "A8" in c.failed_axes


def test_everything_failing_is_not_evaluable_when_all_skipped():
    n = 5
    c = _all_good(
        kinetic_norms=-np.ones(n),
        tracking_overlaps=np.zeros(n),
        relative_gaps=np.zeros(n),
        omega_values=np.array([np.nan] * n),
        observable_shares=np.zeros(n),
        d2_value=float("nan"),
        scale_control_ok=False,
        basis_control_ok=False,
        skip_axes=("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8"),
    )
    assert c.verdict == "NOT_EVALUABLE_V1"


def test_skip_axes_downgrade_rejection_to_not_evaluable():
    """A real failure outside skipped axes still rejects; skipped-only issues
    downgrade to NOT_EVALUABLE instead of rejecting."""
    n = 10
    skipped_only = certify_mode(
        L=6,
        mode_index=0,
        kinetic_norms=np.ones(n),
        tracking_overlaps=np.full(n, 0.99),
        relative_gaps=np.full(n, 0.1),
        omega_values=np.full(n, 1.0),
        observable_shares=np.full(n, 0.3),
        d2_value=float("nan"),
        scale_control_ok=True,
        basis_control_ok=True,
        skip_axes=("A6",),
    )
    assert skipped_only.verdict == "NOT_EVALUABLE_V1"
    assert skipped_only.failed_axes == []  # skipped axes are not hard failures
    assert skipped_only.axes["A6_skipped"]["pass"] is None


def _all_good(**overrides):
    n = 20
    base = dict(
        L=6,
        mode_index=0,
        kinetic_norms=np.ones(n),
        tracking_overlaps=np.full(n, 0.99),
        relative_gaps=np.full(n, 0.1),
        omega_values=np.full(n, 1.0),
        observable_shares=np.full(n, 0.3),
        d2_value=0.95,
        scale_control_ok=True,
        basis_control_ok=True,
    )
    base.update(overrides)
    return certify_mode(**base)


def test_threshold_budget_is_frozen_documentation():
    assert THRESHOLDS["min_tracking_overlap"] == 0.90
    assert THRESHOLDS["max_relative_frequency_span"] == 5e-3
