"""Lock the detected failure, not a false background PASS."""

from pathlib import Path

from ssz_p5.action.onshell import audit_frozen_carrier


def test_frozen_carrier_violation_is_resolved_above_derivative_sensitivity():
    report = audit_frozen_carrier(Path(__file__).resolve().parents[2])
    assert report["status"] == "FAIL_ON_SHELL_IDENTITY"
    assert len(report["routes"]) == 5
    assert report["route_spread"] < 1e-7
    assert 0.07 < min(row["residual"] for row in report["routes"]) < 0.072
    assert abs(report["routes"][-1]["independent_subtraction_residual"]) < 1e-12
    assert report["consistency_checks"]["max_abs_H_minus_2G4"] < 1e-12
