from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / "data/generated/strong_field_transition_2026-09-20"


def _load(name: str):
    return json.loads((DIR / name).read_text())


def test_projected_joint_action_tangent_is_locally_controllable():
    r = _load("TRANSITION_PROJECTED_CONTROLLABILITY.json")
    assert r["status"] == "PASS_LOCAL_TANGENT_CONTROLLABILITY"
    assert r["global_summary"]["metric_vector_tangent_has_one_direction_improving_all_required_L_at_all_sampled_points"]
    assert r["global_summary"]["extended_scalar_identity_tangent_has_one_direction_improving_all_required_L_at_all_sampled_points"]
    for point in r["sample_points"]:
        assert point["simultaneous_all_L"]["with_scalar_identity_tangent"]["all_required_L_gain_positive"]


def test_projected_audit_never_uses_historical_inner_as_absolute_seed():
    r = _load("TRANSITION_PROJECTED_CONTROLLABILITY.json")
    right = r["requested_right_seed_diagnostic"]
    assert right["u_requested"] == 0.712
    assert right["status"] == "NOT_EVALUATED_NO_CERTIFIABLE_ABSOLUTE_ACTION_SEED"


def test_horndeski_a1_survives_onshell_projection_but_is_not_promoted():
    r = _load("TRANSITION_HORNDESKI_A1_PROJECTION.json")
    assert r["status"] == "PASS_A1_SURVIVES_BACKGROUND_NULL_PROJECTION"
    lr = min(r["sample_points"], key=lambda x: abs(x["u_grid"] - 0.7061345733))
    assert lr["projected_fraction"] > 0.5
    assert lr["background_residual_norm"] < 1e-10
    assert not r["promotion"]["a1_bundle_is_finite_transition_solution"]


def test_projection_numerics_are_converged_and_a1_width_warning_is_preserved():
    r = _load("TRANSITION_PROJECTED_CONTROLLABILITY_CONVERGENCE.json")
    assert r["status"] == "PASS_PROJECTED_CONTROLLABILITY_NUMERICS"
    assert r["joint_H_SVT_projected_tangent"]["fd_converged"]
    assert r["joint_H_SVT_projected_tangent"]["relative_spread_scalar_identity_common_gain"] < 1e-4
    assert r["horndeski_a1_null_bundle"]["amplitude_derivative_sign_stable"]
    assert not r["horndeski_a1_null_bundle"]["all_L_sign_robust_to_width"]
    assert not r["promotion"]["a1_alone_promoted_to_transition_solution"]
