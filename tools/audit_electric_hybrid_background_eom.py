#!/usr/bin/env python3
"""Audit the current Electric-Hybrid Central against the full static SVT EOM."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.policy import numerical_policy
from ssz_p5.production.central_action import central_action_inputs
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central
from ssz_p5.production.sources import SOURCE_REGISTRY
from ssz_p5.production.svt_background_eom import evaluate_svt_background, residual_metrics


SETTINGS = ((7, 6), (9, 8), (11, 8), (13, 8))


def audit(root: Path = ROOT):
    policy = numerical_policy()
    current = build_onshell_central(root).action.copy()
    current_mask = (current.u.to_numpy(float) > 0.62) & (current.u.to_numpy(float) < 0.70)

    reference = central_action_inputs(root)
    reference_mask = (reference.u.to_numpy(float) > 0.62) & (reference.u.to_numpy(float) < 0.70)
    stored = pd.read_csv(root / SOURCE_REGISTRY["central_exact_SVT"]["exact_regression"])
    stored_mask = (stored.u.to_numpy(float) > 0.62) & (stored.u.to_numpy(float) < 0.70)

    convergence = []
    for window, degree in SETTINGS:
        cv = evaluate_svt_background(current, window=window, degree=degree)
        rv = evaluate_svt_background(reference, window=window, degree=degree)
        cm = residual_metrics(cv, current_mask)
        rm = residual_metrics(rv, reference_mask)
        convergence.append(
            {
                "window": window,
                "degree": degree,
                "current": cm,
                "archived_exact_reference_direct_replay": rm,
                "scalar_floor_ratio_max": cm["Ephi_direct"]["max_abs"]
                / max(rm["Ephi_direct"]["max_abs"], np.finfo(float).tiny),
                "scalar_floor_ratio_median": cm["Ephi_direct"]["median_abs"]
                / max(rm["Ephi_direct"]["median_abs"], np.finfo(float).tiny),
            }
        )

    chosen = next(row for row in convergence if row["window"] == 9 and row["degree"] == 8)
    c = chosen["current"]
    strict_non_scalar = (
        c["E00"]["max_abs"] < policy["background_residual_abs"]
        and c["E11"]["max_abs"] < policy["background_residual_abs"]
        and c["JA"]["max_abs"] < policy["background_residual_abs"]
    )
    stored_max = float(np.max(np.abs(stored.loc[stored_mask, "Ephi_residual"].to_numpy(float))))
    reference_floor = chosen["archived_exact_reference_direct_replay"]["Ephi_direct"]["max_abs"]
    current_floor = c["Ephi_direct"]["max_abs"]
    scalar_resolution_consistent = bool(current_floor <= reference_floor)

    report = {
        "status": "PASS_NONSCALAR_SCALAR_ALGEBRAIC_CERT_PENDING"
        if strict_non_scalar and scalar_resolution_consistent
        else "FAIL",
        "equation_source": "Heisenberg--Tsujikawa 2018, Eqs. (2.14)--(2.20)",
        "production_domain": "0.62 < u < 0.70",
        "policy_background_residual_abs": policy["background_residual_abs"],
        "strict": {
            "E00": c["E00"]["max_abs"] < policy["background_residual_abs"],
            "E11": c["E11"]["max_abs"] < policy["background_residual_abs"],
            "zero_electric_current_JA": c["JA"]["max_abs"] < policy["background_residual_abs"],
            "scalar_Ephi_machine_precision": False,
        },
        "selected_w9_d8": chosen,
        "archived_reference_stored_Ephi_max_abs": stored_max,
        "scalar_resolution_statement": {
            "current_direct_Jphi_prime_minus_Pphi_max_abs": current_floor,
            "archived_exact_reference_same_direct_evaluator_max_abs": reference_floor,
            "current_over_reference_floor": current_floor / reference_floor,
            "consistent_with_known_direct_differentiation_floor": scalar_resolution_consistent,
            "release_tolerance_not_relaxed": True,
            "remaining_requirement": (
                "reproduce an algebraic/identity scalar-EOM evaluator at the strict "
                "background residual tolerance before ABSOLUTE_FULL_CLOSURE"
            ),
        },
        "derivative_convergence": convergence,
        "absolute_full_closure": False,
    }
    return report


def main():
    out = ROOT / "data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_BACKGROUND_EOM_AUDIT.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    report = audit(ROOT)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
