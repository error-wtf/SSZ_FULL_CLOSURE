#!/usr/bin/env python3
"""Audit the physical 0.70<u<0.715 strong-field transition contract.

This tool does not claim a solved transition.  It freezes the anti-target-fitting
logic: the current hybrid is a left seed, the pure-H core is the right seed,
and the stable inner light ring must be crossed by one regular action.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central, principal_audit
from ssz_p5.production.strong_field_transition import (
    REQUIRED_L,
    U_HISTORICAL_SEAM,
    U_INNER_LIGHT_RING,
    U_LEFT,
    U_RIGHT,
    electric_constitutive_pivot,
    g2xx_kinetic_null_test,
    kinetic_transition_scan,
)
from ssz_p5.production.svt_background_eom import evaluate_svt_background, residual_metrics
from ssz_p5.production.horndeski_tubular import g4xx_background_null_section
from ssz_p5.production.strong_field_transition import compact_transition_bump

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "data/generated/strong_field_transition_2026-09-20"
REPORT = OUTDIR / "STRONG_FIELD_TRANSITION_AUDIT.json"
MODES = OUTDIR / "HYBRID_RIGHT_TRANSITION_MODE_SCAN.csv"


def _old_inner_probe(root: Path) -> dict:
    from ssz_p5.reducer.kinetic_schur import kinetic_schur
    p = root / "data/prestaged/direct41/inner_selected_candidate_41of41.csv"
    d = pd.read_csv(p).sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    # stay away from the extreme core endpoint where the old explicit reducer is known to be ill-conditioned
    mask = (u >= 0.7100) & (u <= 0.7140)
    inds = np.flatnonzero(mask)
    r = kinetic_schur(d, 6)
    eig = np.linalg.eigvalsh((r["K"] + r["K"].swapaxes(1, 2)) / 2.0)
    j = int(inds[np.argmin(eig[inds, 0])])
    left = int(np.argmin(np.abs(u - 0.7100)))
    return {
        "source": str(p.relative_to(root)),
        "L": 6,
        "eig_min_near_u_0p71": float(eig[left, 0]),
        "min_eig_0p710_to_0p714": float(eig[j, 0]),
        "u_of_min": float(u[j]),
        "final_member_eligible": False,
    }


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build, bulk_scans = principal_audit(ROOT)
    action = build.action.sort_values("x").reset_index(drop=True)
    stream = build.direct41.sort_values("x").reset_index(drop=True)

    # Background diagnostic in the right buffer; not a certified production-domain extension.
    bg = evaluate_svt_background(action, window=9, degree=8)
    um = action.u.to_numpy(float)
    right_mask = (um >= U_LEFT) & (um < U_HISTORICAL_SEAM)
    bg_metrics = residual_metrics(bg, right_mask)

    chi = electric_constitutive_pivot(action)
    lr_i = int(np.argmin(np.abs(um - U_INNER_LIGHT_RING)))
    electric_mask = right_mask & (np.abs(action.A0prime.to_numpy(float)) > 1e-9)

    modes, transition = kinetic_transition_scan(stream, REQUIRED_L)
    modes.to_csv(MODES, index=False)

    # Falsification of a tempting one-channel quartic patch.  If its background
    # compensation pivot is singular, it must not be used as a continuation
    # control by itself; at least one additional Horndeski normal direction is
    # required in a joint solve.
    try:
        q4 = compact_transition_bump(action.u.to_numpy(float))
        sec4 = g4xx_background_null_section(action, q4)
        g4_single = {
            "status": "AVAILABLE",
            "min_abs_background_pivot": float(np.min(np.abs(sec4.background_pivot.to_numpy(float)))),
            "max_abs_background_residual": float(np.max(np.abs(sec4[["delta_E00", "delta_E11", "delta_E22"]].to_numpy(float)))),
        }
    except RuntimeError as exc:
        g4_single = {
            "status": "REJECTED_RANK_DEFICIENT_AS_SINGLE_CHANNEL",
            "reason": str(exc),
            "required_next_step": "joint Horndeski normal bundle / rank-resolved background-preserving solve",
        }

    # high-L finite-set scaling is diagnostic only; final closure requires the analytic unreduced principal symbol.
    L = np.array([row["L"] for row in bulk_scans if row["L"] >= 110], float)
    k = np.array([row["min_eig_K"] for row in bulk_scans if row["L"] >= 110], float)
    cr = np.array([row["min_cr2"] for row in bulk_scans if row["L"] >= 110], float)
    pK = float(np.polyfit(np.log(L), np.log(k), 1)[0])
    pR = float(np.polyfit(np.log(L), np.log(cr), 1)[0])

    report = {
        "status": "OPEN_ACTION_CONTINUATION_REQUIRED",
        "scope": "physical right strong-field transition, not historical Inner target replay",
        "physical_domain": [U_LEFT, U_RIGHT],
        "historical_seam_u": U_HISTORICAL_SEAM,
        "historical_seam_is_physical_boundary": False,
        "inner_light_ring_u": U_INNER_LIGHT_RING,
        "left_anchor": "ELECTRIC_HYBRID_ONSHELL_CENTRAL at u=0.70",
        "right_anchor": "punctured pure-Horndeski core at u=0.715",
        "hybrid_bulk_0p62_0p70": {
            "all_required_L_K_and_radial_positive": bool(all(row["pass"] for row in bulk_scans)),
            "scans": bulk_scans,
        },
        "right_buffer_background_diagnostic": bg_metrics,
        "electric_branch_nondegeneracy": {
            "definition": "chi_E = partial J_A / partial A0prime at fixed phi,X using f2FF+4X f2FY",
            "light_ring_u_grid": float(um[lr_i]),
            "light_ring_A0prime": float(action.A0prime.iloc[lr_i]),
            "light_ring_chi_E": float(chi[lr_i]),
            "min_abs_chi_E_where_A0prime_nonzero": float(np.min(np.abs(chi[electric_mask]))),
            "finite": bool(np.all(np.isfinite(chi[electric_mask]))),
        },
        "unmodified_hybrid_right_extension": transition,
        "g2xx_kinetic_direction": g2xx_kinetic_null_test(stream, action),
        "g4xx_single_channel_background_null_test": g4_single,
        "historical_inner_probe": _old_inner_probe(ROOT),
        "high_L_diagnostic": {
            "finite_sample_log_slope_min_eig_K": pK,
            "finite_sample_log_slope_min_cr2": pR,
            "certification": "NOT_A_HIGH_L_PROOF",
            "required_next_gate": "analytic high-L reduction of the unreduced principal symbol / canonically normalized characteristics",
        },
        "solver_contract": {
            "action_first": True,
            "historical_coefficient_profiles_are_targets": False,
            "single_transition_0p70_to_0p715": True,
            "electric_support_required_through_inner_light_ring": True,
            "dehair_only_after_light_ring": True,
            "background_null_controls_preferred": True,
            "mode_tracking_required": True,
            "rank_loss_is_physical_event_not_regularization_target": True,
            "absolute_action_to_41_reemission_required": True,
            "angular_gate_must_regress_exact_SVT_oracle_before_hybrid": True,
            "qnm_boundary_conditions": "regular center + outgoing infinity",
        },
        "promotion": {
            "strong_field_transition": False,
            "global_direct41": False,
            "global_KRGSM": False,
            "same_operator_QNM": False,
            "absolute_full_closure": False,
        },
    }
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    print(f"WROTE {REPORT}")
    print(f"WROTE {MODES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
