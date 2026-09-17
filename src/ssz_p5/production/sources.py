"""Canonical source registry for the locked regional production member.

Codex and production code should use these paths instead of repeatedly searching
by filename.  PRESTAGED entries are regression/handoff aids, not direct-closure
certificates unless explicitly promoted by the final generator.
"""

from __future__ import annotations

SOURCE_REGISTRY = {
    "weak_exterior_H": {
        "coeff_reference": "data/prestaged/direct41/weak_exterior_39of41.csv",
        "role": "PRESTAGED_REGRESSION",
    },
    "outer_same_action_H_SVT": {
        "background": (
            "data/production/"
            "ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv"
        ),
        "coeff_reference": "data/prestaged/direct41/outer_resolved_raw_41of41.csv",
        "role": "BACKGROUND_AUTHORITATIVE_COEFF_PRESTAGED",
    },
    "central_exact_SVT": {
        "lower_jets": (
            "data/authoritative/"
            "ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv"
        ),
        "unreduced_even": (
            "data/production/"
            "ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv"
        ),
        "exact_regression": (
            "data/production/"
            "ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv"
        ),
        "coeff_reference": "data/prestaged/direct41/central_selected_corrected_41of41.csv",
        "odd_reference": "data/prestaged/odd/central_svt_odd_profile.csv",
        "role": "AUTHORITATIVE_EXISTING_MEMBER_WITH_PRESTAGED_NORMALIZATION",
    },
    "inner_same_action_SVT_H": {
        "background": (
            "data/production/"
            "ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv"
        ),
        "coeff_reference": "data/prestaged/direct41/inner_selected_candidate_41of41.csv",
        "role": "BACKGROUND_AUTHORITATIVE_COEFF_PRESTAGED_CANDIDATE",
    },
    "punctured_H_core": {
        "coeff_reference": "data/prestaged/direct41/core_selected_41of41.csv",
        "action_g4xx": (
            "data/production/"
            "ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
        ),
        "action_g5": (
            "data/production/"
            "ssz_p5_F1b_FULL_G5_subcore_onshell_candidate_2026-09-14.csv"
        ),
        "role": "ACTION_AUTHORITATIVE_COEFF_PRESTAGED_REGRESSION",
    },
    "analytic_center": {
        "center_handover": (
            "data/production/"
            "ssz_p5_F1b_FINAL_Cinf_center_to_punctured_handover_2026-09-14.csv"
        ),
        "role": "AUTHORITATIVE_ANALYTIC_CENTER",
    },
}

# The selected central representation combines unreduced coefficients with the
# explicit lower-order completion. Partial action replays must not override it.
CENTRAL_PRODUCTION_SOURCES = {
    **SOURCE_REGISTRY["central_exact_SVT"],
    "corrected_profile": (
        "archive/full_working_snapshot/"
        "ssz_p5_integrable_full_svt_lobe_corrected_profile_2026-09-12.csv"
    ),
}
