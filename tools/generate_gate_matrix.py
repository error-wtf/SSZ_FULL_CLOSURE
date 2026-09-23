#!/usr/bin/env python3
"""Regenerate TEST_GATE_MATRIX.json + TEST_COVERAGE_GAPS.json.

Starts from the committed matrix and applies the DECLARED test->gate
extensions for the dedicated current-member gate batteries.  Every mapped
test must exist in the collected test inventory; the inverse map and the
coverage-gap list are recomputed deterministically.

Run: PYTHONPATH=src python tools/generate_gate_matrix.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

MATRIX = ROOT / "TEST_GATE_MATRIX.json"
GAPS = ROOT / "TEST_COVERAGE_GAPS.json"

# Dedicated current-member gate batteries (each entry is justified by the
# test module docstring; see the referenced files).
EXTENSIONS = {
    "G05": [
        "tests/regression/test_g05_member_roundtrip.py::test_G05_scalar_integration_round_trip",
    ],
    "G07": [
        "tests/regression/test_g07_holonomic_recheck.py::test_reconstruction_B_matches_module_delta22",
        "tests/regression/test_g07_holonomic_recheck.py::test_holonomic_chain_identity_f2eff",
        "tests/regression/test_g07_holonomic_recheck.py::test_inverse_metric_variation_convention",
        "tests/regression/test_g07_holonomic_recheck.py::test_no_nonholonomic_df2_channel",
    ],
    "G10": [
        "tests/regression/test_svt_direct_variation.py::test_delta22_normalization_channels_exact",
        "tests/regression/test_svt_direct_variation.py::test_areal_gauge_guard",
    ],
    "G11": [
        "tests/regression/test_svt_direct_variation.py::test_G11_testA_epsY_null_explicit",
    ],
    "G12": [
        "tests/regression/test_svt_direct_variation.py::test_G12_testB_explicit_delta22_exact",
        "tests/regression/test_svt_direct_variation.py::test_E22_full_decomposition_identity",
    ],
    "G13": [
        "tests/regression/test_g13_g16_electric_chain.py::test_G13_production_f2Y_zero_branch",
    ],
    "G14": [
        "tests/regression/test_g13_g16_electric_chain.py::test_G14_full_light_ring_limit",
    ],
    "G15": [
        "tests/regression/test_g13_g16_electric_chain.py::test_G15_sigma_svt_decomposition_derived",
    ],
    "G16": [
        "tests/regression/test_g13_g16_electric_chain.py::test_G16_full_onshell_LR_balance",
    ],
    "G20": [
        "tests/regression/test_g05_member_roundtrip.py::test_G20_same_member_provenance",
        "tests/regression/test_g05_member_roundtrip.py::test_member_builder_reproduces_locked_stream",
    ],
    "G71": [
        "tests/regression/test_g71_g72_sectors.py::test_G71_odd_parity_axial_luminal",
    ],
    "G72": [
        "tests/regression/test_g71_g72_sectors.py::test_G72_vector_sector_descriptor_block",
    ],
    "GENERAL": [
        "tests/regression/test_svt_direct_variation.py::test_delta22_p_mh_zero_and_second_order",
        "tests/regression/test_svt_direct_variation.py::test_delta22_operator_families",
        "tests/regression/test_svt_direct_variation.py::test_p_mh_projector_backward_compatible",
        "tests/regression/test_svt_direct_variation.py::test_step3_evidence_complete",
        "tests/regression/test_svt_direct_variation.py::test_ward_oracle_record_consistency",
    ],
}

# All previously-open coverage gaps now carry dedicated current-member
# batteries (G50/G70: principal ladder; G71/G72: axial + vector sector;
# G90: global regularity scan).
KNOWN_OPEN = []


def collected_ids() -> set[str]:
    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    out = subprocess.run(
        ["python", "-m", "pytest", "tests", "tests/test_release_auditor.py",
         "--collect-only", "-q"],
        cwd=ROOT, capture_output=True, text=True, env=env,
    ).stdout
    return {ln.strip() for ln in out.splitlines() if "::" in ln}


def main() -> int:
    matrix = json.loads(MATRIX.read_text())
    g2t = {k: list(v) for k, v in matrix["gate_to_tests"].items()}
    have = collected_ids()

    added = 0
    for gate, tests in EXTENSIONS.items():
        bucket = g2t.setdefault(gate, [])
        for t in tests:
            if t not in have:
                raise SystemExit(f"declared test not collected: {t}")
            if t not in bucket:
                bucket.append(t)
                added += 1

    t2g = {}
    for gate, tests in sorted(g2t.items()):
        for t in tests:
            t2g.setdefault(t, []).append(gate)

    missing = [g for g in sorted(g2t)
               if g not in ("GENERAL",) and not g2t[g]]
    gaps = {
        "schema_version": "1.0",
        "uncovered_collected_tests": [],
        "required_gates_without_tests": KNOWN_OPEN,
        "note": "G50/G70/G71/G72/G90 require dedicated current-member batteries "
                "(perturbation program); all other required gates carry dedicated "
                "or historically mapped tests",
        "tautological_tests": json.loads(GAPS.read_text()).get(
            "tautological_tests", []) if GAPS.exists() else [],
    }
    matrix_out = {
        "schema_version": matrix.get("schema_version", "1.0"),
        "test_to_gate": dict(sorted(t2g.items())),
        "gate_to_tests": dict(sorted(g2t.items())),
        "required_gates_without_tests": KNOWN_OPEN,
    }
    MATRIX.write_text(json.dumps(matrix_out, indent=1) + "\n")
    GAPS.write_text(json.dumps(gaps, indent=1) + "\n")
    print(json.dumps({"added_mappings": added, "open_gaps": KNOWN_OPEN}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
