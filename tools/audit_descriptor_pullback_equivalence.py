#!/usr/bin/env python3
"""Audit T^dagger P_descriptor T against the established Strong-H reducer."""
from __future__ import annotations
import json
from pathlib import Path

from ssz_p5.production.hsvt_eps_y import build_region
from ssz_p5.reducer.unreduced_descriptor import compare_descriptor_pullback

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/absolute_attempt_2026-09-19/DESCRIPTOR_PULLBACK_EQUIVALENCE.json"
LS = (6, 12, 20, 42, 110, 420, 1000)


def main() -> int:
    d = build_region(ROOT, "carrier")
    rows = {}
    overall = True
    for L in LS:
        rec = compare_descriptor_pullback(d, float(L))
        rows[str(L)] = rec
        overall &= bool(rec["pass"])
    report = {
        "witness_member": "P5 Horndeski + 0.01 Y, A0prime=0 (rejected globally on-shell; retained for representation audit)",
        "region": "well-conditioned Strong-H carrier",
        "identity": "T^dagger P_generalized_psi_descriptor T == P_established_3field_reducer",
        "derivative_service": "JET9D8 window=9 degree=8",
        "reason_for_operator_pullback": (
            "A naive discretize-then-Schur operation is not a valid equivalence test for the H0 "
            "Lagrange-multiplier DAE. The differential constraint map must be applied before reduction."
        ),
        "multipoles": rows,
        "DESCRIPTOR_TO_KRGSM_PULLBACK": "PASS" if overall else "FAIL",
        "DEEP_CORE_POLICY": "KEEP_DESCRIPTOR_DAE_UNREDUCED_NO_EXPLICIT_D_h1_INVERSE",
        "QNM_STATUS": "REPRESENTATION_VALIDATED; PRODUCTION_QNM_BLOCKED_UNTIL_ONSHELL_ELECTRIC_HYBRID_IS_PROMOTED",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if overall else 2


if __name__ == "__main__":
    raise SystemExit(main())
