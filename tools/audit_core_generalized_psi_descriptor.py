#!/usr/bin/env python3
"""Audit the constraint-preserving generalized-psi descriptor on the HSVT candidate."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

from ssz_p5.production.hsvt_eps_y import build_region
from ssz_p5.reducer.unreduced_descriptor import (
    generalized_psi_descriptor,
    compare_h0_constraint,
    structural_audit,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/absolute_attempt_2026-09-18/CORE_GENERALIZED_PSI_DESCRIPTOR_AUDIT.json"
LS = (6, 12, 20, 42, 110, 420, 1000)


def main() -> int:
    result = {
        "candidate": "P5 Horndeski + 0.01 Y, A0prime=0",
        "derivative_service": "JET9D8 window=9 degree=8",
        "purpose": "Avoid explicit deep-core D_h1 inverse by retaining H0 and h1 as a descriptor constraint pair.",
        "regions": {},
    }
    overall = True
    for region in ("carrier", "core"):
        d = build_region(ROOT, region)
        rr = {"rows": len(d), "multipoles": {}}
        for L in LS:
            P = generalized_psi_descriptor(d, float(L))
            cmp = compare_h0_constraint(d, float(L), P)
            max_scaled = max(
                v["max_scaled"] for v in cmp.values()
                if isinstance(v, dict) and "max_scaled" in v
            )
            row = {
                "status": "PASS" if cmp["pass"] else "FAIL",
                "max_closed_formula_scaled_difference": max_scaled,
                "max_forbidden_H0_term": cmp["structure"]["max_forbidden_abs"],
                "H0_square_max_abs": cmp["H0_square_max_abs"],
                "operator_structure": structural_audit(P),
            }
            rr["multipoles"][str(L)] = row
            overall &= cmp["pass"]
        result["regions"][region] = rr
    result["CORE_GENERALIZED_PSI_DESCRIPTOR"] = "PASS" if overall else "FAIL"
    result["DIRECT_CORE_KRGM"] = "OPEN_DESCRIPTOR_DAE_SPECTRAL_LINEARIZATION"
    result["COUPLED_QNM"] = "BLOCKED_BY_DESCRIPTOR_TO_SPECTRAL_OPERATOR"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if overall else 2


if __name__ == "__main__":
    raise SystemExit(main())
