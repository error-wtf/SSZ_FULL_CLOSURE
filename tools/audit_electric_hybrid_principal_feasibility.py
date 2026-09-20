#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from ssz_p5.production.electric_hybrid_controls import kinetic_feasibility_audit

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "data/generated/absolute_attempt_2026-09-19"
OUT = OUTDIR / "ELECTRIC_HYBRID_PRINCIPAL_RECIPE_AUDIT.json"
STREAM = OUTDIR / "ELECTRIC_HYBRID_PRINCIPAL_RECIPE_REPLAY.csv"

stream, scans = kinetic_feasibility_audit(ROOT)
stream.to_csv(STREAM, index=False)
report = {
    "status": "FEASIBILITY_ONLY_NOT_BACKGROUND_CERTIFIED",
    "scope": "central 0.62<u<0.70 finite-L kinetic/pivot response of the stored primitive-control recipe",
    "validation": {
        "constraint_pivots_nonzero": all(
            row["min_abs_Dh1"] > 25 and row["min_abs_auxiliary_determinant"] > 100
            for row in scans[:3]
        ),
        "high_L_positive_42_110_420_1000": all(row["pass_kinetic"] for row in scans[3:]),
        "all_required_L_positive": all(row["pass_kinetic"] for row in scans),
    },
    "scans": scans,
    "interpretation": (
        "The stored Horndeski primitive controls strongly reduce the old Central kinetic failure "
        "while preserving large constraint-pivot margins, but L=6,12,20 remain negative near "
        "the upper Central edge. Further principal-only shaping is not promoted as closure; the "
        "next control must come from the electric/background action sector."
    ),
}
OUT.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
raise SystemExit(0)
