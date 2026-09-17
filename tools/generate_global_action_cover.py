#!/usr/bin/env python3
"""Generate the authoritative human/machine-readable production action cover."""

from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.production.regions import region_table
from ssz_p5.production.sources import SOURCE_REGISTRY

OUT = ROOT / "data" / "authoritative" / "ssz_p5_GLOBAL_ACTION_COVER_2026-09-17.csv"

rows = []
for spec in region_table():
    key = spec.region.value
    src = SOURCE_REGISTRY.get(key, {})
    rows.append(
        {
            "region": key,
            "u_min": "" if spec.u_min is None else repr(spec.u_min),
            "u_max": "" if spec.u_max is None else repr(spec.u_max),
            "role": src.get("role", "AUTHORITATIVE"),
            "background": src.get("background", ""),
            "unreduced_even": src.get("unreduced_even", ""),
            "coeff_reference": src.get("coeff_reference", ""),
            "odd_reference": src.get("odd_reference", ""),
            "center_handover": src.get("center_handover", ""),
            "action_g4xx": src.get("action_g4xx", ""),
            "action_g5": src.get("action_g5", ""),
            "lower_jets": src.get("lower_jets", ""),
            "notes": spec.notes,
        }
    )
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(OUT)
