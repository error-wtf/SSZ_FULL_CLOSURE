#!/usr/bin/env python3
"""Export and freeze the central representation without a background side-audit."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.production.central_export import export_central

if __name__ == "__main__":
    result = export_central(ROOT, ROOT / "data/generated/central")
    print(json.dumps({k: result[k] for k in (
        "CENTRAL_DIRECT_41", "certificate_scope", "exact_ZK_comparable_slots",
        "exact_ZK_slot_errors", "rows", "independent_action_to_41_verified")}, indent=2))
    raise SystemExit(0 if result["CENTRAL_DIRECT_41"] == "PASS" else 2)
