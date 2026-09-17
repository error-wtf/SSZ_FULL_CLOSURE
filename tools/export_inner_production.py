#!/usr/bin/env python3
"""Normalize the selected inner member and report frozen endpoint agreement."""
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.production.inner_export import export_inner

if __name__ == "__main__":
    report = export_inner(ROOT, ROOT / "data/generated/inner")
    print(json.dumps({k: report[k] for k in ("INNER_DIRECT_41", "SLOT_NORMALIZATION", "finite", "rows", "interface_values", "background_max_abs")}, indent=2))
    raise SystemExit(0 if report["INNER_DIRECT_41"] == "PASS" else 2)
