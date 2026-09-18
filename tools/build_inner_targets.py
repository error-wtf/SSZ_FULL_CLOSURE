#!/usr/bin/env python3
"""Build Inner lower-order targets and their explicit background-null action controls."""

import json
from pathlib import Path

from ssz_p5.production.inner_targets import build_inner_targets

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    report = build_inner_targets(ROOT, ROOT / "data/generated/inner_controls")
    keys = ("status", "central_endpoint_radial_jets", "production_changed")
    print(json.dumps({k: report[k] for k in keys}, indent=2))
