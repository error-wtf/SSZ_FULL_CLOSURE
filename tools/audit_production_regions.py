#!/usr/bin/env python3
"""Verify region assignment before applying sector-specific identities."""

import json
from pathlib import Path

from ssz_p5.production.regions import audit_light_ring_regions


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    report = audit_light_ring_regions(root)
    output = root / "data/diagnostic/PRODUCTION_REGION_ASSIGNMENT.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
