#!/usr/bin/env python3
"""Reproduce the central finite-L necessary gate and independent kinetic check."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ssz_p5.production.kinetic_audit import audit_central_kinetic

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json", type=Path, default=ROOT / "build/REGIONAL_CENTRAL_KINETIC_GATE.json"
    )
    args = parser.parse_args()
    report = audit_central_kinetic(ROOT)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("Central Full-SVT kinetic gate:", report["status"])
    for row in report["scans"]:
        print(f"L={row['L']:4d} min(K)={row['min_eig_K']:.9g} negative rows={row['negative_rows']}")
    print(args.json)
    raise SystemExit(0 if report["status"] == "PASS" else 2)
