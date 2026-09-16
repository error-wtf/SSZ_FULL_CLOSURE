#!/usr/bin/env python3
"""Reproduce the frozen-member contradiction; exit 2 when the identity fails."""

import argparse
import json
from pathlib import Path

from ssz_p5.action.onshell import audit_frozen_carrier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=Path("build/FROZEN_ONSHELL_IDENTITY.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    report = audit_frozen_carrier(root)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(report["status"])
    print(f"r/rs={report['radius']:.16g}")
    for route in report["routes"]:
        print(f"{route['method']}, stride={route.get('stride', 'n/a')}: "
              f"residual={route['residual']:.16g}")
    print(f"Report: {args.json}")
    return 2 if report["status"] == "FAIL_ON_SHELL_IDENTITY" else 0


if __name__ == "__main__":
    raise SystemExit(main())
