#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ssz_p5.observability import build_evidence_index

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description="List machine-readable evidence JSON with hashes/status.")
    parser.add_argument("--status", help="optional status substring filter")
    args = parser.parse_args()
    rows = build_evidence_index(ROOT)
    if args.status:
        needle = args.status.lower()
        rows = [row for row in rows if needle in row["status"].lower()]
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
