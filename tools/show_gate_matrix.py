#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from ssz_p5.observability import build_gate_matrix

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    rows = build_gate_matrix(ROOT)
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
