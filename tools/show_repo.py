#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ssz_p5.observability import build_repo_snapshot, render_snapshot_text

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description="Show the complete current SSZ P5 repository truth view.")
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    parser.add_argument("--write", type=Path, help="also write the rendered view to this path")
    parser.add_argument("--evidence-limit", type=int, default=30)
    args = parser.parse_args()
    snapshot = build_repo_snapshot(ROOT)
    rendered = json.dumps(snapshot, indent=2) + "\n" if args.json else render_snapshot_text(snapshot, evidence_limit=args.evidence_limit)
    print(rendered, end="")
    if args.write:
        path = args.write if args.write.is_absolute() else ROOT / args.write
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
