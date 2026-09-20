#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {".git", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Print the repository file tree with size and optional SHA-256.")
    parser.add_argument("--prefix", default="", help="limit to a relative path prefix, e.g. src or data/generated")
    parser.add_argument("--hash", action="store_true", help="include SHA-256 for every file")
    args = parser.parse_args()
    base = ROOT / args.prefix if args.prefix else ROOT
    if not base.exists():
        raise SystemExit(f"missing path: {base}")
    for path in sorted(p for p in base.rglob("*") if p.is_file()):
        if any(part in EXCLUDED for part in path.parts):
            continue
        rel = path.relative_to(ROOT)
        suffix = f" sha256={sha256(path)}" if args.hash else ""
        print(f"{rel}\t{path.stat().st_size}{suffix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
