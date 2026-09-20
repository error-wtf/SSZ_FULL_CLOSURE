#!/usr/bin/env python3
"""Build or verify the immutable repository inventory (runtime outputs excluded)."""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    "build",
    "dist",
}
EXCLUDED_FILES = {
    "MANIFEST.json",
    "SHA256SUMS",
    "FULL_PIPELINE_REPORT.json",
    "FULL_PIPELINE_REPORT.md",
}


def inventory():
    result = []
    for p in sorted(ROOT.rglob("*")):
        rel = p.relative_to(ROOT)
        if not p.is_file() or any(x in EXCLUDED_DIRS or x.endswith(".egg-info") for x in rel.parts):
            continue
        if len(rel.parts) == 1 and p.name in EXCLUDED_FILES:
            continue
        if p.suffix in {".pyc", ".pyo"}:
            continue
        result.append(
            {
                "path": str(rel),
                "size_bytes": p.stat().st_size,
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
        )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    files = inventory()
    if args.check:
        original = json.loads((ROOT / "MANIFEST.json").read_text())["files"]
        if original != files:
            raise SystemExit("Manifest differs from current files")
        expected = "".join(x["sha256"] + "  " + x["path"] + "\n" for x in files)
        expected += (
            hashlib.sha256((ROOT / "MANIFEST.json").read_bytes()).hexdigest() + "  MANIFEST.json\n"
        )
        if (ROOT / "SHA256SUMS").read_text() != expected:
            raise SystemExit("SHA256SUMS differs from manifest")
        print(f"Manifest and SHA256SUMS: {len(files)} files verified")
        return
    (ROOT / "MANIFEST.json").write_text(
        json.dumps(
            {"release": "checkpoint-2026-09-19", "repository_edition": "electric-hybrid-search", "files": files}, indent=2
        )
        + "\n"
    )
    (ROOT / "SHA256SUMS").write_text(
        "".join(x["sha256"] + "  " + x["path"] + "\n" for x in files)
        + hashlib.sha256((ROOT / "MANIFEST.json").read_bytes()).hexdigest()
        + "  MANIFEST.json\n"
    )
    print(f"Manifest generated: {len(files)} files")


if __name__ == "__main__":
    main()
