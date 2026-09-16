"""Reject historical or explicitly forbidden inputs at the production boundary."""

import fnmatch
import json
from pathlib import Path

from ..config import repo_root


def require_production_input(path: Path, root: Path | None = None) -> Path:
    root = (root or repo_root()).resolve()
    path = path.resolve()
    rel = path.relative_to(root)
    if not path.is_file():
        raise ValueError(f"Missing production input: {rel}")
    if any(
        part in {"archive", "legacy_superseded", "superseded", "diagnostic", "qnm"}
        for part in rel.parts
    ):
        raise ValueError(f"Historical/diagnostic input is not production: {rel}")
    rules = json.loads((root / "PRODUCTION_BLACKLIST.json").read_text())
    if any(fnmatch.fnmatch(path.name, pat) for pat in rules["production_forbidden_patterns"]):
        raise ValueError(f"Blacklisted production input: {rel}")
    return path
