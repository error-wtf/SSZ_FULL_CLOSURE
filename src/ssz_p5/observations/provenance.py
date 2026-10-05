"""Bridge provenance (contract gate 1): every input hash-frozen."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def bridge_provenance_record(root: Path, inputs: list[str]) -> dict:
    """Hash every listed repo-relative input; missing files are recorded,
    never silently skipped."""
    rec = {}
    for rel in inputs:
        p = root / rel
        rec[rel] = sha256_file(p) if p.exists() else "ABSENT"
    return rec


def git_head(root: Path) -> str:
    head = root / ".git/HEAD"
    if not head.exists():
        return "UNKNOWN"
    txt = head.read_text().strip()
    if txt.startswith("ref: "):
        ref = root / ".git" / txt[5:]
        return ref.read_text().strip()[:12] if ref.exists() else "UNKNOWN"
    return txt[:12]


def write_frozen(obj: dict, path: Path) -> str:
    """Write JSON and its .sha256 sidecar; returns the hash."""
    blob = json.dumps(obj, indent=1, sort_keys=True).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blob)
    digest = hashlib.sha256(blob).hexdigest()
    path.with_suffix(path.suffix + ".sha256").write_text(digest + "\n")
    return digest
