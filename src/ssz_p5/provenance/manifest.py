from __future__ import annotations

import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_sha256s(root: Path, sums_file: Path | None = None) -> list[str]:
    sums_file = sums_file or root / "SHA256SUMS"
    errors = []
    for line in sums_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, rel = line.split(None, 1)
        rel = rel.strip().lstrip("*")
        p = root / rel
        if not p.exists():
            errors.append(f"missing:{rel}")
        elif sha256(p) != digest:
            errors.append(f"hash:{rel}")
    return errors
