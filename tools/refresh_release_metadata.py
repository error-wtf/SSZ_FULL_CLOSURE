"""Deterministically stamp MODEL_LOCK.json / EVIDENCE_INDEX.json with the
current HEAD (release-packaging step).  Run BEFORE committing a release:
    PYTHONPATH=src python tools/refresh_release_metadata.py
The provenance preflight test fails with instructions if these fields go
stale relative to HEAD."""
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def head() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()


def refresh() -> dict:
    commit = head()
    changed = {}
    for name in ("MODEL_LOCK.json", "EVIDENCE_INDEX.json"):
        p = ROOT / name
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        if d.get("git_commit") != commit:
            d["git_commit"] = commit
            p.write_text(json.dumps(d, indent=1) + "\n")
            changed[name] = commit
    return changed


if __name__ == "__main__":
    print(refresh() or "already current")
