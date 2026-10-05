"""P3: FROZEN_PIPELINE_RECORD — freeze the observation pipeline.

After P1 (artifact rejection) and P2 (positive recovery) pass, the exact
analysis configuration, code state and GTI policy are frozen into a
hash-bound record.  Any later analysis change requires a NEW record;
silent drift between record and analysis is detectable.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

RECORD_VERSION = "FROZEN_PIPELINE_RECORD_V1"


def freeze_pipeline_record(
    *,
    analysis_config: dict,
    code_files: dict[str, str],   # relpath -> sha256
    git_commit: str,
    gti_policy: dict,
    negative_controls: dict,
    positive_controls: dict,
) -> tuple[dict, str]:
    """Build and hash the frozen record.  Returns (record, sha256)."""
    record = {
        "record_version": RECORD_VERSION,
        "analysis_config": analysis_config,
        "code_files": dict(sorted(code_files.items())),
        "git_commit": git_commit,
        "gti_policy": gti_policy,
        "controls": {
            "negative": negative_controls,
            "positive": positive_controls,
        },
        "frozen": True,
    }
    # hash payload WITHOUT the hash field (verify pops it symmetrically)
    payload = json.dumps(record, sort_keys=True, allow_nan=False).encode()
    digest = hashlib.sha256(payload).hexdigest()
    record["sha256"] = digest
    return record, digest


def verify_pipeline_record(record: dict) -> bool:
    if record.get("record_version") != RECORD_VERSION:
        return False
    if record.get("frozen") is not True:
        return False
    check = dict(record)
    claimed = check.pop("sha256", "")
    payload = json.dumps(check, sort_keys=True, allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest() == claimed


def code_file_hashes(root: Path, relpaths: list[str]) -> dict[str, str]:
    out = {}
    for rel in relpaths:
        p = root / rel
        if not p.exists():
            out[rel] = "ABSENT"
            continue
        out[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out
