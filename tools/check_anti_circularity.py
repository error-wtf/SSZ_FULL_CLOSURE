#!/usr/bin/env python3
"""ANTI_CIRCULARITY_GATE — machine-checkable enforcement of the
Anti-Circular Forward-Prediction Law (docs/ANTI_CIRCULAR_FORWARD_PREDICTION_LAW.md).

Usage:
    python tools/check_anti_circularity.py <comparison_artifact.json>

Exit codes:
    0  GATE PASS (all mandatory fields present, freeze precedes unblinding)
    1  GATE FAIL (missing fields or ordering violation)
    2  usage error

The gate is intentionally dependency-free (stdlib only) so any CI or
consumer (e.g. Spectroscopy-Bridge) can run it.
"""
import json
import sys
from datetime import datetime

MANDATORY = [
    "theory_freeze_commit",
    "theory_freeze_timestamp",
    "observational_unblind_timestamp",
    "data_access_boundary",
    "preregistration_hash",
    "prediction_hash",
]


def parse_ts(s):
    s = s.replace("Z", "+00:00")
    d = datetime.fromisoformat(s)
    if d.tzinfo is None:
        raise ValueError("timestamp must carry an explicit timezone")
    return d


def main():
    if len(sys.argv) != 2:
        print("usage: check_anti_circularity.py <comparison_artifact.json>")
        return 2
    d = json.load(open(sys.argv[1]))

    missing = [k for k in MANDATORY if k not in d or not d[k]]
    if missing:
        print("ANTI_CIRCULARITY_GATE = FAIL")
        print("missing mandatory fields:", ", ".join(missing))
        return 1

    try:
        tf = parse_ts(str(d["theory_freeze_timestamp"]))
        ou = parse_ts(str(d["observational_unblind_timestamp"]))
    except ValueError as e:
        print("ANTI_CIRCULARITY_GATE = FAIL")
        print("timestamp parse error:", e)
        return 1

    if not (tf < ou):
        print("ANTI_CIRCULARITY_GATE = FAIL")
        print(f"ordering violation: theory_freeze_timestamp ({tf.isoformat()}) "
              f"must be strictly BEFORE observational_unblind_timestamp ({ou.isoformat()})")
        return 1

    print("ANTI_CIRCULARITY_GATE = PASS")
    print(f"  theory_freeze_commit           : {d['theory_freeze_commit']}")
    print(f"  theory_freeze_timestamp        : {tf.isoformat()}")
    print(f"  observational_unblind_timestamp: {ou.isoformat()}")
    print(f"  preregistration_hash           : {str(d['preregistration_hash'])[:16]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
