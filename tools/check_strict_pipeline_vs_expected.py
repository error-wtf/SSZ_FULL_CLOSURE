#!/usr/bin/env python3
"""Compare strict-pipeline FAILED checks against the registered expected set.

The strict pipeline on this branch is EXPECTED to fail with exactly the three
documented blockers (data/diagnostic/EXPECTED_STRICT_PIPELINE_FAILURES.json):
  software.pytest                     <- pre-existing G05 member roundtrip numeric
  finite_l.central_regional_kinetic   <- core K_scalar pocket (claim gate)
  absolute_closure.complete_production_chain <- QNM disabled (claim gate)
Same set = reproduction contract holds. Extra/missing = regression or recovery
(either requires updating the registered expectation with justification).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "FULL_PIPELINE_REPORT.json"
EXPECTED = ROOT / "data/diagnostic/EXPECTED_STRICT_PIPELINE_FAILURES.json"


def main() -> int:
    rep = json.loads(REPORT.read_text())
    exp = json.loads(EXPECTED.read_text())
    failed = {f"{c['stage']}.{c['name']}": c
              for c in rep.get("checks", []) if c.get("status") == "FAIL"}
    expected = {e["check"]: e for e in exp["expected_failed_checks"]}
    missing = sorted(set(expected) - set(failed))
    extra = sorted(set(failed) - set(expected))
    ok = not missing and not extra
    print(json.dumps({
        "expected_failures": sorted(expected),
        "observed_failures": sorted(failed),
        "missing_expected_failures": missing,
        "unexpected_new_failures": extra,
        "reproduction_contract_holds": ok,
    }, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
