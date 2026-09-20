#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "REPO_EVIDENCE_REGISTRY.json"
OUTDIR = ROOT / "data/generated/repo_observability"
JSON_OUT = OUTDIR / "ALL_EVIDENCE_RUN.json"
MD_OUT = OUTDIR / "ALL_EVIDENCE_RUN.md"


def _run(entry: dict, timeout: int) -> dict:
    raw = list(entry["command"])
    cmd = [sys.executable, *raw]
    proc = subprocess.run(
        cmd,
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    expected = int(entry.get("expected_exit", 0))
    return {
        "id": entry["id"],
        "category": entry.get("category", "full_check"),
        "scope": entry.get("scope", ""),
        "command": cmd,
        "expected_exit": expected,
        "actual_exit": proc.returncode,
        "expectation_met": proc.returncode == expected,
        "evidence": entry.get("evidence", ""),
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def _render(report: dict) -> str:
    lines = [
        "# SSZ P5 all-evidence run",
        "",
        f"Overall registry result: **{report['status']}**",
        "",
        "| Audit | Category | Expected | Actual | Registry result | Evidence |",
        "| --- | --- | ---: | ---: | --- | --- |",
    ]
    for row in report["runs"]:
        lines.append(
            f"| `{row['id']}` | {row['category']} | {row['expected_exit']} | "
            f"{row['actual_exit']} | {'PASS' if row['expectation_met'] else 'FAIL'} | "
            f"`{row['evidence']}` |"
        )
    lines += [
        "",
        "> Registry PASS means each audit returned its declared expected exit code. "
        "It does not mean Absolute Full Closure; intentional rejection witnesses may expect exit 2.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run registered evidence and classify expected rejection witnesses correctly.")
    parser.add_argument("--full", action="store_true", help="also run pytest and release manifest checks")
    parser.add_argument("--timeout", type=int, default=180, help="per-command timeout in seconds")
    parser.add_argument("--only", action="append", default=[], help="run only named registry id (repeatable)")
    args = parser.parse_args()
    registry = json.loads(REGISTRY.read_text())
    entries = list(registry["audits"])
    if args.full:
        entries.extend(registry.get("full_checks", []))
    if args.only:
        wanted = set(args.only)
        entries = [entry for entry in entries if entry["id"] in wanted]
    runs = []
    for entry in entries:
        print(f"==> {entry['id']}", flush=True)
        try:
            row = _run(entry, args.timeout)
        except subprocess.TimeoutExpired as exc:
            row = {
                "id": entry["id"],
                "category": entry.get("category", "full_check"),
                "scope": entry.get("scope", ""),
                "command": [sys.executable, *entry["command"]],
                "expected_exit": int(entry.get("expected_exit", 0)),
                "actual_exit": "TIMEOUT",
                "expectation_met": False,
                "evidence": entry.get("evidence", ""),
                "stdout_tail": (exc.stdout or "")[-4000:] if isinstance(exc.stdout, str) else "",
                "stderr_tail": (exc.stderr or "")[-4000:] if isinstance(exc.stderr, str) else "",
            }
        runs.append(row)
        print(f"    expected={row['expected_exit']} actual={row['actual_exit']} met={row['expectation_met']}")
    complete_registry = not args.only and not args.full
    if args.full and not args.only:
        complete_registry = True
    met = all(r["expectation_met"] for r in runs)
    report = {
        "schema": 1,
        "absolute_full_closure": "NOT_INFERRED_FROM_THIS_RUN",
        "selection": {
            "full_checks_included": bool(args.full),
            "only": list(args.only),
            "complete_registered_audits": bool(not args.only),
        },
        "status": ("PASS" if met else "FAIL") if not args.only else ("PASS_SELECTED" if met else "FAIL_SELECTED"),
        "runs": runs,
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(report, indent=2) + "\n")
    MD_OUT.write_text(_render(report) + "\n")
    print(JSON_OUT.relative_to(ROOT))
    print(MD_OUT.relative_to(ROOT))
    return 0 if str(report["status"]).startswith("PASS") else 2


if __name__ == "__main__":
    raise SystemExit(main())
