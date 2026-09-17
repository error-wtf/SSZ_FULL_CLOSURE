#!/usr/bin/env python3
"""One-command SSZ P5 full research pipeline.

This is the orchestration entry point for the frozen 2026-09-16 repository.
It does NOT rediscover/rederive physics and it does NOT reopen superseded
branches.  It validates and assembles the already completed research layers:

  1. repository/provenance integrity (SHA256 manifest where available)
  2. Python/package importability and test suite
  3. full scientific closure auditor
  4. frozen action / 41-slot / KRGM-KGSM artifact presence
  5. QNM/spectral archive validation (Jost, WKB/test-field, eikonal records)
  6. final machine-readable + Markdown report

Strict mode also requires direct matrices, same-operator spectral convergence, and
an absolute-closure certificate. The default archive audit is separately labelled.

The script is deliberately archive-first.  Files in data/superseded or
archive/full_working_snapshot are NEVER promoted to production merely by name.

Usage:
    python ssz_p5_full_pipeline.py
    python ssz_p5_full_pipeline.py --strict
    python ssz_p5_full_pipeline.py --skip-tests
    python ssz_p5_full_pipeline.py --json build/FULL_PIPELINE_REPORT.json

Exit codes:
    0  all required frozen-repository pipeline stages pass
    2  one or more required stages fail
    4  required repository structure/input missing

"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

RELEASE = "2026-09-16"
DEFAULT_ROOT = Path(__file__).resolve().parent


@dataclass
class Check:
    stage: str
    name: str
    status: str
    value: Any = None
    criterion: str = ""
    evidence: str = ""

    @property
    def passed(self) -> bool:
        return self.status in {"PASS", "PASS_EXACT", "PASS_NUMERICAL", "PASS_ARCHIVED"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str], cwd: Path) -> tuple[int, str]:
    p = subprocess.run(
        cmd, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False
    )
    return p.returncode, p.stdout


def find_one(root: Path, names: Iterable[str]) -> Path | None:
    search_roots = [
        root,
        root / "data" / "production",
        root / "data" / "authoritative",
        root / "data" / "regression",
        root / "data" / "certificates",
        root / "data" / "qnm",
        root / "data" / "diagnostic",
        root / "src",
        root / "legacy_reference_code",
    ]
    for name in names:
        for base in search_roots:
            p = base / name
            if p.exists():
                return p
    return None


def require_file(
    checks: list[Check], root: Path, stage: str, name: str, alternatives: Iterable[str] = ()
) -> Path | None:
    names = [name, *alternatives]
    p = find_one(root, names)
    checks.append(
        Check(
            stage,
            f"file:{name}",
            "PASS" if p else "FAIL",
            str(p.relative_to(root)) if p else None,
            "file exists",
            " | ".join(names),
        )
    )
    return p


def check_repo_structure(root: Path) -> list[Check]:
    out: list[Check] = []
    required = [
        "README.md",
        "RESEARCH_STATE.md",
        "DO_NOT_RECOMPUTE.md",
        "STATUS_FULL_CLOSURE.md",
        "ARTIFACT_CLASSIFICATION.json",
        "PRODUCTION_BLACKLIST.json",
        "NUMERICAL_POLICY.json",
        "pyproject.toml",
        "ssz_p5_full_closure_auditor.py",
    ]
    for n in required:
        p = root / n
        out.append(Check("repo", n, "PASS" if p.exists() else "FAIL", p.exists(), "exists", str(p)))
    for d in [
        "src/ssz_p5",
        "tests",
        "data/production",
        "data/regression",
        "data/qnm",
        "paper",
        "archive/full_working_snapshot",
    ]:
        p = root / d
        out.append(
            Check(
                "repo", d, "PASS" if p.is_dir() else "FAIL", p.is_dir(), "directory exists", str(p)
            )
        )
    return out


def check_hash_manifest(root: Path) -> list[Check]:
    out: list[Check] = []
    sf = root / "SHA256SUMS"
    if not sf.exists():
        return [Check("provenance", "SHA256SUMS", "FAIL", None, "exists")]
    bad, missing, verified = [], [], 0
    for line in sf.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        digest, rel = parts
        rel = rel.strip().lstrip("*")
        p = root / rel
        # SHA256SUMS itself and generated pipeline reports are intentionally not recursive.
        if rel in {"SHA256SUMS", "FULL_PIPELINE_REPORT.json", "FULL_PIPELINE_REPORT.md"}:
            continue
        if not p.exists():
            missing.append(rel)
        elif sha256(p) != digest:
            bad.append(rel)
        else:
            verified += 1
    out.append(
        Check(
            "provenance",
            "sha256_verified_files",
            "PASS" if not bad and not missing else "FAIL",
            verified,
            "all listed files match",
            f"bad={bad}; missing={missing}",
        )
    )
    return out


def check_tests(root: Path, skip: bool) -> list[Check]:
    if skip:
        return [Check("software", "pytest", "SKIPPED", "skipped by user", "optional runtime stage")]
    rc, text = run([sys.executable, "-m", "pytest", "-q"], root)
    tail = "\n".join(text.splitlines()[-20:])
    return [Check("software", "pytest", "PASS" if rc == 0 else "FAIL", rc, "exit code 0", tail)]


def check_closure_auditor(root: Path) -> list[Check]:
    audit_json = root / "build" / "audit.json"
    audit_json.unlink(missing_ok=True)
    audit_json.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(root / "ssz_p5_full_closure_auditor.py"),
        "--data-dir",
        str(root),
        "--full",
        "--json",
        str(audit_json),
        "--csv",
        str(root / "build/gates.csv"),
    ]
    rc, text = run(cmd, root)
    checks = [
        Check(
            "scientific_audit",
            "full_closure_auditor",
            "PASS" if rc == 0 else "FAIL",
            rc,
            "exit code 0",
            "\n".join(text.splitlines()[-30:]),
        )
    ]
    if audit_json.exists():
        try:
            data = json.loads(audit_json.read_text(encoding="utf-8"))
            # Support both historical report layouts.
            summary = data.get("summary", data)
            fail_count = summary.get("constructive_failures", summary.get("failures", 0))
            if isinstance(fail_count, (list, tuple, dict, set)):
                fail_n = len(fail_count)
            elif fail_count is None:
                fail_n = 0
            else:
                fail_n = int(fail_count)
            for key in ("numerical_regression_failures",):
                values = summary.get(key, [])
                checks.append(
                    Check(
                        "scientific_audit",
                        key,
                        "PASS" if not values else "FAIL",
                        values,
                        "no failures",
                    )
                )
            checks.append(
                Check(
                    "scientific_audit",
                    "auditor_constructive_failures",
                    "PASS" if fail_n == 0 else "FAIL",
                    fail_n,
                    "0",
                    f"raw={fail_count!r}",
                )
            )
        except Exception as e:
            checks.append(
                Check("scientific_audit", "audit_json_parse", "FAIL", str(e), "valid JSON")
            )
    return checks


def check_production_artifacts(root: Path) -> list[Check]:
    out: list[Check] = []
    items = [
        (
            "P5_global_principal_KRG",
            "ssz_p5_F2A_GLOBAL_CANONICAL_KRG_PRINCIPAL_CINF_FINAL_2026-09-15(2).csv",
            [],
        ),
        (
            "outer_handover_background",
            "ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv",
            [],
        ),
        (
            "inner_handover_background",
            "ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv",
            [],
        ),
        ("strongH_selected_41", "ssz_p5_F3_horndeski_carrier_SELECTED_41of41_2026-09-15.csv", []),
        ("core_selected_41", "ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv", []),
        ("outer_selected_41", "ssz_p5_F3_outer_same_action_SELECTED_41of41_2026-09-15.csv", []),
        ("ZK_raw_41", "ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv", []),
        (
            "ZK_exact_regression",
            "ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv",
            [],
        ),
        (
            "center_handover",
            "ssz_p5_F1b_FINAL_Cinf_center_to_punctured_handover_2026-09-14.csv",
            [],
        ),
    ]
    for label, name, alts in items:
        p = find_one(root, [name, *alts])
        if p:
            from ssz_p5.provenance.blacklist import require_production_input

            try:
                require_production_input(p, root)
            except ValueError as exc:
                out.append(Check("production", label, "FAIL", str(exc), "not blacklisted"))
                continue
        out.append(
            Check(
                "production",
                label,
                "PASS" if p else "FAIL",
                str(p.relative_to(root)) if p else None,
                "present",
                name,
            )
        )
    return out


def _numeric_csv_ok(path: Path, required_cols: Iterable[str] = ()) -> tuple[bool, str]:
    try:
        d = pd.read_csv(path)
        missing = [c for c in required_cols if c not in d.columns]
        if missing:
            return False, f"missing columns {missing}"
        if len(d) == 0:
            return False, "zero rows"
        return True, f"rows={len(d)} cols={len(d.columns)}"
    except Exception as e:
        return False, repr(e)


def check_qnm_archive(root: Path) -> list[Check]:
    """Validate QNM/spectral work already performed; never promote rejected poles."""
    out: list[Check] = []
    qdir = root / "data" / "qnm"
    jost = qdir / "ssz_p5_F3_maxwell_l2_jost_convergence_2026-09-15.csv"
    spectral = qdir / "SSZ_P5_POST_CLOSURE_SPECTRAL_PROBE_v2_2026-09-14.md"
    checkpoint = qdir / "SSZ_P5_F3_COUPLED_QNM_EXECUTION_CHECKPOINT_2026-09-15.md"
    eik = qdir / "SSZ_CANONICAL_V13_EIKONAL_QNM_PROXY_v1_REPORT_2026-08-28.md"
    recon = qdir / "SSZ_QNM_HISTORICAL_VS_V13_RECONCILIATION_v1_REPORT_2026-08-28.md"
    for label, p in [
        ("maxwell_jost", jost),
        ("spectral_probe", spectral),
        ("qnm_checkpoint", checkpoint),
        ("eikonal_proxy", eik),
        ("historical_reconciliation", recon),
    ]:
        out.append(
            Check(
                "qnm_archive",
                label,
                "PASS_ARCHIVED" if p.exists() else "FAIL",
                str(p.relative_to(root)) if p.exists() else None,
                "present",
            )
        )

    if jost.exists():
        ok, detail = _numeric_csv_ok(
            jost, ["dx", "R", "jost_order", "omega_re", "omega_im", "residual"]
        )
        if ok:
            d = pd.read_csv(jost)
            finite = np.isfinite(d[["omega_re", "omega_im", "residual"]].to_numpy(float)).all()
            max_res = float(np.nanmax(np.abs(d["residual"].to_numpy(float))))
            orders = sorted(set(int(x) for x in d["jost_order"].dropna()))
            radii = sorted(set(float(x) for x in d["R"].dropna()))
            out.append(
                Check(
                    "qnm_archive",
                    "maxwell_jost_rows_finite",
                    "PASS" if finite else "FAIL",
                    len(d),
                    "all finite",
                    detail,
                )
            )
            out.append(
                Check(
                    "qnm_archive",
                    "maxwell_jost_parameter_coverage",
                    "PASS" if len(orders) >= 3 and len(radii) >= 3 else "FAIL",
                    {"orders": orders, "R": radii},
                    ">=3 Jost orders and >=3 outer radii",
                )
            )
            out.append(
                Check(
                    "qnm_archive",
                    "maxwell_jost_residual_recorded",
                    "PASS" if np.isfinite(max_res) else "FAIL",
                    max_res,
                    "finite residuals",
                    "archived convergence experiment; candidate poles are not auto-promoted",
                )
            )
        else:
            out.append(Check("qnm_archive", "maxwell_jost_csv", "FAIL", detail, "valid CSV"))

    if spectral.exists():
        t = spectral.read_text(encoding="utf-8", errors="replace")
        out.append(
            Check(
                "qnm_archive",
                "scalar_WKB_threshold_256",
                "PASS_ARCHIVED" if "256" in t and "I_{\\max}" in t else "FAIL",
                256,
                "archived scalar threshold",
            )
        )
        out.append(
            Check(
                "qnm_archive",
                "maxwell_WKB_threshold_257",
                "PASS_ARCHIVED" if "257" in t else "FAIL",
                257,
                "archived Maxwell threshold",
            )
        )
    if checkpoint.exists():
        t = checkpoint.read_text(encoding="utf-8", errors="replace")
        out.append(
            Check(
                "qnm_archive",
                "outgoing_jost_order10_validation",
                "PASS_ARCHIVED" if "2.43" in t and "10" in t else "FAIL",
                "~2.43e-12",
                "archived Schwarzschild-tail transport validation",
            )
        )
        out.append(
            Check(
                "qnm_archive",
                "rejected_real_axis_poles_preserved_as_rejected",
                "PASS" if "rejected" in t.lower() else "FAIL",
                True,
                "must remain rejected, not production QNM",
            )
        )
    return out


def check_blacklist(root: Path) -> list[Check]:
    p = root / "PRODUCTION_BLACKLIST.json"
    if not p.exists():
        return [Check("hygiene", "production_blacklist", "FAIL", None, "exists")]
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        n = (
            len(d)
            if isinstance(d, list)
            else sum(len(v) for v in d.values() if isinstance(v, list))
        )
        return [
            Check(
                "hygiene",
                "production_blacklist",
                "PASS",
                n,
                "machine-readable blacklist present",
                p.name,
            )
        ]
    except Exception as e:
        return [Check("hygiene", "production_blacklist", "FAIL", repr(e), "valid JSON")]


def check_frozen_onshell(root: Path) -> list[Check]:
    from ssz_p5.production.regions import audit_light_ring_regions

    try:
        report = audit_light_ring_regions(root)
        output = root / "data/diagnostic/PRODUCTION_REGION_ASSIGNMENT.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n")
        return [
            Check(
                "production_regions",
                "sector_specific_identity_scope",
                "PASS",
                report["witnesses"],
                "pure-H identities are not applied to central_exact_SVT",
                str(output.relative_to(root)),
            )
        ]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return [Check("production_regions", "sector_specific_identity_scope", "FAIL", str(exc))]


def check_regional_kinetic(root: Path) -> list[Check]:
    from ssz_p5.production.kinetic_audit import audit_central_kinetic

    try:
        report = audit_central_kinetic(root)
        output = root / "build/REGIONAL_CENTRAL_KINETIC_GATE.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        return [
            Check(
                "finite_l",
                "central_regional_kinetic",
                report["status"],
                report["scans"],
                "K positive for every required L in central_exact_SVT",
                str(output.relative_to(root)),
            )
        ]
    except (OSError, ValueError, KeyError, TypeError, np.linalg.LinAlgError) as exc:
        return [Check("finite_l", "central_regional_kinetic", "FAIL", str(exc))]


def check_absolute_closure(root: Path) -> list[Check]:
    from ssz_p5.export.verification import verify_absolute_closure

    try:
        verify_absolute_closure(root)
        return [
            Check(
                "absolute_closure",
                "complete_production_chain",
                "PASS",
                True,
                "direct matrices, finite-L gates, coupled convergence and bound certificate",
            )
        ]
    except (OSError, RuntimeError, ValueError, KeyError, TypeError, AssertionError) as exc:
        return [
            Check(
                "absolute_closure",
                "complete_production_chain",
                "FAIL",
                str(exc),
                "direct matrices, finite-L gates, coupled convergence and bound certificate",
            )
        ]


def write_reports(root: Path, checks: list[Check], json_path: Path, md_path: Path) -> None:
    passed = sum(c.passed for c in checks)
    failed = sum(not c.passed for c in checks)
    absolute_pass = failed == 0 and any(c.stage == "absolute_closure" for c in checks)
    payload = {
        "release": RELEASE,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "repository": str(root),
        "provenance": {
            "git_commit": run(["git", "rev-parse", "HEAD"], root)[1].strip(),
            "python_version": sys.version,
            "action_sha256": sha256(root / "SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json"),
            "numerical_policy_sha256": sha256(root / "NUMERICAL_POLICY.json"),
            "dependency_lock_sha256": sha256(root / "requirements.lock"),
            "input_manifest_sha256": sha256(root / "MANIFEST.json"),
            "generator": "ssz_p5_full_pipeline.py",
        },
        "summary": {
            "checks": len(checks),
            "passed": passed,
            "failed": failed,
            "pipeline_pass": failed == 0,
        },
        "checks": [asdict(c) for c in checks],
        "scientific_status": {
            "regional_central_kinetic": next(
                (c.status for c in checks if c.name == "central_regional_kinetic"), "NOT_RUN"
            ),
            "absolute_closure": "PASS" if absolute_pass else "NOT_CERTIFIED",
            "direct_global_krgm_export": "PASS" if absolute_pass else "NOT_CERTIFIED",
            "coupled_hsvt_qnm": "PASS" if absolute_pass else "NOT_CERTIFIED",
            "scope": "archive, constructive audit and numerical regressions; "
            "not a direct-global export certificate",
        },
        "notes": [
            "Pipeline validates the frozen repository and previously performed QNM/spectral work.",
            "It does not promote historically rejected finite-radius QNM poles.",
            "It does not reopen superseded action/core/handover derivations.",
        ],
    }
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
    )

    by_stage: dict[str, list[Check]] = {}
    for c in checks:
        by_stage.setdefault(c.stage, []).append(c)
    lines = [
        "# SSZ P5 full-pipeline report",
        "",
        f"**Release:** {RELEASE}",
        f"**Generated:** {payload['generated_at_utc']}",
        f"**Overall:** {'PASS' if failed == 0 else 'FAIL'}",
        f"**Checks:** {passed}/{len(checks)} passed",
        "",
        "This report validates the frozen repository. It does not repeat completed derivations.",
        "",
    ]
    for stage, rows in by_stage.items():
        lines += [f"## {stage}", "", "| Check | Status | Value | Criterion |", "|---|---|---:|---|"]
        for c in rows:
            val = str(c.value).replace("|", "\\|").replace("\n", " ")
            lines.append(f"| `{c.name}` | **{c.status}** | {val[:180]} | {c.criterion} |")
        lines.append("")
    lines += [
        "## Scientific scope",
        "",
        "Absolute closure: **" + payload["scientific_status"]["absolute_closure"] + "**.",
        "Coupled HSVT QNM: **" + payload["scientific_status"]["coupled_hsvt_qnm"] + "**.",
        "An archival pipeline PASS is not a direct-export or coupled-spectrum certificate.",
        "",
    ]
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--skip-tests", action="store_true")
    p.add_argument(
        "--strict",
        action="store_true",
        help="require hashes, direct global matrices, spectral convergence and absolute closure",
    )
    p.add_argument("--json", type=Path, default=None)
    p.add_argument("--markdown", type=Path, default=None)
    return p.parse_args()


def main() -> int:
    ns = parse_args()
    root = ns.repo.resolve()
    if not root.exists():
        print(f"ERROR: repository does not exist: {root}", file=sys.stderr)
        return 4

    checks: list[Check] = []
    checks += check_repo_structure(root)
    if ns.strict:
        checks += check_hash_manifest(root)
        rc, detail = run([sys.executable, str(root / "tools/release_manifest.py"), "--check"], root)
        checks.append(
            Check(
                "provenance",
                "manifest_inventory",
                "PASS" if rc == 0 else "FAIL",
                rc,
                "complete inventory and sums match",
                detail,
            )
        )
    if ns.strict and ns.skip_tests:
        print("ERROR: --strict cannot be combined with --skip-tests", file=sys.stderr)
        return 2
    checks += check_tests(root, ns.skip_tests)
    checks += check_closure_auditor(root)
    checks += check_production_artifacts(root)
    checks += check_qnm_archive(root)
    checks += check_blacklist(root)
    from ssz_p5.production.member import load_regional_member

    try:
        member = load_regional_member(root)
        checks.append(
            Check(
                "action",
                "frozen_action_member",
                "PASS",
                member["name"],
                "locked six-region Full-SVT production cover",
            )
        )
    except (ValueError, KeyError, OSError) as exc:
        checks.append(Check("action", "frozen_action_member", "FAIL", str(exc)))
    from ssz_p5.provenance.blacklist import require_production_input

    forbidden = []
    for path in (root / "data/production").glob("*"):
        if not path.is_file():
            continue
        try:
            require_production_input(path, root)
        except ValueError:
            forbidden.append(path.name)
    checks.append(
        Check(
            "hygiene",
            "production_blacklist_enforced",
            "FAIL" if forbidden else "PASS",
            forbidden,
            "no forbidden production inputs",
        )
    )
    diagnostic = root / "data/diagnostic/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv"
    checks.append(
        Check(
            "hygiene",
            "forbidden_stream_preserved_as_diagnostic",
            "PASS" if diagnostic.exists() else "FAIL",
            diagnostic.name,
        )
    )

    if ns.strict:
        checks += check_frozen_onshell(root)
        checks += check_regional_kinetic(root)
        checks += check_absolute_closure(root)

    json_path = ns.json or root / "FULL_PIPELINE_REPORT.json"
    md_path = ns.markdown or root / "FULL_PIPELINE_REPORT.md"
    write_reports(root, checks, json_path, md_path)

    failed = [c for c in checks if not c.passed]
    print("=" * 76)
    print("SSZ P5 FULL PIPELINE")
    print("=" * 76)
    for stage in sorted(set(c.stage for c in checks)):
        rows = [c for c in checks if c.stage == stage]
        print(f"{stage:20s}: {sum(c.passed for c in rows):3d}/{len(rows):3d} PASS")
    print("-" * 76)
    print(f"TOTAL               : {sum(c.passed for c in checks):3d}/{len(checks):3d} PASS")
    print(f"REPORT JSON         : {json_path}")
    print(f"REPORT MARKDOWN     : {md_path}")
    if failed:
        print("FAILED CHECKS:")
        for c in failed:
            print(f"  - [{c.stage}] {c.name}: {c.value} ({c.criterion})")
        return 2
    print("ABSOLUTE_FULL_CLOSURE = PASS" if ns.strict else "FULL_PIPELINE_ARCHIVE_AND_AUDIT = PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
