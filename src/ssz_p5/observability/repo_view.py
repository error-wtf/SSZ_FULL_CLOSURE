from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

STATUS_KEYS = (
    "status",
    "CORE_GENERALIZED_PSI_DESCRIPTOR",
    "DESCRIPTOR_TO_KRGSM_PULLBACK",
    "STRONG_H_RADIAL_PENCIL_EQUIVALENCE",
    "DEEP_CORE_DESCRIPTOR_PENCIL",
    "INNER_DIRECT_41",
    "DIRECT_ACTION_REPLAY_COMPLETE",
    "inner_direct_41",
)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return {}
    return value if isinstance(value, dict) else {"value": value}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _first_status(data: dict[str, Any]) -> str:
    for key in STATUS_KEYS:
        if key in data:
            value = data[key]
            if isinstance(value, bool):
                return "PASS" if value else "FAIL"
            return str(value)
    summary = data.get("summary")
    if isinstance(summary, dict):
        for key in ("absolute_full_closure", "FULL_CONSTRUCTIVE_CLOSURE", "DIRECT_GLOBAL_KRGM_EXPORT"):
            if key in summary:
                return str(summary[key])
    return "UNCLASSIFIED"


def _scope(data: dict[str, Any]) -> str:
    for key in ("scope", "region", "identity", "interpretation", "logic"):
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip().replace("\n", " ")[:280]
    return ""


def _working_status(root: Path) -> tuple[dict[str, Any], str]:
    """Return the newest explicit working-status snapshot without rewriting history."""
    for name in (
        "FULL_CLOSURE_WORKING_STATUS_2026-09-21_PROJECTED_CONTROLLABILITY.json",
        "FULL_CLOSURE_WORKING_STATUS_2026-09-20_STRONG_FIELD.json",
        "FULL_CLOSURE_WORKING_STATUS_2026-09-21_PROJECTED_CONTROLLABILITY.json",
    ):
        p = root / name
        if p.exists():
            return _read_json(p), name
    return {}, ""


def build_evidence_index(root: Path) -> list[dict[str, Any]]:
    roots = [root / "data" / "generated", root / "build"]
    paths: list[Path] = []
    for base in roots:
        if base.exists():
            paths.extend(base.rglob("*.json"))
    for name in (
        "ABSOLUTE_CLOSURE_LEDGER.json",
        "CLOSURE_AUDIT_REPORT.json",
        "STRICT_DIRECT_AUDIT_REPORT.json",
        "PACKAGING_STATUS.json",
        "UPDATED_WORKING_STATUS_2026-09-18.json",
        "FULL_CLOSURE_WORKING_STATUS_2026-09-20_STRONG_FIELD.json",
    ):
        p = root / name
        if p.exists():
            paths.append(p)
    rows = []
    seen = set()
    for path in sorted(paths):
        rel = str(path.relative_to(root))
        if rel in seen:
            continue
        seen.add(rel)
        data = _read_json(path)
        rows.append(
            {
                "path": rel,
                "status": _first_status(data),
                "scope": _scope(data),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    return rows


def build_gate_matrix(root: Path) -> list[dict[str, Any]]:
    ledger = _read_json(root / "ABSOLUTE_CLOSURE_LEDGER.json")
    working, working_name = _working_status(root)
    rows: list[dict[str, Any]] = []

    def add(name: str, status: Any, evidence: str = "", note: str = "") -> None:
        rows.append({"gate": name, "status": str(status), "evidence": evidence, "note": note})

    add("ABSOLUTE_FULL_CLOSURE", ledger.get("absolute_full_closure", "UNKNOWN"), "ABSOLUTE_CLOSURE_LEDGER.json")
    add("SOFTWARE_TESTS", ledger.get("software_tests", "UNKNOWN"), "CURRENT_RELEASE_STATUS_2026-09-19.md")
    infra = ledger.get("verified_reusable_infrastructure", {})
    if isinstance(infra, dict):
        for key, value in infra.items():
            add(key.upper(), value, "ABSOLUTE_CLOSURE_LEDGER.json", "reusable infrastructure")
    active = ledger.get("active_search", {})
    if isinstance(active, dict):
        add("ACTIVE_SEARCH", active.get("status", "UNKNOWN"), str(active.get("member", "")))
        pf = active.get("horndeski_principal_feasibility", {})
        if isinstance(pf, dict):
            add("HORNDESKI_PRINCIPAL_FEASIBILITY", pf.get("status", "UNKNOWN"), str(pf.get("audit", "")))
    add("QNM", ledger.get("qnm", "UNKNOWN"), "ABSOLUTE_CLOSURE_LEDGER.json")
    if working:
        central = working.get("central", {})
        transition = working.get("right_strong_field_transition", {})
        downstream = working.get("downstream", {})
        tests = working.get("tests", {}).get("grouped_complete_matrix", {})
        add(
            "SOFTWARE_TESTS_CURRENT",
            f"{tests.get('passed')}/{tests.get('total')}_PASS_GROUPED_COMPLETE_MATRIX" if tests.get("status") == "PASS" else tests.get("status", "UNKNOWN"),
            working_name,
            "current repaired checkpoint; historical release test counts remain preserved separately",
        )
        add(
            "CENTRAL_DIRECT41_NORMALIZATION",
            central.get("historical_direct41_normalization", "UNKNOWN"),
            working_name,
            "scoped reproduction only; not physical closure",
        )
        add(
            "ELECTRIC_HYBRID_BULK",
            central.get("electric_hybrid_bulk_0p62_0p70", "UNKNOWN"),
            "data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_AUDIT.json",
        )
        add(
            "RIGHT_STRONG_FIELD_TRANSITION",
            transition.get("status", "UNKNOWN"),
            "data/generated/strong_field_transition_2026-09-20/STRONG_FIELD_TRANSITION_AUDIT.json",
            "one action solve on 0.70<=u<=0.715; u=0.71 is not a physical seam",
        )
        add("GLOBAL_ABSOLUTE_DIRECT41", downstream.get("global_absolute_direct41", "UNKNOWN"), working_name)
        add("GLOBAL_KRGSM", downstream.get("global_KRGSM", "UNKNOWN"), working_name)
        add("SAME_OPERATOR_QNM", downstream.get("same_operator_QNM", "UNKNOWN"), working_name)
        projected = working.get("local_projected_controllability", {})
        if projected:
            add("PROJECTED_ONSHELL_CONTROLLABILITY", projected.get("status", "UNKNOWN"), "data/generated/strong_field_transition_2026-09-20/TRANSITION_PROJECTED_CONTROLLABILITY.json", "local tangent existence gate; not a finite transition member")

    evidence = {row["path"]: row for row in build_evidence_index(root)}
    preferred = [
        "data/generated/absolute_attempt_2026-09-19/ZERO_VECTOR_LIGHT_RING_AUDIT.json",
        "data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_PRINCIPAL_RECIPE_AUDIT.json",
        "data/generated/absolute_attempt_2026-09-19/DESCRIPTOR_PULLBACK_EQUIVALENCE.json",
        "data/generated/absolute_attempt_2026-09-19/RADIAL_DESCRIPTOR_PENCIL_AUDIT.json",
        "data/generated/absolute_attempt_2026-09-18/CORE_GENERALIZED_PSI_DESCRIPTOR_AUDIT.json",
        "data/generated/inner_y_hessian/INNER_F2_4D_ALGEBRAIC_REACHABILITY.json",
        "data/generated/central/CENTRAL_DIRECT_41_CERTIFICATE.json",
        "data/generated/inner/INNER_DIRECT_41_CERTIFICATE.json",
    ]
    existing = {r["evidence"] for r in rows}
    for path in preferred:
        if path in evidence and path not in existing:
            row = evidence[path]
            add(Path(path).stem, row["status"], path, row["scope"])
    return rows


def build_member_matrix(root: Path) -> list[dict[str, str]]:
    ledger = _read_json(root / "ABSOLUTE_CLOSURE_LEDGER.json")
    working, _working_name = _working_status(root)
    rows: list[dict[str, str]] = []
    if working:
        central = working.get("central", {})
        transition = working.get("right_strong_field_transition", {})
        rows.extend([
            {
                "member": "historical_central_direct41",
                "status": str(central.get("historical_direct41_normalization", "UNKNOWN")),
                "reason": "normalization/reproduction scope only; not final physical member",
            },
            {
                "member": "electric_hybrid_bulk_0p62_0p70",
                "status": str(central.get("electric_hybrid_bulk_0p62_0p70", "UNKNOWN")),
                "reason": "7/7 finite-L K/radial bulk pass; scalar/angular/high-L final gates remain",
            },
            {
                "member": "historical_inner",
                "status": "REJECTED_AS_FINAL_MEMBER" if not transition.get("historical_inner_member_final_production_eligible", True) else "UNKNOWN",
                "reason": "replaced by one action-first strong-field transition on 0.70<=u<=0.715",
            },
        ])
    for key in ("regional_electric_member_2026_09_17", "zero_vector_epsilon_y_member_2026_09_18"):
        data = ledger.get(key, {})
        if isinstance(data, dict):
            reasons = data.get("reasons", data.get("reason", ""))
            if isinstance(reasons, list):
                reasons = "; ".join(str(x) for x in reasons)
            rows.append({"member": key, "status": str(data.get("status", "UNKNOWN")), "reason": str(reasons)})
    active = ledger.get("active_search", {})
    if isinstance(active, dict):
        rows.append(
            {
                "member": str(active.get("member", "active_search")),
                "status": str(active.get("status", "UNKNOWN")),
                "reason": f"next_gate={active.get('next_gate', 'UNKNOWN')}",
            }
        )
    return rows


def build_inventory(root: Path) -> dict[str, Any]:
    excluded = {".git", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache"}
    files = []
    for path in root.rglob("*"):
        if not path.is_file() or any(part in excluded for part in path.parts):
            continue
        files.append(path)
    suffixes = Counter((p.suffix.lower() or "<none>") for p in files)
    tops = Counter(p.relative_to(root).parts[0] for p in files if p.relative_to(root).parts)
    return {
        "total_files": len(files),
        "total_bytes": sum(p.stat().st_size for p in files),
        "by_suffix": dict(sorted(suffixes.items(), key=lambda kv: (-kv[1], kv[0]))),
        "by_top_level": dict(sorted(tops.items(), key=lambda kv: (-kv[1], kv[0]))),
        "important_files": [
            name
            for name in (
                "README.md",
                "RESEARCH_STATE.md",
                "ABSOLUTE_CLOSURE_LEDGER.json",
                "CURRENT_RELEASE_STATUS_2026-09-19.md",
                "CONTINUE_IMPLEMENTATION_2026-09-19.md",
                "MANIFEST.json",
                "SHA256SUMS",
            )
            if (root / name).exists()
        ],
    }


def reproduction_commands() -> list[str]:
    return [
        "python tools/show_repo.py",
        "python tools/show_gate_matrix.py",
        "python tools/show_member_matrix.py",
        "python tools/show_evidence_index.py",
        "python tools/show_repo_tree.py",
        "python tools/show_code_map.py",
        "python tools/run_all_evidence.py",
        "python tools/audit_zero_vector_light_ring.py",
        "python tools/audit_strong_field_transition.py",
        "python tools/audit_electric_hybrid_principal_feasibility.py",
        "python tools/audit_descriptor_pullback_equivalence.py",
        "python tools/audit_radial_descriptor_pencil.py",
        "python tools/audit_core_generalized_psi_descriptor.py",
        "python tools/audit_inner_y_hessian.py",
        "python tools/audit_regional_kinetic.py  # expected FAIL for rejected 2026-09-17 member",
        "python tools/audit_frozen_onshell.py    # expected FAIL for rejected zero-vector witness",
        "pytest -q",
        "python tools/release_manifest.py --check",
    ]


def build_repo_snapshot(root: Path) -> dict[str, Any]:
    ledger = _read_json(root / "ABSOLUTE_CLOSURE_LEDGER.json")
    packaging = _read_json(root / "PACKAGING_STATUS.json")
    working, _working_name = _working_status(root)
    grouped = working.get("tests", {}).get("grouped_complete_matrix", {}) if working else {}
    if grouped.get("status") == "PASS":
        software = f"{grouped.get('passed')}/{grouped.get('total')}_PASS_GROUPED_COMPLETE_MATRIX"
    else:
        software = ledger.get("software_tests", packaging.get("tests", {}).get("status", "UNKNOWN"))
    return {
        "absolute_full_closure": ledger.get("absolute_full_closure", "NOT_CERTIFIED"),
        "software_tests": software,
        "active_search": {
            "status": working.get("right_strong_field_transition", {}).get("status", "UNKNOWN"),
            "member": "RIGHT_STRONG_FIELD_TRANSITION_0p70_0p715",
            "next_gate": "rank-resolved action continuation"
        } if working else ledger.get("active_search", {}),
        "members": build_member_matrix(root),
        "gates": build_gate_matrix(root),
        "evidence": build_evidence_index(root),
        "inventory": build_inventory(root),
        "reproduce": reproduction_commands(),
    }


def _table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    if not rows:
        return "(none)"
    widths = {c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) for c in columns}
    widths = {c: min(w, 72) for c, w in widths.items()}
    header = " | ".join(c.ljust(widths[c]) for c in columns)
    sep = "-+-".join("-" * widths[c] for c in columns)
    body = []
    for row in rows:
        cells = []
        for c in columns:
            value = str(row.get(c, "")).replace("\n", " ")
            if len(value) > widths[c]:
                value = value[: max(0, widths[c] - 1)] + "…"
            cells.append(value.ljust(widths[c]))
        body.append(" | ".join(cells))
    return "\n".join([header, sep, *body])


def render_snapshot_text(snapshot: dict[str, Any], *, evidence_limit: int = 30) -> str:
    out = [
        "SSZ P5 REPOSITORY STATUS",
        "=" * 80,
        f"Absolute Full Closure : {snapshot['absolute_full_closure']}",
        f"Software tests        : {snapshot['software_tests']}",
        "",
        "MEMBERS",
        _table(snapshot["members"], ["member", "status", "reason"]),
        "",
        "GATES",
        _table(snapshot["gates"], ["gate", "status", "evidence"]),
        "",
        "EVIDENCE INDEX",
        _table(snapshot["evidence"][:evidence_limit], ["path", "status", "scope"]),
        "",
        "INVENTORY",
        json.dumps(snapshot["inventory"], indent=2),
        "",
        "REPRODUCTION COMMANDS",
        *(f"  {cmd}" for cmd in snapshot["reproduce"]),
    ]
    return "\n".join(out) + "\n"
