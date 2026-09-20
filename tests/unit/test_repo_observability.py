from __future__ import annotations

import json
from pathlib import Path

from ssz_p5.observability import (
    build_evidence_index,
    build_gate_matrix,
    build_inventory,
    build_member_matrix,
    build_repo_snapshot,
    render_snapshot_text,
)

ROOT = Path(__file__).resolve().parents[2]


def test_repo_snapshot_surfaces_current_truth():
    snapshot = build_repo_snapshot(ROOT)
    assert snapshot["absolute_full_closure"] == "NOT_CERTIFIED"
    assert any(row["status"] == "SEARCH_OPEN" for row in snapshot["members"])
    assert any(row["gate"] == "QNM" and "BLOCKED" in row["status"] for row in snapshot["gates"])
    assert snapshot["inventory"]["total_files"] > 500


def test_evidence_index_hashes_machine_reports():
    rows = build_evidence_index(ROOT)
    by_path = {row["path"]: row for row in rows}
    target = "data/generated/absolute_attempt_2026-09-19/ZERO_VECTOR_LIGHT_RING_AUDIT.json"
    assert target in by_path
    assert len(by_path[target]["sha256"]) == 64
    assert "ZERO_VECTOR" in by_path[target]["status"]


def test_gate_and_member_views_keep_rejections_separate_from_active_search():
    gates = build_gate_matrix(ROOT)
    members = build_member_matrix(ROOT)
    assert any(row["gate"] == "HORNDESKI_PRINCIPAL_FEASIBILITY" for row in gates)
    assert any(row["status"] == "REJECTED_FOR_ABSOLUTE_CLOSURE" for row in members)
    assert any(row["status"] == "REJECTED_AS_GLOBAL_ONSHELL_MEMBER" for row in members)
    assert any(row["status"] == "SEARCH_OPEN" for row in members)


def test_registry_declares_expected_rejection_exit_codes():
    registry = json.loads((ROOT / "REPO_EVIDENCE_REGISTRY.json").read_text())
    expected = {row["id"]: row["expected_exit"] for row in registry["audits"]}
    assert expected["regional_kinetic_rejected_member"] == 2
    assert expected["frozen_onshell_rejected_zero_vector"] == 2
    assert expected["electric_hybrid_principal_feasibility"] == 0
    text = render_snapshot_text(build_repo_snapshot(ROOT), evidence_limit=5)
    assert "ABSOLUTE_FULL_CLOSURE" in text
    assert "SEARCH_OPEN" in text
