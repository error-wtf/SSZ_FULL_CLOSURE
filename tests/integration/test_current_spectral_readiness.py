import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "data/generated/spectral_selection_2026-09-29"
    / "CURRENT_MEMBER_SPECTRAL_READINESS_AUDIT.json"
)


def test_current_member_spectral_readiness_audit():
    subprocess.run(
        [sys.executable, str(ROOT / "tools/audit_current_spectral_selection.py")],
        cwd=ROOT,
        check=True,
    )
    data = json.loads(REPORT.read_text())
    assert data["member_hash_match"] is True
    assert data["principal_status"] == "PRINCIPAL_CHARACTER_SMOOTH"
    assert all(row["pass"] for row in data["principal_audit"])
    assert min(
        row["min_adjacent_matched_overlap"]
        for row in data["principal_mode_tracking"]
    ) > 0.99
    promotion = data["attempted_global_same_action_promotion"]
    assert promotion["historical_or_q2_splicing_performed"] is False
    assert promotion["status"] == "GLOBAL_SAME_ACTION_BUILD_BLOCKED"
    assert data["spectral_weight_status"] == "NOT_YET_EVALUABLE"
