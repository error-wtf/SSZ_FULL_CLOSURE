import json
from pathlib import Path

import numpy as np

from ssz_p5.production.action_handover import flat_step01

ROOT = Path(__file__).resolve().parents[2]
ATTEMPT = ROOT / 'data/generated/absolute_attempt_2026-09-19'


def test_central_single_action_checkpoint_is_positive_for_all_required_L():
    report = json.loads((ATTEMPT / 'ELECTRIC_HYBRID_ONSHELL_CENTRAL_AUDIT.json').read_text())
    assert report['status'] == 'PASS'
    assert len(report['scans']) == 7
    assert all(row['pass'] for row in report['scans'])
    assert all(row['negative_K_rows'] == 0 for row in report['scans'])
    assert all(row['negative_radial_rows'] == 0 for row in report['scans'])


def test_handover_checkpoint_rejects_f4_only_and_retains_rank_six():
    report = json.loads((ATTEMPT / 'ACTION_HANDOVER_CONTROL_SPACE_AUDIT.json').read_text())
    assert report['status'] == 'PENDING_JOINT_BACKGROUND_CONSTRAINED_SOLVE'
    for side in ('outer', 'inner'):
        assert report['one_control_a2_probe'][side]['accepted'] is False
        assert report['compact_six_primitive_rank'][side]['rank'] == 6
    assert report['certification']['joint_background_E00_E11_Ephi_JA_resolve'] is False
    assert report['certification']['direct41_handover_certified'] is False


def test_cinfinity_flat_step_endpoint_and_monotonic_contract():
    t = np.array([-1.0, 0.0, 1e-6, 0.25, 0.5, 0.75, 1.0 - 1e-6, 1.0, 2.0])
    w = flat_step01(t)
    assert w[0] == 0.0 and w[1] == 0.0
    assert w[-2] == 1.0 and w[-1] == 1.0
    assert np.all(np.diff(w) >= 0.0)
    assert abs(w[4] - 0.5) < 1e-15
