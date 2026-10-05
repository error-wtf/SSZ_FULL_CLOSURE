"""Regression: the real frozen N2 artifact cannot open a scientific gate."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
PROFILE = ROOT / 'data/generated/spectral/N2_LUMINAL_BACKGROUND_PROFILE_V2.csv'
V3_RESULT = ROOT / 'data/generated/spectral/N2_LUMINAL_SOLVE_RESULT_V3.json'
V3_PROFILE = ROOT / 'data/generated/spectral/N2_LUMINAL_BACKGROUND_PROFILE_V3.csv'
V3_DECISION = ROOT / 'data/generated/spectral/N3_N4_LUMINAL_DECISION_V3.json'


def test_real_runner_blocks_invalid_background_and_preserves_all_rows(tmp_path):
    before = hashlib.sha256(PROFILE.read_bytes()).hexdigest()
    out = tmp_path / 'decision.json'
    result = subprocess.run(
        [sys.executable, str(ROOT / 'tools/run_n3_n4_luminal_decision.py'),
         '--output', str(out)], cwd=ROOT, env={**os.environ, 'PYTHONPATH': str(ROOT / 'src')},
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    report = json.loads(out.read_text())
    assert report['verdict'] == 'BLOCKED'
    assert report['physical_qnm_claim_allowed'] is False
    assert report['profile_rows'] == report['diagnostic_rows'] == 200
    assert report['n3']['status'] == 'NOT_EVALUABLE'
    assert report['n4']['status'] == 'NOT_EVALUABLE'
    assert report['n4']['executed_L'] == []
    assert report['n4']['requested_L'] == list(range(6, 1001))
    assert not out.with_suffix('.stream41.csv').exists()
    assert hashlib.sha256(PROFILE.read_bytes()).hexdigest() == before


def load_runner():
    import importlib.util
    path = ROOT / 'tools/run_n3_n4_luminal_decision.py'
    spec = importlib.util.spec_from_file_location('n3_audit', path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


def test_kinetic_precheck_rejects_partial_domain_even_if_positive():
    runner = load_runner()
    assert runner.kinetic_precheck([1.0] * 54, expected_rows=200,
                                    certified_error_bound=0.0)['status'] == 'NOT_EVALUABLE'


# ---------------------------------------------------------------------------
# V3 coordinate-convention fail-closed gates (N2_COORDINATE_MISMATCH repair)
# ---------------------------------------------------------------------------

def load_v3_solver():
    import importlib.util
    path = ROOT / 'tools/run_luminal_background_solve_v3.py'
    spec = importlib.util.spec_from_file_location('n2_v3', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_require_dr_convention_rejects_du_columns():
    """A profile whose stored 'phi_r' is actually dphi/du (the V2 defect)
    must be refused: d/dr = -u^2 d/du differs by the factor -u^2."""
    mod = load_v3_solver()
    u = np.linspace(0.71, 100.0, 200)
    phi = 1.0 - 0.15 * np.exp(-(u - 0.71) / 5.0)  # phi(r=1/u), monotone in u
    dphi_du = np.gradient(phi, u)
    with pytest.raises(ValueError, match='d/dr'):
        mod.require_dr_convention(u, phi, dphi_du)          # d/du as-is
    with pytest.raises(ValueError, match='d/dr'):
        mod.require_dr_convention(u, phi, +u**2 * dphi_du)  # wrong sign
    # the chain-rule value is the ONLY accepted column
    mod.require_dr_convention(u, phi, -u**2 * dphi_du)


def test_require_dr_convention_rejects_unrelated_columns():
    mod = load_v3_solver()
    u = np.linspace(0.71, 100.0, 200)
    phi = np.ones(200)
    with pytest.raises(ValueError):
        mod.require_dr_convention(u, phi, np.full(200, -1.28))  # constant != 0


def test_require_static_radial_X_rejects_wrong_convention():
    """X = +phi_u^2/(2f) (V2) or any column violating h phi_r^2 + 2X = 0
    must be refused (N2_WRONG_X repair)."""
    mod = load_v3_solver()
    r = np.linspace(0.01, 1.4, 100)
    h = 0.4 + 0.6 * r
    phi_r = -0.5 * np.ones(100)
    X_correct = -h * phi_r**2 / 2
    mod.require_static_radial_X(h, phi_r, X_correct)
    with pytest.raises(ValueError, match='h\\*phi_r\\^2 \\+ 2X'):
        mod.require_static_radial_X(h, phi_r, +phi_r**2 / 2)   # V2 sign flip
    with pytest.raises(ValueError, match='h\\*phi_r\\^2 \\+ 2X'):
        mod.require_static_radial_X(h, phi_r, X_correct * 1.5)  # arbitrary


def test_v3_solver_rejects_the_frozen_v2_profile_fail_closed():
    """End-to-end: feeding the V2 profile (d/du derivatives, wrong X) into
    the V3 gates raises — no silent consumption of the defective artifact."""
    mod = load_v3_solver()
    d = pd.read_csv(PROFILE).sort_values('u').reset_index(drop=True)
    with pytest.raises(ValueError):
        # V2's phi_p column is d/du (numpy.gradient over u): must be refused
        mod.require_dr_convention(d.u.to_numpy(float),
                                  d.phi.to_numpy(float),
                                  d.phi_p.to_numpy(float))


def test_v3_symbolic_core_if_result_exists():
    """If the V3 run has been executed (artifact present), its recorded
    symbolic facts must include the repaired structure."""
    if not V3_RESULT.exists():
        pytest.skip('V3 result not generated yet')
    res = json.loads(V3_RESULT.read_text())
    facts = res['symbolic_facts']
    assert facts['h_prime_from_E00'] is True
    assert facts['f_prime_from_E11'] is True
    assert facts['E22_bianchi_shadow_identically_zero'] is True
    assert facts['K_identity'] is True
    assert facts['K_on_family_is_J0_squared_form'] is True
    # V2's E11 defect is exactly the missing -h phi'^2 term
    assert facts['V2_E11_missing_term'] == '-h*ph**2'
    assert res['derivative_convention'].startswith('d/dr at r = r_over_rs = 1/u')


def test_v3_profile_states_derivative_convention_fail_closed():
    """The exported background artifact must state its derivative convention
    explicitly (consumers reduce_profile/emitters expect x = r)."""
    if not V3_PROFILE.exists():
        pytest.skip('V3 profile not generated yet')
    d = pd.read_csv(V3_PROFILE)
    assert 'derivative_convention' in d.columns
    conv = d.derivative_convention.iloc[0]
    assert 'd/dr' in conv and '1/u' in conv
    # full domain, no window truncation
    u = d.u.to_numpy(float)
    assert u.min() <= 0.71 + 1e-9 and u.max() >= 100.0 - 1e-9
    # r strictly increasing (r_over_rs column) for the reduced-operator path
    assert np.all(np.diff(d.r_over_rs.to_numpy(float)) > 0)


def test_v3_no_go_and_fail_closed_decision_if_present():
    if not V3_DECISION.exists():
        pytest.skip('V3 decision not generated yet')
    dec = json.loads(V3_DECISION.read_text())
    # fail-closed: no verdict names outside the declared rule set
    assert dec['verdict'] in ('BLOCKED',)
    assert dec['n3']['status'] == 'NOT_EVALUABLE'
    assert dec['n4']['status'] == 'NOT_EVALUABLE'
    assert dec['n4']['physical_qnm_claim_allowed'] is False
    assert dec['n4']['executed_L'] == []
    allowed = {'H1_CONFIRMED', 'H3_PHYSICAL_GHOST', 'NEW_PATHOLOGY_MAP'}
    for rule in dec['n4']['declared_rules_unchanged_echo']:
        assert rule['id'] in allowed
    # the no-go must be quantified, not asserted
    nogo = dec['no_go_certificate']
    assert nogo['incompatibility_factor'] > 100.0
    assert nogo['J0_forced_by_archive_core_sq'] > 0.0
