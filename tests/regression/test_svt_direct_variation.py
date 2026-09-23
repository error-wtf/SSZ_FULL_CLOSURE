"""Regression: direct theta-theta variation of the genuine SVT action.

Pins the authoritative Delta22_SVT derivation (STEP3):
  - exact normalization-channel checks V1/V2/V3 against the frozen repo
    transcriptions E00_full / E11_full / E22_MH_slice,
  - P_MH[Delta22_SVT] == 0 and second-order-EOM property,
  - areal-gauge guard (C = r^2 only AFTER the Euler-Lagrange step),
  - Test B / Test A exactness on the explicit C_bg assembly,
  - projector backward compatibility,
  - recorded Ward-oracle consistency status.
"""
import json
from pathlib import Path

import sympy as sp

from ssz_p5.action import light_ring_identity as lri
from ssz_p5.action import svt_direct_variation as sdv
from ssz_p5.action.symbolic_canon import canonical_simplify


def test_delta22_normalization_channels_exact():
    """V1/V2/V3: the reduced C(r)-metric action reproduces the frozen repo
    transcriptions exactly; the map factors are metric-only."""
    rep = sdv.verify_normalization_channels()
    assert rep['V1_residual'] == 0
    assert rep['V2_residual'] == 0
    assert rep['V3_residual'] == 0
    assert rep['V2_kappa_metric_only'] is True
    assert rep['V3_kappa_metric_only'] is True


def test_delta22_p_mh_zero_and_second_order():
    """P_MH kills the derived correction completely; no third-order field
    channels (h''/phi'''/A0''') appear (second-order-EOM property)."""
    d22 = sdv.build_delta22_svt(verify=False)
    assert sdv.p_mh_delta22_zero(d22) is True
    assert sdv.second_order_field_check(d22) is True
    # electric activation: vanishes identically at A0prime = 0
    assert sp.simplify(d22.subs(lri.ap, 0)) == 0


def test_areal_gauge_guard():
    """The areal gauge C = r^2 must not be imposed before the variation:
    the raw EL channels carry the symbolic C(r); only the final expressions
    are gauge-specialized."""
    ch = sdv.delta22_channels()
    assert ch['E_C_SVT'].has(sdv.CF), 'EL_C must be formed with C(r) symbolic'
    fixed = sdv.to_repo_symbols(ch['E_C_SVT'], areal_gauge=True)
    assert not fixed.has(sdv.CF), 'areal gauge substitution failed'
    unfixed = sdv.to_repo_symbols(ch['E_C_SVT'], areal_gauge=False)
    assert unfixed.has(sdv.CF), 'areal_gauge=False must keep C symbolic'


def test_delta22_operator_families():
    """The correction factorizes by genuine-SVT operator families only."""
    d22 = sdv.build_delta22_svt(verify=False)
    for j in (lri.f3, lri.f3X, lri.f4, lri.f4X, lri.f4XX, lri.tf4):
        assert d22.has(j), f'operator family {j} missing from the derivation'
    # no f2-sector leakage: the f2-family must not appear
    for j in (lri.f2, lri.f2X, lri.f2F, lri.f2Y):
        assert not d22.has(j), f'MH-slice operator {j} leaked into Delta22'


def test_G12_testB_explicit_delta22_exact():
    """Test B on the EXPLICIT C_bg assembly: P_MH[C_bg] - C85|_slots == 0."""
    C_bg = lri.build_general_C_bg(delta22='explicit')
    lhs = lri.P_MH_full(C_bg)
    rhs = lri.C85_undivided().subs(lri.mh_slot_identification())
    assert canonical_simplify(lhs - rhs) == 0


def test_G11_testA_epsY_null_explicit():
    """Test A on the EXPLICIT C_bg assembly: f2Y-channel vanishes at
    A0prime = 0 while staying symbolic off-shell."""
    C_bg = lri.build_general_C_bg(delta22='explicit')
    fy = sp.expand(C_bg.subs(lri.ap, 0)).coeff(lri.f2Y)
    assert fy == 0
    assert sp.expand(C_bg).coeff(lri.f2Y) != 0


def test_E22_full_decomposition_identity():
    """E22_full == P_MH[E22_full] + Delta22_SVT with the derived correction."""
    E22_full = lri.E22_MH_slice() + lri.get_delta22_svt_expr()
    assert sp.simplify(lri.P_MH_full(E22_full) - lri.E22_MH_slice()) == 0


def test_p_mh_projector_backward_compatible():
    """The extended projector (phi-channel jet kills) is a superset: on
    expressions without phi-channel symbols it projects identically to the
    historical kill list."""
    old = lri.build_general_C_bg(delta22='symbol')
    subs_old = {lri.f3: 0, lri.f3X: 0, lri.f4: 0, lri.f4X: 0, lri.f4XX: 0,
                lri.tf4: 0, lri.Delta22_SVT: 0}
    assert sp.simplify(lri.P_MH_full(old) - sp.expand(old.subs(subs_old))) == 0


def test_step3_evidence_complete():
    """The STEP3 evidence artifact must record a fully passing derivation."""
    p = Path(__file__).resolve().parents[2] / 'data' / 'generated' / 'phase2_q2' \
        / 'STEP3_DELTA22_DIRECT_VARIATION.json'
    assert p.exists(), 'STEP3 evidence missing - run tools/derive_svt_reduced_action.py'
    d = json.loads(p.read_text())
    assert d['all_checks_pass'] is True
    for k, v in d['reduced_piece_equality_with_module'].items():
        assert v is True, f'reduced piece {k} not equal to covariant derivation'
    for k in ('V1_residual', 'V2_residual', 'V3_residual',
              'P_MH_delta22_zero', 'second_order_fields'):
        assert d['normalization_channel_checks'][k] is True


def test_ward_oracle_record_consistency():
    """Superseded Ward-oracle record (STEP1 JSON) vs the pinned decomposition
    R = R_el + R_met (light_ring_identity.ward_residual_decomposition):

    * the pure-metric parts agree EXACTLY;
    * the recorded electric part carries a pre-c4-convention normalization:
      recorded_electric - R_el|_{f2Y=0} == +ap^2 f2F hp_r r^3 (h-1)/(4h)
      (term-wise 1/h on the f hp_r channel).

    The recorded residual therefore cannot serve as a definition of E22; it is
    retained as a historical record whose exact deviation this test pins."""
    d = lri.ward_residual_decomposition()
    # recorded STEP1 unmatched residual (verbatim from
    # data/generated/phase2_q2/STEP1_DELTA22_WARD_SPAN_SOLUTION.json)
    r_rec = -lri.r*(lri.ap**2*lri.f*lri.f2F*lri.hp_r*lri.r**2
                    - lri.ap**2*lri.f2F*lri.fp_r*lri.h**2*lri.r**2
                    + 2*lri.f**2*lri.h*lri.hp_r + lri.f**2*lri.hp_r**2*lri.r
                    - lri.fp_r**2*lri.h**2*lri.r)/(4*lri.f*lri.h)
    # 1) pure-metric parts agree exactly
    met_part = -lri.r*(2*lri.f**2*lri.h*lri.hp_r + lri.r*(lri.f**2*lri.hp_r**2
                       - lri.fp_r**2*lri.h**2))/(4*lri.f*lri.h)
    assert sp.simplify(d['R_met'] - met_part) == 0
    # 2) the electric deviation from the pinned decomposition is exactly the
    #    documented pre-c4-convention artifact
    electric_rec = sp.simplify(r_rec - met_part)
    deviation = sp.simplify(electric_rec - d['R_el'].subs(lri.f2Y, 0))
    expected = lri.ap**2*lri.f2F*lri.hp_r*lri.r**3*(lri.h - 1)/(4*lri.h)
    assert sp.simplify(deviation - expected) == 0, \
        'recorded electric part drifted from the registered deviation formula'
    # 3) the pinned decomposition itself remains algebraically exact
    assert sp.simplify(d['R_el'] - lri.c4_JA_channel*lri.JA_channel) == 0
