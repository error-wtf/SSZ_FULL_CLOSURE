"""Regression: KT2023 Maxwell-Horndeski primary-source transcription.

Sources transcribed verbatim (no fitting, no reconstruction):
- Kase & Tsujikawa, arXiv:2301.10362 (PRD 107, 104045 (2023)):
  Eqs. (7), (8), (9), (10), (24), (25), Appendix A Eq. (168), Eq. (85).
- Heisenberg & Tsujikawa, arXiv:1802.07035 (PLB 780, 638 (2018)):
  Eqs. (14), (15) (f-/h-variations, repo production normalization).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from ssz_p5.action import kt_mh_background as kt


def test_s0_slice_diagonal_map():
    """KT2023 component normalization vs HT2018 reduced-action normalization:
    KT_E00 = -HT_E00/(r^2 f) and KT_E11 = +HT_E11/(r^2 f) on the common
    subtheory S0 (Einstein + f2(X,F)).  This is the convention bridge between
    the repo's unsplit evaluator and the primary source."""
    m = kt.s0_slice_diagonal_map()
    assert m['map_holds'], m
    assert m['residual_E00_map'] == 0 and m['residual_E11_map'] == 0


def test_eq85_generating_two_term_identity():
    """C_bg RECIPE (machine-verified, luminal branch): in the project-locked
    luminal Horndeski branch the undivided Eq. (85) expression is an exact
    OFF-SHELL two-term identity  C85 == A*E00 + B*E22  with metric-only
    A = -r f^(3/2)/sqrt(h), B = +r f^(3/2)/sqrt(h).  This is the
    Eq.-85-generating combination the general C_bg assembly must use."""
    p = kt.probe_offshell_two_term_residual(kt.LUMINAL_G4PHI_SECTOR)
    assert p['two_term_identity_exact'] is True
    assert p['full_residual'] == 0
    import sympy as sp
    assert sp.simplify(p['A'] + kt.r * kt.f**sp.Rational(3, 2)
                       / sp.sqrt(kt.h)) == 0
    assert sp.simplify(p['B'] - kt.r * kt.f**sp.Rational(3, 2)
                       / sp.sqrt(kt.h)) == 0


def test_two_term_identity_sector_scope():
    """Sector scope: with general G4(phi,X) (G4X != 0) the metric-only
    two-term identity closes no longer; Eq. (85) is then an on-shell
    statement.  Recorded as scope evidence (project is luminal-locked)."""
    g4x_sector = dict(kt.LUMINAL_G4PHI_SECTOR)
    # re-open the G4X channel only (G4phiX/G4XX stay locked to keep the
    # probe cheap; the identity already fails with G4X alone)
    g4x_sector[kt.G4X] = kt.G4X
    p = kt.probe_offshell_two_term_residual(g4x_sector)
    assert p['two_term_identity_exact'] is False


def test_eq85_onshell_certificate_luminal_branch():
    """PRIMARY-SOURCE CERTIFICATE (machine-verified): in the project-locked
    luminal Horndeski branch (G4 = G4(phi), G3 = 0, G5 = 0) the undivided
    Eq. (85) expression vanishes exactly on the (E00, E22) shell.
    Premises carried in the result: pivot factors f>0, r>0, G4!=0, G4phi!=0,
    h!=0, and 2f - r f' != 0 (ring handled undivided elsewhere)."""
    cert = kt.onshell_eq85_certificate(kt.LUMINAL_G4PHI_SECTOR)
    assert cert['certificate_holds'] is True
    assert cert['undivided_eq85_onshell_residual'] == 0
    # premise tracking: pivot must be algebraically the documented factor
    # product (contains G4, G4phi, h and the light-ring denominator 2f-rf')
    import sympy as sp
    pivot = cert['pivot_D']
    ratio = sp.simplify(pivot * kt.f * kt.r / (kt.G4 * kt.G4phi * kt.h)
                        + (2 * kt.f - kt.r * kt.fp))
    assert ratio == 0, pivot


def test_a4_v8_slot_definitions():
    """KT2023 Appendix A Eq. (168) slot definitions must hold structurally:
    a4 = sqrt(fh)/2 * Hcal and v8 = G2F/(2 sqrt(fh))  (canonical MH v8 slot;
    the permanent mapping is ZK V9 = MH v8 — never the 13-slot V8 column)."""
    import sympy as sp
    assert sp.simplify(kt.a4_kt() - sp.sqrt(kt.f * kt.h) / 2 * kt.hcal()) == 0
    assert sp.simplify(kt.v8_kt() - kt.G2F / (2 * sp.sqrt(kt.f * kt.h))) == 0


def test_ward_span_derivation_status():
    """Delta22_SVT derivation record: the original 'no-go' was a non-holonomic
    jet-space artifact (d(f2)/dr treated as independent).  With the holonomic
    chain rule the Ward+scalar-EOM span closes all derivative channels with
    metric-only coefficients c0..c3 and the E22 map E22_HT = r^2 f * E22_KT;
    ONE unmatched free residual remains (evidence JSON) and must be matched
    before Delta22_SVT is complete.  This test pins the honest status."""
    p = kt.probe_ward_span_no_go()
    assert p["status"] == "SUPERSEDED_NON_HOLONOMIC_ARTIFACT"
    assert "holonomic" in p["df2_coefficients"]["holonomic"]
    assert "E22_HT = r^2 f * E22_KT" in p["correct_route"]
