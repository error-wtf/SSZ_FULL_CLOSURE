"""Regression: structural invariants of the unsplit EOM core (STEP 1)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import sympy as sp
from ssz_p5.action import light_ring_identity as lri


def test_structural_invariants():
    inv = lri.structural_invariants()
    assert inv['svt_terms_carry_ap2']
    assert inv['one_f3free_term']
    assert inv['f3free_is_f2X_channel']
    assert inv['P_MH_leaves_f2X_channel']


def test_P_A_analytic_nonzero():
    r, f, h = lri.r, lri.f, lri.h
    assert sp.simplify(lri.P_A - sp.sqrt(h/f)*r**2) == 0
    # positive symbols => P_A > 0 structurally


def test_epsY_nulltest():
    T = sp.Function('T')
    terms = {'E00': lri.ap**2 * lri.f2Y * T(lri.r, lri.f, lri.h, lri.ph),
             'E11': lri.ap**2 * lri.f2Y * T(lri.r, lri.f, lri.h, lri.ph),
             'E_phi': lri.ap**2 * lri.f2Y * T(lri.r, lri.f, lri.h, lri.ph),
             'E_A': lri.ap**2 * lri.f2Y * T(lri.r, lri.f, lri.h, lri.ph)}
    zeros, all_vanish = lri.epsY_nulltest(terms)
    assert all_vanish and all(zeros.values())


def test_P_MH_leaves_f2X_channel():
    pmh = lri.P_MH(sp.expand(lri.E11 - lri.E00))
    assert pmh.has(lri.f2X) and not pmh.has(lri.f3, lri.f3X, lri.ap)
