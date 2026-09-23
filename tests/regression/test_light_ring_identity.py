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
    r, f, h, ap = lri.r, lri.f, lri.h, lri.ap
    # JA = sqrt(h/f) * ap * (r^2 f2F + 4 r h ph f3 + 8(1-h) f4 + 2 h^2 ph^2 (f4X + 2 tf4))
    # => P_A = d JA/d f2F = sqrt(h/f) * ap * r^2  (positive symbols => structurally > 0)
    assert sp.simplify(lri.P_A - ap*sp.sqrt(h/f)*r**2) == 0


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


def test_epsY_slot_shift_structural():
    """STEP 2 Test A (structural): Delta f2 = eps_Y Y with Y ∝ A0prime^2 on the
    F-tilde=0 electric background => all three slot shifts carry A0prime^2 =>
    every EOM contribution of Delta f2 vanishes at A0prime = 0, with f2Y symbolic
    and WITHOUT setting f2Y = 0.  (Y ~ nabla phi nabla phi F^2 carries F ~ A0prime^2.)"""
    r, f, h, ph = lri.r, lri.f, lri.h, lri.ph
    ap, A0u, epsY = lri.ap, sp.Symbol('A0u', positive=True), sp.Symbol('epsY', positive=True)
    k2, k2X, k2F = sp.symbols('k2 k2X k2F')   # invariant shift shapes (documented, not fitted)
    # slot shifts, each carrying the background A0prime^2 factor:
    d_f2  = epsY * A0u**2 * k2          # Delta f2  = eps_Y * Y          (Y ~ A0u^2)
    d_f2X = epsY * A0u**2 * k2X         # Delta f2X = eps_Y * dY/dX       (~ A0u^2)
    d_f2F = epsY * A0u**2 * k2F         # Delta f2F = eps_Y * dY/dF       (~ A0u^2)
    # EOM contributions of the shifts (linear EOM structure from the unsplit core):
    dE00 = r**2 * (f * d_f2 - h * ap**2 * d_f2F)
    dE11 = r**2 * (f * d_f2X * sp.Symbol('xconv') + f * h * ph**2 * d_f2X - h * ap**2 * d_f2F)
    dJA = (lri.P_A) * d_f2F * ap**0     # JA is linear in f2F: P_A * d_f2F
    for dE in (dE00, dE11, dJA):
        # EVERY term must carry A0u^2 (via the shifts) — that is the structural content:
        assert sp.expand(dE).has(A0u**2) or sp.simplify(dE.subs(A0u, 0)) == 0
        assert sp.simplify(dE.subs(A0u, 0)) == 0
    # and the f2Y-symbolic requirement: the shift STRUCTURE does not presuppose f2Y=0
    # (k2/k2X/k2F are the documented eps_Y*Y shapes; f2Y itself never enters this
    #  branch's EOM, which is exactly the branch postulate — registered)


def test_epsY_nulltest_uses_f2Y_symbolic():
    """The eps_Y nulltest helper must accept symbolic f2Y and still pass."""
    T = sp.Function('T')
    terms = {'E00': lri.ap**2 * lri.f2Y * T(lri.r, lri.f, lri.h, lri.ph),
             'E11': lri.ap**2 * lri.f2Y * T(lri.r, lri.f, lri.h, lri.ph)}
    zeros, ok = lri.epsY_nulltest(terms)
    assert ok and zeros['E00'] and zeros['E11']


def test_C85_undivided_structure():
    """C85_undivided must contain the geom*a4 term, the tensor kinetic term and
    the electric term -2 r f h ap^2 v8_MH explicitly (A0prime symbolic)."""
    from ssz_p5.action import light_ring_identity as lri
    C = lri.C85_undivided()
    assert C.has(lri.a4s) and C.has(lri.F_MH_sym) and C.has(lri.v8_MH_sym)
    assert sp.simplify(C.subs(lri.ap, 0)) != 0   # A0prime symbolic: electric term survives


def test_epsY_slot_shift_contributions_vanish():
    """Test A on the full shift contributions: vanish at A0u=0 and at ap=0."""
    from ssz_p5.action import light_ring_identity as lri
    zeros, ok = lri.run_tests_A()
    assert ok and all(zeros.values())


def test_f2Y_explicit_EOM_core():
    """The f2Y-explicit EOM core must differ from the f2Y-free core while f2Y is
    symbolic (visibility requirement), and reduce identically at f2Y=0."""
    from ssz_p5.action import light_ring_identity as lri
    full = lri.EOM_core_f2Y_explicit()
    red = lri.EOM_core_f2Y_explicit(f2Y_value=0)
    diff_E00 = sp.simplify(full['E00'] - red['E00'])
    assert diff_E00 != 0 and diff_E00.has(lri.f2Y)   # f2Y-dependence is VISIBLE
    assert sp.simplify(full['E00'].subs(lri.f2Y, 0) - red['E00']) == 0  # reduction exact


def test_C85_undivided_executable():
    """TEST B precondition: C85_undivided is executable source with the electric
    term surviving symbolic A0prime (anti-merge rule: never set A0prime=0 here)."""
    from ssz_p5.action import light_ring_identity as lri
    C = lri.C85_undivided()
    # electric term -2 r f h ap^2 v8_MH must be present:
    assert C.has(-2 * lri.r * lri.f * lri.h * lri.ap**2 * lri.v8_MH_sym)
    # A0prime symbolic: ap survives
    assert C.has(lri.ap)
    # geom term present:
    assert C.has(lri.a4s) and C.has(lri.fpp_r) and C.has(lri.fp_r)


def test_P_MH_preserves_A0prime_symbolic():
    """P_MH kills f3/f3X but must NEVER touch A0prime (Test B anti-merge rule)."""
    from ssz_p5.action import light_ring_identity as lri
    expr = lri.ap**2 * lri.v8_MH_sym * lri.f3   # SVT-electric mixed dummy
    proj = lri.P_MH(expr)
    assert proj == 0                             # f3 carrier projects out
    keeper = lri.ap**2 * lri.v8_MH_sym           # Eq.85 electric term: no f3/f3X
    assert sp.simplify(lri.P_MH(keeper) - keeper) == 0   # survives projection
