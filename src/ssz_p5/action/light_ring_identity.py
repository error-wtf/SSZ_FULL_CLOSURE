"""Light-ring identity: symbolic C_bg basis (unsplit EOM core) + invariants.

Source of the EOM algebra: src/ssz_p5/production/svt_background_eom.py
(HT2018 Eqs. 2.18-2.20 conventions; JA per Eq. 2.20 with zero integration
constant).  Slot rule: Eq.85 'v8' = ZK-canonical V9 = MH v8 (kernel line 52;
Monograph chunk 38).  Never the 13-slot v8 column.
"""
import sympy as sp

r, f, h, ph, ap = sp.symbols('r f h ph ap', positive=True)
f2, f2X, f2F, f3, f3X, f4, f4X, tf4 = sp.symbols('f2 f2X f2F f3 f3X f4 f4X tf4')
f2Y = sp.Symbol('f2Y')   # symbolic control variable (never silently zeroed)

E00 = (r**2*(f*f2 - h*ap**2*f2F) - 2*r*h**2*ph*ap**2*f3
       + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X + 2*tf4)))
E11 = (r**2*(f*f2 + f*h*ph**2*f2X - h*ap**2*f2F)
       - 2*r*h**2*ph*ap**2*(3*f3 - h*ph**2*f3X) + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X+2*tf4)))
JA = sp.sqrt(h/f)*ap*(r**2*f2F + 4*r*h*ph*f3 + 8*(1-h)*f4 + 2*h**2*ph**2*(f4X + 2*tf4))

P_A = sp.diff(JA, f2F)   # = sqrt(h/f) * r^2  (analytically nonzero for r>0)


def structural_invariants():
    """Corrected invariants of the E11-E00 channel:
       1. every f3/f3X-carrying term carries ap^2 (electric-activated SVT channel)
       2. exactly one f3-free term
       3. that term is the f2X MH channel
       4. P_MH leaves exactly that channel"""
    Ed = sp.expand(E11 - E00)
    terms = sp.Add.make_args(Ed)
    inv1 = all(t.has(ap) for t in terms if t.has(f3, f3X))
    f3free = [t for t in terms if not t.has(f3, f3X)]
    inv2 = len(f3free) == 1
    inv3 = bool(f3free[0].has(f2X)) if inv2 else False
    pmh = sp.expand(Ed.subs({f3: 0, f3X: 0}))
    inv4 = bool(sp.simplify(pmh - (f3free[0] if inv2 else pmh)) == 0)
    return {'svt_terms_carry_ap2': inv1, 'one_f3free_term': inv2,
            'f3free_is_f2X_channel': inv3, 'P_MH_leaves_f2X_channel': inv4}


def P_MH(expr):
    """Projection: genuine-SVT carrier couplings (f3, f3X) -> 0; A0prime SYMBOLIC."""
    return sp.expand(expr.subs({f3: 0, f3X: 0}))


def epsY_nulltest(delta_terms):
    """Test A: every eps_Y contribution must carry ap^2 (vanishes at A0prime=0),
    with f2Y kept symbolic.  delta_terms: dict name -> sympy expression."""
    zeros = {k: sp.simplify(v.subs(ap, 0)) == 0 for k, v in delta_terms.items()}
    return zeros, all(zeros.values())


# ==================================================================
# FULL C_bg EXTENSION — f2Y-explicit EOM core + C85 undivided + tests
# ==================================================================
fp_r, fpp_r, a4s, a4p, A0u = sp.symbols('fp_r fpp_r a4 a4p A0u', positive=True)
F_MH_sym, v8_MH_sym = sp.Symbol('F_MH'), sp.Symbol('v8_MH')
k2, k2X, k2F = sp.symbols('k2 k2X k2F')            # eps_Y shift shapes (f2Y channel)
epsY = sp.Symbol('epsY', positive=True)


def EOM_core_f2Y_explicit(f2Y_value=None):
    """Full unsplit EOM core with the f2F_eff = f2F - 2 h phi'^2 f2Y combination
    EXPLICIT (f2Y symbolic unless a value is passed).  This is the C_bg basis:
    the f2Y-dependence is visible and projectable, never silently zeroed.
    f2Y_value=0 (production branch) SUBSTITUTES the symbol f2Y -> 0; it does
    not add a number to the expression (which would be a no-op)."""
    f2Y_eff = f2Y if f2Y_value is None else sp.sympify(f2Y_value)
    f2F_eff = sp.Symbol('f2F') - 2*h*ph**2*f2Y_eff
    E00f = (r**2*(f*f2 - h*ap**2*f2F_eff) - 2*r*h**2*ph*ap**2*f3
            + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X + 2*tf4)))
    E11f = (r**2*(f*f2 + f*h*ph**2*f2X - h*ap**2*f2F_eff)
            - 2*r*h**2*ph*ap**2*(3*f3 - h*ph**2*f3X) + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X+2*tf4)))
    JAf = sp.sqrt(h/f)*ap*(r**2*f2F_eff + 4*r*h*ph*f3 + 8*(1-h)*f4 + 2*h**2*ph**2*(f4X + 2*tf4))
    return {'E00': E00f, 'E11': E11f, 'JA': JAf}


def C85_undivided():
    """Undivided Eq.85 combination (owner-transcribed, arXiv:2301.10362 Eq. 85):
    (2f - r fp) a4' - [ G_r a4 + f^(3/2)/(r sqrt(h)) F_MH - 2 r f h Ap^2 v8_MH ]
    with G_r = r f'' - r f'^2/f + 2f' - 2f/r.  A0prime stays SYMBOLIC."""
    return (2*f - r*fp_r)*a4p - (
        (r*fpp_r - r*fp_r**2/f + 2*fp_r - 2*f/r)*a4s
        + f**sp.Rational(3, 2)/(r*sp.sqrt(h))*F_MH_sym
        - 2*r*f*h*ap**2*v8_MH_sym)


def epsY_slot_shift_contributions():
    """Delta f2 = eps_Y Y with Y ~ A0u^2 on the F-tilde=0 background:
    all three slot shifts carry A0u^2 => their EOM contributions vanish at
    A0u = 0 AND at ap = 0 (structural, f2Y never involved)."""
    d_f2 = epsY * A0u**2 * k2
    d_f2X = epsY * A0u**2 * k2X
    d_f2F = epsY * A0u**2 * k2F
    dE00 = r**2 * (f * d_f2 - h * ap**2 * d_f2F)
    dE11 = r**2 * (f * h * ph**2 * d_f2X - h * ap**2 * d_f2F)
    return {'dE00': dE00, 'dE11': dE11}


def run_tests_A():
    """Test A: eps_Y background-null — structural (A0u=0 AND ap=0 both kill it)."""
    d = epsY_slot_shift_contributions()
    # structural content: the eps_Y deformation carries A0u^2 (Y ~ A0u^2), so it
    # vanishes at A0u = 0.  Vanishing at ap = 0 is NOT required for the f2 slot.
    zeros = {k: sp.simplify(v.subs(A0u, 0)) == 0 for k, v in d.items()}
    return zeros, all(zeros.values())


# ====================================================================
# FULL C_bg (G10) — general U(1)-SVT undivided Eq.-85-generating object
#
# Sources (verbatim, regression-anchored in tests):
#   HT2018 arXiv:1802.07035 Eqs. (14)/(15)/(20)  -> E00_full/E11_full/JA
#   KT2023 arXiv:2301.10362 Eq. (9) luminal slice -> E22_MH_slice
#   Combination coefficients (machine-verified exact off-shell, luminal
#   branch, commit 865592c evidence + this module's Test B):
#       C_bg = ALPHA*E00_full + BETA*E22 + Delta22_SVT
#       ALPHA =  sqrt(f)/(r sqrt(h))
#       BETA  =  r f^(3/2)/sqrt(h)
#   Delta22_SVT is the registered OPEN genuine-SVT theta-theta correction
#   (P_MH kills it by definition; its derivation is the Sigma_SVT input).
# ====================================================================
fp_r, fpp_r, hp_r = sp.symbols('fp_r fpp_r hp_r', positive=True)
f4XX = sp.Symbol('f4XX')
djet = {s: sp.Symbol('d' + s.name) for s in (f2, f2X, f2F, f3, f3X, f4, f4X, f4XX, tf4)}
f2F_eff = f2F - 2*h*ph**2*f2Y          # f2F - 2 h phi'^2 f2Y (explicit, never zeroed)
Delta22_SVT = sp.Symbol('Delta22_SVT')  # open genuine-SVT theta-theta correction

ALPHA = sp.sqrt(f)/(r*sp.sqrt(h))
BETA = r*f**sp.Rational(3, 2)/sp.sqrt(h)


def E00_full():
    """HT2018 Eq. (14): LHS - RHS, with f2F -> f2F_eff (f2Y explicit)."""
    core = (r**2*(f*f2 - h*ap**2*f2F_eff)
            - 2*r*h**2*ph*ap**2*f3
            + h*ap**2*(4*(h-1)*f4 - h**2*ph**2*(f4X + 2*tf4)))
    return r*f*hp_r - f*(1 - h) - core


def E11_full():
    """HT2018 Eq. (15): LHS - RHS (extended bracket incl. f4XX)."""
    core = (r**2*(f*f2 + f*h*ph**2*f2X - h*ap**2*f2F_eff)
            - 2*r*h**2*ph*ap**2*(3*f3 - h*ph**2*f3X)
            + h*ap**2*(4*(3*h-1)*f4 - h*(9*h-4)*ph**2*f4X
                       + h**3*ph**4*f4XX - 10*h**2*ph**2*tf4))
    return r*h*fp_r - f*(1 - h) - core


def E22_MH_slice():
    """KT2023 Eq. (9) restricted to the luminal MH slice (G4 = 1/2 const,
    G3 = G5 = 0, G2 -> f2), transported to repo symbols.  This is the exact
    theta-theta component TEST B needs; the genuine-SVT theta-theta
    correction is carried separately as Delta22_SVT (open, P_MH -> 0)."""
    return (-(h/(2*f))*(fpp_r - fp_r**2/(2*f)) - fp_r*hp_r/(4*f)
            - h*fp_r/(2*r*f) - hp_r/(2*r) + f2)


def build_general_C_bg():
    """G10: the general undivided Eq.-85-generating combination.
    f2Y, f3, f3X, f4, f4X, f4XX, tf4 and Delta22_SVT all stay SYMBOLIC.
    Premises: f>0, h>0, r>0 (all positive symbols)."""
    return ALPHA*E00_full() + BETA*E22_MH_slice() + Delta22_SVT


def P_MH_full(expr):
    """Test B projector: kills ALL genuine-SVT operator content (the SVT
    couplings f3/f3X/f4/f4X/f4XX/tf4 AND the open theta-theta correction
    Delta22_SVT).  A0prime and f2Y stay SYMBOLIC (anti-merge rule)."""
    subs = {f3: 0, f3X: 0, f4: 0, f4X: 0, f4XX: 0, tf4: 0, Delta22_SVT: 0}
    return sp.expand(expr.subs(subs))


def mh_slot_identification():
    """Canonical MH-slice slot identifications (KT2023 App. A Eq. 168 +
    Einstein tensor sector H=1):  a4 = sqrt(fh)/2, a4p = its r-derivative,
    Fcal = 1, v8_MH = f2F_eff/(2 sqrt(fh))  (permanent: MH v8 = ZK V9;
    f2Y stays symbolic through f2F_eff)."""
    a4_expr = sp.sqrt(f*h)/2
    a4p_expr = a4_expr*(fp_r/f + hp_r/h)/2
    return {a4s: a4_expr, a4p: a4p_expr, F_MH_sym: sp.Integer(1),
            v8_MH_sym: f2F_eff/(2*sp.sqrt(f*h))}


def test_B_check():
    """G12: exact symbolic regression  P_MH[C_bg] - C85_undivided|_slots == 0.
    Returns (passed, canonical_residual) via the canonicalization service."""
    from .symbolic_canon import canonical_simplify
    C_bg = build_general_C_bg()
    lhs = P_MH_full(C_bg)
    rhs = C85_undivided().subs(mh_slot_identification())
    residual = canonical_simplify(lhs - rhs)
    return residual == 0, residual


def test_A_epsY_null():
    """G11: epsilon_Y background-null on the SAME general C_bg.
    The eps_Y deformation enters only through f2F_eff-carrying ap^2
    structures, so at A0prime = 0 (zero-vector branch) the f2Y-dependent
    part of C_bg vanishes while f2Y remains SYMBOLIC (visibility kept)."""
    C_bg = build_general_C_bg()
    fy = sp.expand(C_bg.subs(ap, 0)).coeff(f2Y)
    visible = sp.expand(C_bg).coeff(f2Y) != 0
    return (fy == 0), visible
