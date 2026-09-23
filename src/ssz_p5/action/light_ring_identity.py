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
