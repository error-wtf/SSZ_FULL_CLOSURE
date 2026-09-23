"""Direct theta-theta variation of the genuine U(1)-SVT action (HT2018).

Authoritative derivation of Delta22_SVT (from commit f4d5ce2 route decision):

    Delta22_SVT = kappa_C * EL_C[L_SVT] |_{C = r^2},   kappa_C = sqrt(h/f),

with kappa_C pinned EXACTLY on the verified MH slice (V3 below) and the
angular metric degree of freedom kept UNFIXED during the variation:

    ds^2 = -f(r) dt^2 + dr^2/h(r) + C(r) dOmega^2,   g_thetatheta = C(r).

The areal gauge C = r^2 is substituted ONLY after the Euler-Lagrange
operator is formed (guard in build_delta22_svt).

Reduced action (Heisenberg & Tsujikawa, arXiv:1802.07035, Eqs. (1)-(13)):

    L_EH  = (1/2) sqrt(-g) R                        [M_pl^2 = 1]
    L_f2  = sqrt(-g) f2(phi, X, F)                  [Y = 4FX substituted;
                                                     f2F_eff = f2F - 2 h phi'^2 f2Y]
    L_3   = sqrt(-g) f3(phi,X) g_rs Ftilde^{mu r} Ftilde^{nu s} nabla_mu nabla_nu phi
    L_4a  = sqrt(-g) (f4X/2 + tf4) Ftilde Ftilde (nabla nabla phi)^2
    L_4b  = sqrt(-g) f4 L^{mu nu ab} F_{mu nu} F_{ab}

The tf3 part of M3 vanishes identically on the static electric background
(structurally verified; see tools/derive_svt_reduced_action.py, which derives
all five reduced pieces from the covariant 4D action and asserts equality
with the transcribed forms below).

Normalization map to the repository E22 convention (exact, metric-only map
factors - verified, never fitted):

    V1: E00_full     == -2 f^(3/2) sqrt(h) * EL_f[L_full]  |_{C=r^2}
    V2: E11_full     == +2 sqrt(f) h^(3/2) * EL_h[L_full]  |_{C=r^2}
    V3: E22_MH_slice == sqrt(h/f)          * EL_C[L_MH]   |_{C=r^2}
"""
from __future__ import annotations

import sympy as sp

from .light_ring_identity import (
    r, f, h, ph, ap, fp_r, fpp_r, hp_r,
    f2, f2X, f2F, f2Y, f3, f3X, f4, f4X, f4XX, tf4,
)

# ------------------------------------------------------------------ symbols
# f, h as canonical repo SYMBOLS (imported); derivation-space FUNCTIONS:
ff, hf, CF = sp.Function('f')(r), sp.Function('h')(r), sp.Function('C')(r)
phif, A0f = sp.Function('phif')(r), sp.Function('A0f')(r)

# fields' radial derivatives (holonomic via sympy Function machinery)
Xe = -hf*sp.diff(phif, r)**2/2
Fe = hf*sp.diff(A0f, r)**2/(2*ff)

# jets as multivariate functions (explicit (phi, X[, F]) arguments -> exact chains)
f2eff_f = sp.Function('f2eff')(phif, Xe, Fe)   # repo convention: Y=4FX substituted,
                                               # f2Y kept symbolic via f2F_eff
f3f = sp.Function('f3')(phif, Xe)
f4f = sp.Function('f4')(phif, Xe)
f4Xf = sp.Function('f4X')(phif, Xe)
tf4f = sp.Function('tf4')(phif)
tf3f = sp.Function('tf3')(phif, Xe)            # vanishes structurally (checked)

# new canonical repo symbols introduced by the direct variation
f2phis = sp.Symbol('f2phi')
f3phis = sp.Symbol('f3phi')
f4phis = sp.Symbol('f4phi')
f4Xphis = sp.Symbol('f4Xphi')   # = repo column 'f4phiX' (d^2 f4 / dphi dX)
tf4phis = sp.Symbol('tf4phi')
phpps = sp.Symbol('phpp')       # phi''
ap2s = sp.Symbol('ap2')         # A0''

_GENUINE_SVT_JETS = (f3, f3X, f4, f4X, f4XX, tf4,
                     f3phis, f4phis, f4Xphis, tf4phis)


# ------------------------------------------------------------------ reduced Lagrangians
def sqrt_g():
    """sqrt(-g)/sin(theta) = sqrt(f/h) * C (angular 4pi convention stripped)."""
    return sp.sqrt(ff/hf)*CF


def ricci_scalar_C():
    """Ricci scalar of ds^2 = -f dt^2 + dr^2/h + C dOmega^2 (machine-derived,
    Bianchi- and Schwarzschild/de Sitter/flat-verified; see
    tools/derive_svt_reduced_action.py)."""
    Cp, Cpp = sp.diff(CF, r), sp.diff(CF, r, 2)
    fp, fpp = sp.diff(ff, r), sp.diff(ff, r, 2)
    hp = sp.diff(hf, r)
    return sp.cancel(-(
        2*CF**2*ff*hf*fpp + CF**2*ff*fp*hp - CF**2*hf*fp**2
        + 4*CF*ff**2*hf*Cpp + 2*CF*ff**2*Cp*hp - 4*CF*ff**2
        + 2*CF*ff*hf*Cp*fp - ff**2*hf*Cp**2
    )/(2*CF**2*ff**2))


def r_thetaphi_C():
    """4 Ftilde^2 R_{theta phi theta phi} angular-free factor: (4C - h C'^2)/4
    is R_{theta phi theta phi}/sin^2(theta) on this metric (C=r^2 gives the
    textbook r^2(1-h))."""
    return (4*CF - hf*sp.diff(CF, r)**2)/4


def reduced_lagrangians():
    """The five reduced pieces (angles stripped, M_pl^2 = 1) on the C(r)-metric."""
    A0p = sp.diff(A0f, r)
    php = sp.diff(phif, r)
    L_EH = sp.cancel(sqrt_g()*ricci_scalar_C()/2)
    L_f2 = sp.cancel(sqrt_g()*f2eff_f)
    # verbatim covariant-tool outputs (tools/derive_svt_reduced_action.py):
    # the C-dependence of the SVT pieces enters only via C', C'^2/C and 1/C -
    # do NOT multiply these by sqrt_g() (which carries a spurious factor C).
    L3 = sp.cancel(sp.sqrt(ff/hf)*f3f*hf**2*A0p**2*sp.diff(CF, r)*php/ff)
    L4a = sp.cancel(sp.sqrt(ff/hf)*(f4Xf/2 + tf4f)*hf**3*A0p**2
                    * sp.diff(CF, r)**2*php**2/(2*CF*ff))
    L4b = sp.cancel((4*CF - hf*sp.diff(CF, r)**2)*f4f*A0p**2
                    / (sp.sqrt(ff/hf)*CF))
    return {'L_EH': L_EH, 'L_f2': L_f2, 'L3': L3, 'L4a': L4a, 'L4b': L4b}


# ------------------------------------------------------------------ EL machinery
def euler_lagrange(expr, q, order):
    """E_q = dL/dq - d/dr dL/dq' + ... + (-1)^order d^order/dr^order dL/dq^(order)."""
    e = sp.diff(expr, q)
    dq = q
    for k in range(1, order + 1):
        dq = sp.diff(dq, r)
        term = sp.diff(expr, dq)
        for _ in range(k):
            term = sp.diff(term, r)
        e = e + (-1)**k*term
    return sp.expand(e)


# ------------------------------------------------- canonicalization to repo symbols
def canonicalize(expr):
    """doit() holonomic chains and map jet partials to canonical repo symbols."""
    e = sp.expand(expr.doit())
    SYM = {
        ('f2eff', 0): f2phis, ('f2eff', 1): f2X,
        ('f2eff', 2): f2F - 2*hf*ph**2*f2Y,          # repo f2F_eff convention
        ('f3', 0): f3phis, ('f3', 1): f3X,
        ('f4', 0): f4phis, ('f4', 1): f4X,
        ('f4X', 0): f4Xphis, ('f4X', 1): f4XX,
        ('tf4', 0): tf4phis,
        ('tf3', 0): sp.Symbol('tf3phi'), ('tf3', 1): sp.Symbol('tf3X'),
    }
    VAL = {'f2eff': f2, 'f3': f3, 'f4': f4, 'f4X': f4X, 'tf4': tf4,
           'tf3': sp.Symbol('tf3')}
    subs_map = {}
    for s in e.atoms(sp.Subs):
        inner = s.expr
        name = inner.expr.func.__name__
        var = s.variables[0]
        idx = None
        for i, a in enumerate(inner.expr.args):
            if a == var:
                idx = i
                break
        if idx is None:
            raise AssertionError(f'unresolvable Subs atom {s}')
        subs_map[s] = SYM[(name, idx)]
    e = e.subs(subs_map)
    for name, val in VAL.items():
        e = e.subs(sp.Function(name)(phif, Xe, Fe), val)
        e = e.subs(sp.Function(name)(phif, Xe), val)
        e = e.subs(sp.Function(name)(phif), val)
    # leftover partials on canonicalized symbols
    for part_sym, base in ((tf4phis, tf4), (f4phis, f4), (f4Xphis, f4X),
                           (f3phis, f3), (f2phis, f2)):
        e = e.subs(sp.Derivative(base, phif), part_sym)
    return sp.expand(e)


def to_repo_symbols(expr, *, areal_gauge=True):
    """Canonicalize; substitute field derivatives; C -> r^2 LAST (areal gauge)."""
    e = canonicalize(expr)
    if areal_gauge:
        e = e.subs(sp.Derivative(CF, (r, 4)), 0)
        e = e.subs(sp.Derivative(CF, (r, 3)), 0)
        e = e.subs(sp.Derivative(CF, (r, 2)), 2)
        e = e.subs(sp.Derivative(CF, r), 2*r)
        e = e.subs(CF, r**2)
    e = e.subs(sp.Derivative(ff, (r, 2)), fpp_r)
    e = e.subs(sp.Derivative(ff, r), fp_r)
    e = e.subs(ff, f)
    e = e.subs(sp.Derivative(hf, (r, 2)), sp.Symbol('hp2'))
    e = e.subs(sp.Derivative(hf, r), hp_r)
    e = e.subs(hf, h)
    e = e.subs(sp.Derivative(phif, (r, 2)), phpps)
    e = e.subs(sp.Derivative(phif, r), ph)
    e = e.subs(sp.Derivative(A0f, (r, 2)), ap2s)
    e = e.subs(sp.Derivative(A0f, r), ap)
    return sp.expand(e)


# ------------------------------------------------------------------ Delta22 core
def delta22_channels():
    """Raw EL channels of the C(r)-metric reduced action (C still symbolic)."""
    L = reduced_lagrangians()
    L_MH = sp.expand(L['L_EH'] + L['L_f2'])
    L_SVT = sp.expand(L['L3'] + L['L4a'] + L['L4b'])
    L_full = sp.expand(L_MH + L_SVT)
    return {
        'E_f': euler_lagrange(L_full, ff, 2),
        'E_h': euler_lagrange(L_full, hf, 2),
        'E_C_MH': euler_lagrange(L_MH, CF, 2),
        'E_C_SVT': euler_lagrange(L_SVT, CF, 2),
    }


def build_delta22_svt(*, verify=True):
    """Delta22_SVT: direct genuine-SVT theta-theta variation, repo E22 norm.

    Definition (authoritative):
        Delta22_SVT = E22_full - P_MH[E22_full]
                    = kappa_C * EL_C[L_SVT] |_{C=r^2},  kappa_C = sqrt(h/f).
    Areal-gauge guard: EL operators are formed with C(r) SYMBOLIC."""
    ch = delta22_channels()
    E_C_SVT_repo = to_repo_symbols(ch['E_C_SVT'], areal_gauge=True)
    assert not E_C_SVT_repo.has(CF), 'areal gauge leaked into the variation'
    delta22 = sp.expand(sp.sqrt(h)/sp.sqrt(f)*E_C_SVT_repo)
    if verify:
        rep = verify_normalization_channels()
        assert rep['V1_residual'] == 0, 'V1 (E00_full channel) failed'
        assert rep['V2_residual'] == 0 and rep['V2_kappa_metric_only'], 'V2 failed'
        assert rep['V3_residual'] == 0 and rep['V3_kappa_metric_only'], 'V3 failed'
        assert p_mh_delta22_zero(delta22), 'P_MH[Delta22_SVT] != 0'
        assert second_order_field_check(delta22), 'third-order channels present'
    return delta22


# ------------------------------------------------------------------ verification battery
def p_mh_delta22_zero(delta22):
    """P_MH[Delta22_SVT] == 0: killing all genuine-SVT operator families
    (values AND phi-channel partials) leaves exactly zero."""
    out = delta22
    for j in _GENUINE_SVT_JETS:
        out = out.subs(j, 0)
    return sp.simplify(out) == 0


def second_order_field_check(expr):
    """Second-order-EOM property: no h''/phi'''/A0''' channels may appear."""
    return not (expr.has(sp.Symbol('hp2')) or expr.has(sp.Symbol('phppp'))
                or expr.has(sp.Symbol('appp')))


def verify_normalization_channels():
    """Exact channel checks against the frozen repository transcriptions.

    V1: kappa_f * EL_f == E00_full     (kappa_f = -2 f^(3/2) sqrt(h))
    V2: kappa_h * EL_h == E11_full     (kappa_h = +2 sqrt(f) h^(3/2))
    V3: kappa_C * EL_C[L_MH] == E22_MH_slice  (kappa_C = sqrt(h/f))
    All kappa are metric-only; residuals must be EXACTLY zero."""
    from . import light_ring_identity as lri
    ch = delta22_channels()
    E_f = to_repo_symbols(ch['E_f'])
    E_h = to_repo_symbols(ch['E_h'])
    E_C_MH = to_repo_symbols(ch['E_C_MH'])
    kappa_f = -2*f**sp.Rational(3, 2)*sp.sqrt(h)
    kappa_h = 2*sp.sqrt(f)*h**sp.Rational(3, 2)
    kappa_C = sp.sqrt(h)/sp.sqrt(f)
    return {
        'V1_residual': sp.simplify(kappa_f*E_f - lri.E00_full()),
        'V2_residual': sp.simplify(kappa_h*E_h - lri.E11_full()),
        'V3_residual': sp.simplify(kappa_C*E_C_MH - lri.E22_MH_slice()),
        'V2_kappa_metric_only': _metric_only(kappa_h),
        'V3_kappa_metric_only': _metric_only(kappa_C),
    }


def _metric_only(expr):
    """True iff expr carries no jet / electric / scalar-field content."""
    jet_syms = (f2, f2X, f2F, f2Y, f2phis, f3, f3X, f3phis, f4, f4X,
                f4XX, f4phis, f4Xphis, tf4, tf4phis, ap2s, phpps)
    return not any(expr.has(s) for s in jet_syms) and not expr.has(ap)
