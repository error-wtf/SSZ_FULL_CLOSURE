"""G07: independent holonomic action-jet re-check (multiple reconstructions).

The holonomic closure of the action jets is validated through THREE mutually
independent derivative reconstructions:

  A. sympy multivariate-function chains (ssz_p5.action.svt_direct_variation
     pipeline; Xe/Fe are genuine Function arguments, d/dr differentiates
     through compositions).

  B. symbol-level chain tables (jets as flat symbols; total derivatives
     applied via the explicit holonomic replacement
     d f3/dr = f3phi phi' + f3X X',  d f4X/dr = f4phiX phi' + f4XX X', ...).

  C. inverse-metric variation convention: varying w = g^thetatheta = 1/C
     instead of C must act on the algebraic sector exactly as
     delta g^thetatheta = -(g^thetatheta)^2 delta g_thetatheta predicts,
     i.e.  EL_w[L |_{C=1/w}] == -(1/w^2) * EL_C[L] |_{C=1/w}.

Reconstructions A and B must produce the SAME Delta22_SVT expression
exactly; any drift between chain implementations breaks G07.
"""
import sympy as sp

from ssz_p5.action import light_ring_identity as lri
from ssz_p5.action import svt_direct_variation as sdv


def _symbol_chain_L_SVT():
    """Reconstruction B: L_SVT3 + L_SVT4a + L_SVT4b with jets as flat repo
    symbols and holonomic chains applied by explicit replacement tables."""
    r, f, h, ph, ap = lri.r, lri.f, lri.h, lri.ph, lri.ap
    fp_r, hp_r = lri.fp_r, lri.hp_r
    f3, f3X, f4, f4X, f4XX, tf4 = lri.f3, lri.f3X, lri.f4, lri.f4X, lri.f4XX, lri.tf4
    f3phi = sp.Symbol('f3phi')
    f4phi = sp.Symbol('f4phi')
    f4Xphi = sp.Symbol('f4Xphi')
    tf4phi = sp.Symbol('tf4phi')
    phpp = sp.Symbol('phpp')
    ap2 = sp.Symbol('ap2')

    C = sp.Function('C')(r)
    Cp = sp.diff(C, r)
    Xp = -(hp_r*ph**2 + 2*h*ph*phpp)/2
    # f2-sector absent from L_SVT; SVT jet chains (holonomic, X = -h ph^2/2):
    chain = {f3: f3phi*ph + f3X*Xp,
             f4: f4phi*ph + f4X*Xp,
             f4X: f4Xphi*ph + f4XX*Xp,
             tf4: tf4phi*ph,
             ap: ap2,
             ph: phpp,
             h: hp_r,
             f: fp_r,
             phpp: 0,
             ap2: 0,
             hp_r: sp.Symbol('hp2'),
             fp_r: 0}

    def d_dr(expr):
        """Manual total derivative: flat symbols carry NO automatic chains -
        the holonomic replacement table IS the reconstruction under test."""
        e = sp.diff(expr, r)
        for s, v in chain.items():
            e = e + sp.diff(expr, s)*v
        return sp.expand(e)

    L3 = sp.sqrt(f/h)*f3*h**2*ap**2*Cp*ph/f
    L4a = sp.sqrt(f/h)*(f4X/2 + tf4)*h**3*ap**2*Cp**2*ph**2/(2*C*f)
    L4b = (4*C - h*Cp**2)*f4*ap**2/(sp.sqrt(f/h)*C)
    L_SVT = sp.expand(L3 + L4a + L4b)

    # EL w.r.t. C (L_SVT has no C'' channel): E = dL/dC - d_dr(dL/dC')
    assert sp.diff(L_SVT, sp.diff(C, r, 2)) == 0
    e = sp.diff(L_SVT, C) - d_dr(sp.diff(L_SVT, Cp))
    e = e.subs(sp.Derivative(C, (r, 3)), 0)
    e = e.subs(sp.Derivative(C, (r, 2)), 2)
    e = e.subs(sp.Derivative(C, r), 2*r)
    e = e.subs(C, r**2)
    e = e.subs(f3, lri.f3).subs(f4, lri.f4).subs(f4X, lri.f4X).subs(tf4, lri.tf4)
    e = e.subs(f3X, lri.f3X).subs(f4XX, lri.f4XX)
    return sp.expand(e*sp.sqrt(h)/sp.sqrt(f))


def test_reconstruction_B_matches_module_delta22():
    """Reconstruction B (flat-symbol chain tables) must reproduce the module's
    function-chain Delta22_SVT EXACTLY (zero residual, canonical service)."""
    from ssz_p5.action.symbolic_canon import canonical_simplify
    d22 = sdv.build_delta22_svt(verify=False)
    d22_B = _symbol_chain_L_SVT()
    assert canonical_simplify(d22 - d22_B) == 0


def test_holonomic_chain_identity_f2eff():
    """The sympy-function chain of the effective f2 jet equals the repo
    holonomic convention df2/dr = f2phi phi' + f2X X' + f2F_eff F' with
    f2F_eff = f2F - 2 h phi'^2 f2Y (svt_background_eom.py chain rule).
    Side A: function-composition chain via the module canonicalization.
    Side B: explicit primitive-derivative table (independent construction)."""
    # side A: function-composition chain, canonicalized to repo symbols
    df2_mapped = sdv.to_repo_symbols(sp.diff(sdv.f2eff_f, sdv.r).doit(),
                                     areal_gauge=False)
    # side B: explicit primitive-derivative table (independent construction)
    f2F_eff = lri.f2F - 2*lri.h*lri.ph**2*lri.f2Y
    Fp = (lri.f*(lri.hp_r*lri.ap**2 + 2*lri.h*lri.ap*sdv.ap2s)
          - lri.fp_r*lri.h*lri.ap**2)/(2*lri.f**2)
    Xp = -(lri.hp_r*lri.ph**2 + 2*lri.h*lri.ph*sdv.phpps)/2
    df2_table = sp.expand(sdv.f2phis*lri.ph + lri.f2X*Xp + f2F_eff*Fp)
    assert sp.simplify(df2_mapped - df2_table) == 0


def test_inverse_metric_variation_convention():
    """Convention pin (prompt Sec. 7): on the ALGEBRAIC sector (L depends on C
    but not C', C''), varying w = g^thetatheta = 1/C acts as
    delta g^thetatheta = -(g^thetatheta)^2 delta g_thetatheta:
        EL_w[L|_{C=1/w}] == -(1/w^2) * EL_C[L]|_{C=1/w}.
    The algebraic C-channel is exactly the f2-sector channel (X, F, Y are
    C-independent on the static electric background)."""
    r = sdv.r
    w = sp.Function('_w_inv')(r)
    L_f2 = sdv.reduced_lagrangians()['L_f2']
    # C-channel (algebraic: no C', C'' dependence)
    assert not sp.diff(L_f2, sp.diff(sdv.CF, r)).has(sdv.CF)
    EL_C = sp.diff(L_f2, sdv.CF)
    # w-channel: substitute C = 1/w, differentiate w.r.t. w
    L_w = L_f2.subs(sdv.CF, 1/w)
    EL_w = sp.diff(L_w, w)
    pred = -1/w**2 * EL_C.subs(sdv.CF, 1/w)
    assert sp.simplify(sp.expand(EL_w - pred)) == 0


def test_no_nonholonomic_df2_channel():
    """Regression for the historical no-go artifact: d(f2)/dr must never be
    treated as an independent chain direction - the holonomic expansion
    absorbs it (df2 = f2phi phi' + f2X X' + f2F_eff F')."""
    mapped = sdv.to_repo_symbols(sp.diff(sdv.f2eff_f, sdv.r).doit(),
                                 areal_gauge=False)
    # after the holonomic mapping no unevaluated jet derivative remains
    assert not mapped.has(sp.Subs)
    assert not mapped.has(sp.Derivative(sdv.phif, sdv.r))
    assert not mapped.has(sp.Derivative(sdv.A0f, sdv.r))
    # and it lives exclusively on the documented channel symbols
    allowed = (sdv.f2phis, lri.f2X, lri.f2F, lri.f2Y, lri.ph, lri.ap,
               lri.hp_r, lri.h, lri.f, lri.fp_r, sdv.phpps, sdv.ap2s)
    assert mapped.free_symbols <= set(allowed), \
        mapped.free_symbols - set(allowed)
