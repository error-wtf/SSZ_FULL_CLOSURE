"""KT2023 Maxwell-Horndeski background core (primary-source transcription).

Sources (verbatim transcriptions):
- Kase & Tsujikawa, arXiv:2301.10362 (PRD 107, 104045 (2023)):
    Eq. (7)  E00,  Eq. (8) E11,  Eq. (9) E22,  Eq. (10) C1..C10,
    Eq. (14) EA0,  Eq. (16) scalar-equation bracket (ON-SHELL rewriting),
    Eq. (24) Hcal, Eq. (25) Fcal, Eq. (26) Gcal,
    Appendix A Eq. (168): a4 = sqrt(fh)/2 * Hcal, v8 = G2F/(2 sqrt(fh)),
    Eq. (85) a4' relation (undivided form in light_ring_identity.C85_undivided).
- Heisenberg & Tsujikawa, arXiv:1802.07035 (PLB 780, 638 (2018)):
    Eq. (12) reduced action, Eqs. (14)/(15)/(18)/(19)/(20) f-, h-variations,
    scalar/vector currents (repo production evaluator conventions).

Established executable facts (regression-anchored here):

1. S0-SLICE DIAGONAL MAP.  On the common subtheory S0 (Einstein + f2(X,F):
   G3=G5=0, G4=Mpl^2 const, genuine SVT jets f3=f3X=tf4=0), the KT2023
   component normalizations and the HT2018 reduced-action variations are
   related by the DIAGONAL map

       KT_E00 = -HT_E00 / (r^2 f)      KT_E11 = +HT_E11 / (r^2 f)

   (verified exactly in test_s0_slice_diagonal_map).

2. EQ.-85-GENERATING COMBINATION (C_bg recipe).  In the project-locked
   luminal Horndeski branch (G4 = G4(phi), G3 = 0, G5 = 0) the undivided
   Eq. (85) expression is an EXACT OFF-SHELL two-term identity,

       C85_undivided == A * KT_E00 + B * KT_E22,
       A = -r f^(3/2)/sqrt(h),   B = +r f^(3/2)/sqrt(h),

   with A, B metric-only (verified exactly in
   test_eq85_generating_two_term_identity).  This is the C_bg-generating
   combination: the general U(1)-SVT C_bg is the SAME combination applied
   to the general theory's (E00, E22); its MH slice is C85_undivided and
   its genuine-SVT remainder is the Sigma_SVT candidate content.

3. SECTOR SCOPE.  Outside the luminal branch (e.g. general G4(phi,X) with
   G4X != 0) the two-term metric-coefficient identity closes no longer;
   Eq. (85) there is an on-shell statement involving additional input.
   The project action family is locked to the luminal branch, so this is
   outside current scope and is recorded as probe evidence.

Premises (Section 32 policy): all denominators appearing (f, sqrt(f*h), r)
are nonzero on the physical branch f>0, h>0, r>0.
"""
from __future__ import annotations

import sympy as sp

# ---------------------------------------------------------------- symbols
r, f, h, ph, app = sp.symbols('r f h ph app', positive=True)
fp, fpp, hp = sp.symbols('fp fpp hp', positive=True)
ppp = sp.Symbol('phpp')          # phi''  (background scalar second derivative)

# Maxwell-Horndeski jets (functions of r; transcription treats them as symbols)
G2, G2X, G2F = sp.symbols('G2 G2X G2F')
G3, G3X, G3phi = sp.symbols('G3 G3X G3phi')
G4, G4X, G4phi, G4phiX, G4XX = sp.symbols('G4 G4X G4phi G4phiX G4XX')
G5, G5X, G5phi, G5phiX, G5XX = sp.symbols('G5 G5X G5phi G5phiX G5XX')
G5phiphi = sp.Symbol('G5phiphi')   # d^2 G5 / dphi^2 (needed in Hcal_phi)
JETS = [G2, G2X, G2F, G3, G3X, G3phi, G4, G4X, G4phi, G4phiX, G4XX,
        G5, G5X, G5phi, G5phiX, G5XX, G5phiphi]
# r-derivatives of the jets (chain-rule bookkeeping)
dG = {s: sp.Symbol('d' + s.name) for s in JETS}

hpp, app2, fppp = sp.symbols('hpp app2 fppp')
CHAIN = {ph: ppp, fp: fpp, hp: hpp, app: app2, fpp: fppp}


def dr(expr):
    """Total radial derivative with full jet chain rule."""
    e = sp.diff(expr, r)
    for s in JETS:
        e = e + sp.diff(expr, s) * dG[s]
    for s, ds in CHAIN.items():
        e = e + sp.diff(expr, s) * ds
    return e


# ------------------------------------------------- KT2023 Eq. (10) C-jets
def c_coefficients():
    C1 = -h**2 * (G3X - 2 * G4phiX) * ph**2 - 2 * G4phi * h
    C2 = 2 * h**3 * (2 * G4XX - G5phiX) * ph**3 - 4 * h**2 * (G4X - G5phi) * ph
    C3 = -h**4 * G5XX * ph**4 + h**2 * G5X * (3 * h - 1) * ph**2
    C4 = (h**2 * (2 * G4XX - G5phiX) * ph**4
          + h * (3 * G5phi - 4 * G4X) * ph**2 - 2 * G4)
    C5 = -sp.Rational(1, 2) * (G5XX * h**3 * ph**5 - h * G5X * (5 * h - 1) * ph**3)
    C6 = h * (G3phi - 2 * G4phi) * ph**2 + G2
    C7 = -2 * h**2 * (2 * G4phiX - G5phi) * ph**3 - 4 * G4phi * h * ph
    C8 = (G5phiX * h**3 * ph**4 - h * (2 * G4X * h - G5phi * h - G5phi) * ph**2
          - 2 * G4 * (h - 1))
    C9 = -h * (G2X - G3phi) * ph**2 - G2
    C10 = (sp.Rational(1, 2) * G5phiX * h**3 * ph**4
           - sp.Rational(1, 2) * h**2 * (2 * G4X - G5phi) * ph**2 - G4 * h)
    return dict(C1=C1, C2=C2, C3=C3, C4=C4, C5=C5,
                C6=C6, C7=C7, C8=C8, C9=C9, C10=C10)


# ------------------------------------------ KT2023 Eqs. (7)/(8)/(9)/(14)
def kt_e00():
    C = c_coefficients()
    return ((C['C1'] + C['C2'] / r + C['C3'] / r**2) * ppp
            + (ph / (2 * h) * C['C1'] + C['C4'] / r + C['C5'] / r**2) * hp
            + C['C6'] + C['C7'] / r + C['C8'] / r**2
            - h / f * G2F * app**2)


def kt_e11():
    C = c_coefficients()
    return (-(ph / (2 * h) * C['C1'] + C['C4'] / r + C['C5'] / r**2) * h * fp / f
            + C['C9'] - 2 * ph / r * C['C1']
            - (ph / (2 * h) * C['C2'] + (h - 1) * C['C4']) / r**2
            + h / f * G2F * app**2)


def kt_e22():
    C = c_coefficients()
    return (( (C['C2'] + ((2 * h - 1) * ph * C['C3'] + 2 * h * C['C5'])
              / (h * ph * r)) * fp / (4 * f) + C['C1'] + C['C2'] / (2 * r) ) * ppp
            + (2 * h * C['C4'] - ph * C['C2']
               + (2 * h * C['C5'] - ph * C['C3']) / r)
            / (4 * f) * (fpp - fp**2 / (2 * f))
            + (C['C4'] + (2 * h * (2 * h + 1) * C['C5'] - ph * C['C3'])
               / (2 * h**2 * r)) * fp * hp / (4 * f)
            + (C['C7'] / 4 + C['C10'] / r) * fp / f
            + (ph / h * C['C1'] + C['C4'] / r) * hp / 2
            + C['C6'] + C['C7'] / (2 * r))


def kt_ea0():
    """Eq. (14): (G2F sqrt(h/f) r^2 A0')' = 0 (undifferentiated current below)."""
    return sp.sqrt(h / f) * r**2 * app * G2F


# ----------------------------- KT2023 Eqs. (24)/(25)/(168) tensor slots
def hcal():
    return 2 * G4 + 2 * h * ph**2 * G4X - h * ph**2 * G5phi - h**2 * ph**3 * G5X / r


def fcal():
    return (2 * G4 + h * ph**2 * G5phi
            - h * ph**2 * (hp * ph / 2 + h * ppp) * G5X)


def gcal():
    return (2 * G4 + 2 * h * ph**2 * G4X
            - h * ph**2 * (G5phi + fp * h * ph * G5X / (2 * f)))


def a4_kt():
    """Appendix A Eq. (168): a4 = sqrt(fh)/2 * Hcal."""
    return sp.sqrt(f * h) / 2 * hcal()


def v8_kt():
    """Appendix A Eq. (168): v8 = G2F/(2 sqrt(fh))  (canonical MH v8 slot)."""
    return G2F / (2 * sp.sqrt(f * h))


def a4_kt_prime():
    """Chain-rule a4' (jets symbolic).  The partial derivatives of Hcal use the
    PARTIAL-derivative jet symbols (G4phi = dG4/dphi, G4phiX = d^2G4/dphi dX,
    ...), NOT radial-derivative symbols; the only new radial jet entering is
    via X' = -(h'phi'^2 + 2 h phi' phi'')/2 and the metric factors."""
    H = hcal()
    H_phi = 2 * G4phi + 2 * h * ph**2 * G4phiX - h * ph**2 * G5phiphi \
        - h**2 * ph**3 * G5phiX / r
    H_X = 2 * h * ph**2 * G4XX - h * ph**2 * G5phiX - h**2 * ph**3 * G5XX / r
    H_h = 2 * ph**2 * G4X - ph**2 * G5phi - 2 * h * ph**3 * G5X / r
    H_r = h**2 * ph**3 * G5X / r**2
    Xp = -(hp * ph**2 + 2 * h * ph * ppp) / 2     # X = -h ph^2 / 2
    a4p_metric = sp.sqrt(f * h) / 2 * (fp / f + hp / h) / 2
    return a4p_metric * H + sp.sqrt(f * h) / 2 * (H_phi * ph + H_X * Xp + H_h * hp + H_r)


def c85_undivided_kt():
    """Undivided Eq. (85): (2f-rf')a4' - [G_r a4 + f^(3/2)Fcal/(r sqrt(h)) - 2rfh ap^2 v8]."""
    G_r = r * fpp - r * fp**2 / f + 2 * fp - 2 * f / r
    return ((2 * f - r * fp) * a4_kt_prime() - G_r * a4_kt()
            - f**sp.Rational(3, 2) / (r * sp.sqrt(h)) * fcal()
            + 2 * r * f * h * app**2 * v8_kt())


# ------------------------------- S0 slice and the diagonal HT map check
def s0_substitutions(mpl_sq=1):
    """Common subtheory S0: Einstein + f2(X,F).  G3=G5=0, G4=mpl_sq/2 const,
    all genuine-SVT jets zero.  Returns the jet-substitution dict."""
    subs = {}
    for s in (G3, G3X, G3phi, G4X, G4phi, G4phiX, G4XX,
              G5, G5X, G5phi, G5phiX, G5XX):
        subs[s] = 0
    subs[G4] = mpl_sq / 2
    for s in JETS:
        subs[dG[s]] = 0
    return subs


def ht_e00_s0(mpl_sq=1):
    """HT2018 Eq. (14) restricted to S0 (f-variation, Mpl^2 = mpl_sq)."""
    f2, f2F = sp.symbols('f2 f2F')
    return (mpl_sq * r * f * hp - mpl_sq * f * (1 - h)
            - r**2 * (f * f2 - h * app**2 * f2F))


def ht_e11_s0(mpl_sq=1):
    """HT2018 Eq. (15) restricted to S0 (h-variation)."""
    f2, f2X, f2F = sp.symbols('f2 f2X f2F')
    return (mpl_sq * r * h * fp - mpl_sq * f * (1 - h)
            - r**2 * (f * f2 + f * h * ph**2 * f2X - h * app**2 * f2F))


def s0_slice_diagonal_map(mpl_sq=1):
    """Verify KT_E00 = -HT_E00/(r^2 f) and KT_E11 = +HT_E11/(r^2 f) on S0,
    with the jet identification G2 -> f2, G2X -> f2X, G2F -> f2F."""
    f2, f2X, f2F = sp.symbols('f2 f2X f2F')
    ident = {G2: f2, G2X: f2X, G2F: f2F}
    e00_kt = kt_e00().subs(s0_substitutions(mpl_sq)).subs(ident)
    e11_kt = kt_e11().subs(s0_substitutions(mpl_sq)).subs(ident)
    r00 = sp.simplify(e00_kt + ht_e00_s0(mpl_sq) / (r**2 * f))
    r11 = sp.simplify(e11_kt - ht_e11_s0(mpl_sq) / (r**2 * f))
    return {'residual_E00_map': sp.simplify(r00),
            'residual_E11_map': sp.simplify(r11),
            'map_holds': r00 == 0 and r11 == 0}


def probe_offshell_two_term_residual(sector=None):
    """C_bg-generating combination evidence.  Returns the metric-only
    coefficients A, B fixed by the app^2-channel and the status of the
    exact off-shell two-term identity C85 == A*E00 + B*E22 in the given
    jet sector.  Luminal branch: holds EXACTLY (C_bg recipe).  Outside the
    luminal branch (G4X != 0) it closes no longer (sector-scope evidence)."""
    subs = dict(sector if sector is not None else LUMINAL_G4PHI_SECTOR)
    E00 = sp.expand(kt_e00().subs(subs))
    E22 = sp.expand(kt_e22().subs(subs))
    C85 = sp.expand(c85_undivided_kt().subs(subs))
    A_app = sp.simplify(sp.cancel(sp.expand(C85).coeff(app**2)
                                  / sp.expand(E00).coeff(app**2)))
    # G2 channel (E00 and E22 both carry C6 = G2 with unit weight):
    B_G2 = -A_app
    full_residual = sp.cancel(sp.together(C85 - A_app * E00 - B_G2 * E22))
    return {'A': sp.simplify(A_app),
            'B': sp.simplify(B_G2),
            'two_term_identity_exact': full_residual == 0,
            'full_residual': sp.factor(full_residual)}


# ------------------------------------------------ Eq. (85) ON-SHELL CERT
LUMINAL_G4PHI_SECTOR = {
    # project-locked Horndeski content: G4 = G4(phi), G3 = 0, G5 = 0
    G3: 0, G3X: 0, G3phi: 0,
    G4X: 0, G4phiX: 0, G4XX: 0,
    G5: 0, G5X: 0, G5phi: 0, G5phiX: 0, G5XX: 0, G5phiphi: 0,
}


def onshell_eq85_certificate(sector=None):
    """Machine certificate of the undivided Eq. (85) in a given jet sector.

    Procedure (primary-source-faithful, no fitting):
      1. transcribe E00 (Eq. 7), E22 (Eq. 9) and the undivided Eq. (85)
         expression (Eqs. 24/25/168 chain rule) with all jets symbolic;
      2. solve E00 = 0, E22 = 0 for the background second derivatives
         (phi'', h') — the KT2023 on-shell usage;
      3. substitute and require the undivided Eq. (85) expression == 0.

    Returns the pivot determinant (premises), the on-shell phi''/h'
    solutions and the residual (0 = certificate holds).  Premises: every
    nonzero pivot factor, i.e. f>0, r>0, G4!=0, G4phi!=0, h!=0 and
    2f - r f' != 0 (Eq. 85 divides by the light-ring denominator; the
    undivided form carries the ring limit — light_ring_identity policy).
    """
    subs = dict(sector or LUMINAL_G4PHI_SECTOR)
    E00 = sp.expand(kt_e00().subs(subs))
    E22 = sp.expand(kt_e22().subs(subs))
    C85 = sp.expand(c85_undivided_kt().subs(subs))
    K00, L00 = E00.coeff(ppp), E00.coeff(hp)
    K22, L22 = E22.coeff(ppp), E22.coeff(hp)
    M00 = sp.expand(E00.subs({ppp: 0, hp: 0}))
    M22 = sp.expand(E22.subs({ppp: 0, hp: 0}))
    D = sp.cancel(sp.together(K00 * L22 - K22 * L00))
    ppp_os = sp.cancel(sp.together((-M00 * L22 + M22 * L00) / D))
    hp_os = sp.cancel(sp.together((-K00 * M22 + K22 * M00) / D))
    residual = sp.simplify(
        sp.cancel(sp.together(C85.subs({ppp: ppp_os, hp: hp_os}))))
    return {'sector': {str(k): v for k, v in subs.items()},
            'pivot_D': sp.factor(D),
            'phi2_onshell': sp.factor(ppp_os),
            'hp_onshell': sp.factor(hp_os),
            'undivided_eq85_onshell_residual': residual,
            'certificate_holds': residual == 0}


# ------------------------------------------------ Ward-span no-go probe
def probe_ward_span_no_go():
    """SUPERSEDED REGISTER (do not use for gating).

    The original version claimed the theta-theta equation cannot be
    reconstructed as E22 = Ward + sum c_i*EOM_i with metric coefficients.
    That claim was an ARTIFACT of a non-holonomic jet space (it treated
    d(f2)/dr as an independent chain direction).  With the holonomic
    chain rule (d f2/dr = f2phi*phi' + f2X*X' + f2F*F' + f2Y*Y', exactly
    the repo convention in svt_background_eom.py:219) the df2 obstruction
    vanishes: in the luminal MH slice all derivative channels (df2F,
    df2X, phpp, hpp, app2, fpp) close with METRIC-ONLY coefficients

        c2 = ap*sqrt(h)*(f r^2 - 1)/(2 r sqrt(f))      (JA' channel)
        c3 = sqrt(f)*h*phi_r*r/(2 sqrt(h))              (scalar-EOM channel)
        c0 = r*(f hp + fp h)/(4 f h)
        c1 = -(-4 f h + f hp r + fp h r)/(4 f h)

    and the E22 normalization map is measured by the fpp channel as
    E22_HT = r^2 f * E22_KT.  ONE unmatched free residual remains
    (see data/generated/phase2_q2/STEP1_DELTA22_WARD_SPAN_SOLUTION.json);
    matching it completes Delta22_SVT.  This probe is kept as the honest
    record of the falsified shortcut claim."""
    return {
        "status": "SUPERSEDED_NON_HOLONOMIC_ARTIFACT",
        "df2_coefficients": {"Ward": "f*r**3/2 (non-holonomic space only)",
                             "holonomic": "absorbed by df2 = f2phi*ph + f2X*Xp + f2F*Fp + f2Y*Yp"},
        "correct_route": "Ward + c0*E00 + c1*E11 + c2*JA' + c3*(Jphi'-Pphi), "
                         "holonomic chains, E22_HT = r^2 f * E22_KT",
        "evidence": "data/generated/phase2_q2/STEP1_DELTA22_WARD_SPAN_SOLUTION.json",
    }
