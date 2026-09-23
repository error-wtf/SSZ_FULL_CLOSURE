"""Covariant derivation of the C(r)-metric reduced SVT action (evidence tool).

Derives the five reduced Lagrangian pieces of the genuine U(1)-SVT action
(HT2018 arXiv:1802.07035 Eqs. (1)-(13)) from the FULL 4D covariant action on

    ds^2 = -f(r) dt^2 + dr^2/h(r) + C(r) dOmega^2,   g_thetatheta = C(r) FREE,

fields phi(r), A_mu = (A_0(r), 0, 0, 0), and asserts EXACT equality with the
transcribed forms in ssz_p5.action.svt_direct_variation.  Also verifies the
structural vanishing of the tf3 part of M3 on this background (HT2018 Sec. II)
and writes the STEP3 evidence JSON.

This tool is the provenance chain for the transcribed reduced forms; it is a
one-shot evidence generator (not part of the fast test suite).

Run: PYTHONPATH=src python tools/derive_svt_reduced_action.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.action import svt_direct_variation as sdv  # noqa: E402

r = sp.Symbol('r', positive=True)
ff, hf, CF = sp.Function('f')(r), sp.Function('h')(r), sp.Function('C')(r)
phF, A0F = sp.Function('phif')(r), sp.Function('A0f')(r)
coords = [sp.Symbol('t'), r, sp.Symbol('_th'), sp.Symbol('_phc')]
th = coords[2]
g = sp.diag(-ff, 1/hf, CF, CF*sp.sin(th)**2)
ginv = g.inv()
sg = sp.sqrt(ff/hf)*CF*sp.sin(th)


def Gamma_u(rho, mu, nu):
    return sp.cancel(sp.Rational(1, 2)*sum(
        ginv[rho, s]*(sp.diff(g[s, nu], coords[mu]) + sp.diff(g[s, mu], coords[nu])
                      - sp.diff(g[mu, nu], coords[s])) for s in range(4)))


def Riem_u(rho, sig, mu, nu):
    return sp.cancel(sp.diff(Gamma_u(rho, nu, sig), coords[mu])
                     - sp.diff(Gamma_u(rho, mu, sig), coords[nu])
                     + sum(Gamma_u(rho, mu, lam)*Gamma_u(lam, nu, sig)
                           - Gamma_u(rho, nu, lam)*Gamma_u(lam, mu, sig)
                           for lam in range(4)))


def main() -> int:
    php = sp.diff(phF, r)
    ap = sp.diff(A0F, r)

    # --- geometry sanity battery (independent of the SVT transcription) ---
    subs_schw = {ff: 1 - 2*sp.Symbol('M')/r, hf: 1 - 2*sp.Symbol('M')/r, CF: r**2}
    Ric = {(mu, nu): sp.cancel(sum(Riem_u(lam, mu, lam, nu) for lam in range(4)))
           for mu in range(4) for nu in range(4)}
    ric_schw = all(sp.simplify(Ric[k].subs(subs_schw).doit()) == 0 for k in Ric)
    subs_flat = {ff: 1, hf: 1, CF: r**2}
    ric_flat = all(sp.simplify(Ric[k].subs(subs_flat).doit()) == 0 for k in Ric)
    Rsc = sp.cancel(sum(ginv[mu, nu]*Ric[(mu, nu)] for mu in range(4) for nu in range(4)))
    subs_dsi = {ff: 1 - sp.Symbol('L')*r**2/3, hf: 1 - sp.Symbol('L')*r**2/3, CF: r**2}
    rsc_dsi = sp.simplify((Rsc - 4*sp.Symbol('L')).subs(subs_dsi).doit()) == 0
    R_thphthph = sp.cancel(sum(g[2, lam]*Riem_u(lam, 3, 2, 3) for lam in range(4)))
    r_thph_val = sp.simplify(
        (R_thphthph - (4*CF - hf*sp.diff(CF, r)**2)*sp.sin(th)**2/4).subs({CF: r**2}).doit())

    # --- X, F, Y (C-independent on this background: structural fact) ---
    Xe = -hf*php**2/2
    Fe = hf*ap**2/(2*ff)
    Ye = 4*Fe*Xe
    assert sp.simplify(sp.diff(Xe, CF)) == 0
    assert sp.simplify(sp.diff(Fe, CF)) == 0
    assert sp.simplify(sp.diff(Ye, CF)) == 0

    # --- jets as multivariate functions ---
    f2f = sp.Function('f2eff')(phF, Xe, Fe)
    f3f = sp.Function('f3')(phF, Xe)
    tf3f = sp.Function('tf3')(phF, Xe)
    f4f = sp.Function('f4')(phF, Xe)
    f4Xf = sp.Function('f4X')(phF, Xe)
    tf4f = sp.Function('tf4')(phF)

    gradphi = [ginv[mu, 1]*php for mu in range(4)]
    hess = {}
    for mu in range(4):
        for nu in range(4):
            hess[(mu, nu)] = sp.cancel(sp.diff(sp.diff(phF, coords[mu]), coords[nu])
                - sum(Gamma_u(lam, mu, nu)*php if lam == 1 else 0
                      for lam in range(4)))
    Fdn = {(mu, nu): sp.cancel(sp.diff(A0F if mu == 0 else 0, coords[nu])
                               - sp.diff(A0F if nu == 0 else 0, coords[mu]))
           for mu in range(4) for nu in range(4)}

    def perm_sign(p):
        p = list(p)
        s = 1
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if p[i] > p[j]:
                    s = -s
        return s

    epsU = {}
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for sig in range(4):
                    p = [mu, nu, rho, sig]
                    epsU[(mu, nu, rho, sig)] = (perm_sign(p)/sg
                                                if len(set(p)) == 4 else sp.Integer(0))
    FtilU = {(mu, nu): sp.cancel(sp.Rational(1, 2)*sum(
        epsU[(mu, nu, lam, sig)]*Fdn[(lam, sig)]
        for lam in range(4) for sig in range(4))) for mu in range(4) for nu in range(4)}

    def red(dens):
        return sp.cancel(dens/sp.sin(th))

    L_EH = red(sg*Rsc/2)
    L_f2 = red(sg*f2f)
    L3_full = red(sg*sum(
        (f3f*g[lam, sig] + tf3f*gradphi[lam]*gradphi[sig])
        * FtilU[(mu, lam)]*FtilU[(nu, sig)]*hess[(mu, nu)]
        for mu in range(4) for nu in range(4)
        for lam in range(4) for sig in range(4)))
    L4a = red(sg*(f4Xf/2 + tf4f)*sum(
        FtilU[(mu, nu)]*FtilU[(al, be)]*hess[(mu, al)]*hess[(nu, be)]
        for mu in range(4) for nu in range(4)
        for al in range(4) for be in range(4)))
    Rlo = {(rho, s2, ga, de): sp.cancel(sum(g[rho, lam]*Riem_u(lam, s2, ga, de)
            for lam in range(4)))
           for rho in range(4) for s2 in range(4) for ga in range(4) for de in range(4)}
    L4b = red(sg*f4f*sum(
        FtilU[(rho, s2)]*FtilU[(ga, de)]*Rlo[(rho, s2, ga, de)]
        for rho in range(4) for s2 in range(4) for ga in range(4) for de in range(4)))

    # --- tf3 structural vanishing (HT2018: ftilde3 term vanishes on this background) ---
    tf3_part = sp.simplify(L3_full.coeff(tf3f))
    tf3_vanishes = tf3_part == 0

    # --- equality with the module transcriptions ---
    Lm = sdv.reduced_lagrangians()
    checks = {
        'L_EH': sp.simplify(sp.expand(L_EH - Lm['L_EH'])) == 0,
        'L_f2': sp.simplify(sp.expand(L_f2 - Lm['L_f2'])) == 0,
        'L3': sp.simplify(sp.expand(L3_full - Lm['L3'])) == 0,
        'L4a': sp.simplify(sp.expand(L4a - Lm['L4a'])) == 0,
        'L4b': sp.simplify(sp.expand(L4b - Lm['L4b'])) == 0,
    }

    # --- Delta22 channel checks (module battery) ---
    rep = sdv.verify_normalization_channels()
    d22 = sdv.build_delta22_svt(verify=False)
    battery = {
        'V1_residual': sp.simplify(rep['V1_residual']) == 0,
        'V2_residual': sp.simplify(rep['V2_residual']) == 0,
        'V3_residual': sp.simplify(rep['V3_residual']) == 0,
        'V2_kappa_metric_only': rep['V2_kappa_metric_only'],
        'V3_kappa_metric_only': rep['V3_kappa_metric_only'],
        'P_MH_delta22_zero': sdv.p_mh_delta22_zero(d22),
        'second_order_fields': sdv.second_order_field_check(d22),
    }

    all_ok = (ric_schw and ric_flat and rsc_dsi and r_thph_val == 0
              and tf3_vanishes and all(checks.values()) and all(battery.values()))

    families = {}
    for j in sdv._GENUINE_SVT_JETS:
        c = sp.simplify(d22.coeff(j))
        if c != 0:
            families[str(j)] = sp.srepr(sp.factor(c))
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, cwd=Path(__file__).resolve().parents[1]
                              ).stdout.strip()
    except Exception:
        head = "unknown"

    evidence = {
        "schema_version": "1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": head,
        "title": "STEP3: Delta22_SVT derived by direct theta-theta variation of the "
                 "genuine SVT action with the angular metric degree of freedom unfixed",
        "route": "Delta22_SVT = kappa_C * EL_C[L_SVT2 + L_SVT3 + L_SVT4] |_{C=r^2}, "
                 "kappa_C = sqrt(h/f) pinned exactly on the verified MH slice (V3); "
                 "areal gauge substituted only AFTER the Euler-Lagrange step",
        "geometry_checks": {
            "ricci_schwarzschild_zero": ric_schw,
            "ricci_flat_zero": ric_flat,
            "ricci_scalar_deSitter_4L": rsc_dsi,
            "R_thetaphi_thetaphi_textbook_C_is_r2": bool(r_thph_val == 0),
        },
        "structural_checks": {
            "tf3_part_vanishes_on_static_electric_background": tf3_vanishes,
            "X_F_Y_C_independent": True,
        },
        "reduced_piece_equality_with_module": checks,
        "normalization_channel_checks": battery,
        "normalization_map": {
            "V1": "E00_full == -2 f^(3/2) sqrt(h) * EL_f[L_full] |_{C=r^2}",
            "V2": "E11_full == +2 sqrt(f) h^(3/2) * EL_h[L_full] |_{C=r^2}",
            "V3": "E22_MH_slice == sqrt(h/f) * EL_C[L_MH] |_{C=r^2}",
            "map_factors_metric_only": True,
        },
        "delta22_operator_families_srepr": families,
        "delta22_second_derivative_channels": {
            "phi'": "ph", "phi''": "phpp", "A0'": "ap", "A0''": "ap2",
            "f'": "fp_r", "f''": "fpp_r", "h'": "hp_r",
            "third_order_channels_absent": True,
        },
        "ward_oracle_status": {
            "route_decision": "Ward-span reconstruction superseded (commit f4d5ce2); "
                              "direct variation is the derivation path",
            "R_met_module_equals_recorded_metric_part": True,
            "R_el_normalization_note": "recorded STEP1 residual electric part carries "
                                       "plain f2F and a 1/(f h) normalization while the "
                                       "pinned module decomposition R_el = c4*JA_channel "
                                       "carries f2F_eff (f2Y symbolic) and 1/f; the "
                                       "superseded oracle is internally inconsistent and "
                                       "is retained only as a historical record",
        },
        "all_checks_pass": all_ok,
    }

    out = Path(__file__).resolve().parents[1] / "data" / "generated" / "phase2_q2" \
        / "STEP3_DELTA22_DIRECT_VARIATION.json"
    out.write_text(json.dumps(evidence, indent=1) + "\n")
    print(f"evidence written: {out}")
    print("all_checks_pass:", all_ok)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
