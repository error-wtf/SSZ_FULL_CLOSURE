#!/usr/bin/env python3
"""N2 V4: luminal background solve with a VARYING G4(phi) law (a1 != 0).

Implements the declared next step of docs/N1_N4_LUMINAL_RESOLVE_SPEC.md
section (f) — the only route to a nontrivial on-shell background after the
V3 no-go (bd46255: constant G4 leaves only the trivial branch).

Law (DECLARED here, first implementation):
    G4 = 1/2 + a1 (phi - 1),  a1 != 0,  G4phi = a1, G4phiphi = 0
    G4X = G4XX = G4phiX = 0,  G3-family = G5-family = 0
    G2 = G2X X with G2X = 1 (X = -h phi'^2/2), G2F = 1, A0' = 0

Architecture (exactly the sympy-verified V3 structure, extended to the
varying law; repo-symbol convention: `ph` is phi', NOT the field value;
the field value enters ONLY through the declared law as symbol PHI):

    E00 = 0  ->  h'   (algebraic, hp-linear, no fp)
    E11 = 0  ->  f'   (algebraic, fp-linear only)
    E22 = 0  ->  phi'' (algebraic after substituting h', f' and
                 f'' = total d/dr of the f' solution — the V3-certified
                 Bianchi-shadow closure; holds for any diffeo-invariant
                 law via paper Eq. (2.16) on the E00/E11 manifold with
                 E_A0 = 0, visually verified on canonical page 4)

Sanity anchors (runtime, fail-closed):
    * a1 = 0 reproduces the V3 exact_rhs (h', f', phi'') to machine
      precision on random states — the constant-law no-go system.
    * eps = 0 (no hair) reproduces the trivial branch exactly.

Family: shooting from the flat end u = 100 (r = r_min = 0.01) with
    f = 1/4, h = 1, phi = 1, phi'(r_min) = eps  (hair amplitude).
    The archive core edge (u = 0.71) is the terminal point.

N3/N4: emission through the VERIFIED chain (quartic_g5zero_primitives +
emit_from_primitives) on certified branches; K_scalar from the emitted
stream; decision per data/diagnostic/N1_N4_DECISION_RULES.json —
UNCHANGED (H1_CONFIRMED / H3_PHYSICAL_GHOST / NEW_PATHOLOGY_MAP).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.action.kt_mh_background import (  # noqa: E402
    G2,
    G2F,
    G2X,
    G3,
    G3X,
    G4,
    G4X,
    G4XX,
    G5,
    G5X,
    G5XX,
    G3phi,
    G4phi,
    G4phiX,
    G5phi,
    G5phiphi,
    G5phiX,
    app,
    kt_e00,
    kt_e11,
    kt_e22,
    ppp,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    f as sf,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    fp as sfp,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    fpp as sfpp,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    h as sh,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    hp as shp,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    ph as sph,
)
from ssz_p5.action.kt_mh_background import (  # noqa: E402
    r as sr,
)

ARCHIVE = ROOT / (
    "data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
)
OUT_DIR = ROOT / "data/generated/spectral"
DERIV_CONVENTION = (
    "d/dr at r = r_over_rs = 1/u (r_s = 1); repo jet-symbol convention: "
    "kt_mh_background `ph` is phi'; the scalar FIELD value enters only via "
    "the law symbol PHI"
)

PHI = sp.Symbol("phi_field")   # scalar FIELD value (law argument)
A1 = sp.Symbol("a1")           # law coupling

LAW = {
    G3: 0, G3X: 0, G3phi: 0,
    G4X: 0, G4phiX: 0, G4XX: 0,
    G5: 0, G5X: 0, G5phi: 0, G5phiX: 0, G5XX: 0, G5phiphi: 0,
    G4: sp.Rational(1, 2) + A1 * (PHI - 1),
    G4phi: A1,
    G2F: 1,
    G2X: 1,
    app: 0,
    G2: -sh * sph**2 / 2,
}

_RED_ARGS = (sr, PHI, sf, sh, sph, A1)


def build_reduction():
    """sympy reduction: closed-form h', f', phi'' for the varying law."""
    E00 = kt_e00().subs(LAW)
    E11 = kt_e11().subs(LAW)
    E22 = kt_e22().subs(LAW)

    # h' from E00 (hp-linear, no fp — structural facts checked)
    coef_hp = sp.diff(E00, shp)
    rest_hp = E00.subs(shp, 0)
    facts = {
        "E00_hp_linear": coef_hp.free_symbols.isdisjoint({shp}),
        "E00_free_of_fp": sfp not in E00.free_symbols,
        "E11_fp_linear_only": E11.free_symbols.isdisjoint({shp, ppp, sfpp}),
    }
    hp_of = sp.cancel(sp.together(-rest_hp / coef_hp))

    # f' from E11
    coef_fp = sp.diff(E11, sfp)
    rest_fp = E11.subs(sfp, 0)
    facts["E11_fp_linear"] = coef_fp.free_symbols.isdisjoint({sfp})
    fp_of = sp.cancel(sp.together(-rest_fp / coef_fp))

    # f'' = total d/dr of the f' solution — EXPLICIT metric chain rule
    # (V3 pattern: d/dsr + d/df*f' + d/dh*h' + d/dph*phi'').  dr_total()
    # is NOT used: it only carries the jet chain.  NO sp.expand anywhere
    # here — expand on these nested radicals is prohibitively slow;
    # cancel/together + subs stay in the sub-second regime.
    dF1_full = (
        sp.diff(fp_of, sr)
        + sp.diff(fp_of, sf) * sfp
        + sp.diff(fp_of, sh) * shp
        + sp.diff(fp_of, sph) * ppp
        + sp.diff(fp_of, PHI) * sph
    )
    dF1 = dF1_full.subs({shp: hp_of, sfp: fp_of})
    A = sp.cancel(sp.together(dF1.subs(ppp, 0)))
    B = sp.cancel(sp.together(sp.diff(dF1, ppp)))
    facts["fpp_linear_in_phpp"] = B.free_symbols.isdisjoint({ppp})

    # phi'' from E22 on the reduced manifold.  hp_of may contain ppp
    # (a1 != 0 coupling); substitute the LINEAR hp = alpha*ppp + beta
    # form explicitly, then E22 must reduce to phpp-linear (Bianchi
    # shadow generalizes).  All operations cancel/together-class only.
    hp_alpha = sp.cancel(sp.together(sp.diff(hp_of, ppp)))
    hp_beta = sp.cancel(sp.together(hp_of.subs(ppp, 0)))
    E22_sub = E22.subs({
        sfp: fp_of, sfpp: A + B * ppp, shp: hp_alpha * ppp + hp_beta,
    })
    ppp_in = ppp in E22_sub.free_symbols
    facts["coupled_E00_E22_had_ppp"] = ppp_in
    coef_pp = sp.diff(E22_sub, ppp)
    rest_pp = E22_sub.subs(ppp, 0)
    facts["E22_phpp_linear_after_reduction"] = (
        coef_pp.free_symbols.isdisjoint({ppp})
    )
    ppp_of = sp.cancel(sp.together(-rest_pp / coef_pp))

    # fail-closed: closed-form solutions must be FREE of the metric-chain
    # symbols except the LEGITIMATE ppp coupling in hp_of (hp = a*ppp+b,
    # composed numerically in rhs).  phpp IS ppp (repo alias).
    chain_syms = {shp, sfp, sfpp, sp.Symbol("hpp"), sp.Symbol("fppp"),
                  sp.Symbol("app2")}
    for nm, ex in (("hp_of", hp_of), ("fp_of", fp_of), ("ppp_of", ppp_of)):
        bad = ex.free_symbols & chain_syms
        facts[f"{nm}_free_of_chain_symbols"] = len(bad) == 0
        if bad:
            print(f"      WARNING {nm} contains chain symbols: {bad}")
    facts["hp_of_couples_only_ppp"] = (
        hp_of.free_symbols - chain_syms
    ).isdisjoint({sp.Symbol("hpp"), sp.Symbol("fppp"), shp, sfp, sfpp}) \
        and ppp in hp_of.free_symbols
    # denominators nonzero on the physical branch (symbolic sanity)
    facts["coef_hp_symbolic_nonzero"] = sp.cancel(
        sp.together(hp_alpha)
    ) != 0
    facts["coef_phpp_symbolic_nonzero"] = sp.cancel(
        sp.together(coef_pp)
    ) != 0

    lf_halpha = sp.lambdify(_RED_ARGS, hp_alpha, "numpy")
    lf_hbeta = sp.lambdify(_RED_ARGS, hp_beta, "numpy")
    lf_fp = sp.lambdify(_RED_ARGS, fp_of, "numpy")
    lf_pp = sp.lambdify(_RED_ARGS, ppp_of, "numpy")
    lf_cd = sp.lambdify(_RED_ARGS, coef_pp, "numpy")  # fail-closed denominator
    return (lf_halpha, lf_hbeta, lf_fp, lf_pp, lf_cd), facts, {
        "hp_of(form)": "hp = hp_alpha * phpp + hp_beta (composed in rhs)",
        "fp_of": str(fp_of), "ppp_of": str(ppp_of),
    }


def rhs_factory(red, a1_val):
    lf_halpha, lf_hbeta, lf_fp, lf_pp, lf_cd = red

    def rhs(r, y):
        f, h, phi, php = y
        if f <= 0 or h <= 0 or not np.isfinite(y).all():
            raise FloatingPointError(f"bad state at r={r}: {y}")
        args = (r, phi, f, h, php, a1_val)
        # coef_pp is proportional to phi' — the ONLY degenerate point is
        # phi' == 0 exactly (the trivial branch, where phi'' == 0 is the
        # solution: E22's rest vanishes at the same order).  Fail-closed
        # elsewhere.
        den = float(np.squeeze(lf_cd(*args)))
        if not np.isfinite(den) or abs(den) < 1e-14:
            if abs(php) < 1e-12:
                hp = float(np.squeeze(lf_hbeta(*args)))
                return [float(np.squeeze(lf_fp(*args))), hp, php, 0.0]
            raise ZeroDivisionError(
                f"E22 phi''-coefficient ~0 at r={r} with phi'={php}"
            )
        ppp = float(np.squeeze(lf_pp(*args)))
        hp = (
            float(np.squeeze(lf_halpha(*args))) * ppp
            + float(np.squeeze(lf_hbeta(*args)))
        )
        return [
            float(np.squeeze(lf_fp(*args))),
            hp,
            php,
            ppp,
        ]

    return rhs


def integrate_branch(rhs, r_from, r_to, eps, *, rtol=1e-12, atol=1e-14):
    def blowup(r, y):
        return min(abs(y[0]) - 1e8, abs(y[1]) - 1e8, abs(y[2]) - 10.0)

    blowup.terminal = True
    blowup.direction = 1.0
    return solve_ivp(
        rhs, (r_from, r_to), [0.25, 1.0, 1.0, eps], method="LSODA",
        rtol=rtol, atol=atol, events=blowup, dense_output=True,
    )


def residuals_on_grid(r, f, h, phi, a1_val):
    """Independent certification: full E00/E11/E22 with jet9d8 derivatives
    (differentiation path, not the integration path)."""
    from ssz_p5.jets.jet9d8 import derivative as jet

    fp = jet(r, f)
    hp = jet(r, h)
    php = jet(r, phi)
    fpp = jet(r, f, 2)
    phpp = jet(r, phi, 2)
    args_syms = (sr, PHI, sf, sh, sph, sfp, sfpp, shp, ppp, A1)
    out = {}
    for name, eq in (("E00", kt_e00), ("E11", kt_e11), ("E22", kt_e22)):
        fn = sp.lambdify(args_syms, eq().subs(LAW), "numpy")
        vals = fn(r, phi, f, h, php, fp, fpp, hp, phpp, a1_val)
        out[name] = np.asarray(vals, float)
    return out


def _stats(v):
    a = np.abs(np.asarray(v, float))
    return {"max": float(a.max()), "median": float(np.median(a))}


def _v3_phpp(r, f, h, php):
    """V3 exact phi'': from (r^2 sqrt(fh) phi')' = 0 =>
    phi'' = -phi' (2/r + f'/(2f) + h'/(2h)) with the V3 f', h'."""
    X = -h * php**2 / 2
    fp = f * (1 - h) / (h * r) - f * X * r / h
    hp = (1 - h) / r + r * X
    return -php * (2 / r + fp / (2 * f) + hp / (2 * h))


def v3_equivalence_check(red, n=5):
    """a1=0 must reproduce the V3 exact_rhs (+ closed phi'') on random
    states: with PHI irrelevant at a1=0, the reduction must collapse to
    the V3 system. Note V3's phi'' is NOT the J0-conservation form at
    arbitrary php — V3 enforces J-conservation; here we compare the
    (f', h') components exactly and phi'' against the J0-consistent
    value for the php that the state carries."""
    rng = np.random.default_rng(20261006)
    worst = 0.0
    for _ in range(n):
        r = float(rng.uniform(0.05, 1.2))
        f = float(rng.uniform(0.1, 0.9))
        h = float(rng.uniform(0.3, 1.0))
        php = float(rng.uniform(-0.2, 0.2))
        rhs = rhs_factory(red, 0.0)
        x = rhs(r, [f, h, 1.0, php])
        X = -h * php**2 / 2
        v3 = [
            f * (1 - h) / (h * r) - f * X * r / h,
            (1 - h) / r + r * X,
            php,
            _v3_phpp(r, f, h, php),
        ]
        err = float(np.max(np.abs(np.array(x) - np.array(v3))))
        worst = max(worst, err)
    return worst


POCKET_LO, POCKET_HI = 0.285, 0.325


def main() -> int:
    t0 = time.time()
    archive = pd.read_csv(ARCHIVE).sort_values("u").reset_index(drop=True)
    u_all = archive.u.to_numpy(float)
    r_all = (1.0 / u_all)[::-1]
    r_min, r_core = float(r_all[0]), float(r_all[-1])
    n_out = len(r_all)

    print("[1/7] sympy reduction of the varying-G4 system ...", flush=True)
    red, facts, exprs = build_reduction()
    for k, v in facts.items():
        print(f"      {k}: {v}", flush=True)
    assert all(facts.values()), "structural reduction facts violated"

    print("[2/7] a1=0 equivalence vs V3 exact_rhs ...", flush=True)
    worst = v3_equivalence_check(red)
    print(f"      max|V4(a1=0) - V3| over random states: {worst:.3e}")
    assert worst < 1e-9, "a1=0 limit does not reproduce V3"

    print("[3/7] eps = 0 -> trivial branch recovery ...", flush=True)
    sol0 = integrate_branch(rhs_factory(red, 0.5), r_min, r_core, 0.0)
    y0 = sol0.sol(r_all)
    recov = {
        "trivial_recovered": bool(
            np.max(np.abs(y0[0] - 0.25)) < 1e-8
            and np.max(np.abs(y0[1] - 1.0)) < 1e-8
            and np.max(np.abs(y0[2] - 1.0)) < 1e-8
        ),
        "max|f-1/4|": float(np.max(np.abs(y0[0] - 0.25))),
        "max|h-1|": float(np.max(np.abs(y0[1] - 1.0))),
    }
    print(f"      trivial recovered: {recov['trivial_recovered']}")
    assert recov["trivial_recovered"]

    print("[4/7] eps-family scan (hair amplitude) ...", flush=True)
    a1_grid = [0.5, -0.5, 2.0, -2.0]
    eps_grid = np.concatenate([
        [0.0],
        np.geomspace(1e-8, 1e-1, 10),
        np.linspace(0.15, 1.5, 10),
        -np.geomspace(1e-8, 1e-1, 10),
        -np.linspace(0.15, 1.5, 10),
    ])
    rows = []
    for a1v in a1_grid:
        rhs = rhs_factory(red, a1v)
        for eps in eps_grid:
            rec = {"a1": a1v, "eps": float(eps)}
            try:
                sol = integrate_branch(rhs, r_min, r_core, float(eps))
            except (FloatingPointError, ZeroDivisionError,
                    np.linalg.LinAlgError, ValueError) as exc:
                rec.update(reached_core=False, reason=f"{type(exc).__name__}")
                rows.append(rec)
                continue
            reached = bool(sol.t[-1] >= r_core * (1 - 1e-9))
            rec["reached_core"] = reached
            if reached:
                fe, he, phie, phpe = sol.y[:, -1]
                rec.update(
                    f_at_core=float(fe), h_at_core=float(he),
                    phi_at_core=float(phie), phip_at_core=float(phpe),
                    h_min=float(sol.y[1].min()), f_min=float(sol.y[0].min()),
                    phi_range=[float(sol.y[2].min()), float(sol.y[2].max())],
                )
            else:
                rec["reason"] = f"blowup_at_r={float(sol.t[-1]):.4f}"
            rows.append(rec)
    fam = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fam_path = OUT_DIR / "N2_LUMINAL_V4_EPS_FAMILY.csv"
    fam.to_csv(fam_path, index=False)
    n_reached = int(fam.reached_core.fillna(False).sum())
    print(f"      {n_reached}/{len(fam)} branches reached the core edge")

    print("[5/7] archive-core matching + certification of reaching "
          "branches ...", flush=True)
    a_fc = float(archive.A_f[0])
    b_hc = float(archive.B_h[0])
    ph_c = float(archive.phi[0])
    cand = fam[fam.reached_core == True].copy()  # noqa: E712
    if len(cand):
        cand["core_mismatch"] = (
            (cand.f_at_core - a_fc) ** 2
            + (cand.h_at_core - b_hc) ** 2
            + (cand.phi_at_core - ph_c) ** 2
        ) ** 0.5
        cand = cand.sort_values("core_mismatch")
    certified = []
    # Interior certification window: jet9d8 8th-degree stencils are
    # one-sided at the domain edges (documented V3 window precedent);
    # the certification excludes the outermost 8 points on each side and
    # says so.  Full-domain residual stats are reported separately.
    MARGIN = 8
    interior = slice(MARGIN, n_out - MARGIN)
    for _, row in cand.head(20).iterrows():
        a1v, eps = float(row.a1), float(row.eps)
        try:
            sol = integrate_branch(rhs_factory(red, a1v), r_min, r_core, eps)
            y = sol.sol(r_all)
            res = residuals_on_grid(r_all, y[0], y[1], y[2], a1v)
            st = {
                k: {
                    "max": float(np.max(np.abs(v))),
                    "max_interior": float(np.max(np.abs(v[interior]))),
                }
                for k, v in res.items()
            }
            ok = (
                st["E00"]["max_interior"] <= 1e-6
                and st["E11"]["max_interior"] <= 1e-6
                and st["E22"]["max_interior"] <= 1e-4
            )
            certified.append({
                "a1": a1v, "eps": eps,
                "converged": bool(ok),
                "residual_stats": st,
                "profile": y,
                "core_mismatch": float(row.core_mismatch),
            })
            print(
                f"      a1={a1v} eps={eps:.3e}: E00 int {st['E00']['max_interior']:.2e} "
                f"(full {st['E00']['max']:.2e}) E11 int {st['E11']['max_interior']:.2e} "
                f"E22 int {st['E22']['max_interior']:.2e} conv={ok} "
                f"mismatch={row.core_mismatch:.4f}",
                flush=True,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"      a1={a1v} eps={eps:.3e}: certification failed: {exc}")
    certified = [c for c in certified if c["converged"]]
    certified.sort(key=lambda c: c["core_mismatch"])

    print("[6/7] N3 emission + direct K_scalar on certified nontrivial "
          "branches ...", flush=True)
    from ssz_p5.jets.jet9d8 import derivative as jet

    emission = []
    for c in certified:
        if abs(c["eps"]) < 1e-14:
            continue
        a1v, eps = c["a1"], c["eps"]
        f_p, h_p, phi_p = c["profile"][0], c["profile"][1], c["profile"][2]
        # ---- direct on-shell K via the V3-certified identity
        #      K = F r/2 (f'/f - h'/h),  F = 2 G4(phi) = 1 + 2 a1 (phi-1)
        fp_g = jet(r_all, f_p)
        hp_g = jet(r_all, h_p)
        F_of_r = 1.0 + 2.0 * a1v * (phi_p - 1.0)
        K_dir = F_of_r * r_all / 2.0 * (fp_g / f_p - hp_g / h_p)
        # stencil-resolution proof: w9 vs w13 (and loose vs tight integration
        # tolerance) must agree far below the |K| scale of the verdict
        fp_w = jet(r_all, f_p, window=13)
        hp_w = jet(r_all, h_p, window=13)
        K_w13 = F_of_r * r_all / 2.0 * (fp_w / f_p - hp_w / h_p)
        k_stencil_diff = float(np.max(np.abs(K_dir - K_w13)))
        band = (r_all >= POCKET_LO) & (r_all <= POCKET_HI)
        entry = {
            "a1": a1v, "eps": eps,
            "K_direct_min": float(np.min(K_dir)),
            "K_direct_max": float(np.max(K_dir)),
            "K_pocket_band_min": float(np.min(K_dir[band])),
            "K_pocket_band_max": float(np.max(K_dir[band])),
            "K_stencil_diff_w9_w13": k_stencil_diff,
            "residual_scale_E00_interior":
                c["residual_stats"]["E00"]["max_interior"],
            "sign_structure": (
                "K<0 everywhere (incl. pocket band)"
                if np.max(K_dir) < 0
                else "K>0 everywhere" if np.min(K_dir) > 0
                else "K sign-changing"
            ),
            "stream41_emitted": False,
        }
        # ---- 41-slot emission through the verified chain (executed)
        try:
            from ssz_p5.coefficients.mh_action_primitives import (
                quartic_g5zero_primitives,
            )
            from ssz_p5.coefficients.mh_general_primitives import (
                emit_from_primitives,
            )

            phi_r = jet(r_all, phi_p)
            feed = pd.DataFrame({
                "x": r_all, "f": f_p, "h": h_p, "phi_r": phi_r,
                "G4": 0.5 + a1v * (phi_p - 1.0),
                "G4X": np.zeros(n_out), "G4XX": np.zeros(n_out),
                "G4phi": np.full(n_out, a1v),
                "G4phiX": np.zeros(n_out), "G3X": np.zeros(n_out),
            })
            prim = quartic_g5zero_primitives(feed)
            c2 = np.sqrt(f_p * h_p) * phi_r * 0.5 * r_all**2
            prim_in = pd.DataFrame({
                "x": r_all, "f": f_p, "h": h_p, "phiprime": phi_r,
                "A0prime": np.zeros(n_out),
                "a1": prim.a1_action.to_numpy(float),
                "c2": c2,
                "c4": prim.c4_action.to_numpy(float),
                "F_tensor": prim.F_tensor_action.to_numpy(float),
                "G_tensor": prim.G_tensor_action.to_numpy(float),
                "H_tensor": prim.H_tensor_action.to_numpy(float),
            })
            stream = emit_from_primitives(prim_in)
            nonfinite = {
                col: int((~np.isfinite(stream[col].to_numpy(float))).sum())
                for col in stream.columns
                if stream[col].dtype.kind == "f"
            }
            entry["stream41_emitted"] = True
            entry["nonfinite_slot_rows"] = {
                k: v for k, v in nonfinite.items() if v
            }
            sname = (
                f"N3_STREAM41_V4_a1_{a1v}_eps_{eps:.2e}.csv"
                .replace("-", "m").replace("+", "p")
            )
            stream.to_csv(OUT_DIR / sname, index=False)
            entry["stream_csv"] = str((OUT_DIR / sname).relative_to(ROOT))
        except Exception as exc:  # noqa: BLE001
            entry["stream41_emitted"] = False
            entry["stream_exception"] = f"{type(exc).__name__}: {exc}"
        emission.append(entry)
        print(
            f"      a1={a1v} eps={eps:.3e}: K in "
            f"[{entry['K_direct_min']:.4e}, {entry['K_direct_max']:.4e}] "
            f"({entry['sign_structure']}), stencil-diff {k_stencil_diff:.1e}, "
            f"stream41={entry['stream41_emitted']}",
            flush=True,
        )

    print("[7/7] N4 decision per UNCHANGED declared rules ...", flush=True)
    rules = json.loads(
        (ROOT / "data/diagnostic/N1_N4_DECISION_RULES.json").read_text()
    )
    decisions = []
    for e in emission:
        kmin_glob = e["K_direct_min"]
        kmin_band = e["K_pocket_band_min"]
        # fail-closed clause: K's error scale = max(stencil diff,
        # background residual scale).  A negative K whose |min| is below
        # that scale is NOT decidable -> NOT_EVALUABLE (no new verdicts).
        k_err = max(e["K_stencil_diff_w9_w13"],
                    e["residual_scale_E00_interior"])
        if abs(kmin_glob) < k_err and abs(kmin_band) < k_err:
            decisions.append({
                "a1": e["a1"], "eps": e["eps"],
                "verdict": "NOT_EVALUABLE",
                "reason": "|K| below error scale (stencil/residual)",
                "k_error_scale": k_err,
            })
            print(f"      a1={e['a1']} eps={e['eps']:.3e}: NOT_EVALUABLE "
                  f"(|K| < {k_err:.1e})")
            continue
        if kmin_band < 0:
            verdict = "H3_PHYSICAL_GHOST"
        elif kmin_glob < 0:
            verdict = "NEW_PATHOLOGY_MAP"
        else:
            verdict = "H1_CANDIDATE_PENDING_FINITE_L_HEALTH"
        decisions.append({
            "a1": e["a1"], "eps": e["eps"], "verdict": verdict,
            "K_scalar_min_global": kmin_glob,
            "K_scalar_min_pocket_band": kmin_band,
            "k_error_scale": k_err,
            "sign_structure": e["sign_structure"],
        })
        print(f"      a1={e['a1']} eps={e['eps']:.3e}: {verdict}")

    result = {
        "audit": "N2_LUMINAL_SOLVE_V4",
        "law": {
            "form": "G4 = 1/2 + a1 (phi - 1)",
            "a1_grid": a1_grid,
            "G4phi": "a1", "G4phiphi": 0,
            "sectors_zero": "G4X-fam, G3-fam, G5-fam",
            "G2": "G2X X, G2X=1", "G2F": 1, "A0prime": 0,
            "declaration_status": (
                "branch decision per spec section (f); decision rules "
                "UNCHANGED (anti-circularity echo attached)"
            ),
        },
        "derivative_convention": DERIV_CONVENTION,
        "structural_facts": facts,
        "reduced_expressions": exprs,
        "v3_equivalence_max_err": worst,
        "eps_zero_recovery": recov,
        "domain": {
            "r_range": [r_min, r_core],
            "u_range": [float(u_all[0]), float(u_all[-1])],
            "points": int(n_out),
        },
        "family_scan": {
            "csv": str(fam_path.relative_to(ROOT)),
            "branches_reached_core": n_reached,
            "total": int(len(fam)),
        },
        "certified_branches": [
            {k: v for k, v in c.items() if k != "profile"}
            for c in certified
        ],
        "n3_emission": emission,
        "n4_decisions": decisions,
        "decision_rules_unchanged": True,
        "decision_rules_echo": rules["decision_rules"],
        "wall_seconds": round(time.time() - t0, 1),
    }

    if certified:
        best = certified[0]
        prof = pd.DataFrame({
            "u": u_all[::-1].copy(),
            "r_over_rs": r_all,
            "f": best["profile"][0], "h": best["profile"][1],
            "phi": best["profile"][2],
            "a1": best["a1"], "eps": best["eps"],
            "derivative_convention": DERIV_CONVENTION,
        })
        prof_path = OUT_DIR / "N2_LUMINAL_BACKGROUND_PROFILE_V4.csv"
        prof.to_csv(prof_path, index=False)
        result["profile_csv"] = str(prof_path.relative_to(ROOT))
        result["best_branch"] = {
            "a1": best["a1"], "eps": best["eps"],
            "core_mismatch": best["core_mismatch"],
            "residual_stats": best["residual_stats"],
        }

    out_path = OUT_DIR / "N2_LUMINAL_SOLVE_RESULT_V4.json"
    out_path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    dec = {
        "audit": "N3_N4_LUMINAL_DECISION_V4",
        "source_result": str(out_path.relative_to(ROOT)),
        "n4_decisions": decisions,
        "decision_rules_echo": rules["decision_rules"],
        "anti_circularity": (
            "rules unchanged from data/diagnostic/N1_N4_DECISION_RULES.json"
        ),
        "law": result["law"],
    }
    (OUT_DIR / "N3_N4_LUMINAL_DECISION_V4.json").write_text(
        json.dumps(dec, indent=1, allow_nan=False) + "\n"
    )
    print(json.dumps({
        "written": [str(out_path.relative_to(ROOT))],
        "certified": len(certified),
        "n4": [d["verdict"] for d in decisions],
        "wall_s": result["wall_seconds"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
