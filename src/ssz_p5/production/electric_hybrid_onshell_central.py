"""Reproducible on-shell electric-hybrid Central action member.

This module reconstructs the post-18-Sep research checkpoint from action data:

  electric A2 re-solve
    -> (f3, f2F) from A2/JA
    -> metric background re-solve (f2, f2X)
    -> holonomic f2 Hessian restore of (v1,v4,c2)
    -> direct Appendix-A lower slots (v5,c3,e3)
    -> background-null Horndeski G2XX lift
    -> one 41-slot Central coefficient stream.

The G2 lift is represented by the genuine local action deformation

    Delta G2(phi,X) = 1/2 q(phi) [X-X_b(phi)]^2,

so Delta G2, Delta G2_X and Delta G2_phi vanish on the background.  Its full
41-slot response used here is therefore the Horndeski primitive response
(c2,c6,e2), not an isolated c2 edit.

This is a Central-member construction only.  It does not certify the two
regional action handovers, angular characteristics, global descriptor, or QNM.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from ..numerics import module
from .full_action_lower import complete_total_action_jets, emit_lower_slots

_REQUIRED_L = (6, 12, 20, 42, 110, 420, 1000)


@dataclass(frozen=True)
class CentralBuild:
    action: pd.DataFrame
    direct41: pd.DataFrame
    pre_lift41: pd.DataFrame
    diagnostics: dict


def _paths(root: Path) -> dict[str, Path]:
    b = Path(root) / "archive/full_working_snapshot"
    return {
        "profile": b / "ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv",
        "raw": b / "ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv",
        "selected": b / "ssz_p5_F2b_central_c3_e3_SELECTED_REPRESENTATIVE_2026-09-15.csv",
        "g2_rank": Path(root)
        / "data/authoritative/ssz_p5_full_horndeski_principal_control_rank_2026-09-12.csv",
    }


def g2xx_profile(u: np.ndarray) -> np.ndarray:
    """Locked background-null radial control recovered from the accepted search."""
    u = np.asarray(u, float)
    return -300.0 - 700.0 * 0.5 * (1.0 + np.tanh((u - 0.66) / 0.01))


def _load_inputs(root: Path):
    p = _paths(root)
    for path in p.values():
        if not path.exists():
            raise FileNotFoundError(path)
    prof = pd.read_csv(p["profile"])
    raw = pd.read_csv(p["raw"])
    selected = pd.read_csv(p["selected"])
    rank = pd.read_csv(p["g2_rank"])
    if not (len(prof) == len(raw) == len(selected)):
        raise ValueError("Central source grids have different lengths")
    return prof, raw, selected, rank


def _electric_a2_resolve(root: Path, prof: pd.DataFrame, raw: pd.DataFrame, selected: pd.DataFrame):
    """Resolve v6, f3, f2F, f2 and f2X on the Central support."""
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")

    r = raw.x.to_numpy(float)
    u = raw.u.to_numpy(float)
    f = raw.f.to_numpy(float)
    h = raw.h.to_numpy(float)
    ph = raw.phiprime.to_numpy(float)
    ap = raw.A0prime.to_numpy(float)
    X = raw.X.to_numpy(float)
    fp = zk.dr(r, f, 1, 9, 8)
    fpp = zk.dr(r, f, 2, 9, 8)
    hp = zk.dr(r, h, 1, 9, 8)
    app = zk.dr(r, ap, 1, 9, 8)

    f4 = prof.f4.to_numpy(float)
    f4x = prof.N4.to_numpy(float)
    a4 = raw.a4.to_numpy(float)
    f3x = prof.f3X_integrated.to_numpy(float)
    f4xx = selected.f4XX_recovered.to_numpy(float)

    band = (u >= 0.61) & (u <= 0.715)
    idx = np.flatnonzero(band)
    order = idx[np.argsort(r[idx])]
    rs = r[order]
    fields = {
        name: PchipInterpolator(rs, arr[order], extrapolate=False)
        for name, arr in (
            ("f", f),
            ("h", h),
            ("ph", ph),
            ("ap", ap),
            ("fp", fp),
            ("fpp", fpp),
            ("hp", hp),
            ("app", app),
            ("f4", f4),
            ("f4x", f4x),
            ("a4", a4),
        )
    }

    def aux(rr: float, vv: float):
        ff = float(fields["f"](rr))
        hh = float(fields["h"](rr))
        pp = float(fields["ph"](rr))
        aa = float(fields["ap"](rr))
        F4 = float(fields["f4"](rr))
        F4x = float(fields["f4x"](rr))
        c = 2.0 * hh**1.5 * aa / (rr * np.sqrt(ff))
        F3 = (vv / c + 4.0 * F4 - hh * pp**2 * F4x) / (rr * pp)
        F2F = -(
            4.0 * rr * hh * pp * F3
            + 8.0 * (1.0 - hh) * F4
            + 2.0 * hh**2 * pp**2 * F4x
        ) / rr**2
        FP = float(fields["fp"](rr))
        v10 = -np.sqrt(ff * hh) / (2.0 * rr) * (
            rr * F2F
            + 2.0 * hh * pp * F3
            + (hh * FP / ff) * (rr * pp * F3 - 4.0 * F4 + hh * pp**2 * F4x)
        )
        a7 = (1.0 - (4.0 * hh * aa**2 / ff) * F4) / (4.0 * rr**2 * np.sqrt(ff * hh))
        return F3, F2F, v10, a7

    def rhs(rr, y):
        vv = float(y[0])
        ff = float(fields["f"](rr))
        hh = float(fields["h"](rr))
        aa = float(fields["ap"](rr))
        FP = float(fields["fp"](rr))
        FPP = float(fields["fpp"](rr))
        HP = float(fields["hp"](rr))
        APP = float(fields["app"](rr))
        A4 = float(fields["a4"](rr))
        _, _, v10, a7 = aux(rr, vv)
        num = A4 * (
            2.0 * ff**2 * (rr * HP + 2.0 * hh)
            + hh * rr**2 * FP**2
            - ff * rr * (rr * FP * HP + 2.0 * hh * (rr * FPP + FP))
        ) - ff * hh * rr**2 * (
            aa * (4.0 * v10 * aa + vv * FP) + ff * vv * APP + 8.0 * a7 * ff**2
        )
        den = ff**2 * hh * rr**2 * aa
        return [num / den]

    i0 = int(np.argmin(np.abs(u - 0.61)))
    r0 = float(r[i0])
    v0 = float(raw.v6.iloc[i0])
    r1 = float(rs.min())
    sol = solve_ivp(
        rhs,
        (r0, r1),
        [v0],
        rtol=2e-11,
        atol=2e-12,
        dense_output=True,
        max_step=2e-4,
    )
    if not sol.success:
        raise RuntimeError(f"electric A2 re-solve failed: {sol.message}")

    v6 = np.full(len(raw), np.nan)
    v6[idx] = sol.sol(r[idx])[0]
    f3 = np.full(len(raw), np.nan)
    f2F = np.full(len(raw), np.nan)
    for i in idx:
        f3[i], f2F[i], _, _ = aux(float(r[i]), float(v6[i]))

    f2old = prof.f2.to_numpy(float)
    f2xold = prof.f2X.to_numpy(float)
    A = r**2 * f
    kappa = h * ph**2
    E00 = (
        r * f * hp
        - (
            f * (1.0 - h)
            + r**2 * (f * f2old - h * ap**2 * f2F)
            - 2.0 * r * h**2 * ph * ap**2 * f3
            + h * ap**2 * (4.0 * (h - 1.0) * f4 - h**2 * ph**2 * f4x)
        )
    )
    E11 = (
        r * h * fp
        - (
            f * (1.0 - h)
            + r**2 * (f * f2old + f * h * ph**2 * f2xold - h * ap**2 * f2F)
            - 2.0 * r * h**2 * ph * ap**2 * (3.0 * f3 - h * ph**2 * f3x)
            + h
            * ap**2
            * (
                4.0 * (3.0 * h - 1.0) * f4
                - h * (9.0 * h - 4.0) * ph**2 * f4x
                + h**3 * ph**4 * f4xx
            )
        )
    )
    f2 = f2old + E00 / A
    f2x = f2xold + (E11 - E00) / (A * kappa)
    return v6, f3, f2F, f2, f2x, sol


def _build_action(root: Path):
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    prof, raw, selected, rank = _load_inputs(root)
    v6, f3, f2F, f2, f2x, sol = _electric_a2_resolve(root, prof, raw, selected)

    uall = raw.u.to_numpy(float)
    ridx = np.flatnonzero((uall >= 0.61) & (uall < 0.71))
    ridx = ridx[np.argsort(raw.x.to_numpy(float)[ridx])]
    rr = raw.x.to_numpy(float)[ridx]

    d = pd.DataFrame(
        {
            "u": uall[ridx],
            "x": rr,
            "phi": raw.phi.to_numpy(float)[ridx],
            "f": raw.f.to_numpy(float)[ridx],
            "h": raw.h.to_numpy(float)[ridx],
            "phiprime": raw.phiprime.to_numpy(float)[ridx],
            "A0prime": raw.A0prime.to_numpy(float)[ridx],
            "X": raw.X.to_numpy(float)[ridx],
            "f2": f2[ridx],
            "f2X": f2x[ridx],
            "f2F": f2F[ridx],
            "f2Y": 0.0,
            "f2XX": prof.HXX.to_numpy(float)[ridx],
            "f2XF": prof.HXF.to_numpy(float)[ridx],
            "f2XY": prof.HXY.to_numpy(float)[ridx],
            "f2FF": prof.HFF.to_numpy(float)[ridx],
            "f2FY": prof.HFY.to_numpy(float)[ridx],
            "f2YY": prof.HYY.to_numpy(float)[ridx],
            "f3": f3[ridx],
            "f3X": prof.f3X_integrated.to_numpy(float)[ridx],
            "f3XX": selected.f3XX_selected.to_numpy(float)[ridx],
            "tf3": prof.tilde_f3.to_numpy(float)[ridx],
            "f4": prof.f4.to_numpy(float)[ridx],
            "f4X": prof.N4.to_numpy(float)[ridx],
            "f4XX": selected.f4XX_recovered.to_numpy(float)[ridx],
            "f4XXX": selected.f4XXX_selected.to_numpy(float)[ridx],
            "tf4": 0.0,
            "v6_A2_resolved": v6[ridx],
        }
    )

    Fbg = d.h.to_numpy(float) * d.A0prime.to_numpy(float) ** 2 / (2.0 * d.f.to_numpy(float))
    Ybg = 4.0 * d.X.to_numpy(float) * Fbg
    Xp = zk.dr(rr, d.X.to_numpy(float), 1, 9, 8)
    Fp = zk.dr(rr, Fbg, 1, 9, 8)
    d["f2phi"] = (
        zk.dr(rr, d.f2.to_numpy(float), 1, 9, 8)
        - d.f2X.to_numpy(float) * Xp
        - d.f2F.to_numpy(float) * Fp
    ) / d.phiprime.to_numpy(float)

    zeros = np.zeros(len(d))
    pre = zk.emit(d, selected_v5=zeros, selected_c3=zeros, selected_e3=zeros, v6_phi_selector="action")
    target = raw.iloc[ridx].reset_index(drop=True)
    target3 = np.column_stack((target.v1, target.v4, target.c2)).astype(float)

    cols = ("f2XX", "f2XF", "f2FF")
    response = []
    for col in cols:
        dd = d.copy()
        dd[col] = dd[col].to_numpy(float) + 1.0
        g = zk.emit(dd, selected_v5=zeros, selected_c3=zeros, selected_e3=zeros, v6_phi_selector="action")
        response.append(
            np.stack(
                (
                    g.v1.to_numpy(float) - pre.v1.to_numpy(float),
                    g.v4.to_numpy(float) - pre.v4.to_numpy(float),
                    g.c2.to_numpy(float) - pre.c2.to_numpy(float),
                ),
                axis=1,
            )
        )

    delta = np.empty((len(d), 3))
    for j in range(len(d)):
        M = np.column_stack([response[k][j] for k in range(3)])
        rhs = target3[j] - np.array((pre.v1.iloc[j], pre.v4.iloc[j], pre.c2.iloc[j]))
        delta[j] = np.linalg.solve(M, rhs)
    for k, col in enumerate(cols):
        d[col] = d[col].to_numpy(float) + delta[:, k]

    d = complete_total_action_jets(d)
    lower, d = emit_lower_slots(d)
    out = zk.emit(
        d,
        selected_v5=lower.v5.to_numpy(float),
        selected_c3=lower.c3.to_numpy(float),
        selected_e3=lower.e3.to_numpy(float),
        v6_phi_selector="action",
    )
    out["a5"] = (
        zk.dr(rr, out.a2.to_numpy(float), 1, 9, 8)
        - zk.dr(rr, out.a1.to_numpy(float), 2, 9, 8)
        - zk.dr(rr, out.A0prime.to_numpy(float) * out.v4.to_numpy(float) / 2.0, 1, 9, 8)
        + out.A0prime.to_numpy(float) * out.v5.to_numpy(float) / 2.0
    )

    # Preserve useful internal diagnostics on the action table.
    d["Fbg_action"] = Fbg
    d["Ybg_action"] = Ybg
    d["principal_delta_f2XX"] = delta[:, 0]
    d["principal_delta_f2XF"] = delta[:, 1]
    d["principal_delta_f2FF"] = delta[:, 2]
    return d, out, target, rank, sol


def _apply_g2xx_lift(action: pd.DataFrame, pre: pd.DataFrame, rank: pd.DataFrame):
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = action.copy()
    out = pre.copy()
    u = d.u.to_numpy(float)
    r = d.x.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    X = d.X.to_numpy(float)
    Xp = zk.dr(r, X, 1, 9, 8)
    q = g2xx_profile(u)

    rank = rank.sort_values("u")
    jac = np.interp(u, rank.u.to_numpy(float), rank.dc2_dG2XX.to_numpy(float))
    dc2 = jac * q

    # Store a true on-background jet representation of
    # 1/2 q(phi)[X-X_b(phi)]^2.
    d["G2_lift"] = 0.0
    d["G2X_lift"] = 0.0
    d["G2phi_lift"] = 0.0
    d["G2XX_lift"] = q
    d["G2Xphi_lift"] = -q * Xp / ph
    d["G2phiX_lift"] = d.G2Xphi_lift
    d["G2phiphi_lift"] = q * (Xp / ph) ** 2
    d["dc2_dG2XX"] = jac
    d["delta_c2_G2_lift"] = dc2

    # Full direct Appendix-A response of the pure G2XX/c2 primitive.
    out["c2"] = out.c2.to_numpy(float) + dc2
    out["c6"] = out.c6.to_numpy(float) - 0.25 * ph * dc2
    out["e2"] = out.e2.to_numpy(float) - dc2 / ph
    return d, out


def build_onshell_central(root: Path) -> CentralBuild:
    """Build the locked Central single-action candidate and its full G2 lift."""
    return _build_onshell_central_cached(str(Path(root).resolve()))


@cache
def _build_onshell_central_cached(root_string: str) -> CentralBuild:
    root = Path(root_string)
    action, pre, target, rank, sol = _build_action(root)
    action, out = _apply_g2xx_lift(action, pre, rank)

    mask = (out.u.to_numpy(float) > 0.62) & (out.u.to_numpy(float) < 0.70)
    # Principal replay is evaluated before the independent G2XX lift.
    replay = {}
    for name in ("v1", "v4", "c2"):
        a = pre.loc[mask, name].to_numpy(float)
        b = target.loc[mask, name].to_numpy(float)
        replay[name] = float(np.max(np.abs(a - b) / np.maximum(1.0, np.abs(b))))

    # Exact chain-rule diagnostics of the stored background-null G2 action jets.
    r = action.x.to_numpy(float)
    ph = action.phiprime.to_numpy(float)
    Xp = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py").dr(
        r, action.X.to_numpy(float), 1, 9, 8
    )
    chain_x = action.G2Xphi_lift.to_numpy(float) * ph + action.G2XX_lift.to_numpy(float) * Xp
    chain_phi = (
        action.G2phiphi_lift.to_numpy(float) * ph
        + action.G2phiX_lift.to_numpy(float) * Xp
    )
    diagnostics = {
        "a2_solve_success": bool(sol.success),
        "a2_nfev": int(sol.nfev),
        "rows": int(len(out)),
        "production_rows": int(mask.sum()),
        "principal_replay_max_scaled": replay,
        "g2_lift": {
            "q_min": float(np.min(action.G2XX_lift)),
            "q_max": float(np.max(action.G2XX_lift)),
            "min_abs_dc2_dG2XX": float(np.min(np.abs(action.dc2_dG2XX[mask]))),
            "max_abs_chain_G2X": float(np.max(np.abs(chain_x[mask]))),
            "max_abs_chain_G2phi": float(np.max(np.abs(chain_phi[mask]))),
        },
    }
    return CentralBuild(action=action, direct41=out, pre_lift41=pre, diagnostics=diagnostics)


def principal_audit(root: Path):
    """Finite-L Central K/G gate from exactly the stream built above."""
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    build = build_onshell_central(root)
    d = build.direct41
    production = (d.u.to_numpy(float) > 0.62) & (d.u.to_numpy(float) < 0.70)
    inds = np.flatnonzero(production)
    rows = []
    for L in _REQUIRED_L:
        a = red.canonical_audit(d, L)
        K = np.asarray(a["K"], float)
        G = np.asarray(a["G"], float)
        Ks = (K + K.swapaxes(1, 2)) / 2.0
        Gs = (G + G.swapaxes(1, 2)) / 2.0
        ke = np.linalg.eigvalsh(Ks)
        radial = np.empty(len(inds))
        for q, i in enumerate(inds):
            w, U = np.linalg.eigh(Ks[i])
            if np.min(w) <= 0:
                radial[q] = -np.inf
            else:
                invsqrt = U @ np.diag(1.0 / np.sqrt(w)) @ U.T
                C = invsqrt @ Gs[i] @ invsqrt
                radial[q] = np.linalg.eigvalsh((C + C.T) / 2.0)[0]
        rows.append(
            {
                "L": int(L),
                "min_eig_K": float(np.min(ke[inds, 0])),
                "negative_K_rows": int(np.sum(ke[inds, 0] <= 0)),
                "min_cr2": float(np.min(radial)),
                "negative_radial_rows": int(np.sum(radial <= 0)),
                "max_R_abs": float(np.max(np.abs(np.asarray(a["R"])[inds]))),
                "max_S_sym": float(np.max(np.abs(np.asarray(a["S"])[inds] + np.asarray(a["S"])[inds].swapaxes(1, 2)))),
                "pass": bool(np.all(ke[inds, 0] > 0) and np.all(radial > 0)),
            }
        )
    return build, rows
