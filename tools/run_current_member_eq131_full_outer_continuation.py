#!/usr/bin/env python3
"""Finite Eq131-compatible continuation across the *entire* native outer tail.

This supersedes the earlier narrow trial, whose active domain stopped at
u=0.708 even though the frozen member extends to u~=0.709985 and reaches its
worst kinetic ghost at the endpoint.

The present diagnostic therefore covers every native row with u>=0.70035.
It modifies only the outer tail, preserves the registered production window,
relinearizes after every accepted step, enforces E00/E11/JA/Eq131 at multiple
collocation points, and verifies K on *all* outer-tail rows for all required L.

It also exports the finite trial action and direct41 stream so downstream
radial/interface audits can be performed without reconstructing the result.

No promotion or production overwrite is performed.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.numerics import module  # noqa:E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa:E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa:E402
from ssz_p5.production.central_action import background_residuals  # noqa:E402

OUTDIR = ROOT / "data/generated/spectral"
REPORT = OUTDIR / "CURRENT_MEMBER_EQ131_FULL_OUTER_CONTINUATION.json"
ACTION_OUT = OUTDIR / "CURRENT_MEMBER_EQ131_FULL_OUTER_CONTINUATION_ACTION.csv"
STREAM_OUT = OUTDIR / "CURRENT_MEMBER_EQ131_FULL_OUTER_CONTINUATION_DIRECT41.csv"

DIRECTIONS = ("f2", "f2X", "f2F", "f3", "f3X", "f4", "f4X", "f4XX")
CENTERS = (0.7013, 0.7030, 0.7048, 0.7062, 0.7076, 0.7088, 0.70955)
WIDTH = 0.00105
LS = (6, 12, 20, 42, 110, 420, 1000)
FD_EPS = 2e-5
MAX_IT = 9
TRUST = 0.06
TAIL_LO = 0.70035


def bump(u, u0):
    z = (np.asarray(u, float) - u0) / WIDTH
    q = np.zeros_like(z)
    m = np.abs(z) < 1.0
    q[m] = np.exp(-1.0 / (1.0 - z[m] ** 2)) / np.exp(-1.0)
    return q


def reemit(action):
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = action.copy().sort_values("x").reset_index(drop=True)
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    h = d.h.to_numpy(float)
    X = d.X.to_numpy(float)
    A = d.A0prime.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    F = h * A * A / (2.0 * f)
    Y = 4.0 * X * F
    Xp = zk.dr(r, X, 1, 9, 8)
    Fp = zk.dr(r, F, 1, 9, 8)
    Yp = zk.dr(r, Y, 1, 9, 8)
    d["f2phi"] = (
        zk.dr(r, d.f2.to_numpy(float), 1, 9, 8)
        - d.f2X.to_numpy(float) * Xp
        - d.f2F.to_numpy(float) * Fp
        - d.f2Y.to_numpy(float) * Yp
    ) / ph
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
        zk.dr(r, out.a2.to_numpy(float), 1, 9, 8)
        - zk.dr(r, out.a1.to_numpy(float), 2, 9, 8)
        - zk.dr(r, out.A0prime.to_numpy(float) * out.v4.to_numpy(float) / 2, 1, 9, 8)
        + out.A0prime.to_numpy(float) * out.v5.to_numpy(float) / 2
    )
    return d, out


def apply_params(base, p):
    d = base.copy()
    u = d.u.to_numpy(float)
    p = np.asarray(p, float).reshape(len(CENTERS), len(DIRECTIONS))
    for ic, c in enumerate(CENTERS):
        sh = bump(u, c)
        for j, name in enumerate(DIRECTIONS):
            d[name] = d[name].to_numpy(float) + p[ic, j] * sh
    return d


def constraint_vector(completed, colloc):
    bg = background_residuals(completed)
    vals = []
    for i in colloc:
        vals.extend((bg["E00"][i], bg["E11"][i], bg["JA"][i], bg["eq131_residual"][i]))
    return np.asarray(vals, float)


def kinetic_vector(stream, red, rows):
    vals, tags = [], []
    for L in LS:
        a = red.canonical_audit(stream, int(L))
        K = np.asarray(a["K"], float)
        Ks = (K + K.swapaxes(1, 2)) / 2.0
        ev = np.linalg.eigvalsh(Ks)[:, 0]
        for i in rows:
            vals.append(ev[i])
            tags.append((L, i))
    return np.asarray(vals, float), tags


def worst_info(vals, tags, u):
    q = int(np.argmin(vals))
    L, i = tags[q]
    return {"min_K": float(vals[q]), "L": int(L), "u": float(u[i]), "row": int(i)}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def radial_outer_summary(stream, red, outer_rows):
    out = {}
    for L in LS:
        a = red.canonical_audit(stream, int(L))
        K = np.asarray(a["K"], float)
        G = np.asarray(a["G"], float)
        Ks = (K + K.swapaxes(1, 2)) / 2
        Gs = (G + G.swapaxes(1, 2)) / 2
        cr = []
        for i in outer_rows:
            w, U = np.linalg.eigh(Ks[i])
            if np.min(w) <= 0:
                continue
            inv = U @ np.diag(1 / np.sqrt(w)) @ U.T
            C = inv @ Gs[i] @ inv
            cr.append(float(np.linalg.eigvalsh((C + C.T) / 2)[0]))
        out[str(L)] = {
            "defined_rows": len(cr),
            "min_cr2_where_K_positive": min(cr) if cr else None,
            "radial_pass_where_defined": bool(cr and min(cr) > 0),
        }
    return out


def main():
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    build = build_onshell_central(ROOT)
    base = build.action.sort_values("x").reset_index(drop=True)
    u = base.u.to_numpy(float)
    umax = float(u.max())

    outer = np.flatnonzero(u >= TAIL_LO)
    if len(outer) < 100:
        raise RuntimeError("insufficient full outer-tail rows")

    colloc = np.array([int(np.argmin(np.abs(u - c))) for c in CENTERS], dtype=int)
    # K objective sampling includes endpoint and a dense uniform tail coverage.
    pick = np.linspace(0, len(outer) - 1, 72).round().astype(int)
    sample = np.unique(np.concatenate([outer[pick], colloc, [outer[-1]]]))

    npar = len(CENTERS) * len(DIRECTIONS)
    p = np.zeros(npar)
    history = []

    comp, stream = reemit(apply_params(base, p))
    E = constraint_vector(comp, colloc)
    kvals, tags = kinetic_vector(stream, red, sample)
    initial_sample = worst_info(kvals, tags, u)
    initial_full_vals, initial_full_tags = kinetic_vector(stream, red, outer)
    initial_full = worst_info(initial_full_vals, initial_full_tags, u)

    for it in range(MAX_IT):
        B = np.zeros((len(E), npar))
        G = np.zeros((len(kvals), npar))
        for j in range(npar):
            pp = p.copy()
            pp[j] += FD_EPS
            cp, sp = reemit(apply_params(base, pp))
            Ep = constraint_vector(cp, colloc)
            kp, _ = kinetic_vector(sp, red, sample)
            B[:, j] = (Ep - E) / FD_EPS
            G[:, j] = (kp - kvals) / FD_EPS

        rn = np.linalg.norm(B, axis=1)
        keep = rn > 1e-12
        Be = B[keep] / rn[keep, None]
        ee = E[keep] / rn[keep]

        obj = np.zeros(npar + 1)
        obj[-1] = -1.0
        Aub = np.column_stack([-G, np.ones(len(kvals))])
        bub = kvals.copy()
        Aeq = np.column_stack([Be, np.zeros(len(Be))]) if len(Be) else None
        beq = -ee if len(Be) else None
        sol = linprog(
            obj,
            A_ub=Aub,
            b_ub=bub,
            A_eq=Aeq,
            b_eq=beq,
            bounds=[(-TRUST, TRUST)] * npar + [(None, None)],
            method="highs",
        )
        if not sol.success:
            history.append({"iteration": it, "status": "LP_FAIL", "message": sol.message})
            break

        dp = sol.x[:-1]
        curmin = float(kvals.min())
        accepted = False
        best = None

        for frac in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125):
            pt = p + frac * dp
            ct, st = reemit(apply_params(base, pt))
            Et = constraint_vector(ct, colloc)
            kt, tagt = kinetic_vector(st, red, sample)
            emax = float(np.max(np.abs(Et)))
            wi = worst_info(kt, tagt, u)
            rec = {
                "iteration": it,
                "fraction": frac,
                "predicted_margin": float(sol.x[-1]),
                "finite_sample_worst": wi,
                "constraint_max_abs_at_collocation": emax,
                "step_max_abs": float(np.max(np.abs(frac * dp))),
            }
            score = wi["min_K"] - 10.0 * emax
            if best is None or score > best[0]:
                best = (score, pt, ct, st, Et, kt, tagt, rec)
            if (
                wi["min_K"] > curmin + max(1e-6, 1e-4 * max(1.0, abs(curmin)))
                and emax < 2e-3
            ):
                p, comp, stream, E, kvals, tags = pt, ct, st, Et, kt, tagt
                rec["status"] = "ACCEPT"
                history.append(rec)
                accepted = True
                break

        if not accepted:
            if best and best[7]["constraint_max_abs_at_collocation"] < 5e-3 and best[7]["finite_sample_worst"]["min_K"] > curmin:
                _, p, comp, stream, E, kvals, tags, rec = best
                rec["status"] = "ACCEPT_BEST"
                history.append(rec)
                accepted = True
            else:
                if best:
                    best[7]["status"] = "STALL"
                    history.append(best[7])
                break

        # Do not stop merely on sampled K; verify entire outer tail first.
        full_vals, full_tags = kinetic_vector(stream, red, outer)
        if float(full_vals.min()) > 0 and float(np.max(np.abs(E))) < 2e-3:
            break

    final_action = apply_params(base, p)
    comp, stream = reemit(final_action)

    full_vals, full_tags = kinetic_vector(stream, red, outer)
    final_full = worst_info(full_vals, full_tags, u)

    bg = background_residuals(comp)
    outer_mask = u >= TAIL_LO
    bgmax = {
        k: float(np.max(np.abs(np.asarray(v)[outer_mask])))
        for k, v in bg.items()
    }

    prod = (u > 0.62) & (u < 0.70)
    prod_shift = 0.0
    for name in DIRECTIONS:
        prod_shift = max(
            prod_shift,
            float(
                np.max(
                    np.abs(
                        final_action[name].to_numpy(float)[prod]
                        - base[name].to_numpy(float)[prod]
                    )
                )
            ),
        )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    comp.to_csv(ACTION_OUT, index=False)
    stream.to_csv(STREAM_OUT, index=False)

    radial = radial_outer_summary(stream, red, outer)

    bg_ok = all(
        bgmax[k] < 2e-3
        for k in ("E00", "E11", "JA", "eq131_residual")
    )
    if final_full["min_K"] > 0 and bg_ok:
        status = "FULL_NATIVE_OUTER_EQ131_K_PASS"
        nxt = "AUDIT_OUTER_RADIAL_G_AND_INTERFACE_THEN_TACKLE_INNER_CR2"
    elif final_full["min_K"] > initial_full["min_K"]:
        status = "FULL_NATIVE_OUTER_EQ131_K_IMPROVED_NOT_CLOSED"
        nxt = "REFINE_OUTER_BASIS_TRUST_AND_NONLINEAR_CORRECTOR"
    else:
        status = "FULL_NATIVE_OUTER_EQ131_K_STALLED"
        nxt = "EXPAND_ACTION_BASIS_OR_NONLINEAR_CORRECTOR"

    report = {
        "scope": "finite Eq131-compatible continuation over every native outer-tail row; no promotion",
        "native_u_max": umax,
        "outer_domain": [TAIL_LO, umax],
        "outer_rows": int(len(outer)),
        "directions": list(DIRECTIONS),
        "centers": list(CENTERS),
        "width": WIDTH,
        "initial_sample_worst": initial_sample,
        "initial_full_outer_worst": initial_full,
        "iterations": history,
        "final_full_outer_worst": final_full,
        "background_max_abs_outer": bgmax,
        "production_action_max_abs_shift": prod_shift,
        "parameter_max_abs": float(np.max(np.abs(p))),
        "nonzero_parameter_count": int(np.sum(np.abs(p) > 1e-10)),
        "outer_radial_summary": radial,
        "artifacts": {
            "action_csv": str(ACTION_OUT.relative_to(ROOT)),
            "action_sha256": sha256(ACTION_OUT),
            "direct41_csv": str(STREAM_OUT.relative_to(ROOT)),
            "direct41_sha256": sha256(STREAM_OUT),
        },
        "status": status,
        "next_target": nxt,
        "guard": (
            "Even FULL_NATIVE_OUTER_EQ131_K_PASS would certify only the current "
            "native outer tail. The known inner-tail radial cr2 failure and the "
            "missing center-to-infinity direct-global KRGM remain separate gates."
        ),
    }
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
