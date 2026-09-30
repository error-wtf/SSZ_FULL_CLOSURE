#!/usr/bin/env python3
"""Local Eq131-compatible kinetic controllability audit at the common ghost.

Tests whether the current Electric action has any *action-level local tangent*
that preserves the four background conditions

    E00 = 0, E11 = 0, JA = 0, Eq131 = 0

to first order while increasing the weakest finite-L kinetic eigenvalue.

This is the natural compatibility test after:
  * the fixed-f4 Eq131 branch develops the common u~0.7013 ghost,
  * boundary-value freedom is absent at u~0.61,
  * f4X drops out of Eq131 after JA reduction,
  * the direct f4-resolved branch remains kinetically unhealthy.

No finite deformation is promoted here.  The result is only a local tangent
existence/nonexistence statement.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.numerics import module  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central  # noqa: E402
from ssz_p5.production.full_action_lower import complete_total_action_jets, emit_lower_slots  # noqa: E402
from ssz_p5.production.central_action import background_residuals  # noqa: E402

OUT = ROOT / "data/generated/spectral/CURRENT_MEMBER_EQ131_PROJECTED_CONTROLLABILITY.json"

DIRECTIONS = ("f2", "f2X", "f2F", "f3", "f3X", "f4", "f4X", "f4XX")
LS = (6, 20, 42, 110, 420, 1000)
POINTS = (0.7013, 0.7030, 0.706135)
EPS = 1e-5


def bump(u, u0, width=7e-4):
    z = (np.asarray(u, float) - float(u0)) / width
    q = np.zeros_like(z)
    m = np.abs(z) < 1
    q[m] = np.exp(-1/(1-z[m]**2)) / np.exp(-1)
    return q


def reemit(action):
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")
    d = action.copy().sort_values("x").reset_index(drop=True)
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float); h = d.h.to_numpy(float)
    X = d.X.to_numpy(float); A = d.A0prime.to_numpy(float)
    ph = d.phiprime.to_numpy(float)
    F = h*A*A/(2*f); Y = 4*X*F
    Xp = zk.dr(r, X, 1, 9, 8); Fp = zk.dr(r, F, 1, 9, 8); Yp = zk.dr(r, Y, 1, 9, 8)
    d["f2phi"] = (
        zk.dr(r, d.f2.to_numpy(float), 1, 9, 8)
        - d.f2X.to_numpy(float)*Xp
        - d.f2F.to_numpy(float)*Fp
        - d.f2Y.to_numpy(float)*Yp
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
        - zk.dr(r, out.A0prime.to_numpy(float)*out.v4.to_numpy(float)/2, 1, 9, 8)
        + out.A0prime.to_numpy(float)*out.v5.to_numpy(float)/2
    )
    return d, out


def evec(completed, i):
    bg = background_residuals(completed)
    return np.array([
        bg["E00"][i],
        bg["E11"][i],
        bg["JA"][i],
        bg["eq131_residual"][i],
    ], float)


def kmin(stream, red, L, i):
    a = red.canonical_audit(stream, int(L))
    K = np.asarray(a["K"], float)
    Ks = (K + K.swapaxes(1,2))/2
    return float(np.linalg.eigvalsh(Ks[i])[0])


def nullspace(B):
    U,S,Vt = np.linalg.svd(B, full_matrices=True)
    tol = max(B.shape)*np.finfo(float).eps*(S[0] if len(S) else 1.0)*100
    r = int(np.sum(S > tol))
    return Vt[r:].T, r, [float(x) for x in S]


def common_gain_lp(B, grads):
    n = B.shape[1]
    # row-normalize equalities
    Beq = B.copy()
    norms = np.linalg.norm(Beq, axis=1)
    keep = norms > 1e-14
    Beq = Beq[keep] / norms[keep,None]
    c = np.zeros(n+1); c[-1] = -1
    Aub=[]; bub=[]
    for L in LS:
        row=np.zeros(n+1)
        row[:n] = -grads[L]
        row[-1] = 1
        Aub.append(row); bub.append(0.0)
    Aeq = np.column_stack([Beq, np.zeros(len(Beq))]) if len(Beq) else None
    beq = np.zeros(len(Beq)) if len(Beq) else None
    sol=linprog(c,A_ub=np.asarray(Aub),b_ub=np.asarray(bub),A_eq=Aeq,b_eq=beq,
                bounds=[(-1,1)]*n+[(None,None)],method="highs")
    if not sol.success:
        return {"success":False,"message":sol.message}
    dp=sol.x[:n]
    gains={str(L):float(grads[L]@dp) for L in LS}
    return {
        "success":True,
        "common_gain":float(sol.x[-1]),
        "all_positive":bool(min(gains.values())>0),
        "gains":gains,
        "direction":{DIRECTIONS[j]:float(dp[j]) for j in range(n) if abs(dp[j])>1e-8},
    }


def point_audit(full_action, u0):
    # local guarded table
    loc = full_action[(full_action.u >= u0-0.003) & (full_action.u <= u0+0.003)].copy()
    loc = loc.sort_values("x").reset_index(drop=True)
    comp0, st0 = reemit(loc)
    u = loc.u.to_numpy(float)
    i = int(np.argmin(np.abs(u-u0)))
    shape = bump(u, u[i])
    E0 = evec(comp0, i)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    base={L:kmin(st0,red,L,i) for L in LS}

    B=np.zeros((4,len(DIRECTIONS)))
    grads={L:np.zeros(len(DIRECTIONS)) for L in LS}

    for j,name in enumerate(DIRECTIONS):
        d=loc.copy()
        d[name]=d[name].to_numpy(float)+EPS*shape
        comp,st=reemit(d)
        B[:,j]=(evec(comp,i)-E0)/EPS
        for L in LS:
            grads[L][j]=(kmin(st,red,L,i)-base[L])/EPS

    N,rankB,S = nullspace(B)
    projected={}
    for L in LS:
        g=grads[L]
        pg=N@(N.T@g) if N.size else np.zeros_like(g)
        projected[str(L)]={
            "lambda_min":base[L],
            "gradient_norm":float(np.linalg.norm(g)),
            "projected_norm":float(np.linalg.norm(pg)),
            "projected_fraction":float(np.linalg.norm(pg)/max(np.linalg.norm(g),1e-300)),
            "structurally_controllable":bool(np.linalg.norm(pg)>1e-10*max(1.0,np.linalg.norm(g))),
        }

    lp=common_gain_lp(B,grads)
    return {
        "u_requested":float(u0),
        "u_grid":float(u[i]),
        "background_at_point":{"E00":float(E0[0]),"E11":float(E0[1]),"JA":float(E0[2]),"Eq131":float(E0[3])},
        "rank_B":rankB,
        "nullity":int(len(DIRECTIONS)-rankB),
        "singular_values_B":S,
        "directions":list(DIRECTIONS),
        "per_L":projected,
        "simultaneous_all_L":lp,
    }


def main():
    build=build_onshell_central(ROOT)
    action=build.action.sort_values("x").reset_index(drop=True)
    points=[point_audit(action,u0) for u0 in POINTS]

    all_structural=all(all(v["structurally_controllable"] for v in p["per_L"].values()) for p in points)
    all_common=all(p["simultaneous_all_L"].get("all_positive",False) for p in points)

    if all_common:
        diagnosis="EQ131_COMPATIBLE_LOCAL_TANGENT_EXISTS_TO_IMPROVE_ALL_TESTED_L"
        nxt="BUILD_PREDICTOR_CORRECTOR_FINITE_EQ131_CONTINUATION"
    elif all_structural:
        diagnosis="PER_L_EQ131_COMPATIBLE_TANGENTS_EXIST_BUT_NO_COMMON_ALL_L_DIRECTION"
        nxt="TEST_MULTIPOINT_MULTI_L_CONTINUATION_OR_ADDITIONAL_ACTION_DIRECTIONS"
    else:
        diagnosis="EQ131_CONSTRAINED_TANGENT_LOSES_KINETIC_CONTROLLABILITY"
        nxt="EXPAND_ACTION_SPACE_OR_REVISIT_EQ131_COMPATIBILITY_ASSUMPTIONS"

    report={
        "scope":"local first-order Eq131-compatible action controllability; no finite member/no promotion",
        "points":points,
        "all_points_per_L_structurally_controllable":all_structural,
        "all_points_common_all_L_positive_direction":all_common,
        "diagnosis":diagnosis,
        "next_target":nxt,
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
