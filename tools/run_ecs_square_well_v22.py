#!/usr/bin/env python3
"""ECS_V2_2 — Exterior Complex Scaling certification on the SAME frozen
half-line square-well root as JOST_V2_2 (structurally independent solver).

Solver family: complex-scaled Hamiltonian eigenvalue pencil  H(theta) psi =
omega^2 psi  on a single complex FD grid (interior unscaled [0, r_match],
rotated tail z = r_match + e^{i theta} (r - r_match)).  Dense/ARPACK
eigenvalues in the lambda = omega^2 plane; NO propagation, NO Wronskian,
NO shared root-finding residual with the Jost side (producer/validator
split per DUAL_SOLVER_CONTRACT_V2).

Geometry (declared BEFORE any production scan):
  |arg(omega_ref)| = atan(1.4820919693/2.2383463024) = 33.53 deg.
  The pole is representable only for theta > 33.53 deg — the historical
  ECS theta=25 failure ("diskrete Kontinuums-Kette, kein Pol") is exactly
  this sector exclusion, not a solver defect.  Production theta: 45 deg,
  stability ladder 40/45/50 deg.

FROZEN GATES (declared in this same commit, before the production run):
  G1 accuracy:      rel_err(omega) <= 0.02 at production resolution
  G2 convergence:   rel_err strictly decreasing across the h-ladder
                    (>= 2 refinements; fitted order recorded, reported as
                    "consistent with r^-p to fitted precision")
  G3 theta-stab:    spread over theta in {40,45,50} <= 2e-3 (rel)
  G4 tail-stab:     spread over tail length in {8,12,16} <= 2e-3 (rel)
  G5 negative ctrl: V = +2.5 (barrier) -> nearest eigenvalue-root at
                    rel distance > 0.1 from REF (cert must FAIL to find it)
  G6 isolation:     nearest other eigenvalue in the lambda-plane >= 0.05
                    from the accepted lambda
  G7 outgoing tail: |psi| on the last 10% of the rotated tail < 0.2 x
                    |psi| near r_match (physical outgoing-decay guard)
Any gate failure -> status ECS_V2_2_NOT_CERTIFIED with per-gate reasons
(fail-closed).  Tolerances are never edited after results exist; a failed
version must be superseded by a NEW declared version.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eig as dense_eig

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
ART = ROOT / "data/generated/spectral"
REF_FILE = ART / "HALF_LINE_SQUARE_WELL_REFERENCE_ROOTS_V2.json"

# ---- frozen problem definition (identical to the Jost V2.2 control) ----
V0, A_INT = -2.5, 1.0          # square well on (0, a), Neumann center
R_MATCH = 6.0                  # unscaled interior ends here
THETA_PROD = 45.0              # deg; sector must contain arg(omega_ref)
H_LADDER = [0.02, 0.01, 0.005]          # interior spacing h
THETA_LADDER = [40.0, 45.0, 50.0]
TAIL_LADDER = [8.0, 12.0, 16.0]        # rotated-tail length
N_TAIL_PER_H = 200             # tail nodes per unit length at h=0.005

# ---- frozen gates (see module docstring) ----
G1_TOL = 0.02
G3_TOL = 2e-3
G4_TOL = 2e-3
G5_MIN_DIST = 0.1
G6_ISO = 0.05
G7_TAIL_RATIO = 0.2


def load_ref() -> complex:
    d = json.loads(REF_FILE.read_text())
    re_, im, _res = d["N_roots"][0]
    return complex(re_, im)


def build_grid(h: float, theta_deg: float, tail_len: float):
    """Interior [0, R_MATCH] with r=a exactly on a node (V2.2-b lesson);
    rotated tail of length tail_len with |spacing| = h."""
    m_per_unit = int(round(1.0 / h))
    assert abs(m_per_unit * h - 1.0) < 1e-12
    n_in = int(round(R_MATCH / h))            # r_match/h must be integer
    assert abs(n_in * h - R_MATCH) < 1e-12
    assert A_INT * m_per_unit % 1 == 0        # a exactly on a node
    r1 = np.arange(0, n_in + 1) * h           # includes 0 and R_MATCH
    n_tail = int(round(tail_len / h))
    s = np.arange(1, n_tail + 1) * h
    th = math.radians(theta_deg)
    z2 = R_MATCH + np.exp(1j * th) * s
    z = np.concatenate([r1.astype(complex), z2])
    return z


def build_h(z: np.ndarray, v0: float, a: float):
    """H = d2/dz2-scaled operator: H psi = (-psi'' + V psi) = lambda psi,
    3-point stencil with the general complex-spacing formula; Neumann
    center (one-sided 2nd order), Dirichlet truncation at z_end."""
    n = len(z)
    main = np.zeros(n, complex)
    lower = np.zeros(n - 1, complex)
    upper = np.zeros(n - 1, complex)
    for i in range(1, n - 1):
        hl = z[i] - z[i - 1]
        hr = z[i + 1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        lower[i - 1] = c * hr
        upper[i] = c * hl
        main[i] = -c * (hl + hr)
    # center row: (-3 psi0 + 4 psi1 - psi2)/(2h) = 0
    h0 = z[1] - z[0]
    main[0] = -3.0 / (2 * h0)
    upper[0] = 4.0 / (2 * h0)
    # NOTE: the psi2 coefficient lands on upper[1]; splice it in:
    # row 0 has entries at cols 0,1,2 -> build rows via lil instead.
    rows = []
    cols = []
    vals = []
    for i in range(1, n - 1):
        hl = z[i] - z[i - 1]
        hr = z[i + 1] - z[i]
        c = 2.0 / (hl * hr * (hl + hr))
        rows += [i, i, i]
        cols += [i - 1, i, i + 1]
        vals += [c * hr, -c * (hl + hr), c * hl]
    h0 = z[1] - z[0]
    h2 = z[2] - z[0]
    rows += [0, 0, 0]
    cols += [0, 1, 2]
    vals += [-3.0 / h2, 4.0 / h2, -1.0 / h2]
    # Dirichlet at the truncated end:
    rows.append(n - 1)
    cols.append(n - 1)
    vals.append(1.0 + 0j)
    D2 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    V = np.where(np.real(z) <= a, v0, 0.0)
    V[-1] = 0.0
    # -psi'' + V psi = lambda psi
    H = (-D2) + sp.diags(V.astype(complex))
    return H.tocsc(), V


def roots_from_eigs(lams, ref_w):
    """lambda -> omega = principal sqrt (Re>0); keep Im(omega)<0 sector."""
    out = []
    for lam in lams:
        w = np.sqrt(complex(lam))
        if w.imag > 0:
            w = -w
        out.append(w)
    return out


def solve_near(H, lam_ref, k=14):
    vals = spla.eigs(H, k=k, sigma=lam_ref, which="LM",
                     return_eigenvectors=True)
    return vals


def main() -> int:
    t0 = time.time()
    REF = load_ref()
    LAM_REF = REF * REF
    out = {
        "audit": "ECS_V2_2_CERTIFICATE",
        "problem": {
            "control": "half-line square well, Neumann center, V=-2.5 on (0,1), outgoing infinity",
            "identical_frozen_inputs_as_jost_v2_2": str(REF_FILE.name),
            "reference_root": [REF.real, REF.imag],
            "lambda_ref": [LAM_REF.real, LAM_REF.imag],
            "arg_omega_ref_deg": math.degrees(math.atan2(
                -REF.imag, REF.real)),
        },
        "method": {
            "family": "exterior complex scaling, complex FD Hamiltonian eigenpencil",
            "independence": ("no propagation, no Wronskian, no shared "
                             "root-finding residual with JOST_V2_2 "
                             "(shares only frozen coefficient inputs)"),
            "sector_note": ("|arg(omega_ref)|=33.53deg -> lambda-arg=-67deg; "
                            "with the continuum rotated DOWN by e^{-2i theta} "
                            "the pole is exposed for theta < 33.5deg "
                            "(historical theta=25 was geometrically correct; "
                            "the earlier 'sector exclusion' reading was "
                            "inverted)"),
        },
        "gates_declared_before_run": {
            "G1_rel_err_max": G1_TOL,
            "G2": "rel_err strictly decreasing across h-ladder",
            "G3_theta_spread_max": G3_TOL,
            "G4_tail_spread_max": G4_TOL,
            "G5_negative_min_rel_dist": G5_MIN_DIST,
            "G6_isolation_min": G6_ISO,
            "G7_tail_ratio_max": G7_TAIL_RATIO,
        },
    }
    gates = {}

    # ---- analytic stencil self-check on the actual complex tail ----
    z = build_grid(0.02, THETA_PROD, 12.0)
    kk = 1.3 - 0.7j
    f = np.exp(1j * kk * z)
    H, Vv = build_h(z, V0, A_INT)
    res = H @ f - (kk**2 + Vv) * f
    # boundary rows (0 and end) are BC rows, exclude from the interior check
    interior_res = np.abs(res[1:-1])
    scale = np.max(np.abs((kk**2 + Vv[1:-1]) * f[1:-1]))
    out["stencil_self_check"] = {
        "max_interior_residual": float(np.max(interior_res)),
        "rel": float(np.max(interior_res) / scale),
    }
    gates["stencil_self_check"] = out["stencil_self_check"]["rel"] < 1e-4

    # ---- dense-vs-ARPACK cross-check at the coarsest grid ----
    zc = build_grid(H_LADDER[0], THETA_PROD, 12.0)
    Hc, _ = build_h(zc, V0, A_INT)
    lam_dense = dense_eig(Hc.toarray(), right=False)
    lam_dense = lam_dense[0] if isinstance(lam_dense, tuple) else lam_dense
    wd = roots_from_eigs(lam_dense, REF)
    cand_d = [w for w in wd
              if abs(w - REF) / abs(REF) < 0.2 and w.real > 0]
    lam_arpack, vecs = solve_near(Hc, LAM_REF)
    wa = roots_from_eigs(lam_arpack, REF)
    cand_a = [w for w in wa if abs(w - REF) / abs(REF) < 0.2 and w.real > 0]
    cross = None
    if cand_d and cand_a:
        bd = min(cand_d, key=lambda w: abs(w - REF))
        ba = min(cand_a, key=lambda w: abs(w - REF))
        cross = abs(bd - ba) / abs(REF)
    out["dense_vs_arpack_cross_check"] = {
        "dense_candidates": [[w.real, w.imag] for w in cand_d],
        "arpack_candidates": [[w.real, w.imag] for w in cand_a],
        "max_rel_disagreement": cross,
    }
    gates["dense_vs_arpack"] = cross is not None and cross < 1e-6

    def ecs_root(h, theta_deg, tail_len):
        z_ = build_grid(h, theta_deg, tail_len)
        H_, _ = build_h(z_, V0, A_INT)
        lam, vecs_ = solve_near(H_, LAM_REF)
        ws = roots_from_eigs(lam, REF)
        ok = [i for i, w in enumerate(ws)
              if w.real > 0 and w.imag < 0
              and abs(w - REF) / abs(REF) < 0.2]
        if not ok:
            return None, (lam, ws, None, z_)
        i_best = min(ok, key=lambda i: abs(ws[i] - REF))
        # G6 isolation in the lambda-plane:
        lam_best = lam[i_best]
        others = [lam[j] for j in range(len(lam)) if j != i_best]
        iso = float(min(abs(o - lam_best) for o in others))
        # G7 outgoing-tail decay of the eigenvector:
        v = vecs_[:, i_best]
        n_tail_zone = int(0.10 * len(z_))
        tail_amp = float(np.mean(np.abs(v[-n_tail_zone:])))
        im_zone = np.abs(np.real(z_) - R_MATCH) < 0.5
        ref_amp = float(np.mean(np.abs(v[im_zone]))) if im_zone.any() else 1.0
        ratio = tail_amp / max(ref_amp, 1e-300)
        return ws[i_best], (lam, ws, (iso, ratio), z_)

    # ---- G2 resolution ladder (production theta/tail) ----
    ladder = []
    for h in H_LADDER:
        w, aux = ecs_root(h, THETA_PROD, 12.0)
        rel = None if w is None else abs(w - REF) / abs(REF)
        ladder.append({"h": h, "omega": None if w is None else [w.real, w.imag],
                       "rel_err": rel})
    out["resolution_ladder"] = ladder
    rels = [e["rel_err"] for e in ladder]
    gates["G1"] = rels[-1] is not None and rels[-1] <= G1_TOL
    gates["G2"] = all(
        rels[i] is not None and rels[i + 1] is not None and rels[i + 1] < rels[i]
        for i in range(len(rels) - 1))
    # fitted order (reported, not gated):
    import numpy as _np
    hs = _np.array([e["h"] for e in ladder])
    es = _np.array(rels)
    if es[-1] is not None and all(e > 0 for e in es):
        p, _c = _np.polyfit(_np.log(hs), _np.log(es), 1)
        out["fitted_convergence_order"] = float(-p)

    # ---- G3 theta ladder at production h ----
    th_roots = []
    for th in THETA_LADDER:
        w, _aux = ecs_root(H_LADDER[-1], th, 12.0)
        th_roots.append(w)
    if all(w is not None for w in th_roots):
        spread = max(abs(a - b) for a in th_roots for b in th_roots) / abs(REF)
        out["theta_ladder"] = {"thetas": THETA_LADDER,
                               "omegas": [[w.real, w.imag] for w in th_roots],
                               "rel_spread": float(spread)}
        gates["G3"] = spread <= G3_TOL
    else:
        out["theta_ladder"] = {"thetas": THETA_LADDER, "omegas": None}
        gates["G3"] = False

    # ---- G4 tail ladder at production h/theta ----
    tl_roots = []
    for tl in TAIL_LADDER:
        w, _aux = ecs_root(H_LADDER[-1], THETA_PROD, tl)
        tl_roots.append(w)
    if all(w is not None for w in tl_roots):
        spread = max(abs(a - b) for a in tl_roots for b in tl_roots) / abs(REF)
        out["tail_ladder"] = {"tail_lengths": TAIL_LADDER,
                              "omegas": [[w.real, w.imag] for w in tl_roots],
                              "rel_spread": float(spread)}
        gates["G4"] = spread <= G4_TOL
    else:
        out["tail_ladder"] = None
        gates["G4"] = False

    # ---- G6/G7 from the production run ----
    w_prod, aux_prod = ecs_root(H_LADDER[-1], THETA_PROD, 12.0)
    extra_prod = aux_prod[2] if aux_prod is not None else None
    if w_prod is not None and extra_prod is not None:
        iso, ratio = extra_prod
        out["production"] = {
            "omega": [w_prod.real, w_prod.imag],
            "rel_err": abs(w_prod - REF) / abs(REF),
            "lambda_isolation": iso,
            "tail_amplitude_ratio": ratio,
        }
        gates["G6"] = iso >= G6_ISO
        gates["G7"] = ratio <= G7_TAIL_RATIO
    else:
        out["production"] = None
        gates["G6"] = False
        gates["G7"] = False

    # ---- G5 negative control: barrier instead of well ----
    zn = build_grid(H_LADDER[-1], THETA_PROD, 12.0)
    Hn, _ = build_h(zn, +2.5, A_INT)
    lam_n, _ = solve_near(Hn, LAM_REF)
    wn = roots_from_eigs(lam_n, REF)
    wn_ok = [w for w in wn if w.real > 0 and w.imag < 0]
    if wn_ok:
        dist = min(abs(w - REF) / abs(REF) for w in wn_ok)
    else:
        dist = float("inf")
    out["negative_control"] = {
        "potential": "+2.5 barrier (V0 sign flip)",
        "n_roots_in_sector": len(wn_ok),
        "min_rel_dist_to_ref": dist,
    }
    gates["G5"] = dist > G5_MIN_DIST

    out["gates"] = {k: bool(v) for k, v in gates.items()}
    all_pass = all(bool(v) for v in gates.values())
    out["status"] = "ECS_V2_2_CERTIFIED" if all_pass else "ECS_V2_2_NOT_CERTIFIED"
    out["failed_gates"] = [k for k, v in gates.items() if not bool(v)]
    out["wall_seconds"] = round(time.time() - t0, 1)

    ART.mkdir(parents=True, exist_ok=True)
    (ART / "ECS_V2_2_CERTIFICATE.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("status", "gates", "failed_gates",
                       "production", "negative_control",
                       "wall_seconds")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
