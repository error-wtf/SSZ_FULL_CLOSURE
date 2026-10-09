#!/usr/bin/env python3
"""M4-RERUN A: corrected ECS tile scan D1/D2/D3 (SSZ_DOMAIN_EXPANSION_V1).

EXACT copy of tools/run_ssz_domain_expansion_v1.py build_fem + ECS path.
ONLY change (the certified M3 bug): ML now accumulates the P1 mass
matrix of K_phys (the correct w^2 weight per Jost B0 = w^2 K - ... and
ECS-discovery C2 = blk(Kc)). Everything else identical: same grid h=0.02,
theta=45deg, same 6x4 sigma lattice + one 3x3 half-spacing refinement,
same residual gate 1e-6, same tiles.
Outputs go to SCRATCH, not to the certified artifact dir.
"""
import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
ART = ROOT / "data/generated/spectral"
OUT = Path("/root/.hermes/cache/scratch/m4_rerun_ecs")
OUT.mkdir(parents=True, exist_ok=True)

BRANCH = {"a1": -0.5, "eps": -0.3, "L": 6}
TILES = {
    "D1": {"re": [1.5, 3.0], "im": [-0.8, -0.05]},
    "D2": {"re": [0.05, 1.5], "im": [-2.0, -0.8]},
    "D3": {"re": [1.5, 3.0], "im": [-2.0, -0.8]},
}
EXPORT = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
BASIS = ART / "V4_CENTER_REGULAR_BASIS_V2.npz"


def w_of(lam):
    w = np.sqrt(lam)
    return -w if w.imag > 0 else w


def in_tile(w, tile):
    return (tile["re"][0] <= w.real <= tile["re"][1]
            and tile["im"][0] <= w.imag <= tile["im"][1])


def main():
    t0 = time.time()
    d3 = np.load(EXPORT)
    cb = np.load(BASIS)
    r = d3["r"]
    K = d3["K_phys"]; G = d3["G_phys"]; S = d3["S_phys"]; M = d3["M_phys"]
    n = len(r); dim = 3
    print(f"SSZ export: {n} nodes, dim {dim}", flush=True)

    def der(A):
        dA = np.empty_like(A)
        dA[1:-1] = (A[2:] - A[:-2]) / (r[2:] - r[:-2])[:, None, None]
        dA[0] = dA[1]; dA[-1] = dA[-2]
        return dA

    Gp = der(G); Sp = der(S)

    def build_fem(theta_deg, h):
        th = math.radians(theta_deg)
        r_max = float(r[-1]); r0 = float(r[0])
        n_in = int(round((r_max - r0) / h))
        n_tail = int(round(12.0 / h))
        ne = n_in + n_tail; nn = ne + 1
        N = nn * dim
        Ks = sp.lil_matrix((N, N), dtype=complex)
        ML = sp.lil_matrix((N, N), dtype=complex)
        MV = sp.lil_matrix((N, N), dtype=complex)
        I3 = np.eye(3)
        for e in range(ne):
            if e < n_in:
                x0 = r0 + e * h
                dz = h
                i0 = min(int((x0 - r0) / h), n - 2)
                t = (x0 - r[i0]) / (r[i0+1] - r[i0])
                def L(Q): return Q[i0]*(1-t) + Q[i0+1]*t
                Km = L(K); Gm = L(G); Sm = L(S); Mm = L(M)
                Gpm = L(Gp); Spm = L(Sp)
            else:
                Km = np.zeros((3,3))
                Gm = I3; Sm = np.zeros((3,3))
                Mm = np.zeros((3,3)); Gpm = np.zeros((3,3)); Spm = np.zeros((3,3))
                dz = h * np.exp(1j * th)
            ke = np.array([[1., -1.], [-1., 1.]]) / dz
            ml = dz * np.array([[1/3, 1/6], [1/6, 1/3]])
            de = np.array([[-0.5, 0.5], [-0.5, 0.5]])
            for a_ in range(dim):
                for b_ in range(dim):
                    g = Gm[a_, b_]
                    m_ab = Mm[a_, b_] + 0.5 * Spm[a_, b_]
                    s_ab = Sm[a_, b_] - Gpm[a_, b_]
                    k_ab = Km[a_, b_]          # <<< FIX: K-assembly for ML
                    for i in range(2):
                        for j in range(2):
                            kv = ke[i, j]; mv = ml[i, j]; dv = de[i, j]
                            if kv != 0 and g != 0:
                                Ks[(e+i)*dim+a_, (e+j)*dim+b_] += kv * g
                            if mv != 0 and m_ab != 0:
                                MV[(e+i)*dim+a_, (e+j)*dim+b_] += mv * m_ab
                            if mv != 0 and k_ab != 0:
                                ML[(e+i)*dim+a_, (e+j)*dim+b_] += mv * k_ab
                            if dv != 0 and s_ab != 0:
                                Ks[(e+i)*dim+a_, (e+j)*dim+b_] += dv * s_ab
        A = Ks + MV
        return A[:-dim, :-dim].tocsc(), ML[:-dim, :-dim].tocsc()

    h = 0.02
    A, L = build_fem(45.0, h)
    ml_nnz = L.nnz
    print(f"ECS pencil dim: {A.shape[0]}; FIXED ML nnz: {ml_nnz} "
          f"(production bug: 0)", flush=True)
    # sanity: the pencil A - w^2 L must be nontrivial now
    assert ml_nnz > 0, "ML still empty - fix did not apply"

    def sigma_lattice(tile, nre=6, nim=4):
        box_re = np.linspace(tile["re"][0] + (tile["re"][1]-tile["re"][0])*0.05,
                             tile["re"][1] - (tile["re"][1]-tile["re"][0])*0.05, nre)
        box_im = np.linspace(tile["im"][1] - (tile["im"][1]-tile["im"][0])*0.08,
                             tile["im"][0] + (tile["im"][1]-tile["im"][0])*0.08, nim)
        return [(float(sr), float(si)) for sr in box_re for si in box_im]

    ecs_by_tile = {}
    for tname, tile in TILES.items():
        cands = []
        seen = set()
        for (sr, si) in sigma_lattice(tile):
            sigma = complex(sr, si)
            try:
                lu = spla.splu((A - sigma * L).tocsc())
            except RuntimeError:
                continue
            OP = spla.LinearOperator((A.shape[0],)*2,
                                     matvec=lambda x: lu.solve(L @ x),
                                     dtype=complex)
            try:
                mu, vecs = spla.eigs(OP, k=20, which="LM", return_eigenvectors=True)
            except Exception:
                continue
            for j in range(mu.size):
                lam = sigma + 1.0 / complex(mu[j])
                w = w_of(lam)
                v = vecs[:, j]
                res = float(np.linalg.norm(A @ v - lam * (L @ v))
                            / max(np.linalg.norm(L @ v), 1e-300))
                if in_tile(w, tile) and res <= 1e-6:
                    key = (round(w.real, 3), round(w.imag, 3))
                    if key in seen:
                        continue
                    seen.add(key)
                    cands.append({"omega_t": [w.real, w.imag], "pencil_resid": res,
                                  "sigma_window": [sr, si]})
                    print(f"  [{tname}] ECS cand: w={w.real:.5f}{w.imag:+.5f}j res={res:.1e}", flush=True)
        refined = list(cands)
        seen2 = {tuple(c["omega_t"]) for c in cands}
        for c in cands:
            wr, wi = c["omega_t"]
            wre = (tile["re"][1]-tile["re"][0])/5/2
            wim = (tile["im"][0]-tile["im"][1])/3/2
            for dr in np.linspace(-wre, wre, 3):
                for di in np.linspace(-wim, wim, 3):
                    if abs(dr) < 1e-12 and abs(di) < 1e-12:
                        continue
                    sigma = complex(wr + dr, wi + di)
                    try:
                        lu = spla.splu((A - sigma * L).tocsc())
                        OP = spla.LinearOperator((A.shape[0],)*2,
                                                 matvec=lambda x: lu.solve(L @ x),
                                                 dtype=complex)
                        mu, vecs = spla.eigs(OP, k=20, which="LM", return_eigenvectors=True)
                    except Exception:
                        continue
                    for j in range(mu.size):
                        lam = sigma + 1.0 / complex(mu[j])
                        w = w_of(lam)
                        v = vecs[:, j]
                        res = float(np.linalg.norm(A @ v - lam * (L @ v))
                                    / max(np.linalg.norm(L @ v), 1e-300))
                        if in_tile(w, tile) and res <= 1e-6:
                            key = (round(w.real, 4), round(w.imag, 4))
                            if key in seen2:
                                continue
                            seen2.add(key)
                            refined.append({"omega_t": [w.real, w.imag],
                                            "pencil_resid": res,
                                            "sigma_window": [sigma.real, sigma.imag],
                                            "refined": True})
                            print(f"  [{tname}] ECS refined: w={w.real:.5f}{w.imag:+.5f}j res={res:.1e}", flush=True)
        ecs_by_tile[tname] = refined
        (OUT / f"ECS_SSZ_{tname}_CANDIDATES_V1_FIXED_ML.json").write_text(
            json.dumps({"audit": f"ECS_SSZ_{tname}_CANDIDATES_V1_FIXED_ML",
                        "tile": tile, "branch": BRANCH,
                        "fix": "ML = P1 mass matrix of K_phys (was: never filled)",
                        "candidates": refined, "n": len(refined)}, indent=1) + "\n")
        print(f"[{tname}] ECS total: {len(refined)}", flush=True)

    summary = {"audit": "M4_RERUN_ECS_D1D3_FIXED_ML",
               "wall_s": round(time.time() - t0, 1),
               "tiles": {t: len(v) for t, v in ecs_by_tile.items()}}
    (OUT / "M4_RERUN_ECS_SUMMARY.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("SUMMARY:", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
