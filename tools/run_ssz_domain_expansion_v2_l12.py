#!/usr/bin/env python3
"""L=12 tile scans — same declared protocol as V1, new operator + 4 tiles."""
import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
ART = ROOT / "data/generated/spectral"

BRANCH = {"a1": -0.5, "eps": -0.3, "L": 12}
TILES = {
    "D0": {"re": [0.05, 1.5], "im": [-0.8, -0.05]},
    "D1": {"re": [1.5, 3.0], "im": [-0.8, -0.05]},
    "D2": {"re": [0.05, 1.5], "im": [-2.0, -0.8]},
    "D3": {"re": [1.5, 3.0], "im": [-2.0, -0.8]},
}
EXPORT = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V4_L12.npz"
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
    print(f"L12 export: {n} nodes, dim {dim}", flush=True)

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
                Km = np.zeros((3,3)); Gm = I3; Sm = np.zeros((3,3))
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
                    for i in range(2):
                        for j in range(2):
                            kv = ke[i, j]; mv = ml[i, j]; dv = de[i, j]
                            if kv != 0 and g != 0:
                                Ks[(e+i)*dim+a_, (e+j)*dim+b_] += kv * g
                            if mv != 0 and m_ab != 0:
                                MV[(e+i)*dim+a_, (e+j)*dim+b_] += mv * m_ab
                            if dv != 0 and s_ab != 0:
                                Ks[(e+i)*dim+a_, (e+j)*dim+b_] += dv * s_ab
        A = Ks + MV
        return A[:-dim, :-dim].tocsc(), ML[:-dim, :-dim].tocsc()

    h = 0.02
    A, L = build_fem(45.0, h)
    print(f"ECS pencil dim: {A.shape[0]}", flush=True)

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
        (ART / f"ECS_SSZ_L12_{tname}_CANDIDATES_V1.json").write_text(
            json.dumps({"audit": f"ECS_SSZ_L12_{tname}_CANDIDATES_V1",
                         "tile": tile, "branch": BRANCH,
                         "candidates": refined, "n": len(refined)}, indent=1) + "\n")
        print(f"[{tname}] ECS total: {len(refined)}", flush=True)

    h_j = 0.02
    xs = np.arange(r[0], r[-1], h_j)

    def coeff_at(x):
        i = min(int((x - r[0]) / h_j), n - 2)
        t = (x - r[i]) / (r[i+1] - r[i])
        def Lq(Q): return Q[i]*(1-t) + Q[i+1]*t
        return Lq(K), Lq(G), Lq(S), Lq(M), Lq(Gp), Lq(Sp)

    def jost_match_residual(w):
        Lmat = np.zeros((3, 3), dtype=complex)
        x = r[0]
        for i in range(len(xs) - 1):
            hh = xs[i+1] - xs[i]
            def f(xx, LL):
                Km, Gm, Sm, Mm, Gpm, Spm = coeff_at(xx)
                return (np.linalg.solve(Gm, (Mm + 0.5*Spm - (w*w)*Km))
                        + np.linalg.solve(Gm, Sm - Gpm) @ LL - LL @ LL)
            try:
                k1 = f(x, Lmat)
                k2 = f(x+hh/2, Lmat+hh/2*k1)
                k3 = f(x+hh/2, Lmat+hh/2*k2)
                k4 = f(x+hh, Lmat+hh*k3)
                Lmat = Lmat + hh/6*(k1+2*k2+2*k3+k4)
            except np.linalg.LinAlgError:
                return None
            x += hh
            if np.max(np.abs(Lmat)) > 1e12:
                return None
        return float(np.linalg.norm(Lmat - 1j * w * np.eye(3)))

    jost_by_tile = {}
    for tname, tile in TILES.items():
        cands = []
        for (sr, si) in sigma_lattice(tile):
            w = complex(sr, si)
            res = jost_match_residual(w)
            if res is not None and res < 1.0:
                cands.append({"omega_t": [w.real, w.imag], "matching_residual": res})
        refined = list(cands)
        seen = {(round(c["omega_t"][0], 3), round(c["omega_t"][1], 3)) for c in cands}
        for c in cands:
            wr, wi = c["omega_t"]
            wre = (tile["re"][1]-tile["re"][0])/5/2
            wim = (tile["im"][0]-tile["im"][1])/3/2
            for dr in np.linspace(-wre, wre, 3):
                for di in np.linspace(-wim, wim, 3):
                    if abs(dr) < 1e-12 and abs(di) < 1e-12:
                        continue
                    w = complex(wr + dr, wi + di)
                    res = jost_match_residual(w)
                    if res is not None and res < 1.0:
                        key = (round(w.real, 3), round(w.imag, 3))
                        if key in seen:
                            continue
                        seen.add(key)
                        refined.append({"omega_t": [w.real, w.imag],
                                         "matching_residual": res, "refined": True})
        best = {}
        for c in sorted(refined, key=lambda c: c["matching_residual"]):
            key = (round(c["omega_t"][0], 3), round(c["omega_t"][1], 3))
            if key not in best:
                best[key] = c
        jost_by_tile[tname] = list(best.values())
        (ART / f"JOST_SSZ_L12_{tname}_CANDIDATES_V1.json").write_text(
            json.dumps({"audit": f"JOST_SSZ_L12_{tname}_CANDIDATES_V1",
                         "tile": tile, "branch": BRANCH,
                         "candidates": jost_by_tile[tname],
                         "n": len(jost_by_tile[tname])}, indent=1) + "\n")
        print(f"[{tname}] Jost total: {len(jost_by_tile[tname])}", flush=True)

    matching = {"audit": "SSZ_DOMAIN_EXPANSION_V2_MATCHING_L12", "branch": BRANCH,
                 "tiles": {}}
    any_common = False
    for tname in TILES:
        ecs_set = {(round(c["omega_t"][0], 2), round(c["omega_t"][1], 2))
                   for c in ecs_by_tile[tname]}
        jost_set = {(round(c["omega_t"][0], 2), round(c["omega_t"][1], 2))
                    for c in jost_by_tile[tname]}
        common = sorted(ecs_set & jost_set)
        matching["tiles"][tname] = {"ecs_n": len(ecs_by_tile[tname]),
                                     "jost_n": len(jost_by_tile[tname]),
                                     "common": common, "common_n": len(common)}
        if common:
            any_common = True
    matching["wall_seconds"] = round(time.time() - t0, 1)
    matching["any_common"] = any_common
    (ART / "SSZ_DOMAIN_EXPANSION_V2_MATCHING_L12.json").write_text(
        json.dumps(matching, indent=1) + "\n")
    print("MATCHING:", json.dumps(matching["tiles"], indent=1))
    print("STATUS:", "COMMON_CANDIDATES_FOUND" if any_common else "NO_COMMON_CANDIDATES")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
