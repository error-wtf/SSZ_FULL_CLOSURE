#!/usr/bin/env python3
"""M4-CHECK: theta-stability of the 5 refined D1 ECS candidates.

Gate logic mirrors ECS_V2_2 certificate G3 (theta spread <= 0.002):
for each candidate w, re-solve the pencil at theta 40/45/50 deg with a
local shift-invert around sigma = w^2 and measure the spread of the
eigenvalue nearest w.
"""
import json
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from pathlib import Path

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
ART = ROOT / "data/generated/spectral"
OUT = Path("/root/.hermes/cache/scratch/m4_rerun_ecs")

BRANCH = {"a1": -0.5, "eps": -0.3, "L": 6}
EXPORT = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
BASIS = ART / "V4_CENTER_REGULAR_BASIS_V2.npz"


def main():
    d3 = np.load(EXPORT)
    cb = np.load(BASIS)
    r = d3["r"]
    K = d3["K_phys"]; G = d3["G_phys"]; S = d3["S_phys"]; M = d3["M_phys"]

    def der(A):
        dA = np.empty_like(A)
        dA[1:-1] = (A[2:] - A[:-2]) / (r[2:] - r[:-2])[:, None, None]
        dA[0] = dA[1]; dA[-1] = dA[-2]
        return dA

    Gp = der(G); Sp = der(S)
    dim = 3

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
                i0 = min(int((x0 - r0) / h), len(r) - 2)
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
                    k_ab = Km[a_, b_]
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

    def w_of(lam):
        w = np.sqrt(lam)
        return -w if w.imag > 0 else w

    cands = json.load(open(OUT / "ECS_SSZ_D1_CANDIDATES_V1_FIXED_ML.json"))
    refined = [c for c in cands["candidates"] if c.get("refined")]
    report = []
    h = 0.02
    pencils = {}
    for th in (40.0, 45.0, 50.0):
        pencils[th] = build_fem(th, h)
        print(f"theta {th}: pencil built", flush=True)
    for c in refined:
        w = complex(*c["omega_t"])
        sigma = w * w
        found = {}
        for th, (A, L) in pencils.items():
            try:
                lu = spla.splu((A - sigma * L).tocsc())
                OP = spla.LinearOperator((A.shape[0],)*2,
                                         matvec=lambda x: lu.solve(L @ x),
                                         dtype=complex)
                mu, _ = spla.eigs(OP, k=20, which="LM", return_eigenvectors=True)
            except Exception as ex:
                print(f"  w={w:.4f} theta={th}: solve failed ({ex})", flush=True)
                continue
            lams = sigma + 1.0 / np.array([complex(m) for m in mu])
            ws = np.array([w_of(l) for l in lams])
            j = int(np.argmin(np.abs(ws - w)))
            if abs(ws[j] - w) < 0.02:
                found[th] = [float(ws[j].real), float(ws[j].imag)]
        if len(found) == 3:
            arr = np.array([found[40], found[45], found[50]])
            spread = float(np.abs(arr - arr.mean(axis=0)).max() / np.abs(arr).max())
            report.append({"w45": c["omega_t"], "found": found,
                           "rel_spread_max": spread, "stable": spread <= 0.002})
            print(f"w={w:.4f}{w.imag:+.4f}j spread={spread:.2e} "
                  f"stable={spread <= 0.002}", flush=True)
        else:
            report.append({"w45": c["omega_t"], "found": found,
                           "stable": False, "note": "not tracked at all thetas"})
            print(f"w={w:.4f}{w.imag:+.4f}j NOT tracked: {sorted(found)}", flush=True)
    (OUT / "M4_THETA_STABILITY_D1.json").write_text(
        json.dumps({"audit": "M4_THETA_STABILITY_D1", "results": report}, indent=1) + "\n")
    n_stable = sum(1 for x in report if x.get("stable"))
    print(f"SUMMARY: {n_stable}/{len(report)} theta-stable", flush=True)


if __name__ == "__main__":
    main()
