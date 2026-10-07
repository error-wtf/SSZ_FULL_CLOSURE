#!/usr/bin/env python3
"""ECS_DISCOVERY_V1 — Solver B: independent exterior-complex-scaling
spectral solve on the frozen (a1=-0.5, eps=-0.3) 3-channel operator.

Independent implementation path vs the Jost solver (matrix pencil on a
discretized grid vs ODE shooting + matching). Consumes the same frozen
operator NPZ.

Pencil (Zhang-Kase convention, see run_jost_discovery_v1.py):
  G Psi'' + A(w,r) Psi' + B(w,r) Psi = 0,
  A(w,r) = i w R + (G' - S),   B(w,r) = w^2 K - i w R - (M + S'/2)
with R = 0 in the exported convention (audit: time-first-odd terms fold
into S; max_time_first_sym_res ~ 0).

Discretize with Dirichlet interior BC at r_in (regularity enforced by
the center basis, P0) and outgoing behavior handled by exterior complex
scaling beyond r_match: r -> r_m + e^{i theta} (r - r_m). The quadratic
eigenvalue problem in w is linearized and solved with scipy.linalg.eig.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig as geig

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
ART = ROOT / "data/generated/spectral"
NPZ = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V2_3DOF.npz"


def deriv_arrays(r, A):
    dA = np.empty_like(A)
    dA[1:-1] = (A[2:] - A[:-2]) / (r[2:] - r[:-2])[:, None, None]
    dA[0] = dA[1]
    dA[-1] = dA[-2]
    return dA


def build_pencil(theta_deg, n_pts, r_in=0.05, r_out=60.0, r_match=30.0):
    """Build the discretized quadratic pencil blocks (C2, C1, C0) with
    complex scaling beyond r_match: C2 w^2 + C1 w + C0, acting on the
    3n-dim field vector."""
    d = np.load(NPZ)
    r_full = d["r"].astype(float)
    theta = np.deg2rad(theta_deg)
    # complex-scaled grid: linear inside, scaled outside
    r_in_s = r_in
    # interior part [r_in, r_match], scaled part r_match..r_out mapped
    n_in = int(n_pts * (r_match - r_in)/(r_out - r_in))
    n_out_ = n_pts - n_in
    r1 = np.linspace(r_in, r_match, n_in, endpoint=False)
    dr_scaled = (r_out - r_match)*np.exp(1j*theta)/n_out_
    # the SCALED coordinate s runs r_match..r_out; physical point is
    # r_m + e^{i theta}(s - r_m)
    r2 = r_match + dr_scaled*np.arange(1, n_out_+1)
    r_grid = np.concatenate([r1, r2])  # complex beyond r_match
    # interpolate operator blocks at |r| (they decay ~ const beyond the
    # window; use the asymptotic values for the scaled region)
    def blocks_at(rr):
        """K,G,S,M at real coordinate rr (clip to domain; asymptote
        constant beyond r[-1])."""
        rq = np.clip(np.real(rr), r_full[0], r_full[-1])
        idx = np.searchsorted(r_full, rq) - 1
        idx = np.clip(idx, 0, len(r_full)-2)
        t = (rq - r_full[idx])/(r_full[idx+1]-r_full[idx])
        def ler(Q):
            return Q[idx]*(1-t)[:, None, None] + Q[idx+1]*t[:, None, None]
        return ler(d["K_phys"]), ler(d["G_phys"]), ler(d["S_phys"]), ler(d["M_phys"])
    K, G, S, M = blocks_at(r_grid)
    Gp = deriv_arrays(np.real(r_grid), G)
    Sp = deriv_arrays(np.real(r_grid), S)
    n = len(r_grid)
    # derivative matrices on the complex-scaled grid:
    # d/dr_phys = (dr_grid/dr_phys)^-1 d/dr_grid; in scaled region the
    # coordinate stretch is e^{i theta}, so D1 -> e^{-i theta} D1 and
    # D2 -> e^{-2i theta} D2 there.
    h = np.diff(r_grid)
    D1 = np.zeros((n, n), complex)
    D2 = np.zeros((n, n), complex)
    for i in range(1, n-1):
        hm, hp = h[i-1], h[i]
        D1[i, i-1] = -hp/(hm*(hm+hp))
        D1[i, i+1] = hm/(hp*(hm+hp))
        D1[i, i] = (hp-hm)/(hm*hp)
        D2[i, i-1] = 2.0/(hm*(hm+hp))
        D2[i, i+1] = 2.0/(hp*(hm+hp))
        D2[i, i] = -2.0/(hm*hp)
    scale = np.where(np.real(r_grid) > r_match,
                      np.exp(-1j*theta), 1.0)
    D1 = D1 * scale[:, None]
    D2 = D2 * (scale[:, None]**2)
    # Dirichlet at r_in (row 0), smooth-continuation row at the end
    D1[0, :] = 0; D1[0, 0] = 1.0
    D2[0, :] = 0
    D1[-1, :] = 0; D1[-1, -1] = 1.0
    D2[-1, :] = 0
    # Kronecker over the 3 channels
    I3 = np.eye(3)
    def blk3(Mb):
        # Mb: (n,3,3) physical blocks -> (3n,3n)
        out = np.zeros((3*n, 3*n), complex)
        for i in range(n):
            out[3*i:3*i+3, 3*i:3*i+3] = Mb[i]
        return out
    D1_3 = np.kron(D1, I3)
    D2_3 = np.kron(D2, I3)
    G_b = blk3(G); A1 = blk3(Gp - S)          # omega-independent Psi'
    M_b = blk3(M + Sp/2.0)                    # -(B part) minus w^2K
    K_b = blk3(K)
    # pencil: (w^2 K) + w*0 + [G D2 + (G'-S) D1 - (M + S'/2)]
    C2 = K_b
    C1 = np.zeros_like(K_b)                   # R = 0 in this export
    C0 = G_b @ D2_3 + A1 @ D1_3 - M_b
    return C2, C1, C0, r_grid


def solve_spectrum(theta_deg, n_pts):
    C2, C1, C0, r_grid = build_pencil(theta_deg, n_pts)
    n3 = C2.shape[0]
    # linearize: x = [Psi; w Psi]
    # [C2 0; 0 I] [Psi'; w Psi'] ... standard companion form:
    # w^2 C2 Psi + w C1 Psi + C0 Psi = 0
    # -> [ -C0  -C1 ; 0 I ] x = w [ C2 0; 0 I ] x
    A_mat = np.zeros((2*n3, 2*n3), complex)
    B_mat = np.zeros((2*n3, 2*n3), complex)
    A_mat[:n3, :n3] = -C0
    A_mat[:n3, n3:] = -C1
    A_mat[n3:, n3:] = np.eye(n3)
    B_mat[:n3, :n3] = C2
    B_mat[n3:, n3:] = np.eye(n3)
    w_all, *_ = geig(A_mat, B_mat)
    return w_all


def main() -> int:
    t0 = time.time()
    out = {"audit": "ECS_DISCOVERY_V1",
            "note": ("independent matrix-pencil path; discovery-only, "
                      "no Jost information used"),
            "runs": []}
    # theta ladder: stability check comes later (C-R4); discovery uses
    # the canonical theta values
    all_specs = {}
    for theta in (15.0, 25.0):
        for n_pts in (700,):
            t1 = time.time()
            w = solve_spectrum(theta, n_pts)
            # keep finite eigenvalues with Im < 0 (decaying) and Re > 0
            keep = (np.abs(w) < 50) & (np.imag(w) < -1e-6) & (np.real(w) > 1e-6)
            ws = w[keep]
            ws = ws[np.argsort(np.real(ws))]
            key = f"theta{int(theta)}_n{n_pts}"
            all_specs[key] = [[round(float(z.real), 6),
                                round(float(z.imag), 6)] for z in ws[:40]]
            out["runs"].append({
                "theta_deg": theta, "n_pts": n_pts,
                "n_kept": int(keep.sum()),
                "wall_seconds": round(time.time()-t1, 1),
            })
            print(f"theta={theta} n={n_pts}: {keep.sum()} kept "
                  f"({time.time()-t1:.0f}s)", flush=True)
    out["spectra"] = all_specs
    out["wall_seconds"] = round(time.time() - t0, 1)
    (ART / "ECS_CANDIDATES_V1.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
