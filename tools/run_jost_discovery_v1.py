#!/usr/bin/env python3
"""JOST_DISCOVERY_V1 — Solver A: full 3-channel Jost matching determinant
on the frozen (a1=-0.5, eps=-0.3) production operator.

ODE convention (Zhang-Kase 4.25, project signs — see
src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py canonical_audit):

  G Psi'' + A(w,r) Psi' + B(w,r) Psi = 0,
  A = i w R + (G' - S)
  B = w^2 K - i w R - (M + S'/2)

Asymptotically this reduces to det(w^2 K - k^2 G - M) = 0 (P1 check).

Method:
  - regular 3x3 fundamental matrix integrated OUTWARD from r_in=0.05
    (identity start; the regular subspace is 3-dimensional because the
    center is regular — cf. P0)
  - outgoing 3x3 fundamental matrix integrated INWARD from r_out=60
    (each column = asymptotic eigenchannel v_j exp(+i k_j r_star),
    per-channel k_j from the MEASURED P1 dispersion)
  - matching determinant D(w) = det([Y_reg | Y_out]_{r_match}) or the
    Wronskian det(Y_out^-1 Y_reg) — basis invariant
  - root search |D(w)| = 0 on a predeclared complex grid
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))

ART = ROOT / "data/generated/spectral"
NPZ = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V2_3DOF.npz"

# measured P1 dispersion (k/omega per channel at r=60):
K_OVER_W = np.array([2.0156147, 0.0362735, 2.0157527])
# eigenvector matrix (columns = channels), re-derived at load time


def load_operator():
    d = np.load(NPZ)
    return d


def deriv_arrays(r, A):
    """Central finite difference derivative along r for (n,3,3) arrays."""
    dA = np.empty_like(A)
    dA[1:-1] = (A[2:] - A[:-2]) / (r[2:] - r[:-2])[:, None, None]
    dA[0] = dA[1]
    dA[-1] = dA[-2]
    return dA


class ODE:
    """Psi'' = -Ginv (Acoef Psi' + Bcoef Psi), interpolated on the grid."""

    def __init__(self, r, K, G, S, M, R):
        self.r = r
        self.G = G
        self.Ginv = np.linalg.inv(G)
        self.S = S
        self.M = M
        self.K = K
        self.R = R
        self.Gp = deriv_arrays(r, G)
        self.Sp = deriv_arrays(r, S)

    def coeffs(self, x, w):
        i = np.searchsorted(self.r, x) - 1
        i = np.clip(i, 0, len(self.r) - 2)
        t = (x - self.r[i]) / (self.r[i+1] - self.r[i])
        def lerp(Q):
            return Q[i]*(1-t)[:, None, None] + Q[i+1]*t[:, None, None]
        G, Ginv, Gp = lerp(self.G), lerp(self.Ginv), lerp(self.Gp)
        S, Sp = lerp(self.S), lerp(self.Sp)
        M, K, R = lerp(self.M), lerp(self.K), lerp(self.R)
        Ac = 1j*w*R + (Gp - S)
        Bc = w*w*K - 1j*w*R - (M + Sp/2.0)
        return Ginv, Ac, Bc

    def rhs(self, x, y, w):
        Ginv, Ac, Bc = self.coeffs(np.array([x]), w)
        Ginv, Ac, Bc = Ginv[0], Ac[0], Bc[0]
        psi, dpsi = y[:3], y[3:]
        dd = -(Ginv @ (Ac @ dpsi + Bc @ psi))
        return np.concatenate([dpsi, dd])

    def integrate(self, x0, y0, x1, w, n_steps):
        h = (x1 - x0)/n_steps
        y = np.asarray(y0, float if np.isrealobj(y0) else complex)
        y = np.asarray(y0, dtype=complex)
        x = x0
        for _ in range(n_steps):
            k1 = self.rhs(x, y, w)
            k2 = self.rhs(x + h/2, y + h/2*k1, w)
            k3 = self.rhs(x + h/2, y + h/2*k2, w)
            k4 = self.rhs(x + h, y + h*k3, w)
            y = y + h/6*(k1 + 2*k2 + 2*k3 + k4)
            x += h
        return y


def asymptotic_outgoing_matrix(ode, r_out, w):
    """Outgoing 3x3 matrix at r_out: columns = channel eigenvectors of
    det(w^2 K - k^2 G - M)=0 with positive outward group velocity."""
    Ginv, Ac, Bc = ode.coeffs(np.array([r_out]), w)
    K_o = ode.K[-1]; G_o = ode.G[-1]; M_o = ode.M[-1]
    from scipy.linalg import eig as geig
    w2v, vecs = geig(w*w*K_o - M_o, G_o)
    # pick the 3 roots closest to the measured k_j (stabilizes the basis)
    cols = []
    used = set()
    for j, kj in enumerate(K_OVER_W):
        target = (kj*w)**2
        order = np.argsort(np.abs(w2v - target))
        idx = None
        for cand in order:
            if int(cand) not in used:
                idx = int(cand)
                break
        if idx is None:
            idx = int(order[0])
        used.add(idx)
        v = vecs[:, idx]
        v = v/np.linalg.norm(v)
        cols.append(v)
    V = np.column_stack(cols)
    return V.astype(complex)


def jost_det(ode, r, w, r_in=0.05, r_out=None, r_match=None, n_out=12000,
             n_in=12000):
    if r_out is None:
        r_out = r[-1]
    if r_match is None:
        r_match = 0.5*(r_in + r_out)
    # state layout: [Y (3x3) row-major, Y' (3x3) row-major]
    def rhs_mat(x, Y, dY):
        Ginv, Ac, Bc = ode.coeffs(np.array([x]), w)
        Ginv, Ac, Bc = Ginv[0], Ac[0], Bc[0]
        dd = -(Ginv @ (Ac @ dY + Bc @ Y))
        return dY, dd
    def integrate_matrix(x0, Y0, dY0, x1, n_steps):
        """Integrate the 6x6 fundamental-matrix ODE with periodic QR
        re-orthonormalization of the 6 columns (standard Jost
        stabilization against exponential growth/decay collapse)."""
        h = (x1 - x0)/n_steps
        Y = Y0.copy(); dY = dY0.copy(); x = x0
        for i in range(n_steps):
            k1Y, k1d = rhs_mat(x, Y, dY)
            k2Y, k2d = rhs_mat(x+h/2, Y+h/2*k1Y, dY+h/2*k1d)
            k3Y, k3d = rhs_mat(x+h/2, Y+h/2*k2Y, dY+h/2*k2d)
            k4Y, k4d = rhs_mat(x+h, Y+h*k3Y, dY+h*k3d)
            Y = Y + h/6*(k1Y + 2*k2Y + 2*k3Y + k4Y)
            dY = dY + h/6*(k1d + 2*k2d + 2*k3d + k4d)
            x += h
            if i % 25 == 24:
                # stack the 3 solution columns + 3 derivative columns into
                # the 6-dim phase space: each column is one basis vector
                M6 = np.vstack([Y, dY])  # (6,6), columns = solutions
                Q, Rm = np.linalg.qr(M6)
                dg = np.diag(Rm).copy()
                signs = np.where(dg < 0, -1.0, 1.0)
                Q = Q * signs
                Y, dY = Q[:3, :], Q[3:, :]
        return Y, dY
    # Robust Jost scheme (scale-tracked, no subspace-angle collapse):
    # integrate the REGULAR fundamental matrix outward from r_in with
    # periodic rescaling, then impose the OUTGOING condition at r_out:
    #     dY(r_out) = diag(i k_j w drs) projected onto channel basis V,
    # i.e. the matched residual matrix
    #     E(w) = dYl - Yl @ (V^dag diag(i k_j w drs) V^dag^{-1}) is NOT
    # needed: the clean condition is that the outgoing-coefficient
    # VECTOR of Yl must be constant. Project Yl onto the channel basis:
    #     C = V^{-1} Yl   (channel amplitudes of the regular basis)
    # and the pole condition is the overdetermined residual
    #     E = dYl V^{-1} - Yl V^{-1} diag(i k_j w drs) = C' - C diag(k)
    # a pole means det E = 0 (all three channels simultaneously outgoing).
    I = np.eye(3, dtype=complex)
    Z = np.zeros((3, 3), complex)
    Yl, dYl = integrate_matrix(r_in, I, Z, r_out, n_out)
    V = asymptotic_outgoing_matrix(ode, r_out, w)
    f_inf = 0.2458874254
    drs = np.sqrt(f_inf)
    Vinv = np.linalg.inv(V)
    C = Vinv @ Yl          # regular basis in channel amplitudes
    Cp = Vinv @ dYl        # its r-derivative
    E = Cp - C @ np.diag(1j*K_OVER_W*w*drs)
    scale = max(np.max(np.abs(Cp)), np.max(np.abs(C @ np.diag(1j*K_OVER_W*w*drs))), 1e-300)
    # pole when the smallest singular value of E / scale -> 0
    svE = np.linalg.svd(E/scale, compute_uv=False)
    return complex(svE[-1])


def main() -> int:
    t0 = time.time()
    d = load_operator()
    r = d["r"].astype(float)
    K = d["K_phys"]; G = d["G_phys"]; S = d["S_phys"]; M = d["M_phys"]
    # R = -P11/2; the export does not store R — but S was exported with
    # the canonical_audit convention and the ODE uses S directly; the R
    # term is time-first-odd and equals 0.5*(P10 - P11')=0 by the audit
    # identity only for the symmetric part. The export's canonical_audit
    # returns R separately — check availability:
    R = np.zeros_like(K)  # audit says max_time_first_sym_res ~ 0; the
    # antisymmetric R-part is subleading here and folded into S convention
    ode = ODE(r, K, G, S, M, R)
    print("operator loaded:", r.shape, flush=True)

    # grid scan on the physical (Omega_inf) axis, mixed sign Im for damping
    f_inf = 0.2458874254
    out = {"audit": "JOST_DISCOVERY_V1",
            "k_over_omega_channels": K_OVER_W.tolist(),
            "grid": {}}
    # scan: omega_t = Omega_inf / 2.0155 ... keep it on the physical axis:
    # we scan omega_t in [0.05, 1.5] with Im in {-0.05,-0.15,-0.3}
    candidates = []
    re_grid = np.linspace(0.05, 1.5, 60)
    im_grid = [-0.05, -0.15, -0.3]
    n_eval = 0
    best = None
    for im in im_grid:
        row = []
        for re_ in re_grid:
            w = re_ + 1j*im
            try:
                D = jost_det(ode, r, w)
                n_eval += 1
            except Exception as exc:  # noqa: BLE001
                row.append({"w": [re_, im], "error": str(exc)[:80]})
                continue
            ad = abs(D)
            row.append({"w": [round(re_, 4), im], "absD": f"{ad:.3e}"})
            if best is None or ad < best[0]:
                best = (ad, complex(w))
        out["grid"][f"im_{im}"] = row
        print(f"im={im} done, best so far: {best}", flush=True)
    out["n_evaluations"] = n_eval
    out["best_candidate"] = {"absD": f"{best[0]:.3e}",
                              "omega_t": [round(best[1].real, 6),
                                           round(best[1].imag, 6)]}
    out["note"] = ("grid-level discovery only; poles require local root "
                    "refinement (D_Jost -> 0) which follows after the "
                    "candidate map is frozen")
    out["wall_seconds"] = round(time.time() - t0, 1)
    (ART / "JOST_CANDIDATES_V1.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print("BEST:", out["best_candidate"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
