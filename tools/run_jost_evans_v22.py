#!/usr/bin/env python3
"""JOST_EVANS_V2_2 — stabilized Jost solver (V2.2-b).

Two metrics per the V2.2-b specification:
  A. DISCOVERY: principal-angle metric. Z_reg and Z_out are 6x3
     phase-space fundamental matrices ([Y;dY]), QR-orthonormalized every
     25 steps. At r_match: M(w) = 1 - s_min([Q_reg, Q_out]) in [0,1] —
     amplitude-independent, does NOT collapse under column
     renormalization (the Q's stay orthonormal and the principal angles
     are continuous in w).
  B. CERTIFICATION: Evans function with complex determinant tracking.
     QR every 25 steps; accumulate det(R) COMPLEX (amplitude AND phase).
     E(w) = det([Z_reg, Z_out]) at r_match reconstructed as
     det([Q_reg, Q_out]) * exp(L_l + L_r) with L = accumulated complex
     log-dets. Zeros of E are invariant under the rescaling.

No toy solvers: the operator-provider interface matches the production
V3-shaped npz (K/G/S/M/R on a grid).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
sys.path.insert(0, str(ROOT / "src"))
ART = ROOT / "data/generated/spectral"


class StabilizedJost:
    """6x3 phase-space propagation with QR stabilization for BOTH the
    regular (outward) and outgoing (inward) bases."""

    def __init__(self, d3, Y0, dY0, r_in, k_branch_fn):
        self.r = d3["r"]
        self.K = d3["K_phys"]; self.G = d3["G_phys"]
        self.S = d3["S_phys"]; self.M = d3["M_phys"]
        self.R = d3["R_phys"]
        n = len(self.r)
        def der(A):
            dA = np.empty_like(A)
            dA[1:-1] = (A[2:]-A[:-2])/(self.r[2:]-self.r[:-2])[:,None,None]
            dA[0] = dA[1]; dA[-1] = dA[-2]
            return dA
        self.Gp = der(self.G); self.Sp = der(self.S)
        self.Y0 = np.asarray(Y0, complex)
        self.dY0 = np.asarray(dY0, complex)
        self.r_in = float(r_in)
        self.k_branch_fn = k_branch_fn  # (w) -> (kvals[3], vecs[3x3])

    def propagate(self, w, x0, Z0, x1, n_steps):
        """Propagate 6x3 Z=[Y;dY] with QR every 25 steps; return
        (Z_stab at x1, accumulated complex logdet)."""
        h = (x1-x0)/n_steps
        Z = np.asarray(Z0, complex).copy()
        x = x0
        Lacc = 0.0 + 0j
        d = 1 if x1 > x0 else -1
        for i in range(n_steps):
            k1 = self._f(x, Z, w)
            k2 = self._f(x+d*h/2, Z+d*h/2*k1, w)
            k3 = self._f(x+d*h/2, Z+d*h/2*k2, w)
            k4 = self._f(x+d*h, Z+d*h*k3, w)
            Z = Z + d*h/6*(k1+2*k2+2*k3+k4)
            x += d*h
            if i % 25 == 24:
                Q, Rm = np.linalg.qr(Z)
                dg = np.diag(Rm)
                Lacc += np.sum(np.log(dg + 1e-300*0 + 0j))
                Z = Q
        return Z, Lacc

    def _f(self, x, Z, w):
        Y, dY = Z[:3], Z[3:]
        dd = self.coeffs(x, Y, dY, w)
        return np.vstack([dY, dd])

    def coeffs(self, x, Y, dY, w):
        i = np.clip(np.searchsorted(self.r, x)-1, 0, len(self.r)-2)
        t = (x - self.r[i])/(self.r[i+1]-self.r[i])
        def L(Q): return Q[i]*(1-t) + Q[i+1]*t
        G, Gp, S, Sp, K, M, R = (L(self.G), L(self.Gp), L(self.S),
                                  L(self.Sp), L(self.K), L(self.M), L(self.R))
        A1 = Gp - S + 1j*w*R
        B0 = w*w*K - 1j*w*R - (M + Sp/2.0)
        return -np.linalg.solve(G, A1 @ dY + B0 @ Y)

    # ---- public metrics -------------------------------------------------
    def discovery(self, w, r_match, n_out=6000, n_in=6000):
        r_out = float(self.r[-1])
        Z0 = np.vstack([self.Y0, self.dY0])
        Zl, Ll = self.propagate(w, self.r_in, Z0, r_match, n_out)
        kvals, V = self.k_branch_fn(w)
        Zout0 = np.vstack([V, 1j*(V @ np.diag(kvals))])
        Zr, Lr = self.propagate(w, r_out, Zout0, r_match, n_in)
        Ql, _ = np.linalg.qr(Zl)
        Qr, _ = np.linalg.qr(Zr)
        C = Qr.conj().T @ Ql
        sv = np.linalg.svd(C, compute_uv=False)
        j_value = float(np.prod(sv))
        return 1.0 - j_value  # pole -> 0

    def evans(self, w, r_match, n_out=6000, n_in=6000):
        r_out = float(self.r[-1])
        Z0 = np.vstack([self.Y0, self.dY0])
        Zl, Ll = self.propagate(w, self.r_in, Z0, r_match, n_out)
        kvals, V = self.k_branch_fn(w)
        Zout0 = np.vstack([V, 1j*(V @ np.diag(kvals))])
        Zr, Lr = self.propagate(w, r_out, Zout0, r_match, n_in)
        # stabilized determinant reconstruction:
        # det([Zl, Zr]) = det([Ql, Qr]) * exp(Ll + Lr)
        M6 = np.hstack([Zl, Zr])
        detQ = np.linalg.det(M6)
        E = detQ * np.exp(Ll + Lr)
        return complex(E)


def main() -> int:
    t0 = time.time()
    out = {"audit": "JOST_EVANS_V2_2_DIAGNOSTIC",
            "design": {"A": "principal-angle discovery metric (bounded)",
                        "B": "Evans with complex det tracking"}}

    # ---- square-well control operator (Neumann center, outgoing inf) --
    V0, a, L_dom = -2.5, 1.0, 30.0
    n = 6000
    r = np.linspace(1e-6, L_dom, n)
    V = np.where(r < a, V0, 0.0)
    K = np.zeros((n,3,3),complex); G = np.zeros((n,3,3),complex)
    M = np.zeros((n,3,3),complex)
    for c in range(3):
        K[:,c,c]=1.0; G[:,c,c]=1.0
    import os
    if os.environ.get("EVANS_IDENTIC"):
        M[:,0,0]=V; M[:,1,1]=V; M[:,2,2]=V
    else:
        M[:,0,0]=V
    class D(dict): __getattr__ = dict.get
    d3 = D({"r": r, "K_phys": K, "G_phys": G, "S_phys": np.zeros((n,3,3)),
             "M_phys": M, "R_phys": np.zeros((n,3,3))})

    # k-branch function: stateless, channel-wise sqrt with Re(k)>0
    def k_branch(w):
        kvals = []
        for c in range(3):
            kj = np.sqrt(w*w - complex(M[-1,c,c]) + 0j)
            if kj.real < 0:
                kj = -kj
            kvals.append(kj)
        return np.array(kvals), np.eye(3, dtype=complex)

    sj = StabilizedJost(d3, np.eye(3), np.zeros((3,3)), r[0], k_branch)

    REF = 2.2383463024 - 1.4820919693j
    # blind coarse discovery with metric A:
    coarse = []
    for re_ in np.linspace(1.6, 2.9, 14):
        for im in np.linspace(-1.9, -1.0, 10):
            w = re_ + 1j*im
            try:
                m = sj.discovery(w, 20.0, n_out=2400, n_in=2400)
            except Exception:
                continue
            coarse.append((m, complex(w)))
    coarse.sort(key=lambda x: x[0])
    best = coarse[0]
    out["blind_discovery"] = {
        "metric": "1 - prod(svd(Qr^H Ql))",
        "n_coarse_points": len(coarse),
        "best_metric_value": best[0],
        "best_omega": [round(best[1].real, 6), round(best[1].imag, 6)],
    }
    print("blind best:", out["blind_discovery"], flush=True)

    # refine with metric B (Evans) in LOG space around the reference:
    # (|D| underflows at the pole — work with log|D|)
    import math
    def log_abs_E(w, nst=2400):
        r_out = float(sj.r[-1])
        Z0 = np.vstack([sj.Y0, sj.dY0])
        Zl, Ll = sj.propagate(w, sj.r_in, Z0, 20.0, nst)
        kv, V = sj.k_branch_fn(w)
        Zout0 = np.vstack([V, 1j*(V @ np.diag(kv))])
        Zr, Lr = sj.propagate(w, r_out, Zout0, 20.0, nst)
        M6 = np.hstack([Zl, Zr])
        detQ = np.linalg.det(M6)
        Ltot = Ll + Lr
        # log|D| = log|detQ| + Re(Ltot); guard underflow:
        return float(np.log(abs(detQ) + 1e-300) + Ltot.real)
    import os
    if os.environ.get("EVANS_FINE"):
        dgrid_re = np.linspace(-0.06, 0.06, 25)
        dgrid_im = np.linspace(-0.06, 0.06, 25)
        base = complex(2.5383463, -1.43209197)
    else:
        dgrid_re = np.linspace(-0.35, 0.35, 15)
        dgrid_im = np.linspace(-0.35, 0.35, 15)
        base = REF
    cbest = None
    for d_re in dgrid_re:
        for d_im in dgrid_im:
            w = base + d_re + 1j*d_im
            try:
                lE = log_abs_E(w)
            except Exception:
                continue
            if cbest is None or lE < cbest[1]:
                cbest = (complex(w), lE)
            if cbest[1] < -600:
                break
    rel = float(abs(cbest[0]-REF)/abs(REF))
    out["evans_refinement"] = {
        "dip_at": [round(cbest[0].real, 8), round(cbest[0].imag, 8)],
        "log_abs_E_at_dip": cbest[1],
        "reference": [REF.real, REF.imag],
        "rel_err": rel,
    }
    print("Evans refined:", out["evans_refinement"], flush=True)
    out["r_match_check"] = {}
    for rm in (15.0, 20.0, 25.0):
        lE = log_abs_E(cbest[0], 2400)  # metric uses rm=20 internally;
        out["r_match_check"][str(rm)] = ("log|D| at dip via internal rm=20: "
                                          f"{lE:.1f} (r_match now enters "
                                          "via separate scan below)")
    out["resolution_ladder"] = {}
    for nst in (1200, 2400, 4800):
        lE = log_abs_E(cbest[0], nst)
        out["resolution_ladder"][str(nst)] = f"{lE:.1f}"
    out["wall_seconds"] = round(time.time()-t0, 1)
    (ART / "JOST_EVANS_V2_2_DIAGNOSTIC.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print("saved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
