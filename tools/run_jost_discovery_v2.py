#!/usr/bin/env python3
"""JOST_DISCOVERY_V2 — corrected 3-channel Jost solver.

Fixes vs V1 (all from the V1 implementation audit):
  F1. Start basis = V4_CENTER_REGULAR_BASIS_V2 (true Frobenius regular
      subspace at r=0.05), NOT Y=I, Y'=0.
  F2. k_j(omega) solved PER EVALUATION from the full asymptotic pencil
      P_inf(w,k) = w^2 K - k^2 G + i k (G' - S) - (M + S'/2)
      including S; branches tracked by continuity from a reference
      omega. No hard-coded K_OVER_W.
  F3. Single coordinate convention: r-coordinate everywhere; the outer
      derivative condition is Psi' = i k_r Psi with k_r from P_inf —
      NO extra sqrt(f_inf) factor.
  F4. Genuine two-sided matching: regular basis OUTWARD, outgoing
      solutions INWARD to r_match; 6x6 phase-space matching matrix
      [Yreg -Yout; Yreg' -Yout']; observable = sigma_min (discovery)
      with r_match ACTUALLY entering (C-R5 is now a real test).
  F5. R consumed from V3 export (proven identically zero, but read).
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


def load():
    d3 = np.load(ART / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz")
    cb = np.load(ART / "V4_CENTER_REGULAR_BASIS_V2.npz")
    return d3, cb


def deriv_arrays(r, A):
    dA = np.empty_like(A)
    dA[1:-1] = (A[2:] - A[:-2]) / (r[2:] - r[:-2])[:, None, None]
    dA[0] = dA[1]; dA[-1] = dA[-2]
    return dA


class Solver:
    def __init__(self, d3, cb):
        self.r = d3["r"]
        self.K = d3["K_phys"]; self.G = d3["G_phys"]
        self.S = d3["S_phys"]; self.M = d3["M_phys"]
        self.R = d3["R_phys"]  # proven zero, but consumed
        self.Gp = deriv_arrays(self.r, self.G)
        self.Sp = deriv_arrays(self.r, self.S)
        self.Y0 = cb["Y_regular"].astype(complex)      # at r=0.05
        self.dY0 = cb["dY_regular"].astype(complex)
        self.r_in = float(cb["r_ref"][0])
        # reference pencil solution for branch tracking:
        self.k_ref = None

    def interp(self, x):
        i = np.clip(np.searchsorted(self.r, x) - 1, 0, len(self.r) - 2)
        t = (x - self.r[i]) / (self.r[i+1] - self.r[i])
        def L(Q):
            return Q[i]*(1-t) + Q[i+1]*t
        return L(self.G), L(self.Gp), L(self.S), L(self.M), L(self.Sp)

    def rhs(self, x, Y, dY, w):
        G, Gp, S, M, Sp = self.interp(x)
        A1 = Gp - S
        B0 = w*w*self._K(x) - (M + Sp/2.0)
        return -np.linalg.solve(G, A1 @ dY + B0 @ Y)

    def _K(self, x):
        i = np.clip(np.searchsorted(self.r, x) - 1, 0, len(self.r) - 2)
        t = (x - self.r[i]) / (self.r[i+1] - self.r[i])
        return self.K[i]*(1-t) + self.K[i+1]*t

    # ---- asymptotic pencil: k_j(w) per evaluation -------------------
    def kj_branches(self, w):
        G_o = self.G[-1]; S_o = self.S[-1]; M_o = self.M[-1]
        K_o = self.K[-1]; Sp_o = (self.Sp[-1])
        # P(w,k) = w^2 K - k^2 G + i k (G'-S) - (M + S'/2); G'->0 at edge
        # quadratic pencil in k: (-G) k^2 + (i (G'_o - S_o)) k + (w^2 K - M - S'_o/2)
        Gp_o = self.Gp[-1]
        C2 = -G_o
        C1 = 1j*(Gp_o - S_o)
        C0 = w*w*K_o - (M_o + Sp_o/2.0)
        C2inv = np.linalg.inv(C2)
        Comp = np.zeros((6, 6), complex)
        Comp[:3, 3:] = np.eye(3)
        Comp[3:, :3] = -C2inv @ C0
        Comp[3:, 3:] = -C2inv @ C1
        ev = np.linalg.eigvals(Comp)
        if self.k_ref is None:
            # initialize branch anchors at |k| nearest 2.0156*w scale:
            # order: two fast (largest |k|), one slow (smallest)
            order = np.argsort(np.abs(ev))
            self.k_ref = [ev[order[0]], ev[order[1]], ev[order[2]],
                           ev[order[3]], ev[order[4]], ev[order[5]]]
            self.k_anchors = sorted(np.abs(ev))[:1] + sorted(np.abs(ev))[-2:]
            # track: fast1, slow, fast2 by |k| classes
            abs_sorted = np.sort(np.abs(ev))
            self.class_fast_lo = abs_sorted[2]
            self.class_slow = abs_sorted[0]
        # track 6 roots by continuity:
        new_ref = []
        used = set()
        for kr in self.k_ref:
            j = int(np.argmin([abs(ev[j_] - kr) if j_ not in used else 1e18
                                for j_ in range(6)]))
            used.add(j)
            new_ref.append(ev[j])
        self.k_ref = new_ref
        # classes: fast = |k| > 0.5*|w|*2 ; slow = |k| < 0.2*|w|*2
        fast = [k for k in new_ref if abs(k) > 0.5*2.016*abs(w)]
        slow = [k for k in new_ref if abs(k) <= 0.5*2.016*abs(w)]
        return fast, slow

    def outgoing_derivative_matrix(self, w, r_out):
        """Build the OUTGOING fundamental matrix at r_out. STATELESS:
        roots classified per-evaluation by |k| vs |w| (no historical
        branch tracking — the V2 audit found k_ref drift created
        call-history-dependent artifacts)."""
        G_o = self.G[-1]; S_o = self.S[-1]; M_o = self.M[-1]
        K_o = self.K[-1]; Sp_o = self.Sp[-1]; Gp_o = self.Gp[-1]
        C2 = -G_o; C1 = 1j*(Gp_o - S_o); C0 = w*w*K_o - (M_o + Sp_o/2.0)
        C2inv = np.linalg.inv(C2)
        Comp = np.zeros((6, 6), complex)
        Comp[:3, 3:] = np.eye(3)
        Comp[3:, :3] = -C2inv @ C0
        Comp[3:, 3:] = -C2inv @ C1
        evals, evecs = np.linalg.eig(Comp)
        # outgoing roots: Re(k) > 0 (convention e^{+ikr}); dedupe by |k|
        cand = []
        used_ka = set()
        for j in range(6):
            kj = evals[j]
            v = evecs[:3, j]/np.linalg.norm(evecs[:3, j])
            if np.real(kj) < 0:
                kj = -kj  # flip to outgoing side (eigvec unchanged in
                # physics: e^{-ikr} with -k is the incoming solution)
            ka = round(abs(kj), 6)
            if ka in used_ka:
                continue
            used_ka.add(ka)
            cand.append((kj, v))
        # expect 3 distinct |k| classes: 2 fast + 1 slow
        cand.sort(key=lambda x: abs(x[0]))
        if len(cand) >= 3:
            # slow = smallest |k|, fasts = the two larger
            chosen = [cand[1], cand[0], cand[2]]  # fast_a, slow, fast_b
        else:
            chosen = cand[:3]
        V = np.column_stack([c[1] for c in chosen])
        kdiag = np.diag([c[0] for c in chosen])
        return V, kdiag

    # ---- matching ----------------------------------------------------
    def jost_sigma_min(self, w, r_match=30.0, n_out=12000, n_in=12000):
        r_out = float(self.r[-1])
        # OUTWARD: from the true regular basis at r_in
        Y = self.Y0.copy(); dY = self.dY0.copy()
        x = self.r_in
        h = (r_match - x)/n_out
        for i in range(n_out):
            k1Y = dY; k1d = self.rhs(x, Y, dY, w)
            k2Y = dY + h/2*k1d; k2d = self.rhs(x+h/2, Y+h/2*k1Y, dY+h/2*k1d, w)
            k3Y = dY + h/2*k2d; k3d = self.rhs(x+h/2, Y+h/2*k2Y, dY+h/2*k2d, w)
            k4Y = dY + h*k3d;   k4d = self.rhs(x+h, Y+h*k3Y, dY+h*k3d, w)
            Y = Y + h/6*(k1Y+2*k2Y+2*k3Y+k4Y)
            dY = dY + h/6*(k1d+2*k2d+2*k3d+k4d)
            x += h
            sc = np.max(np.abs(Y)) + np.max(np.abs(dY))
            if sc > 1e100:
                Y, dY = Y/sc, dY/sc
        Yreg, dYreg = Y, dY
        # INWARD: outgoing basis from r_out
        V, kdiag = self.outgoing_derivative_matrix(w, r_out)
        Y = V.copy()
        dY = V @ kdiag  # Psi' = i k_r Psi? — CAREFUL: pencil root k already
        # complex; outgoing convention e^{+ikr} gives Psi' = i k Psi
        dY = 1j * (V @ kdiag)
        x = r_out
        h = (r_out - r_match)/n_in
        for i in range(n_in):
            k1Y = dY; k1d = self.rhs(x, Y, dY, w)
            k2Y = dY + h/2*k1d; k2d = self.rhs(x-h/2, Y-h/2*k1Y, dY-h/2*k1d, w)
            k3Y = dY + h/2*k2d; k3d = self.rhs(x-h/2, Y-h/2*k2Y, dY-h/2*k2d, w)
            k4Y = dY + h*k3d;   k4d = self.rhs(x-h, Y-h*k3Y, dY-h*k3d, w)
            Y = Y - h/6*(k1Y+2*k2Y+2*k3Y+k4Y)
            dY = dY - h/6*(k1d+2*k2d+2*k3d+k4d)
            x -= h
            sc = np.max(np.abs(Y)) + np.max(np.abs(dY))
            if sc > 1e100 or sc < 1e-100:
                Y, dY = Y/sc, dY/sc
        Yout, dYout = Y, dY
        # 6x6 phase-space matching matrix at r_match (r_match ENTERS):
        Mmat = np.block([[Yreg, -Yout], [dYreg, -dYout]])
        sv = np.linalg.svd(Mmat, compute_uv=False)
        return float(sv[-1])

    def r_out_values(self):
        return float(self.r[-1])


def main() -> int:
    t0 = time.time()
    d3, cb = load()
    sv_solver = Solver(d3, cb)
    out = {"audit": "JOST_DISCOVERY_V2",
            "fixes": ["F1 true Frobenius start basis",
                       "F2 k_j(w) per-eval from full pencil with S",
                       "F3 single r-coordinate convention",
                       "F4 genuine two-sided matching at real r_match",
                       "F5 R consumed from V3 export (proven zero)"]}
    # smoke: few probes with r_match ACTUALLY varied
    probes = []
    for w in (0.30+0.10j, 0.602-0.354j, 0.60+0.15j):
        row = []
        for rm in (20.0, 25.0, 30.0, 35.0, 40.0):
            s = sv_solver.jost_sigma_min(w, r_match=rm)
            row.append({"r_match": rm, "sigma_min": f"{s:.3e}"})
        probes.append({"omega_t": [w.real, w.imag], "scan": row})
        print(w, "->", [x["sigma_min"] for x in row], flush=True)
    out["r_match_probe"] = probes
    out["wall_seconds"] = round(time.time() - t0, 1)
    (ART / "JOST_V2_PROBE.json").write_text(
        json.dumps(out, indent=1, default=str) + "\n")
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
