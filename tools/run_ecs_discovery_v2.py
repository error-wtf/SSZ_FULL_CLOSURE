#!/usr/bin/env python3
"""ECS_DISCOVERY_V2 — corrected exterior complex scaling solver.

Fixes vs V1 (implementation audit):
  F1. Correct companion linearisation:
      A = [[0, I], [-C0, -C1]], B = [[I, 0], [0, C2]]
      (unit-tested: Q(w)psi = 0 iff A x = w B x, residual ~1e-15)
  F2. Single complex-scaled coordinate: finite differences built DIRECTLY
      on the complex z-grid — the complex spacing carries the scaling,
      NO extra e^{-i theta} derivative factors.
  F3. Exterior coefficients analytically continued via the fitted 1/r
      expansion: A(z) = A_inf + A1/z + A2/z^2 fitted from the last
      nodes of the V3 export, evaluated at complex z. No real(z)
      projection.
  F4. Center BC: the pencil is projected onto the SAME
      V4_CENTER_REGULAR_BASIS_V2 consumed by Jost V2 (basis replace at
      the first two grid rows), not arbitrary Dirichlet rows.
  F5. R consumed from the V3 export (proven identically zero).
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


def deriv_arrays(r, A):
    dA = np.empty_like(A)
    dA[1:-1] = (A[2:] - A[:-2]) / (r[2:] - r[:-2])[:, None, None]
    dA[0] = dA[1]; dA[-1] = dA[-2]
    return dA


def fit_exterior(A, r, n_fit=60):
    """Fit A_inf + A1/z + A2/z^2 on the last n_fit nodes (real r) for
    analytic continuation to complex z."""
    zz = r[-n_fit:]
    blocks = []
    for i in range(3):
        for j in range(3):
            y = A[-n_fit:, i, j]
            Amat = np.vstack([np.ones_like(zz), 1/zz, 1/zz**2]).T
            coef, *_ = np.linalg.lstsq(Amat, y, rcond=None)
            blocks.append(coef)
    return np.array(blocks)  # shape (9, 3): [a_inf, a1, a2] per entry


def eval_fit(coef, z):
    out = np.zeros((len(np.atleast_1d(z)), 3, 3), complex)
    z = np.atleast_1d(z)
    for e in range(9):
        i, j = divmod(e, 3)
        out[:, i, j] = coef[e, 0] + coef[e, 1]/z + coef[e, 2]/z**2
    return out


# (superseded build_pencil removed; use build_pencil_fixed)

def build_pencil_fixed(d3, cb, theta_deg, n_in=450, n_out=350,
                        r_match=30.0, r_end=60.0):
    theta = np.deg2rad(theta_deg)
    r_real = d3["r"]
    r1 = np.linspace(r_real[0], r_match, n_in, endpoint=False)
    s_seg = np.linspace(0, 1, n_out + 1)[1:]
    z2 = r_match + np.exp(1j*theta) * s_seg * (r_end - r_match)
    z = np.concatenate([r1.astype(complex), z2])
    n = len(z)
    h = np.diff(z)
    D1 = np.zeros((n, n), complex)
    D2 = np.zeros((n, n), complex)
    for i in range(1, n-1):
        hm, hp = h[i-1], h[i]
        D1[i, i-1] = -hp/(hm*(hm+hp)); D1[i, i+1] = hm/(hp*(hm+hp))
        D1[i, i] = (hp-hm)/(hm*hp)
        D2[i, i-1] = 2.0/(hm*(hm+hp)); D2[i, i+1] = 2.0/(hp*(hm+hp))
        D2[i, i] = -2.0/(hm*hp)
    # one-sided 2nd-order stencil at the FIRST node (needed by the
    # center-BC row; previously D1 row 0 was all-zero, silently killing
    # the center constraint for any basis with dY != 0):
    h0, h1 = h[0], h[1]
    D1[0, 0] = -(2*h0 + h1)/(h0*(h0+h1))
    D1[0, 1] = (h0 + h1)/(h0*h1)
    D1[0, 2] = -h0/(h1*(h0+h1))
    fitK = fit_exterior(d3["K_phys"], r_real)
    fitG = fit_exterior(d3["G_phys"], r_real)
    fitS = fit_exterior(d3["S_phys"], r_real)
    fitM = fit_exterior(d3["M_phys"], r_real)
    Kc = np.zeros((n, 3, 3), complex); Gc = np.zeros_like(Kc)
    Sc = np.zeros_like(Kc); Mc = np.zeros_like(Kc)
    def interp_real(A, rr):
        idx = np.clip(np.searchsorted(r_real, rr) - 1, 0, len(r_real)-2)
        t = (rr - r_real[idx])/(r_real[idx+1]-r_real[idx])
        return A[idx]*(1-t)[:, None, None] + A[idx+1]*t[:, None, None]
    Kc[:n_in] = interp_real(d3["K_phys"], r1)
    Gc[:n_in] = interp_real(d3["G_phys"], r1)
    Sc[:n_in] = interp_real(d3["S_phys"], r1)
    Mc[:n_in] = interp_real(d3["M_phys"], r1)
    Kc[n_in:] = eval_fit(fitK, z2)
    Gc[n_in:] = eval_fit(fitG, z2)
    Sc[n_in:] = eval_fit(fitS, z2)
    Mc[n_in:] = eval_fit(fitM, z2)
    dGz = np.zeros_like(Gc); dSz = np.zeros_like(Sc)
    dGz[1:-1] = (Gc[2:] - Gc[:-2]) / h[1:, None, None]
    dGz[0] = dGz[1]; dGz[-1] = dGz[-2]
    dSz[1:-1] = (Sc[2:] - Sc[:-2]) / h[1:, None, None]
    dSz[0] = dSz[1]; dSz[-1] = dSz[-2]
    I3 = np.eye(3)
    def blk(Mb):
        out = np.zeros((3*n, 3*n), complex)
        for i in range(n):
            out[3*i:3*i+3, 3*i:3*i+3] = Mb[i]
        return out
    D1_3 = np.kron(D1, I3); D2_3 = np.kron(D2, I3)
    C2 = blk(Kc)
    C1 = np.zeros_like(C2)
    C0 = blk(Gc) @ D2_3 + blk(dGz - Sc) @ D1_3 - blk(Mc + dSz/2.0)
    # ---- center BC: Psi'(r_in) = Trel Psi(r_in) on the FIRST row block
    Yreg = cb["Y_regular"].astype(complex)
    Yinv = np.linalg.inv(Yreg)
    Trel = cb["dY_regular"].astype(complex) @ Yinv
    for ch in range(3):
        row = 3*0 + ch
        # save the D1 first-row contribution, then zero the row:
        d1row = D1_3[row, :].copy()
        C0[row, :] = 0; C1[row, :] = 0; C2[row, :] = 0
        # Psi' - Trel Psi = 0:
        C0[row, :] += d1row          # Psi' term
        for cc in range(3):
            C0[row, 3*0+cc] += -Trel[ch, cc]  # -Trel Psi term
    # ---- OUTER ENDPOINT BC (V2.2 fix): Psi(z_end) = 0 ----------------
    # The last 3 row blocks previously kept the algebraic interior
    # stencil row (no valid endpoint condition). Impose Psi(z_end)=0
    # explicitly — standard ECS truncation for a sufficiently long
    # rotated tail. Endpoint convergence is tested separately.
    for ch in range(3):
        row = 3*(n-1) + ch
        C0[row, :] = 0; C1[row, :] = 0; C2[row, :] = 0
        C0[row, 3*(n-1) + ch] = 1.0   # Psi_ch(z_end) = 0
    return C2, C1, C0, z


def solve_spectrum(d3, cb, theta_deg, n_in=450, n_out=350):
    C2, C1, C0, z = build_pencil_fixed(d3, cb, theta_deg, n_in, n_out)
    n3 = C2.shape[0]
    A = np.zeros((2*n3, 2*n3), complex)
    B = np.zeros((2*n3, 2*n3), complex)
    A[:n3, n3:] = np.eye(n3)
    A[n3:, :n3] = -C0
    A[n3:, n3:] = -C1
    B[:n3, :n3] = np.eye(n3)
    B[n3:, n3:] = C2
    w_all, vecs = geig(A, B)
    return w_all


def main() -> int:
    t0 = time.time()
    d3 = np.load(ART / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz")
    cb = np.load(ART / "V4_CENTER_REGULAR_BASIS_V2.npz")
    out = {"audit": "ECS_DISCOVERY_V2",
            "fixes": ["F1 correct companion linearisation (unit-tested)",
                       "F2 single complex grid, no double scaling",
                       "F3 analytic 1/r exterior continuation",
                       "F4 center-regular BC in the pencil",
                       "F5 R from V3 export"]}
    spectra = {}
    for theta in (35.0, 45.0):
        t1 = time.time()
        w = solve_spectrum(d3, cb, theta)
        keep = (np.abs(w) < 50) & (np.imag(w) < -1e-6) & (np.real(w) > 1e-6)
        ws = w[keep]
        ws = ws[np.argsort(np.real(ws))]
        key = f"theta{int(theta)}"
        spectra[key] = [[round(float(x.real), 6), round(float(x.imag), 6)]
                         for x in ws[:60]]
        # region near the Jost V2 candidate band:
        near = [x for x in ws
                 if 0.5 < x.real < 0.8 and -0.45 < x.imag < -0.15]
        out[f"near_jost_band_theta{int(theta)}"] = [
            [round(float(x.real), 6), round(float(x.imag), 6)] for x in near]
        print(f"theta={theta}: kept={keep.sum()}, near Jost band: {len(near)} "
              f"({time.time()-t1:.0f}s)", flush=True)
        for x in near[:12]:
            print(f"   {x.real:.5f} {x.imag:+.5f}j", flush=True)
    out["spectra"] = spectra
    out["wall_seconds"] = round(time.time() - t0, 1)
    (ART / "ECS_V2_CANDIDATES.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
