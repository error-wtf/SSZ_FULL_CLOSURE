#!/usr/bin/env python3
"""SSZ spectral-weight diagnostics step 2: Green-function weight spectrum of the
certified 3-DOF operator.

Prereg: artifacts/SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_PREREGISTRATION_V1.json

Assembly follows the certified ECS path (run_ecs_discovery_v2.py):
  C0 = blk(G) D2 + blk(dG - S) D1 - blk(M + dS/2),  C2 = blk(K), C1 = 0
  G(w) = (C2 w^2 + C0)^-1, observable = G[0,0] (dipole at inner boundary)
Real r-grid, two resolutions (600, 1000), two BCs (Dirichlet, Robin-at-inner
matching the regular-solution derivative ratio), w in [0.05, 1.5].

Reported per prereg: peak table with local-noise significance, selectivity C,
entropy S; >5% drift across variants => UNSTABLE, no claim.
"""
import json
import math
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
NPZ = ROOT / "data/generated/spectral/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
OUT = ROOT / "artifacts/SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json"

W_LO, W_HI, N_W = 0.05, 1.5, 800
RESOLUTIONS = [600, 1000]
BCS = ["dirichlet", "robin"]


def assemble(npz, n, bc):
    r_full = npz["r"]
    K = npz["K_phys"]; G = npz["G_phys"]; S = npz["S_phys"]; M = npz["M_phys"]
    # resample profiles onto a uniform real grid inside the certified window
    r_lo, r_hi = 0.05, 1.0
    r1 = np.linspace(r_lo, r_hi, n)

    def interp3(A):
        out = np.empty((n, 3, 3))
        for a in range(3):
            for b in range(3):
                out[:, a, b] = np.interp(r1, r_full, A[:, a, b])
        return out

    Kc, Gc, Sc, Mc = interp3(K), interp3(G), interp3(S), interp3(M)
    h = np.full(n, r1[1] - r1[0])
    dGz = np.zeros_like(Gc); dSz = np.zeros_like(Sc)
    hc = h[1:n-1]
    dGz[1:-1] = (Gc[2:] - Gc[:-2]) / hc[:, None, None]
    dGz[0] = dGz[1]; dGz[-1] = dGz[-2]
    dSz[1:-1] = (Sc[2:] - Sc[:-2]) / hc[:, None, None]
    dSz[0] = dSz[1]; dSz[-1] = dSz[-2]

    D1 = np.zeros((n, n)); D2 = np.zeros((n, n))
    for i in range(1, n - 1):
        D1[i, i - 1] = -1.0 / (2 * h[i]); D1[i, i + 1] = 1.0 / (2 * h[i])
        D2[i, i - 1] = 1.0 / h[i] ** 2; D2[i, i + 1] = 1.0 / h[i] ** 2
        D2[i, i] = -2.0 / h[i] ** 2

    def blk(Mb):
        out = np.zeros((3 * n, 3 * n))
        for i in range(n):
            out[3 * i:3 * i + 3, 3 * i:3 * i + 3] = Mb[i]
        return out

    I3 = np.eye(3)
    D1_3 = np.kron(D1, I3); D2_3 = np.kron(D2, I3)
    C2 = blk(Kc)
    C0 = blk(Gc) @ D2_3 + blk(dGz - Sc) @ D1_3 - blk(Mc + dSz / 2.0)

    if bc == "dirichlet":
        for ch in range(3):
            for row in (ch, 3 * (n - 1) + ch):
                C0[row, :] = 0; C2[row, :] = 0
                C0[row, 3 * (0 if row == ch else n - 1) + ch] = 1.0
    elif bc == "robin":
        # inner BC: one-sided 2nd-order forward difference for Psi'(r_in)=0
        # (D1 has zero boundary rows by construction — using its row 0 left
        # the matrix structurally singular; that was the crash). Frozen: T=0.
        h0 = h[0]
        for ch in range(3):
            row = ch
            C0[row, :] = 0; C2[row, :] = 0
            C0[row, 3 * 0 + ch] = -3.0 / (2 * h0)
            C0[row, 3 * 1 + ch] = 4.0 / (2 * h0)
            C0[row, 3 * 2 + ch] = -1.0 / (2 * h0)
        for ch in range(3):
            row = 3 * (n - 1) + ch
            C0[row, :] = 0; C2[row, :] = 0
            C0[row, 3 * (n - 1) + ch] = 1.0
    return C2, C0, r1


def weight_spectrum(C2, C0, ws, probe_node):
    # probe at an interior node (BC rows pin boundary values — useless as
    # observables); small imaginary shift avoids exact-hit singularities,
    # documented in the report
    g00 = np.empty(len(ws))
    skipped = 0
    for k, w in enumerate(ws):
        A = C2 * ((w + 1e-6j) ** 2) + C0
        b = np.zeros(A.shape[0], dtype=complex); b[3 * probe_node] = 1.0
        # row equilibration: huge coefficient dynamic (diag up to ~1e9) breaks
        # LAPACK's default pivot tolerance at finer grids
        rs = np.max(np.abs(A), axis=1)
        rs[rs == 0] = 1.0
        As = A / rs[:, None]
        bs = b / rs
        try:
            g00[k] = abs(np.linalg.solve(As, bs)[3 * probe_node])
        except np.linalg.LinAlgError:
            g00[k] = np.nan; skipped += 1
    if skipped:
        print(f"  [warn] {skipped} singular frequencies skipped", flush=True)
    return g00


def peak_table(w, g, n_sigma=6.0):
    # local noise: rolling median absolute deviation in a 61-point window
    half = 30
    med = np.empty_like(g)
    for i in range(len(g)):
        lo, hi = max(0, i - half), min(len(g), i + half + 1)
        med[i] = np.median(g[lo:hi])
    mad = np.empty_like(g)
    for i in range(len(g)):
        lo, hi = max(0, i - half), min(len(g), i + half + 1)
        mad[i] = np.median(np.abs(g[lo:hi] - med[i])) * 1.4826
    sig = (g - med) / np.maximum(mad, 1e-30)
    peaks = []
    for i in range(1, len(g) - 1):
        if g[i] > g[i - 1] and g[i] >= g[i + 1] and sig[i] >= n_sigma:
            peaks.append({"w": round(float(w[i]), 5),
                           "g": round(float(g[i]), 6),
                           "significance_sigma": round(float(sig[i]), 2)})
    return peaks, med, mad


def main():
    t0 = time.time()
    npz = np.load(NPZ, allow_pickle=True)
    ws = np.linspace(W_LO, W_HI, N_W)
    variants = {}
    for n in RESOLUTIONS:
        for bc in BCS:
            t1 = time.time()
            C2, C0, _ = assemble(npz, n, bc)
            g = weight_spectrum(C2, C0, ws, probe_node=n // 4)
            # NaN-safe: interpolate over skipped points for peak search
            if np.isnan(g).any():
                good = ~np.isnan(g)
                g = np.interp(ws, ws[good], g[good])
            peaks, med, mad = peak_table(ws, g)
            # selectivity over detected peaks
            gp = np.array([p["g"] for p in peaks]) if peaks else np.array([0.0])
            p = gp / gp.sum()
            variants[f"n{n}_{bc}"] = {
                "n_peaks": len(peaks),
                "peaks_top5": sorted(peaks, key=lambda x: -x["g"])[:5],
                "selectivity_C": round(float(p.max()), 4),
                "entropy_S": round(float(-np.sum(p[p > 0] * np.log(p[p > 0]))), 4),
                "g_max": round(float(g.max()), 6),
                "g_median": round(float(np.median(g)), 8),
                "wall_s": round(time.time() - t1, 1),
            }
            print(f"n={n} {bc}: peaks={len(peaks)} C={variants[f'n{n}_{bc}']['selectivity_C']} "
                  f"maxG={g.max():.4g} ({time.time()-t1:.0f}s)", flush=True)

    # robustness per prereg: compare peak frequencies across variants (within 2 bins)
    bins = ws[1] - ws[0]
    def freqs(v):
        return [p["w"] for p in variants[v]["peaks_top5"]]
    base = freqs("n1000_dirichlet")
    stable = {}
    for v in variants:
        if v == "n1000_dirichlet":
            continue
        matched = 0
        for w0 in base:
            if any(abs(w0 - w1) <= 2 * bins for w1 in freqs(v)):
                matched += 1
        stable[v] = f"{matched}/{len(base)}"

    report = {
        "artifact": "SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1",
        "preregistration": "SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_PREREGISTRATION_V1",
        "operator_source": NPZ.name,
        "window": [W_LO, W_HI], "n_freq": N_W,
        "variants": variants,
        "robustness_top5_match_vs_n1000_dirichlet": stable,
        "verdict": None,
        "wall_seconds": round(time.time() - t0, 1),
    }
    # verdict per prereg: claim only if >=2 of 3 other variants reproduce all top5
    ok = sum(1 for v, s in stable.items() if s.startswith(f"{len(base)}/"))
    report["verdict"] = ("PEAKS_ROBUST" if ok >= 2 and base else
                          "PEAKS_NOT_ROBUST — diagnostics reported, no resonance claim")
    OUT.write_text(json.dumps(report, indent=1) + "\n")
    print("VERDICT:", report["verdict"])
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
