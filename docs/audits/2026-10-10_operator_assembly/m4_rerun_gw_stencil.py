#!/usr/bin/env python3
"""M4-RERUN B: SSZ spectral-weight diagnostics with CORRECT stencil (2h).

EXACT replica of tools/run_spectral_weight_diagnostics_gw.py assemble +
weight_spectrum. ONLY change: dGz/dSz denominators h -> 2h (the certified
M3 factor-2 bug). Same grid n=400 (fast certification pass, NOT the full
800x2x2 production run), same window, same probe convention.
Output: scratch JSON with production-vs-fixed comparison at n=400.
"""
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
NPZ = ROOT / "data/generated/spectral/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
OUT = Path("/root/.hermes/cache/scratch/m4_rerun_gw")
OUT.mkdir(parents=True, exist_ok=True)

# production run params (from file constants, verified earlier)
W_LO, W_HI = 0.05, 1.5
N_W_CERT = 400          # certification pass (production was 800)
RESOLUTIONS = [400]
BCS = ["dirichlet", "robin"]


def assemble(npz, n, bc, stencil_fix):
    r_full = npz["r"]
    K = npz["K_phys"]; G = npz["G_phys"]; S = npz["S_phys"]; M = npz["M_phys"]
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
    denom = (2.0 if stencil_fix else 1.0) * hc
    dGz[1:-1] = (Gc[2:] - Gc[:-2]) / denom[:, None, None]
    dGz[0] = dGz[1]; dGz[-1] = dGz[-2]
    dSz[1:-1] = (Sc[2:] - Sc[:-2]) / denom[:, None, None]
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
    g00 = np.empty(len(ws))
    for k, w in enumerate(ws):
        A = C2 * ((w + 1e-6j) ** 2) + C0
        b = np.zeros(A.shape[0], dtype=complex); b[3 * probe_node] = 1.0
        rs = np.max(np.abs(A), axis=1); rs[rs == 0] = 1.0
        As = A / rs[:, None]; bs = b / rs
        try:
            g00[k] = abs(np.linalg.solve(As, bs)[3 * probe_node])
        except np.linalg.LinAlgError:
            g00[k] = np.nan
    return g00


def main():
    t0 = time.time()
    npz = np.load(NPZ, allow_pickle=True)
    ws = np.linspace(W_LO, W_HI, N_W_CERT)
    report = {"audit": "M4_RERUN_GW_STENCIL_FIX_N400",
              "note": "certification pass n=400; production was n=600/1000",
              "variants": {}}
    for n in RESOLUTIONS:
        for bc in BCS:
            t1 = time.time()
            C2p, C0p, _ = assemble(npz, n, bc, stencil_fix=False)
            gp = weight_spectrum(C2p, C0p, ws, probe_node=n // 4)
            C2f, C0f, _ = assemble(npz, n, bc, stencil_fix=True)
            gf = weight_spectrum(C2f, C0f, ws, probe_node=n // 4)
            np.save(OUT / f"g_n{n}_{bc}_prod.npy", gp)
            np.save(OUT / f"g_n{n}_{bc}_fixed.npy", gf)
            rel = np.abs(gf - gp) / np.maximum(np.abs(gp), 1e-30)
            report["variants"][f"n{n}_{bc}"] = {
                "prod_g_max": round(float(np.nanmax(gp)), 6),
                "fixed_g_max": round(float(np.nanmax(gf)), 6),
                "prod_g_median": round(float(np.nanmedian(gp)), 8),
                "fixed_g_median": round(float(np.nanmedian(gf)), 8),
                "max_rel_dev": round(float(np.nanmax(rel)), 4),
                "median_rel_dev": round(float(np.nanmedian(rel)), 4),
                "wall_s": round(time.time() - t1, 1),
            }
            print(f"n={n} {bc}: prod_max={np.nanmax(gp):.4g} fixed_max={np.nanmax(gf):.4g} "
                  f"med_rel_dev={np.nanmedian(rel):.3%} ({time.time()-t1:.0f}s)", flush=True)
            # peak tables for both curves (same 6-sigma local-MAD filter)
            for tag, g in (("prod", gp), ("fixed", gf)):
                if np.isnan(g).any():
                    good = ~np.isnan(g)
                    g = np.interp(ws, ws[good], g[good])
                peaks, med, mad = peak_table(ws, g)
                report["variants"][f"n{n}_{bc}"][f"peaks_{tag}"] = {
                    "n": len(peaks), "top5": sorted(peaks, key=lambda x: -x["g"])[:5]}
    report["wall_s"] = round(time.time() - t0, 1)
    (OUT / "M4_RERUN_GW_SUMMARY.json").write_text(json.dumps(report, indent=1) + "\n")
    print("wrote", OUT / "M4_RERUN_GW_SUMMARY.json", flush=True)


if __name__ == "__main__":
    main()
