#!/usr/bin/env python3
"""Global P5 test-field QNM/resonance scan with a complex absorbing layer.

This is an independent global center-to-exterior spectral test on the frozen P5
geometry.  It is deliberately separate from the coupled HSVT QNM certificate.

The radial equations are put on tortoise coordinate r_*:
    [-d^2/dr_*^2 + V_l(r)] psi = omega^2 psi
with
    V_EM = f l(l+1)/r^2
and the minimally coupled scalar potential used by the project.

Regular-center behavior is enforced by Dirichlet at the innermost resolved
center node for l>=1.  Outgoing behavior is approximated with a smooth complex
absorbing potential (CAP) in the weak exterior.  Resolution, CAP onset and CAP
strength are varied.  Roots are only called CAP-stable diagnostics if they stay
within a declared complex-frequency tolerance across all configurations.

This is a real global spectral computation on the P5 geometry, but CAP-stable
roots are not promoted to the coupled HSVT QNM spectrum.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.sparse import diags
from scipy.sparse.linalg import eigs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/qnm/P5_GLOBAL_TESTFIELD_QNM_CAP.json"
GEOM = ROOT / "data/authoritative/ssz_p5_F2A_GLOBAL_CANONICAL_KRG_PRINCIPAL_CINF_FINAL_2026-09-15(2).csv"

KINDS = ("maxwell", "scalar")
LS = (2, 3, 6)
CONFIGS = (
    (700, 0.72, 0.5),
    (900, 0.75, 0.5),
    (1100, 0.75, 0.8),
    (900, 0.78, 0.8),
)
MATCH_TOL = 1.0e-2


def geometry():
    d = pd.read_csv(GEOM, usecols=["x", "f", "h"]).sort_values("x")
    r = d.x.to_numpy(float)
    f = d.f.to_numpy(float)
    h = d.h.to_numpy(float)
    invc = 1.0 / np.sqrt(f * h)
    rs = np.zeros_like(r)
    rs[1:] = np.cumsum(0.5 * (invc[1:] + invc[:-1]) * np.diff(r))
    sq = np.sqrt(f * h)
    dsdr = np.gradient(sq, r, edge_order=2)
    return r, f, h, rs, dsdr


def wkb_guess(kind, ell, r, f, h, rs, dsdr):
    if kind == "maxwell":
        V = f * ell * (ell + 1) / r**2
    else:
        V = f * (ell * (ell + 1) / r**2 + dsdr / (r * np.sqrt(f / h)))
    mask = (r > 1.45) & (r < 1.55)
    ids = np.flatnonzero(mask)
    i = ids[int(np.argmax(V[mask]))]
    cs = CubicSpline(rs, V)
    v0 = float(V[i])
    v2 = float(cs(rs[i], 2))
    z = np.sqrt(v0 - 0.5j * np.sqrt(max(0.0, -2.0 * v2)))
    return z if z.imag <= 0 else -z


def cap_roots(kind, ell, N, onset_frac, eta, sigma, r, f, h, rs, dsdr):
    sr = PchipInterpolator(rs, r)
    sf = PchipInterpolator(rs, f)
    sh = PchipInterpolator(rs, h)
    sd = PchipInterpolator(rs, dsdr)

    t = np.linspace(rs[0], rs[-1], N)
    rr = sr(t)
    ff = sf(t)
    hh = sh(t)
    if kind == "maxwell":
        V = ff * ell * (ell + 1) / rr**2
    else:
        V = ff * (ell * (ell + 1) / rr**2 + sd(t) / (rr * np.sqrt(ff / hh)))

    tt = t[1:-1]
    vv = V[1:-1]
    dx = float(tt[1] - tt[0])
    n = len(tt)
    H = diags(
        [np.full(n - 1, -1 / dx**2), np.full(n, 2 / dx**2) + vv, np.full(n - 1, -1 / dx**2)],
        [-1, 0, 1],
        format="csc",
    ).astype(complex)

    onset = t[0] + onset_frac * (t[-1] - t[0])
    q = np.clip((tt - onset) / (t[-1] - onset), 0, None)
    H += diags(-1j * eta * q**4, 0, format="csc")

    vals = eigs(
        H,
        k=28,
        sigma=sigma**2,
        which="LM",
        return_eigenvectors=False,
        tol=1e-10,
        maxiter=30000,
    )
    w = np.sqrt(vals)
    w = np.asarray([z if z.imag <= 0 else -z for z in w])
    return w[np.argsort(np.abs(w - sigma))]


def common_roots(sets):
    out = []
    for z in sets[0]:
        d = [float(np.min(np.abs(s - z))) for s in sets[1:]]
        if d and max(d) <= MATCH_TOL:
            out.append((z, max(d)))
    out.sort(key=lambda p: p[1])
    return out


def main():
    r, f, h, rs, dsdr = geometry()
    cases = []
    total_growing = 0
    for kind in KINDS:
        for ell in LS:
            sigma = wkb_guess(kind, ell, r, f, h, rs, dsdr)
            sets = []
            cfgs = []
            for N, onset, eta in CONFIGS:
                roots = cap_roots(kind, ell, N, onset, eta, sigma, r, f, h, rs, dsdr)
                sets.append(roots)
                cfgs.append(
                    {
                        "N": N,
                        "cap_onset_fraction": onset,
                        "cap_strength": eta,
                        "nearest_roots": [
                            {"re": float(z.real), "im": float(z.imag)}
                            for z in roots[:12]
                        ],
                    }
                )
            stable = common_roots(sets)
            grow = sum(z.imag > 1e-8 for z, _ in stable)
            total_growing += grow
            cases.append(
                {
                    "field": kind,
                    "ell": ell,
                    "wkb_outer_light_ring_guess": {
                        "re": float(sigma.real),
                        "im": float(sigma.imag),
                    },
                    "configs": cfgs,
                    "cap_stable_root_count": len(stable),
                    "cap_stable_growing_count": grow,
                    "cap_stable_roots": [
                        {
                            "re": float(z.real),
                            "im": float(z.imag),
                            "max_cross_config_delta": float(delta),
                        }
                        for z, delta in stable[:20]
                    ],
                }
            )

    report = {
        "status": "GLOBAL_P5_TESTFIELD_QNM_CAP_SCAN_COMPLETE",
        "geometry": str(GEOM.relative_to(ROOT)),
        "center_to_exterior_x_range": [float(r.min()), float(r.max())],
        "center_to_exterior_rstar_range": [float(rs.min()), float(rs.max())],
        "method": "finite-difference tortoise operator + complex absorbing potential",
        "match_tolerance": MATCH_TOL,
        "cases": cases,
        "cap_stable_growing_modes_in_scanned_cases": int(total_growing),
        "physical_scope": (
            "Global Maxwell/scalar test-field resonance diagnostic on the P5 geometry. "
            "Not the coupled HSVT same-operator QNM certificate."
        ),
        "caution": (
            "CAP roots can include discretized continuum/pseudomodes. A root needs stronger "
            "CAP-trajectory or Jost/ECS convergence before publication as a physical QNM."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": report["status"],
        "stable_growing": report["cap_stable_growing_modes_in_scanned_cases"],
        "summary": {
            f"{c['field']}_l{c['ell']}": {
                "stable": c["cap_stable_root_count"],
                "growing": c["cap_stable_growing_count"],
                "best": c["cap_stable_roots"][:3],
            }
            for c in cases
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
