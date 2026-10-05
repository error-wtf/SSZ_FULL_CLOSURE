#!/usr/bin/env python3
"""NATIVE_MODE_D2_ATLAS: Li-D2/IPR of the native FEM box modes.

Computes the localization diagnostics on the NATIVE radial mode shapes
psi(r) (3047-point FEM grid, Dirichlet box).  These are candidate
structures — the global operator later decides physical status.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.qnm.localization_diagnostics import diagnose_mode  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native-dir", type=Path,
                    default=ROOT / "data/generated/spectral/native_window")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "data/generated/spectral/NATIVE_MODE_D2_ATLAS.json")
    args = ap.parse_args()

    out = {"audit": "NATIVE_MODE_D2_ATLAS_V1",
           "note": ("D2/IPR of the native FEM box modes (Dirichlet, radial "
                    "grid).  Candidate structures — the global operator "
                    "later decides physical status."),
           "modes": []}

    for npz_path in sorted(args.native_dir.glob("L*_native_window_spectroscopy.npz")):
        L = int(npz_path.name.split("_")[0][1:])
        z = np.load(npz_path)
        r, psi, om = z["r"], z["psi"], z["omega"]
        n_r = psi.shape[1]  # radial points axis
        ones_K = np.ones((n_r, 1, 1))
        for m in range(psi.shape[0]):
            for f in range(psi.shape[2]):
                p = psi[m, :, f]
                norm = np.sqrt(np.abs(np.trapezoid(p * p, r)))
                if norm < 1e-12:
                    continue
                p = p / norm
                loc = diagnose_mode(m, p.reshape(-1, 1), ones_K, np.ones(n_r))
                out["modes"].append({
                    "L": L, "mode": m, "field": f,
                    "omega": float(om[m]),
                    "d2": loc.d2_scaling,
                    "ipr_grid": loc.ipr,
                    "class": loc.classification,
                })
        print(f"L={L}: {psi.shape[0]} modes")

    args.output.write_text(json.dumps(out, indent=1) + "\n")
    from collections import Counter
    c = Counter(m["class"] for m in out["modes"])
    d2s = [m["d2"] for m in out["modes"] if m["d2"] == m["d2"]]
    print(f"n mode-fields: {len(out['modes'])} | classes: {dict(c)}")
    print(f"D2 range: {min(d2s):.3f} .. {max(d2s):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
