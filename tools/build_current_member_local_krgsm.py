#!/usr/bin/env python3
"""Build every scientifically available current-member local K/R/G/S/M product.

The frozen current electric member is only defined on u in [0.61,0.71].
Therefore this exporter writes the complete finite-L matrices on the declared
production window 0.62<u<0.70 and labels them LOCAL_STRONG_FIELD.  It never
mislabels them as a center-to-infinity Direct Global KRGM product.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ssz_p5.numerics import module
from ssz_p5.production.electric_hybrid_onshell_central import build_onshell_central

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/spectral_selection/current_member_local_krgm"
MANIFEST = ROOT / "data/generated/phase2_q2/ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"
LS = (6, 12, 20, 42, 110, 420, 1000)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    manifest = json.loads(MANIFEST.read_text())
    build = build_onshell_central(ROOT)
    d = build.direct41.copy().sort_values("x").reset_index(drop=True)
    u = d.u.to_numpy(float)
    keep = (u > 0.62) & (u < 0.70)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")

    OUT.mkdir(parents=True, exist_ok=True)
    artifacts = []
    for L in LS:
        a = red.canonical_audit(d, L)
        path = OUT / f"L{L:04d}.npz"
        np.savez_compressed(
            path,
            x=d.x.to_numpy(float)[keep],
            u=u[keep],
            K=np.asarray(a["K"], float)[keep],
            R=np.asarray(a["R"], float)[keep],
            G=np.asarray(a["G"], float)[keep],
            S=np.asarray(a["S"], float)[keep],
            M=np.asarray(a["M"], float)[keep],
        )
        artifacts.append(
            {"L": L, "path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        )

    report = {
        "schema_version": "1.0",
        "classification": "LOCAL_STRONG_FIELD_KRGSM_EXPORT",
        "global_direct_certificate": False,
        "member_hash": manifest["member_hash"],
        "domain": {"u_min_exclusive": 0.62, "u_max_exclusive": 0.70},
        "L_values": list(LS),
        "artifacts": artifacts,
        "warning": (
            "These matrices are the complete available current-member local "
            "production-window products. They are not center-to-infinity and "
            "must not unlock the QNM gate."
        ),
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
