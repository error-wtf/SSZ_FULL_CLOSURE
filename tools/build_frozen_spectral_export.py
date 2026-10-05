"""Produce the frozen spectral operator export (Producer Contract, SSZ-B1/B2).

Reads the frozen production member, restricts to the declared healthy window
(u in [0.620, 0.700]), reduces the profile operator per L in DEFAULT_L via the
verified reducer, runs the fail-closed finite-L health diagnostics, and writes
a hash-bound npz + provenance JSON for independent bridge validation.

The known nonpositive-Kinetic block at the outer boundary zone (u 0.7013-0.7100)
is explicitly excluded by the window and declared in the provenance record —
it is NOT hidden and the export does NOT claim global kinetic health.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.config import DEFAULT_L, SLOT_NAMES  # noqa: E402
from ssz_p5.geometry.p5 import from_frame  # noqa: E402
from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)
from ssz_p5.provenance.manifest import sha256  # noqa: E402
from ssz_p5.reducer.canonical import reduce_profile  # noqa: E402
from ssz_p5.stability.finite_l import stability_diagnostics  # noqa: E402
from ssz_p5.types import Coefficients41  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
U_MIN, U_MAX = 0.620, 0.700


def main() -> int:
    b = build_onshell_central(ROOT)
    frame = b.direct41.sort_values("x").reset_index(drop=True)
    u_all = 1.0 / frame["x"].to_numpy(float)
    win = (u_all >= U_MIN) & (u_all <= U_MAX)
    sub = frame[win].reset_index(drop=True)
    u_win = 1.0 / sub["x"].to_numpy(float)
    print(f"healthy-window rows: {len(sub)}  u: {u_win.min():.4f} -> {u_win.max():.4f}")

    c = Coefficients41(
        from_frame(sub),
        {s: sub[s].to_numpy(float) for s in SLOT_NAMES},
        "healthy_window",
    )
    health = {}
    ops = {}
    for L in DEFAULT_L:
        op = reduce_profile(c, int(L))
        diag = stability_diagnostics(op)  # fail-closed on nonpositive K
        health[str(int(L))] = diag
        ops[int(L)] = op
        print(
            f"L={int(L)}: min_eig_K={diag['min_eig_K']:.3e} "
            f"min_radial={diag['min_radial']:.3e} pass={diag['pass']}"
        )

    arrays: dict[str, np.ndarray] = {}
    for L in DEFAULT_L:
        op = ops[int(L)]
        prefix = f"{int(L)}_"
        arrays[prefix + "r"] = op.r
        for m in ("K", "R", "G", "S", "M"):
            arrays[prefix + m] = getattr(op, m)
        arrays[prefix + "Dh1"] = op.pivots.Dh1
        arrays[prefix + "DeltaV"] = op.pivots.DeltaV
        arrays[prefix + "pivotA0"] = op.pivots.pivotA0

    out_npz = ROOT / "data/generated/spectral/FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.npz"
    np.savez(out_npz, **arrays)
    npz_sha = sha256(out_npz)
    print("saved:", out_npz, out_npz.stat().st_size, "bytes")

    prov = {
        "export": "FROZEN_SPECTRAL_OPERATOR_EXPORT_V1",
        "producer": "tools/build_frozen_spectral_export.py (SSZ repo)",
        "member_stream": "electric_hybrid_onshell_central.direct41",
        "healthy_window_u": [U_MIN, U_MAX],
        "window_rows": int(len(sub)),
        "L_values": [int(L) for L in DEFAULT_L],
        "npz_sha256": npz_sha,
        "locked_conventions": {
            "reducer": "ssz_p5_profile_operator_reducer_JET9D8_2026-09-16",
            "reducer_sha256": sha256(
                ROOT / "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py"
            ),
            "derivatives": "JET9D8 window=9 degree=8",
            "v12": "-v6/(2h)",
            "coordinate": "x = r/r_s, u = 1/x; radial operator matrices on x-grid",
        },
        "known_limitations": {
            "outer_boundary_nonpositive_K": (
                "u in [0.7013, 0.7100] (x 1.4085-1.4259) shows min eig_K ~ -0.16 "
                "(scale 1.9e4); block contiguous at the domain edge, matches the "
                "documented suspect boundary zone (phi_max > 1). EXCLUDED by the "
                "window, NOT hidden. Global kinetic health is NOT claimed."
            ),
            "sector": "single coupled 3-field sector from the reducer; parity "
            "decomposition not yet exported",
            "bc": "no boundary conditions exported — this is an operator "
            "coefficient export, not a QNM claim",
        },
        "health_L": health,
    }
    out_json = ROOT / (
        "data/generated/spectral/FROZEN_SPECTRAL_OPERATOR_EXPORT_V1.json"
    )
    out_json.write_text(json.dumps(prov, indent=1) + "\n")
    print("saved:", out_json)
    print("npz_sha256:", npz_sha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
