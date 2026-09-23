#!/usr/bin/env python3
"""Freeze the CURRENT same-action electric production member (STEP4).

Builds the locked on-shell Central electric member (build_onshell_central),
exports its full action stream (profiles + jets) and a canonical manifest,
computes the member_hash (SHA256 over the canonical CSV), and updates
MODEL_LOCK.json deterministically.

The historical zero-vector member (SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16,
A0' = 0, epsilon_Y deformation) is NOT reused; it is recorded as
historical-only context.

Run: PYTHONPATH=src python tools/freeze_electric_production_member.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.production.electric_hybrid_onshell_central import (  # noqa: E402
    build_onshell_central,
)

OUT_DIR = ROOT / "data" / "generated" / "phase2_q2"
CSV_PATH = OUT_DIR / "ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv"
MANIFEST_PATH = OUT_DIR / "ELECTRIC_PRODUCTION_MEMBER_CURRENT.json"


def git_head() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                          text=True, cwd=ROOT).stdout.strip()


def main() -> int:
    build = build_onshell_central(ROOT)
    action: pd.DataFrame = build.action.copy()

    # A0 profile by trapezoid integration of A0' on the (monotone-in-x) grid,
    # normalized to A0(r -> max x) = 0 gauge (U(1): only A0' is physical)
    import numpy as np
    x = action.x.to_numpy(float)
    ap = action.A0prime.to_numpy(float)
    dx = np.diff(x)
    A0_raw = np.concatenate([[0.0], np.cumsum(0.5*(ap[1:] + ap[:-1])*dx)])
    action["A0"] = A0_raw - A0_raw[-1]

    diag = build.diagnostics
    cols = ["u", "x", "phi", "f", "h", "phiprime", "A0", "A0prime", "X",
            "Fbg_action", "Ybg_action",
            "f2", "f2X", "f2F", "f2Y", "f2phi", "f2phiphi",
            "f3", "f3X", "f3phi", "f3phiX",
            "f4", "f4X", "f4XX", "f4phi", "f4phiX", "tf4", "tf4phi",
            "G2XX_lift", "G2Xphi_lift", "G2phiphi_lift"]
    missing = [c for c in cols if c not in action.columns]
    assert not missing, f"missing member columns: {missing}"
    stream = action[cols].sort_values("u").reset_index(drop=True)
    csv_text = stream.to_csv(index=True, float_format="%.17e", lineterminator="\n")
    CSV_PATH.write_text(csv_text)
    member_hash = hashlib.sha256(csv_text.encode("utf-8")).hexdigest()

    u = stream.u.to_numpy(float)
    mask_prod = (u > 0.62) & (u < 0.70)
    manifest = {
        "schema_version": "1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_head(),
        "name": "SSZ P5 electric production member (current, same-action, "
                "on-shell Central, Delta22_SVT-derived action)",
        "supersedes": {
            "member": "SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json",
            "reason": "historical ZERO-VECTOR member (A0'=0, eps_Y deformation); "
                      "rejected at the inner light ring (Eq85 requires "
                      "Tensor_F < 0 at H > 0 there); kept as control case only",
        },
        "member_hash": member_hash,
        "member_hash_scope": "sha256 over ELECTRIC_PRODUCTION_MEMBER_CURRENT.csv "
                             "(canonical UTF-8 serialization, index column included)",
        "rows": int(len(stream)),
        "domain": {
            "u_min": float(u.min()), "u_max": float(u.max()),
            "production_window": [0.62, 0.70],
            "production_rows": int(mask_prod.sum()),
            "coordinate": "u = r_s/r (r_s = 1)",
        },
        "interfaces": {
            "inner_u": 0.61, "inner_patch": "strong-field hybrid continuation",
            "outer_u": 0.71, "outer_patch": "outer handover continuation "
                                            "(S_SVT -> 0 across 0.71..0.715)",
            "partition_rule": "S_SVT + T_H = 1 (inner_targets.py), edges pinned 1/0",
        },
        "profiles": ["f", "h", "phi", "A0", "A0prime", "X", "Fbg_action",
                     "Ybg_action"],
        "action_jets": ["f2", "f2X", "f2F", "f2Y", "f2phi", "f2phiphi",
                        "f3", "f3X", "f3phi", "f3phiX",
                        "f4", "f4X", "f4XX", "f4phi", "f4phiX", "tf4", "tf4phi",
                        "G2XX_lift", "G2Xphi_lift", "G2phiphi_lift"],
        "f2Y_branch_rule": "f2Y SYMBOLIC through C_bg/Test A/Test B; production "
                           "branch postulate f2Y = 0 applied only at G13",
        "boundary_provenance": {
            "a2_solve_success": diag["a2_solve_success"],
            "a2_nfev": diag["a2_nfev"],
            "principal_replay_max_scaled": diag["principal_replay_max_scaled"],
            "baseline_streams": [
                "data/generated/absolute_attempt_2026-09-19/"
                "ELECTRIC_HYBRID_CENTRAL_SEARCH_BASELINE.csv",
                "data/generated/absolute_attempt_2026-09-20/"
                "ELECTRIC_HYBRID_PRINCIPAL_LATE_RAMP_FEASIBILITY_C6_CORRECTED.csv",
            ],
            "note": "on-shell Central single-action candidate with full G2 lift; "
                    "u in [0.61, 0.71]",
        },
        "action_definition_status": {
            "C_bg": "ALPHA*E00_full + BETA*E22_MH_slice + BETA*Delta22_SVT with "
                    "Delta22_SVT DERIVED (STEP3, direct theta-theta variation, "
                    "areal gauge after variation)",
            "genuine_svt_theta_theta": "CLOSED",
        },
        "diagnostics": diag,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=1) + "\n")

    # deterministic MODEL_LOCK update
    lock_path = ROOT / "MODEL_LOCK.json"
    lock = json.loads(lock_path.read_text())
    lock["action_member_reference"] = str(MANIFEST_PATH.relative_to(ROOT))
    lock["action_member_sha256"] = member_hash
    lock["action_member_stream"] = str(CSV_PATH.relative_to(ROOT))
    lock["action_member_frozen_at"] = manifest["timestamp"]
    lock["action_member_provenance"] = {
        "build": "build_onshell_central (electric_hybrid_onshell_central)",
        "supersedes_zero_vector_member": True,
        "delta22_svt": "DERIVED (STEP3_DELTA22_DIRECT_VARIATION.json)",
    }
    lock_path.write_text(json.dumps(lock, indent=1) + "\n")

    print(json.dumps({"member_hash": member_hash, "rows": len(stream),
                      "csv": str(CSV_PATH), "manifest": str(MANIFEST_PATH)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
