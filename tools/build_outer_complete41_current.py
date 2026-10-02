#!/usr/bin/env python3
"""Rebuild the complete outer H+SVT 41-slot coefficient stream.

This implements the documented common-action assembly used in the 2026-09-16
outer handover work:

    C_total = C_MH + C_ZK - C_shared

on the unreduced 41-slot action coefficients, followed by the locked lower
selection/recanonicalization.  The SVT part is re-emitted from the authoritative
re-solved outer background plus the documented S-weighted transverse Hessian
closure.  The MH component is the archived direct-action MH cubic stream used
by the original outer construction.  The shared Einstein/Maxwell/scalar
baseline is re-emitted on the identical outer grid and subtracted exactly once.

Scope: this closes the *coefficient-level complete outer assembly/pivot gap*.
It does not by itself prove that the resulting finite-l kinetic operator is
healthy; that is audited immediately downstream.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.config import SLOT_NAMES
from ssz_p5.numerics import module
from ssz_p5.production.regional_coefficients import select_lower

OUTDIR = ROOT / "data/generated/qnm_global_diagnostic"
OUT = OUTDIR / "OUTER_COMPLETE_H_SVT_41.csv"
REPORT = OUTDIR / "OUTER_COMPLETE_H_SVT_41_AUDIT.json"

BG = ROOT / "data/production/ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv"
SVT_CTRL = ROOT / "archive/full_working_snapshot/ssz_p5_outer_svt_repaired_integrable_transition_profile_2026-09-12.csv"
MH41 = ROOT / "archive/full_working_snapshot/ssz_p5_OUTER_MH_CUBIC_SELECTED_41of41_2026-09-16.csv"

HESS = dict(HXX="f2XX", HXF="f2XF", HXY="f2XY", HFF="f2FF", HFY="f2FY", HYY="f2YY")


def digest(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def interp(src: pd.DataFrame, u: np.ndarray, c: str) -> np.ndarray:
    return np.interp(u, src.u.to_numpy(float), src[c].to_numpy(float))


def emit_shared(zk, base: pd.DataFrame, f2F: np.ndarray) -> pd.DataFrame:
    d = base.copy()
    for c in ("f2X","f3","f3X","f3XX","tf3","f4","f4X","f4XX","f4XXX","tf4"):
        d[c] = 0.0
    d["f2F"] = f2F
    d["f2Y"] = 0.0
    for c in HESS.values():
        d[c] = 0.0
    return zk.emit(d, selected_v5=0.0, selected_c3=0.0, selected_e3=0.0)


def main() -> int:
    bg = pd.read_csv(BG).sort_values("u").reset_index(drop=True)
    src = pd.read_csv(SVT_CTRL).sort_values("u").reset_index(drop=True)
    mh = pd.read_csv(MH41).sort_values("u").reset_index(drop=True)
    zk = module("ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py")

    u = bg.u.to_numpy(float)
    r = bg.x.to_numpy(float)
    f = bg.f.to_numpy(float)
    h = bg.h.to_numpy(float)
    X = bg.X.to_numpy(float)
    ph = -np.sqrt(np.maximum(0.0, -2.0 * X / h))

    if len(mh) != len(bg) or np.max(np.abs(mh.u.to_numpy(float)-u)) > 1e-12:
        # historical MH direct stream uses same outer construction grid; interpolate
        # only if a serialization/order difference exists.
        mhi = {"u":u}
        for c in mh.columns:
            if c == "u":
                continue
            if np.issubdtype(mh[c].dtype, np.number):
                mhi[c] = np.interp(u, mh.u.to_numpy(float), mh[c].to_numpy(float))
        mh = pd.DataFrame(mhi)

    base = pd.DataFrame(dict(
        u=u, x=r, phi=bg.phi.to_numpy(float), f=f, h=h,
        phiprime=ph, A0prime=bg.A0prime.to_numpy(float), X=X,
        f2X=bg.f2X.to_numpy(float), f2F=bg.f2F.to_numpy(float), f2Y=0.0,
        f3=bg.f3.to_numpy(float), f3X=bg.f3X.to_numpy(float),
        f3XX=0.0, tf3=0.0,
        f4=bg.f4.to_numpy(float), f4X=bg.N4.to_numpy(float),
        f4XX=0.0, f4XXX=0.0, tf4=0.0,
    ))
    S = bg.S_SVT.to_numpy(float)
    for old, new in HESS.items():
        base[new] = S * interp(src, u, old)
    base["f2YY"] = S * interp(src, u, "HYY")

    svt = zk.emit(base, selected_v5=0.0, selected_c3=0.0, selected_e3=0.0)

    # Shared baseline used by the original construction:
    # Einstein/scalar baseline plus S-weighted Maxwell f2F=1 increment.
    sh0 = emit_shared(zk, base, np.zeros_like(u))
    sh1 = emit_shared(zk, base, np.ones_like(u))
    shared = sh0.copy()
    for k in SLOT_NAMES:
        shared[k] = sh0[k].to_numpy(float) + S * (
            sh1[k].to_numpy(float) - sh0[k].to_numpy(float)
        )

    geom_cols = ["u","x","phi","f","h","phiprime","A0prime"]
    total = pd.DataFrame({
        "u":u, "x":r, "phi":bg.phi.to_numpy(float), "f":f, "h":h,
        "phiprime":ph, "A0prime":bg.A0prime.to_numpy(float),
        "S_SVT":S, "T_H":bg.T_H.to_numpy(float),
    })
    for k in SLOT_NAMES:
        total[k] = (
            mh[k].to_numpy(float)
            + svt[k].to_numpy(float)
            - shared[k].to_numpy(float)
        )
    total = select_lower(total)
    total["production_region"] = "outer_same_action_H_SVT"
    total = total.sort_values("x").reset_index(drop=True)

    # Exact pivot/identity checks.
    DeltaV = 4*total.b1.to_numpy(float)*total.v10.to_numpy(float) - total.v11.to_numpy(float)**2
    v7res = total.v7.to_numpy(float) - total.v2.to_numpy(float)**2/(4*total.v1.to_numpy(float))
    report = {
        "status": "OUTER_COMPLETE_H_SVT_41_ASSEMBLED",
        "scope": "complete unreduced coefficient assembly from action-derived MH + re-emitted SVT - re-emitted shared baseline",
        "rows": int(len(total)),
        "u_range": [float(total.u.min()), float(total.u.max())],
        "source_hashes": {
            "background": digest(BG),
            "svt_controls": digest(SVT_CTRL),
            "mh_direct41": digest(MH41),
        },
        "pivots": {
            "min_abs_v1": float(np.min(np.abs(total.v1))),
            "min_abs_v9": float(np.min(np.abs(total.v9))),
            "min_abs_v10": float(np.min(np.abs(total.v10))),
            "min_abs_DeltaV": float(np.min(np.abs(DeltaV))),
            "nearzero_v1": int(np.sum(np.abs(total.v1)<1e-14)),
            "nearzero_v9": int(np.sum(np.abs(total.v9)<1e-14)),
            "nearzero_v10": int(np.sum(np.abs(total.v10)<1e-14)),
            "nearzero_DeltaV": int(np.sum(np.abs(DeltaV)<1e-14)),
        },
        "max_abs_v7_identity_residual": float(np.max(np.abs(v7res))),
    }
    report["pivot_pass"] = all(
        report["pivots"][k] == 0
        for k in ("nearzero_v1","nearzero_v9","nearzero_v10","nearzero_DeltaV")
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    total.to_csv(OUT, index=False)
    report["output_sha256"] = digest(OUT)
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["pivot_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
