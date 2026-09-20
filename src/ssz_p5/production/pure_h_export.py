"""Direct weak/core Maxwell--Horndeski 41-slot production exports.

These exports do not use the blacklisted historical global selected stream as
an input.  They start from the action-derived F2 regional tables, apply the
locked lower-order representative v5=c3=e3=0, v12=0 on A0'=0 pure-H regions,
and regenerate the accepted holonomic a5 through the frozen 16-Sep identity.
The historical selected stream is used only as an independent regression
oracle.
"""
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import SLOT_NAMES
from ..numerics import module
from ..provenance.manifest import sha256
from .member import MEMBER_FILE, validate_stream_regions

WEAK_SOURCE = "data/regression/ssz_p5_F2_exterior_horndeski_unreduced_39of41_2026-09-14.csv"
CORE_SOURCE = "data/authoritative/ssz_p5_F2_core_punctured_horndeski_unreduced_39of41_2026-09-14.csv"
REGRESSION_STREAM = "archive/full_working_snapshot/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv"
A5_IMPL = "src/ssz_p5_holonomic_a5_closure_2026-09-16.py"


def _scaled_max(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.max(np.abs(a-b)/np.maximum(1.0, np.abs(b))))


def _a5_selected(source: pd.DataFrame, smoothing: float) -> pd.DataFrame:
    # Reuse the frozen accepted implementation; do not fork the identity.
    return module("ssz_p5_holonomic_a5_closure_2026-09-16.py").select(source, smoothing)


def _normalize_pure_h(source: pd.DataFrame, *, region: str, smoothing: float,
                      u_mask) -> pd.DataFrame:
    d = source.copy().reset_index(drop=True)
    if "phiprime" not in d and "phi_r" in d:
        d["phiprime"] = d["phi_r"]
    if "A0prime" not in d:
        d["A0prime"] = 0.0
    if np.max(np.abs(d.A0prime.to_numpy(float))) > 1e-13:
        raise ValueError(f"{region}: pure-H direct export requires A0'=0")

    a5 = _a5_selected(d, smoothing).sort_values("x").reset_index(drop=True)
    s = d.sort_values("x").reset_index(drop=True)
    np.testing.assert_allclose(s.x, a5.x, rtol=0, atol=1e-13)
    s["a5"] = a5.a5_selected.to_numpy(float)
    s["v5"] = 0.0
    s["c3"] = 0.0
    s["e3"] = 0.0
    s["v12"] = 0.0
    s = s.loc[u_mask(s.u.to_numpy(float))].copy().reset_index(drop=True)
    s["region"] = region
    s["source_member"] = MEMBER_FILE
    if not np.isfinite(s[list(SLOT_NAMES)].to_numpy(float)).all():
        raise ValueError(f"{region}: non-finite 41-slot direct export")
    validate_stream_regions(s)
    return s


def _compare_to_regression(root: Path, actual: pd.DataFrame, region: str):
    # Regression only. This stream is explicitly blacklisted as a final global
    # source by the execution contract and is never used to construct actual.
    old = pd.read_csv(root / REGRESSION_STREAM)
    old = old.loc[old.region == region].sort_values("x").reset_index(drop=True)
    act = actual.sort_values("x").reset_index(drop=True)
    if len(old) != len(act):
        raise ValueError(f"{region}: regression row count mismatch {len(act)} != {len(old)}")
    np.testing.assert_allclose(act.x, old.x, rtol=0, atol=1e-13)
    rows=[]
    for name in SLOT_NAMES:
        err=_scaled_max(act[name],old[name])
        rows.append(dict(slot=name,max_scaled_error=err,status="PASS" if err < 1e-10 else "FAIL"))
    return pd.DataFrame(rows)


def _certificate(root: Path, output: Path, *, gate: str, region: str, source_path: str,
                 stream_path: Path, comparison_path: Path, comparison: pd.DataFrame,
                 a5_smoothing: float):
    passed = bool((comparison.status == "PASS").all())
    cert = {
        gate: "PASS" if passed else "FAIL",
        "region": region,
        "scope": "direct action-derived pure-H regional 41-slot export",
        "construction": "F2 action-derived region + v5=c3=e3=0 + v12=0 + frozen holonomic a5",
        "historical_global_stream_used_as_input": False,
        "historical_stream_role": "regression oracle only",
        "a5_smoothing": a5_smoothing,
        "slots": list(SLOT_NAMES),
        "rows": int(len(pd.read_csv(stream_path))),
        "max_regression_scaled_error": float(comparison.max_scaled_error.max()),
        "timestamp": datetime.now(UTC).isoformat(),
        "source": {"path": source_path, "sha256": sha256(root/source_path)},
        "implementation": [
            {"path": "src/ssz_p5/production/pure_h_export.py", "sha256": sha256(root/"src/ssz_p5/production/pure_h_export.py")},
            {"path": A5_IMPL, "sha256": sha256(root/A5_IMPL)},
            {"path": MEMBER_FILE, "sha256": sha256(root/MEMBER_FILE)},
        ],
        "outputs": [
            {"path": str(stream_path.relative_to(root)), "sha256": sha256(stream_path)},
            {"path": str(comparison_path.relative_to(root)), "sha256": sha256(comparison_path)},
        ],
        "absolute_full_closure": "NOT_CERTIFIED",
    }
    return cert


def export_pure_h(root: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    weak_src = pd.read_csv(root/WEAK_SOURCE)
    core_src = pd.read_csv(root/CORE_SOURCE)

    weak = _normalize_pure_h(
        weak_src, region="weak_exterior_H", smoothing=0.0,
        u_mask=lambda u: u < 0.5515230871346237,
    )
    core = _normalize_pure_h(
        core_src, region="punctured_H_core", smoothing=1e-8,
        u_mask=lambda u: u >= 0.715,
    )

    weak_path=output/"weak_exterior_H_DIRECT_41.csv"
    core_path=output/"punctured_H_core_DIRECT_41.csv"
    weak.to_csv(weak_path,index=False); core.to_csv(core_path,index=False)

    wc=_compare_to_regression(root,weak,"weak_exterior_H")
    cc=_compare_to_regression(root,core,"punctured_H_core")
    wcp=output/"WEAK_DIRECT_41_COMPARISON.csv"; ccp=output/"CORE_DIRECT_41_COMPARISON.csv"
    wc.to_csv(wcp,index=False); cc.to_csv(ccp,index=False)

    wcert=_certificate(root,output,gate="WEAK_DIRECT_41",region="weak_exterior_H",
                       source_path=WEAK_SOURCE,stream_path=weak_path,comparison_path=wcp,
                       comparison=wc,a5_smoothing=0.0)
    ccert=_certificate(root,output,gate="CORE_DIRECT_41",region="punctured_H_core",
                       source_path=CORE_SOURCE,stream_path=core_path,comparison_path=ccp,
                       comparison=cc,a5_smoothing=1e-8)
    (output/"WEAK_DIRECT_41_CERTIFICATE.json").write_text(json.dumps(wcert,indent=2)+"\n")
    (output/"CORE_DIRECT_41_CERTIFICATE.json").write_text(json.dumps(ccert,indent=2)+"\n")
    return wcert,ccert
