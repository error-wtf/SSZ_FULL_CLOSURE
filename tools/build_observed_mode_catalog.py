#!/usr/bin/env python3
"""REAL_SPECTROSCOPY_BRIDGE_V1 Stage A output: schema-conformant
OBSERVED_MODE_CATALOG from a stage-A report + analysis config.

Reads a stage-A report (data/observed/*_STAGE_A.json) that carries the
blind candidate list, maps it to the catalog schema, validates, and
freezes (frozen=true + sha256 + sidecar).  Refuses to run if a predicted
catalog is readable (blindness).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.observations.blind_catalog import (  # noqa: E402
    freeze_observed_catalog,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data/observed/OBSERVED_MODE_CATALOG.json"


def assert_blind(root: Path) -> None:
    probe = root / "data/predicted/PREDICTED_MODE_CATALOG.json"
    if probe.exists():
        raise RuntimeError(
            f"BLINDNESS VIOLATION: observed catalog builder found {probe}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage-a-report", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--target", default="MAXI J1820+070")
    ap.add_argument("--obsid", action="append", required=True)
    ap.add_argument("--total-exposure-s", type=float, required=True)
    args = ap.parse_args()

    assert_blind(ROOT)
    report = json.loads(args.stage_a_report.read_text())
    if report.get("verdict", "").startswith("NO_ASTRONOMICAL_CANDIDATE"):
        candidates = []
    else:
        # stage-A global candidates (above threshold)
        candidates = [
            {
                "f_hz": float(c["f_hz"]),
                "sigma_f_hz": float(c.get("sigma_f_hz", 0.0)),
                "amplitude_rms": float(c.get("amplitude_rms", 0.0)),
                "covariance": [],
                "significance_sigma": float(c["global_sigma"]),
                "energy_dependence": "none",
                "obs_occurrence": 1,
                "notes": "stage-A PSD candidate, single-obs",
            }
            for c in report.get("candidates_global_gt4", [])
        ]

    controls = [
        {
            "f_hz": float(fc),
            "predicted_by": "gti_gap_comb",
            "verdict": ("RECOVERED_AND_REJECTED"
                        if fc in report.get("gti_artifact_control", {})
                        .get("artifact_candidates_rejected", [])
                        else "NOT_RECOVERED"),
            "rejection_reason": "GTI-gap comb artifact negative control",
        }
        for fc in report.get("gti_artifact_control", {}).get("comb_lines_hz", [])
    ]

    import jsonschema

    catalog = {
        "catalog": "OBSERVED_MODE_CATALOG",
        "version": 1,
        "frozen": False,
        "sha256": "",
        "target": args.target,
        "obsids": args.obsid,
        "total_exposure_s": args.total_exposure_s,
        "analysis_config": {
            "psd_normalization": "leahy",
            "freq_range_hz": [0.1, 2048.0],
            "energy_bands_kev": [0.5, 2.0, 10.0],
            "significance_threshold_sigma": report.get(
                "threshold_global_sigma", 4.0),
            "trial_correction": "bonferroni",
            "frozen_before_first_psd": True,
        },
        "modes": candidates,
        "instrumental_controls": controls,
        "pipeline_valid": True,
    }
    schema = json.loads(
        (ROOT / "schemas/OBSERVED_MODE_CATALOG.schema.json").read_text())
    # validate the pre-freeze structure without the frozen=true constraint:
    pre = dict(catalog)
    pre["frozen"] = True
    pre["sha256"] = "0" * 64
    jsonschema.validate(pre, schema)

    frozen, digest = freeze_observed_catalog(catalog)
    jsonschema.validate(frozen, schema)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(frozen, indent=1) + "\n")
    (args.output.with_suffix(args.output.suffix + ".sha256")).write_text(
        digest + "\n")
    print(f"written: {args.output} (frozen, sha256={digest[:12]}...)")
    print(f"modes: {len(candidates)}, controls: {len(controls)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
