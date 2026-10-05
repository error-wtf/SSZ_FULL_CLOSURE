#!/usr/bin/env python3
"""REAL_SPECTROSCOPY_BRIDGE_V1 stage C: frozen-catalog comparison.

Third process per contract rule 7.  Refuses unfrozen inputs: both catalogs
must carry a sha256 field that is verified before matching.  Never reads
or re-runs theory; pure comparison plus permutation null.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssz_p5.spectroscopy.mode_match import match_one, null_comparison  # noqa: E402

DEFAULT_POLICY = {"freq_tol_rel": 0.05, "freq_tol_abs_hz": 1.0, "n_sigma": 3.0}


def _sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--observed", type=Path, required=True)
    ap.add_argument("--predicted", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--freq-tol-rel", type=float, default=0.05)
    ap.add_argument("--freq-tol-abs-hz", type=float, default=1.0)
    ap.add_argument("--n-sigma", type=float, default=3.0)
    args = ap.parse_args()

    obs_doc = json.loads(args.observed.read_text())
    pred_doc = json.loads(args.predicted.read_text())

    # freeze verification: sidecar .sha256 files must match
    for p in (args.observed, args.predicted):
        sidecar = p.with_suffix(p.suffix + ".sha256")
        if not sidecar.exists():
            raise SystemExit(f"FROZEN CATALOG REQUIRED: missing {sidecar}")
        expected = sidecar.read_text().split()[0].strip()
        actual = _sha256_file(p)
        if actual != expected:
            raise SystemExit(
                f"FREEZE VIOLATION: {p} sha256 {actual[:12]} != frozen {expected[:12]}")

    policy = dict(
        DEFAULT_POLICY,
        freq_tol_rel=args.freq_tol_rel,
        freq_tol_abs_hz=args.freq_tol_abs_hz,
        n_sigma=args.n_sigma,
    )
    observed = obs_doc.get("candidates") or obs_doc.get("modes") or []
    predicted = pred_doc.get("modes") or []

    matches = [dict(candidate=o, **match_one(o, predicted, policy)) for o in observed]
    null = null_comparison(observed, predicted, policy)

    report = {
        "comparison_version": "FROZEN_CATALOG_COMPARE_V1",
        "observed_source": str(args.observed),
        "predicted_source": str(args.predicted),
        "policy": policy,
        "matches": matches,
        "null_comparison": null,
        "verdict": (
            "MATCH_ABOVE_NULL" if null["p_null"] < 0.05 and null["real_hits"] > 0
            else "NO_MATCH_ABOVE_NULL"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    print(f"written: {args.output}")
    print(json.dumps({"verdict": report["verdict"], "null": null}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
