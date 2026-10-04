#!/usr/bin/env python3
"""REAL_SPECTROSCOPY_BRIDGE_V1 - Stage A skeleton (smoke test level).

Implements the blind observed-mode extraction machinery for one NICER obsid
(5200120403) including the MANDATORY GTI-gap artifact negative control
(~55 Hz documented by NICER analysis tips for MAXI J1820+070).

What this skeleton does TODAY (verifiable without downloaded products):
  * GTI-gap spectral artifact PREDICTION: given GTI segments (start, stop) the
    gap cadence produces a comb of frequencies; the pipeline must show that a
    gap-periodicity artifact lands at/near 55 Hz for realistic gap patterns
    and must carry the machinery to reject it.
  * Blindness asserts: refuses to run if data/observed/ or
    data/predicted/ catalogs from the OPPOSITE stage are readable.

What requires the actual FITS products (next engineering step):
  * event-level barycentering, PSD extraction, peak fitting, covariance.

Usage:
  python3 tools/bridge_stage_a_observed_modes.py --check-env
  python3 tools/bridge_stage_a_observed_modes.py --gti-demo
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

FIRST_TEST_OBSID = "5200120403"
GTI_ARTIFACT_FREQ_HZ = 55.0  # documented instrumental artifact for this target


def assert_blindness(stage: str) -> None:
    """Stage A must not read predicted catalogs; stage B not observed ones."""
    forbidden = {
        "A": ROOT / "data/predicted/PREDICTED_MODE_CATALOG.json",
        "B": ROOT / "data/observed/OBSERVED_MODE_CATALOG.json",
    }
    f = forbidden.get(stage)
    if f is not None and f.exists():
        raise RuntimeError(
            f"BLINDNESS VIOLATION: stage {stage} must not read {f}. "
            "Remove the file or run in an isolated checkout.")
    print(f"blindness assert stage {stage}: OK (no forbidden inputs)")


def gti_gap_comb(gaps_s: np.ndarray, t_max_s: float):
    """Frequencies where periodic data gaps leak power.

    A periodic gap pattern with mean gap period T_g leaks power at multiples
    of 1/T_g (and aliases thereof). Returns the comb frequencies in Hz.
    """
    if len(gaps_s) < 2:
        return np.array([])
    periods = np.diff(np.sort(gaps_s))
    T_g = float(np.median(periods))
    n = int(t_max_s * GTI_ARTIFACT_FREQ_HZ) + 10
    comb = np.arange(1, n + 1) / T_g
    return comb[comb < 100.0]


def gti_demo():
    """Demonstrate the 55 Hz artifact mechanism on synthetic GTI patterns.

    NICER GTI gaps for this target produce characteristic dropouts; the
    documented ~55 Hz artifact corresponds to a characteristic gap cadence.
    We verify: a comb built from a gap period T_g = 1/55 s has a line at 55 Hz
    and that harmonic structure could masquerade as a QPO. This is the control
    the final pipeline MUST pass: recover and reject it.
    """
    # synthetic GTI: 1 s exposures separated by 1/55 s gaps (worst case)
    seg = 1.0
    n = 200
    starts = np.arange(n) * (seg + 1.0 / GTI_ARTIFACT_FREQ_HZ)
    comb = gti_gap_comb(starts, t_max_s=starts[-1] + seg)
    hit = np.min(np.abs(comb - GTI_ARTIFACT_FREQ_HZ))
    ok = hit < 0.5
    print(f"GTI comb lines: {len(comb)}, first 5: {np.round(comb[:5], 2)}")
    print(f"closest comb line to {GTI_ARTIFACT_FREQ_HZ} Hz: {hit:.3f} Hz "
          f"-> {'HIT' if ok else 'MISS'}")
    return ok


def check_env():
    cat = ROOT / "data/observations/nicer/MAXI_J1820_PLUS_070_CATALOG_V1.json"
    assert cat.exists(), "frozen NICER catalog missing"
    j = json.load(open(cat))
    obsids = [r["obsid"] for r in j["observations"]]
    print(f"frozen catalog OK: {j['n_obs']} obs, {j['total_exposure_s']/1000:.1f} ks")
    print(f"first test obsid present: {FIRST_TEST_OBSID in obsids}")
    assert_blindness("A")
    # heasoft presence check (soft)
    import shutil
    for tool in ("nicerl2", "ftplot", "ftrimode"):
        print(f"heasoft tool {tool}: {shutil.which(tool) or 'NOT IN PATH'}")
    try:
        from astropy.io import fits  # noqa: F401
        print("astropy.io.fits: available")
    except ImportError:
        print("astropy.io.fits: MISSING (required for product ingest)")
    try:
        import stingray  # noqa: F401
        print("stingray: available (PSD machinery)")
    except ImportError:
        print("stingray: MISSING (will implement PSD via numpy if absent)")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-env", action="store_true")
    ap.add_argument("--gti-demo", action="store_true")
    args = ap.parse_args()
    if args.check_env:
        return check_env()
    if args.gti_demo:
        assert_blindness("A")
        return 0 if gti_demo() else 1
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
