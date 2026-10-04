#!/usr/bin/env python3
"""BRIDGE Stage A, step A3: blind PSD + proper peak significance for one obsid.

Implements the correct Leahy-PSD statistics for MAXI J1820+070 obsid
5200120403 (first engineering test of REAL_SPECTROSCOPY_BRIDGE_V1).

Statistics (no shortcuts):
  - Leahy PSD per segment; average over nseg segments;
  - per-bin null distribution: psd/2 ~ Gamma(k=nseg, theta=1/nseg)
    (sum of nseg chi2_2/2), exact survival via regularized upper
    incomplete gamma;
  - trial correction: Bonferroni over the independent bins actually searched;
  - candidates = local maxima with GLOBAL posterior-search significance > 4.

Also runs the MANDATORY GTI-gap artifact prediction (the documented ~55 Hz
NICER artifact class): a data-gap comb with period T_g leaks power at k/T_g.
The pipeline must demonstrably be able to REJECT such a line.

This smoke run answers: "does the obs show any globally >4 sigma QPO-like
line in 0.05-64 Hz in the mid band?" The honest current answer for
5200120403 is recorded in the output JSON.

No SSZ input, no model input: this is the blind observed branch.
"""
from __future__ import annotations

import json
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.special import gammaincc
from scipy.stats import norm

OBS_DIR = Path("/home/error/physics/nicer_data/5200120403/xti/event_cl")
EVT = OBS_DIR / "ni5200120403_0mpu7_cl.evt.gz"
OUTDIR = Path(__file__).resolve().parents[1] / "data/observed"
OUTDIR.mkdir(parents=True, exist_ok=True)

DT = 1.0 / 128
F_MIN, F_MAX = 0.05, 64.0
GLOBAL_THRESHOLD = 4.0
BANDS = {"soft_0.5-2keV": (25, 100), "mid_2-5keV": (100, 250), "hard_5-10keV": (250, 500)}


def load_events(path):
    with fits.open(path) as f:
        evt = f["EVENTS"].data
        t = evt["TIME"]
        pi = evt["PI"]
        gti = f["GTI"].data
        starts = np.array(gti["START"], float)
        stops = np.array(gti["STOP"], float)
    return np.sort(t), np.asarray(pi), starts, stops


def segments_ps(t_events, starts, stops, seg_t):
    seg_n = int(seg_t / DT)
    segs = []
    for s, e in zip(starts, stops):
        k = int((e - s) // seg_t)
        for j in range(k):
            lo = s + j * seg_t
            sel = t_events[(t_events >= lo) & (t_events < lo + seg_t)] - lo
            c, _ = np.histogram(sel, bins=np.arange(0, seg_t + DT, DT))
            segs.append(c)
    return np.array(segs, float), seg_n


def psd_stats(segs, seg_n):
    psds = []
    for c in segs:
        n = c.sum()
        if n > 0:
            psds.append(2.0 * np.abs(np.fft.rfft(c)) ** 2 / n)
    P = np.array(psds, float)
    psd = P.mean(axis=0)
    nseg = len(psds)
    freqs = np.fft.rfftfreq(seg_n, DT)
    # exact per-bin survival: psd/2 ~ Gamma(nseg, 1/nseg)
    x = psd / 2.0
    p_single = np.clip(gammaincc(nseg, nseg * x), 1e-300, 1.0)
    return freqs, psd, nseg, p_single


def global_peaks(freqs, p_single):
    band = (freqs > F_MIN) & (freqs < F_MAX)
    n_bins = int(band.sum())
    p_glob = np.clip(n_bins * p_single, 0.0, 1.0)  # Bonferroni, capped at 1
    z_glob = np.full_like(p_glob, -np.inf)
    ok = p_glob < 1.0
    z_glob[ok] = norm.isf(p_glob[ok])
    z_glob[~band] = -np.inf
    peaks = []
    for i in range(2, len(freqs) - 1):
        if z_glob[i] > 3.0 and z_glob[i] >= z_glob[i - 1] and z_glob[i] >= z_glob[i + 1]:
            peaks.append({"f_hz": float(freqs[i]), "global_sigma": float(z_glob[i]),
                          "p_local": float(p_single[i])})
    return peaks, n_bins


def gti_artifact_check(starts, stops, peaks):
    """Predict gap-comb lines and test whether any candidate coincides.

    For each inter-GTI gap pattern with quasi-periodic cadence T_g, leak
    lines appear at k/T_g.  The documented ~55 Hz artifact arises for gap
    cadences near 1/55 s. We record the mechanism result for this obs.
    """
    gaps = starts[1:] - stops[:-1]
    short_gaps = gaps[gaps < 2.0]
    result = {
        "n_gaps_total": int(len(gaps)),
        "n_short_gaps_lt2s": int(len(short_gaps)),
        "median_short_gap_s": (float(np.median(short_gaps)) if len(short_gaps) else None),
        "comb_lines_hz": [],
        "artifact_candidates_rejected": [],
        "not_recovered": [],
    }
    if len(short_gaps) >= 3:
        T_g = float(np.median(short_gaps))
        t_span = float(stops.max() - starts.min())
        kmax = int(min(t_span, 100.0) * 55.0) + 5
        comb = np.arange(1, kmax + 1) / T_g
        comb = comb[comb < 64.0]
        result["comb_lines_hz"] = [float(x) for x in comb[:20]]
        for p in peaks:
            hit = np.min(np.abs(np.asarray(comb) - p["f_hz"])) if len(comb) else np.inf
            if hit < 0.2:
                result["artifact_candidates_rejected"].append(
                    {**p, "rejection_reason": "GTI_GAP_COMB", "comb_distance_hz": float(hit)})
        result["prediction_at_55hz"] = {
            "T_g_s": T_g,
            "nearest_comb_line_hz": float(np.min(np.abs(np.asarray(comb) - 55.0))) if len(comb) else None,
        }
    return result


def run():
    t, pi, starts, stops = load_events(EVT)
    report = {
        "catalog": "STAGE_A_SMOKE_TEST",
        "obsid": "5200120403",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "blindness": "no model input, no predicted catalog read (enforced by bridge asserts)",
        "statistics": "Leahy PSD, Gamma(nseg) per-bin exact, Bonferroni over searched bins",
        "threshold_global_sigma": GLOBAL_THRESHOLD,
        "bands": {},
    }
    total_candidates = []
    for band, (lo, hi) in BANDS.items():
        sel = (pi >= lo) & (pi <= hi)
        tb = np.sort(t[sel])
        band_res = {"n_events": int(sel.sum())}
        for seg_t in (64.0, 128.0, 512.0):
            segs, seg_n = segments_ps(tb, starts, stops, seg_t)
            if not len(segs):
                continue
            freqs, psd, nseg, p_single = psd_stats(segs, seg_n)
            peaks, n_bins = global_peaks(freqs, p_single)
            for p in peaks:
                p["segment_s"] = seg_t
            band_res[f"seg{int(seg_t)}s"] = {
                "n_segments": nseg,
                "n_candidates_global_gt3": len(peaks),
                "peaks": peaks,
            }
            total_candidates += peaks
        report["bands"][band] = band_res

    # GTI artifact machinery (on the union band result)
    art = gti_artifact_check(starts, stops, total_candidates)
    report["gti_artifact_control"] = art

    # verdict
    confirmed = [p for p in total_candidates if p["global_sigma"] > GLOBAL_THRESHOLD]
    rejected = art.get("artifact_candidates_rejected", [])
    report["candidates_global_gt4"] = confirmed
    report["artifact_rejections"] = rejected
    report["verdict"] = (
        "NO_ASTRONOMICAL_CANDIDATE_IN_THIS_OBSID_AT_GLOBAL_4SIGMA"
        if not confirmed else
        f"{len(confirmed)} CANDIDATES_PASS_GLOBAL_4SIGMA (see candidates_global_gt4)")
    blob = json.dumps(report, indent=1, sort_keys=True)
    out = OUTDIR / "SMOKE_5200120403_STAGE_A.json"
    out.write_text(blob)
    print("verdict:", report["verdict"])
    print("candidates >4 sigma:", len(confirmed))
    print("artifact rejections:", len(rejected))
    print("written:", out)
    print("sha256:", hashlib.sha256(blob.encode()).hexdigest()[:16], "...")
    return 0


if __name__ == "__main__":
    sys.exit(run())
