"""Leahy PSD timing with GTI-gap artifact detection (contract rules 5+).

The 55 Hz GTI-gap artifact (NICER analysis tips, ObsID 1200120107,
MAXI J1820+070) is the canonical negative control: a peak produced by
short telemetry gaps must be REJECTED by the continuous-GTI analysis.
"""
from __future__ import annotations

import numpy as np


def leahy_psd(counts: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Leahy-normalized PSD of a (continuous) count series."""
    n = len(counts)
    if n < 4:
        raise ValueError("series too short")
    f = np.fft.rfftfreq(n, d=dt)
    ps = np.abs(np.fft.rfft(counts - counts.mean())) ** 2
    leahy = 2.0 * ps / counts.sum()
    return f[1:], leahy[1:]


def gaps_from_gtis(gti_starts: np.ndarray, gti_stops: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (gap_starts, gap_lengths) between consecutive GTIs."""
    starts = np.asarray(gti_starts, float)
    stops = np.asarray(gti_stops, float)
    gap_starts = stops[:-1]
    gap_lens = starts[1:] - stops[:-1]
    keep = gap_lens > 0
    return gap_starts[keep], gap_lens[keep]


def gap_comb_frequencies(gap_len_s: float) -> float:
    """Characteristic artifact frequency of a periodic gap pattern:
    the gap repetition drives power at 1/gap_length and harmonics."""
    if gap_len_s <= 0:
        raise ValueError("gap length must be positive")
    return 1.0 / gap_len_s


def segment_stability(f_target: float, segments: list[np.ndarray], dt: float,
                      tol_rel: float = 0.05) -> dict:
    """Does the candidate frequency appear in a majority of independent,
    continuous segments?  GTI-gap artifacts vanish when only long continuous
    stretches are analyzed."""
    hits = 0
    used = 0
    for seg in segments:
        if len(seg) < 16:
            continue
        f, p = leahy_psd(seg, dt)
        used += 1
        idx = np.argmin(np.abs(f - f_target))
        local = p[max(0, idx - 2): idx + 3]
        if p[idx] >= local.max() and p[idx] > 2.0:
            hits += 1
    if used == 0:
        return {"n_segments": 0, "hit_fraction": 0.0, "stable": False}
    frac = hits / used
    return {"n_segments": used, "n_hits": hits, "hit_fraction": frac,
            "stable": frac >= 0.6}


def spurious_peak_verdict(f_candidate: float, gti_starts: np.ndarray,
                          gti_stops: np.ndarray, segments_continuous: list[np.ndarray],
                          dt: float, comb_tol_rel: float = 0.10) -> dict:
    """Combined artifact check used by the 55 Hz negative control:
    1) does f_candidate sit on a gap-comb frequency?
    2) does it survive continuous-segment analysis?
    Rejection requires BOTH: comb-proximity AND segment instability."""
    gap_starts, gap_lens = gaps_from_gtis(gti_starts, gti_stops)
    comb_hits = []
    for gl in np.unique(np.round(gap_lens, 9)):
        f_comb = gap_comb_frequencies(float(gl))
        # check candidate and subharmonics
        for h in (1, 2, 3, 4):
            if abs(f_candidate - h * f_comb) / f_candidate <= comb_tol_rel:
                comb_hits.append({"gap_len_s": float(gl), "harmonic": h,
                                  "f_comb": f_comb})
    stab = segment_stability(f_candidate, segments_continuous, dt)
    on_comb = len(comb_hits) > 0
    rejected = on_comb and not stab["stable"]
    return {
        "f_candidate_hz": f_candidate,
        "on_gap_comb": on_comb,
        "comb_hits": comb_hits,
        "segment_stability": stab,
        "verdict": "SPURIOUS_GTIL_GAP_ARTIFACT_REJECTED" if rejected
        else ("SURVIVES_CONTINUOUS_SEGMENTS" if stab["stable"] else "INCONCLUSIVE"),
    }
