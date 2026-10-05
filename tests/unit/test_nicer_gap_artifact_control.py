"""NICER_1200120107_SPURIOUS_55HZ_REJECTED — the mandatory negative control.

The NICER analysis tips document that short GTI gaps (telemetry saturation)
in MAXI J1820+070 ObsID 1200120107 produce a spurious ~55 Hz "QPO".
Reproduced here synthetically: gap-gap periodicity at 1/0.018 s ~ 55.6 Hz
drives comb power; continuous-GTI analysis must reject it.  A REAL signal
(control) must survive the same pipeline.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.observations.timing import (  # noqa: E402
    gap_comb_frequencies,
    spurious_peak_verdict,
)

DT = 1.0 / 512.0
F_ARTIFACT = 1.0 / 0.018  # 55.6 Hz gap periodicity


def _gapped_lightcurve(rng, total_s=64.0, gap_len=0.018, n_gaps=40):
    """Counts with dropouts every ~1/(2*55Hz) creating the comb structure."""
    n = int(total_s / DT)
    t = np.arange(n) * DT
    rate = 50.0 * np.ones(n)
    # periodic dropout pattern at the artifact period:
    period = 2.0 * gap_len
    phase = (t % period) < gap_len
    rate[phase] *= 0.1
    counts = rng.poisson(rate * DT)
    # build GTIs: continuous blocks split by the dropout windows
    gti_starts, gti_stops = [], []
    in_gti = False
    for i in range(n):
        gapped = phase[i]
        if not gapped and not in_gti:
            gti_starts.append(t[i]); in_gti = True
        elif gapped and in_gti:
            gti_stops.append(t[i]); in_gti = False
    if in_gti:
        gti_stops.append(t[-1])
    return counts, np.array(gti_starts), np.array(gti_stops)


def test_gap_comb_frequency_is_55hz():
    assert abs(gap_comb_frequencies(0.018) - 55.56) < 0.5


def test_spurious_55hz_is_rejected():
    """The documented artifact: on-comb AND unstable in continuous segments."""
    rng = np.random.default_rng(7)
    counts, gs, ge = _gapped_lightcurve(rng)
    # continuous segments: cut one long continuous GTI into 8 s segments
    segs = []
    for s0, s1 in zip(gs, ge):
        if s1 - s0 > 8.0:
            k = int((s1 - s0) / 8.0)
            for j in range(k):
                i0 = int((s0 + j * 8.0) / DT)
                i1 = int((s0 + (j + 1) * 8.0) / DT)
                segs.append(counts[i0:i1])
    verdict = spurious_peak_verdict(F_ARTIFACT, gs, ge, segs, DT)
    assert verdict["on_gap_comb"], "55 Hz candidate should sit on the gap comb"
    assert verdict["verdict"] == "SPURIOUS_GTIL_GAP_ARTIFACT_REJECTED"


def test_real_signal_survives_same_pipeline():
    """Control: a genuine 20 Hz oscillation inside continuous GTIs must NOT
    be flagged as artifact."""
    rng = np.random.default_rng(11)
    total_s = 64.0
    n = int(total_s / DT)
    t = np.arange(n) * DT
    rate = 50.0 * (1.0 + 0.5 * np.sin(2 * np.pi * 20.0 * t))
    counts = rng.poisson(rate * DT)
    # two clean GTIs with one mid-gap (unrelated to 20 Hz)
    gs = np.array([0.0, 33.0])
    ge = np.array([32.0, total_s])
    segs = [counts[: int(32.0 / DT)], counts[int(33.0 / DT):]]
    verdict = spurious_peak_verdict(20.0, gs, ge, segs, DT)
    assert not verdict["on_gap_comb"] or verdict["segment_stability"]["stable"]
    assert verdict["verdict"] in ("SURVIVES_CONTINUOUS_SEGMENTS", "INCONCLUSIVE")
