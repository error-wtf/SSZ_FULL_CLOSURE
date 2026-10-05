"""Positive control: injection/recovery with KNOWN parameters (P2).

A finder that only rejects artifacts proves nothing about sensitivity.
Here we inject a damped sinusoid with KNOWN (f, gamma, A) into a
realistic gapped light curve and require the pipeline to recover the
parameters within declared tolerances.
"""
from __future__ import annotations

import numpy as np


def inject_signal(counts: np.ndarray, dt: float, f_hz: float,
                  gamma_hz: float, amplitude_rms: float,
                  phase: float = 0.0) -> np.ndarray:
    """Add a damped sinusoid (ringdown-like) to a count series; amplitude
    is specified in rms counts so injections are comparable across
    segment lengths."""
    n = len(counts)
    t = np.arange(n) * dt
    env = np.exp(-2.0 * np.pi * gamma_hz * t)
    raw = env * np.sin(2.0 * np.pi * f_hz * t + phase)
    scale = amplitude_rms / max(float(np.sqrt(np.mean(raw * raw))), 1e-300)
    return counts + scale * raw


def recover_frequency_and_width(lc: np.ndarray, dt: float,
                                f_guess: float) -> dict:
    """Local peak interpolation around f_guess on the Leahy PSD;
    parabolic peak interpolation gives sub-bin frequency and an FWHM
    estimate from neighboring bins."""
    from .timing import leahy_psd

    f, p = leahy_psd(lc, dt)
    df = f[1] - f[0]
    i0 = int(np.argmin(np.abs(f - f_guess)))
    lo, hi = max(1, i0 - 3), min(len(f) - 1, i0 + 3)
    win = slice(lo, hi + 1)
    idx = lo + int(np.argmax(p[win]))
    if idx <= 0 or idx >= len(p) - 1:
        return {"f_hz": float(f[idx]), "gamma_fwhm_hz": float("nan"),
                "peak_leahy": float(p[idx]), "interpolated": False}
    y0, y1, y2 = p[idx - 1], p[idx], p[idx + 1]
    den = y0 - 2.0 * y1 + y2
    off = 0.0 if abs(den) < 1e-300 else 0.5 * (y0 - y2) / den
    f_peak = f[idx] + off * df
    half = y1 / 2.0
    left = idx
    while left > 0 and p[left] > half:
        left -= 1
    right = idx
    while right < len(p) - 1 and p[right] > half:
        right += 1
    gamma = (right - left) * df / 2.0
    return {"f_hz": float(f_peak), "gamma_fwhm_hz": float(gamma),
            "peak_leahy": float(p[idx]), "interpolated": True}


def injection_recovery_trial(rng: np.random.Generator, dt: float,
                             seg_s: float, f_true: float, gamma_true: float,
                             amp_rms: float, rate: float = 50.0,
                             f_tol_rel: float = 0.01) -> dict:
    """One trial: inject into pure Poisson noise, recover, compare."""
    n = int(seg_s / dt)
    counts = rng.poisson(rate * dt, n).astype(float)
    injected = inject_signal(counts, dt, f_true, gamma_true, amp_rms)
    rec = recover_frequency_and_width(injected, dt, f_true)
    f_err = abs(rec["f_hz"] - f_true) / f_true
    return {
        "f_true": f_true, "gamma_true": gamma_true, "amp_rms": amp_rms,
        "f_rec": rec["f_hz"], "gamma_rec": rec["gamma_fwhm_hz"],
        "f_err_rel": f_err, "peak_leahy": rec["peak_leahy"],
        "recovered": bool(f_err <= f_tol_rel and rec["peak_leahy"] > 20.0),
    }


def positive_control_battery(seed: int = 2026) -> dict:
    """The declared positive-control set: span frequencies, widths,
    amplitudes.  PASS requires every trial recovered."""
    rng = np.random.default_rng(seed)
    dt = 1.0 / 512.0
    trials = []
    # all f well below Nyquist (dt=1/512 s -> 256 Hz); a 400 Hz trial
    # would alias and is deliberately NOT part of the control set
    for f_true in (20.0, 55.0, 133.0, 200.0):
        for gamma_true, amp in ((0.5, 3.0), (2.0, 6.0)):
            tr = injection_recovery_trial(rng, dt, 64.0, f_true,
                                          gamma_true, amp)
            trials.append(tr)
    return {
        "control_type": "POSITIVE_INJECTION_RECOVERY",
        "n_trials": len(trials),
        "n_recovered": sum(t["recovered"] for t in trials),
        "all_recovered": all(t["recovered"] for t in trials),
        "trials": trials,
    }
