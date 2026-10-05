"""Frozen-catalog comparison: OBSERVED x PREDICTED, third-process only.

Contract rule 7: both catalogs arrive hash-frozen; this module never
re-runs theory and never re-detects peaks.  Matching is a pure function
of the two frozen JSON documents plus a fixed tolerance policy.
"""
from __future__ import annotations

import numpy as np


def match_one(observed: dict, predicted: list[dict], policy: dict) -> dict:
    """Best predicted match for one observed candidate.

    policy: {"freq_tol_rel": ..., "freq_tol_abs_hz": ..., "n_sigma": ...}
    A candidate matches when BOTH the relative and the absolute window hit
    (relative alone explodes at high f; absolute alone dies at low f).
    Returns agreement metrics; never mutates inputs.
    """
    f_obs = float(observed["f_hz"])
    sig_f = float(observed.get("sigma_f_hz") or 0.0)
    tol_rel = float(policy["freq_tol_rel"])
    tol_abs = float(policy.get("freq_tol_abs_hz", 1.0))
    n_sigma = float(policy.get("n_sigma", 3.0))
    best = None
    for pred in predicted:
        f_pred = float(pred["f_hz"])
        rel = abs(f_pred - f_obs) / max(abs(f_obs), 1e-300)
        delta = abs(f_pred - f_obs)
        within = (rel <= tol_rel) and (delta <= tol_abs)
        # combined pull if observed sigma available
        pull = delta / sig_f if sig_f > 0 else float("nan")
        if sig_f > 0:
            within = within and (pull <= n_sigma)
        score = (within, -rel)
        if best is None or score > best["_score"]:
            best = {
                "_score": score,
                "predicted_index": pred.get("mode_index"),
                "f_pred_hz": f_pred,
                "freq_rel_dev": rel,
                "within_tolerance": within,
                "pull_sigma": pull,
                "pred_zq": pred.get("Z_q"),
                "pred_ipr": pred.get("ipr"),
                "pred_class": pred.get("localization_class"),
            }
    if best is None:
        return {"matched": False, "reason": "empty predicted catalog"}
    best.pop("_score")
    best["matched"] = best["within_tolerance"]
    return best


def _is_hit(o: dict, p: dict, policy: dict) -> bool:
    """True iff o and p agree within the policy windows (pair-level test)."""
    f_o, f_p = float(o["f_hz"]), float(p["f_hz"])
    sig = float(o.get("sigma_f_hz") or 0.0)
    delta = abs(f_p - f_o)
    ok = (delta / max(abs(f_o), 1e-300) <= float(policy["freq_tol_rel"])
          and delta <= float(policy.get("freq_tol_abs_hz", 1.0)))
    if sig > 0:
        ok = ok and (delta / sig <= float(policy.get("n_sigma", 3.0)))
    return ok


def null_comparison(observed: list[dict], predicted: list[dict],
                    policy: dict, n_draws: int = 2000, seed: int = 12345) -> dict:
    """Explicit null via ASSIGNMENT permutation: each draw randomly reassigns
    predicted modes to observed candidates (without replacement), destroying
    the physical pairing while keeping both frequency sets fixed.  Empirical
    p-value: P(random assignment >= real hits)."""
    rng = np.random.default_rng(seed)
    real_hits = sum(
        1 for o in observed
        if match_one(o, predicted, policy).get("matched"))
    n_obs = min(len(observed), len(predicted))
    count_ge = 0
    for _ in range(n_draws):
        perm = rng.permutation(len(predicted))[:n_obs]
        hits = sum(
            1 for o, pi in zip(observed[:n_obs], perm, strict=False)
            if _is_hit(o, predicted[int(pi)], policy))
        if hits >= real_hits:
            count_ge += 1
    return {
        "real_hits": real_hits,
        "n_observed": len(observed),
        "n_predicted": len(predicted),
        "n_draws": n_draws,
        "p_null": (count_ge + 1) / (n_draws + 1),
        "note": ("assignment-permutation null: random predicted->observed "
                 "reassignment without replacement"),
    }
