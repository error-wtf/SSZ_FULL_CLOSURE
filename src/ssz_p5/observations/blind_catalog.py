"""Blind OBSERVED-mode catalog builder (contract rule 5).

Takes the Stage-A PSD candidates plus analysis config and instrumental
control outcomes and produces a document that CONFORMS to
schemas/OBSERVED_MODE_CATALOG.schema.json.

frozen=true + sha256 are reachable ONLY through freeze_observed_catalog
(freeze is a contract act).  verify_frozen_catalog round-trips the hash.
"""
from __future__ import annotations

import hashlib
import json


def build_observed_catalog(
    *,
    target: str,
    obsids: list[str],
    total_exposure_s: float,
    analysis_config: dict,
    candidates: list[dict],
    instrumental_controls: list[dict],
    pipeline_valid: bool = True,
) -> dict:
    """candidates: dicts with f_hz, sigma_f_hz, amplitude_rms,
    covariance, significance_sigma, energy_dependence, obs_occurrence
    (plus optional gamma_fwhm_hz, notes).  instrumental_controls: dicts
    with f_hz, predicted_by, verdict, rejection_reason."""
    modes = []
    for c in candidates:
        m = {
            "f_hz": float(c["f_hz"]),
            "sigma_f_hz": float(c["sigma_f_hz"]),
            "amplitude_rms": float(c["amplitude_rms"]),
            "covariance": c.get("covariance", []),
            "significance_sigma": float(c["significance_sigma"]),
            "energy_dependence": c.get("energy_dependence", "none"),
            "obs_occurrence": int(c["obs_occurrence"]),
        }
        for opt in ("gamma_fwhm_hz", "notes"):
            if c.get(opt) is not None:
                m[opt] = c[opt]
        modes.append(m)
    return {
        "catalog": "OBSERVED_MODE_CATALOG",
        "version": 1,
        "frozen": False,
        "sha256": "",
        "target": target,
        "obsids": list(obsids),
        "total_exposure_s": float(total_exposure_s),
        "analysis_config": dict(analysis_config),
        "modes": modes,
        "instrumental_controls": [dict(ic) for ic in instrumental_controls],
        "pipeline_valid": pipeline_valid,
    }


def freeze_observed_catalog(doc: dict) -> tuple[dict, str]:
    """Set frozen=true + sha256 over the canonical freeze payload."""
    if doc.get("frozen") is not False:
        raise ValueError("refusing to freeze: catalog is not in unfrozen state")
    frozen = dict(doc)
    frozen["frozen"] = True
    frozen.pop("sha256", None)
    payload = json.dumps(frozen, sort_keys=True, allow_nan=False).encode()
    digest = hashlib.sha256(payload).hexdigest()
    frozen["sha256"] = digest
    return frozen, digest


def verify_frozen_catalog(doc: dict) -> bool:
    if doc.get("frozen") is not True:
        return False
    check = dict(doc)
    claimed = check.pop("sha256", "")
    payload = json.dumps(check, sort_keys=True, allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest() == claimed
