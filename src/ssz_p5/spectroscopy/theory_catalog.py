"""Blind predicted-mode catalog builder (contract rule 6).

Consumes ONLY frozen local-operator outputs; refuses to read anything
under data/observed/ (blindness asserted at open time).

The output conforms to schemas/PREDICTED_MODE_CATALOG.schema.json:
frozen=true requires an explicit freeze step (write_frozen_from_catalog)
that computes sha256 and writes the sidecar — a catalog without that
step is not freezable and not comparable.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .physical_units import FrozenMassPrior, unit_conversion_factor_hz

_STATE_CLASS = {
    "EXTENDED": "extended",
    "LOCALIZED": "localized",
    "MULTIFRACTAL": "multifractal",
    "UNDETERMINED": "undetermined",
}


def assert_blind(root: Path) -> None:
    """Raise if any observed catalog is readable — the theory side must not
    see observed results before the match step."""
    for probe in (
        root / "data/observed/OBSERVED_MODE_CATALOG.json",
        root / "data/observed/SMOKE_5200120403_STAGE_A.json",
    ):
        if probe.exists():
            raise RuntimeError(
                f"BLINDNESS VIOLATION: theory catalog builder found {probe}")


def build_predicted_catalog(
    per_mode_rows: list[dict],
    prior: FrozenMassPrior,
    *,
    member_hash: str,
    healthy_window: dict,
    observables: list[str] | None = None,
) -> dict:
    """per_mode_rows: dicts with omega_bar_re, observable name, u, omega2,
    cluster_residue_weight, ipr, localization_class.  The returned
    document is schema-conformant EXCEPT frozen/sha256, which only the
    explicit freeze step may set (frozen=true is a contract act, not a
    builder default)."""
    conv = unit_conversion_factor_hz(prior.mass_solar)
    modes = []
    for row in per_mode_rows:
        cls = str(row.get("localization_class", "UNDETERMINED")).upper()
        modes.append({
            "observable": str(row["observable"]),
            "u": float(row["u"]),
            "L": int(row.get("L", 0)),
            "omega2": float(row["omega2"]),
            "cluster_residue_weight": float(row.get("cluster_residue_weight", 0.0)),
            "f_hz": conv * float(row["omega_bar_re"]),
            "omega_bar_re": float(row["omega_bar_re"]),
            "ipr": float(row["ipr"]) if row.get("ipr") is not None else None,
            "scaling_dimension": row.get("scaling_dimension"),
            "state_class": _STATE_CLASS.get(cls, "undetermined"),
        })
    return {
        "catalog": "PREDICTED_MODE_CATALOG",
        "version": 1,
        "frozen": False,  # set true only by freeze_predicted_catalog
        "sha256": "",
        "member_hash": member_hash,
        "mass_prior": prior.to_json(),
        "unit_conversion_hz_per_omega_bar": conv,
        "healthy_window": healthy_window,
        "observables": observables or [],
        "modes": modes,
        "global_qnm_included": False,
        "blindness_assert": {
            "did_not_read_observed_catalog": True,
            "enforced_by": (
                "ssz_p5.spectroscopy.theory_catalog.assert_blind at build time"),
        },
    }


def freeze_predicted_catalog(doc: dict) -> tuple[dict, str]:
    """Set frozen=true + sha256 over the canonical freeze payload.  Returns
    (frozen_doc, sha256_hex)."""
    if doc.get("frozen") is not False:
        raise ValueError("refusing to freeze: catalog is not in unfrozen state")
    frozen = dict(doc)
    frozen["frozen"] = True
    frozen.pop("sha256", None)   # hash payload WITHOUT the hash field
    payload = json.dumps(
        frozen, sort_keys=True, allow_nan=False).encode()
    digest = hashlib.sha256(payload).hexdigest()
    frozen["sha256"] = digest
    return frozen, digest


def verify_frozen_catalog(doc: dict) -> bool:
    """Recompute the freeze hash; frozen catalogs must verify or be rejected."""
    if doc.get("frozen") is not True:
        return False
    check = dict(doc)
    claimed = check.pop("sha256", "")
    payload = json.dumps(check, sort_keys=True, allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest() == claimed
