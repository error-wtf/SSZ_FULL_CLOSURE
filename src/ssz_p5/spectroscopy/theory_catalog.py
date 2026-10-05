"""Blind predicted-mode catalog builder (contract rule 6).

Consumes ONLY frozen local-operator outputs; refuses to read anything
under data/observed/ (blindness asserted at open time).
"""
from __future__ import annotations

from pathlib import Path

from .physical_units import FrozenMassPrior, unit_conversion_factor_hz


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
) -> dict:
    """per_mode_rows: dicts with omega_bar_re, Z_q, ipr, localization_class,
    convergence dict.  Frequencies in Hz via the FROZEN prior only."""
    conv = unit_conversion_factor_hz(prior.mass_solar)
    modes = []
    for i, row in enumerate(per_mode_rows):
        modes.append({
            "mode_index": i,
            "f_hz": conv * float(row["omega_bar_re"]),
            "omega_bar_re": float(row["omega_bar_re"]),
            "Z_q": row.get("Z_q"),
            "ipr": row.get("ipr"),
            "localization_class": row.get("localization_class"),
            "convergence": row.get("convergence", {}),
        })
    return {
        "catalog_type": "PREDICTED_MODE_CATALOG",
        "mass_prior": prior.to_json(),
        "unit_conversion_hz_per_omega_bar": conv,
        "n_modes": len(modes),
        "modes": modes,
    }
