"""MODE_CERTIFICATE_V1: per-mode certification axes for local SSZ modes.

A mode earns CERTIFIED_V1 only if it survives ALL of the axes below.
"Perfect mode" means: numerically invariant + physically admissible +
observable — never Im(omega)=0.

Axes (V1, local zone; global-QNM axes are defined but require the
currently blocked global certificate):

  A1  positive kinetic norm        : modal K-norm > 0 at every sample
  A2  tracking continuity          : min metric overlap >= threshold
  A3  eigenvalue-gap integrity     : tracked mode is non-degenerate OR
                                     degeneracy handled cluster-aware
  A4  frequency invariance         : relative frequency span across the
                                     radial window below threshold
  A5  observable residue           : mode carries a non-vanishing share
                                     of at least one frozen observable
  A6  localization determinacy     : Li-style D2 classification computed
                                     and finite (LOCALIZED / EXTENDED /
                                     MULTIFRACTAL all acceptable)
  A7  stencil/grid stability       : scale-control reproducibility flag
                                     from the v3 spectroscopy run
  A8  basis invariance             : basis-resolvent control flag

Global-QNM-only axes (V2, NOT evaluated here while the K_scalar pocket
is unresolved): Jost/ECS pole stability, Im(omega) reproducibility,
boundary-condition independence.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np

CERTIFICATE_VERSION = "MODE_CERTIFICATE_V1"

# Frozen thresholds (V1).  Rationale documented next to each use.
THRESHOLDS = {
    "min_tracking_overlap": 0.90,       # metric overlap along the track
    "min_relative_eigenvalue_gap": 1e-3,  # matches v3 REL_CLUSTER_GAP
    "max_relative_frequency_span": 5e-3,  # dimensionless drift budget
    "min_observable_share": 1e-4,       # visible residue share
    "min_kinetic_norm": 0.0,            # strict positivity
}


@dataclass(frozen=True)
class ModeCertificate:
    version: str
    L: int
    mode_index: int
    axes: dict
    verdict: str  # CERTIFIED_V1 | REJECTED_V1 | NOT_EVALUABLE_V1
    failed_axes: list[str] = field(default_factory=list)

    def to_json(self) -> dict:
        d = asdict(self)
        # schema-conformant field name (schemas/MODE_CERTIFICATE.schema.json)
        d["certificate_version"] = d.pop("version")
        return d


def _check(name: str, ok: bool, detail: dict) -> tuple[bool, dict]:
    return ok, {name: {"pass": bool(ok), **detail}}


def certify_mode(
    *,
    L: int,
    mode_index: int,
    skip_axes: tuple[str, ...] = (),
    kinetic_norms: np.ndarray,
    tracking_overlaps: np.ndarray,
    relative_gaps: np.ndarray,
    omega_values: np.ndarray,
    observable_shares: np.ndarray,
    d2_value: float,
    scale_control_ok: bool,
    basis_control_ok: bool,
) -> ModeCertificate:
    """Evaluate all V1 axes for one tracked mode across the radial window."""
    axes: dict = {}
    failed: list[str] = []

    ok, detail = _check(
        "A1_positive_kinetic_norm",
        bool(np.all(np.asarray(kinetic_norms, float) > THRESHOLDS["min_kinetic_norm"])),
        {"min": float(np.min(kinetic_norms))},
    )
    axes.update(detail)
    if not ok:
        failed.append("A1")

    ok, detail = _check(
        "A2_tracking_continuity",
        bool(np.min(tracking_overlaps) >= THRESHOLDS["min_tracking_overlap"]),
        {"min_overlap": float(np.min(tracking_overlaps))},
    )
    axes.update(detail)
    if not ok:
        failed.append("A2")

    ok, detail = _check(
        "A3_eigenvalue_gap_integrity",
        bool(np.min(relative_gaps) > THRESHOLDS["min_relative_eigenvalue_gap"]),
        {"min_relative_gap": float(np.min(relative_gaps))},
    )
    axes.update(detail)
    if not ok:
        failed.append("A3")

    finite_omega = np.asarray(omega_values, float)
    finite_omega = finite_omega[np.isfinite(finite_omega)]
    if len(finite_omega) >= 2:
        mean_abs = max(abs(float(np.mean(finite_omega))), 1e-300)
        span = float((np.max(finite_omega) - np.min(finite_omega)) / mean_abs)
        ok = span <= THRESHOLDS["max_relative_frequency_span"]
    else:
        span = float("nan")
        ok = False
    ok, detail = _check(
        "A4_frequency_invariance",
        ok,
        {"relative_span": span, "threshold": THRESHOLDS["max_relative_frequency_span"]},
    )
    axes.update(detail)
    if not ok:
        failed.append("A4")

    max_share = float(np.max(observable_shares)) if len(observable_shares) else 0.0
    ok, detail = _check(
        "A5_observable_residue",
        max_share >= THRESHOLDS["min_observable_share"],
        {"max_share": max_share},
    )
    axes.update(detail)
    if not ok:
        failed.append("A5")

    ok, detail = _check(
        "A6_localization_determinacy",
        bool(math.isfinite(d2_value)),
        {"d2": float(d2_value)},
    )
    axes.update(detail)
    if not ok:
        failed.append("A6")

    ok, detail = _check("A7_grid_stability", bool(scale_control_ok), {})
    axes.update(detail)
    if not ok:
        failed.append("A7")

    ok, detail = _check("A8_basis_invariance", bool(basis_control_ok), {})
    axes.update(detail)
    if not ok:
        failed.append("A8")

    for ax in skip_axes:
        axes[f"{ax}_skipped"] = {"pass": None, "note": "not decidable on this input"}
    hard_failed = [f for f in failed if f not in skip_axes]
    if hard_failed:
        verdict = "REJECTED_V1"
    elif failed:
        verdict = "NOT_EVALUABLE_V1"
    else:
        verdict = "CERTIFIED_V1"
    return ModeCertificate(
        version=CERTIFICATE_VERSION,
        L=L,
        mode_index=mode_index,
        axes=axes,
        verdict=verdict,
        failed_axes=hard_failed,
    )
