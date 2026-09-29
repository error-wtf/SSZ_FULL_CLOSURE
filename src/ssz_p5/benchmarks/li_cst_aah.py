"""Curved-spacetime Aubry--Andre--Harper benchmark from Li et al., PRB 108, 094209.

The model combines curvature-modulated hopping
  J_j = J (j/(N-1))^sigma
with the quasiperiodic onsite potential
  V_j = lambda cos(2*pi*phi*j + theta).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh_tridiagonal


PHI_GOLDEN = (np.sqrt(5.0) - 1.0) / 2.0


@dataclass(frozen=True)
class LiResult:
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    ipr: np.ndarray
    critical_site: int


def critical_site(n: int, lam: float, J: float = 1.0, sigma: float = 1.0) -> int:
    if n < 2 or J <= 0 or sigma <= 0 or lam < 0:
        raise ValueError("invalid CST-AAH parameters")
    return int(np.floor((lam / (2.0 * J)) ** (1.0 / sigma) * (n - 1)))


def solve(
    n: int = 2584,
    lam: float = 1.5,
    J: float = 1.0,
    sigma: float = 1.0,
    theta: float = 0.0,
    phi: float = PHI_GOLDEN,
) -> LiResult:
    j = np.arange(1, n + 1, dtype=float)
    diag = lam * np.cos(2.0 * np.pi * phi * j + theta)
    bonds = np.arange(1, n, dtype=float)
    off = J * (bonds / (n - 1.0)) ** sigma
    vals, vecs = eigh_tridiagonal(diag, off, check_finite=True)
    ipr = np.sum(np.abs(vecs) ** 4, axis=0)
    return LiResult(vals, vecs, ipr, critical_site(n, lam, J, sigma))


def sector_probabilities(result: LiResult) -> tuple[np.ndarray, np.ndarray]:
    """Probability of every eigenstate on the analytic left/right sectors."""
    cut = max(0, min(result.eigenvectors.shape[0], result.critical_site))
    p_left = np.sum(np.abs(result.eigenvectors[:cut]) ** 2, axis=0)
    p_right = 1.0 - p_left
    return p_left, p_right


def classify_states(result: LiResult, threshold: float = 0.9) -> dict[str, int]:
    left, right = sector_probabilities(result)
    return {
        "localized_sector": int(np.sum(left >= threshold)),
        "extended_sector": int(np.sum(right >= threshold)),
        "mixed": int(np.sum((left < threshold) & (right < threshold))),
    }
