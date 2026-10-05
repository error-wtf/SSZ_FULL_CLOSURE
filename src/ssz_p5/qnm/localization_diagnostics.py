"""Li-2023-style localization diagnostics for K-normalized SSZ modes.

Reference: Li et al., PhysRevB.108.094209 (2023) — self-consistent
segmentation, IPR/fractal-dimension classification of localized vs
extended vs multifractal states.  Method-transplant only: no claim that
the AAH Hamiltonian equals the KRGSM operator.

For each K-normalized mode psi_n we form the normalized kinetic-weight
distribution

    p_{n,i} = w_i * psi_n^dag K(r_i) psi_n  /  sum_j (same),   sum_i p = 1,

and compute

    IPR_n        = sum_i p_{n,i}^2
    P2_n (inv)   = 1 / IPR                    (participation ratio)
    S_n          = -sum_i p ln p              (Shannon)
    D2_n         = -lim log(IPR) / log(N_eff) (finite-size scaling slope)

Classification follows typical wave-localization practice:
    D2 ~ 1            -> localized
    D2 ~ d (=1 here)  -> extended        (1D radial problem)
    1 < D2 < 1-ish    -> multifractal/intermediate

The scaling dimension is obtained by binning the same mode on nested
sub-grids (Li's segmentation idea): for block size b, p_b(k) = sum of p
within block k; IPR_b = sum_k p_b^2.  D2 from slope of log IPR_b vs
log b via least squares over b in [2, ~N/8].
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class ModeLocalization:
    mode_index: int
    ipr: float
    participation_ratio: float
    shannon_entropy: float
    d2_scaling: float
    classification: str
    n_grid: int


def _clean_p(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    p = p / p.sum()
    # floor to keep log finite for entropies
    floor = 1e-300
    return np.maximum(p, floor)


def mode_probability_weights(psi_n: np.ndarray, K_n: np.ndarray, w_n: np.ndarray) -> np.ndarray:
    """p_{n,i} ∝ w_i * psi† K psi, normalized to 1 (Li-style kinetic weight)."""
    density = np.einsum("ni,nij,nj->n", psi_n, K_n, psi_n)
    p = w_n * density
    total = p.sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError("non-positive modal kinetic weight")
    return p / total


def ipr_from_p(p: np.ndarray) -> float:
    p = _clean_p(p)
    return float(np.sum(p * p))


def shannon_from_p(p: np.ndarray) -> float:
    p = _clean_p(p)
    nz = p[p > 1e-300]
    return float(-np.sum(nz * np.log(nz)))


def d2_from_block_scaling(p: np.ndarray) -> float:
    """D2 via nested-block IPR scaling (Li's segmentation on 1D radial grid)."""
    p = _clean_p(p)
    n = len(p)
    bs, ys = [], []
    for b in (2, 4, 8, 16, 32, 64, 128):
        if b < 2 or b > max(4, n // 8):
            continue
        nb = n // b
        if nb < 2:
            continue
        pb = p[: nb * b].reshape(nb, b).sum(axis=1)
        ipr_b = float(np.sum(pb * pb))
        bs.append(b)
        ys.append(ipr_b)
    if len(bs) < 3:
        return float("nan")
    # log IPR_b = D2 * log b + const (IPR_b grows with block size for extended
    # states: IPR_b ~ b/n) -> D2 = +slope; localized: IPR_b ~ const -> D2 ~ 0.
    A = np.vstack([np.log(np.asarray(bs, float)), np.ones(len(bs))]).T
    slope, _ = np.linalg.lstsq(A, np.log(np.asarray(ys)), rcond=None)[0]
    return float(slope)


def classify(d2: float) -> str:
    """1D radial problem: extended ~ D2≈1, localized ~ D2→0-ish, middle = multifractal."""
    if not np.isfinite(d2):
        return "UNDETERMINED"
    if d2 > 0.75:
        return "EXTENDED"
    if d2 < 0.25:
        return "LOCALIZED"
    return "MULTIFRACTAL"


def diagnose_mode(
    mode_index: int,
    psi_n: np.ndarray,
    K_n: np.ndarray,
    w_n: np.ndarray,
) -> ModeLocalization:
    p = mode_probability_weights(psi_n, K_n, w_n)
    ipr = ipr_from_p(p)
    d2 = d2_from_block_scaling(p)
    return ModeLocalization(
        mode_index=mode_index,
        ipr=ipr,
        participation_ratio=1.0 / ipr if ipr > 0 else float("inf"),
        shannon_entropy=shannon_from_p(p),
        d2_scaling=d2,
        classification=classify(d2),
        n_grid=len(p),
    )


def diagnose_spectrum(psi: np.ndarray, K: np.ndarray, w: np.ndarray) -> list[ModeLocalization]:
    """psi: (k, n, f); K: (n, f, f); w: (n,)."""
    out = []
    for m in range(psi.shape[0]):
        out.append(diagnose_mode(m, psi[m], K, w))
    return out


def to_json(results: list[ModeLocalization]) -> dict:
    return {"modes": [asdict(r) for r in results]}
