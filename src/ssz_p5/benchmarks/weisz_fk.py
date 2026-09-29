"""Frenkel--Kontorova benchmark following Weisz et al., PRB 18, 3275 (1978).

Dimensionless convention:
  spring stiffness alpha = 1
  natural chain spacing a = 1
  beta = V0 Q^2 / (2 alpha)

For an N*lambda = M*a rational approximant there are M atoms in a cell of
length L=M*a and substrate wave number Q=2*pi/lambda.

The q=0 optical residue of a normalized normal mode e_n is its squared overlap
with the normalized uniform displacement vector.  Completeness therefore makes
all residues sum to one.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class FKResult:
    positions: np.ndarray
    omega2: np.ndarray
    modes: np.ndarray
    optical_weights: np.ndarray
    real_ipr: np.ndarray
    k_ipr: np.ndarray


def _energy_grad(y: np.ndarray, *, beta: float, n_periods: int, atoms: int):
    a = 1.0
    lam = atoms * a / n_periods
    q = 2.0 * np.pi / lam
    v0 = 2.0 * beta / q**2
    x = np.arange(atoms, dtype=float) * a + y
    dy = np.roll(y, -1) - y
    energy = 0.5 * np.dot(dy, dy) + v0 * np.cos(q * x).sum()
    grad = 2.0 * y - np.roll(y, 1) - np.roll(y, -1) - v0 * q * np.sin(q * x)
    return float(energy), grad


def equilibrium(beta: float, n_periods: int, atoms: int) -> np.ndarray:
    """Find the lowest deterministic equilibrium among several phase seeds."""
    if beta < 0 or n_periods < 1 or atoms < 2:
        raise ValueError("invalid FK parameters")
    a = 1.0
    lam = atoms * a / n_periods
    seeds = []
    # Common shifts sample the substrate phase; a tiny deterministic ripple
    # helps L-BFGS leave a symmetry saddle if one is encountered.
    for phase in np.linspace(0.0, lam, 12, endpoint=False):
        seed = np.full(atoms, phase)
        seed += 1e-7 * np.sin(2.0 * np.pi * np.arange(atoms) / atoms)
        seeds.append(seed)
    best = None
    for seed in seeds:
        res = minimize(
            lambda yy: _energy_grad(yy, beta=beta, n_periods=n_periods, atoms=atoms),
            seed,
            jac=True,
            method="L-BFGS-B",
            options={"ftol": 1e-14, "gtol": 1e-11, "maxiter": 20000},
        )
        if not res.success and np.linalg.norm(res.jac, ord=np.inf) > 2e-7:
            continue
        if best is None or res.fun < best.fun:
            best = res
    if best is None:
        raise RuntimeError("FK equilibrium solve did not converge")
    return np.arange(atoms, dtype=float) * a + best.x


def solve(beta: float, n_periods: int, atoms: int) -> FKResult:
    x = equilibrium(beta, n_periods, atoms)
    lam = atoms / n_periods
    q = 2.0 * np.pi / lam
    # Hessian of 1/2 sum (x_{j+1}-x_j-a)^2 + V0 sum cos(Q x_j).
    h = np.zeros((atoms, atoms), float)
    np.fill_diagonal(h, 2.0 - 2.0 * beta * np.cos(q * x))
    idx = np.arange(atoms)
    h[idx, (idx + 1) % atoms] = -1.0
    h[idx, (idx - 1) % atoms] = -1.0
    omega2, modes = np.linalg.eigh((h + h.T) / 2.0)
    uniform = np.ones(atoms) / np.sqrt(atoms)
    weights = np.abs(modes.T @ uniform) ** 2
    weights /= weights.sum()
    real_ipr = np.sum(np.abs(modes) ** 4, axis=0)
    fourier = np.fft.fft(modes, axis=0) / np.sqrt(atoms)
    k_ipr = np.sum(np.abs(fourier) ** 4, axis=0).real
    return FKResult(x, omega2, modes, weights, real_ipr, k_ipr)


def local_optical_weights(result: FKResult) -> np.ndarray:
    """W[j,n] = Z_n |psi_n(j)|^2."""
    return np.abs(result.modes) ** 2 * result.optical_weights[None, :]
