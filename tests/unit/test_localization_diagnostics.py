"""Li-style localization diagnostics: analytically known reference states."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.qnm.localization_diagnostics import (  # noqa: E402
    d2_from_block_scaling,
    diagnose_mode,
    ipr_from_p,
)


def _identity_K(n):
    K = np.zeros((n, 1, 1))
    K[:, 0, 0] = 1.0
    return K


def test_delta_peak_is_localized():
    n = 1024
    psi = np.zeros((n, 1))
    psi[n // 2, 0] = 1.0
    r = diagnose_mode(0, psi, _identity_K(n), np.ones(n))
    assert r.ipr == pytest.approx(1.0, abs=1e-12)
    assert r.d2_scaling == pytest.approx(0.0, abs=0.05)
    assert r.classification == "LOCALIZED"


def test_sinus_mode_is_extended():
    n = 1024
    x = np.linspace(0, 1, n)
    psi = np.zeros((n, 1))
    psi[:, 0] = np.sqrt(2.0) * np.sin(np.pi * x)
    r = diagnose_mode(0, psi, _identity_K(n), np.ones(n))
    assert r.ipr == pytest.approx(0.001953, rel=0.6)  # ~ 3/(2n) for sin^2 weighting
    assert r.d2_scaling == pytest.approx(1.0, abs=0.05)
    assert r.classification == "EXTENDED"


def test_narrow_gaussian_intermediate():
    n = 1024
    x = np.linspace(0, 1, n)
    psi = np.zeros((n, 1))
    psi[:, 0] = np.exp(-0.5 * ((x - 0.5) / 0.01) ** 2)
    psi /= np.linalg.norm(psi)
    r = diagnose_mode(0, psi, _identity_K(n), np.ones(n))
    assert r.classification in ("LOCALIZED", "MULTIFRACTAL")
    assert r.ipr < 0.2


def test_ipr_bounds():
    rng = np.random.default_rng(42)
    p = rng.random(500)
    p /= p.sum()
    ipr = ipr_from_p(p)
    assert 1.0 / 500 <= ipr <= 1.0


def test_d2_monotone_in_localization_length():
    """Narrower state -> smaller D2 (sanity of the scaling estimator)."""
    n = 1024
    x = np.linspace(0, 1, n)
    ds = []
    for sigma in (0.005, 0.02, 0.08):
        psi = np.exp(-0.5 * ((x - 0.5) / sigma) ** 2)
        psi /= psi.sum()
        ds.append(d2_from_block_scaling(psi))
    assert ds[0] < ds[1] < ds[2]
