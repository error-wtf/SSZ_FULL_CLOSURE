"""Contract tests for the REAL_SPECTROSCOPY_BRIDGE_V1 core modules."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from ssz_p5.observations.timing import (  # noqa: E402
    gap_comb_frequencies,
)
from ssz_p5.spectroscopy.mode_match import (  # noqa: E402
    match_one,
    null_comparison,
)
from ssz_p5.spectroscopy.physical_units import (  # noqa: E402
    FrozenMassPrior,
    omega_bar_to_hz,
    unit_conversion_factor_hz,
)
from ssz_p5.spectroscopy.theory_catalog import (  # noqa: E402
    assert_blind,
    build_predicted_catalog,
)

POLICY = {"freq_tol_rel": 0.005, "freq_tol_abs_hz": 0.2, "n_sigma": 3.0}


# ---------- physical units ----------

def test_mass_prior_is_hash_frozen():
    p = FrozenMassPrior("test", 10.0, 1.0, "unit test")
    j = p.to_json()
    assert len(j["sha256"]) == 64
    j2 = FrozenMassPrior("test", 10.0, 1.0, "unit test").to_json()
    assert j["sha256"] == j2["sha256"]
    j3 = FrozenMassPrior("test", 11.0, 1.0, "unit test").to_json()
    assert j["sha256"] != j3["sha256"]


def test_omega_conversion_scales_inverse_with_mass():
    f10 = omega_bar_to_hz(1.0, 10.0)
    f20 = omega_bar_to_hz(1.0, 20.0)
    assert f20 == pytest.approx(f10 / 2.0)
    assert unit_conversion_factor_hz(10.0) == pytest.approx(f10)


# ---------- blind theory catalog ----------

def test_blindness_assert_fires_on_observed_files(tmp_path):
    (tmp_path / "data/observed").mkdir(parents=True)
    (tmp_path / "data/observed/OBSERVED_MODE_CATALOG.json").write_text("{}")
    with pytest.raises(RuntimeError, match="BLINDNESS VIOLATION"):
        assert_blind(tmp_path)


def test_predicted_catalog_uses_frozen_prior(tmp_path):
    prior = FrozenMassPrior("unit", 10.0, 1.0, "t")
    rows = [{"omega_bar_re": 0.1, "Z_q": 0.5, "ipr": 0.01,
             "localization_class": "EXTENDED"}]
    cat = build_predicted_catalog(rows, prior)
    assert cat["n_modes"] == 1
    assert cat["modes"][0]["f_hz"] == pytest.approx(
        0.1 * cat["unit_conversion_hz_per_omega_bar"])


# ---------- matching ----------

def test_match_within_and_outside_windows():
    pred = [{"mode_index": 0, "f_hz": 100.05},
            {"mode_index": 1, "f_hz": 130.0}]
    obs = {"f_hz": 100.0, "sigma_f_hz": 0.05}
    m = match_one(obs, pred, POLICY)
    assert m["matched"] and m["predicted_index"] == 0
    obs_far = {"f_hz": 100.0, "sigma_f_hz": 0.05}
    pred_far = [{"mode_index": 0, "f_hz": 101.0}]  # 1% > 0.5% tol
    assert not match_one(obs_far, pred_far, POLICY)["matched"]


def test_null_pvalue_small_for_real_matches():
    obs = [{"f_hz": f, "sigma_f_hz": 0.05} for f in (100.0, 200.0, 300.0)]
    pred = [{"mode_index": i, "f_hz": f} for i, f in enumerate(
        [100.01, 200.02, 299.99, 130.0, 260.0, 390.0])]
    n = null_comparison(obs, pred, POLICY, n_draws=500)
    assert n["real_hits"] == 3
    assert n["p_null"] < 0.05


def test_null_pvalue_flat_for_random():
    rng = np.random.default_rng(3)
    obs = [{"f_hz": float(f), "sigma_f_hz": 0.05} for f in rng.uniform(100, 110, 3)]
    pred = [{"mode_index": i, "f_hz": float(f)} for i, f in enumerate(rng.uniform(200, 300, 8))]
    n = null_comparison(obs, pred, POLICY, n_draws=500)
    assert n["real_hits"] == 0
    assert n["p_null"] > 0.5


# ---------- gap artifact (already covered deeper in its own file) ----------

def test_comb_freq_definition():
    assert gap_comb_frequencies(0.02) == pytest.approx(50.0)
