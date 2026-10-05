"""End-to-end: schema-conformant frozen catalogs -> compare tool verdicts.

Covers the full Stage-A→B→C contract chain with the REAL builders and
the REAL compare tool as subprocess (not a re-implementation).
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ssz_p5.observations.blind_catalog import (  # noqa: E402
    build_observed_catalog,
    freeze_observed_catalog,
)
from ssz_p5.spectroscopy.physical_units import FrozenMassPrior  # noqa: E402
from ssz_p5.spectroscopy.theory_catalog import (  # noqa: E402
    build_predicted_catalog,
    freeze_predicted_catalog,
)

CONFIG = {
    "psd_normalization": "leahy",
    "freq_range_hz": [0.1, 2048.0],
    "energy_bands_kev": [0.5, 2.0, 10.0],
    "significance_threshold_sigma": 4.0,
    "trial_correction": "bonferroni",
    "frozen_before_first_psd": True,
}


def _write_frozen_pair(tmp_path: Path, obs_freqs, pred_freqs):
    obs = build_observed_catalog(
        target="T", obsids=["x"], total_exposure_s=1.0, analysis_config=CONFIG,
        candidates=[{"f_hz": f, "sigma_f_hz": 0.05, "amplitude_rms": 0.01,
                     "significance_sigma": 6.0, "obs_occurrence": 1}
                    for f in obs_freqs],
        instrumental_controls=[])
    frozen_obs, _ = freeze_observed_catalog(obs)
    (tmp_path / "observed.json").write_text(json.dumps(frozen_obs))

    prior = FrozenMassPrior("p", 10.0, 1.0, "e2e")
    rows = [{"omega_bar_re": 1.0 + 0.7 * i, "observable": "q1", "u": 0.65,
             "L": 6, "omega2": 2.0 + i, "cluster_residue_weight": 0.9,
             "ipr": 0.02, "localization_class": "EXTENDED"}
            for i in range(len(pred_freqs))]
    pred = build_predicted_catalog(
        rows, prior, member_hash="0" * 64,
        healthy_window={"u_min": 0.62, "u_max": 0.70, "L_values": [6],
                        "controls_all_pass": True})
    conv = pred["unit_conversion_hz_per_omega_bar"]
    for i, f in enumerate(pred_freqs):
        pred["modes"][i]["omega_bar_re"] = f / conv
        pred["modes"][i]["f_hz"] = f
    frozen_pred, _ = freeze_predicted_catalog(pred)
    (tmp_path / "predicted.json").write_text(json.dumps(frozen_pred))

    for n in ("observed", "predicted"):
        h = hashlib.sha256((tmp_path / f"{n}.json").read_bytes()).hexdigest()
        (tmp_path / f"{n}.json.sha256").write_text(h + "\n")


def _run_compare(tmp_path: Path, out="report.json"):
    return subprocess.run(
        [sys.executable, "tools/compare_frozen_mode_catalogs.py",
         "--observed", str(tmp_path / "observed.json"),
         "--predicted", str(tmp_path / "predicted.json"),
         "--output", str(tmp_path / out)],
        capture_output=True, text=True, cwd=ROOT)


def test_e2e_true_matches_above_null(tmp_path):
    """3 real matches among 8 predicted -> verdict MATCH_ABOVE_NULL."""
    _write_frozen_pair(tmp_path,
                       obs_freqs=[100.0, 200.0, 300.0],
                       pred_freqs=[100.01, 200.02, 299.99,
                                   130.0, 260.0, 390.0, 520.0, 650.0])
    r = _run_compare(tmp_path)
    assert r.returncode == 0, r.stderr
    rep = json.loads((tmp_path / "report.json").read_text())
    assert rep["verdict"] == "MATCH_ABOVE_NULL"
    assert rep["null_comparison"]["real_hits"] == 3
    assert rep["null_comparison"]["p_null"] < 0.05


def test_e2e_random_gives_no_match(tmp_path):
    """Uncorrelated catalogs must NOT produce a match verdict."""
    _write_frozen_pair(tmp_path,
                       obs_freqs=[100.0, 200.0, 300.0],
                       pred_freqs=[530.0, 560.0, 590.0, 620.0, 650.0, 680.0,
                                   710.0, 740.0])
    r = _run_compare(tmp_path)
    assert r.returncode == 0, r.stderr
    rep = json.loads((tmp_path / "report.json").read_text())
    assert rep["verdict"] == "NO_MATCH_ABOVE_NULL"


def test_e2e_freeze_tampering_rejected(tmp_path):
    """Editing a frozen catalog after signing must abort the comparison."""
    _write_frozen_pair(tmp_path,
                       obs_freqs=[100.0], pred_freqs=[100.01])
    doc = json.loads((tmp_path / "observed.json").read_text())
    doc["modes"][0]["f_hz"] = 42.0  # tamper
    (tmp_path / "observed.json").write_text(json.dumps(doc))
    r = _run_compare(tmp_path)
    assert r.returncode != 0
    assert "FREEZE VIOLATION" in r.stderr
