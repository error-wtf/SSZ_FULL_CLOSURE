#!/usr/bin/env python3
"""Build the common Weisz-Li-SSZ spectral-comparison format.

Produces ONE schema that all three model families feed, so cross-model
falsification tests (execution plan stage D) compare dimensionless invariants
rather than pictures:

  SpectralEntry = {
    model: "weisz1978_fk" | "li2023_cst_aah" | "ssz_krgsm_local",
    case: str,
    n_modes, top_residue_fractions (1..3), effective_mode_count,
    spectral_entropy, residue_dynamic_range,
    frequency_separation_vs_residue_separation (Spearman),
    spatial diagnostics where defined: IPR median per region, region split
  }

Inputs (all machine-generated, no hand-entered numbers):
  data/generated/spectral/weisz1978/WEISZ1978_FK_REPRODUCTION.json + CSV
  data/generated/spectral/li2023_cst_aah/LI2023_CSTAAH_REPRODUCTION.json
  data/generated/spectral/SSZ_SPECTROSCOPY_FIXED_OBSERVABLE_V3.json (if present)

Anti-circularity: this script only REFORMATS existing audited outputs. It does
not fit, scale or retune anything.  Missing inputs are recorded as "absent",
never fabricated.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/generated/spectral/COMMON_SPECTRAL_FORMAT_V1.json"


def top_fractions(weights, k=3):
    w = np.sort(np.asarray(weights, float))[::-1]
    s = w.sum()
    return [float(w[:n].sum() / s) for n in range(1, k + 1)]


def entropy(weights):
    w = np.asarray(weights, float)
    w = w / w.sum()
    return float(-(w[w > 0] * np.log(w[w > 0])).sum())


def n_eff(weights):
    w = np.asarray(weights, float)
    w = w / w.sum()
    return float(1.0 / (w**2).sum())


def dynamic_range(weights):
    w = np.sort(np.asarray(weights, float))[::-1]
    w = w / w.sum()
    pos = w[w > 0]
    return float(pos[0] / pos[-1]) if len(pos) > 1 else None


def gap_residue_spearman(omega2, weights):
    """Near-degeneracy vs residue-suppression signature:
    correlate |dlog omega2| between consecutive ranked modes with the residue
    ratio of the same pair."""
    o = np.asarray(omega2, float); w = np.asarray(weights, float)
    m = (o > 0) & (w > 0)
    o, w = o[m], w[m]
    idx = np.argsort(o); o, w = o[idx], w[idx]
    if len(o) < 8:
        return None
    do = np.abs(np.diff(np.log(o)))
    rw = np.log(w[:-1] / w[1:])
    from scipy.stats import spearmanr
    rho, _ = spearmanr(do, rw)
    return float(rho)


def collect_weisz() -> list[dict]:
    j = json.load(open(ROOT / "data/generated/spectral/weisz1978/WEISZ1978_FK_REPRODUCTION.json"))
    entries = []
    for c in j["cases"]:
        entry = {
            "model": "weisz1978_fk",
            "case": f"{c['N_lambda_eq_M_a']} beta={c['beta']}",
            "n_modes": None,
            "top_residue_fractions": top_fractions(c["top_weights"]),
            "effective_mode_count": c["effective_mode_count"],
            "spectral_entropy": c["spectral_entropy"],
            "residue_dynamic_range": dynamic_range(c["top_weights"]) if len(c["top_weights"]) > 1 else None,
            "lowest_omega2": c["lowest_omega2"],
            "hessian_min": c["hessian_min"],
            "spatial": None,
        }
        entries.append(entry)
    # full mode tables from CSV for the gap/residue signature
    df = pd.read_csv(ROOT / "data/generated/spectral/weisz1978/WEISZ1978_FK_MODES.csv")
    for (nl, m, b), g in df.groupby(["N", "M", "beta"]):
        lab = f"{nl} lambda = {m} a beta={b}"
        rho = gap_residue_spearman(g.omega2.to_numpy(), g.q0_residue.to_numpy())
        for e in entries:
            if e["case"] == lab:
                e["gap_residue_spearman"] = rho
                e["n_modes"] = int(len(g))
    return entries


def collect_li() -> list[dict]:
    p = ROOT / "data/generated/spectral/li2023_cst_aah/LI2023_CSTAAH_REPRODUCTION.json"
    if not p.exists():
        return [{"model": "li2023_cst_aah", "case": "ABSENT_RUN_IN_PROGRESS"}]
    j = json.load(open(p))
    st = j["stages"]
    entries = []
    if "fig2_theta_averaged" in st and "4181" in st["fig2_theta_averaged"]:
        d = st["fig2_theta_averaged"]["4181"]
        b_loc = d["btilde_localized_median"]; b_ext = d["btilde_extended_median"]
        entries.append({
            "model": "li2023_cst_aah",
            "case": f"N=4181 sigma=1 lambda=1.5 100theta",
            "jc": d["jc"],
            "frac_states_localized_region": d["frac_states_localized_region"],
            "spatial": {
                "diagnostic": "fractal_dimension_btilde",
                "localized_median": b_loc,
                "extended_median": b_ext,
                "separation": float(b_ext - b_loc),
            },
            "gates": j.get("gates"),
        })
    if "flat_aah_control" in st:
        for tag, d in st["flat_aah_control"].items():
            entries.append({
                "model": "li2023_cst_aah",
                "case": f"flat_control_{tag} lambda={d['lambda']}",
                "spatial": {"diagnostic": "fractal_dimension_btilde",
                            "median": d["btilde_median"],
                            "ipr_median": d["ipr_median"]},
            })
    if "wavepacket_dynamics" in st:
        for name, d in st["wavepacket_dynamics"].items():
            t0, t1 = d["trajectory"][0], d["trajectory"][-1]
            entries.append({
                "model": "li2023_cst_aah",
                "case": f"wavepacket_{name}",
                "participation_growth_ratio": float(t1["participation"] / max(t0["participation"], 1e-300)),
            })
    return entries


def collect_ssz() -> list[dict]:
    p = ROOT / "data/generated/spectral/SSZ_SPECTROSCOPY_FIXED_OBSERVABLE_V3.json"
    if not p.exists():
        return [{"model": "ssz_krgsm_local", "case": "ABSENT_NOT_YET_RERUN_ON_BRANCH"}]
    j = json.load(open(p))
    lp = j["local_principal_spectroscopy_v3"]
    entries = [{
        "model": "ssz_krgsm_local",
        "case": "fixed_observable_v3_status",
        "status": lp["status"],
        "scope": lp["scope"],
        "all_controls_pass": lp["all_controls_pass"],
        "interpretation_guard": j.get("interpretation_guard"),
    }]
    for L, info in lp["per_L"].items():
        for obs_name, obs in info["observables"].items():
            e = {
                "model": "ssz_krgsm_local",
                "case": f"healthy_window L={L} {obs_name}",
                "u_window": [info["u_min"], info["u_max"]],
                "samples": info["samples"],
                "nondegenerate_fraction": info["nondegenerate_sample_fraction"],
                "cluster_top1_range": obs.get("cluster_top1_range"),
                "cluster_top2_range": obs.get("cluster_top2_range"),
                "cluster_neff_range": obs.get("cluster_neff_range"),
                "robust_inversions": obs.get("robust_nondegenerate_inversion_count"),
                # dimensionless class summary: top2 lower bound
                "top2_lower_bound": (min(obs["cluster_top2_range"])
                                     if obs.get("cluster_top2_range") else None),
            }
            entries.append(e)
    return entries


def main():
    entries = collect_weisz() + collect_li() + collect_ssz()
    rep = {
        "format_version": 1,
        "description": "Common dimensionless spectral-comparison schema (Weisz/Li/SSZ)",
        "invariants": ["top_residue_fractions", "effective_mode_count",
                        "spectral_entropy", "residue_dynamic_range",
                        "gap_residue_spearman", "fractal_dimension",
                        "participation_growth_ratio"],
        "anti_circularity": "reformat-only; no fitting; absent != zero",
        "entries": entries,
    }
    OUT.write_text(json.dumps(rep, indent=1))
    n_models = sorted({e["model"] for e in entries})
    print("models present:", n_models)
    for e in entries:
        line = f"  {e['model']:<18} {e.get('case','')[:45]:<45}"
        if e.get("top_residue_fractions"):
            line += f" top3={['%.3f'%x for x in e['top_residue_fractions']]}"
        if e.get("gap_residue_spearman") is not None:
            line += f" rho(gap,res)={e['gap_residue_spearman']:+.2f}"
        print(line)
    print("written:", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
