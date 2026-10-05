#!/usr/bin/env python3
"""V4 ZERTifizierungs-Vervollständigung: alle Branches des eps-Family-
Scans, die den Core-Rand erreichen, mit dem Interior-Fenster zertifizieren
(auf 4 signifikante a1/eps-Kombinationen pro Vorzeichenquadrant verdichtet).

Das V4-Ergebnis (H3) steht auf EINEM zertifizierten Zweig; diese Kampagne
macht daraus einen Zweig-VERBAND: dieselbe sign(K)=-sign(a1*eps)-Struktur
über alle Quadranten, jedes Mitglied einzeln residuen-zertifiziert.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

import importlib.util as _ilu  # noqa: E402

_spec = _ilu.spec_from_file_location(
    "v4s", ROOT / "tools" / "run_luminal_background_solve_v4.py")
_v4 = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_v4)

ARCHIVE = ROOT / (
    "data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv"
)
OUT_DIR = ROOT / "data/generated/spectral"
MARGIN = 8


def main() -> int:
    t0 = time.time()
    archive = pd.read_csv(ARCHIVE).sort_values("u").reset_index(drop=True)
    r_all = (1.0 / archive.u.to_numpy(float))[::-1]
    r_min, r_core = float(r_all[0]), float(r_all[-1])
    n_out = len(r_all)
    interior = slice(MARGIN, n_out - MARGIN)

    red, facts, _e = _v4.build_reduction()
    assert all(facts.values())

    # 8 Zweige: 4 a1-Werte × 2 eps-Vorzeichen, mittlere Amplitude
    grid = []
    for a1v in (0.5, -0.5, 2.0, -2.0):
        for eps in (0.3, -0.3):
            grid.append((a1v, eps))

    rows = []
    for a1v, eps in grid:
        try:
            sol = _v4.integrate_branch(
                _v4.rhs_factory(red, a1v), r_min, r_core, eps)
        except Exception as exc:  # noqa: BLE001
            rows.append({"a1": a1v, "eps": eps, "reached": False,
                         "reason": type(exc).__name__})
            continue
        reached = bool(sol.t[-1] >= r_core * (1 - 1e-9))
        if not reached:
            rows.append({"a1": a1v, "eps": eps, "reached": False,
                         "reason": f"blowup r={float(sol.t[-1]):.3f}"})
            continue
        y = sol.sol(r_all)
        res = _v4.residuals_on_grid(r_all, y[0], y[1], y[2], a1v)
        e00 = float(np.max(np.abs(res["E00"][interior])))
        e11 = float(np.max(np.abs(res["E11"][interior])))
        e22 = float(np.max(np.abs(res["E22"][interior])))
        # frozen thresholds (V4 campaign): E00/E11 2e-6, E22 2e-4 —
        # shooting + jet9d8 discretization floor measured at 1.9e-6
        conv = e00 <= 2e-6 and e11 <= 2e-6 and e22 <= 2e-4

        # K on-shell (V3-certified identity)
        from ssz_p5.jets.jet9d8 import derivative as jet
        fp_g = jet(r_all, y[0])
        hp_g = jet(r_all, y[1])
        F_of_r = 1.0 + 2.0 * a1v * (y[2] - 1.0)
        K = F_of_r * r_all / 2.0 * (fp_g / y[0] - hp_g / y[1])
        fp13 = jet(r_all, y[0], window=13)
        hp13 = jet(r_all, y[1], window=13)
        K13 = F_of_r * r_all / 2.0 * (fp13 / y[0] - hp13 / y[1])
        k_err = float(np.max(np.abs(K - K13)))
        band = (r_all >= 0.285) & (r_all <= 0.325)
        k_err_eff = max(k_err, e00)
        sign_all = ("K<0_all" if np.max(K) < 0
                    else "K>0_all" if np.min(K) > 0 else "sign-change")
        decidable = not (abs(float(np.min(K))) < k_err_eff
                         and abs(float(np.min(K[band]))) < k_err_eff)
        if not decidable:
            verdict = "NOT_EVALUABLE"
        elif float(np.min(K[band])) < 0:
            verdict = "H3_PHYSICAL_GHOST"
        elif float(np.min(K)) < 0:
            verdict = "NEW_PATHOLOGY_MAP"
        else:
            verdict = "H1_CANDIDATE_PENDING_FINITE_L_HEALTH"

        rows.append({
            "a1": a1v, "eps": eps, "reached": True,
            "converged": bool(conv),
            "resid_interior": {"E00": e00, "E11": e11, "E22": e22},
            "K_min": float(np.min(K)), "K_max": float(np.max(K)),
            "K_band_min": float(np.min(K[band])),
            "K_stencil_diff": k_err,
            "sign_structure": sign_all,
            "verdict": verdict,
        })
        print(f"a1={a1v:+.1f} eps={eps:+.1f}: conv={conv} "
              f"K in [{np.min(K):.3e},{np.max(K):.3e}] -> {verdict}",
              flush=True)

    df = pd.DataFrame(rows)
    out_csv = OUT_DIR / "N4_V4_BRANCH_VERBAND.csv"
    df.to_csv(out_csv, index=False)
    converged = df[df.get("converged") == True]  # noqa: E712
    verdicts = converged["verdict"].value_counts().to_dict() if len(converged) else {}
    result = {
        "audit": "N4_V4_BRANCH_VERBAND",
        "purpose": "single-branch H3 verdict extended to an 8-branch band "
                   "across all a1/eps quadrants, each member individually "
                   "residual-certified",
        "branches": rows,
        "n_converged": int(len(converged)),
        "verdict_counts": verdicts,
        "sign_law_consistent": bool(
            len(converged) > 0 and all(
                (r["verdict"] == "H3_PHYSICAL_GHOST") == (r["a1"] * r["eps"] > 0)
                for _, r in converged.iterrows()
            )
        ),
        "decision_rules_unchanged": True,
        "wall_seconds": round(time.time() - t0, 1),
    }
    (OUT_DIR / "N4_V4_BRANCH_VERBAND.json").write_text(
        json.dumps(result, indent=1, allow_nan=False) + "\n")
    print(json.dumps({
        "n_converged": result["n_converged"],
        "verdicts": verdicts,
        "sign_law_consistent": result["sign_law_consistent"],
        "written": str(out_csv.relative_to(ROOT)),
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
