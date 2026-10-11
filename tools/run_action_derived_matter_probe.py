#!/usr/bin/env python3
"""G-E Robustheits-Pruefstand (Lino-Direktive 2026-10-11, Drei-Arbeiten-Runde):

Testet, ob NEW_STRUCTURE_NOT_COLLAPSIBLE (eingefroren Run 11: alpha_best
0.9191907077802026, max resid 0.5642466802128533) unter den drei offenen
Modellierungs-Alternativen aus full-report.md § 6 bestehen bleibt:

  P  Proper-radiale Substratphase:  phi_n = Q * l(r_n) + Q * s(r_n) * u_n
     mit l(r) = (r - R) + (r_s/2) ln(r/R)  =>  dl/dr = s(r)  (analytisch).
     Ersetzt die Mehrskalen-Identifikation phi_n = Q s_n (a n + u_n).
  K  Aktions-Kinetik (Kandidat § 6.1):  Trägheit m s^2/D  statt  m/D^2.
     (Für den g1-Zweig ist s^2/D = (1+xi)^3 = 1/D^3.)
  PK beides kombiniert.

Flache Referenz (ssz=False) ist identisch unter allen Varianten
(D=s=1 => beide Phasenidentifikationen und beide Trägheiten fallen auf m zusammen).

Unabhängige Modenzuordnung (§ 6.2): statt rank-by-rank wird per
Hungarian-Assignment auf der Überlappungsmatrix |V_flat^T V_ssz| der
massgewichteten Symmetrie-Moden gematcht; G-E-alpha und Rest kommen dann
aus den GEMATCHTEN Paaren, nicht aus der Rangordnung.

Gates VOR dem Lauf deklariert. Gate-Ergebnis P3 ist die eigentliche Frage:
G-E robust => Nichtkollapsibilität ist keine Eigenschaft der Modellwahl.
KEIN astrophysikalischer Detektionsclaim.
"""
from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
OUT = ROOT / "artifacts" / "SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json"

_spec = importlib.util.spec_from_file_location(
    "fk_ref", ROOT / "tools" / "run_ssz_fk_coupled_model.py")
fk_ref = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fk_ref)

FKModel = fk_ref.FKModel

# Eingefrorene Run-11-Referenz (artifacts/SSZ_FK_COUPLED_MODEL_V1.json)
FROZEN = {"alpha_best": 0.9191907077802026,
          "max_rel_resid_after_alpha_fit": 0.5642466802128533}

# Testkonfiguration identisch zum eingefrorenen G-E (M=60, R=12, dr=0.05)
GE_KW = dict(M=60, R=12.0, dr=0.05, r_s=4.4, V0=1.0)


class ProbeModel(FKModel):
    """FKModel mit umschaltbarer Phasen-Identifikation und Trägheitsform.

    phase_mode:
      multiscale (eingefroren): phi_n = Q s_n (a n + u_n)
      proper (§ 6.1a):          phi_n = Q [l(r_n) + s_n u_n]
    inertia_mode:
      d2 (eingefroren):         m_n = m / D_n^2
      s2_over_d (§ 6.1):        m_n = m s_n^2 / D_n
    """

    def __init__(self, *args, phase_mode: str = "multiscale",
                 inertia_mode: str = "d2", **kwargs):
        # Geometrie-Slots vor super().__init__ aus den Argumenten ziehen
        # (FKModel speichert R/dr/r_s nicht als Attribute).
        all_kw = dict(zip(("M", "R", "dr", "r_s"), args, strict=False), **kwargs)
        self.R_val = float(all_kw.get("R", 12.0))
        self.dr_val = float(all_kw.get("dr", 0.05))
        self.r_s_val = float(all_kw.get("r_s", 4.4))
        super().__init__(*args, **kwargs)
        self.phase_mode = phase_mode
        self.inertia_mode = inertia_mode
        r = self.R_val + self.dr_val * np.arange(self.M, dtype=float)
        self.l_ell = ((r - self.R_val)
                      + 0.5 * self.r_s_val * np.log(r / self.R_val))
        if inertia_mode == "s2_over_d":
            self.mass = self.m * self.s**2 / self.d
        elif inertia_mode != "d2":
            raise ValueError(inertia_mode)

    def _phase(self, u: np.ndarray) -> np.ndarray:
        n = self.a * np.arange(self.M)
        if self.phase_mode == "multiscale":
            return self.Q * self.s * (n + u)
        if self.phase_mode == "proper":
            return self.Q * (self.l_ell + self.s * u)
        raise ValueError(self.phase_mode)

    def u_pot(self, u: np.ndarray) -> float:
        du = np.diff(u)
        strain = self.s_mid * (self.a + du) - self.a
        u_el = 0.5 * self.kap * float(np.sum(strain * strain))
        u_sub = -self.V0 * float(np.sum(np.cos(self._phase(u))))
        return u_el + u_sub

    def grad(self, u: np.ndarray) -> np.ndarray:
        du = np.diff(u)
        strain = self.s_mid * (self.a + du) - self.a
        g = np.zeros(self.M)
        fel = self.kap * strain * self.s_mid
        g[:-1] -= fel
        g[1:] += fel
        g += self.V0 * self.Q * self.s * np.sin(self._phase(u))
        return g

    def hessian(self, u: np.ndarray) -> np.ndarray:
        H = np.zeros((self.M, self.M))
        k = self.kap * self.s_mid**2
        idx = np.arange(self.M - 1)
        H[idx, idx] += k
        H[idx + 1, idx + 1] += k
        H[idx, idx + 1] -= k
        H[idx + 1, idx] -= k
        hsub = self.V0 * (self.Q * self.s) ** 2 * np.cos(self._phase(u))
        H[np.arange(self.M), np.arange(self.M)] += hsub
        return H


def fd_check(model: ProbeModel, seed: int = 0) -> dict:
    """Zentrale Differenzen gegen grad und hessian (Korrektheit der Ableitungen)."""
    rng = np.random.default_rng(seed)
    u = 0.05 * rng.standard_normal(model.M)
    h = 1e-6
    g = model.grad(u)
    gfd = np.zeros(model.M)
    for i in range(model.M):
        up, um = u.copy(), u.copy()
        up[i] += h
        um[i] -= h
        gfd[i] = (model.u_pot(up) - model.u_pot(um)) / (2.0 * h)
    scale_g = max(float(np.max(np.abs(g))), 1e-12)
    grad_err = float(np.max(np.abs(g - gfd)) / scale_g)

    H = model.hessian(u)
    Hfd = np.zeros_like(H)
    for i in range(model.M):
        up, um = u.copy(), u.copy()
        up[i] += h
        um[i] -= h
        Hfd[:, i] = (model.grad(up) - model.grad(um)) / (2.0 * h)
    scale_H = max(float(np.max(np.abs(H))), 1e-12)
    hess_err = float(np.max(np.abs(H - Hfd)) / scale_H)
    return {"grad_rel_err": grad_err, "hessian_rel_err": hess_err}


def weighted_modes(model: ProbeModel, u_eq: np.ndarray):
    """(w, V) mit V = Orthonormalbasis des massgewichteten Problems."""
    w2, u_modes, _ = model.modes(u_eq)
    sqrt_mass = np.sqrt(model.mass)
    V = sqrt_mass[:, None] * u_modes          # inverse Abbildung von modes()
    return np.sqrt(np.maximum(w2, 0.0)), V


def ge_analysis(model_ssz: ProbeModel, u_s: np.ndarray,
                wf: np.ndarray, Vf: np.ndarray) -> dict:
    """G-E alpha/Rest: dreistufig.
    (a) rank-by-rank (eingefrorene Prozedur),
    (b) Hungarian-Matching auf Einzelvektoren (ill-posed bei Entartung,
        wird nur als Diagnose gefuehrt),
    (c) entartungs-bewusster Unterraum-Matching: Soft-Moden werden zu
        Frequenz-Gruppen geclustert (rel. Toleranz 1e-3); Gruppen werden
        per Hauptwinkel (Projektor-Overlap) zugeordnet; alpha/Rest aus den
        gruppengemittelten Frequenzen.
    """
    ws, Vs = weighted_modes(model_ssz, u_s)
    soft_f = wf > 1e-2 * wf.max()
    soft_s = ws > 1e-2 * ws.max()

    # --- rank-by-rank (eingefroren) ---
    spec_f = np.sort(wf[soft_f])[:20]
    spec_s = np.sort(ws[soft_s])[:20]
    a_rank = float(np.sum(spec_s * spec_f) / np.sum(spec_f**2))
    r_rank = float(np.max(np.abs(spec_s - a_rank * spec_f) / (a_rank * spec_f)))

    # --- Overlap-Matching (Hungarian, maximale Gesamtüberlappung) ---
    Ov = np.abs(Vf.T @ Vs)
    ri, cj = linear_sum_assignment(-Ov)
    # gepaarte soft/soft-Moden, sortiert nach flacher Frequenz, Top 20
    pairs = [(wf[i], ws[j], float(Ov[i, j]))
             for i, j in zip(ri, cj, strict=False)
             if soft_f[i] and soft_s[j]]
    pairs.sort(key=lambda p: p[0])
    top = pairs[:20]
    mf = np.array([p[0] for p in top])
    ms = np.array([p[1] for p in top])
    a_ov = float(np.sum(ms * mf) / np.sum(mf**2))
    r_ov = float(np.max(np.abs(ms - a_ov * mf) / (a_ov * mf)))
    # Zuordnungs-Flip vs. reiner Rangordnung: wie viele gematchte Paare
    # haben NICHT denselben Rang in ihrem jeweiligen soft-Spektrum?
    rank_f = {f: k for k, f in enumerate(np.sort(wf[soft_f]))}
    rank_s = {s: k for k, s in enumerate(np.sort(ws[soft_s]))}
    flips = int(sum(1 for f, s, _ in top if rank_f[f] != rank_s[s]))

    # --- (c) entartungs-bewusstes Gruppen-Matching (Unterraum-Vergleich) ---
    def freq_groups(w: np.ndarray, tol: float = 1e-3):
        """Indizes sortierter Soft-Frequenzen zu Gruppen clustern."""
        order = np.argsort(w)
        groups, cur = [], [order[0]]
        for k in order[1:]:
            if abs(w[k] - w[cur[-1]]) <= tol * w[cur[-1]]:
                cur.append(k)
            else:
                groups.append(cur)
                cur = [k]
        groups.append(cur)
        return groups

    gf = freq_groups(wf[soft_f])
    gs = freq_groups(ws[soft_s])
    # Projektoren der Gruppen im vollen Modenraum
    idx_f = np.flatnonzero(soft_f)
    idx_s = np.flatnonzero(soft_s)
    Pf = [Vf[:, idx_f[g]] @ Vf[:, idx_f[g]].T for g in gf]
    Ps = [Vs[:, idx_s[g]] @ Vs[:, idx_s[g]].T for g in gs]
    # Hauptwinkel pro Paar: singulaerwerte von P_f P_s
    principal = np.array([[np.linalg.svd(A @ B, compute_uv=False)[0]
                           for B in Ps] for A in Pf])
    gri, gjj = linear_sum_assignment(-principal)
    g_pairs = [(float(np.mean(wf[idx_f[gf[i]]])),
                float(np.mean(ws[idx_s[gs[j]]])),
                float(principal[i, j]),
                len(gf[i]), len(gs[j]))
               for i, j in zip(gri, gjj, strict=False)]
    g_pairs.sort(key=lambda p: p[0])
    top_g = [p for p in g_pairs if p[3] == p[4]][:20]  # nur gleiche Multiplizität
    gf_freq = np.array([p[0] for p in top_g])
    gs_freq = np.array([p[1] for p in top_g])
    a_grp = float(np.sum(gs_freq * gf_freq) / np.sum(gf_freq**2))
    r_grp = float(np.max(np.abs(gs_freq - a_grp * gf_freq) / (a_grp * gf_freq)))
    return {
        "rank_by_rank": {"alpha_best": a_rank,
                         "max_rel_resid_after_alpha_fit": r_rank,
                         "verdict": ("NEW_STRUCTURE_NOT_COLLAPSIBLE"
                                     if r_rank >= 1e-5 else
                                     "COLLAPSIBLE_TO_UNIFORM_FACTOR"),
                         "n_pairs": int(len(spec_f))},
        "overlap_matching": {"alpha_best": a_ov,
                             "max_rel_resid_after_alpha_fit": r_ov,
                             "verdict": ("NEW_STRUCTURE_NOT_COLLAPSIBLE"
                                         if r_ov >= 1e-5 else
                                         "COLLAPSIBLE_TO_UNIFORM_FACTOR"),
                             "n_pairs": len(top),
                             "mean_overlap": float(np.mean([p[2] for p in top])),
                             "min_overlap": float(np.min([p[2] for p in top])),
                             "assignment_flips_vs_rank": flips},
        "subspace_group_matching": {
            "alpha_best": a_grp,
            "max_rel_resid_after_alpha_fit": r_grp,
            "verdict": ("NEW_STRUCTURE_NOT_COLLAPSIBLE"
                        if r_grp >= 1e-5 else
                        "COLLAPSIBLE_TO_UNIFORM_FACTOR"),
            "n_groups": len(top_g),
            "min_principal_angle_cos": float(np.min([p[2] for p in top_g])),
            "mean_principal_angle_cos": float(np.mean([p[2] for p in top_g]))},
    }


VARIANTS = [
    ("baseline", dict(phase_mode="multiscale", inertia_mode="d2")),
    ("proper_phase", dict(phase_mode="proper", inertia_mode="d2")),
    ("action_inertia", dict(phase_mode="multiscale", inertia_mode="s2_over_d")),
    ("proper_plus_inertia", dict(phase_mode="proper", inertia_mode="s2_over_d")),
]


def main() -> int:
    t0 = time.time()
    report: dict = {"model": "SSZ_ACTION_DERIVED_MATTER_PROBE_V1",
                    "declared_before_run": True,
                    "claims": "NONE (robustness study of a model-level verdict)",
                    "frozen_reference": FROZEN}
    gates: dict = {}
    results: dict = {}

    # ---------------- flache Referenz (einmal, variantenunabhängig) --------
    fkf = FKModel(ssz=False, **GE_KW)
    uf, _ = fkf.equilibrium()
    wf, Vf = weighted_modes(fkf, uf)

    for name, kw in VARIANTS:
        t_v = time.time()
        model = ProbeModel(ssz=True, **GE_KW, **kw)
        fd = fd_check(model)
        u_eq, ginf = model.equilibrium()
        res = ge_analysis(model, u_eq, wf, Vf)
        res["fd_check"] = fd
        res["grad_inf"] = ginf
        res["wall_seconds"] = round(time.time() - t_v, 1)
        results[name] = res
        print(f"[{name}] rank: a={res['rank_by_rank']['alpha_best']:.6f} "
              f"resid={res['rank_by_rank']['max_rel_resid_after_alpha_fit']:.6f} "
              f"({res['rank_by_rank']['verdict']}) | "
              f"groups: a={res['subspace_group_matching']['alpha_best']:.6f} "
              f"resid={res['subspace_group_matching']['max_rel_resid_after_alpha_fit']:.6f} "
              f"({res['subspace_group_matching']['verdict']}, "
              f"nG={res['subspace_group_matching']['n_groups']}, "
              f"minCos={res['subspace_group_matching']['min_principal_angle_cos']:.3f}) "
              f"| FD g={fd['grad_rel_err']:.2e} H={fd['hessian_rel_err']:.2e} "
              f"| {res['wall_seconds']}s")

    report["variants"] = results

    # ---------------- Gates ------------------------------------------------
    base = results["baseline"]
    d_a = abs(base["rank_by_rank"]["alpha_best"] - FROZEN["alpha_best"])
    d_r = abs(base["rank_by_rank"]["max_rel_resid_after_alpha_fit"]
              - FROZEN["max_rel_resid_after_alpha_fit"])
    gates["P0_reproduce_frozen"] = {
        "delta_alpha": d_a, "delta_resid": d_r,
        "pass": bool(d_a < 1e-6 and d_r < 1e-6)}

    fd_ok = all(v["fd_check"]["grad_rel_err"] < 1e-5
                and v["fd_check"]["hessian_rel_err"] < 1e-5
                for v in results.values())
    gates["P1_fd_consistency"] = {"pass": bool(fd_ok)}

    ov_ok = all(v["subspace_group_matching"]["n_groups"] >= 10
                and v["subspace_group_matching"]["min_principal_angle_cos"] > 0.9
                for v in results.values())
    gates["P2_matching_well_defined"] = {
        "pass": bool(ov_ok),
        "note": ("Gruppen-Matching: >= 10 Gruppen, Hauptwinkel-cos > 0.9. "
                 "Einzelvektor-Hungarian ist bei Frequenzentartung ill-posed "
                 "und wird nur als Diagnose gefuehrt (nicht gegatet).")}

    verdicts = {name: {m: v[m]["verdict"] for m in ("rank_by_rank",
                                                    "overlap_matching",
                                                    "subspace_group_matching")}
                for name, v in results.items()}
    report["verdicts"] = verdicts

    not_collapsible_everywhere = all(
        v[m]["verdict"] == "NEW_STRUCTURE_NOT_COLLAPSIBLE"
        for v in results.values()
        for m in ("rank_by_rank", "overlap_matching",
                  "subspace_group_matching"))
    gates["P3_ge_robustness"] = {
        "pass": bool(not_collapsible_everywhere),
        "verdict": ("GE_ROBUST_AGAINST_PROPER_LENGTH_AND_ACTION_INERTIA"
                    if not_collapsible_everywhere else
                    "GE_MODEL_CHOICE_DEPENDENT — Nichtkollapsibilität ist "
                    "eine Eigenschaft der aktuellen Modellwahl"),
        "note": ("pass heisst: alle 4 Varianten x beide Zuordnungen behalten "
                 "NEW_STRUCTURE_NOT_COLLAPSIBLE. Ein Kollaps einer Variante "
                 "widerlegt die Messung nicht, zeigt aber Modellabhaengigkeit.")}

    report["gates"] = gates
    report["wall_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(report, indent=1) + "\n")
    print("\nGATES:")
    for g, v in gates.items():
        print(f"  {g}: {'PASS' if v['pass'] else 'FAIL'}"
              + (f" — {v.get('verdict')}" if v.get("verdict") else ""))
    print(f"\nwrote {OUT} ({OUT.stat().st_size} B, {report['wall_seconds']}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
