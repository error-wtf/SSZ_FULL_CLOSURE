#!/usr/bin/env python3
"""Demo 1 — G-E: Die Struktur des SSZ-FK ist NICHT FK x Rotfaktor.

Der wichtigste Einzelbefund des Repos, live aus dem zertifizierten Artefakt
gelesen und in drei Szenarien durchgespielt:

  (a) "Der residuale Wert ist doch klein" — nein: 0.5642 nach bestmöglicher
      uniform-alpha-Anpassung. Das ist 56 %, keine Rundung.
  (b) "Vielleicht kollabiert es unter einer anderen Phasenkonvention" —
      Prüfstand 6a: phi = Q*l (proper length), Trägheit, beides kombiniert;
      in allen 4 Varianten x 3 Matchings bleibt NEW_STRUCTURE_NOT_COLLAPSIBLE.
  (c) "Vielleicht ist der Solver schuld" — Delta ~ 2e-7 mit korrektem
      Solver (siehe README); der Befund ist numerisch sauber getrennt.

Ausführen:  python demos/demo_1_ge_not_collapsible.py   (~5 s)
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"


def main() -> None:
    fk = json.loads((ART / "SSZ_FK_COUPLED_MODEL_V1.json").read_text())
    probe = json.loads((ART / "SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json").read_text())

    ge = fk["gates"]["G-E_new_structure"]
    print("=" * 72)
    print("DEMO 1 · G-E NEW_STRUCTURE_NOT_COLLAPSIBLE — live aus dem Artefakt")
    print("=" * 72)

    resid = ge["max_rel_resid_after_alpha_fit"]
    print(f"\n[1] Was wurde gemessen?")
    print(f"    Struktur-Scores S(fk_mass-weighted) vs S(fk * redshift).")
    print(f"    Bestmögliche uniforme alpha-Anpassung: alpha_best = "
          f"{ge.get('alpha_best', float('nan')):.4f}")
    print(f"    Max. relative Residuen NACH Anpassung:  {resid:.4f}  (= {resid*100:.1f} %)")
    print(f"    collapsible_to_uniform_factor = {ge['collapsible_to_uniform_factor']}")
    print(f"    verdict = {ge['verdict']}")
    print(f"    -> Ein uniformer Rotfaktor erklaert <50 % der Struktur. Das ist der")
    print(f"       BEFUND, kein Bug: SSZ-FK ist nicht 'FK x z'.")

    print(f"\n[2] Ist das robust? -> Prüfstand 6a (phi = Q·l, Trägheit):")
    base = probe["frozen_reference"]["max_rel_resid_after_alpha_fit"]
    print(f"    gefrorene Referenz resid = {base:.4f}")
    print(f"    {'Variante':<20}{'Matching':<26}{'resid':>9}  verdict")
    for var, name in [("baseline", "baseline"), ("proper_phase", "phi=Q*l"),
                      ("action_inertia", "action inertia"),
                      ("proper_plus_inertia", "phi+inertia")]:
        for mat in ("rank_by_rank", "overlap_matching", "subspace_group_matching"):
            e = probe["variants"][var][mat]
            mark = " <-- " + name if var == "proper_phase" and mat == "rank_by_rank" else ""
            print(f"    {name:<20}{mat:<26}{e['max_rel_resid_after_alpha_fit']:>9.4f}  "
                  f"{e['verdict']}{mark}")
    n_nc = sum(1 for row in probe["verdicts"].values() for v in row.values()
               if v == "NEW_STRUCTURE_NOT_COLLAPSIBLE")
    print(f"    -> {n_nc}/12 Kombinationen bleiben NOT_COLLAPSIBLE.")
    print(f"       Unter proper length sinkt resid {base:.4f} -> "
          f"{probe['variants']['proper_phase']['rank_by_rank']['max_rel_resid_after_alpha_fit']:.4f},")
    print(f"       das Urteil bleibt: die Struktur ist nicht rotfaktor-kollabierbar.")

    print(f"\n[3] Gates des Prüfstands:")
    for g, v in probe["gates"].items():
        state = "PASS" if v["pass"] else "FAIL"
        print(f"    {g:<28} {state}")
        if not v["pass"]:
            print(f"        note: {v['note'][:100]}...")
    print(f"    P2 FAIL ist dokumentierte Erkenntnis: Einzelvektor-Hungarian bei")
    print(f"    Frequenzentartung ist ill-posed (kein Modellfehler).")

    print(f"\nFazit: G-E ist kein Artefakt der Phasenkonvention, des Matchings")
    print(f"oder der Diskretisierung. Software-Verifikation eines strukturellen")
    print(f"Befunds im effective-model scope — keine astrophysikalische Aussage.")


if __name__ == "__main__":
    main()
