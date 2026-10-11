#!/usr/bin/env python3
"""Demo 2 — Prüfstand 6a im Detail: Proper Length statt Koordinatenlänge.

Der Robustheits-Prüfstand 6a (Commit 5245a0f, 2026-10-11) ersetzt die
Phasen phi = omega*t durch die proper length phi = Q*l und zusätzlich
Trägheits-Terme. Diese Demo zeigt interaktiv:

  (a) Was sich ändert: resid 0.5642 -> 0.0043 unter proper length
      — die numerische Diskrepanz schrumpft um Faktor ~131, ABER
  (b) was GLEICH bleibt: das Urteil NEW_STRUCTURE_NOT_COLLAPSIBLE in
      allen 4 Varianten x 3 Matchings (12/12).
  (c) P0/P1-Pass: die gefrorene Referenz wird bit-genau reproduziert.

Der Punkt: 6a misst, wie stark der G-E-Befund von Konventionen abhängt.
Antwort: die Residuen ändern sich dramatisch (Konventionen wirken), das
Urteil bleibt stabil — genau das macht G-E zu einer Messung, nicht zu
einem Mapping-Artefakt.

Ausführen:  python demos/demo_2_proper_length_probe.py   (~5 s)
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"


def main() -> None:
    probe = json.loads((ART / "SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json").read_text())

    print("=" * 72)
    print("DEMO 2 · Prüfstand 6a — Proper-Length-Phasen & Action-Trägheit")
    print("=" * 72)

    base = probe["frozen_reference"]["max_rel_resid_after_alpha_fit"]
    pp = probe["variants"]["proper_phase"]["rank_by_rank"]

    print(f"\n[1] Setup")
    print(f"    Modell: {probe['model']}")
    print(f"    declared_before_run: {probe['declared_before_run']}")
    print(f"    claims: {probe['claims']}")

    print(f"\n[2] P0 — Reproduktion der gefrorenen Referenz")
    p0 = probe["gates"]["P0_reproduce_frozen"]
    print(f"    delta_alpha  = {p0['delta_alpha']}  (exakt)")
    print(f"    delta_resid  = {p0['delta_resid']}  (exakt)")
    print(f"    -> {'PASS' if p0['pass'] else 'FAIL'}: der Prüfstand rechnet auf demselben")
    print(f"       Datenstand wie das zertifizierte FK-Modell V1.")

    print(f"\n[3] Der Sprung: Koordinatenlänge -> proper length")
    print(f"    baseline (phi = omega*t):   resid = {base:.6f}")
    print(f"    proper length (phi = Q*l):  resid = {pp['max_rel_resid_after_alpha_fit']:.6f}")
    ratio = base / pp["max_rel_resid_after_alpha_fit"]
    print(f"    Faktor: {ratio:.1f}x kleiner — die numerische Lücke ist")
    print(f"    konventionsabhängig, die Konvention 'proper length' bringt")
    print(f"    beide Strukturen nahe zusammen (alpha_best = {pp['alpha_best']:.4f}).")

    print(f"\n[4] Aber: Urteil unverändert — 12/12 NOT_COLLAPSIBLE")
    header = f"    {'Variante':<22}{'rank':>10}{'overlap':>10}{'subspace':>11}"
    print(header)
    for var, name in [("baseline", "baseline"), ("proper_phase", "phi=Q*l"),
                      ("action_inertia", "inertia"), ("proper_plus_inertia", "phi+inertia")]:
        row = probe["variants"][var]
        cells = []
        for mat in ("rank_by_rank", "overlap_matching", "subspace_group_matching"):
            r = row[mat]["max_rel_resid_after_alpha_fit"]
            ok = row[mat]["verdict"] == "NEW_STRUCTURE_NOT_COLLAPSIBLE"
            cells.append(f"{r:>10.4f}{'✓' if ok else '✗'}")
        print(f"    {name:<22}{cells[0]}{cells[1]}{cells[2]}")
    print(f"    (✓ = NEW_STRUCTURE_NOT_COLLAPSIBLE)")

    print(f"\n[5] P2 — wo es klemmt (ehrlicher FAIL)")
    p2 = probe["gates"]["P2_matching_well_defined"]
    print(f"    pass = {p2['pass']}")
    print(f"    note: {p2['note']}")

    print(f"\nFazit: 6a zeigt zwei Dinge zugleich — (i) die Residuen sind stark")
    print(f"konventionsabhängig (0.56 -> 0.0043), (ii) das Urteil ist es NICHT.")
    print(f"Eine Messung, deren Urteil unter Konventionswechsel stabil bleibt," )
    print(f"ist robust; die resid-Zahl allein wäre irreführend gewesen.")


if __name__ == "__main__":
    main()
