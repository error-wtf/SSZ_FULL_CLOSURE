#!/usr/bin/env python3
"""Demo 3 — Ehrliches Negativ: PEAKS_NOT_ROBUST in den Spektralgewichten.

Wissenschaftliche Ehrlichkeit heisst: Negative klar dokumentieren. Diese
Demo zeigt das zentrale ehrliche Negative des Repos live aus dem Artefakt:

  (a) Auf allen 4 Grids (n=600/1000 x Dirichlet/Robin) sind die Spectral
      Weights verschieden — maxG fällt mit der Auflösung 0.05696 -> 0.03417.
  (b) Null robuste Peaks über allen Grids: verdict PEAKS_NOT_ROBUST.
  (c) Was daraus NICHT folgt: irgendein Resonanz-/GW-Claim. Das Artefakt
      existiert als dokumentierter Fehlschlag, nicht als Erfolg.

Ausführen:  python demos/demo_3_peaks_not_robust.py   (~5 s)
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"

GRIDS = [
    ("n600_dirichlet", "n=600 Dirichlet"),
    ("n600_robin", "n=600 Robin"),
    ("n1000_dirichlet", "n=1000 Dirichlet"),
    ("n1000_robin", "n=1000 Robin"),
]


def main() -> None:
    d = json.loads((ART / "SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json").read_text())

    print("=" * 72)
    print("DEMO 3 · Spectral Weight Diagnostics — das ehrliche Negative")
    print("=" * 72)

    print(f"\n[1] Die vier Grids, live aus dem Artefakt")
    print(f"    {'Grid':<20}{'g_max':>10}{'g_median':>11}{'n_peaks':>9}")
    for key, label in GRIDS:
        v = d["variants"][key]
        print(f"    {label:<20}{v['g_max']:>10.5f}{v['g_median']:>11.5f}{v['n_peaks']:>9}")
    g1 = d["variants"]["n600_dirichlet"]["g_max"]
    g2 = d["variants"]["n1000_robin"]["g_max"]
    print(f"\n[2] maxG sinkt mit der Auflösung: {g1:.5f} -> {g2:.5f}")
    print(f"    (Δ = {(g1-g2)/g1*100:.1f} %) — die Gewichte konvergieren NICHT")
    print(f"    gegen eine grid-unabhängige Resonanzstruktur.")

    print(f"\n[3] Verdict des Artefakts: {d['verdict']}")
    print(f"    n_peaks = 0 auf jedem einzelnen Grid und null Peaks, die auf")
    print(f"    allen vieren gleichzeitig über Schwellwert lägen.")

    print(f"\n[4] Warum das im Repo bleibt:")
    print(f"    - Es grenzt den Claim-Raum ab: ohne robuste Peaks kein")
    print(f"      Spektral-Claim gegen GW-Daten. Der Pfad ist verbaut und")
    print(f"      dokumentiert, damit ihn niemand erneut beschreitet.")
    print(f"    - Kontrast: die FK-Gates (Demo 1) bestehen auf SELBEM Operator-")
    print(f"      export — die Maschine kann Peaks finden, findet hier aber keine.")
    print(f"      Damit ist PEAKS_NOT_ROBUST eine Aussage über die Daten, nicht")
    print(f"      über die Pipeline.")

    print(f"\nFazit: ein zertifizierter Fehlschlag als eigenständige Kategorie.")
    print(f"Kein astrophysikalischer Claim; effective-model scope.")


if __name__ == "__main__":
    main()
