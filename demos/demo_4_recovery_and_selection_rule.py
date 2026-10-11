#!/usr/bin/env python3
"""Demo 4 — End-to-End-Kette: blinde Recovery + Selektionsregel B = 0.

Die Pipeline-Kette "Validierte Maschine -> ehrliche Projektion":

  (a) Spectral Response Framework (Stage 3): 3 Moden werden BLIND in ein
      Spektrum injiziert und vom Code gefunden, der die Frequenzen nicht
      kennt — Recovery 3/3, Fehler < 1 Frequenz-Bin, Null-Test über 200
      Rausch-Trials: median max-z = 0.0 (keine falschen Peaks).
  (b) Weisz Response Bridge (Run 6): Dieselbe Maschine auf dem gefrorenen
      SSZ-Operator. Projektion B (ch0,ch2) liefert g_max EXAKT 0.0 —
      keine numerische Nähe zu null, sondern strukturelle Orthogonalität
      (Selektionsregel). A: 0.0930, C: 3e-6.

Der Punkt: die Maschine ist kalibriert (a), und ihre Aussagen auf dem
Operator sind deshalb interpretierbar (b) — inkl. der exakten Null.

Ausführen:  python demos/demo_4_recovery_and_selection_rule.py   (~5 s)
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"


def main() -> None:
    srf = json.loads((ART / "SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json").read_text())
    bridge = json.loads((ART / "SSZ_WEISZ_RESPONSE_BRIDGE_V1.json").read_text())

    print("=" * 72)
    print("DEMO 4 · Blinde Recovery + Selektionsregel — die Maschine und ihre Null")
    print("=" * 72)

    br = srf["detail"]["G2_G3_G4"]["g3_blind_recovery"]
    print(f"\n[1] Blinder Injection/Recovery-Test (Stage 3, G3)")
    print(f"    {'Mode':<7}{'injiziert':>12}{'blind gefunden':>17}{'Fehler':>10}")
    for i, (wi, wr, e) in enumerate(zip(br["injected_ws"], br["recovered_ws"], br["freq_errors"])):
        ok = e < br["bin_width"]
        print(f"    {i+1:<7}{wi:>12.5f}{wr:>17.5f}{e:>10.5f}  {'< bin ✓' if ok else 'ÜBER bin ✗'}")
    print(f"    bin_width = {br['bin_width']} — {br['n_recovered']}/{br['n_injected']} Moden, pass = {br['pass']}")

    null = srf["detail"]["G2_G3_G4"]["g2_noise_null"]
    print(f"\n[2] Null-Test (G2): {null['n_trials']} Rausch-Trials ohne Injektion")
    print(f"    median max-z = {null['median_max_z']}, max max-z = {null['max_max_z']},")
    print(f"    median n_peaks = {null['median_n_peaks']}  -> pass = {null['pass']}")
    print(f"    Die Maschine erfindet keine Peaks.")

    fa = srf["detail"]["G2_G3_G4"]["g4_false_alarm"]
    print(f"\n[3] False-Alarm-Kontrolle (G4): SNR-{fa['injection_snr']:.0f}-Injektion")
    print(f"    andere Peaks über Schwellwert: {fa['n_other_peaks']},")
    print(f"    max SNR anderer Peaks: {fa['max_other_snr']}  -> pass = {fa['pass']}")

    print(f"\n[4] Dieselbe Maschine auf dem gefrorenen Operator (Bridge Run 6)")
    proj = bridge["stage2_ssz_projections"]["projections"]
    print(f"    {'Projektion':<18}{'Kanäle':<14}{'g_max':>12}{'g_median':>12}")
    for key, label, ch in [("A_ch0_ch0", "A", "(ch0,ch0)"),
                           ("B_ch0_ch2", "B", "(ch0,ch2)"),
                           ("C_ch1_ch1_late", "C", "(ch1,ch1 late)")]:
        p = proj[key]
        print(f"    {label:<18}{ch:<14}{p['g_max']:>12.6g}{p['g_median']:>12.6g}")
    print(f"\n[5] Der Punkt: Projektion B ist EXAKT 0.0")
    print(f"    g_max(B) = {proj['B_ch0_ch2']['g_max']} — nicht 1e-9, nicht 'klein',")
    print(f"    sondern exakt null: strukturelle Orthogonalität (Selektionsregel),")
    print(f"    die die validierte Maschine aus (1)-(3) reproduziert.")
    print(f"    A: g_max = {proj['A_ch0_ch0']['g_max']:.4f} — die Projektions-")
    print(f"    abhängigkeit IST das Ergebnis (nicht Unterdrückung durch die Pipeline).")

    print(f"\nFazit: kalibrierte Messkette + interpretierbare Null. Effective-model")
    print(f"scope: Software-Verifikation, kein GW-Detektions-Claim.")


if __name__ == "__main__":
    main()
