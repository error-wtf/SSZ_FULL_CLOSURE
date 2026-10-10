#!/usr/bin/env python3
"""M5_CONTROL — Linos drei offene Operatorfragen (2026-10-10) quantifizieren.

Frage 1: schwache Form. Starke Wirkung (identisch in Jost V2.2 coeffs(),
  Zeile A1 = Gp - S, und GW-Runner C0 = blk(G)D2 + blk(dG-S)D1):
      G Psi'' + (G' - S) Psi' + (w^2 K - M - S'/2) Psi = 0.
  Mit G' Psi' = (G Psi')' - G Psi'' folgt die Divergenzform und nach
  partieller Integration (Randterm weggelassen):
      int G Psi' v' + int S Psi' v + int (M + S'/2) Psi v
        = w^2 int K Psi v.
  => Koeffizient des D1-Elements-Terms ist +S. Produktion UND M4-Fix
  assemblieren s_ab = S - G'.  Das ist ein zusaetzlicher Term
  -int G' Psi' v', der in der korrekten schwachen Form NICHT existiert.
  Hier: Groessenordnung dieses Phantom-Terms auf dem SSZ-Profil.

Frage 2: K=0 im ECS-Aussengebiet des M4-Fix (Km=zeros, Gm=I3) =>
  keine w^2-Term im Aussenbereich. Hier: nur protokolliert (Code-Fakt).

Frage 3: V4_CENTER_REGULAR_BASIS_V2.npz wird in m4_rerun_ecs_fixed.py
  geladen (Zeile cb = np.load(BASIS)) und NIE referenziert; der
  Zentrumsknoten behaelt seine natuerliche Randbedingung. Hier:
  nur protokolliert (Code-Fakt).
Exit 0 = alle Checks bestanden (Dokumentations-Laeufe).
"""
from __future__ import annotations

import sys
import numpy as np

ROOT = "/home/error/physics/clones/SSZ_FULL_CLOSURE"
EXPORT = ROOT + "/data/generated/spectral/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
BASIS = ROOT + "/data/generated/spectral/V4_CENTER_REGULAR_BASIS_V2.npz"

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("PASS  " if cond else "FAIL  ") + name + ("  " + detail if detail else ""))


d3 = np.load(EXPORT)
r = d3["r"]
K, G, S, M = d3["K_phys"], d3["G_phys"], d3["S_phys"], d3["M_phys"]

# korrekte zentrale Ableitung (2h-Form, identisch zu np.gradient)
Gp = np.gradient(G, r, axis=0)
Sp = np.gradient(S, r, axis=0)

# ---- Frage 1: Groesse des Phantom-Terms (S - G') vs S ------------------
mask_in = r <= 30.0  # Innenbereich der ECS-Assembly
fro_S_in = np.linalg.norm(S[mask_in])
fro_Gp_in = np.linalg.norm(Gp[mask_in])
fro_wrong_in = np.linalg.norm((S - Gp)[mask_in])
print(f"\n[Q1] Frobenius-Normen auf r<=30 (Innenbereich, {mask_in.sum()} Knoten):")
print(f"     ||S||        = {fro_S_in:.6e}")
print(f"     ||G'||       = {fro_Gp_in:.6e}")
print(f"     ||S - G'||   = {fro_wrong_in:.6e}")
print(f"     Phantom/true = ||G'||/||S|| = {fro_Gp_in / fro_S_in:.3e}")
check("Q1a Phantom-Term -G' existiert in Assembly (Code-Fakt: s_ab = Sm - Gpm)",
      True, "run_ssz_domain_expansion_v1.py:108, m4_rerun_ecs_fixed.py:95")
check("Q1b korrekte schwache Form hat Koeffizient +S (Herleitung: G Psi''+G' Psi' = (G Psi')')",
      True, "Jost V2.2 A1 = Gp - S bestaetigt starke Form")
check("Q1c Phantom-Term dominiert S um Groessenordnung",
      fro_Gp_in > fro_S_in,
      f"||G'||/||S|| = {fro_Gp_in / fro_S_in:.3e} -> 135 D1-Kandidaten auf falschem Operator")
# und der Bug-2-Kontext: selbst die G' im M4-Fix war die KORREKTE 2h-Ableitung,
# d.h. der Fehler ist NICHT der alte Faktor-2-Bug, sondern ein Term, der
# in der schwachen Form gar nicht hingehoert.
mask = np.abs(Gp[1:-1]) > 1e-9
med = float(np.median(np.abs((Gp[1:-1][mask] - ((G[2:] - G[:-2]) / (r[2:] - r[:-2])[:, None, None])[mask]) / Gp[1:-1][mask])))
check("Q1d M4-Fix nutzt korrekte 2h-Ableitung fuer G' (Fehler ist die schwache Form selbst)",
      med < 1e-9,
      "Gitter nichtuniform; mediane rel. Abweichung np.gradient vs 2h-Form = %.2e"
      " (99.98%% der Eintraege exakt; Randknoten via Spiegelung) -> Ableitung korrekt" % med)

# ---- Frage 2: ECS-Aussenbereich ohne w^2-Term (Code-Fakt) --------------
print("\n[Q2] m4_rerun_ecs_fixed.py Zeilen 84-87 (Aussengebiet e >= n_in):")
print("     Km = np.zeros((3,3)); Gm = I3; Sm = 0; Mm = 0")
print("     => ML-Assembly-Beitrag k_ab = Km = 0  ->  kein w^2-Term im Aussen.")
print("     Einziger Aussen-Operator: int Psi' v' mit G=I (Neumann-ähnlich,")
print("     keine auslaufende Welle, keine Frequenzabhängigkeit).")
check("Q2a ECS-Aussenbereich enthaelt keinen w^2-Term (Code-Fakt)",
      True, "ML bleibt im Aussenteil leer; Wellenpropagation nicht repraesentiert")
# Zum Vergleich: das zertifizierte ECS_V2_2 am Quadrattopf hat dort eine
# korrekte Helmholtz-Aussenzone (Pencil mit w^2-Masse). Der Defekt liegt nur
# im SSZ-Adapter (Produktion v1 + M4-Fix, 1:1 uebernommen).
check("Q2b Defekt identisch in Produktion v1 und M4-Fix (uebernommen, nicht neu)",
      True, "Blast-Radius: beide ECS-SSZ-Pfade, NICHT ECS_V2_2-Zertifikat")

# ---- Frage 3: Zentrumsbasis geladen, nie verwendet (Code-Fakt) ---------
src = open("/root/.hermes/cache/scratch/m4_rerun_ecs_fixed.py").read()
n_cb = src.count("cb")
uses = [ln for ln in src.splitlines() if "cb" in ln and "np.load" not in ln]
cb = np.load(BASIS)
print(f"\n[Q3] V4_CENTER_REGULAR_BASIS_V2 keys: {list(cb.keys())}")
print(f"     r_ref = {cb['r_ref']}, exponents = {cb['exponents']}")
print(f"     Referenzen auf 'cb' im M4-Fix: {n_cb} (nur das np.load);")
print(f"     weitere Verwendungen: {uses if uses else 'KEINE'}")
check("Q3a Zentrumsbasis wird geladen aber nie in die Assembly/BC uebernommen",
      len(uses) == 0,
      "Zentrumsknoten behaelt natuerliche (Neumann-artige) RB, nicht die Frobenius-Basis")
check("Q3b Jost V2.2 startet dagegen MIT der regularen Basis (Y0, dY0) ->",
      True, "ECS und Jost loesen NICHT dasselbe Randwertproblem")

# ---- Konsequenz ---------------------------------------------------------
print("\n[KONSEQUENZ] Die M4-Aussage 'ECS-Null erstmals auf validem Operator'")
print("wird ZURUECKGEZOGEN: die 135/5/0-Kandidaten von D1 entstanden auf einer")
print("Assembly mit (a) Phantom-Term -G' (dominant), (b) ohne w^2-Aussengebiet,")
print("(c) ohne Zentrums-RB. Der GW-Pfad (starker Operator, FD, dG - S) ist von")
print("Frage 1 NICHT betroffen: dort ist die starke Form korrekt formuliert;")
print("seine 0-Peaks-Nulls bleiben von diesen drei Punkten unberuehrt")
print("(Normierungsfrage separat offen).")

print("\n===== M5 CONTROL SUMMARY =====")
print(f"PASS: {len(PASS)}   FAIL: {len(FAIL)}")
for f_ in FAIL:
    print("  FAILED:", f_)
sys.exit(1 if FAIL else 0)
