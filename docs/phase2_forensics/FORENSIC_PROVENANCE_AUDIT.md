# FORENSIC PROVENANCE / REPRESENTATIVE AUDIT — Resolution of the C1/C3 vs historical-K contradiction

Date: 2026-09-22. Scope: forensics only. No new theory, no completion, no
stability optimization, no fitting. Triggered by the owner's identified
contradiction: preregistered C1/C3 streams show min K = -18.522265817831975
while the canon records healthy historical witnesses (K_even_min ~ 10.5415,
c_r values, c_Omega 84.912/643.142, strong-H carrier).

## 1. What the historical numbers actually were (provenance)

- Source: data/production/ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv
  (4000 rows, u in [0.61, 0.715]), recorded 2026-09-16 in
  SSZ_P5_COMPLETE_STATUS_TO_KRGM_QNM (section 16) and the Monograph/HSVT memory.
- Measured quantities:
    K_ZK_exact            min = 10.541549191507217   <-- the famous 10.54
    K_from_e1_identity    min = 10.541549191460803 (residual ~4.6e-11)
    K_even_corrected      min = 0.963915063812136    <-- internal column, lower!
    c_r,T^2 min = 1.0050694018534372; c_r,S^2 min = 0.9037451776188076;
    c_r,V^2 min = 9.458870488366266; cOmega_-^2 min = 84.912; cOmega_+^2 min = 643.142.
- These are REDUCED even-sector ZK-exact combinations on the central exact SVT
  representative (ZK regression 2026-09-13; reducer
  ssz_p5_profile_operator_reducer_JET9D8_2026-09-16), NOT the full finite-L
  Schur minimum eigenvalue of the complete 3x3 reduced kinetic matrix.
- f3_X provenance: explicitly integrated from the scalar background ODE
  (ScalarODEIdentity, A q' + B q + C = 0, q = f3X); lobe scalar residual
  ~1e-12; CSV column f3X_integrated.

## 2. Forensic measurements (this audit)

2.1 Representative identity:
    absolute_reemit's algebraic corrector changes the central member's own
    on-shell data by at most 6.9e-12 (f2), 1.2e-11 (f2X), 0.0 (f2F) scaled.
    => The C1/C3 member stream IS the central member data. No representative
    fork exists in the construction. (Class K "regression error": refuted.)

2.2 f3_X audit (directive section 5):
    The central frame carries nontrivial f3X in [-18.494, 44.728]; scalar
    identity residual on the central frame = 2.79e-4; after the corrector =
    3.43e-4 (recorded; historical lobe level was 1e-12 — small degradation
    from re-solved f2-family + numerical derivatives, NOT a dropped structure).
    C1/C3 retained f3X. (Class E "omitted background-determined jets": refuted.)

2.3 Section-9 decisive replay (historical member, CURRENT emitter + reducer):
    central frame, own on-shell f2 (corrector is a no-op), canonical ZK emitter,
    current kinetic_schur reducer, finite-L minimum eigenvalues:
        L=6: -1.61e-01 ... L=1000: -2.78e+01 (worst near u~0.7099).
    The historical healthy lobe values do NOT survive as FULL finite-L Schur
    minima. This is consistent with — and was already recorded by — the
    2026-09-17 regional kinetic gate (commit 9c85a30): central region FAIL,
    eigenvalues -43.86167 / 0.09070 / 115.42431 at u=0.69, L=6; and
    ALL_EVIDENCE_RUN.json: Central Full-SVT kinetic gate FAIL, L=6 min(K) =
    -692.14. (Class J: historical positive K superseded as a full finite-L
    stability witness by the current reducer's stricter full-matrix gate.)

2.4 Background residual semantics (directive section 0 / 14):
    absolute_reemit's own residuals on the corrected frame: E00 = 4.37e-12,
    E11 = 6.30e-12, JA = 9.90e-15 (machine precision). The Phase-2 diagnostic
    values E00 = 1.04 / E11 = 0.043 were an EVALUATOR INPUT ARTIFACT: they were
    measured on the emitted 41-slot stream, which does not carry the primitive
    columns (f2F etc. default to 0 in the evaluator). BACKGROUND = PASS is
    justified; the Phase-2 caveat note is corrected by this audit. Tolerances
    unchanged.

2.5 C3 stream identity explanation (directive section 8):
    Within the frozen bundle basis B_k = b(u) z^k, EVERY nonzero coefficient
    vector produces a nonzero third derivative (the Gram matrix of the third
    derivatives is positive definite: eigenvalues 1.54e9 ... 4.30e10). Hence
    "minimise ||D'''||^2" has the unique minimiser c = 0 -- the FULL zero
    displacement. Setting the third covariant jet to zero therefore coincides,
    in this frozen basis, with the whole transverse displacement vanishing.
    This is a THEOREM of the frozen basis, not an implementation shortcut;
    it is exactly the nontriviality risk declared in Phase 1.5 and confirmed
    by the owner's rule ("if trivial, do not change -- that is C3 genuine
    prediction").

## 3. Required end result (directive section 15)

    CLASS D — MIXTURE, layer by layer:

    (i)  C1/C3 are valid but are NOT a new/different action: they are the
         central member's own data (machine precision). [layer B resolved:
         no representative error]
    (ii) The historical positive K (10.54) was a reduced even-sector ZK
         combination; as a full finite-L stability witness it was already
         superseded on 2026-09-17 by the full-matrix kinetic gate. [layer J]
    (iii) The old angular values 84.912/643.142 belong to the superseded
         reduced/spline-era angular determination; the action-level direct
         epsilon^3 scan supersedes them (canon supersede list, section 16).
         [layer J, angular]
    (iv) f3_X and other background-determined structure was retained. [layer E
         refuted]
    (v)  No regression/implementation error in the current K measurement.
         [layer K refuted]

## 4. Consequence for the Phase-2 interpretation (correction of record)

The Phase-2 statement "both minimality principles predict the ghost" remains
true as a prediction statement. The Phase-2-adjacent interpretation that this
"demonstrates minimal completion class necessarily fails and new non-minimal
physics is required" is RETRACTED as overstrong:

    - Healthy finite-L K representatives DO exist inside the already-used
      Horndeski/SVT architecture: the strong-H carrier (min eig 1.6051e-07 at
      L=1000, positive at all tested L, positive radial characteristics;
      HSVT_DIRECT_REGION_AUDIT_2026-09-18, "carrier" section, all PASS) and
      the exterior region (9.8899e-07 at L=1000). [Q1 = YES]
    - What is NOT established for those representatives: same-action provenance
      (Q2) and full closure (Q3). The Sept-16 memory itself: "finite-l global
      same-action closure ... not yet passed".
    - Therefore the correct open question is unchanged and sharper: assemble a
      healthy-K representative and the SSZ background into ONE same-action
      member (Q2), then close the remaining gates (Q3). New fundamental
      dynamics is NOT logically forced by the current ghost.

## 5. Historical pass survival ledger (summary table)

| Historical PASS | Status today | Reason |
|---|---|---|
| Lobe K_ZK_exact min 10.54 (2026-09-13/16) | PASS_SCOPE_ONLY | reduced even-sector combination; full finite-L gate supersedes |
| Lobe c_r^2 (T/S/V) | PASS_SCOPE_ONLY | same reduced scope |
| Lobe c_Omega 84.912/643.142 | SUPERSEDED_NUMERICS | reduced-era angular determination; direct eps^3 determinant supersedes (canon section 16 supersede list) |
| Central finite-L kinetic | SUPERSEDED_NUMERICS (2026-09-17, 9c85a30) | full-matrix gate FAILS (L=6,12,20,42) |
| Strong-H carrier finite-L K | PASS_STILL_VALID (scope: carrier representative) | min 1.6e-7 @ L=1000, all-L positive, radial positive |
| Exterior finite-L K | PASS_STILL_VALID (scope) | 9.9e-7 @ L=1000 |
| v12 = -v6/(2h); a5 with -a1''; S-block retained | PASS_STILL_VALID (conventions) | canon supersede list section 16 |
| K-repaired dense member angular | FALSIFICATION WITNESS (direct eps^3) | cminus^2 < 0 on 297/297 |
| C1/C3 blind prediction | PASS_STILL_VALID as prediction | both select the minimal stream; K_GHOST_FAIL is their genuine prediction |
