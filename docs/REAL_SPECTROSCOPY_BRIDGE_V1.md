# REAL_SPECTROSCOPY_BRIDGE_V1 — Implementation Contract (2026-10-05)

Status: ACTIVE WORK ORDER. This document is the binding spec for the bridge
layer between mathematical mode spectroscopy and real observational data.
Nothing here retunes, fits or scales any model. Every stage has explicit
blindness rules and freeze points.

## 0. Context (verified 2026-10-05)

- Working branch: `spectroscopy-real-data-20261004` (forks on
  `qnm-global-final-20261002@5c78ea3`).
- Local principal spectroscopy (v3) passes all controls in the healthy window;
  global coupled QNM claims remain BLOCKED by the core `K_scalar` pocket
  (stencil-, route- and formula-audited; see
  `docs/P1_KSCALAR_ACTION_LEVEL_VERIFICATION_2026-10-05.md`).
- Pinned blind catalogs already in repo (SHA256-frozen):
  - `data/observations/nicer/MAXI_J1820_PLUS_070_CATALOG_V1.json`
  - `data/observations/eso/GALACTIC_CENTER_CATALOG_V1.json`

## 1. Three separate outputs (never merged before stage 3)

```
RAW OBSERVATIONS  --stage A-->  OBSERVED_MODE_CATALOG.json   (blind)
HEALTHY SSZ OP    --stage B-->  PREDICTED_MODE_CATALOG.json (blind)
FROZEN CATALOGS   --stage C-->  statistical comparison only
```

Blindness rules (binding):
- B1: Stage B code must not read stage A outputs. Enforced by CI directory
  separation (`data/observed/` vs `data/predicted/`) + script asserts.
- B2: No frequency may be rescaled after stage C begins. The only shared
  degrees of freedom are calibration/background nuisance parameters, and the
  null model (GR / pure continuum) gets the SAME ones.
- B3: Catalog freeze = SHA256 recorded in `FROZEN_CATALOG_INDEX.json`.
  Any change afterwards = new version, never overwrite.

## 2. Stage A — blind observed-mode catalog (NICER MAXI J1820+070)

### A1. Product ingest
- Download event files (HEASARC, obsids from frozen catalog v1: 20 obs, 25.9 ks).
- Standard screening (`nicerl2` equivalent), GTIs per obs.
- Record per obs: exposure, GTI gaps > 1 s list, barycenter correction version.

### A2. Timing/spectral extraction (model-free)
- Power spectral density per obs (Leahy normalization), frequency range
  0.01–100 Hz, geometric binning.
- Energy-resolved PSDs: 3 bands straddling Fe-K interest (e.g. 0.5–2, 2–5,
  5–10 keV) — band choice frozen BEFORE looking at any PSD.
- Averaged cross spectra for coherence/phase-lag (2 subbands max).
- continua: per-obs PCA-style background estimate via systematic model, NOT
  a fit shared with any model prediction.

### A3. Peak candidate extraction (blind rules)
- Candidate = local maximum with global signifiance > 4 sigma (trial-corrected)
  AND persistence in >= 3 obs or one obs with > 8 sigma.
- For each candidate record: f_n, FWHM/Gamma_n, A_n (rms), covariance matrix,
  significance, per-obs occurrence, energy dependence.
- MANDATORY negative control: the known instrumental ~55 Hz GTI-gap artifact
  (documented by NICER analysis tips for this target) must be RECOVERED by the
  same pipeline and then EXPLICITLY REJECTED with reason `GTI_GAP_ARTIFACT`.
  If the pipeline does not see the 55 Hz peak where GTI gaps predict it, the
  pipeline is broken — block, do not proceed.
- Output: `data/observed/OBSERVED_MODE_CATALOG.json` + SHA256 freeze.

## 3. Stage B — predicted-mode catalog (healthy window only)

- Input: frozen local principal operator (same member hash as v3 run).
- For each fixed observable q (the 5 registered ones) and each radial sample:
  - K-normalized modes, cluster residues, top-k fractions, N_eff, entropy;
  - NEW (Li-type state diagnostics, this work order):
    `p_{n,i} ∝ w_i psi_n^dag(r_i) K(r_i) psi_n(r_i)`,
    `IPR_n = sum_i p_{n,i}^2`, finite-size scaling dimension via L in
    {6,12,20,42,110,420,1000} (existing grid);
    classification: localized / extended / multifractal per mode.
- Global QNM poles: NOT INCLUDED (blocked, physical_qnm_claim_allowed=false).
- Output: `data/predicted/PREDICTED_MODE_CATALOG.json` + SHA256 freeze.

## 4. Stage C — statistical comparison (frozen vs frozen)

- Match statistic: pre-registered before first comparison, written into
  `schemas/COMPARISON_PREREGISTRATION.json`:
  - primary: frequency residual |f_obs - f_pred| / sigma_obs with null model
    (red-noise continuum + detector lines) evaluated with identical nuisance;
  - secondary: residue-class agreement (observed A_n ranking vs predicted
    Z_n ranking), permutation p-value;
  - multiplicity: Benjamini-Hochberg, q < 0.05.
- Verdict classes: NO_MATCH / WEAK / STRONG / FALSE_POSITIVE_CONTROL_FAILED.
- The 55 Hz artifact enters as a SANITY row: pipeline must classify it as
  instrumental. Failure = all results invalid.

## 5. Engineering test order

1. Skeleton schemas (this commit).
2. A1/A2 on the FIRST obsid only (5200120403) — smoke test.
3. A3 blind extraction on that obs; verify 55 Hz GTI prediction machinery.
4. Full 20-obs pass; freeze OBSERVED catalog.
5. B: Li-type diagnostics extension in v3; freeze PREDICTED catalog.
6. C: comparison with pre-registered statistics.

## 6. What this deliberately does NOT do

- No global QNM matching (blocked by core K_scalar pocket).
- No photon-line vs eigenfrequency identity claims (needs transfer function;
  separate later layer `TRANSFER_BRIDGE_V1`).
- No GW ringdown matching (waits for certified Jost/ECS global operator).

---

## UPDATE 2026-10-05 (Commit-Serie): Vertrag implementiert

Die 8 Gates des Vertrags sind jetzt als Pakete/Tools umgesetzt:

| Gate | Umsetzung |
|------|-----------|
| 1 Provenienz | `src/ssz_p5/observations/provenance.py` (sha256 je Input + git head + `write_frozen` mit sidecar) |
| 2 Healthy-domain | unverändert: globaler QNM bleibt hinter `require_direct_krgm_certificate` fail-closed |
| 3 Weisz-Response | unverändert: Weisz-Repro + `compare_ssz_weisz1978_selectivity.py` Pflichtregression |
| 4 Li-Lokalisierung | `src/ssz_p5/qnm/localization_diagnostics.py` (IPR + D2-Blockscaling) |
| 5 Blind observed | `tools/bridge_stage_a_observed_modes.py` + `observations/timing.py` |
| 6 Blind predicted | `src/ssz_p5/spectroscopy/theory_catalog.py` (Blindheits-Assert + FrozenMassPrior) + `physical_units.py` (f = c³/(4πGM)·Re(ω̄), Mass-Prior extern+gehasht) |
| 7 Frozen compare | `tools/compare_frozen_mode_catalogs.py` (sha256-sidecar-Pflicht; Assignment-Permutations-Null; MATCH_ABOVE_NULL/NO_MATCH_ABOVE_NULL) |
| 8 MODE_CERTIFICATE | `src/ssz_p5/qnm/mode_certificate.py` + `tools/run_mode_certificates.py` + `schemas/MODE_CERTIFICATE.schema.json` |

Negativtest (Vertrag): `tests/unit/test_nicer_gap_artifact_control.py` —
reproduziert das dokumentierte 55-Hz-GTI-Lücken-Artefakt (ObsID 1200120107,
MAXI J1820+070) synthetisch: Kandidat sitzt auf der Gap-Kamm-Frequenz UND
verschieden in kontinuierlichen Segmenten → `SPURIOUS_GTIL_GAP_ARTIFACT_REJECTED`;
Kontrollsinal (echte 20-Hz-Oszillation) überlebt dieselbe Pipeline.

Statistik-Verifikation des Null-Tests: 3 echte Matches → p=0.003
MATCH_ABOVE_NULL; unkorrelierter Zufall → p=1.0 NO_MATCH_ABOVE_NULL.
