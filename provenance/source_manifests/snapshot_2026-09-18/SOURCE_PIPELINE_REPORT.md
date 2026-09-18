# SSZ P5 full-pipeline report

**Release:** 2026-09-16
**Generated:** 2026-09-18T09:01:27.349111+00:00
**Overall:** PASS
**Checks:** 45/45 passed

This report validates the frozen repository. It does not repeat completed derivations.

## repo

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `README.md` | **PASS** | True | exists |
| `RESEARCH_STATE.md` | **PASS** | True | exists |
| `DO_NOT_RECOMPUTE.md` | **PASS** | True | exists |
| `STATUS_FULL_CLOSURE.md` | **PASS** | True | exists |
| `ARTIFACT_CLASSIFICATION.json` | **PASS** | True | exists |
| `PRODUCTION_BLACKLIST.json` | **PASS** | True | exists |
| `NUMERICAL_POLICY.json` | **PASS** | True | exists |
| `pyproject.toml` | **PASS** | True | exists |
| `ssz_p5_full_closure_auditor.py` | **PASS** | True | exists |
| `src/ssz_p5` | **PASS** | True | directory exists |
| `tests` | **PASS** | True | directory exists |
| `data/production` | **PASS** | True | directory exists |
| `data/regression` | **PASS** | True | directory exists |
| `data/qnm` | **PASS** | True | directory exists |
| `paper` | **PASS** | True | directory exists |
| `archive/full_working_snapshot` | **PASS** | True | directory exists |

## software

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `pytest` | **PASS** | 0 | exit code 0 |

## scientific_audit

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `full_closure_auditor` | **PASS** | 0 | exit code 0 |
| `numerical_regression_failures` | **PASS** | [] | no failures |
| `auditor_constructive_failures` | **PASS** | 0 | 0 |

## production

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `P5_global_principal_KRG` | **PASS** | data/production/ssz_p5_F2A_GLOBAL_CANONICAL_KRG_PRINCIPAL_CINF_FINAL_2026-09-15(2).csv | present |
| `outer_handover_background` | **PASS** | data/production/ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv | present |
| `inner_handover_background` | **PASS** | data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv | present |
| `strongH_selected_41` | **PASS** | data/production/ssz_p5_F3_horndeski_carrier_SELECTED_41of41_2026-09-15.csv | present |
| `core_selected_41` | **PASS** | data/production/ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv | present |
| `outer_selected_41` | **PASS** | data/production/ssz_p5_F3_outer_same_action_SELECTED_41of41_2026-09-15.csv | present |
| `ZK_raw_41` | **PASS** | data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv | present |
| `ZK_exact_regression` | **PASS** | data/production/ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv | present |
| `center_handover` | **PASS** | data/production/ssz_p5_F1b_FINAL_Cinf_center_to_punctured_handover_2026-09-14.csv | present |

## qnm_archive

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `maxwell_jost` | **PASS_ARCHIVED** | data/qnm/ssz_p5_F3_maxwell_l2_jost_convergence_2026-09-15.csv | present |
| `spectral_probe` | **PASS_ARCHIVED** | data/qnm/SSZ_P5_POST_CLOSURE_SPECTRAL_PROBE_v2_2026-09-14.md | present |
| `qnm_checkpoint` | **PASS_ARCHIVED** | data/qnm/SSZ_P5_F3_COUPLED_QNM_EXECUTION_CHECKPOINT_2026-09-15.md | present |
| `eikonal_proxy` | **PASS_ARCHIVED** | data/qnm/SSZ_CANONICAL_V13_EIKONAL_QNM_PROXY_v1_REPORT_2026-08-28.md | present |
| `historical_reconciliation` | **PASS_ARCHIVED** | data/qnm/SSZ_QNM_HISTORICAL_VS_V13_RECONCILIATION_v1_REPORT_2026-08-28.md | present |
| `maxwell_jost_rows_finite` | **PASS** | 199 | all finite |
| `maxwell_jost_parameter_coverage` | **PASS** | {'orders': [6, 8, 10], 'R': [25.0, 30.0, 35.0, 40.0, 45.0]} | >=3 Jost orders and >=3 outer radii |
| `maxwell_jost_residual_recorded` | **PASS** | 1.2851121504250217e-06 | finite residuals |
| `scalar_WKB_threshold_256` | **PASS_ARCHIVED** | 256 | archived scalar threshold |
| `maxwell_WKB_threshold_257` | **PASS_ARCHIVED** | 257 | archived Maxwell threshold |
| `outgoing_jost_order10_validation` | **PASS_ARCHIVED** | ~2.43e-12 | archived Schwarzschild-tail transport validation |
| `rejected_real_axis_poles_preserved_as_rejected` | **PASS** | True | must remain rejected, not production QNM |

## hygiene

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `production_blacklist` | **PASS** | 8 | machine-readable blacklist present |
| `production_blacklist_enforced` | **PASS** | [] | no forbidden production inputs |
| `forbidden_stream_preserved_as_diagnostic` | **PASS** | ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv |  |

## action

| Check | Status | Value | Criterion |
|---|---|---:|---|
| `frozen_action_member` | **PASS** | SSZ P5 regional Full-SVT production member | locked six-region Full-SVT production cover |

## Scientific scope

Absolute closure: **NOT_CERTIFIED**.
Coupled HSVT QNM: **NOT_CERTIFIED**.
An archival pipeline PASS is not a direct-export or coupled-spectrum certificate.

