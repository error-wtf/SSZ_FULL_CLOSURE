# TEST CATALOG — complete semantic inventory

Machine-generated from TEST_SEMANTICS.json (tools/generate_test_catalog.py).
Every one of the **166 collected tests** is documented with its
scientific claim, category, independence level, member scope and gate.

## Headline report

| Metric | Count |
|---|---|
| Total tests (all PASS) | 166 |
| STRONG (independent) | 40 |
| MEDIUM (partially independent / self-consistency) | 105 |
| TAUTOLOGICAL (change-pins; do not count as physics evidence) | 9 |
| HISTORICAL (zero-vector member scope) | 12 |
| Negative controls | 10 |
| Provenance tests | 25 |

## Categories

| Category | Tests |
|---|---|
| PROVENANCE | 25 |
| SVT_OPERATOR | 25 |
| MH_EQ85 | 18 |
| KINETIC | 10 |
| SOFTWARE_INFRASTRUCTURE | 10 |
| NEGATIVE_CONTROL | 10 |
| SYMBOLIC_IDENTITY | 9 |
| INTERFACE | 9 |
| LIGHT_RING | 8 |
| COORDINATE_CONVENTION | 8 |
| HOLONOMY | 7 |
| MEMBER_IDENTITY | 5 |
| ANGULAR | 5 |
| PERTURBATION_STRUCTURE | 4 |
| VECTOR | 4 |
| BACKGROUND_EOM | 3 |
| REGRESSION_BUG | 2 |
| RADIAL | 2 |
| GLOBAL | 1 |
| HISTORICAL_ONLY | 1 |

## Required gates WITHOUT dedicated tests (coverage-gap backlog)

- G02
- G03
- G04
- G05
- G16
- G50
- G70
- G71
- G72
- G90

## Per-test inventory

### tests/negative/test_production_guards.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_blacklisted_actual_archive_rejected[data/diagnostic/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv]` | PROVENANCE | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | G20 | Production results must derive only from the current single-action emission chain; the stitched SELECTED_41STREAM_V2 historical table must never re-enter production. |
| `test_blacklisted_actual_archive_rejected[archive/full_working_snapshot/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv]` | PROVENANCE | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | G20 | Production results must derive only from the current single-action emission chain; the stitched SELECTED_41STREAM_V2 historical table must never re-enter production. |
| `test_forged_boolean_certificate_rejected` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | A bare boolean pass claim cannot unlock the QNM stage; enabling QNM requires a full hash-bound certificate with artifacts, multipole coverage and gates. |
| `test_corrupt_hash_rejected` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G00 | Every published number is bound to the exact bytes of its input artifact; any modification is detectable. |
| `test_v12_sign` | MEMBER_IDENTITY | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G30 | The selected SVT branch fixes v12 = -v6/(2h); the opposite sign is the blacklisted branch with different vector-sector dynamics. |
| `test_a5_requires_second_derivative` | SYMBOLIC_IDENTITY | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G07 | The a5 identity must keep its genuine -a1'' term; dropping it is the explicitly blacklisted historical bug. |
| `test_jet_refuses_invalid_input[x0-y0]` | HOLONOMY | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G07 | The locked JET9D8 derivative service must abort on invalid input (non-monotone grids, NaN, insufficient stencil) rather than silently produce wrong jets. |
| `test_jet_refuses_invalid_input[x1-y1]` | HOLONOMY | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G07 | The locked JET9D8 derivative service must abort on invalid input (non-monotone grids, NaN, insufficient stencil) rather than silently produce wrong jets. |
| `test_jet_refuses_invalid_input[x2-y2]` | HOLONOMY | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G07 | The locked JET9D8 derivative service must abort on invalid input (non-monotone grids, NaN, insufficient stencil) rather than silently produce wrong jets. |
| `test_schema_count` | PERTURBATION_STRUCTURE | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G30 | The even-parity second-order action has exactly 41 distinct coefficient slots; duplication or loss silently changes the theory being closed. |
### tests/negative/test_mh_action_domain.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_luminal_emitter_rejects_unsupported_action_jet[G4X]` | MH_EQ85 | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G13 | The MH luminal emitter is valid only for G4=G4(phi), G4X=0, G5=0; non-luminal Horndeski jets must abort, not be silently dropped. |
| `test_luminal_emitter_rejects_unsupported_action_jet[G4XX]` | MH_EQ85 | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G13 | The MH luminal emitter is valid only for G4=G4(phi), G4X=0, G5=0; non-luminal Horndeski jets must abort, not be silently dropped. |
| `test_luminal_emitter_rejects_unsupported_action_jet[G4phiX]` | MH_EQ85 | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G13 | The MH luminal emitter is valid only for G4=G4(phi), G4X=0, G5=0; non-luminal Horndeski jets must abort, not be silently dropped. |
| `test_luminal_emitter_rejects_unsupported_action_jet[G5X]` | MH_EQ85 | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G13 | The MH luminal emitter is valid only for G4=G4(phi), G4X=0, G5=0; non-luminal Horndeski jets must abort, not be silently dropped. |
| `test_luminal_emitter_rejects_unsupported_action_jet[G5phi]` | MH_EQ85 | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G13 | The MH luminal emitter is valid only for G4=G4(phi), G4X=0, G5=0; non-luminal Horndeski jets must abort, not be silently dropped. |
| `test_zk_emitter_rejects_horndeski_jets_instead_of_discarding_them[G3X]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G30 | The ZK Appendix-A emitter is derived for the pure U(1)-SVT Einstein baseline; Horndeski jets must abort (they collide with same-named SVT couplings). |
| `test_zk_emitter_rejects_horndeski_jets_instead_of_discarding_them[G4X]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G30 | The ZK Appendix-A emitter is derived for the pure U(1)-SVT Einstein baseline; Horndeski jets must abort (they collide with same-named SVT couplings). |
| `test_zk_emitter_rejects_horndeski_jets_instead_of_discarding_them[G4XX]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G30 | The ZK Appendix-A emitter is derived for the pure U(1)-SVT Einstein baseline; Horndeski jets must abort (they collide with same-named SVT couplings). |
| `test_zk_emitter_rejects_horndeski_jets_instead_of_discarding_them[G5X]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G30 | The ZK Appendix-A emitter is derived for the pure U(1)-SVT Einstein baseline; Horndeski jets must abort (they collide with same-named SVT couplings). |
| `test_zk_emitter_rejects_horndeski_jets_instead_of_discarding_them[G5phi]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G30 | The ZK Appendix-A emitter is derived for the pure U(1)-SVT Einstein baseline; Horndeski jets must abort (they collide with same-named SVT couplings). |
### tests/negative/test_absolute_closure.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_strict_pipeline_requires_direct_krgm` | GLOBAL | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | The absolute-closure verdict is fail-closed: without direct certificate/products/spectral evidence the pipeline must report failure. |
| `test_absolute_full_closure_certificate_cannot_replace_direct_products` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | A hand-written ABSOLUTE_FULL_CLOSURE.json can never substitute recomputed direct matrix products. |
| `test_qnm_same_operator_and_convergence_gate[operator]` | LIGHT_RING | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | QNM evidence is admissible only with the same per-L matrices as the certified products, converged on all axes, and with an Im(omega)>0 instability search. |
| `test_qnm_same_operator_and_convergence_gate[axis]` | LIGHT_RING | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | QNM evidence is admissible only with the same per-L matrices as the certified products, converged on all axes, and with an Im(omega)>0 instability search. |
| `test_qnm_same_operator_and_convergence_gate[one_setting]` | LIGHT_RING | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | QNM evidence is admissible only with the same per-L matrices as the certified products, converged on all axes, and with an Im(omega)>0 instability search. |
| `test_qnm_same_operator_and_convergence_gate[error]` | LIGHT_RING | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | QNM evidence is admissible only with the same per-L matrices as the certified products, converged on all axes, and with an Im(omega)>0 instability search. |
| `test_qnm_same_operator_and_convergence_gate[no_upper_search]` | LIGHT_RING | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | QNM evidence is admissible only with the same per-L matrices as the certified products, converged on all axes, and with an Im(omega)>0 instability search. |
### tests/negative/test_matrix_guards.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_reject_bad_operator[nan]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G32 | The constraint-reduced operator must be well-posed: finite, K/G/M symmetric, S antisymmetric, R block vanishing, constraint pivots nonzero. |
| `test_reject_bad_operator[asymmetric]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G32 | The constraint-reduced operator must be well-posed: finite, K/G/M symmetric, S antisymmetric, R block vanishing, constraint pivots nonzero. |
| `test_reject_bad_operator[symmetric_S]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G32 | The constraint-reduced operator must be well-posed: finite, K/G/M symmetric, S antisymmetric, R block vanishing, constraint pivots nonzero. |
| `test_reject_bad_operator[zero_pivot]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G32 | The constraint-reduced operator must be well-posed: finite, K/G/M symmetric, S antisymmetric, R block vanishing, constraint pivots nonzero. |
| `test_reject_bad_operator[nonzero_R]` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G32 | The constraint-reduced operator must be well-posed: finite, K/G/M symmetric, S antisymmetric, R block vanishing, constraint pivots nonzero. |
| `test_symmetric_whitening_and_negative_kinetic` | KINETIC | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G40 | After symmetric whitening the radial characteristic matrix is positive, and a nonpositive kinetic eigenvalue (ghost) must abort the stability gate. |
### tests/negative/test_certificate_provenance.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_schema_and_hash_guard_accepts_complete_synthetic_fixture` | PROVENANCE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G100 | The QNM certificate guard must be non-vacuous: a genuinely complete, hash-bound certificate is accepted (fails closed without rejecting everything). |
| `test_schema_valid_forgery_is_rejected[artifact_hash]` | PROVENANCE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G100 | Schema conformance alone never enables QNM: hash binding to the exact action member, all gates, full multipole coverage and in-repo artifact paths are required. |
| `test_schema_valid_forgery_is_rejected[action_hash]` | PROVENANCE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G100 | Schema conformance alone never enables QNM: hash binding to the exact action member, all gates, full multipole coverage and in-repo artifact paths are required. |
| `test_schema_valid_forgery_is_rejected[gate]` | PROVENANCE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G100 | Schema conformance alone never enables QNM: hash binding to the exact action member, all gates, full multipole coverage and in-repo artifact paths are required. |
| `test_schema_valid_forgery_is_rejected[multipole]` | PROVENANCE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G100 | Schema conformance alone never enables QNM: hash binding to the exact action member, all gates, full multipole coverage and in-repo artifact paths are required. |
| `test_schema_valid_forgery_is_rejected[escape]` | PROVENANCE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G100 | Schema conformance alone never enables QNM: hash binding to the exact action member, all gates, full multipole coverage and in-repo artifact paths are required. |
### tests/negative/test_blacklist.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_blacklist_is_nonempty` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G20 | The production blacklist must stay populated — it keeps forbidden historical artifacts out of the production input boundary. |
| `test_license_and_authors_metadata` | SOFTWARE_INFRASTRUCTURE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | GENERAL | Release metadata must keep licensing and author attribution consistent (ACSL v1.4, Carmen Casu and Lino Casu). |
### tests/negative/test_qnm_guard.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_qnm_blocked_without_certificate` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G100 | The QNM stage is closed by default: without a direct-global-KRGM certificate any QNM attempt fails loudly. |
### tests/regression/test_light_ring_identity.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_structural_invariants` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G10 | In E11-E00 of the unsplit core every genuine-SVT coupling (f3/f3X) is electric-activated (carries A0'^2), and exactly one f3-free term survives — the MH f2X channel. |
| `test_P_A_analytic_nonzero` | VECTOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | GENERAL | JA depends linearly on f2F with analytically nonzero coefficient P_A = sqrt(h/f)*A0'*r^2, so the electric slot is identifiable in the vector equation. |
| `test_epsY_nulltest` | SVT_OPERATOR | TAUTOLOGICAL | BRANCH_AGNOSTIC | G11 | eps_Y contributions vanish at A0'=0 with f2Y symbolic — validates the nulltest helper on synthetic terms. |
| `test_P_MH_leaves_f2X_channel` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G10 | P_MH of E11-E00 retains exactly the f2X MH piece with no A0' dependence — the projection isolates a purely metric subchannel. |
| `test_epsY_slot_shift_structural` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G11 | Slot shifts Delta f2 = eps_Y*Y (Y ~ A0'^2) propagate the A0'^2 factor into all EOM contributions — background-null at A0'=0 without setting f2Y=0. |
| `test_epsY_nulltest_uses_f2Y_symbolic` | SVT_OPERATOR | TAUTOLOGICAL | BRANCH_AGNOSTIC | G11 | The nulltest helper accepts f2Y as an unevaluated symbol — the coupling is never silently zeroed. |
| `test_C85_undivided_structure` | MH_EQ85 | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G12 | C85_undivided contains the geometric a4 term, F_MH channel and electric term, and does not collapse at A0'=0 (anti-merge rule). |
| `test_epsY_slot_shift_contributions_vanish` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G11 | Module-level Test A: full eps_Y slot-shift EOM contributions vanish identically at A0'=0. |
| `test_f2Y_explicit_EOM_core` | SYMBOLIC_IDENTITY | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G10 | The f2Y-explicit core carries f2F - 2h phi'^2 f2Y visibly; f2Y=0 is an exact substitution, never a silent omission. |
| `test_C85_undivided_executable` | MH_EQ85 | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G12 | The executable undivided Eq. 85 equals an independently transcribed copy of arXiv:2301.10362 Eq. 85 term by term, including the exact electric coefficient +2rfh A0'^2 on v8_MH. |
| `test_P_MH_preserves_A0prime_symbolic` | MH_EQ85 | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G12 | Anti-merge rule: the MH projection kills f3 terms but never touches A0' — the electric term survives projection untouched. |
| `test_G10_C_bg_general_symbolic` | MH_EQ85 | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G10 | C_bg carries every jet, A0', and Delta22_SVT symbolically with machine-verified ALPHA/BETA — no silent zeroing at assembly. |
| `test_G12_MH_projection_exact_eq85` | MH_EQ85 | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G12 | Central identity (Test B): P_MH[C_bg] with canonical slot identifications reproduces the undivided KT2023 Eq. 85 exactly, A0' and f2Y symbolic. |
| `test_G11_epsY_background_null_on_C_bg` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G11 | Test A on the same C_bg: at A0'=0 the f2Y coefficient vanishes structurally while f2Y stays visible in the general object. |
| `test_G12_negative_control_wrong_slot` | NEGATIVE_CONTROL | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G12 | Negative control: a wrong v8 slot (scaled 11/10, the historical V8-vs-V9 confusion class) must break the exact Eq. 85 identity. |
| `test_G15_sigma_svt_operator_decomposition` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G15 | Sigma_SVT = C_bg - P_MH[C_bg] decomposes exactly into the registered operator list (E00-channel f3 and f4/f4X/tf4 terms with ap^2, plus Delta22_SVT at BETA weight). |
| `test_G14_symbol_separation_u_lr_vs_A0_u` | COORDINATE_CONVENTION | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G14 | Symbol-collision guard: u_lr (radial coordinate at the ring) and A0_u (electric derivative) are distinct; the geometric term uses u_lr, the electric term A0_u. |
| `test_ward_residual_electric_channel_JA` | VECTOR | TAUTOLOGICAL | BRANCH_AGNOSTIC | G15 | Registered Ward-residual split: the unmatched electric part is exactly c4*JA (non-derivative current), R_met carries no electric structures — pinning the split the covariant derivation must reproduce. |
### tests/regression/test_kt_mh_background.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_s0_slice_diagonal_map` | SYMBOLIC_IDENTITY | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G10 | Convention bridge: on the S0 slice the KT2023 components and HT2018 variations are related by the exact diagonal map KT_E00 = -HT_E00/(r^2 f), KT_E11 = +HT_E11/(r^2 f). |
| `test_eq85_generating_two_term_identity` | MH_EQ85 | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G10 | Eq.-85-generating recipe: in the luminal branch C85_undivided == A*KT_E00 + B*KT_E22 exactly off-shell with A = -r f^(3/2)/sqrt(h), B = +r f^(3/2)/sqrt(h). |
| `test_two_term_identity_sector_scope` | NEGATIVE_CONTROL | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G10 | Scope evidence: reopening G4X breaks the metric-only two-term identity — Eq. 85 is on-shell outside the luminal lock. |
| `test_eq85_onshell_certificate_luminal_branch` | MH_EQ85 | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G10 | Primary-source certificate: in the luminal branch the undivided Eq. 85 vanishes identically when phi''/h' are solved from E00 = E22 = 0; the pivot carries all premises. |
| `test_a4_v8_slot_definitions` | PROVENANCE | TAUTOLOGICAL | BRANCH_AGNOSTIC | G12 | Canonical slots (KT2023 App. A Eq. 168): a4 = sqrt(fh)/2 * Hcal, v8 = G2F/(2 sqrt(fh)) — MH v8 = ZK V9, never the V8 column. |
| `test_ward_span_derivation_status` | PROVENANCE | TAUTOLOGICAL | BRANCH_AGNOSTIC | GENERAL | Honest derivation record: the Ward-span no-go was a non-holonomic artifact; the holonomic span closes all channels with metric-only coefficients and one open residual. |
### tests/regression/test_full_closure_attempt_20260919.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_central_single_action_checkpoint_is_positive_for_all_required_L` | KINETIC | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | GENERAL | Central single-action checkpoint: kinetic matrix positive and radial sector stable at every required L. |
| `test_handover_checkpoint_rejects_f4_only_and_retains_rank_six` | PERTURBATION_STRUCTURE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | GENERAL | A single f4 control cannot solve the joint background; the six-primitive control space retains rank 6; certifications stay PENDING. |
| `test_cinfinity_flat_step_endpoint_and_monotonic_contract` | SOFTWARE_INFRASTRUCTURE | INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | The C-infinity blending step is exactly flat at/beyond endpoints, monotone in between, midpoint 1/2 — seam gluing injects no spurious derivatives. |
### tests/regression/test_eq85_naming_guard.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_v_slot_mapping_is_documented` | PROVENANCE | TAUTOLOGICAL | BRANCH_AGNOSTIC | G01 | Slot-mapping documentation guard: the kernel documents V9=MHv8 (never the 13-slot V8) together with the neighbor mappings. |
| `test_q_conventions_must_not_be_mixed` | COORDINATE_CONVENTION | TAUTOLOGICAL | BRANCH_AGNOSTIC | G01 | q_r = A0_r^2 v8 and q_u = A0_u^2 v8 differ by the Jacobian factor u^4/r_s^2 — raw comparison is a coordinate-mixing bug; the frozen u^4 constant exposes drift. |
### tests/regression/test_p5_light_ring_zero_vector.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_inner_light_ring_excludes_positive_zero_vector_tensor_branch` | LIGHT_RING | PARTIALLY_INDEPENDENT | ZERO_VECTOR_HISTORICAL | G14 | On the frozen zero-vector carrier, Eq. 85 resolves two rings; at the inner ring the required tensor slot is negative (f/h < -1.6) while the carrier has a4 > 0, H > 0 — the positive zero-vector branch cannot traverse the inner ring; the outer ring requires f/h ~ +1. |
### tests/regression/test_golden_jets.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_shared_jet_matches_frozen_scalar_implementation` | SOFTWARE_INFRASTRUCTURE | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | The shared JET9D8 service reproduces the frozen archived implementation exactly — all derivative-based numerical identities remain the certified ones. |
### tests/regression/test_frozen_onshell_diagnostic.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_frozen_carrier_violation_is_resolved_above_derivative_sensitivity` | REGRESSION_BUG | PARTIALLY_INDEPENDENT | ZERO_VECTOR_HISTORICAL | GENERAL | Locks a KNOWN unresolved violation: the frozen zero-vector carrier fails the Eq. 85 on-shell identity with residual ~ 0.071, far above derivative noise — the contradiction must stay visibly open. |
### tests/regression/test_electric_hybrid_principal_checkpoint.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_principal_search_checkpoint_is_improved_but_not_promoted` | KINETIC | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | GENERAL | Electric hybrid kinetic feasibility: smallest kinetic eigenvalue improves with L (small negatives at L=6/12/20, strictly positive at L >= 42), far from the spurious Newton branch — improved, not solved. |
### tests/regression/test_electric_hybrid_background_eom.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_electric_hybrid_background_eom_checkpoint` | BACKGROUND_EOM | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | GENERAL | The electric hybrid central background satisfies E00 = E11 = JA = 0 within policy tolerance on the production window; the scalar E_phi residual sits at the known direct-differentiation floor (archived exact reference, tolerance not relaxed). |
### tests/test_release_auditor.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_constructive_closure_passes` | SOFTWARE_INFRASTRUCTURE | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | GENERAL | Release meta-guard: the constructive-closure auditor certifies FULL_CONSTRUCTIVE_CLOSURE = PASS over the frozen artifacts. |
### tests/integration/test_constructive_audit.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_constructive_audit_passes` | SOFTWARE_INFRASTRUCTURE | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | GENERAL | Integration duplicate of the release meta-guard with weaker substring assertion. |
### tests/unit/test_regional_production.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_selected_member_is_regional_and_allows_electric_central_background` | MEMBER_IDENTITY | INDEPENDENT | CURRENT_ELECTRIC | G20 | Production is bound to the six-region regional member; A0' != 0 only in genuine-SVT/handover regions, never in pure-Horndeski regions. |
| `test_invalid_coordinates_cannot_silently_select_core[nan]` | COORDINATE_CONVENTION | INDEPENDENT | BRANCH_AGNOSTIC | G20 | NaN/inf/negative radii are rejected, not silently mapped into a region branch. |
| `test_invalid_coordinates_cannot_silently_select_core[inf]` | COORDINATE_CONVENTION | INDEPENDENT | BRANCH_AGNOSTIC | G20 | NaN/inf/negative radii are rejected, not silently mapped into a region branch. |
| `test_invalid_coordinates_cannot_silently_select_core[-1.0]` | COORDINATE_CONVENTION | INDEPENDENT | BRANCH_AGNOSTIC | G20 | NaN/inf/negative radii are rejected, not silently mapped into a region branch. |
| `test_wrong_region_cannot_pass_member_guard` | MEMBER_IDENTITY | INDEPENDENT | CURRENT_ELECTRIC | G20 | Stream region labels must agree with the locked radial cover; self-declared labels cannot override the resolver. |
| `test_changed_member_cover_rejected` | PROVENANCE | INDEPENDENT | CURRENT_ELECTRIC | G20 | The member file's region cover is locked to the resolver; boundary edits invalidate production provenance. |
| `test_shared_baseline_is_subtracted_before_selection` | SVT_OPERATOR | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G30 | Unreduced H+SVT assembly = Horndeski + SVT - shared (inclusion-exclusion), then locked recanonicalization v7 = v2^2/(4v1), v12 = -v6/(2h), pivot v1 nonzero. |
| `test_central_independent_schur_matches_euler_operator` | KINETIC | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G32 | Two independent central reductions (Euler-operator reducer vs Schur elimination) agree on the kinetic matrix, and the strong-field kinetic instability (A0' > 1, negative eigenvalue) is faithfully reproduced, not hidden. |
### tests/unit/test_unreduced_descriptor.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_descriptor_builds_full_eight_field_euler_operator` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G30 | The even-parity EOM assembles as an 8-field descriptor (DAE) Euler operator with correct derivative-order structure and finite principal symbol. |
| `test_local_schur_crosscheck_reports_conditioning` | SVT_OPERATOR | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G32 | Local principal-symbol Schur complement is the sanctioned crosscheck; auxiliary-block conditioning is reported, never hidden. |
| `test_constant_profile_matches_existing_unreduced_symbol` | SVT_OPERATOR | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G30 | The new profile-aware descriptor and the golden unreduced-even kernel have identical Euler-Lagrange Fourier symbols. |
| `test_generalized_psi_descriptor_keeps_h0_constraint_uneliminated` | PERTURBATION_STRUCTURE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G31 | Under the generalized-psi map the H0 Euler row remains a proper constraint-multiplier row (no forbidden derivatives) — avoids the ill-conditioned D_h1 inversion. |
| `test_h0_closed_formula_matches_descriptor_with_same_jet_service` | PERTURBATION_STRUCTURE | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G31 | The descriptor's H0 constraint row equals the closed-form generalized-psi coefficients; the H0^2 term cancels through v7 = v2^2/(4v1). |
| `test_descriptor_pullback_matches_established_reducer_on_real_carrier` | SVT_OPERATOR | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G32 | The descriptor pulled back through the accepted constraint maps reproduces the certified 3-field production reducer at operator level on the real carrier. |
### tests/unit/test_onshell_identity.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_exact_vacuum_controls[0.0-0.0]` | MH_EQ85 | INDEPENDENT | BRANCH_AGNOSTIC | G12 | For exact vacuum solutions (Schwarzschild/Kottler) the KT2023 Eq. 85 identity holds and the independent E00-E22 subtraction vanishes. |
| `test_exact_vacuum_controls[0.4-0.0]` | MH_EQ85 | INDEPENDENT | BRANCH_AGNOSTIC | G12 | For exact vacuum solutions (Schwarzschild/Kottler) the KT2023 Eq. 85 identity holds and the independent E00-E22 subtraction vanishes. |
| `test_exact_vacuum_controls[0.0-0.01]` | MH_EQ85 | INDEPENDENT | BRANCH_AGNOSTIC | G12 | For exact vacuum solutions (Schwarzschild/Kottler) the KT2023 Eq. 85 identity holds and the independent E00-E22 subtraction vanishes. |
| `test_reissner_nordstrom_requires_electric_term` | MH_EQ85 | INDEPENDENT | CURRENT_ELECTRIC | G12 | Reissner-Nordstrom satisfies Eq. 85 only WITH the Maxwell electric term; omitting it leaves the exact analytic residual -f Q^2/r^3. |
| `test_independent_background_subtraction_with_varying_G4` | MH_EQ85 | INDEPENDENT | BRANCH_AGNOSTIC | G12 | Eq. 85 is algebraically equivalent to E00 = E22 for the luminal G5=0 MH sector with varying G4, off-shell, all G2/G3 freedom cancelling. |
| `test_no_division_at_light_ring` | LIGHT_RING | INDEPENDENT | BRANCH_AGNOSTIC | G14 | The undivided Eq. 85 implementation stays exact and finite at 2f - r f' = 0 — never divides by the light-ring denominator. |
### tests/unit/test_central_action.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_background_equations_reissner_nordstrom_control` | BACKGROUND_EOM | INDEPENDENT | CURRENT_ELECTRIC | G06 | The electric central background equations (E00, E11, J_A, undivided Eq. 131) are exactly satisfied by a Reissner-Nordstrom-type solution with conserved charge 0.4. |
| `test_background_checker_detects_action_change` | NEGATIVE_CONTROL | INDEPENDENT | CURRENT_ELECTRIC | G06 | Calibration: perturbing f2 by +0.1 shifts E00 by exactly the analytic amount — the residual validator is falsifiable, not trivially zero. |
| `test_zero_electric_identity_remains_finite` | BACKGROUND_EOM | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | G06 | The undivided Eq. 131 form remains finite at A0' = 0 — never divides by A'. |
| `test_central_published_formula_matches_independent_schur` | KINETIC | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G40 | The published ZK kinetic formula (valid only on the Einstein+genuine-SVT branch) and an independent Schur elimination agree at every tested L; off-branch use fails loudly. |
| `test_central_action_consumes_explicit_inputs` | SVT_OPERATOR | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G30 | The central action replay consumes explicit action inputs: v10 responds to f2F with the exact Appendix-A coefficient; d3 uses the mixed action jet, not the total radial derivative. |
### tests/unit/test_transition_projected_controllability.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_projected_joint_action_tangent_is_locally_controllable` | KINETIC | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G40 | In the transition window the projected joint H+SVT tangent contains a direction simultaneously improving the kinetic minimum eigenvalue for all required L — local controllability of the kinetic obstruction. |
| `test_projected_audit_never_uses_historical_inner_as_absolute_seed` | PROVENANCE | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | G20 | Action-first contract enforced: no certifiable absolute action seed exists right of u = 0.71; the historical Inner stream is never used as an absolute seed. |
| `test_horndeski_a1_survives_onshell_projection_but_is_not_promoted` | KINETIC | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G40 | The Horndeski a1 kinetic-gain direction survives projection onto the background-null space (>50% norm at the inner ring) but is explicitly NOT promoted to a finite transition solution. |
| `test_projection_numerics_are_converged_and_a1_width_warning_is_preserved` | NEGATIVE_CONTROL | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G40 | Tangent numerics converged; the honest warning preserved: a1's all-L gain sign is not robust to sampling width — a1 alone must not be promoted. |
### tests/unit/test_strong_field_transition_contract.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_physical_interval_contains_light_ring_but_not_a_physical_seam` | COORDINATE_CONVENTION | INDEPENDENT | BRANCH_AGNOSTIC | G80 | The strong-field transition interval contains the stable inner light ring; the historical seam u = 0.71 is inside but physically distinct from the ring. |
| `test_transition_coordinate_and_basis_are_flat_at_physical_endpoints` | COORDINATE_CONVENTION | INDEPENDENT | BRANCH_AGNOSTIC | G80 | Transition coordinate maps to s in [0,1]; the C-infinity bump and basis columns vanish exactly at both endpoints — no boundary artifacts injected. |
| `test_mode_tracking_preserves_eigenvector_identity_through_order_exchange` | KINETIC | INDEPENDENT | BRANCH_AGNOSTIC | G40 | At kinetic eigenvalue branch crossings, mode identity is tracked by eigenvector overlap, not sorted eigenvalue labels — the tracked branch keeps its physical eigenvector. |
| `test_first_zero_crossing_is_interpolated` | SOFTWARE_INFRASTRUCTURE | INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | Instability onset is estimated by linear interpolation inside the bracketing interval, not snapped to a grid point. |
### tests/unit/test_repo_observability.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_repo_snapshot_surfaces_current_truth` | SOFTWARE_INFRASTRUCTURE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | GENERAL | The status view reports the honest state: closure NOT_CERTIFIED, member search open, QNM gate BLOCKED — no overstating. |
| `test_evidence_index_hashes_machine_reports` | PROVENANCE | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | G00 | Every machine report is sha256-hashed in the evidence index — tamper-evident, rejection status surfaced. |
| `test_gate_and_member_views_keep_rejections_separate_from_active_search` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G20 | Rejected members stay distinct from the active search; the Horndeski principal-feasibility gate is visible. |
| `test_registry_declares_expected_rejection_exit_codes` | PROVENANCE | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | G00 | Intentional rejection witnesses (expected exit 2) are distinguished from current PASS audits (exit 0) — expected-failures are not broken pipelines. |
### tests/unit/test_import_provenance.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_import_provenance` | PROVENANCE | INDEPENDENT | BRANCH_AGNOSTIC | G00 | All SSZ modules import from the canonical checkout, never from the historical uppercase snapshot, inside the running tree. |
| `test_no_forbidden_path_in_sys_path` | PROVENANCE | INDEPENDENT | BRANCH_AGNOSTIC | G00 | The historical snapshot directory must not be on sys.path — prevents mixed-provenance imports. |
| `test_gate_dependency_enforcement` | SOFTWARE_INFRASTRUCTURE | INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | Fail-closed gating: PASS with a non-PASS required parent certifies BLOCKED_BY_DEPENDENCY transitively; closure verdict derives from certification. |
| `test_release_metadata_current` | PROVENANCE | INDEPENDENT | BRANCH_AGNOSTIC | G01 | MODEL_LOCK/EVIDENCE_INDEX git stamps equal HEAD (or an ancestor with clean tree) — packaged evidence is bound to the code that produced it. |
### tests/unit/test_holonomic_hessian.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_transverse_block_implies_exact_background_null_hessian` | HOLONOMY | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G07 | A background-null f2 deformation has only 3 independent Hessian DOF: nullity forces the mixed/phi Hessian jets from the transverse block (holonomy). |
| `test_first_chain_rule_cannot_be_repaired_by_principal_block` | HOLONOMY | INDEPENDENT | BRANCH_AGNOSTIC | G07 | The first chain-rule row depends only on lower Hessian controls; the principal block cannot repair it — split-control designs are structurally unrepairable. |
| `test_audit_flags_nonholonomic_split_controls` | NEGATIVE_CONTROL | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G07 | The audit flags FAIL for Hessian splits off the holonomic background-null section — INNER_DIRECT_41 is not certifiable from them. |
| `test_appendix_a_lower_response_contains_cross_and_radial_terms` | SVT_OPERATOR | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G30 | The correct ZK Appendix-A f2-Hessian response includes the A0'^2 cross term in c3 and the radial-derivative term in e3 — the historical split-map defect. |
### tests/unit/test_eq47_projection.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_projected_relation_and_literal_residual_are_distinct` | ANGULAR | INDEPENDENT | BRANCH_AGNOSTIC | G60 | The published Eq. 47 shortcut is valid only after the K-to-M projection factor; the literal shortcut differs by a real, tracked R-dependent residual. |
| `test_leading_order_is_scale_invariant_and_not_hardcoded` | ANGULAR | INDEPENDENT | BRANCH_AGNOSTIC | G60 | The leading eps-power (propagating-mode count) is inferred from coefficients, scale-invariant, never hardcoded. |
| `test_reference_action_and_projection_checkpoint` | ANGULAR | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G60 | Central electric member at u~0.65: action-derived angular characteristic satisfies the projected Eq. 46/47 relations; roots match the frozen reference; closure explicitly false. |
| `test_large_root_preserves_cubic_degree` | REGRESSION_BUG | INDEPENDENT | BRANCH_AGNOSTIC | G60 | Root extraction preserves small leading coefficients (huge physical roots ~1e16) — no silent mode loss from polynomial hygiene. |
### tests/unit/test_zk_action_selector.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_default_selector_is_exactly_legacy_holonomic_path` | COORDINATE_CONVENTION | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G07 | The ZK emitter's d3 phi-derivative selector defaults to the holonomic on-curve path (regresses the archived 37/41 table); no silent convention change. |
| `test_action_selector_uses_mixed_action_jets_in_d3_only` | SYMBOLIC_IDENTITY | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G07 | The 'action' selector substitutes the mixed action jet only in d3 (0.5*A0'*dphi(v6)); no other slot depends on the selector. |
| `test_action_selector_rejects_unknown_mode` | NEGATIVE_CONTROL | SELF_CONSISTENCY_ONLY | BRANCH_AGNOSTIC | GENERAL | Unknown selector values fail loudly — no silent defaulting between derivative conventions. |
### tests/unit/test_principal_controls.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_raw_response_matches_actual_emitter` | SYMBOLIC_IDENTITY | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G07 | The closed-form f2-Hessian control response (v1, v4, c2) is exactly what the Appendix-A emitter implements — inversion targets real behavior. |
| `test_raw_inverse_roundtrip_includes_zero_electric_endpoint` | INTERFACE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G80 | The triangular raw response is invertible where its diagonal is nonzero; at A0'=0 rows the vector response vanishes and must be handled by the reduced block, not filled. |
| `test_unreachable_zero_field_vector_target_is_rejected` | NEGATIVE_CONTROL | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G80 | At A0'=0 no f2-Hessian deformation can move v1; demanding one is physically unreachable and rejected. |
### tests/unit/test_holonomic_hessian_y.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_4d_completion_is_background_null` | HOLONOMY | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G07 | A background-null f2(...,Y) deformation has a symmetric 4x4 Hessian obeying the chain rule H@(phi',X',F',Y')=0; the 4D completion from six transverse jets is the unique such Hessian. |
| `test_local_algebraic_matrix_matches_direct_response` | SYMBOLIC_IDENTITY | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G07 | The Appendix-A response in five algebraic channels collapses to a local matrix M(q); closed-form matrix and direct Hessian evaluation are the same physics. |
| `test_y_sector_changes_response_without_breaking_holonomy` | HOLONOMY | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G07 | Y-sector Hessians are genuine electric controls (move v1/v4/c2) while the completed 4D Hessian stays background-null — reachability without breaking holonomy. |
### tests/unit/test_descriptor_pencil.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_local_poly_matrix_differentiates_degree_eight_polynomial` | RADIAL | INDEPENDENT | BRANCH_AGNOSTIC | G100 | The QNM descriptor differentiation matrix is polynomially exact to degree 8 — the linear counterpart of JET9D8. |
| `test_quadratic_linearization_preserves_singular_mass` | KINETIC | INDEPENDENT | BRANCH_AGNOSTIC | G100 | The descriptor linearization keeps the singular kinetic matrix genuinely singular — constraints give infinite eigenvalues, not a fake finite spectrum. |
| `test_semidiscrete_descriptor_keeps_dae_singular_kinetic_rows` | KINETIC | INDEPENDENT | BRANCH_AGNOSTIC | G100 | Semi-discretization preserves node-wise kinetic singularity — the DAE structure survives discretization, no premature Schur complement. |
### tests/unit/test_constraint_stationarity.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_h1_and_vector_constraints_stationary[6]` | SYMBOLIC_IDENTITY | INDEPENDENT | BRANCH_AGNOSTIC | G30 | The h1 and Maxwell-vector constraint maps are stationary points of the unreduced quadratic action: their Euler equations hold identically at L = 6, 42, 1000. |
| `test_h1_and_vector_constraints_stationary[42]` | SYMBOLIC_IDENTITY | INDEPENDENT | BRANCH_AGNOSTIC | G30 | The h1 and Maxwell-vector constraint maps are stationary points of the unreduced quadratic action: their Euler equations hold identically at L = 6, 42, 1000. |
| `test_h1_and_vector_constraints_stationary[1000]` | SYMBOLIC_IDENTITY | INDEPENDENT | BRANCH_AGNOSTIC | G30 | The h1 and Maxwell-vector constraint maps are stationary points of the unreduced quadratic action: their Euler equations hold identically at L = 6, 42, 1000. |
### tests/unit/test_production_regions.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_region_boundaries_exclude_strong_h_shortcut` | MEMBER_IDENTITY | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G20 | The u-partition is locked so the strong-H carrier can never shortcut the inner light-ring band (four regions, half-open boundaries). |
| `test_both_disputed_points_are_central_svt` | LIGHT_RING | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G20 | Both disputed light-ring radii lie in the central exact-SVT lobe with A0' > 1.4 — the ring controversy is resolved inside the electric exact-SVT member. |
### tests/unit/test_lower_order_controls.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_lower_order_inverse_replays_targets_and_reduced_endpoint` | INTERFACE | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G80 | The inner lower-order response is exactly diagonal; inversion realizes targets as action controls; at A0'=0 only the reduced (c3, e3) block acts. |
| `test_nonzero_v5_at_vector_free_endpoint_is_rejected` | NEGATIVE_CONTROL | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G80 | A v5 target at an A0'=0 row is physically unreachable (response ∝ A0') and must fail loudly, not be epsilon-filled. |
### tests/unit/test_jet9d8.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_polynomial_exact_interior` | RADIAL | INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | JET9D8 is exact on polynomial profiles in the interior — every radial jet consumed downstream is trustworthy. |
| `test_profile_service_supports_declared_convergence_stencils` | SOFTWARE_INFRASTRUCTURE | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | The cached vectorized profile-jet service reproduces the reference derivative for all declared convergence stencils (7/6, 9/8, 11/8, 13/8) through order 2. |
### tests/unit/test_inner_targets.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_central_endpoint_jet_uses_only_central_side` | INTERFACE | INDEPENDENT | BRANCH_AGNOSTIC | G80 | Inner handover targets are seeded by one-sided jets from the central side only — no information leaks across the interface. |
| `test_target_preserves_flat_partition_endpoint_jet` | INTERFACE | INDEPENDENT | BRANCH_AGNOSTIC | G80 | The S_SVT/T_H partition multiplies the target linear extension: central jet inside the SVT window, exactly zero in the flat core-side window. |
### tests/unit/test_inner_export.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_endpoint_jet_matches_polynomial` | INTERFACE | INDEPENDENT | BRANCH_AGNOSTIC | G80 | The interface report's endpoint jets measure value/first/second derivatives exactly on polynomials — the seam bookkeeping quantifies what it claims. |
| `test_endpoint_report_detects_lower_order_selection_jump` | INTERFACE | PARTIALLY_INDEPENDENT | BRANCH_AGNOSTIC | G80 | The interface comparison localizes a v5 jump to exactly that slot at radial order 0 with correct magnitude; all other 41 slots zero. |
### tests/unit/test_inner_endpoint_contract.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_endpoint_orders_are_generated_from_contract_only` | INTERFACE | TAUTOLOGICAL | BRANCH_AGNOSTIC | G80 | Endpoint derivative orders encode the radial structure of the lower coefficients (v5/c3: order 1, e3: order 0) — single source of truth. |
| `test_endpoint_equation_count_matches_contract` | INTERFACE | TAUTOLOGICAL | BRANCH_AGNOSTIC | G80 | The frozen contract fixes the endpoint equation count (2 sides x 5 = 10) — no silent enlargement of the constraint system. |
### tests/unit/test_eps_y.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_eps_y_changes_only_v1_v10` | VECTOR | PARTIALLY_INDEPENDENT | ZERO_VECTOR_HISTORICAL | G11 | On the zero-vector branch the eps_Y*Y deformation is exactly background-null and acts only through the vector slots, rescaling v1 and v10 by Z_A = 1 - 2 eps_Y kappa; all other slots bit-identical. |
| `test_eps_y_rejects_electric_background` | NEGATIVE_CONTROL | SELF_CONSISTENCY_ONLY | ZERO_VECTOR_HISTORICAL | G11 | The eps_Y deformation is certified only on the zero-electric branch; applying it to A0' != 0 is refused. |
### tests/unit/test_central_export.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_central_comparison_identifies_only_changed_slot` | MEMBER_IDENTITY | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G20 | The central-stream comparator localizes single-slot drift to exactly that slot and first deviating radius — slot-level member integrity audit. |
| `test_central_comparison_rejects_nonfinite_and_wrong_grid` | NEGATIVE_CONTROL | SELF_CONSISTENCY_ONLY | CURRENT_ELECTRIC | G20 | The comparator is fail-closed: nonfinite coefficients and grid mismatches abort instead of producing silent certificates. |
### tests/unit/test_angular_universal_provenance.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_action_minus_branch_matches_m22_and_cross_identity` | ANGULAR | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G60 | The Laurent-reduced angular mass M22_0 equals the closed-form action shortcut on the action-minus branch; the branch difference obeys the exact cross-identity M22^- - M22^+ = a6 v6 v13/(a4 v9). |
| `test_raw_characteristic_keeps_vector_root_near_luminal` | ANGULAR | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G60 | The raw angular characteristic polynomial retains a vector root at (near) luminal speed — luminality survives the formal reduction without hand repair. |
### tests/unit/test_light_ring_electric.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_light_ring_identity_is_denominator_free_and_electric_control_is_finite` | MH_EQ85 | INDEPENDENT | BRANCH_AGNOSTIC | G14 | The Eq.-85 light-ring balance is evaluated denominator-free (the ring is a regular point) and the required electric control q stays finite at the ring. |
### tests/unit/test_inner_principal.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_principal_targets_match_endpoint_taylor_values` | INTERFACE | INDEPENDENT | BRANCH_AGNOSTIC | G80 | Inner principal (v1, v4, c2) targets are one-sided endpoint Taylor values blended by the stored partition — C-infinity matching, no coefficient interpolation as production. |
### tests/unit/test_hsvt_eps_y_direct.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_direct_hsvt_carrier_is_finite_and_background_null_branch` | VECTOR | PARTIALLY_INDEPENDENT | ZERO_VECTOR_HISTORICAL | G11 | The direct zero-vector HSVT carrier is background-null by construction (A0' = 0, c3 = e3 = v5 = 0) with positive vector kinetic normalization Z_A. |
### tests/unit/test_full_action_lower.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_full_action_lower_direct_v5_and_central_diagnostic` | HISTORICAL_ONLY | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | G07 | Direct one-action lower slots reproduce the historical v5, but the selected c3/e3 are NOT a single-action replay — the audit keeps flagging them as not promotable. |
### tests/unit/test_electric_hybrid_controls.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_electric_hybrid_recipe_replays_corrected_frozen_checkpoint` | PROVENANCE | PARTIALLY_INDEPENDENT | CURRENT_ELECTRIC | GENERAL | The electric-hybrid control recipe replays the frozen corrected checkpoint to 5e-8 — the feasibility study stays auditable (checkpoint is not a production member). |
### tests/unit/test_angular_laurent_series.py

| Test | Category | Independence | Member | Gate | Assertion |
|---|---|---|---|---|---|
| `test_formal_reciprocal` | SOFTWARE_INFRASTRUCTURE | INDEPENDENT | BRANCH_AGNOSTIC | GENERAL | The truncated Laurent-series algebra inverts formally and accurately — the backbone of the large-L angular reducers is arithmetically sound. |

## Tautological tests (change-pins, documented)

- `tests/regression/test_light_ring_identity.py::test_epsY_nulltest`
- `tests/regression/test_light_ring_identity.py::test_epsY_nulltest_uses_f2Y_symbolic`
- `tests/regression/test_light_ring_identity.py::test_ward_residual_electric_channel_JA`
- `tests/regression/test_kt_mh_background.py::test_a4_v8_slot_definitions`
- `tests/regression/test_kt_mh_background.py::test_ward_span_derivation_status`
- `tests/regression/test_eq85_naming_guard.py::test_v_slot_mapping_is_documented`
- `tests/regression/test_eq85_naming_guard.py::test_q_conventions_must_not_be_mixed`
- `tests/unit/test_inner_endpoint_contract.py::test_endpoint_orders_are_generated_from_contract_only`
- `tests/unit/test_inner_endpoint_contract.py::test_endpoint_equation_count_matches_contract`
