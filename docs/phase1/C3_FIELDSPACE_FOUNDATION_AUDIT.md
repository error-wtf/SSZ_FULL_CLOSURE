# C3 FIELD-SPACE FOUNDATION AUDIT (PHASE 1.5)

Status: stability-blind theoretical audit, completed before Phase 2.
Date: 2026-09-22. Scope: C3 only. No K, c_r^2, c_Omega^2, QNM, tangent-sensitivity,
historical stability or PASS/FAIL data was loaded or inspected (see
PHASE1_BLINDING_AUDIT.json; document-level grep check repeated for this audit).

## 1. Field-space metric audit

Manifold: the field space of the canon-closed P5 family, fields (g_mu_nu,
phi_P5, A_mu), restricted to the static spherical electric branch; on-curve
coordinate = arclength s of the background curve.

Chart: (phi_P5, X, F, Y).
  - phi_P5: the segmentation scalar (Monograph anchor "scalar phi=Xi").
  - X: scalar kinetic invariant of the scalar.
  - F: Maxwell invariant F_mu_nu F^mu_nu (electric branch, F-tilde = 0).
  - Y: mixed invariant; on-curve the P5 code uses Y = 4 X F as shorthand
    (strong_field_continuation.py); off-curve Y is an independent coordinate
    of the chart unless a candidate fixes otherwise.

Privilege: operational, not geometric. Every coupling function of the 41-slot
emitter is a function of exactly these invariants (mh_general_primitives.py,
full_action_lower.py); the canon uses no other chart. Under spacetime
diffeomorphisms all four are scalars; under U(1) all four are invariant.
Under field redefinitions they transform but are NOT unique -- the chart is a
representation choice.

Audit result (symbolic, tools runs recorded in the session log):
  - Linear reparametrisation Q = a q: third jets transform linearly,
    gamma_FS transforms tensorially, J_C is a scalar. PASS.
  - Nonlinear reparametrisation Q = q + alpha q^2: the third derivative picks
    up inhomogeneous terms (verified: 6 alpha u' u'' for the one-coordinate
    case). The third jet is NOT a tensor; it is an affine jet-bundle object.
    The naive norm gamma_FS(J^3, J^3) with flat component matrix delta in the
    privileged chart is therefore NOT reparametrisation-invariant. FAIL.

  Repair (declared, not silent): the flat field-space CONNECTION nabla^FS of
  the privileged chart is added to the theory definition. Covariant third jets
  are defined as jets relative to nabla^FS (equivalently: components in any
  nabla^FS-flat chart plus the standard connection correction terms). With
  this declaration the functional is chart-independent by construction.

  Gamma_FS = delta_AB in the privileged chart is thereby a genuine tensor
  (components in one chart + declared connection), but its flatness is a
  postulate of the theory, not a derivation from SSZ.

## 2. Dimensional / jet-scaling audit

  Lambda := r_s (the only scale of the canon background; convention lock
  r_s = 1). In canon units every 41-slot coefficient is dimensionless; a k-th
  jet carries dimension L^(-k); the weight Lambda^(2k) = r_s^(2k) renders each
  summand dimensionless. No other power law is dimensionally consistent with a
  single scale, hence no residual freedom.
  Origin of Lambda: the canon unit convention. NOT stability-derived.
  Minimizer dependence: within a fixed unit system Lambda is a constant, so
  J_C scales by an overall factor under unit changes and the argmin is
  unchanged; the relative weights between jet orders are fixed by the single
  canon scale. No C3_UNDERDETERMINED_SCALE.

## 3. Parametrisation invariance test

  T1 (linear, Q = a q): J_C' - J_C = 0. C3_FIELDSPACE_COVARIANCE (linear): PASS.
  T2 (nonlinear, Q = q + alpha q^2): naive jets FAIL (inhomogeneous terms);
  with the declared nabla^FS: covariant jets transform correctly, J_C scalar.
  Outcome per directive section 3: B.
  C3_FIELDSPACE_METRIC_IS_ADDITIONAL_POSTULATE = true.

## 4. Strict convexity / uniqueness

  Function space: the affine space of admissible covariant third transverse
  jets J^(3 perp)(s) along the curve (real analytic profiles), with the lower
  jets fixed along the whole curve by (a) background EOM JA=E00=E11=0,
  (b) holonomicity chain-rule identities, (c) convention identities
  (v12 = -v6/(2h); corrected a5), (d) regularity and asymptotics.
  Boundary conditions: fixed lower jets everywhere; no endpoint third-jet data
  is declared by the canon.
  Second variation: delta^2 J_C = 2 integral ds gamma_FS(delta J^3, delta J^3)
  > 0 for every nonzero variation (positive-definite metric). Strictly convex;
  the minimiser is unique in the affine admissible set. No null directions
  remain after the constraint projection (U(1) gauge is already eliminated in
  the reduced formalism; residual redundancy: none declared).

  NONTRIVIALITY RISK (declared, not hidden): with no endpoint data for the
  third jets, the unique minimiser is the vanishing covariant third jet.
  This is NOT automatically the falsified C0 member: C0 zeroed the bundle
  displacement amplitudes (a1/c4 tangent directions, polynomial in u),
  whereas C3's zero acts on field-space Taylor jets -- different constraint
  languages, possibly different 41-slot streams. Whether the two streams
  coincide is a mechanical Phase-2 construction question; both outcomes are
  reportable without selection.

## 5. Definition of J^(3 perp)

  Transversality is defined METRIC-FREE: let N_iA(q) be the covectors obtained
  from the transverse field-space derivatives of the three background
  equations (JA, E00, E11) and of the holonomicity identities, evaluated on
  the curve. Then
      J^(3 perp) in ker N(q)   (componentwise along the curve).
  gamma_FS enters ONLY the norm, not the projector. This avoids the hidden
  circularity in which the same metric defines both transverse directions and
  the minimiser. Declared explicitly: with this choice gamma_FS does one job
  (norm), the constraint covectors do the other (directions).

## 6. Euler-Lagrange / quintic statement

  Variation of J_C = integral ds gamma_AB u'''^A u'''^B with fixed endpoint
  2-jets (u, u', u'' fixed at both ends):
    delta J_C = [gamma_AB u'''^B delta u''^A - d/ds(gamma_AB u'''^B) delta u'^A
                 + d^2/ds^2(gamma_AB u'''^B) delta u^A]_0^L
                - integral ds d^3/ds^3(gamma_AB u'''^B) delta u^A.
  Boundary terms vanish for fixed endpoint 2-jets. EL:
      d^3/ds^3( gamma_AB u'''^B ) = 0.
  In the privileged chart gamma_AB = delta_AB is constant, hence
  u^(6)^A = 0 componentwise: every minimiser is a quintic per component --
  the exact analogue of the canonical bridge theorem Xi^(6) = 0 (Variational
  Bridge v2: J[u] = int (u''')^2 dt over H^3 with six endpoint data; unique
  quintic minimiser). The field-space version preserves the order EXACTLY
  because the frozen connection makes gamma constant in the privileged chart;
  for a nontrivial gamma(s) the EL would be the full variable-coefficient
  third-order equation in u''' and the quintic statement would not hold --
  this case is not part of the frozen theory.

## 7. Comparison to C1

  None performed. C1 and C3 remain parallel pre-registered hypotheses. No
  ranking, no survivor selection, no expected-stability reasoning.

## 8. Freeze result

  FINAL CLASSIFICATION:
      C3_NEW_THEORY_HYPOTHESIS_WITH_FIXED_FIELDSPACE_METRIC

  C3_FIELDSPACE_COVARIANCE = PASS (with the declared connection; naive-jet
  covariance FAIL is recorded as the reason the connection is required).
  C3_FIELDSPACE_METRIC_IS_ADDITIONAL_POSTULATE = true.

  Theory definition of C3 (all items enter the theory hash):
    - invariant chart (phi_P5, X, F, Y) with on-curve shorthand Y = 4XF;
    - flat field-space connection nabla^FS of that chart;
    - gamma_FS = delta_AB components in that chart;
    - Lambda = r_s, jet scaling r_s^(2k), canon unit convention r_s = 1;
    - metric-free transversality projector ker N(q) (constraint covectors of
      JA, E00, E11, holonomicity identities);
    - boundary conditions: lower jets fixed everywhere by EOM/holonomy/
      conventions/regularity/asymptotics; no endpoint third-jet data;
    - functional J_C = integral ds r_s^(2k) gamma_FS(J^(3 perp), J^(3 perp)).

  Phase-2 obligations recorded: mechanical stream comparison of the C3 member
  against the C0 witness (nontriviality check, section 4), full blind forward
  replay of C1 and C3, whole-distribution reporting, no post-unblinding
  modification.

## 9. Phase 2 authorisation

  With this audit frozen and hashed in a new preregistration commit, Phase 2
  begins: C1 and C3, each independently through the complete forward chain.
  Do not stop when one fails; do not prefer one when one passes; report both.
