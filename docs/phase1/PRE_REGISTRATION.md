# PRE-REGISTRATION: SSZ COMPLETION ENSEMBLE C0–C6 (PHASE 1)

Status: PRE_REGISTRATION_FROZEN (pending commit)
Date: 2026-09-22
Rule of construction: NO stability information (K, lambda_min(K), c_r^2, c_Omega^2,
angular roots, QNM, published perturbation roots, historical PASS/FAIL of tangent
directions, tangent sensitivity gradients) was used in the definition, selection,
parametrisation or ranking of any principle below. The forbidden artifacts are
listed in PHASE1_BLINDING_AUDIT.json.

## 0. Problem statement

The certified SSZ P5 background fixes the on-curve action data but leaves the
transverse action jets underdetermined. Two concrete completions are already
frozen as falsification witnesses:

  C0  TRANSVERSE_ZERO (c = 0):                    K_GHOST_FAIL
  C0b K-repaired historical member:               ANGULAR_FAIL_ROBUST

Neither result is reused here. Phase 1 asks exactly one question:

  Does SSZ itself supply an independent, covariant completion principle for
  the transverse action jets?

Direction of derivation (binding):

  SSZ principle -> S[g, phi_P5, A] -> Xi(r), D(r)        (forward; required)
  Xi(r), D(r)   -> any S that fits                       (inverse; FORBIDDEN)

## 1. Canonical audit performed before freezing (directive section 11)

### 1.1 Two distinct objects share the symbol "phi"

  phi_gold = (1+sqrt(5))/2 = 1.618033988749895
      Structural constant of the SSZ phenomenology layer. Sources:
      02_FOUNDATIONS/segment_density.md, structural_constants.md, phi_geometry.md,
      formula_compendium.md, SSZ_CANONICAL_FORMULAS_2026.md,
      SSZ_Middle_Bridge_Unified_Dynamics (Xi_s = 1 - exp(-phi/x)).
      Appears as: exponent of the strong-branch Xi; coupling radius
      r_phi = (phi/2) r_s [1 + beta Delta(M)]; frequency grid R = phi^N;
      logarithmic spiral b = ln(phi)/Delta theta.

  phi_P5 = Xi(C)
      The dynamical segmentation scalar of the P5 action (Monograph convention
      anchor, SHA a028d69bd8e6..., chunk 86: "scalar phi=Xi"; C = r_s/r,
      x = r/r_s = 1/u; chi = dXi/dC). NOT the golden ratio. The Monograph lists
      the historical Xi_strong = 1 - exp(-phi_gold r_s/r) as historical
      provenance only ("historical Hermite middle bridge is provenance only and
      is not the P5 field equation").

  RULING (frozen): phi_gold is a fixed structural constant of the SSZ layer and
  is NOT a fit parameter anywhere in this ensemble. phi_gold and phi_P5 are
  never numerically identified. Any completion that needs phi_gold inside a
  covariant action must treat its appearance as NEW_THEORY_HYPOTHESIS (there is
  no canon derivation of the golden ratio inside the P5 action family).

### 1.2 Xi_max contradiction (resolved)

  - SSZ_CANONICAL_FORMULAS_2026.md: "Xi_max = 1.000000 (exactly)".
  - glossary_de.md, phi_geometry.md, full-segmented-memory.md: Xi_max = 0.80171
    = 1 - exp(-phi_gold).

  Resolution per CANONICAL_XI_AND_BLEND_RESOLUTION.md: the min-form
  Xi_sat = min(1 - exp(-phi_gold r_s/r), Xi_max) is a LOCAL metric-pure
  context form. The canonical GLOBAL form is the decay form
  Xi_strong = 1 - exp(-phi_gold r_s/r), which needs no Xi_max (saturation is
  intrinsic). The Xi_max = 1.0 reading is additionally excluded from any global
  channel because min(., 1) tends to 1 at infinity, violating Xi -> 0.

  RULING (frozen): the decay form is the canonical strong branch for this
  ensemble. Xi_max = 1 - exp(-phi_gold) = 0.801712 appears only as the derived
  horizon value. The Xi_max = 1.000 reading is excluded from the P5 channel.

### 1.3 Radial coordinate conventions

  P5: C = r_s/r, x = r/r_s = 1/C, u = C in the dense tooling.
  Middle Bridge: u_MB = GM/(r c^2) = C/2.
  All formulas below are written in C; conversions are documented in
  SSZ_P5_COMPLETE_STATUS_TO_KRGM_QNM_2026-09-16.md.

### 1.4 Construction constraints (directive section 12)

  The following are THEORIES POSTULATES of the SSZ canon and are therefore
  CONSTRUCTION_CONSTRAINTS for every candidate. Their reproduction by any
  member is NOT an independent prediction:

    Xi_weak = r_s/(2r)                                  (weak branch)
    Xi_strong = 1 - exp(-phi_gold r_s/r)                (strong branch, decay form)
    Hermite C2 quintic bridge on x in [1.8, 2.2]        (bridge; J[u]=int (u''')^2 dt)
    phi_gold = 1.618033988749895                        (structural constant)
    U(1) gauge invariance of the vector sector          (Zhang-Kase schema)
    second-order / degenerate field equations           (Horndeski class)
    weak-field GR/PPN limit (beta = gamma = 1)          (canon: exact)
    regularity and asymptotics                          (canon guardrails)

  K, c_r, c_Omega, QNM are candidate INDEPENDENT_PREDICTIONS of any candidate
  member, contingent on none of them having influenced the completion choice.
  This contingency is satisfied by construction below and is re-verified in
  Phase 2.

### 1.5 Canon action-family anchor

  The covariant family is closed and canonically documented:
    - SSZ_P5_FULL_CLOSURE_MONOGRAPH_2026-09-16, chapter 5:
      S = integral d^4x sqrt(-g) (L_H[g, phi] + L_SVT[g, phi, F]),
      hybrid composition rule C_HSVT = C_MH + (C_SVT - C_shared),
      "Reduce(A+B) != Reduce(A)+Reduce(B)".
    - SSZ_P5_Horndeski_SVT_CLOSURE_COMPLETE_2026-09-14: full L2-L5 Horndeski
      realization plus gauge-invariant SVT sector with same-action closure
      theorem for the regular P5 geometry.
    - 41-slot coefficient schema a1..a9, b1..b5, c1..c6, d1..d4, e1..e4, v1..v13
      (Zhang-Kase 13-slot vector block after canonicalisation).
    - Background EOM convention: JA = 0, E00 = 0, E11 = 0.

## 2. The frozen ensemble

### C0 — HISTORICAL BASELINE / TRANSVERSE_ZERO  (control; not modified)

  Archived falsification witness. Member = unmodified background-forced stream
  (c = 0 in the background-null tangent bundle). Frozen status: K_GHOST_FAIL.
  Registered only as control; not re-optimised, not repaired.

### C1 — MINIMAL COVARIANT OPERATOR COMPLETION

  Rule (lexicographic, applied BEFORE any stability evaluation):
  Within the canon-closed Horndeski + U(1)-SVT family, select the covariant
  operator content minimising, in order:
    (1) number of independent operators,
    (2) maximal derivative order,
    (3) analytic/polynomial degree in the field invariants (phi_P5, X, F, Y),
    (4) number of additional free couplings,
  subject to: exact satisfaction of the canonical SSZ background (JA=E00=E11=0
  on the canonical profile), the weak/strong regime structure, U(1) gauge
  invariance, second-order/degenerate equations, regularity, asymptotics.
  NO stability criterion. No coefficient minimisation in an arbitrary basis.
  The member must follow from the background EOM and fixed boundary
  conditions; FAIL_BACKGROUND is a possible and acceptable outcome.

  Classification: CANON_MOTIVATED_HYPOTHESIS (the lexicographic order is a
  methodological choice; the family is canon).
  Expected uniqueness: unique operator set -> unique member, contingent on the
  background solve. Possibility of no admissible operator set is recorded as
  FAIL_BACKGROUND, not fitted around.

### C2 — ANALYTIC SATURATION COMPLETION

  Motivation (canon): the strong branch is Xi_strong = 1 - exp(-phi_gold C):
  exponential saturation is intrinsic. The generator structure D = exp(sigma),
  sigma = -ln(1+Xi), suggests couplings whose off-curve extension is fixed by
  the same analytic generator that produces the on-curve exponential.

  Requirement: couplings defined as covariant functions of the legitimate
  field invariants (phi_P5, X, F, Y) only. Writing Xi(r) into an action is
  forbidden; inverse-engineering from the desired background is forbidden.

  Assessment (theoretical, no stability): infinitely many covariant functions
  of (phi_P5, X, F) reduce to the same on-shell exponential in C; the canon
  contains no criterion singling out one extension; embedding phi_gold into
  the action family would itself be NEW_THEORY_HYPOTHESIS (see 1.1).

  STATUS: UNDERDETERMINED (pre-registered). No member is generated. Not
  decided by stability later. Classification: CANON_MOTIVATED_HYPOTHESIS.

### C3 — VARIATIONAL COMPLETION (frozen functional)

  Motivation (canon): the bridge principle J[u] = integral_0^1 (u'''(t))^2 dt
  with admissible class H^3 and six boundary data; Euler-Lagrange Xi^(6) = 0;
  quintic Hermite is the unique global minimiser (Variational Bridge v2,
  SHA 2c93df70..., section on the canonical construction). The canon confines
  this principle to one radial dimension and warns against naive promotion
  (section 20). C3 is the minimal consistent promotion to the transverse jets.

  Frozen functional: let J^3_T denote the third transverse jet bundle along
  the background curve, i.e. all jet components of the coupling functions
  (in the invariant chart (phi_P5, X, F, Y)) that are NOT fixed by
    (a) the background EOM (JA=E00=E11=0, on-curve),
    (b) holonomicity/chain-rule identities (canon: df3/dr = [f3]_phi phi_r + ...),
    (c) the accepted convention identities (v12 = -v6/(2h); corrected a5),
    (d) regularity and asymptotic conditions.
  Define on the curve, parameterised by arclength s,

      J_C = integral ds  gamma_FS( J^3_T(s), J^3_T(s) ),

  where gamma_FS is the flat field-space metric with unit coordinate norm in
  the canonical invariant chart (phi_P5, X, F, Y), i.e.
  gamma_FS = d phi_P5^2 + d X^2 + d F^2 + d Y^2 restricted to jet indices,
  scaled per jet order by the dimension factor L^(2k) for k-th jets so that
  every summand is dimensionless.

  Justifications (frozen):
    - Covariance: gamma_FS is a genuine tensor on field space; the contraction
      of a jet-bundle tensor with a field-space metric is invariant under
      field redefinitions. The chart is used only for the component statement,
      not for the invariance.
    - Dimensions: jet scaling makes J_C dimensionless; no new dimensionful
      constant enters.
    - Redefinition dependence: none beyond the declared gamma_FS; the choice
      of the flat gamma_FS in the invariant chart is exactly the
      NEW_THEORY_HYPOTHESIS component of C3 and is frozen here.
    - Relation to canon: minimising the squared third jet reproduces, in the
      one-jet case, the bridge theorem's H^3 structure (EL yields polynomial
      jets of bounded degree, the direct analogue of Xi^(6) = 0).

  Unique member: YES by construction (strictly convex quadratic functional in
  the free jet components; minimiser unique in the affine constraint set).
  Classification: CANON_MOTIVATED_HYPOTHESIS with one declared
  NEW_THEORY_HYPOTHESIS component (gamma_FS). No stability information used.

### C4 — MINIMAL ANALYTIC JET COMPLETION

  Requires a basis-invariant meaning of "minimal jet". Such a meaning needs a
  field-space measure, i.e. exactly the object C3 had to postulate. Without
  it, minimality is basis-dependent, and freezing any basis would smuggle the
  falsified TRANSVERSE_ZERO back in as a "principle" (explicitly forbidden).

  STATUS: NOT_WELL_DEFINED (stand-alone). If the frozen gamma_FS of C3 is
  accepted, C4 collapses into a lower-jet variant of C3 and is not a distinct
  principle; it is registered as non-distinct, not as a separate member.
  Classification: DERIVED_FROM_CANON (the negative assessment follows from the
  canon's own redefinition-sensitivity discipline).

### C5 — STRUCTURAL GAUGE / NED COMPLETION

  Canon sources: Middle Bridge sections 21.7-21.9 (SHA 40ba5156... /
  b616c0f9...): NED diagnostics with L_F zeros (x = 2.034986003262,
  2.199789272082; 1.927207549701, 2.133996087105), defect-to-gauge split
  rho = rho_D + rho_G with stress signatures diag(-1,-1,0,0) / diag(-1,-1,1,1),
  rejected single-sector filters. The canon states verbatim that these are
  "action-search filters, not matter identifications" that "reduce the
  candidate space while leaving the fundamental inverse problem open".

  STATUS: FAMILY_NOT_UNIQUE. Free dimensions (documented): the NED function
  L(F, G) family, the constitutive relation, the defect/gauge split. The
  filters exclude sectors (e.g. pure linear radial Maxwell requiring p_t =
  rho exactly) but do not select an action. No member generated.
  Classification: DERIVED_FROM_CANON.

### C6 — SYMMETRY / GENERATOR COMPLETION

  Canon sources: Middle Bridge U(1) chapter ("U(1) = {z in C : |z|=1} ... the
  exponential map sends the additive generator iPhi to the group element
  e^{iPhi}"); JIF Foundation Canonical (phase chain dJ = dPhi/(2pi),
  dJ/dt = f0 D; no unreferenced absolute-phase observable; U(1) without gauge
  boson; no force on light).

  Assessment: the structural parallel sigma = ln D, D = exp(sigma) (R^+
  generator of rescaling clocks) with U(1) = exp(i Phi) (phase generator) is
  formal similarity only. The canon contains no Noether current, no
  connection, and no field-space structure that joins them into an action
  principle; the JIF foundation explicitly refuses absolute-phase physics.

  STATUS: INTERPRETIVE_ONLY_NOT_ACTION_PRINCIPLE. No member generated.
  Classification: DERIVED_FROM_CANON.

## 3. What is frozen here

  - The six principle definitions above, with statuses
      C0  control (archived witness)
      C1  DEFINED, unique member expected (Phase 2 forward replay; FAIL_BACKGROUND possible)
      C2  UNDERDETERMINED
      C3  DEFINED (unique minimiser), one declared NEW_THEORY_HYPOTHESIS component
      C4  NOT_WELL_DEFINED (stand-alone); non-distinct from C3 if gamma_FS accepted
      C5  FAMILY_NOT_UNIQUE
      C6  INTERPRETIVE_ONLY_NOT_ACTION_PRINCIPLE
  - The phi audit rulings (1.1, 1.2) and the construction-constraint list (1.4).
  - The blinding list of forbidden stability artifacts (PHASE1_BLINDING_AUDIT.json).

  Phase 2 (after PRE_REGISTRATION_FROZEN commit): blind forward replay of C1
  and C3 through the full chain (action -> background -> 41-slot -> constraints
  -> K/G/S/M -> finite-L K -> radial -> angular -> odd -> global -> QNM), each
  separately, no member changed after seeing results, whole-distribution
  reporting. C2, C4, C5, C6 generate no members and are reported as-is.
