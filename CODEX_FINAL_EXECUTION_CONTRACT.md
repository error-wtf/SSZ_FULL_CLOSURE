# SSZ P5 — CODEX EXECUTION CONTRACT FOR ABSOLUTE FULL CLOSURE

## 0. Status of this handoff

This repository is a **self-contained handoff** assembled from the complete 2026-09-16 research snapshot plus the locked regional-production correction used on 2026-09-17. The purpose of this file is to prevent repeated research, repeated filesystem searches and repeated validation of already closed subproblems.

The repository has deliberately been pre-staged with the files that previous Codex runs kept treating as missing:

```text
data/authoritative/ssz_p5_F2a_EXPLICIT_GLOBAL_ACTION_COVER_2026-09-15.csv

data/prestaged/direct41/outer_resolved_raw_41of41.csv
data/prestaged/direct41/central_selected_corrected_41of41.csv
data/prestaged/direct41/inner_selected_candidate_41of41.csv
data/prestaged/direct41/core_selected_41of41.csv
data/prestaged/direct41/weak_exterior_39of41.csv

data/prestaged/odd/central_svt_odd_profile.csv
```

These files exist so you do **not** spend time locating or reconstructing historical outputs merely to understand the already completed work. Their roles are encoded in `src/ssz_p5/production/sources.py`.

A PRESTAGED or REGRESSION file is not automatically a Direct-Closure certificate. Use it as a target/reference while implementing or verifying the final direct generator.

---

# 1. Treat the validated project state as input

You are not beginning the physics from first principles.

The project has already gone through repeated derivation, debugging, counterchecks and regression work. Historical files are retained for provenance and therefore contain superseded branches. That archival complexity is not permission to reopen every resolved question.

The operating rule is:

```text
validated project state
+ narrowly scoped implementation
+ regressions
+ final integration tests
```

not:

```text
doubt every prior result
→ derive it again
→ search literature again
→ reopen the member
→ consume the budget
```

A current PASS remains accepted until your changed code causes a relevant current regression to fail.

---

# 2. Authoritative production-region map

Use only `src/ssz_p5/production/regions.py` for production membership.

The locked map is:

```text
u < 0.5515230871346237
    weak_exterior_H

0.5515230871346237 <= u < 0.61
    outer_same_action_H_SVT

0.61 <= u < 0.71
    central_exact_SVT

0.71 <= u < 0.715
    inner_same_action_SVT_H

u >= 0.715
    punctured_H_core

r = 0
    analytic_center
```

Do not infer production membership from filenames, radial coverage of witnesses, or old prose.

Two previously misclassified radii are explicitly regression-tested:

```text
x ≈ 1.5999977950251312 → central_exact_SVT
x ≈ 1.41616064         → central_exact_SVT
```

A pure-Horndeski `A0prime=0` identity is NOT_APPLICABLE there because the selected production member is genuine-SVT with nonzero vector background.

Do not resurrect the old false global contradiction.

---

# 3. Strong-H carrier role is fixed

The strong-H carrier is a:

```text
PRINCIPAL WITNESS
CONTROL WITNESS
HORNDESKI REGRESSION TARGET
```

It is not the global production background through `0.551523... <= u < 0.715`.

The fact that a witness table contains a radius does not make it the selected production member at that radius.

Do not use it to overwrite the outer H/SVT handover or central SVT lobe.

---

# 4. Do not globally impose the later A0prime=0 simplification

A later research branch explored a global Horndeski carrier plus background-null `epsilon_Y Y` with global `A0prime=0`. That branch is not the selected regional Full-SVT production architecture you are closing here.

Do not replace the regional chain

```text
weak H
→ outer H+SVT handover
→ central genuine SVT
→ inner SVT+H handover
→ H core
→ analytic center
```

with a global null-vector Horndeski member.

---

# 5. Precedence when files disagree

Use this precedence order:

```text
1. src/ssz_p5/production/regions.py
2. src/ssz_p5/production/sources.py
3. current re-solved outer/inner production background datasets
4. current central genuine-SVT datasets
5. current action-level core datasets
6. current direct coefficient generators
7. current constraint implementation
8. current profile-aware reducer
9. prestaged/regression datasets
10. witness/control datasets
11. diagnostics
12. historical/superseded files and old prose
```

A filename containing `FINAL`, `COMPLETE`, `MASTER` or `SELECTED` does not override this precedence.

---

# 6. Closed inputs — do not re-derive

Do not reopen these unless a correctly classified downstream production test points directly to them:

```text
P5 geometry
Xi geometry
light-ring locations
regional production assignment
outer background construction
central genuine-SVT member construction
inner background construction
punctured-core construction
analytic-center strategy
41-slot schema
JET9D8 design
corrected a5 convention
corrected v12 convention
selected lower-order member
shared-baseline principle
common constraint architecture
b2 constraint correction
profile-aware reducer architecture
scalar WKB trapping
Maxwell WKB trapping
outgoing Jost asymptotic validation
historical QNM branch reconciliation
```

Do not repeat a derivation simply because an old file shows an earlier convention.

---

# 7. Locked coefficient conventions

The common even-parity unreduced schema is exactly:

```text
a1..a9
b1..b5
c1..c6
d1..d4
e1..e4
v1..v13
```

Total: 41.

The accepted branch uses:

```text
v12 = -v6/(2h)
```

and the corrected holonomic expression

```text
a5 = a2' - a1'' - (A0' v4 / 2)' + A0' v5 / 2
```

For a pure-H region with `A0prime=0`, this reduces to:

```text
a5 = a2' - a1''
```

Do not revisit these signs/formulas.

The current selected lower-order representative remains the encoded one; do not launch another lower-order optimization.

---

# 8. JET9D8 is the production derivative service

Use the established nonuniform local polynomial derivative service:

```text
window = 9
degree = 8
```

Do not replace it with splines or a new smoothing method in production.

Alternative windows/degrees may be used only in the already-defined convergence checks.

---

# 9. Source registry — no filesystem archaeology

The canonical paths are already encoded in:

```text
src/ssz_p5/production/sources.py
```

Use that registry in new production code. If you locate another useful artifact once, add its canonical path to the registry and stop repeatedly searching for it.

Before asking the authors for a file, check the source registry and `data/prestaged/`.

---

# 10. Weak exterior region

Domain:

```text
u < 0.5515230871346237
```

A complete 39/41 weak-exterior regression is already pre-staged at:

```text
data/prestaged/direct41/weak_exterior_39of41.csv
```

The missing lower-order slots are not an invitation to search for a new action. Complete them using the already selected lower-order convention and current corrected formulas.

Use the existing Horndeski/Maxwell-Horndeski emitter and current JET9D8 implementation.

Do not change the weak background.

---

# 11. Outer same-action H/SVT handover

Domain:

```text
0.5515230871346237 <= u < 0.61
```

Authoritative re-solved background:

```text
data/production/ssz_p5_F2_outer_same_action_RESOLVED_background_jets_2026-09-15.csv
```

A full previously generated outer 41-slot stream has already been pre-staged:

```text
data/prestaged/direct41/outer_resolved_raw_41of41.csv
```

This is your direct regression target. Do not ask for the old `ssz_p5_OUTER_RESOLVED_ZK_RAW_41of41_2026-09-16.csv`; its contents are already here.

The entire outer domain starts at `u=0.5515230871346237`. Do not insert the Strong-H witness over `0.551523..0.57`.

Use the existing ZK/SVT action-jet emitter for this electrical H/SVT handover.

The background is already re-solved; do not repeat the background inverse problem.

Known re-solved residuals are at floating-point scale and remain regressions.

---

# 12. Central exact genuine-SVT region

Domain:

```text
0.61 <= u < 0.71
```

Use the existing central Full-SVT member.

Authoritative data already included:

```text
data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv
data/production/ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv
```

A normalized/corrected 41-slot reference is already pre-staged:

```text
data/prestaged/direct41/central_selected_corrected_41of41.csv
```

The central Odd profile is also already pre-staged:

```text
data/prestaged/odd/central_svt_odd_profile.csv
```

Do not search for `ssz_p5_final_svt_lobe_odd_profile_2026-09-12.csv`; the file is already copied here.

Do not replace this region with a Horndeski null-vector member.

Do not rebuild the central lobe.

---

# 13. Inner same-action SVT/H handover

Domain:

```text
0.71 <= u < 0.715
```

Authoritative re-solved background:

```text
data/production/ssz_p5_F2_inner_same_action_RESOLVED_v2_background_jets_2026-09-15.csv
```

A deterministic selected 41-slot candidate is pre-staged:

```text
data/prestaged/direct41/inner_selected_candidate_41of41.csv
```

This candidate exists to give you the expected coefficient stream and endpoint behavior. It is not permission to call the direct higher-jet regeneration closed if that direct provenance has not been produced.

If higher transverse action jets needed for independent regeneration are missing from the background table, reconstruct only those missing higher jets from the existing partition/endpoints/tubular-normal/higher-jet closure machinery.

Do not search for a new background or new member.

Do not change the region map.

---

# 14. Punctured Horndeski core

Domain:

```text
u >= 0.715
```

The selected 41-slot core reference is already pre-staged:

```text
data/prestaged/direct41/core_selected_41of41.csv
```

The action-level core inputs already included are:

```text
data/production/ssz_p5_F1b_G4XX_transverse_core_candidate_2026-09-14.csv
data/production/ssz_p5_F1b_FULL_G5_subcore_onshell_candidate_2026-09-14.csv
```

The known software gap from the research history is narrow: the simple `G4(phi), G5=0` emitter is not a general `G4X/G4XX/G5` Horndeski core emitter.

If the final direct provenance requires a general core emitter, implement/port exactly that action-jet→41-slot functionality and no broader theory redesign.

The new general emitter must reduce to the already validated restricted emitter when the general extra jets are set to their restricted limit.

Use the pre-staged core 41 table as a regression target. Do not back-fit a new action from it.

---

# 15. Analytic center

Use:

```text
data/production/ssz_p5_F1b_FINAL_Cinf_center_to_punctured_handover_2026-09-14.csv
```

and the existing Taylor/analytic center rules.

Do not evaluate formulas at `r=0` that manufacture divisions by `r`, `phi_r`, `X` or `kappa`.

The exact center is an analytic override, not a numerical finite-difference discovery problem.

---

# 16. Global action cover already supplied and reproducible

An older action-cover CSV is already present:

```text
data/authoritative/ssz_p5_F2a_EXPLICIT_GLOBAL_ACTION_COVER_2026-09-15.csv
```

A current cover can be regenerated directly from the locked Python region/source registry by running:

```bash
python tools/generate_global_action_cover.py
```

Do not treat absence of an old byte-identical action-cover file as a blocker.

---

# 17. Build one global direct 41/41 stream

After each regional direct stream is available, assemble exactly one global unreduced coefficient stream according to `regions.py`.

Required output conceptually:

```text
data/authoritative/ssz_p5_GLOBAL_DIRECT_41of41.csv
```

Each row should carry:

```text
radial coordinates
production_region
source/provenance id
all 41 coefficient slots
```

Do not stitch already reduced K/G/M matrices.

Do not use the old `SELECTED_41STREAM_V2` as production; it is blacklisted because its old regional assignment included the Strong-H shortcut.

---

# 18. Interface validation

Check only the required production interfaces:

```text
u = 0.5515230871346237
u = 0.61
u = 0.71
u = 0.715
center limit
```

Use the existing endpoint/continuity regressions.

Do not turn interface validation into a new member search.

Once interface tests pass, freeze/hash the global direct 41 stream and move forward.

---

# 19. Shared baseline rule is final

Do not reduce sectors independently and add the reduced operators.

Correct sequence:

```text
assemble unreduced common action/coefficient system
→ eliminate common auxiliary fields once
→ reduce
```

Where sector subtraction notation is used:

```text
C_total = C_H + (C_SVT - C_shared)
```

Use the existing shared-baseline validation. Do not re-derive the principle.

---

# 20. Constraint system is already implemented

Use the current constraint implementation, including the corrected `b2` contribution.

Do not compare against old constraint implementations unless a current stationarity test fails.

Check the defined pivots (`Dh1`, `DeltaV`, `2 L v9`, etc.) with the analytic-center handling already adopted.

---

# 21. Use the existing profile-aware reducer

Do not build another reducer.

The accepted reduced action convention is:

```text
L2 = dot(Y)^T K dot(Y) - Y'^T G Y' + Y'^T S Y - Y^T M Y
```

with the project canonical `R` convention and the expected matrix symmetries/antisymmetries.

All radial basis derivatives, product-rule terms and integration-by-parts terms must remain included.

---

# 22. Required finite-l production set

At minimum use:

```text
L = 6, 12, 20, 42, 110, 420, 1000
```

The same global member must be used for all multipoles.

Generate K,R,G,S,M artifacts with provenance linking them to the global direct 41 hash and numerical policy.

---

# 23. Finite-l gates

For every required multipole, check the already defined project gates only:

```text
finite matrices
matrix structure
constraint pivots
kinetic positivity where required
physical radial generalized eigenvalues
pure-sector regressions
central exact-SVT regression
high-L M regression
JET9D8 convergence
interface behavior
```

Do not invent new fundamental closure criteria during implementation.

If one gate fails, localize by region/radius/L/coefficient before considering any scientific redesign.

---

# 24. Odd and Even both receive production certificates

Use the supplied/pre-staged central Odd profile and the already established odd-sector conditions. Do not re-derive Odd perturbation theory.

Produce machine-readable current-production certificates for Odd and Even after the direct regional/global products are validated.

---

# 25. QNM begins only after the direct reduced operator is closed

Do not perform more proxy QNM work while Direct KRGSM is incomplete.

When the direct operator is available, use exactly that same K,R,G,S,M basis for the final coupled spectral system.

The repository already contains:

```text
scalar/Maxwell trapping diagnostics
outgoing Jost asymptotic validation
Jost order/radius convergence data
complex-root diagnostic work
historical branch reconciliation
```

Use those as solver/boundary regressions and initial guesses. Do not redo them from scratch.

Real-axis inward shooting of damped modes remains rejected as the final production method.

Prefer extending the existing Jost machinery to a coupled compactified spectral determinant, or use ECS if the existing code makes that substantially cheaper.

Do not compare five different QNM frameworks.

---

# 26. QNM convergence policy

A numerical root is accepted only after the already expected convergence controls are documented, including relevant combinations of:

```text
radial resolution
outer domain / compactification
Jost asymptotic order
root tolerance
basis normalization
branch tracking
```

If ECS is used, include ECS angle/onset convergence.

Explicitly scan the declared production frequency domain for `Im(omega)>0` and record the scanned domain. Do not claim a finite search is a theorem over the entire complex plane.

---

# 27. Failure-handling protocol

If a current mandatory production gate fails, record:

```text
gate
region
radius
L
coefficient/matrix element/eigenvalue
observed value
expected condition
immediate dependency chain
```

Check first for:

```text
wrong production region
wrong source path/column mapping
wrong derivative
wrong sign convention
wrong field basis
missing profile derivative
missing constraint term
missing action jet
double-counted shared baseline
```

Repair the smallest responsible layer.

Only reopen a scientific construction after the failure is shown to be region-correct, formula-correct, numerically converged, basis-correct and independently reproducible.

---

# 28. Test execution discipline

At the beginning, run exactly once:

```bash
python tools/generate_global_action_cover.py
python tools/verify_codex_handoff.py
pytest -q
python ssz_p5_full_pipeline.py --strict
```

During implementation, run only focused tests relevant to changed code.

Run the full suite at these milestones only:

```text
A. all regional direct coefficient paths ready
B. global direct 41 ready
C. direct KRGSM ready
D. coupled spectral solver converged
E. release candidate
```

Do not run the entire historical audit after every small edit.

---

# 29. No file questions unless an actual action-level input is irreducibly absent

You already have:

```text
old action cover
outer resolved raw 41
central corrected 41
inner selected 41 candidate
core selected 41
weak exterior 39/41
central odd profile
re-solved outer/inner backgrounds
central exact SVT data
core action-level G4XX/G5 data
analytic-center data
```

Do not ask Carmen/Lino for those files.

If a genuinely necessary action-level jet is not reconstructible from supplied sources, report exactly:

```text
region
equation requiring it
jet name
files checked
why current higher-jet/tubular machinery cannot reconstruct it
```

A generic "data missing" message is not sufficient.

---

# 30. No new scope

Do not add nonlinear evolution, LIGO analysis, cosmology, UV completion, new P5 variants or new action families.

The target is the defined linear finite-l Full-SVT production/reduction/spectral chain.

---

# 31. Operational meaning of Absolute Full Closure

For this repository, `ABSOLUTE_FULL_CLOSURE` means the defined implementation chain has no remaining gate:

```text
locked regional selected member
→ direct regional coefficients
→ global direct 41/41
→ common constraints
→ direct K,R,G,S,M
→ required Odd/Even finite-l gates
→ coupled spectral system
→ converged production spectral search
→ reproducible clean release
```

It does not mean a proof of all nonlinear/UV/global mathematical properties of the theory.

Do not add new closure criteria after completing the specified ledger.

---

# 32. Required gate ledger

Use a fixed ledger at minimum containing:

```text
PRODUCTION_REGIONS
GLOBAL_ACTION_COVER
WEAK_DIRECT_41
OUTER_DIRECT_41
CENTRAL_DIRECT_41
INNER_DIRECT_41
CORE_DIRECT_41
ANALYTIC_CENTER
INTERFACE_CONTINUITY
GLOBAL_DIRECT_41OF41
SHARED_BASELINE
CONSTRAINT_SYSTEM
CONSTRAINT_PIVOTS
PROFILE_REDUCER
DIRECT_KRGSM_L6
DIRECT_KRGSM_L12
DIRECT_KRGSM_L20
DIRECT_KRGSM_L42
DIRECT_KRGSM_L110
DIRECT_KRGSM_L420
DIRECT_KRGSM_L1000
HIGH_L_M
JET9D8_REGRESSION
PURE_SECTOR_REGRESSIONS
ODD_GLOBAL
EVEN_GLOBAL
COUPLED_RADIAL_OPERATOR
QNM_BOUNDARY_IMPLEMENTATION
QNM_NUMERICAL_CONVERGENCE
QNM_BRANCH_TRACKING
UNSTABLE_MODE_SEARCH
PYTEST
STRICT_PIPELINE
CLEAN_INSTALL
CLEAN_CHECKOUT_PIPELINE
MANIFEST
SHA256
ABSOLUTE_FULL_CLOSURE
```

A historical witness or prestaged reference may not satisfy a Direct gate by existence alone.

---

# 33. Exact forward execution order

Do not deviate unless a current focused gate requires a local repair:

```text
1. Verify this handoff once.
2. Generate current action cover from regions.py.
3. Complete weak direct 41.
4. Regenerate/verify outer direct 41 against the supplied outer reference.
5. Normalize/regenerate central exact SVT direct 41 against supplied reference.
6. Use/regenerate supplied central Odd profile.
7. Complete independent inner direct 41 from resolved inner member; reconstruct only necessary missing higher jets.
8. Implement/port the general Horndeski core emitter only if required for direct core provenance.
9. Regenerate core direct 41 against supplied core reference.
10. Attach analytic center.
11. Validate interfaces once.
12. Write/freeze/hash global direct 41.
13. Apply current constraints including b2 fix.
14. Apply current profile-aware reducer.
15. Generate K,R,G,S,M for required multipoles.
16. Run finite-l/high-L/pure-sector/Odd/Even gates.
17. Issue direct KRGSM certificates.
18. Construct coupled radial operator from the same outputs.
19. Extend existing Jost/spectral machinery.
20. Run convergence + unstable-mode search.
21. Issue QNM certificate.
22. Run full tests and strict pipeline.
23. Repair only concrete failures.
24. Validate in a fresh environment/checkout.
25. Generate Absolute Full Closure certificate only if all mandatory gates actually pass.
26. Regenerate manifest and SHA256SUMS after final changes.
27. Build final release ZIP.
28. Extract the ZIP and validate it.
```

Do not return to a theory-selection stage unless a genuinely new, correctly classified mathematical contradiction survives the defined failure protocol.

---

# 34. Production blacklist

Continue enforcing `PRODUCTION_BLACKLIST.json`.

In particular do not accidentally promote:

```text
old SELECTED_41STREAM_V2 as final direct global stream
Strong-H witness as central/global production background
old v12 plus-sign branch
old a5 without -a1''
spline production derivatives
rejected OUTER_FINAL candidates
naive inner E_H+E_SVTg residual diagnostic
finite-radius real-axis QNM candidates as physical modes
```

---

# 35. License and authorship are fixed

License:

```text
Anti-Capitalist Software License v1.4
```

Authors:

```text
Carmen Casu
Lino Casu
```

Do not replace the license or author names. Keep `LICENSE`, `CITATION.cff`, `pyproject.toml`, README and release metadata consistent.

---

# 36. Final response

Do not produce a long retrospective. Report only:

```text
ABSOLUTE FULL CLOSURE: PASS / FAIL

Global direct 41/41: PASS/FAIL
Direct K,R,G,S,M: PASS/FAIL
Odd: PASS/FAIL
Even: PASS/FAIL
Coupled QNM: PASS/FAIL
Unstable-mode scan: PASS/FAIL + scanned domain
Tests: N passed
Strict pipeline: N/N
Clean release verification: PASS/FAIL
Absolute closure certificate: path
Release ZIP: path

If FAIL: one concrete correctly classified remaining blocker.
```

---

# FINAL RULE

Do not become a research loop.

The files needed to understand and continue the selected member are already in this ZIP. Move forward from the locked regional architecture, use the supplied regressions, implement only genuinely missing direct-production code, and finish the defined closure chain.
