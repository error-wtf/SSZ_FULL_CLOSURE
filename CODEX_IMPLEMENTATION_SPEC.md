# Codex Implementation Specification — SSZ P5 Full Closure Repository

**Release basis:** 2026-09-16  
**Purpose:** remove architectural ambiguity before Codex starts coding.  
**Primary target:** close `DIRECT_GLOBAL_KRGM_EXPORT` without changing the frozen P5 geometry or the production action member.

---

## 0. Read order for Codex

Do not infer architecture from filenames. Read exactly in this order:

1. `SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json`
2. `STATUS_FULL_CLOSURE.md`
3. `CODEX_IMPLEMENTATION_SPEC.md` (this file)
4. `NUMERICAL_POLICY.json`
5. `CODEX_TASK_DAG.json`
6. `README_FOR_CODEX.md`
7. only then the monograph and source modules.

The implementation must preserve the scientific status; it must not redesign the model.

---

## 1. Frozen scientific object

The production model is

\[
S_{\rm prod}=S_{\rm H,P5}+\epsilon_Y\int d^4x\sqrt{-g}\,Y,
\qquad \epsilon_Y=10^{-2},
\]

with

\[
Y=\nabla_\mu\phi\nabla_\nu\phi F^{\mu\alpha}F^\nu{}_{\alpha},
\qquad A_0'(r)=0.
\]

### 1.1 Consequences that are implementation facts, not TODOs

On the zero-vector background:

- `Delta f2 = epsilon_Y * Y` is exactly background-null;
- `Delta E00 = Delta E11 = Delta Ephi = Delta EA = 0`;
- the scalar/metric Horndeski quadratic block is unchanged;
- the genuine-SVT deformation is quadratically nontrivial in the vector sector;
- for constant `epsilon_Y` and `A0prime=0`, its 41-slot effect is **only** a multiplicative change of the Maxwell vector slots `v1` and `v10`.

Define

\[
\kappa=h\phi_r^2=-2X,\qquad
Z_A=1-2\epsilon_Y\kappa.
\]

For the accepted Zhang--Kase conventions and the Maxwell baseline (`f2F=1`),

\[
v_1^{\rm prod}=\frac{r^2}{2}\sqrt{\frac{h}{f}}\,Z_A,
\qquad
v_{10}^{\rm prod}=-\frac{\sqrt{fh}}{2}\,Z_A.
\]

Equivalently, if a regenerated pure-Maxwell/Horndeski 41 stream already has baseline `v1_H` and `v10_H`, apply

```python
kappa = h * phi_r**2
ZA = 1.0 - 2.0 * epsilon_y * kappa
v1  = v1_H  * ZA
v10 = v10_H * ZA
```

on the zero-vector production branch.

For this production member, do **not** import the historical electric-SVT `v2...v13` structure. In particular, for the epsilon-Y deformation at `A0prime=0`, the genuine-SVT delta to all other 41 slots is zero. After complete action-level assembly, enforce the usual canonical identities from the Horndeski/Maxwell baseline (`v2=v3=v4=v5=v6=v7=v8=v11=v12=v13=0` in pure zero-vector Maxwell-Horndeski patches unless a Horndeski emitter explicitly provides otherwise; `v9` remains its baseline coefficient).

The certified global margin is

\[
\min Z_A=0.9872911806592305>0.
\]

This simplification is the preferred production path. The historical electric central-SVT representative is a **regression oracle only**, not the final production action.

---

## 2. No-redesign rules

Codex must not spend tokens exploring alternatives to these decisions:

1. **Geometry:** frozen P5 only.
2. **Scalar patch coordinate:** `phi = Xi` where monotone.
3. **Production vector background:** `A0prime = 0`.
4. **Production genuine-SVT term:** exactly `epsilon_Y * Y`, `epsilon_Y=0.01`.
5. **Derivative service:** JET9D8, window 9, degree 8 unless an analytic derivative is supplied.
6. **v12 branch:** `v12 = -v6/(2*h)`.
7. **a5:** `a5 = a2' - a1'' - (A0' v4/2)' + A0' v5/2`.
8. **Auxiliary vector canonicalization:** after unreduced action addition, `v7=v2**2/(4*v1)`.
9. **Assembly order:** combine unreduced action coefficients first, eliminate common constraints second.
10. **Reduced convention:** retain antisymmetric `S`; final canonical convention has `R=0`.
11. **Center:** analytic Taylor chart controls the exact `r=0` sign; do not divide by `kappa` there.
12. **QNM:** disabled until direct-global-KRGM certificate passes.

---

## 3. Target repository tree

Create exactly this logical structure (minor naming changes require no scientific redesign, but avoid them unless necessary):

```text
ssz-p5/
├── pyproject.toml
├── README.md
├── CITATION.cff
├── LICENSE
├── src/ssz_p5/
│   ├── __init__.py
│   ├── cli.py
│   ├── config.py
│   ├── types.py
│   ├── geometry/
│   │   ├── __init__.py
│   │   ├── p5.py
│   │   ├── center.py
│   │   └── light_rings.py
│   ├── action/
│   │   ├── __init__.py
│   │   ├── model.py
│   │   ├── horndeski.py
│   │   ├── eps_y.py
│   │   ├── partitions.py
│   │   └── patches.py
│   ├── jets/
│   │   ├── __init__.py
│   │   └── jet9d8.py
│   ├── coefficients/
│   │   ├── __init__.py
│   │   ├── schema41.py
│   │   ├── mh.py
│   │   ├── zk.py
│   │   └── assemble.py
│   ├── constraints/
│   │   ├── __init__.py
│   │   └── even.py
│   ├── reducer/
│   │   ├── __init__.py
│   │   ├── differential_operator.py
│   │   └── canonical.py
│   ├── stability/
│   │   ├── __init__.py
│   │   ├── pivots.py
│   │   ├── finite_l.py
│   │   ├── exact_oracles.py
│   │   └── interfaces.py
│   ├── export/
│   │   ├── __init__.py
│   │   ├── coefficients.py
│   │   ├── krgm.py
│   │   └── certificate.py
│   ├── provenance/
│   │   ├── __init__.py
│   │   ├── manifest.py
│   │   └── records.py
│   └── qnm/
│       ├── __init__.py
│       ├── gate.py
│       ├── asymptotics.py
│       └── spectral.py
├── data/
│   ├── authoritative/
│   ├── regression/
│   ├── generated/
│   └── certificates/
├── tests/
│   ├── unit/
│   ├── regression/
│   ├── negative/
│   └── integration/
├── scripts/
└── .github/workflows/
```

Do not keep production code dependent on `/mnt/data` paths. The old source modules may be imported initially only as migration references; production package paths must be relative/configurable.

---

## 4. Required dataclasses and types

Implement these public types. They are API contracts.

```python
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
import numpy as np

Array = np.ndarray

@dataclass(frozen=True)
class P5Background:
    r: Array
    u: Array
    phi: Array
    f: Array
    h: Array
    phi_r: Array
    X: Array
    A0prime: Array
    region: Array | None = None

@dataclass(frozen=True)
class ActionMember:
    release: str
    epsilon_y: float
    vector_branch: Literal["A0prime=0"]
    source_definition: Path

@dataclass(frozen=True)
class Coefficients41:
    background: P5Background
    slots: dict[str, Array]
    provenance_id: str

@dataclass(frozen=True)
class ConstraintPivots:
    Dh1: Array
    DeltaV: Array
    pivotA0: Array

@dataclass(frozen=True)
class ReducedOperator:
    L: int
    r: Array
    K: Array  # shape (N,3,3)
    R: Array  # shape (N,3,3)
    G: Array  # shape (N,3,3)
    S: Array  # shape (N,3,3)
    M: Array  # shape (N,3,3)
    pivots: ConstraintPivots
    provenance_id: str

@dataclass(frozen=True)
class GateResult:
    name: str
    status: Literal["PASS","FAIL","WARN","OPEN","PASS_EXACT","PASS_NUMERICAL"]
    value: object
    criterion: str
    evidence: str

@dataclass(frozen=True)
class DirectKRGMCertificate:
    pass_: bool
    release: str
    action_sha256: str
    coefficient_stream_sha256: str
    L_values: tuple[int, ...]
    gates: tuple[GateResult, ...]
```

### 4.1 Shape invariants

- every background scalar array: `(N,)`;
- every reduced matrix: `(N,3,3)`;
- all arrays in one object use one identical radial ordering;
- production radial ordering must be documented once and never silently reversed;
- before export, `r` must be strictly monotone and unique;
- no implicit pandas index is scientifically meaningful.

---

## 5. Public function contracts

### geometry

```python
def load_frozen_p5(source: Path) -> P5Background: ...
def validate_background(bg: P5Background) -> list[GateResult]: ...
def kappa(bg: P5Background) -> Array: ...
def locate_light_rings(bg: P5Background) -> dict[str, float]: ...
```

`kappa(bg)` must return exactly `bg.h * bg.phi_r**2`; test against `-2*X` where `X` is available.

### jets

```python
def derivative(
    x: Array,
    y: Array,
    order: int = 1,
    *,
    window: int = 9,
    degree: int = 8,
) -> Array: ...
```

Production implementation should be migrated directly from `ssz_p5_higher_jet_closure_2026-09-16.py`, not reinvented with splines.

### action

```python
def load_action_member(path: Path) -> ActionMember: ...
def epsilon_y_margin(bg: P5Background, epsilon_y: float) -> Array: ...
def apply_epsilon_y_deformation(
    coeffs_h: Coefficients41,
    epsilon_y: float,
) -> Coefficients41: ...
```

`apply_epsilon_y_deformation` is intentionally simple for the frozen member:

```python
kap = coeffs_h.background.h * coeffs_h.background.phi_r**2
ZA = 1 - 2 * epsilon_y * kap
out = deep_copy(coeffs_h)
out.slots["v1"]  = coeffs_h.slots["v1"]  * ZA
out.slots["v10"] = coeffs_h.slots["v10"] * ZA
# all other slots unchanged by epsilon-Y on A0prime=0
```

Before accepting this shortcut, assert the baseline is the zero-vector Maxwell/Horndeski branch (`A0prime==0` to tolerance) and run the analytic formula regression against the Zhang--Kase emitter at representative radii.

### coefficients

```python
SLOT_NAMES: tuple[str, ...]  # 41 names, canonical order

def emit_mh_patch(action_patch, bg: P5Background) -> Coefficients41: ...
def emit_zk_regression_patch(action_patch, bg: P5Background) -> Coefficients41: ...
def validate_41_schema(c: Coefficients41) -> list[GateResult]: ...
def canonicalize_aux_vector(c: Coefficients41) -> Coefficients41: ...
```

`emit_zk_regression_patch` is for historical/electric-SVT regression. It is not the production final-member constructor.

### constraints

```python
def constraint_maps(c: Coefficients41, L: int) -> dict[str, Array]: ...
def constraint_pivots(c: Coefficients41, L: int) -> ConstraintPivots: ...
def validate_constraint_maps(c: Coefficients41, L: int) -> list[GateResult]: ...
```

Port the accepted formulas from `ssz_hybrid_full_constraint_maps_JET9D8(1).py` literally before refactoring algebra.

### reducer

```python
def reduce_profile(c: Coefficients41, L: int) -> ReducedOperator: ...
def validate_operator(op: ReducedOperator) -> list[GateResult]: ...
```

Port `ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py` first as a golden implementation. Refactor only after byte/numerical regression tests pass.

### export

```python
def build_global_coefficients(config) -> Coefficients41: ...
def build_global_krgm(config, L_values: tuple[int,...]) -> dict[int, ReducedOperator]: ...
def write_direct_certificate(...)->Path: ...
```

### qnm

```python
def require_direct_krgm_certificate(path: Path) -> DirectKRGMCertificate: ...
```

Every public QNM entry point must call this guard. No bypass flag in production CLI.

---

## 6. Direct-global-KRGM build algorithm

This is the main implementation task. Follow this sequence exactly.

### Phase A — load immutable definitions

1. Load `SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json`.
2. Verify SHA-256 against the release manifest.
3. Verify `epsilon_Y == 0.01` and vector branch is `A0prime=0`.
4. Load P5/action patch authoritative data.

### Phase B — reconstruct the pure Horndeski/Maxwell 41 stream

The final production baseline is **not** the historical electric-SVT global stream.

Construct center-to-infinity pure Horndeski/Maxwell coefficients from the accepted Horndeski action patches:

- weak exterior;
- strong-field carrier through both light rings to `u=0.71`;
- punctured quartic core;
- analytic exact-center chart plus center-to-punctured handover.

Use archived selected 41 tables only as **regression targets** while migrating. The final generated stream must record its action/profile source for every row.

At each patch:

1. regenerate action-dependent slots using the accepted emitter or exact patch formulas;
2. regenerate lower-order derivatives with JET9D8;
3. enforce accepted `a5`;
4. enforce `v12=-v6/(2h)` (trivial where `v6=0`);
5. canonicalize common vector auxiliary completion only after complete patch assembly;
6. verify pure Maxwell-Horndeski exact invariant conditions.

### Phase C — patch interfaces

Do not concatenate by nearest row.

For every interface:

1. work in monotone `phi` when the handover is a real action handover;
2. use the archived flat `C^infinity` patch definition/certificate;
3. re-solve the background-constrained lower jets if the patch construction requires it;
4. regenerate the 41 slots **after** the action patch is fixed;
5. test left/right values and derivative continuity at the order needed by the reducer;
6. include at least four JET9D8 stencil widths around interfaces in a convergence test.

The exact center is a separate analytic chart. Do not ask a finite-difference grid to certify the sign of an `O(r^4)` quantity at `r=0`.

### Phase D — add the production epsilon-Y deformation

Once the pure-H global 41 stream is valid:

```python
kap = h * phi_r**2
ZA = 1 - 2*epsilon_y*kap
assert ZA.min() > 0
v1 *= ZA
v10 *= ZA
```

All other 41 slots remain those of the zero-vector Horndeski/Maxwell baseline for this frozen production member.

Then re-run structural identities and auxiliary canonicalization checks. Do not import historical electric-SVT slot deltas.

### Phase E — reduce

For each `L` in `(6,12,20,42,110,420,1000)`:

1. compute generalized-psi constraint maps;
2. assert all pivots are safely nonzero;
3. build the complete unreduced common action;
4. perform profile-aware `T^dagger P T` reduction retaining all radial product rules;
5. extract canonical `K,R,G,S,M` using the accepted convention;
6. assert no differential-operator key above second order remains;
7. save matrices with one shared radial coordinate and metadata.

### Phase F — certify

The certificate can pass only if every required gate in section 9 passes. Write the certificate last. Never create a provisional `pass:true` file.

---

## 7. Production exports and schemas

### 7.1 Global 41 coefficient export

Preferred file: `data/generated/SSZ_P5_GLOBAL_41_COEFFICIENTS.parquet` (CSV mirror optional).

Required columns:

```text
r,u,phi,f,h,phi_r,X,A0prime,region,patch_id,
[a1..a9],[b1..b5],[c1..c6],[d1..d4],[e1..e4],[v1..v13],
source_action_id,source_sha256,generator_version
```

Never use filename semantics as certification.

### 7.2 KRGM export

One file per L or one HDF5/Zarr container is acceptable. Preferred simple layout:

```text
data/generated/krgm/L0006.npz
...
data/generated/krgm/L1000.npz
```

Each contains:

```text
r,u,K,R,G,S,M,Dh1,DeltaV,pivotA0
```

where each matrix is `(N,3,3)`.

### 7.3 Metadata

Every generated artifact must have a sibling JSON provenance record matching `schemas/artifact_provenance.schema.json`.

---

## 8. Numerical policy

All shared thresholds live in `NUMERICAL_POLICY.json`; production code must import them from one place.

Rules:

- do not scatter magic tolerances across tests;
- distinguish exact algebraic gates, numerical regression gates, and implementation gates;
- exact analytic center certificates are not downgraded because finite differences lose precision;
- relative comparisons use a scale floor, not division by tiny target values;
- interface tests report both absolute and scaled jumps;
- matrix positivity uses the symmetrized matrix `(A+A.T)/2` only after separately testing that antisymmetric contamination is within tolerance;
- generalized radial eigenvalues are computed via a symmetric whitening of positive `K`, not raw `eig(inv(K)@G)`.

Recommended helper:

```python
def scaled_error(a,b, floor=1.0):
    return abs(a-b) / maximum(floor, abs(b))
```

---

## 9. Required direct-export gates

A direct certificate must contain and pass at least these gates for every production L:

### Input/provenance
- action JSON hash matches manifest;
- all authoritative input hashes match;
- no source path is under `legacy_superseded/`;
- no generated row is tagged `CANDIDATE` or `SUPERSEDED`.

### Grid/background
- finite arrays only;
- `f>0`, `h>0`;
- radial grid monotone and duplicate-free;
- `kappa = h*phi_r**2 = -2X` within numerical tolerance where both are available;
- `A0prime=0` on production member within tolerance.

### 41 slots
- all 41 present and finite;
- accepted `a5` identity;
- accepted v12 branch;
- auxiliary H0/vector canonical identity;
- epsilon-Y only modifies `v1,v10` relative to regenerated zero-vector baseline;
- `min(ZA)>0`.

### Constraints
- `min(abs(Dh1)) > pivot_floor`;
- `min(abs(DeltaV)) > pivot_floor`;
- `min(abs(2*L*v9)) > pivot_floor`;
- exact cancellation channels for `delta_phi''` and `h1'` satisfy tolerance.

### Operator structure
- no operator derivative order >2;
- `K` symmetric within tolerance;
- `G` symmetric within tolerance;
- `S+S.T` within tolerance;
- `R` zero within tolerance;
- no NaN/inf in `K,R,G,S,M`.

### Stability
- every eigenvalue of symmetrized `K` strictly positive with configured margin;
- every generalized radial eigenvalue of `(G,K)` strictly positive with configured margin;
- production vector `ZA` strictly positive;
- exact Maxwell-Horndeski invariant oracle passes on H patches;
- analytic center no-ghost coefficient positive;
- historical ZK pure-SVT regression oracle still passes (regression-only, not production path).

### High-L
- mass identities `M11,M12,M13,M22,M23,M33` converge to accepted direct-action formulas;
- corrected `m5_minus` and distinct `m1_plus/m1_minus` routes are tested;
- large-L radial/operator behavior matches archived regression tolerances.

### Interfaces
- no discontinuous production patch jump;
- derivative-sensitive interface quantities converge under JET stencil variation;
- no use of old `SELECTED_41STREAM_V2` as production source.

---

## 10. Negative tests Codex must implement

Tests must prove that the repository rejects known past mistakes.

1. Flip `v12` to plus sign -> test must fail.
2. Remove `-a1''` from `a5` -> finite-l regression must fail.
3. Call reducer separately on H and SVT then add results -> test must fail against combine-before-reduce reference.
4. Load a file from `legacy_superseded/` as production -> provenance guard must fail.
5. Introduce NaN in one slot -> schema/audit must fail.
6. Force one constraint pivot to zero -> reducer must refuse to run.
7. Set epsilon-Y above positivity boundary so `ZA<=0` -> stability test must fail.
8. Try to call QNM without direct certificate -> hard exception.
9. Forge `DIRECT_GLOBAL_KRGM_CERTIFICATE.json` with `pass:true` but wrong hashes -> master audit must fail once repo implementation is complete.
10. Use historical electric-SVT `SELECTED_41STREAM_V2` as final action -> source-policy test must fail.

---

## 11. Exact tests to create

### Unit

```text
test_kappa_identity.py
test_jet9d8_polynomials.py
test_a5_identity.py
test_v12_branch.py
test_eps_y_exact_delta.py
test_aux_vector_canonicalization.py
test_constraint_map_algebra.py
test_canonical_KGSM_extraction.py
```

### Regression

```text
test_mh_emitter_strongH.py
test_zk_emitter_archive.py
test_strongH_finite_l.py
test_highL_mass.py
test_center_taylor_certificate.py
test_interface_stencil_convergence.py
```

### Integration

```text
test_build_global_41.py
test_build_global_krgm.py
test_direct_certificate.py
test_master_auditor_require_direct.py
```

### Negative

Use the ten tests from section 10.

---

## 12. CLI contract

Use one console command, suggested name `ssz-p5`.

```text
ssz-p5 audit [--quick|--full] [--require-direct]
ssz-p5 build coefficients --output ...
ssz-p5 build krgm --L 6 12 20 42 110 420 1000 --output-dir ...
ssz-p5 export certificate --output ...
ssz-p5 plot geometry|stability|operator ...
ssz-p5 qnm ...
```

Rules:

- `qnm` exits nonzero unless a valid direct certificate is supplied;
- `build krgm` never consumes `legacy_superseded` paths;
- `export certificate` reruns required validations; it does not trust cached booleans;
- every CLI command has `--json-log` and records input hashes.

---

## 13. CI contract

### Pull requests

Run:

1. formatting/lint/type checks;
2. unit tests;
3. negative tests;
4. quick constructive auditor;
5. manifest integrity.

### Main branch/nightly

Additionally run:

1. full JET9D8 convergence regressions;
2. all default L direct builds once implemented;
3. global operator stability gates;
4. deterministic artifact hash comparison where appropriate.

### Release tag

Require:

```bash
pytest -q
ssz-p5 audit --full
ssz-p5 build coefficients
ssz-p5 build krgm --L 6 12 20 42 110 420 1000
ssz-p5 export certificate
python ssz_p5_full_closure_auditor.py --data-dir . --require-direct-krgm
```

Only after all return zero may a release tag advertise direct-global-KRGM completion.

---

## 14. QNM implementation boundary

Do not spend tokens implementing a full QNM solver until the direct certificate is green.

After it is green, consume **only** regenerated `ReducedOperator` objects. Do not rebuild perturbation matrices inside the QNM package.

Initial supported methods, in priority order:

1. exterior complex scaling + spectral collocation;
2. compactified Jost/spectral determinant;
3. continued fraction if a stable recurrence is derivable.

Reject real-axis inward shooting as the production method; it was already found susceptible to exponentially amplified incoming contamination.

QNM acceptance must include convergence with radial resolution, domain/compactification, scaling angle (for ECS), basis normalization, and root continuation.

---

## 15. Provenance policy

Every generated file must state:

- release identifier;
- git commit;
- generator module/function;
- Python version;
- dependency lock hash;
- action definition hash;
- input file hashes;
- parameters (`L`, JET window/degree, epsilon-Y, tolerance set);
- output hash;
- status (`authoritative`, `regression`, `diagnostic`, `superseded`).

A status is metadata, never inferred from a filename.

---

## 16. Migration map from handoff scripts

Use these as golden source, then refactor:

| handoff module | target package |
|---|---|
| `ssz_p5_higher_jet_closure_2026-09-16.py` | `jets/jet9d8.py` |
| `ssz_p5_mh_luminal_g4phi_emitter_JET9D8_2026-09-16.py` | `coefficients/mh.py` |
| `ssz_p5_zk_appendixA_emitter_JET9D8_2026-09-16.py` | `coefficients/zk.py` |
| `ssz_hybrid_full_constraint_maps_JET9D8(1).py` | `constraints/even.py` |
| `ssz_hybrid_unreduced_even_kernel.py` | `reducer/differential_operator.py` or `coefficients/assemble.py` |
| `ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py` | `reducer/canonical.py` |
| `ssz_p5_holonomic_a5_closure_2026-09-16.py` | `coefficients/lower_order.py` |
| `ssz_p5_F2b_action_derived_highL_mass_closure_2026-09-15.py` | `stability/high_l.py` |
| `ssz_p5_KGSM_highL_crosscheck_2026-09-16.py` | regression tests |
| `ssz_p5_full_closure_auditor.py` | keep top-level compatibility wrapper; internally call package audit API |

First obtain regression equivalence. Only then simplify code.

---

## 17. Files that are regression-only or forbidden as production truth

### Forbidden production sources

- any `OUTER_FINAL_SAME_ACTION*` Sep-16 variable-G4 transplant candidate;
- any plus-sign `v12` data;
- any old `a5` completion without `-a1''`;
- experimental A2/v6 ODE repair material;
- naive inner `E_H + E_SVTg` sum as full residual;
- deep-core raw direct-K result from rounded historical 41 slots;
- `ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv` as the final production member (it contains the historical electric-SVT patch architecture and an explicitly candidate inner segment).

### Allowed regression oracles

- historical Zhang--Kase electric-SVT exact lobe: validates the ZK emitter and pure-SVT limits;
- old global principal K/R/G export: validates qualitative/principal limits only;
- selected Horndeski 41 tables: migration regression targets, not substitutes for action regeneration;
- previous handover tables: background/action patch certificates and regression targets.

---

## 18. Definition of Done per phase

### D0 — repository boots
- editable install works;
- `ssz-p5 --help` works;
- legacy constructive auditor still passes.

### D1 — golden modules migrated
- JET9D8, emitters, constraint maps, reducer numerically match handoff versions;
- no `/mnt/data` hardcoding remains.

### D2 — pure-H global 41 regenerated
- center-to-infinity pure-H/Maxwell stream generated from frozen patch definitions;
- exact invariant oracles pass;
- interfaces pass continuity/convergence tests.

### D3 — production HSVT 41 generated
- epsilon-Y deformation applied exactly;
- only `v1,v10` change from the pure-H zero-vector baseline;
- `ZA_min` reproduces release value within tolerance.

### D4 — global KRGM generated
- all default L operators exist;
- structure, pivots and positivity pass.

### D5 — direct certificate
- certificate schema validates;
- all hashes match;
- master auditor with `--require-direct-krgm` exits 0.

### D6 — QNM enabled
- only after D5;
- convergence suite exists before any published frequencies.

---

## 19. What Codex should not spend time on

Do not:

- search for another P5 ansatz;
- tune epsilon-Y unless a regression exposes an actual bug;
- revisit the old central electric-SVT action as the production member;
- derive a new reducer architecture;
- replace JET9D8 by splines/automatic smoothing;
- redesign the physical basis `(psi, dphi, V)`;
- “fix” a failed test by weakening tolerances without a documented numerical convergence study;
- implement QNM early;
- infer authority from words like `FINAL` in old filenames.

---

## 20. First 12 Codex commits/tasks

Codex should implement in this order:

1. package skeleton + config + immutable model loader;
2. JET9D8 port + polynomial unit tests;
3. 41-slot schema + dataclasses + provenance types;
4. MH emitter port + regression;
5. ZK emitter port + regression oracle;
6. constraint-map port + pivot tests;
7. profile reducer port + golden regression;
8. pure-H patch loader/regenerator + interface tests;
9. epsilon-Y exact deformation + analytic tests;
10. direct global coefficient exporter;
11. direct global KRGM builder + certificate;
12. CI + release workflow; only then QNM branch.

If a task fails, fix that task before starting the next. Do not parallelize scientifically dependent phases simply to reduce wall-clock time.

---

## 21. Final invariant

The repository is successful only if a clean checkout can reproduce the scientific claim without chat history:

```text
immutable action JSON
        ↓
frozen P5/action patches
        ↓
action-level 41-slot regeneration
        ↓
epsilon-Y exact deformation
        ↓
common constraints
        ↓
profile-aware K,R,G,S,M
        ↓
strict gates + hashes
        ↓
DIRECT_GLOBAL_KRGM_EXPORT = PASS
        ↓
QNM allowed
```

Anything that bypasses an arrow is not the production pipeline.
