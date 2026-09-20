# SSZ P5 — Full Closure Research Repository

**Carmen Casu and Lino Casu · Anti-Capitalist Software License v1.4**

Research sources, action profiles, perturbation code and reproducible verification
for the regional P5 Full-SVT member: weak Horndeski exterior, outer H/SVT
handover, central genuine-SVT lobe, inner SVT/H handover, general Horndeski core,
and an analytic center chart.

**Latest imported checkpoint: September 18, full-action gate.** The updated
software suite passes **96 tests**. The 4D `(phi,X,F,Y)` Hessian completion and
local reachability audit are implemented; the remaining Inner construction needs
additional mixed `f3/f4` action-jet directions. See the
[current status](INTERMEDIATE_STATUS_2026-09-18_FULL_ACTION_GATE.md) and
[integration verification](docs/FULL_ACTION_SNAPSHOT_INTEGRATION_2026-09-18.md).

> **Current result: Absolute Full Closure is not certified.**
> The selected central Full-SVT coefficient member fails the necessary finite-L
> kinetic-positivity gate at L=6,12,20,42. At the interior point u≈0.690004,
> L=6 gives a kinetic eigenvalue ≈−43.8617. An independent Schur calculation,
> stencil/grid checks and basis rescalings reproduce the negative direction.
> This is a statement about this coefficient member, not a no-go theorem for
> P5 geometry. See the [reproducible execution report](docs/REGIONAL_KINETIC_GATE_2026-09-17.md)
> and [machine evidence](data/diagnostic/REGIONAL_CENTRAL_KINETIC_GATE.json).

## Reproduce

Use CPython **3.14** on Linux, from the repository root:

```bash
git clone https://github.com/error-wtf/SSZ_FULL_CLOSURE.git
cd SSZ_FULL_CLOSURE
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
pytest -q
python tools/regenerate_regional_coefficients.py
python tools/audit_regional_kinetic.py
python ssz_p5_full_pipeline.py --strict
```

The regional regeneration command checks the Outer SVT-sector replay and central
normalization. The kinetic audit and strict pipeline currently return **exit 2**:
the frozen central member fails kinetic positivity and the complete global/spectral
certificates are absent. Software tests verify implementations and rejection guards;
their passing does not establish physical closure.

For the separately labelled archive audit, run:

```bash
python ssz_p5_full_pipeline.py
```

The delivered handoff starts at [START_HERE_CODEX.md](START_HERE_CODEX.md), followed
by the [execution contract](CODEX_FINAL_EXECUTION_CONTRACT.md) and
[task DAG](CODEX_TASK_DAG.json). [tools/codex_start.sh](tools/codex_start.sh) performs
the start checks and terminates with the actual strict result. The handoff's
reported archive success must not be confused with a production-closure PASS.

## Production chain and verified scope

```text
Locked regional Full-SVT action / background jets
    → direct unreduced 41-slot regional emission
    → common assembly with one shared-baseline subtraction
    → same-basis K, R, G, S, M
    → finite-L / center / interface / pure-sector / high-L gates
    → coupled spectrum and convergence from that same operator
    → Absolute Full Closure certificate
```

| Layer | Current result | Scope |
| --- | --- | --- |
| Regional source integration | Implemented | Active member, source registry, region resolver and certificate binding |
| Outer raw SVT emission | Regression PASS | 41 slots, max scaled reference difference ≈9.03e-9; full H+SVT assembly is separate |
| Central selected normalization | **41/41 PASS** | 3,809 rows exported with source/output hashes; stored comparable ZK slots reproduce exactly. [Scoped certificate](data/generated/central/CENTRAL_DIRECT_41_CERTIFICATE.json) |
| Inner selected normalization | **41/41 PASS** | 1,400 finite normalized rows; normalization is not yet a single-action Direct-41 certificate |
| Inner action / Direct-41 closure | **OPEN / FAIL** | The old 3D `(phi,X,F)` Hessian exclusion is superseded by the `Y` sector, but the full 4D `f2(phi,X,F,Y)` audit shows that the current independently blended lower/principal target set is not reachable by a background-null `f2` Hessian alone. Full mixed `f3/f4` action-jet controls are now the next gate. |
| Central finite-L kinetic gate | **FAIL for frozen selected member** | Negative kinetic directions remain recorded for L=6,12,20,42; the Inner action work does not supersede this gate |
| Global direct 41 and KRGSM | Not certified | Remaining regional direct paths, center and interfaces are not claimed complete |
| Coupled spectral convergence | Not run | Same-operator global KRGSM and kinetic gates are prerequisites |
| Absolute closure / v1.0.0 | Not released | Strict production evidence remains incomplete |

## September 18 working-snapshot integration

The supplied working snapshot is preserved in
[archive/working_2026-09-18](archive/working_2026-09-18/IMPORT_PROVENANCE.json).
All 63 imported files retain their source hashes. Its historical success reports
remain archive results, not a new global-closure certificate.

The raw `f2FF/f2XF/f2XX` inverse now controls `v1/v4/c2` in the actual ZK
normalization, including rejection of unreachable zero-electric-field targets.
The general MH primitive emitter uses the existing JET9D8 service and reproduces
the imported executable bitwise on the same 3,299-row input. It still consumes
primitive profiles; this is not a complete raw G4/G5 action reconstruction.

```bash
python tools/audit_working_import.py
python tools/build_inner_targets.py
python tools/audit_inner_y_hessian.py
```

The Inner lower-order and principal target maps now have executable local response
inverses, but those split inverses are **not** a common single-action certificate.
The former 3D common-Hessian audit omitted the `Y` dependence of `f2` and is retained
only as historical diagnostic evidence. The corrected 4D audit uses the six transverse
`(XX,XF,XY,FF,FY,YY)` directions and shows a second, more useful limitation: on the
static electric background their local algebraic response to
`(v5,c3,v1,v4,c2)` has maximum rank 4, and the currently independent target blends
leave that column space. Therefore the next production construction must use the
**full action-jet control space**, including mixed `f3/f4` directions, rather than
coefficient-level target matching. See
[the machine-readable reachability audit](data/generated/inner_y_hessian/INNER_F2_4D_ALGEBRAIC_REACHABILITY.json) and
[the continuation contract](CONTINUE_IMPLEMENTATION_2026-09-18.md).

## Remaining production gates

The Inner construction and the Central kinetic failure are independent gates.
Completing Inner/Core does not by itself repair the negative Central directions.
The remaining work is:

1. **Common Inner action:** combine the `f2(phi,X,F,Y)` Hessian with mixed
   `f3/f4` jets in one holonomic action-control space, deriving all six channels
   `v5,c3,e3,v1,v4,c2` together.
2. **Inner Direct-41:** emit the full guard grid, close `a5/v7/v12`, then trim
   and verify Central–Inner and Inner–Core interfaces. Both `INNER_DIRECT_41`
   and `DIRECT_ACTION_REPLAY_COMPLETE` require genuine action replay.
3. **Core/G5:** complete the general quintic action emission and required
   derivative orders for all 41 slots.
4. **Exact center:** construct the regular Taylor chart at `r=0`, with finite
   curvature invariants and compatible coefficient/operator limits.
5. **Global Direct-41:** derive one exterior-to-center stream with common
   conventions and basis, one shared baseline, and smooth action/jet matching.
6. **Central kinetic gate:** test a complete action re-emission of the selected
   member. If the negative mode persists, a changed Central action member must
   pass background, holonomy, interface and stability checks. This is an open
   physical selection problem; software completion alone does not guarantee PASS.
7. **Global KRGSM:** jointly reduce that stream and check constraint pivots,
   kinetic positivity, radial/angular gradients, matrix symmetries, canonical
   `R=0`, sector limits, high-L behavior and interfaces.
8. **Coupled spectrum:** use exactly that operator for converged resonance
   branches and a documented search for unstable poles. Archived test-field
   or WKB values cannot certify this gate.
9. **Release:** reproduce all required certificates and the strict pipeline
   from a fresh extraction, with a verified manifest and SHA256 inventory.
   `ABSOLUTE_FULL_CLOSURE` and `v1.0.0` require strict exit code **0**.

## Selected member and numerical conventions

The active [regional member](SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json)
allows nonzero A0prime in its electric SVT regions. The historical
[epsilon-Y member](SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json) and strong-H carrier
remain scoped alternatives/control witnesses; they do not replace the regional
production action. In particular no global epsilon-Y-only slot modification or
global A0prime=0 condition is imposed on this chain.

- **41 slots:** `a1..a9`, `b1..b5`, `c1..c6`, `d1..d4`, `e1..e4`, `v1..v13`.
- **Derivative service:** JET9D8, window 9 and degree 8.
- **Locked identities:** `v12 = -v6/(2h)` and the corrected `a5` identity,
  including `-a1''`.
- **Assembly:** shared baseline once; eliminate the common constraints after
  unreduced assembly. Reduced sector matrices are not added together.
- **Reduced convention:** retain radial product rules, basis derivatives and
  integration-by-parts terms; verify `S.T = -S` and canonical `R = 0`.
- **Center:** analytic center treatment is required; sampled punctured-core
  rows alone do not certify the exact center.

Tolerances are centralized in [NUMERICAL_POLICY.json](NUMERICAL_POLICY.json).
No production tolerance has been relaxed to obtain a passing result.

## Strict verification and spectral scope

Strict mode checks file integrity, the software tests, existing audit results,
the necessary frozen-carrier on-shell identity, and direct-production closure evidence. It rejects absent direct matrices,
missing multipoles, inconsistent hashes, nonpositive finite-L results, matrices
that do not reproduce from the recorded coefficient stream, and spectral records
bound to a different operator. It also requires convergence records for every
reported branch and an explicitly bounded search in the upper frequency half-plane.
These checks validate the declared evidence; they do not create a spectral solution.

The repository preserves scalar/Maxwell WKB, light-ring/eikonal estimates and
Jost diagnostics under [data/qnm](data/qnm). The previously rejected real-axis
inward-shooting candidates remain diagnostic. No final coupled HSVT frequencies
are claimed in this edition.

## Find the research and code

| Location | Contents |
| --- | --- |
| [paper/](paper/) | Byte-original papers, monograph, LaTeX and figures |
| [src/ssz_p5/](src/ssz_p5/) | Package APIs, JET service, frozen deformation, validation and provenance |
| [src/](src/) | Migrated numerical emitters, constraint maps and profile reducer |
| [data/production/](data/production/) | Curated source/action profiles and historical local witnesses; classification is not a global direct certificate |
| [data/authoritative/](data/authoritative/) | Authoritative source and comparison artifacts |
| [data/regression/](data/regression/) | Reference results used for comparisons |
| [data/diagnostic/](data/diagnostic/) | Excluded streams and current diagnostic findings |
| [tests/](tests/) | Unit, integration, regression and rejection tests |
| [archive/](archive/) | Preserved working snapshots and the earlier import |
| [provenance/source_manifests/](provenance/source_manifests/) | Original supplied inventories |

[PRODUCTION_BLACKLIST.json](PRODUCTION_BLACKLIST.json) controls production
eligibility. A filename containing `FINAL`, `COMPLETE` or `SELECTED` is not a
certificate. The historical papers remain unchanged; current executable behavior
and verified implementation findings are stated separately.

## Integrity, citation and license

```bash
python tools/release_manifest.py --check
```

[MANIFEST.json](MANIFEST.json) and [SHA256SUMS](SHA256SUMS) cover the release
files. Runtime reports are excluded from the immutable inventory. Hash verification
establishes file identity, not scientific correctness.

Use [CITATION.cff](CITATION.cff) for citation metadata. The software license is
[Anti-Capitalist Software License v1.4](LICENSE).

License: Anti-Capitalist Software License v1.4

Authors: Carmen Casu and Lino Casu

Copyright (c) 2026 Carmen Casu and Lino Casu
