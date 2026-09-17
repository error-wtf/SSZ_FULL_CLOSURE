# SSZ P5 — Full Closure Research Repository

**Carmen Casu and Lino Casu · Anti-Capitalist Software License v1.4**

Research sources, action profiles, perturbation code and reproducible verification
for the regional P5 Full-SVT member: weak Horndeski exterior, outer H/SVT
handover, central genuine-SVT lobe, inner SVT/H handover, general Horndeski core,
and an analytic center chart.

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
| Inner selected normalization | **41/41 PASS; interface FAIL** | 1,400 finite rows; selected zero lower-order slots disagree with the frozen central endpoint. [Evidence](docs/INNER_PRODUCTION_EXPORT.md) |
| Central finite-L kinetic gate | **FAIL** | Negative kinetic directions for L=6,12,20,42 |
| Global direct 41 and KRGSM | Not certified | Remaining regional direct paths, center and interfaces are not claimed complete |
| Coupled spectral convergence | Not run | Necessary kinetic gate fails; archived test-field roots are not a substitute |
| Absolute closure / v1.0.0 | Not released | Strict mode rejects the failed gate and absent production evidence |

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
