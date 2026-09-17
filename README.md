# SSZ P5 — Full Closure Research Repository

**Carmen Casu and Lino Casu · Anti-Capitalist Software License v1.4**

Reproducible source, action profiles, perturbation tools and research records for
P5 geometry with a Horndeski–Maxwell carrier and the frozen background-null
scalar–vector–tensor deformation. Original papers and historical calculations
are preserved alongside executable verification code.

> **Current production status: absolute full closure is not certified.**
> The software tests and the archive audit pass. A sector-aware audit now assigns the disputed radii to the genuine-SVT
> production lobe, where `A0prime != 0`; the pure-H null-vector identity is
> therefore `NOT_APPLICABLE` to the global member. The strong-H file remains a
> principal/control witness. Direct global coefficients and coupled spectral
> products are still not certified. See the [region report](data/diagnostic/PRODUCTION_REGION_ASSIGNMENT.json).

## Reproduce

Use CPython **3.14** on Linux and run from the repository root:

```bash
git clone https://github.com/error-wtf/SSZ_FULL_CLOSURE.git
cd SSZ_FULL_CLOSURE
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
pytest -q
python ssz_p5_full_pipeline.py --strict
```

A strict exit code `2` currently reports the frozen-carrier on-shell identity
failure and the missing direct-production closure certificate. Installation instructions and interpretation of the results are in
[REPRODUCE.md](REPRODUCE.md).

To run the separately labelled historical archive and regression audit:

```bash
python ssz_p5_full_pipeline.py
```

Runtime reports are written to `FULL_PIPELINE_REPORT.json`,
`FULL_PIPELINE_REPORT.md`, `build/audit.json` and `build/gates.csv`.

## The production chain

```text
Frozen P5 action / background jets
    → one global unreduced 41-slot stream
    → epsilon_Y correction to v1 and v10
    → common profile-aware constraint elimination
    → K, R, G, S, M at L = 6, 12, 20, 42, 110, 420, 1000
    → finite-L / center / interface / pure-sector / high-L gates
    → coupled spectral calculation using those same matrices
    → convergence and bounded unstable-mode search
    → absolute-closure certificate
```

| Verification layer | Current result | Meaning |
| --- | --- | --- |
| Software tests | PASS | Includes negative certificate tests and direct constraint stationarity |
| Existing archive/audit checks | PASS | Reproduces the defined historical witnesses |
| Frozen carrier, necessary on-shell identity | **FAIL** | Five numerical routes give a nonzero interior residual; no action member has been changed |
| Direct global action → 41-slot regeneration | Not completed | General quartic/quintic core emission is not implemented by the existing restricted emitter |
| Direct global finite-L matrices | Not certified | A complete same-action set of generated matrices is required |
| Coupled spectral convergence | Not certified | Archived test-field calculations cannot substitute for the coupled operator |
| Absolute closure / v1.0.0 | Not released | Strict mode rejects missing or invalid production evidence |

The existing ZK emitter also uses a restricted Einstein/SVT action: its `f4XX`
terms are not Horndeski `G4XX`, and it contains no general `G5` emission.
Both active emitters reject unsupported nonzero Horndeski jets.

The sector assignment is machine readable in [PRODUCTION_REGION_ASSIGNMENT.json](data/diagnostic/PRODUCTION_REGION_ASSIGNMENT.json).
The concrete source-to-emitter comparison is machine readable in
[data/diagnostic/PRODUCTION_SOURCE_CONTRACT.json](data/diagnostic/PRODUCTION_SOURCE_CONTRACT.json).
[Implementation findings](docs/DIRECT_IMPLEMENTATION_FINDINGS.md) record the
verified constraint repair and the remaining emission boundary. These are
implementation findings. The separately documented background inconsistency is
not a computed perturbative instability or a no-go theorem for P5 geometry.

## Frozen model and numerical conventions

The [action member](SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json) fixes
`epsilon_Y = 0.01` and `A0prime = 0`. With `kappa = h * phi_r**2 = -2X`,
`ZA = 1 - 2 * epsilon_Y * kappa` multiplies only the Maxwell slots `v1` and
`v10`. The other 39 slots retain the Horndeski/Maxwell baseline values.

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
