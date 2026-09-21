# SSZ P5 — Full Closure Research Repository

Current Angular checkpoint: [Eq47 projection and dense action audit](CURRENT_WORKING_STATE_2026-09-21_EQ47.md).
The projection hypothesis and genuine-SVT finite-L/Laurent checks pass. Dense
angular positivity fails; the mixed-branch formal leading order remains unresolved.
Continuation beyond u≈0.708 and global QNM are held. Absolute Full Closure is false.
See [machine status](FULL_CLOSURE_WORKING_STATUS_2026-09-21_EQ47.json).

**Carmen Casu and Lino Casu · Anti-Capitalist Software License v1.4**

Research sources, action profiles, perturbation code and reproducible verification
for the regional P5 Full-SVT member: weak Horndeski exterior, outer H/SVT
handover, central genuine-SVT lobe, inner SVT/H handover, general Horndeski core,
and an analytic center chart.

> **Current result: Absolute Full Closure is not certified.**
> Two concrete global members have now been rejected without discarding the
> reusable P5/perturbation infrastructure. The 2026-09-17 regional electric
> member fails Central finite-L kinetic positivity and common-action lower replay.
> The 2026-09-18 global `A0prime=0` epsilon-Y member passes many component tests,
> but a denominator-free Maxwell-Horndeski Eq.85 audit excludes it at the inner
> P5 light ring: with positive tensor H/a4 it requires negative tensor F. The
> active path is therefore an **on-shell electric Horndeski/SVT hybrid**, seeded
> by the machine light-ring audit in
> `data/generated/absolute_attempt_2026-09-19/ZERO_VECTOR_LIGHT_RING_AUDIT.json`.

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

## Repository observability

The repository can print its current truth state directly from machine-readable evidence:

```bash
python tools/show_repo.py
python tools/show_gate_matrix.py
python tools/show_member_matrix.py
python tools/show_evidence_index.py
```

`REPO_EVIDENCE_REGISTRY.json` classifies current audits and intentional rejection
witnesses, including their expected exit codes. Run all registered evidence with:

```bash
python tools/run_all_evidence.py
# add --full to include pytest and the release-manifest check
```

This avoids conflating an expected failing audit of a rejected member with a failure
of the active search, and avoids treating historical PASS reports as current Absolute
Closure. See [REPO_OBSERVABILITY.md](REPO_OBSERVABILITY.md).

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
| P5 geometry / 41-slot / JET9D8 infrastructure | **PASS / retained** | Independent of the final action-member decision |
| 2026-09-17 regional electric member | **Rejected for Absolute Closure** | Central finite-L kinetic FAIL at L=6,12,20,42; selected c3/e3 not a common direct-action replay |
| 2026-09-18 zero-vector epsilon-Y member | **Rejected as global on-shell member** | Component-level finite-L/descriptor PASS evidence retained; inner light-ring Eq.85 forces tensor F/H≈-1.63868 for A0prime=0 |
| Generalized-psi descriptor | **PASS** | H0/h1 constraints retained explicitly; validated on carrier and punctured Core |
| Descriptor pullback equivalence | **PASS** | Same differential operator as established reduced KRGSM on Strong-H chart for L=6,12,20,42,110,420,1000 |
| Singular radial DAE pencil | **PASS_STRUCTURE** | No explicit deep-Core D_h1 inversion |
| Electric hybrid Eq.85 seed | **GENERATED** | Inner light ring needs q=A0prime^2*v8_background_identity≈1.56365 if F=H is used as a local diagnostic target |
| Electric-hybrid Horndeski principal feasibility | **PARTIAL PASS / SEARCH** | Self-contained common-emitter recipe gives L6≈-3.24295, L12≈-0.130994, L20≈-0.0128774; L=42,110,420,1000 positive; constraint pivots remain large. Not background/on-shell certified |
| Active electric hybrid action member | **SEARCH OPEN** | Must solve background/on-shell identities and finite-L K simultaneously from one holonomic action |
| Global same-operator KRGSM | **OPEN** | Requires promoted electric hybrid Direct-41 member |
| Coupled spectral convergence | **BLOCKED** | QNM gate remains closed until same-operator global certificate exists |
| Absolute closure / v1.0.0 | **Not released** | No failed/open gate is relabeled PASS |

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
`data/generated/inner_y_hessian/INNER_F2_4D_ALGEBRAIC_REACHABILITY.json` and
[the continuation contract](CONTINUE_IMPLEMENTATION_2026-09-19.md).

## Selected member and numerical conventions

There is currently **no promoted Absolute-Closure production member**. The
2026-09-17 regional electric member and the 2026-09-18 global zero-vector epsilon-Y
member are retained as rejected research witnesses with useful scoped evidence.
The active search definition is
[SSZ_P5_ELECTRIC_HYBRID_SEARCH_CANDIDATE_2026-09-19.json](SSZ_P5_ELECTRIC_HYBRID_SEARCH_CANDIDATE_2026-09-19.json).
A future promotion must use one holonomic action with nonzero electric support
through the inner light-ring region and must pass all background, Direct-41 and
finite-L gates before QNM is reopened.

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
