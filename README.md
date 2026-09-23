# SSZ P5 — Full Closure

**Can a "segmented spacetime" black-hole geometry — one with an inner
*stable* light ring — be derived from a single, healthy law of gravity?**

*Research repository · Carmen Casu and Lino Casu · Anti-Capitalist Software License v1.4*

---

## The idea in one paragraph

Segmented Spacetime (SSZ) proposes a discrete micro-structure of space and
time. Its "P5" solution predicts a compact object with an unusual feature:
besides the familiar unstable photon orbit it possesses an **inner stable
light ring** — a place where light can circle forever. Structures like that
leave dramatic observational fingerprints (long-lived "echoes" in
gravitational-wave signals), but they are also famously fragile: most
theories that allow them are haunted by ghosts or explosive instabilities.

This repository asks the sharpest possible version of that question:

> Does there exist **one** covariant action — of the general
> scalar–vector–tensor class that contains Maxwell, Horndeski and their
> genuine scalar–vector couplings — that reproduces the P5 background
> **and** passes every stability and consistency test, with no parameter
> ever tuned to make a test pass?

"Full Closure" means exactly that: one action member (one cryptographic
hash) carried through the *entire* chain —

```
background equations → action integrability → light-ring compatibility
→ quadratic perturbation action → constraint elimination
→ no ghosts (K > 0) → radial speeds c_r² > 0 → angular speeds c_Ω² > 0
→ finite-l sectors → odd & vector sectors → interfaces → global regularity
→ quasinormal modes (last)
```

— with every step independently tested and machine-documented.

## Status

| | |
|---|---|
| **Verdict** | `ABSOLUTE_FULL_CLOSURE_PASS = False` — derived programmatically, never asserted |
| **Gates** | 10 of 28 required gates PASS · 18 open · live from `GATE_STATUS.json` |
| **Tests** | 162/162 passing |
| **Head** | main (pushed; live commit via `git rev-parse HEAD`) |

The verdict is computed, not claimed:

```bash
PYTHONPATH=src python -c "
from ssz_p5.closure.gates import full_closure_verdict
print(full_closure_verdict()['ABSOLUTE_FULL_CLOSURE_PASS'])"
```

## Test results

**Software suite — 162/162 PASS** (commit `8eb7f12`, CPython 3.14.4,
numpy 2.5.3, scipy 1.18.1, sympy 1.14.0):

| Suite | Tests | Result | Content |
|---|---|---|---|
| `tests/unit/` | 86 | PASS | evaluators, emitters, jets, provenance guards |
| `tests/negative/` | 42 | PASS | deliberate falsification controls (must fail) |
| `tests/regression/` | 32 | PASS | Eq.85 identities, C_bg chain, KT2023 transcription |
| `tests/integration/` | 1 | PASS | end-to-end pipeline consistency |

Key symbolic regression results (exact, `== 0` via `canonical_simplify`):

| Identity | Result |
|---|---|
| `P_MH[C_bg] − C85_undivided` (Test B, symbolic A0′) | **0** — exact |
| `ΔC_εY` at `A0′ = 0` (Test A, f2Y symbolic) | **0** — exact |
| KT2023 Eq. 85 on-shell certificate (luminal branch) | **0** — exact |
| S0-slice diagonal map `KT_E00 = −HT_E00/(r²f)` | **0** — exact |
| Wrong-slot negative control (v8 → V8 column) | identity **breaks** (as required) |

**Scientific gate status — 10/28 PASS, 1 MARGINAL, 17 NOT_RUN:**

| Gate | Status | Claim |
|---|---|---|
| G00 | **PASS** | repository/environment provenance |
| G01 | **PASS** | model lock (`MODEL_LOCK.json`) |
| G02 | **PASS** | geometry regression (light-ring structure) |
| G03 | **PASS** | background algebraic rank |
| G04 | **PASS** | scalar ODE regularity |
| G05 | **PASS** | scalar integration |
| G06 | **PASS** | independent background residual validation |
| G07 | MARGINAL | holonomic action-jet validation (independent recheck open) |
| G10 | **PASS** | full symbolic `C_bg` |
| G11 | **PASS** | ε_Y null test (Test A) |
| G12 | **PASS** | exact MH Eq.85 projection (Test B) |
| G13 | NOT_RUN | production electric specialization (f2Y = 0) |
| G14 | NOT_RUN | full light-ring limit |
| G15 | NOT_RUN | Σ_SVT operator decomposition |
| G16 | NOT_RUN | full on-shell LR balance |
| G20 | NOT_RUN | same-member provenance gate |
| G30 | NOT_RUN | unreduced common H+SVT quadratic action |
| G31 | NOT_RUN | common constraint rank |
| G32 | NOT_RUN | common constraint elimination |
| G40 | NOT_RUN | kinetic K positivity |
| G50 | NOT_RUN | radial characteristics `c_r²` |
| G60 | NOT_RUN | angular characteristics `c_Ω²` |
| G70 | NOT_RUN | finite-l even sector |
| G71 | NOT_RUN | odd sector |
| G72 | NOT_RUN | vector sector |
| G80 | NOT_RUN | interfaces / patch continuation |
| G90 | NOT_RUN | global regularity |
| G100 | NOT_RUN | QNM / trapping (evaluated LAST) |

`NOT_RUN` is an honest status, not a gap in the bookkeeping: each of these
gates becomes PASS only through its own executable implementation,
independent test and machine-readable evidence.

## Historical test results

The repository keeps the **complete dated trail** of every test campaign —
they are the falsification memory of the project and remain queryable:

| Where | What | Example |
|---|---|---|
| `TEST_STATUS_*.json` | grouped test campaigns per checkpoint era | 2026-09-21 Angular-Laurent: **132/132** (unit 80, negative 42, regression 8, integration 1, release 1) |
| `FULL_CLOSURE_WORKING_STATUS_*.json` | per-era closure state (2026-09-20 → 2026-09-21, six eras) | strong-field, angular oracle, dense branch, EQ47, projected controllability |
| `CURRENT_WORKING_STATE_*.md` | human-readable per-era narratives | same eras |
| `RELEASE_STATUS_*` / `PACKAGING_STATUS_*` | release audits per era | 2026-09-20 strong-field era |
| `docs/checkpoints/` | 8 checkpoint monographs (2026-09-14 → 2026-09-16) | Horndeski/SVT closure-complete, finite-l technical memory, shared-baseline assembly validation, holonomic A5 closure |
| `data/generated/` | 32 evidence directories | Eq.85 light-ring replays, angular EQ47 projection, STEP 0–2 chain, absolute-attempt audits (09-18/19/20) |
| `ABSOLUTE_CLOSURE_LEDGER.json`, `CLOSURE_AUDIT_REPORT.json`, `FULL_PIPELINE_REPORT.json` | cross-era ledgers | every closure attempt with verdict |

Two entries deserve emphasis because they **rejected** concrete members:

- the 2026-09-17 regional electric member (central finite-L kinetic
  positivity failure), and
- the 2026-09-18 global `A0′ = 0` ε_Y member (excluded at the inner light
  ring by the denominator-free Eq.-85 audit).

**Interpretation rule:** historical PASS/FAIL records are regression
references, counterexamples and falsification witnesses **for their own
dated member** — under the same-action rule they are never production
evidence for the current member. The current verdict comes exclusively
from `GATE_STATUS.json` + the gate graph.
```

## What has been proven so far

The heart of the light-ring problem is the undivided Maxwell–Horndeski
relation **Eq. (85)** of Kase & Tsujikawa (PRD 107, 104045 (2023)): its
denominator `2f − r f′` vanishes *exactly at a light ring*, so only its
undivided form can decide whether the P5 ring is compatible with the theory.
This chain is now derived and machine-verified end to end:

1. **Primary sources transcribed verbatim** — Kase–Tsujikawa Eqs. (7)–(10),
   (14), (24), (25), Appendix A; Heisenberg–Tsujikawa Eqs. (14)–(20) — with
   regressions that pin every equation to its source.

2. **The Eq.-85-generating combination was found and verified exactly.**
   In the theory's Maxwell–Horndeski slice,

   ```
   C85_undivided  =  α·E00 + β·E22 + β·Δ22_SVT
   α = +√f/(r√h)          β = +r·f^{3/2}/√h
   ```

   holds **identically** (symbolic `== 0`, not numerically) — derived from
   the Ward/Bianchi structure, not fitted.

3. **Test A** — the historical ε_Y deformation is *background-null* at
   `A0′ = 0`: it vanishes on the zero-vector branch as a structural
   identity, while `f2Y` stays symbolic.

4. **Test B** — the full Maxwell–Horndeski projection of the same `C_bg`
   reproduces the undivided Eq. (85) **exactly**, with the electric support
   term `−2 r f h A0′² v8` intact and `A0′` symbolic. A negative control
   (wrong coefficient slot) breaks the identity — the notorious
   **v8↔V9 slot mapping** can never silently regress again.

5. **Diagnostic values** established under Eq.-85-only assumptions:
   the correctly mapped electric coefficient at the inner ring is
   **v8 = +17.52 > 0** (an earlier −16.82 came from the wrong table
   column), giving the electric channel the right sign to counter the
   positive curvature `W″ ≈ +2.59` of the stable ring.

What remains is concentrated in **one symbolic unknown** — the genuine-SVT
correction `Δ22_SVT` of the angular field equation — followed by the
perturbation-theory gates (K, `c_r²`, `c_Ω²`, finite-l, odd/vector,
interfaces, QNM). The exact next command always lives in
[`CONTINUATION_STATE.json`](data/generated/phase2_q2/CONTINUATION_STATE.json).

## Why so much machinery?

Because this field produces beautiful wrong answers. A wrong coefficient
slot once flipped the sign of the electric support and nearly killed a
viable branch. A stale derivative convention can fake a boundary artifact
into a physics failure. So the repository enforces what the physics
demands:

- **Fail-closed gate graph** — a gate is PASS only with evidence; missing
  evidence means `NOT_RUN`, never PASS. Downstream gates cannot certify on
  broken upstream state.
- **Same-action provenance** — every artifact carries the member hash;
  historical tables can inform but never substitute.
- **One canonical symbolic comparer** — exact algebraic equality, never
  fragile expression-shape matching.
- **No fitting, encoded** — boundary data come from certified
  constructions; the closure pipeline contains no optimizer.

## Getting started

CPython **3.14**, Linux:

```bash
git clone https://github.com/error-wtf/SSZ_FULL_CLOSURE.git
cd SSZ_FULL_CLOSURE
python3.14 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .

pytest -q                                                    # 162 tests
```

The symbolic closure layers run in seconds:

```bash
PYTHONPATH=src pytest tests/regression/test_light_ring_identity.py -v
PYTHONPATH=src pytest tests/regression/test_kt_mh_background.py -v
```

## Repository map

| Path | Content |
|---|---|
| `src/ssz_p5/action/` | Light-ring identity, the `C_bg` construction, KT2023 primary-source core, symbolic canonicalizer |
| `src/ssz_p5/production/` | Unsplit SVT background evaluator, Appendix-A coefficient emitter |
| `src/ssz_p5/closure/` | The fail-closed stage-gate graph and programmatic verdict |
| `tests/` | Unit, regression and scientific suites (independent of the code they verify where possible) |
| `data/generated/` | Machine-readable evidence for every claim |
| `docs/canonical/` | Frozen derivation protocol, canonical facts, supersession ledger |
| `MODEL_LOCK.json` | The one production model / member / convention lock |
| `EVIDENCE_INDEX.json` | Claim → gate → code → test → evidence → commit |
| `GATE_STATUS.json` | Machine-written gate statuses |

Superseded members and forensic audits are kept under dated status files and
`archive/` — never as production dependencies.

## Primary literature

- R. Kase, S. Tsujikawa, *Black hole perturbations in Maxwell–Horndeski
  theories*, PRD **107**, 104045 (2023), [arXiv:2301.10362](https://arxiv.org/abs/2301.10362)
- L. Heisenberg, S. Tsujikawa, *Hairy black hole solutions in U(1)
  gauge-invariant scalar-vector-tensor theories*, PLB **780**, 638 (2018),
  [arXiv:1802.07035](https://arxiv.org/abs/1802.07035)

## License

Anti-Capitalist Software License v1.4 — see [LICENSE](LICENSE).
