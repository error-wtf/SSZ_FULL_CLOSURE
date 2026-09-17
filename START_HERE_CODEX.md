# START HERE, CODEX

This is the self-contained **near-zero-work regional Full-Closure handoff**.

Do not begin by searching the filesystem or re-reading historical theory. The files that previous runs repeatedly treated as missing have already been copied into `data/prestaged/`, and canonical paths are registered in `src/ssz_p5/production/sources.py`.

## First: run the baseline once

```bash
python tools/generate_global_action_cover.py
python tools/verify_codex_handoff.py
pytest -q
python ssz_p5_full_pipeline.py --strict
```

At package creation the expected baseline was:

```text
handoff verification 30/30 PASS
pytest 11/11 PASS
strict archive/audit pipeline 44/44 PASS
```

Do not repeatedly rerun the full baseline while implementing small changes. Use focused tests, then full tests only at the milestones defined in the execution contract.

## Then read exactly these two files

```text
CODEX_FINAL_EXECUTION_CONTRACT.md
CODEX_TASK_DAG.json
```

The contract is authoritative.

## Production architecture — do not reinterpret it

```text
u < 0.5515230871346237          weak_exterior_H
0.5515230871346237 <= u < 0.61 outer_same_action_H_SVT
0.61 <= u < 0.71               central_exact_SVT
0.71 <= u < 0.715              inner_same_action_SVT_H
u >= 0.715                     punctured_H_core
r = 0                          analytic_center
```

The two previously misclassified strong-field points are central genuine-SVT points. Pure-H `A0prime=0` identities are not global gates there.

Do not promote the Strong-H witness to production. Do not globally impose `A0prime=0`. Do not use the historical global epsilon-Y action JSON as the production member for this run. The current member description is:

```text
SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json
```

## Already supplied

You already have the old action cover, outer raw 41, central Odd profile, regional coefficient references, resolved outer/inner backgrounds, central exact-SVT data, core action data, analytic center data and QNM/Jost regression artifacts.

Do not ask the authors for them.

## Forward-only goal

```text
finish/verify missing regional direct generation
→ global Direct 41/41
→ existing b2-corrected constraints
→ existing profile-aware reducer
→ Direct K,R,G,S,M
→ finite-l Odd/Even gates
→ coupled QNM from the SAME operator
→ convergence / unstable-mode search
→ Absolute Full Closure certificate
→ clean verified release ZIP
```

No research loop.
