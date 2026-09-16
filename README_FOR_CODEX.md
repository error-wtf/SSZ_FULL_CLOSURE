# SSZ P5 Full Closure — Codex Handoff

This directory is intended to be sufficient for building a clean, reproducible research repository.

## First command

```bash
python ssz_p5_full_closure_auditor.py --data-dir . --full --json build/audit.json --csv build/gates.csv
```

Expected headline:

```text
FULL_CONSTRUCTIVE_CLOSURE = PASS
DIRECT_GLOBAL_KRGM_EXPORT = OPEN
```

The second line is intentional. Codex should close it by regenerating one center-to-infinity 41-slot stream directly from the action patches and feeding that stream through the accepted profile-aware reducer. Never turn it green by renaming or concatenating historical files.

## Production action

Use `SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json` as the immutable model definition. The genuine-SVT deformation is `Delta f2 = 0.01 Y` on `A0prime=0`. It is exactly background-null and changes the quadratic vector block while leaving the scalar/metric background and their Horndeski quadratic block unchanged.

## Non-negotiable conventions

1. Assemble sectors before eliminating constraints.
2. Accepted branch: `v12 = -v6/(2h)`.
3. Accepted holonomic identity:
   `a5 = a2' - a1'' - (A0' v4 / 2)' + A0' v5 / 2`.
4. Recanonicalize the common auxiliary vector after action-level addition: `v7=v2^2/(4v1)`.
5. Use JET9D8 (`window=9`, polynomial degree `8`) for production radial derivatives unless an analytic derivative is available.
6. Preserve the antisymmetric finite-l `S` block. `R=0` in the accepted canonical convention.
7. Do not treat a numerical deep-core cancellation failure from a rounded historical 41-slot CSV as a physical ghost. Regress the regenerated operator against exact Maxwell-Horndeski invariants.
8. A QNM spectrum is publishable only after the direct-global-KRGM certificate is generated and the auditor passes with `--require-direct-krgm`.

## Directory roles

- `src/`: accepted emitters, constraints, derivative service and reducer.
- `data/authoritative/`: source witnesses and on-shell profiles.
- `data/regression/`: numerical regression targets, not primary model definitions.
- `data/certificates/`: frozen release audit.
- `paper/`: full monograph and source.
- `legacy_superseded/`: known historical traps and diagnostics.
