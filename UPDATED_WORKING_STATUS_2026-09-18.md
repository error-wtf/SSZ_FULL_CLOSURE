# SSZ P5 — Updated Working Closure Checkpoint

**Date:** 2026-09-18
**Source of truth:** `SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json`

## Verified software/repository state

- `pytest`: **86/86 PASS**.
- full non-strict repository pipeline: **45/45 PASS**.
- constructive full-closure auditor: **PASS / exit 0**.
- strict direct-global KRG/M implementation auditor: **OPEN / exit 3** only because `DIRECT_GLOBAL_KRGM_EXPORT` is not yet certified.
- locked conventions preserved: JET9D8, `v12=-v6/(2h)`, holonomic `a5`.

## New closure work in this checkpoint

### Inner lower-order action control — PASS

The lower-order target profiles `(v5,c3,e3)` are no longer edited as coefficient columns.  They are realized through an explicit background-null normal-jet inverse using

- `f2phiF_control`,
- `f2Xphi_control`,
- `f2phiphi_control`.

The re-emitted lower-order map has a maximum scaled replay residual of approximately `1.20e-11`; stored inner background residuals remain at the `1e-15` level.

### Inner principal Hessian action control — PASS

The principal target profiles `(v1,v4,c2)` are now realized by the raw triangular Hessian inverse into `(f2FF,f2XF,f2XX)` and re-emitted through the Appendix-A coefficient map.  Maximum scaled principal replay residual is approximately `3.38e-16`.

At the central interface, the controlled values `v5,c3,e3,v1,v4,c2` match to roughly `1e-12` or better.  At the core interface all controlled values are at machine precision except the vector-free `v1` endpoint, whose residual is approximately `3.4e-9`.

## Remaining strict Inner Direct-41 blocker

`INNER_DIRECT_41` is still **FAIL**, but its scope is now restricted to derivative-sensitive higher-jet slots:

- `a5`,
- `d3`,
- `e4`,
- `v13`,
- small residual in `v9`.

Lower-order control, principal Hessian control, slot schema/normalization, finiteness, and stored background residuals are already passing.  The remaining task is higher mixed/holonomic action-jet replay and endpoint-jet matching; it is no longer a lower-order or principal-control problem.

## Next closure chain

1. Realize the remaining Inner higher mixed/holonomic jets at action level and re-emit `a5,d3,e4,v13,v9`.
2. Promote `INNER_DIRECT_41` to PASS.
3. Regenerate/audit Core Direct-41 including the G5/quintic subcore and analytic-center match.
4. Assemble one direct center-to-infinity 41-slot stream.
5. Reduce the same stream to global finite-l `K,R,G,S,M`.
6. Run coupled same-operator QNM/Jost + unstable-mode scan.
7. Produce strict direct-global certificate and final reproducible release ZIP.

This checkpoint is intentionally **not** labelled `PERFECT_PASS`; no tolerance was relaxed and no failed gate was overwritten.
