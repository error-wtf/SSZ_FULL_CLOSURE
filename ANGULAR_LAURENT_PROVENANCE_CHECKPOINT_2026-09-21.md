# SSZ P5 — Angular Laurent provenance checkpoint — 2026-09-21

This checkpoint extends the angular method infrastructure only.  It does not issue a Dense-branch angular PASS/FAIL and does not advance the strong-field continuation.

## New full provenance pipeline

The common 41-slot action is now propagated formally in `eps=1/L` through the accepted constraint maps, radial product rules, adjoints, and Euler reduction.  The same code exports the leading cubic angular characteristic polynomial and an Eq83-style coefficient map.

The pipeline then compares, without altering the action:

1. finite-L common-action reducer vs formal Laurent reducer;
2. project mass convention vs paper mass convention;
3. action-derived `m5_minus` vs published `m5_plus` comparator;
4. componentwise subleading mass map;
5. raw characteristic polynomial vs published coupled-angular reconstruction from the same extracted coefficients.

## Fresh findings from this safety-lineage replay

- Laurent vs finite-L regression: PASS at sub-`1e-6` scaled error over the certified genuine-SVT interval.
- Exact symbolic source provenance fixes `m5_action = a4*v13 - a6*v6`.
- `M22_minus-M22_plus = a6*v6*v13/(a4*v9)` is reproduced to about `1.8e-13` worst scaled error.
- After the explicit `M_paper=-M_project` convention bridge, `M11_0`, `M12_0`, `M13_0`, `M22_0`, and `M33_0` track the action-derived minus structure tightly.
- The first growing subleading mapping discrepancy is `M23_0` (median about `5.6e-7`, p99 about `6.7e-3`, maximum about `3.7e-2` on `0.62<=u<0.70`).
- The normalized raw characteristic polynomial and the published coupled reconstruction still disagree strongly, so the remaining issue is upstream of root sorting.

At the representative point `u≈0.650015`, this safety-lineage replay gives raw common-action roots approximately

    (-208.724, 0.999637, 549.058)

while the published shortcut reconstruction from the same extracted coefficient set gives

    (1, 168.199, 7721.565).

These numbers are checkpoint-specific diagnostics, not a physical angular verdict.

## Maxwell-Horndeski note

The conversation records a later exact GM-GHS formal-Laurent PASS with characteristic `(1-z)^3`, but the exact input artifact is not present in the attached safety ZIP or discoverable in the available Library search.  This checkpoint therefore provides the exact oracle hook but does not fabricate a machine replay.

## Gate state

- Dense K/R branch: preserved.
- Dense angular gate: OPEN.
- Continuation `0.708 -> 0.715`: HELD.
- Absolute Full Closure: false.

Next: close the componentwise common-action Laurent -> published `tilde K/M` convention mapping, beginning with `M23_0`.
