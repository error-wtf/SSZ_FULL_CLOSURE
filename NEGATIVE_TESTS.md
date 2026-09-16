# Mandatory negative regression suite

These are not optional. The repository is required to prove that known historical mistakes fail loudly.

1. `v12` plus-branch injection must fail.
2. `a5` without `-a1''` must fail finite-l regression.
3. reduce-then-add must disagree with combine-before-reduce and be rejected.
4. `legacy_superseded/` input must be rejected by provenance guard.
5. NaN/Inf injection into any 41 slot must fail.
6. zero `Dh1`, `DeltaV`, or `2Lv9` must abort reduction.
7. choose `epsilon_Y` so `ZA<=0` somewhere: vector stability must fail.
8. QNM call without direct certificate must abort.
9. forged certificate with wrong artifact hash must fail.
10. historical electric-SVT `ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv` must be refused as the final production stream.
11. old variable-G4 transplant `OUTER_FINAL_SAME_ACTION*` must be refused as a production source.
12. deep-core rounded direct-K negativity must not override an exact invariant oracle unless regenerated action-level coefficients also fail.
