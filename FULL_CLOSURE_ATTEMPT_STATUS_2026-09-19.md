# SSZ P5 Full-Closure Attempt — 2026-09-19

This checkpoint deliberately distinguishes **reproducible progress** from a final closure certificate.

## Verified in this archive

- Existing + new regression suite: **118/118 PASS**.
- Reconstructed electric-hybrid Central single-action stream: **PASS**.
- Central exported action jets: **3,809 rows**.
- Central exported Direct-41 stream: **3,809 rows**.
- Finite-L Central kinetic/radial-principal gate: **7/7 PASS**.
- C-infinity interface coordinate system implemented.
- One-control `f4` handover: **rejected by reproducible falsification gate**.
- Compact joint `(a1,c2,c4,F,G,H)` response space: **rank 6/6 on both buffers**.

## Not yet certified

`ABSOLUTE_FULL_CLOSURE = false`.

The direct verification call currently stops at:

    RuntimeError: QNM disabled: invalid or absent direct-global-KRGM certificate: [Errno 2] No such file or directory: '/mnt/data/ssz_fullclosure_attempt_build/SSZ_FULL_CLOSURE/data/certificates/DIRECT_GLOBAL_KRGM_CERTIFICATE.json'

The scientific interface audit independently records:

    status = PENDING_JOINT_BACKGROUND_CONSTRAINED_SOLVE

The next required constructive step is therefore the joint background-constrained SVT + Horndeski solve in `0.61<u<0.62` and `0.70<u<0.71`, followed by Direct-41 re-emission, angular/global KRGSM and same-operator QNM convergence.

## Important files

- `src/ssz_p5/production/electric_hybrid_onshell_central.py`
- `src/ssz_p5/production/action_handover.py`
- `tools/audit_action_handover_control_space.py`
- `data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_ACTION_JETS.csv`
- `data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_41.csv`
- `data/generated/absolute_attempt_2026-09-19/ELECTRIC_HYBRID_ONSHELL_CENTRAL_AUDIT.json`
- `data/generated/absolute_attempt_2026-09-19/ACTION_HANDOVER_CONTROL_SPACE_AUDIT.json`
- `tests/regression/test_full_closure_attempt_20260919.py`
