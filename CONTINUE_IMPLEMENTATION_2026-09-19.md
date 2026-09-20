
## Always inspect the executable truth view first

Before changing the action search, run:

```bash
python tools/show_repo.py
python tools/show_gate_matrix.py
python tools/show_member_matrix.py
```

Use `REPO_EVIDENCE_REGISTRY.json` rather than interpreting historical report
filenames.  Rejected-member audits may intentionally exit 2.

# Continue implementation — electric hybrid checkpoint

This file is the current continuation contract.  Do not restart the project from
historical coefficient dumps.

1. Preserve the verified P5 geometry, 41-slot framework, JET9D8,
   generalized-psi descriptor, descriptor pullback audit and singular DAE pencil.
2. Treat `SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json` as a rejected
   Absolute-Closure member, not as the current target.
3. Treat `SSZ_P5_HSVT_PRODUCTION_CANDIDATE_2026-09-18.json` as a rejected global
   zero-vector member.  Its component-level PASS evidence remains useful.
4. Start the active search from
   `SSZ_P5_ELECTRIC_HYBRID_SEARCH_CANDIDATE_2026-09-19.json` and the machine
   light-ring report/seed in `data/generated/absolute_attempt_2026-09-19/`.
5. The inner light-ring Eq.85 requirement is a hard background gate.  Do not set
   `A0prime=0` globally or force tensor_F=tensor_H by assumption there.
6. Map the Eq.85 electric combination to the chosen SVT action convention before
   interpreting the seed as an action control.  Do not identify similarly named
   coefficient slots without an explicit derivation.
7. Solve the electric background and Horndeski scalar/tensor controls together.
   The old f2-Hessian-only scan and simple Strong-H coefficient transplant are
   known insufficient.
8. Promote a member only after one holonomic action reproduces its background and
   Direct-41 stream and passes finite-L K/radial/angular gates for
   L=6,12,20,42,110,420,1000.
9. Reuse the validated descriptor/DAE machinery on that same member; do not return
   to the ill-conditioned deep-Core explicit `D_h1^{-1}` Schur representation.
10. Coupled QNM remains blocked until a global same-operator certificate exists.

## Added principal-search checkpoint

11. Reproduce the current best conditioned principal search with
    `python tools/audit_electric_hybrid_principal_feasibility.py`.
12. Preserve `ELECTRIC_HYBRID_HORNDESKI_PRIMITIVE_BASIS.csv`,
    `ELECTRIC_HYBRID_CENTRAL_SEARCH_BASELINE.csv`, and
    `ELECTRIC_HYBRID_PRINCIPAL_RECIPE.json`; they replace the external search
    pickle as the continuation state.
13. Do **not** continue unconstrained principal-only Newton tuning.  The stored
    negative L=6/12/20 remnants coexist with wide pivot margins, while a forced
    zero crossing was shown to create a near-singular H0/h1 pivot and huge false
    negative kinetic modes.
14. The next solver must vary the electric/background action controls jointly
    with the Horndeski scalar/tensor controls and re-evaluate the background
    identities before any Direct-41 or finite-L PASS is promoted.
