# Regional Full-SVT execution result — 2026-09-17

The current selected production member is the six-region cover in
`SSZ_P5_REGIONAL_PRODUCTION_MEMBER_2026-09-17.json`. The historical globally
zero-vector epsilon-Y member and the strong-H control witness are not used as
its central background.

## Verified input and implemented steps

The supplied Near-Zero-Work ZIP has SHA-256
`b2c5307256d2e9c2c78cebcd83ce7035c6121ad2c3dd3f679409709754db0021`.
A fresh extraction and editable installation passed. Its original handoff
verification passed 30/30 and pytest passed 11/11. The original archive strict
pipeline returned 43/44, not 44/44, because a module attempted to import the old
absolute `/mnt/data/ssz_p5_higher_jet_closure_2026-09-16.py` path. Integration uses
the existing portable repository loaders and preserves the tested b2 constraint
fix absent from the supplied constraint-map copy.

Implemented in the active repository:

- Regional member and source registry wired into pipeline/certificate validation.
  Electric central and handover backgrounds are accepted in their proper regions;
  the obsolete global A0prime=0 requirement is removed from direct verification.
- Outer pure-SVT 41-slot sector independently regenerated from resolved background
  jets and the existing S-weighted transverse Hessian prescription. Maximum scaled
  difference against the staged raw reference: about 9.03e-9. This is a sector
  regression, not the complete Horndeski+SVT assembly or a global certificate.
  The historical controls' constant continuation below u=0.57 is recorded explicitly.
- Central selected coefficients regenerated from the original unreduced table and
  its existing selected lower-order jets. Encoded nonzero v5, c3 and e3 are retained.
  Maximum scaled difference against the staged corrected reference: about 6.89e-12.
- Common unreduced assembly helper subtracts the shared baseline exactly once,
  requires identical grids/backgrounds, then applies the fixed identities.
- Independent kinetic Schur calculation, finite-L scan and convergence audit.
- The existing JET service now supports its declared 11/8 and 13/8 convergence
  stencils; the production 9/8 path is unchanged. Higher differential orders are
  classified by total derivative order, so mixed orders (2,1) and (2,2) are visible.

## Reproducible necessary gate failure

Region: `central_exact_SVT`.
Coordinates: `u=0.6900037509377344`, `x=r/r_s=1.4492674838955182`.
Multipole parameter: `L=6` (ell=2).
Basis: `(psi, dphi, V)`.
The vector background is nonzero. No Horndeski-only identity is used here.

The kinetic matrix at this interior point has eigenvalues approximately

```text
-43.86167177, 0.09070493, 115.42430929
```

A normalized negative direction is approximately
`(0.04532109, 0.96998501, -0.23890393)`.
It violates the required `K > 0` condition in the actual selected central
coefficient member. The concrete matrix and full-precision values are in
`data/diagnostic/REGIONAL_CENTRAL_KINETIC_GATE.json`.

The independent implementation forms the kinetic quadratic form of the original
41-slot action in `(dot Y, dot Y')`, solves its H1/dA1 two-by-two auxiliary block,
and performs the radial integration by parts. It does not call the constraint-map
or Euler-operator implementation. Both routes agree to better than 5e-8 maximum
scaled error on the audited interior `0.62 < u < 0.70`. The residual dot-Y-prime
quadratic block and antisymmetric mixed kinetic block vanish to floating-point
precision. Constraint pivots are nonzero. Nonsingular constant field rescalings
preserve the negative inertia.

At the same interior point the negative eigenvalue remains between approximately
-43.869 and -43.829 for stride 1/2/4 and stencils 7/6, 9/8, 11/8 and 13/8.
The check also integrates a localized smooth velocity packet directly in the
Schur-reduced kinetic action and compares it with the canonical kinetic integral.
The machine report contains both results, all source hashes and implementation
hashes. This is a robust sign failure, not a claim of arbitrarily precise
continuum eigenvalues from rounded CSV data.

The full required L set is scanned over the selected central interval
`0.61 <= u < 0.71`. The kinetic condition fails for L=6,12,20,42. It passes this
local necessary check for L=110,420,1000; those local passes do not certify radial,
angular, global or spectral stability.

## Consequence and scope

The requested all-PASS production chain cannot be certified while this central
coefficient member is locked. Completing other regions cannot remove a negative
kinetic direction supported wholly inside the central region. No new core search,
background inverse problem or historical QNM calculation is performed in response.

This is not a no-go theorem for P5 geometry, all Full-SVT actions, or all possible
regional members. A changed central action/member would require a separately
specified and validated new selection. This implementation does not silently
replace the frozen central action to obtain PASS.

No global direct-41 certificate, global KRGSM certificate, coupled QNM convergence
certificate, Absolute Full Closure certificate or v1.0.0 closure ZIP is issued.
The remaining global generators are not claimed complete.

## Reproduction

```bash
python tools/regenerate_regional_coefficients.py
python tools/audit_regional_kinetic.py
python ssz_p5_full_pipeline.py --strict
```

The regeneration command passes its stated sector/normalization regression.
The kinetic command returns exit code 2 for the recorded failed necessary gate.
Strict mode also returns 2 and preserves rejection of absent global/spectral
certificates. The default pipeline remains a separately labelled archive audit.
