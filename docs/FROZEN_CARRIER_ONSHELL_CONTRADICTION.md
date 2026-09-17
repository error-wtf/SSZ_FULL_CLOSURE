> **Scope correction (2026-09-17):** This report evaluates a pure-H witness file. It is not a global production-member gate. The sector-aware production assignment places both disputed radii in `central_exact_SVT` with nonzero `A0prime`; see [PRODUCTION_REGION_ASSIGNMENT.json](../data/diagnostic/PRODUCTION_REGION_ASSIGNMENT.json).

# Reproducible inconsistency in the prescribed zero-vector carrier

**Status: FAIL_ON_SHELL_IDENTITY.** This is a necessary background-equation
check on the supplied production input. It is not a coupled-spectrum result,
a perturbative instability, or a no-go theorem for the P5 geometry.

The existing 51-test implementation at `ad02fb0` was the starting point. The
constraint `b2` repair and all frozen source files remain unchanged. The general
core emitter has not been completed: the premise that it is the only remaining
obstacle fails an independent check in the already supplied strong carrier.

## Precisely which assumptions conflict

These simultaneously prescribed inputs are tested:

1. The `f`, `h`, and `H_equals_F_equals_G` functions in
   `data/production/ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv`.
2. `a4 = sqrt(f*h)*H/2`, the common Horndeski coefficient convention.
3. The final action JSON's `A0prime = 0`, with no additional matter sector.
4. The claim that this same carrier is on shell.

The `epsilon_Y * Y` correction has zero background variation at vanishing
field strength: its first variation still contains at least one background
field-strength factor. It cannot supply the missing background equation.

## Necessary equation, before any constraint reduction

[Kase and Tsujikawa, *Black hole perturbations in Maxwell-Horndeski theories*,
arXiv:2301.10362](https://arxiv.org/html/2301.10362) use
`ds² = -f dt² + dr²/h + r² dΩ²`. Their Eq. (85), following from background
Eqs. (7) and (9), requires

```text
(2f - r f') a4'
  = (r f'' - r f'^2/f + 2f' - 2f/r) a4
    + f^(3/2) F/(r sqrt(h)) - 2 r f h A0'^2 v8.
```

Define `R85` as the left side minus the right side. The zero-vector carrier
requires `R85 = 0`. The implementation multiplies through by `2f-rf'`, so no
light-ring division is involved. It also never divides by the scalar derivative.
The chosen witness is an interior point, not a center, endpoint, or transition.
The source uses arbitrary Horndeski functions; this identity does not depend on
the unresolved higher scalar jets controlling `c2`.

## Direct witness

At `r/rs = 1.5999977950251312`, the stored row gives

| Input | Value |
| --- | ---: |
| f | 0.3739125409549743 |
| h | 0.4015751947977254 |
| F = H | 0.8777994828489721 |
| G4 | 0.438899741424486 |
| G4_phi | 0.61633042306571 |

The production JET9D8 service gives

```text
left side  =  0.006453887622634062
right side = -0.06452053578102868
R85        =  0.07097442340366274
```

Every comparison uses the same original row; stride decimation changes the
surrounding stencil, not the physical point:

| Evaluation | R85 |
| --- | ---: |
| JET9D8, original sampling | 0.07097442340366274 |
| JET9D8, stride 2 | 0.07097442586934617 |
| JET9D8, stride 3 | 0.07097442578583507 |
| JET9D8, stride 4 | 0.07097442550915813 |
| Action-jet/scalar chain rule | 0.07097442116687630 |

The last route obtains `phi' = -sqrt(-2X/h)`, `H' = 2 G4_phi phi'`, and
`f'`, `f''` by differentiating the frozen `f=(1+phi)^(-2)` closure. Only `h'`
and `phi''` use the existing JET9D8 service. It does not numerically differentiate
`H` or twice differentiate `f`. The supplied decreasing scalar fixes the sign.
The stored `H=2G4` and scalar closure agree to machine precision.

The route spread is approximately `4.70e-9`, versus a residual of `7.10e-2`.
The report keeps the unchanged `background_residual_abs=1e-10` and separately
reports a conservative diagnostic sensitivity allowance of 100 times this
spread (approximately `4.70e-7`). This is a numerical sensitivity estimate,
**not a rigorous interval-arithmetic error enclosure**. No policy tolerance is
relaxed. The discrepancy is resolved by more than five orders of magnitude
above that allowance.

## Independent background-equation subtraction

For the supplied luminal strong chart `G4=G4(phi), G5=0`, let `H=2G4`.
Substituting in the source's Eq. (10) gives

```text
C2=C3=C5=0, C4=-H, C7=-2hH', C8=-H(h-1), C10=-Hh/2.
```

Subtracting Eq. (9) from Eq. (7) cancels `C1` and `C6`, including all `G2`
and `G3` freedom. On `A0'=0`, the resulting necessary equation is

```text
E00-E22 = -Hh'/(2r) - hH'/r - H(h-1)/r²
          + hHf''/(2f) - hHf'^2/(4f²) + Hf'h'/(4f)
          + hH'f'/(2f) + Hhf'/(2rf) = 0.
```

At the same row it evaluates to `-0.12294499883905527`. Direct algebra gives

```text
R85 = -r f sqrt(f/h) * (E00-E22).
```

The independently coded expressions agree to `2.78e-17` at the witness. This
also explains why choosing missing `G2XX`, `G3XX`, or `c2` values cannot repair
this particular contradiction. A different core outside this interior chart
cannot change the local equation while these prescribed functions stay fixed.

## Controls and reproduction

```bash
python tools/audit_frozen_onshell.py
pytest -q tests/unit/test_onshell_identity.py tests/regression/test_frozen_onshell_diagnostic.py
python ssz_p5_full_pipeline.py --strict
```

The diagnostic exits **2** and writes `build/FROZEN_ONSHELL_IDENTITY.json`.
The committed snapshot at `data/diagnostic/FROZEN_ONSHELL_IDENTITY.json`
records all source SHA256 hashes, the exact row and five evaluation routes.
No selected 41-slot output or reduced matrix is an input.

Analytic tests check Minkowski, Schwarzschild, de Sitter and the correctly
normalized electric Reissner–Nordström solution. For Reissner–Nordström,
omitting its electric term correctly creates a nonzero residual. A varying-G4
analytic off-shell example tests the separate background subtraction. A
light-ring test checks that no division by `2f-rf'` occurs.

Seven new tests pass because they check the implementation and correct failure
detection. They do **not** turn the failed physical identity into a PASS.
Strict mode recomputes the witness rather than trusting the committed report.
The archive auditor remains a separate historical regression.

## Consequence for the requested release

The requested stop condition “a precise reproducible contradiction of the
frozen member” is reached under the specified identification of this CSV as
the final carrier. The narrower claim that the missing general core emitter
alone prevents absolute closure is no longer supported.

No replacement action, electric branch, adjusted metric, or altered H profile
has been selected. The original JSON's historical constructive-PASS label is
preserved as source material; it does not override the new numerical check.
A general Horndeski emitter can still be a useful software component, but cannot
certify this exact prescribed background as on shell. No final coupled QNM
frequencies, absolute PASS certificate, or absolute v1.0.0 ZIP are issued.

## Verification of this implementation checkpoint

- `pytest -q`: **58 passed** (the original 51 plus seven diagnostic/control cases).
- Ruff: **PASS** for the package, tests, diagnostic command, and pipeline.
- `python ssz_p5_full_pipeline.py --strict`: **47/49**, exit **2**.
- The two failing checks are the necessary frozen-carrier on-shell identity
  and the missing complete direct-production certificate.
- Manifest and SHA256 inventory are regenerated after the final source edits.
