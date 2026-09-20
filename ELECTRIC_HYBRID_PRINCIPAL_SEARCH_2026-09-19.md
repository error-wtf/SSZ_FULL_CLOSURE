# Electric-hybrid principal search checkpoint — 2026-09-19

This checkpoint answers one narrow question: is the old Central finite-L kinetic
failure controllable by smooth Horndeski primitive responses while the electric
SVT background is held fixed?

## Result

Yes, strongly, but not completely.  The exact frozen search recipe reduces the
full-grid `0.62<u<0.70` minima to:

| L | min eigenvalue K | u at minimum |
|---:|---:|---:|
| 6 | -3.2429460101 | 0.6999287322 |
| 12 | -0.1309937268 | 0.6999287322 |
| 20 | -0.0128773635 | 0.6999812453 |
| 42 | +0.0012545875 | 0.6999812453 |
| 110 | +0.0004647683 | 0.6999812453 |
| 420 | +3.98849e-5 | 0.6999812453 |
| 1000 | +7.31613e-6 | 0.6999812453 |

The low-L constraint pivots remain wide; at L=6,
`min|Dh1|≈151.06` and `min|auxiliary determinant|≈441.88`.

## What was learned

- A single global-amplitude control already moved the robust u≈0.69 L=6 mode
  from about -43.86 to O(-1).
- Smooth polynomial and C-infinity ramp shapes move all remaining negative modes
  to the same upper-Central edge without harming high L.
- Narrow localized bumps are rejected: derivative jets create new interior fails.
- An unconstrained Newton solve can make selected eigenvalues look nearly positive
  only by driving a constraint pivot toward zero; full-grid K then develops
  O(1e8) negative modes.  That branch is explicitly rejected.
- Further principal-only shaping gives diminishing returns.  The next physically
  meaningful control is the electric/background action sector required already by
  the inner-light-ring Eq.85 audit.

## Reproduce

```bash
python tools/audit_electric_hybrid_principal_feasibility.py
pytest -q tests/unit/test_electric_hybrid_controls.py \
          tests/regression/test_electric_hybrid_principal_checkpoint.py
```

This is **not** a background/on-shell, Direct-41, KRGSM, QNM or Absolute-Full-
Closure certificate.
