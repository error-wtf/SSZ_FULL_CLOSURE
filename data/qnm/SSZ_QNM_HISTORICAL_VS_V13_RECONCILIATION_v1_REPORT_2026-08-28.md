# Historical QNM branch vs canonical v13 — reconciliation audit v1

## Result

The old `qnm_eikonal.py` output and the current v13 strong-field fingerprint are **not the same metric branch**.

The RAG catalog identifies the historical source as

`physics/Segmented-Spacetime-Mass-Projection-Unified-Results/qnm_eikonal.py`

with SHA-256

`da76e035eb6425220bc3b2849e385328fd0049d27004fca06ddb6d1722250de6`

and marks it `canonical: false`.

The archived run reports for 30 M_sun:

- Omega_c = 1386.478000 1/s
- lambda = 1323.476000 1/s

The same archived validation block explicitly declares a different effective metric:

A(U) = 1 - 2U + 2U^2 - (24/5) U^3,
B = 1/A,
U = 1/(2x).

## 1. The historical Omega value is identified exactly

Written in x=r/r_s,

A(x)=1-1/x+1/(2x^2)-3/(5x^3).

Its circular-null condition factorizes:

2A-xA' = (2x-3)(x^2+1)/x^3.

Therefore its unique positive real circular-null radius is exactly

x_ph = 3/2.

At that point

A(3/2) = 0.377777777777778

and

Omega_c r_s/c = sqrt(A)/x = 0.409757531435240.

For 30 M_sun this gives

Omega_c = 1386.478171 1/s,

which reproduces the recorded 1386.478000 1/s to relative error
1.232e-07.

So the historical Omega value is positively identified as a prediction of the old polynomial metric, not the current Xi/Hermite v13 metric.

## 2. The historical lambda value is not metric-consistent

For a static spherical metric with B=1/A, the coordinate-time instability exponent of the unstable null orbit is

lambda^2 = A_m [2A_m-r_m^2 A_m'']/(2 r_m^2).

Applied to the same historical polynomial metric:

lambda r_s/c = 0.492467263365768,

or for 30 M_sun

lambda = 1666.339379 1/s.

But the archived script recorded

lambda = 1323.476000 1/s,

equivalent to

lambda r_s/c = 0.391137971057999.

That is -20.58% below the metric-consistent value.

Therefore the archived lambda cannot be used as the Lyapunov exponent of that polynomial metric without recovering and justifying the exact historical formula.

## 3. It did not numerically 'match GR'

For 30 M_sun the Schwarzschild eikonal rates are

Omega_GR = lambda_GR = 1302.369464 1/s.

The historical recorded Omega is therefore +6.46% from GR,
and the recorded lambda is +1.62% from GR.

Thus the old documentation phrase `match GR's eikonal relations` cannot mean equality of the numerical GR values. At most it describes a GR-style eikonal construction.

## 4. Canonical v13 is a different branch

Current canonical v13:

- x_u = 2.130539516826896
- b_c/r_s = 2.656525211989155
- Omega_c r_s/c = 0.376431586452447
- lambda r_s/c = 1.855517555304344
- lambda/Omega = 4.929229166954
- a_bar = 0.202871476681

Historical polynomial metric, treated consistently:

- x_u = 1.5
- b_c/r_s = 2.440467650459881
- Omega_c r_s/c = 0.409757531435240
- lambda r_s/c = 0.492467263365768
- lambda/Omega = 1.201850425155
- a_bar = 0.832050294338

These are visibly different geometries.

## Scientific disposition

1. The old `qnm_eikonal.py` output is **historical/supporting**, not a canonical-v13 strong-field result.
2. Its Omega value is understood and reproducible from the old polynomial metric.
3. Its lambda value is **superseded/unresolved** as a Lyapunov prediction because it fails the metric-consistency check.
4. The current v13 `Omega_c`, `lambda`, `a_bar` fingerprint remains a metric-level null-geodesic result.
5. Calling it a gravitational-wave QNM spectrum remains conditional until a perturbation sector is derived from the same SSZ action/completion.
