# Canonical v13 eikonal-QNM proxy v1

## Purpose

This is a deliberately **conditional** extension of the corrected static-v13 null-orbit fingerprint. It does not use the superseded historical `qnm_eikonal.py` output.

Assume only the standard eikonal correspondence

`omega_ln ≈ l Omega_c - i (n+1/2) lambda`

for `l >> 1`.

Then the current metric fingerprint implies a sharp prediction for the *shape* of a possible high-l QNM spectrum.

## Metric inputs

GR:
- Omega_c r_s/c = 0.384900179460
- lambda r_s/c = 0.384900179460
- a_bar = Omega/lambda = 1.000000000000

Canonical SSZ v13:
- Omega_c r_s/c = 0.376431586452
- lambda r_s/c = 1.855517555304
- a_bar = Omega/lambda = 0.202871476681

## Mass-independent same-mode discriminator

For the same `(l,n)` and black-hole mass:

- frequency ratio SSZ/GR = **0.977998**
- damping-rate ratio SSZ/GR = **4.820776**
- damping-time ratio SSZ/GR = **0.207435**
- quality-factor ratio SSZ/GR = **0.202871**

Thus the eikonal oscillation frequency is only about -2.20% shifted,
while the amplitude damping rate would be about **4.82 times larger**.

That is the same strong-field pattern already seen geometrically:
small change in orbit frequency, large change in instability.

## Example: 30 M_sun, l=10, n=0

GR:
- f = 2072.785 Hz
- amplitude e-fold time = 1.536 ms
- Q = 10.000
- cycles/e-fold = 3.183

SSZ v13:
- f = 2027.180 Hz
- amplitude e-fold time = 0.319 ms
- Q = 2.029
- cycles/e-fold = 0.646

So under the correspondence, the high-l SSZ mode would oscillate at nearly the same frequency scale but lose amplitude much faster.

## Important limitation

This is **not** an l=2 gravitational-wave prediction.

The current FINAL MASTER still lists the quadratic perturbation action, ghosts/gradients/hyperbolicity/characteristics as open. Until the perturbation equations are derived from the same covariant SSZ completion, the eikonal formula is a metric proxy only.

A future action-derived perturbation sector gives a clean falsification gate:

1. derive the master perturbation equation,
2. take its high-l limit,
3. test whether its potential peak approaches the same corrected unstable null orbit,
4. compare its asymptotic Re(omega) and Im(omega) to this frozen proxy.

If it does not, the null-geodesic/QNM correspondence fails for that SSZ completion and this proxy must not be promoted to a physical QNM prediction.
