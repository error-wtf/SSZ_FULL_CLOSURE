# Current-member spectral preflight — 2026-09-29

Branch: `current-member-spectral-preflight-20260929`

## Purpose

Test the current frozen electric production member with the same finite-L
canonical reducer beyond the certified production window, without promoting the
result to a direct-global KRGM/QNM certificate.

## Rebuilt member

`ELECTRIC_PRODUCTION_MEMBER_CURRENT`

Pinned member hash:

`8bd460ef022a9cdbcc3644abd8aecbfbb910f8e364ac1410378d2641291559cf`

Stored domain: `0.61 <= u <= 0.7099849962490622`

Certified production window: `0.62 < u < 0.70`

Multipoles: `L = 6,12,20,42,110,420,1000`

## Results

The certified production window replays cleanly: K is positive and the minimum
generalized radial characteristic is positive for all seven required L.

The same rebuilt member/reducer does **not** remain kinetically healthy over the
full stored domain. The first inward K=0 crossing occurs at approximately

- L=6:    u = 0.701306416
- L=12:   u = 0.701306406
- L=20:   u = 0.701306382
- L=42:   u = 0.701306242
- L=110:  u = 0.701305147
- L=420:  u = 0.701294238
- L=1000: u = 0.701294091

Thus the loss of K positivity is essentially common across the multipole
ladder and occurs before the inner stable light ring at u ~= 0.7061346.

At the outer light ring u=2/3, interpolated min eig(K) remains positive for all
L. At the inner stable ring, it is negative for all L:

| L | min eig(K), outer ring | min eig(K), inner stable ring |
|---:|---:|---:|
| 6 | 7.9965461e-2 | -6.9107060e-2 |
| 12 | 3.7441255e-2 | -8.9678660e-2 |
| 20 | 1.6894726e-2 | -3.7077410e-1 |
| 42 | 4.3952858e-3 | -6.1098333 |
| 110 | 6.7581796e-4 | -9.8561198 |
| 420 | 4.7168386e-5 | -11.4477796 |
| 1000 | 8.3454838e-6 | -11.7640771 |

The full stored-domain minimum again occurs at the inner endpoint for L=1000:

`min eig(K) = -27.98282633548 at u=0.7099849962490622`.

## Interpretation

This is a diagnostic extension outside the certified production window, not a
new global certificate. It nevertheless establishes a hard local obstruction
for the proposed inner-ring spectral-selection test on the current member:
the exact radial region containing the stable trapping ring is already outside
the positive-K sector of this finite-L operator.

Therefore the current member cannot be used to define a positive kinetic
normalization for a physical Z_n(r) test at the inner stable ring.

The correct chain is:

`production-window principal gates = PASS`

but

`current-member inner stable ring finite-L K = GHOST FAIL`

and independently

`DIRECT_GLOBAL_KRGM_CERTIFICATE = absent`.

So the requested global spectral-selection observable remains unevaluable for
the current member, with an additional upstream reason: the stored same-action
member becomes kinetically indefinite before reaching the inner trapping ring.

## Scope inconsistency found in repository wording

`tests/regression/test_g50_g70_g90_member.py` evaluates G50/G70 K and c_r^2
only on the production mask `0.62 < u < 0.70`, yet its module/test wording says
there is "no finite-l instability anywhere on the member." That wording is
stronger than the executed test. G40 is correctly described elsewhere as
production-window scoped.

The ABSOLUTE/TRUE closure result should therefore be read as closure of the
registered production-window validation corpus, not as a proof of finite-L
kinetic positivity over the whole stored member domain or over a
center-to-infinity global spectral operator.
