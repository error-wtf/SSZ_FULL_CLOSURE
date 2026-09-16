# P5 Post-Closure Spectral Probe — 14 September 2026

## Status

This note evaluates what can be derived directly from the frozen P5 geometry after local same-action closure, without selecting an additional unique global Horndeski–SVT representative. It therefore separates geometric/eikonal and test-field spectral statements from the still representative-dependent coupled tensor–scalar–vector spectrum.

## 1. Outer-ring eikonal frequency scale

At the outer P5 light ring, the frozen closure gives

$$
x_c=\frac{3}{2},\qquad f_c=h_c=\frac13,\qquad W''(u_c)=-2.
$$

The coordinate-time orbital frequency is

$$
\Omega_c r_s=\frac{\sqrt{f_c}}{x_c}=\frac{2}{3\sqrt3}\simeq0.384900179460.
$$

Because the local optical curvature is Schwarzschild-matched, the radial Lyapunov exponent is likewise

$$
\lambda_c r_s=\frac{2}{3\sqrt3}\simeq0.384900179460.
$$

Hence the leading eikonal family is

$$
\omega r_s\simeq \left(\ell+\frac12\right)\Omega_c r_s
-i\left(n+\frac12\right)\lambda_c r_s.
$$

This is a leading eikonal statement, not the final finite-\(\ell\) coupled HSVT spectrum.

## 2. Finite-\(\ell\) massless test-scalar potential

For

$$
ds^2=-f(r)dt^2+\frac{dr^2}{h(r)}+r^2d\Omega^2,
$$

with tortoise coordinate

$$
\frac{dr_*}{dr}=\frac{1}{\sqrt{f h}},
$$

a minimally coupled massless scalar separated as

$$
\Phi=e^{-i\omega t}Y_{\ell m}(\theta,\varphi)\frac{\psi(r)}{r}
$$

obeys

$$
\frac{d^2\psi}{dr_*^2}+\left[\omega^2-V_\ell(r)\right]\psi=0,
$$

with

$$
V_\ell(r)=f(r)\left[
\frac{\ell(\ell+1)}{r^2}
+\frac{1}{r\sqrt{f/h}}\frac{d}{dr}\sqrt{fh}
\right].
$$

The P5 profile produces a local well inside the outer barrier. To test whether this well can support a semiclassical trapped level, define

$$
I(E)=\int_{r_1(E)}^{r_2(E)}\sqrt{E-V_\ell(r)}\,dr_*.
$$

Bohr–Sommerfeld/WKB quantization requires

$$
I(E)=\pi\left(n+\frac12\right).
$$

The maximum well action before the outer barrier opens is

$$
I_{\max}=I(V_{\max}).
$$

A ground state therefore requires

$$
I_{\max}\ge\frac{\pi}{2}.
$$

The scan of the certified P5 profile gives the first crossing at

$$
\boxed{\ell_{\rm threshold}=256},
$$

with

$$
I_{\max}(256)=1.575048979904>\frac{\pi}{2}.
$$

Thus, in this minimally coupled test-scalar proxy,

$$
\boxed{\text{no WKB-trapped ground state for }2\le\ell\le255.}
$$

Representative values are:

| \(\ell\) | \(I_{\max}\) | Ground-state WKB trapping? |
|---:|---:|:---|
| 2 | 0.382223 | No |
| 3 | 0.332322 | No |
| 4 | 0.298897 | No |
| 5 | 0.274376 | No |
| 10 | 0.214727 | No |
| 20 | 0.207652 | No |
| 50 | 0.342698 | No |
| 100 | 0.631571 | No |
| 150 | 0.931625 | No |
| 200 | 1.234506 | No |
| 250 | 1.538522 | No |
| 256 | 1.575049 | Yes |
| 300 | 1.843103 | Yes |
| 500 | 3.063715 | Yes |
| 770 | 4.713417 | Yes |
| 1000 | 6.119231 | Yes |

The important interpretation is that the shallow geometric stable-light-ring well does not automatically imply an extremely long-lived low-multipole test-scalar resonance. The first semiclassical trapped test-scalar level appears only deep in the eikonal regime.

## 3. Why this is encouraging but not yet the physical HSVT QNM proof

The physical P5 theory contains coupled tensor, scalar and vector perturbations. Their finite-\(\ell\) potential matrix contains lower-order mass/mixing terms not fixed by the null potential alone. Therefore the threshold \(\ell\simeq256\) is a geometric/test-field diagnostic, not a proof that the coupled HSVT modes have exactly the same threshold.

The result nevertheless weakens the naive inference

$$
\text{stable light ring}\Longrightarrow\text{dangerous low-}\ell\text{ resonance}.
$$

For the actual P5 well this implication fails for the minimally coupled scalar proxy.

## 4. What remains for the exact coupled spectrum

The closure theorem proves existence of a healthy local same-action family. A unique global finite-\(\ell\) spectrum requires choosing one explicit member of that family and exporting the fully reduced matrices

$$
K(r),\qquad R(r),\qquad G(r),\qquad M(r)
$$

through the complete regular center-to-infinity domain.

The exact next calculation is therefore:

1. select one explicit same-action representative;
2. export the odd and even reduced matrices globally;
3. canonically normalize the kinetic matrix;
4. construct the coupled finite-\(\ell\) radial ODE system;
5. impose regular-center boundary conditions and outgoing conditions at infinity;
6. compute the complex QNM spectrum and search explicitly for \(\operatorname{Im}\omega>0\);
7. continue the spectrum toward large \(\ell\) and compare with the trapped eikonal branch;
8. only after linear spectral stability is established, proceed to nonlinear time-domain evolution.

## 5. Current post-closure verdict

The new calculation gives a materially more favorable picture than the purely geometric statement “there is a stable inner light ring.” The P5 well is shallow enough that the simplest finite-\(\ell\) wave problem does not support a WKB trapped state until approximately \(\ell=256\). The low multipoles that dominate most observational ringdown analyses therefore do not look automatically pathological in this proxy.

This does not close the full tensor–scalar–vector spectral problem, but it removes one of the most immediate low-\(\ell\) concerns and narrows the remaining task to the explicit global representative and its coupled eigenvalue problem.


## 6. Independent Maxwell test-field cross-check

For a test Maxwell field on the same static spherical geometry, the separated radial potential is

$$
V_{\ell}^{(\mathrm{EM})}(r)=f(r)rac{\ell(\ell+1)}{r^2}.
$$

Because $u=r_s/r$, this is exactly

$$
r_s^2 V_{\ell}^{(\mathrm{EM})}=\ell(\ell+1)W(u),
\qquad W(u)=u^2f(u).
$$

Hence the inner-well turning point at the outer-barrier energy is independent of $\ell$, and the maximum WKB action factorizes as

$$
oxed{I_{\max}^{(\mathrm{EM})}(\ell)=\sqrt{\ell(\ell+1)}\,J_{\mathrm{P5}}},
$$

with the purely geometric constant

$$
oxed{J_{\mathrm{P5}}\simeq 0.006114458428}.
$$

The Bohr--Sommerfeld condition $I=\pi(n+1/2)$ therefore yields an explicit threshold law. The first trapped Maxwell ground state occurs at

$$
oxed{\ell_{n=0}=257}.
$$

The first few branch-onset multipoles are

| radial branch $n$ | threshold $\ell$ |
|---:|---:|
| 0 | 257 |
| 1 | 771 |
| 2 | 1284 |
| 3 | 1798 |
| 4 | 2312 |
| 5 | 2826 |
| 6 | 3340 |
| 7 | 3853 |

This is an independent cross-check of the minimally coupled scalar result $\ell_{\rm threshold}\simeq256$. The one-unit difference is caused by the scalar curvature/derivative term absent from the Maxwell potential. The common conclusion is geometric: the P5 stable-light-ring well is too shallow to support a low-multipole semiclassical trapped level in either test-field channel.

## 7. Direct complex-pole shooting: current numerical status

A first Riccati/log-derivative shooting implementation was also attempted for the Maxwell problem, imposing regular-center behavior and an outgoing asymptotic condition. Candidate damped roots with negative imaginary part are readily found at finite outer matching radius, but the first convergence run showed that the naive finite-radius outgoing approximation is not yet sufficiently controlled to label those roots as physical QNMs rather than boundary-placement artifacts. They are therefore **not used as results** here. A production calculation should replace the finite-radius Sommerfeld approximation by a higher-order asymptotic/Jost expansion or compactified exterior integration and demonstrate radial-boundary convergence before quoting complex frequencies.

This failed convergence check is useful: it prevents overclaiming and defines the next numerical implementation step unambiguously.

## 8. Updated post-closure assessment

Two independent geometry-fixed test channels now agree that the stable P5 null well does not generate a semiclassical trapped ground state at the low multipoles relevant to ordinary quadrupolar or octupolar ringdown. The scalar threshold is approximately $\ell=256$ and the exact Maxwell eikonal threshold is $\ell=257$. Stable trapping therefore remains a real asymptotic/high-$\ell$ issue, consistent with the general stable-trapping literature, but the P5 geometry does not exhibit an immediate low-$\ell$ trapped-mode pathology in these probes.

The exact coupled HSVT spectrum remains representative-dependent through the finite-$\ell$ mass/mixing matrix. The next required artifact is therefore one explicit global same-action representative with a center-to-infinity dump of $K,R,G,M$ (and the odd-sector matrices) rather than another local existence argument.
