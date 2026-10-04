# Real Spectroscopy & Certified Modes — Execution Plan (2026-10-04)

## Source of truth

Do **not** use an old ZIP as the working source. The repository state is now:

- `main` = `df09672eaff34f1e32c13341181bfd652976d4e6`
- `qnm-global-final-20261002` is **149 commits ahead of main and 0 behind**
- `spectroscopy-real-data-20261004` is forked from `qnm-global-final-20261002`

The spectroscopy/QNM machinery therefore lives on the spectroscopy branch lineage, not on the current `main` snapshot.

## What is already established

1. Fixed-observable local SSZ spectroscopy exists and includes degeneracy-aware generalized-eigenmode tracking, residue weights and resolvent checks.
2. Weisz-1978 FK spectra have an independent numerical reproduction path.
3. A global coupled QNM diagnostic has been executed, but a **physical coupled QNM claim is not allowed** because the global finite-l operator fails mandatory health gates.
4. The global blocker was localized away from interfaces into the punctured-H core. The archived positive core `K_scalar` oracle is not reproducible from the archived primitive profiles and a stencil-robust negative pocket remains.
5. Therefore: local spectroscopy is scientifically usable; global coupled QNM frequencies are not yet certifiable.

## Definition of a “perfect mode”

A mode is not “perfect” because one eigensolver returns a sharp eigenvalue. A publishable/certified mode must satisfy all of the following:

- same locked operator and member hash;
- positive kinetic and required hyperbolicity/stability gates on the full domain;
- convergence under radial resolution and derivative stencil;
- convergence under independent boundary implementations (Jost/ECS/compactified outgoing);
- stability against window/domain truncation;
- gauge/basis invariance of observable residues;
- degeneracy-aware subspace tracking;
- stable left/right (biorthogonal) eigenvectors for the non-Hermitian problem;
- bounded condition number / pseudospectral sensitivity;
- stable pole location and residue in the retarded Green function;
- successful synthetic-injection recovery.

## Execution ladder

### A. Paper truth set — Weisz

Reproduce every directly testable quantity from the equations/parameter sets in the Weisz papers:

- normal-mode frequencies;
- q=0 optical residues;
- effective visible-mode count;
- high-order near-degenerate residue suppression;
- line-strength hierarchy;
- broadened Green-function spectral curves.

No SSZ parameters enter this stage.

### B. Paper truth set — Zhang/Li curved/quasiperiodic model

Reconstruct the published Hamiltonian exactly from the paper:

- position-dependent hopping / curvature term;
- quasiperiodic onsite term;
- the same irrational/rational approximants and boundary conditions;
- eigenenergies and eigenvectors;
- inverse participation ratio (IPR);
- local density of states (LDOS);
- mobility edge / localized–extended phase boundary;
- wave-packet propagation used as a negative/positive localization control.

The comparison to SSZ is **structural**, not an identity claim.

### C. SSZ local observable spectroscopy

On the registered healthy SSZ window only, compute for fixed observables q:

[
Gv_n=\omega_n^2Kv_n,qquad
Z_n^{(q)}=\frac{|q^T v_n|^2}{2\omega_n},
]

and the direct observable resolvent

[
A_q(\omega,r)=-\frac1\pi\operatorname{Im}
q^T[G-(\omega+i\eta)^2K]^{-1}q.
]

Export raw, machine-readable:

- `omega_n(r,L)`;
- raw and normalized `Z_n(r,L,q)`;
- spectral curves `A_q(omega,r)`;
- mode-tracking overlaps;
- degeneracy clusters;
- effective mode number / entropy;
- radius at every robust residue-order inversion.

### D. Cross-model falsification tests

Do not compare “pictures”. Compare dimensionless observables:

- top-1 / top-2 / top-3 residue fraction;
- effective mode number;
- spectral entropy;
- residue dynamic range;
- near-degenerate frequency separation versus residue separation;
- LDOS concentration;
- IPR and participation ratio where defined;
- rank inversions with radius/position;
- robustness against artificial broadening eta.

The question is not “does SSZ look like Weisz/Li?” but:

> Which spectral-selection invariants are shared, and which falsify the analogy?

### E. Real observational spectroscopy

Use only public, independently archived observations. Candidate empirical layers:

1. **NICER X-ray binaries** — energy-resolved timing/spectral data, QPO/rms/phase-lag structure.
2. **ESO Galactic-center spectroscopy** — redshift/orbital spectral observables around Sgr A*.
3. **ALMA spectral cubes** — frequency-resolved spatial observables for compact-source/accretion environments.
4. **GW public ringdown data** — direct QNM-frequency/damping comparison once the coupled global SSZ operator is certified.

For each dataset:

- download observation metadata first;
- freeze observation IDs and calibration versions;
- keep background/calibration nuisance parameters separate from SSZ parameters;
- forward-model the observable in detector space;
- perform blind/posterior predictive tests;
- compare GR/null and SSZ with the **same nuisance model**.

No fitting of SSZ free functions to each spectrum is allowed.

### F. Global coupled SSZ modes

This stage remains gated.

Current blocker:
- the reconstructed global finite-l operator has a core instability / inconsistent archived positive `K_scalar` oracle;
- interface amplification is not the dominant cause.

Required sequence:

1. reconstruct the core principal oracle directly from covariant action jets;
2. decide whether the x~0.30 negative pocket is physical or an obsolete-oracle inconsistency;
3. regenerate a single direct-global 41-slot same-action stream;
4. require full-domain K/radial/angular/pivot health;
5. build compactified Jost and independent ECS solvers;
6. keep only poles common to both solvers and converged under resolution;
7. compute retarded-Green residues;
8. only then compare physical QNM frequencies with real ringdown data.

## Stop rules

A result is blocked, not “almost passed”, when:

- K <= 0 on the claimed domain;
- radial/angular characteristic gate fails;
- a constraint pivot vanishes;
- a pole moves materially with numerical boundary settings;
- mode identity depends on arbitrary eigenvector rotations inside a degenerate subspace;
- observed-data improvement disappears under the same nuisance model applied to the null model.

## Immediate work order

1. Validate the spectroscopy branch as source-of-truth.
2. Complete the Zhang/Li numerical reproduction next to the existing Weisz reproduction.
3. Regenerate raw SSZ local resolvent/mode tables on the same branch.
4. Build one common comparison table for Weisz, Li and SSZ.
5. Query public observation archives and pin concrete observation IDs.
6. In parallel, resolve the core principal-oracle inconsistency; only this unlocks certified global QNM modes.

