# Reproduce and interpret the current repository

Supported environment: CPython 3.14 on Linux with the pinned dependency lock.

```bash
git clone https://github.com/error-wtf/SSZ_FULL_CLOSURE.git
cd SSZ_FULL_CLOSURE
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
pytest -q
python tools/release_manifest.py --check
python ssz_p5_full_pipeline.py --strict
```

## Expected results in this implementation checkpoint

- Software tests and file-integrity verification succeed.
- `python ssz_p5_full_pipeline.py` reproduces the archive/audit checks.
- `python ssz_p5_full_pipeline.py --strict` exits **2** because the complete
  direct-production and coupled spectral certificates have not been generated,
  and the prescribed carrier fails a necessary on-shell background identity.

The strict failure is intentional evidence of an unfinished production chain;
it must not be ignored when deciding whether to publish an absolute-closure
release. See [implementation findings](docs/DIRECT_IMPLEMENTATION_FINDINGS.md)
for the exact source-to-emitter mismatch. The old archive-only meaning of
`--strict` has been removed.

## What strict mode requires

Strict mode runs the software tests, current manifest/hash checks, the existing
full closure auditor, blacklist checks and archived spectral-record checks. It
additionally requires valid direct coefficient/matrix products for every default
L, verifies their finite-L structure and stability, reproduces the matrices
through the profile reducer, binds spectral results to those same operator
hashes, and requires convergence and absolute-closure evidence.

Runtime output is written to `FULL_PIPELINE_REPORT.json`,
`FULL_PIPELINE_REPORT.md`, `build/audit.json` and `build/gates.csv`.
These outputs and caches are excluded from the immutable release inventory.
`--strict --skip-tests` is rejected.

No archive-only PASS is a final coupled-spectrum certificate. No full-closure
v1.0.0 archive is created while strict mode fails.

## Reproduce the source contract finding

```bash
python tools/audit_production_sources.py
pytest -q tests/unit/test_constraint_stationarity.py tests/negative/test_mh_action_domain.py
```

The JSON source report records the input hashes, sampled radial coverage and
action jets outside the restricted emitter's supported domain. The tests verify
the repaired H1 source against the unreduced Euler equations and ensure that
unsupported Horndeski jets are rejected instead of silently erased.

## After intentional source changes

Regenerate the manifest and SHA256 inventory only after the final file edits:

```bash
python tools/release_manifest.py
python tools/release_manifest.py --check
```

Then test an exported clean tree in a fresh environment. Archived originals,
including the initial import and supplied snapshot manifests, remain immutable.

## Reproduce the frozen-carrier contradiction

```bash
python tools/audit_frozen_onshell.py
pytest -q tests/unit/test_onshell_identity.py tests/regression/test_frozen_onshell_diagnostic.py
```

The diagnostic deliberately exits **2** for the supplied carrier and writes
`build/FROZEN_ONSHELL_IDENTITY.json`. The seven test cases pass because they
verify the identity on analytic controls and verify detection of the frozen
input failure. Passing these tests is not a physical background PASS.

A committed snapshot with source hashes is in
[data/diagnostic/FROZEN_ONSHELL_IDENTITY.json](data/diagnostic/FROZEN_ONSHELL_IDENTITY.json).
[The derivation](docs/FROZEN_CARRIER_ONSHELL_CONTRADICTION.md) identifies the exact
contradictory assumptions, radius, equation, and limits of the numerical check.
