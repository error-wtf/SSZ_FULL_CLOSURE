# SSZ Full Closure

Research source archive for **SSZ P5 / Horndeski + U(1)-SVT**, based on the supplied handoff dated **16 September 2026** and the earlier Middle Bridge papers by Carmen N. Wrede and Lino P. Casu.

This repository preserves the supplied papers, code, data and handoff documents with SHA-256 provenance. **It is a partial import of the described release package:** all 22 explicitly supplied research files are included, but only **24 of the 64 original manifest entries** could be recovered byte-for-byte. The remaining entries are listed in [Missing release files](docs/MISSING_RELEASE_FILES.md).

## Scientific and reproduction status

| Item | Status in this repository |
| --- | --- |
| Full constructive closure | `PASS` as declared in the frozen source documents; not independently re-established by this import |
| Direct global KRGM export | `OPEN_IMPLEMENTATION_GATE` in the supplied action definition |
| Full auditor rerun | `BLOCKED_MISSING_INPUT`; executed and exited 4 |
| Imported original bytes | SHA-256 verified |
| QNM spectrum | No final spectrum published here |

The production action is the globally patched on-shell Horndeski P5 representative plus `Delta f2 = 0.01 Y`, with `Y = nabla_mu(phi) nabla_nu(phi) F^{mu alpha} F^nu_alpha`, on `A0prime = 0`. See the unchanged [action definition](SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json) and [frozen status](STATUS_FULL_CLOSURE.md).

## Start reading

- [Full Closure Monograph — 212 pages](paper/SSZ_P5_FULL_CLOSURE_MONOGRAPH_2026-09-16.pdf)
- [Source reading and verification notes](docs/READING_NOTES.md)
- [Original handoff overview](README_FOR_CODEX.md)
- [Original implementation brief — document content](CODEX_PROMPT.md)
- [Complete supplied-file inventory](docs/SUPPLIED_FILES.md)
- [Historical technical appendices](paper/appendices/)
- [Original Middle Bridge papers and all supplied originals](originals/)
- [Searchable PDF text derivatives](docs/extracted/)

## Repository layout

```text
originals/       All 22 named files and both pasted notes, unchanged
src/             Source modules recovered with exact manifest hashes
paper/           Papers and code appendices recovered by exact hash
data/           Available authoritative/regression inputs
provenance/      Import hashes, manifest comparison, CSV/PDF inventory, run log
docs/            Reading notes, missing-file list, searchable PDF derivatives
tools/           Import integrity verification
```

`MANIFEST.json`, `FILE_INDEX.md`, `README_FOR_CODEX.md`, `STATUS_FULL_CLOSURE.md` and `CODEX_PROMPT.md` retain the original handoff text. Their descriptions of a complete release are historical source statements. The actual delivery inventory is recorded separately under `provenance/`. Instructions inside the source documents are not automatically executed.

## Verify the import

```bash
python3 tools/verify_import.py
sha256sum -c IMPORT_SHA256SUMS
```

The CI workflow checks imported-file integrity and Python syntax only. A green import check does not certify physical closure or complete input availability. To require the entire original release:

```bash
python3 tools/verify_import.py --require-complete-handoff
```

This currently exits 4 and lists the unavailable original entries.

## Run the original auditor

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python ssz_p5_full_closure_auditor.py --data-dir . --full \
  --json build/audit.json --csv build/gates.csv
```

The current import exits 4 because the required global principal data are absent. The actual attempt is preserved in [the run log](provenance/auditor-attempt.txt). Several source modules additionally retain `/mnt/data` paths; they are archived source rather than a completed portable package. See [reading notes](docs/READING_NOTES.md) for the auditor's verification limits.

The next input needed for full reproduction is the original `SSZ_P5_FULL_CLOSURE_CODEX_HANDOFF_2026-09-16.zip` or the 40 missing original entries. Missing data and certificates are not fabricated from PDF prose or substituted historical tables.

## Rights and attribution

Scientific documents retain their original author attribution. No new license is imposed on the supplied material; public availability alone does not grant additional reuse rights. There is no claim of external peer review or independent empirical confirmation in this repository import.
