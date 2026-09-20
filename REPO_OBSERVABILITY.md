# Repository observability

The repository exposes its current scientific/software state through executable views.
Do not infer Absolute Full Closure from filenames or from a single historical report.

## One-command view

```bash
python tools/show_repo.py
```

Machine-readable form:

```bash
python tools/show_repo.py --json
```

The view includes:

- current Absolute-Closure status;
- accepted/rejected/search action members;
- current gate matrix;
- machine-evidence index with SHA-256 hashes;
- repository inventory;
- exact reproduction commands.

## Focused views

```bash
python tools/show_gate_matrix.py
python tools/show_member_matrix.py
python tools/show_evidence_index.py
python tools/show_repo_inventory.py
python tools/show_reproduction_commands.py
```

## Evidence runner

`REPO_EVIDENCE_REGISTRY.json` is the executable registry of audits.  It records the
expected exit code for every audit.  This matters because some rejected-member
witnesses are expected to exit `2`; reproducing that rejection is a registry PASS,
not an Absolute-Closure PASS.

Run the registered scientific evidence:

```bash
python tools/run_all_evidence.py
```

Also run the full software suite and release-manifest verification:

```bash
python tools/run_all_evidence.py --full
```

The runner writes:

- `data/generated/repo_observability/ALL_EVIDENCE_RUN.json`
- `data/generated/repo_observability/ALL_EVIDENCE_RUN.md`

A registry PASS means only that every audit reproduced its declared expected result.
The authoritative Absolute-Closure claim remains `ABSOLUTE_CLOSURE_LEDGER.json`.

## Pre-generated views in checkpoints

The current checkpoint also carries generated, non-authoritative convenience views in
`data/generated/repo_observability/`:

- `REPO_STATUS.txt` / `REPO_STATUS.json`
- `GATE_MATRIX.json`
- `MEMBER_MATRIX.json`
- `EVIDENCE_INDEX.json`
- `INVENTORY.json`
- `REPO_TREE.tsv`
- `CODE_MAP.txt`
- `REPRO_COMMANDS.txt`
- `ALL_EVIDENCE_RUN.json` / `.md`

Regenerate these with the corresponding `tools/show_*` scripts; the authoritative
scientific state remains the underlying ledgers and evidence files.
