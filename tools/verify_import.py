#!/usr/bin/env python3
"""Verify imported bytes; optionally require the complete original handoff."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--require-complete-handoff', action='store_true')
args = p.parse_args()
errors = []
for item in json.loads((root / 'provenance/imported-originals.json').read_text())['files']:
    path = root / item['path']
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
        errors.append(item['path'])
for item in json.loads((root / 'provenance/handoff-manifest-check.json').read_text())['entries']:
    if item['status'] == 'VERIFIED':
        path = root / item['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            errors.append(item['path'])
for path in root.rglob('*.py'):
    if not any(part.startswith('.') for part in path.relative_to(root).parts):
        ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
missing = []
for item in json.loads((root / 'MANIFEST.json').read_text())['files']:
    path = root / item['path']
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
        missing.append(item['path'])
print(f'Import integrity: {"FAIL" if errors else "PASS"}; original handoff entries unavailable: {len(missing)}')
if errors:
    print('\n'.join(errors))
    raise SystemExit(1)
if args.require_complete_handoff and missing:
    print('\n'.join(missing))
    raise SystemExit(4)
