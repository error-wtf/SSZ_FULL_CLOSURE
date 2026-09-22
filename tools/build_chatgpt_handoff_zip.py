#!/usr/bin/env python3
"""Assemble the complete self-contained ChatGPT handoff ZIP (2026-09-22)."""
from __future__ import annotations
import hashlib, json, os, shutil, subprocess, sys, zipfile
from pathlib import Path

ROOT = Path('/home/error/ssz-full-closure')
STAGE = Path('/tmp/opencode/handoff_stage')
ZIP_PATH = Path('/home/error/SSZ_P5_2026-09-22_EQ85_Q2_COMPLETE_HANDOFF.zip')

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

# ---------- curated tree ----------
COPY_DIRS = [
    ('docs/canonical', 'docs/canonical'),
    ('docs/phase1', 'docs/phase1'),
    ('src/ssz_p5', 'src/ssz_p5'),
    ('data/production', 'data/production'),
]
COPY_FILES = [
    'README.md', 'CURRENT_WORKING_STATE_2026-09-21_EQ47.md',
    'FULL_CLOSURE_WORKING_STATUS_2026-09-21_EQ47.json', 'MANIFEST.json', 'SHA256SUMS',
    'requirements.lock',
    'tools/context_guard.py', 'tools/release_manifest.py',
    'tools/audit_dense_angular_direct.py', 'tools/refreeze_dense_baseline.py',
    'tools/build_prospective_member_v2.py', 'tools/phase1_freeze_hashes.py',
    'tools/advance_transition_dense_h_continuation.py',
    'tools/audit_angular_eq47_projection.py', 'tools/run_angular_universal_pipeline.py',
    'tools/audit_zero_vector_light_ring.py',
    'archive/full_working_snapshot/build_ssz_p5_selected_41stream_v2_2026-09-16.py',
    'archive/full_working_snapshot/ssz_p5_holonomic_a5_closure_2026-09-16.py',
    'archive/initial_import/src/ssz_p5_higher_jet_closure_2026-09-16.py',
    'archive/initial_import/originals/SSZ_P5_HSVT_ACTION_MEMBER_2026-09-16.json',
    'archive/initial_import/originals/SSZ_P5_SHARED_BASELINE_ASSEMBLY_VALIDATION_2026-09-16.md',
    'data/production/ssz_p5_horndeski_carrier_through_light_rings_to_core_2026-09-12.csv',
    'data/production/ssz_p5_integrable_svt_unreduced_even_coefficients_2026-09-12.csv',
    'data/production/ssz_p5_integrable_full_svt_lobe_ZK_exact_regression_2026-09-13.csv',
    'data/production/ssz_p5_F3_outer_same_action_SELECTED_41of41_2026-09-15.csv',
    'data/production/ssz_p5_F3_inner_same_action_SELECTED_41of41_CANDIDATE_2026-09-15.csv',
    'data/production/ssz_p5_F3_core_SELECTED_41of41_2026-09-15.csv',
    'data/authoritative/ssz_p5_F2_horndeski_carrier_unreduced_39of41_CORRECTED_2026-09-14.csv',
    'data/diagnostic/ssz_p5_SELECTED_41STREAM_V2_2026-09-16.csv',
    'data/generated/absolute_attempt_2026-09-19/ZERO_VECTOR_LIGHT_RING_AUDIT.json',
    'data/generated/angular_eq47_projection_2026-09-21/EQ47_PROJECTION_AUDIT.json',
    'data/generated/phase2_q2/Q2_INSTRUMENTATION_NOOP_AUDIT.json',
    'data/generated/phase2_q2/Q2_SPLIT_RECOVERY_AUDIT.json',
    'data/generated/phase2_q2/Q2_PRODUCTION_MEMBER_STATUS.json',
    'data/generated/phase2_q2/Q2_SAME_ACTION_COMPATIBILITY_AUDIT.json',
    'data/generated/phase2_q2/Q2_CHAIN_DECOMPOSITION_AUDIT.json',
    'data/generated/phase2_q2/Q2_INTERFACE_TEST.json',
    'data/generated/phase2_q2/Q2_ARCHIVE_RECOVERY_LEDGER.json',
    'data/generated/phase2_q2/HISTORICAL_GENERATOR_REPLAY_AUDIT.json',
    'data/generated/phase2_forensics/HISTORICAL_MEMBER_PROVENANCE_COMPARISON.json',
    'data/generated/phase2/PHASE2_UNBLINDING_LOG.json',
    'data/generated/phase2/PHASE2_RESULTS_C1.json', 'data/generated/phase2/PHASE2_RESULTS_C3.json',
    'data/generated/phase2/PHASE2_K_GATES_C1_C3.json', 'data/generated/phase2/PHASE2_C3_GRAM.json',
    'data/generated/phase2/PHASE2_BACKGROUND_RESIDUALS.json',
    'data/generated/phase2_q2/Q2_SHARED_SECTOR_PROVENANCE.json',
]

def build():
    if STAGE.exists(): shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    for src, dst in COPY_DIRS:
        s, d = ROOT/src, STAGE/dst
        if s.exists(): shutil.copytree(s, d, ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '.venv'))
    for f in COPY_FILES:
        s = ROOT/f
        if s.exists():
            d = STAGE/f; d.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(s, d)
    # pruned heavy generated files: keep only key dense-angular evidence
    dd = STAGE/'data/generated/dense_angular_direct_2026-09-22'
    dd.mkdir(parents=True, exist_ok=True)
    for f in ('DENSE_ANGULAR_DIRECT.csv','DENSE_ANGULAR_DIRECT_SUMMARY.json','DENSE_DIRECT_FINITE_L_CROSS.csv'):
        s = ROOT/'data/generated/dense_angular_direct_2026-09-22'/f
        if s.exists(): shutil.copy2(s, dd/f)
    hist = STAGE/'historical_sources'
    hist.mkdir(exist_ok=True)
    g = ROOT/'archive/full_working_snapshot/build_ssz_p5_selected_41stream_v2_2026-09-16.py'
    shutil.copy2(g, hist/'build_ssz_p5_selected_41stream_v2_2026-09-16.py')
    (hist/'SOURCE_MAP.txt').write_text(
        'build_ssz_p5_selected_41stream_v2_2026-09-16.py\n'
        '  repo path: archive/full_working_snapshot/\n'
        '  entered git: commit 9ef9247 (2026-09-17 initial import); work dated 2026-09-16\n'
        '  era: ERA_ZK_SEP13_16 (see EMITTER_ERA_LEDGER.csv)\n'
        '  role: ORIGINAL GENERATOR of the Sep-16 selected 41-stream (region-priority\n'
        '  stitcher; v5=c3=e3=0 selection; corrected holonomic a5; INTERFACE_JUMPS export)\n'
        '  sha256: '+sha(g)+'\n')
    return len(list(STAGE.rglob('*')))

def security_scan():
    pats = [b'ghp_', b'github_pat_', b'BEGIN OPENSSH PRIVATE KEY', b'BEGIN RSA PRIVATE KEY', b'x-access-token:']
    hits = []
    for p in STAGE.rglob('*'):
        if p.is_file() and p.stat().st_size < 5_000_000:
            try: blob = p.read_bytes()
            except Exception: continue
            for pat in pats:
                if pat in blob: hits.append((str(p.relative_to(STAGE)), pat.decode()))
    (STAGE/'SECURITY_SCAN.json').write_text(json.dumps({'secrets_found': len(hits), 'patterns': [p for _,p in hits]}, indent=1))
    return hits

def write_handoff():
    md = ROOT/'docs'/'handoff'; md.mkdir(parents=True, exist_ok=True)
    h = md/'CHATGPT_HANDOFF_2026-09-22.md'
    h.write_text(Path('/tmp/opencode/handoff_md.md').read_text())
    shutil.copy2(h, STAGE/'CHATGPT_HANDOFF_2026-09-22.md')
    for f in ('REPRODUCE_HANDOFF.md',):
        shutil.copy2(md/f, STAGE/f)

def main():
    n = build()
    hits = security_scan()
    if hits:
        print('SECRETS FOUND — ABORT:', hits); sys.exit(1)
    write_handoff()
    # GIT_STATE
    gs = subprocess.run(['git','-C',str(ROOT),'log','--oneline','-12'],capture_output=True,text=True).stdout
    st = subprocess.run(['git','-C',str(ROOT),'status','--short'],capture_output=True,text=True).stdout
    (STAGE/'GIT_STATE.txt').write_text(
        'branch: main\ncommit: '+(ROOT/'.git'/'refs'/'heads'/'main').read_text().strip()+
        '\nremote: https://github.com/error-wtf/SSZ_FULL_CLOSURE.git\n\nrecent commits:\n'+gs+'\nstatus --short:\n'+(st or '(clean)'))
    # manifest + sha256sums
    files = sorted(p for p in STAGE.rglob('*') if p.is_file())
    # self-exclusion convention: MANIFEST.json cannot contain its own hash; it is
    # verified via ZIP integrity + SHA256SUMS.txt coverage of every other file
    files = [p for p in files if p.name not in ('MANIFEST.json', 'SHA256SUMS.txt')]
    manifest = [{'path': str(p.relative_to(STAGE)), 'size': p.stat().st_size, 'sha256': sha(p)} for p in files]
    (STAGE/'MANIFEST.json').write_text(json.dumps({'files': manifest}, indent=1))
    (STAGE/'SHA256SUMS.txt').write_text(''.join(f"{m['sha256']}  {m['path']}\n" for m in manifest))
    # zip
    if ZIP_PATH.exists(): ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(STAGE.rglob('*')):
            if p.is_file(): z.write(p, p.relative_to(STAGE))
    zsha = sha(ZIP_PATH)
    rel = {'ZIP_PATH': str(ZIP_PATH), 'ZIP_SHA256': zsha, 'FILE_COUNT': len(manifest),
           'ZIP_SIZE_BYTES': ZIP_PATH.stat().st_size,
           'EQ85_STATUS': 'BLOCKED_PENDING_PRIMARY_SOURCE_RETRIEVAL (replay spec included)',
           'EPSILON_Y_MEMBER_STATUS': 'REPLAY_PENDING (neither active nor superseded-as-replayed)',
           'Q2_STATUS': 'OPEN — electric chain; paused on Eq85 replay + C_shared provenance (historical track)',
           'ABSOLUTE_FULL_CLOSURE_STATUS': False}
    Path('/home/error/HANDOFF_RELEASE.json').write_text(json.dumps(rel, indent=1))
    print(json.dumps(rel, indent=1))

if __name__ == '__main__':
    main()
