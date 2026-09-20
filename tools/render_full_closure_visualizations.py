#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import math
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
from PIL import Image

ROOT = Path('/mnt/data/ssz_fullclosure_attempt_build/SSZ_FULL_CLOSURE')
ATTEMPT = ROOT / 'data' / 'generated' / 'absolute_attempt_2026-09-19'
OUT = ATTEMPT / 'visualizations'
OUT.mkdir(parents=True, exist_ok=True)

CENTRAL_AUDIT = json.loads((ATTEMPT / 'ELECTRIC_HYBRID_ONSHELL_CENTRAL_AUDIT.json').read_text())
HANDOVER_AUDIT = json.loads((ATTEMPT / 'ACTION_HANDOVER_CONTROL_SPACE_AUDIT.json').read_text())
BG_AUDIT = json.loads((ATTEMPT / 'ELECTRIC_HYBRID_BACKGROUND_EOM_AUDIT.json').read_text())
FULL_REPORT = json.loads((ROOT / 'FULL_CLOSURE_ATTEMPT_REPORT_2026-09-19.json').read_text())
ACTION = pd.read_csv(ATTEMPT / 'ELECTRIC_HYBRID_ONSHELL_CENTRAL_ACTION_JETS.csv').sort_values('u').reset_index(drop=True)
FORTYONE = pd.read_csv(ATTEMPT / 'ELECTRIC_HYBRID_ONSHELL_CENTRAL_41.csv').sort_values('u').reset_index(drop=True)

# Restrict to production bulk where useful
BULK = ACTION[(ACTION['u'] >= 0.62) & (ACTION['u'] <= 0.70)].copy()
BULK41 = FORTYONE[(FORTYONE['u'] >= 0.62) & (FORTYONE['u'] <= 0.70)].copy()


def savefig(path: Path, dpi: int = 180):
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close()


def human_wrap(s, width=50):
    return '\n'.join(textwrap.wrap(str(s), width=width))


# 1) Status dashboard
fig = plt.figure(figsize=(16, 10))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()

ax.text(0.03, 0.96, 'SSZ P5 — Absolute Full Closure Attempt', fontsize=24, weight='bold', va='top')
ax.text(0.03, 0.925, 'Electric-Hybrid branch · visual status snapshot · 2026-09-19', fontsize=12, va='top')

# Header stats
stats = [
    ('Tests', '118 / 118 PASS'),
    ('Central gate', CENTRAL_AUDIT['status']),
    ('Background EOM', BG_AUDIT['status']),
    ('Interfaces', HANDOVER_AUDIT['status']),
    ('Absolute closure', 'NOT CERTIFIED' if not CENTRAL_AUDIT['absolute_full_closure'] else 'CERTIFIED'),
]

x0 = 0.03
w = 0.18
for i, (k, v) in enumerate(stats):
    x = x0 + i * (w + 0.01)
    patch = FancyBboxPatch((x, 0.82), w, 0.08, boxstyle='round,pad=0.01,rounding_size=0.015',
                           linewidth=1.5, edgecolor='black', facecolor='#f4f4f4')
    ax.add_patch(patch)
    ax.text(x + 0.01, 0.875, k, fontsize=11, weight='bold', va='top')
    ax.text(x + 0.01, 0.835, human_wrap(v, 24), fontsize=10, va='top')

# Roadmap
steps = [
    ('Central on-shell K/G', 'PASS'),
    ('Central background EOM', 'PASS*'),
    ('Scalar algebraic certificate', 'PENDING'),
    ('6D compact control rank', 'PASS'),
    ('Joint HSVT interface solve', 'PENDING'),
    ('Direct-41 handover', 'PENDING'),
    ('Global K/R/G/M', 'PENDING'),
    ('Same-operator QNM', 'PENDING'),
    ('Absolute full closure', 'PENDING'),
]

status_colors = {'PASS': '#7ed957', 'PASS*': '#a3e27f', 'PENDING': '#ffd966', 'FAIL': '#ff6b6b'}
y = 0.68
for i, (label, stat) in enumerate(steps):
    x = 0.05 + i * 0.1
    c = status_colors.get(stat, '#dddddd')
    circ = Circle((x, y), 0.018, facecolor=c, edgecolor='black', lw=1.2)
    ax.add_patch(circ)
    if i < len(steps) - 1:
        ax.plot([x + 0.018, x + 0.082], [y, y], color='black', lw=1.5)
    ax.text(x, y - 0.045, human_wrap(label, 16), ha='center', va='top', fontsize=9)
    ax.text(x, y + 0.03, stat, ha='center', va='bottom', fontsize=9, weight='bold')

# Notes panels
panels = [
    ('Central scan',
     f"Rows: {CENTRAL_AUDIT['diagnostics']['rows']}\n"
     f"Production rows: {CENTRAL_AUDIT['diagnostics']['production_rows']}\n"
     f"A2 solve success: {CENTRAL_AUDIT['diagnostics']['a2_solve_success']}\n"
     f"G2 lift q-range: {CENTRAL_AUDIT['diagnostics']['g2_lift']['q_min']:.3f} … {CENTRAL_AUDIT['diagnostics']['g2_lift']['q_max']:.3f}"),
    ('Principal minima',
     f"L=6: min eig(K)={CENTRAL_AUDIT['scans'][0]['min_eig_K']:.6g}\n"
     f"L=1000: min eig(K)={CENTRAL_AUDIT['scans'][-1]['min_eig_K']:.6g}\n"
     f"L=6: min c_r²={CENTRAL_AUDIT['scans'][0]['min_cr2']:.6g}\n"
     f"L=1000: min c_r²={CENTRAL_AUDIT['scans'][-1]['min_cr2']:.6g}"),
    ('Background EOM',
     f"max |E00| = {BG_AUDIT['selected_w9_d8']['current']['E00']['max_abs']:.3e}\n"
     f"max |E11| = {BG_AUDIT['selected_w9_d8']['current']['E11']['max_abs']:.3e}\n"
     f"max |JA| = {BG_AUDIT['selected_w9_d8']['current']['JA']['max_abs']:.3e}\n"
     f"max |Ephi| = {BG_AUDIT['scalar_resolution_statement']['current_direct_Jphi_prime_minus_Pphi_max_abs']:.3e}"),
    ('Handover',
     f"Outer rank: {HANDOVER_AUDIT['compact_six_primitive_rank']['outer']['rank']}/6\n"
     f"Inner rank: {HANDOVER_AUDIT['compact_six_primitive_rank']['inner']['rank']}/6\n"
     f"Outer κ = {HANDOVER_AUDIT['compact_six_primitive_rank']['outer']['condition_normalized']:.3f}\n"
     f"Inner κ = {HANDOVER_AUDIT['compact_six_primitive_rank']['inner']['condition_normalized']:.3f}"),
]

coords = [(0.05, 0.38), (0.52, 0.38), (0.05, 0.13), (0.52, 0.13)]
for (title, text), (x, y) in zip(panels, coords):
    patch = FancyBboxPatch((x, y), 0.4, 0.18, boxstyle='round,pad=0.015,rounding_size=0.02',
                           linewidth=1.5, edgecolor='black', facecolor='#fafafa')
    ax.add_patch(patch)
    ax.text(x + 0.015, y + 0.155, title, fontsize=13, weight='bold', va='top')
    ax.text(x + 0.015, y + 0.125, text, fontsize=11, va='top', family='monospace')

ax.text(0.03, 0.02,
        'PASS* = non-scalar background equations pass at strict tolerance; scalar algebraic certificate still pending.\n'
        'Next gate: solve joint SVT + (a1,c2,c4,F,G,H) compact controls against full background EOM, then re-emit Direct-41.',
        fontsize=10, va='bottom')
savefig(OUT / 'closure_status_dashboard.png', dpi=180)

# 2) Central profiles static
fig, axs = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
axs = axs.ravel()
axs[0].plot(BULK['u'], BULK['A0prime'])
axs[0].set_title('A0prime across central bulk')
axs[0].set_ylabel('A0prime')
axs[1].plot(BULK['u'], BULK['f4'])
axs[1].set_title('f4 across central bulk')
axs[1].set_ylabel('f4')
axs[2].plot(BULK['u'], BULK['v6_A2_resolved'])
axs[2].set_title('v6 (A2 resolved)')
axs[2].set_ylabel('v6')
axs[2].set_xlabel('u')
axs[3].plot(BULK['u'], BULK['G2XX_lift'])
axs[3].set_title('G2XX lift profile')
axs[3].set_ylabel('G2XX_lift')
axs[3].set_xlabel('u')
for ax_ in axs:
    ax_.grid(True, alpha=0.3)
savefig(OUT / 'central_profiles.png', dpi=180)

# 3) Principal scan
scan = pd.DataFrame(CENTRAL_AUDIT['scans'])
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(scan['L'], scan['min_eig_K'], marker='o', label='min eig(K)')
ax.plot(scan['L'], scan['min_cr2'], marker='s', label='min c_r^2')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('L')
ax.set_ylabel('minimum positive gate value')
ax.set_title('Central electric-hybrid principal scan')
ax.grid(True, which='both', alpha=0.3)
ax.legend()
savefig(OUT / 'principal_scan.png', dpi=180)

# 4) Background residuals static
cur = BG_AUDIT['selected_w9_d8']['current']
ref = BG_AUDIT['selected_w9_d8']['archived_exact_reference_direct_replay']
labels = ['E00', 'E11', 'JA', 'Ephi_direct']
current_vals = [cur[k]['max_abs'] for k in labels]
reference_vals = [
    ref['E00']['max_abs'],
    ref['E11']['max_abs'],
    ref['JA']['max_abs'],
    ref['Ephi_direct']['max_abs'],
]
x = np.arange(len(labels))
width = 0.35
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(x - width/2, current_vals, width, label='current central member')
ax.bar(x + width/2, reference_vals, width, label='archived exact replay floor')
ax.axhline(BG_AUDIT['policy_background_residual_abs'], linestyle='--', linewidth=1.5, label='strict tolerance')
ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(['E00', 'E11', 'JA', 'Ephi'])
ax.set_ylabel('max absolute residual')
ax.set_title('Background EOM residual audit')
ax.grid(True, which='both', axis='y', alpha=0.3)
ax.legend()
savefig(OUT / 'background_eom_residuals.png', dpi=180)

# 5) Handover control / interface static
inner = HANDOVER_AUDIT['compact_six_primitive_rank']['inner']
outer = HANDOVER_AUDIT['compact_six_primitive_rank']['outer']
sv_inner = np.array(inner['singular_values_normalized'])
sv_outer = np.array(outer['singular_values_normalized'])
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
axs[0].plot(np.arange(1, len(sv_outer)+1), sv_outer, marker='o', label='outer')
axs[0].plot(np.arange(1, len(sv_inner)+1), sv_inner, marker='s', label='inner')
axs[0].set_yscale('log')
axs[0].set_xlabel('singular-value index')
axs[0].set_ylabel('normalized singular value')
axs[0].set_title('6D compact control singular spectrum')
axs[0].grid(True, which='both', alpha=0.3)
axs[0].legend()

deltas = [
    HANDOVER_AUDIT['endpoint_action_diagnostics']['outer']['delta_v6'],
    HANDOVER_AUDIT['endpoint_action_diagnostics']['inner']['delta_v6'],
]
axs[1].bar(['outer Δv6', 'inner Δv6'], deltas)
axs[1].axhline(0, color='black', lw=1)
axs[1].set_title('Endpoint mismatch before joint solve')
axs[1].set_ylabel('delta v6 = central - old')
axs[1].grid(True, axis='y', alpha=0.3)
savefig(OUT / 'handover_control_space.png', dpi=180)

# 6) Animated central sweep gif
frames_dir = OUT / '_frames_central_sweep'
frames_dir.mkdir(exist_ok=True)
# Sample frames from bulk to keep size manageable
sample_idx = np.linspace(0, len(BULK)-1, 60).astype(int)
for fi, idx in enumerate(sample_idx):
    row = BULK.iloc[idx]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axs = axs.ravel()
    data_map = [
        ('A0prime', BULK['A0prime']),
        ('f4', BULK['f4']),
        ('v6_A2_resolved', BULK['v6_A2_resolved']),
        ('G2XX_lift', BULK['G2XX_lift']),
    ]
    for ax_, (name, series) in zip(axs, data_map):
        ax_.plot(BULK['u'], series)
        ax_.axvline(row['u'], linewidth=2)
        ax_.set_title(name)
        ax_.grid(True, alpha=0.3)
    axs[2].set_xlabel('u')
    axs[3].set_xlabel('u')
    fig.suptitle(f'Central bulk sweep — frame {fi+1}/60 — u={row["u"]:.6f}')
    p = frames_dir / f'frame_{fi:03d}.png'
    plt.tight_layout()
    plt.savefig(p, dpi=110, bbox_inches='tight')
    plt.close()

imgs = [Image.open(p).convert('P', palette=Image.ADAPTIVE) for p in sorted(frames_dir.glob('frame_*.png'))]
if imgs:
    imgs[0].save(OUT / 'central_bulk_sweep.gif', save_all=True, append_images=imgs[1:], duration=110, loop=0)

# 7) Animated roadmap gif
frames_dir2 = OUT / '_frames_roadmap'
frames_dir2.mkdir(exist_ok=True)
roadmap = steps
for upto in range(1, len(roadmap)+1):
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.text(0.03, 0.92, 'Closure roadmap progression', fontsize=20, weight='bold', va='top')
    y = 0.5
    for i, (label, stat) in enumerate(roadmap):
        x = 0.05 + i * 0.1
        if i < upto:
            c = status_colors.get(stat, '#dddddd')
        else:
            c = '#eeeeee'
        circ = Circle((x, y), 0.03, facecolor=c, edgecolor='black', lw=1.2)
        ax.add_patch(circ)
        if i < len(roadmap) - 1:
            ax.plot([x + 0.03, x + 0.07], [y, y], color='black', lw=1.5)
        ax.text(x, y - 0.08, human_wrap(label, 14), ha='center', va='top', fontsize=9)
        if i < upto:
            ax.text(x, y + 0.06, stat, ha='center', va='bottom', fontsize=9, weight='bold')
    current_label, current_stat = roadmap[upto-1]
    ax.text(0.03, 0.12, f'Current stage: {current_label}  [{current_stat}]', fontsize=13)
    ax.text(0.03, 0.06, human_wrap(FULL_REPORT['next_gate'], 95), fontsize=11)
    p = frames_dir2 / f'frame_{upto:03d}.png'
    plt.savefig(p, dpi=120, bbox_inches='tight')
    plt.close()
imgs = [Image.open(p).convert('P', palette=Image.ADAPTIVE) for p in sorted(frames_dir2.glob('frame_*.png'))]
if imgs:
    imgs[0].save(OUT / 'closure_roadmap.gif', save_all=True, append_images=imgs[1:], duration=550, loop=0)

# 8) Summary markdown
summary = ROOT / 'VISUALIZATION_SUMMARY_2026-09-19.md'
summary.write_text(
    '# Visualization summary — 2026-09-19\n\n'
    'This snapshot packages new static and animated visualizations for the electric-hybrid full-closure attempt.\n\n'
    '## Static figures\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/closure_status_dashboard.png\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/central_profiles.png\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/principal_scan.png\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/background_eom_residuals.png\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/handover_control_space.png\n\n'
    '## Animated figures\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/central_bulk_sweep.gif\n'
    '- data/generated/absolute_attempt_2026-09-19/visualizations/closure_roadmap.gif\n\n'
    '## Generator\n'
    '- tools/render_full_closure_visualizations.py\n'
)

print('Rendered visualizations to', OUT)
for p in sorted(OUT.glob('*')):
    if p.is_file():
        print(p.name, p.stat().st_size)
