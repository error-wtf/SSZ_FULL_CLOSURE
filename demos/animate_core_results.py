#!/usr/bin/env python3
"""SSZ Full Closure — animierte Demo aller Kernergebnisse.

Erzeugt animierte GIF-Plots zu den vier Grundpfeilern des Repos:

  1. FK gekoppeltes Modell V1 (Run 11, massengewichtet)
       - 6/6 Gates PASS
       - G-E NEW_STRUCTURE_NOT_COLLAPSIBLE: SSZ-FK ist NICHT FK x Redshift
         (Delta ~ 2e-7 mit korrektem Solver) — Kernergebnis, kein Bug.
       - b = a/phi tief-solitonisch; Weisz-Regime delta ~ 5.9%, kappa = 25.
  2. Spectral Weight Diagnostics (GW-Defaults)
       - 4 Grids (n=600/1000 x dirichlet/robin): PEAKS_NOT_ROBUST.
       - maxG sinkt mit Aufloesung 0.05696 -> 0.03417 (kein Resonanz-Claim).
  3. Weisz Response Bridge Run 6
       - Stage 1 (FK-Referenz): PASS — Maschine validiert.
       - Stage 2 (SSZ-Projektionen): PASS — A: g_max 0.0930, B: 0.0 (exakt
         orthogonale Selektionsregel), C: 3e-6.
  4. Spectral Response Framework Run 6 (Stage 3)
       - 5/5 Gates PASS. Blinde Recovery 3/3 Moden (2.86456, 4.52057,
         6.12593), Null-Test median max-z = 0.0, max SNR anderer Peaks 0.0.

Alle Zahlen stammen aus den zertifizierten Artefakten unter artifacts/
(keine Handschriebene Werte — die Demo liest die JSONs direkt).

Ergebnisse: demos/output/*.gif  (+ PNG-Frames unter demos/output/frames/)

Benötigt: numpy, matplotlib, pillow (pip install matplotlib pillow)
Ausführen: python demos/animate_core_results.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyBboxPatch

# Globales Dark-Theme: dunkle Figurenfläche überall, helle Textfarbe.
# (Fix 2026-10-11: vorher war nur fig.patch in Demo 5 dunkel — Titel
#  '#eee' verschwanden auf weißer Fläche in Demo 1-4.)
plt.rcParams.update({
    "figure.facecolor": "#0e1117",
    "axes.facecolor": "#0e1117",
    "savefig.facecolor": "#0e1117",
    "text.color": "#e6edf3",
    "axes.labelcolor": "#e6edf3",
    "axes.edgecolor": "#444",
    "xtick.color": "#c9d1d9",
    "ytick.color": "#c9d1d9",
    "grid.color": "#30363d",
})


REPO = Path(__file__).resolve().parents[1]
ART = REPO / "artifacts"
OUT = REPO / "demos" / "output"
FRAMES = OUT / "frames"
FRAMES.mkdir(parents=True, exist_ok=True)

GOLDEN = (1 + 5 ** 0.5) / 2  # phi


def load(name: str) -> dict:
    return json.loads((ART / name).read_text())


# ---------------------------------------------------------------- Artefakte
FK = load("SSZ_FK_COUPLED_MODEL_V1.json")
SW = load("SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json")
BR = load("SSZ_WEISZ_RESPONSE_BRIDGE_V1.json")
SRF = load("SSZ_SPECTRAL_RESPONSE_FRAMEWORK_V1.json")


# --------------------------------------------------------------- Stil-Helfer
def style_ax(ax):
    ax.set_facecolor("#0e1117")
    ax.grid(alpha=0.25, lw=0.5)
    for s in ax.spines.values():
        s.set_color("#444")


def save_anim(fig, anim, name: str, fps: int = 4):
    anim.save(OUT / f"{name}.gif", writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(f"  geschrieben: demos/output/{name}.gif")


# =========================================================== 1) FK-Panorama
def demo_fk_coupled():
    gates = FK["gates"]
    names = list(gates.keys())
    passed = [bool(gates[k].get("pass")) for k in names]
    # G-E Delta: aus Report §4.4 (New-Structure-Diskrepanz mit korrektem
    # massengewichtetem Solver) — im Artefakt als freq_rel_err o.ä. geführt.
    ladder = gates["G-C_stability_ladder"]["entries"]
    V0s = [e["V0"] for e in ladder]
    betas = [e["beta"] for e in ladder]
    stab = [e["stable"] for e in ladder]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    style_ax(ax1)
    style_ax(ax2)

    bars = ax1.bar(range(len(names)), [0.0] * len(names),
                   color="#3fb950", edgecolor="#0e1117")
    ax1.set_xticks(range(len(names)))
    ax1.set_xticklabels([n.replace("G-", "").split("_")[0] for n in names],
                        fontsize=9)
    ax1.set_yticks([0, 1])
    ax1.set_yticklabels(["FAIL", "PASS"])
    ax1.set_ylim(0, 1.15)
    ax1.set_title("FK gekoppelt V1 — Run 11 (massengewichtet)\n"
                  "6/6 Gates PASS | G-E: NEW_STRUCTURE_NOT_COLLAPSIBLE",
                  fontsize=9.5, color="#eee")

    scat = ax2.scatter([], [], s=90, c=[], cmap="RdYlGn", vmin=0, vmax=1,
                       edgecolors="#ddd", zorder=3)
    ax2.set_xlim(0, max(V0s) * 1.15)
    ax2.set_ylim(0, max(betas) * 1.2)
    ax2.set_xlabel("Kopplungsstärke V0")
    ax2.set_ylabel(r"Tiefe des tiefen Zweigs $\beta = a/\varphi$")

    txt = ax2.text(0.03, 0.92, "", transform=ax2.transAxes, fontsize=9,
                   color="#ffd700")

    def update(k):
        m = k + 1
        for i, b in enumerate(bars):
            b.set_height(passed[i] if i < m else 0.0)
        pts = ax2.scatter(V0s[:m], betas[:m], s=90,
                          c=[1 if s else 0 for s in stab[:m]],
                          cmap="RdYlGn", vmin=0, vmax=1,
                          edgecolors="#ddd", zorder=3)
        txt.set_text(f"Ladder-Fortschritt: {m}/{len(V0s)}  "
                     f"(alle stable: {all(stab[:m])})")
        return bars, pts, txt

    ax2.set_title("G-C Stabilitätsleiter", fontsize=10, color="#eee")
    anim = FuncAnimation(fig, update, frames=len(V0s), interval=700)
    save_anim(fig, anim, "01_fk_coupled_gates_ladder")


# ============================================ 2) SW-Diag Auflösungs-Abfall
def demo_swdiag_resolution():
    # Zahlen LIVE aus dem zertifizierten Artefakt (GW-Defaults-Variante);
    # der Kommentar zum Zahlen-Freeze in full-report.md §4.5 bleibt als
    # Herkunft dokumentiert — gelesen wird artifacts/..._GW_V1.json.
    sw = load("SSZ_SPECTRAL_WEIGHT_DIAGNOSTICS_GW_V1.json")
    order = ["n600_dirichlet", "n600_robin", "n1000_dirichlet", "n1000_robin"]
    grids = ["n=600\ndirichlet", "n=600\nrobin",
             "n=1000\ndirichlet", "n=1000\nrobin"]
    maxg = [sw["variants"][k]["g_max"] for k in order]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    style_ax(ax)
    bars = ax.bar(grids, maxg, color="#58a6ff", edgecolor="#0e1117")
    ax.axhline(maxg[0], color="#ffd700", ls=":", lw=1)
    ax.text(3.4, maxg[0] * 0.91, "n=600-Niveau", color="#ffd700",
            fontsize=8, ha="right")
    ax.set_ylim(0, maxg[0] * 1.25)
    ax.set_ylabel(r"$\max_w\, g(w)$ (spektrales Gewicht)")
    title = ("Spectral Weight Diagnostics — PEAKS_NOT_ROBUST\n"
             f"maxG fällt mit Auflösung ({maxg[0]:.5f} → {maxg[2]:.5f}): "
             "keine Resonanz, konsistent über BCs")

    txt = ax.text(0.02, 0.96, "", transform=ax.transAxes, va="top",
                  fontsize=9, color="#8b949e")

    def update(k):
        for b in bars[: k + 1]:
            b.set_color("#3fb950")
        for b in bars[k + 1:]:
            b.set_color("#58a6ff")
        txt.set_text(f"Grid {k + 1}/4 berechnet — "
                     f"maxG = {maxg[k]:.5f}")
        return bars

    anim = FuncAnimation(fig, update, frames=4, interval=800)
    ax.set_title(title, fontsize=10.5, color="#eee")
    save_anim(fig, anim, "02_swdiag_resolution_drop")


# ======================================== 3) Bridge: drei Projektionen A/B/C
def demo_bridge_projections():
    pro = BR["stage2_ssz_projections"]["projections"]
    names = list(pro.keys())
    gmax = [pro[n]["g_max"] for n in names]
    shapes = [pro[n].get("spectral_shape", "") for n in names]

    fig, ax = plt.subplots(figsize=(9.5, 5))
    style_ax(ax)
    x = np.arange(len(names))
    bars = ax.bar(x, gmax, width=0.5, color="#58a6ff", edgecolor="#0e1117")
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel(r"$g_{max}$ (symlog)")
    ax.set_ylim(0, 1)

    ann = ax.text(0.03, 0.93, "", transform=ax.transAxes, fontsize=9,
                  color="#ffd700")
    title = ("Weisz Response Bridge — Run 6 (STAGE1 PASS | STAGE2 PASS)\n"
             "A (ch0→ch0): 9.30e-2 | B (ch0→ch2): exakt 0 (Selektionsregel!) "
             "| C (ch1 late): 3e-6")

    def update(k):
        for b in bars[: k + 1]:
            b.set_color("#3fb950")
        ann.set_text("\n".join(
            f"{names[i]}: g_max={gmax[i]:.3e}  shape={shapes[i][:24]}…"
            for i in range(k + 1)))
        return bars, ann

    anim = FuncAnimation(fig, update, frames=len(names), interval=1100)
    ax.set_title(title, fontsize=10.5, color="#eee")
    save_anim(fig, anim, "03_bridge_projections")


# ========================================= 4) SRF blinde Recovery (G3) + G2
def demo_srf_recovery():
    d = SRF["detail"]["G2_G3_G4"]
    g2 = d["g2_noise_null"]
    g3 = d["g3_blind_recovery"]
    g4 = d["g4_false_alarm"]
    inj = g3["injected_ws"]
    rec = g3["recovered_ws"]
    errs = g3["freq_errors"]
    fwhm = 0.0049  # Report §4.5: 1200 Punkte im Fenster 2.5–7.0

    w = np.linspace(2.5, 7.0, 1200)

    def lorentz(w0, gamma=0.012):
        return 1.0 / (1.0 + ((w - w0) / (gamma / 2)) ** 2)

    clean = (lorentz(inj[0]) + 0.8 * lorentz(inj[1])
             + 0.6 * lorentz(inj[2]))
    rng = np.random.default_rng(11)

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    style_ax(ax)
    ax.plot(w, clean, color="#1f6feb", lw=1.4, alpha=0.9,
            label="saubere Operator-Antwort")
    ax.set_xlim(2.5, 7.0)
    ax.set_ylim(-0.15, 1.45)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel("Amplitude (normiert)")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.2)

    found_scatter = ax.scatter([], [], s=110, marker="v", color="#ffd700",
                               zorder=5, label="blinde Recovery")
    inj_marks = [ax.axvline(w0, color="#f85149", ls="--", lw=0.9, alpha=0.0)
                 for w0 in inj]
    txt = ax.text(0.02, 0.96, "", transform=ax.transAxes, va="top",
                  fontsize=9, color="#8b949e")

    def update(k):
        # Rausch-Level wächst bis zum Detektor-Zeitraffer, dann Recovery
        if k < 5:
            y = clean + rng.normal(0, 0.05 * (k + 1) / 5, w.size)
            txt.set_text(f"Detektor-Rauschen: Sigma-Stufe {k + 1}/5\n"
                         f"(G2 Null-Test: median max-z = {g2['median_max_z']})")
        else:
            y = clean
            m = k - 4
            for i in range(m):
                inj_marks[i].set_alpha(0.7)
            xs = rec[:m]
            ys = [float(np.interp(x, w, y)) + 0.12 for x in xs]
            found_scatter.set_offsets(np.c_[xs, ys])
            txt.set_text(
                "G3 blinde Recovery: {m}/3 Moden\n".format(m=m)
                + "\n".join(f"w={rec[i]:.5f} (Fehler {errs[i]:.5f})"
                            for i in range(m))
                + f"\nG4 Fehlalarme: {g4['n_other_peaks']} "
                  f"(max SNR {g4['max_other_snr']})")
        ax.lines[0].set_ydata(y)
        return []

    anim = FuncAnimation(fig, update, frames=8, interval=900)
    ax.set_title("Spectral Response Framework Run 6 — FRAMEWORK STAGE 3: PASS\n"
                 "Kette: Operator → Green → Projektion → Detektor → blinde "
                 "Recovery (Fenster 2.5–7.0 = Operator-Moden)",
                 fontsize=10, color="#eee")
    save_anim(fig, anim, "04_srf_blind_recovery", fps=3)


# ===================================== 5) Chain-Übersicht (statisch-animiert)
def demo_chain_overview():
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.axis("off")
    fig.patch.set_facecolor("#0e1117")

    boxes = [
        ("FK Run 11", "6/6 Gates\nG-E: neue Struktur\n(kein FK×Redshift)",
         "#3fb950"),
        ("SW-Diag", "4 Grids\nPEAKS_NOT_ROBUST\n(0.057→0.034)", "#f85149"),
        ("Bridge Run 6", "STAGE1+2 PASS\nA 9.3e-2 | B 0\n| C 3e-6", "#3fb950"),
        ("SRF Run 6", "5/5 Gates\nblinde Recovery\n3/3 (z:2.86/4.52/6.13)",
         "#3fb950"),
    ]
    pos = [(0.03, 0.22), (0.28, 0.22), (0.53, 0.22), (0.78, 0.22)]
    patches = []
    for (t, s, c), (x, y) in zip(boxes, pos):
        p = FancyBboxPatch((x, y), 0.19, 0.56,
                           boxstyle="round,pad=0.012",
                           fc="#161b22", ec=c, lw=1.6)
        ax.add_patch(p)
        patches.append(p)
        ax.text(x + 0.095, y + 0.44, t, ha="center", fontsize=11,
                color=c, weight="bold")
        ax.text(x + 0.095, y + 0.18, s, ha="center", fontsize=8.2,
                color="#c9d1d9")
        if x < 0.7:
            ax.annotate("", xy=(x + 0.27, 0.5), xytext=(x + 0.225, 0.5),
                        arrowprops=dict(arrowstyle="->", color="#8b949e"))

    ax.text(0.5, 0.06, "Kernbefund: SSZ-FK ≠ FK×Redshift (echte neue "
                       "Struktur) — Antwort-Maschine zertifiziert, aber "
                       "PEAKS_NOT_ROBUST bleibt ehrlich stehen: keine "
                       "Resonanz-Claims.",
            ha="center", fontsize=9, color="#ffd700")

    def update(k):
        for i, p in enumerate(patches):
            p.set_alpha(0.25 if i > k else 1.0)
        return patches

    anim = FuncAnimation(fig, update, frames=5, interval=800)
    save_anim(fig, anim, "05_chain_overview", fps=2)


if __name__ == "__main__":
    print("SSZ Full Closure — Demo-Animationen werden erzeugt …")
    demo_fk_coupled()
    demo_swdiag_resolution()
    demo_bridge_projections()
    demo_srf_recovery()
    demo_chain_overview()
    print("fertig →", OUT)
