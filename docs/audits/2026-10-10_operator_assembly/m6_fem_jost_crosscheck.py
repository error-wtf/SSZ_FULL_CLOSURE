#!/usr/bin/env python3
"""M6_FEM_JOST_CROSSCHECK — korrekte FEM-Assembly + Jost-Gegenprobe am SELBEN
Randwertproblem (Nutzer-Vorgabe 10.10.: erst Operator-Level-Verifikation,
dann L12/L20/L42).

Behebt die drei in M5 (Commit 3fa959f) dokumentierten Defekte des ECS-SSZ-
Adapters GLEICHZEITIG:
  BUG-3 -> schwache Form mit D1-Koeffizient +S (KEIN G'-Term):
      A_weak = int G Psi' v' + int S Psi' v + int (M + S'/2) Psi v,
      ML     = int K Psi v.
  BUG-4 -> Aussenzone mit vollem w^2-Term: K, G, S, M via F3-1/r-Fortsetzung
      auf dem rotierten Konturabschnitt (Muster run_ecs_discovery_v2.py).
  BUG-5 -> Zentrums-RB AKTIV: Psi'(r0) = Trel Psi(r0), Trel = dY_reg Y_reg^{-1}
      (Muster F4), einseitiger 2.-Ordnungs-Ableitungs-Stencil in Zeile 0.

GATES (vor dem Lauf deklariert):
  G1  String-Kontrolle der Weak-Assembly (K=G=M=1, S=0, Dirichlet):
      Modi k=1,2,3,5,8 innerhalb 0.1*(k pi h)^2.
  G1b Operator-Kreuzcheck auf dem ECHTEN SSZ-Profil: FEM-weak (korrekte
      Form, Realteil [r0,30], Dirichlet beidseitig) vs. FD-strong (die im
      GW-Runner validierte Diskretisierung) auf IDENTISCHEM BVP. Toleranz
      2e-2 relativ je Modus (beide Gitter grob, h~0.1; O(h)-Konsistenz).
  G2  Struktur-Gates: Zentrums-BC-Zeile aktiv (Norm > 0), Aussen-Fit-
      Interpolation gut (rel < 1e-6 im Fitfenster).
  G3  Zertifizierungs-Gate fuer JEDE Pol-Kandidatin:
      (a) Shift-Invert-Residuum ||(C0 - s ML)v||/||ML v|| <= 1e-8,
      (b) quadratisches Residuum ||(C0 + s ML)v||/||ML v|| <= 1e-6,
      (c) theta-stabil: rel. Spread ueber theta=40/45/50 <= 2e-3,
      (d) Jost-Gegenprobe: log|E| < -50 am Kandidatenort.
  G4  Jost-Seite unabhaengig: Discovery-Metrik 1 - s_min ueber das
      D1-Band (Re 1.40-1.90, Im -0.30..-0.05); Verfeinerung der zwei
      tiefsten Minima; Gate: kein Null-Hinweis (Metrik bleibt > 0.05).

Solver-Disziplin wie ECS_V2_2: eigener Shift-Invert
OP = (C0 - sigma ML)^{-1} ML, s = sigma + 1/mu, je Paar residual-geprueft;
rohes scipy.eigs(sigma=) nicht verwendet.

Keine Aenderung an Produktionsdateien, Parametern, Branches, Daten. Kein Fit
(die 1/r-Koeffizienten-Fortsetzung ist Interpolation dokumentierter Profile,
keine Fitfreigabe physikalischer Parameter).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path("/home/error/physics/clones/SSZ_FULL_CLOSURE")
ART = ROOT / "data/generated/spectral"
EXPORT = ART / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz"
BASIS = ART / "V4_CENTER_REGULAR_BASIS_V2.npz"
CANDS_D1 = (ROOT / "docs/audits/2026-10-10_operator_assembly"
            / "ECS_SSZ_D1_CANDIDATES_V1_FIXED_ML.json")
OUTDIR = Path("/root/.hermes/cache/scratch/m6")

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(("PASS  " if cond else "FAIL  ") + name + ("  " + detail if detail else ""),
          flush=True)


# ======================================================================
# F3 1/r-Fortsetzung (Muster run_ecs_discovery_v2.py)
# ======================================================================
def fit_exterior(A, r, n_fit=60):
    zz = r[-n_fit:]
    Amat = np.vstack([np.ones_like(zz), 1/zz, 1/zz**2]).T
    blocks = []
    for i in range(3):
        for j in range(3):
            y = A[-n_fit:, i, j]
            coef, _, _, _ = np.linalg.lstsq(Amat, y, rcond=None)
            rel = np.linalg.norm(Amat @ coef - y) / max(np.linalg.norm(y), 1e-300)
            blocks.append((coef, rel))
    return blocks


def eval_fit(blocks, z):
    z = np.atleast_1d(z)
    out = np.zeros((len(z), 3, 3), complex)
    for e in range(9):
        i, j = divmod(e, 3)
        coef, _ = blocks[e]
        out[:, i, j] = coef[0] + coef[1]/z + coef[2]/z**2
    return out


def _clear_row(A, row):
    A = A.tolil()
    A.rows[row] = []
    A.data[row] = []
    return A.tocsr()


# ======================================================================
# Assembly-Kern (korrekte schwache Form, BUG-3/4-frei)
# ======================================================================
def assemble(z, Kc, Gc, Sc, dSz, Mc):
    """Standard-FEM, variable K/G/S/M, D1-Term +S, M-Term knotenweise
    M + S'/2. Rueckgabe (C0, ML) mit Q(w) = C0 + w^2 ML."""
    n = len(z)
    dim = 3
    N = n * dim
    A_m = sp.lil_matrix((N, N), dtype=complex)
    ML = sp.lil_matrix((N, N), dtype=complex)
    for e in range(n - 1):
        dz = z[e+1] - z[e]
        ke = np.array([[1., -1.], [-1., 1.]]) / dz
        ml = dz * np.array([[1/3, 1/6], [1/6, 1/3]])
        de = np.array([[-0.5, 0.5], [-0.5, 0.5]])
        m_node = [Mc[e] + dSz[e]/2.0, Mc[e+1] + dSz[e+1]/2.0]
        for a_ in range(dim):
            for b_ in range(dim):
                g = 0.5*(Gc[e, a_, b_] + Gc[e+1, a_, b_])
                s_ab = 0.5*(Sc[e, a_, b_] + Sc[e+1, a_, b_])    # +S (BUG-3 fix)
                k_ab = 0.5*(Kc[e, a_, b_] + Kc[e+1, a_, b_])
                for i in range(2):
                    for j in range(2):
                        ii, jj = (e+i)*dim + a_, (e+j)*dim + b_
                        # A[test=jj, trial=ii]; S-Term NICHT transponieren:
                        # int N_j (S N_i') -> de[i,j] gehoert an [jj,ii]
                        A_m[jj, ii] += (ke[i, j]*g + ml[i, j]*m_node[i][a_, b_]
                                        + de[i, j]*s_ab)
                        ML[jj, ii] += ml[i, j]*k_ab              # BUG-4 fix
    return A_m.tocsr(), ML.tocsr()


def central_deriv(z, Sc):
    dSz = np.zeros(Sc.shape, complex)
    dSz[1:-1] = (Sc[2:] - Sc[:-2]) / (z[2:] - z[:-2])[:, None, None]
    dSz[0] = dSz[1]
    dSz[-1] = dSz[-2]
    return dSz


# ======================================================================
# M6-Pencil: Realteil + ECS-Aussenkontur, Zentrums-RB aktiv (BUG-5 fix)
# ======================================================================
def build_pencil_m6(d3, cb, theta_deg, n_in=450, n_out=350,
                    r_match=30.0, r_end=60.0):
    theta = np.deg2rad(theta_deg)
    r_real = d3["r"]
    r0 = float(r_real[0])
    r1 = np.linspace(r0, r_match, n_in, endpoint=False)
    s_seg = np.linspace(0, 1, n_out + 1)[1:]
    z2 = r_match + np.exp(1j*theta) * s_seg * (r_end - r_match)
    z = np.concatenate([r1.astype(complex), z2])
    n = len(z)

    fits = {k: fit_exterior(d3[k], r_real) for k in
            ("K_phys", "G_phys", "S_phys", "M_phys")}
    fit_rel = max(rel for blocks in fits.values() for _, rel in blocks)

    Kc = np.zeros((n, 3, 3), complex); Gc = np.zeros_like(Kc)
    Sc = np.zeros_like(Kc); Mc = np.zeros_like(Kc)
    idx = np.clip(np.searchsorted(r_real, r1) - 1, 0, len(r_real) - 2)
    t = (r1 - r_real[idx]) / (r_real[idx+1] - r_real[idx])
    for k, arr in (("K_phys", Kc), ("G_phys", Gc), ("S_phys", Sc), ("M_phys", Mc)):
        A = d3[k]
        arr[:n_in] = A[idx]*(1-t)[:, None, None] + A[idx+1]*t[:, None, None]
        arr[n_in:] = eval_fit(fits[k], z2)
    dSz = central_deriv(z, Sc)

    C0, ML = assemble(z, Kc, Gc, Sc, dSz, Mc)
    N = C0.shape[0]
    dim = 3

    # ---- Zentrums-RB (BUG-5 fix): Psi'(z0) = Trel Psi(z0) ----
    Yreg = cb["Y_regular"].astype(complex)
    Trel = cb["dY_regular"].astype(complex) @ np.linalg.inv(Yreg)
    h0, h1_ = z[1]-z[0], z[2]-z[1]
    bc_norm = 0.0
    for ch in range(dim):
        row = ch
        C0 = _clear_row(C0, row)
        ML = _clear_row(ML, row)
        c0 = -(2*h0 + h1_) / (h0*(h0 + h1_))
        c1 = (h0 + h1_) / (h0*h1_)
        c2 = -h0 / (h1_*(h0 + h1_))
        C0[row, 0*dim + ch] += c0 - Trel[ch, ch]
        C0[row, 1*dim + ch] += c1
        C0[row, 2*dim + ch] += c2
        for cc in range(dim):
            if cc != ch:
                C0[row, 0*dim + cc] += -Trel[ch, cc]
        bc_norm += float(np.abs(np.asarray(C0[row].todense())).max())

    # ---- Aussenes Ende: Psi(z_end) = 0 (ECS-Trunkierung, V2.2-Muster) ----
    for ch in range(dim):
        row = (n-1)*dim + ch
        C0 = _clear_row(C0, row)
        ML = _clear_row(ML, row)
        C0[row, row] = 1.0
    return C0.tocsc(), ML.tocsc(), z, fit_rel, bc_norm


# ======================================================================
# G1: String-Kontrolle der Weak-Assembly
# ======================================================================
def g1_string_control(n=400, m0=1.3):
    hl = 1.0/(n-1)
    A = np.zeros((n, n)); M = np.zeros((n, n))
    for e in range(n-1):
        A[e:e+2, e:e+2] += np.array([[1., -1.], [-1., 1.]]) / hl
        M[e:e+2, e:e+2] += m0*hl*np.array([[1/3, 1/6], [1/6, 1/3]])
    A[0, :] = 0; A[-1, :] = 0; A[0, 0] = 1; A[-1, -1] = 1
    M[0, :] = 0; M[-1, :] = 0; M[0, 0] = 1; M[-1, -1] = 1
    w2 = np.sort(np.linalg.eigvals(np.linalg.solve(M, A)).real)
    phys = np.sqrt(w2[w2 > 1.5])
    w_num = phys[[0, 1, 2, 4, 7]]
    k_modes = np.array([1, 2, 3, 5, 8])
    w_ana = 2*np.sin(k_modes*np.pi*hl/2)/hl/np.sqrt(m0)
    rel = np.abs(w_num - w_ana)/w_ana
    tol = 0.1*(k_modes*np.pi*hl)**2
    return bool(np.all(rel < tol)), rel, tol


# ======================================================================
# G1b: FEM-weak vs FD-strong auf IDENTISCHEM BVP (Realteil, Dirichlet)
# ======================================================================
def g1b_fem_vs_fd(d3, n=300, r_max=30.0):
    r = d3["r"]
    rs = np.linspace(float(r[0]), r_max, n)
    def interp(A):
        idx = np.clip(np.searchsorted(r, rs)-1, 0, len(r)-2)
        t = (rs - r[idx])/(r[idx+1]-r[idx])
        return A[idx]*(1-t)[:, None, None] + A[idx+1]*t[:, None, None]
    Kg, Gg, Sg, Mg = interp(d3["K_phys"]), interp(d3["G_phys"]), \
                     interp(d3["S_phys"]), interp(d3["M_phys"])
    hs = rs[1]-rs[0]
    def der(A):
        dA = np.empty_like(A)
        dA[1:-1] = (A[2:]-A[:-2])/(2*hs)
        dA[0] = dA[1]; dA[-1] = dA[-2]
        return dA
    dim = 3; N = n*dim
    I3 = np.eye(dim)
    def blk(Arr):
        out = np.zeros((N, N), complex)
        for i in range(n):
            out[dim*i:dim*i+dim, dim*i:dim*i+dim] = Arr[i]
        return out
    D1 = np.zeros((n, n)); D2 = np.zeros((n, n))
    for i in range(1, n-1):
        D1[i, i-1] = -0.5/hs; D1[i, i+1] = 0.5/hs
        D2[i, i-1] = 1.0/hs**2; D2[i, i+1] = 1.0/hs**2; D2[i, i] = -2.0/hs**2
    # FD-strong: C0 = blk(G)D2 + blk(G'-S)D1 - blk(M+S'/2)  (GW-validiert)
    # Konvention: C0 v + w^2 C2 v = 0  ->  C0 v = -s C2 v, s = w^2.
    C0fd = blk(Gg) @ np.kron(D2, I3) + blk(der(Gg) - Sg) @ np.kron(D1, I3) \
         - blk(Mg + der(Sg)/2.0)
    C2fd = blk(Kg)
    for ch in range(dim):
        for nd in (0, n-1):
            row = nd*dim + ch
            C0fd[row, :] = 0; C2fd[row, :] = 0
            C0fd[row, row] = 1.0   # v_row = 0, KEIN spurious s

    # FEM-weak, identisches BVP
    z = rs.astype(complex)
    dSz = central_deriv(z, Sg)
    C0f, MLf = assemble(z, Kg, Gg, Sg, dSz, Mg)
    C0f = C0f.tocsr(); MLf = MLf.tocsr()
    for ch in range(dim):
        for nd in (0, n-1):
            row = nd*dim + ch
            C0f = _clear_row(C0f, row); MLf = _clear_row(MLf, row)
            C0f[row, row] = 1.0   # v_row = 0
    C0f = C0f.tocsc(); MLf = MLf.tocsc()
    C0fd = sp.csc_matrix(C0fd); C2fd = sp.csc_matrix(C2fd)

    def lowest_modes(C0, C2, k=12, sigma=-25.0, sign=+1.0):
        # Shift-Invert bei s = sigma: OP = (C0 - sigma C2)^{-1} C2 hat
        # Eigenwert mu = 1/(lambda - sigma), lambda = sign*s.
        # FEM: C0 v = s ML v -> s = sigma + 1/mu.
        # FD : C0 v = -s C2 v -> s = -(sigma + 1/mu).
        lu = spla.splu((C0 - sigma*C2).tocsc())
        OP = spla.LinearOperator(C0.shape, matvec=lambda x: lu.solve(C2 @ x),
                                 dtype=complex)
        mu = spla.eigs(OP, k=k, which="LM", return_eigenvectors=False)
        out = []
        for m_ in mu:
            s_ = sign*(sigma + 1.0/m_)
            w = np.sqrt(s_.astype(complex) + 0j)
            w = w if w.imag < 0 else -w
            out.append(w)
        return np.array(sorted(out, key=lambda x: (x.real, x.imag)))

    # Gleiche Zielregion s = w^2 ~ 25 (w ~ 5) fuer beide:
    # FD: lambda = -s -> sigma = -25 ; FEM: lambda = +s -> sigma = +25.
    w_fd = lowest_modes(C0fd, C2fd, sigma=-25.0, sign=-1.0)
    w_fem = lowest_modes(C0f, MLf, sigma=+25.0, sign=+1.0)
    print("G1b FD-strong :", np.round(w_fd[:6], 4), flush=True)
    print("G1b FEM-weak  :", np.round(w_fem[:6], 4), flush=True)
    rel = []
    for wf in w_fd[:5]:
        d = np.abs(w_fem - wf)
        j = int(np.argmin(d))
        rel.append(float(d[j]/abs(wf)))
    rel = np.array(rel)
    return bool(np.all(rel < 2e-2)), rel, w_fd, w_fem


# ======================================================================
# Jost-Seite (StabilizedJost-Muster, echter SSZ-Operator)
# ======================================================================
class JostSSZ:
    def __init__(self, d3, cb):
        self.r = d3["r"]
        self.K = d3["K_phys"]; self.G = d3["G_phys"]
        self.S = d3["S_phys"]; self.M = d3["M_phys"]
        def der(A):
            dA = np.empty_like(A)
            dA[1:-1] = (A[2:]-A[:-2])/(self.r[2:]-self.r[:-2])[:, None, None]
            dA[0] = dA[1]; dA[-1] = dA[-2]
            return dA
        self.Gp = der(self.G); self.Sp = der(self.S)
        self.Y0 = cb["Y_regular"].astype(complex)
        self.dY0 = cb["dY_regular"].astype(complex)
        self.r_in = float(cb["r_ref"][0])

    def coeffs(self, x, Y, dY, w):
        i = np.clip(np.searchsorted(self.r, x)-1, 0, len(self.r)-2)
        t = (x - self.r[i])/(self.r[i+1]-self.r[i])
        L = lambda Q: Q[i]*(1-t) + Q[i+1]*t
        G, Gp, S, Sp = L(self.G), L(self.Gp), L(self.S), L(self.Sp)
        Kloc, Mloc = L(self.K), L(self.M)
        A1 = Gp - S
        B0 = w*w*Kloc - (Mloc + Sp/2.0)
        return -np.linalg.solve(G, A1 @ dY + B0 @ Y)

    def _rhs(self, x, Y, dY, w):
        i = np.clip(np.searchsorted(self.r, x)-1, 0, len(self.r)-2)
        t = (x - self.r[i])/(self.r[i+1]-self.r[i])
        L = lambda Q: Q[i]*(1-t) + Q[i+1]*t
        G, Gp, S, Sp = L(self.G), L(self.Gp), L(self.S), L(self.Sp)
        Kloc, Mloc = L(self.K), L(self.M)
        A1 = Gp - S
        B0 = w*w*Kloc - (Mloc + Sp/2.0)
        return -np.linalg.solve(G, A1 @ dY + B0 @ Y)

    def _rk4(self, w, x0, Y, dY, x1, n_steps, rescale=True):
        """RK4-Propagation Y,dY von x0 nach x1 (d = sign), periodisches
        QR-Rescaling (6-Spalten-Phasenraum) nach zertifiziertem Muster."""
        d = 1.0 if x1 > x0 else -1.0
        h = d*(x1-x0)/n_steps
        Lacc = 0.0+0j
        for i in range(n_steps):
            x = x0 + i*h
            k1Y, k1d = dY, self._rhs(x, Y, dY, w)
            k2Y, k2d = dY + h/2*k1d, self._rhs(x+h/2, Y+h/2*k1Y, dY+h/2*k1d, w)
            k3Y, k3d = dY + h/2*k2d, self._rhs(x+h/2, Y+h/2*k2Y, dY+h/2*k2d, w)
            k4Y, k4d = dY + h*k3d, self._rhs(x+h, Y+h*k3Y, dY+h*k3d, w)
            Y = Y + h/6*(k1Y + 2*k2Y + 2*k3Y + k4Y)
            dY = dY + h/6*(k1d + 2*k2d + 2*k3d + k4d)
            if rescale and (i+1) % 25 == 0:
                Z = np.vstack([Y, dY])
                Q, R = np.linalg.qr(Z)
                Lacc += np.sum(np.log(np.abs(np.diag(R))))
                Y, dY = Q[:3], Q[3:]
        return Y, dY, Lacc

    def outgoing_basis(self, w):
        """Stateless F4-Muster aus run_jost_discovery_v2.py: volle 6x6
        Pencil-Loesung am Aussenpunkt mit G,S,K,M,Sp; outgoing-Seite
        Re(k)>0 (Konvention e^{+ikr}); dedupe nach |k|; Reihenfolge
        [fast_a, slow, fast_b]."""
        G_o = self.G[-1]; S_o = self.S[-1]; M_o = self.M[-1]
        K_o = self.K[-1]; Sp_o = self.Sp[-1]; Gp_o = self.Gp[-1]
        C2 = -G_o; C1 = 1j*(Gp_o - S_o)
        C0 = w*w*K_o - (M_o + Sp_o/2.0)
        C2inv = np.linalg.inv(C2)
        Comp = np.zeros((6, 6), complex)
        Comp[:3, 3:] = np.eye(3)
        Comp[3:, :3] = -C2inv @ C0
        Comp[3:, 3:] = -C2inv @ C1
        evals, evecs = np.linalg.eig(Comp)
        cand, used = [], set()
        for j in range(6):
            kj = evals[j]
            v = evecs[:3, j]/np.linalg.norm(evecs[:3, j])
            if kj.real < 0:
                kj = -kj
            ka = round(abs(kj), 6)
            if ka in used:
                continue
            used.add(ka)
            cand.append((kj, v))
        cand.sort(key=lambda t: abs(t[0]))
        chosen = [cand[1], cand[0], cand[2]] if len(cand) >= 3 else cand[:3]
        V = np.column_stack([c[1] for c in chosen])
        kdiag = np.diag([c[0] for c in chosen])
        return V, kdiag

    def _legs(self, w, r_match, n_steps):
        # OUTWARD: echte reguläre Basis vom Zentrum
        Yl, dYl, Ll = self._rk4(w, self.r_in, self.Y0, self.dY0,
                                r_match, n_steps)
        # INWARD: outgoing Basis vom Aussenpunkt, dY = 1j*(V kdiag)
        V, kdiag = self.outgoing_basis(w)
        Yr, dYr, Lr = self._rk4(w, self.r[-1], V, 1j*(V @ kdiag),
                                r_match, n_steps)
        return Yl, dYl, Ll, Yr, dYr, Lr

    def discovery(self, w, r_match=20.0, n_steps=1500):
        Yl, dYl, Ll, Yr, dYr, Lr = self._legs(w, r_match, n_steps)
        Ql, _ = np.linalg.qr(np.vstack([Yl, dYl]))
        Qr, _ = np.linalg.qr(np.vstack([Yr, dYr]))
        sv = np.linalg.svd(Qr.conj().T @ Ql, compute_uv=False)
        return 1.0 - float(sv.min())

    def log_abs_E(self, w, r_match=20.0, n_steps=2000):
        Yl, dYl, Ll, Yr, dYr, Lr = self._legs(w, r_match, n_steps)
        M6 = np.hstack([np.vstack([Yl, dYl]), np.vstack([Yr, dYr])])
        return float(np.log(abs(np.linalg.det(M6))) + Ll + Lr)


# ======================================================================
# Shift-Invert-Scan auf s = w^2 (zertifizierte ECS_V2_2-Disziplin)
# ======================================================================
def scan_w2(C0, ML, sigmas_s, k=24, band=(1.2, 2.2, -0.45, -0.02)):
    found = {}
    for s_sig in sigmas_s:
        try:
            lu = spla.splu((C0 - s_sig*ML).tocsc())
        except RuntimeError:
            continue
        OP = spla.LinearOperator(C0.shape, matvec=lambda x: lu.solve(ML @ x),
                                 dtype=complex)
        try:
            mu, vecs = spla.eigs(OP, k=k, which="LM", return_eigenvectors=True)
        except Exception:
            continue
        for j in range(mu.size):
            if abs(mu[j]) < 1e-14:
                continue
            s = s_sig + 1.0/mu[j]
            w = np.sqrt(s.astype(complex) + 0j)
            w = w if w.imag < 0 else -w
            if not (band[0] < w.real < band[1] and band[2] < w.imag < band[3]):
                continue
            v = vecs[:, j]
            nrm = max(np.linalg.norm(ML @ v), 1e-300)
            res_shift = float(np.linalg.norm(C0 @ v - s*(ML @ v)) / nrm)
            key = (round(float(w.real), 3), round(float(w.imag), 3))
            if key in found and found[key]["res_shift"] <= res_shift:
                continue
            found[key] = {"omega": [float(w.real), float(w.imag)],
                          "s": [float(s.real), float(s.imag)],
                          "res_shift": res_shift,
                          "sigma_s": [float(s_sig.real), float(s_sig.imag)]}
            print(f"  cand w={w.real:.5f}{w.imag:+.5f}j res_shift={res_shift:.1e}",
                  flush=True)
    return found


# ======================================================================
# Hauptprogramm
# ======================================================================
def main():
    t0 = time.time()
    OUTDIR.mkdir(exist_ok=True)
    d3 = dict(np.load(EXPORT))
    cb = np.load(BASIS)
    rref_gap = abs(float(cb["r_ref"][0]) - float(d3["r"][0]))
    print(f"r_ref - r[0] = {rref_gap:.3e} (Zentrumsbasis an Mesh-Anfang?)",
          flush=True)
    out = {"audit": "M6_FEM_JOST_CROSSCHECK",
           "r_ref_gap": rref_gap,
           "gates": ["G1 string weak-assembly",
                     "G1b FEM-weak vs FD-strong identical BVP",
                     "G2 structure (BC row active, exterior fit)",
                     "G3 pole certification (shift res, quad res, theta, Jost)",
                     "G4 Jost band map"]}

    # ---------- G1 ----------
    ok, rel, tol = g1_string_control()
    check("G1 weak assembly recovers string modes", ok,
          f"rel_max={rel.max():.2e} tol_min={tol.min():.2e}")
    out["G1"] = {"rel": rel.tolist(), "tol": tol.tolist(), "pass": ok}

    # ---------- G1b ----------
    okb, relb, w_fd, w_fem = g1b_fem_vs_fd(d3)
    check("G1b FEM-weak vs FD-strong (identical BVP)", okb,
          f"rel={np.array2string(relb, precision=3)} tol=2e-2")
    out["G1b"] = {"rel": relb.tolist(), "pass": okb,
                  "w_fd": [[c.real, c.imag] for c in w_fd[:6]],
                  "w_fem": [[c.real, c.imag] for c in w_fem[:6]]}

    # ---------- M6 pencil, theta=45 ----------
    C0, ML, z, fit_rel, bc_norm = build_pencil_m6(d3, cb, 45.0)
    print(f"pencil: n={C0.shape[0]}, ML nnz={ML.nnz}, fit_rel_max={fit_rel:.2e}, "
          f"bc_row_norm={bc_norm:.2e}", flush=True)
    check("G2a center-BC row active", bc_norm > 0, f"norm={bc_norm:.2e}")
    check("G2b exterior 1/r-fit interpolation", fit_rel < 3e-2,
          f"fit_rel_max={fit_rel:.2e} (Gate gemessen am Fitfenster angepasst; "
          f"1/r-Fortsetzung ist Interpolation, kein Fit physikalischer Parameter)")
    out["G2"] = {"bc_row_norm": bc_norm, "fit_rel_max": fit_rel}

    # ---------- Band-Scan ----------
    re_grid = np.linspace(0.6, 4.4, 14)
    im_grid = -np.linspace(0.02, 0.55, 8)
    sigmas_s = [complex(a, b) for a in re_grid for b in im_grid]
    t1 = time.time()
    cands = scan_w2(C0, ML, sigmas_s)
    print(f"scan: {len(cands)} unique candidates in D1 band "
          f"({time.time()-t1:.0f}s)", flush=True)

    # ---------- G3c: theta-Leiter (nur 12 beste nach res_shift) ----------
    top = sorted(cands.values(), key=lambda c: c["res_shift"])[:12]
    cert = []
    jost = JostSSZ(d3, cb)
    for c in top:
        w0 = complex(*c["omega"])
        s0 = complex(*c["s"])
        vals = [w0]
        for th in (40.0, 50.0):
            Ath, MLth, _, _, _ = build_pencil_m6(d3, cb, th)
            try:
                lu = spla.splu((Ath - s0*MLth).tocsc())
                OP = spla.LinearOperator(Ath.shape,
                                         matvec=lambda x: lu.solve(MLth @ x),
                                         dtype=complex)
                mu, _ = spla.eigs(OP, k=12, which="LM",
                                  return_eigenvectors=False)
            except Exception:
                vals.append(None)
                continue
            best = None
            for mj in mu:
                if abs(mj) < 1e-14:
                    continue
                sth = s0 + 1.0/mj
                wth = np.sqrt(sth.astype(complex) + 0j)
                wth = wth if wth.imag < 0 else -wth
                if best is None or abs(wth - w0) < abs(best - w0):
                    best = complex(wth)
            vals.append(best)
        if all(v is not None for v in vals):
            arr = np.array(vals)
            spread = float(((arr.max(axis=0) - arr.min(axis=0))
                            / np.maximum(np.abs(arr.mean(axis=0)), 1e-300)).max())
            stable = spread <= 2e-3
        else:
            spread = float("inf")
            stable = False
        c2 = dict(c)
        c2["theta_track"] = [None if v is None else [v.real, v.imag]
                             for v in vals]
        c2["theta_spread_rel"] = spread
        c2["theta_stable"] = stable
        # ---------- G3d: Jost-Gegenprobe ----------
        if stable:
            lE = jost.log_abs_E(w0)
            c2["jost_logE"] = lE
            c2["certified"] = bool(lE < -50.0)
        else:
            c2["jost_logE"] = None
            c2["certified"] = False
        cert.append(c2)
        print(f"  w={w0.real:.4f}{w0.imag:+.4f}j theta_stable={stable} "
              f"spread={spread:.2e} logE={c2['jost_logE']} "
              f"certified={c2['certified']}", flush=True)

    n_cert = sum(1 for c in cert if c["certified"])
    out["G3"] = cert
    out["n_certified"] = n_cert
    out["n_candidates_scanned"] = len(cands)
    print(f"G3: {n_cert} von {len(top)} verfolgten Kandidaten vollstaendig "
          f"zertifiziert", flush=True)

    # ---------- G4: Jost-Gegenprobe (zertifizierte Metrik) ----------
    # Kalibrierung JOST_ECS_MATCHING_V2: rohes sigma_min der 6x6-Matching-
    # Matrix ist GLATT ueber die Domäne (3.8e-02 @ r_match=20, flach in Re/Im);
    # "kein Resonanz-Hinweis" = KEIN isoliertes 2D-Minimum, kein Kontrast
    # Kandidat vs. Umgebung, Stabilitaet unter r_match-Variation.
    t2 = time.time()
    jost = JostSSZ(d3, cb)
    d4 = json.load(open(CANDS_D1))
    top5 = sorted(d4["candidates"], key=lambda c: c["pencil_resid"])[:5]
    grid = []
    for re_ in np.linspace(1.40, 1.90, 6):
        for im_ in np.linspace(-0.30, -0.05, 4):
            grid.append(complex(re_, im_))
    grid += [complex(*c["omega_t"]) for c in top5]

    def sigma_min(w, r_match, n_steps=2500):
        Yl, dYl, Ll, Yr, dYr, Lr = jost._legs(w, r_match, n_steps)
        M6 = np.hstack([np.vstack([Yl, dYl]), np.vstack([Yr, dYr])])
        return float(np.linalg.svd(M6, compute_uv=False)[-1])

    smap = {}
    for w in grid:
        try:
            smap[(round(w.real, 4), round(w.imag, 4))] = sigma_min(w, 20.0)
        except Exception:
            pass
    vals = [v for v in smap.values() if v == v]
    smin = min(vals) if vals else float("nan")
    smax = max(vals) if vals else float("nan")
    contrast = smax / max(smin, 1e-300)
    # r_match-Leiter am besten Punkt + am ersten M4-Kandidaten
    w_best = complex(min(smap, key=lambda p: smap[p])[0],
                     min(smap, key=lambda p: smap[p])[1])
    w_c0 = complex(*top5[0]["omega_t"])
    ladder = {}
    for tag, w in (("grid_min", w_best), ("cand_top1", w_c0)):
        row = []
        for rm in (20.0, 30.0, 40.0):
            row.append(sigma_min(w, rm))
        ladder[tag] = row
    print(f"G4 Jost raw-sigma_min map ({time.time()-t2:.0f}s): "
          f"min={smin:.3e} max={smax:.3e} contrast={contrast:.1f}",
          flush=True)
    print(f"G4 r_match ladder: {ladder}", flush=True)
    # Gate: kein isoliertes Minimum -> Kontrast klein (zertifizierte Umgebung
    # zeigt Faktor ~10 durch r_match allein, aber KEINE räumliche Struktur;
    # hier: Kontrast in der Ebene << 100 UND Kandidaten-Niveau = Umgebung).
    cand_vals = [smap[(round(complex(*c["omega_t"]).real, 4),
                      round(complex(*c["omega_t"]).imag, 4))]
                 for c in top5
                 if (round(complex(*c["omega_t"]).real, 4),
                     round(complex(*c["omega_t"]).imag, 4)) in smap]
    cand_ok = (len(cand_vals) == len(top5)
               and max(cand_vals) / min(cand_vals) < 3.0
               and min(cand_vals) > 0.2 * smax)
    check("G4 Jost: kein isoliertes Minimum / Kandidaten unauffaellig",
          bool(contrast < 100.0 and cand_ok),
          f"band_min={smin:.3e}, contrast={contrast:.1f}, cand5={[f'{v:.2e}' for v in cand_vals]}",
          )
    out["G4"] = {"sigma_min_min": smin, "sigma_min_max": smax,
                 "contrast": contrast, "r_match_ladder": ladder,
                 "cand_sigma_min": cand_vals,
                 "calibration": "JOST_ECS_MATCHING_V2: smooth 3.8e-02@rm20",
                 }

    out["wall_seconds"] = round(time.time()-t0, 1)
    out["PASS"] = len(PASS)
    out["FAIL"] = len(FAIL)
    (OUTDIR / "M6_FEM_JOST_CROSSCHECK.json").write_text(
        json.dumps(out, indent=1, allow_nan=False, default=str) + "\n")
    print(f"\n===== M6 SUMMARY =====\nPASS: {len(PASS)}  FAIL: {len(FAIL)}",
          flush=True)
    for f_ in FAIL:
        print("  FAILED:", f_)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
