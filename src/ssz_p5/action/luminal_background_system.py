"""N1: luminal background ODE system (S -> P1 -> K_scalar closure).

Luminal sector: G4 = G4(phi), G3-family = 0, G4X-family = 0, G5-family = 0,
G2 = G2(X, F) with G2F = 1 (Maxwell normalization).

Builds the residual functions E00, E11, E22 from the verified Kase-
Tsujikawa expressions with the full luminal jet substitution, ready for
the N2 relaxation/shooting solve of (f, h, phi)(r).

Convention note: the repo uses u = 1/r and the project asymptotics are
f -> 1/4, h -> 1, phi -> 1 (NOT asymptotically flat).
"""
from __future__ import annotations

import sympy as sp

from .kt_mh_background import JETS, kt_e00, kt_e11, kt_e22

_LUMINAL_ZERO = []
for _s in JETS:
    _n = _s.name
    if _n.startswith("G5") or _n.startswith("G3") or _n in (
            "G4X", "G4XX", "G4phiX"):
        _LUMINAL_ZERO.append((_s, 0))


def luminal_substitutions() -> list[tuple]:
    """The complete luminal jet substitution (G5/G3/G4X families zero)."""
    return list(_LUMINAL_ZERO)


def luminal_equations() -> dict[str, sp.Expr]:
    """E00, E11, E22 with the full luminal substitution applied."""
    zero = luminal_substitutions()
    return {
        "E00": kt_e00().subs(zero),
        "E11": kt_e11().subs(zero),
        "E22": kt_e22().subs(zero),
    }


def system_free_symbols() -> dict[str, list[str]]:
    eqs = luminal_equations()
    return {name: sorted(s.name for s in e.free_symbols) for name, e in eqs.items()}


def scalar_field_equation_hint() -> str:
    """The scalar equation is NOT an independent Einstein component in the
    luminal branch (E00/E11/E22 + Bianchi close the system); N2 solves
    (f, h, phi) with phi via its own equation obtained as the linear
    combination that eliminates app (documented in N1_N4 spec)."""
    return ("E_phi = combination of E00/E11/E22 eliminating app; "
            "see docs/N1_N4_LUMINAL_RESOLVE_SPEC.md section N2")
