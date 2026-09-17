"""Necessary finite-L test of the selected central Full-SVT coefficient member.

A local FAIL prevents global closure. A local PASS never certifies the remaining
regions, interfaces, radial/angular stability, or the coupled spectrum.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DEFAULT_L, SLOT_NAMES
from ..jets.jet9d8 import profile_derivative
from ..numerics import module
from ..provenance.manifest import sha256
from ..reducer.kinetic_schur import kinetic_schur
from .member import MEMBER_FILE, validate_stream_regions
from .regional_coefficients import central_selected, select_lower
from .sources import SOURCE_REGISTRY


def audit_central_kinetic(root: Path):
    root = Path(root)
    source = root / SOURCE_REGISTRY["central_exact_SVT"]["unreduced_even"]
    raw = central_selected(root).sort_values("x").reset_index(drop=True)
    red = module("ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py")
    scans, convergence = [], []
    normalized = central_selected(root).sort_values("x").reset_index(drop=True)
    normalized["region"] = "central_exact_SVT"
    production = (normalized.u >= 0.61) & (normalized.u < 0.71)
    validate_stream_regions(normalized.loc[production])
    for L in DEFAULT_L:
        result = kinetic_schur(normalized, L)
        K = result["K"]
        eig = np.linalg.eigvalsh((K + K.swapaxes(1, 2)) / 2)
        inds = np.flatnonzero(production)
        i = int(inds[np.argmin(eig[inds, 0])])
        scans.append(
            dict(
                L=L,
                rows=len(inds),
                min_eig_K=float(eig[i, 0]),
                x=float(normalized.x.iloc[i]),
                u=float(normalized.u.iloc[i]),
                negative_rows=int(np.sum(eig[inds, 0] <= 0)),
                min_abs_Dh1=float(np.min(abs(result["Dh1"][inds]))),
                min_abs_auxiliary_determinant=float(
                    np.min(abs(result["auxiliary_determinant"][inds]))
                ),
                pass_kinetic=bool(np.all(np.isfinite(K[inds])) and np.all(eig[inds] > 0)),
            )
        )
    for stride, window, degree in [
        (1, 9, 8),
        (2, 9, 8),
        (4, 9, 8),
        (1, 7, 6),
        (1, 11, 8),
        (1, 13, 8),
    ]:
        d = select_lower(raw.iloc[::stride].reset_index(drop=True), window=window, degree=degree)
        result = kinetic_schur(d, 6, window=window, degree=degree)
        replay = red.canonical_audit(d, 6, window, degree)
        i = int(np.argmin(abs(d.u.to_numpy() - 0.69)))
        K = result["K"]
        w, v = np.linalg.eigh((K[i] + K[i].T) / 2)
        # Two nonsingular constant changes of normalization preserve inertia.
        inertia = []
        for scales in ([1.0, 1.0, 1.0], [2.0, 0.5, 3.0], [0.2, 4.0, 1.0]):
            T = np.diag(scales)
            inertia.append(int(np.sum(np.linalg.eigvalsh(T.T @ K[i] @ T) < 0)))
        interior = (d.u > 0.62) & (d.u < 0.70)
        packet = np.exp(-0.5 * ((d.x.to_numpy() - d.x.iloc[i]) / 0.003) ** 2)[:, None] * v[:, 0]
        packet_prime = profile_derivative(d.x.to_numpy(), packet, 1, window, degree)
        z = np.concatenate([packet, packet_prime], axis=1)
        density = np.einsum("ni,nij,nj->n", z, result["quadratic_form"], z)
        canonical_density = np.einsum("ni,nij,nj->n", packet, K, packet)
        energy = float(np.trapezoid(density, d.x))
        canonical_energy = float(np.trapezoid(canonical_density, d.x))
        convergence.append(
            dict(
                stride=stride,
                window=window,
                degree=degree,
                packet_kinetic_integral=energy,
                packet_canonical_integral=canonical_energy,
                x=float(d.x.iloc[i]),
                u=float(d.u.iloc[i]),
                A0prime=float(d.A0prime.iloc[i]),
                eigenvalues=w.tolist(),
                K=K[i].tolist(),
                negative_direction=v[:, 0].tolist(),
                rayleigh=float(v[:, 0] @ K[i] @ v[:, 0]),
                negative_inertia=inertia,
                min_abs_Dh1=float(np.min(abs(result["Dh1"][interior]))),
                replay_scaled_max=float(
                    np.max(
                        abs(K[interior] - replay["K"][interior]) / np.maximum(1, abs(K[interior]))
                    )
                ),
                mixed_second_max=float(np.max(abs(result["mixed_second"][interior]))),
                antisymmetric_cross_max=float(np.max(abs(result["antisymmetric_cross"][interior]))),
            )
        )
    ref = (
        pd.read_csv(root / SOURCE_REGISTRY["central_exact_SVT"]["coeff_reference"])
        .sort_values("x")
        .reset_index(drop=True)
    )
    if not np.array_equal(normalized.x, ref.x):
        raise ValueError("central source/reference grid mismatch")
    error = np.abs(normalized[list(SLOT_NAMES)].to_numpy() - ref[list(SLOT_NAMES)].to_numpy())
    error /= np.maximum(1, np.abs(ref[list(SLOT_NAMES)].to_numpy()))
    validation = dict(
        central_reference_match=bool(error.max() < 1e-7),
        independent_replay_match=all(c["replay_scaled_max"] < 1e-7 for c in convergence),
        no_mixed_kinetic_order=all(
            c["mixed_second_max"] < 1e-10 and c["antisymmetric_cross_max"] < 1e-10
            for c in convergence
        ),
        nonzero_constraint_pivots=all(
            s["min_abs_Dh1"] > 1e-10 and s["min_abs_auxiliary_determinant"] > 1e-10 for s in scans
        ),
        all_required_L_kinetic_positive=all(s["pass_kinetic"] for s in scans),
    )
    status = "PASS" if all(validation.values()) else "FAIL"
    return dict(
        status=status,
        validation=validation,
        scope="necessary central finite-L kinetic gate only",
        region="central_exact_SVT",
        basis=["psi", "dphi", "V"],
        action_sha256=sha256(root / MEMBER_FILE),
        sources=[
            dict(path=str(source.relative_to(root)), sha256=sha256(source)),
            dict(
                path=SOURCE_REGISTRY["central_exact_SVT"]["lower_jets"],
                sha256=sha256(root / SOURCE_REGISTRY["central_exact_SVT"]["lower_jets"]),
            ),
            dict(
                path=SOURCE_REGISTRY["central_exact_SVT"]["coeff_reference"],
                sha256=sha256(root / SOURCE_REGISTRY["central_exact_SVT"]["coeff_reference"]),
            ),
        ],
        implementation=[
            dict(path=name, sha256=sha256(root / name))
            for name in (
                "src/ssz_p5/reducer/kinetic_schur.py",
                "src/ssz_p5/production/regional_coefficients.py",
                "src/ssz_p5/jets/jet9d8.py",
                "src/ssz_p5_profile_operator_reducer_JET9D8_2026-09-16.py",
                "src/ssz_hybrid_full_constraint_maps_JET9D8(1).py",
            )
        ],
        central_reference_max_scaled_error=float(error.max()),
        scans=scans,
        convergence=convergence,
        global_direct_status="NOT_CERTIFIED",
        coupled_spectral_status="NOT_RUN",
        absolute_closure="NOT_CERTIFIED",
        interpretation=(
            "The locked central unreduced coefficient member fails the required kinetic "
            "positivity test. This is not a no-go theorem for P5 geometry or all regional "
            "Full-SVT actions."
        ),
    )
