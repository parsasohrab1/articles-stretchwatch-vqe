"""StretchWatch analysis layer: reference-quality metrics, threshold detection and a reference-state recommendation (FR-03..FR-06).

All quantities are exact for the model Hamiltonian (STO-3G, Jordan-Wigner, PennyLane interleaved spin ordering: qubits 2k, 2k+1 = spatial orbital k).
The quality index and the recommendation rule are heuristics defined here, NOT validated against hardware or larger bases.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from .core import ground_state_in_sector, hamiltonian_matrix, hf_overlap_and_energy


def _bit(n: int, q: int) -> int:
    return n - 1 - q                                    # wire 0 is the most significant bit of the basis index


def excitation_levels(n: int, electrons: int) -> np.ndarray:
    """Excitation level of every computational basis state relative to the HF determinant (electrons in virtual spin orbitals)."""
    virt = sum(1 << _bit(n, q) for q in range(electrons, n))
    return np.array([bin(i & virt).count("1") for i in range(2 ** n)])


def cisd_weight(psi: np.ndarray, n: int, electrons: int) -> float:
    """Probability mass of the ground state inside the HF + singles + doubles space (1.0 means CISD is exact for this state)."""
    lev = excitation_levels(n, electrons)
    return float(np.sum(np.abs(psi[lev <= 2]) ** 2))


def dominant_determinant(psi: np.ndarray) -> tuple:
    i = int(np.argmax(np.abs(psi)))
    return i, float(abs(psi[i]) ** 2)


def reduced_entropy(psi: np.ndarray, n: int, qubits: List[int]) -> float:
    t = psi.reshape([2] * n)
    t = np.moveaxis(t, qubits, list(range(len(qubits)))).reshape(2 ** len(qubits), -1)
    rho = t @ t.conj().T
    ev = np.clip(np.linalg.eigvalsh(rho), 1e-15, 1)
    return float(-np.sum(ev * np.log2(ev)))


def orbital_entropies(psi: np.ndarray, n: int) -> List[float]:
    """Von Neumann entropy (bits) of each spatial orbital = qubit pair (2k, 2k+1); maximum 2 bits."""
    return [reduced_entropy(psi, n, [2 * k, 2 * k + 1]) for k in range(n // 2)]


def start_quality_index(hf_overlap: float, mean_orbital_entropy: float) -> float:
    """Heuristic index in [0, 1]: HF overlap discounted by mean orbital entropy (normalised by its 2-bit maximum)."""
    return float(np.clip(hf_overlap * (1.0 - mean_orbital_entropy / 2.0), 0.0, 1.0))


def recommend_reference(hf_overlap: float, det_overlap: float, cisd: float, threshold: float = 0.9) -> str:
    if hf_overlap >= threshold:
        return "HF"
    if det_overlap > hf_overlap + 0.05 and det_overlap >= threshold:
        return "dominant determinant"
    if cisd >= 0.99:
        return "CISD-type reference"
    return "multireference / adaptive ansatz"


def analyse_point(symbols, r_bohr: float, electrons: int, charge: int = 0, threshold: float = 0.9) -> Dict:
    Hs, n = hamiltonian_matrix(symbols, [0, 0, 0, 0, 0, r_bohr], charge)
    e, psi = ground_state_in_sector(Hs, n, electrons)
    ov, ehf = hf_overlap_and_energy(Hs, n, electrons, psi)
    ent = orbital_entropies(psi, n)
    _, det_ov = dominant_determinant(psi)
    cw = cisd_weight(psi, n, electrons)
    return {"R_bohr": r_bohr, "n_qubits": n, "E_exact": e, "E_HF": ehf, "corr_energy_mHa": 1000 * (ehf - e), "hf_overlap": ov,
            "dominant_det_overlap": det_ov, "cisd_weight": cw, "mean_orbital_entropy": float(np.mean(ent)), "max_orbital_entropy": float(np.max(ent)),
            "start_quality": start_quality_index(ov, float(np.mean(ent))), "recommended_reference": recommend_reference(ov, det_ov, cw, threshold)}


def find_threshold(rs: List[float], overlaps: List[float], threshold: float = 0.9) -> Optional[float]:
    """First R (linear interpolation) at which the HF overlap falls below the threshold; None if it never does, rs[0] if it starts below."""
    if overlaps[0] < threshold:
        return float(rs[0])
    for r0, o0, r1, o1 in zip(rs, overlaps, rs[1:], overlaps[1:]):
        if o0 >= threshold > o1:
            return float(r0 + (o0 - threshold) * (r1 - r0) / (o0 - o1))
    return None
