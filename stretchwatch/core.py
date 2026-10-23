"""StretchWatch: how far does the Hartree-Fock reference stay useful along a bond-stretching curve?

For a molecule at geometry R the qubit Hamiltonian (Jordan-Wigner, STO-3G; PennyLane built-in quantum chemistry) is diagonalised exactly inside the
correct electron-number sector; we report E_exact, the HF energy error, the HF-state overlap |<HF|GS>|^2 (VQE initial-state quality) and the mean
single-spin-orbital entanglement entropy of the ground state (a proxy for the circuit resources an ansatz needs).
"""

from __future__ import annotations

import numpy as np
import pennylane as qml
from scipy.sparse.linalg import eigsh


def hamiltonian_matrix(symbols, coords_bohr, charge=0, basis="sto-3g"):
    H, n = qml.qchem.molecular_hamiltonian(symbols, np.asarray(coords_bohr, dtype=float), charge=charge, basis=basis)
    return _sparse_from_terms(H, n), n


_P = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1.0, -1.0]).astype(complex)}


def _sparse_from_terms(H, n):
    """Build the sparse qubit Hamiltonian from its Pauli terms (avoids PennyLane's sparse_matrix, which needs a matching SciPy)."""
    import scipy.sparse as sp
    coeffs, ops = H.terms()
    total = sp.csr_matrix((2 ** n, 2 ** n), dtype=complex)
    for c, op in zip(coeffs, ops):
        word = ["I"] * n
        for w, name in zip(op.wires, [o.name for o in (op.operands if hasattr(op, "operands") else [op])]):
            word[int(w)] = name[-1] if name.startswith("Pauli") else name
        m = sp.identity(1, dtype=complex, format="csr")
        for ch in word:
            m = sp.kron(m, sp.csr_matrix(_P[ch]), format="csr")
        total = total + complex(c) * m
    return total.tocsr()


def number_operator_diag(n):
    idx = np.arange(2 ** n)
    return np.array([bin(i).count("1") for i in idx], dtype=float)   # wire 0 = most significant bit; popcount is order independent


def ground_state_in_sector(Hs, n, electrons, penalty=50.0):
    N = number_operator_diag(n)
    from scipy.sparse import diags
    Hp = Hs + penalty * diags((N - electrons) ** 2)
    w, v = eigsh(Hp, k=1, which="SA")
    psi = v[:, 0]
    e = float(np.real(np.vdot(psi, Hs @ psi)))
    return e, psi


def hf_index(n, electrons):
    bits = [1] * electrons + [0] * (n - electrons)
    return int("".join(map(str, bits)), 2)


def hf_overlap_and_energy(Hs, n, electrons, psi):
    i = hf_index(n, electrons)
    e_hf = float(np.real(Hs[i, i]))
    return float(abs(psi[i]) ** 2), e_hf


def mean_single_qubit_entropy(psi, n):
    t = psi.reshape([2] * n); s = 0.0
    for q in range(n):
        m = np.moveaxis(t, q, 0).reshape(2, -1)
        rho = m @ m.conj().T
        ev = np.clip(np.linalg.eigvalsh(rho), 1e-15, 1)
        s += float(-np.sum(ev * np.log2(ev)))
    return s / n


def point(symbols, coords_bohr, electrons, charge=0):
    Hs, n = hamiltonian_matrix(symbols, coords_bohr, charge)
    e, psi = ground_state_in_sector(Hs, n, electrons)
    ov, ehf = hf_overlap_and_energy(Hs, n, electrons, psi)
    return {"n_qubits": n, "E_exact": e, "E_HF": ehf, "corr_energy_mHa": 1000 * (ehf - e), "hf_overlap": ov, "mean_qubit_entropy": mean_single_qubit_entropy(psi, n)}
