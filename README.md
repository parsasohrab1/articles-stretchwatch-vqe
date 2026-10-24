# StretchWatch – how good is the Hartree–Fock starting point for VQE along a bond-stretching curve?

For H₂, HeH⁺ and LiH (STO-3G) it computes, at every bond length, the exact ground state (full CI, PySCF), the overlap |⟨HF|FCI⟩|², the entropy of the natural-orbital occupations,
a **quality index** Q = |⟨HF|FCI⟩|² · exp(−S), finds the bond length beyond which HF stops being a good start, and recommends **HF or an alternative reference (CISD)**.

**Maturity: TRL 4** (laboratory validation with exact classical calculations; no quantum hardware or circuit execution, as in [docs/SRS.md](docs/SRS.md)). See [docs/TRL.md](docs/TRL.md).
Moving to TRL 5 needs real VQE runs (for example PennyLane) showing that the index predicts iteration counts or convergence failures, and larger molecules.

## Pipeline
| SRS | What it does | Where |
|-----|--------------|-------|
| FR-01 | Molecule geometry and STO-3G integrals (PySCF) | `chem.build` |
| FR-02 | Exact ground state by full CI | `chem.analyse` |
| FR-03 | Overlap of the HF determinant with the ground state | `chem.analyse` |
| FR-04 | Natural-orbital entropy from the FCI one-body density matrix | `chem.analyse` |
| FR-05 | Scan over bond lengths, CSV | `scan.scan`, `cli.py` |
| FR-06 | Threshold distance (interpolated), HF/CISD recommendation, rank correlation of the index with the HF energy error | `scan.py` |

```bash
pip install -e ".[dev]"
pytest -q
stretchwatch run --molecule H2 --n 12 --out results
stretchwatch run --config molecule.json            # {"molecule": "LiH", "distances": [1.0, 1.6, 2.5, 3.5]}
```

## Lab result (exact calculations; not a quantum-hardware claim)
* H₂: HF energy −1.1168 Ha and FCI −1.1373 Ha at 0.74 Å (STO-3G), as in the literature. |⟨HF|FCI⟩|² falls from 0.99 to 0.50 at dissociation; entropy rises to 2 ln 2; HF stops being a good start (< 0.9) at 1.37 Å.
* LiH: threshold at 2.31 Å. The quality index is perfectly rank-correlated with the HF energy error on the scans (Spearman −1.00).
* HeH⁺ is a negative control: it dissociates into He + H⁺, HF stays exact and no threshold exists.
* CISD is exact for the two-electron molecules and has overlap² ≥ 0.9999 for LiH in this tiny basis, so the CISD result says little about larger systems.
