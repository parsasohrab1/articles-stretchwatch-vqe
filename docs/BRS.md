# Business Requirements Specification (BRS)

**Product:** StretchWatch (articles-stretchwatch-vqe)
*Reconstructed from the README, the requirements document, the TRL assessment and the tests of this repository. Not reviewed by the owner. Do not publish method details covered by filed claims.*

## 1. Business problem
The Hartree-Fock state is the usual starting point of the variational quantum eigensolver (VQE), but it becomes a poor start as a bond is stretched. Users need to know where.

## 2. Business objectives
| ID | Objective |
|----|-----------|
| BO-1 | Compute a quality index of the Hartree-Fock reference along a bond-stretching curve |
| BO-2 | Find the bond length beyond which HF stops being a good start |
| BO-3 | Recommend HF or an alternative reference (CISD) |

## 3. Users and stakeholders (proposed; to be confirmed by the owner)
Quantum-chemistry and quantum-algorithm researchers.

## 4. Scope
In scope: exact classical calculations (STO-3G) for H2, HeH+ and LiH.
Out of scope: quantum hardware, circuit execution, larger basis sets or molecules.

## 5. Constraints and assumptions
- All data are synthetic or public benchmark data; the product makes no clinical, operational or regulatory claim.
- Maturity is TRL 4 (see docs/TRL.md).

## 6. Success measures (to be measured on real data)
Rank correlation of the index with VQE iteration counts or convergence failures on simulators and devices.
