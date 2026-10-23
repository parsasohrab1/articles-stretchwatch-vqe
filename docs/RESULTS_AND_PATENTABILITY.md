# StretchWatch — measured results, prior art, patentability (preliminary, NOT legal advice)
Validation: H2 at 1.4 bohr gives E = -1.1373 Ha (literature FCI/STO-3G), test in `tests/`. Results in `results/curves.json` (STO-3G, bohr):
| molecule | R | HF correlation energy (mHa) | HF overlap | mean qubit entropy (bits) |
|---|---|---|---|---|
| H2 | 1.0 / 1.4 / 2.0 / 3.0 / 4.0 / 6.0 | 13 / 21 / 39 / 100 / 183 / 288 | 0.994 / 0.987 / 0.963 / 0.847 / 0.678 / 0.527 | 0.05 / 0.10 / 0.23 / 0.62 / 0.91 / 1.00 |
| HeH+ | 1.0 ... 6.0 | 9.3 -> 0.0 | 0.997 -> 1.000 (dissociates to He + H+, no static correlation) | 0.03 -> 0.00 |
| LiH | 2.0 / 3.0 / 4.0 / 6.0 / 8.0 | 17 / 20 / 34 / 102 / 495 | 0.982 / 0.975 / 0.935 / 0.674 / 0.008 | 0.05 / 0.06 / 0.12 / 0.34 / 0.39 |
Reading: the HF overlap falls below 0.9 near R = 2.63 bohr for H2 and 4.32 bohr for LiH (linear interpolation on the computed grids); beyond that HF is a poor VQE start and entanglement resources rise. Small basis (STO-3G), 3 molecules, no active-space or noise study — **[to be measured]**. The overlap/entropy numbers are exact for the model Hamiltonian, not estimates.
## Prior art (web search, partial; Google Patents API was blocked)
- Literature already shows HF deviation and overlap loss on stretching (CAFQA, arXiv 2202.12924; arXiv 2409.05939; arXiv 1809.05523 "orthogonality catastrophe"); novelty of the observation is **low**; a tidy reproducible table/tool is the contribution.
- Patents: EP3807804A1, US20200364601A1, US20250259093A1 list HF as an initial/reference state; none seen claiming overlap-based switching of the initial state. Patentability: low.
## Must not be disclosed
Nothing from the 37 applications is used; the attorney must check overlap with the HQCA quantum-classical family before submission.
