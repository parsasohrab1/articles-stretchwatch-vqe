# TRL assessment – StretchWatch

Level reached: **TRL 4 – component/subsystem validated in the laboratory (exact classical calculation).**

| Evidence | Status |
|----------|--------|
| Working prototype: PySCF-based analysis, scan, threshold, recommendation, CLI | done |
| Automated tests: H₂ energies against literature values, variational principle, CISD exact for two electrons, overlap falls and entropy rises with stretching (limit 2 ln 2), HeH⁺ negative control, interpolated threshold, HF/CISD switch, index vs energy error on LiH, CLI | `pytest` – 8 tests, run in CI |
| Validation against known ground truth | full CI is the exact answer in the basis |

Gaps before TRL 5: actual VQE runs on simulators and noisy devices, bigger basis sets and molecules, a rule validated against convergence behaviour, a PennyLane/Qiskit Hamiltonian export.
