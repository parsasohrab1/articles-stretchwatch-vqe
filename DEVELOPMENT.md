# StretchWatch — development status (branch `develop`)
Implemented from the SRS in README.md (first development pass, 2026-10-23):
| SRS item | Status |
|---|---|
| FR-01 geometry -> qubit Hamiltonian | done (PennyLane built-in chemistry, STO-3G, Jordan-Wigner; geometry along z in bohr) |
| FR-02 exact ground state | done (sparse diagonalisation in the correct electron-number sector; H2 at 1.4 bohr = -1.1373 Ha) |
| FR-03 HF state and overlap | done (`stretchwatch.core`, `stretchwatch.analysis`) |
| FR-04 entanglement | done as **spatial-orbital entropy** (qubit pairs 2k,2k+1) in `analysis.py`; **natural-orbital** entropy not implemented |
| FR-05 distance arrays, CSV, JSON input | done: `python -m stretchwatch configs/example.json --out results/example` |
| FR-06 threshold distance | done: linear interpolation of the HF-overlap crossing (default 0.9) |
| Quality index + reference suggestion (section 4) | prototype: `start_quality_index` and `recommend_reference` (HF / dominant determinant / CISD-type / multireference) — **heuristics, not validated** |
| Performance: H2 in 20 points < 1 min | met (measured 2.1 s) |
Tests: `python -m pytest tests -q` (6 pass).
## Deliberate deviation from the SRS appendix
The appendix proposes a synthetic Morse-model data generator. **It is not used**: those curves are invented, not physics. All results come from computed Hamiltonians (not measured data) and are labelled so.
## Results (results/example/*.csv, summary.json)
Threshold at HF overlap 0.9: H2 about 2.63 bohr, LiH about 4.32 bohr, HeH+ never below 0.9 on the computed grid.
## Next steps
Natural-orbital entropies, larger molecules (BeH2, H2O) and active spaces, bases beyond STO-3G, noise, plotting (FR in SRS: "draw curves"), CI check against an independent chemistry code (PySCF) when available, validation of the recommendation rule.
