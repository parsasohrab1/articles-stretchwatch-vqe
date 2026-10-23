"""Validate the computed STO-3G exact curves against EXPERIMENTAL reference constants (NIST WebBook): equilibrium bond length, and H2 dissociation energy.
STO-3G is a minimal basis, so differences of a few percent in r_e and tens of percent in D are expected; the check shows the physics is sane, not that the basis is accurate."""
import json
from pathlib import Path

import numpy as np

from stretchwatch.core import ground_state_in_sector, hamiltonian_matrix

BOHR_PER_A = 1.8897261246
HARTREE_EV = 27.211386


def energy(symbols, r, electrons):
    Hs, n = hamiltonian_matrix(symbols, [0, 0, 0, 0, 0, r], 0)
    return ground_state_in_sector(Hs, n, electrons)[0]


def minimum(symbols, electrons, lo, hi, n=15):
    rs = np.linspace(lo, hi, n)
    es = np.array([energy(symbols, r, electrons) for r in rs])
    i = int(np.argmin(es)); i = min(max(i, 1), n - 2)
    a, b, c = np.polyfit(rs[i - 1:i + 2], es[i - 1:i + 2], 2)
    r0 = -b / (2 * a)
    return float(r0), float(a * r0 ** 2 + b * r0 + c)


def run():
    ref = json.loads((Path(__file__).parent / "data" / "nist_reference.json").read_text())
    out = {"reference": ref["source"]}
    r0, e0 = minimum(["H", "H"], 2, 1.0, 2.0)
    e_far = energy(["H", "H"], 12.0, 2)
    out["H2"] = {"r_e_computed_A": r0 / BOHR_PER_A, "r_e_experiment_A": ref["H2"]["r_e_angstrom"], "r_e_error_percent": 100 * (r0 / BOHR_PER_A / ref["H2"]["r_e_angstrom"] - 1),
                 "D_computed_eV": (e_far - e0) * HARTREE_EV, "D0_experiment_eV": ref["H2"]["D0_cm-1"] * 1.239842e-4, "note": "computed D is De at the STO-3G minimum vs experimental D0 (zero-point energy ~0.27 eV not subtracted)"}
    r1, _ = minimum(["Li", "H"], 4, 2.6, 3.8, 13)
    out["LiH"] = {"r_e_computed_A": r1 / BOHR_PER_A, "r_e_experiment_A": ref["LiH"]["r_e_angstrom"], "r_e_error_percent": 100 * (r1 / BOHR_PER_A / ref["LiH"]["r_e_angstrom"] - 1)}
    return out


if __name__ == "__main__":
    r = run()
    Path("results").mkdir(exist_ok=True)
    Path("results/validation_nist.json").write_text(json.dumps(r, indent=1))
    print(json.dumps(r, indent=1))
