import json
from pathlib import Path

from validate_reference import BOHR_PER_A, minimum


def test_nist_reference_file_and_h2_equilibrium_close_to_experiment():
    ref = json.loads((Path(__file__).resolve().parents[1] / "data" / "nist_reference.json").read_text())
    assert abs(ref["H2"]["r_e_angstrom"] - 0.74144) < 1e-9
    r0, e0 = minimum(["H", "H"], 2, 1.0, 2.0, 9)
    assert abs(r0 / BOHR_PER_A - ref["H2"]["r_e_angstrom"]) / ref["H2"]["r_e_angstrom"] < 0.03 and e0 < -1.13
