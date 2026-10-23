import numpy as np

from stretchwatch.analysis import analyse_point, cisd_weight, excitation_levels, find_threshold, orbital_entropies, recommend_reference, start_quality_index


def test_excitation_levels_and_cisd_weight():
    lev = excitation_levels(4, 2)
    assert lev[0b1100] == 0 and lev[0b0011] == 2 and lev[0b1010] == 1          # HF, double, single (wire 0 = MSB)
    psi = np.zeros(16)
    psi[0b1100] = 0.8
    psi[0b0011] = 0.6
    assert abs(cisd_weight(psi, 4, 2) - 1.0) < 1e-12


def test_orbital_entropy_bounds_and_quality_index():
    hf = np.zeros(16)
    hf[0b1100] = 1.0
    assert max(orbital_entropies(hf, 4)) < 1e-9
    assert start_quality_index(1.0, 0.0) == 1.0 and start_quality_index(0.5, 2.0) == 0.0


def test_threshold_interpolation_and_recommendation():
    assert abs(find_threshold([1, 2, 3], [1.0, 0.95, 0.85], 0.9) - 2.5) < 1e-12
    assert find_threshold([1, 2], [0.99, 0.95], 0.9) is None and find_threshold([1, 2], [0.5, 0.4], 0.9) == 1.0
    assert recommend_reference(0.95, 0.95, 1, 0.9) == "HF" and recommend_reference(0.4, 0.4, 1.0, 0.9) == "CISD-type reference"
    assert recommend_reference(0.4, 0.4, 0.5, 0.9) == "multireference / adaptive ansatz"


def test_h2_equilibrium_and_stretched_behaviour():
    a = analyse_point(["H", "H"], 1.4, 2)
    b = analyse_point(["H", "H"], 6.0, 2)
    assert abs(a["E_exact"] + 1.1373) < 2e-3 and a["recommended_reference"] == "HF" and a["start_quality"] > b["start_quality"]
    assert b["hf_overlap"] < 0.6 and b["cisd_weight"] > 0.999 and b["max_orbital_entropy"] > a["max_orbital_entropy"]
