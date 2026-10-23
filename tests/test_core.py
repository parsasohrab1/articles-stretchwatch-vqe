import numpy as np

from stretchwatch import hf_index, mean_single_qubit_entropy, point


def test_h2_equilibrium_matches_known_sto3g_energy():
    r = point(["H", "H"], [0, 0, 0, 0, 0, 1.4], electrons=2)       # R = 1.4 bohr ~ 0.74 A; literature FCI(STO-3G) ~ -1.1373 Ha
    assert abs(r["E_exact"] - (-1.1373)) < 2e-3 and r["hf_overlap"] > 0.95 and r["corr_energy_mHa"] > 0


def test_helpers():
    assert hf_index(4, 2) == 0b1100
    bell = np.zeros(4); bell[0] = bell[3] = 1 / np.sqrt(2)
    assert abs(mean_single_qubit_entropy(bell, 2) - 1.0) < 1e-9
    prod = np.zeros(4); prod[0] = 1; assert mean_single_qubit_entropy(prod, 2) < 1e-6
