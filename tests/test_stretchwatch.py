import json

import numpy as np
import pytest

from stretchwatch import chem, cli, scan


def test_h2_equilibrium_energies_match_literature():
    d = chem.analyse("H2", 0.74)
    assert d["e_hf"] == pytest.approx(-1.1168, abs=2e-3) and d["e_fci"] == pytest.approx(-1.1373, abs=2e-3)   # STO-3G
    assert d["overlap2"] > 0.98 and d["hf_converged"]


def test_variational_principle_and_cisd_exact_for_two_electrons():
    for R in (0.5, 1.5, 3.0):
        d = chem.analyse("H2", R)
        assert d["e_fci"] <= d["e_hf"] + 1e-9 and d["cisd_overlap2"] > 0.999999 and abs(d["cisd_error"]) < 1e-8


def test_overlap_falls_and_entropy_grows_with_stretching():
    df = scan.scan("H2", [0.5, 0.74, 1.5, 2.5, 4.0])
    assert (np.diff(df.overlap2) < 0).all() and (np.diff(df.entropy) > 0).all() and (np.diff(df.quality) < 0).all()
    assert df.entropy.iloc[-1] == pytest.approx(2 * np.log(2), abs=0.05)         # two fully fractional natural orbitals at dissociation


def test_heh_plus_is_a_negative_control():
    df = scan.scan("HeH+", [0.77, 2.5, 4.0])
    assert df.overlap2.min() > 0.99 and np.isnan(scan.threshold_distance(df, 0.9))   # dissociates into He + H+, HF stays a good start


def test_threshold_distance_is_interpolated_and_ordered():
    df = scan.scan("H2", [0.4, 0.74, 1.0, 1.5, 2.0, 3.0])
    t90, t80 = scan.threshold_distance(df, 0.9), scan.threshold_distance(df, 0.8)
    assert 0.74 < t90 < 1.5 and t90 < t80 and scan.equilibrium(df) == pytest.approx(0.74)


def test_recommend_start_switches_beyond_threshold():
    df = scan.recommend_start(scan.scan("H2", [0.74, 3.0]), 0.9)
    assert list(df.start) == ["HF", "CISD"]


def test_quality_index_tracks_hf_energy_error():
    df = scan.scan("LiH", [1.0, 1.6, 2.5, 3.5, 5.0])
    assert scan.index_vs_error(df) < -0.8


def test_cli_with_config(tmp_path):
    cfg = tmp_path / "c.json"
    cfg.write_text(json.dumps({"molecule": "H2", "distances": [0.74, 2.0]}))
    cli.main(["run", "--config", str(cfg), "--out", str(tmp_path)])
    assert (tmp_path / "H2.csv").exists()
