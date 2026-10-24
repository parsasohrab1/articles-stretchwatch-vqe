import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import chem


def scan(name, distances=None, n=12):
    lo, hi = chem.MOLECULES[name][2]
    Rs = np.asarray(distances if distances is not None else np.linspace(lo, hi, n), float)
    return pd.DataFrame([chem.analyse(name, float(R)) for R in Rs])


def equilibrium(df):
    return float(df.loc[df.e_fci.idxmin(), "R"])


def threshold_distance(df, tau=0.9, key="overlap2"):
    """First bond length beyond the equilibrium at which `key` falls below tau (linear interpolation between grid points); NaN if it never does."""
    d = df[df.R >= equilibrium(df)].sort_values("R")
    v, R = d[key].values, d.R.values
    for i in range(1, len(v)):
        if v[i] < tau <= v[i - 1]:
            return float(R[i - 1] + (R[i] - R[i - 1]) * (v[i - 1] - tau) / (v[i - 1] - v[i]))
    return float("nan") if v[0] >= tau else float(R[0])


def recommend_start(df, tau=0.9):
    """HF where its overlap with the ground state is at least tau, otherwise an alternative reference (CISD)."""
    d = df.copy()
    d["start"] = np.where(d.overlap2 >= tau, "HF", "CISD")
    return d


def index_vs_error(df):
    """Spearman correlation between the starting-point quality index and the HF energy error (should be strongly negative)."""
    return float(spearmanr(df.quality, df.corr_energy).statistic)
