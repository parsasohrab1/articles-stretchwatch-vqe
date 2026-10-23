"""Measured (exact diagonalisation of computed STO-3G Hamiltonians): HF quality and entanglement along stretching curves of H2, HeH+ and LiH."""
import json
from pathlib import Path

import numpy as np

from stretchwatch import point

MOL = {"H2": (["H", "H"], 2, 0, [1.0, 1.4, 2.0, 3.0, 4.0, 6.0]),
       "HeH+": (["He", "H"], 2, 1, [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]),
       "LiH": (["Li", "H"], 4, 0, [2.0, 3.0, 4.0, 6.0, 8.0])}


def run():
    out = {}
    for name, (sym, ne, ch, rs) in MOL.items():
        out[name] = {}
        for r in rs:
            out[name][str(r)] = point(sym, [0, 0, 0, 0, 0, r], ne, ch)
            print(name, r, {k: round(v, 4) for k, v in out[name][str(r)].items()}, flush=True)
    return out


if __name__ == "__main__":
    r = run(); Path("results").mkdir(exist_ok=True); Path("results/curves.json").write_text(json.dumps(r, indent=1))
