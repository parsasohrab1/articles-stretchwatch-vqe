import argparse
import json
import pathlib

from . import scan


def main(argv=None):
    ap = argparse.ArgumentParser(prog="stretchwatch")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="scan bond lengths and judge the Hartree-Fock starting point")
    r.add_argument("--molecule", default="H2", choices=["H2", "HeH+", "LiH"])
    r.add_argument("--n", type=int, default=12, help="number of bond lengths")
    r.add_argument("--config", help='JSON {"molecule": "H2", "distances": [0.5, 0.74, 1.5]} (overrides --molecule/--n)')
    r.add_argument("--tau", type=float, default=0.9, help="HF is a good start while |<HF|FCI>|^2 >= tau")
    r.add_argument("--out", default="results")
    a = ap.parse_args(argv)
    mol, dist = a.molecule, None
    if a.config:
        c = json.load(open(a.config))
        mol, dist = c["molecule"], c.get("distances")
    df = scan.recommend_start(scan.scan(mol, dist, a.n), a.tau)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f"{mol}.csv", index=False)
    print(df[["R", "e_hf", "e_fci", "overlap2", "entropy", "quality", "cisd_overlap2", "start"]].round(4).to_string(index=False))
    print(f"equilibrium {scan.equilibrium(df):.2f} A; HF stops being a good start (overlap^2 < {a.tau}) at {scan.threshold_distance(df, a.tau):.2f} A; "
          f"Spearman(quality, HF energy error) = {scan.index_vs_error(df):.2f}")


if __name__ == "__main__":
    main()
