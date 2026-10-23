"""Command line: python -m stretchwatch configs/example.json  ->  one CSV per molecule plus summary.json (FR-05, FR-06)."""

import argparse
import csv
import json
from pathlib import Path

from .analysis import analyse_point, find_threshold


def run_config(cfg: dict, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    thr = float(cfg.get("overlap_threshold", 0.9))
    summary = {}
    for mol in cfg["molecules"]:
        rows = [analyse_point(mol["symbols"], float(r), int(mol["electrons"]), int(mol.get("charge", 0)), thr) for r in mol["distances_bohr"]]
        with open(out_dir / f"{mol['name']}_curve.csv", "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        summary[mol["name"]] = {"hf_overlap_threshold": thr, "threshold_R_bohr": find_threshold([r["R_bohr"] for r in rows], [r["hf_overlap"] for r in rows], thr),
                                "n_points": len(rows)}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf8")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(prog="stretchwatch")
    ap.add_argument("config")
    ap.add_argument("--out", default="results/run")
    a = ap.parse_args(argv)
    s = run_config(json.loads(Path(a.config).read_text(encoding="utf8")), Path(a.out))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
