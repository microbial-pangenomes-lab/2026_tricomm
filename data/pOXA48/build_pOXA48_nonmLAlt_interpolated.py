#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PLATE_COUNTS = HERE / "plate_counts.tsv"
DEFAULT_OUT = HERE.parent / "model" / "pOXA48_nonmLAlt_interpolated.csv"

PSEUDOCOUNT = 5e-6
SCALE = 1e7

STRAINS = ["LM", "PM", "PL"]
NODRUG_VIALS = ["M0", "M2", "M4", "M6"]
DRUG_VIALS = ["M1", "M3", "M5", "M7"]


def build(plate_counts_path: Path) -> pd.DataFrame:
    df = pd.read_csv(plate_counts_path, sep="\t")

    sub = df[(df["plate"] == "KAN") & (df["replicate"] == "rep2")].copy()
    sub["val"] = sub["CFU/mL"] / SCALE + PSEUDOCOUNT

    def series(vial, strain):
        s = sub[(sub["vial"] == vial) & (sub["strain"] == strain)].sort_values("time")
        return s["time"].values, s["val"].values

    actual_times = sorted(sub["time"].unique())
    midpoints = [(a + b) / 2 for a, b in zip(actual_times[:-1], actual_times[1:])]
    grid = np.array(sorted(set(actual_times) | set(midpoints)))

    out = pd.DataFrame({"time": grid})

    for strain in STRAINS:
        nodrug_matrix = np.array(
            [np.interp(grid, *series(vial, strain)) for vial in NODRUG_VIALS]
        )
        out[f"nodrug_mean_{strain}"] = nodrug_matrix.mean(axis=0)
        for vial in DRUG_VIALS:
            out[f"{vial}_rep2_{strain}"] = np.interp(grid, *series(vial, strain))

    for strain in STRAINS:
        nodrug_matrix = np.array(
            [np.interp(grid, *series(vial, strain)) for vial in NODRUG_VIALS]
        )
        mean = nodrug_matrix.mean(axis=0)
        std = nodrug_matrix.std(axis=0, ddof=1)
        out[f"std_{strain}"] = std
        out[f"lower_{strain}"] = np.minimum(std, mean)

    # match the exact column order of the original file
    col_order = (
        ["time"]
        + [f"nodrug_mean_{s}" for s in STRAINS]
        + [f"{v}_rep2_{s}" for v in DRUG_VIALS for s in STRAINS]
        + [f"std_{s}" for s in STRAINS]
        + [f"lower_{s}" for s in STRAINS]
    )
    return out[col_order]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--plate-counts", type=Path, default=PLATE_COUNTS)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--verify",
        action="store_true",
        help="compare the rebuilt table against --out instead of writing it",
    )
    args = ap.parse_args()

    result = build(args.plate_counts)

    if args.verify:
        target = pd.read_csv(args.out)
        diff = (result - target).abs().max().max()
        print(f"max abs difference vs {args.out}: {diff:.3e}")
        assert diff < 1e-9, "reconstruction does not match target file"
        print("OK: reconstruction matches the existing file.")
    else:
        result.to_csv(args.out, index=False)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
