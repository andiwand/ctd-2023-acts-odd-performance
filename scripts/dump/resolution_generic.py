#!/usr/bin/env python3

import argparse
import pandas as pd

from mycommon.io import read_residuals
from mycommon.agg import agg_resolution


def dump_resolution(x, y, input, output):
    x_range = {
        "eta": (-3, 3),
        "pt": (1, 100),
    }[x]
    x_bins = {
        "eta": 13,
        "pt": 10,
    }[x]

    data = read_residuals(input)

    (x_mid, (std, std_std)) = agg_resolution(x_range, x_bins, data[x], data[f"res_{y}"])

    pd.DataFrame(
        {
            "x": x_mid,
            "y": std,
            "y_err": std_std,
        }
    ).to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("x", choices=["eta", "pt"])
    parser.add_argument("y", choices=["d0", "z0", "qop"])
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    dump_resolution(args.x, args.y, args.input, args.output)
