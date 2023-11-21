#!/usr/bin/env python3

import pandas as pd
import argparse

from mycommon.io import read_pulls
from mycommon.agg import agg_pulls_over_eta


def dump_pulls_over_eta_sausage(input, output, eta_range, eta_bins):
    pull_labels = ["d0", "z0", "t", "phi", "theta", "qop"]

    eta, pulls = read_pulls(input)

    data = {}

    for pull_label, pull in zip(pull_labels, pulls):
        (
            eta_mid,
            (
                mean_binned,
                std_binned,
            ),
        ) = agg_pulls_over_eta(eta_range, eta_bins, eta, pull)

        data["eta"] = eta_mid
        data[pull_label + "_mean"] = mean_binned
        data[pull_label + "_std"] = std_binned

    pd.DataFrame(data).to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--eta-range", nargs=2, type=float, default=(0, 3))
    parser.add_argument("--eta-bins", type=int, default=8)
    args = parser.parse_args()

    dump_pulls_over_eta_sausage(args.input, args.output, args.eta_range, args.eta_bins)
