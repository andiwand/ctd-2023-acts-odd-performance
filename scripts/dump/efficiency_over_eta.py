#!/usr/bin/env python3

import pandas as pd
import argparse

from mycommon.agg import agg_efficiency_over_eta
from mycommon.io import read_track_efficiency


def dump_efficiency_over_eta(input, output):
    eta_range = (-3, 3)
    eta_bins = 30

    eta, track_efficiency = read_track_efficiency(input)

    (
        eta_mid,
        (
            track_efficiency_mean,
            (
                lower_error,
                upper_error,
            ),
        ),
    ) = agg_efficiency_over_eta(eta_range, eta_bins, eta, track_efficiency)

    pd.DataFrame(
        {
            "eta": eta_mid,
            "efficiency": track_efficiency_mean,
            "lower_error": lower_error,
            "upper_error": upper_error,
        }
    ).to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    dump_efficiency_over_eta(args.input, args.output)
