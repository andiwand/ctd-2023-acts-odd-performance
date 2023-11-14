#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse

from mycommon.plot_style import myPlotStyle
from mycommon.label import pt_label
from scipy.stats import binned_statistic


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 25

for file in args.tracksummary:
    tracksummary = uproot.open(file)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(["t_eta", "nHoles"], library="ak"),
        how="outer",
    )
    tracksummary = tracksummary.dropna()

    eta = tracksummary["t_eta"].values
    data = tracksummary["nHoles"].values

    mean, eta_edges, _ = binned_statistic(
        eta, data, bins=eta_bins, range=eta_range, statistic="mean"
    )
    std, _, _ = binned_statistic(
        eta, data, bins=eta_bins, range=eta_range, statistic="std"
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    plt.errorbar(eta_mid, mean, std, fmt="o", label=pt_label(file))

plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
