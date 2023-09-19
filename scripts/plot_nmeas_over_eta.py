#!/usr/bin/env python3

from pathlib import Path
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
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 25

for file in args.tracksummary:
    tracksummary = uproot.open(file)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(["t_eta", "nMeasurements"], library="ak"),
        how="outer",
    )
    tracksummary = tracksummary.dropna()

    eta = tracksummary["t_eta"].values
    data = tracksummary["nMeasurements"].values

    mean, eta_edges, _ = binned_statistic(
        eta, data, bins=eta_bins, range=eta_range, statistic="mean"
    )
    std, _, _ = binned_statistic(
        eta, data, bins=eta_bins, range=eta_range, statistic="std"
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    plt.errorbar(eta_mid, mean, std, fmt="o", label=pt_label(file))

plt.legend()
plt.savefig(Path(__file__).parent.parent / "plots/efficiency_over_eta.png")
plt.show()
