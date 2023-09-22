#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse

from mycommon.plot_style import myPlotStyle
from mycommon.label import pt_label


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
parser.add_argument("--particles", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 25

for file in args.tracksummary:
    tracksummary = uproot.open(file)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(["t_eta"], library="ak"),
        how="outer",
    )
    tracksummary = tracksummary.dropna()

    eta = tracksummary["t_eta"].values

    plt.hist(eta, bins=eta_bins, range=eta_range, histtype="step", label=pt_label(file))

for file in args.particles:
    particles = uproot.open(file)
    particles = ak.to_dataframe(
        particles["particles"].arrays(["eta"], library="ak"),
        how="outer",
    )
    particles = particles.dropna()

    eta = particles["eta"].values

    plt.hist(eta, bins=eta_bins, range=eta_range, histtype="step", label=file)

plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
