#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse

from mycommon.plot_style import myPlotStyle


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (-3, 3)
pt_range = (0, 100)
eta_bins = 11
pt_bins = 11

particles = uproot.open(args.input)
particles = ak.to_dataframe(
    particles["particles"].arrays(["m", "q", "p", "pt", "phi", "eta"], library="ak"),
    how="outer",
).dropna()
particles = particles.query("q != 0")

axes = plt.gcf().subplots(2, 3)

axes[0, 0].hist(particles["eta"], bins=eta_bins, range=eta_range, histtype="step")
axes[0, 0].set_xlabel(r"$\eta$")
axes[0, 0].set_xlim(*eta_range)

axes[0, 1].hist(particles["pt"], bins=pt_bins, range=pt_range, histtype="step")
axes[0, 1].set_xlabel(r"$p_T$")
axes[0, 1].set_xlim(*pt_range)

_, _, _, im = axes[0, 2].hist2d(
    particles["eta"],
    particles["pt"],
    bins=(eta_bins, pt_bins),
    range=(eta_range, pt_range),
    cmap="Oranges",
)
axes[0, 2].set_xlabel(r"$\eta$")
axes[0, 2].set_ylabel(r"$p_T$")
axes[0, 2].set_xlim(*eta_range)
axes[0, 2].set_ylim(*pt_range)
plt.gcf().colorbar(im, ax=axes[0, 2])

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
