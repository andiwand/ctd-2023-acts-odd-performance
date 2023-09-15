#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse
import numpy as np
from scipy.stats import binned_statistic

from common.plot_style import myPlotStyle


def label(file):
    particle_map = {
        "mu": "$\mu$",
        "pi": "$\pi$",
        "e": "$e$",
    }
    split = Path(file).parent.name.split("_")
    #return f"single {particle_map[split[0]]} {split[1].replace('GeV', '')} GeV"
    return f"{split[1].replace('GeV', '')} GeV"


def smoothed_std(data):
    def fit2(data):
        return np.mean(data), np.std(data)

    for _ in range(3):
        m, s = fit2(data)
        data = data[np.abs(data - m) < 3 * s]
    return s


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
args = parser.parse_args()

abs_eta_range = (0, 3)
abs_eta_bins = 15

for file in args.tracksummary:
    tracksummary = uproot.open(file)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(["t_eta", "res_eLOC0_fit"], library="ak"),
        how="outer",
    )
    tracksummary = tracksummary.dropna()

    abs_eta = np.abs(tracksummary["t_eta"].values)
    d0 = tracksummary["res_eLOC0_fit"].values

    resolution_d0_binned, abs_eta_edges, _ = binned_statistic(
        abs_eta, d0, bins=abs_eta_bins, range=abs_eta_range, statistic=smoothed_std
    )
    abs_eta_mid = 0.5 * (abs_eta_edges[:-1] + abs_eta_edges[1:])

    # debugging distribution of eta bins
    """
    digi = np.digitize(abs_eta, abs_eta_edges)
    for i in range(abs_eta_bins):
        f = plt.figure(f"{i} = {abs_eta_edges[i]:.2f} < |eta| < {abs_eta_edges[i+1]:.2f}")
        plt.hist(d0[digi == i], bins=100, range=(-0.4, 0.4), histtype="step")
    """

    plt.figure(0)
    plt.plot(abs_eta_mid, resolution_d0_binned, marker="o", label=label(file))

plt.title("Resolution of $d_0$ over $|\eta|$ for single $\mu$ events")
plt.xlabel("$|\eta|$")
plt.ylabel("$\sigma(d_0)$ [mm]")
plt.xticks(np.linspace(0, 3, 7))
plt.yticks(np.linspace(0, 0.3, 6))
plt.legend()
plt.show()
