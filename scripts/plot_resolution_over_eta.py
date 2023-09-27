#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse
import numpy as np
import pandas as pd
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.stats import smoothed_std
from mycommon.events import split_event_label
from mycommon.label import get_event_variant_label, get_event_type_label
from mycommon.paths import get_event_label_from_path


def get_data(file):
    if str(file).endswith(".root"):
        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(
                ["t_eta", "res_eLOC0_fit"], library="ak"
            ),
            how="outer",
        ).dropna()

        eta = tracksummary["t_eta"].values
        d0 = tracksummary["res_eLOC0_fit"].values

        return eta, d0

    if str(file).endswith(".csv"):
        data = pd.read_csv(file).dropna()

        eta = data["true_eta"].values
        d0 = data["track_res_eLOC0_fit"].values

        return eta, d0

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 31

for file in args.input:
    event_label = get_event_label_from_path(file)
    event, _ = split_event_label(event_label)

    eta, d0 = get_data(file)

    resolution_d0_binned, eta_edges, _ = binned_statistic(
        eta, d0, bins=eta_bins, range=eta_range, statistic=smoothed_std
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    # debugging distribution of eta bins
    """
    digi = np.digitize(eta, eta_edges)
    for i in range(eta_bins):
        f = plt.figure(f"{i} = {eta_edges[i]:.2f} < |eta| < {eta_edges[i+1]:.2f}")
        plt.hist(d0[digi == i], bins=100, range=(-0.4, 0.4), histtype="step")
    """

    plt.plot(
        eta_mid,
        resolution_d0_binned,
        marker="o",
        linestyle="",
        label=get_event_variant_label(event),
    )

plt.title(f"Resolution of $d_0$ over $\eta$ for {get_event_type_label(event)} events")
plt.xlabel("$\eta$")
plt.ylabel("$\sigma(d_0)$ [mm]")
plt.xticks(np.linspace(*eta_range, 10))
plt.yticks(np.linspace(0, 0.3, 6))
plt.xlim(eta_range)
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
