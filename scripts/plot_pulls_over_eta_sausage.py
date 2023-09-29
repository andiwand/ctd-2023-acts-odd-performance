#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import pandas as pd
import argparse
from scipy.stats import binned_statistic

from mycommon.stats import smoothed_mean, smoothed_std


def get_data(file):
    if str(file).endswith(".root"):
        columns = [
            "pull_eLOC0_fit",
            "pull_eLOC1_fit",
            "pull_eT_fit",
            "pull_ePHI_fit",
            "pull_eTHETA_fit",
            "pull_eQOP_fit",
        ]

        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(
                columns + ["t_eta", "nMeasurements"], library="ak"
            ),
            how="outer",
        ).dropna()
        tracksummary = tracksummary.query("nMeasurements >= 10")

        eta = tracksummary["t_eta"].values
        pulls = [tracksummary[col].values for col in columns]

        return eta, pulls

    if str(file).endswith(".csv"):
        columns = [
            "track_pull_eLOC0_fit",
            "track_pull_eLOC1_fit",
            "track_pull_eT_fit",
            "track_pull_ePHI_fit",
            "track_pull_eTHETA_fit",
            "track_pull_eQOP_fit",
        ]

        data = pd.read_csv(file).dropna()

        eta = data["true_eta"].values
        pulls = [data[col].values for col in columns]

        return eta, pulls

    raise ValueError(f"unknown file type: {file}")


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (-3, 3)
pull_range = (-4, 4)
eta_bins = 40
pull_bins = 40

pull_labels = [
    r"$d_0$",
    r"$z_0$",
    r"$t$",
    r"$\phi$",
    r"$\theta$",
    r"$\frac{q}{p}$",
]

eta, pulls = get_data(args.input)

subfigs = plt.gcf().subfigures(2, 3)

for pull_label, pull, subfig in zip(
    pull_labels,
    pulls,
    subfigs.flat,
):
    ax, ax_mean, ax_std = subfig.subplots(
        3,
        1,
        gridspec_kw={"height_ratios": [5, 1, 1], "hspace": 0},
        sharex=True,
    )

    ax.set_title(pull_label)
    ax.set_xlim(eta_range[0], eta_range[1] + 0.2)
    ax.set_ylim(*pull_range)

    h, eta_edges, data_edges, im = ax.hist2d(
        eta,
        pull,
        range=(eta_range, pull_range),
        bins=(eta_bins, pull_bins),
        density=True,
        cmap="Oranges",
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    mean_binned, _, other_digi = binned_statistic(
        eta,
        pull,
        bins=eta_bins,
        range=eta_range,
        statistic=smoothed_mean,
    )
    std_binned, _, _ = binned_statistic(
        eta,
        pull,
        bins=eta_bins,
        range=eta_range,
        statistic=smoothed_std,
    )

    ax.plot(eta_mid, mean_binned, linestyle="-", color="black", label="fit")
    ax.plot(eta_mid, mean_binned - std_binned, linestyle="--", color="black")
    ax.plot(eta_mid, mean_binned + std_binned, linestyle="--", color="black")

    ax_mean.plot(eta_mid, mean_binned, linestyle="-", color="black")
    ax_mean.axhline(0, linestyle="--", color="gray")

    ax_std.plot(eta_mid, std_binned, linestyle="-", color="black")
    ax_std.axhline(1, linestyle="--", color="gray")

    plt.colorbar(im, ax=subfig.get_axes())

plt.supxlabel(r"$|\eta|$")
plt.supylabel(r"pull")
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
