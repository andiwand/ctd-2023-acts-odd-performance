#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse
from scipy.stats import binned_statistic

from mycommon.stats import smoothed_mean, smoothed_std


parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
args = parser.parse_args()

columns = [
    "pull_eLOC0_fit",
    "pull_eLOC1_fit",
    "pull_eT_fit",
    "pull_ePHI_fit",
    "pull_eTHETA_fit",
    "pull_eQOP_fit",
]
eta_range = (-3, 3)
pull_range = (-4, 4)
eta_bins = 40
pull_bins = 40

for file in args.tracksummary:
    tracksummary = uproot.open(file)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(columns + ["t_eta", "nMeasurements"], library="ak"),
        how="outer",
    )
    tracksummary = tracksummary.dropna()
    tracksummary = tracksummary.query("nMeasurements >= 10")

    fig = plt.figure(file, figsize=(8, 6))
    subfigs = fig.subfigures(2, 3)

    for col, subfig in zip(
        columns,
        subfigs.flat,
    ):
        ax, ax_mean, ax_std = subfig.subplots(
            3,
            1,
            gridspec_kw={"height_ratios": [5, 1, 1], "hspace": 0},
            sharex=True,
        )

        ax.set_title(col)
        ax.set_xlim(eta_range[0], eta_range[1]+0.2)
        ax.set_ylim(*pull_range)

        eta = tracksummary["t_eta"].values
        data = tracksummary[col].values

        h, eta_edges, data_edges, im = ax.hist2d(
            eta,
            data,
            range=(eta_range, pull_range),
            bins=(eta_bins, pull_bins),
            density=True,
            cmap="Oranges",
        )
        eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

        mean_binned, _, other_digi = binned_statistic(
            eta,
            data,
            bins=eta_bins,
            range=eta_range,
            statistic=smoothed_mean,
        )
        std_binned, _, _ = binned_statistic(
            eta,
            data,
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

        fig.colorbar(im, ax=subfig.get_axes())

fig.supxlabel(r"$|\eta|$")
fig.supylabel(r"pull")
fig.legend()
fig.savefig(Path(__file__).parent.parent / "plots/pulls_over_eta_sausage.png")
plt.show()
