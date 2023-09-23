#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import uproot
import awkward as ak
import argparse
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.stats import smoothed_mean, smoothed_std
from mycommon.label import particle_label, pt_label, param_label


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

columns = [
    "pull_eLOC0_fit",
    "pull_eLOC1_fit",
    # "pull_eT_fit",
    # "pull_ePHI_fit",
    # "pull_eTHETA_fit",
    "pull_eQOP_fit",
]
eta_range = (0, 3)
pull_range = (-4, 4)
eta_bins = 7

fig = plt.figure(figsize=(16, 8))
axs = fig.subplots(1, 3, sharey=True)

for i, file in enumerate(args.tracksummary):
    tracksummary = uproot.open(file)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(
            columns + ["t_eta", "nMeasurements"], library="ak"
        ),
        how="outer",
    )
    tracksummary = tracksummary.dropna()
    tracksummary = tracksummary.query("nMeasurements >= 10")

    for j, (col, ax) in enumerate(
        zip(
            columns,
            axs.flat,
        )
    ):
        ax.set_title(param_label(col))
        ax.set_xlim(eta_range[0], eta_range[1] + 0.2)
        ax.set_ylim(*pull_range)

        eta = tracksummary["t_eta"].values
        data = tracksummary[col].values

        mean_binned, eta_edges, _ = binned_statistic(
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
        eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

        trans = Affine2D().translate(0.1 * i, 0.0) + ax.transData
        ax.errorbar(
            eta_mid,
            mean_binned,
            yerr=std_binned,
            fmt="o",
            transform=trans,
            label=pt_label(file) if j == 0 else None,
        )

        ax.axhline(0, linestyle="--", color="gray")
        ax.axhline(-1, linestyle="--", color="gray")
        ax.axhline(+1, linestyle="--", color="gray")

fig.suptitle(f"single {particle_label(args.tracksummary[0])} pulls over $\eta$")
fig.supxlabel(r"$\eta$")
fig.supylabel(r"pull")
fig.legend()

if args.output:
    fig.savefig(args.output)
else:
    plt.show()
