#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import uproot
import awkward as ak
import argparse
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.stats import smoothed_mean, smoothed_std
from mycommon.events import split_event_label
from mycommon.label import (
    get_event_variant_label,
    get_event_type_label,
)
from mycommon.paths import get_event_label_from_path


def get_data(file):
    if str(file).endswith(".root"):
        columns = [
            "pull_eLOC0_fit",
            "pull_eLOC1_fit",
            # "pull_eT_fit",
            # "pull_ePHI_fit",
            # "pull_eTHETA_fit",
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
            # "track_pull_eT_fit",
            # "track_pull_ePHI_fit",
            # "track_pull_eTHETA_fit",
            "track_pull_eQOP_fit",
        ]

        data = pd.read_csv(file).dropna()

        eta = data["true_eta"].values
        pulls = [data[col].values for col in columns]

        return eta, pulls

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (0, 3)
pull_range = (-4, 4)
eta_bins = 7

pull_labels = [
    r"$d_0$",
    r"$z_0$",
    # r"$t$",
    # r"$\phi$",
    # r"$\theta$",
    r"$\frac{q}{p}$",
]

fig = plt.figure(figsize=(16, 8))
axs = fig.subplots(1, 3, sharey=True)

for i, file in enumerate(args.tracksummary):
    event_label = get_event_label_from_path(file)
    event, _ = split_event_label(event_label)

    eta, pulls = get_data(file)

    for j, (pull_label, pull, ax) in enumerate(
        zip(
            pull_labels,
            pulls,
            axs.flat,
        )
    ):
        ax.set_title(pull_label)
        ax.set_xlim(eta_range[0], eta_range[1] + 0.2)
        ax.set_ylim(*pull_range)

        mean_binned, eta_edges, _ = binned_statistic(
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
        eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

        trans = Affine2D().translate(0.1 * i, 0.0) + ax.transData
        ax.errorbar(
            eta_mid,
            mean_binned,
            yerr=std_binned,
            fmt="o",
            transform=trans,
            label=get_event_variant_label(event) if j == 0 else None,
        )

        ax.axhline(0, linestyle="--", color="gray")
        ax.axhline(-1, linestyle="--", color="gray")
        ax.axhline(+1, linestyle="--", color="gray")

fig.suptitle(f"{get_event_type_label(event)} pulls over $\eta$")
fig.supxlabel(r"$\eta$")
fig.supylabel(r"pull")
fig.legend()

if args.output:
    fig.savefig(args.output)
else:
    plt.show()
