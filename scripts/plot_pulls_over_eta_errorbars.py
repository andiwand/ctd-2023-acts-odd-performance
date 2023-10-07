#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
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
from mycommon.data import get_pull_data


def plot_pulls_over_eta_errorbars(input, fig):
    eta_range = (0, 3)
    pull_range = (-2, 2)
    eta_bins = 8

    pull_labels = [
        r"$d_0$",
        r"$z_0$",
        # r"$t$",
        # r"$\phi$",
        # r"$\theta$",
        r"$\frac{q}{p}$",
    ]

    axs = fig.subplots(3, 1, sharex=True, sharey=True)

    for i, file in enumerate(input):
        event_label = get_event_label_from_path(file)
        event, _ = split_event_label(event_label)

        eta, pulls = get_pull_data(file)
        pulls = pulls[0], pulls[1], pulls[5]

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+")
    parser.add_argument("--output")
    args = parser.parse_args()

    fig = myPlotStyle()
    plot_pulls_over_eta_errorbars(args.input, fig)

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
