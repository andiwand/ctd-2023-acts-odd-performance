#!/usr/bin/env python3

import matplotlib.pyplot as plt
import argparse
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.stats import robust_mean, robust_std
from mycommon.data import get_pull_data


def plot_pulls_over_eta_sausage(input, fig):
    eta_range = (-3, 3)
    pull_range = (-4, 4)
    eta_bins = 30
    pull_bins = 30

    pull_labels = [
        r"$d_0$",
        r"$z_0$",
        r"$t$",
        r"$\phi$",
        r"$\theta$",
        r"$\frac{q}{p}$",
    ]

    eta, pulls = get_pull_data(input)

    subfigs = fig.subfigures(2, 3)

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
            statistic=robust_mean,
        )
        std_binned, _, _ = binned_statistic(
            eta,
            pull,
            bins=eta_bins,
            range=eta_range,
            statistic=robust_std,
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--output")
    args = parser.parse_args()

    fig = myPlotStyle()
    fig.set_size_inches(16, 10, forward=True)
    plot_pulls_over_eta_sausage(args.input, fig)

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
