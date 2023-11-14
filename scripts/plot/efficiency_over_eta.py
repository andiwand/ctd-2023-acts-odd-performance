#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import argparse

from mycommon.plot_style import myPlotStyle
from mycommon.events import split_event_label
from mycommon.label import (
    get_event_variant_label,
    get_event_label,
    get_event_type_label,
)
from mycommon.paths import get_event_label_from_path, check_same_event_type
from mycommon.io import read_track_efficiency
from mycommon.agg import agg_efficiency_over_eta


def plot_efficiency_over_eta(input, fig, ax):
    eta_range = (0, 3)
    eta_bins = 13

    is_same_event_type = check_same_event_type(input)

    for file in input:
        event_label = get_event_label_from_path(file)
        event, _ = split_event_label(event_label)

        eta, track_efficiency = read_track_efficiency(file)

        (
            eta_mid,
            (
                track_efficiency_mean,
                (
                    lower_error,
                    upper_error,
                ),
            ),
        ) = agg_efficiency_over_eta(eta_range, eta_bins, eta, track_efficiency)
        eta_step = eta_mid[1] - eta_mid[0]

        label = (
            get_event_variant_label(event)
            if is_same_event_type
            else get_event_label(event)
        )
        ax.errorbar(
            x=eta_mid,
            y=track_efficiency_mean,
            yerr=(lower_error, upper_error),
            xerr=eta_step * 0.4,
            fmt="",
            linestyle="",
            label=label,
        )

    ax.axhline(1, linestyle="--", color="gray")

    if is_same_event_type:
        ax.set_title(
            f"Technical efficiency over $\eta$ for {get_event_type_label(event)} events"
        )
    else:
        ax.set_title(f"Technical efficiency over $\eta$")
    ax.set_xlabel("$\eta$")
    ax.set_ylabel("technical efficiency")
    ax.set_xticks(np.linspace(*eta_range, 7))
    ax.set_xlim(eta_range)
    ax.legend()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+")
    parser.add_argument("--output")
    args = parser.parse_args()

    fig = myPlotStyle()
    plot_efficiency_over_eta(args.input, fig, fig.gca())

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
