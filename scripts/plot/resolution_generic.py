#!/usr/bin/env python3

import matplotlib.pyplot as plt
import argparse
import numpy as np

from mycommon.plot_style import myPlotStyle
from mycommon.events import split_event_label
from mycommon.label import (
    get_event_variant_label,
    get_event_type_label,
    get_event_label,
)
from mycommon.paths import get_event_label_from_path, check_same_event_type
from mycommon.io import read_residuals
from mycommon.agg import agg_resolution


def plot_resolution(x, y, input, fig, ax):
    x_range = {
        "eta": (0, 3),
        "pt": (1, 100),
    }[x]
    x_bins = {
        "eta": 13,
        "pt": 10,
    }[x]
    x_label = {
        "eta": r"\eta",
        "pt": r"p_T",
    }[x]
    x_unit = {
        "eta": r"",
        "pt": r" [GeV]",
    }[x]

    y_label = {
        "d0": r"d_0",
        "z0": r"z_0",
        "qop": r"\frac{q}{p}",
    }[y]
    y_unit = {
        "d0": r" [mm]",
        "z0": r" [mm]",
        "qop": r" [$\frac{1}{GeV}$]",
    }[y]

    is_same_event_type = check_same_event_type(input)

    for file in input:
        event_label = get_event_label_from_path(file)
        event, _ = split_event_label(event_label)

        data = read_residuals(file)

        (x_mid, (std, std_std)) = agg_resolution(
            x_range, x_bins, data[x], data[f"res_{y}"]
        )
        x_step = x_mid[1] - x_mid[0]

        label = (
            get_event_variant_label(event)
            if is_same_event_type
            else get_event_label(event)
        )
        ax.errorbar(
            x=x_mid,
            y=std,
            yerr=std_std,
            xerr=x_step * 0.4,
            fmt="",
            linestyle="",
            label=label,
        )

    if is_same_event_type:
        ax.set_title(
            rf"Resolution of ${y_label}$ over ${x_label}$ for {get_event_type_label(event)} events"
        )
    else:
        ax.set_title(rf"Resolution of ${y_label}$ over ${x_label}$")
    ax.set_xlabel(rf"${x_label}${x_unit}")
    ax.set_ylabel(rf"$\sigma_{{{y_label}}}${y_unit}")
    if x == "eta":
        ax.set_xticks(np.linspace(*x_range, 7))
        ax.set_xlim(x_range)
    elif x == "pt":
        ax.set_xticks(np.linspace(0, x_range[1], 11))
        ax.set_xlim(x_range)
    else:
        raise ValueError(f"unknown x: {x}")
    ax.legend()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("x", choices=["eta", "pt"])
    parser.add_argument("y", choices=["d0", "z0", "qop"])
    parser.add_argument("input", nargs="+")
    parser.add_argument("--output")
    args = parser.parse_args()

    fig = myPlotStyle()
    plot_resolution(args.x, args.y, args.input, fig, fig.gca())

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
