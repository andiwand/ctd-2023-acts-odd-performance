#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.events import split_event_label
from mycommon.label import get_event_variant_label, get_event_type_label
from mycommon.paths import get_event_label_from_path
from mycommon.stats import (
    create_clopper_pearson_upper_bounds,
    create_clopper_pearson_lower_bounds,
)


def check_same_event_type(input):
    def get_event_from_path(file):
        event_label = get_event_label_from_path(file)
        event, _ = split_event_label(event_label)
        return event

    event_types = [get_event_from_path(file) for file in input]
    return len(set(event_types)) == 1


def get_data(file):
    if str(file).endswith(".csv"):
        data = pd.read_csv(file)

        eta = data["true_eta"].values
        track_efficiency = data["track_efficiency"].values

        return eta, track_efficiency

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 25

is_same_event_type = check_same_event_type(args.input)

for file in args.input:
    event_label = get_event_label_from_path(file)
    event, _ = split_event_label(event_label)

    eta, track_efficiency = get_data(file)

    track_efficiency_mean, eta_edges, _ = binned_statistic(
        eta,
        track_efficiency,
        bins=eta_bins,
        range=eta_range,
        statistic="mean",
    )
    track_efficiency_upper, _, _ = binned_statistic(
        eta,
        track_efficiency,
        bins=eta_bins,
        range=eta_range,
        statistic=create_clopper_pearson_upper_bounds(),
    )
    track_efficiency_lower, _, _ = binned_statistic(
        eta,
        track_efficiency,
        bins=eta_bins,
        range=eta_range,
        statistic=create_clopper_pearson_lower_bounds(),
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])
    eta_step = eta_edges[1] - eta_edges[0]

    label = (
        get_event_variant_label(event)
        if is_same_event_type
        else get_event_type_label(event)
    )
    plt.errorbar(
        x=eta_mid,
        y=track_efficiency_mean,
        yerr=(
            track_efficiency_mean - track_efficiency_lower,
            track_efficiency_upper - track_efficiency_mean,
        ),
        xerr=eta_step * 0.4,
        fmt="",
        linestyle="",
        label=label,
    )

plt.axhline(1, linestyle="--", color="gray")

if is_same_event_type:
    plt.title(
        f"Technical efficiency over $\eta$ for {get_event_variant_label(event)} events"
    )
else:
    plt.title(f"Technical efficiency over $\eta$")
plt.xlabel("$\eta$")
plt.ylabel("technical efficiency")
plt.xticks(np.linspace(*eta_range, 7))
plt.xlim(eta_range)
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
