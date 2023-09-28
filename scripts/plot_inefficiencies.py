#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

from mycommon.plot_style import myPlotStyle
from mycommon.events import split_event_label
from mycommon.label import get_event_label
from mycommon.paths import get_event_label_from_path


def get_data(file):
    if str(file).endswith(".csv"):
        data = pd.read_csv(file)

        phi = data["true_phi"].values
        eta = data["true_eta"].values
        true_hits = data["true_hits"].values
        track_states = data["track_nStates"].values
        track_hits = data["track_nMeasurements"].values
        track_efficiency = data["track_efficiency"].values

        return phi, eta, track_states, true_hits, track_hits, track_efficiency

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("--output")
args = parser.parse_args()

phi_range = (0, 2 * np.pi)
eta_range = (-3, 3)
hit_range = (0, 25)
phi_bins = 25
eta_bins = 25
hit_bins = 25

event_label = get_event_label_from_path(args.input)
event, _ = split_event_label(event_label)

phi, eta, true_hits, track_states, track_hits, track_efficiency = get_data(args.input)

axes = plt.gcf().subplots(3, 2, sharey=True)

_, _, _, im = axes[0, 0].hist2d(
    phi,
    true_hits,
    weights=1 - track_efficiency,
    bins=(phi_range, hit_bins),
    range=(phi_range, hit_range),
    cmap="Reds",
)
axes[0, 0].set_xlabel(r"$\phi$")
axes[0, 0].set_ylabel(r"# true hits")
axes[0, 0].set_xticks(np.linspace(*phi_range, 7))
axes[0, 0].set_xlim(phi_range)
plt.colorbar(im, ax=axes[0, 0])

axes[1, 0].hist2d(
    phi,
    track_hits,
    weights=1 - track_efficiency,
    bins=(phi_range, hit_bins),
    range=(phi_range, hit_range),
    cmap="Reds",
)
axes[1, 0].set_xlabel(r"$\phi$")
axes[1, 0].set_ylabel(r"# measurements")
axes[1, 0].set_xticks(np.linspace(*phi_range, 7))
axes[1, 0].set_xlim(phi_range)
plt.colorbar(im, ax=axes[1, 0])

axes[2, 0].hist2d(
    phi,
    track_states,
    weights=1 - track_efficiency,
    bins=(phi_range, hit_bins),
    range=(phi_range, hit_range),
    cmap="Reds",
)
axes[2, 0].set_xlabel(r"$\phi$")
axes[2, 0].set_ylabel(r"# tack states")
axes[2, 0].set_xticks(np.linspace(*phi_range, 7))
axes[2, 0].set_xlim(phi_range)
plt.colorbar(im, ax=axes[2, 0])

axes[0, 1].hist2d(
    eta,
    true_hits,
    weights=1 - track_efficiency,
    bins=(eta_bins, hit_bins),
    range=(eta_range, hit_range),
    cmap="Reds",
)
axes[0, 1].set_xlabel(r"$\eta$")
axes[0, 1].set_ylabel(r"# true hits")
axes[0, 1].set_xticks(np.linspace(*eta_range, 7))
axes[0, 1].set_xlim(eta_range)
plt.colorbar(im, ax=axes[0, 1])

axes[1, 1].hist2d(
    eta,
    track_hits,
    weights=1 - track_efficiency,
    bins=(eta_bins, hit_bins),
    range=(eta_range, hit_range),
    cmap="Reds",
)
axes[1, 1].set_xlabel(r"$\eta$")
axes[1, 1].set_ylabel(r"# measurements")
axes[1, 1].set_xticks(np.linspace(*eta_range, 7))
axes[1, 1].set_xlim(eta_range)
plt.colorbar(im, ax=axes[1, 1])

axes[2, 1].hist2d(
    eta,
    track_states,
    weights=1 - track_efficiency,
    bins=(eta_bins, hit_bins),
    range=(eta_range, hit_range),
    cmap="Reds",
)
axes[2, 1].set_xlabel(r"$\eta$")
axes[2, 1].set_ylabel(r"# track states")
axes[2, 1].set_xticks(np.linspace(*eta_range, 7))
axes[2, 1].set_xlim(eta_range)
plt.colorbar(im, ax=axes[2, 1])

plt.suptitle(f"Inefficiency over $\eta$ for {get_event_label(event)}")

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
