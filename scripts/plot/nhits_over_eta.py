#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import pandas as pd
import argparse
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle


def plot_nhits_over_eta(particles, hits, fig, ax):
    eta_range = (-3, 3)
    hits_range = (8, 20)
    eta_bins = 25
    hits_bins = hits_range[1] - hits_range[0]

    particles = ak.to_dataframe(
        uproot.open(particles)["particles"].arrays(
            ["event_id", "particle_id", "eta"], library="ak"
        ),
        how="outer",
    ).dropna()

    hits = ak.to_dataframe(
        uproot.open(hits)["hits"].arrays(["event_id", "particle_id"], library="ak"),
        how="outer",
    ).dropna()

    particle_hits = pd.merge(
        particles, hits, on=["event_id", "particle_id"], how="left"
    )
    particle_hits["hits"] = 1
    particle_hits = particle_hits.groupby(["event_id", "particle_id"]).agg(
        {
            "eta": "first",
            "hits": "sum",
        }
    )
    particle_hits.reset_index(drop=True, inplace=True)

    eta = particle_hits["eta"].values
    hits = particle_hits["hits"].values

    mean, eta_edges, _ = binned_statistic(
        eta, hits, bins=eta_bins, range=eta_range, statistic="mean"
    )
    std, _, _ = binned_statistic(
        eta, hits, bins=eta_bins, range=eta_range, statistic="std"
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])
    eta_step = eta_mid[1] - eta_mid[0]

    ax.errorbar(
        x=eta_mid,
        y=mean,
        yerr=std,
        xerr=eta_step * 0.4,
        fmt="",
        linestyle="",
        color="black",
    )

    _, _, _, im = ax.hist2d(
        x=eta,
        y=hits,
        bins=(eta_bins, hits_bins),
        range=(eta_range, hits_range),
        cmap="Blues",
        density=True,
    )

    fig.colorbar(im, ax=fig.get_axes())
    ax.set_title("Number of hits over $\eta$")
    ax.set_xlabel("$\eta$")
    ax.set_ylabel("number of hits")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("particles")
    parser.add_argument("hits")
    parser.add_argument("--output")
    args = parser.parse_args()
    args = parser.parse_args()

    fig = myPlotStyle()
    plot_nhits_over_eta(args.particles, args.hits, fig, fig.gca())

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
