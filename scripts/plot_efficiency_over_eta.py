#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import numpy as np
import pandas as pd
import argparse
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.label import pt_label, event_type_label


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
parser.add_argument("--particles", required=True, nargs="+")
parser.add_argument("--hits", required=True, nargs="+")
parser.add_argument("--output")
parser.add_argument("--require-number-of-hits", type=int, default=7)
parser.add_argument("--matching-ratio", type=float, default=0.5)
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 25

for tracksummary_file, particles_file, hits_file in zip(
    args.tracksummary, args.particles, args.hits
):
    particles = ak.to_dataframe(
        uproot.open(particles_file)["particles"].arrays(
            ["event_id", "particle_id", "eta", "pt"],
            library="ak",
        ),
        how="outer",
    ).dropna()

    hits = ak.to_dataframe(
        uproot.open(hits_file)["hits"].arrays(
            [
                "event_id",
                "particle_id",
                "volume_id",
                "layer_id",
                "sensitive_id",
                "index",
            ],
            library="ak",
        ),
        how="outer",
    ).dropna()

    particles_hits = pd.merge(
        particles,
        hits,
        how="left",
        left_on=["event_id", "particle_id"],
        right_on=["event_id", "particle_id"],
    )
    particle_efficiency = particles_hits.groupby(["event_id", "particle_id"]).aggregate(
        hit_count=pd.NamedAgg(column="volume_id", aggfunc="count"),
        eta=pd.NamedAgg(column="eta", aggfunc="first"),
    )
    particle_efficiency["efficiency"] = (
        particle_efficiency["hit_count"].values >= args.require_number_of_hits
    ).astype(float)

    tracksummary = ak.to_dataframe(
        uproot.open(tracksummary_file)["tracksummary"].arrays(
            [
                "event_nr",
                "t_eta",
                "t_pT",
                "nMeasurements",
                "majorityParticleId",
                "nMajorityHits",
            ],
            library="ak",
        ),
        how="outer",
    ).dropna()

    track_efficiency = tracksummary.copy()
    track_efficiency["duplicate"] = (
        track_efficiency[["event_nr", "majorityParticleId"]]
        .duplicated(keep="first")
        .astype(float)
    )
    track_efficiency["efficiency"] = (
        (track_efficiency["nMeasurements"].values >= args.require_number_of_hits)
        & (
            track_efficiency["nMajorityHits"].values
            / track_efficiency["nMeasurements"].values
            >= args.matching_ratio
        )
        & (track_efficiency["duplicate"].values == 0)
    ).astype(float)

    particle_efficiency_binned, eta_edges, _ = binned_statistic(
        particle_efficiency["eta"],
        particle_efficiency["efficiency"],
        bins=eta_bins,
        range=eta_range,
        statistic="mean",
    )
    track_efficiency_binned, _, _ = binned_statistic(
        track_efficiency["t_eta"],
        track_efficiency["efficiency"],
        bins=eta_bins,
        range=eta_range,
        statistic="mean",
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    plt.plot(
        eta_mid,
        track_efficiency_binned / particle_efficiency_binned,
        marker="o",
        label=pt_label(tracksummary_file),
    )

plt.axhline(1, linestyle="--", color="gray")

plt.title(f"Efficiency over $\eta$ for {event_type_label(args.tracksummary[0])} events")
plt.xlabel("$\eta$")
plt.ylabel("efficiency")
plt.xticks(np.linspace(*eta_range, 7))
plt.xlim(eta_range)
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
