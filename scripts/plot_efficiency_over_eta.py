#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
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


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
parser.add_argument("--particles", required=True, nargs="+")
parser.add_argument("--hits", required=True, nargs="+")
parser.add_argument("--output")
parser.add_argument("--require-pt", type=float, default=0.99)
parser.add_argument("--require-number-of-hits", type=int, default=7)
parser.add_argument("--matching-ratio", type=float, default=0.5)
args = parser.parse_args()

eta_range = (-3, 3)
eta_bins = 25

for tracksummary_file, particles_file, hits_file in zip(
    args.tracksummary, args.particles, args.hits
):
    event_label = get_event_label_from_path(tracksummary_file)
    event, _ = split_event_label(event_label)

    particles = ak.to_dataframe(
        uproot.open(particles_file)["particles"].arrays(
            ["event_id", "particle_id", "q", "eta", "pt", "vertex_primary"],
            library="ak",
        ),
        how="outer",
    ).dropna()

    # apply truth cuts
    # we only care about the first primary vertex which is the hard scatter for ttbar
    particles = particles[particles["vertex_primary"] == 1]
    # apply pt cut
    particles = particles[particles["pt"] >= args.require_pt]
    # filter for charged particles
    particles = particles[particles["q"] != 0.0]

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
    particles_hits = particles_hits.groupby(["event_id", "particle_id"]).aggregate(
        hit_count=pd.NamedAgg(column="volume_id", aggfunc="count"),
        eta=pd.NamedAgg(column="eta", aggfunc="first"),
        pt=pd.NamedAgg(column="pt", aggfunc="first"),
        vertex_primary=pd.NamedAgg(column="vertex_primary", aggfunc="first"),
    )
    particles_hits.reset_index(inplace=True)

    # calculate true efficiency
    particle_efficiency = particles_hits.copy()
    particle_efficiency["efficiency"] = (
        particle_efficiency["hit_count"].values >= args.require_number_of_hits
    ).astype(float)
    # drop true inefficiencies
    particle_efficiency = particle_efficiency[particle_efficiency["efficiency"] != 0.0]

    tracksummary = ak.to_dataframe(
        uproot.open(tracksummary_file)["tracksummary"].arrays(
            [
                "event_nr",
                "nMeasurements",
                "majorityParticleId",
                "nMajorityHits",
            ],
            library="ak",
        ),
        how="outer",
    ).dropna()

    # apply reco cuts
    # filter for track we care about
    tracksummary = tracksummary[
        (tracksummary["nMeasurements"].values >= args.require_number_of_hits)
    ]

    track_efficiency = pd.merge(
        tracksummary.add_prefix("track_"),
        particle_efficiency.add_prefix("true_"),
        how="right",
        left_on=["track_event_nr", "track_majorityParticleId"],
        right_on=["true_event_id", "true_particle_id"],
    )
    track_efficiency["track_nMeasurements"].fillna(0, inplace=True)

    track_efficiency["track_duplicate"] = (
        track_efficiency[["true_event_id", "true_particle_id"]]
        .duplicated(keep="first")
        .astype(float)
    )
    track_efficiency["track_efficiency"] = (
        (
            track_efficiency["track_nMajorityHits"].values
            / track_efficiency["track_nMeasurements"].values
            >= args.matching_ratio
        )
        & (track_efficiency["track_duplicate"].values == 0.0)
    ).astype(float)

    track_efficiency = track_efficiency.groupby(
        ["true_event_id", "true_particle_id"]
    ).aggregate(
        {
            "true_eta": "first",
            "true_pt": "first",
            "true_hit_count": "first",
            "track_efficiency": "max",
            "track_duplicate": "sum",
        }
    )
    track_efficiency.reset_index(inplace=True)

    track_efficiency_mean, eta_edges, _ = binned_statistic(
        track_efficiency["true_eta"],
        track_efficiency["track_efficiency"],
        bins=eta_bins,
        range=eta_range,
        statistic="mean",
    )
    track_efficiency_upper, _, _ = binned_statistic(
        track_efficiency["true_eta"],
        track_efficiency["track_efficiency"],
        bins=eta_bins,
        range=eta_range,
        statistic=create_clopper_pearson_upper_bounds(),
    )
    track_efficiency_lower, _, _ = binned_statistic(
        track_efficiency["true_eta"],
        track_efficiency["track_efficiency"],
        bins=eta_bins,
        range=eta_range,
        statistic=create_clopper_pearson_lower_bounds(),
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])
    eta_step = eta_edges[1] - eta_edges[0]

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
        label=get_event_variant_label(event),
    )

plt.axhline(1, linestyle="--", color="gray")

plt.title(f"Technical efficiency over $\eta$ for {get_event_type_label(event)} events")
plt.xlabel("$\eta$")
plt.ylabel("technical efficiency")
plt.xticks(np.linspace(*eta_range, 7))
plt.xlim(eta_range)
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
