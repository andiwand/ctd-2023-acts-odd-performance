#!/usr/bin/env python3

import uproot
import awkward as ak
import numpy as np
import pandas as pd
import argparse

import acts

u = acts.UnitConstants


def aggregate_tracks(group):
    if len(group) == 1:
        return group

    best_idx = group[
        group["track_nMeasurements"] == group["track_nMeasurements"].max()
    ]["track_chi2Sum"].idxmin()
    best = group[group.index == best_idx].copy()

    best["track_duplicate"] = group["track_duplicate"].sum()
    best["track_efficiency"] = group["track_efficiency"].max()

    return best


parser = argparse.ArgumentParser()
parser.add_argument("tracksummary")
parser.add_argument("particles")
parser.add_argument("hits")
parser.add_argument("output")
parser.add_argument("--require-primary-vertex", type=int, default=1)
parser.add_argument("--require-pt", type=float, default=1, help="in GeV")
parser.add_argument("--require-max-absz", type=float, default=1.0, help="in meters")
parser.add_argument("--require-max-r", type=float, default=24.0, help="in mm")
parser.add_argument("--require-number-of-hits", type=int, default=7)
parser.add_argument("--matching-ratio", type=float, default=0.5)
args = parser.parse_args()

particles = ak.to_dataframe(
    uproot.open(args.particles)["particles"].arrays(library="ak"),
    how="outer",
).dropna()

# apply truth cuts
# we only care about the first primary vertex which is the hard scatter for ttbar
particles = particles[particles["vertex_primary"] == args.require_primary_vertex]
# apply pt cut
particles = particles[particles["pt"] >= args.require_pt * u.GeV]
# filter for charged particles
particles = particles[particles["q"] != 0.0]
# filter particles originating from the beampipe
particles = particles[particles["vz"].abs() < args.require_max_absz * u.m]
particles = particles[
    np.hypot(particles["vx"], particles["vx"]) < args.require_max_r * u.mm
]

hits = ak.to_dataframe(
    uproot.open(args.hits)["hits"].arrays(library="ak"), how="outer"
).dropna()

particles_hits = pd.merge(
    particles,
    hits,
    how="left",
    left_on=["event_id", "particle_id"],
    right_on=["event_id", "particle_id"],
)
del particles
del hits

particle_efficiency = particles_hits.groupby(["event_id", "particle_id"]).aggregate(
    hits=pd.NamedAgg(column="volume_id", aggfunc="count"),
    q=pd.NamedAgg(column="q", aggfunc="first"),
    phi=pd.NamedAgg(column="phi", aggfunc="first"),
    eta=pd.NamedAgg(column="eta", aggfunc="first"),
    p=pd.NamedAgg(column="p", aggfunc="first"),
    pt=pd.NamedAgg(column="pt", aggfunc="first"),
    vertex_primary=pd.NamedAgg(column="vertex_primary", aggfunc="first"),
)
del particles_hits
particle_efficiency.reset_index(inplace=True)

# calculate true efficiency
particle_efficiency["efficiency"] = (
    particle_efficiency["hits"].values >= args.require_number_of_hits
).astype(int)
# drop true inefficiencies
particle_efficiency = particle_efficiency[particle_efficiency["efficiency"] != 0.0]

tracksummary = ak.to_dataframe(
    uproot.open(args.tracksummary)["tracksummary"].arrays(
        [
            "event_nr",
            "nStates",
            "nMeasurements",
            "nOutliers",
            "nHoles",
            "nSharedHits",
            "chi2Sum",
            "majorityParticleId",
            "nMajorityHits",
            "eLOC0_fit",
            "eLOC1_fit",
            "eT_fit",
            "ePHI_fit",
            "eTHETA_fit",
            "eQOP_fit",
            "res_eLOC0_fit",
            "res_eLOC1_fit",
            "res_eT_fit",
            "res_ePHI_fit",
            "res_eTHETA_fit",
            "res_eQOP_fit",
            "pull_eLOC0_fit",
            "pull_eLOC1_fit",
            "pull_eT_fit",
            "pull_ePHI_fit",
            "pull_eTHETA_fit",
            "pull_eQOP_fit",
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
    particle_efficiency.add_prefix("true_"),
    tracksummary.add_prefix("track_"),
    how="left",
    left_on=["true_event_id", "true_particle_id"],
    right_on=["track_event_nr", "track_majorityParticleId"],
)
del particle_efficiency
del tracksummary
track_efficiency["track_nMeasurements"].fillna(0, inplace=True)

track_efficiency["track_duplicate"] = (
    track_efficiency[["true_event_id", "true_particle_id"]]
    .duplicated(keep="first")
    .astype(int)
)
track_efficiency["track_efficiency"] = (
    (
        track_efficiency["track_nMajorityHits"].values
        / track_efficiency["track_nMeasurements"].values
        >= args.matching_ratio
    )
    & (track_efficiency["track_duplicate"].values == 0)
).astype(int)

track_efficiency = track_efficiency.groupby(
    ["true_event_id", "true_particle_id"]
).apply(aggregate_tracks)
track_efficiency.reset_index(drop=True, inplace=True)

track_efficiency[
    [
        "true_event_id",
        "true_particle_id",
        "true_q",
        "true_phi",
        "true_eta",
        "true_p",
        "true_pt",
        "true_vertex_primary",
        "true_hits",
        "track_duplicate",
        "track_efficiency",
        "track_nStates",
        "track_nMeasurements",
        "track_nOutliers",
        "track_nHoles",
        "track_nSharedHits",
        "track_chi2Sum",
        "track_nMajorityHits",
        "track_eLOC0_fit",
        "track_eLOC1_fit",
        "track_eT_fit",
        "track_ePHI_fit",
        "track_eTHETA_fit",
        "track_eQOP_fit",
        "track_res_eLOC0_fit",
        "track_res_eLOC1_fit",
        "track_res_eT_fit",
        "track_res_ePHI_fit",
        "track_res_eTHETA_fit",
        "track_res_eQOP_fit",
        "track_pull_eLOC0_fit",
        "track_pull_eLOC1_fit",
        "track_pull_eT_fit",
        "track_pull_ePHI_fit",
        "track_pull_eTHETA_fit",
        "track_pull_eQOP_fit",
    ]
].to_csv(args.output, index=False)
