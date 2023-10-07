#!/usr/bin/env python3

import uproot
import awkward as ak
import numpy as np
import pandas as pd
import argparse

import acts

u = acts.UnitConstants


pixel_volumes = [16, 17, 18]
pixel_volumes_first_layers = [10, 2, 2]


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
parser.add_argument("--require-number-of-pixel-hits", type=int, default=3)
parser.add_argument("--require-number-of-layer1-pixel-hits", type=int, default=1)
parser.add_argument("--matching-ratio", type=float, default=0.5)
args = parser.parse_args()

print(f"read particles...")
particles = ak.to_dataframe(
    uproot.open(args.particles)["particles"].arrays(
        [
            "event_id",
            "particle_id",
            "q",
            "phi",
            "eta",
            "p",
            "pt",
            "vertex_primary",
            "vx",
            "vy",
            "vz",
        ],
        library="ak",
    ),
    how="outer",
).dropna()
print(f"{len(particles)} particles read.")

print(f"apply cuts to particles...")
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
    np.hypot(particles["vx"], particles["vy"]) < args.require_max_r * u.mm
]
particles = particles.copy()
print(f"{len(particles)} particles remaining.")

print(f"read hits...")
hits = ak.to_dataframe(
    uproot.open(args.hits)["hits"].arrays(
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
print(f"{len(hits)} hits read.")

print(f"reduce hits to known particles...")
hits = pd.merge(
    particles[["event_id", "particle_id"]],
    hits,
    how="left",
    left_on=["event_id", "particle_id"],
    right_on=["event_id", "particle_id"],
)
hits.reset_index(drop=True, inplace=True)
print(f"{len(hits)} hits left.")

print(f"group hits...")
hits["hits"] = 1
hits["hits_pixel"] = hits["volume_id"].isin(pixel_volumes)
hits["hits_pixel_layer1"] = sum(
    ((hits["volume_id"] == volume_id) & (hits["layer_id"] == layer_id))
    for volume_id, layer_id in zip(pixel_volumes, pixel_volumes_first_layers)
)
hits = hits.groupby(["event_id", "particle_id"]).agg(
    {
        "event_id": "first",
        "particle_id": "first",
        "hits": "sum",
        "hits_pixel": "sum",
        "hits_pixel_layer1": "sum",
    }
)
hits.reset_index(drop=True, inplace=True)

print(f"merge particles and hits...")
particle_efficiency = pd.merge(
    particles,
    hits,
    how="left",
    left_on=["event_id", "particle_id"],
    right_on=["event_id", "particle_id"],
)
particle_efficiency.reset_index(inplace=True)
particles = None
hits = None

print(f"calculate true efficiency and cut...")
# calculate true efficiency
particle_efficiency["efficiency"] = (
    (particle_efficiency["hits"].values >= args.require_number_of_hits)
    & (particle_efficiency["hits_pixel"].values >= args.require_number_of_pixel_hits)
    & (
        particle_efficiency["hits_pixel_layer1"].values
        >= args.require_number_of_layer1_pixel_hits
    )
).astype(int)
# drop true inefficiencies
particle_efficiency = particle_efficiency[particle_efficiency["efficiency"] != 0]
print(f"{len(particle_efficiency)} particles remaining.")

print(f"read tracksummary...")
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
)
# for whatever reason multiTraj_nr and subTraj_nr are not filled properly, creating random index
tracksummary.reset_index(drop=True, inplace=True)
tracksummary["track_nr"] = tracksummary.index
tracksummary = tracksummary.dropna()

print(f"aggregate tracks...")
track_efficiency = pd.merge(
    particle_efficiency[["event_id", "particle_id"]].add_prefix("true_"),
    tracksummary[
        [
            "event_nr",
            "track_nr",
            "nMeasurements",
            "chi2Sum",
            "majorityParticleId",
            "nMajorityHits",
        ]
    ].add_prefix("track_"),
    how="left",
    left_on=["true_event_id", "true_particle_id"],
    right_on=["track_event_nr", "track_majorityParticleId"],
)
track_efficiency["track_nMeasurements"].fillna(0, inplace=True)
track_efficiency.sort_values(
    ["track_nMeasurements", "track_chi2Sum"], ascending=[False, True], inplace=True
)
track_efficiency["track_duplicate"] = 1
aggregate_dict = {column: "first" for column in track_efficiency.columns}
aggregate_dict["track_duplicate"] = "sum"
track_efficiency = track_efficiency.groupby(["true_event_id", "true_particle_id"]).agg(
    aggregate_dict
)
track_efficiency["track_duplicate"] -= 1
track_efficiency.reset_index(drop=True, inplace=True)

print(f"calculate track efficiency...")
track_efficiency["track_efficiency"] = (
    (track_efficiency["track_nMeasurements"].values >= args.require_number_of_hits)
    & (
        track_efficiency["track_nMajorityHits"].values
        / track_efficiency["track_nMeasurements"].values
        >= args.matching_ratio
    )
).astype(int)

print(f"merge tracks into track efficiency...")
track_efficiency = pd.merge(
    track_efficiency[
        [
            "true_event_id",
            "true_particle_id",
            "track_event_nr",
            "track_track_nr",
            "track_duplicate",
            "track_efficiency",
        ]
    ],
    tracksummary.add_prefix("track_"),
    how="left",
    on=["track_event_nr", "track_track_nr"],
)
print(f"merge particles into track efficiency...")
track_efficiency = pd.merge(
    track_efficiency,
    particle_efficiency.add_prefix("true_"),
    how="left",
    on=["true_event_id", "true_particle_id"],
)
track_efficiency.reset_index(inplace=True)

print(f"write {len(track_efficiency)} tracks to csv...")
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
        "true_hits_pixel",
        "true_hits_pixel_layer1",
        "track_track_nr",
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
