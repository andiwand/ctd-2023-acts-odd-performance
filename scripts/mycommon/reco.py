from pathlib import Path
from typing import Optional, Union
from collections import namedtuple

import acts
from acts.examples.reconstruction import (
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedingAlgorithm,
    SeedFinderConfigArg,
    addSeeding,
    TrackSelectorConfig,
    TrackFindingConfig,
    AmbiguityResolutionConfig,
)

u = acts.UnitConstants

from mycommon.events import get_event_type


seedings = [
    "truth_smeared",
    "truth_estimated",
    # "default",
]

RecoConfig = namedtuple(
    "RecoCuts",
    ["track_selector_config", "track_finding_config", "ambi_config"],
    defaults=[TrackSelectorConfig(), TrackFindingConfig(), AmbiguityResolutionConfig()],
)


def list_seedings():
    return list(seedings)


def list_reco_labels():
    return [create_reco_label(seeding) for seeding in seedings]


def create_reco_label(seeding):
    return f"{seeding}"


def split_reco_label(reco_label):
    for seeding in list_seedings():
        if reco_label == create_reco_label(seeding):
            return seeding
    raise ValueError(f"unknown reco label {reco_label}")


def get_reco_config(event, seeding) -> RecoConfig:
    event_type = get_event_type(event)

    if event_type == "single_particles":
        return RecoConfig(
            track_selector_config=TrackSelectorConfig(),
            track_finding_config=TrackFindingConfig(
                chi2CutOff=15.0,
                numMeasurementsCutOff=10,
            ),
            ambi_config=AmbiguityResolutionConfig(
                maximumSharedHits=3,
                maximumIterations=100000000,
                nMeasurementsMin=3,
            ),
        )

    if event_type == "ttbar":
        return RecoConfig(
            track_selector_config=TrackSelectorConfig(
                pt=(0.5 * u.GeV, None),
                absEta=(None, 3.5),
                loc0=(-4.0 * u.mm, 4.0 * u.mm),
                nMeasurementsMin=7,
            ),
            track_finding_config=TrackFindingConfig(
                chi2CutOff=15.0,
                numMeasurementsCutOff=10,
            ),
            ambi_config=AmbiguityResolutionConfig(
                maximumSharedHits=3,
                maximumIterations=100000000,
                nMeasurementsMin=7,
            ),
        )

    raise ValueError(f"unknown event type {event_type}")


def addMySeeding(
    s: acts.examples.Sequencer,
    algortihm: str,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
    geoSelectionConfigFile: str,
    outputDirRoot: Optional[Union[Path, str]] = None,
):
    if algortihm == "truth_smeared":
        seedingAlgorithm = SeedingAlgorithm.TruthSmeared
    elif algortihm == "truth_estimated":
        seedingAlgorithm = SeedingAlgorithm.TruthEstimated
    elif algortihm == "default":
        seedingAlgorithm = SeedingAlgorithm.Default
    else:
        raise ValueError(f"unknown seeding algorithm: {algortihm}")

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        seedingAlgorithm=seedingAlgorithm,
        truthSeedRanges=TruthSeedRanges(
            # using something close to 1 to include for sure
            pt=(0.999 * u.GeV, None),
            eta=(-3.0, 3.0),
            nHits=(3, None),
        ),
        particleSmearingSigmas=ParticleSmearingSigmas(
            d0=20 * u.um,
            d0PtA=30 * u.um,
            d0PtB=0.3 / u.GeV,
            z0=20 * u.um,
            z0PtA=30 * u.um,
            z0PtB=0.3 / u.GeV,
            t0=1 * u.ns,
            phi=0.1 * u.degree,
            theta=0.1 * u.degree,
            pRel=0.01,
        ),
        seedFinderConfigArg=SeedFinderConfigArg(
            r=(33 * u.mm, 200 * u.mm),
            deltaR=(1 * u.mm, 60 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=5,
            radLengthPerSeed=0.1,
            minPt=0.5 * u.GeV,
            impactMax=3 * u.mm,
        ),
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0.1 / u.GeV,
            1 * u.ns,
        ],
        initialVarInflation=[100.0] * 6,
        geoSelectionConfigFile=geoSelectionConfigFile,
        outputDirRoot=outputDirRoot,
    )
