from pathlib import Path
from typing import Optional, Union

import acts
from acts.examples.reconstruction import (
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedingAlgorithm,
    SeedFinderConfigArg,
    addSeeding,
)

u = acts.UnitConstants

seedings = [
    "truth_smeared",
    "truth_estimated",
    # "default",
]


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
        initialVarInflation=[1e3] * 6,
        geoSelectionConfigFile=geoSelectionConfigFile,
        outputDirRoot=outputDirRoot,
    )
