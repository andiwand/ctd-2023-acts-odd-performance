from pathlib import Path
from typing import Optional, Union

import acts
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedingAlgorithm,
)

u = acts.UnitConstants

seedings = [
    "truth_smeared",
]


def list_seedings():
    return list(seedings)


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
    else:
        raise ValueError(f"unknown seeding algorithm: {algortihm}")

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        seedingAlgorithm=seedingAlgorithm,
        truthSeedRanges=TruthSeedRanges(
            nHits=(7, None),
        ),
        particleSmearingSigmas=ParticleSmearingSigmas(
            d0=20 * u.um,
            d0PtA=30 * u.um,
            d0PtB=0.3 / u.GeV,
            z0=20 * u.um,
            z0PtA=30 * u.um,
            z0PtB=0.3 / u.GeV,
            t0=1 * u.ns,
            phi=1 * u.degree,
            theta=1 * u.degree,
            pRel=0.05,
        ),
        initialSigmas=[0.1, 0.1, 0.002, 0.0001, 0.001, 1000],
        initialVarInflation=[1e3] * 6,
        geoSelectionConfigFile=geoSelectionConfigFile,
        outputDirRoot=outputDirRoot,
    )
