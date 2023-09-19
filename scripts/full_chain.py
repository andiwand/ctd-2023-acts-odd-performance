#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
import itertools
from multiprocessing import Pool, freeze_support
import argparse
from typing import Any, Optional, Union
import fnmatch

import acts
from acts.examples.simulation import (
    addParticleGun,
    addPythia8,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedingAlgorithm,
    addKalmanTracks,
    addTrajectoryWriters,
    addCKFTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    VertexFinder,
    addVertexFitting,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants

# ODD configs
geoDir = getOpenDataDetectorDirectory()
materialMap = geoDir / "data/odd-material-maps.root"
digiConfig = geoDir / "config/odd-digi-smearing-config.json"
seedingSel = geoDir / "config/odd-seeding-config.json"
materialDeco = acts.IMaterialDecorator.fromFile(materialMap)

# ODD
detector, trackingGeometry, decorators = getOpenDataDetector(
    geoDir, mdecorator=materialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))


def main():
    # TODO what this
    freeze_support()

    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("--pool-size", default=4)
    parser.add_argument("--filter", default="*")
    args = parser.parse_args()

    events = [
        f"{p}_{pT}GeV" for p, pT in itertools.product(["mu", "pi", "e"], [1, 10, 100])
    ]
    simulations = [
        "fatras",
        "geant4",
    ]
    seedings = [
        "truth_smeared",
    ]

    # TODO fix multiple G4 setups in one process
    # maxtasksperchild=1 because we want to avoid multiple G4 setups in one process
    with Pool(args.pool_size, maxtasksperchild=1) as pool:
        for event, simulation, seeding in itertools.product(
            events, simulations, seedings
        ):
            label = create_label(event, simulation, seeding)

            if not fnmatch.fnmatch(label, args.filter):
                print(f"skipping {label}")
                continue

            outdir = Path(args.outdir) / label
            outdir.mkdir(parents=True, exist_ok=True)

            # TODO number of events will depend on what kind of event we are processing
            events = 100000
            skip = 0

            with tempfile.TemporaryDirectory() as temp:
                # for debugging
                #run_single_particles(outdir, temp, events, skip, event, simulation, seeding)

                pool.apply_async(
                    run_single_particles,
                    args=(outdir, temp, events, skip, event, simulation, seeding),
                )

        pool.close()
        pool.join()


def create_label(event, simulation, seeding):
    return f"{event}_{simulation}_{seeding}"


def run_single_particles(outdir, temp, events, skip, event, simulation, seeding):
    # single thread because we do multiprocessing
    numThreads = 1

    tp = Path(temp)
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=events,
        skip=skip,
        numThreads=numThreads,
        logLevel=acts.logging.INFO,
        trackFpes=False,
    )

    for d in decorators:
        s.addContextDecorator(d)

    eventType = addMyEventGen(
        s,
        event,
        rnd=rnd,
        outputDirRoot=tp,
    )

    addMySimulation(
        s,
        simulation,
        trackingGeometry,
        field,
        rnd=rnd,
        detector=detector,
        preSelectParticles=ParticleSelectorConfig(),
        postSelectParticles=ParticleSelectorConfig(removeSecondaries=True),
        outputDirRoot=tp,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfig,
        rnd=rnd,
        outputDirRoot=tp,
    )

    addMySeeding(
        s,
        seeding,
        trackingGeometry,
        field,
        rnd=rnd,
        outputDirRoot=tp,
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold=100 * u.TeV,
    )
    addTrajectoryWriters(
        s,
        name="kf",
        trajectories="kfTrajectories",
        outputDirRoot=tp,
        writeStates=True,
        writeSummary=True,
        writeCKFperformance=True,
        writeFinderPerformance=False,
        writeFitterPerformance=False,
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        outputDirRoot=tp,
    )

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3,
            maximumIterations=10000,
            nMeasurementsMin=7,
        ),
        outputDirRoot=tp,
    )

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=tp,
    )

    s.run()
    del s

    for stem in [
        "particles",
        "particles_initial",
        "particles_final",
        "hits",
        "measurements",
        "tracksummary_ckf",
        "trackstates_ckf",
        "performance_ckf",
        "tracksummary_ambi",
        "trackstates_ambi",
        "performance_ambi",
        "tracksummary_kf",
        "trackstates_kf",
        "performance_kf",
        "performance_ivf",
    ]:
        perf_file = tp / f"{stem}.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(perf_file, outdir / f"{stem}.root")


def addMyEventGen(
    s: acts.examples.Sequencer,
    event: str,
    rnd: acts.examples.RandomNumbers,
    outputDirRoot: Optional[Union[Path, str]] = None,
):
    vtxGen = acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 1 * u.ns),
    )

    # generate special events on top

    if event.startswith("ttbar_"):
        pu = float(event.split("_")[1])

        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=pu,
            vtxGen=vtxGen,
            rnd=rnd,
        )

        return "ttbar"

    # generate single particles

    if event.startswith("mu_"):
        particle = acts.PdgParticle.eMuon
    elif event.startswith("pi_"):
        particle = acts.PdgParticle.ePionPlus
    elif event.startswith("e_"):
        particle = acts.PdgParticle.eElectron
    else:
        raise ValueError(f"unknown event: {event}")

    pT = float(event.split("_")[1].replace("GeV", ""))

    addParticleGun(
        s,
        ParticleConfig(1, particle, randomizeCharge=True),
        MomentumConfig(pT, pT, transverse=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        PhiConfig(0.0 * u.degree, 360.0 * u.degree),
        vtxGen=vtxGen,
        multiplicity=1,
        rnd=rnd,
        outputDirRoot=outputDirRoot,
    )

    return "single_particles"


def addMySimulation(
    s: acts.examples.Sequencer,
    algorithm: str,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
    detector: Optional[Any] = None,
    inputParticles: str = "particles_input",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    preSelectParticles: Optional[ParticleSelectorConfig] = ParticleSelectorConfig(),
    postSelectParticles: Optional[ParticleSelectorConfig] = None,
    logLevel: Optional[acts.logging.Level] = None,
    **kwargs,
) -> None:
    if algorithm == "fatras":
        addFatras(
            s=s,
            trackingGeometry=trackingGeometry,
            field=field,
            rnd=rnd,
            enableInteractions=True,
            preSelectParticles=preSelectParticles,
            postSelectParticles=postSelectParticles,
            inputParticles=inputParticles,
            outputDirCsv=outputDirCsv,
            outputDirRoot=outputDirRoot,
            logLevel=logLevel,
            **kwargs,
        )
    elif algorithm == "geant4":
        addGeant4(
            s=s,
            detector=detector,
            trackingGeometry=trackingGeometry,
            field=field,
            rnd=rnd,
            inputParticles=inputParticles,
            preSelectParticles=preSelectParticles,
            postSelectParticles=postSelectParticles,
            outputDirCsv=outputDirCsv,
            outputDirRoot=outputDirRoot,
            logLevel=logLevel,
            **kwargs,
        )
    else:
        raise ValueError(f"unknown simulation algorithm: {algorithm}")


def addMySeeding(
    s: acts.examples.Sequencer,
    algortihm: str,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
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
            t0=0.0001 * u.ns,
            phi=1 * u.degree,
            theta=1 * u.degree,
            pRel=0.05,
        ),
        initialSigmas=[0.1, 0.1, 0.002, 0.0001, 0.001, 1000],
        initialVarInflation=[1e3] * 6,
        geoSelectionConfigFile=seedingSel,
        outputDirRoot=outputDirRoot,
    )


if __name__ == "__main__":
    main()
