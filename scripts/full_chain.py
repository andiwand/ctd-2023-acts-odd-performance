#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
import itertools
from multiprocessing import Pool, freeze_support
import argparse
import fnmatch

import acts
from acts.examples.simulation import (
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addKalmanTracks,
    addTrajectoryWriters,
    addCKFTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    VertexFinder,
    addVertexFitting,
)

from mycommon.events import create_event_label, get_number_of_events
from mycommon.detector import get_odd
from mycommon.sim import addMyEventGen, addMySimulation
from mycommon.reco import addMySeeding


u = acts.UnitConstants

detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()


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
    ] + [f"ttbar_{pu}" for pu in [0, 60, 120, 200]]
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
            label = create_event_label(event, simulation, seeding)

            if not fnmatch.fnmatch(label, args.filter):
                print(f"skipping {label}")
                continue

            outdir = Path(args.outdir) / label
            outdir.mkdir(parents=True, exist_ok=True)

            events = get_number_of_events(event)
            skip = 0

            # for debugging
            # run_full_chain(outdir, events, skip, event, simulation, seeding)

            pool.apply_async(
                run_full_chain,
                args=(outdir, events, skip, event, simulation, seeding),
            )

        pool.close()
        pool.join()


def run_full_chain(outdir, events, skip, event, simulation, seeding):
    # single thread because we do multiprocessing
    numThreads = 1

    with tempfile.TemporaryDirectory() as temp:
        tp = Path(temp)
        rnd = acts.examples.RandomNumbers(seed=42)

        s = acts.examples.Sequencer(
            events=events,
            skip=skip,
            numThreads=numThreads,
            trackFpes=False,
        )

        for d in decorators:
            s.addContextDecorator(d)

        addMyEventGen(
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
            geoSelectionConfigFile=seedingSel,
            outputDirRoot=tp,
        )

        addKalmanTracks(
            s,
            trackingGeometry,
            field,
            reverseFilteringMomThreshold=100 * u.TeV,
            # logLevel=acts.logging.VERBOSE,
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
            # logLevel=acts.logging.VERBOSE,
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
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=tp,
        )

        s.run()
        del s

        for stem in [
            # "particles",
            # "particles_initial",
            # "particles_final",
            "hits",
            "measurements",
            "tracksummary_ckf",
            # "trackstates_ckf",
            "performance_ckf",
            "tracksummary_ambi",
            # "trackstates_ambi",
            "performance_ambi",
            "tracksummary_kf",
            # "trackstates_kf",
            "performance_kf",
            "performance_vertexing",
        ]:
            perf_file = tp / f"{stem}.root"
            assert perf_file.exists(), f"Performance file not found: {perf_file}"
            shutil.copy(perf_file, outdir / f"{stem}.root")


if __name__ == "__main__":
    main()
