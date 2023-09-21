#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
import argparse

import acts
from acts.examples.simulation import (
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

from mycommon.events import (
    create_event_label,
    split_event_label,
    get_number_of_events,
    get_event_type,
)
from mycommon.detector import get_odd
from mycommon.reco import addMySeeding


u = acts.UnitConstants

detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("event_label")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads")
    args = parser.parse_args()

    event, simulation = split_event_label(args.event_label)
    event_label = create_event_label(event, simulation)

    indir = Path(args.outdir) / event_label
    outdir = Path(args.outdir) / event_label / "reco"
    outdir.mkdir(parents=True, exist_ok=True)

    events = get_number_of_events(event)
    skip = 0

    with tempfile.TemporaryDirectory() as temp:
        run_reconstruction(args.threads, Path(temp), indir, outdir, events, skip, event)


def run_reconstruction(numThreads, tp, indir, outdir, events, skip, event):
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=events,
        skip=skip,
        numThreads=numThreads,
        trackFpes=False,
    )

    for d in decorators:
        s.addContextDecorator(d)

    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.WARNING,
            particleCollection="particles",
            filePath=str(indir / "particles.root"),
        )
    )
    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.WARNING,
            particleCollection="particles_input",
            filePath=str(indir / "particles.root"),
        )
    )
    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.WARNING,
            particleCollection="particles_selected",
            filePath=str(indir / "particles.root"),
        )
    )

    s.addReader(
        acts.examples.RootSimHitReader(
            level=acts.logging.WARNING,
            simHitCollection="simhits",
            treeName="hits",
            filePath=str(indir / "hits.root"),
        )
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
        "truth_smeared",
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
            nMeasurementsMin=3,
        ),
        outputDirRoot=tp,
    )

    if get_event_type(event) == "ttbar":
        addVertexFitting(
            s,
            field,
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=tp,
        )

    s.run()
    del s

    for stem in [
        "tracksummary_ckf",
        # "trackstates_ckf",
        "performance_ckf",
        "tracksummary_ambi",
        # "trackstates_ambi",
        "performance_ambi",
        "tracksummary_kf",
        # "trackstates_kf",
        "performance_kf",
    ] + (["performance_vertexing"] if get_event_type(event) == "ttbar" else []):
        perf_file = tp / f"{stem}.root"
        assert perf_file.exists(), f"Performance file not found: {perf_file}"
        shutil.copy(perf_file, outdir / f"{stem}.root")


if __name__ == "__main__":
    main()
