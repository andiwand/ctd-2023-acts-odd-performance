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
    addTruthTrackingGsf,
    addTrajectoryWriters,
    TrackFindingConfig,
    addCKFTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    VertexFinder,
    addVertexFitting,
)

from mycommon.events import (
    split_event_label,
    get_number_of_events,
    get_event_type,
)
from mycommon.detector import get_odd
from mycommon.reco import split_reco_label, addMySeeding


u = acts.UnitConstants

detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("event_label")
    parser.add_argument("reco_label")
    parser.add_argument("indir")
    parser.add_argument("outdir")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    event, simulation = split_event_label(args.event_label)
    event_type = get_event_type(event)
    seeding = split_reco_label(args.reco_label)

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    skip = 0
    events = get_number_of_events(event_type)

    with tempfile.TemporaryDirectory() as temp:
        run_reconstruction(args.threads, Path(temp), event, seeding, indir, outdir, skip, events)


def run_reconstruction(numThreads, tp, event, seeding, indir, outdir, skip, events):
    event_type = get_event_type(event)
    is_single_electrons = event.startswith("e_")
    is_ttbar = event_type == "ttbar"

    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=events,
        skip=skip,
        numThreads=numThreads,
        trackFpes=False,
        outputDir=tp,
    )

    for d in decorators:
        s.addContextDecorator(d)

    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.WARNING,
            particleCollection="particles",
            filePath=indir / "particles.root",
        )
    )
    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.WARNING,
            particleCollection="particles_input",
            filePath=indir / "particles.root",
        )
    )
    s.addReader(
        acts.examples.RootParticleReader(
            level=acts.logging.WARNING,
            particleCollection="particles_selected",
            filePath=indir / "particles.root",
        )
    )

    s.addReader(
        acts.examples.RootSimHitReader(
            level=acts.logging.WARNING,
            simHitCollection="simhits",
            treeName="hits",
            filePath=indir / "hits.root",
        )
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfig,
        rnd=rnd,
        # outputDirRoot=tp,
    )

    addMySeeding(
        s,
        seeding,
        trackingGeometry,
        field,
        rnd=rnd,
        geoSelectionConfigFile=seedingSel,
        # outputDirRoot=tp,
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

    if is_single_electrons:
        addTruthTrackingGsf(
            s,
            trackingGeometry,
            field,
        )
        addTrajectoryWriters(
            s,
            name="gsf",
            trajectories="gsf_trajectories",
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
        trackFindingConfig=TrackFindingConfig(
            chi2CutOff=15.0,
            numMeasurementsCutOff=10,
        ),
        # outputDirRoot=tp,
    )

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3,
            maximumIterations=100000,
            nMeasurementsMin=3,
        ),
        outputDirRoot=tp,
    )

    if is_ttbar:
        addVertexFitting(
            s,
            field,
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=tp,
        )

    s.run()
    del s

    outdir.mkdir(parents=True, exist_ok=True)
    for file in (
        [
            "timing.tsv",
            # "measurements.root",
            # "tracksummary_ckf.root",
            # "trackstates_ckf.root",
            # "performance_ckf.root",
            "tracksummary_ambi.root",
            # "trackstates_ambi.root",
            "performance_ambi.root",
            "tracksummary_kf.root",
            # "trackstates_kf.root",
            "performance_kf.root",
        ]
        + (
            [
                "tracksummary_gsf.root",
                # "trackstates_gsf.root",
                "performance_gsf.root",
            ]
            if is_single_electrons
            else []
        )
        + (["performance_vertexing.root"] if is_ttbar else [])
    ):
        source = tp / file
        destination = outdir / file
        assert source.exists(), f"File not found: {source}"
        shutil.copy(source, destination)


if __name__ == "__main__":
    main()
