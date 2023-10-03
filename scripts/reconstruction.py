#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
import argparse

import acts
from acts.examples.simulation import (
    addDigitization,
    addParticleSelection,
    ParticleSelectorConfig,
)
from acts.examples.reconstruction import (
    addTrajectoryWriters,
    addCKFTracks,
    addAmbiguityResolution,
    VertexFinder,
    addVertexFitting,
)

u = acts.UnitConstants

from mycommon.events import (
    split_event_label,
    get_event_type,
)
from mycommon.detector import get_odd
from mycommon.reco import split_reco_label, addMySeeding, get_reco_config

detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("event_label")
    parser.add_argument("reco_label")
    parser.add_argument("indir")
    parser.add_argument("outdir")
    parser.add_argument("--skip", type=int, required=True, help="Skip number of events")
    parser.add_argument("--events", type=int, required=True, help="Number of events")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    event, simulation = split_event_label(args.event_label)
    seeding = split_reco_label(args.reco_label)

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    skip = args.skip
    events = args.events

    with tempfile.TemporaryDirectory() as temp:
        run_reconstruction(
            args.threads, Path(temp), event, seeding, indir, outdir, skip, events
        )


def run_reconstruction(numThreads, tp, event, seeding, indir, outdir, skip, events):
    event_type = get_event_type(event)
    is_single_electrons = event.startswith("e_")
    is_ttbar = event_type == "ttbar"
    is_truth_seeding = seeding.startswith("truth_")

    reco_config = get_reco_config(event, seeding)

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
            particleCollection="particles_input",
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

    addParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            # using something close to 1 to include for sure
            pt=(0.999 * u.GeV, None),
            removeNeutral=True,
        ),
        "particles_input",
        "particles_selected",
    )
    s.addWhiteboardAlias("particles", "particles_selected")

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

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        trackSelectorConfig=reco_config.track_selector_config,
        trackFindingConfig=reco_config.track_finding_config,
        # outputDirRoot=tp,
    )

    addAmbiguityResolution(
        s,
        config=reco_config.ambi_config,
        # outputDirRoot=tp,
    )
    addTrajectoryWriters(
        s,
        name="ambi",
        trajectories="trajectories",
        outputDirRoot=tp,
        writeStates=False,
        writeSummary=True,
        writeCKFperformance=True,
        writeFinderPerformance=False,
        writeFitterPerformance=False,
    )

    # TODO broken; needs fixing in acts examples
    """
    if is_truth_seeding:
        kfOptions = {
            "multipleScattering": True,
            "energyLoss": True,
            "reverseFilteringMomThreshold": float("inf"),
            "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(),
            "level": acts.logging.VERBOSE,
        }
        s.addAlgorithm(
            acts.examples.RefittingAlgorithm(
                acts.logging.VERBOSE,
                inputTracks="ambiTracks",
                outputTracks="kfTracks",
                fit=acts.examples.makeKalmanFitterFunction(trackingGeometry, field, **kfOptions),
            )
        )
        s.addAlgorithm(
            acts.examples.TracksToTrajectories(
                level=acts.logging.INFO,
                inputTracks="kfTracks",
                outputTrajectories="kfTrajectories",
            )
        )
        addTrajectoryWriters(
            s,
            name="kf",
            trajectories="kfTrajectories",
            outputDirRoot=tp,
            writeStates=False,
            writeSummary=True,
            writeCKFperformance=True,
            writeFinderPerformance=False,
            writeFitterPerformance=False,
        )
    """

    # TODO broken; needs fixing in acts examples
    """
    if is_single_electrons and is_truth_seeding:
        gsfOptions = {
            "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
            "maxComponents": 4,
            "abortOnError": False,
            "disableAllMaterialHandling": False,
            "finalReductionMethod": acts.examples.FinalReductionMethod.maxWeight,
            "weightCutoff": 1.0e-4,
            "level": acts.logging.INFO,
        }
        s.addAlgorithm(
            acts.examples.RefittingAlgorithm(
                acts.logging.INFO,
                inputTracks="ambiTracks",
                outputTracks="gsfTracks",
                fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions),
            )
        )
        s.addAlgorithm(
            acts.examples.TracksToTrajectories(
                level=acts.logging.INFO,
                inputTracks="gsfTracks",
                outputTrajectories="gsfTrajectories",
            )
        )
        addTrajectoryWriters(
            s,
            name="gsf",
            trajectories="gsfTrajectories",
            outputDirRoot=tp,
            writeStates=False,
            writeSummary=True,
            writeCKFperformance=True,
            writeFinderPerformance=False,
            writeFitterPerformance=False,
        )
    """

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
        ]
        + (
            [
                # "tracksummary_kf.root",
                # "trackstates_kf.root",
                # "performance_kf.root",
            ]
            if is_truth_seeding
            else []
        )
        + (
            [
                # "tracksummary_gsf.root",
                # "trackstates_gsf.root",
                # "performance_gsf.root",
            ]
            if is_single_electrons and is_truth_seeding
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
