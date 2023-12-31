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
    addTrackWriters,
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("event_label")
    parser.add_argument("reco_label")
    parser.add_argument("indir")
    parser.add_argument("outdir")
    parser.add_argument("--skip", type=int, required=True, help="Skip number of events")
    parser.add_argument("--events", type=int, required=True, help="Number of events")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--output-trackstates", action="store_true")
    args = parser.parse_args()

    event, simulation = split_event_label(args.event_label)
    seeding = split_reco_label(args.reco_label)

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    skip = args.skip
    events = args.events

    with tempfile.TemporaryDirectory() as temp:
        run_reconstruction(
            args.threads,
            Path(temp),
            event,
            seeding,
            indir,
            outdir,
            skip,
            events,
            args.output_trackstates,
        )


def run_reconstruction(
    threads,
    tp,
    event,
    seeding,
    indir,
    outdir,
    skip,
    events,
    output_trackstates=False,
):
    detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()

    event_type = get_event_type(event)
    is_single_electrons = event.startswith("e_")
    is_ttbar = event_type == "ttbar"
    is_truth_seeding = seeding.startswith("truth_")

    reco_config = get_reco_config(event, seeding)

    output_files = []

    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=events,
        skip=skip,
        numThreads=threads,
        trackFpes=False,
        outputDir=tp,
    )
    output_files.append("timing.tsv")

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
    # output_files.append("measurements.root")

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
        ckfConfig=reco_config.ckf_config,
        # outputDirRoot=tp,
        # logLevel=acts.logging.VERBOSE,
    )
    # output_files.append("tracksummary_ckf.root")
    # output_files.append("trackstates_ckf.root")
    # output_files.append("performance_ckf.root")

    addAmbiguityResolution(
        s,
        config=reco_config.ambi_config,
        # outputDirRoot=tp,
    )
    addTrackWriters(
        s,
        name="ambi",
        tracks="tracks",
        outputDirRoot=tp,
        writeStates=output_trackstates,
        writeSummary=True,
        writeCKFperformance=True,
        writeFinderPerformance=False,
        writeFitterPerformance=False,
    )
    output_files.append("tracksummary_ambi.root")
    if output_trackstates:
        output_files.append("trackstates_ambi.root")
    output_files.append("performance_ambi.root")

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
        addTrackWriters(
            s,
            name="kf",
            tracks="kfTracks",
            outputDirRoot=tp,
            writeStates=False,
            writeSummary=True,
            writeCKFperformance=True,
            writeFinderPerformance=False,
            writeFitterPerformance=False,
        )
        output_files.append("tracksummary_kf.root")
        output_files.append("trackstates_kf.root")
        output_files.append("performance_kf.root")
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
        addTrackWriters(
            s,
            name="gsf",
            tracks="gsfTracks",
            outputDirRoot=tp,
            writeStates=False,
            writeSummary=True,
            writeCKFperformance=True,
            writeFinderPerformance=False,
            writeFitterPerformance=False,
        )
        output_files.append("tracksummary_gsf.root")
        output_files.append("trackstates_gsf.root")
        output_files.append("performance_gsf.root")
    """

    if is_ttbar:
        addVertexFitting(
            s,
            field,
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=tp,
        )
        output_files.append("performance_vertexing.root")

    s.run()
    del s

    outdir.mkdir(parents=True, exist_ok=True)
    for file in output_files:
        source = tp / file
        destination = outdir / file
        assert source.exists(), f"File not found: {source}"
        shutil.copy(source, destination)


if __name__ == "__main__":
    main()
