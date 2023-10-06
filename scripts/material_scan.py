#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
import argparse

import acts
import acts.examples
import acts.examples.geant4
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    getG4DetectorConstructionFactory,
)

u = acts.UnitConstants

from mycommon.detector import get_odd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("simulation")
    parser.add_argument("outdir")
    parser.add_argument("--skip", type=int, required=True, help="Skip number of events")
    parser.add_argument("--events", type=int, required=True, help="Number of events")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    simulation = args.simulation
    skip = args.skip
    events = args.events

    with tempfile.TemporaryDirectory() as temp:
        run_material_scan(Path(temp), outdir, events, skip, simulation)

    return 0


def run_material_scan(tp, outdir, events, skip, simulation):
    detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()
    tracks_per_event = 1000
    # single thread because of G4
    numThreads = 1

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

    if simulation == "fatras":
        nav = acts.Navigator(trackingGeometry=trackingGeometry)
        stepper = acts.StraightLineStepper()
        prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

        s.addAlgorithm(
            acts.examples.PropagationAlgorithm(
                propagatorImpl=prop,
                level=acts.logging.INFO,
                randomNumberSvc=rnd,
                ntests=tracks_per_event,
                sterileLogger=True,
                propagationStepCollection="propagation_steps",
                propagationMaterialCollection="material_tracks",
                particleHypothesis=acts.ParticleHypothesis.geantino,
                energyLoss=False,
                multipleScattering=False,
                recordMaterialInteractions=True,
                phiRange=(0 * u.degree, 360 * u.degree),
                etaRange=(-4, 4),
                ptRange=(1.0 * u.GeV, 1.0 * u.GeV),
                d0Sigma=0,
                z0Sigma=50 * u.mm,
                phiSigma=0,
                thetaSigma=0,
                qpSigma=0,
                tSigma=0,
            )
        )
    elif simulation == "geant4":
        addParticleGun(
            s,
            ParticleConfig(num=1, pdg=acts.PdgParticle.eInvalid, charge=0, mass=0),
            MomentumConfig(1.0 * u.GeV, 1.0 * u.GeV, transverse=True),
            EtaConfig(-4.0, 4.0, uniform=True),
            PhiConfig(0.0 * u.degree, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(0, 0, 50 * u.mm, 0),
            ),
            multiplicity=tracks_per_event,
            rnd=rnd,
            outputDirRoot=tp,
        )

        s.addAlgorithm(
            acts.examples.geant4.Geant4MaterialRecording(
                level=acts.logging.INFO,
                inputParticles="particles_input",
                outputMaterialTracks="material_tracks",
                detectorConstructionFactory=getG4DetectorConstructionFactory(detector),
                randomNumbers=rnd,
                excludeMaterials=[],
            )
        )
    else:
        raise ValueError(f"unknown simulation algorithm: {simulation}")

    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            collapseInteractions=True,
            collection="material_tracks",
            filePath=tp / "material_tracks.root",
            level=acts.logging.INFO,
        )
    )

    s.run()
    del s

    outdir.mkdir(parents=True, exist_ok=True)
    for file in [
        "timing.tsv",
        # "particles.root",
        "material_tracks.root",
    ]:
        source = tp / file
        destination = outdir / file
        assert source.exists(), f"File not found: {source}"
        shutil.copy(source, destination)


if __name__ == "__main__":
    main()
