#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
import argparse

import acts
from acts.examples.simulation import (
    ParticleSelectorConfig,
)

from mycommon.events import (
    split_event_label,
    get_event_type,
)
from mycommon.detector import get_odd
from mycommon.sim import addMyEventGen, addMySimulation


u = acts.UnitConstants

detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("event_label")
    parser.add_argument("outdir")
    parser.add_argument("--skip", type=int, required=True, help="Skip number of events")
    parser.add_argument("--events", type=int, required=True, help="Number of events")
    args = parser.parse_args()

    event, simulation = split_event_label(args.event_label)

    outdir = Path(args.outdir)
    skip = args.skip
    events = args.events

    with tempfile.TemporaryDirectory() as temp:
        run_simulation(Path(temp), event, outdir, events, skip, simulation)

    return 0


def run_simulation(tp, event, outdir, events, skip, simulation):
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
        preSelectParticles=ParticleSelectorConfig(
            # these cuts are necessary because of pythia
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
        ),
        postSelectParticles=ParticleSelectorConfig(
            # these cuts should not be necessary for sim
            eta=(-3.0, 3.0),
            pt=(1 * u.GeV, None),
            removeNeutral=True,
        ),
        outputDirRoot=tp,
    )

    s.run()
    del s

    if get_event_type(event) == "ttbar":
        shutil.copy(tp / "pythia8_particles.root", tp / "particles.root")

    outdir.mkdir(parents=True, exist_ok=True)
    for file in [
        "timing.tsv",
        "particles.root",
        "particles_initial.root",
        # "particles_final.root",
        "hits.root",
    ]:
        source = tp / file
        destination = outdir / file
        assert source.exists(), f"File not found: {source}"
        shutil.copy(source, destination)


if __name__ == "__main__":
    main()
