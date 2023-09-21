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
    create_event_label,
    split_event_label,
    get_number_of_events,
    get_event_type,
)
from mycommon.detector import get_odd
from mycommon.sim import addMyEventGen, addMySimulation


u = acts.UnitConstants

detector, trackingGeometry, decorators, field, digiConfig, seedingSel = get_odd()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("event_label")
    args = parser.parse_args()

    event, simulation = split_event_label(args.event_label)
    event_label = create_event_label(event, simulation)

    outdir = Path(args.outdir) / event_label
    outdir.mkdir(parents=True, exist_ok=True)

    events = get_number_of_events(event)
    skip = 0

    with tempfile.TemporaryDirectory() as temp:
        run_simulation(Path(temp), outdir, events, skip, event, simulation)


def run_simulation(tp, outdir, events, skip, event, simulation):
    # single thread because of G4
    numThreads = 1

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

    s.run()
    del s

    if get_event_type(event) == "ttbar":
        shutil.copy(tp / "pythia8_particles.root", tp / "particles.root")

    for stem in [
        "particles",
        # "particles_initial",
        # "particles_final",
        "hits",
    ]:
        perf_file = tp / f"{stem}.root"
        assert perf_file.exists(), f"Performance file not found: {perf_file}"
        shutil.copy(perf_file, outdir / f"{stem}.root")


if __name__ == "__main__":
    main()
