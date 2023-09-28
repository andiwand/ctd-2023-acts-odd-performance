from pathlib import Path
from typing import Any, Optional, Union

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
)

from mycommon.events import get_event_type, get_event_details

u = acts.UnitConstants


simulations = ["fatras", "geant4"]


def list_simulations():
    return simulations


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

    event_type = get_event_type(event)
    event_details = get_event_details(event)

    if event_type == "ttbar":
        pu = int(event.split("_")[1])

        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=pu,
            vtxGen=vtxGen,
            rnd=rnd,
            outputDirRoot=outputDirRoot,
        )

        return

    if event_type == "single_particles":
        particle_label, pt_label = event_details

        if particle_label == "mu":
            particle = acts.PdgParticle.eMuon
        elif particle_label == "pi":
            particle = acts.PdgParticle.ePionPlus
        elif particle_label == "e":
            particle = acts.PdgParticle.eElectron
        else:
            raise ValueError(f"unknown particle label: {particle_label}")

        if isinstance(pt_label, tuple):
            momentum_config = MomentumConfig(pt_label[0], pt_label[1], transverse=True)
        else:
            momentum_config = MomentumConfig(pt_label, pt_label, transverse=True)

        addParticleGun(
            s,
            ParticleConfig(1, particle, randomizeCharge=True),
            momentum_config,
            EtaConfig(-3.0, 3.0, uniform=True),
            PhiConfig(0.0 * u.degree, 360.0 * u.degree),
            vtxGen=vtxGen,
            multiplicity=1,
            rnd=rnd,
            outputDirRoot=outputDirRoot,
        )

    raise ValueError(f"unknown event type: {event_type}")


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
            killVolume=trackingGeometry.worldVolume,
            killAfterTime=25 * u.ns,
            outputDirCsv=outputDirCsv,
            outputDirRoot=outputDirRoot,
            logLevel=logLevel,
        )
    else:
        raise ValueError(f"unknown simulation algorithm: {algorithm}")
