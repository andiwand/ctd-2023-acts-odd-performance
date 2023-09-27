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

    # generate special events on top

    if event.startswith("ttbar_"):
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
