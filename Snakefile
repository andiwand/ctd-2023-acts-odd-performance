import math
from mycommon.events import (
    list_single_particles,
    list_single_particle_pt_labels,
    list_ttbar_pileups,
    list_simulations,
    list_event_labels,
    get_event_type,
    get_number_of_events,
)


SINGLE_PARTICLES = list_single_particles()
PT_VALUES = list_single_particle_pt_labels()
PILEUP = list_ttbar_pileups()
SIMULATIONS = list_simulations()
EVENT_LABELS = list_event_labels()


def get_reco_threads(wildcards):
    event_type = get_event_type(wildcards["event_label"])
    if event_type == "ttbar":
        return int(math.ceil(workflow.cores * 0.5))
    return int(math.ceil(workflow.cores * 0.1))

def get_events_per_slice(event_type):
    if event_type == "single_particles":
        return 10000
    elif event_type == "ttbar":
        return 100
    raise ValueError(f"Unknown event type: {event_type}")

def get_skip_events(event_label):
    event_type = get_event_type(event_label)
    total = get_number_of_events(event_type)
    step = get_events_per_slice(event_type)
    return range(0, total, step), step

def get_simulation_slices(wildcards):
    skip, events = get_skip_events(wildcards["event_label"])
    return expand(
        "data/sim/{event_label}/slices/{skip}_{events}/{prefix}.root",
        event_label=wildcards["event_label"],
        prefix=wildcards["prefix"],
        skip=skip,
        events=events,
    )

def get_all_pt_variants(wildcards):
    return expand(
        "data/reco/{single_particle}_{pt}_{simulation}/truth_matched_tracksummary_ambi.csv",
        single_particle=wildcards.single_particle,
        pt=PT_VALUES,
        simulation=wildcards.simulation
    )


wildcard_constraints:
    event_label="|".join(EVENT_LABELS),
    single_particle="|".join(SINGLE_PARTICLES),
    pt="|".join(PT_VALUES),
    simulation="|".join(SIMULATIONS),
    prefix="particles|particles_initial|hits",
    skip="[0-9]+",
    events="[0-9]+",

rule all:
    input:
        expand("plots/{event_label}/pulls_over_eta_sausage.png", event_label=EVENT_LABELS),
        expand("plots/{event_label}/efficiency_over_eta.png", event_label=EVENT_LABELS),

        expand("plots/{single_particle}_{simulation}/pulls_over_eta_errorbars.png", single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS),
        expand("plots/{single_particle}_{simulation}/resolution_over_eta.png", single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS),

rule all_sim:
    input:
        expand("data/sim/{event_label}/particles.root", event_label=EVENT_LABELS),
        expand("data/sim/{event_label}/particles_initial.root", event_label=EVENT_LABELS),
        expand("data/sim/{event_label}/hits.root", event_label=EVENT_LABELS),

rule simulation:
    input:
        get_simulation_slices,
    output:
        "data/sim/{event_label}/{prefix}.root",
    shell:
        "hadd -f {output} {input}"

rule simulation_slice:
    output:
        "data/sim/{event_label}/slices/{skip}_{events}/particles.root",
        "data/sim/{event_label}/slices/{skip}_{events}/particles_initial.root",
        "data/sim/{event_label}/slices/{skip}_{events}/hits.root",
        "data/sim/{event_label}/slices/{skip}_{events}/stdout.txt",
        "data/sim/{event_label}/slices/{skip}_{events}/stderr.txt",
    shell:
        """
        mkdir -p data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events} || true
        python scripts/simulation.py {wildcards.event_label} --skip {wildcards.skip} --events {wildcards.events} \
          data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/ \
          > data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/stdout.txt \
          2> data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/stderr.txt
        """

rule reconstruction:
    input:
        "data/sim/{event_label}/particles.root",
        "data/sim/{event_label}/particles_initial.root",
        "data/sim/{event_label}/hits.root",
    output:
        "data/reco/{event_label}/measurements.root",
        "data/reco/{event_label}/tracksummary_ckf.root",
        "data/reco/{event_label}/tracksummary_ambi.root",
        "data/reco/{event_label}/stdout.txt",
        "data/reco/{event_label}/stderr.txt",
    threads: get_reco_threads,
    shell:
        """
        mkdir -p data/reco/{wildcards.event_label} || true
        python scripts/reconstruction.py {wildcards.event_label} data/sim/{wildcards.event_label} data/reco/{wildcards.event_label} --threads {threads} \
          > data/reco/{wildcards.event_label}/stdout.txt \
          2> data/reco/{wildcards.event_label}/stderr.txt
        """

rule truth_matching:
    input:
        "data/reco/{event_label}/tracksummary_ambi.root",
        "data/sim/{event_label}/particles.root",
        "data/sim/{event_label}/hits.root",
    output:
        "data/reco/{event_label}/truth_matched_tracksummary_ambi.csv",
    shell:
        """
        python scripts/truth_matching.py {input} {output}
        """

rule plot_pulls_over_eta_sausage:
    input:
        "data/reco/{event_label}/truth_matched_tracksummary_ambi.csv",
    output:
        "plots/{event_label}/pulls_over_eta_sausage.png",
    shell:
        """
        mkdir -p plots/{wildcards.event_label} || true
        python scripts/plot_pulls_over_eta_sausage.py {input} --output {output}
        """

rule plot_pulls_over_eta_errorbars:
    input:
        get_all_pt_variants,
    output:
        "plots/{single_particle}_{simulation}/pulls_over_eta_errorbars.png",
    shell:
        """
        mkdir -p plots/{wildcards.single_particle}_{wildcards.simulation} || true
        python scripts/plot_pulls_over_eta_errorbars.py {input} --output {output}
        """

rule plot_resolution_over_eta:
    input:
        get_all_pt_variants,
    output:
        "plots/{single_particle}_{simulation}/resolution_over_eta.png",
    shell:
        """
        mkdir -p plots/{wildcards.single_particle}_{wildcards.simulation} || true
        python scripts/plot_resolution_over_eta.py {input} --output {output}
        """

rule plot_efficiency_over_eta:
    input:
        "data/reco/{event_label}/truth_matched_tracksummary_ambi.csv",
    output:
        "plots/{event_label}/efficiency_over_eta.png",
    shell:
        """
        mkdir -p plots/{wildcards.event_label} || true
        python scripts/plot_efficiency_over_eta.py {input} --output {output}
        """
