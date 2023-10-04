import math
from mycommon.events import (
    list_single_particles,
    list_single_particle_pt_labels,
    list_single_particle_pt_range_labels,
    list_ttbar_pileups,
    list_simulations,
    list_event_labels,
    get_event_type,
)
from mycommon.reco import (
    list_reco_labels,
)


SINGLE_PARTICLES = list_single_particles()
PT_VALUES = list_single_particle_pt_labels()
PT_RANGES = list_single_particle_pt_range_labels()
PILEUP = list_ttbar_pileups()
SIMULATIONS = list_simulations()
EVENT_LABELS = list_event_labels()
RECO_LABELS = list_reco_labels()
RES_XS = ["eta", "pt"]
RES_YS = ["d0", "z0", "qop"]
MAT_XS = ["eta", "phi"]
MAT_YS = ["l0", "x0"]


def get_reco_threads(wildcards):
    event_type = get_event_type(wildcards["event_label"])
    if event_type == "ttbar":
        return int(math.ceil(workflow.cores * 0.3))
    return int(math.ceil(workflow.cores * 0.1))

def get_number_of_events(wildcards):
    result = None
    event_label = wildcards["event_label"]
    event_type = get_event_type(event_label)
    if event_type == "single_particles":
        result = 200000
    elif event_type == "ttbar":
        result = 1000
    if result is None:
        raise ValueError(f"unknown event type: {event_type}")
    return int(result * config["event_scale"]["event"] * config["event_scale"]["total"])

def get_events_per_slice(wildcards):
    result = None
    event_label = wildcards["event_label"]
    event_type = get_event_type(event_label)
    if event_type == "single_particles":
        result = 10000
    elif event_type == "ttbar":
        raise ValueError(f"Unknown event type: {event_type}")
    return int(result * config["event_scale"]["slice"] * config["event_scale"]["total"])

def get_skip_events(wildcards):
    total = get_number_of_events(wildcards)
    step = get_events_per_slice(wildcards)
    return range(0, total, step), step

def get_simulation_slices(wildcards):
    skip, events = get_skip_events(wildcards)
    return expand(
        "data/sim/{event_label}/slices/{skip}_{events}/{prefix}.root",
        event_label=wildcards["event_label"],
        prefix=wildcards["prefix"],
        skip=skip,
        events=events,
    )

def get_all_pt_variants(wildcards):
    return expand(
        "data/truth_matching/{reco_label}/{single_particle}_{pt}_{simulation}/truth_matched_tracksummary_ambi.csv",
        reco_label=wildcards.reco_label,
        single_particle=wildcards.single_particle,
        pt=PT_VALUES,
        simulation=wildcards.simulation
    )

def get_all_pt_range_variants(wildcards):
    return expand(
        "data/truth_matching/{reco_label}/{single_particle}_{pt_range}_{simulation}/truth_matched_tracksummary_ambi.csv",
        reco_label=wildcards.reco_label,
        single_particle=SINGLE_PARTICLES,
        pt_range=PT_RANGES,
        simulation=wildcards.simulation
    )

def get_all_flavor_variants(wildcards):
    return expand(
        "data/truth_matching/{reco_label}/{single_particle}_{pt}_{simulation}/truth_matched_tracksummary_ambi.csv",
        reco_label=wildcards.reco_label,
        single_particle=SINGLE_PARTICLES,
        pt=wildcards.pt,
        simulation=wildcards.simulation
    )

def get_all_ttbar_variants(wildcards):
    pileup = list(PILEUP)
    pileup.remove(0)

    return expand(
        "data/truth_matching/{reco_label}/ttbar_{pileup}_{simulation}/truth_matched_tracksummary_ambi.csv",
        reco_label=wildcards.reco_label,
        pileup=pileup,
        simulation=wildcards.simulation
    )


configfile: "config.yaml"

wildcard_constraints:
    event_label="|".join(EVENT_LABELS),
    single_particle="|".join(SINGLE_PARTICLES),
    pt="|".join(PT_VALUES),
    pt_range="|".join(PT_RANGES),
    simulation="|".join(SIMULATIONS),
    reco="|".join(RECO_LABELS),
    prefix="particles|particles_initial|hits",
    skip="[0-9]+",
    events="[0-9]+",
    res_x="|".join(RES_XS),
    res_y="|".join(RES_YS),
    mat_x="|".join(MAT_XS),
    mat_y="|".join(MAT_YS),

rule all:
    input:
        expand("plots/material_comparison_{mat_y}_vs_{mat_x}.png", mat_x=MAT_XS, mat_y=MAT_YS),
        "plots/material_comparison.html",

        expand("plots/{reco_label}/{event_label}/pulls_over_eta_sausage.png", reco_label=RECO_LABELS, event_label=EVENT_LABELS),
        expand("plots/{reco_label}/{event_label}/inefficiencies.png", reco_label=RECO_LABELS, event_label=EVENT_LABELS),
        expand("plots/{reco_label}/{event_label}/particles.png", reco_label=RECO_LABELS, event_label=EVENT_LABELS),

        expand("plots/{reco_label}/{single_particle}_{simulation}/pulls_over_eta_errorbars.png", reco_label=RECO_LABELS, single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS),
        expand("plots/{reco_label}/{single_particle}_{simulation}/resolution_{res_y}_over_{res_x}.png", reco_label=RECO_LABELS, single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS, res_x=RES_XS, res_y=RES_YS),
        expand("plots/{reco_label}/{single_particle}_{simulation}/efficiency_over_eta.png", reco_label=RECO_LABELS, single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS),

        expand("plots/{reco_label}/single_particles_{pt}_{simulation}/efficiency_over_eta.png", reco_label=RECO_LABELS, pt=PT_VALUES, simulation=SIMULATIONS),
        expand("plots/{reco_label}/single_particles_{pt_range}_{simulation}/resolution_{res_y}_over_{res_x}.png", reco_label=RECO_LABELS, pt_range=PT_RANGES, simulation=SIMULATIONS, res_x=RES_XS, res_y=RES_YS),

        expand("plots/{reco_label}/ttbar_{simulation}/efficiency_over_eta.png", reco_label=RECO_LABELS, simulation=SIMULATIONS),

rule all_sim:
    input:
        expand("data/sim/material_{simulation}/material_tracks.root", simulation=SIMULATIONS),

        expand("data/sim/{event_label}/particles.root", event_label=EVENT_LABELS),
        expand("data/sim/{event_label}/particles_initial.root", event_label=EVENT_LABELS),
        expand("data/sim/{event_label}/hits.root", event_label=EVENT_LABELS),

rule material_scan:
    output:
        "data/sim/material_{simulation}/material_tracks.root",
    shell:
        # somehow geant4 is crashing when running multiple instances in parallel
        """
        sleep $((RANDOM % 5))
        mkdir -p data/sim/material_{wildcards.simulation} || true
        python scripts/material_scan.py {wildcards.simulation} --skip 0 --events 100 \
          data/sim/material_{wildcards.simulation}/ \
          > data/sim/material_{wildcards.simulation}/stdout.txt \
          2> data/sim/material_{wildcards.simulation}/stderr.txt
        """

rule material_composition:
    input:
        "data/sim/material_{simulation}/material_tracks.root",
    output:
        "data/sim/material_{simulation}/material_composition.root",
    shell:
        """
        ActsAnalysisMaterialComposition \
            -i {input} -o {output} -s \
            --sub-names all beampipe pixel sstrips lstrips solenoid \
            --sub-rmin 0:0:25:200:680:1100 \
            --sub-rmax 2000:25:200:680:1100:2000 \
            --sub-zmin -3200:-3200:-3200:-3200:-3200:-3200 \
            --sub-zmax 3200:3200:3200:3200:3200:3200 \
        """

rule plot_material:
    input:
        expand("data/sim/material_{simulation}/material_tracks.root", simulation=SIMULATIONS),
    output:
        "plots/material_comparison_{mat_y}_vs_{mat_x}.png",
    shell:
        """
        mkdir -p plots || true
        python scripts/plot_material.py {wildcards.mat_x} {wildcards.mat_y} {input} --output {output}
        """

rule histcmp_material:
    input:
        expand("data/sim/material_{simulation}/material_tracks.root", simulation=SIMULATIONS),
    output:
        "plots/material_comparison.html",
    shell:
        """
        mkdir -p plots || true
        histcmp --label-monitored "acts" --label-reference "geant4" --title "ODD material composition" -o {output} {input}
        """

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
    shell:
        # somehow geant4 is crashing when running multiple instances in parallel
        """
        sleep $((RANDOM % 5))
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
        "data/reco/{reco_label}/{event_label}/tracksummary_ambi.root",
    params:
        skip=0,
        events=get_number_of_events,
    threads: get_reco_threads,
    shell:
        """
        mkdir -p data/reco/{wildcards.reco_label}/{wildcards.event_label} || true
        python scripts/reconstruction.py {wildcards.event_label} {wildcards.reco_label} \
          data/sim/{wildcards.event_label} data/reco/{wildcards.reco_label}/{wildcards.event_label} \
          --skip {params.skip} --events {params.events} --threads {threads} \
          > data/reco/{wildcards.reco_label}/{wildcards.event_label}/stdout.txt \
          2> data/reco/{wildcards.reco_label}/{wildcards.event_label}/stderr.txt
        """

rule truth_matching:
    input:
        "data/reco/{reco_label}/{event_label}/tracksummary_ambi.root",
        "data/sim/{event_label}/particles.root",
        "data/sim/{event_label}/hits.root",
    output:
        "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
    shell:
        """
        mkdir -p data/truth_matching/{wildcards.reco_label}/{wildcards.event_label} || true
        python scripts/truth_matching.py {input} {output}
        """

rule plot_pulls_over_eta_sausage:
    input:
        "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
    output:
        "plots/{reco_label}/{event_label}/pulls_over_eta_sausage.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/{wildcards.event_label} || true
        python scripts/plot_pulls_over_eta_sausage.py {input} --output {output}
        """

rule plot_single_particle_pulls_over_eta_errorbars:
    input:
        get_all_pt_variants,
    output:
        "plots/{reco_label}/{single_particle}_{simulation}/pulls_over_eta_errorbars.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/{wildcards.single_particle}_{wildcards.simulation} || true
        python scripts/plot_pulls_over_eta_errorbars.py {input} --output {output}
        """

rule plot_single_particle_resolution:
    input:
        get_all_pt_variants,
    output:
        "plots/{reco_label}/{single_particle}_{simulation}/resolution_{res_y}_over_{res_x}.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/{wildcards.single_particle}_{wildcards.simulation} || true
        python scripts/plot_resolution_generic.py {wildcards.res_x} {wildcards.res_y} {input} --output {output}
        """

rule plot_single_particle_efficiency_over_eta:
    input:
        get_all_pt_variants,
    output:
        "plots/{reco_label}/{single_particle}_{simulation}/efficiency_over_eta.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/{wildcards.single_particle}_{wildcards.simulation} || true
        python scripts/plot_efficiency_over_eta.py {input} --output {output}
        """

rule plot_cross_single_particle_efficiency_over_eta:
    input:
        get_all_flavor_variants,
    output:
        "plots/{reco_label}/single_particles_{pt}_{simulation}/efficiency_over_eta.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/single_particles_{wildcards.pt}_{wildcards.simulation} || true
        python scripts/plot_efficiency_over_eta.py {input} --output {output}
        """

rule plot_ttbar_efficiency_over_eta:
    input:
        get_all_ttbar_variants,
    output:
        "plots/{reco_label}/ttbar_{simulation}/efficiency_over_eta.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/ttbar_{wildcards.simulation} || true
        python scripts/plot_efficiency_over_eta.py {input} --output {output}
        """

rule plot_cross_single_particle_resolution:
    input:
        get_all_pt_range_variants,
    output:
        "plots/{reco_label}/single_particles_{pt_range}_{simulation}/resolution_{res_y}_over_{res_x}.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/single_particles_{wildcards.pt_range}_{wildcards.simulation} || true
        python scripts/plot_resolution_generic.py {wildcards.res_x} {wildcards.res_y} {input} --output {output}
        """

rule plot_inefficiencies:
    input:
        "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
    output:
        "plots/{reco_label}/{event_label}/inefficiencies.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/{wildcards.event_label} || true
        python scripts/plot_inefficiencies.py {input} --output {output}
        """

rule plot_particles:
    input:
        "data/sim/{event_label}/particles.root",
    output:
        "plots/{reco_label}/{event_label}/particles.png",
    shell:
        """
        mkdir -p plots/{wildcards.reco_label}/{wildcards.event_label} || true
        python scripts/plot_particles.py {input} --output {output}
        """
