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
FORMATS = ["pdf"]


def get_reco_threads(wildcards):
    event_label = wildcards["event_label"]
    event_type = get_event_type(event_label)
    reco_label = wildcards["reco_label"]
    if reco_label == "default" and event_type == "ttbar":
        return int(math.ceil(workflow.cores * 0.2))
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
        result = 100
    if result is None:
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
    prefix="particles|hits",
    skip="[0-9]+",
    events="[0-9]+",
    res_x="|".join(RES_XS),
    res_y="|".join(RES_YS),
    mat_x="|".join(MAT_XS),
    mat_y="|".join(MAT_YS),
    format="|".join(FORMATS),

rule all:
    input:
        expand("data/plots/{reco_label}/{event_label}/pulls_over_eta.csv", reco_label=RECO_LABELS, event_label=EVENT_LABELS),
        expand("data/plots/{reco_label}/{event_label}/efficiency_over_eta.csv", reco_label=RECO_LABELS, event_label=EVENT_LABELS),
        expand("data/plots/{reco_label}/{event_label}/resolution_{res_y}_over_{res_x}.csv", reco_label=RECO_LABELS, event_label=EVENT_LABELS, res_x=RES_XS, res_y=RES_YS),

        "data/event_display/ttbar_200_geant4/hits.csv",
        "data/event_display/truth_smeared/ttbar_200_geant4/tracks.csv",

        expand("plots/detector_layout.{format}", simulation=SIMULATIONS, mat_x=MAT_XS, mat_y=MAT_YS, format=FORMATS),

        expand("plots/sim/material_{simulation}_{mat_y}_vs_{mat_x}.{format}", simulation=SIMULATIONS, mat_x=MAT_XS, mat_y=MAT_YS, format=FORMATS),
        "plots/sim/material_comparison.html",

        expand("plots/sim/{event_label}/particles.{format}", reco_label=RECO_LABELS, event_label=EVENT_LABELS, format=FORMATS),
        expand("plots/sim/{event_label}/nhits_over_eta.{format}", reco_label=RECO_LABELS, event_label=EVENT_LABELS, format=FORMATS),

        expand("plots/reco/{reco_label}/{event_label}/pulls_over_eta_sausage.{format}", reco_label=RECO_LABELS, event_label=EVENT_LABELS, format=FORMATS),
        expand("plots/reco/{reco_label}/{event_label}/inefficiencies.{format}", reco_label=RECO_LABELS, event_label=EVENT_LABELS, format=FORMATS),

        expand("plots/reco/{reco_label}/{single_particle}_{simulation}/pulls_over_eta_errorbars.{format}", reco_label=RECO_LABELS, single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS, format=FORMATS),
        expand("plots/reco/{reco_label}/{single_particle}_{simulation}/resolution_{res_y}_over_{res_x}.{format}", reco_label=RECO_LABELS, single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS, res_x=["eta"], res_y=RES_YS, format=FORMATS),
        expand("plots/reco/{reco_label}/{single_particle}_{simulation}/efficiency_over_eta.{format}", reco_label=RECO_LABELS, single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS, format=FORMATS),

        expand("plots/reco/{reco_label}/single_particles_{pt}_{simulation}/efficiency_over_eta.{format}", reco_label=RECO_LABELS, pt=PT_VALUES, simulation=SIMULATIONS, format=FORMATS),
        expand("plots/reco/{reco_label}/single_particles_{pt_range}_{simulation}/resolution_{res_y}_over_{res_x}.{format}", reco_label=RECO_LABELS, pt_range=PT_RANGES, simulation=SIMULATIONS, res_x=RES_XS, res_y=RES_YS, format=FORMATS),

        expand("plots/reco/{reco_label}/ttbar_{simulation}/efficiency_over_eta.{format}", reco_label=RECO_LABELS, simulation=SIMULATIONS, format=FORMATS),

        "plots/final/single_muon_efficiency.pdf",
        "plots/final/single_particle_efficiency.pdf",
        "plots/final/ttbar_efficiency_ts.pdf",
        "plots/final/ttbar_efficiency_te.pdf",
        "plots/final/single_muon_resolution.pdf",
        "plots/final/single_particle_resolution.pdf",
        "plots/final/single_muon_pulls.pdf",

rule all_sim:
    input:
        expand("data/sim/material_{simulation}/material_tracks.root", simulation=SIMULATIONS),

        expand("data/sim/{event_label}/particles.root", event_label=EVENT_LABELS),
        expand("data/sim/{event_label}/hits.root", event_label=EVENT_LABELS),

rule material_scan:
    input:
        script = "scripts/material_scan.py",
    output:
        "data/sim/material_{simulation}/material_tracks.root",
    shell:
        """
        # ugly macos fix
        export ZSH_VERSION=
        source activate.sh

        # somehow geant4 is crashing when running multiple instances in parallel
        sleep $((RANDOM % 5))

        python {input.script} {wildcards.simulation} --skip 0 --events 100 \
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
            --sub-names all inner beampipe pixel sstrips lstrips solenoid \
            --sub-rmin 0:0:0:25:200:680:1140 \
            --sub-rmax 2000:1140:25:200:680:1140:2000 \
            --sub-zmin -3200:-3200:-3200:-3200:-3200:-3200:-3200 \
            --sub-zmax 3200:3200:3200:3200:3200:3200:3200
        """

rule plot_material:
    input:
        "data/sim/material_{simulation}/material_composition.root",
        script = "scripts/plot/material_generic.py",
    output:
        "plots/sim/material_{simulation}_{mat_y}_vs_{mat_x}.{format}",
    shell:
        "python {input.script} {wildcards.mat_x} {wildcards.mat_y} {input} --output {output}"

rule histcmp_material:
    input:
        expand("data/sim/material_{simulation}/material_composition.root", simulation=SIMULATIONS),
    output:
        "plots/sim/material_comparison.html",
    shell:
        """
        histcmp --label-monitored "Acts" --label-reference "Geant4" --title "ODD material composition" -o {output} {input} || true
        """

rule simulation:
    input:
        get_simulation_slices,
    output:
        "data/sim/{event_label}/{prefix}.root",
    shell:
        "hadd -f {output} {input}"

rule simulation_slice:
    input:
        script = "scripts/simulation.py",
    output:
        "data/sim/{event_label}/slices/{skip}_{events}/particles.root",
        "data/sim/{event_label}/slices/{skip}_{events}/hits.root",
    shell:
        """
        # ugly macos fix
        export ZSH_VERSION=
        source activate.sh

        # somehow geant4 is crashing when running multiple instances in parallel
        sleep $((RANDOM % 5))

        python {input.script} {wildcards.event_label} --skip {wildcards.skip} --events {wildcards.events} \
          data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/ \
          > data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/stdout.txt \
          2> data/sim/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/stderr.txt
        """

rule reconstruction:
    input:
        "data/sim/{event_label}/particles.root",
        "data/sim/{event_label}/hits.root",
        script = "scripts/reconstruction.py",
    output:
        "data/reco/{reco_label}/{event_label}/tracksummary_ambi.root",
    params:
        skip=0,
        events=get_number_of_events,
    threads: get_reco_threads
    shell:
        """
        # ugly macos fix
        export ZSH_VERSION=
        source activate.sh
        
        python {input.script} {wildcards.event_label} {wildcards.reco_label} \
          data/sim/{wildcards.event_label} data/reco/{wildcards.reco_label}/{wildcards.event_label} \
          --skip {params.skip} --events {params.events} --threads {threads} \
          > data/reco/{wildcards.reco_label}/{wildcards.event_label}/stdout.txt \
          2> data/reco/{wildcards.reco_label}/{wildcards.event_label}/stderr.txt
        """

rule truth_matching:
    input:
        files = (
            "data/reco/{reco_label}/{event_label}/tracksummary_ambi.root",
            "data/sim/{event_label}/particles.root",
            "data/sim/{event_label}/hits.root",
        ),
        script = "scripts/truth_matching.py",
    output:
        "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
    threads: 6
    shell:
        "python {input.script} {input.files} {output}"

rule dump_pulls_over_eta:
    input:
        file = "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
        script = "scripts/dump/pulls_over_eta.py",
    output:
        "data/plots/{reco_label}/{event_label}/pulls_over_eta.csv",
    shell:
        "python {input.script} {input.file} {output}"

rule dump_efficiency_over_eta:
    input:
        file = "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
        script = "scripts/dump/efficiency_over_eta.py",
    output:
        "data/plots/{reco_label}/{event_label}/efficiency_over_eta.csv",
    shell:
        "python {input.script} {input.file} {output} --eta-range 0 3 --eta-bins 13"

rule dump_resolution:
    input:
        file = "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
        script = "scripts/dump/resolution_generic.py",
    output:
        "data/plots/{reco_label}/{event_label}/resolution_{res_y}_over_{res_x}.csv",
    shell:
        "python {input.script} {wildcards.res_x} {wildcards.res_y} {input.file} {output} --x-bins 13"

rule plot_pulls_over_eta_sausage:
    input:
        file = "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
        script = "scripts/plot/pulls_over_eta_sausage.py",
    output:
        "plots/reco/{reco_label}/{event_label}/pulls_over_eta_sausage.{format}",
    shell:
        "python {input.script} {input.file} --output {output}"

rule plot_single_particle_pulls_over_eta_errorbars:
    input:
        files = get_all_pt_variants,
        script = "scripts/plot/pulls_over_eta_errorbars.py",
    output:
        "plots/reco/{reco_label}/{single_particle}_{simulation}/pulls_over_eta_errorbars.{format}",
    shell:
        "python {input.script} {input.files} --output {output}"

rule plot_single_particle_resolution:
    input:
        files = get_all_pt_variants,
        script = "scripts/plot/resolution_generic.py",
    output:
        "plots/reco/{reco_label}/{single_particle}_{simulation}/resolution_{res_y}_over_{res_x}.{format}",
    shell:
        "python {input.script} {wildcards.res_x} {wildcards.res_y} {input.files} --output {output}"

rule plot_single_particle_efficiency_over_eta:
    input:
        files = get_all_pt_variants,
        script = "scripts/plot/efficiency_over_eta.py",
    output:
        "plots/reco/{reco_label}/{single_particle}_{simulation}/efficiency_over_eta.{format}",
    shell:
        "python {input.script} {input.files} --output {output}"

rule plot_cross_single_particle_efficiency_over_eta:
    input:
        files = get_all_flavor_variants,
        script = "scripts/plot/efficiency_over_eta.py",
    output:
        "plots/reco/{reco_label}/single_particles_{pt}_{simulation}/efficiency_over_eta.{format}",
    shell:
        "python {input.script} {input.files} --output {output}"

rule plot_ttbar_efficiency_over_eta:
    input:
        files = get_all_ttbar_variants,
        script = "scripts/plot/efficiency_over_eta.py",
    output:
        "plots/reco/{reco_label}/ttbar_{simulation}/efficiency_over_eta.{format}",
    shell:
        "python {input.script} {input.files} --output {output}"

rule plot_cross_single_particle_resolution:
    input:
        files = get_all_pt_range_variants,
        script = "scripts/plot/resolution_generic.py",
    output:
        "plots/reco/{reco_label}/single_particles_{pt_range}_{simulation}/resolution_{res_y}_over_{res_x}.{format}",
    shell:
        "python {input.script} {wildcards.res_x} {wildcards.res_y} {input.files} --output {output}"

rule plot_inefficiencies:
    input:
        file = "data/truth_matching/{reco_label}/{event_label}/truth_matched_tracksummary_ambi.csv",
        script = "scripts/plot/inefficiencies.py",
    output:
        "plots/reco/{reco_label}/{event_label}/inefficiencies.{format}",
    shell:
        "python {input.script} {input.file} --output {output}"

rule plot_particles:
    input:
        file = "data/sim/{event_label}/particles.root",
        script = "scripts/plot/particles.py",
    output:
        "plots/sim/{event_label}/particles.{format}",
    shell:
        "python {input.script} {input.file} --output {output}"

rule plot_nhits_over_eta:
    input:
        files = (
            "data/sim/{event_label}/particles.root",
            "data/sim/{event_label}/hits.root",
        ),
        script = "scripts/plot/nhits_over_eta.py",
    output:
        "plots/sim/{event_label}/nhits_over_eta.{format}",
    shell:
        "python {input.script} {input.files} --output {output}"

rule event_display_reco:
    input:
        "data/sim/{event_label}/particles.root",
        "data/sim/{event_label}/hits.root",
        script = "scripts/reconstruction.py",
    output:
        "data/event_display/{reco_label}/{event_label}/trackstates_ambi.root",
    params:
        skip=0,
        events=1,
    shell:
        """
        python {input.script} {wildcards.event_label} {wildcards.reco_label} \
          data/sim/{wildcards.event_label} data/event_display/{wildcards.reco_label}/{wildcards.event_label} \
          --skip {params.skip} --events {params.events} --threads {threads} --output-trackstates \
          > data/event_display/{wildcards.reco_label}/{wildcards.event_label}/stdout.txt \
          2> data/event_display/{wildcards.reco_label}/{wildcards.event_label}/stderr.txt
        """

rule event_display_dump_hits:
    input:
        file = "data/sim/{event_label}/hits.root",
        script = "scripts/dump/hits.py",
    output:
        "data/event_display/{event_label}/hits.csv",
    params:
        skip=0,
    shell:
        "python {input.script} {input.file} {params.skip} {output}"

rule event_display_dump_tracks:
    input:
        file = "data/event_display/{reco_label}/{event_label}/trackstates_ambi.root",
        script = "scripts/dump/tracks.py",
    output:
        "data/event_display/{reco_label}/{event_label}/tracks.csv",
    params:
        skip=0,
    shell:
        "python {input.script} {input.file} {params.skip} {output}"

rule plot_detector_layout:
    input:
        script = "scripts/plot/detector_layout.py",
    output:
        "plots/detector_layout.{format}",
    shell:
        "python {input.script} --output {output}"

rule plot_final_single_muon_efficiency:
    input:
        files = (
            "data/plots/truth_smeared/mu_1GeV_geant4/efficiency_over_eta.csv",
            "data/plots/truth_smeared/mu_10GeV_geant4/efficiency_over_eta.csv",
            "data/plots/truth_smeared/mu_100GeV_geant4/efficiency_over_eta.csv",
        ),
        script = "scripts/plot/final/single_muon_efficiency.py",
    output:
        "plots/final/single_muon_efficiency.pdf",
    shell:
        "python {input.script} {input.files} {output}"

rule plot_final_single_particle_efficiency:
    input:
        files = (
            "data/plots/truth_smeared/mu_10GeV_geant4/efficiency_over_eta.csv",
            "data/plots/truth_smeared/pi_10GeV_geant4/efficiency_over_eta.csv",
            "data/plots/truth_smeared/e_10GeV_geant4/efficiency_over_eta.csv",
        ),
        script = "scripts/plot/final/single_particle_efficiency.py",
    output:
        "plots/final/single_particle_efficiency.pdf",
    shell:
        "python {input.script} {input.files} {output}"

rule plot_final_ttbar_efficiency_ts:
    input:
        files = (
            "data/plots/truth_smeared/ttbar_60_geant4/efficiency_over_eta.csv",
            "data/plots/truth_smeared/ttbar_120_geant4/efficiency_over_eta.csv",
            "data/plots/truth_smeared/ttbar_200_geant4/efficiency_over_eta.csv",
        ),
        script = "scripts/plot/final/ttbar_efficiency_ts.py",
    output:
        "plots/final/ttbar_efficiency_ts.pdf",
    shell:
        "python {input.script} {input.files} {output}"

rule plot_final_ttbar_efficiency_te:
    input:
        files = (
        "data/plots/truth_estimated/ttbar_60_geant4/efficiency_over_eta.csv",
        "data/plots/truth_estimated/ttbar_120_geant4/efficiency_over_eta.csv",
        "data/plots/truth_estimated/ttbar_200_geant4/efficiency_over_eta.csv",
        ),
        script = "scripts/plot/final/ttbar_efficiency_te.py",
    output:
        "plots/final/ttbar_efficiency_te.pdf",
    shell:
        "python {input.script} {input.files} {output}"

rule plot_final_single_muon_resolution:
    input:
        files = (
            "data/plots/truth_smeared/mu_1GeV_geant4/resolution_d0_over_eta.csv",
            "data/plots/truth_smeared/mu_10GeV_geant4/resolution_d0_over_eta.csv",
            "data/plots/truth_smeared/mu_100GeV_geant4/resolution_d0_over_eta.csv",
        ),
        script = "scripts/plot/final/single_muon_resolution.py",
    output:
        "plots/final/single_muon_resolution.pdf",
    shell:
        """
        python {input.script} {input.files} {output}
        """

rule plot_final_single_particle_resolution:
    input:
        files = (
            "data/plots/truth_smeared/mu_1-100GeV_geant4/resolution_z0_over_pt.csv",
            "data/plots/truth_smeared/pi_1-100GeV_geant4/resolution_z0_over_pt.csv",
            "data/plots/truth_smeared/e_1-100GeV_geant4/resolution_z0_over_pt.csv",
        ),
        script = "scripts/plot/final/single_particle_resolution.py",
    output:
        "plots/final/single_particle_resolution.pdf",
    shell:
        """
        python {input.script} {input.files} {output}
        """

rule plot_final_single_muon_pulls:
    input:
        files = (
            "data/plots/truth_smeared/mu_1GeV_geant4/pulls_over_eta.csv",
            "data/plots/truth_smeared/mu_10GeV_geant4/pulls_over_eta.csv",
            "data/plots/truth_smeared/mu_100GeV_geant4/pulls_over_eta.csv",
        ),
        script = "scripts/plot/final/single_muon_pulls.py",
    output:
        "plots/final/single_muon_pulls.pdf",
    shell:
        """
        python {input.script} {input.files} {output}
        """
