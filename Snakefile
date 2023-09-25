from mycommon.events import list_event_labels, get_event_type, get_number_of_events


EVENT_LABELS = list_event_labels()
SINGLE_PARTICLES = ["mu", "pi", "e"]
PT_VALUES = ["1GeV", "10GeV", "100GeV"]
SIMULATIONS = ["fatras", "geant4"]
RECO_THREADS = 4


def get_events_per_slice(event_type):
    if event_type == "single_particles":
        return 10000
    elif event_type == "ttbar":
        return 1
    raise ValueError(f"Unknown event type: {event_type}")

def get_skip_events(event_label):
    event_type = get_event_type(event_label)
    total = get_number_of_events(event_type)
    step = get_events_per_slice(event_type)
    return range(0, total, step), step

def get_simulation_slices(wildcards):
    skip, events = get_skip_events(wildcards["event_label"])
    return expand(
        "data/{event_label}/slices/{skip}_{events}/{prefix}.root",
        event_label=wildcards["event_label"],
        prefix=wildcards["prefix"],
        skip=skip,
        events=events,
    )

def get_all_pt_variants(wildcards):
    return expand(
        "data/{single_particle}_{pt}_{simulation}/reco/tracksummary_ckf.root",
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
        expand("plots/{single_particle}_{simulation}/{pt}/pulls_over_eta_sausage.png", single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS, pt=PT_VALUES),
        expand("plots/{single_particle}_{simulation}/pulls_over_eta_errorbars.png", single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS),
        expand("plots/{single_particle}_{simulation}/resolution_over_eta.png", single_particle=SINGLE_PARTICLES, simulation=SIMULATIONS),

rule all_sim:
    input:
        expand("data/{event_label}/particles.root", event_label=EVENT_LABELS),
        expand("data/{event_label}/particles_initial.root", event_label=EVENT_LABELS),
        expand("data/{event_label}/hits.root", event_label=EVENT_LABELS),

rule simulation:
    input:
        get_simulation_slices,
    output:
        "data/{event_label}/{prefix}.root",
    shell:
        "hadd -f {output} {input}"

rule simulation_slice:
    output:
        "data/{event_label}/slices/{skip}_{events}/particles.root",
        "data/{event_label}/slices/{skip}_{events}/particles_initial.root",
        "data/{event_label}/slices/{skip}_{events}/hits.root",
        "data/{event_label}/slices/{skip}_{events}/stdout.txt",
        "data/{event_label}/slices/{skip}_{events}/stderr.txt",
    shell:
        """
        mkdir -p data/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events} || true
        python scripts/simulation.py {wildcards.event_label} --skip {wildcards.skip} --events {wildcards.events} \
          data/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/ \
          > data/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/stdout.txt \
          2> data/{wildcards.event_label}/slices/{wildcards.skip}_{wildcards.events}/stderr.txt
        """

rule reconstruction:
    input:
        "data/{event_label}/particles.root",
        "data/{event_label}/particles_initial.root",
        "data/{event_label}/hits.root",
    output:
        "data/{event_label}/reco/measurements.root",
        "data/{event_label}/reco/tracksummary_ckf.root",
        "data/{event_label}/reco/stdout.txt",
        "data/{event_label}/reco/stderr.txt",
    shell:
        """
        mkdir -p data/{wildcards.event_label}/reco || true
        python scripts/reconstruction.py {wildcards.event_label} data/{wildcards.event_label} data/{wildcards.event_label}/reco \
          > data/{wildcards.event_label}/reco/stdout.txt \
          2> data/{wildcards.event_label}/reco/stderr.txt
        """

rule plot_pulls_over_eta_sausage:
    input:
        "data/{single_particle}_{pt}_{simulation}/reco/tracksummary_ckf.root",
    output:
        "plots/{single_particle}_{simulation}/{pt}/pulls_over_eta_sausage.png",
    shell:
        """
        mkdir -p plots/{wildcards.single_particle}_{wildcards.simulation}/{wildcards.pt} || true
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
        "data/{single_particle}_{pt}_{simulation}/reco/tracksummary_ckf.root",
        "data/{single_particle}_{pt}_{simulation}/particles.root",
        "data/{single_particle}_{pt}_{simulation}/hits.root",
    output:
        "plots/{single_particle}_{simulation}/{pt}/efficiency_over_eta.png",
    shell:
        """
        mkdir -p plots/{wildcards.single_particle}_{wildcards.simulation}/{wildcards.pt} || true
        python scripts/plot_efficiency_over_eta.py \
          "data/{wildcards.single_particle}_{wildcards.pt}_{wildcards.simulation}/reco/tracksummary_ckf.root" \
          --particles "data/{wildcards.single_particle}_{wildcards.pt}_{wildcards.simulation}/particles.root" \
          --hits "data/{wildcards.single_particle}_{wildcards.pt}_{wildcards.simulation}/hits.root" \
          --output {output}
        """
