from mycommon.events import list_event_labels, get_event_type, get_number_of_events


EVENT_LABELS = list_event_labels()
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
        "data/{event_label}/slices/{prefix}_{skip}_{events}.root",
        event_label=wildcards["event_label"],
        prefix=wildcards["prefix"],
        skip=skip,
        events=events,
    )


rule all:
    input:
        expand("data/{event_label}/reco/tracksummary_ckf.root", event_label=EVENT_LABELS),

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
    wildcard_constraints:
        prefix="particles|particles_initial|hits",
    shell:
        "hadd -f {output} {input}"

rule simulation_slice:
    output:
        "data/{event_label}/slices/particles_{skip}_{events}.root",
        "data/{event_label}/slices/particles_initial_{skip}_{events}.root",
        "data/{event_label}/slices/hits_{skip}_{events}.root",
    shell:
        "python scripts/simulation.py --skip {wildcards.skip} --events {wildcards.events} data/ {wildcards.event_label}"

rule reconstruction:
    input:
        "data/{event_label}/particles.root",
        "data/{event_label}/particles_initial.root",
        "data/{event_label}/hits.root",
    output:
        "data/{event_label}/reco/tracksummary_ckf.root",
    shell:
        "python scripts/reconstruction.py data/ {wildcards.event_label} > data/{wildcards.event_label}/reco/output.txt"
