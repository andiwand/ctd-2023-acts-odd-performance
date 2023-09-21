from mycommon.events import list_event_labels, get_event_type

EVENT_LABELS = list_event_labels()

def get_events(event_type):
    if event_type == "ttbar":
        return 10
    return 200000

def get_events_per_slice(event_type):
    if event_type == "ttbar":
        return 1
    return 10000

def get_skip(event_label):
    event_type = get_event_type(event_label)
    total = get_events(event_type)
    step = get_events_per_slice(event_type)
    return range(0, total, step)

def get_simulation_slices(wildcards):
    skip = get_skip(wildcards["event_label"])
    events = get_events(wildcards["event_label"])
    return (
        expand("data/{event_label}/slices/particles_{skip}_{events}.root", event_label=wildcards["event_label"], skip=skip, events=events) +
        expand("data/{event_label}/slices/hits_{skip}_{events}.root", event_label=wildcards["event_label"], skip=skip, events=events)
    )


rule all:
    input:
        expand("data/{event_label}/reco/tracksummary_ckf.root", event_label=EVENT_LABELS),

rule all_sim:
    input:
        expand("data/{event_label}/particles.root", event_label=EVENT_LABELS),
        expand("data/{event_label}/hits.root", event_label=EVENT_LABELS),

rule simulation:
    input:
        get_simulation_slices,
    output:
        "data/{event_label}/particles.root",
        "data/{event_label}/hits.root",
    shell:
        """
        hadd -f data/{wildcards.event_label}/particles.root data/{wildcards.event_label}/slices/particles_*.root
        hadd -f data/{wildcards.event_label}/hits.root data/{wildcards.event_label}/slices/hits_*.root
        """

rule simulation_slice:
    output:
        "data/{event_label}/slices/particles_{skip}_{events}.root",
        "data/{event_label}/slices/hits_{skip}_{events}.root",
    shell:
        "python scripts/simulation.py --skip {wildcards.skip} --events {wildcards.events} data/ {wildcards.event_label}"

rule reconstruction:
    input:
        "data/{event_label}/particles.root",
        "data/{event_label}/hits.root",
    output:
        "data/{event_label}/reco/tracksummary_ckf.root",
    shell:
        "python scripts/reconstruction.py data/ {wildcards.event_label} > data/{wildcards.event_label}/reco/output.txt"
