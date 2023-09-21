from mycommon.events import list_event_labels

EVENT_LABELS = list_event_labels()

rule all:
    input:
        expand("data/{event_label}/reco/tracksummary_ckf.root", event_label=EVENT_LABELS),

rule all_sim:
    input:
        expand("data/{event_label}/particles.root", event_label=EVENT_LABELS),
        expand("data/{event_label}/hits.root", event_label=EVENT_LABELS),

rule simulation:
    output:
        "data/{event_label}/particles.root",
        "data/{event_label}/hits.root",
    shell:
        "python scripts/simulation.py data/ {wildcards.event_label} > data/{wildcards.event_label}/sim-stdout.txt 2> data/{wildcards.event_label}/sim-stderr.txt"

rule reconstruction:
    input:
        "data/{event_label}/particles.root",
        "data/{event_label}/hits.root",
    output:
        "data/{event_label}/reco/tracksummary_ckf.root",
    shell:
        "python scripts/reconstruction.py data/ {wildcards.event_label} > data/{wildcards.event_label}/reco/reco-stdout.txt 2> data/{wildcards.event_label}/reco/reco-stderr.txt"
