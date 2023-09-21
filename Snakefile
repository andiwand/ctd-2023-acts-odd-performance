from mycommon.events import list_event_labels

DATA_DIR = "data3"
EVENT_LABELS = list_event_labels()

rule all:
    input:
        expand("data3/{event_label}/reco/tracksummary_ckf.root", event_label=EVENT_LABELS),

rule all_simulation:
    input:
        expand("data3/{event_label}/reco/particles.root", event_label=EVENT_LABELS),
        expand("data3/{event_label}/reco/hits.root", event_label=EVENT_LABELS),

rule simulation:
    output:
        "data3/{event_label}/particles.root",
        "data3/{event_label}/hits.root"
    shell:
        "python scripts/simulation.py data3/ {wildcards.event_label} > data3/{wildcards.event_label}/sim-stdout.txt 2> data3/{wildcards.event_label}/sim-stderr.txt"

rule reconstruction:
    input:
        "data3/{event_label}/particles.root",
        "data3/{event_label}/hits.root"
    output:
        "data3/{event_label}/reco/tracksummary_ckf.root"
    shell:
        "python scripts/reconstruction.py data3/ {wildcards.event_label} > data3/{wildcards.event_label}/reco/reco-stdout.txt 2> data3/{wildcards.event_label}/reco/reco-stderr.txt"
