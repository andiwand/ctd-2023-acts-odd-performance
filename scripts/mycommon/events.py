import itertools


events = [
    f"{p}_{pT}GeV" for p, pT in itertools.product(["mu", "pi", "e"], [1, 10, 100])
] + [f"ttbar_{pu}" for pu in [0, 60, 120, 200]]
simulations = [
    "fatras",
    "geant4",
]


def list_events_simulations():
    return list(itertools.product(events, simulations))


def list_event_labels():
    return [
        create_event_label(event, simulation)
        for event, simulation in list_events_simulations()
    ]


def create_event_label(event, simulation):
    return f"{event}_{simulation}"


def split_event_label(event_label):
    for event, simulation in list_events_simulations():
        if event_label == create_event_label(event, simulation):
            return event, simulation
    raise ValueError(f"unknown event label {event_label}")


def get_event_type(event):
    if event.startswith("ttbar_"):
        return "ttbar"
    return "single_particles"


def get_number_of_events(event):
    event_type = get_event_type(event)
    return 200000 if event_type == "single_particles" else 1
