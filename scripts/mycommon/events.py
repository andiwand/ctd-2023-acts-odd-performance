import itertools


single_particles = ["mu", "pi", "e"]
single_particle_pts = [1, 10, 100]
ttbar_pileups = [0, 60, 120, 200]
events = [
    f"{p}_{pT}GeV" for p, pT in itertools.product(single_particles, single_particle_pts)
] + [f"ttbar_{pu}" for pu in ttbar_pileups]
simulations = [
    "fatras",
    "geant4",
]


def list_single_particles():
    return single_particles


def list_single_particle_pt_labels():
    return [f"{pt}GeV" for pt in single_particle_pts]


def list_ttbar_pileups():
    return ttbar_pileups


def list_simulations():
    return simulations


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


def get_number_of_events(event_type):
    if event_type == "single_particles":
        return 200000
    elif event_type == "ttbar":
        return 1000
    raise ValueError(f"Unknown event type: {event_type}")
