import itertools


single_particles = ["mu", "pi", "e"]
single_particle_pts = [1, 10, 100]
single_particle_pt_ranges = ["1-100GeV"]
ttbar_pileups = [0, 60, 120, 200]
simulations = ["fatras", "geant4"]


def list_single_particles():
    return single_particles


def list_single_particle_pt_labels():
    return [f"{pt}GeV" for pt in single_particle_pts]


def list_single_particle_pt_range_labels():
    return single_particle_pt_ranges


def list_ttbar_pileups():
    return ttbar_pileups


def list_events():
    return (
        [
            f"{p}_{pT}GeV"
            for p, pT in itertools.product(single_particles, single_particle_pts)
        ]
        + [
            f"{p}_{pT_range}"
            for p, pT_range in itertools.product(
                single_particles, single_particle_pt_ranges
            )
        ]
        + [f"ttbar_{pu}" for pu in ttbar_pileups]
    )


def list_simulations():
    return simulations


def list_events_simulations():
    return list(itertools.product(list_events(), list_simulations()))


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
    split = event.split("_")
    if split[0] == "ttbar":
        return "ttbar"
    if split[0] in single_particles:
        return "single_particles"
    raise ValueError(f"cannot determine event type: {event}")


def get_event_details(event):
    event_type = get_event_type(event)
    split = event.split("_")
    if event_type == "ttbar":
        pileup = int(split[1])
        return pileup
    if event_type == "single_particles":
        particle, pt = split[0], split[1]
        pt = pt.replace("GeV", "")
        if "-" in pt:
            pt_range = list(map(int, pt.split("-")))
            return particle, (pt_range[0], pt_range[1])
        pt = int(pt)
        return particle, pt
    raise ValueError(f"unknown event type: {event_type}")
