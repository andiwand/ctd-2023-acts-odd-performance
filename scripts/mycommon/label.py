from mycommon.events import get_event_type, get_event_details


def get_event_variant_label(event):
    event_type = get_event_type(event)
    event_details = get_event_details(event)
    if event_type == "ttbar":
        return f"{event_details} pileup"
    if event_type == "single_particles":
        return f"{event_details[1]} GeV"
    raise ValueError(f"unknown event type {event_type}")


def get_single_particle_label(event):
    particle_map = {
        "mu": r"$\mu$",
        "pi": r"$\pi$",
        "e": r"$e$",
    }
    event_type = get_event_type(event)
    assert event_type == "single_particles", f"event type is {event_type}"
    event_details = get_event_details(event)
    return f"{particle_map[event_details[0]]}"


def get_event_type_label(event):
    event_type = get_event_type(event)
    if event_type == "ttbar":
        return r"$t\bar{t}$"
    return f"single {get_single_particle_label(event)}"
