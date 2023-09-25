from mycommon.events import get_event_type, get_event_details


def get_event_variant_label(event):
    event_type = get_event_type(event)
    event_details = get_event_details(event)
    if event_type == "ttbar":
        return f"{event_details[0]} pileup"
    return f"{event_details[1]} GeV"


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


def get_param_label(param):
    param_map = {
        "pull_eLOC0_fit": r"$d_0$",
        "pull_eLOC1_fit": r"$z_0$",
        "pull_eT_fit": r"$t$",
        "pull_ePHI_fit": r"$\phi$",
        "pull_eTHETA_fit": r"$\theta$",
        "pull_eQOP_fit": r"$\frac{q}{p}$",
    }
    return param_map[param]
