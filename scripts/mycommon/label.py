from pathlib import Path


def split_path(file):
    split = Path(file).parent.parent.name.split("_")
    return split[0], split[1]


def pt_label(file):
    split = split_path(file)
    return f"{split[1].replace('GeV', '')} GeV"


def particle_label(file):
    particle_map = {
        "mu": "$\mu$",
        "pi": "$\pi$",
        "e": "$e$",
    }
    split = split_path(file)
    return f"{particle_map[split[0]]}"


def particle_pt_label(file):
    return f"single {particle_label(file)} {pt_label(file)} GeV"


def param_label(param):
    param_map = {
        "pull_eLOC0_fit": r"$d_0$",
        "pull_eLOC1_fit": r"$z_0$",
        "pull_eT_fit": r"$t$",
        "pull_ePHI_fit": r"$\phi$",
        "pull_eTHETA_fit": r"$\theta$",
        "pull_eQOP_fit": r"$\frac{q}{p}$",
    }
    return param_map[param]
