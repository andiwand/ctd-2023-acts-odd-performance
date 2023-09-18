from pathlib import Path


def pt_label(file):
    split = Path(file).parent.name.split("_")
    return f"{split[1].replace('GeV', '')} GeV"


def particle_pt_label(file):
    particle_map = {
        "mu": "$\mu$",
        "pi": "$\pi$",
        "e": "$e$",
    }
    split = Path(file).parent.name.split("_")
    return f"single {particle_map[split[0]]} {split[1].replace('GeV', '')} GeV"


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
