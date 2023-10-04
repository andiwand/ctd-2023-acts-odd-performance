#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse
import numpy as np
import pandas as pd
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.stats import smoothed_std, smoothed_std_std
from mycommon.events import split_event_label
from mycommon.label import get_event_variant_label, get_event_type_label
from mycommon.paths import get_event_label_from_path


def get_data(file):
    if str(file).endswith(".root"):
        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(
                ["t_pT", "t_eta", "res_eLOC0_fit", "res_eLOC1_fit", "res_eQOP_fit"],
                library="ak",
            ),
            how="outer",
        ).dropna()

        pt = tracksummary["t_pT"].values
        eta = tracksummary["t_eta"].values
        res_d0 = tracksummary["res_eLOC0_fit"].values
        res_z0 = tracksummary["res_eLOC1_fit"].values
        res_qop = tracksummary["res_eQOP_fit"].values

        return {
            "pt": pt,
            "eta": eta,
            "res_d0": res_d0,
            "res_z0": res_z0,
            "res_qop": res_qop,
        }

    if str(file).endswith(".csv"):
        data = pd.read_csv(file).dropna()

        pt = data["true_pt"].values
        eta = data["true_eta"].values
        res_d0 = data["track_res_eLOC0_fit"].values
        res_z0 = data["track_res_eLOC1_fit"].values
        res_qop = data["track_res_eQOP_fit"].values

        return {
            "pt": pt,
            "eta": eta,
            "res_d0": res_d0,
            "res_z0": res_z0,
            "res_qop": res_qop,
        }

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("x", choices=["eta", "pt"])
parser.add_argument("y", choices=["d0", "z0", "qop"])
parser.add_argument("input", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

x_range = {
    "eta": (-3, 3),
    "pt": (1, 100),
}[args.x]
x_bins = {
    "eta": 25,
    "pt": 10,
}[args.x]
x_label = {
    "eta": r"$\eta$",
    "pt": r"$p_T$",
}[args.x]
x_unit = {
    "eta": r"",
    "pt": r" [GeV]",
}[args.x]

y_label = {
    "d0": r"$d_0$",
    "z0": r"$z_0$",
    "qop": r"$\frac{q}{p}$",
}[args.y]
y_unit = {
    "d0": r" [mm]",
    "z0": r" [mm]",
    "qop": r" [$\frac{1}{GeV}$]",
}[args.y]

for file in args.input:
    event_label = get_event_label_from_path(file)
    event, _ = split_event_label(event_label)

    data = get_data(file)

    std, x_edges, _ = binned_statistic(
        data[args.x],
        data[f"res_{args.y}"],
        bins=x_bins,
        range=x_range,
        statistic=smoothed_std,
    )
    std_std, _, _ = binned_statistic(
        data[args.x],
        data[f"res_{args.y}"],
        bins=x_bins,
        range=x_range,
        statistic=smoothed_std_std,
    )
    x_mid = 0.5 * (x_edges[:-1] + x_edges[1:])
    x_step = x_mid[1] - x_mid[0]

    plt.errorbar(
        x=x_mid,
        y=std,
        yerr=std_std,
        xerr=x_step * 0.4,
        fmt="",
        linestyle="",
        label=get_event_variant_label(event),
    )

plt.title(
    rf"Resolution of {y_label} over {x_label} for {get_event_type_label(event)} events"
)
plt.xlabel(rf"{x_label}{x_unit}")
plt.ylabel(rf"{y_label}{y_unit}")
if args.x == "eta":
    plt.xticks(np.linspace(*x_range, 7))
    plt.xlim(x_range)
elif args.x == "pt":
    plt.xticks(np.linspace(0, x_range[1], 11))
    plt.xlim(x_range)
else:
    raise ValueError(f"unknown x: {args.x}")
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
