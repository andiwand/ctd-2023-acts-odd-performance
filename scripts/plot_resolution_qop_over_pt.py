#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse
import numpy as np
import pandas as pd
from scipy.stats import binned_statistic

from mycommon.plot_style import myPlotStyle
from mycommon.stats import smoothed_std
from mycommon.events import split_event_label
from mycommon.label import get_event_variant_label, get_event_type_label
from mycommon.paths import get_event_label_from_path


def get_data(file):
    if str(file).endswith(".root"):
        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(["t_pT", "res_eQOP_fit"], library="ak"),
            how="outer",
        ).dropna()

        pt = tracksummary["t_pT"].values
        res_qop = tracksummary["res_eQOP_fit"].values

        return pt, res_qop

    if str(file).endswith(".csv"):
        data = pd.read_csv(file).dropna()

        pt = data["true_pt"].values
        res_qop = data["track_res_eQOP_fit"].values

        return pt, res_qop

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

pt_range = (0, 100)
pt_bins = 31

for file in args.input:
    event_label = get_event_label_from_path(file)
    event, _ = split_event_label(event_label)

    eta, res_qop = get_data(file)

    resolution_qop_binned, pt_edges, _ = binned_statistic(
        eta, res_qop, bins=pt_bins, range=pt_range, statistic=smoothed_std
    )
    pt_mid = 0.5 * (pt_edges[:-1] + pt_edges[1:])

    plt.plot(
        pt_mid,
        resolution_qop_binned,
        marker="o",
        linestyle="",
        label=get_event_variant_label(event),
    )

plt.title(
    rf"Resolution of $\frac{{q}}{{p}}$ over $p_T$ for {get_event_type_label(event)} events"
)
plt.xlabel(r"$p_T$")
plt.ylabel(r"$\sigma(\frac{{q}}{{p}})$ [mm]")
plt.xticks(np.linspace(*pt_range, 11))
plt.xlim(pt_range)
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
