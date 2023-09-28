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
from mycommon.label import get_event_type_label
from mycommon.paths import get_event_label_from_path


def get_data(file):
    if str(file).endswith(".root"):
        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(
                ["t_pT", "t_charge", "t_p", "eQOP_fit"], library="ak"
            ),
            how="outer",
        ).dropna()

        true_pt = tracksummary["t_pT"].values
        true_q = tracksummary["t_charge"].values
        true_p = tracksummary["t_p"].values
        track_qop = tracksummary["eQOP_fit"].values
        res_p = true_q / track_qop - true_p

        return true_pt, res_p

    if str(file).endswith(".csv"):
        data = pd.read_csv(file).dropna()

        true_pt = data["true_pt"].values
        true_q = data["true_q"].values
        true_p = data["true_p"].values
        track_qop = data["track_eQOP_fit"].values
        res_p = true_q / track_qop - true_p

        return true_pt, res_p

    raise ValueError(f"unknown file type: {file}")


myPlotStyle()

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+")
parser.add_argument("--output")
args = parser.parse_args()

pt_range = (1, 100)
pt_bins = 10

for file in args.input:
    event_label = get_event_label_from_path(file)
    event, _ = split_event_label(event_label)

    true_pt, res_p = get_data(file)

    p_std, pt_edges, _ = binned_statistic(
        true_pt, res_p, bins=pt_bins, range=pt_range, statistic=smoothed_std
    )
    p_std_std, pt_edges, _ = binned_statistic(
        true_pt, res_p, bins=pt_bins, range=pt_range, statistic=smoothed_std_std
    )
    pt_mid = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    pt_step = pt_edges[1] - pt_edges[0]

    plt.errorbar(
        x=pt_mid,
        y=p_std,
        yerr=p_std_std,
        xerr=pt_step * 0.4,
        fmt="",
        linestyle="",
        label=get_event_type_label(event),
    )

plt.title(rf"Resolution of $p$ over $p_T$")
plt.xlabel(r"$p_T$")
plt.ylabel(r"$\sigma(p)$ [GeV]")
plt.xticks(np.linspace(0, pt_range[1], 11))
plt.xlim(pt_range)
plt.legend()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
