#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("tracksummary", nargs="+")
args = parser.parse_args()

columns = [
    "pull_eLOC0_fit",
    "pull_eLOC1_fit",
    "pull_eT_fit",
    "pull_ePHI_fit",
    "pull_eTHETA_fit",
    "pull_eQOP_fit",
]
eta_range = (-3, 3)
pull_range = (-4, 4)
std_range = (0, 4)

for f in args.tracksummary:
    tracksummary = uproot.open(f)
    tracksummary = ak.to_dataframe(
        tracksummary["tracksummary"].arrays(columns + ["t_eta"], library="ak"),
        how="outer",
    )
    tracksummary = tracksummary.dropna()

    fig = plt.figure(f, figsize=(8, 6))
    subfigs = fig.subfigures(2, 3)

    for col, subfig in zip(
        columns,
        subfigs.flat,
    ):
        ax, ax_mean, ax_std = subfig.subplots(
            3,
            1,
            gridspec_kw={"height_ratios": [5, 1, 1], "hspace": 0},
            sharex=True,
        )

        ax.set_title(col)
        h, xedges, yedges, im = ax.hist2d(
            tracksummary["t_eta"],
            tracksummary[col],
            range=(eta_range, pull_range),
            bins=(20, 20),
            density=True,
            cmap="Oranges",
        )

        xmiddles = 0.5 * (xedges[:-1] + xedges[1:])
        ymiddles = 0.5 * (yedges[:-1] + yedges[1:])
        ydensity = h / np.sum(h, axis=1)[:, np.newaxis]

        ymean = np.sum(ymiddles * ydensity, axis=1)
        yvar = np.sum((ymiddles - ymean[:, np.newaxis]) ** 2 * ydensity, axis=1)
        ystd = np.sqrt(yvar)

        ax.plot(xmiddles, ymean, linestyle="-", color="black", label="fit")
        ax.plot(xmiddles, ymean - ystd, linestyle="--", color="black")
        ax.plot(xmiddles, ymean + ystd, linestyle="--", color="black")

        ax_mean.plot(xmiddles, ymean, linestyle="-", color="black")
        ax_mean.axhline(0, linestyle="--", color="gray")

        ax_std.plot(xmiddles, ystd, linestyle="-", color="black")
        ax_std.axhline(1, linestyle="--", color="gray")

        fig.colorbar(im, ax=subfig.get_axes())

plt.show()
