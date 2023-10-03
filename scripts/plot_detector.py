#!/usr/bin/env python3

import matplotlib.pyplot as plt
import uproot
import awkward as ak
import numpy as np
import argparse


def etaToSlope(eta):
    return np.tan(2 * np.arctan(np.exp(-eta)))


parser = argparse.ArgumentParser()
parser.add_argument("hits")
parser.add_argument("--more-hits", nargs="+")
args = parser.parse_args()

hits = uproot.open(args.hits)
hits = ak.to_dataframe(hits["hits"].arrays(library="ak"))

fig_xy = plt.figure("XY")
fig_rz = plt.figure("RZ")

mask_xy = (hits["tz"] > -550) & (hits["tz"] < 550)

fig_xy.gca().scatter(hits[mask_xy]["tx"], hits[mask_xy]["ty"], s=1, c="k")
fig_rz.gca().scatter(hits["tz"], np.hypot(hits["tx"], hits["ty"]), s=1, c="k")

for file in args.more_hits:
    morehits = uproot.open(file)
    morehits = ak.to_dataframe(morehits["hits"].arrays(library="ak"))

    fig_xy.gca().scatter(morehits["tx"], morehits["ty"], s=20, label=file)
    fig_rz.gca().scatter(morehits["tz"], np.hypot(morehits["tx"], morehits["ty"]), s=20, label=file)

for eta in [0, 1, 2, 3, 4]:
    slope = etaToSlope(eta)
    fig_rz.gca().axline((0, 0), slope=slope, color="red", linestyle="--")
    fig_rz.gca().axline((0, 0), slope=-slope, color="red", linestyle="--")
fig_rz.gca().set_ylim((0, None))

fig_xy.gca().legend()
fig_rz.gca().legend()

plt.show()
