#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import mplhep
from hist import Hist
import argparse

from mycommon.plot_style import myPlotStyle


plt.rcParams["ytick.right"] = plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["font.size"] = 12.0
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.columnspacing"] = 0.2
plt.rcParams["legend.handletextpad"] = 0.2
plt.rcParams["legend.labelspacing"] = 0.2
plt.rcParams["legend.borderpad"] = 0
plt.rcParams["legend.handlelength"] = 1.0


myPlotStyle()

parser = argparse.ArgumentParser(description="Make material composition plots")
parser.add_argument("x", choices=["phi", "eta"])
parser.add_argument("y", choices=["l0", "x0"])
parser.add_argument("input", help="Input root file with histograms")
parser.add_argument("--output")
args = parser.parse_args()

names = {
    "all": "Full detector",
    "beampipe": "Beam pipe",
    "sstrips": "Short Strips",
    "lstrips": "Long Strips",
    "pixel": "Pixel",
    "solenoid": "Solenoid",
    "ecal": "EM Calorimeter",
}
region_order = ["beampipe", "pixel", "sstrips", "lstrips", "solenoid"]

x_label = {"phi": r"$\phi$", "eta": "$\eta$"}[args.x]
y_label = {"l0": r"$\lambda_0$", "x0": "$X_0$"}[args.y]

hists = []
labels = []

rf = uproot.open(args.input)

all_hist = rf[f"all_{args.y}_vs_{args.x}_all"].to_hist()

for name in region_order:
    key = f"{name}_{args.y}_vs_{args.x}_all"
    if key in rf:
        o = rf[key].to_hist()
        o.axes[0].label = args.x
        hists.append(o)
    else:
        hists.append(Hist(all_hist.axes[0]))
    labels.append(name)

ax = plt.gcf().subplots()
mplhep.histplot(hists, ax=ax, stack=True, histtype="fill", label=labels)
ymin, ymax = ax.get_ylim()
ax.set_xlim(hists[0].axes[0].edges[0], hists[0].axes[0].edges[-1])
ax.set_ylim(top=1.2 * ymax)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.legend(ncol=3)

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
