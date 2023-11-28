"""
originally from https://github.com/paulgessinger
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse

from mycommon.plot_style import myPlotStyle


def line(p1, p2, ax, **kwargs):
    v = np.stack([p1, p2])
    ax.plot(*v.T, **kwargs)


def box(bl, tr, ax, **kwargs):
    bl = np.array(bl)
    tr = np.array(tr)

    c = (bl + tr) / 2
    w = tr[0] - bl[0]
    h = tr[1] - bl[1]

    rect = patches.Rectangle(bl, w, h, **kwargs)
    ax.add_patch(rect)


def draw_eta_lines(
    ax,
    eta_range=(-3, 3),
    n=None,
    s=None,
    fmt="%.0f",
    text_args={},
    rmin=None,
    rmax=None,
    zmin=None,
    zmax=None,
    **kwargs
):
    assert (n is None and s is not None) or (n is not None and s is None)
    eta_min, eta_max = eta_range
    if n is not None:
        eta_vals = np.linspace(eta_min, eta_max, n)
    else:
        eta_vals = np.arange(eta_min, eta_max + s, s)
    thetas = 2 * np.arctan(np.exp(-eta_vals))
    if zmin is None:
        zmin, zmax = ax.get_xlim()
    if rmin is None:
        rmin, rmax = ax.get_ylim()
    z_vals = np.array([zmin if theta > np.pi / 2.0 else zmax for theta in thetas])
    r_vals = z_vals * np.tan(thetas)

    for eta, theta, z_out, r_out in zip(eta_vals, thetas, z_vals, r_vals):
        z_in = 0
        r_in = 0

        ha = "right"
        va = "center"

        kwargs2 = kwargs.copy()
        if eta % 1 == 0:
            kwargs2["color"] = "grey"

        if eta == 0.0:
            ax.text(0, rmax, s=r"$\eta=0$", ha="center", va="bottom", **text_args)
            ax.plot([0, 0], [max(r_in, rmin), rmax], **kwargs2)
            continue

        if r_out > rmax:
            # re-calc z based on r
            z_out = rmax / np.tan(theta)
            r_out = rmax
            va = "bottom"

        if r_out < rmin:
            # would not show anyway
            continue

        if r_in < rmin:
            # re-calc z_in
            z_in = rmin / np.tan(theta)
            r_in = rmin

        ax.plot([z_in, z_out], [r_in, r_out], **kwargs2)
        if eta > 0:
            ha = "left"

        if eta > 0 and eta % 1 == 0:
            ax.text(
                z_out, r_out, s=r"$\eta=%s$" % (fmt % eta), ha=ha, va=va, **text_args
            )


fig = myPlotStyle(figsize=(12, 7))

parser = argparse.ArgumentParser()
parser.add_argument("--output")
args = parser.parse_args()

# set margin on the right side of the plot
fig.tight_layout()
fig.subplots_adjust(right=0.82)

ax = fig.subplots(1, 1)

ax.set_xlim(-3200, 3800)
ax.set_ylim(-20, 1350)

draw_eta_lines(
    ax=ax,
    eta_range=(-3, 3),
    s=0.25,
    color="lightgray",
    rmin=35,
    rmax=1250,
    zmin=-3200,
    zmax=3200,
    linestyle="--",
)


# beampipe
bp_rmin = 23.6
bp_rmax = 24.4
bp_r = (bp_rmin + bp_rmax) / 2
bp_length = 4000

box([-bp_length, bp_rmin], [bp_length, bp_rmax], ax=ax, color="black")
# line([-bp_length, bp_r], [bp_length, bp_r], ax=ax, color="black")


box([-1600, 28], [1600, 200], ax=ax, color="tab:blue", alpha=0.2)

# pixel barrel
module_width = 72
module_gap = 0.5
module_num = 14
dz = module_width * module_num + module_gap * (module_num - 1) / 2
for r in [34, 70, 116, 172]:
    line([-dz / 2, r], [dz / 2, r], ax=ax, color="tab:blue")

# pixel endcap
rmin = 28
rmax = 186
for z in [620, 720, 840, 980, 1120, 1320, 1520]:
    for s in -1, 1:
        #     for r, dr, dz in [
        #         (76, 34, 3.5),
        #         (144, 34, -3.5)
        #     ]:
        #         line([s*(z-dz), r-dr], [s*(z-dz), r+dr], ax=ax, color="tab:blue")
        line([s * z, rmin], [s * z, rmax], ax=ax, color="tab:blue")


box([-3060, 205], [3060, 716], ax=ax, color="tab:red", alpha=0.2)

# sstrip barrel
module_width = 108
module_gap = 0.5
module_num = 21
dz = module_width * module_num + module_gap * (module_num - 1) / 2
for r in [260, 360, 500, 660]:
    line([-dz / 2, r], [dz / 2, r], ax=ax, color="tab:red")

# sstrip endcap
rmin = 210
rmax = 715
for z in [1300, 1550, 1850, 2200, 2550, 2950]:
    for s in -1, 1:
        #     for r, dr, dz in [
        #         (318, 78, 5),
        #         (470, 78, -5),
        #         (622, 78, 5)
        #     ]:
        #         line([s*(z-dz), r-dr], [s*(z-dz), r+dr], ax=ax, color="tab:red")
        line([s * z, rmin], [s * z, rmax], ax=ax, color="tab:red")

box([-3060, 719], [3060, 1105], ax=ax, color="tab:green", alpha=0.2)

# lstrip barrel
module_width = 108
module_gap = 0.5
module_num = 21
dz = module_width * module_num + module_gap * (module_num - 1) / 2
for r in [820, 1020]:
    line([-dz / 2, r], [dz / 2, r], ax=ax, color="tab:green")

# lstrip endcap
rmin = 720
rmax = 1095
for z in [1300, 1600, 1900, 2250, 2600, 3000]:
    for s in -1, 1:
        #     for r, dr, dz in [
        #         (820, 78, 15),
        #         (990, 78, -15),
        #     ]:
        # line([s*(z-dz), r-dr], [s*(z-dz), r+dr], ax=ax, color="tab:green")
        line([s * z, rmin], [s * z, rmax], ax=ax, color="tab:green")

box([-3000, 1160], [3000, 1200], ax=ax, color="tab:purple", alpha=0.7)

ax.text(4700, 15, s="Beam pipe", va="top", ha="center")
ax.text(4700, 80, s="Pixels", ha="center")
ax.text(4700, 420, s="Short Strips", ha="center")
ax.text(4700, 880, s="Long Strips", ha="center")
ax.text(4700, 1200, s="Solenoid", va="bottom", ha="center")

ax.annotate(
    "",
    xy=(3800, 25),
    xytext=(5500, 25),
    arrowprops=dict(
        alpha=0.7, linewidth=1.7, linestyle="--", arrowstyle="-", color="black"
    ),
)
ax.annotate(
    "",
    xy=(3800, 186),
    xytext=(5500, 186),
    arrowprops=dict(
        alpha=0.7, linewidth=1.7, linestyle="--", arrowstyle="-", color="black"
    ),
)
ax.annotate(
    "",
    xy=(3800, 715),
    xytext=(5500, 715),
    arrowprops=dict(
        alpha=0.7, linewidth=1.7, linestyle="--", arrowstyle="-", color="black"
    ),
)
ax.annotate(
    "",
    xy=(3800, 1095),
    xytext=(5500, 1095),
    arrowprops=dict(
        alpha=0.7, linewidth=1.7, linestyle="--", arrowstyle="-", color="black"
    ),
)

ax.set_ylabel("r [mm]")
ax.set_xlabel("z [mm]")

if args.output:
    fig.savefig(args.output)
else:
    plt.show()
