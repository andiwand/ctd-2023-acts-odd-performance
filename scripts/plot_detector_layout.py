"""
provided by https://github.com/paulgessinger
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse


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
    ax, eta_range=(-3, 3), n=None, s=None, fmt="%.2f", text_args={}, rmin=None, **kwargs
):
    assert (n is None and s is not None) or (n is not None and s is None)
    eta_min, eta_max = eta_range
    if n is not None:
        eta_vals = np.linspace(eta_min, eta_max, n)
    else:
        eta_vals = np.arange(eta_min, eta_max + s, s)
    thetas = 2 * np.arctan(np.exp(-eta_vals))
    zmin, zmax = ax.get_xlim()
    if rmin is None:
        rmin, rmax = ax.get_ylim()
    else:
        _, rmax = ax.get_ylim()
    z_vals = np.array([zmin if theta > np.pi / 2.0 else zmax for theta in thetas])
    r_vals = z_vals * np.tan(thetas)

    for eta, theta, z_out, r_out in zip(eta_vals, thetas, z_vals, r_vals):
        z_in = 0
        r_in = 0

        ha = "right"
        va = "center"

        if eta == 0.0:
            ax.text(0, rmax, s=r"$\eta=0$", ha="center", va="bottom", **text_args)
            ax.plot([0, 0], [max(r_in, rmin), rmax], **kwargs)
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

        ax.plot([z_in, z_out], [r_in, r_out], **kwargs)
        if eta > 0:
            ha = "left"

        ax.text(z_out, r_out, s=r"$%s$" % (fmt % eta), ha=ha, va=va, **text_args)


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


parser = argparse.ArgumentParser()
parser.add_argument("--output")
args = parser.parse_args()

fig, ax = plt.subplots(figsize=(10, 7))

ax.set_xlim(-3200, 3200)
ax.set_ylim(-20, 1250)

draw_eta_lines(
    ax=ax, eta_range=(-3, 3), s=0.5, color="lightgray", rmin=35, linestyle="--"
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

ax.text(0, 23, s="Beam pipe", va="top", ha="center")
ax.text(0, 130, s="Pixels", ha="center")
ax.text(0, 420, s="Short Strips", ha="center")
ax.text(0, 900, s="Long Strips", ha="center")
ax.text(0, 1200, s="Solenoid", va="bottom", ha="center")

ax.set_ylabel("r [mm]")
ax.set_xlabel("z [mm]")

fig.tight_layout()

if args.output:
    fig.savefig(args.output)
else:
    plt.show()
