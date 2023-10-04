import matplotlib.pyplot as plt


def myPlotStyle():
    plt.rcParams["ytick.right"] = plt.rcParams["xtick.top"] = True
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["font.size"] = 18.0
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["legend.columnspacing"] = 0.2
    plt.rcParams["legend.handletextpad"] = 0.2
    plt.rcParams["legend.labelspacing"] = 0.2
    plt.rcParams["legend.borderpad"] = 0
    plt.rcParams["legend.handlelength"] = 1.0

    fig = plt.figure(figsize=(16, 9), dpi=100)
    return fig
