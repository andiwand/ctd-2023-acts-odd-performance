import matplotlib
import matplotlib.pyplot as plt


def myPlotStyle():
    font = {
        "size": 18,
    }

    matplotlib.rc("font", **font)

    plt.figure(figsize=(16, 9), dpi=100)
