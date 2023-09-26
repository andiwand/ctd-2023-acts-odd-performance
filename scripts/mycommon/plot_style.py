import matplotlib
import matplotlib.pyplot as plt


def myPlotStyle():
    font = {
        "size": 18,
    }

    matplotlib.rc("font", **font)

    plt.figure(figsize=(12, 8), dpi=80)
