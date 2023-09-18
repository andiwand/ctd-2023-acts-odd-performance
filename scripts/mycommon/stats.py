import numpy as np


def smoothed_std(data):
    return smoothed_mean_std(data)[1]


def smoothed_mean(data):
    return smoothed_mean_std(data)[0]


def smoothed_mean_std(data):
    def fit2(data):
        return np.mean(data), np.std(data)

    for _ in range(3):
        m, s = fit2(data)
        data = data[np.abs(data - np.median(data)) < 3 * s]
    return m, s
