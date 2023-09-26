import numpy as np
import scipy.stats


def smoothed_mean(data):
    return smoothed_mean_std(data)[0]


def smoothed_std(data):
    return smoothed_mean_std(data)[1]


def smoothed_mean_std(data):
    def fit2(data):
        return np.mean(data), np.std(data)

    for _ in range(3):
        m, s = fit2(data)
        data = data[np.abs(data - np.median(data)) < 3 * s]
    return m, s


def clopper_pearson(k, n, alpha=0.32):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """
    p = k / n
    if p == 1.0:
        return p, p, p
    p_upper = np.maximum(
        scipy.stats.beta.ppf(1 - alpha / 2, k + 1, n - k), np.zeros_like(p)
    )
    p_lower = np.minimum(scipy.stats.beta.ppf(alpha / 2, k, n - k + 1), np.ones_like(p))
    return p, p_upper, p_lower


def create_clopper_pearson_upper_bounds(alpha=0.32):
    def interval(data):
        _, p_upper, _ = clopper_pearson(np.sum(data), len(data), alpha)
        return p_upper

    return interval


def create_clopper_pearson_lower_bounds(alpha=0.32):
    def interval(data):
        _, _, p_lower = clopper_pearson(np.sum(data), len(data), alpha)
        return p_lower

    return interval
