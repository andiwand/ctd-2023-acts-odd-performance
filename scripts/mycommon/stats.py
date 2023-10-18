import numpy as np
import scipy.stats
from scipy.optimize import curve_fit


def robust_mean(data):
    (m, s), cov = robust_gauss_fit(data)
    return m


def robust_std(data):
    (m, s), cov = robust_gauss_fit(data)
    return s


def robust_std_std(data):
    (m, s), cov = robust_gauss_fit(data)
    return cov[1, 1] ** 0.5


def robust_gauss_fit_naive(data):
    def fit(data):
        return np.mean(data), np.std(data)

    if len(data) == 0:
        return (0, 0), np.zeros((2, 2))

    for _ in range(3):
        m, s = fit(data)
        data = data[np.abs(data - np.median(data)) < 3 * s]

    return (m, s), np.zeros((2, 2))


def robust_gauss_fit(data):
    def fit(data):
        def gauss(x, m, s):
            return 1 / (s * (2 * np.pi) ** 0.5) * np.exp(-0.5 * ((x - m) / s) ** 2)

        try:
            if len(data) < 20:
                raise ValueError(f"Not enough data to fit a Gaussian: {len(data)}")

            low, mp1s, median, mm1s, high = np.percentile(data, [1, 16, 50, 84, 99])
            hist_range = (low, high)
            bins = 30
            binned, edges = np.histogram(
                data, range=hist_range, bins=bins, density=True
            )
            centers = 0.5 * (edges[1:] + edges[:-1])

            p0 = median, mm1s - mp1s
            params, cov = curve_fit(gauss, centers, binned, p0=p0, maxfev=1000000)

            if cov[1, 1] > 5 * params[1] ** 2:
                raise ValueError(
                    f"Fit failed, covariance too large: {cov[1, 1]} > {5 * params[1] ** 2}"
                )
        except Exception as e:
            print(f"Falling back to naive mean/std. Error: {e}")
            params, cov = (np.mean(data), np.std(data)), np.zeros((2, 2))

        return params, cov

    if len(data) == 0:
        return (0, 0), np.zeros((2, 2))

    for _ in range(3):
        (m, s), cov = fit(data)
        data = data[np.abs(data - m) < 3 * s]

    return (m, s), cov


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
