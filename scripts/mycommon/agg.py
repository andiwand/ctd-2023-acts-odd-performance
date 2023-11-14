from scipy.stats import binned_statistic

from mycommon.stats import (
    create_clopper_pearson_upper_bounds,
    create_clopper_pearson_lower_bounds,
    robust_mean,
    robust_std,
    robust_std_std,
)


def agg_efficiency_over_eta(eta_range, eta_bins, eta, track_efficiency):
    track_efficiency_mean, eta_edges, _ = binned_statistic(
        eta,
        track_efficiency,
        bins=eta_bins,
        range=eta_range,
        statistic="mean",
    )
    track_efficiency_upper, _, _ = binned_statistic(
        eta,
        track_efficiency,
        bins=eta_bins,
        range=eta_range,
        statistic=create_clopper_pearson_upper_bounds(),
    )
    track_efficiency_lower, _, _ = binned_statistic(
        eta,
        track_efficiency,
        bins=eta_bins,
        range=eta_range,
        statistic=create_clopper_pearson_lower_bounds(),
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    return (
        eta_mid,
        (
            track_efficiency_mean,
            (
                track_efficiency_mean - track_efficiency_lower,
                track_efficiency_upper - track_efficiency_mean,
            ),
        ),
    )


def agg_pulls_over_eta(eta_range, eta_bins, eta, pull):
    mean_binned, eta_edges, _ = binned_statistic(
        eta,
        pull,
        bins=eta_bins,
        range=eta_range,
        statistic=robust_mean,
    )
    std_binned, _, _ = binned_statistic(
        eta,
        pull,
        bins=eta_bins,
        range=eta_range,
        statistic=robust_std,
    )
    eta_mid = 0.5 * (eta_edges[:-1] + eta_edges[1:])

    return (
        eta_mid,
        (
            mean_binned,
            std_binned,
        ),
    )


def agg_resolution(x_range, x_bins, x, y):
    std, x_edges, _ = binned_statistic(
        x,
        y,
        bins=x_bins,
        range=x_range,
        statistic=robust_std,
    )
    std_std, _, _ = binned_statistic(
        x,
        y,
        bins=x_bins,
        range=x_range,
        statistic=robust_std_std,
    )
    x_mid = 0.5 * (x_edges[:-1] + x_edges[1:])

    return (
        x_mid,
        (
            std,
            std_std,
        ),
    )
