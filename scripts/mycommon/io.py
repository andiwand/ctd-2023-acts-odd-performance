import uproot
import awkward as ak
import pandas as pd


def read_track_efficiency(file):
    if str(file).endswith(".csv"):
        data = pd.read_csv(file)

        eta = data["true_eta"].values
        track_efficiency = data["track_efficiency"].values

        return eta, track_efficiency

    raise ValueError(f"unknown file type: {file}")


def read_pulls(file):
    if str(file).endswith(".root"):
        columns = [
            "pull_eLOC0_fit",
            "pull_eLOC1_fit",
            "pull_eT_fit",
            "pull_ePHI_fit",
            "pull_eTHETA_fit",
            "pull_eQOP_fit",
        ]

        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(
                columns + ["t_eta", "nMeasurements"], library="ak"
            ),
            how="outer",
        ).dropna()
        tracksummary = tracksummary.query("nMeasurements >= 10")

        eta = tracksummary["t_eta"].values
        pulls = [tracksummary[col].values for col in columns]

        return eta, pulls

    if str(file).endswith(".csv"):
        columns = [
            "track_pull_eLOC0_fit",
            "track_pull_eLOC1_fit",
            "track_pull_eT_fit",
            "track_pull_ePHI_fit",
            "track_pull_eTHETA_fit",
            "track_pull_eQOP_fit",
        ]

        data = pd.read_csv(file).dropna()

        eta = data["true_eta"].values
        pulls = [data[col].values for col in columns]

        return eta, pulls

    raise ValueError(f"unknown file type: {file}")


def read_residuals(file):
    if str(file).endswith(".root"):
        tracksummary = uproot.open(file)
        tracksummary = ak.to_dataframe(
            tracksummary["tracksummary"].arrays(
                ["t_pT", "t_eta", "res_eLOC0_fit", "res_eLOC1_fit", "res_eQOP_fit"],
                library="ak",
            ),
            how="outer",
        ).dropna()

        pt = tracksummary["t_pT"].values
        eta = tracksummary["t_eta"].values
        res_d0 = tracksummary["res_eLOC0_fit"].values
        res_z0 = tracksummary["res_eLOC1_fit"].values
        res_qop = tracksummary["res_eQOP_fit"].values

        return {
            "pt": pt,
            "eta": eta,
            "res_d0": res_d0,
            "res_z0": res_z0,
            "res_qop": res_qop,
        }

    if str(file).endswith(".csv"):
        data = pd.read_csv(file).dropna()

        pt = data["true_pt"].values
        eta = data["true_eta"].values
        res_d0 = data["track_res_eLOC0_fit"].values
        res_z0 = data["track_res_eLOC1_fit"].values
        res_qop = data["track_res_eQOP_fit"].values

        return {
            "pt": pt,
            "eta": eta,
            "res_d0": res_d0,
            "res_z0": res_z0,
            "res_qop": res_qop,
        }

    raise ValueError(f"unknown file type: {file}")
