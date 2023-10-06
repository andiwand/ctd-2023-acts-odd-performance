import uproot
import awkward as ak
import pandas as pd


def get_pull_data(file):
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
