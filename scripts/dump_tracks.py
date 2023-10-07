import argparse
import uproot
import awkward as ak


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("event_nr", type=int)
parser.add_argument("output")
args = parser.parse_args()

trackstates = ak.to_dataframe(
    uproot.open(args.input)["trackstates"].arrays(
        ["event_nr", "multiTraj_nr", "subTraj_nr", "g_x_smt", "g_y_smt", "g_z_smt"],
        library="ak",
    ),
    how="outer",
).dropna()
trackstates = trackstates[trackstates["event_nr"] == int(args.event_nr)]
trackstates[
    ["event_nr", "multiTraj_nr", "subTraj_nr", "g_x_smt", "g_y_smt", "g_z_smt"]
].to_csv(args.output, index=False)
