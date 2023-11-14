import argparse
import uproot
import awkward as ak


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("event_nr", type=int)
parser.add_argument("output")
args = parser.parse_args()

hits = uproot.open(args.input)["hits"].arrays(library="pandas").dropna()
hits = ak.to_dataframe(
    uproot.open(args.input)["hits"].arrays(
        ["event_id", "tx", "ty", "tz"],
        library="ak",
    ),
    how="outer",
).dropna()
hits = hits[hits["event_id"] == int(args.event_nr)]
hits[["event_id", "tx", "ty", "tz"]].to_csv(args.output, index=False)
