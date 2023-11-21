import argparse

import pandas as pd
import ROOT as r

from mycommon.root import getDefaultStyle, createLabel, createLegend, createEfficency


style = getDefaultStyle()

parser = argparse.ArgumentParser()
parser.add_argument("mu_1GeV", help="Path to the 1 GeV muon efficiency csv file")
parser.add_argument("mu_10GeV", help="Path to the 10 GeV muon efficiency csv file")
parser.add_argument("mu_100GeV", help="Path to the 100 GeV muon efficiency csv file")
parser.add_argument("output", help="Path to the output pdf file")
args = parser.parse_args()

files = [args.mu_1GeV, args.mu_10GeV, args.mu_100GeV]
data = [pd.read_csv(file) for file in files]
graphs = [createEfficency(d, s) for d, s in zip(data, style)]

canvas = r.TCanvas("canvas", "", 800, 600)

r.gPad.SetLeftMargin(0.12)
r.gPad.SetRightMargin(0.05)
r.gPad.SetTopMargin(0.05)

graphs[0].GetYaxis().SetRangeUser(0.96, 1.015)
graphs[0].Draw("AP")
for graph in graphs[1:]:
    graph.Draw("Psame")

line = r.TLine(0.0, 1.0, 3.0, 1.0)
line.SetLineStyle(2)
line.Draw()

text = createLabel(event_string="single muons, <#mu>=0", is_ttbar=False)

legend = createLegend()
legend.AddEntry(graphs[0], "p_{T}=1 GeV", "pl")
legend.AddEntry(graphs[1], "p_{T}=10 GeV", "pl")
legend.AddEntry(graphs[2], "p_{T}=100 GeV", "pl")
legend.Draw()

canvas.SaveAs(args.output)
