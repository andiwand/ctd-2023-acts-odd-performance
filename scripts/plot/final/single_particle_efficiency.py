import argparse

import pandas as pd
import ROOT as r

from mycommon.root import getDefaultStyle, createLabel, createLegend, createEfficency


style = getDefaultStyle()

parser = argparse.ArgumentParser()
parser.add_argument("mu", help="Path to the muon efficiency csv file")
parser.add_argument("pi", help="Path to the pion efficiency csv file")
parser.add_argument("e", help="Path to the electron efficiency csv file")
parser.add_argument("output", help="Path to the output pdf file")
args = parser.parse_args()

files = [args.mu, args.pi, args.e]
data = [pd.read_csv(file) for file in files]
graphs = [createEfficency(d, s) for d, s in zip(data, style)]

canvas = r.TCanvas("canvas", "", 800, 600)

r.gPad.SetLeftMargin(0.12)
r.gPad.SetRightMargin(0.05)
r.gPad.SetTopMargin(0.05)

graphs[0].GetYaxis().SetRangeUser(0.75, 1.1)
graphs[0].Draw("AP")
for graph in graphs[1:]:
    graph.Draw("Psame")

line = r.TLine(0.0, 1.0, 3.0, 1.0)
line.SetLineStyle(2)
line.Draw()

text = createLabel(event_string="single particles, <#mu>=0", is_ttbar=False)

legend = createLegend()
legend.AddEntry(graphs[0], "muons", "pl")
legend.AddEntry(graphs[1], "pions", "pl")
legend.AddEntry(graphs[2], "electrons", "pl")
legend.Draw()

canvas.SaveAs(args.output)
