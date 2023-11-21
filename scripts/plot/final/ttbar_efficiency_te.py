import argparse

import pandas as pd
import ROOT as r

from mycommon.root import getDefaultStyle, createLabel, createLegend, createEfficency


style = getDefaultStyle()

parser = argparse.ArgumentParser()
parser.add_argument("ttbar60", help="Path to the ttbar 60 efficiency csv file")
parser.add_argument("ttbar120", help="Path to the ttbar 120 efficiency csv file")
parser.add_argument("ttbar200", help="Path to the ttbar 200 efficiency csv file")
parser.add_argument("output", help="Path to the output pdf file")
args = parser.parse_args()

files = [args.ttbar60, args.ttbar120, args.ttbar200]
data = [pd.read_csv(file) for file in files]
graphs = [createEfficency(d, s) for d, s in zip(data, style)]

canvas = r.TCanvas("canvas", "", 800, 600)

r.gPad.SetLeftMargin(0.12)
r.gPad.SetRightMargin(0.05)
r.gPad.SetTopMargin(0.05)

graphs[0].GetYaxis().SetRangeUser(0.55, 1.15)
graphs[0].Draw("AP")
for graph in graphs[1:]:
    graph.Draw("Psame")

line = r.TLine(0.0, 1.0, 3.0, 1.0)
line.SetLineStyle(2)
line.Draw()

text = createLabel(event_string=r"t#bar{t}", is_ttbar=True)

legend = createLegend()
legend.AddEntry(graphs[0], "<#mu>=60", "pl")
legend.AddEntry(graphs[1], "<#mu>=120", "pl")
legend.AddEntry(graphs[2], "<#mu>=200", "pl")
legend.Draw()

canvas.SaveAs(args.output)
