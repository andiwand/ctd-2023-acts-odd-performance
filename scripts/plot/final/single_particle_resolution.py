import argparse

import pandas as pd
import ROOT as r

from mycommon.root import getDefaultStyle, createLabel, createLegend, createResolution


style = getDefaultStyle()

parser = argparse.ArgumentParser()
parser.add_argument("mu", help="Path to the muon z0 resolution csv file")
parser.add_argument("pi", help="Path to the pion z0 resolution csv file")
parser.add_argument("e", help="Path to the electron z0 resolution csv file")
parser.add_argument("output", help="Path to the output pdf file")
args = parser.parse_args()

files = [args.mu, args.pi, args.e]
data = [pd.read_csv(file) for file in files]
graphs = [createResolution(d, s) for d, s in zip(data, style)]

canvas = r.TCanvas("canvas", "", 800, 600)

r.gPad.SetLeftMargin(0.12)
r.gPad.SetRightMargin(0.05)
r.gPad.SetTopMargin(0.05)

graphs[0].GetXaxis().SetTitle("p_{T} [GeV]")
graphs[0].GetYaxis().SetTitle("#sigma_{z_{0}} [mm]")
graphs[0].GetXaxis().SetRangeUser(0.0, 100.0)
graphs[0].GetYaxis().SetRangeUser(0.0, 0.032)

graphs[0].Draw("AP")
for graph in graphs[1:]:
    graph.Draw("Psame")

#line = r.TLine(0.0, 0.015, 100.0, 0.015)
#line.SetLineStyle(2)
#line.Draw()

text = createLabel(event_string="single particles, <#mu>=0", is_ttbar=False)

legend = createLegend()
legend.AddEntry(graphs[0], "muon", "pl")
legend.AddEntry(graphs[1], "pion", "pl")
legend.AddEntry(graphs[2], "electron", "pl")
legend.Draw()

canvas.SaveAs(args.output)
