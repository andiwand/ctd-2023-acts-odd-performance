import argparse

import pandas as pd
import ROOT as r

from mycommon.root import getDefaultStyle, createLabel, createLegend, createPull


style = getDefaultStyle()
labels = [
    "#frac{d_{0}^{reco}-d_{0}^{true}}{#sigma_{d_{0}}^{reco}}",
    "#frac{z_{0}^{reco}-z_{0}^{true}}{#sigma_{z_{0}}^{reco}}",
    "#frac{q/p^{reco}-q/p^{true}}{#sigma_{q/p}^{reco}}",
]

parser = argparse.ArgumentParser()
parser.add_argument("mu_1GeV", help="Path to the 1 GeV muon pull csv file")
parser.add_argument("mu_10GeV", help="Path to the 10 GeV muon pull csv file")
parser.add_argument("mu_100GeV", help="Path to the 100 GeV muon pull csv file")
parser.add_argument("output", help="Path to the output pdf file")
args = parser.parse_args()

files = [args.mu_1GeV, args.mu_10GeV, args.mu_100GeV]
data = [pd.read_csv(file) for file in files]

canvas = r.TCanvas("canvas", "", 1000, 600)

pads = [
    r.TPad("a", "a", 0.0, 0.6, 0.7, 0.95),
    r.TPad("b", "b", 0.0, 0.3, 0.7, 0.65),
    r.TPad("c", "c", 0.0, 0.0, 0.7, 0.35),
]

all_graphs = []
all_lines = []

for i, (pad, d, title) in enumerate(zip(pads, data, labels)):
    canvas.cd()
    pad.Draw()
    pad.cd()

    pad.SetTopMargin(0.0)
    pad.SetBottomMargin(0.25)
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.01)

    graphs = [
        createPull(d, s, prefix, title) for s, prefix in zip(style, ["d0", "z0", "qop"])
    ]
    all_graphs.append(graphs)

    graphs[0].Draw("AP")
    graphs[1].Draw("Psame")
    graphs[2].Draw("Psame")

    if i < 2:
        graphs[0].GetXaxis().SetLabelOffset(999)
        graphs[0].GetXaxis().SetLabelSize(0)
        graphs[0].GetXaxis().SetTitleSize(0)

    lines = [r.TLine(0.0, y, 3.0, y) for y in [0.0, -1.0, 1.0]]
    all_lines.append(lines)
    for line in lines:
        line.SetLineStyle(2)
        line.Draw()

canvas.cd()

text = createLabel(event_string="single muons, <#mu>=0", is_ttbar=False, x=0.75)

legend = createLegend(x1=0.75, y1=0.60, x2=0.9, y2=0.80)
legend.AddEntry(graphs[0], "p_{T}=1 GeV", "pl")
legend.AddEntry(graphs[1], "p_{T}=10 GeV", "pl")
legend.AddEntry(graphs[2], "p_{T}=100 GeV", "pl")
legend.Draw()

fake_graph = all_graphs[0][0].Clone()
fake_graph.SetMarkerStyle(20)
fake_graph.SetMarkerColor(r.kBlack)
fake_graph.SetMarkerSize(1.5)
fake_graph.SetLineColor(r.kBlack)

fake_legend = createLegend(x1=0.75, y1=0.45, x2=0.9, y2=0.55)
fake_legend.AddEntry(fake_graph, "mean", "p")
fake_legend.AddEntry(fake_graph, "RMS", "e")
fake_legend.Draw()

canvas.SaveAs(args.output)
