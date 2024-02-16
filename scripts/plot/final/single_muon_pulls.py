import argparse

import pandas as pd
import ROOT as r

from mycommon.root import getDefaultStyle, createLabel, createLegend, createPull


style = getDefaultStyle()
labels = [
    "(d_{0}^{reco}-d_{0}^{true})/#sigma_{d_{0}}^{reco}",
    "(z_{0}^{reco}-z_{0}^{true})/#sigma_{z_{0}}^{reco}",
    "(q/p^{reco}-q/p^{true})/#sigma_{q/p}^{reco}",
]

parser = argparse.ArgumentParser()
parser.add_argument("mu_1GeV", help="Path to the 1 GeV muon pull csv file")
parser.add_argument("mu_10GeV", help="Path to the 10 GeV muon pull csv file")
parser.add_argument("mu_100GeV", help="Path to the 100 GeV muon pull csv file")
parser.add_argument("output", help="Path to the output pdf file")
args = parser.parse_args()

files = [args.mu_1GeV, args.mu_10GeV, args.mu_100GeV]
data = [pd.read_csv(file) for file in files]

canvas = r.TCanvas("canvas", "", 800, 700)

canvas.Divide(1, 4)

all_graphs = []
all_lines = []

for i, d, title in zip([2, 3, 4], data, labels):
    canvas.cd(i)

    r.gPad.SetLeftMargin(0.15)
    r.gPad.SetRightMargin(0.02)
    r.gPad.SetTopMargin(0.05)
    r.gPad.SetBottomMargin(0.3)

    graphs = [
        createPull(d, s, prefix, title) for s, prefix in zip(style, ["d0", "z0", "qop"])
    ]
    all_graphs.append(graphs)

    graphs[0].Draw("AP")
    graphs[1].Draw("Psame")
    graphs[2].Draw("Psame")

    lines = [r.TLine(0.0, y, 3.0, y) for y in [0.0, -1.0, 1.0]]
    all_lines.append(lines)
    for line in lines:
        line.SetLineStyle(2)
        line.Draw()

canvas.cd(0)

text = createLabel(event_string="single muons, <#mu>=0", is_ttbar=False, x=0.09)

legend = createLegend()
legend.AddEntry(graphs[0], "p_{T}=1 GeV", "pl")
legend.AddEntry(graphs[1], "p_{T}=10 GeV", "pl")
legend.AddEntry(graphs[2], "p_{T}=100 GeV", "pl")
legend.Draw()

fake_graph = all_graphs[0][0].Clone()
fake_graph.SetMarkerStyle(20)
fake_graph.SetMarkerColor(r.kBlack)
fake_graph.SetMarkerSize(1.5)
fake_graph.SetLineColor(r.kBlack)

fake_legend = createLegend(x1=0.5, y1=0.85, x2=0.7)
fake_legend.AddEntry(fake_graph, "mean", "p")
fake_legend.AddEntry(fake_graph, "RMS", "e")
fake_legend.Draw()

canvas.SaveAs(args.output)
