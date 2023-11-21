import ROOT as r
import numpy as np


def getDefaultStyle():
    return invertDictOfLists(
        {
            "marker": [20, 21, 22],
            "color": [r.kAzure + 1, r.kAzure + 2, r.kAzure + 3],
            "offset": [-0.07, 0.0, 0.07],
            "font": [43] * 3,
            "font_size": [22] * 3,
        }
    )


def invertDictOfLists(d):
    result = []
    keys = list(d.keys())
    for values in zip(*d.values()):
        result.append(dict(zip(keys, values)))
    return result


def createLabel(event_string, is_ttbar, x=0.18):
    text = r.TLatex()
    text.SetNDC()
    text.SetTextFont(43)
    text.SetTextSize(24)
    text.DrawLatex(x, 0.90, "ODD Simulation")
    text.DrawLatex(x, 0.85, event_string)
    if is_ttbar:
        text.DrawLatex(x, 0.80, r"#sqrt{s} = 14 TeV")
    return text


def createLegend(x1=0.68, y1=0.79, x2=1.0, y2=0.93):
    legend = r.TLegend(x1, y1, x2, y2)
    legend.SetTextFont(43)
    legend.SetTextSize(24)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    return legend


def createEfficency(efficiency, style):
    efficiency = efficiency.sort_values("eta")

    eta_error = 0.5 * (efficiency["eta"].values[1] - efficiency["eta"].values[0])
    eta_error = np.full(len(efficiency), eta_error)

    graph = r.TGraphAsymmErrors(
        len(efficiency),
        efficiency["eta"].values,
        efficiency["efficiency"].values,
        eta_error,
        eta_error,
        efficiency["lower_error"].values,
        efficiency["upper_error"].values,
    )

    graph.SetMarkerStyle(style["marker"])
    graph.SetMarkerSize(1.5)
    graph.SetMarkerColor(style["color"])
    graph.SetLineColor(style["color"])
    graph.SetTitle("")

    graph.GetXaxis().SetTitle("|#eta|")
    graph.GetXaxis().SetRangeUser(0.0, 3.0)
    graph.GetXaxis().SetTitleSize(28)
    graph.GetXaxis().SetTitleFont(style["font"])
    graph.GetXaxis().SetTitleOffset(1.0)
    graph.GetXaxis().SetLabelSize(28)
    graph.GetXaxis().SetLabelFont(style["font"])

    graph.GetYaxis().SetTitle("Technical efficiency")
    graph.GetYaxis().SetNdivisions(505)
    graph.GetYaxis().SetTitleSize(28)
    graph.GetYaxis().SetTitleFont(style["font"])
    graph.GetYaxis().SetTitleOffset(1.55)
    graph.GetYaxis().SetLabelSize(28)
    graph.GetYaxis().SetLabelFont(style["font"])

    return graph


def createResolution(resolution, style):
    resolution = resolution.sort_values("x")

    x_error = 0.5 * (resolution["x"].values[1] - resolution["x"].values[0])
    x_error = np.full(len(resolution), x_error)

    graph = r.TGraphErrors(
        len(resolution),
        resolution["x"].values,
        resolution["y"].values,
        x_error,
        resolution["y_err"].values,
    )

    graph.SetMarkerStyle(style["marker"])
    graph.SetMarkerSize(1.5)
    graph.SetMarkerColor(style["color"])
    graph.SetLineColor(style["color"])
    graph.SetTitle("")

    graph.GetXaxis().SetTitleSize(28)
    graph.GetXaxis().SetTitleFont(style["font"])
    graph.GetXaxis().SetTitleOffset(1.0)
    graph.GetXaxis().SetLabelSize(28)
    graph.GetXaxis().SetLabelFont(style["font"])

    graph.GetYaxis().SetNdivisions(505)
    graph.GetYaxis().SetTitleSize(28)
    graph.GetYaxis().SetTitleFont(style["font"])
    graph.GetYaxis().SetTitleOffset(1.55)
    graph.GetYaxis().SetLabelSize(28)
    graph.GetYaxis().SetLabelFont(style["font"])

    return graph


def createPull(pull, style, prefix, title):
    pull = pull.sort_values("eta")

    eta_error = np.zeros(len(pull))

    graph = r.TGraphErrors(
        len(pull),
        pull["eta"].values + style["offset"],
        pull[f"{prefix}_mean"].values,
        eta_error,
        pull[f"{prefix}_std"].values,
    )

    graph.SetMarkerStyle(style["marker"])
    graph.SetMarkerSize(1.5)
    graph.SetMarkerColor(style["color"])
    graph.SetLineColor(style["color"])
    graph.SetTitle("")

    graph.GetXaxis().SetTitle("|#eta|")
    graph.GetXaxis().SetRangeUser(0.0, 3.0)
    graph.GetXaxis().SetTitleFont(style["font"])
    graph.GetXaxis().SetTitleSize(style["font_size"])
    graph.GetXaxis().SetLabelFont(style["font"])
    graph.GetXaxis().SetLabelSize(style["font_size"])

    graph.GetYaxis().SetTitle(title)
    graph.GetYaxis().SetRangeUser(-2.0, 2.0)
    graph.GetYaxis().SetTitleFont(style["font"])
    graph.GetYaxis().SetTitleSize(style["font_size"])
    graph.GetYaxis().SetLabelFont(style["font"])
    graph.GetYaxis().SetLabelSize(style["font_size"])
    graph.GetYaxis().SetTitleOffset(1.5)
    graph.GetYaxis().CenterTitle()
    graph.GetYaxis().SetNdivisions(3)

    return graph
