import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
from matplotlib.lines import Line2D


def determineChannelCharge(prod, decay):
    """
    Determine the electric charge of a channel. If multiple charges are
    possible, returns the lowest charge.
    """
    from .._Higgs.predictions import validProductionFor, validDecayFor, ECharge

    channelCharges = [
        c
        for c in ECharge.__members__
        if validProductionFor(prod, c) and validDecayFor(decay, c)
    ]
    assert len(channelCharges) >= 1
    return ECharge(channelCharges[0])


def limExtent(limDat: dict):
    minLim = min(
        np.min(limDat["analysis"]["limit"]["observed"]),
        np.min(limDat["analysis"]["limit"]["expected"]),
    )
    maxLim = max(
        np.max(limDat["analysis"]["limit"]["observed"]),
        np.max(limDat["analysis"]["limit"]["expected"]),
    )
    return (minLim, maxLim)


def massExtent(limDat: dict):
    return (
        np.min(limDat["analysis"]["grid"]["mass"]),
        np.max(limDat["analysis"]["grid"]["mass"]),
    )


def widthExtent(limDat: dict):
    return (
        np.min(limDat["analysis"]["grid"]["width"]),
        np.max(limDat["analysis"]["grid"]["width"]),
    )


def normalized(limDat: dict):
    return "normalization" in limDat.keys()


def setupPredictions(limDat: dict, setAdditionalRates, cp: str):
    from .. import Predictions
    from .._Higgs import predictions as HP
    from random import randrange
    from math import sqrt

    pred = Predictions()
    p, d = limDat["process"]["channels"][randrange(len(limDat["process"]["channels"]))]
    if determineChannelCharge(p, d) == HP.ECharge.neutral:
        particle = pred.addParticle(HP.NeutralScalar("h", cp))
        if normalized(limDat):
            print("using effective coupling input for SM-normalized limit")
            if limDat["process"]["channels"][0][1] == "inv":
                print("assuming only production is normalized for inv decay")

                def setRate(r):
                    particle.resetChannelRates()
                    HP.effectiveCouplingInput(
                        particle, HP.scaledSMlikeEffCouplings(sqrt(r))
                    )
                    particle.setTotalWidth(0.0)
                    particle.setDecayWidth("directInv", 1.0)
                    setAdditionalRates(particle, r)

            else:

                def setRate(r):
                    particle.resetChannelRates()
                    HP.effectiveCouplingInput(
                        particle, HP.scaledSMlikeEffCouplings(sqrt(r))
                    )
                    setAdditionalRates(particle, r)

        else:
            if d == HP.Decay.none:

                def setRate(r):
                    if limDat["collider"] == HP.Collider.LEP:
                        particle.setNormalizedCxn(limDat["collider"], p, r)
                    else:
                        particle.setCxn(limDat["collider"], p, r)
                    setAdditionalRates(particle, r)

            else:

                def setRate(r):
                    particle.resetChannelRates()
                    particle.setChannelRate(limDat["collider"], p, d, r)
                    setAdditionalRates(particle, r)

    elif determineChannelCharge(p, d) == HP.ECharge.single:
        particle = pred.addParticle(HP.ChargedScalar("h"))

        def setRate(r):
            particle.setChannelRate(limDat["collider"], p, d, r)
            setAdditionalRates(particle, r)

    return (pred, particle, setRate)


def validateChannelLimit(
    lim, fig=None, ax=None, cp: str = "undefined", setAdditionalRates=lambda *args: None
):
    with open(lim.loadedFrom()) as limJson:
        limDat = json.load(limJson)
    pred, particle, setRate = setupPredictions(limDat, setAdditionalRates, cp)

    def testAgainstChannelLimit(mass, rate):
        particle.setMass(mass)
        setRate(rate)
        app = lim.apply(pred)
        if len(app) == 0:
            return (0, 0)
        return (app[0].expRatio(), app[0].obsRatio())

    test = np.frompyfunc(testAgainstChannelLimit, 2, 2)

    rates = np.linspace(*limExtent(limDat))
    masses = np.linspace(*massExtent(limDat))
    XX, YY = np.meshgrid(masses, rates)
    expR, obsR = test(XX, YY)
    if fig is None or ax is None:
        fig, ax = plt.subplots()
    ax.contour(masses, rates, expR, levels=[1], colors=["darkblue"])
    ax.contour(masses, rates, obsR, levels=[1], colors=["firebrick"])
    expLine = ax.plot(
        limDat["analysis"]["grid"]["mass"],
        limDat["analysis"]["limit"]["expected"],
        ls="--",
        c="deepskyblue",
        label="expected",
    )
    obsLine = ax.plot(
        limDat["analysis"]["grid"]["mass"],
        limDat["analysis"]["limit"]["observed"],
        ls="--",
        c="orange",
        label="observed",
    )
    ax.legend(
        handles=[
            Line2D([0], [0], color="darkblue", label="expRatio = 1"),
            Line2D([0], [0], color="firebrick", label="obsRatio = 1"),
            expLine[0],
            obsLine[0],
        ],
        ncol=2,
        title="official                     HB",
    )
    ax.set_xlabel("M [GeV]")
    if not normalized(limDat):
        ax.set_ylabel("rate [pb]")
        ax.set_yscale("log")
    else:
        ax.set_ylabel(r"$\mu$")
    ax.set_title(lim.processDesc(), fontsize="small")
    ax.grid()
    return (fig, ax)


def validateChannelWidthLimit(
    lim,
    cp: str = "undefined",
    setAdditionalRates=lambda *args: None,
    widthsToPlot=None,
    descInTitle=True,
):
    from .._Higgs import predictions as HP

    with open(lim.loadedFrom()) as limJson:
        limDat = json.load(limJson)
    pred, particle, setRate = setupPredictions(limDat, setAdditionalRates, cp)
    relativeWidth = limDat["analysis"]["relativeWidth"]

    rates = np.linspace(*limExtent(limDat))
    masses = np.linspace(*massExtent(limDat))
    XX, YY = np.meshgrid(masses, rates)

    mGrid = limDat["analysis"]["grid"]["mass"]
    wGrid = limDat["analysis"]["grid"]["width"]
    for i, width in enumerate(wGrid):
        if widthsToPlot is None or width in widthsToPlot:

            def testAgainstChannelWidthLimit(mass, rate):
                particle.setMass(mass)
                setRate(rate)
                if relativeWidth:
                    particle.setTotalWidth(width * mass)
                else:
                    particle.setTotalWidth(width)
                app = lim.apply(pred)
                if len(app) == 0:
                    return (0, 0)
                return (app[0].expRatio(), app[0].obsRatio())

            test = np.frompyfunc(testAgainstChannelWidthLimit, 2, 2)

            expR, obsR = test(XX, YY)
            fig, ax = plt.subplots()
            ax.contour(masses, rates, expR, levels=[1], colors=["darkblue"])
            ax.contour(masses, rates, obsR, levels=[1], colors=["firebrick"])
            expLine = ax.plot(
                mGrid,
                limDat["analysis"]["limit"]["expected"][i :: len(wGrid)],
                ls="--",
                c="deepskyblue",
                label="expected",
            )
            obsLine = ax.plot(
                mGrid,
                limDat["analysis"]["limit"]["observed"][i :: len(wGrid)],
                ls="--",
                c="orange",
                label="observed",
            )
            ax.legend(
                handles=[
                    Line2D([0], [0], color="darkblue", label="expRatio = 1"),
                    Line2D([0], [0], color="firebrick", label="obsRatio = 1"),
                    expLine[0],
                    obsLine[0],
                ],
                ncol=2,
                title="official                     HB",
            )
            ax.set_xlabel("M [GeV]")
            if not normalized(limDat):
                ax.set_ylabel("rate [pb]")
                ax.set_yscale("log")
            else:
                ax.set_ylabel(r"$\mu$")
            if descInTitle:
                ax.set_title(
                    lim.processDesc() + " width=" + str(width), fontsize="small"
                )
            else:
                ax.set_title(f"width={width}", fontsize="small")
            ax.grid()
    return


def runLlhLimitOnPlane(lim, benchmarkPlane, massKeys, channels, rateKeys):
    from .. import Predictions
    from .._Higgs import predictions

    assert np.array(rateKeys).shape == (len(massKeys), len(channels))

    hid = ["h" + str(i) for i in range(len(massKeys))]
    pred = Predictions()
    particles = [pred.addParticle(predictions.NeutralScalar(h)) for h in hid]

    def testAgainstLlhLimit(series):
        for p, m, rs in zip(particles, massKeys, rateKeys):
            p.setMass(series[m])
            for r, ch in zip(rs, channels):
                p.setChannelRate(ch, series[r])
        app = lim.apply(pred)
        if len(app) > 0:
            mostSensitive = np.argmax([a.expRatio() for a in app])
            return pd.Series(
                {
                    "expR": app[mostSensitive].expRatio(),
                    "obsR": app[mostSensitive].obsRatio(),
                    "expLlh": app[mostSensitive].expLikelihood(),
                    "obsLlh": app[mostSensitive].obsLikelihood(),
                }
            )
        else:
            return pd.Series({"expR": 0.0, "obsR": 0.0, "expLlh": 0.0, "obsLlh": 0.0})

    return pd.concat(
        [benchmarkPlane, benchmarkPlane.apply(testAgainstLlhLimit, axis=1)], axis=1
    )


def plotLlhLimitInPlane(dfr, xKey, yKey, fig=None, axes=None, cmap="viridis"):
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    if fig is None or axes is None or len(axes) < 2:
        fig, axes = plt.subplots(
            figsize=(10, 4), nrows=1, ncols=2, sharex=True, sharey=True
        )
    for ax, t in zip(axes, ["exp", "obs"]):
        tcplot = ax.tricontour(
            dfr[xKey], dfr[yKey], dfr[t + "R"], levels=[1.0], colors=["w"]
        )
        tcplot.collections[0].set_label(r"$\mathrm{CL}_S >0.05$")
        tcplot = ax.tricontour(
            dfr[xKey], dfr[yKey], dfr[t + "Llh"], levels=[5.99], colors=["k"]
        )
        tcplot.collections[0].set_label(r"$\Delta\chi^2<5.99$")
        cf = ax.tricontourf(
            dfr[xKey],
            dfr[yKey],
            dfr[t + "Llh"],
            levels=np.linspace(0, 20),
            norm=Normalize(0, 20),
            extend="max",
            zorder=0,
            cmap=cmap,
        )
        for c in cf.collections:
            c.set_edgecolor("face")
        ax.grid()
        ax.set_rasterization_zorder(0)
    axes[0].set_title("expected")
    axes[1].set_title("observed")
    sm = ScalarMappable(norm=Normalize(0, 20), cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=axes.flatten(), label=r"$-2\,\mathrm{ln}(L)$", extend="max")
    return fig, axes
