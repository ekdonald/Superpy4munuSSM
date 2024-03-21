from Higgs import predictions as HP
from Higgs import Predictions
import matplotlib.pyplot as plt
import numpy as np


def chisqPlot(
    fig,
    ax,
    x,
    y,
    chisq,
    experiment=None,
    description=None,
    luminosity=None,
    reference=None,
):
    """
    Create a HiggsSignals style `chi^2` plot.

    Parameters
    ----------
     fig: the matplotlib figure to plot onto
     ax: the axis of fig to plot in
     x: gridpoints along the x-axis
     y: gridpoints along the y-axis
     chisq: chi^2 values on the `np.meshgrid(x,y)`
     experiment: experiment name for the label text
     description: description of the measurement for the label text
     luminosity: luminosity in fb-1 for the label text
     reference: a dict of the form `{68: ..., 95: ..., "BF": ...}` where
       `68` and `95`: `(xvalues, yvalues)` of the 68% and 95% CL contours to plot as reference
       `"BF"`: `(x, y)` of the reference best fit point
    """
    from matplotlib.lines import Line2D

    assert chisq.shape == (*x.shape, *y.shape)
    deltaChisq = chisq - np.min(chisq)
    cf = ax.contourf(
        x,
        y,
        deltaChisq,
        vmin=-1.5,
        levels=np.arange(0.0, 6.001, 0.01),
        cmap=plt.cm.get_cmap("GnBu").reversed(),
    )
    for c in cf.collections:
        c.set_rasterized(True)
    fig.colorbar(
        cf,
        ax=ax,
        label=r"$\Delta\chi^2$",
        ticks=np.arange(0, 6.1, 1),
        pad=0,
        fraction=0.1,
    )
    ax.contour(
        x,
        y,
        deltaChisq,
        levels=[2.3, 5.99],
        colors="k",
        linestyles=["-", "--"],
    )

    bfp = np.unravel_index(np.argmin(deltaChisq), deltaChisq.shape)
    ax.plot(x[bfp[1]], y[bfp[0]], marker="*", ls="none", c="k")

    percent = "\%" if plt.rcParams["text.usetex"] else "%"

    hsLegend = ax.legend(
        handles=[
            Line2D([0], [0], color="k", ls="-", label=f"HS 68{percent} CL"),
            Line2D([0], [0], color="k", ls="--", label=f"HS 95{percent} CL"),
            Line2D([0], [0], color="k", ls="none", marker="*", label="HS BFP"),
        ],
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(hsLegend)

    if reference:
        ax.plot(*reference[68], c="tab:gray")
        ax.plot(*reference[95], c="tab:gray", ls="--")
        ax.plot(*reference["BF"], c="tab:gray", marker="*", ls="none")
        refLegend = ax.legend(
            handles=[
                Line2D(
                    [0], [0], color="tab:gray", ls="-", label=f"ref. 68{percent} CL"
                ),
                Line2D(
                    [0], [0], color="tab:gray", ls="--", label=f"ref. 95{percent} CL"
                ),
                Line2D(
                    [0], [0], color="tab:gray", ls="none", marker="*", label="ref. BFP"
                ),
            ],
            loc="upper right",
            frameon=False,
        )
        ax.add_artist(refLegend)

    ax.text(
        0,
        1,
        r"\texttt{HiggsSignals}"
        if plt.rcParams["text.usetex"]
        else r"$\mathrm{HiggsSignals}$",
        horizontalalignment="left",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    exptext = ""
    if experiment is not None:
        if experiment == "LHCComb":
            experiment = "ATLAS + CMS"
        exptext = (
            rf"\textbf{{{str(experiment)}}}"
            if plt.rcParams["text.usetex"]
            else f"{experiment}"
        )
    desctext = f"{description}" if description is not None else ""
    lumitext = (
        rf"(${luminosity:.0f}\,\mathrm{{fb}}^{{-1}}$)" if luminosity is not None else ""
    )

    if exptext or desctext or lumitext:
        ax.text(
            1,
            1,
            " ".join(["using", exptext, desctext, lumitext]),
            horizontalalignment="right",
            verticalalignment="bottom",
            transform=ax.transAxes,
        )

    return fig, ax


def validateMeasurementKappa(
    meas, kappaF, kappaV, description, official={}, fig=None, ax=None
):
    """
    Create a validation plot for the measurement in the fermionic-bosonic kappa plane.

    Parameters
    ----------
     meas: the HS.Measurement to validate
     kappaF: grid values in `kappa_F`
     kappaV: grid values in `kappa_V`
     description: description of the measurement passed to chisqPlot
     official: passed to chisqPlot as reference, see there
     fig: an existing matplotlib figure to use, if None new `fig, ax` will be created
     ax: an existing matplotlib axis to use, if None new `fig, ax` will be created
    """

    def coups(kappaF, kappaV):
        def cgamgam(ctt, cbb, cWW, ctautau):
            c = np.array(
                [
                    662.84,
                    0.18,
                    14731.86,
                    -16.39,
                    -6249.93,
                    77.42,
                    0.21,
                    -17.69,
                    0.40,
                    83.59,
                ]
            )

            return np.sqrt(
                (
                    ctt ** 2 * c[0]
                    + cbb ** 2 * c[1]
                    + cWW ** 2 * c[2]
                    + ctt * cbb * c[3]
                    + ctt * cWW * c[4]
                    + cbb * cWW * c[5]
                    + ctautau ** 2 * c[6]
                    + ctt * ctautau * c[7]
                    + cbb * ctautau * c[8]
                    + ctautau * cWW * c[9]
                )
                / np.sum(c)
            )

        res = HP.scaledSMlikeEffCouplings(kappaF)
        res.ZZ = kappaV
        res.WW = kappaV
        res.gamgam = cgamgam(kappaF, kappaF, kappaV, kappaF)
        return res

    pred = Predictions()
    h = pred.addParticle(HP.NeutralScalar("h", "even"))
    refMass = meas.referenceMass()
    h.setMass(refMass)

    def getChisq(kappaV, kappaF):
        HP.effectiveCouplingInput(
            h, coups(kappaF, kappaV), reference=meas.referenceModel()
        )
        return meas(pred, {})

    testAgainstMeasurement = np.frompyfunc(getChisq, 2, 1)

    XX, YY = np.meshgrid(kappaV, kappaF)
    dat = testAgainstMeasurement(XX, YY)
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    chisqPlot(
        fig,
        ax,
        kappaV,
        kappaF,
        dat,
        experiment=meas.experiment().name,
        description=description,
        luminosity=meas.luminosity(),
        reference=official,
    )
    ax.plot(1, 1, marker="d", c="tab:orange", ls="none", label="SM")
    ax.legend(loc="upper center", frameon=False)
    ax.set_xlabel(r"$\kappa_V$")
    ax.set_ylabel(r"$\kappa_F$")

    return fig, ax


def validateMeasurementModFactors(
    meas, mu1Range, mu2Range, modFacs, description, fig=None, ax=None
):
    """
    Create a validation plot for the measurement in a plane spanned by two modification factors.

    Parameters
    ----------
     meas: the HS.Measurement to validate
     mu1Range: grid values of the first modification factor
     mu2Range: grid values of the second modification factor
     modFacs: a function taking one of each modification factor values and
     returning a dictionary of the corresponding modification factors for all
     the sub-measurements of meas
     description: description of the measurement passed to chisqPlot
     fig: an existing matplotlib figure to use, if None new `fig, ax` will be created
     ax: an existing matplotlib axis to use, if None new `fig, ax` will be created
    """
    pred = Predictions()
    h = pred.addParticle(HP.NeutralScalar("h", "even"))
    h.setMass(meas.referenceMass())
    HP.effectiveCouplingInput(h, HP.smLikeEffCouplings, reference=meas.referenceModel())

    def getChisq(mu1, mu2):
        return meas(pred, modFacs(mu1, mu2))

    testAgainstMeasurement = np.frompyfunc(getChisq, 2, 1)
    XX, YY = np.meshgrid(mu1Range, mu2Range)
    dat = testAgainstMeasurement(XX, YY)

    if fig is None or ax is None:
        fig, ax = plt.subplots()
    chisqPlot(
        fig,
        ax,
        mu1Range,
        mu2Range,
        dat,
        description=description,
        experiment=meas.experiment().name,
        luminosity=meas.luminosity(),
    )
    return fig, ax


def validateMeasurementRates(
    meas, rate1Range, rate2Range, setRates, description, fig=None, ax=None
):
    """
    Create a validation plot for the measurement in a plane of signal rates.

    Parameters
    ----------
     meas: the HS.Measurement to validate
     rate1Range: grid values of the first rate
     rate2Range: grid values of the second rate
     setRates: a function taking a HP.BsmParticle and one of each rate values
     and setting the needed rates on the particle
     description: description of the measurement passed to chisqPlot
     fig: an existing matplotlib figure to use, if None new `fig, ax` will be created
     ax: an existing matplotlib axis to use, if None new `fig, ax` will be created
    """
    pred = Predictions()
    h = pred.addParticle(HP.NeutralScalar("h", "even"))
    h.setMass(meas.referenceMass())

    def getChisq(r1, r2):
        setRates(h, r1, r2)
        return meas(pred, {})

    testAgainstMeasurement = np.frompyfunc(getChisq, 2, 1)
    XX, YY = np.meshgrid(rate1Range, rate2Range)
    dat = testAgainstMeasurement(XX, YY)

    if fig is None or ax is None:
        fig, ax = plt.subplots()
    chisqPlot(
        fig,
        ax,
        rate1Range,
        rate2Range,
        dat,
        description=description,
        experiment=meas.experiment().name,
        luminosity=meas.luminosity(),
    )
    return fig, ax


def validateMassMeasurement2Higgses(
    meas,
    channel,
    m1Range,
    description,
    m2Range=None,
    dm1=1.0,
    dm2=1.0,
    rateRatio=1.0,
    fig=None,
    ax=None,
):
    """
    Create a validation plot for a mass measurement using 2 particles with different masses.

    Parameters
    ----------
     meas: the HS.Measurement to validate
     channel: the channel the mass was measured in
     m1Range: grid values of the first mass
     description: description of the measurement passed to chisqPlot
     m2Range: grid values of the second mass, if None this is equal to m1Range
     dm1: theoretical mass uncertainty in GeV of the first particle
     dm2: theoretical mass uncertainty in GeV of the second particle
     rateRatio: ratio between the rates of the first and the second particle
     fig: an existing matplotlib figure to use, if None new `fig, ax` will be created
     ax: an existing matplotlib axis to use, if None new `fig, ax` will be created
    """

    if m2Range is None:
        m2Range = m1Range
    m1Range = np.array(m1Range, dtype=float)
    m2Range = np.array(m2Range, dtype=float)
    pred = Predictions()
    h1 = pred.addParticle(HP.NeutralScalar("h1", "even"))
    h2 = pred.addParticle(HP.ChargedScalar("h2"))
    h1.setChannelRate(*channel, 1)
    h2.setChannelRate(*channel, 1 / rateRatio)
    h1.setMassUnc(dm1)
    h2.setMassUnc(dm2)

    def getChisq(m1, m2):
        h1.setMass(m1)
        h2.setMass(m2)
        return meas(pred, {})

    testAgainstMeasurement = np.frompyfunc(getChisq, 2, 1)
    XX, YY = np.meshgrid(m1Range, m2Range)
    dat = testAgainstMeasurement(XX, YY)
    dat = np.array(dat, dtype=float)
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    chisqPlot(
        fig,
        ax,
        m1Range,
        m2Range,
        dat,
        experiment=meas.experiment().name,
        description=description,
        luminosity=meas.luminosity(),
    )
    ax.set_ylim(*ax.get_ylim())  # fix ylim
    ax.plot(
        m1Range,
        (1 + rateRatio) * (meas.referenceMass() - m1Range / (1 + 1 / rateRatio)),
        label=r"$\overline{m} = \hat{m}$",
        color="tab:orange",
        ls=":",
    )
    ax.legend(loc="upper right", frameon=False)
    ax.set_xlabel(rf"$m_{{h_1}} (\Delta m_{{h_1}}={dm1}\,\mathrm{{GeV}})$")
    ax.set_ylabel(rf"$m_{{h_2}} (\Delta m_{{h_2}}={dm2}\,\mathrm{{GeV}})$")

    return fig, ax
