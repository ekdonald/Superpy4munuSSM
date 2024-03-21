import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcol
from .ImplementationUtils import printable


def chisqCritBound(ax, x, y, llh, chisqCrit=5.99):
    """
    Plot the contour where llh=chisqCrit and return the set of points that lie
    on it.
    """
    c = ax.tricontour(
        x,
        y,
        llh,
        colors=["tab:gray"],
        levels=[chisqCrit],
    )
    return np.concatenate([path.vertices for path in c.collections[0].get_paths()])


def fitEllipse(points, center):
    """
    Fit an ellipse around the given center through the points.
    """
    X, Y = (points - center).transpose()
    A = np.transpose([X ** 2, 2 * X * Y, Y ** 2])
    b = np.ones_like(X)
    ellipse = np.linalg.lstsq(A, b, rcond=None)[0].squeeze()
    return ellipse


def longAxis(ellipse):
    """
    Return the long axis of an ellipse from the parameters fitted by fitEllipse.
    """
    mat = np.array([[ellipse[0], ellipse[1]], [ellipse[1], ellipse[2]]])
    ev, vec = np.linalg.eig(mat)
    return vec[np.argmax(np.abs(ev))]


def insideEllipse(ellipse, x, y):
    """
    Check if (x,y) lies within the (origin-centered) ellipse.
    """
    return x ** 2 / ellipse[0] ** 2 + y ** 2 / ellipse[1] ** 2 < 1


def constructRegularization(ax, x, y, llh):
    """
    Construct the regularizing ellipse within which the llh will be set to zero
    and plot the steps. See the HB-5 manual for details.
    """
    from operator import truediv

    ax.set_title("original")

    points = chisqCritBound(ax, x, y, llh)

    bfp = np.array([x.iloc[np.argmin(llh)], y.iloc[np.argmin(llh)]])

    ellipse = fitEllipse(points, bfp)
    ax.tricontour(
        x,
        y,
        ellipse[0] * (x - bfp[0]) ** 2
        + 2 * ellipse[1] * (x - bfp[0]) * (y - bfp[1])
        + ellipse[2] * (y - bfp[1]) ** 2,
        levels=[1.0],
        colors=["k"],
    )

    m = -np.abs(truediv(*longAxis(ellipse)))
    ax.plot(np.unique(x), m * (np.unique(x) - bfp[0]) + bfp[1], ls="--", c="k")

    a = max(np.sqrt(bfp[0] ** 2 - bfp[0] * bfp[1] / m), np.unique(x)[1])
    b = max(np.sqrt(-m * bfp[0] * bfp[1] + bfp[1] ** 2), np.unique(x)[1])
    t = np.linspace(0, np.pi / 2)
    ax.plot(
        a * np.cos(t),
        b * np.sin(t),
        c="red",
        label=r"$a = {:.4}, b = {:.4}$".format(a, b),
    )

    cf = ax.tricontourf(
        x,
        y,
        llh,
        norm=mcol.Normalize(0, 20),
        levels=np.linspace(0, 20, 100),
        extend="max",
        cmap="viridis",
        zorder=-1,
    )
    for c in cf.collections:
        c.set_edgecolor("face")

    return (a, b)


def performRegularization(ax, zeroEllipse, dfLlh, xKey, yKey, llhKey, chisqCrit=5.99):
    """
    Regularize the llh given a previously constructed ellipse that defines the
    zero region.
    """
    from scipy.interpolate import RegularGridInterpolator
    from scipy.optimize import root_scalar

    ax.set_title("processed")

    x = np.unique(dfLlh[xKey])
    y = np.unique(dfLlh[yKey])
    XX, YY = np.meshgrid(x, y)
    assert np.allclose(XX.flatten(), dfLlh[xKey])
    assert np.allclose(YY.flatten(), dfLlh[yKey])
    interp = RegularGridInterpolator(
        (x, y), dfLlh[llhKey].values.reshape(XX.shape).transpose()
    )
    assert np.allclose(
        interp((XX.flatten(), YY.flatten())),
        dfLlh[llhKey].values,
    )

    def r2s(p):
        if p.obs > 0:
            angle = np.arctan2(p[yKey], p[xKey])

            def fun(r):
                return interp(r * np.array((np.cos(angle), np.sin(angle)))) - chisqCrit

            sol = root_scalar(
                fun,
                bracket=(
                    p.re,
                    0.99999
                    * min(
                        np.max(dfLlh[xKey]) / max(np.cos(angle), 1e-10),
                        np.max(dfLlh[yKey]) / max(np.sin(angle), 1e-10),
                    ),
                ),
            )
            return sol.root
        return 0.0

    dfLlh.loc[insideEllipse(zeroEllipse, dfLlh[xKey], dfLlh[yKey]), llhKey] = 0
    dfLlh["r0"] = np.sqrt(dfLlh[xKey] ** 2 + dfLlh[yKey] ** 2)
    dfLlh["re"] = (zeroEllipse[0] * zeroEllipse[1] * dfLlh.r0) / np.sqrt(
        zeroEllipse[0] ** 2 * dfLlh[yKey] ** 2 + zeroEllipse[1] ** 2 * dfLlh[xKey] ** 2
    )
    dfLlh["r2s"] = dfLlh[(dfLlh[llhKey] > 0) & (dfLlh[llhKey] < chisqCrit)].apply(
        r2s, axis=1
    )
    dfLlh.loc[(dfLlh[llhKey] > 0) & (dfLlh[llhKey] < chisqCrit), llhKey] = (
        chisqCrit / (dfLlh.r2s - dfLlh.re) ** 2 * (dfLlh.r0 - dfLlh.re) ** 2
    )
    dfLlh.drop(columns=["r0", "re", "r2s"], inplace=True)

    cf = ax.tricontourf(
        dfLlh[xKey],
        dfLlh[yKey],
        dfLlh[llhKey],
        norm=mcol.Normalize(0, 20),
        levels=np.linspace(0, 20, 100),
        extend="both",
        cmap="viridis",
        zorder=-1,
    )
    for c in cf.collections:
        c.set_edgecolor("face")
    ax.tricontour(
        dfLlh[xKey],
        dfLlh[yKey],
        dfLlh[llhKey],
        levels=[chisqCrit],
        colors=["tab:gray"],
    )


def regularizeExclusionLlh(dfLlh, xKey, yKey, llhKey="obs", chisqCrit=5.99):
    """
    Regularize the exclusion likelihood to avoid exclusion of low rates. See the
    HB-5 manual for details on the method.

    Parameters
    ----------
    dfLlh : pandas.DataFrame
        a dataframe containing the values to regularize
    xKey : str
        key for the x-values
    yKey : str
        key for the y-values
    llhKey : str
        key for the likelihood values
    chisqCrit : float
        the chisq cut to use in the construction
    """
    from matplotlib.cm import ScalarMappable

    fig, axes = plt.subplots(figsize=(10, 4), ncols=2, sharex=True, sharey=True)
    zeroEllipse = constructRegularization(
        axes[0], dfLlh[xKey], dfLlh[yKey], dfLlh[llhKey]
    )

    performRegularization(axes[1], zeroEllipse, dfLlh, xKey, yKey, llhKey)

    for ax in axes:
        ax.set_aspect("equal", "box", share=True)
        ax.set_xlim(np.min(dfLlh[xKey]), np.max(dfLlh[xKey]))
        ax.set_ylim(np.min(dfLlh[yKey]), np.max(dfLlh[yKey]))
        ax.grid()
        ax.set_rasterization_zorder(0)
    return fig, axes


def commonRatePlane(
    obs, exp, keys, obsLlhKey="obs", expLlhKey="exp", floatPrecision=14
):
    """
    Match the observed and expected likelihood onto a common grid in the rates.

    Assumes (and verifies) that the original grids are in row-major (Fortran)
    order. The returned grids are in column-major (C) order as needed for the
    likelihood limits.

    Parameters
    ----------
    obs: pandas.DataFrame
        a dataframe containing the observed likelihood and the corresponding grid
    exp: pandas.DataFrame
        a dataframe containing the expected likelihood and the corresponding grid
    keys: list of str
        the keys in obs and exp that correspond to the rates. The ordering of the keys has to match a np.meshgrid.
    obsLlhKey: str
        key of the likelihood column in obs
    expLlhKey: str
        key of the likelihood columns in exp
    floatPrecision: int
        number of significant digits to round to
    """
    from scipy.interpolate import RegularGridInterpolator

    uniqueObsGrids = [printable(np.unique(obs[k]), floatPrecision) for k in keys]
    assert np.all([np.min(x) == 0 for x in uniqueObsGrids])
    uniqueExpGrids = [printable(np.unique(exp[k]), floatPrecision) for k in keys]
    assert np.all([np.min(x) == 0 for x in uniqueExpGrids])
    assert np.all(
        [
            np.allclose(grid.flatten(), obs[k])
            for grid, k in zip(np.meshgrid(*uniqueObsGrids), keys)
        ]
    ), "observed llh is not in row-major order for the given keys"
    assert np.all(
        [
            np.allclose(grid.flatten(), exp[k])
            for grid, k in zip(np.meshgrid(*uniqueExpGrids), keys)
        ]
    ), "expected llh is not in row-major order for the given keys"
    if np.all(
        [np.allclose(*gridpair) for gridpair in zip(uniqueObsGrids, uniqueExpGrids)]
    ):
        obsDat = np.reshape(
            obs[obsLlhKey].values, [len(x) for x in uniqueObsGrids], order="F"
        ).flatten(order="C")
        expDat = np.reshape(
            exp[expLlhKey].values, [len(x) for x in uniqueExpGrids], order="F"
        ).flatten(order="C")
        return {
            "channels": uniqueObsGrids,
            "observed": printable(obsDat, floatPrecision),
            "expected": printable(expDat, floatPrecision),
        }

    else:
        grid = [printable(np.union1d(exp[k], obs[k]), floatPrecision) for k in keys]
        obsI = RegularGridInterpolator(
            uniqueObsGrids,
            np.reshape(
                obs[obsLlhKey].values, [len(x) for x in uniqueObsGrids], order="F"
            ),
            fill_value=np.nan,
            bounds_error=False,
            method="linear",
        )
        obsE = RegularGridInterpolator(
            uniqueObsGrids,
            np.reshape(
                obs[obsLlhKey].values, [len(x) for x in uniqueObsGrids], order="F"
            ),
            method="nearest",
            bounds_error=False,
            fill_value=None,
        )
        assert np.allclose(
            obsI(np.transpose(np.meshgrid(*uniqueObsGrids, indexing="ij"))).flatten(),
            obs[obsLlhKey].values,
        )

        expI = RegularGridInterpolator(
            uniqueExpGrids,
            np.reshape(
                exp[expLlhKey].values, [len(x) for x in uniqueExpGrids], order="F"
            ),
            fill_value=np.nan,
            bounds_error=False,
            method="linear",
        )
        expE = RegularGridInterpolator(
            uniqueExpGrids,
            np.reshape(
                exp[expLlhKey].values, [len(x) for x in uniqueExpGrids], order="F"
            ),
            method="nearest",
            bounds_error=False,
            fill_value=None,
        )
        assert np.allclose(
            expI(np.transpose(np.meshgrid(*uniqueExpGrids, indexing="ij"))).flatten(
                order="C"
            ),
            exp[expLlhKey].values,
        )
        points = np.transpose(np.meshgrid(*grid))
        obsDat = obsI(points).flatten(order="C")
        obsDat[~np.isfinite(obsDat)] = obsE(points).flatten(order="C")[
            ~np.isfinite(obsDat)
        ]
        expDat = expI(points).flatten(order="C")
        expDat[~np.isfinite(expDat)] = expE(points).flatten(order="C")[
            ~np.isfinite(expDat)
        ]
        return {
            "channels": grid,
            "observed": printable(obsDat, floatPrecision),
            "expected": printable(expDat, floatPrecision),
        }
