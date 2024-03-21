"""
Helper functions for the implementation of experimental results and theoretical
predictions.
"""

import numpy as np
from typing import Union
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict


def relevantGridEntries(x, values, rtol: float = 1e-3):
    """
    For values on 1d grid x, return the indices of those values that are not
    simple linear interpolations between the adjacent values.
    """
    from scipy.interpolate import interp1d

    relevantentries = list(range(0, len(x)))
    p = 0
    c = 1
    n = 2
    while n < len(x):
        while n < len(x) and np.allclose(
            interp1d((x[p], x[n]), (values[p], values[n]), axis=0)(x[c]),
            values[c],
            rtol=rtol,
            atol=np.finfo(float).resolution,
        ):
            relevantentries.remove(c)
            c += 1
            n += 1
        p = c
        c = n
        n += 1
    return relevantentries


def fromHB5Table1(
    tablename: str,
    category: str,
):
    address = "https://gitlab.com/higgsbounds/higgsbounds/-/raw/master/data/Expt_tables/{:.3}tables/{}.txt".format(
        category, tablename
    )
    df = read_csv_from_web(
        address, sep=r"\s+", names=["m", "obs", "exp"], comment="#", dtype=np.float64
    )
    rel = relevantGridEntries(df.m, df[["exp", "obs"]].values)
    return df.loc[rel].reset_index(drop=True)


def fromHB5Table2(
    tablename: str,
    category: str,
    p1Name: str = "Mother",
    p2Name: str = "Daughter",
    manuallySelectedM1=None,
    manuallySelectedM2=None,
):
    address = "https://gitlab.com/higgsbounds/higgsbounds/-/raw/master/data/Expt_tables/{:.3}tables/{:.3}tables2/{}_{{}}.txt".format(
        category, category, tablename
    )
    obsRaw = read_csv_from_web(address.format("obs"), sep=r"\s+")
    expRaw = read_csv_from_web(address.format("pred"), sep=r"\s+")
    assert list(obsRaw) == list(expRaw), "matching obs and pred grids"
    for d in (obsRaw, expRaw):
        d.set_axis(d.columns.astype(float), axis="columns", inplace=True)
        d.set_index(-100, inplace=True)
    assert np.allclose(obsRaw.index, expRaw.index), "matching obs and pred grids"

    if manuallySelectedM1 is not None:
        selectedM1 = manuallySelectedM1
    else:
        selectedM1 = obsRaw.index

    if manuallySelectedM2 is not None:
        selectedM2 = manuallySelectedM2
    else:
        selectedM2 = obsRaw.columns[
            np.unique(
                np.concatenate(
                    [
                        relevantGridEntries(
                            obsRaw.columns.values, obsRaw.loc[m1, :].values
                        )
                        for m1 in selectedM1
                    ]
                )
            )
        ]

    if manuallySelectedM1 is None:
        selectedM1 = obsRaw.index[
            np.unique(
                np.concatenate(
                    [
                        relevantGridEntries(obsRaw.index.values, obsRaw[m2].values)
                        for m2 in selectedM2
                    ]
                )
            )
        ]
    m1Vals, m2Vals = np.meshgrid(selectedM1, selectedM2, indexing="ij")

    df = pd.DataFrame(
        data={
            "m{}".format(p1Name): m1Vals.flatten(),
            "m{}".format(p2Name): m2Vals.flatten(),
            "obs": obsRaw.loc[selectedM1, selectedM2].values.flatten(),
            "exp": expRaw.loc[selectedM1, selectedM2].values.flatten(),
        }
    )

    return df


def read_csv_from_web(url: str, **kwargs):
    """
    Read a csv from a url into pandas.
    """
    from requests import get
    from io import StringIO
    from pandas import read_csv

    header = {
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.75 Safari/537.36",
        "X-Requested-With": "XMLHttpRequest",
    }
    response = get(url, headers=header)
    file_object = StringIO(response.content.decode("utf-8").replace("D", "e"))
    return read_csv(file_object, **kwargs)


def readHEPDataCsv(url: str, skip=0, **kwargs):
    """
    Read a HEPData csv into (several) pandas dataframes.
    """
    from requests import get
    from io import StringIO
    from pandas import read_csv

    header = {
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.75 Safari/537.36",
        "X-Requested-With": "XMLHttpRequest",
    }
    response = get(url, headers=header)
    parts = response.content.decode("utf-8").split("\n\n")
    parts.remove("")

    dfs = []
    for p in parts[skip:]:
        try:
            dfs.append(read_csv(StringIO(p), comment="#", **kwargs))
        except ValueError:
            pass
    return dfs


def assembleMetadata(
    inspireIdentifier: Union[str, int],
    source: str,
    idSuffix: int = None,
    luminosity: float = None,
    collider: str = "",
    experiment: str = "",
):
    from .Inspire import getMetadata

    data = getMetadata(inspireIdentifier)
    if idSuffix is not None:
        data["id"] = int(str(data["id"]) + str(idSuffix))
    if luminosity is not None:
        data["luminosity"] = round(luminosity, 2)
    if collider:
        data["collider"] = collider
    if experiment:
        data["experiment"] = experiment
    if data["experiment"] == "":
        raise RuntimeError(
            "Invalid experiment, set the experiment argument to an appropriate value."
        )
    if data["id"] == 0:
        raise RuntimeError("Invalid id of 0.")
    if data["luminosity"] < 0:
        raise RuntimeError(
            "Invalid luminosity, set the luminosity argument to an appropriate value."
        )
    if data["collider"] == "":
        raise RuntimeError(
            "Invalid collider, set the collider argument to an appropriate value."
        )
    data["source"] = source
    return data


def printable(array, precision=14):
    def transformValue(x):
        return float(np.format_float_scientific(x, precision=precision, unique=False))

    return [transformValue(x) for x in array]


def prodPrefix(productionModes: list):
    prodModes = set(productionModes)
    if len(prodModes) == 1:
        return prodModes.pop()
    elif len(prodModes) == 2 and "HW" in prodModes and "HZ" in prodModes:
        return "HV"
    else:
        return "comb"


def pairDecayPrefix(firstDecay: list, secondDecay: list):
    decs1 = set(firstDecay)
    decs2 = set(secondDecay)
    if len(decs1) == 1 and len(decs2) == 1:
        return decs1.pop() + decs2.pop()
    else:
        return "comb"


def getFilename(prefix: str, collider: str, experiment: str, luminosity: float):
    return "{prefix}_{coll}_{exp}_{lumi}.json".format(
        prefix=prefix,
        coll=collider,
        exp=experiment,
        lumi=round(luminosity),
    )


def writeToJson(filename: str, data: dict):
    from pathlib import Path
    import json

    path = Path(filename)
    if path.is_file():
        validExistingFile = True
        with path.open(mode="r") as exFile:
            try:
                oldDat = json.load(exFile)
            except:
                validExistingFile = False
        if validExistingFile:
            if oldDat["id"] != data["id"]:
                raise RuntimeError("Existing file with different id.")

    print("Looks good. Saving this to", filename)
    with path.open(mode="w") as f:
        json.dump(data, f, indent=4)
    return filename


def implementChannelLimit(
    inspireIdentifier: Union[str, int],
    process: dict,
    source: str,
    df: pd.DataFrame,
    idSuffix: int = None,
    massResolution: dict = None,
    prefix: str = "",
    gridRelevanceTolerance: float = 1e-3,
    luminosity: float = None,
    experiment: str = "",
    constraints: Union[list, dict] = [],
    normalization: dict = {},
    collider: str = "",
    acceptances: list = [],
):
    """
    Implements a channel limit.

    All of the dict arguments correspond directly to the json fields of the same
    name. Their formal definition can be found at
    https://higgsbounds.gitlab.io/higgstools/Datafile.html under
    ChannelLimit.

    Parameters
    ----------
     - inspireIdentifier: A unique identifier for the inspire lookup. Either the
                       arxiv id as a string (should work), or the numeric
                       inspire id (will definitely work).
     - process: the process definition, a dict of 1 field:
            {
                "channels": [one or more channels, where each channel has the form [production mode, decay mode]]
            }
     - source: the exact source of the implemented numbers, Fig/Tab number or
               hepdata link
     - df: the dataframe containing the data, required column names are:
        + "m": mass of the particle in GeV
        + "obs": observed 95% CL limit in pb or as a signal strength (if a
                 normalization is given)
        + "exp": expected 95% CL limit in pb or as a signal strength (if a
                 normalization is given)

    Optional Parameters
    -------------------
    - massResolution: the mass resolutions for the particle. Defaults to best
                      possible resolution, but this should usually always be
                      set.
    - normalization: if this is set, the limit is not set on a rate in pb, but
                     on a signal strangth normalized to the rate of the given
                     reference process
    - constraints: constraints that represent model-assumptions made by the
                   analysis.
    - acceptances: acceptances for the individual channels of the analysis

    The following parameters are used to prevent clashes of limitIds or file
    names:
    - idSuffix: a suffix to add to the automatically determined limit id to
                distinguish multiple limits from the same publication
    - prefix: manually set the prefix of the file name. By default this is
              generated from the production modes.

    These parameters can be used to override the inspire based metadata
    detection if something doesn't work out there:
    - luminosity: manually set the integrated luminosity if the automatic
                  detection fails
    - collider: manually set the collider if the automatic detection fails
    """
    data = assembleMetadata(
        inspireIdentifier, source, idSuffix, luminosity, collider, experiment
    )

    data["limitClass"] = "ChannelLimit"
    data["process"] = process
    if constraints:
        if isinstance(constraints, list):
            data["constraints"] = constraints
        elif isinstance(constraints, dict):
            data["constraints"] = [constraints]
    if normalization:
        data["normalization"] = normalization

    rel = relevantGridEntries(
        df.m, df[["exp", "obs"]].values, rtol=gridRelevanceTolerance
    )
    if len(df.m) > 1:
        plt.figure()
        plt.plot(df.m, df.exp, c="C0")
        plt.plot(df.loc[rel].m, df.loc[rel].exp, c="C0", ls="--")
        plt.plot(df.m, df.obs, c="C1")
        plt.plot(df.loc[rel].m, df.loc[rel].obs, c="C1", ls="--")
        plt.show()

    data["analysis"] = OrderedDict()
    if massResolution is not None:
        data["analysis"]["massResolution"] = massResolution
    if acceptances:
        data["analysis"]["acceptances"] = acceptances
        for acc in acceptances:
            for key, val in acc.items():
                if isinstance(val, np.ndarray):
                    acc[key] = printable(val)
    data["analysis"]["grid"] = {"mass": printable(df.loc[rel].m)}
    data["analysis"]["limit"] = {
        "observed": printable(df.loc[rel].obs),
        "expected": printable(df.loc[rel].exp),
    }

    if prefix == "":
        prefix = prodPrefix([x[0] for x in data["process"]["channels"]])
    filename = getFilename(
        prefix,
        data["collider"],
        data["experiment"],
        data["luminosity"],
    )
    return writeToJson(filename, data)


def implementChannelWidthLimit(
    inspireIdentifier: Union[str, int],
    process: dict,
    source: str,
    df: pd.DataFrame,
    massResolution: dict = None,
    normalization: dict = {},
    constraints: Union[list, dict] = [],
    acceptances: list = [],
    idSuffix: int = None,
    prefix: str = "",
    luminosity: float = None,
    collider: str = "",
):
    """
    Implements a width-dependent channel limit.

    All of the dict arguments correspond directly to the json fields of the same
    name. Their formal definition can be found at
    https://higgsbounds.gitlab.io/higgstools/Datafile.html under
    ChannelWidthLimit.

    Parameters
    ----------
     - inspireIdentifier: A unique identifier for the inspire lookup. Either the
                       arxiv id as a string (should work), or the numeric
                       inspire id (will definitely work).
     - process: the process definition, a dict of 1 field:
            {
                "channels": [one or more channels, where each channel has the form [production mode, decay mode]]
            }
     - source: the exact source of the implemented numbers, Fig/Tab number or
               hepdata link
     - df: the dataframe containing the data, required column names are:
        + "m": mass of the particle in GeV
        + "width" or "normWidth": the total width of the particle, either in GeV
                                  ("width") or normalized to "m" ("normWidth")
        + "obs": observed 95% CL limit in pb or as a signal strength (if a
                 normalization is given)
        + "exp": expected 95% CL limit in pb or as a signal strength (if a
                 normalization is given)

    Optional Parameters
    -------------------
    - massResolution: the mass resolutions for the particle. Defaults to best
                      possible resolution, but this should usually always be
                      set.
    - normalization: if this is set, the limit is not set on a rate in pb, but
                     on a signal strangth normalized to the rate of the given
                     reference process
    - constraints: constraints that represent model-assumptions made by the
                   analysis.
    - acceptances: acceptances for the individual channels of the analysis

    The following parameters are used to prevent clashes of limitIds or file
    names:
    - idSuffix: a suffix to add to the automatically determined limit id to
                distinguish multiple limits from the same publication
    - prefix: manually set the prefix of the file name. By default this is
              generated from the production modes.

    These parameters can be used to override the inspire based metadata
    detection if something doesn't work out there:
    - luminosity: manually set the integrated luminosity if the automatic
                  detection fails
    - collider: manually set the collider if the automatic detection fails
    """
    data = assembleMetadata(inspireIdentifier, source, idSuffix, luminosity, collider)

    data["limitClass"] = "ChannelWidthLimit"
    data["process"] = process
    if constraints:
        if isinstance(constraints, list):
            data["constraints"] = constraints
        elif isinstance(constraints, dict):
            data["constraints"] = [constraints]
    if normalization:
        data["normalization"] = normalization
    normalized = "normWidth" in list(df)
    widthKey = "normWidth" if normalized else "width"
    plt.figure()
    for val, group in df.groupby(widthKey):
        plt.plot(group.m, group.exp, label="exp " + str(val))
        plt.plot(group.m, group.obs, label="obs " + str(val))
    plt.show()

    df.sort_values(["m", widthKey], inplace=True)
    mGrid = np.unique(df.m)
    mGrid = np.array(mGrid, dtype=float)
    wGrid = np.unique(df[widthKey])
    wGrid = np.array(wGrid, dtype=float)
    MM, WW = np.meshgrid(mGrid, wGrid, indexing="ij")
    assert np.allclose(MM.flatten(order="C"), np.array(df.m, dtype=float))
    assert np.allclose(WW.flatten(order="C"), np.array(df[widthKey], dtype=float))

    data["analysis"] = OrderedDict()
    if massResolution is not None:
        data["analysis"]["massResolution"] = massResolution
    if acceptances:
        data["analysis"]["acceptances"] = acceptances
        for acc in acceptances:
            for key, val in acc.items():
                if isinstance(val, np.ndarray):
                    acc[key] = printable(val)

    data["analysis"]["grid"] = {"mass": printable(mGrid), "width": printable(wGrid)}
    data["analysis"]["relativeWidth"] = normalized
    data["analysis"]["limit"] = {
        "observed": printable(df.obs),
        "expected": printable(df.exp),
    }

    if prefix == "":
        prefix = prodPrefix([x[0] for x in data["process"]["channels"]])
    filename = getFilename(
        prefix,
        data["collider"],
        data["experiment"],
        data["luminosity"],
    )
    return writeToJson(filename, data)


def implementChainDecayLimit(
    inspireIdentifier: Union[str, int],
    process: dict,
    source: str,
    df: pd.DataFrame,
    massResolution: dict = None,
    constraints: dict = {},
    productionAcceptances: dict = {},
    idSuffix: int = None,
    prefix: str = "",
    luminosity: float = None,
    collider: str = "",
    contourfPlotArgs: dict = None,
):
    """
    Implements a chain decay Limit, i.e. `H1 -> Z/W H2`.

    All of the dict arguments correspond directly to the json fields of the same
    name. Their formal definition can be found at
    https://higgsbounds.gitlab.io/higgstools/Datafile.html under
    ChainDecayLimit.

    Parameters
    ----------
     - inspireIdentifier: A unique identifier for the inspire lookup. Either the
                       arxiv id as a string (should work), or the numeric
                       inspire id (will definitely work).
     - process: the process definition, a dict of 3 fields:
            {
                "production": [production modes of the mother particle],
                "decay": [decay modes of the dauughter particle],
                "chain": "W" or "Z", specifies the involved SM particle
            }
     - source: the exact source of the implemented numbers, Fig/Tab number or
               hepdata link
     - df: the dataframe containing the data, required column names are
        + "mMother":  mass of the mother particle in GeV
        + "mDaughter": mass of the daughter particle in GeV
        + "obs": observed limit in pb
        + "exp": expected limit in pb

    Optional Parameters
    -------------------
    - massResolution: the mass resolutions for mother and daughter particles.
                      Defaults to best possible resolution, but this should
                      usually always be set.
    - constraints: constraints that represent model-assumptions made by the
                   analysis.
    - productionAcceptances: acceptances for the individual production modes of
      the analysis

    The following parameters are used to prevent clashes of limitIds or file
    names:
    - idSuffix: a suffix to add to the automatically determined limit id to
                distinguish multiple limits from the same publication
    - prefix: manually set the prefix of the file name. By default this is
              generated from the production modes.

    These parameters can be used to override the inspire based metadata
    detection if something doesn't work out there:
    - luminosity: manually set the integrated luminosity if the automatic
                  detection fails
    - collider: manually set the collider if the automatic detection fails

    And these parameters can be used to adjust the data validation plot:
     - contourfPlotArgs: dict of kwargs to pass to the contourf call
    """
    from matplotlib.colors import LogNorm

    data = assembleMetadata(inspireIdentifier, source, idSuffix, luminosity, collider)

    data["limitClass"] = "ChainDecayLimit"
    data["process"] = process
    if constraints:
        data["constraints"] = constraints

    df.sort_values(["mMother", "mDaughter"], inplace=True)
    mMotherGrid = np.unique(df.mMother)
    mDaughterGrid = np.unique(df.mDaughter)
    MM, MD = np.meshgrid(mMotherGrid, mDaughterGrid, indexing="ij")
    assert np.allclose(MM.flatten(order="C"), df.mMother)
    assert np.allclose(MD.flatten(order="C"), df.mDaughter)

    if len(mMotherGrid) > 1 and len(mDaughterGrid) > 1:
        fig, axes = plt.subplots(figsize=(8, 4), ncols=2, sharex=True, sharey=True)
        if contourfPlotArgs is None:
            maxLimVal = np.max(
                (np.max(df[df.exp < 1e5].exp), np.max(df[df.obs < 1e5].obs))
            )
            minLimVal = np.min((np.min(df.exp), np.min(df.obs)))

            contourfPlotArgs = {
                "norm": LogNorm(minLimVal, maxLimVal),
                "levels": np.logspace(np.log10(minLimVal), np.log10(maxLimVal)),
            }
        axes[0].contourf(
            mMotherGrid,
            mDaughterGrid,
            df.exp.values.reshape(MM.shape).T,
            **contourfPlotArgs
        )
        cf = axes[1].contourf(
            mMotherGrid,
            mDaughterGrid,
            df.obs.values.reshape(MM.shape).T,
            **contourfPlotArgs
        )
        axes[0].set_title("exp")
        axes[1].set_title("obs")
        axes[0].set_ylabel("mDaughter [GeV]")
        for ax in axes:
            ax.set_xlabel("mMother [GeV]")
        fig.colorbar(cf, ax=axes, label="rate [pb]")
        plt.show()
    elif len(mDaughterGrid) == 1:
        fig, ax = plt.subplots()
        ax.plot(mMotherGrid, df.exp, label="exp")
        ax.plot(mMotherGrid, df.obs, label="obs")
        ax.set_xlabel("mMother [GeV]")
        ax.set_title("mDaughter = {} GeV".format(mDaughterGrid[0]))
        ax.set_yscale("log")
        plt.show()
    elif len(mMotherGrid) == 1:
        fig, ax = plt.subplots()
        ax.plot(mDaughterGrid, df.exp, label="exp")
        ax.plot(mDaughterGrid, df.obs, label="obs")
        ax.set_xlabel("mDaughter [GeV]")
        ax.set_title("mMother = {} GeV".format(mMotherGrid[0]))
        ax.set_yscale("log")
        plt.show()

    data["analysis"] = OrderedDict()
    if massResolution is not None:
        data["analysis"]["massResolution"] = massResolution
    if productionAcceptances:
        data["analysis"]["productionAcceptances"] = productionAcceptances
        for acc in productionAcceptances:
            for key, val in acc.items():
                if isinstance(val, np.ndarray):
                    acc[key] = printable(val)
    data["analysis"]["grid"] = {
        "massMother": printable(mMotherGrid),
        "massDaughter": printable(mDaughterGrid),
    }
    data["analysis"]["limit"] = {
        "observed": printable(df.obs.to_numpy()),
        "expected": printable(df.exp.to_numpy()),
    }

    if prefix == "":
        prefix = prodPrefix(data["process"]["production"])
    filename = getFilename(
        prefix,
        data["collider"],
        data["experiment"],
        data["luminosity"],
    )
    return writeToJson(filename, data)


def implementPairDecayLimit(
    inspireIdentifier: Union[str, int],
    process: dict,
    source: str,
    df: pd.DataFrame,
    massResolution: dict = None,
    constraints: dict = {},
    idSuffix: int = None,
    prefix: str = "",
    luminosity: float = None,
    collider: str = "",
):
    """
    Implements a pair decay limit, e.g. `p p -> H3 -> H1 H2`.

    All of the dict arguments correspond directly to the json fields of the same
    name. Their formal definition can be found at
    https://higgsbounds.gitlab.io/higgstools/Datafile.html under PairDecayLimit.

    Parameters
    ----------
     - inspireIdentifier: A unique identifier for the inspire lookup. Either the
                       arxiv id as a string (should work), or the numeric
                       inspire id (will definitely work).
     - process: the process definition, a dict of 2 fields:
            {
                "production": [production modes of the mother particle]
                "firstDecay": [decay modes of the first daughter particle],
                "secondDecay": [decay modes of the second daughter particle]
            }
     - source: the exact source of the implemented numbers, Fig/Tab number or
               hepdata link
     - df: the dataframe containing the data, required column names are
        + "mMother":  mass of the mother particle in GeV
        + "mDaughter1": mass of the first daughter particle in GeV
        + "mDaughter2": mass of the second daughter particle in GeV, this is
          optional and if it is not given this indicates that the analysis
          assumes equal masses for the two daughter particles
        + "obs": observed limit in pb
        + "exp": expected limit in pb

    Optional Parameters
    -------------------
    - massResolution: the mass resolutions for the two particles. Defaults to
                      best possible resolution, but this should usually always
                      be set.
    - constraints: constraints that represent model-assumptions made by the
                   analysis.

    The following parameters are used to prevent clashes of limitIds or file
    names:
    - idSuffix: a suffix to add to the automatically determined limit id to
                distinguish multiple limits from the same publication
    - prefix: manually set the prefix of the file name. By default this is
              generated from the production modes.

    These parameters can be used to override the inspire based metadata
    detection if something doesn't work out there:
    - luminosity: manually set the integrated luminosity if the automatic
                  detection fails
    - collider: manually set the collider if the automatic detection fails
    """
    from matplotlib.colors import LogNorm

    data = assembleMetadata(inspireIdentifier, source, idSuffix, luminosity, collider)

    equalDaughterMasses = "mDaughter2" not in list(df)
    if "mDaughter" in list(df) and "mDaughter1" not in list(df):
        df.rename(columns={"mDaughter": "mDaughter1"}, inplace=True)

    data["limitClass"] = "PairDecayLimit"
    data["process"] = process
    if constraints:
        data["constraints"] = constraints

    mMotherGrid = np.unique(df.mMother)
    mDaughter1Grid = np.unique(df.mDaughter1)
    if equalDaughterMasses:
        mDaughter2Grid = mDaughter1Grid
    else:
        mDaughter2Grid = np.unique(df.mDaughter2)

    if equalDaughterMasses:
        df.sort_values(["mMother", "mDaughter1"], inplace=True)
        if len(mMotherGrid) > 1 and len(mDaughter1Grid) > 1:
            MM, MD = np.meshgrid(mMotherGrid, mDaughter1Grid, indexing="ij")
            assert np.allclose(MM.flatten(order="C"), df.mMother)
            assert np.allclose(MD.flatten(order="C"), df.mDaughter1)

            maxLimVal = np.max(
                (np.max(df[df.exp < 1e5].exp), np.max(df[df.obs < 1e5].obs))
            )
            minLimVal = np.min((np.min(df.exp), np.min(df.obs)))
            fig, axes = plt.subplots(figsize=(8, 4), ncols=2, sharex=True, sharey=True)
            args = {
                "norm": LogNorm(minLimVal, maxLimVal),
                "levels": np.logspace(np.log10(minLimVal), np.log10(maxLimVal)),
            }
            axes[0].contourf(
                mMotherGrid, mDaughter1Grid, df.exp.values.reshape(MM.shape).T, **args
            )
            cf = axes[1].contourf(
                mMotherGrid, mDaughter1Grid, df.obs.values.reshape(MM.shape).T, **args
            )
            axes[0].set_title("exp")
            axes[1].set_title("obs")
            axes[0].set_ylabel("mDaughter [GeV]")
            for ax in axes:
                ax.set_xlabel("mMother [GeV]")
            fig.colorbar(cf, ax=axes, label="rate [pb]")
            plt.show()
        elif len(mDaughter1Grid) == 1:
            fig, ax = plt.subplots()
            ax.plot(mMotherGrid, df.exp, label="exp")
            ax.plot(mMotherGrid, df.obs, label="obs")
            ax.set_xlabel("mMother [GeV]")
            ax.set_title("mDaughter = {} GeV".format(mDaughter1Grid[0]))
            ax.set_yscale("log")
            plt.show()
        elif len(mMotherGrid) == 1:
            fig, ax = plt.subplots()
            ax.plot(mDaughter1Grid, df.exp, label="exp")
            ax.plot(mDaughter1Grid, df.obs, label="obs")
            ax.set_xlabel("mDaughter [GeV]")
            ax.set_title("mMother = {} GeV".format(mMotherGrid[0]))
            ax.set_yscale("log")
            plt.show()

    data["analysis"] = OrderedDict()
    if equalDaughterMasses:
        data["analysis"]["equalDaughterMasses"] = True
    if massResolution is not None:
        data["analysis"]["massResolution"] = massResolution
    data["analysis"]["grid"] = {
        "massMother": printable(mMotherGrid),
        "massFirstDaughter": printable(mDaughter1Grid),
    }
    if not equalDaughterMasses:
        data["analysis"]["grid"]["massSecondDaughter"] = printable(mDaughter2Grid)
    data["analysis"]["limit"] = {
        "observed": printable(df.obs),
        "expected": printable(df.exp),
    }

    if prefix == "":
        prefix = prodPrefix(data["process"]["production"])
    filename = getFilename(
        prefix,
        data["collider"],
        data["experiment"],
        data["luminosity"],
    )
    return writeToJson(filename, data)


def implementPairProductionLimit(
    inspireIdentifier: Union[str, int],
    process: dict,
    source: str,
    df: pd.DataFrame,
    massResolution: dict = None,
    constraints: dict = {},
    idSuffix: int = None,
    prefix: str = "",
    luminosity: float = None,
    collider: str = "",
):
    """
    Implements a non-resonant pair production Limit, e.g. `p p -> H1 H2`.

    All of the dict arguments correspond directly to the json fields of the same
    name. Their formal definition can be found at
    https://higgsbounds.gitlab.io/higgstools/Datafile.html under
    PairProductionLimit.

    Parameters
    ----------
     - inspireIdentifier: A unique identifier for the inspire lookup. Either the
                       arxiv id as a string (should work), or the numeric
                       inspire id (will definitely work).
     - process: the process definition, a dict of 2 fields:
            {
                "firstDecay": [decay modes of the first particle],
                "secondDecay": [decay modes of the second particle]
            }
     - source: the exact source of the implemented numbers, Fig/Tab number or
               hepdata link
     - df: the dataframe containing the data, required column names are
        + "m1":  mass of the first particle in GeV
        + "m2": mass of the second particle in GeV, this is optional and if it
          is not given this indicates that the analysis assumes equal masses for
          the two particles
        + "obs": observed limit in pb
        + "exp": expected limit in pb


    Optional Parameters
    -------------------
    - massResolution: the mass resolutions for the two particles. Defaults to
                      best possible resolution, but this should usually always
                      be set.
    - constraints: constraints that represent model-assumptions made by the
                   analysis.

    The following parameters are used to prevent clashes of limitIds or file
    names:
    - idSuffix: a suffix to add to the automatically determined limit id to
                distinguish multiple limits from the same publication
    - prefix: manually set the prefix of the file name. By default this is
              generated from the production modes.

    These parameters can be used to override the inspire based metadata
    detection if something doesn't work out there:
    - luminosity: manually set the integrated luminosity if the automatic
                  detection fails
    - collider: manually set the collider if the automatic detection fails
    """
    from matplotlib.colors import LogNorm

    data = assembleMetadata(inspireIdentifier, source, idSuffix, luminosity, collider)

    equalParticleMasses = "m2" not in list(df)
    if "m" in list(df) and "m1" not in list(df):
        df.rename(columns={"m": "m1"}, inplace=True)

    data["limitClass"] = "PairProductionLimit"
    data["process"] = process
    if constraints:
        data["constraints"] = constraints

    m1Grid = np.unique(df.m1)
    if equalParticleMasses:
        m2Grid = m1Grid
    else:
        m2Grid = np.unique(df.m2)

    if equalParticleMasses:
        df.sort_values("m1", inplace=True)
        if len(m1Grid) > 1:
            fig, ax = plt.subplots()
            ax.plot(m1Grid, df.exp, label="exp")
            ax.plot(m1Grid, df.obs, label="obs")
            ax.set_xlabel("m [GeV]")
            ax.set_yscale("log")
            plt.show()
    else:
        df.sort_values(["m1", "m2"], inplace=True)
        MM1, MM2 = np.meshgrid(m1Grid, m2Grid, indexing="ij")
        assert np.allclose(MM1.flatten(order="C"), df.m1)
        assert np.allclose(MM2.flatten(order="C"), df.m2)
        maxLimVal = np.max((np.max(df[df.exp < 1e5].exp), np.max(df[df.obs < 1e5].obs)))
        minLimVal = np.min((np.min(df.exp), np.min(df.obs)))
        fig, axes = plt.subplots(figsize=(8, 4), ncols=2, sharex=True, sharey=True)
        args = {
            "norm": LogNorm(minLimVal, maxLimVal),
            "levels": np.logspace(np.log10(minLimVal), np.log10(maxLimVal)),
        }
        axes[0].contourf(m1Grid, m2Grid, df.exp.values.reshape(MM1.shape).T, **args)
        cf = axes[1].contourf(
            m1Grid, m2Grid, df.obs.values.reshape(MM1.shape).T, **args
        )
        axes[0].set_title("exp")
        axes[1].set_title("obs")
        axes[0].set_ylabel("m_2 [GeV]")
        for ax in axes:
            ax.set_xlabel("m_1 [GeV]")
        fig.colorbar(cf, ax=axes, label="rate [pb]")
        plt.show()

    data["analysis"] = OrderedDict()
    if equalParticleMasses:
        data["analysis"]["equalParticleMasses"] = True
    if massResolution is not None:
        data["analysis"]["massResolution"] = massResolution
    data["analysis"]["grid"] = {
        "massFirstParticle": printable(m1Grid),
    }
    if not equalParticleMasses:
        data["analysis"]["grid"]["massSecondParticle"] = printable(m2Grid)
    data["analysis"]["limit"] = {
        "observed": printable(df.obs),
        "expected": printable(df.exp),
    }

    if prefix == "":
        prefix = "nonres"
    filename = getFilename(
        prefix,
        data["collider"],
        data["experiment"],
        data["luminosity"],
    )
    return writeToJson(filename, data)


def writeCorrelationMatrix(corrMat: pd.DataFrame, requiredKeys: list):
    """Convert a correlation matrix into a minimal dict format that only
    contains the needed non-zero triangular elements."""
    assert set(requiredKeys) == set(
        corrMat
    ), "names in correlation matrix have to match the sub-measurements"
    assert set(corrMat.index) == set(
        corrMat
    ), "correlation matrix has to include all bin names"
    corrMat = corrMat[corrMat.index]
    assert np.allclose(corrMat.to_numpy().T, corrMat.to_numpy()) or np.allclose(
        corrMat.to_numpy(), np.triu(corrMat.to_numpy())
    ), "correlation matrix has to be symmetric or upper triangular"
    names = list(corrMat)
    return {
        names[i]: corrMat.iloc[i, i + 1 :][corrMat.iloc[i, i + 1 :] != 0].to_dict()
        for i in range(len(names) - 1)
        if np.any(corrMat.iloc[i, i + 1 :] != 0)
    }


def implementMeasurement(
    inspireIdentifier: Union[str, int],
    source: str,
    subMeasurements: dict,
    referenceMass: float,
    massResolution: float,
    prefix,
    corrMatExp: pd.DataFrame = None,
    corrMatTheo: pd.DataFrame = None,
    luminosity: float = None,
    experiment: str = "",
    collider: str = "",
    referenceModel: str = "SMHiggsEW",
):
    """
    Implement a HiggsSignals measurement.

    Parameters
    ----------
     - inspireIdentifier: A unique identifier for the inspire lookup. Either the
                       arxiv id as a string (should work), or the numeric
                       inspire id (will definitely work).
     - source: the exact source of the implemented numbers, Fig/Tab number or
               hepdata link
     - subMeasurements: a dict of the sub-measurements, this is directly passed
       into the resulting json. See the datafiles documentation for the
       definition of the fields in the different types of sub-measurements.
     - referenceMass: the reference mass assumed in the measurement
     - massResolution: the (estimated) mass resolution of the measurement used
       for assignment as an absolute +- uncertainty around the referenceMass
     - prefix: the filename prefix of the created json file, should specify the
       decay and/or production mode of the measurement

    Optional Parameters
    -------------------
     - corrMatExp: the correlation matrix of the experimental errors in a
       DataFrame format where both index and column names are equal to the
       subMeasurement keys
     - corrMatTheo: the correlation matrix of the (reference model) theoretical
       errors in a DataFrame format where both index and column names are equal
       to the subMeasurement keys
     - collider: manually set the collider if the automatic detection fails
     - luminosity: manually set the luminosity if automatic detection fails
     - experiment: manually set the experiment if automatic detection fails
     - referenceModel: set the reference model used for the signal strength
       normalization in this measurement
    """
    data = assembleMetadata(
        inspireIdentifier, source, None, luminosity, collider, experiment
    )
    del data["limitClass"]
    data["referenceMass"] = referenceMass
    data["referenceModel"] = referenceModel
    data["massResolution"] = massResolution
    data["subMeasurements"] = subMeasurements
    for name, bin in data["subMeasurements"].items():
        if "obs" in bin:
            data["subMeasurements"][name]["obs"] = printable(bin["obs"])
        if "obsMass" in bin:
            data["subMeasurements"][name]["obsMass"] = printable(bin["obsMass"])
        if "ref" in bin:
            data["subMeasurements"][name]["ref"] = printable(bin["ref"])
        if "channelWeights" in bin:
            data["subMeasurements"][name]["channelWeights"] = printable(
                bin["channelWeights"]
            )

    if corrMatExp is not None or corrMatTheo is not None:
        data["correlations"] = {}

    if corrMatExp is not None:
        data["correlations"]["experimental"] = writeCorrelationMatrix(
            corrMatExp, data["subMeasurements"].keys()
        )

    if corrMatTheo is not None:
        data["correlations"]["theoretical"] = writeCorrelationMatrix(
            corrMatTheo, data["subMeasurements"].keys()
        )

    filename = getFilename(
        prefix,
        data["collider"],
        data["experiment"],
        data["luminosity"],
    )
    return writeToJson(filename, data)
