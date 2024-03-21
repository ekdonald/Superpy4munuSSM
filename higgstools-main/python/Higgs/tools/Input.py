import pandas as pd
from itertools import product, combinations_with_replacement, chain
import numpy as np
import Higgs.predictions as HP
from typing import Union


def predictionsFromDict(
    data: dict,
    neutralIds: list,
    singleChargedIds: list,
    doubleChargedIds: list,
    referenceModels: dict = {},
    massPrefix = "m",
    massUncPrefix = "dm",
    widthPrefix = "w",
    CpPrefix = "CP",
    brPrefix = "br",
    cxnPrefix = "cxn",
    pairCxnPrefix = "paircxn",
    effCPrefix = "effc",
    oddCoupSuffix = "p",
    normCxnSuffix = "norm",
    useExplicitBr: Union[bool, dict] = True,
    calcggH=True,
    calcHgamgam=True
):
    """
    Construct a Higgs.predictions.Predictions object from a set of model
    predictions given in the form of a dictionary of cxns, BRs and couplings.

    If effective couplings of the neutral scalars are provided, the effective
    coupling input is used first, and then any explicitely specified cross
    sections and branching ratios are used to overwrite the values obtained from
    the effective coupling approximation.

    Arguments
    ----------
      - data: the dictionary containing the data. For a particle with ID `"ID"`
        the following quantities can be specified:
            - `m_ID`: **mass of the particle**
            - `dm_ID`: theoretical mass uncertainty
            - `w_ID`: **total width of the particle**
            - `br_ID_XX`: BR(ID->XX) for `"XX" in Higgs.predictions.Decay`
            - `br_ID_ID2_ID3`: BR(ID-> ID2 ID3) where `ID2` and `ID3` are either
              IDs of other particles, or members of
              `Higgs.predictions.ChainDecay`
            - `cxn_ID_XX_COLL`: a production cross section in the
              `Higgs.predictions.Production` mode `"XX"` at the
              `Higgs.predictions.Collider` `"COLL"`
            - `cxn_id_XX_COLL_norm`: a normalized production cross section in
              the `Higgs.predictions.Production` mode `"XX"` at the
              `Higgs.predictions.Collider` `"COLL"`
        Additional recognized fields for neutral particles only are:
            - `CP_ID`: **CP quantum number (only for neutral particles)**
            - `effc_ID_ff_s`: scalar part of the SM-normalized effective
              coupling to the fermions `ff`
            - `effc_ID_ff_p`: pseudoscalar part of the SM-normalized effective
              coupling to the fermions `ff`
            - `effC_ID_VV`: SM-normalized effective coupling to the gauge bosons
              `VV`

        The **highlighted entries** are required. Finally, non-resonant pair
        production cross sections of the particles `ID1` and `ID2` at a
        `Higgs.predictions.Collider` `COLL` can be specified as
        `paircxn_ID1_ID2_COLL`.

     - neutralIds: IDs of the neutral particles
     - singleChargedIds: IDs of the singly charged particles
     - doublyChargedIds: IDs of the doubly charged particles

     Optional Arguments
     ------------------
      - referenceModels: a dictionary that can be used to assign a
        `Higgs.predictions.ReferenceModel` to any particle. If none is specified
        this defaults to `Higgs.predictions.SMHiggsInterp`.
      - useExplicitBr: whether the explicitly given BRs should be used. 
        Alternatively, the BRs are calculated based on the given effective 
        couplings. This can also be controlled for individual Higgs bosons and/or
        decay channels by providing a dictionary (e.g.
        `useExplicitBr = {'h1': False, 'h2': {'WW': False}}`)
     The remaining optional arguments can be used to adjust the prefixes and
     suffixes that make up the keys as described above.
    """

    def effcFor(hj, point):
        realVals = {
            k.split("_")[2]: v
            for k, v in point.items()
            if k.startswith(f"{effCPrefix}_{hj}_")
            and not k.endswith(f"_{oddCoupSuffix}")
        }
        imagVals = {
            k.split("_")[2]: v
            for k, v in point.items()
            if k.startswith(f"{effCPrefix}_{hj}_") and k.endswith(f"_{oddCoupSuffix}")
        }

        return HP.NeutralEffectiveCouplings(
            **{
                k: realVals.get(k, 0.0) + 1j * imagVals[k]
                if k in imagVals
                else realVals.get(k, 0.0)
                for k in np.unique(list(chain(realVals.keys(), imagVals.keys())))
            }
        )

    pred = HP.Predictions()
    for hj in neutralIds:
        try:
            h = pred.addParticle(
                    HP.BsmParticle(
                        hj, charge="neutral", cp=HP.CP(int(data[f"{CpPrefix}_{hj}"]))
                    )
            )
        except KeyError:
            raise ValueError(f"no CP information given for particle {hj}")
        try:
            h.setMass(data[f"{massPrefix}_{hj}"])
        except KeyError:
            raise ValueError(f"no mass given for particle {hj}")
        h.setMassUnc(data.get(f"{massUncPrefix}_{hj}", 0))

        HP.effectiveCouplingInput(
            h,
            effcFor(hj, data),
            referenceModels.get(
                hj,
                HP.ReferenceModel.SMHiggsInterp
            ),
            calcggH=calcggH,
            calcHgamgam=calcHgamgam
        )

        try:
            gamTot = data[f"{widthPrefix}_{hj}"]
        except KeyError:
            raise ValueError(f"no total width given for particle {hj}")

        def set_brs_from_data():
            for br in (b for b in data.keys() if b.startswith(f"{brPrefix}_{hj}_")):
                        finalState = br.split("_")[2:]
                        h.setDecayWidth(*finalState, data[br] * gamTot)

        if type(useExplicitBr) == dict:
            if hj in useExplicitBr.keys():
                if type(useExplicitBr[hj]) == dict:
                    for br in (b for b in data.keys() if b.startswith(f"{brPrefix}_{hj}_")):
                        setExplicitBr = True
                        finalState = br.split("_")[2:]
                        if len(finalState) == 1:
                            if finalState[0] in useExplicitBr[hj].keys():
                                if not useExplicitBr[hj][finalState[0]]: 
                                    setExplicitBr = False
                        if setExplicitBr:
                            h.setDecayWidth(*finalState, data[br] * gamTot)
                elif type(useExplicitBr[hj]) == bool:
                    if useExplicitBr[hj]:
                        set_brs_from_data()
                else:
                    raise ValueError('useExplicitBr does not have suitable format.')
            else:
                set_brs_from_data()
        elif type(useExplicitBr) == bool:
            if useExplicitBr:
                set_brs_from_data()
        else:
            raise ValueError('useExplicitBr does not have suitable format.')


        for cxn in (x for x in data.keys() if x.startswith(f"{cxnPrefix}_{hj}_")):
            chan, coll = cxn.split("_")[2:4]
            if cxn.endswith(f"_{normCxnSuffix}"):
                h.setNormalizedCxn(
                    coll,
                    chan,
                    data[cxn],
                    referenceModels.get(hj, HP.ReferenceModel.SMHiggs if h.mass() > 150 else HP.ReferenceModel.SMHiggsEW),
                )
            else:
                h.setCxn(coll, chan, data[cxn])

        if gamTot > h.totalWidth():
            h.setDecayWidth("unknown", "unknown", gamTot - h.totalWidth())

    for hpj in chain(singleChargedIds, doubleChargedIds):
        if hpj in singleChargedIds:
            h = pred.addParticle(HP.BsmParticle(hpj, charge="single", cp="undefined"))
        else:
            h = pred.addParticle(HP.BsmParticle(hpj, charge="doubly", cp="undefined"))
        h.setMass(data[f"{massPrefix}_{hpj}"])
        h.setMassUnc(data.get(f"{massUncPrefix}_{hpj}", 0))
        h.setTotalWidth(data[f"{widthPrefix}_{hpj}"])
        for br in (b for b in data.keys() if b.startswith(f"{brPrefix}_{hpj}_")):
            finalState = br.split("_")[2:]
            h.setBr(*finalState, data[br])
        for cxn in (x for x in data.keys() if x.startswith(f"{cxnPrefix}_{hpj}_")):
            chan, coll = cxn.split("_")[2:4]
            h.setCxn(coll, chan, data[cxn])

    for paircxn in data.keys():
        if paircxn.startswith(f"{pairCxnPrefix}_"):
            p1, p2, coll = paircxn.split("_")[1:]
            pred.setBsmPairCxn(coll, p1, p2, data[paircxn])
    return pred


def getCPFromEffC(
    p: str,
    data: dict,
    effCPrefix: str = "effc",
    evenCoupSuffix: str = "s",
    oddCoupSuffix: str = "p",
):
    """
    Use the fermionic effective couplings of a particle with ID p that are
    stored in the data dictionary to determine its CP quantum number.
    """
    coups = {k: v for k, v in data.items() if k.startswith(f"{effCPrefix}_{p}_")}
    evenCoups = [v for k, v in coups.items() if k.endswith(f"_{evenCoupSuffix}")]
    oddCoups = [v for k, v in coups.items() if k.endswith(f"_{oddCoupSuffix}")]
    return np.allclose(oddCoups, 0) - np.allclose(evenCoups, 0)


def readHB5SLHA(
    file: str,
    neutralPDGs: list,
    chargedPDGs: list,
    invisibleWidthThreshold: float = 1e-10,
    invisiblePDGs: list = [],
):
    """
    Reads a HiggsBounds-5 compatible SLHA file into a dictionary. The SLHA file
    must contain the blocks HiggsCouplingsBosons and HiggsCouplingsFermions that
    contain the effective couplings of the neutral Higgs bosons as defined in
    the HiggsBounds-5 manual. Additionally, the masses of all considered
    particles have to be present in the MASS block, and a DECAY block has to
    exist for every considered particle and for the top-quark.

    **Uses and required the pylha package.**

    This function is slightly more general than the HB5 SLHA interface, in that
    it supports models other than the MSSM and NMSSM, by using the PDG numbers
    of the neutral and charged particles as input.

    Arguments
    ---------
      - file: filename of the SLHA file to read
      - neutralPDGs: list of the PDG numbers for neutral spin-0 particles
      - chargedPDGs: list of the PDG numbers for charged spin-0 particles

    Optional Arguments
    ------------------
      - invisibleWidthThreshold: any particles with total widths less than this
        quantity (in GeV) are treated as invisible (does not apply for the 
        SM quarks and charged leptons; only used if `inivisible PDGs = []`)
      - invisiblePDGs: list of the PDG numbers for invisible particles 
        If `invisiblePDGs = []`, which is the default, invisibleWidthThreshold is used.

    Returns
    -------
      A dictionary containing all of the data relevant for HiggsTools that could
      be extracted from the SLHA file. See predictionsFromData for a detailed
      description of the dictionary keys.

    """
    import pylha
    from collections import defaultdict

    with open(
        file,
        "r",
    ) as slhafile:
        lines = slhafile.read()
        # avoids issues with SPheno which sometimes outputs a DECAY1L block which is not correctly parsed by pylha
        lines = lines.replace('DECAY1L', 'DECAYOL')
        # avoids issues with FeynHiggs which output all block names with capital letters
        lines = lines.replace('HIGGSCOUPLINGSBOSONS', 'HiggsCouplingsBosons')
        lines = lines.replace('HIGGSCOUPLINGSFERMIONS', 'HiggsCouplingsFermions')
        slha = pylha.load(lines)

    if not "HiggsCouplingsBosons" in slha['BLOCK'].keys() or \
    not "HiggsCouplingsFermions" in slha['BLOCK'].keys() or \
    not "MASS" in slha['BLOCK'].keys() or \
    not "DECAY" in slha.keys():
        raise ValueError('SLHA file needs to contain the blocks: MASS, DECAY, HiggsCouplingsBosons, HiggsCouplingsFermions')

    fermions = {
        **{v + 1: k for v, k in enumerate(["d", "u", "s", "c", "b", "t"])},
        **{11 + v: k for v, k in enumerate(["e", "nu", "mu", "nu", "tau", "nu"])},
    }
    bosons = {21 + v: k for v, k in enumerate(["g", "gam", "Z", "W"])}

    def parseCouplings():
        HiggsCouplingsBosons = slha["BLOCK"]["HiggsCouplingsBosons"]["values"]
        coupsVV = defaultdict(
            float, {tuple(sorted(v[-(v[1]) :])): v[0] for v in HiggsCouplingsBosons}
        )

        HiggsCouplingsFermions = slha["BLOCK"]["HiggsCouplingsFermions"]["values"]
        coupsffEven = defaultdict(
            float, {tuple(sorted(v[-3:])): v[0] for v in HiggsCouplingsFermions}
        )
        coupsffOdd = defaultdict(
            float, {tuple(sorted(v[-3:])): v[1] for v in HiggsCouplingsFermions}
        )
        coups = dict()
        for pdg in neutralPDGs:
            coups.update(
                (f"effc_{pdg}_{k}{k}", coupsVV[tuple(sorted([v, v, pdg]))])
                for v, k in bosons.items()
            )
            coups.update(
                (f"effc_{pdg}_{k}{k}_s", coupsffEven[tuple(sorted([v, v, pdg]))])
                for v, k in fermions.items()
                if "nu" not in k
            )
            coups.update(
                (f"effc_{pdg}_{k}{k}_p", coupsffOdd[tuple(sorted([v, v, pdg]))])
                for v, k in fermions.items()
                if "nu" not in k
            )
            coups[f"CP_{pdg}"] = getCPFromEffC(str(pdg), coups)
        coups.update(
            (
                f"paircxn_{hi}_{hj}_LEP",
                coupsVV[tuple(sorted([23, hi, hj]))] ** 2,
            )
            for hi, hj in combinations_with_replacement(neutralPDGs, 2)
        )
        return coups

    def parseMasses():
        massVals = dict(slha["BLOCK"]["MASS"]["values"])
        return {f"m_{pdg}": massVals[pdg] for pdg in chain(neutralPDGs, chargedPDGs)}

    def parseMassUncs():
        try:
            massUncVals = dict(slha["BLOCK"]["DMASS"]["values"])
            return {
                f"dm_{pdg}": massUncVals.get(pdg, 0.0)
                for pdg in chain(neutralPDGs, chargedPDGs)
            }
        except KeyError:
            return dict()

    def parseDecays():
        try:
            widths = {
                f"w_{pdg}": slha["DECAY"][str(pdg)]["info"][0]
                for pdg in chain(neutralPDGs, chargedPDGs)
            }
        except:
            raise ValueError('SLHA file needs to contain DECAY block for each spin-0 particles defined via neutralIDs and chargedIDs')
        decays = {}
        for p in chain(neutralPDGs, chargedPDGs):
            SMparticles = {**fermions, **bosons}

            decayDat = [
                [*vals[-2:], vals[0]]
                for vals in slha["DECAY"][str(p)]["values"]
                if vals[1] == 2 and len(vals) == 4
            ]

            def toDecay(d1, d2):
                d1 = abs(d1)
                d2 = abs(d2)
                if d1 in SMparticles and d2 in SMparticles:
                    if f"{SMparticles[d1]}{SMparticles[d2]}" in HP.Decay.__members__:
                        return f"br_{p}_{SMparticles[d1]}{SMparticles[d2]}"
                    elif f"{SMparticles[d2]}{SMparticles[d1]}" in HP.Decay.__members__:
                        return f"br_{p}_{SMparticles[d2]}{SMparticles[d1]}"
                    else:
                        return f"unknown_br_{p}_{SMparticles[d1]}_{SMparticles[d2]}"
                elif d1 in SMparticles:
                    return f"br_{p}_{SMparticles[d1]}_{d2}"
                elif d2 in SMparticles:
                    return f"br_{p}_{SMparticles[d2]}_{d1}"
                else:
                    return f"br_{p}_{d1}_{d2}"

            def isInvisible(d1, d2):
                if invisiblePDGs == []:
                    invisibleParticles = [
                        int(k)
                        for k, val in slha["DECAY"].items()
                        if val["info"][0] < invisibleWidthThreshold
                        and (abs(int(k)) not in [1, 2, 3, 4, 5, 6, 11, 13, 15])
                        # last line avoids marking quarks and charged leptons as invisible particles 
                        # this is relevant for SPheno which sets the width of the light leptons and quarks to zero by default
                    ] + [12, 14, 16]
                else:
                    invisibleParticles = invisiblePDGs
                return abs(d1) in invisibleParticles and abs(d2) in invisibleParticles

            decays.update(
                (toDecay(d1, d2), v)
                for d1, d2, v in decayDat
                if not isInvisible(d1, d2)
            )

            if p in neutralPDGs:
                decays[f"br_{p}_directInv"] = sum(
                    v for d1, d2, v in decayDat if isInvisible(d1, d2)
                )

        return {**widths, **decays}

    def parseTopDecays():
        try:
            topDecays = {
                np.max(np.abs(x[-2:])): x[0]
                for x in slha["DECAY"]["6"]["values"]
                if len(x) == 4 and x[1] == 2
            }
        except:
            raise ValueError('SLHA file needs to contain DECAY block for the top quark')

        return {f"cxn_{p}_brtHpb_LHC13": topDecays.get(p, 0) for p in chargedPDGs}

    def parseChargedHiggsCxns():
        hccxns = dict()
        for coll in "LHC8", "LHC13":
            try:
                cxndat = {
                    (tuple(sorted(np.abs([p1, p2])) + [hp])): val
                    for p1, p2, hp, val in slha["BLOCK"][f"ChargedHiggs{coll}"][
                        "values"
                    ]
                }
            except KeyError:
                continue

            for p in chargedPDGs:
                hccxns[f"cxn_Hpmtb_{p}_{coll}"] = cxndat.get((5, 6, p), 0.0)
                hccxns[f"cxn_qqHpm_{p}_{coll}"] = sum(
                    cxndat.get((i, j, p), 0)
                    for i, j in product([1, 2, 3], [2, 3, 4])
                    if i % 2 != j % 2
                )
                hccxns[f"cxn_vbfHpm_{p}_{coll}"] = cxndat.get((1, 1, p), 0.0)
                hccxns[f"cxn_HpmW_{p}_{coll}"] = cxndat.get((0, 24, p), 0.0)
                hccxns[f"cxn_HpmZ_{p}_{coll}"] = cxndat.get((0, 23, p), 0.0)

            for p1, p2 in combinations_with_replacement(chargedPDGs, 2):
                if (0, p1, p2) in cxndat:
                    hccxns[f"paircxn_{p1}_{p2}_{coll}"] = cxndat[(0, p1, p2)]
            for p1, p2 in product(neutralPDGs, chargedPDGs):
                if (0, p1, p2) in cxndat:
                    hccxns[f"paircxn_{p1}_{p2}_{coll}"] = cxndat[(0, p1, p2)]
        return hccxns

    return {
        **parseMasses(),
        **parseMassUncs(),
        **parseDecays(),
        **parseCouplings(),
        **parseTopDecays(),
        **parseChargedHiggsCxns(),
    }


def readHB5Datafiles(prefix: str, neutralIds: list, chargedIds: list):
    """
    Reads input given in the form of the HiggsBounds-5 datafiles into a pandas
    dataframe. Each row of the dataframe can then be used with
    predictionsFromDict as `df.apply(predictionsFromDict, axis=1)`.

    If both effective coupling and hadronic input data is available for the
    given prefix, the effective coupling input is used first, and then any
    explicitely provided cxns and BRs are overwritten.


    Arguments
    ---------
      - prefix: the prefix path to the datafiles, as defined in the
        HiggsBounds-4/5 manuals
      - neutralIds: the IDs to use for the neutral particles, given in the order
        they appear in the datafiles
      - chargedIds: the IDs to use for the charged particles, given in the order
        they appear in the datafiles

    Returns
    -------
      A pandas DataFrame where each row corresponds to a datapoint (row) given
      in the datafiles. See predictionsFromDict for a detailed description of
      the column names.

    """
    try:
        with open(prefix + "MH_GammaTot.dat", "r") as infile:
            nHzero = int((len(infile.readline().split()) - 1) / 2)
    except:
        nHzero = 0
    try:
        with open(prefix + "MHplus_GammaTot.dat", "r") as infile:
            nHplus = int((len(infile.readline().split()) - 1) / 2)
    except:
        nHplus = 0
    if nHzero == 0 and nHplus == 0:
        raise FileNotFoundError(
            f"Neither {prefix}MH_GammaTot.dat nor {prefix}MHplus_GammaTot.dat exist"
        )
    if len(neutralIds) != nHzero:
        raise RuntimeError(
            f"Provided {len(neutralIds)} neutral IDs, but datafiles imply {nHzero} neutral scalars"
        )
    if len(chargedIds) != nHplus:
        raise RuntimeError(
            f"Provided {len(chargedIds)} charged IDs, but datafiles imply {nHplus} charged scalars"
        )

    def symmetricKeys(names):
        return reversed(list(combinations_with_replacement(reversed(names), 2)))

    neut_ff_chans = ["ss", "cc", "bb", "tt", "mumu", "tautau"]
    neut_VV_chans = ["WW", "ZZ", "Zgam", "gamgam", "gg"]

    neut_had_cxns = [
        f"cxn_{hj}_{chan}_{{}}_norm"
        for chan, hj in product(
            (
                "H",
                "ggH",
                "bbH",
                "HW",
                "unused-HZ",
                "vbfH",
                "Htt",
                "tchanHt",
                "schanHt",
                "HtW",
                "qqHZ",
                "ggHZ",
            ),
            neutralIds,
        )
    ]

    charged_had_cxns = (
        [
            f"cxn_{hpj}_{chan}_{{}}"
            for chan, hpj in product(
                (
                    "Hpmtb",
                    "unused-1",
                    "unused-2",
                    "unused-3",
                    "qqHpm",
                    "HpmW",
                    "HpmZ",
                    "vbfHpm",
                ),
                chargedIds,
            )
        ]
        + [f"paircxn_{hpj}_{hpj}_{{}}" for hpj in chargedIds]
        + [f"paircxn_{hpj}_{hi}_{{}}" for hpj, hi in product(chargedIds, neutralIds)]
    )

    keys = {
        "MH_GammaTot": [f"{x}_{hj}" for x, hj in product(("m", "w"), neutralIds)],
        "BR_H_NP": [f"br_{H}_directInv" for H in neutralIds]
        + [
            f"br_{hk}_{hj}_{hi}"
            for hk, (hj, hi) in product(neutralIds, symmetricKeys(neutralIds))
            if hk != hj and hk != hi
        ]
        + [f"br_{hj}_Z_{hi}" for hj, hi in product(neutralIds, repeat=2) if hi != hj]
        + [f"br_{h}_{lfv}" for lfv, h in product(("emu", "etau", "mutau"), neutralIds)]
        + [f"br_{hj}_W_{hpi}" for hj, hpi in product(neutralIds, chargedIds)],
        "effC": [
            key.format(H)
            for key, H in product(
                chain(
                    *(
                        (f"effc_{{}}_{ff}_s", f"effc_{{}}_{ff}_p")
                        for ff in neut_ff_chans
                    ),
                    (f"effc_{{}}_{VV}" for VV in neut_VV_chans),
                ),
                neutralIds,
            )
        ]
        + [f"g{Hi}{Hj}Z" for Hi, Hj in symmetricKeys(neutralIds)],
        "LHC13_1H_hadCS_ratios": [x.format("LHC13") for x in neut_had_cxns],
        "LHC8_1H_hadCS_ratios": [x.format("LHC8") for x in neut_had_cxns],
        "CP_values": [f"CP_{name}" for name in neutralIds],
        "BR_H_OP": [
            f"br_{hi}_{x}"
            for x, hi in product(chain(neut_ff_chans, neut_VV_chans), neutralIds)
        ],
        "LEP_HZ_CS_ratios": [f"cxn_{hi}_eeHZ_LEP_norm" for hi in neutralIds],
        "LEP_H_ff_CS_ratios": [
            f"cxn_{hi}_eeH{ff}_LEP_norm"
            for ff, hi in product(("bb", "tautau"), neutralIds)
        ],
        "LEP_2H_CS_ratios": [
            f"paircxn_{hj}_{hi}_LEP" for hj, hi in symmetricKeys(neutralIds)
        ],
        "MHplus_GammaTot": [f"{x}_{hj}" for x, hj in product(("m", "w"), chargedIds)],
        "BR_t": ["BRt_Wpb"] + [f"cxn_{hp}_brtHpb_LHC13" for hp in chargedIds],
        "BR_Hplus": [
            f"br_{hp}_{x}"
            for x, hp in product(("cs", "cb", "taunu", "tb", "WZ"), chargedIds)
        ]
        + [f"br_{hp}_W_{hi}" for hp, hi in product(chargedIds, neutralIds)],
        "LEP_HpHm_CS_ratios": [f"paircxn_{hpmj}_{hpmj}_LEP" for hpmj in chargedIds],
        "LHC8_Hplus_hadCS": [x.format("LHC8") for x in charged_had_cxns],
        "LHC13_Hplus_hadCS": [x.format("LHC13") for x in charged_had_cxns],
        "MHall_uncertainties": [f"dm_{h}" for h in chain(neutralIds, chargedIds)],
    }

    def readDatafile(name, **kwargs):
        filename = prefix + name + ".dat"
        return pd.read_csv(
            filename, names=keys[name], sep=r"\s+", index_col=0, **kwargs
        )

    df = pd.DataFrame()
    if neutralIds:
        df = pd.concat(
            (df, readDatafile("MH_GammaTot"), readDatafile("BR_H_NP")),
            axis=1,
        )

        effC = True
        try:
            df = df.join(readDatafile("effC"))
        except:
            print(
                f"WARNING: Could not find {prefix}effC.dat, falling back to hadronic input."
            )
            effC = False
            df = pd.concat(
                (
                    df,
                    readDatafile("LHC8_1H_hadCS_ratios"),
                    readDatafile("LHC13_1H_hadCS_ratios"),
                ),
                axis=1,
            )
        try:
            df = df.join(readDatafile("CP_values", dtype=int))
        except:
            if effC:
                for hj in neutralIds:
                    df[f"CP_{hj}"] = df.apply(
                        lambda dat: getCPFromEffC(hj, dat), axis=1
                    )
            else:
                raise

        def effCOptionalRead(name):
            try:
                return readDatafile(name)
            except:
                if effC:
                    return pd.DataFrame()
                else:
                    raise

        df = pd.concat(
            (
                df,
                effCOptionalRead("BR_H_OP"),
                effCOptionalRead("LEP_HZ_CS_ratios"),
                effCOptionalRead("LEP_H_ff_CS_ratios"),
            ),
            axis=1,
        )

        try:
            df = df.join(readDatafile("LEP_2H_CS_ratios"))
        except:
            if effC:
                coupToCxn = {
                    f"g{Hi}{Hj}Z": f"paircxn_{Hi}_{Hj}_LEP"
                    for Hi, Hj in symmetricKeys(neutralIds)
                }
                for coup, cxn in coupToCxn.items():
                    df[cxn] = df[coup] ** 2
            else:
                raise

    if chargedIds:
        df = pd.concat(
            (
                df,
                readDatafile("MHplus_GammaTot"),
                readDatafile("BR_t"),
                readDatafile("BR_Hplus"),
                readDatafile("LEP_HpHm_CS_ratios"),
                readDatafile("LHC8_Hplus_hadCS"),
                readDatafile("LHC13_Hplus_hadCS"),
            ),
            axis=1,
        )

    try:
        df = df.join(readDatafile("MHall_uncertainties"))
    except:
        pass

    return df.drop(columns=[x for x in df.columns if "unused" in x]).fillna(0)
