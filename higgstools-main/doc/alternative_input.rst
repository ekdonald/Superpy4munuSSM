Input via SLHA or data files
============================

As an alternative to providing the input via the HiggsPredictions routine, SLHA or HiggsBound data files (see e.g. Sec. 5.2.1 of `2006.06007 <https://arxiv.org/pdf/2006.06007.pdf>`_) can be used as an input.

The necessary routines are part of the `Input` python subpackage, which can be loaded e.g. via

.. code-block:: python

    import Higgs.tools.Input as hinput

The main function of this package parse predictions contained in a dictionary into a `Higgs.predictions.Predictions` object. The necessary dictionary can be extracted from a SLHA file or a HB5 data file.

.. py:function:: Higgs.tools.Input.predictionsFromDict(data: dict, neutralIds: list, singleChargedIds: list, doubleChargedIds: list, referenceModels: dict = {}, massPrefix="m", massUncPrefix="dm", widthPrefix="w", CpPrefix="CP", brPrefix="br", cxnPrefix="cxn", pairCxnPrefix="paircxn", effCPrefix="effc", oddCoupSuffix="p", normCxnSuffix="norm", useExplicitBr: Union[bool, dict] = True)

    Construct a Higgs.predictions.Predictions object from a set of model
    predictions given in the form of a dictionary of cxns, BRs and couplings.

    If effective couplings of the neutral scalars are provided, the effective
    coupling input is used first, and then any explicitely specified cross
    sections and branching ratios are used to overwrite the values obtained from
    the effective coupling approximation.

    :param data: the dictionary containing the data. For a particle with ID `"ID"`
        the following quantities can be specified: 

        * `m_ID`: **mass of the particle**
        * `dm_ID`: theoretical mass uncertainty 
        * `w_ID`: **total width of the particle** 
        * `br_ID_XX`: BR(ID->XX) for `"XX" in Higgs.predictions.Decay`
        * `br_ID_ID2_ID3`: BR(ID-> ID2 ID3) where `ID2` and `ID3` are either IDs of other particles, or members of `Higgs.predictions.ChainDecay`
        * `cxn_ID_XX_COLL`: a production cross section in the `Higgs.predictions.Production` mode `"XX"` at the `Higgs.predictions.Collider` `"COLL"`
        * `cxn_id_XX_COLL_norm`: a normalized production cross section in the `Higgs.predictions.Production` mode `"XX"` at the `Higgs.predictions.Collider` `"COLL"`

        Additional recognized fields for neutral particles only are:
        
        * `CP_ID`: **CP quantum number (only for neutral particles)**
        * `effc_ID_ff_s`: scalar part of the SM-normalized effective coupling to the fermions `ff`
        * `effc_ID_ff_p`: pseudoscalar part of the SM-normalized effective coupling to the fermions `ff`
        * `effC_ID_VV`: SM-normalized effective coupling to the gauge bosons `VV`
        
        The **highlighted entries** are required. Finally, non-resonant pair
        production cross sections of the particles `ID1` and `ID2` at a
        `Higgs.predictions.Collider` `COLL` can be specified as
        `paircxn_ID1_ID2_COLL`.

    :type data: dict

    :param neutralIds: IDs of the neutral particles
    
    :type neutralIds: list

    :param singleChargedIds: IDs of the singly charged particles
    
    :type singleChargedIds: list

    :param doublyChargedIds: IDs of the doubly charged particles
    
    :type doublyChargedIds: list

    :param referenceModels: a dictionary that can be used to assign a
        `Higgs.predictions.ReferenceModel` to any particle. If none is specified
        this defaults to `Higgs.predictions.SMHiggsInterp`.
    
    :type referenceModels: dict

    :param useExplicitBr: whether the explicitly given BRs should be used. 
        Alternatively, the BRs are calculated based on the given effective 
        couplings. This can also be controlled for individual Higgs bosons and/or
        decay channels by providing a dictionary (e.g.
        `useExplicitBr = {'h1': False, 'h2': {'WW': False}}`)
    
    :type useExplicitBr: Union[bool, dict]
    
    :return: HiggsPredictions object containing all the predictions extracted from the input data.
    :rtype: `Higgs.predictions.Predictions`

The remaining optional arguments can be used to adjust the prefixes and
suffixes that make up the keys as described above.

If no large effects of BSM particles on the SM-like Higgs branching ratios are
expected, we recommend using the effective coupling input for the SM-like Higgs.
Most tools for the calculation of the branching ratios lack the precision level
required to match the experimental sensitivity resulting in chi^2 penalties 
even if the Higgs couplings are completely SM-like. 

If the effective coupling input is used for (some of) the decay channels, 
the SM decay widths (and not the branching ratios) are rescaled. As a result,
explicitly given branching ratios (e.g. for a BSM decay of the SM-like Higgs boson)
might have a slightly different value after the calculation of the other decay
channels via the effective couplings. Due to the same reason, also the total width
might be shifted slightly with respect to the provided width.



SLHA input
^^^^^^^^^^

The SLHA input relies on the function

.. py:function:: Higgs.tools.Input.readHB5SLHA(file: str, neutralPDGs: list, chargedPDGs: list, invisibleWidthThreshold: float = 1e-10, invisiblePDGs: list = [])

    Reads a HiggsBounds-5 compatible SLHA file into a dictionary. The SLHA file
    must contain the blocks HiggsCouplingsBosons and HiggsCouplingsFermions that
    contain the effective couplings of the neutral Higgs bosons as defined in
    the HiggsBounds-5 manual. Additionally, the masses of all considered
    particles have to be present in the MASS block, and a DECAY block has to
    exist for every considered particle and for the top-quark.

    **Uses and requires the pylha package.**

    This function is slightly more general than the HB5 SLHA interface, in that
    it supports models other than the MSSM and NMSSM, by using the PDG numbers
    of the neutral and charged particles as input.

    :param file: filename of the SLHA file to read

    :type file: str

    :param neutralIds: IDs of the neutral particles
    
    :type neutralIds: list

    :param chargedPDGs: IDs of the charged particles
    
    :type chargedPDGs: list

    :param invisibleWidthThreshold: any particles with total widths less than this
        quantity (in GeV) are treated as invisible (does not apply for the 
        SM quarks and charged leptons; only used if `inivisible PDGs = []`)
        
    :type invisibleWidthThreshold: float
    
    :param invisiblePDGs: list of the PDG numbers for invisible particles. If `invisiblePDGs = []`, which is the default, `invisibleWidthThreshold` is used.

    :type invisiblePDGs: list
    
    :return: A dictionary containing all of the data relevant for HiggsTools that could be extracted from the SLHA file. See `predictionsFromData`` for a detailed description of the dictionary keys.
    :rtype: dict


Data file input
^^^^^^^^^^^^^^^

The HB5 data file input relies on the function

.. py:function:: Higgs.tools.Input.readHB5Datafiles(prefix: str, neutralIds: list, chargedIds: list)

    Reads input given in the form of the HiggsBounds-5 datafiles into a pandas
    dataframe. Each row of the dataframe can then be used with
    predictionsFromDict as `df.apply(predictionsFromDict, axis=1)`.

    If both effective coupling and hadronic input data is available for the
    given prefix, the effective coupling input is used first, and then any
    explicitely provided cxns and BRs are overwritten.

    :param prefix: the prefix path to the datafiles, as defined in the HiggsBounds-4/5 manuals

    :type prefix: str

    :param neutralIds: IDs of the neutral particles, given in the order they appear in the datafiles
    
    :type neutralIds: list

    :param chargedPDGs: IDs of the charged particles, given in the order they appear in the datafiles
    
    :type chargedPDGs: list
    
    :return: A pandas DataFrame where each row corresponds to a datapoint (row) given in the datafiles. See predictionsFromDict for a detailed description of the column names.
    :rtype: `pandas.dataframe`