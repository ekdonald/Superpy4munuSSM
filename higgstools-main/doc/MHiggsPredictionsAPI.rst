HiggsPredictions -- Mathematica Interface
=========================================

All Mathematica functions related to HiggsPredictions can be
recognized by the initial letters ``HP``. In contrast to the
C++ and python interfaces the ``Higgs::Predictions`` object is
automatically generated and is referenced internally if needed.

Add and remove particles
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Mathematica
    
    HPAddParticle[id_String, mass_Real, echarge_String, cp_String]

adds a particle with the name ``id``, the mass ``mass``, the :ref:`electrical charge <ECharge>` ``echarge``, and the :ref:`CP character <CP>` ``cp``.

.. code-block:: Mathematica
    
    HPRemoveParticle[id_String]

removes the particle ``id``.

.. code-block:: Mathematica
    
    HPGetParticleIDs[]

lists all particles.

.. code-block:: Mathematica
    
    HPGetCP[id_String]

returns the :ref:`CP character <CP>` of the particle ``id``.

.. code-block:: Mathematica
    
    HPGetCharge[id_String]

returns the :ref:`electrical charge <ECharge>` of the particle ``id``.

.. code-block:: Mathematica
    
    HPSetMass[id_String, value_Real]

sets the mass of the particle ``id`` to ``value``.

.. code-block:: Mathematica
    
    HPGetMass[id_String]

returns the mass of the particle ``id``.

.. code-block:: Mathematica
    
    HPSetMassUnc[id_String, value_Real]

sets the mass uncertainty of the particle ``id`` to ``value``.

.. code-block:: Mathematica
    
    HPGetMassUnc[id_String]

returns the mass uncertainty of the particle ``id``.

Set cross sections and branching ratios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Mathematica
    
    HPSetBrTopWb[value_Real]

sets :math:`\mathrm{BR}(t\to W^+ b)` to ``value`` (default: 1).

.. code-block:: Mathematica
    
    HPSetDecayWidth[id_String, decay_String, value_Real]
    
sets the decay width for the particle ``id`` for the :ref:`decay <Decays>` ``decay`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetDecayWidth[id1_String, chaindecay_String, iddaughter_String, value_Real]
    
sets the decay width for the particle ``id`` for the :ref:`chain decay <ChainDecays>` ``chaindecay`` and the daughter BSM particle ``iddaughter`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetDecayWidth[id_String, id1_String, id2_String, value_Real]
    
sets the decay width for the particle ``id`` to the daughter BSM particles ``id1`` and ``id2`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetBR[id_String, decay_String, value_Real]
    
sets the branching ratio for the particle ``id`` for the :ref:`decay <Decays>` ``decay`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetBR[id1_String, chaindecay_String, iddaughter_String, value_Real]
    
sets the branching ratio for the particle ``id`` for the :ref:`chain decay <ChainDecays>` ``chaindecay`` and the daughter BSM particle ``iddaughter`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetBR[id_String, id1_String, id2_String, value_Real]
    
sets the branching ratio for the particle ``id`` to the daughter BSM particles ``id1`` and ``id2`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetTotalWidth[id_String, value_Real]
    
sets the total width of the particle ``id`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetChannelRate[id_String, coll_String, prod_String, decay_String, value_Real]
    
sets the channel rate of the particle ``id`` at the :ref:`collider <Colliders and Experiments>` ``coll`` for the :ref:`production mode <Productions>` ``prod`` to the value ``value``.

.. code-block:: Mathematica
    
    HPSetBsmPairCxn[coll_String, id1_String, id2_String, value_Real]
    
sets the cross section for non-resonant pair production of BSM particles to the value ``value``, where ``coll`` identifies the :ref:`collider <Colliders and Experiments>`; ``id1``, the first produced BSM particle; and ``id2``, the second produced BSM particle.

.. code-block:: Mathematica
    
    HPSetCoupling[id_String, coup_String, value_Real]
    
sets the :ref:`coupling <Couplings>` ``coup`` of the particle ``id`` to ``value``.

.. code-block:: Mathematica
    
    HPGetCoupling[id_String, coup_String]
    
returns the :ref:`coupling <Couplings>` ``coup`` of the particle ``id``.

Effective coupling input
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Mathematica

    HPEffectiveCouplingInput[id_String]

sets the :ref:`effective couplings for <Effective Couplings>` for the particle ``id``.

The following optional arguments can be given:

* ``uuRe``: :math:`\kappa_u` (type: ``Real``; default: ``0``),
* ``uuIm``: :math:`\tilde\kappa_u` (type: ``Real``; default: ``0``),
* ``ddRe``: :math:`\kappa_d` (type: ``Real``; default: ``0``),
* ``ddIm``: :math:`\tilde\kappa_d` (type: ``Real``; default: ``0``),
* ``ccRe``: :math:`\kappa_c` (type: ``Real``; default: ``0``),
* ``ccIm``: :math:`\tilde\kappa_c` (type: ``Real``; default: ``0``),
* ``ssRe``: :math:`\kappa_s` (type: ``Real``; default: ``0``),
* ``ssIm``: :math:`\tilde\kappa_s` (type: ``Real``; default: ``0``),
* ``bbRe``: :math:`\kappa_b` (type: ``Real``; default: ``0``),
* ``bbIm``: :math:`\tilde\kappa_b` (type: ``Real``; default: ``0``),
* ``ttRe``: :math:`\kappa_t` (type: ``Real``; default: ``0``),
* ``ttIm``: :math:`\tilde\kappa_t` (type: ``Real``; default: ``0``),
* ``eeRe``: :math:`\kappa_e` (type: ``Real``; default: ``0``),
* ``eeIm``: :math:`\tilde\kappa_e` (type: ``Real``; default: ``0``),
* ``mumuRe``: :math:`\kappa_\mu` (type: ``Real``; default: ``0``),
* ``mumuIm``: :math:`\tilde\kappa_\mu` (type: ``Real``; default: ``0``),
* ``tautauRe``: :math:`\kappa_\tau` (type: ``Real``; default: ``0``),
* ``tautauIm``: :math:`\tilde\kappa_\tau` (type: ``Real``; default: ``0``),
* ``WW``: :math:`\kappa_W` (type: ``Real``; default: ``0``),
* ``ZZ``: :math:`\tilde\kappa_{ZZ}` (type: ``Real``; default: ``0``),
* ``Zgam``: :math:`\tilde\kappa_{Z\gamma}` (type: ``Real``; default: ``0``),
* ``gamgam``: :math:`\tilde\kappa_\gamma` (type: ``Real``; default: ``0``),
* ``gg``: :math:`\tilde\kappa_{gg}` (type: ``Real``; default: ``0``),
* ``refModel``: the  :ref:`reference model <Reference Models>` (type: ``String``; default: ``SMHiggsEW``),
* ``ggH``: whether to calculate the ggH cross-section in terms of the effective top and bottom Yukawa couplings or by rescaling the SM-like ggH XS by the squared of the effective gg coupling (type: ``Boolean``; default: ``True``).

The options can be entered e.g. via ``HPeffectiveCouplingInput[h, ttRe -> 1, ttIm -> 0.5]``.

.. code-block:: Mathematica

    HPScaledSMlikeEffCouplings[id_String, scale_Real]

sets the :ref:`effective couplings for <Effective Couplings>` for the particle ``id`` by rescaling the SM couplings by ``scale``.

The following optional arguments can be given:

* ``refModel``: the  :ref:`reference model <Reference Models>` (type: ``String``; default: ``SMHiggsEW``).

.. code-block:: Mathematica

    HPSMLikeEffCouplings[id_String]

sets the :ref:`effective couplings for <Effective Couplings>` for the particle ``id`` to the SM values.

The following optional arguments can be given:

* ``refModel``: the  :ref:`reference model <Reference Models>` (type: ``String``; default: ``SMHiggsEW``).

Retrieve cross sections and branching ratios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: Mathematica
    
    HPGetSMTotalWidth[mass_Real]

returns the total width for a SM Higgs with the mass ``mass``.

.. code-block:: Mathematica
    
    HPGetSMBr[mass_Real, decay_String]

returns the branching ratio for a SM Higgs with the mass ``mass`` for the :ref:`decay <Decays>` ``decay``.

.. code-block:: Mathematica
    
    HPGetBrTopWb[]

returns :math:`\mathrm{BR}(t\to W^+ b)`.

.. code-block:: Mathematica
    
    HPGetSMCxn[mass_Real, collstring_String, prodstring_String]

returns the cross section for a SM Higgs with the mass ``mass`` for the :ref:`collider <Colliders and Experiments>` ``coll`` and the :ref:`production mode <Productions>` ``prod``.

.. code-block:: Mathematica
    
    HPGetCxn[id_String, coll_String, prod_String]

returns the cross section for the particle ``id`` at the :ref:`collider <Colliders and Experiments>` ``coll`` for the :ref:`production mode <Productions>` ``prod``.

.. code-block:: Mathematica
    
    HPGetBR[id_String, decay_String]
    
returns the branching ratio for the particle ``id`` for the :ref:`decay <Decays>` ``decay`` ``decay``.

.. code-block:: Mathematica
    
    HPGetBR[id1_String, chaindecay_String, iddaughter_String, value_Real]
    
returns the branching ratio for the particle ``id`` for the :ref:`chain decay <ChainDecays>` ``chaindecay`` and the daughter BSM particle ``iddaughter``.

.. code-block:: Mathematica
    
    HPGetBR[id_String, id1_String, id2_String, value_Real]
    
returns the branching ratio for the particle ``id`` to the daughter BSM particles ``id1`` and ``id2``.

.. code-block:: Mathematica
    
    HPGetTotalWidth[id_String]
    
returns the total width of the particle ``id``.

.. code-block:: Mathematica
    
    HPGetChannelRate[id_String, coll_String, prod_String, decay_String, value_Real]
    
returns the channel rate of the particle ``id`` at the :ref:`collider <Colliders and Experiments>` ``coll`` for the :ref:`production mode <Productions>` ``prod``.

.. code-block:: Mathematica
    
    HPGetBsmPairCxn[coll_String, id1_String, id2_String]
    
returns the cross section for non-resonant pair production of BSM particles, where ``coll`` identifies the :ref:`collider <Colliders and Experiments>`; ``id1``, the first produced BSM particle; and ``id2``, the second produced BSM particle.