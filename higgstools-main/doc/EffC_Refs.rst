Reference Models and the Effective Coupling Input
-------------------------------------------------

HiggsPredictions also provides access to a lot of tabulated cross section and
branching ratio predictions in different models and setups. These are grouped
into two parts, the implementations of reference models such as the SM Higgs,
and properties computed in the effective coupling approximation, either with
respect to one of the reference model or in a simplified model with a predefined
coupling structure. All of these are documented on this page.

Reference Models
^^^^^^^^^^^^^^^^

HiggsPredictions currently implements two reference models:
::cpp:class:`Higgs::predictions::SMHiggs` and
:cpp:class:`Higgs::predictions::SMHiggsEW` that both model the SM Higgs boson,
but with different higher order corrections included.

.. doxygenclass:: Higgs::predictions::SMHiggs
    :members:

.. doxygenclass:: Higgs::predictions::SMHiggsEW
    :members:

There is also a flag that can be used to select the reference model.

.. doxygenenum:: Higgs::predictions::ReferenceModel


Effective Couplings
^^^^^^^^^^^^^^^^^^^

The effective coupling approximation is a way to parametrize cross sections or
decay widths through their dependence on the couplings of some reference model.
Such parametrizations can then be used to obtain values in any model that is
related to the reference model by a simple rescaling of the individual
couplings. HiggsPredictions uses this approximation both as a convenient input
method for particles that are related to one of the predefined
:cpp:enum:`Higgs::predictions::ReferenceModel` and to provide tabulated cross
sections for some of the more exotic processes that are targeted in BSM searches
for which no dedicated calculation may be available in most models.


Effective Coupling Input for Neutral Higgs-Like Scalars
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
A wide class of BSM models feature one or more neutral scalars with a coupling
structure that is a simple rescaling of the SM Higgs couplings. For such
particle the effective coupling input can be used to easily set all of the
production and decay modes that are present for the SM Higgs by rescaling with
the appropriate effective couplings.

It is highly recommended to use this input method if it is applicable to get a
baseline set of cross sections and branching ratios that you can then expand
upon.

The effective couplings relative to the SM Higgs are defined as follows

.. doxygenstruct:: Higgs::predictions::NeutralEffectiveCouplings
    :members:

These couplings can then be used to set all possible cross sections, branching
ratios, and couplings through

.. doxygenfunction:: Higgs::predictions::effectiveCouplingInput


The SM Higgs boson can be represented by the following values of the effective
couplings. When using these couplings with
:cpp:func:`Higgs::predictions::effectiveCouplingInput` you are guaranteed to
get a particle that has identical cross sections, branching ratios and couplings
to the specified reference model.

.. doxygenvariable:: Higgs::predictions::smLikeEffCouplings

Finally, HiggsPredictions also offers a convenience function for the common case
where all SM couplings are rescaled by the same factor (e.g. in singlet
extensions).

.. doxygenfunction:: Higgs::predictions::scaledSMlikeEffCouplings



Parametrized Cross Sections
"""""""""""""""""""""""""""

HiggsPredictions also provides direct access to parametrizations of production
cross sections in terms of effective couplings. This includes many production
processes that are not present for a SM-like particle. As such they are not set
through the :cpp:func:`Higgs::predictions::effectiveCouplingInput` but must be
used manually to get a value for :cpp:func:`Higgs::predictions::BsmParticle::setCxn`.

In python the functions in this namespace are located in the sub-module
`Higgs.predictions.EffectiveCouplingCxns`.

.. doxygennamespace:: Higgs::predictions::EffectiveCouplingCxns


Parametrized Scaling Factors
""""""""""""""""""""""""""""

For some cross sections, for which high-precision reference model calculations
are available, we provide rescaling factors as a function of the effective
couplings instead. These can be used explicitly with
:cpp:func:`Higgs::predictions::BsmParticle::setNormalizedCxn` but are usually
mostly used internally through
:cpp:func:`Higgs::predictions::effectiveCouplingInput`.

.. doxygennamespace:: Higgs::predictions::EffectiveCouplingRatios
