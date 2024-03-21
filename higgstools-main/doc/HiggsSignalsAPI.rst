HiggsSignals
============

Much like HiggsBounds, the HiggsSignals library has an extremely simple
interface that is defined in C++ by including `Higgs/Signals.hpp`. As for all of
the API, unless otherwise noted the Python and C++ functions and classes are in
1-to-1 correspondence with C++ namespaces mapped to Python modules.

The main functionality of HiggsSignals is provided by the
:cpp:class:`Higgs::signals::Signals` class. The usual use case given some model
predictions in a :cpp:class:`Higgs::predictions::Predictions` object is as
simple as

.. literalinclude:: examples/usingHS.cpp
    :language: c++

or equivalently in python

.. literalinclude:: examples/usingHS.py
    :language: python

.. doxygenclass:: Higgs::signals::Signals
    :members:


Measurements
^^^^^^^^^^^^

All experimental results in HiggsSignals are implemented through the
:cpp:class:`Higgs::signals::Measurement` class. Each instance of the class
typically implements the results of one experimental analysis and directly
corresponds to a single `json` implementation file in the HiggsSignals dataset.

.. doxygenclass:: Higgs::signals::Measurement
    :members:

Each measurement contains one or more
:cpp:class:`Higgs::predictions::SubMeasurement` that represent a single measured
observable. That could be a measured rate, like a signal strength in some
production mode or a single STXS bin, but also a mass measurement or a coupling
measurement. The :cpp:class:`Higgs::predictions::SubMeasurement` represents the
common interface for all of these observables. As a user, you should rarely need
to interact with a :cpp:class:`Higgs::signals::SubMeasurement` directly, but it
is still documented here for completeness.

.. doxygenclass:: Higgs::signals::SubMeasurement
    :members:

Modification Factors
^^^^^^^^^^^^^^^^^^^^

HiggsSignals supports the notion of *modification factors*. A modification
factor is a multiplicative scaling applied to the model prediction for a single
channel (i.e. the rate in a :cpp:enum:`Higgs::predictions::Production` and
:cpp:enum:`Higgs::predictions::Decay` mode at a given
:cpp:enum:`Higgs::predictions::Collider`) in the context of a specific
:cpp:class:`Higgs::signals::SubMeasurement`.

The purpose of these modification factors is to allow capturing differential BSM
effects that go beyond the inclusive rates stored in the
:cpp:class:`Higgs::predictions::Predictions`.

For example, a model could predict that (for whatever reason) radiative jets in
the :math:`gg -> H -> \gamma\gamma` process are harder than they would be in the SM.
Since there are STXS measurements of this process in HiggsSignals this
prediction can be tested by adjusting the modification factors for the relevant
STXS bins. A Python example of this (for just one analysis and only the 1-jet
STXS bins) would look like:

.. literalinclude:: examples/modificationFactors.py
    :language: python

The whole dictionary corresponds to the
:cpp:type:`Higgs::signals::Signals::ModificationFactors`. As noted in its
documentation, the keys at this first level specify the
:cpp:func:`Higgs::signals::Measurement::id`, in this example modification
factors are set for the `ATLAS 2103.06956 analysis
<https://arxiv.org/abs/2103.06956>`_. When evaluating an analysis for which
modification factors are specified, the corresponding dictionary value is passed
on as :cpp:type:`Higgs::signals::Measurement::ModificationFactors`. In this
second level dictionary each key names a
:cpp:type:`Higgs::signals::SubMeasurement` contained in the
:cpp:type:`Higgs::signals::Measurement`. The keys can be obtained from
:cpp:func:`Higgs::signals::Measurement::subMeasurements` or by looking into the
implementation file for the analysis. Again, the array of values for each
sub-measurement name forms the corresponding
:cpp:type:`Higgs::signals::SubMeasurement::ModificationFactors` which are passed
to the corresponding :cpp:class:`Higgs::signals::SubMeasurement` when evaluating
it. Each value in this array is a rescaling factor applied to the respective
channel when evaluating the model predicted rate in the signal
:cpp:class:`Higgs::predictions::ChannelProcess` of the sub-measurement. In the
example there is only one entry in each array since there is only one channel
that contributes to the modified bins. If there are multiple channels, their
order (which the modification factors has to match) can be retrieved through
either the :cpp:func:`Higgs::signals::SubMeasurement::processDesc` (see there)
or by looking into the implementation file.

It is perfectly fine to only specify a subset of the modification factors, like
shown in the example. Any modification factors that are not set will default to
1, i.e. no scaling is done.


Measurement Options
^^^^^^^^^^^^^^^^^^^

The :cpp:struct:`Higgs::signals::MeasurementOptions` contains options that
govern different aspects of the Measurement and Sub-Measurement evaluation and
:math:`\chi^2` computation. The provided default values for all models are
reasonable general choices, but there are certainly models and applications
where different choices could be argued for.

.. doxygenstruct:: Higgs::signals::MeasurementOptions
    :members:

These options rely on the following definitions

.. doxygenenum:: Higgs::signals::PDF

.. doxygenenum:: Higgs::signals::RescaleToRefMass

.. doxygenenum:: Higgs::signals::Correlations



HiggsSignals Errors
^^^^^^^^^^^^^^^^^^^

Finally, if a Measurement file cannot be read the following error is thrown.

.. doxygenclass:: Higgs::signals::InvalidMeasurement
