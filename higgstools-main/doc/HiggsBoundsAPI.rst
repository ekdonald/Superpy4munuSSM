HiggsBounds
===========

The HiggsBounds library has a very simple interface that is defined in C++ by
including `Higgs/Bounds.hpp`. As for all of the API, unless otherwise noted the
Python and C++ functions and classes have the same names where C++ namespaces
correspond to python modules.

The main functionality of the HiggsBounds library is provided by the
:cpp:class:`Higgs::bounds::Bounds` class (also available as
:cpp:class:`Higgs::Bounds` for convenience). The usual use case given some model
`predictions` in a :cpp:class:`Higgs::predictions::Predictions` object is as
simple as

.. literalinclude:: examples/usingHB.cpp
    :language: c++

or equivalently in python

.. literalinclude:: examples/usingHB.py
    :language: python

.. doxygenclass:: Higgs::bounds::Bounds
    :members:


The :cpp:struct:`Higgs::bounds::HBResult` contains all of the results of a
HiggsBounds run. While the main result of interest is usually
:cpp:member:`Higgs::bounds::HBResult::allowed` all of the information on the
limits that this result is based on is included as well.

.. doxygenstruct:: Higgs::bounds::HBResult
    :members:

This allows you to study which limits actually constrain your model by e.g.
storing or printing information about the
:cpp:member:`Higgs::bounds::HBResult::selectedLimits`. For example, you could
store the Inspire-HEP cite keys for all of the limits that constrain your model,
since you probably want to cite those in your paper

.. code-block:: python

    citeKeys = set([x.limit().citeKey() for x in result.selectedLimits if x.obsRatio() > 1])

But you have access to much more information than just that, see the the
:ref:`Limits` section for details.



Limits
^^^^^^

The :cpp:class:`Higgs::bounds::Limit` class provides a unified interface to all
of the different kinds of limits implemented in HiggsBounds. The class offers a
lot of metadata on the limit as well as descriptions of the underlying process.

Those can be particularly useful to learn more about the
:cpp:member:`Higgs::bounds::HBResult::selectedLimits` or all of the loaded
:cpp:func:`Higgs::bounds::Bounds::limits`. For example, to figure out which CMS 13TeV
limits that involve di-photon final states are implemented you could use

.. code-block:: python

    set([l.reference() for l in bounds.limits() if l.experiment()=="CMS
                                                and l.collider()=="LHC13"
                                                and "gamgam" in l.processDesc()])

You can also load single limits from files using
:cpp:func:`Higgs::bounds::Limit::read` and apply them to
:cpp:class:`Higgs::predictions::Predictions` with
:cpp:func:`Higgs::bounds::Limit::apply` allowing you to explicitly study the
effect of individual limits.

.. doxygenclass:: Higgs::bounds::Limit
    :members:
    :membergroups: Constructor Evaluation Metadata


The :cpp:class:`Higgs::bounds::AppliedLimit` class contains a reference to the
underlying limit, and also stores the results of comparing it to model
evaluations.

.. doxygenclass:: Higgs::bounds::AppliedLimit
    :members: limit, obsRatio, expRatio, contributingParticles, obsLikelihood, expLikelihood


Limit Options
^^^^^^^^^^^^^

Finally, the :cpp:struct:`Higgs::bounds::LimitOptions` contains options that
govern how limits are applied to model predictions. It is recommended to keep
the default values for all of these, but the fine-grained control can sometimes
be useful.

.. doxygenstruct:: Higgs::bounds::LimitOptions
    :members:

In this context the :cpp:enum:`Higgs::predictions::MassUncEagerness` appears for the first time. It is defined as follows:

.. doxygenenum:: Higgs::predictions::MassUncEagerness
