The C++ and Python Interfaces
=============================

This documentation provides detailed information on the main user interfaces of
HiggsTools. The documentation is generated from the C++ code, however the
structure of the Python interface is identical, with C++ namespaces
corresponding to Python modules, e.g. the
:cpp:class:`Higgs::predictions::BsmParticle` C++ class is mapped to the
:class:`Higgs.predictions.BsmParticle` Python class.

Additionally, functions in the python interface taking C++ enumeration types as
arguments, support automatic conversion from strings to enumeration values. For
example in

.. code-block:: python

    import Higgs.predictions as HP

    h1 = HP.BsmParticle("h1", "neutral", "even"))
    h2 = HP.BsmParticle("h2", HP.ECharge.neutral, HP.CP.even))

    h1.setCxn("LHC13","ggH",10.)
    h2.setCxn(HP.Collider.LHC13, HP.Production.ggH, 10.)

both `h1` and `h2` are CP-even, neutral scalars with equal gluon fusion cross
sections at the 13TeV LHC.


.. toctree::
    :maxdepth: 2
    :caption: Contents:

    HiggsPredictionsAPI
    HiggsBoundsAPI
    HiggsSignalsAPI
