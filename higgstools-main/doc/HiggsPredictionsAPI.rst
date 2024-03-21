HiggsPredictions
=================


.. toctree::
    :maxdepth: 2
    :caption: Contents:

    Particles
    EffC_Refs
    Basics
    Processes


HiggsPredictions handles model predictions for BSM particles. It both offers a
framework to input externally computed predictions, and contains tabulated
properties for several reference particles that can be used to easily obtain
model predictions using the effective coupling approximation. The whole
interface can be included through the `#include "Higgs/Predictions.hpp"` in C++
and is contained in the `Higgs.predictions` sub-package in python.

All of the model predictions are stored in one object of the
:cpp:class:`Higgs::predictions::Predictions` class (also available through the
convenience typedef :cpp:class:`Higgs::Predictions`).

A typical usage in C++ looks as follows

.. literalinclude:: examples/addingParticles.cpp
    :language: c++

or equivalently in python

.. literalinclude:: examples/addingParticles.py
    :language: python

On top of storing the individual particles,
:cpp:class:`Higgs::predictions::Predictions` also stores and provides access to
properties that do not depend on just one BSM particle, such as non-resonant
pair-production cross sections.

.. doxygenclass:: Higgs::predictions::Predictions
    :members:
