Particles and their Properties
------------------------------

The :cpp:class:`Higgs::predictions::Particle` class provides an interface to all
kinds of particles in HiggsPredictions.

.. doxygenclass:: Higgs::predictions::Particle
    :members:


Model Predictions for BSM Particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :cpp:class:`Higgs::predictions::BsmParticle` is the main class that handles
model input. All model-predicted production and decay rates of a particle can be
specified here. Note that the properties of the SM-like Higgs boson in any model
also has to be implemented through this class, even though this may not
technically be a BSM particle.

.. doxygenclass:: Higgs::predictions::BsmParticle
    :members:


Production Mode, Decay Modes, and Couplings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

HiggsPredictions defines a large set of different production and decay channels
that can be set by the user. In Python all of these values are can be
alternatively represented through a simple string matching the enumerator name
(as used in the previous python examples).

.. _Productions:

.. doxygenenum:: Higgs::predictions::Production

.. _Decays:

.. doxygenenum:: Higgs::predictions::Decay

.. _ChainDecays:

.. doxygenenum:: Higgs::predictions::ChainDecay

.. _Couplings:

.. doxygenenum:: Higgs::predictions::Coupling

Along with all of these definitions come a few functions that can be used to
check if a coupling is valid:

.. doxygenfunction:: Higgs::predictions::validProductionFor

.. doxygenfunction:: Higgs::predictions::validDecayFor

.. doxygenfunction:: Higgs::predictions::validProductionAt

