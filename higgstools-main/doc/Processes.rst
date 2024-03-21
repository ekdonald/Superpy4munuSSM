Processes
---------

The following introduces the different process topologies that are defined in
HiggsPredictions. These are not intended to be used by the user **and do not
even have a corresponding Python interface**. However, the definition of the
different cases and the detailed discussion of symmetry factors may still be of
interest.

The simplest process in HiggsPredictions is the
:cpp:class:`Higgs::predictions::ChannelProcess`. It is used to model any process
with a single particle of interest produced in a SM initial state and decaying
into a SM final state. This covers all single-particle resonance searches in
HiggsBounds and all existing rate measurements in HiggsSignals. Even for mass
and coupling measurements in HiggsSignals a
:cpp:class:`Higgs::predictions::ChanelProcess` is always used for assignment and
weighting of individual particles.

Due to its simplicity and many use cases this process class supports some
additional features compared to the others, such as merging processes together
and splitting them into sub-process.

.. doxygenclass:: Higgs::predictions::ChannelProcess
    :members:


All of the remaining processes involve more than one
:cpp:class:`Higgs::predictions::BsmParticel`. As such they are currently (sadly)
not relevant for HiggsSignals, but are used extensively in HiggsBounds.

The simplest process involving more than one particle is the
:cpp:class:`Higgs::predictions::ChainDecayProcess` involving two particles of
interest, one decaying into the other through a
:cpp:enum:`Higgs::predictions::ChainDecay`.

.. doxygenclass:: Higgs::predictions::ChainDecayProcess
    :members:

The simplest topology involving three particles of interest, one mother particle
decaying into two daughter particles, is modeled by the
:cpp:class:`Higgs::predictions::PairDecayProcess`.

.. doxygenclass:: Higgs::predictions::PairDecayProcess
    :members:

Finally, non-resonant pair production of BSM particles is covered through the
:cpp:class:`Higgs::predictions::PairProductionProcess`.

.. doxygenclass:: Higgs::predictions::PairProductionProcess
    :members:

