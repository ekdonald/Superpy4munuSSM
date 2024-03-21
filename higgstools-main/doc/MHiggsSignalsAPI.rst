HiggsSignals -- Mathematica Interface
=====================================

All Mathematica functions related to HiggsSignals can be
recognized by the initial letters ``HS``. In contrast to the
C++ and python interfaces the ``Higgs::Signals`` object is
automatically generated and is referenced internally if needed.

.. code-block:: Mathematica
    
    HSInitialize[path_String]

initializes HiggsSignals using the measurement data at the given path. 

The optional arguments

* ``theoryMassUncPdf`` (type: ``Real``; default: 0.5)
* ``massSensitiveAssignmentRange`` (type: ``String``; default: ``cautious``)
* ``unassignedMassMeasurementPenalty`` (type: ``String``; default: ``cautious``)
* ``unassignedCouplingMeasurementPenalty`` (type: ``String``; default: ``cautious``)
* ``rescaleToRefMass`` (type: ``String``; default: ``cautious``)

allow to control the behavior of the :math:`\chi^2` calculation. See :ref:`Measurement Options <Measurement Options>` for an explanation of the various options.

.. code-block:: Mathematica
    
    HSGetChisq[]

returns the :math:`\chi^2` value for the defined model predictions.

.. code-block:: Mathematica
    
    HSGetChisqMeasurement[id_Integer]

returns the :math:`\chi^2` value for the measurement with the id ``id``.

.. code-block:: Mathematica
    
    HSGetObservableCount[]

returns the number of observables used in the :math:`\chi^2` calculation.

.. code-block:: Mathematica
    
    HSListMeasurements[]

lists all loaded measurements.

.. code-block:: Mathematica
    
    HSListSubMeasurements[id_Integer]

lists all submeasurements for the measurement with the id ``id``.


