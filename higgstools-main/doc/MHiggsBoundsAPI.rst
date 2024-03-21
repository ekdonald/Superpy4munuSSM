HiggsBounds -- Mathematica Interface
====================================

All Mathematica functions related to HiggsBounds can be
recognized by the initial letters ``HB``. In contrast to the
C++ and python interfaces the ``Higgs::Bounds`` object is
automatically generated and is referenced internally if needed.

.. code-block:: Mathematica
    
    HBInitialize[path_String]

initializes HiggsBounds using the limit data files at the given path. 

The following optional arguments can be given:

* ``applicableResolutionFac`` (type: ``Real``; default: ``0.5``)
* ``clusterMassUnc`` (type: ``String``; default: ``cautious``)
* ``applicableMassunc`` (type: ``String``; default: ``cautious``)
* ``setLimitMassUnc`` (type: ``String``; default: ``cautious``)

See :ref:`Limit Options <Limit Options>` for an explanation of the various options.

.. code-block:: Mathematica
    
    HBRetrieveOptions[]

retrieves options set with ``HBInitialize``.

.. code-block:: Mathematica
    
    HBListLimits[]

returns a list of all implemented experimental limits.

.. code-block:: Mathematica
    
    HBApplyBounds[]

returns whether the parameter point is allowed at 95% C.L..

.. code-block:: Mathematica
    
    HBGetAppliedBounds[]

returns a list of all applied experimental limits.

.. code-block:: Mathematica
    
    HBGetSelectedBounds[]

returns a list of all selected limits.


