The Mathematica interface
=============================

This documentation provides detailed information on the Mathematica user interfaces of
HiggsTools. As a consequence of the structure of the Wolfram language (i.e., no classes),
the Mathematica interfaces differs from the C++ and python interfaces.

Additionally, functions in the python interface taking C++ enumeration types as
arguments, support automatic conversion from strings to enumeration values. The
Mathematica executable is build by specifying 

.. code-block:: bash

    cmake -DHiggsTools_BUILD_MATHEMATICA_INTERFACE=ON ..

during the build process. The resulting executable is located in ``build/wstp/MHiggsTools``
and can be loaded via

.. code-block:: Mathematica

    Install["/path/to/MHiggsTools"];

from within Mathematica. 

At the moment, the Mathematica interface lacks a few minor functionalities of the
C++ and python interfaces. If you are missing a specific feature, please contact
the developers.

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    MHiggsPredictionsAPI
    MHiggsBoundsAPI
    MHiggsSignalsAPI
