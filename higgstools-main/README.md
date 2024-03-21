# HiggsTools

[![pipeline status](https://gitlab.com/higgsbounds/higgstools/badges/develop/pipeline.svg)](https://gitlab.com/higgsbounds/higgstools/-/commits/develop) 
[![coverage report](https://gitlab.com/higgsbounds/higgstools/badges/develop/coverage.svg)](https://gitlab.com/higgsbounds/higgstools/-/commits/develop)

HiggsTools is a complete rewrite and unification of [HiggsBounds-5] and
[HiggsSignals-2] in modern C++. The library includes a fully functional python
interface that provides easy access to all features of the code. Compared to
[HiggsBounds-5] and [HiggsSignals-2] HiggsTools not only features substantial
technical and interface improvements, but also provides a much more robust and
more consistent representation of both model predictions and experimental results.

[**Detailed documentation of the C++, Python, and Mathematica user interfaces is available here**][documentation]

If you encounter any problems please [open an issue][issues] and if you have questions feel free to [write the team][contact].

### The HiggsBounds Collaboration
HiggsTools is developed by the HiggsBounds collaboration.

Current members of the team are Henning Bahl, Thomas Biekötter, Philip Bechtle, Sven Heinemeyer, Cheng Li, Steven Paasch, Georg Weiglein, and Jonas Wittbrodt.

Former members are Oliver Brein, Daniel Dercks, Tobias Klingl, Oscar Stål, Tim Stefaniak, and Karina E. Williams.





## What is HiggsTools?

HiggsTools is a toolbox for comparing a wide class of BSM models to all available experimental results from searches for new (scalar) particles and measurements of the 125GeV Higgs boson at colliders.

HiggsTools is composed of three sub-libraries:

### HiggsPredictions

HiggsPredictions handles the model predictions and user input. This was
previously part of [HiggsBounds-5], but has been separated out for clarity.
HiggsPredictions provides an object oriented interface in both Python and C++
that can be used to set rates for all production and decay modes among the
particles of the model. 

HiggsPredictions also contains a large repository of tabulated cross sections
and branching ratios in reference models, such as the SM Higgs boson, or
parametrized in the effective coupling approximation. This easily accessible
information can allow easily testing many properties of BSM models for which no
dedicated cross section or decay calculators are available.

### HiggsBounds

HiggsBounds takes the model predictions from HiggsPredictions and tests
them against a huge database of experimental results from searches for new
particles. Out of all the implemented limits, HiggsBounds selects the most
sensitive limit to every particle in the model based on the ratio between the
model prediction and the expected limit. The model predictions for all of these
selected limits are then required to lie below the observed limit to obtain an
approximate overall 95% CL exclusion bound. 

The HiggsBounds dataset is contained in a [separate repository][HBDataset].
Detailed information on the implemented experimental results can be found there.

### HiggsSignals

HiggsSignals takes the model predictions from HiggsPredictions and compares them
to a full set of current measurements of the 125GeV Higgs properties. Using all
available correlation information, it computes and returns a $\chi^2$ value that
quantifies the agreement of the model predictions with the measurements.

This HiggsSignals dataset is also contained in a [separate
repository][HSDataSet] with detailed information on the implementation available
there.

## Journal references

  - Henning Bahl, Thomas Biekötter, Sven Heinemeyer, Cheng Li, Steven Paasch Georg Weiglein, 
    Jonas Wittbrodt                                                                   <br/>
    *HiggsTools: BSM scalar phenomenology with new versions of HiggsBounds and HiggsSignals*                                                                     <br/>
    e-Print: [arXiv:2210.09332] [hep-ph]
#### HiggsBounds

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Georg Weiglein, 
    Karina E. Williams                                                                <br/>
    *HiggsBounds: Confronting Arbitrary Higgs Sectors with Exclusion Bounds from
    LEP and the Tevatron*                                                             <br/>
    Comput.Phys.Commun.181:138-167, 2010, e-Print: [arXiv:0811.4169] [hep-ph]

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Georg Weiglein, 
    Karina E. Williams                                                                <br/>
    *HiggsBounds 2.0.0: Confronting Neutral and Charged Higgs Sector Predictions 
    with Exclusion Bounds from LEP and the Tevatron*                                  <br/>
    Comput.Phys.Commun. 182:2605-2631, 2011, e-Print: [arXiv:1102.1898] [hep-ph]

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, 
    Georg Weiglein, Karina Williams                                                   <br/>
    *Recent Developments in HiggsBounds and a Preview of HiggsSignals*                <br/>
    PoS CHARGED2012 (2012) 024, e-Print: [arXiv:1301.2345] [hep-ph]

  - Philip Bechtle, Oliver Brein, Sven Heinemeyer, Oscar Stål, Tim Stefaniak,
    Georg Weiglein, Karina Williams                                                   <br/>
    *HiggsBounds-4: Improved Tests of Extended Higgs Sectors against Exclusion
    Bounds from LEP, the Tevatron and the LHC*                                        <br/>
    Eur.Phys.J C74 (2014) 2693, e-Print: [arXiv:1311.0055] [hep-ph]

  - Philip Bechtle, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, Georg Weiglein        <br/>
    *Applying Exclusion Likelihoods from LHC Searches to Extended Higgs Sectors*      <br/>
    Eur.Phys.J. C75 (2015) no.9, 421, e-Print: [arXiv:1507.06706] [hep-ph]

  - Philip Bechtle, Daniel Dercks, Sven Heinemeyer, Tobias Klingl,  Tim Stefaniak,
    Georg Weiglein, Jonas Wittbrodt                                                   <br/>
    *HiggsBounds-5: Testing Higgs Sectors in the LHC 13 TeV Era*                      <br/>
    Eur.Phys.J.C 80 (2020) 12, 1211, e-Print: [arxiv:2006.06007] [hep-ph]

  - Henning Bahl, Victor Martin Lozano, Tim Stefaniak, Jonas Wittbrodt              <br/>
    *Testing Exotic Scalars with HiggsBounds*                                         <br/>
    Eur.Phys.J.C 82 (2022) 7, 584 e-Print: [arxiv:2109.10366] [hep-ph]

#### HiggsSignals

  - Philip Bechtle, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, Georg Weiglein    <br/>
    "HiggsSignals: Confronting arbitrary Higgs sectors with measurements at the
    Tevatron and the LHC"                                                         <br/>
    Eur.Phys.J. C74 no.2, 2711, 2014, e-Print: [arxiv:1305.1933] [hep-ph]

  - Oscar Stål, Tim Stefaniak  <br/>
    "Constraining extended Higgs sectors with HiggsSignals"                       <br/>
    PoS EPS-HEP2013 314, 2013, e-Print: [arxiv:1310.4039] [hep-ph]

  - Philip Bechtle, Sven Heinemeyer, Oscar Stål, Tim Stefaniak, Georg Weiglein    <br/>
    "Probing the Standard Model with Higgs signal rates from the Tevatron, 
    the LHC and a future ILC"                                                     <br/>
    JHEP 1411 039, 2014, e-Print: [arxiv:1403.1582] [hep-ph]

  - Philip Bechtle, Sven Heinemeyer, Tobias Klingl, Tim Stefaniak, 
    Georg Weiglein, Jonas Wittbrodt                                                <br/>
    "HiggsSignals-2: Probing new physics with precision Higgs measurements 
    in the LHC 13 TeV era"                                                         <br/>
    e-Print: [arxiv:2012.09197] [hep-ph]
## Usage Guide

Compiling HiggsTools requires only a c++17 compliant compiler:
 - gcc >= 9 (gcc-8 does not work)
 - clang >= 5

as well as 
 - CMake >= 3.17 (a working CMake version can be easily installed through
   [pip][pip cmake], if your OS does not offer a sufficiently new package)

Building the python interface also requires python >= 3.5 and the corresponding
development headers (e.g. the `python3-dev` package on ubuntu). All other
dependencies are compile-time only, and are automatically downloaded by CMake.


### Building the C++ Library

The C++ library version of HiggsTools can be built by running e.g.:
```bash
mkdir build && cd build
cmake ..
make
```

The C++ library and the associated header files can also be installed for use
with other codes by running 
```bash
make install
```
You can easily change the install location from the default using the [DESTDIR
variable](https://cmake.org/cmake/help/latest/envvar/DESTDIR.html#envvar:DESTDIR).

Note that installing the C++ library will not install the Python module.


### Building the Python Module

Compilation of the python module is even easier. This can be done by simply
running
```bash
pip install .
```
from within the root folder of the repository. Note that this build applies a
lot of optimizations and may take a few minutes, but does not show any progress
indicators by default. Installing the python interface this way will not install
the C++ interface.

## Building the Mathematica Interface

HiggsTools also provides a functional Mathematica interface. To build this use
```bash
mkdir build && cd build
cmake -DHiggsTools_BUILD_MATHEMATICA_INTERFACE=ON ..
make
```
when compiling the C++ library. This will result in a Mathematica wstp
executable being generated at `build/wstp/MHiggsTools`.

Depending on the specific Mathematica installation, cmake might not be able
to find all required Mathematica files in order to build the interface.
In this case, the user can specify manually the paths to the Mathematica
installation. Depending on the system and the Mathematica version,
different cmake flags have to be set. As an example, on a Linux computer
with an installation of Mathematica version 13.0.1, a complete cmake call
defining all relevant paths manually is
```bash
cmake \
    -DHiggsTools_BUILD_MATHEMATICA_INTERFACE=ON \
    -DMathematica_ROOT_DIR:PATH=/.../mathematica/13.0.1 \
    -DMathematica_WSTP_LIBRARY:FILEPATH=/.../mathematica/13.0.1/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/libWSTP64i4.so \
    -DMathematica_WSTP_INCLUDE_DIR=/.../mathematica/13.0.1/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions \
    -DMathematica_FRONTEND_EXECUTABLE:FILEPATH=/.../mathematica/13.0.1/Executables/mathematica \
    -DMathematica_HOST_ROOT_DIR:PATH=/.../mathematica/13.0.1 \
    -DMathematica_KERNEL_EXECUTABLE:FILEPATH=/.../mathematica/13.0.1/Executables/MathKernel ..
```


### Obtaining the Dataset

After compiling/installing HiggsTools you still need to download the implementation
files for the experimental results. Simply download the respositories for the
[HiggsBounds dataset][HBDataset] and [HiggsSignals dataset][HSDataset] to
convenient locations and you are good to go.


## Examples

### A minimal C++ example: the SM

The 125GeV SM Higgs is implemented as a super simple C++ example and is compiled
automatically. You can try it out by running
```bash
examples/SM /path/to/HBDataSet /path/to/HSDataSet
```
after compiling the C++ library (from the build directory).

### An Artificial Python Example

If you installed the Python module, an artificial minimal example would be. See the [documentation] for details on all the possible production and decay modes that are implemented.


```python
import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS

pred = HP.Predictions() # create the model predictions
bounds = HB.Bounds("/path/to/HBDataset") # load HB dataset
signals = HS.Signals("/path/to/HSDataset") # load HS dataset

# add a SM-like particle
h = pred.addParticle(HP.NeutralScalar("h", "even"))
h.setMass(125.09)
HP.effectiveCouplingInput(h, HP.smLikeEffCouplings)
# evaluate HiggsSignals
chisqSM = signals(pred)

# now give it some lepton-flavor violating decay
# there are very strong limits on this kind of process in HiggsBounds
h.setDecayWidth("emu", 1e-6)

# evaluate HiggsBounds
hbresult = bounds(pred)
print(hbresult)
# evaluate HiggsSignals
chisq = signals(pred)
print(f"HiggsSignals chisq: {chisq} compared to a SM chisq of {chisqSM}")
```


[HiggsBounds-5]: https://gitlab.com/higgsbounds/higgsbounds
[HiggsSignals-2]: https://gitlab.com/higgsbounds/higgssignals
[HBDataset]: https://gitlab.com/higgsbounds/hbdataset
[HSDataset]: https://gitlab.com/higgsbounds/hsdataset
[pip cmake]: https://pypi.org/project/cmake/
[issues]: https://gitlab.com/higgsbounds/higgstools/-/issues
[contact]: mailto:higgstools@desy.de
[documentation]: https://higgsbounds.gitlab.io/higgstools
