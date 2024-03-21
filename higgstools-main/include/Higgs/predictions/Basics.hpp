/**
 * @file Basics.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Enumerations, constants, errors, and other basic utilities used
 * throughout the library
 *
 * @copyright Copyright 2020 by the authors. This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once
#include "Higgs/HiggsTools_export.h"
#include <stdexcept>

namespace Higgs {
namespace predictions {
enum class Production;
enum class Decay;

//! All implemented Colliders.
enum class Collider {
    LEP = 1,    //!< LEP, \f$e^+e^-\f$ at 90-210 GeV
    LHC8 = 8,   //!< LHC run 1, \f$pp\f$ at 8 TeV
    LHC13 = 13, //!< LHC run 2, \f$pp\f$ at 13 TeV
};

//! The type of collider
enum class ColliderType {
    ee, //!< e+e- collider
    pp, //!< hadron collider
};

//! Determine whether the given collider is a hadron or a lepton collider.
HIGGSTOOLS_EXPORT ColliderType classifyCollider(Collider coll) noexcept;

//! Experimental collaborations at the various colliders.
enum class Experiment {
    ATLAS,   //!< the ATLAS detector at the LHC
    CMS,     //!< the CMS detector at the LHC
    LHCComb, //!< combined ATLAS/CMS results
    ALEPH,   //!< the ALEPH detector at LEP
    DELPHI,  //!< the DELPHI detector at LEP
    L3,      //!< the L3 detector at LEP
    OPAL,    //!< the OPAL detector at LEP
    LEPComb, //!< combined ALEPH/DELPHI/L3/OPAL results
};

//! CP quantum numbers
enum class CP {
    even = 1,     //!< CP-even
    odd = -1,     //!< CP-odd
    undefined = 0 //!< not a CP eigenstate or unspecified
};

//! electric charge
enum class ECharge {
    neutral, //!< neutral
    single,  //!< singly charged
    doubly   //!< doubly charged
};

//! Specify the level of eagerness used to perform operations that depend on
//! the theoretical mass uncertainty.
enum class MassUncEagerness {
    cautious, //!< only if certain
    eager,    //!< whenever possible
    ignore    //!< use central values
};

//! Mathematical and physical constants, SM parameters, and global
//! configurations parameters.
namespace constants {
//! \f$ \pi \f$
static constexpr double pi = 3.14159265358979323846;
//! the reduced EW vev \f$ v/\sqrt{2} \f$ in GeV
static constexpr double vEWred = 174.10358025986375;

//! \f$ \mathrm{BR}(Z\to\mathrm{inv}) \f$ from PDG 2020
static constexpr double b_Z_inv = 0.2;
//! on-shell \f$ m_t \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mTop = 172.5;
//! on-shell \f$ m_b \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mBot = 4.92;
//! on-shell \f$ m_c \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mCharm = 1.51;
//! strange quark mass from PDG 2020
static constexpr double mStrange = 93e-3;
//! strange quark mass from PDG 2020
static constexpr double mUp = 2.16e-3;
//! strange quark mass from PDG 2020
static constexpr double mDown = 4.67e-3;
//! estimate of the SM charm-quark Yukawa from the on-shell mass given in the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double smYcharm = mCharm / vEWred;
//! estimate of the SM strange-quark Yukawa from the PDG 2020 current-quark mass
static constexpr double smYstrange = mStrange / vEWred;
//! estimate of the SM up-quark Yukawa from the PDG 2020 current-quark mass
static constexpr double smYup = mUp / vEWred;
//! estimate of the SM down-quark Yukawa from the PDG 2020 current-quark mass
static constexpr double smYdown = mDown / vEWred;

//! on-shell \f$ m_W \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mW = 80.385;
//! \f$ m_e \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mEl = 0.510998928e-3;
//! \f$ m_\mu \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mMu = 105.6583715e-3;
//! \f$ m_\tau \f$ in GeV from the
//! [LHC recommendations](https://cds.cern.ch/record/2047636)
static constexpr double mTau = 1.77682;

//! Rates below this value (in pb) are neglected. For absolute rates this is
//! 0.1zb (zepto barn) which is well below anything thinkable at the moment. For
//! normalized rates this is still tiny, since e.g. SM Higgs rates are never
//! above a few 100 pb.
static constexpr double minimumRate = 1e-10;
} // namespace constants

//! Thrown if invalid or inconsistent user input is detected.
class HIGGSTOOLS_EXPORT InvalidInput : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

//! Thrown when a Channel is encountered that violates basic consistency
//! requirements.
class HIGGSTOOLS_EXPORT InvalidChannel : public std::runtime_error {
  public:
    //! electric charge violation in decay mode
    InvalidChannel(ECharge charge, Decay d) noexcept;
    //! electric charge violation in production mode
    InvalidChannel(ECharge charge, Production p) noexcept;
    //! invalid production mode for collider type
    InvalidChannel(ColliderType collType, Production p) noexcept;
    //! mismatch between charges for production and decay mode
    InvalidChannel(Production p, Decay d) noexcept;
};

} // namespace predictions
} // namespace Higgs
