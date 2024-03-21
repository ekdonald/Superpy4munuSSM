/**
 * @file EffectiveCouplings.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Set particle properties using the effective coupling approximation.
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include <complex>

namespace Higgs {
namespace predictions {

class BsmParticle;

//! Effective coupling of a neutral particle normalized to a SM-like Higgs
//! boson. The real part of the fermionic couplings is the CP-even, and the
//! imaginary part the CP-odd coupling (both normalized to the SM-like CP-even
//! coupling).
struct HIGGSTOOLS_EXPORT NeutralEffectiveCouplings {
    //! \f$\kappa_u + i \kappa_\tilde{u} \f$
    std::complex<double> uu = 0;
    //! \f$\kappa_d + i \kappa_\tilde{d} \f$
    std::complex<double> dd = 0;
    //! \f$\kappa_c + i \kappa_\tilde{c} \f$
    std::complex<double> cc = 0;
    //! \f$\kappa_s + i \kappa_\tilde{s} \f$
    std::complex<double> ss = 0;
    //! \f$\kappa_t + i \kappa_\tilde{t} \f$
    std::complex<double> tt = 0;
    //! \f$\kappa_b + i \kappa_\tilde{b} \f$
    std::complex<double> bb = 0;
    //! \f$\kappa_e + i \kappa_\tilde{e} \f$
    std::complex<double> ee = 0;
    //! \f$\kappa_\mu + i \kappa_\tilde{\mu} \f$
    std::complex<double> mumu = 0;
    //! \f$\kappa_\tau + i \kappa_\tilde{\tau} \f$
    std::complex<double> tautau = 0;

    double WW = 0;     //!< \f$ \kappa_W \f$
    double ZZ = 0;     //!< \f$ \kappa_Z \f$
    double Zgam = 0;   //!< \f$ \kappa_{Z\gamma} \f$
    double gamgam = 0; //!< \f$ \kappa_\gamma \f$
    double gg = 0;     //!< \f$ \kappa_g \f$
};

//! SM-like effective couplings with a global scaling.
HIGGSTOOLS_EXPORT constexpr NeutralEffectiveCouplings
scaledSMlikeEffCouplings(double scale) {
    return NeutralEffectiveCouplings{scale, scale, scale, scale, scale,
                                     scale, scale, scale, scale, scale,
                                     scale, scale, scale, scale};
}

//! values of the effective couplings for a SM-like Higgs boson
constexpr auto smLikeEffCouplings = scaledSMlikeEffCouplings(1.);

/**
 * @brief Set the cxns and BRs of the particle based on the given effective
 * couplings.
 *
 * @throws InvalidInput if the scalar particle or the reference particle is not
 * neutral.
 *
 * @param scalar a neutral scalar particle
 * @param coups the effective couplings
 * @param reference the reference model, either #ReferenceModel::SMHiggs,
 * #ReferenceModel::SMHiggsEW, or #ReferenceModel::SMHiggsInterp.
 * @param calcggH whether to calculate the ggH cross-section in terms of the
 * effective top and bottom Yukawa couplings or by rescaling the SM-like ggH XS
 * by the squared of the effective gg coupling (no effects from colored BSM particles
 * are taken into account).
 * @param calcHgamgam whether to calculate the H->gaga decay width in terms of the
 * effective couplings or by rescaling the SM-like H->gaga decay
 * by the squared of the effective gamgam coupling (no effects from charged BSM particles
 * are taken into account).
 */
HIGGSTOOLS_EXPORT void
effectiveCouplingInput(BsmParticle &scalar,
                       const NeutralEffectiveCouplings &coups,
                       ReferenceModel reference = ReferenceModel::SMHiggsInterp,
                       bool calcggH = true, bool calcHgamgam = true);

/**
 * Functions for calculating derived effective couplings or normalized cross
 * sections from the basic effective couplings.
 *
 * The following effective couplings appear:
 *
 * parameter name | effective coupling
 * ---------------|-------------------
 * `cHWW`         | \f$\kappa_W\f$
 * `cHZZ`         | \f$\kappa_Z\f$
 * `cHtt`         | \f$\kappa_t + i \kappa_\tilde{t}\f$
 * `cHbb`         | \f$\kappa_b + i \kappa_\tilde{b}\f$
 *
 * where the \f$\kappa\f$ are defined in eq (6) of the HB5 manual.
 *
 * The VBF ratios were obtained from a LO calculation with HAWK
 * [1412.5390](https://arxiv.org/abs/1412.5390) and generalized to allow for
 * different WW and ZZ effective couplings. Interference effects were found to
 * be small and neglected. The ratios are clamped for too low and too large
 * masses.
 *
 * The ttH, tH and tWH ratios were obtained from a LO MG5 calculation using the
 * Higgs characterization model (performed by Henning Bahl, thanks!) and
 * generalize Eqs (15-17) of [2007.08542](https://arxiv.org/abs/2007.08542) to
 * masses other than 125 GeV. They were calculated for the 13TeV LHC but are
 * also used for 8TeV if required.
 */
namespace EffectiveCouplingRatios {

//! \f$ \sigma(\mathrm{VBF}\to H)/\sigma(\mathrm{VBF}\to h_\mathrm{SM}) \f$
HIGGSTOOLS_EXPORT double vbfRatio(Collider coll, double mass, double cHZZ,
                                  double cHWW);
//! \f$ \sigma(pp\to t\bar{t}H)/\sigma(pp\to t\bar{t}h_\mathrm{SM}) \f$
HIGGSTOOLS_EXPORT double ttHRatio(Collider coll, double mass,
                                  std::complex<double> cHtt);
//! \f$ \sigma(pp\to tH)/\sigma(pp\to t h_\mathrm{SM}) \f$ t-channel
HIGGSTOOLS_EXPORT double tHtchanRatio(Collider coll, double mass,
                                      std::complex<double> cHtt, double cHWW);
//! \f$ \sigma(pp\to tWH)/\sigma(pp\to tW h_\mathrm{SM}) \f$
HIGGSTOOLS_EXPORT double tWHRatio(Collider coll, double mass,
                                  std::complex<double> cHtt, double cHWW);

} // namespace EffectiveCouplingRatios

/**
 * Functions for cross section predictions in the effective coupling
 * approximation.
 *
 * All of the returned cross section values are in pb.
 *
 * The following effective couplings of a neutral Higgs relative to the SM-Higgs
 * appear
 *
 * parameter name | effective coupling
 * ---------------|-------------------
 * `cHWW`         | \f$\kappa_W\f$
 * `cHZZ`         | \f$\kappa_Z\f$
 * `cHtt`         | \f$\kappa_t + i \kappa_\tilde{t}\f$
 * `cHbb`         | \f$\kappa_b + i \kappa_\tilde{b}\f$
 *
 * where the \f$\kappa\f$ are defined in eq (6) of the HB5 manual. Includes the
 * WH and ZH cross sections including subchannels as described in Section 4.1 of
 * the HB5 manual. The neutral cross sections cover the ranges
 *
 * cross section          |   min  |  max
 * -----------------------|--------|---------
 * ggH                    | 10 GeV | 3000 GeV
 * ppWH                   |  1 GeV | 2950 GeV
 * ppZH, ggZH, qqZH, bbZH |  1 GeV | 4999 GeV
 *
 * for larger masses the cxn is set to 0, for smaller masses the cross section
 * at 1GeV is returned as a conservative estimate.
 *
 * The charged Higgs cross sections in association with top quarks and neutral
 * Higgs bosons are documented in detail below.
 *
 * The production cross sections of a neutral or charged scalar from quarks
 * including those in association with a photon are computed in the generic
 * scalar model defined in 2109.10366 Eqs. (1-2) as a function of the particle
 * mass and the respective (absolute) coupling values `g`. They cover the mass
 * range of [200GeV, 1150GeV] are zero for larger and clamped for smaller
 * masses.
 *
 */
namespace EffectiveCouplingCxns {

//! \f$ gg\to H\f$ at the specified collider
HIGGSTOOLS_EXPORT double ggH(Collider coll, double mass,
                             std::complex<double> cHtt,
                             std::complex<double> cHbb);
//! \f$ pp\to H W^\pm \f$ at the specified collider
HIGGSTOOLS_EXPORT double ppHW(Collider coll, double mass, double cHWW,
                              std::complex<double> cHtt);
//! \f$ pp\to H Z = q\bar{q}/gg/b\bar{b}\to HZ\f$ at the specified collider
HIGGSTOOLS_EXPORT double ppHZ(Collider coll, double mass, double cHZZ,
                              std::complex<double> cHtt,
                              std::complex<double> cHbb);
//! \f$ gg\to H Z \f$ at the specified collider
HIGGSTOOLS_EXPORT double ggHZ(Collider coll, double mass, double cHZZ,
                              std::complex<double> cHtt,
                              std::complex<double> cHbb);
//! \f$ q\bar{q}\to H Z \f$ at the specified collider
HIGGSTOOLS_EXPORT double qqHZ(Collider coll, double mass, double cHZZ,
                              std::complex<double> cHtt);
//! \f$ b\bar{b}\to H Z \f$ at the specified collider
HIGGSTOOLS_EXPORT double bbHZ(Collider coll, double mass,
                              std::complex<double> cHbb);

/**
 * @brief \f$ pp \to H^\pm t (b)\f$
 *
 * 4FS/5FS and intermediate mass range matched results from
 * https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGMSSMCharged as
 * described in detail in the HB5 manual.
 *
 * @param coll the collider, this cxn is only available for the 13TeV LHC and
 * returns 0 otherwise
 * @param mass the charged Higgs mass, the tabulated cross section covers the
 * mass range 145GeV-2TeV. For masses outside this range returns 0. This makes
 * sense also at lower masses, since the process is replaced by on-shell decays
 * of top-quarks which are modeled through
 * #Higgs::predictions::Production::brtHpb.
 * @param cHpmtbR \f$\kappa^\pm_t\f$ as defined in eq (13) of
 * [2006.06007](https://arxiv.org/abs/2006.06007)
 * @param cHpmtbL \f$\kappa^\pm_b\f$ as defined in eq (13) of
 * [2006.06007](https://arxiv.org/abs/2006.06007)
 * @param brtHpb \f$\mathrm{BR}(t\to H^+b)\f$ where \f$H^+\f$ is the charged
 * Higgs for which the cxn is requested. If this is non-zero, it is assumed that
 * this is the only non-SM decay mode of the top-quark
 * @returns double the cross section in pb
 */
HIGGSTOOLS_EXPORT double ppHpmtb(Collider coll, double mass, double cHpmtbR,
                                 double cHpmtbL, double brtHpb);

/**
 * @brief \f$pp\to H^\pm \phi\f$ where \f$\phi\f$ is a neutral scalar particle
 *
 * Approximate (K-factor) NNLO cross sections as described in
 * [2103.07484](https://arxiv.org/abs/2103.07484).
 *
 * @param coll the collider, this cxn is only available for the 13TeV LHC and
 * returns 0 otherwise
 * @param mHpm the charged Higgs mass, the tabulated cross section covers the
 * range 100GeV-500GeV.
 * @param mPhi  the neutral Higgs mass, the tabulated cross section covers the
 * range 10GeV-500GeV.
 * @param cHpmPhiWmp \f$g(H^\pm_i W^\mp h_j/a_j)/\frac{g}{2}\f$ as defined in
 * eqs (10, 11) of [2103.07484](https://arxiv.org/abs/2103.07484) and \f$\phi\f$
 * is identified with either of \f$h_j/a_j\f$.
 * @return double the cross section in pb. Returns 0 for masses outsize the
 * covered ranges.
 */
HIGGSTOOLS_EXPORT double ppHpmPhi(Collider coll, double mHpm, double mPhi,
                                  double cHpmPhiWmp);

//! \f$u\bar{u} \to H \gamma\f$
HIGGSTOOLS_EXPORT double uuHgam(Collider coll, double mH,
                                std::complex<double> guu);
//! \f$d\bar{d} \to H \gamma\f$
HIGGSTOOLS_EXPORT double ddHgam(Collider coll, double mH,
                                std::complex<double> gdd);
//! \f$c\bar{c} \to H \gamma\f$
HIGGSTOOLS_EXPORT double ccHgam(Collider coll, double mH,
                                std::complex<double> gcc);
//! \f$s\bar{s} \to H \gamma\f$
HIGGSTOOLS_EXPORT double ssHgam(Collider coll, double mH,
                                std::complex<double> gss);
//! \f$b\bar{b} \to H \gamma\f$
HIGGSTOOLS_EXPORT double bbHgam(Collider coll, double mH,
                                std::complex<double> gbb);
//! \f$u\bar{c}+\bar{u}c \to H \gamma\f$
HIGGSTOOLS_EXPORT double ucHgam(Collider coll, double mH,
                                std::complex<double> guc);
//! \f$d\bar{s}+\bar{d}s \to H \gamma\f$
HIGGSTOOLS_EXPORT double dsHgam(Collider coll, double mH,
                                std::complex<double> gds);
//! \f$d\bar{b}+\bar{d}b \to H \gamma\f$
HIGGSTOOLS_EXPORT double dbHgam(Collider coll, double mH,
                                std::complex<double> gdb);
//! \f$s\bar{b}+\bar{s}b \to H \gamma\f$
HIGGSTOOLS_EXPORT double sbHgam(Collider coll, double mH,
                                std::complex<double> gsb);
//! \f$u\bar{d} \to H^+ \gamma\f$
HIGGSTOOLS_EXPORT double udHpgam(Collider coll, double mHp, double gLud,
                                 double gRud);
//! \f$c\bar{s} \to H^+ \gamma\f$
HIGGSTOOLS_EXPORT double csHpgam(Collider coll, double mHp, double gLcs,
                                 double gRcs);
//! \f$u\bar{s} \to H^+ \gamma\f$
HIGGSTOOLS_EXPORT double usHpgam(Collider coll, double mHp, double gLus,
                                 double gRus);
//! \f$c\bar{d} \to H^+ \gamma\f$
HIGGSTOOLS_EXPORT double cdHpgam(Collider coll, double mHp, double gLcd,
                                 double gRcd);
//! \f$u\bar{b} \to H^+ \gamma\f$
HIGGSTOOLS_EXPORT double ubHpgam(Collider coll, double mHp, double gLub,
                                 double gRub);
//! \f$c\bar{b} \to H^+ \gamma\f$
HIGGSTOOLS_EXPORT double cbHpgam(Collider coll, double mHp, double gLcb,
                                 double gRcb);
//! \f$\bar{u}d \to H^- \gamma\f$
HIGGSTOOLS_EXPORT double udHmgam(Collider coll, double mHm, double gLud,
                                 double gRud);
//! \f$\bar{c}s \to H^- \gamma\f$
HIGGSTOOLS_EXPORT double csHmgam(Collider coll, double mHm, double gLcs,
                                 double gRcs);
//! \f$\bar{u}s \to H^- \gamma\f$
HIGGSTOOLS_EXPORT double usHmgam(Collider coll, double mHm, double gLus,
                                 double gRus);
//! \f$\bar{c}d \to H^- \gamma\f$
HIGGSTOOLS_EXPORT double cdHmgam(Collider coll, double mHm, double gLcd,
                                 double gRcd);
//! \f$\bar{u}b \to H^- \gamma\f$
HIGGSTOOLS_EXPORT double ubHmgam(Collider coll, double mHm, double gLub,
                                 double gRub);
//! \f$\bar{c}b \to H^- \gamma\f$
HIGGSTOOLS_EXPORT double cbHmgam(Collider coll, double mHm, double gLcb,
                                 double gRcb);
//! \f$u\bar{u} \to H \f$
HIGGSTOOLS_EXPORT double uuH(Collider coll, double mH,
                             std::complex<double> guu);
//! \f$d\bar{d} \to H \f$
HIGGSTOOLS_EXPORT double ddH(Collider coll, double mH,
                             std::complex<double> gdd);
//! \f$c\bar{c} \to H \f$
HIGGSTOOLS_EXPORT double ccH(Collider coll, double mH,
                             std::complex<double> gcc);
//! \f$s\bar{s} \to H \f$
HIGGSTOOLS_EXPORT double ssH(Collider coll, double mH,
                             std::complex<double> gss);
//! \f$u\bar{c}+\bar{u}c \to H \f$
HIGGSTOOLS_EXPORT double ucH(Collider coll, double mH,
                             std::complex<double> guc);
//! \f$d\bar{s}+\bar{d}s \to H \f$
HIGGSTOOLS_EXPORT double dsH(Collider coll, double mH,
                             std::complex<double> gds);
//! \f$d\bar{b}+\bar{d}b \to H \f$
HIGGSTOOLS_EXPORT double dbH(Collider coll, double mH,
                             std::complex<double> gdb);
//! \f$s\bar{b}+\bar{s}b \to H \f$
HIGGSTOOLS_EXPORT double sbH(Collider coll, double mH,
                             std::complex<double> gsb);
//! \f$u\bar{d} \to H^+ \f$
HIGGSTOOLS_EXPORT double udHp(Collider coll, double mHp, double gLud,
                              double gRud);
//! \f$c\bar{s} \to H^+ \f$
HIGGSTOOLS_EXPORT double csHp(Collider coll, double mHp, double gLcs,
                              double gRcs);
//! \f$u\bar{s} \to H^+ \f$
HIGGSTOOLS_EXPORT double usHp(Collider coll, double mHp, double gLus,
                              double gRus);
//! \f$c\bar{d} \to H^+ \f$
HIGGSTOOLS_EXPORT double cdHp(Collider coll, double mHp, double gLcd,
                              double gRcd);
//! \f$u\bar{b} \to H^+ \f$
HIGGSTOOLS_EXPORT double ubHp(Collider coll, double mHp, double gLub,
                              double gRub);
//! \f$c\bar{b} \to H^+ \f$
HIGGSTOOLS_EXPORT double cbHp(Collider coll, double mHp, double gLcb,
                              double gRcb);
//! \f$\bar{u}d \to H^- \f$
HIGGSTOOLS_EXPORT double udHm(Collider coll, double mHm, double gLud,
                              double gRud);
//! \f$\bar{c}s \to H^- \f$
HIGGSTOOLS_EXPORT double csHm(Collider coll, double mHm, double gLcs,
                              double gRcs);
//! \f$\bar{u}s \to H^- \f$
HIGGSTOOLS_EXPORT double usHm(Collider coll, double mHm, double gLus,
                              double gRus);
//! \f$\bar{c}d \to H^- \f$
HIGGSTOOLS_EXPORT double cdHm(Collider coll, double mHm, double gLcd,
                              double gRcd);
//! \f$\bar{u}b \to H^- \f$
HIGGSTOOLS_EXPORT double ubHm(Collider coll, double mHm, double gLub,
                              double gRub);
//! \f$\bar{c}b \to H^- \f$
HIGGSTOOLS_EXPORT double cbHm(Collider coll, double mHm, double gLcb,
                              double gRcb);

} // namespace EffectiveCouplingCxns

} // namespace predictions
} // namespace Higgs
