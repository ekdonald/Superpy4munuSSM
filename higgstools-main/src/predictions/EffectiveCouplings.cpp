#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "predictions/CalcHgamgam.hpp"
#include "predictions/Helpers.hpp"
#include "utilities/ArithmeticArray.hpp"
#include "utilities/LinearInterpolator.hpp"
#include "utilities/Logging.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <magic_enum.hpp>
#include <memory>
#include <utility>
#include <vector>

#include "predictions/data/HpmTables.cpp" // NOLINT(bugprone-suspicious-include)
#include "predictions/data/VHTables.cpp"  // NOLINT(bugprone-suspicious-include)
#include "predictions/data/bbHTables.cpp" // NOLINT(bugprone-suspicious-include)
#include "predictions/data/ggHTables.cpp" // NOLINT(bugprone-suspicious-include)
#include "predictions/data/tHTables.cpp"  // NOLINT(bugprone-suspicious-include)
template <class S>
using RefInterp1d = Higgs::utilities::LinearInterpolator<1, S, true>;

namespace Higgs::predictions {

namespace {
void setHadrEffCCxns(BsmParticle &scalar,
                     const NeutralEffectiveCouplings &coups, Collider coll,
                     ReferenceModel reference, bool calcggH) {
    const auto m = scalar.mass();
    if (calcggH == true) {
        scalar.setNormalizedCxn(
            coll, Production::ggH,
            EffectiveCouplingCxns::ggH(coll, m, coups.tt, coups.bb) /
                EffectiveCouplingCxns::ggH(coll, m, 1., 1.),
            reference);
    } else {
        scalar.setNormalizedCxn(coll, Production::ggH, std::pow(coups.gg, 2),
                                reference);
    }
    scalar.setNormalizedCxn(coll, Production::bbH, norm(coups.bb), reference);
    scalar.setNormalizedCxn(
        coll, Production::vbfH,
        EffectiveCouplingRatios::vbfRatio(coll, m, coups.ZZ, coups.WW),
        reference);

    scalar.setNormalizedCxn(
        coll, Production::HW,
        EffectiveCouplingCxns::ppHW(coll, m, coups.WW, coups.tt) /
            EffectiveCouplingCxns::ppHW(coll, m, 1., 1.),
        reference);

    //! These are the same in both of the valid reference models
    scalar.setCxn(coll, Production::qqHZ,
                  EffectiveCouplingCxns::qqHZ(coll, m, coups.ZZ, coups.tt));
    scalar.setCxn(
        coll, Production::ggHZ,
        EffectiveCouplingCxns::ggHZ(coll, m, coups.ZZ, coups.tt, coups.bb));
    scalar.setCxn(coll, Production::bbHZ,
                  EffectiveCouplingCxns::bbHZ(coll, m, coups.bb));

    scalar.setNormalizedCxn(
        coll, Production::Htt,
        EffectiveCouplingRatios::ttHRatio(coll, m, coups.tt), reference);
    scalar.setNormalizedCxn(coll, Production::schanHt, norm(coups.tt),
                            reference);
    scalar.setNormalizedCxn(
        coll, Production::tchanHt,
        EffectiveCouplingRatios::tHtchanRatio(coll, m, coups.tt, coups.WW),
        reference);
    scalar.setNormalizedCxn(
        coll, Production::HtW,
        EffectiveCouplingRatios::tWHRatio(coll, m, coups.tt, coups.WW),
        reference);

    // this relies on the fact that the CP-even/odd couplings don't matter at
    // tree-level
    scalar.setNormalizedCxn(coll, Production::uuHgam, norm(coups.uu),
                            reference);
    scalar.setNormalizedCxn(coll, Production::ddHgam, norm(coups.dd),
                            reference);
    scalar.setNormalizedCxn(coll, Production::ccHgam, norm(coups.cc),
                            reference);
    scalar.setNormalizedCxn(coll, Production::ssHgam, norm(coups.ss),
                            reference);
    scalar.setNormalizedCxn(coll, Production::bbHgam, norm(coups.bb),
                            reference);

    scalar.setCxn(
        coll, Production::qqH,
        EffectiveCouplingCxns::uuH(coll, m, coups.uu * constants::smYup) +
            EffectiveCouplingCxns::ddH(coll, m, coups.dd * constants::smYdown) +
            EffectiveCouplingCxns::ccH(coll, m,
                                       coups.cc * constants::smYcharm) +
            EffectiveCouplingCxns::ssH(coll, m,
                                       coups.ss * constants::smYstrange));
}

void setEffCDecays(BsmParticle &scalar, const NeutralEffectiveCouplings &coups,
                   std::unique_ptr<Particle> refParticle, bool calcggH,
                   bool calcHgamgam) {
    auto scaledDecay = [&scalar, refWidth = refParticle->totalWidth(),
                        ref = std::move(refParticle)](Decay d, double scaling) {
        scalar.setDecayWidth(d, ref->br(d) * refWidth * scaling);
    };
    const auto m = scalar.mass();
    double beta_t =
        1; // fermion decay rate of CP-odd scalar is suppressed by 1/beta^2 with
           // beta = sqrt(1 - 4*mf**2/m**2); only relevant for H -> tt
    if (m > 2 * constants::mTop) {
        beta_t = sqrt(1 - 4 * std::pow(constants::mTop / m, 2));
    };
    scalar.setTotalWidth(0.);
    scaledDecay(Decay::uu, norm(coups.uu));
    scaledDecay(Decay::dd, norm(coups.dd));
    scaledDecay(Decay::cc, norm(coups.cc));
    scaledDecay(Decay::ss, norm(coups.ss));
    scaledDecay(Decay::tt, std::pow(coups.tt.real(), 2) +
                               std::pow(coups.tt.imag() / beta_t, 2));
    scaledDecay(Decay::bb, norm(coups.bb));

    scaledDecay(Decay::ee, norm(coups.ee));
    scaledDecay(Decay::mumu, norm(coups.mumu));
    scaledDecay(Decay::tautau, norm(coups.tautau));

    scaledDecay(Decay::WW, std::pow(coups.WW, 2));
    scaledDecay(Decay::ZZ, std::pow(coups.ZZ, 2));
    scaledDecay(Decay::Zgam, std::pow(coups.Zgam, 2));
    if (calcHgamgam == true) {
        scaledDecay(Decay::gamgam, kgamma2(m, coups));
    } else {
        scaledDecay(Decay::gamgam, std::pow(coups.gamgam, 2));
    }
    if (calcggH == true) {
        scaledDecay(
            Decay::gg,
            EffectiveCouplingCxns::ggH(Collider::LHC13, m, coups.tt, coups.bb) /
                EffectiveCouplingCxns::ggH(Collider::LHC13, m, 1., 1.));
    } else {
        scaledDecay(Decay::gg, std::pow(coups.gg, 2));
    }
}

void setLepCxns(BsmParticle &scalar, const NeutralEffectiveCouplings &coups) {
    scalar.setNormalizedCxn(Collider::LEP, Production::eeHZ,
                            std::pow(coups.ZZ, 2), ReferenceModel::SMHiggs);
    scalar.setNormalizedCxn(Collider::LEP, Production::eeHbb, norm(coups.bb),
                            ReferenceModel::SMHiggs);
    scalar.setNormalizedCxn(Collider::LEP, Production::eeHtautau,
                            norm(coups.tautau), ReferenceModel::SMHiggs);
}

void setCouplings(BsmParticle &scalar, const NeutralEffectiveCouplings &coups) {
    scalar.setCoupling(Coupling::effCPeTopYuk, coups.tt.real());
    scalar.setCoupling(Coupling::effCPoTopYuk, coups.tt.imag());
    scalar.setCoupling(Coupling::effVV, (coups.ZZ + coups.WW) / 2.);
    scalar.setCoupling(Coupling::alphaCPTauYuk,
                       std::signbit(coups.tautau.real())
                           ? std::arg(-coups.tautau)
                           : std::arg(coups.tautau));

    static_assert(magic_enum::enum_count<Coupling>() == 4,
                  "You added a new Coupling, make sure to set it in the "
                  "effective coupling input, if possible.");
}
} // namespace

void effectiveCouplingInput(BsmParticle &scalar,
                            const NeutralEffectiveCouplings &coups,
                            ReferenceModel reference, bool calcggH,
                            bool calcHgamgam) {
    if (scalar.charge() != ECharge::neutral) {
        throw InvalidInput("Cannot use effectiveCouplingInput with neutral "
                           "couplings on charged particle.");
    }
    scalar.resetChannelRates();
    static_assert(
        magic_enum::enum_count<Collider>() == 3,
        "You changed Higgs::predictions::Collider, consider updating "
        "the collider loop in Higgs::predictions::effectiveCouplingInput");

    auto refParticle = getReference(reference, scalar.mass());
    // LCOV_EXCL_START
    if (refParticle->charge() != ECharge::neutral) {
        throw InvalidInput("Cannot use charged reference particle for neutral "
                           "effectiveCouplingInput");
    }
    // LCOV_EXCL_STOP
    for (auto coll : {Collider::LHC8, Collider::LHC13}) {
        setHadrEffCCxns(scalar, coups, coll, reference, calcggH);
    }
    setLepCxns(scalar, coups);
    setEffCDecays(scalar, coups, std::move(refParticle), calcggH, calcHgamgam);
    setCouplings(scalar, coups);
}

namespace EffectiveCouplingRatios {

// ###################################################### //
// ##################### VBF ratios ##################### //
// ###################################################### //
namespace {
double vbfRatioZZLhc8(double mass) {
    static constexpr auto upperMassBound = 1050.;
    mass = std::min(mass, upperMassBound);
    static constexpr auto am1 = 9.99844541408024e-06;
    static constexpr auto a0 = 0.250107980787387;
    static constexpr auto a1 = 6.52202639729497e-05;
    static constexpr auto a2 = -3.20376357582592e-08;
    return am1 / mass + a0 + a1 * mass + a2 * std::pow(mass, 2);
}

double vbfRatioZZLhc13(double mass) {
    static constexpr auto upperMassBound = 3050.;
    mass = std::min(mass, upperMassBound);
    static constexpr auto a0 = 0.263427271363775;
    static constexpr auto a1 = 2.72728035702793e-05;
    static constexpr auto a2 = -3.89549044847535e-09;
    return a0 + a1 * mass + a2 * std::pow(mass, 2);
}
} // namespace

double vbfRatio(Collider coll, double mass, double cHZZ, double cHWW) {
    switch (coll) {
    case Collider::LHC8:
        return vbfRatioZZLhc8(mass) * std::pow(cHZZ, 2) +
               (1 - vbfRatioZZLhc8(mass)) * std::pow(cHWW, 2);
    case Collider::LHC13:
        return vbfRatioZZLhc13(mass) * std::pow(cHZZ, 2) +
               (1 - vbfRatioZZLhc13(mass)) * std::pow(cHWW, 2);
    case Collider::LEP:
        return 0;
    }
    logger()->error("Unknown collider in EffectiveCouplings::vbfRatio()");
    return 0;
}

// ###################################################### //
// ################## ttH and tH ratios ################# //
// ###################################################### //

double ttHRatio(Collider /* coll */, double mass, std::complex<double> cHtt) {
    static const auto interp = RefInterp1d<double>{tHMassGrid, ttHCoeffs};
    return std::pow(cHtt.real(), 2) + interp({mass}) * std::pow(cHtt.imag(), 2);
}

double tHtchanRatio(Collider /* coll */, double mass, std::complex<double> cHtt,
                    double cHWW) {
    static const auto interp =
        RefInterp1d<decltype(tHCoeffs)::value_type>{tHMassGrid, tHCoeffs};
    const auto [ct2, at2, cV2, cvct] = interp({mass});
    return ct2 * std::pow(cHtt.real(), 2) + at2 * std::pow(cHtt.imag(), 2) +
           cV2 * std::pow(cHWW, 2) + cvct * cHtt.real() * cHWW;
}

double tWHRatio(Collider /* coll */, double mass, std::complex<double> cHtt,
                double cHWW) {
    static const auto interp =
        RefInterp1d<decltype(tWHCoeffs)::value_type>{tHMassGrid, tWHCoeffs};
    const auto [ct2, at2, cV2, cvct] = interp({mass});
    return ct2 * std::pow(cHtt.real(), 2) + at2 * std::pow(cHtt.imag(), 2) +
           cV2 * std::pow(cHWW, 2) + cvct * cHtt.real() * cHWW;
}

} // namespace EffectiveCouplingRatios

namespace EffectiveCouplingCxns {

// ------------------------- ggH ------------------------- //

namespace {
class CxnggH {
    static constexpr auto numberOfChannels = 6U;
    using DataType = utilities::ArithmeticArray<double, numberOfChannels>;
    const RefInterp1d<DataType> interp_;

  public:
    explicit CxnggH(const std::vector<DataType> &data)
        : interp_{gridggh, data} {}

    double operator()(double mass, std::complex<double> cHtt,
                      std::complex<double> cHbb) const {
        if (validMassIn<Clamp::lower>(mass, gridggh.front())) {
            const auto &[tt_e, tt_o, bb_e, bb_o, tb_e, tb_o] = interp_({mass});
            return tt_e * pow(cHtt.real(), 2) + tt_o * pow(cHtt.imag(), 2) +
                   bb_e * pow(cHbb.real(), 2) + bb_o * pow(cHbb.imag(), 2) +
                   tb_e * cHtt.real() * cHbb.real() +
                   tb_o * cHtt.imag() * cHbb.imag();
        }
        return 0;
    }
};
} // namespace

double ggH(Collider coll, double mass, std::complex<double> cHtt,
           std::complex<double> cHbb) {
    switch (coll) {
    case Collider::LHC8:
        static const auto cxn8 = CxnggH{lhc8_ggh};
        return cxn8(mass, cHtt, cHbb);
    case Collider::LHC13:
        static const auto cxn13 = CxnggH{lhc13_ggh};
        return cxn13(mass, cHtt, cHbb);
    case Collider::LEP:
        return 0.;
    }
    logger()->error("Unknown collider in EffectiveCouplings::Cxns::ggH()");
    return 0;
}

// ------------------------- HW ------------------------- //
namespace {
class CxnHW {
    static constexpr auto numberOfChannels = 2U;
    using DataType = utilities::ArithmeticArray<double, numberOfChannels>;
    const RefInterp1d<DataType> interp_;

  public:
    explicit CxnHW(const std::vector<DataType> &data) : interp_{whGrid, data} {}

    double operator()(double mass, double cHWW,
                      std::complex<double> cHtt) const {
        if (validMassIn<Clamp::lower>(mass, whGrid.front())) {
            const auto &[WW, tW] = interp_({mass});
            return WW * pow(cHWW, 2) + tW * cHtt.real() * cHWW;
        }
        return 0;
    }
};
} // namespace

double ppHW(Collider coll, double mass, double cHWW,
            std::complex<double> cHtt) {
    switch (coll) {
    case Collider::LHC8:
        static const auto cxn8 = CxnHW{lhc8WHN2lo};
        return cxn8(mass, cHWW, cHtt);
    case Collider::LHC13:
        static const auto cxn13 = CxnHW{lhc13WHN2lo};
        return cxn13(mass, cHWW, cHtt);
    case Collider::LEP:
        return 0.;
    }
    logger()->error("Unknown collider in EffectiveCouplings::Cxns::ppHW()");
    return 0;
}

// ---------------------- gg -> HZ ---------------------- //

namespace {
class CxnggHZ {
    static constexpr auto nloNumberOfChannels = 9U;
    static constexpr auto nnloNumberOfChannels = 6U;
    using NloType = utilities::ArithmeticArray<double, nloNumberOfChannels>;
    using NnloType = utilities::ArithmeticArray<double, nnloNumberOfChannels>;
    const RefInterp1d<NloType> nlo_;
    const RefInterp1d<NnloType> nnlo_;

    [[nodiscard]] double nlo(double mass, double cHZZ,
                             std::complex<double> cHtt,
                             std::complex<double> cHbb) const {
        const auto &[htht, hbhb, hthb, ZZ, Zht, Zhb, atat, abab, atab] =
            nlo_({mass});
        return htht * pow(cHtt.real(), 2) + hbhb * pow(cHbb.real(), 2) +
               hthb * cHtt.real() * cHbb.real() + ZZ * pow(cHZZ, 2) +
               Zht * cHZZ * cHtt.real() + Zhb * cHZZ * cHbb.real() +
               atat * pow(cHtt.imag(), 2) + abab * pow(cHbb.imag(), 2) +
               atab * cHtt.imag() * cHbb.imag();
    }

    [[nodiscard]] double nnlo(double mass, double cHZZ, double cHtt,
                              double cHbb) const {
        const auto &[htht, hbhb, hthb, ZZ, Zht, Zhb] = nnlo_({mass});
        return htht * pow(cHtt, 2) + hbhb * pow(cHbb, 2) + hthb * cHtt * cHbb +
               ZZ * pow(cHZZ, 2) + Zht * cHZZ * cHtt + Zhb * cHZZ * cHbb;
    }

  public:
    CxnggHZ(const std::vector<NloType> &nloData,
            const std::vector<NnloType> &nnloData)
        : nlo_{zhGrid, nloData}, nnlo_{zhGrid, nnloData} {}

    double operator()(double mass, double cHZZ, std::complex<double> cHtt,
                      std::complex<double> cHbb) const {
        if (validMassIn<Clamp::lower>(mass, zhGrid.front())) {
            const auto nloVal = nlo(mass, cHZZ, cHtt, cHbb);
            if (const auto scalarRef =
                    nlo(mass, cHZZ, cHtt.real(), cHbb.real());
                scalarRef > std::numeric_limits<double>::min()) {
                return nloVal * nnlo(mass, cHZZ, cHtt.real(), cHbb.real()) /
                       scalarRef;
            }
            return nloVal;
        }
        return 0;
    }
};
} // namespace

double ggHZ(Collider coll, double mass, double cHZZ, std::complex<double> cHtt,
            std::complex<double> cHbb) {
    switch (coll) {
    case Collider::LHC8:
        static const auto cxn8 = CxnggHZ{lhc8ggZHNlo, lhc8ggZHN2lo};
        return cxn8(mass, cHZZ, cHtt, cHbb);
    case Collider::LHC13:
        static const auto cxn13 = CxnggHZ{lhc13ggZHNlo, lhc13ggZHN2lo};
        return cxn13(mass, cHZZ, cHtt, cHbb);
    case Collider::LEP:
        return 0.;
    }
    logger()->error("Unknown collider in EffectiveCouplings::Cxns::ggHZ()");
    return 0;
}

// ---------------------- qq -> HZ ---------------------- //

namespace {
class CxnqqHZ {
    static constexpr auto numberOfChannels = 2U;
    using DataType = utilities::ArithmeticArray<double, numberOfChannels>;
    const RefInterp1d<DataType> interp_;

  public:
    explicit CxnqqHZ(const std::vector<DataType> &data)
        : interp_{zhGrid, data} {}

    double operator()(double mass, double cHZZ,
                      std::complex<double> cHtt) const {
        if (validMassIn<Clamp::lower>(mass, zhGrid.front())) {
            const auto &[ZZ, Zht] = interp_({mass});
            return ZZ * pow(cHZZ, 2) + Zht * cHZZ * cHtt.real();
        }
        return 0;
    }
};
} // namespace

double qqHZ(Collider coll, double mass, double cHZZ,
            std::complex<double> cHtt) {
    switch (coll) {
    case Collider::LHC8:
        static const auto cxn8 = CxnqqHZ{lhc8qqZHN2lo};
        return cxn8(mass, cHZZ, cHtt);
    case Collider::LHC13:
        static const auto cxn13 = CxnqqHZ{lhc13qqZHN2lo};
        return cxn13(mass, cHZZ, cHtt);
    case Collider::LEP:
        return 0.;
    }
    logger()->error("Unknown collider in EffectiveCouplings::Cxns::qqHZ()");
    return 0;
}

// ---------------------- bb -> HZ ---------------------- //
double bbHZ(Collider coll, double mass, std::complex<double> cHbb) {
    if (validMassIn<Clamp::lower>(mass, zhGrid.front())) {
        switch (coll) {
        case Collider::LHC8:
            static const auto cxn8 = RefInterp1d<double>{zhGrid, lhc8bbZHN2lo};
            return cxn8({mass}) * norm(cHbb);
        case Collider::LHC13:
            static const auto cxn13 =
                RefInterp1d<double>{zhGrid, lhc13bbZHN2lo};
            return cxn13({mass}) * norm(cHbb);
        case Collider::LEP:
            return 0.;
        }
        logger()->error("Unknown collider in EffectiveCouplings::Cxns::bbHZ()");
    }
    return 0;
}

// ------------------ combined pp -> HZ ----------------- //
double ppHZ(Collider coll, double mass, double cHZZ, std::complex<double> cHtt,
            std::complex<double> cHbb) {
    return ggHZ(coll, mass, cHZZ, cHtt, cHbb) + qqHZ(coll, mass, cHZZ, cHtt) +
           bbHZ(coll, mass, cHbb);
}

// ----------------- charged Higgs cxns ----------------- //
double ppHpmtb(Collider coll, double mass, double cHpmtbR, double cHpmtbL,
               double brtHpb) {
    using DataType = utilities::ArithmeticArray<double, 3>;

    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::none>(mass, gridHpmtb.front())) {
        static const auto cxn = RefInterp1d<DataType>{gridHpmtb, lhc13Hpmtb};
        auto [tt, tb, bb] = cxn({mass});

        double widthFactor = std::pow(std::max(1. - brtHpb, 0.), 2);
        return 2 * widthFactor *
               (tt * std::pow(cHpmtbR, 2) + tb * cHpmtbR * cHpmtbL +
                bb * std::pow(cHpmtbL, 2));
    }
    return 0;
}

double ppHpmPhi(Collider coll, double mHpm, double mPhi, double cHpmPhiWmp) {
    using Interp = utilities::LinearInterpolator<2, double, true>;
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::none>(mPhi, gridHpmPhi.front()) &&
        validMassIn<Clamp::none>(mHpm, gridHpmPhi.back())) {
        static const auto cxn = Interp{gridHpmPhi, lhc13HpmPhi};
        return cxn({mPhi, mHpm}) * std::pow(cHpmPhiWmp, 2);
    }
    return 0;
}

// --------------------------  qq -> phi gam ---------------------------- //
double uuHgam(Collider coll, double mH, std::complex<double> guu) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(guu) *
               (-74088400.0 * std::pow(mH, -3) + 4495430.0 * std::pow(mH, -2) -
                6390.61 / mH + 3.69548 - 0.000795657 * mH);
    }
    return 0.;
}

double ddHgam(Collider coll, double mH, std::complex<double> gdd) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gdd) *
               (-26054700.0 * std::pow(mH, -3) + 977832.0 * std::pow(mH, -2) -
                1895.61 / mH + 1.6094 - 0.000534126 * mH);
    }
    return 0.;
}

double ccHgam(Collider coll, double mH, std::complex<double> gcc) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gcc) *
               (+36385600.0 * std::pow(mH, -3) + 207310.0 * std::pow(mH, -2) -
                680.539 / mH + 0.71179 - 0.000256867 * mH);
    }
    return 0.;
}

double ssHgam(Collider coll, double mH, std::complex<double> gss) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gss) *
               (+15763500.0 * std::pow(mH, -3) + 61422.7 * std::pow(mH, -2) -
                184.656 / mH + 0.172986 - 5.59835e-05 * mH);
    }
    return 0.;
}

double bbHgam(Collider coll, double mH, std::complex<double> gbb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gbb) *
               (+567982.0 * std::pow(mH, -3) + 201805.0 * std::pow(mH, -2) -
                448.806 / mH + 0.406157 - 0.000136612 * mH);
    }
    return 0.;
}

double ucHgam(Collider coll, double mH, std::complex<double> guc) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(guc) *
               (+56527600.0 * std::pow(mH, -3) + 3481620.0 * std::pow(mH, -2) -
                6279.74 / mH + 4.49248 - 0.00119019 * mH);
    }
    return 0.;
}

double dsHgam(Collider coll, double mH, std::complex<double> gds) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gds) *
               (+16885300.0 * std::pow(mH, -3) + 699651.0 * std::pow(mH, -2) -
                1510.12 / mH + 1.28013 - 0.00040275 * mH);
    }
    return 0.;
}

double dbHgam(Collider coll, double mH, std::complex<double> gdb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gdb) *
               (+8137770.0 * std::pow(mH, -3) + 330773.0 * std::pow(mH, -2) -
                711.803 / mH + 0.607085 - 0.000192274 * mH);
    }
    return 0.;
}

double sbHgam(Collider coll, double mH, std::complex<double> gsb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gsb) *
               (-3903840.0 * std::pow(mH, -3) + 279829.0 * std::pow(mH, -2) -
                646.173 / mH + 0.613265 - 0.000215673 * mH);
    }
    return 0.;
}

double udHpgam(Collider coll, double mHp, double gLud, double gRud) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLud, 2) + std::pow(gRud, 2)) *
               (-17338100.0 * std::pow(mHp, -3) + 652027.0 * std::pow(mHp, -2) -
                1223.17 / mHp + 0.98188 - 0.000307538 * mHp);
    }
    return 0.;
}

double csHpgam(Collider coll, double mHp, double gLcs, double gRcs) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLcs, 2) + std::pow(gRcs, 2)) *
               (+5332200.0 * std::pow(mHp, -3) + 52331.9 * std::pow(mHp, -2) -
                136.529 / mHp + 0.125393 - 4.09939e-05 * mHp);
    }
    return 0.;
}

double usHpgam(Collider coll, double mHp, double gLus, double gRus) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLus, 2) + std::pow(gRus, 2)) *
               (+3710340.0 * std::pow(mHp, -3) + 142254.0 * std::pow(mHp, -2) -
                351.939 / mHp + 0.335395 - 0.000116551 * mHp);
    }
    return 0.;
}

double cdHpgam(Collider coll, double mHp, double gLcd, double gRcd) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLcd, 2) + std::pow(gRcd, 2)) *
               (-15052000.0 * std::pow(mHp, -3) + 457490.0 * std::pow(mHp, -2) -
                1006.88 / mHp + 0.917357 - 0.000313157 * mHp);
    }
    return 0.;
}

double ubHpgam(Collider coll, double mHp, double gLub, double gRub) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLub, 2) + std::pow(gRub, 2)) *
               (+2709810.0 * std::pow(mHp, -3) + 72679.9 * std::pow(mHp, -2) -
                190.285 / mHp + 0.182345 - 6.24154e-05 * mHp);
    }
    return 0.;
}

double cbHpgam(Collider coll, double mHp, double gLcb, double gRcb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLcb, 2) + std::pow(gRcb, 2)) *
               (+2847960.0 * std::pow(mHp, -3) + 28145.8 * std::pow(mHp, -2) -
                82.7795 / mHp + 0.0822505 - 2.85173e-05 * mHp);
    }
    return 0.;
}

double udHmgam(Collider coll, double mHm, double gLud, double gRud) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLud, 2) + std::pow(gRud, 2)) *
               (-471225.0 * std::pow(mHm, -3) + 445930.0 * std::pow(mHm, -2) -
                247.076 / mHm - 0.216092 + 0.000173311 * mHm);
    }
    return 0.;
}

double csHmgam(Collider coll, double mHm, double gLcs, double gRcs) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLcs, 2) + std::pow(gRcs, 2)) *
               (+6778240.0 * std::pow(mHm, -3) + 34604.7 * std::pow(mHm, -2) -
                104.975 / mHm + 0.102024 - 3.47606e-05 * mHm);
    }
    return 0.;
}

double usHmgam(Collider coll, double mHm, double gLus, double gRus) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLus, 2) + std::pow(gRus, 2)) *
               (+2837360.0 * std::pow(mHm, -3) + 296278.0 * std::pow(mHm, -2) -
                329.227 / mHm + 0.0805123 + 2.16609e-05 * mHm);
    }
    return 0.;
}

double cdHmgam(Collider coll, double mHm, double gLcd, double gRcd) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLcd, 2) + std::pow(gRcd, 2)) *
               (+2874780.0 * std::pow(mHm, -3) + 148855.0 * std::pow(mHm, -2) -
                395.107 / mHm + 0.39282 - 0.000140121 * mHm);
    }
    return 0.;
}

double ubHmgam(Collider coll, double mHm, double gLub, double gRub) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLub, 2) + std::pow(gRub, 2)) *
               (-24880100.0 * std::pow(mHm, -3) + 396466.0 * std::pow(mHm, -2) -
                815.75 / mHm + 0.775128 - 0.000288148 * mHm);
    }
    return 0.;
}

double cbHmgam(Collider coll, double mHm, double gLcb, double gRcb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLcb, 2) + std::pow(gRcb, 2)) *
               (+3896500.0 * std::pow(mHm, -3) + 19429.3 * std::pow(mHm, -2) -
                57.5982 / mHm + 0.0525333 - 1.6354e-05 * mHm);
    }
    return 0.;
}

// --------------------------  qq -> phi ----------------------------- //
double uuH(Collider coll, double mH, std::complex<double> guu) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(guu) * (+90150900000.0 * std::pow(mH, -3) +
                            41882700.0 * std::pow(mH, -2) - 245875.0 / mH +
                            254.458 - 0.0908537 * mH);
    }
    return 0.;
}

double ddH(Collider coll, double mH, std::complex<double> gdd) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gdd) * (+70783300000.0 * std::pow(mH, -3) -
                            38427600.0 * std::pow(mH, -2) - 42109.9 / mH +
                            58.2871 - 0.0205486 * mH);
    }
    return 0.;
}

double ccH(Collider coll, double mH, std::complex<double> gcc) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gcc) * (+15053800000.0 * std::pow(mH, -3) -
                            43835300.0 * std::pow(mH, -2) + 59437.3 / mH -
                            40.8521 + 0.0113045 * mH);
    }
    return 0.;
}

double ssH(Collider coll, double mH, std::complex<double> gss) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gss) * (+23486200000.0 * std::pow(mH, -3) -
                            66938000.0 * std::pow(mH, -2) + 93341.3 / mH -
                            67.2994 + 0.0196796 * mH);
    }
    return 0.;
}

double ucH(Collider coll, double mH, std::complex<double> guc) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(guc) * (+90959700000.0 * std::pow(mH, -3) -
                            75285200.0 * std::pow(mH, -2) - 37365.5 / mH +
                            87.893 - 0.0371667 * mH);
    }
    return 0.;
}

double dsH(Collider coll, double mH, std::complex<double> gds) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gds) * (+172313000000.0 * std::pow(mH, -3) -
                            175715000.0 * std::pow(mH, -2) - 3836.18 / mH +
                            114.545 - 0.0551281 * mH);
    }
    return 0.;
}

double dbH(Collider coll, double mH, std::complex<double> gdb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gdb) * (+41957900000.0 * std::pow(mH, -3) -
                            18564000.0 * std::pow(mH, -2) - 60170.9 / mH +
                            87.2187 - 0.0351971 * mH);
    }
    return 0.;
}

double sbH(Collider coll, double mH, std::complex<double> gsb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mH, {200., 1150.})) {
        return norm(gsb) * (+26589800000.0 * std::pow(mH, -3) -
                            9680920.0 * std::pow(mH, -2) - 41054.6 / mH +
                            58.3277 - 0.0238107 * mH);
    }
    return 0.;
}

double udHp(Collider coll, double mHp, double gLud, double gRud) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLud, 2) + std::pow(gRud, 2)) *
               (+45284000000.0 * std::pow(mHp, -3) +
                28699600.0 * std::pow(mHp, -2) - 132718.0 / mHp + 132.009 -
                0.0460763 * mHp);
    }
    return 0.;
}

double csHp(Collider coll, double mHp, double gLcs, double gRcs) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLcs, 2) + std::pow(gRcs, 2)) *
               (+8912850000.0 * std::pow(mHp, -3) -
                25080200.0 * std::pow(mHp, -2) + 30539.2 / mHp - 17.502 +
                0.00400584 * mHp);
    }
    return 0.;
}

double usHp(Collider coll, double mHp, double gLus, double gRus) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLus, 2) + std::pow(gRus, 2)) *
               (+40927300000.0 * std::pow(mHp, -3) -
                25000400.0 * std::pow(mHp, -2) - 21620.2 / mHp + 26.8752 -
                0.00701031 * mHp);
    }
    return 0.;
}

double cdHp(Collider coll, double mHp, double gLcd, double gRcd) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLcd, 2) + std::pow(gRcd, 2)) *
               (+14306200000.0 * std::pow(mHp, -3) -
                33913900.0 * std::pow(mHp, -2) + 36824.0 / mHp - 20.178 +
                0.00443296 * mHp);
    }
    return 0.;
}

double ubHp(Collider coll, double mHp, double gLub, double gRub) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLub, 2) + std::pow(gRub, 2)) *
               (+25075400000.0 * std::pow(mHp, -3) -
                2584920.0 * std::pow(mHp, -2) - 60789.7 / mHp + 77.2716 -
                0.0302666 * mHp);
    }
    return 0.;
}

double cbHp(Collider coll, double mHp, double gLcb, double gRcb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHp, {200., 1150.})) {
        return (std::pow(gLcb, 2) + std::pow(gRcb, 2)) *
               (+5221740000.0 * std::pow(mHp, -3) -
                14868400.0 * std::pow(mHp, -2) + 18402.6 / mHp - 10.7985 +
                0.00237351 * mHp);
    }
    return 0.;
}

double udHm(Collider coll, double mHm, double gLud, double gRud) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLud, 2) + std::pow(gRud, 2)) *
               (+33554900000.0 * std::pow(mHm, -3) -
                18198100.0 * std::pow(mHm, -2) - 30947.2 / mHm + 44.7557 -
                0.0172376 * mHm);
    }
    return 0.;
}

double csHm(Collider coll, double mHm, double gLcs, double gRcs) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLcs, 2) + std::pow(gRcs, 2)) *
               (+9772690000.0 * std::pow(mHm, -3) -
                29140300.0 * std::pow(mHm, -2) + 42584.6 / mHm - 32.2039 +
                0.00982901 * mHm);
    }
    return 0.;
}

double usHm(Collider coll, double mHm, double gLus, double gRus) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLus, 2) + std::pow(gRus, 2)) *
               (+16897600000.0 * std::pow(mHm, -3) -
                42210300.0 * std::pow(mHm, -2) + 51821.8 / mHm - 33.1287 +
                0.00863466 * mHm);
    }
    return 0.;
}

double cdHm(Collider coll, double mHm, double gLcd, double gRcd) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLcd, 2) + std::pow(gRcd, 2)) *
               (+22547700000.0 * std::pow(mHm, -3) -
                22707500.0 * std::pow(mHm, -2) - 9033.18 / mHm + 27.6082 -
                0.0128189 * mHm);
    }
    return 0.;
}

double ubHm(Collider coll, double mHm, double gLub, double gRub) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLub, 2) + std::pow(gRub, 2)) *
               (+9811400000.0 * std::pow(mHm, -3) -
                25465300.0 * std::pow(mHm, -2) + 30948.2 / mHm - 19.4331 +
                0.00502026 * mHm);
    }
    return 0.;
}

double cbHm(Collider coll, double mHm, double gLcb, double gRcb) {
    if (coll == Collider::LHC13 &&
        validMassIn<Clamp::lower>(mHm, {200., 1150.})) {
        return (std::pow(gLcb, 2) + std::pow(gRcb, 2)) *
               (+5466740000.0 * std::pow(mHm, -3) -
                16991200.0 * std::pow(mHm, -2) + 24637.7 / mHm - 18.2229 +
                0.00543889 * mHm);
    }
    return 0.;
}

} // namespace EffectiveCouplingCxns
} // namespace Higgs::predictions
