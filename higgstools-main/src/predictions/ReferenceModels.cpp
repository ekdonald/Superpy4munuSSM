#include "Higgs/predictions/ReferenceModels.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "predictions/Helpers.hpp"
#include "utilities/LinearInterpolator.hpp"
#include "utilities/Logging.hpp"
#include <array>
#include <memory>

#include "predictions/data/YR4Brs.cpp"
#include "predictions/data/YR4Cxns.cpp"
#include "predictions/data/tHTables.cpp"

namespace Higgs::predictions {

SMHiggs::SMHiggs(double mass) noexcept : SMHiggs{mass, "hSM"} {}

SMHiggs::SMHiggs(double mass, std::string id) noexcept
    : Particle{std::move(id), CP::even, ECharge::neutral} {
    setMass(mass);
}

std::unique_ptr<Particle> SMHiggs::clone() const {
    return std::make_unique<SMHiggs>(*this);
}

SMHiggsEW::SMHiggsEW(double mass) noexcept : SMHiggs{mass, "hSMEW"} {}

std::unique_ptr<Particle> SMHiggsEW::clone() const {
    return std::make_unique<SMHiggsEW>(*this);
}

SMHiggsInterp::SMHiggsInterp(double mass) noexcept : SMHiggs{mass, "hSMInterp"} {}

std::unique_ptr<Particle> SMHiggsInterp::clone() const {
    return std::make_unique<SMHiggsInterp>(*this);
}
namespace {
using RefInterp1d = Higgs::utilities::LinearInterpolator<1, double, true>;

double ggF(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, bsmGridCxns.front())) {
        switch (coll) {
        case Collider::LHC8:
            static const auto cxn8 = RefInterp1d{bsmGridCxns, lhc8ggF};
            return cxn8({mass});
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{bsmGridCxns, lhc13ggFN2lo};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0;
}

double bbH(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, bsmGridCxns.front())) {
        switch (coll) {
        case Collider::LHC8:
            static const auto cxn8 = RefInterp1d{bsmGridCxns, lhc8bbH};
            return cxn8({mass});
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{bsmGridCxns, lhc13bbH};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0;
}

double vbf(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, bsmGridCxns.front())) {
        switch (coll) {
        case Collider::LHC8:
            static const auto cxn8 = RefInterp1d{bsmGridCxns, lhc8vbf};
            return cxn8({mass});
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{bsmGridCxns, lhc13vbf};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0;
}

double Htt(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, HttGrid.front())) {
        switch (coll) {
        case Collider::LHC8:
            static const auto cxn8 = RefInterp1d{HttGrid, lhc8Htt};
            return cxn8({mass});
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{HttGrid, lhc13Htt};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0.;
}

double tchanHt(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, HtGrid.front())) {
        switch (coll) {
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{HtGrid, lhc13tchanHt};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0.;
}

double schanHt(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, HtGrid.front())) {
        switch (coll) {
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{HtGrid, lhc13schanHt};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0;
}

double HtW(Collider coll, double mass) {
    if (validMassIn<Clamp::lower>(mass, tHMassGrid.front())) {
        switch (coll) {
        case Collider::LHC13:
            static const auto cxn13 = RefInterp1d{tHMassGrid, tWHCxnsLHC13};
            return cxn13({mass});
        default:
            return 0.;
        }
    }
    return 0.;
}
} // namespace

double SMHiggs::cxn(Collider coll, Production p) const noexcept {
    switch (p) {
    default:
        return 0.;
    case Production::H:
        return cxn(coll, Production::ggH) + cxn(coll, Production::bbH) +
               cxn(coll, Production::qqH);
    case Production::ggH:
        return ggF(coll, mass());
    case Production::bbH:
        return bbH(coll, mass());
    case Production::vbfH:
        return vbf(coll, mass());
    case Production::HW:
        return EffectiveCouplingCxns::ppHW(coll, mass(), 1., 1.);
    case Production::HZ:
        return EffectiveCouplingCxns::ppHZ(coll, mass(), 1., 1., 1.);
    case Production::qqHZ:
        return EffectiveCouplingCxns::qqHZ(coll, mass(), 1., 1.);
    case Production::ggHZ:
        return EffectiveCouplingCxns::ggHZ(coll, mass(), 1., 1., 1.);
    case Production::bbHZ:
        return EffectiveCouplingCxns::bbHZ(coll, mass(), 1.);
    case Production::Htt:
        return Htt(coll, mass());
    case Production::tchanHt:
        return tchanHt(coll, mass());
    case Production::schanHt:
        return schanHt(coll, mass());
    case Production::Ht:
        return cxn(coll, Production::tchanHt) + cxn(coll, Production::schanHt);
    case Production::HtW:
        return HtW(coll, mass());
    case Production::qqH:
        return EffectiveCouplingCxns::uuH(coll, mass(), constants::smYup) +
               EffectiveCouplingCxns::ddH(coll, mass(), constants::smYdown) +
               EffectiveCouplingCxns::ccH(coll, mass(), constants::smYcharm) +
               EffectiveCouplingCxns::ssH(coll, mass(), constants::smYstrange);
    case Production::uuHgam:
        return EffectiveCouplingCxns::uuHgam(coll, mass(), constants::smYup);
    case Production::ddHgam:
        return EffectiveCouplingCxns::ddHgam(coll, mass(), constants::smYdown);
    case Production::ccHgam:
        return EffectiveCouplingCxns::ccHgam(coll, mass(), constants::smYcharm);
    case Production::ssHgam:
        return EffectiveCouplingCxns::ssHgam(coll, mass(),
                                             constants::smYstrange);
    case Production::bbHgam:
        return EffectiveCouplingCxns::bbHgam(
            coll, mass(), constants::mBot / constants::vEWred);
    case Production::eeHZ:
    case Production::eeHbb:
    case Production::eeHtautau:
        return coll == Collider::LEP ? 1. : 0.;
    }
}

namespace {
double ggFN3LOLHC13(double mass) {
    if (validMassIn<Clamp::lower>(mass, bsmGridCxns.front())) {
        static const auto cxn13 = RefInterp1d{bsmGridCxns, lhc13ggFN3lo};
        return cxn13({mass});
    }
    return 0;
}
} // namespace

double SMHiggsEW::cxn(Collider coll, Production p) const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::none>(std::get<0>(m), smGridCxns.front())) {
        switch (p) {
        case Production::ggH: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8ggFSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13ggFSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::bbH: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8bbHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13bbHSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::vbfH: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8VBFSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13VBFSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::HW: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8WHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13WHSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::Htt: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8ttHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13ttHSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::tchanHt: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 =
                    RefInterp1d{smGridCxns, lhc8tH_tchanSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 =
                    RefInterp1d{smGridCxns, lhc13tH_tchanSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::schanHt: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 =
                    RefInterp1d{smGridCxns, lhc8tH_schanSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 =
                    RefInterp1d{smGridCxns, lhc13tH_schanSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::HtW: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8tWHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13tWHSMEW};
                return cxn13(m);
            }
            break;
        }
        default:
            break;
        }
    }
    if (p == Production::ggH && coll == Collider::LHC13) {
        return ggFN3LOLHC13(mass());
    }
    return SMHiggs::cxn(coll, p);
}

namespace {
double ggFinterpLHC13(double mass) {
    if (validMassIn<Clamp::lower>(mass, bsmGridCxns.front())) {
        static const auto cxn13 = RefInterp1d{bsmGridCxns, lhc13ggFinterp};
        return cxn13({mass});
    }
    return 0;
}
} // namespace

double SMHiggsInterp::cxn(Collider coll, Production p) const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::none>(std::get<0>(m), smGridCxns.front())) {
        switch (p) {
        case Production::ggH: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8ggFSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13ggFSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::bbH: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8bbHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13bbHSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::vbfH: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8VBFSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13VBFSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::HW: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8WHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13WHSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::Htt: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8ttHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13ttHSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::tchanHt: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 =
                    RefInterp1d{smGridCxns, lhc8tH_tchanSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 =
                    RefInterp1d{smGridCxns, lhc13tH_tchanSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::schanHt: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 =
                    RefInterp1d{smGridCxns, lhc8tH_schanSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 =
                    RefInterp1d{smGridCxns, lhc13tH_schanSMEW};
                return cxn13(m);
            }
            break;
        }
        case Production::HtW: {
            if (coll == Collider::LHC8) {
                static const auto cxn8 = RefInterp1d{smGridCxns, lhc8tWHSMEW};
                return cxn8(m);
            } else if (coll == Collider::LHC13) {
                static const auto cxn13 = RefInterp1d{smGridCxns, lhc13tWHSMEW};
                return cxn13(m);
            }
            break;
        }
        default:
            break;
        }
    }
    if (p == Production::ggH && coll == Collider::LHC13) {
        return ggFinterpLHC13(mass());
    }
    return SMHiggs::cxn(coll, p);
}

double SMHiggs::br(Decay d) const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::none>(std::get<0>(m), bsmGridBrs.front())) {
        switch (d) {
        default:
            return 0.;
        case Decay::cc:
            static const auto cc = RefInterp1d{bsmGridBrs, brcc};
            return cc(m);
        case Decay::ss:
            static const auto ss = RefInterp1d{bsmGridBrs, brss};
            return ss(m);
        case Decay::tt:
            static const auto tt = RefInterp1d{bsmGridBrs, brtt};
            return tt(m);
        case Decay::bb:
            static const auto bb = RefInterp1d{bsmGridBrs, brbb};
            return bb(m);
        case Decay::mumu:
            static const auto mumu = RefInterp1d{bsmGridBrs, brmumu};
            return mumu(m);
        case Decay::tautau:
            static const auto tautau = RefInterp1d{bsmGridBrs, brtautau};
            return tautau(m);
        case Decay::WW:
            static const auto WW = RefInterp1d{bsmGridBrs, brWW};
            return WW(m);
        case Decay::ZZ:
            static const auto ZZ = RefInterp1d{bsmGridBrs, brZZ};
            return ZZ(m);
        case Decay::Zgam:
            static const auto Zgam = RefInterp1d{bsmGridBrs, brZgam};
            return Zgam(m);
        case Decay::gamgam:
            static const auto gamgam = RefInterp1d{bsmGridBrs, brgamgam};
            return gamgam(m);
        case Decay::gg:
            static const auto gg = RefInterp1d{bsmGridBrs, brgg};
            return gg(m);
        case Decay::inv:
            return br(Decay::ZZ) * std::pow(constants::b_Z_inv, 2);
        }
    }
    return 0;
}

double SMHiggsEW::br(Decay d) const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::none>(std::get<0>(m), smGridBrs.front())) {
        switch (d) {
        case Decay::cc:
            static const auto cc = RefInterp1d{smGridBrs, brccSM};
            return cc(m);
        case Decay::ss:
            static const auto ss = RefInterp1d{smGridBrs, brssSM};
            return ss(m);
        case Decay::bb:
            static const auto bb = RefInterp1d{smGridBrs, brbbSM};
            return bb(m);
        case Decay::mumu:
            static const auto mumu = RefInterp1d{smGridBrs, brmumuSM};
            return mumu(m);
        case Decay::tautau:
            static const auto tautau = RefInterp1d{smGridBrs, brtautauSM};
            return tautau(m);
        case Decay::WW:
            static const auto WW = RefInterp1d{smGridBrs, brWWSM};
            return WW(m);
        case Decay::ZZ:
            static const auto ZZ = RefInterp1d{smGridBrs, brZZSM};
            return ZZ(m);
        case Decay::Zgam:
            static const auto Zgam = RefInterp1d{smGridBrs, brZgamSM};
            return Zgam(m);
        case Decay::gamgam:
            static const auto gamgam = RefInterp1d{smGridBrs, brgamgamSM};
            return gamgam(m);
        case Decay::gg:
            static const auto gg = RefInterp1d{smGridBrs, brggSM};
            return gg(m);
        default:
            break;
        }
    }
    return SMHiggs::br(d);
}

double SMHiggsInterp::br(Decay d) const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::none>(std::get<0>(m), smGridBrs.front())) {
        switch (d) {
        case Decay::cc:
            static const auto cc = RefInterp1d{smGridBrs, brccSM};
            return cc(m);
        case Decay::ss:
            static const auto ss = RefInterp1d{smGridBrs, brssSM};
            return ss(m);
        case Decay::bb:
            static const auto bb = RefInterp1d{smGridBrs, brbbSM};
            return bb(m);
        case Decay::mumu:
            static const auto mumu = RefInterp1d{smGridBrs, brmumuSM};
            return mumu(m);
        case Decay::tautau:
            static const auto tautau = RefInterp1d{smGridBrs, brtautauSM};
            return tautau(m);
        case Decay::WW:
            static const auto WW = RefInterp1d{smGridBrs, brWWSM};
            return WW(m);
        case Decay::ZZ:
            static const auto ZZ = RefInterp1d{smGridBrs, brZZSM};
            return ZZ(m);
        case Decay::Zgam:
            static const auto Zgam = RefInterp1d{smGridBrs, brZgamSM};
            return Zgam(m);
        case Decay::gamgam:
            static const auto gamgam = RefInterp1d{smGridBrs, brgamgamSM};
            return gamgam(m);
        case Decay::gg:
            static const auto gg = RefInterp1d{smGridBrs, brggSM};
            return gg(m);
        default:
            break;
        }
    }
    return SMHiggs::br(d);
}

double SMHiggs::totalWidth() const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::both>(std::get<0>(m), bsmGridBrs.front())) {
        static const auto gamtot = RefInterp1d{bsmGridBrs, wtot};
        return gamtot(m);
    }
    return 0.;
}

double SMHiggsEW::totalWidth() const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::both>(std::get<0>(m), smGridBrs.front())) {
        static const auto gamtot = RefInterp1d{smGridBrs, wtotSM};
        return gamtot(m);
    }
    return SMHiggs::totalWidth();
}

double SMHiggsInterp::totalWidth() const noexcept {
    if (auto m = std::array{mass()};
        validMassIn<Clamp::both>(std::get<0>(m), smGridBrs.front())) {
        static const auto gamtot = RefInterp1d{smGridBrs, wtotSM};
        return gamtot(m);
    }
    return SMHiggs::totalWidth();
}

std::optional<double> SMHiggs::coupling(Coupling c) const noexcept {
    switch (c) {
    case Coupling::effCPeTopYuk:
    case Coupling::effVV:
        return 1.;
    case Coupling::effCPoTopYuk:
    case Coupling::alphaCPTauYuk:
        return 0.;
    }
    return std::nullopt; // LCOV_EXCL_LINE
}

std::unique_ptr<Particle> getReference(ReferenceModel model, double mass) {
    switch (model) {
    case ReferenceModel::SMHiggs:
        return std::make_unique<SMHiggs>(mass);
    case ReferenceModel::SMHiggsEW:
        return std::make_unique<SMHiggsEW>(mass);
    case ReferenceModel::SMHiggsInterp:
        return std::make_unique<SMHiggsInterp>(mass);
    }
    throw(std::runtime_error{"Unknown ReferenceModel in getReference()"});
}

} // namespace Higgs::predictions
