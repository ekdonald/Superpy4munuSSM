#include "Higgs.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include <fmt/core.h>
#include <fmt/format.h>
#include <magic_enum.hpp>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <range/v3/core.hpp>
#include <range/v3/view/split.hpp>
#include <range/v3/view/take_last.hpp>

using namespace py::literals;

namespace fmt {
template <> struct formatter<std::complex<double>> : formatter<double> {
    template <typename FormatCtx>
    auto format(const std::complex<double> &x, FormatCtx &ctx)
        -> decltype(ctx.out()) {
        if (x.real() != 0) {
            format_to(ctx.out(), "(");
            formatter<double>::format(x.real(), ctx);
        }
        if (x.real() && x.imag() >= 0)
            format_to(ctx.out(), "+");
        formatter<double>::format(x.imag(), ctx);
        if (x.real() == 0) {
            return format_to(ctx.out(), "j");
        } else {
            return format_to(ctx.out(), "j)");
        }
    }
};
} // namespace fmt

namespace Higgs::predictions {

void bindHiggsPredictions(pybind11::module &higgs) {
    auto hp = higgs.def_submodule(
        "predictions",
        "Stores and helps to obtain theory predictions BSM particles.");

    registerEnum<Collider>(hp);
    registerEnum<Experiment>(hp);
    registerEnum<CP>(hp);
    registerEnum<ECharge>(hp);
    registerEnum<Production>(hp);
    registerEnum<Decay>(hp);
    registerEnum<ChainDecay>(hp);
    registerEnum<MassUncEagerness>(hp);
    registerEnum<ReferenceModel>(hp);
    registerEnum<Coupling>(hp);

    hp.def("validProductionFor", &validProductionFor, "p"_a, "charge"_a);
    hp.def("validDecayFor", &validDecayFor, "d"_a, "charge"_a);
    hp.def("validProductionAt", &validProductionAt, "p"_a, "collider"_a);

    py::class_<Particle>(hp, "Particle")
        .def("id", &Particle::id)
        .def("mass", &Particle::mass)
        .def("setMass", &Particle::setMass)
        .def("massUnc", &Particle::massUnc)
        .def("setMassUnc", &Particle::setMassUnc)
        .def("cp", &Particle::cp)
        .def("charge", &Particle::charge)
        .def(
            "cxn",
            py::overload_cast<Collider, Production>(&Particle::cxn, py::const_),
            "coll"_a, "p"_a)
        .def("br", py::overload_cast<Decay>(&Particle::br, py::const_), "d"_a)
        // this overload manually handles the distinction between chain decays
        // and two BSM decays
        .def(
            "br",
            [](const BsmParticle &self, const std::string &id1,
               const std::string &id2) {
                if (auto chain1 = magic_enum::enum_cast<ChainDecay>(id1);
                    chain1.has_value()) {
                    return self.br(chain1.value(), id2);
                } else if (auto chain2 = magic_enum::enum_cast<ChainDecay>(id2);
                           chain2.has_value()) {
                    return self.br(chain2.value(), id1);
                } else {
                    return self.br(id1, id2);
                }
            },
            "d"_a, "particleId"_a)
        .def("br",
             py::overload_cast<const std::string &, const std::string &>(
                 &Particle::br, py::const_),
             "particleId1"_a, "particleId2"_a)
        .def("totalWidth", &Particle::totalWidth)
        .def("channelRate", &Particle::channelRate, "coll"_a, "p"_a, "d"_a)
        .def("coupling", &Particle::coupling, "c"_a);

    py::class_<BsmParticle, Particle>(hp, "BsmParticle")
        .def(py::init<std::string, ECharge, CP>(), "id"_a, "charge"_a, "cp"_a)
        .def("setCxn", &BsmParticle::setCxn, "coll"_a, "p"_a, "value"_a)
        .def("setNormalizedCxn", &BsmParticle::setNormalizedCxn, "coll"_a,
             "p"_a, "value"_a, "reference"_a)
        .def("setBr", py::overload_cast<Decay, double>(&BsmParticle::setBr),
             "d"_a, "value"_a)
        .def("setBr",
             py::overload_cast<ChainDecay, const std::string &, double>(
                 &BsmParticle::setBr),
             "d"_a, "particleId"_a, "value"_a)
        // this overload manually handles the distinction between chain decays
        // and two BSM decays
        .def(
            "setBr",
            [](BsmParticle &self, const std::string &id1,
               const std::string &id2, double value) {
                if (auto chain1 = magic_enum::enum_cast<ChainDecay>(id1);
                    chain1.has_value()) {
                    return self.setBr(chain1.value(), id2, value);
                } else if (auto chain2 = magic_enum::enum_cast<ChainDecay>(id2);
                           chain2.has_value()) {
                    return self.setBr(chain2.value(), id1, value);
                } else {
                    return self.setBr(id1, id2, value);
                }
            },
            "particleId1"_a, "particleId2"_a, "value"_a)
        .def("setDecayWidth",
             py::overload_cast<Decay, double>(&BsmParticle::setDecayWidth),
             "d"_a, "value"_a)
        .def("setDecayWidth",
             py::overload_cast<ChainDecay, const std::string &, double>(
                 &BsmParticle::setDecayWidth),
             "d"_a, "particleId"_a, "value"_a)
        // this overload manually handles the distinction between chain decays
        // and two BSM decays
        .def(
            "setDecayWidth",
            [](BsmParticle &self, const std::string &id1,
               const std::string &id2, double value) {
                if (auto chain1 = magic_enum::enum_cast<ChainDecay>(id1);
                    chain1.has_value()) {
                    return self.setDecayWidth(chain1.value(), id2, value);
                } else if (auto chain2 = magic_enum::enum_cast<ChainDecay>(id2);
                           chain2.has_value()) {
                    return self.setDecayWidth(chain2.value(), id1, value);
                } else {
                    return self.setDecayWidth(id1, id2, value);
                }
            },
            "particleId1"_a, "particleId2"_a, "value"_a)
        .def(
            "setDecayWidth",
            py::overload_cast<const std::string &, const std::string &, double>(
                &BsmParticle::setDecayWidth),
            "particleId1"_a, "particleId2"_a, "value"_a)
        .def("setTotalWidth", &BsmParticle::setTotalWidth, "value"_a)
        .def("setCoupling", &BsmParticle::setCoupling, "c"_a, "value"_a)
        .def("setChannelRate", &BsmParticle::setChannelRate, "coll"_a, "p"_a,
             "d"_a, "value"_a)
        .def("resetChannelRates", &BsmParticle::resetChannelRates)
        .def("__repr__", [](const BsmParticle &self) {
            return fmt::format(
                "<Higgs.predictions.BsmParticle {0} (CP={1}, e={2}) "
                "with m{0} = {3}+-{4} GeV>",
                self.id(), magic_enum::enum_name(self.cp()), self.charge(),
                self.mass(), self.massUnc());
        });

    hp.def(
        "NeutralScalar",
        [](const std::string &id, CP cp) {
            return BsmParticle(id, ECharge::neutral, cp);
        },
        "id"_a, "cp"_a,
        "Convenience function that returns a new neutral scalar BsmParticle.");

    hp.def(
        "ChargedScalar",
        [](const std::string &id) {
            return BsmParticle(id, ECharge::single, CP::undefined);
        },
        "id"_a,
        "Convenience function that returns a new charged scalar BsmParticle.");

    hp.def(
        "DoublyChargedScalar",
        [](const std::string &id) {
            return BsmParticle(id, ECharge::doubly, CP::undefined);
        },
        "id"_a,
        "Convenience function that returns a new doubly charged scalar "
        "BsmParticle.");

    py::class_<SMHiggs, Particle>(hp, "SMHiggs")
        .def(py::init<double>(), "mass"_a)
        .def("__repr__",
             [](const SMHiggs &h) {
                 return "Higgs.SMHiggs(" + std::to_string(h.mass()) + ")";
             })
        .def("__str__", [](const SMHiggs &h) {
            return "SM-like Higgs with a mass of " + std::to_string(h.mass()) +
                   " GeV";
        });

    py::class_<SMHiggsEW, Particle>(hp, "SMHiggsEW")
        .def(py::init<double>(), "mass"_a)
        .def("__repr__",
             [](const SMHiggsEW &h) {
                 return "Higgs.SMHiggsEW(" + std::to_string(h.mass()) + ")";
             })
        .def("__str__", [](const SMHiggsEW &h) {
            return "SM-like Higgs (incl EW corrections) with a mass of " +
                   std::to_string(h.mass()) + " GeV";
        });

    py::class_<SMHiggsInterp, Particle>(hp, "SMHiggsInterp")
        .def(py::init<double>(), "mass"_a)
        .def("__repr__",
             [](const SMHiggsInterp &h) {
                 return "Higgs.SMHiggsInterp(" + std::to_string(h.mass()) + ")";
             })
        .def("__str__", [](const SMHiggsInterp &h) {
            return "SM-like Higgs (incl EW corrections for low masses) with a "
                   "mass of " +
                   std::to_string(h.mass()) + " GeV";
        });

    py::class_<Predictions>(hp, "Predictions")
        .def(py::init<>())
        .def("addParticle", &Predictions::addParticle, "particle"_a,
             py::return_value_policy::reference_internal)
        .def("particleIds", &Predictions::particleIds)
        .def("particle",
             py::overload_cast<const std::string &>(&Predictions::particle),
             "id"_a, py::return_value_policy::reference_internal)
        .def("particle",
             py::overload_cast<const std::string &>(&Predictions::particle,
                                                    py::const_),
             "id"_a, py::return_value_policy::reference_internal)
        .def("removeParticle", &Predictions::removeParticle)
        .def("setBsmPairCxn", &Predictions::setBsmPairCxn, "collider"_a,
             "id1"_a, "id2"_a, "value"_a)
        .def("bsmPairCxn", &Predictions::bsmPairCxn, "collider"_a, "id1"_a,
             "id2"_a)
        .def("setBrTopWb", &Predictions::setBrTopWb, "value"_a)
        .def("brTopWb", &Predictions::brTopWb);

    py::class_<NeutralEffectiveCouplings>(hp, "NeutralEffectiveCouplings")
        .def(py::init<std::complex<double>, std::complex<double>,
                      std::complex<double>, std::complex<double>,
                      std::complex<double>, std::complex<double>,
                      std::complex<double>, std::complex<double>,
                      std::complex<double>, double, double, double, double,
                      double>(),
             "uu"_a = 0, "dd"_a = 0, "cc"_a = 0, "ss"_a = 0, "tt"_a = 0,
             "bb"_a = 0, "ee"_a = 0, "mumu"_a = 0, "tautau"_a = 0, "WW"_a = 0,
             "ZZ"_a = 0, "Zgam"_a = 0, "gamgam"_a = 0, "gg"_a = 0)
        .def(
            "__repr__",
            [](const NeutralEffectiveCouplings &self) {
                return fmt::format(
                    "Higgs.predictions.NeutralEffectiveCouplings(uu={}, dd={}, "
                    "cc={}, ss={}, tt={}, bb={}, ee={}, mumu={}, tautau={}, "
                    "WW={}, ZZ={}, Zgam={}, gamgam={}, gg={})",
                    self.uu, self.dd, self.cc, self.ss, self.tt, self.bb,
                    self.ee, self.mumu, self.tautau, self.WW, self.ZZ,
                    self.Zgam, self.gamgam, self.gg);
            })
        .def_readwrite("uu", &NeutralEffectiveCouplings::uu)
        .def_readwrite("dd", &NeutralEffectiveCouplings::dd)
        .def_readwrite("cc", &NeutralEffectiveCouplings::cc)
        .def_readwrite("ss", &NeutralEffectiveCouplings::ss)
        .def_readwrite("tt", &NeutralEffectiveCouplings::tt)
        .def_readwrite("bb", &NeutralEffectiveCouplings::bb)
        .def_readwrite("ee", &NeutralEffectiveCouplings::ee)
        .def_readwrite("mumu", &NeutralEffectiveCouplings::mumu)
        .def_readwrite("tautau", &NeutralEffectiveCouplings::tautau)
        .def_readwrite("WW", &NeutralEffectiveCouplings::WW)
        .def_readwrite("ZZ", &NeutralEffectiveCouplings::ZZ)
        .def_readwrite("Zgam", &NeutralEffectiveCouplings::Zgam)
        .def_readwrite("gamgam", &NeutralEffectiveCouplings::gamgam)
        .def_readwrite("gg", &NeutralEffectiveCouplings::gg);

    hp.def("scaledSMlikeEffCouplings", &scaledSMlikeEffCouplings, "scale"_a);
    hp.attr("smLikeEffCouplings") = smLikeEffCouplings;
    hp.def("effectiveCouplingInput", &effectiveCouplingInput, "scalar"_a,
           "coups"_a, "reference"_a = ReferenceModel::SMHiggsInterp,
           "calcggH"_a = true, "calcHgamgam"_a = true);

    auto hpEffCRatios =
        hp.def_submodule("EffectiveCouplingRatios",
                         "Effective coupling parametrisations for cross "
                         "section ratios and derived effective couplings.");

    hpEffCRatios.def("vbfRatio", EffectiveCouplingRatios::vbfRatio,
                     "collider"_a, "mass"_a, "cHZZ"_a, "cHWW"_a);
    hpEffCRatios.def("ttHRatio", EffectiveCouplingRatios::ttHRatio,
                     "collider"_a, "mass"_a, "cHtt"_a);
    hpEffCRatios.def("tHtchanRatio", EffectiveCouplingRatios::tHtchanRatio,
                     "collider"_a, "mass"_a, "cHtt"_a, "cHWW"_a);
    hpEffCRatios.def("tWHRatio", EffectiveCouplingRatios::tWHRatio,
                     "collider"_a, "mass"_a, "cHtt"_a, "cHWW"_a);

    auto hpEffCCxns = hp.def_submodule(
        "EffectiveCouplingCxns",
        "Effective coupling parametrisations for cross cross sections.");

    hpEffCCxns.def("ggH", EffectiveCouplingCxns::ggH, "collider"_a, "mass"_a,
                   "cHtt"_a, "cHbb"_a);
    hpEffCCxns.def("ppHW", EffectiveCouplingCxns::ppHW, "collider"_a, "mass"_a,
                   "cHWW"_a, "cHtt"_a);
    hpEffCCxns.def("ppHZ", EffectiveCouplingCxns::ppHZ, "collider"_a, "mass"_a,
                   "cHZZ"_a, "cHtt"_a, "cHbb"_a);
    hpEffCCxns.def("ggHZ", EffectiveCouplingCxns::ggHZ, "collider"_a, "mass"_a,
                   "cHZZ"_a, "cHtt"_a, "cHbb"_a);
    hpEffCCxns.def("qqHZ", EffectiveCouplingCxns::qqHZ, "collider"_a, "mass"_a,
                   "cHZZ"_a, "cHtt"_a);
    hpEffCCxns.def("bbHZ", EffectiveCouplingCxns::bbHZ, "collider"_a, "mass"_a,
                   "cHbb"_a);
    hpEffCCxns.def("ppHpmtb", EffectiveCouplingCxns::ppHpmtb, "collider"_a,
                   "mass"_a, "cHpmtbR"_a, "cHpmtbL"_a, "brtHpb"_a);
    hpEffCCxns.def("ppHpmPhi", EffectiveCouplingCxns::ppHpmPhi, "collider"_a,
                   "mHpm"_a, "mPhi"_a, "cHpmPhiWmp"_a);
    hpEffCCxns.def("uuHgam", EffectiveCouplingCxns::uuHgam, "collider"_a,
                   "mH"_a, "guu"_a);
    hpEffCCxns.def("ddHgam", EffectiveCouplingCxns::ddHgam, "collider"_a,
                   "mH"_a, "gdd"_a);
    hpEffCCxns.def("ccHgam", EffectiveCouplingCxns::ccHgam, "collider"_a,
                   "mH"_a, "gcc"_a);
    hpEffCCxns.def("ssHgam", EffectiveCouplingCxns::ssHgam, "collider"_a,
                   "mH"_a, "gss"_a);
    hpEffCCxns.def("bbHgam", EffectiveCouplingCxns::bbHgam, "collider"_a,
                   "mH"_a, "gbb"_a);
    hpEffCCxns.def("ucHgam", EffectiveCouplingCxns::ucHgam, "collider"_a,
                   "mH"_a, "guc"_a);
    hpEffCCxns.def("dsHgam", EffectiveCouplingCxns::dsHgam, "collider"_a,
                   "mH"_a, "gds"_a);
    hpEffCCxns.def("dbHgam", EffectiveCouplingCxns::dbHgam, "collider"_a,
                   "mH"_a, "gdb"_a);
    hpEffCCxns.def("dbHgam", EffectiveCouplingCxns::dbHgam, "collider"_a,
                   "mH"_a, "gdb"_a);
    hpEffCCxns.def("sbHgam", EffectiveCouplingCxns::sbHgam, "collider"_a,
                   "mH"_a, "gsb"_a);
    hpEffCCxns.def("udHpgam", EffectiveCouplingCxns::udHpgam, "collider"_a,
                   "mHp"_a, "gLud"_a, "gRud"_a);
    hpEffCCxns.def("csHpgam", EffectiveCouplingCxns::csHpgam, "collider"_a,
                   "mHp"_a, "gLcs"_a, "gRcs"_a);
    hpEffCCxns.def("usHpgam", EffectiveCouplingCxns::usHpgam, "collider"_a,
                   "mHp"_a, "gLus"_a, "gRus"_a);
    hpEffCCxns.def("cdHpgam", EffectiveCouplingCxns::cdHpgam, "collider"_a,
                   "mHp"_a, "gLcd"_a, "gRcd"_a);
    hpEffCCxns.def("ubHpgam", EffectiveCouplingCxns::ubHpgam, "collider"_a,
                   "mHp"_a, "gLub"_a, "gRub"_a);
    hpEffCCxns.def("cbHpgam", EffectiveCouplingCxns::cbHpgam, "collider"_a,
                   "mHp"_a, "gLcb"_a, "gRcb"_a);
    hpEffCCxns.def("udHmgam", EffectiveCouplingCxns::udHmgam, "collider"_a,
                   "mHm"_a, "gLud"_a, "gRud"_a);
    hpEffCCxns.def("csHmgam", EffectiveCouplingCxns::csHmgam, "collider"_a,
                   "mHm"_a, "gLcs"_a, "gRcs"_a);
    hpEffCCxns.def("usHmgam", EffectiveCouplingCxns::usHmgam, "collider"_a,
                   "mHm"_a, "gLus"_a, "gRus"_a);
    hpEffCCxns.def("cdHmgam", EffectiveCouplingCxns::cdHmgam, "collider"_a,
                   "mHm"_a, "gLcd"_a, "gRcd"_a);
    hpEffCCxns.def("ubHmgam", EffectiveCouplingCxns::ubHmgam, "collider"_a,
                   "mHm"_a, "gLub"_a, "gRub"_a);
    hpEffCCxns.def("cbHmgam", EffectiveCouplingCxns::cbHmgam, "collider"_a,
                   "mHm"_a, "gLcb"_a, "gRcb"_a);
    hpEffCCxns.def("uuH", EffectiveCouplingCxns::uuH, "collider"_a, "mH"_a,
                   "guu"_a);
    hpEffCCxns.def("ddH", EffectiveCouplingCxns::ddH, "collider"_a, "mH"_a,
                   "gdd"_a);
    hpEffCCxns.def("ccH", EffectiveCouplingCxns::ccH, "collider"_a, "mH"_a,
                   "gcc"_a);
    hpEffCCxns.def("ssH", EffectiveCouplingCxns::ssH, "collider"_a, "mH"_a,
                   "gss"_a);
    hpEffCCxns.def("ucH", EffectiveCouplingCxns::ucH, "collider"_a, "mH"_a,
                   "guc"_a);
    hpEffCCxns.def("dsH", EffectiveCouplingCxns::dsH, "collider"_a, "mH"_a,
                   "gds"_a);
    hpEffCCxns.def("dbH", EffectiveCouplingCxns::dbH, "collider"_a, "mH"_a,
                   "gdb"_a);
    hpEffCCxns.def("dbH", EffectiveCouplingCxns::dbH, "collider"_a, "mH"_a,
                   "gdb"_a);
    hpEffCCxns.def("sbH", EffectiveCouplingCxns::sbH, "collider"_a, "mH"_a,
                   "gsb"_a);
    hpEffCCxns.def("udHp", EffectiveCouplingCxns::udHp, "collider"_a, "mHp"_a,
                   "gLud"_a, "gRud"_a);
    hpEffCCxns.def("csHp", EffectiveCouplingCxns::csHp, "collider"_a, "mHp"_a,
                   "gLcs"_a, "gRcs"_a);
    hpEffCCxns.def("usHp", EffectiveCouplingCxns::usHp, "collider"_a, "mHp"_a,
                   "gLus"_a, "gRus"_a);
    hpEffCCxns.def("cdHp", EffectiveCouplingCxns::cdHp, "collider"_a, "mHp"_a,
                   "gLcd"_a, "gRcd"_a);
    hpEffCCxns.def("ubHp", EffectiveCouplingCxns::ubHp, "collider"_a, "mHp"_a,
                   "gLub"_a, "gRub"_a);
    hpEffCCxns.def("cbHp", EffectiveCouplingCxns::cbHp, "collider"_a, "mHp"_a,
                   "gLcb"_a, "gRcb"_a);
    hpEffCCxns.def("udHm", EffectiveCouplingCxns::udHm, "collider"_a, "mHm"_a,
                   "gLud"_a, "gRud"_a);
    hpEffCCxns.def("csHm", EffectiveCouplingCxns::csHm, "collider"_a, "mHm"_a,
                   "gLcs"_a, "gRcs"_a);
    hpEffCCxns.def("usHm", EffectiveCouplingCxns::usHm, "collider"_a, "mHm"_a,
                   "gLus"_a, "gRus"_a);
    hpEffCCxns.def("cdHm", EffectiveCouplingCxns::cdHm, "collider"_a, "mHm"_a,
                   "gLcd"_a, "gRcd"_a);
    hpEffCCxns.def("ubHm", EffectiveCouplingCxns::ubHm, "collider"_a, "mHm"_a,
                   "gLub"_a, "gRub"_a);
    hpEffCCxns.def("cbHm", EffectiveCouplingCxns::cbHm, "collider"_a, "mHm"_a,
                   "gLcb"_a, "gRcb"_a);

    auto hpConstants = hp.def_submodule(
        "constants", "Mathematical and physical constants and SM parameters.");

    hpConstants.attr("pi") = py::float_(constants::pi);

    hpConstants.attr("b_Z_inv") = py::float_(constants::b_Z_inv);
    hpConstants.attr("mTop") = py::float_(constants::mTop);
    hpConstants.attr("mBot") = py::float_(constants::mBot);
    hpConstants.attr("mW") = py::float_(constants::mW);
    hpConstants.attr("mMu") = py::float_(constants::mMu);
    hpConstants.attr("mTau") = py::float_(constants::mTau);

    hpConstants.attr("minimumRate") = py::float_(constants::minimumRate);
}

} // namespace Higgs::predictions
