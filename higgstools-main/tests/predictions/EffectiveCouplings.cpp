#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <complex>
#include <functional>
#include <initializer_list>
#include <magic_enum.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/remove.hpp>
#include <string>
#include <string_view>

#include "data/VHTables.cpp"

using namespace std::complex_literals;
using Catch::Approx;
namespace HP = Higgs::predictions;
namespace {
constexpr auto colliders = magic_enum::enum_values<HP::Collider>();
constexpr auto prodChannels = magic_enum::enum_values<HP::Production>();
constexpr auto decChannels = magic_enum::enum_values<HP::Decay>();
} // namespace

TEST_CASE("Neutral EffC input -- cross sections") {
    SECTION("SM-like") {
        auto ref =
            GENERATE(HP::ReferenceModel::SMHiggs, HP::ReferenceModel::SMHiggsEW,
                     HP::ReferenceModel::SMHiggsInterp);
        auto sm =
            HP::BsmParticle("hsm", HP::ECharge::neutral, HP::CP::undefined);
        sm.setMass(125);
        effectiveCouplingInput(sm, HP::smLikeEffCouplings, ref);
        auto smlike =
            HP::BsmParticle("Hsm", HP::ECharge::neutral, HP::CP::undefined);
        smlike.setMass(300);
        effectiveCouplingInput(smlike, HP::smLikeEffCouplings, ref);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                for (auto sref : {std::cref(sm), std::cref(smlike)}) {
                    const auto &s = sref.get();
                    INFO(magic_enum::enum_name(coll)
                         << " " << magic_enum::enum_name(p)
                         << " mass: " << s.mass());
                    if (ref == HP::ReferenceModel::SMHiggs) {
                        CHECK(s.cxn(coll, static_cast<HP::Production>(p)) ==
                              Approx(HP::SMHiggs(s.mass()).cxn(
                                  coll, static_cast<HP::Production>(p))));
                    } else if (ref == HP::ReferenceModel::SMHiggsEW) {
                        CHECK(s.cxn(coll, static_cast<HP::Production>(p)) ==
                              Approx(HP::SMHiggsEW(s.mass()).cxn(
                                  coll, static_cast<HP::Production>(p))));
                    } else if (ref == HP::ReferenceModel::SMHiggsInterp) {
                        CHECK(s.cxn(coll, static_cast<HP::Production>(p)) ==
                              Approx(HP::SMHiggsInterp(s.mass()).cxn(
                                  coll, static_cast<HP::Production>(p))));
                    }
                }
            }
        }
    }

    SECTION("global scaling") {
        auto s =
            HP::BsmParticle("hsm", HP::ECharge::neutral, HP::CP::undefined);
        s.setMass(150);
        auto scaledEffCoups = HP::NeutralEffectiveCouplings{
            2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.};
        effectiveCouplingInput(s, scaledEffCoups, HP::ReferenceModel::SMHiggs);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                INFO(magic_enum::enum_name(coll)
                     << " " << magic_enum::enum_name(p));
                CHECK(s.cxn(coll, static_cast<HP::Production>(p)) ==
                      Approx(4 * HP::SMHiggs(s.mass()).cxn(
                                     coll, static_cast<HP::Production>(p))));
            }
        }
    }

    SECTION("charge mismatch") {
        auto hp = HP::BsmParticle("h+", HP::ECharge::single, HP::CP::undefined);
        REQUIRE_THROWS_AS(effectiveCouplingInput(hp, HP::smLikeEffCouplings),
                          HP::InvalidInput);
    }
}

TEST_CASE("Neutral EffC input --- decays") {
    SECTION("SM-like") {
        auto sm =
            HP::BsmParticle("hsm", HP::ECharge::neutral, HP::CP::undefined);
        sm.setMass(125);
        effectiveCouplingInput(sm, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs);
        auto smlike =
            HP::BsmParticle("Hsm", HP::ECharge::neutral, HP::CP::undefined);
        smlike.setMass(300);
        effectiveCouplingInput(smlike, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs);
        for (auto sref : {std::cref(sm), std::cref(smlike)}) {
            const auto &s = sref.get();
            for (auto d : decChannels) {
                INFO(magic_enum::enum_name(d) << " mass: " << s.mass());
                CHECK(s.br(static_cast<HP::Decay>(d)) ==
                      Approx(
                          HP::SMHiggs(s.mass()).br(static_cast<HP::Decay>(d))));
            }
            CHECK(s.totalWidth() == Approx(HP::SMHiggs(s.mass()).totalWidth()));
        }
    }
    SECTION("global scaling") {
        auto s = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        s.setMass(150);
        auto scaledEffCoups1 = HP::NeutralEffectiveCouplings{
            2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.};
        auto scaledEffCoups2 =
            HP::NeutralEffectiveCouplings{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                          0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        for (const auto &effc : {scaledEffCoups1, scaledEffCoups2}) {
            effectiveCouplingInput(s, effc, HP::ReferenceModel::SMHiggs);
            for (auto d : decChannels) {
                INFO(magic_enum::enum_name(d));
                CHECK(s.br(static_cast<HP::Decay>(d)) ==
                      Approx(
                          HP::SMHiggs(s.mass()).br(static_cast<HP::Decay>(d))));
            }
            CHECK(s.totalWidth() == Approx(effc.ZZ * effc.ZZ *
                                           HP::SMHiggs(s.mass()).totalWidth()));
        }
    }
}

TEST_CASE("neutral effC input --- couplings") {
    auto ref =
        GENERATE(HP::ReferenceModel::SMHiggs, HP::ReferenceModel::SMHiggsEW,
                 HP::ReferenceModel::SMHiggsInterp);
    auto smlike =
        HP::BsmParticle("Hsm", HP::ECharge::neutral, HP::CP::undefined);
    smlike.setMass(400);
    effectiveCouplingInput(smlike, HP::smLikeEffCouplings, ref);
    for (auto c : magic_enum::enum_values<HP::Coupling>()) {
        REQUIRE(smlike.coupling(c) == HP::SMHiggs(10).coupling(c));
    }

    SECTION("tau CP-phase in range") {
        auto coups = HP::smLikeEffCouplings;
        coups.tautau =
            GENERATE(3. + 2.i, -0.4 + 0.6i, 0.5 - 1.1i, -1.4 - 2.2i, 0. + 0i,
                     -0. - 0.i, 1. + 0.i, 1. - 0.i, -1. + 0.i, -1. - 0.i,
                     0. + 1i, 0. - 1i, -0. + 1i, -0. - 1i);
        effectiveCouplingInput(smlike, coups, ref);
        INFO(coups.tautau);
        CHECK(smlike.coupling(HP::Coupling::alphaCPTauYuk) >=
              -HP::constants::pi / 2.);
        CHECK(smlike.coupling(HP::Coupling::alphaCPTauYuk) <=
              HP::constants::pi / 2.);
    }

    SECTION("beta factor for H->tt") {
        auto coups = HP::NeutralEffectiveCouplings();
        coups.tt = 1.;
        effectiveCouplingInput(smlike, coups, ref, false, false);
        auto dw1 = smlike.br(HP::Decay::tt) * smlike.totalWidth();
        coups.tt = 1i;
        effectiveCouplingInput(smlike, coups, ref);
        auto dw2 = smlike.br(HP::Decay::tt) * smlike.totalWidth();
        CHECK(
            dw1 ==
            Approx(dw2 *
                   (1 - 4 * std::pow(HP::constants::mTop / smlike.mass(), 2))));
    }
}

TEST_CASE("VBF ratios") {
    using HP::EffectiveCouplingRatios::vbfRatio;
    SECTION("LHC8") {
        auto mass = GENERATE(10, 125, 1000, 10000);
        REQUIRE(vbfRatio(HP::Collider::LHC8, mass, 1., 1.) == Approx(1.));
        REQUIRE(vbfRatio(HP::Collider::LHC8, mass, 0., 0.) == Approx(0.));
        auto cZZ = GENERATE(-10, -2, -1, 0, 1, 2, 10);
        auto cWW = GENERATE(-10, -2, -1, 0, 1, 2, 10);
        INFO("mass: " << mass << " cZZ: " << cZZ << " cWW: " << cWW);
        CHECK(vbfRatio(HP::Collider::LHC8, mass, cWW, cZZ) >= 0);
    }
    SECTION("LHC13") {
        auto mass = GENERATE(10, 125, 1000, 10000);
        REQUIRE(vbfRatio(HP::Collider::LHC13, mass, 1., 1.) == Approx(1.));
        REQUIRE(vbfRatio(HP::Collider::LHC13, mass, 0., 0.) == Approx(0.));
        auto cZZ = GENERATE(-10, -2, -1, 0, 1, 2, 10);
        auto cWW = GENERATE(-10, -2, -1, 0, 1, 2, 10);
        INFO("mass: " << mass << " cZZ: " << cZZ << " cWW: " << cWW);
        CHECK(vbfRatio(HP::Collider::LHC13, mass, cWW, cZZ) >= 0);
    }

    SECTION("invalid collider") {
        REQUIRE(vbfRatio(HP::Collider::LEP, 125, 1, 1) == 0);
        REQUIRE(vbfRatio(HP::Collider{-1}, 125, 1, 1) == 0);
    }
}

TEST_CASE("ttH ratios") {
    using HP::EffectiveCouplingRatios::ttHRatio;
    SECTION("LHC8") {
        auto mass = GENERATE(10, 125, 1000);
        REQUIRE(ttHRatio(HP::Collider::LHC8, mass, 1.) == Approx(1.));
        REQUIRE(ttHRatio(HP::Collider::LHC8, mass, 0.) == Approx(0.));
        auto stt = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto att = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        INFO("mass: " << mass << " stt: " << stt << " att: " << att);
        CHECK(ttHRatio(HP::Collider::LHC8, mass, {stt, att}) >= 0);
    }
    SECTION("LHC13") {
        auto mass = GENERATE(10, 125, 1000);
        REQUIRE(ttHRatio(HP::Collider::LHC13, mass, 1.) == Approx(1.));
        REQUIRE(ttHRatio(HP::Collider::LHC13, mass, 0.) == Approx(0.));
        auto stt = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto att = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        INFO("mass: " << mass << " stt: " << stt << " att: " << att);
        CHECK(ttHRatio(HP::Collider::LHC13, mass, {stt, att}) >= 0);
    }
}

TEST_CASE("tH-tchan ratios") {
    using HP::EffectiveCouplingRatios::tHtchanRatio;
    SECTION("LHC8") {
        auto mass = GENERATE(10, 125, 1000);
        REQUIRE(tHtchanRatio(HP::Collider::LHC8, mass, 1., 1.) == Approx(1.));
        REQUIRE(tHtchanRatio(HP::Collider::LHC8, mass, 0., 0.) == Approx(0.));
        auto stt = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto att = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto cWW = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        INFO("mass: " << mass << " stt: " << stt << " att: " << att
                      << " cWW: " << cWW);
        CHECK(tHtchanRatio(HP::Collider::LHC8, mass, {stt, att}, cWW) >= 0);
    }
    SECTION("LHC13") {
        auto mass = GENERATE(10, 125, 1000);
        REQUIRE(tHtchanRatio(HP::Collider::LHC13, mass, 1., 1.) == Approx(1.));
        REQUIRE(tHtchanRatio(HP::Collider::LHC13, mass, 0., 0.) == Approx(0.));
        auto stt = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto att = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto cWW = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        INFO("mass: " << mass << " stt: " << stt << " att: " << att
                      << " cWW: " << cWW);
        CHECK(tHtchanRatio(HP::Collider::LHC13, mass, {stt, att}, cWW) >= 0);
    }
}

TEST_CASE("tWH ratios") {
    using HP::EffectiveCouplingRatios::tWHRatio;
    SECTION("LHC8") {
        auto mass = GENERATE(10, 125, 1000);
        REQUIRE(tWHRatio(HP::Collider::LHC8, mass, 1., 1.) == Approx(1.));
        REQUIRE(tWHRatio(HP::Collider::LHC8, mass, 0., 0.) == Approx(0.));
        auto stt = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto att = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto cWW = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        INFO("mass: " << mass << " stt: " << stt << " att: " << att
                      << " cWW: " << cWW);
        CHECK(tWHRatio(HP::Collider::LHC8, mass, {stt, att}, cWW) >= 0);
    }
    SECTION("LHC13") {
        auto mass = GENERATE(10, 125, 1000);
        REQUIRE(tWHRatio(HP::Collider::LHC13, mass, 1., 1.) == Approx(1.));
        REQUIRE(tWHRatio(HP::Collider::LHC13, mass, 0., 0.) == Approx(0.));
        auto stt = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto att = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        auto cWW = GENERATE(-10., -2., -1., 0., 1., 2., 10.);
        INFO("mass: " << mass << " stt: " << stt << " att: " << att
                      << " cWW: " << cWW);
        CHECK(tWHRatio(HP::Collider::LHC13, mass, {stt, att}, cWW) >= 0);
    }
}

TEST_CASE("ggH cxns") {
    SECTION("LHC8") {
        REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC8, 35.,
                                               1. + 1.i, 1. + 1.i) ==
                Approx(784.8956459999999));
    }
    SECTION("LHC13") {
        REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC13, 35.,
                                               1. + 1.i, 1. + 1.i) ==
                Approx(1428.8072739999998));
    }
    SECTION("LEP") {
        REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LEP, 100., 1.,
                                               1.) == 0);
    }
    SECTION("no collider") {
        REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider{5}, 100., 1., 1.) ==
                0);
    }
    SECTION("calcggH") {
        auto smlike =
            HP::BsmParticle("Hsm", HP::ECharge::neutral, HP::CP::undefined);
        smlike.setMass(300);
        effectiveCouplingInput(smlike, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs, false);
        double ggH_XS_calcggH_false =
            smlike.cxn(HP::Collider::LHC13, HP::Production::ggH);
        effectiveCouplingInput(smlike, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs, true);
        double ggH_XS_calcggH_true =
            smlike.cxn(HP::Collider::LHC13, HP::Production::ggH);
        REQUIRE(ggH_XS_calcggH_true == Approx(ggH_XS_calcggH_false));
    }
    SECTION("calcHgamgam") {
        auto smlike =
            HP::BsmParticle("Hsm", HP::ECharge::neutral, HP::CP::undefined);
        smlike.setMass(300);
        effectiveCouplingInput(smlike, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs, false, false);
        double Hgamgam_br_calcHgamgam_false = smlike.br(HP::Decay::gamgam);
        effectiveCouplingInput(smlike, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs, false, true);
        double Hgamgam_br_calcHgamgam_true = smlike.br(HP::Decay::gamgam);
        ;
        REQUIRE(Hgamgam_br_calcHgamgam_true ==
                Approx(Hgamgam_br_calcHgamgam_false));
    }
    SECTION("out of bounds behaviour") {
        SECTION("set to zero for too large masses") {
            REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC8, 1e6,
                                                   1. + 2.i, 2. + 3.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC13, 1e6,
                                                   1. + 2.i, 2. + 3.i) == 0.);
        }
        SECTION("clamped to the closest value for too low masses") {
            REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC8, 1e-6,
                                                   1. + 2.i, 2. + 3.i) ==
                    HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC8, 1.,
                                                   1. + 2.i, 2. + 3.i));
            REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC13, 1e-6,
                                                   1. + 2.i, 2. + 3.i) ==
                    HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC13, 1.,
                                                   1. + 2.i, 2. + 3.i));
        }
        SECTION("vanishes for vanishing couplings") {
            REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC8, 125, 0.,
                                                   0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC13, 125, 0.,
                                                   0.) == 0.);
        }
    }
}

TEST_CASE("WH cxns") {
    SECTION("LHC8") {
        for (const auto &val : testWHLHC8) {
            INFO("mass = " << val.m);
            INFO("ghw = " << val.ghw);
            INFO("ght = " << val.ght);
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC8, val.m,
                                                    val.ghw, val.ght) ==
                    Approx(val.cxn));
        }
    }
    SECTION("LHC13") {
        for (const auto &val : testWHLHC13) {
            INFO("mass = " << val.m);
            INFO("ghw = " << val.ghw);
            INFO("ght = " << val.ght);
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC13, val.m,
                                                    val.ghw, val.ght) ==
                    Approx(val.cxn));
        }
    }
    SECTION("no collider") {
        REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider{5}, 100., 1.,
                                                1.) == 0);
    }
    SECTION("out of bounds behaviour") {
        SECTION("set to zero for too large masses") {
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC8, 1e6, 1.,
                                                    2. + 3.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC13, 1e6,
                                                    1., 2. + 3.i) == 0.);
        }
        SECTION("clamped to the closest value for too low masses") {
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC8, 1e-6,
                                                    1., 2. + 3.i) ==
                    HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC8, 1., 1.,
                                                    2. + 3.i));
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC13, 1e-6,
                                                    1., 2. + 3.i) ==
                    HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC13, 1., 1.,
                                                    2. + 3.i));
        }
        SECTION("vanishes for vanishing couplings") {
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC8, 125, 0.,
                                                    0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ppHW(HP::Collider::LHC13, 125,
                                                    0., 0.) == 0.);
        }
    }
}

TEST_CASE("ZH cxns") {
    SECTION("LHC8") {
        SECTION("gg->ZH") {
            for (const auto &val : testggZHLHC8) {
                INFO("mass = " << val.m);
                INFO("ghz = " << val.ghz);
                INFO("ght = " << val.ght << "+ i " << val.gat);
                INFO("ghb = " << val.ghb << "+ i " << val.gab);
                REQUIRE(HP::EffectiveCouplingCxns::ggHZ(
                            HP::Collider::LHC8, val.m, val.ghz,
                            val.ght + 1i * val.gat,
                            val.ghb + 1i * val.gab) == Approx(val.cxn));
            }
        }
        SECTION("qq->ZH") {
            for (const auto &val : testqqZHLHC8) {
                INFO("mass = " << val.m);
                INFO("ghz = " << val.ghz);
                INFO("ght = " << val.ght);
                REQUIRE(HP::EffectiveCouplingCxns::qqHZ(
                            HP::Collider::LHC8, val.m, val.ghz, val.ght) ==
                        Approx(val.cxn));
            }
        }
        SECTION("bb->ZH") {
            for (const auto &val : testbbZHLHC8) {
                INFO("mass = " << val.m);
                INFO("ghb = " << val.ghb << "+ i " << val.gab);
                REQUIRE(HP::EffectiveCouplingCxns::bbHZ(
                            HP::Collider::LHC8, val.m,
                            val.ghb + 1i * val.gab) == Approx(val.cxn));
            }
        }
        SECTION("pp->ZH") {
            for (const auto &val : testZHLHC8) {
                INFO("mass = " << val.m);
                INFO("ghz = " << val.ghz);
                INFO("ght = " << val.ght << "+ i " << val.gat);
                INFO("ghb = " << val.ghb << "+ i " << val.gab);
                REQUIRE(HP::EffectiveCouplingCxns::ppHZ(
                            HP::Collider::LHC8, val.m, val.ghz,
                            val.ght + 1i * val.gat,
                            val.ghb + 1i * val.gab) == Approx(val.cxn));
            }
        }
    }
    SECTION("LHC13") {
        SECTION("gg->ZH") {
            for (const auto &val : testggZHLHC13) {
                INFO("mass = " << val.m);
                INFO("ghz = " << val.ghz);
                INFO("ght = " << val.ght << "+ i " << val.gat);
                INFO("ghb = " << val.ghb << "+ i " << val.gab);
                REQUIRE(HP::EffectiveCouplingCxns::ggHZ(
                            HP::Collider::LHC13, val.m, val.ghz,
                            val.ght + 1i * val.gat,
                            val.ghb + 1i * val.gab) == Approx(val.cxn));
            }
        }
        SECTION("pure CP-odd special case (due to K-factor)") {
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC13, 200, 0,
                                                    1i, 1i) ==
                    Approx(0.0206313545));
        }
        SECTION("qq->ZH") {
            for (const auto &val : testqqZHLHC13) {
                INFO("mass = " << val.m);
                INFO("ghz = " << val.ghz);
                INFO("ght = " << val.ght);
                REQUIRE(HP::EffectiveCouplingCxns::qqHZ(
                            HP::Collider::LHC13, val.m, val.ghz, val.ght) ==
                        Approx(val.cxn));
            }
        }
        SECTION("bb->ZH") {
            for (const auto &val : testbbZHLHC13) {
                INFO("mass = " << val.m);
                INFO("ghb = " << val.ghb << "+ i " << val.gab);
                REQUIRE(HP::EffectiveCouplingCxns::bbHZ(
                            HP::Collider::LHC13, val.m,
                            val.ghb + 1i * val.gab) == Approx(val.cxn));
            }
        }
        SECTION("pp->ZH") {
            for (const auto &val : testZHLHC13) {
                INFO("mass = " << val.m);
                INFO("ghz = " << val.ghz);
                INFO("ght = " << val.ght << "+ i " << val.gat);
                INFO("ghb = " << val.ghb << "+ i " << val.gab);
                REQUIRE(HP::EffectiveCouplingCxns::ppHZ(
                            HP::Collider::LHC13, val.m, val.ghz,
                            val.ght + 1i * val.gat,
                            val.ghb + 1i * val.gab) == Approx(val.cxn));
            }
        }
    }

    SECTION("no collider") {
        REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider{5}, 100., 1., 1.,
                                                1.) == 0);
        REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider{5}, 100., 1.,
                                                1.) == 0);
        REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider{5}, 100., 1.) ==
                0);
        REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider{5}, 100., 1., 1.,
                                                1.) == 0);
    }

    SECTION("out of bounds behaviour") {
        SECTION("set to zero for too large masses") {
            REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC8, 1e6, 1.,
                                                    2. + 3.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC8, 1e6, 1.,
                                                    2. + 3.i, 4. + 5.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC8, 1e6,
                                                    2. + 3.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC8, 1e6, 1.,
                                                    2. + 3.i, 4. + 5.i) == 0.);

            REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC13, 1e6,
                                                    1., 2. + 3.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC13, 1e6,
                                                    1., 2. + 3.i,
                                                    4. + 5.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC13, 1e6,
                                                    2. + 3.i) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC13, 1e6,
                                                    1., 2. + 3.i,
                                                    4. + 5.i) == 0.);
        }
        SECTION("clamped to the closest value for too low masses") {
            REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC8, 1e-6,
                                                    1., 2. + 3.i) ==
                    HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC8, 1., 1.,
                                                    2. + 3.i));
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC8, 1e-6,
                                                    1., 2. + 3.i, 4. + 5.i) ==
                    HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC8, 1., 1.,
                                                    2. + 3.i, 4. + 5.i));
            REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC8, 1e-6,
                                                    2. + 3.i) ==
                    HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC8, 1.,
                                                    2. + 3.i));
            REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC8, 1e-6,
                                                    1., 2. + 3.i, 4. + 5.i) ==
                    HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC8, 1., 1.,
                                                    2. + 3.i, 4. + 5.i));

            REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC13, 1e-6,
                                                    1., 2. + 3.i) ==
                    HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC13, 1., 1.,
                                                    2. + 3.i));
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC13, 1e-6,
                                                    1., 2. + 3.i, 4. + 5.i) ==
                    HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC13, 1., 1.,
                                                    2. + 3.i, 4. + 5.i));
            REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC13, 1e-6,
                                                    2. + 3.i) ==
                    HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC13, 1.,
                                                    2. + 3.i));
            REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC13, 1e-6,
                                                    1., 2. + 3.i, 4. + 5.i) ==
                    HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC13, 1., 1.,
                                                    2. + 3.i, 4. + 5.i));
        }
        SECTION("vanishes for vanishing couplings") {
            REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC8, 1e6, 0.,
                                                    0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC8, 1e6, 0.,
                                                    0., 0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC8, 1e6,
                                                    0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC8, 1e6, 0.,
                                                    0., 0.) == 0.);

            REQUIRE(HP::EffectiveCouplingCxns::qqHZ(HP::Collider::LHC13, 1e6,
                                                    0., 0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ggHZ(HP::Collider::LHC13, 1e6,
                                                    0., 0., 0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::bbHZ(HP::Collider::LHC13, 1e6,
                                                    0.) == 0.);
            REQUIRE(HP::EffectiveCouplingCxns::ppHZ(HP::Collider::LHC13, 1e6,
                                                    0., 0., 0.) == 0.);
        }
    }
}

TEST_CASE("Hpmtb cross section") {
    using HP::EffectiveCouplingCxns::ppHpmtb;

    SECTION("sanity checks in mass range") {
        auto m = GENERATE(150, 182, 385, 951);
        auto ct = GENERATE(1e-2, 0.2, 0.6, 1, 1.2, 2.3);
        auto cb = GENERATE(1e-2, 0.3, 0.8, 1, 1.9, 5.3);
        // while not physical for most masses, these should still be well
        // behaved
        auto brt = GENERATE(1e-6, 0.1, 0.4);
        REQUIRE(ppHpmtb(HP::Collider::LHC13, m, ct, cb, brt) > 0);
        REQUIRE(ppHpmtb(HP::Collider::LHC13, m, 2 * ct, 2 * cb, brt) ==
                Approx(4 * ppHpmtb(HP::Collider::LHC13, m, ct, cb, brt)));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, m, 0.6 * ct, 0.6 * cb, brt) ==
                Approx(0.36 * ppHpmtb(HP::Collider::LHC13, m, ct, cb, brt)));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, m, 0, 0, brt) == 0);
        REQUIRE(ppHpmtb(HP::Collider::LHC13, m, ct, cb, 1) == 0);
    }

    SECTION("out of bounds") {
        REQUIRE(ppHpmtb(HP::Collider::LHC8, 150, 1, 1, 0) == 0);
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 10, 1, 1, 0) == 0);
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 10, 1, 1, 0) == 0);
    }

    SECTION("comparison with HB5 values") {
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 170.4, 0.8, 1.3, 0.9) ==
                Approx(5.7104019821142930E-002));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 272.64, 0.8, 1.3, 0.9) ==
                Approx(2.2887426804620168E-002));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 374.88, 0.8, 1.3, 0.9) ==
                Approx(1.0219957915073109E-002));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 477.12, 0.8, 1.3, 0.9) ==
                Approx(4.9319906933707625E-003));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 579.36, 0.8, 1.3, 0.9) ==
                Approx(2.5745494441502814E-003));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 681.6, 0.8, 1.3, 0.9) ==
                Approx(1.3950146807622707E-003));
        REQUIRE(ppHpmtb(HP::Collider::LHC13, 783.84, 0.8, 1.3, 0.9) ==
                Approx(7.8918942939299266E-004));
    }
}

TEST_CASE("HpmPhi cross section") {
    using HP::EffectiveCouplingCxns::ppHpmPhi;

    SECTION("sanity checks in mass range") {
        auto mHp = GENERATE(141, 263, 427);
        auto mPhi = GENERATE(26, 182, 322, 491);
        auto c = GENERATE(1e-2, 0.4, 0.9, 1.3, 10);
        INFO(mHp << " " << mPhi << " " << c);
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, mHp, mPhi, c) > 0);
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, mHp, mPhi, 2 * c) ==
                Approx(4 * ppHpmPhi(HP::Collider::LHC13, mHp, mPhi, c)));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, mHp, mPhi, 0.6 * c) ==
                Approx(0.36 * ppHpmPhi(HP::Collider::LHC13, mHp, mPhi, c)));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, mHp, mPhi, 0) == 0);
    }

    SECTION("out of bounds") {
        REQUIRE(ppHpmPhi(HP::Collider::LHC8, 150, 200, 1) == 0);
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 10, 200, 1) == 0);
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 1200, 300, 1) == 0);
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200, 3, 1) == 0);
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200, 1111, 1) == 0);
    }

    SECTION("comparison with python values") {
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 108.0, 13.0, 1) ==
                Approx(2.161983134986043).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200.75, 13.0, 1) ==
                Approx(0.22028794567140134).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 293.5, 13.0, 1) ==
                Approx(0.055308273635297).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 386.25, 13.0, 1) ==
                Approx(0.019640424378485146).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 479.0, 13.0, 1) ==
                Approx(0.00847403491280959).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 108.0, 107.0, 1) ==
                Approx(0.3812701053308828).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200.75, 107.0, 1) ==
                Approx(0.09680372842069261).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 293.5, 107.0, 1) ==
                Approx(0.033219638971781915).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 386.25, 107.0, 1) ==
                Approx(0.013805204937726764).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 479.0, 107.0, 1) ==
                Approx(0.006514286358261084).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 108.0, 201.0, 1) ==
                Approx(0.09635563864675982).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200.75, 201.0, 1) ==
                Approx(0.03837914087702675).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 293.5, 201.0, 1) ==
                Approx(0.0168988663248099).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 386.25, 201.0, 1) ==
                Approx(0.008209497656045834).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 479.0, 201.0, 1) ==
                Approx(0.0042887702064803).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 108.0, 295.0, 1) ==
                Approx(0.032821219712247836).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200.75, 295.0, 1) ==
                Approx(0.016766167251892356).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 293.5, 295.0, 1) ==
                Approx(0.008722870331844141).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 386.25, 295.0, 1) ==
                Approx(0.004725717539573979).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 479.0, 295.0, 1) ==
                Approx(0.0026849162657831063).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 108.0, 389.0, 1) ==
                Approx(0.013592861107744806).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200.75, 389.0, 1) ==
                Approx(0.008019936915620004).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 293.5, 389.0, 1) ==
                Approx(0.004684293466104349).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 386.25, 389.0, 1) ==
                Approx(0.0027698854559177572).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 479.0, 389.0, 1) ==
                Approx(0.0016847909698605812).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 108.0, 483.0, 1) ==
                Approx(0.006349121641106626).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 200.75, 483.0, 1) ==
                Approx(0.004164901133542265).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 293.5, 483.0, 1) ==
                Approx(0.002636781283376125).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 386.25, 483.0, 1) ==
                Approx(0.0016689694799214077).epsilon(1e-2));
        REQUIRE(ppHpmPhi(HP::Collider::LHC13, 479.0, 483.0, 1) ==
                Approx(0.0010614545548502255).epsilon(1e-2));
    }
}

TEST_CASE("qqPhi and qqPhi+photon cxns") {
    auto coll = GENERATE(
        Catch::Generators::from_range(magic_enum::enum_values<HP::Collider>()));

    auto refCoup = std::complex{5., 0.};
    auto equalCoup = GENERATE(std::complex{0., 5.}, std::complex{3., 4.},
                              std::complex{4., 3.});

    double m = GENERATE(1e-5, 20, 250., 600., 1000., 4e5);
    SECTION("neutral scalar") {
        auto cxn = GENERATE(
            HP::EffectiveCouplingCxns::uuHgam,
            HP::EffectiveCouplingCxns::ddHgam,
            HP::EffectiveCouplingCxns::ccHgam,
            HP::EffectiveCouplingCxns::ssHgam,
            HP::EffectiveCouplingCxns::bbHgam,
            HP::EffectiveCouplingCxns::ucHgam,
            HP::EffectiveCouplingCxns::dsHgam,
            HP::EffectiveCouplingCxns::dbHgam,
            HP::EffectiveCouplingCxns::sbHgam, HP::EffectiveCouplingCxns::uuH,
            HP::EffectiveCouplingCxns::ddH, HP::EffectiveCouplingCxns::ccH,
            HP::EffectiveCouplingCxns::ssH, HP::EffectiveCouplingCxns::ucH,
            HP::EffectiveCouplingCxns::dsH, HP::EffectiveCouplingCxns::dbH,
            HP::EffectiveCouplingCxns::sbH);
        if (coll == HP::Collider::LHC13 && m < 1150) {
            CHECK(cxn(coll, m, refCoup) > 0);
        } else {
            CHECK(cxn(coll, m, refCoup) == 0);
        }
        CHECK(cxn(coll, m, 0.) == 0.);
        CHECK(cxn(coll, m, refCoup) == Approx(cxn(coll, m, equalCoup)));
    }
    SECTION("charged scalar") {
        auto cxn = GENERATE(
            HP::EffectiveCouplingCxns::udHpgam,
            HP::EffectiveCouplingCxns::usHpgam,
            HP::EffectiveCouplingCxns::ubHpgam,
            HP::EffectiveCouplingCxns::cdHpgam,
            HP::EffectiveCouplingCxns::csHpgam,
            HP::EffectiveCouplingCxns::cbHpgam,
            HP::EffectiveCouplingCxns::udHmgam,
            HP::EffectiveCouplingCxns::usHmgam,
            HP::EffectiveCouplingCxns::ubHmgam,
            HP::EffectiveCouplingCxns::cdHmgam,
            HP::EffectiveCouplingCxns::csHmgam,
            HP::EffectiveCouplingCxns::cbHmgam, HP::EffectiveCouplingCxns::udHp,
            HP::EffectiveCouplingCxns::usHp, HP::EffectiveCouplingCxns::ubHp,
            HP::EffectiveCouplingCxns::cdHp, HP::EffectiveCouplingCxns::csHp,
            HP::EffectiveCouplingCxns::cbHp, HP::EffectiveCouplingCxns::udHm,
            HP::EffectiveCouplingCxns::usHm, HP::EffectiveCouplingCxns::ubHm,
            HP::EffectiveCouplingCxns::cdHm, HP::EffectiveCouplingCxns::csHm,
            HP::EffectiveCouplingCxns::cbHm);
        if (coll == HP::Collider::LHC13 && m < 1150) {
            CHECK(cxn(coll, m, refCoup.real(), refCoup.imag()) > 0);
        } else {
            CHECK(cxn(coll, m, refCoup.real(), refCoup.imag()) == 0);
        }
        CHECK(cxn(coll, m, refCoup.real(), refCoup.imag()) ==
              Approx(cxn(coll, m, equalCoup.real(), equalCoup.imag())));
        CHECK(cxn(coll, m, 0., 0.) == 0.);
    }
}
