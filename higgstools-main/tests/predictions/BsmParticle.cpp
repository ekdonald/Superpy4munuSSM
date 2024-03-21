#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <cstddef>
#include <magic_enum.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/remove.hpp>
#include <range/v3/view/zip.hpp>
#include <string>
#include <string_view>

using Catch::Approx;

namespace HP = Higgs::predictions;

namespace {

const auto allHadrProds =
    magic_enum::enum_values<HP::Production>() |
    ranges::views::filter(
        [](auto p) { return HP::validProductionAt(p, HP::ColliderType::pp); }) |
    ranges::to<std::vector>;

const auto neutralHadrProds =
    allHadrProds | ranges::views::filter([](auto p) {
        return HP::validProductionFor(p, HP::ECharge::neutral);
    }) |
    ranges::to<std::vector>;

const auto chargedHadrProds =
    allHadrProds | ranges::views::filter([](auto p) {
        return HP::validProductionFor(p, HP::ECharge::single);
    }) |
    ranges::to<std::vector>;

const auto fixChargeHadrProds = std::array{neutralHadrProds, chargedHadrProds};

const auto settableNeutralHadrProds =
    neutralHadrProds | ranges::views::remove(HP::Production::H) |
    ranges::views::remove(HP::Production::HZ) |
    ranges::views::remove(HP::Production::Ht) | ranges::to<std::vector>;

const auto &settableChargedProds = chargedHadrProds;

const auto settableHadrProds =
    std::array{settableNeutralHadrProds, settableChargedProds};

constexpr auto allDecays = magic_enum::enum_values<HP::Decay>();

const auto neutralDecays =
    allDecays | ranges::views::filter([](auto d) {
        return HP::validDecayFor(d, HP::ECharge::neutral);
    }) |
    ranges::to<std::vector>;

const auto chargedDecays = allDecays | ranges::views::filter([](auto d) {
                               return HP::validDecayFor(d, HP::ECharge::single);
                           }) |
                           ranges::to<std::vector>;

const auto fixChargeDecays = std::array{neutralDecays, chargedDecays};

const auto settableNeutralDecays = neutralDecays |
                                   ranges::views::remove(HP::Decay::inv) |
                                   ranges::to<std::vector>;

const auto &settableChargedDecays = chargedDecays;

const auto settableDecays =
    std::array{settableNeutralDecays, settableChargedDecays};
constexpr auto colliders = magic_enum::enum_values<HP::Collider>();
} // namespace

using namespace magic_enum::ostream_operators;

TEST_CASE("BSMParticle properties") {
    for (auto charge : magic_enum::enum_values<HP::ECharge>()) {
        for (auto cp : magic_enum::enum_values<HP::CP>()) {
            auto p = HP::BsmParticle("h", charge, cp);
            p.setMass(200);
            p.setMassUnc(1.);
            REQUIRE(p.id() == "h");
            REQUIRE(p.mass() == Approx(200.));
            REQUIRE(p.cp() == cp);
            REQUIRE(p.charge() == charge);
            REQUIRE(p.massUnc() == Approx(1.));
        }
    }

    SECTION("copying") {
        auto p = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        p.setDecayWidth(HP::Decay::bb, 1.);
        auto p1 = p;
        REQUIRE(p1.id() == p.id());
        REQUIRE(p1.br(HP::Decay::bb) == p.br(HP::Decay::bb));

        HP::Particle *pp = &p;
        REQUIRE(pp->br(HP::Decay::bb) == p.br(HP::Decay::bb));
        auto pp1 = pp->clone();
        REQUIRE(pp1->id() == p.id());
        REQUIRE(pp1->br(HP::Decay::bb) == p.br(HP::Decay::bb));
    }
}

TEST_CASE("BsmParticle rates") {
    auto h0 = HP::BsmParticle("h0", HP::ECharge::neutral, HP::CP::even);
    auto hplus = HP::BsmParticle("h+", HP::ECharge::single, HP::CP::undefined);
    auto testParticles = std::array{h0, hplus};
    for (const auto &[s, cDecs, setDecs, cProds, setProds] :
         ranges::views::zip(testParticles, fixChargeDecays, settableDecays,
                            fixChargeHadrProds, settableHadrProds)) {
        DYNAMIC_SECTION("Default zeros " << s.id()) {
            REQUIRE(s.totalWidth() == Approx(0.));
            for (auto p : allHadrProds) {
                for (auto d : allDecays) {
                    REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) == 0.);
                    REQUIRE(s.cxn(HP::Collider::LHC8, p) == 0.);
                    REQUIRE(s.br(d) == 0.);
                }
            }
        }

        DYNAMIC_SECTION("all of these should be no-ops " << s.id()) {
            s.setCxn(HP::Collider::LHC13, HP::Production::none, 100);
            REQUIRE_NOTHROW(s.setBr(HP::Decay::none, 200));
            REQUIRE_NOTHROW(s.setDecayWidth(HP::Decay::none, -100));
            for (auto p : cProds) {
                s.setChannelRate(HP::Collider::LHC13, p, HP::Decay::none, 200.);
            }
            for (auto d : cDecs) {
                s.setChannelRate(HP::Collider::LHC13, HP::Production::none, d,
                                 300.);
            }

            for (auto p : allHadrProds) {
                for (auto d : allDecays) {
                    INFO(magic_enum::enum_name(p)
                         << "." << magic_enum::enum_name(d));
                    REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) == 0.);
                }
            }
        }

        DYNAMIC_SECTION("total width " << s.id()) {
            s.setTotalWidth(0.8);
            REQUIRE(s.totalWidth() == Approx(0.8));
            REQUIRE_THROWS_AS(s.setTotalWidth(-1.), HP::InvalidInput);
        }
        s.setTotalWidth(1.0);

        DYNAMIC_SECTION("set branching ratios only " << s.id()) {
            for (const auto &[j, d] : ranges::views::enumerate(setDecs)) {
                s.setBr(d, j / 500.);
            }
            for (auto d : cDecs) {
                REQUIRE(s.channelRate(HP::Collider::LHC13, HP::Production::none,
                                      d) == s.br(d));
            }
            for (const auto &[j, d] : ranges::views::enumerate(setDecs)) {
                REQUIRE(s.br(d) == Approx(j / 500.));
            }

            DYNAMIC_SECTION("combined decay modes " << s.id()) {
                CHECK(s.br(HP::Decay::inv) ==
                      Approx(s.br(HP::Decay::directInv) +
                             s.br(HP::Decay::ZZ) *
                                 std::pow(HP::constants::b_Z_inv, 2)));
                if (s.charge() == HP::ECharge::neutral) {
                    s.setBr(HP::Decay::inv, 0.);
                }
                CHECK(s.br(HP::Decay::directInv) == 0.);
            }

            for (auto p :
                 allHadrProds | ranges::views::remove(HP::Production::none)) {
                for (auto d : allDecays) {
                    INFO(magic_enum::enum_name(p)
                         << "." << magic_enum::enum_name(p));
                    if (HP::validProductionAt(p, HP::ColliderType::pp)) {
                        REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) == 0);
                    }
                }
            }
        }

        DYNAMIC_SECTION("set cxn only " << s.id()) {
            for (auto [j, p] : ranges::views::enumerate(setProds)) {
                s.setCxn(HP::Collider::LHC13, p, j);
            }
            for (auto p : allHadrProds) {
                REQUIRE(
                    s.channelRate(HP::Collider::LHC13, p, HP::Decay::none) ==
                    s.cxn(HP::Collider::LHC13, p));
            }
            for (auto [j, p] : ranges::views::enumerate(setProds)) {
                REQUIRE(s.cxn(HP::Collider::LHC13, p) == j);
            }
            DYNAMIC_SECTION("combined production modes " << s.id()) {
                CHECK(s.cxn(HP::Collider::LHC13, HP::Production::H) ==
                      Approx(s.cxn(HP::Collider::LHC13, HP::Production::ggH) +
                             s.cxn(HP::Collider::LHC13, HP::Production::bbH) +
                             s.cxn(HP::Collider::LHC13, HP::Production::qqH)));
                CHECK(s.cxn(HP::Collider::LHC13, HP::Production::HZ) ==
                      Approx(s.cxn(HP::Collider::LHC13, HP::Production::ggHZ) +
                             s.cxn(HP::Collider::LHC13, HP::Production::qqHZ) +
                             s.cxn(HP::Collider::LHC13, HP::Production::bbHZ)));
                CHECK(s.cxn(HP::Collider::LHC13, HP::Production::Ht) ==
                      Approx(
                          s.cxn(HP::Collider::LHC13, HP::Production::schanHt) +
                          s.cxn(HP::Collider::LHC13, HP::Production::tchanHt)));
            }
            for (auto p : allHadrProds) {
                for (auto d :
                     allDecays | ranges::views::remove(HP::Decay::none)) {
                    INFO(magic_enum::enum_name(p)
                         << "." << magic_enum::enum_name(d));
                    if (HP::validProductionAt(p, HP::ColliderType::pp)) {
                        REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) == 0);
                    }
                }
            }
        }
        DYNAMIC_SECTION("set both production and decay channels " << s.id()) {
            for (const auto &[j, d] : ranges::views::enumerate(setDecs)) {
                s.setBr(d, j / 500.);
            }
            for (auto [j, p] : ranges::views::enumerate(setProds)) {
                s.setCxn(HP::Collider::LHC13, p, j);
            }
            for (auto p : allHadrProds) {
                for (auto d : allDecays) {
                    INFO(magic_enum::enum_name(p)
                         << " " << magic_enum::enum_name(d) << " = "
                         << s.cxn(HP::Collider::LHC13, p) << " * " << s.br(d));
                    if (HP::validProductionAt(p, HP::ColliderType::pp)) {
                        if (p == HP::Production::none) {
                            REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                    s.br(d));
                        } else if (d == HP::Decay::none) {
                            REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                    s.cxn(HP::Collider::LHC13, p));
                        } else {
                            REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                    Approx(s.cxn(HP::Collider::LHC13, p) *
                                           s.br(d)));
                        }
                    }
                }
            }
            for (auto [i, p] : ranges::views::enumerate(
                     setProds | ranges::views::remove(HP::Production::none))) {
                for (const auto &[j, d] : ranges::views::enumerate(
                         setDecs | ranges::views::remove(HP::Decay::none))) {
                    REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                            Approx((i + 1) * (j + 1) / 500.));
                }
            }
        }

        DYNAMIC_SECTION("channelRates with nothing set beneath " << s.id()) {
            for (auto [i, p] : ranges::views::enumerate(cProds)) {
                for (const auto &[j, d] : ranges::views::enumerate(cDecs)) {
                    s.setChannelRate(HP::Collider::LHC13, p, d, i * j);
                }
            }
            for (auto [i, p] : ranges::views::enumerate(cProds)) {
                for (const auto &[j, d] : ranges::views::enumerate(cDecs)) {
                    REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) == i * j);
                }
            }
            s.resetChannelRates();
            for (auto p : allHadrProds) {
                for (auto d :
                     allDecays | ranges::views::remove(HP::Decay::none)) {
                    INFO(magic_enum::enum_name(p)
                         << "." << magic_enum::enum_name(d));
                    if (HP::validProductionAt(p, HP::ColliderType::pp)) {
                        REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) == 0);
                    }
                }
            }
        }

        DYNAMIC_SECTION("channelRates with values set beneath " << s.id()) {
            for (const auto &[j, d] : ranges::views::enumerate(setDecs)) {
                s.setBr(d, j / 500.);
            }
            for (auto [j, p] : ranges::views::enumerate(setProds)) {
                s.setCxn(HP::Collider::LHC13, p, j);
            }
            for (auto [i, p] : ranges::views::enumerate(cProds)) {
                for (const auto &[j, d] : ranges::views::enumerate(cDecs)) {
                    INFO(magic_enum::enum_name(p)
                         << "." << magic_enum::enum_name(d));
                    s.setChannelRate(HP::Collider::LHC13, p, d, i * j);
                    if (p == HP::Production::none) {
                        REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                s.br(d));
                    } else if (d == HP::Decay::none) {
                        REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                s.cxn(HP::Collider::LHC13, p));
                    } else {
                        REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                i * j);
                    }
                }
            }
            s.resetChannelRates();
            for (auto p : allHadrProds) {
                for (auto d : allDecays) {
                    INFO(magic_enum::enum_name(p)
                         << " " << magic_enum::enum_name(d) << " = "
                         << s.cxn(HP::Collider::LHC13, p) << " * " << s.br(d));
                    if (HP::validProductionAt(p, HP::ColliderType::pp)) {
                        if (p == HP::Production::none) {
                            REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                    s.br(d));
                        } else if (d == HP::Decay::none) {
                            REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                    s.cxn(HP::Collider::LHC13, p));
                        } else {
                            REQUIRE(s.channelRate(HP::Collider::LHC13, p, d) ==
                                    Approx(s.cxn(HP::Collider::LHC13, p) *
                                           s.br(d)));
                        }
                    }
                }
            }
        }
    }

    SECTION("error cases") {
        REQUIRE_THROWS_AS(
            h0.setCxn(HP::Collider::LHC13, HP::Production::vbfHpm, 1.),
            HP::InvalidChannel);
        REQUIRE_THROWS_AS(h0.setCxn(HP::Collider::LEP, HP::Production::Htt, 1.),
                          HP::InvalidChannel);

        SECTION("combined production modes") {
            REQUIRE_THROWS_AS(
                h0.setCxn(HP::Collider::LHC13, HP::Production::H, 1.),
                HP::InvalidInput);
            REQUIRE_THROWS_AS(
                h0.setCxn(HP::Collider::LHC13, HP::Production::HZ, 1.),
                HP::InvalidInput);
            REQUIRE_THROWS_AS(
                h0.setCxn(HP::Collider::LHC13, HP::Production::Ht, 1.),
                HP::InvalidInput);
        }
    }
}

TEST_CASE("Setting BSMParticle normalized cxns") {
    SECTION("neutral scalar -- SMHiggs") {
        auto ref = GENERATE(HP::ReferenceModel::SMHiggs,
                            HP::ReferenceModel::SMHiggsEW,
                            HP::ReferenceModel::SMHiggsInterp);
        auto s = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        s.setMass(100);
        // set Cxn and overwrite with zero norm cxn
        for (auto [j, p] : ranges::views::enumerate(settableNeutralHadrProds)) {
            s.setCxn(HP::Collider::LHC13, p, j);
            s.setNormalizedCxn(HP::Collider::LHC13, p, 0., ref);
        }
        for (auto p : allHadrProds) {
            REQUIRE(s.cxn(HP::Collider::LHC13, p) == 0);
        }
        for (auto p : settableNeutralHadrProds) {
            s.setNormalizedCxn(HP::Collider::LHC13, p, 1., ref);
        }
        auto refP = HP::getReference(ref, s.mass());
        for (auto p : allHadrProds) {
            INFO(magic_enum::enum_name(p));
            REQUIRE(s.cxn(HP::Collider::LHC13, p) ==
                    Approx(refP->cxn(HP::Collider::LHC13, p)));
        }
    }
}

TEST_CASE("Setting BSMParticle BRs") {
    auto p = HP::BsmParticle("Phi", HP::ECharge::neutral, HP::CP::undefined);
    p.setTotalWidth(10.);
    constexpr auto bb = 0.7;
    constexpr auto ZH1 = 0.2;
    constexpr auto H1H2 = 0.1;
    REQUIRE_NOTHROW(p.setBr(HP::Decay::bb, 0.7));
    REQUIRE_NOTHROW(p.setBr(HP::ChainDecay::Z, "H1", 0.2));
    REQUIRE_NOTHROW(p.setBr("H1", "H2", 0.1));

    CHECK(p.totalWidth() == Approx(10.));
    CHECK(p.br(HP::Decay::bb) == Approx(0.7));
    CHECK(p.br(HP::ChainDecay::Z, "H1") == Approx(0.2));
    CHECK(p.br("H1", "H2") == Approx(0.1));
    CHECK(p.br("H2", "H1") == Approx(0.1));

    SECTION("errors when sum(BRs)>1") {
        REQUIRE_THROWS_AS(p.setBr(HP::Decay::gg, 0.2), HP::InvalidInput);
        REQUIRE_THROWS_AS(p.setBr(HP::ChainDecay::W, "H+", 0.1),
                          HP::InvalidInput);
        REQUIRE_THROWS_AS(p.setBr("H1", "H1", 0.03), HP::InvalidInput);

        auto s =
            HP::BsmParticle("Phi", HP::ECharge::neutral, HP::CP::undefined);
        s.setTotalWidth(2.);
        REQUIRE_THROWS_AS(s.setBr(HP::Decay::ZZ, 2.), HP::InvalidInput);
    }
    SECTION("reducing one BR allows adding another one") {
        REQUIRE_NOTHROW(p.setBr(HP::Decay::bb, 0.5));
        REQUIRE_NOTHROW(p.setBr(HP::Decay::gg, 0.2));
        REQUIRE_NOTHROW(p.setBr(HP::ChainDecay::Z, "H1", 0.1));
        REQUIRE_NOTHROW(p.setBr(HP::ChainDecay::W, "H+", 0.1));
        REQUIRE_NOTHROW(p.setBr("H1", "H2", 0.07));
        REQUIRE_NOTHROW(p.setBr("H1", "H1", 0.03));
        CHECK(p.totalWidth() == Approx(10.));
        CHECK(p.br(HP::Decay::bb) == Approx(0.5));
        CHECK(p.br(HP::Decay::gg) == Approx(0.2));
        CHECK(p.br(HP::ChainDecay::Z, "H1") == Approx(0.1));
        CHECK(p.br(HP::ChainDecay::W, "H+") == Approx(0.1));
        CHECK(p.br("H1", "H2") == Approx(0.07));
        CHECK(p.br("H2", "H1") == Approx(0.07));
        CHECK(p.br("H1", "H1") == Approx(0.03));
    }
    SECTION("changing the total width leaves the BRs unchanges") {
        p.setTotalWidth(1.);
        CHECK(p.br(HP::Decay::bb) == Approx(0.7));
        CHECK(p.br(HP::ChainDecay::Z, "H1") == Approx(0.2));
        CHECK(p.br("H1", "H2") == Approx(0.1));
        CHECK(p.br("H2", "H1") == Approx(0.1));
        p.setTotalWidth(4.14214);
        CHECK(p.br(HP::Decay::bb) == Approx(0.7));
        CHECK(p.br(HP::ChainDecay::Z, "H1") == Approx(0.2));
        CHECK(p.br("H1", "H2") == Approx(0.1));
        CHECK(p.br("H2", "H1") == Approx(0.1));
    }
    SECTION("when the width is zero, all BRs are zero and none can be added") {
        p.setTotalWidth(0.);
        CHECK(p.br(HP::Decay::bb) == 0.);
        CHECK(p.br(HP::ChainDecay::Z, "H1") == 0.);
        CHECK(p.br("H1", "H2") == 0.);
        CHECK(p.br("H2", "H1") == 0.);
        REQUIRE_THROWS_AS(p.setBr(HP::Decay::gg, 0.2), HP::InvalidInput);
        REQUIRE_THROWS_AS(p.setBr(HP::ChainDecay::W, "H+", 0.1),
                          HP::InvalidInput);
        REQUIRE_THROWS_AS(p.setBr("H1", "H1", 0.03), HP::InvalidInput);
    }

    SECTION("error cases") {
        REQUIRE_THROWS_AS(p.setBr(HP::Decay::taunu, 0.), HP::InvalidChannel);
    }
}

TEST_CASE("Setting BSMParticle decay widths") {
    auto p = HP::BsmParticle("Phi", HP::ECharge::neutral, HP::CP::undefined);

    SECTION("invalid widths") {
        REQUIRE_THROWS_AS(p.setDecayWidth(HP::Decay::bb, -0.5),
                          HP::InvalidInput);
        REQUIRE_THROWS_AS(p.setDecayWidth(HP::ChainDecay::Z, "H1", -100.),
                          HP::InvalidInput);
        REQUIRE_THROWS_AS(p.setDecayWidth("X", "H1", -100.), HP::InvalidInput);
    }

    SECTION("tiny negative widths") {
        REQUIRE_NOTHROW(p.setDecayWidth(HP::Decay::bb, -1.0e-11));
        CHECK(p.br(HP::Decay::bb) == Approx(0.0));
    }

    SECTION("adding widths") {
        REQUIRE_NOTHROW(p.setDecayWidth(HP::Decay::bb, 2.0));
        CHECK(p.br(HP::Decay::bb) == Approx(1.0));
        REQUIRE_NOTHROW(p.setDecayWidth(HP::ChainDecay::Z, "H2", 0.5));
        CHECK(p.br(HP::Decay::bb) == Approx(0.8));
        CHECK(p.br(HP::ChainDecay::Z, "H2") == Approx(0.2));
        CHECK(p.totalWidth() == Approx(2.5));
        REQUIRE_NOTHROW(p.setDecayWidth("H1", "H1", 1.5));
        CHECK(p.br(HP::Decay::bb) == Approx(0.5));
        CHECK(p.br(HP::ChainDecay::Z, "H2") == Approx(1 / 8.));
        CHECK(p.br("H1", "H1") == Approx(3 / 8.));
        CHECK(p.totalWidth() == Approx(4.));
    }

    SECTION("widths with prexisting BRs") {
        p.setTotalWidth(3.0);
        p.setBr(HP::Decay::tautau, 2 / 3.);
        p.setBr(HP::ChainDecay::W, "H+", 2 / 9.);
        p.setBr("X+", "H+", 1 / 9.);
        REQUIRE_NOTHROW(p.setDecayWidth(HP::Decay::bb, 2.0));
        CHECK(p.totalWidth() == Approx(5.0));
        CHECK(p.br(HP::Decay::bb) == Approx(2 / 5.));
        CHECK(p.br(HP::Decay::tautau) == Approx(2 / 5.));
        CHECK(p.br(HP::ChainDecay::W, "H+") == Approx(2 / 15.));
        CHECK(p.br("H+", "X+") == Approx(1 / 15.));
        REQUIRE_NOTHROW(p.setDecayWidth("H1", "H2", 1.0));
        CHECK(p.totalWidth() == Approx(6.0));
        CHECK(p.br(HP::Decay::bb) == Approx(1 / 3.));
        CHECK(p.br(HP::Decay::tautau) == Approx(1 / 3.));
        CHECK(p.br(HP::ChainDecay::W, "H+") == Approx(1 / 9.));
        CHECK(p.br("H+", "X+") == Approx(1 / 18.));
        CHECK(p.br("H1", "H2") == Approx(1 / 6.));
    }

    SECTION("reset totalwidth to 0") {
        p.setDecayWidth(HP::Decay::bb, 2.0);
        p.setDecayWidth(HP::ChainDecay::Z, "H2", 0.5);
        p.setDecayWidth("H1", "H1", 1.5);
        CHECK(p.br(HP::Decay::bb) == Approx(0.5));
        CHECK(p.br(HP::ChainDecay::Z, "H2") == Approx(1 / 8.));
        CHECK(p.br("H1", "H1") == Approx(3 / 8.));
        p.setTotalWidth(0.);
        CHECK(p.br(HP::Decay::bb) == 0.);
        CHECK(p.br(HP::ChainDecay::Z, "H2") == 0.);
        CHECK(p.br("H1", "H1") == 0.);
        CHECK(p.totalWidth() == 0.);
    }

    SECTION("error cases") {
        REQUIRE_THROWS_AS(p.setDecayWidth(HP::Decay::taunu, 1.),
                          HP::InvalidChannel);
    }
}

TEST_CASE("Invisible BRs") {
    using HP::BsmParticle, HP::ChainDecay, HP::Decay, HP::constants::b_Z_inv;
    auto pred = Higgs::Predictions{};
    auto &h1 = pred.addParticle(
        BsmParticle{"h1", HP::ECharge::neutral, HP::CP::undefined});
    auto &h2 = pred.addParticle(
        BsmParticle{"h2", HP::ECharge::neutral, HP::CP::undefined});
    auto &h3 = pred.addParticle(
        BsmParticle{"h3", HP::ECharge::neutral, HP::CP::undefined});

    h1.setMass(125);
    h2.setMass(25);
    h3.setMass(340);

    h2.setTotalWidth(1.);
    h3.setTotalWidth(1.);

    h3.setBr("h1", "h2", 0.1);
    h3.setBr("h2", "h2", 0.15);
    h3.setBr("h1", "h1", 0.2);
    h1.setDecayWidth("h2", "h2", 0.25);

    h3.setBr(HP::ChainDecay::Z, "h1", 0.05);
    h3.setBr(HP::ChainDecay::Z, "h2", 0.1);
    h1.setDecayWidth(HP::ChainDecay::Z, "h2", 0.15);

    constexpr auto h1Br = 0.4;
    h1.setDecayWidth(HP::Decay::inv, h1Br);
    h1.setDecayWidth("other", "stuff", 0.2);
    constexpr auto h2Br = 0.5;
    h2.setBr(HP::Decay::inv, h2Br);
    constexpr auto h3Br = 0.3;
    h3.setBr(HP::Decay::inv, h3Br);

    CHECK(h1.br(HP::Decay::inv) ==
          Approx(h1Br + h1.br("h2", "h2") * std::pow(h2Br, 2) +
                 h1.br(HP::ChainDecay::Z, "h2") * b_Z_inv * h2Br));
    CHECK(h2.br(HP::Decay::inv) == h2Br);
    CHECK(h3.br(HP::Decay::inv) ==
          Approx(h3Br + h3.br("h1", "h1") * std::pow(h1.br(HP::Decay::inv), 2) +
                 h3.br("h2", "h2") * std::pow(h2Br, 2) +
                 h3.br("h1", "h2") * h1.br(HP::Decay::inv) * h2Br +
                 h3.br(HP::ChainDecay::Z, "h1") * h1.br(HP::Decay::inv) *
                     b_Z_inv +
                 h3.br(HP::ChainDecay::Z, "h2") * h2Br * b_Z_inv));
}

namespace {
const auto allLepProds =
    magic_enum::enum_values<HP::Production>() |
    ranges::views::filter(
        [](auto p) { return HP::validProductionAt(p, HP::ColliderType::ee); }) |
    ranges::to<std::vector>;

const auto neutralLepProds =
    allLepProds | ranges::views::filter([](auto p) {
        return HP::validProductionFor(p, HP::ECharge::neutral);
    }) |
    ranges::to<std::vector>;
} // namespace

TEST_CASE("LEP cxns") {
    auto s = HP::BsmParticle("Phi", HP::ECharge::neutral, HP::CP::undefined);
    s.setDecayWidth(HP::Decay::gamgam, 1);
    SECTION("Can't set absolute LEP cxns") {
        for (auto p : neutralLepProds) {
            REQUIRE_THROWS_AS(s.setCxn(HP::Collider::LEP, p, 1.),
                              HP::InvalidInput);
        }
    }
    SECTION("Set and retrieve normalized Cxns and channelRates") {
        auto ref = GENERATE(HP::ReferenceModel::SMHiggs,
                            HP::ReferenceModel::SMHiggsEW,
                            HP::ReferenceModel::SMHiggsInterp);
        for (auto [i, p] : neutralLepProds |
                               ranges::views::remove(HP::Production::none) |
                               ranges::views::enumerate) {
            REQUIRE_NOTHROW(s.setNormalizedCxn(HP::Collider::LEP, p, i, ref));
            CHECK(Approx(i) == s.cxn(HP::Collider::LEP, p));
            INFO(s.cxn(HP::Collider::LEP, p) * s.br(HP::Decay::gamgam));
            CHECK(Approx(i) ==
                  s.channelRate(HP::Collider::LEP, p, HP::Decay::gamgam));
            REQUIRE_NOTHROW(s.setChannelRate(HP::Collider::LEP, p,
                                             HP::Decay::gamgam, 10 * i));
            CHECK(Approx(10 * i) ==
                  s.channelRate(HP::Collider::LEP, p, HP::Decay::gamgam));
            REQUIRE(Approx(i) == s.cxn(HP::Collider::LEP, p));
            s.resetChannelRates();
            CHECK(Approx(i) ==
                  s.channelRate(HP::Collider::LEP, p, HP::Decay::gamgam));
        }
    }
}

TEST_CASE("particle context") {

    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
    auto &h2 = pred.addParticle(
        HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));

    h.setMass(100.);
    h2.setMass(10.);
    h.setDecayWidth(HP::Decay::directInv, .5);
    h.setDecayWidth("h2", "h2", .5);
    h2.setDecayWidth(HP::Decay::directInv, 1.);

    REQUIRE(h.br(HP::Decay::directInv) == Approx(0.5));
    REQUIRE(h.br(HP::Decay::inv) == Approx(1.));

    SECTION("copy/move predictions") {
        auto pred2 = pred; // copy construction
        REQUIRE(pred2.particle("h").br(HP::Decay::inv) == Approx(1.));
        pred2 = pred2; // self-assignment
        REQUIRE(pred2.particle("h").br(HP::Decay::inv) == Approx(1.));
        auto pred3 = std::move(pred2); // move construction
        REQUIRE(pred3.particle("h").br(HP::Decay::inv) == Approx(1.));
        pred2 = pred3; // assignment
        REQUIRE(pred2.particle("h").br(HP::Decay::inv) == Approx(1.));
        pred3 = std::move(pred2); // move assignment
        REQUIRE(pred3.particle("h").br(HP::Decay::inv) == Approx(1.));
        pred3 = std::move(pred3); // move self-assignment
        REQUIRE(pred3.particle("h").br(HP::Decay::inv) == Approx(1.));
    }

    SECTION("self assignment works") {
        REQUIRE(h.br(HP::Decay::inv) == Approx(1.));
        h = h;
        CHECK(h.br(HP::Decay::inv) == Approx(1.));
        h = std::move(h);
        CHECK(h.br(HP::Decay::inv) == Approx(1.));
    }

    SECTION("copy/move particles out of context") {
        auto hc = h; // copy construction
        CHECK(hc.br(HP::Decay::inv) == h.br(HP::Decay::directInv));

        hc = h; // assignment
        CHECK(hc.br(HP::Decay::inv) == h.br(HP::Decay::directInv));

        auto pred2 = pred; // make a copy to move from

        auto hc2 = std::move(pred2.particle("h")); // move construction
        CHECK(hc2.br(HP::Decay::inv) == h.br(HP::Decay::directInv));

        pred2 = pred; // make a copy to move from

        hc2 = std::move(pred2.particle("h")); // move assignment
        CHECK(hc2.br(HP::Decay::inv) == h.br(HP::Decay::directInv));
    }

    SECTION("remove and re-add a particle") {
        auto hc = h;
        CHECK(hc.br(HP::Decay::inv) == h.br(HP::Decay::directInv));

        pred.removeParticle("h");
        pred.addParticle(hc);
        REQUIRE(pred.particle("h").br(HP::Decay::inv) == Approx(1.));

        pred.removeParticle("h");
        pred.addParticle(std::move(hc));
        REQUIRE(pred.particle("h").br(HP::Decay::inv) == Approx(1.));
    }
}

TEST_CASE("BSMParticle couplings") {
    auto h = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);

    for (auto [i, c] :
         ranges::views::enumerate(magic_enum::enum_values<HP::Coupling>())) {
        h.setCoupling(c, i);
    }
    for (auto [i, c] :
         ranges::views::enumerate(magic_enum::enum_values<HP::Coupling>())) {
        INFO(magic_enum::enum_name(c));
        CHECK(h.coupling(c) == Approx(i));
    }
    for (auto p : allHadrProds) {
        CHECK(h.cxn(HP::Collider::LHC13, p) == 0.);
    }
    for (auto d : allDecays) {
        CHECK(h.br(d) == Approx(0.));
    }
}
