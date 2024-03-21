#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <magic_enum.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/functional.hpp>
#include <stdexcept>
#include <vector>

namespace HP = Higgs::predictions;
using Catch::Approx;

TEST_CASE("particle IDs") {
    // masses used for validation only
    auto h1 = HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined);
    auto h2 = HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined);
    auto a1 = HP::BsmParticle("A1", HP::ECharge::neutral, HP::CP::undefined);
    auto h17 = HP::BsmParticle("h17", HP::ECharge::neutral, HP::CP::undefined);
    auto hweird =
        HP::BsmParticle("Hey, I'm a Higgs boson!", HP::ECharge::neutral, HP::CP::undefined);

    const auto neutH = std::vector{h1, h2, a1, h17, hweird};

    const auto h1p = HP::BsmParticle("H1+", HP::ECharge::single, HP::CP::undefined);
    const auto superH =
        HP::BsmParticle("Supercharged Higgs", HP::ECharge::single, HP::CP::undefined);
    const auto charH = std::vector{h1p, superH};

    auto pred = Higgs::Predictions{};
    for (const auto &h : neutH) {
        INFO(h.id());
        pred.addParticle(h);
    }
    for (const auto &h : charH) {
        INFO(h.id());
        pred.addParticle(h);
    }
    SECTION("valid constructions and access") {

        auto ids = pred.particleIds();
        REQUIRE(ids.size() == neutH.size() + charH.size());
        for (const auto &h : neutH) {
            INFO(h.id());
            REQUIRE(ranges::find(ids, h.id()) != ids.end());
        }
        for (const auto &h : charH) {
            INFO(h.id());
            REQUIRE(ranges::find(ids, h.id()) != ids.end());
        }
        // remove and re-add
        REQUIRE_NOTHROW(pred.removeParticle(a1.id()));
        REQUIRE_NOTHROW(pred.addParticle(a1));
        REQUIRE_NOTHROW(pred.particle(a1.id()));

        REQUIRE_NOTHROW(pred.removeParticle(h1p.id()));
        REQUIRE_NOTHROW(pred.addParticle(h1p));
        REQUIRE_NOTHROW(pred.particle(h1p.id()));
    }

    SECTION("duplicates") {
        for (const auto &h : neutH) {
            INFO(h.id());
            REQUIRE_THROWS_AS(pred.addParticle(h), HP::InvalidInput);
        }
        for (const auto &h : charH) {
            INFO(h.id());
            REQUIRE_THROWS_AS(pred.addParticle(h), HP::InvalidInput);
        }
        REQUIRE_THROWS_AS(
            pred.addParticle(HP::BsmParticle{h1.id(), HP::ECharge::single, HP::CP::undefined}),
            HP::InvalidInput);
        REQUIRE_THROWS_AS(
            pred.addParticle(HP::BsmParticle{h1p.id(), HP::ECharge::neutral, HP::CP::undefined}),
            HP::InvalidInput);
    }

    SECTION("missing entries") {
        REQUIRE_THROWS_AS(pred.particle("H1"), std::out_of_range);
        REQUIRE_THROWS_AS(pred.particle("h2pp"), std::out_of_range);
        REQUIRE_NOTHROW(pred.removeParticle("H1"));
        REQUIRE_NOTHROW(pred.removeParticle("h2pp"));
    }

    SECTION("vector access") {
        REQUIRE(pred.particles().size() == neutH.size() + charH.size());
        for (const auto &h : neutH) {
            REQUIRE(ranges::find_if(pred.particles(), [&h](const auto &hh) {
                        return h.id() == hh.id();
                    }) != pred.particles().end());
        }
        for (const auto &h : charH) {
            REQUIRE(ranges::find_if(pred.particles(), [&h](const auto &hh) {
                        return h.id() == hh.id();
                    }) != pred.particles().end());
        }
        pred.removeParticle("h1");
        REQUIRE(pred.particles().size() == neutH.size() + charH.size() - 1);
        for (const auto &h : charH) {
            REQUIRE(ranges::find_if(pred.particles(), [&h](const auto &hh) {
                        return h.id() == hh.id();
                    }) != pred.particles().end());
        }
        pred.removeParticle("H1+");
        REQUIRE(pred.particles().size() == neutH.size() + charH.size() - 2);
    }

    SECTION("references") {
        auto pred2 = Higgs::Predictions{};
        auto &h1 =
            pred2.addParticle(HP::BsmParticle{"H1", HP::ECharge::neutral, HP::CP::undefined});
        auto &h2 =
            pred2.addParticle(HP::BsmParticle{"H2", HP::ECharge::neutral, HP::CP::undefined});
        CHECK(pred2.particle("H1").id() == "H1");
        CHECK(pred2.particle("H2").id() == "H2");
        INFO(&h1 << " vs " << &pred2.particle("H1"));

        CHECK(h1.id() == "H1");
        CHECK(h2.id() == "H2");
    }
}

TEST_CASE("pair production") {
    auto pred = Higgs::Predictions();
    auto &h1 = pred.addParticle(HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
    h1.setMass(1);
    auto &h2 = pred.addParticle(HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
    h2.setMass(2);
    auto &h1p = pred.addParticle(HP::BsmParticle("h1p", HP::ECharge::single, HP::CP::undefined));
    h1p.setMass(10);
    auto &h2p = pred.addParticle(HP::BsmParticle("h2p", HP::ECharge::single, HP::CP::undefined));
    h2p.setMass(20);

    SECTION("set and retrieve") {
        for (auto c : magic_enum::enum_values<HP::Collider>()) {
            for (auto &p1 : pred.particles()) {
                for (auto &p2 : pred.particles()) {
                    REQUIRE_NOTHROW(pred.setBsmPairCxn(
                        c, p1.id(), p2.id(),
                        h1.mass() * h2.mass() * static_cast<int>(c)));
                    CHECK(pred.bsmPairCxn(c, p1.id(), p2.id()) ==
                          Approx(h1.mass() * h2.mass() * static_cast<int>(c)));
                    CHECK(pred.bsmPairCxn(c, p2.id(), p1.id()) ==
                          Approx(h1.mass() * h2.mass() * static_cast<int>(c)));
                }
            }
        }
    }

    SECTION("same particle pair production - LHC") {
        REQUIRE_NOTHROW(
            pred.setBsmPairCxn(HP::Collider::LHC13, "h1", "h1", 10));
        REQUIRE(pred.bsmPairCxn(HP::Collider::LHC13, "h1", "h1") == Approx(10));
        REQUIRE(h1.cxn(HP::Collider ::LHC13, HP::Production::pair) ==
                Approx(10));
        h1.setCxn(HP::Collider::LHC13, HP::Production::pair, 20);
        REQUIRE(pred.bsmPairCxn(HP::Collider::LHC13, "h1", "h1") == Approx(20));
        REQUIRE(h1.cxn(HP::Collider ::LHC13, HP::Production::pair) ==
                Approx(20));
    }

    SECTION("same particle pair production - LEP, normalized") {
        REQUIRE_NOTHROW(
            pred.setBsmPairCxn(HP::Collider::LEP, "h1p", "h1p", 10));
        REQUIRE(pred.bsmPairCxn(HP::Collider::LEP, "h1p", "h1p") == Approx(10));
        REQUIRE(h1p.cxn(HP::Collider::LEP, HP::Production::pair) == Approx(10));
        h1p.setNormalizedCxn(HP::Collider::LEP, HP::Production::pair, 20,
                             HP::ReferenceModel::SMHiggs);
        REQUIRE(pred.bsmPairCxn(HP::Collider::LEP, "h1p", "h1p") == Approx(20));
        REQUIRE(h1p.cxn(HP::Collider ::LEP, HP::Production::pair) ==
                Approx(20));
    }

    SECTION("non-existing particles") {
        REQUIRE_THROWS_AS(pred.setBsmPairCxn(HP::Collider::LHC13,
                                             "not-a-particle",
                                             "also-not-a-particle", 10),
                          std::out_of_range);
        REQUIRE(pred.bsmPairCxn(HP::Collider::LHC13, "not-a-particle",
                                "also-not-a-particle") == 0);
        REQUIRE_THROWS_AS(pred.setBsmPairCxn(HP::Collider::LHC13,
                                             "not-a-particle", "not-a-particle",
                                             10),
                          std::out_of_range);
        REQUIRE(pred.bsmPairCxn(HP::Collider::LHC13, "not-a-particle",
                                "not-a-particle") == 0);
    }
}

TEST_CASE("top brs") {
    auto pred = Higgs::Predictions{};

    auto &hp = pred.addParticle(HP::BsmParticle("h+", HP::ECharge::single, HP::CP::undefined));

    hp.setCxn(HP::Collider::LHC13, HP::Production::brtHpb, 0.5);
    REQUIRE(pred.brTopWb() ==
            Approx(1 - hp.cxn(HP::Collider::LHC13, HP::Production::brtHpb)));

    REQUIRE_NOTHROW(pred.setBrTopWb(0.4));
    REQUIRE(pred.brTopWb() == Approx(0.4));

    REQUIRE_THROWS_AS(pred.setBrTopWb(0.7), HP::InvalidInput);
    REQUIRE(pred.brTopWb() == Approx(0.4));

    pred.setBrTopWb(0);
    auto &hp2 = pred.addParticle(HP::BsmParticle("h2+", HP::ECharge::single, HP::CP::undefined));
    hp2.setCxn(HP::Collider::LHC13, HP::Production::brtHpb, 0.4);
    REQUIRE(pred.brTopWb() == Approx(0.1));
    hp2.setCxn(HP::Collider::LHC13, HP::Production::brtHpb, 0.9);
    REQUIRE(pred.brTopWb() == 0);
}
