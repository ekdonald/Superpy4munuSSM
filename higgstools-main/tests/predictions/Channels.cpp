#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Basics.hpp"
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <magic_enum.hpp>
#include <memory>
#include <unordered_map>

using Catch::Approx;
using Higgs::predictions::Production, Higgs::predictions::Decay,
    Higgs::predictions::Collider, Higgs::predictions::ColliderType,
    Higgs::predictions::ECharge;

TEST_CASE("mass resolution") {
    using Higgs::predictions::MassResolution;
    auto res = MassResolution{};
    REQUIRE(res.relative == 0.);
    REQUIRE(res.absolute > 0);
    res = MassResolution{0.1, 10.};
    REQUIRE(res.relative == Approx(0.1));
    REQUIRE(res.absolute == Approx(10.));
}

TEST_CASE("Production and Decay utilities") {
    SECTION("valid production at") {
        CHECK(validProductionAt(Production::none, ColliderType::ee));
        CHECK(validProductionAt(Production::none, ColliderType::pp));
        CHECK(validProductionAt(Production::bbH, ColliderType::pp));
        CHECK_FALSE(validProductionAt(Production::bbH, ColliderType::ee));
        CHECK_FALSE(validProductionAt(Production::bbH, ColliderType{-1}));
    }

    SECTION("valid production for") {
        for (auto c : magic_enum::enum_values<ECharge>())
            CHECK(validProductionFor(Production::none, c));
        CHECK(validProductionFor(Production::ggH, ECharge::neutral));
        CHECK(validProductionFor(Production::HpmW, ECharge::single));
        CHECK_FALSE(validProductionFor(Production::bbH, ECharge::single));
        CHECK_FALSE(validProductionFor(Production::Hpmtb, ECharge::neutral));
        CHECK_FALSE(validProductionFor(Production::Hpmtb, ECharge{-2}));
    }
    SECTION("valid decay for") {
        for (auto c : magic_enum::enum_values<ECharge>())
            CHECK(validDecayFor(Decay::none, c));
        CHECK(validDecayFor(Decay::bb, ECharge::neutral));
        CHECK(validDecayFor(Decay::taunu, ECharge::single));
        CHECK_FALSE(validDecayFor(Decay::gamgam, ECharge::single));
        CHECK_FALSE(validDecayFor(Decay::tb, ECharge::neutral));
        CHECK_FALSE(validDecayFor(Decay::tb, ECharge{-2}));
    }
}
