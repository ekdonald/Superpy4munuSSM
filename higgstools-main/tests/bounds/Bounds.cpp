#include "Higgs/Bounds.hpp"
#include "Higgs/Predictions.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <sstream>

namespace HP = Higgs::predictions;
TEST_CASE("read and apply the testlimits") {

    REQUIRE_THROWS_WITH(Higgs::Bounds{TESTLIMITS_PATH "/invalid"},
                        Catch::Matchers::StartsWith("No valid limit files"));

    auto bounds = Higgs::Bounds{TESTLIMITS_PATH};
    auto pred = Higgs::Predictions{};

    REQUIRE(bounds.limits().size() == 7);

    SECTION("empty predictions") {
        auto res = bounds(pred);
        REQUIRE(res.allowed);
    }

    SECTION("some sensitivity") {
        auto &h = pred.addParticle(HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
        h.setMass(200);
        for (double rate : {1e-2, 1e-1, 1., 10., 100.}) {
            // set rates for the ChannelWidthLimit.json and ChannelLimit.json
            // test limits
            h.setChannelRate(HP::Collider::LHC8, HP::Production::H,
                             HP::Decay::gamgam, rate);
            h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                             HP::Decay::ZZ, rate);
            auto res = bounds(pred);
            REQUIRE(res.selectedLimits.size() == 1);
            CHECK(res.allowed == res.selectedLimits[h.id()].obsRatio() < 1);

            REQUIRE(res.appliedLimits.size() == 3);
            CHECK(res.appliedLimits[0].expRatio() >
                  res.appliedLimits[1].expRatio());
            CHECK(res.appliedLimits[1].expRatio() >
                  res.appliedLimits[2].expRatio());
            auto ss = std::stringstream{};
            ss << res;
            CHECK_THAT(ss.str(),
                       Catch::Matchers::ContainsSubstring("HiggsBounds"));
        }
    }
}
