#include "bounds/limits/WidthLimits.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "utilities/Json.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iomanip>
#include <map>
#include <string>
using namespace std::string_literals;
using Catch::Approx;

namespace HP = Higgs::predictions;
namespace HB = Higgs::bounds;

TEST_CASE("ChannelWidthLimit") {
    auto content = nlohmann::json{};
    content["id"] = 1u;
    content["reference"] = "2020.20200";
    content["citeKey"] = "Aad:1234";
    content["collider"] = HP::Collider::LHC13;
    content["experiment"] = HP::Experiment::ATLAS;
    content["luminosity"] = 30.;
    content["process"] = {{"charge", HP::ECharge::neutral}};
    content["process"]["channels"] = nlohmann::json::array({{"H", "ZZ"}});
    content["analysis"] = {{"relativeWidth", true}};
    content["analysis"]["massResolution"] = {{"absolute", 10},
                                             {"relative", 0.1}};
    content["analysis"]["grid"] = {{"mass", {100, 200, 300}},
                                   {"width", {0., 0.05, 0.3}}};
    content["analysis"]["limit"] = {{"observed", {1, 3, 5, 2, 4, 6, 3, 5, 7}},
                                    {"expected", {2, 4, 6, 3, 5, 7, 4, 6, 8}}};
    {
        INFO(std::setw(4) << content);
        REQUIRE_NOTHROW(HB::ChannelWidthLimit::create(content, "nofile"));
    }
    auto lim = HB::ChannelWidthLimit::create(content, "nofile");
    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined});
    h.setMass(150);
    h.setTotalWidth(10);

    SECTION("applicable") {
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H, HP::Decay::ZZ,
                         2);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 1);
        REQUIRE(res[0].limit()->id() == lim->id());
        CHECK(res[0].expRatio() == Approx(0.4316546763));
        CHECK(res[0].obsRatio() == Approx(0.5504587156));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});
    }

    SECTION("applicable --- normalized") {
        content["normalization"] = {{"reference", "SMHiggs"},
                                    {"process", "signal"}};
        auto limN = HB::ChannelWidthLimit::create(content, "nofile");
        effectiveCouplingInput(h, HP::scaledSMlikeEffCouplings(2.),
                               HP::ReferenceModel::SMHiggs);
        auto res = limN->apply(pred);
        REQUIRE(res.size() == 1);
        CHECK(res[0].limit()->id() == limN->id());
        CHECK(res[0].expRatio() == Approx(1.5885869997));
        CHECK(res[0].obsRatio() == Approx(2.6351139337));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});
    }

    SECTION("applicable --- acceptances") {
        content["analysis"]["acceptances"] = nlohmann::json::parse(R"([
            {"constantAcceptance": 2}
        ])");
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H, HP::Decay::ZZ,
                         2);

        auto limA = HB::ChannelWidthLimit::create(content, "nofile");
        auto res = limA->apply(pred);
        REQUIRE(res.size() == 1);
        CHECK(res[0].limit()->id() == limA->id());
        CHECK(res[0].expRatio() == Approx(2 * 0.4316546763));
        CHECK(res[0].obsRatio() == Approx(2 * 0.5504587156));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});
    }

    SECTION("not applicable --- rate") {
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- width") {
        h.setTotalWidth(100);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }
}
