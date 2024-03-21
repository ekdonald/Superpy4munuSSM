#include "bounds/limits/LikelihoodLimit.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "prettyprint.hpp"
#include "utilities/Json.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <initializer_list>
#include <iomanip>
#include <map>
#include <stdexcept>

using namespace std::string_literals;
using Catch::Approx;
namespace HP = Higgs::predictions;
namespace HB = Higgs::bounds;

TEST_CASE("construct and apply LikelihoodLimit1d") {
    auto id{100U};
    auto reference{"2020.44444"s};
    auto citeKey{"Aad:123456"};
    auto coll{HP::Collider::LHC8};
    auto exp{HP::Experiment::ATLAS};
    auto lumi{30.};
    auto content = nlohmann::json{};
    content["id"] = id;
    content["reference"] = reference;
    content["citeKey"] = citeKey;
    content["collider"] = coll;
    content["experiment"] = exp;
    content["luminosity"] = lumi;

    content["process"] = nlohmann::json::array({});
    content["process"][0] = {{"charge", "neutral"}, {"channels", {}}};
    content["process"][0]["channels"] =
        nlohmann::json::array({{"H", "tautau"}});
    content["analysis"] = {
        {"massResolution", {{"absolute", 10}, {"relative", 0.1}}}};
    content["analysis"]["grid"] = {{"mass", {10, 20, 30}}};
    content["analysis"]["likelihood"] = {
        {"observed", {0., 5.5, 0., 8., 0., 15.}},
        {"expected", {0., 4., 0., 9., 0., 20.}}};
    content["analysis"]["grid"]["channels"] = nlohmann::json::array({{0., 2.}});
    INFO(std::setw(4) << content);

    REQUIRE_NOTHROW(HB::LikelihoodLimit1d::create(content, "nofile"));
    auto lim = HB::LikelihoodLimit1d::create(content, "nofile");
    REQUIRE(lim->id() == id);
    REQUIRE(lim->reference() == reference);
    REQUIRE(lim->citeKey() == citeKey);
    REQUIRE(lim->collider() == coll);
    REQUIRE(lim->experiment() == exp);
    REQUIRE(lim->luminosity() == lumi);

    SECTION("applicable") {
        auto pred = Higgs::Predictions{};
        auto &h = pred.addParticle(HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined});
        h.setMass(15);
        h.setChannelRate(HP::Collider::LHC8, HP::Production::H,
                         HP::Decay::tautau, 0.6);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 1);
        CHECK(res[0].limit()->id() == lim->id());
        CHECK(res[0].expRatio() == Approx(0.5076196547));
        CHECK(res[0].obsRatio() == Approx(0.5204775012));
        CHECK(res[0].expLikelihood() == Approx(1.95));
        CHECK(res[0].obsLikelihood() == Approx(2.025));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});

        SECTION("limit setting eagerness") {
            auto opts = HB::LimitOptions{};
            opts.setLimitMassUnc = HP::MassUncEagerness::eager;
            auto eagerLim =
                HB::LikelihoodLimit1d::create(content, "nofile", opts);

            opts.setLimitMassUnc = HP::MassUncEagerness::ignore;
            auto ignoreLim =
                HB::LikelihoodLimit1d::create(content, "nofile", opts);

            CHECK(eagerLim->apply(pred)[0].expRatio() == res[0].expRatio());
            CHECK(ignoreLim->apply(pred)[0].expRatio() == res[0].expRatio());

            h.setMassUnc(2);
            CHECK(lim->apply(pred)[0].expRatio() < res[0].expRatio());
            CHECK(ignoreLim->apply(pred)[0].expRatio() == res[0].expRatio());
            CHECK(eagerLim->apply(pred)[0].expRatio() > res[0].expRatio());

            SECTION("invalid eagerness default") {
                opts.setLimitMassUnc = HP::MassUncEagerness{4};
                auto invalidLim =
                    HB::LikelihoodLimit1d::create(content, "nofile", opts);
                CHECK(lim->apply(pred)[0].expRatio() ==
                      invalidLim->apply(pred)[0].expRatio());
            }
        }
    }

    SECTION("invalid") {
        SECTION("observed Llh excludes vanishing rates") {
            content["analysis"]["likelihood"]["observed"][0] = 10.;
            content["analysis"]["likelihood"]["observed"][1] = 11.;
            REQUIRE_THROWS_AS(HB::LikelihoodLimit1d::create(content, "nofile"),
                              std::logic_error);
            REQUIRE_THROWS_WITH(
                HB::LikelihoodLimit1d::create(content, "nofile"),
                Catch::Matchers::ContainsSubstring("excludes zero rates"));
        }
        SECTION("expected Llh excludes vanishing rates") {
            content["analysis"]["likelihood"]["expected"][0] = 10.;
            content["analysis"]["likelihood"]["expected"][1] = 11.;
            REQUIRE_THROWS_AS(HB::LikelihoodLimit1d::create(content, "nofile"),
                              std::logic_error);
            REQUIRE_THROWS_WITH(
                HB::LikelihoodLimit1d::create(content, "nofile"),
                Catch::Matchers::ContainsSubstring("excludes zero rates"));
        }
    }

    SECTION("not applicable --- rate") {
        auto pred = Higgs::Predictions{};
        auto h = HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined};
        h.setMass(15);
        pred.addParticle(h);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- mass") {
        auto h = HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined};
        h.setMass(100);
        h.setChannelRate(HP::Collider::LHC8, HP::Production::H,
                         HP::Decay::tautau, 1.5);

        auto pred = Higgs::Predictions{};
        pred.addParticle(std::move(h));
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }
} // namespace Higgs::boundsTEST_CASE("construct and apply LikelihoodLimit1d")

TEST_CASE("construct and apply LikelihoodLimit2d") {
    auto id{100U};
    auto reference{"2020.41414"s};
    auto citeKey{"Aad:123456"};
    auto coll{HP::Collider::LHC13};
    auto exp{HP::Experiment::CMS};
    auto lumi{30.};
    auto content = nlohmann::json{};
    content["id"] = id;
    content["reference"] = reference;
    content["citeKey"] = citeKey;
    content["collider"] = coll;
    content["experiment"] = exp;
    content["luminosity"] = lumi;

    content["process"] = nlohmann::json::array({});
    content["process"][0] = {{"charge", "neutral"}};
    content["process"][0]["channels"] =
        nlohmann::json::array({{"H", "tautau"}});
    content["process"][1] = {{"charge", "neutral"}};
    content["process"][1]["channels"] = nlohmann::json::array({{"H", "bb"}});

    content["analysis"] = {
        {"massResolution", {{"absolute", 10}, {"relative", 0.1}}}};
    content["analysis"]["stackedLlhGrid"] =
        nlohmann::json::parse("[{"
                              "    \"mass\": 10,"
                              "    \"channels\": [[0,2],[0,10]],"
                              "    \"observed\": [0,7,6,15],"
                              "    \"expected\": [0,11,15,9]"
                              "}, {"
                              "    \"mass\": 20,"
                              "    \"channels\": [[0,3],[0,8]],"
                              "    \"observed\": [0,7,6,15],"
                              "    \"expected\": [0,11,15,9]"
                              "}, {"
                              "    \"mass\": 30,"
                              "    \"channels\": [[0,5],[0,4]],"
                              "    \"observed\": [0,7,6,15],"
                              "    \"expected\": [0,11,15,9]"
                              "}]");
    INFO(std::setw(4) << content);

    REQUIRE_NOTHROW(HB::LikelihoodLimit2d::create(content, "nofile"));
    auto lim = HB::LikelihoodLimit2d::create(content, "nofile");
    REQUIRE(lim->id() == id);
    REQUIRE(lim->reference() == reference);
    REQUIRE(lim->citeKey() == citeKey);
    REQUIRE(lim->collider() == coll);
    REQUIRE(lim->experiment() == exp);
    REQUIRE(lim->luminosity() == lumi);

    SECTION("applicable") {
        auto pred = Higgs::Predictions{};
        auto &h = pred.addParticle(HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined});
        h.setMass(15);

        SECTION("both rates") {
            h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                             HP::Decay::tautau, 0.6);
            h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                             HP::Decay::bb, 3.);
            auto res = lim->apply(pred);
            REQUIRE(res.size() == 1);
            CHECK(res[0].limit()->id() == lim->id());
            CHECK(res[0].expRatio() == Approx(1.7318031395));
            CHECK(res[0].obsRatio() == Approx(1.1999417535));
            CHECK(res[0].expLikelihood() == Approx(6.06));
            CHECK(res[0].obsLikelihood() == Approx(4.0275));
            CHECK(res[0].contributingParticles() ==
                  std::vector<std::string>{"h"});
        }
        SECTION("first rate only") {
            h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                             HP::Decay::tautau, 0.6);
            auto res = lim->apply(pred);
            REQUIRE(res.size() == 1);
            CHECK(res[0].limit()->id() == lim->id());
            CHECK(res[0].expRatio() == Approx(0.9761916436));
            CHECK(res[0].obsRatio() == Approx(0.5030925733));
            CHECK(res[0].expLikelihood() == Approx(3.75));
            CHECK(res[0].obsLikelihood() == Approx(1.5));
            CHECK(res[0].contributingParticles() ==
                  std::vector<std::string>{"h"});
        }
        SECTION("second rate only") {
            h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                             HP::Decay::bb, 3.);
            auto res = lim->apply(pred);
            REQUIRE(res.size() == 1);
            CHECK(res[0].limit()->id() == lim->id());
            CHECK(res[0].expRatio() == Approx(0.9664297272));
            CHECK(res[0].obsRatio() == Approx(0.7084619403));
            CHECK(res[0].expLikelihood() == Approx(3.7125));
            CHECK(res[0].obsLikelihood() == Approx(2.3625));
            CHECK(res[0].contributingParticles() ==
                  std::vector<std::string>{"h"});
        }
    }

    SECTION("invalid") {
        SECTION("observed Llh excludes vanishing rates") {
            content["analysis"]["stackedLlhGrid"][0]["observed"][0] = 10.;

            REQUIRE_THROWS_AS(HB::LikelihoodLimit2d::create(content, "nofile"),
                              std::logic_error);
            REQUIRE_THROWS_WITH(
                HB::LikelihoodLimit2d::create(content, "nofile"),
                Catch::Matchers::ContainsSubstring("excludes zero rates"));
        }
        SECTION("expected Llh excludes vanishing rates") {
            content["analysis"]["stackedLlhGrid"][1]["expected"][0] = 10.;
            REQUIRE_THROWS_AS(HB::LikelihoodLimit2d::create(content, "nofile"),
                              std::logic_error);
            REQUIRE_THROWS_WITH(
                HB::LikelihoodLimit2d::create(content, "nofile"),
                Catch::Matchers::ContainsSubstring("excludes zero rates"));
        }
    }

    SECTION("not applicable --- rate") {
        auto h = HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined};
        h.setMass(15);
        auto pred = Higgs::Predictions{};
        pred.addParticle(std::move(h));
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- mass") {
        auto h = HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined};
        h.setMass(100);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                         HP::Decay::tautau, 1.5);

        auto pred = Higgs::Predictions{};
        pred.addParticle(std::move(h));
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }
}
