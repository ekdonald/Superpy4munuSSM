#include "bounds/limits/ProcessLimits.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "prettyprint.hpp"
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

TEST_CASE("constructing a ChannelLimit") {
    auto id{100U};
    auto reference{"XXXX.XXXXX"s};
    auto citeKey{"Aad:123456"};
    auto coll{HP::Collider::LHC13};
    auto exp{HP::Experiment::ATLAS};
    auto lumi{30.};
    auto content = nlohmann::json{};
    content["id"] = id;
    content["reference"] = reference;
    content["citeKey"] = citeKey;
    content["collider"] = coll;
    content["experiment"] = exp;
    content["luminosity"] = lumi;
    content["process"] = {{"charge", HP::ECharge::neutral}};
    content["process"]["channels"] =
        nlohmann::json::array({{"H", "gamgam"}, {"H", "tautau"}});
    content["analysis"] = {
        {"massResolution", {{"absolute", 10}, {"relative", 0.1}}}};
    content["analysis"]["grid"] = {{"mass", {10, 20, 30}}};
    content["analysis"]["limit"] = {{"observed", {1, 3, 5}},
                                    {"expected", {2, 4, 6}}};
    INFO(std::setw(4) << content);

    auto lim = HB::ChannelLimit::create(content, "nofile");
    REQUIRE(lim->id() == id);
    REQUIRE(lim->reference() == reference);
    REQUIRE(lim->citeKey() == citeKey);
    REQUIRE(lim->collider() == coll);
    REQUIRE(lim->experiment() == exp);
    REQUIRE(lim->luminosity() == lumi);
}

TEST_CASE("applying ChannelLimit -- neutral") {
    auto content = nlohmann::json{};
    content["id"] = 1u;
    content["reference"] = "XXXX.XXXXX";
    content["citeKey"] = "Aad:1234";
    content["collider"] = HP::Collider::LHC13;
    content["experiment"] = HP::Experiment::ATLAS;
    content["luminosity"] = 30.;
    content["process"] = {{"charge", HP::ECharge::neutral}};
    content["process"]["channels"] =
        nlohmann::json::array({{"H", "gamgam"}, {"H", "tautau"}});
    content["analysis"] = {
        {"massResolution", {{"absolute", 10}, {"relative", 0.1}}}};
    content["analysis"]["grid"] = {{"mass", {100, 200, 300}}};
    content["analysis"]["limit"] = {{"observed", {1, 3, 5}},
                                    {"expected", {2, 4, 6}}};
    INFO(std::setw(4) << content);
    auto lim = HB::ChannelLimit::create(content, "nofile");

    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::undefined});
    h.setTotalWidth(1.);

    SECTION("applicable") {
        h.setMass(150);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                         HP::Decay::gamgam, 1.5);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 1);
        CHECK(res[0].limit()->id() == lim->id());
        CHECK(res[0].expRatio() == Approx(0.5));
        CHECK(res[0].obsRatio() == Approx(0.75));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});
    }

    SECTION("applicable with acceptances") {
        content["analysis"]["acceptances"] = nlohmann::json::array();
        content["analysis"]["acceptances"][0] = {{"constantAcceptance", 2.}};
        content["analysis"]["acceptances"][1] = {{"constantAcceptance", 3.}};
        auto lima = HB::ChannelLimit::create(content, "nofile");
        h.setMass(150);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                         HP::Decay::gamgam, 1.5);
        auto res = lima->apply(pred);
        REQUIRE(res.size() == 1);
        CHECK(res[0].limit()->id() == lima->id());
        CHECK(res[0].expRatio() == Approx(1.));
        CHECK(res[0].obsRatio() == Approx(1.5));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});
    }

    SECTION("normalized with modelLike") {
        content["normalization"] = {{"reference", "SMHiggs"},
                                    {"process", "signal"}};
        content["constraints"] = nlohmann::json::array(
            {{{"modelLike", "SMHiggs"}, {"process", "signal"}}});
        auto limN = HB::ChannelLimit::create(content, "nofile");
        h.setMass(150);
        effectiveCouplingInput(h, HP::scaledSMlikeEffCouplings(2.),
                               HP::ReferenceModel::SMHiggs);
        auto res = limN->apply(pred);
        REQUIRE(res.size() == 1);
        CHECK(res[0].limit()->id() == limN->id());
        CHECK(res[0].expRatio() == Approx(4 / 3.));
        CHECK(res[0].obsRatio() == Approx(2.));
        CHECK(res[0].contributingParticles() == std::vector<std::string>{"h"});
    }

    SECTION("not applicable --- rate") {
        h.setMass(150);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- mass") {
        h.setMass(10);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                         HP::Decay::gamgam, 1.5);

        auto pred = Higgs::Predictions{};
        pred.addParticle(std::move(h));
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- modelLike signal") {
        content["constraints"] = nlohmann::json::array(
            {{{"modelLike", "SMHiggs"}, {"process", "signal"}}});
        INFO(content);
        auto lim2 = HB::ChannelLimit::create(content, "nofile");
        h.setMass(150);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                         HP::Decay::gamgam, 1.5);
        auto res = lim2->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- modelLike separate proc") {
        content["normalization"] = {{"reference", "SMHiggs"},
                                    {"process", "signal"}};
        content["constraints"] = nlohmann::json::array();
        content["constraints"][0] = {{"modelLike", "SMHiggs"}};
        content["constraints"][0]["process"] = {{"charge", "neutral"}};
        content["constraints"][0]["process"]["channels"] =
            nlohmann::json::array({{"none", "tautau"}, {"none", "mumu"}});
        h.setMass(150);
        effectiveCouplingInput(h, HP::scaledSMlikeEffCouplings(1.));
        auto lim3 = HB::ChannelLimit::create(content, "nofile");
        REQUIRE(lim3->apply(pred).size() == 1);
        h.setBr(HP::Decay::tautau, 0.);
        REQUIRE(lim3->apply(pred).size() == 0);
    }

    SECTION("applicable but irrelevant") {
        h.setMass(150);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::H,
                         HP::Decay::gamgam, 1.5);
        REQUIRE(lim->apply(pred).size() == 1);
        REQUIRE(lim->apply(pred)[0].expRatio() >
                HB::LimitOptions{}.minExpRatio);

        content["analysis"]["limit"]["expected"] = {20000, 40000, 60000};
        auto lim1 = HB::ChannelLimit::create(content, "nofile");
        REQUIRE(lim1->apply(pred).size() == 0);

        auto opts = HB::LimitOptions{};
        opts.minExpRatio = 1e-8;
        auto lim2 = HB::ChannelLimit::create(content, "nofile", opts);
        REQUIRE(lim2->apply(pred).size() == 1);
        REQUIRE(lim2->apply(pred)[0].expRatio() > opts.minExpRatio);
    }
}

TEST_CASE("applying ChannelLimit -- charged") {
    auto content = nlohmann::json{};
    content["id"] = 1u;
    content["reference"] = "XXXX.XXXXX";
    content["citeKey"] = "Aad:1234";
    content["collider"] = HP::Collider::LHC13;
    content["experiment"] = HP::Experiment::ATLAS;
    content["luminosity"] = 30.;
    content["process"] = {{"charge", HP::ECharge::single}};
    content["process"]["channels"] =
        nlohmann::json::array({{"qqHpm", "taunu"}, {"Hpmtb", "tb"}});
    content["analysis"] = {
        {"massResolution", {{"absolute", 10}, {"relative", 0.1}}}};
    content["analysis"]["grid"] = {{"mass", {10, 20, 30}}};
    content["analysis"]["limit"] = {{"observed", {1, 3, 5}},
                                    {"expected", {2, 4, 6}}};
    INFO(std::setw(4) << content);
    auto lim = HB::ChannelLimit::create(content, "nofile");

    SECTION("applicable") {
        auto pred = Higgs::Predictions{};
        auto &h = pred.addParticle(
            HP::BsmParticle{"h", HP::ECharge::single, HP::CP::undefined});
        h.setMass(15);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::qqHpm,
                         HP::Decay::taunu, 1.5);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 1);
        REQUIRE(res[0].limit()->id() == lim->id());
        REQUIRE(res[0].expRatio() == Approx(0.5));
        REQUIRE(res[0].obsRatio() == Approx(0.75));
        REQUIRE(res[0].contributingParticles() ==
                std::vector<std::string>{"h"});
    }

    SECTION("not applicable --- rate") {
        auto pred = Higgs::Predictions{};
        auto &h = pred.addParticle(
            HP::BsmParticle{"h", HP::ECharge::single, HP::CP::undefined});
        h.setMass(15);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }

    SECTION("not applicable --- mass") {
        auto pred = Higgs::Predictions{};
        auto &h = pred.addParticle(
            HP::BsmParticle{"h", HP::ECharge::single, HP::CP::undefined});
        h.setMass(100);
        h.setChannelRate(HP::Collider::LHC13, HP::Production::qqHpm,
                         HP::Decay::taunu, 1.5);
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 0);
    }
}

TEST_CASE("applying ChainDecayLimit") {
    auto content = nlohmann::json{};
    content["id"] = 1u;
    content["reference"] = "XXXX.XXXXX";
    content["citeKey"] = "Aad:1234";
    content["collider"] = HP::Collider::LHC13;
    content["experiment"] = HP::Experiment::ATLAS;
    content["luminosity"] = 30.;

    content["process"] = {{"motherCharge", HP::ECharge::neutral},
                          {"chain", HP::ChainDecay::Z}};
    content["process"]["production"] = nlohmann::json::array({"bbH", "ggH"});
    content["process"]["decay"] = nlohmann::json::array({"gamgam", "ZZ"});

    content["analysis"] = {};
    content["analysis"]["massResolution"] = {
        {"mother", {{"absolute", 10}, {"relative", 0.1}}},
        {"daughter", {{"absolute", 10}, {"relative", 0.1}}}};

    content["analysis"]["grid"] = {{"massMother", {100, 400, 500}},
                                   {"massDaughter", {50, 100, 150}}};
    content["analysis"]["limit"] = {
        {"observed", {10.5, 20.5, 30.5, 11.5, 21.5, 31.5, 12.5, 22.5, 32.5}},
        {"expected", {10, 20, 30, 11, 21, 31, 12, 22, 32}}};

    auto hh = HP::BsmParticle("hh", HP::ECharge::neutral, HP::CP::undefined);
    hh.setMass(320);
    hh.setTotalWidth(1.);
    hh.setCxn(HP::Collider::LHC13, HP::Production::ggH, 30.);
    hh.setBr(HP::ChainDecay::Z, "hl", 0.4);
    hh.setBr(HP::ChainDecay::W, "H+", 0.6);

    auto hl = HP::BsmParticle("hl", HP::ECharge::neutral, HP::CP::undefined);
    hl.setMass(125);
    hl.setTotalWidth(1.);
    hl.setCxn(HP::Collider::LHC13, HP::Production::bbH, 40.);
    hl.setBr(HP::Decay::gamgam, 0.5);
    hl.setBr(HP::ChainDecay::W, "H+", 0.1);

    auto hp = HP::BsmParticle("H+", HP::ECharge::single, HP::CP::undefined);
    hp.setMass(55);
    hp.setTotalWidth(1.);
    hp.setBr(HP::Decay::taunu, 0.25);

    auto pred = Higgs::Predictions{};
    pred.addParticle(std::move(hh));
    pred.addParticle(std::move(hl));
    pred.addParticle(std::move(hp));

    SECTION("Z") {
        {
            INFO(std::setw(4) << content);
            REQUIRE_NOTHROW(HB::ChainDecayLimit::create(content, "nofile"));
        }
        auto lim = HB::ChainDecayLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 1);
        REQUIRE(res[0].limit()->id() == lim->id());
        CHECK(res[0].expRatio() == Approx(0.23316062181));
        CHECK(res[0].obsRatio() == Approx(0.2287166455));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"hh", ">", "hl"});
    }

    SECTION("with production acceptances") {
        content["analysis"]["productionAcceptances"] =
            nlohmann::json::parse(R"([
            {"constantAcceptance": 1e6},
            {"constantAcceptance": 5}
        ])");

        {
            INFO(std::setw(4) << content);
            REQUIRE_NOTHROW(HB::ChainDecayLimit::create(content, "nofile"));
        }
        auto lim = HB::ChainDecayLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 1);
        REQUIRE(res[0].limit()->id() == lim->id());
        CHECK(res[0].expRatio() == Approx(5 * 0.23316062181));
        CHECK(res[0].obsRatio() == Approx(5 * 0.2287166455));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"hh", ">", "hl"});
    }

    SECTION("W") {
        content["process"]["chain"] = HP::ChainDecay::W;
        content["process"]["decay"] = nlohmann::json::array({HP::Decay::taunu});
        auto lim = HB::ChainDecayLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        REQUIRE(res.size() == 2);
        REQUIRE(res[0].limit()->id() == lim->id());
        CHECK(res[0].expRatio() == Approx(0.3835227273));
        CHECK(res[0].obsRatio() == Approx(0.3678474114));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"hh", ">", "H+"});
        CHECK(res[1].expRatio() == Approx(0.0902255639));
        CHECK(res[1].obsRatio() == Approx(0.0863309353));
        CHECK(res[1].contributingParticles() ==
              std::vector<std::string>{"hl", ">", "H+"});
    }
}

TEST_CASE("applying PairDecayLimit") {
    auto content = nlohmann::json{};
    content["id"] = 1u;
    content["reference"] = "XXXX.XXXXX";
    content["citeKey"] = "Aad:1234";
    content["collider"] = HP::Collider::LHC13;
    content["experiment"] = HP::Experiment::ATLAS;
    content["luminosity"] = 30.;

    content["process"] = nlohmann::json::parse(" {"
                                               "     \"collider\": \"LHC13\","
                                               "     \"production\": ["
                                               "         \"ggH\", \"bbH\""
                                               "     ],"
                                               "     \"firstDecay\": ["
                                               "         \"gamgam\", \"mumu\""
                                               "     ],"
                                               "     \"secondDecay\": ["
                                               "         \"bb\", \"mumu\""
                                               "     ]"
                                               " }");

    content["analysis"] = {};
    content["analysis"]["massResolution"] = {
        {"mother", {{"absolute", 10}, {"relative", 0.1}}},
        {"firstDaughter", {{"absolute", 10}, {"relative", 0.1}}},
        {"secondDaughter", {{"absolute", 10}, {"relative", 0.1}}}};

    content["analysis"]["grid"] = {{"massMother", {100, 400, 500}},
                                   {"massFirstDaughter", {50, 100, 150}},
                                   {"massSecondDaughter", {60, 120, 125}}};
    content["analysis"]["limit"] = {
        {"observed",
         {
             10.,  20.,  30.,  11.,  21.,  31.,  12.,  22.,  32.,
             10.5, 20.5, 30.5, 11.5, 21.5, 31.5, 12.5, 22.5, 32.5,
             10.9, 20.9, 30.9, 11.9, 21.9, 31.9, 12.9, 22.9, 32.9,
         }},
        {"expected",
         {
             10.,  20.,  30.,  11.,  21.,  31.,  12.,  22.,  32.,
             10.5, 20.5, 30.5, 11.5, 21.5, 31.5, 12.5, 22.5, 32.5,
             10.9, 20.9, 30.9, 11.9, 21.9, 31.9, 12.9, 22.9, 32.9,
         }}};

    auto hh = HP::BsmParticle("hh", HP::ECharge::neutral, HP::CP::undefined);
    hh.setMass(220);
    hh.setTotalWidth(1.);
    hh.setCxn(HP::Collider::LHC13, HP::Production::ggH, 1.);
    hh.setBr("h1", "h1", 0.4);
    hh.setBr("h1", "h2", 0.3);

    auto h1 = HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined);
    h1.setMass(125);
    h1.setTotalWidth(1.);
    h1.setBr(HP::Decay::gamgam, 0.5);
    h1.setBr(HP::Decay::bb, 0.2);
    h1.setBr(HP::Decay::mumu, 0.1);

    auto h2 = HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined);
    h2.setMass(140);
    h2.setTotalWidth(1.);
    h2.setBr(HP::Decay::mumu, 0.1);
    h2.setBr(HP::Decay::bb, 0.7);

    auto hhp = HP::BsmParticle("H+", HP::ECharge::single, HP::CP::undefined);
    hhp.setMass(300);
    hhp.setTotalWidth(1.);
    hhp.setCxn(HP::Collider::LHC13, HP::Production::vbfHpm, 1.);
    hhp.setBr("h1", "h+", 0.8);
    hhp.setBr("h2", "h+", 0.1);

    auto hlp = HP::BsmParticle("h+", HP::ECharge::single, HP::CP::undefined);
    hlp.setMass(30);
    hlp.setTotalWidth(1.);
    hlp.setBr(HP::Decay::taunu, 0.25);

    auto pred = Higgs::Predictions{};
    pred.addParticle(std::move(hh));
    pred.addParticle(std::move(h1));
    pred.addParticle(std::move(h2));
    pred.addParticle(std::move(hhp));
    pred.addParticle(std::move(hlp));

    SECTION("distinct daughter grids") {
        INFO(std::setw(4) << content);
        REQUIRE_NOTHROW(HB::PairDecayLimit::create(content, "nofile"));
        auto lim = HB::PairDecayLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        for (const auto &r : res) {
            REQUIRE(r.limit()->id() == lim->id());
            CHECK(r.expRatio() ==
                  Approx(r.obsRatio())); // identical limits for simplicity
        }
        REQUIRE(res.size() == 1);
        CHECK(res[0].expRatio() == Approx(0.0046976301));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"hh", ">", "h1", "h2", "+", "h1"});
    }
    SECTION("equal daughter masses") {
        content["analysis"]["equalDaughterMasses"] = true;
        content["analysis"]["grid"].erase("massSecondDaughter");
        content["analysis"]["limit"]["observed"] =
            std::vector{10., 20., 30., 11., 21., 31., 12., 22., 32.};
        content["analysis"]["limit"]["expected"] =
            std::vector{10., 20., 30., 11., 21., 31., 12., 22., 32.};
        INFO(std::setw(4) << content);
        REQUIRE_NOTHROW(HB::PairDecayLimit::create(content, "nofile"));
        auto lim = HB::PairDecayLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        for (const auto &r : res) {
            REQUIRE(r.limit()->id() == lim->id());
            CHECK(r.expRatio() ==
                  Approx(r.obsRatio())); // identical limits for simplicity
            WARN(r.contributingParticles());
        }
        REQUIRE(res.size() == 1);
        CHECK(res[0].expRatio() == Approx(0.0114937343));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"hh", ">", "h1", "h2", "+", "h1", "h2"});
    }
}

TEST_CASE("applying PairProductionLimit") {
    auto content = nlohmann::json{};
    content["id"] = 1u;
    content["reference"] = "XXXX.XXXXX";
    content["citeKey"] = "Aad:1234";
    content["collider"] = HP::Collider::LEP;
    content["experiment"] = HP::Experiment::ATLAS;
    content["luminosity"] = 40.;

    content["process"] = nlohmann::json::parse(R"({
                                                "firstDecay": [
                                                    "gamgam"
                                                ],
                                                "secondDecay": [
                                                    "bb", "mumu"
                                                ]
                                               })");

    content["analysis"] = {};
    content["analysis"]["massResolution"] = {
        {"firstParticle", {{"absolute", 1}, {"relative", 0.}}},
        {"secondParticle", {{"absolute", 10}, {"relative", 0.1}}}};

    content["analysis"]["grid"] = {{"massFirstParticle", {100, 150}},
                                   {"massSecondParticle", {120, 130}}};
    content["analysis"]["limit"] = {{"observed", {10., 30., 10.5, 20.5}},
                                    {"expected", {10., 30., 10.5, 20.5}}};

    auto pred = Higgs::Predictions{};

    auto &h1 = pred.addParticle(
        HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
    h1.setMass(125);
    h1.setTotalWidth(1.);
    h1.setBr(HP::Decay::gamgam, 0.5);
    h1.setBr(HP::Decay::bb, 0.2);
    h1.setBr(HP::Decay::mumu, 0.1);

    auto &h2 = pred.addParticle(
        HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
    h2.setMass(130);
    h2.setTotalWidth(1.);
    h2.setBr(HP::Decay::mumu, 0.1);
    h2.setBr(HP::Decay::bb, 0.7);

    pred.setBsmPairCxn(HP::Collider::LEP, "h1", "h2", 10.);
    pred.setBsmPairCxn(HP::Collider::LEP, "h1", "h1", 5.);
    pred.setBsmPairCxn(HP::Collider::LEP, "h2", "h2", 1e6);

    SECTION("distinct masses") {
        INFO(std::setw(4) << content);
        REQUIRE_NOTHROW(HB::PairProductionLimit::create(content, "nofile"));

        auto lim = HB::PairProductionLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        for (const auto &r : res) {
            REQUIRE(r.limit()->id() == lim->id());
            CHECK(r.expRatio() ==
                  Approx(r.obsRatio())); // identical limits for simplicity
            WARN(r.contributingParticles());
        }
        REQUIRE(res.size() == 1);
        CHECK(res[0].expRatio() == Approx(0.2370225269));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"h1", "+", "h1", "h2"});
    }

    SECTION("equal masses") {
        content["analysis"]["equalParticleMasses"] = true;
        content["analysis"]["grid"].erase("massSecondParticle");
        content["analysis"]["limit"]["observed"] = std::vector{10., 30.};
        content["analysis"]["limit"]["expected"] = std::vector{10., 30.};
        INFO(std::setw(4) << content);
        REQUIRE_NOTHROW(HB::PairProductionLimit::create(content, "nofile"));

        auto lim = HB::PairProductionLimit::create(content, "nofile");
        auto res = lim->apply(pred);
        for (const auto &r : res) {
            REQUIRE(r.limit()->id() == lim->id());
            CHECK(r.expRatio() ==
                  Approx(r.obsRatio())); // identical limits for simplicity
            WARN(r.contributingParticles());
        }
        REQUIRE(res.size() == 1);
        CHECK(res[0].expRatio() == Approx(0.275));
        CHECK(res[0].contributingParticles() ==
              std::vector<std::string>{"h1", "+", "h1", "h2"});
    }
}
