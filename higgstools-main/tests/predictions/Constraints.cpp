#include "predictions/Constraints.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "predictions/Clustering.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "utilities/Json.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iomanip>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/single.hpp>

using Catch::Approx;
namespace HP = Higgs::predictions;

TEST_CASE("Model likeness") {
    auto h1 = HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined);
    h1.setMass(100);
    auto sm = HP::SMHiggs(125.09);

    const auto p1 = HP::ChannelProcess{HP::Collider::LHC13,
                                       {{HP::Production::ggH, HP::Decay::bb},
                                        {HP::Production::vbfH, HP::Decay::bb},
                                        {HP::Production::HW, HP::Decay::bb},
                                        {HP::Production::HZ, HP::Decay::bb},
                                        {HP::Production::Htt, HP::Decay::bb}}};

    const auto pbr =
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::none, HP::Decay::mumu},
                            {HP::Production::none, HP::Decay::tautau}}};

    const auto pcxn =
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::HW, HP::Decay::none},
                            {HP::Production::HZ, HP::Decay::none},
                            {HP::Production::vbfH, HP::Decay::none}}};

    auto contextDummy = Higgs::Predictions{};

    SECTION("SM-like") {

        auto m1 = HP::ModelLikeness(HP::ReferenceModel::SMHiggs, p1);
        auto mbr = HP::ModelLikeness(HP::ReferenceModel::SMHiggs, pbr);
        auto mcxn = HP::ModelLikeness(HP::ReferenceModel::SMHiggs, pcxn);

        REQUIRE_FALSE(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        REQUIRE_FALSE(mbr.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        REQUIRE_FALSE(mcxn.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));

        for (double m : {800., 50., 333.3, 125.09}) {
            sm.setMass(m);
            INFO(p1(sm) << " " << pbr(sm) << " " << pcxn(sm));
            CHECK(m1.check(HP::ParticleSet{sm}, m, contextDummy));
            CHECK(mbr.check(HP::ParticleSet{sm}, m, contextDummy));
            CHECK(mcxn.check(HP::ParticleSet{sm}, m, contextDummy));
        }

        h1.setMass(sm.mass());
        effectiveCouplingInput(h1, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs);
        CHECK(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        CHECK(mbr.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        CHECK(mcxn.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));

        // compare Fig 4 in 1311.0055 (though we use 13TeV and newer cxns,
        // so the numbers are slightly different)
        auto scalegg = [&h1, &sm](double frac) {
            h1.setNormalizedCxn(HP::Collider::LHC13, HP::Production::ggH, frac,
                                HP::ReferenceModel::SMHiggs);
        };
        scalegg(0.86);
        CHECK(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scalegg(0.84);
        CHECK_FALSE(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scalegg(1.18);
        CHECK(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scalegg(1.2);
        CHECK_FALSE(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scalegg(1.);

        auto scaleVV = [&h1, &sm](double frac) {
            h1.setNormalizedCxn(HP::Collider::LHC13, HP::Production::vbfH, frac,
                                HP::ReferenceModel::SMHiggs);
            h1.setNormalizedCxn(HP::Collider::LHC13, HP::Production::HW, frac,
                                HP::ReferenceModel::SMHiggs);
            h1.setNormalizedCxn(HP::Collider::LHC13, HP::Production::qqHZ, frac,
                                HP::ReferenceModel::SMHiggs);
        };
        scaleVV(0.84);
        CHECK(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scaleVV(0.82);
        CHECK_FALSE(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scaleVV(1.18);
        CHECK(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
        scaleVV(1.2);
        CHECK_FALSE(m1.check(HP::ParticleSet{h1}, h1.mass(), contextDummy));
    }

    SECTION("copy and move") {
        auto m = 200.;
        sm.setMass(m);
        const auto m1 = HP::ModelLikeness(HP::ReferenceModel::SMHiggs, p1);
        REQUIRE(m1.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
        auto mcopy{m1};
        REQUIRE(mcopy.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
        auto mmove{std::move(mcopy)};
        REQUIRE(mmove.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
        REQUIRE_NOTHROW(
            mcopy.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
        mcopy = mmove;
        REQUIRE(mcopy.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
        mmove = std::move(mcopy);
        REQUIRE(mmove.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
        REQUIRE_NOTHROW(
            mcopy.check(HP::ParticleSet{sm}, sm.mass(), contextDummy));
    }

    SECTION("json reads") {
        auto data = nlohmann::json::parse(R"({
            "modelLike": "SMHiggs",
            "process": {
                "charge": "neutral",
                "collider": "LHC13",
                "channels": [
                    ["ggH", "gamgam"],
                    ["vbfH", "gamgam"],
                    ["HW","gamgam"],
                    ["HZ","gamgam"],
                    ["Htt","gamgam"]
                ]
            }
        })");

        INFO(std::setw(4) << data);
        REQUIRE_NOTHROW(
            HP::readChannelProcess(data.at("process"), HP::Collider::LHC8));
        CHECK_NOTHROW(HP::ModelLikeness(data, HP::Collider::LHC8));
        CHECK_NOTHROW(HP::ModelLikeness(data, p1));
        data["process"] = "signal";
        CHECK_NOTHROW(HP::ModelLikeness(data, p1));
        CHECK_THROWS(HP::ModelLikeness(data, HP::Collider::LHC8));
    }
}

TEST_CASE("Top decay Consistency") {
    auto pred = Higgs::Predictions{};
    auto &hp = pred.addParticle(
        HP::BsmParticle("Hp", HP::ECharge::single, HP::CP::undefined));
    auto &h = pred.addParticle(
        HP::BsmParticle("H", HP::ECharge::neutral, HP::CP::undefined));
    hp.setMass(150);
    hp.setCxn(HP::Collider::LHC13, HP::Production::brtHpb, 0.5);
    auto con = HP::TopDecayConsistency{std::vector{HP::Production::brtHpb}};
    REQUIRE(
        con.check(pred.particles() | ranges::to<HP::ParticleSet>, 0., pred));
    h.setCxn(HP::Collider::LHC13, HP::Production::brtHu, 0.1);
    REQUIRE_FALSE(
        con.check(pred.particles() | ranges::to<HP::ParticleSet>, 0., pred));
    h.setCxn(HP::Collider::LHC13, HP::Production::brtHu, 0);
    pred.setBrTopWb(0.499);
    REQUIRE(
        con.check(pred.particles() | ranges::to<HP::ParticleSet>, 0., pred));
    pred.setBrTopWb(0.4);
    REQUIRE(pred.brTopWb() == Approx(0.4));
    REQUIRE_FALSE(
        con.check(pred.particles() | ranges::to<HP::ParticleSet>, 0., pred));
    hp.setCxn(HP::Collider::LHC13, HP::Production::brtHpb, 0.6);
    REQUIRE(
        con.check(pred.particles() | ranges::to<HP::ParticleSet>, 0., pred));

    SECTION("json read") {
        auto data = nlohmann::json::parse(R"({
            "topDecayConsistency": ["brtHpb"]
        })");
        REQUIRE_NOTHROW(HP::TopDecayConsistency{data});
    }
}

TEST_CASE("CPValue") {
    auto predDummy = Higgs::Predictions{};
    auto h0 = HP::BsmParticle("h0", HP::ECharge::neutral, HP::CP::undefined);
    auto H = HP::BsmParticle("H", HP::ECharge::neutral, HP::CP::even);
    auto A = HP::BsmParticle("A", HP::ECharge::neutral, HP::CP::odd);
    auto Hp = HP::BsmParticle("H+", HP::ECharge::single, HP::CP::undefined);
    auto con = HP::CPValue{HP::CP::odd};
    REQUIRE_FALSE(con.check(HP::ParticleSet{h0}, 0., predDummy));
    REQUIRE_FALSE(con.check(HP::ParticleSet{H}, 0., predDummy));
    REQUIRE(con.check(HP::ParticleSet{A}, 0., predDummy));
    REQUIRE_FALSE(con.check(HP::ParticleSet{Hp}, 0., predDummy));
    auto cluster = HP::ParticleSet{H, A};
    REQUIRE_FALSE(con.check(cluster, 0., predDummy));

    SECTION("json read") {
        auto data = nlohmann::json::parse(R"({
            "CPValue": "odd"
        })");
        REQUIRE_NOTHROW(HP::CPValue{data});
    }
}

TEST_CASE("lepton Yukawa") {
    auto predDummy = Higgs::Predictions{};
    auto h = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
    auto conE = HP::MumuTautauRatio{HP::CP::even};
    auto conO = HP::MumuTautauRatio{HP::CP::odd};

    REQUIRE_FALSE(conE.check(HP::ParticleSet{h}, 0., predDummy));
    REQUIRE_FALSE(conO.check(HP::ParticleSet{h}, 0., predDummy));

    for (auto m : {5., 10., 20., 100., 500.}) {
        INFO(m);
        h.setTotalWidth(0.);
        h.setMass(m);
        h.setDecayWidth(HP::Decay::tautau, 1.);
        REQUIRE_FALSE(conE.check(HP::ParticleSet{h}, 0., predDummy));
        REQUIRE_FALSE(conO.check(HP::ParticleSet{h}, 0., predDummy));

        SECTION("just the mass ratio") {
            h.setDecayWidth(
                HP::Decay::mumu,
                std::pow(HP::constants::mMu / HP::constants::mTau, 2));
            if (m >= 10) {
                CHECK(conE.check(HP::ParticleSet{h}, 0., predDummy));
                CHECK(conO.check(HP::ParticleSet{h}, 0., predDummy));
            } else {
                CHECK_FALSE(conE.check(HP::ParticleSet{h}, 0., predDummy));
                CHECK_FALSE(conO.check(HP::ParticleSet{h}, 0., predDummy));
            }
        }

        SECTION("tree-level CP-odd ratio") {
            h.setDecayWidth(
                HP::Decay::mumu,
                std::pow(HP::constants::mMu / HP::constants::mTau, 2) /
                    std::sqrt(1 - std::pow(2 * HP::constants::mTau / m, 2)) *
                    h.br(HP::Decay::tautau) * h.totalWidth());
            if (m >= 10) {
                CHECK(conE.check(HP::ParticleSet{h}, 0., predDummy));
            } else {
                CHECK_FALSE(conE.check(HP::ParticleSet{h}, 0., predDummy));
            }
            CHECK(conO.check(HP::ParticleSet{h}, 0., predDummy));
        }

        SECTION("SM-like CP-even") {
            HP::effectiveCouplingInput(h, HP::smLikeEffCouplings);
            CHECK(conE.check(HP::ParticleSet{h}, 0., predDummy));
            if (m >= 10) {
                CHECK(conO.check(HP::ParticleSet{h}, 0., predDummy));
            } else {
                CHECK_FALSE(conO.check(HP::ParticleSet{h}, 0., predDummy));
            }
        }
    }
}

TEST_CASE("top dominated ggH") {
    auto predDummy = Higgs::Predictions{};
    auto h = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
    auto con = HP::TopDominatedHgg{};
    h.setMass(1e5); // no SM reference
    CHECK_FALSE(con.check(HP::ParticleSet{h}, 0., predDummy));
    h.setMass(100);
    HP::effectiveCouplingInput(h, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs, false);
    CHECK(con.check(HP::ParticleSet{h}, 0., predDummy));
    auto coups = HP::smLikeEffCouplings;
    coups.gg *= 1.1;
    HP::effectiveCouplingInput(h, coups, HP::ReferenceModel::SMHiggs, false);
    CHECK_FALSE(con.check(HP::ParticleSet{h}, 0., predDummy));
    coups = HP::NeutralEffectiveCouplings();
    coups.tt = std::complex<double>(0, 1);
    HP::effectiveCouplingInput(h, coups, HP::ReferenceModel::SMHiggs);
    CHECK(con.check(HP::ParticleSet{h}, 0., predDummy));
    h = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::odd);
    auto ref = HP::SMHiggs(h.mass());
    h.setCxn(HP::Collider::LHC13, HP::Production::Htt,
             HP::EffectiveCouplingRatios::ttHRatio(
                 HP::Collider::LHC13, h.mass(), std::complex<double>(0, 1)) *
                 ref.cxn(HP::Collider::LHC13, HP::Production::Htt));
    h.setCxn(HP::Collider::LHC13, HP::Production::ggH,
             HP::EffectiveCouplingCxns::ggH(HP::Collider::LHC13, h.mass(),
                                            std::complex<double>(0, 1), 0.));
    CHECK(con.check(HP::ParticleSet{h}, 0., predDummy));
}

TEST_CASE("read and check constraints") {
    auto data = nlohmann::json::parse(R"([
        {
            "modelLike": "SMHiggs",
            "process": {
                "charge": "neutral",
                "collider": "LHC13",
                "channels": [
                    ["ggH", "gamgam"],
                    ["vbfH", "gamgam"],
                    ["HW","gamgam"],
                    ["HZ","gamgam"],
                    ["Htt","gamgam"]
                ]
            }
        },{
            "topDecayConsistency": ["brtHc"]
        },{
            "CPValue":"odd"
        },{
            "mumuTautauRatio": "even"
        },{
            "topDominatedHgg": true
        }
    ])");
    REQUIRE_NOTHROW(readConstraints(data, HP::Collider::LHC8));
    const auto proc =
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::ggH, HP::Decay::gamgam},
                            {HP::Production::vbfH, HP::Decay::gamgam}}};
    REQUIRE_NOTHROW(readConstraints(data, proc));

    data.push_back(nlohmann::json{});
    REQUIRE_THROWS_AS(readConstraints(data, HP::Collider::LHC8),
                      Higgs::utilities::BadFieldRead);
    REQUIRE_THROWS_AS(readConstraints(data, proc),
                      Higgs::utilities::BadFieldRead);
}

TEST_CASE("Reference rate") {
    const auto p1 =
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::ggH, HP::Decay::gamgam},
                            {HP::Production::vbfH, HP::Decay::gamgam}}};

    auto r1 = HP::ReferenceRate(HP::ReferenceModel::SMHiggs, p1);
    for (auto m : {20, 130, 222, 415}) {
        CHECK(r1(m) == Approx(p1(HP::SMHiggs(m))));
    }

    const auto pcxn =
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::HZ, HP::Decay::none},
                            {HP::Production::vbfH, HP::Decay::none}}};
    auto r2 = HP::ReferenceRate(HP::ReferenceModel::SMHiggs, pcxn);
    for (auto m : {20, 130, 222, 415}) {
        CHECK(r2(m) == Approx(pcxn(HP::SMHiggs(m))));
    }

    SECTION("copy and move") {
        auto m = 200.;
        auto rcopy{r1};
        REQUIRE(rcopy(m) == r1(m));
        auto rmove{std::move(rcopy)};
        REQUIRE(rmove(m) == r1(m));
        REQUIRE_NOTHROW(rcopy(m));
        rcopy = rmove;
        REQUIRE(rcopy(m) == r1(m));
        rmove = std::move(rcopy);
        REQUIRE(rmove(m) == r1(m));
        REQUIRE_NOTHROW(rcopy(m));
    }

    SECTION("json reads") {
        auto data = nlohmann::json::parse(R"({
            "reference": "SMHiggs",
            "process": {
                "charge": "neutral",
                "collider": "LHC13",
                "channels": [
                    ["ggH", "gamgam"],
                    ["vbfH", "gamgam"],
                    ["HW","gamgam"],
                    ["HZ","gamgam"],
                    ["Htt","gamgam"]
                ]
            }
        })");
        REQUIRE(HP::ReferenceRate::read(data, p1));
        data["process"] = "signal";
        REQUIRE(HP::ReferenceRate::read(data, p1));
        REQUIRE_FALSE(HP::ReferenceRate::read(nlohmann::json{}, p1));
    }

    SECTION("reference model selection") {
        for (auto ref : magic_enum::enum_values<HP::ReferenceModel>()) {
            REQUIRE_NOTHROW(HP::getReference(ref, 100));
        }
        REQUIRE_THROWS_WITH(HP::getReference(HP::ReferenceModel{10}, 10),
                            Catch::Matchers::ContainsSubstring(
                                "unknown", Catch::CaseSensitive::No));
    }
}

namespace {
bool equalParticleSet(const HP::ParticleSet &s1, const HP::ParticleSet &s2) {
    return s1.size() == s2.size() &&
           ranges::all_of(ranges::views::zip_with(
                              [](const HP::Particle &p1,
                                 const HP::Particle &p2) { return p1 == p2; },
                              s1, s2),
                          ranges::identity{});
}

bool clusterListContains(const std::vector<HP::ParticleSet> &clusters,
                         const HP::ParticleSet &targetCluster) {
    return ranges::find_if(clusters,
                           [&targetCluster](const HP::ParticleSet &cluster) {
                               return equalParticleSet(cluster, targetCluster);
                           }) != clusters.end();
}
} // namespace

TEST_CASE("cluster constraints") {
    auto h1 = HP::BsmParticle{"h1", HP::ECharge::neutral, HP::CP::even};
    auto a1 = HP::BsmParticle{"a1", HP::ECharge::neutral, HP::CP::odd};
    auto h2 = HP::BsmParticle{"h2", HP::ECharge::neutral, HP::CP::even};
    auto a2 = HP::BsmParticle{"a2", HP::ECharge::neutral, HP::CP::odd};
    auto hp = HP::BsmParticle{"h+", HP::ECharge::single, HP::CP::undefined};
    h1.setMass(100);
    h2.setMass(100);
    a1.setMass(100);
    a2.setMass(100);
    hp.setMass(100);

    const auto candidates = HP::ParticleSet{h1, a1, h2, a2, hp};
    const auto dummyPred = Higgs::Predictions{};
    auto clusters =
        HP::Clustering::cluster(candidates, {}, HP::MassUncEagerness::ignore);

    REQUIRE(clusters.size() == 31);
    HP::onlyValidClusters(
        clusters, std::vector<HP::Constraint>{HP::CPValue{HP::CP::even}},
        dummyPred);
    INFO(fmt::format(
        "{{{}}}",
        fmt::join(
            clusters | ranges::views::transform([](const HP::ParticleSet &c) {
                return fmt::join(
                    c | ranges::views::transform(
                            [](const HP::Particle &p) { return p.id(); }),
                    ", ");
            }),
            "} {")));
    REQUIRE(clusters.size() == 3);
    CHECK(clusterListContains(clusters, HP::ParticleSet{h1}));
    CHECK(clusterListContains(clusters, HP::ParticleSet{h2}));
    CHECK(clusterListContains(clusters, HP::ParticleSet{h1, h2}));
}
