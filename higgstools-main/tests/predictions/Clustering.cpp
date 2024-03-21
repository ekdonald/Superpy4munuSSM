#include "predictions/Clustering.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "prettyprint.hpp" // IWYU pragma: keep
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/zip_with.hpp>

namespace HP = Higgs::predictions;
namespace HPC = Higgs::predictions::Clustering;
using Catch::Approx;
namespace {
bool particleSetContains(const HP::ParticleSet &set, const HP::Particle &elem) {
    return set.find(std::cref(elem)) != set.end();
}

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

TEST_CASE("particles in mass range") {
    auto pred = Higgs::Predictions{};
    auto &h1 = pred.addParticle(HP::BsmParticle{"h1", HP::ECharge::neutral, HP::CP::undefined});
    auto &h2 = pred.addParticle(HP::BsmParticle{"h2", HP::ECharge::neutral, HP::CP::undefined});
    auto &h3 = pred.addParticle(HP::BsmParticle{"h3", HP::ECharge::neutral, HP::CP::undefined});
    auto &h4 = pred.addParticle(HP::BsmParticle{"h4", HP::ECharge::neutral, HP::CP::undefined});
    h1.setMass(1);
    h1.setMassUnc(1.5);
    h2.setMass(2);
    h2.setMassUnc(0.5);
    h3.setMass(3);
    h3.setMassUnc(2);
    h4.setMass(4);

    SECTION("empty") {
        auto emptyPred = Higgs::Predictions{};
        for (auto eag : magic_enum::enum_values<HP::MassUncEagerness>()) {
            auto emptyR = HPC::particlesInMassRange(pred, {10, 20}, eag);
            REQUIRE(emptyR.empty());
            auto emptyP = HPC::particlesInMassRange(emptyPred, {0, 10}, eag);
            REQUIRE(emptyP.empty());
        }
    }

    auto eager = HPC::particlesInMassRange(pred, {1.4, 4.1},
                                           HP::MassUncEagerness::eager);
    REQUIRE(eager.size() == 4);
    CHECK(particleSetContains(eager, h1));
    CHECK(particleSetContains(eager, h2));
    CHECK(particleSetContains(eager, h3));
    CHECK(particleSetContains(eager, h4));

    auto facEager = HPC::particlesInMassRange(pred, {1.4, 4.1},
                                              HP::MassUncEagerness::eager, 0.2);
    REQUIRE(facEager.size() == 3);
    CHECK_FALSE(particleSetContains(facEager, h1));
    CHECK(particleSetContains(facEager, h2));
    CHECK(particleSetContains(facEager, h3));
    CHECK(particleSetContains(facEager, h4));

    auto cautious = HPC::particlesInMassRange(pred, {1.4, 4.1},
                                              HP::MassUncEagerness::cautious);
    REQUIRE(cautious.size() == 2);
    CHECK_FALSE(particleSetContains(cautious, h1));
    CHECK(particleSetContains(cautious, h2));
    CHECK_FALSE(particleSetContains(cautious, h3));
    CHECK(particleSetContains(cautious, h4));

    auto facCautious = HPC::particlesInMassRange(
        pred, {1.4, 4.1}, HP::MassUncEagerness::cautious, 0.5);
    REQUIRE(facCautious.size() == 3);
    CHECK_FALSE(particleSetContains(facCautious, h1));
    CHECK(particleSetContains(facCautious, h2));
    CHECK(particleSetContains(facCautious, h3));
    CHECK(particleSetContains(facCautious, h4));

    auto ignore = HPC::particlesInMassRange(pred, {1.4, 4.1},
                                            HP::MassUncEagerness::ignore);
    REQUIRE(ignore.size() == 3);
    CHECK_FALSE(particleSetContains(ignore, h1));
    CHECK(particleSetContains(ignore, h2));
    CHECK(particleSetContains(ignore, h3));
    CHECK(particleSetContains(ignore, h4));
}

TEST_CASE("filter relevant particles") {
    auto h1 = HP::BsmParticle{"h1", HP::ECharge::neutral, HP::CP::undefined};
    auto h2 = HP::BsmParticle{"h2", HP::ECharge::neutral, HP::CP::undefined};
    auto h3 = HP::BsmParticle{"h3", HP::ECharge::neutral, HP::CP::undefined};
    auto h4 = HP::BsmParticle{"h4", HP::ECharge::neutral, HP::CP::undefined};

    h1.setCxn(HP::Collider::LHC13, HP::Production::ggH, 1.);
    h1.setDecayWidth(HP::Decay::bb, 1);
    h2.setCxn(HP::Collider::LHC13, HP::Production::ggH,
              HP::constants::minimumRate / 10);
    h3.setCxn(HP::Collider::LHC13, HP::Production::Htt, 100.);
    h3.setDecayWidth(HP::Decay::bb, 1);
    h4.setCxn(HP::Collider::LHC13, HP::Production::ggH, 10);

    auto particles = HP::ParticleSet{h1, h2, h3, h4};

    SECTION("1 particle process") {
        const auto proc = [](const HP::Particle &p) {
            return p.cxn(HP::Collider::LHC13, HP::Production::ggH);
        };

        auto candidates = std::array{particles};
        const auto &[ps1] = HPC::filterRelevantParticles(candidates, proc);
        REQUIRE(ps1.size() == 2);
        CHECK(particleSetContains(ps1, h1));
        CHECK(particleSetContains(ps1, h4));
    }

    SECTION("2 particle process") {
        const auto proc = [](const HP::Particle &p1, const HP::Particle &p2) {
            return p1.cxn(HP::Collider::LHC13, HP::Production::ggH) *
                   p2.br(HP::Decay::bb);
        };

        auto candidates = std::array{particles, particles};
        const auto &[ps1, ps2] = HPC::filterRelevantParticles(candidates, proc);
        REQUIRE(ps1.size() == 2);
        CHECK(ps1.find(h1) != ps1.end());
        CHECK(ps1.find(h4) != ps1.end());
        REQUIRE(ps2.size() == 2);
        CHECK(ps2.find(h1) != ps2.end());
        CHECK(ps2.find(h3) != ps2.end());
    }

    SECTION("3 particle process") {
        const auto proc = [](const HP::Particle &p1, const HP::Particle &p2,
                             const HP::Particle &p3) {
            return p1.cxn(HP::Collider::LHC13, HP::Production::ggH) *
                   p3.cxn(HP::Collider::LHC13, HP::Production::Htt) *
                   p2.br(HP::Decay::bb);
        };

        auto candidates = std::array{particles, particles, particles};
        const auto &[ps1, ps2, ps3] =
            HPC::filterRelevantParticles(candidates, proc);
        REQUIRE(ps1.size() == 3);
        CHECK(ps1.find(h1) != ps1.end());
        // the h3 Htt cxn is large enough to make h2 relevant
        CHECK(ps1.find(h2) != ps1.end());
        CHECK(ps1.find(h4) != ps1.end());
        REQUIRE(ps2.size() == 2);
        CHECK(ps2.find(h1) != ps2.end());
        CHECK(ps2.find(h3) != ps2.end());
        REQUIRE(ps3.size() == 1);
        CHECK(ps2.find(h3) != ps2.end());
    }

    SECTION("zero process") {
        const auto proc = [](const auto &...) { return 0; };

        const auto candidates1 = std::array{particles};
        const auto &[ps11] = HPC::filterRelevantParticles(candidates1, proc);
        REQUIRE(ps11.empty());

        const auto candidates2 = std::array{particles, particles};
        const auto &[ps21, ps22] =
            HPC::filterRelevantParticles(candidates2, proc);
        REQUIRE(ps21.empty());
        REQUIRE(ps22.empty());

        const auto candidates3 = std::array{particles, particles, particles};
        const auto &[ps31, ps32, ps33] =
            HPC::filterRelevantParticles(candidates3, proc);
        REQUIRE(ps31.empty());
        REQUIRE(ps32.empty());
        REQUIRE(ps33.empty());

        const auto candidates4 =
            std::array{particles, particles, particles, particles};
        const auto &[ps41, ps42, ps43, ps44] =
            HPC::filterRelevantParticles(candidates4, proc);
        REQUIRE(ps41.empty());
        REQUIRE(ps42.empty());
        REQUIRE(ps43.empty());
        REQUIRE(ps44.empty());
    }
}

TEST_CASE("clustering") {
    auto m = HP::BsmParticle{"m", HP::ECharge::neutral, HP::CP::undefined};
    auto mu = HP::BsmParticle{"mu", HP::ECharge::neutral, HP::CP::undefined};
    auto muu = HP::BsmParticle{"muu", HP::ECharge::neutral, HP::CP::undefined};
    auto m2 = HP::BsmParticle{"m2", HP::ECharge::neutral, HP::CP::undefined};
    auto m3 = HP::BsmParticle{"m3", HP::ECharge::neutral, HP::CP::undefined};
    auto m4 = HP::BsmParticle{"m4", HP::ECharge::neutral, HP::CP::undefined};
    m.setMass(100);
    mu.setMass(100);
    mu.setMassUnc(5);
    muu.setMass(100);
    muu.setMassUnc(20);
    m2.setMass(108);
    m2.setMassUnc(3.5);
    m3.setMass(112);
    m3.setMassUnc(5);
    m4.setMass(200);
    auto candidates = HP::ParticleSet{m, mu, muu, m2, m3, m4};
    SECTION("zero resolution") {
        std::vector<HP::ParticleSet> eagerClusters =
            HPC::cluster(candidates, {}, HP::MassUncEagerness::eager);
        REQUIRE(eagerClusters.size() == 16);
        auto cautiousClusters =
            HPC::cluster(candidates, {}, HP::MassUncEagerness::cautious);
        REQUIRE(cautiousClusters.size() == 6);
        auto centralClusters =
            HPC::cluster(candidates, {}, HP::MassUncEagerness::ignore);
        REQUIRE(centralClusters.size() == 10);

        SECTION("no empty clusters") {
            for (auto &cluster : ranges::views::concat(
                     eagerClusters, cautiousClusters, centralClusters)) {
                CHECK(cluster.size() >= 1);
            }
        }

        SECTION("single particle clusters") {
            for (auto &clusters :
                 {eagerClusters, cautiousClusters, centralClusters}) {
                CHECK(clusterListContains(clusters, HP::ParticleSet{m}));
                CHECK(clusterListContains(clusters, HP::ParticleSet{mu}));
                CHECK(clusterListContains(clusters, HP::ParticleSet{muu}));
                CHECK(clusterListContains(clusters, HP::ParticleSet{m2}));
                CHECK(clusterListContains(clusters, HP::ParticleSet{m3}));
                CHECK(clusterListContains(clusters, HP::ParticleSet{m4}));
            }
        }

        SECTION("same mass clusters") {
            auto sameMassClusters = HPC::cluster(
                HP::ParticleSet{m, mu, muu}, {}, HP::MassUncEagerness::ignore);
            // the -3 is for {m2}, {m3}, and {m4}
            REQUIRE(sameMassClusters.size() == centralClusters.size() - 3);
            for (const auto &clusters : {eagerClusters, centralClusters}) {
                for (const auto &sameMassCluster : sameMassClusters) {
                    CHECK(clusterListContains(clusters, sameMassCluster));
                }
            }

            CHECK_FALSE(
                clusterListContains(cautiousClusters, HP::ParticleSet{m, mu}));
        }

        SECTION("remaining eager clusters") {
            REQUIRE(eagerClusters.size() == centralClusters.size() + 6);
            CHECK(clusterListContains(eagerClusters, HP::ParticleSet{mu, m2}));
            CHECK(clusterListContains(eagerClusters, HP::ParticleSet{muu, m2}));
            CHECK(clusterListContains(eagerClusters,
                                      HP::ParticleSet{mu, muu, m2}));
            CHECK(clusterListContains(eagerClusters, HP::ParticleSet{m3, m2}));
            CHECK(clusterListContains(eagerClusters,
                                      HP::ParticleSet{muu, m2, m3}));
            CHECK(clusterListContains(eagerClusters, HP::ParticleSet{muu, m3}));
        }
    }

    SECTION("finite resolution") {
        std::vector<HP::ParticleSet> eagerClusters =
            HPC::cluster(candidates, {0.05, 5}, HP::MassUncEagerness::eager);
        REQUIRE(eagerClusters.size() == 16 + 16);
        auto cautiousClusters =
            HPC::cluster(candidates, {0.05, 5}, HP::MassUncEagerness::cautious);
        REQUIRE(cautiousClusters.size() == 6 + 1);
        auto centralClusters =
            HPC::cluster(candidates, {0.05, 5}, HP::MassUncEagerness::ignore);
        REQUIRE(centralClusters.size() == 10 + 8);

        SECTION("no empty clusters") {
            for (auto &cluster : ranges::views::concat(
                     eagerClusters, cautiousClusters, centralClusters)) {
                CHECK(cluster.size() >= 1);
            }
        }

        SECTION("all zero resolution clusters") {
            for (const auto &zeroResCluster :
                 HPC::cluster(candidates, {}, HP::MassUncEagerness::eager)) {
                CHECK(clusterListContains(eagerClusters, zeroResCluster));
            }
            for (const auto &zeroResCluster :
                 HPC::cluster(candidates, {}, HP::MassUncEagerness::cautious)) {
                CHECK(clusterListContains(cautiousClusters, zeroResCluster));
            }
            for (const auto &zeroResCluster :
                 HPC::cluster(candidates, {}, HP::MassUncEagerness::ignore)) {
                CHECK(clusterListContains(centralClusters, zeroResCluster));
            }
        }

        SECTION("fully contained in mass resolution") {
            for (auto &clusters :
                 {eagerClusters, cautiousClusters, centralClusters}) {
                CHECK(clusterListContains(clusters, HP::ParticleSet{m, mu}));
            }
        }

        SECTION("central values in unc") {
            auto clustersToCheck = std::vector<HP::ParticleSet>{
                {m, m2},     {mu, m2},     {muu, m2},     {m2, m3},
                {m, mu, m2}, {m, muu, m2}, {mu, muu, m2}, {m, mu, muu, m2}};
            REQUIRE(clustersToCheck.size() == 8);
            for (auto &clusters : {eagerClusters, centralClusters}) {
                for (auto &targetCluster : clustersToCheck) {
                    CHECK(clusterListContains(clusters, targetCluster));
                }
            }
        }

        SECTION("everything except m4 clusters eagerly") {
            auto almostAllClusters =
                HPC::cluster(HP::ParticleSet{m, mu, muu, m2, m3}, {100},
                             HP::MassUncEagerness::ignore);
            // -1 for {m4}
            REQUIRE(almostAllClusters.size() == eagerClusters.size() - 1);
            for (auto &targetCluster : almostAllClusters) {
                CHECK(clusterListContains(eagerClusters, targetCluster));
            }
        }
    }
}

TEST_CASE("Maximal clusters") {
    const auto h1 = HP::BsmParticle{"h1", HP::ECharge::neutral, HP::CP::undefined};
    const auto h2 = HP::BsmParticle{"h2", HP::ECharge::neutral, HP::CP::undefined};
    const auto h3 = HP::BsmParticle{"h3", HP::ECharge::neutral, HP::CP::undefined};
    const auto h4 = HP::BsmParticle{"h4", HP::ECharge::neutral, HP::CP::undefined};

    SECTION("edge cases") {
        auto empty = std::vector<HP::ParticleSet>{};
        HPC::onlyMaximalClusters(empty);
        REQUIRE(empty.empty());

        auto single = std::vector<HP::ParticleSet>{{h1}};
        HPC::onlyMaximalClusters(single);
        REQUIRE(single.size() == 1);
        CHECK(particleSetContains(single[0], h1));
    }

    auto testClusters = std::vector<HP::ParticleSet>{
        {h1},     {h2},     {h3},         {h4},        {h1, h2},
        {h1, h3}, {h2, h4}, {h1, h2, h4}, {h2, h3, h4}};
    HPC::onlyMaximalClusters(testClusters);
    REQUIRE(testClusters.size() == 3);
    CHECK(particleSetContains(testClusters[0], h1));
    CHECK(particleSetContains(testClusters[0], h3));
    CHECK(particleSetContains(testClusters[1], h1));
    CHECK(particleSetContains(testClusters[1], h2));
    CHECK(particleSetContains(testClusters[1], h4));
    CHECK(particleSetContains(testClusters[2], h2));
    CHECK(particleSetContains(testClusters[2], h3));
    CHECK(particleSetContains(testClusters[2], h4));
}

TEST_CASE("perform clustering") {
    auto m = HP::BsmParticle{"m", HP::ECharge::neutral, HP::CP::even};
    auto mu = HP::BsmParticle{"mu", HP::ECharge::neutral, HP::CP::odd};
    auto muu = HP::BsmParticle{"muu", HP::ECharge::neutral, HP::CP::even};
    auto m2 = HP::BsmParticle{"m2", HP::ECharge::neutral, HP::CP::even};
    auto m3 = HP::BsmParticle{"m3", HP::ECharge::neutral, HP::CP::even};
    auto m4 = HP::BsmParticle{"m4", HP::ECharge::neutral, HP::CP::even};
    m.setMass(100);
    mu.setMass(100);
    mu.setMassUnc(5);
    muu.setMass(100);
    muu.setMassUnc(20);
    m2.setMass(108);
    m2.setMassUnc(3.5);
    m3.setMass(112);
    m3.setMassUnc(5);
    m4.setMass(200);
    auto candidates = HP::ParticleSet{m, mu, muu, m2, m3, m4};
    const auto predDummy = Higgs::Predictions{};

    SECTION("1 role") {
        auto cmpClusters =
            HPC::cluster(candidates, {0.05, 5}, HP::MassUncEagerness::eager);
        INFO(fmt::format(
            "all clusters: {{{}}}",
            fmt::join(
                cmpClusters |
                    ranges::views::transform([](const HP::ParticleSet &c) {
                        return fmt::join(c | ranges::views::transform(
                                                 [](const HP::Particle &p) {
                                                     return p.id();
                                                 }),
                                         ", ");
                    }),
                "} {")));
        HP::onlyValidClusters(
            cmpClusters, std::vector<HP::Constraint>{HP::CPValue{HP::CP::even}},
            predDummy);
        INFO(fmt::format(
            "after constraints: {{{}}}",
            fmt::join(
                cmpClusters |
                    ranges::views::transform([](const HP::ParticleSet &c) {
                        return fmt::join(c | ranges::views::transform(
                                                 [](const HP::Particle &p) {
                                                     return p.id();
                                                 }),
                                         ", ");
                    }),
                "} {")));
        HPC::onlyMaximalClusters(cmpClusters);
        INFO(fmt::format(
            "maximal cluster comparison: {{{}}}",
            fmt::join(
                cmpClusters |
                    ranges::views::transform([](const HP::ParticleSet &c) {
                        return fmt::join(c | ranges::views::transform(
                                                 [](const HP::Particle &p) {
                                                     return p.id();
                                                 }),
                                         ", ");
                    }),
                "} {")));
        auto clusters1Role = HPC::performClusteringNew(
            std::array{candidates}, std::array{HP::MassResolution{0.05, 5}},
            HP::MassUncEagerness::eager,
            std::array{std::vector<HP::Constraint>{HP::CPValue{HP::CP::even}}},
            predDummy);
        const auto &[c1] = clusters1Role;
        INFO(fmt::format(
            "{{{}}}",
            fmt::join(
                c1 | ranges::views::transform([](const HP::ParticleSet &c) {
                    return fmt::join(
                        c | ranges::views::transform(
                                [](const HP::Particle &p) { return p.id(); }),
                        ", ");
                }),
                "} {")));
        REQUIRE(c1.size() == 2);
        CHECK(clusterListContains(c1, HP::ParticleSet{m4}));
        CHECK(clusterListContains(c1, HP::ParticleSet{m, muu, m2, m3}));
    }

    SECTION("2 roles, just to be sure") {
        auto clusters2Roles = HPC::performClusteringNew(
            std::array{candidates, candidates},
            std::array{HP::MassResolution{0.05, 5}, HP::MassResolution{}},
            HP::MassUncEagerness::eager,
            std::array{std::vector<HP::Constraint>{HP::CPValue{HP::CP::even}},
                       std::vector<HP::Constraint>{}},
            predDummy);
        const auto &[c1, c2] = clusters2Roles;
        INFO(fmt::format(
            "c1: {{{}}}",
            fmt::join(
                c1 | ranges::views::transform([](const HP::ParticleSet &c) {
                    return fmt::join(
                        c | ranges::views::transform(
                                [](const HP::Particle &p) { return p.id(); }),
                        ", ");
                }),
                "} {")));
        REQUIRE(c1.size() == 2);
        CHECK(clusterListContains(c1, HP::ParticleSet{m4}));
        CHECK(clusterListContains(c1, HP::ParticleSet{m, muu, m2, m3}));
        INFO(fmt::format(
            "c2: {{{}}}",
            fmt::join(
                c2 | ranges::views::transform([](const HP::ParticleSet &c) {
                    return fmt::join(
                        c | ranges::views::transform(
                                [](const HP::Particle &p) { return p.id(); }),
                        ", ");
                }),
                "} {")));
        REQUIRE(c2.size() == 4);
        CHECK(clusterListContains(c2, HP::ParticleSet{m4}));
        CHECK(clusterListContains(c2, HP::ParticleSet{m, mu, muu}));
        CHECK(clusterListContains(c2, HP::ParticleSet{muu, m2, m3}));
        CHECK(clusterListContains(c2, HP::ParticleSet{m2, mu, muu}));
    }
}

namespace {
struct DummyParticle : public HP::Particle {
    double rate_;
    DummyParticle(double mass, double rate)
        : Particle("dummy", HP::CP::undefined, HP::ECharge{}), rate_{rate} {
        setMass(mass);
    }
    double channelRateOr(HP::Collider, HP::Production, HP::Decay,
                         double) const noexcept override {
        return rate_;
    }
    std::unique_ptr<Particle> clone() const override {
        return std::make_unique<DummyParticle>(*this);
    }
};
} // namespace

TEST_CASE("cluster properties") {

    SECTION("rate weighted -- 1 cluster") {
        auto rateFunc = [](const auto &p) {
            return p.channelRate(HP::Collider::LHC8, HP::Production::bbH,
                                 HP::Decay::bb);
        };
        const auto empty = std::vector<DummyParticle>{};
        CHECK(HPC::combinedRate(rateFunc, empty) == 0.);
        CHECK(std::get<0>(HPC::rateWeightedMasses(rateFunc, empty)).mass == 0.);
        CHECK(
            std::get<0>(HPC::rateWeightedMasses(rateFunc, empty)).uncertainty ==
            0.);

        auto testCluster =
            std::vector{DummyParticle{100., 1.}, DummyParticle{140., 1.},
                        DummyParticle{200, 6.}};
        CHECK(HPC::combinedRate(rateFunc, testCluster | ranges::views::take(
                                                            1)) == Approx(1.));
        CHECK(std::get<0>(HPC::rateWeightedMasses(
                              rateFunc, testCluster | ranges::views::take(1)))
                  .mass == 100);

        CHECK(HPC::combinedRate(rateFunc, testCluster | ranges::views::take(
                                                            2)) == Approx(2.));
        CHECK(std::get<0>(HPC::rateWeightedMasses(
                              rateFunc, testCluster | ranges::views::take(2)))
                  .mass == 120);

        CHECK(HPC::combinedRate(rateFunc, testCluster) == Approx(8.));
        CHECK(
            std::get<0>(HPC::rateWeightedMasses(rateFunc, testCluster)).mass ==
            Approx(180));
    }

    SECTION("rate weighted -- 2 clusters") {
        auto rateFunc = [](const HP::Particle &p1, const HP::Particle &p2) {
            return p1.channelRate(HP::Collider::LHC8, HP::Production::bbH,
                                  HP::Decay::bb) *
                   p2.channelRate(HP::Collider::LHC8, HP::Production::bbH,
                                  HP::Decay::bb);
        };
        const auto empty = std::vector<DummyParticle>{};
        CHECK(HPC::combinedRate(rateFunc, empty, empty) == 0.);
        const auto [em1, em2] = HPC::rateWeightedMasses(rateFunc, empty, empty);
        CHECK(em1 == HP::UncertainMass{});
        CHECK(em2 == HP::UncertainMass{});

        auto c1 = std::vector{DummyParticle{400., 1.}, DummyParticle{340., 1.},
                              DummyParticle{500, 3.}};
        auto c2 = std::vector{DummyParticle{100., 0.3},
                              DummyParticle{140., 0.2}, DummyParticle{300, 1.}};

        SECTION("single particles") {
            CHECK(HPC::combinedRate(rateFunc, c1 | ranges::views::take(1),
                                    c2 | ranges::views::take(1)) ==
                  Approx(0.3));
            auto [m1, m2] =
                HPC::rateWeightedMasses(rateFunc, c1 | ranges::views::take(1),
                                        c2 | ranges::views::take(1));
            CHECK(m1.mass == Approx(400.));
            CHECK(m2.mass == Approx(100.));
        }

        SECTION("two particles") {
            CHECK(HPC::combinedRate(rateFunc, c1 | ranges::views::take(2),
                                    c2 | ranges::views::take(2)) == Approx(1.));
            auto [m1, m2] =
                HPC::rateWeightedMasses(rateFunc, c1 | ranges::views::take(2),
                                        c2 | ranges::views::take(2));
            CHECK(m1.mass == Approx(370.));
            CHECK(m2.mass == Approx(116.));
        }

        SECTION("three particles") {
            CHECK(HPC::combinedRate(rateFunc, c1, c2) == Approx(7.5));
            auto [m1, m2] = HPC::rateWeightedMasses(rateFunc, c1, c2);
            CHECK(m1.mass ==
                  Approx((400. * 1. + 340. * 1. + 500. * 3.) * 1.5 / 7.5));
            CHECK(m2.mass ==
                  Approx((100. * 0.3 + 140. * 0.2 + 300. * 1.) * 5. / 7.5));
        }
    }
}
