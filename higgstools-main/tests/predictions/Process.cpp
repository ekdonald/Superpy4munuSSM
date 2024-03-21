#include "Higgs/predictions/Process.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "utilities/Json.hpp"
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <initializer_list>
#include <map>
#include <memory>
#include <range/v3/algorithm/count.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/view.hpp>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
namespace HP = Higgs::predictions;

TEST_CASE("ChannelProcess operations") {
    SECTION("default construction") {
        auto ch = HP::ChannelProcess{};
        REQUIRE(ch.collider() == HP::Collider::LHC13);
        auto h0 = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                          HP::Decay::gamgam, 1.);
        REQUIRE(ch(h0) == 0.);
    }

    const auto ch =
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::ggH, HP::Decay::gamgam},
                            {HP::Production::bbH, HP::Decay::mumu}}};

    SECTION("description") {
        REQUIRE_THAT(ch.to_string(), ContainsSubstring("ggH") &&
                                         ContainsSubstring("bbH") &&
                                         ContainsSubstring("gamgam") &&
                                         ContainsSubstring("mumu") &&
                                         ContainsSubstring("LHC13"));
    }

    SECTION("parsing") {
        const auto j = nlohmann::json::parse(" {"
                                             "     \"channels\": ["
                                             "         [\"ggH\", \"gamgam\"],"
                                             "         [\"bbH\", \"mumu\"]"
                                             "     ]"
                                             " }");
        REQUIRE_NOTHROW(HP::readChannelProcess(j, HP::Collider::LHC8));
        REQUIRE(HP::readChannelProcess(j, HP::Collider::LHC13) == ch);
    }

    SECTION("application to particles") {
        auto h0 = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                          HP::Decay::gamgam, 1.);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::bbH,
                          HP::Decay::mumu, 2.);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                          HP::Decay::mumu, 1000.);
        CHECK(ch(h0) == Approx(3.));
        CHECK(ch(h0, {4., 0.5}) == Approx(5.));
        CHECK(ch(h0, {2., 3.}) == Approx(8.));
        CHECK(ch(h0) == ch(h0, {}));
        CHECK(ch(h0, {3.}) == Approx(5.));
        CHECK(ch(h0, {5., 4., 3., 2., 1.}) == Approx(13.));
    }
    SECTION("merging of processes -- valid") {
        REQUIRE(ch + HP::ChannelProcess{} == ch);
        REQUIRE(HP::ChannelProcess{} + ch == ch);

        auto ch2 = HP::ChannelProcess{HP::Collider::LHC13,
                                      {{HP::Production::HZ, HP::Decay::gamgam},
                                       {HP::Production::HW, HP::Decay::mumu}}};

        auto h0 = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                          HP::Decay::gamgam, 1.);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::bbH,
                          HP::Decay::mumu, 2.);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::HZ,
                          HP::Decay::gamgam, 3.);
        h0.setChannelRate(HP::Collider::LHC13, HP::Production::HW,
                          HP::Decay::mumu, 4.);
        CHECK(ch(h0) == Approx(3.));
        CHECK(ch2(h0) == Approx(3 + 4));

        auto ch3 = ch + ch2;
        CHECK(ch3(h0) == Approx(3 + 3 + 4));
        ch2 += ch3;
        CHECK(ch2(h0) == Approx(3 + 2 * 3 + 2 * 4));
    }
    SECTION("merging of processes -- invalid") {
        const auto ch1 =
            HP::ChannelProcess{HP::Collider::LHC8,
                               {{HP::Production::ggH, HP::Decay::gamgam},
                                {HP::Production::bbH, HP::Decay::mumu}}};
        REQUIRE(ch + ch1 == ch);
        REQUIRE(ch1 + ch == ch1);
    }
    SECTION("subprocesses") {
        const auto subProcs = ch.subprocesses();
        INFO(ch.to_string());
        REQUIRE(subProcs.size() == ch.size());
        REQUIRE(ch == ranges::accumulate(subProcs, HP::ChannelProcess{}));
    }
}

TEST_CASE("ChainDecayProcess operations") {
    SECTION("default construction") {
        auto ch = HP::ChainDecayProcess{};
        REQUIRE(ch.collider() == HP::Collider::LHC13);
        REQUIRE(ch.chain() == HP::ChainDecay{});
        REQUIRE(ch.productionSize() == 0);
        REQUIRE(ch.decaySize() == 0);

        auto hh = HP::BsmParticle("hh", HP::ECharge::neutral, HP::CP::undefined);
        hh.setTotalWidth(1.);
        hh.setCxn(HP::Collider::LHC13, HP::Production::ggH, 1.);
        hh.setBr(HP::ChainDecay::Z, "hl", 0.4);
        auto hl = HP::BsmParticle("hl", HP::ECharge::neutral, HP::CP::undefined);
        hl.setTotalWidth(1.);
        hl.setBr(HP::Decay::gamgam, 0.5);
        REQUIRE(ch(hh, hl) == 0.);
    }

    SECTION("parsing") {
        const auto j = nlohmann::json::parse(" {"
                                             "     \"production\": ["
                                             "         \"ggH\", \"bbH\""
                                             "     ],"
                                             "     \"chain\": \"Z\","
                                             "     \"decay\": ["
                                             "         \"gamgam\""
                                             "     ]"
                                             " }");
        REQUIRE_NOTHROW(HP::readChainDecayProcess(j, HP::Collider::LHC8));
        auto proc = HP::readChainDecayProcess(j, HP::Collider::LHC13);
        REQUIRE(proc.productionSize() == 2);
        REQUIRE(proc.decaySize() == 1);
        CHECK(proc.collider() == HP::Collider::LHC13);
        CHECK(proc.chain() == HP::ChainDecay::Z);
    }

    SECTION("application to particles") {
        auto procZH = HP::ChainDecayProcess(
            HP::ChainDecay::Z, {HP::Production::ggH, HP::Production::bbH},
            {HP::Decay::gamgam, HP::Decay::mumu}, HP::Collider::LHC13);

        auto hh = HP::BsmParticle("hh", HP::ECharge::neutral, HP::CP::undefined);
        hh.setTotalWidth(1.);
        hh.setCxn(HP::Collider::LHC13, HP::Production::ggH, 0.9);
        hh.setCxn(HP::Collider::LHC13, HP::Production::bbH, 0.1);
        hh.setBr(HP::ChainDecay::Z, "hl", 0.4);
        hh.setBr(HP::ChainDecay::W, "H+", 0.6);

        auto hl = HP::BsmParticle("hl", HP::ECharge::neutral, HP::CP::undefined);
        hl.setTotalWidth(1.);
        hl.setCxn(HP::Collider::LHC13, HP::Production::bbH, 2.);
        hl.setBr(HP::Decay::gamgam, 0.5);
        hl.setBr(HP::ChainDecay::W, "H+", 0.1);

        auto hp = HP::BsmParticle("H+", HP::ECharge::single, HP::CP::undefined);
        hp.setTotalWidth(1.);
        hp.setBr(HP::Decay::taunu, 0.25);

        CHECK(procZH(hh, hl) == Approx(1. * 0.4 * 0.5));
        CHECK(procZH(hl, hh) == 0.);

        auto procWH = HP::ChainDecayProcess(
            HP::ChainDecay::W, {HP::Production::ggH, HP::Production::bbH},
            {HP::Decay::taunu}, HP::Collider::LHC13);
        CHECK(procWH(hh, hp) == Approx(1. * 0.6 * 0.25));

        SECTION("production weights") {
            CHECK(procZH(hh, hl, {3, 10, 1e6, 1e6}) ==
                  Approx((3 * 0.9 + 1) * 0.4 * 0.5));
            CHECK(procZH(hh, hl, {}) == procZH(hh, hl));
        }
    }

    SECTION("output") {
        auto proc = HP::ChainDecayProcess(
            HP::ChainDecay::Z, {HP::Production::ggH, HP::Production::bbH},
            {HP::Decay::gamgam, HP::Decay::mumu}, HP::Collider::LHC13);
        REQUIRE(proc.to_string() ==
                "LHC13 [ggH, bbH]->X1->Z(X2->[gamgam, mumu])");
        auto proc2 = HP::ChainDecayProcess(
            HP::ChainDecay::W, {HP::Production::Hpmtb},
            {HP::Decay::gamgam, HP::Decay::mumu}, HP::Collider::LHC13);
        REQUIRE(proc2.to_string() ==
                "LHC13 [Hpmtb]->X1->W(X2->[gamgam, mumu])");
    }
}

TEST_CASE("PairDecayProcess") {
    auto hh = HP::BsmParticle("hh", HP::ECharge::neutral, HP::CP::undefined);
    hh.setTotalWidth(1.);
    hh.setCxn(HP::Collider::LHC13, HP::Production::ggH, 1.);
    hh.setBr("h1", "h1", 0.4);
    hh.setBr("h1", "h2", 0.3);

    auto h1 = HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined);
    h1.setTotalWidth(1.);
    h1.setBr(HP::Decay::bb, 0.5);
    h1.setBr(HP::Decay::gamgam, 0.1);
    h1.setBr(HP::Decay::tautau, 0.4);

    auto h2 = HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined);
    h2.setTotalWidth(1.);
    h2.setBr(HP::Decay::bb, 0.4);
    h2.setBr(HP::Decay::tautau, 0.2);
    h2.setBr(HP::Decay::gamgam, 0.4);

    SECTION("default construction") {
        auto proc = HP::PairDecayProcess{};
        REQUIRE(proc.collider() == HP::Collider::LHC13);

        REQUIRE(proc(hh, h1, h1) == 0);
        REQUIRE(proc(hh, h1, h2) == 0);
        REQUIRE(proc(hh, h2, h2) == 0);
    }

    const auto proc = HP::PairDecayProcess{{HP::Production::ggH},
                                           {HP::Decay::bb, HP::Decay::tautau},
                                           {HP::Decay::bb, HP::Decay::gamgam},
                                           HP::Collider::LHC13};

    const auto procSame = HP::PairDecayProcess{{HP::Production::ggH},
                                               {HP::Decay::bb},
                                               {HP::Decay::bb},
                                               HP::Collider::LHC13};

    const auto proc1 = HP::PairDecayProcess{{HP::Production::ggH},
                                            {HP::Decay::bb},
                                            {HP::Decay::gamgam},
                                            HP::Collider::LHC13};

    const auto procAll = HP::PairDecayProcess{
        {HP::Production::ggH},
        {HP::Decay::bb, HP::Decay::tautau, HP::Decay::gamgam},
        {HP::Decay::bb, HP::Decay::tautau, HP::Decay::gamgam},
        HP::Collider::LHC13};

    SECTION("application to particles") {
        SECTION("process with identical decays, no symmetry factors") {
            CHECK(procSame(hh, h1, h1) ==
                  Approx(hh.cxn(procSame.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h1.id()) *
                         std::pow(h1.br(HP::Decay::bb), 2)));
            CHECK(procSame(hh, h1, h2) ==
                  Approx(hh.cxn(procSame.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h2.id()) * h1.br(HP::Decay::bb) *
                         h2.br(HP::Decay::bb)));
        }

        SECTION("process with distinct decays, symmetry factors for identical "
                "daughter particles") {
            CHECK(proc1(hh, h1, h1) ==
                  Approx(hh.cxn(proc1.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h1.id()) * 2 * h1.br(HP::Decay::bb) *
                         h1.br(HP::Decay::gamgam)));
            CHECK(proc1(hh, h1, h2) ==
                  Approx(hh.cxn(proc1.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h2.id()) * h1.br(HP::Decay::bb) *
                         h2.br(HP::Decay::gamgam)));
        }

        SECTION("process covering all decay modes, BR(h_i h_j -> anything)=1") {
            CHECK(procAll(hh, h1, h1) ==
                  Approx(hh.cxn(procAll.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h1.id())));
            CHECK(procAll(hh, h1, h2) ==
                  Approx(hh.cxn(procAll.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h2.id())));
        }

        SECTION("process with partially distinct decay modes") {
            CHECK(proc(hh, h1, h1) ==
                  Approx(hh.cxn(proc.collider(), HP::Production::ggH) *
                         hh.br(h1.id(), h1.id()) *
                         (std::pow(h1.br(HP::Decay::bb), 2) +
                          2 * h1.br(HP::Decay::bb) * h1.br(HP::Decay::gamgam) +
                          2 * h1.br(HP::Decay::bb) * h1.br(HP::Decay::tautau) +
                          2 * h1.br(HP::Decay::gamgam) *
                              h1.br(HP::Decay::tautau))));
            CHECK(
                proc(hh, h1, h2) ==
                Approx(hh.cxn(proc.collider(), HP::Production::ggH) *
                       hh.br(h1.id(), h2.id()) *
                       (h1.br(HP::Decay::bb) * h2.br(HP::Decay::bb) +
                        h1.br(HP::Decay::bb) * h2.br(HP::Decay::gamgam) +
                        h1.br(HP::Decay::tautau) * h2.br(HP::Decay::bb) +
                        h1.br(HP::Decay::tautau) * h2.br(HP::Decay::gamgam))));
        }
    }

    SECTION("parsing") {
        const auto j = nlohmann::json::parse(R"({
            "production": ["ggH"],
            "firstDecay":["bb", "tautau"],
            "secondDecay":["gamgam","bb"]
        })");
        REQUIRE_NOTHROW(HP::readPairDecayProcess(j, HP::Collider::LHC8));
        auto jproc = HP::readPairDecayProcess(j, HP::Collider::LHC13);
        CHECK(jproc.collider() == HP::Collider::LHC13);
        CHECK(jproc(hh, h1, h1) == proc(hh, h1, h1));
    }
}

TEST_CASE("PairProductionProcess") {
    auto pred = Higgs::Predictions();
    auto &h1 = pred.addParticle(HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
    auto &h2 = pred.addParticle(HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
    pred.setBsmPairCxn(HP::Collider::LHC13, "h1", "h1", 10);
    pred.setBsmPairCxn(HP::Collider::LHC13, "h1", "h2", 20);
    pred.setBsmPairCxn(HP::Collider::LHC13, "h2", "h2", 30);

    h1.setTotalWidth(1.);
    h1.setBr(HP::Decay::bb, 0.5);
    h1.setBr(HP::Decay::gamgam, 0.1);
    h1.setBr(HP::Decay::tautau, 0.4);

    h2.setTotalWidth(1.);
    h2.setBr(HP::Decay::bb, 0.4);
    h2.setBr(HP::Decay::tautau, 0.2);
    h2.setBr(HP::Decay::gamgam, 0.4);

    SECTION("default construction") {
        auto proc = HP::PairProductionProcess{};
        REQUIRE(proc.collider() == HP::Collider::LHC13);

        REQUIRE(proc(pred, h1, h1) == 0);
        REQUIRE(proc(pred, h1, h2) == 0);
        REQUIRE(proc(pred, h2, h2) == 0);
    }

    const auto proc =
        HP::PairProductionProcess{{HP::Decay::bb, HP::Decay::tautau},
                                  {HP::Decay::bb, HP::Decay::gamgam},
                                  HP::Collider::LHC13};

    const auto procSame = HP::PairProductionProcess{
        {HP::Decay::bb}, {HP::Decay::bb}, HP::Collider::LHC13};

    const auto proc1 = HP::PairProductionProcess{
        {HP::Decay::bb}, {HP::Decay::gamgam}, HP::Collider::LHC13};

    const auto procAll = HP::PairProductionProcess{
        {HP::Decay::bb, HP::Decay::tautau, HP::Decay::gamgam},
        {HP::Decay::bb, HP::Decay::tautau, HP::Decay::gamgam},
        HP::Collider::LHC13};

    SECTION("application to particles") {
        SECTION("process with identical decays, no symmetry factors") {
            CHECK(
                procSame(pred, h1, h1) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h1.id(), h1.id()) *
                       std::pow(h1.br(HP::Decay::bb), 2)));
            CHECK(
                procSame(pred, h1, h2) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h1.id(), h2.id()) *
                       h1.br(HP::Decay::bb) * h2.br(HP::Decay::bb)));
            CHECK(
                procSame(pred, h2, h2) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h2.id(), h2.id()) *
                       std::pow(h2.br(HP::Decay::bb), 2)));
        }

        SECTION("process with distinct decays, symmetry factors for identical "
                "initial particles") {
            CHECK(
                proc1(pred, h1, h1) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h1.id(), h1.id()) *
                       2 * h1.br(HP::Decay::bb) * h1.br(HP::Decay::gamgam)));
            CHECK(
                proc1(pred, h1, h2) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h1.id(), h2.id()) *
                       h1.br(HP::Decay::bb) * h2.br(HP::Decay::gamgam)));
            CHECK(
                proc1(pred, h2, h2) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h2.id(), h2.id()) *
                       2 * h2.br(HP::Decay::bb) * h2.br(HP::Decay::gamgam)));
        }

        SECTION("process covering all decay modes, BR(h_i h_j -> anything)=1") {
            CHECK(procAll(pred, h1, h1) ==
                  Approx(h1.cxn(procAll.collider(), HP::Production::pair)));
            CHECK(
                procAll(pred, h1, h2) ==
                Approx(pred.bsmPairCxn(procAll.collider(), h1.id(), h2.id())));
            CHECK(procAll(pred, h2, h2) ==
                  Approx(h2.cxn(procAll.collider(), HP::Production::pair)));
        }

        SECTION("process with partially distinct decay modes") {
            CHECK(
                proc(pred, h1, h1) ==
                Approx(
                    pred.bsmPairCxn(procSame.collider(), h1.id(), h1.id()) *
                    (std::pow(h1.br(HP::Decay::bb), 2) +
                     2 * h1.br(HP::Decay::bb) * h1.br(HP::Decay::gamgam) +
                     2 * h1.br(HP::Decay::bb) * h1.br(HP::Decay::tautau) +
                     2 * h1.br(HP::Decay::gamgam) * h1.br(HP::Decay::tautau))));
            CHECK(
                proc(pred, h1, h2) ==
                Approx(pred.bsmPairCxn(procSame.collider(), h1.id(), h2.id()) *
                       (h1.br(HP::Decay::bb) * h2.br(HP::Decay::bb) +
                        h1.br(HP::Decay::bb) * h2.br(HP::Decay::gamgam) +
                        h1.br(HP::Decay::tautau) * h2.br(HP::Decay::bb) +
                        h1.br(HP::Decay::tautau) * h2.br(HP::Decay::gamgam))));
            CHECK(
                proc(pred, h2, h2) ==
                Approx(
                    pred.bsmPairCxn(procSame.collider(), h2.id(), h2.id()) *
                    (std::pow(h2.br(HP::Decay::bb), 2) +
                     2 * h2.br(HP::Decay::bb) * h2.br(HP::Decay::gamgam) +
                     2 * h2.br(HP::Decay::bb) * h2.br(HP::Decay::tautau) +
                     2 * h2.br(HP::Decay::gamgam) * h2.br(HP::Decay::tautau))));
        }
    }

    SECTION("parsing") {
        const auto j = nlohmann::json::parse(R"({
            "firstDecay":["bb", "tautau"],
            "secondDecay":["gamgam","bb"]
        })");
        REQUIRE_NOTHROW(HP::readPairProductionProcess(j, HP::Collider::LHC8));
        auto jproc = HP::readPairProductionProcess(j, HP::Collider::LHC13);
        CHECK(jproc.collider() == HP::Collider::LHC13);
        CHECK(jproc(pred, h1, h1) == proc(pred, h1, h1));
    }
}

namespace {
struct Identifiable : public HP::Particle {
    Identifiable(std::string &&id)
        : Particle{std::move(id), HP::CP::undefined, HP::ECharge{}} {}
    std::unique_ptr<Particle> clone() const override {
        return std::make_unique<Identifiable>(*this);
    }
};
} // namespace

TEST_CASE("contributing particles") {
    auto h2 = Identifiable{"H2"};
    auto h3 = Identifiable{"H3"};
    auto h4 = Identifiable{"H4"};
    SECTION("single cluster -- ChannelProcess") {
        REQUIRE(HP::ChannelProcess::contributingParticles(HP::ParticleSet{}) ==
                std::vector<std::string>{});

        REQUIRE(HP::ChannelProcess::contributingParticles(
                    HP::ParticleSet{h4}) == std::vector<std::string>{"H4"});
        REQUIRE(HP::ChannelProcess::contributingParticles(HP::ParticleSet{
                    h4, h2}) == std::vector<std::string>{"H2", "H4"});
        REQUIRE(HP::ChannelProcess::contributingParticles(HP::ParticleSet{
                    h3, h2, h4}) == std::vector<std::string>{"H2", "H3", "H4"});
    }

    SECTION("two clusters -- ChainDecayProcess") {
        REQUIRE(HP::ChainDecayProcess::contributingParticles(
                    HP::ParticleSet{}, HP::ParticleSet{}) ==
                std::vector<std::string>{">"});

        REQUIRE(HP::ChainDecayProcess::contributingParticles(
                    HP::ParticleSet{h4}, HP::ParticleSet{h3, h2, h4}) ==
                std::vector<std::string>{"H4", ">", "H2", "H3", "H4"});
    }
    SECTION("three clusters -- PairDecayProcess") {
        REQUIRE(HP::PairDecayProcess::contributingParticles(
                    HP::ParticleSet{}, HP::ParticleSet{}, HP::ParticleSet{}) ==
                std::vector<std::string>{">", "+"});

        REQUIRE(HP::PairDecayProcess::contributingParticles(
                    HP::ParticleSet{h4, h2}, HP::ParticleSet{h3, h2, h4},
                    HP::ParticleSet{h4}) ==
                std::vector<std::string>{"H2", "H4", ">", "H2", "H3", "H4", "+",
                                         "H4"});
    }
}
