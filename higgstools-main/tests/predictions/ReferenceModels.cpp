#include "Higgs/predictions/ReferenceModels.hpp"
#include "Higgs/predictions/Basics.hpp"
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <magic_enum.hpp>
#include <memory>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/remove.hpp>
#include <range/v3/view/transform.hpp>
#include <utility>

namespace HP = Higgs::predictions;
using Catch::Approx;

namespace {
constexpr auto prodChannels = magic_enum::enum_values<HP::Production>();
constexpr auto colliders = magic_enum::enum_values<HP::Collider>();
constexpr auto decays = magic_enum::enum_values<HP::Decay>();
} // namespace

TEMPLATE_TEST_CASE("Reference Higgs cross sections", "", HP::SMHiggs,
                   HP::SMHiggsEW, HP::SMHiggsInterp) {

    constexpr auto isHtCxn = [](HP::Production p) {
        return p == HP::Production::tchanHt || p == HP::Production::schanHt ||
               p == HP::Production::Ht || p == HP::Production::HtW;
    };

    constexpr auto isFV = [](HP::Production p) {
        return p == HP::Production::ucHgam || p == HP::Production::dsHgam ||
               p == HP::Production::dbHgam || p == HP::Production::sbHgam;
    };

    constexpr auto isHgamCxn = [](HP::Production p) {
        return p == HP::Production::uuHgam || p == HP::Production::ddHgam ||
               p == HP::Production::ccHgam || p == HP::Production::ssHgam ||
               p == HP::Production::bbHgam;
    };

    const auto unavailableCxn = [&isHtCxn, &isFV, &isHgamCxn](auto coll,
                                                              auto &p) {
        bool result = false;
        switch (coll) {
        case HP::Collider::LHC8:
            result = p == HP::Production::qqH || p == HP::Production::brtHu ||
                     p == HP::Production::brtHc || isHgamCxn(p) ||
                     !HP::validProductionFor(p, HP::ECharge::neutral) ||
                     !validProductionAt(p, HP::ColliderType::pp);
            if constexpr (std::is_same_v<TestType, HP::SMHiggs>) {
                result |= isHtCxn(p);
            }
            break;
        case HP::Collider::LHC13:
            result = p == HP::Production::brtHu || p == HP::Production::brtHc ||
                     !HP::validProductionFor(p, HP::ECharge::neutral) ||
                     !validProductionAt(p, HP::ColliderType::pp);
            break;
        case HP::Collider::LEP:
            result = !validProductionAt(p, HP::ColliderType::ee);
        }
        result |= p == HP::Production::pair; // not implemented (yet)
        result |= isFV(p);                   // flavor violation
        return result;
    };

    SECTION("extremely light") {
        INFO("for too light masses the cxns should clamp");
        auto v = TestType(1e-4);
        auto vv = TestType(1e-8);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                INFO(magic_enum::enum_name(coll)
                     << " " << magic_enum::enum_name(p));
                CHECK(v.cxn(coll, p) == Approx(vv.cxn(coll, p)));
            }
        }
    }
    SECTION("125") {
        INFO("all h125 reference cxns need to be positive");
        auto h125 = TestType(125);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                INFO(magic_enum::enum_name(coll)
                     << " " << magic_enum::enum_name(p));
                if (p == HP::Production::none || unavailableCxn(coll, p)) {
                    REQUIRE(h125.cxn(coll, p) == 0.);
                } else {
                    CHECK(h125.cxn(coll, p) > 0);
                    CHECK(h125.cxn(coll, p) ==
                          Approx(h125.channelRate(coll, p, HP::Decay::none)));
                }
            }
        }
    }
    SECTION("medium mass") {
        INFO("at intermediate masses, all reference cxns need to be positive");
        auto med = TestType(300);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                if (p == HP::Production::none || unavailableCxn(coll, p) ||
                    // the LHC8 tH cxns in the SMHiggsEW are only available
                    // around 125GeV
                    (coll == HP::Collider::LHC8 && isHtCxn(p))) {
                    REQUIRE(med.cxn(coll, p) == 0.);
                } else {
                    INFO(magic_enum::enum_name(coll)
                         << " " << magic_enum::enum_name(p));
                    CHECK(med.cxn(coll, p) > 0);
                }
            }
        }
    }
    SECTION("extremely heavy") {
        INFO("for too high masses, the reference cxns should be 0");
        auto h = TestType(1e6);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                INFO(magic_enum::enum_name(coll)
                     << " " << magic_enum::enum_name(p));
                if (coll == HP::Collider::LEP && !unavailableCxn(coll, p) &&
                    p != HP::Production::none) {
                    INFO("except for LEP, where everything is normalized and "
                         "remains 1.");
                    CHECK(h.cxn(coll, p) == 1.);
                } else {
                    CHECK(h.cxn(coll, p) == 0.);
                }
            }
        }
    }
}

TEMPLATE_TEST_CASE("SM-like BRs", "", HP::SMHiggs, HP::SMHiggsEW,
                   HP::SMHiggsInterp) {
    for (auto h : {TestType(1e-4), TestType(20), TestType(125), TestType(500),
                   TestType(10000)}) {

        for (auto d : decays) {
            if (d == HP::Decay::none) {
                CHECK(h.br(d) == 0);
            } else {
                CHECK(h.br(d) ==
                      Approx(h.channelRate(HP::Collider::LHC13,
                                           HP::Production::none, d)));
            }
        }
        auto toBr = [&h](auto d) { return h.br(d); };
        auto summedBr =
            ranges::accumulate(decays | ranges::views::remove(HP::Decay::inv) |
                                   ranges::views::transform(toBr),
                               0.);
        if (summedBr > 0) {
            CHECK(summedBr == Approx(1));
        }
        CHECK(
            h.br(HP::Decay::inv) ==
            Approx(h.br(HP::Decay::ZZ) * std::pow(HP::constants::b_Z_inv, 2)));
    }

    SECTION("OOB behaviour") {
        for (auto d : decays) {
            REQUIRE(TestType(1e-6).br(d) == 0); // zero for too low masses
            REQUIRE(TestType(-10).br(d) == 0);  // zero for negative masses
            REQUIRE(TestType(1e6).br(d) == 0);  // zero for too large masses
        }
        REQUIRE(TestType(1e-6).totalWidth() ==
                TestType(1e-3).totalWidth());     // clamp for too low masses
        REQUIRE(TestType(1e-6).totalWidth() > 0); // non-zero for low masses
        REQUIRE(TestType(-10).totalWidth() == 0); // zero for negative masses
        REQUIRE(TestType(1e6).totalWidth() ==
                TestType(1e4).totalWidth());     // clamp for too large masses
        REQUIRE(TestType(1e6).totalWidth() > 0); // non-zero for huge masses
    }
}

TEST_CASE("SMHiggsInterp") {
    SECTION("low mass") {
        double low_mass = 140;
        auto hSMHiggsEW = HP::SMHiggsEW(low_mass);
        auto hSMHiggsInterp = HP::SMHiggsInterp(low_mass);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                INFO(magic_enum::enum_name(coll)
                     << " " << magic_enum::enum_name(p));
                CHECK(hSMHiggsEW.cxn(coll, p) ==
                      Approx(hSMHiggsInterp.cxn(coll, p)));
            }
        }
    }
    SECTION("high mass") {
        double high_mass = 400;
        auto hSMHiggs = HP::SMHiggs(high_mass);
        auto hSMHiggsInterp = HP::SMHiggsInterp(high_mass);
        for (auto coll : colliders) {
            for (auto p : prodChannels) {
                INFO(magic_enum::enum_name(coll)
                     << " " << magic_enum::enum_name(p));
                CHECK(hSMHiggs.cxn(coll, p) ==
                      Approx(hSMHiggsInterp.cxn(coll, p)));
            }
        }
    }
}

TEMPLATE_TEST_CASE("Reference couplings", "", HP::SMHiggs, HP::SMHiggsEW, HP::SMHiggsInterp) {
    auto h = TestType(10);
    for (auto c : magic_enum::enum_values<HP::Coupling>()) {
        if (c == HP::Coupling::effCPeTopYuk || c == HP::Coupling::effVV) {
            CHECK(h.coupling(c) == 1.);
        } else {
            CHECK(h.coupling(c) == 0.);
        }
    }
}

TEMPLATE_TEST_CASE("Clone reference particles", "", HP::SMHiggs,
                   HP::SMHiggsEW, HP::SMHiggsInterp) {
    auto p = TestType{100.};
    HP::Particle *pc = &p;
    auto pcc = pc->clone();
    REQUIRE(pcc->totalWidth() == pc->totalWidth());
}

TEST_CASE("SMHiggs cxns match YR4") {
    SECTION("ggF LHC8") {
        static constexpr auto vals = std::array<std::pair<double, double>, 39>{
            {{1.00000000e+01, 4.61400000e+03}, {8.86842105e+01, 3.63494737e+01},
             {1.67368421e+02, 1.03054737e+01}, {2.46052632e+02, 4.77726316e+00},
             {3.24736842e+02, 3.19973684e+00}, {4.03421053e+02, 2.94292105e+00},
             {4.82105263e+02, 1.53336842e+00}, {5.60789474e+02, 7.58971053e-01},
             {6.39473684e+02, 3.76942105e-01}, {7.18157895e+02, 1.93717368e-01},
             {7.96842105e+02, 1.01515895e-01}, {8.75526316e+02, 5.65161579e-02},
             {9.54210526e+02, 3.15683158e-02}, {1.03289474e+03, 1.85276316e-02},
             {1.11157895e+03, 1.09759368e-02}, {1.19026316e+03, 6.66697368e-03},
             {1.26894737e+03, 4.13794737e-03}, {1.34763158e+03, 2.58083158e-03},
             {1.42631579e+03, 1.66321053e-03}, {1.50500000e+03, 1.06899000e-03},
             {1.58368421e+03, 7.03598947e-04}, {1.66236842e+03, 4.64794737e-04},
             {1.74105263e+03, 3.10382105e-04}, {1.81973684e+03, 2.09805263e-04},
             {1.89842105e+03, 1.41628421e-04}, {1.97710526e+03, 9.75876316e-05},
             {2.05578947e+03, 6.69084211e-05}, {2.13447368e+03, 4.65599474e-05},
             {2.21315789e+03, 3.23557895e-05}, {2.29184211e+03, 2.26257368e-05},
             {2.37052632e+03, 1.59454737e-05}, {2.44921053e+03, 1.11737368e-05},
             {2.52789474e+03, 7.96191579e-06}, {2.60657895e+03, 5.62547368e-06},
             {2.68526316e+03, 4.02180000e-06}, {2.76394737e+03, 2.86816316e-06},
             {2.84263158e+03, 2.04914737e-06}, {2.92131579e+03, 1.47354737e-06},
             {3.00000000e+03, 1.04800000e-06}}};
        for (auto [x, y] : vals) {
            REQUIRE(HP::SMHiggs(x).cxn(HP::Collider::LHC8,
                                       HP::Production::ggH) == Approx(y));
        }
    }
    SECTION("ggF LHC13 NNLO+NNLL") {
        static constexpr auto vals = std::array<std::pair<double, double>, 39>{
            {{1.00000000e+01, 6.99600000e+03}, {8.86842105e+01, 7.71621053e+01},
             {1.67368421e+02, 2.49215789e+01}, {2.46052632e+02, 1.28076316e+01},
             {3.24736842e+02, 9.37163158e+00}, {4.03421053e+02, 9.33126316e+00},
             {4.82105263e+02, 5.25321053e+00}, {5.60789474e+02, 2.79177895e+00},
             {6.39473684e+02, 1.48968421e+00}, {7.18157895e+02, 8.20472105e-01},
             {7.96842105e+02, 4.61144211e-01}, {8.75526316e+02, 2.73858947e-01},
             {9.54210526e+02, 1.63777895e-01}, {1.03289474e+03, 1.02438158e-01},
             {1.11157895e+03, 6.48010526e-02}, {1.19026316e+03, 4.19850000e-02},
             {1.26894737e+03, 2.77814737e-02}, {1.34763158e+03, 1.85267895e-02},
             {1.42631579e+03, 1.27047368e-02}, {1.50500000e+03, 8.73060000e-03},
             {1.58368421e+03, 6.13034737e-03}, {1.66236842e+03, 4.32732632e-03},
             {1.74105263e+03, 3.08972632e-03}, {1.81973684e+03, 2.23234211e-03},
             {1.89842105e+03, 1.61727368e-03}, {1.97710526e+03, 1.19068947e-03},
             {2.05578947e+03, 8.77491579e-04}, {2.13447368e+03, 6.53867895e-04},
             {2.21315789e+03, 4.88668421e-04}, {2.29184211e+03, 3.68288421e-04},
             {2.37052632e+03, 2.79135789e-04}, {2.44921053e+03, 2.11630000e-04},
             {2.52789474e+03, 1.62234737e-04}, {2.60657895e+03, 1.24634211e-04},
             {2.68526316e+03, 9.60702105e-05}, {2.76394737e+03, 7.44742105e-05},
             {2.84263158e+03, 5.77968421e-05}, {2.92131579e+03, 4.51185789e-05},
             {3.00000000e+03, 3.50200000e-05}}};
        for (auto [x, y] : vals) {
            REQUIRE(HP::SMHiggs(x).cxn(HP::Collider::LHC13,
                                       HP::Production::ggH) == Approx(y));
        }
    }
    SECTION("bbH LHC8") {
        static constexpr auto vals = std::array<std::pair<double, double>, 39>{
            {{1.00000000e+01, 6.36600000e+01}, {8.86842105e+01, 6.13778947e-01},
             {1.67368421e+02, 7.18031579e-02}, {2.46052632e+02, 1.59097368e-02},
             {3.24736842e+02, 4.84242105e-03}, {4.03421053e+02, 1.78918421e-03},
             {4.82105263e+02, 7.50905263e-04}, {5.60789474e+02, 3.52687368e-04},
             {6.39473684e+02, 1.74615789e-04}, {7.18157895e+02, 9.14143684e-05},
             {7.96842105e+02, 4.92782105e-05}, {8.75526316e+02, 2.82568421e-05},
             {9.54210526e+02, 1.63267368e-05}, {1.03289474e+03, 9.88289474e-06},
             {1.11157895e+03, 6.03487368e-06}, {1.19026316e+03, 3.77263158e-06},
             {1.26894737e+03, 2.40554737e-06}, {1.34763158e+03, 1.53987895e-06},
             {1.42631579e+03, 1.01542105e-06}, {1.50500000e+03, 6.67250000e-07},
             {1.58368421e+03, 4.48007368e-07}, {1.66236842e+03, 3.01508947e-07},
             {1.74105263e+03, 2.04563158e-07}, {1.81973684e+03, 1.40155263e-07},
             {1.89842105e+03, 9.61671579e-08}, {1.97710526e+03, 6.69074737e-08},
             {2.05578947e+03, 4.65284211e-08}, {2.13447368e+03, 3.26896842e-08},
             {2.21315789e+03, 2.28331579e-08}, {2.29184211e+03, 1.61397895e-08},
             {2.37052632e+03, 1.14483158e-08}, {2.44921053e+03, 8.09374211e-09},
             {2.52789474e+03, 5.73824211e-09}, {2.60657895e+03, 4.10271053e-09},
             {2.68526316e+03, 2.94212632e-09}, {2.76394737e+03, 2.09867895e-09},
             {2.84263158e+03, 1.49957895e-09}, {2.92131579e+03, 1.07367211e-09},
             {3.00000000e+03, 7.64700000e-10}}};
        for (auto [x, y] : vals) {
            REQUIRE(HP::SMHiggs(x).cxn(HP::Collider::LHC8,
                                       HP::Production::bbH) == Approx(y));
        }
    }
    SECTION("bbH LHC13") {
        static constexpr auto vals = std::array<std::pair<double, double>, 39>{
            {{1.00000000e+01, 1.13800000e+02}, {8.86842105e+01, 1.37805263e+00},
             {1.67368421e+02, 1.87584211e-01}, {2.46052632e+02, 4.69539474e-02},
             {3.24736842e+02, 1.59521053e-02}, {4.03421053e+02, 6.49665789e-03},
             {4.82105263e+02, 2.99794737e-03}, {5.60789474e+02, 1.53231579e-03},
             {6.39473684e+02, 8.25315789e-04}, {7.18157895e+02, 4.69609474e-04},
             {7.96842105e+02, 2.75711579e-04}, {8.75526316e+02, 1.70663158e-04},
             {9.54210526e+02, 1.06866737e-04}, {1.03289474e+03, 6.96918421e-05},
             {1.11157895e+03, 4.60520000e-05}, {1.19026316e+03, 3.11026842e-05},
             {1.26894737e+03, 2.14806316e-05}, {1.34763158e+03, 1.48385263e-05},
             {1.42631579e+03, 1.04615789e-05}, {1.50500000e+03, 7.47410000e-06},
             {1.58368421e+03, 5.41032632e-06}, {1.66236842e+03, 3.93694737e-06},
             {1.74105263e+03, 2.89311579e-06}, {1.81973684e+03, 2.14786842e-06},
             {1.89842105e+03, 1.59623158e-06}, {1.97710526e+03, 1.20894211e-06},
             {2.05578947e+03, 9.08293684e-07}, {2.13447368e+03, 6.91612105e-07},
             {2.21315789e+03, 5.30863158e-07}, {2.29184211e+03, 4.07168947e-07},
             {2.37052632e+03, 3.13527368e-07}, {2.44921053e+03, 2.42166316e-07},
             {2.52789474e+03, 1.88549474e-07}, {2.60657895e+03, 1.43955263e-07},
             {2.68526316e+03, 1.14157895e-07}, {2.76394737e+03, 9.04733684e-08},
             {2.84263158e+03, 7.12332632e-08}, {2.92131579e+03, 5.63350526e-08},
             {3.00000000e+03, 4.46700000e-08}}};
        for (auto [x, y] : vals) {
            REQUIRE(HP::SMHiggs(x).cxn(HP::Collider::LHC13,
                                       HP::Production::bbH) == Approx(y));
        }
    }
    SECTION("vbf LHC8") {
        static constexpr auto vals = std::array<std::pair<double, double>, 39>{
            {{1.00000000e+01, 5.47400000e+00}, {8.86842105e+01, 2.33531579e+00},
             {1.67368421e+02, 1.14526316e+00}, {2.46052632e+02, 6.31644737e-01},
             {3.24736842e+02, 3.76663158e-01}, {4.03421053e+02, 2.37018421e-01},
             {4.82105263e+02, 1.55236842e-01}, {5.60789474e+02, 1.05237579e-01},
             {6.39473684e+02, 7.26926316e-02}, {7.18157895e+02, 5.11887895e-02},
             {7.96842105e+02, 3.64136842e-02}, {8.75526316e+02, 2.64486842e-02},
             {9.54210526e+02, 1.92486316e-02}, {1.03289474e+03, 1.42378947e-02},
             {1.11157895e+03, 1.05583368e-02}, {1.19026316e+03, 7.89153158e-03},
             {1.26894737e+03, 5.93613684e-03}, {1.34763158e+03, 4.46592632e-03},
             {1.42631579e+03, 3.39531579e-03}, {1.50500000e+03, 2.57440000e-03},
             {1.58368421e+03, 1.96623158e-03}, {1.66236842e+03, 1.50164211e-03},
             {1.74105263e+03, 1.14886316e-03}, {1.81973684e+03, 8.82105263e-04},
             {1.89842105e+03, 6.75274737e-04}, {1.97710526e+03, 5.20674211e-04},
             {2.05578947e+03, 3.99574737e-04}, {2.13447368e+03, 3.08064737e-04},
             {2.21315789e+03, 2.37052632e-04}, {2.29184211e+03, 1.82421053e-04},
             {2.37052632e+03, 1.40657895e-04}, {2.44921053e+03, 1.07907895e-04},
             {2.52789474e+03, 8.33150526e-05}, {2.60657895e+03, 6.39476316e-05},
             {2.68526316e+03, 4.92111579e-05}, {2.76394737e+03, 3.77756316e-05},
             {2.84263158e+03, 2.89533684e-05}, {2.92131579e+03, 2.22183684e-05},
             {3.00000000e+03, 1.69500000e-05}}};
        for (auto [x, y] : vals) {
            REQUIRE(HP::SMHiggs(x).cxn(HP::Collider::LHC8,
                                       HP::Production::vbfH) == Approx(y));
        }
    }

    SECTION("vbf LHC13") {
        static constexpr auto vals = std::array<std::pair<double, double>, 39>{
            {{1.00000000e+01, 1.12100000e+01}, {8.86842105e+01, 5.32063158e+00},
             {1.67368421e+02, 2.85673684e+00}, {2.46052632e+02, 1.70965789e+00},
             {3.24736842e+02, 1.10100000e+00}, {4.03421053e+02, 7.46197368e-01},
             {4.82105263e+02, 5.25468421e-01}, {5.60789474e+02, 3.82373158e-01},
             {6.39473684e+02, 2.83584211e-01}, {7.18157895e+02, 2.14426316e-01},
             {7.96842105e+02, 1.64050526e-01}, {8.75526316e+02, 1.27789474e-01},
             {9.54210526e+02, 1.00122737e-01}, {1.03289474e+03, 7.95502632e-02},
             {1.11157895e+03, 6.35269474e-02}, {1.19026316e+03, 5.11176842e-02},
             {1.26894737e+03, 4.14112632e-02}, {1.34763158e+03, 3.36497895e-02},
             {1.42631579e+03, 2.75731579e-02}, {1.50500000e+03, 2.26150000e-02},
             {1.58368421e+03, 1.86670526e-02}, {1.66236842e+03, 1.54396842e-02},
             {1.74105263e+03, 1.28109474e-02}, {1.81973684e+03, 1.06668421e-02},
             {1.89842105e+03, 8.88301053e-03}, {1.97710526e+03, 7.43846316e-03},
             {2.05578947e+03, 6.22600000e-03}, {2.13447368e+03, 5.22983158e-03},
             {2.21315789e+03, 4.39657895e-03}, {2.29184211e+03, 3.70120000e-03},
             {2.37052632e+03, 3.12188421e-03}, {2.44921053e+03, 2.63170526e-03},
             {2.52789474e+03, 2.22722105e-03}, {2.60657895e+03, 1.88286842e-03},
             {2.68526316e+03, 1.59569474e-03}, {2.76394737e+03, 1.35222632e-03},
             {2.84263158e+03, 1.14627368e-03}, {2.92131579e+03, 9.73197895e-04},
             {3.00000000e+03, 8.25300000e-04}}};
        for (auto [x, y] : vals) {
            REQUIRE(HP::SMHiggs(x).cxn(HP::Collider::LHC13,
                                       HP::Production::vbfH) == Approx(y));
        }
    }
}
