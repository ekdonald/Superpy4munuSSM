#include "bounds/Likelihood.hpp"
#include "utilities/LinearInterpolator.hpp"
#include <algorithm>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <exception>
#include <initializer_list>
#include <limits>
#include <memory>
#include <vector>

using Catch::Approx;
using Higgs::bounds::Likelihood::expCLs, Higgs::bounds::Likelihood::obsCLs,
    Higgs::bounds::Likelihood::limitsFromLlh;

TEST_CASE("CLs computations") {
    const auto testvals = std::array{0.,    std::numeric_limits<double>::min(),
                                     1e-10, 1e-3,
                                     0.1,   1.,
                                     5.,    10.,
                                     1e3,   1e10,
                                     1e100, std::numeric_limits<double>::max()};

    SECTION("basic sanity checks") {

        for (auto x : testvals) {
            CHECK(expCLs(x) >= 0.);
            for (auto y : testvals) {
                INFO(x << " " << y);
                CHECK(obsCLs(x, y) >= 0.);
                // these are potentially dangerous for the expansion
                CHECK(obsCLs(x, y * 2.99) >= 0.);
            }
        }
    }

    SECTION("asymptotic behaviour") {
        CHECK(expCLs(0.) == Approx(1.));
        CHECK(expCLs(2e3) == Approx(0.));
        CHECK(expCLs(1e6) == 0.);

        for (auto x : testvals) {
            INFO(x);
            CHECK(obsCLs(0., x) == Approx(obsCLs(1e-100, x)));
            CHECK(obsCLs(x, 0.) == Approx(obsCLs(x, 1e-100)));
            CHECK(obsCLs(x, x) == Approx(expCLs(x)));
        }
    }

    SECTION("some test values") {
        CHECK(expCLs(0) == Approx(1.));
        CHECK(expCLs(1) == Approx(3.173105e-1));
        CHECK(expCLs(10) == Approx(1.5654022580025018e-3));

        CHECK(obsCLs(1, 0) == Approx(0.5942867086725301));
        CHECK(obsCLs(2, 1) == Approx(0.24015345570732952));
        CHECK(obsCLs(20, 10) == Approx(8.649796248556664e-4));
        CHECK(obsCLs(1, 2) == Approx(0.21652859987383824));
        CHECK(obsCLs(10, 20) == Approx(1.8458535698646177e-05));
        CHECK(obsCLs(0, 0.001) == Approx(0.9995));
        CHECK(obsCLs(0, 1) == Approx(0.606531));
        CHECK(obsCLs(0, 5.) == Approx(0.082085));
    }
}

TEST_CASE("Likelihood limit") {
    auto exp = [](double s) { return 2 * std::pow(s, 2); };
    auto obs = [](double s) { return 0.3 * std::pow(s, 3); };

    SECTION("equal ratios for equal llhs") {
        auto [e1, e2] = limitsFromLlh(exp, exp);
        REQUIRE(e1 == Approx(e2));
        REQUIRE(expCLs(exp(1 / e1)) == Approx(0.05));
        auto [o1, o2] = limitsFromLlh(obs, obs);
        REQUIRE(o1 == Approx(o2));
        REQUIRE(expCLs(obs(1 / o1)) == Approx(0.05));
    }

    SECTION("test values") {
        auto [e, o] = limitsFromLlh(exp, obs);
        CHECK(expCLs(exp(1 / e)) == Approx(0.05));
        CHECK(obsCLs(exp(1 / o), obs(1 / o)) == Approx(0.05));
        CHECK(e == Approx(1 / 1.3859038234));
        CHECK(o == Approx(1 / 2.1219356246));
    }

    SECTION("different predicted rate") {
        auto [e, o] = limitsFromLlh(exp, obs);

        for (const auto scale : {0.4, 0.8, 1., 2., 4.4, 5.13}) {
            auto expS = [&exp, scale](double s) { return exp(s * scale); };
            auto obsS = [&obs, scale](double s) { return obs(s * scale); };

            auto [eS, oS] = limitsFromLlh(expS, obsS);
            CHECK(e == Approx(eS / scale));
            CHECK(o == Approx(oS / scale));
        }
    }

    SECTION("clamped likelihoods") {
        auto clampE = [&exp](double s) {
            return s < 0.4 ? 100 * exp(s) : 100 * exp(0.4);
        };
        auto clampO = [&obs](double s) {
            return s < 0.4 ? 200 * obs(s) : 200 * obs(0.4);
        };

        REQUIRE(expCLs(clampE(1e3)) < 0.05);
        REQUIRE(obsCLs(clampE(1e3), clampO(1e3)) < 0.05);

        auto [e, o] = limitsFromLlh(clampE, clampO);
        CHECK(e == Approx(7.2155079049));
        CHECK(o == Approx(2.8093271046));
    }

    SECTION("interpolated likelihoods") {
        const auto massG = std::vector{10., 20., 30.};
        const auto rateG = std::vector{0., 2.};
        const auto expLlh = Higgs::utilities::LinearInterpolator(
            std::array{massG, rateG}, {0., 5.5, 3., 8., 10., 15.});
        const auto obsLlh = Higgs::utilities::LinearInterpolator(
            std::array{massG, rateG}, {0., 4., 4., 9., 12., 20.});
        constexpr auto pred = 1.5;
        constexpr auto mass = 15.;

        const auto expInterp = [&](double s) {
            return expLlh(
                {mass, std::clamp(s * pred, rateG.front(), rateG.back())});
        };
        const auto obsInterp = [&](double s) {
            return obsLlh(
                {mass, std::clamp(s * pred, rateG.front(), rateG.back())});
        };

        REQUIRE(expLlh({20., 0.}) == Approx(3.));
        REQUIRE(obsLlh({10., 2.}) == Approx(4.));
        REQUIRE(expCLs(expInterp(1e3)) < 0.05);
        REQUIRE(obsCLs(expInterp(1e3), obsInterp(1e3)) < 0.05);

        auto [e, o] = limitsFromLlh(expInterp, obsInterp);
        CHECK(e == Approx(1.681643));
        CHECK(o == Approx(1.772731));
    }

    SECTION("all scale factors allowed") {
        auto tooLow = [](double s) { return std::min(s * s, 2.); };

        INFO("all s are allowed");
        REQUIRE(expCLs(tooLow(1e10)) > 0.05);
        REQUIRE_NOTHROW(limitsFromLlh(tooLow, tooLow));
        auto [e, o] = limitsFromLlh(tooLow, tooLow);
        CHECK(e == 0);
        CHECK(o == 0);
    }
}
