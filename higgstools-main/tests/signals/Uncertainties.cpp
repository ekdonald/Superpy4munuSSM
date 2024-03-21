
#include "signals/Uncertainties.hpp"
#include "signals/JsonSupport.hpp"
#include "utilities/Json.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;

namespace HS = Higgs::signals;

TEST_CASE("UncertainValue") {
    SECTION("different construction methods") {
        {
            auto v = HS::UncertainValue{1, 0.1, 0.2};
            REQUIRE(v.central() == 1);
            REQUIRE(v.upperUncertainty() == 0.2);
            REQUIRE(v.lowerUncertainty() == 0.1);
            REQUIRE(v.symmetricUncertainty() == Approx(0.15));
            REQUIRE(v == v);
        }

        {
            auto v = HS::UncertainValue::fromInterval(0.9, 1, 1.2);
            REQUIRE(v.central() == 1);
            REQUIRE(v.upperUncertainty() == Approx(0.2));
            REQUIRE(v.lowerUncertainty() == Approx(0.1));
            REQUIRE(v.symmetricUncertainty() == Approx(0.15));
            REQUIRE(v == v);
        }

        {
            auto v = HS::UncertainValue{2};
            REQUIRE(v.central() == 2);
            REQUIRE(v.upperUncertainty() == 0);
            REQUIRE(v.lowerUncertainty() == 0);
            REQUIRE(v.symmetricUncertainty() == 0);
            REQUIRE(v == v);
        }

        {
            auto v = HS::UncertainValue{2};
            REQUIRE(v.central() == 2);
            REQUIRE(v.upperUncertainty() == 0);
            REQUIRE(v.lowerUncertainty() == 0);
            REQUIRE(v.symmetricUncertainty() == 0);
            REQUIRE(v == v);
        }
    }

    SECTION("comparison") {
        REQUIRE(HS::UncertainValue{5, 2} == HS::UncertainValue{5, 2, 2});
        REQUIRE(HS::UncertainValue{3} == HS::UncertainValue{3, 0});
        REQUIRE(HS::UncertainValue{3} == HS::UncertainValue{3, 0, 0});
        REQUIRE(HS::UncertainValue{3} ==
                HS::UncertainValue::fromInterval(3, 3, 3));
    }

    SECTION("error handling") {
        REQUIRE(HS::UncertainValue{-1, -2, -3} == HS::UncertainValue{-1, 2, 3});
        REQUIRE(HS::UncertainValue{-4, -1} == HS::UncertainValue{-4, 1});

        REQUIRE_THROWS_AS(HS::UncertainValue::fromInterval(2, 1, 3),
                          std::domain_error);
        REQUIRE_THROWS_AS(HS::UncertainValue::fromInterval(1, 3, 2),
                          std::domain_error);
    }

    SECTION("json conversion") {
        auto v1 = HS::UncertainValue{1, 0.6, 0.3};
        auto j1 = nlohmann::json(v1);
        INFO(j1);
        REQUIRE(j1 == nlohmann::json::array({0.4, 1., 1.3}));
        auto vv1 = j1.get<HS::UncertainValue>();
        REQUIRE(vv1.central() == v1.central());
        REQUIRE(vv1.lowerUncertainty() == Approx(v1.lowerUncertainty()));
        REQUIRE(vv1.upperUncertainty() == Approx(v1.upperUncertainty()));

        REQUIRE_THROWS_AS(
            nlohmann::json::array({2., 1., 1.5}).get<HS::UncertainValue>(),
            std::domain_error);
    }
}
