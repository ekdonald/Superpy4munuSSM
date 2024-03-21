#include "predictions/Helpers.hpp"
#include "prettyprint.hpp" // IWYU pragma: keep
#include <catch2/catch_test_macros.hpp>
#include <vector>

using Higgs::predictions::Clamp, Higgs::predictions::validMassIn;

TEST_CASE("mass validation") {
    static const auto testGrid = std::vector{2., 3., 4., 5.};
    SECTION("negative mass") {
        auto mass = -1.;
        CHECK_FALSE(validMassIn<Clamp::none>(mass, testGrid));
        REQUIRE(mass == -1.);

        CHECK_FALSE(validMassIn<Clamp::upper>(mass, testGrid));
        REQUIRE(mass == -1.);

        CHECK_FALSE(validMassIn<Clamp::lower>(mass, testGrid));
        REQUIRE(mass == -1.);

        CHECK_FALSE(validMassIn<Clamp::both>(mass, testGrid));
        REQUIRE(mass == -1.);
    }
    SECTION("too low mass") {
        auto mass = 1.;
        CHECK_FALSE(validMassIn<Clamp::none>(mass, testGrid));
        REQUIRE(mass == 1.);

        CHECK_FALSE(validMassIn<Clamp::upper>(mass, testGrid));
        REQUIRE(mass == 1.);

        CHECK(validMassIn<Clamp::lower>(mass, testGrid));
        REQUIRE(mass == testGrid.front());
        mass = 1.;

        CHECK(validMassIn<Clamp::both>(mass, testGrid));
        REQUIRE(mass == testGrid.front());
    }
    SECTION("too high mass") {
        auto mass = 10.;
        CHECK_FALSE(validMassIn<Clamp::none>(mass, testGrid));
        REQUIRE(mass == 10.);

        CHECK(validMassIn<Clamp::upper>(mass, testGrid));
        REQUIRE(mass == testGrid.back());
        mass = 10.;

        CHECK_FALSE(validMassIn<Clamp::lower>(mass, testGrid));
        REQUIRE(mass == 10.);

        CHECK(validMassIn<Clamp::both>(mass, testGrid));
        REQUIRE(mass == testGrid.back());
    }
    SECTION("valid mass") {
        auto mass = 3.5;
        CHECK(validMassIn<Clamp::none>(mass, testGrid));
        REQUIRE(mass == 3.5);

        CHECK(validMassIn<Clamp::upper>(mass, testGrid));
        REQUIRE(mass == 3.5);

        CHECK(validMassIn<Clamp::lower>(mass, testGrid));
        REQUIRE(mass == 3.5);

        CHECK(validMassIn<Clamp::both>(mass, testGrid));
        REQUIRE(mass == 3.5);
    }
}
