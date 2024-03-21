#include "signals/JsonSupport.hpp"
#include "Higgs/signals/Measurement.hpp"
#include "signals/Uncertainties.hpp"
#include "utilities/Json.hpp"
#include <catch2/catch_test_macros.hpp>

namespace HS = Higgs::signals;
TEST_CASE("read correlation matrix") {
    const auto j = nlohmann::json::parse(R"(
            {
                "a":{
                    "b":2,
                    "c":3
                },
                "b":{
                    "c":4
                }
            }
        )");

    SECTION("symmetric dict read") {
        CHECK(HS::readFromSymmNestedDicts(j, "a", "b", 0) == 2);
        CHECK(HS::readFromSymmNestedDicts(j, "a", "b", -1) ==
              HS::readFromSymmNestedDicts(j, "b", "a", -2));

        CHECK(HS::readFromSymmNestedDicts(j, "a", "c", 0) == 3);
        CHECK(HS::readFromSymmNestedDicts(j, "a", "c", -1) ==
              HS::readFromSymmNestedDicts(j, "c", "a", -2));

        CHECK(HS::readFromSymmNestedDicts(j, "b", "c", 0) == 4);
        CHECK(HS::readFromSymmNestedDicts(j, "b", "c", -1) ==
              HS::readFromSymmNestedDicts(j, "c", "b", -2));

        CHECK(HS::readFromSymmNestedDicts(j, "c", "x", -1) == -1);
        CHECK(HS::readFromSymmNestedDicts(j, "x", "a", -2) == -2);
    }

    auto mat = Eigen::MatrixXd(3, 3);
    HS::readCorrelationMatrix(j, {"a", "b", "c"}, mat);
    REQUIRE(mat(0, 1) == 2);
    REQUIRE(mat(0, 2) == 3);
    REQUIRE(mat(1, 0) == mat(0, 1));
    REQUIRE(mat(1, 2) == 4);
    REQUIRE(mat(2, 0) == mat(0, 2));
    REQUIRE(mat(2, 1) == mat(1, 2));

    HS::readCorrelationMatrix(j, {"c", "b", "a"}, mat);
    REQUIRE(mat(0, 1) == 4);
    REQUIRE(mat(0, 2) == 3);
    REQUIRE(mat(1, 0) == mat(0, 1));
    REQUIRE(mat(1, 2) == 2);
    REQUIRE(mat(2, 0) == mat(0, 2));
    REQUIRE(mat(2, 1) == mat(1, 2));
}
