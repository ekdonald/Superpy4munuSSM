#include "utilities/Json.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <exception>
#include <memory>
#include <string>

using namespace std::string_literals;
using Catch::Approx;

TEST_CASE("Json errors") {
    SECTION("Bad enum read") {
        auto e = Higgs::utilities::BadEnumRead("example error");
        REQUIRE(e.what() == "example error"s);

        e = Higgs::utilities::BadEnumRead("example error"s);
        REQUIRE(e.what() == "example error"s);

        std::unique_ptr<std::exception> ep =
            std::make_unique<Higgs::utilities::BadEnumRead>("example error");
        REQUIRE(ep->what() == "example error"s);
    }
    SECTION("BadFileRead") {
        auto e = Higgs::utilities::BadFieldRead("example error");
        REQUIRE(e.what() == "example error"s);

        e = Higgs::utilities::BadFieldRead("example error"s);
        REQUIRE(e.what() == "example error"s);

        std::unique_ptr<std::exception> ep =
            std::make_unique<Higgs::utilities::BadFieldRead>("example error");
        REQUIRE(ep->what() == "example error"s);
    }
}

TEST_CASE("json read helpers") {
    SECTION("readAs") {
        using Higgs::utilities::readAs;
        auto j = nlohmann::json{};
        j["num"] = -1;

        REQUIRE(readAs<int>(j, "num") == -1);
        REQUIRE(readAs<double>(j, "num") == Approx(-1.));
        REQUIRE_THROWS_AS(readAs<std::string>(j, "num"),
                          Higgs::utilities::BadFieldRead);
    }

    SECTION("readIfPresent") {
        using Higgs::utilities::readIfPresent;
        auto j = nlohmann::json{};
        j["num"] = -1;

        REQUIRE(readIfPresent<int>(j, "num") == -1);
        REQUIRE(readIfPresent<double>(j, "num") == Approx(-1.));
        REQUIRE_THROWS_AS(readIfPresent<std::string>(j, "num"),
                          Higgs::utilities::BadFieldRead);
        REQUIRE(readIfPresent<int>(j, "/num") == -1);
        REQUIRE(readIfPresent<int>(j, "foo") == int{});
        REQUIRE(readIfPresent<int>(j, "foo", 10) == 10);
    }

    SECTION("readWithOptDefault") {
        using Higgs::utilities::readWithOptDefault;

        const auto j = nlohmann::json{{"foo", 1}, {"bar", 4}};

        REQUIRE(readWithOptDefault<int>(j, "foo", std::nullopt) == 1);
        REQUIRE(readWithOptDefault<int>(j, "foo", 3) == 1);
        REQUIRE(readWithOptDefault<int>(j, "bar", std::nullopt) == 4);
        REQUIRE(readWithOptDefault<int>(j, "bar", 1) == 4);
        REQUIRE_THROWS_AS(readWithOptDefault<int>(j, "foo2", std::nullopt),
                          nlohmann::json::out_of_range);
        REQUIRE(readWithOptDefault<int>(j, "foo2", 2) == 2);
    }
}
