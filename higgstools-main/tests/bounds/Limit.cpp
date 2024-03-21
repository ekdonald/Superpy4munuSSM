#include "Higgs/bounds/Limit.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <magic_enum.hpp>
#include <stdexcept>

namespace HB = Higgs::bounds;

TEST_CASE("read ChannelLimit") {
    REQUIRE_NOTHROW(HB::Limit::read(TESTLIMITS_PATH "ChannelLimit.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "ChannelLimit.json");
    REQUIRE(lim->id() == 19072749);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}

TEST_CASE("read ChannelWidthLimit") {
    REQUIRE_NOTHROW(HB::Limit::read(TESTLIMITS_PATH "ChannelWidthLimit.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "ChannelWidthLimit.json");
    REQUIRE(lim->id() == 12345678);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}

TEST_CASE("read LikelihoodLimit 1d") {
    REQUIRE_NOTHROW(HB::Limit::read(TESTLIMITS_PATH "LikelihoodLimit1d.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "LikelihoodLimit1d.json");
    REQUIRE(lim->id() == 1412864);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}

TEST_CASE("read LikelihoodLimit 2d") {
    REQUIRE_NOTHROW(HB::Limit::read(TESTLIMITS_PATH "LikelihoodLimit2d.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "LikelihoodLimit2d.json");
    REQUIRE(lim->id() == 14141414);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}

TEST_CASE("no LikelihoodLimit 3d") {
    REQUIRE_THROWS_AS(HB::Limit::read(TESTLIMITS_PATH "LikelihoodLimit3d.json"),
                      std::runtime_error);
}

TEST_CASE("read ChainDecayLimit") {
    REQUIRE_NOTHROW(HB::Limit::read(TESTLIMITS_PATH "ChainDecayLimit.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "ChainDecayLimit.json");
    REQUIRE(lim->id() == 123124125);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}

TEST_CASE("read PairDecayLimit") {
    REQUIRE_NOTHROW(HB::Limit::read(TESTLIMITS_PATH "PairDecayLimit.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "PairDecayLimit.json");
    REQUIRE(lim->id() == 7328721);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}

TEST_CASE("read PairProductionLimit") {
    REQUIRE_NOTHROW(
        HB::Limit::read(TESTLIMITS_PATH "PairProductionLimit.json"));
    auto lim = HB::Limit::read(TESTLIMITS_PATH "PairProductionLimit.json");
    REQUIRE(lim->id() == 123123123);
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(lim->reference()));
    CHECK_THAT(lim->to_string(), Catch::Matchers::ContainsSubstring(std::string(
                                     magic_enum::enum_name(lim->collider()))));
    CHECK_THAT(lim->to_string(),
               Catch::Matchers::ContainsSubstring(
                   std::string(magic_enum::enum_name(lim->experiment()))));
}
