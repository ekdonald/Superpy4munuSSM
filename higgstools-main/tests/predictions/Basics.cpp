#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <magic_enum.hpp>
#include <memory>
#include <string>

namespace HP = Higgs::predictions;
using namespace std::string_literals;

TEST_CASE("collider classification") {
    for (auto c : {HP::Collider::LHC8, HP::Collider::LHC13}) {
        REQUIRE(HP::classifyCollider(c) == HP::ColliderType::pp);
    }
    REQUIRE(HP::classifyCollider(HP::Collider::LEP) == HP::ColliderType::ee);
}

TEST_CASE("HiggsPredictions errors") {
    SECTION("InvalidInput") {
        auto e = HP::InvalidInput("example error");
        REQUIRE(e.what() == "example error"s);

        e = HP::InvalidInput("example error"s);
        REQUIRE(e.what() == "example error"s);

        std::unique_ptr<std::exception> ep =
            std::make_unique<HP::InvalidInput>("example error");
        REQUIRE(ep->what() == "example error"s);
    }

    SECTION("InvalidChannel") {
        using Catch::Matchers::ContainsSubstring;
        REQUIRE_THAT(
            HP::InvalidChannel(HP::ECharge::single, HP::Decay::bb).what(),
            ContainsSubstring(
                std::string(magic_enum::enum_name<HP::Decay::bb>())));
        REQUIRE_THAT(
            HP::InvalidChannel(HP::ECharge::single, HP::Production::H).what(),
            ContainsSubstring(
                std::string(magic_enum::enum_name<HP::Production::H>())));
        REQUIRE_THAT(
            HP::InvalidChannel(HP::ColliderType::pp, HP::Production::eeHbb)
                .what(),
            ContainsSubstring(
                std::string(magic_enum::enum_name<HP::ColliderType::pp>())) &&
                ContainsSubstring(std::string(
                    magic_enum::enum_name<HP::Production::eeHbb>())));
        REQUIRE_THAT(
            HP::InvalidChannel(HP::Production::bbH, HP::Decay::taunu).what(),
            ContainsSubstring(
                std::string(magic_enum::enum_name<HP::Production::bbH>())) &&
                ContainsSubstring(
                    std::string(magic_enum::enum_name<HP::Decay::taunu>())));
    }
}
