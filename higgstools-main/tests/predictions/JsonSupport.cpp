#include "predictions/JsonSupport.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "utilities/Json.hpp"
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <magic_enum.hpp>
#include <map>
#include <memory>
// IWYU pragma: no_forward_declare mpl_::na

using Catch::Approx;

TEMPLATE_TEST_CASE("enum json conversions", "", Higgs::predictions::Collider,
                   Higgs::predictions::ECharge, Higgs::predictions::Experiment,
                   Higgs::predictions::Production, Higgs::predictions::Decay,
                   Higgs::predictions::ChainDecay,
                   Higgs::predictions::Coupling) {
    for (const auto &e : magic_enum::enum_values<TestType>()) {
        const nlohmann::json j = magic_enum::enum_name(e);
        auto convJson = nlohmann::json{};
        auto convE = TestType{};

        to_json(convJson, e);
        REQUIRE(convJson == j);
        from_json(convJson, convE);
        REQUIRE(convE == e);

        from_json(j, convE);
        REQUIRE(convE == e);
        to_json(convJson, convE);
        REQUIRE(convJson == j);

        REQUIRE(j.get<TestType>() == e);

        REQUIRE_THROWS_AS(from_json(nlohmann::json("not an enum name"), convE),
                          Higgs::utilities::BadEnumRead);
    }
}

TEST_CASE("Json conversions") {
    SECTION("massResolution") {
        const auto res = Higgs::predictions::MassResolution{0.2, 20};
        auto convJson = nlohmann::json{};
        auto convRes = Higgs::predictions::MassResolution{};

        to_json(convJson, res);
        REQUIRE(convJson ==
                nlohmann::json{{"absolute", 20}, {"relative", 0.2}});
        from_json(convJson, convRes);
        REQUIRE(res.absolute == convRes.absolute);
        REQUIRE(res.relative == convRes.relative);

        const auto json = nlohmann::json{{"absolute", 10}, {"relative", 0.1}};
        from_json(json, convRes);
        REQUIRE(Approx(10.) == convRes.absolute);
        REQUIRE(Approx(0.1) == convRes.relative);
        to_json(convJson, convRes);
        REQUIRE(convJson == json);
    }
}
