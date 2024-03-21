#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "TestLimit.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace std::string_literals;
namespace fs = std::filesystem;
using Catch::Approx;

TEST_CASE("Applied Limit") {
    auto id{100U};
    auto reference{"XXXX.XXXX"s};
    auto citeKey{"Aad:123456"};
    auto coll{Higgs::predictions::Collider::LHC13};
    auto exp{Higgs::predictions::Experiment::ATLAS};
    auto lumi{30.};
    fs::path file{"/path/.//to////limit/file.yml"};

    auto lim = Higgs::bounds::TestLimit::create(
        id, reference, citeKey, coll, exp, lumi, file.lexically_normal());

    auto particles = std::vector<std::string>{"H1", "H2"};
    auto applied = Higgs::bounds::AppliedLimit{lim, 1.1, 0.9, particles};
    SECTION("properties") {
        REQUIRE(applied.limit()->id() == id);
        REQUIRE(applied.limit()->reference() == reference);
        REQUIRE(applied.limit()->citeKey() == citeKey);
        REQUIRE(applied.limit()->collider() == coll);
        REQUIRE(applied.limit()->experiment() == exp);
        REQUIRE(applied.limit()->luminosity() == lumi);
        REQUIRE(applied.limit()->loadedFrom() ==
                file.lexically_normal().string());
    }

    SECTION("copy and move") {
        auto applied2 = applied;
        REQUIRE(applied2.limit()->id() == id);
        REQUIRE(applied2.limit()->reference() == reference);
        REQUIRE(applied2.limit()->citeKey() == citeKey);
        REQUIRE(applied2.limit()->collider() == coll);
        REQUIRE(applied2.limit()->experiment() == exp);
        REQUIRE(applied2.limit()->luminosity() == lumi);
        REQUIRE(applied2.limit()->loadedFrom() ==
                file.lexically_normal().string());
        auto applied3 = std::move(applied2);
        REQUIRE(applied3.limit()->id() == id);
        REQUIRE(applied3.limit()->reference() == reference);
        REQUIRE(applied3.limit()->citeKey() == citeKey);
        REQUIRE(applied3.limit()->collider() == coll);
        REQUIRE(applied3.limit()->experiment() == exp);
        REQUIRE(applied3.limit()->luminosity() == lumi);
        REQUIRE(applied3.limit()->loadedFrom() ==
                file.lexically_normal().string());
    }

    SECTION("testlimit") {
        REQUIRE(applied.obsRatio() == Approx(1.1));
        REQUIRE(applied.expRatio() == Approx(0.9));
        REQUIRE(applied.contributingParticles() == particles);
    }
}
