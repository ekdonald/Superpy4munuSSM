#include "bounds/limits/BasicLimit.hpp"
#include "TestLimit.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "prettyprint.hpp"             // IWYU pragma: keep
#include "utilities/ArithmeticArray.hpp"
#include "utilities/Json.hpp"
#include "utilities/LinearInterpolator.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iomanip>
#include <magic_enum.hpp>
#include <map>
#include <memory>
#include <string>
#include <string_view>

using namespace std::string_literals;
namespace fs = std::filesystem;
namespace HB = Higgs::bounds;
namespace HP = Higgs::predictions;
namespace HU = Higgs::utilities;
using Catch::Approx;

TEST_CASE("Basic Limit Properties") {
    auto id{100U};
    auto reference{"XXXX.XXXXX"s};
    auto citeKey{"Aad:123456"};
    auto coll{HP::Collider::LHC13};
    auto exp{HP::Experiment::ATLAS};
    auto lumi{30.};
    fs::path file{"/path/.//to////limit/file.yml"};

    auto lim = HB::TestLimit::create(id, reference, citeKey, coll, exp, lumi,
                                     file.lexically_normal());

    REQUIRE(lim->id() == id);
    REQUIRE(lim->reference() == reference);
    REQUIRE(lim->citeKey() == citeKey);
    REQUIRE(lim->collider() == coll);
    REQUIRE(lim->experiment() == exp);
    REQUIRE(lim->luminosity() == lumi);
    REQUIRE(lim->loadedFrom() == file.lexically_normal().string());

    SECTION("json") {
        auto content = nlohmann::json{};
        content["id"] = id;
        content["reference"] = reference;
        content["citeKey"] = citeKey;
        content["collider"] = coll;
        content["experiment"] = exp;
        content["luminosity"] = lumi;

        INFO(std::setw(4) << content);

        REQUIRE(content.at("experiment") == magic_enum::enum_name(exp));
        REQUIRE(content.at("experiment").get<HP::Experiment>() == exp);
        REQUIRE(content.at("collider") == magic_enum::enum_name(coll));
        REQUIRE(content.at("collider").get<HP::Collider>() == coll);

        auto lim2 = HB::TestLimit::create(content, file.lexically_normal());

        REQUIRE(lim2->id() == id);
        REQUIRE(lim2->reference() == reference);
        REQUIRE(lim2->citeKey() == citeKey);
        REQUIRE(lim2->collider() == coll);
        REQUIRE(lim2->experiment() == exp);
        REQUIRE(lim2->luminosity() == lumi);
        REQUIRE(lim2->loadedFrom() == file.lexically_normal().string());
    }
}

TEST_CASE("limit mass ranges") {
    auto opts = HB::LimitOptions{};
    opts.applicableResolutionFac = 2.;
    const auto lim2 = HB::TestLimit(opts);
    const auto lim = HB::TestLimit(HB::LimitOptions{});

    auto extent = std::pair{100., 200.};
    auto massRes = HP::MassResolution{0.05, 3};
    auto range = lim2.applicableMassRange(extent, massRes);
    CHECK(range.first == Approx(84));
    CHECK(range.second == Approx(226));

    const auto extents =
        std::array{extent, std::pair{10., 20.}, std::pair{0., 0.1}};
    const auto massRess = std::array{massRes, HP::MassResolution{0.2, 6}};
    auto ranges = lim.applicableMassRanges(extents, massRess);
    REQUIRE(ranges.size() == 2);
    CHECK(ranges[0].first == Approx(96));
    CHECK(ranges[0].second == Approx(206.5));
    CHECK(ranges[1].first == Approx(6));
    CHECK(ranges[1].second == Approx(25));
}

TEST_CASE("within experimental resolution") {
    auto opts = HB::LimitOptions{};
    opts.applicableMassUnc = HP::MassUncEagerness::eager;
    const auto eager = HB::TestLimit(opts);
    opts.applicableMassUnc = HP::MassUncEagerness::cautious;
    const auto cautious = HB::TestLimit(opts);
    opts.applicableMassUnc = HP::MassUncEagerness::ignore;
    const auto ignore = HB::TestLimit(opts);

    auto m1 = HP::UncertainMass{95, 5};
    auto m2 = HP::UncertainMass{100, 2};

    auto res1 = HP::MassResolution{0., 1.};
    auto res2 = HP::MassResolution{0., 1.5};

    SECTION("invalid") {
        opts.applicableMassUnc = HP::MassUncEagerness{4};
        const auto invalid = HB::TestLimit(opts);
        REQUIRE_FALSE(invalid.withinExpRes(m1, m2, res1, res2));
    }

    REQUIRE(eager.withinExpRes(m1, m2, res1, res2) ==
            eager.withinExpRes(m2, m1, res2, res1));
    CHECK(eager.withinExpRes(m1, m2, res1, res2));
    CHECK_FALSE(eager.withinExpRes(HP::UncertainMass{95, 0}, m2, res1, res2));

    REQUIRE(cautious.withinExpRes(m1, m2, res1, res2) ==
            cautious.withinExpRes(m2, m1, res2, res1));
    CHECK_FALSE(cautious.withinExpRes(m1, m2, res1, res2));
    CHECK(cautious.withinExpRes(m1, m2, HP::MassResolution{0., 11.}, res2));

    REQUIRE(ignore.withinExpRes(m1, m2, res1, res2) ==
            ignore.withinExpRes(m2, m1, res2, res1));
    CHECK_FALSE(ignore.withinExpRes(m1, m2, res1, res2));
    CHECK(ignore.withinExpRes(m1, m2, res1, HP::MassResolution{0., 4.}));
}

TEST_CASE("getting limits with mass uncertainties") {

    auto opts = HB::LimitOptions{};
    opts.setLimitMassUnc = HP::MassUncEagerness::eager;
    const auto eager = HB::TestLimit(opts);
    opts.setLimitMassUnc = HP::MassUncEagerness::cautious;
    const auto cautious = HB::TestLimit(opts);
    opts.setLimitMassUnc = HP::MassUncEagerness::ignore;
    const auto ignore = HB::TestLimit(opts);

    auto dat = HU::LinearInterpolator<1, HU::ArithmeticArray<double, 2>>{
        std::array{std::vector{1., 10.}},
        std::vector{HU::ArithmeticArray{1., 1.},
                    HU::ArithmeticArray{10., 10.}}};

    auto testMasses = std::array{HP::UncertainMass{5, 2}};

    REQUIRE(eager.getLimit(dat, testMasses) <=
            ignore.getLimit(dat, testMasses));
    REQUIRE(ignore.getLimit(dat, testMasses) <=
            cautious.getLimit(dat, testMasses));

    CHECK(ignore.getLimit(dat, testMasses)[0] == Approx(5));
    CHECK(eager.getLimit(dat, testMasses)[0] == Approx(3));
    CHECK(cautious.getLimit(dat, testMasses)[0] == Approx(7));

    SECTION("invalid") {
        opts.setLimitMassUnc = HP::MassUncEagerness{4};
        const auto invalid = HB::TestLimit(opts);
        REQUIRE_NOTHROW(invalid.getLimit(dat, testMasses));
        REQUIRE(invalid.getLimit(dat, testMasses) ==
                cautious.getLimit(dat, testMasses));
    }
}
