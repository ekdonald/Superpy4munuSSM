#include "utilities/Algorithm.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/unreachable_sentinel.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/view/single.hpp>
#include <set>
#include <unordered_map>
#include <vector>

namespace {
struct Dummy {
    double value;
    double weight;
};
} // namespace

namespace HU = Higgs::utilities;
using Catch::Approx;

TEST_CASE("mean") {
    const auto ofValue = [](const Dummy &d) { return d.value; };
    const auto byWeight = [](const Dummy &d) { return d.weight; };
    SECTION("single range") {

        SECTION("Edge Cases") {
            const auto empty = std::vector<Dummy>{};
            CHECK(HU::weightedMean(empty, byWeight, ofValue) == 0.);
            CHECK(HU::mean(empty, ofValue) == 0.);
            const auto zeros = std::vector<Dummy>{{0, 0}, {0, 0}};
            CHECK(HU::weightedMean(zeros, byWeight, ofValue) == 0.);
            CHECK(HU::mean(zeros, ofValue) == 0.);
        }

        SECTION("Normal cases") {
            auto values = std::vector<Dummy>{{1, 3}, {2, 4}, {5, 0}, {2, 1}};
            CHECK(HU::weightedMean(values, byWeight, ofValue) ==
                  Approx(13. / 8.));
            CHECK(HU::mean(values, ofValue) == Approx(10. / 4.));
            values = std::vector<Dummy>{{1, 0.3}, {2, 0.4}, {2, 0.1}};
            CHECK(HU::weightedMean(values, byWeight, ofValue) ==
                  Approx(1.3 / 0.8));
            CHECK(HU::mean(values, ofValue) == Approx(5. / 3.));
        }
    }

    SECTION("two ranges") {
        constexpr auto byWeights = [](const Dummy &d1, const Dummy &d2) {
            return d1.weight * d2.weight;
        };
        auto v1 = std::vector<Dummy>{{1, 3}, {2, 4}, {5, 0}, {2, 1}};
        auto v2 = std::vector<Dummy>{{2, 2}, {1, 4}};

        SECTION("edge Cases") {
            const auto empty = std::vector<Dummy>{};
            const auto zeros = std::vector<Dummy>{{0, 0}, {0, 0}};
            CHECK(HU::weightedMean(v1, empty, byWeights, ofValue) ==
                  std::array{0., 0.});
            CHECK(HU::weightedMean(empty, v2, byWeights, ofValue) ==
                  std::array{0., 0.});
            CHECK(HU::weightedMean(v2, zeros, byWeights, ofValue) ==
                  std::array{0., 0.});
            CHECK(HU::weightedMean(zeros, v1, byWeights, ofValue) ==
                  std::array{0., 0.});
        }

        SECTION("reduces to 1D case due to simple weight function") {
            auto [m1, m2] = HU::weightedMean(v1, v2, byWeights, ofValue);
            CHECK(m1 == Approx(HU::weightedMean(v1, byWeight, ofValue)));
            CHECK(m2 == Approx(HU::weightedMean(v2, byWeight, ofValue)));
        }
    }

    SECTION("three ranges") {
        constexpr auto byWeights = [](const Dummy &d1, const Dummy &d2,
                                      const Dummy &d3) {
            return d1.weight * d2.weight * d3.weight;
        };
        auto v1 = std::vector<Dummy>{{1, 3}, {2, 4}, {5, 0}, {2, 1}};
        auto v2 = std::vector<Dummy>{{2, 2}, {1, 4}};
        auto v3 = std::vector<Dummy>{{3, 2}, {5, 3}};

        SECTION("edge Cases") {
            const auto empty = std::vector<Dummy>{};
            const auto zeros = std::vector<Dummy>{{0, 0}, {0, 0}};
            CHECK(HU::weightedMean(v1, empty, v2, byWeights, ofValue) ==
                  std::array{0., 0., 0.});
            CHECK(HU::weightedMean(empty, v3, v2, byWeights, ofValue) ==
                  std::array{0., 0., 0.});
            CHECK(HU::weightedMean(v2, v1, zeros, byWeights, ofValue) ==
                  std::array{0., 0., 0.});
            CHECK(HU::weightedMean(v3, zeros, v1, byWeights, ofValue) ==
                  std::array{0., 0., 0.});
        }

        SECTION("reduces to 1D case due to simple weight function") {
            auto [m1, m2, m3] =
                HU::weightedMean(v1, v2, v3, byWeights, ofValue);
            CHECK(m1 == Approx(HU::weightedMean(v1, byWeight, ofValue)));
            CHECK(m2 == Approx(HU::weightedMean(v2, byWeight, ofValue)));
            CHECK(m3 == Approx(HU::weightedMean(v3, byWeight, ofValue)));
        }
    }
}

TEST_CASE("get with default") {
    SECTION("map") {
        auto m = std::map<int, int>{};
        m[1] = 10;
        REQUIRE(HU::getWithDefault(m, 1, -1) == 10);
        REQUIRE(HU::getWithDefault(m, 2, -1) == -1);
    }
    SECTION("unordered map") {
        auto m = std::unordered_map<int, int>{};
        m[1] = 10;
        REQUIRE(HU::getWithDefault(m, 1, -1) == 10);
        REQUIRE(HU::getWithDefault(m, 2, -1) == -1);
    }
}

TEST_CASE("erase from set if") {
    auto testset = std::set{1, 2, 3, 4, 5};
    REQUIRE(testset.size() == 5);
    HU::eraseFromSetIf(testset, [](int x) { return x >= 4; });
    REQUIRE(testset.size() == 3);
}
