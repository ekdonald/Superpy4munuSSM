#include "utilities/LinearInterpolator.hpp"
#include "prettyprint.hpp" // IWYU pragma: keep
#include <algorithm>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/algorithm/min.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/range/conversion.hpp>
#include <stdexcept>
#include <vector>

namespace utils = Higgs::utilities;
namespace utilD = Higgs::utilities::Detail;
using Catch::Approx;

TEST_CASE("findAnchor") {
    using Higgs::utilities::Detail::findAnchor;
    SECTION("normal case") {
        const auto edge = std::vector{0., 2., 4.};
        CHECK(findAnchor(edge, 0) == 0);
        CHECK(findAnchor(edge, 1) == 0);
        CHECK(findAnchor(edge, 2) == 1);
        CHECK(findAnchor(edge, 3) == 1);
        CHECK(findAnchor(edge, 4) == 1);
    }

    SECTION("two elements") {
        const auto edge = std::vector{0., 2.};
        CHECK(findAnchor(edge, 0) == 0);
        CHECK(findAnchor(edge, 1) == 0);
        CHECK(findAnchor(edge, 2) == 0);
    }
}

TEST_CASE("RegularGridView") {
    using Higgs::utilities::RegularGridView;
    SECTION("1d") {
        const auto data = std::vector{1, 2, 3, 4};
        const auto edges = std::array{std::vector{1., 2., 3., 4.}};
        const auto grid = RegularGridView{edges, data};
        for (std::size_t i = 0; i != data.size(); ++i) {
            CHECK(grid[{i}] == data[i]);
        }
        CHECK(grid.extent()[0] == std::pair{1., 4.});
        CHECK(grid.clamp(std::array{10.}) == std::array{4.});
        {
            auto located = grid.locate(std::array{2.5});
            CHECK(located.anchor[0] == 1);
            CHECK(located.normDist[0] == 0.5);
        }
        {
            auto located = grid.locate(std::array{-1.});
            CHECK(located.anchor[0] == 0);
            CHECK(located.normDist[0] == 0.);
        }
        {
            auto located = grid.locate(std::array{1e3});
            CHECK(located.anchor[0] == 2);
            CHECK(located.normDist[0] == 1.0);
        }
    }

    SECTION("2d") {
        const auto data =
            std::vector<int>{11, 12, 13, 21, 22, 23, 31, 32, 33, 41, 42, 43};
        const auto edges =
            std::array{std::vector{1., 2., 3.}, std::vector{.1, .2, .3, .4}};
        const auto grid = RegularGridView{edges, data};
        for (std::size_t i = 0; i != data.size(); ++i) {
            CHECK(grid[{i / 4, i % 4}] == data[i]);
        }
        CHECK(grid.extent()[0] == std::pair{1., 3.});
        CHECK(grid.extent()[1] == std::pair{.1, .4});
        CHECK(grid.clamp(std::array{5., -1.}) == std::array{3., .1});
        {
            auto located = grid.locate(std::array{2.5, 0.32});
            CHECK(located.anchor == std::array<std::size_t, 2>{1, 2});
            CHECK(located.normDist[0] == Approx(0.5));
            CHECK(located.normDist[1] == Approx(0.2));
        }
        {
            auto located = grid.locate(std::array{8., -10.});
            CHECK(located.anchor == std::array<std::size_t, 2>{1, 0});
            CHECK(located.normDist == std::array{1.0, 0.});
        }
    }

    SECTION("default constructed") {
        const auto empty1 = RegularGridView<int, 1>{};
        const auto empty2 = RegularGridView<int, 2>{};

        CHECK(empty1.extent()[0] == std::pair<double, double>{});
        CHECK(empty2.extent()[0] == std::pair<double, double>{});
        CHECK(empty2.extent()[1] == std::pair<double, double>{});

        CHECK(empty1[{0}] == 0);
        CHECK(empty2[{0, 0}] == 0);

        CHECK(empty1.clamp(std::array{1.}) == std::array<double, 1>{});
        CHECK(empty2.clamp(std::array{1., -4e3}) == std::array<double, 2>{});

        CHECK(empty1.locate(std::array{4.2}).anchor ==
              std::array<std::size_t, 1>{});
        CHECK(empty1.locate(std::array{4.2}).normDist ==
              std::array<double, 1>{});
        CHECK(empty2.locate(std::array{2., 1.}).anchor ==
              std::array<std::size_t, 2>{});
        CHECK(empty2.locate(std::array{2., 1.}).normDist ==
              std::array<double, 2>{});
    }

    SECTION("Invalid") {
        const auto data = std::vector<int>{1, 2};
        REQUIRE_THROWS_AS(RegularGridView(std::array{std::vector{1.}}, data),
                          std::out_of_range);
        REQUIRE_THROWS_AS(
            RegularGridView(std::array{std::vector{1., 2., 3.}}, data),
            std::out_of_range);
        REQUIRE_THROWS_AS(
            RegularGridView(std::array{std::vector{2., 1.}}, data),
            std::runtime_error);
    }
}

TEST_CASE("RegularGrid") {
    using Higgs::utilities::RegularGrid;
    SECTION("1d") {
        const auto data = std::vector{1, 2, 3, 4};
        const auto edges = std::array{std::vector{1., 2., 3., 4.}};
        const auto grid = RegularGrid{edges, data};
        for (std::size_t i = 0; i != data.size(); ++i) {
            CHECK(grid[{i}] == data[i]);
        }
        CHECK(grid.extent()[0] == std::pair{1., 4.});
        CHECK(grid.clamp(std::array{10.}) == std::array{4.});
        {
            auto located = grid.locate(std::array{2.5});
            CHECK(located.anchor[0] == 1);
            CHECK(located.normDist[0] == 0.5);
        }
        {
            auto located = grid.locate(std::array{-1.});
            CHECK(located.anchor[0] == 0);
            CHECK(located.normDist[0] == 0.);
        }
        {
            auto located = grid.locate(std::array{1e3});
            CHECK(located.anchor[0] == 2);
            CHECK(located.normDist[0] == 1.0);
        }

        SECTION("copy and move") {
            auto grid2 = grid; // copy construction
            REQUIRE(grid2.extent() == grid.extent());
            auto grid3 = std::move(grid2); // move construction
            REQUIRE(grid3.extent() == grid.extent());
            grid2 = grid3; // assignment
            REQUIRE(grid2.extent() == grid.extent());
            grid3 = std::move(grid2); // move assignment
            REQUIRE(grid3.extent() == grid.extent());
            grid3 = grid3; // self assignment
            REQUIRE(grid3.extent() == grid.extent());
            grid3 = std::move(grid3); // self move assignment
            REQUIRE(grid3.extent() == grid.extent());
        }
    }

    SECTION("2d") {
        const auto data =
            std::vector<int>{11, 12, 13, 21, 22, 23, 31, 32, 33, 41, 42, 43};
        const auto edges =
            std::array{std::vector{1., 2., 3.}, std::vector{.1, .2, .3, .4}};
        const auto grid = RegularGrid{edges, data};
        for (std::size_t i = 0; i != data.size(); ++i) {
            CHECK(grid[{i / 4, i % 4}] == data[i]);
        }
        CHECK(grid.extent()[0] == std::pair{1., 3.});
        CHECK(grid.extent()[1] == std::pair{.1, .4});
        CHECK(grid.clamp(std::array{5., -1.}) == std::array{3., .1});
        {
            auto located = grid.locate(std::array{2.5, 0.32});
            CHECK(located.anchor == std::array<std::size_t, 2>{1, 2});
            CHECK(located.normDist[0] == Approx(0.5));
            CHECK(located.normDist[1] == Approx(0.2));
        }
        {
            auto located = grid.locate(std::array{8., -10.});
            CHECK(located.anchor == std::array<std::size_t, 2>{1, 0});
            CHECK(located.normDist == std::array{1.0, 0.});
        }
    }

    SECTION("trivial 2d") {
        const auto data = std::vector<int>{11, 12, 13};
        const auto edges = std::array{std::vector{1., 2., 3.}, std::vector{1.}};
        const auto grid = RegularGrid{edges, data};
        for (std::size_t i = 0; i != data.size(); ++i) {
            CHECK(grid[{i, 0}] == data[i]);
        }
        CHECK(grid.extent()[0] == std::pair{1., 3.});
        CHECK(grid.extent()[1] == std::pair{1., 1.});
        CHECK(grid.clamp(std::array{5., -1.}) == std::array{3., 1.});
        {
            auto located = grid.locate(std::array{2.5, 20.});
            CHECK(located.anchor == std::array<std::size_t, 2>{1, 0});
            CHECK(located.normDist[0] == Approx(0.5));
            CHECK(located.normDist[1] == 0);
        }
        {
            auto located = grid.locate(std::array{8., -10.});
            CHECK(located.anchor == std::array<std::size_t, 2>{1, 0});
            CHECK(located.normDist == std::array{1.0, 0.});
        }
    }

    SECTION("default constructed") {
        const auto empty1 = RegularGrid<int, 1>{};
        const auto empty2 = RegularGrid<int, 2>{};

        CHECK(empty1.extent()[0] == std::pair<double, double>{});
        CHECK(empty2.extent()[0] == std::pair<double, double>{});
        CHECK(empty2.extent()[1] == std::pair<double, double>{});

        CHECK(empty1[{0}] == 0);
        CHECK(empty2[{0, 0}] == 0);

        CHECK(empty1.clamp(std::array{1.}) == std::array<double, 1>{});
        CHECK(empty2.clamp(std::array{1., -4e3}) == std::array<double, 2>{});

        CHECK(empty1.locate(std::array{4.2}).anchor ==
              std::array<std::size_t, 1>{});
        CHECK(empty1.locate(std::array{4.2}).normDist ==
              std::array<double, 1>{});
        CHECK(empty2.locate(std::array{2., 1.}).anchor ==
              std::array<std::size_t, 2>{});
        CHECK(empty2.locate(std::array{2., 1.}).normDist ==
              std::array<double, 2>{});
    }

    SECTION("Invalid") {
        const auto data = std::vector<int>{1, 2};
        REQUIRE_THROWS_AS(RegularGrid(std::array{std::vector{1.}}, data),
                          std::out_of_range);
        REQUIRE_THROWS_AS(
            RegularGrid(std::array{std::vector{1., 2., 3.}}, data),
            std::out_of_range);
        REQUIRE_THROWS_AS(RegularGrid(std::array{std::vector{2., 1.}}, data),
                          std::runtime_error);
    }

    SECTION("copying and moving") {
        auto copyFromTemp = [] {
            auto data = std::vector{1, 2, 3, 4};
            auto edges = std::array{std::vector{1., 2., 3., 4.}};
            auto grid = RegularGrid{edges, data};
            return std::array{grid};
        };
        auto grid2{copyFromTemp()[0]};
        INFO(grid2.extent());
        CHECK(grid2[{0}] == 1);

        grid2 = grid2;
        CHECK(grid2[{0}] == 1);

        grid2 = std::move(grid2);
        CHECK(grid2[{0}] == 1);
    }
}

TEST_CASE("1D interpolation") {
    using Higgs::utilities::LinearInterpolator;
    const auto X = std::vector<double>{0, 2, 4};
    const auto edge = std::array{X};
    auto interp = LinearInterpolator(edge, std::vector{1., 2., 6.});
    CHECK(interp({0}) == Approx(1));
    CHECK(interp({1}) == Approx(1.5));
    CHECK(interp({1.5}) == Approx(1.75));
    CHECK(interp({2}) == Approx(2));
    CHECK(interp({3.5}) == Approx(5));
    REQUIRE(interp({4}) == Approx(6));

    SECTION("nearest neighbor extrapolation") {
        CHECK(interp({-10}) == interp({0}));
        CHECK(interp({1e6}) == interp({4}));
    }

    SECTION("default constructed") {
        auto interp = utils::LinearInterpolator<1>{};
        REQUIRE(interp({10}) == 0.);
        REQUIRE(interp({-1}) == 0.);
    }

    SECTION("trivial second dimension") {
        auto interp = LinearInterpolator(std::array{X, std::vector<double>{0.}},
                                         std::vector{1., 2., 6.});
        CHECK(interp({0, 0}) == Approx(1));
        CHECK(interp({1, 10}) == Approx(1.5));
        CHECK(interp({1.5, -1.}) == Approx(1.75));
        CHECK(interp({2, 0}) == Approx(2));
        CHECK(interp({3.5, 3.5}) == Approx(5));
        REQUIRE(interp({4, -10}) == Approx(6));
    }
}

TEST_CASE("2d interpolation") {
    const auto X = std::vector{0., 1., 2.};
    const auto Y = std::vector{3., 4., 5., 6.};
    const auto edge = std::array{X, Y};
    const auto values = std::vector{27., 48.,  75., 108., 29., 50.,
                                    77., 110., 43., 64.,  91., 124.};
    auto interp = utils::LinearInterpolator(edge, values);
    CHECK(interp({0, 3}) == Approx(27.0));
    CHECK(interp({1.7, 3.9}) == Approx(57.7));
    CHECK(interp({1.1, 3.0}) == Approx(30.4));
    CHECK(interp({0.2, 5.5}) == Approx(91.9));
    CHECK(interp({0.9, 4.1}) == Approx(52.5));
    CHECK(interp({2, 6.0}) == Approx(124.0));

    SECTION("nearest neighbor extrapolation") {
        auto testvals = std::array{-1e6, -1., 0., 2.5, 5., 10., 200.};
        for (auto t1 : testvals) {
            for (auto t2 : testvals) {
                auto clamped =
                    std::array{std::clamp(t1, ranges::min(X), ranges::max(X)),
                               std::clamp(t2, ranges::min(Y), ranges::max(Y))};
                CHECK(interp({t1, t2}) == interp(clamped));
            }
        }
    }
    SECTION("default constructed") {
        auto interp = utils::LinearInterpolator<2>{};
        REQUIRE(interp({10, -400}) == 0.);
        REQUIRE(interp({-1, 13.8}) == 0.);
    }
}

TEST_CASE("Higher dim interpolation") {
    SECTION("3d") {
        const auto X = std::vector{0., 1., 2.};
        const auto Y = std::vector{3., 4., 5., 6.};
        const auto Z = std::vector{8., 9.};
        const auto edges = std::array{X, Y, Z};

        const auto values = std::vector{
            19., 18., 40.,  39.,  67., 66., 100., 99., 21., 20., 42.,  41.,
            69., 68., 102., 101., 35., 34., 56.,  55., 83., 82., 116., 115.};
        auto interp = utils::LinearInterpolator(edges, values);
        CHECK(interp({0, 3, 8.0}) == Approx(19.0));
        CHECK(interp({1.7, 3.9, 8.1}) == Approx(49.6));
        CHECK(interp({1.1, 3.0, 8.6}) == Approx(21.8));
        CHECK(interp({0.2, 5.5, 8.3}) == Approx(83.6));
        CHECK(interp({0.9, 4.1, 8.9}) == Approx(43.599999999999994));
        CHECK(interp({2, 6.0, 9.0}) == Approx(115.0));
        SECTION("nearest neighbor extrapolation") {
            auto testvals = std::array{-1e6, -1., 0., 2.5, 5., 8.3, 10., 200.};
            for (auto t1 : testvals) {
                for (auto t2 : testvals) {
                    for (auto t3 : testvals) {
                        auto clamped = std::array{
                            std::clamp(t1, ranges::min(X), ranges::max(X)),
                            std::clamp(t2, ranges::min(Y), ranges::max(Y)),
                            std::clamp(t3, ranges::min(Z), ranges::max(Z))};
                        CHECK(interp({t1, t2, t3}) == interp(clamped));
                    }
                }
            }
        }
        SECTION("default constructed") {
            auto interp = utils::LinearInterpolator<3>{};
            REQUIRE(interp({10, -400, 2.}) == 0.);
            REQUIRE(interp({-1, 13.8, 8.}) == 0.);
        }
    }
}

TEST_CASE("Unit cell generation") {
    SECTION("2d square") {
        constexpr auto square = utils::Detail::unitCellCorners<2>();

        INFO(square);
        CHECK(std::is_sorted(square.begin(), square.end()));
        CHECK(square[0][0] == 0);
        CHECK(square[0][1] == 0);
        CHECK(square[1][0] == 0);
        CHECK(square[1][1] == 1);
        CHECK(square[2][0] == 1);
        CHECK(square[2][1] == 0);
        CHECK(square[3][0] == 1);
        CHECK(square[3][1] == 1);
    }
    SECTION("3d cube") {
        constexpr auto cube = utils::Detail::unitCellCorners<3>();
        INFO(cube);
        CHECK(std::is_sorted(cube.begin(), cube.end()));
    }

    SECTION("5d whatever") {
        constexpr auto thing = utils::Detail::unitCellCorners<5>();
        INFO(thing);
        CHECK(std::is_sorted(thing.begin(), thing.end()));
    }
}

struct interptest {
    double x1;
    double x2;

    interptest &operator+=(const interptest other) {
        x1 += other.x1;
        x2 += other.x2;
        return *this;
    }

    friend interptest operator+(interptest a, const interptest &b) {
        a += b;
        return a;
    }

    interptest &operator*=(double d) {
        x1 *= d;
        x2 *= d;
        return *this;
    }

    friend interptest operator*(interptest a, double d) {
        a *= d;
        return a;
    }
};

TEST_CASE("Interpolation with structured values") {
    SECTION("1d grid") {
        const auto X = std::vector<double>{0, 2, 4};
        const auto edge = std::array{X};
        const auto values =
            std::vector<interptest>{{1, 0.1}, {2, 0.2}, {6, 0.6}};
        auto interp = utils::LinearInterpolator(edge, values);
        CHECK(interp({0}).x1 == Approx(1));
        CHECK(interp({1}).x1 == Approx(1.5));
        CHECK(interp({1.5}).x1 == Approx(1.75));
        CHECK(interp({2}).x1 == Approx(2));
        CHECK(interp({3.5}).x1 == Approx(5));
        REQUIRE(interp({4}).x1 == Approx(6));

        CHECK(interp({0}).x2 == Approx(0.1));
        CHECK(interp({1}).x2 == Approx(0.15));
        CHECK(interp({1.5}).x2 == Approx(0.175));
        CHECK(interp({2}).x2 == Approx(0.2));
        CHECK(interp({3.5}).x2 == Approx(0.5));
        REQUIRE(interp({4}).x2 == Approx(0.6));
    }

    SECTION("2d grid") {
        const auto X = std::vector{0., 1., 2.};
        const auto Y = std::vector{3., 4., 5., 6.};
        const auto edge = std::array{X, Y};
        const auto values = std::vector<interptest>{
            {27., 0.27}, {48., 0.48}, {75., 0.75}, {108., 1.08},
            {29., 0.29}, {50., 0.5},  {77., 0.77}, {110., 1.1},
            {43., 0.43}, {64., 0.64}, {91., 0.91}, {124., 1.24}};
        auto interp = utils::LinearInterpolator(edge, values);
        CHECK(interp({0, 3}).x1 == Approx(27.0));
        CHECK(interp({1.7, 3.9}).x1 == Approx(57.7));
        CHECK(interp({1.1, 3.0}).x1 == Approx(30.4));
        CHECK(interp({0.2, 5.5}).x1 == Approx(91.9));
        CHECK(interp({0.9, 4.1}).x1 == Approx(52.5));
        CHECK(interp({2, 6.0}).x1 == Approx(124.0));

        CHECK(interp({0, 3}).x2 == Approx(0.270));
        CHECK(interp({1.7, 3.9}).x2 == Approx(0.577));
        CHECK(interp({1.1, 3.0}).x2 == Approx(0.304));
        CHECK(interp({0.2, 5.5}).x2 == Approx(0.919));
        CHECK(interp({0.9, 4.1}).x2 == Approx(0.525));
        CHECK(interp({2, 6.0}).x2 == Approx(1.240));
    }
}

TEST_CASE("Interpolation with reference to data") {
    const std::vector x = {0., 2., 4.};
    const std::array edges = {x};
    std::vector values = {1., 2., 6.};
    const auto interp =
        utils::LinearInterpolator<1, double, true>{edges, values};
    CHECK(interp({0}) == Approx(1));
    CHECK(interp({1}) == Approx(1.5));
    CHECK(interp({1.5}) == Approx(1.75));
    CHECK(interp({2}) == Approx(2));
    CHECK(interp({3.5}) == Approx(5));
    REQUIRE(interp({4}) == Approx(6));

    SECTION("change data values through reference") {
        INFO("never do this, this is just the easiest way to establish that "
             "the stored references work.");
        values = {2., 4., 12.};
        CHECK(interp({0}) == Approx(2));
        CHECK(interp({1}) == Approx(3.));
        CHECK(interp({1.5}) == Approx(3.5));
        CHECK(interp({2}) == Approx(4.));
        CHECK(interp({3.5}) == Approx(10.));
        REQUIRE(interp({4}) == Approx(12.));
    }
}

TEST_CASE("stacked interpolation") {
    auto i1 = utils::LinearInterpolator{std::array{std::vector{2., 5., 10.}},
                                        std::vector{1., 2., 6.}};
    auto i2 = utils::LinearInterpolator{std::array{std::vector{3., 5., 8.}},
                                        std::vector{3., 4., 1.}};
    auto interp = utils::LinearInterpolator<1, utils::LinearInterpolator<1>>(
        std::array{std::vector{0., 2.}}, std::vector{i1, i2});
    REQUIRE(interp(std::array{2., 4.}) == i2({4.}));
    REQUIRE(interp(std::array{0., 3.}) == i1({3.}));
    REQUIRE(interp(std::array{1., 6.}) == Approx(2.9));
}

TEST_CASE("hypercube interpolation") {
    using Higgs::utilities::LinearInterpolator;

    SECTION("1d") {
        const auto X = std::vector<double>{0, 2, 4};
        const auto edge = std::array{X};
        auto interp = LinearInterpolator(edge, std::vector{3., 2., 6.});

        SECTION("one point within") {
            auto corners =
                interp.atAllCorners({1.5}, {3.}) | ranges::to<std::vector>;
            REQUIRE(corners.size() == 2);
            CHECK(corners[0] == interp({1.5}));
            CHECK(corners[1] == interp({3.}));

            auto interior = interp.atAllPointsBetween({1.5}, {3.}) |
                            ranges::to<std::vector>;
            REQUIRE(interior.size() == 1);
            CHECK(interior[0] == interp({2}));

            CHECK(interp.minWithin({1.5}, {3.}) == interior[0]);
            CHECK(interp.maxWithin({1.5}, {3.}) == corners[1]);
        }

        SECTION("no point within") {
            auto corners =
                interp.atAllCorners({1.8}, {0.3}) | ranges::to<std::vector>;
            REQUIRE(corners.size() == 2);
            CHECK(corners[0] == interp({1.8}));
            CHECK(corners[1] == interp({0.3}));
            auto within = interp.atAllPointsBetween({1.8}, {0.3}) |
                          ranges::to<std::vector>;
            REQUIRE(within.empty());

            CHECK(interp.minWithin({1.8}, {0.3}) == corners[0]);
            CHECK(interp.maxWithin({1.8}, {0.3}) == corners[1]);
        }
    }
    SECTION("2d") {
        const auto X = std::vector{0., 1., 2.};
        const auto Y = std::vector{3., 4., 5., 6.};
        const auto edge = std::array{X, Y};
        const auto values = std::vector{27., 48.,  75., 108., 29., 20.,
                                        77., 110., 43., 64.,  91., 124.};
        auto interp = utils::LinearInterpolator(edge, values);

        auto corners = interp.atAllCorners({1.9, 3.1}, {0.1, 5.8}) |
                       ranges::to<std::vector>;
        INFO(corners);
        REQUIRE(corners.size() == 4);
        CHECK(corners[0] == interp({1.9, 3.1}));
        CHECK(corners[1] == interp({1.9, 5.8}));
        CHECK(corners[2] == interp({0.1, 3.1}));
        CHECK(corners[3] == interp({0.1, 5.8}));
        auto interior = interp.atAllPointsBetween({1.9, 3.1}, {0.1, 5.8}) |
                        ranges::to<std::vector>;
        INFO(interior);
        REQUIRE(interior.size() == 2);
        CHECK(interior[0] == interp({1., 4.}));
        CHECK(interior[1] == interp({1., 5.}));

        CHECK(interp.minWithin({1.9, 3.1}, {0.1, 5.8}) == interior[0]);
        CHECK(interp.maxWithin({1.9, 3.1}, {0.1, 5.8}) == corners[1]);
    }
}
