#include "utilities/ArithmeticArray.hpp"
#include <catch2/catch_test_macros.hpp>
#include <memory>

using Higgs::utilities::ArithmeticArray;

TEST_CASE("arithmetic array") {
    const auto x = ArithmeticArray{10, 20, 30};
    const auto y = ArithmeticArray{1, 2, 3};

    SECTION("equality comparisons") {
        REQUIRE(x == x);
        REQUIRE(y == y);
        REQUIRE_FALSE(x != x);
        REQUIRE_FALSE(y != y);
        REQUIRE(x != y);
        REQUIRE_FALSE(x == y);
    }

    SECTION("unary -") { REQUIRE(-x == std::array{-10, -20, -30}); }

    SECTION("compound assignments -- arrays") {
        auto x1 = x;
        x1 += x;
        REQUIRE(x1 == std::array{20, 40, 60});
        x1 -= y;
        REQUIRE(x1 == std::array{19, 38, 57});
        x1 *= y;
        REQUIRE(x1 == std::array{19, 76, 171});
        x1 /= x;
        REQUIRE(x1 == std::array{1, 3, 5});
    }

    SECTION("binary arithmetics -- arrays") {
        REQUIRE(x + y == std::array{11, 22, 33});
        REQUIRE(x - y == std::array{9, 18, 27});
        REQUIRE(x * y == std::array{10, 40, 90});
        REQUIRE(x / y == std::array{10, 10, 10});
    }

    SECTION("compound assignments -- scalars") {
        auto x1 = x;
        x1 *= 2;
        REQUIRE(x1 == x + x);
        x1 /= 2;
        REQUIRE(x1 == x);
        x1 += 3;
        REQUIRE(x1 == std::array{13, 23, 33});
        x1 -= 3;
        REQUIRE(x1 == x);
    }

    SECTION("binary arithmetics -- scalars") {
        CHECK(y + 10 == std::array{11, 12, 13});
        CHECK(10 + x == x + 10);
        CHECK(x - 5 == std::array{5, 15, 25});
        CHECK(5 - y == -(y - 5));
        CHECK(y * 2 == y + y);
        CHECK(2 * y == y * 2);
        CHECK(x / 5 == std::array{2, 4, 6});
        REQUIRE(30 / y == std::array{30, 15, 10});
    }

    SECTION("structured binding") {
        auto [x1, x2, x3] = x;
        REQUIRE(x1 == 10);
        REQUIRE(x2 == 20);
        REQUIRE(x3 == 30);

        auto xx = x;
        auto &[xx1, xx2, xx3] = xx;
        REQUIRE(xx1 == 10);
        REQUIRE(xx2 == 20);
        REQUIRE(xx3 == 30);

        auto getX = [x]() { return x; };
        auto [xxx1, xxx2, xxx3] = getX();
        REQUIRE(xxx1 == 10);
        REQUIRE(xxx2 == 20);
        REQUIRE(xxx3 == 30);
    }

    SECTION("and some checks for doubles") {
        const auto d = ArithmeticArray{1., 1.};
        REQUIRE(d + d - 1. == d);
        REQUIRE(d + 3 == 4 * d);
        REQUIRE(d / (2 * d) == d * 0.5 * d);
    }
}
