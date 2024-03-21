#include "predictions/UncertainMass.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/repeat_n.hpp>

using Catch::Approx;
namespace HP = Higgs::predictions;

TEST_CASE("mass uncertainty") {
    const auto m1 = HP::UncertainMass{100, 10};
    const auto m2 = HP::UncertainMass{200, 5};
    const auto zero = HP::UncertainMass{0, 0};

    SECTION("equality") {
        CHECK(m1 == m1);
        CHECK_FALSE(m1 != m1);
        CHECK_FALSE(m1 == m2);
        CHECK(m1 != m2);
    }

    SECTION("addition") {
        CHECK(m1 + zero == m1);
        auto mt = m1 + m2;
        CHECK(mt.mass == m1.mass + m2.mass);
        CHECK(2 * mt.uncertainty ==
              Approx((m2.mass + m2.uncertainty + m1.mass + m1.uncertainty) -
                     (m2.mass - m2.uncertainty + m1.mass - m1.uncertainty)));
    }

    SECTION("scalar multiplication") {
        CHECK(m1 * 0 == zero);
        CHECK(0 * m1 == zero);

        auto mt = m1 * 3;
        CHECK(mt.mass == 3 * m1.mass);
        CHECK(mt.uncertainty == 3 * m1.uncertainty);
        mt = -1.4 * m1;
        CHECK(mt.mass == -1.4 * m1.mass);
        CHECK(mt.uncertainty == 1.4 * m1.uncertainty);
    }

    SECTION("scalar division") {
        CHECK(m1 / std::numeric_limits<double>::infinity() == zero);
        auto mt = m1 / 3.;
        CHECK(mt.mass == m1.mass / 3.);
        CHECK(mt.uncertainty == m1.uncertainty / 3.);
        CHECK(m1 * 0.5 == m1 / 2.);
        CHECK(m1 / -1. == m1 * -1.);
    }

    SECTION("unary -") { CHECK(-m1 == m1 * -1); }

    SECTION("no-ops") {
        CHECK(m1 + m1 == 2 * m1);
        CHECK((m1 + m1) / 2 == m1);
        CHECK((2 * m1 + 3 * m1) / 5 == m1);
    }

    SECTION("accumulate") {
        CHECK(ranges::accumulate(ranges::views::repeat_n(m1, 2), m2) ==
              m1 + m1 + m2);
    }

    SECTION("lies within") {
        CHECK(m1.liesWithin({0., 200.}, HP::MassUncEagerness::eager));
        CHECK(m1.liesWithin({0., 200.}, HP::MassUncEagerness::cautious));
        CHECK(m1.liesWithin({0., 200.}, HP::MassUncEagerness::ignore));

        CHECK_FALSE(m1.liesWithin({0., 0.}, HP::MassUncEagerness::eager));
        CHECK_FALSE(m1.liesWithin({0., 0.}, HP::MassUncEagerness::cautious));
        CHECK_FALSE(m1.liesWithin({0., 0.}, HP::MassUncEagerness::ignore));

        CHECK(m1.liesWithin({95., 105.}, HP::MassUncEagerness::eager));
        CHECK(m1.liesWithin({95., 105.}, HP::MassUncEagerness::ignore));
        CHECK_FALSE(m1.liesWithin({95., 105.}, HP::MassUncEagerness::cautious));

        CHECK(m1.liesWithin({104., 105.}, HP::MassUncEagerness::eager));
        CHECK_FALSE(m1.liesWithin({104., 105.}, HP::MassUncEagerness::ignore));
        CHECK_FALSE(
            m1.liesWithin({104., 105.}, HP::MassUncEagerness::cautious));

        SECTION("invalid eagerness") {
            CHECK_FALSE(m1.liesWithin({0., 200.}, HP::MassUncEagerness{5}));
        }
    }

    SECTION("formatting") {
        CHECK_THAT(
            fmt::format("{}", m1),
            Catch::Matchers::ContainsSubstring(fmt::format("{}", m1.mass)) &&
                Catch::Matchers::ContainsSubstring(
                    fmt::format("{}", m1.uncertainty)));
    }
}
