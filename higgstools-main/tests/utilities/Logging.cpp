#include "utilities/Logging.hpp"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("getting loggers") {
    SECTION("predictions") {
        auto log = Higgs::predictions::logger();
        REQUIRE_NOTHROW(log->warn("testing the logger"));
    }
    SECTION("bounds") {
        auto log = Higgs::bounds::logger();
        REQUIRE_NOTHROW(log->warn("testing the logger"));
    }
    SECTION("signals") {
        auto log = Higgs::signals::logger();
        REQUIRE_NOTHROW(log->warn("testing the logger"));
    }
}
