#include "Higgs/Signals.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "prettyprint.hpp"
#include "signals/JsonSupport.hpp"
#include "utilities/Format.hpp"
#include "utilities/Json.hpp"
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <range/v3/numeric/accumulate.hpp>

using Catch::Approx;
namespace HS = Higgs::signals;
namespace HP = Higgs::predictions;

TEST_CASE("higgs signals operation") {
    auto signals = HS::Signals(TESTMEASUREMENTS_PATH);

    REQUIRE(signals.measurements().size() == 2);
    REQUIRE(
        signals.observableCount() ==
        ranges::accumulate(signals.measurements(), 0, std::plus{},
                           [](const auto &m) { return m.nSubMeasurements(); }));

    auto pred = HP::Predictions();
    auto &h = pred.addParticle(HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
    h.setMass(125);
    HP::effectiveCouplingInput(h, HP::smLikeEffCouplings);

    REQUIRE(signals(pred) == Approx(ranges::accumulate(
                                 signals.measurements(), 0., std::plus{},
                                 [&pred](const auto &m) { return m(pred); })));
}
