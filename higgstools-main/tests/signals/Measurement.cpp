#include "Higgs/signals/Measurement.hpp"
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
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/map.hpp>

using Catch::Approx;

namespace HP = Higgs::predictions;
namespace HS = Higgs::signals;

TEST_CASE("Measurement read and props") {

    REQUIRE_NOTHROW(HS::Measurement(TESTMEASUREMENTS_PATH "STXS_test.json"));
    auto meas = HS::Measurement(TESTMEASUREMENTS_PATH "STXS_test.json");
    REQUIRE(meas.id() == 200403447);
    REQUIRE(meas.reference() == "2004.03447");
    REQUIRE(meas.citeKey() == "ATLAS:2020rej");
    REQUIRE(meas.collider() == HP::Collider::LHC13);
    REQUIRE(meas.experiment() == HP::Experiment::ATLAS);
    REQUIRE(meas.luminosity() == Approx(139));
    REQUIRE(meas.referenceMass() == Approx(125));
    REQUIRE(meas.loadedFrom() == TESTMEASUREMENTS_PATH "STXS_test.json");
    REQUIRE(meas.nSubMeasurements() == 12);
    REQUIRE(meas.referenceModel() == HP::ReferenceModel::SMHiggs);

    auto measurements = meas.subMeasurements();
    CHECK(measurements.find("gg2H_0j_pTH_High") != measurements.end());
    CHECK(measurements.find("VH_Lep") != measurements.end());

    SECTION("options") {
        REQUIRE(meas.options().rescaleToRefMass ==
                HS::MeasurementOptions{}.rescaleToRefMass);
        auto opts = HS::MeasurementOptions{};
        opts.rescaleToRefMass = HS::RescaleToRefMass::always;
        opts.ignoreTheoryUncertainties = true;
        auto measOpt =
            HS::Measurement{TESTMEASUREMENTS_PATH "STXS_test.json", opts};
        REQUIRE(measOpt.options().rescaleToRefMass == opts.rescaleToRefMass);
    }

    SECTION("copy and move operations") {
        auto measCopy = meas;
        REQUIRE(measCopy.id() == 200403447);

        auto measCopy2 = std::move(measCopy);
        REQUIRE(measCopy2.id() == 200403447);

        measCopy =
            HS::Measurement(TESTMEASUREMENTS_PATH "run1comb_as_in_HS2.json");
        REQUIRE(measCopy.id() == 160602266);

        measCopy = meas;
        REQUIRE(measCopy.id() == 200403447);
    }

    SECTION("evaluation") {
        auto opts = HS::MeasurementOptions{};
        opts.rescaleToRefMass = HS::RescaleToRefMass::withinMassUnc;
        auto meas = HS::Measurement(TESTMEASUREMENTS_PATH "STXS_test.json", opts);
        auto pred = Higgs::Predictions{};
        auto &h = pred.addParticle(
            HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
        h.setMass(125.09);
        HP::effectiveCouplingInput(h, HP::smLikeEffCouplings,
                                   HP::ReferenceModel::SMHiggs);
        REQUIRE(meas(pred, {}) == Approx(7.207238518));
    }
}

TEST_CASE("individual contributions and modification factors") {
    auto options = HS::MeasurementOptions{};
    options.whichCorrelations = HS::Correlations::none;
    auto meas =
        HS::Measurement(TESTMEASUREMENTS_PATH "STXS_test.json", options);

    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
    h.setMass(125.09);
    HP::effectiveCouplingInput(h, HP::smLikeEffCouplings);

    auto contribs = meas.chisqContributions(pred);

    INFO("when correlations are ignored the sum of the contributions is the "
         "total chisq");
    REQUIRE(ranges::accumulate(contribs | ranges::views::values, 0.) ==
            Approx(meas(pred)));

    SECTION("non-trivial modification factors") {
        INFO(contribs);
        auto modFactors = HS::Measurement::ModificationFactors{};
        modFactors["gg2H_pTH_High"] = {0.8};
        modFactors["gg2H_2j"] = {0.5};

        auto modContribs = meas.chisqContributions(pred, modFactors);
        REQUIRE(ranges::accumulate(modContribs | ranges::views::values, 0.) ==
                Approx(meas(pred, modFactors)));

        for (auto key : ranges::views::keys(contribs)) {
            if (key != "gg2H_pTH_High" && key != "gg2H_2j") {
                REQUIRE(modContribs.at(key) == Approx(contribs.at(key)));
            }
        }
        CHECK(modContribs.at("gg2H_pTH_High") > contribs.at("gg2H_pTH_High"));
        CHECK(modContribs.at("gg2H_2j") < contribs.at("gg2H_2j"));
    }
}
