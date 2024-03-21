#include "signals/Measurements/RateMeasurement.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "signals/JsonSupport.hpp"
#include "utilities/Format.hpp"
#include "utilities/Json.hpp"
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/map.hpp>

using Catch::Approx;

namespace HP = Higgs::predictions;
namespace HS = Higgs::signals;

namespace {
template <class SubMeasurementClass, class... Args>
SubMeasurementClass readToValue(Args &&...args) {
    auto ptr = std::shared_ptr<HS::SubMeasurement>{};
    if constexpr (std::is_same<SubMeasurementClass, HS::RateMeasurement>{}) {
        ptr = HS::readRateMeasurement(std::forward<Args>(args)...);
    } else if constexpr (std::is_same<SubMeasurementClass,
                                      HS::MassMeasurement>{}) {
        ptr = HS::readMassMeasurement(std::forward<Args>(args)...);
    }
    return *std::dynamic_pointer_cast<SubMeasurementClass>(ptr);
}
} // namespace

TEST_CASE("Rate measurement") {
    const auto proc = HP::ChannelProcess{
        HP::Collider::LHC13, {{HP::Production::ggH, HP::Decay::ZZ}}};
    const auto obsVal = HS::UncertainValue{0.9, 1.0, 1.1};
    const auto refVal = HS::UncertainValue{0.4, 0.5, 0.6};

    auto rateMeas = HS::RateMeasurement{
        proc, {0.5},  obsVal, std::make_shared<HP::SMHiggs>(125),
        2.5,  refVal, true};
    REQUIRE(rateMeas.observedRate() == obsVal.central());
    REQUIRE(rateMeas.observedRateUnc(HS::Uncertainty::minus) ==
            Approx(obsVal.lowerUncertainty()));
    REQUIRE(rateMeas.observedRateUnc(HS::Uncertainty::plus) ==
            Approx(obsVal.upperUncertainty()));
    REQUIRE(rateMeas.referenceRate() == refVal.central());
    REQUIRE(rateMeas.referenceRateUnc(HS::Uncertainty::minus) ==
            Approx(refVal.lowerUncertainty()));
    REQUIRE(rateMeas.referenceRateUnc(HS::Uncertainty::plus) ==
            Approx(refVal.upperUncertainty()));
    REQUIRE(rateMeas.referenceMass() == Approx(125));
    REQUIRE(rateMeas.massResolution() == Approx(2.5));
    REQUIRE(rateMeas.massSensitive());
    REQUIRE_THAT(rateMeas.processDesc(),
                 Catch::Matchers::ContainsSubstring("ggH") &&
                     Catch::Matchers::ContainsSubstring("LHC13") &&
                     Catch::Matchers::ContainsSubstring("ZZ"));

    SECTION("json read") {
        auto j = nlohmann::json::parse(R"(
            {
                "process": {
                    "channels": [["ggH","ZZ"]]
                },
                "obs": [0.9, 1.0, 1.1],
                "exp": [0.4, 0.5, 0.60],
                "channelWeights": [0.97],
                "massResolution": 2.5,
                "massSensitive": true
            }
        )");
        REQUIRE_NOTHROW(HS::readRateMeasurement(
            j, HP::Collider::LHC8, std::make_shared<HP::SMHiggs>(200), 1.5));
        auto meas = readToValue<HS::RateMeasurement>(
            j, HP::Collider::LHC13, std::make_shared<HP::SMHiggs>(125), 1.5);
        REQUIRE(rateMeas.observedRate() == obsVal.central());
        REQUIRE(rateMeas.observedRateUnc(HS::Uncertainty::minus) ==
                Approx(obsVal.lowerUncertainty()));
        REQUIRE(rateMeas.observedRateUnc(HS::Uncertainty::plus) ==
                Approx(obsVal.upperUncertainty()));
        REQUIRE(rateMeas.referenceRate() == refVal.central());
        REQUIRE(rateMeas.referenceRateUnc(HS::Uncertainty::minus) ==
                Approx(refVal.lowerUncertainty()));
        REQUIRE(rateMeas.referenceRateUnc(HS::Uncertainty::plus) ==
                Approx(refVal.upperUncertainty()));
        REQUIRE(rateMeas.referenceMass() == Approx(125));
        REQUIRE(rateMeas.massResolution() == Approx(2.5));
        REQUIRE(rateMeas.massSensitive() == true);
    }

    SECTION("invalid construction") {
        REQUIRE_THROWS_AS(HS::RateMeasurement(
                              proc, {1., 1.}, obsVal,
                              std::make_shared<HP::SMHiggs>(125), 2.5, refVal),
                          HS::InvalidMeasurement);
        REQUIRE_THROWS_WITH(
            HS::RateMeasurement(proc, {1., 1.}, obsVal,
                                std::make_shared<HP::SMHiggs>(125), 2.5,
                                refVal),
            Catch::Matchers::ContainsSubstring("not match"));
    }

    SECTION("signal strength evaluation") {
        auto p = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
        p.setMass(120);
        p.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                         HP::Decay::ZZ, 10);
        p.setChannelRate(HP::Collider::LHC13, HP::Production::bbH,
                         HP::Decay::WW, 5);
        const auto proc2 =
            HP::ChannelProcess{HP::Collider::LHC13,
                               {{HP::Production::ggH, HP::Decay::ZZ},
                                {HP::Production::bbH, HP::Decay::WW}}};
        auto rateMeas2 = HS::RateMeasurement{
            proc2, {1.1, 0.8},     obsVal, std::make_shared<HP::SMHiggs>(125),
            2.5,   {10., 20., 30.}};

        SECTION("normalization") {
            auto refMass = HP::SMHiggs{rateMeas.referenceMass()};
            auto ref = HP::SMHiggs{p.mass()};

            auto refRate = ref.channelRate(HP::Collider::LHC13,
                                           HP::Production::ggH, HP::Decay::ZZ);
            auto refMassRate = refMass.channelRate(
                HP::Collider::LHC13, HP::Production::ggH, HP::Decay::ZZ);
            auto refRate2 = ref.channelRate(HP::Collider::LHC13,
                                            HP::Production::bbH, HP::Decay::WW);
            auto refMassRate2 = refMass.channelRate(
                HP::Collider::LHC13, HP::Production::bbH, HP::Decay::WW);
            CHECK(
                rateMeas.signalStrength(p, {}, HS::NormalizeAt::particleMass) ==
                Approx(10 / refRate));
            CHECK(rateMeas.signalStrength(p, {},
                                          HS::NormalizeAt::referenceMass) ==
                  Approx(10 / refMassRate));
            CHECK(rateMeas2.signalStrength(p, {},
                                           HS::NormalizeAt::particleMass) ==
                  Approx((10 / refRate * 1.1 * refMassRate +
                          5 / refRate2 * 0.8 * refMassRate2) /
                         (1.1 * refMassRate + 0.8 * refMassRate2)));

            CHECK(rateMeas2.signalStrength(p, {},
                                           HS::NormalizeAt::referenceMass) ==
                  Approx((10 * 1.1 + 5 * 0.8) /
                         (1.1 * refMassRate + 0.8 * refMassRate2)));
        }

        SECTION("modification factors") {
            for (auto norm : magic_enum::enum_values<HS::NormalizeAt>()) {
                for (double m : {-2., 0.5, 1., 3., 100.}) {
                    CHECK(rateMeas.signalStrength(p, {m, 1e3, 1e5}, norm) ==
                          Approx(m * rateMeas.signalStrength(p, {}, norm)));
                    CHECK(rateMeas.signalStrength(p, {m, 0., 1e8}, norm) ==
                          Approx(m * rateMeas.signalStrength(p, {}, norm)));
                    CHECK(rateMeas2.signalStrength(p, {m, m, -1e4}, norm) ==
                          Approx(m * rateMeas2.signalStrength(p, {}, norm)));
                    CHECK(rateMeas2.signalStrength(p, {m, m, 100, 100}, norm) ==
                          Approx(m * rateMeas2.signalStrength(p, {}, norm)));
                }
                auto pModRates = p;
                pModRates.setChannelRate(
                    HP::Collider::LHC13, HP::Production::ggH, HP::Decay::ZZ,
                    0.5 * p.channelRate(HP::Collider::LHC13,
                                        HP::Production::ggH, HP::Decay::ZZ));
                pModRates.setChannelRate(
                    HP::Collider::LHC13, HP::Production::bbH, HP::Decay::WW,
                    3 * p.channelRate(HP::Collider::LHC13, HP::Production::bbH,
                                      HP::Decay::WW));
                CHECK(rateMeas2.signalStrength(p, {0.5, 3.}, norm) ==
                      rateMeas2.signalStrength(pModRates, {}, norm));

                pModRates.setChannelRate(
                    HP::Collider::LHC13, HP::Production::ggH, HP::Decay::ZZ,
                    2 * p.channelRate(HP::Collider::LHC13, HP::Production::ggH,
                                      HP::Decay::ZZ));
                pModRates.setChannelRate(
                    HP::Collider::LHC13, HP::Production::bbH, HP::Decay::WW,
                    1.4 * p.channelRate(HP::Collider::LHC13,
                                        HP::Production::bbH, HP::Decay::WW));
                CHECK(rateMeas2.signalStrength(p, {2., 1.4}, norm) ==
                      Approx(rateMeas2.signalStrength(pModRates, {}, norm)));
            }
        }
    }
}

namespace {
bool particleSetContains(const HP::ParticleSet &set, const HP::Particle &elem) {
    return set.find(std::cref(elem)) != set.end();
}
} // namespace

TEST_CASE("assign particles for rate measurement") {
    auto pred = Higgs::Predictions{};

    auto &h1 = pred.addParticle(
        HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
    auto &h2 = pred.addParticle(
        HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
    auto &h3 = pred.addParticle(
        HP::BsmParticle("h3", HP::ECharge::neutral, HP::CP::undefined));
    auto &h4 = pred.addParticle(
        HP::BsmParticle("h4", HP::ECharge::neutral, HP::CP::undefined));

    h1.setMass(100);
    h1.setMassUnc(10);
    h2.setMass(90);
    h2.setMassUnc(1.);
    h3.setMass(110);
    h3.setMassUnc(4);
    h4.setMass(85);
    h4.setMassUnc(20);

    h1.setChannelRate(HP::Collider::LHC13, HP::Production::ggH, HP::Decay::bb,
                      1.);
    h2.setChannelRate(HP::Collider::LHC13, HP::Production::ggH, HP::Decay::bb,
                      1.);
    h3.setChannelRate(HP::Collider::LHC13, HP::Production::ggH, HP::Decay::bb,
                      1.);
    h4.setChannelRate(HP::Collider::LHC13, HP::Production::ggH, HP::Decay::bb,
                      1.);

    auto meas = HS::RateMeasurement{
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::ggH, HP::Decay::bb}}},
        {1.},
        {1., 1., 1.},
        std::make_shared<HP::SMHiggs>(100.),
        4.,
        {1., 1., 1.}};

    SECTION("box pdf for THU") {
        auto cluster01 = meas.assignParticles(pred, HS::PDF::box, 0.1);
        CHECK(cluster01.size() == 2);
        CHECK(particleSetContains(cluster01, h1));
        CHECK_FALSE(particleSetContains(cluster01, h2));
        CHECK_FALSE(particleSetContains(cluster01, h3));
        CHECK(particleSetContains(cluster01, h4));

        auto cluster1 = meas.assignParticles(pred, HS::PDF::box, 1.);
        CHECK(cluster1.size() == 2);
        CHECK(particleSetContains(cluster1, h1));
        CHECK_FALSE(particleSetContains(cluster1, h2));
        CHECK_FALSE(particleSetContains(cluster1, h3));
        CHECK(particleSetContains(cluster1, h4));

        auto cluster2 = meas.assignParticles(pred, HS::PDF::box, 2);
        CHECK(cluster2.size() == 3);
        CHECK(particleSetContains(cluster2, h1));
        CHECK_FALSE(particleSetContains(cluster2, h2));
        CHECK(particleSetContains(cluster2, h3));
        CHECK(particleSetContains(cluster2, h4));

        auto cluster23 = meas.assignParticles(pred, HS::PDF::box, 2.3);
        CHECK(cluster23.size() == 4);
        CHECK(particleSetContains(cluster23, h1));
        CHECK(particleSetContains(cluster23, h2));
        CHECK(particleSetContains(cluster23, h3));
        CHECK(particleSetContains(cluster23, h4));
    }

    SECTION("gaussian pdf for THU") {
        auto cluster01 = meas.assignParticles(pred, HS::PDF::gaussian, 0.1);
        CHECK(cluster01.size() == 1);
        CHECK(particleSetContains(cluster01, h1));
        CHECK_FALSE(particleSetContains(cluster01, h2));
        CHECK_FALSE(particleSetContains(cluster01, h3));
        CHECK_FALSE(particleSetContains(cluster01, h4));

        auto cluster1 = meas.assignParticles(pred, HS::PDF::gaussian, 1.);
        CHECK(cluster1.size() == 2);
        CHECK(particleSetContains(cluster1, h1));
        CHECK_FALSE(particleSetContains(cluster1, h2));
        CHECK_FALSE(particleSetContains(cluster1, h3));
        CHECK(particleSetContains(cluster1, h4));

        auto cluster23 = meas.assignParticles(pred, HS::PDF::gaussian, 2.3);
        CHECK(cluster23.size() == 3);
        CHECK(particleSetContains(cluster23, h1));
        CHECK_FALSE(particleSetContains(cluster23, h2));
        CHECK(particleSetContains(cluster23, h3));
        CHECK(particleSetContains(cluster23, h4));

        auto cluster25 = meas.assignParticles(pred, HS::PDF::gaussian, 2.5);
        CHECK(cluster25.size() == 4);
        CHECK(particleSetContains(cluster25, h1));
        CHECK(particleSetContains(cluster25, h2));
        CHECK(particleSetContains(cluster25, h3));
        CHECK(particleSetContains(cluster25, h4));
    }
}

TEST_CASE("evaluate rate measurement") {
    auto pred = Higgs::Predictions{};

    auto &h1 = pred.addParticle(
        HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
    auto &h2 = pred.addParticle(
        HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
    auto &h3 = pred.addParticle(
        HP::BsmParticle("h3", HP::ECharge::neutral, HP::CP::undefined));
    auto &h4 = pred.addParticle(
        HP::BsmParticle("h4", HP::ECharge::neutral, HP::CP::undefined));

    h1.setMass(125);
    h1.setMassUnc(10);
    h2.setMass(120);
    h2.setMassUnc(1.);
    h3.setMass(130);
    h3.setMassUnc(4);
    h4.setMass(110);
    h4.setMassUnc(20);

    HP::effectiveCouplingInput(h1, HP::scaledSMlikeEffCouplings(std::sqrt(0.5)),
                               HP::ReferenceModel::SMHiggs);
    auto meas = HS::RateMeasurement{
        HP::ChannelProcess{HP::Collider::LHC13,
                           {{HP::Production::H, HP::Decay::ZZ}}},
        {1.},
        HS::UncertainValue::fromInterval(0.115, 0.170, 0.225),
        std::make_shared<HP::SMHiggs>(125.),
        2.5,
        HS::UncertainValue::fromInterval(0.151, 0.176, 0.201)};

    REQUIRE(meas.signalStrength(h1, {}, HS::NormalizeAt::referenceMass) ==
            Approx(0.5));
    REQUIRE(meas.referenceRate() == 0.176);
    REQUIRE(meas.observedRate() == 0.17);

    auto opts = HS::MeasurementOptions{};
    opts.rescaleToRefMass =
        HS::RescaleToRefMass::always; // that's the HS-2 behaviour
    SECTION("h1 only") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        CHECK(residual == Approx(0.170 - 0.176 / 2.));
        CHECK(obsVariance == Approx(5.5e-2));
        CHECK(refVariance == Approx(0.5 * 2.5e-2));
        CHECK(extraChisq == 0.);

        CHECK(meas.chisq(pred, {}, opts) ==
              Approx(std::pow(residual, 2) /
                         (std::pow(obsVariance, 2) + std::pow(refVariance, 2)) +
                     extraChisq));
        auto modOpts = opts;
        modOpts.ignoreTheoryUncertainties = true;
        CHECK(meas.chisq(pred, {}, modOpts) ==
              Approx(std::pow(residual / obsVariance, 2) + extraChisq));
    }

    HP::effectiveCouplingInput(h2, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs);
    SECTION("h1 & h2, h2 not assigned") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        CHECK(residual == Approx(0.170 - 0.176 / 2.));
        CHECK(obsVariance == Approx(5.5e-2));
        CHECK(refVariance == Approx(0.5 * 2.5e-2));
        CHECK(extraChisq == 0.);
    }

    HP::effectiveCouplingInput(h3, HP::smLikeEffCouplings,
                               HP::ReferenceModel::SMHiggs);
    SECTION("h1,2,3") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        CHECK(residual == Approx(-9.4e-2)); // numbers from HS-2
        CHECK(obsVariance == Approx(5.5e-2));
        CHECK(refVariance == Approx(3.75e-2));
        CHECK(extraChisq == 0.);
    }

    HP::effectiveCouplingInput(h4, HP::scaledSMlikeEffCouplings(0.1),
                               HP::ReferenceModel::SMHiggs);
    SECTION("h1,2,3,4") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        CHECK(residual == Approx(-0.09576)); // numbers from HS-2
        CHECK(obsVariance == Approx(5.5e-2));
        CHECK(refVariance == Approx(3.775e-2));
        CHECK(extraChisq == 0.);
    }
}

TEST_CASE("mass measurement") {
    const auto proc = HP::ChannelProcess{
        HP::Collider::LHC13, {{HP::Production::ggH, HP::Decay::gamgam}}};

    const auto meas =
        HS::MassMeasurement(proc, {1.}, HS::UncertainValue{125.09, 0.24},
                            std::make_shared<HP::SMHiggs>(125.09));

    SECTION("json read") {
        auto j = nlohmann::json::parse(R"(
            {
                "process": {
                    "channels": [["ggH","ZZ"]]
                },
                "obsMass": [122,125,130],
                "channelWeights":[1.1]
            }
        )");
        auto readMeas = HS::readMassMeasurement(
            j, HP::Collider::LHC13, std::make_shared<HP::SMHiggs>(125.));
        REQUIRE(dynamic_cast<HS::MassMeasurement *>(readMeas.get())
                    ->observedMass() == Approx(125));
        REQUIRE(dynamic_cast<HS::MassMeasurement *>(readMeas.get())
                    ->observedMassUnc(HS::Uncertainty::plus) == Approx(5));
    }

    SECTION("single particle") {
        auto pred = HP::Predictions{};
        auto &h = pred.addParticle(
            HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
        h.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                         HP::Decay::gamgam, 1.);

        SECTION("at observed mass, irrelevant ThU") {
            h.setMass(meas.observedMass());
            auto options = HS::MeasurementOptions{};
            options.theoryMassUncPDF =
                GENERATE(HS::PDF::box, HS::PDF::gaussian);
            auto theoUnc = GENERATE(0., 2.5);
            auto [residual, expVar, theoVar, explicitChisq] =
                meas.evaluate(pred, {}, options);
            CHECK(residual == 0.);
            CHECK(expVar == meas.observedMassUnc(HS::Uncertainty::plus));
            CHECK((theoVar == 0. || theoVar == theoUnc));
            CHECK(explicitChisq == 0.);
        }

        SECTION("THU PDFs are irrelevant for 0 THU") {
            h.setMassUnc(0.);

            auto options = HS::MeasurementOptions{};
            options.theoryMassUncPDF =
                GENERATE(HS::PDF::gaussian, HS::PDF::box);
            SECTION("far off, unassigned and clamped with penalty") {
                auto m = GENERATE(0, 100, 143.2, 1000);
                options.massSensitiveAssignmentRange = GENERATE(1, 2, 3);
                options.unassignedMassMeasurementPenalty = GENERATE(1, 2, 3);
                h.setMass(m);
                auto result = meas.evaluate(pred, {}, options);
                REQUIRE(result.refVariance == 0.);
                REQUIRE(result.residual == 0.);
                REQUIRE(
                    result.extraChisq ==
                    Approx(std::pow(options.massSensitiveAssignmentRange, 2) *
                           options.unassignedMassMeasurementPenalty));
                REQUIRE(meas.chisq(pred, {}, options) == result.extraChisq);
            }

            SECTION("assigned in range") {
                auto m = GENERATE(124.7, 125., 125.4, 125.65);
                h.setMass(m);
                auto result = meas.evaluate(pred, {}, options);
                REQUIRE(result.refVariance == 0.);
                REQUIRE(result.residual == Approx(meas.observedMass() - m));
                REQUIRE(result.obsVariance ==
                        (m >= meas.observedMass()
                             ? meas.observedMassUnc(HS::Uncertainty::plus)
                             : meas.observedMassUnc(HS::Uncertainty::minus)));
                REQUIRE(result.extraChisq == 0.);
            }

            SECTION("continuity at border of assignment") {
                options.unassignedMassMeasurementPenalty = 1;
                auto m = meas.observedMass() -
                         options.massSensitiveAssignmentRange *
                             meas.observedMassUnc(HS::Uncertainty::minus);
                h.setMass(m + 1e-6);
                auto chisqAssigned = meas.chisq(pred, {}, options);
                h.setMass(m - 1e-6);
                auto chisqUnassigned = meas.chisq(pred, {}, options);
                CHECK(Approx(chisqAssigned).epsilon(1e-4) == chisqUnassigned);
            }
        }

        SECTION("box THU test values") {
            h.setMassUnc(1.);

            auto options = HS::MeasurementOptions{};
            options.theoryMassUncPDF = HS::PDF::box;
            INFO("compare values to HS2");
            for (auto [m, chisq] : std::vector<std::pair<double, double>>{
                     {123., 18.},
                     {123.5, 6.0434027777778478},
                     {124., 0.14062500000001066},
                     {124.5, 0.},
                     {125., 0.},
                     {125.5, 0.},
                     {126., 0.},
                     {126.5, 2.9184027777777293},
                     {127., 18.}}) {
                INFO(m);
                h.setMass(m);
                CHECK(meas.chisq(pred, {}, options) == Approx(chisq));
            }
        }

        SECTION("gaussian THU test values") {
            h.setMassUnc(1.);

            auto options = HS::MeasurementOptions{};
            options.theoryMassUncPDF = HS::PDF::gaussian;
            INFO("compare values to HS2");
            for (auto [m, chisq] : std::vector<std::pair<double, double>>{
                     {123., 4.1302004538578041},
                     {123.5, 2.3904122541603732},
                     {124., 1.1233925869894170},
                     {124.5, 0.32914145234493569},
                     {125., 7.6588502269294752E-003},
                     {125.5, 0.15894478063539824},
                     {126., 0.78299924357034212},
                     {126.5, 1.8798222390317609},
                     {127., 3.4494137670197058}}) {
                INFO(m);
                h.setMass(m);
                CHECK(meas.chisq(pred, {}, options) == Approx(chisq));
            }
        }
    }

    SECTION("two particles") {

        auto pred = HP::Predictions{};
        auto &h1 = pred.addParticle(
            HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
        h1.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                          HP::Decay::gamgam, 2.);
        auto &h2 = pred.addParticle(
            HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
        h2.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                          HP::Decay::gamgam, 3.);

        h1.setMassUnc(1.);
        h2.setMassUnc(1.);

        SECTION("gaussian THU test values") {

            auto options = HS::MeasurementOptions{};
            // avoids mass-dependence of the reference rates
            options.rescaleToRefMass = HS::RescaleToRefMass::never;
            options.theoryMassUncPDF = HS::PDF::gaussian;
            SECTION("one particle-equivalent cases") {
                for (auto [m, chisq] : std::vector<std::pair<double, double>>{
                         {123., 4.1302004538578041},
                         {123.5, 2.3904122541603732},
                         {124., 1.1233925869894170},
                         {124.5, 0.32914145234493569},
                         {125., 7.6588502269294752E-003},
                         {125.5, 0.15894478063539824},
                         {126., 0.78299924357034212},
                         {126.5, 1.8798222390317609},
                         {127., 3.4494137670197058}}) {
                    h1.setMass(m);
                    for (auto m2 : {121., m, 129.}) {
                        h2.setMass(m2);
                        INFO(m << " " << m2);
                        CHECK(meas.chisq(pred, {}, options) == Approx(chisq));
                    }
                }
            }

            h2.setMassUnc(0.5);

            SECTION("other test values") {
                auto predEff = HP::Predictions{};
                auto &heff = predEff.addParticle(HP::BsmParticle(
                    "heff", HP::ECharge::neutral, HP::CP::undefined));
                heff.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                                    HP::Decay::gamgam, 1.);
                heff.setMassUnc(0.4 * h1.massUnc() + 0.6 * h2.massUnc());

                for (auto [m1, m2, chisq, sepChisq] :
                     std::vector<std::array<double, 4>>{
                         {123, 124, 4.497, 0.444},
                         {123, 126.5, 5.334, 5.327},
                         {124.5, 125.5, 0.44222678017710160,
                          0.44169422454152618},
                         {126, 124, 1.9674072833921132, 1.8449257601074822},
                         {125.5, 125.5, 0.30697589692388100, 0.}}) {
                    h1.setMass(m1);
                    h2.setMass(m2);
                    heff.setMass(0.4 * m1 + 0.6 * m2);
                    CHECK(meas.chisq(pred, {}, options) ==
                          Approx(chisq).epsilon(5e-2));
                    CHECK(meas.chisq(pred, {}, options) -
                              meas.chisq(predEff, {}, options) ==
                          Approx(sepChisq).epsilon(5e-2));
                }
            }
        }
        SECTION("box THU test values") {

            auto options = HS::MeasurementOptions{};
            // avoids mass-dependence of the reference rates
            options.rescaleToRefMass = HS::RescaleToRefMass::never;
            options.theoryMassUncPDF = HS::PDF::box;
            SECTION("one particle-equivalent cases") {
                for (auto [m, chisq] : std::vector<std::pair<double, double>>{
                         {123., 18.},
                         {123.5, 6.0434027777778478},
                         {124., 0.14062500000001066},
                         {124.5, 0.},
                         {125., 0.},
                         {125.5, 0.},
                         {126., 0.},
                         {126.5, 2.9184027777777293},
                         {127., 18.}}) {
                    h1.setMass(m);
                    for (auto m2 : {121., m, 129.}) {
                        h2.setMass(m2);
                        INFO(m << " " << m2);
                        CHECK(meas.chisq(pred, {}, options) == Approx(chisq));
                    }
                }
            }

            h2.setMassUnc(0.5);

            SECTION("other test values") {
                auto predEff = HP::Predictions{};
                auto &heff = predEff.addParticle(HP::BsmParticle(
                    "heff", HP::ECharge::neutral, HP::CP::undefined));
                heff.setChannelRate(HP::Collider::LHC13, HP::Production::ggH,
                                    HP::Decay::gamgam, 1.);
                heff.setMassUnc(0.4 * h1.massUnc() + 0.6 * h2.massUnc());

                for (auto [m1, m2, chisq, sepChisq] :
                     std::vector<std::array<double, 4>>{
                         {124, 124.5, 0.1406, 0.},
                         {124.5, 125.5, 0., 0.},
                         {123.5, 125.5, 1.1906, 1.1906},
                         {126, 126, 0.7765, 0.}}) {
                    h1.setMass(m1);
                    h2.setMass(m2);
                    heff.setMass(0.4 * m1 + 0.6 * m2);
                    INFO(m1 << " " << m2);
                    CHECK(meas.chisq(pred, {}, options) ==
                          Approx(chisq).epsilon(5e-2));
                    CHECK(meas.chisq(pred, {}, options) -
                              meas.chisq(predEff, {}, options) ==
                          Approx(sepChisq).epsilon(5e-2));
                }
            }
        }
    }
}

TEST_CASE("coupling measurement") {
    const auto proc = HP::ChannelProcess{
        HP::Collider::LHC13, {{HP::Production::H, HP::Decay::tautau}}};

    const auto meas = HS::CouplingMeasurement(
        proc, {1.}, HP::Coupling::alphaCPTauYuk,
        HS::UncertainValue::fromInterval(
            -0.349065850398866, -0.0174532925199433, 0.314159265358979),
        std::make_shared<HP::SMHiggs>(125.), 2.5);

    CHECK(meas.observedCoupling() == Approx(-0.0174532925199433));
    CHECK(meas.observedCouplingUnc(HS::Uncertainty::plus) ==
          Approx(0.314159265358979 + 0.0174532925199433));
    CHECK(meas.observedCouplingUnc(HS::Uncertainty::minus) ==
          Approx(0.349065850398866 - 0.0174532925199433));

    auto pred = HP::Predictions{};
    auto &h1 = pred.addParticle(
        HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined));
    // h1 has no coupling set
    h1.setMass(125);
    h1.setChannelRate(HP::Collider::LHC13, HP::Production::H, HP::Decay::tautau,
                      1.);

    auto &h2 = pred.addParticle(
        HP::BsmParticle("h2", HP::ECharge::neutral, HP::CP::undefined));
    h2.setCoupling(HP::Coupling::alphaCPTauYuk, 0.5);
    h2.setMass(500); // h2 is out of mass range
    h2.setChannelRate(HP::Collider::LHC13, HP::Production::H, HP::Decay::tautau,
                      2.);
    auto &h3 = pred.addParticle(
        HP::BsmParticle("h3", HP::ECharge::neutral, HP::CP::undefined));
    h3.setCoupling(HP::Coupling::alphaCPTauYuk, -0.4);
    h3.setMass(125);
    // h3 has no rate in the coupling measurement process

    auto opts = HS::MeasurementOptions{};
    opts.rescaleToRefMass = HS::RescaleToRefMass::never;
    SECTION("unassigned") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        CHECK(residual == 0.);
        CHECK(obsVariance == Approx(0.349065850398866 - 0.0174532925199433));
        CHECK(refVariance == 0.);
        CHECK(extraChisq == opts.unassignedCouplingMeasurementPenalty);
    }

    h1.setCoupling(HP::Coupling::alphaCPTauYuk, 1.);
    SECTION("single particle assigned") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        CHECK(residual == -0.0174532925199433 - 1.);
        CHECK(obsVariance == Approx(0.349065850398866 - 0.0174532925199433));
        CHECK(refVariance == 0.);
        CHECK(extraChisq == 0.);
    }

    h2.setMass(126); // the precise mass value does not matter here
    SECTION("two particles assigned") {
        REQUIRE(meas.assignParticles(pred, HS::PDF::gaussian, 1.).size() == 2);
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        auto avgCoup = (1 + 2 * 0.5) / 3.;
        CHECK(residual == Approx(-0.0174532925199433 - avgCoup));
        CHECK(obsVariance == Approx(0.349065850398866 - 0.0174532925199433));
        CHECK(refVariance == 0.);
        CHECK(extraChisq ==
              Approx((1 * std::pow((1 - avgCoup) / obsVariance, 2) +
                      2 * std::pow((0.5 - avgCoup) / obsVariance, 2)) /
                     3.));
    }

    h3.setChannelRate(HP::Collider::LHC13, HP::Production::H, HP::Decay::tautau,
                      0.1);
    SECTION("all three particles assigned") {
        auto [residual, obsVariance, refVariance, extraChisq] =
            meas.evaluate(pred, {}, opts);
        auto avgCoup = (1 + 2 * 0.5 + 0.1 * (-0.4)) / 3.1;
        CHECK(residual == Approx(-0.0174532925199433 - avgCoup));
        CHECK(obsVariance == Approx(0.349065850398866 - 0.0174532925199433));
        CHECK(refVariance == 0.);
        CHECK(extraChisq ==
              Approx((1 * std::pow((1 - avgCoup) / obsVariance, 2) +
                      2 * std::pow((0.5 - avgCoup) / obsVariance, 2) +
                      0.1 * std::pow((-0.4 - avgCoup) / obsVariance, 2)) /
                     3.1));
    }
}

TEST_CASE("read SubMeasurement") {

    auto json = nlohmann::json::parse(R"([
            {
                "process": {
                    "channels": [["ggH","ZZ"]]
                },
                "obs": [0.9, 1.0, 1.1],
                "exp": [0.4, 0.5, 0.60],
                "channelWeights": [0.97],
                "massResolution": 2.5,
                "massSensitive": true
            },
            {
                "process": {
                    "channels": [["ggH","ZZ"]]
                },
                "obsMass": [122,125,130],
                "channelWeights":[1.1]
            },
            {
                "process": {
                    "channels": [["Htt","ZZ"]]
                },
                "coupling": "effCPoTopYuk",
                "obsCoupling": [-0.3, 0.01, 0.4],
                "channelWeights": [2.0],
                "massResolution": 2.4
            },
            {
                "process": {
                    "channels": [["ggH","ZZ"]]
                },
                "channelWeights":[1.1]
            }
        ])");
    REQUIRE_NOTHROW(HS::readSubMeasurement(
        json[0], HP::Collider::LHC8, std::make_shared<HP::SMHiggs>(120), 4));
    REQUIRE_NOTHROW(HS::readSubMeasurement(
        json[1], HP::Collider::LHC8, std::make_shared<HP::SMHiggs>(120), 4));
    REQUIRE_NOTHROW(HS::readSubMeasurement(
        json[2], HP::Collider::LHC13, std::make_shared<HP::SMHiggs>(125), 4));
    REQUIRE_THROWS_AS(HS::readSubMeasurement(json[3], HP::Collider::LHC8,
                                             std::make_shared<HP::SMHiggs>(120),
                                             4),
                      HS::InvalidMeasurement);
}
