#include "signals/Evaluation.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "prettyprint.hpp"
#include "signals/Measurements/RateMeasurement.hpp"
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <magic_enum.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/map.hpp>

namespace HP = Higgs::predictions;
namespace HS = Higgs::signals;
using Catch::Approx;

TEST_CASE("chisq computation") {

    Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(4, 4);
    Eigen::VectorXd residuals(4);

    CHECK(HS::computeChisq(Eigen::VectorXd::Zero(4), covariance) == 0);

    residuals << -0.1, 0.02, 0.14, -0.3;
    CHECK(HS::computeChisq(residuals, covariance) == residuals.dot(residuals));

    covariance(1, 2) = 0.3;
    covariance(2, 1) = 0.3;
    covariance(1, 3) = -0.1;
    covariance(3, 1) = -0.1;

    CHECK(HS::computeChisq(residuals, covariance) ==
          Approx(residuals.transpose() * covariance.inverse() * residuals));
}

TEST_CASE("compare measurement to old HS") {
    INFO("the measurement uses the SMHiggs without EW corrections as reference "
         "for comparability with HS-2");
    auto opts = HS::MeasurementOptions{};
    opts.rescaleToRefMass =
        HS::RescaleToRefMass::always; // that's the HS-2 behaviour
    auto meas = HS::Measurement(TESTMEASUREMENTS_PATH "STXS_test.json", opts);
    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
    h.setMass(125.09);
    auto hs2mu = 0.9998777;
    HP::effectiveCouplingInput(h,
                               HP::scaledSMlikeEffCouplings(std::sqrt(hs2mu)),
                               HP::ReferenceModel::SMHiggs);

    for (const auto &rm : meas.subMeasurements()) {
        INFO(rm.first);
        CHECK(rm.second->signalStrength(h, {}, HS::NormalizeAt::particleMass) ==
              Approx(hs2mu));
    }

    Eigen::VectorXd residuals = Eigen::VectorXd(meas.nSubMeasurements());
    Eigen::VectorXd modelRateVariances =
        Eigen::VectorXd(meas.nSubMeasurements());
    Eigen::VectorXd referenceRateVariances =
        Eigen::VectorXd(meas.nSubMeasurements());

    auto extraChisq = evaluateSubMeasurements(
        pred, meas, {}, residuals, modelRateVariances, referenceRateVariances);
    REQUIRE(extraChisq == 0.);

    const auto hs2Residuals = std::unordered_map<std::string, double>{
        {"ttH", .025 - .0154},
        {"gg2H_pTH_High", .038 - .015},
        {"qq2Hqq_VBF", .15 - .10759},
        {"qq2Hqq_BSM", .0005 - .00420},
        {"gg2H_1j_pTH_Med", .17 - .11899},
        {"gg2H_1j_pTH_Low", .05 - .17198},
        {"gg2H_2j", .04 - .12698},
        {"gg2H_1j_pTH_High", .009 - .02},
        {"gg2H_0j_pTH_Low", .17 - .17598},
        {"qq2Hqq_VH", .021 - .0138},
        {"gg2H_0j_pTH_High", 0.63 - .54993},
        {"VH_Lep", .022 - .0164}};
    const auto hs2modelVariances =
        std::unordered_map<std::string, double>{{"ttH", 1.7e-2},
                                                {"gg2H_pTH_High", 1.6e-2},
                                                {"qq2Hqq_VBF", 5.2e-2},
                                                {"qq2Hqq_BSM", 7.9e-3},
                                                {"gg2H_1j_pTH_Med", 5e-2},
                                                {"gg2H_1j_pTH_Low", 8e-2},
                                                {"gg2H_2j", 7.5e-2},
                                                {"gg2H_1j_pTH_High", 1.6e-2},
                                                {"gg2H_0j_pTH_Low", 5.5e-2},
                                                {"qq2Hqq_VH", .035},
                                                {"gg2H_0j_pTH_High", 0.11},
                                                {"VH_Lep", 1.8e-2}};

    const auto hs2refVariances =
        std::unordered_map<std::string, double>{{"ttH", 1.3e-3},
                                                {"gg2H_pTH_High", 4e-3},
                                                {"qq2Hqq_VBF", 3.5e-3},
                                                {"qq2Hqq_BSM", 1.8e-4},
                                                {"gg2H_1j_pTH_Med", 1.8e-2},
                                                {"gg2H_1j_pTH_Low", 2.5e-2},
                                                {"gg2H_2j", 2.7e-2},
                                                {"gg2H_1j_pTH_High", 4e-3},
                                                {"gg2H_0j_pTH_Low", 2.5e-2},
                                                {"qq2Hqq_VH", 6e-4},
                                                {"gg2H_0j_pTH_High", 4e-2},
                                                {"VH_Lep", 4e-4}};
    for (auto [i, key] : ranges::views::enumerate(meas.subMeasurements() |
                                                  ranges::views::keys)) {
        INFO(key);
        CHECK(residuals(i) == Approx(hs2Residuals.at(key)).margin(1e-5));
        CHECK(modelRateVariances(i) == Approx(hs2modelVariances.at(key)));
        CHECK(referenceRateVariances(i) ==
              Approx(hs2mu * hs2refVariances.at(key)));
    }
    // agreement not perfect since HS-2 uses rounded experimental correlations
    REQUIRE(meas(pred, {}) == Approx(7.22739420).margin(0.05));
}

TEST_CASE("compare run1 combination measurement to old HS") {
    auto meas =
        HS::Measurement(TESTMEASUREMENTS_PATH "/run1comb_as_in_HS2.json");
    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined));
    h.setMass(125.09);
    h.setCxn(HP::Collider::LHC8, HP::Production::ggH, 18.484981018066406);
    h.setCxn(HP::Collider::LHC8, HP::Production::bbH, 0.20163561889648438);
    h.setCxn(HP::Collider::LHC8, HP::Production::vbfH, 1.6486860534667969);
    h.setCxn(HP::Collider::LHC8, HP::Production::HW, 0.75244766723632817);
    h.setCxn(HP::Collider::LHC8, HP::Production::qqHZ, 0.41478559258036291);
    h.setCxn(HP::Collider::LHC8, HP::Production::Htt, 0.12097043443583963);
    h.setCxn(HP::Collider::LHC8, HP::Production::schanHt,
             1.2113001098632812E-003);
    h.setCxn(HP::Collider::LHC8, HP::Production::tchanHt,
             1.8663001098632809E-002);
    h.setCxn(HP::Collider::LHC8, HP::Production::HtW, 3.5100000000000001E-003);

    h.setTotalWidth(0.1); // doesn't matter
    h.setBr(HP::Decay::gamgam, 2.3091776788432333E-003);
    h.setBr(HP::Decay::WW, 0.20954031578718774);
    h.setBr(HP::Decay::ZZ, 2.6436157910051848E-002);
    h.setBr(HP::Decay::tautau, 6.3339261734749477E-002);
    h.setBr(HP::Decay::bb, 0.58943399070344793);

    for (const auto &rm : meas.subMeasurements()) {
        INFO(rm.first);
        // the ttH numbers in the old HB include EW corrections, so they match
        // much worse here
        if (rm.first.substr(0, 3) != "Htt") {
            CHECK(rm.second->signalStrength(h, {},
                                            HS::NormalizeAt::particleMass) ==
                  Approx(1).margin(1e-2));
        } else {
            CHECK(rm.second->signalStrength(h, {},
                                            HS::NormalizeAt::particleMass) ==
                  Approx(1).margin(0.2));
        }
    }

    { // the whole normalization procedure should be a no-op in this case
        auto rm = dynamic_cast<HS::RateMeasurement *>(
            meas.subMeasurements().at("ggH_gamgam").get());
        CHECK(rm->signalStrength(h, {}, HS::NormalizeAt::particleMass) *
                  rm->referenceRate() ==
              Approx(h.cxn(HP::Collider::LHC8, HP::Production::H) *
                     h.br(HP::Decay::gamgam)));
    }

    { // the whole normalization procedure should be a no-op in this case
        auto rm = dynamic_cast<HS::RateMeasurement *>(
            meas.subMeasurements().at("Htt_WW").get());
        CHECK(rm->signalStrength(h, {}, HS::NormalizeAt::particleMass) *
                  rm->referenceRate() ==
              Approx((h.cxn(HP::Collider::LHC8, HP::Production::Htt) +
                      h.cxn(HP::Collider::LHC8, HP::Production::Ht) +
                      h.cxn(HP::Collider::LHC8, HP::Production::HtW)) *
                     h.br(HP::Decay::WW)));
    }

    Eigen::VectorXd residuals = Eigen::VectorXd(meas.nSubMeasurements());
    Eigen::VectorXd modelRateVariances =
        Eigen::VectorXd(meas.nSubMeasurements());
    Eigen::VectorXd referenceRateVariances =
        Eigen::VectorXd(meas.nSubMeasurements());

    auto extraChisq = evaluateSubMeasurements(
        pred, meas, {}, residuals, modelRateVariances, referenceRateVariances);
    REQUIRE(extraChisq == 0.);

    auto hs2Residuals = std::unordered_map<std::string, double>{
        {"vbfH_tautau", 0.125 - 0.10442655745896462},
        {"ggH_ZZ", 0.58 - 0.49400234821935302},
        {"HW_WW", 0.24 - 0.15766812180603296},
        {"HW_bb", 0.42 - 0.44351823129460893},
        {"HZ_bb", 8.0E-002 - 0.24448872712093778},
        {"ggH_gamgam", 4.8E-002 - 4.3150718031175317E-002},
        {"HZ_gamgam", 5.0E-004 - 9.5781363189233746E-004},
        {"Htt_bb", 8.0E-002 - 8.5087587907782078E-002},
        {"vbfH_ZZ", 3.0E-003 - 4.3584924853548423E-002},
        {"Htt_tautau", -1.5E-002 - 9.1433223836271515E-003},
        {"vbfH_WW", 0.39 - 0.34546619627736491},
        {"Htt_gamgam", 6.4E-004 - 3.3334073338521573E-004},
        {"vbfH_gamgam", 4.6E-003 - 3.8071090340856687E-003},
        {"ggH_WW", 3.5 - 3.9155995511033206},
        {"ggH_tautau", 1.3 - 1.1835965021055166},
        {"HW_tautau", -6.4E-002 - 4.7659479736783465E-002},
        {"HW_gamgam", 7.0E-004 - 1.7375353576797898E-003},
        {"HZ_WW", 0.53 - 8.6914304053265043E-002},
        {"HZ_tautau", 5.8E-002 - 2.6272213212250768E-002},
        {"Htt_WW", 0.14 - 3.0248136892290115E-002},
        {"mh125", 0.}};

    auto hs2modelVariances =
        std::unordered_map<std::string, double>{{"vbfH_tautau", 3.7E-002},
                                                {"ggH_ZZ", 0.16},
                                                {"HW_WW", 0.16},
                                                {"HW_bb", 0.21},
                                                {"HZ_bb", 9E-002},
                                                {"ggH_gamgam", 9.7E-003},
                                                {"HZ_gamgam", 2.9E-003},
                                                {"Htt_bb", 7E-002},
                                                {"vbfH_ZZ", 4.6E-002},
                                                {"Htt_tautau", 3E-002},
                                                {"vbfH_WW", 0.13},
                                                {"Htt_gamgam", 3.8E-004},
                                                {"vbfH_gamgam", 1.8E-003},
                                                {"ggH_WW", 0.7},
                                                {"ggH_tautau", 0.7},
                                                {"HW_tautau", 6.4E-002},
                                                {"HW_gamgam", 2.1E-003},
                                                {"HZ_WW", 0.2},
                                                {"HZ_tautau", 4.7E-002},
                                                {"Htt_WW", 5E-002},
                                                {"mh125", 0.24}};

    for (auto [i, key] : ranges::views::enumerate(meas.subMeasurements() |
                                                  ranges::views::keys)) {
        INFO(key);
        CHECK(residuals(i) == Approx(hs2Residuals.at(key)));
        CHECK(modelRateVariances(i) == Approx(hs2modelVariances.at(key)));
    }
    // permille disagreement from chisq computation
    REQUIRE(meas(pred, {}) == Approx(20.593038620127679).epsilon(5e-3));
}
