#include "Higgs/Signals.hpp"
#include "Higgs/signals/Measurement.hpp"
#include "signals/Evaluation.hpp"
#include "signals/MeasurementImpl.hpp"
#include "utilities/Algorithm.hpp"
#include "utilities/Json.hpp"
#include "utilities/Logging.hpp"
#include <Eigen/Core>
#include <exception>
#include <filesystem>
#include <functional>
#include <memory>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <string_view>
#include <utility>

namespace fs = std::filesystem;

namespace Higgs::signals {

namespace {

std::vector<Measurement>
loadMeasurements(std::string_view dataPath,
                 const MeasurementOptions &measurementOptions) {
    auto log = logger();
    auto measurements = std::vector<Measurement>();
    for (const auto &p : fs::recursive_directory_iterator(dataPath)) {
        if (fs::is_regular_file(p) && p.path().extension() == ".json") {
            const auto filePath = p.path().lexically_normal().string();
            log->trace("reading measurement from {}", filePath);

            try {
                auto measurement = Measurement(filePath, measurementOptions);
                const auto duplicateId = [&measurement](const auto &m) {
                    return m.id() == measurement.id();
                };
                auto found = ranges::find_if(measurements, duplicateId);
                if (found == measurements.end()) {
                    measurements.emplace_back(std::move(measurement));
                } else {
                    log->warn("Duplicate measurement id {} for files 1: {} and "
                              "2: {}, "
                              "skipping 2.",
                              measurement.id(), found->loadedFrom(), filePath);
                }
            } catch (const nlohmann::json::parse_error &) {
                log->warn("Skipping invalid json file {}", filePath);
            } catch (const std::exception &e) {
                log->warn("Skipping {} after read error: {}", filePath,
                          e.what());
            }
        }
    }
    return measurements;
}

} // namespace

struct GlobalCorrelations {
    explicit GlobalCorrelations(std::size_t size)
        : experimental{Eigen::MatrixXd::Identity(size, size)},
          theory{Eigen::MatrixXd::Identity(size, size)} {}

    Eigen::MatrixXd experimental;
    Eigen::MatrixXd theory;
};

Signals::Signals(const std::string &dataPath,
                 const MeasurementOptions &measurementOptions)
    : measurements_{loadMeasurements(dataPath, measurementOptions)},
      nSubMeasurements_{ranges::accumulate(measurements_, std::size_t{0},
                                           std::plus{},
                                           &Measurement::nSubMeasurements)},
      corr_{std::make_shared<GlobalCorrelations>(nSubMeasurements_)} {
    int subMeasCounter = 0;
    for (const auto &meas : measurements_) {
        const auto nSubMeas = meas.nSubMeasurements();
        corr_->experimental.block(subMeasCounter, subMeasCounter, nSubMeas,
                                  nSubMeas) = meas.data_->corrMatExp_;
        corr_->theory.block(subMeasCounter, subMeasCounter, nSubMeas,
                            nSubMeas) = meas.data_->corrMatTheo_;
        subMeasCounter += nSubMeas;
    }
}

double Signals::operator()(
    const Predictions &predictions,
    const std::unordered_map<
        std::size_t, std::unordered_map<std::string, std::vector<double>>>
        &modificationFactors) const {
    auto residuals = Eigen::VectorXd(nSubMeasurements_);
    auto observedVariances = Eigen::VectorXd(nSubMeasurements_);
    auto expectedVariances = Eigen::VectorXd(nSubMeasurements_);

    double extraChisq = 0.;
    int subMeasCounter = 0;
    for (const auto &measurement : measurements_) {
        const auto nSubMeas = measurement.nSubMeasurements();
        auto resBlock = residuals.segment(subMeasCounter, nSubMeas);
        auto modelRateVarBlock =
            observedVariances.segment(subMeasCounter, nSubMeas);
        auto refRateVarBlock =
            expectedVariances.segment(subMeasCounter, nSubMeas);

        extraChisq += evaluateSubMeasurements(
            predictions, measurement,
            utilities::getWithDefault(modificationFactors, measurement.id(),
                                      Measurement::ModificationFactors{}),
            resBlock, modelRateVarBlock, refRateVarBlock);
        subMeasCounter += nSubMeas;
    }

    auto covariance = computeCovariance(corr_->experimental, observedVariances,
                                        corr_->theory, expectedVariances);

    return computeChisq(residuals, covariance) + extraChisq;
}

std::size_t Signals::observableCount() const noexcept {
    return nSubMeasurements_;
}

const std::vector<Measurement> &Signals::measurements() const noexcept {
    return measurements_;
}
} // namespace Higgs::signals
