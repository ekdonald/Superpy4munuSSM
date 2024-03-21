#include "Higgs/signals/Measurement.hpp"
#include "Higgs/Predictions.hpp"
#include "signals/Evaluation.hpp"
#include "signals/MeasurementImpl.hpp"
#include "utilities/Algorithm.hpp"
#include "utilities/Json.hpp"
#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <utility>
#include "utilities/Logging.hpp"

namespace Higgs::signals {

double SubMeasurement::chisq(const predictions::Predictions &predictions,
                             const ModificationFactors &modFacs,
                             const MeasurementOptions &options) const {
    static const auto log = logger();
    auto [residual, obsVariance, theoVariance, extraChisq] =
        evaluate(predictions, modFacs, options);
    // log->info("residual {}", residual);
    // log->info("obsVariance {}", obsVariance);
    // log->info("theoVariance {}", theoVariance);
    // log->info("extraChisq {}", extraChisq);
    if (options.ignoreTheoryUncertainties) {
        return std::pow(residual / obsVariance, 2) + extraChisq;
    } else {
        return std::pow(residual, 2) /
                   (std::pow(obsVariance, 2) + std::pow(theoVariance, 2)) +
               extraChisq;
    }
}

Measurement::Measurement(const std::string &filePath,
                         const MeasurementOptions &options)
    : data_{std::make_shared<MeasurementData>(
          nlohmann::json::parse(std::ifstream{filePath}), filePath, options)} {}

double Measurement::operator()(
    const Predictions &predictions,
    const Measurement::ModificationFactors &modificationFactors) const {
    auto residuals = Eigen::VectorXd(nSubMeasurements());
    auto observedVariances = Eigen::VectorXd(nSubMeasurements());
    auto expectedVariances = Eigen::VectorXd(nSubMeasurements());

    double extraChisq = evaluateSubMeasurements(
        predictions, *this, modificationFactors, residuals, observedVariances,
        expectedVariances);

    auto covariance = computeCovariance(data_->corrMatExp_, observedVariances,
                                        data_->corrMatTheo_, expectedVariances);

    return computeChisq(residuals, covariance) + extraChisq;
}

std::unordered_map<std::string, double> Measurement::chisqContributions(
    const Predictions &predictions,
    const Measurement::ModificationFactors &modificationFactors) const {
    // writing this as a range conversions yields errors on Mac
    auto result = std::unordered_map<std::string, double>{};
    for (const auto &[key, subMeasurement] : subMeasurements()) {
        result.try_emplace(
            key, subMeasurement->chisq(
                     predictions,
                     utilities::getWithDefault(modificationFactors, key, {}),
                     options()));
    }
    return result;
}

std::size_t Measurement::nSubMeasurements() const noexcept {
    return data_->subMeasurements_.size();
}

const MeasurementOptions &Measurement::options() const noexcept {
    return data_->options_;
}
std::size_t Measurement::id() const noexcept { return data_->id_; }
const std::string &Measurement::reference() const noexcept {
    return data_->reference_;
}
const std::string &Measurement::citeKey() const noexcept {
    return data_->citeKey_;
}
predictions::Collider Measurement::collider() const noexcept {
    return data_->collider_;
}
predictions::Experiment Measurement::experiment() const noexcept {
    return data_->experiment_;
}
double Measurement::luminosity() const noexcept { return data_->luminosity_; }
const std::string &Measurement::loadedFrom() const noexcept {
    return data_->loadedFrom_;
}
double Measurement::referenceMass() const noexcept { return data_->refMass_; }
predictions::ReferenceModel Measurement::referenceModel() const noexcept {
    return data_->referenceModel_;
}
const std::unordered_map<std::string, std::shared_ptr<SubMeasurement>> &
Measurement::subMeasurements() const noexcept {
    return data_->subMeasurements_;
}
} // namespace Higgs::signals
