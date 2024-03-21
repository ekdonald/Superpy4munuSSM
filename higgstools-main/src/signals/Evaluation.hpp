#pragma once

#include "Higgs/Predictions.hpp"
#include <Eigen/Core>
#include <string>
#include <unordered_map>
#include <vector>

namespace Higgs::signals {
class Measurement;


double evaluateSubMeasurements(
    const Predictions &predictions, const Measurement &measurement,
    const std::unordered_map<std::string, std::vector<double>>
        &modificationFactors,
    Eigen::Ref<Eigen::VectorXd> residuals,
    Eigen::Ref<Eigen::VectorXd> modelVariances,
    Eigen::Ref<Eigen::VectorXd> referenceVariances);

Eigen::MatrixXd
computeCovariance(const Eigen::MatrixXd &experimentalCorrMat,
                  const Eigen::VectorXd &observedVariances,
                  const Eigen::MatrixXd &theoreticalCorrMat,
                  const Eigen::VectorXd &expectedVariances) noexcept;

Eigen::VectorXd
computeChisqContributions(const Eigen::VectorXd &residuals,
                          const Eigen::MatrixXd &covarianceMatrix) noexcept;

double computeChisq(const Eigen::VectorXd &residuals,
                    const Eigen::MatrixXd &covarianceMatrix) noexcept;

} // namespace Higgs::signals
