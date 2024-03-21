#include "signals/Evaluation.hpp"
#include "Higgs/signals/Measurement.hpp"
#include "utilities/Algorithm.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <memory>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/unreachable_sentinel.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/zip_with.hpp>
#include <unordered_map>
#include <utility>

namespace Higgs::signals {

double evaluateSubMeasurements(
    const Predictions &predictions, const Measurement &measurement,
    const std::unordered_map<std::string, std::vector<double>>
        &modificationFactors,
    Eigen::Ref<Eigen::VectorXd> residuals,
    Eigen::Ref<Eigen::VectorXd> modelVariances,
    Eigen::Ref<Eigen::VectorXd> referenceVariances) {
    auto extraChisq = 0.;
    for (const auto &[i, m] :
         ranges::views::enumerate(measurement.subMeasurements())) {

        const auto rateModFactors =
            utilities::getWithDefault(modificationFactors, m.first, {});

        const auto evaluated = m.second->evaluate(predictions, rateModFactors,
                                                  measurement.options());
        residuals[i] = evaluated.residual;
        modelVariances[i] = evaluated.obsVariance;
        referenceVariances[i] = evaluated.refVariance;
        extraChisq += evaluated.extraChisq;
    }
    return extraChisq;
}

Eigen::MatrixXd
computeCovariance(const Eigen::MatrixXd &experimentalCorrMat,
                  const Eigen::VectorXd &observedVariances,
                  const Eigen::MatrixXd &theoreticalCorrMat,
                  const Eigen::VectorXd &expectedVariances) noexcept {

    return experimentalCorrMat.cwiseProduct(observedVariances *
                                            observedVariances.transpose()) +
           theoreticalCorrMat.cwiseProduct(expectedVariances *
                                           expectedVariances.transpose());
}

Eigen::VectorXd
computeChisqContributions(const Eigen::VectorXd &residuals,
                          const Eigen::MatrixXd &covarianceMatrix) noexcept {
    return covarianceMatrix.llt().solve(residuals).cwiseProduct(residuals);
}

double computeChisq(const Eigen::VectorXd &residuals,
                    const Eigen::MatrixXd &covarianceMatrix) noexcept {
    return computeChisqContributions(residuals, covarianceMatrix).sum();
}

} // namespace Higgs::signals
