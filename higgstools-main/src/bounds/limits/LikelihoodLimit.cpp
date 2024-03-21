#include "bounds/limits/LikelihoodLimit.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "bounds/Likelihood.hpp"
#include "predictions/Clustering.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "predictions/UncertainMass.hpp"
#include "utilities/Json.hpp"
#include "utilities/Logging.hpp"
#include <algorithm>
#include <exception>
#include <iterator>
#include <map>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>
namespace Higgs::bounds {

namespace {

template <std::size_t llhDim>
std::array<double, llhDim + 1> scaleRates(std::array<double, llhDim + 1> point,
                                          double scale) {
    ranges::for_each(point | ranges::views::drop(1),
                     [scale](auto &val) { val *= scale; });
    return point;
}

template <std::size_t llhDim, class Llh>
void validateLlh(const std::vector<double> &masses, const Llh &exp,
                 const Llh &obs) {
    for (const auto m : masses) {
        const auto x0 = scaleRates<llhDim>({m}, 0.);
        if (Likelihood::expCLs(exp(x0)) < Likelihood::limCLs ||
            Likelihood::obsCLs(exp(x0), obs(x0)) < Likelihood::limCLs) {
            throw std::logic_error{"Invalid likelihood: excludes zero rates"};
        }
    }
}

constexpr auto readVector = utilities::readAs<std::vector<double>>;
constexpr auto readMassRes =
    utilities::readIfPresent<predictions::MassResolution>;

std::vector<double> readStackedMasses(const nlohmann::json &j) {
    return j.at("stackedLlhGrid") |
           ranges::views::transform([](const auto &slice) {
               return utilities::readAs<double>(slice, "mass");
           }) |
           ranges::to<std::vector>;
}

utilities::LinearInterpolator<1, utilities::LinearInterpolator<2>>
read2dStackedLlh(const nlohmann::json &j, const std::string &which) {
    auto readMassSlice = [&which](const nlohmann::json &slice) {
        return utilities::LinearInterpolator{
            std::array{readVector(slice, "/channels/0"),
                       readVector(slice, "/channels/1")},
            readVector(slice, which)};
    };
    auto innerInterps = j.at("stackedLlhGrid") |
                        ranges::views::transform(readMassSlice) |
                        ranges::to<std::vector>;
    return utilities::LinearInterpolator{std::array{readStackedMasses(j)},
                                         innerInterps};
}
} // namespace

LikelihoodLimit1d::LikelihoodLimit1d(const nlohmann::json &data,
                                     const std::string &loadedFrom,
                                     const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolutions_{readMassRes(data, "/analysis/massResolution")},
      processes_{predictions::readChannelProcess(data.at("process").at(0),
                                                 collider())},
      obs_{std::array{readVector(data, "/analysis/grid/mass"),
                      readVector(data, "/analysis/grid/channels/0")},
           readVector(data, "/analysis/likelihood/observed")},
      exp_{std::array{readVector(data, "/analysis/grid/mass"),
                      readVector(data, "/analysis/grid/channels/0")},
           readVector(data, "/analysis/likelihood/expected")},
      constraints_{readConstraints(data.value("constraints", nlohmann::json{}),
                                   processes_[0])} {
    validateLlh<llhDim>(readVector(data, "/analysis/grid/mass"), exp_, obs_);
}

LikelihoodLimit2d::LikelihoodLimit2d(const nlohmann::json &data,
                                     const std::string &loadedFrom,
                                     const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolutions_{readMassRes(data, "/analysis/massResolution")},
      processes_{
          predictions::readChannelProcess(data.at("process").at(0), collider()),
          predictions::readChannelProcess(data.at("process").at(1),
                                          collider())},
      mergedProcess_{processes_[0] + processes_[1]}, obs_{read2dStackedLlh(
                                                         data.at("analysis"),
                                                         "observed")},
      exp_{read2dStackedLlh(data.at("analysis"), "expected")},
      constraints_{readConstraints(data.value("constraints", nlohmann::json{}),
                                   mergedProcess_)} {
    validateLlh<llhDim>(readStackedMasses(data.at("analysis")), exp_, obs_);
}

std::string LikelihoodLimit1d::processDesc() const noexcept {
    return fmt::format("1d likelihood {}", processes_.front().to_string());
}
std::string LikelihoodLimit1d::extentDesc() const noexcept {
    return fmt::format("M={}", obs_.extent().front());
}

std::string LikelihoodLimit2d::processDesc() const noexcept {
    return fmt::format("2d likelihood {{{}, {}}}", processes_[0].to_string(),
                       processes_[1].to_string());
}
std::string LikelihoodLimit2d::extentDesc() const noexcept {
    return fmt::format("M={}", obs_.extent().front());
}

namespace {
template <std::size_t gridDim, class ExpLlhFunc>
std::array<double, gridDim>
likelihoodGridPoint(const std::array<double, gridDim> &clusterRates,
                    const predictions::UncertainMass &clusterMass,
                    const ExpLlhFunc &expectedLikelihood,
                    predictions::MassUncEagerness setLimitEagerness) {
    const auto constructGridAt = [clusterMass,
                                  &clusterRates](const auto &massFunc) {
        auto res = clusterRates;
        res[0] = massFunc(clusterMass);
        return res;
    };

    const auto lessExpLlh = [&expectedLikelihood](const auto &x1,
                                                  const auto &x2) {
        return expectedLikelihood(x1) < expectedLikelihood(x2);
    };

    switch (setLimitEagerness) {
    case predictions::MassUncEagerness::ignore:
        return constructGridAt(predictions::centralMass);
    default:
        logger()->error(
            "Unknown MassUncEagerness in "
            "LikelihoodLimit::applyLLhLimit(), defaulting to cautious.");
        [[fallthrough]];
    case predictions::MassUncEagerness::cautious:
        return std::min({constructGridAt(predictions::centralMass),
                         constructGridAt(predictions::lowerMassBound),
                         constructGridAt(predictions::upperMassBound)},
                        lessExpLlh);
    case predictions::MassUncEagerness::eager:
        return std::max({constructGridAt(predictions::centralMass),
                         constructGridAt(predictions::lowerMassBound),
                         constructGridAt(predictions::upperMassBound)},
                        lessExpLlh);
    }
}
} // namespace

std::vector<AppliedLimit>
LikelihoodLimit2d::apply(const predictions::Predictions &prediction) const {
    logger()->trace("Applying likelihood limit {}: {} {}", id(), processDesc(),
                    extentDesc());

    auto relevantParticles = predictions::Clustering::relevantParticles(
        mergedProcess_, prediction,
        applicableMassRanges(obs_.extent(), massResolutions_),
        options().applicableMassUnc);

    auto clustersForRoles = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolutions_, options().clusterMassUnc,
        constraints_, prediction);

    const auto applyToCluster = [this](auto &&c) {
        const auto weightedMass = std::get<0>(
            predictions::Clustering::rateWeightedMasses(mergedProcess_, c));

        auto rates = std::array<double, 1 + llhDim>{};
        ranges::transform(
            processes_, std::next(rates.begin()), [&c](const auto &p) {
                return predictions::Clustering::combinedRate(p, c);
            });

        const auto x = likelihoodGridPoint(rates, weightedMass, exp_,
                                           options().setLimitMassUnc);

        logger()->trace("Cluster of mass and rates {}", x);
        const auto qExp = [&x, this](double s) {
            return exp_(scaleRates<llhDim>(x, s));
        };
        const auto qObs = [&x, this](double s) {
            return obs_(scaleRates<llhDim>(x, s));
        };
        auto [expRatio, obsRatio] = Likelihood::limitsFromLlh(qExp, qObs);
        return AppliedLimit{
            shared_from_this(),
            obsRatio,
            expRatio,
            predictions::ChannelProcess::contributingParticles(c),
            qObs(1),
            qExp(1)};
    };

    return std::get<0>(clustersForRoles) |
           ranges::views::transform(applyToCluster) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}

std::vector<AppliedLimit>
LikelihoodLimit1d::apply(const predictions::Predictions &prediction) const {
    logger()->trace("Applying likelihood limit {}: {} {}", id(), processDesc(),
                    extentDesc());

    auto relevantParticles = predictions::Clustering::relevantParticles(
        std::get<0>(processes_), prediction,
        applicableMassRanges(obs_.extent(), massResolutions_),
        options().applicableMassUnc);

    auto clustersForRoles = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolutions_, options().clusterMassUnc,
        constraints_, prediction);

    const auto applyToCluster = [this](auto &&c) {
        const auto weightedMass =
            std::get<0>(predictions::Clustering::rateWeightedMasses(
                std::get<0>(processes_), c));

        auto rates = std::array<double, 1 + llhDim>{};
        ranges::transform(
            processes_, std::next(rates.begin()), [&c](const auto &p) {
                return predictions::Clustering::combinedRate(p, c);
            });

        const auto x = likelihoodGridPoint(rates, weightedMass, exp_,
                                           options().setLimitMassUnc);

        logger()->trace("Cluster of mass and rates {}", x);
        const auto qExp = [&x, this](double s) {
            return exp_(scaleRates<llhDim>(x, s));
        };
        const auto qObs = [&x, this](double s) {
            return obs_(scaleRates<llhDim>(x, s));
        };
        auto [expRatio, obsRatio] = Likelihood::limitsFromLlh(qExp, qObs);
        return AppliedLimit{
            shared_from_this(),
            obsRatio,
            expRatio,
            predictions::ChannelProcess::contributingParticles(c),
            qObs(1),
            qExp(1)};
    };

    return std::get<0>(clustersForRoles) |
           ranges::views::transform(applyToCluster) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}
} // namespace Higgs::bounds
