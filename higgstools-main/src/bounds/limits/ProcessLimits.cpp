#include "bounds/limits/ProcessLimits.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Process.hpp"
#include "predictions/Clustering.hpp"
#include "predictions/Constraints.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "utilities/Algorithm.hpp"
#include "utilities/ArithmeticArray.hpp"
#include "utilities/Json.hpp"
#include "utilities/Logging.hpp"
#include <array>
#include <exception>
#include <map>
#include <optional>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>
#include <tuple>

namespace Higgs::predictions {
class Particle;
}
namespace Higgs::bounds {
namespace {
constexpr auto readVector = utilities::readAs<std::vector<double>>;
constexpr auto readMassRes =
    utilities::readIfPresent<predictions::MassResolution>;
} // namespace

ChannelLimit::ChannelLimit(const nlohmann::json &data,
                           const std::string &loadedFrom,
                           const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolution_{readMassRes(data, "/analysis/massResolution")},
      process_{predictions::readChannelProcess(data.at("process"), collider())},
      limit_{{readVector(data, "/analysis/grid/mass")},
             zipLimits(readVector(data, "/analysis/limit/expected"),
                       readVector(data, "/analysis/limit/observed"))},
      normalization_{predictions::ReferenceRate::read(
          data.value("normalization", nlohmann::json{}), process_)},
      constraints_{readConstraints(data.value("constraints", nlohmann::json{}),
                                   process_)},
      acceptances_{data.contains("/analysis/acceptances"_json_pointer)
                       ? Acceptances(data["/analysis/acceptances"_json_pointer],
                                     process_.size())
                       : std::optional<Acceptances>{std::nullopt}} {
    if (acceptances_ && normalization_) {
        throw std::runtime_error(
            "Limits with both normalization and non-trivial acceptances are "
            "not yet implemented!");
    }
}

ChainDecayLimit::ChainDecayLimit(const nlohmann::json &data,
                                 const std::string &loadedFrom,
                                 const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolutions_{readMassRes(data, "/analysis/massResolution/mother"),
                       readMassRes(data, "/analysis/massResolution/daughter")},
      process_{
          predictions::readChainDecayProcess(data.at("process"), collider())},
      limit_{{readVector(data, "/analysis/grid/massMother"),
              readVector(data, "/analysis/grid/massDaughter")},
             zipLimits(readVector(data, "/analysis/limit/expected"),
                       readVector(data, "/analysis/limit/observed"))},
      constraints_{
          readConstraints(
              data.value("/constraints/mother"_json_pointer, nlohmann::json{}),
              collider()),
          readConstraints(data.value("/constraints/daughter"_json_pointer,
                                     nlohmann::json{}),
                          collider())},
      productionAcceptances_{
          data.contains("/analysis/productionAcceptances"_json_pointer)
              ? Acceptances(
                    data["/analysis/productionAcceptances"_json_pointer],
                    process_.productionSize())
              : std::optional<Acceptances>{std::nullopt}} {}

PairDecayLimit::PairDecayLimit(const nlohmann::json &data,
                               const std::string &loadedFrom,
                               const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolutions_{
          readMassRes(data, "/analysis/massResolution/mother"),
          readMassRes(data, "/analysis/massResolution/firstDaughter"),
          readMassRes(data, "/analysis/massResolution/secondDaughter")},
      process_{
          predictions::readPairDecayProcess(data.at("process"), collider())},
      equalDaughterMasses_{utilities::readIfPresent<bool>(
          data, "/analysis/equalDaughterMasses", false)},
      limit_{{readVector(data, "/analysis/grid/massMother"),
              readVector(data, "/analysis/grid/massFirstDaughter"),
              equalDaughterMasses_
                  ? std::vector{0.}
                  : readVector(data, "/analysis/grid/massSecondDaughter")},
             zipLimits(readVector(data, "/analysis/limit/expected"),
                       readVector(data, "/analysis/limit/observed"))},
      constraints_{
          readConstraints(
              data.value("/constraints/mother"_json_pointer, nlohmann::json{}),
              collider()),
          readConstraints(data.value("/constraints/firstDaughter"_json_pointer,
                                     nlohmann::json{}),
                          collider()),
          readConstraints(data.value("/constraints/secondDaughter"_json_pointer,
                                     nlohmann::json{}),
                          collider())} {}

PairProductionLimit::PairProductionLimit(const nlohmann::json &data,
                                         const std::string &loadedFrom,
                                         const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolutions_{
          readMassRes(data, "/analysis/massResolution/firstParticle"),
          readMassRes(data, "/analysis/massResolution/secondParticle")},
      process_{predictions::readPairProductionProcess(data.at("process"),
                                                      collider())},
      equalParticleMasses_{utilities::readIfPresent<bool>(
          data, "/analysis/equalParticleMasses", false)},
      limit_{{readVector(data, "/analysis/grid/massFirstParticle"),
              equalParticleMasses_
                  ? std::vector{0.}
                  : readVector(data, "/analysis/grid/massSecondParticle")},
             zipLimits(readVector(data, "/analysis/limit/expected"),
                       readVector(data, "/analysis/limit/observed"))},
      constraints_{
          readConstraints(data.value("/constraints/firstParticle"_json_pointer,
                                     nlohmann::json{}),
                          collider()),
          readConstraints(data.value("/constraints/secondParticle"_json_pointer,
                                     nlohmann::json{}),
                          collider())} {}

std::string ChannelLimit::processDesc() const noexcept {
    return process_.to_string();
}

std::string ChannelLimit::extentDesc() const noexcept {
    return fmt::format("M={}", limit_.extent().front());
}

std::string ChainDecayLimit::processDesc() const noexcept {
    return process_.to_string();
}
std::string ChainDecayLimit::extentDesc() const noexcept {
    auto extent = limit_.extent();
    return fmt::format("M1={}, M2={}", extent[0], extent[1]);
}

std::string PairDecayLimit::processDesc() const noexcept {
    return process_.to_string();
}
std::string PairDecayLimit::extentDesc() const noexcept {
    auto extent = limit_.extent();
    if (equalDaughterMasses_) {
        return fmt::format("M1={}, M2=M3={}", extent[0], extent[1]);
    } else {
        return fmt::format("M1={}, M2={}, M3={}", extent[0], extent[1],
                           extent[2]);
    }
}

std::string PairProductionLimit::processDesc() const noexcept {
    return process_.to_string();
}
std::string PairProductionLimit::extentDesc() const noexcept {
    auto extent = limit_.extent();
    if (equalParticleMasses_) {
        return fmt::format("M1=M2={}", extent[0]);
    } else {
        return fmt::format("M1={}, M2={}", extent[0], extent[1]);
    }
}

std::vector<AppliedLimit>
ChannelLimit::apply(const predictions::Predictions &prediction) const {
    auto log = logger();
    log->trace("Applying limit {}: {} {}", id(), processDesc(), extentDesc());

    const auto particleRate = [this](const predictions::Particle &p) {
        if (acceptances_) {
            return process_(p, acceptances_.value()(p));
        }
        return process_(p);
    };

    const auto appMassRange =
        applicableMassRanges(limit_.extent(), massResolution_);
    auto relevantParticles = predictions::Clustering::relevantParticles(
        process_, prediction, appMassRange, options().applicableMassUnc);
    if (normalization_) {
        utilities::eraseFromSetIf(std::get<0>(relevantParticles),
                                  [this](const predictions::Particle &p) {
                                      return !((*normalization_)(p.mass()) > 0);
                                  });
    }

    const auto [cluster] = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolution_, options().clusterMassUnc,
        constraints_, prediction);

    const auto applyToCluster = [this, &log, &particleRate](auto &&c) {
        const auto mass =
            predictions::Clustering::rateWeightedMasses(particleRate, c);
        auto rate = predictions::Clustering::combinedRate(particleRate, c);
        if (!normalization_) {
            log->trace("Cluster of mass {} has a rate of {} pb",
                       std::get<0>(mass).mass, rate);
        } else {
            // acceptances are not yet implemented for the normalization
            rate /= (*normalization_)(std::get<0>(mass).mass);
            log->trace("Cluster of mass {} has a mu of {}",
                       std::get<0>(mass).mass, rate);
        }
        auto [expRatio, obsRatio] = rate / getLimit(limit_, mass);
        return AppliedLimit{
            shared_from_this(), obsRatio, expRatio,
            predictions::ChannelProcess::contributingParticles(c)};
    };

    return cluster | ranges::views::transform(applyToCluster) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}

std::vector<AppliedLimit>
ChainDecayLimit::apply(const predictions::Predictions &prediction) const {
    auto log = logger();
    log->trace("Applying limit {}: {} {}", id(), processDesc(), extentDesc());

    const auto particleRate = [this](const predictions::Particle &mother,
                                     const predictions::Particle &daughter) {
        if (productionAcceptances_) {
            return process_(mother, daughter,
                            productionAcceptances_.value()(mother));
        }
        return process_(mother, daughter);
    };

    const auto appMassRanges =
        applicableMassRanges(limit_.extent(), massResolutions_);
    const auto relevantParticles = predictions::Clustering::relevantParticles(
        particleRate, prediction, appMassRanges, options().applicableMassUnc);

    const auto clustersForRoles = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolutions_, options().clusterMassUnc,
        constraints_, prediction);

    const auto applyToClusters = [particleRate, this](auto &&clusterSet) {
        auto masses =
            std::apply(predictions::Clustering::weightedMassesFor(particleRate),
                       clusterSet);
        const auto rate = std::apply(
            predictions::Clustering::combinedRateFor(particleRate), clusterSet);
        logger()->trace("Cluster of masses {} has a rate of {}", masses, rate);
        auto [expRatio, obsRatio] = rate / getLimit(limit_, masses);
        return AppliedLimit{
            shared_from_this(), obsRatio, expRatio,
            std::apply(predictions::ChainDecayProcess::contributingParticles,
                       clusterSet)};
    };

    return std::apply(ranges::views::cartesian_product, clustersForRoles) |
           ranges::views::transform(applyToClusters) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}

std::vector<AppliedLimit>
PairDecayLimit::apply(const predictions::Predictions &prediction) const {
    auto log = logger();
    log->trace("Applying limit {}: {} {}", id(), processDesc(), extentDesc());

    auto extent = limit_.extent();
    if (equalDaughterMasses_) {
        extent[2] = extent[1];
    }
    const auto appMassRanges = applicableMassRanges(extent, massResolutions_);
    const auto relevantParticles = predictions::Clustering::relevantParticles(
        process_, prediction, appMassRanges, options().applicableMassUnc);

    const auto clustersForRoles = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolutions_, options().clusterMassUnc,
        constraints_, prediction);

    const auto validDaughterMasses = [this](auto &&clusterSets) {
        if (equalDaughterMasses_) {
            const auto masses =
                std::apply(predictions::Clustering::weightedMassesFor(process_),
                           clusterSets);
            return withinExpRes(std::get<1>(masses), std::get<2>(masses),
                                std::get<1>(massResolutions_),
                                std::get<2>(massResolutions_));
        };
        return true;
    };

    const auto applyToClusters = [this](auto &&clusterSets) {
        auto masses = std::apply(
            predictions::Clustering::weightedMassesFor(process_), clusterSets);
        const auto rate = std::apply(
            predictions::Clustering::combinedRateFor(process_), clusterSets);
        logger()->trace("Cluster of masses {} has a rate of {}", masses, rate);
        auto [expRatio, obsRatio] = rate / getLimit(limit_, masses);
        return AppliedLimit{
            shared_from_this(), obsRatio, expRatio,
            std::apply(predictions::PairDecayProcess::contributingParticles,
                       clusterSets)};
    };

    return std::apply(ranges::views::cartesian_product, clustersForRoles) |
           ranges::views::filter(validDaughterMasses) |
           ranges::views::transform(applyToClusters) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}

std::vector<AppliedLimit>
PairProductionLimit::apply(const predictions::Predictions &prediction) const {
    auto log = logger();
    log->trace("Applying limit {}: {} {}", id(), processDesc(), extentDesc());

    const auto processWithContext =
        [this, &prediction](const predictions::Particle &p1,
                            const predictions::Particle &p2) {
            return process_(prediction, p1, p2);
        };

    auto extent = limit_.extent();
    if (equalParticleMasses_) {
        extent[1] = extent[0];
    }
    const auto appMassRanges = applicableMassRanges(extent, massResolutions_);
    const auto relevantParticles = predictions::Clustering::relevantParticles(
        processWithContext, prediction, appMassRanges,
        options().applicableMassUnc);

    const auto clustersForRoles = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolutions_, options().clusterMassUnc,
        constraints_, prediction);

    const auto validParticleMasses = [this,
                                      &processWithContext](auto &&clusterSets) {
        if (equalParticleMasses_) {
            const auto masses = std::apply(
                predictions::Clustering::weightedMassesFor(processWithContext),
                clusterSets);
            return withinExpRes(std::get<0>(masses), std::get<1>(masses),
                                std::get<0>(massResolutions_),
                                std::get<1>(massResolutions_));
        }
        return true;
    };

    const auto applyToClusters = [this,
                                  &processWithContext](auto &&clusterSets) {
        auto masses = std::apply(
            predictions::Clustering::weightedMassesFor(processWithContext),
            clusterSets);
        const auto rate = std::apply(
            predictions::Clustering::combinedRateFor(processWithContext),
            clusterSets);
        logger()->trace("Cluster of masses {} has a rate of {}", masses, rate);
        auto [expRatio, obsRatio] = rate / getLimit(limit_, masses);
        return AppliedLimit{
            shared_from_this(), obsRatio, expRatio,
            std::apply(
                predictions::PairProductionProcess::contributingParticles,
                clusterSets)};
    };

    return std::apply(ranges::views::cartesian_product, clustersForRoles) |
           ranges::views::filter(validParticleMasses) |
           ranges::views::transform(applyToClusters) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}
} // namespace Higgs::bounds
