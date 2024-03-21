#include "bounds/limits/WidthLimits.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/Clustering.hpp"
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
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>

namespace Higgs::bounds {
namespace {
constexpr auto readVector = utilities::readAs<std::vector<double>>;
constexpr auto readMassRes =
    utilities::readIfPresent<predictions::MassResolution>;
} // namespace

ChannelWidthLimit::ChannelWidthLimit(const nlohmann::json &data,
                                     const std::string &loadedFrom,
                                     const LimitOptions &options)
    : BasicLimit{data, loadedFrom, options},
      massResolution_{readMassRes(data, "/analysis/massResolution")},
      process_{predictions::readChannelProcess(data.at("process"), collider())},
      limit_{{readVector(data, "/analysis/grid/mass"),
              readVector(data, "/analysis/grid/width")},
             zipLimits(readVector(data, "/analysis/limit/expected"),
                       readVector(data, "/analysis/limit/observed"))},
      relativeWidth_{utilities::readAs<bool>(data, "/analysis/relativeWidth")},
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

std::string ChannelWidthLimit::processDesc() const noexcept {
    return process_.to_string();
}

std::string ChannelWidthLimit::extentDesc() const noexcept {
    auto extent = limit_.extent();
    if (relativeWidth_) {
        return fmt::format("M={}, Gam/M={}", extent[0], extent[1]);
    }
    return fmt::format("M={}, Gam={}", extent[0], extent[1]);
}

std::vector<AppliedLimit>
ChannelWidthLimit::apply(const predictions::Predictions &prediction) const {
    auto log = logger();
    log->trace("Applying limit {}: {} {}", id(), processDesc(), extentDesc());

    auto particleRate = [this](const predictions::Particle &p) {
        if (acceptances_) {
            return process_(p, acceptances_.value()(p));
        }
        return process_(p);
    };

    const auto appMassRange =
        applicableMassRanges(limit_.extent(), massResolution_);
    auto relevantParticles = predictions::Clustering::relevantParticles(
        particleRate, prediction, appMassRange, options().applicableMassUnc);

    const auto widthOutOfBounds = [this](const predictions::Particle &p) {
        const auto width =
            relativeWidth_ ? p.totalWidth() / p.mass() : p.totalWidth();
        const auto [minWidth, maxWidth] = limit_.extent().back();
        return !((1 - widthExtentTol) * minWidth <= width &&
                 width <= (1 + widthExtentTol) * maxWidth);
    };

    utilities::eraseFromSetIf(relevantParticles[0], widthOutOfBounds);
    if (normalization_) {
        utilities::eraseFromSetIf(relevantParticles[0],
                                  [this](const predictions::Particle &p) {
                                      return !((*normalization_)(p.mass()) > 0);
                                  });
    }

    const auto [cluster] = predictions::Clustering::performClusteringNew(
        relevantParticles, massResolution_, options().clusterMassUnc,
        constraints_, prediction);

    const auto applyToParticle = [this, &log, &particleRate](auto &&c) {
        const auto mass =
            predictions::Clustering::rateWeightedMasses(particleRate, c);
        auto width =
            predictions::Clustering::rateWeightedWidths(particleRate, c);
        if (relativeWidth_) {
            std::get<0>(width) /= std::get<0>(mass).mass;
        }
        auto rate = predictions::Clustering::combinedRate(particleRate, c);
        if (!normalization_) {
            log->trace("Cluster of mass {} and width {} (relative={}) has a "
                       "rate of {}",
                       std::get<0>(mass).mass, std::get<0>(width),
                       relativeWidth_, rate);
        } else {
            rate /= (*normalization_)(std::get<0>(mass).mass);
            log->trace(
                "Cluster of mass {} and width {} (relative={}) has a mu of {}",
                std::get<0>(mass).mass, std::get<0>(width), relativeWidth_,
                rate);
        }
        auto [expRatio, obsRatio] = rate / getLimit(limit_, mass, width);
        return AppliedLimit{
            shared_from_this(), obsRatio, expRatio,
            predictions::ChannelProcess::contributingParticles(c)};
    };

    return cluster | ranges::views::transform(applyToParticle) |
           filterRelevantAppLimits() | ranges::to<std::vector>;
}
} // namespace Higgs::bounds
