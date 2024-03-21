#include "bounds/limits/BasicLimit.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "predictions/UncertainMass.hpp"
#include "utilities/Json.hpp"
#include <memory>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/zip_with.hpp>
#include <stdexcept>
#include <utilities/Logging.hpp>
#include <utility>

namespace Higgs::bounds {

BasicLimit::BasicLimit(unsigned id, std::string reference, std::string citeKey,
                       predictions::Collider collider,
                       predictions::Experiment experiment, double luminosity,
                       std::string loadedFrom, const LimitOptions &options)
    : id_{id}, reference_{std::move(reference)}, citeKey_{std::move(citeKey)},
      collider_{collider}, experiment_{experiment},
      luminosity_{luminosity}, file_{std::move(loadedFrom)}, options_{options} {
}

BasicLimit::BasicLimit(const nlohmann::json &data,
                       const std::string &loadedFrom,
                       const LimitOptions &options)
    : BasicLimit{utilities::readAs<unsigned>(data, "id"),
                 utilities::readAs<std::string>(data, "reference"),
                 utilities::readAs<std::string>(data, "citeKey"),
                 utilities::readAs<predictions::Collider>(data, "collider"),
                 utilities::readAs<predictions::Experiment>(data, "experiment"),
                 utilities::readAs<double>(data, "luminosity"),
                 loadedFrom,
                 options} {}

unsigned BasicLimit::id() const noexcept { return id_; }

const std::string &BasicLimit::reference() const noexcept { return reference_; }

const std::string &BasicLimit::citeKey() const noexcept { return citeKey_; }

predictions::Collider BasicLimit::collider() const noexcept {
    return collider_;
}

predictions::Experiment BasicLimit::experiment() const noexcept {
    return experiment_;
}

double BasicLimit::luminosity() const noexcept { return luminosity_; }

const std::string &BasicLimit::loadedFrom() const noexcept { return file_; }

const LimitOptions &BasicLimit::options() const noexcept { return options_; }

std::pair<double, double> BasicLimit::applicableMassRange(
    const std::pair<double, double> &extent,
    const predictions::MassResolution &resolution) const noexcept {
    auto massRange = extent;
    massRange.first -=
        (extent.first * resolution.relative + resolution.absolute) *
        options_.applicableResolutionFac;
    massRange.second +=
        (extent.second * resolution.relative + resolution.absolute) *
        options_.applicableResolutionFac;
    return massRange;
}

bool BasicLimit::withinExpRes(
    predictions::UncertainMass m1, predictions::UncertainMass m2,
    const predictions::MassResolution &resM1,
    const predictions::MassResolution &resM2) const noexcept {
    if (m1.mass > m2.mass) {
        return withinExpRes(m2, m1, resM2, resM1);
    }

    auto m1ExpRes = resM1.absolute + resM1.relative * m1.mass;
    auto m2ExpRes = resM2.absolute + resM2.relative * m2.mass;

    switch (options_.applicableMassUnc) {
    case predictions::MassUncEagerness::eager:
        return upperMassBound(m1) + m1ExpRes >= lowerMassBound(m2) - m2ExpRes;
    case predictions::MassUncEagerness::ignore:
        return m1.mass + m1ExpRes >= m2.mass - m2ExpRes;
    case predictions::MassUncEagerness::cautious:
        return lowerMassBound(m1) + m1ExpRes >= upperMassBound(m2) - m2ExpRes;
    }
    logger()->error("Unknown enum value of MassUncEagerness in "
                    "BasicLimit::withinExpRes()");
    return false;
}

std::vector<BasicLimit::ExpObsLim>
BasicLimit::zipLimits(const std::vector<double> &expected,
                      const std::vector<double> &observed) {
    return ranges::views::zip_with(
               [](const auto &e, const auto &o) {
                   return ExpObsLim{{e, o}};
               },
               expected, observed) |
           ranges::to<std::vector>;
    ;
}

} // namespace Higgs::bounds
