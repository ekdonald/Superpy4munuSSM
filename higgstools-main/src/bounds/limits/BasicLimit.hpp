/**
 * @file BasicLimit.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Common limit metadata
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "predictions/UncertainMass.hpp"
#include "utilities/ArithmeticArray.hpp"
#include "utilities/JsonFwd.hpp"
#include "utilities/LinearInterpolator.hpp"
#include "utilities/Logging.hpp"
#include <array>
#include <cstddef>
#include <memory>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/algorithm/min.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/view/filter.hpp>
#include <string>
#include <utility>
#include <vector>

namespace Higgs::predictions {
struct MassResolution;
}
namespace Higgs::bounds {

/**
 * @brief Intermediate limit class that handles the metadata associated with
 * each limit.
 *
 * Since it does not define the apply method this class is still abstract.
 *
 */
class BasicLimit : public Limit {
  public:
    unsigned id() const noexcept override;

    const std::string &reference() const noexcept override;

    const std::string &citeKey() const noexcept override;

    predictions::Collider collider() const noexcept override;

    predictions::Experiment experiment() const noexcept override;

    double luminosity() const noexcept override;

    const std::string &loadedFrom() const noexcept override;

    //! Return the applicable mass ranges from the limit extent, massResolutions
    //! and the applicableResolutionFac option for each role in the limit
    //! process.
    template <std::size_t nParticleRoles, std::size_t nParams>
    std::array<std::pair<double, double>, nParticleRoles> applicableMassRanges(
        const std::array<std::pair<double, double>, nParams> &extents,
        const std::array<predictions::MassResolution, nParticleRoles>
            &massResolutions) const noexcept {
        static_assert(nParticleRoles <= nParams);
        return applicableMassRangesImpl(
            extents, massResolutions,
            std::make_index_sequence<nParticleRoles>{});
    }

    //! Compute the applicable mass range of a limit from it's extent,
    //! massResolution and the applicableResolutionFac option.
    std::pair<double, double> applicableMassRange(
        const std::pair<double, double> &extent,
        const predictions::MassResolution &massResolution) const noexcept;

    //! Check if the two masses with theoretical uncertainties are close to
    //! each other given the corresponding experimental resolutions. Uses
    //! LimitOptions::applicableMassUnc to determine the treatment of the
    //! theoretical uncertainty.
    bool withinExpRes(predictions::UncertainMass m1,
                      predictions::UncertainMass m2,
                      const predictions::MassResolution &resM1,
                      const predictions::MassResolution &resM2) const noexcept;

    //! access to the options for this limit
    const LimitOptions &options() const noexcept override;

    using ExpObsLim = utilities::ArithmeticArray<double, 2>;

    //! Gets the limit from the interolator taking care of the eagerness.
    template <std::size_t nMasses, std::size_t nAddPars = 0>
    ExpObsLim getLimit(
        const utilities::LinearInterpolator<nMasses + nAddPars, ExpObsLim> &lim,
        const std::array<predictions::UncertainMass, nMasses> &masses,
        const std::array<double, nAddPars> &additionalPars = {})
        const noexcept {
        const auto pars = [&additionalPars, &masses](const auto &massFunc) {
            auto res = std::array<double, nMasses + nAddPars>{};
            auto endMasses = ranges::transform(masses, res.begin(), massFunc);
            ranges::copy(additionalPars, endMasses.out);
            return res;
        };
        switch (options_.setLimitMassUnc) {
        case predictions::MassUncEagerness::ignore:
            return lim(pars(predictions::centralMass));
        default:
            logger()->error("Unknown MassUncEagerness in "
                            "BasicLimit::getLimit(), defaulting to cautious.");
            [[fallthrough]];
        case predictions::MassUncEagerness::cautious:
            return lim.maxWithin(pars(predictions::lowerMassBound),
                                 pars(predictions::upperMassBound), std::less{},
                                 getObsLimit);
        case predictions::MassUncEagerness::eager:
            return lim.minWithin(pars(predictions::lowerMassBound),
                                 pars(predictions::upperMassBound), std::less{},
                                 getObsLimit);
        }
    }

    static std::vector<ExpObsLim>
    zipLimits(const std::vector<double> &expected,
              const std::vector<double> &observed);

    //! Returns a filter that removes all appliedLimits with expRatios <=
    //! options.minExpRatio().
    auto filterRelevantAppLimits() const noexcept {
        return ranges::views::filter(
            [min = options_.minExpRatio](const AppliedLimit &appLim) {
                return appLim.expRatio() > min;
            });
    }

  protected:
    //! constructor from explicitely specified metadata
    BasicLimit(unsigned id, std::string reference, std::string citeKey,
               predictions::Collider collider,
               predictions::Experiment experiment, double luminosity,
               std::string loadedFrom, const LimitOptions &options);

    //! constructor that reads (most) metadata from json data
    BasicLimit(const nlohmann::json &data, const std::string &loadedFrom,
               const LimitOptions &options);

  private:
    //! Implementation for applicableMassRanges.
    template <std::size_t nParticleRoles, std::size_t nParams,
              std::size_t... Is>
    auto applicableMassRangesImpl(
        const std::array<std::pair<double, double>, nParams> &extents,
        const std::array<predictions::MassResolution, nParticleRoles>
            &massResolutions,
        std::index_sequence<Is...>) const noexcept {
        return std::array{applicableMassRange(
            std::get<Is>(extents), std::get<Is>(massResolutions))...};
    }
    static constexpr auto getObsLimit = [](const ExpObsLim &lim) {
        return std::get<1>(lim);
    };

    unsigned id_;
    std::string reference_;
    std::string citeKey_;
    predictions::Collider collider_;
    predictions::Experiment experiment_;
    double luminosity_;
    std::string file_;
    LimitOptions options_;
};
} // namespace Higgs::bounds
