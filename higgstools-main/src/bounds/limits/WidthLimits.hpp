/**
 * @file WidthLimits.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Width-dependent limits on a process
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Process.hpp"
#include "bounds/Acceptance.hpp"
#include "bounds/limits/BasicLimit.hpp"
#include "predictions/Constraints.hpp"
#include "utilities/JsonFwd.hpp"
#include "utilities/LinearInterpolator.hpp"
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Higgs::predictions {
class Predictions;
}

namespace Higgs::bounds {

//! A variant of the ChannelLimit with additional total-width dependence.
//! There is no clustering, since width dependent analyses typically assume that
//! there is only one particle.
class ChannelWidthLimit : public BasicLimit {
  public:
    //! factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new ChannelWidthLimit(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &predictions) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

    //! Limit is applied for widths in
    //! `[(1-widthExtentTol) * minWidth, (1 + widthExtentTol) * maxWidth]`
    static constexpr double widthExtentTol = .1;

  protected:
    //! constructor that reads from the json data
    ChannelWidthLimit(const nlohmann::json &data, const std::string &loadedFrom,
                      const LimitOptions &options = LimitOptions{});

  private:
    static constexpr std::size_t nParticleRoles = 1;
    std::array<predictions::MassResolution, nParticleRoles> massResolution_;
    predictions::ChannelProcess process_;
    utilities::LinearInterpolator<2 * nParticleRoles, ExpObsLim> limit_;
    bool relativeWidth_;
    std::optional<predictions::ReferenceRate> normalization_;
    std::array<std::vector<predictions::Constraint>, nParticleRoles>
        constraints_;
    std::optional<Acceptances> acceptances_;
};
} // namespace Higgs::bounds
