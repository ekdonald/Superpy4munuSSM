/**
 * @file LikelihoodLimit.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Limit given as an exclusion Likelihood for one or more processes
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Process.hpp"
#include "bounds/limits/BasicLimit.hpp"
#include "predictions/Constraints.hpp"
#include "utilities/JsonFwd.hpp"
#include "utilities/LinearInterpolator.hpp"
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Higgs::predictions {
class Predictions;
} // namespace Higgs::predictions

namespace Higgs::bounds {

class LikelihoodLimit1d : public BasicLimit {
  public:
    //! this is a 1d likelihood
    constexpr static std::size_t llhDim = 1;
    // !factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new LikelihoodLimit1d(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &prediction) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

  protected:
    LikelihoodLimit1d(const nlohmann::json &data, const std::string &loadedFrom,
                      const LimitOptions &options = LimitOptions{});

  private:
    static constexpr std::size_t nParticles = 1;
    std::array<predictions::MassResolution, nParticles> massResolutions_;
    std::array<predictions::ChannelProcess, llhDim> processes_;
    utilities::LinearInterpolator<1 + llhDim> obs_;
    utilities::LinearInterpolator<1 + llhDim> exp_;
    std::array<std::vector<predictions::Constraint>, nParticles> constraints_;
};

class LikelihoodLimit2d : public BasicLimit {
  public:
    //! this is a 2d likelihood
    constexpr static std::size_t llhDim = 2;
    // !factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new LikelihoodLimit2d(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &prediction) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

  protected:
    LikelihoodLimit2d(const nlohmann::json &data, const std::string &loadedFrom,
                      const LimitOptions &options = LimitOptions{});

  private:
    static constexpr std::size_t nParticles = 1;
    std::array<predictions::MassResolution, nParticles> massResolutions_;
    std::array<predictions::ChannelProcess, 2> processes_;
    predictions::ChannelProcess mergedProcess_;
    utilities::LinearInterpolator<1, utilities::LinearInterpolator<llhDim>>
        obs_;
    utilities::LinearInterpolator<1, utilities::LinearInterpolator<llhDim>>
        exp_;
    std::array<std::vector<predictions::Constraint>, nParticles> constraints_;
};

} // namespace Higgs::bounds
