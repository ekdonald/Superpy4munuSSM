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
} // namespace Higgs::predictions

namespace Higgs::bounds {

class ChannelLimit : public BasicLimit {
  public:
    //! factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new ChannelLimit(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &predictions) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

  protected:
    //! constructor that reads from the json data
    ChannelLimit(const nlohmann::json &data, const std::string &loadedFrom,
                 const LimitOptions &options = LimitOptions{});

  private:
    static constexpr std::size_t nParticleRoles = 1;
    std::array<predictions::MassResolution, nParticleRoles> massResolution_;
    predictions::ChannelProcess process_;
    utilities::LinearInterpolator<nParticleRoles, ExpObsLim> limit_;
    std::optional<predictions::ReferenceRate> normalization_;
    std::array<std::vector<predictions::Constraint>, nParticleRoles>
        constraints_;
    std::optional<Acceptances> acceptances_;
};

class ChainDecayLimit : public BasicLimit {
  public:
    //! factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new ChainDecayLimit(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &predictions) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

  protected:
    //! constructor that reads from the json data
    ChainDecayLimit(const nlohmann::json &data, const std::string &loadedFrom,
                    const LimitOptions &options = LimitOptions{});

  private:
    // { motherParticle, daughterParticle }
    static constexpr std::size_t nParticleRoles = 2;
    std::array<predictions::MassResolution, nParticleRoles> massResolutions_;
    predictions::ChainDecayProcess process_;
    utilities::LinearInterpolator<nParticleRoles, ExpObsLim> limit_;
    std::array<std::vector<predictions::Constraint>, nParticleRoles>
        constraints_;
    std::optional<Acceptances> productionAcceptances_;
};

class PairDecayLimit : public BasicLimit {
  public:
    //! factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new PairDecayLimit(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &prediction) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

  protected:
    //! constructor that reads from the json data
    PairDecayLimit(const nlohmann::json &data, const std::string &loadedFrom,
                   const LimitOptions &options = LimitOptions{});

  private:
    // {mother, first daughter, second daughter }
    static constexpr std::size_t nParticleRoles = 3;
    std::array<predictions::MassResolution, nParticleRoles> massResolutions_;
    predictions::PairDecayProcess process_;
    bool equalDaughterMasses_;
    utilities::LinearInterpolator<nParticleRoles, ExpObsLim> limit_;
    std::array<std::vector<predictions::Constraint>, nParticleRoles>
        constraints_;
};

class PairProductionLimit : public BasicLimit {
  public:
    //! factory method that acts as a constructor
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new PairProductionLimit(std::forward<Args>(args)...));
    } // LCOV_EXCL_LINE

    std::vector<AppliedLimit>
    apply(const predictions::Predictions &prediction) const override;

    std::string processDesc() const noexcept override;
    std::string extentDesc() const noexcept override;

  protected:
    //! constructor that reads from the json data
    PairProductionLimit(const nlohmann::json &data,
                        const std::string &loadedFrom,
                        const LimitOptions &options = LimitOptions{});

  private:
    // { first particle, second particle }
    static constexpr std::size_t nParticleRoles = 2;
    std::array<predictions::MassResolution, nParticleRoles> massResolutions_;
    predictions::PairProductionProcess process_;
    bool equalParticleMasses_;
    utilities::LinearInterpolator<nParticleRoles, ExpObsLim> limit_;
    std::array<std::vector<predictions::Constraint>, nParticleRoles>
        constraints_;
};

} // namespace Higgs::bounds
