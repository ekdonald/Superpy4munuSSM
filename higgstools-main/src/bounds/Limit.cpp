#include "Higgs/bounds/Limit.hpp"
#include "limits/LikelihoodLimit.hpp"
#include "limits/ProcessLimits.hpp"
#include "limits/WidthLimits.hpp"
#include "predictions/FormatSupport.hpp"
#include "utilities/Json.hpp"
#include <fstream>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <utility>

namespace Higgs::bounds {

//! the different kinds of limits supported by HiggsBounds
enum class LimitClass {
    ChannelLimit,        //!< generic limit with one intermediate BSM resonance
    ChannelWidthLimit,   //!< width dependent version of the ChannelLimit
    LikelihoodLimit,     //!< limit given in terms of an exclusion likelihood
    ChainDecayLimit,     //!< limit with two BSM particles as Phi1 > Phi2 + ...
    PairDecayLimit,      //!< limit with three BSM particles as Phi1 > Phi2 Phi3
    PairProductionLimit, //!< limit with two pair-produced BSM particles
};
HIGGSUTILITIES_JSON_CONV_FWD(
    LimitClass) // easiest way to silence -Wmissing-declarations
HIGGSUTILITIES_ENUM_JSON_CONV(LimitClass)

std::shared_ptr<Limit> Limit::read(const std::string &filePath,
                                   const LimitOptions &options) {
    const auto j = nlohmann::json::parse(std::ifstream{filePath});
    switch (utilities::readAs<LimitClass>(j, "limitClass")) {
    case LimitClass::ChannelLimit:
        return ChannelLimit::create(j, filePath, options);
    case LimitClass::ChannelWidthLimit:
        return ChannelWidthLimit::create(j, filePath, options);
    case LimitClass::ChainDecayLimit:
        return ChainDecayLimit::create(j, filePath, options);
    case LimitClass::PairDecayLimit:
        return PairDecayLimit::create(j, filePath, options);
    case LimitClass::PairProductionLimit:
        return PairProductionLimit::create(j, filePath, options);
    case LimitClass::LikelihoodLimit:
        const auto n = j.at("process").size();
        switch (n) {
        case 1UL:
            return LikelihoodLimit1d::create(j, filePath, options);
        case 2UL:
            return LikelihoodLimit2d::create(j, filePath, options);
        default:
            throw std::runtime_error(
                fmt::format("The requested {}-dim LikelihoodLimit is not yet "
                            "implemented.",
                            n));
        }
    }
    throw std::runtime_error("unreacheable"); // LCOV_EXCL_LINE
}

std::string Limit::to_string() const noexcept {
    if (luminosity() > 0) {
        return fmt::format("{} from {} ({} {}fb-1, {})", processDesc(),
                           reference(), experiment(), luminosity(),
                           extentDesc());
    }
    return fmt::format("{} from {} ({}, {})", processDesc(), reference(),
                       experiment(), extentDesc());
}

std::shared_ptr<const Limit> AppliedLimit::limit() const noexcept {
    return limit_;
}
double AppliedLimit::expRatio() const noexcept { return expRatio_; }
double AppliedLimit::obsRatio() const noexcept { return obsRatio_; }
double AppliedLimit::obsLikelihood() const noexcept { return obsLikelihood_; }
double AppliedLimit::expLikelihood() const noexcept { return expLikelihood_; }
const std::vector<std::string> &
AppliedLimit::contributingParticles() const noexcept {
    return contributingParticles_;
}

AppliedLimit::AppliedLimit(std::shared_ptr<const Limit> limit, double obsRatio,
                           double expRatio,
                           std::vector<std::string> contributingParticles,
                           double obsLikelihood, double expLikelihood)
    : limit_{std::move(limit)}, obsRatio_{obsRatio}, expRatio_{expRatio},
      contributingParticles_{std::move(contributingParticles)},
      obsLikelihood_{obsLikelihood}, expLikelihood_{expLikelihood} {}

} // namespace Higgs::bounds
