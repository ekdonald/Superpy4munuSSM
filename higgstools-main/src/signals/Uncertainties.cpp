#include "signals/Uncertainties.hpp"
#include "utilities/Format.hpp"
#include <cmath>
#include <stdexcept>

namespace Higgs::signals {

UncertainValue::UncertainValue(double central) noexcept : central_{central} {}
UncertainValue::UncertainValue(double central, double uncertainty) noexcept
    : central_{central}, lowerUnc_{std::abs(uncertainty)}, upperUnc_{std::abs(
                                                               uncertainty)} {}
UncertainValue::UncertainValue(double central, double lowerUnc,
                               double upperUnc) noexcept
    : central_{central}, lowerUnc_{std::abs(lowerUnc)}, upperUnc_{std::abs(
                                                            upperUnc)} {}
UncertainValue UncertainValue::fromInterval(double lowerBound, double central,
                                            double upperBound) {
    if (lowerBound > central || central > upperBound) {
        throw std::domain_error(fmt::format(
            "Invalid uncertainty interval [{},{}] for central value {}.",
            lowerBound, upperBound, central));
    }
    return UncertainValue{central, central - lowerBound, upperBound - central};
}

double UncertainValue::central() const noexcept { return central_; }
double UncertainValue::symmetricUncertainty() const noexcept {
    return (lowerUnc_ + upperUnc_) / 2.;
}
double UncertainValue::upperUncertainty() const noexcept { return upperUnc_; }
double UncertainValue::lowerUncertainty() const noexcept { return lowerUnc_; }

double UncertainValue::uncertaintyTowards(double target) const noexcept {
    return target >= central_ ? upperUnc_ : lowerUnc_;
}

bool UncertainValue::operator==(const UncertainValue &other) const noexcept {
    return central_ == other.central_ && lowerUnc_ == other.lowerUnc_ &&
           upperUnc_ == other.upperUnc_;
}

} // namespace Higgs::signals
