#include "bounds/Likelihood.hpp"
#include <cmath>
#include <limits>

namespace {
// complement of the CDF for the unit normal distribution
double ccdf(double x) { return std::erfc(x / std::sqrt(2.)) / 2.; }
} // namespace

namespace Higgs::bounds::Likelihood {

double expCLs(double qExp) {
    constexpr auto expCLb = 0.5;
    return ccdf(std::sqrt(qExp)) / expCLb;
}

double obsCLs(double qExp, double qObs) {
    // -- first case of eq. (65) and (66), numerically safe --
    if (qObs <= qExp) {
        return ccdf(std::sqrt(qObs)) / ccdf(std::sqrt(qObs) - std::sqrt(qExp));
    }

    // -- second case --
    // expansion of the ratio for qExp->0
    if (qExp <= std::numeric_limits<double>::min()) {
        return std::exp(-qObs / 2.);
    }

    const auto qq = (qObs - qExp) / (2 * std::sqrt(qExp));
    // handle div by zero if CLb underflows (around qq=35 for double)
    const auto clb = ccdf(qq);
    if (clb <= std::numeric_limits<double>::min()) {
        // series expansion of the ratio for large qq to O(1/qq)
        return std::exp(-qExp / 2. - std::sqrt(qExp) * qq) *
               (1 - std::sqrt(qExp) / qq);
    }
    // otherwise, use second case ratio
    return ccdf((qObs + qExp) / (2 * std::sqrt(qExp))) / clb;
}

} // namespace Higgs::bounds::Likelihood
