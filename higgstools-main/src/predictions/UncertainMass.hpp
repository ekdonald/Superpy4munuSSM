#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "utilities/Format.hpp"
#include <utility>

namespace Higgs::predictions {

//! A mass with a (symmetric) theory uncertainty, `mass +- uncertainty`. The
//! uncertainty is **not** treated as a statistical error and is simply
//! propagated as an interval. Supports addition as well as
//! multiplication and division with a scalar. Multiplication and division of
//! two uncertain masses are not supported as they would lead to
//! asymmetric uncertainties. Subtraction is not supported since interval
//! subtraction can break some algorithms.
struct UncertainMass {
    double mass = 0.;        //!< mass in GeV
    double uncertainty = 0.; //!< uncertainty in GeV

    //! Check if the mass lies inside of an interval.
    bool liesWithin(const std::pair<double, double> &interval,
                    MassUncEagerness eagerness) const noexcept;
};

inline constexpr double centralMass(const UncertainMass &m) { return m.mass; }
inline constexpr double lowerMassBound(const UncertainMass &m) {
    return m.mass - m.uncertainty;
}
inline constexpr double upperMassBound(const UncertainMass &m) {
    return m.mass + m.uncertainty;
}

bool operator==(const UncertainMass &a, const UncertainMass &b) noexcept;
bool operator!=(const UncertainMass &a, const UncertainMass &b) noexcept;

UncertainMass operator-(UncertainMass x) noexcept;

UncertainMass &operator+=(UncertainMass &x,
                          const UncertainMass &other) noexcept;
UncertainMass &operator*=(UncertainMass &x, double s) noexcept;
UncertainMass &operator/=(UncertainMass &x, double s) noexcept;

UncertainMass operator+(UncertainMass a, const UncertainMass &b) noexcept;
UncertainMass operator*(UncertainMass a, double s) noexcept;
UncertainMass operator*(double s, UncertainMass a) noexcept;
UncertainMass operator/(UncertainMass a, double s) noexcept;

} // namespace Higgs::predictions

namespace fmt {
template <>
struct formatter<Higgs::predictions::UncertainMass> // LCOV_EXCL_LINE
    : formatter<std::pair<double, double>> {
    template <typename FormatContext>
    auto format(Higgs::predictions::UncertainMass e, FormatContext &ctx) {
        return formatter<std::pair<double, double>>::format(
            std::pair{e.mass, e.uncertainty}, ctx);
    }
};
} // namespace fmt
