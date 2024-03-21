#include "predictions/UncertainMass.hpp"
#include "utilities/Logging.hpp"
#include <cmath> // IWYU pragma: keep
#include <memory>
#include <tuple>

namespace Higgs::predictions {
bool UncertainMass::liesWithin(const std::pair<double, double> &interval,
                               MassUncEagerness eagerness) const noexcept {
    switch (eagerness) {
    case MassUncEagerness::eager:
        return interval.first <= mass + uncertainty &&
               mass - uncertainty <= interval.second;
    case MassUncEagerness::cautious:
        return interval.first <= mass - uncertainty &&
               mass + uncertainty <= interval.second;
    case MassUncEagerness::ignore:
        return interval.first <= mass && mass <= interval.second;
    }
    logger()->error(
        "Unknown value for MassUncEagerness in UncertainMass::liesWithin()");
    return false;
}

bool operator==(const UncertainMass &a, const UncertainMass &b) noexcept {
    return std::tie(a.mass, a.uncertainty) == std::tie(b.mass, b.uncertainty);
}

bool operator!=(const UncertainMass &a, const UncertainMass &b) noexcept {
    return !(a == b);
}

UncertainMass operator-(UncertainMass x) noexcept {
    x.mass = -x.mass;
    return x;
}

UncertainMass &operator+=(UncertainMass &x,
                          const UncertainMass &other) noexcept {
    x.mass += other.mass;
    x.uncertainty += other.uncertainty;
    return x;
}

UncertainMass &operator*=(UncertainMass &x, double s) noexcept {
    x.mass *= s;
    x.uncertainty *= std::abs(s);
    return x;
}

UncertainMass &operator/=(UncertainMass &x, double s) noexcept {
    x.mass /= s;
    x.uncertainty /= std::abs(s);
    return x;
}

UncertainMass operator+(UncertainMass a, const UncertainMass &b) noexcept {
    return a += b;
}

UncertainMass operator*(UncertainMass a, double s) noexcept { return a *= s; }
UncertainMass operator*(double s, UncertainMass a) noexcept { return a *= s; }
UncertainMass operator/(UncertainMass a, double s) noexcept { return a /= s; }
} // namespace Higgs::predictions
