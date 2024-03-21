/**
 * @file Uncertainties.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Values with (gaussian) uncertainties.
 *
 * @copyright Copyright 2021 by the authors.
 * This file is part of HiggsTools.
 * HiggsTools is released under the GPLv3+.
 */
#pragma once

namespace Higgs::signals {
class UncertainValue {
  public:
    UncertainValue() = default;
    explicit UncertainValue(double central) noexcept;
    UncertainValue(double central, double uncertainty) noexcept;
    UncertainValue(double central, double lowerUnc, double upperUnc) noexcept;
    static UncertainValue fromInterval(double lowerBound, double central,
                                       double upperBound);

    double central() const noexcept;
    double symmetricUncertainty() const noexcept;
    double upperUncertainty() const noexcept;
    double lowerUncertainty() const noexcept;

    double uncertaintyTowards(double target) const noexcept;

    bool operator==(const UncertainValue &other) const noexcept;

  private:
    double central_ = 0;
    double lowerUnc_ = 0;
    double upperUnc_ = 0;
};

} // namespace Higgs::signals
