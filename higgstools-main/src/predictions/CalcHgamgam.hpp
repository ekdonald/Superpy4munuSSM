/**
 * @file calcHgamgam.hpp
 * @author Henning Bahl (henning.bahl@uchicago.edu)
 *
 * @brief Calculation of H -> gamgam signal strength for Higgs::predictions
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsTools.
 * HiggsTools is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/EffectiveCouplings.hpp"

namespace Higgs::predictions {


// main function to calculate H -> gammagamma signal strength in terms of
// effective couplings (based on arXiv:1211.1980, Eqs. (10)-(12))
double kgamma2(double mass, const NeutralEffectiveCouplings &coups) noexcept;

} // namespace Higgs::predictions
