/**
 * @file Helpers.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Miscellaneous utility functions for Higgs::predictions
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "utilities/Logging.hpp"
#include <memory>
#include <vector>
// IWYU pragma: no_include <exception>

namespace Higgs::predictions {

//! how to clamp a value a grid
enum class Clamp {
    none,  //!< don't clamp
    upper, //!< only clamp large values
    lower, //!< only clamp small values
    both   //!< clamp in both directions
};

//! Mass is valid and within bounds of the grid; clamps if requested.
template <Clamp c>
inline bool validMassIn(double &mass, const std::vector<double> &grid) {
    if (mass < 0) {
        static auto log = logger();
        log->warn("Invalid mass value of {} GeV", mass);
        return false;
    }
    if (mass < grid.front()) {
        if constexpr (c == Clamp::lower || c == Clamp::both) {
            mass = grid.front();
        } else {
            return false;
        }
    } else if (mass > grid.back()) {
        if constexpr (c == Clamp::upper || c == Clamp::both) {
            mass = grid.back();
        } else {
            return false;
        }
    }
    return true;
}
} // namespace Higgs::predictions
