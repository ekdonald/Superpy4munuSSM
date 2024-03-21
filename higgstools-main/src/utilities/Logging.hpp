/**
 * @file Logging.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Loggers used throughout the library
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include <fmt/format.h> // IWYU pragma: export
#include <fmt/ranges.h> // IWYU pragma: export
#include <memory>
#include <spdlog/logger.h> // IWYU pragma: export
#include <spdlog/spdlog.h> // IWYU pragma: export

namespace Higgs {
namespace predictions {
//! return the HiggsPredictions logger
std::shared_ptr<spdlog::logger> logger();
} // namespace predictions
namespace bounds {
//! return the HiggsBounds logger
std::shared_ptr<spdlog::logger> logger();
} // namespace bounds
namespace signals {
//! return the HiggsSignals logger
std::shared_ptr<spdlog::logger> logger();
} // namespace signals
} // namespace Higgs
