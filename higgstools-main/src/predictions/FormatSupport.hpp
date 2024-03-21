/**
 * @file FormatSupport.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Specializations of fmt::formatter for various HiggsPredictions types.
 *
 * @copyright Copyright 2021 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "utilities/Format.hpp" // IWYU pragma: export

HIGGSUTILITIES_ENUM_FORMATTER(Higgs::predictions::Experiment)
HIGGSUTILITIES_ENUM_FORMATTER(Higgs::predictions::ColliderType)
HIGGSUTILITIES_ENUM_FORMATTER(Higgs::predictions::Production)
HIGGSUTILITIES_ENUM_FORMATTER(Higgs::predictions::Decay)
HIGGSUTILITIES_ENUM_FORMATTER(Higgs::predictions::Collider)
