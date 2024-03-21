/**
 * @file JsonSupport.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief json support for types in Higgs::predictions
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Process.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "utilities/JsonFwd.hpp"

namespace Higgs::predictions {

HIGGSUTILITIES_JSON_CONV_FWD(Collider)
HIGGSUTILITIES_JSON_CONV_FWD(ECharge)
HIGGSUTILITIES_JSON_CONV_FWD(Experiment)
HIGGSUTILITIES_JSON_CONV_FWD(Production)
HIGGSUTILITIES_JSON_CONV_FWD(Decay)
HIGGSUTILITIES_JSON_CONV_FWD(ChainDecay)
HIGGSUTILITIES_JSON_CONV_FWD(ReferenceModel)
HIGGSUTILITIES_JSON_CONV_FWD(CP)
HIGGSUTILITIES_JSON_CONV_FWD(Coupling)

HIGGSUTILITIES_JSON_CONV_FWD(MassResolution)

/**
 * @brief Read a ChannelProcess from json.
 * @param json the json to read from. This class is represented e.g. as
 * ```json
 * {
 *     "channels": [
 *         ["ggH", "gamgam"],
 *         ["bbH","mumu"],
 *         ...
 *     ]
 * }
 * ```
 * Where the `"channels"` field is read into a vector of pairs `["production
 * mode","decay mode"]`.
 *
 * @param collider the collider at which the process takes place
 * @return ChannelProcess a newly constructed ChannelProcess
 */
ChannelProcess readChannelProcess(const nlohmann::json &json,
                                  Collider collider);

/**
 * @brief Read a ChainDecayProcess from json.
 * @param j the json to read from. This class is represented e.g. as
 * ```json
 * {
 *     "chain": "Z",
 *     "production": ["ggH","bbH"],
 *     "decay": ["bb","tautau"],
 *     "collider": "LHC13"
 * }
 * ```
 * The `production` and `decay` fields are arrays of Production and Decay
 * values, respectively. The `"collider"` field is optional if a corresponding
 * default value is given (see below).
 * @param collider the collider at which the process takes place
 */
ChainDecayProcess readChainDecayProcess(const nlohmann::json &j,
                                        Collider collider);

/**
 * @brief Read a PairDecayProcess from json.
 * @param j the json to read from. This class is represented e.g. as
 * ```json
 * {
 *     "production": [
 *         "vbfH"
 *     ],
 *     "firstDecay": [
 *         "taunu", "munu"
 *     ],
 *     "secondDecay": [
 *         "bb", "mumu"
 *     ]
 * }
 * ```
 * The `production`, `firstDecay` and `secondDecay` fields are arrays of
 * Production and Decay values, respectively.
 * @param collider the collider at which the process takes place
 */
PairDecayProcess readPairDecayProcess(const nlohmann::json &j,
                                      Collider collider);

/**
 * @brief Read a PairProductionProcess from json.
 * @param j the json to read from. This class is represented e.g. as
 * ```json
 * {
 *     "firstDecay": [
 *         "taunu", "munu"
 *     ],
 *     "secondDecay": [
 *         "bb", "mumu"
 *     ]
 * }
 * ```
 * The `firstDecay` and `secondDecay` fields are arrays Decay values.
 * @param collider the collider at which the process takes place
 */
PairProductionProcess readPairProductionProcess(const nlohmann::json &j,
                                                Collider collider);

} // namespace Higgs::predictions
