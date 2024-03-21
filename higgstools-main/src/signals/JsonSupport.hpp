/**
 * @file JsonSupport.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Conversion functions to and from json for HiggsSignals types.
 *
 * @copyright Copyright 2021 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "utilities/JsonFwd.hpp"
#include <Eigen/Core>
#include <memory>
#include <string>
#include <vector>

namespace Higgs::predictions {
class Particle;
}

namespace Higgs::signals {

class SubMeasurement;
class UncertainValue;

HIGGSUTILITIES_JSON_CONV_FWD(UncertainValue)

/**
 * @brief Read a RateMeasurement from json. The json representation is
 * ```json
 * {
 *      "process": {},
 *      "obs": [0.5,1.,1.5],
 *      "exp": [0.7,0.9,1.2],
 *      "channelWeights": [1.,2.],
 *      "massResolution": 2.5,
 *      "massSensitive": true
 * }
 * ```
 * where the process field is read into a predictions::ChannelProcess using
 * predictions::readChannelProcess. The three values in the rate fields
 * represent
 * ```json
 * ["lower uncertainty interval", "central value", "upper uncertainty
 * interval"]
 * ```
 * The `"exp"` field is optional. If it it not present the expected rate
 * will be the rate of the reference particle for the given process and
 * channel weights: `process(*referenceParticle, channelWeights)`. The
 * `"channelWeights"` field is an array of the same size as
 * `/process/channels` that gives multiplicative weights for each channel.
 * It is *optional* and if it is not present all weights default to 1. The
 * `"massResolution"` field is *optional* as well and can be specified to
 * override the `defaultMassResolution`. Finally, the
 * `"massSensitive"` flag marks the rate measurement as mass sensitive, it
 * is *optional* and defaults to `false` if it is not specified.
 *
 * @param json the json data to read from
 * @param collider the collider at which the "process" takes place, passed
 * to predictions::readChannelProcess
 * @param referenceParticle the reference particle at the reference mass for
 * the measurement
 * @param defaultMassResolution the default mass resolution, used if no
 * "massResolution" is given
 * @return pointer to the parsed rate measurement
 */
std::shared_ptr<SubMeasurement>
readRateMeasurement(const nlohmann::json &json, predictions::Collider collider,
                    std::shared_ptr<predictions::Particle> referenceParticle,
                    double defaultMassResolution);

/**
 * @brief Read a MassMeasurement from json. The json representation is
 * ```json
 * {
 *      "process": {},
 *      "obsMass": [124.8, 125, 125.3],
 *      "channelWeights": [1.,2.],
 * }
 * ```
 * where the process field is read into a predictions::ChannelProcess using
 * predictions::readChannelProcess. The three values in the `"obsMass"` field
 * represent
 * ```json
 * ["lower uncertainty interval", "central value", "upper uncertainty interval"]
 * ```
 * The `"channelWeights"` field is an array of the same size as
 * `/process/channels` that gives multiplicative weights for each channel. It is
 * *optional* and if it is not present all weights default to 1.
 *
 * @param json the json data to read from
 * @param collider the collider at which the "process" takes place, passed to
 * predictions::readChannelProcess
 * @param referenceParticle the reference particle at the reference mass for the
 * measurement
 * @return MassMeasurement pointer to the parsed mass measurement
 */
std::shared_ptr<SubMeasurement>
readMassMeasurement(const nlohmann::json &json, predictions::Collider collider,
                    std::shared_ptr<predictions::Particle> referenceParticle);

/**
 * @brief Read a CouplingMeasurement from json. The json representation is
 * ```json
 * {
 *      "process": {},
 *      "coupling": "effCPoTopYuk",
 *      "obsCoupling": [-0.3, 0.01, 0.4],
 *      "channelWeights": [1.,2.],
 *      "massResolution": 2.
 * }
 * ```
 * where the process field is read into a predictions::ChannelProcess using
 * predictions::readChannelProcess and used as the measurementProcess. The
 * `"coupling"` field specifies which coupling was measured. The three values in
 * the `"obsCoupling"` field represent
 * ```json
 * ["lower uncertainty interval", "observed value", "upper uncertainty
 * interval"]
 * ```
 * The `"channelWeights"` field is an array of the same size as
 * `/process/channels` that gives multiplicative weights for each channel. It is
 * *optional* and if it is not present all weights default to 1. The
 * `"massResolution"` field is *optional* as well and can be specified to
 * override the `defaultMassResolution`.
 *
 * @param json the json data to read from
 * @param collider the collider at which the "process" takes place, passed to
 * predictions::readChannelProcess
 * @param referenceParticle the reference particle at the reference mass for the
 * measurement
 * @param defaultMassResolution the default mass resolution, used if no
 * "massResolution" is given
 * @return MassMeasurement pointer to the parsed mass measurement
 */
std::shared_ptr<SubMeasurement> readCouplingMeasurement(
    const nlohmann::json &json, predictions::Collider collider,
    std::shared_ptr<predictions::Particle> referenceParticle,
    double defaultMassResolution);

std::shared_ptr<SubMeasurement>
readSubMeasurement(const nlohmann::json &json, predictions::Collider collider,
                   std::shared_ptr<predictions::Particle> referenceParticle,
                   double defaultMassResolution);

double readFromSymmNestedDicts(const nlohmann::json &json,
                               const std::string &key1, const std::string &key2,
                               double defaultValue);

void readCorrelationMatrix(const nlohmann::json &json,
                           const std::vector<std::string> &binNames,
                           Eigen::Ref<Eigen::MatrixXd> correlationMatrix);

} // namespace Higgs::signals
