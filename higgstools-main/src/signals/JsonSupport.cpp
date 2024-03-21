#include "signals/JsonSupport.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Process.hpp"
#include "Higgs/signals/Measurement.hpp"
#include "predictions/JsonSupport.hpp"
#include "signals/Measurements/RateMeasurement.hpp"
#include "signals/Uncertainties.hpp"
#include "utilities/Format.hpp"
#include "utilities/Json.hpp"
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/unreachable_sentinel.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/zip_with.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Higgs::predictions {
class Particle;
}

namespace Higgs::signals {

void from_json(const nlohmann::json &j, UncertainValue &v) {
    v = UncertainValue::fromInterval(
        j.at(0).get<double>(), j.at(1).get<double>(), j.at(2).get<double>());
}

void to_json(nlohmann::json &j, const UncertainValue &v) {
    j = nlohmann::json::array({v.central() - v.lowerUncertainty(), v.central(),
                               v.central() + v.upperUncertainty()});
}

std::shared_ptr<SubMeasurement>
readRateMeasurement(const nlohmann::json &json, predictions::Collider collider,
                    std::shared_ptr<predictions::Particle> referenceParticle,
                    double defaultMassResolution) {
    auto proc = predictions::readChannelProcess(json.at("process"), collider);
    auto channelWeights = utilities::readIfPresent<std::vector<double>>(
        json, "channelWeights", proc.size(), 1.);
    auto refParticleRate = proc(*referenceParticle, channelWeights);
    return std::make_shared<RateMeasurement>(
        std::move(proc), std::move(channelWeights),
        utilities::readAs<UncertainValue>(json, "obs"),
        std::move(referenceParticle),
        json.value("massResolution", defaultMassResolution),
        json.value("ref", UncertainValue{refParticleRate}),
        json.value("massSensitive", false));
}

std::shared_ptr<SubMeasurement>
readMassMeasurement(const nlohmann::json &json, predictions::Collider collider,
                    std::shared_ptr<predictions::Particle> referenceParticle) {
    auto proc = predictions::readChannelProcess(json.at("process"), collider);
    auto channelWeights = utilities::readIfPresent<std::vector<double>>(
        json, "channelWeights", proc.size(), 1.);
    return std::make_shared<MassMeasurement>(
        std::move(proc), std::move(channelWeights),
        utilities::readAs<UncertainValue>(json, "obsMass"),
        std::move(referenceParticle));
}

std::shared_ptr<SubMeasurement> readCouplingMeasurement(
    const nlohmann::json &json, predictions::Collider collider,
    std::shared_ptr<predictions::Particle> referenceParticle,
    double defaultMassResolution) {
    auto proc = predictions::readChannelProcess(json.at("process"), collider);
    auto channelWeights = utilities::readIfPresent<std::vector<double>>(
        json, "channelWeights", proc.size(), 1.);
    return std::make_shared<CouplingMeasurement>(
        std::move(proc), std::move(channelWeights),
        utilities::readAs<predictions::Coupling>(json, "coupling"),
        utilities::readAs<UncertainValue>(json, "obsCoupling"),
        std::move(referenceParticle),
        json.value("massResolution", defaultMassResolution));
}

std::shared_ptr<SubMeasurement>
readSubMeasurement(const nlohmann::json &json, predictions::Collider collider,
                   std::shared_ptr<predictions::Particle> referenceParticle,
                   double defaultMassResolution) {
    if (json.contains("obs")) {
        return readRateMeasurement(json, collider, referenceParticle,
                                   defaultMassResolution);
    } else if (json.contains("obsMass")) {
        return readMassMeasurement(json, collider, referenceParticle);
    } else if (json.contains("coupling")) {
        return readCouplingMeasurement(json, collider, referenceParticle,
                                       defaultMassResolution);
    }
    throw InvalidMeasurement(
        fmt::format("Invalid sub measurement: {}", json.dump()));
}

double readFromSymmNestedDicts(const nlohmann::json &json,
                               const std::string &key1, const std::string &key2,
                               double defaultValue) {
    if (json.contains(key1) && json[key1].contains(key2)) {
        return json[key1][key2].get<double>();
    } else if (json.contains(key2) && json[key2].contains(key1)) {
        return json[key2][key1].get<double>();
    }
    return defaultValue;
}

void readCorrelationMatrix(const nlohmann::json &json,
                           const std::vector<std::string> &binOrder,
                           Eigen::Ref<Eigen::MatrixXd> correlationMatrix) {
    for (const auto &[i, key1] : ranges::views::enumerate(binOrder)) {
        for (const auto &[j, key2] : ranges::views::enumerate(binOrder)) {
            if (i != j) {
                correlationMatrix(i, j) =
                    readFromSymmNestedDicts(json, key1, key2, 0.);
            }
        }
    }
}

} // namespace Higgs::signals
