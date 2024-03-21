#pragma once
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "Higgs/signals/Measurement.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "utilities/JsonFwd.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

namespace Higgs::predictions {
class Particle;
}
namespace Higgs::signals {
struct MeasurementData {
    MeasurementData(const nlohmann::json &j, std::string loadedFrom,
                    const MeasurementOptions options);

    std::size_t id_;
    std::string reference_;
    std::string citeKey_;
    predictions::Collider collider_;
    predictions::Experiment experiment_;
    double luminosity_;
    std::string loadedFrom_;

    double refMass_;
    predictions::ReferenceModel referenceModel_;
    std::shared_ptr<predictions::Particle> referenceParticle_;
    double massRes_;
    std::unordered_map<std::string, std::shared_ptr<SubMeasurement>>
        subMeasurements_;
    Eigen::MatrixXd corrMatExp_;
    Eigen::MatrixXd corrMatTheo_;
    MeasurementOptions options_;
};
} // namespace Higgs::signals
