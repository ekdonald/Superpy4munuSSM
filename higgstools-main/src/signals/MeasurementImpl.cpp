#include "signals/MeasurementImpl.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "signals/JsonSupport.hpp"
#include "utilities/Json.hpp"
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>
#include <vector>

namespace Higgs::signals {

MeasurementData::MeasurementData(const nlohmann::json &j,
                                 std::string loadedFrom,
                                 const MeasurementOptions options)
    : id_{utilities::readAs<std::size_t>(j, "id")},
      reference_{utilities::readAs<std::string>(j, "reference")},
      citeKey_{utilities::readAs<std::string>(j, "citeKey")},
      collider_{utilities::readAs<predictions::Collider>(j, "collider")},
      experiment_{utilities::readAs<predictions::Experiment>(j, "experiment")},
      luminosity_{utilities::readAs<double>(j, "luminosity")},
      loadedFrom_{loadedFrom}, refMass_{utilities::readAs<double>(
                                   j, "referenceMass")},
      referenceModel_{
          utilities::readAs<predictions::ReferenceModel>(j, "referenceModel")},
      referenceParticle_{getReference(referenceModel_, refMass_)},
      massRes_{utilities::readAs<double>(j, "massResolution")},
      subMeasurements_{}, corrMatExp_{Eigen::MatrixXd::Identity(
                              j.at("subMeasurements").size(),
                              j.at("subMeasurements").size())},
      corrMatTheo_{Eigen::MatrixXd::Identity(j.at("subMeasurements").size(),
                                             j.at("subMeasurements").size())},
      options_{options} {
    for (const auto &[key, jval] : j.at("subMeasurements").items()) {
        subMeasurements_.emplace(
            key,
            readSubMeasurement(jval, collider_, referenceParticle_, massRes_));
    }
    if (j.contains("correlations")) {
        const auto binNames =
            subMeasurements_ | ranges::views::keys | ranges::to<std::vector>;
        if ((options_.whichCorrelations == Correlations::both ||
             options_.whichCorrelations == Correlations::experimentalOnly) &&
            j["correlations"].contains("experimental")) {
            readCorrelationMatrix(j["correlations"]["experimental"], binNames,
                                  corrMatExp_);
        }
        if (options_.ignoreTheoryUncertainties) {
            corrMatTheo_.setZero();
        } else if ((options_.whichCorrelations == Correlations::both ||
                    options_.whichCorrelations == Correlations::theoryOnly) &&
                   j["correlations"].contains("theoretical")) {
            readCorrelationMatrix(j["correlations"]["theoretical"], binNames,
                                  corrMatTheo_);
        }
    }
}

} // namespace Higgs::signals
