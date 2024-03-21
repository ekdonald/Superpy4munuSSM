/**
 * @file RateMeasurement.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief
 *
 * @copyright Copyright 2021 by the authors.
 * This file is part of HiggsTools.
 * HiggsTools is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/Process.hpp"
#include "Higgs/signals/Measurement.hpp"
#include "signals/Uncertainties.hpp"
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace Higgs::predictions {
class Predictions;
}

namespace Higgs::signals {

class CommonMeasurement : public SubMeasurement {
  public:
    CommonMeasurement(
        predictions::ChannelProcess &&process,
        std::vector<double> &&channelWeights,
        std::shared_ptr<predictions::Particle> &&referenceParticle,
        double massResolution);

    /**
     * @brief Compute the model-predicted signal strength of the given particle
     * for this measurement.
     *
     * @param p the particle
     * @param modificationFactors modification factors for the individual
     * channels of the process. Any missing elements are assumed to be 1, any
     * extra elements are ignored.
     * @param normalization whether to use the particle mass or the reference
     * mass for the signal strength normalization
     * @return double signal strength normalized to the reference model rate
     */
    double signalStrength(const predictions::Particle &p,
                          const ModificationFactors &modificationFactors,
                          NormalizeAt normalize) const override;

    predictions::ParticleSet
    assignParticles(const predictions::Predictions &predictions, PDF massPDF,
                    double assignmentMultiplier) const;

    //! The reference mass used in this rate measurement. This is the both the
    //! signal mass assumed in the analysis, and the mass where the
    //! referenceRate() is evaluated.
    double referenceMass() const noexcept;

    //! The mass resolution used for the purpose of assignments.
    double massResolution() const noexcept;

    //! A description of the process this measurement is sensitive to.
    std::string processDesc(bool keepOrder = false) const override;

  protected:
    auto signalStrengthFor(const ModificationFactors &modFacs,
                           RescaleToRefMass rescaling) const {
        return [this, &modFacs,
                rescaling](const predictions::Particle &p) -> double {
            if (rescaling == RescaleToRefMass::always ||
                (rescaling == RescaleToRefMass::withinMassUnc &&
                 std::abs(referenceMass() - p.mass()) <= p.massUnc())) {
                return signalStrength(p, modFacs, NormalizeAt::particleMass);
            }
            return signalStrength(p, modFacs, NormalizeAt::referenceMass);
        };
    }

  private:
    std::vector<double>
    effectiveWeights(double particleMass,
                     const ModificationFactors &modificationFactors,
                     NormalizeAt normalize) const;

    predictions::ChannelProcess process_;
    std::vector<predictions::ChannelProcess> subprocesses_;
    std::vector<double> channelWeights_;
    std::shared_ptr<predictions::Particle> referenceParticle_;
    double massResolution_;
};

class RateMeasurement : public CommonMeasurement {
  public:
    RateMeasurement(predictions::ChannelProcess process,
                    std::vector<double> channelWeights,
                    const UncertainValue &observedRate,
                    std::shared_ptr<predictions::Particle> referenceParticle,
                    double massResolution, const UncertainValue &referenceRate,
                    bool massSensitive = false);

    SubMeasurementEvaluation
    evaluate(const predictions::Predictions &predictions,
             const ModificationFactors &modFacs,
             const MeasurementOptions &options) const override;

    //! The rate observed in this measurement.
    double observedRate() const noexcept;
    //! The uncertainty of the rate observed in this measurement. The argument
    //! selects whether to return the uncertainty for positive or negative
    //! deviations. For symmetric uncertainties both options are identical.
    double observedRateUnc(Uncertainty which) const noexcept;

    //! The rate of the reference model for this measurement.
    double referenceRate() const noexcept;
    //! The theoretical rate uncertainty of the reference model. The argument
    //! selects whether to return the uncertainty for positive or negative
    //! deviations. For symmetric uncertainties both options are identical.
    double referenceRateUnc(Uncertainty uncertainty) const noexcept;

    //! Is this measurement mass sensitive in the sense that its mass resolution
    //! is equal or comparable to the uncertainty of a mass measurement?
    bool massSensitive() const noexcept;

  private:
    UncertainValue obsRate_;
    UncertainValue refRate_;
    bool massSensitive_;
};

class MassMeasurement : public CommonMeasurement {
  public:
    MassMeasurement(predictions::ChannelProcess process,
                    std::vector<double> channelWeights,
                    const UncertainValue &observedMass,
                    std::shared_ptr<predictions::Particle> referenceParticle);

    SubMeasurementEvaluation
    evaluate(const predictions::Predictions &predictions,
             const ModificationFactors &modFacs,
             const MeasurementOptions &options) const override;

    double observedMass() const noexcept;

    double observedMassUnc(Uncertainty uncertainty) const noexcept;

  private:
    UncertainValue obsMass_;
};

class CouplingMeasurement : public CommonMeasurement {
  public:
    CouplingMeasurement(
        predictions::ChannelProcess measurementProcess,
        std::vector<double> channelWeights, predictions::Coupling coupling,
        const UncertainValue &observedCouplingValue,
        std::shared_ptr<predictions::Particle> referenceParticle,
        double massResolution);

    SubMeasurementEvaluation
    evaluate(const predictions::Predictions &predictions,
             const ModificationFactors &modFacs,
             const MeasurementOptions &options) const override;

    double observedCoupling() const noexcept;

    double observedCouplingUnc(Uncertainty uncertainty) const noexcept;

  private:
    predictions::Coupling coup_;
    UncertainValue obsCoup_;
};

} // namespace Higgs::signals
