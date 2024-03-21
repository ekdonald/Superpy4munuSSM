#include "signals/Measurements/RateMeasurement.hpp"
#include "Higgs/Predictions.hpp"
#include "predictions/Clustering.hpp"
#include "signals/Uncertainties.hpp"
#include "utilities/Algorithm.hpp"
#include "utilities/Logging.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <magic_enum.hpp>
#include <optional>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/unreachable_sentinel.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/repeat.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip_with.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace Higgs::signals {

namespace {
void checkWeightCount(const predictions::ChannelProcess &process,
                      const std::vector<double> &channelWeights) {
    if (process.size() != channelWeights.size()) {
        throw InvalidMeasurement(
            fmt::format("Number of weights {} does not match number of "
                        "channels {} in rate measurement",
                        channelWeights.size(), process.size()));
    }
}
} // namespace

CommonMeasurement::CommonMeasurement(
    predictions::ChannelProcess &&process, std::vector<double> &&channelWeights,
    std::shared_ptr<predictions::Particle> &&referenceParticle,
    double massResolution)
    : process_{std::move(process)}, subprocesses_{process_.subprocesses()},
      channelWeights_{std::move(channelWeights)},
      referenceParticle_{std::move(referenceParticle)}, massResolution_{
                                                            massResolution} {
    checkWeightCount(process_, channelWeights_);
}

double CommonMeasurement::signalStrength(
    const predictions::Particle &p,
    const std::vector<double> &modificationFactors,
    NormalizeAt normalize) const {
    static const auto log = logger();
    const auto effWeights =
        effectiveWeights(p.mass(), modificationFactors, normalize);
    log->trace("computing signal strength at {} with weights {}",
               magic_enum::enum_name(normalize), effWeights);
    return process_(p, effWeights) /
           process_(*referenceParticle_, channelWeights_);
}

std::vector<double> CommonMeasurement::effectiveWeights(
    double particleMass, const ModificationFactors &modificationFactors,
    NormalizeAt normalize) const {
    auto paddedModFactors =
        ranges::views::concat(modificationFactors, ranges::views::repeat(1.));
    auto modWeights = ranges::views::zip_with(
        std::multiplies<double>{}, channelWeights_, paddedModFactors);
    const auto reweight = [this, particleMass, normalize](double weight,
                                                          const auto &process) {
        if (normalize == NormalizeAt::particleMass) {
            auto refAtPMass = referenceParticle_->clone();
            refAtPMass->setMass(particleMass);
            if (process(*refAtPMass) > 0) {
                return weight * process(*referenceParticle_) /
                       process(*refAtPMass);
            }
        }
        return weight;
    };
    return ranges::views::zip_with(reweight, modWeights, subprocesses_) |
           ranges::to<std::vector>;
}

double CommonMeasurement::referenceMass() const noexcept {
    return referenceParticle_->mass();
}

double CommonMeasurement::massResolution() const noexcept {
    return massResolution_;
}

std::string CommonMeasurement::processDesc(bool keepOrder) const {
    return process_.to_string(!keepOrder);
}

predictions::ParticleSet
CommonMeasurement::assignParticles(const predictions::Predictions &predictions,
                                   PDF massPDF,
                                   double uncertaintyMultiplier) const {
    const auto deltaM = massResolution() * uncertaintyMultiplier;

    const auto massWindow = std::array{
        std::pair{referenceMass() - deltaM, referenceMass() + deltaM}};

    const auto particleSignalStrength = [this](const predictions::Particle &p) {
        return signalStrength(p, {}, NormalizeAt::referenceMass);
    };

    const double theoUncMultiplier =
        massPDF == PDF::gaussian ? uncertaintyMultiplier : 1.;

    auto [cluster] = predictions::Clustering::relevantParticles(
        particleSignalStrength, predictions, massWindow,
        predictions::MassUncEagerness::eager, theoUncMultiplier);

    if (massPDF == PDF::gaussian) {
        // relevantParticles selects particles with masses in `m +-
        // assignmentMultiplier * (theoUnc + expUnc)` but we only want those in
        // `m +- assignmentMultiplier *sqrt(theoUnc^2 + expUnc^2)`
        const auto outOfGaussianInterval =
            [this, uncertaintyMultiplier](const predictions::Particle &p) {
                return std::abs(p.mass() - referenceMass()) >
                       uncertaintyMultiplier *
                           std::sqrt(std::pow(p.massUnc(), 2) +
                                     std::pow(massResolution(), 2));
            };
        utilities::eraseFromSetIf(cluster, outOfGaussianInterval);
    }
    return cluster;
}

RateMeasurement::RateMeasurement(
    predictions::ChannelProcess process, std::vector<double> channelWeights,
    const UncertainValue &observedRate,
    std::shared_ptr<predictions::Particle> referenceParticle,
    double massResolution, const UncertainValue &referenceRate,
    bool massSensitive)
    : CommonMeasurement{std::move(process), std::move(channelWeights),
                        std::move(referenceParticle), massResolution},
      obsRate_{observedRate}, refRate_{referenceRate}, massSensitive_{
                                                           massSensitive} {}

MassMeasurement::MassMeasurement(
    predictions::ChannelProcess process, std::vector<double> channelWeights,
    const UncertainValue &observedMass,
    std::shared_ptr<predictions::Particle> referenceParticle)
    : CommonMeasurement{std::move(process), std::move(channelWeights),
                        std::move(referenceParticle),
                        observedMass.symmetricUncertainty()},
      obsMass_{observedMass} {}

CouplingMeasurement::CouplingMeasurement(
    predictions::ChannelProcess measuredProcess,
    std::vector<double> channelWeights, predictions::Coupling coupling,
    const UncertainValue &observedCouplingValue,
    std::shared_ptr<predictions::Particle> referenceParticle,
    double massResolution)
    : CommonMeasurement{std::move(measuredProcess), std::move(channelWeights),
                        std::move(referenceParticle), massResolution},
      coup_{coupling}, obsCoup_{observedCouplingValue} {}

SubMeasurementEvaluation
RateMeasurement::evaluate(const Predictions &predictions,
                          const ModificationFactors &modFacs,
                          const MeasurementOptions &options) const {
    static const auto log = logger();
    const auto cluster = assignParticles(
        predictions, options.theoryMassUncPDF,
        massSensitive() ? options.massSensitiveAssignmentRange : 1.);
    log->trace("assigned {} particles to rate measurement", cluster.size());
    const auto mu = signalStrengthFor(modFacs, options.rescaleToRefMass);
    const auto signalStrength =
        predictions::Clustering::combinedRate(mu, cluster);
    log->trace("total signal strength of {} in rate measurement",
               signalStrength);
    const auto modelRate = signalStrength * referenceRate();
    // log->info("observedRate {}", observedRate());
    // log->info("modelRate {}", modelRate);
    // log->info("signalStrength {}", signalStrength);
    // log->info("refRate_.uncertaintyTowards(modelRate) {}", refRate_.uncertaintyTowards(modelRate));

    return {observedRate() - modelRate, obsRate_.uncertaintyTowards(modelRate),
            signalStrength * refRate_.uncertaintyTowards(modelRate)};
}

namespace {
constexpr auto particleMass = [](const predictions::Particle &p) {
    return p.mass();
};
constexpr auto particleMassUnc = [](const predictions::Particle &p) {
    return p.massUnc();
};
constexpr auto zero = [](const auto &) { return 0.; };
//! Separation chisq treating all uncertainties as gaussian errors
template <class Weighting, class Observable, class TheoUnc = decltype(zero)>
double separationChisqGauss(const predictions::ParticleSet &cluster,
                            const Weighting &weightFunc,
                            const Observable &observableFunc, double avgValue,
                            double obsUnc, const TheoUnc &theoUncFunc = zero) {
    if (cluster.size() > 1) {
        const auto sepChisq = [&](const predictions::Particle &p) {
            return std::pow(observableFunc(p) - avgValue, 2) /
                   (std::pow(theoUncFunc(p), 2) + std::pow(obsUnc, 2));
        };
        return utilities::weightedMean(cluster, weightFunc, sepChisq);
    }
    return 0;
}

//! Separation chisq treating the observed uncertainty as a gaussian error and
//! the theory uncertainty as a flat box.
template <class Weighting, class Observable, class TheoUnc = decltype(zero)>
double separationChisqBoxGauss(const predictions::ParticleSet &cluster,
                               const Weighting &weightFunc,
                               const Observable &observableFunc,
                               double avgValue, double obsUnc,
                               const TheoUnc &theoUncFunc = zero) {
    if (cluster.size() > 1) {
        const auto sepChisq = [&](const predictions::Particle &p) {
            return std::pow(std::max(std::abs(observableFunc(p) - avgValue) -
                                         theoUncFunc(p),
                                     0.) /
                                obsUnc,
                            2);
        };
        return utilities::weightedMean(cluster, weightFunc, sepChisq);
    }
    return 0;
}
} // namespace

SubMeasurementEvaluation
MassMeasurement::evaluate(const Predictions &predictions,
                          const ModificationFactors &modFacs,
                          const MeasurementOptions &options) const {
    static const auto log = logger();
    const auto cluster = assignParticles(predictions, options.theoryMassUncPDF,
                                         options.massSensitiveAssignmentRange);
    if (cluster.empty()) {
        log->trace("No particles assigned to mass measurement.");
        return {0, obsMass_.symmetricUncertainty(), 0,
                options.unassignedMassMeasurementPenalty *
                    std::pow(options.massSensitiveAssignmentRange, 2)};
    }
    const auto mu = signalStrengthFor(modFacs, options.rescaleToRefMass);
    const auto [weightedMass] =
        predictions::Clustering::rateWeightedMasses(mu, cluster);
    const auto obsMassUnc = obsMass_.uncertaintyTowards(weightedMass.mass);
    log->trace("evaluating mass chisq at {}+-{} ({} particles) for "
               "measured {}+-{}.",
               weightedMass.mass, weightedMass.uncertainty, cluster.size(),
               obsMass_.central(), obsMassUnc);
    switch (options.theoryMassUncPDF) {
    case PDF::gaussian:
        return {obsMass_.central() - weightedMass.mass, obsMassUnc,
                weightedMass.uncertainty,
                separationChisqGauss(cluster, mu, particleMass,
                                     weightedMass.mass, obsMassUnc,
                                     particleMassUnc)};
    case PDF::box: {
        const auto origResidual = obsMass_.central() - weightedMass.mass;
        const auto reducedResidual =
            std::max(std::abs(origResidual) - weightedMass.uncertainty, 0.);
        return {std::copysign(reducedResidual, origResidual), obsMassUnc, 0.,
                separationChisqBoxGauss(cluster, mu, particleMass,
                                        weightedMass.mass, obsMassUnc,
                                        particleMassUnc)};
    }
    }
    throw std::runtime_error("unreachable"); // LCOV_EXCL_LINE
}

SubMeasurementEvaluation
CouplingMeasurement::evaluate(const predictions::Predictions &predictions,
                              const ModificationFactors &modFacs,
                              const MeasurementOptions &options) const {
    static const auto log = logger();

    auto cluster = assignParticles(predictions, options.theoryMassUncPDF, 1.);
    // only consider particles for which the coupling has been provided
    utilities::eraseFromSetIf(cluster, [this](const predictions::Particle &p) {
        return !p.coupling(coup_).has_value();
    });
    if (cluster.empty()) {
        log->trace("No particles assigned to coupling measurement yielding a "
                   "chisq penalty of {}",
                   options.unassignedCouplingMeasurementPenalty);
        return {0, obsCoup_.symmetricUncertainty(), 0,
                options.unassignedCouplingMeasurementPenalty};
    }
    const auto mu = signalStrengthFor(modFacs, options.rescaleToRefMass);
    const auto couplingValue = [this](const predictions::Particle &p) {
        // we know that the coupling has a value since we discarded all other
        // particles above
        return p.coupling(coup_).value();
    };
    auto weightedCoupling = utilities::weightedMean(cluster, mu, couplingValue);
    auto couplingUnc = obsCoup_.uncertaintyTowards(weightedCoupling);
    log->trace("evaluating coupling chisq at {} ({} particles) for "
               "measured {}+-{}.",
               weightedCoupling, cluster.size(), obsCoup_.central(),
               couplingUnc);
    return {obsCoup_.central() - weightedCoupling, couplingUnc, 0.,
            separationChisqGauss(cluster, mu, couplingValue, weightedCoupling,
                                 couplingUnc)};
}

double RateMeasurement::observedRate() const noexcept {
    return obsRate_.central();
}
double RateMeasurement::referenceRate() const noexcept {
    return refRate_.central();
}
double RateMeasurement::observedRateUnc(Uncertainty which) const noexcept {
    switch (which) {
    case Uncertainty::minus:
        return obsRate_.lowerUncertainty();
    default:
        return obsRate_.upperUncertainty();
    }
}
double RateMeasurement::referenceRateUnc(Uncertainty which) const noexcept {
    switch (which) {
    case Uncertainty::minus:
        return refRate_.lowerUncertainty();
    default:
        return refRate_.upperUncertainty();
    }
}
bool RateMeasurement::massSensitive() const noexcept { return massSensitive_; }

double MassMeasurement::observedMass() const noexcept {
    return obsMass_.central();
}
double MassMeasurement::observedMassUnc(Uncertainty which) const noexcept {
    switch (which) {
    case Uncertainty::minus:
        return obsMass_.lowerUncertainty();
    default:
        return obsMass_.upperUncertainty();
    }
}
double CouplingMeasurement::observedCoupling() const noexcept {
    return obsCoup_.central();
}
double
CouplingMeasurement::observedCouplingUnc(Uncertainty which) const noexcept {
    switch (which) {
    case Uncertainty::minus:
        return obsCoup_.lowerUncertainty();
    default:
        return obsCoup_.upperUncertainty();
    }
}
} // namespace Higgs::signals
