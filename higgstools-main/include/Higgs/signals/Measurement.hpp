/**
 * @file Measurement.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Interfaces to individual HiggsSignals Measurement#s and
 * SubMeasurement%s as well as definitions for all of the LimitOption#s.
 *
 * @copyright Copyright 2022 by the authors. This file is part of HiggsTools.
 * HiggsTools is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Higgs {
namespace predictions {
class Particle;
class Predictions;
} // namespace predictions

namespace signals {

//! Used for specifying the direction of a (potentially asymmetric) uncertainty
enum class Uncertainty {
    plus, //!< uncertainty for positive deviations
    minus //!< uncertainty for negative deviations
};

//! Which correlations to include
enum class Correlations {
    none,             //!< ignore all correlations
    experimentalOnly, //!< only include experimental correlations
    theoryOnly,       //!< only include theoretical correlations
    both              //!< include all available correlation information
};

/**
 * @brief Under what circumstances should a model prediction be rescaled to the
 * reference mass.
 *
 * A model prediction for a rate observables is evaluated from the signal
 * strength \f$\mu\f$ and the reference rate of the observable as
 *
 * \f[\mu \cdot \sigma^\mathrm{observable}_\mathrm{ref}(m_r) =
 * \frac{\sigma^\mathrm{incl}(m)}{\sigma^\mathrm{incl}_\mathrm{ref}(m_\mathrm{norm})}
 * \sigma^\mathrm{observable}_\mathrm{ref}(m_\mathrm{ref}) \f]
 *
 * where \f$m\f$ is the particle mass and \f$m_\mathrm{ref}\f$ is the reference
 * mass. This flag determines whether to choose \f$m_\mathrm{norm} =
 * m_\mathrm{ref}\f$ or \f$m_\mathrm{norm} = m\f$.
 *
 * Somewhat counterintuitively, choosing \f$m_\mathrm{norm} = m_\mathrm{ref}\f$,
 * i.e. normalizing the signal strength to the reference mass, means **not**
 * rescaling to rate to the reference mass. It turns the above equation into
 *
 * \f[ \sigma^\mathrm{incl}(m)
 *  \frac{\sigma^\mathrm{observable}_\mathrm{ref}(m_\mathrm{ref})}{\sigma^\mathrm{incl}_\mathrm{ref}(m_\mathrm{ref})}
 *  = \sigma^\mathrm{incl}(m)\mathcal{A}_\mathrm{ref}(m_r) \f]
 *
 * which assumes the acceptance is constant over the mass range
 * \f$[m_\mathrm{ref},m]\f$ and directly uses the model-predicted rate without
 * any rescaling.
 *
 * Choosing \f$m_\mathrm{norm} = m\f$, i.e. using the signal strength normalized
 * at the particle mass, instead leads to
 *
 * \f[ \sigma^\mathrm{incl}(m)
 * \frac{\sigma^\mathrm{incl}_\mathrm{ref}(m_\mathrm{ref})}{\sigma^\mathrm{incl}_\mathrm{ref}(m)}
 * \mathcal{A}_\mathrm{ref}(m_\mathrm{ref}) \f]
 *
 * which rescales the model-predicted inclusive rate to the reference mass by
 * assuming that the mass dependence of the model rate matches the mass
 * dependence of the reference rate over the mass range
 * \f$[m_\mathrm{ref},m]\f$.  This effectively makes HiggsSignals approximate
 * what the model-predicted rate would be if the particle mass exactly equaled
 * the reference mass. For inclusive rates that are the sums of multiple
 * channels this scaling is performed separately for each channel.
 */
enum class RescaleToRefMass {
    //! always directly use the model-predicted rates
    never,
    //! always rescale the model-predicted rates to the reference mass
    always,
    //! only rescale the model-predicted rate to the reference mass if the
    //! reference mass lies within the theoretical mass uncertainty of the
    //! particle
    withinMassUnc,
};

//! Which probability density function to use
enum class PDF {
    //! Use a flat box PDF. This means treating every value within the
    //! uncertainty as equally likely. This may be appropriate for theory
    //! uncertainties.
    box,
    //! A gaussian PDF. This means treating the uncertainties as gaussian errors
    //! around the central value. This is always used for experimental
    //! uncertainties, but may also be appropriate for theory uncertainties.
    gaussian
};

//! Options that modify the behaviour of Measurement%s and SubMeasurement%s.
struct HIGGSTOOLS_EXPORT MeasurementOptions {
    //! Which probability density function to use for the theoretical mass
    //! uncertainty. A gaussian PDF treats it like a gaussian error and combines
    //! it with the experimental mass uncertainty using error propagation, while
    //! a box PDF treats it as a flat box.
    //!
    //! This option has no impact on the experimental mass uncertainty, which is
    //! always treated as a gaussian error.
    PDF theoryMassUncPDF = PDF::box;

    //! Sets the assignment range for mass measurements and mass-sensitive rate
    //! measurements.
    //!
    //! All particles that lie in range of a mass measurement at this confidence
    //! level are assigned to the measurement. The effect of this parameter
    //! depends on the chosen theoryMassUncPDF. For PDF::box only the
    //! experimental uncertainty of the mass measurement is multiplied, while
    //! for PDF::gaussian case the (root-of-squared-sum) combined experimental
    //! and theoretical uncertainty is multiplied.
    //!
    //! This is also used for other measurements that are marked as
    //! mass-sensitive, i.e. rate measurements in high resolution channels with
    //! identical effect (treating the mass resolution as the experimental
    //! uncertainty).
    //!
    //! If no particles are assigned to a mass measurement the resulting penalty
    //! is the \f$\chi^2\f$ corresponding to this confidence level multiplied
    //! with the #unassignedMassMeasurementPenalty.
    double massSensitiveAssignmentRange = 3.;

    //! Additional multiplicative penalty to apply on the \f$\chi^2\f$ value of
    //! a mass measurement to which no particle could be assigned.
    //! @see massSensitiveAssignmentRange
    double unassignedMassMeasurementPenalty = 2.;

    //! An absolute \f$\chi^2\f$ penalty to apply for a coupling measurement for
    //! which no particle could be assigned.
    double unassignedCouplingMeasurementPenalty = 4.;

    //! Whether to directly use the model predicted rates or rescale them to the
    //! reference mass.
    //! @see RescaleToRefMass
    RescaleToRefMass rescaleToRefMass = RescaleToRefMass::always;

    //! Which types of correlations to include
    Correlations whichCorrelations = Correlations::both;

    //! whether to ignore all reference model theory uncertainties
    bool ignoreTheoryUncertainties = false;
};

//! Quantities obtained from evaluating a SubMeasurement
struct HIGGSTOOLS_EXPORT SubMeasurementEvaluation {
    double residual;        //!< the residual of the model predictions
    double obsVariance;     //!< the observed variance
    double refVariance = 0; //!< the reference model variance
    //! an explicit \f$\chi^2\f$ contribution to add
    double extraChisq = 0;
};

//! Where to normalize a signal strength
enum class NormalizeAt {
    particleMass, //!< normalize to the reference model at the particle mass
    referenceMass //!< normalize to the reference model at the Measurement's
                  //!< reference mass
};

//! A general sub-measurement.
class HIGGSTOOLS_EXPORT SubMeasurement {
  public:
    //! Represents a set of modification factors that rescale the individual
    //! channels of the signal process. They are assumed to be in the order of
    //! the channels that define the process. This order can be retrieved from
    //! processDesc().
    using ModificationFactors = std::vector<double>;

    /**
     * @brief Evaluates the signal strength of a particle in the signal process.
     *
     * @param p the particle
     * @param modificationFactors multiplicative scaling factors to apply to the
     * rates of `p` in the individual channels of the signal process
     * @param normalize at which mass should the denominator of the signal
     * strength be evaluated. @see RescaleToRefMass.
     * @return double the signal strength
     */
    virtual double
    signalStrength(const predictions::Particle &p,
                   const ModificationFactors &modificationFactors,
                   NormalizeAt normalize) const = 0;

    /**
     * @brief Evaluate this measurement for the provided model predictions.
     *
     * @param predictions the model predictions, automatically assigns particles
     * that are relevant to the signal process and lie in the mass range of the
     * measurement.
     * @param modFacs multiplicative scaling factors to apply to the
     * rates of `p` in the individual channels of the signal process
     * @param options measurement options to use
     * @return SubMeasurementEvaluation individual evaluation results
     */
    virtual SubMeasurementEvaluation
    evaluate(const predictions::Predictions &predictions,
             const ModificationFactors &modFacs,
             const MeasurementOptions &options) const = 0;

    /**
     * @brief Compute the \f$\chi^2\f$ of the model predictions for this single
     * SubMeasurement.
     *
     * This cannot include any correlations that may be included when the full
     * Measurement or even the full Higgs::signals::Signals dataset is evaluated
     * instead.
     *
     * @param predictions the model predictions, automatically assigns particles
     * that are relevant to the signal process and lie in the mass range of the
     * measurement.
     * @param modFacs multiplicative scaling factors to apply to the
     * rates of `p` in the individual channels of the signal process
     * @param options measurement options to use
     * @return double \f$\chi^2\f$
     */
    double chisq(const predictions::Predictions &predictions,
                 const ModificationFactors &modFacs,
                 const MeasurementOptions &options) const;

    /**
     * @brief A description of the underlying signal process.
     * @param keepOrder if true, the channels are returned exactly in the order
     * they are defined in, which is the order that the #ModificationFactors
     * have to match. Otherwise channels may be reordered to achieve a more
     * compact representation.
     * @return std::string a description of the signal process
     */
    virtual std::string processDesc(bool keepOrder = false) const = 0;

    //! @internal
    virtual ~SubMeasurement() noexcept = default;
};

struct MeasurementData;

/**
 * @brief A HiggsSignals measurement.
 *
 * A measurement implements the complete set of results from one experimental
 * analysis. The individual results are are stored as SubMeasurement%s and can
 * represent e.g. STXS bins, different signal strengths, or even measurements of
 * masses or couplings. This class also handles any available experimental and
 * theory correlations between these sub-measurements.
 */
class HIGGSTOOLS_EXPORT Measurement {
  public:
    //! Modification factors that can be used to rescale the model predictions
    //! in the SubMeasurement%s. The keys correspond to the names of the
    //! SubMeasurement%s. @see subMeasurements()
    using ModificationFactors =
        std::unordered_map<std::string, SubMeasurement::ModificationFactors>;

    //! Read a measurement from a datafile using the provided options.
    //! @throws InvalidMeasurement if the file is not a valid measurement
    Measurement(const std::string &filePath,
                const MeasurementOptions &options = {});

    //! Provides access to the SubMeasurement%s of this measurement. The keys of
    //! the returned map are the ids of the sub-measurements, which corespond to
    //! the names of the implemented experimental results.
    const std::unordered_map<std::string, std::shared_ptr<SubMeasurement>> &
    subMeasurements() const noexcept;

    //! The number of SubMeasurement%s contained in this measurement.
    std::size_t nSubMeasurements() const noexcept;

    //! The reference mass (in GeV) assumed in this measurement. This is the
    //! same for all of the SubMeasurement%s and is typically the assumed mass
    //! of the signal particle in the experimental analysis.
    double referenceMass() const noexcept;

    //! The reference model used in this measurement for normalization and
    //! rescaling purposes.
    predictions::ReferenceModel referenceModel() const noexcept;

    //! The options that this limit was loaded with.
    const MeasurementOptions &options() const noexcept;

    /**
     * @brief Evaluates the \f$\chi^2\f$ value of this measurement for the given
     * predictions.
     *
     * The \f$\chi^2\f$ is computed taking into account all included (and
     * enabled depending on the MeasurementOptions) correlations between the
     * SubMeasurement%s, but of course ignoring any correlations between
     * different Measurement%s that may be included when evaluating the
     * \f$\chi^2\f$ for the full Higgs::signals::Signals dataset.
     *
     * @param predictions the model predictions for which to evaluate the
     * \f$\chi^2\f$
     * @param modificationFactors modification factors to rescale the model
     * predictions. The keys identify the SubMeasurement%s by name. The
     * SubMeasurement::ModificationFactors are then passed on to the
     * corresponding SubMeasurement%s. Any modification factors that are not
     * provided default to 1 (i.e. no modification).
     * @return double the \f$\chi^2\f$ value
     */
    double
    operator()(const predictions::Predictions &predictions,
               const ModificationFactors &modificationFactors = {}) const;

    /**
     * @brief Evaluates the individual \f$\chi^2\f$ contributions of all of the
     * SubMeasurement%s.
     *
     * Since every SubMeasurement::chisq() is evaluated individually, this
     * necessarily ignores all correlations. Therefore the sum of these
     * contributions usually does not match the \f$\chi^2\f$ obtained by
     * evaluating the full measurement.
     *
     * @param predictions the model predictions for which to evaluate the
     * \f$\chi^2\f$
     * @param modificationFactors modification factors to rescale the model
     * predictions. The keys identify the SubMeasurement%s by name. The
     * SubMeasurement::ModificationFactors are then passed on to the
     * corresponding SubMeasurement%s. Any modification factors that are not
     * provided default to 1 (i.e. no modification).
     * @return std::unordered_map<std::string, double> the \f$\chi^2\f$ of each
     * SubMeasurement by name
     */
    std::unordered_map<std::string, double> chisqContributions(
        const predictions::Predictions &predictions,
        const ModificationFactors &modificationFactors = {}) const;

    // ---- Properties ----

    //! Unique numeric ID for the Limit.
    std::size_t id() const noexcept;

    //! eprint (arXiv or note) reference
    const std::string &reference() const noexcept;

    //! Inspire-HEP bibtex cite key for the limit.
    const std::string &citeKey() const noexcept;

    //! Collider the limit is from.
    predictions::Collider collider() const noexcept;

    //! Experiment the limit is from
    predictions::Experiment experiment() const noexcept;

    //! Luminosity included in the Limit
    double luminosity() const noexcept;

    //! HiggsSignals data file from which this measurement was loaded.
    const std::string &loadedFrom() const noexcept;

  private:
    friend class Signals;
    std::shared_ptr<MeasurementData> data_;
};

//! An error throwsn when an invalid measurement definition is encountered.
class HIGGSTOOLS_EXPORT InvalidMeasurement : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

} // namespace signals
} // namespace Higgs
