/**
 * @file Limit.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Common interface for all HiggsBounds limits
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/predictions/Basics.hpp"
#include <memory>
#include <string>
#include <vector>

namespace Higgs {

namespace predictions {
class Predictions;
enum class Collider;
enum class Experiment;
} // namespace predictions

namespace bounds {
class AppliedLimit;

//! Options that govern the behaviour of Limit::apply().
struct HIGGSTOOLS_EXPORT LimitOptions {
    //! Fraction of the mass resolution used to determine the range of
    //! applicability for the limit. For example, for a mass resolution of ±10
    //! GeV a limit that starts at 100 GeV is applied to particles with masses
    //! ≥95 GeV.
    double applicableResolutionFac = 0.5;

    //! How to cluster particles with mass uncertainties.
    //!   value  | effect
    //! ---------|---------
    //! eager    | cluster as soon as uncertainties + resolution touch
    //! cautious | cluster only if the mass uncertainties overlap entirely
    //! within the resolution ignore   | ignore the mass uncertainties when
    //! clustering
    predictions::MassUncEagerness clusterMassUnc =
        predictions::MassUncEagerness::cautious;

    //! When to apply a limit to particles with mass uncertainties.
    //!   value  | effect
    //! ---------|---------
    //! eager    | if any mass within the uncertainty is within range
    //! cautious | only if the entire uncertainty interval is within range
    //! ignore   | ignore the mass uncertainty
    predictions::MassUncEagerness applicableMassUnc =
        predictions::MassUncEagerness::cautious;

    //! Which mass value of particles (or clusters) with mass uncertainties to
    //! use for the limit extraction. This is currently only exact if the mass
    //! uncertainty is smaller than the mass grid spacing, but still works as an
    //! approximation otherwise.
    //!   value  | effect
    //! ---------|---------
    //! eager    | use the mass that provides the strongest observed limit
    //! cautious | use the mass that provides the weakest observed limit
    //! ignore   | use the central mass value
    predictions::MassUncEagerness setLimitMassUnc =
        predictions::MassUncEagerness::cautious;

    //! Any applied limits with an AppliedLimit::expRatio() <= this value are
    //! considered insensitive and will not be returned.
    double minExpRatio = 1e-4;
};

//! A HiggsBounds limit. All types of limits are accessed through this
//! interface. Provides access to limit metadata and allows applying a limit to
//! a #Higgs::predictions::Predictions.
class HIGGSTOOLS_EXPORT Limit : public std::enable_shared_from_this<Limit> {
  public:
    /** @{ @name Constructor */
    //! Reads a limit from the `filePath`, passes it the `options` and returns
    //! it. This acts as a constructor and is the only way to create Limit%s.
    static std::shared_ptr<Limit> read(const std::string &filePath,
                                       const LimitOptions &options = {});
    //! @}

    /** @{ @name Evaluation */
    /**
     * Apply this limit to the given `predictions`. This function performs the
     * actual work within HiggsBounds. Very roughly it
     *
     * 1. clusters the predictions.particles() according to the mass
     *    resolution of the limit,
     * 2. evaluates the model prediction for each of those clusters (or each
     *    combinations of clusters for more complicated processes),
     * 3. compares those model predictions to the observed and expected limit,
     * 4. returns the results.
     *
     * The reason that the result can be multiple AppliedLimit%s is that one
     * limit may be sensitive to several distinct particles/clusters from the
     * `predictions`, e.g. if it covers a mass range that contains several
     * particles.
     */
    virtual std::vector<AppliedLimit>
    apply(const Higgs::predictions::Predictions &predictions) const = 0;
    //! @}

    //! @{ @name Metadata
    //! Functions that return metadata of the limit and provide descriptions of
    //! its properties.

    //! A unique numeric ID for the Limit.
    virtual unsigned id() const noexcept = 0;

    //! A uniquely identifying reference to the publication this limit is from.
    //! Ideally and arXiv eprint number, otherwise a collaboration preprint
    //! number.
    virtual const std::string &reference() const noexcept = 0;

    //! A description of the underlying process. Specifies involved production
    //! and decay modes as well as any BSM -> BSM decays. The precise format
    //! depends on the type of process.
    virtual std::string processDesc() const noexcept = 0;

    //! A description of the limit extent. Specifies the covered ranges for all
    //! relevant parameters that the limit depends on. This typically involves
    //! one or more masses and potentially decay width(s).
    virtual std::string extentDesc() const noexcept = 0;

    //! The Inspire-HEP bibtex cite key for the publication underlying the
    //! limit.
    virtual const std::string &citeKey() const noexcept = 0;

    //! The collider and center of mass energy at which this limit was obtained.
    //! If data from different energies was combined the highest is given.
    virtual predictions::Collider collider() const noexcept = 0;

    //! The experiment/experimental collaboration that obtained the limit.
    virtual predictions::Experiment experiment() const noexcept = 0;

    //! The amount of data used in the limit as integrated luminosity in
    //! \f$\mathrm{fb}^{-1}\f$. If data from different collider runs was
    //! combined this is the sum.
    virtual double luminosity() const noexcept = 0;

    //! Path to the HiggsBounds limit data file from which this limit was
    //! loaded.
    virtual const std::string &loadedFrom() const noexcept = 0;

    //! Access to the option values of this limit
    virtual const LimitOptions &options() const noexcept = 0;

    //! A full string description of the limit containing most of the
    //! information of the previous functions.
    std::string to_string() const noexcept;
    //! @}

    //! Virtual destructor of abstract base class.
    virtual ~Limit() = default;
};

//! The result of applying a Limit to a Prediction.
class HIGGSTOOLS_EXPORT AppliedLimit {
  public:
    //! Access the underlying limit (e.g. for limit metadata).
    std::shared_ptr<const Limit> limit() const noexcept;
    //! The model-predicted rate divided by the observed 95% C.L. limit.
    double obsRatio() const noexcept;
    //! The model-predicted rate divided by the expected 95% C.L. limit
    double expRatio() const noexcept;
    //! The observed likelihood value of the model prediction. If the limit
    //! includes no likelihood information this is zero.
    double obsLikelihood() const noexcept;
    //! The expected likelihood value of the model prediction. If the limit
    //! includes no likelihood information this is zero.
    double expLikelihood() const noexcept;
    //! The #Higgs::predictions::Particle::id%s of the contributing particles.
    //! The precise format depends on the topology of the process that the limit
    //! is set on. Single entries of `"+"` and `">"` are used to separate the
    //! contributing particles for different roles. E.g.
    //! `{"h4",">","h1","h2","+","h3"}` could be the contributing particles for
    //! a \f$H_i \to H_j H_k\f$ process, where `"h1"` and `"h2"` are close in
    //! mass and are clustered for the role of \f$H_j\f$. See the
    //! @verbatim embed:rst:inline :ref:`Processes` @endverbatim section for details.
    const std::vector<std::string> &contributingParticles() const noexcept;

    //! Construct an applied Limit
    AppliedLimit(std::shared_ptr<const Limit> limit, double obsRatio,
                 double expRatio,
                 std::vector<std::string> contributingParticles,
                 double obsLikelihood = 0, double expLikelihood = 0);
    //! default constructor
    AppliedLimit() = default;

  private:
    std::shared_ptr<const Limit> limit_;
    double obsRatio_;
    double expRatio_;
    std::vector<std::string> contributingParticles_;
    double obsLikelihood_;
    double expLikelihood_;
};

} // namespace bounds
} // namespace Higgs
