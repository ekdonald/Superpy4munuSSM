/**
 * @file ReferenceModels.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Reference Cxns and BRs for a SM-like Higgs boson
 *
 * @copyright Copyright 2022 by the authors.
 * This file is part of HiggsTools.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include <memory>

namespace Higgs {
namespace predictions {

/**
 * @brief A SM-like Higgs boson.
 *
 * The cross sections, BRs and the total width are from the
 * [YR4](https://arxiv.org/abs/1610.07922) unless stated otherwise. Where
 * different calculations are available, they are chosen to be suitable for
 * rescaling to BSM models over a wide mass range. In particular this means that
 * SM EW corrections as well as higher order corrections in the heavy top limit
 * are not included.
 */
class HIGGSTOOLS_EXPORT SMHiggs : public Particle {
  public:
    //! Construct a SM-like Higgs boson with the given mass (in GeV).
    SMHiggs(double mass) noexcept;

    /**
     * @brief Return specified production cross section of this particle.
     *
     * If the mass() is larger than the maximal available tabulated value this
     * returns 0. If mass() is lower than the minimal value this returns the cxn
     * at the lowest available mass, which results in a lower estimate of the
     * correct cross section (and is therfore still useful for HiggsBounds).
     *
     * The VH cxns do not use the YR4 numbers, but use the cxn parametrization
     * from #Higgs::predictions::EffectiveCouplingCxns with SM-like effective
     * couplings instead. The inclusive ZH and WH cxns obtained from this lie in
     * the uncertainty interval of the YR4 numbers, but the implementation
     * additionally provides the sub-channels, which are not separated in the
     * YR4.
     *
     * Since the LEP cxns are defined to be SM-normalized they are always 1 by
     * definition.
     *
     * @param coll the collider for which to return the cross section
     * @param p the process of interest
     * @return double a cross section value in pb (or dimensionless for
     * normalized LEP cxns)
     */
    double cxn(Collider coll, Production p) const noexcept override;

    /**
     * @brief The specified branching ratio of this particle.
     *
     * If the mass() is outsize the tabulated mass range, returns 0.
     *
     * The mass range of the YR4 numbers has been slightly extended down to
     * 10GeV and up to 1TeV using HDECAY. Additionally, HDECAY was used to add
     * the decay into strange quarks.
     *
     * @param d the decay mode of interest
     * @return double a BR
     */
    double br(Decay d) const noexcept override;

    /**
     * @brief The total width of this particle.
     *
     * If the mass() is outside the tabulated mass range, returns the total
     * width value at the closest available mass value. That way this is never
     * 0, which would lead to ill-defined branching ratios.
     *
     * @return double the total width in GeV
     */
    double totalWidth() const noexcept override;

    //! Couplings of this particle.
    std::optional<double> coupling(Coupling c) const noexcept override;

    //! @private
    std::unique_ptr<Particle> clone() const override;

  protected:
    //! A constructor that also changes the id. Useful for subclassing.
    SMHiggs(double mass, std::string id) noexcept;
};

/**
 * @brief The SM Higgs boson at the highest available precision.
 *
 * The cross sections, BRs and the total width are from the
 * [YR4](https://arxiv.org/abs/1610.07922) unless stated otherwise and are
 * chosen to yield the highest available precision including SM EW corrections.
 *
 * Since most of these high precision calculations are only available for a
 * limited mass range around 125GeV, the implementation falls back to the
 * #SMHiggs values if no more precise values are available.
 */
class HIGGSTOOLS_EXPORT SMHiggsEW : public SMHiggs {
  public:
    //! Construct a highest precision SM Higgs boson with the given mass (in
    //! GeV).
    SMHiggsEW(double mass) noexcept;

    /**
     * @brief The specified production cross section of a SM Higgs boson at the
     * hightest precision.
     *
     * If the mass() larger than the maximal available tabulated value this
     * returns 0. If mass() is lower than the minimal value returns the cxn at
     * the lowest available mass.
     *
     * Falls back to the SMHiggs values if no more precise values are available
     * for the given mass/channel.
     *
     * Since STXS measurements can be sensitive to the ZH sub-channels, we do
     * not use the YR4 numbers for ZH, but continue to use the cxn
     * parametrization from #Higgs::predictions::EffectiveCouplingCxns with
     * SM-like effective couplings instead. This means that EW corrections are
     * not included for ZH production and its sub-channels.
     *
     * Since the LEP cxns are defined to be SM-normalized they are always 1 by
     * definition.
     *
     * @param coll the collider for which to return the cross section
     * @param p the process of interest
     * @return double a cross section value in pb (or dimensionless for
     * normalized LEP cxns)
     */
    double cxn(Collider coll, Production p) const noexcept override;

    /**
     * @brief The specified branching ratio of a SM Higgs boson at the highest
     * precision.
     *
     * If the mass() is outsize the tabulated mass range, returns 0. Falls back
     * to the #SMHiggs values if no more precise values are available for the
     * given mass/channel.
     *
     * @param d the decay mode of interest
     * @return double a BR
     */
    double br(Decay d) const noexcept override;

    /**
     * @brief The total width of a SM Higgs boson at the hightest precision.
     *
     * If the mass() is outside the tabulated mass range, returns the total
     * width value at the closest available mass value. That way this is never
     * 0, which would lead to ill-defined branching ratios. Falls back to the
     * #SMHiggs values if no more precise values are available for the given
     * mass.
     *
     * @return double the total width in GeV
     */
    double totalWidth() const noexcept override;

    //! @private
    std::unique_ptr<Particle> clone() const override;
};

/**
 * @brief Interpolation between the SMHiggs and SMHiggsEW reference
 * models. For masses below ~ 150 GeV the SMHiggsEW values are 
 * used (N3LO in the inifinit top-quark mass limit); 
 * for values above ~ 150 GeV, the SMHiggs values
 * (NNLO including finite top-quark mass effects).
 */
class HIGGSTOOLS_EXPORT SMHiggsInterp : public SMHiggs {
  public:
    //! Construct a highest precision SM Higgs boson with the given mass (in
    //! GeV).
    SMHiggsInterp(double mass) noexcept;

    /**
     * @brief The specified production cross section of a SM Higgs boson.
     *
     * If the mass() larger than the maximal available tabulated value this
     * returns 0. If mass() is lower than the minimal value returns the cxn at
     * the lowest available mass.
     *
     * Falls back to the SMHiggs values if no more precise values are available
     * for the given mass/channel.
     *
     * Since STXS measurements can be sensitive to the ZH sub-channels, we do
     * not use the YR4 numbers for ZH, but continue to use the cxn
     * parametrization from #Higgs::predictions::EffectiveCouplingCxns with
     * SM-like effective couplings instead. This means that EW corrections are
     * not included for ZH production and its sub-channels.
     *
     * Since the LEP cxns are defined to be SM-normalized they are always 1 by
     * definition.
     *
     * @param coll the collider for which to return the cross section
     * @param p the process of interest
     * @return double a cross section value in pb (or dimensionless for
     * normalized LEP cxns)
     */
    double cxn(Collider coll, Production p) const noexcept override;

    /**
     * @brief The specified branching ratio of a SM Higgs boson.
     *
     * If the mass() is outsize the tabulated mass range, returns 0. Falls back
     * to the #SMHiggs values if no more precise values are available for the
     * given mass/channel.
     *
     * @param d the decay mode of interest
     * @return double a BR
     */
    double br(Decay d) const noexcept override;

    /**
     * @brief The total width of a SM Higgs boson..
     *
     * If the mass() is outside the tabulated mass range, returns the total
     * width value at the closest available mass value. That way this is never
     * 0, which would lead to ill-defined branching ratios. Falls back to the
     * #SMHiggs values if no more precise values are available for the given
     * mass.
     *
     * @return double the total width in GeV
     */
    double totalWidth() const noexcept override;

    //! @private
    std::unique_ptr<Particle> clone() const override;
};

//! Which reference model to use
enum class ReferenceModel {
    SMHiggs,   //!< #Higgs::predictions::SMHiggs
    SMHiggsEW, //!< #Higgs::predictions::SMHiggsEW
    SMHiggsInterp, //!< #Higgs::predictions::SMHiggsInterp
};

//! Get a reference particle in the specified model at the given mass.
std::unique_ptr<Particle> getReference(ReferenceModel model, double mass);

} // namespace predictions
} // namespace Higgs
