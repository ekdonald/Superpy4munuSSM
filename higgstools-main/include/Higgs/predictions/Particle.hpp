/**
 * @file Particle.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Different kinds of particles
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <string>

namespace Higgs {
namespace predictions {
class Predictions;
class ParticleData;
enum class ReferenceModel;

//! Interface to different kinds of particles.
class HIGGSTOOLS_EXPORT Particle {
  public:
    /** @{ @name IDs, quantum numbers, and masses
     * @brief Basic properties of the particle including quantum numbers and its
     * mass.
     */
    //! A unique id of the particle
    const std::string &id() const noexcept;

    //! The CP quantum number of the particle
    CP cp() const noexcept;
    //! The electric charge of the particle
    ECharge charge() const noexcept;

    //! The mass of the particle in GeV
    double mass() const noexcept;
    //! Set the mass of the particle to the value in GeV
    void setMass(double value) noexcept;
    //! The theoretical mass uncertainty of the particle in GeV
    double massUnc() const noexcept;
    //! Set the theoretical mass uncertainty of the particle in GeV
    void setMassUnc(double value) noexcept;
    //! @}

    /** @{ @name Production and decay rates
     * @brief These functions provide access to cxns, branching ratios, and
     * combined channel rates for different production and decay modes.
     */
    /**
     * @brief A production cross section of this particle.
     *
     * Returns 0 if the production mode does not exist  for any reason (would
     * violate charge, does not make sense for the given collider type, or
     * simply not available/set by the user).
     *
     * @param coll which collider
     * @param p which production mode
     * @return double a cross section. For LHC processes this is an absolute cxn
     * in pb (except for Production::brtHpb). For LEP processes the cxns are
     * always normalized to an appropriate reference model to account for the
     * many different combined collision energies (see the description of the
     * corresponding production modes for details).
     */
    virtual double cxn(Collider coll, Production p) const noexcept;

    /**
     * @brief A branching ratio of this particle into a SM final state.
     *
     * Returns 0 if the decay mode does not exist.
     *
     * @param d which decay mode
     * @return double the branching ratio
     */
    virtual double br(Decay d) const noexcept;

    /**
     * @brief A branching ratio of this particle into a SM + BSM final state.
     *
     * Returns 0 if the decay mode does not exist.
     *
     * @param d which type of chain decay, specifies the SM final state particle
     * @param particleId the BsmParticle.id() of the BSM final state particle
     * @return double the branching ratio
     */
    virtual double br(ChainDecay d,
                      const std::string &particleId) const noexcept;

    /**
     * @brief A branching ratio of this particle into a final state of 2 other
     * BSM particles.
     *
     * @param particleId1 the BsmParticle.id() of the first daughter particle
     * @param particleId2 the BsmParticle.id() of the second daughter particle.
     * @return double the branching ratio
     */
    virtual double br(const std::string &particleId1,
                      const std::string &particleId2) const noexcept;

    /**
     * @brief Return the channelRate for the given channel.
     *
     * This defaults to the narrow width approximation, i.e. `cxn * br` for the
     * specified channel. However, overriding channelRateOr() in a derived
     * Particle allows setting dedicated rates for any channel and returning
     * that one instead of the narrow width approximation result. Channel rates
     * can be explicitely set to e.g. account for interference effects in those
     * channels where they are important (or available) while defaulting to the
     * narrow width approximation everywhere else.
     *
     * @param coll the collider fo the channel
     * @param prod the production mode of the channel
     * @param decay the decay mode of the channel
     * @return double an explicitely set channel rate (as defined in a derived
     * type) or `Particle.cxn(Channel.coll(), Channel.prod()) *
     * Particle.br(Channel.decay())`
     */
    double channelRate(Collider coll, Production prod,
                       Decay decay) const noexcept;
    //! @}

    /** @{ @name Couplings
     * @brief This function provides access to predefined Coupling values of the
     * particle.
     */
    /**
     * @brief Return a coupling value associated to this particle.
     *
     * @param c which coupling
     * @return std::optional<double> value of the coupling. std::nullopt unless
     * the coupling has been set explicitely. If a value is set the unit depends
     * on the coupling that is requested.
     */
    virtual std::optional<double> coupling(Coupling c) const noexcept;
    //! @}

    /** @{ @name Total width
     * @brief And finally, this provides access to the total width.
     */
    //! The total width of a particle in GeV.
    virtual double totalWidth() const noexcept;
    //! @}

    //! @privatesection
    //! polymorphic clone a particle
    virtual std::unique_ptr<Particle> clone() const = 0;

    //! equality comparison of the id and basic particle properties
    bool operator==(const Particle &other) const noexcept;

    //! default destructor
    virtual ~Particle() = default;

  protected:
    /**
     * @brief Implements channelRate functionality for derived Particle%s.
     *
     * By default this simply returns the default value cxnTimesBr, which is set
     * by the channelRate() implementation. Can be overridden in derived classes
     * to implement dedicated channel rates.
     *
     */
    virtual double channelRateOr(Collider coll, Production prod, Decay decay,
                                 double cxnTimesBr) const noexcept;

    //! constructor that sets the basic resonance properties
    Particle(std::string &&id, CP cp, ECharge charge);

    //! @{
    //! copy/move operations that prevent slicing
    Particle &operator=(const Particle &) = default;
    Particle &operator=(Particle &&) = default;
    Particle(const Particle &) = default;
    Particle(Particle &&) = default;
    //! @}

  private:
    std::string id_;
    double mass_;
    double massUnc_;
    CP cp_;
    ECharge charge_;
};

/**
 * @brief A generic Bsm Particle. This class stores the bulk of the user input,
 * including particle massses and quantum numbers, production cross sections,
 * branching ratios, and channel rates (see channelRates()).
 */
class HIGGSTOOLS_EXPORT BsmParticle : public Particle {
  public:
    /** @{ @name Construction */
    //! Constructor that sets the basic properties. The properties can be
    //! accessed through the #Particle interface.
    BsmParticle(std::string id, ECharge charge, CP cp);
    //! @}

    /** @{ @name Setting production cross sections
     * @brief These functions are used to input the model predictions for the
     * particle production cross sections at different colliders. Cross sections
     * set like this can be retrived through the Particle::cxn interface and are
     * also used to construct Particle::channelRate%s if no explicit channel
     * rate is set using #setChannelRate().
     *
     * *For neutral Higgs-like scalars it is highly recommended to use the
     * #effectiveCouplingInput first and then adjust any cross sections for
     * which you have dedicated calculations at a higher precision available.*
     */
    /**
     * Set the production cross section at a collider.
     *
     * Cross sections set like this can be retrived through the Particle::cxn
     * interface and are also used to construct Particle::channelRate%s if no
     * explicit channel rate is set using #setChannelRate().
     *
     * @throws InvalidInput when trying to set a combined production mode.
     * @throws InvalidInput when trying to set LEP cxns.
     *
     * @param coll the collider and CMS energy, LEP production cross sections
     * are always normalized and have to be set through setNormalizedCxn().
     * @param p the production channel
     * @param value the new cross section value in pb
     */
    void setCxn(Collider coll, Production p, double value);

    /**
     * Set the production cross section at a collider through the normalized
     * value to a reference model.
     *
     * It is the users responsibility to check that the requested reference cxn
     * is available, otherwise it will default to zero.
     *
     *
     * @param coll the collider and CMS energy
     * @param p the production channel
     * @param value cross section normalized to corresponding reference cross
     * section.
     * @param reference which reference model to use.
     * For `coll=Collider::LEP` this argument is ignored. All single particle
     * LEP cross sections for a neutral scalars are always assumed to be
     * normalized to the #Higgs::predictions::SMHiggs, while for the neutral and
     * charged pair productions processes at LEP the following reference cross
     * sections are used:
     *  |     charge()     | reference process
     *  |------------------|-------------------
     *  | ECharge::neutral | \f$ e^+ e^- \to Z \to A H \f$
     *  | ECharge::single  | \f$ e^+ e^- \to Z \to H^+ H^- \f$
     * where \f$ H,A,H^\pm \f$ are the non-SM-like CP-even, CP-odd, and charged
     * scalars of a 2HDM where the other CP-even scalar \f$ h \f$ is SM-like
     * (exact alignment limit).
     */
    void setNormalizedCxn(Collider coll, Production p, double value,
                          ReferenceModel reference);

    //! @}

    /** @{ @name Specify branching ratios directly
     * @brief These function can be used to directly set the branching ratio of
     * the particle into different kinds of final states.
     *
     * Branching ratios set like this can be retrived through the Particle::br
     * interface and are also used to construct Particle::channelRate%s if no
     * explicit channel rate is set using #setChannelRate().
     *
     * The #setBr functions do **not** adjust the total width or the other BRs
     * in any way. They can only be used after a non-zero total width has been
     * set using #setTotalWidth. If the sum of branching ratios would exceed 1,
     * they throw an exception.
     *
     * Setting the BRs of a particle through this method would look something
     * like this:
     *
     * @rststar
     * .. literalinclude:: examples/settingBRs.cpp
     *    :language: c++
     * @endrststar
     *
     * or equivalently in python:
     *
     * @rststar
     * .. literalinclude:: examples/settingBRs.py
     *     :language: python
     * @endrststar
     *
     * Using these functions is only recommended if you have computed all
     * branching ratios and the total width using some external decay calculator
     * and just want to transfer them as they are into HiggsPredictions. In most
     * other cases it is easier to set the partial decay widths instead (using
     * #setDecayWidth) which automatically adjusts the total width. Mixed use of
     * the two methods is perfectly fine as well.
     *
     * *For supported particles it is highly recommended to use the
     * #effectiveCouplingInput first and then use the #setDecayWidth functions
     * to adjust or add decay modes.*
     */
    /**
     * Set the total width to the value in GeV. If you want to use the setBr
     * functions you need to set the total width of the particle beforehand.
     * This is not necessary for the setDecayWidth method.
     *
     * Setting the total width to zero will set all BRs to zero. This is the
     * recommended method of "starting over" with clean.
     *
     * @throws InvalidInput for negative values
     */
    void setTotalWidth(double value);

    /**
     * Set the specified branching ratio to the given value.
     * @throws InvalidInput if setting the BR would result in sum(BRs)>1
     */
    void setBr(Decay d, double value);

    /**
     * Set the branching ratio for the specified ChainDecay.
     * @throws InvalidInput if setting the BR would result in sum(BRs)>1
     * @param d the kind of decay chain
     * @param particleId the BSM daughter particle
     * @param value the branching ratio
     */
    void setBr(ChainDecay d, const std::string &particleId, double value);

    /**
     * Set the branching ratio for the decay into two BSM particles.
     * @throws InvalidInput if setting the BR would result in sum(BRs)>1
     * @param particleId1 id of the first daughter particle
     * @param particleId2 id of the second daughter particle, the order of these
     * two arguments is irrelevant
     * @param value the branching ratio
     */
    void setBr(const std::string &particleId1, const std::string &particleId2,
               double value);
    //! @}

    /** @{ @name Setting branching ratios through partial decay widths
     * @brief These functions are used to specify the partial decay widths of
     * the particle into different final states. Setting a partial decay width
     * like this will
     *
     * 1. adjust the total width (taking into account the previous value of the
     *    decay width into the modified final state) thereby also adjusting the
     *    BRs into all other final states
     * 2. set the BR of the modified final state according to the specified
     *    partial decay width and the adjusted total width
     *
     * Branching ratios set like this can be retrived through the Particle::br
     * interface and are also used to construct Particle::channelRate%s if no
     * explicit channel rate is set using #setChannelRate().
     *
     * The following examples demonstrate the use of these functions:
     *
     * @rststar
     * .. literalinclude:: examples/settingDecayWidths.cpp
     *     :language: c++
     * @endrststar
     *
     * or equivalently in python:
     *
     * @rststar
     * .. literalinclude:: examples/settingDecayWidths.py
     *     :language: python
     * @endrststar
     *
     * The examples set a previous total width and BR using #setBr to show how
     * the two methods can be mixed, but this is not necessary. These functions
     * work just fine starting from a particle with `totalWidth()==0`. On the
     * other hand, they work equally well if you have a complicated predefined
     * set of BRs, i.e. from the #effectiveCouplingInput.
     *
     * *For supported particles it is highly recommended to use the
     * #effectiveCouplingInput first and then use these functions to adjust or
     * add decay modes.*
     */
    /**
     * Set the decay width in the given channel to the specified value (in
     * GeV). Automatically adjusts the total width and all other BRs to
     * match.
     * @throws InvalidInput for negative width values
     */
    void setDecayWidth(Decay d, double value);

    /**
     * Set the decay width for the chain decay to the specified value.
     * Automatically adjusts the total width and all other BRs to match.
     * @throws InvalidInput for negative width values
     * @param d the kind of decay chain
     * @param particleId the BSM daughter particle
     * @param value the partial width in GeV
     */
    void setDecayWidth(ChainDecay d, const std::string &particleId,
                       double value);

    /**
     * Set the decay width for the decay into two BSM particles. Automatically
     * adjusts the total width and all other BRs to match.
     * @throws InvalidInput for negative width values
     * @param particleId1 id of the first daughter particle
     * @param particleId2 id of the second daughter particle, the order of these
     * two arguments is irrelevant
     * @param value the partial width in GeV
     */
    void setDecayWidth(const std::string &particleId1,
                       const std::string &particleId2, double value);
    // @}



    /** @{ @name Explicitely setting channel rates
     * @brief By default HiggsPredictions generates model predictions in the
     * narrow width approximation through `cxn*br`. If this is not wanted, it is
     * also possible to explicitely specify the inclusive rates for complete `SM
     * inital state -> BSM resonance -> SM final state` channels. If such a
     * channel rate is set fo a specific channel, it will always be used instead
     * of the corresponding narrow width value.
     */
    //! Set the channelrate for the specified collider channel. If an explicit
    //! channel rate is set, the narrow width approximation is *not* used for
    //! the channel. In contrast to the setBr() and setDecayWidth() functions,
    //! this function **will not try to validate your value in any
    //! way or form**.
    void setChannelRate(Collider coll, Production prod, Decay decay,
                        double value);

    //! Clear any explicitely set channelrates
    void resetChannelRates() noexcept;
    // @}

    /** @{ @name Setting couplings
     * @brief In some cases it is not entirely possible to work on the
     * model-independent level of cross sections and branching ratios. This
     * happens e.g. for scenarios where signal-background interference is
     * important in searches for new particles (e.g. in multi-top final states)
     * or for measurements of CP-phases where the final observable is a
     * specifically defined coupling phase. To handle these cases,
     * HiggsPredictions also allows setting values for specific predefined
     * couplings.
     *
     * **This is not equivalent to the #effectiveCouplingInput()**, but if you
     * use the #effectiveCouplingInput for a supported particle, all of these
     * couplings will be set automatically.
     */
    /**
     * @brief Set a coupling to the given value.
     *
     * **This does not set any cross sections or branching ratios that might be
     * dependent on this coupling.** For that kind of functionality use the
     * Higgs::Predictions::effectiveCouplingInput, which uses a complete set of
     * effective couplings to both set all possible cross sections and branching
     * ratios and store the relevant couplings explicitely through this method.
     *
     * @param c which coupling
     * @param value the value, the unit depends on which coupling is set
     */
    void setCoupling(Coupling c, double value);
    //! @}

    /** @privatesection */
    //! get the production cxn at the given collider
    double cxn(Collider coll, Production p) const noexcept override;

    //! get the total width in GeV
    double totalWidth() const noexcept override;

    //! Access explicitely set couplings
    std::optional<double> coupling(Coupling c) const noexcept override;

    //! get the BR for the given decay
    double br(Decay d) const noexcept override;
    double br(ChainDecay d,
              const std::string &particleId) const noexcept override;
    double br(const std::string &particleId1,
              const std::string &particleId2) const noexcept override;

    //! move/copy operations and destructor (value semantics)
    BsmParticle &operator=(BsmParticle &&other) noexcept;
    BsmParticle(BsmParticle &&other) noexcept;
    BsmParticle &operator=(const BsmParticle &other);
    BsmParticle(const BsmParticle &other);
    std::unique_ptr<Particle> clone() const override;
    ~BsmParticle() noexcept;

    //! @internal
    //! This is necessary so that the Higgs::predictions::Predictions can let a
    //! particle know about the other particles around it.
    friend class ::Higgs::predictions::Predictions;

  protected:
    //! Makes it so that channel rates set with setChannelRate() override the
    //! narrow width approximation.
    double channelRateOr(Collider coll, Production prod, Decay decay,
                         double cxnTimesBr) const noexcept override;

  private:
    std::unique_ptr<ParticleData> data_;
};

//! Ordering predicate for particles.
class ParticleIdOrder {
  public:
    //! Orders two particles by their id.
    bool operator()(const Particle &p1, const Particle &p2) const noexcept {
        return p1.id() < p2.id();
    }
};

//! A reference to a Particle
using ParticleRef = std::reference_wrapper<const Particle>;

//! A set of references to different particles
using ParticleSet = std::set<ParticleRef, ParticleIdOrder>;

} // namespace predictions
} // namespace Higgs
