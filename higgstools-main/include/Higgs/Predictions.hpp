/**
 * @file Predictions.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * Interface of the HiggsPredictions library. Includes all other headers
 * needed for access to all of HiggsPredictions functionality.
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsTools.
 * HiggsTools is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
// IWYU pragma: begin_exports
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
// IWYU pragma: end_exports
#include <list>
#include <memory>
#include <string>
#include <vector>

namespace Higgs {
namespace predictions {
class PredictionsData;

/**
 * The central class of the HiggsPredictions library.
 * Stores and gives access to model predictions for all particles in the model.
 */
class HIGGSTOOLS_EXPORT Predictions {
  public:
    /** @{ @name Particle management
     * @brief The main purpose of the Predictions class is to assemble and store
     * the model predictions for all of the BsmParticle%s in the model. The
     * following functions are responsible for adding, retreiving, and removing
     * such particles.
     */

    //! Construct an empty predictions object. Particles can then be added
    //! using addParticle().
    explicit Predictions() noexcept;

    /**
     * Add the given particle to these predictions.
     *
     * It it highly recommended to immediately add particles to a prediction
     * when you create them and then set their properties through the returned
     * reference.
     *
     * @throws InvalidInput if the ID of the particle is already in use.
     * @param particle the particle to be added
     * @return a reference to the added particle
     */
    BsmParticle &addParticle(BsmParticle particle);

    //! returns a list of the IDs of all contained particles
    std::vector<std::string> particleIds() const noexcept;

    /**
     * Returns a particle by looking up the ID.
     *
     * Do not assign to the returned references and use the member functions of
     * the particle to modify it instead. **Never assign a particle with a
     * different ID to the returned reference.**
     *
     * @throws std::out_of_range if no particle with the ID is found
     * @param id the id to look for
     * @return  the corresponding particle
     */
    BsmParticle &particle(const std::string &id);

    /**
     * Returns a particle by looking up the ID.
     *
     * @throws std::out_of_range if no particle with the ID is found
     * @param id the id to look for
     * @return  the corresponding particle
     */
    const BsmParticle &particle(const std::string &id) const;

    /**
     * Remove a particle based on the ID.
     * @param id the ID to look for, no effect if the id does not exist
     */
    void removeParticle(const std::string &id);
    //!@}

    /** @{ @name Multi-Particle production processes
     * @brief The Predictions object also has to handle all aspects of the model
     * that cannot be associated with a single BsmParticle and stored there.
     * First of all this includes non-resonant pair production of (potentially)
     * distinct particles. This is handled by the following functions.
     */
    /**
     * @brief Access a cross section for non-resonant pair production of BSM
     * particles.
     *
     * Non-resonant in this context means, without contributions from on-shell
     * \f$s\f$-channel BSM particles. Such contributions should be modeled
     * through the production and decay of the intermediate particle instead.
     *
     * If `id1 == id2` this is equivalent to
     * `this->particle(id1).cxn(coll, Production::pair)`.
     *
     * @param coll the collider where the production takes place.
     * @param id1 id of the first produced BSM particle
     * @param id2 id of the second produced BSM particle, the order of id1 and
     * id2 is irrelevant
     * @return double value of the cross section. For the LHC an absolute cross
     * section value in pb. For LEP the cxn is normalized depending on its
     * charge:
     *
     *  |     charge()     | reference process
     *  |------------------|-------------------
     *  | ECharge::neutral | \f$ e^+ e^- \to Z \to A H \f$
     *  | ECharge::single  | \f$ e^+ e^- \to Z \to H^+ H^- \f$
     *
     * where \f$ H,A,H^\pm \f$ are the non-SM-like CP-even, CP-odd, and charged
     * scalars of a 2HDM where the other CP-even scalar \f$ h \f$ is SM-like
     * (exact alignment limit). This normalization is necessary since LEP
     * results combine data at many different CM energies.
     */
    double bsmPairCxn(Collider coll, const std::string &id1,
                      const std::string &id2) const noexcept;

    /**
     * @brief Set a cross section for non-resonant pair production of BSM
     * particles.
     *
     * Non-resonant in this context means, without contributions from on-shell
     * \f$s\f$-channel BSM particles.
     *
     * If `id1 == id2` this is equivalent to
     * `this.particle(id1).setCxn(coll, Production::pair)`.
     *
     * @throws InvalidInput if either id1 or id2 does not name a
     * contained particle
     *
     * @param coll the #Collider where the production takes place.
     * @param id1 id of the first produced BSM particle
     * @param id2 id of the second produced BSM particle, the order of id1 and
     * id2 is irrelevant
     * @param value see the documentation of the return value of bsmPairCxn()
     */
    void setBsmPairCxn(Collider coll, const std::string &id1,
                       const std::string &id2, double value);
    //! @}

    /** @{ @name Top-quark properties
     * @brief Additionally, there are also some properties of the top-quark that
     * can be relevant for charged particles produced in top-decays.
     */
    //! Access \f$\mathrm{BR}(t\to W^+ b)\f$. If not set, this defaults to 1 -
    //! the Production::brtHpb values for all contained charged scalars. In most
    //! models, this is correct and setting this manually is not required.
    double brTopWb() const noexcept;
    //! Explicitely set the \f$\mathrm{BR}(t\to W^+ b)\f$. You should never have
    //! to do this unless your model has BSM top-quark decays into particles
    //! that are not represented by BsmParticle%s. For top-quark decays into
    //! BsmParticle%s set a pseudo-cxn for the Production::brtHpb mode instead.
    //! @throws InvalidInput if setting this would result in a sum(top-BRs)>1
    void setBrTopWb(double value);
    //! @}

    //! @privatesection
    /**
     * Access all contained particles.
     *
     * The order of the particles in the list is unspecified but will not
     * change unless particles are added or removed. This routine is mostly
     * intended for internal use.
     *
     * @return the list of particles
     */
    const std::list<BsmParticle> &particles() const noexcept;

    //! move/copy operations and destructor (value semantics)
    Predictions &operator=(Predictions &&other) noexcept;
    Predictions(Predictions &&other) noexcept;
    Predictions &operator=(const Predictions &other);
    Predictions(const Predictions &other);
    ~Predictions();

  private:
    std::unique_ptr<PredictionsData> data_;
};
} // namespace predictions
using predictions::Predictions;
} // namespace Higgs
