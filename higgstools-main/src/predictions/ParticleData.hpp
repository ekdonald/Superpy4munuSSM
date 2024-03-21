/**
 * @file ParticleData.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Pimpl structs for the Higgs::predictions::Particle classes
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <magic_enum.hpp>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace Higgs::predictions {
class Predictions;

//! Pimpl data for the particle classes
class ParticleData {
  public:
    //! set a cxn
    void setCxn(Collider coll, Production p, double value) noexcept;

    //! get a cxn
    double cxn(Collider coll, Production p) const noexcept;

    //! set a BR
    void setBr(Decay d, double value);

    //! set a BSM BR
    void setBr(std::pair<std::string, std::string> &&particleIds, double value);

    //! set a BR through a decay width
    void setDecayWidth(Decay d, double value);

    //! set a BSM BR through a decay width
    void setDecayWidth(std::pair<std::string, std::string> &&particleIds,
                       double value);

    //! get a BR
    double br(Decay d) const noexcept;

    //! Get the full br into invisible, including chain and pair decays with
    //! other particles lighter than mass (if context is available).
    double fullBrInv(double mass) const;

    //! get a BSM BR
    double
    br(const std::pair<std::string, std::string> &particleIds) const noexcept;

    /**
     * Set the total width to the value in GeV.
     * Setting the total width to zero automatically also zeros all BRs.
     * @throws InvalidInput for negative values
     */
    void setTotalWidth(double value);

    //! total width access
    double totalWidth() const noexcept;

    std::optional<double> coupling(Coupling c) const noexcept;

    void setCoupling(Coupling c, double value);

    //! Set the channelrate for the given collider channel.
    void setChannelRate(Collider coll, Production p, Decay d, double value);

    //! clear any explicitely set channelrates
    void resetChannelRates() noexcept;

    //! the rate for the specified channel, or the default value, if no explicit
    //! channel rate is set
    double channelRateOr(Collider coll, Production p, Decay d,
                         double defaultValue) const noexcept;

    // checks and updates the br to the new value
    void updateBr(double &br, double value);
    // updates a decay width and rescales all BRs
    void updateBrFromDecayWidth(double &br, double width);

    //! Key class which governs access to the setContext method
    class ContextKey {
        friend class ::Higgs::predictions::Predictions;
        ContextKey() = default;
    };

    //! Set a predictions class as context for this particle. Can only
    //! be called by classes with access to the ContextKey (i.e. only from
    //! within Higgs::predictions::Predictions)
    void setContext(const Predictions &context,
                    ContextKey /* unused */) noexcept;

    void clearContext() noexcept;

  private:
    //! total width in GeV
    double wtot_ = 0.;
    //! branching ratios
    std::unordered_map<Decay, double> brs_ = {};
    //! BSM BRs
    std::map<std::pair<std::string, std::string>, double> bsmBrs_ = {};
    //! buffered sum of BRs
    double sumOfBrs_ = 0.;
    //! hadron collider production cross sections in pb
    std::map<std::pair<Collider, Production>, double> cxns_ = {};
    //! couplings
    std::unordered_map<Coupling, double> coups_ = {};

    using CollChannel = std::tuple<Collider, Production, Decay>;

    struct CollChannelHash {
        std::size_t operator()(const CollChannel &collChannel) const noexcept {
            return (static_cast<std::size_t>(std::get<0>(collChannel))
                    << (std::numeric_limits<std::size_t>::digits * 2 / 3)) +
                   (static_cast<std::size_t>(std::get<1>(collChannel))
                    << (std::numeric_limits<std::size_t>::digits / 3)) +
                   static_cast<std::size_t>(std::get<2>(collChannel));
        }
    };

    //! explicitely specified channelrates
    std::unordered_map<CollChannel, double, CollChannelHash> channelRates_ = {};

    //! Allows to reference the Predictions that this particle is a part
    //! of. This is reset to nullptr on copies/moves of the BsmParticle.
    const Predictions *context_ = nullptr;
};

inline auto chainDecayKey(ChainDecay d, const std::string &particleId) {
    return std::pair<std::string, std::string>{magic_enum::enum_name(d),
                                               particleId};
}

inline auto twoBsmDecayKey(const std::string &particleId1,
                           const std::string &particleId2) {
    return std::pair<std::string, std::string>{
        std::minmax(particleId1, particleId2)};
}

} // namespace Higgs::predictions
