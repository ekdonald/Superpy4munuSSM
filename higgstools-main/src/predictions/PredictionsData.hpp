#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Particle.hpp"
#include <list>
#include <map>
#include <string>
#include <utility>

namespace Higgs::predictions {

class PredictionsData {
  public:
    void ensureUnusedId(const std::string &newId) const;

    BsmParticle &addParticle(BsmParticle &&particle);

    void removeParticle(const std::string &id);

    const auto &particles() const noexcept { return particles_; }

    std::list<BsmParticle>::iterator
    findParticle(const std::string &id) noexcept;
    std::list<BsmParticle>::const_iterator
    findParticle(const std::string &id) const noexcept;

    double brTopWb() const noexcept;
    void setBrTopWb(double value);

    void setPairCxn(Collider coll, std::pair<std::string, std::string> &&pairId,
                    double value);
    double pairCxn(Collider coll,
                   std::pair<std::string, std::string> &&pairId) const noexcept;

  private:
    std::list<BsmParticle> particles_;
    double brTopWb_;
    std::map<std::pair<Collider, std::pair<std::string, std::string>>, double>
        pairCxns_;
};
} // namespace Higgs::predictions
