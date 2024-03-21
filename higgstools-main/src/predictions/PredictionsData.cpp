#include "predictions/PredictionsData.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "utilities/Algorithm.hpp"
#include "utilities/Logging.hpp"
#include <exception>
#include <functional>
#include <memory>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <utility>

namespace Higgs::predictions {

void PredictionsData::ensureUnusedId(const std::string &newId) const {
    if (findParticle(newId) != particles_.end()) {
        throw predictions::InvalidInput(
            fmt::format("A scalar with id {} already exists, "
                        "cannot insert duplicate.",
                        newId));
    }
}

std::list<BsmParticle>::iterator
PredictionsData::findParticle(const std::string &id) noexcept {
    return ranges::find(particles_, id,
                        std::mem_fn(&predictions::Particle::id));
}

std::list<BsmParticle>::const_iterator
PredictionsData::findParticle(const std::string &id) const noexcept {
    return ranges::find(particles_, id,
                        std::mem_fn(&predictions::Particle::id));
}

BsmParticle &PredictionsData::addParticle(BsmParticle &&particle) {
    ensureUnusedId(particle.id());
    particles_.push_back(std::move(particle));
    return particles_.back();
}

void PredictionsData::removeParticle(const std::string &id) {
    if (auto loc = findParticle(id); loc != particles_.end()) {
        particles_.erase(loc);
    }
}
void PredictionsData::setPairCxn(Collider coll,
                                 std::pair<std::string, std::string> &&pairId,
                                 double value) {
    pairCxns_[{coll, std::move(pairId)}] = value;
}
double PredictionsData::pairCxn(
    Collider coll,
    std::pair<std::string, std::string> &&pairId) const noexcept {
    return utilities::getWithDefault(pairCxns_, {coll, std::move(pairId)}, 0.);
}

namespace {
template <class ParticleRange>
double sumBrTopBsm(const ParticleRange &particles) {
    constexpr auto getBrtBsm = [](const auto &p) {
        return p.cxn(predictions::Collider::LHC13,
                     predictions::Production::brtHpb) +
               p.cxn(predictions::Collider::LHC13,
                     predictions::Production::brtHc) +
               p.cxn(predictions::Collider::LHC13,
                     predictions::Production::brtHu);
    };
    return ranges::accumulate(particles, 0., std::plus{}, getBrtBsm);
}
} // namespace

double PredictionsData::brTopWb() const noexcept {
    if (brTopWb_ > 0) {
        return brTopWb_;
    }
    const auto brTopBsm = sumBrTopBsm(particles_);
    if (brTopBsm > 1) {
        predictions::logger()->error(
            "The combined BR(t > H^+ b) + BR(t > H q) = {} exceeds 1.",
            brTopBsm);
        return 0.;
    }
    return 1 - brTopBsm;
}

void PredictionsData::setBrTopWb(double value) {
    const auto brTopBsm = sumBrTopBsm(particles_);
    if (value + brTopBsm > 1) {
        throw predictions::InvalidInput(fmt::format(
            "Explicitely setting BR(t > W^+ b) to {} would result in "
            "sum(BR(t))={}>1",
            value, brTopBsm));
    }
    brTopWb_ = value;
}

} // namespace Higgs::predictions
