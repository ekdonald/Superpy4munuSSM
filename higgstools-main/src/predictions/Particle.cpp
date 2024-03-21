#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "predictions/ParticleData.hpp"
#include "utilities/Logging.hpp"
#include <exception>
#include <magic_enum.hpp>
#include <memory>
#include <tuple>
#include <utility>

namespace Higgs::predictions {
Particle::Particle(std::string &&id, CP cp, ECharge charge)
    : id_{std::move(id)}, mass_{0}, massUnc_{0}, cp_{cp}, charge_{charge} {}

const std::string &Particle::id() const noexcept { return id_; }

CP Particle::cp() const noexcept { return cp_; }
ECharge Particle::charge() const noexcept { return charge_; }

double Particle::mass() const noexcept { return mass_; }
void Particle::setMass(double value) noexcept { mass_ = value; }
double Particle::massUnc() const noexcept { return massUnc_; }
void Particle::setMassUnc(double value) noexcept { massUnc_ = value; }

double Particle::totalWidth() const noexcept { return 0; }

double Particle::br(Decay /* d */) const noexcept { return 0.; }
double Particle::br(ChainDecay /*chain*/,
                    const std::string & /*particleId*/) const noexcept {
    return 0.;
}
double Particle::br(const std::string & /*particleId1*/,
                    const std::string & /*particleId2*/) const noexcept {
    return 0.;
}

double Particle::cxn(Collider /* coll */, Production /* p */) const noexcept {
    return 0.;
}

double Particle::channelRate(Collider coll, Production p,
                             Decay d) const noexcept {
    const auto b = br(d);
    if (p == Production::none) {
        return b;
    }
    const auto x = cxn(coll, p);
    if (d == Decay::none) {
        return x;
    }
    return channelRateOr(coll, p, d, x * b);
}

double Particle::channelRateOr(Collider /*coll*/, Production /*p*/, Decay /*d*/,
                               double cxnTimesBr) const noexcept {
    return cxnTimesBr;
}

std::optional<double> Particle::coupling(Coupling /* c */) const noexcept {
    return std::nullopt;
}

bool Particle::operator==(const Particle &other) const noexcept {
    return std::tie(id_, charge_, cp_, mass_, massUnc_) ==
           std::tie(other.id_, other.charge_, other.cp_, other.mass_,
                    other.massUnc_);
}

BsmParticle::BsmParticle(std::string id, ECharge charge, CP cp)
    : Particle{std::move(id), cp, charge},
      data_{std::make_unique<ParticleData>()} {}

BsmParticle::~BsmParticle() noexcept = default;
BsmParticle &BsmParticle::operator=(BsmParticle &&other) noexcept {
    if (this != &other) {
        Particle::operator=(std::move(other));
        data_ = std::move(other.data_);
        data_->clearContext();
    }
    return *this;
}

BsmParticle::BsmParticle(BsmParticle &&other) noexcept
    : Particle{std::move(other)}, data_{std::move(other.data_)} {
    data_->clearContext();
}
BsmParticle &BsmParticle::operator=(const BsmParticle &other) {
    if (this != &other) {
        Particle::operator=(other);
        data_ = std::make_unique<ParticleData>(*other.data_);
        data_->clearContext();
    }
    return *this;
}
BsmParticle::BsmParticle(const BsmParticle &other)
    : Particle{other}, data_{std::make_unique<ParticleData>(*other.data_)} {
    data_->clearContext();
}
std::unique_ptr<Particle> BsmParticle::clone() const {
    return std::make_unique<BsmParticle>(*this);
}

namespace {
void checkProductionMode(Production p, Collider coll, ECharge charge) {
    if (!validProductionFor(p, charge)) {
        throw(InvalidChannel{charge, p});
    }
    if (auto collType = classifyCollider(coll);
        !validProductionAt(p, collType)) {
        throw(InvalidChannel{collType, p});
    }
}
} // namespace

void BsmParticle::setCxn(Collider coll, Production p, double value) {
    checkProductionMode(p, coll, charge());
    if (coll == Collider::LEP) {
        throw(InvalidInput("Can't set absolute LEP cxns due to varying CM "
                           "energies. Use setNormalizedCxn instead."));
    }
    switch (p) {
    default:
        data_->setCxn(coll, p, value);
        return;
    // ---- combined production modes ----
    case Production::H:
        throw(InvalidInput{"Can't set combined production mode H = ggH + bbH + "
                           "qqH, set the underlying modes instead."});
    case Production::HZ:
        throw(InvalidInput{"Can't set combined production mode HZ = qqHZ + "
                           "ggHZ + bbHZ, set the underlying modes instead."});
    case Production::Ht:
        throw(InvalidInput{"Can't set combined production mode Ht = tchanHt + "
                           "schanHt, set the underlying modes instead."});
        // ---- production modes that aren't cross sections ----
    case Production::brtHpb:
    case Production::brtHc:
    case Production::brtHu:
        data_->setCxn(Collider::LHC13, p, value);
        return;
    }
}

void BsmParticle::setNormalizedCxn(Collider coll, Production p, double value,
                                   ReferenceModel reference) {
    checkProductionMode(p, coll, charge());
    switch (coll) {
    default: {
        const auto ref = getReference(reference, mass());
        setCxn(coll, p, value * ref->cxn(coll, p));
        return;
    }
    case Collider::LEP:
        // LEP cxns are always normalized, so just store them as given
        data_->setCxn(coll, p, value);
        return;
    }
}

double BsmParticle::cxn(Collider coll, Production p) const noexcept {
    switch (p) {
    default:
        return data_->cxn(coll, p);
    // ---- combined production modes ----
    case Production::H:
        return data_->cxn(coll, Production::ggH) +
               data_->cxn(coll, Production::bbH) +
               data_->cxn(coll, Production::qqH);
    case Production::HZ:
        return data_->cxn(coll, Production::qqHZ) +
               data_->cxn(coll, Production::ggHZ) +
               data_->cxn(coll, Production::bbHZ);
    case Production::Ht:
        return data_->cxn(coll, Production::tchanHt) +
               data_->cxn(coll, Production::schanHt);
    // ---- production modes that are not cross sections ----
    case Production::brtHpb:
    case Production::brtHc:
    case Production::brtHu:
        return data_->cxn(Collider::LHC13, p);
    }
}

void BsmParticle::setBr(Decay d, double value) {
    if (!validDecayFor(d, charge())) {
        throw(InvalidChannel{charge(), d});
    }
    try {
        switch (d) {
        default:
            data_->setBr(d, value);
            return;
        case Decay::inv:
            logger()->info("Can't set combined decay mode inv, setting "
                           "directInv instead.");
            data_->setBr(Decay::directInv, value);
            return;
        }
    } catch (const InvalidInput &e) {
        logger()->error("Error setting BR->{} of {}", magic_enum::enum_name(d),
                        id());
        throw;
    }
}

void BsmParticle::setBr(ChainDecay d, const std::string &particleId,
                        double value) {
    try {
        data_->setBr(chainDecayKey(d, particleId), value);
    } catch (const InvalidInput &e) {
        logger()->error("Error setting chain decay BR->{} (X = {}) of {}: {}",
                        magic_enum::enum_name(d), particleId, id(), e.what());
        throw;
    }
}

void BsmParticle::setBr(const std::string &particleId1,
                        const std::string &particleId2, double value) {
    try {
        data_->setBr(twoBsmDecayKey(particleId1, particleId2), value);
    } catch (const InvalidInput &e) {
        logger()->error("Error setting chain decay BR->{} {} of {}",
                        particleId1, particleId2, id());
        throw;
    }
}

void BsmParticle::setDecayWidth(Decay d, double value) {
    if (!validDecayFor(d, charge())) {
        throw(InvalidChannel{charge(), d});
    }
    try {
        switch (d) {
        default:
            data_->setDecayWidth(d, value);
            return;
        case Decay::inv:
            logger()->info("Can't set combined decay mode inv, setting "
                           "directInv instead.");
            data_->setDecayWidth(Decay::directInv, value);
            return;
        }
    } catch (const InvalidInput &e) {
        logger()->error("Error setting decay width into {} of {}: {}",
                        magic_enum::enum_name(d), id(), e.what());
        throw;
    }
}

void BsmParticle::setDecayWidth(ChainDecay d, const std::string &particleId,
                                double value) {
    try {
        data_->setDecayWidth(chainDecayKey(d, particleId), value);
    } catch (const InvalidInput &e) {
        logger()->error(
            "Error setting chain decay width into {} (X = {}) of {}: {}",
            magic_enum::enum_name(d), particleId, id(), e.what());
        throw;
    }
}

void BsmParticle::setDecayWidth(const std::string &particleId1,
                                const std::string &particleId2, double value) {
    try {
        data_->setDecayWidth(twoBsmDecayKey(particleId1, particleId2), value);
    } catch (const InvalidInput &e) {
        logger()->error("Error setting chain decay width into {} {} of {}: {}",
                        particleId1, particleId2, id(), e.what());
        throw;
    }
}

double BsmParticle::br(Decay d) const noexcept {
    switch (d) {
    default:
        return data_->br(d);
    case Decay::inv:
        if (charge() == ECharge::neutral) {
            return data_->fullBrInv(mass());
        } else {
            return 0;
        }
    }
}

double BsmParticle::br(ChainDecay d,
                       const std::string &particleId) const noexcept {
    return data_->br(chainDecayKey(d, particleId));
}

double BsmParticle::br(const std::string &particleId1,
                       const std::string &particleId2) const noexcept {
    return data_->br(twoBsmDecayKey(particleId1, particleId2));
}

void BsmParticle::setTotalWidth(double value) {
    try {
        data_->setTotalWidth(value);
    } catch (const InvalidInput &e) {
        logger()->error("Error setting total width of {}", id());
        throw;
    }
}

double BsmParticle::totalWidth() const noexcept { return data_->totalWidth(); }

void BsmParticle::setChannelRate(Collider coll, Production prod, Decay decay,
                                 double value) {
    data_->setChannelRate(coll, prod, decay, value);
}

void BsmParticle::resetChannelRates() noexcept { data_->resetChannelRates(); }

double BsmParticle::channelRateOr(Collider coll, Production prod, Decay decay,
                                  double cxnTimesBr) const noexcept {
    return data_->channelRateOr(coll, prod, decay, cxnTimesBr);
}

std::optional<double> BsmParticle::coupling(Coupling c) const noexcept {
    return data_->coupling(c);
}

void BsmParticle::setCoupling(Coupling c, double value) {
    data_->setCoupling(c, value);
}

} // namespace Higgs::predictions
