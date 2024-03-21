#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "predictions/ParticleData.hpp"
#include "predictions/PredictionsData.hpp"
#include "utilities/Format.hpp"
#include <functional>
#include <memory>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Higgs {

Predictions::~Predictions() = default;
Predictions::Predictions() noexcept
    : data_{std::make_unique<predictions::PredictionsData>()} {}

std::vector<std::string> Predictions::particleIds() const noexcept {
    return data_->particles() |
           ranges::views::transform(std::mem_fn(&predictions::Particle::id)) |
           ranges::to<std::vector>;
} // LCOV_EXCL_LINE

predictions::BsmParticle &
Predictions::addParticle(predictions::BsmParticle particle) {
    auto &p = data_->addParticle(std::move(particle));
    p.data_->setContext(*this, {});
    return p;
}

predictions::BsmParticle &Predictions::particle(const std::string &id) {
    // forward to const overload (see EffCpp Item 3)
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
    return const_cast<predictions::BsmParticle &>(
        std::as_const(*this).particle(id));
}

const predictions::BsmParticle &
Predictions::particle(const std::string &id) const {
    auto loc = data_->findParticle(id);
    if (loc == data_->particles().end()) {
        throw std::out_of_range(
            fmt::format("No particle with the requested id {} exists.", id));
    }
    return *loc;
}

void Predictions::removeParticle(const std::string &id) {
    data_->removeParticle(id);
}

namespace {
auto twoParticleKey(const std::string &id1, const std::string &id2) {
    return id1 < id2 ? std::pair{id1, id2} : std::pair{id2, id1};
}
} // namespace

double Predictions::bsmPairCxn(predictions::Collider coll,
                               const std::string &id1,
                               const std::string &id2) const noexcept {
    if (auto loc1 = data_->findParticle(id1);
        loc1 != data_->particles().end()) {
        if (id1 == id2) {
            return loc1->cxn(coll, predictions::Production::pair);
        }
        if (auto loc2 = data_->findParticle(id2);
            loc2 != data_->particles().end()) {
            return data_->pairCxn(coll, twoParticleKey(id1, id2));
        }
    }
    return 0;
}

void Predictions::setBsmPairCxn(predictions::Collider coll,
                                const std::string &id1, const std::string &id2,
                                double value) {
    if (id1 == id2) {
        if (coll == predictions::Collider::LEP) {
            particle(id1).setNormalizedCxn(
                coll, predictions::Production::pair, value,
                predictions::ReferenceModel::SMHiggs);
        } else {
            particle(id1).setCxn(coll, predictions::Production::pair, value);
        }
    } else {
        // check that both particles exist, throw is not
        particle(id1);
        particle(id2);

        return data_->setPairCxn(coll, twoParticleKey(id1, id2), value);
    }
}

double Predictions::brTopWb() const noexcept { return data_->brTopWb(); }
void Predictions::setBrTopWb(double value) { data_->setBrTopWb(value); }

const std::list<predictions::BsmParticle> &
Predictions::particles() const noexcept {
    return data_->particles();
}

Predictions &Predictions::operator=(Predictions &&other) noexcept {
    data_ = std::move(other.data_);
    for (auto &p : data_->particles()) {
        p.data_->setContext(*this, {});
    }
    return *this;
}
Predictions::Predictions(Predictions &&other) noexcept
    : data_{std::move(other.data_)} {
    for (auto &p : data_->particles()) {
        p.data_->setContext(*this, {});
    }
}
Predictions &Predictions::operator=(const Predictions &other) {
    if (this != &other) {
        data_ = std::make_unique<predictions::PredictionsData>(*other.data_);
        for (auto &p : data_->particles()) {
            p.data_->setContext(*this, {});
        }
    }
    return *this;
}
Predictions::Predictions(const Predictions &other)
    : data_{std::make_unique<predictions::PredictionsData>(*other.data_)} {
    for (auto &p : data_->particles()) {
        p.data_->setContext(*this, {});
    }
}
} // namespace Higgs
