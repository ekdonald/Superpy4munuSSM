#include "Higgs/predictions/Process.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/FormatSupport.hpp" // IWYU pragma: keep
#include "predictions/JsonSupport.hpp"   // IWYU pragma: keep
#include "utilities/Logging.hpp"
#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <range/v3/action/action.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/detail/variant.hpp>
#include <range/v3/functional/arithmetic.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/repeat.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip_with.hpp>
#include <tuple>
#include <utility>
#include <vector>

namespace Higgs::predictions {

// --------- process constructors and properties -------- //
ChannelProcess::ChannelProcess(Collider coll, std::vector<Channel> channels)
    : coll_{coll}, channels_{std::move(channels)} {}

Collider ChannelProcess::collider() const noexcept { return coll_; }

std::string ChannelProcess::to_string(bool compactify) const noexcept {
    if (compactify && channels_.size() > 1) {
        auto prodModes =
            channels_ | ranges::views::keys | ranges::to<std::vector>;
        prodModes |= ranges::actions::sort | ranges::actions::unique;

        auto decModes =
            channels_ | ranges::views::values | ranges::to<std::vector>;
        decModes |= ranges::actions::sort | ranges::actions::unique;
        if (channels_.size() == prodModes.size() * decModes.size()) {
            return fmt::format("{} [{}]>[{}]", collider(),
                               fmt::join(prodModes, ","),
                               fmt::join(decModes, ","));
        }
    }
    auto chans = channels_ | ranges::views::transform([](const auto &c) {
                     return fmt::format("{}>{}", c.first, c.second);
                 });
    return fmt::format("{} [{}]", collider(), fmt::join(chans, " + "));
}

ChannelProcess &ChannelProcess::operator+=(const ChannelProcess &other) {
    if (other.channels_.empty()) {
        return *this;
    }
    if (channels_.empty()) {
        return *this = other;
    }
    if (collider() == other.collider()) {
        channels_.insert(channels_.end(), other.channels_.begin(),
                         other.channels_.end());
    } else {
        logger()->warn(
            "Cannot merge ChannelProcesses, mismatch in colliders ({} vs {})",
            collider(), other.collider());
    }
    return *this;
}

std::vector<ChannelProcess> ChannelProcess::subprocesses() const noexcept {
    const auto constructSubchannel = [coll = collider()](const auto &chan) {
        return ChannelProcess{coll, std::vector{chan}};
    };
    return ranges::views::zip_with(constructSubchannel, channels_) |
           ranges::to<std::vector>;
}

bool ChannelProcess::operator==(const ChannelProcess &other) const noexcept {
    return coll_ == other.coll_ && channels_ == other.channels_;
} // LCOV_EXCL_LINE

ChannelProcess operator+(ChannelProcess lhs, const ChannelProcess &rhs) {
    return lhs += rhs;
}

std::size_t ChannelProcess::size() const noexcept { return channels_.size(); }

ChainDecayProcess::ChainDecayProcess(ChainDecay chain,
                                     std::vector<Production> production,
                                     std::vector<Decay> decay, Collider coll)
    : chain_{chain}, production_{std::move(production)},
      decay_{std::move(decay)}, coll_{coll} {}

ChainDecay ChainDecayProcess::chain() const noexcept { return chain_; }
Collider ChainDecayProcess::collider() const noexcept { return coll_; }

std::string ChainDecayProcess::to_string() const noexcept {
    return fmt::format(
        "{} [{}]->X1->{}(X2->[{}])", coll_, fmt::join(production_, ", "),
        chain_ == ChainDecay::Z ? "Z" : "W", fmt::join(decay_, ", "));
}

std::size_t ChainDecayProcess::productionSize() const noexcept {
    return production_.size();
}
std::size_t ChainDecayProcess::decaySize() const noexcept {
    return decay_.size();
}

namespace {
std::set<std::pair<Decay, Decay>>
generateUnorderedDecayPairs(const std::vector<Decay> &decay1,
                            const std::vector<Decay> &decay2) {
    constexpr auto toUnorderedPair =
        [](const std::tuple<Decay, Decay> &decayPair) {
            const auto &[d1, d2] = decayPair;
            return std::pair<Decay, Decay>{std::minmax(d1, d2)};
        };
    return ranges::views::cartesian_product(decay1, decay2) |
           ranges::views::transform(toUnorderedPair) | ranges::to<std::set>;
}
} // namespace

PairDecayProcess::PairDecayProcess(std::vector<Production> production,
                                   std::vector<Decay> firstDecay,
                                   std::vector<Decay> secondDecay,
                                   Collider coll)
    : production_{std::move(production)}, decay1_{std::move(firstDecay)},
      decay2_{std::move(secondDecay)},
      unorderedDecayPairs_{generateUnorderedDecayPairs(decay1_, decay2_)},
      coll_{coll} {}

Collider PairDecayProcess::collider() const noexcept { return coll_; }

std::string PairDecayProcess::to_string() const noexcept {
    return fmt::format("{} [{}]->X1->(X2->[{}])(X3->[{}])", coll_,
                       fmt::join(production_, ", "), fmt::join(decay1_, ", "),
                       fmt::join(decay2_, ", "));
}

PairProductionProcess::PairProductionProcess(std::vector<Decay> firstDecay,
                                             std::vector<Decay> secondDecay,
                                             Collider coll)
    : decay1_{std::move(firstDecay)}, decay2_{std::move(secondDecay)},
      unorderedDecayPairs_{generateUnorderedDecayPairs(decay1_, decay2_)},
      coll_{coll} {}

Collider PairProductionProcess::collider() const noexcept { return coll_; }

std::string PairProductionProcess::to_string() const noexcept {
    return fmt::format("{} (X1->[{}])(X2->[{}])", coll_,
                       fmt::join(decay1_, ", "), fmt::join(decay2_, ", "));
}

// ---------------- evaluation operators ---------------- //
double ChannelProcess::operator()(const Particle &particle) const noexcept {
    const auto toChannelRate = [&particle, this](const auto &ch) {
        return particle.channelRate(collider(), ch.first, ch.second);
    };
    return ranges::accumulate(channels_, 0., std::plus{}, toChannelRate);
}

double
ChannelProcess::operator()(const Particle &particle,
                           const std::vector<double> &weights) const noexcept {
    const auto toChannelRate = [&particle, this](const auto &ch) {
        return particle.channelRate(collider(), ch.first, ch.second);
    };
    auto paddedWeights =
        ranges::views::concat(weights, ranges::views::repeat(1.));
    return ranges::inner_product(
        channels_ | ranges::views::transform(toChannelRate), paddedWeights, 0.);
}

double ChainDecayProcess::operator()(const Particle &mother,
                                     const Particle &daughter) const noexcept {

    const auto toCxn = [&mother, this](auto p) { return mother.cxn(coll_, p); };
    auto cxn = ranges::accumulate(production_, 0., std::plus{}, toCxn);

    const auto toBr = [&daughter](auto d) { return daughter.br(d); };
    auto decay = ranges::accumulate(decay_, 0., std::plus{}, toBr);

    return cxn * mother.br(chain_, daughter.id()) * decay;
}

double ChainDecayProcess::operator()(
    const Particle &mother, const Particle &daughter,
    const std::vector<double> &productionWeights) const noexcept {

    auto paddedWeights =
        ranges::views::concat(productionWeights, ranges::views::repeat(1.));

    const auto toCxn = [&mother, this](auto p) { return mother.cxn(coll_, p); };
    auto cxn = ranges::inner_product(
        production_ | ranges::views::transform(toCxn), paddedWeights, 0.);

    const auto toBr = [&daughter](auto d) { return daughter.br(d); };
    auto decay = ranges::accumulate(decay_, 0., std::plus{}, toBr);

    return cxn * mother.br(chain_, daughter.id()) * decay;
}

double
PairDecayProcess::operator()(const Particle &mother,
                             const Particle &firstDaughter,
                             const Particle &secondDaughter) const noexcept {

    const auto toCxn = [&mother, this](auto p) { return mother.cxn(coll_, p); };
    auto rate = ranges::accumulate(production_, 0., std::plus{}, toCxn);
    rate *= mother.br(firstDaughter.id(), secondDaughter.id());

    if (firstDaughter.id() != secondDaughter.id() ||
        (mother.charge() == ECharge::neutral &&
         firstDaughter.charge() != ECharge::neutral)) {
        const auto distinctPairBr = [&firstDaughter,
                                     &secondDaughter](const auto &decays) {
            const auto &[d1, d2] = decays;
            return firstDaughter.br(d1) * secondDaughter.br(d2);
        };
        rate *= ranges::accumulate(
            ranges::views::cartesian_product(decay1_, decay2_), 0., std::plus{},
            distinctPairBr);
    } else {
        const auto identicalPairBr = [&firstDaughter](const auto &decays) {
            const auto &[d1, d2] = decays;
            return (d1 == d2 ? 1 : 2) * firstDaughter.br(d1) *
                   firstDaughter.br(d2);
        };
        rate *= ranges::accumulate(unorderedDecayPairs_, 0., std::plus{},
                                   identicalPairBr);
    }
    return rate;
}

double PairProductionProcess::operator()(
    const Predictions &predictions, const Particle &firstParticle,
    const Particle &secondParticle) const noexcept {

    auto rate =
        predictions.bsmPairCxn(coll_, firstParticle.id(), secondParticle.id());

    if (firstParticle.id() != secondParticle.id() ||
        firstParticle.charge() != ECharge::neutral) {
        const auto distinctPairBr = [&firstParticle,
                                     &secondParticle](const auto &decays) {
            const auto &[d1, d2] = decays;
            return firstParticle.br(d1) * secondParticle.br(d2);
        };
        rate *= ranges::accumulate(
            ranges::views::cartesian_product(decay1_, decay2_), 0., std::plus{},
            distinctPairBr);
    } else {
        const auto identicalPairBr = [&firstParticle](const auto &decays) {
            const auto &[d1, d2] = decays;
            return (d1 == d2 ? 1 : 2) * firstParticle.br(d1) *
                   firstParticle.br(d2);
        };
        rate *= ranges::accumulate(unorderedDecayPairs_, 0., std::plus{},
                                   identicalPairBr);
    }

    return rate;
}

// ---------------- contributing particles ---------------- //
namespace {
auto containedParticles(const ParticleSet &c) {
    auto toId = [](const Particle &p) { return p.id(); };
    return c | ranges::views::transform(toId);
}
} // namespace

std::vector<std::string>
ChannelProcess::contributingParticles(const ParticleSet &cluster) {
    return containedParticles(cluster) | ranges::to<std::vector>;
}

std::vector<std::string>
ChainDecayProcess::contributingParticles(const ParticleSet &motherCluster,
                                         const ParticleSet &daughterCluster) {
    return ranges::views::concat(containedParticles(motherCluster),
                                 ranges::views::single(">"),
                                 containedParticles(daughterCluster)) |
           ranges::to<std::vector>;
}

std::vector<std::string> PairDecayProcess::contributingParticles(
    const ParticleSet &motherCluster, const ParticleSet &firstDaughterCluster,
    const ParticleSet &secondDaughterCluster) {
    return ranges::views::concat(containedParticles(motherCluster),
                                 ranges::views::single(">"),
                                 containedParticles(firstDaughterCluster),
                                 ranges::views::single("+"),
                                 containedParticles(secondDaughterCluster)) |
           ranges::to<std::vector>;
}

std::vector<std::string> PairProductionProcess::contributingParticles(
    const ParticleSet &firstParticleCluster,
    const ParticleSet &secondParticleCluster) {
    return ranges::views::concat(containedParticles(firstParticleCluster),
                                 ranges::views::single("+"),
                                 containedParticles(secondParticleCluster)) |
           ranges::to<std::vector>;
}
} // namespace Higgs::predictions
