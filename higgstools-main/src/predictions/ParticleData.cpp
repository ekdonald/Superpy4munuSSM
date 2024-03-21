#include "predictions/ParticleData.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "utilities/Algorithm.hpp"
#include "utilities/Logging.hpp"
#include <exception>
#include <limits>
#include <list>
#include <memory>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/functional/arithmetic.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <tuple>
#include <vector>
#include <math.h>
namespace Higgs::predictions {
void ParticleData::setCxn(Collider coll, Production p, double value) noexcept {
    if (p != Production::none) {
        cxns_[{coll, p}] = value;
    }
}

double ParticleData::cxn(Collider coll, Production p) const noexcept {
    return utilities::getWithDefault(cxns_, {coll, p}, 0.);
}

void ParticleData::updateBr(double &br, double value) {
    if (wtot_ > 0) {
        const auto newSumOfBrs = sumOfBrs_ - br + value;
        if (newSumOfBrs <= 1. + 100. * std::numeric_limits<double>::epsilon()) {
            br = value;
            sumOfBrs_ = newSumOfBrs;
        } else {
            throw InvalidInput("Cannot set BR, sum of BRs would exceed 1.");
        }
    } else {
        throw InvalidInput("Cannot set BR on zero width particle.");
    }
}

void ParticleData::setBr(Decay d, double value) {
    if (d != Decay::none) {
        updateBr(brs_[d], value);
    }
}

void ParticleData::setBr(std::pair<std::string, std::string> &&particleIds,
                         double value) {
    updateBr(bsmBrs_[std::move(particleIds)], value);
}

void ParticleData::updateBrFromDecayWidth(double &br, double width) {
    if (width < 0) {
        if (fabs(width) < 1.0e-10) {
          width = 0.0;
        } else {
          logger()->error("Invalid value {} for partial decay width.", width);
          throw InvalidInput("Invalid negative value for decay width");
        }
    }
    const auto oldWidth = totalWidth();
    const auto newWidth = oldWidth * (1 - br) + width;
    if (newWidth > 0) {
        const auto scaleBr = [scale = oldWidth / newWidth](auto &bsmBR) {
            bsmBR *= scale;
        };
        ranges::for_each(brs_ | ranges::views::values, scaleBr);
        ranges::for_each(bsmBrs_ | ranges::views::values, scaleBr);
        br = width / newWidth;
        sumOfBrs_ = (sumOfBrs_ * oldWidth + width) / newWidth;
    }
    setTotalWidth(newWidth);
}

void ParticleData::setDecayWidth(Decay d, double value) {
    if (d != Decay::none) {
        updateBrFromDecayWidth(brs_[d], value);
    }
}

void ParticleData::setDecayWidth(
    std::pair<std::string, std::string> &&particleIds, double value) {
    updateBrFromDecayWidth(bsmBrs_[std::move(particleIds)], value);
}

double ParticleData::br(Decay d) const noexcept {
    return utilities::getWithDefault(brs_, d, 0.);
}

double ParticleData::br(
    const std::pair<std::string, std::string> &particleIds) const noexcept {
    return utilities::getWithDefault(bsmBrs_, particleIds, 0.);
}

void ParticleData::setTotalWidth(double value) {
    if (value < 0) {
        logger()->error("Invalid value {} for total width", value);
        throw InvalidInput("Invalid value for total width");
    }
    wtot_ = value;
    if (wtot_ < 100. * std::numeric_limits<double>::epsilon()) {
        wtot_ = 0.;
        brs_.clear();
        sumOfBrs_ = 0.;
        bsmBrs_.clear();
    }
}

double ParticleData::totalWidth() const noexcept { return wtot_; }

std::optional<double> ParticleData::coupling(Coupling c) const noexcept {
    if (auto loc = coups_.find(c); loc != coups_.end())
        return loc->second;
    return std::nullopt;
}

void ParticleData::setCoupling(Coupling c, double value) { coups_[c] = value; }

void ParticleData::setChannelRate(Collider coll, Production p, Decay d,
                                  double value) {
    if (p != Production::none && d != Decay::none) {
        channelRates_[CollChannel{coll, p, d}] = value;
    }
}

void ParticleData::resetChannelRates() noexcept { channelRates_.clear(); }

double ParticleData::channelRateOr(Collider coll, Production p, Decay d,
                                   double defaultValue) const noexcept {
    return utilities::getWithDefault(channelRates_, CollChannel{coll, p, d},
                                     defaultValue);
}

double ParticleData::fullBrInv(double mass) const {
    auto b =
        br(Decay::directInv) + br(Decay::ZZ) * std::pow(constants::b_Z_inv, 2);
    if (context_ != nullptr) {
        auto lighterScalars =
            context_->particles() |
            ranges::views::filter([mass](const auto &p) {
                return p.mass() < mass && p.charge() == ECharge::neutral;
            });

        auto lighterInvBrs = lighterScalars |
                             ranges::views::transform([](const auto &p) {
                                 return std::pair{p.id(), p.br(Decay::inv)};
                             }) |
                             ranges::views::filter([](const auto &val) {
                                 return val.second > 0;
                             }) |
                             ranges::to<std::vector>;

        auto chainDecaysToInv =
            lighterInvBrs | ranges::views::transform([this](const auto &val) {
                const auto &[id, brC] = val;
                return br(chainDecayKey(ChainDecay::Z, id)) * brC *
                       constants::b_Z_inv;
            });

        auto pairDecaysToInv =
            ranges::views::cartesian_product(lighterInvBrs, lighterInvBrs) |
            ranges::views::transform([this](const auto &vals) {
                const auto &[v1, v2] = vals;
                const auto &[id1, br1] = v1;
                const auto &[id2, br2] = v2;
                const auto symFac = id1 == id2 ? 1. : 0.5;
                return br(twoBsmDecayKey(id1, id2)) * br1 * br2 * symFac;
            });

        b += ranges::accumulate(chainDecaysToInv, 0.) +
             ranges::accumulate(pairDecaysToInv, 0.);
    }
    return b;
}

void ParticleData::setContext(const Predictions &context,
                              ContextKey /* unused */) noexcept {
    context_ = &context;
}

void ParticleData::clearContext() noexcept { context_ = nullptr; }

} // namespace Higgs::predictions
