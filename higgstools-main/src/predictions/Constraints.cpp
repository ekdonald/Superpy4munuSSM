#include "predictions/Constraints.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "predictions/Clustering.hpp"
#include "predictions/JsonSupport.hpp" // IWYU pragma: keep
#include "utilities/Algorithm.hpp"
#include "utilities/Json.hpp"
#include "utilities/Logging.hpp"
#include <cmath>
#include <functional>
#include <map>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/algorithm/remove_if.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/unreachable_sentinel.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip_with.hpp>
#include <stdexcept>
#include <utility>

namespace Higgs::predictions {

ModelLikeness::ModelLikeness(ReferenceModel model,
                             const ChannelProcess &process)
    : ref_{getReference(model, 0)},
      combinedProc_{process}, subprocs_{process.subprocesses()} {}

ModelLikeness::ModelLikeness(const nlohmann::json &j, Collider collider)
    : ModelLikeness{
          utilities::readAs<ReferenceModel>(j, "modelLike"),
          predictions::readChannelProcess(j.at("process"), collider)} {}

ModelLikeness::ModelLikeness(const nlohmann::json &j,
                             const ChannelProcess &signalProcess)
    : ModelLikeness{utilities::readAs<ReferenceModel>(j, "modelLike"),
                    j.at("process").is_string() && j.at("process") == "signal"
                        ? signalProcess
                        : predictions::readChannelProcess(
                              j.at("process"), signalProcess.collider())} {}

ModelLikeness::ModelLikeness(const ModelLikeness &other)
    : ref_{other.ref_->clone()},
      combinedProc_{other.combinedProc_}, subprocs_{other.subprocs_} {}
ModelLikeness::ModelLikeness(ModelLikeness &&other) noexcept
    : ref_{std::move(other.ref_)},
      combinedProc_{std::move(other.combinedProc_)}, subprocs_{std::move(
                                                         other.subprocs_)} {}
ModelLikeness &ModelLikeness::operator=(ModelLikeness &&other) noexcept {
    ref_ = std::move(other.ref_);
    combinedProc_ = std::move(other.combinedProc_);
    subprocs_ = std::move(other.subprocs_);
    return *this;
}
ModelLikeness &ModelLikeness::operator=(const ModelLikeness &other) {
    auto tmp{other};
    std::swap(*this, tmp);
    return *this;
}

bool ModelLikeness::check(const ParticleSet &cluster, double refMass,
                          const Predictions & /* unused */) const noexcept {
    if (ref_) {
        // locking is needed because the mass of ref_ is modified in check.
        // Constness is ok though, no properties of ref_ are observable
        // outside this function.
        auto lock = std::scoped_lock{mut_};
        ref_->setMass(refMass);
        auto signalStrength = [cluster, this](const auto &proc) mutable {
            return Clustering::combinedRate(proc, cluster) / proc(*ref_);
        };

        auto ci = subprocs_ | ranges::views::transform(signalStrength);

        const auto refWeight = [this](const auto &proc) {
            return proc(*ref_) / combinedProc_(*ref_);
        };
        auto wi = subprocs_ | ranges::views::transform(refWeight);
        const auto weightedDeviation =
            [mu = signalStrength(combinedProc_)](double c, double w) {
                return w > 0 ? w * std::abs(c - mu) / mu : 0;
            };
        return ranges::max(ranges::views::zip_with(weightedDeviation, ci, wi)) <
               modelLikenessCut;
    } else {
        logger()->error("Ref is null in ModelLikeness::check");
        return false;
    }
}

TopDecayConsistency::TopDecayConsistency(
    const std::vector<Production> &modes) noexcept
    : brtBsm{Collider::LHC13,
             modes | ranges::views::transform([](Production p) {
                 return Channel{p, Decay::none};
             }) | ranges::to<std::vector>} {}

TopDecayConsistency::TopDecayConsistency(const nlohmann::json &j)
    : TopDecayConsistency{utilities::readAs<std::vector<Production>>(
          j, "topDecayConsistency")} {}

bool TopDecayConsistency::check(const ParticleSet &cluster,
                                double /* refMass */,
                                const Predictions &predictions) const noexcept {
    return std::abs(predictions.brTopWb() +
                    Clustering::combinedRate(brtBsm, cluster) - 1) <
           topDecayConsistencyCut;
}

CPValue::CPValue(CP cp) noexcept : cp_{cp} {}
CPValue::CPValue(const nlohmann::json &j)
    : CPValue{utilities::readAs<CP>(j, "CPValue")} {}
bool CPValue::check(const ParticleSet &cluster, double /* refMass */,
                    const Predictions & /* predictions */) const noexcept {
    return ranges::all_of(cluster,
                          [this](const Particle &p) { return p.cp() == cp_; });
}

MumuTautauRatio::MumuTautauRatio(CP cp) noexcept : cp_{cp} {}
MumuTautauRatio::MumuTautauRatio(const nlohmann::json &j)
    : MumuTautauRatio{utilities::readAs<CP>(j, "mumuTautauRatio")} {}

bool MumuTautauRatio::check(
    const ParticleSet &cluster, double /* refMass */,
    const Predictions & /* predictions */) const noexcept {
    const auto checkParticle = [this](const Particle &p) {
        if (p.mass() < 2 * constants::mTau) {
            return false;
        }
        auto widthRatio = p.br(Decay::mumu) / p.br(Decay::tautau);
        auto predRatio =
            std::pow(constants::mMu / constants::mTau, 2) /
            std::pow(std::sqrt(1 - std::pow(2 * constants::mTau / p.mass(), 2)),
                     cp_ == CP::even ? 3 : 1);
        return std::abs(widthRatio / predRatio - 1) < maxRelDeviation;
    };
    return ranges::all_of(cluster, checkParticle);
}

bool TopDominatedHgg::check(
    const ParticleSet &cluster, double /* refMass */,
    const Predictions & /* predictions */) const noexcept {
    constexpr auto checkParticle = [](const Particle &p) {
        const auto cxnRatio = [ref = SMHiggs(p.mass()),
                               &p](const Production &prod) {
            return p.cxn(Collider::LHC13, prod) /
                   ref.cxn(Collider::LHC13, prod);
        };
        // if effC are defined, test whether sigma_ggH(effTopC)/sigma_ggH_SM <
        // maxDeviation
        if (p.coupling(Coupling::effCPeTopYuk).has_value() &&
            p.coupling(Coupling::effCPoTopYuk).has_value()) {
            auto effTopC = std::complex<double>(
                p.coupling(Coupling::effCPeTopYuk).value(),
                p.coupling(Coupling::effCPoTopYuk).value());
            return std::abs(
                       std::sqrt(EffectiveCouplingCxns::ggH(
                                     Collider::LHC13, p.mass(), effTopC, 0.) /
                                 EffectiveCouplingCxns::ggH(Collider::LHC13,
                                                            p.mass(), 1., 0.)) -
                       std::sqrt(cxnRatio(Production::ggH))) < maxDeviation;
            // if effC are not defined, test whether sigma_ggH/sigma_ggH_SM is
            // compatible with sigma_ttH/sigma_ttH_SM
        } else {
            // for the CP-odd case, sigma_ggH/sigma_ggH_SM !=
            // sigma_ttH/sigma_ttH_SM even if ggH is top dominated.
            if (p.cp() == CP::odd) {
                auto ref = SMHiggs(p.mass());
                return std::abs(
                           std::sqrt(
                               p.cxn(Collider::LHC13, Production::Htt) /
                               (EffectiveCouplingRatios::ttHRatio(Collider::LHC13, p.mass(),
                                         std::complex<double>(0, 1)) *
                                ref.cxn(Collider::LHC13, Production::Htt))) -
                           std::sqrt(p.cxn(Collider::LHC13, Production::ggH) /
                                     EffectiveCouplingCxns::ggH(
                                         Collider::LHC13, p.mass(),
                                         std::complex<double>(0, 1), 0.))) <
                       maxDeviation;
            } else {
                return std::abs(std::sqrt(cxnRatio(Production::Htt)) -
                                std::sqrt(cxnRatio(Production::ggH))) <
                       maxDeviation;
            };
        };
    };
    return ranges::all_of(cluster, checkParticle);
}

std::vector<Constraint> readConstraints(const nlohmann::json &j,
                                        const ChannelProcess &signalProcess) {
    auto readOneConstraint = [&signalProcess](const auto &jj) -> Constraint {
        if (jj.contains("modelLike")) {
            return ModelLikeness{jj, signalProcess};
        } else if (jj.contains("topDecayConsistency")) {
            return TopDecayConsistency{jj};
        } else if (jj.contains("CPValue")) {
            return CPValue{jj};
        } else if (jj.contains("mumuTautauRatio")) {
            return MumuTautauRatio{jj};
        } else if (jj.contains("topDominatedHgg")) {
            return TopDominatedHgg{};
        }
        throw utilities::BadFieldRead("Unknown constraint in json.");
    };
    return j | ranges::views::transform(readOneConstraint) |
           ranges::to<std::vector>;
}

std::vector<Constraint> readConstraints(const nlohmann::json &j,
                                        Collider collider) {
    auto readOneConstraint = [collider](const auto &jj) -> Constraint {
        if (jj.contains("modelLike")) {
            return ModelLikeness{jj, collider};
        } else if (jj.contains("topDecayConsistency")) {
            return TopDecayConsistency{jj};
        } else if (jj.contains("CPValue")) {
            return CPValue{jj};
        } else if (jj.contains("mumuTautauRatio")) {
            return MumuTautauRatio{jj};
        } else if (jj.contains("topDominatedHgg")) {
            return TopDominatedHgg{};
        }
        throw utilities::BadFieldRead("Unknown constraint in json.");
    };
    return j | ranges::views::transform(readOneConstraint) |
           ranges::to<std::vector>;
}

bool checkConstraints(const std::vector<Constraint> &constraints,
                      const ParticleSet &cluster, double referenceMass,
                      const Predictions &context) {
    auto checkOneConstraint = [&cluster, referenceMass,
                               &context](const auto &con) {
        return con.check(cluster, referenceMass, context);
    };
    return ranges::all_of(constraints,
                          [&checkOneConstraint](const auto &anyCon) {
                              return std::visit(checkOneConstraint, anyCon);
                          });
}

void onlyValidClusters(std::vector<ParticleSet> &clusters,
                       const std::vector<Constraint> &constraints,
                       const Predictions &context) {
    auto failConstraints = [&context,
                            &constraints](const ParticleSet &cluster) {
        constexpr auto mass = [](const Particle &p) { return p.mass(); };
        return !checkConstraints(constraints, cluster,
                                 utilities::mean(cluster, mass), context);
    };
    clusters.erase(ranges::remove_if(clusters, failConstraints),
                   clusters.end());
}

ReferenceRate::ReferenceRate(ReferenceModel model, ChannelProcess process)
    : ref_{getReference(model, 0)}, process_{std::move(process)} {}
ReferenceRate::ReferenceRate(const ReferenceRate &other)
    : ref_{other.ref_->clone()}, process_{other.process_} {}
ReferenceRate::ReferenceRate(ReferenceRate &&other) noexcept
    : ref_{std::move(other.ref_)}, process_{std::move(other.process_)} {}
ReferenceRate &ReferenceRate::operator=(ReferenceRate &&other) noexcept {
    ref_ = std::move(other.ref_);
    process_ = std::move(other.process_);
    return *this;
}
ReferenceRate &ReferenceRate::operator=(const ReferenceRate &other) {
    auto tmp{other};
    std::swap(tmp, *this);
    return *this;
}

double ReferenceRate::operator()(double mass) const noexcept {
    if (ref_) {
        auto lock = std::scoped_lock{mut_};
        ref_->setMass(mass);
        return process_(*ref_);
    } else {
        logger()->error("Ref is null in ReferenceRate::operator()");
        return 0;
    }
}

std::optional<ReferenceRate>
ReferenceRate::read(const nlohmann::json &data,
                    const ChannelProcess &signalProcess) {
    if (data.contains("process") && data.contains("reference")) {
        auto ref = utilities::readAs<ReferenceModel>(data, "reference");
        if (data.at("process").is_string() && data["process"] == "signal") {
            return std::make_optional<ReferenceRate>(ref, signalProcess);
        } else {
            return std::make_optional<ReferenceRate>(
                ref, predictions::readChannelProcess(data["process"],
                                                     signalProcess.collider()));
        }
    }
    return std::nullopt;
}

} // namespace Higgs::predictions
