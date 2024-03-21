#include "predictions/Clustering.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "predictions/UncertainMass.hpp"
#include "range/v3/iterator/basic_iterator.hpp"
#include "range/v3/iterator/diffmax_t.hpp"
#include "range/v3/iterator/operations.hpp"
#include "range/v3/iterator/unreachable_sentinel.hpp"
#include "range/v3/range/conversion.hpp"
#include "range/v3/range/primitives.hpp"
#include "range/v3/utility/get.hpp"
#include "range/v3/view/enumerate.hpp"
#include "range/v3/view/view.hpp"
#include "utilities/Algorithm.hpp"
#include <list>
#include <new>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/algorithm/min.hpp>
#include <range/v3/algorithm/remove_if.hpp>
#include <range/v3/algorithm/set_algorithm.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

namespace Higgs::predictions::Clustering {

ParticleSet particlesInMassRange(const Predictions &prediction,
                                 const std::pair<double, double> &massRange,
                                 MassUncEagerness eagerness,
                                 double massUncMultiplier) {
    constexpr auto toParticleRef = [](const Particle &p) {
        return std::cref(p);
    };
    auto allParticles =
        prediction.particles() | ranges::views::transform(toParticleRef);

    const auto inMassRange = [&massRange, eagerness](const UncertainMass &m) {
        return m.liesWithin(massRange, eagerness);
    };
    const auto uncParticleMasses = [massUncMultiplier](const Particle &p) {
        return UncertainMass{p.mass(), massUncMultiplier * p.massUnc()};
    };
    return allParticles |
           ranges::views::filter(inMassRange, uncParticleMasses) |
           ranges::to<ParticleSet>;
}

namespace {
/**
 * @brief Generates the powerset (the set of all subsets) on a range.
 *
 * For a given range (treated as a set) of size `n` the powerset has \f$ 2^n
 * \f$ elements (it always includes the empty set). E.g. for `{1, 2, 3}` the
 * powerset will contain
 * `{ {}, {1}, {2}, {3}, {1,2}, {1,3}, {2,3}, {1,2,3} }`. The subsets won't
 * necessarily be in this order *but the empty set will always be the first
 * subset*.
 *
 * The implementation uses the isomorphism between the powerset and the
 * binary representations of the numbers from 0 to \f$ 2^n - 1\f$ (see
 * https://en.wikipedia.org/wiki/Power_set#Representing_subsets_as_functions)
 *
 * @tparam R a ranges::forward_range and ranges::viewable_range
 * @param r the given range
 */
auto powerSet(const ParticleSet &particles) {
    const auto powerSetSize = 1UL << ranges::distance(particles);
    const auto subsetAtIndex = [&particles](std::size_t setI) {
        const auto elementPresent = [setI](const auto &elem) -> bool {
            const auto &[elemI, value] = elem;
            return (1ul << elemI) & setI;
        };
        return ranges::views::enumerate(particles) |
               ranges::views::filter(elementPresent) | ranges::views::values;
    };

    return ranges::views::ints(0UL, powerSetSize) |
           ranges::views::transform(subsetAtIndex);
}

template <class CandidateCluster>
bool closeInMass(CandidateCluster &&candidates, MassResolution resolution,
                 MassUncEagerness clusterEagerness) {
    if (ranges::empty(candidates)) {
        return false; // we don't want empty clusters
    }
    if (ranges::distance(candidates) == 1) {
        return true; // a single particle is always a valid cluster
    }

    constexpr auto getMassMinusUnc = ranges::views::transform(
        [](const Particle &p) { return p.mass() - p.massUnc(); });
    constexpr auto getMassPlusUnc = ranges::views::transform(
        [](const Particle &p) { return p.mass() + p.massUnc(); });
    constexpr auto getMass =
        ranges::views::transform([](const Particle &p) { return p.mass(); });

    auto maxMass = 0.;
    auto minMass = 0.;
    switch (clusterEagerness) {
    case MassUncEagerness::eager:
        maxMass = ranges::max(candidates | getMassMinusUnc);
        minMass = ranges::min(candidates | getMassPlusUnc);
        break;
    case MassUncEagerness::cautious:
        maxMass = ranges::max(candidates | getMassPlusUnc);
        minMass = ranges::min(candidates | getMassMinusUnc);
        break;
    case MassUncEagerness::ignore:
        maxMass = ranges::max(candidates | getMass);
        minMass = ranges::min(candidates | getMass);
        break;
    }
    const auto avgMass = utilities::mean(candidates | getMass);
    return maxMass - minMass <=
           avgMass * resolution.relative + resolution.absolute;
}
} // namespace

std::vector<ParticleSet> cluster(const ParticleSet &candidates,
                                 const MassResolution &resolution,
                                 MassUncEagerness clusterEagerness) {
    const auto validMass = [&resolution,
                            clusterEagerness](auto &&clusterCandidate) {
        return closeInMass(clusterCandidate, resolution, clusterEagerness);
    };

    return powerSet(candidates) | ranges::views::filter(validMass) |
           ranges::to<std::vector<ParticleSet>>;
}

void onlyMaximalClusters(std::vector<ParticleSet> &clusters) {
    const auto isSubsetOfAnother = [&clusters](const ParticleSet &cluster) {
        auto clusterIsSubsetOf = [&cluster](const ParticleSet &otherCluster) {
            return cluster.size() < otherCluster.size() &&
                   ranges::includes(otherCluster, cluster, cluster.key_comp());
        };
        return ranges::any_of(clusters, clusterIsSubsetOf);
    };

    clusters.erase(ranges::remove_if(clusters, isSubsetOfAnother),
                   clusters.end());
}

} // namespace Higgs::predictions::Clustering
