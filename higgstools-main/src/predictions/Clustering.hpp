#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/Constraints.hpp"
#include "predictions/UncertainMass.hpp"
#include "utilities/Algorithm.hpp"
#include "utilities/Concepts.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace Higgs::predictions {
class Predictions;
struct MassResolution;

namespace Clustering {

/**
 * @brief Returns a set of those particles from the model prediction that are
 * within the mass range.
 *
 * @param prediction the model prediction from which the BSMParticles are taken
 * @param massRange the mass range to consider
 * @param eagerness how to handle mass uncertainties:
 *   value  | effect
 * ---------|---------
 * eager    | true if any mass within the uncertainty is within range
 * cautious | only true if the entire uncertainty interval is within range
 * ignore   | ignore the mass uncertainty and check the central value
 * @param massUncMultiplier a factor to multiplies the theoretical mass
 * uncertainty
 * @return ParticleSet all particles from the model prediction that are in the
 * mass range
 */
ParticleSet particlesInMassRange(const Predictions &prediction,
                                 const std::pair<double, double> &massRange,
                                 MassUncEagerness eagerness,
                                 double massUncMultiplier = 1.);

//! A function object that filters out the relevant particles for a given
//! process. This class is used like a function through the
//! filterRelevantparticles object.
struct FilterRelevantParticles {
  public:
    /**
     * @brief Given candidate particles and a process, return only those
     * particles that are actually relevant for the process.
     *
     * @tparam Process the type of process to test
     * @tparam nParticleRoles the number of particle roles in the process
     * @param candidates the candidate particles for each role
     * @param process the process to test
     * @return std::array<ParticleSet, nParticleRoles> the relevant particles
     * for each role in the process
     */
    template <class Process, std::size_t nParticleRoles>
    std::array<ParticleSet, nParticleRoles>
    operator()(const std::array<ParticleSet, nParticleRoles> &candidates,
               const Process &process) const noexcept {
        auto positiveRate = [&process](const auto &particleAssignment) {
            return std::apply(process, particleAssignment) >
                   constants::minimumRate;
        };

        auto relevantAssignments =
            std::apply(ranges::views::cartesian_product, candidates) |
            ranges::views::filter(positiveRate);

        return filterRelevantParticlesImpl(
            candidates, std::move(relevantAssignments),
            std::make_index_sequence<nParticleRoles>{});
    }

  private:
    template <std::size_t I>
    static constexpr auto get = ranges::views::transform([](auto &&R) {
        return std::get<I>(R);
    });

    template <class ParticleRange>
    ParticleSet intersection(ParticleSet set,
                             ParticleRange &&range) const noexcept {
        constexpr auto particleId = [](const Particle &p) { return p.id(); };
        for (auto setIt = set.begin(); setIt != set.end();) {
            if (ranges::find(range, particleId(*setIt), particleId) ==
                end(range)) {
                setIt = set.erase(setIt);
            } else {
                ++setIt;
            }
        }
        return set;
    }

    template <std::size_t nParticleRoles, class RelevantAssignments,
              std::size_t... Is>
    auto filterRelevantParticlesImpl(
        const std::array<ParticleSet, nParticleRoles> &candidates,
        RelevantAssignments &&relevantAssignments,
        std::index_sequence<Is...>) const noexcept {

        return std::array{intersection(std::get<Is>(candidates),
                                       relevantAssignments | get<Is>)...};
    }
};

//! Function object to invoke FilterRelevantParticles::operator()()
constexpr inline FilterRelevantParticles filterRelevantParticles{};

//! Obtain the relevant particles for the given process from the model
//! prediction, keeping only those that are within the given massRanges.
template <class Process, std::size_t nParticleRoles>
std::array<ParticleSet, nParticleRoles> relevantParticles(
    const Process &process, const Predictions &prediction,
    const std::array<std::pair<double, double>, nParticleRoles> &massRanges,
    MassUncEagerness applicabilityEagerness, double massUncMultiplier = 1.) {

    auto getCandidateParticles = [&prediction, applicabilityEagerness,
                                  massUncMultiplier](auto &&...massRanges) {
        return std::array{particlesInMassRange(prediction, massRanges,
                                               applicabilityEagerness,
                                               massUncMultiplier)...};
    };
    auto candidates = std::apply(getCandidateParticles, massRanges);

    return filterRelevantParticles(candidates, process);
}

/**
 * @brief Generate all clusters out of the candidate particles, that are within
 * the given mass resolution of each other.
 *
 * @param candidates the candidate particles to cluster
 * @param resolution the mass resolution to use, the relative mass resolution is
 * taken relative to the mean mass of the cluster
 * @param clusterEagerness how to handle mass uncertainties of the involved
 * particles:
 *   value  | effect
 * ---------|---------
 * eager    | cluster as soon as uncertainties + resolution touch
 * cautious | cluster only if the mass uncertainties overlap entirely
 * ignore   | ignore the mass uncertainties when clustering
 * @return std::vector<ParticleSet> a vector of clusters
 */
std::vector<ParticleSet> cluster(const ParticleSet &candidates,
                                 const MassResolution &resolution,
                                 MassUncEagerness clusterEagerness);

//! Remove all clusters from the given vector of clusters that are subclusters
//! of any other cluster in the list.
void onlyMaximalClusters(std::vector<ParticleSet> &clusters);

//! A functor class that fully performs the clustering of the given candidates.
struct PerformClustering {
  public:
    /**
     * @brief Fully performs the clustering for the given candidate particles
     * for any number of particle roles.
     *
     * @tparam nParticleRoles the number of particle roles in the underlying
     * process
     * @param candidates the candidate particles for each role
     * @param massResolutions the mass resolution for each role
     * @param clusterEagerness the eagerness for the clustering, see cluster()
     * @param constraints the constraints on the clusters in each role
     * @param context the context for the checks of the constraints
     * @return std::array<std::vector<ParticleSet>, nParticleRoles>  a list of
     * the maximal clusters fulfilling all constraints for each role
     */
    template <std::size_t nParticleRoles>
    std::array<std::vector<ParticleSet>, nParticleRoles> operator()(
        const std::array<ParticleSet, nParticleRoles> &candidates,
        const std::array<MassResolution, nParticleRoles> &massResolutions,
        MassUncEagerness clusterEagerness,
        const std::array<std::vector<Constraint>, nParticleRoles> &constraints,
        const Predictions &context) const {
        return impl(candidates, massResolutions, clusterEagerness, constraints,
                    context, std::make_index_sequence<nParticleRoles>{});
    }

  private:
    template <std::size_t nParticleRoles, std::size_t... Is>
    std::array<std::vector<ParticleSet>, nParticleRoles>
    impl(const std::array<ParticleSet, nParticleRoles> &candidates,
         const std::array<MassResolution, nParticleRoles> &massResolutions,
         MassUncEagerness clusterEagerness,
         const std::array<std::vector<Constraint>, nParticleRoles> &constraints,
         const Predictions &context, std::index_sequence<Is...>) const {

        const auto maximalClustersWithConstraints =
            [clusterEagerness, &context](
                const ParticleSet &candidates, const MassResolution &resolution,
                const std::vector<Constraint> &constraints) {
                auto clusters =
                    cluster(candidates, resolution, clusterEagerness);
                // the order of the following two lines is crucial
                onlyValidClusters(clusters, constraints, context);
                onlyMaximalClusters(clusters);
                clusters.shrink_to_fit();
                return clusters;
            }; // LCOV_EXCL_LINE

        return std::array{maximalClustersWithConstraints(
            std::get<Is>(candidates), std::get<Is>(massResolutions),
            std::get<Is>(constraints))...};
    }
};

//! Function object to invoke PerformClustering::operator()()
constexpr inline PerformClustering performClusteringNew{};

/**
 * @brief The combined rate of all member particles in the cluster(s).
 * @tparam RateFunc function which takes one particle of each cluster (in
 * order) and returns the corresponding rate
 * @tparam Clusters any number of clusters (in the form of ranges::views)
 * @param rate function object to compute the rate from individual particles
 * @param clusters the clusters
 * @return double the combined rate summed over all possible combinations of
 * particles
 */
template <class RateFunc, class... Clusters>
REQUIRES(std::regular_invocable<RateFunc, ranges::range_value_t<Clusters>...>)
double combinedRate(const RateFunc &rate, Clusters &&...clusters) {
    return utilities::accumulatedCartesianProduct(rate, clusters...);
} // LCOV_EXCL_LINE

//! Convenience function returning a lambda that calls combinedRate for the
//! given rate function on any number of clusters.
template <class RateFunc> auto combinedRateFor(const RateFunc &rate) {
    return
        [&rate](auto &&...clusters) { return combinedRate(rate, clusters...); };
}

/**
 * @brief Weights the masses of the particles in the cluster(s) according to the
 * corresponding rates.
 *
 * @tparam RateFunc function which takes one particle of each cluster (in
 * order) and returns the corresponding rate
 * @tparam Clusters any number of clusters (in the form of ranges::views)
 * @param rate function object to compute the rate from individual particles
 * @param clusters the clusters
 * @return std::array<UncertainMass, sizeof...(Clusters)> the weighted masses of
 * each cluster (in order)
 */
template <class RateFunc, class... Clusters>
REQUIRES(std::regular_invocable<RateFunc, ranges::range_value_t<Clusters>...>)
std::array<UncertainMass, sizeof...(Clusters)>
rateWeightedMasses(const RateFunc &rate, Clusters &&...clusters) {
    static_assert(sizeof...(Clusters) <= 3,
                  "utilities::weightedMean is only implemented for up to three "
                  "ranges. It's straightforward to add more overloads.");
    constexpr auto uncertainMass = [](const Particle &p) {
        return UncertainMass{p.mass(), p.massUnc()};
    };
    if constexpr (sizeof...(Clusters) == 1) {
        return std::array{
            utilities::weightedMean(clusters..., rate, uncertainMass)};
    } else {
        return utilities::weightedMean(clusters..., rate, uncertainMass);
    }
}

//! Convenience function returning a lambda that calls rateWeightedMasses for
//! the given weight function on any number of clusters.
template <class WeightFunc> auto weightedMassesFor(const WeightFunc &weight) {
    return [&weight](auto &&...clusters) {
        return rateWeightedMasses(weight, clusters...);
    };
}

/**
 * @brief Weights the total widths of the particles in the cluster(s) according
 * to the corresponding rates.
 *
 * @tparam RateFunc function which takes one particle of each cluster (in
 * order) and returns the corresponding rate
 * @tparam Clusters any number of clusters (in the form of ranges::views)
 * @param rate function object to compute the rate from individual particles
 * @param clusters the clusters
 * @return std::array<double, sizeof...(Clusters)> the weighted total widths of
 * each cluster (in order)
 */
template <class RateFunc, class... Clusters>
REQUIRES(std::regular_invocable<RateFunc, ranges::range_value_t<Clusters>...>)
std::array<double, sizeof...(Clusters)>
rateWeightedWidths(const RateFunc &rate, Clusters &&...clusters) {
    static_assert(sizeof...(Clusters) <= 3,
                  "utilities::weightedMean is only implemented for up to three "
                  "ranges. It's straightforward to add more overloads.");
    if constexpr (sizeof...(Clusters) == 1) {
        return std::array{utilities::weightedMean(
            clusters..., rate, std::mem_fn(&Particle::totalWidth))};
    } else {
        return utilities::weightedMean(clusters..., rate,
                                       std::mem_fn(&Particle::totalWidth));
    }
}

} // namespace Clustering
} // namespace Higgs::predictions
