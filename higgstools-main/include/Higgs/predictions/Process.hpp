/**
 * @file Process.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Different kinds of processes.
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include <cstddef>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Higgs::predictions {
class Predictions;

//! A channel consists of a production and a decay mode.
using Channel = std::pair<Production, Decay>;

//! A process \f$ X X \to P \to Y Y\f$ for some BSM particle \f$P\f$
//! represented as a sum of channels that specify the production mode \f$ X X
//! \to P\f$ and the decay mode \f$ P \to Y Y\f$ at one given collider.
class ChannelProcess {
  public:
    //! Construct a channel process.
    ChannelProcess(Collider coll, std::vector<Channel> channels);

    //! default constructs an empty channel process.
    ChannelProcess() = default;

    //! The collider at which this process occurs.
    Collider collider() const noexcept;

    //! A string description of the process. If compactify is true the channels
    //! may be reordered for a more compact representation.
    std::string to_string(bool compactify = true) const noexcept;

    //! Returns the rate of the given particle for this process. The rate of a
    //! channel process is simply the sum of the Particle.channelRate#s for the
    //! contained sub-processes.
    double operator()(const Particle &particle) const noexcept;

    //! Returns the rate of the given particle for this process, using the given
    //! weights for each channel. If `weights.size() > size()` the extra
    //! weights are ignored, if `weights.size() < size()` the missing weights
    //! are assumed to be 1.
    double operator()(const Particle &particle,
                      const std::vector<double> &weights) const noexcept;

    //! Returns a vector containing the ids of the particles in the cluster.
    static std::vector<std::string>
    contributingParticles(const ParticleSet &cluster);

    //! Split this process into a list of subprocesses.
    std::vector<ChannelProcess> subprocesses() const noexcept;

    //! Merge the channels of the other process into this one.
    ChannelProcess &operator+=(const ChannelProcess &other);

    //! check if processes are identical (sub-channel order matters)
    bool operator==(const ChannelProcess &other) const noexcept;

    //! Merge two ChannelProcess%es together. Return lhs if the processes can't
    //! be merged.
    friend ChannelProcess operator+(ChannelProcess lhs,
                                    const ChannelProcess &rhs);

    //! Return the number of channels in this process.
    std::size_t size() const noexcept;

  private:
    Collider coll_ = Collider::LHC13;
    std::vector<Channel> channels_;
};

//! A process \f$ XX \to P_1 \to P_2 Y \to ZZ \f$ involving two BSM particles
//! \f$ P_{1,2}\f$ (referred to as mother and daughter particle) with a sum over
//! the production modes \f$ XX \f$ and final states \f$ ZZ \f$.
class ChainDecayProcess {
  public:
    //! Construct a chain decay process from the given chain type, as well as
    //! production and decay modes at a given collider.
    ChainDecayProcess(ChainDecay chain, std::vector<Production> production,
                      std::vector<Decay> decay, Collider coll);

    //! default constructor for an empty process
    ChainDecayProcess() = default;

    //! The kind of chain decay happening
    ChainDecay chain() const noexcept;
    //! The collider at which this process occurs.
    Collider collider() const noexcept;

    //! A string description of the process
    std::string to_string() const noexcept;

    //! Returns the rate of this process for the given particles. Considers all
    //! possible combinations of initial and final states.
    double operator()(const Particle &mother,
                      const Particle &daughter) const noexcept;

    //! Returns the rate of this process for the given particles, using the
    //! given weights for each production mode. If `productionWeights.size() !=
    //! this.productionSize()`, any extra weights are ignored or the missing
    //! weights are assumed to be 1.
    double
    operator()(const Particle &mother, const Particle &daughter,
               const std::vector<double> &productionWeights) const noexcept;

    /// A vector containing the ids of the particles in the clusters.
    /// Formatting for a ChainDecayProcess involving a mother and daughter
    /// cluster. The particles in the mother cluster are separated from
    /// those in the daughter cluster by an entry of `">"`, i.e.
    /// @verbatim embed:rst:leading-slashes
    /// .. code-block:: c++
    ///
    ///       {/* ids of the mothers */ , ">", /* ids of the daughters */}
    /// @endverbatim
    static std::vector<std::string>
    contributingParticles(const ParticleSet &motherCluster,
                          const ParticleSet &daughterCluster);

    //! Returns the number of production modes in this process.
    std::size_t productionSize() const noexcept;

    //! Returns the number of decay modes in this process.
    std::size_t decaySize() const noexcept;

  private:
    ChainDecay chain_;
    std::vector<Production> production_;
    std::vector<Decay> decay_;
    Collider coll_ = Collider::LHC13;
};

//! A process of the form \f$ XX \to P_1 \to P_2 P_3, P_2 \to YY, P_3 \to ZZ
//! \f$, where \f$ P_{1,2,3} \f$ are BSM particles (referred to as mother, first
//! daughter, and second daughter particle, respectively) with a sum over the SM
//! initial state \f$XX\f$ and the SM final states \f$YY, ZZ\f$. For
//! non-resonant pair production without the on-shell BSM resonance, see
//! PairProductionProcess.
class PairDecayProcess {
  public:
    //! Construct a pair decay process where a resonance is produced in the
    //! production modes at the given collider and decays into two daughter
    //! particles that then decay in the modes firstDecay and secondDecay.
    PairDecayProcess(std::vector<Production> production,
                     std::vector<Decay> firstDecay,
                     std::vector<Decay> secondDecay, Collider coll);

    //! default constructor for an empty process
    PairDecayProcess() = default;

    //! The collider at which this process occurs.
    Collider collider() const noexcept;

    //! A string description of the process.
    std::string to_string() const noexcept;

    //! Returns the rate of this process for the given particles. Considers all
    //! possible combinations of mother-particle production modes and
    //! daughter-particle decay modes. Does **not** account for permutations of
    //! the daughter particles (if possible), these are handled at the call site
    //! (typically through a cartesian product). However the process does
    //! include the required symmetry factors for different decay modes of
    //! identical daughter particles. See PairProductionProcess::operator()()
    //! for details on the combinatorics.
    double operator()(const Particle &mother, const Particle &firstDaughter,
                      const Particle &secondDaughter) const noexcept;

    /// A vector containing the ids of the particles in the clusters. Formatting
    /// for a PairDecayProcess involving a mother and two potentially distinct
    /// daughter clusters. Contains one entry of each `">"` and `"+"` to
    /// separate the particles in the different roles, i.e.
    /// @verbatim embed:rst:leading-slashes
    /// .. code-block:: c++
    ///
    ///       {/* ids of the mothers */ , ">", /* ids of the first daughters */,
    ///       "+", /* ids of the second daughters*/}
    /// @endverbatim
    static std::vector<std::string>
    contributingParticles(const ParticleSet &motherCluster,
                          const ParticleSet &firstDaughterCluster,
                          const ParticleSet &secondDaughterCluster);

  private:
    std::vector<Production> production_;
    std::vector<Decay> decay1_;
    std::vector<Decay> decay2_;
    std::set<std::pair<Decay, Decay>> unorderedDecayPairs_;
    Collider coll_ = Collider::LHC13;
};

//! A non-resonant pair production process of the form \f$ XX \to P_1 P_1, P_1
//! \to YY, P_2 \to ZZ \f$, where \f$ P_{1,2} \f$ are BSM particles with a sum
//! over the the SM final states \f$YY, ZZ\f$. For resonant pair production with
//! an on-shell BSM resonance, see PairDecayProcess.
class PairProductionProcess {
  public:
    //! Construct a non-resonant pair production process at the given collider,
    //! where the two produced particles decay into the given decay modes.
    PairProductionProcess(std::vector<Decay> firstDecay,
                          std::vector<Decay> secondDecay, Collider coll);

    //! default constructor for an empty process
    PairProductionProcess() = default;

    //! The collider at which this process occurs.
    Collider collider() const noexcept;

    //! A string description of the process.
    std::string to_string() const noexcept;

    /**
     * Returns the rate of this process for the given particles. If the process
     * cannot occur for the given particles, returns 0.
     *
     * Given two sets of final states \f$A=\{a_1,a_2,\ldots\}\f$ and
     * \f$B=\{b_1,b_2,\ldots\}\f$ the combined decay rate of two particles
     * \f$h_{i,j}\f$ into any combination of those final states is given by
     *
     * \f[
     * \mathrm{Br}(h_i h_j \to A B) =
     * \begin{cases}
     *   \sum_{a\in A}\sum_{b \in B}
     *     \mathrm{Br}(h_i \to a)\mathrm{Br}(h_j \to b)
     *   &i\neq j,\\
     *   \sum_{\{a,b\}\forall a\in A,b\in B} S(\{a,b\})
     *     \mathrm{Br}(h_i \to a)\mathrm{Br}(h_i \to b)
     *   &i=j.
     * \end{cases}
     * \f]
     *
     * The sum in the second case runs over all *unique unordered* pairs
     * \f$\{a,b\}\f$ and the symmetry factor \f$S\f$ is
     *
     * \f[
     * S(\{a,b\})=
     * \begin{cases}
     *   1 & a=b,\\
     *   2 & a\neq b.
     * \end{cases}
     * \f]
     *
     * For example, if \f$A=\{b b, \tau\tau\}\f$ and \f$B=\{bb,\gamma\gamma\}\f$
     * the result would be
     *
     * \f[
     * \mathrm{Br}(h_i h_j \to A B) =
     * \begin{cases}
     *   \mathrm{Br}^i_{bb}\mathrm{Br}^j_{bb}
     *   + \mathrm{Br}^i_{bb}\mathrm{Br}^j_{\gamma\gamma}
     *   + \mathrm{Br}^i_{\tau\tau}\mathrm{Br}^j_{bb}
     *   + \mathrm{Br}^i_{\tau\tau}\mathrm{Br}^j_{\gamma\gamma}
     *   & i\neq j\,\\
     *   {(\mathrm{Br}^i_{bb})}^2
     *   + 2\,\mathrm{Br}^i_{bb}\mathrm{Br}^i_{\gamma\gamma}
     *   + 2\,\mathrm{Br}^i_{\tau\tau}\mathrm{Br}^i_{bb}
     *   + 2\,\mathrm{Br}^i_{\tau\tau}\mathrm{Br}^i_{\gamma\gamma}
     *   & i = j,
     * \end{cases}
     * \f]
     *
     * where \f$\mathrm{Br}^i_{a}=\mathrm{Br}(h_i\to a)\f$.
     * Permutations of \f$h_i\f$ and \f$h_j\f$ when \f$i\neq j\f$ are *not*
     * included here, but are instead accounted for by sums over the
     * corresponding clusters \f$C_1, C_2\f$ at the call sites of this function:
     *
     * \f[
     * \mathrm{Br}(C_1 C_2 \to A B) =
     * \sum_{h_i \in C_1} \sum_{h_j \in C_2}
     * \mathrm{Br}(h_i h_j \to A B).
     * \f]
     *
     * In all of the above charged particle-antiparticle pairs need to be
     * treated as distinct, though they are treated as the same particle
     * everywhere else in the code. To address this, non-resonant pair
     * production of equal id charged particles is always assumed to be neutral
     * overall (i.e. `H+H-` and never `H+H+`).
     *
     * @param predictions the full model predictions, needed to access pair
     * production cross sections of distinct particles
     * @param firstParticle the first produced particle
     * @param secondParticle the second produced particle
     * @return double the rate of this process
     */
    double operator()(const Predictions &predictions,
                      const Particle &firstParticle,
                      const Particle &secondParticle) const noexcept;

    /// A vector containing the ids of the particles in the clusters. Formatting
    /// for a PairProductionProcess involving two produced clusters. Contains
    /// one `"+"` to separate the particles in the different roles, i.e.
    /// @verbatim embed:rst:leading-slashes
    /// .. code-block:: c++
    ///
    ///       {/* ids of first particles */ , "+",
    ///        /* ids of second particles */}
    /// @endverbatim
    static std::vector<std::string>
    contributingParticles(const ParticleSet &firstParticleCluster,
                          const ParticleSet &secondParticleCluster);

  private:
    std::vector<Decay> decay1_;
    std::vector<Decay> decay2_;
    std::set<std::pair<Decay, Decay>> unorderedDecayPairs_;
    Collider coll_ = Collider::LHC13;
};

} // namespace Higgs::predictions
