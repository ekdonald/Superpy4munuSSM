#pragma once

#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/Process.hpp"
#include "utilities/JsonFwd.hpp"
#include <memory>
#include <mutex>
#include <optional>
#include <variant>
#include <vector>

namespace Higgs::predictions {
class Predictions;

//! ModelLikeness constraint that checks that the signal composition is
//! similar to a reference model.
class ModelLikeness {
  public:
    static constexpr double modelLikenessCut = 2e-2;

    //! Checks that the given cluster fulfills the model-likeness constraint at
    //! the given reference mass.
    bool check(const ParticleSet &cluster, double refMass,
               const Predictions & /* unused */) const noexcept;

    static std::vector<ModelLikeness> read(const nlohmann::json &j,
                                           const ChannelProcess &signalProcess);

    static std::vector<ModelLikeness> read(const nlohmann::json &analysis,
                                           Collider collider);

    //! Constructs a new ModelLikeness object that imposes similarity to the
    //! given model for the given process.
    explicit ModelLikeness(ReferenceModel model, const ChannelProcess &process);

    /**
     * @brief Read a ModelLikeness constraint from json.
     *
     * The json layout is e.g.:
     * ```json
     * {
     *    "modelLike": "SMHiggs",
     *    "process": {"...":"..."}
     * }
     * ```
     * where `"modelLike"` can be any ReferenceModel and `"process"` describes a
     * ChannelProcess in the appropriate json format.
     *
     * @param j the json to read from
     * @param collider the collider to use for the process construction
     */
    explicit ModelLikeness(const nlohmann::json &j, Collider collider);

    //! Same as ModelLikeness(const nlohmann::json , Collider) except that if
    //! `"process": "signal"` the signal process is used.
    explicit ModelLikeness(const nlohmann::json &j,
                           const ChannelProcess &signalProcess);

    //! @{
    //! value semantics
    ModelLikeness &operator=(ModelLikeness &&other) noexcept;
    ModelLikeness(ModelLikeness &&other) noexcept;
    ModelLikeness &operator=(const ModelLikeness &other);
    ModelLikeness(const ModelLikeness &other);
    ~ModelLikeness() = default;
    //! @}
  private:
    std::unique_ptr<Particle> ref_;
    ChannelProcess combinedProc_;
    std::vector<ChannelProcess> subprocs_;
    mutable std::mutex mut_ = {};
};

//! Constraint that checks if the top-quark decays into signal particles and Wb
//! sum so approximately one.
class TopDecayConsistency {
  public:
    static constexpr double topDecayConsistencyCut = 2e-2;

    explicit TopDecayConsistency(const std::vector<Production> &modes) noexcept;

    //! Json constructor. The json layout for this class is
    //! ```json
    //! {"topDecayConsistency": ["brtHpb", "brtHq"]}
    //! ```
    //! where the Production modes in the array indicate which top BSM decays
    //! should be considered.
    explicit TopDecayConsistency(const nlohmann::json &j);

    //! Check that the top-quark BRs into `tb` plus all considered decays
    //! involving particles in the cluster sum to 1.
    bool check(const ParticleSet &cluster, double /* refMass */,
               const Predictions &predictions) const noexcept;

  private:
    ChannelProcess brtBsm;
};

//! Constraint on the CP-character of particles.
class CPValue {
  public:
    explicit CPValue(CP cp) noexcept;
    explicit CPValue(const nlohmann::json &j);

    //! True exactly if all particles in the cluster are of the specified
    //! CP-character.
    bool check(const ParticleSet &cluster, double /* refMass */,
               const Predictions & /* predictions */) const noexcept;

  private:
    CP cp_;
};

//! Constraint that checks that the ratio between BR(mumu) and BR(tautau)
//! fulfills the given assumptions. At higher masses this is simply \f$
//! \mathrm{BR}(\tau\tau)/\mathrm{BR}(\mu\mu)\approx \frac{m^2_\tau}{m^2_\mu}\f$
//! independent of the assumed CP, while at BSM particle masses close to
//! \f$2m_\tau\f$ the CP-differences become important and the correct tree-level
//! relation for each case is used.
class MumuTautauRatio {
  public:
    //! chosen such that the SM-Higgs always fulfills this for CP-even
    static constexpr double maxRelDeviation = 5.5e-2;

    explicit MumuTautauRatio(CP cp) noexcept;

    //! Json constructor. The json layout for this class is
    //! ```json
    //! {"MumuTautauRatio": "even"}
    //! ```
    //! where the value specifies the CP-assumption.
    explicit MumuTautauRatio(const nlohmann::json &j);

    //! Check if all particles in the cluster fulfill the assumption on
    //! `BR(mumu)/BR(tautau)`.
    bool check(const ParticleSet &cluster, double /* refMass */,
               const Predictions & /* predictions */) const noexcept;

  private:
    bool checkParticle(const Particle &p) noexcept;
    CP cp_;
};

class TopDominatedHgg {
  public:
    static constexpr double maxDeviation = 5e-2;

    //! True if the coupling of every particle in the cluster to two gluons is
    //! approximately equal to the coupling to top quarks.
    bool check(const ParticleSet &cluster, double /* refMass */,
               const Predictions & /* predictions */) const noexcept;

  private:
    static bool checkParticle(const Particle &p) noexcept;
};

using Constraint = std::variant<TopDecayConsistency, ModelLikeness, CPValue,
                                MumuTautauRatio, TopDominatedHgg>;

//! @{
//! Read an array of constraints from Json. Constraints are identified according
//! to the json layout given in their documentation.
std::vector<Constraint> readConstraints(const nlohmann::json &j,
                                        const ChannelProcess &signalProcess);
std::vector<Constraint> readConstraints(const nlohmann::json &j, Collider coll);
//! @}

//! Check the cluster against all constraints in the vector. The reference mass
//! and context may or may not be used, depending on the specific Constraint%s.
bool checkConstraints(const std::vector<Constraint> &constraints,
                      const ParticleSet &cluster, double referenceMass,
                      const Predictions &context);

//! Remove all clusters from the given vector of clusters that do not fulfill
//! the specified constraints.
void onlyValidClusters(std::vector<ParticleSet> &clusters,
                       const std::vector<Constraint> &constraints,
                       const Predictions &context);

class ReferenceRate {
  public:
    //! return the reference rate for the given mass
    double operator()(double mass) const noexcept;
    //! Construct a rate reference using the given reference model and
    //! process.
    ReferenceRate(ReferenceModel model, ChannelProcess process);

    static std::optional<ReferenceRate>
    read(const nlohmann::json &data, const ChannelProcess &signalProcess);

    ReferenceRate &operator=(ReferenceRate &&other) noexcept;
    ReferenceRate(ReferenceRate &&other) noexcept;
    ReferenceRate &operator=(const ReferenceRate &other);
    ReferenceRate(const ReferenceRate &other);
    ~ReferenceRate() = default;

  private:
    std::unique_ptr<predictions::Particle> ref_;
    predictions::ChannelProcess process_;
    mutable std::mutex mut_ = {};
};
} // namespace Higgs::predictions
