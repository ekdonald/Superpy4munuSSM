/**
 * @file Acceptance.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Non-trivial acceptance times efficiency factors.
 *
 * @copyright Copyright 2021 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/predictions/Channels.hpp"
#include "utilities/JsonFwd.hpp"
#include "utilities/LinearInterpolator.hpp"
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace Higgs::predictions {
class Particle;
}

namespace Higgs::bounds {
//! A constant acceptance.
class ConstantAcceptance {
  public:
    explicit constexpr ConstantAcceptance(double value = 1) : value_{value} {}

    /**
     * @brief Read a new Constant Acceptance object from json
     *
     * @param j the json to read from. The representation is simply
     * ```json
     * {
     *      "constantAcceptance": 1.5
     * }
     * ```
     * where `1.5` is the value.
     */
    explicit ConstantAcceptance(const nlohmann::json &j);

    //! Return the constant acceptance value.
    constexpr double operator()(const predictions::Particle & /* p */) const {
        return value_;
    }

  private:
    double value_;
};

//! An acceptance that depends on the mass of the particle.
class MassDepAcceptance {
  public:
    /**
     * @brief Construct a new mass dependent acceptance object.
     *
     * @throws std::out_of_range if the two grids don't have the same size
     * @param massGrid the mass values corresponding to the acceptance values
     * @param values the acceptance values
     */
    explicit MassDepAcceptance(const std::vector<double> &values,
                               std::vector<double> massGrid);
    /**
     * @brief Read a new mass dependent acceptance from json.
     *
     * @param j the json to read from. The json representation is e.g.
     * ```json
     * {
     *      "massDepAcceptance": [0.6,0.7,0.8],
     *      "massGrid": [100,200,500]
     * }
     * ```
     * where the array in the "massDepAcceptance" field contains the acceptance
     * values.
     */
    explicit MassDepAcceptance(const nlohmann::json &j);

    //! Evaluate the acceptance at the mass of the particle p.
    double operator()(const predictions::Particle &p) const;

  private:
    utilities::LinearInterpolator<1> interp_;
};

class CouplingDepAcceptance;

//! Any kind of Acceptance.
using Acceptance =
    std::variant<ConstantAcceptance, MassDepAcceptance, CouplingDepAcceptance>;

class CouplingDepAcceptance {
  public:
    //! A coupling power may be positive, negative or zero.
    using Power = int;
    //! Defines a Monomial as a product of the contained couplings raised to
    //! their respective powers.
    using Monomial = std::unordered_map<predictions::Coupling, Power>;
    //! Each term contains a monomial and and an acceptance that defines the
    //! coefficient that multiplies this polynomial.
    using Term = std::pair<Monomial, Acceptance>;
    //! A polynomial is the sum of several terms. Since we allow negative
    //! powers, this is technically a Laurent Polynomial.
    using Polynomial = std::vector<Term>;

    /**
     * @brief Construct a new coupling dependent acceptance.
     *
     * The total acceptance is parametrized as a fraction of two Polynomials of
     * the couplings. The acceptance coefficients in each polynomial may be any
     * kind of Acceptance and will be evaluated on the particle. Usually they
     * will be MassDepAcceptance%s, but they can also be ConstantAcceptance%s or
     * even CouplingDepAcceptance%s themselves.
     *
     * @param numerator the Polynomial that defines the coupling dependence and
     * acceptance coefficients of the numerator
     * @param denominator the Polynomial that defines the coupling dependence
     * and acceptance coefficients of the denominator, if this is empty the
     * denominator is assumed to be 1.
     */
    explicit CouplingDepAcceptance(Polynomial numerator,
                                   Polynomial denominator = {});

    /**
     * @brief Read a new coupling dependent acceptance from json.
     *
     * @param j the json to read from. The json representation is e.g.
     * ```json
     * {
     *      "couplingDepAcceptance": [
     *        [
     *          {"effCPeTopYuk": 2, "effCPoTopYuk": 2},
     *          {"constantAcceptance": 3.1}
     *        ],
     *        [
     *          {"effCPeTopYuk": 4},
     *          {
     *            "massDepAcceptance": [0.3,0.7,0.9],
     *            "massGrid": [100,200,500]
     *          }
     *        ]
     *      ],
     *      "denominator": [
     *        [{"effCPeTopYuk": 4}, {"constantAcceptance": 1.1}]
     *      ]
     * }
     * ```
     * where the array in the "couplingDepAcceptance" field defines the
     * numerator polynomial, and the `"denominator"` entry is optional (and
     * default to a denominator of `1.` if absent). Term of the polynomials is
     * represented by a two-element array, where the first element is a
     * dictionary that specifies the coupling powers and the second element is
     * any kind of acceptance in its json representation (this can even by
     * another CouplingDepAcceptance).
     */
    explicit CouplingDepAcceptance(const nlohmann::json &j);

    //! Evaluate the acceptance for the given particle. For both numerator and
    //! denominator polynomial, this evaluates all of the acceptance
    //! coefficients for the particle, multiplies them with their respective
    //! powers of the particle's couplings, and sums the resulting
    //! contributions. Then the numerator is divided by the denominator,
    //! obviously the denominator must never be zero. The only such case that is
    //! guaranteed to be well behaved is if both numerator and denominator are
    //! exactly zero in which case the result is zero.
    double operator()(const predictions::Particle &p) const;

  private:
    static double evaluateCouplingMonomial(const Monomial &monomial,
                                           const predictions::Particle &p);

    Polynomial numerator_;
    Polynomial denominator_;
};

//! A class storing the acceptances corresponding to a process.
class Acceptances {
  public:
    /**
     * @brief Construct a new Acceptances object
     * @throws std::out_of_range if `acceptances.size()!=nChannels`
     * @param acceptances the vector of acceptances to use
     * @param nChannels the number of channels in the corresponding process
     */
    explicit Acceptances(std::vector<Acceptance> acceptances,
                         std::size_t nChannels);
    /**
     * @brief Read acceptances from json.
     * @param j the json to read from. The representation is an array of objects
     * where each object can be read using readAcceptance.
     * @param nChannels the number of channels in the corresponding process
     */
    explicit Acceptances(const nlohmann::json &j, std::size_t nChannels);

    //! Evaluate all acceptances for the given particles.
    std::vector<double> operator()(const predictions::Particle &p) const;

  private:
    std::vector<Acceptance> acceptances_;
};

/**
 * @brief Read an Acceptance from json. The type of acceptance to read is
 * decided based on the presence of a key on the object, that matches the class
 * name (except for the first letter being lowercase). See the json constructors
 * of the individual types for the specific representations.
 *
 * @throws utilities::badFieldRead if the json object cannot be parsed into a
 * known acceptance
 * @param j the json to read from.
 * @return Acceptance an acceptance
 */
Acceptance readAcceptance(const nlohmann::json &j);

} // namespace Higgs::bounds
