/**
 * @file Likelihood.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Exclusion likelihoods and CLs computation
 *
 * @copyright Copyright 2019 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "utilities/Concepts.hpp"
#include "utilities/Logging.hpp"
#include "utilities/RootFinding.hpp"
#include <exception>
#include <memory>
#include <utility>

//! Functions to obtain limits from exclusion likelihood information. The
//! calculation of the CLs from the likelihood relies on the fact that in the
//! conventions of 1007.1727 `qObs`\f$=q_\mu\f$ and in the asymptotic limit
//! `qExp`\f$=\mu^2/\sigma^2\f$.
namespace Higgs::bounds::Likelihood {

/**
 * Expected CLs from likelihood information.
 *
 * Calculated as \f$\mathrm{CL}_{s+b}/\mathrm{CL}_b\f$ using CLsb from eq. (66)
 * of 1007.1727 with \f$\tilde{q}_mu=\f$`qExp` and CLb=0.5.
 */
double expCLs(double qExp);

/**
 * Observed CLs from likelihood information.
 *
 * Calculated using CLsb from eq. (66) of 1007.1727 and CLb from eq. (65) with
 * \f$\mu'=0\f$. Expanded to remove numerical issues when `qExp->0` or
 * `qObs>>qExp`.
 */
double obsCLs(double qExp, double qObs);

//! The CLs cut used for limit setting
inline constexpr auto limCLs = 0.05;

//! expected and observed limit ratios reconstructed from likelihood profiles
struct LikelihoodRatios {
    double expRatio; //!< predicted rate over expected \f$\CL{95}\f$ limit
    double obsRatio; //!< predicted rate over observed \f$\CL{95}\f$ limit
};

/**
 * Find the obsRatio and expRatio from given likelihood profiles.
 *
 * Given the observed and expected likelihood profiles as a function of a
 * scaling factor, this finds and returns the values of the scale factor for
 * which the corresponding \f$\mathrm{CL}_s=0.05\f$.
 *
 * The predicted rate should correspond to a scale factor of 1. The likelihood
 * functions qExp and qObs should increase monotonically as a function of the
 * scale factor (at the very least in the vicinity of `scalefactor=1`).
 *
 * Works by numerically findinf the root of `CLs(llh(scalefactor))-0.05` that is
 * closest to `scalefactor=1`.
 *
 * @tparam ExpLlhFunc a function double(double) taking a scale factor and
 * returning the corresponding likelihood value
 * @tparam ObsLlhFunc a function double(double) taking a scale factor and
 * returning the corresponding likelihood value
 * @param qExp function returning the expected likelihood
 * @param qObs function returning the observed likelihood
 * @return expRatio and obsRatio
 */
template <class ExpLlhFunc, class ObsLlhFunc>
REQUIRES(std::regular_invocable<ExpLlhFunc, double>
             &&std::regular_invocable<ObsLlhFunc, double>)
LikelihoodRatios limitsFromLlh(const ExpLlhFunc &qExp, const ObsLlhFunc &qObs) {
    static const auto log = logger();

    const auto exp = [&qExp](double s) { return expCLs(qExp(s)) - limCLs; };
    const auto obs = [&qExp, &qObs](double s) {
        return obsCLs(qExp(s), qObs(s)) - limCLs;
    };

    constexpr auto maxScaleFactor = 1e3;
    if (exp(maxScaleFactor) > 0 || obs(maxScaleFactor) > 0) {
        log->trace("Skipping likelihood limit since it is far from sensitive "
                   "(exp/obsRatio < {})",
                   1 / maxScaleFactor);
        return {0., 0.};
    }

    constexpr auto guess = 1.;
    constexpr auto bracketingFactor = 2.;
    constexpr auto rising = false;
    static const auto tol = utilities::RootFinding::eps_tolerance<double>{};
    constexpr auto maxIter = 50ul;
    log->trace("Starting expected limit search from llh {} (CLs={})",
               qExp(guess), expCLs(qExp(guess)));
    auto expMaxIter = maxIter;
    auto expScale = std::pair<double, double>{};
    try {
        // NOLINTNEXTLINE(clang-analyzer-core.UndefinedBinaryOperatorResult)
        expScale = utilities::RootFinding::bracket_and_solve_root(
            exp, guess, bracketingFactor, rising, tol, expMaxIter);
    } catch (const std::exception &e) { // LCOV_EXCL_START
        logger()->error(e.what());
        logger()->info("Returning obsratio, expratio=0");
        return {0, 0};
    } // LCOV_EXCL_STOP

    log->trace("Obtained expected limit from llh: 1/expRatio in [{}, {}] after "
               "{} iterations",
               expScale.first, expScale.second, expMaxIter);

    log->trace("Starting observed limit search from llh {} (CLs={})",
               qObs(guess), obsCLs(qExp(guess), qObs(guess)));
    auto obsMaxIter = maxIter;

    // NOLINTNEXTLINE(clang-analyzer-core.UndefinedBinaryOperatorResult)
    auto obsScale = std::pair<double, double>{};
    try {
        obsScale = utilities::RootFinding::bracket_and_solve_root(
            obs, guess, bracketingFactor, rising, tol, obsMaxIter);
    } catch (const std::exception &e) { // LCOV_EXCL_START
        logger()->error(e.what());
        logger()->info("Returning obsratio, expratio=0");
        return {0, 0};
    } // LCOV_EXCL_STOP

    log->trace("Obtained observed limit from llh: 1/obsRatio in [{}, {}] after "
               "{} iterations",
               obsScale.first, obsScale.second, obsMaxIter);

    return {1 / expScale.first, 1 / obsScale.first};
}

} // namespace Higgs::bounds::Likelihood
