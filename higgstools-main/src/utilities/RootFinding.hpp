/**
 * @file RootFinding.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief This file is adapted from <boost/math/tools/toms748_solve.hpp>. All of
 * the logic is identical, only the policies have been removed (with the default
 * behaviour hard-coded) and boost::precision calls replaced by
 * std::numeric_limits.
 *
 * @copyright modified version: Copyright 2019 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 * @copyright original boost version: Copyright John Maddock 2006. Use,
 * modification and distribution are subject to the Boost Software License,
 * Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 */

#include "utilities/Format.hpp"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Higgs::utilities::RootFinding {

template <class T> class eps_tolerance {
  public:
    eps_tolerance() { eps = 4 * std::numeric_limits<T>::epsilon(); }
    bool operator()(const T &a, const T &b) {
        return std::fabs(a - b) <= (eps * std::min(std::fabs(a), std::fabs(b)));
    }

  private:
    T eps;
};

namespace detail {

template <class T> inline int sign(const T &z) {
    return (z == 0) ? 0 : std::signbit(z) ? -1 : 1;
}

template <class F, class T>
void bracket(F f, T &a, T &b, T c, T &fa, T &fb, T &d, T &fd) {
    //
    // Given a point c inside the existing enclosing interval
    // [a, b] sets a = c if f(c) == 0, otherwise finds the new
    // enclosing interval: either [a, c] or [c, b] and sets
    // d and fd to the point that has just been removed from
    // the interval.  In other words d is the third best guess
    // to the root.
    //
    T tol = std::numeric_limits<T>::epsilon() * 2;
    //
    // If the interval [a,b] is very small, or if c is too close
    // to one end of the interval then we need to adjust the
    // location of c accordingly:
    //
    if ((b - a) < 2 * tol * a) {
        c = a + (b - a) / 2;
    } else if (c <= a + fabs(a) * tol) {
        c = a + fabs(a) * tol;
    } else if (c >= b - fabs(b) * tol) {
        c = b - fabs(b) * tol;
    }
    //
    // OK, lets invoke f(c):
    //
    T fc = f(c);
    //
    // if we have a zero then we have an exact solution to the root:
    //
    if (fc == 0) {
        a = c;
        fa = 0;
        d = 0;
        fd = 0;
        return;
    }
    //
    // Non-zero fc, update the interval:
    //
    if (sign(fa) * sign(fc) < 0) {
        d = b;
        fd = fb;
        b = c;
        fb = fc;
    } else {
        d = a;
        fd = fa;
        a = c;
        fa = fc;
    }
}

template <class T> inline T safe_div(T num, T denom, T r) {
    //
    // return num / denom without overflow,
    // return r if overflow would occur.
    //

    if (fabs(denom) < 1) {
        if (fabs(denom * std::numeric_limits<T>::max()) <= fabs(num))
            return r;
    }
    return num / denom;
}

template <class T>
inline T secant_interpolate(const T &a, const T &b, const T &fa, const T &fb) {
    //
    // Performs standard secant interpolation of [a,b] given
    // function evaluations f(a) and f(b).  Performs a bisection
    // if secant interpolation would leave us very close to either
    // a or b.  Rationale: we only call this function when at least
    // one other form of interpolation has already failed, so we know
    // that the function is unlikely to be smooth with a root very
    // close to a or b.
    //

    T tol = std::numeric_limits<T>::epsilon() * 5;
    T c = a - (fa / (fb - fa)) * (b - a);
    if ((c <= a + fabs(a) * tol) || (c >= b - fabs(b) * tol))
        return (a + b) / 2;
    return c;
}

template <class T>
T quadratic_interpolate(const T &a, const T &b, T const &d, const T &fa,
                        const T &fb, T const &fd, unsigned count) {
    //
    // Performs quadratic interpolation to determine the next point,
    // takes count Newton steps to find the location of the
    // quadratic polynomial.
    //
    // Point d must lie outside of the interval [a,b], it is the third
    // best approximation to the root, after a and b.
    //
    // Note: this does not guarantee to find a root
    // inside [a, b], so we fall back to a secant step should
    // the result be out of range.
    //
    // Start by obtaining the coefficients of the quadratic polynomial:
    //
    T B = safe_div(T(fb - fa), T(b - a), std::numeric_limits<T>::max());
    T A = safe_div(T(fd - fb), T(d - b), std::numeric_limits<T>::max());
    A = safe_div(T(A - B), T(d - a), T(0));

    if (A == 0) {
        // failure to determine coefficients, try a secant step:
        return secant_interpolate(a, b, fa, fb);
    }
    //
    // Determine the starting point of the Newton steps:
    //
    T c;
    if (sign(A) * sign(fa) > 0) {
        c = a;
    } else {
        c = b;
    }
    //
    // Take the Newton steps:
    //
    for (unsigned i = 1; i <= count; ++i) {
        // c -= safe_div(B * c, (B + A * (2 * c - a - b)), 1 + c - a);
        c -= safe_div(T(fa + (B + A * (c - b)) * (c - a)),
                      T(B + A * (2 * c - a - b)), T(1 + c - a));
    }
    if ((c <= a) || (c >= b)) {
        // Oops, failure, try a secant step:
        c = secant_interpolate(a, b, fa, fb);
    }
    return c;
}

template <class T>
T cubic_interpolate(const T &a, const T &b, const T &d, const T &e, const T &fa,
                    const T &fb, const T &fd, const T &fe) {
    //
    // Uses inverse cubic interpolation of f(x) at points
    // [a,b,d,e] to obtain an approximate root of f(x).
    // Points d and e lie outside the interval [a,b]
    // and are the third and forth best approximations
    // to the root that we have found so far.
    //
    // Note: this does not guarantee to find a root
    // inside [a, b], so we fall back to quadratic
    // interpolation in case of an erroneous result.
    //
    T q11 = (d - e) * fd / (fe - fd);
    T q21 = (b - d) * fb / (fd - fb);
    T q31 = (a - b) * fa / (fb - fa);
    T d21 = (b - d) * fd / (fd - fb);
    T d31 = (a - b) * fb / (fb - fa);

    T q22 = (d21 - q11) * fb / (fe - fb);
    T q32 = (d31 - q21) * fa / (fd - fa);
    T d32 = (d31 - q21) * fd / (fd - fa);
    T q33 = (d32 - q22) * fa / (fe - fa);
    T c = q31 + q32 + q33 + a;

    if ((c <= a) || (c >= b)) {
        // Out of bounds step, fall back to quadratic interpolation:
        c = quadratic_interpolate(a, b, d, fa, fb, fd, 3);
    }

    return c;
}

} // namespace detail

template <class F, class T, class Tol>
std::pair<T, T> toms748_solve(F f, const T &ax, const T &bx, const T &fax,
                              const T &fbx, Tol tol, std::size_t &max_iter) {
    //
    // Main entry point and logic for Toms Algorithm 748
    // root finder.
    //

    //
    // Sanity check - are we allowed to iterate at all?
    //
    if (max_iter == 0)
        return std::make_pair(ax, bx);

    std::size_t count = max_iter;
    T a, b, fa, fb, c, u, fu, a0, b0, d, fd, e, fe;
    static const T mu = 0.5f;

    // initialise a, b and fa, fb:
    a = ax;
    b = bx;

    fa = fax;
    fb = fbx;

    if (tol(a, b) || (fa == 0) || (fb == 0)) {
        max_iter = 0;
        if (fa == 0)
            b = a;
        else if (fb == 0)
            a = b;
        return std::make_pair(a, b);
    }

    if (detail::sign(fa) * detail::sign(fb) > 0)
        throw std::domain_error("Parameters a and b do not bracket the root");
    // dummy value for fd, e and fe:
    fe = e = fd = 1e5F;

    if (fa != 0) {
        //
        // On the first step we take a secant step:
        //
        c = detail::secant_interpolate(a, b, fa, fb);
        detail::bracket(f, a, b, c, fa, fb, d, fd);
        --count;

        if (count && (fa != 0) && !tol(a, b)) {
            //
            // On the second step we take a quadratic interpolation:
            //
            c = detail::quadratic_interpolate(a, b, d, fa, fb, fd, 2);
            e = d;
            fe = fd;
            detail::bracket(f, a, b, c, fa, fb, d, fd);
            --count;
        }
    }

    while (count && (fa != 0) && !tol(a, b)) {
        // save our brackets:
        a0 = a;
        b0 = b;
        //
        // Starting with the third step taken
        // we can use either quadratic or cubic interpolation.
        // Cubic interpolation requires that all four function values
        // fa, fb, fd, and fe are distinct, should that not be the case
        // then variable prof will get set to true, and we'll end up
        // taking a quadratic step instead.
        //
        T min_diff = std::numeric_limits<T>::min() * 32;
        bool prof = (fabs(fa - fb) < min_diff) || (fabs(fa - fd) < min_diff) ||
                    (fabs(fa - fe) < min_diff) || (fabs(fb - fd) < min_diff) ||
                    (fabs(fb - fe) < min_diff) || (fabs(fd - fe) < min_diff);
        if (prof) {
            c = detail::quadratic_interpolate(a, b, d, fa, fb, fd, 2);
        } else {
            c = detail::cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
        }
        //
        // re-bracket, and check for termination:
        //
        e = d;
        fe = fd;
        detail::bracket(f, a, b, c, fa, fb, d, fd);
        if ((0 == --count) || (fa == 0) || tol(a, b))
            break;
        //
        // Now another interpolated step:
        //
        prof = (fabs(fa - fb) < min_diff) || (fabs(fa - fd) < min_diff) ||
               (fabs(fa - fe) < min_diff) || (fabs(fb - fd) < min_diff) ||
               (fabs(fb - fe) < min_diff) || (fabs(fd - fe) < min_diff);
        if (prof) {
            c = detail::quadratic_interpolate(a, b, d, fa, fb, fd, 3);
        } else {
            c = detail::cubic_interpolate(a, b, d, e, fa, fb, fd, fe);
        }
        //
        // Bracket again, and check termination condition, update e:
        //
        detail::bracket(f, a, b, c, fa, fb, d, fd);
        if ((0 == --count) || (fa == 0) || tol(a, b))
            break;
        //
        // Now we take a double-length secant step:
        //
        if (fabs(fa) < fabs(fb)) {
            u = a;
            fu = fa;
        } else {
            u = b;
            fu = fb;
        }
        c = u - 2 * (fu / (fb - fa)) * (b - a);
        if (fabs(c - u) > (b - a) / 2) {
            c = a + (b - a) / 2;
        }
        //
        // Bracket again, and check termination condition:
        //
        e = d;
        fe = fd;
        detail::bracket(f, a, b, c, fa, fb, d, fd);
        if ((0 == --count) || (fa == 0) || tol(a, b))
            break;
        //
        // And finally... check to see if an additional bisection step is
        // to be taken, we do this if we're not converging fast enough:
        //
        if ((b - a) < mu * (b0 - a0))
            continue;
        //
        // bracket again on a bisection:
        //
        e = d;
        fe = fd;
        detail::bracket(f, a, b, T(a + (b - a) / 2), fa, fb, d, fd);
        --count;
    } // while loop

    max_iter -= count;
    if (fa == 0) {
        b = a;
    } else if (fb == 0) {
        a = b;
    }
    return std::make_pair(a, b);
}

template <class F, class T, class Tol>
std::pair<T, T> bracket_and_solve_root(F f, const T &guess, T factor,
                                       bool rising, Tol tol,
                                       std::size_t &max_iter) {
    //
    // Set up initial brackets:
    //
    T a = guess;
    T b = a;
    T fa = f(a);
    T fb = fa;
    //
    // Set up invocation count:
    //
    std::size_t count = max_iter - 1;

    std::size_t step = 32;

    if ((fa < 0) == (guess < 0 ? !rising : rising)) {
        //
        // Zero is to the right of b, so walk upwards
        // until we find it:
        //
        while (detail::sign(fb) == detail::sign(fa)) {
            if (count == 0)
                throw std::domain_error(fmt::format(
                    "Unable to bracket root, last nearest value was {}", b));
            //
            // Heuristic: normally it's best not to increase the step sizes as
            // we'll just end up with a really wide range to search for the
            // root.  However, if the initial guess was *really* bad then we
            // need to speed up the search otherwise we'll take forever if we're
            // orders of magnitude out.  This happens most often if the guess is
            // a small value (say 1) and the result we're looking for is close
            // to std::numeric_limits<T>::min().
            //
            if ((max_iter - count) % step == 0) {
                factor *= 2;
                if (step > 1)
                    step /= 2;
            }
            //
            // Now go ahead and move our guess by "factor":
            //
            a = b;
            fa = fb;
            b *= factor;
            fb = f(b);
            --count;
        }
    } else {
        //
        // Zero is to the left of a, so walk downwards
        // until we find it:
        //
        while (detail::sign(fb) == detail::sign(fa)) {
            if (fabs(a) < std::numeric_limits<T>::min()) {
                // Escape route just in case the answer is zero!
                max_iter -= count;
                max_iter += 1;
                return a > 0 ? std::make_pair(T(0), T(a))
                             : std::make_pair(T(a), T(0));
            }
            if (count == 0)
                throw std::runtime_error(fmt::format(
                    "Unable to bracket root, last nearest value was {}", a));
            //
            // Heuristic: normally it's best not to increase the step sizes as
            // we'll just end up with a really wide range to search for the
            // root.  However, if the initial guess was *really* bad then we
            // need to speed up the search otherwise we'll take forever if we're
            // orders of magnitude out.  This happens most often if the guess is
            // a small value (say 1) and the result we're looking for is close
            // to std::numeric_limits<T>::min().
            //
            if ((max_iter - count) % step == 0) {
                factor *= 2;
                if (step > 1)
                    step /= 2;
            }
            //
            // Now go ahead and move are guess by "factor":
            //
            b = a;
            fb = fa;
            a /= factor;
            fa = f(a);
            --count;
        }
    }
    max_iter -= count;
    max_iter += 1;
    std::pair<T, T> r =
        toms748_solve(f, (a < 0 ? b : a), (a < 0 ? a : b), (a < 0 ? fb : fa),
                      (a < 0 ? fa : fb), tol, count);
    max_iter += count;
    return r;
}

} // namespace Higgs::utilities::RootFinding
