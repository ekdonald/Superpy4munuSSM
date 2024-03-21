/**
 * @file Algorithm.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Algorithms used throughout the library
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "utilities/Concepts.hpp"
#include <array>
#include <functional>
#include <limits>
#include <range/v3/functional/identity.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/traits.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/transform.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Higgs::utilities {

/**
 * @brief Get an element from a std compliant map (an AssociativeContainer) or a
 * default value it the key does not exist.
 *
 * @tparam Map the AssociativeContainer type
 * @param map the map to read from
 * @param key the key to look for
 * @param defaultValue the default value if the key is not found
 * @return Map::mapped_type the value corresponding to key or the default value
 */
template <class Map>
typename Map::mapped_type
getWithDefault(const Map &map, const typename Map::key_type &key,
               const typename Map::mapped_type &defaultValue) noexcept {
    auto it = map.find(key);
    if (it == map.end())
        return defaultValue;
    return it->second;
}

/**
 * @brief Accumulates the results of calling value on every element of the
 * cartesian product of the ranges.
 *
 * @tparam ValueFunc function that given one element from each range returns the
 * values that should be accumulated
 * @tparam Ranges any number of ranges
 */
template <class ValueFunc, class... Ranges>
REQUIRES(std::regular_invocable<ValueFunc, ranges::range_value_t<Ranges>...>)
auto accumulatedCartesianProduct(ValueFunc value, Ranges &&...ranges) {
    using Result =
        std::invoke_result_t<ValueFunc, ranges::range_value_t<Ranges>...>;
    const auto getValue = [&value](const auto &elements) {
        return std::apply(value, elements);
    };
    return ranges::accumulate(
        ranges::views::cartesian_product(std::forward<Ranges>(ranges)...),
        Result{}, std::plus{}, getValue);
}

/**
 * Calculates the mean over a range. Averages the values returned by
 * `valueFunc(*x)` for `x` in the range.
 * @tparam Range a range
 * @tparam WeightFunc a function that is invocable on the elements of the
 * range
 * @tparam ValueFunc a function that is invocable on the elements of the
 * range
 * @param range the range
 * @param weightFunc function that returns weights as **non-negative**
 * doubles
 * @param valueFunc  function that returns the values
 * @return weighted mean
 */
template <class Range, class ValueFunc = ranges::identity>
REQUIRES(std::ranges::input_range<Range>
             &&std::regular_invocable<ValueFunc, ranges::range_value_t<Range>>)
auto mean(Range &&range, const ValueFunc &value = ranges::identity{}) {
    using Result =
        std::invoke_result_t<ValueFunc, ranges::range_value_t<Range>>;
    const auto updateMean = [&value](const Result &val, const auto &xx) {
        const auto &[n, x] = xx;
        return val + (value(x) - val) / static_cast<double>(n + 1);
    };
    return ranges::accumulate(ranges::views::enumerate(range), Result{},
                              updateMean);
}

/**
 * Calculates the weighted mean over a range. Averages the values returned by
 * `valueFunc(*x)` based on the weights `weightFunc(*x)` for `x` in the range.
 * @tparam Range a range
 * @tparam ValueFunc a function that is invocable on the elements of the range
 * @tparam WeightFunc a function that is invocable on the elements of the range
 * @param range the range
 * @param weight returns the weight of an element as a **non-negative** double
 * @param value returns the value of a range element
 * @return mean
 */
template <class Range, class WeightFunc, class ValueFunc>
REQUIRES(
    std::ranges::input_range<Range>
        &&std::regular_invocable<WeightFunc, ranges::range_value_t<Range>>
            &&std::regular_invocable<ValueFunc, ranges::range_value_t<Range>>)
auto weightedMean(Range &&range, const WeightFunc &weight,
                  const ValueFunc &value) {
    using Result =
        std::invoke_result_t<ValueFunc, ranges::range_value_t<Range>>;

    if (ranges::distance(range) == 1) {
        return value(*range.begin());
    }
    auto addValueAndWeight = [&weight, &value](const auto &prev,
                                               const auto &x) {
        const auto w = weight(x);
        const auto &[prevVal, prevW] = prev;
        return std::pair{prevVal + w * value(x), prevW + w};
    };
    const auto sum = ranges::accumulate(range, std::pair<Result, double>{},
                                        addValueAndWeight);
    if (sum.second < 100 * std::numeric_limits<double>::denorm_min()) {
        return Result{};
    } else {
        return sum.first / sum.second;
    }
}

/**
 * @brief Calculates the weighted mean over two ranges for a weight that
 * depends on both. Returns a two-component array where the elements are the
 * weighted means of the values for each range. The weight for a fixed
 * element e1 in the first range is calculated by summing `weight(e1,e2)`
 * for all `e2` in the other range and vice-versa.
 *
 * @tparam R1 first range
 * @tparam R2 second range
 * @tparam WeightFunc function that is invocable on one element from each
 * range (in order)
 * @tparam ValueFunc a function that is invocable on each element of either
 * range
 * @param r1 the first range
 * @param r2 the second range
 * @param weight returns the weight of a combination of range elements as a
 * **non-negative** double
 * @param value returns the value of a range element
 * @returns a size=2 std::array of the weighted values for each range
 */
template <class R1, class R2, class WeightFunc, class ValueFunc>
REQUIRES(std::regular_invocable<WeightFunc, ranges::range_value_t<R1>,
                                ranges::range_value_t<R2>>
             &&std::regular_invocable<ValueFunc, ranges::range_value_t<R1>>
                 &&std::regular_invocable<ValueFunc, ranges::range_value_t<R2>>)
auto weightedMean(R1 &&r1, R2 &&r2, const WeightFunc &weight,
                  const ValueFunc &value) {
    const auto weight1 = [&r2, &weight](const auto &e) {
        return accumulatedCartesianProduct(
            weight, ranges::views::single(std::cref(e)), r2);
    };
    const auto weight2 = [&r1, &weight](const auto &e) {
        return accumulatedCartesianProduct(weight, r1,
                                           ranges::views::single(std::cref(e)));
    };
    return std::array{weightedMean(r1, weight1, value),
                      weightedMean(r2, weight2, value)};
}

//! Calculates the weighted mean over three ranges for a weight that depends
//! on all three. See the two-range overload for details.
template <class R1, class R2, class R3, class WeightFunc, class ValueFunc>
REQUIRES(
    std::regular_invocable<WeightFunc, ranges::range_value_t<R1>,
                           ranges::range_value_t<R2>, ranges::range_value_t<R3>>
        &&std::regular_invocable<ValueFunc, ranges::range_value_t<R1>>
            &&std::regular_invocable<ValueFunc, ranges::range_value_t<R2>>
                &&std::regular_invocable<ValueFunc, ranges::range_value_t<R3>>)
auto weightedMean(R1 &&r1, R2 &&r2, R3 &&r3, const WeightFunc &weight,
                  const ValueFunc &value) {
    const auto weight1 = [&r2, &r3, &weight](const auto &e) {
        return accumulatedCartesianProduct(
            weight, ranges::views::single(std::cref(e)), r2, r3);
    };
    const auto weight2 = [&r1, &r3, &weight](const auto &e) {
        return accumulatedCartesianProduct(
            weight, r1, ranges::views::single(std::cref(e)), r3);
    };
    const auto weight3 = [&r1, &r2, &weight](const auto &e) {
        return accumulatedCartesianProduct(weight, r1, r2,
                                           ranges::views::single(std::cref(e)));
    };
    return std::array{weightedMean(r1, weight1, value),
                      weightedMean(r2, weight2, value),
                      weightedMean(r3, weight3, value)};
}

template <class Set, class Predicate>
REQUIRES(std::regular_invocable<Predicate, Set::value_type>)
void eraseFromSetIf(Set &set, const Predicate &pred) {
    for (auto i = set.begin(), last = set.end(); i != last;) {
        if (pred(*i)) {
            i = set.erase(i);
        } else {
            ++i;
        }
    }
}

} // namespace Higgs::utilities
