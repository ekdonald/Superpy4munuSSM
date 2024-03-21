/**
 * @file ArithmeticArray.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Array with arithmetic operations
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include <array>
#include <cstddef>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <type_traits>

namespace Higgs::utilities {

/**
 * A statically sized array class with elementwise arithmetic operations.
 *
 * All arithmetic operators work elementwise for both two arrays and an array
 * and a value convertible to T. Inherits from std::array so the full
 * standard library interface is available as well.
 *
 * @tparam T value type
 * @tparam N number of elements
 */
template <class T, std::size_t N>
struct ArithmeticArray : public std::array<T, N> {
    static_assert(std::is_arithmetic<T>::value);

    //! unary -
    ArithmeticArray operator-() const {
        auto res = *this;
        ranges::for_each(res, [](T &t) { t = -t; });
        return res;
    }

    //! @{
    //! compound assignment operators for two arrays
    ArithmeticArray &operator+=(const std::array<T, N> &rhs) & {
        ranges::transform(*this, rhs, this->begin(),
                          [](const T &t, const T &r) { return t + r; });
        return *this;
    }
    ArithmeticArray &operator-=(const std::array<T, N> &rhs) & {
        ranges::transform(*this, rhs, this->begin(),
                          [](const T &t, const T &r) { return t - r; });
        return *this;
    }
    ArithmeticArray &operator*=(const std::array<T, N> &rhs) & {
        ranges::transform(*this, rhs, this->begin(),
                          [](const T &t, const T &r) { return t * r; });
        return *this;
    }
    ArithmeticArray &operator/=(const std::array<T, N> &rhs) & {
        ranges::transform(*this, rhs, this->begin(),
                          [](const T &t, const T &r) { return t / r; });
        return *this;
    }
    //! @}

    //! @{
    //! compound assignment operators for arrays and scalars
    ArithmeticArray &operator+=(const T &s) & {
        ranges::for_each(*this, [&s](T &t) { t += s; });
        return *this;
    }
    ArithmeticArray &operator-=(const T &s) & {
        ranges::for_each(*this, [&s](T &t) { t -= s; });
        return *this;
    }
    ArithmeticArray &operator*=(const T &s) & {
        ranges::for_each(*this, [&s](T &t) { t *= s; });
        return *this;
    }
    ArithmeticArray &operator/=(const T &s) & {
        ranges::for_each(*this, [&s](T &t) { t /= s; });
        return *this;
    }
    //! @}

    //! @{
    //! binary arithmetic operators for two arrays
    friend ArithmeticArray operator+(ArithmeticArray lhs,
                                     std::array<T, N> const &rhs) {
        lhs += rhs;
        return lhs;
    }
    friend ArithmeticArray operator-(ArithmeticArray lhs,
                                     std::array<T, N> const &rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend ArithmeticArray operator*(ArithmeticArray lhs,
                                     std::array<T, N> const &rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend ArithmeticArray operator/(ArithmeticArray lhs,
                                     std::array<T, N> const &rhs) {
        lhs /= rhs;
        return lhs;
    }
    //! @}

    //! @{
    //! binary arithmetic operators for array and scalar
    friend ArithmeticArray operator+(ArithmeticArray lhs, const T &s) {
        lhs += s;
        return lhs;
    }
    friend ArithmeticArray operator-(ArithmeticArray lhs, const T &s) {
        lhs -= s;
        return lhs;
    }
    friend ArithmeticArray operator*(ArithmeticArray lhs, const T &s) {
        lhs *= s;
        return lhs;
    }
    friend ArithmeticArray operator/(ArithmeticArray lhs, const T &s) {
        lhs /= s;
        return lhs;
    }
    friend ArithmeticArray operator+(const T &s, ArithmeticArray rhs) {
        rhs += s;
        return rhs;
    }
    friend ArithmeticArray operator-(const T &s, ArithmeticArray rhs) {
        rhs -= s;
        return -rhs;
    }
    friend ArithmeticArray operator*(const T &s, ArithmeticArray rhs) {
        rhs *= s;
        return rhs;
    }
    friend ArithmeticArray operator/(const T &s, ArithmeticArray rhs) {
        ranges::for_each(rhs, [&s](T &t) { t = s / t; });
        return rhs;
    }
    //! @}

    //! @{
    //! get functions to enable structured bindings
    template <std::size_t I> auto const &get() const & { return (*this)[I]; }
    template <std::size_t I> auto &get() & { return (*this)[I]; }
    template <std::size_t I> auto &&get() && { return (*this)[I]; }
    //! @}
};

//! deduction guide for ArithmeticArray @relates ArithmeticArray
template <class T, class... U>
ArithmeticArray(T, U...) -> ArithmeticArray<T, 1 + sizeof...(U)>;

} // namespace Higgs::utilities

namespace std {
//! specialization to enable structured bindings  @relates ArithmeticArray
template <class T, std::size_t N>
struct tuple_size<Higgs::utilities::ArithmeticArray<T, N>>
    : std::integral_constant<size_t, N> {};
//! specialization to enable structured bindings @relates ArithmeticArray
template <class T, std::size_t N, std::size_t I>
struct tuple_element<I, Higgs::utilities::ArithmeticArray<T, N>> {
    using type = T; //!< all tuple elements are the value_type of the array
};

} // namespace std
