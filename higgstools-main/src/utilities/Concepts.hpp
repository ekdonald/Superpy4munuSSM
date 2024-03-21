/**
 * @file Concepts.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Concepts if C++-20 is available
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#if defined(__cpp_concepts) && defined(__has_include) && __has_include(<concepts>)
#define HIGGSBOUNDS_HAS_CONCEPTS
#endif

#ifdef HIGGSBOUNDS_HAS_CONCEPTS
#include <concepts>
#include <magic_enum.hpp>
#include <ranges>

namespace Higgs::predictions {
class Particle;
}

namespace Higgs::utilities::concepts {
//! A scoped enum/enum class
template <class T> concept ScopedEnum = magic_enum::is_scoped_enum_v<T>;

//! A particle.
template <class T>
concept Particle = std::derived_from<T, Higgs::predictions::Particle>;

//! A range of particles
template <class T>
concept ParticleRange =
    std::ranges::forward_range<T> &&Particle<std::ranges::range_value_t<T>>;

template <class T>
concept EnumWithNone = ScopedEnum<T> &&magic_enum::enum_contains(T::none);

//! A type that can be indexed using operator[]
template <class T> concept Container = requires(T a, typename T::size_type i) {
    typename T::value_type;
    a[i];
};

//! A type that supports addition
template <typename T> concept Addable = requires(T x) {
    { x += x }
    ->std::convertible_to<T>;
    { x + x }
    ->std::convertible_to<T>;
};

//! A type that supports multiplication with another type
template <typename T, typename S>
concept MultiplicableWith = requires(T x, S s) {
    { x *s }
    ->std::convertible_to<T>;
};

} // namespace Higgs::utilities::concepts

//! expands to a requires clause if concepts are supported
#define REQUIRES(...) requires __VA_ARGS__
#else
#define REQUIRES(...)
#endif
