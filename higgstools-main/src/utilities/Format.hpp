#pragma once

#include <fmt/core.h>     // IWYU pragma: export
#include <fmt/format.h>   // IWYU pragma: export
#include <fmt/ranges.h>   // IWYU pragma: export
#include <magic_enum.hpp> // IWYU pragma: export

#define HIGGSUTILITIES_ENUM_FORMATTER(ENUM)                                    \
    namespace fmt {                                                            \
    template <> struct formatter<ENUM> : formatter<string_view> {              \
        template <typename FormatContext>                                      \
        auto format(ENUM e, FormatContext &ctx) {                              \
            return formatter<string_view>::format(magic_enum::enum_name(e),    \
                                                  ctx);                        \
        }                                                                      \
    };                                                                         \
    }
