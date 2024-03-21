#pragma once

#include "utilities/Concepts.hpp"
#include "utilities/JsonFwd.hpp" // IWYU pragma: export
#include <magic_enum.hpp>        // IWYU pragma: export
#include <nlohmann/json.hpp>     // IWYU pragma: export
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

//! define enum<->json conversion functions using the enum field names
#define HIGGSUTILITIES_ENUM_JSON_CONV(ENUM)                                    \
    void to_json(nlohmann::json &j, const ENUM &e) {                           \
        j = magic_enum::enum_name(e);                                          \
    }                                                                          \
    void from_json(const nlohmann::json &j, ENUM &e) {                         \
        auto optE = magic_enum::enum_cast<ENUM>(j.get<std::string_view>());    \
        if (optE) {                                                            \
            e = optE.value();                                                  \
        } else {                                                               \
            throw Higgs::utilities::BadEnumRead(                               \
                "`" + j.get<std::string>() +                                   \
                "` is not a valid value for " #ENUM);                          \
        }                                                                      \
    }

namespace Higgs::utilities {

//! Error reading an enum value from json
class BadEnumRead : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

//! Error reading a field from json
class BadFieldRead : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Read a field from the json as the specified type
 * Utility function to provide better error messages (that include the field
 * name).
 * @throws BadFieldRead if the read failed
 * @tparam T the type to read as
 * @param j the json data
 * @param key the key, if it starts with a `/` it is interpreted as a json
 * pointer
 * @return T the read value
 */
template <typename T>
T readAs(const nlohmann::json &j, const std::string &key) {
    try {
        if (!key.empty() && key.front() == '/')
            return j.at(nlohmann::json::json_pointer{key}).get<T>();
        else
            return j.at(key).get<T>();
    } catch (const nlohmann::json::type_error &e) {
        throw BadFieldRead("error in field `" + key + "`: " + e.what());
    }
}

/**
 * @brief Read the field from json as the specified type or return an
 * object constructed from args.
 *
 * @tparam T the type to read
 * @tparam Args... argument types to the constructor of T
 * @param j the json data
 * @param key the key, if it stars with a `/` it is interpreted as a json
 * pointer
 * @param args arguments to the constructor in case the key is not present
 * @return T the read value or `T(args...)`
 */
template <typename T, class... Args> REQUIRES(std::semiregular<T>)
T readIfPresent(const nlohmann::json &j, const std::string &key,
                Args &&...args) {
    if (j.contains(key) || (!key.empty() && key.front() == '/' &&
                            j.contains(nlohmann::json::json_pointer{key}))) {
        return readAs<T>(j, key);
    }
    return T(std::forward<Args>(args)...);
}
/**
 * @brief Read a field from json with an optionally provided default value.
 *
 * If the field does not exist in the json the default value (if present) is
 * used. Otherwise calls readAs.
 *
 * @tparam T the type to read
 * @param j the json data
 * @param key the key, no json pointer support
 * @param defaultVal a default value if the field is not present.
 * @return T the read value or the default
 */
template <typename T>
T readWithOptDefault(const nlohmann::json &j, const std::string &key,
                     const std::optional<T> &defaultVal) {
    if (defaultVal && !j.contains(key))
        return defaultVal.value();
    return readAs<T>(j, key);
}

} // namespace Higgs::utilities
