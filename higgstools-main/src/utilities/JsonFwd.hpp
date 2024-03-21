#pragma once

#include <nlohmann/json_fwd.hpp> // IWYU pragma: export

//! declare TYPE <-> json conversion functions
#define HIGGSUTILITIES_JSON_CONV_FWD(TYPE)                                     \
    void to_json(nlohmann::json &j, const TYPE &e);                            \
    void from_json(const nlohmann::json &j, TYPE &e);
