#pragma once

#include <magic_enum.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace Higgs {
//! Get the enum type name without any namespaces
template <class Enum> constexpr auto enumTypeName() {
    constexpr auto fullName = magic_enum::enum_type_name<Enum>();
    constexpr auto pos = fullName.find_last_of(':');
    if (pos != std::string_view::npos) {
        return fullName.substr(pos + 1);
    }
    return fullName;
}
// register an Enum type with implicit string conversions in pybind
template <class Enum> void registerEnum(const py::module &m) {
    static constexpr auto typeName = enumTypeName<Enum>();
    auto e = py::enum_<Enum>(m, typeName.data());
    for (const auto &[val, name] : magic_enum::enum_entries<Enum>()) {
        e.value(name.data(), val);
    }
    e.def(py::init([](const std::string &name) {
        auto cast = magic_enum::enum_cast<Enum>(name);
        if (!cast) {
            throw std::invalid_argument(
                name + " is not a valid value for " + std::string(typeName) +
                ", see " + std::string(typeName) + ".__members__().");
        }
        return cast.value();
    }));
    py::implicitly_convertible<std::string, Enum>();
}
} // namespace Higgs
namespace Higgs::predictions {
void bindHiggsPredictions(py::module &higgs);
}

namespace Higgs::bounds {
void bindHiggsBounds(py::module &higgs);
}

namespace Higgs::signals {
void bindHiggsSignals(py::module &higgs);
}
