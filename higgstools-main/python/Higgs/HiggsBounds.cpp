#include "Higgs.hpp"
#include "Higgs/Bounds.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/bounds/Limit.hpp"
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <magic_enum.hpp>
#include <map>
#include <memory>
#include <pybind11/stl.h>
#include <sstream>
#include <string>
#include <vector>

using namespace py::literals;

namespace {
constexpr Higgs::bounds::LimitOptions defaultOptions{};
}
namespace Higgs::bounds {

void bindHiggsBounds(py::module &higgs) {
    auto hb = higgs.def_submodule(
        "bounds", "Exclusion bounds from particle searches at colliders.");

    py::class_<AppliedLimit>(hb, "AppliedLimit")
        .def("limit", &AppliedLimit::limit)
        .def("obsRatio", &AppliedLimit::obsRatio)
        .def("expRatio", &AppliedLimit::expRatio)
        .def("obsLikelihood", &AppliedLimit::obsLikelihood)
        .def("expLikelihood", &AppliedLimit::expLikelihood)
        .def("contributingParticles", &AppliedLimit::contributingParticles)
        .def("__repr__",
             [](const AppliedLimit &self) {
                 return fmt::format(
                     "<Higgs.bounds.AppliedLimit limitId: {}, obsRatio: "
                     "{:5.3f}, expRatio: {:5.3f}, for: {}>",
                     self.limit()->id(), self.obsRatio(), self.expRatio(),
                     self.contributingParticles());
             })
        .def("__str__", [](const AppliedLimit &self) {
            return fmt::format(
                "obsRatio {:5.3f}, expRatio: {:5.3f} for {} with {}",
                self.obsRatio(), self.expRatio(), self.contributingParticles(),
                self.limit()->to_string());
        });

    py::class_<LimitOptions>(hb, "LimitOptions")
        .def(py::init<double, predictions::MassUncEagerness,
                      predictions::MassUncEagerness,
                      predictions::MassUncEagerness, double>(),
             "applicableResolutionFac"_a =
                 defaultOptions.applicableResolutionFac,
             "clusterMassUnc"_a = defaultOptions.clusterMassUnc,
             "applicableMassUnc"_a = defaultOptions.applicableMassUnc,
             "setLimitMassUnc"_a = defaultOptions.setLimitMassUnc,
             "minExpRatio"_a = defaultOptions.minExpRatio)
        .def("__repr__",
             [](const LimitOptions &self) {
                 return fmt::format(
                     "Higgs.bounds.LimitOptions(applicableResolutionFac = {}, "
                     "clusterMassUnc = \"{}\", applicableMassUnc = \"{}\", "
                     "setLimitMassUnc = \"{}\", minExpRatio = {})",
                     self.applicableResolutionFac,
                     magic_enum::enum_name(self.clusterMassUnc),
                     magic_enum::enum_name(self.applicableMassUnc),
                     magic_enum::enum_name(self.setLimitMassUnc),
                     self.minExpRatio);
             })
        .def_readwrite("applicableResolutionFac",
                       &LimitOptions::applicableResolutionFac)
        .def_readwrite("clusterMassUnc", &LimitOptions::clusterMassUnc)
        .def_readwrite("applicableMassUnc", &LimitOptions::applicableMassUnc)
        .def_readwrite("setLimitMassUnc", &LimitOptions::setLimitMassUnc);

    py::class_<Limit, std::shared_ptr<Limit>>(hb, "Limit")
        .def(py::init(&Limit::read), "filePath"_a, "options"_a = LimitOptions{})
        .def("__repr__",
             [](const Limit &l) {
                 return "<Higgs.bounds.Limit " + std::to_string(l.id()) + ": " +
                        l.to_string() + ">";
             })
        .def("__str__", [](const Limit &l) { return l.to_string(); })
        .def("apply", &Limit::apply)
        .def("id", &Limit::id)
        .def("processDesc", &Limit::processDesc)
        .def("extentDesc", &Limit::extentDesc)
        .def("reference", &Limit::reference)
        .def("citeKey", &Limit::citeKey)
        .def("collider", &Limit::collider)
        .def("experiment", &Limit::experiment)
        .def("luminosity", &Limit::luminosity)
        .def("loadedFrom", &Limit::loadedFrom)
        .def("options", &Limit::options);

    py::class_<HBResult>(hb, "HBResult")
        .def_readonly("selectedLimits", &HBResult::selectedLimits)
        .def_readonly("appliedLimits", &HBResult::appliedLimits)
        .def_readonly("allowed", &HBResult::allowed)
        .def("__bool__", [](const HBResult &x) { return x.allowed; })
        .def("__str__",
             [](const HBResult &self) {
                 auto ss = std::ostringstream{};
                 ss << self;
                 return ss.str();
             })
        .def("__repr__", [](const HBResult &self) {
            return fmt::format("<Higgs.bounds.HBResult {}>", self.allowed);
        });

    py::class_<Bounds>(hb, "Bounds")
        .def(py::init<const std::string &, const LimitOptions &>(),
             "dataPath"_a, "limitOptions"_a = LimitOptions{})
        .def("limits", &Bounds::limits)
        .def("__call__", &Bounds::operator());
}

} // namespace Higgs::bounds
