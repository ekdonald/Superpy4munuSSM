#include "Higgs.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/Signals.hpp"
#include "Higgs/signals/Measurement.hpp"
#include <fmt/format.h>
#include <pybind11/stl.h>

using namespace py::literals;

namespace Higgs::signals {

namespace {
constexpr MeasurementOptions defaultOptions{};
}

void bindHiggsSignals(py::module &higgs) {
    auto hs = higgs.def_submodule(
        "signals", "Tests agreement with measurements of scalar resonances.");

    registerEnum<PDF>(hs);
    registerEnum<RescaleToRefMass>(hs);
    registerEnum<Uncertainty>(hs);
    registerEnum<Correlations>(hs);
    registerEnum<NormalizeAt>(hs);

    py::class_<SubMeasurementEvaluation>(hs, "SubMeasurementEvaluation")
        .def_readonly("residual", &SubMeasurementEvaluation::residual)
        .def_readonly("obsVariance", &SubMeasurementEvaluation::obsVariance)
        .def_readonly("refVariance", &SubMeasurementEvaluation::refVariance)
        .def_readonly("extraChisq", &SubMeasurementEvaluation::extraChisq);

    py::class_<SubMeasurement, std::shared_ptr<SubMeasurement>>(
        hs, "SubMeasurement")
        .def("signalStrength", &SubMeasurement::signalStrength)
        .def("chisq", &SubMeasurement::chisq)
        .def("__repr__",
             [](const SubMeasurement &self) {
                 return fmt::format("<Higgs.signals.SubMeasurement of {}>",
                                    self.processDesc(false));
             })
        .def("processDesc", &SubMeasurement::processDesc,
             "keepOrder"_a = false);

    py::class_<MeasurementOptions>(hs, "MeasurementOptions")
        .def(py::init<PDF, double, double, double, RescaleToRefMass,
                      Correlations, bool>(),
             "theoryMassUncPDF"_a = defaultOptions.theoryMassUncPDF,
             "massSensitiveAssignmentRange"_a =
                 defaultOptions.massSensitiveAssignmentRange,
             "unassignedMassMeasurementPenalty"_a =
                 defaultOptions.unassignedMassMeasurementPenalty,
             "unassignedCouplingMeasurementPenalty"_a =
                 defaultOptions.unassignedCouplingMeasurementPenalty,
             "rescaleToRefMass"_a = defaultOptions.rescaleToRefMass,
             "whichCorrelations"_a = defaultOptions.whichCorrelations,
             "ignoreTheoryUncertainties"_a =
                 defaultOptions.ignoreTheoryUncertainties)
        .def("__repr__",
             [](const MeasurementOptions &self) {
                 return fmt::format(
                     "Higgs.signals.MeasurementOptions(massPDF = \"{}\", "
                     "massSensitiveAssignmentRange = {}, "
                     "unassignedMassMeasurementPenalty = {}, "
                     "unassignedCouplingMeasurementPenalty = {}, "
                     "rescaleToRefMass = \"{}\", "
                     "whichCorrelations = \"{}\", "
                     "ignoreTheoryUncertainties = {})",
                     magic_enum::enum_name(self.theoryMassUncPDF),
                     self.massSensitiveAssignmentRange,
                     self.unassignedMassMeasurementPenalty,
                     self.unassignedCouplingMeasurementPenalty,
                     magic_enum::enum_name(self.rescaleToRefMass),
                     magic_enum::enum_name(self.whichCorrelations),
                     self.ignoreTheoryUncertainties ? "True" : "False");
             })
        .def_readwrite("massPDF", &MeasurementOptions::theoryMassUncPDF)
        .def_readwrite("massSensitiveAssignmentRange",
                       &MeasurementOptions::massSensitiveAssignmentRange)
        .def_readwrite("unassignedMassMeasurementPenalty",
                       &MeasurementOptions::unassignedMassMeasurementPenalty)
        .def_readwrite(
            "unassignedCouplingMeasurementPenalty",
            &MeasurementOptions::unassignedCouplingMeasurementPenalty)
        .def_readwrite("rescaleToRefMass",
                       &MeasurementOptions::rescaleToRefMass)
        .def_readwrite("whichCorrelations",
                       &MeasurementOptions::whichCorrelations)
        .def_readwrite("ignoreTheoryUncertainties",
                       &MeasurementOptions::ignoreTheoryUncertainties);

    py::class_<Measurement>(hs, "Measurement")
        .def(py::init<std::string, MeasurementOptions>(), "filePath"_a,
             "options"_a = MeasurementOptions{})
        .def("subMeasurements", &Measurement::subMeasurements,
             py::return_value_policy::reference_internal)
        .def("nSubMeasurements", &Measurement::nSubMeasurements)
        .def("referenceMass", &Measurement::referenceMass)
        .def("referenceModel", &Measurement::referenceModel)
        .def("__call__", &Measurement::operator(), "predictions"_a,
             "modificationFactors"_a = Measurement::ModificationFactors{})
        .def("chisqContributions", &Measurement::chisqContributions,
             "predictions"_a,
             "modificationFactors"_a = Measurement::ModificationFactors{})
        .def("id", &Measurement::id)
        .def("reference", &Measurement::reference)
        .def("citeKey", &Measurement::citeKey)
        .def("collider", &Measurement::collider)
        .def("experiment", &Measurement::experiment)
        .def("luminosity", &Measurement::luminosity)
        .def("loadedFrom", &Measurement::loadedFrom)
        .def("options", &Measurement::options)
        .def("__repr__", [](const Measurement &self) {
            return fmt::format("<Higgs.signals.Measurement {}: {} bin from {}>",
                               self.id(), self.nSubMeasurements(),
                               self.reference());
        });

    py::class_<Signals>(hs, "Signals")
        .def(py::init<std::string, MeasurementOptions>(), "dataPath"_a,
             "measurementOptions"_a = MeasurementOptions{})
        .def("__call__", &Signals::operator(), "predictions"_a,
             "modificationFactors"_a = Signals::ModificationFactors{})
        .def("observableCount", &Signals::observableCount)
        .def("measurements", &Signals::measurements,
             py::return_value_policy::reference_internal)
        .def("__repr__", [](const Signals &self) {
            return fmt::format(
                "<Higgs.Signals with {} observables in {} measurements>",
                self.observableCount(), self.measurements().size());
        });
}

} // namespace Higgs::signals
