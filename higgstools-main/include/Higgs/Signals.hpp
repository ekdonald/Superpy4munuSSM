/**
 * @file Signals.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Interface to the HiggsSignals library. Includes all other headers
 * needed for all of the HiggsSignals functionality.
 *
 * @copyright Copyright 2022 by the authors. This file is part of HiggsTools.
 * HiggsTools is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/Predictions.hpp"
#include "Higgs/signals/Measurement.hpp"
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace Higgs {
namespace signals {

struct GlobalCorrelations;

//! Main class of the HiggsSignals library. This class loads and stores all
//! Measurement%s contained in the HiggsSignals dataset and can be applied to a
//! Higgs::predictions::Prediction to obtain the HiggsSignals \f$\chi^2\f$ of
//! the model predictions.
class HIGGSTOOLS_EXPORT Signals {
  public:
    /**
     * @brief Construct a new Signals object by loading all available
     * Measurements.
     *
     * This reads in all of the Measurements from disk. Make sure to reuse the
     * constructed object whenever possible. If reading fails for any of the
     * json files, those will be skipped and a warning will be logged.
     *
     * @throws std::invalid_argument if no valid measurements are found in the
     * dataPath.
     * @param dataPath The filesystem path to read from. All `.json` files
     * within this folder and its subfolders are read.
     * @param measurementOptions options to pass to all read measurements
     */
    Signals(const std::string &dataPath,
            const MeasurementOptions &measurementOptions = {});

    //! Global modification factors for the entire HiggsSignals dataset. The
    //! numeric key of the map specifies the Measurement::id() of the
    //! measurement that the Measurement::ModificationFactors should be passed
    //! to.
    using ModificationFactors =
        std::unordered_map<std::size_t, Measurement::ModificationFactors>;

    /**
     * @brief Compute the \f$\chi^2\f$ of the given model prediction compared to
     * all of the measurements in the HiggsSignals dataset.
     *
     * This is the main function of HiggsSignals.
     *
     * @param predictions the model predictions to compare to the measurements
     * @param modificationFactors optionally specify modification factors that
     * can separately rescale the model predictions for any Measurement or
     * SubMeasurement and signal process channel. Unspecified modification
     * factors mean that no scaling is performed for the corresponding
     * Measurement or SubMeasurement.
     * @return double the total \f$\chi^2\f$ value
     */
    double
    operator()(const Higgs::predictions::Predictions &predictions,
               const ModificationFactors &modificationFactors = {}) const;

    //! The total number of observables in the loaded measurements
    std::size_t observableCount() const noexcept;

    //! Access the list of loaded measurements.
    const std::vector<Measurement> &measurements() const noexcept;

  private:
    std::vector<Measurement> measurements_;
    std::size_t nSubMeasurements_;
    std::shared_ptr<GlobalCorrelations> corr_;
};
} // namespace signals

using signals::Signals;
} // namespace Higgs
