/**
 * @file Bounds.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * Interface to the HiggsBounds library. Includes all other headers
 * needed for access to all of HiggsBounds functionality.
 *
 * @copyright Copyright 2021 by the authors. This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/bounds/Limit.hpp" // IWYU pragma: export
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>
namespace Higgs {
namespace predictions {
class Predictions;
} // namespace predictions

namespace bounds {

//! The result of running HiggsBounds.
struct HIGGSTOOLS_EXPORT HBResult {
    //! The overall verdict whether the point is allowed at 95% C.L..
    bool allowed;
    //! The selected most sensitive limits for each particle. The overall result
    //! is `allowed` if none of these exclude the model predictions. A limit
    //! excludes the model predictions if AppliedLimit::obsRatio > 1. The
    //! sensitivity is judged based on the AppliedLimit::expRatio().
    std::map<std::string, AppliedLimit> selectedLimits;
    //! All limits that were compared to the model predictions. They are sorted
    //! by decreasing sensitivity, i.e. decreasing AppliedLimit::expRatio().
    //! See also LimitOptions::minExpRatio.
    std::vector<AppliedLimit> appliedLimits;
};

//! Prints a HiggsBounds result. @relates HBResult
HIGGSTOOLS_EXPORT std::ostream &operator<<(std::ostream &os,
                                           const HBResult &res);

//! Main class of the HiggsBounds library. This class loads and stores all
//! available limits from files and can then be applied to
//! #Higgs::predictions::Predictions to compare the model predictions to the
//! available limits.
class HIGGSTOOLS_EXPORT Bounds {
  public:
    /**
     * @brief Construct a new Bounds object by loading all available Limits.
     *
     * This reads in all of the Limits from disk, so this is *pretty slow*.
     * Make sure to reuse the constructed object whenever possible. If reading
     * fails for any of the Limits those will be skipped and a warning will be
     * logged.
     *
     * @throws std::invalid_argument if no valid limits are found in the
     * dataPath.
     * @param dataPath The filesystem path to read from. All ".json" files
     * within this folder and its subfolders are read.
     * @param limitOptions options to pass to all read limits
     */
    Bounds(const std::string &dataPath, const LimitOptions &limitOptions = {});

    /**
     * @brief Check the given predictions against the limits.
     *
     * This is the main function of HiggsBounds. It Limit::apply%'s each of the
     * contained Limit%s to the `predictions`. Out of the resulting
     * AppliedLimit%s the most sensitive limit for each particle in the
     * predictions is selected based on the expected sensitivity. The
     * `predictions` are considered allowed if they are not excluded (based on
     * the observed limit) by any of the selected limits.
     *
     * @param predictions the predictions that are tested against the
     * HiggsBounds limits
     * @return HBResult the result of the HiggsBounds run
     */
    HBResult
    operator()(const Higgs::predictions::Predictions &predictions) const;

    //! Lists all loaded limits.
    const std::vector<std::shared_ptr<Limit>> &limits() const noexcept;

  private:
    std::vector<AppliedLimit>
    applyLimits(const predictions::Predictions &predictions) const noexcept;

    std::vector<std::shared_ptr<Limit>> limits_;
};

} // namespace bounds

using bounds::Bounds;
} // namespace Higgs
