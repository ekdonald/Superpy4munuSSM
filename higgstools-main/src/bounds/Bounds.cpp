#include "Higgs/Bounds.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "utilities/Json.hpp"
#include "utilities/Logging.hpp"
#include <exception>
#include <filesystem>
#include <ostream>
#include <iostream>
#include <fstream>
#include <range/v3/action/action.hpp>
#include <range/v3/action/join.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/dangling.hpp>
#include <range/v3/range_fwd.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

namespace fs = std::filesystem;



namespace Higgs {


namespace {

std::vector<std::shared_ptr<bounds::Limit>>
loadLimits(std::string_view dataPath,
           const bounds::LimitOptions &limitOptions) {
    auto log = bounds::logger();
    auto limits = std::vector<std::shared_ptr<bounds::Limit>>();
    for (const auto &p : fs::recursive_directory_iterator(dataPath)) {
        if (fs::is_regular_file(p) && p.path().extension() == ".json") {
            const auto filePath = p.path().lexically_normal().string();
            log->trace("reading limit from {}", filePath);

            try {
                auto limit = bounds::Limit::read(filePath, limitOptions);
                const auto duplicateId = [&limit](const auto &l) {
                    return l->id() == limit->id();
                };
                auto found = ranges::find_if(limits, duplicateId);
                if (found == limits.end()) {
                    limits.emplace_back(std::move(limit));
                } else {
                    log->warn(
                        "Duplicate limit id {} for files 1: {} and 2: {}, "
                        "skipping 2.",
                        limit->id(), (*found)->loadedFrom(), filePath);
                }
            } catch (const nlohmann::json::parse_error &) {
                log->warn("Skipping invalid json file {}", filePath);
            } catch (const std::exception &e) {
                log->warn("Skipping {} after read error: {}", filePath,
                          e.what());
            }
        }
    }
    return limits;
}

} // namespace

Bounds::Bounds(const std::string &dataPath,
               const bounds::LimitOptions &limitOptions)
    : limits_{loadLimits(dataPath, limitOptions)} {
    if (limits_.empty()) {
        throw std::invalid_argument(fmt::format(
            "No valid limit files found in data path {}", dataPath));
    }
}

const std::vector<std::shared_ptr<bounds::Limit>> &
Bounds::limits() const noexcept {
    return limits_;
}

namespace {
auto contains(const std::string &p) {
    return [&p](const bounds::AppliedLimit &appLim) {
        return ranges::find(appLim.contributingParticles(), p) !=
               appLim.contributingParticles().end();
    };
}
} // namespace

bounds::HBResult
Bounds::operator()(const Higgs::Predictions &predictions) const {
    auto res = bounds::HBResult{};
    res.appliedLimits = applyLimits(predictions);
    for (const auto &p : predictions.particleIds()) {
        auto sensitiveLimit = ranges::find_if(res.appliedLimits, contains(p));
        if (sensitiveLimit != res.appliedLimits.end()) {
            res.selectedLimits[p] = *sensitiveLimit;
        }
    }

    auto areAllowed = [](const auto &lim) { return lim.obsRatio() < 1.; };
    res.allowed =
        ranges::all_of(ranges::views::values(res.selectedLimits), areAllowed);
    return res;
}

std::vector<bounds::AppliedLimit>
Bounds::applyLimits(const Higgs::Predictions &predictions) const noexcept {
    auto applyLimit = [&predictions](const auto &l) {
        return l->apply(predictions);
    };
    auto byExpRatio = [](const auto &a, const auto &b) {
        return a.expRatio() > b.expRatio();
    };
    return limits_ | ranges::views::transform(applyLimit) |
           ranges::actions::join | ranges::actions::sort(byExpRatio);
}

//std::ostream &bounds::operator<<(std::ostream &os,
//                                 const bounds::HBResult &res) {
//    os << "HiggsBounds result: " << (res.allowed ? "allowed" : "excluded")
//       << "\n"
//          "\tparticle | obsRatio | expRatio | selected limit description\n"
//          "\t---------|----------|----------|---------------------------\n";
//    for (const auto &[p, lim] : res.selectedLimits) {
//        os << fmt::format("\t {: ^7} | {: ^ 8.3f} | {: ^ 8.3f} | {}\n", p,
//                          lim.obsRatio(), lim.expRatio(),
//                          lim.limit()->to_string());
//    }
//    return os;
//    
// }



// By DEK  std::ofstream  --- on 30/01/2024

  std::ofstream MyFile("hboutput.dat");
  std::ostream & bounds::operator<<(std::ostream &os,
                                 const bounds::HBResult &res) {
  // write in file
  MyFile << "# particle  obsRatio  expRatio  status\n";
          // "\t--------------------------------------------------------\n";    
  for (const auto &[p, lim] : res.selectedLimits) {
         MyFile << fmt::format(" {: ^7}  {: ^ 8.3f}  {: ^ 8.3f}  {}\n", p,
                          lim.obsRatio(), lim.expRatio(),
                           (res.allowed ? "allowed" : "excluded") );
    }
  MyFile.close();
  
  
  // write the same in std
  os << "HiggsBounds result: " << (res.allowed ? "allowed" : "excluded")
         << "\n"
          "\tparticle | obsRatio | expRatio | selected limit description\n"
          "\t---------|----------|----------|---------------------------\n";    
  for (const auto &[p, lim] : res.selectedLimits) {
         os << fmt::format("\t {: ^7} | {: ^ 8.3f} | {: ^ 8.3f} | {}\n", p,
                          lim.obsRatio(), lim.expRatio(),
                          lim.limit()->to_string());
    }
    
  return os;
  }



} // namespace Higgs








