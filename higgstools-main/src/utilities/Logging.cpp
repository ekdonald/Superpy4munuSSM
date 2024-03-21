#include "Logging.hpp"
#include <spdlog/cfg/env.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <string>

namespace Higgs {

namespace {

bool initializeLogger(const std::string &name) {
    auto logger = spdlog::stdout_color_mt(name);
    spdlog::cfg::load_env_levels();
    return logger != nullptr;
}

} // namespace

std::shared_ptr<spdlog::logger> predictions::logger() {
    constexpr auto name = "HiggsPredictions";
    [[maybe_unused]] static const auto initialized = initializeLogger(name);
    return spdlog::get(name);
}
std::shared_ptr<spdlog::logger> bounds::logger() {
    constexpr auto name = "HiggsBounds";
    [[maybe_unused]] static const auto initialized = initializeLogger(name);
    return spdlog::get(name);
}
std::shared_ptr<spdlog::logger> signals::logger() {
    constexpr auto name = "HiggsSignals";
    [[maybe_unused]] static const auto initialized = initializeLogger(name);
    return spdlog::get(name);
}

} // namespace Higgs
