#pragma once

#include "Higgs/Predictions.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "bounds/limits/BasicLimit.hpp"
#include <vector>

namespace Higgs {
namespace predictions {
class Predictions;
}

namespace bounds {
class TestLimit : public Higgs::bounds::BasicLimit {
  public:
    template <typename... Args> static auto create(Args &&...args) {
        return std::shared_ptr<Limit>(
            new TestLimit(std::forward<Args>(args)...));
    }

    std::vector<Higgs::bounds::AppliedLimit>
    apply(const Higgs::predictions::Predictions &prediction) const override {
        return {};
    }

    std::string processDesc() const noexcept override { return ""; }
    std::string extentDesc() const noexcept override { return ""; }

    template <typename... Args>
    explicit TestLimit(Args &&...args)
        : BasicLimit{std::forward<Args>(args)..., LimitOptions{}} {}

    explicit TestLimit(LimitOptions options)
        : BasicLimit{0, "", "", {}, {}, 0, {}, options} {}
};
} // namespace bounds
} // namespace Higgs
