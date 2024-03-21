#include "Higgs/predictions/Basics.hpp"
#include "predictions/FormatSupport.hpp" // IWYU pragma: keep
#include "utilities/Logging.hpp"
#include <magic_enum.hpp>
#include <memory>

namespace Higgs::predictions {

ColliderType classifyCollider(Collider coll) noexcept {
    static_assert(magic_enum::enum_count<Collider>() == 3,
                  "You have changed Higgs::predictions::Collider and need to "
                  "update Higgs::predictions::classifyCollider as well!");
    if (coll == Collider::LEP) {
        return ColliderType::ee;
    }
    return ColliderType::pp;
}

InvalidChannel::InvalidChannel(ECharge charge, Decay d) noexcept
    : std::runtime_error{fmt::format(
          "Invalid decay mode {} for charge {} particle", d, charge)} {
    logger()->error(what());
}

InvalidChannel::InvalidChannel(ECharge charge, Production p) noexcept
    : std::runtime_error{fmt::format(
          "Invalid production mode {} for charge {} particle", p, charge)} {
    logger()->error(what());
}

InvalidChannel::InvalidChannel(ColliderType collType, Production p) noexcept
    : std::runtime_error{fmt::format(
          "Invalid production mode {} at {}-collider", p, collType)} {
    logger()->error(what());
}

InvalidChannel::InvalidChannel(Production p, Decay d) noexcept
    : std::runtime_error{
          fmt::format("Inconsistent charge between production mode {} "
                      "and decay mode {}",
                      p, d)} {
    logger()->error(what());
}

} // namespace Higgs::predictions
