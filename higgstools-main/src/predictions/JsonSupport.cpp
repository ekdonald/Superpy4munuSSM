#include "predictions/JsonSupport.hpp"
#include "utilities/Json.hpp"
#include <magic_enum.hpp> // IWYU pragma: keep
#include <stdexcept>
#include <vector>

namespace Higgs::predictions {
HIGGSUTILITIES_ENUM_JSON_CONV(Collider)
HIGGSUTILITIES_ENUM_JSON_CONV(ECharge)
HIGGSUTILITIES_ENUM_JSON_CONV(Experiment)
HIGGSUTILITIES_ENUM_JSON_CONV(Production)
HIGGSUTILITIES_ENUM_JSON_CONV(Decay)
HIGGSUTILITIES_ENUM_JSON_CONV(ChainDecay)
HIGGSUTILITIES_ENUM_JSON_CONV(ReferenceModel)
HIGGSUTILITIES_ENUM_JSON_CONV(CP)
HIGGSUTILITIES_ENUM_JSON_CONV(Coupling)

void from_json(const nlohmann::json &j, MassResolution &m) {
    j.at("absolute").get_to(m.absolute);
    j.at("relative").get_to(m.relative);
}

void to_json(nlohmann::json &j, const MassResolution &m) {
    j["absolute"] = m.absolute;
    j["relative"] = m.relative;
}

ChannelProcess readChannelProcess(const nlohmann::json &json,
                                  Collider collider) {
    return ChannelProcess{
        collider, utilities::readAs<std::vector<Channel>>(json, "channels")};
}

ChainDecayProcess readChainDecayProcess(const nlohmann::json &j,
                                        Collider collider) {
    return ChainDecayProcess{
        utilities::readAs<ChainDecay>(j, "chain"),
        utilities::readAs<std::vector<Production>>(j, "production"),
        utilities::readAs<std::vector<Decay>>(j, "decay"), collider};
}

PairDecayProcess readPairDecayProcess(const nlohmann::json &j,
                                      Collider collider) {
    return PairDecayProcess{
        utilities::readAs<std::vector<Production>>(j, "production"),
        utilities::readAs<std::vector<Decay>>(j, "firstDecay"),
        utilities::readAs<std::vector<Decay>>(j, "secondDecay"), collider};
}

PairProductionProcess readPairProductionProcess(const nlohmann::json &j,
                                                Collider collider) {
    return PairProductionProcess{
        utilities::readAs<std::vector<Decay>>(j, "firstDecay"),
        utilities::readAs<std::vector<Decay>>(j, "secondDecay"), collider};
}
} // namespace Higgs::predictions
