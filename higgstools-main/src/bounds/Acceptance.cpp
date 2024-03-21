#include "bounds/Acceptance.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "predictions/JsonSupport.hpp"
#include "utilities/Format.hpp"
#include "utilities/Json.hpp"
#include <algorithm>
#include <map>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/utility/get.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <stdexcept>
#include <utility>

namespace Higgs::bounds {

ConstantAcceptance::ConstantAcceptance(const nlohmann::json &j)
    : ConstantAcceptance{utilities::readAs<double>(j, "constantAcceptance")} {}

MassDepAcceptance::MassDepAcceptance(const std::vector<double> &values,
                                     std::vector<double> massGrid)
    : interp_{{std::move(massGrid)}, values} {}

MassDepAcceptance::MassDepAcceptance(const nlohmann::json &j)
    : MassDepAcceptance::MassDepAcceptance{
          utilities::readAs<std::vector<double>>(j, "massDepAcceptance"),
          utilities::readAs<std::vector<double>>(j, "massGrid")} {}

double MassDepAcceptance::operator()(const predictions::Particle &p) const {
    return interp_({p.mass()});
}

CouplingDepAcceptance::CouplingDepAcceptance(
    CouplingDepAcceptance::Polynomial numerator,
    CouplingDepAcceptance::Polynomial denominator)
    : numerator_{std::move(numerator)}, denominator_{std::move(denominator)} {}

namespace {

auto readMonomial(const nlohmann::json &jmonomial) {
    // the json.items() proxy does not work well with ranges
    auto result = CouplingDepAcceptance::Monomial{};
    for (auto [key, val] : jmonomial.items()) {
        auto coup = magic_enum::enum_cast<predictions::Coupling>(key);
        if (!coup) {
            throw utilities::BadEnumRead("Invalid value \"" + key +
                                         "\" for Higgs::predictions::Coupling");
        }
        result.emplace(coup.value(), val.get<int>());
    }
    return result;
}

constexpr auto readTerm = [](const nlohmann::json &jterm) {
    return CouplingDepAcceptance::Term{readMonomial(jterm.at(0)),
                                       readAcceptance(jterm.at(1))};
}; // LCOV_EXCL_LINE

auto readPolynomialAcceptance(const nlohmann::json &j) {
    return j | ranges::views::transform(readTerm) |
           ranges::to<CouplingDepAcceptance::Polynomial>;
}

} // namespace

CouplingDepAcceptance::CouplingDepAcceptance(const nlohmann::json &j)
    : numerator_{readPolynomialAcceptance(j.at("couplingDepAcceptance"))},
      denominator_{j.contains("denominator")
                       ? readPolynomialAcceptance(j.at("denominator"))
                       : Polynomial{}} {}

double CouplingDepAcceptance::operator()(const predictions::Particle &p) const {
    auto evaluateAcceptanceTerm = [&p](const Term &term) {
        auto &[monomial, acc] = term;
        auto evaluateAcceptance = [&p](const auto &a) { return a(p); };
        return evaluateCouplingMonomial(monomial, p) *
               std::visit(evaluateAcceptance, acc);
    };
    auto result =
        ranges::accumulate(numerator_, 0., std::plus{}, evaluateAcceptanceTerm);
    if (!denominator_.empty() && result > 0) {
        result /= ranges::accumulate(denominator_, 0., std::plus{},
                                     evaluateAcceptanceTerm);
    }
    return result;
}

double CouplingDepAcceptance::evaluateCouplingMonomial(
    const CouplingDepAcceptance::Monomial &monomial,
    const predictions::Particle &p) {
    auto evaluateCouplingPower = [&p](const auto &couplingPower) {
        return std::pow(p.coupling(couplingPower.first).value_or(0.),
                        couplingPower.second);
    };
    return ranges::accumulate(monomial, 1., std::multiplies{},
                              evaluateCouplingPower);
}

Acceptances::Acceptances(std::vector<Acceptance> acceptances,
                         std::size_t nChannels)
    : acceptances_{std::move(acceptances)} {
    if (acceptances_.size() != nChannels) {
        throw std::out_of_range(fmt::format(
            "Number of acceptances {} does not match number of channels {}",
            acceptances_.size(), nChannels));
    }
}

Acceptances::Acceptances(const nlohmann::json &j, std::size_t nChannels)
    : Acceptances{j | ranges::views::transform(readAcceptance) |
                      ranges::to<std::vector>,
                  nChannels} {}

std::vector<double>
Acceptances::operator()(const predictions::Particle &p) const {
    auto evaluateAcceptance = [&p](const auto &a) { return a(p); };
    return acceptances_ |
           ranges::views::transform([&evaluateAcceptance](const auto &acc) {
               return std::visit(evaluateAcceptance, acc);
           }) |
           ranges::to<std::vector>;
}

Acceptance readAcceptance(const nlohmann::json &j) {
    if (j.contains("constantAcceptance")) {
        return ConstantAcceptance{j};
    } else if (j.contains("massDepAcceptance")) {
        return MassDepAcceptance{j};
    } else if (j.contains("couplingDepAcceptance")) {
        return CouplingDepAcceptance{j};
    }
    throw utilities::BadFieldRead("Unknown acceptance in json.");
}

} // namespace Higgs::bounds
