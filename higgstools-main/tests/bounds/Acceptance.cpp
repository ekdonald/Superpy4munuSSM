#include "bounds/Acceptance.hpp"
#include "Higgs/predictions/Channels.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "predictions/JsonSupport.hpp"
#include "utilities/Json.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <string>

using Catch::Approx;
namespace HB = Higgs::bounds;
namespace HP = Higgs::predictions;
namespace HU = Higgs::utilities;

TEST_CASE("Acceptance -- constant") {
    const auto acc1 = HB::ConstantAcceptance();
    const auto acc2 = HB::ConstantAcceptance(2);

    const auto p1 = HP::BsmParticle("p1", HP::ECharge::neutral, HP::CP::undefined);
    const auto p2 = HP::SMHiggs(200);

    CHECK(acc1(p1) == 1.);
    CHECK(acc2(p1) == 2.);
    CHECK(acc1(p2) == 1.);
    CHECK(acc2(p2) == 2.);

    SECTION("json read") {
        const auto j = nlohmann::json::parse(R"(
            {
                "constantAcceptance": 1.5
            }
        )");
        auto accj = HB::ConstantAcceptance(j);
        CHECK(accj(p1) == 1.5);
        CHECK(accj(p2) == 1.5);
        auto accj2 = HB::readAcceptance(j);
        REQUIRE(std::holds_alternative<HB::ConstantAcceptance>(accj2));
        CHECK(std::get<HB::ConstantAcceptance>(accj2)(p1) == 1.5);
        CHECK(std::get<HB::ConstantAcceptance>(accj2)(p2) == 1.5);
    }
}

TEST_CASE("Acceptance -- mass dependent") {
    const auto acc = HB::MassDepAcceptance({1.5, 2.5}, {100, 200});

    const auto p50 = HP::SMHiggs(50);
    const auto p100 = HP::SMHiggs(100);
    const auto p150 = HP::SMHiggs(150);
    const auto p200 = HP::SMHiggs(200);
    const auto p300 = HP::SMHiggs(300);

    CHECK(acc(p50) == acc(p100));
    CHECK(acc(p100) == 1.5);
    CHECK(acc(p150) == Approx(2.));
    CHECK(acc(p200) == 2.5);
    CHECK(acc(p200) == acc(p300));

    REQUIRE_THROWS_AS(HB::MassDepAcceptance({1., 2., 3.}, {100., 200.}),
                      std::out_of_range);

    SECTION("json read") {
        const auto j = nlohmann::json::parse(R"(
            {
                "massDepAcceptance": [0.2, 0.6],
                "massGrid": [100,200]
            }
        )");
        auto accj = HB::MassDepAcceptance(j);
        CHECK(accj(p50) == accj(p100));
        CHECK(accj(p100) == 0.2);
        CHECK(accj(p150) == Approx(0.4));
        CHECK(accj(p200) == 0.6);
        CHECK(accj(p200) == accj(p300));

        auto accj2 = HB::readAcceptance(j);
        REQUIRE(std::holds_alternative<HB::MassDepAcceptance>(accj2));
        CHECK(std::get<HB::MassDepAcceptance>(accj2)(p50) ==
              std::get<HB::MassDepAcceptance>(accj2)(p100));
        CHECK(std::get<HB::MassDepAcceptance>(accj2)(p100) == 0.2);
        CHECK(std::get<HB::MassDepAcceptance>(accj2)(p150) == Approx(0.4));
        CHECK(std::get<HB::MassDepAcceptance>(accj2)(p200) == 0.6);
        CHECK(std::get<HB::MassDepAcceptance>(accj2)(p200) ==
              std::get<HB::MassDepAcceptance>(accj2)(p300));
    }
}

TEST_CASE("Acceptance -- coupling dependent") {
    const auto cte4 =
        HB::CouplingDepAcceptance::Monomial{{HP::Coupling::effCPeTopYuk, 4}};
    const auto cte2cto2 = HB::CouplingDepAcceptance::Monomial{
        {HP::Coupling::effCPeTopYuk, 2}, {HP::Coupling::effCPoTopYuk, 2}};
    const auto cto4 =
        HB::CouplingDepAcceptance::Monomial{{HP::Coupling::effCPoTopYuk, 4}};

    const auto massDepAcc1 = HB::MassDepAcceptance({1.5, 2.5}, {100, 200});

    const auto poly1 = HB::CouplingDepAcceptance::Polynomial{
        {cte4, HB::MassDepAcceptance({1.5, 2.5}, {100, 200})},
        {cte2cto2, HB::MassDepAcceptance({1.3, 0.4, 0.9}, {100, 150, 200})},
        {cto4, HB::ConstantAcceptance(0.4)}};

    const auto poly2 = HB::CouplingDepAcceptance::Polynomial{
        {cte4, HB::MassDepAcceptance({2.2, 1.4}, {100, 200})}};

    const auto poly3 = HB::CouplingDepAcceptance::Polynomial{
        {cto4, HB::CouplingDepAcceptance{poly2}}};

    auto p = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
    p.setMass(125);
    p.setCoupling(HP::Coupling::effCPeTopYuk, 2.);
    p.setCoupling(HP::Coupling::effCPoTopYuk, 0.3);

    SECTION("numerator == denominator") {
        const auto idAcc1 = HB::CouplingDepAcceptance{poly1, poly1};
        const auto idAcc2 = HB::CouplingDepAcceptance{poly2, poly2};
        const auto idAcc3 = HB::CouplingDepAcceptance{poly3, poly3};
        REQUIRE(idAcc1(p) == Approx(1.));
        REQUIRE(idAcc2(p) == Approx(1.));
        REQUIRE(idAcc3(p) == Approx(1.));
    }

    const auto acc1 = HB::CouplingDepAcceptance{poly1};
    const auto acc2 = HB::CouplingDepAcceptance{poly2};
    const auto acc12 = HB::CouplingDepAcceptance{poly1, poly2};
    const auto acc21 = HB::CouplingDepAcceptance{poly2, poly1};
    const auto acc31 = HB::CouplingDepAcceptance{poly3, poly1};

    CHECK(acc1(p) == Approx(1.75 * std::pow(2., 4.) +
                            1.7 / 2. * std::pow(2, 2) * std::pow(0.3, 2) +
                            0.4 * std::pow(0.3, 4)));
    CHECK(acc2(p) == Approx(2. * std::pow(2., 4)));
    CHECK(acc12(p) == Approx(acc1(p) / acc2(p)));
    CHECK(acc21(p) == Approx(acc2(p) / acc1(p)));
    CHECK(acc31(p) == Approx(std::pow(0.3, 4) * std::pow(2, 4) * 2. / acc1(p)));

    SECTION("some edge cases with zeros") {
        const auto zeroPoly = HB::CouplingDepAcceptance::Polynomial{
            {cte4, HB::ConstantAcceptance{0.}}};
        CHECK(HB::CouplingDepAcceptance{zeroPoly, poly2}(p) == 0.);
        CHECK(HB::CouplingDepAcceptance{zeroPoly, poly1}(p) == 0.);
        CHECK(HB::CouplingDepAcceptance{zeroPoly}(p) == 0.);

        CHECK(HB::CouplingDepAcceptance{zeroPoly, zeroPoly}(p) == 0.);
    }

    SECTION("json read") {
        auto j = nlohmann::json::parse(R"(
        {
            "couplingDepAcceptance": [
                [
                    {"effCPeTopYuk": 2, "effCPoTopYuk": 2},
                    {"constantAcceptance": 3.1}
                ],
                [
                    {"effCPeTopYuk": 4},
                    {
                        "massDepAcceptance": [0.3,0.7,0.9],
                        "massGrid": [100,200,500]
                    }
                ],
                [
                    {"effCPoTopYuk": 4},
                    {
                        "couplingDepAcceptance": [
                            [{"alphaCPTauYuk": 4}, {"constantAcceptance": 1.1}]
                        ]
                    }
                ]
            ]
        })");

        auto readAcc1 = HB::CouplingDepAcceptance(j);
        auto readAcc2 =
            std::get<HB::CouplingDepAcceptance>(HB::readAcceptance(j));
        REQUIRE(readAcc1(p) == Approx(0.4 * std::pow(2., 4.) +
                                      3.1 * std::pow(0.3, 2) * std::pow(2, 2)));
        REQUIRE(readAcc1(p) == readAcc2(p));
        j["denominator"] = nlohmann::json::parse(R"([
                [{"effCPeTopYuk": 4}, {"constantAcceptance": 1.1}]
            ])");
        auto readAcc3 = HB::CouplingDepAcceptance(j);
        REQUIRE(readAcc3(p) == Approx(readAcc1(p) / (1.1 * std::pow(2, 4))));

        SECTION("invalid") {
            j["denominator"] = nlohmann::json::parse(R"([
                [{"notacoupling": 4}, {"constantAcceptance": 1.1}]
            ])");
            REQUIRE_THROWS_AS(HB::CouplingDepAcceptance(j), HU::BadEnumRead);
        }
    }
}

TEST_CASE("Acceptances") {
    const auto accConst = HB::ConstantAcceptance();
    const auto accMass = HB::MassDepAcceptance({1.5, 2.5}, {100, 200});

    REQUIRE_THROWS_AS(HB::Acceptances({accConst, accMass}, 3),
                      std::out_of_range);
    auto accs = HB::Acceptances{{accConst, accMass}, 2};

    const auto p = HP::SMHiggs(170);
    auto res = accs(p);
    CHECK(res[0] == accConst(p));
    CHECK(res[1] == accMass(p));

    SECTION("json read") {
        const auto j = nlohmann::json::parse(R"([
            {
                "constantAcceptance": 1.5
            },
            {
                "massDepAcceptance": [0.2, 0.6],
                "massGrid": [100,200]
            }
        ])");
        REQUIRE_THROWS_AS(HB::Acceptances(j, 1), std::out_of_range);
        auto accsj = HB::Acceptances(j, 2);
        const auto p = HP::SMHiggs(110);
        auto resj = accsj(p);
        CHECK(resj[0] == HB::ConstantAcceptance(j[0])(p));
        CHECK(resj[1] == HB::MassDepAcceptance(j[1])(p));

        const auto jInvalid = nlohmann::json::parse(R"([
            {"not an acceptance":{}}
        ])");
        REQUIRE_THROWS_AS(HB::Acceptances(jInvalid, 1), HU::BadFieldRead);
    }
}
