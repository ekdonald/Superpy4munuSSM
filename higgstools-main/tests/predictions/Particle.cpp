#include "Higgs/predictions/Particle.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <string>

using Catch::Approx;
using namespace std::string_literals;
namespace HP = Higgs::predictions;

TEST_CASE("Particle value semantics") {

    auto h1 = HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::undefined);
    h1.setMass(100);
    auto h2 = HP::BsmParticle("h2", HP::ECharge::single, HP::CP::undefined);
    h2.setMass(200);
    SECTION("copying") {
        auto s1{h1};
        auto s2(h2);
        CHECK(s1.id() == h1.id());
        CHECK(s2.id() == h2.id());
        s1 = s2;
        CHECK(s1.id() == s2.id());
        CHECK(s1.id() == h2.id());
        s2 = h1;
        CHECK(s2.id() == h1.id());
    }

    SECTION("moves") {
        auto s1{h1};
        auto s2(h2);
        auto s3{std::move(s1)};
        CHECK(s3.id() == h1.id());
        s1 = std::move(s2);
        CHECK(s1.id() == h2.id());
    }

    SECTION("copying deep copies") {
        auto s1 = h1;
        s1.setTotalWidth(100.);
        auto s2 = s1;
        REQUIRE(s2.totalWidth() == Approx(100.));
        s1.setTotalWidth(50.);
        REQUIRE(s1.totalWidth() == Approx(50.));
        REQUIRE(s2.totalWidth() == Approx(100.));
    }
}

class DummyParticle : public HP::Particle {
  public:
    DummyParticle()
        : Particle{"dummy", HP::CP::undefined, HP::ECharge::neutral} {}
    std::unique_ptr<Particle> clone() const {
        return std::make_unique<DummyParticle>(*this);
    }
};

TEST_CASE("Particle interface") {

    auto n = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::undefined);
    n.setMass(100);
    auto d = DummyParticle();

    SECTION("everything callable and zero") {
        for (const HP::Particle &x :
             {std::cref<HP::Particle>(n), std::cref<HP::Particle>(d)}) {
            CHECK(x.cxn(HP::Collider{}, HP::Production{}) == 0);
            CHECK(x.br(HP::Decay{}) == 0.);
            CHECK(x.br(HP::ChainDecay{}, "") == 0.);
            CHECK(x.br("", "") == 0.);
            CHECK(x.channelRate(HP::Collider{}, HP::Production{},
                                HP::Decay{}) == 0.);
            CHECK(x.totalWidth() == 0.);
            CHECK_FALSE(x.coupling(HP::Coupling{}).has_value());
        }
    }
}
