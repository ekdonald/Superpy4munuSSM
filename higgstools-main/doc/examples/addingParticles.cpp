#include "Higgs/Predictions.hpp"
namespace HP = Higgs::predictions;

int main() {
    auto pred = Higgs::Predictions();
    auto &h1 = pred.addParticle(
        HP::BsmParticle("h1", HP::ECharge::neutral, HP::CP::even));
    auto &hc = pred.addParticle(HP::BsmParticle("H+", HP::ECharge::single, HP::CP::undefined));
    // set properties on h1 and hc
}
