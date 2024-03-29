#include "Higgs/Predictions.hpp"
namespace HP = Higgs::predictions;

int main() {
    auto h = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::even);

    h.setTotalWidth(0.1);
    h.setBr(HP::Decay::bb, 1.);
    // BR(h->bb)=100%
    h.setDecayWidth(HP::Decay::tautau, 0.2);
    // totalWidth() = 1.2, BR(h->bb)~83.4%, BR(h->tautau)~16.7%
    h.setDecayWidth(HP::Decay::gamgam, 0.01);
    // totalWidth() = 1.21, BR(h->bb)~82.7%, BR(h->tautau)~16.5%,
    // BR(h->gamgam)~0.8%
    h.setDecayWidth(HP::Decay::tautau, 0.4);
    // totalWidth() = 1.41, BR(h->bb)~70.9%, BR(h->tautau)~28.4%,
    // BR(h->gamgam)~0.7%
}
