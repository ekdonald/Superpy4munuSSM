#include "Higgs/Predictions.hpp"
namespace HP = Higgs::predictions;

int main() {
    auto h = HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::even);

    h.setTotalWidth(0.3);
    h.setBr(HP::Decay::bb, 0.8);
    h.setBr(HP::Decay::tautau, 0.17);
    h.setBr(HP::Decay::gamgam, 3e-2);
    // BR(h->bb)=80%, BR(h->tautau)=17%, BR(h->gamgam)=3%

    try {
        h.setBr(HP::Decay::WW, 0.11);
    } catch (const HP::InvalidInput &) {
        // this will be an error, since sum(BRs) is already 100%
    }
}
