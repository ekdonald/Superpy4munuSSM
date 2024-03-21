#include "Higgs/Bounds.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/Signals.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace HP = Higgs::predictions;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0]
                  << " /path/to/HBDataSet /path/to/HSDataSet" << std::endl;
        return 1;
    }

    auto pred = Higgs::Predictions{};
    auto &h = pred.addParticle(
        HP::BsmParticle{"h", HP::ECharge::neutral, HP::CP::even});
    h.setMass(125.09);
    effectiveCouplingInput(h, HP::smLikeEffCouplings,
                           HP::ReferenceModel::SMHiggsEW);

    const auto bounds = Higgs::Bounds(argv[1]);
    const auto resultHB = bounds(pred);
    std::cout << resultHB << std::endl;
    std::cout << "All applied limits: obsRatio (expRatio)\n";
    for (const auto &al : resultHB.appliedLimits) {
        std::cout << al.limit()->id() << " " << al.limit()->processDesc()
                  << ": " << al.obsRatio() << " (" << al.expRatio() << ")"
                  << std::endl;
    }

    const auto signals = Higgs::Signals(argv[2]);
    auto resultHS = signals(pred);
    std::cout << "\n HiggsSignals chisq: " << resultHS << " from "
              << signals.observableCount() << " observables" << std::endl;
}
