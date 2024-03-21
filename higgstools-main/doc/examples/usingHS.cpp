#include "Higgs/Predictions.hpp"
#include "Higgs/Signals.hpp"
#include <iostream>

int main() {
    auto signals = Higgs::Signals{"/path/to/HSDataSet"};

    auto pred = Higgs::Predictions{};
    // set model predictions on pred

    auto hsChisq = signals(pred);
    std::cout << "HS chisq = " << hsChisq << std::endl;
}
