#include "Higgs/Bounds.hpp"
#include "Higgs/Predictions.hpp"
#include <iostream>

int main() {
    auto bounds = Higgs::Bounds{"/path/to/HBDataSet"};

    auto pred = Higgs::Predictions{};
    // set model predictions on pred

    auto hbResult = bounds(pred);
    std::cout << hbResult << std::endl;
}
