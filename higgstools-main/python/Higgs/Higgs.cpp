#include "Higgs.hpp"

PYBIND11_MODULE(_Higgs, h) {
    // TODO: h.doc() =

    Higgs::predictions::bindHiggsPredictions(h);
    Higgs::bounds::bindHiggsBounds(h);
    Higgs::signals::bindHiggsSignals(h);
}
