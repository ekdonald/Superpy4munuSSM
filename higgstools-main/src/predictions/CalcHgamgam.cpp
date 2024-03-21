#include "predictions/CalcHgamgam.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "utilities/Logging.hpp"
#include <complex>

namespace Higgs::predictions {

namespace {

// helper loop function
inline std::complex<double> fLoop(double tau) {
    if (tau <= 1) {
        return pow(std::asin(sqrt(tau)), 2);
    } else {
        return -1 / 4. *
               pow(std::log((1 + sqrt(1 - 1 / tau)) / (1 - sqrt(1 - 1 / tau))) -
                       std::complex<double>(0, 1) * constants::pi,
                   2);
    }
}

// fermion loop function for CP-even scalar
inline std::complex<double> H12(double tau) {
    return ((tau - 1) * fLoop(tau) + tau) / pow(tau, 2);
}

// fermion loop function for CP-odd scalar
inline std::complex<double> A12(double tau) { return fLoop(tau) / tau; }

// boson loop function for CP-even scalar
inline std::complex<double> H1(double tau) {
    return (3 * (2 * tau - 1) * fLoop(tau) + 3 * tau + 2 * pow(tau, 2)) /
           (2 * pow(tau, 2));
}
} // namespace

double kgamma2(double mass, const NeutralEffectiveCouplings &coups) noexcept {
    auto tauF = [](double m, double mf) {
        return pow(m, 2) / (4 * pow(mf, 2));
    };
    double taut = tauF(mass, constants::mTop);
    double taub = tauF(mass, constants::mBot);
    double tauc = tauF(mass, constants::mCharm);
    double taus = tauF(mass, constants::mStrange);
    double tauu = tauF(mass, constants::mUp);
    double taud = tauF(mass, constants::mDown);
    double tautau = tauF(mass, constants::mTau);
    double taumu = tauF(mass, constants::mMu);
    double taue = tauF(mass, constants::mEl);
    double tauW = tauF(mass, constants::mW);

    double BSMcontr =
        pow(abs((coups.tt.real() * H12(taut) + coups.cc.real() * H12(tauc) +
                 coups.uu.real() * H12(tauu)) *
                    (4 / 3.) +
                (coups.bb.real() * H12(taub) + coups.ss.real() * H12(taus) +
                 coups.dd.real() * H12(taud)) *
                    (1 / 3.) +
                coups.tautau.real() * H12(tautau) +
                coups.mumu.real() * H12(taumu) + coups.ee.real() * H12(taue) -
                coups.WW * H1(tauW)),
            2) +
        pow(abs((coups.tt.imag() * A12(taut) + coups.cc.imag() * A12(tauc) +
                 coups.uu.imag() * A12(tauu)) *
                    (4 / 3.) +
                (coups.bb.imag() * A12(taub) + coups.ss.imag() * A12(taus) +
                 coups.dd.imag() * A12(taud)) *
                    (1 / 3.) +
                coups.mumu.imag() * A12(taumu) +
                coups.mumu.imag() * A12(taumu) + coups.ee.imag() * A12(taue)),
            2);
    double SMcontr = pow(abs((H12(taut) + H12(tauc) + H12(tauu)) * (4 / 3.) +
                             (H12(taub) + H12(taus) + H12(taud)) * (1 / 3.) +
                             H12(tautau) + H12(taumu) + H12(taue) - H1(tauW)),
                         2);

    return BSMcontr / SMcontr;
}

} // namespace Higgs::predictions
