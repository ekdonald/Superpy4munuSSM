#include "Higgs/Bounds.hpp"
#include "Higgs/Predictions.hpp"
#include "Higgs/Signals.hpp"
#include "Higgs/bounds/Limit.hpp"
#include "Higgs/predictions/Basics.hpp"
#include "Higgs/predictions/EffectiveCouplings.hpp"
#include "Higgs/predictions/Particle.hpp"
#include "Higgs/predictions/ReferenceModels.hpp"
#include "utilities/Format.hpp"
#include <exception>
#include <iostream>
#include <magic_enum.hpp>
#include <memory>
#include <string>
#include <vector>
#include <wstp.h>

namespace {

////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////

Higgs::Predictions pred{};

std::unique_ptr<Higgs::Bounds> bounds{};
Higgs::bounds::LimitOptions HBoptions{};
Higgs::bounds::HBResult resultHB;

std::unique_ptr<Higgs::Signals> signals{};
Higgs::signals::MeasurementOptions HSoptions{};

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////

int putReal(double d) { return WSPutReal(stdlink, d); }

int putString(const char *s) { return WSPutString(stdlink, s); }

int putFunction(const char *s, int argc) {
    return WSPutFunction(stdlink, s, argc);
}

int putSymbol(const char *s) { return WSPutSymbol(stdlink, s); }

int putInteger(int i) { return WSPutInteger(stdlink, i); }

void putRule(const char *s) {
    putFunction("Rule", 2);
    putSymbol(s);
}

void putSymbolRule(const char *s, const char *i) {
    putRule(s);
    putSymbol(i);
}

void putIntRule(const char *s, int i) {
    putRule(s);
    putInteger(i);
}

void putRealRule(const char *s, double i) {
    putRule(s);
    putReal(i);
}

void putStringRule(const char *s, const char *i) {
    putRule(s);
    putString(i);
}

void putListRule(const char *s, int i) {
    putRule(s);
    putFunction("List", i);
}

void putError(const std::string &error) {
    putFunction("CompoundExpression", 2);
    putFunction("ToExpression", 1);
    putString(("Print[\"HB error: " + error + "\"]").c_str());
    putSymbol("$Failed");
}

inline const char *const BoolToString(bool b) { return b ? "True" : "False"; }

bool StringToBool(std::string const &s) {
    if (s == "True") {
        return true;
    } else if (s == "False") {
        return false;
    } else {
        putError("Can not cast " + s + " to string. Please enter either 'True' or 'False'.");
    }
    return s == "True"; 
}

template <class Enum> Enum enumFromName(const std::string_view enumName) {
    auto opt = magic_enum::enum_cast<Enum>(enumName);
    if (opt) {
        return opt.value();
    }
    throw(std::runtime_error(
        fmt::format("{} is not a valid value for {}, options are: [{}]",
                    enumName, magic_enum::enum_type_name<Enum>(),
                    fmt::join(magic_enum::enum_names<Enum>(), ", "))));
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// HiggsPredictions functions
////////////////////////////////////////////////////////////////////////////////

void mHPSetBrTopWb(double value) {
    try {
        pred.setBrTopWb(value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetBrTopWb() {
    try {
        putReal(pred.brTopWb());
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetSMCxn(double mass, const char *collstring, const char *prodstring) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        auto prod = enumFromName<Higgs::predictions::Production>(prodstring);
        auto SMHiggs = Higgs::predictions::SMHiggs(mass);
        putReal(SMHiggs.cxn(coll, prod));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetSMBr(double mass, const char *decaystring) {
    try {
        auto decay = enumFromName<Higgs::predictions::Decay>(decaystring);
        auto SMHiggs = Higgs::predictions::SMHiggs(mass);
        putReal(SMHiggs.br(decay));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetSMTotalWidth(double mass) {
    try {
        auto SMHiggs = Higgs::predictions::SMHiggs(mass);
        putReal(SMHiggs.totalWidth());
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPAddParticle(const char *id, double mass, const char *echargestring,
                    const char *cpstring) {
    try {
        auto echarge = enumFromName<Higgs::predictions::ECharge>(echargestring);
        auto cp = enumFromName<Higgs::predictions::CP>(cpstring);
        auto &part =
            pred.addParticle(Higgs::predictions::BsmParticle{id, echarge, cp});
        part.setMass(mass);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPRemoveParticle(const char *id) {
    try {
        pred.removeParticle(id);
    } catch (const std::exception &e) {
        putError(e.what());
    }
    putSymbol("Null");
}

void mHPGetParticleIDs() {
    putFunction("List", pred.particleIds().size());
    for (const auto &id : pred.particleIds()) {
        putFunction("List", 2);
        putString(id.c_str());
        putReal(pred.particle(id).mass());
    }
}

void mHPGetCP(const char *id) {
    try {
        putString(magic_enum::enum_name(pred.particle(id).cp()).data());
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetCharge(const char *id) {
    try {
        putString(magic_enum::enum_name(pred.particle(id).charge()).data());
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetMass(const char *id, double value) {
    try {
        pred.particle(id).setMass(value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetMass(const char *id) {
    try {
        putReal(pred.particle(id).mass());
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetMassUnc(const char *id, double value) {
    try {
        pred.particle(id).setMassUnc(value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetMassUnc(const char *id) {
    try {
        putReal(pred.particle(id).massUnc());
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetCxn(const char *id, const char *collstring, const char *prodstring,
               double value) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        auto prod = enumFromName<Higgs::predictions::Production>(prodstring);
        pred.particle(id).setCxn(coll, prod, value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetNormalizedCxn(const char *id, const char *collstring,
                         const char *prodstring, double value) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        auto prod = enumFromName<Higgs::predictions::Production>(prodstring);
        pred.particle(id).setNormalizedCxn(
            coll, prod, value, Higgs::predictions::ReferenceModel::SMHiggs);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetCxn(const char *id, const char *collstring, const char *prodstring) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        auto prod = enumFromName<Higgs::predictions::Production>(prodstring);
        WSPutReal(stdlink, pred.particle(id).cxn(coll, prod));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetDecayWidth1(const char *id, const char *decaystring, double value) {
    try {
        auto decay = enumFromName<Higgs::predictions::Decay>(decaystring);
        pred.particle(id).setDecayWidth(decay, value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetDecayWidth2(const char *id, const char *decaystring1,
                       const char *decaystring2, double value) {
    try {
        auto chaindecay =
            magic_enum::enum_cast<Higgs::predictions::ChainDecay>(decaystring1);
        if (chaindecay) {
            pred.particle(id).setDecayWidth(chaindecay.value(), decaystring2,
                                            value);
        } else {
            pred.particle(id).setDecayWidth(decaystring1, decaystring2, value);
        }
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetBR1(const char *id, const char *decaystring, double value) {
    try {
        auto decay = enumFromName<Higgs::predictions::Decay>(decaystring);
        pred.particle(id).setBr(decay, value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetBR1(const char *id, const char *decaystring) {
    try {
        auto decay = enumFromName<Higgs::predictions::Decay>(decaystring);
        putReal(pred.particle(id).br(decay));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetBR2(const char *id, const char *decaystring1,
               const char *decaystring2, double value) {
    try {
        auto chaindecay =
            magic_enum::enum_cast<Higgs::predictions::ChainDecay>(decaystring1);
        if (chaindecay) {
            pred.particle(id).setBr(chaindecay.value(), decaystring2, value);
        } else {
            pred.particle(id).setBr(decaystring1, decaystring2, value);
        }
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetBR2(const char *id, const char *decaystring1,
               const char *decaystring2) {
    try {
        auto chaindecay =
            magic_enum::enum_cast<Higgs::predictions::ChainDecay>(decaystring1);
        if (chaindecay) {
            putReal(pred.particle(id).br(chaindecay.value(), decaystring2));
        } else {
            putReal(pred.particle(id).br(decaystring1, decaystring2));
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPEffectiveCouplingInput(
    const char *id, double uu_re, double uu_im, double dd_re, double dd_im,
    double cc_re, double cc_im, double ss_re, double ss_im, double tt_re,
    double tt_im, double bb_re, double bb_im, double ee_re, double ee_im,
    double mumu_re, double mumu_im, double tautau_re, double tautau_im,
    double WW, double ZZ, double Zgam, double gamgam, double gg,
    const char *refModel, const char *calcggH, const char *calcHgamgam) {
    try {
        auto effC = Higgs::predictions::NeutralEffectiveCouplings();
        effC.uu = std::complex<double>(uu_re, uu_im);
        effC.dd = std::complex<double>(dd_re, dd_im);
        effC.cc = std::complex<double>(cc_re, cc_im);
        effC.ss = std::complex<double>(ss_re, ss_im);
        effC.tt = std::complex<double>(tt_re, tt_im);
        effC.bb = std::complex<double>(bb_re, bb_im);
        effC.ee = std::complex<double>(ee_re, ee_im);
        effC.mumu = std::complex<double>(mumu_re, mumu_im);
        effC.tautau = std::complex<double>(tautau_re, tautau_im);
        effC.WW = WW;
        effC.ZZ = ZZ;
        effC.Zgam = Zgam;
        effC.gamgam = gamgam;
        effC.gg = gg;
        Higgs::predictions::effectiveCouplingInput(
            pred.particle(id), effC,
            enumFromName<Higgs::predictions::ReferenceModel>(refModel),
            StringToBool(calcggH), StringToBool(calcHgamgam));
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetCoupling(const char *id, const char *coupstring, double value) {
    try {
        auto coup = enumFromName<Higgs::predictions::Coupling>(coupstring);
        pred.particle(id).setCoupling(coup, value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetCoupling(const char *id, const char *coupstring) {
    try {
        auto coup = enumFromName<Higgs::predictions::Coupling>(coupstring);
        auto value = pred.particle(id).coupling(coup);
        if (value) {
            putReal(*value);
        } else {
            putSymbol("Null");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetTotalWidth(const char *id, double value) {
    try {
        pred.particle(id).setTotalWidth(value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetTotalWidth(const char *id) {
    try {
        putReal(pred.particle(id).totalWidth());
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetChannelRate(const char *id, const char *collstring,
                       const char *prodstring, const char *decaystring,
                       double value) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        auto prod = enumFromName<Higgs::predictions::Production>(prodstring);
        auto decay = enumFromName<Higgs::predictions::Decay>(decaystring);
        pred.particle(id).setChannelRate(coll, prod, decay, value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetChannelRate(const char *id, const char *collstring,
                       const char *prodstring, const char *decaystring) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        auto prod = enumFromName<Higgs::predictions::Production>(prodstring);
        auto decay = enumFromName<Higgs::predictions::Decay>(decaystring);
        putReal(pred.particle(id).channelRate(coll, prod, decay));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPResetChannelRates(const char *id) {
    try {
        pred.particle(id).resetChannelRates();
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPSetBsmPairCxn(const char *collstring, const char *id1, const char *id2,
                      double value) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        pred.setBsmPairCxn(coll, id1, id2, value);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHPGetBsmPairCxn(const char *collstring, const char *id1,
                      const char *id2) {
    try {
        auto coll = enumFromName<Higgs::predictions::Collider>(collstring);
        putReal(pred.bsmPairCxn(coll, id1, id2));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

////////////////////////////////////////////////////////////////////////////////
// HiggsBounds functions
////////////////////////////////////////////////////////////////////////////////

void mHBInitialize(const char *path, double res,
                   const char *clustermassuncstring,
                   const char *applicablemassuncstring,
                   const char *setlimitmassuncstring) {
    try {
        auto clustermassunc =
            enumFromName<Higgs::predictions::MassUncEagerness>(
                clustermassuncstring);
        auto applicablemassunc =
            enumFromName<Higgs::predictions::MassUncEagerness>(
                applicablemassuncstring);
        auto setlimitmassunc =
            enumFromName<Higgs::predictions::MassUncEagerness>(
                setlimitmassuncstring);

        HBoptions.applicableResolutionFac = res;
        HBoptions.clusterMassUnc = clustermassunc;
        HBoptions.applicableMassUnc = applicablemassunc;
        HBoptions.setLimitMassUnc = setlimitmassunc;

        bounds = std::make_unique<Higgs::Bounds>(path, HBoptions);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHBRetrieveOptions() {
    try {
        putFunction("List", 4);
        putRealRule("applicableResolutionFac",
                    HBoptions.applicableResolutionFac);
        putStringRule("clusterMassUnc",
                      magic_enum::enum_name(HBoptions.clusterMassUnc).data());
        putStringRule(
            "applicableMassUnc",
            magic_enum::enum_name(HBoptions.applicableMassUnc).data());
        putStringRule("setLimitMassUnc",
                      magic_enum::enum_name(HBoptions.setLimitMassUnc).data());
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHBListLimits() {
    try {
        if (bounds != nullptr) {
            putFunction("List", (*bounds).limits().size());
            for (const auto &lim : (*bounds).limits()) {
                putFunction("List", 9);

                putIntRule("id", (*lim).id());
                putStringRule("reference", (*lim).reference().c_str());
                putStringRule("process", (*lim).processDesc().c_str());
                putStringRule("limitExtent", (*lim).extentDesc().c_str());
                putStringRule("citeKey", (*lim).citeKey().c_str());
                putStringRule("collider",
                              magic_enum::enum_name((*lim).collider()).data());
                putStringRule(
                    "experiment",
                    magic_enum::enum_name((*lim).experiment()).data());
                putRealRule("luminosity", (*lim).luminosity());
                putStringRule("loadedfrom", (*lim).loadedFrom().c_str());
            }
        } else {
            putError("Need to call HBInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHBApplyBounds() {
    try {
        if (bounds != nullptr) {
            resultHB = (*bounds)(pred);
            putSymbol(BoolToString(resultHB.allowed));
        } else {
            putError("Need to call HBInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHBGetAppliedBounds() {
    try {
        putFunction("List", (resultHB.appliedLimits).size());
        for (const auto &al : resultHB.appliedLimits) {
            putFunction("List", 7);

            putIntRule("limitID", al.limit()->id());
            putStringRule("process", al.limit()->processDesc().c_str());
            putRealRule("obsRatio", al.obsRatio());
            putRealRule("expRatio", al.expRatio());
            putRealRule("obsLikelihood", al.obsLikelihood());
            putRealRule("expLikelihood", al.expLikelihood());

            auto contrParts = al.contributingParticles();
            putListRule("contributingParticles", size(contrParts));
            for (const auto &id : contrParts) {
                putString(id.c_str());
            }
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHBGetSelectedBounds() {
    try {
        putFunction("List", (resultHB.selectedLimits).size());
        for (const auto &[p, lim] : resultHB.selectedLimits) {
            putFunction("List", 4);

            putStringRule("particleID", p.c_str());
            putRealRule("obsRatio", lim.obsRatio());
            putRealRule("expRatio", lim.expRatio());
            putStringRule("limitInfo", lim.limit()->to_string().c_str());
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

////////////////////////////////////////////////////////////////////////////////
// HiggsSignals functions
////////////////////////////////////////////////////////////////////////////////

void mHSInitialize(const char *path, const char *theoryMassUncPDFString,
                   double massSensitiveAssignmentRange,
                   double unassignedMassMeasurementPenalty,
                   double unassignedCouplingMeasurementPenalty,
                   const char *rescaleToRefMassString,
                   const char *whichCorrelationsString,
                   const char *ignoreTheoryUncertainties) {
    try {
        auto theoryMassUncPDF =
            enumFromName<Higgs::signals::PDF>(theoryMassUncPDFString);
        auto rescaleToRefMass = enumFromName<Higgs::signals::RescaleToRefMass>(
            rescaleToRefMassString);
        auto whichCorrelations =
            enumFromName<Higgs::signals::Correlations>(whichCorrelationsString);

        HSoptions.theoryMassUncPDF = theoryMassUncPDF;
        HSoptions.massSensitiveAssignmentRange = massSensitiveAssignmentRange;
        HSoptions.unassignedMassMeasurementPenalty =
            unassignedMassMeasurementPenalty;
        HSoptions.unassignedCouplingMeasurementPenalty =
            unassignedCouplingMeasurementPenalty;
        HSoptions.rescaleToRefMass = rescaleToRefMass;
        HSoptions.whichCorrelations = whichCorrelations;
        HSoptions.ignoreTheoryUncertainties =
            StringToBool(ignoreTheoryUncertainties);

        signals = std::make_unique<Higgs::Signals>(path, HSoptions);
        putSymbol("Null");
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHSRetrieveOptions() {
    try {
        putFunction("List", 7);
        putStringRule("theoryMassUncPdf",
                      magic_enum::enum_name(HSoptions.theoryMassUncPDF).data());
        putRealRule("massSensitiveAssignmentRange",
                    HSoptions.massSensitiveAssignmentRange);
        putRealRule("unassignedMassMeasurementPenalty",
                    HSoptions.unassignedMassMeasurementPenalty);
        putRealRule("unassignedCouplingMeasurementPenalty",
                    HSoptions.unassignedCouplingMeasurementPenalty);
        putStringRule("rescaleToRefMass",
                      magic_enum::enum_name(HSoptions.rescaleToRefMass).data());
        putStringRule(
            "whichCorrelations",
            magic_enum::enum_name(HSoptions.whichCorrelations).data());
        putSymbolRule("ignoreTheoryUncertainties",
                      BoolToString(HSoptions.ignoreTheoryUncertainties));
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHSGetChisq() {
    try {
        if (signals != nullptr) {
            putReal((*signals)(pred));
        } else {
            putError("Need to call HSInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHSGetChisqMeasurement(int id) {
    try {
        if (signals != nullptr) {
            for (const auto m : (*signals).measurements()) {
                if (m.id() == id) {
                    putReal(m(pred));
                }
            }
        } else {
            putError("Need to call HSInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHSGetObservableCount() {
    try {
        if (signals != nullptr) {
            putReal((*signals).observableCount());
        } else {
            putError("Need to call HSInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHSListMeasurements() {
    try {
        if (signals != nullptr) {
            putFunction("List", (*signals).measurements().size());
            for (const auto m : (*signals).measurements()) {
                putFunction("List", 10);

                putIntRule("id", m.id());
                putIntRule("Nsubmeasurements", m.nSubMeasurements());
                putRealRule("refMass", m.referenceMass());
                putStringRule("refModel",
                              magic_enum::enum_name(m.referenceModel()).data());
                putStringRule("reference", m.reference().c_str());
                putStringRule("citeKey", m.citeKey().c_str());
                putStringRule("collider",
                              magic_enum::enum_name(m.collider()).data());
                putStringRule("experiment",
                              magic_enum::enum_name(m.experiment()).data());
                putRealRule("luminosity", m.luminosity());
                putStringRule("loadedfrom", m.loadedFrom().c_str());
            }
        } else {
            putError("Need to call HSInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

void mHSListSubMeasurements(int id) {
    try {
        if (signals != nullptr) {
            for (const auto m : (*signals).measurements()) {
                if (m.id() == id) {
                    putFunction("List", m.subMeasurements().size());
                    for (const auto subm : m.subMeasurements()) {
                        putString(subm.first.c_str());
                    }
                }
            }
        } else {
            putError("Need to call HSInitialize first");
        }
    } catch (const std::exception &e) {
        putError(e.what());
    }
}

////////////////////////////////////////////////////////////////////////////////
// WSTP internals
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) { return WSMain(argc, argv); }
