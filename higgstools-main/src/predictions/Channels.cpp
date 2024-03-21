#include "Higgs/predictions/Channels.hpp"
#include <array>
#include <magic_enum.hpp>

using magic_enum::enum_contains;
namespace Higgs::predictions {

namespace {

#define DECAY_MODE(m) m = static_cast<int>(Decay::m)

//! Neutral scalar decay modes.
//! This is a subset of the Decay modes.
enum class NeutralDecay {
    DECAY_MODE(none),
    DECAY_MODE(uu),
    DECAY_MODE(dd),
    DECAY_MODE(cc),
    DECAY_MODE(ss),
    DECAY_MODE(tt),
    DECAY_MODE(bb),
    DECAY_MODE(ee),
    DECAY_MODE(mumu),
    DECAY_MODE(tautau),
    DECAY_MODE(WW),
    DECAY_MODE(ZZ),
    DECAY_MODE(Zgam),
    DECAY_MODE(gamgam),
    DECAY_MODE(gg),
    DECAY_MODE(directInv),
    DECAY_MODE(inv),
    DECAY_MODE(emu),
    DECAY_MODE(etau),
    DECAY_MODE(mutau),
    DECAY_MODE(uc),
    DECAY_MODE(ds),
    DECAY_MODE(db),
    DECAY_MODE(sb),
};

//! Charged scalar decay modes.
//! This is a subset of the Decay modes.
enum class ChargedDecay {
    DECAY_MODE(none),
    DECAY_MODE(ud),
    DECAY_MODE(us),
    DECAY_MODE(ub),
    DECAY_MODE(cd),
    DECAY_MODE(cs),
    DECAY_MODE(cb),
    DECAY_MODE(tb),
    DECAY_MODE(enu),
    DECAY_MODE(munu),
    DECAY_MODE(taunu),
    DECAY_MODE(WZ),
    DECAY_MODE(Wgam),
};

//! Charged scalar decay modes.
//! This is a subset of the Decay modes.
enum class DoublyChargedDecay {
    DECAY_MODE(none),
    DECAY_MODE(WWsamesign),
    DECAY_MODE(eesamesign),
    DECAY_MODE(mumusamesign),
    DECAY_MODE(tautausamesign),
    DECAY_MODE(emusamesign),
    DECAY_MODE(etausamesign),
    DECAY_MODE(mutausamesign),
};

#undef DECAY_MODE
constexpr bool allDecaysAccountedFor() {
    bool check = true;
    for (auto x : magic_enum::enum_values<Decay>()) {
        check = check && (enum_contains(static_cast<NeutralDecay>(x)) ||
                          enum_contains(static_cast<ChargedDecay>(x)) ||
                          enum_contains(static_cast<DoublyChargedDecay>(x)));
    }
    return check;
}
} // namespace

bool validDecayFor(Decay d, ECharge charge) {
    static_assert(allDecaysAccountedFor(),
                  "All Decays must be either ChargedDecays, DoublyChargedDecay "
                  "or NeutralDecays. "
                  "If you just added a new decay mode please also add it to "
                  "the appropriate enum in Channels.cpp.");
    switch (charge) {
    case ECharge::neutral:
        return magic_enum::enum_contains(static_cast<NeutralDecay>(d));
    case ECharge::single:
        return magic_enum::enum_contains(static_cast<ChargedDecay>(d));
    case ECharge::doubly:
        return magic_enum::enum_contains(static_cast<DoublyChargedDecay>(d));
    }
    return false;
}

namespace {
#define PRODUCTION_MODE(m) m = static_cast<int>(Production::m)
//! Neutral scalar production modes at hadron colliders.
//! This is a subset of the Production modes.
enum class HadrNeutralProduction {
    PRODUCTION_MODE(none),
    PRODUCTION_MODE(pair),
    PRODUCTION_MODE(vbfH),
    PRODUCTION_MODE(HW),
    PRODUCTION_MODE(Htt),
    PRODUCTION_MODE(ggH),
    PRODUCTION_MODE(bbH),
    PRODUCTION_MODE(tchanHt),
    PRODUCTION_MODE(schanHt),
    PRODUCTION_MODE(qqHZ),
    PRODUCTION_MODE(ggHZ),
    PRODUCTION_MODE(bbHZ),
    PRODUCTION_MODE(HtW),
    PRODUCTION_MODE(qqH),
    PRODUCTION_MODE(brtHc),
    PRODUCTION_MODE(brtHu),

    PRODUCTION_MODE(uuHgam),
    PRODUCTION_MODE(ddHgam),
    PRODUCTION_MODE(ccHgam),
    PRODUCTION_MODE(ssHgam),
    PRODUCTION_MODE(bbHgam),
    PRODUCTION_MODE(ucHgam),
    PRODUCTION_MODE(dsHgam),
    PRODUCTION_MODE(dbHgam),
    PRODUCTION_MODE(sbHgam),

    PRODUCTION_MODE(H),
    PRODUCTION_MODE(HZ),
    PRODUCTION_MODE(Ht),
};

enum class LepNeutralProduction {
    PRODUCTION_MODE(none),
    PRODUCTION_MODE(pair),
    PRODUCTION_MODE(eeHZ),
    PRODUCTION_MODE(eeHbb),
    PRODUCTION_MODE(eeHtautau),
};

//! Charged scalar production modes at hadron colliders.
//! This is a subset of the Production modes.
enum class HadrChargedProduction {
    PRODUCTION_MODE(none),
    PRODUCTION_MODE(pair),
    PRODUCTION_MODE(qqHpm),
    PRODUCTION_MODE(vbfHpm),
    PRODUCTION_MODE(Hpmtb),
    PRODUCTION_MODE(HpmW),
    PRODUCTION_MODE(HpmZ),
    PRODUCTION_MODE(brtHpb),

    PRODUCTION_MODE(udHpgam),
    PRODUCTION_MODE(usHpgam),
    PRODUCTION_MODE(ubHpgam),
    PRODUCTION_MODE(cdHpgam),
    PRODUCTION_MODE(csHpgam),
    PRODUCTION_MODE(cbHpgam),
    PRODUCTION_MODE(udHmgam),
    PRODUCTION_MODE(usHmgam),
    PRODUCTION_MODE(ubHmgam),
    PRODUCTION_MODE(cdHmgam),
    PRODUCTION_MODE(csHmgam),
    PRODUCTION_MODE(cbHmgam),
};

//! Doubly charged scalar production modes at hadron colliders.
//! This is a subset of the Production modes.
enum class HadrDoublyChargedProduction {
    PRODUCTION_MODE(none),
    PRODUCTION_MODE(pair),
};

#undef PRODUCTION_MODE
constexpr bool allProdsAccountedFor() {
    bool check = true;
    for (auto x : magic_enum::enum_values<Production>()) {
        check = check &&
                (enum_contains(static_cast<HadrNeutralProduction>(x)) ||
                 enum_contains(static_cast<HadrChargedProduction>(x)) ||
                 enum_contains(static_cast<HadrDoublyChargedProduction>(x)) ||
                 enum_contains(static_cast<LepNeutralProduction>(x)));
    }
    return check;
}
} // namespace

bool validProductionFor(Production p, ECharge charge) {
    static_assert(
        allProdsAccountedFor(),
        "All production modes must be HadrNeutral, LepNeutral, "
        "HadrCharged, or LepCharged. If you just added a new decay mode, "
        "please also add it to the appropriate enum in Channels.cpp.");
    switch (charge) {
    case ECharge::neutral:
        return magic_enum::enum_contains(
                   static_cast<HadrNeutralProduction>(p)) ||
               magic_enum::enum_contains(static_cast<LepNeutralProduction>(p));
    case ECharge::single:
        return magic_enum::enum_contains(static_cast<HadrChargedProduction>(p));
    case ECharge::doubly:
        return magic_enum::enum_contains(
            static_cast<HadrDoublyChargedProduction>(p));
    }
    return false;
}

bool validProductionAt(Production p, ColliderType collType) {
    static_assert(magic_enum::enum_count<ColliderType>() == 2,
                  "You have changed the ColliderType and also need to update "
                  "the logic in Channels.cpp");
    switch (collType) {
    case ColliderType::pp:
        return enum_contains(static_cast<HadrNeutralProduction>(p)) ||
               enum_contains(static_cast<HadrChargedProduction>(p)) ||
               enum_contains(static_cast<HadrDoublyChargedProduction>(p));
    case ColliderType::ee:
        return enum_contains(static_cast<LepNeutralProduction>(p));
    }
    return false;
}

} // namespace Higgs::predictions
