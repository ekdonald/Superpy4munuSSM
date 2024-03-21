/**
 * @file Channels.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Channels and Processes modeled by HiggsPredictions
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "Higgs/HiggsTools_export.h"
#include "Higgs/predictions/Basics.hpp"
#include <limits>
namespace Higgs {
namespace predictions {

//! All implemented particle production modes.
enum class Production {
    //! no production mode. The associated Particle::cxn() for this is always
    //! zero. If used in Particle::channelRate() the production process is
    //! ignored and the Particle::br() of the channel decay mode is returned
    //! instead.
    none,
    //! Non-resonant pair production of this particle
    pair,
    /**<
     * @rststar
     * Production of a neutral particle :math:`H`
     * """"""""""""""""""""""""""""""""""""""""""
     * @endrststar
     */
    ggH,     //!< \f$ gg \to H \f$
    bbH,     //!< \f$ b\bar{b} \to H \f$ (5FS) or \f$ gg\to H b\bar{b} \f$ (4FS)
    vbfH,    //!< \f$ pp \to qq H \f$ in VBF
    HW,      //!< \f$ pp \to W^\pm H \f$
    qqHZ,    //!< \f$ q\bar{q} \to Z H \f$
    ggHZ,    //!< \f$ gg \to Z H \f$
    bbHZ,    //!< \f$ b\bar{b} \to Z H \f$
    Htt,     //!< \f$ pp \to t\bar{t} H \f$
    tchanHt, //!< \f$ pp \to H t \f$ in the \f$t\f$-channel
    schanHt, //!< \f$ pp \to H t \f$ in the \f$s\f$-channel
    HtW,     //!< \f$ pp \to H t W \f$
    qqH,     //!< \f$ q\bar{q} \to H \f$ where \f$q\in\{u,d,c,s\}\f$
    eeHZ,    //!< \f$ e^+e^- \to Z H \f$
    eeHbb,   //!< \f$ e^+ e^- \to b\bar{b} H \f$
    eeHtautau, //!< \f$ e^+ e^- \to \tau^+\tau^- H \f$
    //! @brief \f$ pp\to t\bar{t}, t\to H c\f$:
    //! This does not refer to the cxn but to \f$\mathrm{BR}(t\to H c)\f$.
    brtHc,
    //! @brief \f$ pp\to t\bar{t}, t\to H u\f$:
    //! This does not refer to the cxn but to \f$\mathrm{BR}(t\to H u)\f$.
    brtHu,
    //! @brief \f$ pp \to H \f$: `ggH + bbH + qqH`
    //! Combined production mode, can't be set directly.
    H,
    //! @brief \f$ pp \to HZ \f$: `qqHZ + ggHZ + bbHZ`
    //! Combined production mode, can't be set directly.
    HZ,
    //! @brief \f$ pp \to Ht \f$: `tchanHt + schanHt`
    //! Combined production mode, can't be set directly.
    Ht,
    uuHgam, //!< \f$ u\bar{u}\to H \gamma\f$
    ddHgam, //!< \f$ d\bar{d}\to H \gamma\f$
    ccHgam, //!< \f$ c\bar{c}\to H \gamma\f$
    ssHgam, //!< \f$ s\bar{s}\to H \gamma\f$
    bbHgam, //!< \f$ b\bar{b}\to H \gamma\f$
    ucHgam, //!< \f$ u\bar{c}+\bar{u}c\to H \gamma\f$
    dsHgam, //!< \f$ d\bar{s}+\bar{d}s\to H \gamma\f$
    dbHgam, //!< \f$ d\bar{b}+\bar{d}b\to H \gamma\f$
    sbHgam, //!< \f$ s\bar{b}+\bar{s}b\to H \gamma\f$

    /**<
     * @rststar
     * Production of a singly charged particle :math:`H^\pm`
     * """""""""""""""""""""""""""""""""""""""""""""""""""""
     * @endrststar
     */
    Hpmtb,  //!< \f$ gb\to H^\pm t \f$ (5FS) or \f$ gg \to t b H^\pm \f$ (4FS)
    qqHpm,  //!< \f$ q\bar{q'} \to H^\pm \f$ where \f$q,q'\in\{u,d,c,s\}\f$
    vbfHpm, //!< \f$ pp \to qq H^\pm \f$ in VBF
    HpmW,   //!< \f$ pp \to H^\pm W^\mp \f$
    HpmZ,   //!< \f$ pp \to H^\pm Z \f$
    //! @brief \f$ pp\to t\bar{t}, t\to H^+ b\f$:
    //! This does not refer to the cxn but to \f$\mathrm{BR}(t\to H^+ b)\f$.
    brtHpb,
    udHpgam, //!< \f$ u\bar{d} \to H^+ \gamma\f$
    usHpgam, //!< \f$ u\bar{s} \to H^+ \gamma\f$
    ubHpgam, //!< \f$ u\bar{b} \to H^+ \gamma\f$
    cdHpgam, //!< \f$ c\bar{d} \to H^+ \gamma\f$
    csHpgam, //!< \f$ c\bar{s} \to H^+ \gamma\f$
    cbHpgam, //!< \f$ c\bar{b} \to H^+ \gamma\f$
    udHmgam, //!< \f$ \bar{u}d \to H^- \gamma\f$
    usHmgam, //!< \f$ \bar{u}s \to H^- \gamma\f$
    ubHmgam, //!< \f$ \bar{u}b \to H^- \gamma\f$
    cdHmgam, //!< \f$ \bar{c}d \to H^- \gamma\f$
    csHmgam, //!< \f$ \bar{c}s \to H^- \gamma\f$
    cbHmgam, //!< \f$ \bar{c}b \to H^- \gamma\f$
};

//! All implemented decay modes.
enum class Decay {
    //! no decay mode. The associated Particle.br() for this is always
    //! zero, but the value is useful in Channel%s, where the
    //! Particle.channelRate() of a Channel with Channel.decay()==Decay::none
    //! is simply the Particle.cxn() of Channel.prod().
    none,
    /**<
     * @rststar
     * Elementary neutral final states
     * """""""""""""""""""""""""""""""
     * @endrststar
     */
    uu,        //!< \f$ u\bar{u} \f$
    dd,        //!< \f$ d\bar{d} \f$
    cc,        //!< \f$ c\bar{c} \f$
    ss,        //!< \f$ s\bar{s} \f$
    tt,        //!< \f$ t\bar{t} \f$
    bb,        //!< \f$ b\bar{b} \f$
    ee,        //!< \f$ e^+e^- \f$
    mumu,      //!< \f$ \mu^+\mu^- \f$
    tautau,    //!< \f$ \tau^+\tau^- \f$
    WW,        //!< \f$ W^+W^- \f$
    ZZ,        //!< \f$ ZZ \f$
    Zgam,      //!< \f$ Z \gamma \f$
    gamgam,    //!< \f$ \gamma\gamma \f$
    gg,        //!< \f$ gg \f$
    directInv, //!< direct decays into invisible
    emu,       //!< \f$ e^\pm \mu^\mp \f$
    etau,      //!< \f$ e^\pm \tau^\mp \f$
    mutau,     //!< \f$ \mu^\pm \tau^\mp \f$
    uc, //!< \f$ u\bar{c} + \bar{u}c \f$
    ds, //!< \f$ d\bar{s} + \bar{d}s \f$
    db, //!< \f$ d\bar{b} + \bar{d}b \f$
    sb, //!< \f$ s\bar{b} + \bar{s}b \f$
    /**<
     * @rststar
     * Combined neutral final states
     * """""""""""""""""""""""""""""
     * @endrststar
     */
    inv, //!< full decays into invisible. Includes Decay::directInv as well as
         //!< \f$ ZZ\to 4\nu \f$, \f$ Z X \to \nu\nu + \mathrm{inv} \f$, and
         //!< \f$ X1 X2 \to \mathrm{inv}\f$, where they invisible decays of the BSM
         //!< particles \f$X\f$ are themselves the full invisible BR.

    /**<
     * @rststar
     * Singly charged final states
     * """""""""""""""""""""""""""
     * @endrststar
     */
    ud,    //!< \f$u\bar{d} \f$
    us,    //!< \f$u\bar{s} \f$
    ub,    //!< \f$u\bar{b} \f$
    cd,    //!< \f$c\bar{d} \f$
    cs,    //!< \f$c\bar{s} \f$
    cb,    //!< \f$c\bar{b} \f$
    tb,    //!< \f$t\bar{b} \f$
    enu,   //!< \f$e^+\nu_e \f$
    munu,  //!< \f$\mu^+\nu_\mu \f$
    taunu, //!< \f$\tau^+\nu_\tau \f$
    WZ,    //!< \f$W^+Z \f$
    Wgam,  //!< \f$W^+ \gamma \f$

    /**<
     * @rststar
     * Doubly charged final states
     * """""""""""""""""""""""""""
     * @endrststar
     */
    WWsamesign,     //!< \f$ W^\pm W^\pm \f$
    eesamesign,     //!< \f$ e^\pm e^\pm \f$
    mumusamesign,   //!< \f$ \mu^\pm \mu^\pm \f$
    tautausamesign, //!< \f$ \tau^\pm \tau^\pm \f$
    emusamesign,    //!< \f$ e^\pm \mu^\pm \f$
    etausamesign,   //!< \f$ e^\pm \tau^\pm \f$
    mutausamesign,  //!< \f$ \mu^\pm \tau^\pm \f$

};

//! Chain decays with final states containing one BSM particle
enum class ChainDecay {
    Z, //!< \f$ \phi_i \to \phi_j Z \f$
    W, //!< \f$ \phi_i \to \phi_j W^\pm \f$, where the charged difference
       //!< between \f$\phi_{i,j}\f$ is one
};

//! Couplings that some searches or measurements are sensitive to
enum class Coupling {
    effCPeTopYuk, //!< tree-level SM-normalized CP-even top Yukawa coupling
    effCPoTopYuk, //!< tree-level SM-normalized CP-odd top Yukawa coupling
    effVV,        //!< tree-level SM-normalized WW/ZZ coupling
    //! the CP-phase of the tree-level SM-normalized tau Yukawa coupling as
    //! defined in eq. (2) of 2110.04836. The overall sign of
    //! \f$\kappa + i\tilde \kappa\f$ is ignored such that the value lies
    //! in [-pi/2, pi/2].
    alphaCPTauYuk,
};

//! Check if the production mode p is valid for a particle of the given
//! charge. Production::none is always valid.
HIGGSTOOLS_EXPORT bool validProductionFor(Production p, ECharge charge);

//! Check if the decay mode d is valid for a particle of the given charge.
//! Decay::none is always valid.
HIGGSTOOLS_EXPORT bool validDecayFor(Decay d, ECharge charge);

//! Check if the production mode p is valid at the given type of collider.
//! Production::none is always valid.
HIGGSTOOLS_EXPORT bool validProductionAt(Production prod,
                                         ColliderType collType);

//! specified relative and absolute mass resolution
struct HIGGSTOOLS_EXPORT MassResolution {
    //! relative mass resolution
    double relative = 0.;
    //! absolute mass resolution, tiny non-zero default value to account for
    //! numerical inaccuracies
    double absolute = 1e3 * std::numeric_limits<double>::min();
};

} // namespace predictions
} // namespace Higgs
