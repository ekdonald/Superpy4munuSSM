import Higgs.predictions as HP
import Higgs.signals as HS

# setup inclusive model predictions with a single SM-like Higgs
pred = HP.Predictions()
h = pred.addParticle(HP.BsmParticle("h", charge="neutral", cp="even"))
h.setMass(125.09)
HP.effectiveCouplingInput(h, HP.smLikeEffCouplings)

signals = HS.Signals("path/to/HSDataset")

modFacs = {
    210306956: {
        "ggH_1J_low_pTH": [0.8],
        "ggH_1J_med_pTH": [0.9],
        "ggH_1J_high_pTH": [1.2],
    }
}

chisqMod = signals(pred, modificationFactors=modFacs)
