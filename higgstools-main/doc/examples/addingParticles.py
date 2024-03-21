import Higgs.predictions as HP

pred = HP.Predictions()
h1 = pred.addParticle(HP.BsmParticle("h1", "neutral", "even"))
hc = pred.addParticle(HP.BsmParticle("hc", "single"))
# set properties on h1 and hc
