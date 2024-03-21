import Higgs.predictions as HP

h = HP.BsmParticle("h", charge="neutral", cp="even")

h.setTotalWidth(0.3)
h.setBr("bb", 0.8)
h.setBr("tautau", 0.17)
h.setBr("gamgam", 3e-2)
# BR(h->bb)=80%, BR(h->tautau)=17%, BR(h->gamgam)=3%

try:
    h.setBr("WW", 0.11)
except:
    print("this will be an error, since sum(BRs) is already 100%")
