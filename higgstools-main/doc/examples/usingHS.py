import Higgs.predictions as HP
import Higgs.signals as HS


signals = HS.Signals("path/to/HSDataSet")

pred = HP.Predictions()
# set model predictions on pred

hsChisq = signals(pred)
print(hsChisq)
