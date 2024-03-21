import Higgs.predictions as HP
import Higgs.bounds as HB


bounds = HB.Bounds("path/to/HBDataSet")

pred = HP.Predictions()
# set model predictions on pred

hbResult = bounds(pred)
print(hbResult)
