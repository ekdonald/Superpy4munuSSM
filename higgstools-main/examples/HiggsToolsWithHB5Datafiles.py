#!/usr/bin/env python3

from sys import argv
import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS
from Higgs.tools.Input import readHB5Datafiles, predictionsFromDict

if len(argv) not in [6, 8]:
    print(
        f"Usage: {argv[0]} datafile/prefix/path nHzero nHplus path/to/HBDataset path/to/HSDataset [calcggH] [calcHgamgam]"
    )
    if len(argv) == 7:
        print(
            f"       Default arguments calcggH and calcHgamgam have to be specified together."
        )
    exit(1)

bounds = HB.Bounds(argv[4])
signals = HS.Signals(argv[5])

nHzero = int(argv[2])
nHplus = int(argv[3])

try:
    if argv[6] == "True":
        calcggH = True
    elif argv[6] == "False":
        calcggH = False
    else:
        print("Invalid value for optional argument calcggH. Must be 'True' or 'False' (default).")
        exit(1)
    if argv[7] == "True":
        calcHgamgam = True
    elif argv[7] == "False":
        calcHgamgam = False
    else:
        print("Invalid value for optional argument calcHgamgam. Must be 'True' or 'False' (default).")
        exit(1)
except IndexError:
    calcggH = True
    calcHgamgam = True

neutralIds = [f"H{i+1}" for i in range(nHzero)]
chargedIds = [f"H{i+1}+" for i in range(nHplus)]

print(f"-> reading input from datafiles with prefix {argv[1]}")

df = readHB5Datafiles(argv[1], neutralIds, chargedIds)
# you can easily set additional input by adding columns to df (see the
# documentation of `predictionsFromDict`)


def runHBHS(datapoint):
    pred = predictionsFromDict(
        datapoint,
        neutralIds,
        chargedIds,
        [],
        calcggH=calcggH,
        calcHgamgam=calcHgamgam)
    # alternatively you can set additional input directly on pred in here
    return (bounds(pred), signals(pred))

print(f"-> running HiggsTools")

df[["HBResult", "HSChisq"]] = df.apply(runHBHS, axis=1, result_type="expand")

result = df[[f"m_{x}" for x in neutralIds + chargedIds]].copy()

result["HSchisq"] = df.HSChisq

result["HBresult"] = df.HBResult.apply(lambda x: int(x.allowed))
for h in neutralIds + chargedIds:

    def ifAnyLimitWasSelected(func):
        def inner(res):
            if h in res.selectedLimits.keys():
                return func(res.selectedLimits[h])
            return None

        return inner

    result[f"selLim_{h}"] = df.HBResult.map(
        ifAnyLimitWasSelected(lambda selLim: selLim.limit().reference())
    )
    result[f"selLim_{h}_obsRatio"] = df.HBResult.map(
        ifAnyLimitWasSelected(lambda selLim: selLim.obsRatio())
    )
    result[f"selLim_{h}_contributing"] = df.HBResult.map(
        ifAnyLimitWasSelected(lambda selLim: " ".join(selLim.contributingParticles()))
    )


outfile = f"{argv[1]}HiggsTools_results.tsv"

print(f"-> writing output to {outfile}")
result.to_csv(outfile, sep="\t")
