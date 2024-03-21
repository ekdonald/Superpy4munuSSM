#########################################################################
#                                                                       #
#    L i k e l i h o o d                                                #
#                                                                       #
#########################################################################

# External modules.
import os
import math
import scipy
import subprocess
import pyslha213mod
import tempfile
import pandas as pd
import numpy as NP
import sys
# import StringIO
from io import StringIO

import shutil			# ADDED

# SuperPy modules.
import Priors
import Cube
#import fastlim_superpy


#################################################################################################

# Read the InputFile to see which subprograms and constraints will be taken into account

def ReadInputFile(filename, start):								#ADDED
        """ Read from InputFile if MultiNest has to take into account certain constraints 

        Arguments:
        filename -- File to be read from.
        start -- String of beginning of line to read.


        Returns:
        Parameter retrieved from the file.

        """

        aux = 0

        for line in open(filename, 'r'):
            if line.lstrip().startswith(start):
                aux = line.split(start[-1])[-1]
                parameter = aux[:-1]
                if parameter not in ("True","False"):
                    print ("ERROR: Check InputFile, write True or False in:", start)
                    sys.exit()
                return parameter

        if aux == 0:
            print ("ERROR: Check InputFile, there is no:", start)
            sys.exit()



# NEUTRALINOS LSP?
          
NeutralinosLSPOK = ReadInputFile('InputFile','NeutralinosLSPOK = ')


# SUBPROGRAMS

##MicrOmegas          
MicrOmegasOK = ReadInputFile('InputFile','MicrOmegasOK = ')
constraintRelicDensity = ReadInputFile('InputFile','constraintRelicDensity = ')

#constraintRelicDensityUpperLimit = ReadInputFile('InputFile','constraintRelicDensityUpperLimit = ')
addcubeoh2 = ReadInputFile('InputFile','addcubeoh2 = ')
addSISDpn = ReadInputFile('InputFile','addSISDpn = ')


## HiggsSignals & HiggsBounds
HiggsSignalsOK = ReadInputFile('InputFile','HiggsSignalsOK = ')
constraintHiggsBS = ReadInputFile('InputFile','constraintHiggsBS = ')


# HiggsSignals & HiggsBounds
HiggsToolsOK = ReadInputFile('InputFile','HiggsToolsOK = ')
constraintHiggsBS = ReadInputFile('InputFile','constraintHiggsBS = ')


# SModelS
SModelSOK = ReadInputFile('InputFile','SModelSOK = ')
#constraintSModelS = ReadInputFile('InputFile','constraintSModelS = ')


# DDCalc
DDCalcOK = ReadInputFile('InputFile','DDCalcOK = ')

DDCalcXENON1T_2017 = ReadInputFile('InputFile','DDCalcXENON1T_2017 = ')
DDCalcCRESST_2017 = ReadInputFile('InputFile','DDCalcCRESST_2017 = ')
DDCalcDARWIN = ReadInputFile('InputFile','DDCalcDARWIN = ')
DDCalcPICO_60_2017 = ReadInputFile('InputFile','DDCalcPICO_60_2017 = ')
DDCalcXENON1T_2018 = ReadInputFile('InputFile','DDCalcXENON1T_2018 = ')

DoDDCalcXENON1T_2017 = ReadInputFile('InputFile','DoDDCalcXENON1T_2017 = ')
DoDDCalcCRESST_2017 = ReadInputFile('InputFile','DoDDCalcCRESST_2017 = ')
DoDDCalcDARWIN = ReadInputFile('InputFile','DoDDCalcDARWIN = ')
DoDDCalcPICO_60_2017 = ReadInputFile('InputFile','DoDDCalcPICO_60_2017 = ')
DoDDCalcXENON1T_2018 = ReadInputFile('InputFile','DoDDCalcXENON1T_2018 = ')

StorepvalueDDCalc = ReadInputFile('InputFile','StorepvalueDDCalc = ')


# Neutrino Physics, i.e. mixing angles and masses
NeutrinoPhysicsOK = ReadInputFile('InputFile','NeutrinoPhysicsOK = ')
constraintNeutrinoPhysics = ReadInputFile('InputFile','constraintNeutrinoPhysics = ')
constraintSumNeutrinos = ReadInputFile('InputFile','constraintSumNeutrinos = ')


# Several observables calculated byb SPheno
SPhenoObservablesOK = ReadInputFile('InputFile','SPhenoObservablesOK = ')
constraintSPhenoObservables = ReadInputFile('InputFile','constraintSPhenoObservables = ')


# AUXILIARY FILES

DeleteOK = ReadInputFile('InputFile','DeleteOK = ')

DeleteSLHAinput = ReadInputFile('InputFile','DeleteSLHAinput = ')
DeleteSPheno = ReadInputFile('InputFile','DeleteSPheno = ')
DeleteMicrOmegas = ReadInputFile('InputFile','DeleteMicrOmegas = ')
DeleteHiggsSignals = ReadInputFile('InputFile','DeleteHiggsSignals = ')
DeleteHiggsTools = ReadInputFile('InputFile','DeleteHiggsTools = ')
DeleteSModelS = ReadInputFile('InputFile','DeleteSModelS = ')
DeleteDDCalc = ReadInputFile('InputFile','DeleteDDCalc = ')


# INPUT PARAMETER OPTIONS

InputParamFormat = Priors.ReadInputFileParam('InputFileParam','InputParamFormat = ')   # Priors

if InputParamFormat == 'D':
    Ydaux11 = Priors.ReadInputFileParam('InputFileParam','Ydaux11 = ')
    Ydaux22 = Priors.ReadInputFileParam('InputFileParam','Ydaux22 = ')
    Ydaux33 = Priors.ReadInputFileParam('InputFileParam','Ydaux33 = ')

    Yuaux11 = Priors.ReadInputFileParam('InputFileParam','Yuaux11 = ')
    Yuaux22 = Priors.ReadInputFileParam('InputFileParam','Yuaux22 = ')
    Yuaux33 = Priors.ReadInputFileParam('InputFileParam','Yuaux33 = ')

    Yeaux11 = Priors.ReadInputFileParam('InputFileParam','Yeaux11 = ')
    Yeaux22 = Priors.ReadInputFileParam('InputFileParam','Yeaux22 = ')
    Yeaux33 = Priors.ReadInputFileParam('InputFileParam','Yeaux33 = ')


# SHOW IN TERMINAL OPTIONS
ShowPredictions = ReadInputFile('InputFile','ShowPredictions = ')
ShowTotalLoglike = ReadInputFile('InputFile','ShowTotalLoglike = ')
ShowSPhenoProcess = ReadInputFile('InputFile','ShowSPhenoProcess = ')


#################################################################################################

def myloglike(cube, ndim, nparams):
    """ MultiNest callback function.
    Calculate the model's predictions (by calling external programs),
    saves the extra parameters in the cube, and returns log likelihood to Multinest.

    Arguments:
    cube -- Unit hypercube from which model parameters are mapped.
    ndim -- Number of model paramers.
    nparams -- Total number of parameters in the cube.

    Returns: The total log likelihood for the model parameter
    point under consideration.

    """
    # Set up constraints class.
    Constraints = MUNUSSMConstraintTracker()							#MODIFIED

    # Copy cube to constraints, so it can work out predictions etc.
    for i, name in enumerate(sorted(Priors.MUNUSSMModelTracker().param.keys(), key=str.lower)):		#MODIFIED
        Constraints.param[name] = cube[i]

    # Set predictions and loglikes.
    Constraints.SetPredictions()
    Constraints.SetLogLike()

    # Copy constraints to cube.
    for name in sorted(Constraints.constraint.keys(), key=str.lower):
        Cube.AddCube(
            cube,
            Constraints.constraint[name].theory,
            Cube.label,
            name)

    # Copy associated chi2s to cube. Better to print chi2 than loglike,
    # beceause MultiNest prints chi2 and they can be treated in the
    # same way when plotting.
    for name in sorted(Constraints.constraint.keys(), key=str.lower):
        Cube.AddCube(cube, -
                     2 *
                     Constraints.constraint[name].loglike, Cube.label, 'chi2:' +
                     name)

    # Copy SLHA masses to the cube.
    for key in Constraints.masses:
        Cube.AddCube(
            cube,
            Constraints.masses[key],
            Cube.label,
            'Mass:' +
            str(key))

    # Print-out cube for debugging.
    if ShowPredictions == 'True':                           # ADDED option in InputFile
        print ('Predictions:')
        for label, param in zip(Cube.label.items(), cube):
            print (label, param)
    if ShowTotalLoglike == 'True':                          # ADDED option in InputFile
        print ('Total loglike', Constraints.loglike)

    # Start cube count from 0 again.
    Cube.AddCube.reset()

    # Print an info file for the cube.
    # Decorator insures this is printed only once.
    # Point must be physical, else the info will be incomplete.
    if Constraints.physical:
        Cube.PrintInfo()

    # Return the log likelihood to MultiNest.
    return Constraints.loglike


#########################################################################

# This class evalulates likelihoods and the model's predictions etc.

class MUNUSSMConstraintTracker:							#MODIFIED

    """ Contains SPhenomunuSSM3G constraints and functions to
    evaluate the model's predictions. """
    
    def __init__(self):
        """ Sets-up the model. """
        # Dictionary of model's parameters.
        self.param = {}
        # Dictionary of constraint information.
        self.constraint = {}
        # Whether the model has an acceptable mass spectrum, EWSB etc.
        self.physical = True
        # The total loglike associated with model.
        self.loglike = -1e101
        # List SLHA masses.
        self.masses = []
        #self.HTfilesroot = ''

        # You set the constraint values here. You must make sure that
        # you have computed all the relevant theoretical values etc.

        # Naturalness priors.
        # Implemented in modified SOFTSUSY.
        # Although it is a prior, easier to put here in the likelihood.
        #self.constraint['Natural'] = ExternalConstraint()				#MODIFIED

        # LEP, Higgs and Tevatron data.
        # Implemented in HiggsSignals.
        # http://higgsbounds.hepforge.org/
        if constraintHiggsBS == 'True':
            self.constraint['Higgs'] = ExternalConstraint()				#MODIFIED
            self.constraint['Higgspvalue'] = ExternalConstraint()				#MODIFIED
            self.constraint['HiggsBounds'] = ExternalConstraint()				#MODIFIED

        # LHC direct searches.
        # Implemented in Fast-Lim.
        # http://fastlim.web.cern.ch/fastlim/
        #self.constraint['LHC'] = ExternalConstraint()					#MODIFIED

        # Interpolate lower bound on (m0,m12) plane.
        # https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CombinedSummaryPlots/SUSY/ATLAS_SUSY_MSUGRA/ATLAS_SUSY_MSUGRA.png
        # ATLAS-CONF-2013-047
        #self.constraint['LHC_interp'] = InterpolateLowerConstraint(			#MODIFIED
        #    'atlas_m0m12.dat')

        # Relic density of neutralinos.
        # Planck.
        # http://arxiv.org/pdf/1303.5076v2.pdf
        #self.constraint['oh2'] = GaussConstraintFractionalTau(				#MODIFIED
        #    0.1199,
        #    0.0027,
        #    0.1)
        # https://arxiv.org/abs/1502.01589 Table 4 TT,TE,EE + lowP + lensing + (BAO + JLA + H0) [68 % limits] There are less constraining values in the paper
        if constraintRelicDensity == 'True':
            if constraintRelicDensityUpperLimit == 'True':
                self.constraint['oh2'] = UpperConstraint(0.1188)
            else:
                self.constraint['oh2'] = GaussConstraintFractionalTau(0.1188,0.11,0.1)#(0.1188,0.0010,0.1) # CAREFUL WITH THIS!

        if constraintRelicDensity == 'False' and addcubeoh2 == 'True':
            self.constraint['oh2'] = ExternalConstraint()

        if addSISDpn == 'True':
            self.constraint['SIproton'] = ExternalConstraint()
            self.constraint['SIneutron'] = ExternalConstraint()
            self.constraint['SDproton'] = ExternalConstraint()
            self.constraint['SDneutron'] = ExternalConstraint()
            
        # Direct Detection constraints.
        # Implemented in DDCalc.
        # http://ddcalc.hepforge.org/
        if DDCalcXENON1T_2017 == 'True':
            self.constraint['DDCalcXENON1T_2017'] = ExternalConstraint()				#MODIFIED
            if StorepvalueDDCalc == 'True':
                self.constraint['pvalueDDCalcXENON1T_2017'] = ExternalConstraint()
        if DDCalcCRESST_2017 == 'True':
            self.constraint['DDCalcCRESST_2017'] = ExternalConstraint()				#MODIFIED
            if StorepvalueDDCalc == 'True':
                self.constraint['pvalueDDCalcCRESST_2017'] = ExternalConstraint()
        if DDCalcDARWIN == 'True':
            self.constraint['DDCalcDARWIN'] = ExternalConstraint()				#MODIFIED
            if StorepvalueDDCalc == 'True':
                self.constraint['pvalueDDCalcDARWIN'] = ExternalConstraint()
        if DDCalcPICO_60_2017 == 'True':
            self.constraint['DDCalcPICO_60_2017'] = ExternalConstraint()				#MODIFIED
            if StorepvalueDDCalc == 'True':
                self.constraint['pvalueDDCalcPICO_60_2017'] = ExternalConstraint()
        if DDCalcXENON1T_2018 == 'True':
            self.constraint['DDCalcXENON1T_2018'] = ExternalConstraint()				#MODIFIED
            if StorepvalueDDCalc == 'True':
                self.constraint['pvalueDDCalcXENON1T_2018'] = ExternalConstraint()


        # Neutrino Physics. Mass differences and mixing angles
        # https://arxiv.org/abs/1703.04471
        if constraintNeutrinoPhysics == 'True':
            self.constraint['dm2'] = GaussConstraint(7.37e-5, 0.165e-5, 0.)		        	#MODIFIED
            self.constraint['Dm2'] = GaussConstraint(2.525e-3, 0.036e-3, 0.)		    	#MODIFIED

            self.constraint['sin2theta12'] = GaussConstraint(2.97e-1, 0.165e-1, 0.)			#MODIFIED
            self.constraint['sin2theta13'] = GaussConstraint(2.15e-2, 0.07e-2, 0.)			#MODIFIED
            self.constraint['sin2theta23'] = GaussConstraint(4.25e-1, 0.18e-1, 0.)			#MODIFIED

        # Neutrino Physics. Sum of neutrino masses
        # https://arxiv.org/pdf/1502.01589.pdf Planck 95% C.L. pag. 41
        if constraintSumNeutrinos == 'True':
            self.constraint['sumneutrinos'] = UpperConstraint(0.23)			#MODIFIED


        # Selected SPheno Observables, taken from the SPheno Output (some old observables calculated with SuperISO)
        if constraintSPhenoObservables == 'True':
        # Anomalous magnetic moment of muon.
        # PDG.
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-g-2-muon-anom-mag-moment.pdf      
        # http://pdg.lbl.gov/2016/reviews/rpp2016-rev-g-2-muon-anom-mag-moment.pdf UPDATED (same values) eq.(15):a_mu^exp-a_mu^SM
            self.constraint['gm2'] = GaussConstraint(28.8e-10, 8.0e-10, 1e-10)

        # BR(b -> s gamma).
        # PDG.
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-b-meson-prod-decay.pdf
        # http://pdg.lbl.gov/2016/reviews/rpp2016-rev-b-meson-prod-decay.pdf            UPDATED (same values)
        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf REVIEW OF PARTICLE PHYSICS 2016    UPDATED pag 1144
            self.constraint['bsg'] = GaussConstraint(3.43e-4, 0.22e-4, 0.21e-4)

        # BR(Bd -> mu mu).
        # PDG.
        # http://pdg.lbl.gov/2016/reviews/rpp2016-rev-b-meson-prod-decay.pdf
        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf REVIEW OF PARTICLE PHYSICS 2016    ADDED pag 1145
            self.constraint['bdmumu'] = GaussConstraint(3.9e-10, 1.6e-10, 0.)

        # BR(Bs -> mu mu).
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/BranchingRatio.action?parCode=S086&desig=15
        # http://pdg.lbl.gov/2016/reviews/rpp2016-rev-b-meson-prod-decay.pdf            UPDATED
        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf REVIEW OF PARTICLE PHYSICS 2016    UPDATED pag 1145
        #    self.constraint['bsmumu'] = GaussConstraintFractionalTau(3.1e-9, 0.7e-9, 0.14)
            self.constraint['bsmumu'] = GaussConstraint(2.8e-9, 0.7e-9, 0.)

        # BR(b -> tau nu).
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/BranchingRatio.action?parCode=S041&desig=18
        # https://arxiv.org/pdf/1611.00231.pdf    pag 9, table 2 Belle y Babar,   2016
        # https://arxiv.org/pdf/1503.05613.pdf      pag 3, current world average 2015
        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf REVIEW OF PARTICLE PHYSICS 2016    UPDATED pag 1144
            self.constraint['btaunu'] = GaussConstraint(1.14e-4, 0.27e-4, 0.38e-4)	

        # BR(mu->e gamma)
        # PDG.
        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf REVIEW OF PARTICLE PHYSICS 2016 pag 718
        # http://pdg.lbl.gov/2017/tables/rpp2017-sum-leptons.pdf             UPDATED
            self.constraint['mueg'] = UpperConstraint(4.2e-13)			

        # BR(mu->3e)
        # PDG.
        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf REVIEW OF PARTICLE PHYSICS 2016 pag 718
        # http://pdg.lbl.gov/2017/tables/rpp2017-sum-leptons.pdf             UPDATED
            self.constraint['mu3e'] = UpperConstraint(1.0e-12)			

        # W-boson mass.
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/DataBlock.action?node=S043
        #self.constraint['mw'] = GaussConstraint(80.385, 0.015, 0.015)			#MODIFIED

        # Leptonic sin eff theta, effective weak-mixing angle.
        # PDG.
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-standard-model.pdf
        #self.constraint['sineff'] = GaussConstraint(0.23126, 0.00005, 15e-5)		#MODIFIED

        # delta MBs.
        # PDG.
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-b-bar-mixing.pdf
        #self.constraint['deltaMb'] = GaussConstraint(17.761, 0.022, 2.4)		#MODIFIED

    def SetPredictions(self):
        """ Run the auxilliary programs for a particular model
        point to find the model's preddictions.
        Arguments:

        Returns:

        """

        #print ("... len(MUNUSSMConstraintTracker().constraint) = ", len(MUNUSSMConstraintTracker().constraint))

        print ("Subprograms:")
        
        
        # Call SPheno and read the mass spectrum.
        print ("Calling SPheno...")

        self.HTrootname = os.path.basename(tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=False).name)
        print ("self.HTrootname = ", self.HTrootname)

        self.munussm(self.HTrootname)
        print('readslha self.SLHA', os.path.basename(self.SLHA))
        self.readslha(self.SLHA)


        # Call auxillary programs if physical.
        if self.physical:
            if MicrOmegasOK == 'True':
                print ("Calling micrOMEGAs...")
                self.micromegas(self.SLHA)

        if self.physical:
            if HiggsSignalsOK == 'True':
                print ("Calling HiggsSignals...")
                self.higgssignals(self.SLHA)

        if self.physical:
            if HiggsToolsOK == 'True':
                print ("Calling higgstools...")
                self.higgstools(self.SLHA, self.HTrootname)

        if self.physical:
            if DDCalcOK == 'True':
                print ("Calling DDCalc...")
                self.ddcalc(self.microfile)

        if self.physical:
            if SModelSOK == 'True':
                print ("Calling SModelS...")
                self.smodels('inputFiles/slha/modSPhenoOutput.slha')  

        if self.physical:
            if NeutrinoPhysicsOK == 'True':
                print ("Calling Neutrino Physics routine...")
                self.neutrinophysics()             

        if self.physical:
            if SPhenoObservablesOK == 'True':
                print ("Calling Selected SPheno Obsevables...")
                self.sphenoobservables()  

        if DeleteOK == 'True':
            print ("Deleting selected auxiliary files generated by the subprograms...")
            self.delete()             
            print ('NOTE: Some empty files may remain in tmp folder')

        if DeleteOK == 'True':
            # remove ht input files from hbout3g directoory 
            try:
                del_file_patterns('../pyMUNUSSM/temFiles/', self.HTrootname)
            except OSError:
                #pass
                print ("SPHENO files inside not deleted")


    def munussm(self, HTrootname):									#ADDED    RUNS SPHENO
        """Call SPheno to calculate spectrum and decay tables for the MUNUSSM.
        Arguments:

        Returns:

        """

        # Create SLHA input file.
        self.SLHAIN = self.writeslha(self.param)
    

        print('test position 1 self.SLHAIN = ', self.SLHAIN)

        self.SLHA = RunProgramSPheno(
	    './bin/SPhenomunuSSM3G',
            '../SPhenomunuSSM',
            [self.SLHAIN, HTrootname])
        print('test position 2')

        self.CheckProgram(self.SLHA, ["ERROR"])
        
        #print('test position 3 self.SLHA', self.SLHA)



    def writeslha(self, param, MZ=9.11876000e+01):
        """ Write an SLHA input file, SLHAIN, for a given
        parameter point.

        Arguments:
        param -- Dictionary of model parameters.
        MZ -- Value of Z-boson mass.

        Return:
        Name of SLHA input file.

        """
        SLHAIN = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', mode = "w", delete=False)
        # Pass information in SLHA format. NB that you could easily add more
        # parameters to the EXTPAR section, to relax the CNMSSM.
        SLHAIN.write(
            """
Block MODSEL    #  
1 0             #  1/0: High/low scale input 
2 1             # Boundary Condition  
6 1             # Generation Mixing 
12 %s           # Renormalization scale 
Block SMINPUTS    # Standard Model inputs 
1 %s  # alpha_em^-1(MZ)^MSbar
2 1.166370E-05    # G_F,Fermi constant 
3 %s    # alpha_s(MZ) SM MSbar 
4 %s    # Z-boson pole mass 
5 %s    # m_b(mb) SM MSbar 
6 %s    # m_top(pole) 
7 1.776690E+00    # m_tau(pole) 
Block MINPAR      # Input parameters 
3   %s    # TanBeta
Block EXTPAR     # Input parameters 
65  %s           # vR1Input
66  %s           # vR2Input
67  %s           # vR3Input
200 %s           # vL1Input
201 %s           # vL2Input
202 %s           # vL3Input
Block SPhenoInput   # SPheno specific input 
  1 -1              # error level 
  2  0              # SPA conventions 
  6  1              #
  7  0              # Skip 2-loop Higgs corrections 
  8  3              # Method used for two-loop calculation 
  9  1              # Gaugeless limit used at two-loop 
 10  0              # safe-mode used at two-loop 
 11 1               # calculate branching ratios 
 13 0               # 3-Body decays: none (0), fermion (1), scalar (2), both (3) 
 14 0               # Run couplings to scale of decaying particle 
 12 1.000E-30       # write only branching ratios larger than this value 
 15 1.000E-30       # write only decay if width larger than this value 
 31 -1              # fixed GUT scale (-1: dynamical GUT scale) 
 32 0               # Strict unification 
 34 1.000E-03       # Precision of mass calculation 
 35 30               # Maximal number of iterations
 36 2               # Minimal number of iterations
 37 1               # Set Yukawa scheme  
 38 2               # 1- or 2-Loop RGEs 
 50 1               # Majorana phases: use only positive masses 
 51 0               # Write Output in CKM basis 
 52 0               # Write spectrum in case of tachyonic states 
 55 1               # Calculate one loop masses 
 57 1               # Calculate low energy constraints 
 65 1               # Solution tadpole equation 
 70 0               #
 75 0               # Write WHIZARD files 
 76 0               # Write HiggsBounds file 
 86 0               # Maximal width to be counted as invisible in Higgs decays; -1: only LSP 
510 0               # Write tree level values for tadpole solutions 
515 0               # Write parameter values at GUT scale 
520 0               # Write effective Higgs couplings (HiggsBounds blocks) 
525 0               # Write loop contributions to diphoton decay of Higgs 
530 0               # Write Blocks for Vevacious 
Block MSD2IN    #  
1 1   %s         # md2(1,1)
1 2   %s         # md2(1,2)
1 3   %s         # md2(1,3)
2 1   %s         # md2(2,1)
2 2   %s         # md2(2,2)
2 3   %s         # md2(2,3)
3 1   %s         # md2(3,1)
3 2   %s         # md2(3,2)
3 3   %s         # md2(3,3)
Block MSE2IN    #  
1 1   %s         # me2(1,1)
1 2   %s         # me2(1,2)
1 3   %s         # me2(1,3)
2 1   %s         # me2(2,1)
2 2   %s         # me2(2,2)
2 3   %s         # me2(2,3)
3 1   %s         # me2(3,1)
3 2   %s         # me2(3,2)
3 3   %s         # me2(3,3)
Block MSL2IN    #  
1 1   %s         # ml2(1,1)
1 2   %s         # ml2(1,2)
1 3   %s         # ml2(1,3)
2 1   %s         # ml2(2,1)
2 2   %s         # ml2(2,2)
2 3   %s         # ml2(2,3)
3 1   %s         # ml2(3,1)
3 2   %s         # ml2(3,2)
3 3   %s         # ml2(3,3)
Block RVM2LH1IN    #  
1   %s         # mlHd2(1)
2   %s         # mlHd2(2)
3   %s         # mlHd2(3)
Block MSQ2IN    #  
1 1   %s         # mq2(1,1)
1 2   %s         # mq2(1,2)
1 3   %s         # mq2(1,3)
2 1   %s         # mq2(2,1)
2 2   %s         # mq2(2,2)
2 3   %s         # mq2(2,3)
3 1   %s         # mq2(3,1)
3 2   %s         # mq2(3,2)
3 3   %s         # mq2(3,3)
Block MSU2IN    #  
1 1   %s         # mu2(1,1)
1 2   %s         # mu2(1,2)
1 3   %s         # mu2(1,3)
2 1   %s         # mu2(2,1)
2 2   %s         # mu2(2,2)
2 3   %s         # mu2(2,3)
3 1   %s         # mu2(3,1)
3 2   %s         # mu2(3,2)
3 3   %s         # mu2(3,3)
Block MV2IN    #  
1 1   %s         # mv2(1,1)
1 2   %s         # mv2(1,2)
1 3   %s         # mv2(1,3)
2 1   %s         # mv2(2,1)
2 2   %s         # mv2(2,2)
2 3   %s         # mv2(2,3)
3 1   %s         # mv2(3,1)
3 2   %s         # mv2(3,2)
3 3   %s         # mv2(3,3)
Block YVIN    #  
1 1   %s         # Yv(1,1) 
1 2   %s         # Yv(1,2)
1 3   %s         # Yv(1,3)
2 1   %s         # Yv(2,1)
2 2   %s         # Yv(2,2)
2 3   %s         # Yv(2,3)
3 1   %s         # Yv(3,1)
3 2   %s         # Yv(3,2)
3 3   %s         # Yv(3,3)
Block LAMIN    #  
1   %s           # lam(1)
2   %s           # lam(2)
3   %s           # lam(3)
Block KAPIN    #  
1 1 1   %s         # kap(1,1,1)
1 1 2   %s         # kap(1,1,2)
1 1 3   %s         # kap(1,1,3)
1 2 1   %s         # kap(1,2,1)
1 2 2   %s         # kap(1,2,2)
1 2 3   %s         # kap(1,2,3)
1 3 1   %s         # kap(1,3,1)
1 3 2   %s         # kap(1,3,2)
1 3 3   %s         # kap(1,3,3)
2 1 1   %s         # kap(2,1,1)
2 1 2   %s         # kap(2,1,2)
2 1 3   %s         # kap(2,1,3)
2 2 1   %s         # kap(2,2,1)
2 2 2   %s         # kap(2,2,2)
2 2 3   %s         # kap(2,2,3)
2 3 1   %s         # kap(2,3,1)
2 3 2   %s         # kap(2,3,2)
2 3 3   %s         # kap(2,3,3)
3 1 1   %s         # kap(3,1,1)
3 1 2   %s         # kap(3,1,2)
3 1 3   %s         # kap(3,1,3)
3 2 1   %s         # kap(3,2,1)
3 2 2   %s         # kap(3,2,2)
3 2 3   %s         # kap(3,2,3)
3 3 1   %s         # kap(3,3,1)
3 3 2   %s         # kap(3,3,2)
3 3 3   %s         # kap(3,3,3)
Block TDIN    #  
1 1   %s         # Td(1,1)
1 2   %s         # Td(1,2)
1 3   %s         # Td(1,3)
2 1   %s         # Td(2,1)
2 2   %s         # Td(2,2)
2 3   %s         # Td(2,3)
3 1   %s         # Td(3,1)
3 2   %s         # Td(3,2)
3 3   %s         # Td(3,3)
Block TEIN    # 
1 1   %s         # Te(1,1)
1 2   %s         # Te(1,2)
1 3   %s         # Te(1,3)
2 1   %s         # Te(2,1)
2 2   %s         # Te(2,2)
2 3   %s         # Te(2,3)
3 1   %s         # Te(3,1)
3 2   %s         # Te(3,2)
3 3   %s         # Te(3,3)
Block TUIN    #  
1 1   %s         # Tu(1,1)
1 2   %s         # Tu(1,2)
1 3   %s         # Tu(1,3)
2 1   %s         # Tu(2,1)
2 2   %s         # Tu(2,2)
2 3   %s         # Tu(2,3)
3 1   %s         # Tu(3,1)
3 2   %s         # Tu(3,2)
3 3   %s         # Tu(3,3)
Block TVIN    #  
1 1   %s         # Tv(1,1)
1 2   %s         # Tv(1,2)
1 3   %s         # Tv(1,3)
2 1   %s         # Tv(2,1)
2 2   %s         # Tv(2,2)
2 3   %s         # Tv(2,3)
3 1   %s         # Tv(3,1)
3 2   %s         # Tv(3,2)
3 3   %s         # Tv(3,3)
Block TLAMIN   #  
1   %s         # Tlam(1)
2   %s         # Tlam(2)
3   %s         # Tlam(3)
Block TKIN    #  
1 1 1   %s         # Tk(1,1,1)
1 1 2   %s         # Tk(1,1,2)
1 1 3   %s         # Tk(1,1,3)
1 2 1   %s         # Tk(1,2,1)
1 2 2   %s         # Tk(1,2,2)
1 2 3   %s         # Tk(1,2,3)
1 3 1   %s         # Tk(1,3,1)
1 3 2   %s         # Tk(1,3,2)
1 3 3   %s         # Tk(1,3,3)
2 1 1   %s         # Tk(2,1,1)
2 1 2   %s         # Tk(2,1,2)
2 1 3   %s         # Tk(2,1,3)
2 2 1   %s         # Tk(2,2,1)
2 2 2   %s         # Tk(2,2,2)
2 2 3   %s         # Tk(2,2,3)
2 3 1   %s         # Tk(2,3,1)
2 3 2   %s         # Tk(2,3,2)
2 3 3   %s         # Tk(2,3,3)
3 1 1   %s         # Tk(3,1,1)
3 1 2   %s         # Tk(3,1,2)
3 1 3   %s         # Tk(3,1,3)
3 2 1   %s         # Tk(3,2,1)
3 2 2   %s         # Tk(3,2,2)
3 2 3   %s         # Tk(3,2,3)
3 3 1   %s         # Tk(3,3,1)
3 3 2   %s         # Tk(3,3,2)
3 3 3   %s         # Tk(3,3,3)
Block MSOFTIN       #  
1   %s         # M1
3   %s         # M3
2   %s         # M2
 """ %
            (
                param['renormscale'], 
                param['invalpha'],
                param['alphas'], 
                MZ, 
                param['mb'], 
                param['mt'], 
                param['tanbeta'], 

                param['vR1'],
                param['vR1'] if InputParamFormat in {'B','A'} else  param['vR2'],
                param['vR1'] if InputParamFormat in {'B','A'} else  param['vR3'],
                param['vL1'],
                param['vL2'],
                param['vL3'],

                param['m2d11'], 
                0.0, 
                0.0, 
                0.0, 
                param['m2d22'] if InputParamFormat in {'B','A'} else param['m2d11'], 
                0.0, 
                0.0, 
                0.0, 
                param['m2d33'] if InputParamFormat in {'B','A'} else param['m2d11'], 

                param['m2e11'], 
                0.0, 
                0.0, 
                0.0, 
                param['m2e22'] if InputParamFormat in {'B','A'} else param['m2e11'], 
                0.0, 
                0.0, 
                0.0, 
                param['m2e33'] if InputParamFormat in {'B','A'} else param['m2e11'], 

                0.0, 
                0.0, 
                0.0, 
                0.0, 
                0.0,  
                0.0, 
                0.0, 
                0.0, 
                0.0, 

                0.0, 
                0.0, 
                0.0, 

                param['m2Q11'], 
                0.0, 
                0.0, 
                0.0,
                param['m2Q22'] if InputParamFormat in {'B','A'} else param['m2Q11'], 
                0.0, 
                0.0, 
                0.0,
                param['m2Q33'] if InputParamFormat in {'B','A'} else param['m2Q11'], 

                param['m2u11'], 
                0.0, 
                0.0, 
                0.0,
                param['m2u22'] if InputParamFormat in {'B','A'} else param['m2u11'], 
                0.0, 
                0.0, 
                0.0,
                param['m2u33'] if InputParamFormat in {'B','A'} else param['m2u11'], 

                0.0, 
                0.0, 
                0.0, 
                0.0, 
                0.0,  
                0.0, 
                0.0, 
                0.0, 
                0.0, 

                param['Yv11'], 
                param['Yv12'] if InputParamFormat == 'A' else 0.0, 
                param['Yv13'] if InputParamFormat == 'A' else 0.0, 
                param['Yv21'] if InputParamFormat == 'A' else 0.0, 
                param['Yv22'], 
                param['Yv23'] if InputParamFormat == 'A' else 0.0, 
                param['Yv31'] if InputParamFormat == 'A' else 0.0, 
                param['Yv32'] if InputParamFormat == 'A' else 0.0, 
                param['Yv33'], 

                param['lam1'], 
                param['lam2'], 
                param['lam3'], 

                param['kap111'], 
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0, 
                0.0,
                0.0,
                0.0,
                param['kap222'],
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0, 
                0.0,
                0.0,
                0.0,
                param['kap333'],

                0.0 if InputParamFormat in {'B','A'} else param['Td11'],
                0.0 if InputParamFormat in {'B','A'} else param['Td12'], 
                0.0 if InputParamFormat in {'B','A'} else param['Td13'], 
                0.0 if InputParamFormat in {'B','A'} else param['Td21'], 
                0.0 if InputParamFormat in {'B','A'} else param['Td22'], 
                0.0 if InputParamFormat in {'B','A'} else param['Td23'], 
                0.0 if InputParamFormat in {'B','A'} else param['Td31'], 
                0.0 if InputParamFormat in {'B','A'} else param['Td32'], 
                param['Td33'], 

                0.0 if InputParamFormat in {'B','A'} else param['Te11'],
                0.0 if InputParamFormat in {'B','A'} else param['Te12'], 
                0.0 if InputParamFormat in {'B','A'} else param['Te13'], 
                0.0 if InputParamFormat in {'B','A'} else param['Te21'], 
                0.0 if InputParamFormat in {'B','A'} else param['Te22'], 
                0.0 if InputParamFormat in {'B','A'} else param['Te23'], 
                0.0 if InputParamFormat in {'B','A'} else param['Te31'], 
                0.0 if InputParamFormat in {'B','A'} else param['Te32'], 
                param['Te33'], 

                0.0 if InputParamFormat in {'B','A'} else param['Tu11'],
                0.0 if InputParamFormat in {'B','A'} else param['Tu12'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tu13'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tu21'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tu22'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tu23'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tu31'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tu32'], 
                param['Tu33'], 

                param['Tv11'],
                0.0 if InputParamFormat in {'B','A'} else param['Tv12'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tv13'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tv21'],
                param['Tv22'] if InputParamFormat in {'B','A'} else param['Tv22'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tv23'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tv31'], 
                0.0 if InputParamFormat in {'B','A'} else param['Tv32'],
                param['Tv33'] if InputParamFormat in {'B','A'} else param['Tv33'], 

                param['Tlam1'], 
                param['Tlam1'] if InputParamFormat in {'B','A'} else param['Tlam2'],
                param['Tlam1'] if InputParamFormat in {'B','A'} else param['Tlam3'],

                param['Tk111'],
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0, 
                0.0,
                0.0,
                0.0,
                param['Tk222'] if InputParamFormat == 'A' else param['Tk111'],
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0, 
                0.0,
                0.0,
                0.0,
                param['Tk333'] if InputParamFormat == 'A' else param['Tk111'],

                param['M1'], 
                param['M3'] if InputParamFormat in {'B','A'} else (param['M1']*6), 
                param['M2'] if InputParamFormat in {'B','A'} else (param['M1']*2), 
             ))

        # Close file so that output is flushed.
        SLHAIN.close()
        #print ("SLHA input file:", SLHAIN.name)
        return SLHAIN.name



    # TO DO --- ADAPT TO MUNUSSM ---------------- FEB 12, 2024

    def readslha(self, input_file):						#MODIFIED
        """ Read an SLHA file with PySLHA.
        Populates masses, mu and neutralino mixings.
        If the point is unphysical,
        populate the mass blocks with zeros.

        Arguments:
        input_file -- Name of input file.

        Returns:

        """
        #print("Test before try except --- self.blocks, self.decays", input_file)
        #self.blocks, self.decays = pyslha213mod.readSLHAFile(input_file)

        #print('entering readSLHA --- Likelihood') 
        try:
            # Read the blocks in the SLHA file.
	        #print("Test before self.blocks, self.decays", input_file)
	        self.blocks, self.decays = pyslha213mod.readSLHAFile(input_file)
	        #print("Test after self.blocks, self.decays", input_file)
	    
        except Exception as e:
            # With expected running, shouldn't get any problems. But best
            # to be defensive. A missing mass block would cause an
            # exception, but e.g. stau LSP would not.
            self.physical = False
            print ('Caught trouble in the SLHA file:', e)

            # Still need to return data of correct length.			#MODIFIED
            self.masses = [0] * 59   # Number of particles masses + 1 (entry in the spheno output file)
            #self.mu = 0.                                           # COMMENTED. there is mu_eff as input parameter

            # NB neutralino mixing matrix is 10 by 10.                
            self.neutralinoRe = [0] * 100                             # OPTION 1: 100 (entire matrix)
            self.neutralinoIm = [0] * 100  

            # NB Higgs Scalar-Mixing-Matrix is 8 by 8: ZH.
            self.scalarhiggs = [0] * 64

            # NB Higgs Pseudo-Scalar-Mixing-Matrix is 8 by 8: ZA.
            self.pseudoscalarhiggs = [0] * 64

            # Pick out Neutrinos Mixing-Matrix: UV in SARAH {{FvL,conj[FvR]},{Fvm,UV}} , UVMIX in SPheno
            # Values to calculate observables in the Neutrino Sector, see def neutrinophysics.
            # Entries selected to get Normal Hierarchy
            self.UV31Re = 0.
            self.UV21Re = 0.
            self.UV32Re = 0.

            self.UV31Im = 0.
            self.UV21Im = 0.
            self.UV32Im = 0.

            # Pick out SPheno observables from the blocks: SPhenoLowEnergy, FlavorKitQFV, FlavorKitLFV
            # Values to constraint observables
            self.gm2 = 0.
            self.bsg = 0.
            self.btaunu = 0.
            self.TanBeta = 0.
            self.bdmumu = 0.
            self.bsmumu = 0.
            self.mueg = 0.
            self.mu3e = 0.
            self.LSPparticle = 0.

        else:
            # Save the mass block.
            self.masses = self.blocks['MASS'].entries
            
            ## Pick out neutralino mixing.
            #self.neutralinoRe = self.blocks['UVMIX'].entries      
            #self.neutralinoIm = self.blocks['IMUVMIX'].entries  

            ## Pick out Higgs Scalar-Mixing-Matrix: ZH 
            #self.scalarhiggs = self.blocks['SCALARMIX'].entries

            ## Pick out Higgs Pseudo-Scalar-Mixing-Matrix: ZA 
            #self.pseudoscalarhiggs = self.blocks['PSEUDOSCALARMIX'].entries     


            # Pick out Neutrinos Mixing-Matrix: UV in SARAH {{FvL,conj[FvR]},{Fvm,UV}} , UVMIX in SPheno
            # Values to calculate observables in the Neutrino Sector, see def neutrinophysics.
            # Entries selected to get Normal Hierarchy

            self.UV31Re = self.blocks['UVMIX'].entries[3,1]
            self.UV21Re = self.blocks['UVMIX'].entries[2,1] 
            self.UV32Re = self.blocks['UVMIX'].entries[3,2] 

            self.UV31Im = self.blocks['IMUVMIX'].entries[3,1] 
            self.UV21Im = self.blocks['IMUVMIX'].entries[2,1]
            self.UV32Im = self.blocks['IMUVMIX'].entries[3,2] 


            # Pick out SPheno observables from the blocks: SPhenoLowEnergy, FlavorKitQFV, FlavorKitLFV
            # Values to constraint observables
            self.gm2 = self.blocks['SPHENOLOWENERGY'].entries[21]
            self.bsg = self.blocks['FLAVORKITQFV'].entries[200]       # BR(B->X_s gamma)
            self.btaunu  = self.blocks['FLAVORKITQFV'].entries[502]   # BR(B->tau nu)
            self.TanBeta = self.blocks['MINPAR'].entries[3]           # BR(B->tau nu) 
            self.bdmumu = self.blocks['FLAVORKITQFV'].entries[4004]   # BR(B^0_d->mu mu)
            self.bsmumu = self.blocks['FLAVORKITQFV'].entries[4006]   # BR(B^0_s->mu mu)
            self.mueg   = self.blocks['FLAVORKITLFV'].entries[701]    # BR(mu->e gamma)
            self.mu3e   = self.blocks['FLAVORKITLFV'].entries[901]    # BR(mu->3e)



    def micromegas(self, input_file):				# NOT YET ADAPTED TO MUNUSSM
        """Call micrOMEGAs to obtain DM predictions
        for model.

        Arguments:
        input_file -- Name of input file.

        Returns:

        """
        self.microfile = RunProgram(
            './OmegaMultiNestDD',
            '../micromegas_4.3.4/SPhenomunuSSM',
            input_file)
        self.CheckProgram(self.microfile, ["error"])
        if self.physical:
            calcoh2 = self.ReadParameter(                       # MODIFIED
                self.microfile,
                'Xf=',
                split='Omega h^2=')
            getLSPmass = self.ReadParameter(				    # MODIFIED
                self.microfile,
                'Dark matter candidate',
                split='mass=')
  
            if constraintRelicDensity == 'True':                # ADDED
                self.constraint['oh2'].theory = calcoh2
            #else:
            #    print 'oh2', calcoh2

            if constraintRelicDensity == 'False' and addcubeoh2 == 'True':  # ADDED
                self.constraint['oh2'].theory = calcoh2
                self.constraint['oh2'].loglike = 0.

            # Spin-Independent WIMP-proton cross-section in pb
            sigmapSI = self.ReadParameterMiddle(				# MODIFIED
                self.microfile,' proton  SI',split='  SD ')
            # Spin-Independent WIMP-neutron cross-section in pb
            sigmanSI = self.ReadParameterMiddle(				# MODIFIED
                self.microfile,' neutron SI',split='  SD ')
            # Spin-Dependent WIMP-proton cross-section in pb
            sigmapSD = self.ReadParameter(				        # MODIFIED
                self.microfile,' proton  SI',split='  SD ')
            # Spin-Dependent WIMP-neutron cross-section in pb
            sigmanSD = self.ReadParameter(			        	# MODIFIED
                self.microfile,' neutron SI',split='  SD ')


            if self.LSPparticle == 1000022:
                LSPname = 'Chi_1'
            elif self.LSPparticle == 1000012:
                LSPname = 'SvRe_1'
            elif self.LSPparticle == 3000012:
                LSPname = 'SvIm_1'
            else:
                LSPname = 'SLHA # ' + str(self.LSPparticle)

            #print ' LSP', LSPname, ' ; mass ', getLSPmass, 'GeV'
            #print ' oh2 ', calcoh2
            #print ' sigmapSI ', sigmapSI, 'pb (*E-36 to cm2)'
            #print ' sigmanSI ', sigmanSI, 'pb'
            #print ' sigmapSD ', sigmapSD, 'pb'
            #print ' sigmanSD ', sigmanSD, 'pb'

            if addSISDpn == 'True':  # ADDED
                self.constraint['SIproton'].theory = sigmapSI
                self.constraint['SIproton'].loglike = 0.

                self.constraint['SIneutron'].theory = sigmanSI
                self.constraint['SIneutron'].loglike = 0.

                self.constraint['SDproton'].theory = sigmapSD
                self.constraint['SDproton'].loglike = 0.

                self.constraint['SDneutron'].theory = sigmanSD
                self.constraint['SDneutron'].loglike = 0.



    def higgstools(self, input_file, HTrootname):
        """Call Higgstools to find whether points is excluded by
        LEP, Tevatron, LHC Higgs searches.

        Arguments:
        input_file -- Name of input file.

        Returns:

        """

        # by DEK ... change '/home/dkpatcha/anaconda3/bin/python3' to your python executable
        self.higgstoolsfile = RunProgramHiggsTools(
            'python3', 
            '../higgstools-main/buildpy',
            ['HiggsToolsWithHB5Datafiles.py', HTrootname]
            )
        self.CheckProgram(self.higgstoolsfile, ["error"])
        
        
        if self.physical:
            print () 
            print ('self.higgstoolsfile = ', self.higgstoolsfile) 
            print () 
            cols_to_use = ["m_H1",   "m_H2",   "m_H3",   "m_H4",   "m_H5",   "m_H6",   "m_H7", 
                           "m_H8",   "m_H9",   "m_H10",   "m_H11",   "m_H12",   "m_H13",   "m_H14",
                           "m_H15",   "m_H1+",   "m_H2+",   "m_H3+",   "m_H4+",   "m_H5+",   "m_H6+", 
                           "m_H7+",   "HSchisq",   "HBresult"]

            dataout = pd.read_csv(self.higgstoolsfile, sep='\t', usecols=cols_to_use)

            Higgscalc = -0.5 * dataout["HSchisq"].iloc[0]    
            HiggsBounds = dataout["HBresult"].iloc[0]  # dataout.HBresult
            
            if constraintHiggsBS == 'True':         # ADDED notice its a .loglike
                self.constraint['Higgs'].loglike = Higgscalc
                self.constraint['Higgspvalue'].theory = 0.0  # not computed in this version
                self.constraint['Higgspvalue'].loglike = 0.0
                self.constraint['HiggsBounds'].theory = HiggsBounds
                if float(HiggsBounds) == 1.0:
                    self.constraint['HiggsBounds'].loglike = -0.0
                else:
                    self.constraint['HiggsBounds'].loglike = -1e101
            else:
                print ('chi2:Higgs ', -2.0*Higgscalc )             # as said before
                print ('HBresult (1: allowed, 0: excluded) ', HiggsBounds)

            #print ('  ')
            print (' chi2:Higgs = ', -2.0*Higgscalc)            # as said before
            print (' HBresult (1: allowed, 0: excluded) =', HiggsBounds)
            #print ('  ')



    def higgssignals(self, input_file): # NOT YET ADAPTED TO MUNUSSM
        """Call HiggsSignals to find whether points is excluded by
        LEP, Tevatron, LHC Higgs searches.

        Arguments:
        input_file -- Name of input file.

        Returns:

        """
        self.higgssignalsfile = RunProgram(
            './SuperPy',
            '../HiggsSignals-2.2.3beta/example_programs',
            input_file)
        self.CheckProgram(self.higgssignalsfile, ["error"])
        if self.physical:
            # -0.5* because you have a -2* in the CUBE definition above
            Higgscalc = -0.5 * self.ReadParameter(self.higgssignalsfile, 'chi^2 (total peak+STXS) =')     
            # p-value given by HS
            Higgspvaluecalc = self.ReadParameter(self.higgssignalsfile, 'Probability peak+STXS',split=') = ') 
            HiggsBounds = self.ReadParameter(self.higgssignalsfile, 'HB result (1: allowed,',split=') = ')
            if constraintHiggsBS == 'True':                    # ADDED notice its a .loglike
                self.constraint['Higgs'].loglike = Higgscalc
                self.constraint['Higgspvalue'].theory = Higgspvaluecalc
                self.constraint['Higgspvalue'].loglike = 0.0
                self.constraint['HiggsBounds'].theory = HiggsBounds
                if float(HiggsBounds) == 1.0:
                    self.constraint['HiggsBounds'].loglike = -0.0
                else:
                    self.constraint['HiggsBounds'].loglike = -1e101
            else:
                print ('chi2:Higgs ', -2 *Higgscalc)              # as said before
                print ('p-value:HiggsHS ', Higgspvaluecalc)
                print ('HBresult (1: allowed) ', HiggsBounds)

            print (' chi2:Higgs ', -2 *Higgscalc)            # as said before
            print (' p-value:HiggsHS ', Higgspvaluecalc)
            print (' HBresult (1: allowed) ', HiggsBounds)


    def ddcalc(self, input_file):
        """Call DDCalc to get SI and SD constraints

        Arguments:
        input_file -- Name of input file.

        Returns:

        """

        if self.physical:    # READ SI and SD cross sections from the file generated by MicrOmegas 
            # Relic Density calculated by MicrOmegas
            calcoh2 = self.ReadParameter(                       # MODIFIED
                self.microfile,'Xf=',split='Omega h^2=')
            # LSP mass in GeV
            getLSPmass = self.ReadParameter(				    # MODIFIED
                self.microfile,'Dark matter candidate',split='mass=')
            # Spin-Independent WIMP-proton cross-section in pb
            sigmapSI = self.ReadParameterMiddle(				# MODIFIED
                self.microfile,' proton  SI',split='  SD ')
            # Spin-Independent WIMP-neutron cross-section in pb
            sigmanSI = self.ReadParameterMiddle(				# MODIFIED
                self.microfile,' neutron SI',split='  SD ')
            # Spin-Dependent WIMP-proton cross-section in pb
            sigmapSD = self.ReadParameter(				        # MODIFIED
                self.microfile,' proton  SI',split='  SD ')
            # Spin-Dependent WIMP-neutron cross-section in pb
            sigmanSD = self.ReadParameter(			        	# MODIFIED
                self.microfile,' neutron SI',split='  SD ')

            # normalized: as the DM relic density calculated by micromegas can be different from the measured value, the constraints will change. This is needed to be informed to DDCalc: using an effective cross section, or a different local DM density
            #sigSIp_normalized = sigmapSI*(calcoh2/0.1188)
            #sigSIn_normalized = sigmanSI*(calcoh2/0.1188)
            #sigSDp_normalized = sigmapSD*(calcoh2/0.1188)
            #sigSDn_normalized = sigmanSD*(calcoh2/0.1188)
            localrho_normalized = 0.3*(calcoh2/0.1188)

            DDCalcAux = open('../DDCalc_v2.0.0/multinestpoint.dat', 'w+')   # Create a file with the relevant info for DDCalc
            DDCalcAux.write('%s\n' % getLSPmass)
            #DDCalcAux.write('%s\n' % sigSIp_normalized)
            #DDCalcAux.write('%s\n' % sigSIn_normalized)
            #DDCalcAux.write('%s\n' % sigSDp_normalized)
            #DDCalcAux.write('%s' % sigSDn_normalized)
            DDCalcAux.write('%s\n' % sigmapSI)
            DDCalcAux.write('%s\n' % sigmanSI)
            DDCalcAux.write('%s\n' % sigmapSD)
            DDCalcAux.write('%s\n' % sigmanSD)
            DDCalcAux.write('%s' % localrho_normalized)
            DDCalcAux.close()


            # RUNS DDCalc for every efficiency file, and creates a new file, "self.DDCalcfile" with all the outputs

            input_string = './DDCalc_MultiNest'

            self.DDCalcfile = RunProgramDDCalc(
                input_string,
                '../DDCalc_v2.0.0')
            self.CheckProgram(self.DDCalcfile, ["error"])


            # READ the loglikelihood and p-value calculated 

            Log_LDDCalcXe17 = self.ReadParameter(self.DDCalcfile,'XENON1T_2017 Log(likelihood):',split=':')
            pvalDDCalcXe17 = self.ReadParameter(self.DDCalcfile,'XENON1T_2017 p-value(no bg subtraction):',split=':')

            Log_LDDCalcCR17 = self.ReadParameter(self.DDCalcfile,'CRESST_2017 Log(likelihood):',split=':')
            pvalDDCalcCR17 = self.ReadParameter(self.DDCalcfile,'CRESST_2017 p-value(no bg subtraction):',split=':')

            Log_LDDCalcDRW = self.ReadParameter(self.DDCalcfile,'DARWIN Log(likelihood):',split=':')
            pvalDDCalcDRW = self.ReadParameter(self.DDCalcfile,'DARWIN p-value(no bg subtraction):',split=':')

            Log_LDDCalcP6017 = self.ReadParameter(self.DDCalcfile,'PICO_60_2017 Log(likelihood):',split=':')
            pvalDDCalcP6017 = self.ReadParameter(self.DDCalcfile,'PICO_60_2017 p-value(no bg subtraction):',split=':')

            Log_LDDCalcXe18 = self.ReadParameter(self.DDCalcfile,'XENON1T_2018 Log(likelihood):',split=':')
            pvalDDCalcXe18 = self.ReadParameter(self.DDCalcfile,'XENON1T_2018 p-value(no bg subtraction):',split=':')


            # ADD the loglikelihood and p-value to MultiNest

            if DDCalcXENON1T_2017 == 'True':                    # ADDED # Adds the loglike to the cube
                print (' XENON1T_2017 loglike ', Log_LDDCalcXe17)
                print (' XENON1T_2017 p-value ', pvalDDCalcXe17)
                self.constraint['DDCalcXENON1T_2017'].loglike = Log_LDDCalcXe17
                if StorepvalueDDCalc == 'True':          # Adds the p-value to the saved values, but not to the cube
                    self.constraint['pvalueDDCalcXENON1T_2017'].theory = pvalDDCalcXe17
                    self.constraint['pvalueDDCalcXENON1T_2017'].loglike = 0.
            else:
                if DoDDCalcXENON1T_2017 == 'True':
                    print (' XENON1T_2017 loglike ', Log_LDDCalcXe17)
                    print (' XENON1T_2017 p-value ', pvalDDCalcXe17)


            if DDCalcCRESST_2017 == 'True':                    # ADDED # Adds the loglike to the cube
                print (' CRESST_2017 loglike ', Log_LDDCalcCR17)
                print (' CRESST_2017 p-value ', pvalDDCalcCR17)
                self.constraint['DDCalcCRESST_2017'].loglike = Log_LDDCalcCR17
                if StorepvalueDDCalc == 'True':          # Adds the p-value to the saved values, but not to the cube
                    self.constraint['pvalueDDCalcCRESST_2017'].theory = pvalDDCalcCR17
                    self.constraint['pvalueDDCalcCRESST_2017'].loglike = 0.
            else:
                if DoDDCalcCRESST_2017 == 'True':
                    print (' CRESST_2017 loglike ', Log_LDDCalcCR17)
                    print (' CRESST_2017 p-value ', pvalDDCalcCR17)


            if DDCalcDARWIN == 'True':                    # ADDED # Adds the loglike to the cube
                print (' DARWIN loglike ', Log_LDDCalcDRW)
                print (' DARWIN p-value ', pvalDDCalcDRW)
                self.constraint['DDCalcDARWIN'].loglike = Log_LDDCalcDRW
                if StorepvalueDDCalc == 'True':          # Adds the p-value to the saved values, but not to the cube
                    self.constraint['pvalueDDCalcDARWIN'].theory = pvalDDCalcDRW
                    self.constraint['pvalueDDCalcDARWIN'].loglike = 0.
            else:
                if DoDDCalcDARWIN == 'True':
                    print (' DARWIN loglike ', Log_LDDCalcDRW)
                    print (' DARWIN p-value ', pvalDDCalcDRW)


            if DDCalcPICO_60_2017 == 'True':                    # ADDED # Adds the loglike to the cube
                print (' PICO_60_2017 loglike ', Log_LDDCalcP6017)
                print (' PICO_60_2017 p-value ', pvalDDCalcP6017)
                self.constraint['DDCalcPICO_60_2017'].loglike = Log_LDDCalcP6017
                if StorepvalueDDCalc == 'True':          # Adds the p-value to the saved values, but not to the cube
                    self.constraint['pvalueDDCalcPICO_60_2017'].theory = pvalDDCalcP6017
                    self.constraint['pvalueDDCalcPICO_60_2017'].loglike = 0.
            else:
                if DoDDCalcPICO_60_2017 == 'True':
                    print (' PICO_60_2017 loglike ', Log_LDDCalcP6017)
                    print (' PICO_60_2017 p-value ', pvalDDCalcP6017)


            if DDCalcXENON1T_2018 == 'True':                    # ADDED # Adds the loglike to the cube
                print (' XENON1T_2018 loglike ', Log_LDDCalcXe18)
                print (' XENON1T_2018 p-value ', pvalDDCalcXe18)
                self.constraint['DDCalcXENON1T_2018'].loglike = Log_LDDCalcXe18
                if StorepvalueDDCalc == 'True':          # Adds the p-value to the saved values, but not to the cube
                    self.constraint['pvalueDDCalcXENON1T_2018'].theory = pvalDDCalcXe18
                    self.constraint['pvalueDDCalcXENON1T_2018'].loglike = 0.
            else:
                if DoDDCalcXENON1T_2018 == 'True':
                    print (' XENON1T_2018 loglike ', Log_LDDCalcXe18)
                    print (' XENON1T_2018 p-value ', pvalDDCalcXe18)




    def smodels(self, input_file):    # NOT YET ADAPTED TO MUNUSSM
        """Call SModelS to find whether points is excluded LHC searches.

        Arguments:
        input_file -- Name of input file.

        Returns:

        """
        if self.physical:
            os.chdir('../../smodels-v1.1.1')                    # REMEMBER TO CORRECT THIS PATH

                # Generate the cross section. The XSECTION blocks in the SLHA file
            xseccomputer ='./smodelsTools.py xseccomputer -p -f '
            #input_file ='inputFiles/slha/modSPhenoOutput.slha'
            xseccomputer_file = ' > xseccomputerPRINT'
            generateXS= xseccomputer + input_file + xseccomputer_file

            os.system(generateXS)

                # Check if everything is right in the SLHA file                 # GOTTA ADD IFs everywhere
            slhacheck ='./smodelsTools.py slhachecker -m 0.001 -s 0.01 -f '
            checker_file = ' > checker'
            checkSLHA= slhacheck + input_file + checker_file

            os.system(checkSLHA)


            analysisLHC ='./runSModelS.py -f '
            param_ini = ' -p parameters.ini -o ./ -v warning'
            analysis_file = ' > SModelSOutput_MultiNest'
            command4= analysisLHC + input_file + param_ini + analysis_file

            os.system(command4)


            os.chdir('../superpy-munussm-v.0.0/SPhenomunuSSM')          # REMEMBER TO CORRECT THIS PATH



    def neutrinophysics(self):          # ADDED
        """Routine to calculate neutrino physics observables. Experimental constraints applied. Change settings in "InputFile"

        Arguments:
        Input file sets

        Returns:

        """
        GeVtoeV = 1000000000        # GeV to eV factor

        if self.physical:

            # m_2^2 - m_1^2
            calcdm2 = (self.masses[14] * GeVtoeV)**2 - (self.masses[12] * GeVtoeV)**2 
            # m_3^2 - (m_2^2 + m_1^2)/2
            calcDm2 = (self.masses[16] * GeVtoeV)**2 - (((self.masses[14] * GeVtoeV)**2 + (self.masses[12] * GeVtoeV)**2)/2) 

            FAC1 = NP.abs(NP.sqrt(self.UV31Re**2 + self.UV31Im**2))  
            FAC2 = NP.abs(NP.sqrt(1. - FAC1*FAC1))
            FAC3 = NP.abs(NP.sqrt(self.UV21Re**2 + self.UV21Im**2)) 
            FAC4 = NP.abs(NP.sqrt(self.UV32Re**2 + self.UV32Im**2)) 

            calcsin2theta13 = FAC1**2
            calcsin2theta12 = (FAC3/FAC2)**2
            calcsin2theta23 = (FAC4/FAC2)**2

            '''
            ! ----------- Neutrino masses and mixing --------
            FAC1 = NP.abs(NP.sqrt(self.UV31Re**2 + self.UV31Im**2))  # abs(sqrt(Real(UV(3,1),dp)**2 + Aimag(UV(3,1))**2))
            FAC2 = NP.abs(NP.sqrt(1. - FAC1*FAC1))
            FAC3 = NP.abs(NP.sqrt(self.UV21Re**2 + self.UV21Im**2))  # abs(sqrt(Real(UV(2,1),dp)**2 + Aimag(UV(2,1))**2))
            FAC4 = NP.abs(NP.sqrt(self.UV32Re**2 + self.UV32Im**2))  # abs(sqrt(Real(UV(3,2),dp)**2 + Aimag(UV(3,2))**2))
            sinsqth13 = FAC1**2
            sinsqth12 = (FAC3/FAC2)**2
            sinsqth23 = (FAC4/FAC2)**2
            '''
            calcsumneutrinos = (math.fabs(self.masses[12]) + math.fabs(self.masses[14]) + math.fabs(self.masses[16])) * GeVtoeV

            if constraintNeutrinoPhysics == 'True':           # ADDED
                self.constraint['dm2'].theory = calcdm2
                self.constraint['Dm2'].theory = calcDm2
                self.constraint['sin2theta12'].theory = calcsin2theta12
                self.constraint['sin2theta13'].theory = calcsin2theta13
                self.constraint['sin2theta23'].theory = calcsin2theta23

                print (' dm2', calcdm2)
                print (' Dm2', calcDm2)
                print (' sin2theta12', calcsin2theta12)
                print (' sin2theta13', calcsin2theta13)
                print (' sin2theta23', calcsin2theta23)

            else:
                print (' dm2', calcdm2)
                print (' Dm2', calcDm2)
                print (' sin2theta12', calcsin2theta12)
                print (' sin2theta13', calcsin2theta13)
                print (' sin2theta23', calcsin2theta23)


            if constraintSumNeutrinos == 'True':           # ADDED
                self.constraint['sumneutrinos'].theory = calcsumneutrinos
            else:
                print (' sumneutrinos', calcsumneutrinos)



    def sphenoobservables(self):                                  # ADDED
        """Routine to get the Selected Spheno observables. Experimental constraints applied. Change settings in "InputFile"

        Arguments:
        Input file sets

        Returns:

        """
        if self.physical:

            calcgm2 = self.gm2
            calcbsg = self.bsg
#            calcbmunu = self.bmunu
#            calcbtaunu = self.btaunu
#### ALTERNATIVE TO self.btaunu  ####
# Ratio BR / BR_SM      #https://arxiv.org/pdf/hep-ph/0605012.pdf
# BR        # http://pdg.lbl.gov/2017/download/rpp2016-Chin.Phys.C.40.100001.pdf     pag 1144
# BR_SM  #https://arxiv.org/pdf/1503.05613.pdf    https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.110.131801
            calcbtaunu = (1-(((5.27931**2)/(self.masses[37]**2))*((self.TanBeta**2)/(1+(0.01*self.TanBeta)))))**2 * 0.000075
            calcbdmumu = self.bdmumu
            calcbsmumu = self.bsmumu
            calcmueg = self.mueg
            calcmu3e = self.mu3e

            if constraintSPhenoObservables == 'True':           # ADDED
                self.constraint['gm2'].theory = calcgm2
                self.constraint['bsg'].theory = calcbsg
#                self.constraint['bmunu'].theory = calcbmunu
                self.constraint['btaunu'].theory = calcbtaunu
                self.constraint['bdmumu'].theory = calcbdmumu
                self.constraint['bsmumu'].theory = calcbsmumu
                self.constraint['mueg'].theory = calcmueg
                self.constraint['mu3e'].theory = calcmu3e

            else:
                print (' gm2', calcgm2)
                print (' bsg', calcbsg)
#                print( 'bmunu', calcbmunu)
                print (' btaunu', calcbtaunu)
                print (' bdmumu', calcbdmumu)
                print (' bsmumu', calcbsmumu)
                print (' mueg', calcmueg)
                print (' mu3e', calcmu3e)



    def delete(self):
        """Deletes all the auxiliary files generated by the subprograms. Change settings in "InputFile"

        Arguments:
        Input file sets

        Returns:

        """
        if self.physical:
            if DeleteSLHAinput == 'True':                  # Delete SLHA input auxiliary file
                os.remove(self.SLHAIN)
                print (' SLHA input deleted')
            if DeleteSPheno == 'True':                     # Delete SPheno auxiliary file
                os.remove(self.SLHA)
                print (' SPheno output file deleted')
            if MicrOmegasOK == 'True':
                if DeleteMicrOmegas == 'True':              # Delete MicrOmegas auxiliary file
                    os.remove(self.microfile)
                    print (' micrOmegas output file deleted')
            if HiggsSignalsOK == 'True':
                if DeleteHiggsSignals == 'True':            # Delete HiggsSignals auxiliary file
                    os.remove(self.higgssignalsfile)
                    print (' HS-HB output file deleted')

            if HiggsToolsOK == 'True':
                if DeleteHiggsTools == 'True':              # Delete HiggsSignals auxiliary file
                    os.remove(self.higgstoolsfile)
                    print (' Higgstools output file deleted')

            if SModelSOK == 'True':
                if DeleteSModelS == 'True':                  # Delete SModelS auxiliary file
                    print (' SModelS not activated yet')
            if DDCalcOK == 'True':
                if DeleteDDCalc == 'True':                   # Delete DDCalc auxiliary file
                    os.remove(self.DDCalcfile)
                    print (' DDCalc output file deleted')
             
        else:
            if DeleteSLHAinput == 'True':                    # Delete SLHA input auxiliary file
                os.remove(self.SLHAIN)
                print (' SLHA input deleted')

            if DeleteSPheno == 'True':                    # Delete SPheno auxiliary file
                try:                                      # Check if the file was generated
                    os.remove(self.SLHA)
                    print (' SPheno output file deleted')
                except:
                    print(" SPheno output file was not generated")

            if MicrOmegasOK == 'True':
                if DeleteMicrOmegas == 'True':                  # Delete MicrOmegas auxiliary file
                    try:                                        # Check if the file was generated
                        os.remove(self.microfile)
                        print (' micrOmegas output file deleted')
                    except:
                        print(" micrOmegas output file was not generated")

            if HiggsSignalsOK == 'True':
                if DeleteHiggsSignals == 'True':                # Delete HiggsSignals auxiliary file
                    try:                                        # Check if the file was generated
                        os.remove(self.higgssignalsfile)
                        print (' HS-HB output file deleted')
                    except:
                        print(" HS-HB output file was not generated")

            if DDCalcOK == 'True':
                if DeleteDDCalc == 'True':                    # Delete DDCalc auxiliary file
                    try:                                      # Check if the file was generated
                        os.remove(self.DDCalcfile)
                        print (' DDCalc output file deleted')
                    except:
                        print(" DDCalc output file was not generated")



    def SetLogLike(self):
        """ Loops over the constraints, calculates each log likelihood
        and sums the log likelihoods.
        Arguments:

        Returns:

        """
        if self.physical:
            # Set the loglikelhoods from the constraints's SetLogLike
            # functions.
            self.loglike = 0.
            for name in self.constraint.keys():
                if self.constraint[name].apply:
                    self.loglike = self.loglike + \
                        self.constraint[name].SetLogLike()
        else:
            # If the point is unphysical, all loglikelihoods = biggest possible.
            # Logzero is -1e101.
            for name in self.constraint.keys():
                self.constraint[name].loglike = -1e101
            self.loglike = -1e101

    def ReadParameter(self, filename, start, split=None):
        """ Read a parameter from a file. The parameter is right after the string "split".

        Arguments:
        filename -- File to be read from.
        start -- String of beginning of line to read.
        split -- string on which to split the line. The parameter to be read must be right after this string.

        Returns:
        Parameter retrieved from the file.

        """
        if split is None:
            split = start[-1]

        start = start.strip()
        try:
            for line in open(filename, 'r'):
                if line.lstrip().startswith(start):
                    parameter = float(line.split(split)[-1])
                    if math.isnan(parameter):
                        self.physical = False
                        return 999.
                    else:
                        return parameter
        except:
            self.physical = False
            return 999.

        else:
            self.physical = False
            return 999.



    def ReadParameterMiddle(self, filename, start, split=None):
        """ Read a parameter from a file. The parameter is between the strings "start" and "split".

        Arguments:
        filename -- File to be read from.
        start -- String of beginning of line to read. The parameter to be read must be right after this string.
        split -- string on which to split the line. The parameter to be read must be right before this string.

        Returns:
        Parameter retrieved from the file.

        """
        if split is None:
            split = start[-1]

        aux = len(start)
        start = start.strip()
        try:
            for line in open(filename, 'r'):
                if line.lstrip().startswith(start):
                    parameter = float(line.split(split)[0][aux:])
                    if math.isnan(parameter):
                        self.physical = False
                        return 999.
                    else:
                        return parameter
        except:
            self.physical = False
            return 999.

        else:
            self.physical = False
            return 999.





    def CheckProgram(self, filename, errors):
        """ Check whether a program contained errors.

        Arguments:
        filename -- File to be read from.
        errors -- Keywords that indicate an error.

        """
        if filename is None:
            print ("Error - no output returned.", filename)
            self.physical = False
            return

        for line in open(filename, 'r'):
            for word in errors:
                if word in line:
                    print( "Error in program output.")
                    print( line)
                    self.physical = False
                    return

#########################################################################

# These classes are various likelihoods for experimental measurements.
# I have included Gaussians, error functions etc. They all have a SetLogLike
# function, from which the loglikelhood is calculated.


@Cube.memoize
class GaussConstraint:

    """ A Gaussian likelihood constraint. """

    def __init__(self, mu, sigma, tau=0, apply=True):
        """ Initializes a Gaussian constraint.
        Arguments:
        mu -- Experimental mean of measurement.
        sigma -- Experimental variance of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.mu = mu
        self.sigma = sigma
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Caclulate Gaussian likelihood from model's predictions and
        data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        self.loglike = - math.pow(self.mu - self.theory, 2) /\
            (2 * (math.pow(self.sigma, 2) +
                  math.pow(self.tau, 2)))
        return self.loglike


@Cube.memoize
class GaussConstraintFractionalTau:

    """ A Gaussian likelihood constraint, with a fractional
    theoretical error.
    """

    def __init__(self, mu, sigma, tau, apply=True):
        """ Initializes a Gaussian constraint, with a fractional
        theoretical error.
        Arguments:
        mu -- Experimental mean of measurement.
        sigma -- Experimental variance of measurement.
        tau -- Fractional theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.mu = mu
        self.sigma = sigma
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Caclulates Gaussian likelihood from model's predictions and
        data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # NB that here theory error = tau * theoretical value.
        self.loglike = - math.pow(self.mu - self.theory, 2) /\
            (2 * (math.pow(self.sigma, 2) +
                  math.pow(self.tau * self.theory, 2)))
        return self.loglike


@Cube.memoize
class UpperConstraint:

    """ Upper bound, Gaussian error function constraint. """

    def __init__(self, limit, tau=None, apply=True):
        """ Initializes an upper bound constraint .
        Arguments:
        limit -- Upper bound of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.limit = limit
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Calcualte upper bound likelihood with Gaussian error
        function from model's predictions and data.
        Arguments:

        Returns: The loglike from this constraint.

        """

        if self.tau:
          # NB that erfc is a complementary error function.
          try:
              self.loglike = math.log(
                  0.5 * scipy.special.erfc((self.theory - self.limit) / (math.sqrt(2) * self.tau)))
          # Sometimes we are trying to take log zero.
          except ValueError:
              self.loglike = -1e101

        else:
          self.loglike = (self.theory > self.limit) * -1e101

        return self.loglike


@Cube.memoize
class UpperConstraintFractionalTau:

    """ Upper bound, Gaussian error function constraint, with fractional
    theory error.
    """

    def __init__(self, limit, tau, apply=True):
        """ Initializes an upper bound constraint .
        Arguments:
        limit -- Upper bound of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.limit = limit
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Calcualte upper bound likelihood with Gaussian error
        function from model's predictions and data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # NB that erfc is a complementary error function.
        # NB that here theory error = tau * theoretical value.
        try:
            self.loglike = math.log(0.5 *
                                    scipy.special.erfc((self.theory -
                                                        self.limit) /
                                                       (math.sqrt(2) *
                                                        self.tau *
                                                        self.theory)))
        # Sometimes we are trying to take log zero.
        except ValueError as ZeroDivisionError:
            self.loglike = -1e101
        return self.loglike


@Cube.memoize
class LowerConstraint:

    """ Lower bound, Gaussian error function constraint. """

    def __init__(self, limit, tau=None, apply=True):
        """ Initializes a lower bound constraint.
        Arguments:
        limit -- Lower bound of measurement.
        tau -- Theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.limit = limit
        self.tau = tau
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ Calcualte upper bound likelihood with Gaussian error
        function from model's predictions and data.
        Arguments:

        Returns: The loglike from this constraint.

        """
        if self.tau:
          # NB that erf is an error function.
          try:
              self.loglike = math.log(0.5 * (1 + scipy.special.erf((
                  self.theory - self.limit) /
                  (math.sqrt(2) * self.tau))))
          # Sometimes we are trying to take log zero.
          except ValueError:
              self.loglike = -1e101
        else:
          self.loglike = (self.theory < self.limit) * -1e101

        return self.loglike


@Cube.memoize
class InterpolateUpperConstraint:

    """ Interpolate am upper 2D limit from a data file. """

    def __init__(self, file, tau=None, apply=True):
        """ Initializes a 2D upper bound constraint.
        The constraint is on the y-co-ordinate, and is calculated as
        a function of x.

        Arguments:
        file -- Name of the data file.
        tau -- Fractional theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.file = os.path.abspath(file)
        self.tau = tau
        self.apply = apply

        # The x and y theory values.
        self.theoryx = 0.
        self.theory = 999.

        # The derived limit on the y-parameter, which is a function of x.
        self.limit = 0.

        self.loglike = -1e101
        self.arg = self.args

        # Import the data from file. Do this here so
        # we only do it once.
        # The data should be x, y.
        self.data = NP.loadtxt(self.file)

    def SetLogLike(self):
        """ Calcualte lower bound with interpolation
        function, then apply an error function.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # Interpolate the y-value of the limit
        # corresponding to the theory x-value.
        self.limit = NP.interp(self.theoryx, self.data[0],
                               self.data[1], left=None, right=None)


        if self.tau:

          # Now calcualte likelihood with Gaussian error function - erf.
          like = 0.5 - 0.5 * \
              scipy.special.erf((self.theory - self.limit) / (2. ** 0.5 * self.tau * self.theory))
          try:
              self.loglike = math.log(like)
          # Sometimes we are trying to take log zero.
          except ValueError:
              self.loglike = -1e101

        else:

          self.loglike = (self.theory > self.limit) * -1e101

        return self.loglike


@Cube.memoize
class InterpolateLowerConstraint:

    """ Interpolate a lower 2D limit from a data file. """

    def __init__(self, file, tau=None, apply=True):
        """ Initializes a 2D lower bound constraint.
        The constraint is on the y-co-ordinate, and is calculated as
        a function of x.

        Arguments:
        file -- Name of the data file.
        tau -- Fractional theoretical error in calculation.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.file = os.path.abspath(file)
        self.tau = tau
        self.apply = apply

        # The x and y theory values.
        self.theoryx = 0.
        self.theory = 999.

        # The derived limit on the y-parameter, which is a function of x.
        self.limit = 0.

        self.loglike = -1e101
        self.arg = self.args

        # Import the data from file. Do this here so
        # we only do it once.
        # The data should be x, y.
        self.data = NP.loadtxt(self.file)

    def SetLogLike(self):
        """ Calcualte lower bound with interpolation
        function, then apply an error function.
        Arguments:

        Returns: The loglike from this constraint.

        """
        # Interpolate the y-value of the limit
        # corresponding to the theory x-value.
        self.limit = NP.interp(self.theoryx, self.data[0],
                               self.data[1], left=None, right=None)

        if self.tau:

          # Now calcualte likelihood with Gaussian error function - erf.
          like = 0.5 + 0.5 * \
              scipy.special.erf((self.theory - self.limit) / (2. ** 0.5 * self.tau * self.theory))
          try:
              self.loglike = math.log(like)
          # Sometimes we are trying to take log zero.
          except ValueError:
              self.loglike = -1e101

        else:

          self.loglike = (self.theory > self.limit) * -1e101

        return self.loglike


@Cube.memoize
class LikeMapConstraint:

    """ Interpolate likelihood data file. """

    def __init__(self, file, apply=True):
        """ Interpolates a likelihood from a data file
        Arguments:
        file -- Name of the data file.
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.file = os.path.abspath(file)
        self.apply = apply
        # The x parameter.
        self.theory = 0.
        # The y parameter.
        self.theoryy = 0.
        self.loglike = -1e101
        self.arg = self.args
        # Load file now so only loaded once.
        self.data = NP.loadtxt(self.file)

    def SetLogLike(self):
        """ Calcualte lower bound with interpolation
        function, then apply an error function.
        Arguments:

        Returns: The loglike from this constraint.

        """
        try:
            # Interpolate the likelihood from the datafile.
            self.loglike = math.log(mlab.griddata(
                self.data[0], self.data[1], self.data[2],
                NP.array([self.theory]), NP.array([self.theoryy]),
                interp='nn'))

        # Sometimes we are trying to take log zero.
        except ValueError:
            self.loglike = -1e101

        return self.loglike


@Cube.memoize
class ExternalConstraint:

    """ External contribution to likelihood - the likelihood is read from an external
    program.
    """

    def __init__(self, apply=True):
        """ Dummy constraint with no contribution to the likelhoood.
        Arguments:
        apply -- Whether to apply this constraint.

        Returns:

        """
        self.apply = apply

        self.theory = 0.
        self.loglike = -1e101
        self.arg = self.args

    def SetLogLike(self):
        """ External constraint - nothing to compute.
        Arguments:

        Returns: The loglike from this constraint.

        """
        return self.loglike

#########################################################################

# These are functions for manipulating programs and datafiles.


def RunProgram(executable, path, arguments, shell=False):
    """ Call a program
    and save its output to a file.

    Arguments:
    executable -- Name of executable file.
    path -- Path to executable.
    arguments -- Command-line arguments for executable.
    shell -- Whether to run through a shell, sometimes useufl.

    Returns:
    Name of file containing program output.

    """

    command = [executable, arguments]
    outputfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=False)
    errorfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=True)

    try:
        # Run program.
        program = subprocess.Popen(
            command,
            stdout=outputfile.fileno(),
            stderr=errorfile.fileno(), shell=shell,
            cwd=os.path.abspath(path))
        # Wait until it has finished!
        program.wait()

        # Return filename.
        print ("Output file name:", outputfile.name)
        return outputfile.name
    except:
        import sys
        print ('Error running program.', sys.exc_info()[0])
        os.remove(outputfile.name)
        return None


    
def RunProgramSPheno(executable, path, arguments, shell=False):
    """ Call Spheno
    and save its output to a file.

    Arguments:
    executable -- Name of executable file.
    path -- Path to executable.
    arguments -- Command-line arguments for executable.
    shell -- Whether to run through a shell, sometimes useufl.

    Returns:
    Name of file containing program output.

    """
    
    cwdir = os.getcwd()
    pardir = os.path.dirname(cwdir)

    # create HT input file rootname
    HTfilesroot = arguments[1] 
    #print('arguments SPHENO = ', arguments)
    #print('HTfilesroot arguments = ', HTfilesroot, arguments)

    command = [executable, arguments[0], HTfilesroot, pardir]
    outputfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=False)
    errorfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=False)
    #print('command position = ', command)


    if ShowSPhenoProcess == 'False':
        processfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=True)
   
    try:
        # Run program.	# stdout=None, because it means what SPheno writes in the screen, not the program's output
	    program = subprocess.Popen(
            command,
            cwd=os.path.abspath(path),
            stdout=processfile.fileno() if ShowSPhenoProcess == 'False' else None, 
            #stdout=None,
            stderr=errorfile.fileno(), shell=shell)
	    #print('command position 2 os.path.abspath(path) = ', os.path.abspath(path))
	    # Wait until it has finished!
	    program.wait()
	    #print('command position 3 errorfile outputfile ', errorfile.name, outputfile.name)
		
        # Copy SPheno output to the temporary file folder to use it as input for the other external programs
	    dir_SPheno_output = '../SPhenomunuSSM/output/'+HTfilesroot+'SPheno.spc.munuSSM3G'
	    temp_out_dir = outputfile.name
	    #print('command position 4 0  ... dir_SPheno_output, temp_out_dir = ', dir_SPheno_output, temp_out_dir)

	    #print('command position 4')
	    shutil.move(dir_SPheno_output, temp_out_dir)
	    #print('command position 5 ... dir_SPheno_output, temp_out_dir = ', dir_SPheno_output, temp_out_dir)

        # TODO
        # remove ht input files from hbout3g directoory 

        # Return filename.
	    print ("Output file name:", outputfile.name)
	    return outputfile.name
    except:
        import sys
        print ('Error running program.', sys.exc_info()[0])
        os.remove(outputfile.name)
        return None


def RunProgramHiggsTools(executable, path, arguments, shell=False):
    """ Call HiggsTools
    and save its output to a file.

    Arguments:
    executable -- Name of executable file.
    path -- Path to executable.
    arguments -- Command-line arguments for executable.
    shell -- Whether to run through a shell, sometimes useufl.

    Returns:
    Name of file containing program output.

    """
    
    
    
    HTfilesroot = arguments[1] 
    #print('HTfilesroot HIGGSTOOLS = ', HTfilesroot)

    # 15 neutral CP even and CP odd scalars, and 7 charged scalars
    htinname = 'HBout3G/'+HTfilesroot+'-'
    command = [executable, arguments[0], htinname, '15', '7',
     '../data/hbdataset-master', 
     '../data/hsdataset-main',
     'True', 'True', HTfilesroot]

    outputfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=False)
    errorfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=True)
    #print('errorfile HT = ', errorfile.name)
   
    '''
    # Run program. # stdout=None, because it means what HT writes in the screen, not the program's output
    program = subprocess.Popen(
             command,
             cwd=os.path.abspath(path),
             stdout=outputfile.fileno() if ShowSPhenoProcess == 'False' else None, 
             #stdout=None,
             stderr=errorfile.fileno(), shell=shell)
    print('command position HT os.path.abspath(path) = ', os.path.abspath(path))
    # Wait until it has finished!
    program.wait()
    # print('command position 3')
    return outputfile.name
    '''

    try:
        # Run program. # stdout=None, because it means what HT writes in the screen, not the program's output
	    program = subprocess.Popen(
            command,
            cwd=os.path.abspath(path),
            stdout=outputfile.fileno() if ShowSPhenoProcess == 'False' else None,  # ADDED option in InputFile
            # stdout=None,
            stderr=errorfile.fileno(), shell=shell)
	    #print('command position HT os.path.abspath(path) = ', os.path.abspath(path))
	    # Wait until it has finished!
	    program.wait()
	    #print('HT command position 3')
		
        # Copy HT output to the temporary file folder to use it as input for the other external programs
	    dir_HT_output = '../higgstools-main/buildpy/'+htinname+'HiggsTools_results.tsv' 
	    # os.path.basename(inputFilepath) 
	    temp_out_dir = outputfile.name
	    #print('command position 4 HT', temp_out_dir, dir_HT_output)
	    shutil.move(dir_HT_output, temp_out_dir)
	    #print('command position 5 HT')

        # Return filename.
	    #print('HT outputfile.name = ', outputfile.name)
	    #print ("Output file name:", outputfile.name, os.path.abspath(path)+'/HBout3G/', HTfilesroot)
        
        # TODO
        # remove ht input files from hbout3g directoory 
	    try:
	        del_file_patterns(os.path.abspath(path)+'/HBout3G/', HTfilesroot)
	    except OSError:
	        #pass
	        print ("HT files inside HBout3G/ not deleted")

	    return outputfile.name
    except:
        import sys
        print ('Error running program.', sys.exc_info()[0])
        os.remove(outputfile.name)
        return None



def RunProgramDDCalc(executable, path, shell=True):
    """ Call a program
    and save its output to a file.

    Arguments:
    executable -- Name of executable file.
    path -- Path to executable.
    arguments -- Command-line arguments for executable.
    shell -- Whether to run through a shell, sometimes useufl.

    Returns:
    Name of file containing program output.

    """

    outputfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=False)
    command = executable
    errorfile = tempfile.NamedTemporaryFile(dir='../pyMUNUSSM/temFiles/', delete=True)

    try:
        # Run program.
        program = subprocess.Popen(
            command,
            stdout=outputfile.fileno(),
            stderr=errorfile.fileno(), shell=shell,
            cwd=os.path.abspath(path))
        # Wait until it has finished!
        program.wait()

        # Return filename.
        #print ("Output file name:", outputfile.name)
        return outputfile.name
    except:
        import sys
        print ('Error running program.', sys.exc_info()[0])
        os.remove(outputfile.name)
        return None

def del_file_with_patterns(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def del_file_patterns(dir, pattern):
	import glob
	# Getting All Files List
	fileList = glob.glob(dir+'/'+pattern+'*', recursive=True)
	# Remove all files one by one
	for file in fileList:
		try:
		    os.remove(file)
		except OSError:
		    print("Error while deleting file")
