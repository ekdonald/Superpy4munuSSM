#########################################################################
#                                                                       #
#    P r i o r s                                                        #
#    Priors for the NMUNUSSM are set here.                          #
#                                                                       #
#########################################################################

# External modules.
import math
import scipy
from scipy import special

# Superpy modules.
import numpy as NP
import Cube


################################################################################################

# Read the InputFile to see the parameters and parameters options (aka diagonal or not)

def ReadInputFileParam(filename, start):								#ADDED
        """ Read from InputFileParam to set parameters options 

        Arguments:
        filename -- File to be read from.
        start -- String of beginning of line to read.


        Returns:
        Parameter retrieved from the file.

        """

        aux = 0

        for line in open(filename, 'r'):
            if line.lstrip().startswith(start):
                aux = line.split(start)
                parameter = aux[-1][:-1]
                return parameter

        if aux == 0:
            print ("ERROR: Check InputFile, there is no:", start_)
            sys.exit()

# PARAMETER SET UP (diagonal matrix or not?)
InputParamFormat = ReadInputFileParam('InputFileParam','InputParamFormat = ')

# SHOW IN TERMINAL OPTIONS
ShowParameters = ReadInputFileParam('InputFile','ShowParameters = ')


################################################################################################

def myprior(cube, ndim, nparams):
    """ MultiNest callback function. Sets the model parameters
    from the unit hypercube.
    Arguments:
    cube -- Unit hypercube from which model parameters are mapped.
    ndim -- Number of model paramers.
    nparams -- Total number of parameters in the cube .

    Returns:

    """
    # Set-up the paramters - MUNUSSM model.					#MODIFIED!
    Model = MUNUSSMModelTracker()

    # Set the model parameters from the cube.
    Model.SetParams(cube)

    # Copy them to the cube so that they are printed, in alphabetical order.
    for name in sorted(Model.param.keys(), key=str.lower):
        Cube.AddCube(cube, Model.param[name].value, Cube.label, name)

    # Print-out cube model parameters for debugging.
    print ('*****************************************************')
    if ShowParameters == 'True':                           # ADDED option in InputFile
        print ('Parameters:')
        for i in range(Cube.AddCube.count() + 1):  # Count begins at 0 - need +1.
            print (Cube.label[i], cube[i])


#########################################################################

# Model class - this is where everything is stored.
# This class contains all the parameters for the model,
# and the functions to evalulate them.


#------ BY DONALD --- STOPPED HERE AT 03:AM ON FEB 13, 2024

class MUNUSSMModelTracker:								# TO BE MODIFIED ON FEB 13, 2024

    """ Contains NMUNUSSM model parameters, priors and functions to
    evaluate the parameters. """

    def __init__(self):
        """ Sets-up the model.
        Specify the priors for the model and nuisance parameters.
        """

        self.param = {}  # Holds model's parameters.

        # First, create a file in the HT folder to tell it the number of free model parameters

        Nparam = ReadInputFileParam('InputFileParam','Nparam'+InputParamFormat+' = ')

        NparamAux = open('../higgstools-main/buildpy/Nparam.dat', 'w+') 
        NparamAux.write('%s' % Nparam)
        NparamAux.close()

        # Set the priors for each parameter. NB Gaussian priors are in a way
        # infinite, so you don't need to specify upper and lower bounds.
        # NB If you use a log parameter, be sure its always strictly > 0.
        # The parameters are in GeV.


	    # Renormalization scale		# TP1MOD
        self.param['renormscale'] = eval(ReadInputFileParam('InputFileParam','renormscale'+InputParamFormat+' = '))

        # Soft masses for gauginos.
        self.param['M1'] = eval(ReadInputFileParam('InputFileParam','M1'+InputParamFormat+' = '))
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['M2'] = eval(ReadInputFileParam('InputFileParam','M2'+InputParamFormat+' = '))
            self.param['M3'] = eval(ReadInputFileParam('InputFileParam','M3'+InputParamFormat+' = '))


        #Ratio of Higgs VEVs, defined at EW scale.
        self.param['tanbeta'] = eval(ReadInputFileParam('InputFileParam','tanbeta'+InputParamFormat+' = '))

        # right vevs.
        self.param['vR1'] = eval(ReadInputFileParam('InputFileParam','vR1'+InputParamFormat+' = '))
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['vR2'] = self.param['vR1']
            self.param['vR3'] = self.param['vR1'] 
        else:
            self.param['vR2'] = eval(ReadInputFileParam('InputFileParam','vR2'+InputParamFormat+' = '))
            self.param['vR3'] = eval(ReadInputFileParam('InputFileParam','vR3'+InputParamFormat+' = '))

        # left vevs.
        self.param['vL1'] = eval(ReadInputFileParam('InputFileParam','vL1'+InputParamFormat+' = '))
        self.param['vL2'] = eval(ReadInputFileParam('InputFileParam','vL2'+InputParamFormat+' = '))
        self.param['vL3'] = eval(ReadInputFileParam('InputFileParam','vL3'+InputParamFormat+' = '))
	
	    #Couplings in the superpotential and soft-breaking counter-parts

        #lambda_i vR Hu Hd in MUNUSSM superpotential.
        self.param['lam1'] = eval(ReadInputFileParam('InputFileParam','lambda1'+InputParamFormat+' = '))
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['lam2'] = self.param['lam1']
            self.param['lam3'] = self.param['lam1']  
        else:
            self.param['lam2'] = eval(ReadInputFileParam('InputFileParam','lambda2'+InputParamFormat+' = '))
            self.param['lam3'] = eval(ReadInputFileParam('InputFileParam','lambda3'+InputParamFormat+' = '))


        # Tlambda_i vR Hu Hd MUNUSSM soft-breaking term (Tlambda = Alambda * lambda).
        self.param['Tlam1'] = eval(ReadInputFileParam('InputFileParam','Tlambda1'+InputParamFormat+' = ')) 
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Tlam2'] = self.param['Tlam1']
            self.param['Tlam3'] = self.param['Tlam1']
        else:
            self.param['Tlam2'] = eval(ReadInputFileParam('InputFileParam','Tlambda2'+InputParamFormat+' = '))
            self.param['Tlam3'] = eval(ReadInputFileParam('InputFileParam','Tlambda3'+InputParamFormat+' = ')) 


        # kappa_ijk vR vR vR in MUNUSSM superpotential. Only diagonal terms are non-zero
        self.param['kap111'] = eval(ReadInputFileParam('InputFileParam','kappa1'+InputParamFormat+' = '))
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['kap222'] = self.param['kap111']
            self.param['kap333'] = self.param['kap111']
        else:
            #self.param['kap111'] = eval(ReadInputFileParam('InputFileParam','kappa1'+InputParamFormat+' = '))
            self.param['kap222'] = eval(ReadInputFileParam('InputFileParam','kappa2'+InputParamFormat+' = '))
            self.param['kap333'] = eval(ReadInputFileParam('InputFileParam','kappa3'+InputParamFormat+' = '))


        # Tkappa_ijk vR vR vR in MUNUSSM soft-breaking term (Tkappa=Akappa*kappa). Off-diagonal are zero
        self.param['Tk111'] = eval(ReadInputFileParam('InputFileParam','Tkappa1'+InputParamFormat+' = '))
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Tk222'] = self.param['Tk111'] 
            self.param['Tk333'] = self.param['Tk111'] 
        else:
            self.param['Tk222'] = eval(ReadInputFileParam('InputFileParam','Tkappa2'+InputParamFormat+' = '))
            self.param['Tk333'] = eval(ReadInputFileParam('InputFileParam','Tkappa3'+InputParamFormat+' = '))


	    #Yukawa matrix, neutrinos

        # Y_v		Y_v L H_u N term in the superpotential
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Yv11'] = eval(ReadInputFileParam('InputFileParam','Yv11'+InputParamFormat+' = '))
            self.param['Yv22'] = eval(ReadInputFileParam('InputFileParam','Yv22'+InputParamFormat+' = '))
            self.param['Yv33'] = eval(ReadInputFileParam('InputFileParam','Yv33'+InputParamFormat+' = '))
        else:
            self.param['Yv11'] = eval(ReadInputFileParam('InputFileParam','Yv11'+InputParamFormat+' = '))
            self.param['Yv12'] = eval(ReadInputFileParam('InputFileParam','Yv12'+InputParamFormat+' = '))
            self.param['Yv13'] = eval(ReadInputFileParam('InputFileParam','Yv13'+InputParamFormat+' = '))
            self.param['Yv21'] = eval(ReadInputFileParam('InputFileParam','Yv21'+InputParamFormat+' = '))
            self.param['Yv22'] = eval(ReadInputFileParam('InputFileParam','Yv22'+InputParamFormat+' = '))
            self.param['Yv23'] = eval(ReadInputFileParam('InputFileParam','Yv23'+InputParamFormat+' = '))
            self.param['Yv31'] = eval(ReadInputFileParam('InputFileParam','Yv31'+InputParamFormat+' = '))
            self.param['Yv32'] = eval(ReadInputFileParam('InputFileParam','Yv32'+InputParamFormat+' = '))
            self.param['Yv33'] = eval(ReadInputFileParam('InputFileParam','Yv33'+InputParamFormat+' = '))


	    #Soft-breaking masses for all the particles

        # m^2_d.
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['m2d11'] = eval(ReadInputFileParam('InputFileParam','m2d11'+InputParamFormat+' = '))
            self.param['m2d22'] = eval(ReadInputFileParam('InputFileParam','m2d22'+InputParamFormat+' = '))
            self.param['m2d33'] = eval(ReadInputFileParam('InputFileParam','m2d33'+InputParamFormat+' = '))
        else:
            self.param['m2d11'] = eval(ReadInputFileParam('InputFileParam','m2d11'+InputParamFormat+' = '))
            self.param['m2d22'] = eval(ReadInputFileParam('InputFileParam','m2d22'+InputParamFormat+' = '))
            self.param['m2d33'] = eval(ReadInputFileParam('InputFileParam','m2d33'+InputParamFormat+' = '))
            self.param['m2d12'] = eval(ReadInputFileParam('InputFileParam','m2d12'+InputParamFormat+' = '))
            self.param['m2d13'] = eval(ReadInputFileParam('InputFileParam','m2d13'+InputParamFormat+' = '))
            self.param['m2d21'] = eval(ReadInputFileParam('InputFileParam','m2d21'+InputParamFormat+' = '))
            self.param['m2d23'] = eval(ReadInputFileParam('InputFileParam','m2d23'+InputParamFormat+' = '))
            self.param['m2d31'] = eval(ReadInputFileParam('InputFileParam','m2d31'+InputParamFormat+' = '))
            self.param['m2d32'] = eval(ReadInputFileParam('InputFileParam','m2d32'+InputParamFormat+' = '))


        # m^2_u.
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['m2u22'] = eval(ReadInputFileParam('InputFileParam','m2u22'+InputParamFormat+' = '))
            self.param['m2u33'] = eval(ReadInputFileParam('InputFileParam','m2u33'+InputParamFormat+' = '))
            self.param['m2u11'] = eval(ReadInputFileParam('InputFileParam','m2u11'+InputParamFormat+' = '))
        else:
            self.param['m2u22'] = eval(ReadInputFileParam('InputFileParam','m2u22'+InputParamFormat+' = '))
            self.param['m2u33'] = eval(ReadInputFileParam('InputFileParam','m2u33'+InputParamFormat+' = '))
            self.param['m2u11'] = eval(ReadInputFileParam('InputFileParam','m2u11'+InputParamFormat+' = '))
            self.param['m2u12'] = eval(ReadInputFileParam('InputFileParam','m2u12'+InputParamFormat+' = '))
            self.param['m2u13'] = eval(ReadInputFileParam('InputFileParam','m2u13'+InputParamFormat+' = '))
            self.param['m2u21'] = eval(ReadInputFileParam('InputFileParam','m2u21'+InputParamFormat+' = '))
            self.param['m2u23'] = eval(ReadInputFileParam('InputFileParam','m2u23'+InputParamFormat+' = '))
            self.param['m2u31'] = eval(ReadInputFileParam('InputFileParam','m2u31'+InputParamFormat+' = '))
            self.param['m2u32'] = eval(ReadInputFileParam('InputFileParam','m2u32'+InputParamFormat+' = '))


        # m^2_Q.
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['m2Q11'] = eval(ReadInputFileParam('InputFileParam','m2Q11'+InputParamFormat+' = '))
            self.param['m2Q22'] = eval(ReadInputFileParam('InputFileParam','m2Q22'+InputParamFormat+' = '))
            self.param['m2Q33'] = eval(ReadInputFileParam('InputFileParam','m2Q33'+InputParamFormat+' = '))
        else:
            self.param['m2Q11'] = eval(ReadInputFileParam('InputFileParam','m2Q11'+InputParamFormat+' = '))
            self.param['m2Q22'] = eval(ReadInputFileParam('InputFileParam','m2Q22'+InputParamFormat+' = '))
            self.param['m2Q33'] = eval(ReadInputFileParam('InputFileParam','m2Q33'+InputParamFormat+' = '))
            self.param['m2Q12'] = eval(ReadInputFileParam('InputFileParam','m2Q12'+InputParamFormat+' = '))
            self.param['m2Q13'] = eval(ReadInputFileParam('InputFileParam','m2Q13'+InputParamFormat+' = '))
            self.param['m2Q21'] = eval(ReadInputFileParam('InputFileParam','m2Q21'+InputParamFormat+' = '))
            self.param['m2Q23'] = eval(ReadInputFileParam('InputFileParam','m2Q23'+InputParamFormat+' = '))
            self.param['m2Q31'] = eval(ReadInputFileParam('InputFileParam','m2Q31'+InputParamFormat+' = '))
            self.param['m2Q32'] = eval(ReadInputFileParam('InputFileParam','m2Q32'+InputParamFormat+' = '))

        # mlHd2. --- always set to zero.
        #self.param['mlHd21'] = eval(ReadInputFileParam('InputFileParam','mlHd21'+InputParamFormat+' = '))
        #self.param['mlHd22'] = eval(ReadInputFileParam('InputFileParam','mlHd22'+InputParamFormat+' = '))
        #self.param['mlHd23'] = eval(ReadInputFileParam('InputFileParam','mlHd23'+InputParamFormat+' = '))


        # m^2_e.
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['m2e11'] = eval(ReadInputFileParam('InputFileParam','m2e11'+InputParamFormat+' = '))
            self.param['m2e22'] = eval(ReadInputFileParam('InputFileParam','m2e22'+InputParamFormat+' = '))
            self.param['m2e33'] = eval(ReadInputFileParam('InputFileParam','m2e33'+InputParamFormat+' = '))
        else:
            self.param['m2e11'] = eval(ReadInputFileParam('InputFileParam','m2e11'+InputParamFormat+' = '))
            self.param['m2e22'] = eval(ReadInputFileParam('InputFileParam','m2e22'+InputParamFormat+' = '))
            self.param['m2e33'] = eval(ReadInputFileParam('InputFileParam','m2e33'+InputParamFormat+' = '))
            self.param['m2e12'] = eval(ReadInputFileParam('InputFileParam','m2e12'+InputParamFormat+' = '))
            self.param['m2e13'] = eval(ReadInputFileParam('InputFileParam','m2e13'+InputParamFormat+' = '))
            self.param['m2e21'] = eval(ReadInputFileParam('InputFileParam','m2e21'+InputParamFormat+' = '))
            self.param['m2e23'] = eval(ReadInputFileParam('InputFileParam','m2e23'+InputParamFormat+' = '))
            self.param['m2e31'] = eval(ReadInputFileParam('InputFileParam','m2e31'+InputParamFormat+' = '))
            self.param['m2e32'] = eval(ReadInputFileParam('InputFileParam','m2e32'+InputParamFormat+' = '))


        ## m^2_v. solved from tadpoles --- no need to provide them
        # self.param['m2v11'] = eval(ReadInputFileParam('InputFileParam','m2v11'+InputParamFormat+' = '))
        # self.param['m2v22'] = eval(ReadInputFileParam('InputFileParam','m2v22'+InputParamFormat+' = '))
        # self.param['m2v33'] = eval(ReadInputFileParam('InputFileParam','m2v33'+InputParamFormat+' = '))
        # self.param['m2v12'] = eval(ReadInputFileParam('InputFileParam','m2v12'+InputParamFormat+' = '))
        # self.param['m2v13'] = eval(ReadInputFileParam('InputFileParam','m2v13'+InputParamFormat+' = '))
        # self.param['m2v21'] = eval(ReadInputFileParam('InputFileParam','m2v21'+InputParamFormat+' = '))
        # self.param['m2v23'] = eval(ReadInputFileParam('InputFileParam','m2v23'+InputParamFormat+' = '))
        # self.param['m2v31'] = eval(ReadInputFileParam('InputFileParam','m2v31'+InputParamFormat+' = '))
        # self.param['m2v32'] = eval(ReadInputFileParam('InputFileParam','m2v32'+InputParamFormat+' = '))


        ## m^2_L. solved from tadpoles --- no need to provide them
        # self.param['m2L11'] = eval(ReadInputFileParam('InputFileParam','m2L11'+InputParamFormat+' = '))
        # self.param['m2L22'] = eval(ReadInputFileParam('InputFileParam','m2L22'+InputParamFormat+' = '))
        # self.param['m2L33'] = eval(ReadInputFileParam('InputFileParam','m2L33'+InputParamFormat+' = '))
        # self.param['m2L12'] = eval(ReadInputFileParam('InputFileParam','m2L12'+InputParamFormat+' = '))
        # self.param['m2L13'] = eval(ReadInputFileParam('InputFileParam','m2L13'+InputParamFormat+' = '))
        # self.param['m2L21'] = eval(ReadInputFileParam('InputFileParam','m2L21'+InputParamFormat+' = '))
        # self.param['m2L23'] = eval(ReadInputFileParam('InputFileParam','m2L23'+InputParamFormat+' = '))
        # self.param['m2L31'] = eval(ReadInputFileParam('InputFileParam','m2L31'+InputParamFormat+' = '))
        # self.param['m2L32'] = eval(ReadInputFileParam('InputFileParam','m2L32'+InputParamFormat+' = '))


	    #Other Soft-breaking terms

        # Td_ij  # Td_ij H_d Q d in MUNUSSM soft-breaking term (Td = Ad * Yd).
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Td33'] = eval(ReadInputFileParam('InputFileParam','Td33'+InputParamFormat+' = '))
        else:
            self.param['Td11'] = eval(ReadInputFileParam('InputFileParam','Td11'+InputParamFormat+' = '))
            self.param['Td22'] = eval(ReadInputFileParam('InputFileParam','Td22'+InputParamFormat+' = '))
            self.param['Td12'] = eval(ReadInputFileParam('InputFileParam','Td12'+InputParamFormat+' = '))
            self.param['Td13'] = eval(ReadInputFileParam('InputFileParam','Td13'+InputParamFormat+' = '))
            self.param['Td21'] = eval(ReadInputFileParam('InputFileParam','Td21'+InputParamFormat+' = '))
            self.param['Td23'] = eval(ReadInputFileParam('InputFileParam','Td23'+InputParamFormat+' = '))
            self.param['Td31'] = eval(ReadInputFileParam('InputFileParam','Td31'+InputParamFormat+' = '))
            self.param['Td32'] = eval(ReadInputFileParam('InputFileParam','Td32'+InputParamFormat+' = '))
            self.param['Td33'] = eval(ReadInputFileParam('InputFileParam','Td33'+InputParamFormat+' = '))


        # Tu_ij  # Tu_ij H_u Q u in NMUNUSSM soft-breaking term (Tu = Au * Yu).
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Tu33'] = eval(ReadInputFileParam('InputFileParam','Tu33'+InputParamFormat+' = '))
        else:
            self.param['Tu11'] = eval(ReadInputFileParam('InputFileParam','Tu11'+InputParamFormat+' = '))
            self.param['Tu22'] = eval(ReadInputFileParam('InputFileParam','Tu22'+InputParamFormat+' = '))
            self.param['Tu12'] = eval(ReadInputFileParam('InputFileParam','Tu12'+InputParamFormat+' = '))
            self.param['Tu13'] = eval(ReadInputFileParam('InputFileParam','Tu13'+InputParamFormat+' = '))
            self.param['Tu21'] = eval(ReadInputFileParam('InputFileParam','Tu21'+InputParamFormat+' = '))
            self.param['Tu23'] = eval(ReadInputFileParam('InputFileParam','Tu23'+InputParamFormat+' = '))
            self.param['Tu31'] = eval(ReadInputFileParam('InputFileParam','Tu31'+InputParamFormat+' = '))
            self.param['Tu32'] = eval(ReadInputFileParam('InputFileParam','Tu32'+InputParamFormat+' = '))
            self.param['Tu33'] = eval(ReadInputFileParam('InputFileParam','Tu33'+InputParamFormat+' = '))


        # Te_ij  # Te_ij H_d L e in NMUNUSSM soft-breaking term (Te = Ae * Ye).
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Te33'] = eval(ReadInputFileParam('InputFileParam','Te33'+InputParamFormat+' = '))
        else:
            self.param['Te11'] = eval(ReadInputFileParam('InputFileParam','Te11'+InputParamFormat+' = '))
            self.param['Te22'] = eval(ReadInputFileParam('InputFileParam','Te22'+InputParamFormat+' = '))
            self.param['Te12'] = eval(ReadInputFileParam('InputFileParam','Te12'+InputParamFormat+' = '))
            self.param['Te13'] = eval(ReadInputFileParam('InputFileParam','Te13'+InputParamFormat+' = '))
            self.param['Te21'] = eval(ReadInputFileParam('InputFileParam','Te21'+InputParamFormat+' = '))
            self.param['Te23'] = eval(ReadInputFileParam('InputFileParam','Te23'+InputParamFormat+' = '))
            self.param['Te31'] = eval(ReadInputFileParam('InputFileParam','Te31'+InputParamFormat+' = '))
            self.param['Te32'] = eval(ReadInputFileParam('InputFileParam','Te32'+InputParamFormat+' = '))
            self.param['Te33'] = eval(ReadInputFileParam('InputFileParam','Te33'+InputParamFormat+' = '))


        # Tv_ij  # Tv_ij H_u L vR in NMUNUSSM soft-breaking term (Tv = Av * Yv).
        self.param['Tv11'] = eval(ReadInputFileParam('InputFileParam','Tv11'+InputParamFormat+' = '))
        if InputParamFormat == 'A' or InputParamFormat == 'B':
            self.param['Tv22'] = eval(ReadInputFileParam('InputFileParam','Tv22'+InputParamFormat+' = '))
            self.param['Tv33'] = eval(ReadInputFileParam('InputFileParam','Tv33'+InputParamFormat+' = '))
        else:
            self.param['Tv11'] = eval(ReadInputFileParam('InputFileParam','Tv11'+InputParamFormat+' = '))
            self.param['Tv22'] = eval(ReadInputFileParam('InputFileParam','Tv22'+InputParamFormat+' = '))
            self.param['Tv33'] = eval(ReadInputFileParam('InputFileParam','Tv33'+InputParamFormat+' = '))
            self.param['Tv12'] = eval(ReadInputFileParam('InputFileParam','Tv12'+InputParamFormat+' = '))
            self.param['Tv13'] = eval(ReadInputFileParam('InputFileParam','Tv13'+InputParamFormat+' = '))
            self.param['Tv21'] = eval(ReadInputFileParam('InputFileParam','Tv21'+InputParamFormat+' = '))
            self.param['Tv23'] = eval(ReadInputFileParam('InputFileParam','Tv23'+InputParamFormat+' = '))
            self.param['Tv31'] = eval(ReadInputFileParam('InputFileParam','Tv31'+InputParamFormat+' = '))
            self.param['Tv32'] = eval(ReadInputFileParam('InputFileParam','Tv32'+InputParamFormat+' = '))


        # Nusisance parameters.

        # Bottom MS-bar mass.
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/DataBlock.action?node=Q005M
        self.param['mb'] = eval(ReadInputFileParam('InputFileParam','mb'+InputParamFormat+' = '))

        # Top pole mass.
        # PDG.
        # http://pdg8.lbl.gov/rpp2014v1/pdgLive/DataBlock.action?node=Q007T7
        self.param['mt'] = eval(ReadInputFileParam('InputFileParam','mt'+InputParamFormat+' = '))

        # Strong coupling.
        # PDG
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-qcd.pdf
        self.param['alphas'] = eval(ReadInputFileParam('InputFileParam','alphas'+InputParamFormat+' = '))

        # Reciprocal of EM coupling at MZ.
        # PDG.
        # http://pdg.lbl.gov/2014/reviews/rpp2014-rev-standard-model.pd
        self.param['invalpha'] = eval(ReadInputFileParam('InputFileParam','invalpha'+InputParamFormat+' = '))


    def SetParams(self, x):
        """ Sets the parameters from the unit hypercube.
        Arguments:
        x -- The whole unit hypercube vector.

        Returns:

        """
        for i, name in enumerate(self.param.keys()):
            self.param[name].SetValue(x[i])

#########################################################################

# These classes are all the types of prior one might use, they all have a
# SetValue function, which sets the value of parameter from MultiNest's unit
# hypercube. i.e. all priors are mapped from a Uniform(0,1) distribution.


@Cube.memoize
class LogParameter:

    """ A parameter with a log prior. """

    def __init__(self, min, max):
        """ Initializes a parameter with a log prior.
        Arguments:
        min -- Minimum of prior.
        max -- Maximum of prior.

        Returns:

        """
        self.minimum = min
        self.maximum = max
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Calculates value of parameter with log prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Log the ranges, map from 0-1 to min-max in log space,
        # then exponentiate.
        self.value = math.pow(10, math.log10(self.minimum) + x *
                              (math.log10(self.maximum) -
                               math.log10(self.minimum)))


@Cube.memoize
class NegLogParameter:                                  # ADDED

    """ A parameter with a log prior, for negative values. """

    def __init__(self, min, max):
        """ Initializes a parameter with a log prior.
        Arguments:
        min -- Minimum of prior.
        max -- Maximum of prior.

        Returns:

        """
        self.minimum = -max
        self.maximum = -min
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Calculates value of parameter with log prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Log the ranges, map from 0-1 to min-max in log space,
        # then exponentiate.
        self.value = -math.pow(10, math.log10(self.minimum) + x *
                              (math.log10(self.maximum) -
                               math.log10(self.minimum)))


@Cube.memoize
class LinearParameter:

    """ A parameter with a linear prior. """

    def __init__(self, min, max):
        """ Initializes a parameter with a linear prior .
        Arguments:
        min -- Minimum of prior.
        max -- Maximum of prior.

        Returns:

        """
        self.minimum = min
        self.maximum = max
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Calculates value of parameter with a linear prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Map from 0-1 to min-max.
        self.value = self.minimum + x * (self.maximum - self.minimum)


@Cube.memoize
class FixedParameter:

    """ A parameter with a fixed value. """

    def __init__(self, fixed):
        """ Initializes a parameter with a fixed value.
        You might, for example, want to fix sgn mu in the CMSSM.
        Arguments:
        fixed -- Fixed value of parameters.

        Returns:

        """
        self.fix = fixed
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Sets value of parameter with a fixed prior, ignore unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube, ignored.

        Returns:

        """
        # NB that although x isn't used, it should stil be here
        # so that this class can be used like the others.
        self.value = self.fix
        self.arg = self.args


@Cube.memoize
class SignParameter:

    """ A parameter that is either +/- 1. """

    def __init__(self):
        """ Initializes a parameter that is either +/-1.
        For example, sign mu in the CMSSM.
        Arguments:

        Returns:

        """
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Sets value of parameter to +/- 1, with equal probability.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # You could change this to favour a particular sign.
        if x > 0.5:
            self.value = +1
        else:
            self.value = -1


@Cube.memoize
class GaussParameter:

    """ A parameter with a Gaussian prior """

    def __init__(self, mu, sigma):
        """ Initializes a parameter with a Gaussian prior.
        Arguments:
        mu -- Mean of Gaussian.
        sigma -- Variance of Gaussian.

        Returns:

        """
        self.mu = mu
        self.sigma = sigma
        self.value = 0
        self.arg = self.args

    def SetValue(self, x):
        """ Sets value of parameter with a Gaussian prior from unit
        hypercube.
        Arguments:
        x -- One element of unit hypercube.

        Returns:

        """
        # Map from 0-1 to a Gaussian! erfcinv is a the inverse of the
        # complementary Gaussian error function.
        self.value = self.mu + self.sigma * math.sqrt(2) *\
            scipy.special.erfcinv(2 * (1 - x))
