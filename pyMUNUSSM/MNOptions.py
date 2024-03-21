#########################################################################
#                                                                       #
#    M u l t i n e s t     O p t i o n s                                #
#                                                                       #
#########################################################################

# External modules.
import os

# Internal modules.
import Priors  # For setting callback function for setting number of paramters.
import Likelihood  # For setting callback function and number of constraints.

# Make chains directory for MultiNest.
if not os.path.exists("chains"):
    os.mkdir("chains")

# File names - best to convert to absolute path here.
outputfiles_basename = os.path.abspath('./chains/munuSSM')		

# These are the settings that are passed to MultiNest.
# Loglike callback function.
LogLikelihood = Likelihood.myloglike
# Prior callback function.
Prior = Priors.myprior
# Number of model parameters.
n_dims = len(Priors.MUNUSSMModelTracker().param)	
# Total number of parameters in cube.
# This is model parameters + constraints + chi2 from constraints (2) + sparticle/particle masses (58)
# sparticle/particle masses 
n_params = n_dims + 2 * \
    len(Likelihood.MUNUSSMConstraintTracker().constraint) + 2 + 60 		
# Which parameters to mode cluster, first n.
n_clustering_params = 2
# Should be a vector of 1s and 0s - 1 = wrapped, 0 = unwrapped.
wrapped_params = None
# Whether to separate modes.
multimodal = True
const_efficiency_mode = False
n_live_points = 1000
evidence_tolerance = 1.0
sampling_efficiency = 2.0
n_iter_before_update = 1
null_log_evidence = -1e+90
# For memory allocation only, if memory exceeded, MultiNest halts.
max_modes = 5
# Seed the random number generator, if -1 seeded from system clock.
seed = -1
# Whether to write information to screen.
verbose = True
# Resume from unfinished run.
resume = True
# Extra unused parameter, could be used to pass info to MultiNest.
context = 0
# Importance sampling.
IS = True
