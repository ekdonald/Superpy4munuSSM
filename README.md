# ABOUT THIS SUPERPY:

    S u p e r P y

    Finds regions of SUSY models' parameter spaces that are in
    best agreement with experimental data.
    Project: SuperPy.
    Author: Andrew Fowlie, KBFI, Tallinn.
    Version: 1.1
    Date: 05/14
    https://github.com/andrewfowlie/superpy


   # ADAPTATION TO MUNUSSM  
    Author:  D. E. Kpatcha 
    Date: 2024
    Updated: subprograms
    New: more subprograms (SPhenomunussm is the main one), tools, model, easier and better 
		 input formats, etc

    The modifications w.r.t the original SuperPy are listed in MUNUSSM/MODIFICATIONS


# QUICK STARTUP: INSTALLATION

    1. The first step is to clean 
		a) use "make clean"  to clean DDCalc_v, MultiNest-master, and SPhenomunussm
          (warning ... higgstools not cleanned when using make clean)
		b) or use "make cleanall"  to clean all the subprograms (if installing for the first time)
        c) you may want to clean higgstools only use 
           "make clean_higgstools" for this purpose       



    2. Install now
		"make" or "make all"   # maybe you will need to adjust something (hopefully not)
                               # you need python 2.7 or later prob, and libraries, etc
                               # tested with python 3.10


    2. Inside ../superpy-munussm-v*/pyMUNUSSM/

    2.a) in the terminal:
        export LD_LIBRARY_PATH=/your/path/here/superpy-master-new/MultiNest-master/lib:$LD_LIBRARY_PATH
        in every new shell you want to run the program, and run 
        "python SuperPy.py"

    2.b) add once
        export LD_LIBRARY_PATH=/your/path/here/superpy-master-new/MultiNest-master/lib:$LD_LIBRARY_PATH
        at the end of your .bashrc, 
		and run "python SuperPy.py"
		
    3.)  BEFORE RUNNING YOU NEED TO PROVIDE DATA TO HIGGSTOOLS
   
        a) download the latest HiggsBounds dataset from: https://gitlab.com/higgsbounds/hbdataset

        b) download the latest HiggsSignals dataset from: https://gitlab.com/higgsbounds/hsdataset
     
    Then (unzip if necessary) put the two datasets ("hbdataset" and "hsdataset") inside
    
    "./higgstools-main/data/" folder



# QUICK STARTUP: RUNNING THE CODE

    1.) run 
        "./run.sh"
        mini script that contains the "export ..." command and "python SuperPy.py"


    And it will start!

    2. Output in pyMUNUSSM/chains


# QUICK INPUTS:


    There are 2 InputFiles into the folder "pyMUNUSSM"

    1. InputFile. 
        Set True or False if you want to use a particular subprogram, or constraint, or to show some 
        process in the shell, or to keep/delete auxiliary files, etc
        More details explained in that file.

    2. InputFileParam.
        Set the parameter values (with linear, gaussian, log, neglog, etc distributions) for this model 
        SPhenomunuSSM.
        Set the amount of parameters (GUT relation or not, diagonal matrix or not, etc)
        More details explained in that file.

    3. The amount of points for Multinest set in "MNOptions.py" 
        n_live_points

	
# CURRENT IMPLEMENTED SUBPROGRAMS CONSTRAINTS:

    0. spheno for munussm

    1. higgstools (https://gitlab.com/higgsbounds/higgstools)
       software for comparing a wide class of BSM models to all available experimental results 
       from searches for new (scalar) particles and measurements of the 125 GeV Higgs boson at 
       colliders. We use it to compare the munussm predictions to the
       experimental observations.

    2. neutrino physics (See Likelihood.py file) 
       squared masses difference, mixing angles, and cosmological constraint on the sum of the masses
    
    3. spheno observables: for example (can be extended to more observables ... feel free to 
       chose/add what you need. See Likelihood.py file)
       g-2 of muon  
       b -> s gamma 
       B -> mu nu 
       B -> tau nu 
       Bd -> mu mu 
       Bs -> mu mu
       mu -> e gamma 
       mu -> eee

# NOTES
SEE pyMUNUSSM/MODIFICATIONS file for more info on how munussm and other subprogramms have been interfaced





