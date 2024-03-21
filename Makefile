#########################################################################
#                                                                       #
#     M a k e f i l e
#                                                                       #
#########################################################################


# SuperPy - compile the codes. You might need to tweak individual makefiles
# for your own system. You will need LAPACK, gcc compilers and python and python-dev.

default: all

all: ddcalc higgstools pymultinest multinest python spheno run


ddcalc:
	cd ./DDCalc_v*; make all
	

higgstools:
	rm -rf ./higgstools-main/buildpy/HBout3G/*
	# See the readme for proper configuration (you may want to build only python or C++)
	# Bulding python library (you may want to install it locally. check the HT installation instructions)
	cd higgstools-main/ && pip3 install . ;
	# Bulding C++ library
	# cd higgstools-main/ && mkdir -p build && cd build/ && cmake .. && make;
	# make running bash file executable
	cd higgstools-main/buildpy/  &&  chmod +x run_htools.sh


#micromegas:
	# TODO
	#make -C micromegas_*
	#make -C micromegas_*/NMSSMplusRHN main=OmegaMultiNestDD.cpp


pymultinest:
	#cd PyMultiNest-master/; python3 setup.py install
	cd PyMultiNest-master/ && pip3 install pymultinest


multinest:
	rm -rf MultiNest-master/build/* || true;
	# mkdir -p MultiNest-master/build;
	cd MultiNest-master/build/ &&  cmake .. && make
	# You probably will need to copy this line in the hidden file ./bashrc, with the correct path
	# or copy this line in every shell you use
	# export LD_LIBRARY_PATH=./superpy-munussm-v.0.0/MultiNest-master/lib/:$LD_LIBRARY_PATH


run:
	rm -rf ./pyMUNUSSM/temFiles/*
	cd pyMUNUSSM/; chmod +x run.sh


# Build the Python libraries. You might need to sudo these commands.
# Also, if you have your own machine, rather than a networked machine,
# you might want to install things globally, rather than locally, by
# removing the --user argument.


python:
	cd ./pyslha-*; python3 setup.py install --user
	cd ./PyMultiNest-master; python3 setup.py install --user


spheno:
	cd SPhenomunuSSM; make 


clean:
	cd ./DDCalc_v*; make clean
	-rm -rf ./MultiNest-master/build/*
	make -C SPhenomunuSSM cleanall
	rm -rf ./higgstools-main/buildpy/HBout3G/*
	rm -rf ./pyMUNUSSM/temFiles/*
	cd pyMUNUSSM/; chmod -x run.sh
	

clean_higgstools:
	rm -rf ./higgstools-main/buildpy/HBout3G/*
	cd higgstools-main/ && rm -rf HiggsTools.egg-info _skbuild

	
cleanall:	
	cd ./DDCalc_v*; make clean
	cd higgstools-main/ && rm -rf HiggsTools.egg-info  _skbuild
	-rm -rf ./MultiNest-master/build/*
	make -C SPhenomunuSSM cleanall
	rm -rf ./higgstools-main/buildpy/HBout3G/*
	rm -rf ./pyMUNUSSM/temFiles/*
	cd pyMUNUSSM/; chmod -x run.sh
