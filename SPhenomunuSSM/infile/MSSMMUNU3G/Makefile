modname = munuSSM3G  # as an example  
PreDef = -DGENERATIONMIXING 
# setting various paths  
InDir = ../include
Mdir = ${InDir}
name = ../lib/libSPhenomunuSSM3G.a
 
# check if SARAH module and SPheno are compatibel  
minV=330.00 
cVersion =$(shell expr $(version) \>= $(minV))
#  
# options for various compilers  
#  
# Default Compiler  
F90=gfortran
comp= -c -O -module ${Mdir} -I${InDir}  
LFlagsB= -O  
# Intels ifort,debug modus  
ifeq (${F90},ifortg)  
F90=ifort  
comp= -c -g -module ${Mdir} -I${InDir}  
LFlagsB= -g  
endif  
# gfortran  
ifeq (${F90},gfortran)  
#comp= -c -g -ffree-line-length-none -J${Mdir} -I${InDir}  
comp= -c -O0 -ffree-line-length-none -g -I${Mdir} -J../include -I${InDir}
LFlagsB= -g  
endif  
# g95  
ifeq (${F90},g95)  
comp= -c -O -fmod=${Mdir} -I${InDir}  
LFlagsB= -O  
endif  
# Lahey F95 compiler  
ifeq (${F90},lf95)  
comp=-c -O -M ${Mdir} -I${InDir}  
LFlagsB=-O  
endif  
# NAG f95/2003  
ifeq (${F90},nagfor)  
comp= -c -O -mdir ${Mdir} -I${InDir}  
LFlagsB= -O -DONLYDOUBLE -mdir ${MDir} -I${InDir}  
endif   
.SUFFIXES : .o .ps .f90 .F90 .a 
bin/SPhenomunuSSM3G:
ifeq (${cVersion},1)
	 cd ../src ; ${MAKE} F90=${F90} 
	 ${MAKE} F90=${F90} ${name} 
	 ${MAKE} F90=${F90} SPhenomunuSSM3G.o 
	 ${F90} -o SPhenomunuSSM3G ${LFlagsB} SPhenomunuSSM3G.o ../lib/libSPhenomunuSSM3G.a ../lib/libSPheno.a
	 mv SPhenomunuSSM3G ../bin
	 rm SPhenomunuSSM3G.o  
${name}:  ${name}(Model_Data_munuSSM3G.o)  \
 ${name}(RGEs_munuSSM3G.o)   \
 ${name}(Couplings_munuSSM3G.o) ${name}(SusyMasses_munuSSM3G.o) ${name}(TadpoleEquations_munuSSM3G.o) \
 ${name}(LoopCouplings_munuSSM3G.o) ${name}(CouplingsForDecays_munuSSM3G.o) \
 ${name}(SusyDecays_munuSSM3G.o) \
 ${name}(Glu_munuSSM3G.o)  ${name}(Sd_munuSSM3G.o)  ${name}(Su_munuSSM3G.o)  ${name}(Cha_munuSSM3G.o)  ${name}(Chi_munuSSM3G.o)  \
 ${name}(EffPotFunctions.o) ${name}(DerivativesEffPotFunctions.o) ${name}(EffectivePotential_munuSSM3G.o) \
 ${name}(2LPoleFunctions.o) ${name}(2LPole_munuSSM3G.o) \
 ${name}(BranchingRatios_munuSSM3G.o) ${name}(LoopMasses_munuSSM3G.o) \
 ${name}(HiggsCS_munuSSM3G.o) ${name}(RunSM_munuSSM3G.o) \
 ${name}(LowEnergy_munuSSM3G.o) \
${name}(FlavorKit_LFV_munuSSM3G.o) ${name}(FlavorKit_QFV_munuSSM3G.o) ${name}(FlavorKit_Observables_munuSSM3G.o)\
 ${name}(SugraRuns_munuSSM3G.o)  ${name}(InputOutput_munuSSM3G.o) 
else 
	 @echo -------------------------------------------------------------------  
	 @echo ERROR:  
	 @echo The installed SPheno is version not compatibel with this module 
	 @echo Please, upgrade at least to SPheno version 3.3.0.  
	 @echo The current SPheno version can be downloaded from 
	 @echo http://www.hepforge.org/downloads/spheno 
	 @echo ------------------------------------------------------------------- 
endif 
clean: 
	 rm -f *.o *~ */*.o */*~
cleanall: 
	 rm -f bin/SPheno3 lib/*.a *~ */*.o */*~ include/*
#
# Suffix rules
#
.f90.a:
	 ${F90} ${comp} $< 
	 ar -ruc $@ $*.o
	 rm -f $*.o
.F90.a:
	 ${F90} ${comp} ${PreDef} $< 
	 ar -ruc $@ $*.o
	 rm -f $*.o
.f90.o:
	 ${F90} ${comp} $<
.f90.ps:
	 a2ps -o $*.ps $<
.h.ps:
	 a2ps -o $*.ps $<
