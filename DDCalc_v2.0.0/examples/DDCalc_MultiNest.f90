!#######################################################################
! DDCALC PROGRAM (FORTRAN) -> MODIFIED TO USE MULTINEST
! 
! For various types of DM-nucleon interactions, as well as for various 
! direct detection experiments, it computes the expected signal rates,
! any by comparing to the observed number of events, also the 
! log(likelihood) value for each of the examples.
! In order to convert those likelihood values into an upper bound on
! e.g. the scattering cross section, please have a look at
! DDCalc_exclusionC.cpp.
! 
! Run:
!   ./DDCalc_MultiNest
!
!   Created by
!   Chris Savage    Nordita/Utah/SavageArchitectures
!   Pat Scott       Imperial College
!   Martin White    Adelaide
!   Felix Kahlhoefer   RWTH Aachen
!   Sebastian Wild     DESY
!   ddcalc@projects.hepforge.org
! 
!   Modified by
!   Andres Daniel Perez
!
!#######################################################################

PROGRAM DDCalc_MultiNest

  ! The module must be loaded via a 'USE DDCalc' statement
  ! in any routine that is to make use of DDCalc routines.
  ! To use anything but the default analysis also requires
  ! explicitly loading the DDExperiments module.
  USE DDCalc
  USE DDExperiments
  USE DDRates


  IMPLICIT NONE
  INTEGER :: i_bin,type
  REAL*8 :: mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD, rho_normalized
  REAL*8 :: DM_spin, fp, fn, ap, an
  double precision    :: mDM_Data, sigmap_SI_Data, sigman_SI_Data, sigmap_SD_Data, sigman_SD_Data, rho_normalized_Data

  ! These three sets of indices refer to instance of the three types that
  ! are the bedrock of DDCalc.  (Almost) every calculation needs to be
  ! provided with an instance of each of these to do its job.  Passing the
  ! index of one of them to DDCalc tells it which one in its internal cache
  ! to use for the calculation requested. You can create as many
  ! different instances as you want, corresponding to e.g. different
  ! detectors/analyses, WIMP models and DM halo models; the factory
  ! funcions create the instances in DDCalc and return you the index of
  ! the resulting object.
  TYPE(WIMPStruct) :: WIMP
  TYPE(HaloStruct) :: Halo
  TYPE(DetectorStruct) :: Detector

  open(unit = 100, file = 'multinestpoint.dat', status = 'old', action = 'read')
    read(100,*) mDM_Data
    read(100,*) sigmap_SI_Data
    read(100,*) sigman_SI_Data
    read(100,*) sigmap_SD_Data
    read(100,*) sigman_SD_Data
    read(100,*) rho_normalized_Data

  rho_normalized = rho_normalized_Data				! Local dark matter density [GeV/cm^3], normalized by the actual value of the relic density calculated by micromegas

  ! Initialise a DM Halo object to default values.
  Halo = DDCalc_InitHalo()
    
  ! Initialise a WIMP objects to default values.
  WIMP = DDCalc_InitWIMP()

  ! Optionally set the Standard Halo Model parameters to values different from the default choices:
  !     rho     Local dark matter density [GeV/cm^3]
  !     vrot    Local disk rotation speed [km/s]
  !     v0      Maxwell-Boltzmann most probable speed [km/s]
  !     vesc    Galactic escape speed [km/s]
  CALL DDCalc_SetSHM(Halo, rho_normalized,235d0,235d0,550d0)

  ! Optional: Use a tabulated velocity integral instead of the Standard Halo Model one
  !  CALL DDCalc_SetHalo(Halo,rho=0.3d0,g_file='Halos/gvmin_VL2_shell.txt',g_column=5,Nvmin=1000)


  ! ****************************************************************************************************************
  ! Example 1: XENON1T (2018) analysis, with standard SI/SD interactions specified
  !            by WIMP-nucleon cross sections.
  Detector = XENON1T_2018_Init()		    ! Initalize the XENON1T_2018 detector.
  mDM = mDM_Data                           	! DM Mass in GeV.
  sigmap_SI = sigmap_SI_Data				! SI WIMP-proton cross section in pb.
  sigman_SI = sigman_SI_Data				! SI WIMP-neutron cross section in pb.
						                    !   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = sigmap_SD_Data				! SD WIMP-proton cross section in pb.
  sigman_SD = sigman_SD_Data				! SD WIMP-neutron cross section in pb.
  CALL DDCalc_SetWIMP_msigma(WIMP, mDM, &
	sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Experiment 1: Mixed SI and SD interactions at XENON1T (2018)'
  WRITE (*,*) '           specified by the WIMP-nucleon cross sections.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'sigmap_SI = ', sigmap_SI, ' pb'
  WRITE (*,*) 'sigman_SI = ', sigman_SI, ' pb'
  WRITE (*,*) 'sigmap_SD = ', sigmap_SD, ' pb'
  WRITE (*,*) 'sigman_SD = ', sigman_SD, ' pb'
  WRITE (*,*) 'rho_normalized = ', rho_normalized, ' GeV/cm^3'
  WRITE (*,*)
  WRITE (*,*) 'XENON1T_2018 Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'XENON1T_2018 Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'XENON1T_2018 Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'XENON1T_2018 Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 3
  WRITE (*,*) 'XENON1T_2018 p-value(no bg subtraction):            ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 4
  WRITE (*,*) 'XENON1T_2018 p-value(background subtraction):       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************


  ! ****************************************************************************************************************
  ! Example 2: XENON1T (2017) analysis, with standard SI/SD interactions specified
  !            by WIMP-nucleon cross sections.
  Detector = XENON1T_2017_Init()		    ! Initalize the XENON1T_2017 detector.
  mDM = mDM_Data                           	! DM Mass in GeV.
  sigmap_SI = sigmap_SI_Data				! SI WIMP-proton cross section in pb.
  sigman_SI = sigman_SI_Data				! SI WIMP-neutron cross section in pb.
						                    !   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = sigmap_SD_Data				! SD WIMP-proton cross section in pb.
  sigman_SD = sigman_SD_Data				! SD WIMP-neutron cross section in pb.
  CALL DDCalc_SetWIMP_msigma(WIMP, mDM, &
	sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Experiment 2: Mixed SI and SD interactions at XENON1T (2017)'
  WRITE (*,*) '           specified by the WIMP-nucleon cross sections.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'sigmap_SI = ', sigmap_SI, ' pb'
  WRITE (*,*) 'sigman_SI = ', sigman_SI, ' pb'
  WRITE (*,*) 'sigmap_SD = ', sigmap_SD, ' pb'
  WRITE (*,*) 'sigman_SD = ', sigman_SD, ' pb'
  WRITE (*,*) 'rho_normalized = ', rho_normalized, ' GeV/cm^3'
  WRITE (*,*)
  WRITE (*,*) 'XENON1T_2017 Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'XENON1T_2017 Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'XENON1T_2017 Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'XENON1T_2017 Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 3
  WRITE (*,*) 'XENON1T_2017 p-value(no bg subtraction):            ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 4
  WRITE (*,*) 'XENON1T_2017 p-value(background subtraction):       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************


  ! ****************************************************************************************************************
  ! Example 3: PICO-60 (2017) analysis, with standard SI/SD interactions specified
  !            by WIMP-nucleon cross sections.
  Detector = PICO_60_2017_Init()        	! Initalize the PICO_60_2017 detector.
  mDM = mDM_Data                           	! DM Mass in GeV.
  sigmap_SI = sigmap_SI_Data				! SI WIMP-proton cross section in pb.
  sigman_SI = sigman_SI_Data				! SI WIMP-neutron cross section in pb.
						                    !   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = sigmap_SD_Data				! SD WIMP-proton cross section in pb.
  sigman_SD = sigman_SD_Data				! SD WIMP-neutron cross section in pb.
  CALL DDCalc_SetWIMP_msigma(WIMP, mDM, &
	sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Experiment 3: Mixed SI and SD interactions at PICO-60 (2017)'
  WRITE (*,*) '           specified by the WIMP-nucleon cross sections.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'sigmap_SI = ', sigmap_SI, ' pb'
  WRITE (*,*) 'sigman_SI = ', sigman_SI, ' pb'
  WRITE (*,*) 'sigmap_SD = ', sigmap_SD, ' pb'
  WRITE (*,*) 'sigman_SD = ', sigman_SD, ' pb'
  WRITE (*,*) 'rho_normalized = ', rho_normalized, ' GeV/cm^3'
  WRITE (*,*)
  WRITE (*,*) 'PICO_60_2017 Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'PICO_60_2017 Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'PICO_60_2017 Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'PICO_60_2017 Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 3
  WRITE (*,*) 'PICO_60_2017 p-value(no bg subtraction):            ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 4
  WRITE (*,*) 'PICO_60_2017 p-value(background subtraction):       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************



  ! ****************************************************************************************************************
  ! Example 4: CRESST (2017) analysis, with standard SI/SD interactions specified
  !            by WIMP-nucleon cross sections.
  Detector = CRESST_II_Init()           	! Initalize the CRESST_II detector.
  mDM = mDM_Data                           	! DM Mass in GeV.
  sigmap_SI = sigmap_SI_Data				! SI WIMP-proton cross section in pb.
  sigman_SI = sigman_SI_Data				! SI WIMP-neutron cross section in pb.
						                    !   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = sigmap_SD_Data				! SD WIMP-proton cross section in pb.
  sigman_SD = sigman_SD_Data				! SD WIMP-neutron cross section in pb.
  CALL DDCalc_SetWIMP_msigma(WIMP, mDM, &
	sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Experiment 4: Mixed SI and SD interactions at CRESST_II'
  WRITE (*,*) '           specified by the WIMP-nucleon cross sections.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'sigmap_SI = ', sigmap_SI, ' pb'
  WRITE (*,*) 'sigman_SI = ', sigman_SI, ' pb'
  WRITE (*,*) 'sigmap_SD = ', sigmap_SD, ' pb'
  WRITE (*,*) 'sigman_SD = ', sigman_SD, ' pb'
  WRITE (*,*) 'rho_normalized = ', rho_normalized, ' GeV/cm^3'
  WRITE (*,*)
  WRITE (*,*) 'CRESST_2017 Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'CRESST_2017 Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'CRESST_2017 Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'CRESST_2017 Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 3
  WRITE (*,*) 'CRESST_2017 p-value(no bg subtraction):            ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 4
  WRITE (*,*) 'CRESST_2017 p-value(background subtraction):       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************
    
     

  ! ****************************************************************************************************************
  ! Example 5: DARWIN analysis, with standard SI/SD interactions specified
  !            by WIMP-nucleon cross sections.
  Detector = DARWIN_Init()           	    ! Initalize the DARWIN detector.
  mDM = mDM_Data                           	! DM Mass in GeV.
  sigmap_SI = sigmap_SI_Data				! SI WIMP-proton cross section in pb.
  sigman_SI = sigman_SI_Data				! SI WIMP-neutron cross section in pb.
						                    !   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = sigmap_SD_Data				! SD WIMP-proton cross section in pb.
  sigman_SD = sigman_SD_Data				! SD WIMP-neutron cross section in pb.
  CALL DDCalc_SetWIMP_msigma(WIMP, mDM, &
	sigmap_SI, sigman_SI, sigmap_SD, sigman_SD)
  CALL DDCalc_CalcRates(Detector,WIMP,Halo)	! This performs the actual calculation of the rates.

  WRITE (*,*) '******************************************************************'
  WRITE (*,*) 'Experiment 5: Mixed SI and SD interactions at DARWIN'
  WRITE (*,*) '           specified by the WIMP-nucleon cross sections.'
  WRITE (*,*) 'mDM       = ', mDM, ' GeV'
  WRITE (*,*) 'sigmap_SI = ', sigmap_SI, ' pb'
  WRITE (*,*) 'sigman_SI = ', sigman_SI, ' pb'
  WRITE (*,*) 'sigmap_SD = ', sigmap_SD, ' pb'
  WRITE (*,*) 'sigman_SD = ', sigman_SD, ' pb'
  WRITE (*,*) 'rho_normalized = ', rho_normalized, ' GeV/cm^3'
  WRITE (*,*)
  WRITE (*,*) 'DARWIN Expected number of signal events:      ', DDCalc_Signal(Detector)
  WRITE (*,*) 'DARWIN Observed number of events:             ', DDCalc_Events(Detector)
  WRITE (*,*) 'DARWIN Expected number of background events:  ', DDCalc_Background(Detector)
  WRITE (*,*) 'DARWIN Log(likelihood):                       ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 3
  WRITE (*,*) 'DARWIN p-value(no bg subtraction):            ', DDCalc_LogLikelihood(Detector)
Detector%StatisticFlag = 4
  WRITE (*,*) 'DARWIN p-value(background subtraction):       ', DDCalc_LogLikelihood(Detector)
  WRITE (*,*) '******************************************************************'
  WRITE (*,*)
  ! ****************************************************************************************************************


  
END PROGRAM



!#######################################################################
! END OF FILE
!#######################################################################

