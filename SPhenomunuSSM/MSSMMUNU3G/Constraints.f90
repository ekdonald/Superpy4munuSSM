! File created at 14:42 on 26.8.2015   
! Last edited by Donald 0n 02/11/2015
! ----------------------------------------------------------------------  
 
 
subroutine constraints(PROB)
 
Use InputOutput_munuSSM3G
Use Model_Data_munuSSM3G

!! Use BranchingRatios_munuSSM3G !! by DK on 03-02-16
  
Implicit None

!!***********************************************************************
!      PROB(1) =/= 0   # Warning: selectron/smuon too light 
!      PROB(2) =/= 0   # Warning: stau too light 
!      PROB(3) =/= 0   # Warning: charged Higgs excluded by LEP
!      PROB(4) =/= 0   # Warning: chargino mass excluded by LEP
!      PROB(5) =/= 0   # Warning: squark/gluino too light
!      PROB(6) =/= 0   # Warning: scalar top is too light 
!      PROB(7) =/= 0   # Warning: scalar bottom is too light 
!      PROB(9) =/= 0   # Warning: Z decay too large, excluded
!      PROB(8) =/= 0   # Warning: 80% of the neutral or Charged Goldstone is not from Hu or Hd
!      PROB(10) =/= 0  # Warning: Total SM-like Higgs decay total width is too large, > 22MeV
!      PROB(11) =/= 0  # Warning: excluded by e e -> Z h independent flavor
!      PROB(12) =/= 0  # Warning: excluded by e e -> Z h, h -> b b
!      PROB(13) =/= 0  # Warning: excluded by e e -> Z h, h -> j j
!      PROB(14) =/= 0  # Warning: excluded by e e -> Z h, h -> tau tau
!      PROB(15) =/= 0  # Warning: excluded by e e -> Z h, h -> 2 photons 
!      PROB(16) =/= 0  # Warning: excluded by e e -> P Z independent flavor
!      PROB(17) =/= 0  # Warning: excluded by P -> b b
!      PROB(18) =/= 0  # Warning: excluded by P -> j j
!      PROB(19) =/= 0  # Warning: excluded by P -> tau tau
!      PROB(20) =/= 0  # Warning: excluded by P -> 2 photons 
!      PROB(21) =/= 0  # Warning: excluded by ee -> hZ-> ZAA/Zhh -> Z4bs 
!      PROB(22) =/= 0  # Warning: excluded by ee -> hZ, h-> hh/AA -> 4taus
!      PROB(41) =/= 0  # Warning: excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)
!      PROB(23) =/= 0  # Warning: excluded by ee -> hZ, h-> hh/ AA -> 2bs 2taus
!      PROB(24) =/= 0  # Warning: excluded by ee -> Zh -> ZAA -> Z + light pairs
!      PROB(25) =/= 0  # Warning: excluded by ee -> hA -> 4bs
!      PROB(26) =/= 0  # Warning: excluded by ee -> hA -> 4taus
!      PROB(27) =/= 0  # Warning: excluded by ee -> hA -> 2bs 2taus
!      PROB(28) =/= 0  # Warning: excluded by ee -> hA -> AAA -> 6bs
!      PROB(29) =/= 0  # Warning: excluded by ee -> hA -> AAA -> 6taus
!      PROB(30) =/= 0  # Warning: excluded by top -> b H+, H+ -> c s (CDF, D0)
!      PROB(31) =/= 0  # Warning: excluded by top -> b H+, H+ -> tau nu_tau (D0)
!      PROB(32) =/= 0  # Warning: excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
!      PROB(51) =/= 0  # Warning: excluded by H/A->tautau
!      PROB(52) =/= 0  # Warning: 80% of 125 GeV Higgs composition is not from Hu
!      PROB(53) =/= 0  # Warning: no Higgs in the 122-128 GeV range 
!      PROB(54) =/= 0  # Warning: excluded H_125->AA->4mu (CMS)
!      PROB(55) =/= 0  # Warning: 80% of the 125 Higgs composition is not a mixture of Hu & Hd
!      PROB(60) =/= 0  # Warning: excluded by stop -> b l sneutrino
!      PROB(61) =/= 0  # Warning: excluded by stop -> neutralino c
!      PROB(62) =/= 0  # Warning: excluded by sbottom -> neutralino b
!      PROB(64) =/= 0  # Warning: excluded by ggF->H/A->gamgam (65GeV < M < 122GeV, ATLAS)
!      PROB(65) =/= 0  # Warning: excluded by Upsilon(1S) -> A gamma
!      PROB(66) =/= 0  # Warning: excluded by eta_b(1S) mass measurement


!!*********** Declaration for scan and constraints ***********************

Complex(dp), save :: lam1INMin,lam1INMax,lam2INMin,lam2INMax,lam3INMin,lam3INMax
Real(dp)::lam1INMin_Re,lam1INMin_Im,lam1INMax_Re,lam1INMax_Im,lam2INMin_Re, &
& lam2INMin_Im,lam2INMax_Re,lam2INMax_Im,lam3INMin_Re,lam3INMin_Im, &
& lam3INMax_Re,lam3INMax_Im
!
Complex(dp), save:: Tlam1INMin,Tlam1INMax,Tlam2INMin,Tlam2INMax,Tlam3INMin,Tlam3INMax
Real(dp)::Tlam1INMin_Re,Tlam1INMin_Im,Tlam1INMax_Re,Tlam1INMax_Im,Tlam2INMin_Re, &
& Tlam2INMin_Im,Tlam2INMax_Re,Tlam2INMax_Im,Tlam3INMin_Re,Tlam3INMin_Im, &
& Tlam3INMax_Re,Tlam3INMax_Im
!
Complex(dp), save :: kap111INMax,kap112INMax,kap113INMax,kap121INMax,kap122INMax, &
& kap123INMax,kap131INMax,kap132INMax,kap133INMax,kap211INMax, &
& kap212INMax,kap213INMax,kap221INMax,kap222INMax,kap223INMax,kap231INMax, &
& kap232INMax,kap233INMax,kap311INMax,kap312INMax,kap313INMax, &
& kap321INMax,kap322INMax,kap323INMax,kap331INMax,kap332INMax,kap333INMax
!!
Complex(dp), save:: kap111INMin,kap112INMin,kap113INMin,kap121INMin,kap122INMin, &
& kap123INMin,kap131INMin,kap132INMin,kap133INMin,kap211INMin, &
& kap212INMin,kap213INMin,kap221INMin,kap222INMin,kap223INMin,kap231INMin, &
& kap232INMin,kap233INMin,kap311INMin,kap312INMin,kap313INMin, &
& kap321INMin,kap322INMin,kap323INMin,kap331INMin,kap332INMin,kap333INMin
!
Real(dp):: kap111INMax_Re,kap111INMax_Im,kap112INMax_Re, &
& kap112INMax_Im,kap113INMax_Re,kap113INMax_Im, &
& kap121INMax_Re,kap121INMax_Im,kap122INMax_Re,kap122INMax_Im,kap123INMax_Re, &
& kap123INMax_Im,kap131INMax_Re,kap131INMax_Im,kap132INMax_Re,kap132INMax_Im, &
& kap133INMax_Re,kap133INMax_Im,kap211INMax_Re,kap211INMax_Im,kap212INMax_Re, &
& kap212INMax_Im,kap213INMax_Re,kap213INMax_Im,kap221INMax_Re,kap221INMax_Im, &
& kap222INMax_Re,kap222INMax_Im,kap223INMax_Re,kap223INMax_Im,kap231INMax_Re, &
& kap231INMax_Im,kap232INMax_Re,kap232INMax_Im,kap233INMax_Re,kap233INMax_Im, &
& kap311INMax_Re,kap311INMax_Im,kap312INMax_Re,kap312INMax_Im,kap313INMax_Re, &
& kap313INMax_Im,kap321INMax_Re,kap321INMax_Im,kap322INMax_Re,kap322INMax_Im, &
& kap323INMax_Re,kap323INMax_Im,kap331INMax_Re,kap331INMax_Im,kap332INMax_Re, &
& kap332INMax_Im,kap333INMax_Re,kap333INMax_Im
!!
Real(dp):: kap111INMin_Re,kap111INMin_Im,kap112INMin_Re, &
& kap112INMin_Im,kap113INMin_Re,kap113INMin_Im, &
& kap121INMin_Re,kap121INMin_Im,kap122INMin_Re,kap122INMin_Im,kap123INMin_Re, &
& kap123INMin_Im,kap131INMin_Re,kap131INMin_Im,kap132INMin_Re,kap132INMin_Im, &
& kap133INMin_Re,kap133INMin_Im,kap211INMin_Re,kap211INMin_Im,kap212INMin_Re, &
& kap212INMin_Im,kap213INMin_Re,kap213INMin_Im,kap221INMin_Re,kap221INMin_Im, &
& kap222INMin_Re,kap222INMin_Im,kap223INMin_Re,kap223INMin_Im,kap231INMin_Re, &
& kap231INMin_Im,kap232INMin_Re,kap232INMin_Im,kap233INMin_Re,kap233INMin_Im, &
& kap311INMin_Re,kap311INMin_Im,kap312INMin_Re,kap312INMin_Im,kap313INMin_Re, &
& kap313INMin_Im,kap321INMin_Re,kap321INMin_Im,kap322INMin_Re,kap322INMin_Im, &
& kap323INMin_Re,kap323INMin_Im,kap331INMin_Re,kap331INMin_Im,kap332INMin_Re, &
& kap332INMin_Im,kap333INMin_Re,kap333INMin_Im
!
Complex(dp), save :: Tkap111INMax,Tkap112INMax,Tkap113INMax,Tkap121INMax,Tkap122INMax, &
& Tkap123INMax,Tkap131INMax,Tkap132INMax,Tkap133INMax,Tkap211INMax, &
& Tkap212INMax,Tkap213INMax,Tkap221INMax,Tkap222INMax,Tkap223INMax,Tkap231INMax, &
& Tkap232INMax,Tkap233INMax,Tkap311INMax,Tkap312INMax,Tkap313INMax, &
& Tkap321INMax,Tkap322INMax,Tkap323INMax,Tkap331INMax,Tkap332INMax,Tkap333INMax
!!
Complex(dp), save :: Tkap111INMin,Tkap112INMin,Tkap113INMin,Tkap121INMin,Tkap122INMin, &
& Tkap123INMin,Tkap131INMin,Tkap132INMin,Tkap133INMin,Tkap211INMin, &
& Tkap212INMin,Tkap213INMin,Tkap221INMin,Tkap222INMin,Tkap223INMin,Tkap231INMin, &
& Tkap232INMin,Tkap233INMin,Tkap311INMin,Tkap312INMin,Tkap313INMin, &
& Tkap321INMin,Tkap322INMin,Tkap323INMin,Tkap331INMin,Tkap332INMin,Tkap333INMin
!
Real(dp):: Tkap111INMax_Re,Tkap111INMax_Im,Tkap112INMax_Re, &
& Tkap112INMax_Im,Tkap113INMax_Re,Tkap113INMax_Im, &
& Tkap121INMax_Re,Tkap121INMax_Im,Tkap122INMax_Re,Tkap122INMax_Im,Tkap123INMax_Re, &
& Tkap123INMax_Im,Tkap131INMax_Re,Tkap131INMax_Im,Tkap132INMax_Re,Tkap132INMax_Im, &
& Tkap133INMax_Re,Tkap133INMax_Im,Tkap211INMax_Re,Tkap211INMax_Im,Tkap212INMax_Re, &
& Tkap212INMax_Im,Tkap213INMax_Re,Tkap213INMax_Im,Tkap221INMax_Re,Tkap221INMax_Im, &
& Tkap222INMax_Re,Tkap222INMax_Im,Tkap223INMax_Re,Tkap223INMax_Im,Tkap231INMax_Re, &
& Tkap231INMax_Im,Tkap232INMax_Re,Tkap232INMax_Im,Tkap233INMax_Re,Tkap233INMax_Im, &
& Tkap311INMax_Re,Tkap311INMax_Im,Tkap312INMax_Re,Tkap312INMax_Im,Tkap313INMax_Re, &
& Tkap313INMax_Im,Tkap321INMax_Re,Tkap321INMax_Im,Tkap322INMax_Re,Tkap322INMax_Im, &
& Tkap323INMax_Re,Tkap323INMax_Im,Tkap331INMax_Re,Tkap331INMax_Im,Tkap332INMax_Re, &
& Tkap332INMax_Im,Tkap333INMax_Re,Tkap333INMax_Im
!!
Real(dp):: Tkap111INMin_Re,Tkap111INMin_Im,Tkap112INMin_Re, &
& Tkap112INMin_Im,Tkap113INMin_Re,Tkap113INMin_Im, &
& Tkap121INMin_Re,Tkap121INMin_Im,Tkap122INMin_Re,Tkap122INMin_Im,Tkap123INMin_Re, &
& Tkap123INMin_Im,Tkap131INMin_Re,Tkap131INMin_Im,Tkap132INMin_Re,Tkap132INMin_Im, &
& Tkap133INMin_Re,Tkap133INMin_Im,Tkap211INMin_Re,Tkap211INMin_Im,Tkap212INMin_Re, &
& Tkap212INMin_Im,Tkap213INMin_Re,Tkap213INMin_Im,Tkap221INMin_Re,Tkap221INMin_Im, &
& Tkap222INMin_Re,Tkap222INMin_Im,Tkap223INMin_Re,Tkap223INMin_Im,Tkap231INMin_Re, &
& Tkap231INMin_Im,Tkap232INMin_Re,Tkap232INMin_Im,Tkap233INMin_Re,Tkap233INMin_Im, &
& Tkap311INMin_Re,Tkap311INMin_Im,Tkap312INMin_Re,Tkap312INMin_Im,Tkap313INMin_Re, &
& Tkap313INMin_Im,Tkap321INMin_Re,Tkap321INMin_Im,Tkap322INMin_Re,Tkap322INMin_Im, &
& Tkap323INMin_Re,Tkap323INMin_Im,Tkap331INMin_Re,Tkap331INMin_Im,Tkap332INMin_Re, &
& Tkap332INMin_Im,Tkap333INMin_Re,Tkap333INMin_Im
!
Complex(dp), save :: Yv11INMax,Yv12INMax,Yv13INMax,Yv21INMax,Yv22INMax, &
& Yv23INMax,Yv31INMax,Yv32INMax,Yv33INMax
Complex(dp):: Yv11INMin,Yv12INMin,Yv13INMin,Yv21INMin,Yv22INMin, &
& Yv23INMin,Yv31INMin,Yv32INMin,Yv33INMin
!
Real(dp):: Yv11INMax_Re,Yv11INMax_Im,Yv12INMax_Re,Yv12INMax_Im, &
& Yv13INMax_Re,Yv13INMax_Im,Yv21INMax_Re,Yv21INMax_Im,Yv22INMax_Re, &
& Yv22INMax_Im,Yv23INMax_Re,Yv23INMax_Im,Yv31INMax_Re,Yv31INMax_Im, &
& Yv32INMax_Re,Yv32INMax_Im,Yv33INMax_Re,Yv33INMax_Im
Real(dp):: Yv11INMin_Re,Yv11INMin_Im,Yv12INMin_Re,Yv12INMin_Im, &
& Yv13INMin_Re,Yv13INMin_Im,Yv21INMin_Re,Yv21INMin_Im,Yv22INMin_Re, &
& Yv22INMin_Im,Yv23INMin_Re,Yv23INMin_Im,Yv31INMin_Re,Yv31INMin_Im, &
& Yv32INMin_Re,Yv32INMin_Im,Yv33INMin_Re,Yv33INMin_Im
!
Complex(dp), save :: Tv11INMax,Tv12INMax,Tv13INMax,Tv21INMax,Tv22INMax, &
& Tv23INMax,Tv31INMax,Tv32INMax,Tv33INMax
Complex(dp), save :: Tv11INMin,Tv12INMin,Tv13INMin,Tv21INMin,Tv22INMin, &
& Tv23INMin,Tv31INMin,Tv32INMin,Tv33INMin
!
Real(dp):: Tv11INMax_Re,Tv11INMax_Im,Tv12INMax_Re,Tv12INMax_Im, &
& Tv13INMax_Re,Tv13INMax_Im,Tv21INMax_Re,Tv21INMax_Im,Tv22INMax_Re, &
& Tv22INMax_Im,Tv23INMax_Re,Tv23INMax_Im,Tv31INMax_Re,Tv31INMax_Im, &
& Tv32INMax_Re,Tv32INMax_Im,Tv33INMax_Re,Tv33INMax_Im
Real(dp):: Tv11INMin_Re,Tv11INMin_Im,Tv12INMin_Re,Tv12INMin_Im, &
& Tv13INMin_Re,Tv13INMin_Im,Tv21INMin_Re,Tv21INMin_Im,Tv22INMin_Re, &
& Tv22INMin_Im,Tv23INMin_Re,Tv23INMin_Im,Tv31INMin_Re,Tv31INMin_Im, &
& Tv32INMin_Re,Tv32INMin_Im,Tv33INMin_Re,Tv33INMin_Im
!
Complex(dp), save :: Tu11INMax,Tu12INMax,Tu13INMax,Tu21INMax,Tu22INMax, &
& Tu23INMax,Tu31INMax,Tu32INMax,Tu33INMax
Complex(dp), save :: Tu11INMin,Tu12INMin,Tu13INMin,Tu21INMin,Tu22INMin, &
& Tu23INMin,Tu31INMin,Tu32INMin,Tu33INMin
!
Real(dp):: Tu11INMax_Re,Tu11INMax_Im,Tu12INMax_Re,Tu12INMax_Im, &
& Tu13INMax_Re,Tu13INMax_Im,Tu21INMax_Re,Tu21INMax_Im,Tu22INMax_Re, &
& Tu22INMax_Im,Tu23INMax_Re,Tu23INMax_Im,Tu31INMax_Re,Tu31INMax_Im, &
& Tu32INMax_Re,Tu32INMax_Im,Tu33INMax_Re,Tu33INMax_Im
Real(dp):: Tu11INMin_Re,Tu11INMin_Im,Tu12INMin_Re,Tu12INMin_Im, &
& Tu13INMin_Re,Tu13INMin_Im,Tu21INMin_Re,Tu21INMin_Im,Tu22INMin_Re, &
& Tu22INMin_Im,Tu23INMin_Re,Tu23INMin_Im,Tu31INMin_Re,Tu31INMin_Im, &
& Tu32INMin_Re,Tu32INMin_Im,Tu33INMin_Re,Tu33INMin_Im
!
Real(dp), save ::vL1INMin,vL1INMax,vL2INMin,vL2INMax,vL3INMin,vL3INMax, &
& vR1INMin,vR1INMax,vR2INMin,vR2INMax,vR3INMin,vR3INMax, &
& tbINMin,tbINMax

Real(dp) :: msw2,mf1,mf2,mf3,mf4,mf5,mf6,mvu,mvd,smu 
integer :: Nmax, n, iEL, iMU, iTA, l, NscanPt, Nslep, iQ, ios,Nglusq,i,J,iH,jQ,kZ, K
integer,parameter :: seed = 86456 !scan DK 10.09.2015
Real(dp) :: MEL,MER,MMR,MML,MLL,MTT,MRR,Mhiggs,Mhiggs1,Mhiggs2,MSLMIN,MSTMIN, MTR, MTL
Real(dp) :: seLmix, smuLmix, stauLmix, seRmix, smuRmix, stauRmix, MMIN, MGLMIN
Real(dp) :: suLmix,suRmix,scLmix,scRmix,sdLmix,sdRmix,ssLmix,ssRmix,chDmix,chUmix
Real(dp) :: MDL,MDR,MGL,MUL,MUR,MSQMIN,TestSqGlu,MCHAMIN,TestChargino, MCMIN
Real(dp) :: CMASS, TestChHiggs, MCHu, MCHd
Real(dp) :: MBL,MTopL,MBR,MTopR, sbLmix, sbRmix, stLmix, stRmix, Testsb, Testst, MSTopMIN, MSBMIN
Real(dp) :: newGAMMAZsp,newGAMMAZchichi, FactorZsp, FactorZchichi,newGAMMAZ
complex(dp) :: OLnnz1,OLnnz2,OLnnz3,OLnnz,ORnnz1,ORnnz2,ORnnz3,ORnnz
complex(dp) :: Ospz1,Ospz2,Ospz3,Ospz
!!
 INTEGER, save :: NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
 INTEGER, save :: NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
 INTEGER, save :: NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
 INTEGER, save :: Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
 INTEGER, save :: Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
 INTEGER, save :: Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
 INTEGER, save :: Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
 INTEGER, save :: Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
 INTEGER, save :: Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
 INTEGER, save :: Nstblsn,Nstnc,Nsbnb,Nglsq, iZ, jZ ,kQ,iHp, iH1, iH2,iH125,iH126,iH127, N_lightPj,JPj

 DOUBLE PRECISION hZjj(1000,2), hZind(1000,2), hZinv(1000,2), hZll(1000, 2), hZbb(1000,2), hZgg(1000, 2)
 DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3),AAZ2b2tau(10000,3),cccc02(1000,2),cccc04(100,2)
 DOUBLE PRECISION cccc05(100,2), cccc06(100,2), cccc08(100,2), cccc1(100,2)
 DOUBLE PRECISION ccgg02(1000,2), ccgg04(100,2), ccgg05(100,2), ccgg06(100,2), ccgg08(100,2), ccgg1(100,2)
 DOUBLE PRECISION cctt02(1000,2), cctt04(100,2), cctt05(100,2), cctt06(100,2), cctt08(100,2), cctt1(100,2)
 DOUBLE PRECISION gggg02(1000,2), gggg04(100,2), gggg05(100,2), gggg06(100,2), gggg08(100,2), gggg1(100,2)
 DOUBLE PRECISION tttt02(1000,2), tttt04(100,2), tttt05(100,2), tttt06(100,2), tttt08(100,2), tttt1(100,2)
 DOUBLE PRECISION ttgg02(1000,2), ttgg04(100,2), ttgg05(100,2), ttgg06(100,2), ttgg08(100,2), ttgg1(100,2)
 DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3),hA2b2tau(10000,3),hA2tau2b(10000,3),AAA6b(10000,3),AAA6tau(10000,3)
 DOUBLE PRECISION stblsn(1000,2),stnc(1000,2),sbnb(1000,2)

 REAL(dp):: PROB(1000), BR4b, hu, hd, h1, h2, R, BRmax, Rmax, vtotSM2, vmunu2,tempMH,Rtest1,Rtest2,Rtest3
 REAL(dp):: PROBtest(100),GAMMAHtot,BR4bmax,BRPPmax,BRLLmax,BRWWmax,BRZZmax,BRBBmax,GammaHtotmax,MNL1,MNL2,MNL3,MNL4
 REAL(dp):: Massmax,Massmin,brtoplim,brtopcs,brtoptau,brtopa1
 COMPLEX(dp) :: FactorHww 

 DOUBLE PRECISION XIN(6),HM1(6),HM2(6),HM3(6),HM4(6),HM5(6)
 DOUBLE PRECISION DMA,HMIN,HMAX
 DOUBLE PRECISION minf(9),msup(9),binf(9),bsup(9)
 DOUBLE PRECISION mh4m(4),mh4p(4),br4m(4),br4p(4),brtau(8)
 DOUBLE PRECISION mh7(5),br7(5),br9(5),ceff,MAtest
 DOUBLE PRECISION alphaSamY,alphaSamZ
 DOUBLE PRECISION MY, MBQM,ALSMY,BRYMUMU,ZZ_Pj,AP,RC,DELTA,PI_Pj,ALEM0,RmaxPj,ALEMMVZ
 DOUBLE PRECISION M0Pj,GYEE,RETA,GAM2,XXPj,DPj,MEMAX,YMAX,FMAX,FPj,gPj,UUPj,VVPj,YYPj

 DATA XIN/.2d0,.25d0,.4d0,.5d0,.7d0,1d0/
 DATA HM1/70d0,94d0,99.7d0,102.3d0,105.2d0,108.4d0/
 DATA HM2/70d0,93d0,98.6d0,101.6d0,105d0,108d0/
 DATA HM3/70d0,92.3d0,98.3d0,100.3d0,104.6d0,107.6d0/
 DATA HM4/70d0,90.8d0,98d0,100.1d0,104d0,107d0/
 DATA HM5/70d0,86.7d0,97d0,99d0,103.5d0,106d0/
 DATA minf/60d0,70d0,80d0,100d0,110d0,120d0,130d0,140d0,150d0/
 DATA msup/70d0,80d0,100d0,110d0,120d0,130d0,140d0,150d0,155d0/
 DATA binf/.09d0,0d0,.21d0,.21d0,.15d0,.12d0,.08d0,.1d0,.20d0/
 DATA bsup/.12d0,0d0,.21d0,.15d0,.12d0,.08d0,.1d0,.13d0,.19d0/
 DATA mh4m/85d0,90d0,100d0,120d0/
 DATA mh4p/90d0,100d0,120d0,140d0/
 DATA br4m/.5d0,.33d0,.27d0,.35d0/
 DATA br4p/.33d0,.27d0,.35d0,.52d0/
 DATA mh7/90d0,100d0,120d0,140d0,160d0/
 DATA br7/.13d0,.1d0,.13d0,.2d0,.3d0/
 DATA br9/.11d0,.08d0,.09d0,.14d0,.21d0/
 DATA brtau/.16d0,.15d0,.16,.17d0,.175,.18d0,.19d0,.18d0/


  INTEGER I1,J1,NX,JBAR,JBARbb,NY
  PARAMETER(NX=18)
  PARAMETER(NY=6)

  DOUBLE PRECISION XX1(NY),LSIGBR(NY),LIMIT
  DOUBLE PRECISION D1,D2,CJ2,BRHTOAA
  DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz,MHcen

  DATA XX1/.25d0,.5d0,.75d0,1.d0,2.d0,3.d0/
  DATA LSIGBR/3.8d-3,3.6d-3,4.1d-3,4.2d-3,4.4d-3,4.6d-3/


  DOUBLE PRECISION HMAS(NX),XSM(NX),XSMbb(NX),LCMS(NX),LCMSbb(NX)
  DOUBLE PRECISION LATLASgg(NX),LATLASbb(NX)
  DOUBLE PRECISION MH(15),XSMH(15),LCMSH(15),SIG(15),LATLASH(15)
  DOUBLE PRECISION XSMHbb(15),LCMSHbb(15),SIGbb(15),LATLASHbb(15)
  DOUBLE PRECISION DEL,SIGTOT,MBAR,LCMSMB,LATLASMB
  DOUBLE PRECISION SIGTOTbb,MBARbb,LCMSMBbb,LATLASMBbb

      DATA HMAS/90d0,100d0,120d0,140d0,160d0,180d0,200d0,250d0,300d0,  &
     & 350d0,400d0,450d0,500d0,600d0,700d0,800d0,900d0,1000d0/

! SM Higgs ggF prod. cross sect. at 8 TeV from 
! https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
! CERNYellowReportPageAt8TeV#gluon_gluon_Fusion_Process
      DATA XSM/36.32d0,29.68d0,20.8d0,15.42d0,11.96d0,8.98d0,7.081d0, &
     & 4.783d0,3.594d0,3.401d0,2.921d0,2.002d0,1.283d0,.523d0,.229d0, &
     & .1097d0,.0571d0,.032d0/

! SM Higgs bbH prod. cross sect. at 8 TeV from 
! https://twiki.cern.ch/twiki/bin/view/LHCPhysics/
! /CrossSectionsFigures#MSSM_WG_plots (estimated)
      DATA XSMbb/0.56d0,0.42d0,0.25d0,0.15d0,8.7d-2,5.1d-2,3.2d-2, &
     & 1.2d-2,4.5d-3,2.8d-3,1.4d-3,8.2d-4,4.9d-4,2.6d-4,1.2d-4,  &
     & 5.2d-5,2.4d-5,1.2d-5/

! Upper limit on ggF->H->tautau (8 TeV) from CMS-PAS-HIG-13-021, Table 7
      DATA LCMS/50.2d0,31.3d0,7.38d0,2.27d0,.845d0,.549d0,.517d0,.315d0, &
     & .15d0,.112d0,.103d0,.607d-1,.385d-1,.193d-1,.143d-1,.115d-1, &
     & .923d-2,.865d-2/

! Upper limit on bbH->tautau (8 TeV) from CMS-PAS-HIG-13-021, Table 8
      DATA LCMSbb/6.03d0,4.14d0,1.76d0,1.25d0,.814d0,.659d0,.553d0, &
     & .217d0,.975d-1,.638d-1,.613d-1,.431d-1,.320d-1,.203d-1,.173d-1,.166d-1, &
     & .146d-1,.133d-1/ 

! Upper limit on ggF->H->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASgg/29.1d0,24.0d0,5.25d0,2.02d0,1.39d0,1.00d0,.794d0, &
     &  .281d0,.127d0,.112d0,.773d-1,.400d-1,.240d-1,.177d-1,.127d-1, &
     &  .993d-2,.840d-2,.735d-2/

! Upper limit on bbH->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASbb/6.32d0,6.32d0,2.73d0,1.27d0,.966d0,.606d0,.393d0, &
     &  .305d0,.116d0,.101d0,.656d-1,.363d-1,.238d-1,.159d-1,.117d-1, &
     &  .943d-2,.785d-2,.716d-2/


 DOUBLE PRECISION ggHgg(100,4),dummy(100,4),SMXS(100,6)
 Integer :: NggHgg,NSMXS 


!!*************************************************************** *************
!!*****************************************************************************

Real(dp) :: epsI=0.00001_dp, deltaM = 0.000001_dp 
Real(dp) :: mGut = -1._dp, ratioWoM = 0._dp
Integer :: kont 
 
Integer,Parameter :: p_max=100
Real(dp) :: Ecms(p_max),Pm(p_max),Pp(p_max), dt, tz, Qin, gSM(11) 
Real(dp) :: vev, sinw2
Logical :: ISR(p_max)=.False.
Logical :: CalcTBD
Real(dp) :: ae,amu,atau,EDMe,EDMmu,EDMtau,dRho,BrBsGamma,ratioBsGamma,BrDmunu,ratioDmunu,         & 
& BrDsmunu,ratioDsmunu,BrDstaunu,ratioDstaunu,BrBmunu,ratioBmunu,BrBtaunu,               & 
& ratioBtaunu,BrKmunu,ratioKmunu,RK,RKSM,muEgamma,tauEgamma,tauMuGamma,CRmuEAl,          & 
& CRmuETi,CRmuESr,CRmuESb,CRmuEAu,CRmuEPb,BRmuTo3e,BRtauTo3e,BRtauTo3mu,BRtauToemumu,    & 
& BRtauTomuee,BRtauToemumu2,BRtauTomuee2,BrZtoMuE,BrZtoTauE,BrZtoTauMu,BrhtoMuE,         & 
& BrhtoTauE,BrhtoTauMu,DeltaMBs,ratioDeltaMBs,DeltaMBq,ratioDeltaMBq,BrTautoEPi,         & 
& BrTautoEEta,BrTautoEEtap,BrTautoMuPi,BrTautoMuEta,BrTautoMuEtap,BrB0dEE,               & 
& ratioB0dEE,BrB0sEE,ratioB0sEE,BrB0dMuMu,ratioB0dMuMu,BrB0sMuMu,ratioB0sMuMu,           & 
& BrB0dTauTau,ratioB0dTauTau,BrB0sTauTau,ratioB0sTauTau,BrBtoSEE,ratioBtoSEE,            & 
& BrBtoSMuMu,ratioBtoSMuMu,BrBtoKmumu,ratioBtoKmumu,BrBtoSnunu,ratioBtoSnunu,            & 
& BrBtoDnunu,ratioBtoDnunu,BrKptoPipnunu,ratioKptoPipnunu,BrKltoPinunu,ratioKltoPinunu,  & 
& DelMK,ratioDelMK,epsK,ratioepsK

  character(8)  :: date
  character(10) :: time
  character(5)  :: zone

 CHARACTER(256) :: EXPCON_PATH

Real(dp) :: CLEOTAU, CLEOMU

logical, save :: first = .true.

logical :: debug = .false.

!!!**************************************************************
Real(dp) :: mvR(3),mneuvR(3),mneuBino(3),mneuWino(3),mBino,mWino
Integer :: ip


!! if(debug) then   !!! Not to use this subroutine


!Real(dp) :: DeltaM_squarks, DeltaMSu(6), DeltaMSd(6)
!!!**************************************************************

!The EXPCON_PATH variable is set:

      !CALL getenv('EXPCON_PATH',EXPCON_PATH)
      !if(EXPCON_PATH.eq.' ')  EXPCON_PATH='../EXPCON'

 CALL getcwd(EXPCON_PATH)

EXPCON_PATH=trim(EXPCON_PATH)//'/contributions/SPheno-4.0.0/data/'

!!-------------- Neutral 125 GeV Higgs (122 < Mh < 128 GeV)-------------


  call date_and_time(DATE=date,ZONE=zone)
  call date_and_time(TIME=time)

!--------------------------------------------------------------
!!---------- Writting constraints tests output ----------------
!--------------------------------------------------------------


 OPEN(unit = 17, file = trim(EXPCON_PATH)//'OutputConstraints.dat') 

 WRITE(17,*) "# ***** Test For LEP+LHC+TEVATRON Contraints ******"
 WRITE(17,*) "        SPheno # Spectrum calculator "
 WRITE(17,*) "        3.3.6      # Version number "
 write(17,*)  "       Started: ", date(7:8),".",date(5:6),".",date(1:4)," , ",time(1:2),":",time(3:4),":",time(5:6)
 WRITE(17,*) "      "

!
 MHmin = 122d0
 MHmax = 128d0
DO iEL = 1, 8
 If( (Mhh(iEL).gt. MHmin) .and. (Mhh(iEL).lt. MHmax) ) THEN  
   iH125 = iH125 + 1
!   WRITE(17,*)  ZH(iEL,1),' ',  ZH(iEL,2),' ', ZH(iEL,1)**2 + ZH(iEL,2)**2   !! DK
!   WRITE(17,*) iEL, '  ',  Mhh(iEL),' ', BRhh(iEL,184)  !! DK
!   WRITE(17,*)  '  ' !! DK

    PROBtest(1) = DDIM(1d0, (ZH(iEL,1)**2 + ZH(iEL,2)**2 )/ 7.d-1) 
   If( PROBtest(1).eq. 0d0 ) THEN      !! 80% of Higgs composition must be a mixture of Hu & Hd 
     iH126 = iH126 +1
      PROBtest(2) = DDIM(1d0, DABS(ZH(iEL,2))/ 8.d-1)
       If( PROBtest(2).eq. 0d0) THEN !! 80% of the true 125 GeV Higgs composition must be Hu 
!         iH127 = iH127 +1             
         iH = iEL  
         Mhiggs = Mhh(iH)
       Endif
   EndIf
 Endif
 ENDDO


 WRITE(17,*) ' Neutral Scalars, pseudoscalars masses ' 
 DO  iEL = 2, 8
 WRITE(17,*) iEL, '  ',  Mhh(iEL),  MAh(iEL) 
 ENDDO

 WRITE(17,*) '  ' 
 WRITE(17,*) ' Charged Higgs, charginos masses '
 DO  iEL = 2, 8
  WRITE(17,*) iEL, '  ',  MHpm(iEL), Mcha(iEL)
 ENDDO
 
 WRITE(17,*) '  ' 
 WRITE(17,*) ' Neutalinos  '
 DO  iEL = 1, 10
  WRITE(17,*) iEL, '  ',  Mchi(iEL)
 ENDDO

 OPEN(unit = 23, file = trim(EXPCON_PATH)//'MhForKapLam.dat')
 OPEN(unit = 24, file = trim(EXPCON_PATH)//'KappaLambda.dat')

!!-----------------------------------------------------------------------------------
  hd = 1d0/DSQRT(2d0*DSQRT(2d0)*(1d0+Real(TanBeta,dp)**2)*G_F) ! = vd/sqrt(2)
  hu = hd*Real(TanBeta,dp)                                     ! = vu/sqrt(2)
!!******************************************** Test DK
 mBino = Real(M1,dp)
 mWino = Real(M2,dp)
 mvR(1)= 2* Real(kap(1,1,1),dp)*vR(1)
 mvR(2)= 2* Real(kap(2,2,2),dp)*vR(2)
 mvR(3)= 2* Real(kap(3,3,3),dp)*vR(3)

Do ip =1, 3
mneuBino(ip) = (vL(ip)*Real(g2,dp))**2 / mBino
mneuWino(ip) = (vL(ip)*Real(g2,dp))**2 / mWino
mneuvR(ip)   = (Real(Yv(ip,ip),dp)*hu)**2 / mvR(ip)
Enddo
  WRITE(23,*)  MChi(1), MChi(2), MChi(3)
  WRITE(23,*)  mBino, mWino
  WRITE(23,*)  mneuBino(1), mneuBino(2), mneuBino(3)
  WRITE(23,*)  mneuWino(1), mneuWino(2), mneuWino(3)
  WRITE(23,*)  mvR(1), mvR(2), mvR(3)
  WRITE(23,*)  mneuvR(1), mneuvR(2), mneuvR(3)
  WRITE(23,*)  ' '
  WRITE(23,*)  ' '
!!******************************************** Test DK

 OPEN(unit = 20, file = trim(EXPCON_PATH)//'MhForvLYv.dat')
  WRITE(20,*) Mhiggs, vL(1), Real(Yv(1,1), dp)


!!-----------------------------------------------------
!init
PROB = 0.d0

!!----------- Test on chargino mass -------------------
!!-----------------------------------------------------
 PROB(4) = DDIM(1d0,DABS(MCha(4))/103.5)

!!-----------------------------------------------------

!!--- set debug = .false. in order to not use my subroutine

!!init
!!PROB = 0.d0

if(first) then 


!-----------------------------------------------------------------
!--- Reading the scan input files --------------
!-----------------------------------------------
!Open(unit = 98, file = trim(EXPCON_PATH)//'Rinput_data', status = 'old')!for real
!Open(unit = 99, file = trim(EXPCON_PATH)//'Cinput_data', status = 'old')!for complex

!

  call date_and_time(DATE=date,ZONE=zone)
  call date_and_time(TIME=time)

!--------------------------------------------------------------
!!---------- Writting constraints tests output ----------------
!--------------------------------------------------------------


if(debug) then
 OPEN(unit = 17, file = trim(EXPCON_PATH)//'OutputConstraints.dat') 

 WRITE(17,*) "# ***** Test For LEP+LHC+TEVATRON Contraints ******"
 WRITE(17,*) "        SPheno # Spectrum calculator "
 WRITE(17,*) "        3.3.6      # Version number "
 write(17,*)  "       Started: ", date(7:8),".",date(5:6),".",date(1:4)," , ",time(1:2),":",time(3:4),":",time(5:6)
 WRITE(17,*) "      "
!endif

!!------------------------------------------------------------
!!-------- Openning Mhh in (Kappa, lamba) plane -------
!------------------------------------------------------
!OPEN(unit = 20, file = '../MhForKapLam1.dat')
!OPEN(unit = 21, file = '../MhForKapLam2.dat')
!OPEN(unit = 22, file = '../MhForKapLam3.dat')

!only first time so I have to define static variables a la NMSSMBayes-like
 OPEN(unit = 23, file = trim(EXPCON_PATH)//'MhForKapLam.dat')
 OPEN(unit = 24, file = trim(EXPCON_PATH)//'KappaLambda.dat')

! OPEN(unit = 20, file = trim(EXPCON_PATH)//'SusyMassPrecision.dat')

!! -----------------------------------------------------------
!!---------- Reading InputConstraintsMin.dat file -------
!!-------------------------------------------------------

!-------------------------------------------------------
! Reading from LEP+Tevaton+lhc constraints files 
!-------------------------------------------------------

OPEN(unit = 205, file = trim(EXPCON_PATH)//'EXPCON/hZinv.dat', status='old') !
NhZinv = 0 
DO iQ=1, 1000! NhZinv
  READ(205,*,IOSTAT=ios) hZinv(iQ,1), hZinv(iQ,2)
  IF (ios /= 0) EXIT
  NhZinv = NhZinv +1
ENDDO
!!-------------------------------------------------------------------
OPEN(unit = 206, file = trim(EXPCON_PATH)//'EXPCON/hZind.dat', status='old') !
NhZind = 0 
DO iQ=1, 1000
  READ(206,*,IOSTAT=ios) hZind(iQ,1), hZind(iQ,2)
  IF (ios /= 0) EXIT
  NhZind = NhZind +1
ENDDO
!!-------------------------------------------------------------------
NhZbb = 0
OPEN(unit = 207, file = trim(EXPCON_PATH)//'EXPCON/hZbb.dat', status='old') !
DO iQ=1, 1000
  READ(207,*,IOSTAT=ios) hZbb(iQ,1), hZbb(iQ,2)
  IF (ios /= 0) EXIT
  NhZbb = NhZbb+1
ENDDO

!------------------------------------------------------- h -> 2jets
OPEN(unit = 208, file = trim(EXPCON_PATH)//'EXPCON/hZjj.dat', status='old') !
NhZjj = 0 
DO iQ=1, 1000
  READ(208,*,IOSTAT=ios) hZjj(iQ,1), hZjj(iQ,2)
  IF (ios /= 0) EXIT
  NhZjj = NhZjj+1
ENDDO
!!----------------------------------------------------- h -> tau tau
OPEN(unit = 209, file = trim(EXPCON_PATH)//'EXPCON/hZll.dat', status='old') !
NhZll = 0 
DO iQ=1, 1000
  READ(209,*,IOSTAT=ios) hZll(iQ,1), hZll(iQ,2)
  IF (ios /= 0) EXIT
  NhZll = NhZll+1
ENDDO
!!----------------------------------------------------- h -> 2 photons
NhZgg = 0 
OPEN(unit = 210, file = trim(EXPCON_PATH)//'EXPCON/hZgg.dat', status='old') !
DO iQ=1, 100
  READ(210,*,IOSTAT=ios) hZgg(iQ,1), hZgg(iQ,2)
  IF (ios /= 0) EXIT
  NhZgg =  NhZgg +1
ENDDO

!!----------------------------------------------------- h -> AA -> 4b
NAAZ4b = 0 
OPEN(unit = 211, file = trim(EXPCON_PATH)//'EXPCON/AAZ4b.dat', status='old') !
DO iQ=1, 5000
  READ(211,*,IOSTAT=ios) AAZ4b(iQ,1), AAZ4b(iQ,2), AAZ4b(iQ,3)
  IF (ios /= 0) EXIT
  NAAZ4b = NAAZ4b +1
ENDDO
!!----------------------------------------------------- h -> AA -> 4tau
NAAZ4tau = 0 
OPEN(unit = 212, file = trim(EXPCON_PATH)//'EXPCON/AAZ4tau.dat', status='old') !
DO iQ=1, 2000
  READ(212,*,IOSTAT=ios) AAZ4tau(iQ,1), AAZ4tau(iQ,2), AAZ4tau(iQ,3)
  IF (ios /= 0) EXIT
  NAAZ4tau = NAAZ4tau +1
ENDDO
!!----------------------------------------- h-> AA -> 2b+2tau
NAAZ2b2tau = 0 
OPEN(unit = 213, file = trim(EXPCON_PATH)//'EXPCON/AAZ2b2tau.dat', status='old') !
DO iQ=1, 3000
  READ(213,*,IOSTAT=ios) AAZ2b2tau(iQ,1), AAZ2b2tau(iQ,2), AAZ2b2tau(iQ,3)
  IF (ios /= 0) EXIT
  NAAZ2b2tau = NAAZ2b2tau +1
ENDDO
!!------------------------------------------- h -> AA -> 4c
Ncccc02 =  0 
OPEN(unit = 214, file = trim(EXPCON_PATH)//'EXPCON/cccc02.dat', status='old') !
DO iQ=1, 50
  READ(214,*,IOSTAT=ios) cccc02(iQ,1), cccc02(iQ,2)
  IF (ios /= 0) EXIT
  Ncccc02 = Ncccc02 + 1
ENDDO
!!-------------------------------------------------
Ncccc04 = 0
OPEN(unit = 215, file = trim(EXPCON_PATH)//'EXPCON/cccc04.dat', status='old') !
DO iQ=1, 30
  READ(215,*,IOSTAT=ios) cccc04(iQ,1), cccc04(iQ,2)
  IF (ios /= 0) EXIT
  Ncccc04 = Ncccc04 + 1
ENDDO
!!-------------------------------------------------
Ncccc05 = 0
OPEN(unit = 216, file = trim(EXPCON_PATH)//'EXPCON/cccc05.dat', status='old') !
DO iQ=1, 30
  READ(216,*,IOSTAT=ios) cccc05(iQ,1), cccc05(iQ,2)
  IF (ios /= 0) EXIT
  Ncccc05 = Ncccc05 + 1
ENDDO
!!-------------------------------------------------
Ncccc06 = 0 
OPEN(unit = 217, file = trim(EXPCON_PATH)//'EXPCON/cccc06.dat', status='old') !
DO iQ=1, 30
  READ(217,*,IOSTAT=ios) cccc06(iQ,1), cccc06(iQ,2)
  IF (ios /= 0) EXIT
  Ncccc06 = Ncccc06 +1
ENDDO
!!-------------------------------------------------
Ncccc08 = 0
OPEN(unit = 218, file = trim(EXPCON_PATH)//'EXPCON/cccc08.dat', status='old') !
DO iQ=1, 30
  READ(218,*,IOSTAT=ios) cccc08(iQ,1), cccc08(iQ,2)
  IF (ios /= 0) EXIT
  Ncccc08 = Ncccc08 +1
ENDDO
!!-------------------------------------------------
Ncccc1 = 0
OPEN(unit = 219, file = trim(EXPCON_PATH)//'EXPCON/cccc1.dat', status='old') !
DO iQ=1, 30
  READ(219,*,IOSTAT=ios) cccc1(iQ,1), cccc1(iQ,2)
  IF (ios /= 0) EXIT
  Ncccc1 = Ncccc1 +1
ENDDO

!!-------------------------------------------------
Nccgg02 = 0
OPEN(unit = 220, file = trim(EXPCON_PATH)//'EXPCON/ccgg02.dat', status='old') !
DO iQ=1, 1000
  READ(220,*,IOSTAT=ios) ccgg02(iQ,1), ccgg02(iQ,2)
  IF (ios /= 0) EXIT
  Nccgg02 = Nccgg02 +1
ENDDO

!!-------------------------------------------------
Nccgg04 = 0
OPEN(unit = 221, file = trim(EXPCON_PATH)//'EXPCON/ccgg04.dat', status='old') !
DO iQ=1, 1000
  READ(221,*,IOSTAT=ios) ccgg04(iQ,1), ccgg04(iQ,2)
  IF (ios /= 0) EXIT
  Nccgg04 = Nccgg04 +1
ENDDO

!!-------------------------------------------------
Nccgg05 = 0
OPEN(unit = 222, file = trim(EXPCON_PATH)//'EXPCON/ccgg05.dat', status='old') !
DO iQ=1, 1000
  READ(222,*,IOSTAT=ios) ccgg05(iQ,1), ccgg05(iQ,2)
  IF (ios /= 0) EXIT
  Nccgg05 = Nccgg05 +1
ENDDO

!!-------------------------------------------------
Nccgg06 = 0 
OPEN(unit = 223, file = trim(EXPCON_PATH)//'EXPCON/ccgg06.dat', status='old') !
DO iQ=1, 1000
  READ(223,*,IOSTAT=ios) ccgg06(iQ,1), ccgg06(iQ,2)
  IF (ios /= 0) EXIT
  Nccgg06 = Nccgg06 +1
ENDDO

!!-------------------------------------------------
Nccgg08 = 0
OPEN(unit = 224, file = trim(EXPCON_PATH)//'EXPCON/ccgg08.dat', status='old') !
DO iQ=1, 1000
  READ(224,*,IOSTAT=ios) ccgg08(iQ,1), ccgg08(iQ,2)
  IF (ios /= 0) EXIT
  Nccgg08 = Nccgg08 +1
ENDDO

!!-------------------------------------------------
Nccgg1 = 0
OPEN(unit = 225, file = trim(EXPCON_PATH)//'EXPCON/ccgg1.dat', status='old') !
DO iQ=1, 1000
  READ(225,*,IOSTAT=ios) ccgg1(iQ,1), ccgg1(iQ,2)
  IF (ios /= 0) EXIT
  Nccgg1 = Nccgg1 +1
ENDDO

!!-------------------------------------------
Ncctt02 =  0 
OPEN(unit = 226, file = trim(EXPCON_PATH)//'EXPCON/cctt02.dat', status='old') !
DO iQ=1, 300
  READ(226,*,IOSTAT=ios) cctt02(iQ,1), cctt02(iQ,2)
  IF (ios /= 0) EXIT
  Ncctt02 = Ncctt02 + 1
ENDDO
!!-------------------------------------------------
Ncctt04 = 0 
OPEN(unit = 227, file = trim(EXPCON_PATH)//'EXPCON/cctt04.dat', status='old') !
DO iQ=1, 300
  READ(227,*,IOSTAT=ios) cctt04(iQ,1), cctt04(iQ,2)
  IF (ios /= 0) EXIT
  Ncctt04 = Ncctt04 + 1
ENDDO
!!-------------------------------------------------
Ncctt05 = 0 
OPEN(unit = 228, file = trim(EXPCON_PATH)//'EXPCON/cctt05.dat', status='old') !
DO iQ=1, 300
  READ(228,*,IOSTAT=ios) cctt05(iQ,1), cctt05(iQ,2)
  IF (ios /= 0) EXIT
  Ncctt05 = Ncctt05 + 1
ENDDO
!!-------------------------------------------------
Ncctt06 = 0 ! Read ATLAS upper limit
! ggHgg(IQ,1): Higgs mass
! ggHgg(IQ,4): upper limit in fb

OPEN(unit = 229, file = trim(EXPCON_PATH)//'EXPCON/cctt06.dat', status='old') !
DO iQ=1, 300
  READ(229,*,IOSTAT=ios) cctt06(iQ,1), cctt06(iQ,2)
  IF (ios /= 0) EXIT
  Ncctt06 = Ncctt06 +1
ENDDO
!!-------------------------------------------------
Ncctt08 = 0 
OPEN(unit = 230, file = trim(EXPCON_PATH)//'EXPCON/cctt08.dat', status='old') !
DO iQ=1, 300
  READ(230,*,IOSTAT=ios) cctt08(iQ,1), cctt08(iQ,2)
  IF (ios /= 0) EXIT
  Ncctt08 = Ncctt08 +1
ENDDO
!!-------------------------------------------------
Ncctt1 = 0 
OPEN(unit = 231, file = trim(EXPCON_PATH)//'EXPCON/cctt1.dat', status='old') !
DO iQ=1, 300
  READ(231,*,IOSTAT=ios) cctt1(iQ,1), cctt1(iQ,2)
  IF (ios /= 0) EXIT
  Ncctt1 = Ncctt1 +1
ENDDO
!!----------------------------------------------------
Ngggg02 =  0 
OPEN(unit = 232, file = trim(EXPCON_PATH)//'EXPCON/gggg02.dat', status='old') !
DO iQ=1, 300
  READ(232,*,IOSTAT=ios) gggg02(iQ,1), gggg02(iQ,2)
  IF (ios /= 0) EXIT
  Ngggg02 = Ngggg02 + 1
ENDDO
!!-------------------------------------------------
Ngggg04 = 0 
OPEN(unit = 233, file = trim(EXPCON_PATH)//'EXPCON/gggg04.dat', status='old') !
DO iQ=1, 300
  READ(233,*,IOSTAT=ios) gggg04(iQ,1), gggg04(iQ,2)
  IF (ios /= 0) EXIT
  Ngggg04 = Ngggg04 + 1
ENDDO
!!-------------------------------------------------
Ngggg05 = 0 
OPEN(unit = 234, file = trim(EXPCON_PATH)//'EXPCON/gggg05.dat', status='old') !
DO iQ=1, 300
  READ(234,*,IOSTAT=ios) gggg05(iQ,1), gggg05(iQ,2)
  IF (ios /= 0) EXIT
  Ngggg05 = Ngggg05 + 1
ENDDO
!!-------------------------------------------------
Ngggg06 = 0
OPEN(unit = 235, file = trim(EXPCON_PATH)//'EXPCON/gggg06.dat', status='old') !
DO iQ=1, 300
  READ(235,*,IOSTAT=ios) gggg06(iQ,1), gggg06(iQ,2)
  IF (ios /= 0) EXIT
  Ngggg06 = Ngggg06 +1
ENDDO
!!-------------------------------------------------
Ngggg08 = 0 
OPEN(unit = 236, file = trim(EXPCON_PATH)//'EXPCON/gggg08.dat', status='old') !
DO iQ=1, 300
  READ(236,*,IOSTAT=ios) gggg08(iQ,1), gggg08(iQ,2)
  IF (ios /= 0) EXIT
  Ngggg08 = Ngggg08 +1
ENDDO
!!-------------------------------------------------
Ngggg1 = 0 
OPEN(unit = 237, file = trim(EXPCON_PATH)//'EXPCON/gggg1.dat', status='old') !
DO iQ=1, 300
  READ(237,*,IOSTAT=ios) gggg1(iQ,1), gggg1(iQ,2)
  IF (ios /= 0) EXIT
  Ngggg1 = Ngggg1 +1
ENDDO
!!-------------------------------------------
Nttgg02 = 0 
OPEN(unit = 238, file = trim(EXPCON_PATH)//'EXPCON/ttgg02.dat', status='old') !
DO iQ=1, 300
  READ(238,*,IOSTAT=ios) ttgg02(iQ,1), ttgg02(iQ,2)
  IF (ios /= 0) EXIT
  Nttgg02 = Nttgg02 + 1
ENDDO
!!-------------------------------------------------
Nttgg04 = 0 
OPEN(unit = 239, file = trim(EXPCON_PATH)//'EXPCON/ttgg04.dat', status='old') !
DO iQ=1, 300
  READ(239,*,IOSTAT=ios) ttgg04(iQ,1), ttgg04(iQ,2)
  IF (ios /= 0) EXIT
  Nttgg04 = Nttgg04 + 1
ENDDO
!!-------------------------------------------------
Nttgg05 = 0 
OPEN(unit = 240, file = trim(EXPCON_PATH)//'EXPCON/ttgg05.dat', status='old') !
DO iQ=1, 300
  READ(240,*,IOSTAT=ios) ttgg05(iQ,1), ttgg05(iQ,2)
  IF (ios /= 0) EXIT
  Nttgg05 = Nttgg05 + 1
ENDDO
!!-------------------------------------------------
Nttgg06 = 0 
OPEN(unit = 241, file = trim(EXPCON_PATH)//'EXPCON/ttgg06.dat', status='old') !
DO iQ=1, 300
  READ(241,*,IOSTAT=ios) ttgg06(iQ,1), ttgg06(iQ,2)
  IF (ios /= 0) EXIT
  Nttgg06 = Nttgg06 +1 
ENDDO
!!-------------------------------------------------
Nttgg08 = 0
OPEN(unit = 242, file = trim(EXPCON_PATH)//'EXPCON/ttgg08.dat', status='old') !
DO iQ=1, 300
  READ(242,*,IOSTAT=ios) ttgg08(iQ,1), ttgg08(iQ,2)
  IF (ios /= 0) EXIT
  Nttgg08 = Nttgg08 +1
ENDDO
!!-------------------------------------------------
Nttgg1 = 0 
OPEN(unit = 243, file = trim(EXPCON_PATH)//'EXPCON/ttgg1.dat', status='old') !
DO iQ=1, 300
  READ(243,*,IOSTAT=ios) ttgg1(iQ,1), ttgg1(iQ,2)
  IF (ios /= 0) EXIT
  Nttgg1 = Nttgg1 +1
ENDDO
!!-------------------------------------------
Ntttt02 =  0 
OPEN(unit = 245, file = trim(EXPCON_PATH)//'EXPCON/tttt02.dat', status='old') !
DO iQ=1, 300
  READ(245,*,IOSTAT=ios) tttt02(iQ,1), tttt02(iQ,2)
  IF (ios /= 0) EXIT
  Ntttt02 = Ntttt02 + 1
ENDDO
!!-------------------------------------------------
Ntttt04 = 0 
OPEN(unit = 246, file = trim(EXPCON_PATH)//'EXPCON/tttt04.dat', status='old') !
DO iQ=1, 300
  READ(246,*,IOSTAT=ios) tttt04(iQ,1), tttt04(iQ,2)
  IF (ios /= 0) EXIT
  Ntttt04 = Ntttt04 + 1
ENDDO
!!-------------------------------------------------
Ntttt05 = 0
OPEN(unit = 247, file = trim(EXPCON_PATH)//'EXPCON/tttt05.dat', status='old') !
DO iQ=1, 300
  READ(247,*,IOSTAT=ios) tttt05(iQ,1), tttt05(iQ,2)
  IF (ios /= 0) EXIT
  Ntttt05 = Ntttt05 + 1
ENDDO
!!-------------------------------------------------
Ntttt06 = 0 
OPEN(unit = 248, file = trim(EXPCON_PATH)//'EXPCON/tttt06.dat', status='old') !
DO iQ=1, 300
  READ(248,*,IOSTAT=ios) tttt06(iQ,1), tttt06(iQ,2)
  IF (ios /= 0) EXIT
  Ntttt06 = Ntttt06 +1
ENDDO
!!-------------------------------------------------
Ntttt08 = 0 
OPEN(unit = 249, file = trim(EXPCON_PATH)//'EXPCON/tttt08.dat', status='old') !
DO iQ=1, 300
  READ(249,*,IOSTAT=ios) tttt08(iQ,1), tttt08(iQ,2)
  IF (ios /= 0) EXIT
  Ntttt08 = Ntttt08 +1
ENDDO
!!-------------------------------------------------
Ntttt1 = 0
OPEN(unit = 250, file = trim(EXPCON_PATH)//'EXPCON/tttt1.dat', status='old') !
DO iQ=1, 300
  READ(250,*,IOSTAT=ios) tttt1(iQ,1), tttt1(iQ,2)
  IF (ios /= 0) EXIT
  Ntttt1 = Ntttt1 +1
ENDDO
!!-------------------------------------------------------
NhA4b = 0 
OPEN(unit = 260, file = trim(EXPCON_PATH)//'EXPCON/hA4b.dat', status='old') !
DO iQ=1, 10000
  READ(260,*,IOSTAT=ios) hA4b(iQ,1), hA4b(iQ,2), hA4b(iQ,3)
  IF (ios /= 0) EXIT
  NhA4b = NhA4b +1
ENDDO
!!----------------------------------------------------------
NhA4tau = 0 
OPEN(unit = 261, file = trim(EXPCON_PATH)//'EXPCON/hA4tau.dat', status='old') !
DO iQ=1, 10000
  READ(261,*,IOSTAT=ios) hA4tau(iQ,1), hA4tau(iQ,2), hA4tau(iQ,3)
  IF (ios /= 0) EXIT
  NhA4tau = NhA4tau +1
ENDDO
!!----------------------------------------------------------
NhA2b2tau = 0 
OPEN(unit = 262, file = trim(EXPCON_PATH)//'EXPCON/hA2b2tau.dat', status='old') !
DO iQ=1, 10000
  READ(262,*,IOSTAT=ios) hA2b2tau(iQ,1), hA2b2tau(iQ,2), hA2b2tau(iQ,3)
  IF (ios /= 0) EXIT
  NhA2b2tau = NhA2b2tau +1
ENDDO
!!----------------------------------------------------------
NhA2tau2b = 0 
OPEN(unit = 263, file = trim(EXPCON_PATH)//'EXPCON/hA2tau2b.dat', status='old') !
DO iQ=1, 10000
  READ(263,*,IOSTAT=ios) hA2tau2b(iQ,1),hA2tau2b(iQ,2), hA2tau2b(iQ,3)
  IF (ios /= 0) EXIT
  NhA2tau2b = NhA2tau2b +1
ENDDO
!!----------------------------------------------------------
NAAA6b = 0 
OPEN(unit = 264, file = trim(EXPCON_PATH)//'EXPCON/AAA6b.dat', status='old') !
DO iQ=1, 10000
  READ(264,*,IOSTAT=ios) AAA6b(iQ,1), AAA6b(iQ,2), AAA6b(iQ,3)
  IF (ios /= 0) EXIT
  NAAA6b = NAAA6b +1
ENDDO
!!----------------------------------------------------------
NAAA6tau = 0 
OPEN(unit = 265, file = trim(EXPCON_PATH)//'EXPCON/AAA6tau.dat', status='old') !
DO iQ=1, 10000
  READ(265,*,IOSTAT=ios) AAA6tau(iQ,1), AAA6tau(iQ,2), AAA6tau(iQ,3)
  IF (ios /= 0) EXIT
  NAAA6tau = NAAA6tau +1
ENDDO

!!----------------------------------------------------------
Nstblsn = 0 
OPEN(unit = 275, file = trim(EXPCON_PATH)//'EXPCON/stblsn.dat', status='old') !
DO iQ=1, 10000
  READ(275,*,IOSTAT=ios) stblsn(iQ,1), stblsn(iQ,2)
  IF (ios /= 0) EXIT
  Nstblsn = Nstblsn+1
ENDDO
!!----------------------------------------------------------
Nstnc = 0 
OPEN(unit = 276, file = trim(EXPCON_PATH)//'EXPCON/stnc.dat', status='old') !
DO iQ=1, 10000
  READ(276,*,IOSTAT=ios) stnc(iQ,1), stnc(iQ,2)
  IF (ios /= 0) EXIT
  Nstnc = Nstnc+1
ENDDO
!!----------------------------------------------------------
Nsbnb = 0 
OPEN(unit = 277, file = trim(EXPCON_PATH)//'EXPCON/sbnb.dat', status='old') !
DO iQ=1, 10000
  READ(277,*,IOSTAT=ios) sbnb(iQ,1), sbnb(iQ,2)
  IF (ios /= 0) EXIT
  Nsbnb = Nsbnb+1
ENDDO

!!--------------------------------------------------------------
! Read ATLAS upper limit
! ggHgg(IQ,1): Higgs mass
! ggHgg(IQ,4): upper limit in fb
NggHgg = 0 
OPEN(unit = 284, file = trim(EXPCON_PATH)//'EXPCON/ggHgg.dat', status='old') !
DO iQ=1, 100
  READ(284,*,IOSTAT=ios) ggHgg(iQ,1), ggHgg(iQ,2), ggHgg(iQ,3), ggHgg(iQ,4)
  IF (ios /= 0) EXIT
  NggHgg = NggHgg +1
ENDDO

! Read SM Higgs ggF production cross section (60-122 GeV)
! SMXS(I,1): Higgs mass
! SMXS.65-122.dat:  ggF production cross section in pb
! SMXS(I,2): ggF production cross section in fb
NSMXS = 0 
OPEN(unit = 285, file = trim(EXPCON_PATH)//'EXPCON/SMXS.65-122.dat', status='old') !
DO iQ=1, 100
  READ(285,*,IOSTAT=ios) SMXS(iQ,1), SMXS(iQ,2), SMXS(iQ,3), SMXS(iQ,4), SMXS(iQ,5), SMXS(iQ,6)
  SMXS(iQ,2)=1.d3*SMXS(iQ,2)
  IF (ios /= 0) EXIT
  NSMXS = NSMXS +1
ENDDO
!!-------------------------------------------------------------

endif !end first

first = .false.
  
!!******************************************************
!!********** Starting the scan *************************
!!******************************************************



!WRITE(24,*) real(kapIN(1,1,1), dp), real(lamIN(1), dp)   !! DK on 29-02-16

!!**************************************************************
!!--- Starting constraints tests -------------------------------
!!**************************************************************

!!--------------------------------------------------------------
!!--------------- Test on charged slepton masses ---------------
!!--------------------------------------------------------------
iH = 0
iH1 = 0
iH125 = 0
iH126 = 0
iH127 = 0
Mhiggs = 0.d0
! PROB(55) = 5.d0
!!******************
DO iEL = 1, 8

   if(DABS(ZP(iEL,3)).GE. 0.71 ) then ! picking up eL~ 
    MEL = MHpm(iEL)
   endif
   if(DABS(ZP(iEL,4)).GE. 0.71 ) then ! picking up muL~ 
    MML = MHpm(iEL)
   endif
!!---------------------------------------------
   if(DABS(ZP(iEL,6)).GE. 0.71) then  ! picking up eR~
    MER = MHpm(iEL)
   endif
   if(DABS(ZP(iEL,7)).GE. 0.71) then ! picking up muR~ 
    MMR = MHpm(iEL)
   endif
!!-------------------------------------------------------------------
   if(DABS(ZP(iEL,5)).GE. 0.71) then ! picking up tauL~ 
    MTL = MHpm(iEL)
   endif
   if(DABS(ZP(iEL,8)).GE. 0.71) then ! picking up tauRL~ 
    MTR = MHpm(iEL)
   endif 
!!---------------------------------------------------------------------
   if(DABS(ZP(iEL,1)).GE. 0.71) then   ! picking up charged higgs 
    MCHd = MHpm(iEL) 
   endif
   if(DABS(ZP(iEL,2)).GE. 0.71) then 
    MCHu = MHpm(iEL)
   endif 

!!---------------------------------------------------------------------
!! Looking for the true ~nueL and ~numuL out of 8 scalars and 7 pseudo
 If(DABS(ZH(iEL,3)).GE. 0.71) then   ! > 50%
  MNL1 = Mhh(iEL) 
 Endif
 If(DABS(ZH(iEL,4)).GE. 0.71) then   ! > 50%
  MNL2 = Mhh(iEL)
 Endif 
 If(DABS(ZA(iEL,3)).GE. 0.71) then   ! > 50%
  MNL3 = MAh(iEL) 
 Endif
 If(DABS(ZA(iEL,4)).GE. 0.71) then   ! > 50%
  MNL4 = MAh(iEL)
 Endif

!!-------------- Neutral 125 GeV Higgs (122 < Mh < 128 GeV)-------------
 MHmin = 122d0
 MHmax = 128d0
 If( (Mhh(iEL).gt. MHmin) .and. (Mhh(iEL).lt. MHmax) ) THEN  
   iH125 = iH125 + 1
   WRITE(17,*)  ZH(iEL,1),' ',  ZH(iEL,2),' ', ZH(iEL,1)**2 + ZH(iEL,2)**2   !! DK
   WRITE(17,*) iEL, '  ',  Mhh(iEL),' ', BRhh(iEL,184)  !! DK
   WRITE(17,*)  '  ' !! DK

    PROBtest(1) = DDIM(1d0, (ZH(iEL,1)**2 + ZH(iEL,2)**2 )/ 7.5d-1) 
   If( PROBtest(1).eq. 0d0 ) THEN      !! 80% of Higgs composition must be a mixture of Hu & Hd 
     iH126 = iH126 +1
!      PROBtest(2) = DDIM(1d0, DABS(ZH(iEL,2))/ 8.9d-1)
!       If( PROBtest(2).eq. 0d0) THEN !! 80% of the true 125 GeV Higgs composition must be Hu 
         iH127 = iH127 +1             
         iH = iEL  
         Mhiggs = Mhh(iH)
!       Endif
   EndIf
 Endif


ENDDO
!!******************

If(iH125 .eq. 0) then
 PROB(53) = 5.d0
endif
If(iH126 .eq. 0) then
 PROB(55) = 5.d0
endif
If(iH127 .eq. 0) then
 PROB(52) = 5.d0
endif


!!-----------------------------------------------------------------------
!! ---- Filling histograms for Mhh plot in (kappa, lambda) plane --------
!!-----------------------------------------------------------------------

 !WRITE(21,*) real(kapIN(2,2,2),dp),real(lamIN(2),dp), Mhiggs
 !WRITE(22,*) real(kapIN(3,3,3),dp),real(lamIN(3),dp), Mhiggs
! WRITE(17,*)  Mhh(1),Mhh(2),Mhh(3),Mhh(4),Mhh(5),Mhh(6)
 !WRITE(23,*)  BRAh(iEL,4),gThh(1),gThh(4),BRhh(iEL,1),Real(ratioGG(2),dp)
 !WRITE(23,*)  iH125, iH126, iH127, iH, Mhiggs
 !WRITE(23,*)  gTSd(3),gTSu(3)
!!---------------------------------------------------
!!----- test on stau and charged sleptons ----------- 
!!---------------------------------------------------

 MRR = MIN(MER, MMR) ! choosing the MIN of R-handed Sleptons
 MLL = MIN(MEL, MML) ! choosing the MIN of L-handed Sleptons
 MTT = MIN(MTR, MTL) ! choosing the MIN of R&L handed taus
 
! PROB(1) = DDIM(1d0,MIN(MRR,MLL)/MSLMIN) 
! PROB(2) = DDIM(1d0,MTT/MSTMIN)
 
!!-------------------------------------------------------------------
!!----------- Test on charged Higgs mass (the one from Hu,Hd)  ------
!!-------------------------------------------------------------------
 CMASS = MAX(MCHu, MCHd) ! The min would concide with the Goldstone
 PROB(3) = DDIM(1d0,CMASS/78.6)
 
!!-----------------------------------------------------
!!----------- Test on chargino mass -------------------
!!-----------------------------------------------------
 PROB(4) = DDIM(1d0,DABS(MCha(4))/103.5)

!!-----------------------------------------------------
!!----------- Test on gluino/squark masses ------------
!!-----------------------------------------------------

DO iQ = 1, 6

  If(dabs(real(ZU(iQ,1))).GT. 0.71)  Then  ! looking for the correct uL~, > 50% of ~uL
    MUL = MSu(iQ)
  EndIf
  If(dabs(real(ZU(iQ,2))).GT. 0.71)  Then  ! looking for the correct cL~ 
    MUL = MIN(MSu(iQ),MUL)
  EndIf
  If(dabs(real(ZD(iQ,1))).GT. 0.71)  Then  ! looking for the correct  dL~
   MDL = MSd(iQ)
  EndIf
  If(dabs(real(ZD(iQ,2))).GT. 0.71)  Then  !  looking for the correct sL~
   MDL = MIN(MSd(iQ),MDL)
  EndIf
!!
  If(dabs(real(ZU(iQ,4))).GT. 0.71)  Then  !  looking for the correct uR~
   MUR = MSu(iQ)
  EndIf
  If(dabs(real(ZU(iQ,5))).GT. 0.71)  Then  !  looking for the correct cR~
    MUR = MIN(MSu(iQ),MUR)
  EndIf
  If(dabs(real(ZD(iQ,4))).GT. 0.71)  Then  !  looking for the correct dR~
   MDR = MSd(iQ)
  EndIf
  If(dabs(real(ZD(iQ,5))).GT. 0.71)  Then  !  looking for the correct sR~
   MDR = MIN(MSd(iQ),MDR)
  EndIf
!!
  If(dabs(real(ZD(iQ,3))).GT. 0.71)  Then  !  looking for the correct bL~
   MBL = MSd(iQ)
  EndIf
  If(dabs(real(ZD(iQ,6))).GT. 0.71)  Then  !  looking for the correct bR~
   MBR = MSd(iQ)
  EndIf
!!
  If(dabs(real(ZU(iQ,3))).GT. 0.71)  Then  !  looking for the correct tL~
   MTopL = MSu(iQ)
  EndIf
  If(dabs(real(ZU(iQ,6))).GT. 0.71)  Then  !  looking for the correct tR~
   MTopR = MSu(iQ)
  EndIf
ENDDO

!!*****************************************************************************
!! Precision of masses computations
!Write(20,*) '# -- Contains info about precisons on susy masses calculations --'
!DO iQ = 1, 6

!DeltaMSd(iQ) = (MSd_1L(iQ)-MSd(iQ))/MSd_1L(iQ)
!DeltaMSu(iQ) = (MSu_1L(iQ)-MSu(iQ))/MSu_1L(iQ)
!DeltaM_squarks = DeltaM_squarks + DeltaMSu(iQ) + DeltaMSd(iQ)
!Write(20,*) 'Squarks' , DeltaMSd(iQ), DeltaMSu(iQ)
!Write(20,*) 'Squarks' , DeltaMSd(iQ), DeltaMSu(iQ)
!ENDDO


!!-----------------------------------------------------------
!!-----------------------------------------------------------
!  WRITE(23,*) MCHAMIN,MTopR,MTopL,MIN(MTopR,MTopL)


!------------------------------------------------------------ 
!  PROB(5) = DDIM(1d0,MIN(MUR,MUL,MDR,MDL)/MSQMIN) +DDIM(1d0,DABS(MGlu)/MGLMIN)
!  PROB(6) = DDIM(1d0,MIN(MBR,MBL)/MSBMIN)
!  PROB(7) = DDIM(1d0,MIN(MTopR,MTopL)/MSTopMIN)

!!*********************************************************************************************
!!------------- Test on stop (lightest) -> b l sneutrino -- ! Switching them because these analysis are

!!------------------------------------------------------------
!!------------- Test on stop (lightest) -> b l sneutrino -- ! Switching them because these analysis are
! model dependens.
!!------------------------------------------------------------
! I=1
! DO WHILE(stblsn(I,1) .LE. MIN(MTopR,MTopL) .AND. I.LT.Nstblsn)
!  I=I+1
! ENDDO

! IF(I.GT.1 .AND. MIN(MTopR,MTopL).LT.stblsn(Nstblsn,1))THEN
!  MMIN = stblsn(I-1,2)+((MIN(MTopR,MTopL)- stblsn(I-1,1))/(stblsn(I,1)-stblsn(I-1,1)))*(stblsn(I,2)-stblsn(I-1,2))
!  PROB(60) = DDIM(1d0,MNL1/MMIN)+DDIM(1d0,MNL2/MMIN)+DDIM(1d0,MNL3/MMIN)+DDIM(1d0,MNL4/MMIN)
! ENDIF

!!-----------------------------------------------------------------
!!------ Test on stop (lightest) -> neutralino(lightest) + c ------
!!-----------------------------------------------------------------
! I=1
! DO WHILE(stnc(I,1).LE.DABS(MChi(4)) .AND. I.LT.Nstnc)
!  I=I+1
! ENDDO

! IF(I.GT.1 .AND. DABS(MChi(4)).LT.stnc(Nstnc,1))THEN
!  MMIN = stnc(I-1,2)+((DABS(MChi(4))-stnc(I-1,1))/(stnc(I,1)-stnc(I-1,1)))*(stnc(I,2)-stnc(I-1,2))
!  PROB(61)=DDIM(1d0, MIN(MTopR,MTopL)/MMIN)
! ENDIF

!!-----------------------------------------------------------------
!!------ Test on sbottom (lightest) -> neutralino(lightest) b -----
!!-----------------------------------------------------------------
! I=1
! DO WHILE(sbnb(I,1).LE.MIN(MBR,MBL).AND. I.LT.Nsbnb)
!  I=I+1
! ENDDO
!IF(I.GT.1 .AND. MIN(MBR,MBL) .LT.sbnb(Nsbnb,1))THEN
!  MMIN = sbnb(I-1,2)+((MIN(MBR,MBL)-sbnb(I-1,1))/(sbnb(I,1)-sbnb(I-1,1)))*(sbnb(I,2)-sbnb(I-1,2))
!  PROB(62)= DDIM(1d0,DABS(MChi(4))/MMIN)
! ENDIF

!!***************************************************************
!!---------------------------------------------------------
!!------ Checking for the correct Goldstone bosons --------
!!---------------------------------------------------------

!! At least 80% of the Goldstones comme from Hu and Hd
  PROB(8) = DDIM(1.d0, (ZA(1,1)**2 + ZA(1,2)**2 )/ 8.0d-1)+ DDIM(1.d0, (ZP(1,1)**2 + ZP(1,2)**2 )/ 8.0d-1)

!!---------------------------------------------------------------------------
!!----- Constraint on new Z decay in munuSSM (ref: arXiv:1403.3675v3) -------
!!---------------------------------------------------------------------------
! BR(A->tau tau) = BR(jQ,99)
! BR(A->b b)     = BR(jQ,175)
! BR(h->tau tau) = BR(iQ,108)
! BR(h->b b)     = BR(iQ,184)


!! Here the computations are done assuming that ~chi0 are vR-like and h/A are ~vc-like 

!! Start computing Z -> ~chi0 ~chi0 width
 OLnnz1 = (0.d0,0.d0)
 OLnnz2 = (0.d0,0.d0)
 OLnnz3 = (0.d0,0.d0)
 ORnnz1 = (0.d0,0.d0)
 ORnnz2 = (0.d0,0.d0)
 ORnnz3 = (0.d0,0.d0)
 OLnnz  = (0.d0,0.d0)
 ORnnz  = (0.d0,0.d0)
 newGAMMAZchichi = 0.d0
 FactorZchichi = (Real(g2,dp))**2 / (48*pi*MVZ**5)
 Do iZ = 4, 6    ! 7 if one adds ~B
  Do jZ = iZ, 6  ! 7 
   if(MChi(iZ)+ MChi(jZ).le. MVZ) then
     If(((real(UV(iZ,5)).gt. 7.1d-1).or.(real(UV(iZ,6)).gt. 7.1d-1).or.(real(UV(iZ,7)).gt. 7.1d-1)) &
    & .and.((real(UV(jZ,5)).gt. 7.1d-1).or.(real(UV(jZ,6)).gt. 7.1d-1).or.(real(UV(jZ,7)).gt. 7.1d-1)) ) Then
!! The last "IF condition" checks the vR-like of neutralinos

       OLnnz1 = OLnnz1 + ( UV(iZ,3)*DCONJG(UV(jZ,3)) - UV(iZ,4)*DCONJG(UV(jZ,4)) &
    &  + UV(iZ,8)*DCONJG(UV(jZ,8)) )/(2*sqrt(0.769)) 

       OLnnz2 = OLnnz2 + ( UV(iZ,3)*DCONJG(UV(jZ,3)) - UV(iZ,4)*DCONJG(UV(jZ,4)) &
    &  + UV(iZ,9)*DCONJG(UV(jZ,9)) )/(2*sqrt(0.769))

       OLnnz3 = OLnnz3 + ( UV(iZ,3)*DCONJG(UV(jZ,3)) - UV(iZ,4)*DCONJG(UV(jZ,4)) &
    &  + UV(iZ,10)*DCONJG(UV(jZ,10)))/(2*sqrt(0.769))
!
      ORnnz1 = ORnnz1 + (- UV(jZ,3)*DCONJG(UV(iZ,3)) + UV(jZ,4)*DCONJG(UV(iZ,4)) &
    &  - UV(jZ,8)*DCONJG(UV(iZ,8)) )/(2*sqrt(0.769))

       ORnnz2 = ORnnz2 + (- UV(jZ,3)*DCONJG(UV(iZ,3)) + UV(jZ,4)*DCONJG(UV(iZ,4)) &
    &  - UV(jZ,9)*DCONJG(UV(iZ,9)) )/(2*sqrt(0.769))

       ORnnz3 = ORnnz3 + (- UV(jZ,3)*DCONJG(UV(iZ,3)) + UV(jZ,4)*DCONJG(UV(iZ,4)) &
    & - UV(jZ,10)*DCONJG(UV(iZ,10)) )/(2*sqrt(0.769))
!
    OLnnz = OLnnz + OLnnz1 + OLnnz2 + OLnnz3
    ORnnz = ORnnz + ORnnz1 + ORnnz2 + ORnnz3
!
    newGAMMAZchichi = newGAMMAZchichi + FactorZchichi*sqrt((MVZ**2 - MChi(jZ)**2 &
   &  - MChi(jZ)**2)**2 - 4*MChi(iZ)**2 * MChi(jZ)**2)* &
   & ( (OLnnz*DCONJG(OLnnz) + ORnnz*DCONJG(ORnnz))*(2*MVZ**4 &
   &  - MVZ**2 * (MChi(iZ)**2 + MChi(jZ)**2)-(MChi(iZ)**2 - MChi(jZ)**2)**2) &
   &  + 12*Real(DCONJG(OLnnz)*ORnnz)*MChi(jZ)*MChi(jZ)*MVZ**2)
   Endif
   else
    newGAMMAZchichi = 0.d0
   endif 
  Enddo
 Enddo

!! Start computing Z -> h A
 Ospz1 = (0.d0,0.d0)
 Ospz2 = (0.d0,0.d0)
 Ospz3 = (0.d0,0.d0)
 Ospz  = 0.d0
 newGAMMAZsp = 0.d0
 FactorZsp = (Real(g2,dp))**2 / (192*pi*0.769*MVZ**5)

 Do iZ = 1, 3
  Do jZ = 2, 4 
   if(Mhh(iZ)+ MAh(jZ).le. MVZ) then
     If(((real(ZH(iZ,6)).gt. 7.1d-1).or.(real(ZH(iZ,7)).gt. 7.1d-1).or.(real(ZH(iZ,8)).gt. 7.1d-1)) &
    & .and.((real(ZA(jZ,3)).gt. 7.1d-1).or.(real(ZA(jZ,4)).gt. 7.1d-1).or.(real(ZA(jZ,5)).gt. 7.1d-1)) ) Then
!! The last "IF condition" checks the ~vc-like of h/A

    Ospz1 = Ospz1 +  ZH(iZ,1)*ZA(jZ,1) - ZH(iZ,2)*ZA(jZ,2) + ZH(iZ,6)*ZA(jZ,6) ! check if ZA(3 or 6)
    Ospz2 = Ospz2 +  ZH(iZ,1)*ZA(jZ,1) - ZH(iZ,2)*ZA(jZ,2) + ZH(iZ,7)*ZA(jZ,7) ! check if ZA(3 or 7)
    Ospz3 = Ospz3 +  ZH(iZ,1)*ZA(jZ,1) - ZH(iZ,2)*ZA(jZ,2) + ZH(iZ,8)*ZA(jZ,8) ! check if ZA(3 or 8)
    Ospz  =  Ospz + Ospz1 + Ospz2 + Ospz3
!
    newGAMMAZsp =newGAMMAZsp+Ospz*DCONJG(Ospz)*FactorZsp*sqrt((MVZ**2-Mhh(iZ)**2 & 
    & -MAh(jZ)**2)**2)*(MVZ**4-2*(MVZ**2)*(Mhh(iZ)**2+MAh(jZ)**2)+(Mhh(iZ)**2 - MAh(jZ)**2)**2 )
    Endif
    else 
     newGAMMAZsp = 0.d0
   endif 
  Enddo
 Enddo
!

 newGAMMAZ = newGAMMAZsp + newGAMMAZchichi  !! total new Z decays width

!! Here  at 1sigma δGamma = 3.4 MeV can accomodate total new decay of Z (ref: arXiv:1403.3675v3). 
!! We use this value to constraint new Z decays  in munuSSM

 PROB(9) = DDIM(newGAMMAZ/3.4d-3, 1.d0)
!!*****************************************************************
!!-----------------------------------------------------------------
!!------- Constraint on the total Higgs width in the munuSSM ------
!!-----------------------------------------------------------------
 GammaHtotmax = 22.d-3
!rruiz: problem with gThh since iH can't be 0 due to array dim.
 IF(iH.GE.1) THEN
 PROB(10) = DDIM(gThh(iH)/22.d-3, 1.d0)
 ENDIF
!!-------------------------------------------------------------
!!*************************************************************

!!-----------------------------------------------------------------------------------
  hd = 1d0/DSQRT(2d0*DSQRT(2d0)*(1d0+Real(TanBeta,dp)**2)*G_F) ! = vd/sqrt(2)
  hu = hd*Real(TanBeta,dp)                                     ! = vu/sqrt(2)
!!******************************************** Test DK
 mBino = Real(M1,dp)
 mWino = Real(M2,dp)
 mvR(1)= 2* Real(kap(1,1,1),dp)*vR(1)
 mvR(2)= 2* Real(kap(2,2,2),dp)*vR(2)
 mvR(3)= 2* Real(kap(3,3,3),dp)*vR(3)

Do ip =1, 3
mneuBino(ip) = (vL(ip)*Real(g2,dp))**2 / mBino
mneuWino(ip) = (vL(ip)*Real(g2,dp))**2 / mWino
mneuvR(ip)   = (Real(Yv(ip,ip),dp)*hu)**2 / mvR(ip)
Enddo
  WRITE(23,*)  MChi(1), MChi(2), MChi(3)
  WRITE(23,*)  mBino, mWino
  WRITE(23,*)  mneuBino(1), mneuBino(2), mneuBino(3)
  WRITE(23,*)  mneuWino(1), mneuWino(2), mneuWino(3)
  WRITE(23,*)  mvR(1), mvR(2), mvR(3)
  WRITE(23,*)  mneuvR(1), mneuvR(2), mneuvR(3)
  WRITE(23,*)  ' '
  WRITE(23,*)  ' '
!!******************************************** Test DK

!!-------------------------------------------------------------------
!!------- Constraints on new scalars decays in munuSSM from LEP -----
!!-------------------------------------------------------------------
! BR(A->tau tau) = BR(jQ,99)
! BR(A->b b)     = BR(jQ,175)
! BR(h->tau tau) = BR(iQ,108)
! BR(h->b b)     = BR(iQ,184)

  vtotSM2 = hd**2 + hu**2 + vL(1)**2 + vL(2)**2 + vL(3)**2

DO iQ= 1, 8  !! loop over the 8 scalars 

  vmunu2  = ZH(iQ,1)*hd + ZH(iQ,2)*hu + vL(1)*ZH(iQ,6) + vL(2)*ZH(iQ,7) + vL(3)*ZH(iQ,8)
 
!rruiz: mistake here iQ has to be larger or equal than 1
  IF(iQ.GT.1 .AND. DABS(Mhh(iQ)-Mhh(iQ-1)).LT.3d0)THEN 
   tempMH = Mhh(iQ-1)
   Mhh(iQ) = tempMH
   R = R + vmunu2**2/vtotSM2
  ELSE
   R = vmunu2**2 / vtotSM2
  ENDIF

!!------- ee-> Z h independent flavour ------------------------------
 IF(Mhh(iQ)+MVZ.LT. 209.d0)THEN

 Rmax=1d0
 J=1
  DO WHILE(hZind(J,1).LE.Mhh(iQ) .AND. J.LT.NhZind)
     J=J+1
  ENDDO
  IF(J.GT.1 .AND. Mhh(iQ).LT.hZind(NhZind,1)) THEN
    Rmax = hZind(J-1,2)+ ((Mhh(iQ)-hZind(J-1,1))/(hZind(J,1) &
     &  -hZind(J-1,1)))*(hZind(J,2)-hZind(J-1,2))
  ENDIF 
PROB(11) = PROB(11) + DDIM(R/Rmax, 1.d0)

!!-------ee-> Z h -> Z b b ------------------------------------------
  Rmax=1d0
  J=1
  DO WHILE(hZbb(J,1).LE.Mhh(iQ) .AND. J.LT.NhZbb)
  J=J+1
  ENDDO
 IF(J.GT.1 .AND. Mhh(iQ).LT.hZbb(NhZbb,1)) THEN
   Rmax = hZbb(J-1,2)+((Mhh(iQ)-hZbb(J-1,1))/(hZbb(J,1)-hZbb(J-1,1))) &
  &   *(hZbb(J,2)-hZbb(J-1,2))
  ENDIF
  PROB(12) = PROB(12) + DDIM(R*BRhh(iQ,184)/Rmax, 1.d0)

!!  WRITE(17,*)  R, Rmax, BRhh(iQ,184)

!!-------ee-> Z h -> Z j j ------------------------------------------
 Rmax=1d0
 J=1
 DO WHILE(hZjj(J,1).LE.Mhh(iQ) .AND. J.LT.NhZjj)
 J=J+1
 ENDDO

 IF(J.GT.1 .AND. Mhh(iQ).LT.hZjj(NhZjj,1)) THEN
  Rmax = hZjj(J-1,2)+((Mhh(iQ)-hZjj(J-1,1))/(hZjj(J,1)-hZjj(J-1,1))) &
  &  *(hZjj(J,2)-hZjj(J-1,2))
 ENDIF
 PROB(13) = PROB(13) + DDIM(R*(BRhh(iQ,2)+BRhh(iQ,180)+BRhh(iQ,184)+ BRhh(iQ,189))/Rmax, 1.d0)

!
!!------- ee-> Z h -> Z tau tau -------------------------------------
  Rmax=1d0
  J=1
  DO WHILE(hZll(J,1).LE.Mhh(iQ) .AND. J.LT.NhZll)
  J=J+1
  ENDDO
  IF(J.GT.1 .AND. Mhh(iQ).LT.hZll(NhZll,1)) THEN
   Rmax=hZll(J-1,2)+((Mhh(iQ)-hZll(J-1,1))/(hZll(J,1)-hZll(J-1,1))) &
   &     *(hZll(J,2)-hZll(J-1,2))
  ENDIF
  PROB(14) = PROB(14) + DDIM(R*BRhh(iQ,108)/Rmax, 1.d0)

!!-------  ee-> Z h -> Z Gammma Gamma ------------------------
  Rmax=1d0
  J=1
  DO WHILE(hZgg(J,1).LE.Mhh(iQ) .AND. J.LT.NhZgg)
  J=J+1
  ENDDO

  IF(J.GT.1 .AND. Mhh(iQ).LT.hZgg(NhZgg,1)) THEN
   Rmax=hZgg(J-1,2)+ ((Mhh(iQ)-hZgg(J-1,1))/(hZgg(J,1)-hZgg(J-1,1))) &
   &   *(hZgg(J,2)-hZgg(J-1,2))
  ENDIF
  PROB(15) = PROB(15) + DDIM(R*BRhh(iQ,1)/Rmax, 1.d0)

!!-------------------------------------------------------
!!----------------  ee -> Z h, h-> AA/hh -> 4bs ---------
!!-------------------------------------------------------

  IF(Mhh(iQ).GT. 2d0*Mhh(1))THEN !! h -> h1h1 -> 4bs
    Rmax=1d0
    DO kQ=1,NAAZ4b
! rruiz: it looks like KB is 10001 larger than 10000 which is the max number
     IF(AINT(Mhh(iQ)).EQ.AAZ4b(kQ,1) .AND. AINT(Mhh(1)).EQ.AAZ4b(kQ,2))THEN
     Rmax=AAZ4b(kQ,3)
     GOTO 15
     ENDIF
    ENDDO
15  PROB(21)=PROB(21)+DDIM(R*BRHHH(iQ,1)*BRhh(1,184)**2/Rmax,1d0)
  ENDIF

 IF(Mhh(iQ).GT. 2d0*Mhh(1))THEN !! h -> h1 h1 -> 4taus
   Rmax=1d0
   DO kQ=1,NAAZ4tau
    IF(AINT(Mhh(iQ)).EQ.AAZ4tau(kQ,1) .AND.AINT(Mhh(1)).EQ.AAZ4tau(kQ,2))THEN
     Rmax=AAZ4tau(kQ,3)
     GOTO 35
    ENDIF
   ENDDO
35  PROB(22)=PROB(22)+DDIM(R*BRHHH(iQ,1)*BRhh(1,108)**2/Rmax,1d0)
 ENDIF

 IF(Mhh(iQ).GT. 2d0*Mhh(1)) THEN  !! h -> h1 h1
   Rmax=1d0
   DO kQ=1,NAAZ2b2tau
    IF(AINT(Mhh(iQ)).EQ.AAZ2b2tau(kQ,1) .AND. AINT(Mhh(1)).EQ.AAZ2b2tau(kQ,2)) THEN
    Rmax=AAZ2b2tau(kQ,3)
    GOTO 22
    ENDIF
   ENDDO
22   PROB(23)=PROB(23) + DDIM(R*BRHHH(iQ,1)*2d0*BRhh(1,184)*BRhh(1,108)/Rmax,1d0)
 ENDIF


  DO jQ = 2, 8

  IF(Mhh(iQ).GT. 2d0*MAh(jQ))THEN !! h -> AA -> 4bs
    Rmax=1d0
    DO kQ=1,NAAZ4b
    IF(AINT(Mhh(iQ)).EQ.AAZ4b(kQ,1) .AND. AINT(MAh(jQ)).EQ.AAZ4b(kQ,2))THEN
      Rmax=AAZ4b(kQ,3)
      GOTO 10
    ENDIF
    ENDDO
10   PROB(21)=PROB(21)+DDIM(R*BRHAA(iQ,jQ)*BRAh(jQ,175)**2/Rmax,1d0)
!    Write(17,*) '1- BRAh(',jQ,'175)=',BRAh(jQ,175), 'BRHAA(',iQ,jQ,')=',BRHAA(iQ,jQ)
  ENDIF

  IF(Mhh(iQ).GT. 2d0*Mhh(jQ))THEN !! h -> hh -> 4bs
    Rmax=1d0
    DO kQ=1,NAAZ4b
    IF(AINT(Mhh(iQ)).EQ.AAZ4b(kQ,1) .AND. AINT(Mhh(jQ)).EQ.AAZ4b(kQ,2))THEN
     Rmax=AAZ4b(kQ,3)
     GOTO 20
    ENDIF
    ENDDO
20  PROB(21)=PROB(21)+DDIM(R*BRHHH(iQ,jQ)*BRhh(jQ,184)**2/Rmax,1d0)

  ENDIF

!!---------------------------------------------------------
!!----------------- ee -> hZ, h -> AA/hh -> 4taus  --------
!!---------------------------------------------------------
 IF(Mhh(iQ).GT. 2d0*MAh(jQ))THEN  !! h -> AA -> 4taus
   Rmax=1d0
   DO kQ=1,NAAZ4tau
    IF(AINT(Mhh(iQ)).EQ.AAZ4tau(kQ,1) .AND.AINT(MAh(jQ)).EQ.AAZ4tau(kQ,2))THEN
     Rmax=AAZ4tau(kQ,3)
     GOTO 30
    ENDIF
   ENDDO
!!       Apply only for MAh < 9.4GeV where A <-> eta_b mixing is absent
30 IF(MAh(jQ).LE. 9.4d0) THEN
    PROB(22)=PROB(22)+DDIM(R*BRHAA(iQ,jQ)*BRAh(jQ,99)**2/Rmax,1d0)
   ENDIF
 ENDIF      

 IF(Mhh(iQ).GT. 2d0*Mhh(jQ))THEN !! h -> h h -> 4taus
   Rmax=1d0
   DO kQ=1,NAAZ4tau
    IF(AINT(Mhh(iQ)).EQ.AAZ4tau(kQ,1) .AND.AINT(Mhh(jQ)).EQ.AAZ4tau(kQ,2))THEN
    Rmax=AAZ4tau(kQ,3)
    GOTO 40
    ENDIF
   ENDDO
40  PROB(22)=PROB(22)+DDIM(R*BRHHH(iQ,jQ)*BRhh(jQ,108)**2/Rmax,1d0)
!    Write(17,*) '4- BRHHH(',iQ,jQ,')=',BRHHH(iQ,jQ), 'BRhh(',jQ,'108)=',BRhh(jQ,108)
 ENDIF

!!----------------------------------------------------------------
!!----- ee -> hZ, h -> AA/hh -> 4taus (new ALEPH analysis) -------
!!----------------------------------------------------------------
       RMAX=1d0
       IF(Mhh(iQ).GE.70d0)THEN
        IF(MAh(jQ).GE.4d0.AND.MAh(jQ).LT.6d0)THEN
             DMA=(6d0-MAh(jQ))/2d0
             K=0
2001     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2001
         ENDIF
        ELSEIF(MAh(jQ).GE.6d0.AND.MAh(jQ).LT.8d0)THEN
         DMA=(8d0-MAh(jQ))/2d0
              K=0
2002     K=K+1
         HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2002
         ENDIF
        ELSEIF(MAh(jQ).GE.8d0.AND.MAh(jQ).LT.10d0)THEN
              DMA=(10d0-MAh(jQ))/2d0
              K=0
2003     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2003
         ENDIF
        ELSEIF(MAh(jQ).GE.10d0.AND.MAh(jQ).LE.12d0)THEN
              DMA=(12d0-MAh(jQ))/2d0
         K=0
2004     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2004
         ENDIF
        ENDIF
       ENDIF

!! Apply only for MA < 9.4GeV where A <-> eta_b mixing is absent:
       IF(MAh(jQ).LE.9.4d0) PROB(41)= PROB(41)+DDIM(R*BRHAA(iQ,jQ)*BRAh(jQ,99)**2/Rmax,1d0)

       RMAX=1d0
       IF(Mhh(iQ).GE.70d0)THEN
        IF(Mhh(1).GE.4d0.AND.Mhh(1).LT.6d0)THEN
             DMA=(6d0-Mhh(1))/2d0
             K=0
2005     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2005
         ENDIF
        ELSEIF(Mhh(1).GE.6d0.AND.Mhh(1).LT.8d0)THEN
             DMA=(8d0-Mhh(jQ-1))/2d0
             K=0
2006     K=K+1
             HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2006
         ENDIF
        ELSEIF(Mhh(1).GE.8d0.AND.Mhh(1).LT.10d0)THEN
             DMA=(10d0-Mhh(jQ))/2d0
             K=0
2007     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2007
         ENDIF
        ELSEIF(Mhh(1).GE.10d0.AND.Mhh(1).LE.12d0)THEN
             DMA=(12d0-Mhh(jQ))/2d0
             K=0
2008     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2008
         ENDIF
        ENDIF
       ENDIF
       PROB(41)= PROB(41)+DDIM(R*BRHHH(iQ,1)*BRhh(1,108)**2/Rmax,1d0)
!! 
       RMAX=1d0
       IF(Mhh(iQ).GE.70d0)THEN
        IF(Mhh(jQ).GE.4d0.AND.Mhh(jQ).LT.6d0)THEN
             DMA=(6d0-Mhh(jQ))/2d0
             K=0
2009     K=K+1
             HMIN=HM2(K)+(HM1(K)-HM2(K))*DMA
             HMAX=HM2(K+1)+(HM1(K+1)-HM2(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2009
         ENDIF
        ELSEIF(Mhh(jQ).GE.6d0.AND.Mhh(jQ).LT.8d0)THEN
             DMA=(8d0-Mhh(jQ-1))/2d0
             K=0
2010     K=K+1
             HMIN=HM3(K)+(HM2(K)-HM3(K))*DMA
             HMAX=HM3(K+1)+(HM2(K+1)-HM3(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2010
         ENDIF
        ELSEIF(Mhh(jQ).GE.8d0.AND.Mhh(jQ).LT.10d0)THEN
             DMA=(10d0-Mhh(jQ))/2d0
             K=0
2011     K=K+1
             HMIN=HM4(K)+(HM3(K)-HM4(K))*DMA
             HMAX=HM4(K+1)+(HM3(K+1)-HM4(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2011
         ENDIF
        ELSEIF(Mhh(jQ).GE.10d0.AND.Mhh(jQ).LE.12d0)THEN
             DMA=(12d0-Mhh(jQ))/2d0
             K=0
2012     K=K+1
             HMIN=HM5(K)+(HM4(K)-HM5(K))*DMA
             HMAX=HM5(K+1)+(HM4(K+1)-HM5(K+1))*DMA
             IF(Mhh(iQ).LE.HMAX)THEN
              RMAX=XIN(K)+(XIN(K+1)-XIN(K))*(Mhh(iQ)-HMIN)/(HMAX-HMIN)
             ELSEIF(K.LT.5)THEN
          GOTO 2012
         ENDIF
        ENDIF
       ENDIF
       PROB(41)=PROB(41)+DDIM(R*BRHHH(iQ,jQ)*BRhh(jQ,108)**2/Rmax,1d0)

!!------------------------------------------------------
!!----------- ee -> hZ, h -> AA/hh -> 2bs 2taus --------
!!------------------------------------------------------

 IF(Mhh(iQ).GT.2d0*MAh(jQ))THEN
   Rmax=1d0
   DO kQ =1,NAAZ2b2tau
    IF(AINT(Mhh(iQ)).EQ.AAZ2b2tau(kQ,1) .AND. AINT(MAh(jQ)).EQ.AAZ2b2tau(kQ,2))THEN
    Rmax=AAZ2b2tau(kQ,3)
    GOTO 21
    ENDIF
   ENDDO
!! Apply only for MAh < 9.4GeV or MAh > 10.5GeV without A <-> eta_b mixing
21  IF(MAh(jQ).LE. 9.4d0 .OR. MAh(jQ).GE. 10.5d0) THEN
      PROB(23)= PROB(23) + DDIM(R*BRHAA(iQ,jQ)*2d0*BRAh(jQ,99)*BRAh(jQ,175)/Rmax,1d0)
    ENDIF
 ENDIF

 IF(Mhh(iQ).GT. 2d0*Mhh(jQ)) THEN  !! h -> h h
   Rmax=1d0
   DO kQ=1,NAAZ2b2tau
    IF(AINT(Mhh(iQ)).EQ.AAZ2b2tau(kQ,1) .AND. AINT(Mhh(jQ)).EQ.AAZ2b2tau(kQ,2)) THEN
    Rmax=AAZ2b2tau(kQ,3)
    GOTO 23
    ENDIF
   ENDDO
23   PROB(23)=PROB(23) + DDIM(R*BRHHH(iQ,jQ)*2d0*BRhh(jQ,184)*BRhh(jQ,108)/Rmax,1d0)
 ENDIF

!!-----------------------------------------------------------
!!----------  ee -> Z h -> Z AA -> Z + light pairs   --------
!!-----------------------------------------------------------

 IF((Mhh(iQ).LT.2d0*MAh(jQ)).OR.(Mhh(iQ).LT.40d0).OR.(Mhh(iQ).GT.90d0).OR.(MAh(jQ).LT.2d0).OR.(MAh(jQ).GT.12d0))GOTO 752

!!-------------------  AA -> cccc
  ceff=R*BRHAA(iQ,jQ)*BRAh(jQ,180)**2
  IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 102
  IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 202
  IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 302
  IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 402
  IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 502
  IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 602
  GOTO 702

102  Rmax=.2d0
     kQ=1
  DO WHILE(cccc02(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Ncccc02)
    kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT. cccc02(Ncccc02,1))  THEN
    MAtest=cccc02(kQ-1,2)+((Mhh(iQ)-cccc02(kQ-1,1))/ (cccc02(kQ,1)-cccc02(kQ-1,1))) &
 &  *(cccc02(kQ,2)-cccc02(kQ-1,2))
  ENDIF
   PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
   GOTO 702

202  Rmax=.4d0
     kQ=1
  DO WHILE(cccc04(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Ncccc04)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1.AND.Mhh(iQ).LT.cccc04(Ncccc04,1)) THEN
   MAtest=cccc04(kQ-1,2)+((Mhh(iQ)-cccc04(kQ-1,1))/(cccc04(kQ,1)-cccc04(kQ-1,1))) &
 &  *(cccc04(kQ,2)-cccc04(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 702

302  Rmax=.5d0
     kQ=1
  DO WHILE(cccc05(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Ncccc05)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.cccc05(Ncccc05,1)) THEN
   MAtest=cccc05(kQ-1,2)+((Mhh(iQ)-cccc05(kQ-1,1))/(cccc05(kQ,1)-cccc05(kQ-1,1))) &
 &  *(cccc05(kQ,2)-cccc05(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 702

402  Rmax=.6d0
     kQ=1
  DO WHILE(cccc06(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Ncccc06)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.cccc06(Ncccc06,1)) THEN
   MAtest=cccc06(kQ-1,2)+((Mhh(iQ)-cccc06(kQ-1,1))/(cccc06(kQ,1)-cccc06(kQ-1,1))) &
 &  *(cccc06(kQ,2)-cccc06(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 702

502  Rmax=.8d0
     kQ=1
  DO WHILE(cccc08(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Ncccc08)
   kQ=kQ+1
  ENDDO
   MAtest=0d0
   IF(kQ.GT.1 .AND. Mhh(iQ).LT.cccc08(Ncccc08,1)) THEN
    MAtest=cccc08(kQ-1,2)+((Mhh(iQ)-cccc08(kQ-1,1))/(cccc08(kQ,1)-cccc08(kQ-1,1))) &
 &   *(cccc08(kQ,2)-cccc08(kQ-1,2))
   ENDIF
   PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
   GOTO 702

602  Rmax=1d0
     kQ=1
  DO WHILE(cccc1(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Ncccc1)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.cccc1(Ncccc1,1)) THEN
   MAtest=cccc1(kQ-1,2)+(Mhh(iQ)-cccc1(kQ-1,1))/(cccc1(kQ,1)-cccc1(kQ-1,1)) &
 &  *(cccc1(kQ,2)-cccc1(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
      
702 CONTINUE

!--------------------  AA -> ccjj 
  ceff=R*BRHAA(iQ,jQ)*BRAh(jQ,180)*(BRAh(jQ,2)+BRAh(jQ,171))*2d0
  IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 112
  IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 212
  IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 312
  IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 412
  IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 512
  IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 612
  GOTO 712

112  Rmax=.2d0
     kQ=1
  DO WHILE(ccgg02(kQ,1).LE.Mhh(iQ).AND. kQ.LT.Nccgg02)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.ccgg02(Nccgg02,1)) THEN
   MAtest = ccgg02(kQ-1,2)+((Mhh(iQ)-ccgg02(kQ-1,1))/(ccgg02(kQ,1)-ccgg02(kQ-1,1)) )&
 &   *(ccgg02(kQ,2)-ccgg02(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 712

212   Rmax=.4d0
       kQ=1
  DO WHILE(ccgg04(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Nccgg04)
    kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.ccgg04(Nccgg04,1)) THEN
   MAtest = ccgg04(kQ-1,2)+((Mhh(iQ)-ccgg04(kQ-1,1))/(ccgg04(kQ,1)-ccgg04(kQ-1,1))) &
 &   *(ccgg04(kQ,2)-ccgg04(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 712

312   Rmax=.5d0
       kQ=1
  DO WHILE(ccgg05(kQ,1).LE.Mhh(iQ).AND. kQ.LT.Nccgg05)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.ccgg05(Nccgg05,1)) THEN
   MAtest = ccgg05(kQ-1,2)+((Mhh(iQ)-ccgg05(kQ-1,1))/(ccgg05(kQ,1)-ccgg05(kQ-1,1))) &
 &  *(ccgg05(KQ,2)-ccgg05(KQ-1,2)) 
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 712

412   Rmax=.6d0
      kQ=1
  DO WHILE(ccgg06(kQ,1).LE.Mhh(iQ).AND. kQ.LT.Nccgg06)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.ccgg06(Nccgg06,1)) THEN 
   MAtest = ccgg06(kQ-1,2)+((Mhh(iQ)-ccgg06(kQ-1,1))/ (ccgg06(kQ,1)-ccgg06(kQ-1,1))) &
 &   *(ccgg06(kQ,2)-ccgg06(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 712

512   Rmax=.8d0
      kQ=1
  DO WHILE(ccgg08(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Nccgg08)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.ccgg08(Nccgg08,1)) THEN
   MAtest = ccgg08(kQ-1,2)+((Mhh(iQ)-ccgg08(kQ-1,1))/ (ccgg08(kQ,1)-ccgg08(kQ-1,1))) &
 &  *(ccgg08(kQ,2)-ccgg08(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 712

612   Rmax=1d0
      kQ=1
  DO WHILE(ccgg1(kQ,1).LE.Mhh(iQ) .AND. kQ.LT.Nccgg1)
   kQ=kQ+1
  ENDDO
  MAtest=0d0
  IF(kQ.GT.1 .AND. Mhh(iQ).LT.ccgg1(Nccgg1,1)) THEN
  MAtest=ccgg1(kQ-1,2)+((Mhh(iQ)-ccgg1(kQ-1,1))/ (ccgg1(kQ,1)-ccgg1(kQ-1,1))) &
 &  *(ccgg1(kQ,2)-ccgg1(kQ-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
    
712   CONTINUE

!----------------------  AA -> cctautau  
  ceff=R*BRHAA(iQ,jQ)*BRAh(jQ, 180)*BRAh(jQ,99)*2d0
  IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 122
  IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 222
  IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 322
  IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 422
  IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 522
  IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 622
  GOTO 722

122   Rmax=.2d0
      J=1
  DO WHILE(cctt02(J,1).LE.Mhh(iQ) .AND. J.LT.Ncctt02)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.cctt02(Ncctt02,1)) THEN
   MAtest=cctt02(J-1,2)+((Mhh(iQ)-cctt02(J-1,1))/(cctt02(J,1)-cctt02(J-1,1))) &
 &      *(cctt02(J,2)-cctt02(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 722

222   Rmax=.4d0
      J=1
  DO WHILE(cctt04(J,1).LE.Mhh(iQ) .AND. J.LT.Ncctt04)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.cctt04(Ncctt04,1)) THEN
  MAtest=cctt04(J-1,2)+((Mhh(iQ)-cctt04(J-1,1))/(cctt04(J,1)-cctt04(J-1,1))) &
 &     *(cctt04(J,2)-cctt04(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 722

322   Rmax=.5d0
      J=1
  DO WHILE(cctt05(J,1).LE.Mhh(iQ) .AND. J.LT.Ncctt05)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.cctt05(Ncctt05,1)) THEN
   MAtest=cctt05(J-1,2)+((Mhh(iQ)-cctt05(J-1,1))/(cctt05(J,1)-cctt05(J-1,1))) &
 &      *(cctt05(J,2)-cctt05(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 722

422   Rmax=.6d0
      J=1
  DO WHILE(cctt06(J,1).LE.Mhh(iQ) .AND. J.LT.Ncctt06)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.cctt06(Ncctt06,1)) THEN
   MAtest=cctt06(J-1,2)+((Mhh(iQ)-cctt06(J-1,1))/(cctt06(J,1)-cctt06(J-1,1))) &
&     *(cctt06(J,2)-cctt06(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 722

522   Rmax=.8d0
      J=1
  DO WHILE(cctt08(J,1).LE.Mhh(iQ).AND. J.LT.Ncctt08)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.cctt08(Ncctt08,1)) THEN
   MAtest=cctt08(J-1,2)+((Mhh(iQ)-cctt08(J-1,1))/(cctt08(J,1)-cctt08(J-1,1))) &
&     *(cctt08(J,2)-cctt08(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 722

622   Rmax=1d0
      J=1
  DO WHILE(cctt1(J,1).LE.Mhh(iQ) .AND. J.LT.Ncctt1)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.cctt1(Ncctt1,1)) THEN
   MAtest=cctt1(J-1,2)+((Mhh(iQ)-cctt1(J-1,1))/(cctt1(J,1)-cctt1(J-1,1))) &
&      *(cctt1(J,2)-cctt1(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
    
 722   CONTINUE

!!-----------------------  AA -> jjjj
  ceff=R*BRHAA(iQ,jQ)*(BRAh(jQ,2)+BRAh(jQ,171))**2
  IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 132
  IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 232
  IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 332
  IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 432
  IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 532
  IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 632
  GOTO 732

132   Rmax=.2d0
      J=1
  DO WHILE(gggg02(J,1).LE.Mhh(iQ) .AND. J.LT.Ngggg02)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.gggg02(Ngggg02,1)) THEN
   MAtest=gggg02(J-1,2)+((Mhh(iQ)-gggg02(J-1,1))/(gggg02(J,1)-gggg02(J-1,1))) &
&     *(gggg02(J,2)-gggg02(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 732

232   Rmax=.4d0
     J=1
  DO WHILE(gggg04(J,1).LE.Mhh(iQ) .AND. J.LT.Ngggg04)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.gggg04(Ngggg04,1)) THEN
   MAtest=gggg04(J-1,2)+((Mhh(iQ)-gggg04(J-1,1))/(gggg04(J,1)-gggg04(J-1,1))) &
 &     *(gggg04(J,2)-gggg04(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 732

332   Rmax=.5d0
     J=1
  DO WHILE(gggg05(J,1).LE.Mhh(iQ) .AND. J.LT.Ngggg05) 
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.gggg05(Ngggg05,1)) THEN
   MAtest=gggg05(J-1,2)+((Mhh(iQ)-gggg05(J-1,1))/ (gggg05(J,1)-gggg05(J-1,1))) &
&      *(gggg05(J,2)-gggg05(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 732

432   Rmax=.6d0
    J=1
  DO WHILE(gggg06(J,1).LE.Mhh(iQ) .AND. J.LT.Ngggg06)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.gggg06(Ngggg06,1)) THEN
   MAtest=gggg06(J-1,2)+((Mhh(iQ)-gggg06(J-1,1))/(gggg06(J,1)-gggg06(J-1,1))) &
&      *(gggg06(J,2)-gggg06(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 732

532   Rmax=.8d0
    J=1
  DO WHILE(gggg08(J,1).LE.Mhh(iQ).AND. J.LT.Ngggg08)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.gggg08(Ngggg08,1)) THEN
   MAtest=gggg08(J-1,2)+((Mhh(iQ)-gggg08(J-1,1))/(gggg08(J,1)-gggg08(J-1,1))) &
&     *(gggg08(J,2)-gggg08(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 732

632   Rmax=1d0
    J=1
  DO WHILE(gggg1(J,1).LE.Mhh(iQ) .AND. J.LT.Ngggg1)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.gggg1(Ngggg1,1)) THEN
   MAtest=gggg1(J-1,2)+((Mhh(iQ)-gggg1(J-1,1))/(gggg1(J,1)-gggg1(J-1,1))) &
&     *(gggg1(J,2)-gggg1(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
      
732   CONTINUE

!!-------------------  AA -> tautau jj
   ceff=R*BRHAA(iQ,jQ)*BRAh(jQ,99)*(BRAh(jQ,2)+BRAh(jQ,171))*2d0
  IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 142
  IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 242
  IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 342
  IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 442
  IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 542
  IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 642
  GOTO 742

142   Rmax=.2d0
     J=1
  DO WHILE(ttgg02(J,1).LE.Mhh(jQ).AND. J.LT.Nttgg02)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(jQ).LT.ttgg02(Nttgg02,1)) THEN
   MAtest=ttgg02(J-1,2)+((Mhh(jQ)-ttgg02(J-1,1))/(ttgg02(J,1)-ttgg02(J-1,1))) & 
&     *(ttgg02(J,2)-ttgg02(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 742

242   Rmax=.4d0
     J=1
  DO WHILE(ttgg04(J,1).LE.Mhh(iQ).AND. J.LT.Nttgg04)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.ttgg04(Nttgg04,1)) THEN
   MAtest=ttgg04(J-1,2)+((Mhh(iQ)-ttgg04(J-1,1))/(ttgg04(J,1)-ttgg04(J-1,1)) )&
&     *(ttgg04(J,2)-ttgg04(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 742

342   Rmax=.5d0
    J=1
  DO WHILE(ttgg05(J,1).LE.Mhh(iQ) .AND. J.LT.Nttgg05)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.ttgg05(Nttgg05,1)) THEN
   MAtest=ttgg05(J-1,2)+((Mhh(iQ)-ttgg05(J-1,1))/ (ttgg05(J,1)-ttgg05(J-1,1))) &
&     *(ttgg05(J,2)-ttgg05(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 742

442   Rmax=.6d0
    J=1
  DO WHILE(ttgg06(J,1).LE.Mhh(iQ).AND. J.LT.Nttgg06)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.ttgg06(Nttgg06,1)) THEN
   MAtest=ttgg06(J-1,2)+((Mhh(iQ)-ttgg06(J-1,1))/ (ttgg06(J,1)-ttgg06(J-1,1))) &
&     *(ttgg06(J,2)-ttgg06(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 742

542   Rmax=.8d0
   J=1
  DO WHILE(ttgg08(J,1).LE.Mhh(iQ) .AND. J.LT.Nttgg08)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.ttgg08(Nttgg08,1)) THEN
   MAtest=ttgg08(J-1,2)+((Mhh(iQ)-ttgg08(J-1,1))/(ttgg08(J,1)-ttgg08(J-1,1))) &
&     *(ttgg08(J,2)-ttgg08(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 742

642   Rmax=1d0
    J=1
  DO WHILE(ttgg1(J,1).LE.Mhh(iQ) .AND. J.LT.Nttgg1)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.ttgg1(Nttgg1,1)) THEN
   MAtest=ttgg1(J-1,2)+((Mhh(iQ)-ttgg1(J-1,1))/(ttgg1(J,1)-ttgg1(J-1,1))) &
&     *(ttgg1(J,2)-ttgg1(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
      
742   CONTINUE

!!----------------------   AA -> tautautautau
  ceff=R*BRHAA(iQ,jQ)*BRAh(jQ,99)**2
  IF(ceff.GE.0.0.AND.ceff.LT.0.2)GOTO 152
  IF(ceff.GE.0.2.AND.ceff.LT.0.4)GOTO 252
  IF(ceff.GE.0.4.AND.ceff.LT.0.5)GOTO 352
  IF(ceff.GE.0.5.AND.ceff.LT.0.6)GOTO 452
  IF(ceff.GE.0.6.AND.ceff.LT.0.8)GOTO 552
  IF(ceff.GE.0.8.AND.ceff.LT.1.0)GOTO 652
  GOTO 752

152   Rmax=.2d0
     J=1
  DO WHILE(tttt02(J,1).LE.Mhh(iQ).AND. J.LT.Ntttt02)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.tttt02(Ntttt02,1)) THEN
   MAtest=tttt02(J-1,2)+((Mhh(iQ)-tttt02(J-1,1))/(tttt02(J,1)-tttt02(J-1,1))) &
&    *(tttt02(J,2)-tttt02(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 752

252   Rmax=.4d0
    J=1
  DO WHILE(tttt04(J,1).LE.Mhh(iQ) .AND. J.LT.Ntttt04)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.tttt04(Ntttt04,1)) THEN
   MAtest=tttt04(J-1,2)+((Mhh(iQ)-tttt04(J-1,1))/ (tttt04(J,1)-tttt04(J-1,1))) &
&    *(tttt04(J,2)-tttt04(J-1,2))
  ENDIF
  PROB(19)=PROB(19)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 752

352   Rmax=.5d0
    J=1
  DO WHILE(tttt05(J,1).LE.Mhh(iQ).AND. J.LT.Ntttt05)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.tttt05(Ntttt05,1)) THEN
    MAtest=tttt05(J-1,2)+((Mhh(iQ)-tttt05(J-1,1))/(tttt05(J,1)-tttt05(J-1,1))) &
&      *(tttt05(J,2)-tttt05(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 752

452   Rmax=.6d0
    J=1
  DO WHILE(tttt06(J,1).LE.Mhh(iQ).AND. J.LT.Ntttt06)
  J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.tttt06(Ntttt06,1)) THEN
   MAtest=tttt06(J-1,2)+((Mhh(iQ)-tttt06(J-1,1))/ (tttt06(J,1)-tttt06(J-1,1))) &
&     *(tttt06(J,2)-tttt06(J-1,2))
  ENDIF   
  PROB(24)=PROB(24)+10d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 752

552   Rmax=.8d0
    J=1
  DO WHILE(tttt08(J,1).LE.Mhh(iQ) .AND. J.LT.Ntttt08)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.tttt08(Ntttt08,1)) THEN
   MAtest=tttt08(J-1,2)+((Mhh(iQ)-tttt08(J-1,1))/(tttt08(J,1)-tttt08(J-1,1))) &
&     *(tttt08(J,2)-tttt08(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
  GOTO 752

652   Rmax=1d0
   J=1
  DO WHILE(tttt1(J,1).LE.Mhh(iQ).AND. J.LT.Ntttt1)
   J=J+1
  ENDDO
  MAtest=0d0
  IF(J.GT.1 .AND. Mhh(iQ).LT.tttt1(Ntttt1,1)) THEN
   MAtest=tttt1(J-1,2)+((Mhh(iQ)-tttt1(J-1,1))/(tttt1(J,1)-tttt1(J-1,1))) &
&    *(tttt1(J,2)-tttt1(J-1,2))
  ENDIF
  PROB(24)=PROB(24)+5d0*(Rmax-ceff)*DDIM(MAtest/MAh(jQ),1d0)
    
752   CONTINUE

ENDDO !! jQ loop
ENDIF !! mh + mvz .lt. 209
ENDDO !! iQ loop

!!-----------------------------------------------------------
!!---------------- Associated production e e -> h A ---------
!!-----------------------------------------------------------

DO iQ = 1, 8
DO jQ = 2, 8

R =  (1+Real(TanBeta,dp)**2)**2*ZA(jQ,2)**2*(ZH(iQ,2)*hd - ZH(iQ,1)*hu)**2/(hd**2+hu**2+vL(1)**2+vL(2)**2+vL(3)**2)

!!------------  ee -> hA -> 4bs
 IF(Mhh(iQ)+MAh(jQ).LT. 209.d0 )THEN

  Massmin=MIN(Mhh(iQ),MAh(jQ))
  Massmax=MAX(Mhh(iQ),MAh(jQ))

  Rmax=1d0
  DO K=1,NhA4b
   IF(AINT(Massmax).EQ.hA4b(K,1) .AND. AINT(Massmin).EQ.hA4b(K,2))THEN
    Rmax = hA4b(K,3)
    GOTO 1003
   ENDIF
  ENDDO
1003     PROB(25)=PROB(25)+DDIM(R*BRhh(iQ,184)*BRAh(jQ,175)/Rmax,1d0)

!!----------  ee -> hA -> 4taus
  Rmax=1d0
  DO K=1,NhA4tau
   IF(AINT(Massmax).EQ.hA4tau(K,1) .AND. AINT(Massmin).EQ.hA4tau(K,2))THEN
    Rmax=hA4tau(K,3)
    GOTO 1004
    ENDIF
  ENDDO
1004     PROB(26)=PROB(26)+DDIM(R*BRhh(iQ,108)*BRAh(jQ,99)/Rmax,1d0)

!!---------  ee -> hA -> 2b 2taus

  Rmax=1d0
  IF(MAh(jQ).GT.Mhh(iQ))THEN
   DO K=1,NhA2b2tau
    IF(AINT(Massmax).EQ.hA2b2tau(K,1) .AND. AINT(Massmin).EQ.hA2b2tau(K,2))THEN
    Rmax=hA2b2tau(K,3)
    GOTO 1005
    ENDIF
   ENDDO
1005      PROB(27)=PROB(27)+DDIM(R*BRhh(iQ,108)*BRAh(jQ,175)/Rmax,1d0)
  ENDIF

  Rmax=1d0
  IF(Mhh(iQ).GT.MAh(jQ))THEN
   DO K=1,NhA2b2tau
    IF(AINT(Massmax).EQ.hA2b2tau(K,1) .AND. AINT(Massmin).EQ.hA2b2tau(K,2))THEN
     Rmax=hA2b2tau(K,3)
     GOTO 1006
    ENDIF
   ENDDO
1006      PROB(27)=PROB(27)+DDIM(R*BRhh(iQ,184)*BRAh(jQ,99)/Rmax,1d0)
  ENDIF

  Rmax=1d0
  IF(MAh(jQ).GT.Mhh(iQ))THEN
   DO K=1,NhA2tau2b
    IF(AINT(Massmax).EQ.hA2tau2b(K,1) .AND. AINT(Massmin).EQ.hA2tau2b(K,2))THEN
    Rmax=hA2tau2b(K,3)
    GOTO 1007
    ENDIF
   ENDDO
1007      PROB(27)=PROB(27)+DDIM(R*BRhh(iQ,184)*BRAh(jQ,99)/Rmax,1d0)
  ENDIF

  Rmax=1d0
  IF(Mhh(iQ).GT.MAh(jQ))THEN
   DO K=1,NhA2tau2b
    IF(AINT(Massmax).EQ.hA2tau2b(K,1) .AND. AINT(Massmin).EQ.hA2tau2b(K,2))THEN
    Rmax=hA2tau2b(K,3)
    GOTO 1008
    ENDIF
   ENDDO
1008      PROB(27)=PROB(27)+DDIM(R*BRhh(iQ,108)*BRAh(jQ,175)/Rmax,1d0)
  ENDIF

!!-----------  ee -> hA -> AAA -> 6bs
  IF(Mhh(iQ).GT.2d0*MAh(jQ))THEN
   Rmax=1d0
   DO K=1,NAAA6b
    IF(AINT(Mhh(iQ)).EQ.AAA6b(K,1) .AND. AINT(MAh(jQ)).EQ.AAA6b(K,2))THEN
    Rmax=AAA6b(K,3)
    GOTO 1009
    ENDIF
   ENDDO
1009      PROB(28)=PROB(28)+DDIM(R*BRHAA(iQ,jQ)*BRAh(jQ,175)**3/Rmax,1d0)

!!----------  ee -> hA -> AAA -> 6taus
   Rmax=1d0
   DO K=1,NAAA6tau
    IF(AINT(Mhh(iQ)).EQ.AAA6tau(K,1) .AND. AINT(MAh(jQ)).EQ.AAA6tau(K,2))THEN
    Rmax=AAA6tau(K,3)
    GOTO 1010
    ENDIF
   ENDDO
1010     PROB(29)=PROB(29)+DDIM(R*BRHAA(iQ,jQ)*BRAh(jQ,99)**3/Rmax,1d0)
  ENDIF

ENDIF !! MH+MA.LT. 209
ENDDO !! enddo over jQ = 2, 8 
ENDDO !! enddo over iQ = 1, 8 


!!*****************************************************************************
  DO iHp = 2, 8  !! loop over all 7 charged higgs

!!------------------------------------------------------------------------
!!--- top -> H+ b, H+ -> c s (from CDF, 0907.1269, and D0, 0908.1811) ----
!!------------------------------------------------------------------------    
!BRFu(3,iHp+94) = BR(t-> H+ b)
!BRHpm(iHp,111) = BR(H+-> c s)
    
  brtopcs  = BRFu(3,iHp+94)*BRHpm(iHp,111)        !!      BR_tHb(iHp)*BR_Hcs(iHp) 
  PROB(30) = 0d0
      
  DO I=1,9
   IF((minf(i).LE.MHpm(iHp)).AND.(MHpm(iHp).LE.msup(i))) THEN
   brtoplim = binf(i)+((MHpm(iHp)-minf(i))*(bsup(i)-binf(i)))/(msup(i)-minf(i))
   PROB(30) = PROB(30)+DDIM(brtopcs/brtoplim,1d0)
   ENDIF
  ENDDO

!   Write(17,*) 'BR_Hcs(iHp)=',BR_Hcs(iHp), 'BRHpm(iHp,111)=',BRHpm(iHp,111)
!!------------------------------------------------------------
!!--- top -> H+ b; H+ -> tau nu_tau (from D0, 0908.1811) -----
!!------------------------------------------------------------
!BRHpm(iHp,69)  = BR(H+-> tau nu_tau)

  brtoptau = BRFu(3,iHp+94)*BRHpm(iHp,69)  ! BR_tHb(iHp)*BR_Htaunu(iHp),
  PROB(31) = 0d0
  
  DO I=1,7
   IF((minf(i+2).LE.MHpm(iHp)).AND.(MHpm(iHp).LE.msup(i+2))) THEN
   brtoplim=brtau(i)+((MHpm(iHp)-minf(i+2))*(brtau(i+1)-brtau(i)))/(msup(i+2)-minf(i+2))
   PROB(31)=PROB(31)+DDIM(brtoptau/brtoplim,1d0)
   ENDIF
  ENDDO

!!---------------------------------------------------------------------
!!--- top -> H+ b; H+ -> W+ A; A -> 2taus (from CDF Note 10104) ---
!!---------------------------------------------------------------------
! BRHpm(iHp,jQ+48) = BR(H+-> W+ A)
! BRAh(jQ,99) = BR(A -> 2 taus)
! MA < 2Mb

  DO jQ = 2, 8 !! over A
   brtopa1 = BRFu(3,iHp+94)*BRHpm(iHp,jQ+48)*BRAh(jQ,99)
   PROB(32)=0d0

   IF(MAh(jQ).LE.4d0) then     
    DO I=1,4
     IF((mh4m(i).LE.MHpm(iHp)).AND.(MHpm(iHp).LE.mh4p(i))) THEN
     brtoplim = br4m(i)+(MHpm(iHp)-mh4m(i))*(br4p(i)-br4m(i))/(mh4p(i)-mh4m(i))
     PROB(32) = PROB(32)+DDIM(brtopa1/brtoplim,1d0)
     ENDIF
    ENDDO
   ENDIF

   IF((MAh(jQ).gt.4d0).and.(MAh(jQ).LE.7d0)) then     
    DO I=1,3
     IF((mh7(i).LE.MHpm(iHp)).AND.(MHpm(iHp).LE.mh7(i+1))) THEN
     brtoplim=br4p(i)+(MAh(jQ)-4d0)*(br7(i)-br4p(i))/3d0 +(br4p(i+1)-br4p(i) &
&    +(br7(i+1)-br7(i)-br4p(i+1)+br4p(i))*(MAh(jQ)-4d0)/3d0) &
&    *(MHpm(iHp)-mh7(i))/(mh7(i+1)-mh7(i))
     PROB(32)= PROB(32)+DDIM(brtopa1/brtoplim,1d0)
     ENDIF
    ENDDO
   ENDIF

   IF((MAh(jQ).gt.7d0).and.(MAh(jQ).LE.9d0)) then     
    DO I=1,4
     IF((mh7(i).LE.MHpm(iHp)).AND.(MHpm(iHp).LE.mh7(i+1))) THEN
     brtoplim=br7(i)+(MAh(jQ)-7d0)*(br9(i)-br7(i))/2d0 +(br7(i+1)-br7(i)  &
&    +(br9(i+1)-br9(i)-br7(i+1)+br7(i))*(MAh(jQ)-7d0)/2d0) &
&    *(MHpm(iHp)-mh7(i))/(mh7(i+1)-mh7(i))
     PROB(32) = PROB(32)+DDIM(brtopa1/brtoplim,1d0)
     ENDIF
    ENDDO
   ENDIF
 
   ENDDO
 ENDDO


!!******************************************************************
!!** Constraints from ggF/bb->H/A->tautau from CMS-PAS-HIG-13-021 **
!!------------------------------------------------------------------
!! Loop over all 15 Higgs (CP odd + even)

DO I=1,15

!Write(23,*) ' I =' , I

  IF(I.LE.8) MH(I)=Mhh(I)
  IF(I.GE.9) MH(I)=MAh(I-7)
  J=1
  DO WHILE(HMAS(J).LE.MH(I) .AND. J.LT.NX)
  J=J+1
  ENDDO
  IF(J.GE.2 .AND. MH(I).LT.HMAS(NX)) THEN
   XSMH(I)=0d0
   LCMSH(I)=1d10
   XSMHbb(I)=0d0
   LCMSHbb(I)=1d10

! SM Higgs ggF prod. cross sect.:
   XSMH(I)=XSM(J-1)+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(XSM(J)-XSM(J-1))

! ggF Signal cross section*BR:
   IF(I.LE.8) SIG(I)= Real(ratioGG(I),dp)**2*BRhh(I,108)*XSMH(I)
   IF(I.GE.9) SIG(I)= Real(ratioPGG(I-7),dp)**2*BRAh(I-7,99)*XSMH(I)

! CMS ggF limit:
   LCMSH(I)= LCMS(J-1)+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(LCMS(J)-LCMS(J-1))

   PROB(51)= PROB(51)+DDIM(1d0,LCMSH(I)/SIG(I))

!Write(23,*) 'CMS ggF limit,   ggFSM*BR'
!Write(23,*) 'LCMSH(I), SIG(I)',LCMSH(I),SIG(I)
!Write(23,*) '  ' 

! ATLAS ggF limit:
   LATLASH(I)= LATLASgg(J-1)+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(LATLASgg(J)-LATLASgg(J-1))

! Correct for jump in Fig.7 at MA=200 GeV: J=8, 
! modif. LATLASgg(J-1)=LATLASgg(7)=.96D0 and not .794d0:
   IF(J.EQ.8) THEN
   LATLASH(I)=.96D0+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(LATLASgg(J)-.96D0)
   ENDIF
   PROB(51)= PROB(51)+DDIM(1D0,LATLASH(I)/SIG(I))
 
!Write(23,*) 'ATLAS ggF limit,    ggFSM*BR' 
!Write(23,*) ' LATLASH(I) , SIG(I)', LATLASH(I) , SIG(I)  !!DK
!Write(23,*) '  ' 

!!***********************************************
!!------- SM Higgs bbH prod. cross sect.:
  XSMHbb(I)= XSMbb(J-1)+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(XSMbb(J)-XSMbb(J-1))

!!------- bbH Signal cross section*BR:
  IF(I.LE.8) SIGbb(I) = rHB_S_S_Fd(I,3)**2*XSMHbb(I)*BRhh(I,108)
  IF(I.GE.9) SIGbb(I) = rHB_P_P_Fd(I-7,3)**2*XSMHbb(I)*BRAh(I-7,99)

!!------- CMS Hbb limit:
   LCMSHbb(I)= LCMSbb(J-1)+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(LCMSbb(J)-LCMSbb(J-1))
   PROB(51)= PROB(51)+DDIM(1d0,LCMSHbb(I)/SIGbb(I))

!Write(23,*) 'CMS bbH limit,    bbHSM*BR' 
!Write(23,*) ' LCMSHbb(I), SIGbb(I)', LCMSHbb(I), SIGbb(I)  !!DK
!Write(23,*) '  ' 

!!------- ATLAS Hbb limit:
   LATLASHbb(I) = LATLASbb(J-1)+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(LATLASbb(J)-LATLASbb(J-1))

!!------- Correct for jump in Fig.7 at MA=200 GeV: J=8, 
!!------- modif. LATLASbb(J-1)=LATLASbb(7)=.858D0 and not .393d0:
 IF(J.EQ.8) THEN
 LATLASHbb(I) = .858d0+((MH(I)-HMAS(J-1))/(HMAS(J)-HMAS(J-1)))*(LATLASbb(J)-.858d0)
 ENDIF
 PROB(51) = PROB(51)+DDIM(1d0,LATLASHbb(I)/SIGbb(I))

!Write(23,*) 'ATLAS bbH limit,   bbHSM*BR'
!Write(23,*) ' LATLASHbb(I), SIGbb(I)', LATLASHbb(I), SIGbb(I)  !!DK
!Write(23,*) '  ' 

!!***********************************************
!---- Combine signal rates of 2 Higgses ---------
!!-----------------------------------------------
IF(I .GT. 1) THEN  !! BY DK on 25.02.16

 DO I1=1,I-1

!  Write(23,*) ' I1 = ', I1 

   J1=1
   DO WHILE(HMAS(J1).LE.MH(I1) .AND. J1.LT.NX)
    J1 = J1+1
   ENDDO
 IF(MH(I1).GE.HMAS(1) .AND. MH(I1).LT.HMAS(NX)) THEN
! Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
   MBAR = (SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
   JBAR = 1
   DO WHILE(HMAS(JBAR).LE.MBAR.AND.JBAR.LT.NX)
   JBAR = JBAR+1
   ENDDO
   MBARbb = (SIGbb(I)*MH(I)+SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
   JBARbb = 1
   DO WHILE(HMAS(JBARbb).LE.MBARbb.AND.JBARbb.LT.NX)
   JBARbb = JBARbb+1
   ENDDO
! DEL=mass difference divided by a (small) resolution squared:
! [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0; otherwise the combined signal rate is small]
   DEL = ((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2
! Estimate of the combined ggF signal rates:
   SIGTOT = SIG(I)+SIG(I1)-SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
! Continue only if SIGTOT > 0 and 90<MBAR<1000 and |MH(I)-MH(I1)|/MBAR<0.20:
   IF(SIGTOT.GT.0D0.AND.MBAR.GE.HMAS(1).AND.MBAR.LT.HMAS(NX).AND.dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN
! CMS ggF limit at MBAR:
     LCMSMB = LCMS(JBAR-1)+((MBAR-HMAS(JBAR-1))/(HMAS(JBAR)-HMAS(JBAR-1)))*(LCMS(JBAR)-LCMS(JBAR-1))
     PROB(51) = PROB(51)+DDIM(1d0,LCMSMB/SIGTOT)
!Write(23,*) 'CMS ggF limit,    bbHSM*BR,  combined rate'
!Write(23,*) 'LCMSMB, SIGTOT', LCMSMB, SIGTOT  !!DK
!Write(23,*) '  ' 

! ATLAS ggF limit at MBAR:
       LATLASMB = LATLASgg(JBAR-1)+((MBAR-HMAS(JBAR-1))/(HMAS(JBAR)-HMAS(JBAR-1)))*(LATLASgg(JBAR)-LATLASgg(JBAR-1))
! Correct for jump in Fig.7 at MA=200 GeV: JBAR=8, modif. LATLASgg(JBAR-1)=LATLASgg(7)=.96D0 and not .794d0:
   IF(J.EQ.8) THEN
     LATLASMB = .96D0+((MBAR-HMAS(JBAR-1))/(HMAS(JBAR)-HMAS(JBAR-1)))*(LATLASgg(JBAR)-.96D0)
   ENDIF
    PROB(51)= PROB(51)+DDIM(1D0,LATLASMB/SIGTOT)

   ENDIF

!Write(23,*) 'ATLAS ggH limit,    bbHSM*BR'
!Write(23,*) 'LATLASMB, SIGTOT', LATLASMB, SIGTOT  !!DK
!Write(23,*) '  ' 

! Estimate of the combined bbH signal rates:
   SIGTOTbb = SIGbb(I)+SIGbb(I1)-SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
! Continue only if SIGTOTbb > 0 and 90<MBARbb<1000  and |MH(I)-MH(I1)|/MBARbb<0.20: 
   IF(SIGTOTbb.GT.0D0.AND.MBARbb.GE.HMAS(1).AND. MBARbb.LT.HMAS(NX).AND.dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN
! CMS bbH limit at MBARbb:
    LCMSMBbb = LCMSbb(JBARbb-1)+((MBARbb-HMAS(JBARbb-1))/(HMAS(JBARbb)-HMAS(JBARbb-1)))*(LCMSbb(JBARbb)-LCMSbb(JBARbb-1))
    PROB(51) = PROB(51)+DDIM(1d0,LCMSMBbb/SIGTOTbb)

!Write(23,*) 'LCMSMBbb, SIGTOTbb', LCMSMBbb, SIGTOTbb  !!DK
!Write(23,*) '  ' 

! ATLAS bbH limit at MBARbb:
    LATLASMBbb = LATLASbb(JBARbb-1)+((MBARbb-HMAS(JBARbb-1))/(HMAS(JBARbb)-HMAS(JBARbb-1)))*(LATLASbb(JBARbb)-LATLASbb(JBARbb-1))
! Correct for jump in Fig.7 at MA=200 GeV: JBARbb=8, modif. LATLASbb(JBARbb-1)=LATLASbb(7)=.858D0 and not .393d0:
    IF(J.EQ.8) THEN
     LATLASMBbb = .858D0+((MBARbb-HMAS(JBARbb-1))/(HMAS(JBARbb)-HMAS(JBARbb-1)))*(LATLASbb(JBARbb)-.858D0)
    ENDIF
     PROB(51) = PROB(51)+DDIM(1D0,LATLASMBbb/SIGTOTbb)
    ENDIF
!Write(23,*) 'LATLASMBbb, SIGTOTbb', LATLASMBbb, SIGTOTbb  !!DK
!Write(23,*) '  ' 
!Write(23,*) '  ' 
!Write(23,*) '  ' 
   ENDIF
  ENDDO   

ENDIF !! BY DK on 25.02.16  
    
 ENDIF
ENDDO

!!--------------------------------------------------------------
!!--- CMS_AA_4MU,  Sparticle searches at LHC (NMSPEC) ----------
!!--------------------------------------------------------------
!!  iH correspond to the position of the correct 125 GeV higgs
 IF(iH.GE.1) THEN
DO jQ = 2, 8
 J=1
 DO WHILE(MAh(jQ) .GE. XX1(1) .AND. XX1(J) .LE. MAh(jQ) .AND. J.LT.NY)
  J=J+1
 ENDDO
 IF(J.GT.1 .AND. MAh(jQ).LT.XX1(NY))THEN
  LIMIT = LIMIT + LSIGBR(J-1)+ ((MAh(jQ)-XX1(J-1))/(XX1(J)-XX1(J-1)))*(LSIGBR(J)-LSIGBR(J-1))
  PROB(54) = DDIM(Real(ratioGG(iH),dp)**2*BRHAA(iH,jQ)*BRAh(jQ,93)**2/LIMIT,1d0)
 ENDIF
ENDDO
 ENDIF

!!*************************************************************
!! Constraints from ggF->H/A->gamgam from ATLAS-CONF-2014-031, M_H/A < 122
!!*************************************************************
! Read ATLAS upper limit
! ggHgg(IQ,1): Higgs mass
! ggHgg(IQ,4): upper limit in fb

! rruiz put above
!NggHgg = 0 
!OPEN(unit = 274, file = trim(EXPCON_PATH)//'EXPCON/ggHgg.dat', status='old') !
!DO iQ=1, 100
!  READ(274,*,IOSTAT=ios) ggHgg(iQ,1), ggHgg(iQ,2), ggHgg(iQ,3), ggHgg(iQ,4)
!  IF (ios /= 0) EXIT
!  NggHgg = NggHgg +1
!ENDDO

! Read SM Higgs ggF production cross section (60-122 GeV)
! SMXS(I,1): Higgs mass
! SMXS.65-122.dat:  ggF production cross section in pb
! SMXS(I,2): ggF production cross section in fb
!NSMXS = 0 
!OPEN(unit = 275, file = trim(EXPCON_PATH)//'EXPCON/SMXS.65-122.dat', status='old') !
!DO iQ=1, 100
!  READ(275,*,IOSTAT=ios) SMXS(iQ,1), SMXS(iQ,2), SMXS(iQ,3), SMXS(iQ,4), SMXS(i!Q,5), SMXS(iQ,6)
!  SMXS(iQ,2)=1.d3*SMXS(iQ,2)
!  IF (ios /= 0) EXIT
!  NSMXS = NSMXS +1
!ENDDO

!! Loop over 15 Higgses
   DO I=1,15
     IF(I.LE.8) MH(I)=Mhh(I)
     IF(I.GE.9) MH(I)=MAh(I-7)
      J=1
      DO WHILE(SMXS(J,1).LE.MH(I) .AND. J.LT.70)
        J=J+1
      ENDDO

     IF(J.GE.2 .AND. MH(I).LT.122d0) THEN
      XSMH(I)=0d0
      LATLASH(I)=1d10
! SM Higgs ggF prod. cross sect.:
      XSMH(I)=SMXS(J-1,2)+((MH(I)-SMXS(J-1,1))/(SMXS(J,1)-SMXS(J-1,1)))*(SMXS(J,2)-SMXS(J-1,2))

! ggF Signal cross section*BR(H->gamgam):
   IF(I.LE.8) SIG(I)= Real(ratioGG(I),dp)**2*BRhh(I,2)*XSMH(I)
   IF(I.GE.9) SIG(I)= Real(ratioPGG(I-7),dp)**2*BRAh(I-7,2)*XSMH(I)

! ATLAS limit:
      LATLASH(I)=ggHgg(J-1,4)+((MH(I)-SMXS(J-1,1))/(SMXS(J,1)-SMXS(J-1,1)))*(ggHgg(J,4)-ggHgg(J-1,4))
      PROB(64)=PROB(64)+DDIM(1d0,LATLASH(I)/SIG(I))

     ENDIF
    ENDDO

!!*************************************************************
!!*************************************************************
!!-----------------Light A Physics-----------------------------
alphaSamZ = AlphaSDR(9115.35d-2,MSd,MSu,MGlu,MFd,MFu)
alphaSamY = AlphaSDR(946.d-2,MSd,MSu,MGlu,MFd,MFu)
PI_Pj=4d0*DATAN(1d0) !! PI
ALEMMVZ = alpha_MSbar(mZ,mW)
gPj=(Real(g1,dp)+Real(g2,dp))/2d0
N_lightPj = 0
JPj = 0

!!*******************************************************
 DO jQ = 2, 8

!------ Test on Upsilon(1S) -> A gamma (from CLEO)
!! Available only for MA < MY

      MY=9.46d0 ! Upsilon(1S) mass
      MBQM=4.9d0 ! b quark mass in quark models
      ALSMY=alphaSamY ! alpha_s at MY, 2 loop calculation
      BRYMUMU=2.48d-2 ! BR(Upsilon(1S) -> mu mu)
      IF(MAh(jQ).LT.MY)THEN

      ZZ_Pj=1d0-MAh(jQ)**2/MY**2 ! energy fraction of the photon
      AP=6d0*ZZ_Pj**.2d0 ! Nason function for QCD corrections
      RC=1d0+4d0*ALSMY/(3d0*PI_Pj)*(4d0-AP) ! QCD corrections
      DELTA=1.2d0**2/MBQM**2 ! function for rel. corrections
      RC=RC* & ! relativistic corrections (for MA<~8.8 GeV)
     &  (MY**2-MAh(jQ)**2)**2/(4d0*MBQM**2-MAh(jQ)**2)**2*(1d0- &
     &  DELTA/3d0*(36d0*MBQM**2+MAh(jQ)**2)/(4d0*MBQM**2-MAh(jQ)**2))

      RC=MAX(RC,1d-6)

      RMAXPj=0d0
      RMAXPj=DSQRT(2d0)*PI_Pj*ALEMMVZ*CLEOTAU(MAh(jQ))/(G_F*MBQM**2*ZZ_Pj*RC*BRYMUMU)
      IF(RMAX.NE.0d0) PROB(65)=DDIM(rHB_P_P_Fd(jQ,3)**2*BRAh(jQ,99)/RMAX,1d0)

      RMAXPj=0d0
      RMAXPj=DSQRT(2d0)*PI_Pj*ALEMMVZ*CLEOMU(MAh(jQ))/(G_F*MBQM**2*ZZ_Pj*RC*BRYMUMU)
      IF(RMAX.NE.0d0) PROB(65)=DDIM(rHB_P_P_Fd(jQ,3)**2*BRAh(jQ,93)/RMAX,1d0)

      IF(gTAh(jQ).GT.1d-2) PROB(65) = -PROB(65)

!!   Write(17,*) 'PROB(80)= ',PROB(80)
     ENDIF
     
      IF( (9.4d0 .le. MAh(jQ)).and. (MAh(jQ).le. 10.5d0)) THEN
       N_lightPj = N_lightPj + 1
       JPj = jQ
      ENDIF

 ENDDO

!!*************************************************************
!!----- Test on etab(1S) mass difference (BABAR - theory) O
   
     IF(N_lightPj .eq. 1) THEN  !! constraint avaible only one A satisfies 9.4<mA<10.5 GeV

      M0Pj = 9.389d0 ! etab(1S) mass
      GYEE = 1.34d-6 ! Gamma(Upsilon(1S) -> e+ e-)
      RETA = GYEE*9d0*MY**2/(4d0*ALEMMVZ**2)*(1d0+16d0*ALSMY/(3d0*PI)) ! radial wave fun. at the origin
       ! Resolution of the 3rd degree eq. for the limit on Xd
      GAM2=(M0Pj*2d-2)**2
      XXPj = MAh(JPj)**2-M0Pj**2
      IF(MAh(JPj).LT.M0Pj)THEN
       MEMAX = M0Pj-3.d-2
      ELSE
       MEMAX = M0Pj+4.d-2
      ENDIF
      YMAX = MEMAX**2-M0Pj**2
      FMAX = XXPj*YMAX*(1d0+GAM2/(XXPj+YMAX)**2)

      DPj = XXPj**2-GAM2/27d0
      IF(DPj.LT.0d0)THEN
       UUPj = 2d0*DSQRT(GAM2)/DSQRT(3d0)
       VVPj = -3d0*DSQRT(3d0)*XXPj/DSQRT(GAM2)
       YYPj = UUPj*DCOS(DACOS(VVPj)/3d0 + 4d0*PI/3d0)
       FPj  = XXPj*YYPj*(1d0+GAM2/(XXPj+YYPj)**2)
       FMAX = MAX(FMAX,FPj)
      ENDIF

      RMAX=8d0*PI*MZ**2/(3d0*gPj*RETA*M0Pj**3)*FMAX
      PROB(66)=DDIM(rHB_P_P_Fd(jQ,3)**2/RMAX,1d0)
 
      ENDIF

!!*************************************************************
!!*************************************************************

!!*************************************************************

!if(debug) then

If(PROB(53).ne.0.d0) then
   Write(17,*) '# Warning: no Higgs in the 122-128 GeV range'
Else
  If(PROB(55).ne.0.d0) then
   Write(17,*) '# Warning: 80% of the 125 GeV Higgs composition is not a mixture of Hu & Hd'
  else
     If(PROB(52).ne.0.d0) then
     Write(17,*) '# Warning: 80% of the 125 GeV Higgs composition is not from Hu'
     Endif
  Endif
Endif

!IF(PROB(1).NE.0d0) Then
!  WRITE(17,*) '# Warning: selectron or smuon too light '
!Endif
!IF(PROB(2).NE.0d0) then
! WRITE(17,*) '# Warning: stau too light '
!Endif
IF(PROB(3).NE.0d0) then
 WRITE(17,*) '# Warning: charged Higgs mass excluded by LEP '
Endif
IF(PROB(4).NE.0.d0) THEN
  WRITE(17,*) '# Warning: chargino mass excluded by LEP'     
ENDIF
!IF(PROB(5).NE.0d0) THEN
!    WRITE(17,*) '# Warning: squark or gluino too light '
!ENDIF 
!IF(PROB(6).NE.0d0) THEN
!    WRITE(17,*) '# Warning: scalar bottom is too light '
!ENDIF 
!IF(PROB(7).NE.0d0) THEN
!    WRITE(17,*) '# Warning: scalar top is too light '
!ENDIF 
If (PROB(8).NE.0.d0) Then
    Write(17,*) '# Warning: 80% of the neutral or Charged Goldstone is not from Hu or Hd'
Endif
If(PROB(9).NE.0.d0) Then
 Write(17,*) '# Warning: Z decay width too large, excluded'
Endif
If(PROB(10).ne.0.d0) then
  Write(17,*) '# Warning: Total SM-like Higgs total decay width is too large, > 22MeV(CMS)'
Endif
If(PROB(11).ne.0.d0) then
   Write(17,*) '# Warning: excluded by e e -> Z h '
Endif
If(PROB(12).ne.0.d0) then
   Write(17,*) '# Warning: excluded by e e -> Z h, h -> b b'
Endif
If(PROB(13).ne.0.d0) then
   Write(17,*) '# Warning: excluded by e e -> Z h, h -> j j'
Endif
If(PROB(14).ne.0.d0) then
   Write(17,*) '# Warning: excluded by e e -> Z h, h -> tau tau '
Endif
If(PROB(15).ne.0.d0) then
   Write(17,*) '# Warning: excluded by e e -> Z h, h -> 2 photons '
Endif
If(PROB(21).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> Z h, h-> hh/ AA -> 4bs '
Endif
If(PROB(22).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> hZ, h-> hh/AA -> 4taus'
Endif
If(PROB(41).ne.0d0) then
   Write(17,*) '#  excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)'
Endif
If(PROB(23).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> Z h, h-> hh/ AA -> 2bs 2taus '
Endif
If(PROB(24).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> Z h -> ZAA -> Z + light pairs '
Endif
If(PROB(25).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> hA -> 4bs'
Endif
If(PROB(26).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> hA -> 4taus'
Endif
If(PROB(27).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> hA -> 2bs 2taus'
Endif
If(PROB(28).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> hA -> AAA -> 6bs'
Endif
If(PROB(29).ne.0.d0) then
   Write(17,*) '# Warning: excluded by ee -> hA -> AAA -> 6taus'
Endif
If(PROB(30).ne.0.d0) then
   Write(17,*) '# Warning: excluded by top -> b H+, H+ -> c s (CDF, D0)'
Endif
If(PROB(31).ne.0.d0) then
   Write(17,*) '# Warning: excluded by top -> b H+, H+ -> tau nu_tau (D0)'
Endif
If(PROB(32).ne.0.d0) then
   Write(17,*) '# Warning: excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)'
Endif
If(PROB(51).ne.0.d0) then
   Write(17,*) '# Warning: excluded by H/A->tautau'
Endif
If(PROB(54).ne.0.d0) then
   Write(17,*) '# Warning: excluded H_125->AA->4mu (CMS)'
Endif

!IF(PROB(60).NE.0d0) THEN
!   WRITE(17,*) '# Warning: excluded by stop -> b l sneutrino'
!ENDIF
!IF(PROB(61).NE.0d0) THEN
!  WRITE(17,*) '# Warning: excluded by stop -> neutralino c'
!ENDIF
!IF(PROB(62).NE.0d0) THEN
!  WRITE(17,*) '# Warning: excluded by sbottom -> neutralino b'
!ENDIF

If(PROB(64).ne.0.d0) then
  Write(17,*) '# Warning: excluded by ggF->H/A->gamgam (65GeV < M < 122GeV, ATLAS)'
Endif

If(PROB(65).ne.0.d0) then
  Write(17,*) '# Warning: excluded by Upsilon(1S) -> A gamma'
Endif

If(PROB(66).ne.0.d0) then
  Write(17,*) '# Warning: excluded by eta_b(1S) mass measurement'
Endif



!!********************************************************


!999 continue

!write(20,*) ' '
!write(21,*) ' '
!write(22,*) ' '
write(17,*) ' '

!endif

! ENDDO !! end for number of scan points


  call date_and_time(DATE=date,ZONE=zone)
  call date_and_time(TIME=time)
  Write(17,*) "Finished: " , date(7:8),".",date(5:6),".",date(1:4)," , ",time(1:2),":",time(3:4),":",time(5:6)
write(17,*) ' '
write(17,*) ' '

ENDif!!(debug)

End subroutine constraints

!!**********************************************************
DOUBLE PRECISION FUNCTION CLEOTAU(MX)

!  CLEO constraints on BR(Y -> A gamma)*BR(A -> tau tau)

 IMPLICIT NONE

  INTEGER I,N
   PARAMETER (N=17)
   DOUBLE PRECISION MX,X(N),M(N)

   DATA M/3.75d0,4.25d0,4.75d0,5.1d0,5.8d0,6.15d0,6.6d0,7d0, &
  &      7.4d0,7.6d0,8d0,8.25d0,8.6d0,9d0,9.25d0,9.35d0,9.41d0/
   DATA X/2.9d-5,2.5d-5,2d-5,2.3d-5,5.1d-5,2.5d-5,2.5d-5,2.7d-5,  &
  &     4.5d-5,3.7d-5,2.7d-5,7.2d-5,6.8d-5,8.6d-5,2.1d-4,2.85d-4, &
  &     4.75d-4/

  CLEOTAU = 0d0

  IF(MX.LT.M(1).OR.MX.GT.M(N))RETURN

  DO I=2,N
    IF(MX.LT.M(I))THEN
      CLEOTAU=X(I-1)+(MX-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
    RETURN
    ENDIF
  ENDDO

END function  CLEOTAU
!!*************************************************************

      DOUBLE PRECISION FUNCTION CLEOMU(MX)

!  CLEO constraints on BR(Y -> A gamma)*BR(A -> mu mu)

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION MX,X(N),M(N)

      DATA M/.25d0,3.75d0/
      DATA X/9.d-6,9.d-6/

      CLEOMU=0d0

      IF(MX.LT.M(1).OR.MX.GT.M(N))RETURN

      DO I=2,N
       IF(MX.LT.M(I))THEN
        CLEOMU=X(I-1)+(MX-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
        RETURN
       ENDIF
      ENDDO

      END function cleomu
!!*************************************************************











