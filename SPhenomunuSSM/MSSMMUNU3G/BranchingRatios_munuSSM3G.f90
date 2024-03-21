! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:25 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module BranchingRatios_munuSSM3G 
 
Use Control 
Use Couplings_munuSSM3G 
Use Model_Data_munuSSM3G 
Use LoopCouplings_munuSSM3G 
Use Glu3Decays_munuSSM3G 
Use Sd3Decays_munuSSM3G 
Use Su3Decays_munuSSM3G 
Use Cha3Decays_munuSSM3G 
Use Chi3Decays_munuSSM3G 
Use SUSYDecays_munuSSM3G 
 
Contains 
 
Subroutine CalculateBR(CTBD,fac3,epsI,deltaM,kont,MAh,MAh2,MCha,MCha2,MChi,           & 
& MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,              & 
& MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,             & 
& vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,             & 
& mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,gPSd,gTSd,BRSd,gPSu,gTSu,BRSu,gPhh,gThh,           & 
& BRhh,gPAh,gTAh,BRAh,gPHpm,gTHpm,BRHpm,gPGlu,gTGlu,BRGlu,gPFu,gTFu,BRFu,gPCha,          & 
& gTCha,BRCha,gPChi,gTChi,BRChi)

Real(dp), Intent(in) :: epsI, deltaM, fac3 
Integer, Intent(inout) :: kont 
Logical, Intent(in) :: CTBD 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Real(dp),Intent(inout) :: gPSd(6,2880),gTSd(6),BRSd(6,2880),gPSu(6,2880),gTSu(6),BRSu(6,2880),gPhh(8,359),      & 
& gThh(8),BRhh(8,359),gPAh(8,356),gTAh(8),BRAh(8,356),gPHpm(8,223),gTHpm(8),             & 
& BRHpm(8,223),gPGlu(1,261),gTGlu,BRGlu(1,261),gPFu(3,168),gTFu(3),BRFu(3,168),          & 
& gPCha(5,1010),gTCha(5),BRCha(5,1010),gPChi(10,1729),gTChi(10),BRChi(10,1729)

Complex(dp) :: cplHiggsPP(8),cplHiggsGG(8),cplPseudoHiggsPP(8),cplPseudoHiggsGG(8),cplHiggsZZvirt(8),& 
& cplHiggsWWvirt(8)

Real(dp) :: gGluFdcFdChi(1,3,3,10),gGluFdcFucCha(1,3,3,5),gGluFucFuChi(1,3,3,10),gSdAhChaFu(6,8,5,3),& 
& gSdAhChiFd(6,8,10,3),gSdAhFdGlu(6,8,3,1),gSdSuChaChi(6,6,5,10),gSdChaFdcHpm(6,5,3,8),  & 
& gSdhhChaFu(6,8,5,3),gSdChaGluSu(6,5,1,6),gSdSdChacCha(6,6,5,5),gSdSdChiChi(6,6,10,10), & 
& gSdhhChiFd(6,8,10,3),gSdHpmChiFu(6,8,10,3),gSdChiGluSd(6,10,1,6),gSdFdFdcSd(6,3,3,6),  & 
& gSdFuFdcSu(6,3,3,6),gSdHpmFdcCha(6,8,3,5),gSdSdFdcFd(6,6,3,3),gSdSuFdcFu(6,6,3,3),     & 
& gSdSdFucFu(6,6,3,3),gSdhhFdGlu(6,8,3,1),gSdHpmFuGlu(6,8,3,1),gSdGluGluSd(6,1,1,6),     & 
& gSuAhChiFu(6,8,10,3),gSuAhFdcCha(6,8,3,5),gSuAhFuGlu(6,8,3,1),gSuSuChiChi(6,6,10,10),  & 
& gSucHpmChiFd(6,8,10,3),gSuhhChiFu(6,8,10,3),gSuChiGluSu(6,10,1,6),gSuSdChicCha(6,6,10,5),& 
& gSuFdFucSd(6,3,3,6),gSuhhFdcCha(6,8,3,5),gSuSuFdcFd(6,6,3,3),gSucHpmChaFu(6,8,5,3),    & 
& gSuFuFucSu(6,3,3,6),gSucChaFuHpm(6,5,3,8),gSuSdFucFd(6,6,3,3),gSuSuFucFu(6,6,3,3),     & 
& gSucHpmFdGlu(6,8,3,1),gSuhhFuGlu(6,8,3,1),gSuGluGluSu(6,1,1,6),gSuGluSdcCha(6,1,6,5),  & 
& gSuSuChacCha(6,6,5,5),gChaChacChaCha(5,5,5,5),gChaChaChiChi(5,5,10,10),gChaChacFdFd(5,5,3,3),& 
& gChaChacFuFu(5,5,3,3),gChaChicFuFd(5,10,3,3),gChaFdcFuGlu(5,3,3,1),gChiChicChaCha(10,10,5,5),& 
& gChiChiChiChi(10,10,10,10),gChiChicFdFd(10,10,3,3),gChiChicFuFu(10,10,3,3),            & 
& gChiChacFdFu(10,5,3,3),gChiFdcFdGlu(10,3,3,1),gChiFucFuGlu(10,3,3,1)

Complex(dp) :: coup 
Real(dp) :: vev 
Real(dp) :: gTVZ,gTVWm

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateBR'
 
!Write(*,*) "Calculating branching ratios and decay widths" 
gTVWm = gamW 
gTVZ = gamZ 
gPSd = 0._dp 
gTSd = 0._dp 
BRSd = 0._dp 
If (.Not.CTBD) Then 
Call SdTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,            & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,             & 
& pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,              & 
& Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,              & 
& M2,M3,vd,vu,vL,vR,gPSd(:,1:192),gTSd,BRSd(:,1:192))

If (Enable3BDecaysS) Then 
If (MaxVal(gTSd).Lt.MaxVal(fac3*Abs(MSd))) Then 
Call SdThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gSdAhChaFu,gSdAhChiFd,gSdAhFdGlu,gSdSuChaChi,gSdChaFdcHpm,              & 
& gSdhhChaFu,gSdChaGluSu,gSdSdChacCha,gSdSdChiChi,gSdhhChiFd,gSdHpmChiFu,gSdChiGluSd,    & 
& gSdFdFdcSd,gSdFuFdcSu,gSdHpmFdcCha,gSdSdFdcFd,gSdSuFdcFu,gSdSdFucFu,gSdhhFdGlu,        & 
& gSdHpmFuGlu,gSdGluGluSd,epsI,deltaM,.False.,gTSd,gPSd(:,193:2880),BRSd(:,193:2880))

Else 
Call SdThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gSdAhChaFu,gSdAhChiFd,gSdAhFdGlu,gSdSuChaChi,gSdChaFdcHpm,              & 
& gSdhhChaFu,gSdChaGluSu,gSdSdChacCha,gSdSdChiChi,gSdhhChiFd,gSdHpmChiFu,gSdChiGluSd,    & 
& gSdFdFdcSd,gSdFuFdcSu,gSdHpmFdcCha,gSdSdFdcFd,gSdSuFdcFu,gSdSdFucFu,gSdhhFdGlu,        & 
& gSdHpmFuGlu,gSdGluGluSd,epsI,deltaM,.True.,gTSd,gPSd(:,193:2880),BRSd(:,193:2880))

End If 

End If 
Else 
Call SdThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gSdAhChaFu,gSdAhChiFd,gSdAhFdGlu,gSdSuChaChi,gSdChaFdcHpm,              & 
& gSdhhChaFu,gSdChaGluSu,gSdSdChacCha,gSdSdChiChi,gSdhhChiFd,gSdHpmChiFu,gSdChiGluSd,    & 
& gSdFdFdcSd,gSdFuFdcSu,gSdHpmFdcCha,gSdSdFdcFd,gSdSuFdcFu,gSdSdFucFu,gSdhhFdGlu,        & 
& gSdHpmFuGlu,gSdGluGluSd,epsI,deltaM,.False.,gTSd,gPSd(:,193:2880),BRSd(:,193:2880))

End If 
Do i1=1,6
gTSd(i1) =Sum(gPSd(i1,:)) 
If (gTSd(i1).Gt.0._dp) BRSd(i1,: ) =gPSd(i1,:)/gTSd(i1) 
End Do 
 

gPSu = 0._dp 
gTSu = 0._dp 
BRSu = 0._dp 
If (.Not.CTBD) Then 
Call SuTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,            & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,             & 
& pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,              & 
& Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,              & 
& M2,M3,vd,vu,vL,vR,gPSu(:,1:192),gTSu,BRSu(:,1:192))

If (Enable3BDecaysS) Then 
If (MaxVal(gTSu).Lt.MaxVal(fac3*Abs(MSu))) Then 
Call SuThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gSuAhChiFu,gSuAhFdcCha,gSuAhFuGlu,gSuSuChiChi,gSucHpmChiFd,             & 
& gSuhhChiFu,gSuChiGluSu,gSuSdChicCha,gSuFdFucSd,gSuhhFdcCha,gSuSuFdcFd,gSucHpmChaFu,    & 
& gSuFuFucSu,gSucChaFuHpm,gSuSdFucFd,gSuSuFucFu,gSucHpmFdGlu,gSuhhFuGlu,gSuGluGluSu,     & 
& gSuGluSdcCha,gSuSuChacCha,epsI,deltaM,.False.,gTSu,gPSu(:,193:2880),BRSu(:,193:2880))

Else 
Call SuThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gSuAhChiFu,gSuAhFdcCha,gSuAhFuGlu,gSuSuChiChi,gSucHpmChiFd,             & 
& gSuhhChiFu,gSuChiGluSu,gSuSdChicCha,gSuFdFucSd,gSuhhFdcCha,gSuSuFdcFd,gSucHpmChaFu,    & 
& gSuFuFucSu,gSucChaFuHpm,gSuSdFucFd,gSuSuFucFu,gSucHpmFdGlu,gSuhhFuGlu,gSuGluGluSu,     & 
& gSuGluSdcCha,gSuSuChacCha,epsI,deltaM,.True.,gTSu,gPSu(:,193:2880),BRSu(:,193:2880))

End If 

End If 
Else 
Call SuThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gSuAhChiFu,gSuAhFdcCha,gSuAhFuGlu,gSuSuChiChi,gSucHpmChiFd,             & 
& gSuhhChiFu,gSuChiGluSu,gSuSdChicCha,gSuFdFucSd,gSuhhFdcCha,gSuSuFdcFd,gSucHpmChaFu,    & 
& gSuFuFucSu,gSucChaFuHpm,gSuSdFucFd,gSuSuFucFu,gSucHpmFdGlu,gSuhhFuGlu,gSuGluGluSu,     & 
& gSuGluSdcCha,gSuSuChacCha,epsI,deltaM,.False.,gTSu,gPSu(:,193:2880),BRSu(:,193:2880))

End If 
Do i1=1,6
gTSu(i1) =Sum(gPSu(i1,:)) 
If (gTSu(i1).Gt.0._dp) BRSu(i1,: ) =gPSu(i1,:)/gTSu(i1) 
End Do 
 

gPhh = 0._dp 
gThh = 0._dp 
BRhh = 0._dp 

Call hhTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,            & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,             & 
& pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,              & 
& Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,              & 
& M2,M3,vd,vu,vL,vR,gPhh,gThh,BRhh)


gPAh = 0._dp 
gTAh = 0._dp 
BRAh = 0._dp 
Call AhTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,            & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,             & 
& pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,              & 
& Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,              & 
& M2,M3,vd,vu,vL,vR,gPAh,gTAh,BRAh)


! Set Goldstone Widhts 
gTAh(1)=gTVZ


gPHpm = 0._dp 
gTHpm = 0._dp 
BRHpm = 0._dp 
Call HpmTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,               & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,              & 
& MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,             & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,vd,vu,vL,vR,gPHpm,gTHpm,BRHpm)

! Set Goldstone Widhts 
gTHpm(1)=gTVWm


gPGlu = 0._dp 
gTGlu = 0._dp 
BRGlu = 0._dp 
If (.Not.CTBD) Then 
Call GluTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,               & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,              & 
& MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,             & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,vd,vu,vL,vR,gPGlu(:,1:36),gTGlu,BRGlu(:,1:36))

If (Enable3BDecaysF) Then 
If (gTGlu.Lt.fac3*Abs(MGlu)) Then 
Call GluThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTSd,gTSu,gGluFdcFdChi,gGluFdcFucCha,gGluFucFuChi,epsI,deltaM,          & 
& .False.,gTGlu,gPGlu(:,37:261),BRGlu(:,37:261))

Else 
Call GluThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTSd,gTSu,gGluFdcFdChi,gGluFdcFucCha,gGluFucFuChi,epsI,deltaM,          & 
& .True.,gTGlu,gPGlu(:,37:261),BRGlu(:,37:261))

End If 

End If 
Else 
Call GluThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTSd,gTSu,gGluFdcFdChi,gGluFdcFucCha,gGluFucFuChi,epsI,deltaM,          & 
& .False.,gTGlu,gPGlu(:,37:261),BRGlu(:,37:261))

End If 
Do i1=1,1
gTGlu =Sum(gPGlu(i1,:)) 
If (gTGlu.Gt.0._dp) BRGlu(i1,: ) =gPGlu(i1,:)/gTGlu 
End Do 
 

gPFu = 0._dp 
gTFu = 0._dp 
BRFu = 0._dp 
Call FuTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,            & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,             & 
& pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,              & 
& Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,              & 
& M2,M3,vd,vu,vL,vR,gPFu,gTFu,BRFu)

gPCha = 0._dp 
gTCha = 0._dp 
BRCha = 0._dp 
If (.Not.CTBD) Then 
Call ChaTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,               & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,              & 
& MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,             & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,vd,vu,vL,vR,gPCha(:,1:196),gTCha,BRCha(:,1:196))

If (Enable3BDecaysF) Then 
If (MaxVal(gTCha).Lt.MaxVal(fac3*Abs(MCha))) Then 
Call ChaThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChaChacChaCha,gChaChaChiChi,      & 
& gChaChacFdFd,gChaChacFuFu,gChaChicFuFd,gChaFdcFuGlu,epsI,deltaM,.False.,               & 
& gTCha,gPCha(:,197:1010),BRCha(:,197:1010))

Else 
Call ChaThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChaChacChaCha,gChaChaChiChi,      & 
& gChaChacFdFd,gChaChacFuFu,gChaChicFuFd,gChaFdcFuGlu,epsI,deltaM,.True.,gTCha,          & 
& gPCha(:,197:1010),BRCha(:,197:1010))

End If 
 
End If 
Else 
Call ChaThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChaChacChaCha,gChaChaChiChi,      & 
& gChaChacFdFd,gChaChacFuFu,gChaChicFuFd,gChaFdcFuGlu,epsI,deltaM,.False.,               & 
& gTCha,gPCha(:,197:1010),BRCha(:,197:1010))

End If 
Do i1=1,5
gTCha(i1) =Sum(gPCha(i1,:)) 
If (gTCha(i1).Gt.0._dp) BRCha(i1,: ) =gPCha(i1,:)/gTCha(i1) 
End Do 
 
gPChi = 0._dp 
gTChi = 0._dp 
BRChi = 0._dp 
If (.Not.CTBD) Then 
Call ChiTwoBodyDecay(-1,DeltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,               & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,              & 
& MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,             & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,vd,vu,vL,vR,gPChi(:,1:236),gTChi,BRChi(:,1:236))


If (Enable3BDecaysF) Then 
If (MaxVal(gTChi).Lt.MaxVal(fac3*Abs(MChi))) Then 
Call ChiThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChiChicChaCha,gChiChiChiChi,      & 
& gChiChicFdFd,gChiChicFuFu,gChiChacFdFu,gChiFdcFdGlu,gChiFucFuGlu,epsI,deltaM,          & 
& .False.,gTChi,gPChi(:,237:1729),BRChi(:,237:1729))

Else 
Call ChiThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChiChicChaCha,gChiChiChiChi,      & 
& gChiChicFdFd,gChiChicFuFu,gChiChacFdFu,gChiFdcFdGlu,gChiFucFuGlu,epsI,deltaM,          & 
& .True.,gTChi,gPChi(:,237:1729),BRChi(:,237:1729))

End If 
 
End If 
Else 
Call ChiThreeBodyDecay(-1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChiChicChaCha,gChiChiChiChi,      & 
& gChiChicFdFd,gChiChicFuFu,gChiChacFdFu,gChiFdcFdGlu,gChiFucFuGlu,epsI,deltaM,          & 
& .False.,gTChi,gPChi(:,237:1729),BRChi(:,237:1729))

End If 
Do i1=1,10
gTChi(i1) =Sum(gPChi(i1,:)) 
If (gTChi(i1).Gt.0._dp) BRChi(i1,: ) =gPChi(i1,:)/gTChi(i1) 
End Do 
 

Iname = Iname - 1 
 
End Subroutine CalculateBR 
End Module BranchingRatios_munuSSM3G 
 
