! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:25 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module Chi3Decays_munuSSM3G 
 
Use Control 
Use CouplingsForDecays_munuSSM3G 
Use ThreeBodyPhaseSpace 
 
Contains 
 
Subroutine ChiThreeBodyDecay(n_in,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,            & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,              & 
& MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,             & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,vd,vu,vL,vR,gTAh,gThh,gTHpm,gTSd,gTSu,gTVWm,gTVZ,gChiChicChaCha,              & 
& gChiChiChiChi,gChiChicFdFd,gChiChicFuFu,gChiChacFdFu,gChiFdcFdGlu,gChiFucFuGlu,        & 
& epsI,deltaM,CheckRealStates,gT,gPartial,BR)

Implicit None 
 
Real(dp),Intent(in) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Complex(dp) :: cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),  & 
& cplcChaChaVZL(5,5),cplcChaChaVZR(5,5),cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),   & 
& cplcChaChiVWmL(5,10),cplcChaChiVWmR(5,10),cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),     & 
& cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),       & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),               & 
& cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),           & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),         & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),           & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),cplChaFucSdL(5,3,6),& 
& cplChaFucSdR(5,3,6),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),cplChiChacVWmL(10,5),& 
& cplChiChacVWmR(10,5),cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),& 
& cplChiChihhR(10,10,8),cplChiChiVZL(10,10),cplChiChiVZR(10,10),cplChiFdcSdL(10,3,6),    & 
& cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),cplGluFdcSdL(3,6),      & 
& cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6)

Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(inout) :: gChiChicChaCha(10,10,5,5),gChiChiChiChi(10,10,10,10),gChiChicFdFd(10,10,3,3),         & 
& gChiChicFuFu(10,10,3,3),gChiChacFdFu(10,5,3,3),gChiFdcFdGlu(10,3,3,1),gChiFucFuGlu(10,3,3,1)

Real(dp),Intent(in) :: gTAh(8),gThh(8),gTHpm(8),gTSd(6),gTSu(6),gTVWm,gTVZ

Real(dp) :: gChiChicChaChai(10,5,5),gChiChiChiChii(10,10,10),gChiChicFdFdi(10,3,3),               & 
& gChiChicFuFui(10,3,3),gChiChacFdFui(5,3,3),gChiFdcFdGlui(3,3,1),gChiFucFuGlui(3,3,1)

Real(dp) :: gTAhtemp(8),gThhtemp(8),gTHpmtemp(8),gTSdtemp(6),gTSutemp(6),gTVWmtemp,gTVZtemp
Integer :: NVs,NVst,NSs,NVVst,NVVss,NVSss,NVSst,NSSss,NSSst
Complex(dp), Allocatable :: IntegralVVst(:,:),IntegralVSss(:,:),IntegralVSst(:,:),IntegralSSss(:,:)               & 
& ,IntegralSSst(:,:)
Real(dp), Allocatable :: IntegralVs(:,:),IntegralVst(:,:),IntegralSs(:,:),IntegralVVss(:,:)
Real(dp), Intent(inout), Optional :: BR(:,:), gPartial(:,:) 
Real(dp), Intent(inout) :: gT(:) 
Integer, Intent(in) :: n_in 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Intent(in) ::  CheckRealStates 
Integer :: i_start, i_end, i_run, n_out, n_length, gt1, gt2, gt3, i1 
Logical :: check 
Iname = Iname +1 
NameOfUnit(Iname) = 'ChiThreeBodyDecay' 
 
Allocate( IntegralVs(25000,8) ) 
Allocate( IntegralVst(25000,12) ) 
Allocate( IntegralSs(500000,10) ) 
Allocate( IntegralVVst(25000,12) ) 
Allocate( IntegralVVss(500000,12) ) 
Allocate( IntegralVSss(500000,12) ) 
Allocate( IntegralVSst(500000,16) ) 
Allocate( IntegralSSss(500000,12) ) 
Allocate( IntegralSSst(500000,16) ) 

 
If (CheckRealStates) Then 
gTAhtemp = 0._dp 
gThhtemp = 0._dp 
gTHpmtemp = 0._dp 
gTSdtemp = 0._dp 
gTSutemp = 0._dp 
gTVWmtemp = 0._dp 
gTVZtemp = 0._dp 
Else 
gTAhtemp = gTAh 
gThhtemp = gThh 
gTHpmtemp = gTHpm 
gTSdtemp = gTSd 
gTSutemp = gTSu 
gTVWmtemp = gTVWm 
gTVZtemp = gTVZ 
End If 
 
check=CheckRealStates 

 
If (n_in.Lt.0) Then 
 i_start = 1 
 i_end = 10 
 Else If ( (n_in.Ge.1).And.(n_in.Le. 10) ) Then 
 i_start = n_in 
 i_end = n_in 
Else 
 If (ErrorLevel.Ge.-1) Then 
   Write (ErrCan, *) 'Problem in subroutine'//NameOfUnit(Iname) 
   Write (ErrCan, *) 'Value of n_in out of range, (n_in,10) = ',n_in,10 
 End If 
 
 If (ErrorLevel.Gt.0) then
   Call TerminateProgram 
   return
 endif 

 If (Present(BR)) BR = 0._dp 
 Iname = Iname - 1 
 Return 
End If 
 
Do i_run = i_start, i_end 
 
Call CouplingsFor_Chi_decays_3B(MChi(i_run),i_run,MAh,MAh2,MCha,MCha2,MChi,           & 
& MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,              & 
& MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,             & 
& g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,            & 
& me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,          & 
& cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,cplcChaChiHpmR,               & 
& cplcChaChiVWmL,cplcChaChiVWmR,cplcFdChaSuL,cplcFdChaSuR,cplcFdChiSdL,cplcFdChiSdR,     & 
& cplcFdFdAhL,cplcFdFdAhR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,cplcFdFdVZR,               & 
& cplcFdFuHpmL,cplcFdFuHpmR,cplcFdFuVWmL,cplcFdFuVWmR,cplcFdGluSdL,cplcFdGluSdR,         & 
& cplcFuChiSuL,cplcFuChiSuR,cplcFuFuAhL,cplcFuFuAhR,cplcFuFuhhL,cplcFuFuhhR,             & 
& cplcFuFuVZL,cplcFuFuVZR,cplcFuGluSuL,cplcFuGluSuR,cplChaFucSdL,cplChaFucSdR,           & 
& cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,cplChiChacVWmR,cplChiChiAhL,              & 
& cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,cplChiFdcSdL,         & 
& cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,         & 
& cplGluFucSuR,deltaM)

IntegralVs = 0._dp 
NVs = 0  
IntegralVst = 0._dp 
NVst = 0  
IntegralSs = 0._dp 
NSs = 0  
IntegralVVst = 0._dp 
NVVst = 0  
IntegralVVss = 0._dp 
NVVss = 0  
IntegralVSss = 0._dp 
NVSss = 0  
IntegralVSst = 0._dp 
NVSst = 0  
IntegralSSss = 0._dp 
NSSss = 0  
IntegralSSst = 0._dp 
NSSst = 0  

 
gChiChicChaChai = 0._dp 
Call ChiToChicChaCha(i_run,MChi,MCha,MVZ,MVWm,MHpm,MAh,Mhh,cplcChaChaAhL,             & 
& cplcChaChaAhR,cplcChaChahhL,cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,  & 
& cplcChaChiHpmR,cplcChaChiVWmL,cplcChaChiVWmR,cplChiChacHpmL,cplChiChacHpmR,            & 
& cplChiChacVWmL,cplChiChacVWmR,cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,     & 
& cplChiChiVZL,cplChiChiVZR,IntegralSs,IntegralSSss,IntegralSSst,IntegralVs,             & 
& IntegralVSss,IntegralVSst,IntegralVVss,IntegralVVst,NSs,NSSss,NSSst,NVs,               & 
& NVSss,NVSst,NVVss,NVVst,gTVZtemp,gTVWmtemp,gTHpmtemp,gTAhtemp,gThhtemp,deltaM,         & 
& epsI,check,gChiChicChaChai)

gChiChicChaCha(i_run,:,:,:) = gChiChicChaChai 
gT(i_run) = gT(i_run) + Sum(gChiChicChaChai) 
 
gChiChiChiChii = 0._dp 
Call ChiToChiChiChi(i_run,MChi,MVZ,MAh,Mhh,cplChiChiAhL,cplChiChiAhR,cplChiChihhL,    & 
& cplChiChihhR,cplChiChiVZL,cplChiChiVZR,IntegralSs,IntegralSSss,IntegralSSst,           & 
& IntegralVs,IntegralVSss,IntegralVSst,IntegralVVss,IntegralVVst,NSs,NSSss,              & 
& NSSst,NVs,NVSss,NVSst,NVVss,NVVst,gTVZtemp,gTAhtemp,gThhtemp,deltaM,epsI,              & 
& check,gChiChiChiChii)

gChiChiChiChi(i_run,:,:,:) = gChiChiChiChii 
gT(i_run) = gT(i_run) + Sum(gChiChiChiChii) 
 
gChiChicFdFdi = 0._dp 
Call ChiToChicFdFd(i_run,MChi,MFd,MVZ,MSd,MAh,Mhh,cplcFdChiSdL,cplcFdChiSdR,          & 
& cplcFdFdAhL,cplcFdFdAhR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,cplcFdFdVZR,               & 
& cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,         & 
& cplChiFdcSdL,cplChiFdcSdR,IntegralSs,IntegralSSss,IntegralSSst,IntegralVs,             & 
& IntegralVSss,IntegralVSst,IntegralVVss,NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,          & 
& gTVZtemp,gTSdtemp,gTAhtemp,gThhtemp,deltaM,epsI,check,gChiChicFdFdi)

gChiChicFdFd(i_run,:,:,:) = gChiChicFdFdi 
gT(i_run) = gT(i_run) + Sum(gChiChicFdFdi) 
 
gChiChicFuFui = 0._dp 
Call ChiToChicFuFu(i_run,MChi,MFu,MVZ,MSu,MAh,Mhh,cplcFuChiSuL,cplcFuChiSuR,          & 
& cplcFuFuAhL,cplcFuFuAhR,cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,         & 
& cplChiFucSuL,cplChiFucSuR,IntegralSs,IntegralSSss,IntegralSSst,IntegralVs,             & 
& IntegralVSss,IntegralVSst,IntegralVVss,NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,          & 
& gTVZtemp,gTSutemp,gTAhtemp,gThhtemp,deltaM,epsI,check,gChiChicFuFui)

gChiChicFuFu(i_run,:,:,:) = gChiChicFuFui 
gT(i_run) = gT(i_run) + Sum(gChiChicFuFui) 
 
gChiChacFdFui = 0._dp 
Call ChiToChacFdFu(i_run,MCha,MFd,MFu,MVWm,MHpm,MSu,MSd,MChi,cplcFdChaSuL,            & 
& cplcFdChaSuR,cplcFdChiSdL,cplcFdChiSdR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdFuVWmL,         & 
& cplcFdFuVWmR,cplChaFucSdL,cplChaFucSdR,cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,   & 
& cplChiChacVWmR,cplChiFucSuL,cplChiFucSuR,IntegralSs,IntegralSSss,IntegralSSst,         & 
& IntegralVs,IntegralVSss,IntegralVSst,IntegralVVss,NSs,NSSss,NSSst,NVs,NVSss,           & 
& NVSst,NVVss,gTVWmtemp,gTHpmtemp,gTSutemp,gTSdtemp,deltaM,epsI,check,gChiChacFdFui)

gChiChacFdFu(i_run,:,:,:) = gChiChacFdFui 
gT(i_run) = gT(i_run) + 2._dp*Sum(gChiChacFdFui) 
 
gChiFdcFdGlui = 0._dp 
Call ChiToFdcFdGlu(i_run,MFd,MGlu,MSd,MChi,cplcFdChiSdL,cplcFdChiSdR,cplcFdGluSdL,    & 
& cplcFdGluSdR,cplChiFdcSdL,cplChiFdcSdR,cplGluFdcSdL,cplGluFdcSdR,IntegralSs,           & 
& IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSdtemp,deltaM,epsI,check,gChiFdcFdGlui)

gChiFdcFdGlu(i_run,:,:,:) = gChiFdcFdGlui 
gT(i_run) = gT(i_run) + Sum(gChiFdcFdGlui) 
 
gChiFucFuGlui = 0._dp 
Call ChiToFucFuGlu(i_run,MFu,MGlu,MSu,MChi,cplcFuChiSuL,cplcFuChiSuR,cplcFuGluSuL,    & 
& cplcFuGluSuR,cplChiFucSuL,cplChiFucSuR,cplGluFucSuL,cplGluFucSuR,IntegralSs,           & 
& IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSutemp,deltaM,epsI,check,gChiFucFuGlui)

gChiFucFuGlu(i_run,:,:,:) = gChiFucFuGlui 
gT(i_run) = gT(i_run) + Sum(gChiFucFuGlui) 
 
End Do 
 

If (Present(gPartial)) Then
Do i1 = i_start, i_end 
 
n_length=1
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
gPartial(i1,n_length)= gChiChicChaCha(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
gPartial(i1,n_length)= gChiChiChiChi(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
gPartial(i1,n_length)= gChiChicFdFd(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
gPartial(i1,n_length)= gChiChicFuFu(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
gPartial(i1,n_length)= 2._dp*gChiChacFdFu(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,1
gPartial(i1,n_length)= gChiFdcFdGlu(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,1
gPartial(i1,n_length)= gChiFucFuGlu(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
If (Present(BR).And.(gT(i1).Gt.0._dp)) Then 
BR(i1,:)=gPartial(i1,:)/gT(i1)
Else If (Present(BR)) Then
BR(i1,:)=0._dp
End If
 
End Do 
End if 
Deallocate( IntegralVs ) 
Deallocate( IntegralVst ) 
Deallocate( IntegralSs ) 
Deallocate( IntegralVVst ) 
Deallocate( IntegralVVss ) 
Deallocate( IntegralVSss ) 
Deallocate( IntegralVSst ) 
Deallocate( IntegralSSss ) 
Deallocate( IntegralSSst ) 
Iname = Iname - 1 
 
End Subroutine ChiThreeBodyDecay
 
 
Subroutine ChiToChicChaCha(iIN,MChi,MCha,MVZ,MVWm,MHpm,MAh,Mhh,cplcChaChaAhL,         & 
& cplcChaChaAhR,cplcChaChahhL,cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,  & 
& cplcChaChiHpmR,cplcChaChiVWmL,cplcChaChiVWmR,cplChiChacHpmL,cplChiChacHpmR,            & 
& cplChiChacVWmL,cplChiChacVWmR,cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,     & 
& cplChiChiVZL,cplChiChiVZR,IntegralSs,IntegralSSss,IntegralSSst,IntegralVs,             & 
& IntegralVSss,IntegralVSst,IntegralVVss,IntegralVVst,NSs,NSSss,NSSst,NVs,               & 
& NVSss,NVSst,NVVss,NVVst,gTVZ,gTVWm,gTHpm,gTAh,gThh,deltaM,epsI,check,g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MChi(10),MCha(5),MVZ,MVWm,MHpm(8),MAh(8),Mhh(8)

Complex(dp),Intent(in) :: cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),  & 
& cplcChaChaVZL(5,5),cplcChaChaVZR(5,5),cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),   & 
& cplcChaChiVWmL(5,10),cplcChaChiVWmR(5,10),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),& 
& cplChiChacVWmL(10,5),cplChiChacVWmR(10,5),cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8), & 
& cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),cplChiChiVZL(10,10),cplChiChiVZR(10,10)

Real(dp),Intent(inout) :: IntegralSs(500000,10),IntegralVs(25000,8),IntegralVVss(500000,12)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16),IntegralVSss(500000,12),              & 
& IntegralVSst(500000,16),IntegralVVst(25000,12)

Real(dp),Intent(inout) :: gTVZ,gTVWm,gTHpm(8),gTAh(8),gThh(8)

Integer, Intent(inout) :: NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,NVVst
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 1225 
Allocate( gSum(10,5,5, Isum) ) 
Allocate( Contribution(10,5,5, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1, iIN-1
      Do gt2=1,5
        Do gt3=1,5
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MCha(gt3))+Abs(MCha(gt2))+Abs(MChi(gt1)))) Then 
!-------------- 
!  VZ 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVZ 
Boson2(2) =gTVZ 
 
Boson4(1) = MVZ 
Boson4(2) =gTVZ 
Boson4(3) = MVZ 
Boson4(4) =gTVZ 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup(2) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup(1) = Conjg(cplChiChiVZR(iIN,gt1)) 
coup(4) = Conjg(cplcChaChaVZL(gt2,gt3)) 
coup(3) = Conjg(cplcChaChaVZR(gt2,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ'



!-------------- 
!  conj[VWm] 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVWm 
Boson2(2) =gTVWm 
 
Boson4(1) = MVWm 
Boson4(2) =gTVWm 
Boson4(3) = MVWm 
Boson4(4) =gTVWm 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MCha(gt3) 
mass(3) = -MCha(gt2) 
mass(4) = MChi(gt1) 
 
coup(2) = Conjg(cplChiChacVWmL(iIN,gt3)) 
coup(1) = Conjg(cplChiChacVWmR(iIN,gt3)) 
coup(4) = Conjg(cplcChaChiVWmL(gt2,gt1)) 
coup(3) = Conjg(cplcChaChiVWmR(gt2,gt1))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[VWm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm]'



!-------------- 
!  VWm 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVWm 
Boson2(2) =gTVWm 
 
Boson4(1) = MVWm 
Boson4(2) =gTVWm 
Boson4(3) = MVWm 
Boson4(4) =gTVWm 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MCha(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MCha(gt3) 
 
coup(2) = Conjg(cplcChaChiVWmL(gt2,iIN)) 
coup(1) = Conjg(cplcChaChiVWmR(gt2,iIN)) 
coup(4) = Conjg(cplChiChacVWmL(gt1,gt3)) 
coup(3) = Conjg(cplChiChacVWmR(gt1,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VWm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='VWm'



!-------------- 
!  conj[Hpm] 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MHpm(i1) 
Boson2(2) =gTHpm(i1) 
 
Boson4(1) = MHpm(i1) 
Boson4(2) =gTHpm(i1) 
Boson4(3) = MHpm(i1) 
Boson4(4) =gTHpm(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MCha(gt3) 
mass(3) = -MCha(gt2) 
mass(4) = MChi(gt1) 
 
coup(2) = Conjg(cplChiChacHpmL(iIN,gt3,i1)) 
coup(1) = Conjg(cplChiChacHpmR(iIN,gt3,i1)) 
coup(4) = Conjg(cplcChaChiHpmL(gt2,gt1,i1)) 
coup(3) = Conjg(cplcChaChiHpmR(gt2,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm]'
      End Do 



!-------------- 
!  Hpm 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MHpm(i1) 
Boson2(2) =gTHpm(i1) 
 
Boson4(1) = MHpm(i1) 
Boson4(2) =gTHpm(i1) 
Boson4(3) = MHpm(i1) 
Boson4(4) =gTHpm(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MCha(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MCha(gt3) 
 
coup(2) = Conjg(cplcChaChiHpmL(gt2,iIN,i1)) 
coup(1) = Conjg(cplcChaChiHpmR(gt2,iIN,i1)) 
coup(4) = Conjg(cplChiChacHpmL(gt1,gt3,i1)) 
coup(3) = Conjg(cplChiChacHpmR(gt1,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Hpm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Hpm'
      End Do 



!-------------- 
!  Ah 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MAh(i1) 
Boson2(2) =gTAh(i1) 
 
Boson4(1) = MAh(i1) 
Boson4(2) =gTAh(i1) 
Boson4(3) = MAh(i1) 
Boson4(4) =gTAh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup(2) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChiAhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcChaChaAhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcChaChaAhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah'
      End Do 



!-------------- 
!  hh 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = Mhh(i1) 
Boson2(2) =gThh(i1) 
 
Boson4(1) = Mhh(i1) 
Boson4(2) =gThh(i1) 
Boson4(3) = Mhh(i1) 
Boson4(4) =gThh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup(2) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChihhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcChaChahhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcChaChahhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='hh'
      End Do 



!-------------- 
!  VZ, conj[VWm] 
!-------------- 
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MVWm 
Boson4(4) = gTVWm 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(4) = Conjg(cplChiChacVWmL(iIN,gt3)) 
coup2(3) = Conjg(cplChiChacVWmR(iIN,gt3))  
coup2(5) = cplcChaChaVZL(gt2,gt3) 
coup2(6) = cplcChaChaVZR(gt2,gt3) 
coup2(8) = Conjg(cplcChaChiVWmL(gt2,gt1)) 
coup2(7) = Conjg(cplcChaChiVWmR(gt2,gt1)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ,conj[VWm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,conj[VWm]'



!-------------- 
!  VZ, VWm 
!-------------- 
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MVWm 
Boson4(4) = gTVWm 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(4) = Conjg(cplcChaChiVWmL(gt2,iIN)) 
coup2(3) = Conjg(cplcChaChiVWmR(gt2,iIN))  
coup2(5) = cplcChaChaVZL(gt2,gt3) 
coup2(6) = cplcChaChaVZR(gt2,gt3) 
coup2(8) = Conjg(cplChiChacVWmL(gt1,gt3)) 
coup2(7) = Conjg(cplChiChacVWmR(gt1,gt3)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ,VWm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,VWm'



!-------------- 
!  VZ, conj[Hpm] 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChacHpmL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChacHpmR(iIN,gt3,i2))  
coup2(5) = cplcChaChaVZL(gt2,gt3) 
coup2(6) = cplcChaChaVZR(gt2,gt3) 
coup2(7) = Conjg(cplcChaChiHpmL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcChaChiHpmR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ,conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,conj[Hpm]'
      End Do 



!-------------- 
!  VZ, Hpm 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplcChaChiHpmL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcChaChiHpmR(gt2,iIN,i2))  
coup2(5) = cplcChaChaVZL(gt2,gt3) 
coup2(6) = cplcChaChaVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChacHpmL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChacHpmR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ,Hpm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Hpm'
      End Do 



!-------------- 
!  VZ, Ah 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcChaChaVZL(gt2,gt3) 
coup2(6) = cplcChaChaVZR(gt2,gt3) 
coup2(7) = Conjg(cplcChaChaAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChaAhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Ah'
      End Do 



!-------------- 
!  VZ, hh 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcChaChaVZL(gt2,gt3) 
coup2(6) = cplcChaChaVZR(gt2,gt3) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VZ,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,hh'
      End Do 



!-------------- 
!  conj[VWm], VWm 
!-------------- 
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MVWm 
Boson4(4) = gTVWm 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(4) = MCha(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChacVWmL(iIN,gt3) 
coup2(2) = cplChiChacVWmR(iIN,gt3) 
coup2(4) = Conjg(cplcChaChiVWmL(gt2,iIN)) 
coup2(3) = Conjg(cplcChaChiVWmR(gt2,iIN))  
coup2(5) = cplcChaChiVWmL(gt2,gt1) 
coup2(6) = cplcChaChiVWmR(gt2,gt1) 
coup2(8) = Conjg(cplChiChacVWmL(gt1,gt3)) 
coup2(7) = Conjg(cplChiChacVWmR(gt1,gt3)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[VWm],VWm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],VWm'



!-------------- 
!  conj[VWm], conj[Hpm] 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt3) 
mass(3) = -MCha(gt2) 
mass(4) = MChi(gt1) 
 
coup2(1) = cplChiChacVWmL(iIN,gt3) 
coup2(2) = cplChiChacVWmR(iIN,gt3) 
coup2(3) = Conjg(cplChiChacHpmL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChacHpmR(iIN,gt3,i2))  
coup2(5) = cplcChaChiVWmL(gt2,gt1) 
coup2(6) = cplcChaChiVWmR(gt2,gt1) 
coup2(7) = Conjg(cplcChaChiHpmL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcChaChiHpmR(gt2,gt1,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[VWm],conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],conj[Hpm]'
      End Do 



!-------------- 
!  conj[VWm], Hpm 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(4) = MCha(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChacVWmL(iIN,gt3) 
coup2(2) = cplChiChacVWmR(iIN,gt3) 
coup2(3) = Conjg(cplcChaChiHpmL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcChaChiHpmR(gt2,iIN,i2))  
coup2(5) = cplcChaChiVWmL(gt2,gt1) 
coup2(6) = cplcChaChiVWmR(gt2,gt1) 
coup2(7) = Conjg(cplChiChacHpmL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChacHpmR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[VWm],Hpm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],Hpm'
      End Do 



!-------------- 
!  conj[VWm], Ah 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt3) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChacVWmL(iIN,gt3) 
coup2(2) = cplChiChacVWmR(iIN,gt3) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcChaChiVWmL(gt2,gt1) 
coup2(6) = cplcChaChiVWmR(gt2,gt1) 
coup2(7) = Conjg(cplcChaChaAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChaAhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[VWm],Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],Ah'
      End Do 



!-------------- 
!  conj[VWm], hh 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt3) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChacVWmL(iIN,gt3) 
coup2(2) = cplChiChacVWmR(iIN,gt3) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcChaChiVWmL(gt2,gt1) 
coup2(6) = cplcChaChiVWmR(gt2,gt1) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[VWm],hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],hh'
      End Do 



!-------------- 
!  VWm, conj[Hpm] 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt3) 
mass(4) = MCha(gt2) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplcChaChiVWmL(gt2,iIN) 
coup2(2) = cplcChaChiVWmR(gt2,iIN) 
coup2(3) = Conjg(cplChiChacHpmL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChacHpmR(iIN,gt3,i2))  
coup2(5) = cplChiChacVWmL(gt1,gt3) 
coup2(6) = cplChiChacVWmR(gt1,gt3) 
coup2(7) = Conjg(cplcChaChiHpmL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcChaChiHpmR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VWm,conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VWm,conj[Hpm]'
      End Do 



!-------------- 
!  VWm, Hpm 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplcChaChiVWmL(gt2,iIN) 
coup2(2) = cplcChaChiVWmR(gt2,iIN) 
coup2(3) = Conjg(cplcChaChiHpmL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcChaChiHpmR(gt2,iIN,i2))  
coup2(5) = cplChiChacVWmL(gt1,gt3) 
coup2(6) = cplChiChacVWmR(gt1,gt3) 
coup2(7) = Conjg(cplChiChacHpmL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChacHpmR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VWm,Hpm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VWm,Hpm'
      End Do 



!-------------- 
!  VWm, Ah 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt2) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplcChaChiVWmL(gt2,iIN) 
coup2(2) = cplcChaChiVWmR(gt2,iIN) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChacVWmL(gt1,gt3) 
coup2(6) = cplChiChacVWmR(gt1,gt3) 
coup2(7) = Conjg(cplcChaChaAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChaAhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VWm,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VWm,Ah'
      End Do 



!-------------- 
!  VWm, hh 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt2) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplcChaChiVWmL(gt2,iIN) 
coup2(2) = cplcChaChiVWmR(gt2,iIN) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChacVWmL(gt1,gt3) 
coup2(6) = cplChiChacVWmR(gt1,gt3) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: VWm,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VWm,hh'
      End Do 



!-------------- 
!  conj[Hpm], conj[Hpm] 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt3) 
mass(3) = -MCha(gt2) 
mass(4) = MChi(gt1) 
 
coup2(1) = cplChiChacHpmL(iIN,gt3,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChacHpmL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChacHpmR(iIN,gt3,i2))  
coup2(5) = cplcChaChiHpmL(gt2,gt1,i1) 
coup2(6) = cplcChaChiHpmR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcChaChiHpmL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcChaChiHpmR(gt2,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[Hpm],conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],conj[Hpm]'
        End Do 
      End Do 



!-------------- 
!  conj[Hpm], Hpm 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(4) = MCha(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChacHpmL(iIN,gt3,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt3,i1) 
coup2(3) = Conjg(cplcChaChiHpmL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcChaChiHpmR(gt2,iIN,i2))  
coup2(5) = cplcChaChiHpmL(gt2,gt1,i1) 
coup2(6) = cplcChaChiHpmR(gt2,gt1,i1) 
coup2(7) = Conjg(cplChiChacHpmL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChacHpmR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[Hpm],Hpm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],Hpm'
        End Do 
      End Do 



!-------------- 
!  conj[Hpm], Ah 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt3) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChacHpmL(iIN,gt3,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcChaChiHpmL(gt2,gt1,i1) 
coup2(6) = cplcChaChiHpmR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcChaChaAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChaAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[Hpm],Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],Ah'
        End Do 
      End Do 



!-------------- 
!  conj[Hpm], hh 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt3) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplChiChacHpmL(iIN,gt3,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcChaChiHpmL(gt2,gt1,i1) 
coup2(6) = cplcChaChiHpmR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: conj[Hpm],hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],hh'
        End Do 
      End Do 



!-------------- 
!  Hpm, Hpm 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplcChaChiHpmL(gt2,iIN,i1) 
coup2(2) = cplcChaChiHpmR(gt2,iIN,i1) 
coup2(3) = Conjg(cplcChaChiHpmL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcChaChiHpmR(gt2,iIN,i2))  
coup2(5) = cplChiChacHpmL(gt1,gt3,i1) 
coup2(6) = cplChiChacHpmR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChacHpmL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChacHpmR(gt1,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Hpm,Hpm" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Hpm,Hpm'
        End Do 
      End Do 



!-------------- 
!  Hpm, Ah 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt2) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplcChaChiHpmL(gt2,iIN,i1) 
coup2(2) = cplcChaChiHpmR(gt2,iIN,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChacHpmL(gt1,gt3,i1) 
coup2(6) = cplChiChacHpmR(gt1,gt3,i1) 
coup2(7) = Conjg(cplcChaChaAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChaAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Hpm,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Hpm,Ah'
        End Do 
      End Do 



!-------------- 
!  Hpm, hh 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MCha(gt2) 
mass(3) = -MCha(gt2) 
 
coup2(1) = cplcChaChiHpmL(gt2,iIN,i1) 
coup2(2) = cplcChaChiHpmR(gt2,iIN,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChacHpmL(gt1,gt3,i1) 
coup2(6) = cplChiChacHpmR(gt1,gt3,i1) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Hpm,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Hpm,hh'
        End Do 
      End Do 



!-------------- 
!  Ah, Ah 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcChaChaAhL(gt2,gt3,i1) 
coup2(6) = cplcChaChaAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcChaChaAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChaAhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Ah,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,Ah'
        End Do 
      End Do 



!-------------- 
!  Ah, hh 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcChaChaAhL(gt2,gt3,i1) 
coup2(6) = cplcChaChaAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: Ah,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,hh'
        End Do 
      End Do 



!-------------- 
!  hh, hh 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = Mhh(i1) 
Boson4(2) = gThh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MCha(gt2) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcChaChahhL(gt2,gt3,i1) 
coup2(6) = cplcChaChahhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcChaChahhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcChaChahhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cCha Cha Propagator: hh,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='hh,hh'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,gt3,:)= 0._dp  
End If 
       End Do 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1, iIN-1
      Do gt2=1,5
        Do gt3=1,5
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:1225))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToChicChaCha 
 
 
Subroutine ChiToChiChiChi(iIN,MChi,MVZ,MAh,Mhh,cplChiChiAhL,cplChiChiAhR,             & 
& cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,IntegralSs,IntegralSSss,           & 
& IntegralSSst,IntegralVs,IntegralVSss,IntegralVSst,IntegralVVss,IntegralVVst,           & 
& NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,NVVst,gTVZ,gTAh,gThh,deltaM,epsI,check,          & 
& g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MChi(10),MVZ,MAh(8),Mhh(8)

Complex(dp),Intent(in) :: cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),& 
& cplChiChiVZL(10,10),cplChiChiVZR(10,10)

Real(dp),Intent(inout) :: IntegralSs(500000,10),IntegralVs(25000,8),IntegralVVss(500000,12)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16),IntegralVSss(500000,12),              & 
& IntegralVSst(500000,16),IntegralVVst(25000,12)

Real(dp),Intent(inout) :: gTVZ,gTAh(8),gThh(8)

Integer, Intent(inout) :: NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,NVVst
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 289 
Allocate( gSum(10,10,10, Isum) ) 
Allocate( Contribution(10,10,10, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1, iIN-1
      Do gt2=gt1, iIN-1
        Do gt3=gt2, iIN-1
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MChi(gt3))+Abs(MChi(gt2))+Abs(MChi(gt1)))) Then 
!-------------- 
!  VZ 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVZ 
Boson2(2) =gTVZ 
 
Boson4(1) = MVZ 
Boson4(2) =gTVZ 
Boson4(3) = MVZ 
Boson4(4) =gTVZ 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup(1) = Conjg(cplChiChiVZR(iIN,gt1)) 
coup(4) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup(3) = Conjg(cplChiChiVZR(gt2,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
If (gt1.Eq.gt2) Then 
resR=resR/2._dp 
End If
resR= 1*resR ! color factor 
resS = resS + resR 
 
 mass(2) = MChi(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplChiChiVZL(iIN,gt2)) 
coup(1) = Conjg(cplChiChiVZR(iIN,gt2)) 
coup(4) = Conjg(cplChiChiVZL(gt1,gt3)) 
coup(3) = Conjg(cplChiChiVZR(gt1,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
If (gt1.Eq.gt2) Then 
resR=resR/2._dp 
End If
resR= 1*resR ! color factor 
resS = resS + resR 
 
 mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(4) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup2(3) = Conjg(cplChiChiVZR(iIN,gt1))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(8) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup2(7) = Conjg(cplChiChiVZR(gt2,gt3)) 
Call IntegrateGaugeSS(Boson4, mass, coup2, deltaM, epsI,IntegralVVss,NVVss, resR, check) 
If (resR.ne.resR) resR = 0._dp
resC = -2._dp*resR 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(4) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup2(3) = Conjg(cplChiChiVZR(iIN,gt1))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(8) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup2(7) = Conjg(cplChiChiVZR(gt2,gt3)) 
Call IntegrateGaugeSS(Boson4, mass, coup2, deltaM, epsI,IntegralVVss,NVVss, resR, check) 
If (resR.ne.resR) resR = 0._dp
resC = -2._dp*resR 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt3) 
coup2(2) = cplChiChiVZR(iIN,gt3) 
coup2(4) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup2(3) = Conjg(cplChiChiVZR(iIN,gt1))  
coup2(5) = cplChiChiVZL(gt1,gt2) 
coup2(6) = cplChiChiVZR(gt1,gt2) 
coup2(8) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup2(7) = Conjg(cplChiChiVZR(gt2,gt3)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt3) 
coup2(2) = cplChiChiVZR(iIN,gt3) 
coup2(4) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup2(3) = Conjg(cplChiChiVZR(iIN,gt1))  
coup2(5) = cplChiChiVZL(gt1,gt2) 
coup2(6) = cplChiChiVZR(gt1,gt2) 
coup2(8) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup2(7) = Conjg(cplChiChiVZR(gt2,gt3)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt2) 
coup2(2) = cplChiChiVZR(iIN,gt2) 
coup2(4) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup2(3) = Conjg(cplChiChiVZR(iIN,gt1))  
coup2(5) = cplChiChiVZL(gt1,gt3) 
coup2(6) = cplChiChiVZR(gt1,gt3) 
coup2(8) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup2(7) = Conjg(cplChiChiVZR(gt2,gt3)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt2) 
coup2(2) = cplChiChiVZR(iIN,gt2) 
coup2(4) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup2(3) = Conjg(cplChiChiVZR(iIN,gt1))  
coup2(5) = cplChiChiVZL(gt1,gt3) 
coup2(6) = cplChiChiVZR(gt1,gt3) 
coup2(8) = Conjg(cplChiChiVZL(gt2,gt3)) 
coup2(7) = Conjg(cplChiChiVZR(gt2,gt3)) 
Call IntegrateGaugeST(Boson4, mass, coup2, deltaM, epsI,IntegralVVst,NVVst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: VZ" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ'



!-------------- 
!  Ah 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MAh(i1) 
Boson2(2) =gTAh(i1) 
 
Boson4(1) = MAh(i1) 
Boson4(2) =gTAh(i1) 
Boson4(3) = MAh(i1) 
Boson4(4) =gTAh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChiAhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplChiChiAhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
If (gt1.Eq.gt2) Then 
resR=resR/2._dp 
End If
resR= 1*resR ! color factor 
resS = resS + resR 
 
 mass(2) = MChi(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplChiChiAhL(iIN,gt2,i1)) 
coup(1) = Conjg(cplChiChiAhR(iIN,gt2,i1)) 
coup(4) = Conjg(cplChiChiAhL(gt1,gt3,i1)) 
coup(3) = Conjg(cplChiChiAhR(gt1,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
If (gt1.Eq.gt2) Then 
resR=resR/2._dp 
End If
resR= 1*resR ! color factor 
resS = resS + resR 
 
 mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i1))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i1)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i1))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i1)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt3,i1) 
coup2(2) = cplChiChiAhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i1))  
coup2(5) = cplChiChiAhL(gt1,gt2,i1) 
coup2(6) = cplChiChiAhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt3,i1) 
coup2(2) = cplChiChiAhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i1))  
coup2(5) = cplChiChiAhL(gt1,gt2,i1) 
coup2(6) = cplChiChiAhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt2,i1) 
coup2(2) = cplChiChiAhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i1))  
coup2(5) = cplChiChiAhL(gt1,gt3,i1) 
coup2(6) = cplChiChiAhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt2,i1) 
coup2(2) = cplChiChiAhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i1))  
coup2(5) = cplChiChiAhL(gt1,gt3,i1) 
coup2(6) = cplChiChiAhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah'
      End Do 



!-------------- 
!  hh 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = Mhh(i1) 
Boson2(2) =gThh(i1) 
 
Boson4(1) = Mhh(i1) 
Boson4(2) =gThh(i1) 
Boson4(3) = Mhh(i1) 
Boson4(4) =gThh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChihhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplChiChihhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
If (gt1.Eq.gt2) Then 
resR=resR/2._dp 
End If
resR= 1*resR ! color factor 
resS = resS + resR 
 
 mass(2) = MChi(gt2) 
mass(3) = -MChi(gt1) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplChiChihhL(iIN,gt2,i1)) 
coup(1) = Conjg(cplChiChihhR(iIN,gt2,i1)) 
coup(4) = Conjg(cplChiChihhL(gt1,gt3,i1)) 
coup(3) = Conjg(cplChiChihhR(gt1,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
If (gt1.Eq.gt2) Then 
resR=resR/2._dp 
End If
resR= 1*resR ! color factor 
resS = resS + resR 
 
 mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i1))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i1)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i1))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i1)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChihhL(iIN,gt3,i1) 
coup2(2) = cplChiChihhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i1))  
coup2(5) = cplChiChihhL(gt1,gt2,i1) 
coup2(6) = cplChiChihhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChihhL(iIN,gt3,i1) 
coup2(2) = cplChiChihhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i1))  
coup2(5) = cplChiChihhL(gt1,gt2,i1) 
coup2(6) = cplChiChihhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt2,i1) 
coup2(2) = cplChiChihhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i1))  
coup2(5) = cplChiChihhL(gt1,gt3,i1) 
coup2(6) = cplChiChihhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt2,i1) 
coup2(2) = cplChiChihhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i1))  
coup2(5) = cplChiChihhL(gt1,gt3,i1) 
coup2(6) = cplChiChihhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i1)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i1)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='hh'
      End Do 



!-------------- 
!  VZ, Ah 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt2,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt3,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt3,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt2,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt2,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt3) 
coup2(2) = cplChiChiVZR(iIN,gt3) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChiVZL(gt1,gt2) 
coup2(6) = cplChiChiVZR(gt1,gt2) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiVZL(iIN,gt3) 
coup2(2) = cplChiChiVZR(iIN,gt3) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt2,i2))  
coup2(5) = cplChiChiVZL(gt1,gt2) 
coup2(6) = cplChiChiVZR(gt1,gt2) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt2) 
coup2(2) = cplChiChiVZR(iIN,gt2) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChiVZL(gt1,gt3) 
coup2(6) = cplChiChiVZR(gt1,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiVZL(iIN,gt2) 
coup2(2) = cplChiChiVZR(iIN,gt2) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt3,i2))  
coup2(5) = cplChiChiVZL(gt1,gt3) 
coup2(6) = cplChiChiVZR(gt1,gt3) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt2,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: VZ,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Ah'
      End Do 



!-------------- 
!  VZ, hh 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt2,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChiVZL(gt2,gt3) 
coup2(6) = cplChiChiVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt3) 
coup2(2) = cplChiChiVZR(iIN,gt3) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChiVZL(gt1,gt2) 
coup2(6) = cplChiChiVZR(gt1,gt2) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiVZL(iIN,gt3) 
coup2(2) = cplChiChiVZR(iIN,gt3) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChiVZL(gt1,gt2) 
coup2(6) = cplChiChiVZR(gt1,gt2) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt2) 
coup2(2) = cplChiChiVZR(iIN,gt2) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChiVZL(gt1,gt3) 
coup2(6) = cplChiChiVZR(gt1,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiVZL(iIN,gt2) 
coup2(2) = cplChiChiVZR(iIN,gt2) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChiVZL(gt1,gt3) 
coup2(6) = cplChiChiVZR(gt1,gt3) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt2,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: VZ,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,hh'
      End Do 



!-------------- 
!  Ah, Ah 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt2,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt3,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt3,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt2,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt2,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt3,i1) 
coup2(2) = cplChiChiAhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChiAhL(gt1,gt2,i1) 
coup2(6) = cplChiChiAhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiAhL(iIN,gt3,i1) 
coup2(2) = cplChiChiAhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt2,i2))  
coup2(5) = cplChiChiAhL(gt1,gt2,i1) 
coup2(6) = cplChiChiAhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt2,i1) 
coup2(2) = cplChiChiAhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiChiAhL(gt1,gt3,i1) 
coup2(6) = cplChiChiAhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiAhL(iIN,gt2,i1) 
coup2(2) = cplChiChiAhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt3,i2))  
coup2(5) = cplChiChiAhL(gt1,gt3,i1) 
coup2(6) = cplChiChiAhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChiAhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChiAhR(gt1,gt2,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: Ah,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,Ah'
        End Do 
      End Do 



!-------------- 
!  Ah, hh 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt2,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChiAhL(gt2,gt3,i1) 
coup2(6) = cplChiChiAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChiAhL(iIN,gt3,i1) 
coup2(2) = cplChiChiAhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChiAhL(gt1,gt2,i1) 
coup2(6) = cplChiChiAhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiAhL(iIN,gt3,i1) 
coup2(2) = cplChiChiAhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChiAhL(gt1,gt2,i1) 
coup2(6) = cplChiChiAhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt2,i1) 
coup2(2) = cplChiChiAhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChiAhL(gt1,gt3,i1) 
coup2(6) = cplChiChiAhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChiAhL(iIN,gt2,i1) 
coup2(2) = cplChiChiAhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChiAhL(gt1,gt3,i1) 
coup2(6) = cplChiChiAhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt2,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: Ah,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,hh'
        End Do 
      End Do 



!-------------- 
!  hh, hh 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = Mhh(i1) 
Boson4(2) = gThh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MChi(gt2) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt2,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChihhL(gt2,gt3,i1) 
coup2(6) = cplChiChihhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt2) 
 
coup2(1) = cplChiChihhL(iIN,gt3,i1) 
coup2(2) = cplChiChihhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChihhL(gt1,gt2,i1) 
coup2(6) = cplChiChihhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt2) 
mass(4) = MChi(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChihhL(iIN,gt3,i1) 
coup2(2) = cplChiChihhR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt2,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt2,i2))  
coup2(5) = cplChiChihhL(gt1,gt2,i1) 
coup2(6) = cplChiChihhR(gt1,gt2,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt1) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt2,i1) 
coup2(2) = cplChiChihhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiChihhL(gt1,gt3,i1) 
coup2(6) = cplChiChihhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
mass(2) = MChi(gt3) 
mass(4) = MChi(gt2) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiChihhL(iIN,gt2,i1) 
coup2(2) = cplChiChihhR(iIN,gt2,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt3,i2))  
coup2(5) = cplChiChihhL(gt1,gt3,i1) 
coup2(6) = cplChiChihhR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiChihhL(gt1,gt2,i2)) 
coup2(8) = Conjg(cplChiChihhR(gt1,gt2,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1*resC ! color factor 
If ((gt1.Eq.gt2).And.(gt2.Eq.gt3)) Then 
resC=resC/6._dp
Else If ((gt1.Eq.gt2).Or.(gt2.Eq.gt3).Or.(gt1.Eq.gt3)) Then
resC=resC/2._dp 
End If
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi Chi Chi Propagator: hh,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='hh,hh'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,gt3,:)= 0._dp  
End If 
       End Do 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1, iIN-1
      Do gt2=gt1, iIN-1
        Do gt3=gt2, iIN-1
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:289))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToChiChiChi 
 
 
Subroutine ChiToChicFdFd(iIN,MChi,MFd,MVZ,MSd,MAh,Mhh,cplcFdChiSdL,cplcFdChiSdR,      & 
& cplcFdFdAhL,cplcFdFdAhR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,cplcFdFdVZR,               & 
& cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,         & 
& cplChiFdcSdL,cplChiFdcSdR,IntegralSs,IntegralSSss,IntegralSSst,IntegralVs,             & 
& IntegralVSss,IntegralVSst,IntegralVVss,NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,          & 
& gTVZ,gTSd,gTAh,gThh,deltaM,epsI,check,g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MChi(10),MFd(3),MVZ,MSd(6),MAh(8),Mhh(8)

Complex(dp),Intent(in) :: cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),      & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),               & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),& 
& cplChiChiVZL(10,10),cplChiChiVZR(10,10),cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6)

Real(dp),Intent(inout) :: IntegralSs(500000,10),IntegralVs(25000,8),IntegralVVss(500000,12)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16),IntegralVSss(500000,12),              & 
& IntegralVSst(500000,16)

Real(dp),Intent(inout) :: gTVZ,gTSd(6),gTAh(8),gThh(8)

Integer, Intent(inout) :: NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 841 
Allocate( gSum(10,3,3, Isum) ) 
Allocate( Contribution(10,3,3, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1, iIN-1
      Do gt2=1,3
        Do gt3=1,3
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MFd(gt3))+Abs(MFd(gt2))+Abs(MChi(gt1)))) Then 
!-------------- 
!  VZ 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVZ 
Boson2(2) =gTVZ 
 
Boson4(1) = MVZ 
Boson4(2) =gTVZ 
Boson4(3) = MVZ 
Boson4(4) =gTVZ 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup(2) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup(1) = Conjg(cplChiChiVZR(iIN,gt1)) 
coup(4) = Conjg(cplcFdFdVZL(gt2,gt3)) 
coup(3) = Conjg(cplcFdFdVZR(gt2,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: VZ" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ'



!-------------- 
!  conj[Sd] 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSd(i1) 
Boson2(2) =gTSd(i1) 
 
Boson4(1) = MSd(i1) 
Boson4(2) =gTSd(i1) 
Boson4(3) = MSd(i1) 
Boson4(4) =gTSd(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFd(gt3) 
mass(3) = -MFd(gt2) 
mass(4) = MChi(gt1) 
 
coup(2) = Conjg(cplChiFdcSdL(iIN,gt3,i1)) 
coup(1) = Conjg(cplChiFdcSdR(iIN,gt3,i1)) 
coup(4) = Conjg(cplcFdChiSdL(gt2,gt1,i1)) 
coup(3) = Conjg(cplcFdChiSdR(gt2,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: conj[Sd]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Sd]'
      End Do 



!-------------- 
!  Sd 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSd(i1) 
Boson2(2) =gTSd(i1) 
 
Boson4(1) = MSd(i1) 
Boson4(2) =gTSd(i1) 
Boson4(3) = MSd(i1) 
Boson4(4) =gTSd(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFd(gt2) 
mass(3) = -MFd(gt3) 
mass(4) = MChi(gt1) 
 
coup(2) = Conjg(cplcFdChiSdL(gt2,iIN,i1)) 
coup(1) = Conjg(cplcFdChiSdR(gt2,iIN,i1)) 
coup(4) = Conjg(cplChiFdcSdL(gt1,gt3,i1)) 
coup(3) = Conjg(cplChiFdcSdR(gt1,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Sd'
      End Do 



!-------------- 
!  Ah 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MAh(i1) 
Boson2(2) =gTAh(i1) 
 
Boson4(1) = MAh(i1) 
Boson4(2) =gTAh(i1) 
Boson4(3) = MAh(i1) 
Boson4(4) =gTAh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup(2) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChiAhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFdFdAhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFdFdAhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah'
      End Do 



!-------------- 
!  hh 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = Mhh(i1) 
Boson2(2) =gThh(i1) 
 
Boson4(1) = Mhh(i1) 
Boson4(2) =gThh(i1) 
Boson4(3) = Mhh(i1) 
Boson4(4) =gThh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup(2) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChihhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFdFdhhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFdFdhhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='hh'
      End Do 



!-------------- 
!  VZ, conj[Sd] 
!-------------- 
  Do i2=1,6
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiFdcSdL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFdcSdR(iIN,gt3,i2))  
coup2(5) = cplcFdFdVZL(gt2,gt3) 
coup2(6) = cplcFdFdVZR(gt2,gt3) 
coup2(7) = Conjg(cplcFdChiSdL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFdChiSdR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: VZ,conj[Sd]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,conj[Sd]'
      End Do 



!-------------- 
!  VZ, Sd 
!-------------- 
  Do i2=1,6
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplcFdFdVZL(gt2,gt3) 
coup2(6) = cplcFdFdVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiFdcSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiFdcSdR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: VZ,Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Sd'
      End Do 



!-------------- 
!  VZ, Ah 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcFdFdVZL(gt2,gt3) 
coup2(6) = cplcFdFdVZR(gt2,gt3) 
coup2(7) = Conjg(cplcFdFdAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdAhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: VZ,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Ah'
      End Do 



!-------------- 
!  VZ, hh 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFdFdVZL(gt2,gt3) 
coup2(6) = cplcFdFdVZR(gt2,gt3) 
coup2(7) = Conjg(cplcFdFdhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdhhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: VZ,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,hh'
      End Do 



!-------------- 
!  conj[Sd], conj[Sd] 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt3) 
mass(3) = -MFd(gt2) 
mass(4) = MChi(gt1) 
 
coup2(1) = cplChiFdcSdL(iIN,gt3,i1) 
coup2(2) = cplChiFdcSdR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiFdcSdL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFdcSdR(iIN,gt3,i2))  
coup2(5) = cplcFdChiSdL(gt2,gt1,i1) 
coup2(6) = cplcFdChiSdR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFdChiSdL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFdChiSdR(gt2,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: conj[Sd],conj[Sd]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Sd],conj[Sd]'
        End Do 
      End Do 



!-------------- 
!  conj[Sd], Sd 
!-------------- 
Do i1=1,6
  Do i2=1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(4) = MFd(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiFdcSdL(iIN,gt3,i1) 
coup2(2) = cplChiFdcSdR(iIN,gt3,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplcFdChiSdL(gt2,gt1,i1) 
coup2(6) = cplcFdChiSdR(gt2,gt1,i1) 
coup2(7) = Conjg(cplChiFdcSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiFdcSdR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: conj[Sd],Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Sd],Sd'
        End Do 
      End Do 



!-------------- 
!  conj[Sd], Ah 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFd(gt3) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplChiFdcSdL(iIN,gt3,i1) 
coup2(2) = cplChiFdcSdR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcFdChiSdL(gt2,gt1,i1) 
coup2(6) = cplcFdChiSdR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFdFdAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: conj[Sd],Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Sd],Ah'
        End Do 
      End Do 



!-------------- 
!  conj[Sd], hh 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFd(gt3) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplChiFdcSdL(iIN,gt3,i1) 
coup2(2) = cplChiFdcSdR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFdChiSdL(gt2,gt1,i1) 
coup2(6) = cplcFdChiSdR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFdFdhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdhhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: conj[Sd],hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Sd],hh'
        End Do 
      End Do 



!-------------- 
!  Sd, Sd 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(3) = -MFd(gt3) 
mass(4) = MChi(gt1) 
 
coup2(1) = cplcFdChiSdL(gt2,iIN,i1) 
coup2(2) = cplcFdChiSdR(gt2,iIN,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplChiFdcSdL(gt1,gt3,i1) 
coup2(6) = cplChiFdcSdR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiFdcSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiFdcSdR(gt1,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Sd,Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Sd,Sd'
        End Do 
      End Do 



!-------------- 
!  Sd, Ah 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFd(gt2) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplcFdChiSdL(gt2,iIN,i1) 
coup2(2) = cplcFdChiSdR(gt2,iIN,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiFdcSdL(gt1,gt3,i1) 
coup2(6) = cplChiFdcSdR(gt1,gt3,i1) 
coup2(7) = Conjg(cplcFdFdAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Sd,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Sd,Ah'
        End Do 
      End Do 



!-------------- 
!  Sd, hh 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFd(gt2) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplcFdChiSdL(gt2,iIN,i1) 
coup2(2) = cplcFdChiSdR(gt2,iIN,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiFdcSdL(gt1,gt3,i1) 
coup2(6) = cplChiFdcSdR(gt1,gt3,i1) 
coup2(7) = Conjg(cplcFdFdhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdhhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Sd,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Sd,hh'
        End Do 
      End Do 



!-------------- 
!  Ah, Ah 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcFdFdAhL(gt2,gt3,i1) 
coup2(6) = cplcFdFdAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFdFdAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdAhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Ah,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,Ah'
        End Do 
      End Do 



!-------------- 
!  Ah, hh 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFdFdAhL(gt2,gt3,i1) 
coup2(6) = cplcFdFdAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFdFdhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdhhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: Ah,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,hh'
        End Do 
      End Do 



!-------------- 
!  hh, hh 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = Mhh(i1) 
Boson4(2) = gThh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFd(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFdFdhhL(gt2,gt3,i1) 
coup2(6) = cplcFdFdhhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFdFdhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFdhhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFd Fd Propagator: hh,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='hh,hh'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,gt3,:)= 0._dp  
End If 
       End Do 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1, iIN-1
      Do gt2=1,3
        Do gt3=1,3
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:841))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToChicFdFd 
 
 
Subroutine ChiToChicFuFu(iIN,MChi,MFu,MVZ,MSu,MAh,Mhh,cplcFuChiSuL,cplcFuChiSuR,      & 
& cplcFuFuAhL,cplcFuFuAhR,cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,         & 
& cplChiFucSuL,cplChiFucSuR,IntegralSs,IntegralSSss,IntegralSSst,IntegralVs,             & 
& IntegralVSss,IntegralVSst,IntegralVVss,NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss,          & 
& gTVZ,gTSu,gTAh,gThh,deltaM,epsI,check,g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MChi(10),MFu(3),MVZ,MSu(6),MAh(8),Mhh(8)

Complex(dp),Intent(in) :: cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),      & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),               & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),& 
& cplChiChiVZL(10,10),cplChiChiVZR(10,10),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6)

Real(dp),Intent(inout) :: IntegralSs(500000,10),IntegralVs(25000,8),IntegralVVss(500000,12)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16),IntegralVSss(500000,12),              & 
& IntegralVSst(500000,16)

Real(dp),Intent(inout) :: gTVZ,gTSu(6),gTAh(8),gThh(8)

Integer, Intent(inout) :: NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 841 
Allocate( gSum(10,3,3, Isum) ) 
Allocate( Contribution(10,3,3, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1, iIN-1
      Do gt2=1,3
        Do gt3=1,3
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MFu(gt3))+Abs(MFu(gt2))+Abs(MChi(gt1)))) Then 
!-------------- 
!  VZ 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVZ 
Boson2(2) =gTVZ 
 
Boson4(1) = MVZ 
Boson4(2) =gTVZ 
Boson4(3) = MVZ 
Boson4(4) =gTVZ 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup(2) = Conjg(cplChiChiVZL(iIN,gt1)) 
coup(1) = Conjg(cplChiChiVZR(iIN,gt1)) 
coup(4) = Conjg(cplcFuFuVZL(gt2,gt3)) 
coup(3) = Conjg(cplcFuFuVZR(gt2,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: VZ" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ'



!-------------- 
!  conj[Su] 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSu(i1) 
Boson2(2) =gTSu(i1) 
 
Boson4(1) = MSu(i1) 
Boson4(2) =gTSu(i1) 
Boson4(3) = MSu(i1) 
Boson4(4) =gTSu(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFu(gt3) 
mass(3) = -MFu(gt2) 
mass(4) = MChi(gt1) 
 
coup(2) = Conjg(cplChiFucSuL(iIN,gt3,i1)) 
coup(1) = Conjg(cplChiFucSuR(iIN,gt3,i1)) 
coup(4) = Conjg(cplcFuChiSuL(gt2,gt1,i1)) 
coup(3) = Conjg(cplcFuChiSuR(gt2,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su]'
      End Do 



!-------------- 
!  Su 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSu(i1) 
Boson2(2) =gTSu(i1) 
 
Boson4(1) = MSu(i1) 
Boson4(2) =gTSu(i1) 
Boson4(3) = MSu(i1) 
Boson4(4) =gTSu(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFu(gt2) 
mass(3) = -MFu(gt3) 
mass(4) = MChi(gt1) 
 
coup(2) = Conjg(cplcFuChiSuL(gt2,iIN,i1)) 
coup(1) = Conjg(cplcFuChiSuR(gt2,iIN,i1)) 
coup(4) = Conjg(cplChiFucSuL(gt1,gt3,i1)) 
coup(3) = Conjg(cplChiFucSuR(gt1,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Su'
      End Do 



!-------------- 
!  Ah 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MAh(i1) 
Boson2(2) =gTAh(i1) 
 
Boson4(1) = MAh(i1) 
Boson4(2) =gTAh(i1) 
Boson4(3) = MAh(i1) 
Boson4(4) =gTAh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup(2) = Conjg(cplChiChiAhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChiAhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFuFuAhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFuFuAhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah'
      End Do 



!-------------- 
!  hh 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = Mhh(i1) 
Boson2(2) =gThh(i1) 
 
Boson4(1) = Mhh(i1) 
Boson4(2) =gThh(i1) 
Boson4(3) = Mhh(i1) 
Boson4(4) =gThh(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup(2) = Conjg(cplChiChihhL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChihhR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFuFuhhL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFuFuhhR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='hh'
      End Do 



!-------------- 
!  VZ, conj[Su] 
!-------------- 
  Do i2=1,6
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt3) 
mass(4) = MChi(gt1) 
mass(3) = -MFu(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiFucSuL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFucSuR(iIN,gt3,i2))  
coup2(5) = cplcFuFuVZL(gt2,gt3) 
coup2(6) = cplcFuFuVZR(gt2,gt3) 
coup2(7) = Conjg(cplcFuChiSuL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFuChiSuR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: VZ,conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,conj[Su]'
      End Do 



!-------------- 
!  VZ, Su 
!-------------- 
  Do i2=1,6
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt2) 
mass(4) = MChi(gt1) 
mass(3) = -MFu(gt2) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplcFuChiSuL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFuChiSuR(gt2,iIN,i2))  
coup2(5) = cplcFuFuVZL(gt2,gt3) 
coup2(6) = cplcFuFuVZR(gt2,gt3) 
coup2(7) = Conjg(cplChiFucSuL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiFucSuR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: VZ,Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Su'
      End Do 



!-------------- 
!  VZ, Ah 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcFuFuVZL(gt2,gt3) 
coup2(6) = cplcFuFuVZR(gt2,gt3) 
coup2(7) = Conjg(cplcFuFuAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuAhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: VZ,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,Ah'
      End Do 



!-------------- 
!  VZ, hh 
!-------------- 
  Do i2=1,8
Boson4(1) = MVZ 
Boson4(2) = gTVZ 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChiVZL(iIN,gt1) 
coup2(2) = cplChiChiVZR(iIN,gt1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFuFuVZL(gt2,gt3) 
coup2(6) = cplcFuFuVZR(gt2,gt3) 
coup2(7) = Conjg(cplcFuFuhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuhhR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: VZ,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='VZ,hh'
      End Do 



!-------------- 
!  conj[Su], conj[Su] 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt3) 
mass(3) = -MFu(gt2) 
mass(4) = MChi(gt1) 
 
coup2(1) = cplChiFucSuL(iIN,gt3,i1) 
coup2(2) = cplChiFucSuR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiFucSuL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFucSuR(iIN,gt3,i2))  
coup2(5) = cplcFuChiSuL(gt2,gt1,i1) 
coup2(6) = cplcFuChiSuR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFuChiSuL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFuChiSuR(gt2,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: conj[Su],conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su],conj[Su]'
        End Do 
      End Do 



!-------------- 
!  conj[Su], Su 
!-------------- 
Do i1=1,6
  Do i2=1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt2) 
mass(4) = MFu(gt3) 
mass(3) = -MChi(gt1) 
 
coup2(1) = cplChiFucSuL(iIN,gt3,i1) 
coup2(2) = cplChiFucSuR(iIN,gt3,i1) 
coup2(3) = Conjg(cplcFuChiSuL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFuChiSuR(gt2,iIN,i2))  
coup2(5) = cplcFuChiSuL(gt2,gt1,i1) 
coup2(6) = cplcFuChiSuR(gt2,gt1,i1) 
coup2(7) = Conjg(cplChiFucSuL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiFucSuR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: conj[Su],Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su],Su'
        End Do 
      End Do 



!-------------- 
!  conj[Su], Ah 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFu(gt3) 
mass(3) = -MFu(gt2) 
 
coup2(1) = cplChiFucSuL(iIN,gt3,i1) 
coup2(2) = cplChiFucSuR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcFuChiSuL(gt2,gt1,i1) 
coup2(6) = cplcFuChiSuR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFuFuAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: conj[Su],Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su],Ah'
        End Do 
      End Do 



!-------------- 
!  conj[Su], hh 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFu(gt3) 
mass(3) = -MFu(gt2) 
 
coup2(1) = cplChiFucSuL(iIN,gt3,i1) 
coup2(2) = cplChiFucSuR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFuChiSuL(gt2,gt1,i1) 
coup2(6) = cplcFuChiSuR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFuFuhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuhhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: conj[Su],hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su],hh'
        End Do 
      End Do 



!-------------- 
!  Su, Su 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt2) 
mass(3) = -MFu(gt3) 
mass(4) = MChi(gt1) 
 
coup2(1) = cplcFuChiSuL(gt2,iIN,i1) 
coup2(2) = cplcFuChiSuR(gt2,iIN,i1) 
coup2(3) = Conjg(cplcFuChiSuL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFuChiSuR(gt2,iIN,i2))  
coup2(5) = cplChiFucSuL(gt1,gt3,i1) 
coup2(6) = cplChiFucSuR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChiFucSuL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChiFucSuR(gt1,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Su,Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Su,Su'
        End Do 
      End Do 



!-------------- 
!  Su, Ah 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFu(gt2) 
mass(3) = -MFu(gt2) 
 
coup2(1) = cplcFuChiSuL(gt2,iIN,i1) 
coup2(2) = cplcFuChiSuR(gt2,iIN,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplChiFucSuL(gt1,gt3,i1) 
coup2(6) = cplChiFucSuR(gt1,gt3,i1) 
coup2(7) = Conjg(cplcFuFuAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuAhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Su,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Su,Ah'
        End Do 
      End Do 



!-------------- 
!  Su, hh 
!-------------- 
Do i1=1,6
  Do i2=1,8
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(4) = MFu(gt2) 
mass(3) = -MFu(gt2) 
 
coup2(1) = cplcFuChiSuL(gt2,iIN,i1) 
coup2(2) = cplcFuChiSuR(gt2,iIN,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplChiFucSuL(gt1,gt3,i1) 
coup2(6) = cplChiFucSuR(gt1,gt3,i1) 
coup2(7) = Conjg(cplcFuFuhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuhhR(gt2,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Su,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Su,hh'
        End Do 
      End Do 



!-------------- 
!  Ah, Ah 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = MAh(i2) 
Boson4(4) = gTAh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChiAhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChiAhR(iIN,gt1,i2))  
coup2(5) = cplcFuFuAhL(gt2,gt3,i1) 
coup2(6) = cplcFuFuAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFuFuAhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuAhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Ah,Ah" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,Ah'
        End Do 
      End Do 



!-------------- 
!  Ah, hh 
!-------------- 
Do i1=1,8
  Do i2=1,8
Boson4(1) = MAh(i1) 
Boson4(2) = gTAh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChiAhL(iIN,gt1,i1) 
coup2(2) = cplChiChiAhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFuFuAhL(gt2,gt3,i1) 
coup2(6) = cplcFuFuAhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFuFuhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuhhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: Ah,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Ah,hh'
        End Do 
      End Do 



!-------------- 
!  hh, hh 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = Mhh(i1) 
Boson4(2) = gThh(i1) 
Boson4(3) = Mhh(i2) 
Boson4(4) = gThh(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MChi(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChihhL(iIN,gt1,i1) 
coup2(2) = cplChiChihhR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChihhL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChihhR(iIN,gt1,i2))  
coup2(5) = cplcFuFuhhL(gt2,gt3,i1) 
coup2(6) = cplcFuFuhhR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFuFuhhL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuFuhhR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Chi cFu Fu Propagator: hh,hh" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='hh,hh'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,gt3,:)= 0._dp  
End If 
       End Do 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1, iIN-1
      Do gt2=1,3
        Do gt3=1,3
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:841))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToChicFuFu 
 
 
Subroutine ChiToChacFdFu(iIN,MCha,MFd,MFu,MVWm,MHpm,MSu,MSd,MChi,cplcFdChaSuL,        & 
& cplcFdChaSuR,cplcFdChiSdL,cplcFdChiSdR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdFuVWmL,         & 
& cplcFdFuVWmR,cplChaFucSdL,cplChaFucSdR,cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,   & 
& cplChiChacVWmR,cplChiFucSuL,cplChiFucSuR,IntegralSs,IntegralSSss,IntegralSSst,         & 
& IntegralVs,IntegralVSss,IntegralVSst,IntegralVVss,NSs,NSSss,NSSst,NVs,NVSss,           & 
& NVSst,NVVss,gTVWm,gTHpm,gTSu,gTSd,deltaM,epsI,check,g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MCha(5),MFd(3),MFu(3),MVWm,MHpm(8),MSu(6),MSd(6),MChi(10)

Complex(dp),Intent(in) :: cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),    & 
& cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),           & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8), & 
& cplChiChacVWmL(10,5),cplChiChacVWmR(10,5),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6)

Real(dp),Intent(inout) :: IntegralSs(500000,10),IntegralVs(25000,8),IntegralVVss(500000,12)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16),IntegralVSss(500000,12),              & 
& IntegralVSst(500000,16)

Real(dp),Intent(inout) :: gTVWm,gTHpm(8),gTSu(6),gTSd(6)

Integer, Intent(inout) :: NSs,NSSss,NSSst,NVs,NVSss,NVSst,NVVss
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 441 
Allocate( gSum(5,3,3, Isum) ) 
Allocate( Contribution(5,3,3, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1,5
      Do gt2=1,3
        Do gt3=1,3
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MFu(gt3))+Abs(MFd(gt2))+Abs(MCha(gt1)))) Then 
!-------------- 
!  conj[VWm] 
!-------------- 
Isum = Isum + 1 
Boson2(1) = MVWm 
Boson2(2) =gTVWm 
 
Boson4(1) = MVWm 
Boson4(2) =gTVWm 
Boson4(3) = MVWm 
Boson4(4) =gTVWm 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MCha(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFu(gt3) 
 
coup(2) = Conjg(cplChiChacVWmL(iIN,gt1)) 
coup(1) = Conjg(cplChiChacVWmR(iIN,gt1)) 
coup(4) = Conjg(cplcFdFuVWmL(gt2,gt3)) 
coup(3) = Conjg(cplcFdFuVWmR(gt2,gt3))
Call IntegrateGaugeSS(Boson2,mass,coup,deltaM,epsI,IntegralVs,NVs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[VWm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm]'



!-------------- 
!  conj[Hpm] 
!-------------- 
Do i1=1,8
Isum = Isum + 1 
Boson2(1) = MHpm(i1) 
Boson2(2) =gTHpm(i1) 
 
Boson4(1) = MHpm(i1) 
Boson4(2) =gTHpm(i1) 
Boson4(3) = MHpm(i1) 
Boson4(4) =gTHpm(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MCha(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFu(gt3) 
 
coup(2) = Conjg(cplChiChacHpmL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiChacHpmR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFdFuHpmL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFdFuHpmR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm]'
      End Do 



!-------------- 
!  conj[Su] 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSu(i1) 
Boson2(2) =gTSu(i1) 
 
Boson4(1) = MSu(i1) 
Boson4(2) =gTSu(i1) 
Boson4(3) = MSu(i1) 
Boson4(4) =gTSu(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFu(gt3) 
mass(3) = -MFd(gt2) 
mass(4) = MCha(gt1) 
 
coup(2) = Conjg(cplChiFucSuL(iIN,gt3,i1)) 
coup(1) = Conjg(cplChiFucSuR(iIN,gt3,i1)) 
coup(4) = Conjg(cplcFdChaSuL(gt2,gt1,i1)) 
coup(3) = Conjg(cplcFdChaSuR(gt2,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su]'
      End Do 



!-------------- 
!  Sd 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSd(i1) 
Boson2(2) =gTSd(i1) 
 
Boson4(1) = MSd(i1) 
Boson4(2) =gTSd(i1) 
Boson4(3) = MSd(i1) 
Boson4(4) =gTSd(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFd(gt2) 
mass(3) = -MFu(gt3) 
mass(4) = MCha(gt1) 
 
coup(2) = Conjg(cplcFdChiSdL(gt2,iIN,i1)) 
coup(1) = Conjg(cplcFdChiSdR(gt2,iIN,i1)) 
coup(4) = Conjg(cplChaFucSdL(gt1,gt3,i1)) 
coup(3) = Conjg(cplChaFucSdR(gt1,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 3*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp
Else 
gSum(gt1,gt2,gt3,Isum)=resD
End If 
Contribution(gt1,gt2,gt3,Isum)='Sd'
      End Do 



!-------------- 
!  conj[VWm], conj[Hpm] 
!-------------- 
  Do i2=1,8
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChacVWmL(iIN,gt1) 
coup2(2) = cplChiChacVWmR(iIN,gt1) 
coup2(3) = Conjg(cplChiChacHpmL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChacHpmR(iIN,gt1,i2))  
coup2(5) = cplcFdFuVWmL(gt2,gt3) 
coup2(6) = cplcFdFuVWmR(gt2,gt3) 
coup2(7) = Conjg(cplcFdFuHpmL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFuHpmR(gt2,gt3,i2)) 
Call IntegrateGaugeSscalarS(Boson4, mass, coup2, deltaM, epsI,IntegralVSss,NVSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[VWm],conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],conj[Hpm]'
      End Do 



!-------------- 
!  conj[VWm], conj[Su] 
!-------------- 
  Do i2=1,6
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt3) 
mass(4) = MCha(gt1) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplChiChacVWmL(iIN,gt1) 
coup2(2) = cplChiChacVWmR(iIN,gt1) 
coup2(3) = Conjg(cplChiFucSuL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFucSuR(iIN,gt3,i2))  
coup2(5) = cplcFdFuVWmL(gt2,gt3) 
coup2(6) = cplcFdFuVWmR(gt2,gt3) 
coup2(7) = Conjg(cplcFdChaSuL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFdChaSuR(gt2,gt1,i2)) 
coupT = coup2(5) 
coup2(5) = coup2(6) 
coup2(6) = coupT 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[VWm],conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],conj[Su]'
      End Do 



!-------------- 
!  conj[VWm], Sd 
!-------------- 
  Do i2=1,6
Boson4(1) = MVWm 
Boson4(2) = gTVWm 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(4) = MCha(gt1) 
mass(3) = -MFu(gt3) 
 
coup2(1) = cplChiChacVWmL(iIN,gt1) 
coup2(2) = cplChiChacVWmR(iIN,gt1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplcFdFuVWmL(gt2,gt3) 
coup2(6) = cplcFdFuVWmR(gt2,gt3) 
coup2(7) = Conjg(cplChaFucSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChaFucSdR(gt1,gt3,i2)) 
Call IntegrateGaugeSscalarT(Boson4, mass, coup2, deltaM, epsI,IntegralVSst,NVSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[VWm],Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[VWm],Sd'
      End Do 



!-------------- 
!  conj[Hpm], conj[Hpm] 
!-------------- 
Do i1=1,7
  Do i2=i1+1,8
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MHpm(i2) 
Boson4(4) = gTHpm(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MCha(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MFu(gt3) 
 
coup2(1) = cplChiChacHpmL(iIN,gt1,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiChacHpmL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiChacHpmR(iIN,gt1,i2))  
coup2(5) = cplcFdFuHpmL(gt2,gt3,i1) 
coup2(6) = cplcFdFuHpmR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFdFuHpmL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdFuHpmR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Hpm],conj[Hpm]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],conj[Hpm]'
        End Do 
      End Do 



!-------------- 
!  conj[Hpm], conj[Su] 
!-------------- 
Do i1=1,8
  Do i2=1,6
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt3) 
mass(4) = MCha(gt1) 
mass(3) = -MFd(gt2) 
 
coup2(1) = cplChiChacHpmL(iIN,gt1,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiFucSuL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFucSuR(iIN,gt3,i2))  
coup2(5) = cplcFdFuHpmL(gt2,gt3,i1) 
coup2(6) = cplcFdFuHpmR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFdChaSuL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFdChaSuR(gt2,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Hpm],conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],conj[Su]'
        End Do 
      End Do 



!-------------- 
!  conj[Hpm], Sd 
!-------------- 
Do i1=1,8
  Do i2=1,6
Boson4(1) = MHpm(i1) 
Boson4(2) = gTHpm(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(4) = MCha(gt1) 
mass(3) = -MFu(gt3) 
 
coup2(1) = cplChiChacHpmL(iIN,gt1,i1) 
coup2(2) = cplChiChacHpmR(iIN,gt1,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplcFdFuHpmL(gt2,gt3,i1) 
coup2(6) = cplcFdFuHpmR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChaFucSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChaFucSdR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Hpm],Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Hpm],Sd'
        End Do 
      End Do 



!-------------- 
!  conj[Su], conj[Su] 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt3) 
mass(3) = -MFd(gt2) 
mass(4) = MCha(gt1) 
 
coup2(1) = cplChiFucSuL(iIN,gt3,i1) 
coup2(2) = cplChiFucSuR(iIN,gt3,i1) 
coup2(3) = Conjg(cplChiFucSuL(iIN,gt3,i2)) 
coup2(4) = Conjg(cplChiFucSuR(iIN,gt3,i2))  
coup2(5) = cplcFdChaSuL(gt2,gt1,i1) 
coup2(6) = cplcFdChaSuR(gt2,gt1,i1) 
coup2(7) = Conjg(cplcFdChaSuL(gt2,gt1,i2)) 
coup2(8) = Conjg(cplcFdChaSuR(gt2,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Su],conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su],conj[Su]'
        End Do 
      End Do 



!-------------- 
!  conj[Su], Sd 
!-------------- 
Do i1=1,6
  Do i2=1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(4) = MFu(gt3) 
mass(3) = -MCha(gt1) 
 
coup2(1) = cplChiFucSuL(iIN,gt3,i1) 
coup2(2) = cplChiFucSuR(iIN,gt3,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplcFdChaSuL(gt2,gt1,i1) 
coup2(6) = cplcFdChaSuR(gt2,gt1,i1) 
coup2(7) = Conjg(cplChaFucSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChaFucSdR(gt1,gt3,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: conj[Su],Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Su],Sd'
        End Do 
      End Do 



!-------------- 
!  Sd, Sd 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(3) = -MFu(gt3) 
mass(4) = MCha(gt1) 
 
coup2(1) = cplcFdChiSdL(gt2,iIN,i1) 
coup2(2) = cplcFdChiSdR(gt2,iIN,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplChaFucSdL(gt1,gt3,i1) 
coup2(6) = cplChaFucSdR(gt1,gt3,i1) 
coup2(7) = Conjg(cplChaFucSdL(gt1,gt3,i2)) 
coup2(8) = Conjg(cplChaFucSdR(gt1,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 3*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Cha cFd Fu Propagator: Sd,Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='Sd,Sd'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,gt3,:)= 0._dp  
End If 
       End Do 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1,5
      Do gt2=1,3
        Do gt3=1,3
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:441))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToChacFdFu 
 
 
Subroutine ChiToFdcFdGlu(iIN,MFd,MGlu,MSd,MChi,cplcFdChiSdL,cplcFdChiSdR,             & 
& cplcFdGluSdL,cplcFdGluSdR,cplChiFdcSdL,cplChiFdcSdR,cplGluFdcSdL,cplGluFdcSdR,         & 
& IntegralSs,IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSd,deltaM,epsI,check,           & 
& g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MFd(3),MGlu,MSd(6),MChi(10)

Complex(dp),Intent(in) :: cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),        & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6)

Real(dp),Intent(inout) :: IntegralSs(500000,10)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16)

Real(dp),Intent(inout) :: gTSd(6)

Integer, Intent(inout) :: NSs,NSSss,NSSst
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 144 
Allocate( gSum(3,3,1, Isum) ) 
Allocate( Contribution(3,3,1, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1,3
      Do gt2=1,3
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MGlu)+Abs(MFd(gt2))+Abs(MFd(gt1)))) Then 
!-------------- 
!  conj[Sd] 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSd(i1) 
Boson2(2) =gTSd(i1) 
 
Boson4(1) = MSd(i1) 
Boson4(2) =gTSd(i1) 
Boson4(3) = MSd(i1) 
Boson4(4) =gTSd(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFd(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MGlu 
 
coup(2) = Conjg(cplChiFdcSdL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiFdcSdR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFdGluSdL(gt2,i1)) 
coup(3) = Conjg(cplcFdGluSdR(gt2,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 4*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fd cFd Glu Propagator: conj[Sd]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp
Else 
gSum(gt1,gt2,1,Isum)=resD
End If 
Contribution(gt1,gt2,1,Isum)='conj[Sd]'
      End Do 



!-------------- 
!  Sd 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSd(i1) 
Boson2(2) =gTSd(i1) 
 
Boson4(1) = MSd(i1) 
Boson4(2) =gTSd(i1) 
Boson4(3) = MSd(i1) 
Boson4(4) =gTSd(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFd(gt2) 
mass(3) = -MGlu 
mass(4) = MFd(gt1) 
 
coup(2) = Conjg(cplcFdChiSdL(gt2,iIN,i1)) 
coup(1) = Conjg(cplcFdChiSdR(gt2,iIN,i1)) 
coup(4) = Conjg(cplGluFdcSdL(gt1,i1)) 
coup(3) = Conjg(cplGluFdcSdR(gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 4*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fd cFd Glu Propagator: Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp
Else 
gSum(gt1,gt2,1,Isum)=resD
End If 
Contribution(gt1,gt2,1,Isum)='Sd'
      End Do 



!-------------- 
!  conj[Sd], conj[Sd] 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt1) 
mass(3) = -MFd(gt2) 
mass(4) = MGlu 
 
coup2(1) = cplChiFdcSdL(iIN,gt1,i1) 
coup2(2) = cplChiFdcSdR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiFdcSdL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiFdcSdR(iIN,gt1,i2))  
coup2(5) = cplcFdGluSdL(gt2,i1) 
coup2(6) = cplcFdGluSdR(gt2,i1) 
coup2(7) = Conjg(cplcFdGluSdL(gt2,i2)) 
coup2(8) = Conjg(cplcFdGluSdR(gt2,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 4*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fd cFd Glu Propagator: conj[Sd],conj[Sd]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp  
Else 
gSum(gt1,gt2,1,Isum)= resS  
End If 
Contribution(gt1,gt2,1,Isum)='conj[Sd],conj[Sd]'
        End Do 
      End Do 



!-------------- 
!  conj[Sd], Sd 
!-------------- 
Do i1=1,6
  Do i2=1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(4) = MFd(gt1) 
mass(3) = -MGlu 
 
coup2(1) = cplChiFdcSdL(iIN,gt1,i1) 
coup2(2) = cplChiFdcSdR(iIN,gt1,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplcFdGluSdL(gt2,i1) 
coup2(6) = cplcFdGluSdR(gt2,i1) 
coup2(7) = Conjg(cplGluFdcSdL(gt1,i2)) 
coup2(8) = Conjg(cplGluFdcSdR(gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 4*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fd cFd Glu Propagator: conj[Sd],Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp  
Else 
gSum(gt1,gt2,1,Isum)= resS  
End If 
Contribution(gt1,gt2,1,Isum)='conj[Sd],Sd'
        End Do 
      End Do 



!-------------- 
!  Sd, Sd 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSd(i2) 
Boson4(4) = gTSd(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFd(gt2) 
mass(3) = -MGlu 
mass(4) = MFd(gt1) 
 
coup2(1) = cplcFdChiSdL(gt2,iIN,i1) 
coup2(2) = cplcFdChiSdR(gt2,iIN,i1) 
coup2(3) = Conjg(cplcFdChiSdL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFdChiSdR(gt2,iIN,i2))  
coup2(5) = cplGluFdcSdL(gt1,i1) 
coup2(6) = cplGluFdcSdR(gt1,i1) 
coup2(7) = Conjg(cplGluFdcSdL(gt1,i2)) 
coup2(8) = Conjg(cplGluFdcSdR(gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 4*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fd cFd Glu Propagator: Sd,Sd" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp  
Else 
gSum(gt1,gt2,1,Isum)= resS  
End If 
Contribution(gt1,gt2,1,Isum)='Sd,Sd'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,1,:)= 0._dp  
End If 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1,3
      Do gt2=1,3
g(gt1,gt2,1)=Sum(gSum(gt1,gt2,1,1:144))
If (g(gt1,gt2,1).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,1)=0._dp
End If
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToFdcFdGlu 
 
 
Subroutine ChiToFucFuGlu(iIN,MFu,MGlu,MSu,MChi,cplcFuChiSuL,cplcFuChiSuR,             & 
& cplcFuGluSuL,cplcFuGluSuR,cplChiFucSuL,cplChiFucSuR,cplGluFucSuL,cplGluFucSuR,         & 
& IntegralSs,IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSu,deltaM,epsI,check,           & 
& g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MFu(3),MGlu,MSu(6),MChi(10)

Complex(dp),Intent(in) :: cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),        & 
& cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6)

Real(dp),Intent(inout) :: IntegralSs(500000,10)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16)

Real(dp),Intent(inout) :: gTSu(6)

Integer, Intent(inout) :: NSs,NSSss,NSSst
Real(dp),Intent(inout)::g(:,:,:) 
Logical, Intent(in) :: check 
Integer, Intent(in) :: iIN 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Optional :: WriteContributions 
Integer :: i1,i2,gt1,gt2,gt3, Isum 
Real(dp) :: resR,  res1, res2, resD, m_in 
Complex(dp) :: resC, resS 
Real(dp), Allocatable :: gSum(:,:,:,:) 
Character(len=20), Allocatable :: Contribution(:,:,:,:) 
Real(dp) :: Boson2(2), mass(4),  Boson4(4) 
Complex(dp) :: coup(4), coup2(8),coupT 
mass(1) = MChi(iIN) 
 
Isum = 144 
Allocate( gSum(3,3,1, Isum) ) 
Allocate( Contribution(3,3,1, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1,3
      Do gt2=1,3
Isum = 0 
 
If(Abs(MChi(iIN)).gt.(Abs(MGlu)+Abs(MFu(gt2))+Abs(MFu(gt1)))) Then 
!-------------- 
!  conj[Su] 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSu(i1) 
Boson2(2) =gTSu(i1) 
 
Boson4(1) = MSu(i1) 
Boson4(2) =gTSu(i1) 
Boson4(3) = MSu(i1) 
Boson4(4) =gTSu(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFu(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MGlu 
 
coup(2) = Conjg(cplChiFucSuL(iIN,gt1,i1)) 
coup(1) = Conjg(cplChiFucSuR(iIN,gt1,i1)) 
coup(4) = Conjg(cplcFuGluSuL(gt2,i1)) 
coup(3) = Conjg(cplcFuGluSuR(gt2,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 4*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fu cFu Glu Propagator: conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp
Else 
gSum(gt1,gt2,1,Isum)=resD
End If 
Contribution(gt1,gt2,1,Isum)='conj[Su]'
      End Do 



!-------------- 
!  Su 
!-------------- 
Do i1=1,6
Isum = Isum + 1 
Boson2(1) = MSu(i1) 
Boson2(2) =gTSu(i1) 
 
Boson4(1) = MSu(i1) 
Boson4(2) =gTSu(i1) 
Boson4(3) = MSu(i1) 
Boson4(4) =gTSu(i1) 
 
resS=0._dp 
resD=0._dp 
 
mass(2) = MFu(gt2) 
mass(3) = -MGlu 
mass(4) = MFu(gt1) 
 
coup(2) = Conjg(cplcFuChiSuL(gt2,iIN,i1)) 
coup(1) = Conjg(cplcFuChiSuR(gt2,iIN,i1)) 
coup(4) = Conjg(cplGluFucSuL(gt1,i1)) 
coup(3) = Conjg(cplGluFucSuR(gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 4*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fu cFu Glu Propagator: Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp
Else 
gSum(gt1,gt2,1,Isum)=resD
End If 
Contribution(gt1,gt2,1,Isum)='Su'
      End Do 



!-------------- 
!  conj[Su], conj[Su] 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt1) 
mass(3) = -MFu(gt2) 
mass(4) = MGlu 
 
coup2(1) = cplChiFucSuL(iIN,gt1,i1) 
coup2(2) = cplChiFucSuR(iIN,gt1,i1) 
coup2(3) = Conjg(cplChiFucSuL(iIN,gt1,i2)) 
coup2(4) = Conjg(cplChiFucSuR(iIN,gt1,i2))  
coup2(5) = cplcFuGluSuL(gt2,i1) 
coup2(6) = cplcFuGluSuR(gt2,i1) 
coup2(7) = Conjg(cplcFuGluSuL(gt2,i2)) 
coup2(8) = Conjg(cplcFuGluSuR(gt2,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 4*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fu cFu Glu Propagator: conj[Su],conj[Su]" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp  
Else 
gSum(gt1,gt2,1,Isum)= resS  
End If 
Contribution(gt1,gt2,1,Isum)='conj[Su],conj[Su]'
        End Do 
      End Do 



!-------------- 
!  conj[Su], Su 
!-------------- 
Do i1=1,6
  Do i2=1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt2) 
mass(4) = MFu(gt1) 
mass(3) = -MGlu 
 
coup2(1) = cplChiFucSuL(iIN,gt1,i1) 
coup2(2) = cplChiFucSuR(iIN,gt1,i1) 
coup2(3) = Conjg(cplcFuChiSuL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFuChiSuR(gt2,iIN,i2))  
coup2(5) = cplcFuGluSuL(gt2,i1) 
coup2(6) = cplcFuGluSuR(gt2,i1) 
coup2(7) = Conjg(cplGluFucSuL(gt1,i2)) 
coup2(8) = Conjg(cplGluFucSuR(gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 4*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fu cFu Glu Propagator: conj[Su],Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp  
Else 
gSum(gt1,gt2,1,Isum)= resS  
End If 
Contribution(gt1,gt2,1,Isum)='conj[Su],Su'
        End Do 
      End Do 



!-------------- 
!  Su, Su 
!-------------- 
Do i1=1,5
  Do i2=i1+1,6
Boson4(1) = MSu(i1) 
Boson4(2) = gTSu(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt2) 
mass(3) = -MGlu 
mass(4) = MFu(gt1) 
 
coup2(1) = cplcFuChiSuL(gt2,iIN,i1) 
coup2(2) = cplcFuChiSuR(gt2,iIN,i1) 
coup2(3) = Conjg(cplcFuChiSuL(gt2,iIN,i2)) 
coup2(4) = Conjg(cplcFuChiSuR(gt2,iIN,i2))  
coup2(5) = cplGluFucSuL(gt1,i1) 
coup2(6) = cplGluFucSuR(gt1,i1) 
coup2(7) = Conjg(cplGluFucSuL(gt1,i2)) 
coup2(8) = Conjg(cplGluFucSuR(gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 4*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Chi->Fu cFu Glu Propagator: Su,Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,1,Isum)= 0._dp  
Else 
gSum(gt1,gt2,1,Isum)= resS  
End If 
Contribution(gt1,gt2,1,Isum)='Su,Su'
        End Do 
      End Do 



Else 
gSum(gt1,gt2,1,:)= 0._dp  
End If 
     End Do 
   End Do 
!---------- 
!Summing 
!---------- 
g=0._dp 
    Do gt1=1,3
      Do gt2=1,3
g(gt1,gt2,1)=Sum(gSum(gt1,gt2,1,1:144))
If (g(gt1,gt2,1).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,1)=0._dp
End If
     End Do 
   End Do 
  g = oo512pi3 / Abs(MChi(iIN))**3*g
End Subroutine ChiToFucFuGlu 
 
 
End Module Chi3Decays_munuSSM3G 
 
