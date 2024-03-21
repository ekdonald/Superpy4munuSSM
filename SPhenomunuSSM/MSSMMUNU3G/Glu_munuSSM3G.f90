! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:24 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module Glu3Decays_munuSSM3G 
 
Use Control 
Use CouplingsForDecays_munuSSM3G 
Use ThreeBodyPhaseSpace 
 
Contains 
 
Subroutine GluThreeBodyDecay(n_in,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,            & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,              & 
& MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,Yd,Ye,             & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,vd,vu,vL,vR,gTSd,gTSu,gGluFdcFdChi,gGluFdcFucCha,gGluFucFuChi,epsI,           & 
& deltaM,CheckRealStates,gT,gPartial,BR)

Implicit None 
 
Real(dp),Intent(in) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Complex(dp) :: cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),  & 
& cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),         & 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),         & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6)

Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(inout) :: gGluFdcFdChi(1,3,3,10),gGluFdcFucCha(1,3,3,5),gGluFucFuChi(1,3,3,10)

Real(dp),Intent(in) :: gTSd(6),gTSu(6)

Real(dp) :: gGluFdcFdChii(3,3,10),gGluFdcFucChai(3,3,5),gGluFucFuChii(3,3,10)

Real(dp) :: gTSdtemp(6),gTSutemp(6)
Integer :: NVs,NVst,NSs,NVVst,NVVss,NVSss,NVSst,NSSss,NSSst
Complex(dp), Allocatable :: IntegralVVst(:,:),IntegralVSss(:,:),IntegralVSst(:,:),IntegralSSss(:,:)               & 
& ,IntegralSSst(:,:)
Real(dp), Allocatable :: IntegralVs(:,:),IntegralVst(:,:),IntegralSs(:,:),IntegralVVss(:,:)
Real(dp), Intent(inout), Optional :: BR(:,:), gPartial(:,:) 
Real(dp), Intent(inout) :: gT 
Integer, Intent(in) :: n_in 
Real(dp), Intent(in) :: epsI, deltaM 
Logical, Intent(in) ::  CheckRealStates 
Integer :: i_start, i_end, i_run, n_out, n_length, gt1, gt2, gt3, i1 
Logical :: check 
Iname = Iname +1 
NameOfUnit(Iname) = 'GluThreeBodyDecay' 
 
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
gTSdtemp = 0._dp 
gTSutemp = 0._dp 
Else 
gTSdtemp = gTSd 
gTSutemp = gTSu 
End If 
 
check=CheckRealStates 

 
If (n_in.Lt.0) Then 
 i_start = 1 
 i_end = 1 
 Else If ( (n_in.Ge.1).And.(n_in.Le. 1) ) Then 
 i_start = n_in 
 i_end = n_in 
Else 
 If (ErrorLevel.Ge.-1) Then 
   Write (ErrCan, *) 'Problem in subroutine'//NameOfUnit(Iname) 
   Write (ErrCan, *) 'Value of n_in out of range, (n_in,1) = ',n_in,1 
 End If 
 
 If (ErrorLevel.Gt.0) then
   Call TerminateProgram 
   return
 endif 
 If (Present(BR)) BR = 0._dp 
 Iname = Iname - 1 
 Return 
End If 
 
i_run = 1 
 
Call CouplingsFor_Glu_decays_3B(MGlu,i_run,MAh,MAh2,MCha,MCha2,MChi,MChi2,            & 
& MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,               & 
& g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,               & 
& me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,cplcChacFuSdL,cplcChacFuSdR,cplcChaFdcSuL,          & 
& cplcChaFdcSuR,cplcFdChiSdL,cplcFdChiSdR,cplcFdGluSdL,cplcFdGluSdR,cplcFuChiSuL,        & 
& cplcFuChiSuR,cplcFuGluSuL,cplcFuGluSuR,cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,         & 
& cplChiFucSuR,cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,deltaM)

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

 
gGluFdcFdChii = 0._dp 
Call GluToFdcFdChi(i_run,MFd,MChi,MSd,MGlu,cplcFdChiSdL,cplcFdChiSdR,cplcFdGluSdL,    & 
& cplcFdGluSdR,cplChiFdcSdL,cplChiFdcSdR,cplGluFdcSdL,cplGluFdcSdR,IntegralSs,           & 
& IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSdtemp,deltaM,epsI,check,gGluFdcFdChii)

gGluFdcFdChi(i_run,:,:,:) = gGluFdcFdChii 
gT = gT + Sum(gGluFdcFdChii) 
 
gGluFdcFucChai = 0._dp 
Call GluToFdcFucCha(i_run,MFd,MFu,MCha,MSd,MSu,MGlu,cplcChacFuSdL,cplcChacFuSdR,      & 
& cplcChaFdcSuL,cplcChaFdcSuR,cplcFuGluSuL,cplcFuGluSuR,cplGluFdcSdL,cplGluFdcSdR,       & 
& IntegralSs,IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSdtemp,gTSutemp,deltaM,         & 
& epsI,check,gGluFdcFucChai)

gGluFdcFucCha(i_run,:,:,:) = gGluFdcFucChai 
gT = gT + 2._dp*Sum(gGluFdcFucChai) 
 
gGluFucFuChii = 0._dp 
Call GluToFucFuChi(i_run,MFu,MChi,MSu,MGlu,cplcFuChiSuL,cplcFuChiSuR,cplcFuGluSuL,    & 
& cplcFuGluSuR,cplChiFucSuL,cplChiFucSuR,cplGluFucSuL,cplGluFucSuR,IntegralSs,           & 
& IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSutemp,deltaM,epsI,check,gGluFucFuChii)

gGluFucFuChi(i_run,:,:,:) = gGluFucFuChii 
gT = gT + Sum(gGluFucFuChii) 
 
If (Present(gPartial)) Then
Do i1 = i_start, i_end 
 
n_length=1
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,10
gPartial(i1,n_length)= gGluFdcFdChi(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,5
gPartial(i1,n_length)= 2._dp*gGluFdcFucCha(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,10
gPartial(i1,n_length)= gGluFucFuChi(i1,gt1,gt2,gt3)
n_length=n_length+1
     End Do 
  End Do 
End Do 
If (Present(BR).And.(gT.Gt.0._dp)) Then 
BR(i1,:)=gPartial(i1,:)/gT
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
 
End Subroutine GluThreeBodyDecay
 
 
Subroutine GluToFdcFdChi(iIN,MFd,MChi,MSd,MGlu,cplcFdChiSdL,cplcFdChiSdR,             & 
& cplcFdGluSdL,cplcFdGluSdR,cplChiFdcSdL,cplChiFdcSdR,cplGluFdcSdL,cplGluFdcSdR,         & 
& IntegralSs,IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSd,deltaM,epsI,check,           & 
& g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MFd(3),MChi(10),MSd(6),MGlu

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
mass(1) = MGlu 
 
Isum = 144 
Allocate( gSum(3,3,10, Isum) ) 
Allocate( Contribution(3,3,10, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1,3
      Do gt2=1,3
        Do gt3=1,10
Isum = 0 
 
If(Abs(MGlu).gt.(Abs(MChi(gt3))+Abs(MFd(gt2))+Abs(MFd(gt1)))) Then 
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
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplGluFdcSdL(gt1,i1)) 
coup(1) = Conjg(cplGluFdcSdR(gt1,i1)) 
coup(4) = Conjg(cplcFdChiSdL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFdChiSdR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1._dp/2._dp*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFd Chi Propagator: conj[Sd]" 
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
mass(3) = -MFd(gt1) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplcFdGluSdL(gt2,i1)) 
coup(1) = Conjg(cplcFdGluSdR(gt2,i1)) 
coup(4) = Conjg(cplChiFdcSdL(gt3,gt1,i1)) 
coup(3) = Conjg(cplChiFdcSdR(gt3,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1._dp/2._dp*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFd Chi Propagator: Sd" 
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
mass(4) = MChi(gt3) 
 
coup2(1) = cplGluFdcSdL(gt1,i1) 
coup2(2) = cplGluFdcSdR(gt1,i1) 
coup2(3) = Conjg(cplGluFdcSdL(gt1,i2)) 
coup2(4) = Conjg(cplGluFdcSdR(gt1,i2))  
coup2(5) = cplcFdChiSdL(gt2,gt3,i1) 
coup2(6) = cplcFdChiSdR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFdChiSdL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFdChiSdR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFd Chi Propagator: conj[Sd],conj[Sd]" 
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
mass(4) = MFd(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplGluFdcSdL(gt1,i1) 
coup2(2) = cplGluFdcSdR(gt1,i1) 
coup2(3) = Conjg(cplcFdGluSdL(gt2,i2)) 
coup2(4) = Conjg(cplcFdGluSdR(gt2,i2))  
coup2(5) = cplcFdChiSdL(gt2,gt3,i1) 
coup2(6) = cplcFdChiSdR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiFdcSdL(gt3,gt1,i2)) 
coup2(8) = Conjg(cplChiFdcSdR(gt3,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFd Chi Propagator: conj[Sd],Sd" 
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
mass(3) = -MFd(gt1) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplcFdGluSdL(gt2,i1) 
coup2(2) = cplcFdGluSdR(gt2,i1) 
coup2(3) = Conjg(cplcFdGluSdL(gt2,i2)) 
coup2(4) = Conjg(cplcFdGluSdR(gt2,i2))  
coup2(5) = cplChiFdcSdL(gt3,gt1,i1) 
coup2(6) = cplChiFdcSdR(gt3,gt1,i1) 
coup2(7) = Conjg(cplChiFdcSdL(gt3,gt1,i2)) 
coup2(8) = Conjg(cplChiFdcSdR(gt3,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFd Chi Propagator: Sd,Sd" 
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
    Do gt1=1,3
      Do gt2=1,3
        Do gt3=1,10
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:144))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MGlu)**3*g
End Subroutine GluToFdcFdChi 
 
 
Subroutine GluToFdcFucCha(iIN,MFd,MFu,MCha,MSd,MSu,MGlu,cplcChacFuSdL,cplcChacFuSdR,  & 
& cplcChaFdcSuL,cplcChaFdcSuR,cplcFuGluSuL,cplcFuGluSuR,cplGluFdcSdL,cplGluFdcSdR,       & 
& IntegralSs,IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSd,gTSu,deltaM,epsI,            & 
& check,g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MFd(3),MFu(3),MCha(5),MSd(6),MSu(6),MGlu

Complex(dp),Intent(in) :: cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),  & 
& cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6)

Real(dp),Intent(inout) :: IntegralSs(500000,10)

Complex(dp),Intent(inout) :: IntegralSSss(500000,12),IntegralSSst(500000,16)

Real(dp),Intent(inout) :: gTSd(6),gTSu(6)

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
mass(1) = MGlu 
 
Isum = 144 
Allocate( gSum(3,3,5, Isum) ) 
Allocate( Contribution(3,3,5, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1,3
      Do gt2=1,3
        Do gt3=1,5
Isum = 0 
 
If(Abs(MGlu).gt.(Abs(MCha(gt3))+Abs(MFu(gt2))+Abs(MFd(gt1)))) Then 
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
mass(3) = -MFu(gt2) 
mass(4) = MCha(gt3) 
 
coup(2) = Conjg(cplGluFdcSdL(gt1,i1)) 
coup(1) = Conjg(cplGluFdcSdR(gt1,i1)) 
coup(4) = Conjg(cplcChacFuSdL(gt3,gt2,i1)) 
coup(3) = Conjg(cplcChacFuSdR(gt3,gt2,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1._dp/2._dp*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFu cCha Propagator: conj[Sd]" 
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
mass(3) = -MCha(gt3) 
mass(4) = MFd(gt1) 
 
coup(2) = Conjg(cplcFuGluSuL(gt2,i1)) 
coup(1) = Conjg(cplcFuGluSuR(gt2,i1)) 
coup(4) = Conjg(cplcChaFdcSuL(gt3,gt1,i1)) 
coup(3) = Conjg(cplcChaFdcSuR(gt3,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1._dp/2._dp*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFu cCha Propagator: Su" 
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
mass(3) = -MFu(gt2) 
mass(4) = MCha(gt3) 
 
coup2(1) = cplGluFdcSdL(gt1,i1) 
coup2(2) = cplGluFdcSdR(gt1,i1) 
coup2(3) = Conjg(cplGluFdcSdL(gt1,i2)) 
coup2(4) = Conjg(cplGluFdcSdR(gt1,i2))  
coup2(5) = cplcChacFuSdL(gt3,gt2,i1) 
coup2(6) = cplcChacFuSdR(gt3,gt2,i1) 
coup2(7) = Conjg(cplcChacFuSdL(gt3,gt2,i2)) 
coup2(8) = Conjg(cplcChacFuSdR(gt3,gt2,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFu cCha Propagator: conj[Sd],conj[Sd]" 
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
!  conj[Sd], Su 
!-------------- 
Do i1=1,6
  Do i2=1,6
Boson4(1) = MSd(i1) 
Boson4(2) = gTSd(i1) 
Boson4(3) = MSu(i2) 
Boson4(4) = gTSu(i2) 
Isum = Isum + 1 
 
resS = 0._dp 
mass(2) = MFu(gt2) 
mass(4) = MFd(gt1) 
mass(3) = -MCha(gt3) 
 
coup2(1) = cplGluFdcSdL(gt1,i1) 
coup2(2) = cplGluFdcSdR(gt1,i1) 
coup2(3) = Conjg(cplcFuGluSuL(gt2,i2)) 
coup2(4) = Conjg(cplcFuGluSuR(gt2,i2))  
coup2(5) = cplcChacFuSdL(gt3,gt2,i1) 
coup2(6) = cplcChacFuSdR(gt3,gt2,i1) 
coup2(7) = Conjg(cplcChaFdcSuL(gt3,gt1,i2)) 
coup2(8) = Conjg(cplcChaFdcSuR(gt3,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFu cCha Propagator: conj[Sd],Su" 
Write(*,*)  "M_in: ",m_in 
Write(*,*)  "mass: ",mass 
Write(*,*)  "coup: ",coup 
gSum(gt1,gt2,gt3,Isum)= 0._dp  
Else 
gSum(gt1,gt2,gt3,Isum)= resS  
End If 
Contribution(gt1,gt2,gt3,Isum)='conj[Sd],Su'
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
mass(3) = -MCha(gt3) 
mass(4) = MFd(gt1) 
 
coup2(1) = cplcFuGluSuL(gt2,i1) 
coup2(2) = cplcFuGluSuR(gt2,i1) 
coup2(3) = Conjg(cplcFuGluSuL(gt2,i2)) 
coup2(4) = Conjg(cplcFuGluSuR(gt2,i2))  
coup2(5) = cplcChaFdcSuL(gt3,gt1,i1) 
coup2(6) = cplcChaFdcSuR(gt3,gt1,i1) 
coup2(7) = Conjg(cplcChaFdcSuL(gt3,gt1,i2)) 
coup2(8) = Conjg(cplcChaFdcSuR(gt3,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fd cFu cCha Propagator: Su,Su" 
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
    Do gt1=1,3
      Do gt2=1,3
        Do gt3=1,5
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:144))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MGlu)**3*g
End Subroutine GluToFdcFucCha 
 
 
Subroutine GluToFucFuChi(iIN,MFu,MChi,MSu,MGlu,cplcFuChiSuL,cplcFuChiSuR,             & 
& cplcFuGluSuL,cplcFuGluSuR,cplChiFucSuL,cplChiFucSuR,cplGluFucSuL,cplGluFucSuR,         & 
& IntegralSs,IntegralSSss,IntegralSSst,NSs,NSSss,NSSst,gTSu,deltaM,epsI,check,           & 
& g,WriteContributions)

Implicit None 
 
Real(dp),Intent(in) :: MFu(3),MChi(10),MSu(6),MGlu

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
mass(1) = MGlu 
 
Isum = 144 
Allocate( gSum(3,3,10, Isum) ) 
Allocate( Contribution(3,3,10, Isum) ) 
gSum = 0._dp  
Contribution = ' ' 
 
Isum = 0 
 
    Do gt1=1,3
      Do gt2=1,3
        Do gt3=1,10
Isum = 0 
 
If(Abs(MGlu).gt.(Abs(MChi(gt3))+Abs(MFu(gt2))+Abs(MFu(gt1)))) Then 
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
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplGluFucSuL(gt1,i1)) 
coup(1) = Conjg(cplGluFucSuR(gt1,i1)) 
coup(4) = Conjg(cplcFuChiSuL(gt2,gt3,i1)) 
coup(3) = Conjg(cplcFuChiSuR(gt2,gt3,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1._dp/2._dp*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fu cFu Chi Propagator: conj[Su]" 
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
mass(3) = -MFu(gt1) 
mass(4) = MChi(gt3) 
 
coup(2) = Conjg(cplcFuGluSuL(gt2,i1)) 
coup(1) = Conjg(cplcFuGluSuR(gt2,i1)) 
coup(4) = Conjg(cplChiFucSuL(gt3,gt1,i1)) 
coup(3) = Conjg(cplChiFucSuR(gt3,gt1,i1))
Call IntegrateScalarSS(Boson2,mass,coup,deltaM,epsI,IntegralSs,NSs,resR, check) 
resR= 1._dp/2._dp*resR ! color factor 
resS = resS + resR 
 
 resD = resD + resS 
If (resD.ne.resD) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fu cFu Chi Propagator: Su" 
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
mass(4) = MChi(gt3) 
 
coup2(1) = cplGluFucSuL(gt1,i1) 
coup2(2) = cplGluFucSuR(gt1,i1) 
coup2(3) = Conjg(cplGluFucSuL(gt1,i2)) 
coup2(4) = Conjg(cplGluFucSuR(gt1,i2))  
coup2(5) = cplcFuChiSuL(gt2,gt3,i1) 
coup2(6) = cplcFuChiSuR(gt2,gt3,i1) 
coup2(7) = Conjg(cplcFuChiSuL(gt2,gt3,i2)) 
coup2(8) = Conjg(cplcFuChiSuR(gt2,gt3,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fu cFu Chi Propagator: conj[Su],conj[Su]" 
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
mass(4) = MFu(gt1) 
mass(3) = -MChi(gt3) 
 
coup2(1) = cplGluFucSuL(gt1,i1) 
coup2(2) = cplGluFucSuR(gt1,i1) 
coup2(3) = Conjg(cplcFuGluSuL(gt2,i2)) 
coup2(4) = Conjg(cplcFuGluSuR(gt2,i2))  
coup2(5) = cplcFuChiSuL(gt2,gt3,i1) 
coup2(6) = cplcFuChiSuR(gt2,gt3,i1) 
coup2(7) = Conjg(cplChiFucSuL(gt3,gt1,i2)) 
coup2(8) = Conjg(cplChiFucSuR(gt3,gt1,i2)) 
Call IntegrateScalarST(Boson4, mass, coup2, deltaM, epsI,IntegralSSst,NSSst, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = -2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fu cFu Chi Propagator: conj[Su],Su" 
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
mass(3) = -MFu(gt1) 
mass(4) = MChi(gt3) 
 
coup2(1) = cplcFuGluSuL(gt2,i1) 
coup2(2) = cplcFuGluSuR(gt2,i1) 
coup2(3) = Conjg(cplcFuGluSuL(gt2,i2)) 
coup2(4) = Conjg(cplcFuGluSuR(gt2,i2))  
coup2(5) = cplChiFucSuL(gt3,gt1,i1) 
coup2(6) = cplChiFucSuR(gt3,gt1,i1) 
coup2(7) = Conjg(cplChiFucSuL(gt3,gt1,i2)) 
coup2(8) = Conjg(cplChiFucSuR(gt3,gt1,i2)) 
Call IntegrateScalarS1S2(Boson4, mass, coup2, deltaM, epsI,IntegralSSss,NSSss, resC, check) 
If (resC.ne.resC) resC = 0._dp
resC = 2._dp*resC 
resC= 1._dp/2._dp*resC ! color factor 
resS = resS + resC 
If (resS.ne.resS) Then 
Write(*,*) "NaN appearing in the following diagrams: " 
Write(*,*) "Glu->Fu cFu Chi Propagator: Su,Su" 
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
    Do gt1=1,3
      Do gt2=1,3
        Do gt3=1,10
g(gt1,gt2,gt3)=Sum(gSum(gt1,gt2,gt3,1:144))
If (g(gt1,gt2,gt3).Lt.0._dp) Then
  Write (ErrCan,*)'Error in Subroutine'//NameOfUnit(Iname)
  g(gt1,gt2,gt3)=0._dp
End If
       End Do 
     End Do 
   End Do 
  g = oo512pi3 / Abs(MGlu)**3*g
End Subroutine GluToFucFuChi 
 
 
End Module Glu3Decays_munuSSM3G 
 
