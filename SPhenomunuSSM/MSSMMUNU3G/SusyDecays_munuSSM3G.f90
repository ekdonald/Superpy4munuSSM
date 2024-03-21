! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:24 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module SusyDecays_munuSSM3G
 
Use Control 
Use DecayFunctions 
Use LoopCouplings_munuSSM3G 
Use CouplingsForDecays_munuSSM3G 
Use Model_Data_munuSSM3G 
Use Mathematics, Only: Li2 
 
 Contains 
 
  
Subroutine SdTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,             & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplAhSdcSd(8,6,6),cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplChiFdcSdL(10,3,6),       & 
& cplChiFdcSdR(10,3,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplhhSdcSd(8,6,6),            & 
& cplHpmSucSd(8,6,6),cplSdcSdVZ(6,6),cplSucSdVWm(6,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'SdTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 6
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.6) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,6

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,6

  End If 
  If (ErrorLevel.Gt.0) then
     Call TerminateProgram
     return
  endif 
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = MSd(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Sd_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,cplAhSdcSd,cplChaFucSdL,cplChaFucSdR,cplChiFdcSdL,          & 
& cplChiFdcSdR,cplGluFdcSdL,cplGluFdcSdR,cplhhSdcSd,cplHpmSucSd,cplSdcSdVZ,              & 
& cplSucSdVWm,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Sd, Ah
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=2, 8
m1out = MSd(gt1)
m2out = MAh(gt2)
coup = cplAhSdcSd(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fu, Cha
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 5
m1out = MFu(gt1)
m2out = MCha(gt2)
coupL = cplChaFucSdL(gt2,gt1,i1)
coupR = cplChaFucSdR(gt2,gt1,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fd, Chi
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 10
m1out = MFd(gt1)
m2out = MChi(gt2)
coupL = cplChiFdcSdL(gt2,gt1,i1)
coupR = cplChiFdcSdR(gt2,gt1,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Glu, Fd
! ----------------------------------------------

 
  Do gt2=1, 3
m1out = MGlu
m2out = MFd(gt2)
coupL = cplGluFdcSdL(gt2,i1)
coupR = cplGluFdcSdR(gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 4._dp/3._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
End Do 
 

 
! ----------------------------------------------
! Sd, hh
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 8
m1out = MSd(gt1)
m2out = Mhh(gt2)
coup = cplhhSdcSd(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Su, Hpm
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=2, 8
m1out = MSu(gt1)
m2out = MHpm(gt2)
coup = cplHpmSucSd(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Sd, VZ
! ----------------------------------------------

 
Do gt1= 1, 6
m1out = MSd(gt1)
m2out = MVZ
coup = cplSdcSdVZ(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Su, VWm
! ----------------------------------------------

 
Do gt1= 1, 6
m1out = MSu(gt1)
m2out = MVWm
coup = cplSucSdVWm(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
End Subroutine SdTwoBodyDecay
 
 
Subroutine SuTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,             & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplAhSucSu(8,6,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),cplcChaFdcSuL(5,3,6),     & 
& cplcChaFdcSuR(5,3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6),cplhhSucSu(8,6,6),            & 
& cplSdcHpmcSu(6,8,6),cplSdcSucVWm(6,6),cplSucSuVZ(6,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'SuTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 6
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.6) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,6

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,6

  End If 
  If (ErrorLevel.Gt.0) then
    Call TerminateProgram 
    return
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = MSu(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Su_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,cplAhSucSu,cplChiFucSuL,cplChiFucSuR,cplcChaFdcSuL,         & 
& cplcChaFdcSuR,cplGluFucSuL,cplGluFucSuR,cplhhSucSu,cplSdcHpmcSu,cplSdcSucVWm,          & 
& cplSucSuVZ,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Su, Ah
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=2, 8
m1out = MSu(gt1)
m2out = MAh(gt2)
coup = cplAhSucSu(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fu, Chi
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 10
m1out = MFu(gt1)
m2out = MChi(gt2)
coupL = cplChiFucSuL(gt2,gt1,i1)
coupR = cplChiFucSuR(gt2,gt1,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Cha], Fd
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=1, 3
m1out = MCha(gt1)
m2out = MFd(gt2)
coupL = cplcChaFdcSuL(gt1,gt2,i1)
coupR = cplcChaFdcSuR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Glu, Fu
! ----------------------------------------------

 
  Do gt2=1, 3
m1out = MGlu
m2out = MFu(gt2)
coupL = cplGluFucSuL(gt2,i1)
coupR = cplGluFucSuR(gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 4._dp/3._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
End Do 
 

 
! ----------------------------------------------
! Su, hh
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 8
m1out = MSu(gt1)
m2out = Mhh(gt2)
coup = cplhhSucSu(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! conj[Hpm], Sd
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2=1, 6
m1out = MHpm(gt1)
m2out = MSd(gt2)
coup = cplSdcHpmcSu(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Sd, conj[VWm]
! ----------------------------------------------

 
Do gt1= 1, 6
m1out = MSd(gt1)
m2out = MVWm
coup = cplSdcSucVWm(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Su, VZ
! ----------------------------------------------

 
Do gt1= 1, 6
m1out = MSu(gt1)
m2out = MVZ
coup = cplSucSuVZ(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
End Subroutine SuTwoBodyDecay
 
 
Subroutine hhTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,             & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplHiggsPP(8),cplHiggsGG(8),cplHiggsZZvirt(8),cplHiggsWWvirt(8),cplAhAhhh(8,8,8),     & 
& cplAhhhhh(8,8,8),cplAhhhVZ(8,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),             & 
& cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),     & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),            & 
& cplhhHpmcVWm(8,8),cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplhhcVWmVWm(8),cplhhVZVZ(8)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Real(dp) :: alpha3 
Iname = Iname + 1 
NameOfUnit(Iname) = 'hhTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 8
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.8) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,8

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,8

  End If 
  If (ErrorLevel.Gt.0) then
    Call TerminateProgram 
    return 
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = Mhh(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_hh_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,cplHiggsPP,cplHiggsGG,cplHiggsZZvirt,cplHiggsWWvirt,        & 
& cplAhAhhh,cplAhhhhh,cplAhhhVZ,cplcChaChahhL,cplcChaChahhR,cplChiChihhL,cplChiChihhR,   & 
& cplcFdFdhhL,cplcFdFdhhR,cplcFuFuhhL,cplcFuFuhhR,cplhhhhhh,cplhhHpmcHpm,cplhhHpmcVWm,   & 
& cplhhSdcSd,cplhhSucSu,cplhhcVWmVWm,cplhhVZVZ,deltaM)

!alpha3 = AlphaSDR(m_in,MSd,MSu,MGlu,MFd,MFu) 
alpha3 = g3running**2/(4._dp*Pi) 
i_count = 1 

 
! ----------------------------------------------
! VP, VP
! ----------------------------------------------

 
m1out = 0.
m2out = 0.
gam = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cplHiggsPP(i1))**2 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! VG, VG
! ----------------------------------------------

 
m1out = 0.
m2out = 0.
gam = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cplHiggsGG(i1))**2 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! VZ, VZ
! ----------------------------------------------

 
m1out = MVZ
m2out = MVZ
If (m_in.le.2._dp*m1out) Then 
coupReal = cplHiggsZZvirt(i1)
Call ScalarToVectorBosonsVR(m_in,m1out,coupReal,gam) 
Else 
gam = 0._dp 
End if 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! VWm, VWm
! ----------------------------------------------

 
m1out = MVWm
m2out = MVWm
If (m_in.le.2._dp*m1out) Then 
coupReal = cplHiggsWWvirt(i1)
Call ScalarToVectorBosonsVR(m_in,m1out,coupReal,gam) 
Else 
gam = 0._dp 
End if 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! Ah, Ah
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2= gt1, 8
m1out = MAh(gt1)
m2out = MAh(gt2)
coup = cplAhAhhh(gt1,gt2,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
If (gt1.ne.gt2) gam = 2._dp*gam 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If (gt1.eq.gt2) Then 
  BRHAA(i1,gt1) = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! hh, Ah
! ----------------------------------------------

 
Do gt1= 1, 8
  Do gt2=2, 8
m1out = Mhh(gt1)
m2out = MAh(gt2)
coup = cplAhhhhh(gt2,i1,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Ah, VZ
! ----------------------------------------------

 
Do gt1= 2, 8
m1out = MAh(gt1)
m2out = MVZ
coup = cplAhhhVZ(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! bar[Cha], Cha
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=1, 5
m1out = MCha(gt1)
m2out = MCha(gt2)
coupL = cplcChaChahhL(gt1,gt2,i1)
coupR = cplcChaChahhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Chi, Chi
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2= gt1, 10
m1out = MChi(gt1)
m2out = MChi(gt2)
coupL = cplChiChihhL(gt1,gt2,i1)
coupR = cplChiChihhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
If (gt1.ne.gt2) gam = 2._dp*gam 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Fd], Fd
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 3
m1out = MFd(gt1)
m2out = MFd(gt2)
coupL = cplcFdFdhhL(gt1,gt2,i1)
coupR = cplcFdFdhhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gam = gam * FFqcd(m1out,m_in,alpha3) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Fu], Fu
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 3
m1out = MFu(gt1)
m2out = MFu(gt2)
coupL = cplcFuFuhhL(gt1,gt2,i1)
coupR = cplcFuFuhhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gam = gam * FFqcd(m1out,m_in,alpha3) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! hh, hh
! ----------------------------------------------

 
Do gt1= 1, 8
  Do gt2= gt1, 8
m1out = Mhh(gt1)
m2out = Mhh(gt2)
coup = cplhhhhhh(i1,gt1,gt2)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
If (gt1.ne.gt2) gam = 2._dp*gam 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If (gt1.eq.gt2) Then 
  BRHHH(i1,gt1) = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! conj[Hpm], Hpm
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2=2, 8
m1out = MHpm(gt1)
m2out = MHpm(gt2)
coup = cplhhHpmcHpm(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Hpm, conj[VWm]
! ----------------------------------------------

 
Do gt1= 2, 8
m1out = MHpm(gt1)
m2out = MVWm
coup = cplhhHpmcVWm(i1,gt1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 2*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! conj[Sd], Sd
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 6
m1out = MSd(gt1)
m2out = MSd(gt2)
coup = cplhhSdcSd(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! conj[Su], Su
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 6
m1out = MSu(gt1)
m2out = MSu(gt2)
coup = cplhhSucSu(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! VWm, conj[VWm]
! ----------------------------------------------

 
m1out = MVWm
m2out = MVWm
coup = cplhhcVWmVWm(i1)
Call ScalarToTwoVectorBosons(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! VZ, VZ
! ----------------------------------------------

 
m1out = MVZ
m2out = MVZ
coup = cplhhVZVZ(i1)
Call ScalarToTwoVectorBosons(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  BRHHH(i1,:) = BRHHH(i1,:)/gT(i1) 
  BRHAA(i1,:) = BRHAA(i1,:)/gT(i1) 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
Contains 
 
  Real(dp) Function FFqcd(mf, mA, alpha_s)
  Implicit None
   Real(dp) , Intent(in) :: mf, mA, alpha_s
   Real(dp) :: fac, beta, beta2, ratio, R_beta_1, Ln_R_beta_1, Ln_beta

   FFqcd = 0._dp
   ratio = mf / mA
   If (ratio.Ge.0.5_dp) Return ! decay is kinematically forbitten

   If (ratio.Ge.0.495_dp) Return ! Coloumb singularity

    beta2 = 1._dp - 4._dp * ratio**2
    beta = Sqrt(beta2)
    
    R_beta_1 = (1. - beta) / (1._dp + beta)
    Ln_beta = Log(beta)
    Ln_R_beta_1 = Log(R_beta_1)

    fac = (3._dp + 34._dp * beta2 - 13._dp * beta**4) / (16._dp * beta**3)  &
      &     * (-Ln_R_beta_1)                                                &
      & + 0.375_dp * (7._dp - beta2) - 3._dp * Log(4._dp/(1._dp - beta**2)) &
      & - 4._dp * Ln_beta                                                   &
      & + (1._dp + beta**2)                                                 &
      &       * ( 4._dp * Li2(R_beta_1) + 2._dp * Li2(- R_beta_1)           &
      &         + Ln_R_beta_1 * ( 3._dp * Log(2._dp/(1._dp + beta))         &
      &                         + 2._dp * Ln_beta )  ) / beta

    fac = fac - 3._dp * Log(ratio)  ! absorb large logarithms in mass

    FFqcd = 1._dp + 5._dp * alpha_s * fac * oo3pi 

  End  Function FFqcd
End Subroutine hhTwoBodyDecay
 
 
Subroutine AhTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,             & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplPseudoHiggsPP(8),cplPseudoHiggsGG(8),cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),            & 
& cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8), & 
& cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),           & 
& cplAhhhhh(8,8,8),cplAhhhVZ(8,8),cplAhHpmcHpm(8,8,8),cplAhHpmcVWm(8,8),cplAhSdcSd(8,6,6),& 
& cplAhSucSu(8,6,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Real(dp) :: alpha3 
Iname = Iname + 1 
NameOfUnit(Iname) = 'AhTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 2 
  i_end = 8
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.8) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,8

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,8

  End If 
  If (ErrorLevel.Gt.0) then
     Call TerminateProgram
     return 
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = MAh(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Ah_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,cplPseudoHiggsPP,cplPseudoHiggsGG,cplAhAhAh,cplAhAhhh,      & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplAhhhhh,cplAhhhVZ,cplAhHpmcHpm,cplAhHpmcVWm,cplAhSdcSd,      & 
& cplAhSucSu,deltaM)

!alpha3 = AlphaSDR(m_in,MSd,MSu,MGlu,MFd,MFu) 
alpha3 = g3running**2/(4._dp*Pi) 
i_count = 1 

 
! ----------------------------------------------
! VP, VP
! ----------------------------------------------

 
m1out = 0.
m2out = 0.
gam = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cplPseudoHiggsPP(i1))**2 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! VG, VG
! ----------------------------------------------

 
m1out = 0.
m2out = 0.
gam = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cplPseudoHiggsGG(i1))**2 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 

 
! ----------------------------------------------
! Ah, Ah
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2= gt1, 8
m1out = MAh(gt1)
m2out = MAh(gt2)
coup = cplAhAhAh(i1,gt1,gt2)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
If (gt1.ne.gt2) gam = 2._dp*gam 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If (gt1.eq.gt2) Then 
  BRAAA(i1,gt1) = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! hh, Ah
! ----------------------------------------------

 
Do gt1= 1, 8
  Do gt2=2, 8
m1out = Mhh(gt1)
m2out = MAh(gt2)
coup = cplAhAhhh(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Cha], Cha
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=1, 5
m1out = MCha(gt1)
m2out = MCha(gt2)
coupL = cplcChaChaAhL(gt1,gt2,i1)
coupR = cplcChaChaAhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Chi, Chi
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2= gt1, 10
m1out = MChi(gt1)
m2out = MChi(gt2)
coupL = cplChiChiAhL(gt1,gt2,i1)
coupR = cplChiChiAhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
If (gt1.ne.gt2) gam = 2._dp*gam 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Fd], Fd
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 3
m1out = MFd(gt1)
m2out = MFd(gt2)
coupL = cplcFdFdAhL(gt1,gt2,i1)
coupR = cplcFdFdAhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gam = gam * FFqcd(m1out,m_in,alpha3) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Fu], Fu
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 3
m1out = MFu(gt1)
m2out = MFu(gt2)
coupL = cplcFuFuAhL(gt1,gt2,i1)
coupR = cplcFuFuAhR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gam = gam * FFqcd(m1out,m_in,alpha3) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! hh, hh
! ----------------------------------------------

 
Do gt1= 1, 8
  Do gt2= gt1, 8
m1out = Mhh(gt1)
m2out = Mhh(gt2)
coup = cplAhhhhh(i1,gt1,gt2)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
If (gt1.ne.gt2) gam = 2._dp*gam 
gPartial(i1,i_count) = 1._dp/2._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If (gt1.eq.gt2) Then 
  BRAHH(i1,gt1) = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! hh, VZ
! ----------------------------------------------

 
Do gt1= 1, 8
m1out = Mhh(gt1)
m2out = MVZ
coup = cplAhhhVZ(i1,gt1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! conj[Hpm], Hpm
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2=2, 8
m1out = MHpm(gt1)
m2out = MHpm(gt2)
coup = cplAhHpmcHpm(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Hpm, conj[VWm]
! ----------------------------------------------

 
Do gt1= 2, 8
m1out = MHpm(gt1)
m2out = MVWm
coup = cplAhHpmcVWm(i1,gt1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 2*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! conj[Sd], Sd
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 6
m1out = MSd(gt1)
m2out = MSd(gt2)
coup = cplAhSdcSd(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! conj[Su], Su
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 6
m1out = MSu(gt1)
m2out = MSu(gt2)
coup = cplAhSucSu(i1,gt2,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
  BRAHH(i1,:) = BRAHH(i1,:)/gT(i1) 
  BRAAA(i1,:) = BRAAA(i1,:)/gT(i1) 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
Contains 
 
 Real(dp) Function FFqcd(mf, mA, alpha_s)
  implicit none
   Real(dp) , Intent(in) :: mf, mA, alpha_s
   Real(dp) :: fac, beta, beta2, ratio, R_beta_1, Ln_R_beta_1, Ln_beta

   FFqcd = 0._dp
   ratio = mf / mA
   if (ratio.ge.0.5_dp) return ! decay is kinematically forbitten

   if (ratio.ge.0.495_dp) return ! Coloumb singularity

    beta2 = 1._dp - 4._dp * ratio**2
    beta = Sqrt(beta2)
    
    R_beta_1 = (1. - beta) / (1._dp + beta)
    Ln_beta = Log(beta)
    Ln_R_beta_1 = Log(R_beta_1)

    fac = (19._dp + 2._dp * beta2 + 3._dp * beta**4) / (16._dp * beta)      &
      &     * (-Ln_R_beta_1)                                                &
      & + 0.375_dp * (7._dp - beta2) - 3._dp * Log(4._dp/(1._dp - beta**2)) &
      & - 4._dp * Ln_beta                                                   &
      & + (1._dp + beta**2)                                                 &
      &       * ( 4._dp * Li2(R_beta_1) + 2._dp * Li2(- R_beta_1)           &
      &         + Ln_R_beta_1 * ( 3._dp * Log(2._dp/(1._dp + beta))         &
      &                         + 2._dp * Ln_beta )  ) / beta
    fac =  fac - 3._dp * Log(ratio)  ! absorb large logarithms in mass

    FFqcd = 1._dp + 5._dp * alpha_s * fac * oo3pi 

  end  Function FFqcd
End Subroutine AhTwoBodyDecay
 
 
Subroutine HpmTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplAhHpmcHpm(8,8,8),cplAhcHpmVWm(8,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),  & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplhhHpmcHpm(8,8,8),cplhhcHpmVWm(8,8),       & 
& cplHpmcHpmVZ(8,8),cplSdcHpmcSu(6,8,6),cplcHpmVWmVZ(8)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'HpmTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 2 
  i_end = 8
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.8) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,8

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,8

  End If 
  If (ErrorLevel.Gt.0) then
     Call TerminateProgram 
     return
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = MHpm(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Hpm_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,               & 
& MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,               & 
& g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,               & 
& me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,cplAhHpmcHpm,cplAhcHpmVWm,cplChiChacHpmL,           & 
& cplChiChacHpmR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplhhHpmcHpm,cplhhcHpmVWm,cplHpmcHpmVZ,     & 
& cplSdcHpmcSu,cplcHpmVWmVZ,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Hpm, Ah
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2=2, 8
m1out = MHpm(gt1)
m2out = MAh(gt2)
coup = cplAhHpmcHpm(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Ah, VWm
! ----------------------------------------------

 
Do gt1= 2, 8
m1out = MAh(gt1)
m2out = MVWm
coup = cplAhcHpmVWm(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Chi, Cha
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2=1, 5
m1out = MChi(gt1)
m2out = MCha(gt2)
coupL = cplChiChacHpmL(gt1,gt2,i1)
coupR = cplChiChacHpmR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If ((gt1.eq.gt2).and.(gt1.eq.3)) Then 
  BR_Htaunu(i1) = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Fu], Fd
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 3
m1out = MFu(gt1)
m2out = MFd(gt2)
coupL = cplcFuFdcHpmL(gt1,gt2,i1)
coupR = cplcFuFdcHpmR(gt1,gt2,i1)
Call ScalarToTwoFermions(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If ((gt1.eq.2).and.(gt2.eq.3)) Then 
  BR_Hcb(i1) = gPartial(i1,i_count) 
End if 
If ((gt1.eq.2).and.(gt2.eq.2)) Then 
  BR_Hcs(i1) = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Hpm, hh
! ----------------------------------------------

 
Do gt1= 2, 8
  Do gt2=1, 8
m1out = MHpm(gt1)
m2out = Mhh(gt2)
coup = cplhhHpmcHpm(gt2,gt1,i1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! hh, VWm
! ----------------------------------------------

 
Do gt1= 1, 8
m1out = Mhh(gt1)
m2out = MVWm
coup = cplhhcHpmVWm(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Hpm, VZ
! ----------------------------------------------

 
Do gt1= 2, 8
m1out = MHpm(gt1)
m2out = MVZ
coup = cplHpmcHpmVZ(gt1,i1)
Call ScalarToScalarVectorBoson(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! conj[Su], Sd
! ----------------------------------------------

 
Do gt1= 1, 6
  Do gt2=1, 6
m1out = MSu(gt1)
m2out = MSd(gt2)
coup = cplSdcHpmcSu(gt2,i1,gt1)
Call ScalarToTwoScalars(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! VWm, VZ
! ----------------------------------------------

 
m1out = MVWm
m2out = MVZ
coup = cplcHpmVWmVZ(i1)
Call ScalarToTwoVectorBosons(m_in,m1out,m2out,coup,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
End Subroutine HpmTwoBodyDecay
 
 
Subroutine GluTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'GluTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 1 
  gT = 0._dp 
  gPartial = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,1

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,1

  End If 
  If (ErrorLevel.Gt.0) then 
    Call TerminateProgram 
    return
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
i1=1 
m_in = MGlu 
Call CouplingsFor_Glu_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,               & 
& MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,               & 
& g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,               & 
& me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,             & 
& cplGluFucSuR,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Fd, conj[Sd]
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 6
m1out = MFd(gt1)
m2out = MSd(gt2)
coupL = cplGluFdcSdL(gt1,gt2)
coupR = cplGluFdcSdR(gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(1,i_count) = 1*gam 
gT = gT + gPartial(1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fu, conj[Su]
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 6
m1out = MFu(gt1)
m2out = MSu(gt2)
coupL = cplGluFucSuL(gt1,gt2)
coupR = cplGluFucSuR(gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(1,i_count) = 1*gam 
gT = gT + gPartial(1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
If ((Present(BR)).And.(gT.Eq.0)) Then 
  BR(1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(1,:) = gPartial(1,:)/gT 
End if 
 
Iname = Iname - 1 
 
End Subroutine GluTwoBodyDecay
 
 
Subroutine FuTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,             & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),      & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplcFuFdcVWmL(3,3),cplcFuFdcVWmR(3,3),       & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),               & 
& cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'FuTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 3
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.3) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,3

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,3

  End If 
  If (ErrorLevel.Gt.0) then
     Call TerminateProgram
     return 
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = MFu(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Fu_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,cplcFuFuAhL,cplcFuFuAhR,cplcFuChiSuL,cplcFuChiSuR,          & 
& cplcFuFdcHpmL,cplcFuFdcHpmR,cplcFuFdcVWmL,cplcFuFdcVWmR,cplcFuFuhhL,cplcFuFuhhR,       & 
& cplcFuFuVZL,cplcFuFuVZR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Fu, Ah
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=2, 8
m1out = MFu(gt1)
m2out = MAh(gt2)
coupL = cplcFuFuAhL(i1,gt1,gt2)
coupR = cplcFuFuAhR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Chi, Su
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2=1, 6
m1out = MChi(gt1)
m2out = MSu(gt2)
coupL = cplcFuChiSuL(i1,gt1,gt2)
coupR = cplcFuChiSuR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fd, conj[Hpm]
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=2, 8
m1out = MFd(gt1)
m2out = MHpm(gt2)
coupL = cplcFuFdcHpmL(i1,gt1,gt2)
coupR = cplcFuFdcHpmR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If (gt1.eq.3) Then 
  BR_tHb = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fd, conj[VWm]
! ----------------------------------------------

 
Do gt1= 1, 3
m1out = MFd(gt1)
m2out = MVWm
coupL = cplcFuFdcVWmL(i1,gt1)
coupR = cplcFuFdcVWmR(i1,gt1)
Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
If (gt1.eq.3) Then 
  BR_tWb = gPartial(i1,i_count) 
End if 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Fu, hh
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 8
m1out = MFu(gt1)
m2out = Mhh(gt2)
coupL = cplcFuFuhhL(i1,gt1,gt2)
coupR = cplcFuFuhhR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Fu, VZ
! ----------------------------------------------

 
Do gt1= 1, 3
m1out = MFu(gt1)
m2out = MVZ
coupL = cplcFuFuVZL(i1,gt1)
coupR = cplcFuFuVZR(i1,gt1)
Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Glu, Su
! ----------------------------------------------

 
  Do gt2=1, 6
m1out = MGlu
m2out = MSu(gt2)
coupL = cplcFuGluSuL(i1,gt2)
coupR = cplcFuGluSuR(i1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 4._dp/3._dp*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
End Do 
 

 
! ----------------------------------------------
! bar[Cha], Sd
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=1, 6
m1out = MCha(gt1)
m2out = MSd(gt2)
coupL = cplcChacFuSdL(gt1,i1,gt2)
coupR = cplcChacFuSdR(gt1,i1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
End Subroutine FuTwoBodyDecay
 
 
Subroutine ChaTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),  & 
& cplcChaChaVZL(5,5),cplcChaChaVZR(5,5),cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),   & 
& cplcChaChiVWmL(5,10),cplcChaChiVWmR(5,10),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),   & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'ChaTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 5
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.5) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,5

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,5

  End If 
  If (ErrorLevel.Gt.0) then
    Call TerminateProgram
    return 
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 
Do i1=i_start,i_end 
m_in = MCha(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Cha_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,               & 
& MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,               & 
& g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,               & 
& me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,          & 
& cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,cplcChaChiHpmR,               & 
& cplcChaChiVWmL,cplcChaChiVWmR,cplcChaFdcSuL,cplcChaFdcSuR,cplcChacFuSdL,               & 
& cplcChacFuSdR,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Cha, Ah
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=2, 8
m1out = MCha(gt1)
m2out = MAh(gt2)
coupL = cplcChaChaAhL(i1,gt1,gt2)
coupR = cplcChaChaAhR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Cha, hh
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=1, 8
m1out = MCha(gt1)
m2out = Mhh(gt2)
coupL = cplcChaChahhL(i1,gt1,gt2)
coupR = cplcChaChahhR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Cha, VZ
! ----------------------------------------------

 
Do gt1= 1, 5
m1out = MCha(gt1)
m2out = MVZ
coupL = cplcChaChaVZL(i1,gt1)
coupR = cplcChaChaVZR(i1,gt1)
Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Chi, Hpm
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2=2, 8
m1out = MChi(gt1)
m2out = MHpm(gt2)
coupL = cplcChaChiHpmL(i1,gt1,gt2)
coupR = cplcChaChiHpmR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! Chi, VWm
! ----------------------------------------------

 
Do gt1= 1, 10
m1out = MChi(gt1)
m2out = MVWm
coupL = cplcChaChiVWmL(i1,gt1)
coupR = cplcChaChiVWmR(i1,gt1)
Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Fd, conj[Su]
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 6
m1out = MFd(gt1)
m2out = MSu(gt2)
coupL = cplcChaFdcSuL(i1,gt1,gt2)
coupR = cplcChaFdcSuR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 

 
! ----------------------------------------------
! bar[Fu], Sd
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 6
m1out = MFu(gt1)
m2out = MSd(gt2)
coupL = cplcChacFuSdL(i1,gt1,gt2)
coupR = cplcChacFuSdR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 3*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 
End Do 
 
Iname = Iname - 1 
 
End Subroutine ChaTwoBodyDecay
 
 
Subroutine ChiTwoBodyDecay(i_in,deltaM,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,            & 
& MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,             & 
& MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,vd,vu,vL,vR,gPartial,gT,BR)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3),MAh(8),MAh2(8),MCha(5),MCha2(5),        & 
& MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),            & 
& MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),         & 
& ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),       & 
& M1,M2,M3,pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),             & 
& ZUL(3,3),ZUR(3,3),ZW(2,2)

Complex(dp) :: cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),& 
& cplChiChacVWmL(10,5),cplChiChacVWmR(10,5),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8), & 
& cplChiChiVZL(10,10),cplChiChiVZR(10,10),cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),     & 
& cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6)

Integer, Intent(in) :: i_in 
Real(dp), Intent(inout) :: gPartial(:,:), gT(:) 
Real(dp), Intent(in) :: deltaM 
Real(dp), Optional, Intent(inout) :: BR(:,:) 
Integer :: i1, i2, i3, i4, i_start, i_end, i_count, gt1, gt2, gt3, gt4 
Real(dp) :: gam, m_in, m1out, m2out, coupReal 
Complex(dp) :: coupC, coupR, coupL, coup 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'ChiTwoBodyDecay'
 
If (i_in.Lt.0) Then 
  i_start = 1 
  i_end = 10
  gT = 0._dp 
  gPartial = 0._dp 
Else If ( (i_in.Ge.1).And.(i_in.Le.10) ) Then 
  i_start = i_in 
  i_end = i_in 
  gT(i_in) = 0._dp 
  gPartial(i_in,:) = 0._dp 
Else 
  If (ErrorLevel.Ge.-1) Then 
     Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,10

     Write(*,*) 'Problem in subroutine '//NameOfUnit(Iname) 
     Write(*,*) 'Value of i_in out of range, (i_in,i_max) = ', i_in,10

  End If 
  If (ErrorLevel.Gt.0) then
     Call TerminateProgram 
     return
  endif
  If (Present(BR)) BR = 0._dp 
  Iname = Iname -1 
  Return 
End If 
 

Do i1=i_start,i_end 
m_in = MChi(i1) 
If (m_in.Eq.0._dp) Cycle 
Call CouplingsFor_Chi_decays_2B(m_in,i1,MAh,MAh2,MCha,MCha2,MChi,MChi2,               & 
& MFd,MFd2,MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,g1,               & 
& g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,               & 
& me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,cplChiChiAhL,cplChiChiAhR,cplChiChacHpmL,           & 
& cplChiChacHpmR,cplChiChacVWmL,cplChiChacVWmR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,   & 
& cplChiChiVZR,cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,deltaM)

i_count = 1 

 
! ----------------------------------------------
! Chi, Ah
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2=2, 8
m1out = MChi(gt1)
m2out = MAh(gt2)
coupL = cplChiChiAhL(i1,gt1,gt2)
coupR = cplChiChiAhR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
 
! ----------------------------------------------
! Cha, conj[Hpm]
! ----------------------------------------------

 
Do gt1= 1, 5
  Do gt2=2, 8
m1out = MCha(gt1)
m2out = MHpm(gt2)
coupL = cplChiChacHpmL(i1,gt1,gt2)
coupR = cplChiChacHpmR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 2*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
 
! ----------------------------------------------
! Cha, conj[VWm]
! ----------------------------------------------

 
Do gt1= 1, 5
m1out = MCha(gt1)
m2out = MVWm
coupL = cplChiChacVWmL(i1,gt1)
coupR = cplChiChacVWmR(i1,gt1)
Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 2*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 


! ----------------------------------------------
! Chi, hh
! ----------------------------------------------

 
Do gt1= 1, 10
  Do gt2=1, 8
m1out = MChi(gt1)
m2out = Mhh(gt2)
coupL = cplChiChihhL(i1,gt1,gt2)
coupR = cplChiChihhR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
 
! ----------------------------------------------
! Chi, VZ
! ----------------------------------------------

 
Do gt1= 1, 10
m1out = MChi(gt1)
m2out = MVZ
coupL = cplChiChiVZL(i1,gt1)
coupR = cplChiChiVZR(i1,gt1)
Call FermionToFermionVectorBoson(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 1*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 

 
! ----------------------------------------------
! Fd, conj[Sd]
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 6
m1out = MFd(gt1)
m2out = MSd(gt2)
coupL = cplChiFdcSdL(i1,gt1,gt2)
coupR = cplChiFdcSdR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 6*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
 
! ----------------------------------------------
! Fu, conj[Su]
! ----------------------------------------------

 
Do gt1= 1, 3
  Do gt2=1, 6
m1out = MFu(gt1)
m2out = MSu(gt2)
coupL = cplChiFucSuL(i1,gt1,gt2)
coupR = cplChiFucSuR(i1,gt1,gt2)
Call FermionToFermionScalar(m_in,m1out,m2out,coupL,coupR,gam) 
gPartial(i1,i_count) = 6*gam 
gT(i1) = gT(i1) + gPartial(i1,i_count) 
i_count = i_count +1 
  End Do 
End Do 
 
If ((Present(BR)).And.(gT(i1).Eq.0)) Then 
  BR(i1,:) = 0._dp 
Else If (Present(BR)) Then 
  BR(i1,:) = gPartial(i1,:)/gT(i1) 
End if 
 

End Do 
 

Iname = Iname - 1 
 
End Subroutine ChiTwoBodyDecay
 
 
Subroutine ScalarToTwoVectorbosonsNew(mS,mV1,mV2,coup,width)
  implicit none
   real(dp), intent(in) :: mS,mV1,mV2
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: coup

   real(dp) :: mSsq,mV1sq, mV2sq,kappa

   if ( abs(mS).le.( abs(mV1)+abs(mV2) ) ) then
    width = 0._dp

   elseif (mV1.eq.0._dp) then
    write(ErrCan,*) 'Server warning, in subroutine ScalarToTwoVectorbosons'
    write(ErrCan,*) 'm_V1 = 0, setting width to 0'
    width = 0._dp

   elseif (mV2.eq.0._dp) then
    write(ErrCan,*) 'Server warning, in subroutine ScalarToTwoVectorbosons'
    write(ErrCan,*) 'm_V2 = 0, setting width to 0'
    width = 0._dp


   elseif (Abs(coup).eq.0._dp) then
    width = 0._dp

   else
    mSsq = mS * mS
    mV1sq = mV1**2
    mV2sq = mV2**2
    kappa = Sqrt( (mSsq-mV1sq-mV2sq)**2 - 4._dp * mV1sq*mV2sq )/(mS**3)
    width = coup**2/(4._dp*mV1sq*mV2sq)*(mV1sq**2 + 10._dp*mV1sq*mV2sq - &
             & 2._dp*(mV1sq+mV2sq)*mSsq + mV2sq**2 +mSsq**2)
    width = oo16Pi*width*kappa

   endif

  End Subroutine ScalarToTwoVectorbosonsNew

Subroutine FermionToFermionVectorBosonMassless(mFin,mFout,mV,coupL,coupR,width)
  implicit none
   real(dp), intent(in) :: mFin,mFout,mV
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: coupL, coupR


   if ( abs(mFin).le.( abs(mFout)+abs(mV) ) ) then
    width = 0._dp

   else

    width = 0.5_dp*oo16Pi*(Abs(coupL)**2 + Abs(coupR)**2)*mFin**3


   endif

  End Subroutine FermionToFermionVectorBosonMassless

End Module SusyDecays_munuSSM3G 
 
