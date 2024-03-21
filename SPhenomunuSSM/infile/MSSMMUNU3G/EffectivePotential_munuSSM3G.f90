Module EffectivePotential_munuSSM3G 
 
Use Control 
Use Couplings_munuSSM3G 
Use LoopFunctions 
Use Mathematics 
Use MathematicsQP 
Use Model_Data_munuSSM3G 
Use StandardModel 
Use SusyMasses_munuSSM3G 
Use EffPotFunctions
Use DerivativesEffPotFunctions
 
Contains 
 
Subroutine CalculateCorrectionsEffPot(ti_ep2L,pi_ep2L,vd,vu,vR,vL,g1,g2,              & 
& g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,              & 
& mv2,mlHd2,M1,M2,M3,kont)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Integer , Intent(inout):: kont 
Integer :: iout 
Real(dp) :: err,h_start,vevs(8) 
Real(dp), Intent(out) :: ti_ep2L(8)  
Real(dp), Intent(out) :: pi_ep2L(8,8)


err2L = 0._dp 
If (.not.PurelyNumericalEffPot) Then 
epsM = 1.0E-8_dp 
epsD = 1.0E-8_dp 
! 2nd deriv. also calculates the 1st deriv. of V
Call SecondDerivativeEffPot2Loop(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,            & 
& Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,kont,             & 
& pi_ep2L,ti_ep2L)

Else 
epsM = 1.0E-6_dp 
epsD = 1.0E-6_dp 
vevs = (/vd,vu,vR,vL/) 
! Calculate 1st (ti_ep) and 2nd derivatives (pi_ep)
ti_ep2L(1) = partialDiff_Ridders(EffPotFunction2Loop,vevs,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(2) = partialDiff_Ridders(EffPotFunction2Loop,vevs,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(3) = partialDiff_Ridders(EffPotFunction2Loop,vevs,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(4) = partialDiff_Ridders(EffPotFunction2Loop,vevs,4,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(5) = partialDiff_Ridders(EffPotFunction2Loop,vevs,5,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(6) = partialDiff_Ridders(EffPotFunction2Loop,vevs,6,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(7) = partialDiff_Ridders(EffPotFunction2Loop,vevs,7,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
ti_ep2L(8) = partialDiff_Ridders(EffPotFunction2Loop,vevs,8,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(1,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,1,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(2,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,2,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(2,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,2,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(3,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,3,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(3,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,3,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(3,3) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,3,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(4,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,4,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(4,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,4,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(4,3) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,4,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(4,4) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,4,4,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(5,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,5,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(5,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,5,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(5,3) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,5,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(5,4) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,5,4,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(5,5) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,5,5,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(6,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,6,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(6,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,6,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(6,3) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,6,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(6,4) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,6,4,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(6,5) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,6,5,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(6,6) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,6,6,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,3) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,4) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,4,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,5) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,5,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,6) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,6,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(7,7) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,7,7,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,1) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,1,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,2) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,2,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,3) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,3,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,4) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,4,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,5) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,5,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,6) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,6,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,7) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,7,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(8,8) = partialDiffXY_Ridders(EffPotFunction2Loop,vevs,8,8,8,err,h_start,iout) 
If (err.gt.err2L) err2L = err 
pi_ep2L(1,2)=pi_ep2L(2,1)
pi_ep2L(1,3)=pi_ep2L(3,1)
pi_ep2L(2,3)=pi_ep2L(3,2)
pi_ep2L(1,4)=pi_ep2L(4,1)
pi_ep2L(2,4)=pi_ep2L(4,2)
pi_ep2L(3,4)=pi_ep2L(4,3)
pi_ep2L(1,5)=pi_ep2L(5,1)
pi_ep2L(2,5)=pi_ep2L(5,2)
pi_ep2L(3,5)=pi_ep2L(5,3)
pi_ep2L(4,5)=pi_ep2L(5,4)
pi_ep2L(1,6)=pi_ep2L(6,1)
pi_ep2L(2,6)=pi_ep2L(6,2)
pi_ep2L(3,6)=pi_ep2L(6,3)
pi_ep2L(4,6)=pi_ep2L(6,4)
pi_ep2L(5,6)=pi_ep2L(6,5)
pi_ep2L(1,7)=pi_ep2L(7,1)
pi_ep2L(2,7)=pi_ep2L(7,2)
pi_ep2L(3,7)=pi_ep2L(7,3)
pi_ep2L(4,7)=pi_ep2L(7,4)
pi_ep2L(5,7)=pi_ep2L(7,5)
pi_ep2L(6,7)=pi_ep2L(7,6)
pi_ep2L(1,8)=pi_ep2L(8,1)
pi_ep2L(2,8)=pi_ep2L(8,2)
pi_ep2L(3,8)=pi_ep2L(8,3)
pi_ep2L(4,8)=pi_ep2L(8,4)
pi_ep2L(5,8)=pi_ep2L(8,5)
pi_ep2L(6,8)=pi_ep2L(8,6)
pi_ep2L(7,8)=pi_ep2L(8,7)
End If 
Contains 

Real(dp) Function EffPotFunction(vevs) 
  ! a wrapping function to be passed to numerical differentiation 
  Implicit None 
  Real(dp), Intent(in) :: vevs(8) 
  Real(dp) :: effPot 
Call CalculateEffPot(vevs(1),vevs(2),vevs(3:5),vevs(6:8),g1,g2,g3,Yd,Ye,              & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,kont,effPot)

  EffPotFunction = effPot 
  End Function 

Real(dp) Function EffPotFunction2Loop(vevs) 
  ! a wrapping function to be passed to numerical differentiation 
  Implicit None 
  Real(dp), Intent(in) :: vevs(8) 
  Real(dp) :: effPot2L 
Call CalculateEffPot2Loop(vevs(1),vevs(2),vevs(3:5),vevs(6:8),g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,kont,effPot2L)

  EffPotFunction2Loop = effPot2L 
  End Function 

End subroutine CalculateCorrectionsEffPot 


Subroutine CalculateEffPot(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,               & 
& Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,kont,effPot)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Integer, Intent(inout):: kont
Integer :: i 
Real(dp) :: effpot,Qscale,result,temp


Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

Qscale= getRenormalizationScale()
result=0._dp
temp=0._dp
! sum over real scalars (color *[2 if complex]) 
Do i=1,6
temp=temp+(3*2) * h_SMartin(MSd2(i),Qscale) ! Sd
End Do
Do i=1,6
temp=temp+(3*2) * h_SMartin(MSu2(i),Qscale) ! Su
End Do
Do i=1,8
temp=temp+          h_SMartin(Mhh2(i),Qscale) ! hh
End Do
Do i=2,8
temp=temp+          h_SMartin(MAh2(i),Qscale) ! Ah
End Do
Do i=2,8
temp=temp+(2)    * h_SMartin(MHpm2(i),Qscale) ! Hpm
End Do
result=(+1)*temp !scalars

temp=0._dp
! sum over two-component fermions (*color [*2 if Dirac Fermion]) 
Do i=1,5
temp=temp+          h_SMartin(MCha2(i),Qscale)*2 ! Cha
End Do
Do i=1,10
temp=temp+          h_SMartin(MChi2(i),Qscale)*1 ! Chi
End Do
Do i=1,3
temp=temp+(3)   * h_SMartin(MFd2(i),Qscale)*2 ! Fd
End Do
Do i=1,3
temp=temp+(3)   * h_SMartin(MFu2(i),Qscale)*2 ! Fu
End Do
temp=temp+(8)   * h_SMartin(MGlu2,Qscale)*1 ! Glu
result=result+(-2)*temp ! fermions

temp=0._dp
! sum over real vectors (color *[2 if complex]) 
temp=temp+(8)   * h_SMartin(0._dp,Qscale) ! VG
temp=temp+          h_SMartin(0._dp,Qscale) ! VP
temp=temp+          h_SMartin(MVZ2,Qscale) ! VZ
temp=temp+(2)    * h_SMartin(MVWm2,Qscale) ! VWm
result=result+(+3)*temp ! vectors

effPot = result * oo16pi2

End Subroutine CalculateEffPot 



 
 
Subroutine CalculateEffPot2Loop(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,             & 
& Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,kont,effPot2L)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Integer, Intent(inout):: kont
Real(dp), Intent(out) :: effpot2L
Integer :: i,i1,i2,i3,includeGhosts,NrContr 
Integer :: NrContr1,NrContr2 !nr of contributing diagrams
Real(dp) :: Qscale,result,colorfactor,coeff,coeffbar
Complex(dp) :: temp,coupx,coupxbar,coup1,coup2,coup1L,coup1R,coup2L,coup2R
Complex(dp) :: dcoupx,dcoupxbar,dcoup1,dcoup2,dcoup1L,dcoup1R,dcoup2L,dcoup2R
Complex(dp) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

Real(dp) :: results1(33),results2(17)


Call TreeMassesEffPot(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,               & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,.True.,kont)

Call CouplingsForEffPot3(lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZA,ZH,Ye,Te,               & 
& ZP,Yd,Td,ZD,Yu,Tu,ZU,g3,ZER,ZEL,UV,ZDL,ZDR,ZUL,ZUR,pG,cplAhAhAh,cplAhAhhh,             & 
& cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,cplhhSdcSd,        & 
& cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplcChaChaAhL,     & 
& cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,cplcFuFuAhL,           & 
& cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,cplcChaChahhL,     & 
& cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,cplChiFdcSdL,        & 
& cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,cplcFdChiSdL,     & 
& cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,cplcFdFdhhL,          & 
& cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplGluFucSuL,      & 
& cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdGluSdL,           & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,cplcFdFdVGL,        & 
& cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR)

Call CouplingsForEffPot4(lam,Yv,kap,ZA,ZH,Ye,ZP,Yd,ZD,Yu,ZU,g3,cplAhAhAhAh,           & 
& cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,cplhhhhHpmcHpm,       & 
& cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,cplHpmSucHpmcSu,           & 
& cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)


Qscale = getRenormalizationScale()
result=0._dp
results1=0._dp
results2=0._dp
temp=0._dp
! ----- Topology1 (sunrise): diagrams w. 3 Particles and 2 Vertices

! ----- diagrams of type SSS, 11 ------ 
! ---- Ah,Ah,Ah ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhAhAh(i1,i2,i3)
coup2 = cplAhAhAh(i1,i2,i3)
colorfactor=1
temp=temp+colorfactor*1._dp/12._dp*abs(coup1)**2*Fep_SSS(MAh2(i1),MAh2(i2),MAh2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Ah, Ah]' 
    End Do
  End Do
End Do
results1(1)=temp
! ---- Ah,Ah,hh ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhAhhh(i1,i2,i3)
coup2 = cplAhAhhh(i1,i2,i3)
colorfactor=1
temp=temp+colorfactor*0.25_dp*abs(coup1)**2*Fep_SSS(MAh2(i1),MAh2(i2),Mhh2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Ah, hh]' 
    End Do
  End Do
End Do
results1(2)=temp
! ---- Ah,hh,hh ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhhhhh(i1,i2,i3)
coup2 = cplAhhhhh(i1,i2,i3)
colorfactor=1
temp=temp+colorfactor*0.25_dp*abs(coup1)**2*Fep_SSS(MAh2(i1),Mhh2(i2),Mhh2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, hh, hh]' 
    End Do
  End Do
End Do
results1(3)=temp
! ---- Ah,Hpm,conj[Hpm] ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhHpmcHpm(i1,i2,i3)
coup2 = cplAhHpmcHpm(i1,i3,i2)
colorfactor=1
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSS(MAh2(i1),MHpm2(i2),MHpm2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Hpm, conj[Hpm]]' 
    End Do
  End Do
End Do
results1(4)=temp
! ---- Ah,Sd,conj[Sd] ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplAhSdcSd(i1,i2,i3)
coup2 = cplAhSdcSd(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSS(MAh2(i1),MSd2(i2),MSd2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Sd, conj[Sd]]' 
    End Do
  End Do
End Do
results1(5)=temp
! ---- Ah,Su,conj[Su] ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplAhSucSu(i1,i2,i3)
coup2 = cplAhSucSu(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSS(MAh2(i1),MSu2(i2),MSu2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Su, conj[Su]]' 
    End Do
  End Do
End Do
results1(6)=temp
! ---- hh,hh,hh ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplhhhhhh(i1,i2,i3)
coup2 = cplhhhhhh(i1,i2,i3)
colorfactor=1
temp=temp+colorfactor*1._dp/12._dp*abs(coup1)**2*Fep_SSS(Mhh2(i1),Mhh2(i2),Mhh2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, hh, hh]' 
    End Do
  End Do
End Do
results1(7)=temp
! ---- hh,Hpm,conj[Hpm] ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplhhHpmcHpm(i1,i2,i3)
coup2 = cplhhHpmcHpm(i1,i3,i2)
colorfactor=1
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSS(Mhh2(i1),MHpm2(i2),MHpm2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, Hpm, conj[Hpm]]' 
    End Do
  End Do
End Do
results1(8)=temp
! ---- hh,Sd,conj[Sd] ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplhhSdcSd(i1,i2,i3)
coup2 = cplhhSdcSd(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSS(Mhh2(i1),MSd2(i2),MSd2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, Sd, conj[Sd]]' 
    End Do
  End Do
End Do
results1(9)=temp
! ---- hh,Su,conj[Su] ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplhhSucSu(i1,i2,i3)
coup2 = cplhhSucSu(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSS(Mhh2(i1),MSu2(i2),MSu2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, Su, conj[Su]]' 
    End Do
  End Do
End Do
results1(10)=temp
! ---- Sd,conj[Hpm],conj[Su] ----
temp=0._dp
Do i1=1,6
 Do i2=1,8
    Do i3=1,6
coup1 = cplSdcHpmcSu(i1,i2,i3)
coup2 = cplHpmSucSd(i2,i3,i1)
colorfactor=3
temp=temp+colorfactor*1._dp*abs(coup1)**2*Fep_SSS(MSd2(i1),MHpm2(i2),MSu2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Sd, conj[Hpm], conj[Su]]' 
    End Do
  End Do
End Do
results1(11)=temp
! ----- diagrams of type FFS, 16 ------ 
! ---- Ah,Cha,bar[Cha] ----
temp=0._dp
Do i1=1,8
 Do i2=1,5
    Do i3=1,5
coup1L = cplcChaChaAhL(i3,i2,i1)
coup1R = cplcChaChaAhR(i3,i2,i1)
coup2L = cplcChaChaAhL(i2,i3,i1)
coup2R = cplcChaChaAhR(i2,i3,i1)
colorfactor=1
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MCha2(i3),MCha2(i2),MAh2(i1),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Cha, bar[Cha]]' 
temp=temp+colorfactor*0.5_dp*2*Real(coup1L*conjg(coup1R),dp)*MCha(i2)*MCha(i3)*Fep_FFSbar(MCha2(i3),MCha2(i2),MAh2(i1),Qscale)
    End Do
  End Do
End Do
results1(12)=temp
! ---- Ah,Chi,Chi ----
temp=0._dp
Do i1=1,8
 Do i2=1,10
    Do i3=1,10
coup1L = cplChiChiAhL(i2,i3,i1)
coup1R = cplChiChiAhR(i2,i3,i1)
coup2L = cplChiChiAhL(i2,i3,i1)
coup2R = cplChiChiAhR(i2,i3,i1)
colorfactor=1
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2)*Fep_FFS(MChi2(i3),MChi2(i2),MAh2(i1),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Chi, Chi]' 
temp=temp+colorfactor*0.5_dp*Real(coup1L**2,dp)*MChi(i2)*MChi(i3)*Fep_FFSbar(MChi2(i3),MChi2(i2),MAh2(i1),Qscale)
    End Do
  End Do
End Do
results1(13)=temp
! ---- Ah,Fd,bar[Fd] ----
temp=0._dp
Do i1=1,8
 Do i2=1,3
    Do i3=1,3
coup1L = cplcFdFdAhL(i3,i2,i1)
coup1R = cplcFdFdAhR(i3,i2,i1)
coup2L = cplcFdFdAhL(i2,i3,i1)
coup2R = cplcFdFdAhR(i2,i3,i1)
colorfactor=3
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFd2(i3),MFd2(i2),MAh2(i1),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Fd, bar[Fd]]' 
temp=temp+colorfactor*0.5_dp*2*Real(coup1L*conjg(coup1R),dp)*MFd(i2)*MFd(i3)*Fep_FFSbar(MFd2(i3),MFd2(i2),MAh2(i1),Qscale)
    End Do
  End Do
End Do
results1(14)=temp
! ---- Ah,Fu,bar[Fu] ----
temp=0._dp
Do i1=1,8
 Do i2=1,3
    Do i3=1,3
coup1L = cplcFuFuAhL(i3,i2,i1)
coup1R = cplcFuFuAhR(i3,i2,i1)
coup2L = cplcFuFuAhL(i2,i3,i1)
coup2R = cplcFuFuAhR(i2,i3,i1)
colorfactor=3
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFu2(i3),MFu2(i2),MAh2(i1),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Fu, bar[Fu]]' 
temp=temp+colorfactor*0.5_dp*2*Real(coup1L*conjg(coup1R),dp)*MFu(i2)*MFu(i3)*Fep_FFSbar(MFu2(i3),MFu2(i2),MAh2(i1),Qscale)
    End Do
  End Do
End Do
results1(15)=temp
! ---- Cha,hh,bar[Cha] ----
temp=0._dp
Do i1=1,5
 Do i2=1,8
    Do i3=1,5
coup1L = cplcChaChahhL(i3,i1,i2)
coup1R = cplcChaChahhR(i3,i1,i2)
coup2L = cplcChaChahhL(i1,i3,i2)
coup2R = cplcChaChahhR(i1,i3,i2)
colorfactor=1
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MCha2(i3),MCha2(i1),Mhh2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Cha, hh, bar[Cha]]' 
temp=temp+colorfactor*0.5_dp*2*Real(coup1L*conjg(coup1R),dp)*MCha(i1)*MCha(i3)*Fep_FFSbar(MCha2(i3),MCha2(i1),Mhh2(i2),Qscale)
    End Do
  End Do
End Do
results1(16)=temp
! ---- Chi,Chi,hh ----
temp=0._dp
Do i1=1,10
 Do i2=1,10
    Do i3=1,8
coup1L = cplChiChihhL(i1,i2,i3)
coup1R = cplChiChihhR(i1,i2,i3)
coup2L = cplChiChihhL(i1,i2,i3)
coup2R = cplChiChihhR(i1,i2,i3)
colorfactor=1
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2)*Fep_FFS(MChi2(i2),MChi2(i1),Mhh2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Chi, hh]' 
temp=temp+colorfactor*0.5_dp*Real(coup1L**2,dp)*MChi(i1)*MChi(i2)*Fep_FFSbar(MChi2(i2),MChi2(i1),Mhh2(i3),Qscale)
    End Do
  End Do
End Do
results1(17)=temp
! ---- Chi,Hpm,bar[Cha] ----
temp=0._dp
Do i1=1,10
 Do i2=1,8
    Do i3=1,5
coup1L = cplcChaChiHpmL(i3,i1,i2)
coup1R = cplcChaChiHpmR(i3,i1,i2)
coup2L = cplChiChacHpmL(i1,i3,i2)
coup2R = cplChiChacHpmR(i1,i3,i2)
colorfactor=1
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MCha2(i3),MChi2(i1),MHpm2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Hpm, bar[Cha]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MChi(i1)*MCha(i3)*Fep_FFSbar(MCha2(i3),MChi2(i1),MHpm2(i2),Qscale)
    End Do
  End Do
End Do
results1(18)=temp
! ---- Chi,Sd,bar[Fd] ----
temp=0._dp
Do i1=1,10
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFdChiSdL(i3,i1,i2)
coup1R = cplcFdChiSdR(i3,i1,i2)
coup2L = cplChiFdcSdL(i1,i3,i2)
coup2R = cplChiFdcSdR(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFd2(i3),MChi2(i1),MSd2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Sd, bar[Fd]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MChi(i1)*MFd(i3)*Fep_FFSbar(MFd2(i3),MChi2(i1),MSd2(i2),Qscale)
    End Do
  End Do
End Do
results1(19)=temp
! ---- Chi,Su,bar[Fu] ----
temp=0._dp
Do i1=1,10
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFuChiSuL(i3,i1,i2)
coup1R = cplcFuChiSuR(i3,i1,i2)
coup2L = cplChiFucSuL(i1,i3,i2)
coup2R = cplChiFucSuR(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFu2(i3),MChi2(i1),MSu2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Su, bar[Fu]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MChi(i1)*MFu(i3)*Fep_FFSbar(MFu2(i3),MChi2(i1),MSu2(i2),Qscale)
    End Do
  End Do
End Do
results1(20)=temp
! ---- Fd,hh,bar[Fd] ----
temp=0._dp
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFdFdhhL(i3,i1,i2)
coup1R = cplcFdFdhhR(i3,i1,i2)
coup2L = cplcFdFdhhL(i1,i3,i2)
coup2R = cplcFdFdhhR(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFd2(i3),MFd2(i1),Mhh2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fd, hh, bar[Fd]]' 
temp=temp+colorfactor*0.5_dp*2*Real(coup1L*conjg(coup1R),dp)*MFd(i1)*MFd(i3)*Fep_FFSbar(MFd2(i3),MFd2(i1),Mhh2(i2),Qscale)
    End Do
  End Do
End Do
results1(21)=temp
! ---- Fd,bar[Cha],conj[Su] ----
temp=0._dp
Do i1=1,3
 Do i2=1,5
    Do i3=1,6
coup1L = cplcChaFdcSuL(i2,i1,i3)
coup1R = cplcChaFdcSuR(i2,i1,i3)
coup2L = cplcFdChaSuL(i1,i2,i3)
coup2R = cplcFdChaSuR(i1,i2,i3)
colorfactor=3
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MCha2(i2),MFd2(i1),MSu2(i3),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fd, bar[Cha], conj[Su]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MFd(i1)*MCha(i2)*Fep_FFSbar(MCha2(i2),MFd2(i1),MSu2(i3),Qscale)
    End Do
  End Do
End Do
results1(22)=temp
! ---- Fu,hh,bar[Fu] ----
temp=0._dp
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFuFuhhL(i3,i1,i2)
coup1R = cplcFuFuhhR(i3,i1,i2)
coup2L = cplcFuFuhhL(i1,i3,i2)
coup2R = cplcFuFuhhR(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFu2(i3),MFu2(i1),Mhh2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fu, hh, bar[Fu]]' 
temp=temp+colorfactor*0.5_dp*2*Real(coup1L*conjg(coup1R),dp)*MFu(i1)*MFu(i3)*Fep_FFSbar(MFu2(i3),MFu2(i1),Mhh2(i2),Qscale)
    End Do
  End Do
End Do
results1(23)=temp
! ---- Fu,Hpm,bar[Fd] ----
temp=0._dp
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFdFuHpmL(i3,i1,i2)
coup1R = cplcFdFuHpmR(i3,i1,i2)
coup2L = cplcFuFdcHpmL(i1,i3,i2)
coup2R = cplcFuFdcHpmR(i1,i3,i2)
colorfactor=3
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFd2(i3),MFu2(i1),MHpm2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fu, Hpm, bar[Fd]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MFu(i1)*MFd(i3)*Fep_FFSbar(MFd2(i3),MFu2(i1),MHpm2(i2),Qscale)
    End Do
  End Do
End Do
results1(24)=temp
! ---- Glu,Sd,bar[Fd] ----
temp=0._dp
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFdGluSdL(i3,i2)
coup1R = cplcFdGluSdR(i3,i2)
coup2L = cplGluFdcSdL(i3,i2)
coup2R = cplGluFdcSdR(i3,i2)
colorfactor=4
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFd2(i3),MGlu2,MSd2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Glu, Sd, bar[Fd]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MGlu*MFd(i3)*Fep_FFSbar(MFd2(i3),MGlu2,MSd2(i2),Qscale)
    End Do
  End Do
results1(25)=temp
! ---- Glu,Su,bar[Fu] ----
temp=0._dp
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFuGluSuL(i3,i2)
coup1R = cplcFuGluSuR(i3,i2)
coup2L = cplGluFucSuL(i3,i2)
coup2R = cplGluFucSuR(i3,i2)
colorfactor=4
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFu2(i3),MGlu2,MSu2(i2),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Glu, Su, bar[Fu]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MGlu*MFu(i3)*Fep_FFSbar(MFu2(i3),MGlu2,MSu2(i2),Qscale)
    End Do
  End Do
results1(26)=temp
! ---- Sd,bar[Cha],bar[Fu] ----
temp=0._dp
Do i1=1,6
 Do i2=1,5
    Do i3=1,3
coup1L = cplcChacFuSdL(i2,i3,i1)
coup1R = cplcChacFuSdR(i2,i3,i1)
coup2L = cplChaFucSdL(i2,i3,i1)
coup2R = cplChaFucSdR(i2,i3,i1)
colorfactor=3
temp=temp+colorfactor*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFS(MFu2(i3),MCha2(i2),MSd2(i1),Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Sd, bar[Cha], bar[Fu]]' 
temp=temp+colorfactor*2*Real(coup1L*conjg(coup1R),dp)*MCha(i2)*MFu(i3)*Fep_FFSbar(MFu2(i3),MCha2(i2),MSd2(i1),Qscale)
    End Do
  End Do
End Do
results1(27)=temp
! ----- diagrams of type SSV, 2 ------ 
! ---- Sd,VG,conj[Sd] ----
temp=0._dp
Do i1=1,6
    Do i3=1,6
coup1 = cplSdcSdVG(i1,i3)
coup2 = cplSdcSdVG(i3,i1)
colorfactor=4
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSV(MSd2(i3),MSd2(i1),0._dp,Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSV C[Sd, VG, conj[Sd]]' 
    End Do
End Do
results1(28)=temp
! ---- Su,VG,conj[Su] ----
temp=0._dp
Do i1=1,6
    Do i3=1,6
coup1 = cplSucSuVG(i1,i3)
coup2 = cplSucSuVG(i3,i1)
colorfactor=4
temp=temp+colorfactor*0.5_dp*abs(coup1)**2*Fep_SSV(MSu2(i3),MSu2(i1),0._dp,Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSV C[Su, VG, conj[Su]]' 
    End Do
End Do
results1(29)=temp
! ----- diagrams of type FFV, 3 ------ 
! ---- Fd,VG,bar[Fd] ----
temp=0._dp
Do i1=1,3
    Do i3=1,3
coup1L = cplcFdFdVGL(i3,i1)
coup1R = cplcFdFdVGR(i3,i1)
coup2L = cplcFdFdVGL(i1,i3)
coup2R = cplcFdFdVGR(i1,i3)
colorfactor=4
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFV(MFd2(i3),MFd2(i1),0._dp,Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFV C[Fd, VG, bar[Fd]]' 
temp=temp+colorfactor*0.5_dp*2*Real(-coup1L*conjg(coup1R),dp)*MFd(i1)*MFd(i3)*Fep_FFVbar(MFd2(i3),MFd2(i1),0._dp,Qscale)
    End Do
End Do
results1(30)=temp
! ---- Fu,VG,bar[Fu] ----
temp=0._dp
Do i1=1,3
    Do i3=1,3
coup1L = cplcFuFuVGL(i3,i1)
coup1R = cplcFuFuVGR(i3,i1)
coup2L = cplcFuFuVGL(i1,i3)
coup2R = cplcFuFuVGR(i1,i3)
colorfactor=4
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2+abs(coup1R)**2)*Fep_FFV(MFu2(i3),MFu2(i1),0._dp,Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFV C[Fu, VG, bar[Fu]]' 
temp=temp+colorfactor*0.5_dp*2*Real(-coup1L*conjg(coup1R),dp)*MFu(i1)*MFu(i3)*Fep_FFVbar(MFu2(i3),MFu2(i1),0._dp,Qscale)
    End Do
End Do
results1(31)=temp
! ---- Glu,Glu,VG ----
temp=0._dp
coup1L = cplGluGluVGL
coup1R = cplGluGluVGR
coup2L = cplGluGluVGL
coup2R = cplGluGluVGR
colorfactor=24
temp=temp+colorfactor*0.5_dp*(abs(coup1L)**2)*Fep_FFV(MGlu2,MGlu2,0._dp,Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFV C[Glu, Glu, VG]' 
temp=temp+colorfactor*0.5_dp*Real(coup1L**2,dp)*MGlu*MGlu*Fep_FFVbar(MGlu2,MGlu2,0._dp,Qscale)
results1(32)=temp
! ----- diagrams of type VVV, 1 ------ 
! ---- VG,VG,VG ----
temp=0._dp
coup1 = cplVGVGVG
coup2 = cplVGVGVG
colorfactor=24
temp=temp+colorfactor*1._dp/12._dp*(coup1)**2*Fep_gauge(0._dp,0._dp,0._dp,Qscale)
 if (.not.(temp.eq.temp))  write(*,*) 'NaN at VVV C[VG, VG, VG]' 
results1(33)=temp
! ----- Topology2: diagrams w. 2 Particles and 1 Vertex

! ----- diagrams of type SS, 15 ------ 
! ---- Ah,Ah ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhAhAh(i1,i1,i2,i2)
temp=temp+(-1._dp/8._dp)*(-coup1)*Fep_SS(MAh2(i1),MAh2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Ah, Ah]' 
  End Do
End Do
results2(1)=temp
! ---- Ah,hh ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhhhhh(i1,i1,i2,i2)
temp=temp+(-0.25_dp)*(-coup1)*Fep_SS(MAh2(i1),Mhh2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, hh, hh]' 
  End Do
End Do
results2(2)=temp
! ---- Ah,Hpm ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhHpmcHpm(i1,i1,i2,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(MAh2(i1),MHpm2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Hpm, conj[Hpm]]' 
  End Do
End Do
results2(3)=temp
! ---- Ah,Sd ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
coup1 = cplAhAhSdcSd(i1,i1,i2,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(MAh2(i1),MSd2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Sd, conj[Sd]]' 
  End Do
End Do
results2(4)=temp
! ---- Ah,Su ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
coup1 = cplAhAhSucSu(i1,i1,i2,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(MAh2(i1),MSu2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Su, conj[Su]]' 
  End Do
End Do
results2(5)=temp
! ---- hh,hh ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
coup1 = cplhhhhhhhh(i1,i1,i2,i2)
temp=temp+(-1._dp/8._dp)*(-coup1)*Fep_SS(Mhh2(i1),Mhh2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, hh, hh]' 
  End Do
End Do
results2(6)=temp
! ---- hh,Hpm ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
coup1 = cplhhhhHpmcHpm(i1,i1,i2,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(Mhh2(i1),MHpm2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, Hpm, conj[Hpm]]' 
  End Do
End Do
results2(7)=temp
! ---- hh,Sd ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
coup1 = cplhhhhSdcSd(i1,i1,i2,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(Mhh2(i1),MSd2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, Sd, conj[Sd]]' 
  End Do
End Do
results2(8)=temp
! ---- hh,Su ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
coup1 = cplhhhhSucSu(i1,i1,i2,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(Mhh2(i1),MSu2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, Su, conj[Su]]' 
  End Do
End Do
results2(9)=temp
! ---- Hpm,Hpm ----
temp=0._dp
Do i1=1,8
 Do i2=1,8
coup1 = cplHpmHpmcHpmcHpm(i1,i2,i1,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(MHpm2(i1),MHpm2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Hpm, Hpm, conj[Hpm], conj[Hpm]]' 
  End Do
End Do
results2(10)=temp
! ---- Hpm,Sd ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
coup1 = cplHpmSdcHpmcSd(i1,i2,i1,i2)
temp=temp+(-1._dp)*(-coup1)*Fep_SS(MHpm2(i1),MSd2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Hpm, Sd, conj[Hpm], conj[Sd]]' 
  End Do
End Do
results2(11)=temp
! ---- Hpm,Su ----
temp=0._dp
Do i1=1,8
 Do i2=1,6
coup1 = cplHpmSucHpmcSu(i1,i2,i1,i2)
temp=temp+(-1._dp)*(-coup1)*Fep_SS(MHpm2(i1),MSu2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Hpm, Su, conj[Hpm], conj[Su]]' 
  End Do
End Do
results2(12)=temp
! ---- Sd,Sd ----
temp=0._dp
Do i1=1,6
 Do i2=1,6
coup1 = cplSdSdcSdcSd(i1,i2,i1,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(MSd2(i1),MSd2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Sd, Sd, conj[Sd], conj[Sd]]' 
  End Do
End Do
results2(13)=temp
! ---- Sd,Su ----
temp=0._dp
Do i1=1,6
 Do i2=1,6
coup1 = cplSdSucSdcSu(i1,i2,i1,i2)
temp=temp+(-1._dp)*(-coup1)*Fep_SS(MSd2(i1),MSu2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Sd, Su, conj[Sd], conj[Su]]' 
  End Do
End Do
results2(14)=temp
! ---- Su,Su ----
temp=0._dp
Do i1=1,6
 Do i2=1,6
coup1 = cplSuSucSucSu(i1,i2,i1,i2)
temp=temp+(-0.5_dp)*(-coup1)*Fep_SS(MSu2(i1),MSu2(i2),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Su, Su, conj[Su], conj[Su]]' 
  End Do
End Do
results2(15)=temp
! ----- diagrams of type VS, 2 ------ 
! ---- Sd,VG ----
temp=0._dp
Do i1=1,6
coup1 = cplSdcSdVGVG(i1,i1)
temp=temp+0.25_dp*coup1*Fep_VS(0._dp,MSd2(i1),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at VS C[Sd, VG, VG, conj[Sd]]' 
End Do
results2(16)=temp
! ---- Su,VG ----
temp=0._dp
Do i1=1,6
coup1 = cplSucSuVGVG(i1,i1)
temp=temp+0.25_dp*coup1*Fep_VS(0._dp,MSu2(i1),Qscale)
if (.not.(temp.eq.temp))  write(*,*) 'NaN at VS C[Su, VG, VG, conj[Su]]' 
End Do
results2(17)=temp
result = sum(results1)+sum(results2) 
effPot2L = result * oo16pi2 * oo16Pi2
End Subroutine CalculateEffPot2Loop


Subroutine SecondDerivativeEffPot2Loop(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,kont,pi2L,ti2L)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Integer, Intent(inout):: kont
Real(dp), Intent(out) :: pi2L(8,8),ti2L(8)
Integer :: i,i1,i2,i3,NrContr 
Integer :: iv1, iv2 
Integer :: NrContr1,NrContr2 !nr of contributing diagrams
Real(dp) :: Q2,colorfactor,coeff,coeffbar
Complex(dp) :: result,result_ti,temp,temp_ti,temp_tj,tempbar,tempbar_ti,tempbar_tj
Complex(dp) :: coup1,coup2,coup1L,coup1R,coup2L,coup2R,coupx,coupxbar
Complex(dp) :: Di_coup1,Di_coup2,Di_coup1L,Di_coup1R,Di_coup2L,Di_coup2R,Di_coupx,Di_coupxbar
Complex(dp) :: Dj_coup1,Dj_coup2,Dj_coup1L,Dj_coup1R,Dj_coup2L,Dj_coup2R,Dj_coupx,Dj_coupxbar
Complex(dp) :: DDcoup1,DDcoup2,DDcoup1L,DDcoup1R,DDcoup2L,DDcoup2R,DDcoupx,DDcoupxbar
Complex(dp) :: results1(33),results2(17)
Complex(dp) :: results1_ti(33),results2_ti(17)
Real(dp) :: gout(112198) 
Complex(dp) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

Complex(dp) :: dcplAhAhAh(8,8,8,8),dcplAhAhhh(8,8,8,8),dcplAhhhhh(8,8,8,8),dcplAhHpmcHpm(8,8,8,8),   & 
& dcplAhSdcSd(8,6,6,8),dcplAhSucSu(8,6,6,8),dcplhhhhhh(8,8,8,8),dcplhhHpmcHpm(8,8,8,8),  & 
& dcplhhSdcSd(8,6,6,8),dcplhhSucSu(8,6,6,8),dcplHpmSucSd(8,6,6,8),dcplSdcHpmcSu(6,8,6,8),& 
& dcplSdcSdVG(6,6,8),dcplSucSuVG(6,6,8),dcplVGVGVG(8),dcplcChaChaAhL(5,5,8,8),           & 
& dcplcChaChaAhR(5,5,8,8),dcplChiChiAhL(10,10,8,8),dcplChiChiAhR(10,10,8,8),             & 
& dcplcFdFdAhL(3,3,8,8),dcplcFdFdAhR(3,3,8,8),dcplcFuFuAhL(3,3,8,8),dcplcFuFuAhR(3,3,8,8),& 
& dcplChiChacHpmL(10,5,8,8),dcplChiChacHpmR(10,5,8,8),dcplChaFucSdL(5,3,6,8),            & 
& dcplChaFucSdR(5,3,6,8),dcplcChaChahhL(5,5,8,8),dcplcChaChahhR(5,5,8,8),dcplcFdChaSuL(3,5,6,8),& 
& dcplcFdChaSuR(3,5,6,8),dcplChiChihhL(10,10,8,8),dcplChiChihhR(10,10,8,8),              & 
& dcplChiFdcSdL(10,3,6,8),dcplChiFdcSdR(10,3,6,8),dcplChiFucSuL(10,3,6,8),               & 
& dcplChiFucSuR(10,3,6,8),dcplcChaChiHpmL(5,10,8,8),dcplcChaChiHpmR(5,10,8,8),           & 
& dcplcFdChiSdL(3,10,6,8),dcplcFdChiSdR(3,10,6,8),dcplcFuChiSuL(3,10,6,8),               & 
& dcplcFuChiSuR(3,10,6,8),dcplGluFdcSdL(3,6,8),dcplGluFdcSdR(3,6,8),dcplcFdFdhhL(3,3,8,8),& 
& dcplcFdFdhhR(3,3,8,8),dcplcChaFdcSuL(5,3,6,8),dcplcChaFdcSuR(5,3,6,8),dcplcFuFdcHpmL(3,3,8,8),& 
& dcplcFuFdcHpmR(3,3,8,8),dcplGluFucSuL(3,6,8),dcplGluFucSuR(3,6,8),dcplcFuFuhhL(3,3,8,8),& 
& dcplcFuFuhhR(3,3,8,8),dcplcFdFuHpmL(3,3,8,8),dcplcFdFuHpmR(3,3,8,8),dcplcFdGluSdL(3,6,8),& 
& dcplcFdGluSdR(3,6,8),dcplcFuGluSuL(3,6,8),dcplcFuGluSuR(3,6,8),dcplcChacFuSdL(5,3,6,8),& 
& dcplcChacFuSdR(5,3,6,8),dcplcFdFdVGL(3,3,8),dcplcFdFdVGR(3,3,8),dcplcFuFuVGL(3,3,8),   & 
& dcplcFuFuVGR(3,3,8),dcplGluGluVGL(8),dcplGluGluVGR(8)

Complex(dp) :: dcplAhAhAhAh(8,8,8,8,8),dcplAhAhhhhh(8,8,8,8,8),dcplAhAhHpmcHpm(8,8,8,8,8),           & 
& dcplAhAhSdcSd(8,8,6,6,8),dcplAhAhSucSu(8,8,6,6,8),dcplhhhhhhhh(8,8,8,8,8),             & 
& dcplhhhhHpmcHpm(8,8,8,8,8),dcplhhhhSdcSd(8,8,6,6,8),dcplhhhhSucSu(8,8,6,6,8),          & 
& dcplHpmHpmcHpmcHpm(8,8,8,8,8),dcplHpmSdcHpmcSd(8,6,8,6,8),dcplHpmSucHpmcSu(8,6,8,6,8), & 
& dcplSdSdcSdcSd(6,6,6,6,8),dcplSdSucSdcSu(6,6,6,6,8),dcplSuSucSucSu(6,6,6,6,8),         & 
& dcplSdcSdVGVG(6,6,8),dcplSucSuVGVG(6,6,8)

Complex(dp) :: ddcplAhAhAh(8,8,8,8,8),ddcplAhAhhh(8,8,8,8,8),ddcplAhhhhh(8,8,8,8,8),ddcplAhHpmcHpm(8,8,8,8,8),& 
& ddcplAhSdcSd(8,6,6,8,8),ddcplAhSucSu(8,6,6,8,8),ddcplhhhhhh(8,8,8,8,8),ddcplhhHpmcHpm(8,8,8,8,8),& 
& ddcplhhSdcSd(8,6,6,8,8),ddcplhhSucSu(8,6,6,8,8),ddcplHpmSucSd(8,6,6,8,8),              & 
& ddcplSdcHpmcSu(6,8,6,8,8),ddcplSdcSdVG(6,6,8,8),ddcplSucSuVG(6,6,8,8),ddcplVGVGVG(8,8),& 
& ddcplcChaChaAhL(5,5,8,8,8),ddcplcChaChaAhR(5,5,8,8,8),ddcplChiChiAhL(10,10,8,8,8),     & 
& ddcplChiChiAhR(10,10,8,8,8),ddcplcFdFdAhL(3,3,8,8,8),ddcplcFdFdAhR(3,3,8,8,8),         & 
& ddcplcFuFuAhL(3,3,8,8,8),ddcplcFuFuAhR(3,3,8,8,8),ddcplChiChacHpmL(10,5,8,8,8),        & 
& ddcplChiChacHpmR(10,5,8,8,8),ddcplChaFucSdL(5,3,6,8,8),ddcplChaFucSdR(5,3,6,8,8),      & 
& ddcplcChaChahhL(5,5,8,8,8),ddcplcChaChahhR(5,5,8,8,8),ddcplcFdChaSuL(3,5,6,8,8),       & 
& ddcplcFdChaSuR(3,5,6,8,8),ddcplChiChihhL(10,10,8,8,8),ddcplChiChihhR(10,10,8,8,8),     & 
& ddcplChiFdcSdL(10,3,6,8,8),ddcplChiFdcSdR(10,3,6,8,8),ddcplChiFucSuL(10,3,6,8,8),      & 
& ddcplChiFucSuR(10,3,6,8,8),ddcplcChaChiHpmL(5,10,8,8,8),ddcplcChaChiHpmR(5,10,8,8,8),  & 
& ddcplcFdChiSdL(3,10,6,8,8),ddcplcFdChiSdR(3,10,6,8,8),ddcplcFuChiSuL(3,10,6,8,8),      & 
& ddcplcFuChiSuR(3,10,6,8,8),ddcplGluFdcSdL(3,6,8,8),ddcplGluFdcSdR(3,6,8,8),            & 
& ddcplcFdFdhhL(3,3,8,8,8),ddcplcFdFdhhR(3,3,8,8,8),ddcplcChaFdcSuL(5,3,6,8,8),          & 
& ddcplcChaFdcSuR(5,3,6,8,8),ddcplcFuFdcHpmL(3,3,8,8,8),ddcplcFuFdcHpmR(3,3,8,8,8),      & 
& ddcplGluFucSuL(3,6,8,8),ddcplGluFucSuR(3,6,8,8),ddcplcFuFuhhL(3,3,8,8,8),              & 
& ddcplcFuFuhhR(3,3,8,8,8),ddcplcFdFuHpmL(3,3,8,8,8),ddcplcFdFuHpmR(3,3,8,8,8),          & 
& ddcplcFdGluSdL(3,6,8,8),ddcplcFdGluSdR(3,6,8,8),ddcplcFuGluSuL(3,6,8,8),               & 
& ddcplcFuGluSuR(3,6,8,8),ddcplcChacFuSdL(5,3,6,8,8),ddcplcChacFuSdR(5,3,6,8,8),         & 
& ddcplcFdFdVGL(3,3,8,8),ddcplcFdFdVGR(3,3,8,8),ddcplcFuFuVGL(3,3,8,8),ddcplcFuFuVGR(3,3,8,8),& 
& ddcplGluGluVGL(8,8),ddcplGluGluVGR(8,8)

Complex(dp) :: ddcplAhAhAhAh(8,8,8,8,8,8),ddcplAhAhhhhh(8,8,8,8,8,8),ddcplAhAhHpmcHpm(8,8,8,8,8,8),  & 
& ddcplAhAhSdcSd(8,8,6,6,8,8),ddcplAhAhSucSu(8,8,6,6,8,8),ddcplhhhhhhhh(8,8,8,8,8,8),    & 
& ddcplhhhhHpmcHpm(8,8,8,8,8,8),ddcplhhhhSdcSd(8,8,6,6,8,8),ddcplhhhhSucSu(8,8,6,6,8,8), & 
& ddcplHpmHpmcHpmcHpm(8,8,8,8,8,8),ddcplHpmSdcHpmcSd(8,6,8,6,8,8),ddcplHpmSucHpmcSu(8,6,8,6,8,8),& 
& ddcplSdSdcSdcSd(6,6,6,6,8,8),ddcplSdSucSdcSu(6,6,6,6,8,8),ddcplSuSucSucSu(6,6,6,6,8,8),& 
& ddcplSdcSdVGVG(6,6,8,8),ddcplSucSuVGVG(6,6,8,8)

Real(dp) :: MSd(6),MSd2(6),MSu(6),MSu2(6),Mhh(8),Mhh2(8),MAh(8),MAh2(8),MHpm(8),MHpm2(8),         & 
& MChi(10),MChi2(10),MCha(5),MCha2(5),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,          & 
& MVZ,MVZ2,MVWm,MVWm2

Complex(dp) :: dMSd(6,8),dMSd2(6,8),dMSu(6,8),dMSu2(6,8),dMhh(8,8),dMhh2(8,8),dMAh(8,8),             & 
& dMAh2(8,8),dMHpm(8,8),dMHpm2(8,8),dMChi(10,8),dMChi2(10,8),dMCha(5,8),dMCha2(5,8),     & 
& dMFd(3,8),dMFd2(3,8),dMFu(3,8),dMFu2(3,8),dMGlu(1,8),dMGlu2(1,8),dMVZ(1,8),            & 
& dMVZ2(1,8),dMVWm(1,8),dMVWm2(1,8)

Complex(dp) :: ddMSd(6,8,8),ddMSd2(6,8,8),ddMSu(6,8,8),ddMSu2(6,8,8),ddMhh(8,8,8),ddMhh2(8,8,8),     & 
& ddMAh(8,8,8),ddMAh2(8,8,8),ddMHpm(8,8,8),ddMHpm2(8,8,8),ddMChi(10,8,8),ddMChi2(10,8,8),& 
& ddMCha(5,8,8),ddMCha2(5,8,8),ddMFd(3,8,8),ddMFd2(3,8,8),ddMFu(3,8,8),ddMFu2(3,8,8),    & 
& ddMGlu(1,8,8),ddMGlu2(1,8,8),ddMVZ(1,8,8),ddMVZ2(1,8,8),ddMVWm(1,8,8),ddMVWm2(1,8,8)

!! ------------------------------------------------- 
!! Calculate masses, couplings and their derivatives 
!! ------------------------------------------------- 

Do i1=1,8
Call FirstDerivativeMassesCoups(i1,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,              & 
& kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,gout)

Call GToMassesCoups(gout,MSd,MSd2,MSu,MSu2,Mhh,Mhh2,MAh,MAh2,MHpm,MHpm2,              & 
& MChi,MChi2,MCha,MCha2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MVZ,MVZ2,MVWm,MVWm2,cplAhAhAh,      & 
& cplAhAhhh,cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,         & 
& cplhhSdcSd,cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,        & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,       & 
& cplcChaChahhL,cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,       & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,     & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplcFdFdhhL,cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplGluFucSuL,cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,       & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR,             & 
& cplAhAhAhAh,cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,          & 
& cplhhhhHpmcHpm,cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,            & 
& cplHpmSucHpmcSu,cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)

dMSd(:,i1) = MSd
dMSd2(:,i1) = MSd2
dMSu(:,i1) = MSu
dMSu2(:,i1) = MSu2
dMhh(:,i1) = Mhh
dMhh2(:,i1) = Mhh2
dMAh(:,i1) = MAh
dMAh2(:,i1) = MAh2
dMHpm(:,i1) = MHpm
dMHpm2(:,i1) = MHpm2
dMChi(1:,i1) = MChi
dMChi2(1:,i1) = MChi2
dMCha(:,i1) = MCha
dMCha2(:,i1) = MCha2
dMFd(:,i1) = MFd
dMFd2(:,i1) = MFd2
dMFu(:,i1) = MFu
dMFu2(:,i1) = MFu2
dMGlu(:,i1) = MGlu
dMGlu2(:,i1) = MGlu2
dMVZ(:,i1) = MVZ
dMVZ2(:,i1) = MVZ2
dMVWm(:,i1) = MVWm
dMVWm2(:,i1) = MVWm2
dcplAhAhAh(:,:,:,i1) = cplAhAhAh
dcplAhAhhh(:,:,:,i1) = cplAhAhhh
dcplAhhhhh(:,:,:,i1) = cplAhhhhh
dcplAhHpmcHpm(:,:,:,i1) = cplAhHpmcHpm
dcplAhSdcSd(:,:,:,i1) = cplAhSdcSd
dcplAhSucSu(:,:,:,i1) = cplAhSucSu
dcplhhhhhh(:,:,:,i1) = cplhhhhhh
dcplhhHpmcHpm(:,:,:,i1) = cplhhHpmcHpm
dcplhhSdcSd(:,:,:,i1) = cplhhSdcSd
dcplhhSucSu(:,:,:,i1) = cplhhSucSu
dcplHpmSucSd(:,:,:,i1) = cplHpmSucSd
dcplSdcHpmcSu(:,:,:,i1) = cplSdcHpmcSu
dcplSdcSdVG(:,:,i1) = cplSdcSdVG
dcplSucSuVG(:,:,i1) = cplSucSuVG
dcplVGVGVG(i1) = cplVGVGVG
dcplcChaChaAhL(:,:,:,i1) = cplcChaChaAhL
dcplcChaChaAhR(:,:,:,i1) = cplcChaChaAhR
dcplChiChiAhL(1:,1:,:,i1) = cplChiChiAhL
dcplChiChiAhR(1:,1:,:,i1) = cplChiChiAhR
dcplcFdFdAhL(:,:,:,i1) = cplcFdFdAhL
dcplcFdFdAhR(:,:,:,i1) = cplcFdFdAhR
dcplcFuFuAhL(:,:,:,i1) = cplcFuFuAhL
dcplcFuFuAhR(:,:,:,i1) = cplcFuFuAhR
dcplChiChacHpmL(1:,:,:,i1) = cplChiChacHpmL
dcplChiChacHpmR(1:,:,:,i1) = cplChiChacHpmR
dcplChaFucSdL(:,:,:,i1) = cplChaFucSdL
dcplChaFucSdR(:,:,:,i1) = cplChaFucSdR
dcplcChaChahhL(:,:,:,i1) = cplcChaChahhL
dcplcChaChahhR(:,:,:,i1) = cplcChaChahhR
dcplcFdChaSuL(:,:,:,i1) = cplcFdChaSuL
dcplcFdChaSuR(:,:,:,i1) = cplcFdChaSuR
dcplChiChihhL(1:,1:,:,i1) = cplChiChihhL
dcplChiChihhR(1:,1:,:,i1) = cplChiChihhR
dcplChiFdcSdL(1:,:,:,i1) = cplChiFdcSdL
dcplChiFdcSdR(1:,:,:,i1) = cplChiFdcSdR
dcplChiFucSuL(1:,:,:,i1) = cplChiFucSuL
dcplChiFucSuR(1:,:,:,i1) = cplChiFucSuR
dcplcChaChiHpmL(:,1:,:,i1) = cplcChaChiHpmL
dcplcChaChiHpmR(:,1:,:,i1) = cplcChaChiHpmR
dcplcFdChiSdL(:,1:,:,i1) = cplcFdChiSdL
dcplcFdChiSdR(:,1:,:,i1) = cplcFdChiSdR
dcplcFuChiSuL(:,1:,:,i1) = cplcFuChiSuL
dcplcFuChiSuR(:,1:,:,i1) = cplcFuChiSuR
dcplGluFdcSdL(:,:,i1) = cplGluFdcSdL
dcplGluFdcSdR(:,:,i1) = cplGluFdcSdR
dcplcFdFdhhL(:,:,:,i1) = cplcFdFdhhL
dcplcFdFdhhR(:,:,:,i1) = cplcFdFdhhR
dcplcChaFdcSuL(:,:,:,i1) = cplcChaFdcSuL
dcplcChaFdcSuR(:,:,:,i1) = cplcChaFdcSuR
dcplcFuFdcHpmL(:,:,:,i1) = cplcFuFdcHpmL
dcplcFuFdcHpmR(:,:,:,i1) = cplcFuFdcHpmR
dcplGluFucSuL(:,:,i1) = cplGluFucSuL
dcplGluFucSuR(:,:,i1) = cplGluFucSuR
dcplcFuFuhhL(:,:,:,i1) = cplcFuFuhhL
dcplcFuFuhhR(:,:,:,i1) = cplcFuFuhhR
dcplcFdFuHpmL(:,:,:,i1) = cplcFdFuHpmL
dcplcFdFuHpmR(:,:,:,i1) = cplcFdFuHpmR
dcplcFdGluSdL(:,:,i1) = cplcFdGluSdL
dcplcFdGluSdR(:,:,i1) = cplcFdGluSdR
dcplcFuGluSuL(:,:,i1) = cplcFuGluSuL
dcplcFuGluSuR(:,:,i1) = cplcFuGluSuR
dcplcChacFuSdL(:,:,:,i1) = cplcChacFuSdL
dcplcChacFuSdR(:,:,:,i1) = cplcChacFuSdR
dcplcFdFdVGL(:,:,i1) = cplcFdFdVGL
dcplcFdFdVGR(:,:,i1) = cplcFdFdVGR
dcplcFuFuVGL(:,:,i1) = cplcFuFuVGL
dcplcFuFuVGR(:,:,i1) = cplcFuFuVGR
dcplGluGluVGL(i1) = cplGluGluVGL
dcplGluGluVGR(i1) = cplGluGluVGR
dcplAhAhAhAh(:,:,:,:,i1) = cplAhAhAhAh
dcplAhAhhhhh(:,:,:,:,i1) = cplAhAhhhhh
dcplAhAhHpmcHpm(:,:,:,:,i1) = cplAhAhHpmcHpm
dcplAhAhSdcSd(:,:,:,:,i1) = cplAhAhSdcSd
dcplAhAhSucSu(:,:,:,:,i1) = cplAhAhSucSu
dcplhhhhhhhh(:,:,:,:,i1) = cplhhhhhhhh
dcplhhhhHpmcHpm(:,:,:,:,i1) = cplhhhhHpmcHpm
dcplhhhhSdcSd(:,:,:,:,i1) = cplhhhhSdcSd
dcplhhhhSucSu(:,:,:,:,i1) = cplhhhhSucSu
dcplHpmHpmcHpmcHpm(:,:,:,:,i1) = cplHpmHpmcHpmcHpm
dcplHpmSdcHpmcSd(:,:,:,:,i1) = cplHpmSdcHpmcSd
dcplHpmSucHpmcSu(:,:,:,:,i1) = cplHpmSucHpmcSu
dcplSdSdcSdcSd(:,:,:,:,i1) = cplSdSdcSdcSd
dcplSdSucSdcSu(:,:,:,:,i1) = cplSdSucSdcSu
dcplSuSucSucSu(:,:,:,:,i1) = cplSuSucSucSu
dcplSdcSdVGVG(:,:,i1) = cplSdcSdVGVG
dcplSucSuVGVG(:,:,i1) = cplSucSuVGVG
End Do 
 
Do i1=1,8
  Do i2=i1,8
Call SecondDerivativeMassesCoups(i1,i2,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,gout)

Call GToMassesCoups(gout,MSd,MSd2,MSu,MSu2,Mhh,Mhh2,MAh,MAh2,MHpm,MHpm2,              & 
& MChi,MChi2,MCha,MCha2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MVZ,MVZ2,MVWm,MVWm2,cplAhAhAh,      & 
& cplAhAhhh,cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,         & 
& cplhhSdcSd,cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,        & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,       & 
& cplcChaChahhL,cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,       & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,     & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplcFdFdhhL,cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplGluFucSuL,cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,       & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR,             & 
& cplAhAhAhAh,cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,          & 
& cplhhhhHpmcHpm,cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,            & 
& cplHpmSucHpmcSu,cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)

ddMSd(:,i1,i2) = MSd
ddMSd2(:,i1,i2) = MSd2
ddMSu(:,i1,i2) = MSu
ddMSu2(:,i1,i2) = MSu2
ddMhh(:,i1,i2) = Mhh
ddMhh2(:,i1,i2) = Mhh2
ddMAh(:,i1,i2) = MAh
ddMAh2(:,i1,i2) = MAh2
ddMHpm(:,i1,i2) = MHpm
ddMHpm2(:,i1,i2) = MHpm2
ddMChi(1:,i1,i2) = MChi
ddMChi2(1:,i1,i2) = MChi2
ddMCha(:,i1,i2) = MCha
ddMCha2(:,i1,i2) = MCha2
ddMFd(:,i1,i2) = MFd
ddMFd2(:,i1,i2) = MFd2
ddMFu(:,i1,i2) = MFu
ddMFu2(:,i1,i2) = MFu2
ddMGlu(:,i1,i2) = MGlu
ddMGlu2(:,i1,i2) = MGlu2
ddMVZ(:,i1,i2) = MVZ
ddMVZ2(:,i1,i2) = MVZ2
ddMVWm(:,i1,i2) = MVWm
ddMVWm2(:,i1,i2) = MVWm2
ddcplAhAhAh(:,:,:,i1,i2) = cplAhAhAh
ddcplAhAhhh(:,:,:,i1,i2) = cplAhAhhh
ddcplAhhhhh(:,:,:,i1,i2) = cplAhhhhh
ddcplAhHpmcHpm(:,:,:,i1,i2) = cplAhHpmcHpm
ddcplAhSdcSd(:,:,:,i1,i2) = cplAhSdcSd
ddcplAhSucSu(:,:,:,i1,i2) = cplAhSucSu
ddcplhhhhhh(:,:,:,i1,i2) = cplhhhhhh
ddcplhhHpmcHpm(:,:,:,i1,i2) = cplhhHpmcHpm
ddcplhhSdcSd(:,:,:,i1,i2) = cplhhSdcSd
ddcplhhSucSu(:,:,:,i1,i2) = cplhhSucSu
ddcplHpmSucSd(:,:,:,i1,i2) = cplHpmSucSd
ddcplSdcHpmcSu(:,:,:,i1,i2) = cplSdcHpmcSu
ddcplSdcSdVG(:,:,i1,i2) = cplSdcSdVG
ddcplSucSuVG(:,:,i1,i2) = cplSucSuVG
ddcplVGVGVG(i1,i2) = cplVGVGVG
ddcplcChaChaAhL(:,:,:,i1,i2) = cplcChaChaAhL
ddcplcChaChaAhR(:,:,:,i1,i2) = cplcChaChaAhR
ddcplChiChiAhL(1:,1:,:,i1,i2) = cplChiChiAhL
ddcplChiChiAhR(1:,1:,:,i1,i2) = cplChiChiAhR
ddcplcFdFdAhL(:,:,:,i1,i2) = cplcFdFdAhL
ddcplcFdFdAhR(:,:,:,i1,i2) = cplcFdFdAhR
ddcplcFuFuAhL(:,:,:,i1,i2) = cplcFuFuAhL
ddcplcFuFuAhR(:,:,:,i1,i2) = cplcFuFuAhR
ddcplChiChacHpmL(1:,:,:,i1,i2) = cplChiChacHpmL
ddcplChiChacHpmR(1:,:,:,i1,i2) = cplChiChacHpmR
ddcplChaFucSdL(:,:,:,i1,i2) = cplChaFucSdL
ddcplChaFucSdR(:,:,:,i1,i2) = cplChaFucSdR
ddcplcChaChahhL(:,:,:,i1,i2) = cplcChaChahhL
ddcplcChaChahhR(:,:,:,i1,i2) = cplcChaChahhR
ddcplcFdChaSuL(:,:,:,i1,i2) = cplcFdChaSuL
ddcplcFdChaSuR(:,:,:,i1,i2) = cplcFdChaSuR
ddcplChiChihhL(1:,1:,:,i1,i2) = cplChiChihhL
ddcplChiChihhR(1:,1:,:,i1,i2) = cplChiChihhR
ddcplChiFdcSdL(1:,:,:,i1,i2) = cplChiFdcSdL
ddcplChiFdcSdR(1:,:,:,i1,i2) = cplChiFdcSdR
ddcplChiFucSuL(1:,:,:,i1,i2) = cplChiFucSuL
ddcplChiFucSuR(1:,:,:,i1,i2) = cplChiFucSuR
ddcplcChaChiHpmL(:,1:,:,i1,i2) = cplcChaChiHpmL
ddcplcChaChiHpmR(:,1:,:,i1,i2) = cplcChaChiHpmR
ddcplcFdChiSdL(:,1:,:,i1,i2) = cplcFdChiSdL
ddcplcFdChiSdR(:,1:,:,i1,i2) = cplcFdChiSdR
ddcplcFuChiSuL(:,1:,:,i1,i2) = cplcFuChiSuL
ddcplcFuChiSuR(:,1:,:,i1,i2) = cplcFuChiSuR
ddcplGluFdcSdL(:,:,i1,i2) = cplGluFdcSdL
ddcplGluFdcSdR(:,:,i1,i2) = cplGluFdcSdR
ddcplcFdFdhhL(:,:,:,i1,i2) = cplcFdFdhhL
ddcplcFdFdhhR(:,:,:,i1,i2) = cplcFdFdhhR
ddcplcChaFdcSuL(:,:,:,i1,i2) = cplcChaFdcSuL
ddcplcChaFdcSuR(:,:,:,i1,i2) = cplcChaFdcSuR
ddcplcFuFdcHpmL(:,:,:,i1,i2) = cplcFuFdcHpmL
ddcplcFuFdcHpmR(:,:,:,i1,i2) = cplcFuFdcHpmR
ddcplGluFucSuL(:,:,i1,i2) = cplGluFucSuL
ddcplGluFucSuR(:,:,i1,i2) = cplGluFucSuR
ddcplcFuFuhhL(:,:,:,i1,i2) = cplcFuFuhhL
ddcplcFuFuhhR(:,:,:,i1,i2) = cplcFuFuhhR
ddcplcFdFuHpmL(:,:,:,i1,i2) = cplcFdFuHpmL
ddcplcFdFuHpmR(:,:,:,i1,i2) = cplcFdFuHpmR
ddcplcFdGluSdL(:,:,i1,i2) = cplcFdGluSdL
ddcplcFdGluSdR(:,:,i1,i2) = cplcFdGluSdR
ddcplcFuGluSuL(:,:,i1,i2) = cplcFuGluSuL
ddcplcFuGluSuR(:,:,i1,i2) = cplcFuGluSuR
ddcplcChacFuSdL(:,:,:,i1,i2) = cplcChacFuSdL
ddcplcChacFuSdR(:,:,:,i1,i2) = cplcChacFuSdR
ddcplcFdFdVGL(:,:,i1,i2) = cplcFdFdVGL
ddcplcFdFdVGR(:,:,i1,i2) = cplcFdFdVGR
ddcplcFuFuVGL(:,:,i1,i2) = cplcFuFuVGL
ddcplcFuFuVGR(:,:,i1,i2) = cplcFuFuVGR
ddcplGluGluVGL(i1,i2) = cplGluGluVGL
ddcplGluGluVGR(i1,i2) = cplGluGluVGR
ddcplAhAhAhAh(:,:,:,:,i1,i2) = cplAhAhAhAh
ddcplAhAhhhhh(:,:,:,:,i1,i2) = cplAhAhhhhh
ddcplAhAhHpmcHpm(:,:,:,:,i1,i2) = cplAhAhHpmcHpm
ddcplAhAhSdcSd(:,:,:,:,i1,i2) = cplAhAhSdcSd
ddcplAhAhSucSu(:,:,:,:,i1,i2) = cplAhAhSucSu
ddcplhhhhhhhh(:,:,:,:,i1,i2) = cplhhhhhhhh
ddcplhhhhHpmcHpm(:,:,:,:,i1,i2) = cplhhhhHpmcHpm
ddcplhhhhSdcSd(:,:,:,:,i1,i2) = cplhhhhSdcSd
ddcplhhhhSucSu(:,:,:,:,i1,i2) = cplhhhhSucSu
ddcplHpmHpmcHpmcHpm(:,:,:,:,i1,i2) = cplHpmHpmcHpmcHpm
ddcplHpmSdcHpmcSd(:,:,:,:,i1,i2) = cplHpmSdcHpmcSd
ddcplHpmSucHpmcSu(:,:,:,:,i1,i2) = cplHpmSucHpmcSu
ddcplSdSdcSdcSd(:,:,:,:,i1,i2) = cplSdSdcSdcSd
ddcplSdSucSdcSu(:,:,:,:,i1,i2) = cplSdSucSdcSu
ddcplSuSucSucSu(:,:,:,:,i1,i2) = cplSuSucSucSu
ddcplSdcSdVGVG(:,:,i1,i2) = cplSdcSdVGVG
ddcplSucSuVGVG(:,:,i1,i2) = cplSucSuVGVG
  End Do 
 End Do 
 
Call TreeMassesEffPot(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,               & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,.True.,kont)

Call CouplingsForEffPot3(lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZA,ZH,Ye,Te,               & 
& ZP,Yd,Td,ZD,Yu,Tu,ZU,g3,ZER,ZEL,UV,ZDL,ZDR,ZUL,ZUR,pG,cplAhAhAh,cplAhAhhh,             & 
& cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,cplhhSdcSd,        & 
& cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplcChaChaAhL,     & 
& cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,cplcFuFuAhL,           & 
& cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,cplcChaChahhL,     & 
& cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,cplChiFdcSdL,        & 
& cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,cplcFdChiSdL,     & 
& cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,cplcFdFdhhL,          & 
& cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplGluFucSuL,      & 
& cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdGluSdL,           & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,cplcFdFdVGL,        & 
& cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR)

Call CouplingsForEffPot4(lam,Yv,kap,ZA,ZH,Ye,ZP,Yd,ZD,Yu,ZU,g3,cplAhAhAhAh,           & 
& cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,cplhhhhHpmcHpm,       & 
& cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,cplHpmSucHpmcSu,           & 
& cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)



!! ------------------------------------------------- 
!! Calculate derivative of effective potential      
!! ------------------------------------------------- 



Q2 = getRenormalizationScale()
Do iv1=1,8
  Do iv2=iv1,8
    result = ZeroC
    result_ti = ZeroC
    results1 = ZeroC
    results1_ti = ZeroC
    results2 = ZeroC
    results2_ti = ZeroC


! ----- Topology1 (sunrise): diagrams w. 3 Particles and 2 Vertices


! ----- diagrams of type SSS, 11 ------ 

! ---- Ah,Ah,Ah ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhAhAh(i1,i2,i3)
coup2 = cplAhAhAh(i1,i2,i3)
Di_coup1 = dcplAhAhAh(i1,i2,i3,iv1)
Dj_coup1 = dcplAhAhAh(i1,i2,i3,iv2)
DDcoup1 = ddcplAhAhAh(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MAh2(i1),MAh2(i2),MAh2(i3),dMAh2(i1,iv1)            & 
& ,dMAh2(i2,iv1),dMAh2(i3,iv1),dMAh2(i1,iv2),dMAh2(i2,iv2),dMAh2(i3,iv2),ddMAh2(i1,iv1,iv2)& 
& ,ddMAh2(i2,iv1,iv2),ddMAh2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 1._dp/12._dp
colorfactor = 1
results1(1)=results1(1) + coeff*colorfactor*temp
results1_ti(1)=results1_ti(1) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(1).eq.results1(1)))  write(*,*) 'NaN at SSS C[Ah, Ah, Ah]' 
! ---- Ah,Ah,hh ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhAhhh(i1,i2,i3)
coup2 = cplAhAhhh(i1,i2,i3)
Di_coup1 = dcplAhAhhh(i1,i2,i3,iv1)
Dj_coup1 = dcplAhAhhh(i1,i2,i3,iv2)
DDcoup1 = ddcplAhAhhh(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MAh2(i1),MAh2(i2),Mhh2(i3),dMAh2(i1,iv1)            & 
& ,dMAh2(i2,iv1),dMhh2(i3,iv1),dMAh2(i1,iv2),dMAh2(i2,iv2),dMhh2(i3,iv2),ddMAh2(i1,iv1,iv2)& 
& ,ddMAh2(i2,iv1,iv2),ddMhh2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.25_dp
colorfactor = 1
results1(2)=results1(2) + coeff*colorfactor*temp
results1_ti(2)=results1_ti(2) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(2).eq.results1(2)))  write(*,*) 'NaN at SSS C[Ah, Ah, hh]' 
! ---- Ah,hh,hh ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhhhhh(i1,i2,i3)
coup2 = cplAhhhhh(i1,i2,i3)
Di_coup1 = dcplAhhhhh(i1,i2,i3,iv1)
Dj_coup1 = dcplAhhhhh(i1,i2,i3,iv2)
DDcoup1 = ddcplAhhhhh(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MAh2(i1),Mhh2(i2),Mhh2(i3),dMAh2(i1,iv1)            & 
& ,dMhh2(i2,iv1),dMhh2(i3,iv1),dMAh2(i1,iv2),dMhh2(i2,iv2),dMhh2(i3,iv2),ddMAh2(i1,iv1,iv2)& 
& ,ddMhh2(i2,iv1,iv2),ddMhh2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.25_dp
colorfactor = 1
results1(3)=results1(3) + coeff*colorfactor*temp
results1_ti(3)=results1_ti(3) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(3).eq.results1(3)))  write(*,*) 'NaN at SSS C[Ah, hh, hh]' 
! ---- Ah,Hpm,conj[Hpm] ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhHpmcHpm(i1,i2,i3)
coup2 = cplAhHpmcHpm(i1,i3,i2)
Di_coup1 = dcplAhHpmcHpm(i1,i2,i3,iv1)
Dj_coup1 = dcplAhHpmcHpm(i1,i2,i3,iv2)
DDcoup1 = ddcplAhHpmcHpm(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MAh2(i1),MHpm2(i2),MHpm2(i3),dMAh2(i1,iv1)          & 
& ,dMHpm2(i2,iv1),dMHpm2(i3,iv1),dMAh2(i1,iv2),dMHpm2(i2,iv2),dMHpm2(i3,iv2)             & 
& ,ddMAh2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),ddMHpm2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 1
results1(4)=results1(4) + coeff*colorfactor*temp
results1_ti(4)=results1_ti(4) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(4).eq.results1(4)))  write(*,*) 'NaN at SSS C[Ah, Hpm, conj[Hpm]]' 
! ---- Ah,Sd,conj[Sd] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplAhSdcSd(i1,i2,i3)
coup2 = cplAhSdcSd(i1,i3,i2)
Di_coup1 = dcplAhSdcSd(i1,i2,i3,iv1)
Dj_coup1 = dcplAhSdcSd(i1,i2,i3,iv2)
DDcoup1 = ddcplAhSdcSd(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MAh2(i1),MSd2(i2),MSd2(i3),dMAh2(i1,iv1)            & 
& ,dMSd2(i2,iv1),dMSd2(i3,iv1),dMAh2(i1,iv2),dMSd2(i2,iv2),dMSd2(i3,iv2),ddMAh2(i1,iv1,iv2)& 
& ,ddMSd2(i2,iv1,iv2),ddMSd2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 3
results1(5)=results1(5) + coeff*colorfactor*temp
results1_ti(5)=results1_ti(5) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(5).eq.results1(5)))  write(*,*) 'NaN at SSS C[Ah, Sd, conj[Sd]]' 
! ---- Ah,Su,conj[Su] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplAhSucSu(i1,i2,i3)
coup2 = cplAhSucSu(i1,i3,i2)
Di_coup1 = dcplAhSucSu(i1,i2,i3,iv1)
Dj_coup1 = dcplAhSucSu(i1,i2,i3,iv2)
DDcoup1 = ddcplAhSucSu(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MAh2(i1),MSu2(i2),MSu2(i3),dMAh2(i1,iv1)            & 
& ,dMSu2(i2,iv1),dMSu2(i3,iv1),dMAh2(i1,iv2),dMSu2(i2,iv2),dMSu2(i3,iv2),ddMAh2(i1,iv1,iv2)& 
& ,ddMSu2(i2,iv1,iv2),ddMSu2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 3
results1(6)=results1(6) + coeff*colorfactor*temp
results1_ti(6)=results1_ti(6) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(6).eq.results1(6)))  write(*,*) 'NaN at SSS C[Ah, Su, conj[Su]]' 
! ---- hh,hh,hh ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplhhhhhh(i1,i2,i3)
coup2 = cplhhhhhh(i1,i2,i3)
Di_coup1 = dcplhhhhhh(i1,i2,i3,iv1)
Dj_coup1 = dcplhhhhhh(i1,i2,i3,iv2)
DDcoup1 = ddcplhhhhhh(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(Mhh2(i1),Mhh2(i2),Mhh2(i3),dMhh2(i1,iv1)            & 
& ,dMhh2(i2,iv1),dMhh2(i3,iv1),dMhh2(i1,iv2),dMhh2(i2,iv2),dMhh2(i3,iv2),ddMhh2(i1,iv1,iv2)& 
& ,ddMhh2(i2,iv1,iv2),ddMhh2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 1._dp/12._dp
colorfactor = 1
results1(7)=results1(7) + coeff*colorfactor*temp
results1_ti(7)=results1_ti(7) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(7).eq.results1(7)))  write(*,*) 'NaN at SSS C[hh, hh, hh]' 
! ---- hh,Hpm,conj[Hpm] ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplhhHpmcHpm(i1,i2,i3)
coup2 = cplhhHpmcHpm(i1,i3,i2)
Di_coup1 = dcplhhHpmcHpm(i1,i2,i3,iv1)
Dj_coup1 = dcplhhHpmcHpm(i1,i2,i3,iv2)
DDcoup1 = ddcplhhHpmcHpm(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(Mhh2(i1),MHpm2(i2),MHpm2(i3),dMhh2(i1,iv1)          & 
& ,dMHpm2(i2,iv1),dMHpm2(i3,iv1),dMhh2(i1,iv2),dMHpm2(i2,iv2),dMHpm2(i3,iv2)             & 
& ,ddMhh2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),ddMHpm2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 1
results1(8)=results1(8) + coeff*colorfactor*temp
results1_ti(8)=results1_ti(8) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(8).eq.results1(8)))  write(*,*) 'NaN at SSS C[hh, Hpm, conj[Hpm]]' 
! ---- hh,Sd,conj[Sd] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplhhSdcSd(i1,i2,i3)
coup2 = cplhhSdcSd(i1,i3,i2)
Di_coup1 = dcplhhSdcSd(i1,i2,i3,iv1)
Dj_coup1 = dcplhhSdcSd(i1,i2,i3,iv2)
DDcoup1 = ddcplhhSdcSd(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(Mhh2(i1),MSd2(i2),MSd2(i3),dMhh2(i1,iv1)            & 
& ,dMSd2(i2,iv1),dMSd2(i3,iv1),dMhh2(i1,iv2),dMSd2(i2,iv2),dMSd2(i3,iv2),ddMhh2(i1,iv1,iv2)& 
& ,ddMSd2(i2,iv1,iv2),ddMSd2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 3
results1(9)=results1(9) + coeff*colorfactor*temp
results1_ti(9)=results1_ti(9) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(9).eq.results1(9)))  write(*,*) 'NaN at SSS C[hh, Sd, conj[Sd]]' 
! ---- hh,Su,conj[Su] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplhhSucSu(i1,i2,i3)
coup2 = cplhhSucSu(i1,i3,i2)
Di_coup1 = dcplhhSucSu(i1,i2,i3,iv1)
Dj_coup1 = dcplhhSucSu(i1,i2,i3,iv2)
DDcoup1 = ddcplhhSucSu(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(Mhh2(i1),MSu2(i2),MSu2(i3),dMhh2(i1,iv1)            & 
& ,dMSu2(i2,iv1),dMSu2(i3,iv1),dMhh2(i1,iv2),dMSu2(i2,iv2),dMSu2(i3,iv2),ddMhh2(i1,iv1,iv2)& 
& ,ddMSu2(i2,iv1,iv2),ddMSu2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 3
results1(10)=results1(10) + coeff*colorfactor*temp
results1_ti(10)=results1_ti(10) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(10).eq.results1(10)))  write(*,*) 'NaN at SSS C[hh, Su, conj[Su]]' 
! ---- Sd,conj[Hpm],conj[Su] ----
Do i1=1,6
 Do i2=1,8
    Do i3=1,6
coup1 = cplSdcHpmcSu(i1,i2,i3)
coup2 = cplHpmSucSd(i2,i3,i1)
Di_coup1 = dcplSdcHpmcSu(i1,i2,i3,iv1)
Dj_coup1 = dcplSdcHpmcSu(i1,i2,i3,iv2)
DDcoup1 = ddcplSdcHpmcSu(i1,i2,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MSd2(i1),MHpm2(i2),MSu2(i3),dMSd2(i1,iv1)           & 
& ,dMHpm2(i2,iv1),dMSu2(i3,iv1),dMSd2(i1,iv2),dMHpm2(i2,iv2),dMSu2(i3,iv2)               & 
& ,ddMSd2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),ddMSu2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'SSS   ',Q2,temp,temp_ti,temp_tj)
coeff = 1._dp
colorfactor = 3
results1(11)=results1(11) + coeff*colorfactor*temp
results1_ti(11)=results1_ti(11) + coeff*colorfactor*temp_ti
    End Do
  End Do
End Do
if (.not.(results1(11).eq.results1(11)))  write(*,*) 'NaN at SSS C[Sd, conj[Hpm], conj[Su]]' 

! ----- diagrams of type FFS, 16 ------ 

! ---- Ah,Cha,bar[Cha] ----
Do i1=1,8
 Do i2=1,5
    Do i3=1,5
coup1L = cplcChaChaAhL(i3,i2,i1)
coup1R = cplcChaChaAhR(i3,i2,i1)
coup2L = cplcChaChaAhL(i2,i3,i1)
coup2R = cplcChaChaAhR(i2,i3,i1)
Di_coup1L = dcplcChaChaAhL(i3,i2,i1,iv1)
Di_coup1R = dcplcChaChaAhR(i3,i2,i1,iv1)
Dj_coup1L = dcplcChaChaAhL(i3,i2,i1,iv2)
Dj_coup1R = dcplcChaChaAhR(i3,i2,i1,iv2)
DDcoup1L = ddcplcChaChaAhL(i3,i2,i1,iv1,iv2)
DDcoup1R = ddcplcChaChaAhR(i3,i2,i1,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MCha2(i3),MCha2(i2),MAh2(i1),dMCha2(i3,iv1)         & 
& ,dMCha2(i2,iv1),dMAh2(i1,iv1),dMCha2(i3,iv2),dMCha2(i2,iv2),dMAh2(i1,iv2)              & 
& ,ddMCha2(i3,iv1,iv2),ddMCha2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MCha2(i3),MCha2(i2),MAh2(i1),dMCha2(i3,iv1)         & 
& ,dMCha2(i2,iv1),dMAh2(i1,iv1),dMCha2(i3,iv2),dMCha2(i2,iv2),dMAh2(i1,iv2)              & 
& ,ddMCha2(i3,iv1,iv2),ddMCha2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 1
results1(12)=results1(12) + coeff*colorfactor*temp
results1_ti(12)=results1_ti(12) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(12)=results1(12) + coeffbar*colorfactor*tempbar
results1_ti(12)=results1_ti(12) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(12).eq.results1(12)))  write(*,*) 'NaN at FFS C[Ah, Cha, bar[Cha]]' 
! ---- Ah,Chi,Chi ----
Do i1=1,8
 Do i2=1,10
    Do i3=1,10
coup1L = cplChiChiAhL(i2,i3,i1)
coup1R = cplChiChiAhR(i2,i3,i1)
coup2L = cplChiChiAhL(i2,i3,i1)
coup2R = cplChiChiAhR(i2,i3,i1)
Di_coup1L = dcplChiChiAhL(i2,i3,i1,iv1)
Di_coup1R = dcplChiChiAhR(i2,i3,i1,iv1)
Dj_coup1L = dcplChiChiAhL(i2,i3,i1,iv2)
Dj_coup1R = dcplChiChiAhR(i2,i3,i1,iv2)
DDcoup1L = ddcplChiChiAhL(i2,i3,i1,iv1,iv2)
DDcoup1R = ddcplChiChiAhR(i2,i3,i1,iv1,iv2)
coupx = abs(coup1L)**2
Di_coupx=Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) 
Dj_coupx=Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L)  
Call SecondDerivativeVeff_sunrise(MChi2(i3),MChi2(i2),MAh2(i1),dMChi2(i3,iv1)         & 
& ,dMChi2(i2,iv1),dMAh2(i1,iv1),dMChi2(i3,iv2),dMChi2(i2,iv2),dMAh2(i1,iv2)              & 
& ,ddMChi2(i3,iv1,iv2),ddMChi2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = Real(coup1L**2,dp) 
Di_coupxbar = Real(2*Di_coup1L*coup1L,dp) 
Dj_coupxbar = Real(2*Dj_coup1L*coup1L,dp) 
DDcoupxbar = Real(2*DDcoup1L*coup1L + 2*Di_coup1L*Dj_coup1L,dp) 
Call SecondDerivativeVeff_sunrise(MChi2(i3),MChi2(i2),MAh2(i1),dMChi2(i3,iv1)         & 
& ,dMChi2(i2,iv1),dMAh2(i1,iv1),dMChi2(i3,iv2),dMChi2(i2,iv2),dMAh2(i1,iv2)              & 
& ,ddMChi2(i3,iv1,iv2),ddMChi2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 1
results1(13)=results1(13) + coeff*colorfactor*temp
results1_ti(13)=results1_ti(13) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(13)=results1(13) + coeffbar*colorfactor*tempbar
results1_ti(13)=results1_ti(13) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(13).eq.results1(13)))  write(*,*) 'NaN at FFS C[Ah, Chi, Chi]' 
! ---- Ah,Fd,bar[Fd] ----
Do i1=1,8
 Do i2=1,3
    Do i3=1,3
coup1L = cplcFdFdAhL(i3,i2,i1)
coup1R = cplcFdFdAhR(i3,i2,i1)
coup2L = cplcFdFdAhL(i2,i3,i1)
coup2R = cplcFdFdAhR(i2,i3,i1)
Di_coup1L = dcplcFdFdAhL(i3,i2,i1,iv1)
Di_coup1R = dcplcFdFdAhR(i3,i2,i1,iv1)
Dj_coup1L = dcplcFdFdAhL(i3,i2,i1,iv2)
Dj_coup1R = dcplcFdFdAhR(i3,i2,i1,iv2)
DDcoup1L = ddcplcFdFdAhL(i3,i2,i1,iv1,iv2)
DDcoup1R = ddcplcFdFdAhR(i3,i2,i1,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFd2(i2),MAh2(i1),dMFd2(i3,iv1)            & 
& ,dMFd2(i2,iv1),dMAh2(i1,iv1),dMFd2(i3,iv2),dMFd2(i2,iv2),dMAh2(i1,iv2),ddMFd2(i3,iv1,iv2)& 
& ,ddMFd2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFd2(i2),MAh2(i1),dMFd2(i3,iv1)            & 
& ,dMFd2(i2,iv1),dMAh2(i1,iv1),dMFd2(i3,iv2),dMFd2(i2,iv2),dMAh2(i1,iv2),ddMFd2(i3,iv1,iv2)& 
& ,ddMFd2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 3
results1(14)=results1(14) + coeff*colorfactor*temp
results1_ti(14)=results1_ti(14) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(14)=results1(14) + coeffbar*colorfactor*tempbar
results1_ti(14)=results1_ti(14) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(14).eq.results1(14)))  write(*,*) 'NaN at FFS C[Ah, Fd, bar[Fd]]' 
! ---- Ah,Fu,bar[Fu] ----
Do i1=1,8
 Do i2=1,3
    Do i3=1,3
coup1L = cplcFuFuAhL(i3,i2,i1)
coup1R = cplcFuFuAhR(i3,i2,i1)
coup2L = cplcFuFuAhL(i2,i3,i1)
coup2R = cplcFuFuAhR(i2,i3,i1)
Di_coup1L = dcplcFuFuAhL(i3,i2,i1,iv1)
Di_coup1R = dcplcFuFuAhR(i3,i2,i1,iv1)
Dj_coup1L = dcplcFuFuAhL(i3,i2,i1,iv2)
Dj_coup1R = dcplcFuFuAhR(i3,i2,i1,iv2)
DDcoup1L = ddcplcFuFuAhL(i3,i2,i1,iv1,iv2)
DDcoup1R = ddcplcFuFuAhR(i3,i2,i1,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFu2(i3),MFu2(i2),MAh2(i1),dMFu2(i3,iv1)            & 
& ,dMFu2(i2,iv1),dMAh2(i1,iv1),dMFu2(i3,iv2),dMFu2(i2,iv2),dMAh2(i1,iv2),ddMFu2(i3,iv1,iv2)& 
& ,ddMFu2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFu2(i3),MFu2(i2),MAh2(i1),dMFu2(i3,iv1)            & 
& ,dMFu2(i2,iv1),dMAh2(i1,iv1),dMFu2(i3,iv2),dMFu2(i2,iv2),dMAh2(i1,iv2),ddMFu2(i3,iv1,iv2)& 
& ,ddMFu2(i2,iv1,iv2),ddMAh2(i1,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 3
results1(15)=results1(15) + coeff*colorfactor*temp
results1_ti(15)=results1_ti(15) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(15)=results1(15) + coeffbar*colorfactor*tempbar
results1_ti(15)=results1_ti(15) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(15).eq.results1(15)))  write(*,*) 'NaN at FFS C[Ah, Fu, bar[Fu]]' 
! ---- Cha,hh,bar[Cha] ----
Do i1=1,5
 Do i2=1,8
    Do i3=1,5
coup1L = cplcChaChahhL(i3,i1,i2)
coup1R = cplcChaChahhR(i3,i1,i2)
coup2L = cplcChaChahhL(i1,i3,i2)
coup2R = cplcChaChahhR(i1,i3,i2)
Di_coup1L = dcplcChaChahhL(i3,i1,i2,iv1)
Di_coup1R = dcplcChaChahhR(i3,i1,i2,iv1)
Dj_coup1L = dcplcChaChahhL(i3,i1,i2,iv2)
Dj_coup1R = dcplcChaChahhR(i3,i1,i2,iv2)
DDcoup1L = ddcplcChaChahhL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcChaChahhR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MCha2(i3),MCha2(i1),Mhh2(i2),dMCha2(i3,iv1)         & 
& ,dMCha2(i1,iv1),dMhh2(i2,iv1),dMCha2(i3,iv2),dMCha2(i1,iv2),dMhh2(i2,iv2)              & 
& ,ddMCha2(i3,iv1,iv2),ddMCha2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MCha2(i3),MCha2(i1),Mhh2(i2),dMCha2(i3,iv1)         & 
& ,dMCha2(i1,iv1),dMhh2(i2,iv1),dMCha2(i3,iv2),dMCha2(i1,iv2),dMhh2(i2,iv2)              & 
& ,ddMCha2(i3,iv1,iv2),ddMCha2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 1
results1(16)=results1(16) + coeff*colorfactor*temp
results1_ti(16)=results1_ti(16) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(16)=results1(16) + coeffbar*colorfactor*tempbar
results1_ti(16)=results1_ti(16) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(16).eq.results1(16)))  write(*,*) 'NaN at FFS C[Cha, hh, bar[Cha]]' 
! ---- Chi,Chi,hh ----
Do i1=1,10
 Do i2=1,10
    Do i3=1,8
coup1L = cplChiChihhL(i1,i2,i3)
coup1R = cplChiChihhR(i1,i2,i3)
coup2L = cplChiChihhL(i1,i2,i3)
coup2R = cplChiChihhR(i1,i2,i3)
Di_coup1L = dcplChiChihhL(i1,i2,i3,iv1)
Di_coup1R = dcplChiChihhR(i1,i2,i3,iv1)
Dj_coup1L = dcplChiChihhL(i1,i2,i3,iv2)
Dj_coup1R = dcplChiChihhR(i1,i2,i3,iv2)
DDcoup1L = ddcplChiChihhL(i1,i2,i3,iv1,iv2)
DDcoup1R = ddcplChiChihhR(i1,i2,i3,iv1,iv2)
coupx = abs(coup1L)**2
Di_coupx=Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) 
Dj_coupx=Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L)  
Call SecondDerivativeVeff_sunrise(MChi2(i2),MChi2(i1),Mhh2(i3),dMChi2(i2,iv1)         & 
& ,dMChi2(i1,iv1),dMhh2(i3,iv1),dMChi2(i2,iv2),dMChi2(i1,iv2),dMhh2(i3,iv2)              & 
& ,ddMChi2(i2,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMhh2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = Real(coup1L**2,dp) 
Di_coupxbar = Real(2*Di_coup1L*coup1L,dp) 
Dj_coupxbar = Real(2*Dj_coup1L*coup1L,dp) 
DDcoupxbar = Real(2*DDcoup1L*coup1L + 2*Di_coup1L*Dj_coup1L,dp) 
Call SecondDerivativeVeff_sunrise(MChi2(i2),MChi2(i1),Mhh2(i3),dMChi2(i2,iv1)         & 
& ,dMChi2(i1,iv1),dMhh2(i3,iv1),dMChi2(i2,iv2),dMChi2(i1,iv2),dMhh2(i3,iv2)              & 
& ,ddMChi2(i2,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMhh2(i3,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 1
results1(17)=results1(17) + coeff*colorfactor*temp
results1_ti(17)=results1_ti(17) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(17)=results1(17) + coeffbar*colorfactor*tempbar
results1_ti(17)=results1_ti(17) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(17).eq.results1(17)))  write(*,*) 'NaN at FFS C[Chi, Chi, hh]' 
! ---- Chi,Hpm,bar[Cha] ----
Do i1=1,10
 Do i2=1,8
    Do i3=1,5
coup1L = cplcChaChiHpmL(i3,i1,i2)
coup1R = cplcChaChiHpmR(i3,i1,i2)
coup2L = cplChiChacHpmL(i1,i3,i2)
coup2R = cplChiChacHpmR(i1,i3,i2)
Di_coup1L = dcplcChaChiHpmL(i3,i1,i2,iv1)
Di_coup1R = dcplcChaChiHpmR(i3,i1,i2,iv1)
Dj_coup1L = dcplcChaChiHpmL(i3,i1,i2,iv2)
Dj_coup1R = dcplcChaChiHpmR(i3,i1,i2,iv2)
DDcoup1L = ddcplcChaChiHpmL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcChaChiHpmR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MCha2(i3),MChi2(i1),MHpm2(i2),dMCha2(i3,iv1)        & 
& ,dMChi2(i1,iv1),dMHpm2(i2,iv1),dMCha2(i3,iv2),dMChi2(i1,iv2),dMHpm2(i2,iv2)            & 
& ,ddMCha2(i3,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MCha2(i3),MChi2(i1),MHpm2(i2),dMCha2(i3,iv1)        & 
& ,dMChi2(i1,iv1),dMHpm2(i2,iv1),dMCha2(i3,iv2),dMChi2(i1,iv2),dMHpm2(i2,iv2)            & 
& ,ddMCha2(i3,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 1
results1(18)=results1(18) + coeff*colorfactor*temp
results1_ti(18)=results1_ti(18) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(18)=results1(18) + coeffbar*colorfactor*tempbar
results1_ti(18)=results1_ti(18) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(18).eq.results1(18)))  write(*,*) 'NaN at FFS C[Chi, Hpm, bar[Cha]]' 
! ---- Chi,Sd,bar[Fd] ----
Do i1=1,10
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFdChiSdL(i3,i1,i2)
coup1R = cplcFdChiSdR(i3,i1,i2)
coup2L = cplChiFdcSdL(i1,i3,i2)
coup2R = cplChiFdcSdR(i1,i3,i2)
Di_coup1L = dcplcFdChiSdL(i3,i1,i2,iv1)
Di_coup1R = dcplcFdChiSdR(i3,i1,i2,iv1)
Dj_coup1L = dcplcFdChiSdL(i3,i1,i2,iv2)
Dj_coup1R = dcplcFdChiSdR(i3,i1,i2,iv2)
DDcoup1L = ddcplcFdChiSdL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcFdChiSdR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFd2(i3),MChi2(i1),MSd2(i2),dMFd2(i3,iv1)           & 
& ,dMChi2(i1,iv1),dMSd2(i2,iv1),dMFd2(i3,iv2),dMChi2(i1,iv2),dMSd2(i2,iv2)               & 
& ,ddMFd2(i3,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMSd2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFd2(i3),MChi2(i1),MSd2(i2),dMFd2(i3,iv1)           & 
& ,dMChi2(i1,iv1),dMSd2(i2,iv1),dMFd2(i3,iv2),dMChi2(i1,iv2),dMSd2(i2,iv2)               & 
& ,ddMFd2(i3,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMSd2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 3
results1(19)=results1(19) + coeff*colorfactor*temp
results1_ti(19)=results1_ti(19) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(19)=results1(19) + coeffbar*colorfactor*tempbar
results1_ti(19)=results1_ti(19) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(19).eq.results1(19)))  write(*,*) 'NaN at FFS C[Chi, Sd, bar[Fd]]' 
! ---- Chi,Su,bar[Fu] ----
Do i1=1,10
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFuChiSuL(i3,i1,i2)
coup1R = cplcFuChiSuR(i3,i1,i2)
coup2L = cplChiFucSuL(i1,i3,i2)
coup2R = cplChiFucSuR(i1,i3,i2)
Di_coup1L = dcplcFuChiSuL(i3,i1,i2,iv1)
Di_coup1R = dcplcFuChiSuR(i3,i1,i2,iv1)
Dj_coup1L = dcplcFuChiSuL(i3,i1,i2,iv2)
Dj_coup1R = dcplcFuChiSuR(i3,i1,i2,iv2)
DDcoup1L = ddcplcFuChiSuL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcFuChiSuR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFu2(i3),MChi2(i1),MSu2(i2),dMFu2(i3,iv1)           & 
& ,dMChi2(i1,iv1),dMSu2(i2,iv1),dMFu2(i3,iv2),dMChi2(i1,iv2),dMSu2(i2,iv2)               & 
& ,ddMFu2(i3,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFu2(i3),MChi2(i1),MSu2(i2),dMFu2(i3,iv1)           & 
& ,dMChi2(i1,iv1),dMSu2(i2,iv1),dMFu2(i3,iv2),dMChi2(i1,iv2),dMSu2(i2,iv2)               & 
& ,ddMFu2(i3,iv1,iv2),ddMChi2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 3
results1(20)=results1(20) + coeff*colorfactor*temp
results1_ti(20)=results1_ti(20) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(20)=results1(20) + coeffbar*colorfactor*tempbar
results1_ti(20)=results1_ti(20) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(20).eq.results1(20)))  write(*,*) 'NaN at FFS C[Chi, Su, bar[Fu]]' 
! ---- Fd,hh,bar[Fd] ----
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFdFdhhL(i3,i1,i2)
coup1R = cplcFdFdhhR(i3,i1,i2)
coup2L = cplcFdFdhhL(i1,i3,i2)
coup2R = cplcFdFdhhR(i1,i3,i2)
Di_coup1L = dcplcFdFdhhL(i3,i1,i2,iv1)
Di_coup1R = dcplcFdFdhhR(i3,i1,i2,iv1)
Dj_coup1L = dcplcFdFdhhL(i3,i1,i2,iv2)
Dj_coup1R = dcplcFdFdhhR(i3,i1,i2,iv2)
DDcoup1L = ddcplcFdFdhhL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcFdFdhhR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),Mhh2(i2),dMFd2(i3,iv1)            & 
& ,dMFd2(i1,iv1),dMhh2(i2,iv1),dMFd2(i3,iv2),dMFd2(i1,iv2),dMhh2(i2,iv2),ddMFd2(i3,iv1,iv2)& 
& ,ddMFd2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),Mhh2(i2),dMFd2(i3,iv1)            & 
& ,dMFd2(i1,iv1),dMhh2(i2,iv1),dMFd2(i3,iv2),dMFd2(i1,iv2),dMhh2(i2,iv2),ddMFd2(i3,iv1,iv2)& 
& ,ddMFd2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 3
results1(21)=results1(21) + coeff*colorfactor*temp
results1_ti(21)=results1_ti(21) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(21)=results1(21) + coeffbar*colorfactor*tempbar
results1_ti(21)=results1_ti(21) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(21).eq.results1(21)))  write(*,*) 'NaN at FFS C[Fd, hh, bar[Fd]]' 
! ---- Fd,bar[Cha],conj[Su] ----
Do i1=1,3
 Do i2=1,5
    Do i3=1,6
coup1L = cplcChaFdcSuL(i2,i1,i3)
coup1R = cplcChaFdcSuR(i2,i1,i3)
coup2L = cplcFdChaSuL(i1,i2,i3)
coup2R = cplcFdChaSuR(i1,i2,i3)
Di_coup1L = dcplcChaFdcSuL(i2,i1,i3,iv1)
Di_coup1R = dcplcChaFdcSuR(i2,i1,i3,iv1)
Dj_coup1L = dcplcChaFdcSuL(i2,i1,i3,iv2)
Dj_coup1R = dcplcChaFdcSuR(i2,i1,i3,iv2)
DDcoup1L = ddcplcChaFdcSuL(i2,i1,i3,iv1,iv2)
DDcoup1R = ddcplcChaFdcSuR(i2,i1,i3,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MCha2(i2),MFd2(i1),MSu2(i3),dMCha2(i2,iv1)          & 
& ,dMFd2(i1,iv1),dMSu2(i3,iv1),dMCha2(i2,iv2),dMFd2(i1,iv2),dMSu2(i3,iv2),ddMCha2(i2,iv1,iv2)& 
& ,ddMFd2(i1,iv1,iv2),ddMSu2(i3,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MCha2(i2),MFd2(i1),MSu2(i3),dMCha2(i2,iv1)          & 
& ,dMFd2(i1,iv1),dMSu2(i3,iv1),dMCha2(i2,iv2),dMFd2(i1,iv2),dMSu2(i3,iv2),ddMCha2(i2,iv1,iv2)& 
& ,ddMFd2(i1,iv1,iv2),ddMSu2(i3,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 3
results1(22)=results1(22) + coeff*colorfactor*temp
results1_ti(22)=results1_ti(22) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(22)=results1(22) + coeffbar*colorfactor*tempbar
results1_ti(22)=results1_ti(22) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(22).eq.results1(22)))  write(*,*) 'NaN at FFS C[Fd, bar[Cha], conj[Su]]' 
! ---- Fu,hh,bar[Fu] ----
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFuFuhhL(i3,i1,i2)
coup1R = cplcFuFuhhR(i3,i1,i2)
coup2L = cplcFuFuhhL(i1,i3,i2)
coup2R = cplcFuFuhhR(i1,i3,i2)
Di_coup1L = dcplcFuFuhhL(i3,i1,i2,iv1)
Di_coup1R = dcplcFuFuhhR(i3,i1,i2,iv1)
Dj_coup1L = dcplcFuFuhhL(i3,i1,i2,iv2)
Dj_coup1R = dcplcFuFuhhR(i3,i1,i2,iv2)
DDcoup1L = ddcplcFuFuhhL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcFuFuhhR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),Mhh2(i2),dMFu2(i3,iv1)            & 
& ,dMFu2(i1,iv1),dMhh2(i2,iv1),dMFu2(i3,iv2),dMFu2(i1,iv2),dMhh2(i2,iv2),ddMFu2(i3,iv1,iv2)& 
& ,ddMFu2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),Mhh2(i2),dMFu2(i3,iv1)            & 
& ,dMFu2(i1,iv1),dMhh2(i2,iv1),dMFu2(i3,iv2),dMFu2(i1,iv2),dMhh2(i2,iv2),ddMFu2(i3,iv1,iv2)& 
& ,ddMFu2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 3
results1(23)=results1(23) + coeff*colorfactor*temp
results1_ti(23)=results1_ti(23) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(23)=results1(23) + coeffbar*colorfactor*tempbar
results1_ti(23)=results1_ti(23) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(23).eq.results1(23)))  write(*,*) 'NaN at FFS C[Fu, hh, bar[Fu]]' 
! ---- Fu,Hpm,bar[Fd] ----
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFdFuHpmL(i3,i1,i2)
coup1R = cplcFdFuHpmR(i3,i1,i2)
coup2L = cplcFuFdcHpmL(i1,i3,i2)
coup2R = cplcFuFdcHpmR(i1,i3,i2)
Di_coup1L = dcplcFdFuHpmL(i3,i1,i2,iv1)
Di_coup1R = dcplcFdFuHpmR(i3,i1,i2,iv1)
Dj_coup1L = dcplcFdFuHpmL(i3,i1,i2,iv2)
Dj_coup1R = dcplcFdFuHpmR(i3,i1,i2,iv2)
DDcoup1L = ddcplcFdFuHpmL(i3,i1,i2,iv1,iv2)
DDcoup1R = ddcplcFdFuHpmR(i3,i1,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFu2(i1),MHpm2(i2),dMFd2(i3,iv1)           & 
& ,dMFu2(i1,iv1),dMHpm2(i2,iv1),dMFd2(i3,iv2),dMFu2(i1,iv2),dMHpm2(i2,iv2)               & 
& ,ddMFd2(i3,iv1,iv2),ddMFu2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFu2(i1),MHpm2(i2),dMFd2(i3,iv1)           & 
& ,dMFu2(i1,iv1),dMHpm2(i2,iv1),dMFd2(i3,iv2),dMFu2(i1,iv2),dMHpm2(i2,iv2)               & 
& ,ddMFd2(i3,iv1,iv2),ddMFu2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 3
results1(24)=results1(24) + coeff*colorfactor*temp
results1_ti(24)=results1_ti(24) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(24)=results1(24) + coeffbar*colorfactor*tempbar
results1_ti(24)=results1_ti(24) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(24).eq.results1(24)))  write(*,*) 'NaN at FFS C[Fu, Hpm, bar[Fd]]' 
! ---- Glu,Sd,bar[Fd] ----
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFdGluSdL(i3,i2)
coup1R = cplcFdGluSdR(i3,i2)
coup2L = cplGluFdcSdL(i3,i2)
coup2R = cplGluFdcSdR(i3,i2)
Di_coup1L = dcplcFdGluSdL(i3,i2,iv1)
Di_coup1R = dcplcFdGluSdR(i3,i2,iv1)
Dj_coup1L = dcplcFdGluSdL(i3,i2,iv2)
Dj_coup1R = dcplcFdGluSdR(i3,i2,iv2)
DDcoup1L = ddcplcFdGluSdL(i3,i2,iv1,iv2)
DDcoup1R = ddcplcFdGluSdR(i3,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFd2(i3),MGlu2,MSd2(i2),dMFd2(i3,iv1)               & 
& ,dMGlu2(1,iv1),dMSd2(i2,iv1),dMFd2(i3,iv2),dMGlu2(1,iv2),dMSd2(i2,iv2),ddMFd2(i3,iv1,iv2)& 
& ,ddMGlu2(1,iv1,iv2),ddMSd2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFd2(i3),MGlu2,MSd2(i2),dMFd2(i3,iv1)               & 
& ,dMGlu2(1,iv1),dMSd2(i2,iv1),dMFd2(i3,iv2),dMGlu2(1,iv2),dMSd2(i2,iv2),ddMFd2(i3,iv1,iv2)& 
& ,ddMGlu2(1,iv1,iv2),ddMSd2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 4
results1(25)=results1(25) + coeff*colorfactor*temp
results1_ti(25)=results1_ti(25) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(25)=results1(25) + coeffbar*colorfactor*tempbar
results1_ti(25)=results1_ti(25) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
if (.not.(results1(25).eq.results1(25)))  write(*,*) 'NaN at FFS C[Glu, Sd, bar[Fd]]' 
! ---- Glu,Su,bar[Fu] ----
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFuGluSuL(i3,i2)
coup1R = cplcFuGluSuR(i3,i2)
coup2L = cplGluFucSuL(i3,i2)
coup2R = cplGluFucSuR(i3,i2)
Di_coup1L = dcplcFuGluSuL(i3,i2,iv1)
Di_coup1R = dcplcFuGluSuR(i3,i2,iv1)
Dj_coup1L = dcplcFuGluSuL(i3,i2,iv2)
Dj_coup1R = dcplcFuGluSuR(i3,i2,iv2)
DDcoup1L = ddcplcFuGluSuL(i3,i2,iv1,iv2)
DDcoup1R = ddcplcFuGluSuR(i3,i2,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFu2(i3),MGlu2,MSu2(i2),dMFu2(i3,iv1)               & 
& ,dMGlu2(1,iv1),dMSu2(i2,iv1),dMFu2(i3,iv2),dMGlu2(1,iv2),dMSu2(i2,iv2),ddMFu2(i3,iv1,iv2)& 
& ,ddMGlu2(1,iv1,iv2),ddMSu2(i2,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFu2(i3),MGlu2,MSu2(i2),dMFu2(i3,iv1)               & 
& ,dMGlu2(1,iv1),dMSu2(i2,iv1),dMFu2(i3,iv2),dMGlu2(1,iv2),dMSu2(i2,iv2),ddMFu2(i3,iv1,iv2)& 
& ,ddMGlu2(1,iv1,iv2),ddMSu2(i2,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 4
results1(26)=results1(26) + coeff*colorfactor*temp
results1_ti(26)=results1_ti(26) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(26)=results1(26) + coeffbar*colorfactor*tempbar
results1_ti(26)=results1_ti(26) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
if (.not.(results1(26).eq.results1(26)))  write(*,*) 'NaN at FFS C[Glu, Su, bar[Fu]]' 
! ---- Sd,bar[Cha],bar[Fu] ----
Do i1=1,6
 Do i2=1,5
    Do i3=1,3
coup1L = cplcChacFuSdL(i2,i3,i1)
coup1R = cplcChacFuSdR(i2,i3,i1)
coup2L = cplChaFucSdL(i2,i3,i1)
coup2R = cplChaFucSdR(i2,i3,i1)
Di_coup1L = dcplcChacFuSdL(i2,i3,i1,iv1)
Di_coup1R = dcplcChacFuSdR(i2,i3,i1,iv1)
Dj_coup1L = dcplcChacFuSdL(i2,i3,i1,iv2)
Dj_coup1R = dcplcChacFuSdR(i2,i3,i1,iv2)
DDcoup1L = ddcplcChacFuSdL(i2,i3,i1,iv1,iv2)
DDcoup1R = ddcplcChacFuSdR(i2,i3,i1,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFu2(i3),MCha2(i2),MSd2(i1),dMFu2(i3,iv1)           & 
& ,dMCha2(i2,iv1),dMSd2(i1,iv1),dMFu2(i3,iv2),dMCha2(i2,iv2),dMSd2(i1,iv2)               & 
& ,ddMFu2(i3,iv1,iv2),ddMCha2(i2,iv1,iv2),ddMSd2(i1,iv1,iv2),coupx,Di_coupx,Dj_coupx,DDcoupx,'FFS   ',Q2,temp,temp_ti,temp_tj)
coupxbar = 2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = 2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = 2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = 2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFu2(i3),MCha2(i2),MSd2(i1),dMFu2(i3,iv1)           & 
& ,dMCha2(i2,iv1),dMSd2(i1,iv1),dMFu2(i3,iv2),dMCha2(i2,iv2),dMSd2(i1,iv2)               & 
& ,ddMFu2(i3,iv1,iv2),ddMCha2(i2,iv1,iv2),ddMSd2(i1,iv1,iv2),coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFSbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 1._dp
colorfactor = 3
results1(27)=results1(27) + coeff*colorfactor*temp
results1_ti(27)=results1_ti(27) + coeff*colorfactor*temp_ti
coeffbar = 1._dp
results1(27)=results1(27) + coeffbar*colorfactor*tempbar
results1_ti(27)=results1_ti(27) + coeffbar*colorfactor*tempbar_ti
    End Do
  End Do
End Do
if (.not.(results1(27).eq.results1(27)))  write(*,*) 'NaN at FFS C[Sd, bar[Cha], bar[Fu]]' 

! ----- diagrams of type SSV, 2 ------ 

! ---- Sd,VG,conj[Sd] ----
Do i1=1,6
    Do i3=1,6
coup1 = cplSdcSdVG(i1,i3)
coup2 = cplSdcSdVG(i3,i1)
Di_coup1 = dcplSdcSdVG(i1,i3,iv1)
Dj_coup1 = dcplSdcSdVG(i1,i3,iv2)
DDcoup1 = ddcplSdcSdVG(i1,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MSd2(i3),MSd2(i1),0._dp,dMSd2(i3,iv1)               & 
& ,dMSd2(i1,iv1),ZeroC,dMSd2(i3,iv2),dMSd2(i1,iv2),ZeroC,ddMSd2(i3,iv1,iv2)              & 
& ,ddMSd2(i1,iv1,iv2),ZeroC,coupx,Di_coupx,Dj_coupx,DDcoupx,'SSV   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 4
results1(28)=results1(28) + coeff*colorfactor*temp
results1_ti(28)=results1_ti(28) + coeff*colorfactor*temp_ti
    End Do
End Do
if (.not.(results1(28).eq.results1(28)))  write(*,*) 'NaN at SSV C[Sd, VG, conj[Sd]]' 
! ---- Su,VG,conj[Su] ----
Do i1=1,6
    Do i3=1,6
coup1 = cplSucSuVG(i1,i3)
coup2 = cplSucSuVG(i3,i1)
Di_coup1 = dcplSucSuVG(i1,i3,iv1)
Dj_coup1 = dcplSucSuVG(i1,i3,iv2)
DDcoup1 = ddcplSucSuVG(i1,i3,iv1,iv2)
coupx=abs(coup1)**2 
Di_coupx=Di_coup1*conjg(coup1)+coup1*conjg(Di_coup1) 
Dj_coupx=Dj_coup1*conjg(coup1)+coup1*conjg(Dj_coup1) 
DDcoupx = DDcoup1*conjg(coup1)+coup1*conjg(DDcoup1) & 
& + Di_coup1*conjg(Dj_coup1)+Dj_coup1*conjg(Di_coup1)  
Call SecondDerivativeVeff_sunrise(MSu2(i3),MSu2(i1),0._dp,dMSu2(i3,iv1)               & 
& ,dMSu2(i1,iv1),ZeroC,dMSu2(i3,iv2),dMSu2(i1,iv2),ZeroC,ddMSu2(i3,iv1,iv2)              & 
& ,ddMSu2(i1,iv1,iv2),ZeroC,coupx,Di_coupx,Dj_coupx,DDcoupx,'SSV   ',Q2,temp,temp_ti,temp_tj)
coeff = 0.5_dp
colorfactor = 4
results1(29)=results1(29) + coeff*colorfactor*temp
results1_ti(29)=results1_ti(29) + coeff*colorfactor*temp_ti
    End Do
End Do
if (.not.(results1(29).eq.results1(29)))  write(*,*) 'NaN at SSV C[Su, VG, conj[Su]]' 

! ----- diagrams of type FFV, 3 ------ 

! ---- Fd,VG,bar[Fd] ----
Do i1=1,3
    Do i3=1,3
coup1L = cplcFdFdVGL(i3,i1)
coup1R = cplcFdFdVGR(i3,i1)
coup2L = cplcFdFdVGL(i1,i3)
coup2R = cplcFdFdVGR(i1,i3)
Di_coup1L = dcplcFdFdVGL(i3,i1,iv1)
Di_coup1R = dcplcFdFdVGR(i3,i1,iv1)
Dj_coup1L = dcplcFdFdVGL(i3,i1,iv2)
Dj_coup1R = dcplcFdFdVGR(i3,i1,iv2)
DDcoup1L = ddcplcFdFdVGL(i3,i1,iv1,iv2)
DDcoup1R = ddcplcFdFdVGR(i3,i1,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),0._dp,dMFd2(i3,iv1)               & 
& ,dMFd2(i1,iv1),ZeroC,dMFd2(i3,iv2),dMFd2(i1,iv2),ZeroC,ddMFd2(i3,iv1,iv2)              & 
& ,ddMFd2(i1,iv1,iv2),ZeroC,coupx,Di_coupx,Dj_coupx,DDcoupx,'FFV   ',Q2,temp,temp_ti,temp_tj)
coupxbar = -2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = -2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = -2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = -2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),0._dp,dMFd2(i3,iv1)               & 
& ,dMFd2(i1,iv1),ZeroC,dMFd2(i3,iv2),dMFd2(i1,iv2),ZeroC,ddMFd2(i3,iv1,iv2)              & 
& ,ddMFd2(i1,iv1,iv2),ZeroC,coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFVbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 4
results1(30)=results1(30) + coeff*colorfactor*temp
results1_ti(30)=results1_ti(30) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(30)=results1(30) + coeffbar*colorfactor*tempbar
results1_ti(30)=results1_ti(30) + coeffbar*colorfactor*tempbar_ti
    End Do
End Do
if (.not.(results1(30).eq.results1(30)))  write(*,*) 'NaN at FFV C[Fd, VG, bar[Fd]]' 
! ---- Fu,VG,bar[Fu] ----
Do i1=1,3
    Do i3=1,3
coup1L = cplcFuFuVGL(i3,i1)
coup1R = cplcFuFuVGR(i3,i1)
coup2L = cplcFuFuVGL(i1,i3)
coup2R = cplcFuFuVGR(i1,i3)
Di_coup1L = dcplcFuFuVGL(i3,i1,iv1)
Di_coup1R = dcplcFuFuVGR(i3,i1,iv1)
Dj_coup1L = dcplcFuFuVGL(i3,i1,iv2)
Dj_coup1R = dcplcFuFuVGR(i3,i1,iv2)
DDcoup1L = ddcplcFuFuVGL(i3,i1,iv1,iv2)
DDcoup1R = ddcplcFuFuVGR(i3,i1,iv1,iv2)
coupx = (abs(coup1L)**2 + abs(coup1R)**2) 
Di_coupx = Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) &
& + Di_coup1R*conjg(coup1R)+coup1R*conjg(Di_coup1R) 
Dj_coupx = Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) &
& + Dj_coup1R*conjg(coup1R)+coup1R*conjg(Dj_coup1R) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L) & 
& + DDcoup1R*conjg(coup1R)+coup1R*conjg(DDcoup1R) & 
& + Di_coup1R*conjg(Dj_coup1R)+Dj_coup1R*conjg(Di_coup1R)  
Call SecondDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),0._dp,dMFu2(i3,iv1)               & 
& ,dMFu2(i1,iv1),ZeroC,dMFu2(i3,iv2),dMFu2(i1,iv2),ZeroC,ddMFu2(i3,iv1,iv2)              & 
& ,ddMFu2(i1,iv1,iv2),ZeroC,coupx,Di_coupx,Dj_coupx,DDcoupx,'FFV   ',Q2,temp,temp_ti,temp_tj)
coupxbar = -2*Real(coup1L*conjg(coup1R),dp) 
Di_coupxbar = -2*Real(Di_coup1L*conjg(coup1R)+coup1L*conjg(Di_coup1R),dp) 
Dj_coupxbar = -2*Real(Dj_coup1L*conjg(coup1R)+coup1L*conjg(Dj_coup1R),dp) 
DDcoupxbar = -2*Real(DDcoup1L*conjg(coup1R)+coup1L*conjg(DDcoup1R) &  
&          + Di_coup1L*conjg(Dj_coup1R)+Dj_coup1L*conjg(Di_coup1R)  ,dp) 
Call SecondDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),0._dp,dMFu2(i3,iv1)               & 
& ,dMFu2(i1,iv1),ZeroC,dMFu2(i3,iv2),dMFu2(i1,iv2),ZeroC,ddMFu2(i3,iv1,iv2)              & 
& ,ddMFu2(i1,iv1,iv2),ZeroC,coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFVbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 4
results1(31)=results1(31) + coeff*colorfactor*temp
results1_ti(31)=results1_ti(31) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(31)=results1(31) + coeffbar*colorfactor*tempbar
results1_ti(31)=results1_ti(31) + coeffbar*colorfactor*tempbar_ti
    End Do
End Do
if (.not.(results1(31).eq.results1(31)))  write(*,*) 'NaN at FFV C[Fu, VG, bar[Fu]]' 
! ---- Glu,Glu,VG ----
coup1L = cplGluGluVGL
coup1R = cplGluGluVGR
coup2L = cplGluGluVGL
coup2R = cplGluGluVGR
Di_coup1L = dcplGluGluVGL(iv1)
Di_coup1R = dcplGluGluVGR(iv1)
Dj_coup1L = dcplGluGluVGL(iv2)
Dj_coup1R = dcplGluGluVGR(iv2)
DDcoup1L = ddcplGluGluVGL(iv1,iv2)
DDcoup1R = ddcplGluGluVGR(iv1,iv2)
coupx = abs(coup1L)**2
Di_coupx=Di_coup1L*conjg(coup1L)+coup1L*conjg(Di_coup1L) 
Dj_coupx=Dj_coup1L*conjg(coup1L)+coup1L*conjg(Dj_coup1L) 
DDcoupx = DDcoup1L*conjg(coup1L)+coup1L*conjg(DDcoup1L) & 
& + Di_coup1L*conjg(Dj_coup1L)+Dj_coup1L*conjg(Di_coup1L)  
Call SecondDerivativeVeff_sunrise(MGlu2,MGlu2,0._dp,dMGlu2(1,iv1),dMGlu2(1,iv1)       & 
& ,ZeroC,dMGlu2(1,iv2),dMGlu2(1,iv2),ZeroC,ddMGlu2(1,iv1,iv2),ddMGlu2(1,iv1,iv2)         & 
& ,ZeroC,coupx,Di_coupx,Dj_coupx,DDcoupx,'FFV   ',Q2,temp,temp_ti,temp_tj)
coupxbar = Real(coup1L**2,dp) 
Di_coupxbar = Real(2*Di_coup1L*coup1L,dp) 
Dj_coupxbar = Real(2*Dj_coup1L*coup1L,dp) 
DDcoupxbar = Real(2*DDcoup1L*coup1L + 2*Di_coup1L*Dj_coup1L,dp) 
Call SecondDerivativeVeff_sunrise(MGlu2,MGlu2,0._dp,dMGlu2(1,iv1),dMGlu2(1,iv1)       & 
& ,ZeroC,dMGlu2(1,iv2),dMGlu2(1,iv2),ZeroC,ddMGlu2(1,iv1,iv2),ddMGlu2(1,iv1,iv2)         & 
& ,ZeroC,coupxbar,Di_coupxbar,Dj_coupxbar,DDcoupxbar,'FFVbar',Q2,tempbar,tempbar_ti,tempbar_tj)
coeff = 0.5_dp
colorfactor = 24
results1(32)=results1(32) + coeff*colorfactor*temp
results1_ti(32)=results1_ti(32) + coeff*colorfactor*temp_ti
coeffbar = 0.5_dp
results1(32)=results1(32) + coeffbar*colorfactor*tempbar
results1_ti(32)=results1_ti(32) + coeffbar*colorfactor*tempbar_ti
if (.not.(results1(32).eq.results1(32)))  write(*,*) 'NaN at FFV C[Glu, Glu, VG]' 

! ----- diagrams of type VVV, 1 ------ 

! ---- VG,VG,VG ----
coup1 = cplVGVGVG
coup2 = cplVGVGVG
Di_coup1 = dcplVGVGVG(iv1)
Dj_coup1 = dcplVGVGVG(iv2)
DDcoup1 = ddcplVGVGVG(iv1,iv2)
coeff = 0.000
colorfactor = 24
results1(33)=results1(33) + coeff*colorfactor*temp
results1_ti(33)=results1_ti(33) + coeff*colorfactor*temp_ti
if (.not.(results1(33).eq.results1(33)))  write(*,*) 'NaN at VVV C[VG, VG, VG]' 
! ----- Topology2: diagrams w. 2 Particles and 1 Vertex


! ----- diagrams of type SS, 15 ------ 

! ---- Ah,Ah ----
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhAhAh(i1,i1,i2,i2)
Di_coup1 = dcplAhAhAhAh(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplAhAhAhAh(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplAhAhAhAh(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MAh2(i1),MAh2(i2),dMAh2(i1,iv1),dMAh2(i2,iv1)         & 
& ,dMAh2(i1,iv2),dMAh2(i2,iv2),ddMAh2(i1,iv1,iv2),ddMAh2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-1._dp/8._dp)
results2(1)=results2(1) + coeff*temp
results2_ti(1)=results2_ti(1) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(1).eq.results2(1)))  write(*,*) 'NaN at SS C[Ah, Ah, Ah, Ah]' 
! ---- Ah,hh ----
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhhhhh(i1,i1,i2,i2)
Di_coup1 = dcplAhAhhhhh(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplAhAhhhhh(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplAhAhhhhh(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MAh2(i1),Mhh2(i2),dMAh2(i1,iv1),dMhh2(i2,iv1)         & 
& ,dMAh2(i1,iv2),dMhh2(i2,iv2),ddMAh2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.25_dp)
results2(2)=results2(2) + coeff*temp
results2_ti(2)=results2_ti(2) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(2).eq.results2(2)))  write(*,*) 'NaN at SS C[Ah, Ah, hh, hh]' 
! ---- Ah,Hpm ----
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhHpmcHpm(i1,i1,i2,i2)
Di_coup1 = dcplAhAhHpmcHpm(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplAhAhHpmcHpm(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplAhAhHpmcHpm(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MAh2(i1),MHpm2(i2),dMAh2(i1,iv1),dMHpm2(i2,iv1)       & 
& ,dMAh2(i1,iv2),dMHpm2(i2,iv2),ddMAh2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(3)=results2(3) + coeff*temp
results2_ti(3)=results2_ti(3) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(3).eq.results2(3)))  write(*,*) 'NaN at SS C[Ah, Ah, Hpm, conj[Hpm]]' 
! ---- Ah,Sd ----
Do i1=1,8
 Do i2=1,6
coup1 = cplAhAhSdcSd(i1,i1,i2,i2)
Di_coup1 = dcplAhAhSdcSd(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplAhAhSdcSd(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplAhAhSdcSd(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MAh2(i1),MSd2(i2),dMAh2(i1,iv1),dMSd2(i2,iv1)         & 
& ,dMAh2(i1,iv2),dMSd2(i2,iv2),ddMAh2(i1,iv1,iv2),ddMSd2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(4)=results2(4) + coeff*temp
results2_ti(4)=results2_ti(4) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(4).eq.results2(4)))  write(*,*) 'NaN at SS C[Ah, Ah, Sd, conj[Sd]]' 
! ---- Ah,Su ----
Do i1=1,8
 Do i2=1,6
coup1 = cplAhAhSucSu(i1,i1,i2,i2)
Di_coup1 = dcplAhAhSucSu(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplAhAhSucSu(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplAhAhSucSu(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MAh2(i1),MSu2(i2),dMAh2(i1,iv1),dMSu2(i2,iv1)         & 
& ,dMAh2(i1,iv2),dMSu2(i2,iv2),ddMAh2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(5)=results2(5) + coeff*temp
results2_ti(5)=results2_ti(5) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(5).eq.results2(5)))  write(*,*) 'NaN at SS C[Ah, Ah, Su, conj[Su]]' 
! ---- hh,hh ----
Do i1=1,8
 Do i2=1,8
coup1 = cplhhhhhhhh(i1,i1,i2,i2)
Di_coup1 = dcplhhhhhhhh(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplhhhhhhhh(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplhhhhhhhh(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(Mhh2(i1),Mhh2(i2),dMhh2(i1,iv1),dMhh2(i2,iv1)         & 
& ,dMhh2(i1,iv2),dMhh2(i2,iv2),ddMhh2(i1,iv1,iv2),ddMhh2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-1._dp/8._dp)
results2(6)=results2(6) + coeff*temp
results2_ti(6)=results2_ti(6) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(6).eq.results2(6)))  write(*,*) 'NaN at SS C[hh, hh, hh, hh]' 
! ---- hh,Hpm ----
Do i1=1,8
 Do i2=1,8
coup1 = cplhhhhHpmcHpm(i1,i1,i2,i2)
Di_coup1 = dcplhhhhHpmcHpm(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplhhhhHpmcHpm(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplhhhhHpmcHpm(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(Mhh2(i1),MHpm2(i2),dMhh2(i1,iv1),dMHpm2(i2,iv1)       & 
& ,dMhh2(i1,iv2),dMHpm2(i2,iv2),ddMhh2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(7)=results2(7) + coeff*temp
results2_ti(7)=results2_ti(7) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(7).eq.results2(7)))  write(*,*) 'NaN at SS C[hh, hh, Hpm, conj[Hpm]]' 
! ---- hh,Sd ----
Do i1=1,8
 Do i2=1,6
coup1 = cplhhhhSdcSd(i1,i1,i2,i2)
Di_coup1 = dcplhhhhSdcSd(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplhhhhSdcSd(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplhhhhSdcSd(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(Mhh2(i1),MSd2(i2),dMhh2(i1,iv1),dMSd2(i2,iv1)         & 
& ,dMhh2(i1,iv2),dMSd2(i2,iv2),ddMhh2(i1,iv1,iv2),ddMSd2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(8)=results2(8) + coeff*temp
results2_ti(8)=results2_ti(8) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(8).eq.results2(8)))  write(*,*) 'NaN at SS C[hh, hh, Sd, conj[Sd]]' 
! ---- hh,Su ----
Do i1=1,8
 Do i2=1,6
coup1 = cplhhhhSucSu(i1,i1,i2,i2)
Di_coup1 = dcplhhhhSucSu(i1,i1,i2,i2,iv1)
Dj_coup1 = dcplhhhhSucSu(i1,i1,i2,i2,iv2)
DDcoup1 = ddcplhhhhSucSu(i1,i1,i2,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(Mhh2(i1),MSu2(i2),dMhh2(i1,iv1),dMSu2(i2,iv1)         & 
& ,dMhh2(i1,iv2),dMSu2(i2,iv2),ddMhh2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(9)=results2(9) + coeff*temp
results2_ti(9)=results2_ti(9) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(9).eq.results2(9)))  write(*,*) 'NaN at SS C[hh, hh, Su, conj[Su]]' 
! ---- Hpm,Hpm ----
Do i1=1,8
 Do i2=1,8
coup1 = cplHpmHpmcHpmcHpm(i1,i2,i1,i2)
Di_coup1 = dcplHpmHpmcHpmcHpm(i1,i2,i1,i2,iv1)
Dj_coup1 = dcplHpmHpmcHpmcHpm(i1,i2,i1,i2,iv2)
DDcoup1 = ddcplHpmHpmcHpmcHpm(i1,i2,i1,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MHpm2(i1),MHpm2(i2),dMHpm2(i1,iv1),dMHpm2(i2,iv1)     & 
& ,dMHpm2(i1,iv2),dMHpm2(i2,iv2),ddMHpm2(i1,iv1,iv2),ddMHpm2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(10)=results2(10) + coeff*temp
results2_ti(10)=results2_ti(10) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(10).eq.results2(10)))  write(*,*) 'NaN at SS C[Hpm, Hpm, conj[Hpm], conj[Hpm]]' 
! ---- Hpm,Sd ----
Do i1=1,8
 Do i2=1,6
coup1 = cplHpmSdcHpmcSd(i1,i2,i1,i2)
Di_coup1 = dcplHpmSdcHpmcSd(i1,i2,i1,i2,iv1)
Dj_coup1 = dcplHpmSdcHpmcSd(i1,i2,i1,i2,iv2)
DDcoup1 = ddcplHpmSdcHpmcSd(i1,i2,i1,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MHpm2(i1),MSd2(i2),dMHpm2(i1,iv1),dMSd2(i2,iv1)       & 
& ,dMHpm2(i1,iv2),dMSd2(i2,iv2),ddMHpm2(i1,iv1,iv2),ddMSd2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-1._dp)
results2(11)=results2(11) + coeff*temp
results2_ti(11)=results2_ti(11) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(11).eq.results2(11)))  write(*,*) 'NaN at SS C[Hpm, Sd, conj[Hpm], conj[Sd]]' 
! ---- Hpm,Su ----
Do i1=1,8
 Do i2=1,6
coup1 = cplHpmSucHpmcSu(i1,i2,i1,i2)
Di_coup1 = dcplHpmSucHpmcSu(i1,i2,i1,i2,iv1)
Dj_coup1 = dcplHpmSucHpmcSu(i1,i2,i1,i2,iv2)
DDcoup1 = ddcplHpmSucHpmcSu(i1,i2,i1,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MHpm2(i1),MSu2(i2),dMHpm2(i1,iv1),dMSu2(i2,iv1)       & 
& ,dMHpm2(i1,iv2),dMSu2(i2,iv2),ddMHpm2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-1._dp)
results2(12)=results2(12) + coeff*temp
results2_ti(12)=results2_ti(12) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(12).eq.results2(12)))  write(*,*) 'NaN at SS C[Hpm, Su, conj[Hpm], conj[Su]]' 
! ---- Sd,Sd ----
Do i1=1,6
 Do i2=1,6
coup1 = cplSdSdcSdcSd(i1,i2,i1,i2)
Di_coup1 = dcplSdSdcSdcSd(i1,i2,i1,i2,iv1)
Dj_coup1 = dcplSdSdcSdcSd(i1,i2,i1,i2,iv2)
DDcoup1 = ddcplSdSdcSdcSd(i1,i2,i1,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MSd2(i1),MSd2(i2),dMSd2(i1,iv1),dMSd2(i2,iv1)         & 
& ,dMSd2(i1,iv2),dMSd2(i2,iv2),ddMSd2(i1,iv1,iv2),ddMSd2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(13)=results2(13) + coeff*temp
results2_ti(13)=results2_ti(13) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(13).eq.results2(13)))  write(*,*) 'NaN at SS C[Sd, Sd, conj[Sd], conj[Sd]]' 
! ---- Sd,Su ----
Do i1=1,6
 Do i2=1,6
coup1 = cplSdSucSdcSu(i1,i2,i1,i2)
Di_coup1 = dcplSdSucSdcSu(i1,i2,i1,i2,iv1)
Dj_coup1 = dcplSdSucSdcSu(i1,i2,i1,i2,iv2)
DDcoup1 = ddcplSdSucSdcSu(i1,i2,i1,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MSd2(i1),MSu2(i2),dMSd2(i1,iv1),dMSu2(i2,iv1)         & 
& ,dMSd2(i1,iv2),dMSu2(i2,iv2),ddMSd2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-1._dp)
results2(14)=results2(14) + coeff*temp
results2_ti(14)=results2_ti(14) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(14).eq.results2(14)))  write(*,*) 'NaN at SS C[Sd, Su, conj[Sd], conj[Su]]' 
! ---- Su,Su ----
Do i1=1,6
 Do i2=1,6
coup1 = cplSuSucSucSu(i1,i2,i1,i2)
Di_coup1 = dcplSuSucSucSu(i1,i2,i1,i2,iv1)
Dj_coup1 = dcplSuSucSucSu(i1,i2,i1,i2,iv2)
DDcoup1 = ddcplSuSucSucSu(i1,i2,i1,i2,iv1,iv2)
Call SecondDerivativeVeff_balls(MSu2(i1),MSu2(i2),dMSu2(i1,iv1),dMSu2(i2,iv1)         & 
& ,dMSu2(i1,iv2),dMSu2(i2,iv2),ddMSu2(i1,iv1,iv2),ddMSu2(i2,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'SS',Q2,temp,temp_ti,temp_tj)
coeff = (-0.5_dp)
results2(15)=results2(15) + coeff*temp
results2_ti(15)=results2_ti(15) + coeff*temp_ti
  End Do
End Do
if (.not.(results2(15).eq.results2(15)))  write(*,*) 'NaN at SS C[Su, Su, conj[Su], conj[Su]]' 

! ----- diagrams of type VS, 2 ------ 

! ---- Sd,VG ----
Do i1=1,6
coup1 = cplSdcSdVGVG(i1,i1)
Di_coup1 = dcplSdcSdVGVG(i1,i1,iv1)
Dj_coup1 = dcplSdcSdVGVG(i1,i1,iv2)
DDcoup1 = ddcplSdcSdVGVG(i1,i1,iv1,iv2)
Call SecondDerivativeVeff_balls(0._dp,MSd2(i1),ZeroC,dMSd2(i1,iv1),ZeroC,dMSd2(i1,iv2)& 
& ,ZeroC,ddMSd2(i1,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'VS',Q2,temp,temp_ti,temp_tj)
coeff = 0._dp
results2(16)=results2(16) + coeff*temp
results2_ti(16)=results2_ti(16) + coeff*temp_ti
End Do
if (.not.(results2(16).eq.results2(16)))  write(*,*) 'NaN at VS C[Sd, VG, VG, conj[Sd]]' 
! ---- Su,VG ----
Do i1=1,6
coup1 = cplSucSuVGVG(i1,i1)
Di_coup1 = dcplSucSuVGVG(i1,i1,iv1)
Dj_coup1 = dcplSucSuVGVG(i1,i1,iv2)
DDcoup1 = ddcplSucSuVGVG(i1,i1,iv1,iv2)
Call SecondDerivativeVeff_balls(0._dp,MSu2(i1),ZeroC,dMSu2(i1,iv1),ZeroC,dMSu2(i1,iv2)& 
& ,ZeroC,ddMSu2(i1,iv1,iv2),coup1,Di_coup1,Dj_coup1,DDcoup1,'VS',Q2,temp,temp_ti,temp_tj)
coeff = 0._dp
results2(17)=results2(17) + coeff*temp
results2_ti(17)=results2_ti(17) + coeff*temp_ti
End Do
if (.not.(results2(17).eq.results2(17)))  write(*,*) 'NaN at VS C[Su, VG, VG, conj[Su]]' 

  result = sum(results1)+sum(results2) ! 2nd deriv. of V
  result_ti = sum(results1_ti)+sum(results2_ti) ! 1st deriv. of V
  pi2L(iv1,iv2) = oo16pi2**2 * Real(result,dp) 
  End Do 
  ti2L(iv1) = oo16pi2**2 * Real(result_ti,dp) 
End Do 
Do iv1=1,8
  Do iv2=1,iv1-1
  pi2L(iv1,iv2) = pi2L(iv2,iv1) 
  End Do 
End Do 
End Subroutine SecondDerivativeEffPot2Loop 





Subroutine FirstDerivativeEffPot2Loop(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,              & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,kont,ti2L)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Integer, Intent(inout):: kont
Real(dp), Intent(out) :: ti2L(8)
Integer :: i,i1,i2,i3,includeGhosts,NrContr 
Integer :: iv1, iv2 
Integer :: NrContr1,NrContr2 !nr of contributing diagrams
Real(dp) :: Q2,colorfactor,coeff,coeffbar
Complex(dp) :: result,temp,tempbar
Complex(dp) :: coup1,coup2,coup1L,coup1R,coup2L,coup2R,coupx,coupxbar
Complex(dp) :: dcoup1,dcoup2,dcoup1L,dcoup1R,dcoup2L,dcoup2R,dcoupx,dcoupxbar
Real(dp) :: gout(112198) 
Real(dp) :: results1(33),results2(17)
Complex(dp) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

Complex(dp) :: dcplAhAhAh(8,8,8,8),dcplAhAhhh(8,8,8,8),dcplAhhhhh(8,8,8,8),dcplAhHpmcHpm(8,8,8,8),   & 
& dcplAhSdcSd(8,6,6,8),dcplAhSucSu(8,6,6,8),dcplhhhhhh(8,8,8,8),dcplhhHpmcHpm(8,8,8,8),  & 
& dcplhhSdcSd(8,6,6,8),dcplhhSucSu(8,6,6,8),dcplHpmSucSd(8,6,6,8),dcplSdcHpmcSu(6,8,6,8),& 
& dcplSdcSdVG(6,6,8),dcplSucSuVG(6,6,8),dcplVGVGVG(8),dcplcChaChaAhL(5,5,8,8),           & 
& dcplcChaChaAhR(5,5,8,8),dcplChiChiAhL(10,10,8,8),dcplChiChiAhR(10,10,8,8),             & 
& dcplcFdFdAhL(3,3,8,8),dcplcFdFdAhR(3,3,8,8),dcplcFuFuAhL(3,3,8,8),dcplcFuFuAhR(3,3,8,8),& 
& dcplChiChacHpmL(10,5,8,8),dcplChiChacHpmR(10,5,8,8),dcplChaFucSdL(5,3,6,8),            & 
& dcplChaFucSdR(5,3,6,8),dcplcChaChahhL(5,5,8,8),dcplcChaChahhR(5,5,8,8),dcplcFdChaSuL(3,5,6,8),& 
& dcplcFdChaSuR(3,5,6,8),dcplChiChihhL(10,10,8,8),dcplChiChihhR(10,10,8,8),              & 
& dcplChiFdcSdL(10,3,6,8),dcplChiFdcSdR(10,3,6,8),dcplChiFucSuL(10,3,6,8),               & 
& dcplChiFucSuR(10,3,6,8),dcplcChaChiHpmL(5,10,8,8),dcplcChaChiHpmR(5,10,8,8),           & 
& dcplcFdChiSdL(3,10,6,8),dcplcFdChiSdR(3,10,6,8),dcplcFuChiSuL(3,10,6,8),               & 
& dcplcFuChiSuR(3,10,6,8),dcplGluFdcSdL(3,6,8),dcplGluFdcSdR(3,6,8),dcplcFdFdhhL(3,3,8,8),& 
& dcplcFdFdhhR(3,3,8,8),dcplcChaFdcSuL(5,3,6,8),dcplcChaFdcSuR(5,3,6,8),dcplcFuFdcHpmL(3,3,8,8),& 
& dcplcFuFdcHpmR(3,3,8,8),dcplGluFucSuL(3,6,8),dcplGluFucSuR(3,6,8),dcplcFuFuhhL(3,3,8,8),& 
& dcplcFuFuhhR(3,3,8,8),dcplcFdFuHpmL(3,3,8,8),dcplcFdFuHpmR(3,3,8,8),dcplcFdGluSdL(3,6,8),& 
& dcplcFdGluSdR(3,6,8),dcplcFuGluSuL(3,6,8),dcplcFuGluSuR(3,6,8),dcplcChacFuSdL(5,3,6,8),& 
& dcplcChacFuSdR(5,3,6,8),dcplcFdFdVGL(3,3,8),dcplcFdFdVGR(3,3,8),dcplcFuFuVGL(3,3,8),   & 
& dcplcFuFuVGR(3,3,8),dcplGluGluVGL(8),dcplGluGluVGR(8)

Complex(dp) :: dcplAhAhAhAh(8,8,8,8,8),dcplAhAhhhhh(8,8,8,8,8),dcplAhAhHpmcHpm(8,8,8,8,8),           & 
& dcplAhAhSdcSd(8,8,6,6,8),dcplAhAhSucSu(8,8,6,6,8),dcplhhhhhhhh(8,8,8,8,8),             & 
& dcplhhhhHpmcHpm(8,8,8,8,8),dcplhhhhSdcSd(8,8,6,6,8),dcplhhhhSucSu(8,8,6,6,8),          & 
& dcplHpmHpmcHpmcHpm(8,8,8,8,8),dcplHpmSdcHpmcSd(8,6,8,6,8),dcplHpmSucHpmcSu(8,6,8,6,8), & 
& dcplSdSdcSdcSd(6,6,6,6,8),dcplSdSucSdcSu(6,6,6,6,8),dcplSuSucSucSu(6,6,6,6,8),         & 
& dcplSdcSdVGVG(6,6,8),dcplSucSuVGVG(6,6,8)

Real(dp) :: MSd(6),MSd2(6),MSu(6),MSu2(6),Mhh(8),Mhh2(8),MAh(8),MAh2(8),MHpm(8),MHpm2(8),         & 
& MChi(10),MChi2(10),MCha(5),MCha2(5),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,          & 
& MVZ,MVZ2,MVWm,MVWm2

Complex(dp) :: dMSd(6,8),dMSd2(6,8),dMSu(6,8),dMSu2(6,8),dMhh(8,8),dMhh2(8,8),dMAh(8,8),             & 
& dMAh2(8,8),dMHpm(8,8),dMHpm2(8,8),dMChi(10,8),dMChi2(10,8),dMCha(5,8),dMCha2(5,8),     & 
& dMFd(3,8),dMFd2(3,8),dMFu(3,8),dMFu2(3,8),dMGlu(1,8),dMGlu2(1,8),dMVZ(1,8),            & 
& dMVZ2(1,8),dMVWm(1,8),dMVWm2(1,8)

!! ------------------------------------------------- 
!! Calculate masses, couplings and their derivatives 
!! ------------------------------------------------- 

Do i1=1,8
Call FirstDerivativeMassesCoups(i1,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,              & 
& kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,gout)

Call GToMassesCoups(gout,MSd,MSd2,MSu,MSu2,Mhh,Mhh2,MAh,MAh2,MHpm,MHpm2,              & 
& MChi,MChi2,MCha,MCha2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MVZ,MVZ2,MVWm,MVWm2,cplAhAhAh,      & 
& cplAhAhhh,cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,         & 
& cplhhSdcSd,cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,        & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,       & 
& cplcChaChahhL,cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,       & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,     & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplcFdFdhhL,cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplGluFucSuL,cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,       & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR,             & 
& cplAhAhAhAh,cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,          & 
& cplhhhhHpmcHpm,cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,            & 
& cplHpmSucHpmcSu,cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)

dMSd(:,i1) = MSd
dMSd2(:,i1) = MSd2
dMSu(:,i1) = MSu
dMSu2(:,i1) = MSu2
dMhh(:,i1) = Mhh
dMhh2(:,i1) = Mhh2
dMAh(:,i1) = MAh
dMAh2(:,i1) = MAh2
dMHpm(:,i1) = MHpm
dMHpm2(:,i1) = MHpm2
dMChi(1:,i1) = MChi
dMChi2(1:,i1) = MChi2
dMCha(:,i1) = MCha
dMCha2(:,i1) = MCha2
dMFd(:,i1) = MFd
dMFd2(:,i1) = MFd2
dMFu(:,i1) = MFu
dMFu2(:,i1) = MFu2
dMGlu(:,i1) = MGlu
dMGlu2(:,i1) = MGlu2
dMVZ(:,i1) = MVZ
dMVZ2(:,i1) = MVZ2
dMVWm(:,i1) = MVWm
dMVWm2(:,i1) = MVWm2
dcplAhAhAh(:,:,:,i1) = cplAhAhAh
dcplAhAhhh(:,:,:,i1) = cplAhAhhh
dcplAhhhhh(:,:,:,i1) = cplAhhhhh
dcplAhHpmcHpm(:,:,:,i1) = cplAhHpmcHpm
dcplAhSdcSd(:,:,:,i1) = cplAhSdcSd
dcplAhSucSu(:,:,:,i1) = cplAhSucSu
dcplhhhhhh(:,:,:,i1) = cplhhhhhh
dcplhhHpmcHpm(:,:,:,i1) = cplhhHpmcHpm
dcplhhSdcSd(:,:,:,i1) = cplhhSdcSd
dcplhhSucSu(:,:,:,i1) = cplhhSucSu
dcplHpmSucSd(:,:,:,i1) = cplHpmSucSd
dcplSdcHpmcSu(:,:,:,i1) = cplSdcHpmcSu
dcplSdcSdVG(:,:,i1) = cplSdcSdVG
dcplSucSuVG(:,:,i1) = cplSucSuVG
dcplVGVGVG(i1) = cplVGVGVG
dcplcChaChaAhL(:,:,:,i1) = cplcChaChaAhL
dcplcChaChaAhR(:,:,:,i1) = cplcChaChaAhR
dcplChiChiAhL(1:,1:,:,i1) = cplChiChiAhL
dcplChiChiAhR(1:,1:,:,i1) = cplChiChiAhR
dcplcFdFdAhL(:,:,:,i1) = cplcFdFdAhL
dcplcFdFdAhR(:,:,:,i1) = cplcFdFdAhR
dcplcFuFuAhL(:,:,:,i1) = cplcFuFuAhL
dcplcFuFuAhR(:,:,:,i1) = cplcFuFuAhR
dcplChiChacHpmL(1:,:,:,i1) = cplChiChacHpmL
dcplChiChacHpmR(1:,:,:,i1) = cplChiChacHpmR
dcplChaFucSdL(:,:,:,i1) = cplChaFucSdL
dcplChaFucSdR(:,:,:,i1) = cplChaFucSdR
dcplcChaChahhL(:,:,:,i1) = cplcChaChahhL
dcplcChaChahhR(:,:,:,i1) = cplcChaChahhR
dcplcFdChaSuL(:,:,:,i1) = cplcFdChaSuL
dcplcFdChaSuR(:,:,:,i1) = cplcFdChaSuR
dcplChiChihhL(1:,1:,:,i1) = cplChiChihhL
dcplChiChihhR(1:,1:,:,i1) = cplChiChihhR
dcplChiFdcSdL(1:,:,:,i1) = cplChiFdcSdL
dcplChiFdcSdR(1:,:,:,i1) = cplChiFdcSdR
dcplChiFucSuL(1:,:,:,i1) = cplChiFucSuL
dcplChiFucSuR(1:,:,:,i1) = cplChiFucSuR
dcplcChaChiHpmL(:,1:,:,i1) = cplcChaChiHpmL
dcplcChaChiHpmR(:,1:,:,i1) = cplcChaChiHpmR
dcplcFdChiSdL(:,1:,:,i1) = cplcFdChiSdL
dcplcFdChiSdR(:,1:,:,i1) = cplcFdChiSdR
dcplcFuChiSuL(:,1:,:,i1) = cplcFuChiSuL
dcplcFuChiSuR(:,1:,:,i1) = cplcFuChiSuR
dcplGluFdcSdL(:,:,i1) = cplGluFdcSdL
dcplGluFdcSdR(:,:,i1) = cplGluFdcSdR
dcplcFdFdhhL(:,:,:,i1) = cplcFdFdhhL
dcplcFdFdhhR(:,:,:,i1) = cplcFdFdhhR
dcplcChaFdcSuL(:,:,:,i1) = cplcChaFdcSuL
dcplcChaFdcSuR(:,:,:,i1) = cplcChaFdcSuR
dcplcFuFdcHpmL(:,:,:,i1) = cplcFuFdcHpmL
dcplcFuFdcHpmR(:,:,:,i1) = cplcFuFdcHpmR
dcplGluFucSuL(:,:,i1) = cplGluFucSuL
dcplGluFucSuR(:,:,i1) = cplGluFucSuR
dcplcFuFuhhL(:,:,:,i1) = cplcFuFuhhL
dcplcFuFuhhR(:,:,:,i1) = cplcFuFuhhR
dcplcFdFuHpmL(:,:,:,i1) = cplcFdFuHpmL
dcplcFdFuHpmR(:,:,:,i1) = cplcFdFuHpmR
dcplcFdGluSdL(:,:,i1) = cplcFdGluSdL
dcplcFdGluSdR(:,:,i1) = cplcFdGluSdR
dcplcFuGluSuL(:,:,i1) = cplcFuGluSuL
dcplcFuGluSuR(:,:,i1) = cplcFuGluSuR
dcplcChacFuSdL(:,:,:,i1) = cplcChacFuSdL
dcplcChacFuSdR(:,:,:,i1) = cplcChacFuSdR
dcplcFdFdVGL(:,:,i1) = cplcFdFdVGL
dcplcFdFdVGR(:,:,i1) = cplcFdFdVGR
dcplcFuFuVGL(:,:,i1) = cplcFuFuVGL
dcplcFuFuVGR(:,:,i1) = cplcFuFuVGR
dcplGluGluVGL(i1) = cplGluGluVGL
dcplGluGluVGR(i1) = cplGluGluVGR
dcplAhAhAhAh(:,:,:,:,i1) = cplAhAhAhAh
dcplAhAhhhhh(:,:,:,:,i1) = cplAhAhhhhh
dcplAhAhHpmcHpm(:,:,:,:,i1) = cplAhAhHpmcHpm
dcplAhAhSdcSd(:,:,:,:,i1) = cplAhAhSdcSd
dcplAhAhSucSu(:,:,:,:,i1) = cplAhAhSucSu
dcplhhhhhhhh(:,:,:,:,i1) = cplhhhhhhhh
dcplhhhhHpmcHpm(:,:,:,:,i1) = cplhhhhHpmcHpm
dcplhhhhSdcSd(:,:,:,:,i1) = cplhhhhSdcSd
dcplhhhhSucSu(:,:,:,:,i1) = cplhhhhSucSu
dcplHpmHpmcHpmcHpm(:,:,:,:,i1) = cplHpmHpmcHpmcHpm
dcplHpmSdcHpmcSd(:,:,:,:,i1) = cplHpmSdcHpmcSd
dcplHpmSucHpmcSu(:,:,:,:,i1) = cplHpmSucHpmcSu
dcplSdSdcSdcSd(:,:,:,:,i1) = cplSdSdcSdcSd
dcplSdSucSdcSu(:,:,:,:,i1) = cplSdSucSdcSu
dcplSuSucSucSu(:,:,:,:,i1) = cplSuSucSucSu
dcplSdcSdVGVG(:,:,i1) = cplSdcSdVGVG
dcplSucSuVGVG(:,:,i1) = cplSucSuVGVG
End Do 
 
Call TreeMassesEffPot(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,               & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,.True.,kont)

Call CouplingsForEffPot3(lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZA,ZH,Ye,Te,               & 
& ZP,Yd,Td,ZD,Yu,Tu,ZU,g3,ZER,ZEL,UV,ZDL,ZDR,ZUL,ZUR,pG,cplAhAhAh,cplAhAhhh,             & 
& cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,cplhhSdcSd,        & 
& cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplcChaChaAhL,     & 
& cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,cplcFuFuAhL,           & 
& cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,cplcChaChahhL,     & 
& cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,cplChiFdcSdL,        & 
& cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,cplcFdChiSdL,     & 
& cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,cplcFdFdhhL,          & 
& cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplGluFucSuL,      & 
& cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdGluSdL,           & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,cplcFdFdVGL,        & 
& cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR)

Call CouplingsForEffPot4(lam,Yv,kap,ZA,ZH,Ye,ZP,Yd,ZD,Yu,ZU,g3,cplAhAhAhAh,           & 
& cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,cplhhhhHpmcHpm,       & 
& cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,cplHpmSucHpmcSu,           & 
& cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)



!! ------------------------------------------------- 
!! Calculate derivative of effective potential      
!! ------------------------------------------------- 



Q2 = getRenormalizationScale()
ti2L = 0._dp
result = ZeroC
results1 = ZeroC
results2 = ZeroC
Do iv1=1,8
! ----- Topology1 (sunrise): diagrams w. 3 Particles and 2 Vertices

! ----- diagrams of type SSS, 11 ------ 
! ---- Ah,Ah,Ah ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhAhAh(i1,i2,i3)
coup2 = cplAhAhAh(i1,i2,i3)
dcoup1 = dcplAhAhAh(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MAh2(i1),MAh2(i2),MAh2(i3),dMAh2(i1,iv1)             & 
& ,dMAh2(i2,iv1),dMAh2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 1._dp/12._dp
colorfactor = 1
results1(1)=results1(1) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Ah, Ah]' 
    End Do
  End Do
End Do
! ---- Ah,Ah,hh ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhAhhh(i1,i2,i3)
coup2 = cplAhAhhh(i1,i2,i3)
dcoup1 = dcplAhAhhh(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MAh2(i1),MAh2(i2),Mhh2(i3),dMAh2(i1,iv1)             & 
& ,dMAh2(i2,iv1),dMhh2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.25_dp
colorfactor = 1
results1(2)=results1(2) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Ah, hh]' 
    End Do
  End Do
End Do
! ---- Ah,hh,hh ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhhhhh(i1,i2,i3)
coup2 = cplAhhhhh(i1,i2,i3)
dcoup1 = dcplAhhhhh(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MAh2(i1),Mhh2(i2),Mhh2(i3),dMAh2(i1,iv1)             & 
& ,dMhh2(i2,iv1),dMhh2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.25_dp
colorfactor = 1
results1(3)=results1(3) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, hh, hh]' 
    End Do
  End Do
End Do
! ---- Ah,Hpm,conj[Hpm] ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplAhHpmcHpm(i1,i2,i3)
coup2 = cplAhHpmcHpm(i1,i3,i2)
dcoup1 = dcplAhHpmcHpm(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MAh2(i1),MHpm2(i2),MHpm2(i3),dMAh2(i1,iv1)           & 
& ,dMHpm2(i2,iv1),dMHpm2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 1
results1(4)=results1(4) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Hpm, conj[Hpm]]' 
    End Do
  End Do
End Do
! ---- Ah,Sd,conj[Sd] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplAhSdcSd(i1,i2,i3)
coup2 = cplAhSdcSd(i1,i3,i2)
dcoup1 = dcplAhSdcSd(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MAh2(i1),MSd2(i2),MSd2(i3),dMAh2(i1,iv1)             & 
& ,dMSd2(i2,iv1),dMSd2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 3
results1(5)=results1(5) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Sd, conj[Sd]]' 
    End Do
  End Do
End Do
! ---- Ah,Su,conj[Su] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplAhSucSu(i1,i2,i3)
coup2 = cplAhSucSu(i1,i3,i2)
dcoup1 = dcplAhSucSu(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MAh2(i1),MSu2(i2),MSu2(i3),dMAh2(i1,iv1)             & 
& ,dMSu2(i2,iv1),dMSu2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 3
results1(6)=results1(6) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Ah, Su, conj[Su]]' 
    End Do
  End Do
End Do
! ---- hh,hh,hh ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplhhhhhh(i1,i2,i3)
coup2 = cplhhhhhh(i1,i2,i3)
dcoup1 = dcplhhhhhh(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(Mhh2(i1),Mhh2(i2),Mhh2(i3),dMhh2(i1,iv1)             & 
& ,dMhh2(i2,iv1),dMhh2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 1._dp/12._dp
colorfactor = 1
results1(7)=results1(7) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, hh, hh]' 
    End Do
  End Do
End Do
! ---- hh,Hpm,conj[Hpm] ----
Do i1=1,8
 Do i2=1,8
    Do i3=1,8
coup1 = cplhhHpmcHpm(i1,i2,i3)
coup2 = cplhhHpmcHpm(i1,i3,i2)
dcoup1 = dcplhhHpmcHpm(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(Mhh2(i1),MHpm2(i2),MHpm2(i3),dMhh2(i1,iv1)           & 
& ,dMHpm2(i2,iv1),dMHpm2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 1
results1(8)=results1(8) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, Hpm, conj[Hpm]]' 
    End Do
  End Do
End Do
! ---- hh,Sd,conj[Sd] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplhhSdcSd(i1,i2,i3)
coup2 = cplhhSdcSd(i1,i3,i2)
dcoup1 = dcplhhSdcSd(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(Mhh2(i1),MSd2(i2),MSd2(i3),dMhh2(i1,iv1)             & 
& ,dMSd2(i2,iv1),dMSd2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 3
results1(9)=results1(9) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, Sd, conj[Sd]]' 
    End Do
  End Do
End Do
! ---- hh,Su,conj[Su] ----
Do i1=1,8
 Do i2=1,6
    Do i3=1,6
coup1 = cplhhSucSu(i1,i2,i3)
coup2 = cplhhSucSu(i1,i3,i2)
dcoup1 = dcplhhSucSu(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(Mhh2(i1),MSu2(i2),MSu2(i3),dMhh2(i1,iv1)             & 
& ,dMSu2(i2,iv1),dMSu2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 3
results1(10)=results1(10) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[hh, Su, conj[Su]]' 
    End Do
  End Do
End Do
! ---- Sd,conj[Hpm],conj[Su] ----
Do i1=1,6
 Do i2=1,8
    Do i3=1,6
coup1 = cplSdcHpmcSu(i1,i2,i3)
coup2 = cplHpmSucSd(i2,i3,i1)
dcoup1 = dcplSdcHpmcSu(i1,i2,i3,iv1)
coupx=abs(coup1)**2 
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MSd2(i1),MHpm2(i2),MSu2(i3),dMSd2(i1,iv1)            & 
& ,dMHpm2(i2,iv1),dMSu2(i3,iv1),coupx,dcoupx,'SSS   ',Q2,temp)
coeff = 1._dp
colorfactor = 3
results1(11)=results1(11) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSS C[Sd, conj[Hpm], conj[Su]]' 
    End Do
  End Do
End Do
! ----- diagrams of type FFS, 16 ------ 
! ---- Ah,Cha,bar[Cha] ----
Do i1=1,8
 Do i2=1,5
    Do i3=1,5
coup1L = cplcChaChaAhL(i3,i2,i1)
coup1R = cplcChaChaAhR(i3,i2,i1)
coup2L = cplcChaChaAhL(i2,i3,i1)
coup2R = cplcChaChaAhR(i2,i3,i1)
dcoup1L = dcplcChaChaAhL(i3,i2,i1,iv1)
dcoup1R = dcplcChaChaAhR(i3,i2,i1,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MCha2(i3),MCha2(i2),MAh2(i1),dMCha2(i3,iv1)          & 
& ,dMCha2(i2,iv1),dMAh2(i1,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MCha2(i3),MCha2(i2),MAh2(i1),dMCha2(i3,iv1)          & 
& ,dMCha2(i2,iv1),dMAh2(i1,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 1
results1(12)=results1(12) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Cha, bar[Cha]]' 
coeffbar = 0.5_dp
results1(12)=results1(12) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Ah,Chi,Chi ----
Do i1=1,8
 Do i2=1,10
    Do i3=1,10
coup1L = cplChiChiAhL(i2,i3,i1)
coup1R = cplChiChiAhR(i2,i3,i1)
coup2L = cplChiChiAhL(i2,i3,i1)
coup2R = cplChiChiAhR(i2,i3,i1)
dcoup1L = dcplChiChiAhL(i2,i3,i1,iv1)
dcoup1R = dcplChiChiAhR(i2,i3,i1,iv1)
coupx=(abs(coup1L)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L) 
coupxbar=Real(coup1L**2,dp) 
dcoupxbar=Real(2*dcoup1L*coup1L,dp) 
Call FirstDerivativeVeff_sunrise(MChi2(i3),MChi2(i2),MAh2(i1),dMChi2(i3,iv1)          & 
& ,dMChi2(i2,iv1),dMAh2(i1,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MChi2(i3),MChi2(i2),MAh2(i1),dMChi2(i3,iv1)          & 
& ,dMChi2(i2,iv1),dMAh2(i1,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 1
results1(13)=results1(13) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Chi, Chi]' 
coeffbar = 0.5_dp
results1(13)=results1(13) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Ah,Fd,bar[Fd] ----
Do i1=1,8
 Do i2=1,3
    Do i3=1,3
coup1L = cplcFdFdAhL(i3,i2,i1)
coup1R = cplcFdFdAhR(i3,i2,i1)
coup2L = cplcFdFdAhL(i2,i3,i1)
coup2R = cplcFdFdAhR(i2,i3,i1)
dcoup1L = dcplcFdFdAhL(i3,i2,i1,iv1)
dcoup1R = dcplcFdFdAhR(i3,i2,i1,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFd2(i2),MAh2(i1),dMFd2(i3,iv1)             & 
& ,dMFd2(i2,iv1),dMAh2(i1,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFd2(i2),MAh2(i1),dMFd2(i3,iv1)             & 
& ,dMFd2(i2,iv1),dMAh2(i1,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 3
results1(14)=results1(14) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Fd, bar[Fd]]' 
coeffbar = 0.5_dp
results1(14)=results1(14) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Ah,Fu,bar[Fu] ----
Do i1=1,8
 Do i2=1,3
    Do i3=1,3
coup1L = cplcFuFuAhL(i3,i2,i1)
coup1R = cplcFuFuAhR(i3,i2,i1)
coup2L = cplcFuFuAhL(i2,i3,i1)
coup2R = cplcFuFuAhR(i2,i3,i1)
dcoup1L = dcplcFuFuAhL(i3,i2,i1,iv1)
dcoup1R = dcplcFuFuAhR(i3,i2,i1,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFu2(i3),MFu2(i2),MAh2(i1),dMFu2(i3,iv1)             & 
& ,dMFu2(i2,iv1),dMAh2(i1,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFu2(i3),MFu2(i2),MAh2(i1),dMFu2(i3,iv1)             & 
& ,dMFu2(i2,iv1),dMAh2(i1,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 3
results1(15)=results1(15) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Ah, Fu, bar[Fu]]' 
coeffbar = 0.5_dp
results1(15)=results1(15) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Cha,hh,bar[Cha] ----
Do i1=1,5
 Do i2=1,8
    Do i3=1,5
coup1L = cplcChaChahhL(i3,i1,i2)
coup1R = cplcChaChahhR(i3,i1,i2)
coup2L = cplcChaChahhL(i1,i3,i2)
coup2R = cplcChaChahhR(i1,i3,i2)
dcoup1L = dcplcChaChahhL(i3,i1,i2,iv1)
dcoup1R = dcplcChaChahhR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MCha2(i3),MCha2(i1),Mhh2(i2),dMCha2(i3,iv1)          & 
& ,dMCha2(i1,iv1),dMhh2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MCha2(i3),MCha2(i1),Mhh2(i2),dMCha2(i3,iv1)          & 
& ,dMCha2(i1,iv1),dMhh2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 1
results1(16)=results1(16) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Cha, hh, bar[Cha]]' 
coeffbar = 0.5_dp
results1(16)=results1(16) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Chi,Chi,hh ----
Do i1=1,10
 Do i2=1,10
    Do i3=1,8
coup1L = cplChiChihhL(i1,i2,i3)
coup1R = cplChiChihhR(i1,i2,i3)
coup2L = cplChiChihhL(i1,i2,i3)
coup2R = cplChiChihhR(i1,i2,i3)
dcoup1L = dcplChiChihhL(i1,i2,i3,iv1)
dcoup1R = dcplChiChihhR(i1,i2,i3,iv1)
coupx=(abs(coup1L)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L) 
coupxbar=Real(coup1L**2,dp) 
dcoupxbar=Real(2*dcoup1L*coup1L,dp) 
Call FirstDerivativeVeff_sunrise(MChi2(i2),MChi2(i1),Mhh2(i3),dMChi2(i2,iv1)          & 
& ,dMChi2(i1,iv1),dMhh2(i3,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MChi2(i2),MChi2(i1),Mhh2(i3),dMChi2(i2,iv1)          & 
& ,dMChi2(i1,iv1),dMhh2(i3,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 1
results1(17)=results1(17) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Chi, hh]' 
coeffbar = 0.5_dp
results1(17)=results1(17) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Chi,Hpm,bar[Cha] ----
Do i1=1,10
 Do i2=1,8
    Do i3=1,5
coup1L = cplcChaChiHpmL(i3,i1,i2)
coup1R = cplcChaChiHpmR(i3,i1,i2)
coup2L = cplChiChacHpmL(i1,i3,i2)
coup2R = cplChiChacHpmR(i1,i3,i2)
dcoup1L = dcplcChaChiHpmL(i3,i1,i2,iv1)
dcoup1R = dcplcChaChiHpmR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MCha2(i3),MChi2(i1),MHpm2(i2),dMCha2(i3,iv1)         & 
& ,dMChi2(i1,iv1),dMHpm2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MCha2(i3),MChi2(i1),MHpm2(i2),dMCha2(i3,iv1)         & 
& ,dMChi2(i1,iv1),dMHpm2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 1
results1(18)=results1(18) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Hpm, bar[Cha]]' 
coeffbar = 1._dp
results1(18)=results1(18) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Chi,Sd,bar[Fd] ----
Do i1=1,10
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFdChiSdL(i3,i1,i2)
coup1R = cplcFdChiSdR(i3,i1,i2)
coup2L = cplChiFdcSdL(i1,i3,i2)
coup2R = cplChiFdcSdR(i1,i3,i2)
dcoup1L = dcplcFdChiSdL(i3,i1,i2,iv1)
dcoup1R = dcplcFdChiSdR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFd2(i3),MChi2(i1),MSd2(i2),dMFd2(i3,iv1)            & 
& ,dMChi2(i1,iv1),dMSd2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFd2(i3),MChi2(i1),MSd2(i2),dMFd2(i3,iv1)            & 
& ,dMChi2(i1,iv1),dMSd2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 3
results1(19)=results1(19) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Sd, bar[Fd]]' 
coeffbar = 1._dp
results1(19)=results1(19) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Chi,Su,bar[Fu] ----
Do i1=1,10
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFuChiSuL(i3,i1,i2)
coup1R = cplcFuChiSuR(i3,i1,i2)
coup2L = cplChiFucSuL(i1,i3,i2)
coup2R = cplChiFucSuR(i1,i3,i2)
dcoup1L = dcplcFuChiSuL(i3,i1,i2,iv1)
dcoup1R = dcplcFuChiSuR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFu2(i3),MChi2(i1),MSu2(i2),dMFu2(i3,iv1)            & 
& ,dMChi2(i1,iv1),dMSu2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFu2(i3),MChi2(i1),MSu2(i2),dMFu2(i3,iv1)            & 
& ,dMChi2(i1,iv1),dMSu2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 3
results1(20)=results1(20) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Chi, Su, bar[Fu]]' 
coeffbar = 1._dp
results1(20)=results1(20) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Fd,hh,bar[Fd] ----
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFdFdhhL(i3,i1,i2)
coup1R = cplcFdFdhhR(i3,i1,i2)
coup2L = cplcFdFdhhL(i1,i3,i2)
coup2R = cplcFdFdhhR(i1,i3,i2)
dcoup1L = dcplcFdFdhhL(i3,i1,i2,iv1)
dcoup1R = dcplcFdFdhhR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),Mhh2(i2),dMFd2(i3,iv1)             & 
& ,dMFd2(i1,iv1),dMhh2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),Mhh2(i2),dMFd2(i3,iv1)             & 
& ,dMFd2(i1,iv1),dMhh2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 3
results1(21)=results1(21) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fd, hh, bar[Fd]]' 
coeffbar = 0.5_dp
results1(21)=results1(21) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Fd,bar[Cha],conj[Su] ----
Do i1=1,3
 Do i2=1,5
    Do i3=1,6
coup1L = cplcChaFdcSuL(i2,i1,i3)
coup1R = cplcChaFdcSuR(i2,i1,i3)
coup2L = cplcFdChaSuL(i1,i2,i3)
coup2R = cplcFdChaSuR(i1,i2,i3)
dcoup1L = dcplcChaFdcSuL(i2,i1,i3,iv1)
dcoup1R = dcplcChaFdcSuR(i2,i1,i3,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MCha2(i2),MFd2(i1),MSu2(i3),dMCha2(i2,iv1)           & 
& ,dMFd2(i1,iv1),dMSu2(i3,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MCha2(i2),MFd2(i1),MSu2(i3),dMCha2(i2,iv1)           & 
& ,dMFd2(i1,iv1),dMSu2(i3,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 3
results1(22)=results1(22) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fd, bar[Cha], conj[Su]]' 
coeffbar = 1._dp
results1(22)=results1(22) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Fu,hh,bar[Fu] ----
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFuFuhhL(i3,i1,i2)
coup1R = cplcFuFuhhR(i3,i1,i2)
coup2L = cplcFuFuhhL(i1,i3,i2)
coup2R = cplcFuFuhhR(i1,i3,i2)
dcoup1L = dcplcFuFuhhL(i3,i1,i2,iv1)
dcoup1R = dcplcFuFuhhR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),Mhh2(i2),dMFu2(i3,iv1)             & 
& ,dMFu2(i1,iv1),dMhh2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),Mhh2(i2),dMFu2(i3,iv1)             & 
& ,dMFu2(i1,iv1),dMhh2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 0.5_dp
colorfactor = 3
results1(23)=results1(23) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fu, hh, bar[Fu]]' 
coeffbar = 0.5_dp
results1(23)=results1(23) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Fu,Hpm,bar[Fd] ----
Do i1=1,3
 Do i2=1,8
    Do i3=1,3
coup1L = cplcFdFuHpmL(i3,i1,i2)
coup1R = cplcFdFuHpmR(i3,i1,i2)
coup2L = cplcFuFdcHpmL(i1,i3,i2)
coup2R = cplcFuFdcHpmR(i1,i3,i2)
dcoup1L = dcplcFdFuHpmL(i3,i1,i2,iv1)
dcoup1R = dcplcFdFuHpmR(i3,i1,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFu2(i1),MHpm2(i2),dMFd2(i3,iv1)            & 
& ,dMFu2(i1,iv1),dMHpm2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFu2(i1),MHpm2(i2),dMFd2(i3,iv1)            & 
& ,dMFu2(i1,iv1),dMHpm2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 3
results1(24)=results1(24) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Fu, Hpm, bar[Fd]]' 
coeffbar = 1._dp
results1(24)=results1(24) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ---- Glu,Sd,bar[Fd] ----
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFdGluSdL(i3,i2)
coup1R = cplcFdGluSdR(i3,i2)
coup2L = cplGluFdcSdL(i3,i2)
coup2R = cplGluFdcSdR(i3,i2)
dcoup1L = dcplcFdGluSdL(i3,i2,iv1)
dcoup1R = dcplcFdGluSdR(i3,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFd2(i3),MGlu2,MSd2(i2),dMFd2(i3,iv1),dMGlu2(1,iv1)  & 
& ,dMSd2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFd2(i3),MGlu2,MSd2(i2),dMFd2(i3,iv1),dMGlu2(1,iv1)  & 
& ,dMSd2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 4
results1(25)=results1(25) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Glu, Sd, bar[Fd]]' 
coeffbar = 1._dp
results1(25)=results1(25) + coeffbar*colorfactor*tempbar
    End Do
  End Do
! ---- Glu,Su,bar[Fu] ----
 Do i2=1,6
    Do i3=1,3
coup1L = cplcFuGluSuL(i3,i2)
coup1R = cplcFuGluSuR(i3,i2)
coup2L = cplGluFucSuL(i3,i2)
coup2R = cplGluFucSuR(i3,i2)
dcoup1L = dcplcFuGluSuL(i3,i2,iv1)
dcoup1R = dcplcFuGluSuR(i3,i2,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFu2(i3),MGlu2,MSu2(i2),dMFu2(i3,iv1),dMGlu2(1,iv1)  & 
& ,dMSu2(i2,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFu2(i3),MGlu2,MSu2(i2),dMFu2(i3,iv1),dMGlu2(1,iv1)  & 
& ,dMSu2(i2,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 4
results1(26)=results1(26) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Glu, Su, bar[Fu]]' 
coeffbar = 1._dp
results1(26)=results1(26) + coeffbar*colorfactor*tempbar
    End Do
  End Do
! ---- Sd,bar[Cha],bar[Fu] ----
Do i1=1,6
 Do i2=1,5
    Do i3=1,3
coup1L = cplcChacFuSdL(i2,i3,i1)
coup1R = cplcChacFuSdR(i2,i3,i1)
coup2L = cplChaFucSdL(i2,i3,i1)
coup2R = cplChaFucSdR(i2,i3,i1)
dcoup1L = dcplcChacFuSdL(i2,i3,i1,iv1)
dcoup1R = dcplcChacFuSdR(i2,i3,i1,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFu2(i3),MCha2(i2),MSd2(i1),dMFu2(i3,iv1)            & 
& ,dMCha2(i2,iv1),dMSd2(i1,iv1),coupx,dcoupx,'FFS   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFu2(i3),MCha2(i2),MSd2(i1),dMFu2(i3,iv1)            & 
& ,dMCha2(i2,iv1),dMSd2(i1,iv1),coupxbar,dcoupxbar,'FFSbar',Q2,tempbar)
coeff = 1._dp
colorfactor = 3
results1(27)=results1(27) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFS C[Sd, bar[Cha], bar[Fu]]' 
coeffbar = 1._dp
results1(27)=results1(27) + coeffbar*colorfactor*tempbar
    End Do
  End Do
End Do
! ----- diagrams of type SSV, 2 ------ 
! ---- Sd,VG,conj[Sd] ----
Do i1=1,6
    Do i3=1,6
coup1 = cplSdcSdVG(i1,i3)
coup2 = cplSdcSdVG(i3,i1)
dcoup1 = dcplSdcSdVG(i1,i3,iv1)
coupx=abs(coup1)**2
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MSd2(i3),MSd2(i1),0._dp,dMSd2(i3,iv1),dMSd2(i1,iv1)  & 
& ,ZeroC,coupx,dcoupx,'SSV   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 4
results1(28)=results1(28) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSV C[Sd, VG, conj[Sd]]' 
    End Do
End Do
! ---- Su,VG,conj[Su] ----
Do i1=1,6
    Do i3=1,6
coup1 = cplSucSuVG(i1,i3)
coup2 = cplSucSuVG(i3,i1)
dcoup1 = dcplSucSuVG(i1,i3,iv1)
coupx=abs(coup1)**2
dcoupx=dcoup1*conjg(coup1)+coup1*conjg(dcoup1) 
Call FirstDerivativeVeff_sunrise(MSu2(i3),MSu2(i1),0._dp,dMSu2(i3,iv1),dMSu2(i1,iv1)  & 
& ,ZeroC,coupx,dcoupx,'SSV   ',Q2,temp)
coeff = 0.5_dp
colorfactor = 4
results1(29)=results1(29) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SSV C[Su, VG, conj[Su]]' 
    End Do
End Do
! ----- diagrams of type FFV, 3 ------ 
! ---- Fd,VG,bar[Fd] ----
Do i1=1,3
    Do i3=1,3
coup1L = cplcFdFdVGL(i3,i1)
coup1R = cplcFdFdVGR(i3,i1)
coup2L = cplcFdFdVGL(i1,i3)
coup2R = cplcFdFdVGR(i1,i3)
dcoup1L = dcplcFdFdVGL(i3,i1,iv1)
dcoup1R = dcplcFdFdVGR(i3,i1,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=-2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=-2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),0._dp,dMFd2(i3,iv1),dMFd2(i1,iv1)  & 
& ,ZeroC,coupx,dcoupx,'FFV   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFd2(i3),MFd2(i1),0._dp,dMFd2(i3,iv1),dMFd2(i1,iv1)  & 
& ,ZeroC,coupxbar,dcoupxbar,'FFVbar',Q2,temp)
coeff = 0.5_dp
colorfactor = 4
results1(30)=results1(30) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFV C[Fd, VG, bar[Fd]]' 
coeffbar = 0.5_dp
results1(30)=results1(30) + coeffbar*colorfactor*tempbar
    End Do
End Do
! ---- Fu,VG,bar[Fu] ----
Do i1=1,3
    Do i3=1,3
coup1L = cplcFuFuVGL(i3,i1)
coup1R = cplcFuFuVGR(i3,i1)
coup2L = cplcFuFuVGL(i1,i3)
coup2R = cplcFuFuVGR(i1,i3)
dcoup1L = dcplcFuFuVGL(i3,i1,iv1)
dcoup1R = dcplcFuFuVGR(i3,i1,iv1)
coupx=(abs(coup1L)**2+abs(coup1R)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L)+dcoup1R*conjg(coup1R)+coup1R*conjg(dcoup1R) 
coupxbar=-2*Real(coup1L*conjg(coup1R),dp) 
dcoupxbar=-2*Real(dcoup1L*conjg(coup1R)+coup1L*conjg(dcoup1R),dp) 
Call FirstDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),0._dp,dMFu2(i3,iv1),dMFu2(i1,iv1)  & 
& ,ZeroC,coupx,dcoupx,'FFV   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MFu2(i3),MFu2(i1),0._dp,dMFu2(i3,iv1),dMFu2(i1,iv1)  & 
& ,ZeroC,coupxbar,dcoupxbar,'FFVbar',Q2,temp)
coeff = 0.5_dp
colorfactor = 4
results1(31)=results1(31) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFV C[Fu, VG, bar[Fu]]' 
coeffbar = 0.5_dp
results1(31)=results1(31) + coeffbar*colorfactor*tempbar
    End Do
End Do
! ---- Glu,Glu,VG ----
coup1L = cplGluGluVGL
coup1R = cplGluGluVGR
coup2L = cplGluGluVGL
coup2R = cplGluGluVGR
dcoup1L = dcplGluGluVGL(iv1)
dcoup1R = dcplGluGluVGR(iv1)
coupx=(abs(coup1L)**2) 
dcoupx=dcoup1L*conjg(coup1L)+coup1L*conjg(dcoup1L) 
coupxbar=Real(coup1L**2,dp) 
dcoupxbar=Real(2*dcoup1L*coup1L,dp) 
Call FirstDerivativeVeff_sunrise(MGlu2,MGlu2,0._dp,dMGlu2(1,iv1),dMGlu2(1,iv1)        & 
& ,ZeroC,coupx,dcoupx,'FFV   ',Q2,temp)
Call FirstDerivativeVeff_sunrise(MGlu2,MGlu2,0._dp,dMGlu2(1,iv1),dMGlu2(1,iv1)        & 
& ,ZeroC,coupxbar,dcoupxbar,'FFVbar',Q2,temp)
coeff = 0.5_dp
colorfactor = 24
results1(32)=results1(32) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at FFV C[Glu, Glu, VG]' 
coeffbar = 0.5_dp
results1(32)=results1(32) + coeffbar*colorfactor*tempbar
! ----- diagrams of type VVV, 1 ------ 
! ---- VG,VG,VG ----
coup1 = cplVGVGVG
coup2 = cplVGVGVG
dcoup1 = dcplVGVGVG(iv1)
coeff = 0.0000
colorfactor = 24
results1(33)=results1(33) + coeff*colorfactor*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at VVV C[VG, VG, VG]' 
! ----- Topology2: diagrams w. 2 Particles and 1 Vertex

! ----- diagrams of type SS, 15 ------ 
! ---- Ah,Ah ----
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhAhAh(i1,i1,i2,i2)
dcoup1 = dcplAhAhAhAh(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(MAh2(i1),MAh2(i2),dMAh2(i1,iv1),dMAh2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-1._dp/8._dp)
results2(1)=results2(1) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Ah, Ah]' 
  End Do
End Do
! ---- Ah,hh ----
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhhhhh(i1,i1,i2,i2)
dcoup1 = dcplAhAhhhhh(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(MAh2(i1),Mhh2(i2),dMAh2(i1,iv1),dMhh2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.25_dp)
results2(2)=results2(2) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, hh, hh]' 
  End Do
End Do
! ---- Ah,Hpm ----
Do i1=1,8
 Do i2=1,8
coup1 = cplAhAhHpmcHpm(i1,i1,i2,i2)
dcoup1 = dcplAhAhHpmcHpm(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(MAh2(i1),MHpm2(i2),dMAh2(i1,iv1),dMHpm2(i2,iv1)        & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(3)=results2(3) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Hpm, conj[Hpm]]' 
  End Do
End Do
! ---- Ah,Sd ----
Do i1=1,8
 Do i2=1,6
coup1 = cplAhAhSdcSd(i1,i1,i2,i2)
dcoup1 = dcplAhAhSdcSd(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(MAh2(i1),MSd2(i2),dMAh2(i1,iv1),dMSd2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(4)=results2(4) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Sd, conj[Sd]]' 
  End Do
End Do
! ---- Ah,Su ----
Do i1=1,8
 Do i2=1,6
coup1 = cplAhAhSucSu(i1,i1,i2,i2)
dcoup1 = dcplAhAhSucSu(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(MAh2(i1),MSu2(i2),dMAh2(i1,iv1),dMSu2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(5)=results2(5) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Ah, Ah, Su, conj[Su]]' 
  End Do
End Do
! ---- hh,hh ----
Do i1=1,8
 Do i2=1,8
coup1 = cplhhhhhhhh(i1,i1,i2,i2)
dcoup1 = dcplhhhhhhhh(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(Mhh2(i1),Mhh2(i2),dMhh2(i1,iv1),dMhh2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-1._dp/8._dp)
results2(6)=results2(6) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, hh, hh]' 
  End Do
End Do
! ---- hh,Hpm ----
Do i1=1,8
 Do i2=1,8
coup1 = cplhhhhHpmcHpm(i1,i1,i2,i2)
dcoup1 = dcplhhhhHpmcHpm(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(Mhh2(i1),MHpm2(i2),dMhh2(i1,iv1),dMHpm2(i2,iv1)        & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(7)=results2(7) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, Hpm, conj[Hpm]]' 
  End Do
End Do
! ---- hh,Sd ----
Do i1=1,8
 Do i2=1,6
coup1 = cplhhhhSdcSd(i1,i1,i2,i2)
dcoup1 = dcplhhhhSdcSd(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(Mhh2(i1),MSd2(i2),dMhh2(i1,iv1),dMSd2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(8)=results2(8) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, Sd, conj[Sd]]' 
  End Do
End Do
! ---- hh,Su ----
Do i1=1,8
 Do i2=1,6
coup1 = cplhhhhSucSu(i1,i1,i2,i2)
dcoup1 = dcplhhhhSucSu(i1,i1,i2,i2,iv1)
Call FirstDerivativeVeff_balls(Mhh2(i1),MSu2(i2),dMhh2(i1,iv1),dMSu2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(9)=results2(9) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[hh, hh, Su, conj[Su]]' 
  End Do
End Do
! ---- Hpm,Hpm ----
Do i1=1,8
 Do i2=1,8
coup1 = cplHpmHpmcHpmcHpm(i1,i2,i1,i2)
dcoup1 = dcplHpmHpmcHpmcHpm(i1,i2,i1,i2,iv1)
Call FirstDerivativeVeff_balls(MHpm2(i1),MHpm2(i2),dMHpm2(i1,iv1),dMHpm2(i2,iv1)      & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(10)=results2(10) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Hpm, Hpm, conj[Hpm], conj[Hpm]]' 
  End Do
End Do
! ---- Hpm,Sd ----
Do i1=1,8
 Do i2=1,6
coup1 = cplHpmSdcHpmcSd(i1,i2,i1,i2)
dcoup1 = dcplHpmSdcHpmcSd(i1,i2,i1,i2,iv1)
Call FirstDerivativeVeff_balls(MHpm2(i1),MSd2(i2),dMHpm2(i1,iv1),dMSd2(i2,iv1)        & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-1._dp)
results2(11)=results2(11) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Hpm, Sd, conj[Hpm], conj[Sd]]' 
  End Do
End Do
! ---- Hpm,Su ----
Do i1=1,8
 Do i2=1,6
coup1 = cplHpmSucHpmcSu(i1,i2,i1,i2)
dcoup1 = dcplHpmSucHpmcSu(i1,i2,i1,i2,iv1)
Call FirstDerivativeVeff_balls(MHpm2(i1),MSu2(i2),dMHpm2(i1,iv1),dMSu2(i2,iv1)        & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-1._dp)
results2(12)=results2(12) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Hpm, Su, conj[Hpm], conj[Su]]' 
  End Do
End Do
! ---- Sd,Sd ----
Do i1=1,6
 Do i2=1,6
coup1 = cplSdSdcSdcSd(i1,i2,i1,i2)
dcoup1 = dcplSdSdcSdcSd(i1,i2,i1,i2,iv1)
Call FirstDerivativeVeff_balls(MSd2(i1),MSd2(i2),dMSd2(i1,iv1),dMSd2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(13)=results2(13) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Sd, Sd, conj[Sd], conj[Sd]]' 
  End Do
End Do
! ---- Sd,Su ----
Do i1=1,6
 Do i2=1,6
coup1 = cplSdSucSdcSu(i1,i2,i1,i2)
dcoup1 = dcplSdSucSdcSu(i1,i2,i1,i2,iv1)
Call FirstDerivativeVeff_balls(MSd2(i1),MSu2(i2),dMSd2(i1,iv1),dMSu2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-1._dp)
results2(14)=results2(14) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Sd, Su, conj[Sd], conj[Su]]' 
  End Do
End Do
! ---- Su,Su ----
Do i1=1,6
 Do i2=1,6
coup1 = cplSuSucSucSu(i1,i2,i1,i2)
dcoup1 = dcplSuSucSucSu(i1,i2,i1,i2,iv1)
Call FirstDerivativeVeff_balls(MSu2(i1),MSu2(i2),dMSu2(i1,iv1),dMSu2(i2,iv1)          & 
& ,coup1,dcoup1,'SS',Q2,temp)
coeff = (-0.5_dp)
results2(15)=results2(15) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at SS C[Su, Su, conj[Su], conj[Su]]' 
  End Do
End Do
! ----- diagrams of type VS, 2 ------ 
! ---- Sd,VG ----
Do i1=1,6
coup1 = cplSdcSdVGVG(i1,i1)
dcoup1 = dcplSdcSdVGVG(i1,i1,iv1)
Call FirstDerivativeVeff_balls(0._dp,MSd2(i1),ZeroC,dMSd2(i1,iv1),coup1,dcoup1,'VS',Q2,temp)
coeff = 0._dp
results2(16)=results2(16) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at VS C[Sd, VG, VG, conj[Sd]]' 
End Do
! ---- Su,VG ----
Do i1=1,6
coup1 = cplSucSuVGVG(i1,i1)
dcoup1 = dcplSucSuVGVG(i1,i1,iv1)
Call FirstDerivativeVeff_balls(0._dp,MSu2(i1),ZeroC,dMSu2(i1,iv1),coup1,dcoup1,'VS',Q2,temp)
coeff = 0._dp
results2(17)=results2(17) + coeff*temp
if (.not.(temp.eq.temp))  write(*,*) 'NaN at VS C[Su, VG, VG, conj[Su]]' 
End Do
result = sum(results1)+sum(results2) 


 
ti2L(iv1) = oo16pi2**2 * Real(result,dp) 
End Do 
End Subroutine FirstDerivativeEffPot2Loop 


Subroutine SecondDerivativeMassesCoups(i1,i2,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,              & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,gout)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Real(dp), Intent(out) :: gout(:) 
Real(dp) :: err, vevs(8) 
Integer :: iout 
Integer, Intent(in) :: i1,i2 
vevs = (/vd,vu,vR,vL/) 
gout = partialDiffXY_RiddersMulDim(AllMassesCouplings,112198,vevs,i1,i2,8,err) 
If (err.gt.err2L) err2L = err 
End Subroutine SecondDerivativeMassesCoups

Subroutine FirstDerivativeMassesCoups(i1,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,              & 
& Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,              & 
& M2,M3,gout)

Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Real(dp), Intent(out) :: gout(:) 
Real(dp) :: err, vevs(8) 
Integer :: iout 
Integer, Intent(in) :: i1 
vevs = (/vd,vu,vR,vL/) 
gout = partialDiff_RiddersMulDim(AllMassesCouplings,112198,vevs,i1,8,err) 
If (err.gt.err2L) err2L = err 
End Subroutine FirstDerivativeMassesCoups

Subroutine AllMassesCouplings(vevs,gout) 
Implicit None 
Real(dp), Intent(out) :: gout(112198) 
Real(dp), Intent(in) :: vevs(8) 
Integer :: kont 
Real(dp) :: vd,vu,vR(3),vL(3)

Complex(dp) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

vd=vevs(1) 
vu=vevs(2) 
vR=vevs(3:5)
vL=vevs(6:8)
Call WrapperForDerivatives(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,               & 
& Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,gout)

End Subroutine 

Subroutine WrapperForDerivatives(vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,            & 
& Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,gout)

Implicit None 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Real(dp), Intent(out) :: gout(:) 
Integer :: kont 
Real(dp) :: MSd(6),MSd2(6),MSu(6),MSu2(6),Mhh(8),Mhh2(8),MAh(8),MAh2(8),MHpm(8),MHpm2(8),         & 
& MChi(10),MChi2(10),MCha(5),MCha2(5),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,          & 
& MVZ,MVZ2,MVWm,MVWm2

Complex(dp) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

Call TreeMassesEffPot(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,               & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,.True.,kont)

Call CouplingsForEffPot3(lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZA,ZH,Ye,Te,               & 
& ZP,Yd,Td,ZD,Yu,Tu,ZU,g3,ZER,ZEL,UV,ZDL,ZDR,ZUL,ZUR,pG,cplAhAhAh,cplAhAhhh,             & 
& cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,cplhhSdcSd,        & 
& cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplcChaChaAhL,     & 
& cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,cplcFuFuAhL,           & 
& cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,cplcChaChahhL,     & 
& cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,cplChiFdcSdL,        & 
& cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,cplcFdChiSdL,     & 
& cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,cplcFdFdhhL,          & 
& cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplGluFucSuL,      & 
& cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdGluSdL,           & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,cplcFdFdVGL,        & 
& cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR)

Call CouplingsForEffPot4(lam,Yv,kap,ZA,ZH,Ye,ZP,Yd,ZD,Yu,ZU,g3,cplAhAhAhAh,           & 
& cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,cplhhhhHpmcHpm,       & 
& cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,cplHpmSucHpmcSu,           & 
& cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)

Call MassesCoupsToG(MSd,MSd2,MSu,MSu2,Mhh,Mhh2,MAh,MAh2,MHpm,MHpm2,MChi,              & 
& MChi2,MCha,MCha2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MVZ,MVZ2,MVWm,MVWm2,cplAhAhAh,           & 
& cplAhAhhh,cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,         & 
& cplhhSdcSd,cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,        & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,       & 
& cplcChaChahhL,cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,       & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,     & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplcFdFdhhL,cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplGluFucSuL,cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,       & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR,             & 
& cplAhAhAhAh,cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,          & 
& cplhhhhHpmcHpm,cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,            & 
& cplHpmSucHpmcSu,cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG,gout)

End Subroutine WrapperForDerivatives

Subroutine MassesCoupsToG(MSd,MSd2,MSu,MSu2,Mhh,Mhh2,MAh,MAh2,MHpm,MHpm2,             & 
& MChi,MChi2,MCha,MCha2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MVZ,MVZ2,MVWm,MVWm2,cplAhAhAh,      & 
& cplAhAhhh,cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,         & 
& cplhhSdcSd,cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,        & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,       & 
& cplcChaChahhL,cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,       & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,     & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplcFdFdhhL,cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplGluFucSuL,cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,       & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR,             & 
& cplAhAhAhAh,cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,          & 
& cplhhhhHpmcHpm,cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,            & 
& cplHpmSucHpmcSu,cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG,g)

Implicit None 
Real(dp), Intent(out) :: g(:) 
Integer :: i1,i2,i3,i4, sumI
Real(dp),Intent(in) :: MSd(6),MSd2(6),MSu(6),MSu2(6),Mhh(8),Mhh2(8),MAh(8),MAh2(8),MHpm(8),MHpm2(8),         & 
& MChi(10),MChi2(10),MCha(5),MCha2(5),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,          & 
& MVZ,MVZ2,MVWm,MVWm2

Complex(dp),Intent(in) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp),Intent(in) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

g(1:6) = MSd
g(7:12) = MSd2
g(13:18) = MSu
g(19:24) = MSu2
g(25:32) = Mhh
g(33:40) = Mhh2
g(41:48) = MAh
g(49:56) = MAh2
g(57:64) = MHpm
g(65:72) = MHpm2
g(73:82) = MChi
g(83:92) = MChi2
g(93:97) = MCha
g(98:102) = MCha2
g(103:105) = MFd
g(106:108) = MFd2
g(109:111) = MFu
g(112:114) = MFu2
g(115) = MGlu
g(116) = MGlu2
g(117) = MVZ
g(118) = MVZ2
g(119) = MVWm
g(120) = MVWm2
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
g(SumI+121) = Real(cplAhAhAh(i1,i2,i3), dp) 
g(SumI+122) = Aimag(cplAhAhAh(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
g(SumI+1145) = Real(cplAhAhhh(i1,i2,i3), dp) 
g(SumI+1146) = Aimag(cplAhAhhh(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
g(SumI+2169) = Real(cplAhhhhh(i1,i2,i3), dp) 
g(SumI+2170) = Aimag(cplAhhhhh(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
g(SumI+3193) = Real(cplAhHpmcHpm(i1,i2,i3), dp) 
g(SumI+3194) = Aimag(cplAhHpmcHpm(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
g(SumI+4217) = Real(cplAhSdcSd(i1,i2,i3), dp) 
g(SumI+4218) = Aimag(cplAhSdcSd(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
g(SumI+4793) = Real(cplAhSucSu(i1,i2,i3), dp) 
g(SumI+4794) = Aimag(cplAhSucSu(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
g(SumI+5369) = Real(cplhhhhhh(i1,i2,i3), dp) 
g(SumI+5370) = Aimag(cplhhhhhh(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
g(SumI+6393) = Real(cplhhHpmcHpm(i1,i2,i3), dp) 
g(SumI+6394) = Aimag(cplhhHpmcHpm(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
g(SumI+7417) = Real(cplhhSdcSd(i1,i2,i3), dp) 
g(SumI+7418) = Aimag(cplhhSdcSd(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
g(SumI+7993) = Real(cplhhSucSu(i1,i2,i3), dp) 
g(SumI+7994) = Aimag(cplhhSucSu(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
g(SumI+8569) = Real(cplHpmSucSd(i1,i2,i3), dp) 
g(SumI+8570) = Aimag(cplHpmSucSd(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,8
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*48
SumI = SumI*2 
g(SumI+9145) = Real(cplSdcHpmcSu(i1,i2,i3), dp) 
g(SumI+9146) = Aimag(cplSdcHpmcSu(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+9721) = Real(cplSdcSdVG(i1,i2), dp) 
g(SumI+9722) = Aimag(cplSdcSdVG(i1,i2)) 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+9793) = Real(cplSucSuVG(i1,i2), dp) 
g(SumI+9794) = Aimag(cplSucSuVG(i1,i2)) 
End Do 
End Do 

g(9865) = Real(cplVGVGVG,dp)  
g(9866) = Aimag(cplVGVGVG)  
Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
g(SumI+9867) = Real(cplcChaChaAhL(i1,i2,i3), dp) 
g(SumI+9868) = Aimag(cplcChaChaAhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
g(SumI+10267) = Real(cplcChaChaAhR(i1,i2,i3), dp) 
g(SumI+10268) = Aimag(cplcChaChaAhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
g(SumI+10667) = Real(cplChiChiAhL(i1,i2,i3), dp) 
g(SumI+10668) = Aimag(cplChiChiAhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
g(SumI+12267) = Real(cplChiChiAhR(i1,i2,i3), dp) 
g(SumI+12268) = Aimag(cplChiChiAhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+13867) = Real(cplcFdFdAhL(i1,i2,i3), dp) 
g(SumI+13868) = Aimag(cplcFdFdAhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+14011) = Real(cplcFdFdAhR(i1,i2,i3), dp) 
g(SumI+14012) = Aimag(cplcFdFdAhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+14155) = Real(cplcFuFuAhL(i1,i2,i3), dp) 
g(SumI+14156) = Aimag(cplcFuFuAhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+14299) = Real(cplcFuFuAhR(i1,i2,i3), dp) 
g(SumI+14300) = Aimag(cplcFuFuAhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
g(SumI+14443) = Real(cplChiChacHpmL(i1,i2,i3), dp) 
g(SumI+14444) = Aimag(cplChiChacHpmL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
g(SumI+15243) = Real(cplChiChacHpmR(i1,i2,i3), dp) 
g(SumI+15244) = Aimag(cplChiChacHpmR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+16043) = Real(cplChaFucSdL(i1,i2,i3), dp) 
g(SumI+16044) = Aimag(cplChaFucSdL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+16223) = Real(cplChaFucSdR(i1,i2,i3), dp) 
g(SumI+16224) = Aimag(cplChaFucSdR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
g(SumI+16403) = Real(cplcChaChahhL(i1,i2,i3), dp) 
g(SumI+16404) = Aimag(cplcChaChahhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
g(SumI+16803) = Real(cplcChaChahhR(i1,i2,i3), dp) 
g(SumI+16804) = Aimag(cplcChaChahhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,5
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*30
SumI = SumI*2 
g(SumI+17203) = Real(cplcFdChaSuL(i1,i2,i3), dp) 
g(SumI+17204) = Aimag(cplcFdChaSuL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,5
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*30
SumI = SumI*2 
g(SumI+17383) = Real(cplcFdChaSuR(i1,i2,i3), dp) 
g(SumI+17384) = Aimag(cplcFdChaSuR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
g(SumI+17563) = Real(cplChiChihhL(i1,i2,i3), dp) 
g(SumI+17564) = Aimag(cplChiChihhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
g(SumI+19163) = Real(cplChiChihhR(i1,i2,i3), dp) 
g(SumI+19164) = Aimag(cplChiChihhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+20763) = Real(cplChiFdcSdL(i1,i2,i3), dp) 
g(SumI+20764) = Aimag(cplChiFdcSdL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+21123) = Real(cplChiFdcSdR(i1,i2,i3), dp) 
g(SumI+21124) = Aimag(cplChiFdcSdR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+21483) = Real(cplChiFucSuL(i1,i2,i3), dp) 
g(SumI+21484) = Aimag(cplChiFucSuL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+21843) = Real(cplChiFucSuR(i1,i2,i3), dp) 
g(SumI+21844) = Aimag(cplChiFucSuR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
g(SumI+22203) = Real(cplcChaChiHpmL(i1,i2,i3), dp) 
g(SumI+22204) = Aimag(cplcChaChiHpmL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
g(SumI+23003) = Real(cplcChaChiHpmR(i1,i2,i3), dp) 
g(SumI+23004) = Aimag(cplcChaChiHpmR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
g(SumI+23803) = Real(cplcFdChiSdL(i1,i2,i3), dp) 
g(SumI+23804) = Aimag(cplcFdChiSdL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
g(SumI+24163) = Real(cplcFdChiSdR(i1,i2,i3), dp) 
g(SumI+24164) = Aimag(cplcFdChiSdR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
g(SumI+24523) = Real(cplcFuChiSuL(i1,i2,i3), dp) 
g(SumI+24524) = Aimag(cplcFuChiSuL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
g(SumI+24883) = Real(cplcFuChiSuR(i1,i2,i3), dp) 
g(SumI+24884) = Aimag(cplcFuChiSuR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+25243) = Real(cplGluFdcSdL(i1,i2), dp) 
g(SumI+25244) = Aimag(cplGluFdcSdL(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+25279) = Real(cplGluFdcSdR(i1,i2), dp) 
g(SumI+25280) = Aimag(cplGluFdcSdR(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+25315) = Real(cplcFdFdhhL(i1,i2,i3), dp) 
g(SumI+25316) = Aimag(cplcFdFdhhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+25459) = Real(cplcFdFdhhR(i1,i2,i3), dp) 
g(SumI+25460) = Aimag(cplcFdFdhhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+25603) = Real(cplcChaFdcSuL(i1,i2,i3), dp) 
g(SumI+25604) = Aimag(cplcChaFdcSuL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+25783) = Real(cplcChaFdcSuR(i1,i2,i3), dp) 
g(SumI+25784) = Aimag(cplcChaFdcSuR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+25963) = Real(cplcFuFdcHpmL(i1,i2,i3), dp) 
g(SumI+25964) = Aimag(cplcFuFdcHpmL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+26107) = Real(cplcFuFdcHpmR(i1,i2,i3), dp) 
g(SumI+26108) = Aimag(cplcFuFdcHpmR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+26251) = Real(cplGluFucSuL(i1,i2), dp) 
g(SumI+26252) = Aimag(cplGluFucSuL(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+26287) = Real(cplGluFucSuR(i1,i2), dp) 
g(SumI+26288) = Aimag(cplGluFucSuR(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+26323) = Real(cplcFuFuhhL(i1,i2,i3), dp) 
g(SumI+26324) = Aimag(cplcFuFuhhL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+26467) = Real(cplcFuFuhhR(i1,i2,i3), dp) 
g(SumI+26468) = Aimag(cplcFuFuhhR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+26611) = Real(cplcFdFuHpmL(i1,i2,i3), dp) 
g(SumI+26612) = Aimag(cplcFdFuHpmL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
g(SumI+26755) = Real(cplcFdFuHpmR(i1,i2,i3), dp) 
g(SumI+26756) = Aimag(cplcFdFuHpmR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+26899) = Real(cplcFdGluSdL(i1,i2), dp) 
g(SumI+26900) = Aimag(cplcFdGluSdL(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+26935) = Real(cplcFdGluSdR(i1,i2), dp) 
g(SumI+26936) = Aimag(cplcFdGluSdR(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+26971) = Real(cplcFuGluSuL(i1,i2), dp) 
g(SumI+26972) = Aimag(cplcFuGluSuL(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+27007) = Real(cplcFuGluSuR(i1,i2), dp) 
g(SumI+27008) = Aimag(cplcFuGluSuR(i1,i2)) 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+27043) = Real(cplcChacFuSdL(i1,i2,i3), dp) 
g(SumI+27044) = Aimag(cplcChacFuSdL(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
g(SumI+27223) = Real(cplcChacFuSdR(i1,i2,i3), dp) 
g(SumI+27224) = Aimag(cplcChacFuSdR(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+27403) = Real(cplcFdFdVGL(i1,i2), dp) 
g(SumI+27404) = Aimag(cplcFdFdVGL(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+27421) = Real(cplcFdFdVGR(i1,i2), dp) 
g(SumI+27422) = Aimag(cplcFdFdVGR(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+27439) = Real(cplcFuFuVGL(i1,i2), dp) 
g(SumI+27440) = Aimag(cplcFuFuVGL(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+27457) = Real(cplcFuFuVGR(i1,i2), dp) 
g(SumI+27458) = Aimag(cplcFuFuVGR(i1,i2)) 
End Do 
End Do 

g(27475) = Real(cplGluGluVGL,dp)  
g(27476) = Aimag(cplGluGluVGL)  
g(27477) = Real(cplGluGluVGR,dp)  
g(27478) = Aimag(cplGluGluVGR)  
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
g(SumI+27479) = Real(cplAhAhAhAh(i1,i2,i3,i4), dp) 
g(SumI+27480) = Aimag(cplAhAhAhAh(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
g(SumI+35671) = Real(cplAhAhhhhh(i1,i2,i3,i4), dp) 
g(SumI+35672) = Aimag(cplAhAhhhhh(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
g(SumI+43863) = Real(cplAhAhHpmcHpm(i1,i2,i3,i4), dp) 
g(SumI+43864) = Aimag(cplAhAhHpmcHpm(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
g(SumI+52055) = Real(cplAhAhSdcSd(i1,i2,i3,i4), dp) 
g(SumI+52056) = Aimag(cplAhAhSdcSd(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
g(SumI+56663) = Real(cplAhAhSucSu(i1,i2,i3,i4), dp) 
g(SumI+56664) = Aimag(cplAhAhSucSu(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
g(SumI+61271) = Real(cplhhhhhhhh(i1,i2,i3,i4), dp) 
g(SumI+61272) = Aimag(cplhhhhhhhh(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
g(SumI+69463) = Real(cplhhhhHpmcHpm(i1,i2,i3,i4), dp) 
g(SumI+69464) = Aimag(cplhhhhHpmcHpm(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
g(SumI+77655) = Real(cplhhhhSdcSd(i1,i2,i3,i4), dp) 
g(SumI+77656) = Aimag(cplhhhhSdcSd(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
g(SumI+82263) = Real(cplhhhhSucSu(i1,i2,i3,i4), dp) 
g(SumI+82264) = Aimag(cplhhhhSucSu(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
g(SumI+86871) = Real(cplHpmHpmcHpmcHpm(i1,i2,i3,i4), dp) 
g(SumI+86872) = Aimag(cplHpmHpmcHpmcHpm(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,8
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*48 + (i1-1)*288
SumI = SumI*2 
g(SumI+95063) = Real(cplHpmSdcHpmcSd(i1,i2,i3,i4), dp) 
g(SumI+95064) = Aimag(cplHpmSdcHpmcSd(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,8
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*48 + (i1-1)*288
SumI = SumI*2 
g(SumI+99671) = Real(cplHpmSucHpmcSu(i1,i2,i3,i4), dp) 
g(SumI+99672) = Aimag(cplHpmSucHpmcSu(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*216
SumI = SumI*2 
g(SumI+104279) = Real(cplSdSdcSdcSd(i1,i2,i3,i4), dp) 
g(SumI+104280) = Aimag(cplSdSdcSdcSd(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*216
SumI = SumI*2 
g(SumI+106871) = Real(cplSdSucSdcSu(i1,i2,i3,i4), dp) 
g(SumI+106872) = Aimag(cplSdSucSdcSu(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*216
SumI = SumI*2 
g(SumI+109463) = Real(cplSuSucSucSu(i1,i2,i3,i4), dp) 
g(SumI+109464) = Aimag(cplSuSucSucSu(i1,i2,i3,i4)) 
End Do 
End Do 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+112055) = Real(cplSdcSdVGVG(i1,i2), dp) 
g(SumI+112056) = Aimag(cplSdcSdVGVG(i1,i2)) 
End Do 
End Do 

Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
g(SumI+112127) = Real(cplSucSuVGVG(i1,i2), dp) 
g(SumI+112128) = Aimag(cplSucSuVGVG(i1,i2)) 
End Do 
End Do 

End Subroutine MassesCoupsToG

Subroutine GToMassesCoups(g,MSd,MSd2,MSu,MSu2,Mhh,Mhh2,MAh,MAh2,MHpm,MHpm2,           & 
& MChi,MChi2,MCha,MCha2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MVZ,MVZ2,MVWm,MVWm2,cplAhAhAh,      & 
& cplAhAhhh,cplAhhhhh,cplAhHpmcHpm,cplAhSdcSd,cplAhSucSu,cplhhhhhh,cplhhHpmcHpm,         & 
& cplhhSdcSd,cplhhSucSu,cplHpmSucSd,cplSdcHpmcSu,cplSdcSdVG,cplSucSuVG,cplVGVGVG,        & 
& cplcChaChaAhL,cplcChaChaAhR,cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,         & 
& cplcFuFuAhL,cplcFuFuAhR,cplChiChacHpmL,cplChiChacHpmR,cplChaFucSdL,cplChaFucSdR,       & 
& cplcChaChahhL,cplcChaChahhR,cplcFdChaSuL,cplcFdChaSuR,cplChiChihhL,cplChiChihhR,       & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplcChaChiHpmL,cplcChaChiHpmR,     & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFuChiSuL,cplcFuChiSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplcFdFdhhL,cplcFdFdhhR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplGluFucSuL,cplGluFucSuR,cplcFuFuhhL,cplcFuFuhhR,cplcFdFuHpmL,cplcFdFuHpmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,       & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplGluGluVGL,cplGluGluVGR,             & 
& cplAhAhAhAh,cplAhAhhhhh,cplAhAhHpmcHpm,cplAhAhSdcSd,cplAhAhSucSu,cplhhhhhhhh,          & 
& cplhhhhHpmcHpm,cplhhhhSdcSd,cplhhhhSucSu,cplHpmHpmcHpmcHpm,cplHpmSdcHpmcSd,            & 
& cplHpmSucHpmcSu,cplSdSdcSdcSd,cplSdSucSdcSu,cplSuSucSucSu,cplSdcSdVGVG,cplSucSuVGVG)

Implicit None 
Real(dp), Intent(in) :: g(:) 
Integer :: i1,i2,i3,i4, sumI
Real(dp),Intent(inout) :: MSd(6),MSd2(6),MSu(6),MSu2(6),Mhh(8),Mhh2(8),MAh(8),MAh2(8),MHpm(8),MHpm2(8),         & 
& MChi(10),MChi2(10),MCha(5),MCha2(5),MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,          & 
& MVZ,MVZ2,MVWm,MVWm2

Complex(dp),Intent(inout) :: cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),cplAhhhhh(8,8,8),cplAhHpmcHpm(8,8,8),               & 
& cplAhSdcSd(8,6,6),cplAhSucSu(8,6,6),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),              & 
& cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplHpmSucSd(8,6,6),cplSdcHpmcSu(6,8,6),            & 
& cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),   & 
& cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),     & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),   & 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),     & 
& cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),& 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),         & 
& cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),       & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplGluFucSuL(3,6),cplGluFucSuR(3,6),         & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),         & 
& cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),               & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),           & 
& cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplGluGluVGL,cplGluGluVGR

Complex(dp),Intent(inout) :: cplAhAhAhAh(8,8,8,8),cplAhAhhhhh(8,8,8,8),cplAhAhHpmcHpm(8,8,8,8),cplAhAhSdcSd(8,8,6,6),& 
& cplAhAhSucSu(8,8,6,6),cplhhhhhhhh(8,8,8,8),cplhhhhHpmcHpm(8,8,8,8),cplhhhhSdcSd(8,8,6,6),& 
& cplhhhhSucSu(8,8,6,6),cplHpmHpmcHpmcHpm(8,8,8,8),cplHpmSdcHpmcSd(8,6,8,6),             & 
& cplHpmSucHpmcSu(8,6,8,6),cplSdSdcSdcSd(6,6,6,6),cplSdSucSdcSu(6,6,6,6),cplSuSucSucSu(6,6,6,6),& 
& cplSdcSdVGVG(6,6),cplSucSuVGVG(6,6)

MSd=g(1:6) 
MSd2=g(7:12) 
MSu=g(13:18) 
MSu2=g(19:24) 
Mhh=g(25:32) 
Mhh2=g(33:40) 
MAh=g(41:48) 
MAh2=g(49:56) 
MHpm=g(57:64) 
MHpm2=g(65:72) 
MChi=g(73:82) 
MChi2=g(83:92) 
MCha=g(93:97) 
MCha2=g(98:102) 
MFd=g(103:105) 
MFd2=g(106:108) 
MFu=g(109:111) 
MFu2=g(112:114) 
MGlu=g(115) 
MGlu2=g(116) 
MVZ=g(117) 
MVZ2=g(118) 
MVWm=g(119) 
MVWm2=g(120) 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
cplAhAhAh(i1,i2,i3) = Cmplx( g(SumI+121), g(SumI+122), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
cplAhAhhh(i1,i2,i3) = Cmplx( g(SumI+1145), g(SumI+1146), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
cplAhhhhh(i1,i2,i3) = Cmplx( g(SumI+2169), g(SumI+2170), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
cplAhHpmcHpm(i1,i2,i3) = Cmplx( g(SumI+3193), g(SumI+3194), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
cplAhSdcSd(i1,i2,i3) = Cmplx( g(SumI+4217), g(SumI+4218), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
cplAhSucSu(i1,i2,i3) = Cmplx( g(SumI+4793), g(SumI+4794), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
cplhhhhhh(i1,i2,i3) = Cmplx( g(SumI+5369), g(SumI+5370), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*64
SumI = SumI*2 
cplhhHpmcHpm(i1,i2,i3) = Cmplx( g(SumI+6393), g(SumI+6394), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
cplhhSdcSd(i1,i2,i3) = Cmplx( g(SumI+7417), g(SumI+7418), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
cplhhSucSu(i1,i2,i3) = Cmplx( g(SumI+7993), g(SumI+7994), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*36
SumI = SumI*2 
cplHpmSucSd(i1,i2,i3) = Cmplx( g(SumI+8569), g(SumI+8570), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,8
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*48
SumI = SumI*2 
cplSdcHpmcSu(i1,i2,i3) = Cmplx( g(SumI+9145), g(SumI+9146), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplSdcSdVG(i1,i2) = Cmplx( g(SumI+9721), g(SumI+9722), dp) 
End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplSucSuVG(i1,i2) = Cmplx( g(SumI+9793), g(SumI+9794), dp) 
End Do 
 End Do 
 
cplVGVGVG= Cmplx(g(9865),g(9866),dp) 
Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
cplcChaChaAhL(i1,i2,i3) = Cmplx( g(SumI+9867), g(SumI+9868), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
cplcChaChaAhR(i1,i2,i3) = Cmplx( g(SumI+10267), g(SumI+10268), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
cplChiChiAhL(i1,i2,i3) = Cmplx( g(SumI+10667), g(SumI+10668), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
cplChiChiAhR(i1,i2,i3) = Cmplx( g(SumI+12267), g(SumI+12268), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFdFdAhL(i1,i2,i3) = Cmplx( g(SumI+13867), g(SumI+13868), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFdFdAhR(i1,i2,i3) = Cmplx( g(SumI+14011), g(SumI+14012), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFuFuAhL(i1,i2,i3) = Cmplx( g(SumI+14155), g(SumI+14156), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFuFuAhR(i1,i2,i3) = Cmplx( g(SumI+14299), g(SumI+14300), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
cplChiChacHpmL(i1,i2,i3) = Cmplx( g(SumI+14443), g(SumI+14444), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
cplChiChacHpmR(i1,i2,i3) = Cmplx( g(SumI+15243), g(SumI+15244), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplChaFucSdL(i1,i2,i3) = Cmplx( g(SumI+16043), g(SumI+16044), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplChaFucSdR(i1,i2,i3) = Cmplx( g(SumI+16223), g(SumI+16224), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
cplcChaChahhL(i1,i2,i3) = Cmplx( g(SumI+16403), g(SumI+16404), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,5
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*40
SumI = SumI*2 
cplcChaChahhR(i1,i2,i3) = Cmplx( g(SumI+16803), g(SumI+16804), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,5
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*30
SumI = SumI*2 
cplcFdChaSuL(i1,i2,i3) = Cmplx( g(SumI+17203), g(SumI+17204), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,5
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*30
SumI = SumI*2 
cplcFdChaSuR(i1,i2,i3) = Cmplx( g(SumI+17383), g(SumI+17384), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
cplChiChihhL(i1,i2,i3) = Cmplx( g(SumI+17563), g(SumI+17564), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
cplChiChihhR(i1,i2,i3) = Cmplx( g(SumI+19163), g(SumI+19164), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplChiFdcSdL(i1,i2,i3) = Cmplx( g(SumI+20763), g(SumI+20764), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplChiFdcSdR(i1,i2,i3) = Cmplx( g(SumI+21123), g(SumI+21124), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplChiFucSuL(i1,i2,i3) = Cmplx( g(SumI+21483), g(SumI+21484), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,10
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplChiFucSuR(i1,i2,i3) = Cmplx( g(SumI+21843), g(SumI+21844), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
cplcChaChiHpmL(i1,i2,i3) = Cmplx( g(SumI+22203), g(SumI+22204), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,10
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*80
SumI = SumI*2 
cplcChaChiHpmR(i1,i2,i3) = Cmplx( g(SumI+23003), g(SumI+23004), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
cplcFdChiSdL(i1,i2,i3) = Cmplx( g(SumI+23803), g(SumI+23804), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
cplcFdChiSdR(i1,i2,i3) = Cmplx( g(SumI+24163), g(SumI+24164), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
cplcFuChiSuL(i1,i2,i3) = Cmplx( g(SumI+24523), g(SumI+24524), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,10
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*60
SumI = SumI*2 
cplcFuChiSuR(i1,i2,i3) = Cmplx( g(SumI+24883), g(SumI+24884), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplGluFdcSdL(i1,i2) = Cmplx( g(SumI+25243), g(SumI+25244), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplGluFdcSdR(i1,i2) = Cmplx( g(SumI+25279), g(SumI+25280), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFdFdhhL(i1,i2,i3) = Cmplx( g(SumI+25315), g(SumI+25316), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFdFdhhR(i1,i2,i3) = Cmplx( g(SumI+25459), g(SumI+25460), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplcChaFdcSuL(i1,i2,i3) = Cmplx( g(SumI+25603), g(SumI+25604), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplcChaFdcSuR(i1,i2,i3) = Cmplx( g(SumI+25783), g(SumI+25784), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFuFdcHpmL(i1,i2,i3) = Cmplx( g(SumI+25963), g(SumI+25964), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFuFdcHpmR(i1,i2,i3) = Cmplx( g(SumI+26107), g(SumI+26108), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplGluFucSuL(i1,i2) = Cmplx( g(SumI+26251), g(SumI+26252), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplGluFucSuR(i1,i2) = Cmplx( g(SumI+26287), g(SumI+26288), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFuFuhhL(i1,i2,i3) = Cmplx( g(SumI+26323), g(SumI+26324), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFuFuhhR(i1,i2,i3) = Cmplx( g(SumI+26467), g(SumI+26468), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFdFuHpmL(i1,i2,i3) = Cmplx( g(SumI+26611), g(SumI+26612), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,8
SumI = (i3-1) + (i2-1)*8 + (i1-1)*24
SumI = SumI*2 
cplcFdFuHpmR(i1,i2,i3) = Cmplx( g(SumI+26755), g(SumI+26756), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplcFdGluSdL(i1,i2) = Cmplx( g(SumI+26899), g(SumI+26900), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplcFdGluSdR(i1,i2) = Cmplx( g(SumI+26935), g(SumI+26936), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplcFuGluSuL(i1,i2) = Cmplx( g(SumI+26971), g(SumI+26972), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplcFuGluSuR(i1,i2) = Cmplx( g(SumI+27007), g(SumI+27008), dp) 
End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplcChacFuSdL(i1,i2,i3) = Cmplx( g(SumI+27043), g(SumI+27044), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,5
Do i2 = 1,3
Do i3 = 1,6
SumI = (i3-1) + (i2-1)*6 + (i1-1)*18
SumI = SumI*2 
cplcChacFuSdR(i1,i2,i3) = Cmplx( g(SumI+27223), g(SumI+27224), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
cplcFdFdVGL(i1,i2) = Cmplx( g(SumI+27403), g(SumI+27404), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
cplcFdFdVGR(i1,i2) = Cmplx( g(SumI+27421), g(SumI+27422), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
cplcFuFuVGL(i1,i2) = Cmplx( g(SumI+27439), g(SumI+27440), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
cplcFuFuVGR(i1,i2) = Cmplx( g(SumI+27457), g(SumI+27458), dp) 
End Do 
 End Do 
 
cplGluGluVGL= Cmplx(g(27475),g(27476),dp) 
cplGluGluVGR= Cmplx(g(27477),g(27478),dp) 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
cplAhAhAhAh(i1,i2,i3,i4) = Cmplx( g(SumI+27479), g(SumI+27480), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
cplAhAhhhhh(i1,i2,i3,i4) = Cmplx( g(SumI+35671), g(SumI+35672), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
cplAhAhHpmcHpm(i1,i2,i3,i4) = Cmplx( g(SumI+43863), g(SumI+43864), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
cplAhAhSdcSd(i1,i2,i3,i4) = Cmplx( g(SumI+52055), g(SumI+52056), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
cplAhAhSucSu(i1,i2,i3,i4) = Cmplx( g(SumI+56663), g(SumI+56664), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
cplhhhhhhhh(i1,i2,i3,i4) = Cmplx( g(SumI+61271), g(SumI+61272), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
cplhhhhHpmcHpm(i1,i2,i3,i4) = Cmplx( g(SumI+69463), g(SumI+69464), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
cplhhhhSdcSd(i1,i2,i3,i4) = Cmplx( g(SumI+77655), g(SumI+77656), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*288
SumI = SumI*2 
cplhhhhSucSu(i1,i2,i3,i4) = Cmplx( g(SumI+82263), g(SumI+82264), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,8
Do i3 = 1,8
Do i4 = 1,8
SumI = (i4-1) + (i3-1)*8 + (i2-1)*64 + (i1-1)*512
SumI = SumI*2 
cplHpmHpmcHpmcHpm(i1,i2,i3,i4) = Cmplx( g(SumI+86871), g(SumI+86872), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,8
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*48 + (i1-1)*288
SumI = SumI*2 
cplHpmSdcHpmcSd(i1,i2,i3,i4) = Cmplx( g(SumI+95063), g(SumI+95064), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,8
Do i2 = 1,6
Do i3 = 1,8
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*48 + (i1-1)*288
SumI = SumI*2 
cplHpmSucHpmcSu(i1,i2,i3,i4) = Cmplx( g(SumI+99671), g(SumI+99672), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*216
SumI = SumI*2 
cplSdSdcSdcSd(i1,i2,i3,i4) = Cmplx( g(SumI+104279), g(SumI+104280), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*216
SumI = SumI*2 
cplSdSucSdcSu(i1,i2,i3,i4) = Cmplx( g(SumI+106871), g(SumI+106872), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
Do i3 = 1,6
Do i4 = 1,6
SumI = (i4-1) + (i3-1)*6 + (i2-1)*36 + (i1-1)*216
SumI = SumI*2 
cplSuSucSucSu(i1,i2,i3,i4) = Cmplx( g(SumI+109463), g(SumI+109464), dp) 
End Do 
 End Do 
 End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplSdcSdVGVG(i1,i2) = Cmplx( g(SumI+112055), g(SumI+112056), dp) 
End Do 
 End Do 
 
Do i1 = 1,6
Do i2 = 1,6
SumI = (i2-1) + (i1-1)*6
SumI = SumI*2 
cplSucSuVGVG(i1,i2) = Cmplx( g(SumI+112127), g(SumI+112128), dp) 
End Do 
 End Do 
 
End Subroutine GToMassesCoups

End Module EffectivePotential_munuSSM3G 
 
