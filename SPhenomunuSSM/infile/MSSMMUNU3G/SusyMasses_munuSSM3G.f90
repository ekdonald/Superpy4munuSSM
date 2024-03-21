! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:20 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module SusyMasses_munuSSM3G 
 
Use Control 
Use Mathematics 
Use MathematicsQP 
Use Model_Data_munuSSM3G 
!Use StandardModel 
 
 
Logical :: SignOfMassChanged =.False.  
Logical :: SignOfMuChanged =.False.  
Contains 
 
Subroutine TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,               & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,               & 
& Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,             & 
& mlHd2,M1,M2,M3,GenerationMixing,kont)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(out) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(out) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Logical, Intent(in) :: GenerationMixing 
Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,kontSave

character(len=20) :: varcheck
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'TreeMassesmunuSSM'
 
kont = 0 
Call CalculateVPVZ(g1,g2,vd,vu,vL,ZZ,MVZ,MVZ2,TW,kont)

Call CalculateVWm(g2,vd,vu,vL,ZW,MVWm,MVWm2,kont)

MGlu = M3 
pG = Abs(MGlu)/MGlu
pG = Sqrt(pG)
MGlu = Abs(MGlu) 
MGlu2 = MGlu**2 
Call CalculateMSd(g1,g2,Yd,Td,lam,mq2,md2,vd,vu,vL,vR,ZD,MSd,MSd2,kont)

Call CalculateMSu(g1,g2,lam,Yv,Yu,Tu,mq2,mu2,vd,vu,vL,vR,ZU,MSu,MSu2,kont)

Call CalculateMhh(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,mHu2,mv2,vd,             & 
& vu,vL,vR,ZH,Mhh,Mhh2,kont)

kontSave = kont 
Call CalculateMAh(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,mHu2,mv2,vd,             & 
& vu,vL,vR,TW,ZA,MAh,MAh2,kont)

kont = kontSave 
kontSave = kont 
Call CalculateMHpm(g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,ml2,mlHd2,mHd2,mHu2,me2,            & 
& vd,vu,vL,vR,ZP,MHpm,MHpm2,kont)

kont = kontSave 
Call CalculateMChi(g1,g2,lam,Yv,kap,M1,M2,vd,vu,vL,vR,UV,MChi,kont)

MChi2 = MChi**2 
Call CalculateMCha(g2,Ye,lam,Yv,M2,vd,vu,vL,vR,ZER,ZEL,MCha,kont)

MCha2 = MCha**2 
Call CalculateMFd(Yd,vd,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
Call CalculateMFu(Yu,vu,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 
 
 Call SortGoldstones(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,            & 
& MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,             & 
& ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,kont)

!!-----------------------------------------------------------
!Write(*,*) "Tree Mhh", Mhh(1), Mhh(2), Mhh(3), Mhh(4), Mhh(5)
!Write(*,*) "Tree MAh", MAh(1), MAh(2), MAh(3), MAh(4), MAh(5)
!Write(*,*) "Tree MHpm", MHpm(1), MHpm(2), MHpm(3), MHpm(4), MHpm(5)
!Write(*,*) "Tree Mch1", Mchi(1), Mchi(3), Mchi(4), Mchi(5), Mchi(6)
!Write(*,*) "Tree Mcha", Mcha(1), Mcha(2), Mcha(3), Mcha(4), Mcha(5)
!Write(*,*) "  "
!!-----------------------------------------------------------

varcheck = "no"


If ((HighScaleModel.Eq."LOW").and.(.not.SUSYrunningFromMZ)) Then 
 If (SignOfMassChanged) Then  
 If (.Not.IgnoreNegativeMasses) Then 
  Write(*,*) " Stopping 1A calculation because of negative mass squared." 

   !Open(375, action='write',file= 'checktachyons.dat',status='unknown')
   !write(375,*) "yes"
   !close(375)
   varcheck = "yes"

  Call TerminateProgram 
  return
 Else 
  SignOfMassChanged= .False. 
  kont=0  
 End If 
End If 
If (SignOfMuChanged) Then 
 If (.Not.IgnoreMuSignFlip) Then 
  Write(*,*) " Stopping 2B calculation because of negative mass squared in tadpoles." 
  Call TerminateProgram 
  return
 Else 
  SignOfMuChanged= .False. 
  kont=0 
 End If 
End If 
End if 
 
!! -------------------------------- 
! Setting Goldstone masses 
! -------------------------------- 
 
TW = ACos(Abs(ZZ(1,1)))
Iname = Iname - 1 


   Open(375, action='write',file= 'checktachyons.dat',status='unknown')
   write(375,*) varcheck
   close(375)

 
End Subroutine  TreeMasses 
 
 
Subroutine TreeMassesEffPot(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,              & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,             & 
& pG,TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,               & 
& g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,              & 
& mv2,mlHd2,M1,M2,M3,GenerationMixing,kont)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(out) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(out) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Logical, Intent(in) :: GenerationMixing 
Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,kontSave 
Iname = Iname + 1 
NameOfUnit(Iname) = 'TreeMassesmunuSSM'
 
kont = 0 
Call CalculateVPVZEffPot(g1,g2,vd,vu,vL,ZZ,MVZ,MVZ2,TW,kont)

Call CalculateVWmEffPot(g2,vd,vu,vL,ZW,MVWm,MVWm2,kont)

MGlu = M3 
pG = Abs(MGlu)/MGlu
pG = Sqrt(pG)
MGlu = Abs(MGlu) 
MGlu2 = MGlu**2 
Call CalculateMSdEffPot(g1,g2,Yd,Td,lam,mq2,md2,vd,vu,vL,vR,ZD,MSd,MSd2,kont)

Call CalculateMSuEffPot(g1,g2,lam,Yv,Yu,Tu,mq2,mu2,vd,vu,vL,vR,ZU,MSu,MSu2,kont)

Call CalculateMhhEffPot(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,mHu2,              & 
& mv2,vd,vu,vL,vR,ZH,Mhh,Mhh2,kont)

kontSave = kont 
Call CalculateMAhEffPot(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,mHu2,              & 
& mv2,vd,vu,vL,vR,TW,ZA,MAh,MAh2,kont)

kont = kontSave 
kontSave = kont 
Call CalculateMHpmEffPot(g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,ml2,mlHd2,mHd2,               & 
& mHu2,me2,vd,vu,vL,vR,ZP,MHpm,MHpm2,kont)

kont = kontSave 
Call CalculateMChiEffPot(g1,g2,lam,Yv,kap,M1,M2,vd,vu,vL,vR,UV,MChi,kont)

MChi2 = MChi**2 
Call CalculateMChaEffPot(g2,Ye,lam,Yv,M2,vd,vu,vL,vR,ZER,ZEL,MCha,kont)

MCha2 = MCha**2 
Call CalculateMFdEffPot(Yd,vd,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
Call CalculateMFuEffPot(Yu,vu,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 

 
 If ((HighScaleModel.Eq."LOW").and.(.not.SUSYrunningFromMZ)) Then 
 If (SignOfMassChanged) Then  
 If (.Not.IgnoreNegativeMasses) Then 
  Write(*,*) " Stopping 3C calculation because of negative mass squared." 
  Call TerminateProgram 
  return
 Else 
  SignOfMassChanged= .False. 
  kont=0  
 End If 
End If 
If (SignOfMuChanged) Then 
 If (.Not.IgnoreMuSignFlip) Then 
  Write(*,*) " Stopping 4D calculation because of negative mass squared in tadpoles." 
  Call TerminateProgram 
  return
 Else 
  SignOfMuChanged= .False. 
  kont=0 
 End If 
End If 
End if 
 Iname = Iname - 1 
 
End Subroutine  TreeMassesEffPot 
 
 
Subroutine CalculateMSd(g1,g2,Yd,Td,lam,mq2,md2,vd,vu,vL,vR,ZD,MSd,MSd2,kont)

Real(dp), Intent(in) :: g1,g2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: Yd(3,3),Td(3,3),lam(3),mq2(3,3),md2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MSd(6), MSd2(6) 
Complex(dp), Intent(out) :: ZD(6,6) 
 
Complex(dp) :: mat(6,6)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMSd'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)-(g1**2*vd**2)/24._dp
mat(1,1) = mat(1,1)-(g2**2*vd**2)/8._dp
mat(1,1) = mat(1,1)+(g1**2*vu**2)/24._dp
mat(1,1) = mat(1,1)+(g2**2*vu**2)/8._dp
mat(1,1) = mat(1,1)+mq2(1,1)
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+mq2(1,2)
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,2))/2._dp
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+mq2(1,3)
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,3))/2._dp
End Do 
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+(vd*Conjg(Td(1,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(Yd(1,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+(vd*Conjg(Td(2,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(Yd(2,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+(vd*Conjg(Td(3,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu*Conjg(Yd(3,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)-(g1**2*vd**2)/24._dp
mat(2,2) = mat(2,2)-(g2**2*vd**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vu**2)/24._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+mq2(2,2)
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vd**2*Conjg(Yd(j1,2))*Yd(j1,2))/2._dp
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+mq2(2,3)
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd**2*Conjg(Yd(j1,2))*Yd(j1,3))/2._dp
End Do 
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(vd*Conjg(Td(1,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vu*Conjg(Yd(1,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(vd*Conjg(Td(2,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vu*Conjg(Yd(2,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(2,6) = 0._dp 
mat(2,6) = mat(2,6)+(vd*Conjg(Td(3,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vu*Conjg(Yd(3,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)-(g1**2*vd**2)/24._dp
mat(3,3) = mat(3,3)-(g2**2*vd**2)/8._dp
mat(3,3) = mat(3,3)+(g1**2*vu**2)/24._dp
mat(3,3) = mat(3,3)+(g2**2*vu**2)/8._dp
mat(3,3) = mat(3,3)+mq2(3,3)
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd**2*Conjg(Yd(j1,3))*Yd(j1,3))/2._dp
End Do 
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(vd*Conjg(Td(1,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yd(1,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(vd*Conjg(Td(2,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yd(2,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vd*Conjg(Td(3,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yd(3,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)-(g1**2*vd**2)/12._dp
mat(4,4) = mat(4,4)+(g1**2*vu**2)/12._dp
mat(4,4) = mat(4,4)+md2(1,1)
Do j1 = 1,3
mat(4,4) = mat(4,4)-(g1**2*vL(j1)**2)/12._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd**2*Conjg(Yd(1,j1))*Yd(1,j1))/2._dp
End Do 
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+md2(1,2)
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd**2*Conjg(Yd(2,j1))*Yd(1,j1))/2._dp
End Do 
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+md2(1,3)
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vd**2*Conjg(Yd(3,j1))*Yd(1,j1))/2._dp
End Do 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)-(g1**2*vd**2)/12._dp
mat(5,5) = mat(5,5)+(g1**2*vu**2)/12._dp
mat(5,5) = mat(5,5)+md2(2,2)
Do j1 = 1,3
mat(5,5) = mat(5,5)-(g1**2*vL(j1)**2)/12._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd**2*Conjg(Yd(2,j1))*Yd(2,j1))/2._dp
End Do 
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+md2(2,3)
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vd**2*Conjg(Yd(3,j1))*Yd(2,j1))/2._dp
End Do 
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)-(g1**2*vd**2)/12._dp
mat(6,6) = mat(6,6)+(g1**2*vu**2)/12._dp
mat(6,6) = mat(6,6)+md2(3,3)
Do j1 = 1,3
mat(6,6) = mat(6,6)-(g1**2*vL(j1)**2)/12._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vd**2*Conjg(Yd(3,j1))*Yd(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,MSd2,ZD,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,6
  If (Abs(MSd2(i1)).Le.MaxMassNumericalZero) MSd2(i1) = 1.E-10_dp 
  If (MSd2(i1).ne.MSd2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (MSd2(i1).Ge.0._dp) Then 
  MSd(i1)=Sqrt(MSd2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MSd2(i1) 
    End If 
  MSd = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MSd2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MSd2(i1) 
  MSd2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMSd 

Subroutine CalculateMSu(g1,g2,lam,Yv,Yu,Tu,mq2,mu2,vd,vu,vL,vR,ZU,MSu,MSu2,kont)

Real(dp), Intent(in) :: g1,g2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: lam(3),Yv(3,3),Yu(3,3),Tu(3,3),mq2(3,3),mu2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MSu(6), MSu2(6) 
Complex(dp), Intent(out) :: ZU(6,6) 
 
Complex(dp) :: mat(6,6)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMSu'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)-(g1**2*vd**2)/24._dp
mat(1,1) = mat(1,1)+(g2**2*vd**2)/8._dp
mat(1,1) = mat(1,1)+(g1**2*vu**2)/24._dp
mat(1,1) = mat(1,1)-(g2**2*vu**2)/8._dp
mat(1,1) = mat(1,1)+mq2(1,1)
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+mq2(1,2)
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,2))/2._dp
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+mq2(1,3)
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,3))/2._dp
End Do 
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+(vu*Conjg(Tu(1,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vd*Conjg(Yu(1,1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(Conjg(Yu(1,1))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+(vu*Conjg(Tu(2,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vd*Conjg(Yu(2,1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(Conjg(Yu(2,1))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+(vu*Conjg(Tu(3,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vd*Conjg(Yu(3,1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,6) = mat(1,6)+(Conjg(Yu(3,1))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)-(g1**2*vd**2)/24._dp
mat(2,2) = mat(2,2)+(g2**2*vd**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vu**2)/24._dp
mat(2,2) = mat(2,2)-(g2**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+mq2(2,2)
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vu**2*Conjg(Yu(j1,2))*Yu(j1,2))/2._dp
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+mq2(2,3)
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu**2*Conjg(Yu(j1,2))*Yu(j1,3))/2._dp
End Do 
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(vu*Conjg(Tu(1,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(Yu(1,2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yu(1,2))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(vu*Conjg(Tu(2,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(Yu(2,2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yu(2,2))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,6) = 0._dp 
mat(2,6) = mat(2,6)+(vu*Conjg(Tu(3,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*Conjg(Yu(3,2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yu(3,2))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)-(g1**2*vd**2)/24._dp
mat(3,3) = mat(3,3)+(g2**2*vd**2)/8._dp
mat(3,3) = mat(3,3)+(g1**2*vu**2)/24._dp
mat(3,3) = mat(3,3)-(g2**2*vu**2)/8._dp
mat(3,3) = mat(3,3)+mq2(3,3)
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu**2*Conjg(Yu(j1,3))*Yu(j1,3))/2._dp
End Do 
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(vu*Conjg(Tu(1,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yu(1,3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yu(1,3))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(vu*Conjg(Tu(2,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yu(2,3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yu(2,3))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vu*Conjg(Tu(3,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yu(3,3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yu(3,3))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+(g1**2*vd**2)/6._dp
mat(4,4) = mat(4,4)-(g1**2*vu**2)/6._dp
mat(4,4) = mat(4,4)+mu2(1,1)
Do j1 = 1,3
mat(4,4) = mat(4,4)+(g1**2*vL(j1)**2)/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu**2*Conjg(Yu(1,j1))*Yu(1,j1))/2._dp
End Do 
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+mu2(1,2)
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yu(2,j1))*Yu(1,j1))/2._dp
End Do 
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+mu2(1,3)
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu**2*Conjg(Yu(3,j1))*Yu(1,j1))/2._dp
End Do 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+(g1**2*vd**2)/6._dp
mat(5,5) = mat(5,5)-(g1**2*vu**2)/6._dp
mat(5,5) = mat(5,5)+mu2(2,2)
Do j1 = 1,3
mat(5,5) = mat(5,5)+(g1**2*vL(j1)**2)/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu**2*Conjg(Yu(2,j1))*Yu(2,j1))/2._dp
End Do 
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+mu2(2,3)
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu**2*Conjg(Yu(3,j1))*Yu(2,j1))/2._dp
End Do 
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)+(g1**2*vd**2)/6._dp
mat(6,6) = mat(6,6)-(g1**2*vu**2)/6._dp
mat(6,6) = mat(6,6)+mu2(3,3)
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g1**2*vL(j1)**2)/6._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vu**2*Conjg(Yu(3,j1))*Yu(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,MSu2,ZU,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,6
  If (Abs(MSu2(i1)).Le.MaxMassNumericalZero) MSu2(i1) = 1.E-10_dp 
  If (MSu2(i1).ne.MSu2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return 
    End If 
  If (MSu2(i1).Ge.0._dp) Then 
  MSu(i1)=Sqrt(MSu2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MSu2(i1) 
    End If 
  MSu = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MSu2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MSu2(i1) 
  MSu2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMSu 

Subroutine CalculateMhh(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,mHu2,              & 
& mv2,vd,vu,vL,vR,ZH,Mhh,Mhh2,kont)

Real(dp), Intent(in) :: g1,g2,mlHd2(3),mHd2,mHu2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: lam(3),Tlam(3),Yv(3,3),Tv(3,3),kap(3,3,3),Tk(3,3,3),ml2(3,3),mv2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: Mhh(8), Mhh2(8) 
Real(dp), Intent(out) :: ZH(8,8) 
 
Real(dp) :: mat(8,8)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMhh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+mHd2
mat(1,1) = mat(1,1)+(3*g1**2*vd**2)/8._dp
mat(1,1) = mat(1,1)+(3*g2**2*vd**2)/8._dp
mat(1,1) = mat(1,1)-(g1**2*vu**2)/8._dp
mat(1,1) = mat(1,1)-(g2**2*vu**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(1,1) = mat(1,1)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vu**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g1**2*vd*vu)/4._dp
mat(1,2) = mat(1,2)-(g2**2*vd*vu)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(Tlam(j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)+vd*vu*Conjg(lam(j1))*lam(j1)
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vR(j1)*Tlam(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vu*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vu*Conjg(Yv(j2,j1))*vL(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)-(vu*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vd*Conjg(lam(j1))*vR(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vd*Conjg(lam(1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,3) = mat(1,3)-(vu*Tlam(1))/(2._dp*sqrt(2._dp))
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)-(vu*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vd*Conjg(lam(j1))*vR(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vd*Conjg(lam(2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,4) = mat(1,4)-(vu*Tlam(2))/(2._dp*sqrt(2._dp))
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)-(vu*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vd*Conjg(lam(j1))*vR(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vd*Conjg(lam(3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,5) = mat(1,5)-(vu*Tlam(3))/(2._dp*sqrt(2._dp))
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+mlHd2(1)
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(lam(j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(Yv(1,j1))*lam(j1))/4._dp
End Do 
mat(1,6) = mat(1,6)+(g1**2*vd*vL(1))/4._dp
mat(1,6) = mat(1,6)+(g2**2*vd*vL(1))/4._dp
mat(1,7) = 0._dp 
mat(1,7) = mat(1,7)+mlHd2(2)
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(lam(j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(Yv(2,j1))*lam(j1))/4._dp
End Do 
mat(1,7) = mat(1,7)+(g1**2*vd*vL(2))/4._dp
mat(1,7) = mat(1,7)+(g2**2*vd*vL(2))/4._dp
mat(1,8) = 0._dp 
mat(1,8) = mat(1,8)+mlHd2(3)
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(lam(j1))*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(Yv(3,j1))*lam(j1))/4._dp
End Do 
mat(1,8) = mat(1,8)+(g1**2*vd*vL(3))/4._dp
mat(1,8) = mat(1,8)+(g2**2*vd*vL(3))/4._dp
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mHu2
mat(2,2) = mat(2,2)-(g1**2*vd**2)/8._dp
mat(2,2) = mat(2,2)-(g2**2*vd**2)/8._dp
mat(2,2) = mat(2,2)+(3*g1**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+(3*g2**2*vu**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(2,2) = mat(2,2)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vd**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(Yv(j2,j1))*vL(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j1,j3))*vR(j2)*vR(j3)*Yv(j1,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j3,j1))*vL(j2)*vL(j3)*Yv(j2,j1))/2._dp
End Do 
End Do 
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-(vd*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Tv(j1,1))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(lam(j1))*vR(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(lam(1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vL(j1)*Tv(j1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(Yv(j1,j2))*vR(j2)*Yv(j1,1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(Yv(j1,1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(1,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(1,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j1,1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j1,j3,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j3,1,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j3,j1,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
End Do 
mat(2,3) = mat(2,3)-(vd*Tlam(1))/(2._dp*sqrt(2._dp))
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)-(vd*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Tv(j1,2))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(lam(j1))*vR(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(lam(2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vL(j1)*Tv(j1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(Yv(j1,j2))*vR(j2)*Yv(j1,2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(Yv(j1,2))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(2,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(2,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j1,2,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j1,j3,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j3,2,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j3,j1,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
End Do 
mat(2,4) = mat(2,4)-(vd*Tlam(2))/(2._dp*sqrt(2._dp))
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)-(vd*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Tv(j1,3))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(lam(j1))*vR(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(lam(3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vL(j1)*Tv(j1,3))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(Yv(j1,j2))*vR(j2)*Yv(j1,3))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(Yv(j1,3))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(3,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(3,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j1,3,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j1,j3,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j3,3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j3,j1,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
End Do 
mat(2,5) = mat(2,5)-(vd*Tlam(3))/(2._dp*sqrt(2._dp))
mat(2,6) = 0._dp 
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Tv(1,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*vu*Conjg(lam(j1))*Yv(1,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*vu*Conjg(Yv(1,j1))*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)+(vR(j1)*Tv(1,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(vu*Conjg(Yv(1,j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,6) = mat(2,6)-(g1**2*vu*vL(1))/4._dp
mat(2,6) = mat(2,6)-(g2**2*vu*vL(1))/4._dp
mat(2,7) = 0._dp 
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Tv(2,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vd*vu*Conjg(lam(j1))*Yv(2,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vd*vu*Conjg(Yv(2,j1))*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vR(j1)*Tv(2,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vu*Conjg(Yv(2,j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,7) = mat(2,7)-(g1**2*vu*vL(2))/4._dp
mat(2,7) = mat(2,7)-(g2**2*vu*vL(2))/4._dp
mat(2,8) = 0._dp 
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Tv(3,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vd*vu*Conjg(lam(j1))*Yv(3,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vd*vu*Conjg(Yv(3,j1))*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)+(vR(j1)*Tv(3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(vu*Conjg(Yv(3,j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,8) = mat(2,8)-(g1**2*vu*vL(3))/4._dp
mat(2,8) = mat(2,8)-(g2**2*vu*vL(3))/4._dp
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+mv2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Tk(1,1,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Tk(1,j1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Tk(j1,1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(lam(j1))*kap(1,1,j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(lam(j1))*kap(1,j1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(lam(j1))*kap(j1,1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(kap(1,1,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(kap(1,j1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(kap(j1,1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vR(j1)*Tk(1,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vR(j1)*Tk(1,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vR(j1)*Tk(j1,1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(kap(1,1,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(kap(1,j1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(kap(j1,1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,1,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,1))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(3,3) = mat(3,3)+(vd**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,3) = mat(3,3)+(vu**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+mv2(1,2)/2._dp
mat(3,4) = mat(3,4)+mv2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(1,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(1,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(2,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(2,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(j1,1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(j1,2,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(1,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(1,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(2,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(2,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(j1,1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(j1,2,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(1,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(1,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(2,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(2,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(j1,1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(j1,2,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(1,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(1,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(2,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(2,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(j1,1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(j1,2,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(1,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(1,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(2,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(2,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(j1,1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(j1,2,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+mv2(1,3)/2._dp
mat(3,5) = mat(3,5)+mv2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(1,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(1,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(3,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(3,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(j1,1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(j1,3,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(1,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(1,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(3,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(3,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(j1,1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(j1,3,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(1,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(1,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(3,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(3,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(j1,1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(j1,3,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(1,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(1,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(3,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(3,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(j1,1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(j1,3,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(1,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(1,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(3,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(3,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(j1,1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(j1,3,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vu*Conjg(Tv(1,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(lam(1))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yv(1,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yv(1,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(1,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,6) = mat(3,6)+(vu*Tv(1,1))/(2._dp*sqrt(2._dp))
mat(3,7) = 0._dp 
mat(3,7) = mat(3,7)+(vu*Conjg(Tv(2,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(lam(1))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(Yv(2,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(Yv(2,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(2,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,7) = mat(3,7)+(vu*Tv(2,1))/(2._dp*sqrt(2._dp))
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vu*Conjg(Tv(3,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(lam(1))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(Yv(3,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(Yv(3,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(3,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,8) = mat(3,8)+(vu*Tv(3,1))/(2._dp*sqrt(2._dp))
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+mv2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,2))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Tk(2,2,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Tk(2,j1,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Tk(j1,2,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(lam(j1))*kap(2,2,j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(lam(j1))*kap(2,j1,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(lam(j1))*kap(j1,2,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(kap(2,2,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(kap(2,j1,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(kap(j1,2,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vR(j1)*Tk(2,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vR(j1)*Tk(2,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vR(j1)*Tk(j1,2,2))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(kap(2,2,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(kap(2,j1,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(kap(j1,2,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,2))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,2))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(4,4) = mat(4,4)+(vd**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,4) = mat(4,4)+(vu**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+mv2(2,3)/2._dp
mat(4,5) = mat(4,5)+mv2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(2,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(2,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(3,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(3,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(j1,2,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(j1,3,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(2,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(2,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(3,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(3,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(j1,2,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(j1,3,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(2,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(2,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(3,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(3,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(j1,2,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(j1,3,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(2,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(2,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(3,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(3,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(j1,2,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(j1,3,2))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(2,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(2,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(3,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(3,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(j1,2,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(j1,3,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,2))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+(vu*Conjg(Tv(1,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(lam(2))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(Yv(1,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(Yv(1,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(1,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,6) = mat(4,6)+(vu*Tv(1,2))/(2._dp*sqrt(2._dp))
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(vu*Conjg(Tv(2,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(lam(2))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(Yv(2,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(Yv(2,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(2,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,7) = mat(4,7)+(vu*Tv(2,2))/(2._dp*sqrt(2._dp))
mat(4,8) = 0._dp 
mat(4,8) = mat(4,8)+(vu*Conjg(Tv(3,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(lam(2))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(Yv(3,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(Yv(3,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(3,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,8) = mat(4,8)+(vu*Tv(3,2))/(2._dp*sqrt(2._dp))
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+mv2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,3))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Tk(3,3,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Tk(3,j1,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Tk(j1,3,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(lam(j1))*kap(3,3,j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(lam(j1))*kap(3,j1,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(lam(j1))*kap(j1,3,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(kap(3,3,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(kap(3,j1,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(kap(j1,3,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*Tk(3,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*Tk(3,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*Tk(j1,3,3))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(kap(3,3,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(kap(3,j1,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(kap(j1,3,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,3,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,3))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,3))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(5,5) = mat(5,5)+(vd**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,5) = mat(5,5)+(vu**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(vu*Conjg(Tv(1,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(lam(3))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(Yv(1,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(Yv(1,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(1,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,6) = mat(5,6)+(vu*Tv(1,3))/(2._dp*sqrt(2._dp))
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)+(vu*Conjg(Tv(2,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(lam(3))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(Yv(2,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(Yv(2,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(2,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,7) = mat(5,7)+(vu*Tv(2,3))/(2._dp*sqrt(2._dp))
mat(5,8) = 0._dp 
mat(5,8) = mat(5,8)+(vu*Conjg(Tv(3,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(lam(3))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(Yv(3,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(Yv(3,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(3,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,8) = mat(5,8)+(vu*Tv(3,3))/(2._dp*sqrt(2._dp))
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)+(g1**2*vd**2)/8._dp
mat(6,6) = mat(6,6)+(g2**2*vd**2)/8._dp
mat(6,6) = mat(6,6)-(g1**2*vu**2)/8._dp
mat(6,6) = mat(6,6)-(g2**2*vu**2)/8._dp
mat(6,6) = mat(6,6)+ml2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(6,6) = mat(6,6)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vu**2*Conjg(Yv(1,j1))*Yv(1,j1))/2._dp
End Do 
mat(6,6) = mat(6,6)+(g1**2*vL(1)**2)/4._dp
mat(6,6) = mat(6,6)+(g2**2*vL(1)**2)/4._dp
mat(6,7) = 0._dp 
mat(6,7) = mat(6,7)+ml2(1,2)/2._dp
mat(6,7) = mat(6,7)+ml2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(2,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(1,j1))*Yv(2,j1))/4._dp
End Do 
mat(6,7) = mat(6,7)+(g1**2*vL(1)*vL(2))/4._dp
mat(6,7) = mat(6,7)+(g2**2*vL(1)*vL(2))/4._dp
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)+ml2(1,3)/2._dp
mat(6,8) = mat(6,8)+ml2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(3,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(1,j1))*Yv(3,j1))/4._dp
End Do 
mat(6,8) = mat(6,8)+(g1**2*vL(1)*vL(3))/4._dp
mat(6,8) = mat(6,8)+(g2**2*vL(1)*vL(3))/4._dp
mat(7,7) = 0._dp 
mat(7,7) = mat(7,7)+(g1**2*vd**2)/8._dp
mat(7,7) = mat(7,7)+(g2**2*vd**2)/8._dp
mat(7,7) = mat(7,7)-(g1**2*vu**2)/8._dp
mat(7,7) = mat(7,7)-(g2**2*vu**2)/8._dp
mat(7,7) = mat(7,7)+ml2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(7,7) = mat(7,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(vu**2*Conjg(Yv(2,j1))*Yv(2,j1))/2._dp
End Do 
mat(7,7) = mat(7,7)+(g1**2*vL(2)**2)/4._dp
mat(7,7) = mat(7,7)+(g2**2*vL(2)**2)/4._dp
mat(7,8) = 0._dp 
mat(7,8) = mat(7,8)+ml2(2,3)/2._dp
mat(7,8) = mat(7,8)+ml2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(3,j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(2,j1))*Yv(3,j1))/4._dp
End Do 
mat(7,8) = mat(7,8)+(g1**2*vL(2)*vL(3))/4._dp
mat(7,8) = mat(7,8)+(g2**2*vL(2)*vL(3))/4._dp
mat(8,8) = 0._dp 
mat(8,8) = mat(8,8)+(g1**2*vd**2)/8._dp
mat(8,8) = mat(8,8)+(g2**2*vd**2)/8._dp
mat(8,8) = mat(8,8)-(g1**2*vu**2)/8._dp
mat(8,8) = mat(8,8)-(g2**2*vu**2)/8._dp
mat(8,8) = mat(8,8)+ml2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(8,8) = mat(8,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(vu**2*Conjg(Yv(3,j1))*Yv(3,j1))/2._dp
End Do 
mat(8,8) = mat(8,8)+(g1**2*vL(3)**2)/4._dp
mat(8,8) = mat(8,8)+(g2**2*vL(3)**2)/4._dp

 
 Do i1=2,8
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,Mhh2,ZH,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,8
  If (Abs(Mhh2(i1)).Le.MaxMassNumericalZero) Mhh2(i1) = 1.E-10_dp 
  If (Mhh2(i1).ne.Mhh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Mhh2(i1).Ge.0._dp) Then 
  Mhh(i1)=Sqrt(Mhh2(i1) ) 
!  Write(*,*)  'Good ', i1,Mhh2(i1), ErrorLevel, mHd2, mHu2
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Mhh2(i1) 
    End If 
  Mhh = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Mhh2(i1) 
    Write(*,*)  'Bad ', i1,Mhh2(i1), ErrorLevel, mHd2, mHu2
    ! Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
    ! Write(*,*) 'in the calculation of the masses' 
    ! Write(*,*) 'occurred a negative mass squared!' 
    ! Write(*,*) i1,Mhh2(i1) 
  Mhh2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMhh 

Subroutine CalculateMAh(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,mHu2,              & 
& mv2,vd,vu,vL,vR,TW,ZA,MAh,MAh2,kont)

Real(dp), Intent(in) :: g1,g2,mlHd2(3),mHd2,mHu2,vd,vu,vL(3),vR(3),TW

Complex(dp), Intent(in) :: lam(3),Tlam(3),Yv(3,3),Tv(3,3),kap(3,3,3),Tk(3,3,3),ml2(3,3),mv2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MAh(8), MAh2(8) 
Real(dp), Intent(out) :: ZA(8,8) 
 
Real(dp) :: mat(8,8)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMAh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+mHd2
mat(1,1) = mat(1,1)+(g1**2*vd**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vd**2)/8._dp
mat(1,1) = mat(1,1)-(g1**2*vu**2)/8._dp
mat(1,1) = mat(1,1)-(g2**2*vu**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vd**2*Cos(TW)**2*RXiZ)/4._dp
mat(1,1) = mat(1,1)+(g1*g2*vd**2*Cos(TW)*RXiZ*Sin(TW))/2._dp
mat(1,1) = mat(1,1)+(g1**2*vd**2*RXiZ*Sin(TW)**2)/4._dp
Do j1 = 1,3
Do j2 = 1,3
mat(1,1) = mat(1,1)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vu**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g2**2*vd*vu*Cos(TW)**2*RXiZ)/4._dp
mat(1,2) = mat(1,2)-(g1*g2*vd*vu*Cos(TW)*RXiZ*Sin(TW))/2._dp
mat(1,2) = mat(1,2)-(g1**2*vd*vu*RXiZ*Sin(TW)**2)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(Tlam(j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vR(j1)*Tlam(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)-(vu*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)+(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,3) = mat(1,3)-(vu*Tlam(1))/(2._dp*sqrt(2._dp))
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)-(vu*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)+(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,4) = mat(1,4)-(vu*Tlam(2))/(2._dp*sqrt(2._dp))
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)-(vu*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)+(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,5) = mat(1,5)-(vu*Tlam(3))/(2._dp*sqrt(2._dp))
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+mlHd2(1)
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(lam(j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(Yv(1,j1))*lam(j1))/4._dp
End Do 
mat(1,6) = mat(1,6)+(g2**2*vd*Cos(TW)**2*RXiZ*vL(1))/4._dp
mat(1,6) = mat(1,6)+(g1*g2*vd*Cos(TW)*RXiZ*Sin(TW)*vL(1))/2._dp
mat(1,6) = mat(1,6)+(g1**2*vd*RXiZ*Sin(TW)**2*vL(1))/4._dp
mat(1,7) = 0._dp 
mat(1,7) = mat(1,7)+mlHd2(2)
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(lam(j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(Yv(2,j1))*lam(j1))/4._dp
End Do 
mat(1,7) = mat(1,7)+(g2**2*vd*Cos(TW)**2*RXiZ*vL(2))/4._dp
mat(1,7) = mat(1,7)+(g1*g2*vd*Cos(TW)*RXiZ*Sin(TW)*vL(2))/2._dp
mat(1,7) = mat(1,7)+(g1**2*vd*RXiZ*Sin(TW)**2*vL(2))/4._dp
mat(1,8) = 0._dp 
mat(1,8) = mat(1,8)+mlHd2(3)
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(lam(j1))*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(Yv(3,j1))*lam(j1))/4._dp
End Do 
mat(1,8) = mat(1,8)+(g2**2*vd*Cos(TW)**2*RXiZ*vL(3))/4._dp
mat(1,8) = mat(1,8)+(g1*g2*vd*Cos(TW)*RXiZ*Sin(TW)*vL(3))/2._dp
mat(1,8) = mat(1,8)+(g1**2*vd*RXiZ*Sin(TW)**2*vL(3))/4._dp
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mHu2
mat(2,2) = mat(2,2)-(g1**2*vd**2)/8._dp
mat(2,2) = mat(2,2)-(g2**2*vd**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2*Cos(TW)**2*RXiZ)/4._dp
mat(2,2) = mat(2,2)+(g1*g2*vu**2*Cos(TW)*RXiZ*Sin(TW))/2._dp
mat(2,2) = mat(2,2)+(g1**2*vu**2*RXiZ*Sin(TW)**2)/4._dp
Do j1 = 1,3
Do j2 = 1,3
mat(2,2) = mat(2,2)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vd**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(Yv(j2,j1))*vL(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j1,j3))*vR(j2)*vR(j3)*Yv(j1,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j3,j1))*vL(j2)*vL(j3)*Yv(j2,j1))/2._dp
End Do 
End Do 
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-(vd*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Tv(j1,1))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vL(j1)*Tv(j1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(1,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(1,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j1,1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j1,j3,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j3,1,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j3,j1,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
End Do 
mat(2,3) = mat(2,3)-(vd*Tlam(1))/(2._dp*sqrt(2._dp))
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)-(vd*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Tv(j1,2))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vL(j1)*Tv(j1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(2,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(2,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j1,2,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j1,j3,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j3,2,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j3,j1,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
End Do 
mat(2,4) = mat(2,4)-(vd*Tlam(2))/(2._dp*sqrt(2._dp))
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)-(vd*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Tv(j1,3))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vL(j1)*Tv(j1,3))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(3,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(3,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j1,3,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j1,j3,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j3,3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j3,j1,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
End Do 
mat(2,5) = mat(2,5)-(vd*Tlam(3))/(2._dp*sqrt(2._dp))
mat(2,6) = 0._dp 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Tv(1,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vR(j1)*Tv(1,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,6) = mat(2,6)-(g2**2*vu*Cos(TW)**2*RXiZ*vL(1))/4._dp
mat(2,6) = mat(2,6)-(g1*g2*vu*Cos(TW)*RXiZ*Sin(TW)*vL(1))/2._dp
mat(2,6) = mat(2,6)-(g1**2*vu*RXiZ*Sin(TW)**2*vL(1))/4._dp
mat(2,7) = 0._dp 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Tv(2,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vR(j1)*Tv(2,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,7) = mat(2,7)-(g2**2*vu*Cos(TW)**2*RXiZ*vL(2))/4._dp
mat(2,7) = mat(2,7)-(g1*g2*vu*Cos(TW)*RXiZ*Sin(TW)*vL(2))/2._dp
mat(2,7) = mat(2,7)-(g1**2*vu*RXiZ*Sin(TW)**2*vL(2))/4._dp
mat(2,8) = 0._dp 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Tv(3,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vR(j1)*Tv(3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,8) = mat(2,8)-(g2**2*vu*Cos(TW)**2*RXiZ*vL(3))/4._dp
mat(2,8) = mat(2,8)-(g1*g2*vu*Cos(TW)*RXiZ*Sin(TW)*vL(3))/2._dp
mat(2,8) = mat(2,8)-(g1**2*vu*RXiZ*Sin(TW)**2*vL(3))/4._dp
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+mv2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(Tk(1,1,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(Tk(1,j1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(Tk(j1,1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(lam(j1))*kap(1,1,j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(lam(j1))*kap(1,j1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(lam(j1))*kap(j1,1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(kap(1,1,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(kap(1,j1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(kap(j1,1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vR(j1)*Tk(1,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vR(j1)*Tk(1,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vR(j1)*Tk(j1,1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(kap(1,1,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(kap(1,j1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(kap(j1,1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,1,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,1))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(3,3) = mat(3,3)+(vd**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,3) = mat(3,3)+(vu**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+mv2(1,2)/2._dp
mat(3,4) = mat(3,4)+mv2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(1,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(1,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(2,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(2,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(j1,1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(j1,2,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(1,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(1,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(2,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(2,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(j1,1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(j1,2,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(1,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(1,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(2,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(2,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(j1,1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(j1,2,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(1,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(1,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(2,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(2,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(j1,1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(j1,2,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(1,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(1,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(2,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(2,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(j1,1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(j1,2,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+mv2(1,3)/2._dp
mat(3,5) = mat(3,5)+mv2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(1,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(1,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(3,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(3,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(j1,1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(j1,3,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(1,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(1,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(3,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(3,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(j1,1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(j1,3,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(1,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(1,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(3,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(3,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(j1,1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(j1,3,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(1,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(1,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(3,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(3,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(j1,1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(j1,3,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(1,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(1,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(3,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(3,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(j1,1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(j1,3,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vu*Conjg(Tv(1,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)-(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vd*Conjg(lam(1))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vd*Conjg(Yv(1,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yv(1,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(1,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,6) = mat(3,6)+(vu*Tv(1,1))/(2._dp*sqrt(2._dp))
mat(3,7) = 0._dp 
mat(3,7) = mat(3,7)+(vu*Conjg(Tv(2,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)-(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vd*Conjg(lam(1))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vd*Conjg(Yv(2,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(Yv(2,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(2,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,7) = mat(3,7)+(vu*Tv(2,1))/(2._dp*sqrt(2._dp))
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vu*Conjg(Tv(3,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)-(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vd*Conjg(lam(1))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vd*Conjg(Yv(3,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(Yv(3,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(3,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,8) = mat(3,8)+(vu*Tv(3,1))/(2._dp*sqrt(2._dp))
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+mv2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,2))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(Tk(2,2,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(Tk(2,j1,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(Tk(j1,2,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(lam(j1))*kap(2,2,j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(lam(j1))*kap(2,j1,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(lam(j1))*kap(j1,2,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(kap(2,2,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(kap(2,j1,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(kap(j1,2,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vR(j1)*Tk(2,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vR(j1)*Tk(2,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vR(j1)*Tk(j1,2,2))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(kap(2,2,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(kap(2,j1,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(kap(j1,2,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,2))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,2))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(4,4) = mat(4,4)+(vd**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,4) = mat(4,4)+(vu**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+mv2(2,3)/2._dp
mat(4,5) = mat(4,5)+mv2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(2,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(2,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(3,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(3,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(j1,2,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(j1,3,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(2,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(2,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(3,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(3,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(j1,2,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(j1,3,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(2,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(2,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(3,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(3,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(j1,2,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(j1,3,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(2,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(2,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(3,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(3,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(j1,2,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(j1,3,2))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(2,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(2,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(3,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(3,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(j1,2,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(j1,3,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,2))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+(vu*Conjg(Tv(1,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)-(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vd*Conjg(lam(2))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vd*Conjg(Yv(1,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(Yv(1,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(1,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,6) = mat(4,6)+(vu*Tv(1,2))/(2._dp*sqrt(2._dp))
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(vu*Conjg(Tv(2,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)-(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vd*Conjg(lam(2))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vd*Conjg(Yv(2,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(Yv(2,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(2,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,7) = mat(4,7)+(vu*Tv(2,2))/(2._dp*sqrt(2._dp))
mat(4,8) = 0._dp 
mat(4,8) = mat(4,8)+(vu*Conjg(Tv(3,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)-(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vd*Conjg(lam(2))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vd*Conjg(Yv(3,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(Yv(3,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(3,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,8) = mat(4,8)+(vu*Tv(3,2))/(2._dp*sqrt(2._dp))
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+mv2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,3))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(Tk(3,3,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(Tk(3,j1,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(Tk(j1,3,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(lam(j1))*kap(3,3,j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(lam(j1))*kap(3,j1,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(lam(j1))*kap(j1,3,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(kap(3,3,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(kap(3,j1,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(kap(j1,3,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vR(j1)*Tk(3,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vR(j1)*Tk(3,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vR(j1)*Tk(j1,3,3))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(kap(3,3,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(kap(3,j1,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(kap(j1,3,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,3,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,3))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,3))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(5,5) = mat(5,5)+(vd**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,5) = mat(5,5)+(vu**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(vu*Conjg(Tv(1,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)-(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vd*Conjg(lam(3))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vd*Conjg(Yv(1,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(Yv(1,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(1,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,6) = mat(5,6)+(vu*Tv(1,3))/(2._dp*sqrt(2._dp))
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)+(vu*Conjg(Tv(2,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)-(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vd*Conjg(lam(3))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vd*Conjg(Yv(2,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(Yv(2,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(2,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,7) = mat(5,7)+(vu*Tv(2,3))/(2._dp*sqrt(2._dp))
mat(5,8) = 0._dp 
mat(5,8) = mat(5,8)+(vu*Conjg(Tv(3,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)-(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vd*Conjg(lam(3))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vd*Conjg(Yv(3,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(Yv(3,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(3,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,8) = mat(5,8)+(vu*Tv(3,3))/(2._dp*sqrt(2._dp))
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)+(g1**2*vd**2)/8._dp
mat(6,6) = mat(6,6)+(g2**2*vd**2)/8._dp
mat(6,6) = mat(6,6)-(g1**2*vu**2)/8._dp
mat(6,6) = mat(6,6)-(g2**2*vu**2)/8._dp
mat(6,6) = mat(6,6)+ml2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(6,6) = mat(6,6)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vu**2*Conjg(Yv(1,j1))*Yv(1,j1))/2._dp
End Do 
mat(6,6) = mat(6,6)+(g2**2*Cos(TW)**2*RXiZ*vL(1)**2)/4._dp
mat(6,6) = mat(6,6)+(g1*g2*Cos(TW)*RXiZ*Sin(TW)*vL(1)**2)/2._dp
mat(6,6) = mat(6,6)+(g1**2*RXiZ*Sin(TW)**2*vL(1)**2)/4._dp
mat(6,7) = 0._dp 
mat(6,7) = mat(6,7)+ml2(1,2)/2._dp
mat(6,7) = mat(6,7)+ml2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(2,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(1,j1))*Yv(2,j1))/4._dp
End Do 
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)+ml2(1,3)/2._dp
mat(6,8) = mat(6,8)+ml2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(3,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(1,j1))*Yv(3,j1))/4._dp
End Do 
mat(7,7) = 0._dp 
mat(7,7) = mat(7,7)+(g1**2*vd**2)/8._dp
mat(7,7) = mat(7,7)+(g2**2*vd**2)/8._dp
mat(7,7) = mat(7,7)-(g1**2*vu**2)/8._dp
mat(7,7) = mat(7,7)-(g2**2*vu**2)/8._dp
mat(7,7) = mat(7,7)+ml2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(7,7) = mat(7,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(vu**2*Conjg(Yv(2,j1))*Yv(2,j1))/2._dp
End Do 
mat(7,7) = mat(7,7)+(g2**2*Cos(TW)**2*RXiZ*vL(2)**2)/4._dp
mat(7,7) = mat(7,7)+(g1*g2*Cos(TW)*RXiZ*Sin(TW)*vL(2)**2)/2._dp
mat(7,7) = mat(7,7)+(g1**2*RXiZ*Sin(TW)**2*vL(2)**2)/4._dp
mat(7,8) = 0._dp 
mat(7,8) = mat(7,8)+ml2(2,3)/2._dp
mat(7,8) = mat(7,8)+ml2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(3,j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(2,j1))*Yv(3,j1))/4._dp
End Do 
mat(8,8) = 0._dp 
mat(8,8) = mat(8,8)+(g1**2*vd**2)/8._dp
mat(8,8) = mat(8,8)+(g2**2*vd**2)/8._dp
mat(8,8) = mat(8,8)-(g1**2*vu**2)/8._dp
mat(8,8) = mat(8,8)-(g2**2*vu**2)/8._dp
mat(8,8) = mat(8,8)+ml2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(8,8) = mat(8,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(vu**2*Conjg(Yv(3,j1))*Yv(3,j1))/2._dp
End Do 
mat(8,8) = mat(8,8)+(g2**2*Cos(TW)**2*RXiZ*vL(3)**2)/4._dp
mat(8,8) = mat(8,8)+(g1*g2*Cos(TW)*RXiZ*Sin(TW)*vL(3)**2)/2._dp
mat(8,8) = mat(8,8)+(g1**2*RXiZ*Sin(TW)**2*vL(3)**2)/4._dp

 
 Do i1=2,8
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MAh2,ZA,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,8
  If (Abs(MAh2(i1)).Le.MaxMassNumericalZero) MAh2(i1) = 1.E-10_dp 
  If (MAh2(i1).ne.MAh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (MAh2(i1).Ge.0._dp) Then 
  MAh(i1)=Sqrt(MAh2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MAh2(i1) 
    End If 
  MAh = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MAh2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MAh2(i1) 
  MAh2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMAh 

Subroutine CalculateMHpm(g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,ml2,mlHd2,mHd2,               & 
& mHu2,me2,vd,vu,vL,vR,ZP,MHpm,MHpm2,kont)

Real(dp), Intent(in) :: g1,g2,mlHd2(3),mHd2,mHu2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: Ye(3,3),Te(3,3),lam(3),Tlam(3),Yv(3,3),Tv(3,3),kap(3,3,3),ml2(3,3),me2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MHpm(8), MHpm2(8) 
Real(dp), Intent(out) :: ZP(8,8) 
 
Real(dp) :: mat(8,8)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMHpm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+mHd2
mat(1,1) = mat(1,1)+(g1**2*vd**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vd**2)/8._dp
mat(1,1) = mat(1,1)-(g1**2*vu**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vu**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vd**2*RXiWm)/4._dp
Do j1 = 1,3
Do j2 = 1,3
mat(1,1) = mat(1,1)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,1) = mat(1,1)+(Conjg(Ye(j1,j3))*vL(j2)*vL(j3)*Ye(j1,j2))/2._dp
End Do 
End Do 
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(g2**2*vd*vu)/4._dp
mat(1,2) = mat(1,2)-(g2**2*vd*vu*RXiWm)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(Tlam(j1))*vR(j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vd*vu*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vu*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/6._dp
End Do 
End Do 
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+mlHd2(1)
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vd*Conjg(Ye(j1,j2))*vL(j2)*Ye(j1,1))/2._dp
End Do 
End Do 
mat(1,3) = mat(1,3)+(g2**2*vd*vL(1))/4._dp
mat(1,3) = mat(1,3)+(g2**2*vd*RXiWm*vL(1))/4._dp
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+mlHd2(2)
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vd*Conjg(Ye(j1,j2))*vL(j2)*Ye(j1,2))/2._dp
End Do 
End Do 
mat(1,4) = mat(1,4)+(g2**2*vd*vL(2))/4._dp
mat(1,4) = mat(1,4)+(g2**2*vd*RXiWm*vL(2))/4._dp
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+mlHd2(3)
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vd*Conjg(Ye(j1,j2))*vL(j2)*Ye(j1,3))/2._dp
End Do 
End Do 
mat(1,5) = mat(1,5)+(g2**2*vd*vL(3))/4._dp
mat(1,5) = mat(1,5)+(g2**2*vd*RXiWm*vL(3))/4._dp
mat(1,6) = 0._dp 
Do j1 = 1,3
mat(1,6) = mat(1,6)-((Conjg(Te(1,j1))*vL(j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu*Conjg(Ye(1,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,7) = 0._dp 
Do j1 = 1,3
mat(1,7) = mat(1,7)-((Conjg(Te(2,j1))*vL(j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu*Conjg(Ye(2,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,8) = 0._dp 
Do j1 = 1,3
mat(1,8) = mat(1,8)-((Conjg(Te(3,j1))*vL(j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu*Conjg(Ye(3,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mHu2
mat(2,2) = mat(2,2)-(g1**2*vd**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vd**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2*RXiWm)/4._dp
Do j1 = 1,3
Do j2 = 1,3
mat(2,2) = mat(2,2)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j1,j3))*vR(j2)*vR(j3)*Yv(j1,j2))/2._dp
End Do 
End Do 
End Do 
mat(2,3) = 0._dp 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*vu*Conjg(lam(j1))*Yv(1,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)-((vR(j1)*Tv(1,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(1,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(1,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(1,j1))/6._dp
End Do 
End Do 
End Do 
mat(2,3) = mat(2,3)+(g2**2*vu*vL(1))/4._dp
mat(2,3) = mat(2,3)-(g2**2*vu*RXiWm*vL(1))/4._dp
mat(2,4) = 0._dp 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*vu*Conjg(lam(j1))*Yv(2,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)-((vR(j1)*Tv(2,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(2,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(2,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(2,j1))/6._dp
End Do 
End Do 
End Do 
mat(2,4) = mat(2,4)+(g2**2*vu*vL(2))/4._dp
mat(2,4) = mat(2,4)-(g2**2*vu*RXiWm*vL(2))/4._dp
mat(2,5) = 0._dp 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*vu*Conjg(lam(j1))*Yv(3,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)-((vR(j1)*Tv(3,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(3,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(3,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(3,j1))/6._dp
End Do 
End Do 
End Do 
mat(2,5) = mat(2,5)+(g2**2*vu*vL(3))/4._dp
mat(2,5) = mat(2,5)-(g2**2*vu*RXiWm*vL(3))/4._dp
mat(2,6) = 0._dp 
Do j1 = 1,3
Do j2 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Ye(1,j2))*vL(j2)*vR(j1)*lam(j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*Conjg(Ye(1,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,7) = 0._dp 
Do j1 = 1,3
Do j2 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Ye(2,j2))*vL(j2)*vR(j1)*lam(j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vd*Conjg(Ye(2,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,8) = 0._dp 
Do j1 = 1,3
Do j2 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Ye(3,j2))*vL(j2)*vR(j1)*lam(j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vd*Conjg(Ye(3,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(g1**2*vd**2)/8._dp
mat(3,3) = mat(3,3)-(g2**2*vd**2)/8._dp
mat(3,3) = mat(3,3)-(g1**2*vu**2)/8._dp
mat(3,3) = mat(3,3)+(g2**2*vu**2)/8._dp
mat(3,3) = mat(3,3)+ml2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,1))/2._dp
End Do 
mat(3,3) = mat(3,3)+(g2**2*vL(1)**2)/4._dp
mat(3,3) = mat(3,3)+(g2**2*RXiWm*vL(1)**2)/4._dp
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+ml2(1,2)
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,2))/2._dp
End Do 
mat(3,4) = mat(3,4)+(g2**2*vL(1)*vL(2))/4._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+ml2(1,3)
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,3))/2._dp
End Do 
mat(3,5) = mat(3,5)+(g2**2*vL(1)*vL(3))/4._dp
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vd*Conjg(Te(1,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Ye(1,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,7) = 0._dp 
mat(3,7) = mat(3,7)+(vd*Conjg(Te(2,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Ye(2,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vd*Conjg(Te(3,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Ye(3,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+(g1**2*vd**2)/8._dp
mat(4,4) = mat(4,4)-(g2**2*vd**2)/8._dp
mat(4,4) = mat(4,4)-(g1**2*vu**2)/8._dp
mat(4,4) = mat(4,4)+(g2**2*vu**2)/8._dp
mat(4,4) = mat(4,4)+ml2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd**2*Conjg(Ye(j1,2))*Ye(j1,2))/2._dp
End Do 
mat(4,4) = mat(4,4)+(g2**2*vL(2)**2)/4._dp
mat(4,4) = mat(4,4)+(g2**2*RXiWm*vL(2)**2)/4._dp
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+ml2(2,3)
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd**2*Conjg(Ye(j1,2))*Ye(j1,3))/2._dp
End Do 
mat(4,5) = mat(4,5)+(g2**2*vL(2)*vL(3))/4._dp
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+(vd*Conjg(Te(1,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Ye(1,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(vd*Conjg(Te(2,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Ye(2,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,8) = 0._dp 
mat(4,8) = mat(4,8)+(vd*Conjg(Te(3,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Ye(3,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+(g1**2*vd**2)/8._dp
mat(5,5) = mat(5,5)-(g2**2*vd**2)/8._dp
mat(5,5) = mat(5,5)-(g1**2*vu**2)/8._dp
mat(5,5) = mat(5,5)+(g2**2*vu**2)/8._dp
mat(5,5) = mat(5,5)+ml2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd**2*Conjg(Ye(j1,3))*Ye(j1,3))/2._dp
End Do 
mat(5,5) = mat(5,5)+(g2**2*vL(3)**2)/4._dp
mat(5,5) = mat(5,5)+(g2**2*RXiWm*vL(3)**2)/4._dp
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(vd*Conjg(Te(1,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Ye(1,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)+(vd*Conjg(Te(2,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Ye(2,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(5,8) = 0._dp 
mat(5,8) = mat(5,8)+(vd*Conjg(Te(3,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Ye(3,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)-(g1**2*vd**2)/4._dp
mat(6,6) = mat(6,6)+(g1**2*vu**2)/4._dp
mat(6,6) = mat(6,6)+me2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(6,6) = mat(6,6)+(Conjg(Ye(1,j2))*vL(j1)*vL(j2)*Ye(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)-(g1**2*vL(j1)**2)/4._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vd**2*Conjg(Ye(1,j1))*Ye(1,j1))/2._dp
End Do 
mat(6,7) = 0._dp 
mat(6,7) = mat(6,7)+me2(1,2)
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Ye(2,j2))*vL(j1)*vL(j2)*Ye(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vd**2*Conjg(Ye(2,j1))*Ye(1,j1))/2._dp
End Do 
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)+me2(1,3)
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Ye(3,j2))*vL(j1)*vL(j2)*Ye(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vd**2*Conjg(Ye(3,j1))*Ye(1,j1))/2._dp
End Do 
mat(7,7) = 0._dp 
mat(7,7) = mat(7,7)-(g1**2*vd**2)/4._dp
mat(7,7) = mat(7,7)+(g1**2*vu**2)/4._dp
mat(7,7) = mat(7,7)+me2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(7,7) = mat(7,7)+(Conjg(Ye(2,j2))*vL(j1)*vL(j2)*Ye(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)-(g1**2*vL(j1)**2)/4._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(vd**2*Conjg(Ye(2,j1))*Ye(2,j1))/2._dp
End Do 
mat(7,8) = 0._dp 
mat(7,8) = mat(7,8)+me2(2,3)
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Ye(3,j2))*vL(j1)*vL(j2)*Ye(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vd**2*Conjg(Ye(3,j1))*Ye(2,j1))/2._dp
End Do 
mat(8,8) = 0._dp 
mat(8,8) = mat(8,8)-(g1**2*vd**2)/4._dp
mat(8,8) = mat(8,8)+(g1**2*vu**2)/4._dp
mat(8,8) = mat(8,8)+me2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(8,8) = mat(8,8)+(Conjg(Ye(3,j2))*vL(j1)*vL(j2)*Ye(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)-(g1**2*vL(j1)**2)/4._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(vd**2*Conjg(Ye(3,j1))*Ye(3,j1))/2._dp
End Do 

Do i1=2,8
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
!  print*,i1,mat(i1,:)
End do 
 
Call EigenSystem(mat,MHpm2,ZP,ierr,test) 
  
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 

Do i1=1,8
  If (Abs(MHpm2(i1)).Le.MaxMassNumericalZero) MHpm2(i1) = 1.E-10_dp 
  If (MHpm2(i1).ne.MHpm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (MHpm2(i1).Ge.0._dp) Then 
  MHpm(i1)=Sqrt(MHpm2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MHpm2(i1) 
    End If 
  MHpm = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MHpm2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MHpm2(i1) 
  MHpm2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMHpm 

Subroutine CalculateMChi(g1,g2,lam,Yv,kap,M1,M2,vd,vu,vL,vR,UV,MChi,kont)

Real(dp) ,Intent(in) :: g1,g2,vd,vu,vL(3),vR(3)

Complex(dp) ,Intent(in) :: lam(3),Yv(3,3),kap(3,3,3),M1,M2

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MChi(10) 
Complex(dp), Intent(out) ::  UV(10,10) 
                              
Complex(dp) :: mat(10,10), mat2(10,10), phaseM, E10(10) 

Real(dp) :: UVa(10,10), test(2), eig(10) 

Real(dp) :: MChitemp(10) 
Complex(dp) :: UVtemp(10,10) 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMChi'
 
mat(1,1) = 0._dp 
mat(1,2) = 0._dp 
mat(1,3) = 0._dp 
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)-(g1*vL(1))/2._dp
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+(g2*vL(1))/2._dp
mat(1,6) = 0._dp 
mat(1,7) = 0._dp 
Do j1 = 1,3
mat(1,7) = mat(1,7)+(vR(j1)*Yv(1,j1))/sqrt(2._dp)
End Do 
mat(1,8) = 0._dp 
mat(1,8) = mat(1,8)+(vu*Yv(1,1))/sqrt(2._dp)
mat(1,9) = 0._dp 
mat(1,9) = mat(1,9)+(vu*Yv(1,2))/sqrt(2._dp)
mat(1,10) = 0._dp 
mat(1,10) = mat(1,10)+(vu*Yv(1,3))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,3) = 0._dp 
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)-(g1*vL(2))/2._dp
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(g2*vL(2))/2._dp
mat(2,6) = 0._dp 
mat(2,7) = 0._dp 
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vR(j1)*Yv(2,j1))/sqrt(2._dp)
End Do 
mat(2,8) = 0._dp 
mat(2,8) = mat(2,8)+(vu*Yv(2,1))/sqrt(2._dp)
mat(2,9) = 0._dp 
mat(2,9) = mat(2,9)+(vu*Yv(2,2))/sqrt(2._dp)
mat(2,10) = 0._dp 
mat(2,10) = mat(2,10)+(vu*Yv(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)-(g1*vL(3))/2._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(g2*vL(3))/2._dp
mat(3,6) = 0._dp 
mat(3,7) = 0._dp 
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vR(j1)*Yv(3,j1))/sqrt(2._dp)
End Do 
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vu*Yv(3,1))/sqrt(2._dp)
mat(3,9) = 0._dp 
mat(3,9) = mat(3,9)+(vu*Yv(3,2))/sqrt(2._dp)
mat(3,10) = 0._dp 
mat(3,10) = mat(3,10)+(vu*Yv(3,3))/sqrt(2._dp)
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+M1
mat(4,5) = 0._dp 
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)-(g1*vd)/2._dp
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(g1*vu)/2._dp
mat(4,8) = 0._dp 
mat(4,9) = 0._dp 
mat(4,10) = 0._dp 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+M2
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(g2*vd)/2._dp
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)-(g2*vu)/2._dp
mat(5,8) = 0._dp 
mat(5,9) = 0._dp 
mat(5,10) = 0._dp 
mat(6,6) = 0._dp 
mat(6,7) = 0._dp 
Do j1 = 1,3
mat(6,7) = mat(6,7)-((vR(j1)*lam(j1))/sqrt(2._dp))
End Do 
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)-((vu*lam(1))/sqrt(2._dp))
mat(6,9) = 0._dp 
mat(6,9) = mat(6,9)-((vu*lam(2))/sqrt(2._dp))
mat(6,10) = 0._dp 
mat(6,10) = mat(6,10)-((vu*lam(3))/sqrt(2._dp))
mat(7,7) = 0._dp 
mat(7,8) = 0._dp 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vL(j1)*Yv(j1,1))/sqrt(2._dp)
End Do 
mat(7,8) = mat(7,8)-((vd*lam(1))/sqrt(2._dp))
mat(7,9) = 0._dp 
Do j1 = 1,3
mat(7,9) = mat(7,9)+(vL(j1)*Yv(j1,2))/sqrt(2._dp)
End Do 
mat(7,9) = mat(7,9)-((vd*lam(2))/sqrt(2._dp))
mat(7,10) = 0._dp 
Do j1 = 1,3
mat(7,10) = mat(7,10)+(vL(j1)*Yv(j1,3))/sqrt(2._dp)
End Do 
mat(7,10) = mat(7,10)-((vd*lam(3))/sqrt(2._dp))
mat(8,8) = 0._dp 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(sqrt(2._dp)*vR(j1)*kap(1,1,j1))/3._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(sqrt(2._dp)*vR(j1)*kap(1,j1,1))/3._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(sqrt(2._dp)*vR(j1)*kap(j1,1,1))/3._dp
End Do 
mat(8,9) = 0._dp 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(1,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(1,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(2,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(2,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(j1,1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(j1,2,1))/(3._dp*sqrt(2._dp))
End Do 
mat(8,10) = 0._dp 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(1,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(1,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(3,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(3,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(j1,1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(j1,3,1))/(3._dp*sqrt(2._dp))
End Do 
mat(9,9) = 0._dp 
Do j1 = 1,3
mat(9,9) = mat(9,9)+(sqrt(2._dp)*vR(j1)*kap(2,2,j1))/3._dp
End Do 
Do j1 = 1,3
mat(9,9) = mat(9,9)+(sqrt(2._dp)*vR(j1)*kap(2,j1,2))/3._dp
End Do 
Do j1 = 1,3
mat(9,9) = mat(9,9)+(sqrt(2._dp)*vR(j1)*kap(j1,2,2))/3._dp
End Do 
mat(9,10) = 0._dp 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(2,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(2,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(3,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(3,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(j1,2,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(j1,3,2))/(3._dp*sqrt(2._dp))
End Do 
mat(10,10) = 0._dp 
Do j1 = 1,3
mat(10,10) = mat(10,10)+(sqrt(2._dp)*vR(j1)*kap(3,3,j1))/3._dp
End Do 
Do j1 = 1,3
mat(10,10) = mat(10,10)+(sqrt(2._dp)*vR(j1)*kap(3,j1,3))/3._dp
End Do 
Do j1 = 1,3
mat(10,10) = mat(10,10)+(sqrt(2._dp)*vR(j1)*kap(j1,3,3))/3._dp
End Do 

 
 Do i1=2,10
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
If (Maxval(Abs(Aimag(mat))).Eq.0._dp) Then 
Call EigenSystemQP(Real(mat,dp),Eig,UVa,ierr,test) 
 
   Do i1=1,10
   If ((Eig(i1).Lt.0._dp).or.(Abs(eig(i1)).lt.1E-15)) Then 
    MChi(i1) = - Eig(i1) 
    UV(i1,:) = (0._dp,1._dp)*UVa(i1,:) 
   Else 
    MChi(i1) = Eig(i1) 
    UV(i1,:) = UVa(i1,:)
    End If 
   End Do 
 
Do i1=1,9
  Do i2=i1+1,10
    If (MChi(i1).Gt.MChi(i2)) Then 
      Eig(1) = MChi(i1) 
      MChi(i1) = MChi(i2) 
      MChi(i2) =  Eig(1) 
      E10 = UV(i1,:) 
      UV(i1,:) = UV(i2,:) 
      UV(i2,:) = E10
    End If 
   End Do 
End Do 
 
Else 
 
mat2 = Matmul( Transpose(Conjg( mat) ), mat ) 
Call EigensystemQP(mat2, Eig, UV, ierr, test) 
mat2 = Matmul( Conjg(UV), Matmul( mat, Transpose( Conjg(UV)))) 
Do i1=1,10
  If (Eig(i1).ne.Eig(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  UV(i1,:)= phaseM * UV(i1,:) 
End if 
  If ((Abs(Eig(i1)).Le.MaxMassNumericalZero).and.(Eig(i1).lt.0._dp)) Eig(i1) = Abs(Eig(i1))+1.E-10_dp 
  If (Eig(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      !Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      !Write(*,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Call TerminateProgram
      return 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Eig(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,Eig(i1) 
  Eig(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
MChi = Sqrt( Eig ) 
 
End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  !Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMChi, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


If ((Abs(UV(1,2)).gt.Abs(UV(2,1))).And.(MChi(1).lt.1.0E-15_dp).And.(MChi(2).lt.1.0E-15_dp)) Then 

   MChitemp = MChi 

   UVtemp = UV 

   UV(1,:) = UVtemp(2,:) 

   UV(2,:) = UVtemp(1,:) 

   MChi(1) = MChitemp(2) 

   MChi(2) = MChitemp(1) 

End If 
 
Iname = Iname - 1 
 
End Subroutine CalculateMChi 

Subroutine CalculateMCha(g2,Ye,lam,Yv,M2,vd,vu,vL,vR,ZER,ZEL,MCha,kont)

Real(dp),Intent(in) :: g2,vd,vu,vL(3),vR(3)

Complex(dp),Intent(in) :: Ye(3,3),lam(3),Yv(3,3),M2

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MCha(5) 
 Complex(dp), Intent(out) :: ZER(5,5), ZEL(5,5) 
 
 Complex(dp) :: mat(5,5)=0._dp, mat2(5,5)=0._dp, phaseM 

Complex(dp) :: ZER2(5,5), ZEL2(5,5) 
 
 Real(dp) :: ZER1(5,5), ZEL1(5,5), test(2), MCha2(5) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMCha'
 
MCha = 0._dp 
ZER = 0._dp 
ZEL = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(vd*Ye(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vd*Ye(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(vd*Ye(3,1))/sqrt(2._dp)
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+(g2*vL(1))/sqrt(2._dp)
mat(1,5) = 0._dp 
Do j1 = 1,3
mat(1,5) = mat(1,5)-((vR(j1)*Yv(1,j1))/sqrt(2._dp))
End Do 
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vd*Ye(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(vd*Ye(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(vd*Ye(3,2))/sqrt(2._dp)
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(g2*vL(2))/sqrt(2._dp)
mat(2,5) = 0._dp 
Do j1 = 1,3
mat(2,5) = mat(2,5)-((vR(j1)*Yv(2,j1))/sqrt(2._dp))
End Do 
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(vd*Ye(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(vd*Ye(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(vd*Ye(3,3))/sqrt(2._dp)
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(g2*vL(3))/sqrt(2._dp)
mat(3,5) = 0._dp 
Do j1 = 1,3
mat(3,5) = mat(3,5)-((vR(j1)*Yv(3,j1))/sqrt(2._dp))
End Do 
mat(4,1) = 0._dp 
mat(4,2) = 0._dp 
mat(4,3) = 0._dp 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+M2
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+(g2*vu)/sqrt(2._dp)
mat(5,1) = 0._dp 
Do j1 = 1,3
mat(5,1) = mat(5,1)-((vL(j1)*Ye(1,j1))/sqrt(2._dp))
End Do 
mat(5,2) = 0._dp 
Do j1 = 1,3
mat(5,2) = mat(5,2)-((vL(j1)*Ye(2,j1))/sqrt(2._dp))
End Do 
mat(5,3) = 0._dp 
Do j1 = 1,3
mat(5,3) = mat(5,3)-((vL(j1)*Ye(3,j1))/sqrt(2._dp))
End Do 
mat(5,4) = 0._dp 
mat(5,4) = mat(5,4)+(g2*vd)/sqrt(2._dp)
mat(5,5) = 0._dp 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*lam(j1))/sqrt(2._dp)
End Do 

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2,ZEL1,ierr,test) 
ZEL2 = ZEL1 
Else 
Call EigenSystem(mat2,MCha2,ZEL2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MCha2,ZER1,ierr,test) 
                  
                  
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MCha2,ZER2,ierr,test) 
 
 
End If 
ZER2 = Conjg(ZER2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZER2),mat),Transpose( Conjg(ZEL2))) 
Do i1=1,5
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZEL2(i1,:) = phaseM *ZEL2(i1,:) 
 End if 
End Do 
 
Do i1=1,5
If (Abs(ZEL2(i1,i1)).gt.0._dp) Then 
phaseM = ZEL2(i1,i1) / Abs(ZEL2(i1,i1)) 
ZEL2(i1,:) = Conjg(phaseM) *ZEL2(i1,:) 
 ZER2(i1,:) = phaseM *ZER2(i1,:) 
 End if 
  If (MCha2(i1).ne.MCha2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Abs(MCha2(i1)).Le.MaxMassNumericalZero) MCha2(i1) = Abs(MCha2(i1))+1.E-10_dp 
  If (MCha2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MCha2(i1) 
      !Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      !Write(*,*) 'a mass squarred is negative: ',i1,MCha2(i1) 
      Call TerminateProgram
      return 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MCha2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MCha2(i1) 
  MCha2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMCha, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MCha = Sqrt( MCha2 ) 
ZER = ZER2 
ZEL = ZEL2 
Iname = Iname - 1 
 
End Subroutine CalculateMCha 

Subroutine CalculateMFd(Yd,vd,ZDL,ZDR,MFd,kont)

Real(dp),Intent(in) :: vd

Complex(dp),Intent(in) :: Yd(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFd(3) 
 Complex(dp), Intent(out) :: ZDL(3,3), ZDR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZDL2(3,3), ZDR2(3,3) 
 
 Real(dp) :: ZDL1(3,3), ZDR1(3,3), test(2), MFd2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFd'
 
MFd = 0._dp 
ZDL = 0._dp 
ZDR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(vd*Yd(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vd*Yd(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(vd*Yd(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vd*Yd(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(vd*Yd(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(vd*Yd(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(vd*Yd(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(vd*Yd(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(vd*Yd(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2,ZDR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFd2,ZDL1,ierr,test) 
                  
                  
ZDL2 = ZDL1 
Else 
Call EigenSystem(mat2,MFd2,ZDL2,ierr,test) 
 
 
End If 
ZDL2 = Conjg(ZDL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZDL2),mat),Transpose( Conjg(ZDR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZDR2(i1,:) = phaseM *ZDR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZDR2(i1,i1)).gt.0._dp) Then 
phaseM = ZDR2(i1,i1) / Abs(ZDR2(i1,i1)) 
ZDR2(i1,:) = Conjg(phaseM) *ZDR2(i1,:) 
 ZDL2(i1,:) = phaseM *ZDL2(i1,:) 
 End if 
  If (MFd2(i1).ne.MFd2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram
      return 
    End If 
  If (Abs(MFd2(i1)).Le.MaxMassNumericalZero) MFd2(i1) = Abs(MFd2(i1))+1.E-10_dp 
  If (MFd2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MFd2(i1) 
      !Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      !Write(*,*) 'a mass squarred is negative: ',i1,MFd2(i1) 
      Call TerminateProgram 
      return
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MFd2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MFd2(i1) 
  MFd2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFd, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFd = Sqrt( MFd2 ) 
ZDL = ZDL2 
ZDR = ZDR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFd 

Subroutine CalculateMFu(Yu,vu,ZUL,ZUR,MFu,kont)

Real(dp),Intent(in) :: vu

Complex(dp),Intent(in) :: Yu(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFu(3) 
 Complex(dp), Intent(out) :: ZUL(3,3), ZUR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZUL2(3,3), ZUR2(3,3) 
 
 Real(dp) :: ZUL1(3,3), ZUR1(3,3), test(2), MFu2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFu'
 
MFu = 0._dp 
ZUL = 0._dp 
ZUR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(vu*Yu(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vu*Yu(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(vu*Yu(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vu*Yu(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(vu*Yu(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(vu*Yu(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(vu*Yu(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(vu*Yu(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(vu*Yu(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2,ZUR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFu2,ZUL1,ierr,test) 
                  
                  
ZUL2 = ZUL1 
Else 
Call EigenSystem(mat2,MFu2,ZUL2,ierr,test) 
 
 
End If 
ZUL2 = Conjg(ZUL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZUL2),mat),Transpose( Conjg(ZUR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZUR2(i1,:) = phaseM *ZUR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZUR2(i1,i1)).gt.0._dp) Then 
phaseM = ZUR2(i1,i1) / Abs(ZUR2(i1,i1)) 
ZUR2(i1,:) = Conjg(phaseM) *ZUR2(i1,:) 
 ZUL2(i1,:) = phaseM *ZUL2(i1,:) 
 End if 
  If (MFu2(i1).ne.MFu2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Abs(MFu2(i1)).Le.MaxMassNumericalZero) MFu2(i1) = Abs(MFu2(i1))+1.E-10_dp 
  If (MFu2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MFu2(i1) 
      !Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      !Write(*,*) 'a mass squarred is negative: ',i1,MFu2(i1) 
      Call TerminateProgram 
      return
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MFu2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,MFu2(i1) 
  MFu2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFu, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFu = Sqrt( MFu2 ) 
ZUL = ZUL2 
ZUR = ZUR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFu 

Subroutine CalculateVPVZ(g1,g2,vd,vu,vL,ZZ,MVZ,MVZ2,TW,kont)

Real(dp), Intent(in) :: g1,g2,vd,vu,vL(3)

Real(dp), Intent(out) :: TW

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVZ, MVZ2
Real(dp) :: VPVZ2(2),VPVZ(2)  

Real(dp), Intent(out) :: ZZ(2,2) 
 
Real(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVPVZ'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g1**2*vd**2)/4._dp
mat(1,1) = mat(1,1)+(g1**2*vu**2)/4._dp
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/4._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g1*g2*vd**2)/4._dp
mat(1,2) = mat(1,2)-(g1*g2*vu**2)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)-(g1*g2*vL(j1)**2)/4._dp
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*vd**2)/4._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2)/4._dp
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/4._dp
End Do 

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,VPVZ2,ZZ,ierr,test) 
 
 
ZZ = Transpose(ZZ) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(VPVZ2(i1)).Le.1.E-10_dp*(Maxval(VPVZ2))) VPVZ2(i1) = 1.E-10_dp 
  If (VPVZ2(i1).ne.VPVZ2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (VPVZ2(i1).Ge.0._dp) Then 
  VPVZ(i1) =Sqrt(VPVZ2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,VPVZ2(i1) 
    End If 
  VPVZ(i1)= 1._dp 
  VPVZ2(i1)= 1._dp  
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,VPVZ2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,VPVZ2(i1) 
  VPVZ(i1)= 1._dp 
  VPVZ2(i1) = 1._dp  
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
MVZ = VPVZ(2) 
MVZ2 = VPVZ2(2) 
TW = ACos(Abs(ZZ(1,1)))

 Iname = Iname - 1 
 
End Subroutine CalculateVPVZ 

Subroutine CalculateVWm(g2,vd,vu,vL,ZW,MVWm,MVWm2,kont)

Real(dp), Intent(in) :: g2,vd,vu,vL(3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVWm, MVWm2
Real(dp) :: VWm2(2),VWm(2)  

Complex(dp), Intent(out) :: ZW(2,2) 
 
Complex(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVWm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g2**2*vd**2)/4._dp
mat(1,1) = mat(1,1)+(g2**2*vu**2)/4._dp
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/4._dp
End Do 
mat(1,2) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*vd**2)/4._dp
mat(2,2) = mat(2,2)+(g2**2*vu**2)/4._dp
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/4._dp
End Do 

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,VWm2,ZW,ierr,test) 
 
 
ZW = Transpose(ZW) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(VWm2(i1)).Le.1.E-10_dp*(Maxval(VWm2))) VWm2(i1) = 1.E-10_dp 
  If (VWm2(i1).ne.VWm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (VWm2(i1).Ge.0._dp) Then 
  VWm(i1) =Sqrt(VWm2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,VWm2(i1) 
    End If 
  VWm(i1)= 1._dp 
  VWm2(i1)= 1._dp  
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,VWm2(i1) 
     !Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     !Write(*,*) 'in the calculation of the masses' 
     !Write(*,*) 'occurred a negative mass squared!' 
     !Write(*,*) i1,VWm2(i1) 
  VWm(i1)= 1._dp 
  VWm2(i1) = 1._dp  
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
MVWm = VWm(1) 
MVWm2 = VWm2(1) 

 Iname = Iname - 1 
 
End Subroutine CalculateVWm 

Subroutine CalculateMSdEffPot(g1,g2,Yd,Td,lam,mq2,md2,vd,vu,vL,vR,ZD,MSd,             & 
& MSd2,kont)

Real(dp), Intent(in) :: g1,g2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: Yd(3,3),Td(3,3),lam(3),mq2(3,3),md2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MSd(6), MSd2(6) 
Complex(dp), Intent(out) :: ZD(6,6) 
 
Complex(dp) :: mat(6,6)  

Real(dp) :: MSd2temp(6), Q2 
Complex(dp) :: ZDtemp(6,6), ZDtemp2(6,6) 
 
Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMSd'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)-(g1**2*vdFix**2)/24._dp
mat(1,1) = mat(1,1)-(g2**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)+(g1**2*vuFix**2)/24._dp
mat(1,1) = mat(1,1)+(g2**2*vuFix**2)/8._dp
mat(1,1) = mat(1,1)+mq2(1,1)
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+mq2(1,2)
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,2))/2._dp
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+mq2(1,3)
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,3))/2._dp
End Do 
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+(vd*Conjg(Td(1,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(Yd(1,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+(vd*Conjg(Td(2,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(Yd(2,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+(vd*Conjg(Td(3,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu*Conjg(Yd(3,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)-(g1**2*vdFix**2)/24._dp
mat(2,2) = mat(2,2)-(g2**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vuFix**2)/24._dp
mat(2,2) = mat(2,2)+(g2**2*vuFix**2)/8._dp
mat(2,2) = mat(2,2)+mq2(2,2)
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vd**2*Conjg(Yd(j1,2))*Yd(j1,2))/2._dp
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+mq2(2,3)
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd**2*Conjg(Yd(j1,2))*Yd(j1,3))/2._dp
End Do 
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(vd*Conjg(Td(1,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vu*Conjg(Yd(1,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(vd*Conjg(Td(2,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vu*Conjg(Yd(2,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(2,6) = 0._dp 
mat(2,6) = mat(2,6)+(vd*Conjg(Td(3,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vu*Conjg(Yd(3,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)-(g1**2*vdFix**2)/24._dp
mat(3,3) = mat(3,3)-(g2**2*vdFix**2)/8._dp
mat(3,3) = mat(3,3)+(g1**2*vuFix**2)/24._dp
mat(3,3) = mat(3,3)+(g2**2*vuFix**2)/8._dp
mat(3,3) = mat(3,3)+mq2(3,3)
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd**2*Conjg(Yd(j1,3))*Yd(j1,3))/2._dp
End Do 
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(vd*Conjg(Td(1,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yd(1,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(vd*Conjg(Td(2,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yd(2,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vd*Conjg(Td(3,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yd(3,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)-(g1**2*vdFix**2)/12._dp
mat(4,4) = mat(4,4)+(g1**2*vuFix**2)/12._dp
mat(4,4) = mat(4,4)+md2(1,1)
Do j1 = 1,3
mat(4,4) = mat(4,4)-(g1**2*vL(j1)**2)/12._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd**2*Conjg(Yd(1,j1))*Yd(1,j1))/2._dp
End Do 
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+md2(1,2)
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd**2*Conjg(Yd(2,j1))*Yd(1,j1))/2._dp
End Do 
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+md2(1,3)
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vd**2*Conjg(Yd(3,j1))*Yd(1,j1))/2._dp
End Do 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)-(g1**2*vdFix**2)/12._dp
mat(5,5) = mat(5,5)+(g1**2*vuFix**2)/12._dp
mat(5,5) = mat(5,5)+md2(2,2)
Do j1 = 1,3
mat(5,5) = mat(5,5)-(g1**2*vL(j1)**2)/12._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd**2*Conjg(Yd(2,j1))*Yd(2,j1))/2._dp
End Do 
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+md2(2,3)
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vd**2*Conjg(Yd(3,j1))*Yd(2,j1))/2._dp
End Do 
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)-(g1**2*vdFix**2)/12._dp
mat(6,6) = mat(6,6)+(g1**2*vuFix**2)/12._dp
mat(6,6) = mat(6,6)+md2(3,3)
Do j1 = 1,3
mat(6,6) = mat(6,6)-(g1**2*vL(j1)**2)/12._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vd**2*Conjg(Yd(3,j1))*Yd(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,MSd2,ZD,ierr,test) 
 
 
! Fix order
  ZDtemp2=ZD
Do i1=1,6
  pos=Maxloc(Abs(ZDtemp2(i1,:)),1)
  ZDtemp(pos,:)=ZD(i1,:)
  MSd2temp(pos)=MSd2(i1)
  ZDtemp2(:,pos)=0._dp
End do
  MSd2 = MSd2temp
  ZD = ZDtemp
! Fix phases
Do i1=1,6
  pos=Maxloc(Abs(ZD(i1,:)),1)
  If (Real(ZD(i1,pos),dp).lt.0._dp) Then
    ZD(i1,:)=-ZD(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,6
  If (MSd2(i1).ne.MSd2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (MSd2(i1).Ge.0._dp) Then 
  MSd(i1)=Sqrt(MSd2(i1) ) 
  Else 
  MSd = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMSdEffPot 

Subroutine CalculateMSuEffPot(g1,g2,lam,Yv,Yu,Tu,mq2,mu2,vd,vu,vL,vR,ZU,              & 
& MSu,MSu2,kont)

Real(dp), Intent(in) :: g1,g2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: lam(3),Yv(3,3),Yu(3,3),Tu(3,3),mq2(3,3),mu2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MSu(6), MSu2(6) 
Complex(dp), Intent(out) :: ZU(6,6) 
 
Complex(dp) :: mat(6,6)  

Real(dp) :: MSu2temp(6), Q2 
Complex(dp) :: ZUtemp(6,6), ZUtemp2(6,6) 
 
Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMSu'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)-(g1**2*vdFix**2)/24._dp
mat(1,1) = mat(1,1)+(g2**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)+(g1**2*vuFix**2)/24._dp
mat(1,1) = mat(1,1)-(g2**2*vuFix**2)/8._dp
mat(1,1) = mat(1,1)+mq2(1,1)
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+mq2(1,2)
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,2))/2._dp
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+mq2(1,3)
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,3))/2._dp
End Do 
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+(vu*Conjg(Tu(1,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vd*Conjg(Yu(1,1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(Conjg(Yu(1,1))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+(vu*Conjg(Tu(2,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vd*Conjg(Yu(2,1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(Conjg(Yu(2,1))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+(vu*Conjg(Tu(3,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vd*Conjg(Yu(3,1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,6) = mat(1,6)+(Conjg(Yu(3,1))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)-(g1**2*vdFix**2)/24._dp
mat(2,2) = mat(2,2)+(g2**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vuFix**2)/24._dp
mat(2,2) = mat(2,2)-(g2**2*vuFix**2)/8._dp
mat(2,2) = mat(2,2)+mq2(2,2)
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vu**2*Conjg(Yu(j1,2))*Yu(j1,2))/2._dp
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+mq2(2,3)
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu**2*Conjg(Yu(j1,2))*Yu(j1,3))/2._dp
End Do 
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(vu*Conjg(Tu(1,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(Yu(1,2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yu(1,2))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(vu*Conjg(Tu(2,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(Yu(2,2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yu(2,2))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,6) = 0._dp 
mat(2,6) = mat(2,6)+(vu*Conjg(Tu(3,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*Conjg(Yu(3,2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yu(3,2))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)-(g1**2*vdFix**2)/24._dp
mat(3,3) = mat(3,3)+(g2**2*vdFix**2)/8._dp
mat(3,3) = mat(3,3)+(g1**2*vuFix**2)/24._dp
mat(3,3) = mat(3,3)-(g2**2*vuFix**2)/8._dp
mat(3,3) = mat(3,3)+mq2(3,3)
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g1**2*vL(j1)**2)/24._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu**2*Conjg(Yu(j1,3))*Yu(j1,3))/2._dp
End Do 
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(vu*Conjg(Tu(1,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yu(1,3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yu(1,3))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(vu*Conjg(Tu(2,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yu(2,3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yu(2,3))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vu*Conjg(Tu(3,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yu(3,3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yu(3,3))*vL(j1)*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+(g1**2*vdFix**2)/6._dp
mat(4,4) = mat(4,4)-(g1**2*vuFix**2)/6._dp
mat(4,4) = mat(4,4)+mu2(1,1)
Do j1 = 1,3
mat(4,4) = mat(4,4)+(g1**2*vL(j1)**2)/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu**2*Conjg(Yu(1,j1))*Yu(1,j1))/2._dp
End Do 
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+mu2(1,2)
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yu(2,j1))*Yu(1,j1))/2._dp
End Do 
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+mu2(1,3)
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu**2*Conjg(Yu(3,j1))*Yu(1,j1))/2._dp
End Do 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+(g1**2*vdFix**2)/6._dp
mat(5,5) = mat(5,5)-(g1**2*vuFix**2)/6._dp
mat(5,5) = mat(5,5)+mu2(2,2)
Do j1 = 1,3
mat(5,5) = mat(5,5)+(g1**2*vL(j1)**2)/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu**2*Conjg(Yu(2,j1))*Yu(2,j1))/2._dp
End Do 
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+mu2(2,3)
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu**2*Conjg(Yu(3,j1))*Yu(2,j1))/2._dp
End Do 
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)+(g1**2*vdFix**2)/6._dp
mat(6,6) = mat(6,6)-(g1**2*vuFix**2)/6._dp
mat(6,6) = mat(6,6)+mu2(3,3)
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g1**2*vL(j1)**2)/6._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vu**2*Conjg(Yu(3,j1))*Yu(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,MSu2,ZU,ierr,test) 
 
 
! Fix order
  ZUtemp2=ZU
Do i1=1,6
  pos=Maxloc(Abs(ZUtemp2(i1,:)),1)
  ZUtemp(pos,:)=ZU(i1,:)
  MSu2temp(pos)=MSu2(i1)
  ZUtemp2(:,pos)=0._dp
End do
  MSu2 = MSu2temp
  ZU = ZUtemp
! Fix phases
Do i1=1,6
  pos=Maxloc(Abs(ZU(i1,:)),1)
  If (Real(ZU(i1,pos),dp).lt.0._dp) Then
    ZU(i1,:)=-ZU(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,6
  If (MSu2(i1).ne.MSu2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram
      return 
    End If 
  If (MSu2(i1).Ge.0._dp) Then 
  MSu(i1)=Sqrt(MSu2(i1) ) 
  Else 
  MSu = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMSuEffPot 

Subroutine CalculateMhhEffPot(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,             & 
& mHu2,mv2,vd,vu,vL,vR,ZH,Mhh,Mhh2,kont)

Real(dp), Intent(in) :: g1,g2,mlHd2(3),mHd2,mHu2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: lam(3),Tlam(3),Yv(3,3),Tv(3,3),kap(3,3,3),Tk(3,3,3),ml2(3,3),mv2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: Mhh(8), Mhh2(8) 
Real(dp), Intent(out) :: ZH(8,8) 
 
Real(dp) :: mat(8,8)  

Real(dp) :: Mhh2temp(8), Q2 
Real(dp) :: ZHtemp(8,8),ZHtemp2(8,8) 
 
Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMhh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+mHd2
mat(1,1) = mat(1,1)+(3*g1**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)+(3*g2**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)-(g1**2*vuFix**2)/8._dp
mat(1,1) = mat(1,1)-(g2**2*vuFix**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(1,1) = mat(1,1)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vu**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g1**2*vdFix*vuFix)/4._dp
mat(1,2) = mat(1,2)-(g2**2*vdFix*vuFix)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(Tlam(j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)+vd*vu*Conjg(lam(j1))*lam(j1)
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vR(j1)*Tlam(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vu*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vu*Conjg(Yv(j2,j1))*vL(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)-(vu*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vd*Conjg(lam(j1))*vR(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vd*Conjg(lam(1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,3) = mat(1,3)-(vu*Tlam(1))/(2._dp*sqrt(2._dp))
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)-(vu*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vd*Conjg(lam(j1))*vR(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vd*Conjg(lam(2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,4) = mat(1,4)-(vu*Tlam(2))/(2._dp*sqrt(2._dp))
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)-(vu*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vd*Conjg(lam(j1))*vR(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vd*Conjg(lam(3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,5) = mat(1,5)-(vu*Tlam(3))/(2._dp*sqrt(2._dp))
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+mlHd2(1)
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(lam(j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(Yv(1,j1))*lam(j1))/4._dp
End Do 
mat(1,6) = mat(1,6)+(g1**2*vd*vL(1))/4._dp
mat(1,6) = mat(1,6)+(g2**2*vd*vL(1))/4._dp
mat(1,7) = 0._dp 
mat(1,7) = mat(1,7)+mlHd2(2)
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(lam(j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(Yv(2,j1))*lam(j1))/4._dp
End Do 
mat(1,7) = mat(1,7)+(g1**2*vd*vL(2))/4._dp
mat(1,7) = mat(1,7)+(g2**2*vd*vL(2))/4._dp
mat(1,8) = 0._dp 
mat(1,8) = mat(1,8)+mlHd2(3)
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(lam(j1))*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(Yv(3,j1))*lam(j1))/4._dp
End Do 
mat(1,8) = mat(1,8)+(g1**2*vd*vL(3))/4._dp
mat(1,8) = mat(1,8)+(g2**2*vd*vL(3))/4._dp
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mHu2
mat(2,2) = mat(2,2)-(g1**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)-(g2**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)+(3*g1**2*vuFix**2)/8._dp
mat(2,2) = mat(2,2)+(3*g2**2*vuFix**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(2,2) = mat(2,2)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vd**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(Yv(j2,j1))*vL(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j1,j3))*vR(j2)*vR(j3)*Yv(j1,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j3,j1))*vL(j2)*vL(j3)*Yv(j2,j1))/2._dp
End Do 
End Do 
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-(vd*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Tv(j1,1))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(lam(j1))*vR(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(lam(1))*vR(j1)*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vL(j1)*Tv(j1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(Yv(j1,j2))*vR(j2)*Yv(j1,1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vu*Conjg(Yv(j1,1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vd*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(1,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(1,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j1,1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j1,j3,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j3,1,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(kap(j3,j1,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
End Do 
mat(2,3) = mat(2,3)-(vd*Tlam(1))/(2._dp*sqrt(2._dp))
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)-(vd*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Tv(j1,2))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(lam(j1))*vR(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(lam(2))*vR(j1)*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vL(j1)*Tv(j1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(Yv(j1,j2))*vR(j2)*Yv(j1,2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vu*Conjg(Yv(j1,2))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vd*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(2,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(2,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j1,2,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j1,j3,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j3,2,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(kap(j3,j1,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
End Do 
mat(2,4) = mat(2,4)-(vd*Tlam(2))/(2._dp*sqrt(2._dp))
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)-(vd*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Tv(j1,3))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(lam(j1))*vR(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(lam(3))*vR(j1)*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vL(j1)*Tv(j1,3))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(Yv(j1,j2))*vR(j2)*Yv(j1,3))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vu*Conjg(Yv(j1,3))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vd*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(3,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(3,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j1,3,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j1,j3,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j3,3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(kap(j3,j1,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
End Do 
mat(2,5) = mat(2,5)-(vd*Tlam(3))/(2._dp*sqrt(2._dp))
mat(2,6) = 0._dp 
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Tv(1,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*vu*Conjg(lam(j1))*Yv(1,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*vu*Conjg(Yv(1,j1))*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)+(vR(j1)*Tv(1,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(vu*Conjg(Yv(1,j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)+(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,6) = mat(2,6)-(g1**2*vu*vL(1))/4._dp
mat(2,6) = mat(2,6)-(g2**2*vu*vL(1))/4._dp
mat(2,7) = 0._dp 
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Tv(2,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vd*vu*Conjg(lam(j1))*Yv(2,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vd*vu*Conjg(Yv(2,j1))*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vR(j1)*Tv(2,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vu*Conjg(Yv(2,j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)+(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,7) = mat(2,7)-(g1**2*vu*vL(2))/4._dp
mat(2,7) = mat(2,7)-(g2**2*vu*vL(2))/4._dp
mat(2,8) = 0._dp 
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Tv(3,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vd*vu*Conjg(lam(j1))*Yv(3,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vd*vu*Conjg(Yv(3,j1))*lam(j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)+(vR(j1)*Tv(3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(vu*Conjg(Yv(3,j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)+(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,8) = mat(2,8)-(g1**2*vu*vL(3))/4._dp
mat(2,8) = mat(2,8)-(g2**2*vu*vL(3))/4._dp
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+mv2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Tk(1,1,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Tk(1,j1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Tk(j1,1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(lam(j1))*kap(1,1,j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(lam(j1))*kap(1,j1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(lam(j1))*kap(j1,1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(kap(1,1,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(kap(1,j1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*vu*Conjg(kap(j1,1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vR(j1)*Tk(1,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vR(j1)*Tk(1,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vR(j1)*Tk(j1,1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(kap(1,1,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(kap(1,j1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(kap(j1,1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,1,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,1))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(3,3) = mat(3,3)+(vd**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,3) = mat(3,3)+(vu**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+mv2(1,2)/2._dp
mat(3,4) = mat(3,4)+mv2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(1,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(1,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(2,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(2,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(j1,1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Tk(j1,2,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(1,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(1,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(2,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(2,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(j1,1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(lam(j1))*kap(j1,2,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(1,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(1,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(2,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(2,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(j1,1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*vu*Conjg(kap(j1,2,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(1,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(1,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(2,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(2,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(j1,1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vR(j1)*Tk(j1,2,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(1,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(1,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(2,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(2,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(j1,1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(kap(j1,2,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+mv2(1,3)/2._dp
mat(3,5) = mat(3,5)+mv2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(1,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(1,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(3,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(3,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(j1,1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Tk(j1,3,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(1,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(1,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(3,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(3,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(j1,1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(lam(j1))*kap(j1,3,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(1,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(1,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(3,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(3,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(j1,1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*vu*Conjg(kap(j1,3,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(1,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(1,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(3,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(3,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(j1,1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vR(j1)*Tk(j1,3,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(1,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(1,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(3,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(3,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(j1,1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(kap(j1,3,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vu*Conjg(Tv(1,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(lam(1))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yv(1,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yv(1,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(1,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,6) = mat(3,6)+(vu*Tv(1,1))/(2._dp*sqrt(2._dp))
mat(3,7) = 0._dp 
mat(3,7) = mat(3,7)+(vu*Conjg(Tv(2,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(lam(1))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(Yv(2,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(Yv(2,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(2,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,7) = mat(3,7)+(vu*Tv(2,1))/(2._dp*sqrt(2._dp))
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vu*Conjg(Tv(3,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(lam(1))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(Yv(3,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(Yv(3,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(3,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,8) = mat(3,8)+(vu*Tv(3,1))/(2._dp*sqrt(2._dp))
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+mv2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,2))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Tk(2,2,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Tk(2,j1,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Tk(j1,2,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(lam(j1))*kap(2,2,j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(lam(j1))*kap(2,j1,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(lam(j1))*kap(j1,2,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(kap(2,2,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(kap(2,j1,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*vu*Conjg(kap(j1,2,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vR(j1)*Tk(2,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vR(j1)*Tk(2,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vR(j1)*Tk(j1,2,2))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(kap(2,2,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(kap(2,j1,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(kap(j1,2,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,2))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,2))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(4,4) = mat(4,4)+(vd**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,4) = mat(4,4)+(vu**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+mv2(2,3)/2._dp
mat(4,5) = mat(4,5)+mv2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(2,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(2,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(3,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(3,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(j1,2,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Tk(j1,3,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(2,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(2,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(3,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(3,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(j1,2,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(lam(j1))*kap(j1,3,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(2,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(2,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(3,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(3,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(j1,2,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*vu*Conjg(kap(j1,3,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(2,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(2,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(3,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(3,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(j1,2,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vR(j1)*Tk(j1,3,2))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(2,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(2,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(3,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(3,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(j1,2,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(kap(j1,3,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,2))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+(vu*Conjg(Tv(1,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(lam(2))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(Yv(1,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(Yv(1,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(1,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,6) = mat(4,6)+(vu*Tv(1,2))/(2._dp*sqrt(2._dp))
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(vu*Conjg(Tv(2,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(lam(2))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(Yv(2,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(Yv(2,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(2,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,7) = mat(4,7)+(vu*Tv(2,2))/(2._dp*sqrt(2._dp))
mat(4,8) = 0._dp 
mat(4,8) = mat(4,8)+(vu*Conjg(Tv(3,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(lam(2))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(Yv(3,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(Yv(3,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(3,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,8) = mat(4,8)+(vu*Tv(3,2))/(2._dp*sqrt(2._dp))
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+mv2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,3))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Tk(3,3,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Tk(3,j1,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Tk(j1,3,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(lam(j1))*kap(3,3,j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(lam(j1))*kap(3,j1,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(lam(j1))*kap(j1,3,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(kap(3,3,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(kap(3,j1,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*vu*Conjg(kap(j1,3,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*Tk(3,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*Tk(3,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*Tk(j1,3,3))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(kap(3,3,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(kap(3,j1,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(kap(j1,3,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,3,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,3))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,3))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(5,5) = mat(5,5)+(vd**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,5) = mat(5,5)+(vu**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(vu*Conjg(Tv(1,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(lam(3))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(Yv(1,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(Yv(1,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(1,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,6) = mat(5,6)+(vu*Tv(1,3))/(2._dp*sqrt(2._dp))
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)+(vu*Conjg(Tv(2,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(lam(3))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(Yv(2,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(Yv(2,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(2,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,7) = mat(5,7)+(vu*Tv(2,3))/(2._dp*sqrt(2._dp))
mat(5,8) = 0._dp 
mat(5,8) = mat(5,8)+(vu*Conjg(Tv(3,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(lam(3))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(Yv(3,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(Yv(3,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(3,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,8) = mat(5,8)+(vu*Tv(3,3))/(2._dp*sqrt(2._dp))
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)+(g1**2*vdFix**2)/8._dp
mat(6,6) = mat(6,6)+(g2**2*vdFix**2)/8._dp
mat(6,6) = mat(6,6)-(g1**2*vuFix**2)/8._dp
mat(6,6) = mat(6,6)-(g2**2*vuFix**2)/8._dp
mat(6,6) = mat(6,6)+ml2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(6,6) = mat(6,6)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vu**2*Conjg(Yv(1,j1))*Yv(1,j1))/2._dp
End Do 
mat(6,6) = mat(6,6)+(g1**2*vL(1)**2)/4._dp
mat(6,6) = mat(6,6)+(g2**2*vL(1)**2)/4._dp
mat(6,7) = 0._dp 
mat(6,7) = mat(6,7)+ml2(1,2)/2._dp
mat(6,7) = mat(6,7)+ml2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(2,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(1,j1))*Yv(2,j1))/4._dp
End Do 
mat(6,7) = mat(6,7)+(g1**2*vL(1)*vL(2))/4._dp
mat(6,7) = mat(6,7)+(g2**2*vL(1)*vL(2))/4._dp
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)+ml2(1,3)/2._dp
mat(6,8) = mat(6,8)+ml2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(3,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(1,j1))*Yv(3,j1))/4._dp
End Do 
mat(6,8) = mat(6,8)+(g1**2*vL(1)*vL(3))/4._dp
mat(6,8) = mat(6,8)+(g2**2*vL(1)*vL(3))/4._dp
mat(7,7) = 0._dp 
mat(7,7) = mat(7,7)+(g1**2*vdFix**2)/8._dp
mat(7,7) = mat(7,7)+(g2**2*vdFix**2)/8._dp
mat(7,7) = mat(7,7)-(g1**2*vuFix**2)/8._dp
mat(7,7) = mat(7,7)-(g2**2*vuFix**2)/8._dp
mat(7,7) = mat(7,7)+ml2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(7,7) = mat(7,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(vu**2*Conjg(Yv(2,j1))*Yv(2,j1))/2._dp
End Do 
mat(7,7) = mat(7,7)+(g1**2*vL(2)**2)/4._dp
mat(7,7) = mat(7,7)+(g2**2*vL(2)**2)/4._dp
mat(7,8) = 0._dp 
mat(7,8) = mat(7,8)+ml2(2,3)/2._dp
mat(7,8) = mat(7,8)+ml2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(3,j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(2,j1))*Yv(3,j1))/4._dp
End Do 
mat(7,8) = mat(7,8)+(g1**2*vL(2)*vL(3))/4._dp
mat(7,8) = mat(7,8)+(g2**2*vL(2)*vL(3))/4._dp
mat(8,8) = 0._dp 
mat(8,8) = mat(8,8)+(g1**2*vdFix**2)/8._dp
mat(8,8) = mat(8,8)+(g2**2*vdFix**2)/8._dp
mat(8,8) = mat(8,8)-(g1**2*vuFix**2)/8._dp
mat(8,8) = mat(8,8)-(g2**2*vuFix**2)/8._dp
mat(8,8) = mat(8,8)+ml2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(8,8) = mat(8,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(vu**2*Conjg(Yv(3,j1))*Yv(3,j1))/2._dp
End Do 
mat(8,8) = mat(8,8)+(g1**2*vL(3)**2)/4._dp
mat(8,8) = mat(8,8)+(g2**2*vL(3)**2)/4._dp

 
 Do i1=2,8
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,Mhh2,ZH,ierr,test) 
 
 
! Fix phases
Do i1=1,8
  pos=Maxloc(Abs(ZH(i1,:)),1)
  If (Real(ZH(i1,pos),dp).lt.0._dp) Then
    ZH(i1,:)=-ZH(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,8
  If (Mhh2(i1).ne.Mhh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Mhh2(i1).Ge.0._dp) Then 
  Mhh(i1)=Sqrt(Mhh2(i1) ) 
  Else 
  Mhh = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMhhEffPot 

Subroutine CalculateMAhEffPot(g1,g2,lam,Tlam,Yv,Tv,kap,Tk,ml2,mlHd2,mHd2,             & 
& mHu2,mv2,vd,vu,vL,vR,TW,ZA,MAh,MAh2,kont)

Real(dp), Intent(in) :: g1,g2,mlHd2(3),mHd2,mHu2,vd,vu,vL(3),vR(3),TW

Complex(dp), Intent(in) :: lam(3),Tlam(3),Yv(3,3),Tv(3,3),kap(3,3,3),Tk(3,3,3),ml2(3,3),mv2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MAh(8), MAh2(8) 
Real(dp), Intent(out) :: ZA(8,8) 
 
Real(dp) :: ZAFIX(8,8) 
 
Real(dp) :: mat(8,8)  

Real(dp) ::  test(2), Q2 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMAh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+mHd2
mat(1,1) = mat(1,1)+(g1**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)-(g1**2*vuFix**2)/8._dp
mat(1,1) = mat(1,1)-(g2**2*vuFix**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(1,1) = mat(1,1)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(vu**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
mat(1,2) = 0._dp 
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(Tlam(j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vR(j1)*Tlam(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*lam(j1))/12._dp
End Do 
End Do 
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)-(vu*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)+(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)+(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)+(vu*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,3) = mat(1,3)-(vu*Tlam(1))/(2._dp*sqrt(2._dp))
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)-(vu*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)+(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)+(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)+(vu*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,4) = mat(1,4)-(vu*Tlam(2))/(2._dp*sqrt(2._dp))
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)-(vu*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)+(Conjg(lam(j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)+(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*lam(j1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*lam(3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)+(vu*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
mat(1,5) = mat(1,5)-(vu*Tlam(3))/(2._dp*sqrt(2._dp))
mat(1,6) = 0._dp 
mat(1,6) = mat(1,6)+mlHd2(1)
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,6) = mat(1,6)-(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(lam(j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu**2*Conjg(Yv(1,j1))*lam(j1))/4._dp
End Do 
mat(1,7) = 0._dp 
mat(1,7) = mat(1,7)+mlHd2(2)
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,7) = mat(1,7)-(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(lam(j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu**2*Conjg(Yv(2,j1))*lam(j1))/4._dp
End Do 
mat(1,8) = 0._dp 
mat(1,8) = mat(1,8)+mlHd2(3)
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(1,8) = mat(1,8)-(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*lam(j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(lam(j1))*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu**2*Conjg(Yv(3,j1))*lam(j1))/4._dp
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mHu2
mat(2,2) = mat(2,2)-(g1**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)-(g2**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vuFix**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vuFix**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(2,2) = mat(2,2)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(vd**2*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)-(vd*Conjg(Yv(j2,j1))*vL(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j1,j3))*vR(j2)*vR(j3)*Yv(j1,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j3,j1))*vL(j2)*vL(j3)*Yv(j2,j1))/2._dp
End Do 
End Do 
End Do 
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-(vd*Conjg(Tlam(1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,3) = mat(2,3)+(Conjg(Tv(j1,1))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vL(j1)*Tv(j1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(1,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(1,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j1,1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j1,j2,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j2,1,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*Conjg(kap(j2,j1,1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(1,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(1,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j1,1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j1,j3,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j3,1,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j3,j1,1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
End Do 
mat(2,3) = mat(2,3)-(vd*Tlam(1))/(2._dp*sqrt(2._dp))
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)-(vd*Conjg(Tlam(2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,4) = mat(2,4)+(Conjg(Tv(j1,2))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vL(j1)*Tv(j1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(2,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(2,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j1,2,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j1,j2,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j2,2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*Conjg(kap(j2,j1,2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(2,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(2,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j1,2,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j1,j3,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j3,2,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j3,j1,2))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
End Do 
mat(2,4) = mat(2,4)-(vd*Tlam(2))/(2._dp*sqrt(2._dp))
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)-(vd*Conjg(Tlam(3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
mat(2,5) = mat(2,5)+(Conjg(Tv(j1,3))*vL(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vL(j1)*Tv(j1,3))/(2._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(lam(j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(3,j1,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(3,j2,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j1,3,j2))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j1,j2,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j2,3,j1))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*Conjg(kap(j2,j1,3))*vR(j2)*lam(j1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(3,j1,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(3,j3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j1,3,j3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j1,j3,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j3,3,j1))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j3,j1,3))*vL(j2)*vR(j3)*Yv(j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(Yv(j3,j1))*vL(j3)*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
End Do 
mat(2,5) = mat(2,5)-(vd*Tlam(3))/(2._dp*sqrt(2._dp))
mat(2,6) = 0._dp 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Tv(1,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vR(j1)*Tv(1,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(1,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Yv(1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,7) = 0._dp 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Tv(2,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vR(j1)*Tv(2,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(2,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Yv(2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(2,8) = 0._dp 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Tv(3,j1))*vR(j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vR(j1)*Tv(3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(3,j1))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/12._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Yv(3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/12._dp
End Do 
End Do 
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+mv2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(Tk(1,1,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(Tk(1,j1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(Tk(j1,1,1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(lam(j1))*kap(1,1,j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(lam(j1))*kap(1,j1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(lam(j1))*kap(j1,1,1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(1))/2._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(kap(1,1,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(kap(1,j1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd*vu*Conjg(kap(j1,1,1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vR(j1)*Tk(1,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vR(j1)*Tk(1,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vR(j1)*Tk(j1,1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(kap(1,1,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(kap(1,j1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(kap(j1,1,1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,1,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,1))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(1,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,1,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(1,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,3) = mat(3,3)-(Conjg(kap(j1,1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(3,3) = mat(3,3)+(vd**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,3) = mat(3,3)+(vu**2*Conjg(lam(1))*lam(1))/2._dp
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+mv2(1,2)/2._dp
mat(3,4) = mat(3,4)+mv2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(1,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(1,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(2,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(2,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(j1,1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(Tk(j1,2,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(1,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(1,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(2,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(2,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(j1,1,2))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(lam(j1))*kap(j1,2,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(1,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(1,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(2,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(2,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(j1,1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd*vu*Conjg(kap(j1,2,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(1,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(1,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(2,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(2,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(j1,1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vR(j1)*Tk(j1,2,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(1,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(1,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(2,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(2,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(j1,1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(kap(j1,2,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(1,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(2,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,4) = mat(3,4)-(Conjg(kap(j1,2,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(2))*lam(1))/4._dp
mat(3,4) = mat(3,4)+(vd**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,4) = mat(3,4)+(vu**2*Conjg(lam(1))*lam(2))/4._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+mv2(1,3)/2._dp
mat(3,5) = mat(3,5)+mv2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(j2,1))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(1,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(1,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(3,1,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(3,j1,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(j1,1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(Tk(j1,3,1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vu**2*Conjg(Yv(j1,1))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(lam(1))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(1,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(1,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(3,1,j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(3,j1,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(j1,1,3))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(lam(j1))*kap(j1,3,1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vd*Conjg(Yv(j1,1))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(1,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(1,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(3,1,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(3,j1,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(j1,1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd*vu*Conjg(kap(j1,3,1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(1,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(1,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(3,1,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(3,j1,1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(j1,1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vR(j1)*Tk(j1,3,1))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(1,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(1,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(3,1,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(3,j1,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(j1,1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(kap(j1,3,1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(1,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,1))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(1,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,1,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(1,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j1,j3,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)+(Conjg(kap(j3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(1,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,1,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(3,j1,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(3,5) = mat(3,5)-(Conjg(kap(j1,3,1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(3))*lam(1))/4._dp
mat(3,5) = mat(3,5)+(vd**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,5) = mat(3,5)+(vu**2*Conjg(lam(1))*lam(3))/4._dp
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vu*Conjg(Tv(1,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,6) = mat(3,6)-(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vd*Conjg(lam(1))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)+(vd*Conjg(Yv(1,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vd*Conjg(Yv(1,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)+(Conjg(Yv(1,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,6) = mat(3,6)+(vu*Tv(1,1))/(2._dp*sqrt(2._dp))
mat(3,7) = 0._dp 
mat(3,7) = mat(3,7)+(vu*Conjg(Tv(2,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,7) = mat(3,7)-(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vd*Conjg(lam(1))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vd*Conjg(Yv(2,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vd*Conjg(Yv(2,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)+(Conjg(Yv(2,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,7) = mat(3,7)+(vu*Tv(2,1))/(2._dp*sqrt(2._dp))
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vu*Conjg(Tv(3,1)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)-(Conjg(Yv(j2,1))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(3,8) = mat(3,8)-(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vd*Conjg(lam(1))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)+(vd*Conjg(Yv(3,j1))*vR(j1)*lam(1))/4._dp
End Do 
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vd*Conjg(Yv(3,1))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,1))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(1,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(1,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j1,1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j1,j2,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j2,1,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(kap(j2,j1,1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)+(Conjg(Yv(3,1))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(1,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,1))/12._dp
End Do 
End Do 
mat(3,8) = mat(3,8)+(vu*Tv(3,1))/(2._dp*sqrt(2._dp))
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+mv2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,2))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(Tk(2,2,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(Tk(2,j1,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(Tk(j1,2,2))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(lam(j1))*kap(2,2,j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(lam(j1))*kap(2,j1,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(lam(j1))*kap(j1,2,2))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(2))/2._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(kap(2,2,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(kap(2,j1,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd*vu*Conjg(kap(j1,2,2))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vR(j1)*Tk(2,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vR(j1)*Tk(2,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vR(j1)*Tk(j1,2,2))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(kap(2,2,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(kap(2,j1,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(kap(j1,2,2))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,2))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,2))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(2,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,2,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(2,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,4) = mat(4,4)-(Conjg(kap(j1,2,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(4,4) = mat(4,4)+(vd**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,4) = mat(4,4)+(vu**2*Conjg(lam(2))*lam(2))/2._dp
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+mv2(2,3)/2._dp
mat(4,5) = mat(4,5)+mv2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(j2,2))*vL(j1)*vL(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(2,3,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(2,j1,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(3,2,j1))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(3,j1,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(j1,2,3))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(Tk(j1,3,2))*vR(j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vu**2*Conjg(Yv(j1,2))*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(lam(2))*vL(j1)*Yv(j1,3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(2,3,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(2,j1,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(3,2,j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(3,j1,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(j1,2,3))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(lam(j1))*kap(j1,3,2))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vd*Conjg(Yv(j1,2))*vL(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(2,3,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(2,j1,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(3,2,j1))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(3,j1,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(j1,2,3))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd*vu*Conjg(kap(j1,3,2))*lam(j1))/12._dp
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(2,3,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(2,j1,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(3,2,j1))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(3,j1,2))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(j1,2,3))/(6._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vR(j1)*Tk(j1,3,2))/(6._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(2,3,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(2,j1,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(3,2,j1))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(3,j1,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(j1,2,3))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(kap(j1,3,2))*vL(j2)*Yv(j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(2,j1,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,2))/12._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(2,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j1,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(3,j2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,2,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j1,j2,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j1,j3,j2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,2,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,2))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(2,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,2,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j1,j3,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)+(Conjg(kap(j3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j1,j3))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(2,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,2,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(3,j1,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,2,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(4,5) = mat(4,5)-(Conjg(kap(j1,3,2))*vR(j2)*vR(j3)*kap(j2,j3,j1))/36._dp
End Do 
End Do 
End Do 
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(3))*lam(2))/4._dp
mat(4,5) = mat(4,5)+(vd**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,5) = mat(4,5)+(vu**2*Conjg(lam(2))*lam(3))/4._dp
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+(vu*Conjg(Tv(1,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,6) = mat(4,6)-(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vd*Conjg(lam(2))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)+(vd*Conjg(Yv(1,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vd*Conjg(Yv(1,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)+(Conjg(Yv(1,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,6) = mat(4,6)+(vu*Tv(1,2))/(2._dp*sqrt(2._dp))
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(vu*Conjg(Tv(2,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,7) = mat(4,7)-(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vd*Conjg(lam(2))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)+(vd*Conjg(Yv(2,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vd*Conjg(Yv(2,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)+(Conjg(Yv(2,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,7) = mat(4,7)+(vu*Tv(2,2))/(2._dp*sqrt(2._dp))
mat(4,8) = 0._dp 
mat(4,8) = mat(4,8)+(vu*Conjg(Tv(3,2)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)-(Conjg(Yv(j2,2))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(4,8) = mat(4,8)-(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,2))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vd*Conjg(lam(2))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)+(vd*Conjg(Yv(3,j1))*vR(j1)*lam(2))/4._dp
End Do 
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vd*Conjg(Yv(3,2))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(2,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(2,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j1,2,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j1,j2,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j2,2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(kap(j2,j1,2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)+(Conjg(Yv(3,2))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(2,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,2,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,2))/12._dp
End Do 
End Do 
mat(4,8) = mat(4,8)+(vu*Tv(3,2))/(2._dp*sqrt(2._dp))
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+mv2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Yv(j2,3))*vL(j1)*vL(j2)*Yv(j1,3))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(Tk(3,3,j1))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(Tk(3,j1,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(Tk(j1,3,3))*vR(j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vu**2*Conjg(Yv(j1,3))*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(lam(3))*vL(j1)*Yv(j1,3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(lam(j1))*kap(3,3,j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(lam(j1))*kap(3,j1,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(lam(j1))*kap(j1,3,3))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vd*Conjg(Yv(j1,3))*vL(j1)*lam(3))/2._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(kap(3,3,j1))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(kap(3,j1,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd*vu*Conjg(kap(j1,3,3))*lam(j1))/6._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vR(j1)*Tk(3,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vR(j1)*Tk(3,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vR(j1)*Tk(j1,3,3))/(3._dp*sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(kap(3,3,j1))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(kap(3,j1,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(kap(j1,3,3))*vL(j2)*Yv(j2,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,3,j1))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(3,j1,3))/6._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*kap(j1,3,3))/6._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j1,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(3,j2,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j1,j2,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j1,j3,j2))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j1,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(3,j3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,3,j3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j1,j3,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)+(Conjg(kap(j3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j1,j3))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,3,j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(3,j1,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(5,5) = mat(5,5)-(Conjg(kap(j1,3,3))*vR(j2)*vR(j3)*kap(j2,j3,j1))/18._dp
End Do 
End Do 
End Do 
mat(5,5) = mat(5,5)+(vd**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,5) = mat(5,5)+(vu**2*Conjg(lam(3))*lam(3))/2._dp
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(vu*Conjg(Tv(1,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,6) = mat(5,6)-(Conjg(Yv(1,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(lam(j1))*vR(j1)*Yv(1,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vd*Conjg(lam(3))*vR(j1)*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)+(vd*Conjg(Yv(1,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vd*Conjg(Yv(1,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(1,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(1,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)+(Conjg(Yv(1,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Yv(1,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,6) = mat(5,6)+(vu*Tv(1,3))/(2._dp*sqrt(2._dp))
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)+(vu*Conjg(Tv(2,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,7) = mat(5,7)-(Conjg(Yv(2,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(lam(j1))*vR(j1)*Yv(2,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vd*Conjg(lam(3))*vR(j1)*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)+(vd*Conjg(Yv(2,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vd*Conjg(Yv(2,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(2,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)+(Conjg(Yv(2,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Yv(2,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,7) = mat(5,7)+(vu*Tv(2,3))/(2._dp*sqrt(2._dp))
mat(5,8) = 0._dp 
mat(5,8) = mat(5,8)+(vu*Conjg(Tv(3,3)))/(2._dp*sqrt(2._dp))
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)-(Conjg(Yv(j2,3))*vL(j2)*vR(j1)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(5,8) = mat(5,8)-(Conjg(Yv(3,j2))*vL(j1)*vR(j2)*Yv(j1,3))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(lam(j1))*vR(j1)*Yv(3,3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vd*Conjg(lam(3))*vR(j1)*Yv(3,j1))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)+(vd*Conjg(Yv(3,j1))*vR(j1)*lam(3))/4._dp
End Do 
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vd*Conjg(Yv(3,3))*vR(j1)*lam(j1))/4._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(j1,j2))*vL(j1)*vR(j2)*Yv(3,3))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(3,j1,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(3,j2,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j1,3,j2))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j1,j2,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j2,3,j1))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(kap(j2,j1,3))*vR(j2)*Yv(3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)+(Conjg(Yv(3,3))*vL(j1)*vR(j2)*Yv(j1,j2))/4._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j1,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(3,j2,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,3,j2))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j1,j2,3))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,3,j1))/12._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Yv(3,j1))*vR(j2)*kap(j2,j1,3))/12._dp
End Do 
End Do 
mat(5,8) = mat(5,8)+(vu*Tv(3,3))/(2._dp*sqrt(2._dp))
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)+(g1**2*vdFix**2)/8._dp
mat(6,6) = mat(6,6)+(g2**2*vdFix**2)/8._dp
mat(6,6) = mat(6,6)-(g1**2*vuFix**2)/8._dp
mat(6,6) = mat(6,6)-(g2**2*vuFix**2)/8._dp
mat(6,6) = mat(6,6)+ml2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(6,6) = mat(6,6)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vu**2*Conjg(Yv(1,j1))*Yv(1,j1))/2._dp
End Do 
mat(6,7) = 0._dp 
mat(6,7) = mat(6,7)+ml2(1,2)/2._dp
mat(6,7) = mat(6,7)+ml2(2,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(2,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vu**2*Conjg(Yv(1,j1))*Yv(2,j1))/4._dp
End Do 
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)+ml2(1,3)/2._dp
mat(6,8) = mat(6,8)+ml2(3,1)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(1,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(3,j1))*Yv(1,j1))/4._dp
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vu**2*Conjg(Yv(1,j1))*Yv(3,j1))/4._dp
End Do 
mat(7,7) = 0._dp 
mat(7,7) = mat(7,7)+(g1**2*vdFix**2)/8._dp
mat(7,7) = mat(7,7)+(g2**2*vdFix**2)/8._dp
mat(7,7) = mat(7,7)-(g1**2*vuFix**2)/8._dp
mat(7,7) = mat(7,7)-(g2**2*vuFix**2)/8._dp
mat(7,7) = mat(7,7)+ml2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(7,7) = mat(7,7)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(vu**2*Conjg(Yv(2,j1))*Yv(2,j1))/2._dp
End Do 
mat(7,8) = 0._dp 
mat(7,8) = mat(7,8)+ml2(2,3)/2._dp
mat(7,8) = mat(7,8)+ml2(3,2)/2._dp
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(2,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(3,j1))/4._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(3,j1))*Yv(2,j1))/4._dp
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vu**2*Conjg(Yv(2,j1))*Yv(3,j1))/4._dp
End Do 
mat(8,8) = 0._dp 
mat(8,8) = mat(8,8)+(g1**2*vdFix**2)/8._dp
mat(8,8) = mat(8,8)+(g2**2*vdFix**2)/8._dp
mat(8,8) = mat(8,8)-(g1**2*vuFix**2)/8._dp
mat(8,8) = mat(8,8)-(g2**2*vuFix**2)/8._dp
mat(8,8) = mat(8,8)+ml2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(8,8) = mat(8,8)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(vu**2*Conjg(Yv(3,j1))*Yv(3,j1))/2._dp
End Do 

 
 Do i1=2,8
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MAh2,ZA,ierr,test) 
 
 
! Fix phases
Do i1=1,8
  pos=Maxloc(Abs(ZA(i1,:)),1)
  If (Real(ZA(i1,pos),dp).lt.0._dp) Then
    ZA(i1,:)=-ZA(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,8
  If (MAh2(i1).ne.MAh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram
      return 
    End If 
  If (MAh2(i1).Ge.0._dp) Then 
  MAh(i1)=Sqrt(MAh2(i1) ) 
  Else 
  MAh = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMAhEffPot 

Subroutine CalculateMHpmEffPot(g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,ml2,mlHd2,              & 
& mHd2,mHu2,me2,vd,vu,vL,vR,ZP,MHpm,MHpm2,kont)

Real(dp), Intent(in) :: g1,g2,mlHd2(3),mHd2,mHu2,vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: Ye(3,3),Te(3,3),lam(3),Tlam(3),Yv(3,3),Tv(3,3),kap(3,3,3),ml2(3,3),me2(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MHpm(8), MHpm2(8) 
Real(dp), Intent(out) :: ZP(8,8) 
 
Real(dp) :: ZPFIX(8,8) 
 
Real(dp) :: mat(8,8)  

Real(dp) ::  test(2), Q2 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMHpm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+mHd2
mat(1,1) = mat(1,1)+(g1**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vdFix**2)/8._dp
mat(1,1) = mat(1,1)-(g1**2*vuFix**2)/8._dp
mat(1,1) = mat(1,1)+(g2**2*vuFix**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(1,1) = mat(1,1)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(1,1) = mat(1,1)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,1) = mat(1,1)+(Conjg(Ye(j1,j3))*vL(j2)*vL(j3)*Ye(j1,j2))/2._dp
End Do 
End Do 
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(g2**2*vdFix*vuFix)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(Tlam(j1))*vR(j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
mat(1,2) = mat(1,2)-(vd*vu*Conjg(lam(j1))*lam(j1))/2._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(vu*Conjg(lam(j1))*vL(j2)*Yv(j2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j1,j3,j2))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j1,j3))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(1,2) = mat(1,2)+(Conjg(lam(j1))*vR(j2)*vR(j3)*kap(j2,j3,j1))/6._dp
End Do 
End Do 
End Do 
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+mlHd2(1)
Do j1 = 1,3
Do j2 = 1,3
mat(1,3) = mat(1,3)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,3) = mat(1,3)-(vd*Conjg(Ye(j1,j2))*vL(j2)*Ye(j1,1))/2._dp
End Do 
End Do 
mat(1,3) = mat(1,3)+(g2**2*vd*vL(1))/4._dp
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+mlHd2(2)
Do j1 = 1,3
Do j2 = 1,3
mat(1,4) = mat(1,4)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,4) = mat(1,4)-(vd*Conjg(Ye(j1,j2))*vL(j2)*Ye(j1,2))/2._dp
End Do 
End Do 
mat(1,4) = mat(1,4)+(g2**2*vd*vL(2))/4._dp
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+mlHd2(3)
Do j1 = 1,3
Do j2 = 1,3
mat(1,5) = mat(1,5)-(Conjg(lam(j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,5) = mat(1,5)-(vd*Conjg(Ye(j1,j2))*vL(j2)*Ye(j1,3))/2._dp
End Do 
End Do 
mat(1,5) = mat(1,5)+(g2**2*vd*vL(3))/4._dp
mat(1,6) = 0._dp 
Do j1 = 1,3
mat(1,6) = mat(1,6)-((Conjg(Te(1,j1))*vL(j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,6) = mat(1,6)-(vu*Conjg(Ye(1,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,7) = 0._dp 
Do j1 = 1,3
mat(1,7) = mat(1,7)-((Conjg(Te(2,j1))*vL(j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,7) = mat(1,7)-(vu*Conjg(Ye(2,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(1,8) = 0._dp 
Do j1 = 1,3
mat(1,8) = mat(1,8)-((Conjg(Te(3,j1))*vL(j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(1,8) = mat(1,8)-(vu*Conjg(Ye(3,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mHu2
mat(2,2) = mat(2,2)-(g1**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vdFix**2)/8._dp
mat(2,2) = mat(2,2)+(g1**2*vuFix**2)/8._dp
mat(2,2) = mat(2,2)+(g2**2*vuFix**2)/8._dp
Do j1 = 1,3
Do j2 = 1,3
mat(2,2) = mat(2,2)+(Conjg(lam(j2))*vR(j1)*vR(j2)*lam(j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)-(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/8._dp
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,2) = mat(2,2)+(Conjg(Yv(j1,j3))*vR(j2)*vR(j3)*Yv(j1,j2))/2._dp
End Do 
End Do 
End Do 
mat(2,3) = 0._dp 
Do j1 = 1,3
mat(2,3) = mat(2,3)+(vd*vu*Conjg(lam(j1))*Yv(1,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,3) = mat(2,3)-((vR(j1)*Tv(1,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(1,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(1,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,3) = mat(2,3)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(1,j1))/6._dp
End Do 
End Do 
End Do 
mat(2,3) = mat(2,3)+(g2**2*vu*vL(1))/4._dp
mat(2,4) = 0._dp 
Do j1 = 1,3
mat(2,4) = mat(2,4)+(vd*vu*Conjg(lam(j1))*Yv(2,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,4) = mat(2,4)-((vR(j1)*Tv(2,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(2,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(2,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,4) = mat(2,4)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(2,j1))/6._dp
End Do 
End Do 
End Do 
mat(2,4) = mat(2,4)+(g2**2*vu*vL(2))/4._dp
mat(2,5) = 0._dp 
Do j1 = 1,3
mat(2,5) = mat(2,5)+(vd*vu*Conjg(lam(j1))*Yv(3,j1))/2._dp
End Do 
Do j1 = 1,3
mat(2,5) = mat(2,5)-((vR(j1)*Tv(3,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(vu*Conjg(Yv(j2,j1))*vL(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j1,j3,j2))*vR(j2)*vR(j3)*Yv(3,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j2,j1,j3))*vR(j2)*vR(j3)*Yv(3,j1))/6._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
mat(2,5) = mat(2,5)-(Conjg(kap(j2,j3,j1))*vR(j2)*vR(j3)*Yv(3,j1))/6._dp
End Do 
End Do 
End Do 
mat(2,5) = mat(2,5)+(g2**2*vu*vL(3))/4._dp
mat(2,6) = 0._dp 
Do j1 = 1,3
Do j2 = 1,3
mat(2,6) = mat(2,6)-(Conjg(Ye(1,j2))*vL(j2)*vR(j1)*lam(j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,6) = mat(2,6)-(vd*Conjg(Ye(1,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,7) = 0._dp 
Do j1 = 1,3
Do j2 = 1,3
mat(2,7) = mat(2,7)-(Conjg(Ye(2,j2))*vL(j2)*vR(j1)*lam(j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,7) = mat(2,7)-(vd*Conjg(Ye(2,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(2,8) = 0._dp 
Do j1 = 1,3
Do j2 = 1,3
mat(2,8) = mat(2,8)-(Conjg(Ye(3,j2))*vL(j2)*vR(j1)*lam(j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
mat(2,8) = mat(2,8)-(vd*Conjg(Ye(3,j1))*vR(j2)*Yv(j1,j2))/2._dp
End Do 
End Do 
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(g1**2*vdFix**2)/8._dp
mat(3,3) = mat(3,3)-(g2**2*vdFix**2)/8._dp
mat(3,3) = mat(3,3)-(g1**2*vuFix**2)/8._dp
mat(3,3) = mat(3,3)+(g2**2*vuFix**2)/8._dp
mat(3,3) = mat(3,3)+ml2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(3,3) = mat(3,3)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(3,3) = mat(3,3)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,1))/2._dp
End Do 
mat(3,3) = mat(3,3)+(g2**2*vL(1)**2)/4._dp
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+ml2(1,2)
Do j1 = 1,3
Do j2 = 1,3
mat(3,4) = mat(3,4)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,4) = mat(3,4)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,2))/2._dp
End Do 
mat(3,4) = mat(3,4)+(g2**2*vL(1)*vL(2))/4._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+ml2(1,3)
Do j1 = 1,3
Do j2 = 1,3
mat(3,5) = mat(3,5)+(Conjg(Yv(1,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(3,5) = mat(3,5)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,3))/2._dp
End Do 
mat(3,5) = mat(3,5)+(g2**2*vL(1)*vL(3))/4._dp
mat(3,6) = 0._dp 
mat(3,6) = mat(3,6)+(vd*Conjg(Te(1,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,6) = mat(3,6)-(vu*Conjg(Ye(1,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,7) = 0._dp 
mat(3,7) = mat(3,7)+(vd*Conjg(Te(2,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,7) = mat(3,7)-(vu*Conjg(Ye(2,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vd*Conjg(Te(3,1)))/sqrt(2._dp)
Do j1 = 1,3
mat(3,8) = mat(3,8)-(vu*Conjg(Ye(3,1))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+(g1**2*vdFix**2)/8._dp
mat(4,4) = mat(4,4)-(g2**2*vdFix**2)/8._dp
mat(4,4) = mat(4,4)-(g1**2*vuFix**2)/8._dp
mat(4,4) = mat(4,4)+(g2**2*vuFix**2)/8._dp
mat(4,4) = mat(4,4)+ml2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(4,4) = mat(4,4)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(4,4) = mat(4,4)+(vd**2*Conjg(Ye(j1,2))*Ye(j1,2))/2._dp
End Do 
mat(4,4) = mat(4,4)+(g2**2*vL(2)**2)/4._dp
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+ml2(2,3)
Do j1 = 1,3
Do j2 = 1,3
mat(4,5) = mat(4,5)+(Conjg(Yv(2,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(4,5) = mat(4,5)+(vd**2*Conjg(Ye(j1,2))*Ye(j1,3))/2._dp
End Do 
mat(4,5) = mat(4,5)+(g2**2*vL(2)*vL(3))/4._dp
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)+(vd*Conjg(Te(1,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(4,6) = mat(4,6)-(vu*Conjg(Ye(1,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(vd*Conjg(Te(2,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(4,7) = mat(4,7)-(vu*Conjg(Ye(2,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(4,8) = 0._dp 
mat(4,8) = mat(4,8)+(vd*Conjg(Te(3,2)))/sqrt(2._dp)
Do j1 = 1,3
mat(4,8) = mat(4,8)-(vu*Conjg(Ye(3,2))*vR(j1)*lam(j1))/2._dp
End Do 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+(g1**2*vdFix**2)/8._dp
mat(5,5) = mat(5,5)-(g2**2*vdFix**2)/8._dp
mat(5,5) = mat(5,5)-(g1**2*vuFix**2)/8._dp
mat(5,5) = mat(5,5)+(g2**2*vuFix**2)/8._dp
mat(5,5) = mat(5,5)+ml2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(5,5) = mat(5,5)+(Conjg(Yv(3,j2))*vR(j1)*vR(j2)*Yv(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(g1**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)-(g2**2*vL(j1)**2)/8._dp
End Do 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vd**2*Conjg(Ye(j1,3))*Ye(j1,3))/2._dp
End Do 
mat(5,5) = mat(5,5)+(g2**2*vL(3)**2)/4._dp
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(vd*Conjg(Te(1,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(5,6) = mat(5,6)-(vu*Conjg(Ye(1,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)+(vd*Conjg(Te(2,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(5,7) = mat(5,7)-(vu*Conjg(Ye(2,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(5,8) = 0._dp 
mat(5,8) = mat(5,8)+(vd*Conjg(Te(3,3)))/sqrt(2._dp)
Do j1 = 1,3
mat(5,8) = mat(5,8)-(vu*Conjg(Ye(3,3))*vR(j1)*lam(j1))/2._dp
End Do 
mat(6,6) = 0._dp 
mat(6,6) = mat(6,6)-(g1**2*vdFix**2)/4._dp
mat(6,6) = mat(6,6)+(g1**2*vuFix**2)/4._dp
mat(6,6) = mat(6,6)+me2(1,1)
Do j1 = 1,3
Do j2 = 1,3
mat(6,6) = mat(6,6)+(Conjg(Ye(1,j2))*vL(j1)*vL(j2)*Ye(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)-(g1**2*vL(j1)**2)/4._dp
End Do 
Do j1 = 1,3
mat(6,6) = mat(6,6)+(vd**2*Conjg(Ye(1,j1))*Ye(1,j1))/2._dp
End Do 
mat(6,7) = 0._dp 
mat(6,7) = mat(6,7)+me2(1,2)
Do j1 = 1,3
Do j2 = 1,3
mat(6,7) = mat(6,7)+(Conjg(Ye(2,j2))*vL(j1)*vL(j2)*Ye(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,7) = mat(6,7)+(vd**2*Conjg(Ye(2,j1))*Ye(1,j1))/2._dp
End Do 
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)+me2(1,3)
Do j1 = 1,3
Do j2 = 1,3
mat(6,8) = mat(6,8)+(Conjg(Ye(3,j2))*vL(j1)*vL(j2)*Ye(1,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(6,8) = mat(6,8)+(vd**2*Conjg(Ye(3,j1))*Ye(1,j1))/2._dp
End Do 
mat(7,7) = 0._dp 
mat(7,7) = mat(7,7)-(g1**2*vdFix**2)/4._dp
mat(7,7) = mat(7,7)+(g1**2*vuFix**2)/4._dp
mat(7,7) = mat(7,7)+me2(2,2)
Do j1 = 1,3
Do j2 = 1,3
mat(7,7) = mat(7,7)+(Conjg(Ye(2,j2))*vL(j1)*vL(j2)*Ye(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)-(g1**2*vL(j1)**2)/4._dp
End Do 
Do j1 = 1,3
mat(7,7) = mat(7,7)+(vd**2*Conjg(Ye(2,j1))*Ye(2,j1))/2._dp
End Do 
mat(7,8) = 0._dp 
mat(7,8) = mat(7,8)+me2(2,3)
Do j1 = 1,3
Do j2 = 1,3
mat(7,8) = mat(7,8)+(Conjg(Ye(3,j2))*vL(j1)*vL(j2)*Ye(2,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vd**2*Conjg(Ye(3,j1))*Ye(2,j1))/2._dp
End Do 
mat(8,8) = 0._dp 
mat(8,8) = mat(8,8)-(g1**2*vdFix**2)/4._dp
mat(8,8) = mat(8,8)+(g1**2*vuFix**2)/4._dp
mat(8,8) = mat(8,8)+me2(3,3)
Do j1 = 1,3
Do j2 = 1,3
mat(8,8) = mat(8,8)+(Conjg(Ye(3,j2))*vL(j1)*vL(j2)*Ye(3,j1))/2._dp
End Do 
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)-(g1**2*vL(j1)**2)/4._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(vd**2*Conjg(Ye(3,j1))*Ye(3,j1))/2._dp
End Do 

 
 Do i1=2,8
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MHpm2,ZP,ierr,test) 
 
 
! Fix phases
Do i1=1,8
  pos=Maxloc(Abs(ZP(i1,:)),1)
  If (Real(ZP(i1,pos),dp).lt.0._dp) Then
    ZP(i1,:)=-ZP(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,8
  If (MHpm2(i1).ne.MHpm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (MHpm2(i1).Ge.0._dp) Then 
  MHpm(i1)=Sqrt(MHpm2(i1) ) 
  Else 
  MHpm = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMHpmEffPot 

Subroutine CalculateMChiEffPot(g1,g2,lam,Yv,kap,M1,M2,vd,vu,vL,vR,UV,MChi,kont)

Real(dp) ,Intent(in) :: g1,g2,vd,vu,vL(3),vR(3)

Complex(dp) ,Intent(in) :: lam(3),Yv(3,3),kap(3,3,3),M1,M2

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MChi(10) 
Complex(dp), Intent(out) ::  UV(10,10) 
                              
Complex(dp) :: mat(10,10), mat2(10,10), phaseM, E10(10) 

Real(dp) :: UVa(10,10), test(2), eig(10) 

Real(dp) :: MChitemp(10) 
Complex(dp) :: UVtemp(10,10) 
 
Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMChi'
 
mat(1,1) = 0._dp 
mat(1,2) = 0._dp 
mat(1,3) = 0._dp 
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)-(g1*vL(1))/2._dp
mat(1,5) = 0._dp 
mat(1,5) = mat(1,5)+(g2*vL(1))/2._dp
mat(1,6) = 0._dp 
mat(1,7) = 0._dp 
Do j1 = 1,3
mat(1,7) = mat(1,7)+(vR(j1)*Yv(1,j1))/sqrt(2._dp)
End Do 
mat(1,8) = 0._dp 
mat(1,8) = mat(1,8)+(vu*Yv(1,1))/sqrt(2._dp)
mat(1,9) = 0._dp 
mat(1,9) = mat(1,9)+(vu*Yv(1,2))/sqrt(2._dp)
mat(1,10) = 0._dp 
mat(1,10) = mat(1,10)+(vu*Yv(1,3))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,3) = 0._dp 
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)-(g1*vL(2))/2._dp
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(g2*vL(2))/2._dp
mat(2,6) = 0._dp 
mat(2,7) = 0._dp 
Do j1 = 1,3
mat(2,7) = mat(2,7)+(vR(j1)*Yv(2,j1))/sqrt(2._dp)
End Do 
mat(2,8) = 0._dp 
mat(2,8) = mat(2,8)+(vu*Yv(2,1))/sqrt(2._dp)
mat(2,9) = 0._dp 
mat(2,9) = mat(2,9)+(vu*Yv(2,2))/sqrt(2._dp)
mat(2,10) = 0._dp 
mat(2,10) = mat(2,10)+(vu*Yv(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)-(g1*vL(3))/2._dp
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(g2*vL(3))/2._dp
mat(3,6) = 0._dp 
mat(3,7) = 0._dp 
Do j1 = 1,3
mat(3,7) = mat(3,7)+(vR(j1)*Yv(3,j1))/sqrt(2._dp)
End Do 
mat(3,8) = 0._dp 
mat(3,8) = mat(3,8)+(vu*Yv(3,1))/sqrt(2._dp)
mat(3,9) = 0._dp 
mat(3,9) = mat(3,9)+(vu*Yv(3,2))/sqrt(2._dp)
mat(3,10) = 0._dp 
mat(3,10) = mat(3,10)+(vu*Yv(3,3))/sqrt(2._dp)
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+M1
mat(4,5) = 0._dp 
mat(4,6) = 0._dp 
mat(4,6) = mat(4,6)-(g1*vd)/2._dp
mat(4,7) = 0._dp 
mat(4,7) = mat(4,7)+(g1*vu)/2._dp
mat(4,8) = 0._dp 
mat(4,9) = 0._dp 
mat(4,10) = 0._dp 
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+M2
mat(5,6) = 0._dp 
mat(5,6) = mat(5,6)+(g2*vd)/2._dp
mat(5,7) = 0._dp 
mat(5,7) = mat(5,7)-(g2*vu)/2._dp
mat(5,8) = 0._dp 
mat(5,9) = 0._dp 
mat(5,10) = 0._dp 
mat(6,6) = 0._dp 
mat(6,7) = 0._dp 
Do j1 = 1,3
mat(6,7) = mat(6,7)-((vR(j1)*lam(j1))/sqrt(2._dp))
End Do 
mat(6,8) = 0._dp 
mat(6,8) = mat(6,8)-((vu*lam(1))/sqrt(2._dp))
mat(6,9) = 0._dp 
mat(6,9) = mat(6,9)-((vu*lam(2))/sqrt(2._dp))
mat(6,10) = 0._dp 
mat(6,10) = mat(6,10)-((vu*lam(3))/sqrt(2._dp))
mat(7,7) = 0._dp 
mat(7,8) = 0._dp 
Do j1 = 1,3
mat(7,8) = mat(7,8)+(vL(j1)*Yv(j1,1))/sqrt(2._dp)
End Do 
mat(7,8) = mat(7,8)-((vd*lam(1))/sqrt(2._dp))
mat(7,9) = 0._dp 
Do j1 = 1,3
mat(7,9) = mat(7,9)+(vL(j1)*Yv(j1,2))/sqrt(2._dp)
End Do 
mat(7,9) = mat(7,9)-((vd*lam(2))/sqrt(2._dp))
mat(7,10) = 0._dp 
Do j1 = 1,3
mat(7,10) = mat(7,10)+(vL(j1)*Yv(j1,3))/sqrt(2._dp)
End Do 
mat(7,10) = mat(7,10)-((vd*lam(3))/sqrt(2._dp))
mat(8,8) = 0._dp 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(sqrt(2._dp)*vR(j1)*kap(1,1,j1))/3._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(sqrt(2._dp)*vR(j1)*kap(1,j1,1))/3._dp
End Do 
Do j1 = 1,3
mat(8,8) = mat(8,8)+(sqrt(2._dp)*vR(j1)*kap(j1,1,1))/3._dp
End Do 
mat(8,9) = 0._dp 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(1,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(1,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(2,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(2,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(j1,1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,9) = mat(8,9)+(vR(j1)*kap(j1,2,1))/(3._dp*sqrt(2._dp))
End Do 
mat(8,10) = 0._dp 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(1,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(1,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(3,1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(3,j1,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(j1,1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(8,10) = mat(8,10)+(vR(j1)*kap(j1,3,1))/(3._dp*sqrt(2._dp))
End Do 
mat(9,9) = 0._dp 
Do j1 = 1,3
mat(9,9) = mat(9,9)+(sqrt(2._dp)*vR(j1)*kap(2,2,j1))/3._dp
End Do 
Do j1 = 1,3
mat(9,9) = mat(9,9)+(sqrt(2._dp)*vR(j1)*kap(2,j1,2))/3._dp
End Do 
Do j1 = 1,3
mat(9,9) = mat(9,9)+(sqrt(2._dp)*vR(j1)*kap(j1,2,2))/3._dp
End Do 
mat(9,10) = 0._dp 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(2,3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(2,j1,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(3,2,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(3,j1,2))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(j1,2,3))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
mat(9,10) = mat(9,10)+(vR(j1)*kap(j1,3,2))/(3._dp*sqrt(2._dp))
End Do 
mat(10,10) = 0._dp 
Do j1 = 1,3
mat(10,10) = mat(10,10)+(sqrt(2._dp)*vR(j1)*kap(3,3,j1))/3._dp
End Do 
Do j1 = 1,3
mat(10,10) = mat(10,10)+(sqrt(2._dp)*vR(j1)*kap(3,j1,3))/3._dp
End Do 
Do j1 = 1,3
mat(10,10) = mat(10,10)+(sqrt(2._dp)*vR(j1)*kap(j1,3,3))/3._dp
End Do 

 
 Do i1=2,10
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
If (Maxval(Abs(Aimag(mat))).Eq.0._dp) Then 
Call EigenSystemQP(Real(mat,dp),Eig,UVa,ierr,test) 
 
   Do i1=1,10
   If ((Eig(i1).Lt.0._dp).or.(Abs(eig(i1)).lt.1E-15)) Then 
    MChi(i1) = - Eig(i1) 
    UV(i1,:) = (0._dp,1._dp)*UVa(i1,:) 
   Else 
    MChi(i1) = Eig(i1) 
    UV(i1,:) = UVa(i1,:)
    End If 
   End Do 
 
  Do i1=1,10
   pos=Maxloc(Abs(UV(i1,:)),1) 
   If (Abs(Real(UV(i1,pos),dp)).gt.Abs(Aimag(UV(i1,pos)))) Then 
      If (Real(UV(i1,pos),dp).lt.0._dp) Then 
        UV(i1,:)=-UV(i1,:) 
       End If 
    Else 
      If (Aimag(UV(i1,pos)).lt.0._dp) Then 
        UV(i1,:)=-UV(i1,:) 
      End If 
    End If 
 End Do 
 
Do i1=1,9
  Do i2=i1+1,10
    If (MChi(i1).Gt.MChi(i2)) Then 
      Eig(1) = MChi(i1) 
      MChi(i1) = MChi(i2) 
      MChi(i2) =  Eig(1) 
      E10 = UV(i1,:) 
      UV(i1,:) = UV(i2,:) 
      UV(i2,:) = E10
    End If 
   End Do 
End Do 
 
Else 
 
mat2 = Matmul( Transpose(Conjg( mat) ), mat ) 
Call EigensystemQP(mat2, Eig, UV, ierr, test) 
mat2 = Matmul( Conjg(UV), Matmul( mat, Transpose( Conjg(UV)))) 
Do i1=1,10
  If (Eig(i1).ne.Eig(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  UV(i1,:)= phaseM * UV(i1,:) 
End if 
  If ((Abs(Eig(i1)).Le.MaxMassNumericalZero).and.(Eig(i1).lt.0._dp)) Eig(i1) = Abs(Eig(i1))+1.E-10_dp 
  If (Eig(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
MChi = Sqrt( Eig ) 
 
End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMChi, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


If ((Abs(UV(1,2)).gt.Abs(UV(2,1))).And.(MChi(1).lt.1.0E-15_dp).And.(MChi(2).lt.1.0E-15_dp)) Then 

   MChitemp = MChi 

   UVtemp = UV 

   UV(1,:) = UVtemp(2,:) 

   UV(2,:) = UVtemp(1,:) 

   MChi(1) = MChitemp(2) 

   MChi(2) = MChitemp(1) 

End If 
 
Iname = Iname - 1 
 
End Subroutine CalculateMChiEffPot 

Subroutine CalculateMChaEffPot(g2,Ye,lam,Yv,M2,vd,vu,vL,vR,ZER,ZEL,MCha,kont)

Real(dp),Intent(in) :: g2,vd,vu,vL(3),vR(3)

Complex(dp),Intent(in) :: Ye(3,3),lam(3),Yv(3,3),M2

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MCha(5) 
 Complex(dp), Intent(out) :: ZER(5,5), ZEL(5,5) 
 
 Complex(dp) :: mat(5,5)=0._dp, mat2(5,5)=0._dp, phaseM 

Complex(dp) :: ZER2(5,5), ZEL2(5,5) 
 
 Real(dp) :: ZER1(5,5), ZEL1(5,5), test(2), MCha2(5) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMCha'
 
MCha = 0._dp 
ZER = 0._dp 
ZEL = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(vd*Ye(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vd*Ye(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(vd*Ye(3,1))/sqrt(2._dp)
mat(1,4) = 0._dp 
mat(1,4) = mat(1,4)+(g2*vL(1))/sqrt(2._dp)
mat(1,5) = 0._dp 
Do j1 = 1,3
mat(1,5) = mat(1,5)-((vR(j1)*Yv(1,j1))/sqrt(2._dp))
End Do 
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vd*Ye(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(vd*Ye(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(vd*Ye(3,2))/sqrt(2._dp)
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(g2*vL(2))/sqrt(2._dp)
mat(2,5) = 0._dp 
Do j1 = 1,3
mat(2,5) = mat(2,5)-((vR(j1)*Yv(2,j1))/sqrt(2._dp))
End Do 
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(vd*Ye(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(vd*Ye(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(vd*Ye(3,3))/sqrt(2._dp)
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(g2*vL(3))/sqrt(2._dp)
mat(3,5) = 0._dp 
Do j1 = 1,3
mat(3,5) = mat(3,5)-((vR(j1)*Yv(3,j1))/sqrt(2._dp))
End Do 
mat(4,1) = 0._dp 
mat(4,2) = 0._dp 
mat(4,3) = 0._dp 
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+M2
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+(g2*vu)/sqrt(2._dp)
mat(5,1) = 0._dp 
Do j1 = 1,3
mat(5,1) = mat(5,1)-((vL(j1)*Ye(1,j1))/sqrt(2._dp))
End Do 
mat(5,2) = 0._dp 
Do j1 = 1,3
mat(5,2) = mat(5,2)-((vL(j1)*Ye(2,j1))/sqrt(2._dp))
End Do 
mat(5,3) = 0._dp 
Do j1 = 1,3
mat(5,3) = mat(5,3)-((vL(j1)*Ye(3,j1))/sqrt(2._dp))
End Do 
mat(5,4) = 0._dp 
mat(5,4) = mat(5,4)+(g2*vd)/sqrt(2._dp)
mat(5,5) = 0._dp 
Do j1 = 1,3
mat(5,5) = mat(5,5)+(vR(j1)*lam(j1))/sqrt(2._dp)
End Do 

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2,ZEL1,ierr,test) 
ZEL2 = ZEL1 
Else 
Call EigenSystem(mat2,MCha2,ZEL2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MCha2,ZER1,ierr,test) 
                  
                  
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MCha2,ZER2,ierr,test) 
 
 
End If 
ZER2 = Conjg(ZER2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZER2),mat),Transpose( Conjg(ZEL2))) 
Do i1=1,5
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZEL2(i1,:) = phaseM *ZEL2(i1,:) 
 End if 
End Do 
 
Do i1=1,5
If (Abs(ZEL2(i1,i1)).gt.0._dp) Then 
phaseM = ZEL2(i1,i1) / Abs(ZEL2(i1,i1)) 
ZEL2(i1,:) = Conjg(phaseM) *ZEL2(i1,:) 
 ZER2(i1,:) = phaseM *ZER2(i1,:) 
 End if 
  If (MCha2(i1).ne.MCha2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Abs(MCha2(i1)).Le.MaxMassNumericalZero) MCha2(i1) = Abs(MCha2(i1))+1.E-10_dp 
  If (MCha2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMCha, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MCha = Sqrt( MCha2 ) 
ZER = ZER2 
ZEL = ZEL2 
Iname = Iname - 1 
 
End Subroutine CalculateMChaEffPot 

Subroutine CalculateMFdEffPot(Yd,vd,ZDL,ZDR,MFd,kont)

Real(dp),Intent(in) :: vd

Complex(dp),Intent(in) :: Yd(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFd(3) 
 Complex(dp), Intent(out) :: ZDL(3,3), ZDR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZDL2(3,3), ZDR2(3,3) 
 
 Real(dp) :: ZDL1(3,3), ZDR1(3,3), test(2), MFd2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFd'
 
MFd = 0._dp 
ZDL = 0._dp 
ZDR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(vd*Yd(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vd*Yd(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(vd*Yd(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vd*Yd(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(vd*Yd(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(vd*Yd(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(vd*Yd(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(vd*Yd(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(vd*Yd(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2,ZDR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFd2,ZDL1,ierr,test) 
                  
                  
ZDL2 = ZDL1 
Else 
Call EigenSystem(mat2,MFd2,ZDL2,ierr,test) 
 
 
End If 
ZDL2 = Conjg(ZDL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZDL2),mat),Transpose( Conjg(ZDR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZDR2(i1,:) = phaseM *ZDR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZDR2(i1,i1)).gt.0._dp) Then 
phaseM = ZDR2(i1,i1) / Abs(ZDR2(i1,i1)) 
ZDR2(i1,:) = Conjg(phaseM) *ZDR2(i1,:) 
 ZDL2(i1,:) = phaseM *ZDL2(i1,:) 
 End if 
  If (MFd2(i1).ne.MFd2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Abs(MFd2(i1)).Le.MaxMassNumericalZero) MFd2(i1) = Abs(MFd2(i1))+1.E-10_dp 
  If (MFd2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFd, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFd = Sqrt( MFd2 ) 
ZDL = ZDL2 
ZDR = ZDR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFdEffPot 

Subroutine CalculateMFuEffPot(Yu,vu,ZUL,ZUR,MFu,kont)

Real(dp),Intent(in) :: vu

Complex(dp),Intent(in) :: Yu(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFu(3) 
 Complex(dp), Intent(out) :: ZUL(3,3), ZUR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZUL2(3,3), ZUR2(3,3) 
 
 Real(dp) :: ZUL1(3,3), ZUR1(3,3), test(2), MFu2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFu'
 
MFu = 0._dp 
ZUL = 0._dp 
ZUR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(vu*Yu(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vu*Yu(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(vu*Yu(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vu*Yu(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(vu*Yu(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(vu*Yu(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(vu*Yu(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(vu*Yu(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(vu*Yu(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2,ZUR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFu2,ZUL1,ierr,test) 
                  
                  
ZUL2 = ZUL1 
Else 
Call EigenSystem(mat2,MFu2,ZUL2,ierr,test) 
 
 
End If 
ZUL2 = Conjg(ZUL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZUL2),mat),Transpose( Conjg(ZUR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZUR2(i1,:) = phaseM *ZUR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZUR2(i1,i1)).gt.0._dp) Then 
phaseM = ZUR2(i1,i1) / Abs(ZUR2(i1,i1)) 
ZUR2(i1,:) = Conjg(phaseM) *ZUR2(i1,:) 
 ZUL2(i1,:) = phaseM *ZUL2(i1,:) 
 End if 
  If (MFu2(i1).ne.MFu2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (Abs(MFu2(i1)).Le.MaxMassNumericalZero) MFu2(i1) = Abs(MFu2(i1))+1.E-10_dp 
  If (MFu2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFu, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFu = Sqrt( MFu2 ) 
ZUL = ZUL2 
ZUR = ZUR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFuEffPot 

Subroutine CalculateVPVZEffPot(g1,g2,vd,vu,vL,ZZ,MVZ,MVZ2,TW,kont)

Real(dp), Intent(in) :: g1,g2,vd,vu,vL(3)

Real(dp), Intent(out) :: TW

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVZ, MVZ2
Real(dp) :: VPVZ2(2),VPVZ(2)  

Real(dp), Intent(out) :: ZZ(2,2) 
 
Real(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVPVZ'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g1**2*vdFix**2)/4._dp
mat(1,1) = mat(1,1)+(g1**2*vuFix**2)/4._dp
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g1**2*vL(j1)**2)/4._dp
End Do 
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g1*g2*vdFix**2)/4._dp
mat(1,2) = mat(1,2)-(g1*g2*vuFix**2)/4._dp
Do j1 = 1,3
mat(1,2) = mat(1,2)-(g1*g2*vL(j1)**2)/4._dp
End Do 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*vdFix**2)/4._dp
mat(2,2) = mat(2,2)+(g2**2*vuFix**2)/4._dp
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/4._dp
End Do 

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,VPVZ2,ZZ,ierr,test) 
 
 
ZZ = Transpose(ZZ) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(VPVZ2(i1)).Le.1.E-10_dp*(Maxval(VPVZ2))) VPVZ2(i1) = 1.E-10_dp 
  If (VPVZ2(i1).ne.VPVZ2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (VPVZ2(i1).Ge.0._dp) Then 
  VPVZ(i1) =Sqrt(VPVZ2(i1) ) 
  Else 
  VPVZ(i1)= 1._dp 
  VPVZ(i1)= 1._dp 
! kont = -104 
 End if 
End Do 
 
MVZ = VPVZ(2) 
MVZ2 = VPVZ2(2) 

 Iname = Iname - 1 
 
End Subroutine CalculateVPVZEffPot 

Subroutine CalculateVWmEffPot(g2,vd,vu,vL,ZW,MVWm,MVWm2,kont)

Real(dp), Intent(in) :: g2,vd,vu,vL(3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVWm, MVWm2
Real(dp) :: VWm2(2),VWm(2)  

Complex(dp), Intent(out) :: ZW(2,2) 
 
Complex(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVWm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g2**2*vdFix**2)/4._dp
mat(1,1) = mat(1,1)+(g2**2*vuFix**2)/4._dp
Do j1 = 1,3
mat(1,1) = mat(1,1)+(g2**2*vL(j1)**2)/4._dp
End Do 
mat(1,2) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*vdFix**2)/4._dp
mat(2,2) = mat(2,2)+(g2**2*vuFix**2)/4._dp
Do j1 = 1,3
mat(2,2) = mat(2,2)+(g2**2*vL(j1)**2)/4._dp
End Do 

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,VWm2,ZW,ierr,test) 
 
 
ZW = Transpose(ZW) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
    return
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(VWm2(i1)).Le.1.E-10_dp*(Maxval(VWm2))) VWm2(i1) = 1.E-10_dp 
  If (VWm2(i1).ne.VWm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
      return
    End If 
  If (VWm2(i1).Ge.0._dp) Then 
  VWm(i1) =Sqrt(VWm2(i1) ) 
  Else 
  VWm(i1)= 1._dp 
  VWm(i1)= 1._dp 
! kont = -104 
 End if 
End Do 
 
MVWm = VWm(1) 
MVWm2 = VWm2(1) 

 Iname = Iname - 1 
 
End Subroutine CalculateVWmEffPot 

Subroutine TreeMassesSM(MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MVWm,MVWm2,           & 
& MVZ,MVZ2,TW,ZER,ZEL,ZDL,ZDR,UV,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,               & 
& lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,             & 
& M1,M2,M3,GenerationMixing,kont)

Implicit None 
 
Real(dp),Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Real(dp),Intent(out) :: MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),MFu2(3),MVWm,               & 
& MVWm2,MVZ,MVZ2,TW,ZZ(2,2)

Complex(dp),Intent(out) :: ZER(5,5),ZEL(5,5),ZDL(3,3),ZDR(3,3),UV(10,10),ZUL(3,3),ZUR(3,3),ZW(2,2)

Real(dp),Intent(in) :: vd,vu,vR(3),vL(3)

Logical, Intent(in) :: GenerationMixing 
Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,kontSave 
Iname = Iname + 1 
NameOfUnit(Iname) = 'TreeMassesmunuSSM'
 
kont = 0 
Call CalculateMChi(g1,g2,lam,Yv,kap,M1,M2,vd,vu,vL,vR,UV,MChi,kont)

MChi2 = MChi**2 
Call CalculateMCha(g2,Ye,lam,Yv,M2,vd,vu,vL,vR,ZER,ZEL,MCha,kont)

MCha2 = MCha**2 
Call CalculateMFd(Yd,vd,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
Call CalculateMFu(Yu,vu,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 

 
 Call CalculateVPVZ(g1,g2,vd,vu,vL,ZZ,MVZ,MVZ2,TW,kont)

Call CalculateVWm(g2,vd,vu,vL,ZW,MVWm,MVWm2,kont)

Iname = Iname - 1 
 
End Subroutine  TreeMassesSM 
 
 
Subroutine SortGoldstones(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,           & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,               & 
& TW,ZER,ZEL,ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,kont)

Real(dp),Intent(inout) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp),Intent(inout) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Integer, Intent(inout) :: kont 
Integer :: i1, i2, pos 
Real(dp) :: MAhtemp(8) 
Complex(dp) :: ZAhtemp(8,8) 
Real(dp) :: MHpmtemp(8) 
Complex(dp) :: ZHpmtemp(8,8) 


pos = MinLoc(Abs(MAh2-MVZ2*RXiZ),1) 
If (pos.ne.1) Then 
  MAhtemp = MAh2 
  ZAhtemp = ZA 
  MAh2(1) = MAhtemp(pos) 
  ZA(1,:) = ZAhtemp(pos,:) 
  MAh2(pos) = MAhtemp(1) 
  ZA(pos,:) = ZAhtemp(1,:) 
End if 

 ! Reorder the physical states by their mass 
Do i1=2,8
 pos = Minloc(MAh2(i1:8),1) + i1 -1  
If (pos.ne.i1) Then 
  MAhtemp = MAh2 
  ZAhtemp = ZA 
  MAh2(i1) = MAhtemp(pos) 
  ZA(i1,:) = ZAhtemp(pos,:) 
  MAh2(pos) = MAhtemp(i1) 
  ZA(pos,:) = ZAhtemp(i1,:) 
End if 
End do 
MAh = sqrt(MAh2) 

 
 
pos = MinLoc(Abs(MHpm2-MVWm2*RXiWm),1) 
If (pos.ne.1) Then 
  MHpmtemp = MHpm2 
  ZHpmtemp = ZP 
  MHpm2(1) = MHpmtemp(pos) 
  ZP(1,:) = ZHpmtemp(pos,:) 
  MHpm2(pos) = MHpmtemp(1) 
  ZP(pos,:) = ZHpmtemp(1,:) 
End if 

 ! Reorder the physical states by their mass 
Do i1=2,8
 pos = Minloc(MHpm2(i1:8),1) + i1 -1  
If (pos.ne.i1) Then 
  MHpmtemp = MHpm2 
  ZHpmtemp = ZP 
  MHpm2(i1) = MHpmtemp(pos) 
  ZP(i1,:) = ZHpmtemp(pos,:) 
  MHpm2(pos) = MHpmtemp(i1) 
  ZP(pos,:) = ZHpmtemp(i1,:) 
End if 
End do 
MHpm = sqrt(MHpm2) 

 
 
End subroutine SortGoldstones 


End Module SusyMasses_munuSSM3G 
 
