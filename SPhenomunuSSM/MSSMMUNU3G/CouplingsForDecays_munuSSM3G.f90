! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:42 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module CouplingsForDecays_munuSSM3G
 
Use Control 
Use Model_Data_munuSSM3G 
Use RGEs_munuSSM3G 
Use Couplings_munuSSM3G 
Use LoopCouplings_munuSSM3G 
Use Tadpoles_munuSSM3G 
 Use SusyMasses_munuSSM3G 
Use Mathematics, Only: CompareMatrices, Adjungate 
 
Use StandardModel 
Contains 
 
 
 
Subroutine CouplingsFor_Sd_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplAhSdcSd,cplChaFucSdL,cplChaFucSdR,cplChiFdcSdL,cplChiFdcSdR,cplGluFdcSdL,           & 
& cplGluFdcSdR,cplhhSdcSd,cplHpmSucSd,cplSdcSdVZ,cplSucSdVWm,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplAhSdcSd(8,6,6),cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplChiFdcSdL(10,3,6),       & 
& cplChiFdcSdR(10,3,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplhhSdcSd(8,6,6),            & 
& cplHpmSucSd(8,6,6),cplSdcSdVZ(6,6),cplSucSdVWm(6,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Sd_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSdcSdT(gt1,gt2,gt3,Yd,Td,lam,vu,vR,ZD,ZA,cplAhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSdcSdT(gt1,gt2,gt3,g1,g2,Yd,Td,lam,vd,vu,vL,vR,ZD,ZH,cplhhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplHpmSucSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingHpmSucSdT(gt1,gt2,gt3,g2,Yd,Td,Ye,lam,Yv,Yu,Tu,vd,vu,vL,vR,              & 
& ZD,ZU,ZP,cplHpmSucSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcSdVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSdcSdVZT(gt1,gt2,g1,g2,ZD,TW,cplSdcSdVZ(gt1,gt2))

 End Do 
End Do 


cplSucSdVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSucSdVWmT(gt1,gt2,g2,ZD,ZU,cplSucSdVWm(gt1,gt2))

 End Do 
End Do 


cplChaFucSdL = 0._dp 
cplChaFucSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChaFucSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplChaFucSdL(gt1,gt2,gt3)& 
& ,cplChaFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Sd_decays_2B
 
Subroutine CouplingsFor_Su_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplAhSucSu,cplChiFucSuL,cplChiFucSuR,cplcChaFdcSuL,cplcChaFdcSuR,cplGluFucSuL,         & 
& cplGluFucSuR,cplhhSucSu,cplSdcHpmcSu,cplSdcSucVWm,cplSucSuVZ,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplAhSucSu(8,6,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),cplcChaFdcSuL(5,3,6),     & 
& cplcChaFdcSuR(5,3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6),cplhhSucSu(8,6,6),            & 
& cplSdcHpmcSu(6,8,6),cplSdcSucVWm(6,6),cplSucSuVZ(6,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Su_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSucSuT(gt1,gt2,gt3,lam,Yv,Yu,Tu,vd,vL,vR,ZU,ZA,cplAhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSucSuT(gt1,gt2,gt3,g1,g2,lam,Yv,Yu,Tu,vd,vu,vL,vR,ZU,ZH,               & 
& cplhhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcHpmcSu = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 8
  Do gt3 = 1, 6
Call CouplingSdcHpmcSuT(gt1,gt2,gt3,g2,Yd,Td,Ye,lam,Yv,Yu,Tu,vd,vu,vL,vR,             & 
& ZD,ZU,ZP,cplSdcHpmcSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcSucVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSdcSucVWmT(gt1,gt2,g2,ZD,ZU,cplSdcSucVWm(gt1,gt2))

 End Do 
End Do 


cplSucSuVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSucSuVZT(gt1,gt2,g1,g2,ZU,TW,cplSucSuVZ(gt1,gt2))

 End Do 
End Do 


cplChiFucSuL = 0._dp 
cplChiFucSuR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFucSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplChiFucSuL(gt1,gt2,gt3)   & 
& ,cplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaFdcSuL = 0._dp 
cplcChaFdcSuR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChaFdcSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcChaFdcSuL(gt1,gt2,gt3)& 
& ,cplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFucSuL = 0._dp 
cplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFucSuT(gt2,gt3,g3,pG,ZU,ZUL,ZUR,cplGluFucSuL(gt2,gt3),cplGluFucSuR(gt2,gt3))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Su_decays_2B
 
Subroutine CouplingsFor_hh_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplHiggsPP,cplHiggsGG,cplHiggsZZvirt,cplHiggsWWvirt,cplAhAhhh,cplAhhhhh,               & 
& cplAhhhVZ,cplcChaChahhL,cplcChaChahhR,cplChiChihhL,cplChiChihhR,cplcFdFdhhL,           & 
& cplcFdFdhhR,cplcFuFuhhL,cplcFuFuhhR,cplhhhhhh,cplhhHpmcHpm,cplhhHpmcVWm,               & 
& cplhhSdcSd,cplhhSucSu,cplhhcVWmVWm,cplhhVZVZ,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplHiggsPP(8),cplHiggsGG(8),cplHiggsZZvirt(8),cplHiggsWWvirt(8),cplAhAhhh(8,8,8),     & 
& cplAhhhhh(8,8,8),cplAhhhVZ(8,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),             & 
& cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),     & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplhhhhhh(8,8,8),cplhhHpmcHpm(8,8,8),            & 
& cplhhHpmcVWm(8,8),cplhhSdcSd(8,6,6),cplhhSucSu(8,6,6),cplhhcVWmVWm(8),cplhhVZVZ(8)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Complex(dp) :: ratCha(5),ratFd(3),ratFu(3),ratHpm(8),ratSd(6),ratSu(6),ratVWm

Complex(dp) :: ratPCha(5),ratPFd(3),ratPFu(3),ratPHpm(8),ratPSd(6),ratPSu(6),ratPVWm

Complex(dp) :: coup 
Real(dp) :: vev, rq, rsq 
Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_hh_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
cplAhAhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhAhhhT(gt1,gt2,gt3,g1,g2,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,             & 
& ZH,ZA,cplAhAhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhhhhhT(gt1,gt2,gt3,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZH,ZA,             & 
& cplAhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhhhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplinghhhhhhT(gt1,gt2,gt3,g1,g2,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,             & 
& ZH,cplhhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhHpmcHpm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplinghhHpmcHpmT(gt1,gt2,gt3,g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,vd,vu,             & 
& vL,vR,ZH,ZP,cplhhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSdcSdT(gt1,gt2,gt3,g1,g2,Yd,Td,lam,vd,vu,vL,vR,ZD,ZH,cplhhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSucSuT(gt1,gt2,gt3,g1,g2,lam,Yv,Yu,Tu,vd,vu,vL,vR,ZU,ZH,               & 
& cplhhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhVZ = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhhhVZT(gt1,gt2,g1,g2,ZH,ZA,TW,cplAhhhVZ(gt1,gt2))

 End Do 
End Do 


cplhhHpmcVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplinghhHpmcVWmT(gt1,gt2,g2,ZH,ZP,cplhhHpmcVWm(gt1,gt2))

 End Do 
End Do 


cplhhcVWmVWm = 0._dp 
Do gt1 = 1, 8
Call CouplinghhcVWmVWmT(gt1,g2,vd,vu,vL,ZH,cplhhcVWmVWm(gt1))

End Do 


cplhhVZVZ = 0._dp 
Do gt1 = 1, 8
Call CouplinghhVZVZT(gt1,g1,g2,vd,vu,vL,ZH,TW,cplhhVZVZ(gt1))

End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdhhL = 0._dp 
cplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdhhT(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,cplcFdFdhhL(gt1,gt2,gt3)              & 
& ,cplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


vev = Sqrt(vd**2 + vu**2 + vL(1)**2 + vL(2)**2 + vL(3)**2)
cplHiggsWWvirt = cplhhcVWmVWm/vev 
cplHiggsZZvirt = cplhhVZVZ*Sqrt(7._dp/12._dp-10._dp/9._dp*Sin(TW)**2+40._dp/27._dp*Sin(TW)**4)/vev 
 

!----------------------------------------------------
! Scalar Higgs coupling ratios 
!----------------------------------------------------
 
Do i2=1, 5
ratCha(i2) = cplcChaChahhL(i2,i2,i1)*1._dp*vev/MCha(i2) 
End Do 
Do i2=1, 3
ratFd(i2) = cplcFdFdhhL(i2,i2,i1)*1._dp*vev/MFd(i2) 
End Do 
Do i2=1, 3
ratFu(i2) = cplcFuFuhhL(i2,i2,i1)*1._dp*vev/MFu(i2) 
End Do 
Do i2=1, 8
ratHpm(i2) = 0.5_dp*cplhhHpmcHpm(i1,i2,i2)*vev/MHpm2(i2) 
End Do 
Do i2=1, 6
ratSd(i2) = 0.5_dp*cplhhSdcSd(i1,i2,i2)*vev/MSd2(i2) 
End Do 
Do i2=1, 6
ratSu(i2) = 0.5_dp*cplhhSucSu(i1,i2,i2)*vev/MSu2(i2) 
End Do 
ratVWm = -0.5_dp*cplhhcVWmVWm(i1)*vev/MVWm2 
rq = 1._dp - oo4pi2*g3**2 
rsq = 1._dp + 2._dp*oo3pi2*g3**2
Call CoupHiggsToPhoton(m_in,i1,ratCha,ratFd,ratFu,ratHpm,ratSd,ratSu,ratVWm,          & 
& MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,rq,rsq,coup)

cplHiggsPP(i1) = coup*Alpha 
CoupHPP(i1) = coup 
Call CoupHiggsToPhotonSM(m_in,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,rq,coup)

ratioPP(i1) = Abs(cplHiggsPP(i1)/(coup*Alpha))**2 
rq = rsq 
Call CoupHiggsToGluon(m_in,i1,ratFd,ratFu,ratSd,ratSu,MFd,MFu,MSd,MSu,rq,rsq,coup)

cplHiggsGG(i1) = coup*AlphaS_MZ 
CoupHGG(i1) = coup 
Call CoupHiggsToGluonSM(m_in,MFd,MFu,MSd,MSu,rq,coup)

coup = coup*Sqrt(1._dp + 12._dp*oo48pi2*(95._dp/4._dp - 7._dp/6._dp*NFlav(m_in)) *g3**2) 
cplHiggsGG(i1) = cplHiggsGG(i1)*Sqrt(1._dp + 12._dp*oo48pi2*(95._dp/4._dp - 7._dp/6._dp*NFlav(m_in)) *g3**2) 
ratioGG(i1) = Abs(cplHiggsGG(i1)/(coup*AlphaS_MZ))**2 
!----------------------------------------------------
! Coupling ratios for HiggsBounds 
!----------------------------------------------------
 
Do i2=1, 5
rHB_S_S_Cha(i1,i2) = Abs((cplcChaChahhL(i2,i2,i1)+cplcChaChahhR(i2,i2,i1))*vev/(2._dp*MCha(i2)))**2 
rHB_S_P_Cha(i1,i2) = Abs((cplcChaChahhL(i2,i2,i1)-cplcChaChahhR(i2,i2,i1))*vev/(2._dp*MCha(i2)))**2 
End Do 
Do i2=1, 3
rHB_S_S_Fu(i1,i2) = Abs((cplcFuFuhhL(i2,i2,i1)+cplcFuFuhhR(i2,i2,i1))*vev/(2._dp*MFu(i2)))**2 
rHB_S_P_Fu(i1,i2) = Abs((cplcFuFuhhL(i2,i2,i1)-cplcFuFuhhR(i2,i2,i1))*vev/(2._dp*MFu(i2)))**2 
End Do 
Do i2=1, 3
rHB_S_S_Fd(i1,i2) = Abs((cplcFdFdhhL(i2,i2,i1)+cplcFdFdhhR(i2,i2,i1))*vev/(2._dp*MFd(i2)))**2 
rHB_S_P_Fd(i1,i2) = Abs((cplcFdFdhhL(i2,i2,i1)-cplcFdFdhhR(i2,i2,i1))*vev/(2._dp*MFd(i2)))**2 
End Do 
rHB_S_VZ(i1) = Abs(-0.5_dp*cplhhVZVZ(i1)*vev/MVZ2)**2 
rHB_S_VWm(i1) = Abs(-0.5_dp*cplhhcVWmVWm(i1)*vev/MVWm2)**2 
Do i2=1, 10
rHB_S_S_Chi(i1,i2) = Abs((cplChiChihhL(i2,i2,i1)+cplChiChihhR(i2,i2,i1))*vev/(2._dp*MChi(i2)))**2 
rHB_S_P_Chi(i1,i2) = Abs((cplChiChihhL(i2,i2,i1)-cplChiChihhR(i2,i2,i1))*vev/(2._dp*MChi(i2)))**2 
End Do 
If (i1.eq.1) Then 
CPL_A_H_Z = Abs(cplAhhhVZ**2/(g2**2/(cos(TW)*4._dp)))
CPL_H_H_Z = 0._dp 
End if 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhAhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhAhhhT(gt1,gt2,gt3,g1,g2,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,             & 
& ZH,ZA,cplAhAhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhhhhhT(gt1,gt2,gt3,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZH,ZA,             & 
& cplAhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhhhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplinghhhhhhT(gt1,gt2,gt3,g1,g2,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,             & 
& ZH,cplhhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhHpmcHpm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplinghhHpmcHpmT(gt1,gt2,gt3,g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,vd,vu,             & 
& vL,vR,ZH,ZP,cplhhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSdcSdT(gt1,gt2,gt3,g1,g2,Yd,Td,lam,vd,vu,vL,vR,ZD,ZH,cplhhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSucSuT(gt1,gt2,gt3,g1,g2,lam,Yv,Yu,Tu,vd,vu,vL,vR,ZU,ZH,               & 
& cplhhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhVZ = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhhhVZT(gt1,gt2,g1,g2,ZH,ZA,TW,cplAhhhVZ(gt1,gt2))

 End Do 
End Do 


cplhhHpmcVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplinghhHpmcVWmT(gt1,gt2,g2,ZH,ZP,cplhhHpmcVWm(gt1,gt2))

 End Do 
End Do 


cplhhcVWmVWm = 0._dp 
Do gt1 = 1, 8
Call CouplinghhcVWmVWmT(gt1,g2,vd,vu,vL,ZH,cplhhcVWmVWm(gt1))

End Do 


cplhhVZVZ = 0._dp 
Do gt1 = 1, 8
Call CouplinghhVZVZT(gt1,g1,g2,vd,vu,vL,ZH,TW,cplhhVZVZ(gt1))

End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdhhL = 0._dp 
cplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdhhT(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,cplcFdFdhhL(gt1,gt2,gt3)              & 
& ,cplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_hh_decays_2B
 
Subroutine CouplingsFor_Ah_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplPseudoHiggsPP,cplPseudoHiggsGG,cplAhAhAh,cplAhAhhh,cplcChaChaAhL,cplcChaChaAhR,     & 
& cplChiChiAhL,cplChiChiAhR,cplcFdFdAhL,cplcFdFdAhR,cplcFuFuAhL,cplcFuFuAhR,             & 
& cplAhhhhh,cplAhhhVZ,cplAhHpmcHpm,cplAhHpmcVWm,cplAhSdcSd,cplAhSucSu,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplPseudoHiggsPP(8),cplPseudoHiggsGG(8),cplAhAhAh(8,8,8),cplAhAhhh(8,8,8),            & 
& cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8), & 
& cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),           & 
& cplAhhhhh(8,8,8),cplAhhhVZ(8,8),cplAhHpmcHpm(8,8,8),cplAhHpmcVWm(8,8),cplAhSdcSd(8,6,6),& 
& cplAhSucSu(8,6,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Complex(dp) :: ratCha(5),ratFd(3),ratFu(3),ratHpm(8),ratSd(6),ratSu(6),ratVWm

Complex(dp) :: ratPCha(5),ratPFd(3),ratPFu(3),ratPHpm(8),ratPSd(6),ratPSu(6),ratPVWm

Complex(dp) :: coup 
Real(dp) :: vev, rq, rsq 
Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Ah_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
cplAhAhAh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhAhAhT(gt1,gt2,gt3,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZA,cplAhAhAh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhAhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhAhhhT(gt1,gt2,gt3,g1,g2,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,             & 
& ZH,ZA,cplAhAhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhhhhhT(gt1,gt2,gt3,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZH,ZA,             & 
& cplAhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhHpmcHpm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhHpmcHpmT(gt1,gt2,gt3,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,vd,vu,vL,             & 
& vR,ZA,ZP,cplAhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSdcSdT(gt1,gt2,gt3,Yd,Td,lam,vu,vR,ZD,ZA,cplAhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSucSuT(gt1,gt2,gt3,lam,Yv,Yu,Tu,vd,vL,vR,ZU,ZA,cplAhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhVZ = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhhhVZT(gt1,gt2,g1,g2,ZH,ZA,TW,cplAhhhVZ(gt1,gt2))

 End Do 
End Do 


cplAhHpmcVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhHpmcVWmT(gt1,gt2,g2,ZA,ZP,cplAhHpmcVWm(gt1,gt2))

 End Do 
End Do 


cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdAhL = 0._dp 
cplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdAhT(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,cplcFdFdAhL(gt1,gt2,gt3)              & 
& ,cplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


vev = Sqrt(vd**2 + vu**2 + vL(1)**2 + vL(2)**2 + vL(3)**2)
!----------------------------------------------------
! Pseudo Scalar Higgs coupling ratios 
!----------------------------------------------------
 
Do i2=1, 5
ratPCha(i2) = cplcChaChaAhL(i2,i2,i1)*1._dp*vev/MCha(i2) 
End Do 
Do i2=1, 3
ratPFd(i2) = cplcFdFdAhL(i2,i2,i1)*1._dp*vev/MFd(i2) 
End Do 
Do i2=1, 3
ratPFu(i2) = cplcFuFuAhL(i2,i2,i1)*1._dp*vev/MFu(i2) 
End Do 
Do i2=1, 8
ratPHpm(i2) = 0.5_dp*cplAhHpmcHpm(i1,i2,i2)*vev/MHpm2(i2) 
End Do 
Do i2=1, 6
ratPSd(i2) = 0.5_dp*cplAhSdcSd(i1,i2,i2)*vev/MSd2(i2) 
End Do 
Do i2=1, 6
ratPSu(i2) = 0.5_dp*cplAhSucSu(i1,i2,i2)*vev/MSu2(i2) 
End Do 
ratPVWm = 0._dp 
Call CoupPseudoHiggsToPhoton(m_in,i1,ratPCha,ratPFd,ratPFu,ratPHpm,ratPSd,            & 
& ratPSu,ratPVWm,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,coup)

cplPseudoHiggsPP(i1) = 2._dp*coup*Alpha 
CoupAPP(i1) = 2._dp*coup 
Call CoupPseudoHiggsToPhotonSM(m_in,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,coup)

ratioPPP(i1) = Abs(cplPseudoHiggsPP(i1)/(2._dp*coup*oo4pi*(1._dp-mW2/mZ2)*g2**2))**2 
Call CoupPseudoHiggsToGluon(Mhh(i1),i1,ratPFd,ratPFu,ratPSd,ratPSu,MFd,               & 
& MFu,MSd,MSu,coup)

cplPseudoHiggsGG(i1) = 2._dp*coup*AlphaS_MZ 
CoupAGG(i1) = 2._dp*coup 
Call CoupPseudoHiggsToGluonSM(m_in,MFd,MFu,MSd,MSu,coup)

ratioPGG(i1) = Abs(cplPseudoHiggsGG(i1)/(2._dp*coup*oo4pi*g3**2))**2 

!----------------------------------------------------
! Coupling ratios for HiggsBounds 
!----------------------------------------------------
 
Do i2=1, 5
rHB_P_S_Cha(i1,i2) = 1._dp*Abs((cplcChaChaAhL(i2,i2,i1)+cplcChaChaAhR(i2,i2,i1))*vev/(2._dp*MCha(i2)))**2 
rHB_P_P_Cha(i1,i2) = 1._dp*Abs((cplcChaChaAhL(i2,i2,i1)-cplcChaChaAhR(i2,i2,i1))*vev/(2._dp*MCha(i2)))**2 
End Do 
Do i2=1, 3
rHB_P_S_Fu(i1,i2) = 1._dp*Abs((cplcFuFuAhL(i2,i2,i1)+cplcFuFuAhR(i2,i2,i1))*vev/(2._dp*MFu(i2)))**2 
rHB_P_P_Fu(i1,i2) = 1._dp*Abs((cplcFuFuAhL(i2,i2,i1)-cplcFuFuAhR(i2,i2,i1))*vev/(2._dp*MFu(i2)))**2 
End Do 
Do i2=1, 3
rHB_P_S_Fd(i1,i2) = 1._dp*Abs((cplcFdFdAhL(i2,i2,i1)+cplcFdFdAhR(i2,i2,i1))*vev/(2._dp*MFd(i2)))**2 
rHB_P_P_Fd(i1,i2) = 1._dp*Abs((cplcFdFdAhL(i2,i2,i1)-cplcFdFdAhR(i2,i2,i1))*vev/(2._dp*MFd(i2)))**2 
End Do 
Do i2=1, 10
rHB_P_S_Chi(i1,i2) = 1._dp*Abs((cplChiChiAhL(i2,i2,i1)+cplChiChiAhR(i2,i2,i1))*vev/(2._dp*MChi(i2)))**2 
rHB_P_P_Chi(i1,i2) = 1._dp*Abs((cplChiChiAhL(i2,i2,i1)-cplChiChiAhR(i2,i2,i1))*vev/(2._dp*MChi(i2)))**2 
End Do 
If (i1.eq.2) Then 
CPL_A_A_Z = 0._dp 
End if 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhAhAh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhAhAhT(gt1,gt2,gt3,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZA,cplAhAhAh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhAhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhAhhhT(gt1,gt2,gt3,g1,g2,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,             & 
& ZH,ZA,cplAhAhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhhh = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhhhhhT(gt1,gt2,gt3,lam,Tlam,Yv,Tv,kap,Tk,vd,vu,vL,vR,ZH,ZA,             & 
& cplAhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhHpmcHpm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhHpmcHpmT(gt1,gt2,gt3,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,vd,vu,vL,             & 
& vR,ZA,ZP,cplAhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSdcSdT(gt1,gt2,gt3,Yd,Td,lam,vu,vR,ZD,ZA,cplAhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSucSuT(gt1,gt2,gt3,lam,Yv,Yu,Tu,vd,vL,vR,ZU,ZA,cplAhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhhhVZ = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhhhVZT(gt1,gt2,g1,g2,ZH,ZA,TW,cplAhhhVZ(gt1,gt2))

 End Do 
End Do 


cplAhHpmcVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhHpmcVWmT(gt1,gt2,g2,ZA,ZP,cplAhHpmcVWm(gt1,gt2))

 End Do 
End Do 


cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdAhL = 0._dp 
cplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdAhT(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,cplcFdFdAhL(gt1,gt2,gt3)              & 
& ,cplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Ah_decays_2B
 
Subroutine CouplingsFor_Hpm_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplAhHpmcHpm,cplAhcHpmVWm,cplChiChacHpmL,cplChiChacHpmR,cplcFuFdcHpmL,cplcFuFdcHpmR,   & 
& cplhhHpmcHpm,cplhhcHpmVWm,cplHpmcHpmVZ,cplSdcHpmcSu,cplcHpmVWmVZ,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplAhHpmcHpm(8,8,8),cplAhcHpmVWm(8,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),  & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplhhHpmcHpm(8,8,8),cplhhcHpmVWm(8,8),       & 
& cplHpmcHpmVZ(8,8),cplSdcHpmcSu(6,8,6),cplcHpmVWmVZ(8)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Hpm_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhHpmcHpm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplingAhHpmcHpmT(gt1,gt2,gt3,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,vd,vu,vL,             & 
& vR,ZA,ZP,cplAhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhHpmcHpm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
  Do gt3 = 1, 8
Call CouplinghhHpmcHpmT(gt1,gt2,gt3,g1,g2,Ye,Te,lam,Tlam,Yv,Tv,kap,vd,vu,             & 
& vL,vR,ZH,ZP,cplhhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcHpmcSu = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 8
  Do gt3 = 1, 6
Call CouplingSdcHpmcSuT(gt1,gt2,gt3,g2,Yd,Td,Ye,lam,Yv,Yu,Tu,vd,vu,vL,vR,             & 
& ZD,ZU,ZP,cplSdcHpmcSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplAhcHpmVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhcHpmVWmT(gt1,gt2,g2,ZA,ZP,cplAhcHpmVWm(gt1,gt2))

 End Do 
End Do 


cplhhcHpmVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplinghhcHpmVWmT(gt1,gt2,g2,ZH,ZP,cplhhcHpmVWm(gt1,gt2))

 End Do 
End Do 


cplHpmcHpmVZ = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingHpmcHpmVZT(gt1,gt2,g1,g2,ZP,TW,cplHpmcHpmVZ(gt1,gt2))

 End Do 
End Do 


cplcHpmVWmVZ = 0._dp 
Do gt1 = 1, 8
Call CouplingcHpmVWmVZT(gt1,g1,g2,vd,vu,vL,ZP,TW,cplcHpmVWmVZ(gt1))

End Do 


cplChiChacHpmL = 0._dp 
cplChiChacHpmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingChiChacHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplChiChacHpmL(gt1,gt2,gt3)& 
& ,cplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFdcHpmL = 0._dp 
cplcFuFdcHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFdcHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFuFdcHpmL(gt1,gt2,gt3)& 
& ,cplcFuFdcHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Hpm_decays_2B
 
Subroutine CouplingsFor_Glu_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Glu_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


cplGluFucSuL = 0._dp 
cplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFucSuT(gt2,gt3,g3,pG,ZU,ZUL,ZUR,cplGluFucSuL(gt2,gt3),cplGluFucSuR(gt2,gt3))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Glu_decays_2B
 
Subroutine CouplingsFor_Fu_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplcFuFuAhL,cplcFuFuAhR,cplcFuChiSuL,cplcFuChiSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,         & 
& cplcFuFdcVWmL,cplcFuFdcVWmR,cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,           & 
& cplcFuGluSuL,cplcFuGluSuR,cplcChacFuSdL,cplcChacFuSdR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),      & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplcFuFdcVWmL(3,3),cplcFuFdcVWmR(3,3),       & 
& cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),               & 
& cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Fu_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuChiSuL = 0._dp 
cplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFuChiSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplcFuChiSuL(gt1,gt2,gt3)   & 
& ,cplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFdcHpmL = 0._dp 
cplcFuFdcHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFdcHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFuFdcHpmL(gt1,gt2,gt3)& 
& ,cplcFuFdcHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuGluSuL = 0._dp 
cplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFuGluSuT(gt1,gt3,g3,pG,ZU,ZUL,ZUR,cplcFuGluSuL(gt1,gt3),cplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


cplcChacFuSdL = 0._dp 
cplcChacFuSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChacFuSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplcChacFuSdL(gt1,gt2,gt3)& 
& ,cplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFdcVWmL = 0._dp 
cplcFuFdcVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFdcVWmT(gt1,gt2,g2,ZDL,ZUL,cplcFuFdcVWmL(gt1,gt2),cplcFuFdcVWmR(gt1,gt2))

 End Do 
End Do 


cplcFuFuVZL = 0._dp 
cplcFuFuVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFuVZT(gt1,gt2,g1,g2,TW,cplcFuFuVZL(gt1,gt2),cplcFuFuVZR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Fu_decays_2B
 
Subroutine CouplingsFor_Cha_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,   & 
& cplcChaChiHpmL,cplcChaChiHpmR,cplcChaChiVWmL,cplcChaChiVWmR,cplcChaFdcSuL,             & 
& cplcChaFdcSuR,cplcChacFuSdL,cplcChacFuSdR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),  & 
& cplcChaChaVZL(5,5),cplcChaChaVZR(5,5),cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),   & 
& cplcChaChiVWmL(5,10),cplcChaChiVWmR(5,10),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),   & 
& cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Cha_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChiHpmL = 0._dp 
cplcChaChiHpmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingcChaChiHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplcChaChiHpmL(gt1,gt2,gt3)& 
& ,cplcChaChiHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaFdcSuL = 0._dp 
cplcChaFdcSuR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChaFdcSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcChaFdcSuL(gt1,gt2,gt3)& 
& ,cplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChacFuSdL = 0._dp 
cplcChacFuSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChacFuSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplcChacFuSdL(gt1,gt2,gt3)& 
& ,cplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChaVZL = 0._dp 
cplcChaChaVZR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
Call CouplingcChaChaVZT(gt1,gt2,g1,g2,ZER,ZEL,TW,cplcChaChaVZL(gt1,gt2)               & 
& ,cplcChaChaVZR(gt1,gt2))

 End Do 
End Do 


cplcChaChiVWmL = 0._dp 
cplcChaChiVWmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
Call CouplingcChaChiVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplcChaChiVWmL(gt1,gt2),               & 
& cplcChaChiVWmR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Cha_decays_2B
 
Subroutine CouplingsFor_Chi_decays_2B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplChiChiAhL,cplChiChiAhR,cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,cplChiChacVWmR, & 
& cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,cplChiFdcSdL,cplChiFdcSdR,         & 
& cplChiFucSuL,cplChiFucSuR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),& 
& cplChiChacVWmL(10,5),cplChiChacVWmR(10,5),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8), & 
& cplChiChiVZL(10,10),cplChiChiVZR(10,10),cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),     & 
& cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Chi_2B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacHpmL = 0._dp 
cplChiChacHpmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingChiChacHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplChiChacHpmL(gt1,gt2,gt3)& 
& ,cplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFucSuL = 0._dp 
cplChiFucSuR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFucSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplChiFucSuL(gt1,gt2,gt3)   & 
& ,cplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacVWmL = 0._dp 
cplChiChacVWmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
Call CouplingChiChacVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplChiChacVWmL(gt1,gt2),               & 
& cplChiChacVWmR(gt1,gt2))

 End Do 
End Do 


cplChiChiVZL = 0._dp 
cplChiChiVZR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
Call CouplingChiChiVZT(gt1,gt2,g1,g2,UV,TW,cplChiChiVZL(gt1,gt2),cplChiChiVZR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Chi_decays_2B
 
Subroutine CouplingsFor_Glu_decays_3B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplcChacFuSdL,cplcChacFuSdR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFdChiSdL,cplcFdChiSdR,     & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuChiSuL,cplcFuChiSuR,cplcFuGluSuL,cplcFuGluSuR,         & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplGluFucSuL,cplGluFucSuR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),  & 
& cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),         & 
& cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),         & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Glu_3B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFucSuL = 0._dp 
cplChiFucSuR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFucSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplChiFucSuL(gt1,gt2,gt3)   & 
& ,cplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChiSdL = 0._dp 
cplcFdChiSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFdChiSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplcFdChiSdL(gt1,gt2,gt3)   & 
& ,cplcFdChiSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuChiSuL = 0._dp 
cplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFuChiSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplcFuChiSuL(gt1,gt2,gt3)   & 
& ,cplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


cplcChaFdcSuL = 0._dp 
cplcChaFdcSuR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChaFdcSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcChaFdcSuL(gt1,gt2,gt3)& 
& ,cplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFucSuL = 0._dp 
cplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFucSuT(gt2,gt3,g3,pG,ZU,ZUL,ZUR,cplGluFucSuL(gt2,gt3),cplGluFucSuR(gt2,gt3))

 End Do 
End Do 


cplcFdGluSdL = 0._dp 
cplcFdGluSdR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFdGluSdT(gt1,gt3,g3,pG,ZD,ZDL,ZDR,cplcFdGluSdL(gt1,gt3),cplcFdGluSdR(gt1,gt3))

 End Do 
End Do 


cplcFuGluSuL = 0._dp 
cplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFuGluSuT(gt1,gt3,g3,pG,ZU,ZUL,ZUR,cplcFuGluSuL(gt1,gt3),cplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


cplcChacFuSdL = 0._dp 
cplcChacFuSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChacFuSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplcChacFuSdL(gt1,gt2,gt3)& 
& ,cplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Glu_decays_3B
 
Subroutine CouplingsFor_Sd_decays_3B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplAhSdcSd,cplcChacFuSdL,cplcChacFuSdR,cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,      & 
& cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,cplcChaChiHpmR,               & 
& cplcChaFdcSuL,cplcChaFdcSuR,cplcFdChaSuL,cplcFdChaSuR,cplcFdChiSdL,cplcFdChiSdR,       & 
& cplcFdFdAhL,cplcFdFdAhR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,cplcFdFdVZR,               & 
& cplcFdFuHpmL,cplcFdFuHpmR,cplcFdGluSdL,cplcFdGluSdR,cplcFuChiSuL,cplcFuChiSuR,         & 
& cplcFuFdcHpmL,cplcFuFdcHpmR,cplcFuFdcVWmL,cplcFuFdcVWmR,cplcFuFuAhL,cplcFuFuAhR,       & 
& cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,cplcFuGluSuL,cplcFuGluSuR,             & 
& cplChaFucSdL,cplChaFucSdR,cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,cplChiChacVWmR, & 
& cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,         & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplGluFucSuL,cplGluFucSuR,cplhhSdcSd,cplHpmSucSd,cplSdcSdVZ,cplSucSdVWm,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplAhSdcSd(8,6,6),cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaChaAhL(5,5,8),     & 
& cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),cplcChaChaVZL(5,5),     & 
& cplcChaChaVZR(5,5),cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcChaFdcSuL(5,3,6), & 
& cplcChaFdcSuR(5,3,6),cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplcFdChiSdL(3,10,6),     & 
& cplcFdChiSdR(3,10,6),cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),cplcFdFdhhL(3,3,8),         & 
& cplcFdFdhhR(3,3,8),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFdFuHpmL(3,3,8),              & 
& cplcFdFuHpmR(3,3,8),cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuChiSuL(3,10,6),          & 
& cplcFuChiSuR(3,10,6),cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplcFuFdcVWmL(3,3),     & 
& cplcFuFdcVWmR(3,3),cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplcFuFuhhL(3,3,8),           & 
& cplcFuFuhhR(3,3,8),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),& 
& cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8), & 
& cplChiChacVWmL(10,5),cplChiChacVWmR(10,5),cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8), & 
& cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),cplChiChiVZL(10,10),cplChiChiVZR(10,10),   & 
& cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),   & 
& cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6),               & 
& cplhhSdcSd(8,6,6),cplHpmSucSd(8,6,6),cplSdcSdVZ(6,6),cplSucSdVWm(6,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Sd_3B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSdcSdT(gt1,gt2,gt3,Yd,Td,lam,vu,vR,ZD,ZA,cplAhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSdcSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSdcSdT(gt1,gt2,gt3,g1,g2,Yd,Td,lam,vd,vu,vL,vR,ZD,ZH,cplhhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplHpmSucSd = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingHpmSucSdT(gt1,gt2,gt3,g2,Yd,Td,Ye,lam,Yv,Yu,Tu,vd,vu,vL,vR,              & 
& ZD,ZU,ZP,cplHpmSucSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcSdVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSdcSdVZT(gt1,gt2,g1,g2,ZD,TW,cplSdcSdVZ(gt1,gt2))

 End Do 
End Do 


cplSucSdVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSucSdVWmT(gt1,gt2,g2,ZD,ZU,cplSucSdVWm(gt1,gt2))

 End Do 
End Do 


cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdAhL = 0._dp 
cplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdAhT(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,cplcFdFdAhL(gt1,gt2,gt3)              & 
& ,cplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacHpmL = 0._dp 
cplChiChacHpmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingChiChacHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplChiChacHpmL(gt1,gt2,gt3)& 
& ,cplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChaFucSdL = 0._dp 
cplChaFucSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChaFucSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplChaFucSdL(gt1,gt2,gt3)& 
& ,cplChaFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChaSuL = 0._dp 
cplcFdChaSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 5
  Do gt3 = 1, 6
Call CouplingcFdChaSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcFdChaSuL(gt1,gt2,gt3)& 
& ,cplcFdChaSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFucSuL = 0._dp 
cplChiFucSuR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFucSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplChiFucSuL(gt1,gt2,gt3)   & 
& ,cplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChiHpmL = 0._dp 
cplcChaChiHpmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingcChaChiHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplcChaChiHpmL(gt1,gt2,gt3)& 
& ,cplcChaChiHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChiSdL = 0._dp 
cplcFdChiSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFdChiSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplcFdChiSdL(gt1,gt2,gt3)   & 
& ,cplcFdChiSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuChiSuL = 0._dp 
cplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFuChiSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplcFuChiSuL(gt1,gt2,gt3)   & 
& ,cplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


cplcFdFdhhL = 0._dp 
cplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdhhT(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,cplcFdFdhhL(gt1,gt2,gt3)              & 
& ,cplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaFdcSuL = 0._dp 
cplcChaFdcSuR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChaFdcSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcChaFdcSuL(gt1,gt2,gt3)& 
& ,cplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFdcHpmL = 0._dp 
cplcFuFdcHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFdcHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFuFdcHpmL(gt1,gt2,gt3)& 
& ,cplcFuFdcHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFucSuL = 0._dp 
cplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFucSuT(gt2,gt3,g3,pG,ZU,ZUL,ZUR,cplGluFucSuL(gt2,gt3),cplGluFucSuR(gt2,gt3))

 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFuHpmL = 0._dp 
cplcFdFuHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFuHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFdFuHpmL(gt1,gt2,gt3) & 
& ,cplcFdFuHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdGluSdL = 0._dp 
cplcFdGluSdR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFdGluSdT(gt1,gt3,g3,pG,ZD,ZDL,ZDR,cplcFdGluSdL(gt1,gt3),cplcFdGluSdR(gt1,gt3))

 End Do 
End Do 


cplcFuGluSuL = 0._dp 
cplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFuGluSuT(gt1,gt3,g3,pG,ZU,ZUL,ZUR,cplcFuGluSuL(gt1,gt3),cplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


cplcChacFuSdL = 0._dp 
cplcChacFuSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChacFuSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplcChacFuSdL(gt1,gt2,gt3)& 
& ,cplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacVWmL = 0._dp 
cplChiChacVWmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
Call CouplingChiChacVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplChiChacVWmL(gt1,gt2),               & 
& cplChiChacVWmR(gt1,gt2))

 End Do 
End Do 


cplcChaChaVZL = 0._dp 
cplcChaChaVZR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
Call CouplingcChaChaVZT(gt1,gt2,g1,g2,ZER,ZEL,TW,cplcChaChaVZL(gt1,gt2)               & 
& ,cplcChaChaVZR(gt1,gt2))

 End Do 
End Do 


cplChiChiVZL = 0._dp 
cplChiChiVZR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
Call CouplingChiChiVZT(gt1,gt2,g1,g2,UV,TW,cplChiChiVZL(gt1,gt2),cplChiChiVZR(gt1,gt2))

 End Do 
End Do 


cplcFdFdVZL = 0._dp 
cplcFdFdVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFdVZT(gt1,gt2,g1,g2,TW,cplcFdFdVZL(gt1,gt2),cplcFdFdVZR(gt1,gt2))

 End Do 
End Do 


cplcFuFdcVWmL = 0._dp 
cplcFuFdcVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFdcVWmT(gt1,gt2,g2,ZDL,ZUL,cplcFuFdcVWmL(gt1,gt2),cplcFuFdcVWmR(gt1,gt2))

 End Do 
End Do 


cplcFuFuVZL = 0._dp 
cplcFuFuVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFuVZT(gt1,gt2,g1,g2,TW,cplcFuFuVZL(gt1,gt2),cplcFuFuVZR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Sd_decays_3B
 
Subroutine CouplingsFor_Su_decays_3B(m_in,i1,MAhinput,MAh2input,MChainput,            & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplAhSucSu,cplcChacFuSdL,cplcChacFuSdR,cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,      & 
& cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,cplcChaChiHpmR,               & 
& cplcChaChiVWmL,cplcChaChiVWmR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFdChaSuL,cplcFdChaSuR,   & 
& cplcFdChiSdL,cplcFdChiSdR,cplcFdFdAhL,cplcFdFdAhR,cplcFdFdhhL,cplcFdFdhhR,             & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdFuVWmL,cplcFdFuVWmR,           & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuChiSuL,cplcFuChiSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplcFuFuAhL,cplcFuFuAhR,cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplcFuGluSuL,cplcFuGluSuR,cplChaFucSdL,cplChaFucSdR,cplChiChacHpmL,cplChiChacHpmR,     & 
& cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,         & 
& cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplGluFdcSdL,cplGluFdcSdR,         & 
& cplGluFucSuL,cplGluFucSuR,cplhhSucSu,cplSdcHpmcSu,cplSdcSucVWm,cplSucSuVZ,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplAhSucSu(8,6,6),cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaChaAhL(5,5,8),     & 
& cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),cplcChaChaVZL(5,5),     & 
& cplcChaChaVZR(5,5),cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcChaChiVWmL(5,10), & 
& cplcChaChiVWmR(5,10),cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),cplcFdChaSuL(3,5,6),    & 
& cplcFdChaSuR(3,5,6),cplcFdChiSdL(3,10,6),cplcFdChiSdR(3,10,6),cplcFdFdAhL(3,3,8),      & 
& cplcFdFdAhR(3,3,8),cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),cplcFdFdVZL(3,3),             & 
& cplcFdFdVZR(3,3),cplcFdFuHpmL(3,3,8),cplcFdFuHpmR(3,3,8),cplcFdFuVWmL(3,3),            & 
& cplcFdFuVWmR(3,3),cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuChiSuL(3,10,6),            & 
& cplcFuChiSuR(3,10,6),cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplcFuFuAhL(3,3,8),     & 
& cplcFuFuAhR(3,3,8),cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),cplcFuFuVZL(3,3),             & 
& cplcFuFuVZR(3,3),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),cplChaFucSdL(5,3,6),              & 
& cplChaFucSdR(5,3,6),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),cplChiChiAhL(10,10,8),& 
& cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),cplChiChiVZL(10,10), & 
& cplChiChiVZR(10,10),cplChiFdcSdL(10,3,6),cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),    & 
& cplChiFucSuR(10,3,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),            & 
& cplGluFucSuR(3,6),cplhhSucSu(8,6,6),cplSdcHpmcSu(6,8,6),cplSdcSucVWm(6,6),             & 
& cplSucSuVZ(6,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Su_3B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplAhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplingAhSucSuT(gt1,gt2,gt3,lam,Yv,Yu,Tu,vd,vL,vR,ZU,ZA,cplAhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplhhSucSu = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CouplinghhSucSuT(gt1,gt2,gt3,g1,g2,lam,Yv,Yu,Tu,vd,vu,vL,vR,ZU,ZH,               & 
& cplhhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcHpmcSu = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 8
  Do gt3 = 1, 6
Call CouplingSdcHpmcSuT(gt1,gt2,gt3,g2,Yd,Td,Ye,lam,Yv,Yu,Tu,vd,vu,vL,vR,             & 
& ZD,ZU,ZP,cplSdcHpmcSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplSdcSucVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSdcSucVWmT(gt1,gt2,g2,ZD,ZU,cplSdcSucVWm(gt1,gt2))

 End Do 
End Do 


cplSucSuVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSucSuVZT(gt1,gt2,g1,g2,ZU,TW,cplSucSuVZ(gt1,gt2))

 End Do 
End Do 


cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdAhL = 0._dp 
cplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdAhT(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,cplcFdFdAhL(gt1,gt2,gt3)              & 
& ,cplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacHpmL = 0._dp 
cplChiChacHpmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingChiChacHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplChiChacHpmL(gt1,gt2,gt3)& 
& ,cplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChaFucSdL = 0._dp 
cplChaFucSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChaFucSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplChaFucSdL(gt1,gt2,gt3)& 
& ,cplChaFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChaSuL = 0._dp 
cplcFdChaSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 5
  Do gt3 = 1, 6
Call CouplingcFdChaSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcFdChaSuL(gt1,gt2,gt3)& 
& ,cplcFdChaSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFucSuL = 0._dp 
cplChiFucSuR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFucSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplChiFucSuL(gt1,gt2,gt3)   & 
& ,cplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChiHpmL = 0._dp 
cplcChaChiHpmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingcChaChiHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplcChaChiHpmL(gt1,gt2,gt3)& 
& ,cplcChaChiHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChiSdL = 0._dp 
cplcFdChiSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFdChiSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplcFdChiSdL(gt1,gt2,gt3)   & 
& ,cplcFdChiSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuChiSuL = 0._dp 
cplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFuChiSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplcFuChiSuL(gt1,gt2,gt3)   & 
& ,cplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


cplcFdFdhhL = 0._dp 
cplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdhhT(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,cplcFdFdhhL(gt1,gt2,gt3)              & 
& ,cplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaFdcSuL = 0._dp 
cplcChaFdcSuR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChaFdcSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcChaFdcSuL(gt1,gt2,gt3)& 
& ,cplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFdcHpmL = 0._dp 
cplcFuFdcHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFdcHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFuFdcHpmL(gt1,gt2,gt3)& 
& ,cplcFuFdcHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFucSuL = 0._dp 
cplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFucSuT(gt2,gt3,g3,pG,ZU,ZUL,ZUR,cplGluFucSuL(gt2,gt3),cplGluFucSuR(gt2,gt3))

 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFuHpmL = 0._dp 
cplcFdFuHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFuHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFdFuHpmL(gt1,gt2,gt3) & 
& ,cplcFdFuHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdGluSdL = 0._dp 
cplcFdGluSdR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFdGluSdT(gt1,gt3,g3,pG,ZD,ZDL,ZDR,cplcFdGluSdL(gt1,gt3),cplcFdGluSdR(gt1,gt3))

 End Do 
End Do 


cplcFuGluSuL = 0._dp 
cplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFuGluSuT(gt1,gt3,g3,pG,ZU,ZUL,ZUR,cplcFuGluSuL(gt1,gt3),cplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


cplcChacFuSdL = 0._dp 
cplcChacFuSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChacFuSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplcChacFuSdL(gt1,gt2,gt3)& 
& ,cplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChaVZL = 0._dp 
cplcChaChaVZR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
Call CouplingcChaChaVZT(gt1,gt2,g1,g2,ZER,ZEL,TW,cplcChaChaVZL(gt1,gt2)               & 
& ,cplcChaChaVZR(gt1,gt2))

 End Do 
End Do 


cplChiChiVZL = 0._dp 
cplChiChiVZR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
Call CouplingChiChiVZT(gt1,gt2,g1,g2,UV,TW,cplChiChiVZL(gt1,gt2),cplChiChiVZR(gt1,gt2))

 End Do 
End Do 


cplcChaChiVWmL = 0._dp 
cplcChaChiVWmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
Call CouplingcChaChiVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplcChaChiVWmL(gt1,gt2),               & 
& cplcChaChiVWmR(gt1,gt2))

 End Do 
End Do 


cplcFdFdVZL = 0._dp 
cplcFdFdVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFdVZT(gt1,gt2,g1,g2,TW,cplcFdFdVZL(gt1,gt2),cplcFdFdVZR(gt1,gt2))

 End Do 
End Do 


cplcFdFuVWmL = 0._dp 
cplcFdFuVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFuVWmT(gt1,gt2,g2,ZDL,ZUL,cplcFdFuVWmL(gt1,gt2),cplcFdFuVWmR(gt1,gt2))

 End Do 
End Do 


cplcFuFuVZL = 0._dp 
cplcFuFuVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFuVZT(gt1,gt2,g1,g2,TW,cplcFuFuVZL(gt1,gt2),cplcFuFuVZR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Su_decays_3B
 
Subroutine CouplingsFor_Cha_decays_3B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplcChacFuSdL,cplcChacFuSdR,cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,cplcChaChahhR,   & 
& cplcChaChaVZL,cplcChaChaVZR,cplcChaChiHpmL,cplcChaChiHpmR,cplcChaChiVWmL,              & 
& cplcChaChiVWmR,cplcChaFdcSuL,cplcChaFdcSuR,cplcFdChaSuL,cplcFdChaSuR,cplcFdFdAhL,      & 
& cplcFdFdAhR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,cplcFdFdVZR,cplcFuChiSuL,              & 
& cplcFuChiSuR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplcFuFdcVWmL,cplcFuFdcVWmR,cplcFuFuAhL,      & 
& cplcFuFuAhR,cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,cplcFuGluSuL,              & 
& cplcFuGluSuR,cplChaFucSdL,cplChaFucSdR,cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,   & 
& cplChiChacVWmR,cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,       & 
& cplChiChiVZR,cplChiFdcSdL,cplChiFdcSdR,cplGluFdcSdL,cplGluFdcSdR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),  & 
& cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),cplcChaChaVZL(5,5),cplcChaChaVZR(5,5),       & 
& cplcChaChiHpmL(5,10,8),cplcChaChiHpmR(5,10,8),cplcChaChiVWmL(5,10),cplcChaChiVWmR(5,10),& 
& cplcChaFdcSuL(5,3,6),cplcChaFdcSuR(5,3,6),cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),     & 
& cplcFdFdAhL(3,3,8),cplcFdFdAhR(3,3,8),cplcFdFdhhL(3,3,8),cplcFdFdhhR(3,3,8),           & 
& cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFuChiSuL(3,10,6),cplcFuChiSuR(3,10,6),           & 
& cplcFuFdcHpmL(3,3,8),cplcFuFdcHpmR(3,3,8),cplcFuFdcVWmL(3,3),cplcFuFdcVWmR(3,3),       & 
& cplcFuFuAhL(3,3,8),cplcFuFuAhR(3,3,8),cplcFuFuhhL(3,3,8),cplcFuFuhhR(3,3,8),           & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcFuGluSuL(3,6),cplcFuGluSuR(3,6),cplChaFucSdL(5,3,6),& 
& cplChaFucSdR(5,3,6),cplChiChacHpmL(10,5,8),cplChiChacHpmR(10,5,8),cplChiChacVWmL(10,5),& 
& cplChiChacVWmR(10,5),cplChiChiAhL(10,10,8),cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),& 
& cplChiChihhR(10,10,8),cplChiChiVZL(10,10),cplChiChiVZR(10,10),cplChiFdcSdL(10,3,6),    & 
& cplChiFdcSdR(10,3,6),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6)

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Cha_3B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdAhL = 0._dp 
cplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdAhT(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,cplcFdFdAhL(gt1,gt2,gt3)              & 
& ,cplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacHpmL = 0._dp 
cplChiChacHpmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingChiChacHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplChiChacHpmL(gt1,gt2,gt3)& 
& ,cplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChaFucSdL = 0._dp 
cplChaFucSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChaFucSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplChaFucSdL(gt1,gt2,gt3)& 
& ,cplChaFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChaSuL = 0._dp 
cplcFdChaSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 5
  Do gt3 = 1, 6
Call CouplingcFdChaSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcFdChaSuL(gt1,gt2,gt3)& 
& ,cplcFdChaSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChiHpmL = 0._dp 
cplcChaChiHpmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingcChaChiHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplcChaChiHpmL(gt1,gt2,gt3)& 
& ,cplcChaChiHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuChiSuL = 0._dp 
cplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFuChiSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplcFuChiSuL(gt1,gt2,gt3)   & 
& ,cplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


cplcFdFdhhL = 0._dp 
cplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdhhT(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,cplcFdFdhhL(gt1,gt2,gt3)              & 
& ,cplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaFdcSuL = 0._dp 
cplcChaFdcSuR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChaFdcSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcChaFdcSuL(gt1,gt2,gt3)& 
& ,cplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFdcHpmL = 0._dp 
cplcFuFdcHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFdcHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFuFdcHpmL(gt1,gt2,gt3)& 
& ,cplcFuFdcHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuGluSuL = 0._dp 
cplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFuGluSuT(gt1,gt3,g3,pG,ZU,ZUL,ZUR,cplcFuGluSuL(gt1,gt3),cplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


cplcChacFuSdL = 0._dp 
cplcChacFuSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingcChacFuSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplcChacFuSdL(gt1,gt2,gt3)& 
& ,cplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacVWmL = 0._dp 
cplChiChacVWmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
Call CouplingChiChacVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplChiChacVWmL(gt1,gt2),               & 
& cplChiChacVWmR(gt1,gt2))

 End Do 
End Do 


cplcChaChaVZL = 0._dp 
cplcChaChaVZR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
Call CouplingcChaChaVZT(gt1,gt2,g1,g2,ZER,ZEL,TW,cplcChaChaVZL(gt1,gt2)               & 
& ,cplcChaChaVZR(gt1,gt2))

 End Do 
End Do 


cplChiChiVZL = 0._dp 
cplChiChiVZR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
Call CouplingChiChiVZT(gt1,gt2,g1,g2,UV,TW,cplChiChiVZL(gt1,gt2),cplChiChiVZR(gt1,gt2))

 End Do 
End Do 


cplcChaChiVWmL = 0._dp 
cplcChaChiVWmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
Call CouplingcChaChiVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplcChaChiVWmL(gt1,gt2),               & 
& cplcChaChiVWmR(gt1,gt2))

 End Do 
End Do 


cplcFdFdVZL = 0._dp 
cplcFdFdVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFdVZT(gt1,gt2,g1,g2,TW,cplcFdFdVZL(gt1,gt2),cplcFdFdVZR(gt1,gt2))

 End Do 
End Do 


cplcFuFdcVWmL = 0._dp 
cplcFuFdcVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFdcVWmT(gt1,gt2,g2,ZDL,ZUL,cplcFuFdcVWmL(gt1,gt2),cplcFuFdcVWmR(gt1,gt2))

 End Do 
End Do 


cplcFuFuVZL = 0._dp 
cplcFuFuVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFuVZT(gt1,gt2,g1,g2,TW,cplcFuFuVZL(gt1,gt2),cplcFuFuVZR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Cha_decays_3B
 
Subroutine CouplingsFor_Chi_decays_3B(m_in,i1,MAhinput,MAh2input,MChainput,           & 
& MCha2input,MChiinput,MChi2input,MFdinput,MFd2input,MFuinput,MFu2input,MGluinput,       & 
& MGlu2input,Mhhinput,Mhh2input,MHpminput,MHpm2input,MSdinput,MSd2input,MSuinput,        & 
& MSu2input,MVWminput,MVWm2input,MVZinput,MVZ2input,pGinput,TWinput,ZERinput,            & 
& ZELinput,ZAinput,ZDinput,ZDLinput,ZDRinput,ZHinput,UVinput,ZPinput,ZUinput,            & 
& ZULinput,ZURinput,ZWinput,ZZinput,g1input,g2input,g3input,Ydinput,Yeinput,             & 
& laminput,Yvinput,Yuinput,kapinput,Tdinput,Teinput,Tlaminput,Tvinput,Tuinput,           & 
& Tkinput,mq2input,ml2input,mHd2input,mHu2input,md2input,mu2input,me2input,              & 
& mv2input,mlHd2input,M1input,M2input,M3input,vdinput,vuinput,vLinput,vRinput,           & 
& cplcChaChaAhL,cplcChaChaAhR,cplcChaChahhL,cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,   & 
& cplcChaChiHpmL,cplcChaChiHpmR,cplcChaChiVWmL,cplcChaChiVWmR,cplcFdChaSuL,              & 
& cplcFdChaSuR,cplcFdChiSdL,cplcFdChiSdR,cplcFdFdAhL,cplcFdFdAhR,cplcFdFdhhL,            & 
& cplcFdFdhhR,cplcFdFdVZL,cplcFdFdVZR,cplcFdFuHpmL,cplcFdFuHpmR,cplcFdFuVWmL,            & 
& cplcFdFuVWmR,cplcFdGluSdL,cplcFdGluSdR,cplcFuChiSuL,cplcFuChiSuR,cplcFuFuAhL,          & 
& cplcFuFuAhR,cplcFuFuhhL,cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,cplcFuGluSuL,              & 
& cplcFuGluSuR,cplChaFucSdL,cplChaFucSdR,cplChiChacHpmL,cplChiChacHpmR,cplChiChacVWmL,   & 
& cplChiChacVWmR,cplChiChiAhL,cplChiChiAhR,cplChiChihhL,cplChiChihhR,cplChiChiVZL,       & 
& cplChiChiVZR,cplChiFdcSdL,cplChiFdcSdR,cplChiFucSuL,cplChiFucSuR,cplGluFdcSdL,         & 
& cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,deltaM)

Implicit None 

Real(dp), Intent(in) :: m_in 
Real(dp), Intent(in) :: deltaM 
Integer, Intent(in) :: i1 
Real(dp),Intent(in) :: g1input,g2input,g3input,mHd2input,mHu2input,mlHd2input(3),vdinput,vuinput,            & 
& vLinput(3),vRinput(3)

Complex(dp),Intent(in) :: Ydinput(3,3),Yeinput(3,3),laminput(3),Yvinput(3,3),Yuinput(3,3),kapinput(3,3,3),      & 
& Tdinput(3,3),Teinput(3,3),Tlaminput(3),Tvinput(3,3),Tuinput(3,3),Tkinput(3,3,3),       & 
& mq2input(3,3),ml2input(3,3),md2input(3,3),mu2input(3,3),me2input(3,3),mv2input(3,3),   & 
& M1input,M2input,M3input

Real(dp),Intent(in) :: MAhinput(8),MAh2input(8),MChainput(5),MCha2input(5),MChiinput(10),MChi2input(10),     & 
& MFdinput(3),MFd2input(3),MFuinput(3),MFu2input(3),MGluinput,MGlu2input,Mhhinput(8),    & 
& Mhh2input(8),MHpminput(8),MHpm2input(8),MSdinput(6),MSd2input(6),MSuinput(6),          & 
& MSu2input(6),MVWminput,MVWm2input,MVZinput,MVZ2input,TWinput,ZAinput(8,8),             & 
& ZHinput(8,8),ZPinput(8,8),ZZinput(2,2)

Complex(dp),Intent(in) :: pGinput,ZERinput(5,5),ZELinput(5,5),ZDinput(6,6),ZDLinput(3,3),ZDRinput(3,3),         & 
& UVinput(10,10),ZUinput(6,6),ZULinput(3,3),ZURinput(3,3),ZWinput(2,2)

Real(dp) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Complex(dp),Intent(out) :: cplcChaChaAhL(5,5,8),cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),  & 
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

Real(dp) ::  g1D(394) 
Integer :: i2, i3, gt1, gt2, gt3, kont 
Real(dp) :: MAh(8),MAh2(8),MCha(5),MCha2(5),MChi(10),MChi2(10),MFd(3),MFd2(3),MFu(3),             & 
& MFu2(3),MGlu,MGlu2,Mhh(8),Mhh2(8),MHpm(8),MHpm2(8),MSd(6),MSd2(6),MSu(6),              & 
& MSu2(6),MVWm,MVWm2,MVZ,MVZ2,TW,ZA(8,8),ZH(8,8),ZP(8,8),ZZ(2,2)

Complex(dp) :: pG,ZER(5,5),ZEL(5,5),ZD(6,6),ZDL(3,3),ZDR(3,3),UV(10,10),ZU(6,6),ZUL(3,3),            & 
& ZUR(3,3),ZW(2,2)

Real(dp) :: gSM(11), sinW2, dt, tz, Qin 
Iname = Iname + 1 
NameOfUnit(Iname) = 'Couplings_Chi_3B'
 
sinW2=1._dp-mW2/mZ2 
g1 = g1input 
g2 = g2input 
g3 = g3input 
Yd = Ydinput 
Ye = Yeinput 
lam = laminput 
Yv = Yvinput 
Yu = Yuinput 
kap = kapinput 
Td = Tdinput 
Te = Teinput 
Tlam = Tlaminput 
Tv = Tvinput 
Tu = Tuinput 
Tk = Tkinput 
mq2 = mq2input 
ml2 = ml2input 
mHd2 = mHd2input 
mHu2 = mHu2input 
md2 = md2input 
mu2 = mu2input 
me2 = me2input 
mv2 = mv2input 
mlHd2 = mlHd2input 
M1 = M1input 
M2 = M2input 
M3 = M3input 
vd = vdinput 
vu = vuinput 
vL = vLinput 
vR = vRinput 
Qin=sqrt(getRenormalizationScale()) 

 
 ! --- GUT normalize gauge couplings --- 
g1 = Sqrt(5._dp/3._dp)*g1 
! ----------------------- 
 
Call ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,               & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g1D)

If ((m_in.le.Qin).and.(RunningCouplingsDecays)) Then 
  tz=Log(m_in/Qin) 
  If (m_in.le.mz) tz=Log(mz/Qin)  
  dt=tz/50._dp 
  Call odeint(g1D,394,0._dp,tz,deltaM,dt,0._dp,rge394,kont)

End if 
Call GToParameters394(g1D,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,              & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)


 
 ! --- Remove GUT-normalization of gauge couplings --- 
g1 = Sqrt(3._dp/5._dp)*g1 
! ----------------------- 
 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,             & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,(/ ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC, ZeroC /))

! --- Calculate running tree-level masses for loop induced couplings and Quark mixing matrices --- 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,          & 
& Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,pG,TW,ZER,ZEL,               & 
& ZA,ZD,ZDL,ZDR,ZH,UV,ZP,ZU,ZUL,ZUR,ZW,ZZ,vd,vu,vR,vL,g1,g2,g3,Yd,Ye,lam,Yv,             & 
& Yu,kap,Td,Te,Tlam,Tv,Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,              & 
& M3,.True.,kont)

ZH = ZHinput 
ZA = ZAinput 
! --- Use the 1-loop mixing matrices calculated at M_SUSY in the vertices --- 
ZA = ZAinput 
ZD = ZDinput 
ZH = ZHinput 
UV = UVinput 
ZP = ZPinput 
ZU = ZUinput 
ZW = ZWinput 
ZZ = ZZinput 
cplcChaChaAhL = 0._dp 
cplcChaChaAhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChaAhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZA,ZER,ZEL,cplcChaChaAhL(gt1,gt2,gt3)& 
& ,cplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChiAhL = 0._dp 
cplChiChiAhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChiAhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZA,UV,cplChiChiAhL(gt1,gt2,gt3)   & 
& ,cplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFdAhL = 0._dp 
cplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdAhT(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,cplcFdFdAhL(gt1,gt2,gt3)              & 
& ,cplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuFuAhL = 0._dp 
cplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuAhT(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,cplcFuFuAhL(gt1,gt2,gt3)              & 
& ,cplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChacHpmL = 0._dp 
cplChiChacHpmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingChiChacHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplChiChacHpmL(gt1,gt2,gt3)& 
& ,cplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChaFucSdL = 0._dp 
cplChaFucSdR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChaFucSdT(gt1,gt2,gt3,g2,Yd,Yu,ZD,ZER,ZEL,ZUL,ZUR,cplChaFucSdL(gt1,gt2,gt3)& 
& ,cplChaFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChahhL = 0._dp 
cplcChaChahhR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
  Do gt3 = 1, 8
Call CouplingcChaChahhT(gt1,gt2,gt3,g2,Ye,lam,Yv,ZH,ZER,ZEL,cplcChaChahhL(gt1,gt2,gt3)& 
& ,cplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChaSuL = 0._dp 
cplcFdChaSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 5
  Do gt3 = 1, 6
Call CouplingcFdChaSuT(gt1,gt2,gt3,g2,Yd,Yu,ZU,ZER,ZEL,ZDL,ZDR,cplcFdChaSuL(gt1,gt2,gt3)& 
& ,cplcFdChaSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiChihhL = 0._dp 
cplChiChihhR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingChiChihhT(gt1,gt2,gt3,g1,g2,lam,Yv,kap,ZH,UV,cplChiChihhL(gt1,gt2,gt3)   & 
& ,cplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFdcSdL = 0._dp 
cplChiFdcSdR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFdcSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplChiFdcSdL(gt1,gt2,gt3)   & 
& ,cplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplChiFucSuL = 0._dp 
cplChiFucSuR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CouplingChiFucSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplChiFucSuL(gt1,gt2,gt3)   & 
& ,cplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcChaChiHpmL = 0._dp 
cplcChaChiHpmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
  Do gt3 = 1, 8
Call CouplingcChaChiHpmT(gt1,gt2,gt3,g1,g2,Ye,lam,Yv,ZP,UV,ZER,ZEL,cplcChaChiHpmL(gt1,gt2,gt3)& 
& ,cplcChaChiHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdChiSdL = 0._dp 
cplcFdChiSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFdChiSdT(gt1,gt2,gt3,g1,g2,Yd,ZD,UV,ZDL,ZDR,cplcFdChiSdL(gt1,gt2,gt3)   & 
& ,cplcFdChiSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFuChiSuL = 0._dp 
cplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 10
  Do gt3 = 1, 6
Call CouplingcFuChiSuT(gt1,gt2,gt3,g1,g2,Yu,ZU,UV,ZUL,ZUR,cplcFuChiSuL(gt1,gt2,gt3)   & 
& ,cplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFdcSdL = 0._dp 
cplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFdcSdT(gt2,gt3,g3,pG,ZD,ZDL,ZDR,cplGluFdcSdL(gt2,gt3),cplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


cplcFdFdhhL = 0._dp 
cplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFdhhT(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,cplcFdFdhhL(gt1,gt2,gt3)              & 
& ,cplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplGluFucSuL = 0._dp 
cplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CouplingGluFucSuT(gt2,gt3,g3,pG,ZU,ZUL,ZUR,cplGluFucSuL(gt2,gt3),cplGluFucSuR(gt2,gt3))

 End Do 
End Do 


cplcFuFuhhL = 0._dp 
cplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFuFuhhT(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,cplcFuFuhhL(gt1,gt2,gt3)              & 
& ,cplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdFuHpmL = 0._dp 
cplcFdFuHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 8
Call CouplingcFdFuHpmT(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,cplcFdFuHpmL(gt1,gt2,gt3) & 
& ,cplcFdFuHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


cplcFdGluSdL = 0._dp 
cplcFdGluSdR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFdGluSdT(gt1,gt3,g3,pG,ZD,ZDL,ZDR,cplcFdGluSdL(gt1,gt3),cplcFdGluSdR(gt1,gt3))

 End Do 
End Do 


cplcFuGluSuL = 0._dp 
cplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CouplingcFuGluSuT(gt1,gt3,g3,pG,ZU,ZUL,ZUR,cplcFuGluSuL(gt1,gt3),cplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


cplChiChacVWmL = 0._dp 
cplChiChacVWmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
Call CouplingChiChacVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplChiChacVWmL(gt1,gt2),               & 
& cplChiChacVWmR(gt1,gt2))

 End Do 
End Do 


cplcChaChaVZL = 0._dp 
cplcChaChaVZR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 5
Call CouplingcChaChaVZT(gt1,gt2,g1,g2,ZER,ZEL,TW,cplcChaChaVZL(gt1,gt2)               & 
& ,cplcChaChaVZR(gt1,gt2))

 End Do 
End Do 


cplChiChiVZL = 0._dp 
cplChiChiVZR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 10
Call CouplingChiChiVZT(gt1,gt2,g1,g2,UV,TW,cplChiChiVZL(gt1,gt2),cplChiChiVZR(gt1,gt2))

 End Do 
End Do 


cplcChaChiVWmL = 0._dp 
cplcChaChiVWmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
Call CouplingcChaChiVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplcChaChiVWmL(gt1,gt2),               & 
& cplcChaChiVWmR(gt1,gt2))

 End Do 
End Do 


cplcFdFdVZL = 0._dp 
cplcFdFdVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFdVZT(gt1,gt2,g1,g2,TW,cplcFdFdVZL(gt1,gt2),cplcFdFdVZR(gt1,gt2))

 End Do 
End Do 


cplcFdFuVWmL = 0._dp 
cplcFdFuVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFuVWmT(gt1,gt2,g2,ZDL,ZUL,cplcFdFuVWmL(gt1,gt2),cplcFdFuVWmR(gt1,gt2))

 End Do 
End Do 


cplcFuFuVZL = 0._dp 
cplcFuFuVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFuFuVZT(gt1,gt2,g1,g2,TW,cplcFuFuVZL(gt1,gt2),cplcFuFuVZR(gt1,gt2))

 End Do 
End Do 


Iname = Iname - 1 
 
End subroutine CouplingsFor_Chi_decays_3B
 
Function NFlav(m_in) 
Implicit None 
Real(dp), Intent(in) :: m_in 
Real(dp) :: NFlav 
If (m_in.lt.mf_d(3)) Then 
  NFlav = 4._dp 
Else If (m_in.lt.mf_u(3)) Then 
  NFlav = 5._dp 
Else 
  NFlav = 6._dp 
End if 
End Function

Subroutine RunSM(scale_out,deltaM,tb,g1,g2,g3,Yu, Yd, Ye, vd, vu) 
Implicit None
Real(dp), Intent(in) :: scale_out,deltaM, tb
Real(dp), Intent(out) :: g1, g2, g3, vd, vu
Complex(dp), Intent(out) :: Yu(3,3), Yd(3,3), Ye(3,3)
Real(dp) :: dt, gSM(14), gSM2(2), gSM3(3), mtopMS,  sinw2, vev, tz, alphaStop 
Integer :: kont

Yd = 0._dp
Ye = 0._dp
Yu = 0._dp

If (.not.RunningTopMZ) Then

! Calculating alpha_S(m_top)
gSM2(1)=sqrt(Alpha_mZ*4*Pi) 
gSM2(2)=sqrt(AlphaS_mZ*4*Pi) 


tz=Log(sqrt(mz2)/mf_u(3)) 
dt=tz/50._dp 
Call odeint(gSM2,2,tz,0._dp,deltaM,dt,0._dp,RGEAlphaS,kont)



alphaStop = gSM2(2)**2/4._dp/Pi



! m_top^pole to m_top^MS(m_top) 

mtopMS = mf_u(3)*(1._dp - 4._dp/3._dp*alphaStop/Pi)


! Running m_top^MS(m_top) to M_Z 

gSM3(1)=gSM2(1) 
gSM3(2)=gSM2(2)
gSM3(3)=mtopMS

tz=Log(sqrt(mz2)/mf_u(3)) 
dt=tz/50._dp 
Call odeint(gSM3,3,0._dp,tz,deltaM,dt,0._dp,RGEtop,kont)


mf_u_mz_running = gSM3(3)


RunningTopMZ = .True.

End if

! Starting values at MZ

gSM(1)=sqrt(Alpha_mZ*4*Pi) 
gSM(2)=sqrt(AlphaS_mZ*4*Pi) 
gSM(3)= 0.486E-03_dp ! mf_l_mz(1) 
gSM(4)= 0.10272 !mf_l_mz(2) 
gSM(5)= 1.74624 !mf_l_mz(3) 
gSM(6)= 1.27E-03_dp ! mf_u_mz(1) 
gSM(7)= 0.619  ! mf_u_mz(2) 
gSM(8)= mf_u_mz_running ! m_top 
gSM(9)= 2.9E-03_dp !mf_d_mz(1) 
gSM(10)= 0.055 !mf_d_mz(2) 
gSM(11)= 2.85 ! mf_d_mz(3) 
 

! To get the running sin(w2) 
SinW2 = 0.22290_dp
gSM(12) = 5._dp/3._dp*Alpha_MZ/(1-sinW2)
gSM(13) = Alpha_MZ/Sinw2
gSM(14) = AlphaS_mZ

  nUp =2._dp 
  nDown =3._dp 
  nLep =3._dp 
 

If (scale_out.gt.sqrt(mz2)) Then

 ! From M_Z to Min(M_top,scale_out) 
 If (scale_out.gt.mf_u(3)) Then 
  tz=Log(sqrt(mz2)/mf_u(3)) 
  dt=tz/50._dp 
 Else 
  tz=Log(sqrt(mz2)/scale_out) 
  dt=tz/50._dp 
 End if 

  Call odeint(gSM,14,tz,0._dp,deltaM,dt,0._dp,rge11,kont)


 ! From M_top to M_SUSY if M_top < M_SUSY 
 If (scale_out.gt.mf_u(3)) Then 
  tz=Log(mf_u(3)/scale_out) 
  dt=tz/50._dp 
  nUp =3._dp 
  Call odeint(gSM,14,tz,0._dp,deltaM,dt,0._dp,rge11,kont)
 End if 
Else

 ! From M_Z down to scale_out
  tz=Log(scale_out/sqrt(mz2)) 
  dt=tz/50._dp 
  Call odeint(gSM,14,0._dp,tz,deltaM,dt,0._dp,rge11_SMa,kont)

End if

! Calculating Couplings 

 sinW2=1._dp-mW2/mZ2 
 vev=Sqrt(mZ2*(1._dp-sinW2)*SinW2/(gSM(1)**2/4._dp))
 vd=vev/Sqrt(1._dp+tb**2)
 vu=tb*vd
 
Yd(1,1) =gSM(9)*sqrt(2._dp)/vd 
Yd(2,2) =gSM(10)*sqrt(2._dp)/vd 
Yd(3,3) =gSM(11)*sqrt(2._dp)/vd 
Ye(1,1) =gSM(3)*sqrt(2._dp)/vd 
Ye(2,2)=gSM(4)*sqrt(2._dp)/vd 
Ye(3,3)=gSM(5)*sqrt(2._dp)/vd 
Yu(1,1)=gSM(6)*sqrt(2._dp)/vu 
Yu(2,2)=gSM(7)*sqrt(2._dp)/vu 
Yu(3,3)=gSM(8)*sqrt(2._dp)/vu 


g3 =gSM(2) 
g3running=gSM(2) 

g1 = sqrt(gSM(12)*4._dp*Pi*3._dp/5._dp)
g2 = sqrt(gSM(13)*4._dp*Pi)
! g3 = sqrt(gSM(3)*4._dp*Pi)

sinw2 = g1**2/(g1**2 + g2**2)

!g2=gSM(1)/sqrt(sinW2) 
!g1 = g2*Sqrt(sinW2/(1._dp-sinW2)) 

If (GenerationMixing) Then 
If (TransposedYukawa) Then ! check, if superpotential is Yu Hu u q  or Yu Hu q u
 Yu= Matmul(Transpose(CKM),Transpose(Yu))
Else 
 Yu=Transpose(Matmul(Transpose(CKM),Transpose(Yu)))
End if 
End If


End Subroutine RunSM


Subroutine RunSMohdm(scale_out,deltaM,g1,g2,g3,Yu, Yd, Ye, v) 
Implicit None
Real(dp), Intent(in) :: scale_out,deltaM
Real(dp), Intent(out) :: g1, g2, g3, v
Complex(dp), Intent(out) :: Yu(3,3), Yd(3,3), Ye(3,3)
Real(dp) :: dt, gSM(14), gSM2(2), gSM3(3), mtopMS,  sinw2, vev, tz, alphaStop 
Integer :: kont

Yd = 0._dp
Ye = 0._dp
Yu = 0._dp

If (.not.RunningTopMZ) Then

! Calculating alpha_S(m_top)
gSM2(1)=sqrt(Alpha_mZ*4*Pi) 
gSM2(2)=sqrt(AlphaS_mZ*4*Pi) 


tz=Log(sqrt(mz2)/mf_u(3)) 
dt=tz/50._dp 
Call odeint(gSM2,2,tz,0._dp,deltaM,dt,0._dp,RGEAlphaS,kont)



alphaStop = gSM2(2)**2/4._dp/Pi



! m_top^pole to m_top^MS(m_top) 

mtopMS = mf_u(3)*(1._dp - 4._dp/3._dp*alphaStop/Pi)


! Running m_top^MS(m_top) to M_Z 

gSM3(1)=gSM2(1) 
gSM3(2)=gSM2(2)
gSM3(3)=mtopMS

tz=Log(sqrt(mz2)/mf_u(3)) 
dt=tz/50._dp 
Call odeint(gSM3,3,0._dp,tz,deltaM,dt,0._dp,RGEtop,kont)


mf_u_mz_running = gSM3(3)


RunningTopMZ = .True.

End if

! Starting values at MZ

gSM(1)=sqrt(Alpha_mZ*4*Pi) 
gSM(2)=sqrt(AlphaS_mZ*4*Pi) 
gSM(3)= 0.486E-03_dp ! mf_l_mz(1) 
gSM(4)= 0.10272 !mf_l_mz(2) 
gSM(5)= 1.74624 !mf_l_mz(3) 
gSM(6)= 1.27E-03_dp ! mf_u_mz(1) 
gSM(7)= 0.619  ! mf_u_mz(2) 
gSM(8)= mf_u_mz_running ! m_top 
gSM(9)= 2.9E-03_dp !mf_d_mz(1) 
gSM(10)= 0.055 !mf_d_mz(2) 
gSM(11)= 2.85 ! mf_d_mz(3) 
 

! To get the running sin(w2) 
SinW2 = 0.22290_dp
gSM(12) = 5._dp/3._dp*Alpha_MZ/(1-sinW2)
gSM(13) = Alpha_MZ/Sinw2
gSM(14) = AlphaS_mZ

  nUp =2._dp 
  nDown =3._dp 
  nLep =3._dp 
 

If (scale_out.gt.sqrt(mz2)) Then

 ! From M_Z to Min(M_top,scale_out) 
 If (scale_out.gt.mf_u(3)) Then 
  tz=Log(sqrt(mz2)/mf_u(3)) 
  dt=tz/50._dp 
 Else 
  tz=Log(sqrt(mz2)/scale_out) 
  dt=tz/50._dp 
 End if 

  Call odeint(gSM,14,tz,0._dp,deltaM,dt,0._dp,rge11,kont)


 ! From M_top to M_SUSY if M_top < M_SUSY 
 If (scale_out.gt.mf_u(3)) Then 
  tz=Log(mf_u(3)/scale_out) 
  dt=tz/50._dp 
  nUp =3._dp 
  Call odeint(gSM,14,tz,0._dp,deltaM,dt,0._dp,rge11,kont)
 End if 
Else

 ! From M_Z down to scale_out
  tz=Log(scale_out/sqrt(mz2)) 
  dt=tz/50._dp 
  Call odeint(gSM,14,0._dp,tz,deltaM,dt,0._dp,rge11_SMa,kont)

End if

! Calculating Couplings 

 sinW2=1._dp-mW2/mZ2 
 vev=Sqrt(mZ2*(1._dp-sinW2)*SinW2/(gSM(1)**2/4._dp))
 v = vev
 
Yd(1,1) =gSM(9)*sqrt(2._dp)/v 
Yd(2,2) =gSM(10)*sqrt(2._dp)/v 
Yd(3,3) =gSM(11)*sqrt(2._dp)/v 
Ye(1,1) =gSM(3)*sqrt(2._dp)/v 
Ye(2,2)=gSM(4)*sqrt(2._dp)/v 
Ye(3,3)=gSM(5)*sqrt(2._dp)/v 
Yu(1,1)=gSM(6)*sqrt(2._dp)/v 
Yu(2,2)=gSM(7)*sqrt(2._dp)/v 
Yu(3,3)=gSM(8)*sqrt(2._dp)/v 


g3 =gSM(2) 
g3running=gSM(2) 

g1 = sqrt(gSM(12)*4._dp*Pi*3._dp/5._dp)
g2 = sqrt(gSM(13)*4._dp*Pi)
! g3 = sqrt(gSM(3)*4._dp*Pi)

sinw2 = g1**2/(g1**2 + g2**2)

g2=gSM(1)/sqrt(sinW2) 
g1 = g2*Sqrt(sinW2/(1._dp-sinW2)) 

If (GenerationMixing) Then 
If (TransposedYukawa) Then ! check, if superpotential is Yu Hu u q  or Yu Hu q u
 Yu= Matmul(Transpose(CKM),Transpose(Yu))
Else 
 Yu=Transpose(Matmul(Transpose(CKM),Transpose(Yu)))
End if 
End If


End Subroutine RunSMohdm
End Module CouplingsForDecays_munuSSM3G
