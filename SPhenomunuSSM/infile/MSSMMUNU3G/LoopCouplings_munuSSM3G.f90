! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:24 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module LoopCouplings_munuSSM3G 
 
Use Control 
Use Couplings_munuSSM3G 
Use Mathematics 
Use LoopFunctions 
Use StandardModel 
 
 Contains 
 
Real(dp) Function Alpha_MSbar(Q,mW,mt) 
Implicit None 
Real(dp),Intent(in)::Q,mW 
Real(dp),Intent(in),Optional::mt 
Real(dp)::DeltaAlpha 
If (MZ_input) Then 
Alpha_MSbar=Alpha_mZ_MS 
If (Present(mt)) Then 
DeltaAlpha=-8._dp*Log(Q/mt)/(9._dp*Pi) 
Alpha_MSbar=Alpha_MSbar/(1._dp+DeltaAlpha*alpha) 
End If 
Else 
DeltaAlpha=3.5_dp*Log(Q/mW)/Pi+0.5_dp*oo3pi 
If (Present(mt)) DeltaAlpha=DeltaAlpha-8._dp*Log(Q/mt)/(9._dp*Pi) 
Alpha_MSbar=Alpha/(1._dp-Delta_Alpha_Lepton-Delta_Alpha_Hadron& 
&+DeltaAlpha*alpha) 
Alpha_MZ_MS=Alpha_MSbar 
End If 
End Function Alpha_MSbar
 
 
Real(dp) Function AlphaEwDR(Q,MVWm,MSd,MSu,MHpm,MCha,MFd,MFu) 
 
Real(dp),Intent(in)::Q,MVWm,MSd(6),MSu(6),MHpm(8),MCha(5),MFd(3),MFu(3)
Integer::i1 
Real(dp)::DeltaAlpha 
If (MZ_input) then 
DeltaAlpha=1._dp-Alpha/Alpha_MZ_MS! MSbar value^=mW+light fermions 
DeltaAlpha=DeltaAlpha+alpha/(6._dp*Pi)! conversion to DRbar 
If (MVWm.gt.Q) Then 
DeltaAlpha=DeltaAlpha-1._dp/3._dp*Log(MVWm/ Q)*Alpha/(2._dp*pi) 
End if 
Do i1=1,6
If (MSd(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha-1._dp/9._dp*Log(MSd(i1)/ Q)*Alpha/(2._dp*pi)  
End if 
End Do 
Do i1=1,6
If (MSu(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha-4._dp/9._dp*Log(MSu(i1)/ Q)*Alpha/(2._dp*pi)  
End if 
End Do 
Do i1=2,8
If (MHpm(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha-1._dp/3._dp*Log(MHpm(i1)/ Q)*Alpha/(2._dp*pi)  
End if 
End Do 
Do i1=1,5
If (MCha(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha-4._dp/3._dp*Log(MCha(i1)/ Q)*Alpha/(2._dp*pi)  
End if 
End Do 
Do i1=1,3
If (MFd(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha-4._dp/9._dp*Log(MFd(i1)/ Q)*Alpha/(2._dp*pi)  
End if 
End Do 
Do i1=1,3
If (MFu(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha-16._dp/9._dp*Log(MFu(i1)/ Q)*Alpha/(2._dp*pi)  
End if 
End Do 
Else 
DeltaAlpha=7._dp*Log(Q/mW) !+16._dp*Log(mf_u(3)/Q)/9._dp 
If (MVWm.gt.Q) Then 
DeltaAlpha=DeltaAlpha+1._dp/3._dp*Log(MVWm/ Q)
End if 
Do i1=1,6
If (MSd(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha+1._dp/9._dp*Log(MSd(i1)/ Q) 
End if 
End Do 
Do i1=1,6
If (MSu(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha+4._dp/9._dp*Log(MSu(i1)/ Q) 
End if 
End Do 
Do i1=2,8
If (MHpm(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha+1._dp/3._dp*Log(MHpm(i1)/ Q) 
End if 
End Do 
Do i1=1,5
If (MCha(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha+4._dp/3._dp*Log(MCha(i1)/ Q) 
End if 
End Do 
Do i1=1,3
If (MFd(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha+4._dp/9._dp*Log(MFd(i1)/ Q) 
End if 
End Do 
Do i1=1,3
If (MFu(i1).gt.Q) Then 
DeltaAlpha=DeltaAlpha+16._dp/9._dp*Log(MFu(i1)/ Q) 
End if 
End Do 
DeltaAlpha=Delta_Alpha_Lepton+Delta_Alpha_Hadron& 
    &-alpha*DeltaAlpha/(2._dp*Pi) 
End If 
 
AlphaEwDR=Alpha/(1._dp-DeltaAlpha) 
 
End Function AlphaEwDR 
 
 
Real(dp) Function AlphaSDR(Q,MSd,MSu,MGlu,MFd,MFu) 
Real(dp),Intent(in)::Q,MSd(6),MSu(6),MGlu,MFd(3),MFu(3) 
Integer::i1 
Real(dp)::DeltaAlpha 
DeltaAlpha = 0.5_dp !- 2._dp*Log(sqrt(mf_u2(3))/Q)/3._dp 
Do i1=1,6
 If (Abs(MSd(i1)/ Q).gt.1._dp) Then 
  DeltaAlpha=DeltaAlpha-1._dp/6._dp*Log(MSd(i1)/ Q) 
 End If 
End Do 
Do i1=1,6
 If (Abs(MSu(i1)/ Q).gt.1._dp) Then 
  DeltaAlpha=DeltaAlpha-1._dp/6._dp*Log(MSu(i1)/ Q) 
 End If 
End Do 
DeltaAlpha=DeltaAlpha-2*Log(MGlu/ Q) 
Do i1=1,3
 If (Abs(MFd(i1)/ Q).gt.1._dp) Then 
  DeltaAlpha=DeltaAlpha-2._dp/3._dp*Log(MFd(i1)/ Q) 
 End If 
End Do 
Do i1=1,3
 If (Abs(MFu(i1)/ Q).gt.1._dp) Then 
  DeltaAlpha=DeltaAlpha-2._dp/3._dp*Log(MFu(i1)/ Q) 
 End If 
End Do 
DeltaAlpha=AlphaS_mZ*DeltaAlpha/(2._dp*Pi) 
   AlphaSDR=AlphaS_mZ/(1._dp-DeltaAlpha)
 
End Function AlphaSDR 
Real(dp) Function AlphaEW_T(AlphaEW_In, Q,MVWm,MSd,MSu,MHpm,MCha,MFd,MFu) 
 
Real(dp),Intent(in)::AlphaEW_In, Q,MVWm,MSd(6),MSu(6),MHpm(8),MCha(5),MFd(3),MFu(3)
Integer::i1 
Real(dp)::DeltaAlpha 
DeltaAlpha=1._dp/6._dp 
Do i1=1,6
DeltaAlpha=DeltaAlpha+1._dp/9._dp*Log(MSd(i1)/ Q) 
End Do 
Do i1=1,6
DeltaAlpha=DeltaAlpha+4._dp/9._dp*Log(MSu(i1)/ Q) 
End Do 
Do i1=2,8
DeltaAlpha=DeltaAlpha+1._dp/3._dp*Log(MHpm(i1)/ Q) 
End Do 
DeltaAlpha=-alpha_in*DeltaAlpha/(2._dp*Pi) 
AlphaEW_T=Alpha_in/(1._dp-DeltaAlpha) 
 
End Function AlphaEW_T 
 
 
Real(dp) Function AlphaS_T(AlphaS_In, Q,MSd,MSu,MGlu,MFd,MFu) 
Real(dp),Intent(in):: AlphaS_In, Q,MSd(6),MSu(6),MGlu,MFd(3),MFu(3) 
Integer::i1 
Real(dp)::DeltaAlpha 
DeltaAlpha = 0.5_dp 
Do i1=1,6
 If (Abs(MSd(i1)/mf_u(3)).gt.1._dp) Then 
DeltaAlpha=DeltaAlpha-1._dp/6._dp*Log(MSd(i1)/ Q) 
 End If 
End Do 
Do i1=1,6
 If (Abs(MSu(i1)/mf_u(3)).gt.1._dp) Then 
DeltaAlpha=DeltaAlpha-1._dp/6._dp*Log(MSu(i1)/ Q) 
 End If 
End Do 
DeltaAlpha=DeltaAlpha-2*Log(MGlu/ Q) 
Do i1=1,3
 If (Abs(MFd(i1)/mf_u(3)).gt.1._dp) Then 
DeltaAlpha=DeltaAlpha-2._dp/3._dp*Log(MFd(i1)/ Q) 
 End If 
End Do 
Do i1=1,3
 If (Abs(MFu(i1)/mf_u(3)).gt.1._dp) Then 
DeltaAlpha=DeltaAlpha-2._dp/3._dp*Log(MFu(i1)/ Q) 
 End If 
End Do 
DeltaAlpha=AlphaS_In*DeltaAlpha/(2._dp*Pi) 
   AlphaS_T=AlphaS_In/(1._dp-DeltaAlpha)
 
End Function AlphaS_T



Subroutine DeltaVB(sinW2,sinW2_dr,rho,MAh,MCha,MChi,MFd,MFu,Mhh,MHpm,MSd,             & 
& MSu,MVWm,g1,g2,ZER,ZEL,vd,vL,vu,Yd,Ye,Yu,Yv,ZA,ZD,ZDL,ZDR,ZH,UVinput,ZP,               & 
& ZU,ZUL,ZUR,kap,lam,res)

Implicit None 
Real(dp),Intent(in) :: MAh(8),MCha(5),MChi(10),MFd(3),MFu(3),Mhh(8),MHpm(8),MSd(6),MSu(6),MVWm,              & 
& g1,g2,vd,vL(3),vu,ZA(8,8),ZH(8,8),ZP(8,8)

Complex(dp),Intent(in) :: ZER(5,5),ZEL(5,5),Yd(3,3),Ye(3,3),Yu(3,3),Yv(3,3),ZD(6,6),ZDL(3,3),ZDR(3,3),          & 
& UVinput(10,10),ZU(6,6),ZUL(3,3),ZUR(3,3),kap(3,3,3),lam(3)

Real(dp) :: MAh2(8),MCha2(5),MChi2(10),MFd2(3),MFu2(3),Mhh2(8),MHpm2(8),MSd2(6),MSu2(6),MVWm2

Complex(dp) :: cplAhcHpmVWm(8,8),cplcChacFuSdL(5,3,6),cplcChacFuSdR(5,3,6),cplcChaChaAhL(5,5,8),     & 
& cplcChaChaAhR(5,5,8),cplcChaChahhL(5,5,8),cplcChaChahhR(5,5,8),cplcChaChiHpmL(5,10,8), & 
& cplcChaChiHpmR(5,10,8),cplcChaChiVWmL(5,10),cplcChaChiVWmR(5,10),cplcChaFdcSuL(5,3,6), & 
& cplcChaFdcSuR(5,3,6),cplcFdChaSuL(3,5,6),cplcFdChaSuR(3,5,6),cplcFdChiSdL(3,10,6),     & 
& cplcFdChiSdR(3,10,6),cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),cplcFuChiSuL(3,10,6),         & 
& cplcFuChiSuR(3,10,6),cplChaFucSdL(5,3,6),cplChaFucSdR(5,3,6),cplChiChacHpmL(10,5,8),   & 
& cplChiChacHpmR(10,5,8),cplChiChacVWmL(10,5),cplChiChacVWmR(10,5),cplChiChiAhL(10,10,8),& 
& cplChiChiAhR(10,10,8),cplChiChihhL(10,10,8),cplChiChihhR(10,10,8),cplChiFdcSdL(10,3,6),& 
& cplChiFdcSdR(10,3,6),cplChiFucSuL(10,3,6),cplChiFucSuR(10,3,6),cplhhcHpmVWm(8,8),      & 
& cplhhcVWmVWm(8),cplSucSdVWm(6,6)

Complex(dp) :: UV(10,10) 
Integer :: i1,i2,i3,i4,gt1,gt2,gt3,gt4 
Real(dp), Intent(in) :: sinW2,sinW2_Dr, rho 
Real(dp), Intent(out) :: res 
 
Complex(dp) :: sumI, coup1L,coup1R,coup2L,coup2R,coup3L,coup3R,coup3, coup4L,coup4R, teil 
Complex(dp) :: D27m2, D0m2, vertex, phase 
Real(dp) :: cosW2, cosW2_Dr, chargefactor 
Iname = Iname+1
NameOfUnit(Iname) = "DeltaVB" 
MAh2 = MAh**2 
MCha2 = MCha**2 
MChi2 = MChi**2 
MFd2 = MFd**2 
MFu2 = MFu**2 
Mhh2 = Mhh**2 
MHpm2 = MHpm**2 
MSd2 = MSd**2 
MSu2 = MSu**2 
MVWm2 = MVWm**2 

 
 ! Fix neutrino phases 
 
Do i1=1,10
  phase=(1._dp,0._dp) 
  If (Abs(UVinput(i1,i1)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,i1))/Abs(UVinput(i1,i1)) 
     Else If (Abs(UVinput(i1,1)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,1))/Abs(UVinput(i1,1)) 
     Else If (Abs(UVinput(i1,2)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,2))/Abs(UVinput(i1,2)) 
     Else If (Abs(UVinput(i1,3)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,3))/Abs(UVinput(i1,3)) 
     Else If (Abs(UVinput(i1,4)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,4))/Abs(UVinput(i1,4)) 
     Else If (Abs(UVinput(i1,5)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,5))/Abs(UVinput(i1,5)) 
     Else If (Abs(UVinput(i1,6)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,6))/Abs(UVinput(i1,6)) 
     Else If (Abs(UVinput(i1,7)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,7))/Abs(UVinput(i1,7)) 
     Else If (Abs(UVinput(i1,8)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,8))/Abs(UVinput(i1,8)) 
     Else If (Abs(UVinput(i1,9)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,9))/Abs(UVinput(i1,9)) 
     Else If (Abs(UVinput(i1,10)).Ne.0._dp) Then 
     phase=Conjg(UVinput(i1,10))/Abs(UVinput(i1,10)) 
  End If 
 UV(i1,:) = UVinput(i1,:)*phase 
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


cplAhcHpmVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplingAhcHpmVWmT(gt1,gt2,g2,ZA,ZP,cplAhcHpmVWm(gt1,gt2))

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


cplChiChacVWmL = 0._dp 
cplChiChacVWmR = 0._dp 
Do gt1 = 1, 10
 Do gt2 = 1, 5
Call CouplingChiChacVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplChiChacVWmL(gt1,gt2),               & 
& cplChiChacVWmR(gt1,gt2))

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


cplcChaChiVWmL = 0._dp 
cplcChaChiVWmR = 0._dp 
Do gt1 = 1, 5
 Do gt2 = 1, 10
Call CouplingcChaChiVWmT(gt1,gt2,g2,UV,ZER,ZEL,cplcChaChiVWmL(gt1,gt2),               & 
& cplcChaChiVWmR(gt1,gt2))

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


cplcFdFuVWmL = 0._dp 
cplcFdFuVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CouplingcFdFuVWmT(gt1,gt2,g2,ZDL,ZUL,cplcFdFuVWmL(gt1,gt2),cplcFdFuVWmR(gt1,gt2))

 End Do 
End Do 


cplhhcHpmVWm = 0._dp 
Do gt1 = 1, 8
 Do gt2 = 1, 8
Call CouplinghhcHpmVWmT(gt1,gt2,g2,ZH,ZP,cplhhcHpmVWm(gt1,gt2))

 End Do 
End Do 


cplhhcVWmVWm = 0._dp 
Do gt1 = 1, 8
Call CouplinghhcVWmVWmT(gt1,g2,vd,vu,vL,ZH,cplhhcVWmVWm(gt1))

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


cplSucSdVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CouplingSucSdVWmT(gt1,gt2,g2,ZD,ZU,cplSucSdVWm(gt1,gt2))

 End Do 
End Do 


!-------------------------- 
!SM Contributions 
!-------------------------- 
cosW2 = 1._dp - sinW2 
cosW2_DR = 1._dp - sinW2_DR 
sumI = 6._dp & 
  & + Log(cosW2)*(3.5_dp - 2.5_dp *sinW2   & 
  & - sinW2_DR*(5._dp - 1.5_dp*cosW2/cosW2_DR))/sinW2   
res = sumI*g2**2*rho 
 
 
If (IncludeBSMdeltaVB) Then 
!-------------------------- 
!BSM  Contributions 
!-------------------------- 
teil = 0._dp 
 
Do gt1=1,10 
 Do gt2=1,10 
sumI = 0._dp 
 
If (mf_l2(2).gt. 0.5_dp*MChi2(gt2)) Then 
chargefactor = 1 
Do i1=1,8
  Do i2=1,10
If ((MAh2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i2,i1)
coup1R = cplChiChiAhR(gt1,i2,i1)
coup2R = Conjg(cplChiChiAhL(gt2,i2,i1))
coup2L = Conjg(cplChiChiAhR(gt2,i2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MChi2(i2),MAh2(i1))  
End if 
   End Do
  End Do




chargefactor = 1 
Do i1=1,5
  Do i2=1,8
If ((MCha2(i1).gt.mf_l2(2)).Or.(MHpm2(i2).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i1,i2)
coup1R = cplChiChacHpmR(gt1,i1,i2)
coup2R = Conjg(cplChiChacHpmL(gt2,i1,i2))
coup2L = Conjg(cplChiChacHpmR(gt2,i1,i2))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MCha2(i1),MHpm2(i2))  
End if 
   End Do
  End Do




chargefactor = 1 
Do i1=1,10
  Do i2=1,8
If ((MChi2(i1).gt.mf_l2(2)).Or.(Mhh2(i2).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i1,i2)
coup1R = cplChiChihhR(gt1,i1,i2)
coup2R = Conjg(cplChiChihhL(gt2,i1,i2))
coup2L = Conjg(cplChiChihhR(gt2,i1,i2))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MChi2(i1),Mhh2(i2))  
End if 
   End Do
  End Do




chargefactor = 3 
Do i1=1,3
  Do i2=1,6
If ((MFd2(i1).gt.mf_l2(2)).Or.(MSd2(i2).gt.mf_l2(2))) Then
coup1L = cplChiFdcSdL(gt1,i1,i2)
coup1R = cplChiFdcSdR(gt1,i1,i2)
coup2R = Conjg(cplChiFdcSdL(gt2,i1,i2))
coup2L = Conjg(cplChiFdcSdR(gt2,i1,i2))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MFd2(i1),MSd2(i2))  
End if 
   End Do
  End Do




chargefactor = 3 
Do i1=1,3
  Do i2=1,6
If ((MFu2(i1).gt.mf_l2(2)).Or.(MSu2(i2).gt.mf_l2(2))) Then
coup1L = cplChiFucSuL(gt1,i1,i2)
coup1R = cplChiFucSuR(gt1,i1,i2)
coup2R = Conjg(cplChiFucSuL(gt2,i1,i2))
coup2L = Conjg(cplChiFucSuR(gt2,i1,i2))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MFu2(i1),MSu2(i2))  
End if 
   End Do
  End Do




chargefactor = 1 
Do i1=1,8
  Do i2=1,5
If ((MHpm2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2))) Then
coup1L = cplcChaChiHpmL(i2,gt1,i1)
coup1R = cplcChaChiHpmR(i2,gt1,i1)
coup2R = Conjg(cplcChaChiHpmL(i2,gt2,i1))
coup2L = Conjg(cplcChaChiHpmR(i2,gt2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MCha2(i2),MHpm2(i1))  
End if 
   End Do
  End Do




chargefactor = 3 
Do i1=1,6
  Do i2=1,3
If ((MSd2(i1).gt.mf_l2(2)).Or.(MFd2(i2).gt.mf_l2(2))) Then
coup1L = cplcFdChiSdL(i2,gt1,i1)
coup1R = cplcFdChiSdR(i2,gt1,i1)
coup2R = Conjg(cplcFdChiSdL(i2,gt2,i1))
coup2L = Conjg(cplcFdChiSdR(i2,gt2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MFd2(i2),MSd2(i1))  
End if 
   End Do
  End Do




chargefactor = 3 
Do i1=1,6
  Do i2=1,3
If ((MSu2(i1).gt.mf_l2(2)).Or.(MFu2(i2).gt.mf_l2(2))) Then
coup1L = cplcFuChiSuL(i2,gt1,i1)
coup1R = cplcFuChiSuR(i2,gt1,i1)
coup2R = Conjg(cplcFuChiSuL(i2,gt2,i1))
coup2L = Conjg(cplcFuChiSuR(i2,gt2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MFu2(i2),MSu2(i1))  
End if 
   End Do
  End Do




res = res + sumI*(UV(gt2,1)+UV(gt2,2)) 
End if 
End Do 
 
End Do 
 
Do gt1=1,2 
Do  gt2=1,5 
sumI = 0._dp 
 
chargefactor = 1 
Do i1=1,8
  Do i2=1,5
If ((MAh2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2))) Then
coup1L = cplcChaChaAhL(i2,gt1,i1)
coup1R = cplcChaChaAhR(i2,gt1,i1)
coup2R = Conjg(cplcChaChaAhL(i2,gt2,i1))
coup2L = Conjg(cplcChaChaAhR(i2,gt2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MCha2(i2),MAh2(i1))  
End if 
   End Do
  End Do




chargefactor = 1 
Do i1=1,10
  Do i2=1,8
If ((MChi2(i1).gt.mf_l2(2)).Or.(MHpm2(i2).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(i1,gt1,i2)
coup1R = cplChiChacHpmR(i1,gt1,i2)
coup2R = Conjg(cplChiChacHpmL(i1,gt2,i2))
coup2L = Conjg(cplChiChacHpmR(i1,gt2,i2))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MChi2(i1),MHpm2(i2))  
End if 
   End Do
  End Do




chargefactor = 3 
Do i1=1,3
  Do i2=1,6
If ((MFu2(i1).gt.mf_l2(2)).Or.(MSd2(i2).gt.mf_l2(2))) Then
coup1L = cplChaFucSdL(gt1,i1,i2)
coup1R = cplChaFucSdR(gt1,i1,i2)
coup2R = Conjg(cplChaFucSdL(gt2,i1,i2))
coup2L = Conjg(cplChaFucSdR(gt2,i1,i2))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MFu2(i1),MSd2(i2))  
End if 
   End Do
  End Do




chargefactor = 1 
Do i1=1,8
  Do i2=1,5
If ((Mhh2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2))) Then
coup1L = cplcChaChahhL(i2,gt1,i1)
coup1R = cplcChaChahhR(i2,gt1,i1)
coup2R = Conjg(cplcChaChahhL(i2,gt2,i1))
coup2L = Conjg(cplcChaChahhR(i2,gt2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MCha2(i2),Mhh2(i1))  
End if 
   End Do
  End Do




chargefactor = 3 
Do i1=1,6
  Do i2=1,3
If ((MSu2(i1).gt.mf_l2(2)).Or.(MFd2(i2).gt.mf_l2(2))) Then
coup1L = cplcFdChaSuL(i2,gt1,i1)
coup1R = cplcFdChaSuR(i2,gt1,i1)
coup2R = Conjg(cplcFdChaSuL(i2,gt2,i1))
coup2L = Conjg(cplcFdChaSuR(i2,gt2,i1))
sumI = sumI + chargefactor*0.5_dp*coup1L*coup2R*B1(0._dp,MFd2(i2),MSu2(i1))  
End if 
   End Do
  End Do




res = res + sumI 
End Do 
 
End Do 
 
vertex = 0._dp 
 
Do gt1=1,10 
 Do gt2=1,2 
chargefactor = 1 
Do i1= 1,8
  Do i2= 1,10
   Do i3= 1,5
  If ((MAh2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2)).Or.(MCha2(i3).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i2,i1)
coup1R = cplChiChiAhR(gt1,i2,i1)
coup2L = cplcChaChaAhL(gt2,i3,i1)
coup2R = cplcChaChaAhR(gt2,i3,i1)
coup3L = -cplcChaChiVWmR(i3,i2)
coup3R = -cplcChaChiVWmL(i3,i2)
vertex = vertex + chargefactor*(coup1L*coup2R*(-sqrt2*coup3R*MChi(i2)*MCha(i3)& 
& *C0_3m(MAh2(i1),MChi2(i2),MCha2(i3)) + oosqrt2*coup3L* & 
& (B0(0._dp,MChi2(i2),MCha2(i3))-0.5_dp +MAh2(i1)*C0_3m(MAh2(i1),MChi2(i2),MCha2(i3))))) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,10
  Do i2= 1,8
   Do i3= 1,8
  If ((MChi2(i1).gt.mf_l2(2)).Or.(MAh2(i2).gt.mf_l2(2)).Or.(MHpm2(i3).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i1,i2)
coup1R = cplChiChiAhR(gt1,i1,i2)
coup2L = cplcChaChiHpmL(gt2,i1,i3)
coup2R = cplcChaChiHpmR(gt2,i1,i3)
coup3 = cplAhcHpmVWm(i2,i3)
vertex = vertex + chargefactor*(0.5_dp*sqrt2*coup1L*coup2R*coup3*(MChi2(i1)*C0_3m(MChi2(i1),MAh2(i2),MHpm2(i3)) + B0(0._dp,MAh2(i2),MHpm2(i3)) +0.5_dp)) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,10
  Do i2= 1,8
   Do i3= 1,8
  If ((MChi2(i1).gt.mf_l2(2)).Or.(Mhh2(i2).gt.mf_l2(2)).Or.(MHpm2(i3).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i1,i2)
coup1R = cplChiChihhR(gt1,i1,i2)
coup2L = cplcChaChiHpmL(gt2,i1,i3)
coup2R = cplcChaChiHpmR(gt2,i1,i3)
coup3 = cplhhcHpmVWm(i2,i3)
vertex = vertex + chargefactor*(0.5_dp*sqrt2*coup1L*coup2R*coup3*(MChi2(i1)*C0_3m(MChi2(i1),Mhh2(i2),MHpm2(i3)) + B0(0._dp,Mhh2(i2),MHpm2(i3)) +0.5_dp)) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,10
  Do i2= 1,8
  If ((MChi2(i1).gt.mf_l2(2)).Or.(Mhh2(i2).gt.mf_l2(2)).Or.(MVWm2.gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i1,i2)
coup1R = cplChiChihhR(gt1,i1,i2)
coup2L = cplcChaChiVWmL(gt2,i1)
coup2R = cplcChaChiVWmR(gt2,i1)
coup3 = cplhhcVWmVWm(i2)
End if 
   End Do
  End Do


chargefactor = 3 
Do i1= 1,3
  Do i2= 1,6
   Do i3= 1,6
  If ((MFu2(i1).gt.mf_l2(2)).Or.(MSu2(i2).gt.mf_l2(2)).Or.(MSd2(i3).gt.mf_l2(2))) Then
coup1L = cplChiFucSuL(gt1,i1,i2)
coup1R = cplChiFucSuR(gt1,i1,i2)
coup2L = cplcChacFuSdL(gt2,i1,i3)
coup2R = cplcChacFuSdR(gt2,i1,i3)
coup3 = cplSucSdVWm(i2,i3)
vertex = vertex + chargefactor*(0.5_dp*sqrt2*coup1L*coup2R*coup3*(MFu2(i1)*C0_3m(MFu2(i1),MSu2(i2),MSd2(i3)) + B0(0._dp,MSu2(i2),MSd2(i3)) +0.5_dp)) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,8
  Do i2= 1,10
   Do i3= 1,5
  If ((Mhh2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2)).Or.(MCha2(i3).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i2,i1)
coup1R = cplChiChihhR(gt1,i2,i1)
coup2L = cplcChaChahhL(gt2,i3,i1)
coup2R = cplcChaChahhR(gt2,i3,i1)
coup3L = -cplcChaChiVWmR(i3,i2)
coup3R = -cplcChaChiVWmL(i3,i2)
vertex = vertex + chargefactor*(coup1L*coup2R*(-sqrt2*coup3R*MChi(i2)*MCha(i3)& 
& *C0_3m(Mhh2(i1),MChi2(i2),MCha2(i3)) + oosqrt2*coup3L* & 
& (B0(0._dp,MChi2(i2),MCha2(i3))-0.5_dp +Mhh2(i1)*C0_3m(Mhh2(i1),MChi2(i2),MCha2(i3))))) 
End if 
   End Do
  End Do
End Do


chargefactor = 3 
Do i1= 1,6
  Do i2= 1,3
   Do i3= 1,3
  If ((MSu2(i1).gt.mf_l2(2)).Or.(MFu2(i2).gt.mf_l2(2)).Or.(MFd2(i3).gt.mf_l2(2))) Then
coup1L = cplcFuChiSuL(i2,gt1,i1)
coup1R = cplcFuChiSuR(i2,gt1,i1)
coup2L = cplcChaFdcSuL(gt2,i3,i1)
coup2R = cplcChaFdcSuR(gt2,i3,i1)
coup3L = -cplcFdFuVWmR(i3,i2)
coup3R = -cplcFdFuVWmL(i3,i2)
vertex = vertex + chargefactor*(coup1L*coup2R*(-sqrt2*coup3R*MFu(i2)*MFd(i3)& 
& *C0_3m(MSu2(i1),MFu2(i2),MFd2(i3)) + oosqrt2*coup3L* & 
& (B0(0._dp,MFu2(i2),MFd2(i3))-0.5_dp +MSu2(i1)*C0_3m(MSu2(i1),MFu2(i2),MFd2(i3))))) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,5
  Do i2= 1,8
   Do i3= 1,8
  If ((MCha2(i1).gt.mf_l2(2)).Or.(MHpm2(i2).gt.mf_l2(2)).Or.(MAh2(i3).gt.mf_l2(2))) Then
coup1L = cplcChaChiHpmL(i1,gt1,i2)
coup1R = cplcChaChiHpmR(i1,gt1,i2)
coup2L = cplcChaChaAhL(gt2,i1,i3)
coup2R = cplcChaChaAhR(gt2,i1,i3)
coup3 = -cplAhcHpmVWm(i3,i2)
vertex = vertex + chargefactor*(0.5_dp*sqrt2*coup1L*coup2R*coup3*(MCha2(i1)*C0_3m(MCha2(i1),MHpm2(i2),MAh2(i3)) + B0(0._dp,MHpm2(i2),MAh2(i3)) +0.5_dp)) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,5
  Do i2= 1,8
   Do i3= 1,8
  If ((MCha2(i1).gt.mf_l2(2)).Or.(MHpm2(i2).gt.mf_l2(2)).Or.(Mhh2(i3).gt.mf_l2(2))) Then
coup1L = cplcChaChiHpmL(i1,gt1,i2)
coup1R = cplcChaChiHpmR(i1,gt1,i2)
coup2L = cplcChaChahhL(gt2,i1,i3)
coup2R = cplcChaChahhR(gt2,i1,i3)
coup3 = -cplhhcHpmVWm(i3,i2)
vertex = vertex + chargefactor*(0.5_dp*sqrt2*coup1L*coup2R*coup3*(MCha2(i1)*C0_3m(MCha2(i1),MHpm2(i2),Mhh2(i3)) + B0(0._dp,MHpm2(i2),Mhh2(i3)) +0.5_dp)) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,5
   Do i3= 1,8
  If ((MCha2(i1).gt.mf_l2(2)).Or.(MVWm2.gt.mf_l2(2)).Or.(Mhh2(i3).gt.mf_l2(2))) Then
coup1L = cplcChaChiVWmL(i1,gt1)
coup1R = cplcChaChiVWmR(i1,gt1)
coup2L = cplcChaChahhL(gt2,i1,i3)
coup2R = cplcChaChahhR(gt2,i1,i3)
coup3 = cplhhcVWmVWm(i3)
End if 
   End Do
End Do


chargefactor = 3 
Do i1= 1,3
  Do i2= 1,6
   Do i3= 1,6
  If ((MFd2(i1).gt.mf_l2(2)).Or.(MSd2(i2).gt.mf_l2(2)).Or.(MSu2(i3).gt.mf_l2(2))) Then
coup1L = cplcFdChiSdL(i1,gt1,i2)
coup1R = cplcFdChiSdR(i1,gt1,i2)
coup2L = cplcChaFdcSuL(gt2,i1,i3)
coup2R = cplcChaFdcSuR(gt2,i1,i3)
coup3 = -cplSucSdVWm(i3,i2)
vertex = vertex + chargefactor*(0.5_dp*sqrt2*coup1L*coup2R*coup3*(MFd2(i1)*C0_3m(MFd2(i1),MSd2(i2),MSu2(i3)) + B0(0._dp,MSd2(i2),MSu2(i3)) +0.5_dp)) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
Do i1= 1,8
  Do i2= 1,5
   Do i3= 1,10
  If ((MHpm2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2)).Or.(MChi2(i3).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i2,i1)
coup1R = cplChiChacHpmR(gt1,i2,i1)
coup2L = cplcChaChiHpmL(gt2,i3,i1)
coup2R = cplcChaChiHpmR(gt2,i3,i1)
coup3L = cplcChaChiVWmL(i2,i3)
coup3R = cplcChaChiVWmR(i2,i3)
vertex = vertex + chargefactor*(coup1L*coup2R*(-sqrt2*coup3R*MCha(i2)*MChi(i3)& 
& *C0_3m(MHpm2(i1),MCha2(i2),MChi2(i3)) + oosqrt2*coup3L* & 
& (B0(0._dp,MCha2(i2),MChi2(i3))-0.5_dp +MHpm2(i1)*C0_3m(MHpm2(i1),MCha2(i2),MChi2(i3))))) 
End if 
   End Do
  End Do
End Do


chargefactor = 3 
Do i1= 1,6
  Do i2= 1,3
   Do i3= 1,3
  If ((MSd2(i1).gt.mf_l2(2)).Or.(MFd2(i2).gt.mf_l2(2)).Or.(MFu2(i3).gt.mf_l2(2))) Then
coup1L = cplChiFdcSdL(gt1,i2,i1)
coup1R = cplChiFdcSdR(gt1,i2,i1)
coup2L = cplcChacFuSdL(gt2,i3,i1)
coup2R = cplcChacFuSdR(gt2,i3,i1)
coup3L = cplcFdFuVWmL(i2,i3)
coup3R = cplcFdFuVWmR(i2,i3)
vertex = vertex + chargefactor*(coup1L*coup2R*(-sqrt2*coup3R*MFd(i2)*MFu(i3)& 
& *C0_3m(MSd2(i1),MFd2(i2),MFu2(i3)) + oosqrt2*coup3L* & 
& (B0(0._dp,MFd2(i2),MFu2(i3))-0.5_dp +MSd2(i1)*C0_3m(MSd2(i1),MFd2(i2),MFu2(i3))))) 
End if 
   End Do
  End Do
End Do


chargefactor = 1 
  Do i2= 1,5
   Do i3= 1,10
  If ((MVWm2.gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2)).Or.(MChi2(i3).gt.mf_l2(2))) Then
coup1L = cplChiChacVWmL(gt1,i2)
coup1R = cplChiChacVWmR(gt1,i2)
coup2L = cplcChaChiVWmL(gt2,i3)
coup2R = cplcChaChiVWmR(gt2,i3)
coup3L = cplcChaChiVWmL(i2,i3)
coup3R = cplcChaChiVWmR(i2,i3)
End if 
  End Do
End Do


 End Do 
 
End Do 
 
res = res + vertex/g2 
Do gt1=1,10 
 Do gt2=1,10 
gt3 = 2 
gt4 = 1 
! Cha,Hpm,Chi,Hpm
chargefactor = 1 
Do i1=1,5
  Do i2=1,8
    Do i3=1,10
      Do i4=1,8
If ((MCha2(i1).gt.mf_l2(2)).Or.(MHpm2(i2).gt.mf_l2(2)).Or.(MChi2(i3).gt.mf_l2(2)).Or.(MHpm2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i1,i4)
coup1R = cplChiChacHpmR(gt1,i1,i4)
coup2L = cplcChaChiHpmL(i1,gt2,i2)
coup2R = cplcChaChiHpmR(i1,gt2,i2)
coup3L = cplChiChacHpmL(i3,gt3,i2)
coup3R = cplChiChacHpmR(i3,gt3,i2)
coup4L = cplcChaChiHpmL(gt4,i3,i4)
coup4R = cplcChaChiHpmR(gt4,i3,i4)
D27m2 = D27_Bagger(MHpm2(i2),MHpm2(i4),MCha2(i1),MChi2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Chi,Ah,bar[Cha],Ah
chargefactor = 1 
Do i1=1,10
  Do i2=1,8
    Do i3=1,5
      Do i4=1,8
If ((MChi2(i1).gt.mf_l2(2)).Or.(MAh2(i2).gt.mf_l2(2)).Or.(MCha2(i3).gt.mf_l2(2)).Or.(MAh2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i1,i4)
coup1R = cplChiChiAhR(gt1,i1,i4)
coup2L = cplChiChiAhL(gt2,i1,i2)
coup2R = cplChiChiAhR(gt2,i1,i2)
coup3L = cplcChaChaAhL(i3,gt3,i2)
coup3R = cplcChaChaAhR(i3,gt3,i2)
coup4L = cplcChaChaAhL(gt4,i3,i4)
coup4R = cplcChaChaAhR(gt4,i3,i4)
D27m2 = D27_Bagger(MAh2(i2),MAh2(i4),MChi2(i1),MCha2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Chi,Ah,bar[Cha],hh
chargefactor = 1 
Do i1=1,10
  Do i2=1,8
    Do i3=1,5
      Do i4=1,8
If ((MChi2(i1).gt.mf_l2(2)).Or.(MAh2(i2).gt.mf_l2(2)).Or.(MCha2(i3).gt.mf_l2(2)).Or.(Mhh2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i1,i4)
coup1R = cplChiChihhR(gt1,i1,i4)
coup2L = cplChiChiAhL(gt2,i1,i2)
coup2R = cplChiChiAhR(gt2,i1,i2)
coup3L = cplcChaChaAhL(i3,gt3,i2)
coup3R = cplcChaChaAhR(i3,gt3,i2)
coup4L = cplcChaChahhL(gt4,i3,i4)
coup4R = cplcChaChahhR(gt4,i3,i4)
D27m2 = D27_Bagger(MAh2(i2),Mhh2(i4),MChi2(i1),MCha2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Chi,hh,bar[Cha],Ah
chargefactor = 1 
Do i1=1,10
  Do i2=1,8
    Do i3=1,5
      Do i4=1,8
If ((MChi2(i1).gt.mf_l2(2)).Or.(Mhh2(i2).gt.mf_l2(2)).Or.(MCha2(i3).gt.mf_l2(2)).Or.(MAh2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i1,i4)
coup1R = cplChiChiAhR(gt1,i1,i4)
coup2L = cplChiChihhL(gt2,i1,i2)
coup2R = cplChiChihhR(gt2,i1,i2)
coup3L = cplcChaChahhL(i3,gt3,i2)
coup3R = cplcChaChahhR(i3,gt3,i2)
coup4L = cplcChaChaAhL(gt4,i3,i4)
coup4R = cplcChaChaAhR(gt4,i3,i4)
D27m2 = D27_Bagger(Mhh2(i2),MAh2(i4),MChi2(i1),MCha2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Chi,hh,bar[Cha],hh
chargefactor = 1 
Do i1=1,10
  Do i2=1,8
    Do i3=1,5
      Do i4=1,8
If ((MChi2(i1).gt.mf_l2(2)).Or.(Mhh2(i2).gt.mf_l2(2)).Or.(MCha2(i3).gt.mf_l2(2)).Or.(Mhh2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i1,i4)
coup1R = cplChiChihhR(gt1,i1,i4)
coup2L = cplChiChihhL(gt2,i1,i2)
coup2R = cplChiChihhR(gt2,i1,i2)
coup3L = cplcChaChahhL(i3,gt3,i2)
coup3R = cplcChaChahhR(i3,gt3,i2)
coup4L = cplcChaChahhL(gt4,i3,i4)
coup4R = cplcChaChahhR(gt4,i3,i4)
D27m2 = D27_Bagger(Mhh2(i2),Mhh2(i4),MChi2(i1),MCha2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Fd,Sd,Fu,Sd
chargefactor = 3 
Do i1=1,3
  Do i2=1,6
    Do i3=1,3
      Do i4=1,6
If ((MFd2(i1).gt.mf_l2(2)).Or.(MSd2(i2).gt.mf_l2(2)).Or.(MFu2(i3).gt.mf_l2(2)).Or.(MSd2(i4).gt.mf_l2(2))) Then
coup1L = cplChiFdcSdL(gt1,i1,i4)
coup1R = cplChiFdcSdR(gt1,i1,i4)
coup2L = cplcFdChiSdL(i1,gt2,i2)
coup2R = cplcFdChiSdR(i1,gt2,i2)
coup3L = cplChaFucSdL(gt3,i3,i2)
coup3R = cplChaFucSdR(gt3,i3,i2)
coup4L = cplcChacFuSdL(gt4,i3,i4)
coup4R = cplcChacFuSdR(gt4,i3,i4)
D27m2 = D27_Bagger(MSd2(i2),MSd2(i4),MFd2(i1),MFu2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! bar[Fu],conj[Su],bar[Fd],conj[Su]
chargefactor = 3 
Do i1=1,3
  Do i2=1,6
    Do i3=1,3
      Do i4=1,6
If ((MFu2(i1).gt.mf_l2(2)).Or.(MSu2(i2).gt.mf_l2(2)).Or.(MFd2(i3).gt.mf_l2(2)).Or.(MSu2(i4).gt.mf_l2(2))) Then
coup1L = cplcFuChiSuL(i1,gt1,i4)
coup1R = cplcFuChiSuR(i1,gt1,i4)
coup2L = cplChiFucSuL(gt2,i1,i2)
coup2R = cplChiFucSuR(gt2,i1,i2)
coup3L = cplcFdChaSuL(i3,gt3,i2)
coup3R = cplcFdChaSuR(i3,gt3,i2)
coup4L = cplcChaFdcSuL(gt4,i3,i4)
coup4R = cplcChaFdcSuR(gt4,i3,i4)
D27m2 = D27_Bagger(MSu2(i2),MSu2(i4),MFu2(i1),MFd2(i3))
If(Real(D27m2,dp).eq.Real(D27m2,dp)) Then 
teil = teil + D27m2*chargefactor*coup1L*coup2R*coup3L*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Ah,Chi,Hpm,Chi
chargefactor = 1 
Do i1=1,8
  Do i2=1,10
    Do i3=1,8
      Do i4=1,10
If ((MAh2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2)).Or.(MHpm2(i3).gt.mf_l2(2)).Or.(MChi2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i4,i1)
coup1R = cplChiChiAhR(gt1,i4,i1)
coup2L = cplChiChiAhL(gt2,i2,i1)
coup2R = cplChiChiAhR(gt2,i2,i1)
coup3L = cplcChaChiHpmL(gt4,i2,i3)
coup3R = cplcChaChiHpmR(gt4,i2,i3)
coup4L = cplChiChacHpmL(i4,gt3,i3)
coup4R = cplChiChacHpmR(i4,gt3,i3)
D0m2 = D0_Bagger(MAh2(i1),MHpm2(i3),MChi2(i2),MChi2(i4))*MChi(i2)*MChi(i4) 
D27m2 = D27_Bagger(MAh2(i1),MHpm2(i3),MChi2(i2),MChi2(i4))
If ((Real(D27m2,dp).eq.Real(D27m2,dp)).And.(Real(D0m2,dp).eq.Real(D0m2,dp))) Then 
teil = teil + 0.5_dp*chargefactor*D27m2*coup1L*coup2R*coup3L*coup4R+D0m2*coup1L*coup2L*coup3R*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! hh,Chi,Hpm,Chi
chargefactor = 1 
Do i1=1,8
  Do i2=1,10
    Do i3=1,8
      Do i4=1,10
If ((Mhh2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2)).Or.(MHpm2(i3).gt.mf_l2(2)).Or.(MChi2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i4,i1)
coup1R = cplChiChihhR(gt1,i4,i1)
coup2L = cplChiChihhL(gt2,i2,i1)
coup2R = cplChiChihhR(gt2,i2,i1)
coup3L = cplcChaChiHpmL(gt4,i2,i3)
coup3R = cplcChaChiHpmR(gt4,i2,i3)
coup4L = cplChiChacHpmL(i4,gt3,i3)
coup4R = cplChiChacHpmR(i4,gt3,i3)
D0m2 = D0_Bagger(Mhh2(i1),MHpm2(i3),MChi2(i2),MChi2(i4))*MChi(i2)*MChi(i4) 
D27m2 = D27_Bagger(Mhh2(i1),MHpm2(i3),MChi2(i2),MChi2(i4))
If ((Real(D27m2,dp).eq.Real(D27m2,dp)).And.(Real(D0m2,dp).eq.Real(D0m2,dp))) Then 
teil = teil + 0.5_dp*chargefactor*D27m2*coup1L*coup2R*coup3L*coup4R+D0m2*coup1L*coup2L*coup3R*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Su,Fu,Sd,Fu
chargefactor = 3 
Do i1=1,6
  Do i2=1,3
    Do i3=1,6
      Do i4=1,3
If ((MSu2(i1).gt.mf_l2(2)).Or.(MFu2(i2).gt.mf_l2(2)).Or.(MSd2(i3).gt.mf_l2(2)).Or.(MFu2(i4).gt.mf_l2(2))) Then
coup1L = cplcFuChiSuL(i4,gt1,i1)
coup1R = cplcFuChiSuR(i4,gt1,i1)
coup2L = cplChiFucSuL(gt2,i2,i1)
coup2R = cplChiFucSuR(gt2,i2,i1)
coup3L = cplcChacFuSdL(gt4,i2,i3)
coup3R = cplcChacFuSdR(gt4,i2,i3)
coup4L = cplChaFucSdL(gt3,i4,i3)
coup4R = cplChaFucSdR(gt3,i4,i3)
D0m2 = D0_Bagger(MSu2(i1),MSd2(i3),MFu2(i2),MFu2(i4))*MFu(i2)*MFu(i4) 
D27m2 = D27_Bagger(MSu2(i1),MSd2(i3),MFu2(i2),MFu2(i4))
If ((Real(D27m2,dp).eq.Real(D27m2,dp)).And.(Real(D0m2,dp).eq.Real(D0m2,dp))) Then 
teil = teil + 0.5_dp*chargefactor*D27m2*coup1L*coup2R*coup3L*coup4R+D0m2*coup1L*coup2L*coup3R*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! conj[Hpm],bar[Cha],Ah,bar[Cha]
chargefactor = 1 
Do i1=1,8
  Do i2=1,5
    Do i3=1,8
      Do i4=1,5
If ((MHpm2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2)).Or.(MAh2(i3).gt.mf_l2(2)).Or.(MCha2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i4,i1)
coup1R = cplChiChacHpmR(gt1,i4,i1)
coup2L = cplcChaChiHpmL(i2,gt2,i1)
coup2R = cplcChaChiHpmR(i2,gt2,i1)
coup3L = cplcChaChaAhL(gt4,i2,i3)
coup3R = cplcChaChaAhR(gt4,i2,i3)
coup4L = cplcChaChaAhL(i4,gt3,i3)
coup4R = cplcChaChaAhR(i4,gt3,i3)
D0m2 = D0_Bagger(MHpm2(i1),MAh2(i3),MCha2(i2),MCha2(i4))*MCha(i2)*MCha(i4) 
D27m2 = D27_Bagger(MHpm2(i1),MAh2(i3),MCha2(i2),MCha2(i4))
If ((Real(D27m2,dp).eq.Real(D27m2,dp)).And.(Real(D0m2,dp).eq.Real(D0m2,dp))) Then 
teil = teil + 0.5_dp*chargefactor*D27m2*coup1L*coup2R*coup3L*coup4R+D0m2*coup1L*coup2L*coup3R*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! conj[Hpm],bar[Cha],hh,bar[Cha]
chargefactor = 1 
Do i1=1,8
  Do i2=1,5
    Do i3=1,8
      Do i4=1,5
If ((MHpm2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2)).Or.(Mhh2(i3).gt.mf_l2(2)).Or.(MCha2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i4,i1)
coup1R = cplChiChacHpmR(gt1,i4,i1)
coup2L = cplcChaChiHpmL(i2,gt2,i1)
coup2R = cplcChaChiHpmR(i2,gt2,i1)
coup3L = cplcChaChahhL(gt4,i2,i3)
coup3R = cplcChaChahhR(gt4,i2,i3)
coup4L = cplcChaChahhL(i4,gt3,i3)
coup4R = cplcChaChahhR(i4,gt3,i3)
D0m2 = D0_Bagger(MHpm2(i1),Mhh2(i3),MCha2(i2),MCha2(i4))*MCha(i2)*MCha(i4) 
D27m2 = D27_Bagger(MHpm2(i1),Mhh2(i3),MCha2(i2),MCha2(i4))
If ((Real(D27m2,dp).eq.Real(D27m2,dp)).And.(Real(D0m2,dp).eq.Real(D0m2,dp))) Then 
teil = teil + 0.5_dp*chargefactor*D27m2*coup1L*coup2R*coup3L*coup4R+D0m2*coup1L*coup2L*coup3R*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! conj[Sd],bar[Fd],conj[Su],bar[Fd]
chargefactor = 3 
Do i1=1,6
  Do i2=1,3
    Do i3=1,6
      Do i4=1,3
If ((MSd2(i1).gt.mf_l2(2)).Or.(MFd2(i2).gt.mf_l2(2)).Or.(MSu2(i3).gt.mf_l2(2)).Or.(MFd2(i4).gt.mf_l2(2))) Then
coup1L = cplChiFdcSdL(gt1,i4,i1)
coup1R = cplChiFdcSdR(gt1,i4,i1)
coup2L = cplcFdChiSdL(i2,gt2,i1)
coup2R = cplcFdChiSdR(i2,gt2,i1)
coup3L = cplcChaFdcSuL(gt4,i2,i3)
coup3R = cplcChaFdcSuR(gt4,i2,i3)
coup4L = cplcFdChaSuL(i4,gt3,i3)
coup4R = cplcFdChaSuR(i4,gt3,i3)
D0m2 = D0_Bagger(MSd2(i1),MSu2(i3),MFd2(i2),MFd2(i4))*MFd(i2)*MFd(i4) 
D27m2 = D27_Bagger(MSd2(i1),MSu2(i3),MFd2(i2),MFd2(i4))
If ((Real(D27m2,dp).eq.Real(D27m2,dp)).And.(Real(D0m2,dp).eq.Real(D0m2,dp))) Then 
teil = teil + 0.5_dp*chargefactor*D27m2*coup1L*coup2R*coup3L*coup4R+D0m2*coup1L*coup2L*coup3R*coup4R 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Ah,Cha,Hpm,Chi
chargefactor = 1 
Do i1=1,8
  Do i2=1,5
    Do i3=1,8
      Do i4=1,10
If ((MAh2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2)).Or.(MHpm2(i3).gt.mf_l2(2)).Or.(MChi2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChiAhL(gt1,i4,i1)
coup1R = cplChiChiAhR(gt1,i4,i1)
coup2L = cplcChaChaAhL(gt4,i2,i1)
coup2R = cplcChaChaAhR(gt4,i2,i1)
coup3L = cplcChaChiHpmL(i2,gt2,i3)
coup3R = cplcChaChiHpmR(i2,gt2,i3)
coup4L = cplChiChacHpmL(i4,gt3,i3)
coup4R = cplChiChacHpmR(i4,gt3,i3)
D0m2 = D0_Bagger(MAh2(i1),MHpm2(i3),MCha2(i2),MChi2(i4))*MCha(i2)*MChi(i4) 
If (Real(D0m2,dp).eq.Real(D0m2,dp)) Then 
teil = teil + 0.5_dp*chargefactor*D0m2*coup1L*coup2R*coup3R*coup4L 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! hh,Cha,Hpm,Chi
chargefactor = 1 
Do i1=1,8
  Do i2=1,5
    Do i3=1,8
      Do i4=1,10
If ((Mhh2(i1).gt.mf_l2(2)).Or.(MCha2(i2).gt.mf_l2(2)).Or.(MHpm2(i3).gt.mf_l2(2)).Or.(MChi2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChihhL(gt1,i4,i1)
coup1R = cplChiChihhR(gt1,i4,i1)
coup2L = cplcChaChahhL(gt4,i2,i1)
coup2R = cplcChaChahhR(gt4,i2,i1)
coup3L = cplcChaChiHpmL(i2,gt2,i3)
coup3R = cplcChaChiHpmR(i2,gt2,i3)
coup4L = cplChiChacHpmL(i4,gt3,i3)
coup4R = cplChiChacHpmR(i4,gt3,i3)
D0m2 = D0_Bagger(Mhh2(i1),MHpm2(i3),MCha2(i2),MChi2(i4))*MCha(i2)*MChi(i4) 
If (Real(D0m2,dp).eq.Real(D0m2,dp)) Then 
teil = teil + 0.5_dp*chargefactor*D0m2*coup1L*coup2R*coup3R*coup4L 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! Su,Fd,Sd,Fu
chargefactor = 3 
Do i1=1,6
  Do i2=1,3
    Do i3=1,6
      Do i4=1,3
If ((MSu2(i1).gt.mf_l2(2)).Or.(MFd2(i2).gt.mf_l2(2)).Or.(MSd2(i3).gt.mf_l2(2)).Or.(MFu2(i4).gt.mf_l2(2))) Then
coup1L = cplcFuChiSuL(i4,gt1,i1)
coup1R = cplcFuChiSuR(i4,gt1,i1)
coup2L = cplcChaFdcSuL(gt4,i2,i1)
coup2R = cplcChaFdcSuR(gt4,i2,i1)
coup3L = cplcFdChiSdL(i2,gt2,i3)
coup3R = cplcFdChiSdR(i2,gt2,i3)
coup4L = cplChaFucSdL(gt3,i4,i3)
coup4R = cplChaFucSdR(gt3,i4,i3)
D0m2 = D0_Bagger(MSu2(i1),MSd2(i3),MFd2(i2),MFu2(i4))*MFd(i2)*MFu(i4) 
If (Real(D0m2,dp).eq.Real(D0m2,dp)) Then 
teil = teil + 0.5_dp*chargefactor*D0m2*coup1L*coup2R*coup3R*coup4L 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! conj[Hpm],Chi,Ah,bar[Cha]
chargefactor = 1 
Do i1=1,8
  Do i2=1,10
    Do i3=1,8
      Do i4=1,5
If ((MHpm2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2)).Or.(MAh2(i3).gt.mf_l2(2)).Or.(MCha2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i4,i1)
coup1R = cplChiChacHpmR(gt1,i4,i1)
coup2L = cplcChaChiHpmL(gt4,i2,i1)
coup2R = cplcChaChiHpmR(gt4,i2,i1)
coup3L = cplChiChiAhL(gt2,i2,i3)
coup3R = cplChiChiAhR(gt2,i2,i3)
coup4L = cplcChaChaAhL(i4,gt3,i3)
coup4R = cplcChaChaAhR(i4,gt3,i3)
D0m2 = D0_Bagger(MHpm2(i1),MAh2(i3),MChi2(i2),MCha2(i4))*MChi(i2)*MCha(i4) 
If (Real(D0m2,dp).eq.Real(D0m2,dp)) Then 
teil = teil + 0.5_dp*chargefactor*D0m2*coup1L*coup2R*coup3R*coup4L 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! conj[Hpm],Chi,hh,bar[Cha]
chargefactor = 1 
Do i1=1,8
  Do i2=1,10
    Do i3=1,8
      Do i4=1,5
If ((MHpm2(i1).gt.mf_l2(2)).Or.(MChi2(i2).gt.mf_l2(2)).Or.(Mhh2(i3).gt.mf_l2(2)).Or.(MCha2(i4).gt.mf_l2(2))) Then
coup1L = cplChiChacHpmL(gt1,i4,i1)
coup1R = cplChiChacHpmR(gt1,i4,i1)
coup2L = cplcChaChiHpmL(gt4,i2,i1)
coup2R = cplcChaChiHpmR(gt4,i2,i1)
coup3L = cplChiChihhL(gt2,i2,i3)
coup3R = cplChiChihhR(gt2,i2,i3)
coup4L = cplcChaChahhL(i4,gt3,i3)
coup4R = cplcChaChahhR(i4,gt3,i3)
D0m2 = D0_Bagger(MHpm2(i1),Mhh2(i3),MChi2(i2),MCha2(i4))*MChi(i2)*MCha(i4) 
If (Real(D0m2,dp).eq.Real(D0m2,dp)) Then 
teil = teil + 0.5_dp*chargefactor*D0m2*coup1L*coup2R*coup3R*coup4L 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


 ! conj[Sd],bar[Fu],conj[Su],bar[Fd]
chargefactor = 3 
Do i1=1,6
  Do i2=1,3
    Do i3=1,6
      Do i4=1,3
If ((MSd2(i1).gt.mf_l2(2)).Or.(MFu2(i2).gt.mf_l2(2)).Or.(MSu2(i3).gt.mf_l2(2)).Or.(MFd2(i4).gt.mf_l2(2))) Then
coup1L = cplChiFdcSdL(gt1,i4,i1)
coup1R = cplChiFdcSdR(gt1,i4,i1)
coup2L = cplcChacFuSdL(gt4,i2,i1)
coup2R = cplcChacFuSdR(gt4,i2,i1)
coup3L = cplChiFucSuL(gt2,i2,i3)
coup3R = cplChiFucSuR(gt2,i2,i3)
coup4L = cplcFdChaSuL(i4,gt3,i3)
coup4R = cplcFdChaSuR(i4,gt3,i3)
D0m2 = D0_Bagger(MSd2(i1),MSu2(i3),MFu2(i2),MFd2(i4))*MFu(i2)*MFd(i4) 
If (Real(D0m2,dp).eq.Real(D0m2,dp)) Then 
teil = teil + 0.5_dp*chargefactor*D0m2*coup1L*coup2R*coup3R*coup4L 
End if
End if 
    End Do 
   End Do 
  End Do 
End Do 


  End Do 
 
End Do 
 

 
sumI = -2._dp*cosW2_DR*mz2*Real(teil,dp)/g2**2 
res = res + SumI 
End if ! BSM part 
res = res*oo16pi2 
Iname = Iname-1
End subroutine DeltaVB 
 
 
Subroutine CoupHiggsToPhoton(mHiggs,inG,ratCha,ratFd,ratFu,ratHpm,ratSd,              & 
& ratSu,ratVWm,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,rq,rsq,coup)

Implicit None 
Complex(dp),Intent(in) :: ratCha(5),ratFd(3),ratFu(3),ratHpm(8),ratSd(6),ratSu(6),ratVWm

Real(dp),Intent(in) :: MCha(5),MFd(3),MFu(3),MHpm(8),MVWm,MSd(6),MSu(6)

Integer, Intent(in) :: inG 
Real(dp), Intent(in) :: mHiggs, rq, rsq 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
coup = coup + 1._dp*1*(-1)**2*ratVWm*A_one(mH2p/MVWm**2)
HPPloopVWm(inG) = HPPloopVWm(inG) + 1._dp*1*(-1)**2*ratVWm*A_one(mH2p/MVWm**2)
Do i1 = 1 , 6
coup = coup + rsq*3*(-1._dp/3._dp)**2*ratSd(i1)*A_zero(mH2p/MSd(i1)**2)
HPPloopSd(i1,inG) = HPPloopSd(i1,inG) + rsq*3*(-1._dp/3._dp)**2*ratSd(i1)*A_zero(mH2p/MSd(i1)**2)
End Do 
Do i1 = 1 , 6
coup = coup + rsq*3*(2._dp/3._dp)**2*ratSu(i1)*A_zero(mH2p/MSu(i1)**2)
HPPloopSu(i1,inG) = HPPloopSu(i1,inG) + rsq*3*(2._dp/3._dp)**2*ratSu(i1)*A_zero(mH2p/MSu(i1)**2)
End Do 
Do i1 = 2 , 8
coup = coup + 1._dp*1*(-1)**2*ratHpm(i1)*A_zero(mH2p/MHpm(i1)**2)
HPPloopHpm(i1,inG) = HPPloopHpm(i1,inG) + 1._dp*1*(-1)**2*ratHpm(i1)*A_zero(mH2p/MHpm(i1)**2)
End Do 
Do i1 = 1 , 5
coup = coup + 1._dp*1*(-1)**2*ratCha(i1)*A_onehalf(mH2p/MCha(i1)**2)
HPPloopCha(i1,inG) = HPPloopCha(i1,inG) + 1._dp*1*(-1)**2*ratCha(i1)*A_onehalf(mH2p/MCha(i1)**2)
End Do 
Do i1 = 1 , 3
If (MFd(i1).gt.mHiggs) Then 
coup = coup + rq*3*(-1._dp/3._dp)**2*ratFd(i1)*A_onehalf(mH2p/MFd(i1)**2)
HPPloopFd(i1,inG) = HPPloopFd(i1,inG) + rq*3*(-1._dp/3._dp)**2*ratFd(i1)*A_onehalf(mH2p/MFd(i1)**2)
Else 
coup = coup + 3*(-1._dp/3._dp)**2*ratFd(i1)*A_onehalf(mH2p/MFd(i1)**2)
HPPloopFd(i1,inG) = HPPloopFd(i1,inG) + 3*(-1._dp/3._dp)**2*ratFd(i1)*A_onehalf(mH2p/MFd(i1)**2)
End if 
End Do 
Do i1 = 1 , 3
If (MFu(i1).gt.mHiggs) Then 
coup = coup + rq*3*(2._dp/3._dp)**2*ratFu(i1)*A_onehalf(mH2p/MFu(i1)**2)
HPPloopFu(i1,inG) = HPPloopFu(i1,inG) + rq*3*(2._dp/3._dp)**2*ratFu(i1)*A_onehalf(mH2p/MFu(i1)**2)
Else 
coup = coup + 3*(2._dp/3._dp)**2*ratFu(i1)*A_onehalf(mH2p/MFu(i1)**2)
HPPloopFu(i1,inG) = HPPloopFu(i1,inG) + 3*(2._dp/3._dp)**2*ratFu(i1)*A_onehalf(mH2p/MFu(i1)**2)
End if 
End Do 
End Subroutine CoupHiggsToPhoton

Subroutine CoupHiggsToGluon(mHiggs,inG,ratFd,ratFu,ratSd,ratSu,MFd,MFu,               & 
& MSd,MSu,rq,rsq,coup)

Implicit None 
Complex(dp),Intent(in) :: ratFd(3),ratFu(3),ratSd(6),ratSu(6)

Real(dp),Intent(in) :: MFd(3),MFu(3),MSd(6),MSu(6)

Integer, Intent(in) :: inG 
Real(dp), Intent(in) :: mHiggs, rq, rsq 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
Do i1 = 1 , 6
coup = coup + rsq*1*ratSd(i1)*A_zero(mH2p/MSd(i1)**2)
End Do 
Do i1 = 1 , 6
coup = coup + rsq*1*ratSu(i1)*A_zero(mH2p/MSu(i1)**2)
End Do 
Do i1 = 1 , 3
coup = coup + rq*1*ratFd(i1)*A_onehalf(mH2p/MFd(i1)**2)
End Do 
Do i1 = 1 , 3
coup = coup + rq*1*ratFu(i1)*A_onehalf(mH2p/MFu(i1)**2)
End Do 
coup = 0.75_dp*coup 
End Subroutine CoupHiggsToGluon



Subroutine CoupHiggsToPhotonSM(mHiggs,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,rq,coup)

Implicit None 
Real(dp),Intent(in) :: MCha(5),MFd(3),MFu(3),MHpm(8),MVWm,MSd(6),MSu(6)

Real(dp), Intent(in) :: mHiggs, rq 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
coup = coup + 1._dp*1*(-1)**2*A_one(mH2p/MVWm**2)
Do i1 =1, 3 
coup = coup + 1._dp*1*(-1)**2*A_onehalf(mH2p/MCha(i1)**2)
End Do 
Do i1 =1, 3 
If (MFd(i1).gt.mHiggs) Then 
coup = coup + rq*3*(-1._dp/3._dp)**2*A_onehalf(mH2p/MFd(i1)**2)
Else 
coup = coup + 3*(-1._dp/3._dp)**2*A_onehalf(mH2p/MFd(i1)**2)
End if 
End Do 
Do i1 =1, 3 
If (MFu(i1).gt.mHiggs) Then 
coup = coup + rq*3*(2._dp/3._dp)**2*A_onehalf(mH2p/MFu(i1)**2)
Else 
coup = coup + 3*(2._dp/3._dp)**2*A_onehalf(mH2p/MFu(i1)**2)
End if 
End Do 
End Subroutine CoupHiggsToPhotonSM 



Subroutine CoupHiggsToGluonSM(mHiggs,MFd,MFu,MSd,MSu,rq,coup)

Implicit None 
Real(dp),Intent(in) :: MFd(3),MFu(3),MSd(6),MSu(6)

Real(dp), Intent(in) :: mHiggs, rq 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
Do i1 =1, 3 
coup = coup + rq*1*A_onehalf(mH2p/MFd(i1)**2)
End Do 
Do i1 =1, 3 
coup = coup + rq*1*A_onehalf(mH2p/MFu(i1)**2)
End Do 
coup = 0.75_dp*coup 
End Subroutine CoupHiggsToGluonSM 

Subroutine CoupPseudoHiggsToPhoton(mHiggs,inG,ratCha,ratFd,ratFu,ratHpm,              & 
& ratSd,ratSu,ratVWm,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,coup)

Implicit None 
Complex(dp),Intent(in) :: ratCha(5),ratFd(3),ratFu(3),ratHpm(8),ratSd(6),ratSu(6),ratVWm

Real(dp),Intent(in) :: MCha(5),MFd(3),MFu(3),MHpm(8),MVWm,MSd(6),MSu(6)

Real(dp), Intent(in) :: mHiggs 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Integer, Intent(in) :: inG 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
Do i1 =1, 5
coup = coup + 1*(-1)**2*ratCha(i1)*AP_onehalf(mH2p/MCha(i1)**2)
End Do 
Do i1 =1, 3
coup = coup + 3*(-1._dp/3._dp)**2*ratFd(i1)*AP_onehalf(mH2p/MFd(i1)**2)
End Do 
Do i1 =1, 3
coup = coup + 3*(2._dp/3._dp)**2*ratFu(i1)*AP_onehalf(mH2p/MFu(i1)**2)
End Do 
End Subroutine CoupPseudoHiggsToPhoton

Subroutine CoupPseudoHiggsToGluon(mHiggs,inG,ratFd,ratFu,ratSd,ratSu,MFd,             & 
& MFu,MSd,MSu,coup)

Implicit None 
Complex(dp),Intent(in) :: ratFd(3),ratFu(3),ratSd(6),ratSu(6)

Real(dp),Intent(in) :: MFd(3),MFu(3),MSd(6),MSu(6)

Real(dp), Intent(in) :: mHiggs 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Integer, Intent(in) :: inG 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
Do i1 =1, 3
coup = coup + 1*ratFd(i1)*AP_onehalf(mH2p/MFd(i1)**2)
End Do 
Do i1 =1, 3
coup = coup + 1*ratFu(i1)*AP_onehalf(mH2p/MFu(i1)**2)
End Do 
End Subroutine CoupPseudoHiggsToGluon

Subroutine CoupPseudoHiggsToPhotonSM(mHiggs,MCha,MFd,MFu,MHpm,MVWm,MSd,MSu,coup)

Implicit None 
Real(dp),Intent(in) :: MCha(5),MFd(3),MFu(3),MHpm(8),MVWm,MSd(6),MSu(6)

Real(dp), Intent(in) :: mHiggs 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
Do i1 =1, 5
coup = coup + 1*(-1)**2**AP_onehalf(mH2p/MCha(i1)**2)
End Do 
Do i1 =1, 3
coup = coup + 3*(-1._dp/3._dp)**2**AP_onehalf(mH2p/MFd(i1)**2)
End Do 
Do i1 =1, 3
coup = coup + 3*(2._dp/3._dp)**2**AP_onehalf(mH2p/MFu(i1)**2)
End Do 
End Subroutine CoupPseudoHiggsToPhotonSM 

Subroutine CoupPseudoHiggsToGluonSM(mHiggs,MFd,MFu,MSd,MSu,coup)

Implicit None 
Real(dp),Intent(in) :: MFd(3),MFu(3),MSd(6),MSu(6)

Real(dp), Intent(in) :: mHiggs 
Complex(dp), Intent(out) :: coup 
Integer :: i1 
Real(dp) :: Mh2p 
mH2p = 0.25_dp*mHiggs**2 
coup = 0._dp 
 
Do i1 =1, 3
coup = coup + 1*AP_onehalf(mH2p/MFd(i1)**2)
End Do 
Do i1 =1, 3
coup = coup + 1*AP_onehalf(mH2p/MFu(i1)**2)
End Do 
End Subroutine CoupPseudoHiggsToGluonSM 

End Module LoopCouplings_munuSSM3G 
 
