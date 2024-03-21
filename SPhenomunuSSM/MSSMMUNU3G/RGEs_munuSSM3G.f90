! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:14 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module RGEs_munuSSM3G 
 
Use Control 
Use Model_Data_munuSSM3G 
Use Mathematics 
 
Logical,Private,Save::OnlyDiagonal

Real(dp),Parameter::id3R(3,3)=& 
   & Reshape(Source=(/& 
   & 1,0,0,& 
 &0,1,0,& 
 &0,0,1& 
 &/),shape=(/3,3/)) 
Contains 


Subroutine GToParameters386(g,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,             & 
& Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3)

Implicit None 
Real(dp), Intent(in) :: g(386) 
Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp),Intent(out) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Integer i1, i2, i3, i4, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'GToParameters386' 
 
g1= g(1) 
g2= g(2) 
g3= g(3) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yd(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Ye(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
lam(i1) = Cmplx( g(SumI+40), g(SumI+41), dp) 
End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yv(i1,i2) = Cmplx( g(SumI+46), g(SumI+47), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yu(i1,i2) = Cmplx( g(SumI+64), g(SumI+65), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
kap(i1,i2,i3) = Cmplx( g(SumI+82), g(SumI+83), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Td(i1,i2) = Cmplx( g(SumI+136), g(SumI+137), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Te(i1,i2) = Cmplx( g(SumI+154), g(SumI+155), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
Tlam(i1) = Cmplx( g(SumI+172), g(SumI+173), dp) 
End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Tv(i1,i2) = Cmplx( g(SumI+178), g(SumI+179), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Tu(i1,i2) = Cmplx( g(SumI+196), g(SumI+197), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
Tk(i1,i2,i3) = Cmplx( g(SumI+214), g(SumI+215), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mq2(i1,i2) = Cmplx( g(SumI+268), g(SumI+269), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
ml2(i1,i2) = Cmplx( g(SumI+286), g(SumI+287), dp) 
End Do 
 End Do 
 
mHd2= g(304) 
mHu2= g(305) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
md2(i1,i2) = Cmplx( g(SumI+306), g(SumI+307), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mu2(i1,i2) = Cmplx( g(SumI+324), g(SumI+325), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
me2(i1,i2) = Cmplx( g(SumI+342), g(SumI+343), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mv2(i1,i2) = Cmplx( g(SumI+360), g(SumI+361), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
mlHd2(i1) =  g(SumI+378) 
End Do 
 
M1= Cmplx(g(381),g(382),dp) 
M2= Cmplx(g(383),g(384),dp) 
M3= Cmplx(g(385),g(386),dp) 
Do i1=1,386 
If (g(i1).ne.g(i1)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Write(*,*) "At position ", i1 
 Call TerminateProgram 
 return
End if 
End do 
Iname = Iname - 1 
 
End Subroutine GToParameters386

Subroutine ParametersToG386(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,               & 
& Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,g)

Implicit None 
Real(dp), Intent(out) :: g(386) 
Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3)

Complex(dp), Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Integer i1, i2, i3, i4, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'ParametersToG386' 
 
g(1) = g1  
g(2) = g2  
g(3) = g3  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+4) = Real(Yd(i1,i2), dp) 
g(SumI+5) = Aimag(Yd(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+22) = Real(Ye(i1,i2), dp) 
g(SumI+23) = Aimag(Ye(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
g(SumI+40) = Real(lam(i1), dp) 
g(SumI+41) = Aimag(lam(i1)) 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+46) = Real(Yv(i1,i2), dp) 
g(SumI+47) = Aimag(Yv(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+64) = Real(Yu(i1,i2), dp) 
g(SumI+65) = Aimag(Yu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
g(SumI+82) = Real(kap(i1,i2,i3), dp) 
g(SumI+83) = Aimag(kap(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+136) = Real(Td(i1,i2), dp) 
g(SumI+137) = Aimag(Td(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+154) = Real(Te(i1,i2), dp) 
g(SumI+155) = Aimag(Te(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
g(SumI+172) = Real(Tlam(i1), dp) 
g(SumI+173) = Aimag(Tlam(i1)) 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+178) = Real(Tv(i1,i2), dp) 
g(SumI+179) = Aimag(Tv(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+196) = Real(Tu(i1,i2), dp) 
g(SumI+197) = Aimag(Tu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
g(SumI+214) = Real(Tk(i1,i2,i3), dp) 
g(SumI+215) = Aimag(Tk(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+268) = Real(mq2(i1,i2), dp) 
g(SumI+269) = Aimag(mq2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+286) = Real(ml2(i1,i2), dp) 
g(SumI+287) = Aimag(ml2(i1,i2)) 
End Do 
End Do 

g(304) = mHd2  
g(305) = mHu2  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+306) = Real(md2(i1,i2), dp) 
g(SumI+307) = Aimag(md2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+324) = Real(mu2(i1,i2), dp) 
g(SumI+325) = Aimag(mu2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+342) = Real(me2(i1,i2), dp) 
g(SumI+343) = Aimag(me2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+360) = Real(mv2(i1,i2), dp) 
g(SumI+361) = Aimag(mv2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
g(SumI+378) = mlHd2(i1) 
End Do 

g(381) = Real(M1,dp)  
g(382) = Aimag(M1)  
g(383) = Real(M2,dp)  
g(384) = Aimag(M2)  
g(385) = Real(M3,dp)  
g(386) = Aimag(M3)  
Iname = Iname - 1 
 
End Subroutine ParametersToG386

Subroutine rge386(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2,i3,i4 
Integer :: j1,j2,j3,j4,j5,j6,j7 
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,         & 
& Dg3,mHd2,betamHd21,betamHd22,DmHd2,mHu2,betamHu21,betamHu22,DmHu2,mlHd2(3)             & 
& ,betamlHd21(3),betamlHd22(3),DmlHd2(3)
Complex(dp) :: Yd(3,3),betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3)          & 
& ,betaYe1(3,3),betaYe2(3,3),DYe(3,3),adjYe(3,3),lam(3),betalam1(3),betalam2(3)          & 
& ,Dlam(3),Yv(3,3),betaYv1(3,3),betaYv2(3,3),DYv(3,3),adjYv(3,3),Yu(3,3),betaYu1(3,3)    & 
& ,betaYu2(3,3),DYu(3,3),adjYu(3,3),kap(3,3,3),betakap1(3,3,3),betakap2(3,3,3)           & 
& ,Dkap(3,3,3),Td(3,3),betaTd1(3,3),betaTd2(3,3),DTd(3,3),adjTd(3,3),Te(3,3)             & 
& ,betaTe1(3,3),betaTe2(3,3),DTe(3,3),adjTe(3,3),Tlam(3),betaTlam1(3),betaTlam2(3)       & 
& ,DTlam(3),Tv(3,3),betaTv1(3,3),betaTv2(3,3),DTv(3,3),adjTv(3,3),Tu(3,3),               & 
& betaTu1(3,3),betaTu2(3,3),DTu(3,3),adjTu(3,3),Tk(3,3,3),betaTk1(3,3,3),betaTk2(3,3,3)  & 
& ,DTk(3,3,3),mq2(3,3),betamq21(3,3),betamq22(3,3),Dmq2(3,3),ml2(3,3),betaml21(3,3)      & 
& ,betaml22(3,3),Dml2(3,3),md2(3,3),betamd21(3,3),betamd22(3,3),Dmd2(3,3),               & 
& mu2(3,3),betamu21(3,3),betamu22(3,3),Dmu2(3,3),me2(3,3),betame21(3,3),betame22(3,3)    & 
& ,Dme2(3,3),mv2(3,3),betamv21(3,3),betamv22(3,3),Dmv2(3,3),M1,betaM11,betaM12,          & 
& DM1,M2,betaM21,betaM22,DM2,M3,betaM31,betaM32,DM3
Real(dp) :: Tr1(3),Tr2(3),Tr3(3) 
Real(dp) :: Tr2U1(3,3) 
Real(dp) :: AbsM1,AbsM2,AbsM3
Complex(dp) :: md2Yd(3,3),me2Ye(3,3),ml2adjYe(3,3),ml2CYv(3,3),mq2adjYd(3,3),mq2adjYu(3,3),          & 
& mu2Yu(3,3),mv2lam(3),mv2TpYv(3,3),Ydmq2(3,3),YdadjYd(3,3),Yeml2(3,3),YeadjYe(3,3),     & 
& Yumq2(3,3),YuadjYu(3,3),YvadjYv(3,3),YvClam(3),adjYdmd2(3,3),adjYdYd(3,3),             & 
& adjYdTd(3,3),adjYeme2(3,3),adjYeYe(3,3),adjYeTe(3,3),adjYumu2(3,3),adjYuYu(3,3),       & 
& adjYuTu(3,3),adjYvmlHd2(3),adjYvYv(3,3),adjYvCml2(3,3),adjYvTv(3,3),adjTdTd(3,3),      & 
& adjTeTe(3,3),adjTuTu(3,3),adjkap1mv2(3,3),adjkap1lam(3),adjkap2mv2(3,3),               & 
& adjkap2lam(3),adjkap3mv2(3,3),adjkap3lam(3),Cmv2adjYv(3,3),Cmv2adjkap1(3,3),           & 
& Cmv2adjkap2(3,3),Cmv2adjkap3(3,3),Cmv2Clam(3),CYeYv(3,3),CYeTv(3,3),CYvmv2(3,3),       & 
& CYvlam(3),CYvTlam(3),CYvTpYv(3,3),CYvTpTv(3,3),CTdTpTd(3,3),CTeTpTe(3,3),              & 
& CTuTpTu(3,3),CTvTlam(3),CTvTpTv(3,3),Ckap1Tpkap1(3,3),Ckap1TpTk1(3,3),Ckap2Tpkap2(3,3),& 
& Ckap2TpTk2(3,3),Ckap3Tpkap3(3,3),Ckap3TpTk3(3,3),TdadjTd(3,3),TeadjTe(3,3),            & 
& TuadjTu(3,3),TvClam(3),TpYvml2(3,3),TpYvCYv(3,3),TpTvCTv(3,3),kap1adjkap1(3,3),        & 
& kap1Ckap1(3,3),kap1Ckap2(3,3),kap1Ckap3(3,3),kap2adjkap2(3,3),kap2Ckap1(3,3),          & 
& kap2Ckap2(3,3),kap2Ckap3(3,3),kap3adjkap3(3,3),kap3Ckap1(3,3),kap3Ckap2(3,3),          & 
& kap3Ckap3(3,3),Tk1adjTk1(3,3),Tk2adjTk2(3,3),Tk3adjTk3(3,3),md2YdadjYd(3,3),           & 
& me2YeadjYe(3,3),ml2adjYeYe(3,3),ml2CYvlam(3),ml2CYvTpYv(3,3),mq2adjYdYd(3,3),          & 
& mq2adjYuYu(3,3),mu2YuadjYu(3,3),mv2TpYvCYv(3,3),mv2kap1adjkap1(3,3),mv2kap2adjkap2(3,3),& 
& mv2kap3adjkap3(3,3),Ydmq2adjYd(3,3),YdadjYdmd2(3,3),YdadjYdYd(3,3),YdadjYdTd(3,3),     & 
& YdadjYuYu(3,3),YdadjYuTu(3,3),Yeml2adjYe(3,3),YeadjYeme2(3,3),YeadjYeYe(3,3),          & 
& YeadjYeTe(3,3),YeCYvTpYv(3,3),YeCYvTpTv(3,3),Yumq2adjYu(3,3),YuadjYdYd(3,3),           & 
& YuadjYdTd(3,3),YuadjYumu2(3,3),YuadjYuYu(3,3),YuadjYuTu(3,3),YvadjYvYv(3,3),           & 
& YvadjYvCml2(3,3),YvadjYvTv(3,3),YvCmv2adjYv(3,3),YvCkap1Tpkap1(3,3),YvCkap1TpTk1(3,3), & 
& YvCkap2Tpkap2(3,3),YvCkap2TpTk2(3,3),YvCkap3Tpkap3(3,3),YvCkap3TpTk3(3,3),             & 
& adjYdmd2Yd(3,3),adjYdYdmq2(3,3),adjYeme2Ye(3,3),adjYeYeml2(3,3),adjYumu2Yu(3,3),       & 
& adjYuYumq2(3,3),CYvmv2lam(3),CYvmv2TpYv(3,3),CYvTpYvml2(3,3),TdadjYdYd(3,3),           & 
& TdadjYuYu(3,3),TeadjYeYe(3,3),TeCYvTpYv(3,3),TuadjYdYd(3,3),TuadjYuYu(3,3),            & 
& TvadjYvYv(3,3),TpYeCYeYv(3,3),TpYeCYeTv(3,3),TpYvml2CYv(3,3),TpYvCYvmv2(3,3),          & 
& TpYvCYvlam(3),TpYvCYvTlam(3),TpTeCYeYv(3,3),TpTvCYvlam(3),kap1adjkap1mv2(3,3),         & 
& kap1Cmv2adjkap1(3,3),kap2adjkap2mv2(3,3),kap2Cmv2adjkap2(3,3),kap3adjkap3mv2(3,3),     & 
& kap3Cmv2adjkap3(3,3),Tk1adjkap1lam(3),Tk2adjkap2lam(3),Tk3adjkap3lam(3)

Complex(dp) :: mv2kap1(3,3),mv2kap2(3,3),mv2kap3(3,3),YdadjYu(3,3),YdadjTd(3,3),YdadjTu(3,3),        & 
& YeadjTe(3,3),YeCYv(3,3),YeCTv(3,3),YuadjYd(3,3),YuadjTd(3,3),YuadjTu(3,3),             & 
& Yvadjkap1(3,3),Yvadjkap2(3,3),Yvadjkap3(3,3),YvadjTk1(3,3),YvadjTk2(3,3),              & 
& YvadjTk3(3,3),YvCTlam(3),YvCkap1(3,3),YvCkap2(3,3),YvCkap3(3,3),adjYdCmd2(3,3),        & 
& adjYeCme2(3,3),adjYuCmu2(3,3),adjTdYd(3,3),adjTeYe(3,3),adjTuYu(3,3),adjkap1Tlam(3),   & 
& adjkap1TpYv(3,3),adjkap1Tpkap1(3,3),adjkap1Tpkap2(3,3),adjkap1Tpkap3(3,3),             & 
& adjkap1kap1(3,3),adjkap1kap2(3,3),adjkap1kap3(3,3),adjkap1Tk1(3,3),adjkap1Tk2(3,3),    & 
& adjkap1Tk3(3,3),adjkap2Tlam(3),adjkap2TpYv(3,3),adjkap2Tpkap1(3,3),adjkap2Tpkap2(3,3), & 
& adjkap2Tpkap3(3,3),adjkap2kap1(3,3),adjkap2kap2(3,3),adjkap2kap3(3,3),adjkap2Tk1(3,3), & 
& adjkap2Tk2(3,3),adjkap2Tk3(3,3),adjkap3Tlam(3),adjkap3TpYv(3,3),adjkap3Tpkap1(3,3),    & 
& adjkap3Tpkap2(3,3),adjkap3Tpkap3(3,3),adjkap3kap1(3,3),adjkap3kap2(3,3),               & 
& adjkap3kap3(3,3),adjkap3Tk1(3,3),adjkap3Tk2(3,3),adjkap3Tk3(3,3),adjTk1lam(3),         & 
& adjTk1TpYv(3,3),adjTk1Tk1(3,3),adjTk2lam(3),adjTk2TpYv(3,3),adjTk2Tk2(3,3),            & 
& adjTk3lam(3),adjTk3TpYv(3,3),adjTk3Tk3(3,3),Cme2CYe(3,3),Cml2adjYe(3,3),               & 
& Cmq2adjYd(3,3),Cmq2adjYu(3,3),Cmv2Ckap1(3,3),Cmv2Ckap2(3,3),Cmv2Ckap3(3,3),            & 
& CYeCml2(3,3),CYvTpkap1(3,3),CYvTpkap2(3,3),CYvTpkap3(3,3),CYvTpTk1(3,3),               & 
& CYvTpTk2(3,3),CYvTpTk3(3,3),CTdTpYd(3,3),CTeTv(3,3),CTeTpYe(3,3),CTuTpYu(3,3),         & 
& CTvlam(3),CTvTpYv(3,3),Ckap1mv2(3,3),Ckap1lam(3),Ckap1Tlam(3),Ckap1TpTv(3,3),          & 
& Ckap1Tpkap2(3,3),Ckap1Tpkap3(3,3),Ckap1TpTk2(3,3),Ckap1TpTk3(3,3),Ckap1kap1(3,3),      & 
& Ckap1kap2(3,3),Ckap1kap3(3,3),Ckap1Tk1(3,3),Ckap1Tk2(3,3),Ckap1Tk3(3,3),               & 
& Ckap2mv2(3,3),Ckap2lam(3),Ckap2Tlam(3),Ckap2TpTv(3,3),Ckap2Tpkap1(3,3),Ckap2Tpkap3(3,3),& 
& Ckap2TpTk1(3,3),Ckap2TpTk3(3,3),Ckap2kap1(3,3),Ckap2kap2(3,3),Ckap2kap3(3,3),          & 
& Ckap2Tk1(3,3),Ckap2Tk2(3,3),Ckap2Tk3(3,3),Ckap3mv2(3,3),Ckap3lam(3),Ckap3Tlam(3),      & 
& Ckap3TpTv(3,3),Ckap3Tpkap1(3,3),Ckap3Tpkap2(3,3),Ckap3TpTk1(3,3),Ckap3TpTk2(3,3),      & 
& Ckap3kap1(3,3),Ckap3kap2(3,3),Ckap3kap3(3,3),Ckap3Tk1(3,3),Ckap3Tk2(3,3),              & 
& Ckap3Tk3(3,3),CTk1lam(3),CTk1Tpkap1(3,3),CTk1Tpkap2(3,3),CTk1Tpkap3(3,3),              & 
& CTk2lam(3),CTk2Tpkap1(3,3),CTk2Tpkap2(3,3),CTk2Tpkap3(3,3),CTk3lam(3),CTk3Tpkap1(3,3), & 
& CTk3Tpkap2(3,3),CTk3Tpkap3(3,3),TdadjYd(3,3),TdadjYu(3,3),TdadjTu(3,3),TeadjYe(3,3),   & 
& TeCYv(3,3),TeCTv(3,3),TuadjYd(3,3),TuadjYu(3,3),TuadjTd(3,3),TvadjYv(3,3),             & 
& TvadjTv(3,3),TvCTlam(3),TpYeCYe(3,3),TpYvmlHd2(3),TpYvadjYe(3,3),TpYvadjTe(3,3),       & 
& TpYvCTv(3,3),TpTeCTe(3,3),TpTvadjYe(3,3),TpTvadjTe(3,3),TpTvCYv(3,3),Tpkap1adjYv(3,3), & 
& Tpkap1Clam(3),Tpkap1Ckap1(3,3),Tpkap1Ckap2(3,3),Tpkap1Ckap3(3,3),Tpkap2adjYv(3,3),     & 
& Tpkap2Clam(3),Tpkap2Ckap1(3,3),Tpkap2Ckap2(3,3),Tpkap2Ckap3(3,3),Tpkap3adjYv(3,3),     & 
& Tpkap3Clam(3),Tpkap3Ckap1(3,3),Tpkap3Ckap2(3,3),Tpkap3Ckap3(3,3),TpTk1adjYv(3,3),      & 
& TpTk1Clam(3),TpTk2adjYv(3,3),TpTk2Clam(3),TpTk3adjYv(3,3),TpTk3Clam(3),kap1adjYv(3,3), & 
& kap1adjkap2(3,3),kap1adjkap3(3,3),kap1adjTk1(3,3),kap1adjTk2(3,3),kap1adjTk3(3,3),     & 
& kap1Clam(3),kap2adjYv(3,3),kap2adjkap1(3,3),kap2adjkap3(3,3),kap2adjTk1(3,3)

Complex(dp) :: kap2adjTk2(3,3),kap2adjTk3(3,3),kap2Clam(3),kap3adjYv(3,3),kap3adjkap1(3,3),           & 
& kap3adjkap2(3,3),kap3adjTk1(3,3),kap3adjTk2(3,3),kap3adjTk3(3,3),kap3Clam(3),          & 
& Tk1adjTv(3,3),Tk1Clam(3),Tk2adjTv(3,3),Tk2Clam(3),Tk3adjTv(3,3),Tk3Clam(3),            & 
& md2YdadjYu(3,3),me2YeCYv(3,3),ml2YvadjYv(3,3),mu2YuadjYd(3,3),mv2TpYvadjYe(3,3),       & 
& mv2kap1adjkap2(3,3),mv2kap1adjkap3(3,3),mv2kap1Clam(3),mv2kap1Ckap1(3,3),              & 
& mv2kap1Ckap2(3,3),mv2kap1Ckap3(3,3),mv2kap2adjkap1(3,3),mv2kap2adjkap3(3,3),           & 
& mv2kap2Clam(3),mv2kap2Ckap1(3,3),mv2kap2Ckap2(3,3),mv2kap2Ckap3(3,3),mv2kap3adjkap1(3,3),& 
& mv2kap3adjkap2(3,3),mv2kap3Clam(3),mv2kap3Ckap1(3,3),mv2kap3Ckap2(3,3),mv2kap3Ckap3(3,3),& 
& Ydmq2adjYu(3,3),YdadjYdCmd2(3,3),YdadjYumu2(3,3),YdadjTdTd(3,3),YdCmq2adjYd(3,3),      & 
& Yeml2CYv(3,3),YeadjYeCme2(3,3),YeadjTeTe(3,3),YeCml2adjYe(3,3),YeCYvmv2(3,3),          & 
& YeCYvlam(3),YeCYvTlam(3),YeCTvTlam(3),Yumq2adjYd(3,3),YuadjYdmd2(3,3),YuadjYuCmu2(3,3),& 
& YuadjTuTu(3,3),YuCmq2adjYu(3,3),YvadjYvmlHd2(3),Yvadjkap1mv2(3,3),Yvadjkap1kap1(3,3),  & 
& Yvadjkap1kap2(3,3),Yvadjkap1kap3(3,3),Yvadjkap1Tk1(3,3),Yvadjkap1Tk2(3,3),             & 
& Yvadjkap1Tk3(3,3),Yvadjkap2mv2(3,3),Yvadjkap2kap1(3,3),Yvadjkap2kap2(3,3),             & 
& Yvadjkap2kap3(3,3),Yvadjkap2Tk1(3,3),Yvadjkap2Tk2(3,3),Yvadjkap2Tk3(3,3),              & 
& Yvadjkap3mv2(3,3),Yvadjkap3kap1(3,3),Yvadjkap3kap2(3,3),Yvadjkap3kap3(3,3),            & 
& Yvadjkap3Tk1(3,3),Yvadjkap3Tk2(3,3),Yvadjkap3Tk3(3,3),YvCmv2adjkap1(3,3),              & 
& YvCmv2adjkap2(3,3),YvCmv2adjkap3(3,3),YvCmv2Clam(3),YvCkap1lam(3),YvCkap1kap2(3,3),    & 
& YvCkap1kap3(3,3),YvCkap1Tk1(3,3),YvCkap1Tk2(3,3),YvCkap1Tk3(3,3),YvCkap2lam(3),        & 
& YvCkap2kap1(3,3),YvCkap2kap3(3,3),YvCkap2Tk1(3,3),YvCkap2Tk2(3,3),YvCkap2Tk3(3,3),     & 
& YvCkap3lam(3),YvCkap3kap1(3,3),YvCkap3kap2(3,3),YvCkap3Tk1(3,3),YvCkap3Tk2(3,3),       & 
& YvCkap3Tk3(3,3),adjYdYdadjYd(3,3),adjYdYdadjYu(3,3),adjYdYdadjTd(3,3),adjYdYdadjTu(3,3),& 
& adjYdTdadjYd(3,3),adjYdTdadjYu(3,3),adjYdTdadjTd(3,3),adjYdTdadjTu(3,3),               & 
& adjYeYeadjYe(3,3),adjYeYeadjTe(3,3),adjYeYeCYv(3,3),adjYeYeCTv(3,3),adjYeTeadjYe(3,3), & 
& adjYeTeadjTe(3,3),adjYeTeCYv(3,3),adjYeTeCTv(3,3),adjYuYuadjYd(3,3),adjYuYuadjYu(3,3), & 
& adjYuYuadjTd(3,3),adjYuYuadjTu(3,3),adjYuTuadjYd(3,3),adjYuTuadjYu(3,3),               & 
& adjYuTuadjTd(3,3),adjYuTuadjTu(3,3),adjYvYvadjYv(3,3),adjYvYvadjkap1(3,3),             & 
& adjYvYvadjkap2(3,3),adjYvYvadjkap3(3,3),adjYvYvadjTk1(3,3),adjYvYvadjTk2(3,3),         & 
& adjYvYvadjTk3(3,3),adjYvYvClam(3),adjYvYvCTlam(3),adjYvYvCkap1(3,3),adjYvYvCkap2(3,3), & 
& adjYvYvCkap3(3,3),adjYvTvadjYv(3,3),adjYvTvadjTv(3,3),adjYvTvClam(3),adjYvTvCTlam(3),  & 
& adjYvTpYeCYe(3,3),adjYvTpTeCTe(3,3),adjTdYdadjYd(3,3),adjTdYdadjYu(3,3),               & 
& adjTdTdadjYd(3,3),adjTdTdadjYu(3,3),adjTeYeadjYe(3,3),adjTeYeCYv(3,3),adjTeTeadjYe(3,3),& 
& adjTeTeCYv(3,3),adjTuYuadjYd(3,3),adjTuYuadjYu(3,3),adjTuTuadjYd(3,3),adjTuTuadjYu(3,3),& 
& adjTvYvClam(3),adjTvTvadjYv(3,3),adjTvTvClam(3),adjkap1mv2lam(3),adjkap1mv2kap1(3,3),  & 
& adjkap1mv2kap2(3,3),adjkap1mv2kap3(3,3),adjkap1TpYvml2(3,3),adjkap1kap1Clam(3),        & 
& adjkap1kap2Clam(3),adjkap1kap3Clam(3),adjkap2mv2lam(3),adjkap2mv2kap1(3,3),            & 
& adjkap2mv2kap2(3,3),adjkap2mv2kap3(3,3),adjkap2TpYvml2(3,3),adjkap2kap1Clam(3)

Complex(dp) :: adjkap2kap2Clam(3),adjkap2kap3Clam(3),adjkap3mv2lam(3),adjkap3mv2kap1(3,3),            & 
& adjkap3mv2kap2(3,3),adjkap3mv2kap3(3,3),adjkap3TpYvml2(3,3),adjkap3kap1Clam(3),        & 
& adjkap3kap2Clam(3),adjkap3kap3Clam(3),Cml2YvadjYv(3,3),Cml2Yvadjkap1(3,3),             & 
& Cml2Yvadjkap2(3,3),Cml2Yvadjkap3(3,3),Cml2YvClam(3),Cml2TpYeCYe(3,3),Cmv2Ckap1kap1(3,3),& 
& Cmv2Ckap1kap2(3,3),Cmv2Ckap1kap3(3,3),Cmv2Ckap2kap1(3,3),Cmv2Ckap2kap2(3,3),           & 
& Cmv2Ckap2kap3(3,3),Cmv2Ckap3kap1(3,3),Cmv2Ckap3kap2(3,3),Cmv2Ckap3kap3(3,3),           & 
& CYeYvClam(3),CYeTvClam(3),CYvTpYvadjYe(3,3),CYvTpYvadjTe(3,3),CYvTpYvCYv(3,3),         & 
& CYvTpYvCTv(3,3),CYvTpTvadjTe(3,3),CYvTpTvCTv(3,3),CYvTpkap1Ckap1(3,3),CYvTpkap1Ckap2(3,3),& 
& CYvTpkap1Ckap3(3,3),CYvTpkap2Ckap1(3,3),CYvTpkap2Ckap2(3,3),CYvTpkap2Ckap3(3,3),       & 
& CYvTpkap3Ckap1(3,3),CYvTpkap3Ckap2(3,3),CYvTpkap3Ckap3(3,3),CTvTpYvadjYe(3,3),         & 
& CTvTpYvCYv(3,3),CTvTpTvadjYe(3,3),CTvTpTvCYv(3,3),Ckap1mv2lam(3),Ckap1TpYvCYv(3,3),    & 
& Ckap1TpYvCTv(3,3),Ckap1Tpkap1adjYv(3,3),Ckap1Tpkap1Clam(3),Ckap1TpTk1adjYv(3,3),       & 
& Ckap1kap1adjYv(3,3),Ckap1kap1Ckap1(3,3),Ckap1kap1Ckap2(3,3),Ckap1kap1Ckap3(3,3),       & 
& Ckap1kap2Ckap1(3,3),Ckap1kap2Ckap2(3,3),Ckap1kap2Ckap3(3,3),Ckap1kap3Ckap1(3,3),       & 
& Ckap1kap3Ckap2(3,3),Ckap1kap3Ckap3(3,3),Ckap1Tk1adjTv(3,3),Ckap2mv2lam(3),             & 
& Ckap2TpYvCYv(3,3),Ckap2TpYvCTv(3,3),Ckap2Tpkap2adjYv(3,3),Ckap2Tpkap2Clam(3),          & 
& Ckap2TpTk2adjYv(3,3),Ckap2kap1Ckap1(3,3),Ckap2kap1Ckap2(3,3),Ckap2kap1Ckap3(3,3),      & 
& Ckap2kap2adjYv(3,3),Ckap2kap2Ckap1(3,3),Ckap2kap2Ckap2(3,3),Ckap2kap2Ckap3(3,3),       & 
& Ckap2kap3Ckap1(3,3),Ckap2kap3Ckap2(3,3),Ckap2kap3Ckap3(3,3),Ckap2Tk2adjTv(3,3),        & 
& Ckap3mv2lam(3),Ckap3TpYvCYv(3,3),Ckap3TpYvCTv(3,3),Ckap3Tpkap3adjYv(3,3),              & 
& Ckap3Tpkap3Clam(3),Ckap3TpTk3adjYv(3,3),Ckap3kap1Ckap1(3,3),Ckap3kap1Ckap2(3,3),       & 
& Ckap3kap1Ckap3(3,3),Ckap3kap2Ckap1(3,3),Ckap3kap2Ckap2(3,3),Ckap3kap2Ckap3(3,3),       & 
& Ckap3kap3adjYv(3,3),Ckap3kap3Ckap1(3,3),Ckap3kap3Ckap2(3,3),Ckap3kap3Ckap3(3,3),       & 
& Ckap3Tk3adjTv(3,3),CTk1TpTk1adjYv(3,3),CTk2TpTk2adjYv(3,3),CTk3TpTk3adjYv(3,3),        & 
& TdadjTdYd(3,3),TeadjTeYe(3,3),TeCYvlam(3),TeCTvlam(3),TeCTvTpYv(3,3),TuadjTuYu(3,3),   & 
& Tvadjkap1kap1(3,3),Tvadjkap1kap2(3,3),Tvadjkap1kap3(3,3),Tvadjkap2kap1(3,3),           & 
& Tvadjkap2kap2(3,3),Tvadjkap2kap3(3,3),Tvadjkap3kap1(3,3),Tvadjkap3kap2(3,3),           & 
& Tvadjkap3kap3(3,3),TvCkap1kap1(3,3),TvCkap1kap2(3,3),TvCkap1kap3(3,3),TvCkap2kap1(3,3),& 
& TvCkap2kap2(3,3),TvCkap2kap3(3,3),TvCkap3kap1(3,3),TvCkap3kap2(3,3),TvCkap3kap3(3,3),  & 
& TpYeCme2CYe(3,3),TpYeCYeCml2(3,3),TpYeCTeTv(3,3),TpYvml2adjYe(3,3),TpYvadjYeme2(3,3),  & 
& TpYvadjYeYe(3,3),TpYvadjYeTe(3,3),TpYvCYvTpYv(3,3),TpYvCYvTpTv(3,3),TpYvCYvTpkap1(3,3),& 
& TpYvCYvTpkap2(3,3),TpYvCYvTpkap3(3,3),TpYvCYvTpTk1(3,3),TpYvCYvTpTk2(3,3),             & 
& TpYvCYvTpTk3(3,3),TpYvCTvTlam(3),TpYvCTvTpTv(3,3),TpTvadjYeYe(3,3),TpTvCYvTpYv(3,3),   & 
& TpTvCTvlam(3),TpTvCTvTpYv(3,3),Tpkap1adjYvYv(3,3),Tpkap1adjkap1lam(3),Tpkap1adjkap1Tlam(3),& 
& Tpkap1adjkap2lam(3),Tpkap1adjkap2Tlam(3),Tpkap1adjkap3lam(3),Tpkap1adjkap3Tlam(3),     & 
& Tpkap1Cmv2adjYv(3,3),Tpkap1Cmv2adjkap1(3,3),Tpkap1Cmv2adjkap2(3,3),Tpkap1Cmv2adjkap3(3,3),& 
& Tpkap1Cmv2Clam(3),Tpkap1Ckap1lam(3),Tpkap1Ckap1Tlam(3),Tpkap1Ckap1TpTv(3,3)

Complex(dp) :: Tpkap1Ckap1TpTk1(3,3),Tpkap1Ckap1TpTk2(3,3),Tpkap1Ckap1TpTk3(3,3),Tpkap1Ckap2lam(3),   & 
& Tpkap1Ckap2Tlam(3),Tpkap1Ckap2TpTv(3,3),Tpkap1Ckap2TpTk1(3,3),Tpkap1Ckap2TpTk2(3,3),   & 
& Tpkap1Ckap2TpTk3(3,3),Tpkap1Ckap3lam(3),Tpkap1Ckap3Tlam(3),Tpkap1Ckap3TpTv(3,3),       & 
& Tpkap1Ckap3TpTk1(3,3),Tpkap1Ckap3TpTk2(3,3),Tpkap1Ckap3TpTk3(3,3),Tpkap2adjYvYv(3,3),  & 
& Tpkap2adjkap1lam(3),Tpkap2adjkap1Tlam(3),Tpkap2adjkap2lam(3),Tpkap2adjkap2Tlam(3),     & 
& Tpkap2adjkap3lam(3),Tpkap2adjkap3Tlam(3),Tpkap2Cmv2adjYv(3,3),Tpkap2Cmv2adjkap1(3,3),  & 
& Tpkap2Cmv2adjkap2(3,3),Tpkap2Cmv2adjkap3(3,3),Tpkap2Cmv2Clam(3),Tpkap2Ckap1lam(3),     & 
& Tpkap2Ckap1Tlam(3),Tpkap2Ckap1TpTv(3,3),Tpkap2Ckap1TpTk1(3,3),Tpkap2Ckap1TpTk2(3,3),   & 
& Tpkap2Ckap1TpTk3(3,3),Tpkap2Ckap2lam(3),Tpkap2Ckap2Tlam(3),Tpkap2Ckap2TpTv(3,3),       & 
& Tpkap2Ckap2TpTk1(3,3),Tpkap2Ckap2TpTk2(3,3),Tpkap2Ckap2TpTk3(3,3),Tpkap2Ckap3lam(3),   & 
& Tpkap2Ckap3Tlam(3),Tpkap2Ckap3TpTv(3,3),Tpkap2Ckap3TpTk1(3,3),Tpkap2Ckap3TpTk2(3,3),   & 
& Tpkap2Ckap3TpTk3(3,3),Tpkap3adjYvYv(3,3),Tpkap3adjkap1lam(3),Tpkap3adjkap1Tlam(3),     & 
& Tpkap3adjkap2lam(3),Tpkap3adjkap2Tlam(3),Tpkap3adjkap3lam(3),Tpkap3adjkap3Tlam(3),     & 
& Tpkap3Cmv2adjYv(3,3),Tpkap3Cmv2adjkap1(3,3),Tpkap3Cmv2adjkap2(3,3),Tpkap3Cmv2adjkap3(3,3),& 
& Tpkap3Cmv2Clam(3),Tpkap3Ckap1lam(3),Tpkap3Ckap1Tlam(3),Tpkap3Ckap1TpTv(3,3),           & 
& Tpkap3Ckap1TpTk1(3,3),Tpkap3Ckap1TpTk2(3,3),Tpkap3Ckap1TpTk3(3,3),Tpkap3Ckap2lam(3),   & 
& Tpkap3Ckap2Tlam(3),Tpkap3Ckap2TpTv(3,3),Tpkap3Ckap2TpTk1(3,3),Tpkap3Ckap2TpTk2(3,3),   & 
& Tpkap3Ckap2TpTk3(3,3),Tpkap3Ckap3lam(3),Tpkap3Ckap3Tlam(3),Tpkap3Ckap3TpTv(3,3),       & 
& Tpkap3Ckap3TpTk1(3,3),Tpkap3Ckap3TpTk2(3,3),Tpkap3Ckap3TpTk3(3,3),TpTk1adjYvYv(3,3),   & 
& TpTk1adjkap1lam(3),TpTk1adjkap2lam(3),TpTk1adjkap3lam(3),TpTk1Ckap1lam(3),             & 
& TpTk1Ckap2lam(3),TpTk1Ckap3lam(3),TpTk2adjYvYv(3,3),TpTk2adjkap1lam(3),TpTk2adjkap2lam(3),& 
& TpTk2adjkap3lam(3),TpTk2Ckap1lam(3),TpTk2Ckap2lam(3),TpTk2Ckap3lam(3),TpTk3adjYvYv(3,3),& 
& TpTk3adjkap1lam(3),TpTk3adjkap2lam(3),TpTk3adjkap3lam(3),TpTk3Ckap1lam(3),             & 
& TpTk3Ckap2lam(3),TpTk3Ckap3lam(3),kap1adjYvCml2(3,3),kap1adjYvTv(3,3),kap1adjkap1lam(3),& 
& kap1adjkap1TpYv(3,3),kap1adjkap1Tpkap1(3,3),kap1adjkap1Tpkap2(3,3),kap1adjkap1Tpkap3(3,3),& 
& kap1Cmv2Clam(3),kap1Ckap1lam(3),kap1Ckap1kap1(3,3),kap1Ckap1kap2(3,3),kap1Ckap1kap3(3,3),& 
& kap1Ckap1Tk1(3,3),kap1Ckap1Tk2(3,3),kap1Ckap1Tk3(3,3),kap1Ckap2lam(3),kap1Ckap2kap1(3,3),& 
& kap1Ckap2kap2(3,3),kap1Ckap2kap3(3,3),kap1Ckap2Tk1(3,3),kap1Ckap2Tk2(3,3),             & 
& kap1Ckap2Tk3(3,3),kap1Ckap3lam(3),kap1Ckap3kap1(3,3),kap1Ckap3kap2(3,3),               & 
& kap1Ckap3kap3(3,3),kap1Ckap3Tk1(3,3),kap1Ckap3Tk2(3,3),kap1Ckap3Tk3(3,3),              & 
& kap2adjYvCml2(3,3),kap2adjYvTv(3,3),kap2adjkap2lam(3),kap2adjkap2TpYv(3,3),            & 
& kap2adjkap2Tpkap1(3,3),kap2adjkap2Tpkap2(3,3),kap2adjkap2Tpkap3(3,3),kap2Cmv2Clam(3),  & 
& kap2Ckap1lam(3),kap2Ckap1kap1(3,3),kap2Ckap1kap2(3,3),kap2Ckap1kap3(3,3),              & 
& kap2Ckap1Tk1(3,3),kap2Ckap1Tk2(3,3),kap2Ckap1Tk3(3,3),kap2Ckap2lam(3),kap2Ckap2kap1(3,3),& 
& kap2Ckap2kap2(3,3),kap2Ckap2kap3(3,3),kap2Ckap2Tk1(3,3),kap2Ckap2Tk2(3,3),             & 
& kap2Ckap2Tk3(3,3),kap2Ckap3lam(3),kap2Ckap3kap1(3,3),kap2Ckap3kap2(3,3),               & 
& kap2Ckap3kap3(3,3),kap2Ckap3Tk1(3,3),kap2Ckap3Tk2(3,3),kap2Ckap3Tk3(3,3)

Complex(dp) :: kap3adjYvCml2(3,3),kap3adjYvTv(3,3),kap3adjkap3lam(3),kap3adjkap3TpYv(3,3),            & 
& kap3adjkap3Tpkap1(3,3),kap3adjkap3Tpkap2(3,3),kap3adjkap3Tpkap3(3,3),kap3Cmv2Clam(3),  & 
& kap3Ckap1lam(3),kap3Ckap1kap1(3,3),kap3Ckap1kap2(3,3),kap3Ckap1kap3(3,3),              & 
& kap3Ckap1Tk1(3,3),kap3Ckap1Tk2(3,3),kap3Ckap1Tk3(3,3),kap3Ckap2lam(3),kap3Ckap2kap1(3,3),& 
& kap3Ckap2kap2(3,3),kap3Ckap2kap3(3,3),kap3Ckap2Tk1(3,3),kap3Ckap2Tk2(3,3),             & 
& kap3Ckap2Tk3(3,3),kap3Ckap3lam(3),kap3Ckap3kap1(3,3),kap3Ckap3kap2(3,3),               & 
& kap3Ckap3kap3(3,3),kap3Ckap3Tk1(3,3),kap3Ckap3Tk2(3,3),kap3Ckap3Tk3(3,3),              & 
& Tk1adjkap1TpYv(3,3),Tk1adjTk1lam(3),Tk1adjTk1TpYv(3,3),Tk1Ckap1kap1(3,3),              & 
& Tk1Ckap1kap2(3,3),Tk1Ckap1kap3(3,3),Tk1Ckap2kap1(3,3),Tk1Ckap2kap2(3,3),               & 
& Tk1Ckap2kap3(3,3),Tk1Ckap3kap1(3,3),Tk1Ckap3kap2(3,3),Tk1Ckap3kap3(3,3),               & 
& Tk2adjkap2TpYv(3,3),Tk2adjTk2lam(3),Tk2adjTk2TpYv(3,3),Tk2Ckap1kap1(3,3),              & 
& Tk2Ckap1kap2(3,3),Tk2Ckap1kap3(3,3),Tk2Ckap2kap1(3,3),Tk2Ckap2kap2(3,3),               & 
& Tk2Ckap2kap3(3,3),Tk2Ckap3kap1(3,3),Tk2Ckap3kap2(3,3),Tk2Ckap3kap3(3,3),               & 
& Tk3adjkap3TpYv(3,3),Tk3adjTk3lam(3),Tk3adjTk3TpYv(3,3),Tk3Ckap1kap1(3,3),              & 
& Tk3Ckap1kap2(3,3),Tk3Ckap1kap3(3,3),Tk3Ckap2kap1(3,3),Tk3Ckap2kap2(3,3),               & 
& Tk3Ckap2kap3(3,3),Tk3Ckap3kap1(3,3),Tk3Ckap3kap2(3,3),Tk3Ckap3kap3(3,3),               & 
& md2YdadjYdYd(3,3),me2YeadjYeYe(3,3),me2YeCYvlam(3),ml2adjYeYeadjYe(3,3),               & 
& ml2adjYeYeCYv(3,3),ml2CYvTpYvadjYe(3,3),ml2CYvTpYvCYv(3,3),mq2adjYdYdadjYd(3,3),       & 
& mq2adjYdYdadjYu(3,3),mq2adjYuYuadjYd(3,3),mq2adjYuYuadjYu(3,3),mu2YuadjYuYu(3,3),      & 
& mv2TpYvCYvlam(3),mv2TpYvCYvTpYv(3,3),mv2kap1adjkap1TpYv(3,3),mv2kap2adjkap2TpYv(3,3),  & 
& mv2kap3adjkap3TpYv(3,3),Ydmq2adjYdYd(3,3),YdadjYdmd2Yd(3,3),YdadjYdYdmq2(3,3),         & 
& YdadjYdYdadjYd(3,3),YdadjYdTdadjYd(3,3),YdadjYdTdadjTd(3,3),YdadjYuYuadjYd(3,3),       & 
& YdadjYuTuadjYd(3,3),YdadjYuTuadjTd(3,3),YdadjTdTdadjYd(3,3),YdadjTuTuadjYd(3,3),       & 
& Yeml2adjYeYe(3,3),Yeml2CYvlam(3),YeadjYeme2Ye(3,3),YeadjYeYeml2(3,3),YeadjYeYeadjYe(3,3),& 
& YeadjYeTeadjYe(3,3),YeadjYeTeadjTe(3,3),YeadjTeTeadjYe(3,3),YeCYvmv2lam(3),            & 
& YeCYvTpYvadjYe(3,3),YeCYvTpTvadjTe(3,3),YeCTvTpTvadjYe(3,3),Yumq2adjYuYu(3,3),         & 
& YuadjYdYdadjYu(3,3),YuadjYdTdadjYu(3,3),YuadjYdTdadjTu(3,3),YuadjYumu2Yu(3,3),         & 
& YuadjYuYumq2(3,3),YuadjYuYuadjYu(3,3),YuadjYuTuadjYu(3,3),YuadjYuTuadjTu(3,3),         & 
& YuadjTdTdadjYu(3,3),YuadjTuTuadjYu(3,3),YvadjYvYvadjYv(3,3),YvadjYvYvClam(3),          & 
& YvadjYvTvadjYv(3,3),YvadjYvTvadjTv(3,3),YvadjYvTvClam(3),YvadjYvTpYeCYe(3,3),          & 
& YvadjYvTpTeCTe(3,3),YvadjTvTvadjYv(3,3),Yvadjkap1mv2kap1(3,3),Yvadjkap1mv2kap2(3,3),   & 
& Yvadjkap1mv2kap3(3,3),Yvadjkap2mv2kap1(3,3),Yvadjkap2mv2kap2(3,3),Yvadjkap2mv2kap3(3,3),& 
& Yvadjkap3mv2kap1(3,3),Yvadjkap3mv2kap2(3,3),Yvadjkap3mv2kap3(3,3),YvCmv2Ckap1kap1(3,3),& 
& YvCmv2Ckap1kap2(3,3),YvCmv2Ckap1kap3(3,3),YvCmv2Ckap2kap1(3,3),YvCmv2Ckap2kap2(3,3),   & 
& YvCmv2Ckap2kap3(3,3),YvCmv2Ckap3kap1(3,3),YvCmv2Ckap3kap2(3,3),YvCmv2Ckap3kap3(3,3),   & 
& YvCkap1Tpkap1adjYv(3,3),YvCkap1Tpkap1Clam(3),YvCkap1TpTk1adjYv(3,3),YvCkap1kap1adjYv(3,3),& 
& YvCkap1Tk1adjTv(3,3),YvCkap2Tpkap2adjYv(3,3),YvCkap2Tpkap2Clam(3),YvCkap2TpTk2adjYv(3,3)

Complex(dp) :: YvCkap2kap2adjYv(3,3),YvCkap2Tk2adjTv(3,3),YvCkap3Tpkap3adjYv(3,3),YvCkap3Tpkap3Clam(3),& 
& YvCkap3TpTk3adjYv(3,3),YvCkap3kap3adjYv(3,3),YvCkap3Tk3adjTv(3,3),YvCTk1TpTk1adjYv(3,3),& 
& YvCTk2TpTk2adjYv(3,3),YvCTk3TpTk3adjYv(3,3),adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3), & 
& adjYdYdmq2adjYd(3,3),adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3),adjYdYdadjYdYd(3,3),    & 
& adjYdYdadjYdTd(3,3),adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3),adjYdYdadjYuTu(3,3),      & 
& adjYdYdadjTdTd(3,3),adjYdTdadjYdYd(3,3),adjYdTdadjYuYu(3,3),adjYdTdadjTdYd(3,3),       & 
& adjYeme2YeadjYe(3,3),adjYeme2YeCYv(3,3),adjYeYeml2adjYe(3,3),adjYeYeml2CYv(3,3),       & 
& adjYeYeadjYeme2(3,3),adjYeYeadjYeYe(3,3),adjYeYeadjYeTe(3,3),adjYeYeadjTeTe(3,3),      & 
& adjYeYeCYvmv2(3,3),adjYeYeCYvlam(3),adjYeYeCTvTlam(3),adjYeTeadjYeYe(3,3),             & 
& adjYeTeadjTeYe(3,3),adjYeTeCYvlam(3),adjYeTeCYvTpYv(3,3),adjYeTeCTvlam(3),             & 
& adjYeTeCTvTpYv(3,3),adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),    & 
& adjYuYumq2adjYu(3,3),adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdTd(3,3),     & 
& adjYuYuadjYumu2(3,3),adjYuYuadjYuYu(3,3),adjYuYuadjYuTu(3,3),adjYuYuadjTuTu(3,3),      & 
& adjYuTuadjYdYd(3,3),adjYuTuadjYuYu(3,3),adjYuTuadjTuYu(3,3),adjYvYvadjYvmlHd2(3),      & 
& adjYvYvadjYvYv(3,3),adjYvYvadjYvTv(3,3),adjYvYvadjkap1mv2(3,3),adjYvYvadjkap2mv2(3,3), & 
& adjYvYvadjkap3mv2(3,3),adjYvYvCmv2adjkap1(3,3),adjYvYvCmv2adjkap2(3,3),adjYvYvCmv2adjkap3(3,3),& 
& adjYvYvCmv2Clam(3),adjYvYvCkap1lam(3),adjYvYvCkap1kap2(3,3),adjYvYvCkap1kap3(3,3),     & 
& adjYvYvCkap2lam(3),adjYvYvCkap2kap1(3,3),adjYvYvCkap2kap3(3,3),adjYvYvCkap3lam(3),     & 
& adjYvYvCkap3kap1(3,3),adjYvYvCkap3kap2(3,3),adjYvCml2YvadjYv(3,3),adjYvCml2Yvadjkap1(3,3),& 
& adjYvCml2Yvadjkap2(3,3),adjYvCml2Yvadjkap3(3,3),adjYvCml2YvClam(3),adjYvCml2TpYeCYe(3,3),& 
& adjYvTvadjYvYv(3,3),adjYvTvadjkap1kap1(3,3),adjYvTvadjkap1kap2(3,3),adjYvTvadjkap1kap3(3,3),& 
& adjYvTvadjkap2kap1(3,3),adjYvTvadjkap2kap2(3,3),adjYvTvadjkap2kap3(3,3),               & 
& adjYvTvadjkap3kap1(3,3),adjYvTvadjkap3kap2(3,3),adjYvTvadjkap3kap3(3,3),               & 
& adjYvTpYeCme2CYe(3,3),adjYvTpYeCYeYv(3,3),adjYvTpYeCYeCml2(3,3),adjYvTpYeCYeTv(3,3),   & 
& adjYvTpYeCTeTv(3,3),adjYvTpTeCYeYv(3,3),adjTdYdadjYdTd(3,3),adjTdTdadjYdYd(3,3),       & 
& adjTeYeadjYeTe(3,3),adjTeYeCYvTlam(3),adjTeTeadjYeYe(3,3),adjTeTeCYvlam(3),            & 
& adjTuYuadjYuTu(3,3),adjTuTuadjYuYu(3,3),adjTvYvCkap1Tk2(3,3),adjTvYvCkap1Tk3(3,3),     & 
& adjTvYvCkap2Tk1(3,3),adjTvYvCkap2Tk3(3,3),adjTvYvCkap3Tk1(3,3),adjTvYvCkap3Tk2(3,3),   & 
& adjkap1Tpkap1Cmv2Clam(3),adjkap1Tpkap2Cmv2Clam(3),adjkap1Tpkap3Cmv2Clam(3),            & 
& adjkap2Tpkap1Cmv2Clam(3),adjkap2Tpkap2Cmv2Clam(3),adjkap2Tpkap3Cmv2Clam(3),            & 
& adjkap3Tpkap1Cmv2Clam(3),adjkap3Tpkap2Cmv2Clam(3),adjkap3Tpkap3Cmv2Clam(3),            & 
& Cml2YvCkap1kap2(3,3),Cml2YvCkap1kap3(3,3),Cml2YvCkap2kap1(3,3),Cml2YvCkap2kap3(3,3),   & 
& Cml2YvCkap3kap1(3,3),Cml2YvCkap3kap2(3,3),Cmv2adjYvYvadjYv(3,3),Cmv2adjYvYvadjkap1(3,3),& 
& Cmv2adjYvYvadjkap2(3,3),Cmv2adjYvYvadjkap3(3,3),Cmv2adjYvYvClam(3),Cmv2adjYvTpYeCYe(3,3),& 
& Cmv2Ckap1kap1adjYv(3,3),Cmv2Ckap2kap2adjYv(3,3),Cmv2Ckap3kap3adjYv(3,3),               & 
& CYeTpYeCYeYv(3,3),CYeTpYeCYeTv(3,3),CYeTpTeCYeYv(3,3),CYvmv2TpYvadjYe(3,3),            & 
& CYvmv2TpYvCYv(3,3),CYvmv2kap1Ckap1(3,3),CYvmv2kap1Ckap2(3,3),CYvmv2kap1Ckap3(3,3)

Complex(dp) :: CYvmv2kap2Ckap1(3,3),CYvmv2kap2Ckap2(3,3),CYvmv2kap2Ckap3(3,3),CYvmv2kap3Ckap1(3,3),   & 
& CYvmv2kap3Ckap2(3,3),CYvmv2kap3Ckap3(3,3),CYvTpYvml2adjYe(3,3),CYvTpYvml2CYv(3,3),     & 
& CYvTpYvadjYeme2(3,3),CYvTpYvadjYeYe(3,3),CYvTpYvadjYeTe(3,3),CYvTpYvCYvmv2(3,3),       & 
& CYvTpYvCYvlam(3),CYvTpYvCYvTlam(3),CYvTpYvCYvTpYv(3,3),CYvTpYvCYvTpTv(3,3),            & 
& CYvTpYvCTvTlam(3),CYvTpYvCTvTpTv(3,3),CYvTpTvadjYeYe(3,3),CYvTpTvCYvlam(3),            & 
& CYvTpTvCYvTpYv(3,3),CYvTpTvCTvlam(3),CYvTpTvCTvTpYv(3,3),CYvkap1adjkap1lam(3),         & 
& CYvkap1adjkap1TpYv(3,3),CYvkap2adjkap2lam(3),CYvkap2adjkap2TpYv(3,3),CYvkap3adjkap3lam(3),& 
& CYvkap3adjkap3TpYv(3,3),CYvTk1adjkap1lam(3),CYvTk1adjkap1TpYv(3,3),CYvTk1adjTk1lam(3), & 
& CYvTk1adjTk1TpYv(3,3),CYvTk2adjkap2lam(3),CYvTk2adjkap2TpYv(3,3),CYvTk2adjTk2lam(3),   & 
& CYvTk2adjTk2TpYv(3,3),CYvTk3adjkap3lam(3),CYvTk3adjkap3TpYv(3,3),CYvTk3adjTk3lam(3),   & 
& CYvTk3adjTk3TpYv(3,3),CTvTpYvCYvTlam(3),CTvTpYvCYvTpTv(3,3),CTvTpTvCYvlam(3),          & 
& CTvTpTvCYvTpYv(3,3),Ckap1TpYvCYvTpTk1(3,3),Ckap1Tpkap1adjYvYv(3,3),Ckap1Tpkap1Cmv2adjYv(3,3),& 
& Ckap1Tpkap1Cmv2adjkap1(3,3),Ckap1Tpkap1Cmv2adjkap2(3,3),Ckap1Tpkap1Cmv2adjkap3(3,3),   & 
& Ckap1TpTk1adjYvYv(3,3),Ckap1kap1adjYvCml2(3,3),Ckap2TpYvCYvTpTk2(3,3),Ckap2Tpkap2adjYvYv(3,3),& 
& Ckap2Tpkap2Cmv2adjYv(3,3),Ckap2Tpkap2Cmv2adjkap1(3,3),Ckap2Tpkap2Cmv2adjkap2(3,3),     & 
& Ckap2Tpkap2Cmv2adjkap3(3,3),Ckap2TpTk2adjYvYv(3,3),Ckap2kap2adjYvCml2(3,3),            & 
& Ckap3TpYvCYvTpTk3(3,3),Ckap3Tpkap3adjYvYv(3,3),Ckap3Tpkap3Cmv2adjYv(3,3),              & 
& Ckap3Tpkap3Cmv2adjkap1(3,3),Ckap3Tpkap3Cmv2adjkap2(3,3),Ckap3Tpkap3Cmv2adjkap3(3,3),   & 
& Ckap3TpTk3adjYvYv(3,3),Ckap3kap3adjYvCml2(3,3),TdadjYdYdadjTd(3,3),TdadjYuYuadjTd(3,3),& 
& TdadjTdYdadjYd(3,3),TdadjTuYuadjYd(3,3),TeadjYeYeadjTe(3,3),TeadjTeYeadjYe(3,3),       & 
& TeCYvTpYvadjTe(3,3),TeCTvTpYvadjYe(3,3),TuadjYdYdadjTu(3,3),TuadjYuYuadjTu(3,3),       & 
& TuadjTdYdadjYu(3,3),TuadjTuYuadjYu(3,3),TvadjYvYvClam(3),TpYeCYeYvClam(3),             & 
& TpYeCYeTvClam(3),TpYvml2CYvlam(3),TpYvml2CYvTpYv(3,3),TpYvadjYeYeCYv(3,3),             & 
& TpYvadjYeTeCYv(3,3),TpYvadjYeTeCTv(3,3),TpYvadjTeTeCYv(3,3),TpYvCYvmv2lam(3),          & 
& TpYvCYvmv2TpYv(3,3),TpYvCYvTpYvml2(3,3),TpYvCYvTpYvCYv(3,3),TpYvCYvTpTvCTv(3,3),       & 
& TpYvCTvTpTvCYv(3,3),TpTeCYeYvClam(3),TpTvadjYeYeCTv(3,3),TpTvadjTeYeCYv(3,3),          & 
& TpTvCYvTpYvCYv(3,3),TpTvCYvTpYvCTv(3,3),TpTvCTvTpYvCYv(3,3),Tpkap1adjkap1mv2lam(3),    & 
& Tpkap1adjkap2mv2lam(3),Tpkap1adjkap3mv2lam(3),Tpkap2adjkap1mv2lam(3),Tpkap2adjkap2mv2lam(3),& 
& Tpkap2adjkap3mv2lam(3),Tpkap3adjkap1mv2lam(3),Tpkap3adjkap2mv2lam(3),Tpkap3adjkap3mv2lam(3),& 
& kap1adjYvYvadjkap1(3,3),kap1adjYvYvCkap1(3,3),kap1adjYvYvCkap2(3,3),kap1adjYvYvCkap3(3,3),& 
& kap1adjkap1TpYvml2(3,3),kap1Cmv2Ckap1kap1(3,3),kap1Cmv2Ckap1kap2(3,3),kap1Cmv2Ckap1kap3(3,3),& 
& kap1Cmv2Ckap2kap1(3,3),kap1Cmv2Ckap2kap2(3,3),kap1Cmv2Ckap2kap3(3,3),kap1Cmv2Ckap3kap1(3,3),& 
& kap1Cmv2Ckap3kap2(3,3),kap1Cmv2Ckap3kap3(3,3),kap2adjYvYvadjkap2(3,3),kap2adjYvYvCkap1(3,3),& 
& kap2adjYvYvCkap2(3,3),kap2adjYvYvCkap3(3,3),kap2adjkap2TpYvml2(3,3),kap2Cmv2Ckap1kap1(3,3),& 
& kap2Cmv2Ckap1kap2(3,3),kap2Cmv2Ckap1kap3(3,3),kap2Cmv2Ckap2kap1(3,3),kap2Cmv2Ckap2kap2(3,3),& 
& kap2Cmv2Ckap2kap3(3,3),kap2Cmv2Ckap3kap1(3,3),kap2Cmv2Ckap3kap2(3,3),kap2Cmv2Ckap3kap3(3,3),& 
& kap3adjYvYvadjkap3(3,3),kap3adjYvYvCkap1(3,3),kap3adjYvYvCkap2(3,3),kap3adjYvYvCkap3(3,3)

Complex(dp) :: kap3adjkap3TpYvml2(3,3),kap3Cmv2Ckap1kap1(3,3),kap3Cmv2Ckap1kap2(3,3),kap3Cmv2Ckap1kap3(3,3),& 
& kap3Cmv2Ckap2kap1(3,3),kap3Cmv2Ckap2kap2(3,3),kap3Cmv2Ckap2kap3(3,3),kap3Cmv2Ckap3kap1(3,3),& 
& kap3Cmv2Ckap3kap2(3,3),kap3Cmv2Ckap3kap3(3,3),Tk1adjYvYvadjTk1(3,3),Tk1adjYvYvCkap1(3,3),& 
& Tk2adjYvYvadjTk2(3,3),Tk2adjYvYvCkap2(3,3),Tk3adjYvYvadjTk3(3,3),Tk3adjYvYvCkap3(3,3), & 
& md2YdadjYdYdadjYd(3,3),md2YdadjYuYuadjYd(3,3),me2YeadjYeYeadjYe(3,3),me2YeCYvTpYvadjYe(3,3),& 
& ml2adjYeYeadjYeYe(3,3),ml2adjYeYeCYvlam(3),ml2CYvTpYvCYvlam(3),ml2CYvTpYvCYvTpYv(3,3), & 
& mq2adjYdYdadjYdYd(3,3),mq2adjYdYdadjYuYu(3,3),mq2adjYuYuadjYdYd(3,3),mq2adjYuYuadjYuYu(3,3),& 
& mu2YuadjYdYdadjYu(3,3),mu2YuadjYuYuadjYu(3,3),mv2TpYvadjYeYeCYv(3,3),mv2TpYvCYvTpYvCYv(3,3),& 
& mv2kap1adjYvYvadjkap1(3,3),mv2kap2adjYvYvadjkap2(3,3),mv2kap3adjYvYvadjkap3(3,3),      & 
& Ydmq2adjYdYdadjYd(3,3),Ydmq2adjYuYuadjYd(3,3),YdadjYdmd2YdadjYd(3,3),YdadjYdYdmq2adjYd(3,3),& 
& YdadjYdYdadjYdmd2(3,3),YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdTd(3,3),YdadjYdTdadjYdYd(3,3),& 
& YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),YdadjYuYuadjYdmd2(3,3),YdadjYuYuadjYdYd(3,3),& 
& YdadjYuYuadjYdTd(3,3),YdadjYuYuadjYuYu(3,3),YdadjYuYuadjYuTu(3,3),YdadjYuTuadjYdYd(3,3),& 
& YdadjYuTuadjYuYu(3,3),Yeml2adjYeYeadjYe(3,3),Yeml2CYvTpYvadjYe(3,3),YeadjYeme2YeadjYe(3,3),& 
& YeadjYeYeml2adjYe(3,3),YeadjYeYeadjYeme2(3,3),YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeTe(3,3),& 
& YeadjYeTeadjYeYe(3,3),YeCYvmv2TpYvadjYe(3,3),YeCYvTpYvml2adjYe(3,3),YeCYvTpYvadjYeme2(3,3),& 
& YeCYvTpYvadjYeYe(3,3),YeCYvTpYvadjYeTe(3,3),YeCYvTpYvCYvTpYv(3,3),YeCYvTpYvCYvTpTv(3,3),& 
& YeCYvTpTvadjYeYe(3,3),YeCYvTpTvCYvTpYv(3,3),YeCYvkap1adjkap1TpYv(3,3),YeCYvkap2adjkap2TpYv(3,3),& 
& YeCYvkap3adjkap3TpYv(3,3),YeCYvTk1adjkap1TpYv(3,3),YeCYvTk2adjkap2TpYv(3,3),           & 
& YeCYvTk3adjkap3TpYv(3,3),Yumq2adjYdYdadjYu(3,3),Yumq2adjYuYuadjYu(3,3),YuadjYdmd2YdadjYu(3,3),& 
& YuadjYdYdmq2adjYu(3,3),YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYdTd(3,3),YuadjYdYdadjYumu2(3,3),& 
& YuadjYdYdadjYuYu(3,3),YuadjYdYdadjYuTu(3,3),YuadjYdTdadjYdYd(3,3),YuadjYdTdadjYuYu(3,3),& 
& YuadjYumu2YuadjYu(3,3),YuadjYuYumq2adjYu(3,3),YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),& 
& YuadjYuYuadjYuTu(3,3),YuadjYuTuadjYuYu(3,3),YvadjYvYvadjYvYv(3,3),YvadjYvYvadjYvTv(3,3),& 
& YvadjYvCml2YvadjYv(3,3),YvadjYvCml2TpYeCYe(3,3),YvadjYvTvadjYvYv(3,3),YvadjYvTpYeCme2CYe(3,3),& 
& YvadjYvTpYeCYeYv(3,3),YvadjYvTpYeCYeCml2(3,3),YvadjYvTpYeCYeTv(3,3),YvadjYvTpTeCYeYv(3,3),& 
& YvCmv2adjYvYvadjYv(3,3),YvCmv2adjYvTpYeCYe(3,3),YvCmv2Ckap1kap1adjYv(3,3),             & 
& YvCmv2Ckap2kap2adjYv(3,3),YvCmv2Ckap3kap3adjYv(3,3),YvCkap1TpYvCYvTpTk1(3,3),          & 
& YvCkap1Tpkap1adjYvYv(3,3),YvCkap1Tpkap1Cmv2adjYv(3,3),YvCkap1TpTk1adjYvYv(3,3),        & 
& YvCkap1kap1adjYvCml2(3,3),YvCkap2TpYvCYvTpTk2(3,3),YvCkap2Tpkap2adjYvYv(3,3),          & 
& YvCkap2Tpkap2Cmv2adjYv(3,3),YvCkap2TpTk2adjYvYv(3,3),YvCkap2kap2adjYvCml2(3,3),        & 
& YvCkap3TpYvCYvTpTk3(3,3),YvCkap3Tpkap3adjYvYv(3,3),YvCkap3Tpkap3Cmv2adjYv(3,3),        & 
& YvCkap3TpTk3adjYvYv(3,3),YvCkap3kap3adjYvCml2(3,3),adjYdmd2YdadjYdYd(3,3),             & 
& adjYdYdmq2adjYdYd(3,3),adjYdYdadjYdmd2Yd(3,3),adjYdYdadjYdYdmq2(3,3),adjYeme2YeadjYeYe(3,3),& 
& adjYeme2YeCYvlam(3),adjYeYeml2adjYeYe(3,3),adjYeYeml2CYvlam(3),adjYeYeadjYeme2Ye(3,3), & 
& adjYeYeadjYeYeml2(3,3),adjYeYeCYvmv2lam(3),adjYumu2YuadjYuYu(3,3),adjYuYumq2adjYuYu(3,3),& 
& adjYuYuadjYumu2Yu(3,3),adjYuYuadjYuYumq2(3,3),adjYvYvadjkap1mv2kap2(3,3)

Complex(dp) :: adjYvYvadjkap1mv2kap3(3,3),adjYvYvadjkap2mv2kap1(3,3),adjYvYvadjkap2mv2kap3(3,3),      & 
& adjYvYvadjkap3mv2kap1(3,3),adjYvYvadjkap3mv2kap2(3,3),adjYvYvCmv2Ckap1kap2(3,3),       & 
& adjYvYvCmv2Ckap1kap3(3,3),adjYvYvCmv2Ckap2kap1(3,3),adjYvYvCmv2Ckap2kap3(3,3),         & 
& adjYvYvCmv2Ckap3kap1(3,3),adjYvYvCmv2Ckap3kap2(3,3),adjYvCml2YvCkap1kap2(3,3),         & 
& adjYvCml2YvCkap1kap3(3,3),adjYvCml2YvCkap2kap1(3,3),adjYvCml2YvCkap2kap3(3,3),         & 
& adjYvCml2YvCkap3kap1(3,3),adjYvCml2YvCkap3kap2(3,3),CYvmv2TpYvCYvlam(3),               & 
& CYvmv2TpYvCYvTpYv(3,3),CYvmv2kap1adjkap1TpYv(3,3),CYvmv2kap2adjkap2TpYv(3,3),          & 
& CYvmv2kap3adjkap3TpYv(3,3),CYvTpYvml2CYvlam(3),CYvTpYvml2CYvTpYv(3,3),CYvTpYvCYvmv2lam(3),& 
& CYvTpYvCYvmv2TpYv(3,3),CYvTpYvCYvTpYvml2(3,3),CYvkap1adjkap1TpYvml2(3,3),              & 
& CYvkap2adjkap2TpYvml2(3,3),CYvkap3adjkap3TpYvml2(3,3),TdadjYdYdadjYdYd(3,3),           & 
& TdadjYuYuadjYdYd(3,3),TdadjYuYuadjYuYu(3,3),TeadjYeYeadjYeYe(3,3),TeCYvTpYvadjYeYe(3,3),& 
& TeCYvTpYvCYvTpYv(3,3),TuadjYdYdadjYdYd(3,3),TuadjYdYdadjYuYu(3,3),TuadjYuYuadjYuYu(3,3),& 
& TvadjYvYvadjYvYv(3,3),TvadjYvTpYeCYeYv(3,3),TpYeCYeTpYeCYeYv(3,3),TpYeCYeTpYeCYeTv(3,3),& 
& TpYeCYeTpTeCYeYv(3,3),TpYvml2adjYeYeCYv(3,3),TpYvml2CYvTpYvCYv(3,3),TpYvadjYeme2YeCYv(3,3),& 
& TpYvadjYeYeml2CYv(3,3),TpYvadjYeYeCYvmv2(3,3),TpYvadjYeYeCYvlam(3),TpYvadjYeTeCYvlam(3),& 
& TpYvCYvmv2TpYvCYv(3,3),TpYvCYvTpYvml2CYv(3,3),TpYvCYvTpYvCYvmv2(3,3),TpYvCYvTpYvCYvlam(3),& 
& TpYvCYvTpYvCYvTlam(3),TpYvCYvTpTvCYvlam(3),TpYvCYvkap1adjkap1lam(3),TpYvCYvkap2adjkap2lam(3),& 
& TpYvCYvkap3adjkap3lam(3),TpYvCYvTk1adjkap1lam(3),TpYvCYvTk2adjkap2lam(3),              & 
& TpYvCYvTk3adjkap3lam(3),TpTeCYeTpYeCYeYv(3,3),TpTvadjYeYeCYvlam(3),TpTvCYvTpYvCYvlam(3),& 
& kap1adjYvYvadjkap1mv2(3,3),kap1adjYvYvCmv2adjkap1(3,3),kap1Cmv2adjYvYvadjkap1(3,3),    & 
& kap1Ckap1Tpkap1Cmv2adjkap1(3,3),kap1Ckap2Tpkap2Cmv2adjkap1(3,3),kap1Ckap3Tpkap3Cmv2adjkap1(3,3),& 
& kap2adjYvYvadjkap2mv2(3,3),kap2adjYvYvCmv2adjkap2(3,3),kap2Cmv2adjYvYvadjkap2(3,3),    & 
& kap2Ckap1Tpkap1Cmv2adjkap2(3,3),kap2Ckap2Tpkap2Cmv2adjkap2(3,3),kap2Ckap3Tpkap3Cmv2adjkap2(3,3),& 
& kap3adjYvYvadjkap3mv2(3,3),kap3adjYvYvCmv2adjkap3(3,3),kap3Cmv2adjYvYvadjkap3(3,3),    & 
& kap3Ckap1Tpkap1Cmv2adjkap3(3,3),kap3Ckap2Tpkap2Cmv2adjkap3(3,3),kap3Ckap3Tpkap3Cmv2adjkap3(3,3),& 
& Tk1adjYvYvCkap1lam(3),Tk2adjYvYvCkap2lam(3),Tk3adjYvYvCkap3lam(3)

Complex(dp) :: Trmd2,Trme2,Trml2,Trmq2,Trmu2,TrYdadjYd,TrYeadjYe,TrYuadjYu,TrYvadjYv,TradjYdTd,      & 
& TradjYeTe,TradjYuTu,TradjYvTv,TrCTdTpTd,TrCTeTpTe,TrCTuTpTu,TrCTvTpTv,Trmd2YdadjYd,    & 
& Trme2YeadjYe,Trml2adjYeYe,Trmq2adjYdYd,Trmq2adjYuYu,Trmu2YuadjYu,TrYvadjYvCml2,        & 
& TrYvCmv2adjYv

Complex(dp) :: TradjTk1Tk1,TradjTk2Tk2,TradjTk3Tk3,TrCTdTpYd,TrCTeTpYe,TrCTuTpYu,TrCTvTpYv,          & 
& TrCkap1Tpkap1,TrCkap1Tpkap2,TrCkap1Tpkap3,TrCkap1kap1,TrCkap1kap2,TrCkap1kap3,         & 
& TrCkap1Tk1,TrCkap1Tk2,TrCkap1Tk3,TrCkap2Tpkap1,TrCkap2Tpkap2,TrCkap2Tpkap3,            & 
& TrCkap2kap1,TrCkap2kap2,TrCkap2kap3,TrCkap2Tk1,TrCkap2Tk2,TrCkap2Tk3,TrCkap3Tpkap1,    & 
& TrCkap3Tpkap2,TrCkap3Tpkap3,TrCkap3kap1,TrCkap3kap2,TrCkap3kap3,TrCkap3Tk1,            & 
& TrCkap3Tk2,TrCkap3Tk3,TrCTk1Tpkap1,TrCTk1Tpkap2,TrCTk1Tpkap3,TrCTk2Tpkap1,             & 
& TrCTk2Tpkap2,TrCTk2Tpkap3,TrCTk3Tpkap1,TrCTk3Tpkap2,TrCTk3Tpkap3,Trml2YvadjYv,         & 
& Trmv2kap1adjkap1,Trmv2kap1adjkap2,Trmv2kap1adjkap3,Trmv2kap2adjkap1,Trmv2kap2adjkap2,  & 
& Trmv2kap2adjkap3,Trmv2kap3adjkap1,Trmv2kap3adjkap2,Trmv2kap3adjkap3,TrYdadjYdCmd2,     & 
& TrYdCmq2adjYd,TrYeadjYeCme2,TrYeCml2adjYe,TrYuadjYuCmu2,TrYuCmq2adjYu,TrYdadjYdYdadjYd,& 
& TrYdadjYdTdadjYd,TrYdadjYdTdadjTd,TrYdadjYuYuadjYd,TrYdadjYuTuadjYd,TrYdadjYuTuadjTd,  & 
& TrYdadjTdTdadjYd,TrYdadjTuTuadjYd,TrYeadjYeYeadjYe,TrYeadjYeTeadjYe,TrYeadjYeTeadjTe,  & 
& TrYeadjTeTeadjYe,TrYeCTvTpTvadjYe,TrYuadjYdTdadjYu,TrYuadjYdTdadjTu,TrYuadjYuYuadjYu,  & 
& TrYuadjYuTuadjYu,TrYuadjYuTuadjTu,TrYuadjTdTdadjYu,TrYuadjTuTuadjYu,TrYvadjYvYvadjYv,  & 
& TrYvadjYvTvadjYv,TrYvadjYvTvadjTv,TrYvadjYvTpYeCYe,TrYvadjYvTpTeCTe,TrYvadjTvTvadjYv,  & 
& TrYvCkap1Tpkap1adjYv,TrYvCkap1TpTk1adjYv,TrYvCkap1kap1adjYv,TrYvCkap1Tk1adjTv,         & 
& TrYvCkap2Tpkap2adjYv,TrYvCkap2TpTk2adjYv,TrYvCkap2kap2adjYv,TrYvCkap2Tk2adjTv,         & 
& TrYvCkap3Tpkap3adjYv,TrYvCkap3TpTk3adjYv,TrYvCkap3kap3adjYv,TrYvCkap3Tk3adjTv,         & 
& TrYvCTk1TpTk1adjYv,TrYvCTk2TpTk2adjYv,TrYvCTk3TpTk3adjYv,TradjYeTeCYvTpYv,             & 
& TradjYeTeCTvTpYv,TradjYvTvadjkap1kap1,TradjYvTvadjkap2kap2,TradjYvTvadjkap3kap3,       & 
& TradjYvTpYeCYeTv,TradjYvTpYeCTeTv,Trmd2YdadjYdYdadjYd,Trmd2YdadjYuYuadjYd,             & 
& Trme2YeadjYeYeadjYe,Trml2adjYeYeadjYeYe,Trmq2adjYdYdadjYdYd,Trmq2adjYdYdadjYuYu,       & 
& Trmq2adjYuYuadjYdYd,Trmq2adjYuYuadjYuYu,Trmu2YuadjYdYdadjYu,Trmu2YuadjYuYuadjYu,       & 
& Trmv2kap1adjYvYvadjkap1,Trmv2kap2adjYvYvadjkap2,Trmv2kap3adjYvYvadjkap3,               & 
& TrYvadjYvCml2YvadjYv,TrYvadjYvCml2TpYeCYe,TrYvadjYvTpYeCme2CYe,TrYvadjYvTpYeCYeCml2,   & 
& TrYvCmv2adjYvYvadjYv,TrYvCmv2adjYvTpYeCYe,TrYvCmv2Ckap1kap1adjYv,TrYvCmv2Ckap2kap2adjYv,& 
& TrYvCmv2Ckap3kap3adjYv,TrYvCkap1Tpkap1Cmv2adjYv,TrYvCkap1kap1adjYvCml2,TrYvCkap2Tpkap2Cmv2adjYv,& 
& TrYvCkap2kap2adjYvCml2,TrYvCkap3Tpkap3Cmv2adjYv,TrYvCkap3kap3adjYvCml2

Real(dp) :: SPlamxxClam,SPlamxxadjYvmlHd2,SPlamxxCmv2Clam,SPClamxxTlam,SPTlamxxCTlam

Real(dp) :: SPmlHd2xxadjYeYeCYvlam,SPlamxxCTlam,SPlamxxadjkap1lam,SPlamxxadjkap2lam,              & 
& SPlamxxadjkap3lam,SPlamxxadjYvYvClam,SPlamxxadjTvTvClam,SPlamxxadjYvYvadjYvmlHd2,      & 
& SPlamxxadjYvYvCmv2Clam,SPlamxxCmv2adjYvYvClam,SPClamxxTpYvmlHd2,SPClamxxTpTvCYvlam,    & 
& SPClamxxTpYvml2CYvlam,SPTpkap1Clamxxadjkap1lam,SPTpkap2Clamxxadjkap2lam,               & 
& SPTpkap3Clamxxadjkap3lam,SPTpTk1ClamxxadjTk1lam,SPTpTk2ClamxxadjTk2lam,SPTpTk3ClamxxadjTk3lam,& 
& SPadjYvYvClamxxTlam,SPadjYvYvCTlamxxTlam,SPadjTvYvClamxxTlam,SPTpTvCYvlamxxCTlam

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Complex(dp) :: sqrt3ov5,ooSqrt15

Real(dp) :: g1p4,g1p5,g2p4,g2p5,g3p4,g3p5

Complex(dp) :: SPlamxxClamp2

Complex(dp) :: Dylami1Clami2(3,3),Dymv2lami1Clami2(3,3),DyYvClami1lami2(3,3),Dylami1adjYvmlHd2i2(3,3),& 
& Dylami1Cmv2Clami2(3,3),DyCYvlami1mlHd2i2(3,3),DyTvClami1lami2(3,3),DyTvi11kap1Ckap1i21(3,3),& 
& DyTvi11kap2Ckap1i22(3,3),DyTvi11kap3Ckap1i23(3,3),DyTvi12kap1Ckap2i21(3,3),            & 
& DyTvi12kap2Ckap2i22(3,3),DyTvi12kap3Ckap2i23(3,3),DyTvi13kap1Ckap3i21(3,3),            & 
& DyTvi13kap2Ckap3i22(3,3),DyTvi13kap3Ckap3i23(3,3),DyTlami1CTlami2(3,3),DyYvClami1Tlami2(3,3)

Complex(dp) :: Dylami1CTlami2(3,3),DyCYvlami1YvClami2(3,3),DyCYvTlami1YvCTlami2(3,3),DyCTvTlami1YvClami2(3,3),& 
& DyCTvlami1TvClami2(3,3),DyCYvlami1TvCTlami2(3,3),Dykap1Clami1Cmv2Ckap1i21(3,3),        & 
& Dykap1Clami1Cmv2Ckap2i21(3,3),Dykap1Clami1Cmv2Ckap3i21(3,3),Dykap1Clami1Ckap1mv2i21(3,3),& 
& Dykap1Clami1Ckap2mv2i21(3,3),Dykap1Clami1Ckap3mv2i21(3,3),DyYvCkap1i11kap1Clami2(3,3), & 
& DyYvCkap2i11kap1Clami2(3,3),DyYvCkap3i11kap1Clami2(3,3),Dykap2Clami1Cmv2Ckap1i22(3,3), & 
& Dykap2Clami1Cmv2Ckap2i22(3,3),Dykap2Clami1Cmv2Ckap3i22(3,3),Dykap2Clami1Ckap1mv2i22(3,3),& 
& Dykap2Clami1Ckap2mv2i22(3,3),Dykap2Clami1Ckap3mv2i22(3,3),DyYvCkap1i12kap2Clami2(3,3), & 
& DyYvCkap2i12kap2Clami2(3,3),DyYvCkap3i12kap2Clami2(3,3),Dykap3Clami1Cmv2Ckap1i23(3,3), & 
& Dykap3Clami1Cmv2Ckap2i23(3,3),Dykap3Clami1Cmv2Ckap3i23(3,3),Dykap3Clami1Ckap1mv2i23(3,3),& 
& Dykap3Clami1Ckap2mv2i23(3,3),Dykap3Clami1Ckap3mv2i23(3,3),DyYvCkap1i13kap3Clami2(3,3), & 
& DyYvCkap2i13kap3Clami2(3,3),DyYvCkap3i13kap3Clami2(3,3),DyTk1Clami1CTk1lami2(3,3),     & 
& DyTk2Clami1CTk2lami2(3,3),DyTk3Clami1CTk3lami2(3,3),Dyml2CYvlami1YvClami2(3,3),        & 
& Dymv2kap1Clami1Ckap1lami2(3,3),Dymv2kap2Clami1Ckap2lami2(3,3),Dymv2kap3Clami1Ckap3lami2(3,3),& 
& DyYeCYvlami1YvClami2(3,3),DyYeCYvlami1TvClami2(3,3),DyYeCYvTlami1YvClami2(3,3),        & 
& DyCYvlami1YvadjYvmlHd2i2(3,3),DyYvadjkap1kap1i11adjYvTv1i2(3,3),DyCYvi11Yvadjkap1kap1i21(3,3),& 
& Dyml2CYvi11Yvadjkap1kap1i21(3,3),DyTeCYvi11Yvadjkap1kap1i21(3,3),DyYvadjkap1kap2i11adjYvTv2i2(3,3),& 
& DyCYvi12Yvadjkap1kap2i21(3,3),Dyml2CYvi12Yvadjkap1kap2i21(3,3),DyTeCYvi12Yvadjkap1kap2i21(3,3),& 
& DyYvadjkap1kap3i11adjYvTv3i2(3,3),DyCYvi13Yvadjkap1kap3i21(3,3),Dyml2CYvi13Yvadjkap1kap3i21(3,3),& 
& DyTeCYvi13Yvadjkap1kap3i21(3,3),DyYvadjkap1Tk1i11lami2(3,3),DyYvadjkap1Tk2i11lami2(3,3),& 
& DyYvadjkap1Tk3i11lami2(3,3),DyYvadjkap2kap1i12adjYvTv1i2(3,3),DyCYvi11Yvadjkap2kap1i22(3,3),& 
& Dyml2CYvi11Yvadjkap2kap1i22(3,3),DyTeCYvi11Yvadjkap2kap1i22(3,3),DyYvadjkap2kap2i12adjYvTv2i2(3,3),& 
& DyCYvi12Yvadjkap2kap2i22(3,3),Dyml2CYvi12Yvadjkap2kap2i22(3,3),DyTeCYvi12Yvadjkap2kap2i22(3,3),& 
& DyYvadjkap2kap3i12adjYvTv3i2(3,3),DyCYvi13Yvadjkap2kap3i22(3,3),Dyml2CYvi13Yvadjkap2kap3i22(3,3),& 
& DyTeCYvi13Yvadjkap2kap3i22(3,3),DyYvadjkap2Tk1i12lami2(3,3),DyYvadjkap2Tk2i12lami2(3,3),& 
& DyYvadjkap2Tk3i12lami2(3,3),DyYvadjkap3kap1i13adjYvTv1i2(3,3),DyCYvi11Yvadjkap3kap1i23(3,3),& 
& Dyml2CYvi11Yvadjkap3kap1i23(3,3),DyTeCYvi11Yvadjkap3kap1i23(3,3),DyYvadjkap3kap2i13adjYvTv2i2(3,3),& 
& DyCYvi12Yvadjkap3kap2i23(3,3),Dyml2CYvi12Yvadjkap3kap2i23(3,3),DyTeCYvi12Yvadjkap3kap2i23(3,3),& 
& DyYvadjkap3kap3i13adjYvTv3i2(3,3),DyCYvi13Yvadjkap3kap3i23(3,3),Dyml2CYvi13Yvadjkap3kap3i23(3,3),& 
& DyTeCYvi13Yvadjkap3kap3i23(3,3),DyYvadjkap3Tk1i13lami2(3,3),DyYvadjkap3Tk2i13lami2(3,3),& 
& DyYvadjkap3Tk3i13lami2(3,3),DyCYvlami1YvCmv2Clami2(3,3),DyYvCkap1lami1Tk1Clami2(3,3),  & 
& DyCTvi11YvCkap1Tk1i21(3,3),DyCTvi12YvCkap1Tk2i21(3,3),DyCTvi13YvCkap1Tk3i21(3,3),      & 
& DyYvCkap2lami1Tk2Clami2(3,3),DyCTvi11YvCkap2Tk1i22(3,3),DyCTvi12YvCkap2Tk2i22(3,3),    & 
& DyCTvi13YvCkap2Tk3i22(3,3),DyYvCkap3lami1Tk3Clami2(3,3),DyCTvi11YvCkap3Tk1i23(3,3),    & 
& DyCTvi12YvCkap3Tk2i23(3,3),DyCTvi13YvCkap3Tk3i23(3,3),Dylami1adjYvYvClami2(3,3),       & 
& Dymv2lami1adjYvYvClami2(3,3),Dylami1adjYvTvCTlami2(3,3),Dylami1adjTvTvClami2(3,3),     & 
& DyCYvlami1Cml2YvClami2(3,3),DyCYvmv2lami1YvClami2(3,3),DyCYvTpkap1Ckap1i11mlHd2i2(3,3),& 
& DyCYvTpkap1Ckap2i11mlHd2i2(3,3),DyCYvTpkap1Ckap3i11mlHd2i2(3,3),DyCYvTpkap2Ckap1i12mlHd2i2(3,3)

Complex(dp) :: DyCYvTpkap2Ckap2i12mlHd2i2(3,3),DyCYvTpkap2Ckap3i12mlHd2i2(3,3),DyCYvTpkap3Ckap1i13mlHd2i2(3,3),& 
& DyCYvTpkap3Ckap2i13mlHd2i2(3,3),DyCYvTpkap3Ckap3i13mlHd2i2(3,3),Dykap1Clami1Ckap1mv2lami2(3,3),& 
& Dymv2kap1i11Ckap1kap1Ckap1i21(3,3),Dymv2kap1i12Ckap1kap1Ckap2i21(3,3),Dymv2kap1i13Ckap1kap1Ckap3i21(3,3),& 
& Dymv2kap1i11Ckap1kap2Ckap1i22(3,3),Dymv2kap1i12Ckap1kap2Ckap2i22(3,3),Dymv2kap1i13Ckap1kap2Ckap3i22(3,3),& 
& Dymv2kap1i11Ckap1kap3Ckap1i23(3,3),Dymv2kap1i12Ckap1kap3Ckap2i23(3,3),Dymv2kap1i13Ckap1kap3Ckap3i23(3,3),& 
& Dykap2Clami1Ckap2mv2lami2(3,3),Dymv2kap2i11Ckap2kap1Ckap1i21(3,3),Dymv2kap2i12Ckap2kap1Ckap2i21(3,3),& 
& Dymv2kap2i13Ckap2kap1Ckap3i21(3,3),Dymv2kap2i11Ckap2kap2Ckap1i22(3,3),Dymv2kap2i12Ckap2kap2Ckap2i22(3,3),& 
& Dymv2kap2i13Ckap2kap2Ckap3i22(3,3),Dymv2kap2i11Ckap2kap3Ckap1i23(3,3),Dymv2kap2i12Ckap2kap3Ckap2i23(3,3),& 
& Dymv2kap2i13Ckap2kap3Ckap3i23(3,3),Dykap3Clami1Ckap3mv2lami2(3,3),Dymv2kap3i11Ckap3kap1Ckap1i21(3,3),& 
& Dymv2kap3i12Ckap3kap1Ckap2i21(3,3),Dymv2kap3i13Ckap3kap1Ckap3i21(3,3),Dymv2kap3i11Ckap3kap2Ckap1i22(3,3),& 
& Dymv2kap3i12Ckap3kap2Ckap2i22(3,3),Dymv2kap3i13Ckap3kap2Ckap3i22(3,3),Dymv2kap3i11Ckap3kap3Ckap1i23(3,3),& 
& Dymv2kap3i12Ckap3kap3Ckap2i23(3,3),Dymv2kap3i13Ckap3kap3Ckap3i23(3,3),DyTeCYvlami1YvClami2(3,3),& 
& DyTvCkap1kap1i11TpYvCYvi21(3,3),DyYeCYvi11TvCkap1kap1i21(3,3),DyTvCkap1kap2i11TpYvCYvi22(3,3),& 
& DyYeCYvi12TvCkap1kap2i21(3,3),DyTvCkap1kap3i11TpYvCYvi23(3,3),DyYeCYvi13TvCkap1kap3i21(3,3),& 
& DyTvCkap2kap1i12TpYvCYvi21(3,3),DyYeCYvi11TvCkap2kap1i22(3,3),DyTvCkap2kap2i12TpYvCYvi22(3,3),& 
& DyYeCYvi12TvCkap2kap2i22(3,3),DyTvCkap2kap3i12TpYvCYvi23(3,3),DyYeCYvi13TvCkap2kap3i22(3,3),& 
& DyTvCkap3kap1i13TpYvCYvi21(3,3),DyYeCYvi11TvCkap3kap1i23(3,3),DyTvCkap3kap2i13TpYvCYvi22(3,3),& 
& DyYeCYvi12TvCkap3kap2i23(3,3),DyTvCkap3kap3i13TpYvCYvi23(3,3),DyYeCYvi13TvCkap3kap3i23(3,3),& 
& DyTpYvCYvlami1Clami2(3,3),DyTpYvCYvlami1adjYvmlHd2i2(3,3),DyTpYvCYvlami1Cmv2Clami2(3,3),& 
& DyYvClami1TpYvCYvlami2(3,3),DyTvClami1TpYvCYvlami2(3,3),DyTpYvCYvTlami1CTlami2(3,3),   & 
& DyYvClami1TpYvCYvTlami2(3,3),DyTpYvCYvTpkap11i1Ckap1mv2i21(3,3),DyYvCkap1i11TpYvCYvTpkap11i2(3,3),& 
& DyTpYvCYvTpkap12i1Ckap2mv2i21(3,3),DyYvCkap2i11TpYvCYvTpkap12i2(3,3),DyTpYvCYvTpkap13i1Ckap3mv2i21(3,3),& 
& DyYvCkap3i11TpYvCYvTpkap13i2(3,3),DyTpYvCYvTpkap21i1Ckap1mv2i22(3,3),DyYvCkap1i12TpYvCYvTpkap21i2(3,3),& 
& DyTpYvCYvTpkap22i1Ckap2mv2i22(3,3),DyYvCkap2i12TpYvCYvTpkap22i2(3,3),DyTpYvCYvTpkap23i1Ckap3mv2i22(3,3),& 
& DyYvCkap3i12TpYvCYvTpkap23i2(3,3),DyTpYvCYvTpkap31i1Ckap1mv2i23(3,3),DyYvCkap1i13TpYvCYvTpkap31i2(3,3),& 
& DyTpYvCYvTpkap32i1Ckap2mv2i23(3,3),DyYvCkap2i13TpYvCYvTpkap32i2(3,3),DyTpYvCYvTpkap33i1Ckap3mv2i23(3,3),& 
& DyYvCkap3i13TpYvCYvTpkap33i2(3,3),DyTpYvCTvTlami1Clami2(3,3),DyTpTvCYvlami1CTlami2(3,3),& 
& DyYvClami1TpTvCYvlami2(3,3),DyTpTvCTvlami1Clami2(3,3),DyTpkap1Ckap1TpTv1i1lami2(3,3),  & 
& DyTpkap1Ckap2TpTv2i1lami2(3,3),DyTpkap1Ckap3TpTv3i1lami2(3,3),DyTpkap2Ckap1TpTv1i1lami2(3,3),& 
& DyTpkap2Ckap2TpTv2i1lami2(3,3),DyTpkap2Ckap3TpTv3i1lami2(3,3),DyTpkap3Ckap1TpTv1i1lami2(3,3),& 
& DyTpkap3Ckap2TpTv2i1lami2(3,3),DyTpkap3Ckap3TpTv3i1lami2(3,3),DyYvCkap1i11kap1adjYvTvi21(3,3),& 
& DyYvCkap2i11kap1adjYvTvi22(3,3),DyYvCkap3i11kap1adjYvTvi23(3,3),DyYvCkap1i11kap1adjkap1Tpkap11i2(3,3),& 
& DyYvCkap2i11kap1adjkap1Tpkap12i2(3,3),DyYvCkap3i11kap1adjkap1Tpkap13i2(3,3),           & 
& DyYvCkap1i12kap1adjkap1Tpkap21i2(3,3),DyYvCkap2i12kap1adjkap1Tpkap22i2(3,3),           & 
& DyYvCkap3i12kap1adjkap1Tpkap23i2(3,3),DyYvCkap1i13kap1adjkap1Tpkap31i2(3,3),           & 
& DyYvCkap2i13kap1adjkap1Tpkap32i2(3,3),DyYvCkap3i13kap1adjkap1Tpkap33i2(3,3),           & 
& Dykap1Cmv2Clami1Ckap1i21(3,3),Dykap1Cmv2Clami1Ckap2i21(3,3),Dykap1Cmv2Clami1Ckap3i21(3,3)

Complex(dp) :: Dykap1Ckap1kap1i11Cmv2Ckap1i21(3,3),Dykap1Ckap1kap1i11Ckap1mv2i21(3,3),Dykap1Ckap1kap2i11Cmv2Ckap2i21(3,3),& 
& Dykap1Ckap1kap2i11Ckap2mv2i21(3,3),Dykap1Ckap1kap3i11Cmv2Ckap3i21(3,3),Dykap1Ckap1kap3i11Ckap3mv2i21(3,3),& 
& DyYvCkap1i11kap1Ckap1Tk1i21(3,3),DyYvCkap2i11kap1Ckap1Tk2i21(3,3),DyYvCkap3i11kap1Ckap1Tk3i21(3,3),& 
& Dykap1Ckap2kap1i12Cmv2Ckap1i21(3,3),Dykap1Ckap2kap1i12Ckap1mv2i21(3,3),Dykap1Ckap2kap2i12Cmv2Ckap2i21(3,3),& 
& Dykap1Ckap2kap2i12Ckap2mv2i21(3,3),Dykap1Ckap2kap3i12Cmv2Ckap3i21(3,3),Dykap1Ckap2kap3i12Ckap3mv2i21(3,3),& 
& DyYvCkap1i11kap1Ckap2Tk1i22(3,3),DyYvCkap2i11kap1Ckap2Tk2i22(3,3),DyYvCkap3i11kap1Ckap2Tk3i22(3,3),& 
& Dykap1Ckap3kap1i13Cmv2Ckap1i21(3,3),Dykap1Ckap3kap1i13Ckap1mv2i21(3,3),Dykap1Ckap3kap2i13Cmv2Ckap2i21(3,3),& 
& Dykap1Ckap3kap2i13Ckap2mv2i21(3,3),Dykap1Ckap3kap3i13Cmv2Ckap3i21(3,3),Dykap1Ckap3kap3i13Ckap3mv2i21(3,3),& 
& DyYvCkap1i11kap1Ckap3Tk1i23(3,3),DyYvCkap2i11kap1Ckap3Tk2i23(3,3),DyYvCkap3i11kap1Ckap3Tk3i23(3,3),& 
& DyYvCkap1i12kap2adjYvTvi21(3,3),DyYvCkap2i12kap2adjYvTvi22(3,3),DyYvCkap3i12kap2adjYvTvi23(3,3),& 
& DyYvCkap1i11kap2adjkap2Tpkap11i2(3,3),DyYvCkap2i11kap2adjkap2Tpkap12i2(3,3),           & 
& DyYvCkap3i11kap2adjkap2Tpkap13i2(3,3),DyYvCkap1i12kap2adjkap2Tpkap21i2(3,3),           & 
& DyYvCkap2i12kap2adjkap2Tpkap22i2(3,3),DyYvCkap3i12kap2adjkap2Tpkap23i2(3,3),           & 
& DyYvCkap1i13kap2adjkap2Tpkap31i2(3,3),DyYvCkap2i13kap2adjkap2Tpkap32i2(3,3),           & 
& DyYvCkap3i13kap2adjkap2Tpkap33i2(3,3),Dykap2Cmv2Clami1Ckap1i22(3,3),Dykap2Cmv2Clami1Ckap2i22(3,3),& 
& Dykap2Cmv2Clami1Ckap3i22(3,3),Dykap2Ckap1kap1i11Cmv2Ckap1i22(3,3),Dykap2Ckap1kap1i11Ckap1mv2i22(3,3),& 
& Dykap2Ckap1kap2i11Cmv2Ckap2i22(3,3),Dykap2Ckap1kap2i11Ckap2mv2i22(3,3),Dykap2Ckap1kap3i11Cmv2Ckap3i22(3,3),& 
& Dykap2Ckap1kap3i11Ckap3mv2i22(3,3),DyYvCkap1i12kap2Ckap1Tk1i21(3,3),DyYvCkap2i12kap2Ckap1Tk2i21(3,3),& 
& DyYvCkap3i12kap2Ckap1Tk3i21(3,3),Dykap2Ckap2kap1i12Cmv2Ckap1i22(3,3),Dykap2Ckap2kap1i12Ckap1mv2i22(3,3),& 
& Dykap2Ckap2kap2i12Cmv2Ckap2i22(3,3),Dykap2Ckap2kap2i12Ckap2mv2i22(3,3),Dykap2Ckap2kap3i12Cmv2Ckap3i22(3,3),& 
& Dykap2Ckap2kap3i12Ckap3mv2i22(3,3),DyYvCkap1i12kap2Ckap2Tk1i22(3,3),DyYvCkap2i12kap2Ckap2Tk2i22(3,3),& 
& DyYvCkap3i12kap2Ckap2Tk3i22(3,3),Dykap2Ckap3kap1i13Cmv2Ckap1i22(3,3),Dykap2Ckap3kap1i13Ckap1mv2i22(3,3),& 
& Dykap2Ckap3kap2i13Cmv2Ckap2i22(3,3),Dykap2Ckap3kap2i13Ckap2mv2i22(3,3),Dykap2Ckap3kap3i13Cmv2Ckap3i22(3,3),& 
& Dykap2Ckap3kap3i13Ckap3mv2i22(3,3),DyYvCkap1i12kap2Ckap3Tk1i23(3,3),DyYvCkap2i12kap2Ckap3Tk2i23(3,3),& 
& DyYvCkap3i12kap2Ckap3Tk3i23(3,3),DyYvCkap1i13kap3adjYvTvi21(3,3),DyYvCkap2i13kap3adjYvTvi22(3,3),& 
& DyYvCkap3i13kap3adjYvTvi23(3,3),DyYvCkap1i11kap3adjkap3Tpkap11i2(3,3),DyYvCkap2i11kap3adjkap3Tpkap12i2(3,3),& 
& DyYvCkap3i11kap3adjkap3Tpkap13i2(3,3),DyYvCkap1i12kap3adjkap3Tpkap21i2(3,3),           & 
& DyYvCkap2i12kap3adjkap3Tpkap22i2(3,3),DyYvCkap3i12kap3adjkap3Tpkap23i2(3,3),           & 
& DyYvCkap1i13kap3adjkap3Tpkap31i2(3,3),DyYvCkap2i13kap3adjkap3Tpkap32i2(3,3),           & 
& DyYvCkap3i13kap3adjkap3Tpkap33i2(3,3),Dykap3Cmv2Clami1Ckap1i23(3,3),Dykap3Cmv2Clami1Ckap2i23(3,3),& 
& Dykap3Cmv2Clami1Ckap3i23(3,3),Dykap3Ckap1kap1i11Cmv2Ckap1i23(3,3),Dykap3Ckap1kap1i11Ckap1mv2i23(3,3),& 
& Dykap3Ckap1kap2i11Cmv2Ckap2i23(3,3),Dykap3Ckap1kap2i11Ckap2mv2i23(3,3),Dykap3Ckap1kap3i11Cmv2Ckap3i23(3,3),& 
& Dykap3Ckap1kap3i11Ckap3mv2i23(3,3),DyYvCkap1i13kap3Ckap1Tk1i21(3,3),DyYvCkap2i13kap3Ckap1Tk2i21(3,3),& 
& DyYvCkap3i13kap3Ckap1Tk3i21(3,3),Dykap3Ckap2kap1i12Cmv2Ckap1i23(3,3),Dykap3Ckap2kap1i12Ckap1mv2i23(3,3),& 
& Dykap3Ckap2kap2i12Cmv2Ckap2i23(3,3),Dykap3Ckap2kap2i12Ckap2mv2i23(3,3),Dykap3Ckap2kap3i12Cmv2Ckap3i23(3,3),& 
& Dykap3Ckap2kap3i12Ckap3mv2i23(3,3),DyYvCkap1i13kap3Ckap2Tk1i22(3,3),DyYvCkap2i13kap3Ckap2Tk2i22(3,3),& 
& DyYvCkap3i13kap3Ckap2Tk3i22(3,3),Dykap3Ckap3kap1i13Cmv2Ckap1i23(3,3),Dykap3Ckap3kap1i13Ckap1mv2i23(3,3)

Complex(dp) :: Dykap3Ckap3kap2i13Cmv2Ckap2i23(3,3),Dykap3Ckap3kap2i13Ckap2mv2i23(3,3),Dykap3Ckap3kap3i13Cmv2Ckap3i23(3,3),& 
& Dykap3Ckap3kap3i13Ckap3mv2i23(3,3),DyYvCkap1i13kap3Ckap3Tk1i23(3,3),DyYvCkap2i13kap3Ckap3Tk2i23(3,3),& 
& DyYvCkap3i13kap3Ckap3Tk3i23(3,3),DyYvCkap1i11Tk1Ckap1kap1i21(3,3),DyYvCkap1i12Tk1Ckap1kap2i21(3,3),& 
& DyYvCkap1i13Tk1Ckap1kap3i21(3,3),DyYvCkap1i11Tk1Ckap2kap1i22(3,3),DyYvCkap1i12Tk1Ckap2kap2i22(3,3),& 
& DyYvCkap1i13Tk1Ckap2kap3i22(3,3),DyYvCkap1i11Tk1Ckap3kap1i23(3,3),DyYvCkap1i12Tk1Ckap3kap2i23(3,3),& 
& DyYvCkap1i13Tk1Ckap3kap3i23(3,3),DyYvCkap2i11Tk2Ckap1kap1i21(3,3),DyYvCkap2i12Tk2Ckap1kap2i21(3,3),& 
& DyYvCkap2i13Tk2Ckap1kap3i21(3,3),DyYvCkap2i11Tk2Ckap2kap1i22(3,3),DyYvCkap2i12Tk2Ckap2kap2i22(3,3),& 
& DyYvCkap2i13Tk2Ckap2kap3i22(3,3),DyYvCkap2i11Tk2Ckap3kap1i23(3,3),DyYvCkap2i12Tk2Ckap3kap2i23(3,3),& 
& DyYvCkap2i13Tk2Ckap3kap3i23(3,3),DyYvCkap3i11Tk3Ckap1kap1i21(3,3),DyYvCkap3i12Tk3Ckap1kap2i21(3,3),& 
& DyYvCkap3i13Tk3Ckap1kap3i21(3,3),DyYvCkap3i11Tk3Ckap2kap1i22(3,3),DyYvCkap3i12Tk3Ckap2kap2i22(3,3),& 
& DyYvCkap3i13Tk3Ckap2kap3i22(3,3),DyYvCkap3i11Tk3Ckap3kap1i23(3,3),DyYvCkap3i12Tk3Ckap3kap2i23(3,3),& 
& DyYvCkap3i13Tk3Ckap3kap3i23(3,3),Dymv2TpYvCYvlami1Clami2(3,3),DyYvadjYvYvClami1lami2(3,3),& 
& DyYvadjYvTvClami1lami2(3,3),DyCYvi11Yvadjkap1mv2kap1i21(3,3),DyCYvi12Yvadjkap1mv2kap2i21(3,3),& 
& DyCYvi13Yvadjkap1mv2kap3i21(3,3),DyCYvi11Yvadjkap2mv2kap1i22(3,3),DyCYvi12Yvadjkap2mv2kap2i22(3,3),& 
& DyCYvi13Yvadjkap2mv2kap3i22(3,3),DyCYvi11Yvadjkap3mv2kap1i23(3,3),DyCYvi12Yvadjkap3mv2kap2i23(3,3),& 
& DyCYvi13Yvadjkap3mv2kap3i23(3,3),DyCYvi11YvCmv2Ckap1kap1i21(3,3),DyCYvi12YvCmv2Ckap1kap2i21(3,3),& 
& DyCYvi13YvCmv2Ckap1kap3i21(3,3),DyCYvi11YvCmv2Ckap2kap1i22(3,3),DyCYvi12YvCmv2Ckap2kap2i22(3,3),& 
& DyCYvi13YvCmv2Ckap2kap3i22(3,3),DyCYvi11YvCmv2Ckap3kap1i23(3,3),DyCYvi12YvCmv2Ckap3kap2i23(3,3),& 
& DyCYvi13YvCmv2Ckap3kap3i23(3,3),DyYvCkap1Tpkap1Clami1lami2(3,3),DyYvCkap2Tpkap2Clami1lami2(3,3),& 
& DyYvCkap3Tpkap3Clami1lami2(3,3),DyadjYeYeCYvlami1mlHd2i2(3,3),Dylami1adjYvYvadjYvmlHd2i2(3,3),& 
& Dylami1adjYvYvCmv2Clami2(3,3),Dylami1adjYvCml2YvClami2(3,3),Dylami1Cmv2adjYvYvClami2(3,3),& 
& DyCYvTpYvCYvlami1mlHd2i2(3,3),DyTvadjYvYvClami1lami2(3,3),DyTpYeCYeYvClami1lami2(3,3), & 
& DyTpYeCYeTvClami1lami2(3,3),DyTpYvml2CYvlami1Clami2(3,3),DyTpYvCYvmv2lami1Clami2(3,3), & 
& DyTpTeCYeYvClami1lami2(3,3),Dykap1Cmv2Ckap1kap1i11Ckap1i21(3,3),Dykap1Cmv2Ckap1kap2i11Ckap2i21(3,3),& 
& Dykap1Cmv2Ckap1kap3i11Ckap3i21(3,3),Dykap1Cmv2Ckap2kap1i12Ckap1i21(3,3),               & 
& Dykap1Cmv2Ckap2kap2i12Ckap2i21(3,3),Dykap1Cmv2Ckap2kap3i12Ckap3i21(3,3),               & 
& Dykap1Cmv2Ckap3kap1i13Ckap1i21(3,3),Dykap1Cmv2Ckap3kap2i13Ckap2i21(3,3),               & 
& Dykap1Cmv2Ckap3kap3i13Ckap3i21(3,3),Dykap2Cmv2Ckap1kap1i11Ckap1i22(3,3),               & 
& Dykap2Cmv2Ckap1kap2i11Ckap2i22(3,3),Dykap2Cmv2Ckap1kap3i11Ckap3i22(3,3),               & 
& Dykap2Cmv2Ckap2kap1i12Ckap1i22(3,3),Dykap2Cmv2Ckap2kap2i12Ckap2i22(3,3),               & 
& Dykap2Cmv2Ckap2kap3i12Ckap3i22(3,3),Dykap2Cmv2Ckap3kap1i13Ckap1i22(3,3),               & 
& Dykap2Cmv2Ckap3kap2i13Ckap2i22(3,3),Dykap2Cmv2Ckap3kap3i13Ckap3i22(3,3),               & 
& Dykap3Cmv2Ckap1kap1i11Ckap1i23(3,3),Dykap3Cmv2Ckap1kap2i11Ckap2i23(3,3),               & 
& Dykap3Cmv2Ckap1kap3i11Ckap3i23(3,3),Dykap3Cmv2Ckap2kap1i12Ckap1i23(3,3),               & 
& Dykap3Cmv2Ckap2kap2i12Ckap2i23(3,3),Dykap3Cmv2Ckap2kap3i12Ckap3i23(3,3),               & 
& Dykap3Cmv2Ckap3kap1i13Ckap1i23(3,3),Dykap3Cmv2Ckap3kap2i13Ckap2i23(3,3),               & 
& Dykap3Cmv2Ckap3kap3i13Ckap3i23(3,3),DyTvi11kap1Ckap1i22(3,3),DyTvi11kap1Ckap1i23(3,3)

Complex(dp) :: DyTvi11kap2Ckap1i21(3,3),DyTvi11kap2Ckap1i23(3,3),DyTvi11kap3Ckap1i21(3,3),            & 
& DyTvi11kap3Ckap1i22(3,3),DyTvi11kap1Ckap1lami2(3,3),DyTvi11kap2Ckap1lami2(3,3),        & 
& DyTvi11kap3Ckap1lami2(3,3),DyTvi11kap1adjYvYvCkap1i21(3,3),DyTvi11kap2adjYvYvCkap1i22(3,3),& 
& DyTvi11kap3adjYvYvCkap1i23(3,3),DyTvi12kap1Ckap2i22(3,3),DyTvi12kap1Ckap2i23(3,3),     & 
& DyTvi12kap2Ckap2i21(3,3),DyTvi12kap2Ckap2i23(3,3),DyTvi12kap3Ckap2i21(3,3),            & 
& DyTvi12kap3Ckap2i22(3,3),DyTvi12kap1Ckap2lami2(3,3),DyTvi12kap2Ckap2lami2(3,3),        & 
& DyTvi12kap3Ckap2lami2(3,3),DyTvi12kap1adjYvYvCkap2i21(3,3),DyTvi12kap2adjYvYvCkap2i22(3,3),& 
& DyTvi12kap3adjYvYvCkap2i23(3,3),DyTvi13kap1Ckap3i22(3,3),DyTvi13kap1Ckap3i23(3,3),     & 
& DyTvi13kap2Ckap3i21(3,3),DyTvi13kap2Ckap3i23(3,3),DyTvi13kap3Ckap3i21(3,3),            & 
& DyTvi13kap3Ckap3i22(3,3),DyTvi13kap1Ckap3lami2(3,3),DyTvi13kap2Ckap3lami2(3,3),        & 
& DyTvi13kap3Ckap3lami2(3,3),DyTvi13kap1adjYvYvCkap3i21(3,3),DyTvi13kap2adjYvYvCkap3i22(3,3),& 
& DyTvi13kap3adjYvYvCkap3i23(3,3),DyCYvi11Tvi21(3,3),DyCYvi12Tvi21(3,3),DyCYvi13Tvi21(3,3),& 
& DyCTvi11Tvi21(3,3),DyCTvi12Tvi21(3,3),DyCTvi13Tvi21(3,3),DyCYvi11Tvi22(3,3),           & 
& DyCYvi12Tvi22(3,3),DyCYvi13Tvi22(3,3),DyCTvi11Tvi22(3,3),DyCTvi12Tvi22(3,3),           & 
& DyCTvi13Tvi22(3,3),DyCYvi11Tvi23(3,3),DyCYvi12Tvi23(3,3),DyCYvi13Tvi23(3,3),           & 
& DyCTvi11Tvi23(3,3),DyCTvi12Tvi23(3,3),DyCTvi13Tvi23(3,3),DyTlami1Clami2(3,3),          & 
& DyTlami1adjYvYvCTlami2(3,3),DyTlami1adjTvYvClami2(3,3),DyYvadjkap1kap1i11Tlami2(3,3),  & 
& DyYvadjkap1kap2i11Tlami2(3,3),DyYvadjkap1kap3i11Tlami2(3,3),DyYvadjkap2kap1i12Tlami2(3,3),& 
& DyYvadjkap2kap2i12Tlami2(3,3),DyYvadjkap2kap3i12Tlami2(3,3),DyYvadjkap3kap1i13Tlami2(3,3),& 
& DyYvadjkap3kap2i13Tlami2(3,3),DyYvadjkap3kap3i13Tlami2(3,3),DyYvadjYvYvClami1Tlami2(3,3),& 
& DyTpYeCYeYvClami1Tlami2(3,3),Dykap1i11Ckap1TpYvCYvi21(3,3),Dykap1i11adjYvCml2Yvadjkap11i2(3,3),& 
& Dykap1i12Ckap2TpYvCYvi21(3,3),Dykap1i12adjYvCml2Yvadjkap21i2(3,3),Dykap1i13Ckap3TpYvCYvi21(3,3),& 
& Dykap1i13adjYvCml2Yvadjkap31i2(3,3),Dykap2i11Ckap1TpYvCYvi22(3,3),Dykap2i11adjYvCml2Yvadjkap12i2(3,3),& 
& Dykap2i12Ckap2TpYvCYvi22(3,3),Dykap2i12adjYvCml2Yvadjkap22i2(3,3),Dykap2i13Ckap3TpYvCYvi22(3,3),& 
& Dykap2i13adjYvCml2Yvadjkap32i2(3,3),Dykap3i11Ckap1TpYvCYvi23(3,3),Dykap3i11adjYvCml2Yvadjkap13i2(3,3),& 
& Dykap3i12Ckap2TpYvCYvi23(3,3),Dykap3i12adjYvCml2Yvadjkap23i2(3,3),Dykap3i13Ckap3TpYvCYvi23(3,3),& 
& Dykap3i13adjYvCml2Yvadjkap33i2(3,3),DyTk1i11Ckap1i21(3,3),DyTk1i11Ckap1i22(3,3),       & 
& DyTk1i11Ckap1i23(3,3),DyTk1i11CTk1i21(3,3),DyTk1i11CTk1i22(3,3),DyTk1i11CTk1i23(3,3),  & 
& DyTk1i11Ckap1lami2(3,3),DyTk1i11Ckap1TpYvCTvi21(3,3),DyTk1i12Ckap2i21(3,3),            & 
& DyTk1i12Ckap2i22(3,3),DyTk1i12Ckap2i23(3,3),DyTk1i12CTk1i21(3,3),DyTk1i12CTk1i22(3,3), & 
& DyTk1i12CTk1i23(3,3),DyTk1i12Ckap2lami2(3,3),DyTk1i12Ckap2TpYvCTvi21(3,3),             & 
& DyTk1i13Ckap3i21(3,3),DyTk1i13Ckap3i22(3,3),DyTk1i13Ckap3i23(3,3),DyTk1i13CTk1i21(3,3),& 
& DyTk1i13CTk1i22(3,3),DyTk1i13CTk1i23(3,3),DyTk1i13Ckap3lami2(3,3),DyTk1i13Ckap3TpYvCTvi21(3,3),& 
& DyTk2i11Ckap1i21(3,3),DyTk2i11Ckap1i22(3,3),DyTk2i11Ckap1i23(3,3),DyTk2i11CTk2i21(3,3),& 
& DyTk2i11CTk2i22(3,3),DyTk2i11CTk2i23(3,3),DyTk2i11Ckap1lami2(3,3),DyTk2i11Ckap1TpYvCTvi22(3,3),& 
& DyTk2i12Ckap2i21(3,3),DyTk2i12Ckap2i22(3,3),DyTk2i12Ckap2i23(3,3),DyTk2i12CTk2i21(3,3),& 
& DyTk2i12CTk2i22(3,3),DyTk2i12CTk2i23(3,3),DyTk2i12Ckap2lami2(3,3),DyTk2i12Ckap2TpYvCTvi22(3,3)

Complex(dp) :: DyTk2i13Ckap3i21(3,3),DyTk2i13Ckap3i22(3,3),DyTk2i13Ckap3i23(3,3),DyTk2i13CTk2i21(3,3),& 
& DyTk2i13CTk2i22(3,3),DyTk2i13CTk2i23(3,3),DyTk2i13Ckap3lami2(3,3),DyTk2i13Ckap3TpYvCTvi22(3,3),& 
& DyTk3i11Ckap1i21(3,3),DyTk3i11Ckap1i22(3,3),DyTk3i11Ckap1i23(3,3),DyTk3i11CTk3i21(3,3),& 
& DyTk3i11CTk3i22(3,3),DyTk3i11CTk3i23(3,3),DyTk3i11Ckap1lami2(3,3),DyTk3i11Ckap1TpYvCTvi23(3,3),& 
& DyTk3i12Ckap2i21(3,3),DyTk3i12Ckap2i22(3,3),DyTk3i12Ckap2i23(3,3),DyTk3i12CTk3i21(3,3),& 
& DyTk3i12CTk3i22(3,3),DyTk3i12CTk3i23(3,3),DyTk3i12Ckap2lami2(3,3),DyTk3i12Ckap2TpYvCTvi23(3,3),& 
& DyTk3i13Ckap3i21(3,3),DyTk3i13Ckap3i22(3,3),DyTk3i13Ckap3i23(3,3),DyTk3i13CTk3i21(3,3),& 
& DyTk3i13CTk3i22(3,3),DyTk3i13CTk3i23(3,3),DyTk3i13Ckap3lami2(3,3),DyTk3i13Ckap3TpYvCTvi23(3,3)

Complex(dp) :: kap1(3,3), adjkap1(3,3) 
Complex(dp) :: kap2(3,3), adjkap2(3,3) 
Complex(dp) :: kap3(3,3), adjkap3(3,3) 
Complex(dp) :: Tk1(3,3), adjTk1(3,3) 
Complex(dp) :: Tk2(3,3), adjTk2(3,3) 
Complex(dp) :: Tk3(3,3), adjTk3(3,3) 
Iname = Iname +1 
NameOfUnit(Iname) = 'rge386' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters386(gy,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,               & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3)

kap1=kap(:,:,1) 
Call Adjungate(kap1,adjkap1) 
kap2=kap(:,:,2) 
Call Adjungate(kap2,adjkap2) 
kap3=kap(:,:,3) 
Call Adjungate(kap3,adjkap3) 
Tk1=Tk(:,:,1) 
Call Adjungate(Tk1,adjTk1) 
Tk2=Tk(:,:,2) 
Call Adjungate(Tk2,adjTk2) 
Tk3=Tk(:,:,3) 
Call Adjungate(Tk3,adjTk3) 
AbsM1 = Conjg(M1)*M1
AbsM2 = Conjg(M2)*M2
AbsM3 = Conjg(M3)*M3
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yv,adjYv)
Call Adjungate(Yu,adjYu)
Call Adjungate(Td,adjTd)
Call Adjungate(Te,adjTe)
Call Adjungate(Tv,adjTv)
Call Adjungate(Tu,adjTu)
 md2Yd = Matmul(md2,Yd) 
 me2Ye = Matmul(me2,Ye) 
 ml2adjYe = Matmul(ml2,adjYe) 
 ml2CYv = Matmul(ml2,Conjg(Yv)) 
 mq2adjYd = Matmul(mq2,adjYd) 
 mq2adjYu = Matmul(mq2,adjYu) 
 mu2Yu = Matmul(mu2,Yu) 
 mv2lam = Matmul(mv2,lam) 
 mv2TpYv = Matmul(mv2,Transpose(Yv)) 
 Ydmq2 = Matmul(Yd,mq2) 
 YdadjYd = Matmul(Yd,adjYd) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 Yeml2 = Matmul(Ye,ml2) 
 YeadjYe = Matmul(Ye,adjYe) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 Yumq2 = Matmul(Yu,mq2) 
 YuadjYu = Matmul(Yu,adjYu) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 YvadjYv = Matmul(Yv,adjYv) 
Forall(i2=1:3)  YvadjYv(i2,i2) =  Real(YvadjYv(i2,i2),dp) 
 YvClam = Matmul(Yv,Conjg(lam)) 
 adjYdmd2 = Matmul(adjYd,md2) 
 adjYdYd = Matmul(adjYd,Yd) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdTd = Matmul(adjYd,Td) 
 adjYeme2 = Matmul(adjYe,me2) 
 adjYeYe = Matmul(adjYe,Ye) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYeTe = Matmul(adjYe,Te) 
 adjYumu2 = Matmul(adjYu,mu2) 
 adjYuYu = Matmul(adjYu,Yu) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYuTu = Matmul(adjYu,Tu) 
 adjYvmlHd2 = Matmul(adjYv,mlHd2) 
 adjYvYv = Matmul(adjYv,Yv) 
Forall(i2=1:3)  adjYvYv(i2,i2) =  Real(adjYvYv(i2,i2),dp) 
 adjYvCml2 = Matmul(adjYv,Conjg(ml2)) 
 adjYvTv = Matmul(adjYv,Tv) 
 adjTdTd = Matmul(adjTd,Td) 
 adjTeTe = Matmul(adjTe,Te) 
 adjTuTu = Matmul(adjTu,Tu) 
 adjkap1mv2 = Matmul(adjkap1,mv2) 
 adjkap1lam = Matmul(adjkap1,lam) 
 adjkap2mv2 = Matmul(adjkap2,mv2) 
 adjkap2lam = Matmul(adjkap2,lam) 
 adjkap3mv2 = Matmul(adjkap3,mv2) 
 adjkap3lam = Matmul(adjkap3,lam) 
 Cmv2adjYv = Matmul(Conjg(mv2),adjYv) 
 Cmv2adjkap1 = Matmul(Conjg(mv2),adjkap1) 
 Cmv2adjkap2 = Matmul(Conjg(mv2),adjkap2) 
 Cmv2adjkap3 = Matmul(Conjg(mv2),adjkap3) 
 Cmv2Clam = Matmul(Conjg(mv2),Conjg(lam)) 
 CYeYv = Matmul(Conjg(Ye),Yv) 
 CYeTv = Matmul(Conjg(Ye),Tv) 
 CYvmv2 = Matmul(Conjg(Yv),mv2) 
 CYvlam = Matmul(Conjg(Yv),lam) 
 CYvTlam = Matmul(Conjg(Yv),Tlam) 
 CYvTpYv = Matmul(Conjg(Yv),Transpose(Yv)) 
Forall(i2=1:3)  CYvTpYv(i2,i2) =  Real(CYvTpYv(i2,i2),dp) 
 CYvTpTv = Matmul(Conjg(Yv),Transpose(Tv)) 
 CTdTpTd = Matmul(Conjg(Td),Transpose(Td)) 
Forall(i2=1:3)  CTdTpTd(i2,i2) =  Real(CTdTpTd(i2,i2),dp) 
 CTeTpTe = Matmul(Conjg(Te),Transpose(Te)) 
Forall(i2=1:3)  CTeTpTe(i2,i2) =  Real(CTeTpTe(i2,i2),dp) 
 CTuTpTu = Matmul(Conjg(Tu),Transpose(Tu)) 
Forall(i2=1:3)  CTuTpTu(i2,i2) =  Real(CTuTpTu(i2,i2),dp) 
 CTvTlam = Matmul(Conjg(Tv),Tlam) 
 CTvTpTv = Matmul(Conjg(Tv),Transpose(Tv)) 
Forall(i2=1:3)  CTvTpTv(i2,i2) =  Real(CTvTpTv(i2,i2),dp) 
 Ckap1Tpkap1 = Matmul(Conjg(kap1),Transpose(kap1)) 
Forall(i2=1:3)  Ckap1Tpkap1(i2,i2) =  Real(Ckap1Tpkap1(i2,i2),dp) 
 Ckap1TpTk1 = Matmul(Conjg(kap1),Transpose(Tk1)) 
 Ckap2Tpkap2 = Matmul(Conjg(kap2),Transpose(kap2)) 
Forall(i2=1:3)  Ckap2Tpkap2(i2,i2) =  Real(Ckap2Tpkap2(i2,i2),dp) 
 Ckap2TpTk2 = Matmul(Conjg(kap2),Transpose(Tk2)) 
 Ckap3Tpkap3 = Matmul(Conjg(kap3),Transpose(kap3)) 
Forall(i2=1:3)  Ckap3Tpkap3(i2,i2) =  Real(Ckap3Tpkap3(i2,i2),dp) 
 Ckap3TpTk3 = Matmul(Conjg(kap3),Transpose(Tk3)) 
 TdadjTd = Matmul(Td,adjTd) 
 TeadjTe = Matmul(Te,adjTe) 
 TuadjTu = Matmul(Tu,adjTu) 
 TvClam = Matmul(Tv,Conjg(lam)) 
 TpYvml2 = Matmul(Transpose(Yv),ml2) 
 TpYvCYv = Matmul(Transpose(Yv),Conjg(Yv)) 
Forall(i2=1:3)  TpYvCYv(i2,i2) =  Real(TpYvCYv(i2,i2),dp) 
 TpTvCTv = Matmul(Transpose(Tv),Conjg(Tv)) 
Forall(i2=1:3)  TpTvCTv(i2,i2) =  Real(TpTvCTv(i2,i2),dp) 
 kap1adjkap1 = Matmul(kap1,adjkap1) 
 kap1Ckap1 = Matmul(kap1,Conjg(kap1)) 
 kap1Ckap2 = Matmul(kap1,Conjg(kap2)) 
 kap1Ckap3 = Matmul(kap1,Conjg(kap3)) 
 kap2adjkap2 = Matmul(kap2,adjkap2) 
 kap2Ckap1 = Matmul(kap2,Conjg(kap1)) 
 kap2Ckap2 = Matmul(kap2,Conjg(kap2)) 
 kap2Ckap3 = Matmul(kap2,Conjg(kap3)) 
 kap3adjkap3 = Matmul(kap3,adjkap3) 
 kap3Ckap1 = Matmul(kap3,Conjg(kap1)) 
 kap3Ckap2 = Matmul(kap3,Conjg(kap2)) 
 kap3Ckap3 = Matmul(kap3,Conjg(kap3)) 
 Tk1adjTk1 = Matmul(Tk1,adjTk1) 
 Tk2adjTk2 = Matmul(Tk2,adjTk2) 
 Tk3adjTk3 = Matmul(Tk3,adjTk3) 
 md2YdadjYd = Matmul(md2,YdadjYd) 
 me2YeadjYe = Matmul(me2,YeadjYe) 
 ml2adjYeYe = Matmul(ml2,adjYeYe) 
 ml2CYvlam = Matmul(ml2,CYvlam) 
 ml2CYvTpYv = Matmul(ml2,CYvTpYv) 
 mq2adjYdYd = Matmul(mq2,adjYdYd) 
 mq2adjYuYu = Matmul(mq2,adjYuYu) 
 mu2YuadjYu = Matmul(mu2,YuadjYu) 
 mv2TpYvCYv = Matmul(mv2,TpYvCYv) 
 mv2kap1adjkap1 = Matmul(mv2,kap1adjkap1) 
 mv2kap2adjkap2 = Matmul(mv2,kap2adjkap2) 
 mv2kap3adjkap3 = Matmul(mv2,kap3adjkap3) 
 Ydmq2adjYd = Matmul(Yd,mq2adjYd) 
Forall(i2=1:3)  Ydmq2adjYd(i2,i2) =  Real(Ydmq2adjYd(i2,i2),dp) 
 YdadjYdmd2 = Matmul(Yd,adjYdmd2) 
 YdadjYdYd = Matmul(Yd,adjYdYd) 
 YdadjYdTd = Matmul(Yd,adjYdTd) 
 YdadjYuYu = Matmul(Yd,adjYuYu) 
 YdadjYuTu = Matmul(Yd,adjYuTu) 
 Yeml2adjYe = Matmul(Ye,ml2adjYe) 
Forall(i2=1:3)  Yeml2adjYe(i2,i2) =  Real(Yeml2adjYe(i2,i2),dp) 
 YeadjYeme2 = Matmul(Ye,adjYeme2) 
 YeadjYeYe = Matmul(Ye,adjYeYe) 
 YeadjYeTe = Matmul(Ye,adjYeTe) 
 YeCYvTpYv = Matmul(Ye,CYvTpYv) 
 YeCYvTpTv = Matmul(Ye,CYvTpTv) 
 Yumq2adjYu = Matmul(Yu,mq2adjYu) 
Forall(i2=1:3)  Yumq2adjYu(i2,i2) =  Real(Yumq2adjYu(i2,i2),dp) 
 YuadjYdYd = Matmul(Yu,adjYdYd) 
 YuadjYdTd = Matmul(Yu,adjYdTd) 
 YuadjYumu2 = Matmul(Yu,adjYumu2) 
 YuadjYuYu = Matmul(Yu,adjYuYu) 
 YuadjYuTu = Matmul(Yu,adjYuTu) 
 YvadjYvYv = Matmul(Yv,adjYvYv) 
 YvadjYvCml2 = Matmul(Yv,adjYvCml2) 
 YvadjYvTv = Matmul(Yv,adjYvTv) 
 YvCmv2adjYv = Matmul(Yv,Cmv2adjYv) 
Forall(i2=1:3)  YvCmv2adjYv(i2,i2) =  Real(YvCmv2adjYv(i2,i2),dp) 
 YvCkap1Tpkap1 = Matmul(Yv,Ckap1Tpkap1) 
 YvCkap1TpTk1 = Matmul(Yv,Ckap1TpTk1) 
 YvCkap2Tpkap2 = Matmul(Yv,Ckap2Tpkap2) 
 YvCkap2TpTk2 = Matmul(Yv,Ckap2TpTk2) 
 YvCkap3Tpkap3 = Matmul(Yv,Ckap3Tpkap3) 
 YvCkap3TpTk3 = Matmul(Yv,Ckap3TpTk3) 
 adjYdmd2Yd = Matmul(adjYd,md2Yd) 
Forall(i2=1:3)  adjYdmd2Yd(i2,i2) =  Real(adjYdmd2Yd(i2,i2),dp) 
 adjYdYdmq2 = Matmul(adjYd,Ydmq2) 
 adjYeme2Ye = Matmul(adjYe,me2Ye) 
Forall(i2=1:3)  adjYeme2Ye(i2,i2) =  Real(adjYeme2Ye(i2,i2),dp) 
 adjYeYeml2 = Matmul(adjYe,Yeml2) 
 adjYumu2Yu = Matmul(adjYu,mu2Yu) 
Forall(i2=1:3)  adjYumu2Yu(i2,i2) =  Real(adjYumu2Yu(i2,i2),dp) 
 adjYuYumq2 = Matmul(adjYu,Yumq2) 
 CYvmv2lam = Matmul(Conjg(Yv),mv2lam) 
 CYvmv2TpYv = Matmul(Conjg(Yv),mv2TpYv) 
Forall(i2=1:3)  CYvmv2TpYv(i2,i2) =  Real(CYvmv2TpYv(i2,i2),dp) 
 CYvTpYvml2 = Matmul(Conjg(Yv),TpYvml2) 
 TdadjYdYd = Matmul(Td,adjYdYd) 
 TdadjYuYu = Matmul(Td,adjYuYu) 
 TeadjYeYe = Matmul(Te,adjYeYe) 
 TeCYvTpYv = Matmul(Te,CYvTpYv) 
 TuadjYdYd = Matmul(Tu,adjYdYd) 
 TuadjYuYu = Matmul(Tu,adjYuYu) 
 TvadjYvYv = Matmul(Tv,adjYvYv) 
 TpYeCYeYv = Matmul(Transpose(Ye),CYeYv) 
 TpYeCYeTv = Matmul(Transpose(Ye),CYeTv) 
 TpYvml2CYv = Matmul(Transpose(Yv),ml2CYv) 
Forall(i2=1:3)  TpYvml2CYv(i2,i2) =  Real(TpYvml2CYv(i2,i2),dp) 
 TpYvCYvmv2 = Matmul(Transpose(Yv),CYvmv2) 
 TpYvCYvlam = Matmul(Transpose(Yv),CYvlam) 
 TpYvCYvTlam = Matmul(Transpose(Yv),CYvTlam) 
 TpTeCYeYv = Matmul(Transpose(Te),CYeYv) 
 TpTvCYvlam = Matmul(Transpose(Tv),CYvlam) 
 kap1adjkap1mv2 = Matmul(kap1,adjkap1mv2) 
 kap1Cmv2adjkap1 = Matmul(kap1,Cmv2adjkap1) 
 kap2adjkap2mv2 = Matmul(kap2,adjkap2mv2) 
 kap2Cmv2adjkap2 = Matmul(kap2,Cmv2adjkap2) 
 kap3adjkap3mv2 = Matmul(kap3,adjkap3mv2) 
 kap3Cmv2adjkap3 = Matmul(kap3,Cmv2adjkap3) 
 Tk1adjkap1lam = Matmul(Tk1,adjkap1lam) 
 Tk2adjkap2lam = Matmul(Tk2,adjkap2lam) 
 Tk3adjkap3lam = Matmul(Tk3,adjkap3lam) 
 Trmd2 = Real(cTrace(md2),dp) 
 Trme2 = Real(cTrace(me2),dp) 
 Trml2 = Real(cTrace(ml2),dp) 
 Trmq2 = Real(cTrace(mq2),dp) 
 Trmu2 = Real(cTrace(mu2),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYvadjYv = Real(cTrace(YvadjYv),dp) 
 TradjYdTd = Real(cTrace(adjYdTd),dp) 
 TradjYeTe = Real(cTrace(adjYeTe),dp) 
 TradjYuTu = Real(cTrace(adjYuTu),dp) 
 TradjYvTv = Real(cTrace(adjYvTv),dp) 
 TrCTdTpTd = Real(cTrace(CTdTpTd),dp) 
 TrCTeTpTe = Real(cTrace(CTeTpTe),dp) 
 TrCTuTpTu = Real(cTrace(CTuTpTu),dp) 
 TrCTvTpTv = Real(cTrace(CTvTpTv),dp) 
 Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp) 
 Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp) 
 Trml2adjYeYe = Real(cTrace(ml2adjYeYe),dp) 
 Trmq2adjYdYd = Real(cTrace(mq2adjYdYd),dp) 
 Trmq2adjYuYu = Real(cTrace(mq2adjYuYu),dp) 
 Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp) 
 TrYvadjYvCml2 = Real(cTrace(YvadjYvCml2),dp) 
 TrYvCmv2adjYv = Real(cTrace(YvCmv2adjYv),dp) 
 SPlamxxClam = DOT_PRODUCT(Conjg(lam),Conjg(lam)) 
 SPlamxxadjYvmlHd2 = DOT_PRODUCT(Conjg(lam),adjYvmlHd2) 
 SPlamxxCmv2Clam = DOT_PRODUCT(Conjg(lam),Cmv2Clam) 
 SPClamxxTlam = DOT_PRODUCT(Conjg(Conjg(lam)),Tlam) 
 SPTlamxxCTlam = DOT_PRODUCT(Conjg(Tlam),Conjg(Tlam)) 
 sqrt3ov5 =Sqrt(3._dp/5._dp) 
 ooSqrt15 =1._dp/sqrt(15._dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 
 g1p4 =g1**4 
 g1p5 =g1**5 
 g2p4 =g2**4 
 g2p5 =g2**5 
 g3p4 =g3**4 
 g3p5 =g3**5 
 SPlamxxClamp2 =SPlamxxClam**2 
Do i1=1,3
  Do i2=1,3
Dylami1Clami2(i1,i2) = Conjg(lam(i2))*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2lami1Clami2(i1,i2) = Conjg(lam(i2))*mv2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1lami2(i1,i2) = YvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvmlHd2i2(i1,i2) = adjYvmlHd2(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1Cmv2Clami2(i1,i2) = Cmv2Clam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1mlHd2i2(i1,i2) = CYvlam(i1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvClami1lami2(i1,i2) = TvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1i21(i1,i2) = kap1Ckap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1i22(i1,i2) = kap2Ckap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1i23(i1,i2) = kap3Ckap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2i21(i1,i2) = kap1Ckap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2i22(i1,i2) = kap2Ckap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2i23(i1,i2) = kap3Ckap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3i21(i1,i2) = kap1Ckap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3i22(i1,i2) = kap2Ckap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3i23(i1,i2) = kap3Ckap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1CTlami2(i1,i2) = Conjg(Tlam(i2))*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1Tlami2(i1,i2) = YvClam(i1)*Tlam(i2) 
  End Do 
End Do 


If (TwoLoopRGE) Then 
 mv2kap1 = Matmul(mv2,kap1) 
 mv2kap2 = Matmul(mv2,kap2) 
 mv2kap3 = Matmul(mv2,kap3) 
 YdadjYu = Matmul(Yd,adjYu) 
 YdadjTd = Matmul(Yd,adjTd) 
 YdadjTu = Matmul(Yd,adjTu) 
 YeadjTe = Matmul(Ye,adjTe) 
 YeCYv = Matmul(Ye,Conjg(Yv)) 
 YeCTv = Matmul(Ye,Conjg(Tv)) 
 YuadjYd = Matmul(Yu,adjYd) 
 YuadjTd = Matmul(Yu,adjTd) 
 YuadjTu = Matmul(Yu,adjTu) 
 Yvadjkap1 = Matmul(Yv,adjkap1) 
 Yvadjkap2 = Matmul(Yv,adjkap2) 
 Yvadjkap3 = Matmul(Yv,adjkap3) 
 YvadjTk1 = Matmul(Yv,adjTk1) 
 YvadjTk2 = Matmul(Yv,adjTk2) 
 YvadjTk3 = Matmul(Yv,adjTk3) 
 YvCTlam = Matmul(Yv,Conjg(Tlam)) 
 YvCkap1 = Matmul(Yv,Conjg(kap1)) 
 YvCkap2 = Matmul(Yv,Conjg(kap2)) 
 YvCkap3 = Matmul(Yv,Conjg(kap3)) 
 adjYdCmd2 = Matmul(adjYd,Conjg(md2)) 
 adjYeCme2 = Matmul(adjYe,Conjg(me2)) 
 adjYuCmu2 = Matmul(adjYu,Conjg(mu2)) 
 adjTdYd = Matmul(adjTd,Yd) 
 adjTeYe = Matmul(adjTe,Ye) 
 adjTuYu = Matmul(adjTu,Yu) 
 adjkap1Tlam = Matmul(adjkap1,Tlam) 
 adjkap1TpYv = Matmul(adjkap1,Transpose(Yv)) 
 adjkap1Tpkap1 = Matmul(adjkap1,Transpose(kap1)) 
 adjkap1Tpkap2 = Matmul(adjkap1,Transpose(kap2)) 
 adjkap1Tpkap3 = Matmul(adjkap1,Transpose(kap3)) 
 adjkap1kap1 = Matmul(adjkap1,kap1) 
 adjkap1kap2 = Matmul(adjkap1,kap2) 
 adjkap1kap3 = Matmul(adjkap1,kap3) 
 adjkap1Tk1 = Matmul(adjkap1,Tk1) 
 adjkap1Tk2 = Matmul(adjkap1,Tk2) 
 adjkap1Tk3 = Matmul(adjkap1,Tk3) 
 adjkap2Tlam = Matmul(adjkap2,Tlam) 
 adjkap2TpYv = Matmul(adjkap2,Transpose(Yv)) 
 adjkap2Tpkap1 = Matmul(adjkap2,Transpose(kap1)) 
 adjkap2Tpkap2 = Matmul(adjkap2,Transpose(kap2)) 
 adjkap2Tpkap3 = Matmul(adjkap2,Transpose(kap3)) 
 adjkap2kap1 = Matmul(adjkap2,kap1) 
 adjkap2kap2 = Matmul(adjkap2,kap2) 
 adjkap2kap3 = Matmul(adjkap2,kap3) 
 adjkap2Tk1 = Matmul(adjkap2,Tk1) 
 adjkap2Tk2 = Matmul(adjkap2,Tk2) 
 adjkap2Tk3 = Matmul(adjkap2,Tk3) 
 adjkap3Tlam = Matmul(adjkap3,Tlam) 
 adjkap3TpYv = Matmul(adjkap3,Transpose(Yv)) 
 adjkap3Tpkap1 = Matmul(adjkap3,Transpose(kap1)) 
 adjkap3Tpkap2 = Matmul(adjkap3,Transpose(kap2)) 
 adjkap3Tpkap3 = Matmul(adjkap3,Transpose(kap3)) 
 adjkap3kap1 = Matmul(adjkap3,kap1) 
 adjkap3kap2 = Matmul(adjkap3,kap2) 
 adjkap3kap3 = Matmul(adjkap3,kap3) 
 adjkap3Tk1 = Matmul(adjkap3,Tk1) 
 adjkap3Tk2 = Matmul(adjkap3,Tk2) 
 adjkap3Tk3 = Matmul(adjkap3,Tk3) 
 adjTk1lam = Matmul(adjTk1,lam) 
 adjTk1TpYv = Matmul(adjTk1,Transpose(Yv)) 
 adjTk1Tk1 = Matmul(adjTk1,Tk1) 
 adjTk2lam = Matmul(adjTk2,lam) 
 adjTk2TpYv = Matmul(adjTk2,Transpose(Yv)) 
 adjTk2Tk2 = Matmul(adjTk2,Tk2) 
 adjTk3lam = Matmul(adjTk3,lam) 
 adjTk3TpYv = Matmul(adjTk3,Transpose(Yv)) 
 adjTk3Tk3 = Matmul(adjTk3,Tk3) 
 Cme2CYe = Matmul(Conjg(me2),Conjg(Ye)) 
 Cml2adjYe = Matmul(Conjg(ml2),adjYe) 
 Cmq2adjYd = Matmul(Conjg(mq2),adjYd) 
 Cmq2adjYu = Matmul(Conjg(mq2),adjYu) 
 Cmv2Ckap1 = Matmul(Conjg(mv2),Conjg(kap1)) 
 Cmv2Ckap2 = Matmul(Conjg(mv2),Conjg(kap2)) 
 Cmv2Ckap3 = Matmul(Conjg(mv2),Conjg(kap3)) 
 CYeCml2 = Matmul(Conjg(Ye),Conjg(ml2)) 
 CYvTpkap1 = Matmul(Conjg(Yv),Transpose(kap1)) 
 CYvTpkap2 = Matmul(Conjg(Yv),Transpose(kap2)) 
 CYvTpkap3 = Matmul(Conjg(Yv),Transpose(kap3)) 
 CYvTpTk1 = Matmul(Conjg(Yv),Transpose(Tk1)) 
 CYvTpTk2 = Matmul(Conjg(Yv),Transpose(Tk2)) 
 CYvTpTk3 = Matmul(Conjg(Yv),Transpose(Tk3)) 
 CTdTpYd = Matmul(Conjg(Td),Transpose(Yd)) 
 CTeTv = Matmul(Conjg(Te),Tv) 
 CTeTpYe = Matmul(Conjg(Te),Transpose(Ye)) 
 CTuTpYu = Matmul(Conjg(Tu),Transpose(Yu)) 
 CTvlam = Matmul(Conjg(Tv),lam) 
 CTvTpYv = Matmul(Conjg(Tv),Transpose(Yv)) 
 Ckap1mv2 = Matmul(Conjg(kap1),mv2) 
 Ckap1lam = Matmul(Conjg(kap1),lam) 
 Ckap1Tlam = Matmul(Conjg(kap1),Tlam) 
 Ckap1TpTv = Matmul(Conjg(kap1),Transpose(Tv)) 
 Ckap1Tpkap2 = Matmul(Conjg(kap1),Transpose(kap2)) 
 Ckap1Tpkap3 = Matmul(Conjg(kap1),Transpose(kap3)) 
 Ckap1TpTk2 = Matmul(Conjg(kap1),Transpose(Tk2)) 
 Ckap1TpTk3 = Matmul(Conjg(kap1),Transpose(Tk3)) 
 Ckap1kap1 = Matmul(Conjg(kap1),kap1) 
 Ckap1kap2 = Matmul(Conjg(kap1),kap2) 
 Ckap1kap3 = Matmul(Conjg(kap1),kap3) 
 Ckap1Tk1 = Matmul(Conjg(kap1),Tk1) 
 Ckap1Tk2 = Matmul(Conjg(kap1),Tk2) 
 Ckap1Tk3 = Matmul(Conjg(kap1),Tk3) 
 Ckap2mv2 = Matmul(Conjg(kap2),mv2) 
 Ckap2lam = Matmul(Conjg(kap2),lam) 
 Ckap2Tlam = Matmul(Conjg(kap2),Tlam) 
 Ckap2TpTv = Matmul(Conjg(kap2),Transpose(Tv)) 
 Ckap2Tpkap1 = Matmul(Conjg(kap2),Transpose(kap1)) 
 Ckap2Tpkap3 = Matmul(Conjg(kap2),Transpose(kap3)) 
 Ckap2TpTk1 = Matmul(Conjg(kap2),Transpose(Tk1)) 
 Ckap2TpTk3 = Matmul(Conjg(kap2),Transpose(Tk3)) 
 Ckap2kap1 = Matmul(Conjg(kap2),kap1) 
 Ckap2kap2 = Matmul(Conjg(kap2),kap2) 
 Ckap2kap3 = Matmul(Conjg(kap2),kap3) 
 Ckap2Tk1 = Matmul(Conjg(kap2),Tk1) 
 Ckap2Tk2 = Matmul(Conjg(kap2),Tk2) 
 Ckap2Tk3 = Matmul(Conjg(kap2),Tk3) 
 Ckap3mv2 = Matmul(Conjg(kap3),mv2) 
 Ckap3lam = Matmul(Conjg(kap3),lam) 
 Ckap3Tlam = Matmul(Conjg(kap3),Tlam) 
 Ckap3TpTv = Matmul(Conjg(kap3),Transpose(Tv)) 
 Ckap3Tpkap1 = Matmul(Conjg(kap3),Transpose(kap1)) 
 Ckap3Tpkap2 = Matmul(Conjg(kap3),Transpose(kap2)) 
 Ckap3TpTk1 = Matmul(Conjg(kap3),Transpose(Tk1)) 
 Ckap3TpTk2 = Matmul(Conjg(kap3),Transpose(Tk2)) 
 Ckap3kap1 = Matmul(Conjg(kap3),kap1) 
 Ckap3kap2 = Matmul(Conjg(kap3),kap2) 
 Ckap3kap3 = Matmul(Conjg(kap3),kap3) 
 Ckap3Tk1 = Matmul(Conjg(kap3),Tk1) 
 Ckap3Tk2 = Matmul(Conjg(kap3),Tk2) 
 Ckap3Tk3 = Matmul(Conjg(kap3),Tk3) 
 CTk1lam = Matmul(Conjg(Tk1),lam) 
 CTk1Tpkap1 = Matmul(Conjg(Tk1),Transpose(kap1)) 
 CTk1Tpkap2 = Matmul(Conjg(Tk1),Transpose(kap2)) 
 CTk1Tpkap3 = Matmul(Conjg(Tk1),Transpose(kap3)) 
 CTk2lam = Matmul(Conjg(Tk2),lam) 
 CTk2Tpkap1 = Matmul(Conjg(Tk2),Transpose(kap1)) 
 CTk2Tpkap2 = Matmul(Conjg(Tk2),Transpose(kap2)) 
 CTk2Tpkap3 = Matmul(Conjg(Tk2),Transpose(kap3)) 
 CTk3lam = Matmul(Conjg(Tk3),lam) 
 CTk3Tpkap1 = Matmul(Conjg(Tk3),Transpose(kap1)) 
 CTk3Tpkap2 = Matmul(Conjg(Tk3),Transpose(kap2)) 
 CTk3Tpkap3 = Matmul(Conjg(Tk3),Transpose(kap3)) 
 TdadjYd = Matmul(Td,adjYd) 
 TdadjYu = Matmul(Td,adjYu) 
 TdadjTu = Matmul(Td,adjTu) 
 TeadjYe = Matmul(Te,adjYe) 
 TeCYv = Matmul(Te,Conjg(Yv)) 
 TeCTv = Matmul(Te,Conjg(Tv)) 
 TuadjYd = Matmul(Tu,adjYd) 
 TuadjYu = Matmul(Tu,adjYu) 
 TuadjTd = Matmul(Tu,adjTd) 
 TvadjYv = Matmul(Tv,adjYv) 
 TvadjTv = Matmul(Tv,adjTv) 
 TvCTlam = Matmul(Tv,Conjg(Tlam)) 
 TpYeCYe = Matmul(Transpose(Ye),Conjg(Ye)) 
Forall(i2=1:3)  TpYeCYe(i2,i2) =  Real(TpYeCYe(i2,i2),dp) 
 TpYvmlHd2 = Matmul(Transpose(Yv),mlHd2) 
 TpYvadjYe = Matmul(Transpose(Yv),adjYe) 
 TpYvadjTe = Matmul(Transpose(Yv),adjTe) 
 TpYvCTv = Matmul(Transpose(Yv),Conjg(Tv)) 
 TpTeCTe = Matmul(Transpose(Te),Conjg(Te)) 
Forall(i2=1:3)  TpTeCTe(i2,i2) =  Real(TpTeCTe(i2,i2),dp) 
 TpTvadjYe = Matmul(Transpose(Tv),adjYe) 
 TpTvadjTe = Matmul(Transpose(Tv),adjTe) 
 TpTvCYv = Matmul(Transpose(Tv),Conjg(Yv)) 
 Tpkap1adjYv = Matmul(Transpose(kap1),adjYv) 
 Tpkap1Clam = Matmul(Transpose(kap1),Conjg(lam)) 
 Tpkap1Ckap1 = Matmul(Transpose(kap1),Conjg(kap1)) 
Forall(i2=1:3)  Tpkap1Ckap1(i2,i2) =  Real(Tpkap1Ckap1(i2,i2),dp) 
 Tpkap1Ckap2 = Matmul(Transpose(kap1),Conjg(kap2)) 
 Tpkap1Ckap3 = Matmul(Transpose(kap1),Conjg(kap3)) 
 Tpkap2adjYv = Matmul(Transpose(kap2),adjYv) 
 Tpkap2Clam = Matmul(Transpose(kap2),Conjg(lam)) 
 Tpkap2Ckap1 = Matmul(Transpose(kap2),Conjg(kap1)) 
 Tpkap2Ckap2 = Matmul(Transpose(kap2),Conjg(kap2)) 
Forall(i2=1:3)  Tpkap2Ckap2(i2,i2) =  Real(Tpkap2Ckap2(i2,i2),dp) 
 Tpkap2Ckap3 = Matmul(Transpose(kap2),Conjg(kap3)) 
 Tpkap3adjYv = Matmul(Transpose(kap3),adjYv) 
 Tpkap3Clam = Matmul(Transpose(kap3),Conjg(lam)) 
 Tpkap3Ckap1 = Matmul(Transpose(kap3),Conjg(kap1)) 
 Tpkap3Ckap2 = Matmul(Transpose(kap3),Conjg(kap2)) 
 Tpkap3Ckap3 = Matmul(Transpose(kap3),Conjg(kap3)) 
Forall(i2=1:3)  Tpkap3Ckap3(i2,i2) =  Real(Tpkap3Ckap3(i2,i2),dp) 
 TpTk1adjYv = Matmul(Transpose(Tk1),adjYv) 
 TpTk1Clam = Matmul(Transpose(Tk1),Conjg(lam)) 
 TpTk2adjYv = Matmul(Transpose(Tk2),adjYv) 
 TpTk2Clam = Matmul(Transpose(Tk2),Conjg(lam)) 
 TpTk3adjYv = Matmul(Transpose(Tk3),adjYv) 
 TpTk3Clam = Matmul(Transpose(Tk3),Conjg(lam)) 
 kap1adjYv = Matmul(kap1,adjYv) 
 kap1adjkap2 = Matmul(kap1,adjkap2) 
 kap1adjkap3 = Matmul(kap1,adjkap3) 
 kap1adjTk1 = Matmul(kap1,adjTk1) 
 kap1adjTk2 = Matmul(kap1,adjTk2) 
 kap1adjTk3 = Matmul(kap1,adjTk3) 
 kap1Clam = Matmul(kap1,Conjg(lam)) 
 kap2adjYv = Matmul(kap2,adjYv) 
 kap2adjkap1 = Matmul(kap2,adjkap1) 
 kap2adjkap3 = Matmul(kap2,adjkap3) 
 kap2adjTk1 = Matmul(kap2,adjTk1) 
 kap2adjTk2 = Matmul(kap2,adjTk2) 
 kap2adjTk3 = Matmul(kap2,adjTk3) 
 kap2Clam = Matmul(kap2,Conjg(lam)) 
 kap3adjYv = Matmul(kap3,adjYv) 
 kap3adjkap1 = Matmul(kap3,adjkap1) 
 kap3adjkap2 = Matmul(kap3,adjkap2) 
 kap3adjTk1 = Matmul(kap3,adjTk1) 
 kap3adjTk2 = Matmul(kap3,adjTk2) 
 kap3adjTk3 = Matmul(kap3,adjTk3) 
 kap3Clam = Matmul(kap3,Conjg(lam)) 
 Tk1adjTv = Matmul(Tk1,adjTv) 
 Tk1Clam = Matmul(Tk1,Conjg(lam)) 
 Tk2adjTv = Matmul(Tk2,adjTv) 
 Tk2Clam = Matmul(Tk2,Conjg(lam)) 
 Tk3adjTv = Matmul(Tk3,adjTv) 
 Tk3Clam = Matmul(Tk3,Conjg(lam)) 
 md2YdadjYu = Matmul(md2,YdadjYu) 
 me2YeCYv = Matmul(me2,YeCYv) 
 ml2YvadjYv = Matmul(ml2,YvadjYv) 
 mu2YuadjYd = Matmul(mu2,YuadjYd) 
 mv2TpYvadjYe = Matmul(mv2,TpYvadjYe) 
 mv2kap1adjkap2 = Matmul(mv2,kap1adjkap2) 
 mv2kap1adjkap3 = Matmul(mv2,kap1adjkap3) 
 mv2kap1Clam = Matmul(mv2,kap1Clam) 
 mv2kap1Ckap1 = Matmul(mv2,kap1Ckap1) 
 mv2kap1Ckap2 = Matmul(mv2,kap1Ckap2) 
 mv2kap1Ckap3 = Matmul(mv2,kap1Ckap3) 
 mv2kap2adjkap1 = Matmul(mv2,kap2adjkap1) 
 mv2kap2adjkap3 = Matmul(mv2,kap2adjkap3) 
 mv2kap2Clam = Matmul(mv2,kap2Clam) 
 mv2kap2Ckap1 = Matmul(mv2,kap2Ckap1) 
 mv2kap2Ckap2 = Matmul(mv2,kap2Ckap2) 
 mv2kap2Ckap3 = Matmul(mv2,kap2Ckap3) 
 mv2kap3adjkap1 = Matmul(mv2,kap3adjkap1) 
 mv2kap3adjkap2 = Matmul(mv2,kap3adjkap2) 
 mv2kap3Clam = Matmul(mv2,kap3Clam) 
 mv2kap3Ckap1 = Matmul(mv2,kap3Ckap1) 
 mv2kap3Ckap2 = Matmul(mv2,kap3Ckap2) 
 mv2kap3Ckap3 = Matmul(mv2,kap3Ckap3) 
 Ydmq2adjYu = Matmul(Yd,mq2adjYu) 
 YdadjYdCmd2 = Matmul(Yd,adjYdCmd2) 
 YdadjYumu2 = Matmul(Yd,adjYumu2) 
 YdadjTdTd = Matmul(Yd,adjTdTd) 
 YdCmq2adjYd = Matmul(Yd,Cmq2adjYd) 
Forall(i2=1:3)  YdCmq2adjYd(i2,i2) =  Real(YdCmq2adjYd(i2,i2),dp) 
 Yeml2CYv = Matmul(Ye,ml2CYv) 
 YeadjYeCme2 = Matmul(Ye,adjYeCme2) 
 YeadjTeTe = Matmul(Ye,adjTeTe) 
 YeCml2adjYe = Matmul(Ye,Cml2adjYe) 
Forall(i2=1:3)  YeCml2adjYe(i2,i2) =  Real(YeCml2adjYe(i2,i2),dp) 
 YeCYvmv2 = Matmul(Ye,CYvmv2) 
 YeCYvlam = Matmul(Ye,CYvlam) 
 YeCYvTlam = Matmul(Ye,CYvTlam) 
 YeCTvTlam = Matmul(Ye,CTvTlam) 
 Yumq2adjYd = Matmul(Yu,mq2adjYd) 
 YuadjYdmd2 = Matmul(Yu,adjYdmd2) 
 YuadjYuCmu2 = Matmul(Yu,adjYuCmu2) 
 YuadjTuTu = Matmul(Yu,adjTuTu) 
 YuCmq2adjYu = Matmul(Yu,Cmq2adjYu) 
Forall(i2=1:3)  YuCmq2adjYu(i2,i2) =  Real(YuCmq2adjYu(i2,i2),dp) 
 YvadjYvmlHd2 = Matmul(Yv,adjYvmlHd2) 
 Yvadjkap1mv2 = Matmul(Yv,adjkap1mv2) 
 Yvadjkap1kap1 = Matmul(Yv,adjkap1kap1) 
 Yvadjkap1kap2 = Matmul(Yv,adjkap1kap2) 
 Yvadjkap1kap3 = Matmul(Yv,adjkap1kap3) 
 Yvadjkap1Tk1 = Matmul(Yv,adjkap1Tk1) 
 Yvadjkap1Tk2 = Matmul(Yv,adjkap1Tk2) 
 Yvadjkap1Tk3 = Matmul(Yv,adjkap1Tk3) 
 Yvadjkap2mv2 = Matmul(Yv,adjkap2mv2) 
 Yvadjkap2kap1 = Matmul(Yv,adjkap2kap1) 
 Yvadjkap2kap2 = Matmul(Yv,adjkap2kap2) 
 Yvadjkap2kap3 = Matmul(Yv,adjkap2kap3) 
 Yvadjkap2Tk1 = Matmul(Yv,adjkap2Tk1) 
 Yvadjkap2Tk2 = Matmul(Yv,adjkap2Tk2) 
 Yvadjkap2Tk3 = Matmul(Yv,adjkap2Tk3) 
 Yvadjkap3mv2 = Matmul(Yv,adjkap3mv2) 
 Yvadjkap3kap1 = Matmul(Yv,adjkap3kap1) 
 Yvadjkap3kap2 = Matmul(Yv,adjkap3kap2) 
 Yvadjkap3kap3 = Matmul(Yv,adjkap3kap3) 
 Yvadjkap3Tk1 = Matmul(Yv,adjkap3Tk1) 
 Yvadjkap3Tk2 = Matmul(Yv,adjkap3Tk2) 
 Yvadjkap3Tk3 = Matmul(Yv,adjkap3Tk3) 
 YvCmv2adjkap1 = Matmul(Yv,Cmv2adjkap1) 
 YvCmv2adjkap2 = Matmul(Yv,Cmv2adjkap2) 
 YvCmv2adjkap3 = Matmul(Yv,Cmv2adjkap3) 
 YvCmv2Clam = Matmul(Yv,Cmv2Clam) 
 YvCkap1lam = Matmul(Yv,Ckap1lam) 
 YvCkap1kap2 = Matmul(Yv,Ckap1kap2) 
 YvCkap1kap3 = Matmul(Yv,Ckap1kap3) 
 YvCkap1Tk1 = Matmul(Yv,Ckap1Tk1) 
 YvCkap1Tk2 = Matmul(Yv,Ckap1Tk2) 
 YvCkap1Tk3 = Matmul(Yv,Ckap1Tk3) 
 YvCkap2lam = Matmul(Yv,Ckap2lam) 
 YvCkap2kap1 = Matmul(Yv,Ckap2kap1) 
 YvCkap2kap3 = Matmul(Yv,Ckap2kap3) 
 YvCkap2Tk1 = Matmul(Yv,Ckap2Tk1) 
 YvCkap2Tk2 = Matmul(Yv,Ckap2Tk2) 
 YvCkap2Tk3 = Matmul(Yv,Ckap2Tk3) 
 YvCkap3lam = Matmul(Yv,Ckap3lam) 
 YvCkap3kap1 = Matmul(Yv,Ckap3kap1) 
 YvCkap3kap2 = Matmul(Yv,Ckap3kap2) 
 YvCkap3Tk1 = Matmul(Yv,Ckap3Tk1) 
 YvCkap3Tk2 = Matmul(Yv,Ckap3Tk2) 
 YvCkap3Tk3 = Matmul(Yv,Ckap3Tk3) 
 adjYdYdadjYd = Matmul(adjYd,YdadjYd) 
 adjYdYdadjYu = Matmul(adjYd,YdadjYu) 
 adjYdYdadjTd = Matmul(adjYd,YdadjTd) 
 adjYdYdadjTu = Matmul(adjYd,YdadjTu) 
 adjYdTdadjYd = Matmul(adjYd,TdadjYd) 
 adjYdTdadjYu = Matmul(adjYd,TdadjYu) 
 adjYdTdadjTd = Matmul(adjYd,TdadjTd) 
 adjYdTdadjTu = Matmul(adjYd,TdadjTu) 
 adjYeYeadjYe = Matmul(adjYe,YeadjYe) 
 adjYeYeadjTe = Matmul(adjYe,YeadjTe) 
 adjYeYeCYv = Matmul(adjYe,YeCYv) 
 adjYeYeCTv = Matmul(adjYe,YeCTv) 
 adjYeTeadjYe = Matmul(adjYe,TeadjYe) 
 adjYeTeadjTe = Matmul(adjYe,TeadjTe) 
 adjYeTeCYv = Matmul(adjYe,TeCYv) 
 adjYeTeCTv = Matmul(adjYe,TeCTv) 
 adjYuYuadjYd = Matmul(adjYu,YuadjYd) 
 adjYuYuadjYu = Matmul(adjYu,YuadjYu) 
 adjYuYuadjTd = Matmul(adjYu,YuadjTd) 
 adjYuYuadjTu = Matmul(adjYu,YuadjTu) 
 adjYuTuadjYd = Matmul(adjYu,TuadjYd) 
 adjYuTuadjYu = Matmul(adjYu,TuadjYu) 
 adjYuTuadjTd = Matmul(adjYu,TuadjTd) 
 adjYuTuadjTu = Matmul(adjYu,TuadjTu) 
 adjYvYvadjYv = Matmul(adjYv,YvadjYv) 
 adjYvYvadjkap1 = Matmul(adjYv,Yvadjkap1) 
 adjYvYvadjkap2 = Matmul(adjYv,Yvadjkap2) 
 adjYvYvadjkap3 = Matmul(adjYv,Yvadjkap3) 
 adjYvYvadjTk1 = Matmul(adjYv,YvadjTk1) 
 adjYvYvadjTk2 = Matmul(adjYv,YvadjTk2) 
 adjYvYvadjTk3 = Matmul(adjYv,YvadjTk3) 
 adjYvYvClam = Matmul(adjYv,YvClam) 
 adjYvYvCTlam = Matmul(adjYv,YvCTlam) 
 adjYvYvCkap1 = Matmul(adjYv,YvCkap1) 
 adjYvYvCkap2 = Matmul(adjYv,YvCkap2) 
 adjYvYvCkap3 = Matmul(adjYv,YvCkap3) 
 adjYvTvadjYv = Matmul(adjYv,TvadjYv) 
 adjYvTvadjTv = Matmul(adjYv,TvadjTv) 
 adjYvTvClam = Matmul(adjYv,TvClam) 
 adjYvTvCTlam = Matmul(adjYv,TvCTlam) 
 adjYvTpYeCYe = Matmul(adjYv,TpYeCYe) 
 adjYvTpTeCTe = Matmul(adjYv,TpTeCTe) 
 adjTdYdadjYd = Matmul(adjTd,YdadjYd) 
 adjTdYdadjYu = Matmul(adjTd,YdadjYu) 
 adjTdTdadjYd = Matmul(adjTd,TdadjYd) 
 adjTdTdadjYu = Matmul(adjTd,TdadjYu) 
 adjTeYeadjYe = Matmul(adjTe,YeadjYe) 
 adjTeYeCYv = Matmul(adjTe,YeCYv) 
 adjTeTeadjYe = Matmul(adjTe,TeadjYe) 
 adjTeTeCYv = Matmul(adjTe,TeCYv) 
 adjTuYuadjYd = Matmul(adjTu,YuadjYd) 
 adjTuYuadjYu = Matmul(adjTu,YuadjYu) 
 adjTuTuadjYd = Matmul(adjTu,TuadjYd) 
 adjTuTuadjYu = Matmul(adjTu,TuadjYu) 
 adjTvYvClam = Matmul(adjTv,YvClam) 
 adjTvTvadjYv = Matmul(adjTv,TvadjYv) 
 adjTvTvClam = Matmul(adjTv,TvClam) 
 adjkap1mv2lam = Matmul(adjkap1,mv2lam) 
 adjkap1mv2kap1 = Matmul(adjkap1,mv2kap1) 
 adjkap1mv2kap2 = Matmul(adjkap1,mv2kap2) 
 adjkap1mv2kap3 = Matmul(adjkap1,mv2kap3) 
 adjkap1TpYvml2 = Matmul(adjkap1,TpYvml2) 
 adjkap1kap1Clam = Matmul(adjkap1,kap1Clam) 
 adjkap1kap2Clam = Matmul(adjkap1,kap2Clam) 
 adjkap1kap3Clam = Matmul(adjkap1,kap3Clam) 
 adjkap2mv2lam = Matmul(adjkap2,mv2lam) 
 adjkap2mv2kap1 = Matmul(adjkap2,mv2kap1) 
 adjkap2mv2kap2 = Matmul(adjkap2,mv2kap2) 
 adjkap2mv2kap3 = Matmul(adjkap2,mv2kap3) 
 adjkap2TpYvml2 = Matmul(adjkap2,TpYvml2) 
 adjkap2kap1Clam = Matmul(adjkap2,kap1Clam) 
 adjkap2kap2Clam = Matmul(adjkap2,kap2Clam) 
 adjkap2kap3Clam = Matmul(adjkap2,kap3Clam) 
 adjkap3mv2lam = Matmul(adjkap3,mv2lam) 
 adjkap3mv2kap1 = Matmul(adjkap3,mv2kap1) 
 adjkap3mv2kap2 = Matmul(adjkap3,mv2kap2) 
 adjkap3mv2kap3 = Matmul(adjkap3,mv2kap3) 
 adjkap3TpYvml2 = Matmul(adjkap3,TpYvml2) 
 adjkap3kap1Clam = Matmul(adjkap3,kap1Clam) 
 adjkap3kap2Clam = Matmul(adjkap3,kap2Clam) 
 adjkap3kap3Clam = Matmul(adjkap3,kap3Clam) 
 Cml2YvadjYv = Matmul(Conjg(ml2),YvadjYv) 
 Cml2Yvadjkap1 = Matmul(Conjg(ml2),Yvadjkap1) 
 Cml2Yvadjkap2 = Matmul(Conjg(ml2),Yvadjkap2) 
 Cml2Yvadjkap3 = Matmul(Conjg(ml2),Yvadjkap3) 
 Cml2YvClam = Matmul(Conjg(ml2),YvClam) 
 Cml2TpYeCYe = Matmul(Conjg(ml2),TpYeCYe) 
 Cmv2Ckap1kap1 = Matmul(Conjg(mv2),Ckap1kap1) 
 Cmv2Ckap1kap2 = Matmul(Conjg(mv2),Ckap1kap2) 
 Cmv2Ckap1kap3 = Matmul(Conjg(mv2),Ckap1kap3) 
 Cmv2Ckap2kap1 = Matmul(Conjg(mv2),Ckap2kap1) 
 Cmv2Ckap2kap2 = Matmul(Conjg(mv2),Ckap2kap2) 
 Cmv2Ckap2kap3 = Matmul(Conjg(mv2),Ckap2kap3) 
 Cmv2Ckap3kap1 = Matmul(Conjg(mv2),Ckap3kap1) 
 Cmv2Ckap3kap2 = Matmul(Conjg(mv2),Ckap3kap2) 
 Cmv2Ckap3kap3 = Matmul(Conjg(mv2),Ckap3kap3) 
 CYeYvClam = Matmul(Conjg(Ye),YvClam) 
 CYeTvClam = Matmul(Conjg(Ye),TvClam) 
 CYvTpYvadjYe = Matmul(Conjg(Yv),TpYvadjYe) 
 CYvTpYvadjTe = Matmul(Conjg(Yv),TpYvadjTe) 
 CYvTpYvCYv = Matmul(Conjg(Yv),TpYvCYv) 
 CYvTpYvCTv = Matmul(Conjg(Yv),TpYvCTv) 
 CYvTpTvadjTe = Matmul(Conjg(Yv),TpTvadjTe) 
 CYvTpTvCTv = Matmul(Conjg(Yv),TpTvCTv) 
 CYvTpkap1Ckap1 = Matmul(Conjg(Yv),Tpkap1Ckap1) 
 CYvTpkap1Ckap2 = Matmul(Conjg(Yv),Tpkap1Ckap2) 
 CYvTpkap1Ckap3 = Matmul(Conjg(Yv),Tpkap1Ckap3) 
 CYvTpkap2Ckap1 = Matmul(Conjg(Yv),Tpkap2Ckap1) 
 CYvTpkap2Ckap2 = Matmul(Conjg(Yv),Tpkap2Ckap2) 
 CYvTpkap2Ckap3 = Matmul(Conjg(Yv),Tpkap2Ckap3) 
 CYvTpkap3Ckap1 = Matmul(Conjg(Yv),Tpkap3Ckap1) 
 CYvTpkap3Ckap2 = Matmul(Conjg(Yv),Tpkap3Ckap2) 
 CYvTpkap3Ckap3 = Matmul(Conjg(Yv),Tpkap3Ckap3) 
 CTvTpYvadjYe = Matmul(Conjg(Tv),TpYvadjYe) 
 CTvTpYvCYv = Matmul(Conjg(Tv),TpYvCYv) 
 CTvTpTvadjYe = Matmul(Conjg(Tv),TpTvadjYe) 
 CTvTpTvCYv = Matmul(Conjg(Tv),TpTvCYv) 
 Ckap1mv2lam = Matmul(Conjg(kap1),mv2lam) 
 Ckap1TpYvCYv = Matmul(Conjg(kap1),TpYvCYv) 
 Ckap1TpYvCTv = Matmul(Conjg(kap1),TpYvCTv) 
 Ckap1Tpkap1adjYv = Matmul(Conjg(kap1),Tpkap1adjYv) 
 Ckap1Tpkap1Clam = Matmul(Conjg(kap1),Tpkap1Clam) 
 Ckap1TpTk1adjYv = Matmul(Conjg(kap1),TpTk1adjYv) 
 Ckap1kap1adjYv = Matmul(Conjg(kap1),kap1adjYv) 
 Ckap1kap1Ckap1 = Matmul(Conjg(kap1),kap1Ckap1) 
 Ckap1kap1Ckap2 = Matmul(Conjg(kap1),kap1Ckap2) 
 Ckap1kap1Ckap3 = Matmul(Conjg(kap1),kap1Ckap3) 
 Ckap1kap2Ckap1 = Matmul(Conjg(kap1),kap2Ckap1) 
 Ckap1kap2Ckap2 = Matmul(Conjg(kap1),kap2Ckap2) 
 Ckap1kap2Ckap3 = Matmul(Conjg(kap1),kap2Ckap3) 
 Ckap1kap3Ckap1 = Matmul(Conjg(kap1),kap3Ckap1) 
 Ckap1kap3Ckap2 = Matmul(Conjg(kap1),kap3Ckap2) 
 Ckap1kap3Ckap3 = Matmul(Conjg(kap1),kap3Ckap3) 
 Ckap1Tk1adjTv = Matmul(Conjg(kap1),Tk1adjTv) 
 Ckap2mv2lam = Matmul(Conjg(kap2),mv2lam) 
 Ckap2TpYvCYv = Matmul(Conjg(kap2),TpYvCYv) 
 Ckap2TpYvCTv = Matmul(Conjg(kap2),TpYvCTv) 
 Ckap2Tpkap2adjYv = Matmul(Conjg(kap2),Tpkap2adjYv) 
 Ckap2Tpkap2Clam = Matmul(Conjg(kap2),Tpkap2Clam) 
 Ckap2TpTk2adjYv = Matmul(Conjg(kap2),TpTk2adjYv) 
 Ckap2kap1Ckap1 = Matmul(Conjg(kap2),kap1Ckap1) 
 Ckap2kap1Ckap2 = Matmul(Conjg(kap2),kap1Ckap2) 
 Ckap2kap1Ckap3 = Matmul(Conjg(kap2),kap1Ckap3) 
 Ckap2kap2adjYv = Matmul(Conjg(kap2),kap2adjYv) 
 Ckap2kap2Ckap1 = Matmul(Conjg(kap2),kap2Ckap1) 
 Ckap2kap2Ckap2 = Matmul(Conjg(kap2),kap2Ckap2) 
 Ckap2kap2Ckap3 = Matmul(Conjg(kap2),kap2Ckap3) 
 Ckap2kap3Ckap1 = Matmul(Conjg(kap2),kap3Ckap1) 
 Ckap2kap3Ckap2 = Matmul(Conjg(kap2),kap3Ckap2) 
 Ckap2kap3Ckap3 = Matmul(Conjg(kap2),kap3Ckap3) 
 Ckap2Tk2adjTv = Matmul(Conjg(kap2),Tk2adjTv) 
 Ckap3mv2lam = Matmul(Conjg(kap3),mv2lam) 
 Ckap3TpYvCYv = Matmul(Conjg(kap3),TpYvCYv) 
 Ckap3TpYvCTv = Matmul(Conjg(kap3),TpYvCTv) 
 Ckap3Tpkap3adjYv = Matmul(Conjg(kap3),Tpkap3adjYv) 
 Ckap3Tpkap3Clam = Matmul(Conjg(kap3),Tpkap3Clam) 
 Ckap3TpTk3adjYv = Matmul(Conjg(kap3),TpTk3adjYv) 
 Ckap3kap1Ckap1 = Matmul(Conjg(kap3),kap1Ckap1) 
 Ckap3kap1Ckap2 = Matmul(Conjg(kap3),kap1Ckap2) 
 Ckap3kap1Ckap3 = Matmul(Conjg(kap3),kap1Ckap3) 
 Ckap3kap2Ckap1 = Matmul(Conjg(kap3),kap2Ckap1) 
 Ckap3kap2Ckap2 = Matmul(Conjg(kap3),kap2Ckap2) 
 Ckap3kap2Ckap3 = Matmul(Conjg(kap3),kap2Ckap3) 
 Ckap3kap3adjYv = Matmul(Conjg(kap3),kap3adjYv) 
 Ckap3kap3Ckap1 = Matmul(Conjg(kap3),kap3Ckap1) 
 Ckap3kap3Ckap2 = Matmul(Conjg(kap3),kap3Ckap2) 
 Ckap3kap3Ckap3 = Matmul(Conjg(kap3),kap3Ckap3) 
 Ckap3Tk3adjTv = Matmul(Conjg(kap3),Tk3adjTv) 
 CTk1TpTk1adjYv = Matmul(Conjg(Tk1),TpTk1adjYv) 
 CTk2TpTk2adjYv = Matmul(Conjg(Tk2),TpTk2adjYv) 
 CTk3TpTk3adjYv = Matmul(Conjg(Tk3),TpTk3adjYv) 
 TdadjTdYd = Matmul(Td,adjTdYd) 
 TeadjTeYe = Matmul(Te,adjTeYe) 
 TeCYvlam = Matmul(Te,CYvlam) 
 TeCTvlam = Matmul(Te,CTvlam) 
 TeCTvTpYv = Matmul(Te,CTvTpYv) 
 TuadjTuYu = Matmul(Tu,adjTuYu) 
 Tvadjkap1kap1 = Matmul(Tv,adjkap1kap1) 
 Tvadjkap1kap2 = Matmul(Tv,adjkap1kap2) 
 Tvadjkap1kap3 = Matmul(Tv,adjkap1kap3) 
 Tvadjkap2kap1 = Matmul(Tv,adjkap2kap1) 
 Tvadjkap2kap2 = Matmul(Tv,adjkap2kap2) 
 Tvadjkap2kap3 = Matmul(Tv,adjkap2kap3) 
 Tvadjkap3kap1 = Matmul(Tv,adjkap3kap1) 
 Tvadjkap3kap2 = Matmul(Tv,adjkap3kap2) 
 Tvadjkap3kap3 = Matmul(Tv,adjkap3kap3) 
 TvCkap1kap1 = Matmul(Tv,Ckap1kap1) 
 TvCkap1kap2 = Matmul(Tv,Ckap1kap2) 
 TvCkap1kap3 = Matmul(Tv,Ckap1kap3) 
 TvCkap2kap1 = Matmul(Tv,Ckap2kap1) 
 TvCkap2kap2 = Matmul(Tv,Ckap2kap2) 
 TvCkap2kap3 = Matmul(Tv,Ckap2kap3) 
 TvCkap3kap1 = Matmul(Tv,Ckap3kap1) 
 TvCkap3kap2 = Matmul(Tv,Ckap3kap2) 
 TvCkap3kap3 = Matmul(Tv,Ckap3kap3) 
 TpYeCme2CYe = Matmul(Transpose(Ye),Cme2CYe) 
Forall(i2=1:3)  TpYeCme2CYe(i2,i2) =  Real(TpYeCme2CYe(i2,i2),dp) 
 TpYeCYeCml2 = Matmul(Transpose(Ye),CYeCml2) 
 TpYeCTeTv = Matmul(Transpose(Ye),CTeTv) 
 TpYvml2adjYe = Matmul(Transpose(Yv),ml2adjYe) 
 TpYvadjYeme2 = Matmul(Transpose(Yv),adjYeme2) 
 TpYvadjYeYe = Matmul(Transpose(Yv),adjYeYe) 
 TpYvadjYeTe = Matmul(Transpose(Yv),adjYeTe) 
 TpYvCYvTpYv = Matmul(Transpose(Yv),CYvTpYv) 
 TpYvCYvTpTv = Matmul(Transpose(Yv),CYvTpTv) 
 TpYvCYvTpkap1 = Matmul(Transpose(Yv),CYvTpkap1) 
 TpYvCYvTpkap2 = Matmul(Transpose(Yv),CYvTpkap2) 
 TpYvCYvTpkap3 = Matmul(Transpose(Yv),CYvTpkap3) 
 TpYvCYvTpTk1 = Matmul(Transpose(Yv),CYvTpTk1) 
 TpYvCYvTpTk2 = Matmul(Transpose(Yv),CYvTpTk2) 
 TpYvCYvTpTk3 = Matmul(Transpose(Yv),CYvTpTk3) 
 TpYvCTvTlam = Matmul(Transpose(Yv),CTvTlam) 
 TpYvCTvTpTv = Matmul(Transpose(Yv),CTvTpTv) 
 TpTvadjYeYe = Matmul(Transpose(Tv),adjYeYe) 
 TpTvCYvTpYv = Matmul(Transpose(Tv),CYvTpYv) 
 TpTvCTvlam = Matmul(Transpose(Tv),CTvlam) 
 TpTvCTvTpYv = Matmul(Transpose(Tv),CTvTpYv) 
 Tpkap1adjYvYv = Matmul(Transpose(kap1),adjYvYv) 
 Tpkap1adjkap1lam = Matmul(Transpose(kap1),adjkap1lam) 
 Tpkap1adjkap1Tlam = Matmul(Transpose(kap1),adjkap1Tlam) 
 Tpkap1adjkap2lam = Matmul(Transpose(kap1),adjkap2lam) 
 Tpkap1adjkap2Tlam = Matmul(Transpose(kap1),adjkap2Tlam) 
 Tpkap1adjkap3lam = Matmul(Transpose(kap1),adjkap3lam) 
 Tpkap1adjkap3Tlam = Matmul(Transpose(kap1),adjkap3Tlam) 
 Tpkap1Cmv2adjYv = Matmul(Transpose(kap1),Cmv2adjYv) 
 Tpkap1Cmv2adjkap1 = Matmul(Transpose(kap1),Cmv2adjkap1) 
 Tpkap1Cmv2adjkap2 = Matmul(Transpose(kap1),Cmv2adjkap2) 
 Tpkap1Cmv2adjkap3 = Matmul(Transpose(kap1),Cmv2adjkap3) 
 Tpkap1Cmv2Clam = Matmul(Transpose(kap1),Cmv2Clam) 
 Tpkap1Ckap1lam = Matmul(Transpose(kap1),Ckap1lam) 
 Tpkap1Ckap1Tlam = Matmul(Transpose(kap1),Ckap1Tlam) 
 Tpkap1Ckap1TpTv = Matmul(Transpose(kap1),Ckap1TpTv) 
 Tpkap1Ckap1TpTk1 = Matmul(Transpose(kap1),Ckap1TpTk1) 
 Tpkap1Ckap1TpTk2 = Matmul(Transpose(kap1),Ckap1TpTk2) 
 Tpkap1Ckap1TpTk3 = Matmul(Transpose(kap1),Ckap1TpTk3) 
 Tpkap1Ckap2lam = Matmul(Transpose(kap1),Ckap2lam) 
 Tpkap1Ckap2Tlam = Matmul(Transpose(kap1),Ckap2Tlam) 
 Tpkap1Ckap2TpTv = Matmul(Transpose(kap1),Ckap2TpTv) 
 Tpkap1Ckap2TpTk1 = Matmul(Transpose(kap1),Ckap2TpTk1) 
 Tpkap1Ckap2TpTk2 = Matmul(Transpose(kap1),Ckap2TpTk2) 
 Tpkap1Ckap2TpTk3 = Matmul(Transpose(kap1),Ckap2TpTk3) 
 Tpkap1Ckap3lam = Matmul(Transpose(kap1),Ckap3lam) 
 Tpkap1Ckap3Tlam = Matmul(Transpose(kap1),Ckap3Tlam) 
 Tpkap1Ckap3TpTv = Matmul(Transpose(kap1),Ckap3TpTv) 
 Tpkap1Ckap3TpTk1 = Matmul(Transpose(kap1),Ckap3TpTk1) 
 Tpkap1Ckap3TpTk2 = Matmul(Transpose(kap1),Ckap3TpTk2) 
 Tpkap1Ckap3TpTk3 = Matmul(Transpose(kap1),Ckap3TpTk3) 
 Tpkap2adjYvYv = Matmul(Transpose(kap2),adjYvYv) 
 Tpkap2adjkap1lam = Matmul(Transpose(kap2),adjkap1lam) 
 Tpkap2adjkap1Tlam = Matmul(Transpose(kap2),adjkap1Tlam) 
 Tpkap2adjkap2lam = Matmul(Transpose(kap2),adjkap2lam) 
 Tpkap2adjkap2Tlam = Matmul(Transpose(kap2),adjkap2Tlam) 
 Tpkap2adjkap3lam = Matmul(Transpose(kap2),adjkap3lam) 
 Tpkap2adjkap3Tlam = Matmul(Transpose(kap2),adjkap3Tlam) 
 Tpkap2Cmv2adjYv = Matmul(Transpose(kap2),Cmv2adjYv) 
 Tpkap2Cmv2adjkap1 = Matmul(Transpose(kap2),Cmv2adjkap1) 
 Tpkap2Cmv2adjkap2 = Matmul(Transpose(kap2),Cmv2adjkap2) 
 Tpkap2Cmv2adjkap3 = Matmul(Transpose(kap2),Cmv2adjkap3) 
 Tpkap2Cmv2Clam = Matmul(Transpose(kap2),Cmv2Clam) 
 Tpkap2Ckap1lam = Matmul(Transpose(kap2),Ckap1lam) 
 Tpkap2Ckap1Tlam = Matmul(Transpose(kap2),Ckap1Tlam) 
 Tpkap2Ckap1TpTv = Matmul(Transpose(kap2),Ckap1TpTv) 
 Tpkap2Ckap1TpTk1 = Matmul(Transpose(kap2),Ckap1TpTk1) 
 Tpkap2Ckap1TpTk2 = Matmul(Transpose(kap2),Ckap1TpTk2) 
 Tpkap2Ckap1TpTk3 = Matmul(Transpose(kap2),Ckap1TpTk3) 
 Tpkap2Ckap2lam = Matmul(Transpose(kap2),Ckap2lam) 
 Tpkap2Ckap2Tlam = Matmul(Transpose(kap2),Ckap2Tlam) 
 Tpkap2Ckap2TpTv = Matmul(Transpose(kap2),Ckap2TpTv) 
 Tpkap2Ckap2TpTk1 = Matmul(Transpose(kap2),Ckap2TpTk1) 
 Tpkap2Ckap2TpTk2 = Matmul(Transpose(kap2),Ckap2TpTk2) 
 Tpkap2Ckap2TpTk3 = Matmul(Transpose(kap2),Ckap2TpTk3) 
 Tpkap2Ckap3lam = Matmul(Transpose(kap2),Ckap3lam) 
 Tpkap2Ckap3Tlam = Matmul(Transpose(kap2),Ckap3Tlam) 
 Tpkap2Ckap3TpTv = Matmul(Transpose(kap2),Ckap3TpTv) 
 Tpkap2Ckap3TpTk1 = Matmul(Transpose(kap2),Ckap3TpTk1) 
 Tpkap2Ckap3TpTk2 = Matmul(Transpose(kap2),Ckap3TpTk2) 
 Tpkap2Ckap3TpTk3 = Matmul(Transpose(kap2),Ckap3TpTk3) 
 Tpkap3adjYvYv = Matmul(Transpose(kap3),adjYvYv) 
 Tpkap3adjkap1lam = Matmul(Transpose(kap3),adjkap1lam) 
 Tpkap3adjkap1Tlam = Matmul(Transpose(kap3),adjkap1Tlam) 
 Tpkap3adjkap2lam = Matmul(Transpose(kap3),adjkap2lam) 
 Tpkap3adjkap2Tlam = Matmul(Transpose(kap3),adjkap2Tlam) 
 Tpkap3adjkap3lam = Matmul(Transpose(kap3),adjkap3lam) 
 Tpkap3adjkap3Tlam = Matmul(Transpose(kap3),adjkap3Tlam) 
 Tpkap3Cmv2adjYv = Matmul(Transpose(kap3),Cmv2adjYv) 
 Tpkap3Cmv2adjkap1 = Matmul(Transpose(kap3),Cmv2adjkap1) 
 Tpkap3Cmv2adjkap2 = Matmul(Transpose(kap3),Cmv2adjkap2) 
 Tpkap3Cmv2adjkap3 = Matmul(Transpose(kap3),Cmv2adjkap3) 
 Tpkap3Cmv2Clam = Matmul(Transpose(kap3),Cmv2Clam) 
 Tpkap3Ckap1lam = Matmul(Transpose(kap3),Ckap1lam) 
 Tpkap3Ckap1Tlam = Matmul(Transpose(kap3),Ckap1Tlam) 
 Tpkap3Ckap1TpTv = Matmul(Transpose(kap3),Ckap1TpTv) 
 Tpkap3Ckap1TpTk1 = Matmul(Transpose(kap3),Ckap1TpTk1) 
 Tpkap3Ckap1TpTk2 = Matmul(Transpose(kap3),Ckap1TpTk2) 
 Tpkap3Ckap1TpTk3 = Matmul(Transpose(kap3),Ckap1TpTk3) 
 Tpkap3Ckap2lam = Matmul(Transpose(kap3),Ckap2lam) 
 Tpkap3Ckap2Tlam = Matmul(Transpose(kap3),Ckap2Tlam) 
 Tpkap3Ckap2TpTv = Matmul(Transpose(kap3),Ckap2TpTv) 
 Tpkap3Ckap2TpTk1 = Matmul(Transpose(kap3),Ckap2TpTk1) 
 Tpkap3Ckap2TpTk2 = Matmul(Transpose(kap3),Ckap2TpTk2) 
 Tpkap3Ckap2TpTk3 = Matmul(Transpose(kap3),Ckap2TpTk3) 
 Tpkap3Ckap3lam = Matmul(Transpose(kap3),Ckap3lam) 
 Tpkap3Ckap3Tlam = Matmul(Transpose(kap3),Ckap3Tlam) 
 Tpkap3Ckap3TpTv = Matmul(Transpose(kap3),Ckap3TpTv) 
 Tpkap3Ckap3TpTk1 = Matmul(Transpose(kap3),Ckap3TpTk1) 
 Tpkap3Ckap3TpTk2 = Matmul(Transpose(kap3),Ckap3TpTk2) 
 Tpkap3Ckap3TpTk3 = Matmul(Transpose(kap3),Ckap3TpTk3) 
 TpTk1adjYvYv = Matmul(Transpose(Tk1),adjYvYv) 
 TpTk1adjkap1lam = Matmul(Transpose(Tk1),adjkap1lam) 
 TpTk1adjkap2lam = Matmul(Transpose(Tk1),adjkap2lam) 
 TpTk1adjkap3lam = Matmul(Transpose(Tk1),adjkap3lam) 
 TpTk1Ckap1lam = Matmul(Transpose(Tk1),Ckap1lam) 
 TpTk1Ckap2lam = Matmul(Transpose(Tk1),Ckap2lam) 
 TpTk1Ckap3lam = Matmul(Transpose(Tk1),Ckap3lam) 
 TpTk2adjYvYv = Matmul(Transpose(Tk2),adjYvYv) 
 TpTk2adjkap1lam = Matmul(Transpose(Tk2),adjkap1lam) 
 TpTk2adjkap2lam = Matmul(Transpose(Tk2),adjkap2lam) 
 TpTk2adjkap3lam = Matmul(Transpose(Tk2),adjkap3lam) 
 TpTk2Ckap1lam = Matmul(Transpose(Tk2),Ckap1lam) 
 TpTk2Ckap2lam = Matmul(Transpose(Tk2),Ckap2lam) 
 TpTk2Ckap3lam = Matmul(Transpose(Tk2),Ckap3lam) 
 TpTk3adjYvYv = Matmul(Transpose(Tk3),adjYvYv) 
 TpTk3adjkap1lam = Matmul(Transpose(Tk3),adjkap1lam) 
 TpTk3adjkap2lam = Matmul(Transpose(Tk3),adjkap2lam) 
 TpTk3adjkap3lam = Matmul(Transpose(Tk3),adjkap3lam) 
 TpTk3Ckap1lam = Matmul(Transpose(Tk3),Ckap1lam) 
 TpTk3Ckap2lam = Matmul(Transpose(Tk3),Ckap2lam) 
 TpTk3Ckap3lam = Matmul(Transpose(Tk3),Ckap3lam) 
 kap1adjYvCml2 = Matmul(kap1,adjYvCml2) 
 kap1adjYvTv = Matmul(kap1,adjYvTv) 
 kap1adjkap1lam = Matmul(kap1,adjkap1lam) 
 kap1adjkap1TpYv = Matmul(kap1,adjkap1TpYv) 
 kap1adjkap1Tpkap1 = Matmul(kap1,adjkap1Tpkap1) 
 kap1adjkap1Tpkap2 = Matmul(kap1,adjkap1Tpkap2) 
 kap1adjkap1Tpkap3 = Matmul(kap1,adjkap1Tpkap3) 
 kap1Cmv2Clam = Matmul(kap1,Cmv2Clam) 
 kap1Ckap1lam = Matmul(kap1,Ckap1lam) 
 kap1Ckap1kap1 = Matmul(kap1,Ckap1kap1) 
 kap1Ckap1kap2 = Matmul(kap1,Ckap1kap2) 
 kap1Ckap1kap3 = Matmul(kap1,Ckap1kap3) 
 kap1Ckap1Tk1 = Matmul(kap1,Ckap1Tk1) 
 kap1Ckap1Tk2 = Matmul(kap1,Ckap1Tk2) 
 kap1Ckap1Tk3 = Matmul(kap1,Ckap1Tk3) 
 kap1Ckap2lam = Matmul(kap1,Ckap2lam) 
 kap1Ckap2kap1 = Matmul(kap1,Ckap2kap1) 
 kap1Ckap2kap2 = Matmul(kap1,Ckap2kap2) 
 kap1Ckap2kap3 = Matmul(kap1,Ckap2kap3) 
 kap1Ckap2Tk1 = Matmul(kap1,Ckap2Tk1) 
 kap1Ckap2Tk2 = Matmul(kap1,Ckap2Tk2) 
 kap1Ckap2Tk3 = Matmul(kap1,Ckap2Tk3) 
 kap1Ckap3lam = Matmul(kap1,Ckap3lam) 
 kap1Ckap3kap1 = Matmul(kap1,Ckap3kap1) 
 kap1Ckap3kap2 = Matmul(kap1,Ckap3kap2) 
 kap1Ckap3kap3 = Matmul(kap1,Ckap3kap3) 
 kap1Ckap3Tk1 = Matmul(kap1,Ckap3Tk1) 
 kap1Ckap3Tk2 = Matmul(kap1,Ckap3Tk2) 
 kap1Ckap3Tk3 = Matmul(kap1,Ckap3Tk3) 
 kap2adjYvCml2 = Matmul(kap2,adjYvCml2) 
 kap2adjYvTv = Matmul(kap2,adjYvTv) 
 kap2adjkap2lam = Matmul(kap2,adjkap2lam) 
 kap2adjkap2TpYv = Matmul(kap2,adjkap2TpYv) 
 kap2adjkap2Tpkap1 = Matmul(kap2,adjkap2Tpkap1) 
 kap2adjkap2Tpkap2 = Matmul(kap2,adjkap2Tpkap2) 
 kap2adjkap2Tpkap3 = Matmul(kap2,adjkap2Tpkap3) 
 kap2Cmv2Clam = Matmul(kap2,Cmv2Clam) 
 kap2Ckap1lam = Matmul(kap2,Ckap1lam) 
 kap2Ckap1kap1 = Matmul(kap2,Ckap1kap1) 
 kap2Ckap1kap2 = Matmul(kap2,Ckap1kap2) 
 kap2Ckap1kap3 = Matmul(kap2,Ckap1kap3) 
 kap2Ckap1Tk1 = Matmul(kap2,Ckap1Tk1) 
 kap2Ckap1Tk2 = Matmul(kap2,Ckap1Tk2) 
 kap2Ckap1Tk3 = Matmul(kap2,Ckap1Tk3) 
 kap2Ckap2lam = Matmul(kap2,Ckap2lam) 
 kap2Ckap2kap1 = Matmul(kap2,Ckap2kap1) 
 kap2Ckap2kap2 = Matmul(kap2,Ckap2kap2) 
 kap2Ckap2kap3 = Matmul(kap2,Ckap2kap3) 
 kap2Ckap2Tk1 = Matmul(kap2,Ckap2Tk1) 
 kap2Ckap2Tk2 = Matmul(kap2,Ckap2Tk2) 
 kap2Ckap2Tk3 = Matmul(kap2,Ckap2Tk3) 
 kap2Ckap3lam = Matmul(kap2,Ckap3lam) 
 kap2Ckap3kap1 = Matmul(kap2,Ckap3kap1) 
 kap2Ckap3kap2 = Matmul(kap2,Ckap3kap2) 
 kap2Ckap3kap3 = Matmul(kap2,Ckap3kap3) 
 kap2Ckap3Tk1 = Matmul(kap2,Ckap3Tk1) 
 kap2Ckap3Tk2 = Matmul(kap2,Ckap3Tk2) 
 kap2Ckap3Tk3 = Matmul(kap2,Ckap3Tk3) 
 kap3adjYvCml2 = Matmul(kap3,adjYvCml2) 
 kap3adjYvTv = Matmul(kap3,adjYvTv) 
 kap3adjkap3lam = Matmul(kap3,adjkap3lam) 
 kap3adjkap3TpYv = Matmul(kap3,adjkap3TpYv) 
 kap3adjkap3Tpkap1 = Matmul(kap3,adjkap3Tpkap1) 
 kap3adjkap3Tpkap2 = Matmul(kap3,adjkap3Tpkap2) 
 kap3adjkap3Tpkap3 = Matmul(kap3,adjkap3Tpkap3) 
 kap3Cmv2Clam = Matmul(kap3,Cmv2Clam) 
 kap3Ckap1lam = Matmul(kap3,Ckap1lam) 
 kap3Ckap1kap1 = Matmul(kap3,Ckap1kap1) 
 kap3Ckap1kap2 = Matmul(kap3,Ckap1kap2) 
 kap3Ckap1kap3 = Matmul(kap3,Ckap1kap3) 
 kap3Ckap1Tk1 = Matmul(kap3,Ckap1Tk1) 
 kap3Ckap1Tk2 = Matmul(kap3,Ckap1Tk2) 
 kap3Ckap1Tk3 = Matmul(kap3,Ckap1Tk3) 
 kap3Ckap2lam = Matmul(kap3,Ckap2lam) 
 kap3Ckap2kap1 = Matmul(kap3,Ckap2kap1) 
 kap3Ckap2kap2 = Matmul(kap3,Ckap2kap2) 
 kap3Ckap2kap3 = Matmul(kap3,Ckap2kap3) 
 kap3Ckap2Tk1 = Matmul(kap3,Ckap2Tk1) 
 kap3Ckap2Tk2 = Matmul(kap3,Ckap2Tk2) 
 kap3Ckap2Tk3 = Matmul(kap3,Ckap2Tk3) 
 kap3Ckap3lam = Matmul(kap3,Ckap3lam) 
 kap3Ckap3kap1 = Matmul(kap3,Ckap3kap1) 
 kap3Ckap3kap2 = Matmul(kap3,Ckap3kap2) 
 kap3Ckap3kap3 = Matmul(kap3,Ckap3kap3) 
 kap3Ckap3Tk1 = Matmul(kap3,Ckap3Tk1) 
 kap3Ckap3Tk2 = Matmul(kap3,Ckap3Tk2) 
 kap3Ckap3Tk3 = Matmul(kap3,Ckap3Tk3) 
 Tk1adjkap1TpYv = Matmul(Tk1,adjkap1TpYv) 
 Tk1adjTk1lam = Matmul(Tk1,adjTk1lam) 
 Tk1adjTk1TpYv = Matmul(Tk1,adjTk1TpYv) 
 Tk1Ckap1kap1 = Matmul(Tk1,Ckap1kap1) 
 Tk1Ckap1kap2 = Matmul(Tk1,Ckap1kap2) 
 Tk1Ckap1kap3 = Matmul(Tk1,Ckap1kap3) 
 Tk1Ckap2kap1 = Matmul(Tk1,Ckap2kap1) 
 Tk1Ckap2kap2 = Matmul(Tk1,Ckap2kap2) 
 Tk1Ckap2kap3 = Matmul(Tk1,Ckap2kap3) 
 Tk1Ckap3kap1 = Matmul(Tk1,Ckap3kap1) 
 Tk1Ckap3kap2 = Matmul(Tk1,Ckap3kap2) 
 Tk1Ckap3kap3 = Matmul(Tk1,Ckap3kap3) 
 Tk2adjkap2TpYv = Matmul(Tk2,adjkap2TpYv) 
 Tk2adjTk2lam = Matmul(Tk2,adjTk2lam) 
 Tk2adjTk2TpYv = Matmul(Tk2,adjTk2TpYv) 
 Tk2Ckap1kap1 = Matmul(Tk2,Ckap1kap1) 
 Tk2Ckap1kap2 = Matmul(Tk2,Ckap1kap2) 
 Tk2Ckap1kap3 = Matmul(Tk2,Ckap1kap3) 
 Tk2Ckap2kap1 = Matmul(Tk2,Ckap2kap1) 
 Tk2Ckap2kap2 = Matmul(Tk2,Ckap2kap2) 
 Tk2Ckap2kap3 = Matmul(Tk2,Ckap2kap3) 
 Tk2Ckap3kap1 = Matmul(Tk2,Ckap3kap1) 
 Tk2Ckap3kap2 = Matmul(Tk2,Ckap3kap2) 
 Tk2Ckap3kap3 = Matmul(Tk2,Ckap3kap3) 
 Tk3adjkap3TpYv = Matmul(Tk3,adjkap3TpYv) 
 Tk3adjTk3lam = Matmul(Tk3,adjTk3lam) 
 Tk3adjTk3TpYv = Matmul(Tk3,adjTk3TpYv) 
 Tk3Ckap1kap1 = Matmul(Tk3,Ckap1kap1) 
 Tk3Ckap1kap2 = Matmul(Tk3,Ckap1kap2) 
 Tk3Ckap1kap3 = Matmul(Tk3,Ckap1kap3) 
 Tk3Ckap2kap1 = Matmul(Tk3,Ckap2kap1) 
 Tk3Ckap2kap2 = Matmul(Tk3,Ckap2kap2) 
 Tk3Ckap2kap3 = Matmul(Tk3,Ckap2kap3) 
 Tk3Ckap3kap1 = Matmul(Tk3,Ckap3kap1) 
 Tk3Ckap3kap2 = Matmul(Tk3,Ckap3kap2) 
 Tk3Ckap3kap3 = Matmul(Tk3,Ckap3kap3) 
 md2YdadjYdYd = Matmul(md2,YdadjYdYd) 
 me2YeadjYeYe = Matmul(me2,YeadjYeYe) 
 me2YeCYvlam = Matmul(me2,YeCYvlam) 
 ml2adjYeYeadjYe = Matmul(ml2,adjYeYeadjYe) 
 ml2adjYeYeCYv = Matmul(ml2,adjYeYeCYv) 
 ml2CYvTpYvadjYe = Matmul(ml2,CYvTpYvadjYe) 
 ml2CYvTpYvCYv = Matmul(ml2,CYvTpYvCYv) 
 mq2adjYdYdadjYd = Matmul(mq2,adjYdYdadjYd) 
 mq2adjYdYdadjYu = Matmul(mq2,adjYdYdadjYu) 
 mq2adjYuYuadjYd = Matmul(mq2,adjYuYuadjYd) 
 mq2adjYuYuadjYu = Matmul(mq2,adjYuYuadjYu) 
 mu2YuadjYuYu = Matmul(mu2,YuadjYuYu) 
 mv2TpYvCYvlam = Matmul(mv2,TpYvCYvlam) 
 mv2TpYvCYvTpYv = Matmul(mv2,TpYvCYvTpYv) 
 mv2kap1adjkap1TpYv = Matmul(mv2,kap1adjkap1TpYv) 
 mv2kap2adjkap2TpYv = Matmul(mv2,kap2adjkap2TpYv) 
 mv2kap3adjkap3TpYv = Matmul(mv2,kap3adjkap3TpYv) 
 Ydmq2adjYdYd = Matmul(Yd,mq2adjYdYd) 
 YdadjYdmd2Yd = Matmul(Yd,adjYdmd2Yd) 
 YdadjYdYdmq2 = Matmul(Yd,adjYdYdmq2) 
 YdadjYdYdadjYd = Matmul(Yd,adjYdYdadjYd) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTdadjYd = Matmul(Yd,adjYdTdadjYd) 
 YdadjYdTdadjTd = Matmul(Yd,adjYdTdadjTd) 
 YdadjYuYuadjYd = Matmul(Yd,adjYuYuadjYd) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YdadjYuTuadjYd = Matmul(Yd,adjYuTuadjYd) 
 YdadjYuTuadjTd = Matmul(Yd,adjYuTuadjTd) 
 YdadjTdTdadjYd = Matmul(Yd,adjTdTdadjYd) 
 YdadjTuTuadjYd = Matmul(Yd,adjTuTuadjYd) 
 Yeml2adjYeYe = Matmul(Ye,ml2adjYeYe) 
 Yeml2CYvlam = Matmul(Ye,ml2CYvlam) 
 YeadjYeme2Ye = Matmul(Ye,adjYeme2Ye) 
 YeadjYeYeml2 = Matmul(Ye,adjYeYeml2) 
 YeadjYeYeadjYe = Matmul(Ye,adjYeYeadjYe) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYeTeadjYe = Matmul(Ye,adjYeTeadjYe) 
 YeadjYeTeadjTe = Matmul(Ye,adjYeTeadjTe) 
 YeadjTeTeadjYe = Matmul(Ye,adjTeTeadjYe) 
 YeCYvmv2lam = Matmul(Ye,CYvmv2lam) 
 YeCYvTpYvadjYe = Matmul(Ye,CYvTpYvadjYe) 
Forall(i2=1:3)  YeCYvTpYvadjYe(i2,i2) =  Real(YeCYvTpYvadjYe(i2,i2),dp) 
 YeCYvTpTvadjTe = Matmul(Ye,CYvTpTvadjTe) 
 YeCTvTpTvadjYe = Matmul(Ye,CTvTpTvadjYe) 
Forall(i2=1:3)  YeCTvTpTvadjYe(i2,i2) =  Real(YeCTvTpTvadjYe(i2,i2),dp) 
 Yumq2adjYuYu = Matmul(Yu,mq2adjYuYu) 
 YuadjYdYdadjYu = Matmul(Yu,adjYdYdadjYu) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYdTdadjYu = Matmul(Yu,adjYdTdadjYu) 
 YuadjYdTdadjTu = Matmul(Yu,adjYdTdadjTu) 
 YuadjYumu2Yu = Matmul(Yu,adjYumu2Yu) 
 YuadjYuYumq2 = Matmul(Yu,adjYuYumq2) 
 YuadjYuYuadjYu = Matmul(Yu,adjYuYuadjYu) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YuadjYuTuadjYu = Matmul(Yu,adjYuTuadjYu) 
 YuadjYuTuadjTu = Matmul(Yu,adjYuTuadjTu) 
 YuadjTdTdadjYu = Matmul(Yu,adjTdTdadjYu) 
 YuadjTuTuadjYu = Matmul(Yu,adjTuTuadjYu) 
 YvadjYvYvadjYv = Matmul(Yv,adjYvYvadjYv) 
Forall(i2=1:3)  YvadjYvYvadjYv(i2,i2) =  Real(YvadjYvYvadjYv(i2,i2),dp) 
 YvadjYvYvClam = Matmul(Yv,adjYvYvClam) 
 YvadjYvTvadjYv = Matmul(Yv,adjYvTvadjYv) 
 YvadjYvTvadjTv = Matmul(Yv,adjYvTvadjTv) 
 YvadjYvTvClam = Matmul(Yv,adjYvTvClam) 
 YvadjYvTpYeCYe = Matmul(Yv,adjYvTpYeCYe) 
 YvadjYvTpTeCTe = Matmul(Yv,adjYvTpTeCTe) 
 YvadjTvTvadjYv = Matmul(Yv,adjTvTvadjYv) 
 Yvadjkap1mv2kap1 = Matmul(Yv,adjkap1mv2kap1) 
 Yvadjkap1mv2kap2 = Matmul(Yv,adjkap1mv2kap2) 
 Yvadjkap1mv2kap3 = Matmul(Yv,adjkap1mv2kap3) 
 Yvadjkap2mv2kap1 = Matmul(Yv,adjkap2mv2kap1) 
 Yvadjkap2mv2kap2 = Matmul(Yv,adjkap2mv2kap2) 
 Yvadjkap2mv2kap3 = Matmul(Yv,adjkap2mv2kap3) 
 Yvadjkap3mv2kap1 = Matmul(Yv,adjkap3mv2kap1) 
 Yvadjkap3mv2kap2 = Matmul(Yv,adjkap3mv2kap2) 
 Yvadjkap3mv2kap3 = Matmul(Yv,adjkap3mv2kap3) 
 YvCmv2Ckap1kap1 = Matmul(Yv,Cmv2Ckap1kap1) 
 YvCmv2Ckap1kap2 = Matmul(Yv,Cmv2Ckap1kap2) 
 YvCmv2Ckap1kap3 = Matmul(Yv,Cmv2Ckap1kap3) 
 YvCmv2Ckap2kap1 = Matmul(Yv,Cmv2Ckap2kap1) 
 YvCmv2Ckap2kap2 = Matmul(Yv,Cmv2Ckap2kap2) 
 YvCmv2Ckap2kap3 = Matmul(Yv,Cmv2Ckap2kap3) 
 YvCmv2Ckap3kap1 = Matmul(Yv,Cmv2Ckap3kap1) 
 YvCmv2Ckap3kap2 = Matmul(Yv,Cmv2Ckap3kap2) 
 YvCmv2Ckap3kap3 = Matmul(Yv,Cmv2Ckap3kap3) 
 YvCkap1Tpkap1adjYv = Matmul(Yv,Ckap1Tpkap1adjYv) 
Forall(i2=1:3)  YvCkap1Tpkap1adjYv(i2,i2) =  Real(YvCkap1Tpkap1adjYv(i2,i2),dp) 
 YvCkap1Tpkap1Clam = Matmul(Yv,Ckap1Tpkap1Clam) 
 YvCkap1TpTk1adjYv = Matmul(Yv,Ckap1TpTk1adjYv) 
 YvCkap1kap1adjYv = Matmul(Yv,Ckap1kap1adjYv) 
 YvCkap1Tk1adjTv = Matmul(Yv,Ckap1Tk1adjTv) 
 YvCkap2Tpkap2adjYv = Matmul(Yv,Ckap2Tpkap2adjYv) 
Forall(i2=1:3)  YvCkap2Tpkap2adjYv(i2,i2) =  Real(YvCkap2Tpkap2adjYv(i2,i2),dp) 
 YvCkap2Tpkap2Clam = Matmul(Yv,Ckap2Tpkap2Clam) 
 YvCkap2TpTk2adjYv = Matmul(Yv,Ckap2TpTk2adjYv) 
 YvCkap2kap2adjYv = Matmul(Yv,Ckap2kap2adjYv) 
 YvCkap2Tk2adjTv = Matmul(Yv,Ckap2Tk2adjTv) 
 YvCkap3Tpkap3adjYv = Matmul(Yv,Ckap3Tpkap3adjYv) 
Forall(i2=1:3)  YvCkap3Tpkap3adjYv(i2,i2) =  Real(YvCkap3Tpkap3adjYv(i2,i2),dp) 
 YvCkap3Tpkap3Clam = Matmul(Yv,Ckap3Tpkap3Clam) 
 YvCkap3TpTk3adjYv = Matmul(Yv,Ckap3TpTk3adjYv) 
 YvCkap3kap3adjYv = Matmul(Yv,Ckap3kap3adjYv) 
 YvCkap3Tk3adjTv = Matmul(Yv,Ckap3Tk3adjTv) 
 YvCTk1TpTk1adjYv = Matmul(Yv,CTk1TpTk1adjYv) 
Forall(i2=1:3)  YvCTk1TpTk1adjYv(i2,i2) =  Real(YvCTk1TpTk1adjYv(i2,i2),dp) 
 YvCTk2TpTk2adjYv = Matmul(Yv,CTk2TpTk2adjYv) 
Forall(i2=1:3)  YvCTk2TpTk2adjYv(i2,i2) =  Real(YvCTk2TpTk2adjYv(i2,i2),dp) 
 YvCTk3TpTk3adjYv = Matmul(Yv,CTk3TpTk3adjYv) 
Forall(i2=1:3)  YvCTk3TpTk3adjYv(i2,i2) =  Real(YvCTk3TpTk3adjYv(i2,i2),dp) 
 adjYdmd2YdadjYd = Matmul(adjYd,md2YdadjYd) 
 adjYdmd2YdadjYu = Matmul(adjYd,md2YdadjYu) 
 adjYdYdmq2adjYd = Matmul(adjYd,Ydmq2adjYd) 
 adjYdYdmq2adjYu = Matmul(adjYd,Ydmq2adjYu) 
 adjYdYdadjYdmd2 = Matmul(adjYd,YdadjYdmd2) 
 adjYdYdadjYdYd = Matmul(adjYd,YdadjYdYd) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdTd = Matmul(adjYd,YdadjYdTd) 
 adjYdYdadjYumu2 = Matmul(adjYd,YdadjYumu2) 
 adjYdYdadjYuYu = Matmul(adjYd,YdadjYuYu) 
 adjYdYdadjYuTu = Matmul(adjYd,YdadjYuTu) 
 adjYdYdadjTdTd = Matmul(adjYd,YdadjTdTd) 
 adjYdTdadjYdYd = Matmul(adjYd,TdadjYdYd) 
 adjYdTdadjYuYu = Matmul(adjYd,TdadjYuYu) 
 adjYdTdadjTdYd = Matmul(adjYd,TdadjTdYd) 
 adjYeme2YeadjYe = Matmul(adjYe,me2YeadjYe) 
 adjYeme2YeCYv = Matmul(adjYe,me2YeCYv) 
 adjYeYeml2adjYe = Matmul(adjYe,Yeml2adjYe) 
 adjYeYeml2CYv = Matmul(adjYe,Yeml2CYv) 
 adjYeYeadjYeme2 = Matmul(adjYe,YeadjYeme2) 
 adjYeYeadjYeYe = Matmul(adjYe,YeadjYeYe) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYeTe = Matmul(adjYe,YeadjYeTe) 
 adjYeYeadjTeTe = Matmul(adjYe,YeadjTeTe) 
 adjYeYeCYvmv2 = Matmul(adjYe,YeCYvmv2) 
 adjYeYeCYvlam = Matmul(adjYe,YeCYvlam) 
 adjYeYeCTvTlam = Matmul(adjYe,YeCTvTlam) 
 adjYeTeadjYeYe = Matmul(adjYe,TeadjYeYe) 
 adjYeTeadjTeYe = Matmul(adjYe,TeadjTeYe) 
 adjYeTeCYvlam = Matmul(adjYe,TeCYvlam) 
 adjYeTeCYvTpYv = Matmul(adjYe,TeCYvTpYv) 
 adjYeTeCTvlam = Matmul(adjYe,TeCTvlam) 
 adjYeTeCTvTpYv = Matmul(adjYe,TeCTvTpYv) 
 adjYumu2YuadjYd = Matmul(adjYu,mu2YuadjYd) 
 adjYumu2YuadjYu = Matmul(adjYu,mu2YuadjYu) 
 adjYuYumq2adjYd = Matmul(adjYu,Yumq2adjYd) 
 adjYuYumq2adjYu = Matmul(adjYu,Yumq2adjYu) 
 adjYuYuadjYdmd2 = Matmul(adjYu,YuadjYdmd2) 
 adjYuYuadjYdYd = Matmul(adjYu,YuadjYdYd) 
 adjYuYuadjYdTd = Matmul(adjYu,YuadjYdTd) 
 adjYuYuadjYumu2 = Matmul(adjYu,YuadjYumu2) 
 adjYuYuadjYuYu = Matmul(adjYu,YuadjYuYu) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYuYuadjYuTu = Matmul(adjYu,YuadjYuTu) 
 adjYuYuadjTuTu = Matmul(adjYu,YuadjTuTu) 
 adjYuTuadjYdYd = Matmul(adjYu,TuadjYdYd) 
 adjYuTuadjYuYu = Matmul(adjYu,TuadjYuYu) 
 adjYuTuadjTuYu = Matmul(adjYu,TuadjTuYu) 
 adjYvYvadjYvmlHd2 = Matmul(adjYv,YvadjYvmlHd2) 
 adjYvYvadjYvYv = Matmul(adjYv,YvadjYvYv) 
Forall(i2=1:3)  adjYvYvadjYvYv(i2,i2) =  Real(adjYvYvadjYvYv(i2,i2),dp) 
 adjYvYvadjYvTv = Matmul(adjYv,YvadjYvTv) 
 adjYvYvadjkap1mv2 = Matmul(adjYv,Yvadjkap1mv2) 
 adjYvYvadjkap2mv2 = Matmul(adjYv,Yvadjkap2mv2) 
 adjYvYvadjkap3mv2 = Matmul(adjYv,Yvadjkap3mv2) 
 adjYvYvCmv2adjkap1 = Matmul(adjYv,YvCmv2adjkap1) 
 adjYvYvCmv2adjkap2 = Matmul(adjYv,YvCmv2adjkap2) 
 adjYvYvCmv2adjkap3 = Matmul(adjYv,YvCmv2adjkap3) 
 adjYvYvCmv2Clam = Matmul(adjYv,YvCmv2Clam) 
 adjYvYvCkap1lam = Matmul(adjYv,YvCkap1lam) 
 adjYvYvCkap1kap2 = Matmul(adjYv,YvCkap1kap2) 
 adjYvYvCkap1kap3 = Matmul(adjYv,YvCkap1kap3) 
 adjYvYvCkap2lam = Matmul(adjYv,YvCkap2lam) 
 adjYvYvCkap2kap1 = Matmul(adjYv,YvCkap2kap1) 
 adjYvYvCkap2kap3 = Matmul(adjYv,YvCkap2kap3) 
 adjYvYvCkap3lam = Matmul(adjYv,YvCkap3lam) 
 adjYvYvCkap3kap1 = Matmul(adjYv,YvCkap3kap1) 
 adjYvYvCkap3kap2 = Matmul(adjYv,YvCkap3kap2) 
 adjYvCml2YvadjYv = Matmul(adjYv,Cml2YvadjYv) 
 adjYvCml2Yvadjkap1 = Matmul(adjYv,Cml2Yvadjkap1) 
 adjYvCml2Yvadjkap2 = Matmul(adjYv,Cml2Yvadjkap2) 
 adjYvCml2Yvadjkap3 = Matmul(adjYv,Cml2Yvadjkap3) 
 adjYvCml2YvClam = Matmul(adjYv,Cml2YvClam) 
 adjYvCml2TpYeCYe = Matmul(adjYv,Cml2TpYeCYe) 
 adjYvTvadjYvYv = Matmul(adjYv,TvadjYvYv) 
 adjYvTvadjkap1kap1 = Matmul(adjYv,Tvadjkap1kap1) 
 adjYvTvadjkap1kap2 = Matmul(adjYv,Tvadjkap1kap2) 
 adjYvTvadjkap1kap3 = Matmul(adjYv,Tvadjkap1kap3) 
 adjYvTvadjkap2kap1 = Matmul(adjYv,Tvadjkap2kap1) 
 adjYvTvadjkap2kap2 = Matmul(adjYv,Tvadjkap2kap2) 
 adjYvTvadjkap2kap3 = Matmul(adjYv,Tvadjkap2kap3) 
 adjYvTvadjkap3kap1 = Matmul(adjYv,Tvadjkap3kap1) 
 adjYvTvadjkap3kap2 = Matmul(adjYv,Tvadjkap3kap2) 
 adjYvTvadjkap3kap3 = Matmul(adjYv,Tvadjkap3kap3) 
 adjYvTpYeCme2CYe = Matmul(adjYv,TpYeCme2CYe) 
 adjYvTpYeCYeYv = Matmul(adjYv,TpYeCYeYv) 
Forall(i2=1:3)  adjYvTpYeCYeYv(i2,i2) =  Real(adjYvTpYeCYeYv(i2,i2),dp) 
 adjYvTpYeCYeCml2 = Matmul(adjYv,TpYeCYeCml2) 
 adjYvTpYeCYeTv = Matmul(adjYv,TpYeCYeTv) 
 adjYvTpYeCTeTv = Matmul(adjYv,TpYeCTeTv) 
 adjYvTpTeCYeYv = Matmul(adjYv,TpTeCYeYv) 
 adjTdYdadjYdTd = Matmul(adjTd,YdadjYdTd) 
 adjTdTdadjYdYd = Matmul(adjTd,TdadjYdYd) 
 adjTeYeadjYeTe = Matmul(adjTe,YeadjYeTe) 
 adjTeYeCYvTlam = Matmul(adjTe,YeCYvTlam) 
 adjTeTeadjYeYe = Matmul(adjTe,TeadjYeYe) 
 adjTeTeCYvlam = Matmul(adjTe,TeCYvlam) 
 adjTuYuadjYuTu = Matmul(adjTu,YuadjYuTu) 
 adjTuTuadjYuYu = Matmul(adjTu,TuadjYuYu) 
 adjTvYvCkap1Tk2 = Matmul(adjTv,YvCkap1Tk2) 
 adjTvYvCkap1Tk3 = Matmul(adjTv,YvCkap1Tk3) 
 adjTvYvCkap2Tk1 = Matmul(adjTv,YvCkap2Tk1) 
 adjTvYvCkap2Tk3 = Matmul(adjTv,YvCkap2Tk3) 
 adjTvYvCkap3Tk1 = Matmul(adjTv,YvCkap3Tk1) 
 adjTvYvCkap3Tk2 = Matmul(adjTv,YvCkap3Tk2) 
 adjkap1Tpkap1Cmv2Clam = Matmul(adjkap1,Tpkap1Cmv2Clam) 
 adjkap1Tpkap2Cmv2Clam = Matmul(adjkap1,Tpkap2Cmv2Clam) 
 adjkap1Tpkap3Cmv2Clam = Matmul(adjkap1,Tpkap3Cmv2Clam) 
 adjkap2Tpkap1Cmv2Clam = Matmul(adjkap2,Tpkap1Cmv2Clam) 
 adjkap2Tpkap2Cmv2Clam = Matmul(adjkap2,Tpkap2Cmv2Clam) 
 adjkap2Tpkap3Cmv2Clam = Matmul(adjkap2,Tpkap3Cmv2Clam) 
 adjkap3Tpkap1Cmv2Clam = Matmul(adjkap3,Tpkap1Cmv2Clam) 
 adjkap3Tpkap2Cmv2Clam = Matmul(adjkap3,Tpkap2Cmv2Clam) 
 adjkap3Tpkap3Cmv2Clam = Matmul(adjkap3,Tpkap3Cmv2Clam) 
 Cml2YvCkap1kap2 = Matmul(Conjg(ml2),YvCkap1kap2) 
 Cml2YvCkap1kap3 = Matmul(Conjg(ml2),YvCkap1kap3) 
 Cml2YvCkap2kap1 = Matmul(Conjg(ml2),YvCkap2kap1) 
 Cml2YvCkap2kap3 = Matmul(Conjg(ml2),YvCkap2kap3) 
 Cml2YvCkap3kap1 = Matmul(Conjg(ml2),YvCkap3kap1) 
 Cml2YvCkap3kap2 = Matmul(Conjg(ml2),YvCkap3kap2) 
 Cmv2adjYvYvadjYv = Matmul(Conjg(mv2),adjYvYvadjYv) 
 Cmv2adjYvYvadjkap1 = Matmul(Conjg(mv2),adjYvYvadjkap1) 
 Cmv2adjYvYvadjkap2 = Matmul(Conjg(mv2),adjYvYvadjkap2) 
 Cmv2adjYvYvadjkap3 = Matmul(Conjg(mv2),adjYvYvadjkap3) 
 Cmv2adjYvYvClam = Matmul(Conjg(mv2),adjYvYvClam) 
 Cmv2adjYvTpYeCYe = Matmul(Conjg(mv2),adjYvTpYeCYe) 
 Cmv2Ckap1kap1adjYv = Matmul(Conjg(mv2),Ckap1kap1adjYv) 
 Cmv2Ckap2kap2adjYv = Matmul(Conjg(mv2),Ckap2kap2adjYv) 
 Cmv2Ckap3kap3adjYv = Matmul(Conjg(mv2),Ckap3kap3adjYv) 
 CYeTpYeCYeYv = Matmul(Conjg(Ye),TpYeCYeYv) 
 CYeTpYeCYeTv = Matmul(Conjg(Ye),TpYeCYeTv) 
 CYeTpTeCYeYv = Matmul(Conjg(Ye),TpTeCYeYv) 
 CYvmv2TpYvadjYe = Matmul(Conjg(Yv),mv2TpYvadjYe) 
 CYvmv2TpYvCYv = Matmul(Conjg(Yv),mv2TpYvCYv) 
 CYvmv2kap1Ckap1 = Matmul(Conjg(Yv),mv2kap1Ckap1) 
 CYvmv2kap1Ckap2 = Matmul(Conjg(Yv),mv2kap1Ckap2) 
 CYvmv2kap1Ckap3 = Matmul(Conjg(Yv),mv2kap1Ckap3) 
 CYvmv2kap2Ckap1 = Matmul(Conjg(Yv),mv2kap2Ckap1) 
 CYvmv2kap2Ckap2 = Matmul(Conjg(Yv),mv2kap2Ckap2) 
 CYvmv2kap2Ckap3 = Matmul(Conjg(Yv),mv2kap2Ckap3) 
 CYvmv2kap3Ckap1 = Matmul(Conjg(Yv),mv2kap3Ckap1) 
 CYvmv2kap3Ckap2 = Matmul(Conjg(Yv),mv2kap3Ckap2) 
 CYvmv2kap3Ckap3 = Matmul(Conjg(Yv),mv2kap3Ckap3) 
 CYvTpYvml2adjYe = Matmul(Conjg(Yv),TpYvml2adjYe) 
 CYvTpYvml2CYv = Matmul(Conjg(Yv),TpYvml2CYv) 
 CYvTpYvadjYeme2 = Matmul(Conjg(Yv),TpYvadjYeme2) 
 CYvTpYvadjYeYe = Matmul(Conjg(Yv),TpYvadjYeYe) 
 CYvTpYvadjYeTe = Matmul(Conjg(Yv),TpYvadjYeTe) 
 CYvTpYvCYvmv2 = Matmul(Conjg(Yv),TpYvCYvmv2) 
 CYvTpYvCYvlam = Matmul(Conjg(Yv),TpYvCYvlam) 
 CYvTpYvCYvTlam = Matmul(Conjg(Yv),TpYvCYvTlam) 
 CYvTpYvCYvTpYv = Matmul(Conjg(Yv),TpYvCYvTpYv) 
Forall(i2=1:3)  CYvTpYvCYvTpYv(i2,i2) =  Real(CYvTpYvCYvTpYv(i2,i2),dp) 
 CYvTpYvCYvTpTv = Matmul(Conjg(Yv),TpYvCYvTpTv) 
 CYvTpYvCTvTlam = Matmul(Conjg(Yv),TpYvCTvTlam) 
 CYvTpYvCTvTpTv = Matmul(Conjg(Yv),TpYvCTvTpTv) 
 CYvTpTvadjYeYe = Matmul(Conjg(Yv),TpTvadjYeYe) 
 CYvTpTvCYvlam = Matmul(Conjg(Yv),TpTvCYvlam) 
 CYvTpTvCYvTpYv = Matmul(Conjg(Yv),TpTvCYvTpYv) 
 CYvTpTvCTvlam = Matmul(Conjg(Yv),TpTvCTvlam) 
 CYvTpTvCTvTpYv = Matmul(Conjg(Yv),TpTvCTvTpYv) 
Forall(i2=1:3)  CYvTpTvCTvTpYv(i2,i2) =  Real(CYvTpTvCTvTpYv(i2,i2),dp) 
 CYvkap1adjkap1lam = Matmul(Conjg(Yv),kap1adjkap1lam) 
 CYvkap1adjkap1TpYv = Matmul(Conjg(Yv),kap1adjkap1TpYv) 
 CYvkap2adjkap2lam = Matmul(Conjg(Yv),kap2adjkap2lam) 
 CYvkap2adjkap2TpYv = Matmul(Conjg(Yv),kap2adjkap2TpYv) 
 CYvkap3adjkap3lam = Matmul(Conjg(Yv),kap3adjkap3lam) 
 CYvkap3adjkap3TpYv = Matmul(Conjg(Yv),kap3adjkap3TpYv) 
 CYvTk1adjkap1lam = Matmul(Conjg(Yv),Tk1adjkap1lam) 
 CYvTk1adjkap1TpYv = Matmul(Conjg(Yv),Tk1adjkap1TpYv) 
 CYvTk1adjTk1lam = Matmul(Conjg(Yv),Tk1adjTk1lam) 
 CYvTk1adjTk1TpYv = Matmul(Conjg(Yv),Tk1adjTk1TpYv) 
 CYvTk2adjkap2lam = Matmul(Conjg(Yv),Tk2adjkap2lam) 
 CYvTk2adjkap2TpYv = Matmul(Conjg(Yv),Tk2adjkap2TpYv) 
 CYvTk2adjTk2lam = Matmul(Conjg(Yv),Tk2adjTk2lam) 
 CYvTk2adjTk2TpYv = Matmul(Conjg(Yv),Tk2adjTk2TpYv) 
 CYvTk3adjkap3lam = Matmul(Conjg(Yv),Tk3adjkap3lam) 
 CYvTk3adjkap3TpYv = Matmul(Conjg(Yv),Tk3adjkap3TpYv) 
 CYvTk3adjTk3lam = Matmul(Conjg(Yv),Tk3adjTk3lam) 
 CYvTk3adjTk3TpYv = Matmul(Conjg(Yv),Tk3adjTk3TpYv) 
 CTvTpYvCYvTlam = Matmul(Conjg(Tv),TpYvCYvTlam) 
 CTvTpYvCYvTpTv = Matmul(Conjg(Tv),TpYvCYvTpTv) 
Forall(i2=1:3)  CTvTpYvCYvTpTv(i2,i2) =  Real(CTvTpYvCYvTpTv(i2,i2),dp) 
 CTvTpTvCYvlam = Matmul(Conjg(Tv),TpTvCYvlam) 
 CTvTpTvCYvTpYv = Matmul(Conjg(Tv),TpTvCYvTpYv) 
 Ckap1TpYvCYvTpTk1 = Matmul(Conjg(kap1),TpYvCYvTpTk1) 
 Ckap1Tpkap1adjYvYv = Matmul(Conjg(kap1),Tpkap1adjYvYv) 
 Ckap1Tpkap1Cmv2adjYv = Matmul(Conjg(kap1),Tpkap1Cmv2adjYv) 
 Ckap1Tpkap1Cmv2adjkap1 = Matmul(Conjg(kap1),Tpkap1Cmv2adjkap1) 
 Ckap1Tpkap1Cmv2adjkap2 = Matmul(Conjg(kap1),Tpkap1Cmv2adjkap2) 
 Ckap1Tpkap1Cmv2adjkap3 = Matmul(Conjg(kap1),Tpkap1Cmv2adjkap3) 
 Ckap1TpTk1adjYvYv = Matmul(Conjg(kap1),TpTk1adjYvYv) 
 Ckap1kap1adjYvCml2 = Matmul(Conjg(kap1),kap1adjYvCml2) 
 Ckap2TpYvCYvTpTk2 = Matmul(Conjg(kap2),TpYvCYvTpTk2) 
 Ckap2Tpkap2adjYvYv = Matmul(Conjg(kap2),Tpkap2adjYvYv) 
 Ckap2Tpkap2Cmv2adjYv = Matmul(Conjg(kap2),Tpkap2Cmv2adjYv) 
 Ckap2Tpkap2Cmv2adjkap1 = Matmul(Conjg(kap2),Tpkap2Cmv2adjkap1) 
 Ckap2Tpkap2Cmv2adjkap2 = Matmul(Conjg(kap2),Tpkap2Cmv2adjkap2) 
 Ckap2Tpkap2Cmv2adjkap3 = Matmul(Conjg(kap2),Tpkap2Cmv2adjkap3) 
 Ckap2TpTk2adjYvYv = Matmul(Conjg(kap2),TpTk2adjYvYv) 
 Ckap2kap2adjYvCml2 = Matmul(Conjg(kap2),kap2adjYvCml2) 
 Ckap3TpYvCYvTpTk3 = Matmul(Conjg(kap3),TpYvCYvTpTk3) 
 Ckap3Tpkap3adjYvYv = Matmul(Conjg(kap3),Tpkap3adjYvYv) 
 Ckap3Tpkap3Cmv2adjYv = Matmul(Conjg(kap3),Tpkap3Cmv2adjYv) 
 Ckap3Tpkap3Cmv2adjkap1 = Matmul(Conjg(kap3),Tpkap3Cmv2adjkap1) 
 Ckap3Tpkap3Cmv2adjkap2 = Matmul(Conjg(kap3),Tpkap3Cmv2adjkap2) 
 Ckap3Tpkap3Cmv2adjkap3 = Matmul(Conjg(kap3),Tpkap3Cmv2adjkap3) 
 Ckap3TpTk3adjYvYv = Matmul(Conjg(kap3),TpTk3adjYvYv) 
 Ckap3kap3adjYvCml2 = Matmul(Conjg(kap3),kap3adjYvCml2) 
 TdadjYdYdadjTd = Matmul(Td,adjYdYdadjTd) 
 TdadjYuYuadjTd = Matmul(Td,adjYuYuadjTd) 
 TdadjTdYdadjYd = Matmul(Td,adjTdYdadjYd) 
 TdadjTuYuadjYd = Matmul(Td,adjTuYuadjYd) 
 TeadjYeYeadjTe = Matmul(Te,adjYeYeadjTe) 
 TeadjTeYeadjYe = Matmul(Te,adjTeYeadjYe) 
 TeCYvTpYvadjTe = Matmul(Te,CYvTpYvadjTe) 
 TeCTvTpYvadjYe = Matmul(Te,CTvTpYvadjYe) 
 TuadjYdYdadjTu = Matmul(Tu,adjYdYdadjTu) 
 TuadjYuYuadjTu = Matmul(Tu,adjYuYuadjTu) 
 TuadjTdYdadjYu = Matmul(Tu,adjTdYdadjYu) 
 TuadjTuYuadjYu = Matmul(Tu,adjTuYuadjYu) 
 TvadjYvYvClam = Matmul(Tv,adjYvYvClam) 
 TpYeCYeYvClam = Matmul(Transpose(Ye),CYeYvClam) 
 TpYeCYeTvClam = Matmul(Transpose(Ye),CYeTvClam) 
 TpYvml2CYvlam = Matmul(Transpose(Yv),ml2CYvlam) 
 TpYvml2CYvTpYv = Matmul(Transpose(Yv),ml2CYvTpYv) 
 TpYvadjYeYeCYv = Matmul(Transpose(Yv),adjYeYeCYv) 
Forall(i2=1:3)  TpYvadjYeYeCYv(i2,i2) =  Real(TpYvadjYeYeCYv(i2,i2),dp) 
 TpYvadjYeTeCYv = Matmul(Transpose(Yv),adjYeTeCYv) 
 TpYvadjYeTeCTv = Matmul(Transpose(Yv),adjYeTeCTv) 
 TpYvadjTeTeCYv = Matmul(Transpose(Yv),adjTeTeCYv) 
 TpYvCYvmv2lam = Matmul(Transpose(Yv),CYvmv2lam) 
 TpYvCYvmv2TpYv = Matmul(Transpose(Yv),CYvmv2TpYv) 
 TpYvCYvTpYvml2 = Matmul(Transpose(Yv),CYvTpYvml2) 
 TpYvCYvTpYvCYv = Matmul(Transpose(Yv),CYvTpYvCYv) 
Forall(i2=1:3)  TpYvCYvTpYvCYv(i2,i2) =  Real(TpYvCYvTpYvCYv(i2,i2),dp) 
 TpYvCYvTpTvCTv = Matmul(Transpose(Yv),CYvTpTvCTv) 
 TpYvCTvTpTvCYv = Matmul(Transpose(Yv),CTvTpTvCYv) 
Forall(i2=1:3)  TpYvCTvTpTvCYv(i2,i2) =  Real(TpYvCTvTpTvCYv(i2,i2),dp) 
 TpTeCYeYvClam = Matmul(Transpose(Te),CYeYvClam) 
 TpTvadjYeYeCTv = Matmul(Transpose(Tv),adjYeYeCTv) 
Forall(i2=1:3)  TpTvadjYeYeCTv(i2,i2) =  Real(TpTvadjYeYeCTv(i2,i2),dp) 
 TpTvadjTeYeCYv = Matmul(Transpose(Tv),adjTeYeCYv) 
 TpTvCYvTpYvCYv = Matmul(Transpose(Tv),CYvTpYvCYv) 
 TpTvCYvTpYvCTv = Matmul(Transpose(Tv),CYvTpYvCTv) 
Forall(i2=1:3)  TpTvCYvTpYvCTv(i2,i2) =  Real(TpTvCYvTpYvCTv(i2,i2),dp) 
 TpTvCTvTpYvCYv = Matmul(Transpose(Tv),CTvTpYvCYv) 
 Tpkap1adjkap1mv2lam = Matmul(Transpose(kap1),adjkap1mv2lam) 
 Tpkap1adjkap2mv2lam = Matmul(Transpose(kap1),adjkap2mv2lam) 
 Tpkap1adjkap3mv2lam = Matmul(Transpose(kap1),adjkap3mv2lam) 
 Tpkap2adjkap1mv2lam = Matmul(Transpose(kap2),adjkap1mv2lam) 
 Tpkap2adjkap2mv2lam = Matmul(Transpose(kap2),adjkap2mv2lam) 
 Tpkap2adjkap3mv2lam = Matmul(Transpose(kap2),adjkap3mv2lam) 
 Tpkap3adjkap1mv2lam = Matmul(Transpose(kap3),adjkap1mv2lam) 
 Tpkap3adjkap2mv2lam = Matmul(Transpose(kap3),adjkap2mv2lam) 
 Tpkap3adjkap3mv2lam = Matmul(Transpose(kap3),adjkap3mv2lam) 
 kap1adjYvYvadjkap1 = Matmul(kap1,adjYvYvadjkap1) 
 kap1adjYvYvCkap1 = Matmul(kap1,adjYvYvCkap1) 
 kap1adjYvYvCkap2 = Matmul(kap1,adjYvYvCkap2) 
 kap1adjYvYvCkap3 = Matmul(kap1,adjYvYvCkap3) 
 kap1adjkap1TpYvml2 = Matmul(kap1,adjkap1TpYvml2) 
 kap1Cmv2Ckap1kap1 = Matmul(kap1,Cmv2Ckap1kap1) 
 kap1Cmv2Ckap1kap2 = Matmul(kap1,Cmv2Ckap1kap2) 
 kap1Cmv2Ckap1kap3 = Matmul(kap1,Cmv2Ckap1kap3) 
 kap1Cmv2Ckap2kap1 = Matmul(kap1,Cmv2Ckap2kap1) 
 kap1Cmv2Ckap2kap2 = Matmul(kap1,Cmv2Ckap2kap2) 
 kap1Cmv2Ckap2kap3 = Matmul(kap1,Cmv2Ckap2kap3) 
 kap1Cmv2Ckap3kap1 = Matmul(kap1,Cmv2Ckap3kap1) 
 kap1Cmv2Ckap3kap2 = Matmul(kap1,Cmv2Ckap3kap2) 
 kap1Cmv2Ckap3kap3 = Matmul(kap1,Cmv2Ckap3kap3) 
 kap2adjYvYvadjkap2 = Matmul(kap2,adjYvYvadjkap2) 
 kap2adjYvYvCkap1 = Matmul(kap2,adjYvYvCkap1) 
 kap2adjYvYvCkap2 = Matmul(kap2,adjYvYvCkap2) 
 kap2adjYvYvCkap3 = Matmul(kap2,adjYvYvCkap3) 
 kap2adjkap2TpYvml2 = Matmul(kap2,adjkap2TpYvml2) 
 kap2Cmv2Ckap1kap1 = Matmul(kap2,Cmv2Ckap1kap1) 
 kap2Cmv2Ckap1kap2 = Matmul(kap2,Cmv2Ckap1kap2) 
 kap2Cmv2Ckap1kap3 = Matmul(kap2,Cmv2Ckap1kap3) 
 kap2Cmv2Ckap2kap1 = Matmul(kap2,Cmv2Ckap2kap1) 
 kap2Cmv2Ckap2kap2 = Matmul(kap2,Cmv2Ckap2kap2) 
 kap2Cmv2Ckap2kap3 = Matmul(kap2,Cmv2Ckap2kap3) 
 kap2Cmv2Ckap3kap1 = Matmul(kap2,Cmv2Ckap3kap1) 
 kap2Cmv2Ckap3kap2 = Matmul(kap2,Cmv2Ckap3kap2) 
 kap2Cmv2Ckap3kap3 = Matmul(kap2,Cmv2Ckap3kap3) 
 kap3adjYvYvadjkap3 = Matmul(kap3,adjYvYvadjkap3) 
 kap3adjYvYvCkap1 = Matmul(kap3,adjYvYvCkap1) 
 kap3adjYvYvCkap2 = Matmul(kap3,adjYvYvCkap2) 
 kap3adjYvYvCkap3 = Matmul(kap3,adjYvYvCkap3) 
 kap3adjkap3TpYvml2 = Matmul(kap3,adjkap3TpYvml2) 
 kap3Cmv2Ckap1kap1 = Matmul(kap3,Cmv2Ckap1kap1) 
 kap3Cmv2Ckap1kap2 = Matmul(kap3,Cmv2Ckap1kap2) 
 kap3Cmv2Ckap1kap3 = Matmul(kap3,Cmv2Ckap1kap3) 
 kap3Cmv2Ckap2kap1 = Matmul(kap3,Cmv2Ckap2kap1) 
 kap3Cmv2Ckap2kap2 = Matmul(kap3,Cmv2Ckap2kap2) 
 kap3Cmv2Ckap2kap3 = Matmul(kap3,Cmv2Ckap2kap3) 
 kap3Cmv2Ckap3kap1 = Matmul(kap3,Cmv2Ckap3kap1) 
 kap3Cmv2Ckap3kap2 = Matmul(kap3,Cmv2Ckap3kap2) 
 kap3Cmv2Ckap3kap3 = Matmul(kap3,Cmv2Ckap3kap3) 
 Tk1adjYvYvadjTk1 = Matmul(Tk1,adjYvYvadjTk1) 
 Tk1adjYvYvCkap1 = Matmul(Tk1,adjYvYvCkap1) 
 Tk2adjYvYvadjTk2 = Matmul(Tk2,adjYvYvadjTk2) 
 Tk2adjYvYvCkap2 = Matmul(Tk2,adjYvYvCkap2) 
 Tk3adjYvYvadjTk3 = Matmul(Tk3,adjYvYvadjTk3) 
 Tk3adjYvYvCkap3 = Matmul(Tk3,adjYvYvCkap3) 
 md2YdadjYdYdadjYd = Matmul(md2,YdadjYdYdadjYd) 
 md2YdadjYuYuadjYd = Matmul(md2,YdadjYuYuadjYd) 
 me2YeadjYeYeadjYe = Matmul(me2,YeadjYeYeadjYe) 
 me2YeCYvTpYvadjYe = Matmul(me2,YeCYvTpYvadjYe) 
 ml2adjYeYeadjYeYe = Matmul(ml2,adjYeYeadjYeYe) 
 ml2adjYeYeCYvlam = Matmul(ml2,adjYeYeCYvlam) 
 ml2CYvTpYvCYvlam = Matmul(ml2,CYvTpYvCYvlam) 
 ml2CYvTpYvCYvTpYv = Matmul(ml2,CYvTpYvCYvTpYv) 
 mq2adjYdYdadjYdYd = Matmul(mq2,adjYdYdadjYdYd) 
 mq2adjYdYdadjYuYu = Matmul(mq2,adjYdYdadjYuYu) 
 mq2adjYuYuadjYdYd = Matmul(mq2,adjYuYuadjYdYd) 
 mq2adjYuYuadjYuYu = Matmul(mq2,adjYuYuadjYuYu) 
 mu2YuadjYdYdadjYu = Matmul(mu2,YuadjYdYdadjYu) 
 mu2YuadjYuYuadjYu = Matmul(mu2,YuadjYuYuadjYu) 
 mv2TpYvadjYeYeCYv = Matmul(mv2,TpYvadjYeYeCYv) 
 mv2TpYvCYvTpYvCYv = Matmul(mv2,TpYvCYvTpYvCYv) 
 mv2kap1adjYvYvadjkap1 = Matmul(mv2,kap1adjYvYvadjkap1) 
 mv2kap2adjYvYvadjkap2 = Matmul(mv2,kap2adjYvYvadjkap2) 
 mv2kap3adjYvYvadjkap3 = Matmul(mv2,kap3adjYvYvadjkap3) 
 Ydmq2adjYdYdadjYd = Matmul(Yd,mq2adjYdYdadjYd) 
 Ydmq2adjYuYuadjYd = Matmul(Yd,mq2adjYuYuadjYd) 
 YdadjYdmd2YdadjYd = Matmul(Yd,adjYdmd2YdadjYd) 
Forall(i2=1:3)  YdadjYdmd2YdadjYd(i2,i2) =  Real(YdadjYdmd2YdadjYd(i2,i2),dp) 
 YdadjYdYdmq2adjYd = Matmul(Yd,adjYdYdmq2adjYd) 
 YdadjYdYdadjYdmd2 = Matmul(Yd,adjYdYdadjYdmd2) 
 YdadjYdYdadjYdYd = Matmul(Yd,adjYdYdadjYdYd) 
 YdadjYdYdadjYdTd = Matmul(Yd,adjYdYdadjYdTd) 
 YdadjYdTdadjYdYd = Matmul(Yd,adjYdTdadjYdYd) 
 YdadjYumu2YuadjYd = Matmul(Yd,adjYumu2YuadjYd) 
Forall(i2=1:3)  YdadjYumu2YuadjYd(i2,i2) =  Real(YdadjYumu2YuadjYd(i2,i2),dp) 
 YdadjYuYumq2adjYd = Matmul(Yd,adjYuYumq2adjYd) 
 YdadjYuYuadjYdmd2 = Matmul(Yd,adjYuYuadjYdmd2) 
 YdadjYuYuadjYdYd = Matmul(Yd,adjYuYuadjYdYd) 
 YdadjYuYuadjYdTd = Matmul(Yd,adjYuYuadjYdTd) 
 YdadjYuYuadjYuYu = Matmul(Yd,adjYuYuadjYuYu) 
 YdadjYuYuadjYuTu = Matmul(Yd,adjYuYuadjYuTu) 
 YdadjYuTuadjYdYd = Matmul(Yd,adjYuTuadjYdYd) 
 YdadjYuTuadjYuYu = Matmul(Yd,adjYuTuadjYuYu) 
 Yeml2adjYeYeadjYe = Matmul(Ye,ml2adjYeYeadjYe) 
 Yeml2CYvTpYvadjYe = Matmul(Ye,ml2CYvTpYvadjYe) 
 YeadjYeme2YeadjYe = Matmul(Ye,adjYeme2YeadjYe) 
Forall(i2=1:3)  YeadjYeme2YeadjYe(i2,i2) =  Real(YeadjYeme2YeadjYe(i2,i2),dp) 
 YeadjYeYeml2adjYe = Matmul(Ye,adjYeYeml2adjYe) 
 YeadjYeYeadjYeme2 = Matmul(Ye,adjYeYeadjYeme2) 
 YeadjYeYeadjYeYe = Matmul(Ye,adjYeYeadjYeYe) 
 YeadjYeYeadjYeTe = Matmul(Ye,adjYeYeadjYeTe) 
 YeadjYeTeadjYeYe = Matmul(Ye,adjYeTeadjYeYe) 
 YeCYvmv2TpYvadjYe = Matmul(Ye,CYvmv2TpYvadjYe) 
Forall(i2=1:3)  YeCYvmv2TpYvadjYe(i2,i2) =  Real(YeCYvmv2TpYvadjYe(i2,i2),dp) 
 YeCYvTpYvml2adjYe = Matmul(Ye,CYvTpYvml2adjYe) 
 YeCYvTpYvadjYeme2 = Matmul(Ye,CYvTpYvadjYeme2) 
 YeCYvTpYvadjYeYe = Matmul(Ye,CYvTpYvadjYeYe) 
 YeCYvTpYvadjYeTe = Matmul(Ye,CYvTpYvadjYeTe) 
 YeCYvTpYvCYvTpYv = Matmul(Ye,CYvTpYvCYvTpYv) 
 YeCYvTpYvCYvTpTv = Matmul(Ye,CYvTpYvCYvTpTv) 
 YeCYvTpTvadjYeYe = Matmul(Ye,CYvTpTvadjYeYe) 
 YeCYvTpTvCYvTpYv = Matmul(Ye,CYvTpTvCYvTpYv) 
 YeCYvkap1adjkap1TpYv = Matmul(Ye,CYvkap1adjkap1TpYv) 
 YeCYvkap2adjkap2TpYv = Matmul(Ye,CYvkap2adjkap2TpYv) 
 YeCYvkap3adjkap3TpYv = Matmul(Ye,CYvkap3adjkap3TpYv) 
 YeCYvTk1adjkap1TpYv = Matmul(Ye,CYvTk1adjkap1TpYv) 
 YeCYvTk2adjkap2TpYv = Matmul(Ye,CYvTk2adjkap2TpYv) 
 YeCYvTk3adjkap3TpYv = Matmul(Ye,CYvTk3adjkap3TpYv) 
 Yumq2adjYdYdadjYu = Matmul(Yu,mq2adjYdYdadjYu) 
 Yumq2adjYuYuadjYu = Matmul(Yu,mq2adjYuYuadjYu) 
 YuadjYdmd2YdadjYu = Matmul(Yu,adjYdmd2YdadjYu) 
Forall(i2=1:3)  YuadjYdmd2YdadjYu(i2,i2) =  Real(YuadjYdmd2YdadjYu(i2,i2),dp) 
 YuadjYdYdmq2adjYu = Matmul(Yu,adjYdYdmq2adjYu) 
 YuadjYdYdadjYdYd = Matmul(Yu,adjYdYdadjYdYd) 
 YuadjYdYdadjYdTd = Matmul(Yu,adjYdYdadjYdTd) 
 YuadjYdYdadjYumu2 = Matmul(Yu,adjYdYdadjYumu2) 
 YuadjYdYdadjYuYu = Matmul(Yu,adjYdYdadjYuYu) 
 YuadjYdYdadjYuTu = Matmul(Yu,adjYdYdadjYuTu) 
 YuadjYdTdadjYdYd = Matmul(Yu,adjYdTdadjYdYd) 
 YuadjYdTdadjYuYu = Matmul(Yu,adjYdTdadjYuYu) 
 YuadjYumu2YuadjYu = Matmul(Yu,adjYumu2YuadjYu) 
Forall(i2=1:3)  YuadjYumu2YuadjYu(i2,i2) =  Real(YuadjYumu2YuadjYu(i2,i2),dp) 
 YuadjYuYumq2adjYu = Matmul(Yu,adjYuYumq2adjYu) 
 YuadjYuYuadjYumu2 = Matmul(Yu,adjYuYuadjYumu2) 
 YuadjYuYuadjYuYu = Matmul(Yu,adjYuYuadjYuYu) 
 YuadjYuYuadjYuTu = Matmul(Yu,adjYuYuadjYuTu) 
 YuadjYuTuadjYuYu = Matmul(Yu,adjYuTuadjYuYu) 
 YvadjYvYvadjYvYv = Matmul(Yv,adjYvYvadjYvYv) 
 YvadjYvYvadjYvTv = Matmul(Yv,adjYvYvadjYvTv) 
 YvadjYvCml2YvadjYv = Matmul(Yv,adjYvCml2YvadjYv) 
Forall(i2=1:3)  YvadjYvCml2YvadjYv(i2,i2) =  Real(YvadjYvCml2YvadjYv(i2,i2),dp) 
 YvadjYvCml2TpYeCYe = Matmul(Yv,adjYvCml2TpYeCYe) 
 YvadjYvTvadjYvYv = Matmul(Yv,adjYvTvadjYvYv) 
 YvadjYvTpYeCme2CYe = Matmul(Yv,adjYvTpYeCme2CYe) 
 YvadjYvTpYeCYeYv = Matmul(Yv,adjYvTpYeCYeYv) 
 YvadjYvTpYeCYeCml2 = Matmul(Yv,adjYvTpYeCYeCml2) 
 YvadjYvTpYeCYeTv = Matmul(Yv,adjYvTpYeCYeTv) 
 YvadjYvTpTeCYeYv = Matmul(Yv,adjYvTpTeCYeYv) 
 YvCmv2adjYvYvadjYv = Matmul(Yv,Cmv2adjYvYvadjYv) 
 YvCmv2adjYvTpYeCYe = Matmul(Yv,Cmv2adjYvTpYeCYe) 
 YvCmv2Ckap1kap1adjYv = Matmul(Yv,Cmv2Ckap1kap1adjYv) 
 YvCmv2Ckap2kap2adjYv = Matmul(Yv,Cmv2Ckap2kap2adjYv) 
 YvCmv2Ckap3kap3adjYv = Matmul(Yv,Cmv2Ckap3kap3adjYv) 
 YvCkap1TpYvCYvTpTk1 = Matmul(Yv,Ckap1TpYvCYvTpTk1) 
 YvCkap1Tpkap1adjYvYv = Matmul(Yv,Ckap1Tpkap1adjYvYv) 
 YvCkap1Tpkap1Cmv2adjYv = Matmul(Yv,Ckap1Tpkap1Cmv2adjYv) 
 YvCkap1TpTk1adjYvYv = Matmul(Yv,Ckap1TpTk1adjYvYv) 
 YvCkap1kap1adjYvCml2 = Matmul(Yv,Ckap1kap1adjYvCml2) 
 YvCkap2TpYvCYvTpTk2 = Matmul(Yv,Ckap2TpYvCYvTpTk2) 
 YvCkap2Tpkap2adjYvYv = Matmul(Yv,Ckap2Tpkap2adjYvYv) 
 YvCkap2Tpkap2Cmv2adjYv = Matmul(Yv,Ckap2Tpkap2Cmv2adjYv) 
 YvCkap2TpTk2adjYvYv = Matmul(Yv,Ckap2TpTk2adjYvYv) 
 YvCkap2kap2adjYvCml2 = Matmul(Yv,Ckap2kap2adjYvCml2) 
 YvCkap3TpYvCYvTpTk3 = Matmul(Yv,Ckap3TpYvCYvTpTk3) 
 YvCkap3Tpkap3adjYvYv = Matmul(Yv,Ckap3Tpkap3adjYvYv) 
 YvCkap3Tpkap3Cmv2adjYv = Matmul(Yv,Ckap3Tpkap3Cmv2adjYv) 
 YvCkap3TpTk3adjYvYv = Matmul(Yv,Ckap3TpTk3adjYvYv) 
 YvCkap3kap3adjYvCml2 = Matmul(Yv,Ckap3kap3adjYvCml2) 
 adjYdmd2YdadjYdYd = Matmul(adjYd,md2YdadjYdYd) 
 adjYdYdmq2adjYdYd = Matmul(adjYd,Ydmq2adjYdYd) 
Forall(i2=1:3)  adjYdYdmq2adjYdYd(i2,i2) =  Real(adjYdYdmq2adjYdYd(i2,i2),dp) 
 adjYdYdadjYdmd2Yd = Matmul(adjYd,YdadjYdmd2Yd) 
 adjYdYdadjYdYdmq2 = Matmul(adjYd,YdadjYdYdmq2) 
 adjYeme2YeadjYeYe = Matmul(adjYe,me2YeadjYeYe) 
 adjYeme2YeCYvlam = Matmul(adjYe,me2YeCYvlam) 
 adjYeYeml2adjYeYe = Matmul(adjYe,Yeml2adjYeYe) 
Forall(i2=1:3)  adjYeYeml2adjYeYe(i2,i2) =  Real(adjYeYeml2adjYeYe(i2,i2),dp) 
 adjYeYeml2CYvlam = Matmul(adjYe,Yeml2CYvlam) 
 adjYeYeadjYeme2Ye = Matmul(adjYe,YeadjYeme2Ye) 
 adjYeYeadjYeYeml2 = Matmul(adjYe,YeadjYeYeml2) 
 adjYeYeCYvmv2lam = Matmul(adjYe,YeCYvmv2lam) 
 adjYumu2YuadjYuYu = Matmul(adjYu,mu2YuadjYuYu) 
 adjYuYumq2adjYuYu = Matmul(adjYu,Yumq2adjYuYu) 
Forall(i2=1:3)  adjYuYumq2adjYuYu(i2,i2) =  Real(adjYuYumq2adjYuYu(i2,i2),dp) 
 adjYuYuadjYumu2Yu = Matmul(adjYu,YuadjYumu2Yu) 
 adjYuYuadjYuYumq2 = Matmul(adjYu,YuadjYuYumq2) 
 adjYvYvadjkap1mv2kap2 = Matmul(adjYv,Yvadjkap1mv2kap2) 
 adjYvYvadjkap1mv2kap3 = Matmul(adjYv,Yvadjkap1mv2kap3) 
 adjYvYvadjkap2mv2kap1 = Matmul(adjYv,Yvadjkap2mv2kap1) 
 adjYvYvadjkap2mv2kap3 = Matmul(adjYv,Yvadjkap2mv2kap3) 
 adjYvYvadjkap3mv2kap1 = Matmul(adjYv,Yvadjkap3mv2kap1) 
 adjYvYvadjkap3mv2kap2 = Matmul(adjYv,Yvadjkap3mv2kap2) 
 adjYvYvCmv2Ckap1kap2 = Matmul(adjYv,YvCmv2Ckap1kap2) 
 adjYvYvCmv2Ckap1kap3 = Matmul(adjYv,YvCmv2Ckap1kap3) 
 adjYvYvCmv2Ckap2kap1 = Matmul(adjYv,YvCmv2Ckap2kap1) 
 adjYvYvCmv2Ckap2kap3 = Matmul(adjYv,YvCmv2Ckap2kap3) 
 adjYvYvCmv2Ckap3kap1 = Matmul(adjYv,YvCmv2Ckap3kap1) 
 adjYvYvCmv2Ckap3kap2 = Matmul(adjYv,YvCmv2Ckap3kap2) 
 adjYvCml2YvCkap1kap2 = Matmul(adjYv,Cml2YvCkap1kap2) 
 adjYvCml2YvCkap1kap3 = Matmul(adjYv,Cml2YvCkap1kap3) 
 adjYvCml2YvCkap2kap1 = Matmul(adjYv,Cml2YvCkap2kap1) 
 adjYvCml2YvCkap2kap3 = Matmul(adjYv,Cml2YvCkap2kap3) 
 adjYvCml2YvCkap3kap1 = Matmul(adjYv,Cml2YvCkap3kap1) 
 adjYvCml2YvCkap3kap2 = Matmul(adjYv,Cml2YvCkap3kap2) 
 CYvmv2TpYvCYvlam = Matmul(Conjg(Yv),mv2TpYvCYvlam) 
 CYvmv2TpYvCYvTpYv = Matmul(Conjg(Yv),mv2TpYvCYvTpYv) 
 CYvmv2kap1adjkap1TpYv = Matmul(Conjg(Yv),mv2kap1adjkap1TpYv) 
 CYvmv2kap2adjkap2TpYv = Matmul(Conjg(Yv),mv2kap2adjkap2TpYv) 
 CYvmv2kap3adjkap3TpYv = Matmul(Conjg(Yv),mv2kap3adjkap3TpYv) 
 CYvTpYvml2CYvlam = Matmul(Conjg(Yv),TpYvml2CYvlam) 
 CYvTpYvml2CYvTpYv = Matmul(Conjg(Yv),TpYvml2CYvTpYv) 
Forall(i2=1:3)  CYvTpYvml2CYvTpYv(i2,i2) =  Real(CYvTpYvml2CYvTpYv(i2,i2),dp) 
 CYvTpYvCYvmv2lam = Matmul(Conjg(Yv),TpYvCYvmv2lam) 
 CYvTpYvCYvmv2TpYv = Matmul(Conjg(Yv),TpYvCYvmv2TpYv) 
 CYvTpYvCYvTpYvml2 = Matmul(Conjg(Yv),TpYvCYvTpYvml2) 
 CYvkap1adjkap1TpYvml2 = Matmul(Conjg(Yv),kap1adjkap1TpYvml2) 
 CYvkap2adjkap2TpYvml2 = Matmul(Conjg(Yv),kap2adjkap2TpYvml2) 
 CYvkap3adjkap3TpYvml2 = Matmul(Conjg(Yv),kap3adjkap3TpYvml2) 
 TdadjYdYdadjYdYd = Matmul(Td,adjYdYdadjYdYd) 
 TdadjYuYuadjYdYd = Matmul(Td,adjYuYuadjYdYd) 
 TdadjYuYuadjYuYu = Matmul(Td,adjYuYuadjYuYu) 
 TeadjYeYeadjYeYe = Matmul(Te,adjYeYeadjYeYe) 
 TeCYvTpYvadjYeYe = Matmul(Te,CYvTpYvadjYeYe) 
 TeCYvTpYvCYvTpYv = Matmul(Te,CYvTpYvCYvTpYv) 
 TuadjYdYdadjYdYd = Matmul(Tu,adjYdYdadjYdYd) 
 TuadjYdYdadjYuYu = Matmul(Tu,adjYdYdadjYuYu) 
 TuadjYuYuadjYuYu = Matmul(Tu,adjYuYuadjYuYu) 
 TvadjYvYvadjYvYv = Matmul(Tv,adjYvYvadjYvYv) 
 TvadjYvTpYeCYeYv = Matmul(Tv,adjYvTpYeCYeYv) 
 TpYeCYeTpYeCYeYv = Matmul(Transpose(Ye),CYeTpYeCYeYv) 
 TpYeCYeTpYeCYeTv = Matmul(Transpose(Ye),CYeTpYeCYeTv) 
 TpYeCYeTpTeCYeYv = Matmul(Transpose(Ye),CYeTpTeCYeYv) 
 TpYvml2adjYeYeCYv = Matmul(Transpose(Yv),ml2adjYeYeCYv) 
 TpYvml2CYvTpYvCYv = Matmul(Transpose(Yv),ml2CYvTpYvCYv) 
 TpYvadjYeme2YeCYv = Matmul(Transpose(Yv),adjYeme2YeCYv) 
Forall(i2=1:3)  TpYvadjYeme2YeCYv(i2,i2) =  Real(TpYvadjYeme2YeCYv(i2,i2),dp) 
 TpYvadjYeYeml2CYv = Matmul(Transpose(Yv),adjYeYeml2CYv) 
 TpYvadjYeYeCYvmv2 = Matmul(Transpose(Yv),adjYeYeCYvmv2) 
 TpYvadjYeYeCYvlam = Matmul(Transpose(Yv),adjYeYeCYvlam) 
 TpYvadjYeTeCYvlam = Matmul(Transpose(Yv),adjYeTeCYvlam) 
 TpYvCYvmv2TpYvCYv = Matmul(Transpose(Yv),CYvmv2TpYvCYv) 
Forall(i2=1:3)  TpYvCYvmv2TpYvCYv(i2,i2) =  Real(TpYvCYvmv2TpYvCYv(i2,i2),dp) 
 TpYvCYvTpYvml2CYv = Matmul(Transpose(Yv),CYvTpYvml2CYv) 
 TpYvCYvTpYvCYvmv2 = Matmul(Transpose(Yv),CYvTpYvCYvmv2) 
 TpYvCYvTpYvCYvlam = Matmul(Transpose(Yv),CYvTpYvCYvlam) 
 TpYvCYvTpYvCYvTlam = Matmul(Transpose(Yv),CYvTpYvCYvTlam) 
 TpYvCYvTpTvCYvlam = Matmul(Transpose(Yv),CYvTpTvCYvlam) 
 TpYvCYvkap1adjkap1lam = Matmul(Transpose(Yv),CYvkap1adjkap1lam) 
 TpYvCYvkap2adjkap2lam = Matmul(Transpose(Yv),CYvkap2adjkap2lam) 
 TpYvCYvkap3adjkap3lam = Matmul(Transpose(Yv),CYvkap3adjkap3lam) 
 TpYvCYvTk1adjkap1lam = Matmul(Transpose(Yv),CYvTk1adjkap1lam) 
 TpYvCYvTk2adjkap2lam = Matmul(Transpose(Yv),CYvTk2adjkap2lam) 
 TpYvCYvTk3adjkap3lam = Matmul(Transpose(Yv),CYvTk3adjkap3lam) 
 TpTeCYeTpYeCYeYv = Matmul(Transpose(Te),CYeTpYeCYeYv) 
 TpTvadjYeYeCYvlam = Matmul(Transpose(Tv),adjYeYeCYvlam) 
 TpTvCYvTpYvCYvlam = Matmul(Transpose(Tv),CYvTpYvCYvlam) 
 kap1adjYvYvadjkap1mv2 = Matmul(kap1,adjYvYvadjkap1mv2) 
 kap1adjYvYvCmv2adjkap1 = Matmul(kap1,adjYvYvCmv2adjkap1) 
 kap1Cmv2adjYvYvadjkap1 = Matmul(kap1,Cmv2adjYvYvadjkap1) 
 kap1Ckap1Tpkap1Cmv2adjkap1 = Matmul(kap1,Ckap1Tpkap1Cmv2adjkap1) 
 kap1Ckap2Tpkap2Cmv2adjkap1 = Matmul(kap1,Ckap2Tpkap2Cmv2adjkap1) 
 kap1Ckap3Tpkap3Cmv2adjkap1 = Matmul(kap1,Ckap3Tpkap3Cmv2adjkap1) 
 kap2adjYvYvadjkap2mv2 = Matmul(kap2,adjYvYvadjkap2mv2) 
 kap2adjYvYvCmv2adjkap2 = Matmul(kap2,adjYvYvCmv2adjkap2) 
 kap2Cmv2adjYvYvadjkap2 = Matmul(kap2,Cmv2adjYvYvadjkap2) 
 kap2Ckap1Tpkap1Cmv2adjkap2 = Matmul(kap2,Ckap1Tpkap1Cmv2adjkap2) 
 kap2Ckap2Tpkap2Cmv2adjkap2 = Matmul(kap2,Ckap2Tpkap2Cmv2adjkap2) 
 kap2Ckap3Tpkap3Cmv2adjkap2 = Matmul(kap2,Ckap3Tpkap3Cmv2adjkap2) 
 kap3adjYvYvadjkap3mv2 = Matmul(kap3,adjYvYvadjkap3mv2) 
 kap3adjYvYvCmv2adjkap3 = Matmul(kap3,adjYvYvCmv2adjkap3) 
 kap3Cmv2adjYvYvadjkap3 = Matmul(kap3,Cmv2adjYvYvadjkap3) 
 kap3Ckap1Tpkap1Cmv2adjkap3 = Matmul(kap3,Ckap1Tpkap1Cmv2adjkap3) 
 kap3Ckap2Tpkap2Cmv2adjkap3 = Matmul(kap3,Ckap2Tpkap2Cmv2adjkap3) 
 kap3Ckap3Tpkap3Cmv2adjkap3 = Matmul(kap3,Ckap3Tpkap3Cmv2adjkap3) 
 Tk1adjYvYvCkap1lam = Matmul(Tk1,adjYvYvCkap1lam) 
 Tk2adjYvYvCkap2lam = Matmul(Tk2,adjYvYvCkap2lam) 
 Tk3adjYvYvCkap3lam = Matmul(Tk3,adjYvYvCkap3lam) 
 TradjTk1Tk1 = cTrace(adjTk1Tk1) 
 TradjTk2Tk2 = cTrace(adjTk2Tk2) 
 TradjTk3Tk3 = cTrace(adjTk3Tk3) 
 TrCTdTpYd = cTrace(CTdTpYd) 
 TrCTeTpYe = cTrace(CTeTpYe) 
 TrCTuTpYu = cTrace(CTuTpYu) 
 TrCTvTpYv = cTrace(CTvTpYv) 
 TrCkap1Tpkap1 = cTrace(Ckap1Tpkap1) 
 TrCkap1Tpkap2 = cTrace(Ckap1Tpkap2) 
 TrCkap1Tpkap3 = cTrace(Ckap1Tpkap3) 
 TrCkap1kap1 = cTrace(Ckap1kap1) 
 TrCkap1kap2 = cTrace(Ckap1kap2) 
 TrCkap1kap3 = cTrace(Ckap1kap3) 
 TrCkap1Tk1 = cTrace(Ckap1Tk1) 
 TrCkap1Tk2 = cTrace(Ckap1Tk2) 
 TrCkap1Tk3 = cTrace(Ckap1Tk3) 
 TrCkap2Tpkap1 = cTrace(Ckap2Tpkap1) 
 TrCkap2Tpkap2 = cTrace(Ckap2Tpkap2) 
 TrCkap2Tpkap3 = cTrace(Ckap2Tpkap3) 
 TrCkap2kap1 = cTrace(Ckap2kap1) 
 TrCkap2kap2 = cTrace(Ckap2kap2) 
 TrCkap2kap3 = cTrace(Ckap2kap3) 
 TrCkap2Tk1 = cTrace(Ckap2Tk1) 
 TrCkap2Tk2 = cTrace(Ckap2Tk2) 
 TrCkap2Tk3 = cTrace(Ckap2Tk3) 
 TrCkap3Tpkap1 = cTrace(Ckap3Tpkap1) 
 TrCkap3Tpkap2 = cTrace(Ckap3Tpkap2) 
 TrCkap3Tpkap3 = cTrace(Ckap3Tpkap3) 
 TrCkap3kap1 = cTrace(Ckap3kap1) 
 TrCkap3kap2 = cTrace(Ckap3kap2) 
 TrCkap3kap3 = cTrace(Ckap3kap3) 
 TrCkap3Tk1 = cTrace(Ckap3Tk1) 
 TrCkap3Tk2 = cTrace(Ckap3Tk2) 
 TrCkap3Tk3 = cTrace(Ckap3Tk3) 
 TrCTk1Tpkap1 = cTrace(CTk1Tpkap1) 
 TrCTk1Tpkap2 = cTrace(CTk1Tpkap2) 
 TrCTk1Tpkap3 = cTrace(CTk1Tpkap3) 
 TrCTk2Tpkap1 = cTrace(CTk2Tpkap1) 
 TrCTk2Tpkap2 = cTrace(CTk2Tpkap2) 
 TrCTk2Tpkap3 = cTrace(CTk2Tpkap3) 
 TrCTk3Tpkap1 = cTrace(CTk3Tpkap1) 
 TrCTk3Tpkap2 = cTrace(CTk3Tpkap2) 
 TrCTk3Tpkap3 = cTrace(CTk3Tpkap3) 
 Trml2YvadjYv = cTrace(ml2YvadjYv) 
 Trmv2kap1adjkap1 = cTrace(mv2kap1adjkap1) 
 Trmv2kap1adjkap2 = cTrace(mv2kap1adjkap2) 
 Trmv2kap1adjkap3 = cTrace(mv2kap1adjkap3) 
 Trmv2kap2adjkap1 = cTrace(mv2kap2adjkap1) 
 Trmv2kap2adjkap2 = cTrace(mv2kap2adjkap2) 
 Trmv2kap2adjkap3 = cTrace(mv2kap2adjkap3) 
 Trmv2kap3adjkap1 = cTrace(mv2kap3adjkap1) 
 Trmv2kap3adjkap2 = cTrace(mv2kap3adjkap2) 
 Trmv2kap3adjkap3 = cTrace(mv2kap3adjkap3) 
 TrYdadjYdCmd2 = cTrace(YdadjYdCmd2) 
 TrYdCmq2adjYd = cTrace(YdCmq2adjYd) 
 TrYeadjYeCme2 = cTrace(YeadjYeCme2) 
 TrYeCml2adjYe = cTrace(YeCml2adjYe) 
 TrYuadjYuCmu2 = cTrace(YuadjYuCmu2) 
 TrYuCmq2adjYu = cTrace(YuCmq2adjYu) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdTdadjYd = cTrace(YdadjYdTdadjYd) 
 TrYdadjYdTdadjTd = cTrace(YdadjYdTdadjTd) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYdadjYuTuadjYd = cTrace(YdadjYuTuadjYd) 
 TrYdadjYuTuadjTd = cTrace(YdadjYuTuadjTd) 
 TrYdadjTdTdadjYd = cTrace(YdadjTdTdadjYd) 
 TrYdadjTuTuadjYd = cTrace(YdadjTuTuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYeTeadjYe = cTrace(YeadjYeTeadjYe) 
 TrYeadjYeTeadjTe = cTrace(YeadjYeTeadjTe) 
 TrYeadjTeTeadjYe = cTrace(YeadjTeTeadjYe) 
 TrYeCTvTpTvadjYe = cTrace(YeCTvTpTvadjYe) 
 TrYuadjYdTdadjYu = cTrace(YuadjYdTdadjYu) 
 TrYuadjYdTdadjTu = cTrace(YuadjYdTdadjTu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYuadjYuTuadjYu = cTrace(YuadjYuTuadjYu) 
 TrYuadjYuTuadjTu = cTrace(YuadjYuTuadjTu) 
 TrYuadjTdTdadjYu = cTrace(YuadjTdTdadjYu) 
 TrYuadjTuTuadjYu = cTrace(YuadjTuTuadjYu) 
 TrYvadjYvYvadjYv = cTrace(YvadjYvYvadjYv) 
 TrYvadjYvTvadjYv = cTrace(YvadjYvTvadjYv) 
 TrYvadjYvTvadjTv = cTrace(YvadjYvTvadjTv) 
 TrYvadjYvTpYeCYe = cTrace(YvadjYvTpYeCYe) 
 TrYvadjYvTpTeCTe = cTrace(YvadjYvTpTeCTe) 
 TrYvadjTvTvadjYv = cTrace(YvadjTvTvadjYv) 
 TrYvCkap1Tpkap1adjYv = cTrace(YvCkap1Tpkap1adjYv) 
 TrYvCkap1TpTk1adjYv = cTrace(YvCkap1TpTk1adjYv) 
 TrYvCkap1kap1adjYv = cTrace(YvCkap1kap1adjYv) 
 TrYvCkap1Tk1adjTv = cTrace(YvCkap1Tk1adjTv) 
 TrYvCkap2Tpkap2adjYv = cTrace(YvCkap2Tpkap2adjYv) 
 TrYvCkap2TpTk2adjYv = cTrace(YvCkap2TpTk2adjYv) 
 TrYvCkap2kap2adjYv = cTrace(YvCkap2kap2adjYv) 
 TrYvCkap2Tk2adjTv = cTrace(YvCkap2Tk2adjTv) 
 TrYvCkap3Tpkap3adjYv = cTrace(YvCkap3Tpkap3adjYv) 
 TrYvCkap3TpTk3adjYv = cTrace(YvCkap3TpTk3adjYv) 
 TrYvCkap3kap3adjYv = cTrace(YvCkap3kap3adjYv) 
 TrYvCkap3Tk3adjTv = cTrace(YvCkap3Tk3adjTv) 
 TrYvCTk1TpTk1adjYv = cTrace(YvCTk1TpTk1adjYv) 
 TrYvCTk2TpTk2adjYv = cTrace(YvCTk2TpTk2adjYv) 
 TrYvCTk3TpTk3adjYv = cTrace(YvCTk3TpTk3adjYv) 
 TradjYeTeCYvTpYv = cTrace(adjYeTeCYvTpYv) 
 TradjYeTeCTvTpYv = cTrace(adjYeTeCTvTpYv) 
 TradjYvTvadjkap1kap1 = cTrace(adjYvTvadjkap1kap1) 
 TradjYvTvadjkap2kap2 = cTrace(adjYvTvadjkap2kap2) 
 TradjYvTvadjkap3kap3 = cTrace(adjYvTvadjkap3kap3) 
 TradjYvTpYeCYeTv = cTrace(adjYvTpYeCYeTv) 
 TradjYvTpYeCTeTv = cTrace(adjYvTpYeCTeTv) 
 Trmd2YdadjYdYdadjYd = cTrace(md2YdadjYdYdadjYd) 
 Trmd2YdadjYuYuadjYd = cTrace(md2YdadjYuYuadjYd) 
 Trme2YeadjYeYeadjYe = cTrace(me2YeadjYeYeadjYe) 
 Trml2adjYeYeadjYeYe = cTrace(ml2adjYeYeadjYeYe) 
 Trmq2adjYdYdadjYdYd = cTrace(mq2adjYdYdadjYdYd) 
 Trmq2adjYdYdadjYuYu = cTrace(mq2adjYdYdadjYuYu) 
 Trmq2adjYuYuadjYdYd = cTrace(mq2adjYuYuadjYdYd) 
 Trmq2adjYuYuadjYuYu = cTrace(mq2adjYuYuadjYuYu) 
 Trmu2YuadjYdYdadjYu = cTrace(mu2YuadjYdYdadjYu) 
 Trmu2YuadjYuYuadjYu = cTrace(mu2YuadjYuYuadjYu) 
 Trmv2kap1adjYvYvadjkap1 = cTrace(mv2kap1adjYvYvadjkap1) 
 Trmv2kap2adjYvYvadjkap2 = cTrace(mv2kap2adjYvYvadjkap2) 
 Trmv2kap3adjYvYvadjkap3 = cTrace(mv2kap3adjYvYvadjkap3) 
 TrYvadjYvCml2YvadjYv = cTrace(YvadjYvCml2YvadjYv) 
 TrYvadjYvCml2TpYeCYe = cTrace(YvadjYvCml2TpYeCYe) 
 TrYvadjYvTpYeCme2CYe = cTrace(YvadjYvTpYeCme2CYe) 
 TrYvadjYvTpYeCYeCml2 = cTrace(YvadjYvTpYeCYeCml2) 
 TrYvCmv2adjYvYvadjYv = cTrace(YvCmv2adjYvYvadjYv) 
 TrYvCmv2adjYvTpYeCYe = cTrace(YvCmv2adjYvTpYeCYe) 
 TrYvCmv2Ckap1kap1adjYv = cTrace(YvCmv2Ckap1kap1adjYv) 
 TrYvCmv2Ckap2kap2adjYv = cTrace(YvCmv2Ckap2kap2adjYv) 
 TrYvCmv2Ckap3kap3adjYv = cTrace(YvCmv2Ckap3kap3adjYv) 
 TrYvCkap1Tpkap1Cmv2adjYv = cTrace(YvCkap1Tpkap1Cmv2adjYv) 
 TrYvCkap1kap1adjYvCml2 = cTrace(YvCkap1kap1adjYvCml2) 
 TrYvCkap2Tpkap2Cmv2adjYv = cTrace(YvCkap2Tpkap2Cmv2adjYv) 
 TrYvCkap2kap2adjYvCml2 = cTrace(YvCkap2kap2adjYvCml2) 
 TrYvCkap3Tpkap3Cmv2adjYv = cTrace(YvCkap3Tpkap3Cmv2adjYv) 
 TrYvCkap3kap3adjYvCml2 = cTrace(YvCkap3kap3adjYvCml2) 
 SPmlHd2xxadjYeYeCYvlam = DOT_PRODUCT(mlHd2,adjYeYeCYvlam) 
 SPlamxxCTlam = DOT_PRODUCT(Conjg(lam),Conjg(Tlam)) 
 SPlamxxadjkap1lam = DOT_PRODUCT(Conjg(lam),adjkap1lam) 
 SPlamxxadjkap2lam = DOT_PRODUCT(Conjg(lam),adjkap2lam) 
 SPlamxxadjkap3lam = DOT_PRODUCT(Conjg(lam),adjkap3lam) 
 SPlamxxadjYvYvClam = DOT_PRODUCT(Conjg(lam),adjYvYvClam) 
 SPlamxxadjTvTvClam = DOT_PRODUCT(Conjg(lam),adjTvTvClam) 
 SPlamxxadjYvYvadjYvmlHd2 = DOT_PRODUCT(Conjg(lam),adjYvYvadjYvmlHd2) 
 SPlamxxadjYvYvCmv2Clam = DOT_PRODUCT(Conjg(lam),adjYvYvCmv2Clam) 
 SPlamxxCmv2adjYvYvClam = DOT_PRODUCT(Conjg(lam),Cmv2adjYvYvClam) 
 SPClamxxTpYvmlHd2 = DOT_PRODUCT(Conjg(Conjg(lam)),TpYvmlHd2) 
 SPClamxxTpTvCYvlam = DOT_PRODUCT(Conjg(Conjg(lam)),TpTvCYvlam) 
 SPClamxxTpYvml2CYvlam = DOT_PRODUCT(Conjg(Conjg(lam)),TpYvml2CYvlam) 
 SPTpkap1Clamxxadjkap1lam = DOT_PRODUCT(Conjg(Tpkap1Clam),adjkap1lam) 
 SPTpkap2Clamxxadjkap2lam = DOT_PRODUCT(Conjg(Tpkap2Clam),adjkap2lam) 
 SPTpkap3Clamxxadjkap3lam = DOT_PRODUCT(Conjg(Tpkap3Clam),adjkap3lam) 
 SPTpTk1ClamxxadjTk1lam = DOT_PRODUCT(Conjg(TpTk1Clam),adjTk1lam) 
 SPTpTk2ClamxxadjTk2lam = DOT_PRODUCT(Conjg(TpTk2Clam),adjTk2lam) 
 SPTpTk3ClamxxadjTk3lam = DOT_PRODUCT(Conjg(TpTk3Clam),adjTk3lam) 
 SPadjYvYvClamxxTlam = DOT_PRODUCT(Conjg(adjYvYvClam),Tlam) 
 SPadjYvYvCTlamxxTlam = DOT_PRODUCT(Conjg(adjYvYvCTlam),Tlam) 
 SPadjTvYvClamxxTlam = DOT_PRODUCT(Conjg(adjTvYvClam),Tlam) 
 SPTpTvCYvlamxxCTlam = DOT_PRODUCT(Conjg(TpTvCYvlam),Conjg(Tlam)) 
 g1p4 =g1**4 
 g1p5 =g1**5 
 g2p4 =g2**4 
 g2p5 =g2**5 
 g3p4 =g3**4 
 g3p5 =g3**5 
 SPlamxxClamp2 =SPlamxxClam**2 
Do i1=1,3
  Do i2=1,3
Dylami1CTlami2(i1,i2) = Conjg(Tlam(i2))*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1YvClami2(i1,i2) = CYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTlami1YvCTlami2(i1,i2) = CYvTlam(i1)*YvCTlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvTlami1YvClami2(i1,i2) = CTvTlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvlami1TvClami2(i1,i2) = CTvlam(i1)*TvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1TvCTlami2(i1,i2) = CYvlam(i1)*TvCTlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Clami2(i1,i2) = kap1Clam(i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Clami2(i1,i2) = kap1Clam(i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Clami2(i1,i2) = kap1Clam(i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Clami2(i1,i2) = kap2Clam(i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Clami2(i1,i2) = kap2Clam(i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Clami2(i1,i2) = kap2Clam(i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Clami2(i1,i2) = kap3Clam(i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Clami2(i1,i2) = kap3Clam(i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Clami2(i1,i2) = kap3Clam(i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1Clami1CTk1lami2(i1,i2) = CTk1lam(i2)*Tk1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2Clami1CTk2lami2(i1,i2) = CTk2lam(i2)*Tk2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3Clami1CTk3lami2(i1,i2) = CTk3lam(i2)*Tk3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvlami1YvClami2(i1,i2) = ml2CYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1Clami1Ckap1lami2(i1,i2) = Ckap1lam(i2)*mv2kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2Clami1Ckap2lami2(i1,i2) = Ckap2lam(i2)*mv2kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3Clami1Ckap3lami2(i1,i2) = Ckap3lam(i2)*mv2kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvlami1YvClami2(i1,i2) = YeCYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvlami1TvClami2(i1,i2) = TvClam(i2)*YeCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvTlami1YvClami2(i1,i2) = YeCYvTlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1YvadjYvmlHd2i2(i1,i2) = CYvlam(i1)*YvadjYvmlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap1i11adjYvTv1i2(i1,i2) = adjYvTv(1,i2)*Yvadjkap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap1kap1i21(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi11Yvadjkap1kap1i21(i1,i2) = ml2CYv(i1,1)*Yvadjkap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi11Yvadjkap1kap1i21(i1,i2) = TeCYv(i1,1)*Yvadjkap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap2i11adjYvTv2i2(i1,i2) = adjYvTv(2,i2)*Yvadjkap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap1kap2i21(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi12Yvadjkap1kap2i21(i1,i2) = ml2CYv(i1,2)*Yvadjkap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi12Yvadjkap1kap2i21(i1,i2) = TeCYv(i1,2)*Yvadjkap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap3i11adjYvTv3i2(i1,i2) = adjYvTv(3,i2)*Yvadjkap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap1kap3i21(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi13Yvadjkap1kap3i21(i1,i2) = ml2CYv(i1,3)*Yvadjkap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi13Yvadjkap1kap3i21(i1,i2) = TeCYv(i1,3)*Yvadjkap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1Tk1i11lami2(i1,i2) = Yvadjkap1Tk1(i1,1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1Tk2i11lami2(i1,i2) = Yvadjkap1Tk2(i1,1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1Tk3i11lami2(i1,i2) = Yvadjkap1Tk3(i1,1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap1i12adjYvTv1i2(i1,i2) = adjYvTv(1,i2)*Yvadjkap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap2kap1i22(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi11Yvadjkap2kap1i22(i1,i2) = ml2CYv(i1,1)*Yvadjkap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi11Yvadjkap2kap1i22(i1,i2) = TeCYv(i1,1)*Yvadjkap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap2i12adjYvTv2i2(i1,i2) = adjYvTv(2,i2)*Yvadjkap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap2kap2i22(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi12Yvadjkap2kap2i22(i1,i2) = ml2CYv(i1,2)*Yvadjkap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi12Yvadjkap2kap2i22(i1,i2) = TeCYv(i1,2)*Yvadjkap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap3i12adjYvTv3i2(i1,i2) = adjYvTv(3,i2)*Yvadjkap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap2kap3i22(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi13Yvadjkap2kap3i22(i1,i2) = ml2CYv(i1,3)*Yvadjkap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi13Yvadjkap2kap3i22(i1,i2) = TeCYv(i1,3)*Yvadjkap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2Tk1i12lami2(i1,i2) = Yvadjkap2Tk1(i1,2)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2Tk2i12lami2(i1,i2) = Yvadjkap2Tk2(i1,2)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2Tk3i12lami2(i1,i2) = Yvadjkap2Tk3(i1,2)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap1i13adjYvTv1i2(i1,i2) = adjYvTv(1,i2)*Yvadjkap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap3kap1i23(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi11Yvadjkap3kap1i23(i1,i2) = ml2CYv(i1,1)*Yvadjkap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi11Yvadjkap3kap1i23(i1,i2) = TeCYv(i1,1)*Yvadjkap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap2i13adjYvTv2i2(i1,i2) = adjYvTv(2,i2)*Yvadjkap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap3kap2i23(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi12Yvadjkap3kap2i23(i1,i2) = ml2CYv(i1,2)*Yvadjkap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi12Yvadjkap3kap2i23(i1,i2) = TeCYv(i1,2)*Yvadjkap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap3i13adjYvTv3i2(i1,i2) = adjYvTv(3,i2)*Yvadjkap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap3kap3i23(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi13Yvadjkap3kap3i23(i1,i2) = ml2CYv(i1,3)*Yvadjkap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi13Yvadjkap3kap3i23(i1,i2) = TeCYv(i1,3)*Yvadjkap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3Tk1i13lami2(i1,i2) = Yvadjkap3Tk1(i1,3)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3Tk2i13lami2(i1,i2) = Yvadjkap3Tk2(i1,3)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3Tk3i13lami2(i1,i2) = Yvadjkap3Tk3(i1,3)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1YvCmv2Clami2(i1,i2) = CYvlam(i1)*YvCmv2Clam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1lami1Tk1Clami2(i1,i2) = Tk1Clam(i2)*YvCkap1lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11YvCkap1Tk1i21(i1,i2) = Conjg(Tv(i1,1))*YvCkap1Tk1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12YvCkap1Tk2i21(i1,i2) = Conjg(Tv(i1,2))*YvCkap1Tk2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13YvCkap1Tk3i21(i1,i2) = Conjg(Tv(i1,3))*YvCkap1Tk3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2lami1Tk2Clami2(i1,i2) = Tk2Clam(i2)*YvCkap2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11YvCkap2Tk1i22(i1,i2) = Conjg(Tv(i1,1))*YvCkap2Tk1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12YvCkap2Tk2i22(i1,i2) = Conjg(Tv(i1,2))*YvCkap2Tk2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13YvCkap2Tk3i22(i1,i2) = Conjg(Tv(i1,3))*YvCkap2Tk3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3lami1Tk3Clami2(i1,i2) = Tk3Clam(i2)*YvCkap3lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11YvCkap3Tk1i23(i1,i2) = Conjg(Tv(i1,1))*YvCkap3Tk1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12YvCkap3Tk2i23(i1,i2) = Conjg(Tv(i1,2))*YvCkap3Tk2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13YvCkap3Tk3i23(i1,i2) = Conjg(Tv(i1,3))*YvCkap3Tk3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvYvClami2(i1,i2) = adjYvYvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2lami1adjYvYvClami2(i1,i2) = adjYvYvClam(i2)*mv2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvTvCTlami2(i1,i2) = adjYvTvCTlam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjTvTvClami2(i1,i2) = adjTvTvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1Cml2YvClami2(i1,i2) = Cml2YvClam(i2)*CYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvmv2lami1YvClami2(i1,i2) = CYvmv2lam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap1Ckap1i11mlHd2i2(i1,i2) = CYvTpkap1Ckap1(i1,1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap1Ckap2i11mlHd2i2(i1,i2) = CYvTpkap1Ckap2(i1,1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap1Ckap3i11mlHd2i2(i1,i2) = CYvTpkap1Ckap3(i1,1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap2Ckap1i12mlHd2i2(i1,i2) = CYvTpkap2Ckap1(i1,2)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap2Ckap2i12mlHd2i2(i1,i2) = CYvTpkap2Ckap2(i1,2)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap2Ckap3i12mlHd2i2(i1,i2) = CYvTpkap2Ckap3(i1,2)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap3Ckap1i13mlHd2i2(i1,i2) = CYvTpkap3Ckap1(i1,3)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap3Ckap2i13mlHd2i2(i1,i2) = CYvTpkap3Ckap2(i1,3)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap3Ckap3i13mlHd2i2(i1,i2) = CYvTpkap3Ckap3(i1,3)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap1mv2lami2(i1,i2) = Ckap1mv2lam(i2)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i11Ckap1kap1Ckap1i21(i1,i2) = Ckap1kap1Ckap1(i2,1)*mv2kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i12Ckap1kap1Ckap2i21(i1,i2) = Ckap1kap1Ckap2(i2,1)*mv2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i13Ckap1kap1Ckap3i21(i1,i2) = Ckap1kap1Ckap3(i2,1)*mv2kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i11Ckap1kap2Ckap1i22(i1,i2) = Ckap1kap2Ckap1(i2,2)*mv2kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i12Ckap1kap2Ckap2i22(i1,i2) = Ckap1kap2Ckap2(i2,2)*mv2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i13Ckap1kap2Ckap3i22(i1,i2) = Ckap1kap2Ckap3(i2,2)*mv2kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i11Ckap1kap3Ckap1i23(i1,i2) = Ckap1kap3Ckap1(i2,3)*mv2kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i12Ckap1kap3Ckap2i23(i1,i2) = Ckap1kap3Ckap2(i2,3)*mv2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i13Ckap1kap3Ckap3i23(i1,i2) = Ckap1kap3Ckap3(i2,3)*mv2kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap2mv2lami2(i1,i2) = Ckap2mv2lam(i2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i11Ckap2kap1Ckap1i21(i1,i2) = Ckap2kap1Ckap1(i2,1)*mv2kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i12Ckap2kap1Ckap2i21(i1,i2) = Ckap2kap1Ckap2(i2,1)*mv2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i13Ckap2kap1Ckap3i21(i1,i2) = Ckap2kap1Ckap3(i2,1)*mv2kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i11Ckap2kap2Ckap1i22(i1,i2) = Ckap2kap2Ckap1(i2,2)*mv2kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i12Ckap2kap2Ckap2i22(i1,i2) = Ckap2kap2Ckap2(i2,2)*mv2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i13Ckap2kap2Ckap3i22(i1,i2) = Ckap2kap2Ckap3(i2,2)*mv2kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i11Ckap2kap3Ckap1i23(i1,i2) = Ckap2kap3Ckap1(i2,3)*mv2kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i12Ckap2kap3Ckap2i23(i1,i2) = Ckap2kap3Ckap2(i2,3)*mv2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i13Ckap2kap3Ckap3i23(i1,i2) = Ckap2kap3Ckap3(i2,3)*mv2kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap3mv2lami2(i1,i2) = Ckap3mv2lam(i2)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i11Ckap3kap1Ckap1i21(i1,i2) = Ckap3kap1Ckap1(i2,1)*mv2kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i12Ckap3kap1Ckap2i21(i1,i2) = Ckap3kap1Ckap2(i2,1)*mv2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i13Ckap3kap1Ckap3i21(i1,i2) = Ckap3kap1Ckap3(i2,1)*mv2kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i11Ckap3kap2Ckap1i22(i1,i2) = Ckap3kap2Ckap1(i2,2)*mv2kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i12Ckap3kap2Ckap2i22(i1,i2) = Ckap3kap2Ckap2(i2,2)*mv2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i13Ckap3kap2Ckap3i22(i1,i2) = Ckap3kap2Ckap3(i2,2)*mv2kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i11Ckap3kap3Ckap1i23(i1,i2) = Ckap3kap3Ckap1(i2,3)*mv2kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i12Ckap3kap3Ckap2i23(i1,i2) = Ckap3kap3Ckap2(i2,3)*mv2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i13Ckap3kap3Ckap3i23(i1,i2) = Ckap3kap3Ckap3(i2,3)*mv2kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvlami1YvClami2(i1,i2) = TeCYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap1kap1i11TpYvCYvi21(i1,i2) = TpYvCYv(i2,1)*TvCkap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi11TvCkap1kap1i21(i1,i2) = TvCkap1kap1(i2,1)*YeCYv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap1kap2i11TpYvCYvi22(i1,i2) = TpYvCYv(i2,2)*TvCkap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi12TvCkap1kap2i21(i1,i2) = TvCkap1kap2(i2,1)*YeCYv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap1kap3i11TpYvCYvi23(i1,i2) = TpYvCYv(i2,3)*TvCkap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi13TvCkap1kap3i21(i1,i2) = TvCkap1kap3(i2,1)*YeCYv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap2kap1i12TpYvCYvi21(i1,i2) = TpYvCYv(i2,1)*TvCkap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi11TvCkap2kap1i22(i1,i2) = TvCkap2kap1(i2,2)*YeCYv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap2kap2i12TpYvCYvi22(i1,i2) = TpYvCYv(i2,2)*TvCkap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi12TvCkap2kap2i22(i1,i2) = TvCkap2kap2(i2,2)*YeCYv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap2kap3i12TpYvCYvi23(i1,i2) = TpYvCYv(i2,3)*TvCkap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi13TvCkap2kap3i22(i1,i2) = TvCkap2kap3(i2,2)*YeCYv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap3kap1i13TpYvCYvi21(i1,i2) = TpYvCYv(i2,1)*TvCkap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi11TvCkap3kap1i23(i1,i2) = TvCkap3kap1(i2,3)*YeCYv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap3kap2i13TpYvCYvi22(i1,i2) = TpYvCYv(i2,2)*TvCkap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi12TvCkap3kap2i23(i1,i2) = TvCkap3kap2(i2,3)*YeCYv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap3kap3i13TpYvCYvi23(i1,i2) = TpYvCYv(i2,3)*TvCkap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi13TvCkap3kap3i23(i1,i2) = TvCkap3kap3(i2,3)*YeCYv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvlami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvlami1adjYvmlHd2i2(i1,i2) = adjYvmlHd2(i2)*TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvlami1Cmv2Clami2(i1,i2) = Cmv2Clam(i2)*TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1TpYvCYvlami2(i1,i2) = TpYvCYvlam(i2)*YvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvClami1TpYvCYvlami2(i1,i2) = TpYvCYvlam(i2)*TvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTlami1CTlami2(i1,i2) = Conjg(Tlam(i2))*TpYvCYvTlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1TpYvCYvTlami2(i1,i2) = TpYvCYvTlam(i2)*YvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap11i1Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*TpYvCYvTpkap1(1,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11TpYvCYvTpkap11i2(i1,i2) = TpYvCYvTpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap12i1Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*TpYvCYvTpkap1(2,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11TpYvCYvTpkap12i2(i1,i2) = TpYvCYvTpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap13i1Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*TpYvCYvTpkap1(3,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11TpYvCYvTpkap13i2(i1,i2) = TpYvCYvTpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap21i1Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*TpYvCYvTpkap2(1,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12TpYvCYvTpkap21i2(i1,i2) = TpYvCYvTpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap22i1Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*TpYvCYvTpkap2(2,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12TpYvCYvTpkap22i2(i1,i2) = TpYvCYvTpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap23i1Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*TpYvCYvTpkap2(3,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12TpYvCYvTpkap23i2(i1,i2) = TpYvCYvTpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap31i1Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*TpYvCYvTpkap3(1,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13TpYvCYvTpkap31i2(i1,i2) = TpYvCYvTpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap32i1Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*TpYvCYvTpkap3(2,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13TpYvCYvTpkap32i2(i1,i2) = TpYvCYvTpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap33i1Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*TpYvCYvTpkap3(3,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13TpYvCYvTpkap33i2(i1,i2) = TpYvCYvTpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCTvTlami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvCTvTlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpTvCYvlami1CTlami2(i1,i2) = Conjg(Tlam(i2))*TpTvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1TpTvCYvlami2(i1,i2) = TpTvCYvlam(i2)*YvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpTvCTvlami1Clami2(i1,i2) = Conjg(lam(i2))*TpTvCTvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap1Ckap1TpTv1i1lami2(i1,i2) = Tpkap1Ckap1TpTv(1,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap1Ckap2TpTv2i1lami2(i1,i2) = Tpkap1Ckap2TpTv(2,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap1Ckap3TpTv3i1lami2(i1,i2) = Tpkap1Ckap3TpTv(3,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap2Ckap1TpTv1i1lami2(i1,i2) = Tpkap2Ckap1TpTv(1,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap2Ckap2TpTv2i1lami2(i1,i2) = Tpkap2Ckap2TpTv(2,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap2Ckap3TpTv3i1lami2(i1,i2) = Tpkap2Ckap3TpTv(3,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap3Ckap1TpTv1i1lami2(i1,i2) = Tpkap3Ckap1TpTv(1,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap3Ckap2TpTv2i1lami2(i1,i2) = Tpkap3Ckap2TpTv(2,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap3Ckap3TpTv3i1lami2(i1,i2) = Tpkap3Ckap3TpTv(3,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1adjYvTvi21(i1,i2) = kap1adjYvTv(i2,1)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1adjYvTvi22(i1,i2) = kap1adjYvTv(i2,2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1adjYvTvi23(i1,i2) = kap1adjYvTv(i2,3)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1adjkap1Tpkap11i2(i1,i2) = kap1adjkap1Tpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1adjkap1Tpkap12i2(i1,i2) = kap1adjkap1Tpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1adjkap1Tpkap13i2(i1,i2) = kap1adjkap1Tpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap1adjkap1Tpkap21i2(i1,i2) = kap1adjkap1Tpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap1adjkap1Tpkap22i2(i1,i2) = kap1adjkap1Tpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap1adjkap1Tpkap23i2(i1,i2) = kap1adjkap1Tpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap1adjkap1Tpkap31i2(i1,i2) = kap1adjkap1Tpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap1adjkap1Tpkap32i2(i1,i2) = kap1adjkap1Tpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap1adjkap1Tpkap33i2(i1,i2) = kap1adjkap1Tpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Clami1Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Clami1Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Clami1Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap1i11Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap1i11Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap2i11Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap2i11Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap3i11Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap3i11Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Ckap1Tk1i21(i1,i2) = kap1Ckap1Tk1(i2,1)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Ckap1Tk2i21(i1,i2) = kap1Ckap1Tk2(i2,1)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Ckap1Tk3i21(i1,i2) = kap1Ckap1Tk3(i2,1)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap1i12Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap1i12Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap2i12Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap2i12Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap3i12Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap3i12Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Ckap2Tk1i22(i1,i2) = kap1Ckap2Tk1(i2,2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Ckap2Tk2i22(i1,i2) = kap1Ckap2Tk2(i2,2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Ckap2Tk3i22(i1,i2) = kap1Ckap2Tk3(i2,2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap1i13Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap1i13Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap2i13Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap2i13Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap3i13Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap3i13Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Ckap3Tk1i23(i1,i2) = kap1Ckap3Tk1(i2,3)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Ckap3Tk2i23(i1,i2) = kap1Ckap3Tk2(i2,3)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Ckap3Tk3i23(i1,i2) = kap1Ckap3Tk3(i2,3)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2adjYvTvi21(i1,i2) = kap2adjYvTv(i2,1)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2adjYvTvi22(i1,i2) = kap2adjYvTv(i2,2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2adjYvTvi23(i1,i2) = kap2adjYvTv(i2,3)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap2adjkap2Tpkap11i2(i1,i2) = kap2adjkap2Tpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap2adjkap2Tpkap12i2(i1,i2) = kap2adjkap2Tpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap2adjkap2Tpkap13i2(i1,i2) = kap2adjkap2Tpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2adjkap2Tpkap21i2(i1,i2) = kap2adjkap2Tpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2adjkap2Tpkap22i2(i1,i2) = kap2adjkap2Tpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2adjkap2Tpkap23i2(i1,i2) = kap2adjkap2Tpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap2adjkap2Tpkap31i2(i1,i2) = kap2adjkap2Tpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap2adjkap2Tpkap32i2(i1,i2) = kap2adjkap2Tpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap2adjkap2Tpkap33i2(i1,i2) = kap2adjkap2Tpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Clami1Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Clami1Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Clami1Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap1i11Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap1i11Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap2i11Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap2i11Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap3i11Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap3i11Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Ckap1Tk1i21(i1,i2) = kap2Ckap1Tk1(i2,1)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Ckap1Tk2i21(i1,i2) = kap2Ckap1Tk2(i2,1)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Ckap1Tk3i21(i1,i2) = kap2Ckap1Tk3(i2,1)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap1i12Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap1i12Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap2i12Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap2i12Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap3i12Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap3i12Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Ckap2Tk1i22(i1,i2) = kap2Ckap2Tk1(i2,2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Ckap2Tk2i22(i1,i2) = kap2Ckap2Tk2(i2,2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Ckap2Tk3i22(i1,i2) = kap2Ckap2Tk3(i2,2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap1i13Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap1i13Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap2i13Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap2i13Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap3i13Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap3i13Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Ckap3Tk1i23(i1,i2) = kap2Ckap3Tk1(i2,3)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Ckap3Tk2i23(i1,i2) = kap2Ckap3Tk2(i2,3)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Ckap3Tk3i23(i1,i2) = kap2Ckap3Tk3(i2,3)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3adjYvTvi21(i1,i2) = kap3adjYvTv(i2,1)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3adjYvTvi22(i1,i2) = kap3adjYvTv(i2,2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3adjYvTvi23(i1,i2) = kap3adjYvTv(i2,3)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap3adjkap3Tpkap11i2(i1,i2) = kap3adjkap3Tpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap3adjkap3Tpkap12i2(i1,i2) = kap3adjkap3Tpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap3adjkap3Tpkap13i2(i1,i2) = kap3adjkap3Tpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap3adjkap3Tpkap21i2(i1,i2) = kap3adjkap3Tpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap3adjkap3Tpkap22i2(i1,i2) = kap3adjkap3Tpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap3adjkap3Tpkap23i2(i1,i2) = kap3adjkap3Tpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3adjkap3Tpkap31i2(i1,i2) = kap3adjkap3Tpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3adjkap3Tpkap32i2(i1,i2) = kap3adjkap3Tpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3adjkap3Tpkap33i2(i1,i2) = kap3adjkap3Tpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Clami1Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Clami1Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Clami1Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap1i11Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap1i11Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap2i11Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap2i11Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap3i11Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap3i11Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Ckap1Tk1i21(i1,i2) = kap3Ckap1Tk1(i2,1)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Ckap1Tk2i21(i1,i2) = kap3Ckap1Tk2(i2,1)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Ckap1Tk3i21(i1,i2) = kap3Ckap1Tk3(i2,1)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap1i12Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap1i12Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap2i12Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap2i12Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap3i12Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap3i12Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Ckap2Tk1i22(i1,i2) = kap3Ckap2Tk1(i2,2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Ckap2Tk2i22(i1,i2) = kap3Ckap2Tk2(i2,2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Ckap2Tk3i22(i1,i2) = kap3Ckap2Tk3(i2,2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap1i13Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap1i13Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap2i13Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap2i13Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap3i13Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap3i13Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Ckap3Tk1i23(i1,i2) = kap3Ckap3Tk1(i2,3)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Ckap3Tk2i23(i1,i2) = kap3Ckap3Tk2(i2,3)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Ckap3Tk3i23(i1,i2) = kap3Ckap3Tk3(i2,3)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11Tk1Ckap1kap1i21(i1,i2) = Tk1Ckap1kap1(i2,1)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12Tk1Ckap1kap2i21(i1,i2) = Tk1Ckap1kap2(i2,1)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13Tk1Ckap1kap3i21(i1,i2) = Tk1Ckap1kap3(i2,1)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11Tk1Ckap2kap1i22(i1,i2) = Tk1Ckap2kap1(i2,2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12Tk1Ckap2kap2i22(i1,i2) = Tk1Ckap2kap2(i2,2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13Tk1Ckap2kap3i22(i1,i2) = Tk1Ckap2kap3(i2,2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11Tk1Ckap3kap1i23(i1,i2) = Tk1Ckap3kap1(i2,3)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12Tk1Ckap3kap2i23(i1,i2) = Tk1Ckap3kap2(i2,3)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13Tk1Ckap3kap3i23(i1,i2) = Tk1Ckap3kap3(i2,3)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11Tk2Ckap1kap1i21(i1,i2) = Tk2Ckap1kap1(i2,1)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12Tk2Ckap1kap2i21(i1,i2) = Tk2Ckap1kap2(i2,1)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13Tk2Ckap1kap3i21(i1,i2) = Tk2Ckap1kap3(i2,1)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11Tk2Ckap2kap1i22(i1,i2) = Tk2Ckap2kap1(i2,2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12Tk2Ckap2kap2i22(i1,i2) = Tk2Ckap2kap2(i2,2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13Tk2Ckap2kap3i22(i1,i2) = Tk2Ckap2kap3(i2,2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11Tk2Ckap3kap1i23(i1,i2) = Tk2Ckap3kap1(i2,3)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12Tk2Ckap3kap2i23(i1,i2) = Tk2Ckap3kap2(i2,3)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13Tk2Ckap3kap3i23(i1,i2) = Tk2Ckap3kap3(i2,3)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11Tk3Ckap1kap1i21(i1,i2) = Tk3Ckap1kap1(i2,1)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12Tk3Ckap1kap2i21(i1,i2) = Tk3Ckap1kap2(i2,1)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13Tk3Ckap1kap3i21(i1,i2) = Tk3Ckap1kap3(i2,1)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11Tk3Ckap2kap1i22(i1,i2) = Tk3Ckap2kap1(i2,2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12Tk3Ckap2kap2i22(i1,i2) = Tk3Ckap2kap2(i2,2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13Tk3Ckap2kap3i22(i1,i2) = Tk3Ckap2kap3(i2,2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11Tk3Ckap3kap1i23(i1,i2) = Tk3Ckap3kap1(i2,3)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12Tk3Ckap3kap2i23(i1,i2) = Tk3Ckap3kap2(i2,3)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13Tk3Ckap3kap3i23(i1,i2) = Tk3Ckap3kap3(i2,3)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2TpYvCYvlami1Clami2(i1,i2) = Conjg(lam(i2))*mv2TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjYvYvClami1lami2(i1,i2) = YvadjYvYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjYvTvClami1lami2(i1,i2) = YvadjYvTvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap1mv2kap1i21(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap1mv2kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap1mv2kap2i21(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap1mv2kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap1mv2kap3i21(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap1mv2kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap2mv2kap1i22(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap2mv2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap2mv2kap2i22(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap2mv2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap2mv2kap3i22(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap2mv2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap3mv2kap1i23(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap3mv2kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap3mv2kap2i23(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap3mv2kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap3mv2kap3i23(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap3mv2kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11YvCmv2Ckap1kap1i21(i1,i2) = Conjg(Yv(i1,1))*YvCmv2Ckap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12YvCmv2Ckap1kap2i21(i1,i2) = Conjg(Yv(i1,2))*YvCmv2Ckap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13YvCmv2Ckap1kap3i21(i1,i2) = Conjg(Yv(i1,3))*YvCmv2Ckap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11YvCmv2Ckap2kap1i22(i1,i2) = Conjg(Yv(i1,1))*YvCmv2Ckap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12YvCmv2Ckap2kap2i22(i1,i2) = Conjg(Yv(i1,2))*YvCmv2Ckap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13YvCmv2Ckap2kap3i22(i1,i2) = Conjg(Yv(i1,3))*YvCmv2Ckap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11YvCmv2Ckap3kap1i23(i1,i2) = Conjg(Yv(i1,1))*YvCmv2Ckap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12YvCmv2Ckap3kap2i23(i1,i2) = Conjg(Yv(i1,2))*YvCmv2Ckap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13YvCmv2Ckap3kap3i23(i1,i2) = Conjg(Yv(i1,3))*YvCmv2Ckap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1Tpkap1Clami1lami2(i1,i2) = YvCkap1Tpkap1Clam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2Tpkap2Clami1lami2(i1,i2) = YvCkap2Tpkap2Clam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3Tpkap3Clami1lami2(i1,i2) = YvCkap3Tpkap3Clam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyadjYeYeCYvlami1mlHd2i2(i1,i2) = adjYeYeCYvlam(i1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvYvadjYvmlHd2i2(i1,i2) = adjYvYvadjYvmlHd2(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvYvCmv2Clami2(i1,i2) = adjYvYvCmv2Clam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvCml2YvClami2(i1,i2) = adjYvCml2YvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1Cmv2adjYvYvClami2(i1,i2) = Cmv2adjYvYvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpYvCYvlami1mlHd2i2(i1,i2) = CYvTpYvCYvlam(i1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvadjYvYvClami1lami2(i1,i2) = TvadjYvYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYeCYeYvClami1lami2(i1,i2) = TpYeCYeYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYeCYeTvClami1lami2(i1,i2) = TpYeCYeTvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvml2CYvlami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvml2CYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvmv2lami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvCYvmv2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpTeCYeYvClami1lami2(i1,i2) = TpTeCYeYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap1kap1i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap1kap2i11Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap1kap3i11Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap2kap1i12Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap2kap2i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap2kap3i12Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap3kap1i13Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap3kap2i13Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap3kap3i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap1kap1i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap1kap2i11Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap1kap3i11Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap2kap1i12Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap2kap2i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap2kap3i12Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap3kap1i13Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap3kap2i13Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap3kap3i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap1kap1i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap1kap2i11Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap1kap3i11Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap2kap1i12Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap2kap2i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap2kap3i12Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap3kap1i13Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap3kap2i13Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap3kap3i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1i22(i1,i2) = kap1Ckap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1i23(i1,i2) = kap1Ckap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1i21(i1,i2) = kap2Ckap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1i23(i1,i2) = kap2Ckap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1i21(i1,i2) = kap3Ckap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1i22(i1,i2) = kap3Ckap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1lami2(i1,i2) = kap1Ckap1lam(i2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1lami2(i1,i2) = kap2Ckap1lam(i2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1lami2(i1,i2) = kap3Ckap1lam(i2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1adjYvYvCkap1i21(i1,i2) = kap1adjYvYvCkap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2adjYvYvCkap1i22(i1,i2) = kap2adjYvYvCkap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3adjYvYvCkap1i23(i1,i2) = kap3adjYvYvCkap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2i22(i1,i2) = kap1Ckap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2i23(i1,i2) = kap1Ckap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2i21(i1,i2) = kap2Ckap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2i23(i1,i2) = kap2Ckap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2i21(i1,i2) = kap3Ckap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2i22(i1,i2) = kap3Ckap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2lami2(i1,i2) = kap1Ckap2lam(i2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2lami2(i1,i2) = kap2Ckap2lam(i2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2lami2(i1,i2) = kap3Ckap2lam(i2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1adjYvYvCkap2i21(i1,i2) = kap1adjYvYvCkap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2adjYvYvCkap2i22(i1,i2) = kap2adjYvYvCkap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3adjYvYvCkap2i23(i1,i2) = kap3adjYvYvCkap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3i22(i1,i2) = kap1Ckap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3i23(i1,i2) = kap1Ckap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3i21(i1,i2) = kap2Ckap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3i23(i1,i2) = kap2Ckap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3i21(i1,i2) = kap3Ckap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3i22(i1,i2) = kap3Ckap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3lami2(i1,i2) = kap1Ckap3lam(i2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3lami2(i1,i2) = kap2Ckap3lam(i2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3lami2(i1,i2) = kap3Ckap3lam(i2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1adjYvYvCkap3i21(i1,i2) = kap1adjYvYvCkap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2adjYvYvCkap3i22(i1,i2) = kap2adjYvYvCkap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3adjYvYvCkap3i23(i1,i2) = kap3adjYvYvCkap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Tvi21(i1,i2) = Conjg(Yv(i1,1))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Tvi21(i1,i2) = Conjg(Yv(i1,2))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Tvi21(i1,i2) = Conjg(Yv(i1,3))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11Tvi21(i1,i2) = Conjg(Tv(i1,1))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12Tvi21(i1,i2) = Conjg(Tv(i1,2))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13Tvi21(i1,i2) = Conjg(Tv(i1,3))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Tvi22(i1,i2) = Conjg(Yv(i1,1))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Tvi22(i1,i2) = Conjg(Yv(i1,2))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Tvi22(i1,i2) = Conjg(Yv(i1,3))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11Tvi22(i1,i2) = Conjg(Tv(i1,1))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12Tvi22(i1,i2) = Conjg(Tv(i1,2))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13Tvi22(i1,i2) = Conjg(Tv(i1,3))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Tvi23(i1,i2) = Conjg(Yv(i1,1))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Tvi23(i1,i2) = Conjg(Yv(i1,2))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Tvi23(i1,i2) = Conjg(Yv(i1,3))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11Tvi23(i1,i2) = Conjg(Tv(i1,1))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12Tvi23(i1,i2) = Conjg(Tv(i1,2))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13Tvi23(i1,i2) = Conjg(Tv(i1,3))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1Clami2(i1,i2) = Conjg(lam(i2))*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1adjYvYvCTlami2(i1,i2) = adjYvYvCTlam(i2)*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1adjTvYvClami2(i1,i2) = adjTvYvClam(i2)*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap1i11Tlami2(i1,i2) = Yvadjkap1kap1(i1,1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap2i11Tlami2(i1,i2) = Yvadjkap1kap2(i1,1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap3i11Tlami2(i1,i2) = Yvadjkap1kap3(i1,1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap1i12Tlami2(i1,i2) = Yvadjkap2kap1(i1,2)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap2i12Tlami2(i1,i2) = Yvadjkap2kap2(i1,2)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap3i12Tlami2(i1,i2) = Yvadjkap2kap3(i1,2)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap1i13Tlami2(i1,i2) = Yvadjkap3kap1(i1,3)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap2i13Tlami2(i1,i2) = Yvadjkap3kap2(i1,3)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap3i13Tlami2(i1,i2) = Yvadjkap3kap3(i1,3)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjYvYvClami1Tlami2(i1,i2) = YvadjYvYvClam(i1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYeCYeYvClami1Tlami2(i1,i2) = TpYeCYeYvClam(i1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i11Ckap1TpYvCYvi21(i1,i2) = Ckap1TpYvCYv(i2,1)*kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i11adjYvCml2Yvadjkap11i2(i1,i2) = adjYvCml2Yvadjkap1(1,i2)*kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i12Ckap2TpYvCYvi21(i1,i2) = Ckap2TpYvCYv(i2,1)*kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i12adjYvCml2Yvadjkap21i2(i1,i2) = adjYvCml2Yvadjkap2(1,i2)*kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i13Ckap3TpYvCYvi21(i1,i2) = Ckap3TpYvCYv(i2,1)*kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i13adjYvCml2Yvadjkap31i2(i1,i2) = adjYvCml2Yvadjkap3(1,i2)*kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i11Ckap1TpYvCYvi22(i1,i2) = Ckap1TpYvCYv(i2,2)*kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i11adjYvCml2Yvadjkap12i2(i1,i2) = adjYvCml2Yvadjkap1(2,i2)*kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i12Ckap2TpYvCYvi22(i1,i2) = Ckap2TpYvCYv(i2,2)*kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i12adjYvCml2Yvadjkap22i2(i1,i2) = adjYvCml2Yvadjkap2(2,i2)*kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i13Ckap3TpYvCYvi22(i1,i2) = Ckap3TpYvCYv(i2,2)*kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i13adjYvCml2Yvadjkap32i2(i1,i2) = adjYvCml2Yvadjkap3(2,i2)*kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i11Ckap1TpYvCYvi23(i1,i2) = Ckap1TpYvCYv(i2,3)*kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i11adjYvCml2Yvadjkap13i2(i1,i2) = adjYvCml2Yvadjkap1(3,i2)*kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i12Ckap2TpYvCYvi23(i1,i2) = Ckap2TpYvCYv(i2,3)*kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i12adjYvCml2Yvadjkap23i2(i1,i2) = adjYvCml2Yvadjkap2(3,i2)*kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i13Ckap3TpYvCYvi23(i1,i2) = Ckap3TpYvCYv(i2,3)*kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i13adjYvCml2Yvadjkap33i2(i1,i2) = adjYvCml2Yvadjkap3(3,i2)*kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11CTk1i21(i1,i2) = Conjg(Tk1(i2,1))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11CTk1i22(i1,i2) = Conjg(Tk1(i2,2))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11CTk1i23(i1,i2) = Conjg(Tk1(i2,3))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1lami2(i1,i2) = Ckap1lam(i2)*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1TpYvCTvi21(i1,i2) = Ckap1TpYvCTv(i2,1)*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12CTk1i21(i1,i2) = Conjg(Tk1(i2,1))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12CTk1i22(i1,i2) = Conjg(Tk1(i2,2))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12CTk1i23(i1,i2) = Conjg(Tk1(i2,3))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2lami2(i1,i2) = Ckap2lam(i2)*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2TpYvCTvi21(i1,i2) = Ckap2TpYvCTv(i2,1)*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13CTk1i21(i1,i2) = Conjg(Tk1(i2,1))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13CTk1i22(i1,i2) = Conjg(Tk1(i2,2))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13CTk1i23(i1,i2) = Conjg(Tk1(i2,3))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3lami2(i1,i2) = Ckap3lam(i2)*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3TpYvCTvi21(i1,i2) = Ckap3TpYvCTv(i2,1)*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11CTk2i21(i1,i2) = Conjg(Tk2(i2,1))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11CTk2i22(i1,i2) = Conjg(Tk2(i2,2))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11CTk2i23(i1,i2) = Conjg(Tk2(i2,3))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1lami2(i1,i2) = Ckap1lam(i2)*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1TpYvCTvi22(i1,i2) = Ckap1TpYvCTv(i2,2)*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12CTk2i21(i1,i2) = Conjg(Tk2(i2,1))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12CTk2i22(i1,i2) = Conjg(Tk2(i2,2))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12CTk2i23(i1,i2) = Conjg(Tk2(i2,3))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2lami2(i1,i2) = Ckap2lam(i2)*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2TpYvCTvi22(i1,i2) = Ckap2TpYvCTv(i2,2)*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13CTk2i21(i1,i2) = Conjg(Tk2(i2,1))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13CTk2i22(i1,i2) = Conjg(Tk2(i2,2))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13CTk2i23(i1,i2) = Conjg(Tk2(i2,3))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3lami2(i1,i2) = Ckap3lam(i2)*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3TpYvCTvi22(i1,i2) = Ckap3TpYvCTv(i2,2)*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11CTk3i21(i1,i2) = Conjg(Tk3(i2,1))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11CTk3i22(i1,i2) = Conjg(Tk3(i2,2))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11CTk3i23(i1,i2) = Conjg(Tk3(i2,3))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1lami2(i1,i2) = Ckap1lam(i2)*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1TpYvCTvi23(i1,i2) = Ckap1TpYvCTv(i2,3)*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12CTk3i21(i1,i2) = Conjg(Tk3(i2,1))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12CTk3i22(i1,i2) = Conjg(Tk3(i2,2))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12CTk3i23(i1,i2) = Conjg(Tk3(i2,3))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2lami2(i1,i2) = Ckap2lam(i2)*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2TpYvCTvi23(i1,i2) = Ckap2TpYvCTv(i2,3)*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13CTk3i21(i1,i2) = Conjg(Tk3(i2,1))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13CTk3i22(i1,i2) = Conjg(Tk3(i2,2))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13CTk3i23(i1,i2) = Conjg(Tk3(i2,3))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3lami2(i1,i2) = Ckap3lam(i2)*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3TpYvCTvi23(i1,i2) = Ckap3TpYvCTv(i2,3)*Tk3(i1,3) 
  End Do 
End Do 
End If 
 
 
Tr1(1) = g1*sqrt3ov5*(-1._dp*(mHd2) + mHu2 + Trmd2 + Trme2 - Trml2 + Trmq2 -          & 
&  2._dp*(Trmu2))

If (TwoLoopRGE) Then 
Tr2U1(1, 1) = (g1p2*(3._dp*(mHd2) + 3._dp*(mHu2) + 2._dp*(Trmd2) + 6._dp*(Trme2)      & 
&  + 3._dp*(Trml2) + Trmq2 + 8._dp*(Trmu2)))/10._dp

Tr3(1) = (g1*ooSqrt15*(-9*g1p2*mHd2 - 45*g2p2*mHd2 + 9*g1p2*mHu2 + 45*g2p2*mHu2 -     & 
&  30._dp*(SPClamxxTpYvmlHd2) + 30*(mHd2 - mHu2)*SPlamxxClam + 4*g1p2*Trmd2 +            & 
&  80*g3p2*Trmd2 + 36*g1p2*Trme2 - 9*g1p2*Trml2 - 45*g2p2*Trml2 + 30._dp*(Trml2YvadjYv)  & 
&  + g1p2*Trmq2 + 45*g2p2*Trmq2 + 80*g3p2*Trmq2 - 32*g1p2*Trmu2 - 160*g3p2*Trmu2 +       & 
&  90*mHd2*TrYdadjYd - 60._dp*(TrYdadjYdCmd2) - 30._dp*(TrYdCmq2adjYd) + 30*mHd2*TrYeadjYe -& 
&  60._dp*(TrYeadjYeCme2) + 30._dp*(TrYeCml2adjYe) - 90*mHu2*TrYuadjYu + 120._dp*(TrYuadjYuCmu2)& 
&  - 30._dp*(TrYuCmq2adjYu) - 30*mHu2*TrYvadjYv))/20._dp

Tr2(2) = (mHd2 + mHu2 + Trml2 + 3._dp*(Trmq2))/2._dp

Tr2(3) = (Trmd2 + 2._dp*(Trmq2) + Trmu2)/2._dp

End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = 33._dp*(g1p3)/5._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (199._dp*(g1p5) + 5*g1p3*(27._dp*(g2p2) + 88._dp*(g3p2) - 6._dp*(SPlamxxClam) -       & 
&  14._dp*(TrYdadjYd) - 18._dp*(TrYeadjYe) - 26._dp*(TrYuadjYu) - 6._dp*(TrYvadjYv)))/25._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = g2p3

 
 
If (TwoLoopRGE) Then 
betag22 = (9*g1p2*g2p3)/5._dp + 25._dp*(g2p5) + 2*g2p3*(12._dp*(g3p2) - SPlamxxClam -           & 
&  3._dp*(TrYdadjYd) - TrYeadjYe - 3._dp*(TrYuadjYu) - TrYvadjYv)

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = -3._dp*(g3p3)

 
 
If (TwoLoopRGE) Then 
betag32 = (11*g1p2*g3p3)/5._dp + 9*g2p2*g3p3 + 14._dp*(g3p5) - 4*g3p3*TrYdadjYd -               & 
&  4*g3p3*TrYuadjYu

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = (-7._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)/3._dp +               & 
&  SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd) + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = ((287._dp*(g1p4) + 2*g1p2*(45._dp*(g2p2) + 40._dp*(g3p2) - 18._dp*(TrYdadjYd) +       & 
&  54._dp*(TrYeadjYe)) + 5*(135._dp*(g2p4) + 2*(72*g2p2*g3p2 - 16._dp*(g3p4) -           & 
&  9*(2._dp*(SPlamxxadjYvYvClam) + 3._dp*(SPlamxxClamp2) + 2._dp*(SPTpkap1Clamxxadjkap1lam) +& 
&  2._dp*(SPTpkap2Clamxxadjkap2lam) + 2._dp*(SPTpkap3Clamxxadjkap3lam) - 16*g3p2*TrYdadjYd +& 
&  9._dp*(TrYdadjYdYdadjYd) + 3._dp*(TrYdadjYuYuadjYd) + 3._dp*(TrYeadjYeYeadjYe) +      & 
&  3*SPlamxxClam*TrYuadjYu + SPlamxxClam*TrYvadjYv + TrYvadjYvTpYeCYe))))*Yd)/90._dp +   & 
&  (4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 3*(SPlamxxClam + 3._dp*(TrYdadjYd) +             & 
&  TrYeadjYe))*YdadjYdYd - 4._dp*(YdadjYdYdadjYdYd) + (4*g1p2*YdadjYuYu)/5._dp -         & 
&  SPlamxxClam*YdadjYuYu - 3*TrYuadjYu*YdadjYuYu - TrYvadjYv*YdadjYuYu - 2._dp*(YdadjYuYuadjYdYd) -& 
&  2._dp*(YdadjYuYuadjYuYu)

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
Call Chop(DYd) 

!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (-9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + SPlamxxClam + 3._dp*(TrYdadjYd)      & 
&  + TrYeadjYe)*Ye + 3._dp*(YeadjYeYe) + YeCYvTpYv

 
 
If (TwoLoopRGE) Then 
betaYe2 = -2._dp*(DyYeCYvlami1YvClami2) + ((135._dp*(g1p4) + 2*g1p2*(9._dp*(g2p2) -             & 
&  2._dp*(TrYdadjYd) + 6._dp*(TrYeadjYe)) + 5*(15._dp*(g2p4) - 2*(2._dp*(SPlamxxadjYvYvClam) +& 
&  3._dp*(SPlamxxClamp2) + 2._dp*(SPTpkap1Clamxxadjkap1lam) + 2._dp*(SPTpkap2Clamxxadjkap2lam) +& 
&  2._dp*(SPTpkap3Clamxxadjkap3lam) - 16*g3p2*TrYdadjYd + 9._dp*(TrYdadjYdYdadjYd) +     & 
&  3._dp*(TrYdadjYuYuadjYd) + 3._dp*(TrYeadjYeYeadjYe) + 3*SPlamxxClam*TrYuadjYu +       & 
&  SPlamxxClam*TrYvadjYv + TrYvadjYvTpYeCYe)))*Ye)/10._dp + 6*g2p2*YeadjYeYe -           & 
&  3*SPlamxxClam*YeadjYeYe - 9*TrYdadjYd*YeadjYeYe - 3*TrYeadjYe*YeadjYeYe -             & 
&  4._dp*(YeadjYeYeadjYeYe) - 2._dp*(YeCYvkap1adjkap1TpYv) - 2._dp*(YeCYvkap2adjkap2TpYv) -& 
&  2._dp*(YeCYvkap3adjkap3TpYv) - SPlamxxClam*YeCYvTpYv - 3*TrYuadjYu*YeCYvTpYv -        & 
&  TrYvadjYv*YeCYvTpYv - 2._dp*(YeCYvTpYvadjYeYe) - 2._dp*(YeCYvTpYvCYvTpYv)

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
Call Chop(DYe) 

!-------------------- 
! lam 
!-------------------- 
 
Do i1 = 1,3
betalam1(i1) = 3*TpYvCYvlam(i1) + 2*kap1Ckap1(i1,1)*lam(1) + 2*kap2Ckap1(i1,2)*lam(1) +              & 
&  2*kap3Ckap1(i1,3)*lam(1) + 2*kap1Ckap2(i1,1)*lam(2) + 2*kap2Ckap2(i1,2)*lam(2) +      & 
&  2*kap3Ckap2(i1,3)*lam(2) + 2*kap1Ckap3(i1,1)*lam(3) + 2*kap2Ckap3(i1,2)*lam(3) +      & 
&  2*kap3Ckap3(i1,3)*lam(3) - (3*g1p2*lam(i1))/5._dp - 3*g2p2*lam(i1) + 4*SPlamxxClam*lam(i1) +& 
&  3*TrYdadjYd*lam(i1) + TrYeadjYe*lam(i1) + 3*TrYuadjYu*lam(i1) + TrYvadjYv*lam(i1)
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
betalam2(i1) = -8*adjkap1lam(1)*kap1adjkap1Tpkap1(1,i1) - 8*adjkap2lam(1)*kap1adjkap1Tpkap1(2,       & 
& i1) - 8*adjkap3lam(1)*kap1adjkap1Tpkap1(3,i1) - 8*adjkap1lam(2)*kap1adjkap1Tpkap2(1,   & 
& i1) - 8*adjkap2lam(2)*kap1adjkap1Tpkap2(2,i1) - 8*adjkap3lam(2)*kap1adjkap1Tpkap2(3,   & 
& i1) - 8*adjkap1lam(3)*kap1adjkap1Tpkap3(1,i1) - 8*adjkap2lam(3)*kap1adjkap1Tpkap3(2,   & 
& i1) - 8*adjkap3lam(3)*kap1adjkap1Tpkap3(3,i1) - 8*adjkap1lam(1)*kap2adjkap2Tpkap1(1,   & 
& i1) - 8*adjkap2lam(1)*kap2adjkap2Tpkap1(2,i1) - 8*adjkap3lam(1)*kap2adjkap2Tpkap1(3,   & 
& i1) - 8*adjkap1lam(2)*kap2adjkap2Tpkap2(1,i1) - 8*adjkap2lam(2)*kap2adjkap2Tpkap2(2,   & 
& i1) - 8*adjkap3lam(2)*kap2adjkap2Tpkap2(3,i1) - 8*adjkap1lam(3)*kap2adjkap2Tpkap3(1,   & 
& i1) - 8*adjkap2lam(3)*kap2adjkap2Tpkap3(2,i1) - 8*adjkap3lam(3)*kap2adjkap2Tpkap3(3,   & 
& i1) - 8*adjkap1lam(1)*kap3adjkap3Tpkap1(1,i1) - 8*adjkap2lam(1)*kap3adjkap3Tpkap1(2,   & 
& i1) - 8*adjkap3lam(1)*kap3adjkap3Tpkap1(3,i1) - 8*adjkap1lam(2)*kap3adjkap3Tpkap2(1,   & 
& i1) - 8*adjkap2lam(2)*kap3adjkap3Tpkap2(2,i1) - 8*adjkap3lam(2)*kap3adjkap3Tpkap2(3,   & 
& i1) - 8*adjkap1lam(3)*kap3adjkap3Tpkap3(1,i1) - 8*adjkap2lam(3)*kap3adjkap3Tpkap3(2,   & 
& i1) - 8*adjkap3lam(3)*kap3adjkap3Tpkap3(3,i1) - TpYvadjYeYeCYvlam(i1) - 2*TpYvCYvkap1adjkap1lam(i1) -& 
&  2*TpYvCYvkap2adjkap2lam(i1) - 2*TpYvCYvkap3adjkap3lam(i1) + (6*g1p2*TpYvCYvlam(i1))/5._dp +& 
&  6*g2p2*TpYvCYvlam(i1) - 7*SPlamxxClam*TpYvCYvlam(i1) - 9*TrYuadjYu*TpYvCYvlam(i1) -   & 
&  3*TrYvadjYv*TpYvCYvlam(i1) - 8*adjkap1lam(1)*TpYvCYvTpkap1(1,i1) - 8*adjkap2lam(1)*TpYvCYvTpkap1(2,& 
& i1) - 8*adjkap3lam(1)*TpYvCYvTpkap1(3,i1) - 8*adjkap1lam(2)*TpYvCYvTpkap2(1,           & 
& i1) - 8*adjkap2lam(2)*TpYvCYvTpkap2(2,i1) - 8*adjkap3lam(2)*TpYvCYvTpkap2(3,           & 
& i1) - 8*adjkap1lam(3)*TpYvCYvTpkap3(1,i1) - 8*adjkap2lam(3)*TpYvCYvTpkap3(2,           & 
& i1) - 8*adjkap3lam(3)*TpYvCYvTpkap3(3,i1) - 4*TpYvCYvTpYvCYvlam(i1) - 8*adjkap1lam(1)*kap1Clam(i1)*lam(1)  
betalam2(i1) =  betalam2(i1)- 8*adjkap1lam(2)*kap2Clam(i1)*lam(1) - 8*adjkap1lam(3)*kap3Clam(i1)*lam(1) -           & 
&  8*adjkap2lam(1)*kap1Clam(i1)*lam(2) - 8*adjkap2lam(2)*kap2Clam(i1)*lam(2) -           & 
&  8*adjkap2lam(3)*kap3Clam(i1)*lam(2) - 8*adjkap3lam(1)*kap1Clam(i1)*lam(3) -           & 
&  8*adjkap3lam(2)*kap2Clam(i1)*lam(3) - 8*adjkap3lam(3)*kap3Clam(i1)*lam(3) +           & 
&  (207*g1p4*lam(i1))/50._dp + (9*g1p2*g2p2*lam(i1))/5._dp + (15*g2p4*lam(i1))/2._dp -   & 
&  10*SPlamxxadjYvYvClam*lam(i1) + (6*g1p2*SPlamxxClam*lam(i1))/5._dp + 6*g2p2*SPlamxxClam*lam(i1) -& 
&  10*SPlamxxClamp2*lam(i1) - 4*SPTpkap1Clamxxadjkap1lam*lam(i1) - 4*SPTpkap2Clamxxadjkap2lam*lam(i1) -& 
&  4*SPTpkap3Clamxxadjkap3lam*lam(i1) - (2*g1p2*TrYdadjYd*lam(i1))/5._dp +               & 
&  16*g3p2*TrYdadjYd*lam(i1) - 9*SPlamxxClam*TrYdadjYd*lam(i1) - 9*TrYdadjYdYdadjYd*lam(i1) -& 
&  6*TrYdadjYuYuadjYd*lam(i1) + (6*g1p2*TrYeadjYe*lam(i1))/5._dp - 3*SPlamxxClam*TrYeadjYe*lam(i1) -& 
&  3*TrYeadjYeYeadjYe*lam(i1) + (4*g1p2*TrYuadjYu*lam(i1))/5._dp + 16*g3p2*TrYuadjYu*lam(i1) -& 
&  9*SPlamxxClam*TrYuadjYu*lam(i1) - 9*TrYuadjYuYuadjYu*lam(i1) - 3*SPlamxxClam*TrYvadjYv*lam(i1) -& 
&  2*TrYvadjYvTpYeCYe*lam(i1) - 3*TrYvadjYvYvadjYv*lam(i1) - 2*TrYvCkap1Tpkap1adjYv*lam(i1) -& 
&  2*TrYvCkap2Tpkap2adjYv*lam(i1) - 2*TrYvCkap3Tpkap3adjYv*lam(i1)
End Do

 
Dlam = oo16pi2*( betalam1 + oo16pi2 * betalam2 ) 

 
Else 
Dlam = oo16pi2* betalam1 
End If 
 
 
Call Chop(Dlam) 

!-------------------- 
! Yv 
!-------------------- 
 
betaYv1  = 3._dp*(DyYvClami1lami2) + TpYeCYeYv - (3*g1p2*Yv)/5._dp - 3*g2p2*Yv +      & 
&  SPlamxxClam*Yv + 3*TrYuadjYu*Yv + TrYvadjYv*Yv + 3._dp*(YvadjYvYv) + 2._dp*(YvCkap1Tpkap1)& 
&  + 2._dp*(YvCkap2Tpkap2) + 2._dp*(YvCkap3Tpkap3)

 
 
If (TwoLoopRGE) Then 
betaYv2 = DyTpYeCYeYvClami1lami2 - 4._dp*(DyYvadjYvYvClami1lami2) - 8._dp*(DyYvCkap1i11kap1adjkap1Tpkap11i2) -& 
&  8._dp*(DyYvCkap1i11kap2adjkap2Tpkap11i2) - 8._dp*(DyYvCkap1i11kap3adjkap3Tpkap11i2) - & 
&  8._dp*(DyYvCkap1i11TpYvCYvTpkap11i2) - 8._dp*(DyYvCkap1i12kap1adjkap1Tpkap21i2) -     & 
&  8._dp*(DyYvCkap1i12kap2adjkap2Tpkap21i2) - 8._dp*(DyYvCkap1i12kap3adjkap3Tpkap21i2) - & 
&  8._dp*(DyYvCkap1i12TpYvCYvTpkap21i2) - 8._dp*(DyYvCkap1i13kap1adjkap1Tpkap31i2) -     & 
&  8._dp*(DyYvCkap1i13kap2adjkap2Tpkap31i2) - 8._dp*(DyYvCkap1i13kap3adjkap3Tpkap31i2) - & 
&  8._dp*(DyYvCkap1i13TpYvCYvTpkap31i2) - 2._dp*(DyYvCkap1Tpkap1Clami1lami2) -           & 
&  8._dp*(DyYvCkap2i11kap1adjkap1Tpkap12i2) - 8._dp*(DyYvCkap2i11kap2adjkap2Tpkap12i2) - & 
&  8._dp*(DyYvCkap2i11kap3adjkap3Tpkap12i2) - 8._dp*(DyYvCkap2i11TpYvCYvTpkap12i2) -     & 
&  8._dp*(DyYvCkap2i12kap1adjkap1Tpkap22i2) - 8._dp*(DyYvCkap2i12kap2adjkap2Tpkap22i2) - & 
&  8._dp*(DyYvCkap2i12kap3adjkap3Tpkap22i2) - 8._dp*(DyYvCkap2i12TpYvCYvTpkap22i2) -     & 
&  8._dp*(DyYvCkap2i13kap1adjkap1Tpkap32i2) - 8._dp*(DyYvCkap2i13kap2adjkap2Tpkap32i2) - & 
&  8._dp*(DyYvCkap2i13kap3adjkap3Tpkap32i2) - 8._dp*(DyYvCkap2i13TpYvCYvTpkap32i2) -     & 
&  2._dp*(DyYvCkap2Tpkap2Clami1lami2) - 8._dp*(DyYvCkap3i11kap1adjkap1Tpkap13i2) -       & 
&  8._dp*(DyYvCkap3i11kap2adjkap2Tpkap13i2) - 8._dp*(DyYvCkap3i11kap3adjkap3Tpkap13i2) - & 
&  8._dp*(DyYvCkap3i11TpYvCYvTpkap13i2) - 8._dp*(DyYvCkap3i12kap1adjkap1Tpkap23i2) -     & 
&  8._dp*(DyYvCkap3i12kap2adjkap2Tpkap23i2) - 8._dp*(DyYvCkap3i12kap3adjkap3Tpkap23i2) - & 
&  8._dp*(DyYvCkap3i12TpYvCYvTpkap23i2) - 8._dp*(DyYvCkap3i13kap1adjkap1Tpkap33i2) -     & 
&  8._dp*(DyYvCkap3i13kap2adjkap2Tpkap33i2) - 8._dp*(DyYvCkap3i13kap3adjkap3Tpkap33i2) - & 
&  8._dp*(DyYvCkap3i13TpYvCYvTpkap33i2) - 2._dp*(DyYvCkap3Tpkap3Clami1lami2) -           & 
&  4._dp*(DyYvClami1TpYvCYvlami2) + (6*DyYvClami1lami2*g1p2)/5._dp + 6*DyYvClami1lami2*g2p2  
betaYv2 =  betaYv2- 7*DyYvClami1lami2*SPlamxxClam - 2._dp*(TpYeCYeTpYeCYeYv) + (6*g1p2*TpYeCYeYv)/5._dp - & 
&  SPlamxxClam*TpYeCYeYv - 6*DyYvClami1lami2*TrYdadjYd - 3*TpYeCYeYv*TrYdadjYd -         & 
&  2*DyYvClami1lami2*TrYeadjYe - TpYeCYeYv*TrYeadjYe - 9*DyYvClami1lami2*TrYuadjYu -     & 
&  3*DyYvClami1lami2*TrYvadjYv + (207*g1p4*Yv)/50._dp + (9*g1p2*g2p2*Yv)/5._dp +         & 
&  (15*g2p4*Yv)/2._dp - 6*SPlamxxadjYvYvClam*Yv - 3*SPlamxxClamp2*Yv - 2*SPTpkap1Clamxxadjkap1lam*Yv -& 
&  2*SPTpkap2Clamxxadjkap2lam*Yv - 2*SPTpkap3Clamxxadjkap3lam*Yv - 3*SPlamxxClam*TrYdadjYd*Yv -& 
&  3*TrYdadjYuYuadjYd*Yv - SPlamxxClam*TrYeadjYe*Yv + (4*g1p2*TrYuadjYu*Yv)/5._dp +      & 
&  16*g3p2*TrYuadjYu*Yv - 9*TrYuadjYuYuadjYu*Yv - TrYvadjYvTpYeCYe*Yv - 3*TrYvadjYvYvadjYv*Yv -& 
&  2*TrYvCkap1Tpkap1adjYv*Yv - 2*TrYvCkap2Tpkap2adjYv*Yv - 2*TrYvCkap3Tpkap3adjYv*Yv -   & 
&  2._dp*(YvadjYvTpYeCYeYv) + (6*g1p2*YvadjYvYv)/5._dp + 6*g2p2*YvadjYvYv -              & 
&  3*SPlamxxClam*YvadjYvYv - 9*TrYuadjYu*YvadjYvYv - 3*TrYvadjYv*YvadjYvYv -             & 
&  4._dp*(YvadjYvYvadjYvYv) - 2._dp*(YvCkap1Tpkap1adjYvYv) - 2._dp*(YvCkap2Tpkap2adjYvYv) -& 
&  2._dp*(YvCkap3Tpkap3adjYvYv) - 8*DyYvCkap1i11kap1Clami2*lam(1) - 8*DyYvCkap1i12kap2Clami2*lam(1) -& 
&  8*DyYvCkap1i13kap3Clami2*lam(1) - 8*DyYvCkap2i11kap1Clami2*lam(2) - 8*DyYvCkap2i12kap2Clami2*lam(2) -& 
&  8*DyYvCkap2i13kap3Clami2*lam(2) - 8*DyYvCkap3i11kap1Clami2*lam(3) - 8*DyYvCkap3i12kap2Clami2*lam(3) -& 
&  8*DyYvCkap3i13kap3Clami2*lam(3)

 
DYv = oo16pi2*( betaYv1 + oo16pi2 * betaYv2 ) 

 
Else 
DYv = oo16pi2* betaYv1 
End If 
 
 
Call Chop(DYv) 

!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = (-13._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)/3._dp +              & 
&  SPlamxxClam + 3._dp*(TrYuadjYu) + TrYvadjYv)*Yu + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = ((2743._dp*(g1p4) + 5*(g1p2*(90._dp*(g2p2) + 272._dp*(g3p2) + 72._dp*(TrYuadjYu)) +   & 
&  5*(135._dp*(g2p4) + 2*(72*g2p2*g3p2 - 16._dp*(g3p4) - 9*(6._dp*(SPlamxxadjYvYvClam) + & 
&  3._dp*(SPlamxxClamp2) + 2._dp*(SPTpkap1Clamxxadjkap1lam) + 2._dp*(SPTpkap2Clamxxadjkap2lam) +& 
&  2._dp*(SPTpkap3Clamxxadjkap3lam) + 3*SPlamxxClam*TrYdadjYd + 3._dp*(TrYdadjYuYuadjYd) +& 
&  SPlamxxClam*TrYeadjYe - 16*g3p2*TrYuadjYu + 9._dp*(TrYuadjYuYuadjYu) + TrYvadjYvTpYeCYe +& 
&  3._dp*(TrYvadjYvYvadjYv) + 2._dp*(TrYvCkap1Tpkap1adjYv) + 2._dp*(TrYvCkap2Tpkap2adjYv) +& 
&  2._dp*(TrYvCkap3Tpkap3adjYv))))))*Yu)/450._dp + ((2._dp*(g1p2) - 5*(SPlamxxClam +     & 
&  3._dp*(TrYdadjYd) + TrYeadjYe))*YuadjYdYd)/5._dp - 2._dp*(YuadjYdYdadjYdYd) -         & 
&  2._dp*(YuadjYdYdadjYuYu) + (2*g1p2*YuadjYuYu)/5._dp + 6*g2p2*YuadjYuYu -              & 
&  3*SPlamxxClam*YuadjYuYu - 9*TrYuadjYu*YuadjYuYu - 3*TrYvadjYv*YuadjYuYu -             & 
&  4._dp*(YuadjYuYuadjYuYu)

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
Call Chop(DYu) 

!-------------------- 
! kap 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betakap1(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap1(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap1(i1,i2,3) = 0
End Do
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betakap2(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap2(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap2(i1,i2,3) = 0
End Do
End Do

 
Dkap = oo16pi2*( betakap1 + oo16pi2 * betakap2 ) 

 
Else 
Dkap = oo16pi2* betakap1 
End If 
 
 
!-------------------- 
! Td 
!-------------------- 
 
betaTd1  = 5._dp*(TdadjYdYd) + TdadjYuYu + (14*g1p2*M1*Yd)/15._dp + (32*g3p2*M3*Yd)   & 
& /3._dp + 6*g2p2*M2*Yd + 2*SPClamxxTlam*Yd + 6*TradjYdTd*Yd + 2*TradjYeTe*Yd +          & 
&  4._dp*(YdadjYdTd) + 2._dp*(YdadjYuTu) - (7*g1p2*Td)/15._dp - 3*g2p2*Td -              & 
&  (16*g3p2*Td)/3._dp + SPlamxxClam*Td + 3*TrYdadjYd*Td + TrYeadjYe*Td

 
 
If (TwoLoopRGE) Then 
betaTd2 = -6._dp*(TdadjYdYdadjYdYd) + (4*g1p2*TdadjYuYu)/5._dp - SPlamxxClam*TdadjYuYu -        & 
&  4._dp*(TdadjYuYuadjYdYd) - 2._dp*(TdadjYuYuadjYuYu) + TdadjYdYd*(6._dp*(g1p2)/5._dp + & 
&  12._dp*(g2p2) - 5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)) - 3*TdadjYuYu*TrYuadjYu -& 
&  TdadjYuYu*TrYvadjYv - (574*g1p4*M1*Yd)/45._dp - 2*g1p2*g2p2*M1*Yd - (16*g1p2*g3p2*M1*Yd)/9._dp -& 
&  (16*g1p2*g3p2*M3*Yd)/9._dp - 16*g2p2*g3p2*M3*Yd + (64*g3p4*M3*Yd)/9._dp -             & 
&  2*g1p2*g2p2*M2*Yd - 30*g2p4*M2*Yd - 16*g2p2*g3p2*M2*Yd - 4*SPadjYvYvClamxxTlam*Yd -   & 
&  4*SPClamxxTpTvCYvlam*Yd - 12*SPClamxxTlam*SPlamxxClam*Yd - (4*g1p2*TradjYdTd*Yd)/5._dp +& 
&  32*g3p2*TradjYdTd*Yd + (12*g1p2*TradjYeTe*Yd)/5._dp - 2*TradjYeTeCYvTpYv*Yd -         & 
&  6*SPlamxxClam*TradjYuTu*Yd - 2*TradjYvTpYeCYeTv*Yd - 2*SPlamxxClam*TradjYvTv*Yd +     & 
&  (4*g1p2*M1*TrYdadjYd*Yd)/5._dp - 32*g3p2*M3*TrYdadjYd*Yd - 36*TrYdadjYdTdadjYd*Yd -   & 
&  6*TrYdadjYuTuadjYd*Yd - (12*g1p2*M1*TrYeadjYe*Yd)/5._dp - 12*TrYeadjYeTeadjYe*Yd -    & 
&  6*TrYuadjYdTdadjYu*Yd - 6*SPClamxxTlam*TrYuadjYu*Yd - 2*SPClamxxTlam*TrYvadjYv*Yd +   & 
&  (6*g1p2*YdadjYdTd)/5._dp + 6*g2p2*YdadjYdTd - 4*SPlamxxClam*YdadjYdTd -               & 
&  12*TrYdadjYd*YdadjYdTd - 4*TrYeadjYe*YdadjYdTd - 8._dp*(YdadjYdTdadjYdYd) -           & 
&  (8*g1p2*M1*YdadjYdYd)/5._dp - 12*g2p2*M2*YdadjYdYd - 6*SPClamxxTlam*YdadjYdYd -       & 
&  18*TradjYdTd*YdadjYdYd - 6*TradjYeTe*YdadjYdYd - 6._dp*(YdadjYdYdadjYdTd) +           & 
&  (8*g1p2*YdadjYuTu)/5._dp - 2*SPlamxxClam*YdadjYuTu - 6*TrYuadjYu*YdadjYuTu -          & 
&  2*TrYvadjYv*YdadjYuTu - 4._dp*(YdadjYuTuadjYdYd) - 4._dp*(YdadjYuTuadjYuYu) -         & 
&  (8*g1p2*M1*YdadjYuYu)/5._dp - 2*SPClamxxTlam*YdadjYuYu - 6*TradjYuTu*YdadjYuYu -      & 
&  2*TradjYvTv*YdadjYuYu - 2._dp*(YdadjYuYuadjYdTd) - 4._dp*(YdadjYuYuadjYuTu) +         & 
&  (287*g1p4*Td)/90._dp + g1p2*g2p2*Td + (15*g2p4*Td)/2._dp + (8*g1p2*g3p2*Td)/9._dp  
betaTd2 =  betaTd2+ 8*g2p2*g3p2*Td - (16*g3p4*Td)/9._dp - 2*SPlamxxadjYvYvClam*Td - 3*SPlamxxClamp2*Td -  & 
&  (2*g1p2*TrYdadjYd*Td)/5._dp + 16*g3p2*TrYdadjYd*Td - 9*TrYdadjYdYdadjYd*Td -          & 
&  3*TrYdadjYuYuadjYd*Td + (6*g1p2*TrYeadjYe*Td)/5._dp - 3*TrYeadjYeYeadjYe*Td -         & 
&  3*SPlamxxClam*TrYuadjYu*Td - SPlamxxClam*TrYvadjYv*Td - TrYvadjYvTpYeCYe*Td -         & 
&  2*Conjg(lam(1))*Td*Tpkap1Ckap1lam(1) - 4*Yd*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) -        & 
&  2*Conjg(lam(1))*Td*Tpkap1Ckap2lam(2) - 4*Yd*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) -        & 
&  2*Conjg(lam(1))*Td*Tpkap1Ckap3lam(3) - 4*Yd*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) -        & 
&  2*Conjg(lam(2))*Td*Tpkap2Ckap1lam(1) - 4*Yd*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) -        & 
&  2*Conjg(lam(2))*Td*Tpkap2Ckap2lam(2) - 4*Yd*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) -        & 
&  2*Conjg(lam(2))*Td*Tpkap2Ckap3lam(3) - 4*Yd*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) -        & 
&  2*Conjg(lam(3))*Td*Tpkap3Ckap1lam(1) - 4*Yd*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) -        & 
&  2*Conjg(lam(3))*Td*Tpkap3Ckap2lam(2) - 4*Yd*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) -        & 
&  2*Conjg(lam(3))*Td*Tpkap3Ckap3lam(3) - 4*Yd*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) -        & 
&  4*Yd*Conjg(lam(1))*TpTk1Ckap1lam(1) - 4*Yd*Conjg(lam(1))*TpTk1Ckap2lam(2) -           & 
&  4*Yd*Conjg(lam(1))*TpTk1Ckap3lam(3) - 4*Yd*Conjg(lam(2))*TpTk2Ckap1lam(1) -           & 
&  4*Yd*Conjg(lam(2))*TpTk2Ckap2lam(2) - 4*Yd*Conjg(lam(2))*TpTk2Ckap3lam(3) -           & 
&  4*Yd*Conjg(lam(3))*TpTk3Ckap1lam(1) - 4*Yd*Conjg(lam(3))*TpTk3Ckap2lam(2) -           & 
&  4*Yd*Conjg(lam(3))*TpTk3Ckap3lam(3)

 
DTd = oo16pi2*( betaTd1 + oo16pi2 * betaTd2 ) 

 
Else 
DTd = oo16pi2* betaTd1 
End If 
 
 
Call Chop(DTd) 

!-------------------- 
! Te 
!-------------------- 
 
betaTe1  = 5._dp*(TeadjYeYe) + TeCYvTpYv + (18*g1p2*M1*Ye)/5._dp + 6*g2p2*M2*Ye +     & 
&  2*SPClamxxTlam*Ye + 6*TradjYdTd*Ye + 2*TradjYeTe*Ye + 4._dp*(YeadjYeTe)               & 
&  + 2._dp*(YeCYvTpTv) - (9*g1p2*Te)/5._dp - 3*g2p2*Te + SPlamxxClam*Te + 3*TrYdadjYd*Te +& 
&  TrYeadjYe*Te

 
 
If (TwoLoopRGE) Then 
betaTe2 = -2._dp*(DyTeCYvi11Yvadjkap1kap1i21) - 2._dp*(DyTeCYvi11Yvadjkap2kap1i22) -            & 
&  2._dp*(DyTeCYvi11Yvadjkap3kap1i23) - 2._dp*(DyTeCYvi12Yvadjkap1kap2i21) -             & 
&  2._dp*(DyTeCYvi12Yvadjkap2kap2i22) - 2._dp*(DyTeCYvi12Yvadjkap3kap2i23) -             & 
&  2._dp*(DyTeCYvi13Yvadjkap1kap3i21) - 2._dp*(DyTeCYvi13Yvadjkap2kap3i22) -             & 
&  2._dp*(DyTeCYvi13Yvadjkap3kap3i23) - 2._dp*(DyTeCYvlami1YvClami2) - 4._dp*(DyYeCYvi11TvCkap1kap1i21) -& 
&  4._dp*(DyYeCYvi11TvCkap2kap1i22) - 4._dp*(DyYeCYvi11TvCkap3kap1i23) - 4._dp*(DyYeCYvi12TvCkap1kap2i21) -& 
&  4._dp*(DyYeCYvi12TvCkap2kap2i22) - 4._dp*(DyYeCYvi12TvCkap3kap2i23) - 4._dp*(DyYeCYvi13TvCkap1kap3i21) -& 
&  4._dp*(DyYeCYvi13TvCkap2kap3i22) - 4._dp*(DyYeCYvi13TvCkap3kap3i23) - 4._dp*(DyYeCYvlami1TvClami2) -& 
&  4._dp*(DyYeCYvTlami1YvClami2) - (6*g1p2*TeadjYeYe)/5._dp + 12*g2p2*TeadjYeYe -        & 
&  5*SPlamxxClam*TeadjYeYe - 6._dp*(TeadjYeYeadjYeYe) - SPlamxxClam*TeCYvTpYv -          & 
&  4._dp*(TeCYvTpYvadjYeYe) - 2._dp*(TeCYvTpYvCYvTpYv) - 15*TeadjYeYe*TrYdadjYd -        & 
&  5*TeadjYeYe*TrYeadjYe - 3*TeCYvTpYv*TrYuadjYu - TeCYvTpYv*TrYvadjYv - 54*g1p4*M1*Ye - & 
&  (18*g1p2*g2p2*M1*Ye)/5._dp - (18*g1p2*g2p2*M2*Ye)/5._dp - 30*g2p4*M2*Ye -             & 
&  4*SPadjYvYvClamxxTlam*Ye - 4*SPClamxxTpTvCYvlam*Ye - 12*SPClamxxTlam*SPlamxxClam*Ye - & 
&  (4*g1p2*TradjYdTd*Ye)/5._dp + 32*g3p2*TradjYdTd*Ye + (12*g1p2*TradjYeTe*Ye)/5._dp -   & 
&  2*TradjYeTeCYvTpYv*Ye - 6*SPlamxxClam*TradjYuTu*Ye - 2*TradjYvTpYeCYeTv*Ye -          & 
&  2*SPlamxxClam*TradjYvTv*Ye + (4*g1p2*M1*TrYdadjYd*Ye)/5._dp - 32*g3p2*M3*TrYdadjYd*Ye -& 
&  36*TrYdadjYdTdadjYd*Ye - 6*TrYdadjYuTuadjYd*Ye - (12*g1p2*M1*TrYeadjYe*Ye)/5._dp -    & 
&  12*TrYeadjYeTeadjYe*Ye - 6*TrYuadjYdTdadjYu*Ye - 6*SPClamxxTlam*TrYuadjYu*Ye -        & 
&  2*SPClamxxTlam*TrYvadjYv*Ye + (6*g1p2*YeadjYeTe)/5._dp + 6*g2p2*YeadjYeTe -           & 
&  4*SPlamxxClam*YeadjYeTe - 12*TrYdadjYd*YeadjYeTe - 4*TrYeadjYe*YeadjYeTe  
betaTe2 =  betaTe2- 8._dp*(YeadjYeTeadjYeYe) - 12*g2p2*M2*YeadjYeYe - 6*SPClamxxTlam*YeadjYeYe -          & 
&  18*TradjYdTd*YeadjYeYe - 6*TradjYeTe*YeadjYeYe - 6._dp*(YeadjYeYeadjYeTe) -           & 
&  4._dp*(YeCYvTk1adjkap1TpYv) - 4._dp*(YeCYvTk2adjkap2TpYv) - 4._dp*(YeCYvTk3adjkap3TpYv) -& 
&  2*SPlamxxClam*YeCYvTpTv - 6*TrYuadjYu*YeCYvTpTv - 2*TrYvadjYv*YeCYvTpTv -             & 
&  4._dp*(YeCYvTpTvadjYeYe) - 4._dp*(YeCYvTpTvCYvTpYv) - 2*SPClamxxTlam*YeCYvTpYv -      & 
&  6*TradjYuTu*YeCYvTpYv - 2*TradjYvTv*YeCYvTpYv - 2._dp*(YeCYvTpYvadjYeTe) -            & 
&  4._dp*(YeCYvTpYvCYvTpTv) + (27*g1p4*Te)/2._dp + (9*g1p2*g2p2*Te)/5._dp +              & 
&  (15*g2p4*Te)/2._dp - 2*SPlamxxadjYvYvClam*Te - 3*SPlamxxClamp2*Te - (2*g1p2*TrYdadjYd*Te)/5._dp +& 
&  16*g3p2*TrYdadjYd*Te - 9*TrYdadjYdYdadjYd*Te - 3*TrYdadjYuYuadjYd*Te + (6*g1p2*TrYeadjYe*Te)/5._dp -& 
&  3*TrYeadjYeYeadjYe*Te - 3*SPlamxxClam*TrYuadjYu*Te - SPlamxxClam*TrYvadjYv*Te -       & 
&  TrYvadjYvTpYeCYe*Te - 2*Conjg(lam(1))*Te*Tpkap1Ckap1lam(1) - 4*Ye*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) -& 
&  2*Conjg(lam(1))*Te*Tpkap1Ckap2lam(2) - 4*Ye*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) -        & 
&  2*Conjg(lam(1))*Te*Tpkap1Ckap3lam(3) - 4*Ye*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) -        & 
&  2*Conjg(lam(2))*Te*Tpkap2Ckap1lam(1) - 4*Ye*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) -        & 
&  2*Conjg(lam(2))*Te*Tpkap2Ckap2lam(2) - 4*Ye*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) -        & 
&  2*Conjg(lam(2))*Te*Tpkap2Ckap3lam(3) - 4*Ye*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) -        & 
&  2*Conjg(lam(3))*Te*Tpkap3Ckap1lam(1) - 4*Ye*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) -        & 
&  2*Conjg(lam(3))*Te*Tpkap3Ckap2lam(2) - 4*Ye*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) -        & 
&  2*Conjg(lam(3))*Te*Tpkap3Ckap3lam(3) - 4*Ye*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) -        & 
&  4*Ye*Conjg(lam(1))*TpTk1Ckap1lam(1) - 4*Ye*Conjg(lam(1))*TpTk1Ckap2lam(2) -           & 
&  4*Ye*Conjg(lam(1))*TpTk1Ckap3lam(3) - 4*Ye*Conjg(lam(2))*TpTk2Ckap1lam(1)  
betaTe2 =  betaTe2- 4*Ye*Conjg(lam(2))*TpTk2Ckap2lam(2) - 4*Ye*Conjg(lam(2))*TpTk2Ckap3lam(3) -           & 
&  4*Ye*Conjg(lam(3))*TpTk3Ckap1lam(1) - 4*Ye*Conjg(lam(3))*TpTk3Ckap2lam(2) -           & 
&  4*Ye*Conjg(lam(3))*TpTk3Ckap3lam(3)

 
DTe = oo16pi2*( betaTe1 + oo16pi2 * betaTe2 ) 

 
Else 
DTe = oo16pi2* betaTe1 
End If 
 
 
Call Chop(DTe) 

!-------------------- 
! Tlam 
!-------------------- 
 
Do i1 = 1,3
betaTlam1(i1) = 4*Tk1adjkap1lam(i1) + 4*Tk2adjkap2lam(i1) + 4*Tk3adjkap3lam(i1) + 5*TpTvCYvlam(i1) +  & 
&  4*TpYvCYvTlam(i1) + (6*g1p2*M1*lam(i1))/5._dp + 6*g2p2*M2*lam(i1) + 6*SPClamxxTlam*lam(i1) +& 
&  6*TradjYdTd*lam(i1) + 2*TradjYeTe*lam(i1) + 6*TradjYuTu*lam(i1) + 2*TradjYvTv*lam(i1) +& 
&  2*kap1Ckap1(i1,1)*Tlam(1) + 2*kap2Ckap1(i1,2)*Tlam(1) + 2*kap3Ckap1(i1,               & 
& 3)*Tlam(1) + 2*kap1Ckap2(i1,1)*Tlam(2) + 2*kap2Ckap2(i1,2)*Tlam(2) + 2*kap3Ckap2(i1,   & 
& 3)*Tlam(2) + 2*kap1Ckap3(i1,1)*Tlam(3) + 2*kap2Ckap3(i1,2)*Tlam(3) + 2*kap3Ckap3(i1,   & 
& 3)*Tlam(3) - (3*g1p2*Tlam(i1))/5._dp - 3*g2p2*Tlam(i1) + 6*SPlamxxClam*Tlam(i1) +      & 
&  3*TrYdadjYd*Tlam(i1) + TrYeadjYe*Tlam(i1) + 3*TrYuadjYu*Tlam(i1) + TrYvadjYv*Tlam(i1)
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
betaTlam2(i1) = -16*TradjYvTv*kap1Ckap1lam(i1) - 16*TrCkap1Tk1*kap1Ckap1lam(i1) - 16*Ckap2Tk1(1,      & 
& 2)*kap1Ckap1lam(i1) - 16*Ckap3Tk1(1,3)*kap1Ckap1lam(i1) - 16*TrCkap1Tk2*kap1Ckap2lam(i1) -& 
&  16*adjYvTv(1,2)*kap1Ckap2lam(i1) - 16*Ckap2Tk2(1,2)*kap1Ckap2lam(i1) - 16*Ckap3Tk2(1, & 
& 3)*kap1Ckap2lam(i1) - 16*TrCkap1Tk3*kap1Ckap3lam(i1) - 16*adjYvTv(1,3)*kap1Ckap3lam(i1) -& 
&  16*Ckap2Tk3(1,2)*kap1Ckap3lam(i1) - 16*Ckap3Tk3(1,3)*kap1Ckap3lam(i1) -               & 
&  16*TrCkap2Tk1*kap2Ckap1lam(i1) - 16*adjYvTv(2,1)*kap2Ckap1lam(i1) - 16*Ckap1Tk1(2,    & 
& 1)*kap2Ckap1lam(i1) - 16*Ckap3Tk1(2,3)*kap2Ckap1lam(i1) - 16*TradjYvTv*kap2Ckap2lam(i1) -& 
&  16*TrCkap2Tk2*kap2Ckap2lam(i1) - 16*Ckap1Tk2(2,1)*kap2Ckap2lam(i1) - 16*Ckap3Tk2(2,   & 
& 3)*kap2Ckap2lam(i1) - 16*TrCkap2Tk3*kap2Ckap3lam(i1) - 16*adjYvTv(2,3)*kap2Ckap3lam(i1) -& 
&  16*Ckap1Tk3(2,1)*kap2Ckap3lam(i1) - 16*Ckap3Tk3(2,3)*kap2Ckap3lam(i1) -               & 
&  16*TrCkap3Tk1*kap3Ckap1lam(i1) - 16*adjYvTv(3,1)*kap3Ckap1lam(i1) - 16*Ckap1Tk1(3,    & 
& 1)*kap3Ckap1lam(i1) - 16*Ckap2Tk1(3,2)*kap3Ckap1lam(i1) - 16*TrCkap3Tk2*kap3Ckap2lam(i1) -& 
&  16*adjYvTv(3,2)*kap3Ckap2lam(i1) - 16*Ckap1Tk2(3,1)*kap3Ckap2lam(i1) - 16*Ckap2Tk2(3, & 
& 2)*kap3Ckap2lam(i1) - 16*TradjYvTv*kap3Ckap3lam(i1) - 16*TrCkap3Tk3*kap3Ckap3lam(i1) - & 
&  16*Ckap1Tk3(3,1)*kap3Ckap3lam(i1) - 16*Ckap2Tk3(3,2)*kap3Ckap3lam(i1) -               & 
&  16*Tk1adjYvYvCkap1lam(i1) - 16*SPlamxxadjkap1lam*Tk1Clam(i1) - 16*Tk2adjYvYvCkap2lam(i1) -& 
&  16*SPlamxxadjkap2lam*Tk2Clam(i1) - 16*Tk3adjYvYvCkap3lam(i1) - 16*SPlamxxadjkap3lam*Tk3Clam(i1) -& 
&  2*adjYvTv(1,i1)*Tpkap1Ckap1lam(1) - 16*Ckap1lam(1)*Tpkap1Ckap1TpTk1(1,i1) -           & 
&  16*Ckap2lam(1)*Tpkap1Ckap1TpTk2(1,i1) - 16*Ckap3lam(1)*Tpkap1Ckap1TpTk3(1,            & 
& i1) - 2*adjYvTv(1,i1)*Tpkap1Ckap2lam(2) - 16*Ckap1lam(1)*Tpkap1Ckap2TpTk1(2,           & 
& i1) - 16*Ckap2lam(1)*Tpkap1Ckap2TpTk2(2,i1) - 16*Ckap3lam(1)*Tpkap1Ckap2TpTk3(2,       & 
& i1) - 2*adjYvTv(1,i1)*Tpkap1Ckap3lam(3) - 16*Ckap1lam(1)*Tpkap1Ckap3TpTk1(3,           & 
& i1) - 16*Ckap2lam(1)*Tpkap1Ckap3TpTk2(3,i1) - 16*Ckap3lam(1)*Tpkap1Ckap3TpTk3(3,       & 
& i1) - 2*adjYvTv(2,i1)*Tpkap2Ckap1lam(1) - 16*Ckap1lam(2)*Tpkap2Ckap1TpTk1(1,           & 
& i1) - 16*Ckap2lam(2)*Tpkap2Ckap1TpTk2(1,i1) - 16*Ckap3lam(2)*Tpkap2Ckap1TpTk3(1,       & 
& i1) - 2*adjYvTv(2,i1)*Tpkap2Ckap2lam(2) - 16*Ckap1lam(2)*Tpkap2Ckap2TpTk1(2,           & 
& i1) - 16*Ckap2lam(2)*Tpkap2Ckap2TpTk2(2,i1) - 16*Ckap3lam(2)*Tpkap2Ckap2TpTk3(2,       & 
& i1) - 2*adjYvTv(2,i1)*Tpkap2Ckap3lam(3) - 16*Ckap1lam(2)*Tpkap2Ckap3TpTk1(3,           & 
& i1) - 16*Ckap2lam(2)*Tpkap2Ckap3TpTk2(3,i1) - 16*Ckap3lam(2)*Tpkap2Ckap3TpTk3(3,       & 
& i1) - 2*adjYvTv(3,i1)*Tpkap3Ckap1lam(1) - 16*Ckap1lam(3)*Tpkap3Ckap1TpTk1(1,           & 
& i1) - 16*Ckap2lam(3)*Tpkap3Ckap1TpTk2(1,i1) - 16*Ckap3lam(3)*Tpkap3Ckap1TpTk3(1,       & 
& i1) - 2*adjYvTv(3,i1)*Tpkap3Ckap2lam(2) - 16*Ckap1lam(3)*Tpkap3Ckap2TpTk1(2,           & 
& i1) - 16*Ckap2lam(3)*Tpkap3Ckap2TpTk2(2,i1) - 16*Ckap3lam(3)*Tpkap3Ckap2TpTk3(2,       & 
& i1) - 2*adjYvTv(3,i1)*Tpkap3Ckap3lam(3) - 16*Ckap1lam(3)*Tpkap3Ckap3TpTk1(3,           & 
& i1) - 16*Ckap2lam(3)*Tpkap3Ckap3TpTk2(3,i1) - 16*Ckap3lam(3)*Tpkap3Ckap3TpTk3(3,       & 
& i1) - 3*TpTvadjYeYeCYvlam(i1) + (12*g1p2*TpTvCYvlam(i1))/5._dp + 12*g2p2*TpTvCYvlam(i1)  
betaTlam2(i1) =  betaTlam2(i1)- 11*SPlamxxClam*TpTvCYvlam(i1) - 15*TrYuadjYu*TpTvCYvlam(i1) - 5*TrYvadjYv*TpTvCYvlam(i1) -& 
&  6*TpTvCYvTpYvCYvlam(i1) - 2*TpYvadjYeTeCYvlam(i1) - 4*Tpkap1adjkap1Tlam(1)*TpYvCYv(i1,& 
& 1) - 4*Tpkap1adjkap2Tlam(2)*TpYvCYv(i1,1) - 4*Tpkap1adjkap3Tlam(3)*TpYvCYv(i1,         & 
& 1) - 4*Tpkap2adjkap1Tlam(1)*TpYvCYv(i1,2) - 4*Tpkap2adjkap2Tlam(2)*TpYvCYv(i1,         & 
& 2) - 4*Tpkap2adjkap3Tlam(3)*TpYvCYv(i1,2) - 4*Tpkap3adjkap1Tlam(1)*TpYvCYv(i1,         & 
& 3) - 4*Tpkap3adjkap2Tlam(2)*TpYvCYv(i1,3) - 4*Tpkap3adjkap3Tlam(3)*TpYvCYv(i1,         & 
& 3) - (12*g1p2*M1*TpYvCYvlam(i1))/5._dp - 12*g2p2*M2*TpYvCYvlam(i1) - 12*SPClamxxTlam*TpYvCYvlam(i1) -& 
&  18*TradjYuTu*TpYvCYvlam(i1) - 6*TradjYvTv*TpYvCYvlam(i1) - 4*TpYvCYvTk1adjkap1lam(i1) -& 
&  4*TpYvCYvTk2adjkap2lam(i1) - 4*TpYvCYvTk3adjkap3lam(i1) + (6*g1p2*TpYvCYvTlam(i1))/5._dp +& 
&  6*g2p2*TpYvCYvTlam(i1) - 12*SPlamxxClam*TpYvCYvTlam(i1) - 12*TrYuadjYu*TpYvCYvTlam(i1) -& 
&  4*TrYvadjYv*TpYvCYvTlam(i1) - 8*TpYvCYvTpTvCYvlam(i1) - 6*TpYvCYvTpYvCYvTlam(i1) -    & 
&  (414*g1p4*M1*lam(i1))/25._dp - (18*g1p2*g2p2*M1*lam(i1))/5._dp - (18*g1p2*g2p2*M2*lam(i1))/5._dp -& 
&  30*g2p4*M2*lam(i1) - 18*SPadjYvYvClamxxTlam*lam(i1) + (6*g1p2*SPClamxxTlam*lam(i1))/5._dp +& 
&  6*g2p2*SPClamxxTlam*lam(i1) - 20*SPClamxxTpTvCYvlam*lam(i1) - (12*g1p2*M1*SPlamxxClam*lam(i1))/5._dp -& 
&  12*g2p2*M2*SPlamxxClam*lam(i1) - 36*SPClamxxTlam*SPlamxxClam*lam(i1) - (4*g1p2*TradjYdTd*lam(i1))/5._dp +& 
&  32*g3p2*TradjYdTd*lam(i1) - 18*SPlamxxClam*TradjYdTd*lam(i1) + (12*g1p2*TradjYeTe*lam(i1))/5._dp -& 
&  6*SPlamxxClam*TradjYeTe*lam(i1) - 4*TradjYeTeCYvTpYv*lam(i1) + (8*g1p2*TradjYuTu*lam(i1))/5._dp +& 
&  32*g3p2*TradjYuTu*lam(i1) - 18*SPlamxxClam*TradjYuTu*lam(i1) - 4*TradjYvTpYeCYeTv*lam(i1) -& 
&  6*SPlamxxClam*TradjYvTv*lam(i1) - 4*TradjYvTvadjkap1kap1*lam(i1) - 4*TradjYvTvadjkap2kap2*lam(i1) -& 
&  4*TradjYvTvadjkap3kap3*lam(i1) + (4*g1p2*M1*TrYdadjYd*lam(i1))/5._dp - 32*g3p2*M3*TrYdadjYd*lam(i1) -& 
&  12*SPClamxxTlam*TrYdadjYd*lam(i1) - 36*TrYdadjYdTdadjYd*lam(i1) - 12*TrYdadjYuTuadjYd*lam(i1)  
betaTlam2(i1) =  betaTlam2(i1)- (12*g1p2*M1*TrYeadjYe*lam(i1))/5._dp - 4*SPClamxxTlam*TrYeadjYe*lam(i1) -             & 
&  12*TrYeadjYeTeadjYe*lam(i1) - 12*TrYuadjYdTdadjYu*lam(i1) - (8*g1p2*M1*TrYuadjYu*lam(i1))/5._dp -& 
&  32*g3p2*M3*TrYuadjYu*lam(i1) - 12*SPClamxxTlam*TrYuadjYu*lam(i1) - 36*TrYuadjYuTuadjYu*lam(i1) -& 
&  4*SPClamxxTlam*TrYvadjYv*lam(i1) - 12*TrYvadjYvTvadjYv*lam(i1) - 4*TrYvCkap1TpTk1adjYv*lam(i1) -& 
&  4*TrYvCkap2TpTk2adjYv*lam(i1) - 4*TrYvCkap3TpTk3adjYv*lam(i1) - 4*adjYvTvadjkap1kap2(2,& 
& 1)*lam(i1) - 4*adjYvTvadjkap1kap3(3,1)*lam(i1) - 4*adjYvTvadjkap2kap1(1,               & 
& 2)*lam(i1) - 4*adjYvTvadjkap2kap3(3,2)*lam(i1) - 4*adjYvTvadjkap3kap1(1,               & 
& 3)*lam(i1) - 4*adjYvTvadjkap3kap2(2,3)*lam(i1) - 8*Conjg(lam(1))*Tpkap1Ckap1Tlam(1)*lam(i1) -& 
&  8*Conjg(lam(1))*Tpkap1Ckap2Tlam(2)*lam(i1) - 8*Conjg(lam(1))*Tpkap1Ckap3Tlam(3)*lam(i1) -& 
&  8*Conjg(lam(2))*Tpkap2Ckap1Tlam(1)*lam(i1) - 8*Conjg(lam(2))*Tpkap2Ckap2Tlam(2)*lam(i1) -& 
&  8*Conjg(lam(2))*Tpkap2Ckap3Tlam(3)*lam(i1) - 8*Conjg(lam(3))*Tpkap3Ckap1Tlam(1)*lam(i1) -& 
&  8*Conjg(lam(3))*Tpkap3Ckap2Tlam(2)*lam(i1) - 8*Conjg(lam(3))*Tpkap3Ckap3Tlam(3)*lam(i1) -& 
&  8*Conjg(lam(1))*TpTk1Ckap1lam(1)*lam(i1) - 8*Conjg(lam(1))*TpTk1Ckap2lam(2)*lam(i1) - & 
&  8*Conjg(lam(1))*TpTk1Ckap3lam(3)*lam(i1) - 8*Conjg(lam(2))*TpTk2Ckap1lam(1)*lam(i1) - & 
&  8*Conjg(lam(2))*TpTk2Ckap2lam(2)*lam(i1) - 8*Conjg(lam(2))*TpTk2Ckap3lam(3)*lam(i1) - & 
&  8*Conjg(lam(3))*TpTk3Ckap1lam(1)*lam(i1) - 8*Conjg(lam(3))*TpTk3Ckap2lam(2)*lam(i1) - & 
&  8*Conjg(lam(3))*TpTk3Ckap3lam(3)*lam(i1) - 8*kap1adjYvYvCkap1(i1,1)*Tlam(1) -         & 
&  8*TrCkap1kap1*kap1Ckap1(i1,1)*Tlam(1) - 8*Ckap2kap1(1,2)*kap1Ckap1(i1,1)*Tlam(1) -    & 
&  8*Ckap3kap1(1,3)*kap1Ckap1(i1,1)*Tlam(1) - 8*TrCkap1kap2*kap1Ckap1(i1,2)*Tlam(1) -    & 
&  8*Ckap2kap2(1,2)*kap1Ckap1(i1,2)*Tlam(1) - 8*Ckap3kap2(1,3)*kap1Ckap1(i1,             & 
& 2)*Tlam(1) - 8*TrCkap1kap3*kap1Ckap1(i1,3)*Tlam(1) - 8*Ckap2kap3(1,2)*kap1Ckap1(i1,    & 
& 3)*Tlam(1) - 8*Ckap3kap3(1,3)*kap1Ckap1(i1,3)*Tlam(1) - 24*Conjg(lam(1))*kap1Ckap1lam(i1)*Tlam(1)  
betaTlam2(i1) =  betaTlam2(i1)- 8*kap2adjYvYvCkap1(i1,2)*Tlam(1) - 8*TrCkap2kap1*kap2Ckap1(i1,1)*Tlam(1) -            & 
&  8*Ckap1kap1(2,1)*kap2Ckap1(i1,1)*Tlam(1) - 8*Ckap3kap1(2,3)*kap2Ckap1(i1,             & 
& 1)*Tlam(1) - 8*TrCkap2kap2*kap2Ckap1(i1,2)*Tlam(1) - 8*Ckap1kap2(2,1)*kap2Ckap1(i1,    & 
& 2)*Tlam(1) - 8*Ckap3kap2(2,3)*kap2Ckap1(i1,2)*Tlam(1) - 8*TrCkap2kap3*kap2Ckap1(i1,    & 
& 3)*Tlam(1) - 8*Ckap1kap3(2,1)*kap2Ckap1(i1,3)*Tlam(1) - 8*Ckap3kap3(2,3)*kap2Ckap1(i1, & 
& 3)*Tlam(1) - 24*Conjg(lam(2))*kap2Ckap1lam(i1)*Tlam(1) - 8*kap3adjYvYvCkap1(i1,        & 
& 3)*Tlam(1) - 8*TrCkap3kap1*kap3Ckap1(i1,1)*Tlam(1) - 8*Ckap1kap1(3,1)*kap3Ckap1(i1,    & 
& 1)*Tlam(1) - 8*Ckap2kap1(3,2)*kap3Ckap1(i1,1)*Tlam(1) - 8*TrCkap3kap2*kap3Ckap1(i1,    & 
& 2)*Tlam(1) - 8*Ckap1kap2(3,1)*kap3Ckap1(i1,2)*Tlam(1) - 8*Ckap2kap2(3,2)*kap3Ckap1(i1, & 
& 2)*Tlam(1) - 8*TrCkap3kap3*kap3Ckap1(i1,3)*Tlam(1) - 8*Ckap1kap3(3,1)*kap3Ckap1(i1,    & 
& 3)*Tlam(1) - 8*Ckap2kap3(3,2)*kap3Ckap1(i1,3)*Tlam(1) - 24*Conjg(lam(3))*kap3Ckap1lam(i1)*Tlam(1) -& 
&  8*kap1adjYvYvCkap2(i1,1)*Tlam(2) - 8*TrCkap1kap1*kap1Ckap2(i1,1)*Tlam(2) -            & 
&  8*Ckap2kap1(1,2)*kap1Ckap2(i1,1)*Tlam(2) - 8*Ckap3kap1(1,3)*kap1Ckap2(i1,             & 
& 1)*Tlam(2) - 8*TrCkap1kap2*kap1Ckap2(i1,2)*Tlam(2) - 8*Ckap2kap2(1,2)*kap1Ckap2(i1,    & 
& 2)*Tlam(2) - 8*Ckap3kap2(1,3)*kap1Ckap2(i1,2)*Tlam(2) - 8*TrCkap1kap3*kap1Ckap2(i1,    & 
& 3)*Tlam(2) - 8*Ckap2kap3(1,2)*kap1Ckap2(i1,3)*Tlam(2) - 8*Ckap3kap3(1,3)*kap1Ckap2(i1, & 
& 3)*Tlam(2) - 24*Conjg(lam(1))*kap1Ckap2lam(i1)*Tlam(2) - 8*kap2adjYvYvCkap2(i1,        & 
& 2)*Tlam(2) - 8*TrCkap2kap1*kap2Ckap2(i1,1)*Tlam(2) - 8*Ckap1kap1(2,1)*kap2Ckap2(i1,    & 
& 1)*Tlam(2) - 8*Ckap3kap1(2,3)*kap2Ckap2(i1,1)*Tlam(2) - 8*TrCkap2kap2*kap2Ckap2(i1,    & 
& 2)*Tlam(2) - 8*Ckap1kap2(2,1)*kap2Ckap2(i1,2)*Tlam(2) - 8*Ckap3kap2(2,3)*kap2Ckap2(i1, & 
& 2)*Tlam(2) - 8*TrCkap2kap3*kap2Ckap2(i1,3)*Tlam(2) - 8*Ckap1kap3(2,1)*kap2Ckap2(i1,    & 
& 3)*Tlam(2) - 8*Ckap3kap3(2,3)*kap2Ckap2(i1,3)*Tlam(2) - 24*Conjg(lam(2))*kap2Ckap2lam(i1)*Tlam(2)  
betaTlam2(i1) =  betaTlam2(i1)- 8*kap3adjYvYvCkap2(i1,3)*Tlam(2) - 8*TrCkap3kap1*kap3Ckap2(i1,1)*Tlam(2) -            & 
&  8*Ckap1kap1(3,1)*kap3Ckap2(i1,1)*Tlam(2) - 8*Ckap2kap1(3,2)*kap3Ckap2(i1,             & 
& 1)*Tlam(2) - 8*TrCkap3kap2*kap3Ckap2(i1,2)*Tlam(2) - 8*Ckap1kap2(3,1)*kap3Ckap2(i1,    & 
& 2)*Tlam(2) - 8*Ckap2kap2(3,2)*kap3Ckap2(i1,2)*Tlam(2) - 8*TrCkap3kap3*kap3Ckap2(i1,    & 
& 3)*Tlam(2) - 8*Ckap1kap3(3,1)*kap3Ckap2(i1,3)*Tlam(2) - 8*Ckap2kap3(3,2)*kap3Ckap2(i1, & 
& 3)*Tlam(2) - 24*Conjg(lam(3))*kap3Ckap2lam(i1)*Tlam(2) - 8*kap1adjYvYvCkap3(i1,        & 
& 1)*Tlam(3) - 8*TrCkap1kap1*kap1Ckap3(i1,1)*Tlam(3) - 8*Ckap2kap1(1,2)*kap1Ckap3(i1,    & 
& 1)*Tlam(3) - 8*Ckap3kap1(1,3)*kap1Ckap3(i1,1)*Tlam(3) - 8*TrCkap1kap2*kap1Ckap3(i1,    & 
& 2)*Tlam(3) - 8*Ckap2kap2(1,2)*kap1Ckap3(i1,2)*Tlam(3) - 8*Ckap3kap2(1,3)*kap1Ckap3(i1, & 
& 2)*Tlam(3) - 8*TrCkap1kap3*kap1Ckap3(i1,3)*Tlam(3) - 8*Ckap2kap3(1,2)*kap1Ckap3(i1,    & 
& 3)*Tlam(3) - 8*Ckap3kap3(1,3)*kap1Ckap3(i1,3)*Tlam(3) - 24*Conjg(lam(1))*kap1Ckap3lam(i1)*Tlam(3) -& 
&  8*kap2adjYvYvCkap3(i1,2)*Tlam(3) - 8*TrCkap2kap1*kap2Ckap3(i1,1)*Tlam(3) -            & 
&  8*Ckap1kap1(2,1)*kap2Ckap3(i1,1)*Tlam(3) - 8*Ckap3kap1(2,3)*kap2Ckap3(i1,             & 
& 1)*Tlam(3) - 8*TrCkap2kap2*kap2Ckap3(i1,2)*Tlam(3) - 8*Ckap1kap2(2,1)*kap2Ckap3(i1,    & 
& 2)*Tlam(3) - 8*Ckap3kap2(2,3)*kap2Ckap3(i1,2)*Tlam(3) - 8*TrCkap2kap3*kap2Ckap3(i1,    & 
& 3)*Tlam(3) - 8*Ckap1kap3(2,1)*kap2Ckap3(i1,3)*Tlam(3) - 8*Ckap3kap3(2,3)*kap2Ckap3(i1, & 
& 3)*Tlam(3) - 24*Conjg(lam(2))*kap2Ckap3lam(i1)*Tlam(3) - 8*kap3adjYvYvCkap3(i1,        & 
& 3)*Tlam(3) - 8*TrCkap3kap1*kap3Ckap3(i1,1)*Tlam(3) - 8*Ckap1kap1(3,1)*kap3Ckap3(i1,    & 
& 1)*Tlam(3) - 8*Ckap2kap1(3,2)*kap3Ckap3(i1,1)*Tlam(3) - 8*TrCkap3kap2*kap3Ckap3(i1,    & 
& 2)*Tlam(3) - 8*Ckap1kap2(3,1)*kap3Ckap3(i1,2)*Tlam(3) - 8*Ckap2kap2(3,2)*kap3Ckap3(i1, & 
& 2)*Tlam(3) - 8*TrCkap3kap3*kap3Ckap3(i1,3)*Tlam(3) - 8*Ckap1kap3(3,1)*kap3Ckap3(i1,    & 
& 3)*Tlam(3) - 8*Ckap2kap3(3,2)*kap3Ckap3(i1,3)*Tlam(3) - 24*Conjg(lam(3))*kap3Ckap3lam(i1)*Tlam(3)  
betaTlam2(i1) =  betaTlam2(i1)+ (207*g1p4*Tlam(i1))/50._dp + (9*g1p2*g2p2*Tlam(i1))/5._dp + (15*g2p4*Tlam(i1))/2._dp -& 
&  12*SPlamxxadjYvYvClam*Tlam(i1) + (12*g1p2*SPlamxxClam*Tlam(i1))/5._dp +               & 
&  12*g2p2*SPlamxxClam*Tlam(i1) - 14*SPlamxxClamp2*Tlam(i1) - (2*g1p2*TrYdadjYd*Tlam(i1))/5._dp +& 
&  16*g3p2*TrYdadjYd*Tlam(i1) - 15*SPlamxxClam*TrYdadjYd*Tlam(i1) - 9*TrYdadjYdYdadjYd*Tlam(i1) -& 
&  6*TrYdadjYuYuadjYd*Tlam(i1) + (6*g1p2*TrYeadjYe*Tlam(i1))/5._dp - 5*SPlamxxClam*TrYeadjYe*Tlam(i1) -& 
&  3*TrYeadjYeYeadjYe*Tlam(i1) + (4*g1p2*TrYuadjYu*Tlam(i1))/5._dp + 16*g3p2*TrYuadjYu*Tlam(i1) -& 
&  15*SPlamxxClam*TrYuadjYu*Tlam(i1) - 9*TrYuadjYuYuadjYu*Tlam(i1) - 5*SPlamxxClam*TrYvadjYv*Tlam(i1) -& 
&  2*TrYvadjYvTpYeCYe*Tlam(i1) - 3*TrYvadjYvYvadjYv*Tlam(i1) - 2*TrYvCkap1kap1adjYv*Tlam(i1) -& 
&  2*TrYvCkap2kap2adjYv*Tlam(i1) - 2*TrYvCkap3kap3adjYv*Tlam(i1) - 2*adjYvYvCkap1kap2(2, & 
& 1)*Tlam(i1) - 2*adjYvYvCkap1kap3(3,1)*Tlam(i1) - 2*adjYvYvCkap2kap1(1,2)*Tlam(i1) -    & 
&  2*adjYvYvCkap2kap3(3,2)*Tlam(i1) - 2*adjYvYvCkap3kap1(1,3)*Tlam(i1) - 2*adjYvYvCkap3kap2(2,& 
& 3)*Tlam(i1) - 4*Conjg(lam(1))*Tpkap1Ckap1lam(1)*Tlam(i1) - 4*Conjg(lam(1))*Tpkap1Ckap2lam(2)*Tlam(i1) -& 
&  4*Conjg(lam(1))*Tpkap1Ckap3lam(3)*Tlam(i1) - 4*Conjg(lam(2))*Tpkap2Ckap1lam(1)*Tlam(i1) -& 
&  4*Conjg(lam(2))*Tpkap2Ckap2lam(2)*Tlam(i1) - 4*Conjg(lam(2))*Tpkap2Ckap3lam(3)*Tlam(i1) -& 
&  4*Conjg(lam(3))*Tpkap3Ckap1lam(1)*Tlam(i1) - 4*Conjg(lam(3))*Tpkap3Ckap2lam(2)*Tlam(i1) -& 
&  4*Conjg(lam(3))*Tpkap3Ckap3lam(3)*Tlam(i1)
End Do

 
DTlam = oo16pi2*( betaTlam1 + oo16pi2 * betaTlam2 ) 

 
Else 
DTlam = oo16pi2* betaTlam1 
End If 
 
 
Call Chop(DTlam) 

!-------------------- 
! Tv 
!-------------------- 
 
betaTv1  = 4._dp*(DyTvClami1lami2) + 2._dp*(DyTvi11kap1Ckap1i21) + 2._dp*(DyTvi11kap2Ckap1i22)& 
&  + 2._dp*(DyTvi11kap3Ckap1i23) + 2._dp*(DyTvi12kap1Ckap2i21) + 2._dp*(DyTvi12kap2Ckap2i22)& 
&  + 2._dp*(DyTvi12kap3Ckap2i23) + 2._dp*(DyTvi13kap1Ckap3i21) + 2._dp*(DyTvi13kap2Ckap3i22)& 
&  + 2._dp*(DyTvi13kap3Ckap3i23) + 5._dp*(DyYvClami1Tlami2) + 2._dp*(TpTeCYeYv)          & 
&  + TpYeCYeTv + 4._dp*(TvadjYvYv) + (6*g1p2*M1*Yv)/5._dp + 6*g2p2*M2*Yv +               & 
&  2*SPClamxxTlam*Yv + 6*TradjYuTu*Yv + 2*TradjYvTv*Yv + 5._dp*(YvadjYvTv)               & 
&  + 4._dp*(YvCkap1TpTk1) + 4._dp*(YvCkap2TpTk2) + 4._dp*(YvCkap3TpTk3) - (3*g1p2*Tv)    & 
& /5._dp - 3*g2p2*Tv + SPlamxxClam*Tv + 3*TrYuadjYu*Tv + TrYvadjYv*Tv

 
 
If (TwoLoopRGE) Then 
betaTv2 = 2._dp*(DyTpTeCYeYvClami1lami2) + 2._dp*(DyTpYeCYeTvClami1lami2) + DyTpYeCYeYvClami1Tlami2 -& 
&  6._dp*(DyTvadjYvYvClami1lami2) - 4._dp*(DyTvCkap1kap1i11TpYvCYvi21) - 4._dp*(DyTvCkap1kap2i11TpYvCYvi22) -& 
&  4._dp*(DyTvCkap1kap3i11TpYvCYvi23) - 4._dp*(DyTvCkap2kap1i12TpYvCYvi21) -             & 
&  4._dp*(DyTvCkap2kap2i12TpYvCYvi22) - 4._dp*(DyTvCkap2kap3i12TpYvCYvi23) -             & 
&  4._dp*(DyTvCkap3kap1i13TpYvCYvi21) - 4._dp*(DyTvCkap3kap2i13TpYvCYvi22) -             & 
&  4._dp*(DyTvCkap3kap3i13TpYvCYvi23) - 6._dp*(DyTvClami1TpYvCYvlami2) - 8._dp*(DyTvi11kap1adjYvYvCkap1i21) -& 
&  8._dp*(DyTvi11kap2adjYvYvCkap1i22) - 8._dp*(DyTvi11kap3adjYvYvCkap1i23) -             & 
&  8._dp*(DyTvi12kap1adjYvYvCkap2i21) - 8._dp*(DyTvi12kap2adjYvYvCkap2i22) -             & 
&  8._dp*(DyTvi12kap3adjYvYvCkap2i23) - 8._dp*(DyTvi13kap1adjYvYvCkap3i21) -             & 
&  8._dp*(DyTvi13kap2adjYvYvCkap3i22) - 8._dp*(DyTvi13kap3adjYvYvCkap3i23) -             & 
&  2._dp*(DyYvadjkap1kap1i11adjYvTv1i2) - 2._dp*(DyYvadjkap1kap2i11adjYvTv2i2) -         & 
&  2._dp*(DyYvadjkap1kap3i11adjYvTv3i2) - 2._dp*(DyYvadjkap2kap1i12adjYvTv1i2) -         & 
&  2._dp*(DyYvadjkap2kap2i12adjYvTv2i2) - 2._dp*(DyYvadjkap2kap3i12adjYvTv3i2) -         & 
&  2._dp*(DyYvadjkap3kap1i13adjYvTv1i2) - 2._dp*(DyYvadjkap3kap2i13adjYvTv2i2) -         & 
&  2._dp*(DyYvadjkap3kap3i13adjYvTv3i2) - 8._dp*(DyYvadjYvTvClami1lami2) -               & 
&  6._dp*(DyYvadjYvYvClami1Tlami2) - 16._dp*(DyYvCkap1i11kap1adjYvTvi21) -               & 
&  16._dp*(DyYvCkap1i11kap1Ckap1Tk1i21) - 16._dp*(DyYvCkap1i11kap1Ckap2Tk1i22) -         & 
&  16._dp*(DyYvCkap1i11kap1Ckap3Tk1i23) - 16._dp*(DyYvCkap1i11Tk1Ckap1kap1i21) -         & 
&  16._dp*(DyYvCkap1i11Tk1Ckap2kap1i22) - 16._dp*(DyYvCkap1i11Tk1Ckap3kap1i23) -         & 
&  16._dp*(DyYvCkap1i12kap2adjYvTvi21) - 16._dp*(DyYvCkap1i12kap2Ckap1Tk1i21) -          & 
&  16._dp*(DyYvCkap1i12kap2Ckap2Tk1i22) - 16._dp*(DyYvCkap1i12kap2Ckap3Tk1i23)  
betaTv2 =  betaTv2- 16._dp*(DyYvCkap1i12Tk1Ckap1kap2i21) - 16._dp*(DyYvCkap1i12Tk1Ckap2kap2i22) -         & 
&  16._dp*(DyYvCkap1i12Tk1Ckap3kap2i23) - 16._dp*(DyYvCkap1i13kap3adjYvTvi21) -          & 
&  16._dp*(DyYvCkap1i13kap3Ckap1Tk1i21) - 16._dp*(DyYvCkap1i13kap3Ckap2Tk1i22) -         & 
&  16._dp*(DyYvCkap1i13kap3Ckap3Tk1i23) - 16._dp*(DyYvCkap1i13Tk1Ckap1kap3i21) -         & 
&  16._dp*(DyYvCkap1i13Tk1Ckap2kap3i22) - 16._dp*(DyYvCkap1i13Tk1Ckap3kap3i23) -         & 
&  16._dp*(DyYvCkap1lami1Tk1Clami2) - 16._dp*(DyYvCkap2i11kap1adjYvTvi22) -              & 
&  16._dp*(DyYvCkap2i11kap1Ckap1Tk2i21) - 16._dp*(DyYvCkap2i11kap1Ckap2Tk2i22) -         & 
&  16._dp*(DyYvCkap2i11kap1Ckap3Tk2i23) - 16._dp*(DyYvCkap2i11Tk2Ckap1kap1i21) -         & 
&  16._dp*(DyYvCkap2i11Tk2Ckap2kap1i22) - 16._dp*(DyYvCkap2i11Tk2Ckap3kap1i23) -         & 
&  16._dp*(DyYvCkap2i12kap2adjYvTvi22) - 16._dp*(DyYvCkap2i12kap2Ckap1Tk2i21) -          & 
&  16._dp*(DyYvCkap2i12kap2Ckap2Tk2i22) - 16._dp*(DyYvCkap2i12kap2Ckap3Tk2i23) -         & 
&  16._dp*(DyYvCkap2i12Tk2Ckap1kap2i21) - 16._dp*(DyYvCkap2i12Tk2Ckap2kap2i22) -         & 
&  16._dp*(DyYvCkap2i12Tk2Ckap3kap2i23) - 16._dp*(DyYvCkap2i13kap3adjYvTvi22) -          & 
&  16._dp*(DyYvCkap2i13kap3Ckap1Tk2i21) - 16._dp*(DyYvCkap2i13kap3Ckap2Tk2i22) -         & 
&  16._dp*(DyYvCkap2i13kap3Ckap3Tk2i23) - 16._dp*(DyYvCkap2i13Tk2Ckap1kap3i21) -         & 
&  16._dp*(DyYvCkap2i13Tk2Ckap2kap3i22) - 16._dp*(DyYvCkap2i13Tk2Ckap3kap3i23) -         & 
&  16._dp*(DyYvCkap2lami1Tk2Clami2) - 16._dp*(DyYvCkap3i11kap1adjYvTvi23) -              & 
&  16._dp*(DyYvCkap3i11kap1Ckap1Tk3i21) - 16._dp*(DyYvCkap3i11kap1Ckap2Tk3i22) -         & 
&  16._dp*(DyYvCkap3i11kap1Ckap3Tk3i23) - 16._dp*(DyYvCkap3i11Tk3Ckap1kap1i21) -         & 
&  16._dp*(DyYvCkap3i11Tk3Ckap2kap1i22) - 16._dp*(DyYvCkap3i11Tk3Ckap3kap1i23) -         & 
&  16._dp*(DyYvCkap3i12kap2adjYvTvi23) - 16._dp*(DyYvCkap3i12kap2Ckap1Tk3i21)  
betaTv2 =  betaTv2- 16._dp*(DyYvCkap3i12kap2Ckap2Tk3i22) - 16._dp*(DyYvCkap3i12kap2Ckap3Tk3i23) -         & 
&  16._dp*(DyYvCkap3i12Tk3Ckap1kap2i21) - 16._dp*(DyYvCkap3i12Tk3Ckap2kap2i22) -         & 
&  16._dp*(DyYvCkap3i12Tk3Ckap3kap2i23) - 16._dp*(DyYvCkap3i13kap3adjYvTvi23) -          & 
&  16._dp*(DyYvCkap3i13kap3Ckap1Tk3i21) - 16._dp*(DyYvCkap3i13kap3Ckap2Tk3i22) -         & 
&  16._dp*(DyYvCkap3i13kap3Ckap3Tk3i23) - 16._dp*(DyYvCkap3i13Tk3Ckap1kap3i21) -         & 
&  16._dp*(DyYvCkap3i13Tk3Ckap2kap3i22) - 16._dp*(DyYvCkap3i13Tk3Ckap3kap3i23) -         & 
&  16._dp*(DyYvCkap3lami1Tk3Clami2) - 6._dp*(DyYvClami1TpTvCYvlami2) - 8._dp*(DyYvClami1TpYvCYvTlami2) +& 
&  (6*DyTvClami1lami2*g1p2)/5._dp + (12*DyYvClami1Tlami2*g1p2)/5._dp + 6*DyTvClami1lami2*g2p2 +& 
&  12*DyYvClami1Tlami2*g2p2 - (12*DyYvClami1lami2*g1p2*M1)/5._dp - 12*DyYvClami1lami2*g2p2*M2 -& 
&  14*DyYvClami1lami2*SPClamxxTlam - 10*DyTvClami1lami2*SPlamxxClam - 11*DyYvClami1Tlami2*SPlamxxClam -& 
&  4._dp*(TpTeCYeTpYeCYeYv) + (12*g1p2*TpTeCYeYv)/5._dp - 2*SPlamxxClam*TpTeCYeYv -      & 
&  4._dp*(TpYeCYeTpTeCYeYv) - 2._dp*(TpYeCYeTpYeCYeTv) + (6*g1p2*TpYeCYeTv)/5._dp -      & 
&  SPlamxxClam*TpYeCYeTv - (12*g1p2*M1*TpYeCYeYv)/5._dp - 2*SPClamxxTlam*TpYeCYeYv -     & 
&  12*DyYvClami1lami2*TradjYdTd - 6*TpYeCYeYv*TradjYdTd - 4*DyYvClami1lami2*TradjYeTe -  & 
&  2*TpYeCYeYv*TradjYeTe - 18*DyYvClami1lami2*TradjYuTu - 6*DyYvClami1lami2*TradjYvTv -  & 
&  8*DyTvi11kap1Ckap1i21*TrCkap1kap1 - 8*DyTvi12kap1Ckap2i21*TrCkap1kap1 -               & 
&  8*DyTvi13kap1Ckap3i21*TrCkap1kap1 - 8*DyTvi11kap1Ckap1i22*TrCkap1kap2 -               & 
&  8*DyTvi12kap1Ckap2i22*TrCkap1kap2 - 8*DyTvi13kap1Ckap3i22*TrCkap1kap2 -               & 
&  8*DyTvi11kap1Ckap1i23*TrCkap1kap3 - 8*DyTvi12kap1Ckap2i23*TrCkap1kap3 -               & 
&  8*DyTvi13kap1Ckap3i23*TrCkap1kap3 - 8*DyTvi11kap2Ckap1i21*TrCkap2kap1 -               & 
&  8*DyTvi12kap2Ckap2i21*TrCkap2kap1 - 8*DyTvi13kap2Ckap3i21*TrCkap2kap1  
betaTv2 =  betaTv2- 8*DyTvi11kap2Ckap1i22*TrCkap2kap2 - 8*DyTvi12kap2Ckap2i22*TrCkap2kap2 -               & 
&  8*DyTvi13kap2Ckap3i22*TrCkap2kap2 - 8*DyTvi11kap2Ckap1i23*TrCkap2kap3 -               & 
&  8*DyTvi12kap2Ckap2i23*TrCkap2kap3 - 8*DyTvi13kap2Ckap3i23*TrCkap2kap3 -               & 
&  8*DyTvi11kap3Ckap1i21*TrCkap3kap1 - 8*DyTvi12kap3Ckap2i21*TrCkap3kap1 -               & 
&  8*DyTvi13kap3Ckap3i21*TrCkap3kap1 - 8*DyTvi11kap3Ckap1i22*TrCkap3kap2 -               & 
&  8*DyTvi12kap3Ckap2i22*TrCkap3kap2 - 8*DyTvi13kap3Ckap3i22*TrCkap3kap2 -               & 
&  8*DyTvi11kap3Ckap1i23*TrCkap3kap3 - 8*DyTvi12kap3Ckap2i23*TrCkap3kap3 -               & 
&  8*DyTvi13kap3Ckap3i23*TrCkap3kap3 - 6*DyTvClami1lami2*TrYdadjYd - 12*DyYvClami1Tlami2*TrYdadjYd -& 
&  6*TpTeCYeYv*TrYdadjYd - 3*TpYeCYeTv*TrYdadjYd - 2*DyTvClami1lami2*TrYeadjYe -         & 
&  4*DyYvClami1Tlami2*TrYeadjYe - 2*TpTeCYeYv*TrYeadjYe - TpYeCYeTv*TrYeadjYe -          & 
&  12*DyTvClami1lami2*TrYuadjYu - 15*DyYvClami1Tlami2*TrYuadjYu - 4*DyTvClami1lami2*TrYvadjYv -& 
&  5*DyYvClami1Tlami2*TrYvadjYv - 2._dp*(TvadjYvTpYeCYeYv) + (6*g1p2*TvadjYvYv)/5._dp +  & 
&  6*g2p2*TvadjYvYv - 4*SPlamxxClam*TvadjYvYv - 12*TrYuadjYu*TvadjYvYv - 4*TrYvadjYv*TvadjYvYv -& 
&  6._dp*(TvadjYvYvadjYvYv) - (414*g1p4*M1*Yv)/25._dp - (18*g1p2*g2p2*M1*Yv)/5._dp -     & 
&  (18*g1p2*g2p2*M2*Yv)/5._dp - 30*g2p4*M2*Yv - 12*SPadjYvYvClamxxTlam*Yv -              & 
&  12*SPClamxxTpTvCYvlam*Yv - 12*SPClamxxTlam*SPlamxxClam*Yv - 6*SPlamxxClam*TradjYdTd*Yv -& 
&  2*SPlamxxClam*TradjYeTe*Yv - 2*TradjYeTeCYvTpYv*Yv + (8*g1p2*TradjYuTu*Yv)/5._dp +    & 
&  32*g3p2*TradjYuTu*Yv - 2*TradjYvTpYeCYeTv*Yv - 4*TradjYvTvadjkap1kap1*Yv -            & 
&  4*TradjYvTvadjkap2kap2*Yv - 4*TradjYvTvadjkap3kap3*Yv - 6*SPClamxxTlam*TrYdadjYd*Yv - & 
&  6*TrYdadjYuTuadjYd*Yv - 2*SPClamxxTlam*TrYeadjYe*Yv - 6*TrYuadjYdTdadjYu*Yv -         & 
&  (8*g1p2*M1*TrYuadjYu*Yv)/5._dp - 32*g3p2*M3*TrYuadjYu*Yv - 36*TrYuadjYuTuadjYu*Yv  
betaTv2 =  betaTv2- 12*TrYvadjYvTvadjYv*Yv - 4*TrYvCkap1TpTk1adjYv*Yv - 4*TrYvCkap2TpTk2adjYv*Yv -        & 
&  4*TrYvCkap3TpTk3adjYv*Yv - 4._dp*(YvadjYvTpTeCYeYv) - 4._dp*(YvadjYvTpYeCYeTv) +      & 
&  (12*g1p2*YvadjYvTv)/5._dp + 12*g2p2*YvadjYvTv - 5*SPlamxxClam*YvadjYvTv -             & 
&  15*TrYuadjYu*YvadjYvTv - 5*TrYvadjYv*YvadjYvTv - 8._dp*(YvadjYvTvadjYvYv) -           & 
&  (12*g1p2*M1*YvadjYvYv)/5._dp - 12*g2p2*M2*YvadjYvYv - 6*SPClamxxTlam*YvadjYvYv -      & 
&  18*TradjYuTu*YvadjYvYv - 6*TradjYvTv*YvadjYvYv - 6._dp*(YvadjYvYvadjYvTv) -           & 
&  4._dp*(YvCkap1TpTk1adjYvYv) - 16._dp*(YvCkap1TpYvCYvTpTk1) - 4._dp*(YvCkap2TpTk2adjYvYv) -& 
&  16._dp*(YvCkap2TpYvCYvTpTk2) - 4._dp*(YvCkap3TpTk3adjYvYv) - 16._dp*(YvCkap3TpYvCYvTpTk3) -& 
&  4*Yv*adjYvTvadjkap1kap2(2,1) - 4*Yv*adjYvTvadjkap1kap3(3,1) - 4*Yv*adjYvTvadjkap2kap1(1,& 
& 2) - 4*Yv*adjYvTvadjkap2kap3(3,2) - 4*Yv*adjYvTvadjkap3kap1(1,3) - 4*Yv*adjYvTvadjkap3kap2(2,& 
& 3) - 8*DyTvi11kap2Ckap1i21*Ckap1kap1(2,1) - 8*DyTvi12kap2Ckap2i21*Ckap1kap1(2,         & 
& 1) - 8*DyTvi13kap2Ckap3i21*Ckap1kap1(2,1) - 8*DyTvi11kap3Ckap1i21*Ckap1kap1(3,         & 
& 1) - 8*DyTvi12kap3Ckap2i21*Ckap1kap1(3,1) - 8*DyTvi13kap3Ckap3i21*Ckap1kap1(3,         & 
& 1) - 8*DyTvi11kap2Ckap1i22*Ckap1kap2(2,1) - 8*DyTvi12kap2Ckap2i22*Ckap1kap2(2,         & 
& 1) - 8*DyTvi13kap2Ckap3i22*Ckap1kap2(2,1) - 8*DyTvi11kap3Ckap1i22*Ckap1kap2(3,         & 
& 1) - 8*DyTvi12kap3Ckap2i22*Ckap1kap2(3,1) - 8*DyTvi13kap3Ckap3i22*Ckap1kap2(3,         & 
& 1) - 8*DyTvi11kap2Ckap1i23*Ckap1kap3(2,1) - 8*DyTvi12kap2Ckap2i23*Ckap1kap3(2,         & 
& 1) - 8*DyTvi13kap2Ckap3i23*Ckap1kap3(2,1) - 8*DyTvi11kap3Ckap1i23*Ckap1kap3(3,         & 
& 1) - 8*DyTvi12kap3Ckap2i23*Ckap1kap3(3,1) - 8*DyTvi13kap3Ckap3i23*Ckap1kap3(3,         & 
& 1) - 8*DyTvi11kap1Ckap1i21*Ckap2kap1(1,2) - 8*DyTvi12kap1Ckap2i21*Ckap2kap1(1,         & 
& 2) - 8*DyTvi13kap1Ckap3i21*Ckap2kap1(1,2) - 8*DyTvi11kap3Ckap1i21*Ckap2kap1(3,         & 
& 2) - 8*DyTvi12kap3Ckap2i21*Ckap2kap1(3,2) - 8*DyTvi13kap3Ckap3i21*Ckap2kap1(3,         & 
& 2) - 8*DyTvi11kap1Ckap1i22*Ckap2kap2(1,2) - 8*DyTvi12kap1Ckap2i22*Ckap2kap2(1,         & 
& 2) - 8*DyTvi13kap1Ckap3i22*Ckap2kap2(1,2) - 8*DyTvi11kap3Ckap1i22*Ckap2kap2(3,         & 
& 2) - 8*DyTvi12kap3Ckap2i22*Ckap2kap2(3,2) - 8*DyTvi13kap3Ckap3i22*Ckap2kap2(3,         & 
& 2) - 8*DyTvi11kap1Ckap1i23*Ckap2kap3(1,2) - 8*DyTvi12kap1Ckap2i23*Ckap2kap3(1,         & 
& 2) - 8*DyTvi13kap1Ckap3i23*Ckap2kap3(1,2) - 8*DyTvi11kap3Ckap1i23*Ckap2kap3(3,         & 
& 2) - 8*DyTvi12kap3Ckap2i23*Ckap2kap3(3,2) - 8*DyTvi13kap3Ckap3i23*Ckap2kap3(3,         & 
& 2) - 8*DyTvi11kap1Ckap1i21*Ckap3kap1(1,3) - 8*DyTvi12kap1Ckap2i21*Ckap3kap1(1,         & 
& 3) - 8*DyTvi13kap1Ckap3i21*Ckap3kap1(1,3) - 8*DyTvi11kap2Ckap1i21*Ckap3kap1(2,         & 
& 3) - 8*DyTvi12kap2Ckap2i21*Ckap3kap1(2,3) - 8*DyTvi13kap2Ckap3i21*Ckap3kap1(2,         & 
& 3) - 8*DyTvi11kap1Ckap1i22*Ckap3kap2(1,3) - 8*DyTvi12kap1Ckap2i22*Ckap3kap2(1,         & 
& 3) - 8*DyTvi13kap1Ckap3i22*Ckap3kap2(1,3) - 8*DyTvi11kap2Ckap1i22*Ckap3kap2(2,         & 
& 3) - 8*DyTvi12kap2Ckap2i22*Ckap3kap2(2,3) - 8*DyTvi13kap2Ckap3i22*Ckap3kap2(2,         & 
& 3) - 8*DyTvi11kap1Ckap1i23*Ckap3kap3(1,3) - 8*DyTvi12kap1Ckap2i23*Ckap3kap3(1,         & 
& 3) - 8*DyTvi13kap1Ckap3i23*Ckap3kap3(1,3) - 8*DyTvi11kap2Ckap1i23*Ckap3kap3(2,         & 
& 3) - 8*DyTvi12kap2Ckap2i23*Ckap3kap3(2,3) - 8*DyTvi13kap2Ckap3i23*Ckap3kap3(2,         & 
& 3) - 4*DyTpkap1Ckap1TpTv1i1lami2*Conjg(lam(1)) - 4*DyTpkap1Ckap2TpTv2i1lami2*Conjg(lam(1))  
betaTv2 =  betaTv2- 4*DyTpkap1Ckap3TpTv3i1lami2*Conjg(lam(1)) - 8*DyTvi11kap1Ckap1lami2*Conjg(lam(1)) -   & 
&  8*DyTvi12kap1Ckap2lami2*Conjg(lam(1)) - 8*DyTvi13kap1Ckap3lami2*Conjg(lam(1)) -       & 
&  2*DyYvadjkap1kap1i11Tlami2*Conjg(lam(1)) - 4*DyYvadjkap1Tk1i11lami2*Conjg(lam(1)) -   & 
&  2*DyYvadjkap2kap1i12Tlami2*Conjg(lam(1)) - 4*DyYvadjkap2Tk1i12lami2*Conjg(lam(1)) -   & 
&  2*DyYvadjkap3kap1i13Tlami2*Conjg(lam(1)) - 4*DyYvadjkap3Tk1i13lami2*Conjg(lam(1)) -   & 
&  4*DyTpkap2Ckap1TpTv1i1lami2*Conjg(lam(2)) - 4*DyTpkap2Ckap2TpTv2i1lami2*Conjg(lam(2)) -& 
&  4*DyTpkap2Ckap3TpTv3i1lami2*Conjg(lam(2)) - 8*DyTvi11kap2Ckap1lami2*Conjg(lam(2)) -   & 
&  8*DyTvi12kap2Ckap2lami2*Conjg(lam(2)) - 8*DyTvi13kap2Ckap3lami2*Conjg(lam(2)) -       & 
&  2*DyYvadjkap1kap2i11Tlami2*Conjg(lam(2)) - 4*DyYvadjkap1Tk2i11lami2*Conjg(lam(2)) -   & 
&  2*DyYvadjkap2kap2i12Tlami2*Conjg(lam(2)) - 4*DyYvadjkap2Tk2i12lami2*Conjg(lam(2)) -   & 
&  2*DyYvadjkap3kap2i13Tlami2*Conjg(lam(2)) - 4*DyYvadjkap3Tk2i13lami2*Conjg(lam(2)) -   & 
&  4*DyTpkap3Ckap1TpTv1i1lami2*Conjg(lam(3)) - 4*DyTpkap3Ckap2TpTv2i1lami2*Conjg(lam(3)) -& 
&  4*DyTpkap3Ckap3TpTv3i1lami2*Conjg(lam(3)) - 8*DyTvi11kap3Ckap1lami2*Conjg(lam(3)) -   & 
&  8*DyTvi12kap3Ckap2lami2*Conjg(lam(3)) - 8*DyTvi13kap3Ckap3lami2*Conjg(lam(3)) -       & 
&  2*DyYvadjkap1kap3i11Tlami2*Conjg(lam(3)) - 4*DyYvadjkap1Tk3i11lami2*Conjg(lam(3)) -   & 
&  2*DyYvadjkap2kap3i12Tlami2*Conjg(lam(3)) - 4*DyYvadjkap2Tk3i12lami2*Conjg(lam(3)) -   & 
&  2*DyYvadjkap3kap3i13Tlami2*Conjg(lam(3)) - 4*DyYvadjkap3Tk3i13lami2*Conjg(lam(3)) +   & 
&  (207*g1p4*Tv)/50._dp + (9*g1p2*g2p2*Tv)/5._dp + (15*g2p4*Tv)/2._dp - 6*SPlamxxadjYvYvClam*Tv -& 
&  3*SPlamxxClamp2*Tv - 3*SPlamxxClam*TrYdadjYd*Tv - 3*TrYdadjYuYuadjYd*Tv -             & 
&  SPlamxxClam*TrYeadjYe*Tv + (4*g1p2*TrYuadjYu*Tv)/5._dp + 16*g3p2*TrYuadjYu*Tv -       & 
&  9*TrYuadjYuYuadjYu*Tv - TrYvadjYvTpYeCYe*Tv - 3*TrYvadjYvYvadjYv*Tv - 2*TrYvCkap1kap1adjYv*Tv  
betaTv2 =  betaTv2- 2*TrYvCkap2kap2adjYv*Tv - 2*TrYvCkap3kap3adjYv*Tv - 2*adjYvYvCkap1kap2(2,             & 
& 1)*Tv - 2*adjYvYvCkap1kap3(3,1)*Tv - 2*adjYvYvCkap2kap1(1,2)*Tv - 2*adjYvYvCkap2kap3(3,& 
& 2)*Tv - 2*adjYvYvCkap3kap1(1,3)*Tv - 2*adjYvYvCkap3kap2(2,3)*Tv - 2*Conjg(lam(1))*Tv*Tpkap1Ckap1lam(1) -& 
&  4*Yv*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) - 2*Conjg(lam(1))*Tv*Tpkap1Ckap2lam(2) -        & 
&  4*Yv*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) - 2*Conjg(lam(1))*Tv*Tpkap1Ckap3lam(3) -        & 
&  4*Yv*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) - 2*Conjg(lam(2))*Tv*Tpkap2Ckap1lam(1) -        & 
&  4*Yv*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) - 2*Conjg(lam(2))*Tv*Tpkap2Ckap2lam(2) -        & 
&  4*Yv*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) - 2*Conjg(lam(2))*Tv*Tpkap2Ckap3lam(3) -        & 
&  4*Yv*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) - 2*Conjg(lam(3))*Tv*Tpkap3Ckap1lam(1) -        & 
&  4*Yv*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) - 2*Conjg(lam(3))*Tv*Tpkap3Ckap2lam(2) -        & 
&  4*Yv*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) - 2*Conjg(lam(3))*Tv*Tpkap3Ckap3lam(3) -        & 
&  4*Yv*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) - 4*Yv*Conjg(lam(1))*TpTk1Ckap1lam(1) -         & 
&  4*Yv*Conjg(lam(1))*TpTk1Ckap2lam(2) - 4*Yv*Conjg(lam(1))*TpTk1Ckap3lam(3) -           & 
&  4*Yv*Conjg(lam(2))*TpTk2Ckap1lam(1) - 4*Yv*Conjg(lam(2))*TpTk2Ckap2lam(2) -           & 
&  4*Yv*Conjg(lam(2))*TpTk2Ckap3lam(3) - 4*Yv*Conjg(lam(3))*TpTk3Ckap1lam(1) -           & 
&  4*Yv*Conjg(lam(3))*TpTk3Ckap2lam(2) - 4*Yv*Conjg(lam(3))*TpTk3Ckap3lam(3) -           & 
&  16*DyYvCkap1i11kap1Clami2*Tlam(1) - 16*DyYvCkap1i12kap2Clami2*Tlam(1) -               & 
&  16*DyYvCkap1i13kap3Clami2*Tlam(1) - 16*DyYvCkap2i11kap1Clami2*Tlam(2) -               & 
&  16*DyYvCkap2i12kap2Clami2*Tlam(2) - 16*DyYvCkap2i13kap3Clami2*Tlam(2) -               & 
&  16*DyYvCkap3i11kap1Clami2*Tlam(3) - 16*DyYvCkap3i12kap2Clami2*Tlam(3) -               & 
&  16*DyYvCkap3i13kap3Clami2*Tlam(3)

 
DTv = oo16pi2*( betaTv1 + oo16pi2 * betaTv2 ) 

 
Else 
DTv = oo16pi2* betaTv1 
End If 
 
 
Call Chop(DTv) 

!-------------------- 
! Tu 
!-------------------- 
 
betaTu1  = TuadjYdYd + 5._dp*(TuadjYuYu) + (26*g1p2*M1*Yu)/15._dp + (32*g3p2*M3*Yu)   & 
& /3._dp + 6*g2p2*M2*Yu + 2*SPClamxxTlam*Yu + 6*TradjYuTu*Yu + 2*TradjYvTv*Yu +          & 
&  2._dp*(YuadjYdTd) + 4._dp*(YuadjYuTu) - (13*g1p2*Tu)/15._dp - 3*g2p2*Tu -             & 
&  (16*g3p2*Tu)/3._dp + SPlamxxClam*Tu + 3*TrYuadjYu*Tu + TrYvadjYv*Tu

 
 
If (TwoLoopRGE) Then 
betaTu2 = ((2._dp*(g1p2) - 5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe))*TuadjYdYd)/5._dp -  & 
&  2._dp*(TuadjYdYdadjYdYd) - 4._dp*(TuadjYdYdadjYuYu) + 12*g2p2*TuadjYuYu -             & 
&  5*SPlamxxClam*TuadjYuYu - 15*TrYuadjYu*TuadjYuYu - 5*TrYvadjYv*TuadjYuYu -            & 
&  6._dp*(TuadjYuYuadjYuYu) - (5486*g1p4*M1*Yu)/225._dp - 2*g1p2*g2p2*M1*Yu -            & 
&  (272*g1p2*g3p2*M1*Yu)/45._dp - (272*g1p2*g3p2*M3*Yu)/45._dp - 16*g2p2*g3p2*M3*Yu +    & 
&  (64*g3p4*M3*Yu)/9._dp - 2*g1p2*g2p2*M2*Yu - 30*g2p4*M2*Yu - 16*g2p2*g3p2*M2*Yu -      & 
&  12*SPadjYvYvClamxxTlam*Yu - 12*SPClamxxTpTvCYvlam*Yu - 12*SPClamxxTlam*SPlamxxClam*Yu -& 
&  6*SPlamxxClam*TradjYdTd*Yu - 2*SPlamxxClam*TradjYeTe*Yu - 2*TradjYeTeCYvTpYv*Yu +     & 
&  (8*g1p2*TradjYuTu*Yu)/5._dp + 32*g3p2*TradjYuTu*Yu - 2*TradjYvTpYeCYeTv*Yu -          & 
&  4*TradjYvTvadjkap1kap1*Yu - 4*TradjYvTvadjkap2kap2*Yu - 4*TradjYvTvadjkap3kap3*Yu -   & 
&  6*SPClamxxTlam*TrYdadjYd*Yu - 6*TrYdadjYuTuadjYd*Yu - 2*SPClamxxTlam*TrYeadjYe*Yu -   & 
&  6*TrYuadjYdTdadjYu*Yu - (8*g1p2*M1*TrYuadjYu*Yu)/5._dp - 32*g3p2*M3*TrYuadjYu*Yu -    & 
&  36*TrYuadjYuTuadjYu*Yu - 12*TrYvadjYvTvadjYv*Yu - 4*TrYvCkap1TpTk1adjYv*Yu -          & 
&  4*TrYvCkap2TpTk2adjYv*Yu - 4*TrYvCkap3TpTk3adjYv*Yu + (4*g1p2*YuadjYdTd)/5._dp -      & 
&  2*SPlamxxClam*YuadjYdTd - 6*TrYdadjYd*YuadjYdTd - 2*TrYeadjYe*YuadjYdTd -             & 
&  4._dp*(YuadjYdTdadjYdYd) - 4._dp*(YuadjYdTdadjYuYu) - (4*g1p2*M1*YuadjYdYd)/5._dp -   & 
&  2*SPClamxxTlam*YuadjYdYd - 6*TradjYdTd*YuadjYdYd - 2*TradjYeTe*YuadjYdYd -            & 
&  4._dp*(YuadjYdYdadjYdTd) - 2._dp*(YuadjYdYdadjYuTu) + (6*g1p2*YuadjYuTu)/5._dp +      & 
&  6*g2p2*YuadjYuTu - 4*SPlamxxClam*YuadjYuTu - 12*TrYuadjYu*YuadjYuTu - 4*TrYvadjYv*YuadjYuTu -& 
&  8._dp*(YuadjYuTuadjYuYu) - (4*g1p2*M1*YuadjYuYu)/5._dp - 12*g2p2*M2*YuadjYuYu -       & 
&  6*SPClamxxTlam*YuadjYuYu - 18*TradjYuTu*YuadjYuYu - 6*TradjYvTv*YuadjYuYu  
betaTu2 =  betaTu2- 6._dp*(YuadjYuYuadjYuTu) - 4*Yu*adjYvTvadjkap1kap2(2,1) - 4*Yu*adjYvTvadjkap1kap3(3,  & 
& 1) - 4*Yu*adjYvTvadjkap2kap1(1,2) - 4*Yu*adjYvTvadjkap2kap3(3,2) - 4*Yu*adjYvTvadjkap3kap1(1,& 
& 3) - 4*Yu*adjYvTvadjkap3kap2(2,3) + (2743*g1p4*Tu)/450._dp + g1p2*g2p2*Tu +            & 
&  (15*g2p4*Tu)/2._dp + (136*g1p2*g3p2*Tu)/45._dp + 8*g2p2*g3p2*Tu - (16*g3p4*Tu)/9._dp -& 
&  6*SPlamxxadjYvYvClam*Tu - 3*SPlamxxClamp2*Tu - 3*SPlamxxClam*TrYdadjYd*Tu -           & 
&  3*TrYdadjYuYuadjYd*Tu - SPlamxxClam*TrYeadjYe*Tu + (4*g1p2*TrYuadjYu*Tu)/5._dp +      & 
&  16*g3p2*TrYuadjYu*Tu - 9*TrYuadjYuYuadjYu*Tu - TrYvadjYvTpYeCYe*Tu - 3*TrYvadjYvYvadjYv*Tu -& 
&  2*TrYvCkap1kap1adjYv*Tu - 2*TrYvCkap2kap2adjYv*Tu - 2*TrYvCkap3kap3adjYv*Tu -         & 
&  2*adjYvYvCkap1kap2(2,1)*Tu - 2*adjYvYvCkap1kap3(3,1)*Tu - 2*adjYvYvCkap2kap1(1,       & 
& 2)*Tu - 2*adjYvYvCkap2kap3(3,2)*Tu - 2*adjYvYvCkap3kap1(1,3)*Tu - 2*adjYvYvCkap3kap2(2,& 
& 3)*Tu - 2*Conjg(lam(1))*Tu*Tpkap1Ckap1lam(1) - 4*Yu*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) - & 
&  2*Conjg(lam(1))*Tu*Tpkap1Ckap2lam(2) - 4*Yu*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) -        & 
&  2*Conjg(lam(1))*Tu*Tpkap1Ckap3lam(3) - 4*Yu*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) -        & 
&  2*Conjg(lam(2))*Tu*Tpkap2Ckap1lam(1) - 4*Yu*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) -        & 
&  2*Conjg(lam(2))*Tu*Tpkap2Ckap2lam(2) - 4*Yu*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) -        & 
&  2*Conjg(lam(2))*Tu*Tpkap2Ckap3lam(3) - 4*Yu*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) -        & 
&  2*Conjg(lam(3))*Tu*Tpkap3Ckap1lam(1) - 4*Yu*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) -        & 
&  2*Conjg(lam(3))*Tu*Tpkap3Ckap2lam(2) - 4*Yu*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) -        & 
&  2*Conjg(lam(3))*Tu*Tpkap3Ckap3lam(3) - 4*Yu*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) -        & 
&  4*Yu*Conjg(lam(1))*TpTk1Ckap1lam(1) - 4*Yu*Conjg(lam(1))*TpTk1Ckap2lam(2) -           & 
&  4*Yu*Conjg(lam(1))*TpTk1Ckap3lam(3) - 4*Yu*Conjg(lam(2))*TpTk2Ckap1lam(1)  
betaTu2 =  betaTu2- 4*Yu*Conjg(lam(2))*TpTk2Ckap2lam(2) - 4*Yu*Conjg(lam(2))*TpTk2Ckap3lam(3) -           & 
&  4*Yu*Conjg(lam(3))*TpTk3Ckap1lam(1) - 4*Yu*Conjg(lam(3))*TpTk3Ckap2lam(2) -           & 
&  4*Yu*Conjg(lam(3))*TpTk3Ckap3lam(3)

 
DTu = oo16pi2*( betaTu1 + oo16pi2 * betaTu2 ) 

 
Else 
DTu = oo16pi2* betaTu1 
End If 
 
 
Call Chop(DTu) 

!-------------------- 
! Tk 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaTk1(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk1(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk1(i1,i2,3) = 0
End Do
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaTk2(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk2(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk2(i1,i2,3) = 0
End Do
End Do

 
DTk = oo16pi2*( betaTk1 + oo16pi2 * betaTk2 ) 

 
Else 
DTk = oo16pi2* betaTk1 
End If 
 
 
!-------------------- 
! mq2 
!-------------------- 
 
betamq21  = 2._dp*(adjTdTd) + 2._dp*(adjTuTu) + 2._dp*(adjYdmd2Yd) + adjYdYdmq2 +     & 
&  2._dp*(adjYumu2Yu) + adjYuYumq2 - (2*AbsM1*g1p2*id3R)/15._dp - 6*AbsM2*g2p2*id3R -    & 
&  (32*AbsM3*g3p2*id3R)/3._dp + 2*adjYdYd*mHd2 + 2*adjYuYu*mHu2 + mq2adjYdYd +           & 
&  mq2adjYuYu + g1*id3R*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamq22 = -4._dp*(adjTdTdadjYdYd) - 4._dp*(adjTdYdadjYdTd) - 4._dp*(adjTuTuadjYuYu) -           & 
&  4._dp*(adjTuYuadjYuTu) - 4._dp*(adjYdmd2YdadjYdYd) - 4._dp*(adjYdTdadjTdYd) -         & 
&  4._dp*(adjYdYdadjTdTd) - 4._dp*(adjYdYdadjYdmd2Yd) - 2._dp*(adjYdYdadjYdYdmq2) -      & 
&  4._dp*(adjYdYdmq2adjYdYd) - 4._dp*(adjYumu2YuadjYuYu) - 4._dp*(adjYuTuadjTuYu) -      & 
&  4._dp*(adjYuYuadjTuTu) - 4._dp*(adjYuYuadjYumu2Yu) - 2._dp*(adjYuYuadjYuYumq2) -      & 
&  4._dp*(adjYuYumq2adjYuYu) + (8*adjTuTu*g1p2)/5._dp + (4*adjYdmd2Yd*g1p2)/5._dp +      & 
&  (8*AbsM1*adjYdYd*g1p2)/5._dp + (2*adjYdYdmq2*g1p2)/5._dp + (8*adjYumu2Yu*g1p2)/5._dp +& 
&  (16*AbsM1*adjYuYu*g1p2)/5._dp + (4*adjYuYumq2*g1p2)/5._dp + (199*AbsM1*g1p4*id3R)/75._dp +& 
&  (2*AbsM1*g1p2*g2p2*id3R)/5._dp + (2*AbsM2*g1p2*g2p2*id3R)/5._dp + 33*AbsM2*g2p4*id3R +& 
&  (32*AbsM1*g1p2*g3p2*id3R)/45._dp + (32*AbsM3*g1p2*g3p2*id3R)/45._dp + 32*AbsM2*g2p2*g3p2*id3R +& 
&  32*AbsM3*g2p2*g3p2*id3R - (128*AbsM3*g3p4*id3R)/3._dp - (4*adjTdYd*g1p2*M1)/5._dp -   & 
&  (8*adjTuYu*g1p2*M1)/5._dp - 8*adjYdYdadjYdYd*mHd2 + (4*adjYdYd*g1p2*mHd2)/5._dp -     & 
&  8*adjYuYuadjYuYu*mHu2 + (8*adjYuYu*g1p2*mHu2)/5._dp + (2*g1p2*mq2adjYdYd)/5._dp -     & 
&  2._dp*(mq2adjYdYdadjYdYd) + (4*g1p2*mq2adjYuYu)/5._dp - 2._dp*(mq2adjYuYuadjYuYu) -   & 
&  2*adjTdYd*SPClamxxTlam - 2*adjTuYu*SPClamxxTlam + 2*adjYdYd*SPlamxxadjYvmlHd2 +       & 
&  2*adjYuYu*SPlamxxadjYvmlHd2 - 2*adjTuTu*SPlamxxClam - 2*adjYdmd2Yd*SPlamxxClam -      & 
&  adjYdYdmq2*SPlamxxClam - 2*adjYumu2Yu*SPlamxxClam - adjYuYumq2*SPlamxxClam -          & 
&  4*adjYdYd*mHd2*SPlamxxClam - 2*adjYuYu*mHd2*SPlamxxClam - 2*adjYdYd*mHu2*SPlamxxClam -& 
&  4*adjYuYu*mHu2*SPlamxxClam - mq2adjYdYd*SPlamxxClam - mq2adjYuYu*SPlamxxClam -        & 
&  2*adjYdYd*SPlamxxCmv2Clam - 2*adjYuYu*SPlamxxCmv2Clam - 2*adjYdTd*SPlamxxCTlam -      & 
&  2*adjYuTu*SPlamxxCTlam - 2*adjYdYd*SPTlamxxCTlam - 2*adjYuYu*SPTlamxxCTlam  
betamq22 =  betamq22- 6*adjTdYd*TradjYdTd - 2*adjTdYd*TradjYeTe - 6*adjTuYu*TradjYuTu - 2*adjTuYu*TradjYvTv -& 
&  6*adjYdYd*TrCTdTpTd - 6*adjYdTd*TrCTdTpYd - 2*adjYdYd*TrCTeTpTe - 2*adjYdTd*TrCTeTpYe -& 
&  6*adjYuYu*TrCTuTpTu - 6*adjYuTu*TrCTuTpYu - 2*adjYuYu*TrCTvTpTv - 2*adjYuTu*TrCTvTpYv -& 
&  6*adjYdYd*Trmd2YdadjYd - 2*adjYdYd*Trme2YeadjYe - 2*adjYdYd*Trml2adjYeYe -            & 
&  6*adjYdYd*Trmq2adjYdYd - 6*adjYuYu*Trmq2adjYuYu - 6*adjYuYu*Trmu2YuadjYu -            & 
&  6*adjYdmd2Yd*TrYdadjYd - 3*adjYdYdmq2*TrYdadjYd - 12*adjYdYd*mHd2*TrYdadjYd -         & 
&  3*mq2adjYdYd*TrYdadjYd - 2*adjYdmd2Yd*TrYeadjYe - adjYdYdmq2*TrYeadjYe -              & 
&  4*adjYdYd*mHd2*TrYeadjYe - mq2adjYdYd*TrYeadjYe + (2*adjTdTd*(2._dp*(g1p2) -          & 
&  5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)))/5._dp - 6*adjTuTu*TrYuadjYu -       & 
&  6*adjYumu2Yu*TrYuadjYu - 3*adjYuYumq2*TrYuadjYu - 12*adjYuYu*mHu2*TrYuadjYu -         & 
&  3*mq2adjYuYu*TrYuadjYu - 2*adjTuTu*TrYvadjYv - 2*adjYumu2Yu*TrYvadjYv -               & 
&  adjYuYumq2*TrYvadjYv - 4*adjYuYu*mHu2*TrYvadjYv - mq2adjYuYu*TrYvadjYv -              & 
&  2*adjYuYu*TrYvadjYvCml2 - 2*adjYuYu*TrYvCmv2adjYv - (4*adjYdTd*g1p2*Conjg(M1))/5._dp -& 
&  (8*adjYuTu*g1p2*Conjg(M1))/5._dp + (16*g1p2*g3p2*id3R*M3*Conjg(M1))/45._dp +          & 
&  (g1p2*g2p2*id3R*M2*Conjg(M1))/5._dp + (16*g1p2*g3p2*id3R*M1*Conjg(M3))/45._dp +       & 
&  16*g2p2*g3p2*id3R*M2*Conjg(M3) + (g1p2*g2p2*id3R*M1*Conjg(M2))/5._dp + 16*g2p2*g3p2*id3R*M3*Conjg(M2) +& 
&  6*g2p4*id3R*Tr2(2) + (32*g3p4*id3R*Tr2(3))/3._dp + (2*g1p2*id3R*Tr2U1(1,              & 
& 1))/15._dp + 4*g1*id3R*ooSqrt15*Tr3(1)

 
Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 ) 

 
Else 
Dmq2 = oo16pi2* betamq21 
End If 
 
 
Call Chop(Dmq2) 

Forall(i1=1:3) Dmq2(i1,i1) =  Real(Dmq2(i1,i1),dp) 
Dmq2 = 0.5_dp*(Dmq2+ Conjg(Transpose(Dmq2)) ) 
!-------------------- 
! ml2 
!-------------------- 
 
betaml21  = 2._dp*(adjTeTe) + 2._dp*(adjYeme2Ye) + adjYeYeml2 + 2._dp*(CTvTpTv)       & 
&  + 2._dp*(CYvmv2TpYv) + CYvTpYvml2 - DyCYvlami1mlHd2i2 - (6*AbsM1*g1p2*id3R)           & 
& /5._dp - 6*AbsM2*g2p2*id3R + 2*adjYeYe*mHd2 + 2*CYvTpYv*mHu2 + ml2adjYeYe +            & 
&  ml2CYvTpYv - g1*id3R*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betaml22 = -4._dp*(adjTeTeadjYeYe) - 4._dp*(adjTeYeadjYeTe) - 4._dp*(adjYeme2YeadjYeYe) -        & 
&  4._dp*(adjYeTeadjTeYe) - 4._dp*(adjYeYeadjTeTe) - 4._dp*(adjYeYeadjYeme2Ye) -         & 
&  2._dp*(adjYeYeadjYeYeml2) - 4._dp*(adjYeYeml2adjYeYe) - 4._dp*(CTvTpTvCYvTpYv) -      & 
&  4._dp*(CTvTpYvCYvTpTv) - 2._dp*(CYvkap1adjkap1TpYvml2) - 2._dp*(CYvkap2adjkap2TpYvml2) -& 
&  2._dp*(CYvkap3adjkap3TpYvml2) - 4._dp*(CYvmv2kap1adjkap1TpYv) - 4._dp*(CYvmv2kap2adjkap2TpYv) -& 
&  4._dp*(CYvmv2kap3adjkap3TpYv) - 4._dp*(CYvmv2TpYvCYvTpYv) - 4._dp*(CYvTk1adjTk1TpYv) -& 
&  4._dp*(CYvTk2adjTk2TpYv) - 4._dp*(CYvTk3adjTk3TpYv) - 4._dp*(CYvTpTvCTvTpYv) -        & 
&  4._dp*(CYvTpYvCTvTpTv) - 4._dp*(CYvTpYvCYvmv2TpYv) - 2._dp*(CYvTpYvCYvTpYvml2) -      & 
&  4._dp*(CYvTpYvml2CYvTpYv) - DyadjYeYeCYvlami1mlHd2i2 - 4._dp*(DyCTvi11YvCkap1Tk1i21) -& 
&  4._dp*(DyCTvi11YvCkap2Tk1i22) - 4._dp*(DyCTvi11YvCkap3Tk1i23) - 4._dp*(DyCTvi12YvCkap1Tk2i21) -& 
&  4._dp*(DyCTvi12YvCkap2Tk2i22) - 4._dp*(DyCTvi12YvCkap3Tk2i23) - 4._dp*(DyCTvi13YvCkap1Tk3i21) -& 
&  4._dp*(DyCTvi13YvCkap2Tk3i22) - 4._dp*(DyCTvi13YvCkap3Tk3i23) - 4._dp*(DyCTvlami1TvClami2) -& 
&  4._dp*(DyCTvTlami1YvClami2) - 8._dp*(DyCYvi11Yvadjkap1mv2kap1i21) - 8._dp*(DyCYvi11Yvadjkap2mv2kap1i22) -& 
&  8._dp*(DyCYvi11Yvadjkap3mv2kap1i23) - 4._dp*(DyCYvi11YvCmv2Ckap1kap1i21) -            & 
&  4._dp*(DyCYvi11YvCmv2Ckap2kap1i22) - 4._dp*(DyCYvi11YvCmv2Ckap3kap1i23) -             & 
&  8._dp*(DyCYvi12Yvadjkap1mv2kap2i21) - 8._dp*(DyCYvi12Yvadjkap2mv2kap2i22) -           & 
&  8._dp*(DyCYvi12Yvadjkap3mv2kap2i23) - 4._dp*(DyCYvi12YvCmv2Ckap1kap2i21) -            & 
&  4._dp*(DyCYvi12YvCmv2Ckap2kap2i22) - 4._dp*(DyCYvi12YvCmv2Ckap3kap2i23) -             & 
&  8._dp*(DyCYvi13Yvadjkap1mv2kap3i21) - 8._dp*(DyCYvi13Yvadjkap2mv2kap3i22) -           & 
&  8._dp*(DyCYvi13Yvadjkap3mv2kap3i23) - 4._dp*(DyCYvi13YvCmv2Ckap1kap3i21) -            & 
&  4._dp*(DyCYvi13YvCmv2Ckap2kap3i22) - 4._dp*(DyCYvi13YvCmv2Ckap3kap3i23)  
betaml22 =  betaml22- 2._dp*(DyCYvlami1Cml2YvClami2) - 4._dp*(DyCYvlami1TvCTlami2) + 4._dp*(DyCYvlami1YvadjYvmlHd2i2) -& 
&  4._dp*(DyCYvlami1YvCmv2Clami2) - 4._dp*(DyCYvmv2lami1YvClami2) - 4._dp*(DyCYvTlami1YvCTlami2) +& 
&  2._dp*(DyCYvTpYvCYvlami1mlHd2i2) - 2._dp*(Dyml2CYvi11Yvadjkap1kap1i21) -              & 
&  2._dp*(Dyml2CYvi11Yvadjkap2kap1i22) - 2._dp*(Dyml2CYvi11Yvadjkap3kap1i23) -           & 
&  2._dp*(Dyml2CYvi12Yvadjkap1kap2i21) - 2._dp*(Dyml2CYvi12Yvadjkap2kap2i22) -           & 
&  2._dp*(Dyml2CYvi12Yvadjkap3kap2i23) - 2._dp*(Dyml2CYvi13Yvadjkap1kap3i21) -           & 
&  2._dp*(Dyml2CYvi13Yvadjkap2kap3i22) - 2._dp*(Dyml2CYvi13Yvadjkap3kap3i23) -           & 
&  2._dp*(Dyml2CYvlami1YvClami2) + (12*adjYeme2Ye*g1p2)/5._dp + (24*AbsM1*adjYeYe*g1p2)/5._dp +& 
&  (6*adjYeYeml2*g1p2)/5._dp + (621*AbsM1*g1p4*id3R)/25._dp + (18*AbsM1*g1p2*g2p2*id3R)/5._dp +& 
&  (18*AbsM2*g1p2*g2p2*id3R)/5._dp + 33*AbsM2*g2p4*id3R - (12*adjTeYe*g1p2*M1)/5._dp -   & 
&  8*adjYeYeadjYeYe*mHd2 - 4*DyCYvlami1YvClami2*mHd2 + (12*adjYeYe*g1p2*mHd2)/5._dp -    & 
&  8*CYvTpYvCYvTpYv*mHu2 - 4*DyCYvi11Yvadjkap1kap1i21*mHu2 - 4*DyCYvi11Yvadjkap2kap1i22*mHu2 -& 
&  4*DyCYvi11Yvadjkap3kap1i23*mHu2 - 4*DyCYvi12Yvadjkap1kap2i21*mHu2 - 4*DyCYvi12Yvadjkap2kap2i22*mHu2 -& 
&  4*DyCYvi12Yvadjkap3kap2i23*mHu2 - 4*DyCYvi13Yvadjkap1kap3i21*mHu2 - 4*DyCYvi13Yvadjkap2kap3i22*mHu2 -& 
&  4*DyCYvi13Yvadjkap3kap3i23*mHu2 - 8*DyCYvlami1YvClami2*mHu2 + (6*g1p2*ml2adjYeYe)/5._dp -& 
&  2._dp*(ml2adjYeYeadjYeYe) - 2._dp*(ml2CYvTpYvCYvTpYv) - 2*adjTeYe*SPClamxxTlam -      & 
&  2*CTvTpYv*SPClamxxTlam + 2*adjYeYe*SPlamxxadjYvmlHd2 + 2*CYvTpYv*SPlamxxadjYvmlHd2 -  & 
&  2*adjYeme2Ye*SPlamxxClam - adjYeYeml2*SPlamxxClam - 2*CTvTpTv*SPlamxxClam -           & 
&  2*CYvmv2TpYv*SPlamxxClam - CYvTpYvml2*SPlamxxClam + 3*DyCYvlami1mlHd2i2*SPlamxxClam - & 
&  4*adjYeYe*mHd2*SPlamxxClam - 2*CYvTpYv*mHd2*SPlamxxClam - 2*adjYeYe*mHu2*SPlamxxClam -& 
&  4*CYvTpYv*mHu2*SPlamxxClam - ml2adjYeYe*SPlamxxClam - ml2CYvTpYv*SPlamxxClam  
betaml22 =  betaml22- 2*adjYeYe*SPlamxxCmv2Clam - 2*CYvTpYv*SPlamxxCmv2Clam - 2*adjYeTe*SPlamxxCTlam -      & 
&  2*CYvTpTv*SPlamxxCTlam - 2*adjYeYe*SPTlamxxCTlam - 2*CYvTpYv*SPTlamxxCTlam -          & 
&  6*adjTeYe*TradjYdTd - 2*adjTeYe*TradjYeTe - 6*CTvTpYv*TradjYuTu - 2*CTvTpYv*TradjYvTv -& 
&  4*DyCTvi11Tvi21*TrCkap1Tpkap1 - 4*DyCTvi12Tvi21*TrCkap1Tpkap2 - 4*DyCTvi13Tvi21*TrCkap1Tpkap3 -& 
&  4*DyCTvi11Tvi22*TrCkap2Tpkap1 - 4*DyCTvi12Tvi22*TrCkap2Tpkap2 - 4*DyCTvi13Tvi22*TrCkap2Tpkap3 -& 
&  4*DyCTvi11Tvi23*TrCkap3Tpkap1 - 4*DyCTvi12Tvi23*TrCkap3Tpkap2 - 4*DyCTvi13Tvi23*TrCkap3Tpkap3 -& 
&  6*adjYeYe*TrCTdTpTd - 6*adjYeTe*TrCTdTpYd - 2*adjYeYe*TrCTeTpTe - 2*adjYeTe*TrCTeTpYe -& 
&  4*DyCYvi11Tvi21*TrCTk1Tpkap1 - 4*DyCYvi12Tvi21*TrCTk1Tpkap2 - 4*DyCYvi13Tvi21*TrCTk1Tpkap3 -& 
&  4*DyCYvi11Tvi22*TrCTk2Tpkap1 - 4*DyCYvi12Tvi22*TrCTk2Tpkap2 - 4*DyCYvi13Tvi22*TrCTk2Tpkap3 -& 
&  4*DyCYvi11Tvi23*TrCTk3Tpkap1 - 4*DyCYvi12Tvi23*TrCTk3Tpkap2 - 4*DyCYvi13Tvi23*TrCTk3Tpkap3 -& 
&  6*CYvTpYv*TrCTuTpTu - 6*CYvTpTv*TrCTuTpYu - 2*CYvTpYv*TrCTvTpTv - 2*CYvTpTv*TrCTvTpYv -& 
&  6*adjYeYe*Trmd2YdadjYd - 2*adjYeYe*Trme2YeadjYe - 2*adjYeYe*Trml2adjYeYe -            & 
&  6*adjYeYe*Trmq2adjYdYd - 6*CYvTpYv*Trmq2adjYuYu - 6*CYvTpYv*Trmu2YuadjYu -            & 
&  6*adjYeme2Ye*TrYdadjYd - 3*adjYeYeml2*TrYdadjYd - 12*adjYeYe*mHd2*TrYdadjYd -         & 
&  3*ml2adjYeYe*TrYdadjYd - 2*adjYeme2Ye*TrYeadjYe - adjYeYeml2*TrYeadjYe -              & 
&  4*adjYeYe*mHd2*TrYeadjYe - ml2adjYeYe*TrYeadjYe + (2*adjTeTe*(6._dp*(g1p2) -          & 
&  5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)))/5._dp - 6*CTvTpTv*TrYuadjYu -       & 
&  6*CYvmv2TpYv*TrYuadjYu - 3*CYvTpYvml2*TrYuadjYu + 3*DyCYvlami1mlHd2i2*TrYuadjYu -     & 
&  12*CYvTpYv*mHu2*TrYuadjYu - 3*ml2CYvTpYv*TrYuadjYu - 2*CTvTpTv*TrYvadjYv -            & 
&  2*CYvmv2TpYv*TrYvadjYv - CYvTpYvml2*TrYvadjYv + DyCYvlami1mlHd2i2*TrYvadjYv -         & 
&  4*CYvTpYv*mHu2*TrYvadjYv - ml2CYvTpYv*TrYvadjYv - 2*CYvTpYv*TrYvadjYvCml2  
betaml22 =  betaml22- 2*CYvTpYv*TrYvCmv2adjYv - (12*adjYeTe*g1p2*Conjg(M1))/5._dp + (9*g1p2*g2p2*id3R*M2*Conjg(M1))/5._dp +& 
&  (9*g1p2*g2p2*id3R*M1*Conjg(M2))/5._dp + 6*g2p4*id3R*Tr2(2) + (6*g1p2*id3R*Tr2U1(1,    & 
& 1))/5._dp - 4*g1*id3R*sqrt3ov5*Tr3(1) + 2*DyCYvTpkap1Ckap1i11mlHd2i2*lam(1) +          & 
&  2*DyCYvTpkap2Ckap1i12mlHd2i2*lam(1) + 2*DyCYvTpkap3Ckap1i13mlHd2i2*lam(1) +           & 
&  2*DyCYvTpkap1Ckap2i11mlHd2i2*lam(2) + 2*DyCYvTpkap2Ckap2i12mlHd2i2*lam(2) +           & 
&  2*DyCYvTpkap3Ckap2i13mlHd2i2*lam(2) + 2*DyCYvTpkap1Ckap3i11mlHd2i2*lam(3) +           & 
&  2*DyCYvTpkap2Ckap3i12mlHd2i2*lam(3) + 2*DyCYvTpkap3Ckap3i13mlHd2i2*lam(3)

 
Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 ) 

 
Else 
Dml2 = oo16pi2* betaml21 
End If 
 
 
Call Chop(Dml2) 

Forall(i1=1:3) Dml2(i1,i1) =  Real(Dml2(i1,i1),dp) 
Dml2 = 0.5_dp*(Dml2+ Conjg(Transpose(Dml2)) ) 
!-------------------- 
! mHd2 
!-------------------- 
 
betamHd21  = (-6*AbsM1*g1p2)/5._dp - 6*AbsM2*g2p2 - SPlamxxadjYvmlHd2 +               & 
&  2*mHd2*SPlamxxClam + 2*mHu2*SPlamxxClam + 2._dp*(SPlamxxCmv2Clam) + 2._dp*(SPTlamxxCTlam)& 
&  + 6._dp*(TrCTdTpTd) + 2._dp*(TrCTeTpTe) + 6._dp*(Trmd2YdadjYd) + 2._dp*(Trme2YeadjYe) & 
&  + 2._dp*(Trml2adjYeYe) + 6._dp*(Trmq2adjYdYd) + 6*mHd2*TrYdadjYd + 2*mHd2*TrYeadjYe - & 
&  g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHd22 = (621*AbsM1*g1p4)/25._dp + (18*AbsM1*g1p2*g2p2)/5._dp + (18*AbsM2*g1p2*g2p2)/5._dp +   & 
&  33*AbsM2*g2p4 - 4._dp*(SPadjTvYvClamxxTlam) - 4._dp*(SPadjYvYvCTlamxxTlam) -          & 
&  4._dp*(SPClamxxTpYvml2CYvlam) - 4._dp*(SPlamxxadjTvTvClam) + 2._dp*(SPlamxxadjYvYvadjYvmlHd2) -& 
&  4*mHd2*SPlamxxadjYvYvClam - 8*mHu2*SPlamxxadjYvYvClam - 4._dp*(SPlamxxadjYvYvCmv2Clam) +& 
&  9*SPlamxxadjYvmlHd2*SPlamxxClam - 12*mHd2*SPlamxxClamp2 - 12*mHu2*SPlamxxClamp2 -     & 
&  4._dp*(SPlamxxCmv2adjYvYvClam) - 12*SPlamxxClam*SPlamxxCmv2Clam - 12*SPClamxxTlam*SPlamxxCTlam +& 
&  SPmlHd2xxadjYeYeCYvlam - 12*SPlamxxClam*SPTlamxxCTlam - 4._dp*(SPTpTk1ClamxxadjTk1lam) -& 
&  4._dp*(SPTpTk2ClamxxadjTk2lam) - 4._dp*(SPTpTk3ClamxxadjTk3lam) - 4._dp*(SPTpTvCYvlamxxCTlam) -& 
&  2._dp*(TradjYeTeCTvTpYv) - 6*SPlamxxCTlam*TradjYuTu - 2._dp*(TradjYvTpYeCTeTv) -      & 
&  2*SPlamxxCTlam*TradjYvTv - (4*g1p2*TrCTdTpTd)/5._dp + 32*g3p2*TrCTdTpTd +             & 
&  (4*g1p2*M1*TrCTdTpYd)/5._dp - 32*g3p2*M3*TrCTdTpYd + (12*g1p2*TrCTeTpTe)/5._dp -      & 
&  (12*g1p2*M1*TrCTeTpYe)/5._dp - 6*SPlamxxClam*TrCTuTpTu - 6*SPClamxxTlam*TrCTuTpYu -   & 
&  2*SPlamxxClam*TrCTvTpTv - 2*SPClamxxTlam*TrCTvTpYv - (4*g1p2*Trmd2YdadjYd)/5._dp +    & 
&  32*g3p2*Trmd2YdadjYd - 36._dp*(Trmd2YdadjYdYdadjYd) - 6._dp*(Trmd2YdadjYuYuadjYd) +   & 
&  (12*g1p2*Trme2YeadjYe)/5._dp - 12._dp*(Trme2YeadjYeYeadjYe) + (12*g1p2*Trml2adjYeYe)/5._dp -& 
&  12._dp*(Trml2adjYeYeadjYeYe) - (4*g1p2*Trmq2adjYdYd)/5._dp + 32*g3p2*Trmq2adjYdYd -   & 
&  36._dp*(Trmq2adjYdYdadjYdYd) - 6._dp*(Trmq2adjYdYdadjYuYu) - 6*SPlamxxClam*Trmq2adjYuYu -& 
&  6._dp*(Trmq2adjYuYuadjYdYd) - 6._dp*(Trmu2YuadjYdYdadjYu) - 6*SPlamxxClam*Trmu2YuadjYu -& 
&  36._dp*(TrYdadjTdTdadjYd) - 6._dp*(TrYdadjTuTuadjYd) - (8*AbsM1*g1p2*TrYdadjYd)/5._dp +& 
&  64*AbsM3*g3p2*TrYdadjYd - (4*g1p2*mHd2*TrYdadjYd)/5._dp + 32*g3p2*mHd2*TrYdadjYd -    & 
&  36._dp*(TrYdadjYdTdadjTd) - 36*mHd2*TrYdadjYdYdadjYd - 6._dp*(TrYdadjYuTuadjTd)  
betamHd22 =  betamHd22- 6*mHd2*TrYdadjYuYuadjYd - 6*mHu2*TrYdadjYuYuadjYd - 12._dp*(TrYeadjTeTeadjYe) +       & 
&  (24*AbsM1*g1p2*TrYeadjYe)/5._dp + (12*g1p2*mHd2*TrYeadjYe)/5._dp - 12._dp*(TrYeadjYeTeadjTe) -& 
&  12*mHd2*TrYeadjYeYeadjYe - 2._dp*(TrYeCTvTpTvadjYe) - 6._dp*(TrYuadjTdTdadjYu) -      & 
&  6._dp*(TrYuadjYdTdadjTu) + 3*SPlamxxadjYvmlHd2*TrYuadjYu - 6*mHd2*SPlamxxClam*TrYuadjYu -& 
&  12*mHu2*SPlamxxClam*TrYuadjYu - 6*SPlamxxCmv2Clam*TrYuadjYu - 6*SPTlamxxCTlam*TrYuadjYu +& 
&  SPlamxxadjYvmlHd2*TrYvadjYv - 2*mHd2*SPlamxxClam*TrYvadjYv - 4*mHu2*SPlamxxClam*TrYvadjYv -& 
&  2*SPlamxxCmv2Clam*TrYvadjYv - 2*SPTlamxxCTlam*TrYvadjYv - 2*SPlamxxClam*TrYvadjYvCml2 -& 
&  2._dp*(TrYvadjYvCml2TpYeCYe) - 2._dp*(TrYvadjYvTpTeCTe) - 2._dp*(TrYvadjYvTpYeCme2CYe) -& 
&  2*mHd2*TrYvadjYvTpYeCYe - 2*mHu2*TrYvadjYvTpYeCYe - 2._dp*(TrYvadjYvTpYeCYeCml2) -    & 
&  2*SPlamxxClam*TrYvCmv2adjYv - 2._dp*(TrYvCmv2adjYvTpYeCYe) + (g1p2*(9*g2p2*M2 +       & 
&  4*(TradjYdTd - 3._dp*(TradjYeTe)))*Conjg(M1))/5._dp - 32*g3p2*TradjYdTd*Conjg(M3) +   & 
&  (9*g1p2*g2p2*M1*Conjg(M2))/5._dp + 2*adjYvmlHd2(1)*Tpkap1adjkap1lam(1) -              & 
&  4*Conjg(lam(1))*Tpkap1adjkap1mv2lam(1) + 2*adjYvmlHd2(1)*Tpkap1adjkap2lam(2) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap2mv2lam(2) + 2*adjYvmlHd2(1)*Tpkap1adjkap3lam(3) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap3mv2lam(3) - 2*mHd2*Conjg(lam(1))*Tpkap1Ckap1lam(1) -     & 
&  4*mHu2*Conjg(lam(1))*Tpkap1Ckap1lam(1) - 2*mHd2*Conjg(lam(1))*Tpkap1Ckap2lam(2) -     & 
&  4*mHu2*Conjg(lam(1))*Tpkap1Ckap2lam(2) - 2*mHd2*Conjg(lam(1))*Tpkap1Ckap3lam(3) -     & 
&  4*mHu2*Conjg(lam(1))*Tpkap1Ckap3lam(3) + 2*adjYvmlHd2(2)*Tpkap2adjkap1lam(1) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap1mv2lam(1) + 2*adjYvmlHd2(2)*Tpkap2adjkap2lam(2) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap2mv2lam(2) + 2*adjYvmlHd2(2)*Tpkap2adjkap3lam(3) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap3mv2lam(3) - 2*mHd2*Conjg(lam(2))*Tpkap2Ckap1lam(1)  
betamHd22 =  betamHd22- 4*mHu2*Conjg(lam(2))*Tpkap2Ckap1lam(1) - 2*mHd2*Conjg(lam(2))*Tpkap2Ckap2lam(2) -     & 
&  4*mHu2*Conjg(lam(2))*Tpkap2Ckap2lam(2) - 2*mHd2*Conjg(lam(2))*Tpkap2Ckap3lam(3) -     & 
&  4*mHu2*Conjg(lam(2))*Tpkap2Ckap3lam(3) + 2*adjYvmlHd2(3)*Tpkap3adjkap1lam(1) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap1mv2lam(1) + 2*adjYvmlHd2(3)*Tpkap3adjkap2lam(2) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap2mv2lam(2) + 2*adjYvmlHd2(3)*Tpkap3adjkap3lam(3) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap3mv2lam(3) - 2*mHd2*Conjg(lam(3))*Tpkap3Ckap1lam(1) -     & 
&  4*mHu2*Conjg(lam(3))*Tpkap3Ckap1lam(1) - 2*mHd2*Conjg(lam(3))*Tpkap3Ckap2lam(2) -     & 
&  4*mHu2*Conjg(lam(3))*Tpkap3Ckap2lam(2) - 2*mHd2*Conjg(lam(3))*Tpkap3Ckap3lam(3) -     & 
&  4*mHu2*Conjg(lam(3))*Tpkap3Ckap3lam(3) - 4*Conjg(Tlam(1))*TpTk1adjkap1lam(1) -        & 
&  4*Conjg(Tlam(1))*TpTk1adjkap2lam(2) - 4*Conjg(Tlam(1))*TpTk1adjkap3lam(3) -           & 
&  4*Conjg(Tlam(2))*TpTk2adjkap1lam(1) - 4*Conjg(Tlam(2))*TpTk2adjkap2lam(2) -           & 
&  4*Conjg(Tlam(2))*TpTk2adjkap3lam(3) - 4*Conjg(Tlam(3))*TpTk3adjkap1lam(1) -           & 
&  4*Conjg(Tlam(3))*TpTk3adjkap2lam(2) - 4*Conjg(Tlam(3))*TpTk3adjkap3lam(3) +           & 
&  6*g2p4*Tr2(2) + (6*g1p2*Tr2U1(1,1))/5._dp - 4*g1*sqrt3ov5*Tr3(1) - 2*mHd2*adjkap1kap1Clam(1)*lam(1) -& 
&  2*mHd2*adjkap1kap2Clam(2)*lam(1) - 2*mHd2*adjkap1kap3Clam(3)*lam(1) - 4*adjkap1Tpkap1Cmv2Clam(1)*lam(1) -& 
&  4*adjkap1Tpkap2Cmv2Clam(2)*lam(1) - 4*adjkap1Tpkap3Cmv2Clam(3)*lam(1) -               & 
&  8*Trmv2kap1adjkap1*Conjg(lam(1))*lam(1) - 8*Trmv2kap2adjkap1*Conjg(lam(2))*lam(1) -   & 
&  8*Trmv2kap3adjkap1*Conjg(lam(3))*lam(1) - 2*mHd2*adjkap2kap1Clam(1)*lam(2) -          & 
&  2*mHd2*adjkap2kap2Clam(2)*lam(2) - 2*mHd2*adjkap2kap3Clam(3)*lam(2) - 4*adjkap2Tpkap1Cmv2Clam(1)*lam(2) -& 
&  4*adjkap2Tpkap2Cmv2Clam(2)*lam(2) - 4*adjkap2Tpkap3Cmv2Clam(3)*lam(2) -               & 
&  8*Trmv2kap1adjkap2*Conjg(lam(1))*lam(2) - 8*Trmv2kap2adjkap2*Conjg(lam(2))*lam(2)  
betamHd22 =  betamHd22- 8*Trmv2kap3adjkap2*Conjg(lam(3))*lam(2) - 2*mHd2*adjkap3kap1Clam(1)*lam(3) -          & 
&  2*mHd2*adjkap3kap2Clam(2)*lam(3) - 2*mHd2*adjkap3kap3Clam(3)*lam(3) - 4*adjkap3Tpkap1Cmv2Clam(1)*lam(3) -& 
&  4*adjkap3Tpkap2Cmv2Clam(2)*lam(3) - 4*adjkap3Tpkap3Cmv2Clam(3)*lam(3) -               & 
&  8*Trmv2kap1adjkap3*Conjg(lam(1))*lam(3) - 8*Trmv2kap2adjkap3*Conjg(lam(2))*lam(3) -   & 
&  8*Trmv2kap3adjkap3*Conjg(lam(3))*lam(3) - 4*TrCTk1Tpkap1*Conjg(lam(1))*Tlam(1) -      & 
&  4*TrCTk1Tpkap2*Conjg(lam(2))*Tlam(1) - 4*TrCTk1Tpkap3*Conjg(lam(3))*Tlam(1) -         & 
&  4*TrCkap1Tpkap1*Conjg(Tlam(1))*Tlam(1) - 4*TrCkap1Tpkap2*Conjg(Tlam(2))*Tlam(1) -     & 
&  4*TrCkap1Tpkap3*Conjg(Tlam(3))*Tlam(1) - 4*TrCTk2Tpkap1*Conjg(lam(1))*Tlam(2) -       & 
&  4*TrCTk2Tpkap2*Conjg(lam(2))*Tlam(2) - 4*TrCTk2Tpkap3*Conjg(lam(3))*Tlam(2) -         & 
&  4*TrCkap2Tpkap1*Conjg(Tlam(1))*Tlam(2) - 4*TrCkap2Tpkap2*Conjg(Tlam(2))*Tlam(2) -     & 
&  4*TrCkap2Tpkap3*Conjg(Tlam(3))*Tlam(2) - 4*TrCTk3Tpkap1*Conjg(lam(1))*Tlam(3) -       & 
&  4*TrCTk3Tpkap2*Conjg(lam(2))*Tlam(3) - 4*TrCTk3Tpkap3*Conjg(lam(3))*Tlam(3) -         & 
&  4*TrCkap3Tpkap1*Conjg(Tlam(1))*Tlam(3) - 4*TrCkap3Tpkap2*Conjg(Tlam(2))*Tlam(3) -     & 
&  4*TrCkap3Tpkap3*Conjg(Tlam(3))*Tlam(3)

 
DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 ) 

 
Else 
DmHd2 = oo16pi2* betamHd21 
End If 
 
 
!-------------------- 
! mHu2 
!-------------------- 
 
betamHu21  = (-6*AbsM1*g1p2)/5._dp - 6*AbsM2*g2p2 - 2._dp*(SPlamxxadjYvmlHd2)         & 
&  + 2*mHd2*SPlamxxClam + 2*mHu2*SPlamxxClam + 2._dp*(SPlamxxCmv2Clam) + 2._dp*(SPTlamxxCTlam)& 
&  + 6._dp*(TrCTuTpTu) + 2._dp*(TrCTvTpTv) + 6._dp*(Trmq2adjYuYu) + 6._dp*(Trmu2YuadjYu) & 
&  + 6*mHu2*TrYuadjYu + 2*mHu2*TrYvadjYv + 2._dp*(TrYvadjYvCml2) + 2._dp*(TrYvCmv2adjYv) & 
&  + g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHu22 = (621*AbsM1*g1p4)/25._dp + (18*AbsM1*g1p2*g2p2)/5._dp + (18*AbsM2*g1p2*g2p2)/5._dp +   & 
&  33*AbsM2*g2p4 - 12._dp*(SPadjTvYvClamxxTlam) - 12._dp*(SPadjYvYvCTlamxxTlam) -        & 
&  12._dp*(SPClamxxTpYvml2CYvlam) - 12._dp*(SPlamxxadjTvTvClam) + 12._dp*(SPlamxxadjYvYvadjYvmlHd2) -& 
&  12*mHd2*SPlamxxadjYvYvClam - 24*mHu2*SPlamxxadjYvYvClam - 12._dp*(SPlamxxadjYvYvCmv2Clam) +& 
&  12*SPlamxxadjYvmlHd2*SPlamxxClam - 12*mHd2*SPlamxxClamp2 - 12*mHu2*SPlamxxClamp2 -    & 
&  12._dp*(SPlamxxCmv2adjYvYvClam) - 12*SPlamxxClam*SPlamxxCmv2Clam - 12*SPClamxxTlam*SPlamxxCTlam -& 
&  12*SPlamxxClam*SPTlamxxCTlam - 4._dp*(SPTpTk1ClamxxadjTk1lam) - 4._dp*(SPTpTk2ClamxxadjTk2lam) -& 
&  4._dp*(SPTpTk3ClamxxadjTk3lam) - 12._dp*(SPTpTvCYvlamxxCTlam) - 6*SPlamxxCTlam*TradjYdTd -& 
&  2*SPlamxxCTlam*TradjYeTe - 2._dp*(TradjYeTeCTvTpYv) - 2._dp*(TradjYvTpYeCTeTv) -      & 
&  6*SPlamxxClam*TrCTdTpTd - 6*SPClamxxTlam*TrCTdTpYd - 2*SPlamxxClam*TrCTeTpTe -        & 
&  2*SPClamxxTlam*TrCTeTpYe - 4*TradjYvTv*TrCTk1Tpkap1 - 4*TradjYvTv*TrCTk2Tpkap2 -      & 
&  4*TradjYvTv*TrCTk3Tpkap3 + (8*g1p2*TrCTuTpTu)/5._dp + 32*g3p2*TrCTuTpTu -             & 
&  (8*g1p2*M1*TrCTuTpYu)/5._dp - 32*g3p2*M3*TrCTuTpYu - 4*TrCkap1Tpkap1*TrCTvTpTv -      & 
&  4*TrCkap2Tpkap2*TrCTvTpTv - 4*TrCkap3Tpkap3*TrCTvTpTv - 6*SPlamxxClam*Trmd2YdadjYd -  & 
&  6._dp*(Trmd2YdadjYuYuadjYd) - 2*SPlamxxClam*Trme2YeadjYe - 2*SPlamxxClam*Trml2adjYeYe -& 
&  6*SPlamxxClam*Trmq2adjYdYd - 6._dp*(Trmq2adjYdYdadjYuYu) + (8*g1p2*Trmq2adjYuYu)/5._dp +& 
&  32*g3p2*Trmq2adjYuYu - 6._dp*(Trmq2adjYuYuadjYdYd) - 36._dp*(Trmq2adjYuYuadjYuYu) -   & 
&  6._dp*(Trmu2YuadjYdYdadjYu) + (8*g1p2*Trmu2YuadjYu)/5._dp + 32*g3p2*Trmu2YuadjYu -    & 
&  36._dp*(Trmu2YuadjYuYuadjYu) - 8._dp*(Trmv2kap1adjYvYvadjkap1) - 8._dp*(Trmv2kap2adjYvYvadjkap2) -& 
&  8._dp*(Trmv2kap3adjYvYvadjkap3) - 6._dp*(TrYdadjTuTuadjYd) + 6*SPlamxxadjYvmlHd2*TrYdadjYd -& 
&  12*mHd2*SPlamxxClam*TrYdadjYd - 6*mHu2*SPlamxxClam*TrYdadjYd - 6*SPlamxxCmv2Clam*TrYdadjYd  
betamHu22 =  betamHu22- 6*SPTlamxxCTlam*TrYdadjYd - 6._dp*(TrYdadjYuTuadjTd) - 6*mHd2*TrYdadjYuYuadjYd -      & 
&  6*mHu2*TrYdadjYuYuadjYd + 2*SPlamxxadjYvmlHd2*TrYeadjYe - 4*mHd2*SPlamxxClam*TrYeadjYe -& 
&  2*mHu2*SPlamxxClam*TrYeadjYe - 2*SPlamxxCmv2Clam*TrYeadjYe - 2*SPTlamxxCTlam*TrYeadjYe -& 
&  2._dp*(TrYeCTvTpTvadjYe) - 6._dp*(TrYuadjTdTdadjYu) - 36._dp*(TrYuadjTuTuadjYu) -     & 
&  6._dp*(TrYuadjYdTdadjTu) + (16*AbsM1*g1p2*TrYuadjYu)/5._dp + 64*AbsM3*g3p2*TrYuadjYu +& 
&  (8*g1p2*mHu2*TrYuadjYu)/5._dp + 32*g3p2*mHu2*TrYuadjYu - 36._dp*(TrYuadjYuTuadjTu) -  & 
&  36*mHu2*TrYuadjYuYuadjYu - 12._dp*(TrYvadjTvTvadjYv) - 2._dp*(TrYvadjYvCml2TpYeCYe) - & 
&  12._dp*(TrYvadjYvCml2YvadjYv) - 2._dp*(TrYvadjYvTpTeCTe) - 2._dp*(TrYvadjYvTpYeCme2CYe) -& 
&  2*mHd2*TrYvadjYvTpYeCYe - 2*mHu2*TrYvadjYvTpYeCYe - 2._dp*(TrYvadjYvTpYeCYeCml2) -    & 
&  12._dp*(TrYvadjYvTvadjTv) - 12*mHu2*TrYvadjYvYvadjYv - 2*mHu2*TrYvCkap1kap1adjYv -    & 
&  4._dp*(TrYvCkap1kap1adjYvCml2) - 4._dp*(TrYvCkap1Tk1adjTv) - 2*mHu2*TrYvCkap1Tpkap1adjYv -& 
&  4._dp*(TrYvCkap1Tpkap1Cmv2adjYv) - 2*mHu2*TrYvCkap2kap2adjYv - 4._dp*(TrYvCkap2kap2adjYvCml2) -& 
&  4._dp*(TrYvCkap2Tk2adjTv) - 2*mHu2*TrYvCkap2Tpkap2adjYv - 4._dp*(TrYvCkap2Tpkap2Cmv2adjYv) -& 
&  2*mHu2*TrYvCkap3kap3adjYv - 4._dp*(TrYvCkap3kap3adjYvCml2) - 4._dp*(TrYvCkap3Tk3adjTv) -& 
&  2*mHu2*TrYvCkap3Tpkap3adjYv - 4._dp*(TrYvCkap3Tpkap3Cmv2adjYv) - 2._dp*(TrYvCmv2adjYvTpYeCYe) -& 
&  12._dp*(TrYvCmv2adjYvYvadjYv) - 4._dp*(TrYvCmv2Ckap1kap1adjYv) - 4._dp*(TrYvCmv2Ckap2kap2adjYv) -& 
&  4._dp*(TrYvCmv2Ckap3kap3adjYv) - 4._dp*(TrYvCTk1TpTk1adjYv) - 4._dp*(TrYvCTk2TpTk2adjYv) -& 
&  4._dp*(TrYvCTk3TpTk3adjYv) - 4*adjTvYvCkap1Tk2(2,1) - 4*adjTvYvCkap1Tk3(3,            & 
& 1) - 4*adjTvYvCkap2Tk1(1,2) - 4*adjTvYvCkap2Tk3(3,2) - 4*adjTvYvCkap3Tk1(1,            & 
& 3) - 4*adjTvYvCkap3Tk2(2,3) - 4*adjYvCml2YvCkap1kap2(2,1) - 4*adjYvCml2YvCkap1kap3(3,  & 
& 1) - 4*adjYvCml2YvCkap2kap1(1,2) - 4*adjYvCml2YvCkap2kap3(3,2) - 4*adjYvCml2YvCkap3kap1(1,& 
& 3) - 4*adjYvCml2YvCkap3kap2(2,3) - 4*TrCTk2Tpkap1*adjYvTv(1,2) - 4*TrCTk3Tpkap1*adjYvTv(1,& 
& 3) - 4*TrCTk1Tpkap2*adjYvTv(2,1) - 4*TrCTk3Tpkap2*adjYvTv(2,3) - 4*TrCTk1Tpkap3*adjYvTv(3,& 
& 1) - 4*TrCTk2Tpkap3*adjYvTv(3,2) - 8*adjYvYvadjkap1mv2kap2(2,1) - 8*adjYvYvadjkap1mv2kap3(3,& 
& 1) - 8*adjYvYvadjkap2mv2kap1(1,2) - 8*adjYvYvadjkap2mv2kap3(3,2) - 8*adjYvYvadjkap3mv2kap1(1,& 
& 3) - 8*adjYvYvadjkap3mv2kap2(2,3) - 2*mHu2*adjYvYvCkap1kap2(2,1) - 2*mHu2*adjYvYvCkap1kap3(3,& 
& 1) - 2*mHu2*adjYvYvCkap2kap1(1,2) - 2*mHu2*adjYvYvCkap2kap3(3,2) - 2*mHu2*adjYvYvCkap3kap1(1,& 
& 3) - 2*mHu2*adjYvYvCkap3kap2(2,3) - 4*adjYvYvCmv2Ckap1kap2(2,1) - 4*adjYvYvCmv2Ckap1kap3(3,& 
& 1) - 4*adjYvYvCmv2Ckap2kap1(1,2) - 4*adjYvYvCmv2Ckap2kap3(3,2) - 4*adjYvYvCmv2Ckap3kap1(1,& 
& 3) - 4*adjYvYvCmv2Ckap3kap2(2,3) + (g1p2*(9*g2p2*M2 - 8._dp*(TradjYuTu))*Conjg(M1))/5._dp  
betamHu22 =  betamHu22- 32*g3p2*TradjYuTu*Conjg(M3) + (9*g1p2*g2p2*M1*Conjg(M2))/5._dp + 4*adjYvmlHd2(1)*Tpkap1adjkap1lam(1) -& 
&  4*Conjg(lam(1))*Tpkap1adjkap1mv2lam(1) + 4*adjYvmlHd2(1)*Tpkap1adjkap2lam(2) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap2mv2lam(2) + 4*adjYvmlHd2(1)*Tpkap1adjkap3lam(3) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap3mv2lam(3) - 4*mHd2*Conjg(lam(1))*Tpkap1Ckap1lam(1) -     & 
&  2*mHu2*Conjg(lam(1))*Tpkap1Ckap1lam(1) - 4*mHd2*Conjg(lam(1))*Tpkap1Ckap2lam(2) -     & 
&  2*mHu2*Conjg(lam(1))*Tpkap1Ckap2lam(2) - 4*mHd2*Conjg(lam(1))*Tpkap1Ckap3lam(3) -     & 
&  2*mHu2*Conjg(lam(1))*Tpkap1Ckap3lam(3) + 4*adjYvmlHd2(2)*Tpkap2adjkap1lam(1) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap1mv2lam(1) + 4*adjYvmlHd2(2)*Tpkap2adjkap2lam(2) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap2mv2lam(2) + 4*adjYvmlHd2(2)*Tpkap2adjkap3lam(3) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap3mv2lam(3) - 4*mHd2*Conjg(lam(2))*Tpkap2Ckap1lam(1) -     & 
&  2*mHu2*Conjg(lam(2))*Tpkap2Ckap1lam(1) - 4*mHd2*Conjg(lam(2))*Tpkap2Ckap2lam(2) -     & 
&  2*mHu2*Conjg(lam(2))*Tpkap2Ckap2lam(2) - 4*mHd2*Conjg(lam(2))*Tpkap2Ckap3lam(3) -     & 
&  2*mHu2*Conjg(lam(2))*Tpkap2Ckap3lam(3) + 4*adjYvmlHd2(3)*Tpkap3adjkap1lam(1) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap1mv2lam(1) + 4*adjYvmlHd2(3)*Tpkap3adjkap2lam(2) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap2mv2lam(2) + 4*adjYvmlHd2(3)*Tpkap3adjkap3lam(3) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap3mv2lam(3) - 4*mHd2*Conjg(lam(3))*Tpkap3Ckap1lam(1) -     & 
&  2*mHu2*Conjg(lam(3))*Tpkap3Ckap1lam(1) - 4*mHd2*Conjg(lam(3))*Tpkap3Ckap2lam(2) -     & 
&  2*mHu2*Conjg(lam(3))*Tpkap3Ckap2lam(2) - 4*mHd2*Conjg(lam(3))*Tpkap3Ckap3lam(3) -     & 
&  2*mHu2*Conjg(lam(3))*Tpkap3Ckap3lam(3) - 4*Conjg(Tlam(1))*TpTk1adjkap1lam(1) -        & 
&  4*Conjg(Tlam(1))*TpTk1adjkap2lam(2) - 4*Conjg(Tlam(1))*TpTk1adjkap3lam(3) -           & 
&  4*Conjg(Tlam(2))*TpTk2adjkap1lam(1) - 4*Conjg(Tlam(2))*TpTk2adjkap2lam(2)  
betamHu22 =  betamHu22- 4*Conjg(Tlam(2))*TpTk2adjkap3lam(3) - 4*Conjg(Tlam(3))*TpTk3adjkap1lam(1) -           & 
&  4*Conjg(Tlam(3))*TpTk3adjkap2lam(2) - 4*Conjg(Tlam(3))*TpTk3adjkap3lam(3) -           & 
&  4*TrCkap1Tpkap2*TpTvCTv(1,2) - 4*TrCkap1Tpkap3*TpTvCTv(1,3) - 4*TrCkap2Tpkap1*TpTvCTv(2,& 
& 1) - 4*TrCkap2Tpkap3*TpTvCTv(2,3) - 4*TrCkap3Tpkap1*TpTvCTv(3,1) - 4*TrCkap3Tpkap2*TpTvCTv(3,& 
& 2) + 6*g2p4*Tr2(2) + (6*g1p2*Tr2U1(1,1))/5._dp + 4*g1*sqrt3ov5*Tr3(1) - 2*mHu2*adjkap1kap1Clam(1)*lam(1) -& 
&  2*mHu2*adjkap1kap2Clam(2)*lam(1) - 2*mHu2*adjkap1kap3Clam(3)*lam(1) - 4*adjkap1Tpkap1Cmv2Clam(1)*lam(1) -& 
&  4*adjkap1Tpkap2Cmv2Clam(2)*lam(1) - 4*adjkap1Tpkap3Cmv2Clam(3)*lam(1) -               & 
&  8*Trmv2kap1adjkap1*Conjg(lam(1))*lam(1) - 8*Trmv2kap2adjkap1*Conjg(lam(2))*lam(1) -   & 
&  8*Trmv2kap3adjkap1*Conjg(lam(3))*lam(1) - 2*mHu2*adjkap2kap1Clam(1)*lam(2) -          & 
&  2*mHu2*adjkap2kap2Clam(2)*lam(2) - 2*mHu2*adjkap2kap3Clam(3)*lam(2) - 4*adjkap2Tpkap1Cmv2Clam(1)*lam(2) -& 
&  4*adjkap2Tpkap2Cmv2Clam(2)*lam(2) - 4*adjkap2Tpkap3Cmv2Clam(3)*lam(2) -               & 
&  8*Trmv2kap1adjkap2*Conjg(lam(1))*lam(2) - 8*Trmv2kap2adjkap2*Conjg(lam(2))*lam(2) -   & 
&  8*Trmv2kap3adjkap2*Conjg(lam(3))*lam(2) - 2*mHu2*adjkap3kap1Clam(1)*lam(3) -          & 
&  2*mHu2*adjkap3kap2Clam(2)*lam(3) - 2*mHu2*adjkap3kap3Clam(3)*lam(3) - 4*adjkap3Tpkap1Cmv2Clam(1)*lam(3) -& 
&  4*adjkap3Tpkap2Cmv2Clam(2)*lam(3) - 4*adjkap3Tpkap3Cmv2Clam(3)*lam(3) -               & 
&  8*Trmv2kap1adjkap3*Conjg(lam(1))*lam(3) - 8*Trmv2kap2adjkap3*Conjg(lam(2))*lam(3) -   & 
&  8*Trmv2kap3adjkap3*Conjg(lam(3))*lam(3) - 4*TrCTk1Tpkap1*Conjg(lam(1))*Tlam(1) -      & 
&  4*TrCTk1Tpkap2*Conjg(lam(2))*Tlam(1) - 4*TrCTk1Tpkap3*Conjg(lam(3))*Tlam(1) -         & 
&  4*TrCkap1Tpkap1*Conjg(Tlam(1))*Tlam(1) - 4*TrCkap1Tpkap2*Conjg(Tlam(2))*Tlam(1) -     & 
&  4*TrCkap1Tpkap3*Conjg(Tlam(3))*Tlam(1) - 4*TrCTk2Tpkap1*Conjg(lam(1))*Tlam(2) -       & 
&  4*TrCTk2Tpkap2*Conjg(lam(2))*Tlam(2) - 4*TrCTk2Tpkap3*Conjg(lam(3))*Tlam(2)  
betamHu22 =  betamHu22- 4*TrCkap2Tpkap1*Conjg(Tlam(1))*Tlam(2) - 4*TrCkap2Tpkap2*Conjg(Tlam(2))*Tlam(2) -     & 
&  4*TrCkap2Tpkap3*Conjg(Tlam(3))*Tlam(2) - 4*TrCTk3Tpkap1*Conjg(lam(1))*Tlam(3) -       & 
&  4*TrCTk3Tpkap2*Conjg(lam(2))*Tlam(3) - 4*TrCTk3Tpkap3*Conjg(lam(3))*Tlam(3) -         & 
&  4*TrCkap3Tpkap1*Conjg(Tlam(1))*Tlam(3) - 4*TrCkap3Tpkap2*Conjg(Tlam(2))*Tlam(3) -     & 
&  4*TrCkap3Tpkap3*Conjg(Tlam(3))*Tlam(3)

 
DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 ) 

 
Else 
DmHu2 = oo16pi2* betamHu21 
End If 
 
 
!-------------------- 
! md2 
!-------------------- 
 
betamd21  = 2*(md2YdadjYd + 2._dp*(TdadjTd) + 2*mHd2*YdadjYd + YdadjYdmd2 +           & 
&  2._dp*(Ydmq2adjYd)) - (2*id3R*(4*AbsM1*g1p2 + 80*AbsM3*g3p2 - 15*g1*ooSqrt15*Tr1(1)))/15._dp

 
 
If (TwoLoopRGE) Then 
betamd22 = md2YdadjYd*(2._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 2*(SPlamxxClam + 3._dp*(TrYdadjYd) +  & 
&  TrYeadjYe)) - (2*(5._dp*(md2YdadjYdYdadjYd) + 5._dp*(md2YdadjYuYuadjYd) -             & 
&  2*g1p2*TdadjTd - 30*g2p2*TdadjTd + 10*SPlamxxClam*TdadjTd + 10._dp*(TdadjTdYdadjYd) + & 
&  10._dp*(TdadjTuYuadjYd) + 10*SPlamxxCTlam*TdadjYd + 10._dp*(TdadjYdYdadjTd) +         & 
&  10._dp*(TdadjYuYuadjTd) + 30*TdadjYd*TrCTdTpYd + 10*TdadjYd*TrCTeTpYe +               & 
&  30*TdadjTd*TrYdadjYd + 10*TdadjTd*TrYeadjYe + 2*g1p2*M1*YdadjTd + 30*g2p2*M2*YdadjTd +& 
&  10*SPClamxxTlam*YdadjTd + 30*TradjYdTd*YdadjTd + 10*TradjYeTe*YdadjTd +               & 
&  10._dp*(YdadjTdTdadjYd) + 10._dp*(YdadjTuTuadjYd) - 4*AbsM1*g1p2*YdadjYd -            & 
&  60*AbsM2*g2p2*YdadjYd - 2*g1p2*mHd2*YdadjYd - 30*g2p2*mHd2*YdadjYd - 10*SPlamxxadjYvmlHd2*YdadjYd +& 
&  20*mHd2*SPlamxxClam*YdadjYd + 10*mHu2*SPlamxxClam*YdadjYd + 10*SPlamxxCmv2Clam*YdadjYd +& 
&  10*SPTlamxxCTlam*YdadjYd + 30*TrCTdTpTd*YdadjYd + 10*TrCTeTpTe*YdadjYd +              & 
&  30*Trmd2YdadjYd*YdadjYd + 10*Trme2YeadjYe*YdadjYd + 10*Trml2adjYeYe*YdadjYd +         & 
&  30*Trmq2adjYdYd*YdadjYd + 60*mHd2*TrYdadjYd*YdadjYd + 20*mHd2*TrYeadjYe*YdadjYd -     & 
&  g1p2*YdadjYdmd2 - 15*g2p2*YdadjYdmd2 + 5*SPlamxxClam*YdadjYdmd2 + 15*TrYdadjYd*YdadjYdmd2 +& 
&  5*TrYeadjYe*YdadjYdmd2 + 10._dp*(YdadjYdmd2YdadjYd) + 10._dp*(YdadjYdTdadjTd) +       & 
&  20*mHd2*YdadjYdYdadjYd + 5._dp*(YdadjYdYdadjYdmd2) + 10._dp*(YdadjYdYdmq2adjYd) +     & 
&  10._dp*(YdadjYumu2YuadjYd) + 10._dp*(YdadjYuTuadjTd) + 10*mHd2*YdadjYuYuadjYd +       & 
&  10*mHu2*YdadjYuYuadjYd + 5._dp*(YdadjYuYuadjYdmd2) + 10._dp*(YdadjYuYumq2adjYd) -     & 
&  2*g1p2*Ydmq2adjYd - 30*g2p2*Ydmq2adjYd + 10*SPlamxxClam*Ydmq2adjYd + 30*TrYdadjYd*Ydmq2adjYd +& 
&  10*TrYeadjYe*Ydmq2adjYd + 10._dp*(Ydmq2adjYdYdadjYd) + 10._dp*(Ydmq2adjYuYuadjYd) +   & 
&  2*g1p2*TdadjYd*Conjg(M1) + 30*g2p2*TdadjYd*Conjg(M2)))/5._dp + (8*id3R*(303*AbsM1*g1p4 +& 
&  80*AbsM1*g1p2*g3p2 + 80*AbsM3*g1p2*g3p2 - 1200*AbsM3*g3p4 + 40*g1p2*g3p2*M3*Conjg(M1) +& 
&  40*g1p2*g3p2*M1*Conjg(M3) + 300*g3p4*Tr2(3) + 15*g1p2*Tr2U1(1,1) + 225*g1*ooSqrt15*Tr3(1)))/225._dp

 
Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 ) 

 
Else 
Dmd2 = oo16pi2* betamd21 
End If 
 
 
Call Chop(Dmd2) 

Forall(i1=1:3) Dmd2(i1,i1) =  Real(Dmd2(i1,i1),dp) 
Dmd2 = 0.5_dp*(Dmd2+ Conjg(Transpose(Dmd2)) ) 
!-------------------- 
! mu2 
!-------------------- 
 
betamu21  = 2*(mu2YuadjYu + 2._dp*(TuadjTu) + 2*mHu2*YuadjYu + YuadjYumu2 +           & 
&  2._dp*(Yumq2adjYu)) - (4*id3R*(8*AbsM1*g1p2 + 40*AbsM3*g3p2 + 15*g1*ooSqrt15*Tr1(1)))/15._dp

 
 
If (TwoLoopRGE) Then 
betamu22 = (-2*(5._dp*(mu2YuadjYdYdadjYu) + 5._dp*(mu2YuadjYuYuadjYu) + mu2YuadjYu*(g1p2 +       & 
&  5*(-3._dp*(g2p2) + SPlamxxClam + 3._dp*(TrYuadjYu) + TrYvadjYv)) + 10._dp*(TuadjTdYdadjYu) +& 
&  2*g1p2*TuadjTu - 30*g2p2*TuadjTu + 10*SPlamxxClam*TuadjTu + 30*TrYuadjYu*TuadjTu +    & 
&  10*TrYvadjYv*TuadjTu + 10._dp*(TuadjTuYuadjYu) + 10._dp*(TuadjYdYdadjTu) +            & 
&  10*SPlamxxCTlam*TuadjYu + 30*TrCTuTpYu*TuadjYu + 10*TrCTvTpYv*TuadjYu +               & 
&  10._dp*(TuadjYuYuadjTu) + 10._dp*(YuadjTdTdadjYu) - 2*g1p2*M1*YuadjTu +               & 
&  30*g2p2*M2*YuadjTu + 10*SPClamxxTlam*YuadjTu + 30*TradjYuTu*YuadjTu + 10*TradjYvTv*YuadjTu +& 
&  10._dp*(YuadjTuTuadjYu) + 10._dp*(YuadjYdmd2YdadjYu) + 10._dp*(YuadjYdTdadjTu) +      & 
&  10*mHd2*YuadjYdYdadjYu + 10*mHu2*YuadjYdYdadjYu + 5._dp*(YuadjYdYdadjYumu2) +         & 
&  10._dp*(YuadjYdYdmq2adjYu) + 4*AbsM1*g1p2*YuadjYu - 60*AbsM2*g2p2*YuadjYu +           & 
&  2*g1p2*mHu2*YuadjYu - 30*g2p2*mHu2*YuadjYu - 10*SPlamxxadjYvmlHd2*YuadjYu +           & 
&  10*mHd2*SPlamxxClam*YuadjYu + 20*mHu2*SPlamxxClam*YuadjYu + 10*SPlamxxCmv2Clam*YuadjYu +& 
&  10*SPTlamxxCTlam*YuadjYu + 30*TrCTuTpTu*YuadjYu + 10*TrCTvTpTv*YuadjYu +              & 
&  30*Trmq2adjYuYu*YuadjYu + 30*Trmu2YuadjYu*YuadjYu + 60*mHu2*TrYuadjYu*YuadjYu +       & 
&  20*mHu2*TrYvadjYv*YuadjYu + 10*TrYvadjYvCml2*YuadjYu + 10*TrYvCmv2adjYv*YuadjYu +     & 
&  g1p2*YuadjYumu2 - 15*g2p2*YuadjYumu2 + 5*SPlamxxClam*YuadjYumu2 + 15*TrYuadjYu*YuadjYumu2 +& 
&  5*TrYvadjYv*YuadjYumu2 + 10._dp*(YuadjYumu2YuadjYu) + 10._dp*(YuadjYuTuadjTu) +       & 
&  20*mHu2*YuadjYuYuadjYu + 5._dp*(YuadjYuYuadjYumu2) + 10._dp*(YuadjYuYumq2adjYu) +     & 
&  10._dp*(Yumq2adjYdYdadjYu) + 2*g1p2*Yumq2adjYu - 30*g2p2*Yumq2adjYu + 10*SPlamxxClam*Yumq2adjYu +& 
&  30*TrYuadjYu*Yumq2adjYu + 10*TrYvadjYv*Yumq2adjYu + 10._dp*(Yumq2adjYuYuadjYu) -      & 
&  2*g1p2*TuadjYu*Conjg(M1) + 30*g2p2*TuadjYu*Conjg(M2)))/5._dp + (16*id3R*(642*AbsM1*g1p4 +& 
&  160*AbsM1*g1p2*g3p2 + 160*AbsM3*g1p2*g3p2 - 600*AbsM3*g3p4 + 80*g1p2*g3p2*M3*Conjg(M1) +& 
&  80*g1p2*g3p2*M1*Conjg(M3) + 150*g3p4*Tr2(3) + 30*g1p2*Tr2U1(1,1) - 225*g1*ooSqrt15*Tr3(1)))/225._dp

 
Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 ) 

 
Else 
Dmu2 = oo16pi2* betamu21 
End If 
 
 
Call Chop(Dmu2) 

Forall(i1=1:3) Dmu2(i1,i1) =  Real(Dmu2(i1,i1),dp) 
Dmu2 = 0.5_dp*(Dmu2+ Conjg(Transpose(Dmu2)) ) 
!-------------------- 
! me2 
!-------------------- 
 
betame21  = 2*(me2YeadjYe + 2._dp*(TeadjTe) + 2*mHd2*YeadjYe + YeadjYeme2 +           & 
&  2._dp*(Yeml2adjYe)) + id3R*((-24*AbsM1*g1p2)/5._dp + 2*g1*sqrt3ov5*Tr1(1))

 
 
If (TwoLoopRGE) Then 
betame22 = (-2*(5._dp*(me2YeadjYeYeadjYe) + 5._dp*(me2YeCYvTpYvadjYe) + 6*g1p2*TeadjTe -         & 
&  30*g2p2*TeadjTe + 10*SPlamxxClam*TeadjTe + 10._dp*(TeadjTeYeadjYe) + 10*SPlamxxCTlam*TeadjYe +& 
&  10._dp*(TeadjYeYeadjTe) + 10._dp*(TeCTvTpYvadjYe) + 10._dp*(TeCYvTpYvadjTe) +         & 
&  30*TeadjYe*TrCTdTpYd + 10*TeadjYe*TrCTeTpYe + 30*TeadjTe*TrYdadjYd + 10*TeadjTe*TrYeadjYe +& 
&  me2YeadjYe*(3._dp*(g1p2) + 5*(-3._dp*(g2p2) + SPlamxxClam + 3._dp*(TrYdadjYd) +       & 
&  TrYeadjYe)) - 6*g1p2*M1*YeadjTe + 30*g2p2*M2*YeadjTe + 10*SPClamxxTlam*YeadjTe +      & 
&  30*TradjYdTd*YeadjTe + 10*TradjYeTe*YeadjTe + 10._dp*(YeadjTeTeadjYe) +               & 
&  12*AbsM1*g1p2*YeadjYe - 60*AbsM2*g2p2*YeadjYe + 6*g1p2*mHd2*YeadjYe - 30*g2p2*mHd2*YeadjYe -& 
&  10*SPlamxxadjYvmlHd2*YeadjYe + 20*mHd2*SPlamxxClam*YeadjYe + 10*mHu2*SPlamxxClam*YeadjYe +& 
&  10*SPlamxxCmv2Clam*YeadjYe + 10*SPTlamxxCTlam*YeadjYe + 30*TrCTdTpTd*YeadjYe +        & 
&  10*TrCTeTpTe*YeadjYe + 30*Trmd2YdadjYd*YeadjYe + 10*Trme2YeadjYe*YeadjYe +            & 
&  10*Trml2adjYeYe*YeadjYe + 30*Trmq2adjYdYd*YeadjYe + 60*mHd2*TrYdadjYd*YeadjYe +       & 
&  20*mHd2*TrYeadjYe*YeadjYe + 3*g1p2*YeadjYeme2 - 15*g2p2*YeadjYeme2 + 5*SPlamxxClam*YeadjYeme2 +& 
&  15*TrYdadjYd*YeadjYeme2 + 5*TrYeadjYe*YeadjYeme2 + 10._dp*(YeadjYeme2YeadjYe) +       & 
&  10._dp*(YeadjYeTeadjTe) + 20*mHd2*YeadjYeYeadjYe + 5._dp*(YeadjYeYeadjYeme2) +        & 
&  10._dp*(YeadjYeYeml2adjYe) + 10._dp*(YeCTvTpTvadjYe) + 10._dp*(YeCYvmv2TpYvadjYe) +   & 
&  10._dp*(YeCYvTpTvadjTe) + 10*mHd2*YeCYvTpYvadjYe + 10*mHu2*YeCYvTpYvadjYe +           & 
&  5._dp*(YeCYvTpYvadjYeme2) + 10._dp*(YeCYvTpYvml2adjYe) + 6*g1p2*Yeml2adjYe -          & 
&  30*g2p2*Yeml2adjYe + 10*SPlamxxClam*Yeml2adjYe + 30*TrYdadjYd*Yeml2adjYe +            & 
&  10*TrYeadjYe*Yeml2adjYe + 10._dp*(Yeml2adjYeYeadjYe) + 10._dp*(Yeml2CYvTpYvadjYe) -   & 
&  6*g1p2*TeadjYe*Conjg(M1) + 30*g2p2*TeadjYe*Conjg(M2)))/5._dp + (8*id3R*(351*AbsM1*g1p4 +& 
&  15*g1p2*Tr2U1(1,1) + 25*g1*sqrt3ov5*Tr3(1)))/25._dp

 
Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 ) 

 
Else 
Dme2 = oo16pi2* betame21 
End If 
 
 
Call Chop(Dme2) 

Forall(i1=1:3) Dme2(i1,i1) =  Real(Dme2(i1,i1),dp) 
Dme2 = 0.5_dp*(Dme2+ Conjg(Transpose(Dme2)) ) 
!-------------------- 
! mv2 
!-------------------- 
 
betamv21  = -4._dp*(Dylami1adjYvmlHd2i2) + 2._dp*(Dylami1Cmv2Clami2) + 2._dp*(Dymv2lami1Clami2)& 
&  + 4._dp*(DyTlami1CTlami2) + 2._dp*(kap1adjkap1mv2) + 8._dp*(kap1Cmv2adjkap1)          & 
&  + 2._dp*(kap2adjkap2mv2) + 8._dp*(kap2Cmv2adjkap2) + 2._dp*(kap3adjkap3mv2)           & 
&  + 8._dp*(kap3Cmv2adjkap3) + 4*Dylami1Clami2*mHd2 + 4*Dylami1Clami2*mHu2 +             & 
&  2._dp*(mv2kap1adjkap1) + 2._dp*(mv2kap2adjkap2) + 2._dp*(mv2kap3adjkap3)              & 
&  + 2._dp*(mv2TpYvCYv) + 4._dp*(Tk1adjTk1) + 4._dp*(Tk2adjTk2) + 4._dp*(Tk3adjTk3)      & 
&  + 4._dp*(TpTvCTv) + 4*mHu2*TpYvCYv + 2._dp*(TpYvCYvmv2) + 4._dp*(TpYvml2CYv)

 
 
If (TwoLoopRGE) Then 
betamv22 = -16._dp*(Dykap1Ckap1kap1i11Ckap1mv2i21) - 8._dp*(Dykap1Ckap1kap1i11Cmv2Ckap1i21) -    & 
&  16._dp*(Dykap1Ckap1kap2i11Ckap2mv2i21) - 8._dp*(Dykap1Ckap1kap2i11Cmv2Ckap2i21) -     & 
&  16._dp*(Dykap1Ckap1kap3i11Ckap3mv2i21) - 8._dp*(Dykap1Ckap1kap3i11Cmv2Ckap3i21) -     & 
&  16._dp*(Dykap1Ckap2kap1i12Ckap1mv2i21) - 8._dp*(Dykap1Ckap2kap1i12Cmv2Ckap1i21) -     & 
&  16._dp*(Dykap1Ckap2kap2i12Ckap2mv2i21) - 8._dp*(Dykap1Ckap2kap2i12Cmv2Ckap2i21) -     & 
&  16._dp*(Dykap1Ckap2kap3i12Ckap3mv2i21) - 8._dp*(Dykap1Ckap2kap3i12Cmv2Ckap3i21) -     & 
&  16._dp*(Dykap1Ckap3kap1i13Ckap1mv2i21) - 8._dp*(Dykap1Ckap3kap1i13Cmv2Ckap1i21) -     & 
&  16._dp*(Dykap1Ckap3kap2i13Ckap2mv2i21) - 8._dp*(Dykap1Ckap3kap2i13Cmv2Ckap2i21) -     & 
&  16._dp*(Dykap1Ckap3kap3i13Ckap3mv2i21) - 8._dp*(Dykap1Ckap3kap3i13Cmv2Ckap3i21) -     & 
&  16._dp*(Dykap1Clami1Ckap1mv2lami2) - 16._dp*(Dykap1Cmv2Ckap1kap1i11Ckap1i21) -        & 
&  16._dp*(Dykap1Cmv2Ckap1kap2i11Ckap2i21) - 16._dp*(Dykap1Cmv2Ckap1kap3i11Ckap3i21) -   & 
&  16._dp*(Dykap1Cmv2Ckap2kap1i12Ckap1i21) - 16._dp*(Dykap1Cmv2Ckap2kap2i12Ckap2i21) -   & 
&  16._dp*(Dykap1Cmv2Ckap2kap3i12Ckap3i21) - 16._dp*(Dykap1Cmv2Ckap3kap1i13Ckap1i21) -   & 
&  16._dp*(Dykap1Cmv2Ckap3kap2i13Ckap2i21) - 16._dp*(Dykap1Cmv2Ckap3kap3i13Ckap3i21) -   & 
&  16._dp*(Dykap1i11adjYvCml2Yvadjkap11i2) - 16._dp*(Dykap1i12adjYvCml2Yvadjkap21i2) -   & 
&  16._dp*(Dykap1i13adjYvCml2Yvadjkap31i2) - 16._dp*(Dykap2Ckap1kap1i11Ckap1mv2i22) -    & 
&  8._dp*(Dykap2Ckap1kap1i11Cmv2Ckap1i22) - 16._dp*(Dykap2Ckap1kap2i11Ckap2mv2i22) -     & 
&  8._dp*(Dykap2Ckap1kap2i11Cmv2Ckap2i22) - 16._dp*(Dykap2Ckap1kap3i11Ckap3mv2i22) -     & 
&  8._dp*(Dykap2Ckap1kap3i11Cmv2Ckap3i22) - 16._dp*(Dykap2Ckap2kap1i12Ckap1mv2i22) -     & 
&  8._dp*(Dykap2Ckap2kap1i12Cmv2Ckap1i22) - 16._dp*(Dykap2Ckap2kap2i12Ckap2mv2i22) -     & 
&  8._dp*(Dykap2Ckap2kap2i12Cmv2Ckap2i22) - 16._dp*(Dykap2Ckap2kap3i12Ckap3mv2i22)  
betamv22 =  betamv22- 8._dp*(Dykap2Ckap2kap3i12Cmv2Ckap3i22) - 16._dp*(Dykap2Ckap3kap1i13Ckap1mv2i22) -     & 
&  8._dp*(Dykap2Ckap3kap1i13Cmv2Ckap1i22) - 16._dp*(Dykap2Ckap3kap2i13Ckap2mv2i22) -     & 
&  8._dp*(Dykap2Ckap3kap2i13Cmv2Ckap2i22) - 16._dp*(Dykap2Ckap3kap3i13Ckap3mv2i22) -     & 
&  8._dp*(Dykap2Ckap3kap3i13Cmv2Ckap3i22) - 16._dp*(Dykap2Clami1Ckap2mv2lami2) -         & 
&  16._dp*(Dykap2Cmv2Ckap1kap1i11Ckap1i22) - 16._dp*(Dykap2Cmv2Ckap1kap2i11Ckap2i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap1kap3i11Ckap3i22) - 16._dp*(Dykap2Cmv2Ckap2kap1i12Ckap1i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap2kap2i12Ckap2i22) - 16._dp*(Dykap2Cmv2Ckap2kap3i12Ckap3i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap3kap1i13Ckap1i22) - 16._dp*(Dykap2Cmv2Ckap3kap2i13Ckap2i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap3kap3i13Ckap3i22) - 16._dp*(Dykap2i11adjYvCml2Yvadjkap12i2) -   & 
&  16._dp*(Dykap2i12adjYvCml2Yvadjkap22i2) - 16._dp*(Dykap2i13adjYvCml2Yvadjkap32i2) -   & 
&  16._dp*(Dykap3Ckap1kap1i11Ckap1mv2i23) - 8._dp*(Dykap3Ckap1kap1i11Cmv2Ckap1i23) -     & 
&  16._dp*(Dykap3Ckap1kap2i11Ckap2mv2i23) - 8._dp*(Dykap3Ckap1kap2i11Cmv2Ckap2i23) -     & 
&  16._dp*(Dykap3Ckap1kap3i11Ckap3mv2i23) - 8._dp*(Dykap3Ckap1kap3i11Cmv2Ckap3i23) -     & 
&  16._dp*(Dykap3Ckap2kap1i12Ckap1mv2i23) - 8._dp*(Dykap3Ckap2kap1i12Cmv2Ckap1i23) -     & 
&  16._dp*(Dykap3Ckap2kap2i12Ckap2mv2i23) - 8._dp*(Dykap3Ckap2kap2i12Cmv2Ckap2i23) -     & 
&  16._dp*(Dykap3Ckap2kap3i12Ckap3mv2i23) - 8._dp*(Dykap3Ckap2kap3i12Cmv2Ckap3i23) -     & 
&  16._dp*(Dykap3Ckap3kap1i13Ckap1mv2i23) - 8._dp*(Dykap3Ckap3kap1i13Cmv2Ckap1i23) -     & 
&  16._dp*(Dykap3Ckap3kap2i13Ckap2mv2i23) - 8._dp*(Dykap3Ckap3kap2i13Cmv2Ckap2i23) -     & 
&  16._dp*(Dykap3Ckap3kap3i13Ckap3mv2i23) - 8._dp*(Dykap3Ckap3kap3i13Cmv2Ckap3i23) -     & 
&  16._dp*(Dykap3Clami1Ckap3mv2lami2) - 16._dp*(Dykap3Cmv2Ckap1kap1i11Ckap1i23) -        & 
&  16._dp*(Dykap3Cmv2Ckap1kap2i11Ckap2i23) - 16._dp*(Dykap3Cmv2Ckap1kap3i11Ckap3i23)  
betamv22 =  betamv22- 16._dp*(Dykap3Cmv2Ckap2kap1i12Ckap1i23) - 16._dp*(Dykap3Cmv2Ckap2kap2i12Ckap2i23) -   & 
&  16._dp*(Dykap3Cmv2Ckap2kap3i12Ckap3i23) - 16._dp*(Dykap3Cmv2Ckap3kap1i13Ckap1i23) -   & 
&  16._dp*(Dykap3Cmv2Ckap3kap2i13Ckap2i23) - 16._dp*(Dykap3Cmv2Ckap3kap3i13Ckap3i23) -   & 
&  16._dp*(Dykap3i11adjYvCml2Yvadjkap13i2) - 16._dp*(Dykap3i12adjYvCml2Yvadjkap23i2) -   & 
&  16._dp*(Dykap3i13adjYvCml2Yvadjkap33i2) - 4._dp*(Dylami1adjTvTvClami2) -              & 
&  4._dp*(Dylami1adjYvCml2YvClami2) - 4._dp*(Dylami1adjYvTvCTlami2) + 4._dp*(Dylami1adjYvYvadjYvmlHd2i2) -& 
&  4._dp*(Dylami1adjYvYvCmv2Clami2) - 2._dp*(Dylami1Cmv2adjYvYvClami2) - 8._dp*(Dymv2kap1Clami1Ckap1lami2) -& 
&  8._dp*(Dymv2kap1i11Ckap1kap1Ckap1i21) - 8._dp*(Dymv2kap1i11Ckap1kap2Ckap1i22) -       & 
&  8._dp*(Dymv2kap1i11Ckap1kap3Ckap1i23) - 8._dp*(Dymv2kap1i12Ckap1kap1Ckap2i21) -       & 
&  8._dp*(Dymv2kap1i12Ckap1kap2Ckap2i22) - 8._dp*(Dymv2kap1i12Ckap1kap3Ckap2i23) -       & 
&  8._dp*(Dymv2kap1i13Ckap1kap1Ckap3i21) - 8._dp*(Dymv2kap1i13Ckap1kap2Ckap3i22) -       & 
&  8._dp*(Dymv2kap1i13Ckap1kap3Ckap3i23) - 8._dp*(Dymv2kap2Clami1Ckap2lami2) -           & 
&  8._dp*(Dymv2kap2i11Ckap2kap1Ckap1i21) - 8._dp*(Dymv2kap2i11Ckap2kap2Ckap1i22) -       & 
&  8._dp*(Dymv2kap2i11Ckap2kap3Ckap1i23) - 8._dp*(Dymv2kap2i12Ckap2kap1Ckap2i21) -       & 
&  8._dp*(Dymv2kap2i12Ckap2kap2Ckap2i22) - 8._dp*(Dymv2kap2i12Ckap2kap3Ckap2i23) -       & 
&  8._dp*(Dymv2kap2i13Ckap2kap1Ckap3i21) - 8._dp*(Dymv2kap2i13Ckap2kap2Ckap3i22) -       & 
&  8._dp*(Dymv2kap2i13Ckap2kap3Ckap3i23) - 8._dp*(Dymv2kap3Clami1Ckap3lami2) -           & 
&  8._dp*(Dymv2kap3i11Ckap3kap1Ckap1i21) - 8._dp*(Dymv2kap3i11Ckap3kap2Ckap1i22) -       & 
&  8._dp*(Dymv2kap3i11Ckap3kap3Ckap1i23) - 8._dp*(Dymv2kap3i12Ckap3kap1Ckap2i21) -       & 
&  8._dp*(Dymv2kap3i12Ckap3kap2Ckap2i22) - 8._dp*(Dymv2kap3i12Ckap3kap3Ckap2i23) -       & 
&  8._dp*(Dymv2kap3i13Ckap3kap1Ckap3i21) - 8._dp*(Dymv2kap3i13Ckap3kap2Ckap3i22)  
betamv22 =  betamv22- 8._dp*(Dymv2kap3i13Ckap3kap3Ckap3i23) - 2._dp*(Dymv2lami1adjYvYvClami2) -             & 
&  2._dp*(Dymv2TpYvCYvlami1Clami2) - 16._dp*(DyTk1Clami1CTk1lami2) - 16._dp*(DyTk1i11Ckap1TpYvCTvi21) -& 
&  16._dp*(DyTk1i12Ckap2TpYvCTvi21) - 16._dp*(DyTk1i13Ckap3TpYvCTvi21) - 16._dp*(DyTk2Clami1CTk2lami2) -& 
&  16._dp*(DyTk2i11Ckap1TpYvCTvi22) - 16._dp*(DyTk2i12Ckap2TpYvCTvi22) - 16._dp*(DyTk2i13Ckap3TpYvCTvi22) -& 
&  16._dp*(DyTk3Clami1CTk3lami2) - 16._dp*(DyTk3i11Ckap1TpYvCTvi23) - 16._dp*(DyTk3i12Ckap2TpYvCTvi23) -& 
&  16._dp*(DyTk3i13Ckap3TpYvCTvi23) - 4._dp*(DyTlami1adjTvYvClami2) - 4._dp*(DyTlami1adjYvYvCTlami2) -& 
&  4._dp*(DyTpTvCTvlami1Clami2) - 4._dp*(DyTpTvCYvlami1CTlami2) - 4._dp*(DyTpYvCTvTlami1Clami2) +& 
&  4._dp*(DyTpYvCYvlami1adjYvmlHd2i2) - 2._dp*(DyTpYvCYvlami1Cmv2Clami2) -               & 
&  4._dp*(DyTpYvCYvmv2lami1Clami2) - 4._dp*(DyTpYvCYvTlami1CTlami2) - 16._dp*(DyTpYvCYvTpkap11i1Ckap1mv2i21) -& 
&  16._dp*(DyTpYvCYvTpkap12i1Ckap2mv2i21) - 16._dp*(DyTpYvCYvTpkap13i1Ckap3mv2i21) -     & 
&  16._dp*(DyTpYvCYvTpkap21i1Ckap1mv2i22) - 16._dp*(DyTpYvCYvTpkap22i1Ckap2mv2i22) -     & 
&  16._dp*(DyTpYvCYvTpkap23i1Ckap3mv2i22) - 16._dp*(DyTpYvCYvTpkap31i1Ckap1mv2i23) -     & 
&  16._dp*(DyTpYvCYvTpkap32i1Ckap2mv2i23) - 16._dp*(DyTpYvCYvTpkap33i1Ckap3mv2i23) -     & 
&  4._dp*(DyTpYvml2CYvlami1Clami2) - (12*Dylami1adjYvmlHd2i2*g1p2)/5._dp +               & 
&  (24*AbsM1*Dylami1Clami2*g1p2)/5._dp + (6*Dylami1Cmv2Clami2*g1p2)/5._dp +              & 
&  (6*Dymv2lami1Clami2*g1p2)/5._dp + (12*DyTlami1CTlami2*g1p2)/5._dp - 12*Dylami1adjYvmlHd2i2*g2p2 +& 
&  24*AbsM2*Dylami1Clami2*g2p2 + 6*Dylami1Cmv2Clami2*g2p2 + 6*Dymv2lami1Clami2*g2p2 +    & 
&  12*DyTlami1CTlami2*g2p2 - 8._dp*(kap1adjYvYvadjkap1mv2) - 16._dp*(kap1adjYvYvCmv2adjkap1) -& 
&  16._dp*(kap1Ckap1Tpkap1Cmv2adjkap1) - 16._dp*(kap1Ckap2Tpkap2Cmv2adjkap1) -           & 
&  16._dp*(kap1Ckap3Tpkap3Cmv2adjkap1) - 16._dp*(kap1Cmv2adjYvYvadjkap1) -               & 
&  8._dp*(kap2adjYvYvadjkap2mv2) - 16._dp*(kap2adjYvYvCmv2adjkap2) - 16._dp*(kap2Ckap1Tpkap1Cmv2adjkap2)  
betamv22 =  betamv22- 16._dp*(kap2Ckap2Tpkap2Cmv2adjkap2) - 16._dp*(kap2Ckap3Tpkap3Cmv2adjkap2) -           & 
&  16._dp*(kap2Cmv2adjYvYvadjkap2) - 8._dp*(kap3adjYvYvadjkap3mv2) - 16._dp*(kap3adjYvYvCmv2adjkap3) -& 
&  16._dp*(kap3Ckap1Tpkap1Cmv2adjkap3) - 16._dp*(kap3Ckap2Tpkap2Cmv2adjkap3) -           & 
&  16._dp*(kap3Ckap3Tpkap3Cmv2adjkap3) - 16._dp*(kap3Cmv2adjYvYvadjkap3) -               & 
&  (12*Dylami1CTlami2*g1p2*M1)/5._dp - 12*Dylami1CTlami2*g2p2*M2 - 4*Dylami1adjYvYvClami2*mHd2 -& 
&  4*DyTpYvCYvlami1Clami2*mHd2 + (12*Dylami1Clami2*g1p2*mHd2)/5._dp + 12*Dylami1Clami2*g2p2*mHd2 -& 
&  16*Dykap1i11Ckap1TpYvCYvi21*mHu2 - 16*Dykap1i12Ckap2TpYvCYvi21*mHu2 - 16*Dykap1i13Ckap3TpYvCYvi21*mHu2 -& 
&  16*Dykap2i11Ckap1TpYvCYvi22*mHu2 - 16*Dykap2i12Ckap2TpYvCYvi22*mHu2 - 16*Dykap2i13Ckap3TpYvCYvi22*mHu2 -& 
&  16*Dykap3i11Ckap1TpYvCYvi23*mHu2 - 16*Dykap3i12Ckap2TpYvCYvi23*mHu2 - 16*Dykap3i13Ckap3TpYvCYvi23*mHu2 -& 
&  8*Dylami1adjYvYvClami2*mHu2 - 8*DyTpYvCYvlami1Clami2*mHu2 + (12*Dylami1Clami2*g1p2*mHu2)/5._dp +& 
&  12*Dylami1Clami2*g2p2*mHu2 - 8._dp*(mv2kap1adjYvYvadjkap1) - 8._dp*(mv2kap2adjYvYvadjkap2) -& 
&  8._dp*(mv2kap3adjYvYvadjkap3) - 2._dp*(mv2TpYvadjYeYeCYv) + (6*g1p2*mv2TpYvCYv)/5._dp +& 
&  6*g2p2*mv2TpYvCYv - 2._dp*(mv2TpYvCYvTpYvCYv) - 8*Dylami1CTlami2*SPClamxxTlam +       & 
&  8*Dylami1Clami2*SPlamxxadjYvmlHd2 + 8*Dylami1adjYvmlHd2i2*SPlamxxClam -               & 
&  4*Dylami1Cmv2Clami2*SPlamxxClam - 4*Dymv2lami1Clami2*SPlamxxClam - 8*DyTlami1CTlami2*SPlamxxClam -& 
&  16*Dylami1Clami2*mHd2*SPlamxxClam - 16*Dylami1Clami2*mHu2*SPlamxxClam -               & 
&  2*mv2TpYvCYv*SPlamxxClam - 8*Dylami1Clami2*SPlamxxCmv2Clam - 8*DyTlami1Clami2*SPlamxxCTlam -& 
&  8*Dylami1Clami2*SPTlamxxCTlam - 16._dp*(Tk1adjYvYvadjTk1) - 16._dp*(Tk2adjYvYvadjTk2) -& 
&  16._dp*(Tk3adjYvYvadjTk3) - 4._dp*(TpTvadjTeYeCYv) - 4._dp*(TpTvadjYeYeCTv) +         & 
&  (12*g1p2*TpTvCTv)/5._dp + 12*g2p2*TpTvCTv - 4*SPlamxxClam*TpTvCTv - 4._dp*(TpTvCTvTpYvCYv) -& 
&  4*SPlamxxCTlam*TpTvCYv - 4._dp*(TpTvCYvTpYvCTv) - 4._dp*(TpYvadjTeTeCYv)  
betamv22 =  betamv22- 4._dp*(TpYvadjYeme2YeCYv) - 4._dp*(TpYvadjYeTeCTv) - 4*mHd2*TpYvadjYeYeCYv -          & 
&  4*mHu2*TpYvadjYeYeCYv - 2._dp*(TpYvadjYeYeCYvmv2) - 4._dp*(TpYvadjYeYeml2CYv) -       & 
&  (12*g1p2*M1*TpYvCTv)/5._dp - 12*g2p2*M2*TpYvCTv - 4*SPClamxxTlam*TpYvCTv -            & 
&  4._dp*(TpYvCTvTpTvCYv) + (24*AbsM1*g1p2*TpYvCYv)/5._dp + 24*AbsM2*g2p2*TpYvCYv +      & 
&  (12*g1p2*mHu2*TpYvCYv)/5._dp + 12*g2p2*mHu2*TpYvCYv + 4*SPlamxxadjYvmlHd2*TpYvCYv -   & 
&  4*mHd2*SPlamxxClam*TpYvCYv - 8*mHu2*SPlamxxClam*TpYvCYv - 4*SPlamxxCmv2Clam*TpYvCYv - & 
&  4*SPTlamxxCTlam*TpYvCYv + (6*g1p2*TpYvCYvmv2)/5._dp + 6*g2p2*TpYvCYvmv2 -             & 
&  2*SPlamxxClam*TpYvCYvmv2 - 4._dp*(TpYvCYvmv2TpYvCYv) - 4._dp*(TpYvCYvTpTvCTv) -       & 
&  8*mHu2*TpYvCYvTpYvCYv - 2._dp*(TpYvCYvTpYvCYvmv2) - 4._dp*(TpYvCYvTpYvml2CYv) -       & 
&  4._dp*(TpYvml2adjYeYeCYv) + (12*g1p2*TpYvml2CYv)/5._dp + 12*g2p2*TpYvml2CYv -         & 
&  4*SPlamxxClam*TpYvml2CYv - 4._dp*(TpYvml2CYvTpYvCYv) - 16*kap1adjkap1*TradjTk1Tk1 -   & 
&  16*kap2adjkap2*TradjTk1Tk1 - 16*kap3adjkap3*TradjTk1Tk1 - 16*kap1adjkap1*TradjTk2Tk2 -& 
&  16*kap2adjkap2*TradjTk2Tk2 - 16*kap3adjkap3*TradjTk2Tk2 - 16*kap1adjkap1*TradjTk3Tk3 -& 
&  16*kap2adjkap2*TradjTk3Tk3 - 16*kap3adjkap3*TradjTk3Tk3 - 12*Dylami1CTlami2*TradjYdTd -& 
&  4*Dylami1CTlami2*TradjYeTe - 12*Dylami1CTlami2*TradjYuTu - 12*TpYvCTv*TradjYuTu -     & 
&  4*Dylami1CTlami2*TradjYvTv - 16*kap1adjTk1*TradjYvTv - 16*kap2adjTk2*TradjYvTv -      & 
&  16*kap3adjTk3*TradjYvTv - 4*TpYvCTv*TradjYvTv - 16*kap1adjTk1*TrCkap1Tk1 -            & 
&  16*kap1adjTk2*TrCkap1Tk2 - 16*kap1adjTk3*TrCkap1Tk3 - 16*DyTk1i11CTk1i21*TrCkap1Tpkap1 -& 
&  16*DyTk2i11CTk2i21*TrCkap1Tpkap1 - 16*DyTk3i11CTk3i21*TrCkap1Tpkap1 - 16*DyTk1i11CTk1i22*TrCkap1Tpkap2 -& 
&  16*DyTk2i11CTk2i22*TrCkap1Tpkap2 - 16*DyTk3i11CTk3i22*TrCkap1Tpkap2 - 16*DyTk1i11CTk1i23*TrCkap1Tpkap3 -& 
&  16*DyTk2i11CTk2i23*TrCkap1Tpkap3 - 16*DyTk3i11CTk3i23*TrCkap1Tpkap3 - 16*kap2adjTk1*TrCkap2Tk1  
betamv22 =  betamv22- 16*kap2adjTk2*TrCkap2Tk2 - 16*kap2adjTk3*TrCkap2Tk3 - 16*DyTk1i12CTk1i21*TrCkap2Tpkap1 -& 
&  16*DyTk2i12CTk2i21*TrCkap2Tpkap1 - 16*DyTk3i12CTk3i21*TrCkap2Tpkap1 - 16*DyTk1i12CTk1i22*TrCkap2Tpkap2 -& 
&  16*DyTk2i12CTk2i22*TrCkap2Tpkap2 - 16*DyTk3i12CTk3i22*TrCkap2Tpkap2 - 16*DyTk1i12CTk1i23*TrCkap2Tpkap3 -& 
&  16*DyTk2i12CTk2i23*TrCkap2Tpkap3 - 16*DyTk3i12CTk3i23*TrCkap2Tpkap3 - 16*kap3adjTk1*TrCkap3Tk1 -& 
&  16*kap3adjTk2*TrCkap3Tk2 - 16*kap3adjTk3*TrCkap3Tk3 - 16*DyTk1i13CTk1i21*TrCkap3Tpkap1 -& 
&  16*DyTk2i13CTk2i21*TrCkap3Tpkap1 - 16*DyTk3i13CTk3i21*TrCkap3Tpkap1 - 16*DyTk1i13CTk1i22*TrCkap3Tpkap2 -& 
&  16*DyTk2i13CTk2i22*TrCkap3Tpkap2 - 16*DyTk3i13CTk3i22*TrCkap3Tpkap2 - 16*DyTk1i13CTk1i23*TrCkap3Tpkap3 -& 
&  16*DyTk2i13CTk2i23*TrCkap3Tpkap3 - 16*DyTk3i13CTk3i23*TrCkap3Tpkap3 - 12*Dylami1Clami2*TrCTdTpTd -& 
&  12*DyTlami1Clami2*TrCTdTpYd - 4*Dylami1Clami2*TrCTeTpTe - 4*DyTlami1Clami2*TrCTeTpYe -& 
&  16*DyTk1i11Ckap1i21*TrCTk1Tpkap1 - 16*DyTk1i12Ckap2i21*TrCTk1Tpkap1 - 16*DyTk1i13Ckap3i21*TrCTk1Tpkap1 -& 
&  16*DyTk1i11Ckap1i22*TrCTk1Tpkap2 - 16*DyTk1i12Ckap2i22*TrCTk1Tpkap2 - 16*DyTk1i13Ckap3i22*TrCTk1Tpkap2 -& 
&  16*DyTk1i11Ckap1i23*TrCTk1Tpkap3 - 16*DyTk1i12Ckap2i23*TrCTk1Tpkap3 - 16*DyTk1i13Ckap3i23*TrCTk1Tpkap3 -& 
&  16*DyTk2i11Ckap1i21*TrCTk2Tpkap1 - 16*DyTk2i12Ckap2i21*TrCTk2Tpkap1 - 16*DyTk2i13Ckap3i21*TrCTk2Tpkap1 -& 
&  16*DyTk2i11Ckap1i22*TrCTk2Tpkap2 - 16*DyTk2i12Ckap2i22*TrCTk2Tpkap2 - 16*DyTk2i13Ckap3i22*TrCTk2Tpkap2 -& 
&  16*DyTk2i11Ckap1i23*TrCTk2Tpkap3 - 16*DyTk2i12Ckap2i23*TrCTk2Tpkap3 - 16*DyTk2i13Ckap3i23*TrCTk2Tpkap3 -& 
&  16*DyTk3i11Ckap1i21*TrCTk3Tpkap1 - 16*DyTk3i12Ckap2i21*TrCTk3Tpkap1 - 16*DyTk3i13Ckap3i21*TrCTk3Tpkap1 -& 
&  16*DyTk3i11Ckap1i22*TrCTk3Tpkap2 - 16*DyTk3i12Ckap2i22*TrCTk3Tpkap2 - 16*DyTk3i13Ckap3i22*TrCTk3Tpkap2 -& 
&  16*DyTk3i11Ckap1i23*TrCTk3Tpkap3 - 16*DyTk3i12Ckap2i23*TrCTk3Tpkap3 - 16*DyTk3i13Ckap3i23*TrCTk3Tpkap3 -& 
&  12*Dylami1Clami2*TrCTuTpTu - 12*TpYvCYv*TrCTuTpTu - 12*DyTlami1Clami2*TrCTuTpYu -     & 
&  12*TpTvCYv*TrCTuTpYu - 4*Dylami1Clami2*TrCTvTpTv - 16*kap1adjkap1*TrCTvTpTv -         & 
&  16*kap2adjkap2*TrCTvTpTv - 16*kap3adjkap3*TrCTvTpTv - 4*TpYvCYv*TrCTvTpTv  
betamv22 =  betamv22- 4*DyTlami1Clami2*TrCTvTpYv - 4*TpTvCYv*TrCTvTpYv - 12*Dylami1Clami2*Trmd2YdadjYd -    & 
&  4*Dylami1Clami2*Trme2YeadjYe - 4*Dylami1Clami2*Trml2adjYeYe - 12*Dylami1Clami2*Trmq2adjYdYd -& 
&  12*Dylami1Clami2*Trmq2adjYuYu - 12*TpYvCYv*Trmq2adjYuYu - 12*Dylami1Clami2*Trmu2YuadjYu -& 
&  12*TpYvCYv*Trmu2YuadjYu - 32*kap1adjkap1*Trmv2kap1adjkap1 - 32*kap2adjkap1*Trmv2kap1adjkap2 -& 
&  32*kap3adjkap1*Trmv2kap1adjkap3 - 32*kap1adjkap2*Trmv2kap2adjkap1 - 32*kap2adjkap2*Trmv2kap2adjkap2 -& 
&  32*kap3adjkap2*Trmv2kap2adjkap3 - 32*kap1adjkap3*Trmv2kap3adjkap1 - 32*kap2adjkap3*Trmv2kap3adjkap2 -& 
&  32*kap3adjkap3*Trmv2kap3adjkap3 + 12*Dylami1adjYvmlHd2i2*TrYdadjYd - 6*Dylami1Cmv2Clami2*TrYdadjYd -& 
&  6*Dymv2lami1Clami2*TrYdadjYd - 12*DyTlami1CTlami2*TrYdadjYd - 24*Dylami1Clami2*mHd2*TrYdadjYd -& 
&  12*Dylami1Clami2*mHu2*TrYdadjYd + 4*Dylami1adjYvmlHd2i2*TrYeadjYe - 2*Dylami1Cmv2Clami2*TrYeadjYe -& 
&  2*Dymv2lami1Clami2*TrYeadjYe - 4*DyTlami1CTlami2*TrYeadjYe - 8*Dylami1Clami2*mHd2*TrYeadjYe -& 
&  4*Dylami1Clami2*mHu2*TrYeadjYe + 12*Dylami1adjYvmlHd2i2*TrYuadjYu - 6*Dylami1Cmv2Clami2*TrYuadjYu -& 
&  6*Dymv2lami1Clami2*TrYuadjYu - 12*DyTlami1CTlami2*TrYuadjYu - 12*Dylami1Clami2*mHd2*TrYuadjYu -& 
&  24*Dylami1Clami2*mHu2*TrYuadjYu - 6*mv2TpYvCYv*TrYuadjYu - 12*TpTvCTv*TrYuadjYu -     & 
&  24*mHu2*TpYvCYv*TrYuadjYu - 6*TpYvCYvmv2*TrYuadjYu - 12*TpYvml2CYv*TrYuadjYu +        & 
&  4*Dylami1adjYvmlHd2i2*TrYvadjYv - 2*Dylami1Cmv2Clami2*TrYvadjYv - 2*Dymv2lami1Clami2*TrYvadjYv -& 
&  4*DyTlami1CTlami2*TrYvadjYv - 4*Dylami1Clami2*mHd2*TrYvadjYv - 8*Dylami1Clami2*mHu2*TrYvadjYv -& 
&  2*mv2TpYvCYv*TrYvadjYv - 4*TpTvCTv*TrYvadjYv - 8*mHu2*TpYvCYv*TrYvadjYv -             & 
&  2*TpYvCYvmv2*TrYvadjYv - 4*TpYvml2CYv*TrYvadjYv - 4*Dylami1Clami2*TrYvadjYvCml2 -     & 
&  4*TpYvCYv*TrYvadjYvCml2 - 4*Dylami1Clami2*TrYvCmv2adjYv - 4*TpYvCYv*TrYvCmv2adjYv -   & 
&  32*kap2adjkap1*adjkap1mv2kap1(2,1) - 32*kap3adjkap1*adjkap1mv2kap1(3,1) -             & 
&  32*kap2adjkap2*adjkap1mv2kap2(2,1) - 32*kap3adjkap2*adjkap1mv2kap2(3,1)  
betamv22 =  betamv22- 32*kap2adjkap3*adjkap1mv2kap3(2,1) - 32*kap3adjkap3*adjkap1mv2kap3(3,1) -             & 
&  32*kap1adjkap1*adjkap2mv2kap1(1,2) - 32*kap3adjkap1*adjkap2mv2kap1(3,2) -             & 
&  32*kap1adjkap2*adjkap2mv2kap2(1,2) - 32*kap3adjkap2*adjkap2mv2kap2(3,2) -             & 
&  32*kap1adjkap3*adjkap2mv2kap3(1,2) - 32*kap3adjkap3*adjkap2mv2kap3(3,2) -             & 
&  32*kap1adjkap1*adjkap3mv2kap1(1,3) - 32*kap2adjkap1*adjkap3mv2kap1(2,3) -             & 
&  32*kap1adjkap2*adjkap3mv2kap2(1,3) - 32*kap2adjkap2*adjkap3mv2kap2(2,3) -             & 
&  32*kap1adjkap3*adjkap3mv2kap3(1,3) - 32*kap2adjkap3*adjkap3mv2kap3(2,3) -             & 
&  16*kap1adjTk2*adjYvTv(1,2) - 16*kap1adjTk3*adjYvTv(1,3) - 16*kap2adjTk1*adjYvTv(2,    & 
& 1) - 16*kap2adjTk3*adjYvTv(2,3) - 16*kap3adjTk1*adjYvTv(3,1) - 16*kap3adjTk2*adjYvTv(3,& 
& 2) - 16*kap2adjTk1*Ckap1Tk1(2,1) - 16*kap3adjTk1*Ckap1Tk1(3,1) - 16*kap2adjTk2*Ckap1Tk2(2,& 
& 1) - 16*kap3adjTk2*Ckap1Tk2(3,1) - 16*kap2adjTk3*Ckap1Tk3(2,1) - 16*kap3adjTk3*Ckap1Tk3(3,& 
& 1) - 16*kap1adjTk1*Ckap2Tk1(1,2) - 16*kap3adjTk1*Ckap2Tk1(3,2) - 16*kap1adjTk2*Ckap2Tk2(1,& 
& 2) - 16*kap3adjTk2*Ckap2Tk2(3,2) - 16*kap1adjTk3*Ckap2Tk3(1,2) - 16*kap3adjTk3*Ckap2Tk3(3,& 
& 2) - 16*kap1adjTk1*Ckap3Tk1(1,3) - 16*kap2adjTk1*Ckap3Tk1(2,3) - 16*kap1adjTk2*Ckap3Tk2(1,& 
& 3) - 16*kap2adjTk2*Ckap3Tk2(2,3) - 16*kap1adjTk3*Ckap3Tk3(1,3) - 16*kap2adjTk3*Ckap3Tk3(2,& 
& 3) - (12*DyTlami1Clami2*g1p2*Conjg(M1))/5._dp - (12*g1p2*TpTvCYv*Conjg(M1))/5._dp -    & 
&  12*DyTlami1Clami2*g2p2*Conjg(M2) - 12*g2p2*TpTvCYv*Conjg(M2) - 16*DyTk1i11Ckap1lami2*Conjg(Tlam(1)) -& 
&  16*DyTk1i12Ckap2lami2*Conjg(Tlam(1)) - 16*DyTk1i13Ckap3lami2*Conjg(Tlam(1)) -         & 
&  16*DyTk2i11Ckap1lami2*Conjg(Tlam(2)) - 16*DyTk2i12Ckap2lami2*Conjg(Tlam(2)) -         & 
&  16*DyTk2i13Ckap3lami2*Conjg(Tlam(2)) - 16*DyTk3i11Ckap1lami2*Conjg(Tlam(3)) -         & 
&  16*DyTk3i12Ckap2lami2*Conjg(Tlam(3)) - 16*DyTk3i13Ckap3lami2*Conjg(Tlam(3))  
betamv22 =  betamv22- 16*kap2adjkap1*Tk1adjTk1(1,2) - 16*kap3adjkap1*Tk1adjTk1(1,3) - 16*kap1adjkap2*Tk1adjTk1(2,& 
& 1) - 16*kap3adjkap2*Tk1adjTk1(2,3) - 16*kap1adjkap3*Tk1adjTk1(3,1) - 16*kap2adjkap3*Tk1adjTk1(3,& 
& 2) - 16*kap2adjkap1*Tk2adjTk2(1,2) - 16*kap3adjkap1*Tk2adjTk2(1,3) - 16*kap1adjkap2*Tk2adjTk2(2,& 
& 1) - 16*kap3adjkap2*Tk2adjTk2(2,3) - 16*kap1adjkap3*Tk2adjTk2(3,1) - 16*kap2adjkap3*Tk2adjTk2(3,& 
& 2) - 16*kap2adjkap1*Tk3adjTk3(1,2) - 16*kap3adjkap1*Tk3adjTk3(1,3) - 16*kap1adjkap2*Tk3adjTk3(2,& 
& 1) - 16*kap3adjkap2*Tk3adjTk3(2,3) - 16*kap1adjkap3*Tk3adjTk3(3,1) - 16*kap2adjkap3*Tk3adjTk3(3,& 
& 2) - 16*kap2adjkap1*TpTvCTv(1,2) - 16*kap3adjkap1*TpTvCTv(1,3) - 16*kap1adjkap2*TpTvCTv(2,& 
& 1) - 16*kap3adjkap2*TpTvCTv(2,3) - 16*kap1adjkap3*TpTvCTv(3,1) - 16*kap2adjkap3*TpTvCTv(3,& 
& 2) - 16*Dykap1Clami1Ckap1mv2i21*lam(1) - 8*Dykap1Clami1Cmv2Ckap1i21*lam(1) -           & 
&  16*Dykap1Cmv2Clami1Ckap1i21*lam(1) - 16*Dykap2Clami1Ckap1mv2i22*lam(1) -              & 
&  8*Dykap2Clami1Cmv2Ckap1i22*lam(1) - 16*Dykap2Cmv2Clami1Ckap1i22*lam(1) -              & 
&  16*Dykap3Clami1Ckap1mv2i23*lam(1) - 8*Dykap3Clami1Cmv2Ckap1i23*lam(1) -               & 
&  16*Dykap3Cmv2Clami1Ckap1i23*lam(1) + 16*kap1adjkap1*adjYvmlHd2(1)*lam(1) +            & 
&  16*kap2adjkap1*adjYvmlHd2(2)*lam(1) + 16*kap3adjkap1*adjYvmlHd2(3)*lam(1) -           & 
&  16*kap1adjkap1*mHd2*Conjg(lam(1))*lam(1) - 16*kap1adjkap1*mHu2*Conjg(lam(1))*lam(1) - & 
&  16*kap2adjkap1*mHd2*Conjg(lam(2))*lam(1) - 16*kap2adjkap1*mHu2*Conjg(lam(2))*lam(1) - & 
&  16*kap3adjkap1*mHd2*Conjg(lam(3))*lam(1) - 16*kap3adjkap1*mHu2*Conjg(lam(3))*lam(1) - & 
&  16*Dykap1Clami1Ckap2mv2i21*lam(2) - 8*Dykap1Clami1Cmv2Ckap2i21*lam(2) -               & 
&  16*Dykap1Cmv2Clami1Ckap2i21*lam(2) - 16*Dykap2Clami1Ckap2mv2i22*lam(2) -              & 
&  8*Dykap2Clami1Cmv2Ckap2i22*lam(2) - 16*Dykap2Cmv2Clami1Ckap2i22*lam(2) -              & 
&  16*Dykap3Clami1Ckap2mv2i23*lam(2) - 8*Dykap3Clami1Cmv2Ckap2i23*lam(2)  
betamv22 =  betamv22- 16*Dykap3Cmv2Clami1Ckap2i23*lam(2) + 16*kap1adjkap2*adjYvmlHd2(1)*lam(2) +            & 
&  16*kap2adjkap2*adjYvmlHd2(2)*lam(2) + 16*kap3adjkap2*adjYvmlHd2(3)*lam(2) -           & 
&  16*kap1adjkap2*mHd2*Conjg(lam(1))*lam(2) - 16*kap1adjkap2*mHu2*Conjg(lam(1))*lam(2) - & 
&  16*kap2adjkap2*mHd2*Conjg(lam(2))*lam(2) - 16*kap2adjkap2*mHu2*Conjg(lam(2))*lam(2) - & 
&  16*kap3adjkap2*mHd2*Conjg(lam(3))*lam(2) - 16*kap3adjkap2*mHu2*Conjg(lam(3))*lam(2) - & 
&  16*Dykap1Clami1Ckap3mv2i21*lam(3) - 8*Dykap1Clami1Cmv2Ckap3i21*lam(3) -               & 
&  16*Dykap1Cmv2Clami1Ckap3i21*lam(3) - 16*Dykap2Clami1Ckap3mv2i22*lam(3) -              & 
&  8*Dykap2Clami1Cmv2Ckap3i22*lam(3) - 16*Dykap2Cmv2Clami1Ckap3i22*lam(3) -              & 
&  16*Dykap3Clami1Ckap3mv2i23*lam(3) - 8*Dykap3Clami1Cmv2Ckap3i23*lam(3) -               & 
&  16*Dykap3Cmv2Clami1Ckap3i23*lam(3) + 16*kap1adjkap3*adjYvmlHd2(1)*lam(3) +            & 
&  16*kap2adjkap3*adjYvmlHd2(2)*lam(3) + 16*kap3adjkap3*adjYvmlHd2(3)*lam(3) -           & 
&  16*kap1adjkap3*mHd2*Conjg(lam(1))*lam(3) - 16*kap1adjkap3*mHu2*Conjg(lam(1))*lam(3) - & 
&  16*kap2adjkap3*mHd2*Conjg(lam(2))*lam(3) - 16*kap2adjkap3*mHu2*Conjg(lam(2))*lam(3) - & 
&  16*kap3adjkap3*mHd2*Conjg(lam(3))*lam(3) - 16*kap3adjkap3*mHu2*Conjg(lam(3))*lam(3) - & 
&  16*kap1adjTk1*Conjg(lam(1))*Tlam(1) - 16*kap2adjTk1*Conjg(lam(2))*Tlam(1) -           & 
&  16*kap3adjTk1*Conjg(lam(3))*Tlam(1) - 16*kap1adjkap1*Conjg(Tlam(1))*Tlam(1) -         & 
&  16*kap2adjkap1*Conjg(Tlam(2))*Tlam(1) - 16*kap3adjkap1*Conjg(Tlam(3))*Tlam(1) -       & 
&  16*kap1adjTk2*Conjg(lam(1))*Tlam(2) - 16*kap2adjTk2*Conjg(lam(2))*Tlam(2) -           & 
&  16*kap3adjTk2*Conjg(lam(3))*Tlam(2) - 16*kap1adjkap2*Conjg(Tlam(1))*Tlam(2) -         & 
&  16*kap2adjkap2*Conjg(Tlam(2))*Tlam(2) - 16*kap3adjkap2*Conjg(Tlam(3))*Tlam(2) -       & 
&  16*kap1adjTk3*Conjg(lam(1))*Tlam(3) - 16*kap2adjTk3*Conjg(lam(2))*Tlam(3)  
betamv22 =  betamv22- 16*kap3adjTk3*Conjg(lam(3))*Tlam(3) - 16*kap1adjkap3*Conjg(Tlam(1))*Tlam(3) -         & 
&  16*kap2adjkap3*Conjg(Tlam(2))*Tlam(3) - 16*kap3adjkap3*Conjg(Tlam(3))*Tlam(3)

 
Dmv2 = oo16pi2*( betamv21 + oo16pi2 * betamv22 ) 

 
Else 
Dmv2 = oo16pi2* betamv21 
End If 
 
 
Call Chop(Dmv2) 

Forall(i1=1:3) Dmv2(i1,i1) =  Real(Dmv2(i1,i1),dp) 
Dmv2 = 0.5_dp*(Dmv2+ Conjg(Transpose(Dmv2)) ) 
!-------------------- 
! mlHd2 
!-------------------- 
 
betamlHd21  = -2._dp*(CTvTlam) - 2._dp*(CYvmv2lam) - CYvlam*(mHd2 + 2._dp*(mHu2))     & 
&  - ml2CYvlam

 
 
If (TwoLoopRGE) Then 
betamlHd22 = -2._dp*(adjTeTeCYvlam) - 2._dp*(adjTeYeCYvTlam) - 2._dp*(adjYeme2YeCYvlam) -          & 
&  2._dp*(adjYeTeCTvlam) - 2._dp*(adjYeYeCTvTlam) - 2._dp*(adjYeYeCYvmv2lam) -           & 
&  2._dp*(adjYeYeml2CYvlam) + 4._dp*(CTvTpTvCYvlam) + 4._dp*(CTvTpYvCYvTlam) +           & 
&  4._dp*(CYvmv2TpYvCYvlam) + 4._dp*(CYvTk1adjTk1lam) + 4._dp*(CYvTk2adjTk2lam) +        & 
&  4._dp*(CYvTk3adjTk3lam) + 4._dp*(CYvTpTvCTvlam) + 4._dp*(CYvTpYvCTvTlam) +            & 
&  4._dp*(CYvTpYvCYvmv2lam) + 4._dp*(CYvTpYvml2CYvlam) - 3*adjYeYeCYvlam*mHd2 +          & 
&  2*CYvTpYvCYvlam*mHd2 - 2*adjYeYeCYvlam*mHu2 + 8*CYvTpYvCYvlam*mHu2 - ml2adjYeYeCYvlam +& 
&  2._dp*(ml2CYvTpYvCYvlam) + 6*CTvlam*SPClamxxTlam - 6*CYvlam*SPlamxxadjYvmlHd2 +       & 
&  6*CTvTlam*SPlamxxClam + 6*CYvmv2lam*SPlamxxClam + 9*CYvlam*mHd2*SPlamxxClam +         & 
&  12*CYvlam*mHu2*SPlamxxClam + 3*ml2CYvlam*SPlamxxClam + 6*CYvlam*SPlamxxCmv2Clam +     & 
&  6*CYvTlam*SPlamxxCTlam + 6*CYvlam*SPTlamxxCTlam + 6*CTvlam*TradjYuTu + 2*CTvlam*TradjYvTv +& 
&  6*CYvlam*TrCTuTpTu + 6*CYvTlam*TrCTuTpYu + 2*CYvlam*TrCTvTpTv + 2*CYvTlam*TrCTvTpYv + & 
&  6*CYvlam*Trmq2adjYuYu + 6*CYvlam*Trmu2YuadjYu + 6*CTvTlam*TrYuadjYu + 6*CYvmv2lam*TrYuadjYu +& 
&  3*CYvlam*mHd2*TrYuadjYu + 12*CYvlam*mHu2*TrYuadjYu + 3*ml2CYvlam*TrYuadjYu +          & 
&  2*CTvTlam*TrYvadjYv + 2*CYvmv2lam*TrYvadjYv + CYvlam*mHd2*TrYvadjYv + 4*CYvlam*mHu2*TrYvadjYv +& 
&  ml2CYvlam*TrYvadjYv + 2*CYvlam*TrYvadjYvCml2 + 2*CYvlam*TrYvCmv2adjYv +               & 
&  4*Conjg(Yv(i1,1))*Tpkap1adjkap1mv2lam(1) + 4*Conjg(Yv(i1,1))*Tpkap1adjkap2mv2lam(2) + & 
&  4*Conjg(Yv(i1,1))*Tpkap1adjkap3mv2lam(3) + 4*mHu2*Conjg(Yv(i1,1))*Tpkap1Ckap1lam(1) + & 
&  2*ml2CYv(i1,1)*Tpkap1Ckap1lam(1) + 4*mHu2*Conjg(Yv(i1,1))*Tpkap1Ckap2lam(2) +         & 
&  2*ml2CYv(i1,1)*Tpkap1Ckap2lam(2) + 4*mHu2*Conjg(Yv(i1,1))*Tpkap1Ckap3lam(3) +         & 
&  2*ml2CYv(i1,1)*Tpkap1Ckap3lam(3) + 4*Conjg(Yv(i1,2))*Tpkap2adjkap1mv2lam(1)  
betamlHd22 =  betamlHd22+ 4*Conjg(Yv(i1,2))*Tpkap2adjkap2mv2lam(2) + 4*Conjg(Yv(i1,2))*Tpkap2adjkap3mv2lam(3) + & 
&  4*mHu2*Conjg(Yv(i1,2))*Tpkap2Ckap1lam(1) + 2*ml2CYv(i1,2)*Tpkap2Ckap1lam(1) +         & 
&  4*mHu2*Conjg(Yv(i1,2))*Tpkap2Ckap2lam(2) + 2*ml2CYv(i1,2)*Tpkap2Ckap2lam(2) +         & 
&  4*mHu2*Conjg(Yv(i1,2))*Tpkap2Ckap3lam(3) + 2*ml2CYv(i1,2)*Tpkap2Ckap3lam(3) +         & 
&  4*Conjg(Yv(i1,3))*Tpkap3adjkap1mv2lam(1) + 4*Conjg(Yv(i1,3))*Tpkap3adjkap2mv2lam(2) + & 
&  4*Conjg(Yv(i1,3))*Tpkap3adjkap3mv2lam(3) + 4*mHu2*Conjg(Yv(i1,3))*Tpkap3Ckap1lam(1) + & 
&  2*ml2CYv(i1,3)*Tpkap3Ckap1lam(1) + 4*mHu2*Conjg(Yv(i1,3))*Tpkap3Ckap2lam(2) +         & 
&  2*ml2CYv(i1,3)*Tpkap3Ckap2lam(2) + 4*mHu2*Conjg(Yv(i1,3))*Tpkap3Ckap3lam(3) +         & 
&  2*ml2CYv(i1,3)*Tpkap3Ckap3lam(3) + 4*Conjg(Tv(i1,1))*TpTk1adjkap1lam(1) +             & 
&  4*Conjg(Tv(i1,1))*TpTk1adjkap2lam(2) + 4*Conjg(Tv(i1,1))*TpTk1adjkap3lam(3) +         & 
&  4*Conjg(Tv(i1,2))*TpTk2adjkap1lam(1) + 4*Conjg(Tv(i1,2))*TpTk2adjkap2lam(2) +         & 
&  4*Conjg(Tv(i1,2))*TpTk2adjkap3lam(3) + 4*Conjg(Tv(i1,3))*TpTk3adjkap1lam(1) +         & 
&  4*Conjg(Tv(i1,3))*TpTk3adjkap2lam(2) + 4*Conjg(Tv(i1,3))*TpTk3adjkap3lam(3) +         & 
&  8*Trmv2kap1adjkap1*Conjg(Yv(i1,1))*lam(1) + 8*Trmv2kap2adjkap1*Conjg(Yv(i1,           & 
& 2))*lam(1) + 8*Trmv2kap3adjkap1*Conjg(Yv(i1,3))*lam(1) + 4*CYvmv2kap1Ckap1(i1,         & 
& 1)*lam(1) + 4*CYvmv2kap2Ckap1(i1,2)*lam(1) + 4*CYvmv2kap3Ckap1(i1,3)*lam(1) +          & 
&  2*mHd2*CYvTpkap1Ckap1(i1,1)*lam(1) + 2*mHd2*CYvTpkap2Ckap1(i1,2)*lam(1) +             & 
&  2*mHd2*CYvTpkap3Ckap1(i1,3)*lam(1) + 8*Trmv2kap1adjkap2*Conjg(Yv(i1,1))*lam(2) +      & 
&  8*Trmv2kap2adjkap2*Conjg(Yv(i1,2))*lam(2) + 8*Trmv2kap3adjkap2*Conjg(Yv(i1,           & 
& 3))*lam(2) + 4*CYvmv2kap1Ckap2(i1,1)*lam(2) + 4*CYvmv2kap2Ckap2(i1,2)*lam(2) +         & 
&  4*CYvmv2kap3Ckap2(i1,3)*lam(2) + 2*mHd2*CYvTpkap1Ckap2(i1,1)*lam(2) + 2*mHd2*CYvTpkap2Ckap2(i1,& 
& 2)*lam(2) + 2*mHd2*CYvTpkap3Ckap2(i1,3)*lam(2) + 8*Trmv2kap1adjkap3*Conjg(Yv(i1,       & 
& 1))*lam(3) + 8*Trmv2kap2adjkap3*Conjg(Yv(i1,2))*lam(3) + 8*Trmv2kap3adjkap3*Conjg(Yv(i1,& 
& 3))*lam(3) + 4*CYvmv2kap1Ckap3(i1,1)*lam(3) + 4*CYvmv2kap2Ckap3(i1,2)*lam(3)  
betamlHd22 =  betamlHd22+ 4*CYvmv2kap3Ckap3(i1,3)*lam(3) + 2*mHd2*CYvTpkap1Ckap3(i1,1)*lam(3) + 2*mHd2*CYvTpkap2Ckap3(i1,& 
& 2)*lam(3) + 2*mHd2*CYvTpkap3Ckap3(i1,3)*lam(3) + 4*TrCTk1Tpkap1*Conjg(Yv(i1,           & 
& 1))*Tlam(1) + 4*TrCTk1Tpkap2*Conjg(Yv(i1,2))*Tlam(1) + 4*TrCTk1Tpkap3*Conjg(Yv(i1,     & 
& 3))*Tlam(1) + 4*TrCkap1Tpkap1*Conjg(Tv(i1,1))*Tlam(1) + 4*TrCkap1Tpkap2*Conjg(Tv(i1,   & 
& 2))*Tlam(1) + 4*TrCkap1Tpkap3*Conjg(Tv(i1,3))*Tlam(1) + 4*TrCTk2Tpkap1*Conjg(Yv(i1,    & 
& 1))*Tlam(2) + 4*TrCTk2Tpkap2*Conjg(Yv(i1,2))*Tlam(2) + 4*TrCTk2Tpkap3*Conjg(Yv(i1,     & 
& 3))*Tlam(2) + 4*TrCkap2Tpkap1*Conjg(Tv(i1,1))*Tlam(2) + 4*TrCkap2Tpkap2*Conjg(Tv(i1,   & 
& 2))*Tlam(2) + 4*TrCkap2Tpkap3*Conjg(Tv(i1,3))*Tlam(2) + 4*TrCTk3Tpkap1*Conjg(Yv(i1,    & 
& 1))*Tlam(3) + 4*TrCTk3Tpkap2*Conjg(Yv(i1,2))*Tlam(3) + 4*TrCTk3Tpkap3*Conjg(Yv(i1,     & 
& 3))*Tlam(3) + 4*TrCkap3Tpkap1*Conjg(Tv(i1,1))*Tlam(3) + 4*TrCkap3Tpkap2*Conjg(Tv(i1,   & 
& 2))*Tlam(3) + 4*TrCkap3Tpkap3*Conjg(Tv(i1,3))*Tlam(3)

 
DmlHd2 = oo16pi2*( betamlHd21 + oo16pi2 * betamlHd22 ) 

 
Else 
DmlHd2 = oo16pi2* betamlHd21 
End If 
 
 
!-------------------- 
! M1 
!-------------------- 
 
betaM11  = (66*g1p2*M1)/5._dp

 
 
If (TwoLoopRGE) Then 
betaM12 = (2*(398*g1p4*M1 + 5*g1p2*(27*g2p2*(M1 + M2) + 2*(44*g3p2*(M1 + M3) + 3._dp*(SPClamxxTlam) -& 
&  3*M1*SPlamxxClam + 7._dp*(TradjYdTd) + 9._dp*(TradjYeTe) + 13._dp*(TradjYuTu) +       & 
&  3._dp*(TradjYvTv) - 7*M1*TrYdadjYd - 9*M1*TrYeadjYe - 13*M1*TrYuadjYu -               & 
&  3*M1*TrYvadjYv))))/25._dp

 
DM1 = oo16pi2*( betaM11 + oo16pi2 * betaM12 ) 

 
Else 
DM1 = oo16pi2* betaM11 
End If 
 
 
Call Chop(DM1) 

!-------------------- 
! M2 
!-------------------- 
 
betaM21  = 2*g2p2*M2

 
 
If (TwoLoopRGE) Then 
betaM22 = (18*g1p2*g2p2*(M1 + M2))/5._dp + 4*(25*g2p4*M2 + g2p2*(12*g3p2*(M3 + M2) +            & 
&  SPClamxxTlam - M2*SPlamxxClam + 3._dp*(TradjYdTd) + TradjYeTe + 3._dp*(TradjYuTu) +   & 
&  TradjYvTv - 3*M2*TrYdadjYd - M2*TrYeadjYe - 3*M2*TrYuadjYu - M2*TrYvadjYv))

 
DM2 = oo16pi2*( betaM21 + oo16pi2 * betaM22 ) 

 
Else 
DM2 = oo16pi2* betaM21 
End If 
 
 
Call Chop(DM2) 

!-------------------- 
! M3 
!-------------------- 
 
betaM31  = -6*g3p2*M3

 
 
If (TwoLoopRGE) Then 
betaM32 = (22*g1p2*g3p2*(M1 + M3))/5._dp + 2*(9*g2p2*g3p2*(M3 + M2) + 4*(7*g3p4*M3 +            & 
&  g3p2*(TradjYdTd + TradjYuTu - M3*(TrYdadjYd + TrYuadjYu))))

 
DM3 = oo16pi2*( betaM31 + oo16pi2 * betaM32 ) 

 
Else 
DM3 = oo16pi2* betaM31 
End If 
 
 
Call Chop(DM3) 

Call ParametersToG386(Dg1,Dg2,Dg3,DYd,DYe,Dlam,DYv,DYu,Dkap,DTd,DTe,DTlam,            & 
& DTv,DTu,DTk,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,Dmv2,DmlHd2,DM1,DM2,DM3,f)

Iname = Iname - 1 
 
End Subroutine rge386  

Subroutine GToParameters394(g,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,             & 
& Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)

Implicit None 
Real(dp), Intent(in) :: g(394) 
Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp),Intent(out) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Integer i1, i2, i3, i4, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'GToParameters394' 
 
g1= g(1) 
g2= g(2) 
g3= g(3) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yd(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Ye(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
lam(i1) = Cmplx( g(SumI+40), g(SumI+41), dp) 
End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yv(i1,i2) = Cmplx( g(SumI+46), g(SumI+47), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yu(i1,i2) = Cmplx( g(SumI+64), g(SumI+65), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
kap(i1,i2,i3) = Cmplx( g(SumI+82), g(SumI+83), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Td(i1,i2) = Cmplx( g(SumI+136), g(SumI+137), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Te(i1,i2) = Cmplx( g(SumI+154), g(SumI+155), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
Tlam(i1) = Cmplx( g(SumI+172), g(SumI+173), dp) 
End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Tv(i1,i2) = Cmplx( g(SumI+178), g(SumI+179), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Tu(i1,i2) = Cmplx( g(SumI+196), g(SumI+197), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
Tk(i1,i2,i3) = Cmplx( g(SumI+214), g(SumI+215), dp) 
End Do 
 End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mq2(i1,i2) = Cmplx( g(SumI+268), g(SumI+269), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
ml2(i1,i2) = Cmplx( g(SumI+286), g(SumI+287), dp) 
End Do 
 End Do 
 
mHd2= g(304) 
mHu2= g(305) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
md2(i1,i2) = Cmplx( g(SumI+306), g(SumI+307), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mu2(i1,i2) = Cmplx( g(SumI+324), g(SumI+325), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
me2(i1,i2) = Cmplx( g(SumI+342), g(SumI+343), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mv2(i1,i2) = Cmplx( g(SumI+360), g(SumI+361), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
mlHd2(i1) =  g(SumI+378) 
End Do 
 
M1= Cmplx(g(381),g(382),dp) 
M2= Cmplx(g(383),g(384),dp) 
M3= Cmplx(g(385),g(386),dp) 
vd= g(387) 
vu= g(388) 
Do i1 = 1,3
SumI = (i1-1) 
vL(i1) =  g(SumI+389) 
End Do 
 
Do i1 = 1,3
SumI = (i1-1) 
vR(i1) =  g(SumI+392) 
End Do 
 
Do i1=1,394 
If (g(i1).ne.g(i1)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Write(*,*) "At position ", i1 
 Call TerminateProgram 
    return
End if 
End do 
Iname = Iname - 1 
 
End Subroutine GToParameters394

Subroutine ParametersToG394(g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,               & 
& Tu,Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR,g)

Implicit None 
Real(dp), Intent(out) :: g(394) 
Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2,mlHd2(3),vd,vu,vL(3),vR(3)

Complex(dp), Intent(in) :: Yd(3,3),Ye(3,3),lam(3),Yv(3,3),Yu(3,3),kap(3,3,3),Td(3,3),Te(3,3),Tlam(3),            & 
& Tv(3,3),Tu(3,3),Tk(3,3,3),mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),mv2(3,3),M1,M2,M3

Integer i1, i2, i3, i4, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'ParametersToG394' 
 
g(1) = g1  
g(2) = g2  
g(3) = g3  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+4) = Real(Yd(i1,i2), dp) 
g(SumI+5) = Aimag(Yd(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+22) = Real(Ye(i1,i2), dp) 
g(SumI+23) = Aimag(Ye(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
g(SumI+40) = Real(lam(i1), dp) 
g(SumI+41) = Aimag(lam(i1)) 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+46) = Real(Yv(i1,i2), dp) 
g(SumI+47) = Aimag(Yv(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+64) = Real(Yu(i1,i2), dp) 
g(SumI+65) = Aimag(Yu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
g(SumI+82) = Real(kap(i1,i2,i3), dp) 
g(SumI+83) = Aimag(kap(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+136) = Real(Td(i1,i2), dp) 
g(SumI+137) = Aimag(Td(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+154) = Real(Te(i1,i2), dp) 
g(SumI+155) = Aimag(Te(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
SumI = SumI*2 
g(SumI+172) = Real(Tlam(i1), dp) 
g(SumI+173) = Aimag(Tlam(i1)) 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+178) = Real(Tv(i1,i2), dp) 
g(SumI+179) = Aimag(Tv(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+196) = Real(Tu(i1,i2), dp) 
g(SumI+197) = Aimag(Tu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
Do i3 = 1,3
SumI = (i3-1) + (i2-1)*3 + (i1-1)*9
SumI = SumI*2 
g(SumI+214) = Real(Tk(i1,i2,i3), dp) 
g(SumI+215) = Aimag(Tk(i1,i2,i3)) 
End Do 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+268) = Real(mq2(i1,i2), dp) 
g(SumI+269) = Aimag(mq2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+286) = Real(ml2(i1,i2), dp) 
g(SumI+287) = Aimag(ml2(i1,i2)) 
End Do 
End Do 

g(304) = mHd2  
g(305) = mHu2  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+306) = Real(md2(i1,i2), dp) 
g(SumI+307) = Aimag(md2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+324) = Real(mu2(i1,i2), dp) 
g(SumI+325) = Aimag(mu2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+342) = Real(me2(i1,i2), dp) 
g(SumI+343) = Aimag(me2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+360) = Real(mv2(i1,i2), dp) 
g(SumI+361) = Aimag(mv2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
g(SumI+378) = mlHd2(i1) 
End Do 

g(381) = Real(M1,dp)  
g(382) = Aimag(M1)  
g(383) = Real(M2,dp)  
g(384) = Aimag(M2)  
g(385) = Real(M3,dp)  
g(386) = Aimag(M3)  
g(387) = vd  
g(388) = vu  
Do i1 = 1,3
SumI = (i1-1) 
g(SumI+389) = vL(i1) 
End Do 

Do i1 = 1,3
SumI = (i1-1) 
g(SumI+392) = vR(i1) 
End Do 

Iname = Iname - 1 
 
End Subroutine ParametersToG394

Subroutine rge394(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2,i3,i4 
Integer :: j1,j2,j3,j4,j5,j6,j7 
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,         & 
& Dg3,mHd2,betamHd21,betamHd22,DmHd2,mHu2,betamHu21,betamHu22,DmHu2,mlHd2(3)             & 
& ,betamlHd21(3),betamlHd22(3),DmlHd2(3),vd,betavd1,betavd2,Dvd,vu,betavu1,              & 
& betavu2,Dvu,vL(3),betavL1(3),betavL2(3),DvL(3),vR(3),betavR1(3),betavR2(3),DvR(3)
Complex(dp) :: Yd(3,3),betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3)          & 
& ,betaYe1(3,3),betaYe2(3,3),DYe(3,3),adjYe(3,3),lam(3),betalam1(3),betalam2(3)          & 
& ,Dlam(3),Yv(3,3),betaYv1(3,3),betaYv2(3,3),DYv(3,3),adjYv(3,3),Yu(3,3),betaYu1(3,3)    & 
& ,betaYu2(3,3),DYu(3,3),adjYu(3,3),kap(3,3,3),betakap1(3,3,3),betakap2(3,3,3)           & 
& ,Dkap(3,3,3),Td(3,3),betaTd1(3,3),betaTd2(3,3),DTd(3,3),adjTd(3,3),Te(3,3)             & 
& ,betaTe1(3,3),betaTe2(3,3),DTe(3,3),adjTe(3,3),Tlam(3),betaTlam1(3),betaTlam2(3)       & 
& ,DTlam(3),Tv(3,3),betaTv1(3,3),betaTv2(3,3),DTv(3,3),adjTv(3,3),Tu(3,3),               & 
& betaTu1(3,3),betaTu2(3,3),DTu(3,3),adjTu(3,3),Tk(3,3,3),betaTk1(3,3,3),betaTk2(3,3,3)  & 
& ,DTk(3,3,3),mq2(3,3),betamq21(3,3),betamq22(3,3),Dmq2(3,3),ml2(3,3),betaml21(3,3)      & 
& ,betaml22(3,3),Dml2(3,3),md2(3,3),betamd21(3,3),betamd22(3,3),Dmd2(3,3),               & 
& mu2(3,3),betamu21(3,3),betamu22(3,3),Dmu2(3,3),me2(3,3),betame21(3,3),betame22(3,3)    & 
& ,Dme2(3,3),mv2(3,3),betamv21(3,3),betamv22(3,3),Dmv2(3,3),M1,betaM11,betaM12,          & 
& DM1,M2,betaM21,betaM22,DM2,M3,betaM31,betaM32,DM3
Real(dp) :: Tr1(3),Tr2(3),Tr3(3) 
Real(dp) :: Tr2U1(3,3) 
Real(dp) :: AbsM1,AbsM2,AbsM3
Complex(dp) :: md2Yd(3,3),me2Ye(3,3),ml2adjYe(3,3),ml2CYv(3,3),mq2adjYd(3,3),mq2adjYu(3,3),          & 
& mu2Yu(3,3),mv2lam(3),mv2TpYv(3,3),Ydmq2(3,3),YdadjYd(3,3),Yeml2(3,3),YeadjYe(3,3),     & 
& Yumq2(3,3),YuadjYu(3,3),YvadjYv(3,3),YvClam(3),adjYdmd2(3,3),adjYdYd(3,3),             & 
& adjYdTd(3,3),adjYeme2(3,3),adjYeYe(3,3),adjYeTe(3,3),adjYumu2(3,3),adjYuYu(3,3),       & 
& adjYuTu(3,3),adjYvmlHd2(3),adjYvvL(3),adjYvYv(3,3),adjYvCml2(3,3),adjYvTv(3,3),        & 
& adjTdTd(3,3),adjTeTe(3,3),adjTuTu(3,3),adjkap1mv2(3,3),adjkap1lam(3),adjkap2mv2(3,3),  & 
& adjkap2lam(3),adjkap3mv2(3,3),adjkap3lam(3),Cmv2adjYv(3,3),Cmv2adjkap1(3,3),           & 
& Cmv2adjkap2(3,3),Cmv2adjkap3(3,3),Cmv2Clam(3),CYevL(3),CYeYv(3,3),CYeTv(3,3),          & 
& CYvmv2(3,3),CYvvR(3),CYvlam(3),CYvTlam(3),CYvTpYv(3,3),CYvTpTv(3,3),CTdTpTd(3,3),      & 
& CTeTpTe(3,3),CTuTpTu(3,3),CTvTlam(3),CTvTpTv(3,3),Ckap1vR(3),Ckap1Tpkap1(3,3),         & 
& Ckap1TpTk1(3,3),Ckap2vR(3),Ckap2Tpkap2(3,3),Ckap2TpTk2(3,3),Ckap3vR(3),Ckap3Tpkap3(3,3),& 
& Ckap3TpTk3(3,3),TdadjTd(3,3),TeadjTe(3,3),TuadjTu(3,3),TvClam(3),TpYvml2(3,3),         & 
& TpYvCYv(3,3),TpTvCTv(3,3),kap1adjkap1(3,3),kap1Ckap1(3,3),kap1Ckap2(3,3),              & 
& kap1Ckap3(3,3),kap2adjkap2(3,3),kap2Ckap1(3,3),kap2Ckap2(3,3),kap2Ckap3(3,3),          & 
& kap3adjkap3(3,3),kap3Ckap1(3,3),kap3Ckap2(3,3),kap3Ckap3(3,3),Tk1adjTk1(3,3),          & 
& Tk2adjTk2(3,3),Tk3adjTk3(3,3),md2YdadjYd(3,3),me2YeadjYe(3,3),ml2adjYeYe(3,3),         & 
& ml2CYvlam(3),ml2CYvTpYv(3,3),mq2adjYdYd(3,3),mq2adjYuYu(3,3),mu2YuadjYu(3,3),          & 
& mv2TpYvCYv(3,3),mv2kap1adjkap1(3,3),mv2kap2adjkap2(3,3),mv2kap3adjkap3(3,3),           & 
& Ydmq2adjYd(3,3),YdadjYdmd2(3,3),YdadjYdYd(3,3),YdadjYdTd(3,3),YdadjYuYu(3,3),          & 
& YdadjYuTu(3,3),Yeml2adjYe(3,3),YeadjYeme2(3,3),YeadjYeYe(3,3),YeadjYeTe(3,3),          & 
& YeCYvTpYv(3,3),YeCYvTpTv(3,3),Yumq2adjYu(3,3),YuadjYdYd(3,3),YuadjYdTd(3,3),           & 
& YuadjYumu2(3,3),YuadjYuYu(3,3),YuadjYuTu(3,3),YvadjYvvL(3),YvadjYvYv(3,3),             & 
& YvadjYvCml2(3,3),YvadjYvTv(3,3),YvCmv2adjYv(3,3),YvCkap1Tpkap1(3,3),YvCkap1TpTk1(3,3), & 
& YvCkap2Tpkap2(3,3),YvCkap2TpTk2(3,3),YvCkap3Tpkap3(3,3),YvCkap3TpTk3(3,3),             & 
& adjYdmd2Yd(3,3),adjYdYdmq2(3,3),adjYeme2Ye(3,3),adjYeYeml2(3,3),adjYumu2Yu(3,3),       & 
& adjYuYumq2(3,3),CYvmv2lam(3),CYvmv2TpYv(3,3),CYvTpYvml2(3,3),TdadjYdYd(3,3),           & 
& TdadjYuYu(3,3),TeadjYeYe(3,3),TeCYvTpYv(3,3),TuadjYdYd(3,3),TuadjYuYu(3,3),            & 
& TvadjYvYv(3,3),TpYeCYevL(3),TpYeCYeYv(3,3),TpYeCYeTv(3,3),TpYvml2CYv(3,3),             & 
& TpYvCYvmv2(3,3),TpYvCYvvR(3),TpYvCYvlam(3),TpYvCYvTlam(3),TpTeCYeYv(3,3),              & 
& TpTvCYvlam(3),kap1adjkap1mv2(3,3),kap1Cmv2adjkap1(3,3),kap1Ckap1vR(3),kap2adjkap2mv2(3,3),& 
& kap2Cmv2adjkap2(3,3),kap2Ckap2vR(3),kap3adjkap3mv2(3,3),kap3Cmv2adjkap3(3,3),          & 
& kap3Ckap3vR(3),Tk1adjkap1lam(3),Tk2adjkap2lam(3),Tk3adjkap3lam(3)

Complex(dp) :: mv2kap1(3,3),mv2kap2(3,3),mv2kap3(3,3),YdadjYu(3,3),YdadjTd(3,3),YdadjTu(3,3),        & 
& YevL(3),YeadjTe(3,3),YeCYv(3,3),YeCTv(3,3),YuadjYd(3,3),YuadjTd(3,3),YuadjTu(3,3),     & 
& Yvadjkap1(3,3),Yvadjkap2(3,3),Yvadjkap3(3,3),YvadjTk1(3,3),YvadjTk2(3,3),              & 
& YvadjTk3(3,3),YvCTlam(3),YvCkap1(3,3),YvCkap2(3,3),YvCkap3(3,3),adjYdCmd2(3,3),        & 
& adjYeCme2(3,3),adjYuCmu2(3,3),adjTdYd(3,3),adjTeYe(3,3),adjTuYu(3,3),adjkap1vR(3),     & 
& adjkap1Tlam(3),adjkap1TpYv(3,3),adjkap1Tpkap1(3,3),adjkap1Tpkap2(3,3),adjkap1Tpkap3(3,3),& 
& adjkap1kap1(3,3),adjkap1kap2(3,3),adjkap1kap3(3,3),adjkap1Tk1(3,3),adjkap1Tk2(3,3),    & 
& adjkap1Tk3(3,3),adjkap2vR(3),adjkap2Tlam(3),adjkap2TpYv(3,3),adjkap2Tpkap1(3,3),       & 
& adjkap2Tpkap2(3,3),adjkap2Tpkap3(3,3),adjkap2kap1(3,3),adjkap2kap2(3,3),               & 
& adjkap2kap3(3,3),adjkap2Tk1(3,3),adjkap2Tk2(3,3),adjkap2Tk3(3,3),adjkap3vR(3),         & 
& adjkap3Tlam(3),adjkap3TpYv(3,3),adjkap3Tpkap1(3,3),adjkap3Tpkap2(3,3),adjkap3Tpkap3(3,3),& 
& adjkap3kap1(3,3),adjkap3kap2(3,3),adjkap3kap3(3,3),adjkap3Tk1(3,3),adjkap3Tk2(3,3),    & 
& adjkap3Tk3(3,3),adjTk1lam(3),adjTk1TpYv(3,3),adjTk1Tk1(3,3),adjTk2lam(3),              & 
& adjTk2TpYv(3,3),adjTk2Tk2(3,3),adjTk3lam(3),adjTk3TpYv(3,3),adjTk3Tk3(3,3),            & 
& Cme2CYe(3,3),Cml2adjYe(3,3),Cmq2adjYd(3,3),Cmq2adjYu(3,3),Cmv2Ckap1(3,3),              & 
& Cmv2Ckap2(3,3),Cmv2Ckap3(3,3),CYeCml2(3,3),CYvTpkap1(3,3),CYvTpkap2(3,3),              & 
& CYvTpkap3(3,3),CYvTpTk1(3,3),CYvTpTk2(3,3),CYvTpTk3(3,3),CTdTpYd(3,3),CTeTv(3,3),      & 
& CTeTpYe(3,3),CTuTpYu(3,3),CTvlam(3),CTvTpYv(3,3),Ckap1mv2(3,3),Ckap1lam(3),            & 
& Ckap1Tlam(3),Ckap1TpTv(3,3),Ckap1Tpkap2(3,3),Ckap1Tpkap3(3,3),Ckap1TpTk2(3,3),         & 
& Ckap1TpTk3(3,3),Ckap1kap1(3,3),Ckap1kap2(3,3),Ckap1kap3(3,3),Ckap1Tk1(3,3),            & 
& Ckap1Tk2(3,3),Ckap1Tk3(3,3),Ckap2mv2(3,3),Ckap2lam(3),Ckap2Tlam(3),Ckap2TpTv(3,3),     & 
& Ckap2Tpkap1(3,3),Ckap2Tpkap3(3,3),Ckap2TpTk1(3,3),Ckap2TpTk3(3,3),Ckap2kap1(3,3),      & 
& Ckap2kap2(3,3),Ckap2kap3(3,3),Ckap2Tk1(3,3),Ckap2Tk2(3,3),Ckap2Tk3(3,3),               & 
& Ckap3mv2(3,3),Ckap3lam(3),Ckap3Tlam(3),Ckap3TpTv(3,3),Ckap3Tpkap1(3,3),Ckap3Tpkap2(3,3),& 
& Ckap3TpTk1(3,3),Ckap3TpTk2(3,3),Ckap3kap1(3,3),Ckap3kap2(3,3),Ckap3kap3(3,3),          & 
& Ckap3Tk1(3,3),Ckap3Tk2(3,3),Ckap3Tk3(3,3),CTk1lam(3),CTk1Tpkap1(3,3),CTk1Tpkap2(3,3),  & 
& CTk1Tpkap3(3,3),CTk2lam(3),CTk2Tpkap1(3,3),CTk2Tpkap2(3,3),CTk2Tpkap3(3,3),            & 
& CTk3lam(3),CTk3Tpkap1(3,3),CTk3Tpkap2(3,3),CTk3Tpkap3(3,3),TdadjYd(3,3),               & 
& TdadjYu(3,3),TdadjTu(3,3),TeadjYe(3,3),TeCYv(3,3),TeCTv(3,3),TuadjYd(3,3),             & 
& TuadjYu(3,3),TuadjTd(3,3),TvadjYv(3,3),TvadjTv(3,3),TvCTlam(3),TpYeCYe(3,3),           & 
& TpYvmlHd2(3),TpYvvL(3),TpYvadjYe(3,3),TpYvadjTe(3,3),TpYvCTv(3,3),TpTeCTe(3,3),        & 
& TpTvadjYe(3,3),TpTvadjTe(3,3),TpTvCYv(3,3),Tpkap1adjYv(3,3),Tpkap1Clam(3),             & 
& Tpkap1Ckap1(3,3),Tpkap1Ckap2(3,3),Tpkap1Ckap3(3,3),Tpkap2adjYv(3,3),Tpkap2Clam(3),     & 
& Tpkap2Ckap1(3,3),Tpkap2Ckap2(3,3),Tpkap2Ckap3(3,3),Tpkap3adjYv(3,3),Tpkap3Clam(3),     & 
& Tpkap3Ckap1(3,3),Tpkap3Ckap2(3,3),Tpkap3Ckap3(3,3),TpTk1adjYv(3,3),TpTk1Clam(3),       & 
& TpTk2adjYv(3,3),TpTk2Clam(3),TpTk3adjYv(3,3),TpTk3Clam(3),kap1adjYv(3,3),              & 
& kap1adjkap2(3,3),kap1adjkap3(3,3),kap1adjTk1(3,3),kap1adjTk2(3,3),kap1adjTk3(3,3)

Complex(dp) :: kap1Clam(3),kap2adjYv(3,3),kap2adjkap1(3,3),kap2adjkap3(3,3),kap2adjTk1(3,3),          & 
& kap2adjTk2(3,3),kap2adjTk3(3,3),kap2Clam(3),kap3adjYv(3,3),kap3adjkap1(3,3),           & 
& kap3adjkap2(3,3),kap3adjTk1(3,3),kap3adjTk2(3,3),kap3adjTk3(3,3),kap3Clam(3),          & 
& Tk1adjTv(3,3),Tk1Clam(3),Tk2adjTv(3,3),Tk2Clam(3),Tk3adjTv(3,3),Tk3Clam(3),            & 
& md2YdadjYu(3,3),me2YeCYv(3,3),ml2YvadjYv(3,3),mu2YuadjYd(3,3),mv2TpYvadjYe(3,3),       & 
& mv2kap1adjkap2(3,3),mv2kap1adjkap3(3,3),mv2kap1Clam(3),mv2kap1Ckap1(3,3),              & 
& mv2kap1Ckap2(3,3),mv2kap1Ckap3(3,3),mv2kap2adjkap1(3,3),mv2kap2adjkap3(3,3),           & 
& mv2kap2Clam(3),mv2kap2Ckap1(3,3),mv2kap2Ckap2(3,3),mv2kap2Ckap3(3,3),mv2kap3adjkap1(3,3),& 
& mv2kap3adjkap2(3,3),mv2kap3Clam(3),mv2kap3Ckap1(3,3),mv2kap3Ckap2(3,3),mv2kap3Ckap3(3,3),& 
& Ydmq2adjYu(3,3),YdadjYdCmd2(3,3),YdadjYumu2(3,3),YdadjTdTd(3,3),YdCmq2adjYd(3,3),      & 
& Yeml2CYv(3,3),YeadjYeCme2(3,3),YeadjTeTe(3,3),YeCml2adjYe(3,3),YeCYvmv2(3,3),          & 
& YeCYvvR(3),YeCYvlam(3),YeCYvTlam(3),YeCTvTlam(3),Yumq2adjYd(3,3),YuadjYdmd2(3,3),      & 
& YuadjYuCmu2(3,3),YuadjTuTu(3,3),YuCmq2adjYu(3,3),YvadjYvmlHd2(3),Yvadjkap1mv2(3,3),    & 
& Yvadjkap1kap1(3,3),Yvadjkap1kap2(3,3),Yvadjkap1kap3(3,3),Yvadjkap1Tk1(3,3),            & 
& Yvadjkap1Tk2(3,3),Yvadjkap1Tk3(3,3),Yvadjkap2mv2(3,3),Yvadjkap2kap1(3,3),              & 
& Yvadjkap2kap2(3,3),Yvadjkap2kap3(3,3),Yvadjkap2Tk1(3,3),Yvadjkap2Tk2(3,3),             & 
& Yvadjkap2Tk3(3,3),Yvadjkap3mv2(3,3),Yvadjkap3kap1(3,3),Yvadjkap3kap2(3,3),             & 
& Yvadjkap3kap3(3,3),Yvadjkap3Tk1(3,3),Yvadjkap3Tk2(3,3),Yvadjkap3Tk3(3,3),              & 
& YvCmv2adjkap1(3,3),YvCmv2adjkap2(3,3),YvCmv2adjkap3(3,3),YvCmv2Clam(3),YvCkap1lam(3),  & 
& YvCkap1kap2(3,3),YvCkap1kap3(3,3),YvCkap1Tk1(3,3),YvCkap1Tk2(3,3),YvCkap1Tk3(3,3),     & 
& YvCkap2lam(3),YvCkap2kap1(3,3),YvCkap2kap3(3,3),YvCkap2Tk1(3,3),YvCkap2Tk2(3,3),       & 
& YvCkap2Tk3(3,3),YvCkap3lam(3),YvCkap3kap1(3,3),YvCkap3kap2(3,3),YvCkap3Tk1(3,3),       & 
& YvCkap3Tk2(3,3),YvCkap3Tk3(3,3),adjYdYdadjYd(3,3),adjYdYdadjYu(3,3),adjYdYdadjTd(3,3), & 
& adjYdYdadjTu(3,3),adjYdTdadjYd(3,3),adjYdTdadjYu(3,3),adjYdTdadjTd(3,3),               & 
& adjYdTdadjTu(3,3),adjYeYevL(3),adjYeYeadjYe(3,3),adjYeYeadjTe(3,3),adjYeYeCYv(3,3),    & 
& adjYeYeCTv(3,3),adjYeTeadjYe(3,3),adjYeTeadjTe(3,3),adjYeTeCYv(3,3),adjYeTeCTv(3,3),   & 
& adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYuYuadjTd(3,3),adjYuYuadjTu(3,3),               & 
& adjYuTuadjYd(3,3),adjYuTuadjYu(3,3),adjYuTuadjTd(3,3),adjYuTuadjTu(3,3),               & 
& adjYvYvadjYv(3,3),adjYvYvadjkap1(3,3),adjYvYvadjkap2(3,3),adjYvYvadjkap3(3,3),         & 
& adjYvYvadjTk1(3,3),adjYvYvadjTk2(3,3),adjYvYvadjTk3(3,3),adjYvYvClam(3),               & 
& adjYvYvCTlam(3),adjYvYvCkap1(3,3),adjYvYvCkap2(3,3),adjYvYvCkap3(3,3),adjYvTvadjYv(3,3),& 
& adjYvTvadjTv(3,3),adjYvTvClam(3),adjYvTvCTlam(3),adjYvTpYeCYe(3,3),adjYvTpTeCTe(3,3),  & 
& adjTdYdadjYd(3,3),adjTdYdadjYu(3,3),adjTdTdadjYd(3,3),adjTdTdadjYu(3,3),               & 
& adjTeYeadjYe(3,3),adjTeYeCYv(3,3),adjTeTeadjYe(3,3),adjTeTeCYv(3,3),adjTuYuadjYd(3,3), & 
& adjTuYuadjYu(3,3),adjTuTuadjYd(3,3),adjTuTuadjYu(3,3),adjTvYvClam(3),adjTvTvadjYv(3,3),& 
& adjTvTvClam(3),adjkap1mv2lam(3),adjkap1mv2kap1(3,3),adjkap1mv2kap2(3,3),               & 
& adjkap1mv2kap3(3,3),adjkap1TpYvml2(3,3),adjkap1kap1Clam(3),adjkap1kap2Clam(3)

Complex(dp) :: adjkap1kap3Clam(3),adjkap2mv2lam(3),adjkap2mv2kap1(3,3),adjkap2mv2kap2(3,3),           & 
& adjkap2mv2kap3(3,3),adjkap2TpYvml2(3,3),adjkap2kap1Clam(3),adjkap2kap2Clam(3),         & 
& adjkap2kap3Clam(3),adjkap3mv2lam(3),adjkap3mv2kap1(3,3),adjkap3mv2kap2(3,3),           & 
& adjkap3mv2kap3(3,3),adjkap3TpYvml2(3,3),adjkap3kap1Clam(3),adjkap3kap2Clam(3),         & 
& adjkap3kap3Clam(3),Cml2YvadjYv(3,3),Cml2Yvadjkap1(3,3),Cml2Yvadjkap2(3,3),             & 
& Cml2Yvadjkap3(3,3),Cml2YvClam(3),Cml2TpYeCYe(3,3),Cmv2Ckap1kap1(3,3),Cmv2Ckap1kap2(3,3),& 
& Cmv2Ckap1kap3(3,3),Cmv2Ckap2kap1(3,3),Cmv2Ckap2kap2(3,3),Cmv2Ckap2kap3(3,3),           & 
& Cmv2Ckap3kap1(3,3),Cmv2Ckap3kap2(3,3),Cmv2Ckap3kap3(3,3),CYeYvClam(3),CYeTvClam(3),    & 
& CYvTpYvvL(3),CYvTpYvadjYe(3,3),CYvTpYvadjTe(3,3),CYvTpYvCYv(3,3),CYvTpYvCTv(3,3),      & 
& CYvTpTvadjTe(3,3),CYvTpTvCTv(3,3),CYvTpkap1Ckap1(3,3),CYvTpkap1Ckap2(3,3),             & 
& CYvTpkap1Ckap3(3,3),CYvTpkap2Ckap1(3,3),CYvTpkap2Ckap2(3,3),CYvTpkap2Ckap3(3,3),       & 
& CYvTpkap3Ckap1(3,3),CYvTpkap3Ckap2(3,3),CYvTpkap3Ckap3(3,3),CTvTpYvadjYe(3,3),         & 
& CTvTpYvCYv(3,3),CTvTpTvadjYe(3,3),CTvTpTvCYv(3,3),Ckap1mv2lam(3),Ckap1TpYvCYv(3,3),    & 
& Ckap1TpYvCTv(3,3),Ckap1Tpkap1adjYv(3,3),Ckap1Tpkap1Clam(3),Ckap1TpTk1adjYv(3,3),       & 
& Ckap1kap1adjYv(3,3),Ckap1kap1Ckap1(3,3),Ckap1kap1Ckap2(3,3),Ckap1kap1Ckap3(3,3),       & 
& Ckap1kap2Ckap1(3,3),Ckap1kap2Ckap2(3,3),Ckap1kap2Ckap3(3,3),Ckap1kap3Ckap1(3,3),       & 
& Ckap1kap3Ckap2(3,3),Ckap1kap3Ckap3(3,3),Ckap1Tk1adjTv(3,3),Ckap2mv2lam(3),             & 
& Ckap2TpYvCYv(3,3),Ckap2TpYvCTv(3,3),Ckap2Tpkap2adjYv(3,3),Ckap2Tpkap2Clam(3),          & 
& Ckap2TpTk2adjYv(3,3),Ckap2kap1Ckap1(3,3),Ckap2kap1Ckap2(3,3),Ckap2kap1Ckap3(3,3),      & 
& Ckap2kap2adjYv(3,3),Ckap2kap2Ckap1(3,3),Ckap2kap2Ckap2(3,3),Ckap2kap2Ckap3(3,3),       & 
& Ckap2kap3Ckap1(3,3),Ckap2kap3Ckap2(3,3),Ckap2kap3Ckap3(3,3),Ckap2Tk2adjTv(3,3),        & 
& Ckap3mv2lam(3),Ckap3TpYvCYv(3,3),Ckap3TpYvCTv(3,3),Ckap3Tpkap3adjYv(3,3),              & 
& Ckap3Tpkap3Clam(3),Ckap3TpTk3adjYv(3,3),Ckap3kap1Ckap1(3,3),Ckap3kap1Ckap2(3,3),       & 
& Ckap3kap1Ckap3(3,3),Ckap3kap2Ckap1(3,3),Ckap3kap2Ckap2(3,3),Ckap3kap2Ckap3(3,3),       & 
& Ckap3kap3adjYv(3,3),Ckap3kap3Ckap1(3,3),Ckap3kap3Ckap2(3,3),Ckap3kap3Ckap3(3,3),       & 
& Ckap3Tk3adjTv(3,3),CTk1TpTk1adjYv(3,3),CTk2TpTk2adjYv(3,3),CTk3TpTk3adjYv(3,3),        & 
& TdadjTdYd(3,3),TeadjTeYe(3,3),TeCYvlam(3),TeCTvlam(3),TeCTvTpYv(3,3),TuadjTuYu(3,3),   & 
& Tvadjkap1kap1(3,3),Tvadjkap1kap2(3,3),Tvadjkap1kap3(3,3),Tvadjkap2kap1(3,3),           & 
& Tvadjkap2kap2(3,3),Tvadjkap2kap3(3,3),Tvadjkap3kap1(3,3),Tvadjkap3kap2(3,3),           & 
& Tvadjkap3kap3(3,3),TvCkap1kap1(3,3),TvCkap1kap2(3,3),TvCkap1kap3(3,3),TvCkap2kap1(3,3),& 
& TvCkap2kap2(3,3),TvCkap2kap3(3,3),TvCkap3kap1(3,3),TvCkap3kap2(3,3),TvCkap3kap3(3,3),  & 
& TpYeCme2CYe(3,3),TpYeCYeCml2(3,3),TpYeCTeTv(3,3),TpYvml2adjYe(3,3),TpYvadjYeme2(3,3),  & 
& TpYvadjYeYe(3,3),TpYvadjYeTe(3,3),TpYvCYvTpYv(3,3),TpYvCYvTpTv(3,3),TpYvCYvTpkap1(3,3),& 
& TpYvCYvTpkap2(3,3),TpYvCYvTpkap3(3,3),TpYvCYvTpTk1(3,3),TpYvCYvTpTk2(3,3),             & 
& TpYvCYvTpTk3(3,3),TpYvCTvTlam(3),TpYvCTvTpTv(3,3),TpTvadjYeYe(3,3),TpTvCYvTpYv(3,3),   & 
& TpTvCTvlam(3),TpTvCTvTpYv(3,3),Tpkap1adjYvvL(3),Tpkap1adjYvYv(3,3),Tpkap1adjkap1lam(3),& 
& Tpkap1adjkap1Tlam(3),Tpkap1adjkap2lam(3),Tpkap1adjkap2Tlam(3),Tpkap1adjkap3lam(3)

Complex(dp) :: Tpkap1adjkap3Tlam(3),Tpkap1Cmv2adjYv(3,3),Tpkap1Cmv2adjkap1(3,3),Tpkap1Cmv2adjkap2(3,3),& 
& Tpkap1Cmv2adjkap3(3,3),Tpkap1Cmv2Clam(3),Tpkap1Ckap1lam(3),Tpkap1Ckap1Tlam(3),         & 
& Tpkap1Ckap1TpTv(3,3),Tpkap1Ckap1TpTk1(3,3),Tpkap1Ckap1TpTk2(3,3),Tpkap1Ckap1TpTk3(3,3),& 
& Tpkap1Ckap2lam(3),Tpkap1Ckap2Tlam(3),Tpkap1Ckap2TpTv(3,3),Tpkap1Ckap2TpTk1(3,3),       & 
& Tpkap1Ckap2TpTk2(3,3),Tpkap1Ckap2TpTk3(3,3),Tpkap1Ckap3lam(3),Tpkap1Ckap3Tlam(3),      & 
& Tpkap1Ckap3TpTv(3,3),Tpkap1Ckap3TpTk1(3,3),Tpkap1Ckap3TpTk2(3,3),Tpkap1Ckap3TpTk3(3,3),& 
& Tpkap2adjYvvL(3),Tpkap2adjYvYv(3,3),Tpkap2adjkap1lam(3),Tpkap2adjkap1Tlam(3),          & 
& Tpkap2adjkap2lam(3),Tpkap2adjkap2Tlam(3),Tpkap2adjkap3lam(3),Tpkap2adjkap3Tlam(3),     & 
& Tpkap2Cmv2adjYv(3,3),Tpkap2Cmv2adjkap1(3,3),Tpkap2Cmv2adjkap2(3,3),Tpkap2Cmv2adjkap3(3,3),& 
& Tpkap2Cmv2Clam(3),Tpkap2Ckap1lam(3),Tpkap2Ckap1Tlam(3),Tpkap2Ckap1TpTv(3,3),           & 
& Tpkap2Ckap1TpTk1(3,3),Tpkap2Ckap1TpTk2(3,3),Tpkap2Ckap1TpTk3(3,3),Tpkap2Ckap2lam(3),   & 
& Tpkap2Ckap2Tlam(3),Tpkap2Ckap2TpTv(3,3),Tpkap2Ckap2TpTk1(3,3),Tpkap2Ckap2TpTk2(3,3),   & 
& Tpkap2Ckap2TpTk3(3,3),Tpkap2Ckap3lam(3),Tpkap2Ckap3Tlam(3),Tpkap2Ckap3TpTv(3,3),       & 
& Tpkap2Ckap3TpTk1(3,3),Tpkap2Ckap3TpTk2(3,3),Tpkap2Ckap3TpTk3(3,3),Tpkap3adjYvvL(3),    & 
& Tpkap3adjYvYv(3,3),Tpkap3adjkap1lam(3),Tpkap3adjkap1Tlam(3),Tpkap3adjkap2lam(3),       & 
& Tpkap3adjkap2Tlam(3),Tpkap3adjkap3lam(3),Tpkap3adjkap3Tlam(3),Tpkap3Cmv2adjYv(3,3),    & 
& Tpkap3Cmv2adjkap1(3,3),Tpkap3Cmv2adjkap2(3,3),Tpkap3Cmv2adjkap3(3,3),Tpkap3Cmv2Clam(3),& 
& Tpkap3Ckap1lam(3),Tpkap3Ckap1Tlam(3),Tpkap3Ckap1TpTv(3,3),Tpkap3Ckap1TpTk1(3,3),       & 
& Tpkap3Ckap1TpTk2(3,3),Tpkap3Ckap1TpTk3(3,3),Tpkap3Ckap2lam(3),Tpkap3Ckap2Tlam(3),      & 
& Tpkap3Ckap2TpTv(3,3),Tpkap3Ckap2TpTk1(3,3),Tpkap3Ckap2TpTk2(3,3),Tpkap3Ckap2TpTk3(3,3),& 
& Tpkap3Ckap3lam(3),Tpkap3Ckap3Tlam(3),Tpkap3Ckap3TpTv(3,3),Tpkap3Ckap3TpTk1(3,3),       & 
& Tpkap3Ckap3TpTk2(3,3),Tpkap3Ckap3TpTk3(3,3),TpTk1adjYvYv(3,3),TpTk1adjkap1lam(3),      & 
& TpTk1adjkap2lam(3),TpTk1adjkap3lam(3),TpTk1Ckap1lam(3),TpTk1Ckap2lam(3),               & 
& TpTk1Ckap3lam(3),TpTk2adjYvYv(3,3),TpTk2adjkap1lam(3),TpTk2adjkap2lam(3),              & 
& TpTk2adjkap3lam(3),TpTk2Ckap1lam(3),TpTk2Ckap2lam(3),TpTk2Ckap3lam(3),TpTk3adjYvYv(3,3),& 
& TpTk3adjkap1lam(3),TpTk3adjkap2lam(3),TpTk3adjkap3lam(3),TpTk3Ckap1lam(3),             & 
& TpTk3Ckap2lam(3),TpTk3Ckap3lam(3),kap1adjYvCml2(3,3),kap1adjYvTv(3,3),kap1adjkap1lam(3),& 
& kap1adjkap1TpYv(3,3),kap1adjkap1Tpkap1(3,3),kap1adjkap1Tpkap2(3,3),kap1adjkap1Tpkap3(3,3),& 
& kap1Cmv2Clam(3),kap1Ckap1lam(3),kap1Ckap1kap1(3,3),kap1Ckap1kap2(3,3),kap1Ckap1kap3(3,3),& 
& kap1Ckap1Tk1(3,3),kap1Ckap1Tk2(3,3),kap1Ckap1Tk3(3,3),kap1Ckap2lam(3),kap1Ckap2kap1(3,3),& 
& kap1Ckap2kap2(3,3),kap1Ckap2kap3(3,3),kap1Ckap2Tk1(3,3),kap1Ckap2Tk2(3,3),             & 
& kap1Ckap2Tk3(3,3),kap1Ckap3lam(3),kap1Ckap3kap1(3,3),kap1Ckap3kap2(3,3),               & 
& kap1Ckap3kap3(3,3),kap1Ckap3Tk1(3,3),kap1Ckap3Tk2(3,3),kap1Ckap3Tk3(3,3),              & 
& kap2adjYvCml2(3,3),kap2adjYvTv(3,3),kap2adjkap2lam(3),kap2adjkap2TpYv(3,3),            & 
& kap2adjkap2Tpkap1(3,3),kap2adjkap2Tpkap2(3,3),kap2adjkap2Tpkap3(3,3),kap2Cmv2Clam(3),  & 
& kap2Ckap1lam(3),kap2Ckap1kap1(3,3),kap2Ckap1kap2(3,3),kap2Ckap1kap3(3,3),              & 
& kap2Ckap1Tk1(3,3),kap2Ckap1Tk2(3,3),kap2Ckap1Tk3(3,3),kap2Ckap2lam(3),kap2Ckap2kap1(3,3)

Complex(dp) :: kap2Ckap2kap2(3,3),kap2Ckap2kap3(3,3),kap2Ckap2Tk1(3,3),kap2Ckap2Tk2(3,3),             & 
& kap2Ckap2Tk3(3,3),kap2Ckap3lam(3),kap2Ckap3kap1(3,3),kap2Ckap3kap2(3,3),               & 
& kap2Ckap3kap3(3,3),kap2Ckap3Tk1(3,3),kap2Ckap3Tk2(3,3),kap2Ckap3Tk3(3,3),              & 
& kap3adjYvCml2(3,3),kap3adjYvTv(3,3),kap3adjkap3lam(3),kap3adjkap3TpYv(3,3),            & 
& kap3adjkap3Tpkap1(3,3),kap3adjkap3Tpkap2(3,3),kap3adjkap3Tpkap3(3,3),kap3Cmv2Clam(3),  & 
& kap3Ckap1lam(3),kap3Ckap1kap1(3,3),kap3Ckap1kap2(3,3),kap3Ckap1kap3(3,3),              & 
& kap3Ckap1Tk1(3,3),kap3Ckap1Tk2(3,3),kap3Ckap1Tk3(3,3),kap3Ckap2lam(3),kap3Ckap2kap1(3,3),& 
& kap3Ckap2kap2(3,3),kap3Ckap2kap3(3,3),kap3Ckap2Tk1(3,3),kap3Ckap2Tk2(3,3),             & 
& kap3Ckap2Tk3(3,3),kap3Ckap3lam(3),kap3Ckap3kap1(3,3),kap3Ckap3kap2(3,3),               & 
& kap3Ckap3kap3(3,3),kap3Ckap3Tk1(3,3),kap3Ckap3Tk2(3,3),kap3Ckap3Tk3(3,3),              & 
& Tk1adjkap1TpYv(3,3),Tk1adjTk1lam(3),Tk1adjTk1TpYv(3,3),Tk1Ckap1kap1(3,3),              & 
& Tk1Ckap1kap2(3,3),Tk1Ckap1kap3(3,3),Tk1Ckap2kap1(3,3),Tk1Ckap2kap2(3,3),               & 
& Tk1Ckap2kap3(3,3),Tk1Ckap3kap1(3,3),Tk1Ckap3kap2(3,3),Tk1Ckap3kap3(3,3),               & 
& Tk2adjkap2TpYv(3,3),Tk2adjTk2lam(3),Tk2adjTk2TpYv(3,3),Tk2Ckap1kap1(3,3),              & 
& Tk2Ckap1kap2(3,3),Tk2Ckap1kap3(3,3),Tk2Ckap2kap1(3,3),Tk2Ckap2kap2(3,3),               & 
& Tk2Ckap2kap3(3,3),Tk2Ckap3kap1(3,3),Tk2Ckap3kap2(3,3),Tk2Ckap3kap3(3,3),               & 
& Tk3adjkap3TpYv(3,3),Tk3adjTk3lam(3),Tk3adjTk3TpYv(3,3),Tk3Ckap1kap1(3,3),              & 
& Tk3Ckap1kap2(3,3),Tk3Ckap1kap3(3,3),Tk3Ckap2kap1(3,3),Tk3Ckap2kap2(3,3),               & 
& Tk3Ckap2kap3(3,3),Tk3Ckap3kap1(3,3),Tk3Ckap3kap2(3,3),Tk3Ckap3kap3(3,3),               & 
& md2YdadjYdYd(3,3),me2YeadjYeYe(3,3),me2YeCYvlam(3),ml2adjYeYeadjYe(3,3),               & 
& ml2adjYeYeCYv(3,3),ml2CYvTpYvadjYe(3,3),ml2CYvTpYvCYv(3,3),mq2adjYdYdadjYd(3,3),       & 
& mq2adjYdYdadjYu(3,3),mq2adjYuYuadjYd(3,3),mq2adjYuYuadjYu(3,3),mu2YuadjYuYu(3,3),      & 
& mv2TpYvCYvlam(3),mv2TpYvCYvTpYv(3,3),mv2kap1adjkap1TpYv(3,3),mv2kap2adjkap2TpYv(3,3),  & 
& mv2kap3adjkap3TpYv(3,3),Ydmq2adjYdYd(3,3),YdadjYdmd2Yd(3,3),YdadjYdYdmq2(3,3),         & 
& YdadjYdYdadjYd(3,3),YdadjYdTdadjYd(3,3),YdadjYdTdadjTd(3,3),YdadjYuYuadjYd(3,3),       & 
& YdadjYuTuadjYd(3,3),YdadjYuTuadjTd(3,3),YdadjTdTdadjYd(3,3),YdadjTuTuadjYd(3,3),       & 
& Yeml2adjYeYe(3,3),Yeml2CYvlam(3),YeadjYeme2Ye(3,3),YeadjYeYeml2(3,3),YeadjYeYeadjYe(3,3),& 
& YeadjYeTeadjYe(3,3),YeadjYeTeadjTe(3,3),YeadjTeTeadjYe(3,3),YeCYvmv2lam(3),            & 
& YeCYvTpYvadjYe(3,3),YeCYvTpTvadjTe(3,3),YeCTvTpTvadjYe(3,3),Yumq2adjYuYu(3,3),         & 
& YuadjYdYdadjYu(3,3),YuadjYdTdadjYu(3,3),YuadjYdTdadjTu(3,3),YuadjYumu2Yu(3,3),         & 
& YuadjYuYumq2(3,3),YuadjYuYuadjYu(3,3),YuadjYuTuadjYu(3,3),YuadjYuTuadjTu(3,3),         & 
& YuadjTdTdadjYu(3,3),YuadjTuTuadjYu(3,3),YvadjYvYvadjYv(3,3),YvadjYvYvClam(3),          & 
& YvadjYvTvadjYv(3,3),YvadjYvTvadjTv(3,3),YvadjYvTvClam(3),YvadjYvTpYeCYe(3,3),          & 
& YvadjYvTpTeCTe(3,3),YvadjTvTvadjYv(3,3),Yvadjkap1mv2kap1(3,3),Yvadjkap1mv2kap2(3,3),   & 
& Yvadjkap1mv2kap3(3,3),Yvadjkap2mv2kap1(3,3),Yvadjkap2mv2kap2(3,3),Yvadjkap2mv2kap3(3,3),& 
& Yvadjkap3mv2kap1(3,3),Yvadjkap3mv2kap2(3,3),Yvadjkap3mv2kap3(3,3),YvCmv2Ckap1kap1(3,3),& 
& YvCmv2Ckap1kap2(3,3),YvCmv2Ckap1kap3(3,3),YvCmv2Ckap2kap1(3,3),YvCmv2Ckap2kap2(3,3)

Complex(dp) :: YvCmv2Ckap2kap3(3,3),YvCmv2Ckap3kap1(3,3),YvCmv2Ckap3kap2(3,3),YvCmv2Ckap3kap3(3,3),   & 
& YvCkap1Tpkap1adjYv(3,3),YvCkap1Tpkap1Clam(3),YvCkap1TpTk1adjYv(3,3),YvCkap1kap1adjYv(3,3),& 
& YvCkap1Tk1adjTv(3,3),YvCkap2Tpkap2adjYv(3,3),YvCkap2Tpkap2Clam(3),YvCkap2TpTk2adjYv(3,3),& 
& YvCkap2kap2adjYv(3,3),YvCkap2Tk2adjTv(3,3),YvCkap3Tpkap3adjYv(3,3),YvCkap3Tpkap3Clam(3),& 
& YvCkap3TpTk3adjYv(3,3),YvCkap3kap3adjYv(3,3),YvCkap3Tk3adjTv(3,3),YvCTk1TpTk1adjYv(3,3),& 
& YvCTk2TpTk2adjYv(3,3),YvCTk3TpTk3adjYv(3,3),adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3), & 
& adjYdYdmq2adjYd(3,3),adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3),adjYdYdadjYdYd(3,3),    & 
& adjYdYdadjYdTd(3,3),adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3),adjYdYdadjYuTu(3,3),      & 
& adjYdYdadjTdTd(3,3),adjYdTdadjYdYd(3,3),adjYdTdadjYuYu(3,3),adjYdTdadjTdYd(3,3),       & 
& adjYeme2YeadjYe(3,3),adjYeme2YeCYv(3,3),adjYeYeml2adjYe(3,3),adjYeYeml2CYv(3,3),       & 
& adjYeYeadjYeme2(3,3),adjYeYeadjYeYe(3,3),adjYeYeadjYeTe(3,3),adjYeYeadjTeTe(3,3),      & 
& adjYeYeCYvmv2(3,3),adjYeYeCYvvR(3),adjYeYeCYvlam(3),adjYeYeCTvTlam(3),adjYeTeadjYeYe(3,3),& 
& adjYeTeadjTeYe(3,3),adjYeTeCYvlam(3),adjYeTeCYvTpYv(3,3),adjYeTeCTvlam(3),             & 
& adjYeTeCTvTpYv(3,3),adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),    & 
& adjYuYumq2adjYu(3,3),adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdTd(3,3),     & 
& adjYuYuadjYumu2(3,3),adjYuYuadjYuYu(3,3),adjYuYuadjYuTu(3,3),adjYuYuadjTuTu(3,3),      & 
& adjYuTuadjYdYd(3,3),adjYuTuadjYuYu(3,3),adjYuTuadjTuYu(3,3),adjYvYvadjYvmlHd2(3),      & 
& adjYvYvadjYvvL(3),adjYvYvadjYvYv(3,3),adjYvYvadjYvTv(3,3),adjYvYvadjkap1mv2(3,3),      & 
& adjYvYvadjkap2mv2(3,3),adjYvYvadjkap3mv2(3,3),adjYvYvCmv2adjkap1(3,3),adjYvYvCmv2adjkap2(3,3),& 
& adjYvYvCmv2adjkap3(3,3),adjYvYvCmv2Clam(3),adjYvYvCkap1lam(3),adjYvYvCkap1kap2(3,3),   & 
& adjYvYvCkap1kap3(3,3),adjYvYvCkap2lam(3),adjYvYvCkap2kap1(3,3),adjYvYvCkap2kap3(3,3),  & 
& adjYvYvCkap3lam(3),adjYvYvCkap3kap1(3,3),adjYvYvCkap3kap2(3,3),adjYvCml2YvadjYv(3,3),  & 
& adjYvCml2Yvadjkap1(3,3),adjYvCml2Yvadjkap2(3,3),adjYvCml2Yvadjkap3(3,3),               & 
& adjYvCml2YvClam(3),adjYvCml2TpYeCYe(3,3),adjYvTvadjYvYv(3,3),adjYvTvadjkap1kap1(3,3),  & 
& adjYvTvadjkap1kap2(3,3),adjYvTvadjkap1kap3(3,3),adjYvTvadjkap2kap1(3,3),               & 
& adjYvTvadjkap2kap2(3,3),adjYvTvadjkap2kap3(3,3),adjYvTvadjkap3kap1(3,3),               & 
& adjYvTvadjkap3kap2(3,3),adjYvTvadjkap3kap3(3,3),adjYvTpYeCme2CYe(3,3),adjYvTpYeCYeYv(3,3),& 
& adjYvTpYeCYeCml2(3,3),adjYvTpYeCYeTv(3,3),adjYvTpYeCTeTv(3,3),adjYvTpTeCYeYv(3,3),     & 
& adjTdYdadjYdTd(3,3),adjTdTdadjYdYd(3,3),adjTeYeadjYeTe(3,3),adjTeYeCYvTlam(3),         & 
& adjTeTeadjYeYe(3,3),adjTeTeCYvlam(3),adjTuYuadjYuTu(3,3),adjTuTuadjYuYu(3,3),          & 
& adjTvYvCkap1Tk2(3,3),adjTvYvCkap1Tk3(3,3),adjTvYvCkap2Tk1(3,3),adjTvYvCkap2Tk3(3,3),   & 
& adjTvYvCkap3Tk1(3,3),adjTvYvCkap3Tk2(3,3),adjkap1Tpkap1Cmv2Clam(3),adjkap1Tpkap2Cmv2Clam(3),& 
& adjkap1Tpkap3Cmv2Clam(3),adjkap2Tpkap1Cmv2Clam(3),adjkap2Tpkap2Cmv2Clam(3),            & 
& adjkap2Tpkap3Cmv2Clam(3),adjkap3Tpkap1Cmv2Clam(3),adjkap3Tpkap2Cmv2Clam(3),            & 
& adjkap3Tpkap3Cmv2Clam(3),Cml2YvCkap1kap2(3,3),Cml2YvCkap1kap3(3,3),Cml2YvCkap2kap1(3,3),& 
& Cml2YvCkap2kap3(3,3),Cml2YvCkap3kap1(3,3),Cml2YvCkap3kap2(3,3),Cmv2adjYvYvadjYv(3,3),  & 
& Cmv2adjYvYvadjkap1(3,3),Cmv2adjYvYvadjkap2(3,3),Cmv2adjYvYvadjkap3(3,3)

Complex(dp) :: Cmv2adjYvYvClam(3),Cmv2adjYvTpYeCYe(3,3),Cmv2Ckap1kap1adjYv(3,3),Cmv2Ckap2kap2adjYv(3,3),& 
& Cmv2Ckap3kap3adjYv(3,3),CYeTpYeCYevL(3),CYeTpYeCYeYv(3,3),CYeTpYeCYeTv(3,3),           & 
& CYeTpTeCYeYv(3,3),CYvmv2TpYvadjYe(3,3),CYvmv2TpYvCYv(3,3),CYvmv2kap1Ckap1(3,3),        & 
& CYvmv2kap1Ckap2(3,3),CYvmv2kap1Ckap3(3,3),CYvmv2kap2Ckap1(3,3),CYvmv2kap2Ckap2(3,3),   & 
& CYvmv2kap2Ckap3(3,3),CYvmv2kap3Ckap1(3,3),CYvmv2kap3Ckap2(3,3),CYvmv2kap3Ckap3(3,3),   & 
& CYvTpYvml2adjYe(3,3),CYvTpYvml2CYv(3,3),CYvTpYvadjYeme2(3,3),CYvTpYvadjYeYe(3,3),      & 
& CYvTpYvadjYeTe(3,3),CYvTpYvCYvmv2(3,3),CYvTpYvCYvvR(3),CYvTpYvCYvlam(3),               & 
& CYvTpYvCYvTlam(3),CYvTpYvCYvTpYv(3,3),CYvTpYvCYvTpTv(3,3),CYvTpYvCTvTlam(3),           & 
& CYvTpYvCTvTpTv(3,3),CYvTpTvadjYeYe(3,3),CYvTpTvCYvlam(3),CYvTpTvCYvTpYv(3,3),          & 
& CYvTpTvCTvlam(3),CYvTpTvCTvTpYv(3,3),CYvkap1adjkap1lam(3),CYvkap1adjkap1TpYv(3,3),     & 
& CYvkap2adjkap2lam(3),CYvkap2adjkap2TpYv(3,3),CYvkap3adjkap3lam(3),CYvkap3adjkap3TpYv(3,3),& 
& CYvTk1adjkap1lam(3),CYvTk1adjkap1TpYv(3,3),CYvTk1adjTk1lam(3),CYvTk1adjTk1TpYv(3,3),   & 
& CYvTk2adjkap2lam(3),CYvTk2adjkap2TpYv(3,3),CYvTk2adjTk2lam(3),CYvTk2adjTk2TpYv(3,3),   & 
& CYvTk3adjkap3lam(3),CYvTk3adjkap3TpYv(3,3),CYvTk3adjTk3lam(3),CYvTk3adjTk3TpYv(3,3),   & 
& CTvTpYvCYvTlam(3),CTvTpYvCYvTpTv(3,3),CTvTpTvCYvlam(3),CTvTpTvCYvTpYv(3,3),            & 
& Ckap1TpYvCYvTpTk1(3,3),Ckap1Tpkap1adjYvvL(3),Ckap1Tpkap1adjYvYv(3,3),Ckap1Tpkap1Cmv2adjYv(3,3),& 
& Ckap1Tpkap1Cmv2adjkap1(3,3),Ckap1Tpkap1Cmv2adjkap2(3,3),Ckap1Tpkap1Cmv2adjkap3(3,3),   & 
& Ckap1TpTk1adjYvYv(3,3),Ckap1kap1adjYvCml2(3,3),Ckap2TpYvCYvTpTk2(3,3),Ckap2Tpkap2adjYvvL(3),& 
& Ckap2Tpkap2adjYvYv(3,3),Ckap2Tpkap2Cmv2adjYv(3,3),Ckap2Tpkap2Cmv2adjkap1(3,3),         & 
& Ckap2Tpkap2Cmv2adjkap2(3,3),Ckap2Tpkap2Cmv2adjkap3(3,3),Ckap2TpTk2adjYvYv(3,3),        & 
& Ckap2kap2adjYvCml2(3,3),Ckap3TpYvCYvTpTk3(3,3),Ckap3Tpkap3adjYvvL(3),Ckap3Tpkap3adjYvYv(3,3),& 
& Ckap3Tpkap3Cmv2adjYv(3,3),Ckap3Tpkap3Cmv2adjkap1(3,3),Ckap3Tpkap3Cmv2adjkap2(3,3),     & 
& Ckap3Tpkap3Cmv2adjkap3(3,3),Ckap3TpTk3adjYvYv(3,3),Ckap3kap3adjYvCml2(3,3),            & 
& TdadjYdYdadjTd(3,3),TdadjYuYuadjTd(3,3),TdadjTdYdadjYd(3,3),TdadjTuYuadjYd(3,3),       & 
& TeadjYeYeadjTe(3,3),TeadjTeYeadjYe(3,3),TeCYvTpYvadjTe(3,3),TeCTvTpYvadjYe(3,3),       & 
& TuadjYdYdadjTu(3,3),TuadjYuYuadjTu(3,3),TuadjTdYdadjYu(3,3),TuadjTuYuadjYu(3,3),       & 
& TvadjYvYvClam(3),TpYeCYeYvClam(3),TpYeCYeTvClam(3),TpYvml2CYvlam(3),TpYvml2CYvTpYv(3,3),& 
& TpYvadjYeYeCYv(3,3),TpYvadjYeTeCYv(3,3),TpYvadjYeTeCTv(3,3),TpYvadjTeTeCYv(3,3),       & 
& TpYvCYvmv2lam(3),TpYvCYvmv2TpYv(3,3),TpYvCYvTpYvml2(3,3),TpYvCYvTpYvCYv(3,3),          & 
& TpYvCYvTpTvCTv(3,3),TpYvCTvTpTvCYv(3,3),TpTeCYeYvClam(3),TpTvadjYeYeCTv(3,3),          & 
& TpTvadjTeYeCYv(3,3),TpTvCYvTpYvCYv(3,3),TpTvCYvTpYvCTv(3,3),TpTvCTvTpYvCYv(3,3),       & 
& Tpkap1adjkap1mv2lam(3),Tpkap1adjkap2mv2lam(3),Tpkap1adjkap3mv2lam(3),Tpkap2adjkap1mv2lam(3),& 
& Tpkap2adjkap2mv2lam(3),Tpkap2adjkap3mv2lam(3),Tpkap3adjkap1mv2lam(3),Tpkap3adjkap2mv2lam(3),& 
& Tpkap3adjkap3mv2lam(3),kap1adjYvYvadjkap1(3,3),kap1adjYvYvCkap1(3,3),kap1adjYvYvCkap2(3,3),& 
& kap1adjYvYvCkap3(3,3),kap1adjkap1TpYvml2(3,3),kap1Cmv2Ckap1kap1(3,3),kap1Cmv2Ckap1kap2(3,3),& 
& kap1Cmv2Ckap1kap3(3,3),kap1Cmv2Ckap2kap1(3,3),kap1Cmv2Ckap2kap2(3,3),kap1Cmv2Ckap2kap3(3,3),& 
& kap1Cmv2Ckap3kap1(3,3),kap1Cmv2Ckap3kap2(3,3),kap1Cmv2Ckap3kap3(3,3),kap2adjYvYvadjkap2(3,3)

Complex(dp) :: kap2adjYvYvCkap1(3,3),kap2adjYvYvCkap2(3,3),kap2adjYvYvCkap3(3,3),kap2adjkap2TpYvml2(3,3),& 
& kap2Cmv2Ckap1kap1(3,3),kap2Cmv2Ckap1kap2(3,3),kap2Cmv2Ckap1kap3(3,3),kap2Cmv2Ckap2kap1(3,3),& 
& kap2Cmv2Ckap2kap2(3,3),kap2Cmv2Ckap2kap3(3,3),kap2Cmv2Ckap3kap1(3,3),kap2Cmv2Ckap3kap2(3,3),& 
& kap2Cmv2Ckap3kap3(3,3),kap3adjYvYvadjkap3(3,3),kap3adjYvYvCkap1(3,3),kap3adjYvYvCkap2(3,3),& 
& kap3adjYvYvCkap3(3,3),kap3adjkap3TpYvml2(3,3),kap3Cmv2Ckap1kap1(3,3),kap3Cmv2Ckap1kap2(3,3),& 
& kap3Cmv2Ckap1kap3(3,3),kap3Cmv2Ckap2kap1(3,3),kap3Cmv2Ckap2kap2(3,3),kap3Cmv2Ckap2kap3(3,3),& 
& kap3Cmv2Ckap3kap1(3,3),kap3Cmv2Ckap3kap2(3,3),kap3Cmv2Ckap3kap3(3,3),Tk1adjYvYvadjTk1(3,3),& 
& Tk1adjYvYvCkap1(3,3),Tk2adjYvYvadjTk2(3,3),Tk2adjYvYvCkap2(3,3),Tk3adjYvYvadjTk3(3,3), & 
& Tk3adjYvYvCkap3(3,3),md2YdadjYdYdadjYd(3,3),md2YdadjYuYuadjYd(3,3),me2YeadjYeYeadjYe(3,3),& 
& me2YeCYvTpYvadjYe(3,3),ml2adjYeYeadjYeYe(3,3),ml2adjYeYeCYvlam(3),ml2CYvTpYvCYvlam(3), & 
& ml2CYvTpYvCYvTpYv(3,3),mq2adjYdYdadjYdYd(3,3),mq2adjYdYdadjYuYu(3,3),mq2adjYuYuadjYdYd(3,3),& 
& mq2adjYuYuadjYuYu(3,3),mu2YuadjYdYdadjYu(3,3),mu2YuadjYuYuadjYu(3,3),mv2TpYvadjYeYeCYv(3,3),& 
& mv2TpYvCYvTpYvCYv(3,3),mv2kap1adjYvYvadjkap1(3,3),mv2kap2adjYvYvadjkap2(3,3),          & 
& mv2kap3adjYvYvadjkap3(3,3),Ydmq2adjYdYdadjYd(3,3),Ydmq2adjYuYuadjYd(3,3),              & 
& YdadjYdmd2YdadjYd(3,3),YdadjYdYdmq2adjYd(3,3),YdadjYdYdadjYdmd2(3,3),YdadjYdYdadjYdYd(3,3),& 
& YdadjYdYdadjYdTd(3,3),YdadjYdTdadjYdYd(3,3),YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),& 
& YdadjYuYuadjYdmd2(3,3),YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYdTd(3,3),YdadjYuYuadjYuYu(3,3),& 
& YdadjYuYuadjYuTu(3,3),YdadjYuTuadjYdYd(3,3),YdadjYuTuadjYuYu(3,3),Yeml2adjYeYeadjYe(3,3),& 
& Yeml2CYvTpYvadjYe(3,3),YeadjYeme2YeadjYe(3,3),YeadjYeYeml2adjYe(3,3),YeadjYeYeadjYeme2(3,3),& 
& YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeTe(3,3),YeadjYeTeadjYeYe(3,3),YeCYvmv2TpYvadjYe(3,3),& 
& YeCYvTpYvml2adjYe(3,3),YeCYvTpYvadjYeme2(3,3),YeCYvTpYvadjYeYe(3,3),YeCYvTpYvadjYeTe(3,3),& 
& YeCYvTpYvCYvTpYv(3,3),YeCYvTpYvCYvTpTv(3,3),YeCYvTpTvadjYeYe(3,3),YeCYvTpTvCYvTpYv(3,3),& 
& YeCYvkap1adjkap1TpYv(3,3),YeCYvkap2adjkap2TpYv(3,3),YeCYvkap3adjkap3TpYv(3,3),         & 
& YeCYvTk1adjkap1TpYv(3,3),YeCYvTk2adjkap2TpYv(3,3),YeCYvTk3adjkap3TpYv(3,3),            & 
& Yumq2adjYdYdadjYu(3,3),Yumq2adjYuYuadjYu(3,3),YuadjYdmd2YdadjYu(3,3),YuadjYdYdmq2adjYu(3,3),& 
& YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYdTd(3,3),YuadjYdYdadjYumu2(3,3),YuadjYdYdadjYuYu(3,3),& 
& YuadjYdYdadjYuTu(3,3),YuadjYdTdadjYdYd(3,3),YuadjYdTdadjYuYu(3,3),YuadjYumu2YuadjYu(3,3),& 
& YuadjYuYumq2adjYu(3,3),YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),YuadjYuYuadjYuTu(3,3),& 
& YuadjYuTuadjYuYu(3,3),YvadjYvYvadjYvvL(3),YvadjYvYvadjYvYv(3,3),YvadjYvYvadjYvTv(3,3), & 
& YvadjYvCml2YvadjYv(3,3),YvadjYvCml2TpYeCYe(3,3),YvadjYvTvadjYvYv(3,3),YvadjYvTpYeCme2CYe(3,3),& 
& YvadjYvTpYeCYeYv(3,3),YvadjYvTpYeCYeCml2(3,3),YvadjYvTpYeCYeTv(3,3),YvadjYvTpTeCYeYv(3,3),& 
& YvCmv2adjYvYvadjYv(3,3),YvCmv2adjYvTpYeCYe(3,3),YvCmv2Ckap1kap1adjYv(3,3),             & 
& YvCmv2Ckap2kap2adjYv(3,3),YvCmv2Ckap3kap3adjYv(3,3),YvCkap1TpYvCYvTpTk1(3,3),          & 
& YvCkap1Tpkap1adjYvvL(3),YvCkap1Tpkap1adjYvYv(3,3),YvCkap1Tpkap1Cmv2adjYv(3,3),         & 
& YvCkap1TpTk1adjYvYv(3,3),YvCkap1kap1adjYvCml2(3,3),YvCkap2TpYvCYvTpTk2(3,3),           & 
& YvCkap2Tpkap2adjYvvL(3),YvCkap2Tpkap2adjYvYv(3,3),YvCkap2Tpkap2Cmv2adjYv(3,3),         & 
& YvCkap2TpTk2adjYvYv(3,3),YvCkap2kap2adjYvCml2(3,3),YvCkap3TpYvCYvTpTk3(3,3)

Complex(dp) :: YvCkap3Tpkap3adjYvvL(3),YvCkap3Tpkap3adjYvYv(3,3),YvCkap3Tpkap3Cmv2adjYv(3,3),         & 
& YvCkap3TpTk3adjYvYv(3,3),YvCkap3kap3adjYvCml2(3,3),adjYdmd2YdadjYdYd(3,3),             & 
& adjYdYdmq2adjYdYd(3,3),adjYdYdadjYdmd2Yd(3,3),adjYdYdadjYdYdmq2(3,3),adjYeme2YeadjYeYe(3,3),& 
& adjYeme2YeCYvlam(3),adjYeYeml2adjYeYe(3,3),adjYeYeml2CYvlam(3),adjYeYeadjYeme2Ye(3,3), & 
& adjYeYeadjYeYeml2(3,3),adjYeYeCYvmv2lam(3),adjYumu2YuadjYuYu(3,3),adjYuYumq2adjYuYu(3,3),& 
& adjYuYuadjYumu2Yu(3,3),adjYuYuadjYuYumq2(3,3),adjYvYvadjkap1mv2kap2(3,3),              & 
& adjYvYvadjkap1mv2kap3(3,3),adjYvYvadjkap2mv2kap1(3,3),adjYvYvadjkap2mv2kap3(3,3),      & 
& adjYvYvadjkap3mv2kap1(3,3),adjYvYvadjkap3mv2kap2(3,3),adjYvYvCmv2Ckap1kap2(3,3),       & 
& adjYvYvCmv2Ckap1kap3(3,3),adjYvYvCmv2Ckap2kap1(3,3),adjYvYvCmv2Ckap2kap3(3,3),         & 
& adjYvYvCmv2Ckap3kap1(3,3),adjYvYvCmv2Ckap3kap2(3,3),adjYvCml2YvCkap1kap2(3,3),         & 
& adjYvCml2YvCkap1kap3(3,3),adjYvCml2YvCkap2kap1(3,3),adjYvCml2YvCkap2kap3(3,3),         & 
& adjYvCml2YvCkap3kap1(3,3),adjYvCml2YvCkap3kap2(3,3),CYvmv2TpYvCYvlam(3),               & 
& CYvmv2TpYvCYvTpYv(3,3),CYvmv2kap1adjkap1TpYv(3,3),CYvmv2kap2adjkap2TpYv(3,3),          & 
& CYvmv2kap3adjkap3TpYv(3,3),CYvTpYvml2CYvlam(3),CYvTpYvml2CYvTpYv(3,3),CYvTpYvCYvmv2lam(3),& 
& CYvTpYvCYvmv2TpYv(3,3),CYvTpYvCYvTpYvml2(3,3),CYvkap1adjkap1TpYvml2(3,3),              & 
& CYvkap2adjkap2TpYvml2(3,3),CYvkap3adjkap3TpYvml2(3,3),TdadjYdYdadjYdYd(3,3),           & 
& TdadjYuYuadjYdYd(3,3),TdadjYuYuadjYuYu(3,3),TeadjYeYeadjYeYe(3,3),TeCYvTpYvadjYeYe(3,3),& 
& TeCYvTpYvCYvTpYv(3,3),TuadjYdYdadjYdYd(3,3),TuadjYdYdadjYuYu(3,3),TuadjYuYuadjYuYu(3,3),& 
& TvadjYvYvadjYvYv(3,3),TvadjYvTpYeCYeYv(3,3),TpYeCYeTpYeCYevL(3),TpYeCYeTpYeCYeYv(3,3), & 
& TpYeCYeTpYeCYeTv(3,3),TpYeCYeTpTeCYeYv(3,3),TpYvml2adjYeYeCYv(3,3),TpYvml2CYvTpYvCYv(3,3),& 
& TpYvadjYeme2YeCYv(3,3),TpYvadjYeYeml2CYv(3,3),TpYvadjYeYeCYvmv2(3,3),TpYvadjYeYeCYvvR(3),& 
& TpYvadjYeYeCYvlam(3),TpYvadjYeTeCYvlam(3),TpYvCYvmv2TpYvCYv(3,3),TpYvCYvTpYvml2CYv(3,3),& 
& TpYvCYvTpYvCYvmv2(3,3),TpYvCYvTpYvCYvvR(3),TpYvCYvTpYvCYvlam(3),TpYvCYvTpYvCYvTlam(3), & 
& TpYvCYvTpTvCYvlam(3),TpYvCYvkap1adjkap1lam(3),TpYvCYvkap2adjkap2lam(3),TpYvCYvkap3adjkap3lam(3),& 
& TpYvCYvTk1adjkap1lam(3),TpYvCYvTk2adjkap2lam(3),TpYvCYvTk3adjkap3lam(3),               & 
& TpTeCYeTpYeCYeYv(3,3),TpTvadjYeYeCYvlam(3),TpTvCYvTpYvCYvlam(3),kap1adjYvYvadjkap1mv2(3,3),& 
& kap1adjYvYvCmv2adjkap1(3,3),kap1Cmv2adjYvYvadjkap1(3,3),kap1Ckap1Tpkap1Cmv2adjkap1(3,3),& 
& kap1Ckap2Tpkap2Cmv2adjkap1(3,3),kap1Ckap3Tpkap3Cmv2adjkap1(3,3),kap2adjYvYvadjkap2mv2(3,3),& 
& kap2adjYvYvCmv2adjkap2(3,3),kap2Cmv2adjYvYvadjkap2(3,3),kap2Ckap1Tpkap1Cmv2adjkap2(3,3),& 
& kap2Ckap2Tpkap2Cmv2adjkap2(3,3),kap2Ckap3Tpkap3Cmv2adjkap2(3,3),kap3adjYvYvadjkap3mv2(3,3),& 
& kap3adjYvYvCmv2adjkap3(3,3),kap3Cmv2adjYvYvadjkap3(3,3),kap3Ckap1Tpkap1Cmv2adjkap3(3,3),& 
& kap3Ckap2Tpkap2Cmv2adjkap3(3,3),kap3Ckap3Tpkap3Cmv2adjkap3(3,3),Tk1adjYvYvCkap1lam(3), & 
& Tk2adjYvYvCkap2lam(3),Tk3adjYvYvCkap3lam(3)

Complex(dp) :: Trmd2,Trme2,Trml2,Trmq2,Trmu2,TrYdadjYd,TrYeadjYe,TrYuadjYu,TrYvadjYv,TradjYdTd,      & 
& TradjYeTe,TradjYuTu,TradjYvTv,TrCTdTpTd,TrCTeTpTe,TrCTuTpTu,TrCTvTpTv,Trmd2YdadjYd,    & 
& Trme2YeadjYe,Trml2adjYeYe,Trmq2adjYdYd,Trmq2adjYuYu,Trmu2YuadjYu,TrYvadjYvCml2,        & 
& TrYvCmv2adjYv

Complex(dp) :: TradjTk1Tk1,TradjTk2Tk2,TradjTk3Tk3,TrCTdTpYd,TrCTeTpYe,TrCTuTpYu,TrCTvTpYv,          & 
& TrCkap1Tpkap1,TrCkap1Tpkap2,TrCkap1Tpkap3,TrCkap1kap1,TrCkap1kap2,TrCkap1kap3,         & 
& TrCkap1Tk1,TrCkap1Tk2,TrCkap1Tk3,TrCkap2Tpkap1,TrCkap2Tpkap2,TrCkap2Tpkap3,            & 
& TrCkap2kap1,TrCkap2kap2,TrCkap2kap3,TrCkap2Tk1,TrCkap2Tk2,TrCkap2Tk3,TrCkap3Tpkap1,    & 
& TrCkap3Tpkap2,TrCkap3Tpkap3,TrCkap3kap1,TrCkap3kap2,TrCkap3kap3,TrCkap3Tk1,            & 
& TrCkap3Tk2,TrCkap3Tk3,TrCTk1Tpkap1,TrCTk1Tpkap2,TrCTk1Tpkap3,TrCTk2Tpkap1,             & 
& TrCTk2Tpkap2,TrCTk2Tpkap3,TrCTk3Tpkap1,TrCTk3Tpkap2,TrCTk3Tpkap3,Trml2YvadjYv,         & 
& Trmv2kap1adjkap1,Trmv2kap1adjkap2,Trmv2kap1adjkap3,Trmv2kap2adjkap1,Trmv2kap2adjkap2,  & 
& Trmv2kap2adjkap3,Trmv2kap3adjkap1,Trmv2kap3adjkap2,Trmv2kap3adjkap3,TrYdadjYdCmd2,     & 
& TrYdCmq2adjYd,TrYeadjYeCme2,TrYeCml2adjYe,TrYuadjYuCmu2,TrYuCmq2adjYu,TrYdadjYdYdadjYd,& 
& TrYdadjYdTdadjYd,TrYdadjYdTdadjTd,TrYdadjYuYuadjYd,TrYdadjYuTuadjYd,TrYdadjYuTuadjTd,  & 
& TrYdadjTdTdadjYd,TrYdadjTuTuadjYd,TrYeadjYeYeadjYe,TrYeadjYeTeadjYe,TrYeadjYeTeadjTe,  & 
& TrYeadjTeTeadjYe,TrYeCTvTpTvadjYe,TrYuadjYdTdadjYu,TrYuadjYdTdadjTu,TrYuadjYuYuadjYu,  & 
& TrYuadjYuTuadjYu,TrYuadjYuTuadjTu,TrYuadjTdTdadjYu,TrYuadjTuTuadjYu,TrYvadjYvYvadjYv,  & 
& TrYvadjYvTvadjYv,TrYvadjYvTvadjTv,TrYvadjYvTpYeCYe,TrYvadjYvTpTeCTe,TrYvadjTvTvadjYv,  & 
& TrYvCkap1Tpkap1adjYv,TrYvCkap1TpTk1adjYv,TrYvCkap1kap1adjYv,TrYvCkap1Tk1adjTv,         & 
& TrYvCkap2Tpkap2adjYv,TrYvCkap2TpTk2adjYv,TrYvCkap2kap2adjYv,TrYvCkap2Tk2adjTv,         & 
& TrYvCkap3Tpkap3adjYv,TrYvCkap3TpTk3adjYv,TrYvCkap3kap3adjYv,TrYvCkap3Tk3adjTv,         & 
& TrYvCTk1TpTk1adjYv,TrYvCTk2TpTk2adjYv,TrYvCTk3TpTk3adjYv,TradjYeTeCYvTpYv,             & 
& TradjYeTeCTvTpYv,TradjYvTvadjkap1kap1,TradjYvTvadjkap2kap2,TradjYvTvadjkap3kap3,       & 
& TradjYvTpYeCYeTv,TradjYvTpYeCTeTv,Trmd2YdadjYdYdadjYd,Trmd2YdadjYuYuadjYd,             & 
& Trme2YeadjYeYeadjYe,Trml2adjYeYeadjYeYe,Trmq2adjYdYdadjYdYd,Trmq2adjYdYdadjYuYu,       & 
& Trmq2adjYuYuadjYdYd,Trmq2adjYuYuadjYuYu,Trmu2YuadjYdYdadjYu,Trmu2YuadjYuYuadjYu,       & 
& Trmv2kap1adjYvYvadjkap1,Trmv2kap2adjYvYvadjkap2,Trmv2kap3adjYvYvadjkap3,               & 
& TrYvadjYvCml2YvadjYv,TrYvadjYvCml2TpYeCYe,TrYvadjYvTpYeCme2CYe,TrYvadjYvTpYeCYeCml2,   & 
& TrYvCmv2adjYvYvadjYv,TrYvCmv2adjYvTpYeCYe,TrYvCmv2Ckap1kap1adjYv,TrYvCmv2Ckap2kap2adjYv,& 
& TrYvCmv2Ckap3kap3adjYv,TrYvCkap1Tpkap1Cmv2adjYv,TrYvCkap1kap1adjYvCml2,TrYvCkap2Tpkap2Cmv2adjYv,& 
& TrYvCkap2kap2adjYvCml2,TrYvCkap3Tpkap3Cmv2adjYv,TrYvCkap3kap3adjYvCml2

Real(dp) :: SPvRxxClam,SPlamxxClam,SPlamxxadjYvmlHd2,SPlamxxadjYvvL,SPlamxxCmv2Clam,              & 
& SPClamxxTlam,SPTlamxxCTlam

Real(dp) :: SPmlHd2xxadjYeYeCYvlam,SPvLxxadjYeYeCYvlam,SPvRxxadjYvYvClam,SPlamxxCTlam,            & 
& SPlamxxadjkap1lam,SPlamxxadjkap2lam,SPlamxxadjkap3lam,SPlamxxadjYvYvClam,              & 
& SPlamxxadjTvTvClam,SPlamxxadjYvYvadjYvmlHd2,SPlamxxadjYvYvadjYvvL,SPlamxxadjYvYvCmv2Clam,& 
& SPlamxxCmv2adjYvYvClam,SPClamxxTpYvmlHd2,SPClamxxTpTvCYvlam,SPClamxxTpYvml2CYvlam,     & 
& SPadjkap1lamxxTpkap1adjYvvL,SPadjkap2lamxxTpkap2adjYvvL,SPadjkap3lamxxTpkap3adjYvvL,   & 
& SPTpkap1Clamxxadjkap1lam,SPTpkap2Clamxxadjkap2lam,SPTpkap3Clamxxadjkap3lam,            & 
& SPTpTk1ClamxxadjTk1lam,SPTpTk2ClamxxadjTk2lam,SPTpTk3ClamxxadjTk3lam,SPadjYvYvClamxxTlam,& 
& SPadjYvYvCTlamxxTlam,SPadjTvYvClamxxTlam,SPTpTvCYvlamxxCTlam

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Complex(dp) :: sqrt3ov5,ooSqrt15

Real(dp) :: g1p4,g1p5,g2p4,g2p5,g3p4,g3p5

Complex(dp) :: Xip2,SPlamxxClamp2

Complex(dp) :: Dylami1Clami2(3,3),Dymv2lami1Clami2(3,3),DyYvClami1lami2(3,3),Dylami1adjYvmlHd2i2(3,3),& 
& Dylami1Cmv2Clami2(3,3),DyCYvlami1mlHd2i2(3,3),DyTvClami1lami2(3,3),DyTvi11kap1Ckap1i21(3,3),& 
& DyTvi11kap2Ckap1i22(3,3),DyTvi11kap3Ckap1i23(3,3),DyTvi12kap1Ckap2i21(3,3),            & 
& DyTvi12kap2Ckap2i22(3,3),DyTvi12kap3Ckap2i23(3,3),DyTvi13kap1Ckap3i21(3,3),            & 
& DyTvi13kap2Ckap3i22(3,3),DyTvi13kap3Ckap3i23(3,3),DyTlami1CTlami2(3,3),DyYvClami1Tlami2(3,3)

Complex(dp) :: Dylami1CTlami2(3,3),DyCYvlami1YvClami2(3,3),DyCYvTlami1YvCTlami2(3,3),DyCTvTlami1YvClami2(3,3),& 
& DyCTvlami1TvClami2(3,3),DyCYvlami1TvCTlami2(3,3),Dykap1Clami1Cmv2Ckap1i21(3,3),        & 
& Dykap1Clami1Cmv2Ckap2i21(3,3),Dykap1Clami1Cmv2Ckap3i21(3,3),Dykap1Clami1Ckap1mv2i21(3,3),& 
& Dykap1Clami1Ckap2mv2i21(3,3),Dykap1Clami1Ckap3mv2i21(3,3),DyYvCkap1i11kap1Clami2(3,3), & 
& DyYvCkap2i11kap1Clami2(3,3),DyYvCkap3i11kap1Clami2(3,3),Dykap2Clami1Cmv2Ckap1i22(3,3), & 
& Dykap2Clami1Cmv2Ckap2i22(3,3),Dykap2Clami1Cmv2Ckap3i22(3,3),Dykap2Clami1Ckap1mv2i22(3,3),& 
& Dykap2Clami1Ckap2mv2i22(3,3),Dykap2Clami1Ckap3mv2i22(3,3),DyYvCkap1i12kap2Clami2(3,3), & 
& DyYvCkap2i12kap2Clami2(3,3),DyYvCkap3i12kap2Clami2(3,3),Dykap3Clami1Cmv2Ckap1i23(3,3), & 
& Dykap3Clami1Cmv2Ckap2i23(3,3),Dykap3Clami1Cmv2Ckap3i23(3,3),Dykap3Clami1Ckap1mv2i23(3,3),& 
& Dykap3Clami1Ckap2mv2i23(3,3),Dykap3Clami1Ckap3mv2i23(3,3),DyYvCkap1i13kap3Clami2(3,3), & 
& DyYvCkap2i13kap3Clami2(3,3),DyYvCkap3i13kap3Clami2(3,3),DyTk1Clami1CTk1lami2(3,3),     & 
& DyTk2Clami1CTk2lami2(3,3),DyTk3Clami1CTk3lami2(3,3),Dyml2CYvlami1YvClami2(3,3),        & 
& Dymv2kap1Clami1Ckap1lami2(3,3),Dymv2kap2Clami1Ckap2lami2(3,3),Dymv2kap3Clami1Ckap3lami2(3,3),& 
& DyYeCYvlami1YvClami2(3,3),DyYeCYvlami1TvClami2(3,3),DyYeCYvTlami1YvClami2(3,3),        & 
& DyCYvlami1YvadjYvmlHd2i2(3,3),DyYvadjkap1kap1i11adjYvTv1i2(3,3),DyCYvi11Yvadjkap1kap1i21(3,3),& 
& Dyml2CYvi11Yvadjkap1kap1i21(3,3),DyTeCYvi11Yvadjkap1kap1i21(3,3),DyYvadjkap1kap2i11adjYvTv2i2(3,3),& 
& DyCYvi12Yvadjkap1kap2i21(3,3),Dyml2CYvi12Yvadjkap1kap2i21(3,3),DyTeCYvi12Yvadjkap1kap2i21(3,3),& 
& DyYvadjkap1kap3i11adjYvTv3i2(3,3),DyCYvi13Yvadjkap1kap3i21(3,3),Dyml2CYvi13Yvadjkap1kap3i21(3,3),& 
& DyTeCYvi13Yvadjkap1kap3i21(3,3),DyYvadjkap1Tk1i11lami2(3,3),DyYvadjkap1Tk2i11lami2(3,3),& 
& DyYvadjkap1Tk3i11lami2(3,3),DyYvadjkap2kap1i12adjYvTv1i2(3,3),DyCYvi11Yvadjkap2kap1i22(3,3),& 
& Dyml2CYvi11Yvadjkap2kap1i22(3,3),DyTeCYvi11Yvadjkap2kap1i22(3,3),DyYvadjkap2kap2i12adjYvTv2i2(3,3),& 
& DyCYvi12Yvadjkap2kap2i22(3,3),Dyml2CYvi12Yvadjkap2kap2i22(3,3),DyTeCYvi12Yvadjkap2kap2i22(3,3),& 
& DyYvadjkap2kap3i12adjYvTv3i2(3,3),DyCYvi13Yvadjkap2kap3i22(3,3),Dyml2CYvi13Yvadjkap2kap3i22(3,3),& 
& DyTeCYvi13Yvadjkap2kap3i22(3,3),DyYvadjkap2Tk1i12lami2(3,3),DyYvadjkap2Tk2i12lami2(3,3),& 
& DyYvadjkap2Tk3i12lami2(3,3),DyYvadjkap3kap1i13adjYvTv1i2(3,3),DyCYvi11Yvadjkap3kap1i23(3,3),& 
& Dyml2CYvi11Yvadjkap3kap1i23(3,3),DyTeCYvi11Yvadjkap3kap1i23(3,3),DyYvadjkap3kap2i13adjYvTv2i2(3,3),& 
& DyCYvi12Yvadjkap3kap2i23(3,3),Dyml2CYvi12Yvadjkap3kap2i23(3,3),DyTeCYvi12Yvadjkap3kap2i23(3,3),& 
& DyYvadjkap3kap3i13adjYvTv3i2(3,3),DyCYvi13Yvadjkap3kap3i23(3,3),Dyml2CYvi13Yvadjkap3kap3i23(3,3),& 
& DyTeCYvi13Yvadjkap3kap3i23(3,3),DyYvadjkap3Tk1i13lami2(3,3),DyYvadjkap3Tk2i13lami2(3,3),& 
& DyYvadjkap3Tk3i13lami2(3,3),DyCYvlami1YvCmv2Clami2(3,3),DyYvCkap1lami1Tk1Clami2(3,3),  & 
& DyCTvi11YvCkap1Tk1i21(3,3),DyCTvi12YvCkap1Tk2i21(3,3),DyCTvi13YvCkap1Tk3i21(3,3),      & 
& DyYvCkap2lami1Tk2Clami2(3,3),DyCTvi11YvCkap2Tk1i22(3,3),DyCTvi12YvCkap2Tk2i22(3,3),    & 
& DyCTvi13YvCkap2Tk3i22(3,3),DyYvCkap3lami1Tk3Clami2(3,3),DyCTvi11YvCkap3Tk1i23(3,3),    & 
& DyCTvi12YvCkap3Tk2i23(3,3),DyCTvi13YvCkap3Tk3i23(3,3),Dylami1adjYvYvClami2(3,3),       & 
& Dymv2lami1adjYvYvClami2(3,3),Dylami1adjYvTvCTlami2(3,3),Dylami1adjTvTvClami2(3,3),     & 
& DyCYvlami1Cml2YvClami2(3,3),DyCYvmv2lami1YvClami2(3,3),DyCYvTpkap1Ckap1i11mlHd2i2(3,3),& 
& DyCYvTpkap1Ckap2i11mlHd2i2(3,3),DyCYvTpkap1Ckap3i11mlHd2i2(3,3),DyCYvTpkap2Ckap1i12mlHd2i2(3,3)

Complex(dp) :: DyCYvTpkap2Ckap2i12mlHd2i2(3,3),DyCYvTpkap2Ckap3i12mlHd2i2(3,3),DyCYvTpkap3Ckap1i13mlHd2i2(3,3),& 
& DyCYvTpkap3Ckap2i13mlHd2i2(3,3),DyCYvTpkap3Ckap3i13mlHd2i2(3,3),Dykap1Clami1Ckap1mv2lami2(3,3),& 
& Dymv2kap1i11Ckap1kap1Ckap1i21(3,3),Dymv2kap1i12Ckap1kap1Ckap2i21(3,3),Dymv2kap1i13Ckap1kap1Ckap3i21(3,3),& 
& Dymv2kap1i11Ckap1kap2Ckap1i22(3,3),Dymv2kap1i12Ckap1kap2Ckap2i22(3,3),Dymv2kap1i13Ckap1kap2Ckap3i22(3,3),& 
& Dymv2kap1i11Ckap1kap3Ckap1i23(3,3),Dymv2kap1i12Ckap1kap3Ckap2i23(3,3),Dymv2kap1i13Ckap1kap3Ckap3i23(3,3),& 
& Dykap2Clami1Ckap2mv2lami2(3,3),Dymv2kap2i11Ckap2kap1Ckap1i21(3,3),Dymv2kap2i12Ckap2kap1Ckap2i21(3,3),& 
& Dymv2kap2i13Ckap2kap1Ckap3i21(3,3),Dymv2kap2i11Ckap2kap2Ckap1i22(3,3),Dymv2kap2i12Ckap2kap2Ckap2i22(3,3),& 
& Dymv2kap2i13Ckap2kap2Ckap3i22(3,3),Dymv2kap2i11Ckap2kap3Ckap1i23(3,3),Dymv2kap2i12Ckap2kap3Ckap2i23(3,3),& 
& Dymv2kap2i13Ckap2kap3Ckap3i23(3,3),Dykap3Clami1Ckap3mv2lami2(3,3),Dymv2kap3i11Ckap3kap1Ckap1i21(3,3),& 
& Dymv2kap3i12Ckap3kap1Ckap2i21(3,3),Dymv2kap3i13Ckap3kap1Ckap3i21(3,3),Dymv2kap3i11Ckap3kap2Ckap1i22(3,3),& 
& Dymv2kap3i12Ckap3kap2Ckap2i22(3,3),Dymv2kap3i13Ckap3kap2Ckap3i22(3,3),Dymv2kap3i11Ckap3kap3Ckap1i23(3,3),& 
& Dymv2kap3i12Ckap3kap3Ckap2i23(3,3),Dymv2kap3i13Ckap3kap3Ckap3i23(3,3),DyTeCYvlami1YvClami2(3,3),& 
& DyTvCkap1kap1i11TpYvCYvi21(3,3),DyYeCYvi11TvCkap1kap1i21(3,3),DyTvCkap1kap2i11TpYvCYvi22(3,3),& 
& DyYeCYvi12TvCkap1kap2i21(3,3),DyTvCkap1kap3i11TpYvCYvi23(3,3),DyYeCYvi13TvCkap1kap3i21(3,3),& 
& DyTvCkap2kap1i12TpYvCYvi21(3,3),DyYeCYvi11TvCkap2kap1i22(3,3),DyTvCkap2kap2i12TpYvCYvi22(3,3),& 
& DyYeCYvi12TvCkap2kap2i22(3,3),DyTvCkap2kap3i12TpYvCYvi23(3,3),DyYeCYvi13TvCkap2kap3i22(3,3),& 
& DyTvCkap3kap1i13TpYvCYvi21(3,3),DyYeCYvi11TvCkap3kap1i23(3,3),DyTvCkap3kap2i13TpYvCYvi22(3,3),& 
& DyYeCYvi12TvCkap3kap2i23(3,3),DyTvCkap3kap3i13TpYvCYvi23(3,3),DyYeCYvi13TvCkap3kap3i23(3,3),& 
& DyTpYvCYvlami1Clami2(3,3),DyTpYvCYvlami1adjYvmlHd2i2(3,3),DyTpYvCYvlami1Cmv2Clami2(3,3),& 
& DyYvClami1TpYvCYvlami2(3,3),DyTvClami1TpYvCYvlami2(3,3),DyTpYvCYvTlami1CTlami2(3,3),   & 
& DyYvClami1TpYvCYvTlami2(3,3),DyTpYvCYvTpkap11i1Ckap1mv2i21(3,3),DyYvCkap1i11TpYvCYvTpkap11i2(3,3),& 
& DyTpYvCYvTpkap12i1Ckap2mv2i21(3,3),DyYvCkap2i11TpYvCYvTpkap12i2(3,3),DyTpYvCYvTpkap13i1Ckap3mv2i21(3,3),& 
& DyYvCkap3i11TpYvCYvTpkap13i2(3,3),DyTpYvCYvTpkap21i1Ckap1mv2i22(3,3),DyYvCkap1i12TpYvCYvTpkap21i2(3,3),& 
& DyTpYvCYvTpkap22i1Ckap2mv2i22(3,3),DyYvCkap2i12TpYvCYvTpkap22i2(3,3),DyTpYvCYvTpkap23i1Ckap3mv2i22(3,3),& 
& DyYvCkap3i12TpYvCYvTpkap23i2(3,3),DyTpYvCYvTpkap31i1Ckap1mv2i23(3,3),DyYvCkap1i13TpYvCYvTpkap31i2(3,3),& 
& DyTpYvCYvTpkap32i1Ckap2mv2i23(3,3),DyYvCkap2i13TpYvCYvTpkap32i2(3,3),DyTpYvCYvTpkap33i1Ckap3mv2i23(3,3),& 
& DyYvCkap3i13TpYvCYvTpkap33i2(3,3),DyTpYvCTvTlami1Clami2(3,3),DyTpTvCYvlami1CTlami2(3,3),& 
& DyYvClami1TpTvCYvlami2(3,3),DyTpTvCTvlami1Clami2(3,3),DyTpkap1Ckap1TpTv1i1lami2(3,3),  & 
& DyTpkap1Ckap2TpTv2i1lami2(3,3),DyTpkap1Ckap3TpTv3i1lami2(3,3),DyTpkap2Ckap1TpTv1i1lami2(3,3),& 
& DyTpkap2Ckap2TpTv2i1lami2(3,3),DyTpkap2Ckap3TpTv3i1lami2(3,3),DyTpkap3Ckap1TpTv1i1lami2(3,3),& 
& DyTpkap3Ckap2TpTv2i1lami2(3,3),DyTpkap3Ckap3TpTv3i1lami2(3,3),DyYvCkap1i11kap1adjYvTvi21(3,3),& 
& DyYvCkap2i11kap1adjYvTvi22(3,3),DyYvCkap3i11kap1adjYvTvi23(3,3),DyYvCkap1i11kap1adjkap1Tpkap11i2(3,3),& 
& DyYvCkap2i11kap1adjkap1Tpkap12i2(3,3),DyYvCkap3i11kap1adjkap1Tpkap13i2(3,3),           & 
& DyYvCkap1i12kap1adjkap1Tpkap21i2(3,3),DyYvCkap2i12kap1adjkap1Tpkap22i2(3,3),           & 
& DyYvCkap3i12kap1adjkap1Tpkap23i2(3,3),DyYvCkap1i13kap1adjkap1Tpkap31i2(3,3),           & 
& DyYvCkap2i13kap1adjkap1Tpkap32i2(3,3),DyYvCkap3i13kap1adjkap1Tpkap33i2(3,3),           & 
& Dykap1Cmv2Clami1Ckap1i21(3,3),Dykap1Cmv2Clami1Ckap2i21(3,3),Dykap1Cmv2Clami1Ckap3i21(3,3)

Complex(dp) :: Dykap1Ckap1kap1i11Cmv2Ckap1i21(3,3),Dykap1Ckap1kap1i11Ckap1mv2i21(3,3),Dykap1Ckap1kap2i11Cmv2Ckap2i21(3,3),& 
& Dykap1Ckap1kap2i11Ckap2mv2i21(3,3),Dykap1Ckap1kap3i11Cmv2Ckap3i21(3,3),Dykap1Ckap1kap3i11Ckap3mv2i21(3,3),& 
& DyYvCkap1i11kap1Ckap1Tk1i21(3,3),DyYvCkap2i11kap1Ckap1Tk2i21(3,3),DyYvCkap3i11kap1Ckap1Tk3i21(3,3),& 
& Dykap1Ckap2kap1i12Cmv2Ckap1i21(3,3),Dykap1Ckap2kap1i12Ckap1mv2i21(3,3),Dykap1Ckap2kap2i12Cmv2Ckap2i21(3,3),& 
& Dykap1Ckap2kap2i12Ckap2mv2i21(3,3),Dykap1Ckap2kap3i12Cmv2Ckap3i21(3,3),Dykap1Ckap2kap3i12Ckap3mv2i21(3,3),& 
& DyYvCkap1i11kap1Ckap2Tk1i22(3,3),DyYvCkap2i11kap1Ckap2Tk2i22(3,3),DyYvCkap3i11kap1Ckap2Tk3i22(3,3),& 
& Dykap1Ckap3kap1i13Cmv2Ckap1i21(3,3),Dykap1Ckap3kap1i13Ckap1mv2i21(3,3),Dykap1Ckap3kap2i13Cmv2Ckap2i21(3,3),& 
& Dykap1Ckap3kap2i13Ckap2mv2i21(3,3),Dykap1Ckap3kap3i13Cmv2Ckap3i21(3,3),Dykap1Ckap3kap3i13Ckap3mv2i21(3,3),& 
& DyYvCkap1i11kap1Ckap3Tk1i23(3,3),DyYvCkap2i11kap1Ckap3Tk2i23(3,3),DyYvCkap3i11kap1Ckap3Tk3i23(3,3),& 
& DyYvCkap1i12kap2adjYvTvi21(3,3),DyYvCkap2i12kap2adjYvTvi22(3,3),DyYvCkap3i12kap2adjYvTvi23(3,3),& 
& DyYvCkap1i11kap2adjkap2Tpkap11i2(3,3),DyYvCkap2i11kap2adjkap2Tpkap12i2(3,3),           & 
& DyYvCkap3i11kap2adjkap2Tpkap13i2(3,3),DyYvCkap1i12kap2adjkap2Tpkap21i2(3,3),           & 
& DyYvCkap2i12kap2adjkap2Tpkap22i2(3,3),DyYvCkap3i12kap2adjkap2Tpkap23i2(3,3),           & 
& DyYvCkap1i13kap2adjkap2Tpkap31i2(3,3),DyYvCkap2i13kap2adjkap2Tpkap32i2(3,3),           & 
& DyYvCkap3i13kap2adjkap2Tpkap33i2(3,3),Dykap2Cmv2Clami1Ckap1i22(3,3),Dykap2Cmv2Clami1Ckap2i22(3,3),& 
& Dykap2Cmv2Clami1Ckap3i22(3,3),Dykap2Ckap1kap1i11Cmv2Ckap1i22(3,3),Dykap2Ckap1kap1i11Ckap1mv2i22(3,3),& 
& Dykap2Ckap1kap2i11Cmv2Ckap2i22(3,3),Dykap2Ckap1kap2i11Ckap2mv2i22(3,3),Dykap2Ckap1kap3i11Cmv2Ckap3i22(3,3),& 
& Dykap2Ckap1kap3i11Ckap3mv2i22(3,3),DyYvCkap1i12kap2Ckap1Tk1i21(3,3),DyYvCkap2i12kap2Ckap1Tk2i21(3,3),& 
& DyYvCkap3i12kap2Ckap1Tk3i21(3,3),Dykap2Ckap2kap1i12Cmv2Ckap1i22(3,3),Dykap2Ckap2kap1i12Ckap1mv2i22(3,3),& 
& Dykap2Ckap2kap2i12Cmv2Ckap2i22(3,3),Dykap2Ckap2kap2i12Ckap2mv2i22(3,3),Dykap2Ckap2kap3i12Cmv2Ckap3i22(3,3),& 
& Dykap2Ckap2kap3i12Ckap3mv2i22(3,3),DyYvCkap1i12kap2Ckap2Tk1i22(3,3),DyYvCkap2i12kap2Ckap2Tk2i22(3,3),& 
& DyYvCkap3i12kap2Ckap2Tk3i22(3,3),Dykap2Ckap3kap1i13Cmv2Ckap1i22(3,3),Dykap2Ckap3kap1i13Ckap1mv2i22(3,3),& 
& Dykap2Ckap3kap2i13Cmv2Ckap2i22(3,3),Dykap2Ckap3kap2i13Ckap2mv2i22(3,3),Dykap2Ckap3kap3i13Cmv2Ckap3i22(3,3),& 
& Dykap2Ckap3kap3i13Ckap3mv2i22(3,3),DyYvCkap1i12kap2Ckap3Tk1i23(3,3),DyYvCkap2i12kap2Ckap3Tk2i23(3,3),& 
& DyYvCkap3i12kap2Ckap3Tk3i23(3,3),DyYvCkap1i13kap3adjYvTvi21(3,3),DyYvCkap2i13kap3adjYvTvi22(3,3),& 
& DyYvCkap3i13kap3adjYvTvi23(3,3),DyYvCkap1i11kap3adjkap3Tpkap11i2(3,3),DyYvCkap2i11kap3adjkap3Tpkap12i2(3,3),& 
& DyYvCkap3i11kap3adjkap3Tpkap13i2(3,3),DyYvCkap1i12kap3adjkap3Tpkap21i2(3,3),           & 
& DyYvCkap2i12kap3adjkap3Tpkap22i2(3,3),DyYvCkap3i12kap3adjkap3Tpkap23i2(3,3),           & 
& DyYvCkap1i13kap3adjkap3Tpkap31i2(3,3),DyYvCkap2i13kap3adjkap3Tpkap32i2(3,3),           & 
& DyYvCkap3i13kap3adjkap3Tpkap33i2(3,3),Dykap3Cmv2Clami1Ckap1i23(3,3),Dykap3Cmv2Clami1Ckap2i23(3,3),& 
& Dykap3Cmv2Clami1Ckap3i23(3,3),Dykap3Ckap1kap1i11Cmv2Ckap1i23(3,3),Dykap3Ckap1kap1i11Ckap1mv2i23(3,3),& 
& Dykap3Ckap1kap2i11Cmv2Ckap2i23(3,3),Dykap3Ckap1kap2i11Ckap2mv2i23(3,3),Dykap3Ckap1kap3i11Cmv2Ckap3i23(3,3),& 
& Dykap3Ckap1kap3i11Ckap3mv2i23(3,3),DyYvCkap1i13kap3Ckap1Tk1i21(3,3),DyYvCkap2i13kap3Ckap1Tk2i21(3,3),& 
& DyYvCkap3i13kap3Ckap1Tk3i21(3,3),Dykap3Ckap2kap1i12Cmv2Ckap1i23(3,3),Dykap3Ckap2kap1i12Ckap1mv2i23(3,3),& 
& Dykap3Ckap2kap2i12Cmv2Ckap2i23(3,3),Dykap3Ckap2kap2i12Ckap2mv2i23(3,3),Dykap3Ckap2kap3i12Cmv2Ckap3i23(3,3),& 
& Dykap3Ckap2kap3i12Ckap3mv2i23(3,3),DyYvCkap1i13kap3Ckap2Tk1i22(3,3),DyYvCkap2i13kap3Ckap2Tk2i22(3,3),& 
& DyYvCkap3i13kap3Ckap2Tk3i22(3,3),Dykap3Ckap3kap1i13Cmv2Ckap1i23(3,3),Dykap3Ckap3kap1i13Ckap1mv2i23(3,3)

Complex(dp) :: Dykap3Ckap3kap2i13Cmv2Ckap2i23(3,3),Dykap3Ckap3kap2i13Ckap2mv2i23(3,3),Dykap3Ckap3kap3i13Cmv2Ckap3i23(3,3),& 
& Dykap3Ckap3kap3i13Ckap3mv2i23(3,3),DyYvCkap1i13kap3Ckap3Tk1i23(3,3),DyYvCkap2i13kap3Ckap3Tk2i23(3,3),& 
& DyYvCkap3i13kap3Ckap3Tk3i23(3,3),DyYvCkap1i11Tk1Ckap1kap1i21(3,3),DyYvCkap1i12Tk1Ckap1kap2i21(3,3),& 
& DyYvCkap1i13Tk1Ckap1kap3i21(3,3),DyYvCkap1i11Tk1Ckap2kap1i22(3,3),DyYvCkap1i12Tk1Ckap2kap2i22(3,3),& 
& DyYvCkap1i13Tk1Ckap2kap3i22(3,3),DyYvCkap1i11Tk1Ckap3kap1i23(3,3),DyYvCkap1i12Tk1Ckap3kap2i23(3,3),& 
& DyYvCkap1i13Tk1Ckap3kap3i23(3,3),DyYvCkap2i11Tk2Ckap1kap1i21(3,3),DyYvCkap2i12Tk2Ckap1kap2i21(3,3),& 
& DyYvCkap2i13Tk2Ckap1kap3i21(3,3),DyYvCkap2i11Tk2Ckap2kap1i22(3,3),DyYvCkap2i12Tk2Ckap2kap2i22(3,3),& 
& DyYvCkap2i13Tk2Ckap2kap3i22(3,3),DyYvCkap2i11Tk2Ckap3kap1i23(3,3),DyYvCkap2i12Tk2Ckap3kap2i23(3,3),& 
& DyYvCkap2i13Tk2Ckap3kap3i23(3,3),DyYvCkap3i11Tk3Ckap1kap1i21(3,3),DyYvCkap3i12Tk3Ckap1kap2i21(3,3),& 
& DyYvCkap3i13Tk3Ckap1kap3i21(3,3),DyYvCkap3i11Tk3Ckap2kap1i22(3,3),DyYvCkap3i12Tk3Ckap2kap2i22(3,3),& 
& DyYvCkap3i13Tk3Ckap2kap3i22(3,3),DyYvCkap3i11Tk3Ckap3kap1i23(3,3),DyYvCkap3i12Tk3Ckap3kap2i23(3,3),& 
& DyYvCkap3i13Tk3Ckap3kap3i23(3,3),Dymv2TpYvCYvlami1Clami2(3,3),DyYvadjYvYvClami1lami2(3,3),& 
& DyYvadjYvTvClami1lami2(3,3),DyCYvi11Yvadjkap1mv2kap1i21(3,3),DyCYvi12Yvadjkap1mv2kap2i21(3,3),& 
& DyCYvi13Yvadjkap1mv2kap3i21(3,3),DyCYvi11Yvadjkap2mv2kap1i22(3,3),DyCYvi12Yvadjkap2mv2kap2i22(3,3),& 
& DyCYvi13Yvadjkap2mv2kap3i22(3,3),DyCYvi11Yvadjkap3mv2kap1i23(3,3),DyCYvi12Yvadjkap3mv2kap2i23(3,3),& 
& DyCYvi13Yvadjkap3mv2kap3i23(3,3),DyCYvi11YvCmv2Ckap1kap1i21(3,3),DyCYvi12YvCmv2Ckap1kap2i21(3,3),& 
& DyCYvi13YvCmv2Ckap1kap3i21(3,3),DyCYvi11YvCmv2Ckap2kap1i22(3,3),DyCYvi12YvCmv2Ckap2kap2i22(3,3),& 
& DyCYvi13YvCmv2Ckap2kap3i22(3,3),DyCYvi11YvCmv2Ckap3kap1i23(3,3),DyCYvi12YvCmv2Ckap3kap2i23(3,3),& 
& DyCYvi13YvCmv2Ckap3kap3i23(3,3),DyYvCkap1Tpkap1Clami1lami2(3,3),DyYvCkap2Tpkap2Clami1lami2(3,3),& 
& DyYvCkap3Tpkap3Clami1lami2(3,3),DyadjYeYeCYvlami1mlHd2i2(3,3),Dylami1adjYvYvadjYvmlHd2i2(3,3),& 
& Dylami1adjYvYvCmv2Clami2(3,3),Dylami1adjYvCml2YvClami2(3,3),Dylami1Cmv2adjYvYvClami2(3,3),& 
& DyCYvTpYvCYvlami1mlHd2i2(3,3),DyTvadjYvYvClami1lami2(3,3),DyTpYeCYeYvClami1lami2(3,3), & 
& DyTpYeCYeTvClami1lami2(3,3),DyTpYvml2CYvlami1Clami2(3,3),DyTpYvCYvmv2lami1Clami2(3,3), & 
& DyTpTeCYeYvClami1lami2(3,3),Dykap1Cmv2Ckap1kap1i11Ckap1i21(3,3),Dykap1Cmv2Ckap1kap2i11Ckap2i21(3,3),& 
& Dykap1Cmv2Ckap1kap3i11Ckap3i21(3,3),Dykap1Cmv2Ckap2kap1i12Ckap1i21(3,3),               & 
& Dykap1Cmv2Ckap2kap2i12Ckap2i21(3,3),Dykap1Cmv2Ckap2kap3i12Ckap3i21(3,3),               & 
& Dykap1Cmv2Ckap3kap1i13Ckap1i21(3,3),Dykap1Cmv2Ckap3kap2i13Ckap2i21(3,3),               & 
& Dykap1Cmv2Ckap3kap3i13Ckap3i21(3,3),Dykap2Cmv2Ckap1kap1i11Ckap1i22(3,3),               & 
& Dykap2Cmv2Ckap1kap2i11Ckap2i22(3,3),Dykap2Cmv2Ckap1kap3i11Ckap3i22(3,3),               & 
& Dykap2Cmv2Ckap2kap1i12Ckap1i22(3,3),Dykap2Cmv2Ckap2kap2i12Ckap2i22(3,3),               & 
& Dykap2Cmv2Ckap2kap3i12Ckap3i22(3,3),Dykap2Cmv2Ckap3kap1i13Ckap1i22(3,3),               & 
& Dykap2Cmv2Ckap3kap2i13Ckap2i22(3,3),Dykap2Cmv2Ckap3kap3i13Ckap3i22(3,3),               & 
& Dykap3Cmv2Ckap1kap1i11Ckap1i23(3,3),Dykap3Cmv2Ckap1kap2i11Ckap2i23(3,3),               & 
& Dykap3Cmv2Ckap1kap3i11Ckap3i23(3,3),Dykap3Cmv2Ckap2kap1i12Ckap1i23(3,3),               & 
& Dykap3Cmv2Ckap2kap2i12Ckap2i23(3,3),Dykap3Cmv2Ckap2kap3i12Ckap3i23(3,3),               & 
& Dykap3Cmv2Ckap3kap1i13Ckap1i23(3,3),Dykap3Cmv2Ckap3kap2i13Ckap2i23(3,3),               & 
& Dykap3Cmv2Ckap3kap3i13Ckap3i23(3,3),DyTvi11kap1Ckap1i22(3,3),DyTvi11kap1Ckap1i23(3,3)

Complex(dp) :: DyTvi11kap2Ckap1i21(3,3),DyTvi11kap2Ckap1i23(3,3),DyTvi11kap3Ckap1i21(3,3),            & 
& DyTvi11kap3Ckap1i22(3,3),DyTvi11kap1Ckap1lami2(3,3),DyTvi11kap2Ckap1lami2(3,3),        & 
& DyTvi11kap3Ckap1lami2(3,3),DyTvi11kap1adjYvYvCkap1i21(3,3),DyTvi11kap2adjYvYvCkap1i22(3,3),& 
& DyTvi11kap3adjYvYvCkap1i23(3,3),DyTvi12kap1Ckap2i22(3,3),DyTvi12kap1Ckap2i23(3,3),     & 
& DyTvi12kap2Ckap2i21(3,3),DyTvi12kap2Ckap2i23(3,3),DyTvi12kap3Ckap2i21(3,3),            & 
& DyTvi12kap3Ckap2i22(3,3),DyTvi12kap1Ckap2lami2(3,3),DyTvi12kap2Ckap2lami2(3,3),        & 
& DyTvi12kap3Ckap2lami2(3,3),DyTvi12kap1adjYvYvCkap2i21(3,3),DyTvi12kap2adjYvYvCkap2i22(3,3),& 
& DyTvi12kap3adjYvYvCkap2i23(3,3),DyTvi13kap1Ckap3i22(3,3),DyTvi13kap1Ckap3i23(3,3),     & 
& DyTvi13kap2Ckap3i21(3,3),DyTvi13kap2Ckap3i23(3,3),DyTvi13kap3Ckap3i21(3,3),            & 
& DyTvi13kap3Ckap3i22(3,3),DyTvi13kap1Ckap3lami2(3,3),DyTvi13kap2Ckap3lami2(3,3),        & 
& DyTvi13kap3Ckap3lami2(3,3),DyTvi13kap1adjYvYvCkap3i21(3,3),DyTvi13kap2adjYvYvCkap3i22(3,3),& 
& DyTvi13kap3adjYvYvCkap3i23(3,3),DyCYvi11Tvi21(3,3),DyCYvi12Tvi21(3,3),DyCYvi13Tvi21(3,3),& 
& DyCTvi11Tvi21(3,3),DyCTvi12Tvi21(3,3),DyCTvi13Tvi21(3,3),DyCYvi11Tvi22(3,3),           & 
& DyCYvi12Tvi22(3,3),DyCYvi13Tvi22(3,3),DyCTvi11Tvi22(3,3),DyCTvi12Tvi22(3,3),           & 
& DyCTvi13Tvi22(3,3),DyCYvi11Tvi23(3,3),DyCYvi12Tvi23(3,3),DyCYvi13Tvi23(3,3),           & 
& DyCTvi11Tvi23(3,3),DyCTvi12Tvi23(3,3),DyCTvi13Tvi23(3,3),DyTlami1Clami2(3,3),          & 
& DyTlami1adjYvYvCTlami2(3,3),DyTlami1adjTvYvClami2(3,3),DyYvadjkap1kap1i11Tlami2(3,3),  & 
& DyYvadjkap1kap2i11Tlami2(3,3),DyYvadjkap1kap3i11Tlami2(3,3),DyYvadjkap2kap1i12Tlami2(3,3),& 
& DyYvadjkap2kap2i12Tlami2(3,3),DyYvadjkap2kap3i12Tlami2(3,3),DyYvadjkap3kap1i13Tlami2(3,3),& 
& DyYvadjkap3kap2i13Tlami2(3,3),DyYvadjkap3kap3i13Tlami2(3,3),DyYvadjYvYvClami1Tlami2(3,3),& 
& DyTpYeCYeYvClami1Tlami2(3,3),Dykap1i11Ckap1TpYvCYvi21(3,3),Dykap1i11adjYvCml2Yvadjkap11i2(3,3),& 
& Dykap1i12Ckap2TpYvCYvi21(3,3),Dykap1i12adjYvCml2Yvadjkap21i2(3,3),Dykap1i13Ckap3TpYvCYvi21(3,3),& 
& Dykap1i13adjYvCml2Yvadjkap31i2(3,3),Dykap2i11Ckap1TpYvCYvi22(3,3),Dykap2i11adjYvCml2Yvadjkap12i2(3,3),& 
& Dykap2i12Ckap2TpYvCYvi22(3,3),Dykap2i12adjYvCml2Yvadjkap22i2(3,3),Dykap2i13Ckap3TpYvCYvi22(3,3),& 
& Dykap2i13adjYvCml2Yvadjkap32i2(3,3),Dykap3i11Ckap1TpYvCYvi23(3,3),Dykap3i11adjYvCml2Yvadjkap13i2(3,3),& 
& Dykap3i12Ckap2TpYvCYvi23(3,3),Dykap3i12adjYvCml2Yvadjkap23i2(3,3),Dykap3i13Ckap3TpYvCYvi23(3,3),& 
& Dykap3i13adjYvCml2Yvadjkap33i2(3,3),DyTk1i11Ckap1i21(3,3),DyTk1i11Ckap1i22(3,3),       & 
& DyTk1i11Ckap1i23(3,3),DyTk1i11CTk1i21(3,3),DyTk1i11CTk1i22(3,3),DyTk1i11CTk1i23(3,3),  & 
& DyTk1i11Ckap1lami2(3,3),DyTk1i11Ckap1TpYvCTvi21(3,3),DyTk1i12Ckap2i21(3,3),            & 
& DyTk1i12Ckap2i22(3,3),DyTk1i12Ckap2i23(3,3),DyTk1i12CTk1i21(3,3),DyTk1i12CTk1i22(3,3), & 
& DyTk1i12CTk1i23(3,3),DyTk1i12Ckap2lami2(3,3),DyTk1i12Ckap2TpYvCTvi21(3,3),             & 
& DyTk1i13Ckap3i21(3,3),DyTk1i13Ckap3i22(3,3),DyTk1i13Ckap3i23(3,3),DyTk1i13CTk1i21(3,3),& 
& DyTk1i13CTk1i22(3,3),DyTk1i13CTk1i23(3,3),DyTk1i13Ckap3lami2(3,3),DyTk1i13Ckap3TpYvCTvi21(3,3),& 
& DyTk2i11Ckap1i21(3,3),DyTk2i11Ckap1i22(3,3),DyTk2i11Ckap1i23(3,3),DyTk2i11CTk2i21(3,3),& 
& DyTk2i11CTk2i22(3,3),DyTk2i11CTk2i23(3,3),DyTk2i11Ckap1lami2(3,3),DyTk2i11Ckap1TpYvCTvi22(3,3),& 
& DyTk2i12Ckap2i21(3,3),DyTk2i12Ckap2i22(3,3),DyTk2i12Ckap2i23(3,3),DyTk2i12CTk2i21(3,3),& 
& DyTk2i12CTk2i22(3,3),DyTk2i12CTk2i23(3,3),DyTk2i12Ckap2lami2(3,3),DyTk2i12Ckap2TpYvCTvi22(3,3)

Complex(dp) :: DyTk2i13Ckap3i21(3,3),DyTk2i13Ckap3i22(3,3),DyTk2i13Ckap3i23(3,3),DyTk2i13CTk2i21(3,3),& 
& DyTk2i13CTk2i22(3,3),DyTk2i13CTk2i23(3,3),DyTk2i13Ckap3lami2(3,3),DyTk2i13Ckap3TpYvCTvi22(3,3),& 
& DyTk3i11Ckap1i21(3,3),DyTk3i11Ckap1i22(3,3),DyTk3i11Ckap1i23(3,3),DyTk3i11CTk3i21(3,3),& 
& DyTk3i11CTk3i22(3,3),DyTk3i11CTk3i23(3,3),DyTk3i11Ckap1lami2(3,3),DyTk3i11Ckap1TpYvCTvi23(3,3),& 
& DyTk3i12Ckap2i21(3,3),DyTk3i12Ckap2i22(3,3),DyTk3i12Ckap2i23(3,3),DyTk3i12CTk3i21(3,3),& 
& DyTk3i12CTk3i22(3,3),DyTk3i12CTk3i23(3,3),DyTk3i12Ckap2lami2(3,3),DyTk3i12Ckap2TpYvCTvi23(3,3),& 
& DyTk3i13Ckap3i21(3,3),DyTk3i13Ckap3i22(3,3),DyTk3i13Ckap3i23(3,3),DyTk3i13CTk3i21(3,3),& 
& DyTk3i13CTk3i22(3,3),DyTk3i13CTk3i23(3,3),DyTk3i13Ckap3lami2(3,3),DyTk3i13Ckap3TpYvCTvi23(3,3)

Complex(dp) :: kap1(3,3), adjkap1(3,3) 
Complex(dp) :: kap2(3,3), adjkap2(3,3) 
Complex(dp) :: kap3(3,3), adjkap3(3,3) 
Complex(dp) :: Tk1(3,3), adjTk1(3,3) 
Complex(dp) :: Tk2(3,3), adjTk2(3,3) 
Complex(dp) :: Tk3(3,3), adjTk3(3,3) 
Iname = Iname +1 
NameOfUnit(Iname) = 'rge394' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters394(gy,g1,g2,g3,Yd,Ye,lam,Yv,Yu,kap,Td,Te,Tlam,Tv,Tu,               & 
& Tk,mq2,ml2,mHd2,mHu2,md2,mu2,me2,mv2,mlHd2,M1,M2,M3,vd,vu,vL,vR)

kap1=kap(:,:,1) 
Call Adjungate(kap1,adjkap1) 
kap2=kap(:,:,2) 
Call Adjungate(kap2,adjkap2) 
kap3=kap(:,:,3) 
Call Adjungate(kap3,adjkap3) 
Tk1=Tk(:,:,1) 
Call Adjungate(Tk1,adjTk1) 
Tk2=Tk(:,:,2) 
Call Adjungate(Tk2,adjTk2) 
Tk3=Tk(:,:,3) 
Call Adjungate(Tk3,adjTk3) 
AbsM1 = Conjg(M1)*M1
AbsM2 = Conjg(M2)*M2
AbsM3 = Conjg(M3)*M3
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yv,adjYv)
Call Adjungate(Yu,adjYu)
Call Adjungate(Td,adjTd)
Call Adjungate(Te,adjTe)
Call Adjungate(Tv,adjTv)
Call Adjungate(Tu,adjTu)
 md2Yd = Matmul(md2,Yd) 
 me2Ye = Matmul(me2,Ye) 
 ml2adjYe = Matmul(ml2,adjYe) 
 ml2CYv = Matmul(ml2,Conjg(Yv)) 
 mq2adjYd = Matmul(mq2,adjYd) 
 mq2adjYu = Matmul(mq2,adjYu) 
 mu2Yu = Matmul(mu2,Yu) 
 mv2lam = Matmul(mv2,lam) 
 mv2TpYv = Matmul(mv2,Transpose(Yv)) 
 Ydmq2 = Matmul(Yd,mq2) 
 YdadjYd = Matmul(Yd,adjYd) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 Yeml2 = Matmul(Ye,ml2) 
 YeadjYe = Matmul(Ye,adjYe) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 Yumq2 = Matmul(Yu,mq2) 
 YuadjYu = Matmul(Yu,adjYu) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 YvadjYv = Matmul(Yv,adjYv) 
Forall(i2=1:3)  YvadjYv(i2,i2) =  Real(YvadjYv(i2,i2),dp) 
 YvClam = Matmul(Yv,Conjg(lam)) 
 adjYdmd2 = Matmul(adjYd,md2) 
 adjYdYd = Matmul(adjYd,Yd) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdTd = Matmul(adjYd,Td) 
 adjYeme2 = Matmul(adjYe,me2) 
 adjYeYe = Matmul(adjYe,Ye) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYeTe = Matmul(adjYe,Te) 
 adjYumu2 = Matmul(adjYu,mu2) 
 adjYuYu = Matmul(adjYu,Yu) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYuTu = Matmul(adjYu,Tu) 
 adjYvmlHd2 = Matmul(adjYv,mlHd2) 
 adjYvvL = Matmul(adjYv,vL) 
 adjYvYv = Matmul(adjYv,Yv) 
Forall(i2=1:3)  adjYvYv(i2,i2) =  Real(adjYvYv(i2,i2),dp) 
 adjYvCml2 = Matmul(adjYv,Conjg(ml2)) 
 adjYvTv = Matmul(adjYv,Tv) 
 adjTdTd = Matmul(adjTd,Td) 
 adjTeTe = Matmul(adjTe,Te) 
 adjTuTu = Matmul(adjTu,Tu) 
 adjkap1mv2 = Matmul(adjkap1,mv2) 
 adjkap1lam = Matmul(adjkap1,lam) 
 adjkap2mv2 = Matmul(adjkap2,mv2) 
 adjkap2lam = Matmul(adjkap2,lam) 
 adjkap3mv2 = Matmul(adjkap3,mv2) 
 adjkap3lam = Matmul(adjkap3,lam) 
 Cmv2adjYv = Matmul(Conjg(mv2),adjYv) 
 Cmv2adjkap1 = Matmul(Conjg(mv2),adjkap1) 
 Cmv2adjkap2 = Matmul(Conjg(mv2),adjkap2) 
 Cmv2adjkap3 = Matmul(Conjg(mv2),adjkap3) 
 Cmv2Clam = Matmul(Conjg(mv2),Conjg(lam)) 
 CYevL = Matmul(Conjg(Ye),vL) 
 CYeYv = Matmul(Conjg(Ye),Yv) 
 CYeTv = Matmul(Conjg(Ye),Tv) 
 CYvmv2 = Matmul(Conjg(Yv),mv2) 
 CYvvR = Matmul(Conjg(Yv),vR) 
 CYvlam = Matmul(Conjg(Yv),lam) 
 CYvTlam = Matmul(Conjg(Yv),Tlam) 
 CYvTpYv = Matmul(Conjg(Yv),Transpose(Yv)) 
Forall(i2=1:3)  CYvTpYv(i2,i2) =  Real(CYvTpYv(i2,i2),dp) 
 CYvTpTv = Matmul(Conjg(Yv),Transpose(Tv)) 
 CTdTpTd = Matmul(Conjg(Td),Transpose(Td)) 
Forall(i2=1:3)  CTdTpTd(i2,i2) =  Real(CTdTpTd(i2,i2),dp) 
 CTeTpTe = Matmul(Conjg(Te),Transpose(Te)) 
Forall(i2=1:3)  CTeTpTe(i2,i2) =  Real(CTeTpTe(i2,i2),dp) 
 CTuTpTu = Matmul(Conjg(Tu),Transpose(Tu)) 
Forall(i2=1:3)  CTuTpTu(i2,i2) =  Real(CTuTpTu(i2,i2),dp) 
 CTvTlam = Matmul(Conjg(Tv),Tlam) 
 CTvTpTv = Matmul(Conjg(Tv),Transpose(Tv)) 
Forall(i2=1:3)  CTvTpTv(i2,i2) =  Real(CTvTpTv(i2,i2),dp) 
 Ckap1vR = Matmul(Conjg(kap1),vR) 
 Ckap1Tpkap1 = Matmul(Conjg(kap1),Transpose(kap1)) 
Forall(i2=1:3)  Ckap1Tpkap1(i2,i2) =  Real(Ckap1Tpkap1(i2,i2),dp) 
 Ckap1TpTk1 = Matmul(Conjg(kap1),Transpose(Tk1)) 
 Ckap2vR = Matmul(Conjg(kap2),vR) 
 Ckap2Tpkap2 = Matmul(Conjg(kap2),Transpose(kap2)) 
Forall(i2=1:3)  Ckap2Tpkap2(i2,i2) =  Real(Ckap2Tpkap2(i2,i2),dp) 
 Ckap2TpTk2 = Matmul(Conjg(kap2),Transpose(Tk2)) 
 Ckap3vR = Matmul(Conjg(kap3),vR) 
 Ckap3Tpkap3 = Matmul(Conjg(kap3),Transpose(kap3)) 
Forall(i2=1:3)  Ckap3Tpkap3(i2,i2) =  Real(Ckap3Tpkap3(i2,i2),dp) 
 Ckap3TpTk3 = Matmul(Conjg(kap3),Transpose(Tk3)) 
 TdadjTd = Matmul(Td,adjTd) 
 TeadjTe = Matmul(Te,adjTe) 
 TuadjTu = Matmul(Tu,adjTu) 
 TvClam = Matmul(Tv,Conjg(lam)) 
 TpYvml2 = Matmul(Transpose(Yv),ml2) 
 TpYvCYv = Matmul(Transpose(Yv),Conjg(Yv)) 
Forall(i2=1:3)  TpYvCYv(i2,i2) =  Real(TpYvCYv(i2,i2),dp) 
 TpTvCTv = Matmul(Transpose(Tv),Conjg(Tv)) 
Forall(i2=1:3)  TpTvCTv(i2,i2) =  Real(TpTvCTv(i2,i2),dp) 
 kap1adjkap1 = Matmul(kap1,adjkap1) 
 kap1Ckap1 = Matmul(kap1,Conjg(kap1)) 
 kap1Ckap2 = Matmul(kap1,Conjg(kap2)) 
 kap1Ckap3 = Matmul(kap1,Conjg(kap3)) 
 kap2adjkap2 = Matmul(kap2,adjkap2) 
 kap2Ckap1 = Matmul(kap2,Conjg(kap1)) 
 kap2Ckap2 = Matmul(kap2,Conjg(kap2)) 
 kap2Ckap3 = Matmul(kap2,Conjg(kap3)) 
 kap3adjkap3 = Matmul(kap3,adjkap3) 
 kap3Ckap1 = Matmul(kap3,Conjg(kap1)) 
 kap3Ckap2 = Matmul(kap3,Conjg(kap2)) 
 kap3Ckap3 = Matmul(kap3,Conjg(kap3)) 
 Tk1adjTk1 = Matmul(Tk1,adjTk1) 
 Tk2adjTk2 = Matmul(Tk2,adjTk2) 
 Tk3adjTk3 = Matmul(Tk3,adjTk3) 
 md2YdadjYd = Matmul(md2,YdadjYd) 
 me2YeadjYe = Matmul(me2,YeadjYe) 
 ml2adjYeYe = Matmul(ml2,adjYeYe) 
 ml2CYvlam = Matmul(ml2,CYvlam) 
 ml2CYvTpYv = Matmul(ml2,CYvTpYv) 
 mq2adjYdYd = Matmul(mq2,adjYdYd) 
 mq2adjYuYu = Matmul(mq2,adjYuYu) 
 mu2YuadjYu = Matmul(mu2,YuadjYu) 
 mv2TpYvCYv = Matmul(mv2,TpYvCYv) 
 mv2kap1adjkap1 = Matmul(mv2,kap1adjkap1) 
 mv2kap2adjkap2 = Matmul(mv2,kap2adjkap2) 
 mv2kap3adjkap3 = Matmul(mv2,kap3adjkap3) 
 Ydmq2adjYd = Matmul(Yd,mq2adjYd) 
Forall(i2=1:3)  Ydmq2adjYd(i2,i2) =  Real(Ydmq2adjYd(i2,i2),dp) 
 YdadjYdmd2 = Matmul(Yd,adjYdmd2) 
 YdadjYdYd = Matmul(Yd,adjYdYd) 
 YdadjYdTd = Matmul(Yd,adjYdTd) 
 YdadjYuYu = Matmul(Yd,adjYuYu) 
 YdadjYuTu = Matmul(Yd,adjYuTu) 
 Yeml2adjYe = Matmul(Ye,ml2adjYe) 
Forall(i2=1:3)  Yeml2adjYe(i2,i2) =  Real(Yeml2adjYe(i2,i2),dp) 
 YeadjYeme2 = Matmul(Ye,adjYeme2) 
 YeadjYeYe = Matmul(Ye,adjYeYe) 
 YeadjYeTe = Matmul(Ye,adjYeTe) 
 YeCYvTpYv = Matmul(Ye,CYvTpYv) 
 YeCYvTpTv = Matmul(Ye,CYvTpTv) 
 Yumq2adjYu = Matmul(Yu,mq2adjYu) 
Forall(i2=1:3)  Yumq2adjYu(i2,i2) =  Real(Yumq2adjYu(i2,i2),dp) 
 YuadjYdYd = Matmul(Yu,adjYdYd) 
 YuadjYdTd = Matmul(Yu,adjYdTd) 
 YuadjYumu2 = Matmul(Yu,adjYumu2) 
 YuadjYuYu = Matmul(Yu,adjYuYu) 
 YuadjYuTu = Matmul(Yu,adjYuTu) 
 YvadjYvvL = Matmul(Yv,adjYvvL) 
 YvadjYvYv = Matmul(Yv,adjYvYv) 
 YvadjYvCml2 = Matmul(Yv,adjYvCml2) 
 YvadjYvTv = Matmul(Yv,adjYvTv) 
 YvCmv2adjYv = Matmul(Yv,Cmv2adjYv) 
Forall(i2=1:3)  YvCmv2adjYv(i2,i2) =  Real(YvCmv2adjYv(i2,i2),dp) 
 YvCkap1Tpkap1 = Matmul(Yv,Ckap1Tpkap1) 
 YvCkap1TpTk1 = Matmul(Yv,Ckap1TpTk1) 
 YvCkap2Tpkap2 = Matmul(Yv,Ckap2Tpkap2) 
 YvCkap2TpTk2 = Matmul(Yv,Ckap2TpTk2) 
 YvCkap3Tpkap3 = Matmul(Yv,Ckap3Tpkap3) 
 YvCkap3TpTk3 = Matmul(Yv,Ckap3TpTk3) 
 adjYdmd2Yd = Matmul(adjYd,md2Yd) 
Forall(i2=1:3)  adjYdmd2Yd(i2,i2) =  Real(adjYdmd2Yd(i2,i2),dp) 
 adjYdYdmq2 = Matmul(adjYd,Ydmq2) 
 adjYeme2Ye = Matmul(adjYe,me2Ye) 
Forall(i2=1:3)  adjYeme2Ye(i2,i2) =  Real(adjYeme2Ye(i2,i2),dp) 
 adjYeYeml2 = Matmul(adjYe,Yeml2) 
 adjYumu2Yu = Matmul(adjYu,mu2Yu) 
Forall(i2=1:3)  adjYumu2Yu(i2,i2) =  Real(adjYumu2Yu(i2,i2),dp) 
 adjYuYumq2 = Matmul(adjYu,Yumq2) 
 CYvmv2lam = Matmul(Conjg(Yv),mv2lam) 
 CYvmv2TpYv = Matmul(Conjg(Yv),mv2TpYv) 
Forall(i2=1:3)  CYvmv2TpYv(i2,i2) =  Real(CYvmv2TpYv(i2,i2),dp) 
 CYvTpYvml2 = Matmul(Conjg(Yv),TpYvml2) 
 TdadjYdYd = Matmul(Td,adjYdYd) 
 TdadjYuYu = Matmul(Td,adjYuYu) 
 TeadjYeYe = Matmul(Te,adjYeYe) 
 TeCYvTpYv = Matmul(Te,CYvTpYv) 
 TuadjYdYd = Matmul(Tu,adjYdYd) 
 TuadjYuYu = Matmul(Tu,adjYuYu) 
 TvadjYvYv = Matmul(Tv,adjYvYv) 
 TpYeCYevL = Matmul(Transpose(Ye),CYevL) 
 TpYeCYeYv = Matmul(Transpose(Ye),CYeYv) 
 TpYeCYeTv = Matmul(Transpose(Ye),CYeTv) 
 TpYvml2CYv = Matmul(Transpose(Yv),ml2CYv) 
Forall(i2=1:3)  TpYvml2CYv(i2,i2) =  Real(TpYvml2CYv(i2,i2),dp) 
 TpYvCYvmv2 = Matmul(Transpose(Yv),CYvmv2) 
 TpYvCYvvR = Matmul(Transpose(Yv),CYvvR) 
 TpYvCYvlam = Matmul(Transpose(Yv),CYvlam) 
 TpYvCYvTlam = Matmul(Transpose(Yv),CYvTlam) 
 TpTeCYeYv = Matmul(Transpose(Te),CYeYv) 
 TpTvCYvlam = Matmul(Transpose(Tv),CYvlam) 
 kap1adjkap1mv2 = Matmul(kap1,adjkap1mv2) 
 kap1Cmv2adjkap1 = Matmul(kap1,Cmv2adjkap1) 
 kap1Ckap1vR = Matmul(kap1,Ckap1vR) 
 kap2adjkap2mv2 = Matmul(kap2,adjkap2mv2) 
 kap2Cmv2adjkap2 = Matmul(kap2,Cmv2adjkap2) 
 kap2Ckap2vR = Matmul(kap2,Ckap2vR) 
 kap3adjkap3mv2 = Matmul(kap3,adjkap3mv2) 
 kap3Cmv2adjkap3 = Matmul(kap3,Cmv2adjkap3) 
 kap3Ckap3vR = Matmul(kap3,Ckap3vR) 
 Tk1adjkap1lam = Matmul(Tk1,adjkap1lam) 
 Tk2adjkap2lam = Matmul(Tk2,adjkap2lam) 
 Tk3adjkap3lam = Matmul(Tk3,adjkap3lam) 
 Trmd2 = Real(cTrace(md2),dp) 
 Trme2 = Real(cTrace(me2),dp) 
 Trml2 = Real(cTrace(ml2),dp) 
 Trmq2 = Real(cTrace(mq2),dp) 
 Trmu2 = Real(cTrace(mu2),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYvadjYv = Real(cTrace(YvadjYv),dp) 
 TradjYdTd = Real(cTrace(adjYdTd),dp) 
 TradjYeTe = Real(cTrace(adjYeTe),dp) 
 TradjYuTu = Real(cTrace(adjYuTu),dp) 
 TradjYvTv = Real(cTrace(adjYvTv),dp) 
 TrCTdTpTd = Real(cTrace(CTdTpTd),dp) 
 TrCTeTpTe = Real(cTrace(CTeTpTe),dp) 
 TrCTuTpTu = Real(cTrace(CTuTpTu),dp) 
 TrCTvTpTv = Real(cTrace(CTvTpTv),dp) 
 Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp) 
 Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp) 
 Trml2adjYeYe = Real(cTrace(ml2adjYeYe),dp) 
 Trmq2adjYdYd = Real(cTrace(mq2adjYdYd),dp) 
 Trmq2adjYuYu = Real(cTrace(mq2adjYuYu),dp) 
 Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp) 
 TrYvadjYvCml2 = Real(cTrace(YvadjYvCml2),dp) 
 TrYvCmv2adjYv = Real(cTrace(YvCmv2adjYv),dp) 
 SPvRxxClam = DOT_PRODUCT(vR,Conjg(lam)) 
 SPlamxxClam = DOT_PRODUCT(Conjg(lam),Conjg(lam)) 
 SPlamxxadjYvmlHd2 = DOT_PRODUCT(Conjg(lam),adjYvmlHd2) 
 SPlamxxadjYvvL = DOT_PRODUCT(Conjg(lam),adjYvvL) 
 SPlamxxCmv2Clam = DOT_PRODUCT(Conjg(lam),Cmv2Clam) 
 SPClamxxTlam = DOT_PRODUCT(Conjg(Conjg(lam)),Tlam) 
 SPTlamxxCTlam = DOT_PRODUCT(Conjg(Tlam),Conjg(Tlam)) 
 sqrt3ov5 =Sqrt(3._dp/5._dp) 
 ooSqrt15 =1._dp/sqrt(15._dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 
 g1p4 =g1**4 
 g1p5 =g1**5 
 g2p4 =g2**4 
 g2p5 =g2**5 
 g3p4 =g3**4 
 g3p5 =g3**5 
 Xip2 =Xi**2 
 SPlamxxClamp2 =SPlamxxClam**2 
Do i1=1,3
  Do i2=1,3
Dylami1Clami2(i1,i2) = Conjg(lam(i2))*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2lami1Clami2(i1,i2) = Conjg(lam(i2))*mv2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1lami2(i1,i2) = YvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvmlHd2i2(i1,i2) = adjYvmlHd2(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1Cmv2Clami2(i1,i2) = Cmv2Clam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1mlHd2i2(i1,i2) = CYvlam(i1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvClami1lami2(i1,i2) = TvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1i21(i1,i2) = kap1Ckap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1i22(i1,i2) = kap2Ckap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1i23(i1,i2) = kap3Ckap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2i21(i1,i2) = kap1Ckap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2i22(i1,i2) = kap2Ckap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2i23(i1,i2) = kap3Ckap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3i21(i1,i2) = kap1Ckap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3i22(i1,i2) = kap2Ckap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3i23(i1,i2) = kap3Ckap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1CTlami2(i1,i2) = Conjg(Tlam(i2))*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1Tlami2(i1,i2) = YvClam(i1)*Tlam(i2) 
  End Do 
End Do 


If (TwoLoopRGE) Then 
 mv2kap1 = Matmul(mv2,kap1) 
 mv2kap2 = Matmul(mv2,kap2) 
 mv2kap3 = Matmul(mv2,kap3) 
 YdadjYu = Matmul(Yd,adjYu) 
 YdadjTd = Matmul(Yd,adjTd) 
 YdadjTu = Matmul(Yd,adjTu) 
 YevL = Matmul(Ye,vL) 
 YeadjTe = Matmul(Ye,adjTe) 
 YeCYv = Matmul(Ye,Conjg(Yv)) 
 YeCTv = Matmul(Ye,Conjg(Tv)) 
 YuadjYd = Matmul(Yu,adjYd) 
 YuadjTd = Matmul(Yu,adjTd) 
 YuadjTu = Matmul(Yu,adjTu) 
 Yvadjkap1 = Matmul(Yv,adjkap1) 
 Yvadjkap2 = Matmul(Yv,adjkap2) 
 Yvadjkap3 = Matmul(Yv,adjkap3) 
 YvadjTk1 = Matmul(Yv,adjTk1) 
 YvadjTk2 = Matmul(Yv,adjTk2) 
 YvadjTk3 = Matmul(Yv,adjTk3) 
 YvCTlam = Matmul(Yv,Conjg(Tlam)) 
 YvCkap1 = Matmul(Yv,Conjg(kap1)) 
 YvCkap2 = Matmul(Yv,Conjg(kap2)) 
 YvCkap3 = Matmul(Yv,Conjg(kap3)) 
 adjYdCmd2 = Matmul(adjYd,Conjg(md2)) 
 adjYeCme2 = Matmul(adjYe,Conjg(me2)) 
 adjYuCmu2 = Matmul(adjYu,Conjg(mu2)) 
 adjTdYd = Matmul(adjTd,Yd) 
 adjTeYe = Matmul(adjTe,Ye) 
 adjTuYu = Matmul(adjTu,Yu) 
 adjkap1vR = Matmul(adjkap1,vR) 
 adjkap1Tlam = Matmul(adjkap1,Tlam) 
 adjkap1TpYv = Matmul(adjkap1,Transpose(Yv)) 
 adjkap1Tpkap1 = Matmul(adjkap1,Transpose(kap1)) 
 adjkap1Tpkap2 = Matmul(adjkap1,Transpose(kap2)) 
 adjkap1Tpkap3 = Matmul(adjkap1,Transpose(kap3)) 
 adjkap1kap1 = Matmul(adjkap1,kap1) 
 adjkap1kap2 = Matmul(adjkap1,kap2) 
 adjkap1kap3 = Matmul(adjkap1,kap3) 
 adjkap1Tk1 = Matmul(adjkap1,Tk1) 
 adjkap1Tk2 = Matmul(adjkap1,Tk2) 
 adjkap1Tk3 = Matmul(adjkap1,Tk3) 
 adjkap2vR = Matmul(adjkap2,vR) 
 adjkap2Tlam = Matmul(adjkap2,Tlam) 
 adjkap2TpYv = Matmul(adjkap2,Transpose(Yv)) 
 adjkap2Tpkap1 = Matmul(adjkap2,Transpose(kap1)) 
 adjkap2Tpkap2 = Matmul(adjkap2,Transpose(kap2)) 
 adjkap2Tpkap3 = Matmul(adjkap2,Transpose(kap3)) 
 adjkap2kap1 = Matmul(adjkap2,kap1) 
 adjkap2kap2 = Matmul(adjkap2,kap2) 
 adjkap2kap3 = Matmul(adjkap2,kap3) 
 adjkap2Tk1 = Matmul(adjkap2,Tk1) 
 adjkap2Tk2 = Matmul(adjkap2,Tk2) 
 adjkap2Tk3 = Matmul(adjkap2,Tk3) 
 adjkap3vR = Matmul(adjkap3,vR) 
 adjkap3Tlam = Matmul(adjkap3,Tlam) 
 adjkap3TpYv = Matmul(adjkap3,Transpose(Yv)) 
 adjkap3Tpkap1 = Matmul(adjkap3,Transpose(kap1)) 
 adjkap3Tpkap2 = Matmul(adjkap3,Transpose(kap2)) 
 adjkap3Tpkap3 = Matmul(adjkap3,Transpose(kap3)) 
 adjkap3kap1 = Matmul(adjkap3,kap1) 
 adjkap3kap2 = Matmul(adjkap3,kap2) 
 adjkap3kap3 = Matmul(adjkap3,kap3) 
 adjkap3Tk1 = Matmul(adjkap3,Tk1) 
 adjkap3Tk2 = Matmul(adjkap3,Tk2) 
 adjkap3Tk3 = Matmul(adjkap3,Tk3) 
 adjTk1lam = Matmul(adjTk1,lam) 
 adjTk1TpYv = Matmul(adjTk1,Transpose(Yv)) 
 adjTk1Tk1 = Matmul(adjTk1,Tk1) 
 adjTk2lam = Matmul(adjTk2,lam) 
 adjTk2TpYv = Matmul(adjTk2,Transpose(Yv)) 
 adjTk2Tk2 = Matmul(adjTk2,Tk2) 
 adjTk3lam = Matmul(adjTk3,lam) 
 adjTk3TpYv = Matmul(adjTk3,Transpose(Yv)) 
 adjTk3Tk3 = Matmul(adjTk3,Tk3) 
 Cme2CYe = Matmul(Conjg(me2),Conjg(Ye)) 
 Cml2adjYe = Matmul(Conjg(ml2),adjYe) 
 Cmq2adjYd = Matmul(Conjg(mq2),adjYd) 
 Cmq2adjYu = Matmul(Conjg(mq2),adjYu) 
 Cmv2Ckap1 = Matmul(Conjg(mv2),Conjg(kap1)) 
 Cmv2Ckap2 = Matmul(Conjg(mv2),Conjg(kap2)) 
 Cmv2Ckap3 = Matmul(Conjg(mv2),Conjg(kap3)) 
 CYeCml2 = Matmul(Conjg(Ye),Conjg(ml2)) 
 CYvTpkap1 = Matmul(Conjg(Yv),Transpose(kap1)) 
 CYvTpkap2 = Matmul(Conjg(Yv),Transpose(kap2)) 
 CYvTpkap3 = Matmul(Conjg(Yv),Transpose(kap3)) 
 CYvTpTk1 = Matmul(Conjg(Yv),Transpose(Tk1)) 
 CYvTpTk2 = Matmul(Conjg(Yv),Transpose(Tk2)) 
 CYvTpTk3 = Matmul(Conjg(Yv),Transpose(Tk3)) 
 CTdTpYd = Matmul(Conjg(Td),Transpose(Yd)) 
 CTeTv = Matmul(Conjg(Te),Tv) 
 CTeTpYe = Matmul(Conjg(Te),Transpose(Ye)) 
 CTuTpYu = Matmul(Conjg(Tu),Transpose(Yu)) 
 CTvlam = Matmul(Conjg(Tv),lam) 
 CTvTpYv = Matmul(Conjg(Tv),Transpose(Yv)) 
 Ckap1mv2 = Matmul(Conjg(kap1),mv2) 
 Ckap1lam = Matmul(Conjg(kap1),lam) 
 Ckap1Tlam = Matmul(Conjg(kap1),Tlam) 
 Ckap1TpTv = Matmul(Conjg(kap1),Transpose(Tv)) 
 Ckap1Tpkap2 = Matmul(Conjg(kap1),Transpose(kap2)) 
 Ckap1Tpkap3 = Matmul(Conjg(kap1),Transpose(kap3)) 
 Ckap1TpTk2 = Matmul(Conjg(kap1),Transpose(Tk2)) 
 Ckap1TpTk3 = Matmul(Conjg(kap1),Transpose(Tk3)) 
 Ckap1kap1 = Matmul(Conjg(kap1),kap1) 
 Ckap1kap2 = Matmul(Conjg(kap1),kap2) 
 Ckap1kap3 = Matmul(Conjg(kap1),kap3) 
 Ckap1Tk1 = Matmul(Conjg(kap1),Tk1) 
 Ckap1Tk2 = Matmul(Conjg(kap1),Tk2) 
 Ckap1Tk3 = Matmul(Conjg(kap1),Tk3) 
 Ckap2mv2 = Matmul(Conjg(kap2),mv2) 
 Ckap2lam = Matmul(Conjg(kap2),lam) 
 Ckap2Tlam = Matmul(Conjg(kap2),Tlam) 
 Ckap2TpTv = Matmul(Conjg(kap2),Transpose(Tv)) 
 Ckap2Tpkap1 = Matmul(Conjg(kap2),Transpose(kap1)) 
 Ckap2Tpkap3 = Matmul(Conjg(kap2),Transpose(kap3)) 
 Ckap2TpTk1 = Matmul(Conjg(kap2),Transpose(Tk1)) 
 Ckap2TpTk3 = Matmul(Conjg(kap2),Transpose(Tk3)) 
 Ckap2kap1 = Matmul(Conjg(kap2),kap1) 
 Ckap2kap2 = Matmul(Conjg(kap2),kap2) 
 Ckap2kap3 = Matmul(Conjg(kap2),kap3) 
 Ckap2Tk1 = Matmul(Conjg(kap2),Tk1) 
 Ckap2Tk2 = Matmul(Conjg(kap2),Tk2) 
 Ckap2Tk3 = Matmul(Conjg(kap2),Tk3) 
 Ckap3mv2 = Matmul(Conjg(kap3),mv2) 
 Ckap3lam = Matmul(Conjg(kap3),lam) 
 Ckap3Tlam = Matmul(Conjg(kap3),Tlam) 
 Ckap3TpTv = Matmul(Conjg(kap3),Transpose(Tv)) 
 Ckap3Tpkap1 = Matmul(Conjg(kap3),Transpose(kap1)) 
 Ckap3Tpkap2 = Matmul(Conjg(kap3),Transpose(kap2)) 
 Ckap3TpTk1 = Matmul(Conjg(kap3),Transpose(Tk1)) 
 Ckap3TpTk2 = Matmul(Conjg(kap3),Transpose(Tk2)) 
 Ckap3kap1 = Matmul(Conjg(kap3),kap1) 
 Ckap3kap2 = Matmul(Conjg(kap3),kap2) 
 Ckap3kap3 = Matmul(Conjg(kap3),kap3) 
 Ckap3Tk1 = Matmul(Conjg(kap3),Tk1) 
 Ckap3Tk2 = Matmul(Conjg(kap3),Tk2) 
 Ckap3Tk3 = Matmul(Conjg(kap3),Tk3) 
 CTk1lam = Matmul(Conjg(Tk1),lam) 
 CTk1Tpkap1 = Matmul(Conjg(Tk1),Transpose(kap1)) 
 CTk1Tpkap2 = Matmul(Conjg(Tk1),Transpose(kap2)) 
 CTk1Tpkap3 = Matmul(Conjg(Tk1),Transpose(kap3)) 
 CTk2lam = Matmul(Conjg(Tk2),lam) 
 CTk2Tpkap1 = Matmul(Conjg(Tk2),Transpose(kap1)) 
 CTk2Tpkap2 = Matmul(Conjg(Tk2),Transpose(kap2)) 
 CTk2Tpkap3 = Matmul(Conjg(Tk2),Transpose(kap3)) 
 CTk3lam = Matmul(Conjg(Tk3),lam) 
 CTk3Tpkap1 = Matmul(Conjg(Tk3),Transpose(kap1)) 
 CTk3Tpkap2 = Matmul(Conjg(Tk3),Transpose(kap2)) 
 CTk3Tpkap3 = Matmul(Conjg(Tk3),Transpose(kap3)) 
 TdadjYd = Matmul(Td,adjYd) 
 TdadjYu = Matmul(Td,adjYu) 
 TdadjTu = Matmul(Td,adjTu) 
 TeadjYe = Matmul(Te,adjYe) 
 TeCYv = Matmul(Te,Conjg(Yv)) 
 TeCTv = Matmul(Te,Conjg(Tv)) 
 TuadjYd = Matmul(Tu,adjYd) 
 TuadjYu = Matmul(Tu,adjYu) 
 TuadjTd = Matmul(Tu,adjTd) 
 TvadjYv = Matmul(Tv,adjYv) 
 TvadjTv = Matmul(Tv,adjTv) 
 TvCTlam = Matmul(Tv,Conjg(Tlam)) 
 TpYeCYe = Matmul(Transpose(Ye),Conjg(Ye)) 
Forall(i2=1:3)  TpYeCYe(i2,i2) =  Real(TpYeCYe(i2,i2),dp) 
 TpYvmlHd2 = Matmul(Transpose(Yv),mlHd2) 
 TpYvvL = Matmul(Transpose(Yv),vL) 
 TpYvadjYe = Matmul(Transpose(Yv),adjYe) 
 TpYvadjTe = Matmul(Transpose(Yv),adjTe) 
 TpYvCTv = Matmul(Transpose(Yv),Conjg(Tv)) 
 TpTeCTe = Matmul(Transpose(Te),Conjg(Te)) 
Forall(i2=1:3)  TpTeCTe(i2,i2) =  Real(TpTeCTe(i2,i2),dp) 
 TpTvadjYe = Matmul(Transpose(Tv),adjYe) 
 TpTvadjTe = Matmul(Transpose(Tv),adjTe) 
 TpTvCYv = Matmul(Transpose(Tv),Conjg(Yv)) 
 Tpkap1adjYv = Matmul(Transpose(kap1),adjYv) 
 Tpkap1Clam = Matmul(Transpose(kap1),Conjg(lam)) 
 Tpkap1Ckap1 = Matmul(Transpose(kap1),Conjg(kap1)) 
Forall(i2=1:3)  Tpkap1Ckap1(i2,i2) =  Real(Tpkap1Ckap1(i2,i2),dp) 
 Tpkap1Ckap2 = Matmul(Transpose(kap1),Conjg(kap2)) 
 Tpkap1Ckap3 = Matmul(Transpose(kap1),Conjg(kap3)) 
 Tpkap2adjYv = Matmul(Transpose(kap2),adjYv) 
 Tpkap2Clam = Matmul(Transpose(kap2),Conjg(lam)) 
 Tpkap2Ckap1 = Matmul(Transpose(kap2),Conjg(kap1)) 
 Tpkap2Ckap2 = Matmul(Transpose(kap2),Conjg(kap2)) 
Forall(i2=1:3)  Tpkap2Ckap2(i2,i2) =  Real(Tpkap2Ckap2(i2,i2),dp) 
 Tpkap2Ckap3 = Matmul(Transpose(kap2),Conjg(kap3)) 
 Tpkap3adjYv = Matmul(Transpose(kap3),adjYv) 
 Tpkap3Clam = Matmul(Transpose(kap3),Conjg(lam)) 
 Tpkap3Ckap1 = Matmul(Transpose(kap3),Conjg(kap1)) 
 Tpkap3Ckap2 = Matmul(Transpose(kap3),Conjg(kap2)) 
 Tpkap3Ckap3 = Matmul(Transpose(kap3),Conjg(kap3)) 
Forall(i2=1:3)  Tpkap3Ckap3(i2,i2) =  Real(Tpkap3Ckap3(i2,i2),dp) 
 TpTk1adjYv = Matmul(Transpose(Tk1),adjYv) 
 TpTk1Clam = Matmul(Transpose(Tk1),Conjg(lam)) 
 TpTk2adjYv = Matmul(Transpose(Tk2),adjYv) 
 TpTk2Clam = Matmul(Transpose(Tk2),Conjg(lam)) 
 TpTk3adjYv = Matmul(Transpose(Tk3),adjYv) 
 TpTk3Clam = Matmul(Transpose(Tk3),Conjg(lam)) 
 kap1adjYv = Matmul(kap1,adjYv) 
 kap1adjkap2 = Matmul(kap1,adjkap2) 
 kap1adjkap3 = Matmul(kap1,adjkap3) 
 kap1adjTk1 = Matmul(kap1,adjTk1) 
 kap1adjTk2 = Matmul(kap1,adjTk2) 
 kap1adjTk3 = Matmul(kap1,adjTk3) 
 kap1Clam = Matmul(kap1,Conjg(lam)) 
 kap2adjYv = Matmul(kap2,adjYv) 
 kap2adjkap1 = Matmul(kap2,adjkap1) 
 kap2adjkap3 = Matmul(kap2,adjkap3) 
 kap2adjTk1 = Matmul(kap2,adjTk1) 
 kap2adjTk2 = Matmul(kap2,adjTk2) 
 kap2adjTk3 = Matmul(kap2,adjTk3) 
 kap2Clam = Matmul(kap2,Conjg(lam)) 
 kap3adjYv = Matmul(kap3,adjYv) 
 kap3adjkap1 = Matmul(kap3,adjkap1) 
 kap3adjkap2 = Matmul(kap3,adjkap2) 
 kap3adjTk1 = Matmul(kap3,adjTk1) 
 kap3adjTk2 = Matmul(kap3,adjTk2) 
 kap3adjTk3 = Matmul(kap3,adjTk3) 
 kap3Clam = Matmul(kap3,Conjg(lam)) 
 Tk1adjTv = Matmul(Tk1,adjTv) 
 Tk1Clam = Matmul(Tk1,Conjg(lam)) 
 Tk2adjTv = Matmul(Tk2,adjTv) 
 Tk2Clam = Matmul(Tk2,Conjg(lam)) 
 Tk3adjTv = Matmul(Tk3,adjTv) 
 Tk3Clam = Matmul(Tk3,Conjg(lam)) 
 md2YdadjYu = Matmul(md2,YdadjYu) 
 me2YeCYv = Matmul(me2,YeCYv) 
 ml2YvadjYv = Matmul(ml2,YvadjYv) 
 mu2YuadjYd = Matmul(mu2,YuadjYd) 
 mv2TpYvadjYe = Matmul(mv2,TpYvadjYe) 
 mv2kap1adjkap2 = Matmul(mv2,kap1adjkap2) 
 mv2kap1adjkap3 = Matmul(mv2,kap1adjkap3) 
 mv2kap1Clam = Matmul(mv2,kap1Clam) 
 mv2kap1Ckap1 = Matmul(mv2,kap1Ckap1) 
 mv2kap1Ckap2 = Matmul(mv2,kap1Ckap2) 
 mv2kap1Ckap3 = Matmul(mv2,kap1Ckap3) 
 mv2kap2adjkap1 = Matmul(mv2,kap2adjkap1) 
 mv2kap2adjkap3 = Matmul(mv2,kap2adjkap3) 
 mv2kap2Clam = Matmul(mv2,kap2Clam) 
 mv2kap2Ckap1 = Matmul(mv2,kap2Ckap1) 
 mv2kap2Ckap2 = Matmul(mv2,kap2Ckap2) 
 mv2kap2Ckap3 = Matmul(mv2,kap2Ckap3) 
 mv2kap3adjkap1 = Matmul(mv2,kap3adjkap1) 
 mv2kap3adjkap2 = Matmul(mv2,kap3adjkap2) 
 mv2kap3Clam = Matmul(mv2,kap3Clam) 
 mv2kap3Ckap1 = Matmul(mv2,kap3Ckap1) 
 mv2kap3Ckap2 = Matmul(mv2,kap3Ckap2) 
 mv2kap3Ckap3 = Matmul(mv2,kap3Ckap3) 
 Ydmq2adjYu = Matmul(Yd,mq2adjYu) 
 YdadjYdCmd2 = Matmul(Yd,adjYdCmd2) 
 YdadjYumu2 = Matmul(Yd,adjYumu2) 
 YdadjTdTd = Matmul(Yd,adjTdTd) 
 YdCmq2adjYd = Matmul(Yd,Cmq2adjYd) 
Forall(i2=1:3)  YdCmq2adjYd(i2,i2) =  Real(YdCmq2adjYd(i2,i2),dp) 
 Yeml2CYv = Matmul(Ye,ml2CYv) 
 YeadjYeCme2 = Matmul(Ye,adjYeCme2) 
 YeadjTeTe = Matmul(Ye,adjTeTe) 
 YeCml2adjYe = Matmul(Ye,Cml2adjYe) 
Forall(i2=1:3)  YeCml2adjYe(i2,i2) =  Real(YeCml2adjYe(i2,i2),dp) 
 YeCYvmv2 = Matmul(Ye,CYvmv2) 
 YeCYvvR = Matmul(Ye,CYvvR) 
 YeCYvlam = Matmul(Ye,CYvlam) 
 YeCYvTlam = Matmul(Ye,CYvTlam) 
 YeCTvTlam = Matmul(Ye,CTvTlam) 
 Yumq2adjYd = Matmul(Yu,mq2adjYd) 
 YuadjYdmd2 = Matmul(Yu,adjYdmd2) 
 YuadjYuCmu2 = Matmul(Yu,adjYuCmu2) 
 YuadjTuTu = Matmul(Yu,adjTuTu) 
 YuCmq2adjYu = Matmul(Yu,Cmq2adjYu) 
Forall(i2=1:3)  YuCmq2adjYu(i2,i2) =  Real(YuCmq2adjYu(i2,i2),dp) 
 YvadjYvmlHd2 = Matmul(Yv,adjYvmlHd2) 
 Yvadjkap1mv2 = Matmul(Yv,adjkap1mv2) 
 Yvadjkap1kap1 = Matmul(Yv,adjkap1kap1) 
 Yvadjkap1kap2 = Matmul(Yv,adjkap1kap2) 
 Yvadjkap1kap3 = Matmul(Yv,adjkap1kap3) 
 Yvadjkap1Tk1 = Matmul(Yv,adjkap1Tk1) 
 Yvadjkap1Tk2 = Matmul(Yv,adjkap1Tk2) 
 Yvadjkap1Tk3 = Matmul(Yv,adjkap1Tk3) 
 Yvadjkap2mv2 = Matmul(Yv,adjkap2mv2) 
 Yvadjkap2kap1 = Matmul(Yv,adjkap2kap1) 
 Yvadjkap2kap2 = Matmul(Yv,adjkap2kap2) 
 Yvadjkap2kap3 = Matmul(Yv,adjkap2kap3) 
 Yvadjkap2Tk1 = Matmul(Yv,adjkap2Tk1) 
 Yvadjkap2Tk2 = Matmul(Yv,adjkap2Tk2) 
 Yvadjkap2Tk3 = Matmul(Yv,adjkap2Tk3) 
 Yvadjkap3mv2 = Matmul(Yv,adjkap3mv2) 
 Yvadjkap3kap1 = Matmul(Yv,adjkap3kap1) 
 Yvadjkap3kap2 = Matmul(Yv,adjkap3kap2) 
 Yvadjkap3kap3 = Matmul(Yv,adjkap3kap3) 
 Yvadjkap3Tk1 = Matmul(Yv,adjkap3Tk1) 
 Yvadjkap3Tk2 = Matmul(Yv,adjkap3Tk2) 
 Yvadjkap3Tk3 = Matmul(Yv,adjkap3Tk3) 
 YvCmv2adjkap1 = Matmul(Yv,Cmv2adjkap1) 
 YvCmv2adjkap2 = Matmul(Yv,Cmv2adjkap2) 
 YvCmv2adjkap3 = Matmul(Yv,Cmv2adjkap3) 
 YvCmv2Clam = Matmul(Yv,Cmv2Clam) 
 YvCkap1lam = Matmul(Yv,Ckap1lam) 
 YvCkap1kap2 = Matmul(Yv,Ckap1kap2) 
 YvCkap1kap3 = Matmul(Yv,Ckap1kap3) 
 YvCkap1Tk1 = Matmul(Yv,Ckap1Tk1) 
 YvCkap1Tk2 = Matmul(Yv,Ckap1Tk2) 
 YvCkap1Tk3 = Matmul(Yv,Ckap1Tk3) 
 YvCkap2lam = Matmul(Yv,Ckap2lam) 
 YvCkap2kap1 = Matmul(Yv,Ckap2kap1) 
 YvCkap2kap3 = Matmul(Yv,Ckap2kap3) 
 YvCkap2Tk1 = Matmul(Yv,Ckap2Tk1) 
 YvCkap2Tk2 = Matmul(Yv,Ckap2Tk2) 
 YvCkap2Tk3 = Matmul(Yv,Ckap2Tk3) 
 YvCkap3lam = Matmul(Yv,Ckap3lam) 
 YvCkap3kap1 = Matmul(Yv,Ckap3kap1) 
 YvCkap3kap2 = Matmul(Yv,Ckap3kap2) 
 YvCkap3Tk1 = Matmul(Yv,Ckap3Tk1) 
 YvCkap3Tk2 = Matmul(Yv,Ckap3Tk2) 
 YvCkap3Tk3 = Matmul(Yv,Ckap3Tk3) 
 adjYdYdadjYd = Matmul(adjYd,YdadjYd) 
 adjYdYdadjYu = Matmul(adjYd,YdadjYu) 
 adjYdYdadjTd = Matmul(adjYd,YdadjTd) 
 adjYdYdadjTu = Matmul(adjYd,YdadjTu) 
 adjYdTdadjYd = Matmul(adjYd,TdadjYd) 
 adjYdTdadjYu = Matmul(adjYd,TdadjYu) 
 adjYdTdadjTd = Matmul(adjYd,TdadjTd) 
 adjYdTdadjTu = Matmul(adjYd,TdadjTu) 
 adjYeYevL = Matmul(adjYe,YevL) 
 adjYeYeadjYe = Matmul(adjYe,YeadjYe) 
 adjYeYeadjTe = Matmul(adjYe,YeadjTe) 
 adjYeYeCYv = Matmul(adjYe,YeCYv) 
 adjYeYeCTv = Matmul(adjYe,YeCTv) 
 adjYeTeadjYe = Matmul(adjYe,TeadjYe) 
 adjYeTeadjTe = Matmul(adjYe,TeadjTe) 
 adjYeTeCYv = Matmul(adjYe,TeCYv) 
 adjYeTeCTv = Matmul(adjYe,TeCTv) 
 adjYuYuadjYd = Matmul(adjYu,YuadjYd) 
 adjYuYuadjYu = Matmul(adjYu,YuadjYu) 
 adjYuYuadjTd = Matmul(adjYu,YuadjTd) 
 adjYuYuadjTu = Matmul(adjYu,YuadjTu) 
 adjYuTuadjYd = Matmul(adjYu,TuadjYd) 
 adjYuTuadjYu = Matmul(adjYu,TuadjYu) 
 adjYuTuadjTd = Matmul(adjYu,TuadjTd) 
 adjYuTuadjTu = Matmul(adjYu,TuadjTu) 
 adjYvYvadjYv = Matmul(adjYv,YvadjYv) 
 adjYvYvadjkap1 = Matmul(adjYv,Yvadjkap1) 
 adjYvYvadjkap2 = Matmul(adjYv,Yvadjkap2) 
 adjYvYvadjkap3 = Matmul(adjYv,Yvadjkap3) 
 adjYvYvadjTk1 = Matmul(adjYv,YvadjTk1) 
 adjYvYvadjTk2 = Matmul(adjYv,YvadjTk2) 
 adjYvYvadjTk3 = Matmul(adjYv,YvadjTk3) 
 adjYvYvClam = Matmul(adjYv,YvClam) 
 adjYvYvCTlam = Matmul(adjYv,YvCTlam) 
 adjYvYvCkap1 = Matmul(adjYv,YvCkap1) 
 adjYvYvCkap2 = Matmul(adjYv,YvCkap2) 
 adjYvYvCkap3 = Matmul(adjYv,YvCkap3) 
 adjYvTvadjYv = Matmul(adjYv,TvadjYv) 
 adjYvTvadjTv = Matmul(adjYv,TvadjTv) 
 adjYvTvClam = Matmul(adjYv,TvClam) 
 adjYvTvCTlam = Matmul(adjYv,TvCTlam) 
 adjYvTpYeCYe = Matmul(adjYv,TpYeCYe) 
 adjYvTpTeCTe = Matmul(adjYv,TpTeCTe) 
 adjTdYdadjYd = Matmul(adjTd,YdadjYd) 
 adjTdYdadjYu = Matmul(adjTd,YdadjYu) 
 adjTdTdadjYd = Matmul(adjTd,TdadjYd) 
 adjTdTdadjYu = Matmul(adjTd,TdadjYu) 
 adjTeYeadjYe = Matmul(adjTe,YeadjYe) 
 adjTeYeCYv = Matmul(adjTe,YeCYv) 
 adjTeTeadjYe = Matmul(adjTe,TeadjYe) 
 adjTeTeCYv = Matmul(adjTe,TeCYv) 
 adjTuYuadjYd = Matmul(adjTu,YuadjYd) 
 adjTuYuadjYu = Matmul(adjTu,YuadjYu) 
 adjTuTuadjYd = Matmul(adjTu,TuadjYd) 
 adjTuTuadjYu = Matmul(adjTu,TuadjYu) 
 adjTvYvClam = Matmul(adjTv,YvClam) 
 adjTvTvadjYv = Matmul(adjTv,TvadjYv) 
 adjTvTvClam = Matmul(adjTv,TvClam) 
 adjkap1mv2lam = Matmul(adjkap1,mv2lam) 
 adjkap1mv2kap1 = Matmul(adjkap1,mv2kap1) 
 adjkap1mv2kap2 = Matmul(adjkap1,mv2kap2) 
 adjkap1mv2kap3 = Matmul(adjkap1,mv2kap3) 
 adjkap1TpYvml2 = Matmul(adjkap1,TpYvml2) 
 adjkap1kap1Clam = Matmul(adjkap1,kap1Clam) 
 adjkap1kap2Clam = Matmul(adjkap1,kap2Clam) 
 adjkap1kap3Clam = Matmul(adjkap1,kap3Clam) 
 adjkap2mv2lam = Matmul(adjkap2,mv2lam) 
 adjkap2mv2kap1 = Matmul(adjkap2,mv2kap1) 
 adjkap2mv2kap2 = Matmul(adjkap2,mv2kap2) 
 adjkap2mv2kap3 = Matmul(adjkap2,mv2kap3) 
 adjkap2TpYvml2 = Matmul(adjkap2,TpYvml2) 
 adjkap2kap1Clam = Matmul(adjkap2,kap1Clam) 
 adjkap2kap2Clam = Matmul(adjkap2,kap2Clam) 
 adjkap2kap3Clam = Matmul(adjkap2,kap3Clam) 
 adjkap3mv2lam = Matmul(adjkap3,mv2lam) 
 adjkap3mv2kap1 = Matmul(adjkap3,mv2kap1) 
 adjkap3mv2kap2 = Matmul(adjkap3,mv2kap2) 
 adjkap3mv2kap3 = Matmul(adjkap3,mv2kap3) 
 adjkap3TpYvml2 = Matmul(adjkap3,TpYvml2) 
 adjkap3kap1Clam = Matmul(adjkap3,kap1Clam) 
 adjkap3kap2Clam = Matmul(adjkap3,kap2Clam) 
 adjkap3kap3Clam = Matmul(adjkap3,kap3Clam) 
 Cml2YvadjYv = Matmul(Conjg(ml2),YvadjYv) 
 Cml2Yvadjkap1 = Matmul(Conjg(ml2),Yvadjkap1) 
 Cml2Yvadjkap2 = Matmul(Conjg(ml2),Yvadjkap2) 
 Cml2Yvadjkap3 = Matmul(Conjg(ml2),Yvadjkap3) 
 Cml2YvClam = Matmul(Conjg(ml2),YvClam) 
 Cml2TpYeCYe = Matmul(Conjg(ml2),TpYeCYe) 
 Cmv2Ckap1kap1 = Matmul(Conjg(mv2),Ckap1kap1) 
 Cmv2Ckap1kap2 = Matmul(Conjg(mv2),Ckap1kap2) 
 Cmv2Ckap1kap3 = Matmul(Conjg(mv2),Ckap1kap3) 
 Cmv2Ckap2kap1 = Matmul(Conjg(mv2),Ckap2kap1) 
 Cmv2Ckap2kap2 = Matmul(Conjg(mv2),Ckap2kap2) 
 Cmv2Ckap2kap3 = Matmul(Conjg(mv2),Ckap2kap3) 
 Cmv2Ckap3kap1 = Matmul(Conjg(mv2),Ckap3kap1) 
 Cmv2Ckap3kap2 = Matmul(Conjg(mv2),Ckap3kap2) 
 Cmv2Ckap3kap3 = Matmul(Conjg(mv2),Ckap3kap3) 
 CYeYvClam = Matmul(Conjg(Ye),YvClam) 
 CYeTvClam = Matmul(Conjg(Ye),TvClam) 
 CYvTpYvvL = Matmul(Conjg(Yv),TpYvvL) 
 CYvTpYvadjYe = Matmul(Conjg(Yv),TpYvadjYe) 
 CYvTpYvadjTe = Matmul(Conjg(Yv),TpYvadjTe) 
 CYvTpYvCYv = Matmul(Conjg(Yv),TpYvCYv) 
 CYvTpYvCTv = Matmul(Conjg(Yv),TpYvCTv) 
 CYvTpTvadjTe = Matmul(Conjg(Yv),TpTvadjTe) 
 CYvTpTvCTv = Matmul(Conjg(Yv),TpTvCTv) 
 CYvTpkap1Ckap1 = Matmul(Conjg(Yv),Tpkap1Ckap1) 
 CYvTpkap1Ckap2 = Matmul(Conjg(Yv),Tpkap1Ckap2) 
 CYvTpkap1Ckap3 = Matmul(Conjg(Yv),Tpkap1Ckap3) 
 CYvTpkap2Ckap1 = Matmul(Conjg(Yv),Tpkap2Ckap1) 
 CYvTpkap2Ckap2 = Matmul(Conjg(Yv),Tpkap2Ckap2) 
 CYvTpkap2Ckap3 = Matmul(Conjg(Yv),Tpkap2Ckap3) 
 CYvTpkap3Ckap1 = Matmul(Conjg(Yv),Tpkap3Ckap1) 
 CYvTpkap3Ckap2 = Matmul(Conjg(Yv),Tpkap3Ckap2) 
 CYvTpkap3Ckap3 = Matmul(Conjg(Yv),Tpkap3Ckap3) 
 CTvTpYvadjYe = Matmul(Conjg(Tv),TpYvadjYe) 
 CTvTpYvCYv = Matmul(Conjg(Tv),TpYvCYv) 
 CTvTpTvadjYe = Matmul(Conjg(Tv),TpTvadjYe) 
 CTvTpTvCYv = Matmul(Conjg(Tv),TpTvCYv) 
 Ckap1mv2lam = Matmul(Conjg(kap1),mv2lam) 
 Ckap1TpYvCYv = Matmul(Conjg(kap1),TpYvCYv) 
 Ckap1TpYvCTv = Matmul(Conjg(kap1),TpYvCTv) 
 Ckap1Tpkap1adjYv = Matmul(Conjg(kap1),Tpkap1adjYv) 
 Ckap1Tpkap1Clam = Matmul(Conjg(kap1),Tpkap1Clam) 
 Ckap1TpTk1adjYv = Matmul(Conjg(kap1),TpTk1adjYv) 
 Ckap1kap1adjYv = Matmul(Conjg(kap1),kap1adjYv) 
 Ckap1kap1Ckap1 = Matmul(Conjg(kap1),kap1Ckap1) 
 Ckap1kap1Ckap2 = Matmul(Conjg(kap1),kap1Ckap2) 
 Ckap1kap1Ckap3 = Matmul(Conjg(kap1),kap1Ckap3) 
 Ckap1kap2Ckap1 = Matmul(Conjg(kap1),kap2Ckap1) 
 Ckap1kap2Ckap2 = Matmul(Conjg(kap1),kap2Ckap2) 
 Ckap1kap2Ckap3 = Matmul(Conjg(kap1),kap2Ckap3) 
 Ckap1kap3Ckap1 = Matmul(Conjg(kap1),kap3Ckap1) 
 Ckap1kap3Ckap2 = Matmul(Conjg(kap1),kap3Ckap2) 
 Ckap1kap3Ckap3 = Matmul(Conjg(kap1),kap3Ckap3) 
 Ckap1Tk1adjTv = Matmul(Conjg(kap1),Tk1adjTv) 
 Ckap2mv2lam = Matmul(Conjg(kap2),mv2lam) 
 Ckap2TpYvCYv = Matmul(Conjg(kap2),TpYvCYv) 
 Ckap2TpYvCTv = Matmul(Conjg(kap2),TpYvCTv) 
 Ckap2Tpkap2adjYv = Matmul(Conjg(kap2),Tpkap2adjYv) 
 Ckap2Tpkap2Clam = Matmul(Conjg(kap2),Tpkap2Clam) 
 Ckap2TpTk2adjYv = Matmul(Conjg(kap2),TpTk2adjYv) 
 Ckap2kap1Ckap1 = Matmul(Conjg(kap2),kap1Ckap1) 
 Ckap2kap1Ckap2 = Matmul(Conjg(kap2),kap1Ckap2) 
 Ckap2kap1Ckap3 = Matmul(Conjg(kap2),kap1Ckap3) 
 Ckap2kap2adjYv = Matmul(Conjg(kap2),kap2adjYv) 
 Ckap2kap2Ckap1 = Matmul(Conjg(kap2),kap2Ckap1) 
 Ckap2kap2Ckap2 = Matmul(Conjg(kap2),kap2Ckap2) 
 Ckap2kap2Ckap3 = Matmul(Conjg(kap2),kap2Ckap3) 
 Ckap2kap3Ckap1 = Matmul(Conjg(kap2),kap3Ckap1) 
 Ckap2kap3Ckap2 = Matmul(Conjg(kap2),kap3Ckap2) 
 Ckap2kap3Ckap3 = Matmul(Conjg(kap2),kap3Ckap3) 
 Ckap2Tk2adjTv = Matmul(Conjg(kap2),Tk2adjTv) 
 Ckap3mv2lam = Matmul(Conjg(kap3),mv2lam) 
 Ckap3TpYvCYv = Matmul(Conjg(kap3),TpYvCYv) 
 Ckap3TpYvCTv = Matmul(Conjg(kap3),TpYvCTv) 
 Ckap3Tpkap3adjYv = Matmul(Conjg(kap3),Tpkap3adjYv) 
 Ckap3Tpkap3Clam = Matmul(Conjg(kap3),Tpkap3Clam) 
 Ckap3TpTk3adjYv = Matmul(Conjg(kap3),TpTk3adjYv) 
 Ckap3kap1Ckap1 = Matmul(Conjg(kap3),kap1Ckap1) 
 Ckap3kap1Ckap2 = Matmul(Conjg(kap3),kap1Ckap2) 
 Ckap3kap1Ckap3 = Matmul(Conjg(kap3),kap1Ckap3) 
 Ckap3kap2Ckap1 = Matmul(Conjg(kap3),kap2Ckap1) 
 Ckap3kap2Ckap2 = Matmul(Conjg(kap3),kap2Ckap2) 
 Ckap3kap2Ckap3 = Matmul(Conjg(kap3),kap2Ckap3) 
 Ckap3kap3adjYv = Matmul(Conjg(kap3),kap3adjYv) 
 Ckap3kap3Ckap1 = Matmul(Conjg(kap3),kap3Ckap1) 
 Ckap3kap3Ckap2 = Matmul(Conjg(kap3),kap3Ckap2) 
 Ckap3kap3Ckap3 = Matmul(Conjg(kap3),kap3Ckap3) 
 Ckap3Tk3adjTv = Matmul(Conjg(kap3),Tk3adjTv) 
 CTk1TpTk1adjYv = Matmul(Conjg(Tk1),TpTk1adjYv) 
 CTk2TpTk2adjYv = Matmul(Conjg(Tk2),TpTk2adjYv) 
 CTk3TpTk3adjYv = Matmul(Conjg(Tk3),TpTk3adjYv) 
 TdadjTdYd = Matmul(Td,adjTdYd) 
 TeadjTeYe = Matmul(Te,adjTeYe) 
 TeCYvlam = Matmul(Te,CYvlam) 
 TeCTvlam = Matmul(Te,CTvlam) 
 TeCTvTpYv = Matmul(Te,CTvTpYv) 
 TuadjTuYu = Matmul(Tu,adjTuYu) 
 Tvadjkap1kap1 = Matmul(Tv,adjkap1kap1) 
 Tvadjkap1kap2 = Matmul(Tv,adjkap1kap2) 
 Tvadjkap1kap3 = Matmul(Tv,adjkap1kap3) 
 Tvadjkap2kap1 = Matmul(Tv,adjkap2kap1) 
 Tvadjkap2kap2 = Matmul(Tv,adjkap2kap2) 
 Tvadjkap2kap3 = Matmul(Tv,adjkap2kap3) 
 Tvadjkap3kap1 = Matmul(Tv,adjkap3kap1) 
 Tvadjkap3kap2 = Matmul(Tv,adjkap3kap2) 
 Tvadjkap3kap3 = Matmul(Tv,adjkap3kap3) 
 TvCkap1kap1 = Matmul(Tv,Ckap1kap1) 
 TvCkap1kap2 = Matmul(Tv,Ckap1kap2) 
 TvCkap1kap3 = Matmul(Tv,Ckap1kap3) 
 TvCkap2kap1 = Matmul(Tv,Ckap2kap1) 
 TvCkap2kap2 = Matmul(Tv,Ckap2kap2) 
 TvCkap2kap3 = Matmul(Tv,Ckap2kap3) 
 TvCkap3kap1 = Matmul(Tv,Ckap3kap1) 
 TvCkap3kap2 = Matmul(Tv,Ckap3kap2) 
 TvCkap3kap3 = Matmul(Tv,Ckap3kap3) 
 TpYeCme2CYe = Matmul(Transpose(Ye),Cme2CYe) 
Forall(i2=1:3)  TpYeCme2CYe(i2,i2) =  Real(TpYeCme2CYe(i2,i2),dp) 
 TpYeCYeCml2 = Matmul(Transpose(Ye),CYeCml2) 
 TpYeCTeTv = Matmul(Transpose(Ye),CTeTv) 
 TpYvml2adjYe = Matmul(Transpose(Yv),ml2adjYe) 
 TpYvadjYeme2 = Matmul(Transpose(Yv),adjYeme2) 
 TpYvadjYeYe = Matmul(Transpose(Yv),adjYeYe) 
 TpYvadjYeTe = Matmul(Transpose(Yv),adjYeTe) 
 TpYvCYvTpYv = Matmul(Transpose(Yv),CYvTpYv) 
 TpYvCYvTpTv = Matmul(Transpose(Yv),CYvTpTv) 
 TpYvCYvTpkap1 = Matmul(Transpose(Yv),CYvTpkap1) 
 TpYvCYvTpkap2 = Matmul(Transpose(Yv),CYvTpkap2) 
 TpYvCYvTpkap3 = Matmul(Transpose(Yv),CYvTpkap3) 
 TpYvCYvTpTk1 = Matmul(Transpose(Yv),CYvTpTk1) 
 TpYvCYvTpTk2 = Matmul(Transpose(Yv),CYvTpTk2) 
 TpYvCYvTpTk3 = Matmul(Transpose(Yv),CYvTpTk3) 
 TpYvCTvTlam = Matmul(Transpose(Yv),CTvTlam) 
 TpYvCTvTpTv = Matmul(Transpose(Yv),CTvTpTv) 
 TpTvadjYeYe = Matmul(Transpose(Tv),adjYeYe) 
 TpTvCYvTpYv = Matmul(Transpose(Tv),CYvTpYv) 
 TpTvCTvlam = Matmul(Transpose(Tv),CTvlam) 
 TpTvCTvTpYv = Matmul(Transpose(Tv),CTvTpYv) 
 Tpkap1adjYvvL = Matmul(Transpose(kap1),adjYvvL) 
 Tpkap1adjYvYv = Matmul(Transpose(kap1),adjYvYv) 
 Tpkap1adjkap1lam = Matmul(Transpose(kap1),adjkap1lam) 
 Tpkap1adjkap1Tlam = Matmul(Transpose(kap1),adjkap1Tlam) 
 Tpkap1adjkap2lam = Matmul(Transpose(kap1),adjkap2lam) 
 Tpkap1adjkap2Tlam = Matmul(Transpose(kap1),adjkap2Tlam) 
 Tpkap1adjkap3lam = Matmul(Transpose(kap1),adjkap3lam) 
 Tpkap1adjkap3Tlam = Matmul(Transpose(kap1),adjkap3Tlam) 
 Tpkap1Cmv2adjYv = Matmul(Transpose(kap1),Cmv2adjYv) 
 Tpkap1Cmv2adjkap1 = Matmul(Transpose(kap1),Cmv2adjkap1) 
 Tpkap1Cmv2adjkap2 = Matmul(Transpose(kap1),Cmv2adjkap2) 
 Tpkap1Cmv2adjkap3 = Matmul(Transpose(kap1),Cmv2adjkap3) 
 Tpkap1Cmv2Clam = Matmul(Transpose(kap1),Cmv2Clam) 
 Tpkap1Ckap1lam = Matmul(Transpose(kap1),Ckap1lam) 
 Tpkap1Ckap1Tlam = Matmul(Transpose(kap1),Ckap1Tlam) 
 Tpkap1Ckap1TpTv = Matmul(Transpose(kap1),Ckap1TpTv) 
 Tpkap1Ckap1TpTk1 = Matmul(Transpose(kap1),Ckap1TpTk1) 
 Tpkap1Ckap1TpTk2 = Matmul(Transpose(kap1),Ckap1TpTk2) 
 Tpkap1Ckap1TpTk3 = Matmul(Transpose(kap1),Ckap1TpTk3) 
 Tpkap1Ckap2lam = Matmul(Transpose(kap1),Ckap2lam) 
 Tpkap1Ckap2Tlam = Matmul(Transpose(kap1),Ckap2Tlam) 
 Tpkap1Ckap2TpTv = Matmul(Transpose(kap1),Ckap2TpTv) 
 Tpkap1Ckap2TpTk1 = Matmul(Transpose(kap1),Ckap2TpTk1) 
 Tpkap1Ckap2TpTk2 = Matmul(Transpose(kap1),Ckap2TpTk2) 
 Tpkap1Ckap2TpTk3 = Matmul(Transpose(kap1),Ckap2TpTk3) 
 Tpkap1Ckap3lam = Matmul(Transpose(kap1),Ckap3lam) 
 Tpkap1Ckap3Tlam = Matmul(Transpose(kap1),Ckap3Tlam) 
 Tpkap1Ckap3TpTv = Matmul(Transpose(kap1),Ckap3TpTv) 
 Tpkap1Ckap3TpTk1 = Matmul(Transpose(kap1),Ckap3TpTk1) 
 Tpkap1Ckap3TpTk2 = Matmul(Transpose(kap1),Ckap3TpTk2) 
 Tpkap1Ckap3TpTk3 = Matmul(Transpose(kap1),Ckap3TpTk3) 
 Tpkap2adjYvvL = Matmul(Transpose(kap2),adjYvvL) 
 Tpkap2adjYvYv = Matmul(Transpose(kap2),adjYvYv) 
 Tpkap2adjkap1lam = Matmul(Transpose(kap2),adjkap1lam) 
 Tpkap2adjkap1Tlam = Matmul(Transpose(kap2),adjkap1Tlam) 
 Tpkap2adjkap2lam = Matmul(Transpose(kap2),adjkap2lam) 
 Tpkap2adjkap2Tlam = Matmul(Transpose(kap2),adjkap2Tlam) 
 Tpkap2adjkap3lam = Matmul(Transpose(kap2),adjkap3lam) 
 Tpkap2adjkap3Tlam = Matmul(Transpose(kap2),adjkap3Tlam) 
 Tpkap2Cmv2adjYv = Matmul(Transpose(kap2),Cmv2adjYv) 
 Tpkap2Cmv2adjkap1 = Matmul(Transpose(kap2),Cmv2adjkap1) 
 Tpkap2Cmv2adjkap2 = Matmul(Transpose(kap2),Cmv2adjkap2) 
 Tpkap2Cmv2adjkap3 = Matmul(Transpose(kap2),Cmv2adjkap3) 
 Tpkap2Cmv2Clam = Matmul(Transpose(kap2),Cmv2Clam) 
 Tpkap2Ckap1lam = Matmul(Transpose(kap2),Ckap1lam) 
 Tpkap2Ckap1Tlam = Matmul(Transpose(kap2),Ckap1Tlam) 
 Tpkap2Ckap1TpTv = Matmul(Transpose(kap2),Ckap1TpTv) 
 Tpkap2Ckap1TpTk1 = Matmul(Transpose(kap2),Ckap1TpTk1) 
 Tpkap2Ckap1TpTk2 = Matmul(Transpose(kap2),Ckap1TpTk2) 
 Tpkap2Ckap1TpTk3 = Matmul(Transpose(kap2),Ckap1TpTk3) 
 Tpkap2Ckap2lam = Matmul(Transpose(kap2),Ckap2lam) 
 Tpkap2Ckap2Tlam = Matmul(Transpose(kap2),Ckap2Tlam) 
 Tpkap2Ckap2TpTv = Matmul(Transpose(kap2),Ckap2TpTv) 
 Tpkap2Ckap2TpTk1 = Matmul(Transpose(kap2),Ckap2TpTk1) 
 Tpkap2Ckap2TpTk2 = Matmul(Transpose(kap2),Ckap2TpTk2) 
 Tpkap2Ckap2TpTk3 = Matmul(Transpose(kap2),Ckap2TpTk3) 
 Tpkap2Ckap3lam = Matmul(Transpose(kap2),Ckap3lam) 
 Tpkap2Ckap3Tlam = Matmul(Transpose(kap2),Ckap3Tlam) 
 Tpkap2Ckap3TpTv = Matmul(Transpose(kap2),Ckap3TpTv) 
 Tpkap2Ckap3TpTk1 = Matmul(Transpose(kap2),Ckap3TpTk1) 
 Tpkap2Ckap3TpTk2 = Matmul(Transpose(kap2),Ckap3TpTk2) 
 Tpkap2Ckap3TpTk3 = Matmul(Transpose(kap2),Ckap3TpTk3) 
 Tpkap3adjYvvL = Matmul(Transpose(kap3),adjYvvL) 
 Tpkap3adjYvYv = Matmul(Transpose(kap3),adjYvYv) 
 Tpkap3adjkap1lam = Matmul(Transpose(kap3),adjkap1lam) 
 Tpkap3adjkap1Tlam = Matmul(Transpose(kap3),adjkap1Tlam) 
 Tpkap3adjkap2lam = Matmul(Transpose(kap3),adjkap2lam) 
 Tpkap3adjkap2Tlam = Matmul(Transpose(kap3),adjkap2Tlam) 
 Tpkap3adjkap3lam = Matmul(Transpose(kap3),adjkap3lam) 
 Tpkap3adjkap3Tlam = Matmul(Transpose(kap3),adjkap3Tlam) 
 Tpkap3Cmv2adjYv = Matmul(Transpose(kap3),Cmv2adjYv) 
 Tpkap3Cmv2adjkap1 = Matmul(Transpose(kap3),Cmv2adjkap1) 
 Tpkap3Cmv2adjkap2 = Matmul(Transpose(kap3),Cmv2adjkap2) 
 Tpkap3Cmv2adjkap3 = Matmul(Transpose(kap3),Cmv2adjkap3) 
 Tpkap3Cmv2Clam = Matmul(Transpose(kap3),Cmv2Clam) 
 Tpkap3Ckap1lam = Matmul(Transpose(kap3),Ckap1lam) 
 Tpkap3Ckap1Tlam = Matmul(Transpose(kap3),Ckap1Tlam) 
 Tpkap3Ckap1TpTv = Matmul(Transpose(kap3),Ckap1TpTv) 
 Tpkap3Ckap1TpTk1 = Matmul(Transpose(kap3),Ckap1TpTk1) 
 Tpkap3Ckap1TpTk2 = Matmul(Transpose(kap3),Ckap1TpTk2) 
 Tpkap3Ckap1TpTk3 = Matmul(Transpose(kap3),Ckap1TpTk3) 
 Tpkap3Ckap2lam = Matmul(Transpose(kap3),Ckap2lam) 
 Tpkap3Ckap2Tlam = Matmul(Transpose(kap3),Ckap2Tlam) 
 Tpkap3Ckap2TpTv = Matmul(Transpose(kap3),Ckap2TpTv) 
 Tpkap3Ckap2TpTk1 = Matmul(Transpose(kap3),Ckap2TpTk1) 
 Tpkap3Ckap2TpTk2 = Matmul(Transpose(kap3),Ckap2TpTk2) 
 Tpkap3Ckap2TpTk3 = Matmul(Transpose(kap3),Ckap2TpTk3) 
 Tpkap3Ckap3lam = Matmul(Transpose(kap3),Ckap3lam) 
 Tpkap3Ckap3Tlam = Matmul(Transpose(kap3),Ckap3Tlam) 
 Tpkap3Ckap3TpTv = Matmul(Transpose(kap3),Ckap3TpTv) 
 Tpkap3Ckap3TpTk1 = Matmul(Transpose(kap3),Ckap3TpTk1) 
 Tpkap3Ckap3TpTk2 = Matmul(Transpose(kap3),Ckap3TpTk2) 
 Tpkap3Ckap3TpTk3 = Matmul(Transpose(kap3),Ckap3TpTk3) 
 TpTk1adjYvYv = Matmul(Transpose(Tk1),adjYvYv) 
 TpTk1adjkap1lam = Matmul(Transpose(Tk1),adjkap1lam) 
 TpTk1adjkap2lam = Matmul(Transpose(Tk1),adjkap2lam) 
 TpTk1adjkap3lam = Matmul(Transpose(Tk1),adjkap3lam) 
 TpTk1Ckap1lam = Matmul(Transpose(Tk1),Ckap1lam) 
 TpTk1Ckap2lam = Matmul(Transpose(Tk1),Ckap2lam) 
 TpTk1Ckap3lam = Matmul(Transpose(Tk1),Ckap3lam) 
 TpTk2adjYvYv = Matmul(Transpose(Tk2),adjYvYv) 
 TpTk2adjkap1lam = Matmul(Transpose(Tk2),adjkap1lam) 
 TpTk2adjkap2lam = Matmul(Transpose(Tk2),adjkap2lam) 
 TpTk2adjkap3lam = Matmul(Transpose(Tk2),adjkap3lam) 
 TpTk2Ckap1lam = Matmul(Transpose(Tk2),Ckap1lam) 
 TpTk2Ckap2lam = Matmul(Transpose(Tk2),Ckap2lam) 
 TpTk2Ckap3lam = Matmul(Transpose(Tk2),Ckap3lam) 
 TpTk3adjYvYv = Matmul(Transpose(Tk3),adjYvYv) 
 TpTk3adjkap1lam = Matmul(Transpose(Tk3),adjkap1lam) 
 TpTk3adjkap2lam = Matmul(Transpose(Tk3),adjkap2lam) 
 TpTk3adjkap3lam = Matmul(Transpose(Tk3),adjkap3lam) 
 TpTk3Ckap1lam = Matmul(Transpose(Tk3),Ckap1lam) 
 TpTk3Ckap2lam = Matmul(Transpose(Tk3),Ckap2lam) 
 TpTk3Ckap3lam = Matmul(Transpose(Tk3),Ckap3lam) 
 kap1adjYvCml2 = Matmul(kap1,adjYvCml2) 
 kap1adjYvTv = Matmul(kap1,adjYvTv) 
 kap1adjkap1lam = Matmul(kap1,adjkap1lam) 
 kap1adjkap1TpYv = Matmul(kap1,adjkap1TpYv) 
 kap1adjkap1Tpkap1 = Matmul(kap1,adjkap1Tpkap1) 
 kap1adjkap1Tpkap2 = Matmul(kap1,adjkap1Tpkap2) 
 kap1adjkap1Tpkap3 = Matmul(kap1,adjkap1Tpkap3) 
 kap1Cmv2Clam = Matmul(kap1,Cmv2Clam) 
 kap1Ckap1lam = Matmul(kap1,Ckap1lam) 
 kap1Ckap1kap1 = Matmul(kap1,Ckap1kap1) 
 kap1Ckap1kap2 = Matmul(kap1,Ckap1kap2) 
 kap1Ckap1kap3 = Matmul(kap1,Ckap1kap3) 
 kap1Ckap1Tk1 = Matmul(kap1,Ckap1Tk1) 
 kap1Ckap1Tk2 = Matmul(kap1,Ckap1Tk2) 
 kap1Ckap1Tk3 = Matmul(kap1,Ckap1Tk3) 
 kap1Ckap2lam = Matmul(kap1,Ckap2lam) 
 kap1Ckap2kap1 = Matmul(kap1,Ckap2kap1) 
 kap1Ckap2kap2 = Matmul(kap1,Ckap2kap2) 
 kap1Ckap2kap3 = Matmul(kap1,Ckap2kap3) 
 kap1Ckap2Tk1 = Matmul(kap1,Ckap2Tk1) 
 kap1Ckap2Tk2 = Matmul(kap1,Ckap2Tk2) 
 kap1Ckap2Tk3 = Matmul(kap1,Ckap2Tk3) 
 kap1Ckap3lam = Matmul(kap1,Ckap3lam) 
 kap1Ckap3kap1 = Matmul(kap1,Ckap3kap1) 
 kap1Ckap3kap2 = Matmul(kap1,Ckap3kap2) 
 kap1Ckap3kap3 = Matmul(kap1,Ckap3kap3) 
 kap1Ckap3Tk1 = Matmul(kap1,Ckap3Tk1) 
 kap1Ckap3Tk2 = Matmul(kap1,Ckap3Tk2) 
 kap1Ckap3Tk3 = Matmul(kap1,Ckap3Tk3) 
 kap2adjYvCml2 = Matmul(kap2,adjYvCml2) 
 kap2adjYvTv = Matmul(kap2,adjYvTv) 
 kap2adjkap2lam = Matmul(kap2,adjkap2lam) 
 kap2adjkap2TpYv = Matmul(kap2,adjkap2TpYv) 
 kap2adjkap2Tpkap1 = Matmul(kap2,adjkap2Tpkap1) 
 kap2adjkap2Tpkap2 = Matmul(kap2,adjkap2Tpkap2) 
 kap2adjkap2Tpkap3 = Matmul(kap2,adjkap2Tpkap3) 
 kap2Cmv2Clam = Matmul(kap2,Cmv2Clam) 
 kap2Ckap1lam = Matmul(kap2,Ckap1lam) 
 kap2Ckap1kap1 = Matmul(kap2,Ckap1kap1) 
 kap2Ckap1kap2 = Matmul(kap2,Ckap1kap2) 
 kap2Ckap1kap3 = Matmul(kap2,Ckap1kap3) 
 kap2Ckap1Tk1 = Matmul(kap2,Ckap1Tk1) 
 kap2Ckap1Tk2 = Matmul(kap2,Ckap1Tk2) 
 kap2Ckap1Tk3 = Matmul(kap2,Ckap1Tk3) 
 kap2Ckap2lam = Matmul(kap2,Ckap2lam) 
 kap2Ckap2kap1 = Matmul(kap2,Ckap2kap1) 
 kap2Ckap2kap2 = Matmul(kap2,Ckap2kap2) 
 kap2Ckap2kap3 = Matmul(kap2,Ckap2kap3) 
 kap2Ckap2Tk1 = Matmul(kap2,Ckap2Tk1) 
 kap2Ckap2Tk2 = Matmul(kap2,Ckap2Tk2) 
 kap2Ckap2Tk3 = Matmul(kap2,Ckap2Tk3) 
 kap2Ckap3lam = Matmul(kap2,Ckap3lam) 
 kap2Ckap3kap1 = Matmul(kap2,Ckap3kap1) 
 kap2Ckap3kap2 = Matmul(kap2,Ckap3kap2) 
 kap2Ckap3kap3 = Matmul(kap2,Ckap3kap3) 
 kap2Ckap3Tk1 = Matmul(kap2,Ckap3Tk1) 
 kap2Ckap3Tk2 = Matmul(kap2,Ckap3Tk2) 
 kap2Ckap3Tk3 = Matmul(kap2,Ckap3Tk3) 
 kap3adjYvCml2 = Matmul(kap3,adjYvCml2) 
 kap3adjYvTv = Matmul(kap3,adjYvTv) 
 kap3adjkap3lam = Matmul(kap3,adjkap3lam) 
 kap3adjkap3TpYv = Matmul(kap3,adjkap3TpYv) 
 kap3adjkap3Tpkap1 = Matmul(kap3,adjkap3Tpkap1) 
 kap3adjkap3Tpkap2 = Matmul(kap3,adjkap3Tpkap2) 
 kap3adjkap3Tpkap3 = Matmul(kap3,adjkap3Tpkap3) 
 kap3Cmv2Clam = Matmul(kap3,Cmv2Clam) 
 kap3Ckap1lam = Matmul(kap3,Ckap1lam) 
 kap3Ckap1kap1 = Matmul(kap3,Ckap1kap1) 
 kap3Ckap1kap2 = Matmul(kap3,Ckap1kap2) 
 kap3Ckap1kap3 = Matmul(kap3,Ckap1kap3) 
 kap3Ckap1Tk1 = Matmul(kap3,Ckap1Tk1) 
 kap3Ckap1Tk2 = Matmul(kap3,Ckap1Tk2) 
 kap3Ckap1Tk3 = Matmul(kap3,Ckap1Tk3) 
 kap3Ckap2lam = Matmul(kap3,Ckap2lam) 
 kap3Ckap2kap1 = Matmul(kap3,Ckap2kap1) 
 kap3Ckap2kap2 = Matmul(kap3,Ckap2kap2) 
 kap3Ckap2kap3 = Matmul(kap3,Ckap2kap3) 
 kap3Ckap2Tk1 = Matmul(kap3,Ckap2Tk1) 
 kap3Ckap2Tk2 = Matmul(kap3,Ckap2Tk2) 
 kap3Ckap2Tk3 = Matmul(kap3,Ckap2Tk3) 
 kap3Ckap3lam = Matmul(kap3,Ckap3lam) 
 kap3Ckap3kap1 = Matmul(kap3,Ckap3kap1) 
 kap3Ckap3kap2 = Matmul(kap3,Ckap3kap2) 
 kap3Ckap3kap3 = Matmul(kap3,Ckap3kap3) 
 kap3Ckap3Tk1 = Matmul(kap3,Ckap3Tk1) 
 kap3Ckap3Tk2 = Matmul(kap3,Ckap3Tk2) 
 kap3Ckap3Tk3 = Matmul(kap3,Ckap3Tk3) 
 Tk1adjkap1TpYv = Matmul(Tk1,adjkap1TpYv) 
 Tk1adjTk1lam = Matmul(Tk1,adjTk1lam) 
 Tk1adjTk1TpYv = Matmul(Tk1,adjTk1TpYv) 
 Tk1Ckap1kap1 = Matmul(Tk1,Ckap1kap1) 
 Tk1Ckap1kap2 = Matmul(Tk1,Ckap1kap2) 
 Tk1Ckap1kap3 = Matmul(Tk1,Ckap1kap3) 
 Tk1Ckap2kap1 = Matmul(Tk1,Ckap2kap1) 
 Tk1Ckap2kap2 = Matmul(Tk1,Ckap2kap2) 
 Tk1Ckap2kap3 = Matmul(Tk1,Ckap2kap3) 
 Tk1Ckap3kap1 = Matmul(Tk1,Ckap3kap1) 
 Tk1Ckap3kap2 = Matmul(Tk1,Ckap3kap2) 
 Tk1Ckap3kap3 = Matmul(Tk1,Ckap3kap3) 
 Tk2adjkap2TpYv = Matmul(Tk2,adjkap2TpYv) 
 Tk2adjTk2lam = Matmul(Tk2,adjTk2lam) 
 Tk2adjTk2TpYv = Matmul(Tk2,adjTk2TpYv) 
 Tk2Ckap1kap1 = Matmul(Tk2,Ckap1kap1) 
 Tk2Ckap1kap2 = Matmul(Tk2,Ckap1kap2) 
 Tk2Ckap1kap3 = Matmul(Tk2,Ckap1kap3) 
 Tk2Ckap2kap1 = Matmul(Tk2,Ckap2kap1) 
 Tk2Ckap2kap2 = Matmul(Tk2,Ckap2kap2) 
 Tk2Ckap2kap3 = Matmul(Tk2,Ckap2kap3) 
 Tk2Ckap3kap1 = Matmul(Tk2,Ckap3kap1) 
 Tk2Ckap3kap2 = Matmul(Tk2,Ckap3kap2) 
 Tk2Ckap3kap3 = Matmul(Tk2,Ckap3kap3) 
 Tk3adjkap3TpYv = Matmul(Tk3,adjkap3TpYv) 
 Tk3adjTk3lam = Matmul(Tk3,adjTk3lam) 
 Tk3adjTk3TpYv = Matmul(Tk3,adjTk3TpYv) 
 Tk3Ckap1kap1 = Matmul(Tk3,Ckap1kap1) 
 Tk3Ckap1kap2 = Matmul(Tk3,Ckap1kap2) 
 Tk3Ckap1kap3 = Matmul(Tk3,Ckap1kap3) 
 Tk3Ckap2kap1 = Matmul(Tk3,Ckap2kap1) 
 Tk3Ckap2kap2 = Matmul(Tk3,Ckap2kap2) 
 Tk3Ckap2kap3 = Matmul(Tk3,Ckap2kap3) 
 Tk3Ckap3kap1 = Matmul(Tk3,Ckap3kap1) 
 Tk3Ckap3kap2 = Matmul(Tk3,Ckap3kap2) 
 Tk3Ckap3kap3 = Matmul(Tk3,Ckap3kap3) 
 md2YdadjYdYd = Matmul(md2,YdadjYdYd) 
 me2YeadjYeYe = Matmul(me2,YeadjYeYe) 
 me2YeCYvlam = Matmul(me2,YeCYvlam) 
 ml2adjYeYeadjYe = Matmul(ml2,adjYeYeadjYe) 
 ml2adjYeYeCYv = Matmul(ml2,adjYeYeCYv) 
 ml2CYvTpYvadjYe = Matmul(ml2,CYvTpYvadjYe) 
 ml2CYvTpYvCYv = Matmul(ml2,CYvTpYvCYv) 
 mq2adjYdYdadjYd = Matmul(mq2,adjYdYdadjYd) 
 mq2adjYdYdadjYu = Matmul(mq2,adjYdYdadjYu) 
 mq2adjYuYuadjYd = Matmul(mq2,adjYuYuadjYd) 
 mq2adjYuYuadjYu = Matmul(mq2,adjYuYuadjYu) 
 mu2YuadjYuYu = Matmul(mu2,YuadjYuYu) 
 mv2TpYvCYvlam = Matmul(mv2,TpYvCYvlam) 
 mv2TpYvCYvTpYv = Matmul(mv2,TpYvCYvTpYv) 
 mv2kap1adjkap1TpYv = Matmul(mv2,kap1adjkap1TpYv) 
 mv2kap2adjkap2TpYv = Matmul(mv2,kap2adjkap2TpYv) 
 mv2kap3adjkap3TpYv = Matmul(mv2,kap3adjkap3TpYv) 
 Ydmq2adjYdYd = Matmul(Yd,mq2adjYdYd) 
 YdadjYdmd2Yd = Matmul(Yd,adjYdmd2Yd) 
 YdadjYdYdmq2 = Matmul(Yd,adjYdYdmq2) 
 YdadjYdYdadjYd = Matmul(Yd,adjYdYdadjYd) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTdadjYd = Matmul(Yd,adjYdTdadjYd) 
 YdadjYdTdadjTd = Matmul(Yd,adjYdTdadjTd) 
 YdadjYuYuadjYd = Matmul(Yd,adjYuYuadjYd) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YdadjYuTuadjYd = Matmul(Yd,adjYuTuadjYd) 
 YdadjYuTuadjTd = Matmul(Yd,adjYuTuadjTd) 
 YdadjTdTdadjYd = Matmul(Yd,adjTdTdadjYd) 
 YdadjTuTuadjYd = Matmul(Yd,adjTuTuadjYd) 
 Yeml2adjYeYe = Matmul(Ye,ml2adjYeYe) 
 Yeml2CYvlam = Matmul(Ye,ml2CYvlam) 
 YeadjYeme2Ye = Matmul(Ye,adjYeme2Ye) 
 YeadjYeYeml2 = Matmul(Ye,adjYeYeml2) 
 YeadjYeYeadjYe = Matmul(Ye,adjYeYeadjYe) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYeTeadjYe = Matmul(Ye,adjYeTeadjYe) 
 YeadjYeTeadjTe = Matmul(Ye,adjYeTeadjTe) 
 YeadjTeTeadjYe = Matmul(Ye,adjTeTeadjYe) 
 YeCYvmv2lam = Matmul(Ye,CYvmv2lam) 
 YeCYvTpYvadjYe = Matmul(Ye,CYvTpYvadjYe) 
Forall(i2=1:3)  YeCYvTpYvadjYe(i2,i2) =  Real(YeCYvTpYvadjYe(i2,i2),dp) 
 YeCYvTpTvadjTe = Matmul(Ye,CYvTpTvadjTe) 
 YeCTvTpTvadjYe = Matmul(Ye,CTvTpTvadjYe) 
Forall(i2=1:3)  YeCTvTpTvadjYe(i2,i2) =  Real(YeCTvTpTvadjYe(i2,i2),dp) 
 Yumq2adjYuYu = Matmul(Yu,mq2adjYuYu) 
 YuadjYdYdadjYu = Matmul(Yu,adjYdYdadjYu) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYdTdadjYu = Matmul(Yu,adjYdTdadjYu) 
 YuadjYdTdadjTu = Matmul(Yu,adjYdTdadjTu) 
 YuadjYumu2Yu = Matmul(Yu,adjYumu2Yu) 
 YuadjYuYumq2 = Matmul(Yu,adjYuYumq2) 
 YuadjYuYuadjYu = Matmul(Yu,adjYuYuadjYu) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YuadjYuTuadjYu = Matmul(Yu,adjYuTuadjYu) 
 YuadjYuTuadjTu = Matmul(Yu,adjYuTuadjTu) 
 YuadjTdTdadjYu = Matmul(Yu,adjTdTdadjYu) 
 YuadjTuTuadjYu = Matmul(Yu,adjTuTuadjYu) 
 YvadjYvYvadjYv = Matmul(Yv,adjYvYvadjYv) 
Forall(i2=1:3)  YvadjYvYvadjYv(i2,i2) =  Real(YvadjYvYvadjYv(i2,i2),dp) 
 YvadjYvYvClam = Matmul(Yv,adjYvYvClam) 
 YvadjYvTvadjYv = Matmul(Yv,adjYvTvadjYv) 
 YvadjYvTvadjTv = Matmul(Yv,adjYvTvadjTv) 
 YvadjYvTvClam = Matmul(Yv,adjYvTvClam) 
 YvadjYvTpYeCYe = Matmul(Yv,adjYvTpYeCYe) 
 YvadjYvTpTeCTe = Matmul(Yv,adjYvTpTeCTe) 
 YvadjTvTvadjYv = Matmul(Yv,adjTvTvadjYv) 
 Yvadjkap1mv2kap1 = Matmul(Yv,adjkap1mv2kap1) 
 Yvadjkap1mv2kap2 = Matmul(Yv,adjkap1mv2kap2) 
 Yvadjkap1mv2kap3 = Matmul(Yv,adjkap1mv2kap3) 
 Yvadjkap2mv2kap1 = Matmul(Yv,adjkap2mv2kap1) 
 Yvadjkap2mv2kap2 = Matmul(Yv,adjkap2mv2kap2) 
 Yvadjkap2mv2kap3 = Matmul(Yv,adjkap2mv2kap3) 
 Yvadjkap3mv2kap1 = Matmul(Yv,adjkap3mv2kap1) 
 Yvadjkap3mv2kap2 = Matmul(Yv,adjkap3mv2kap2) 
 Yvadjkap3mv2kap3 = Matmul(Yv,adjkap3mv2kap3) 
 YvCmv2Ckap1kap1 = Matmul(Yv,Cmv2Ckap1kap1) 
 YvCmv2Ckap1kap2 = Matmul(Yv,Cmv2Ckap1kap2) 
 YvCmv2Ckap1kap3 = Matmul(Yv,Cmv2Ckap1kap3) 
 YvCmv2Ckap2kap1 = Matmul(Yv,Cmv2Ckap2kap1) 
 YvCmv2Ckap2kap2 = Matmul(Yv,Cmv2Ckap2kap2) 
 YvCmv2Ckap2kap3 = Matmul(Yv,Cmv2Ckap2kap3) 
 YvCmv2Ckap3kap1 = Matmul(Yv,Cmv2Ckap3kap1) 
 YvCmv2Ckap3kap2 = Matmul(Yv,Cmv2Ckap3kap2) 
 YvCmv2Ckap3kap3 = Matmul(Yv,Cmv2Ckap3kap3) 
 YvCkap1Tpkap1adjYv = Matmul(Yv,Ckap1Tpkap1adjYv) 
Forall(i2=1:3)  YvCkap1Tpkap1adjYv(i2,i2) =  Real(YvCkap1Tpkap1adjYv(i2,i2),dp) 
 YvCkap1Tpkap1Clam = Matmul(Yv,Ckap1Tpkap1Clam) 
 YvCkap1TpTk1adjYv = Matmul(Yv,Ckap1TpTk1adjYv) 
 YvCkap1kap1adjYv = Matmul(Yv,Ckap1kap1adjYv) 
 YvCkap1Tk1adjTv = Matmul(Yv,Ckap1Tk1adjTv) 
 YvCkap2Tpkap2adjYv = Matmul(Yv,Ckap2Tpkap2adjYv) 
Forall(i2=1:3)  YvCkap2Tpkap2adjYv(i2,i2) =  Real(YvCkap2Tpkap2adjYv(i2,i2),dp) 
 YvCkap2Tpkap2Clam = Matmul(Yv,Ckap2Tpkap2Clam) 
 YvCkap2TpTk2adjYv = Matmul(Yv,Ckap2TpTk2adjYv) 
 YvCkap2kap2adjYv = Matmul(Yv,Ckap2kap2adjYv) 
 YvCkap2Tk2adjTv = Matmul(Yv,Ckap2Tk2adjTv) 
 YvCkap3Tpkap3adjYv = Matmul(Yv,Ckap3Tpkap3adjYv) 
Forall(i2=1:3)  YvCkap3Tpkap3adjYv(i2,i2) =  Real(YvCkap3Tpkap3adjYv(i2,i2),dp) 
 YvCkap3Tpkap3Clam = Matmul(Yv,Ckap3Tpkap3Clam) 
 YvCkap3TpTk3adjYv = Matmul(Yv,Ckap3TpTk3adjYv) 
 YvCkap3kap3adjYv = Matmul(Yv,Ckap3kap3adjYv) 
 YvCkap3Tk3adjTv = Matmul(Yv,Ckap3Tk3adjTv) 
 YvCTk1TpTk1adjYv = Matmul(Yv,CTk1TpTk1adjYv) 
Forall(i2=1:3)  YvCTk1TpTk1adjYv(i2,i2) =  Real(YvCTk1TpTk1adjYv(i2,i2),dp) 
 YvCTk2TpTk2adjYv = Matmul(Yv,CTk2TpTk2adjYv) 
Forall(i2=1:3)  YvCTk2TpTk2adjYv(i2,i2) =  Real(YvCTk2TpTk2adjYv(i2,i2),dp) 
 YvCTk3TpTk3adjYv = Matmul(Yv,CTk3TpTk3adjYv) 
Forall(i2=1:3)  YvCTk3TpTk3adjYv(i2,i2) =  Real(YvCTk3TpTk3adjYv(i2,i2),dp) 
 adjYdmd2YdadjYd = Matmul(adjYd,md2YdadjYd) 
 adjYdmd2YdadjYu = Matmul(adjYd,md2YdadjYu) 
 adjYdYdmq2adjYd = Matmul(adjYd,Ydmq2adjYd) 
 adjYdYdmq2adjYu = Matmul(adjYd,Ydmq2adjYu) 
 adjYdYdadjYdmd2 = Matmul(adjYd,YdadjYdmd2) 
 adjYdYdadjYdYd = Matmul(adjYd,YdadjYdYd) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdTd = Matmul(adjYd,YdadjYdTd) 
 adjYdYdadjYumu2 = Matmul(adjYd,YdadjYumu2) 
 adjYdYdadjYuYu = Matmul(adjYd,YdadjYuYu) 
 adjYdYdadjYuTu = Matmul(adjYd,YdadjYuTu) 
 adjYdYdadjTdTd = Matmul(adjYd,YdadjTdTd) 
 adjYdTdadjYdYd = Matmul(adjYd,TdadjYdYd) 
 adjYdTdadjYuYu = Matmul(adjYd,TdadjYuYu) 
 adjYdTdadjTdYd = Matmul(adjYd,TdadjTdYd) 
 adjYeme2YeadjYe = Matmul(adjYe,me2YeadjYe) 
 adjYeme2YeCYv = Matmul(adjYe,me2YeCYv) 
 adjYeYeml2adjYe = Matmul(adjYe,Yeml2adjYe) 
 adjYeYeml2CYv = Matmul(adjYe,Yeml2CYv) 
 adjYeYeadjYeme2 = Matmul(adjYe,YeadjYeme2) 
 adjYeYeadjYeYe = Matmul(adjYe,YeadjYeYe) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYeTe = Matmul(adjYe,YeadjYeTe) 
 adjYeYeadjTeTe = Matmul(adjYe,YeadjTeTe) 
 adjYeYeCYvmv2 = Matmul(adjYe,YeCYvmv2) 
 adjYeYeCYvvR = Matmul(adjYe,YeCYvvR) 
 adjYeYeCYvlam = Matmul(adjYe,YeCYvlam) 
 adjYeYeCTvTlam = Matmul(adjYe,YeCTvTlam) 
 adjYeTeadjYeYe = Matmul(adjYe,TeadjYeYe) 
 adjYeTeadjTeYe = Matmul(adjYe,TeadjTeYe) 
 adjYeTeCYvlam = Matmul(adjYe,TeCYvlam) 
 adjYeTeCYvTpYv = Matmul(adjYe,TeCYvTpYv) 
 adjYeTeCTvlam = Matmul(adjYe,TeCTvlam) 
 adjYeTeCTvTpYv = Matmul(adjYe,TeCTvTpYv) 
 adjYumu2YuadjYd = Matmul(adjYu,mu2YuadjYd) 
 adjYumu2YuadjYu = Matmul(adjYu,mu2YuadjYu) 
 adjYuYumq2adjYd = Matmul(adjYu,Yumq2adjYd) 
 adjYuYumq2adjYu = Matmul(adjYu,Yumq2adjYu) 
 adjYuYuadjYdmd2 = Matmul(adjYu,YuadjYdmd2) 
 adjYuYuadjYdYd = Matmul(adjYu,YuadjYdYd) 
 adjYuYuadjYdTd = Matmul(adjYu,YuadjYdTd) 
 adjYuYuadjYumu2 = Matmul(adjYu,YuadjYumu2) 
 adjYuYuadjYuYu = Matmul(adjYu,YuadjYuYu) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYuYuadjYuTu = Matmul(adjYu,YuadjYuTu) 
 adjYuYuadjTuTu = Matmul(adjYu,YuadjTuTu) 
 adjYuTuadjYdYd = Matmul(adjYu,TuadjYdYd) 
 adjYuTuadjYuYu = Matmul(adjYu,TuadjYuYu) 
 adjYuTuadjTuYu = Matmul(adjYu,TuadjTuYu) 
 adjYvYvadjYvmlHd2 = Matmul(adjYv,YvadjYvmlHd2) 
 adjYvYvadjYvvL = Matmul(adjYv,YvadjYvvL) 
 adjYvYvadjYvYv = Matmul(adjYv,YvadjYvYv) 
Forall(i2=1:3)  adjYvYvadjYvYv(i2,i2) =  Real(adjYvYvadjYvYv(i2,i2),dp) 
 adjYvYvadjYvTv = Matmul(adjYv,YvadjYvTv) 
 adjYvYvadjkap1mv2 = Matmul(adjYv,Yvadjkap1mv2) 
 adjYvYvadjkap2mv2 = Matmul(adjYv,Yvadjkap2mv2) 
 adjYvYvadjkap3mv2 = Matmul(adjYv,Yvadjkap3mv2) 
 adjYvYvCmv2adjkap1 = Matmul(adjYv,YvCmv2adjkap1) 
 adjYvYvCmv2adjkap2 = Matmul(adjYv,YvCmv2adjkap2) 
 adjYvYvCmv2adjkap3 = Matmul(adjYv,YvCmv2adjkap3) 
 adjYvYvCmv2Clam = Matmul(adjYv,YvCmv2Clam) 
 adjYvYvCkap1lam = Matmul(adjYv,YvCkap1lam) 
 adjYvYvCkap1kap2 = Matmul(adjYv,YvCkap1kap2) 
 adjYvYvCkap1kap3 = Matmul(adjYv,YvCkap1kap3) 
 adjYvYvCkap2lam = Matmul(adjYv,YvCkap2lam) 
 adjYvYvCkap2kap1 = Matmul(adjYv,YvCkap2kap1) 
 adjYvYvCkap2kap3 = Matmul(adjYv,YvCkap2kap3) 
 adjYvYvCkap3lam = Matmul(adjYv,YvCkap3lam) 
 adjYvYvCkap3kap1 = Matmul(adjYv,YvCkap3kap1) 
 adjYvYvCkap3kap2 = Matmul(adjYv,YvCkap3kap2) 
 adjYvCml2YvadjYv = Matmul(adjYv,Cml2YvadjYv) 
 adjYvCml2Yvadjkap1 = Matmul(adjYv,Cml2Yvadjkap1) 
 adjYvCml2Yvadjkap2 = Matmul(adjYv,Cml2Yvadjkap2) 
 adjYvCml2Yvadjkap3 = Matmul(adjYv,Cml2Yvadjkap3) 
 adjYvCml2YvClam = Matmul(adjYv,Cml2YvClam) 
 adjYvCml2TpYeCYe = Matmul(adjYv,Cml2TpYeCYe) 
 adjYvTvadjYvYv = Matmul(adjYv,TvadjYvYv) 
 adjYvTvadjkap1kap1 = Matmul(adjYv,Tvadjkap1kap1) 
 adjYvTvadjkap1kap2 = Matmul(adjYv,Tvadjkap1kap2) 
 adjYvTvadjkap1kap3 = Matmul(adjYv,Tvadjkap1kap3) 
 adjYvTvadjkap2kap1 = Matmul(adjYv,Tvadjkap2kap1) 
 adjYvTvadjkap2kap2 = Matmul(adjYv,Tvadjkap2kap2) 
 adjYvTvadjkap2kap3 = Matmul(adjYv,Tvadjkap2kap3) 
 adjYvTvadjkap3kap1 = Matmul(adjYv,Tvadjkap3kap1) 
 adjYvTvadjkap3kap2 = Matmul(adjYv,Tvadjkap3kap2) 
 adjYvTvadjkap3kap3 = Matmul(adjYv,Tvadjkap3kap3) 
 adjYvTpYeCme2CYe = Matmul(adjYv,TpYeCme2CYe) 
 adjYvTpYeCYeYv = Matmul(adjYv,TpYeCYeYv) 
Forall(i2=1:3)  adjYvTpYeCYeYv(i2,i2) =  Real(adjYvTpYeCYeYv(i2,i2),dp) 
 adjYvTpYeCYeCml2 = Matmul(adjYv,TpYeCYeCml2) 
 adjYvTpYeCYeTv = Matmul(adjYv,TpYeCYeTv) 
 adjYvTpYeCTeTv = Matmul(adjYv,TpYeCTeTv) 
 adjYvTpTeCYeYv = Matmul(adjYv,TpTeCYeYv) 
 adjTdYdadjYdTd = Matmul(adjTd,YdadjYdTd) 
 adjTdTdadjYdYd = Matmul(adjTd,TdadjYdYd) 
 adjTeYeadjYeTe = Matmul(adjTe,YeadjYeTe) 
 adjTeYeCYvTlam = Matmul(adjTe,YeCYvTlam) 
 adjTeTeadjYeYe = Matmul(adjTe,TeadjYeYe) 
 adjTeTeCYvlam = Matmul(adjTe,TeCYvlam) 
 adjTuYuadjYuTu = Matmul(adjTu,YuadjYuTu) 
 adjTuTuadjYuYu = Matmul(adjTu,TuadjYuYu) 
 adjTvYvCkap1Tk2 = Matmul(adjTv,YvCkap1Tk2) 
 adjTvYvCkap1Tk3 = Matmul(adjTv,YvCkap1Tk3) 
 adjTvYvCkap2Tk1 = Matmul(adjTv,YvCkap2Tk1) 
 adjTvYvCkap2Tk3 = Matmul(adjTv,YvCkap2Tk3) 
 adjTvYvCkap3Tk1 = Matmul(adjTv,YvCkap3Tk1) 
 adjTvYvCkap3Tk2 = Matmul(adjTv,YvCkap3Tk2) 
 adjkap1Tpkap1Cmv2Clam = Matmul(adjkap1,Tpkap1Cmv2Clam) 
 adjkap1Tpkap2Cmv2Clam = Matmul(adjkap1,Tpkap2Cmv2Clam) 
 adjkap1Tpkap3Cmv2Clam = Matmul(adjkap1,Tpkap3Cmv2Clam) 
 adjkap2Tpkap1Cmv2Clam = Matmul(adjkap2,Tpkap1Cmv2Clam) 
 adjkap2Tpkap2Cmv2Clam = Matmul(adjkap2,Tpkap2Cmv2Clam) 
 adjkap2Tpkap3Cmv2Clam = Matmul(adjkap2,Tpkap3Cmv2Clam) 
 adjkap3Tpkap1Cmv2Clam = Matmul(adjkap3,Tpkap1Cmv2Clam) 
 adjkap3Tpkap2Cmv2Clam = Matmul(adjkap3,Tpkap2Cmv2Clam) 
 adjkap3Tpkap3Cmv2Clam = Matmul(adjkap3,Tpkap3Cmv2Clam) 
 Cml2YvCkap1kap2 = Matmul(Conjg(ml2),YvCkap1kap2) 
 Cml2YvCkap1kap3 = Matmul(Conjg(ml2),YvCkap1kap3) 
 Cml2YvCkap2kap1 = Matmul(Conjg(ml2),YvCkap2kap1) 
 Cml2YvCkap2kap3 = Matmul(Conjg(ml2),YvCkap2kap3) 
 Cml2YvCkap3kap1 = Matmul(Conjg(ml2),YvCkap3kap1) 
 Cml2YvCkap3kap2 = Matmul(Conjg(ml2),YvCkap3kap2) 
 Cmv2adjYvYvadjYv = Matmul(Conjg(mv2),adjYvYvadjYv) 
 Cmv2adjYvYvadjkap1 = Matmul(Conjg(mv2),adjYvYvadjkap1) 
 Cmv2adjYvYvadjkap2 = Matmul(Conjg(mv2),adjYvYvadjkap2) 
 Cmv2adjYvYvadjkap3 = Matmul(Conjg(mv2),adjYvYvadjkap3) 
 Cmv2adjYvYvClam = Matmul(Conjg(mv2),adjYvYvClam) 
 Cmv2adjYvTpYeCYe = Matmul(Conjg(mv2),adjYvTpYeCYe) 
 Cmv2Ckap1kap1adjYv = Matmul(Conjg(mv2),Ckap1kap1adjYv) 
 Cmv2Ckap2kap2adjYv = Matmul(Conjg(mv2),Ckap2kap2adjYv) 
 Cmv2Ckap3kap3adjYv = Matmul(Conjg(mv2),Ckap3kap3adjYv) 
 CYeTpYeCYevL = Matmul(Conjg(Ye),TpYeCYevL) 
 CYeTpYeCYeYv = Matmul(Conjg(Ye),TpYeCYeYv) 
 CYeTpYeCYeTv = Matmul(Conjg(Ye),TpYeCYeTv) 
 CYeTpTeCYeYv = Matmul(Conjg(Ye),TpTeCYeYv) 
 CYvmv2TpYvadjYe = Matmul(Conjg(Yv),mv2TpYvadjYe) 
 CYvmv2TpYvCYv = Matmul(Conjg(Yv),mv2TpYvCYv) 
 CYvmv2kap1Ckap1 = Matmul(Conjg(Yv),mv2kap1Ckap1) 
 CYvmv2kap1Ckap2 = Matmul(Conjg(Yv),mv2kap1Ckap2) 
 CYvmv2kap1Ckap3 = Matmul(Conjg(Yv),mv2kap1Ckap3) 
 CYvmv2kap2Ckap1 = Matmul(Conjg(Yv),mv2kap2Ckap1) 
 CYvmv2kap2Ckap2 = Matmul(Conjg(Yv),mv2kap2Ckap2) 
 CYvmv2kap2Ckap3 = Matmul(Conjg(Yv),mv2kap2Ckap3) 
 CYvmv2kap3Ckap1 = Matmul(Conjg(Yv),mv2kap3Ckap1) 
 CYvmv2kap3Ckap2 = Matmul(Conjg(Yv),mv2kap3Ckap2) 
 CYvmv2kap3Ckap3 = Matmul(Conjg(Yv),mv2kap3Ckap3) 
 CYvTpYvml2adjYe = Matmul(Conjg(Yv),TpYvml2adjYe) 
 CYvTpYvml2CYv = Matmul(Conjg(Yv),TpYvml2CYv) 
 CYvTpYvadjYeme2 = Matmul(Conjg(Yv),TpYvadjYeme2) 
 CYvTpYvadjYeYe = Matmul(Conjg(Yv),TpYvadjYeYe) 
 CYvTpYvadjYeTe = Matmul(Conjg(Yv),TpYvadjYeTe) 
 CYvTpYvCYvmv2 = Matmul(Conjg(Yv),TpYvCYvmv2) 
 CYvTpYvCYvvR = Matmul(Conjg(Yv),TpYvCYvvR) 
 CYvTpYvCYvlam = Matmul(Conjg(Yv),TpYvCYvlam) 
 CYvTpYvCYvTlam = Matmul(Conjg(Yv),TpYvCYvTlam) 
 CYvTpYvCYvTpYv = Matmul(Conjg(Yv),TpYvCYvTpYv) 
Forall(i2=1:3)  CYvTpYvCYvTpYv(i2,i2) =  Real(CYvTpYvCYvTpYv(i2,i2),dp) 
 CYvTpYvCYvTpTv = Matmul(Conjg(Yv),TpYvCYvTpTv) 
 CYvTpYvCTvTlam = Matmul(Conjg(Yv),TpYvCTvTlam) 
 CYvTpYvCTvTpTv = Matmul(Conjg(Yv),TpYvCTvTpTv) 
 CYvTpTvadjYeYe = Matmul(Conjg(Yv),TpTvadjYeYe) 
 CYvTpTvCYvlam = Matmul(Conjg(Yv),TpTvCYvlam) 
 CYvTpTvCYvTpYv = Matmul(Conjg(Yv),TpTvCYvTpYv) 
 CYvTpTvCTvlam = Matmul(Conjg(Yv),TpTvCTvlam) 
 CYvTpTvCTvTpYv = Matmul(Conjg(Yv),TpTvCTvTpYv) 
Forall(i2=1:3)  CYvTpTvCTvTpYv(i2,i2) =  Real(CYvTpTvCTvTpYv(i2,i2),dp) 
 CYvkap1adjkap1lam = Matmul(Conjg(Yv),kap1adjkap1lam) 
 CYvkap1adjkap1TpYv = Matmul(Conjg(Yv),kap1adjkap1TpYv) 
 CYvkap2adjkap2lam = Matmul(Conjg(Yv),kap2adjkap2lam) 
 CYvkap2adjkap2TpYv = Matmul(Conjg(Yv),kap2adjkap2TpYv) 
 CYvkap3adjkap3lam = Matmul(Conjg(Yv),kap3adjkap3lam) 
 CYvkap3adjkap3TpYv = Matmul(Conjg(Yv),kap3adjkap3TpYv) 
 CYvTk1adjkap1lam = Matmul(Conjg(Yv),Tk1adjkap1lam) 
 CYvTk1adjkap1TpYv = Matmul(Conjg(Yv),Tk1adjkap1TpYv) 
 CYvTk1adjTk1lam = Matmul(Conjg(Yv),Tk1adjTk1lam) 
 CYvTk1adjTk1TpYv = Matmul(Conjg(Yv),Tk1adjTk1TpYv) 
 CYvTk2adjkap2lam = Matmul(Conjg(Yv),Tk2adjkap2lam) 
 CYvTk2adjkap2TpYv = Matmul(Conjg(Yv),Tk2adjkap2TpYv) 
 CYvTk2adjTk2lam = Matmul(Conjg(Yv),Tk2adjTk2lam) 
 CYvTk2adjTk2TpYv = Matmul(Conjg(Yv),Tk2adjTk2TpYv) 
 CYvTk3adjkap3lam = Matmul(Conjg(Yv),Tk3adjkap3lam) 
 CYvTk3adjkap3TpYv = Matmul(Conjg(Yv),Tk3adjkap3TpYv) 
 CYvTk3adjTk3lam = Matmul(Conjg(Yv),Tk3adjTk3lam) 
 CYvTk3adjTk3TpYv = Matmul(Conjg(Yv),Tk3adjTk3TpYv) 
 CTvTpYvCYvTlam = Matmul(Conjg(Tv),TpYvCYvTlam) 
 CTvTpYvCYvTpTv = Matmul(Conjg(Tv),TpYvCYvTpTv) 
Forall(i2=1:3)  CTvTpYvCYvTpTv(i2,i2) =  Real(CTvTpYvCYvTpTv(i2,i2),dp) 
 CTvTpTvCYvlam = Matmul(Conjg(Tv),TpTvCYvlam) 
 CTvTpTvCYvTpYv = Matmul(Conjg(Tv),TpTvCYvTpYv) 
 Ckap1TpYvCYvTpTk1 = Matmul(Conjg(kap1),TpYvCYvTpTk1) 
 Ckap1Tpkap1adjYvvL = Matmul(Conjg(kap1),Tpkap1adjYvvL) 
 Ckap1Tpkap1adjYvYv = Matmul(Conjg(kap1),Tpkap1adjYvYv) 
 Ckap1Tpkap1Cmv2adjYv = Matmul(Conjg(kap1),Tpkap1Cmv2adjYv) 
 Ckap1Tpkap1Cmv2adjkap1 = Matmul(Conjg(kap1),Tpkap1Cmv2adjkap1) 
 Ckap1Tpkap1Cmv2adjkap2 = Matmul(Conjg(kap1),Tpkap1Cmv2adjkap2) 
 Ckap1Tpkap1Cmv2adjkap3 = Matmul(Conjg(kap1),Tpkap1Cmv2adjkap3) 
 Ckap1TpTk1adjYvYv = Matmul(Conjg(kap1),TpTk1adjYvYv) 
 Ckap1kap1adjYvCml2 = Matmul(Conjg(kap1),kap1adjYvCml2) 
 Ckap2TpYvCYvTpTk2 = Matmul(Conjg(kap2),TpYvCYvTpTk2) 
 Ckap2Tpkap2adjYvvL = Matmul(Conjg(kap2),Tpkap2adjYvvL) 
 Ckap2Tpkap2adjYvYv = Matmul(Conjg(kap2),Tpkap2adjYvYv) 
 Ckap2Tpkap2Cmv2adjYv = Matmul(Conjg(kap2),Tpkap2Cmv2adjYv) 
 Ckap2Tpkap2Cmv2adjkap1 = Matmul(Conjg(kap2),Tpkap2Cmv2adjkap1) 
 Ckap2Tpkap2Cmv2adjkap2 = Matmul(Conjg(kap2),Tpkap2Cmv2adjkap2) 
 Ckap2Tpkap2Cmv2adjkap3 = Matmul(Conjg(kap2),Tpkap2Cmv2adjkap3) 
 Ckap2TpTk2adjYvYv = Matmul(Conjg(kap2),TpTk2adjYvYv) 
 Ckap2kap2adjYvCml2 = Matmul(Conjg(kap2),kap2adjYvCml2) 
 Ckap3TpYvCYvTpTk3 = Matmul(Conjg(kap3),TpYvCYvTpTk3) 
 Ckap3Tpkap3adjYvvL = Matmul(Conjg(kap3),Tpkap3adjYvvL) 
 Ckap3Tpkap3adjYvYv = Matmul(Conjg(kap3),Tpkap3adjYvYv) 
 Ckap3Tpkap3Cmv2adjYv = Matmul(Conjg(kap3),Tpkap3Cmv2adjYv) 
 Ckap3Tpkap3Cmv2adjkap1 = Matmul(Conjg(kap3),Tpkap3Cmv2adjkap1) 
 Ckap3Tpkap3Cmv2adjkap2 = Matmul(Conjg(kap3),Tpkap3Cmv2adjkap2) 
 Ckap3Tpkap3Cmv2adjkap3 = Matmul(Conjg(kap3),Tpkap3Cmv2adjkap3) 
 Ckap3TpTk3adjYvYv = Matmul(Conjg(kap3),TpTk3adjYvYv) 
 Ckap3kap3adjYvCml2 = Matmul(Conjg(kap3),kap3adjYvCml2) 
 TdadjYdYdadjTd = Matmul(Td,adjYdYdadjTd) 
 TdadjYuYuadjTd = Matmul(Td,adjYuYuadjTd) 
 TdadjTdYdadjYd = Matmul(Td,adjTdYdadjYd) 
 TdadjTuYuadjYd = Matmul(Td,adjTuYuadjYd) 
 TeadjYeYeadjTe = Matmul(Te,adjYeYeadjTe) 
 TeadjTeYeadjYe = Matmul(Te,adjTeYeadjYe) 
 TeCYvTpYvadjTe = Matmul(Te,CYvTpYvadjTe) 
 TeCTvTpYvadjYe = Matmul(Te,CTvTpYvadjYe) 
 TuadjYdYdadjTu = Matmul(Tu,adjYdYdadjTu) 
 TuadjYuYuadjTu = Matmul(Tu,adjYuYuadjTu) 
 TuadjTdYdadjYu = Matmul(Tu,adjTdYdadjYu) 
 TuadjTuYuadjYu = Matmul(Tu,adjTuYuadjYu) 
 TvadjYvYvClam = Matmul(Tv,adjYvYvClam) 
 TpYeCYeYvClam = Matmul(Transpose(Ye),CYeYvClam) 
 TpYeCYeTvClam = Matmul(Transpose(Ye),CYeTvClam) 
 TpYvml2CYvlam = Matmul(Transpose(Yv),ml2CYvlam) 
 TpYvml2CYvTpYv = Matmul(Transpose(Yv),ml2CYvTpYv) 
 TpYvadjYeYeCYv = Matmul(Transpose(Yv),adjYeYeCYv) 
Forall(i2=1:3)  TpYvadjYeYeCYv(i2,i2) =  Real(TpYvadjYeYeCYv(i2,i2),dp) 
 TpYvadjYeTeCYv = Matmul(Transpose(Yv),adjYeTeCYv) 
 TpYvadjYeTeCTv = Matmul(Transpose(Yv),adjYeTeCTv) 
 TpYvadjTeTeCYv = Matmul(Transpose(Yv),adjTeTeCYv) 
 TpYvCYvmv2lam = Matmul(Transpose(Yv),CYvmv2lam) 
 TpYvCYvmv2TpYv = Matmul(Transpose(Yv),CYvmv2TpYv) 
 TpYvCYvTpYvml2 = Matmul(Transpose(Yv),CYvTpYvml2) 
 TpYvCYvTpYvCYv = Matmul(Transpose(Yv),CYvTpYvCYv) 
Forall(i2=1:3)  TpYvCYvTpYvCYv(i2,i2) =  Real(TpYvCYvTpYvCYv(i2,i2),dp) 
 TpYvCYvTpTvCTv = Matmul(Transpose(Yv),CYvTpTvCTv) 
 TpYvCTvTpTvCYv = Matmul(Transpose(Yv),CTvTpTvCYv) 
Forall(i2=1:3)  TpYvCTvTpTvCYv(i2,i2) =  Real(TpYvCTvTpTvCYv(i2,i2),dp) 
 TpTeCYeYvClam = Matmul(Transpose(Te),CYeYvClam) 
 TpTvadjYeYeCTv = Matmul(Transpose(Tv),adjYeYeCTv) 
Forall(i2=1:3)  TpTvadjYeYeCTv(i2,i2) =  Real(TpTvadjYeYeCTv(i2,i2),dp) 
 TpTvadjTeYeCYv = Matmul(Transpose(Tv),adjTeYeCYv) 
 TpTvCYvTpYvCYv = Matmul(Transpose(Tv),CYvTpYvCYv) 
 TpTvCYvTpYvCTv = Matmul(Transpose(Tv),CYvTpYvCTv) 
Forall(i2=1:3)  TpTvCYvTpYvCTv(i2,i2) =  Real(TpTvCYvTpYvCTv(i2,i2),dp) 
 TpTvCTvTpYvCYv = Matmul(Transpose(Tv),CTvTpYvCYv) 
 Tpkap1adjkap1mv2lam = Matmul(Transpose(kap1),adjkap1mv2lam) 
 Tpkap1adjkap2mv2lam = Matmul(Transpose(kap1),adjkap2mv2lam) 
 Tpkap1adjkap3mv2lam = Matmul(Transpose(kap1),adjkap3mv2lam) 
 Tpkap2adjkap1mv2lam = Matmul(Transpose(kap2),adjkap1mv2lam) 
 Tpkap2adjkap2mv2lam = Matmul(Transpose(kap2),adjkap2mv2lam) 
 Tpkap2adjkap3mv2lam = Matmul(Transpose(kap2),adjkap3mv2lam) 
 Tpkap3adjkap1mv2lam = Matmul(Transpose(kap3),adjkap1mv2lam) 
 Tpkap3adjkap2mv2lam = Matmul(Transpose(kap3),adjkap2mv2lam) 
 Tpkap3adjkap3mv2lam = Matmul(Transpose(kap3),adjkap3mv2lam) 
 kap1adjYvYvadjkap1 = Matmul(kap1,adjYvYvadjkap1) 
 kap1adjYvYvCkap1 = Matmul(kap1,adjYvYvCkap1) 
 kap1adjYvYvCkap2 = Matmul(kap1,adjYvYvCkap2) 
 kap1adjYvYvCkap3 = Matmul(kap1,adjYvYvCkap3) 
 kap1adjkap1TpYvml2 = Matmul(kap1,adjkap1TpYvml2) 
 kap1Cmv2Ckap1kap1 = Matmul(kap1,Cmv2Ckap1kap1) 
 kap1Cmv2Ckap1kap2 = Matmul(kap1,Cmv2Ckap1kap2) 
 kap1Cmv2Ckap1kap3 = Matmul(kap1,Cmv2Ckap1kap3) 
 kap1Cmv2Ckap2kap1 = Matmul(kap1,Cmv2Ckap2kap1) 
 kap1Cmv2Ckap2kap2 = Matmul(kap1,Cmv2Ckap2kap2) 
 kap1Cmv2Ckap2kap3 = Matmul(kap1,Cmv2Ckap2kap3) 
 kap1Cmv2Ckap3kap1 = Matmul(kap1,Cmv2Ckap3kap1) 
 kap1Cmv2Ckap3kap2 = Matmul(kap1,Cmv2Ckap3kap2) 
 kap1Cmv2Ckap3kap3 = Matmul(kap1,Cmv2Ckap3kap3) 
 kap2adjYvYvadjkap2 = Matmul(kap2,adjYvYvadjkap2) 
 kap2adjYvYvCkap1 = Matmul(kap2,adjYvYvCkap1) 
 kap2adjYvYvCkap2 = Matmul(kap2,adjYvYvCkap2) 
 kap2adjYvYvCkap3 = Matmul(kap2,adjYvYvCkap3) 
 kap2adjkap2TpYvml2 = Matmul(kap2,adjkap2TpYvml2) 
 kap2Cmv2Ckap1kap1 = Matmul(kap2,Cmv2Ckap1kap1) 
 kap2Cmv2Ckap1kap2 = Matmul(kap2,Cmv2Ckap1kap2) 
 kap2Cmv2Ckap1kap3 = Matmul(kap2,Cmv2Ckap1kap3) 
 kap2Cmv2Ckap2kap1 = Matmul(kap2,Cmv2Ckap2kap1) 
 kap2Cmv2Ckap2kap2 = Matmul(kap2,Cmv2Ckap2kap2) 
 kap2Cmv2Ckap2kap3 = Matmul(kap2,Cmv2Ckap2kap3) 
 kap2Cmv2Ckap3kap1 = Matmul(kap2,Cmv2Ckap3kap1) 
 kap2Cmv2Ckap3kap2 = Matmul(kap2,Cmv2Ckap3kap2) 
 kap2Cmv2Ckap3kap3 = Matmul(kap2,Cmv2Ckap3kap3) 
 kap3adjYvYvadjkap3 = Matmul(kap3,adjYvYvadjkap3) 
 kap3adjYvYvCkap1 = Matmul(kap3,adjYvYvCkap1) 
 kap3adjYvYvCkap2 = Matmul(kap3,adjYvYvCkap2) 
 kap3adjYvYvCkap3 = Matmul(kap3,adjYvYvCkap3) 
 kap3adjkap3TpYvml2 = Matmul(kap3,adjkap3TpYvml2) 
 kap3Cmv2Ckap1kap1 = Matmul(kap3,Cmv2Ckap1kap1) 
 kap3Cmv2Ckap1kap2 = Matmul(kap3,Cmv2Ckap1kap2) 
 kap3Cmv2Ckap1kap3 = Matmul(kap3,Cmv2Ckap1kap3) 
 kap3Cmv2Ckap2kap1 = Matmul(kap3,Cmv2Ckap2kap1) 
 kap3Cmv2Ckap2kap2 = Matmul(kap3,Cmv2Ckap2kap2) 
 kap3Cmv2Ckap2kap3 = Matmul(kap3,Cmv2Ckap2kap3) 
 kap3Cmv2Ckap3kap1 = Matmul(kap3,Cmv2Ckap3kap1) 
 kap3Cmv2Ckap3kap2 = Matmul(kap3,Cmv2Ckap3kap2) 
 kap3Cmv2Ckap3kap3 = Matmul(kap3,Cmv2Ckap3kap3) 
 Tk1adjYvYvadjTk1 = Matmul(Tk1,adjYvYvadjTk1) 
 Tk1adjYvYvCkap1 = Matmul(Tk1,adjYvYvCkap1) 
 Tk2adjYvYvadjTk2 = Matmul(Tk2,adjYvYvadjTk2) 
 Tk2adjYvYvCkap2 = Matmul(Tk2,adjYvYvCkap2) 
 Tk3adjYvYvadjTk3 = Matmul(Tk3,adjYvYvadjTk3) 
 Tk3adjYvYvCkap3 = Matmul(Tk3,adjYvYvCkap3) 
 md2YdadjYdYdadjYd = Matmul(md2,YdadjYdYdadjYd) 
 md2YdadjYuYuadjYd = Matmul(md2,YdadjYuYuadjYd) 
 me2YeadjYeYeadjYe = Matmul(me2,YeadjYeYeadjYe) 
 me2YeCYvTpYvadjYe = Matmul(me2,YeCYvTpYvadjYe) 
 ml2adjYeYeadjYeYe = Matmul(ml2,adjYeYeadjYeYe) 
 ml2adjYeYeCYvlam = Matmul(ml2,adjYeYeCYvlam) 
 ml2CYvTpYvCYvlam = Matmul(ml2,CYvTpYvCYvlam) 
 ml2CYvTpYvCYvTpYv = Matmul(ml2,CYvTpYvCYvTpYv) 
 mq2adjYdYdadjYdYd = Matmul(mq2,adjYdYdadjYdYd) 
 mq2adjYdYdadjYuYu = Matmul(mq2,adjYdYdadjYuYu) 
 mq2adjYuYuadjYdYd = Matmul(mq2,adjYuYuadjYdYd) 
 mq2adjYuYuadjYuYu = Matmul(mq2,adjYuYuadjYuYu) 
 mu2YuadjYdYdadjYu = Matmul(mu2,YuadjYdYdadjYu) 
 mu2YuadjYuYuadjYu = Matmul(mu2,YuadjYuYuadjYu) 
 mv2TpYvadjYeYeCYv = Matmul(mv2,TpYvadjYeYeCYv) 
 mv2TpYvCYvTpYvCYv = Matmul(mv2,TpYvCYvTpYvCYv) 
 mv2kap1adjYvYvadjkap1 = Matmul(mv2,kap1adjYvYvadjkap1) 
 mv2kap2adjYvYvadjkap2 = Matmul(mv2,kap2adjYvYvadjkap2) 
 mv2kap3adjYvYvadjkap3 = Matmul(mv2,kap3adjYvYvadjkap3) 
 Ydmq2adjYdYdadjYd = Matmul(Yd,mq2adjYdYdadjYd) 
 Ydmq2adjYuYuadjYd = Matmul(Yd,mq2adjYuYuadjYd) 
 YdadjYdmd2YdadjYd = Matmul(Yd,adjYdmd2YdadjYd) 
Forall(i2=1:3)  YdadjYdmd2YdadjYd(i2,i2) =  Real(YdadjYdmd2YdadjYd(i2,i2),dp) 
 YdadjYdYdmq2adjYd = Matmul(Yd,adjYdYdmq2adjYd) 
 YdadjYdYdadjYdmd2 = Matmul(Yd,adjYdYdadjYdmd2) 
 YdadjYdYdadjYdYd = Matmul(Yd,adjYdYdadjYdYd) 
 YdadjYdYdadjYdTd = Matmul(Yd,adjYdYdadjYdTd) 
 YdadjYdTdadjYdYd = Matmul(Yd,adjYdTdadjYdYd) 
 YdadjYumu2YuadjYd = Matmul(Yd,adjYumu2YuadjYd) 
Forall(i2=1:3)  YdadjYumu2YuadjYd(i2,i2) =  Real(YdadjYumu2YuadjYd(i2,i2),dp) 
 YdadjYuYumq2adjYd = Matmul(Yd,adjYuYumq2adjYd) 
 YdadjYuYuadjYdmd2 = Matmul(Yd,adjYuYuadjYdmd2) 
 YdadjYuYuadjYdYd = Matmul(Yd,adjYuYuadjYdYd) 
 YdadjYuYuadjYdTd = Matmul(Yd,adjYuYuadjYdTd) 
 YdadjYuYuadjYuYu = Matmul(Yd,adjYuYuadjYuYu) 
 YdadjYuYuadjYuTu = Matmul(Yd,adjYuYuadjYuTu) 
 YdadjYuTuadjYdYd = Matmul(Yd,adjYuTuadjYdYd) 
 YdadjYuTuadjYuYu = Matmul(Yd,adjYuTuadjYuYu) 
 Yeml2adjYeYeadjYe = Matmul(Ye,ml2adjYeYeadjYe) 
 Yeml2CYvTpYvadjYe = Matmul(Ye,ml2CYvTpYvadjYe) 
 YeadjYeme2YeadjYe = Matmul(Ye,adjYeme2YeadjYe) 
Forall(i2=1:3)  YeadjYeme2YeadjYe(i2,i2) =  Real(YeadjYeme2YeadjYe(i2,i2),dp) 
 YeadjYeYeml2adjYe = Matmul(Ye,adjYeYeml2adjYe) 
 YeadjYeYeadjYeme2 = Matmul(Ye,adjYeYeadjYeme2) 
 YeadjYeYeadjYeYe = Matmul(Ye,adjYeYeadjYeYe) 
 YeadjYeYeadjYeTe = Matmul(Ye,adjYeYeadjYeTe) 
 YeadjYeTeadjYeYe = Matmul(Ye,adjYeTeadjYeYe) 
 YeCYvmv2TpYvadjYe = Matmul(Ye,CYvmv2TpYvadjYe) 
Forall(i2=1:3)  YeCYvmv2TpYvadjYe(i2,i2) =  Real(YeCYvmv2TpYvadjYe(i2,i2),dp) 
 YeCYvTpYvml2adjYe = Matmul(Ye,CYvTpYvml2adjYe) 
 YeCYvTpYvadjYeme2 = Matmul(Ye,CYvTpYvadjYeme2) 
 YeCYvTpYvadjYeYe = Matmul(Ye,CYvTpYvadjYeYe) 
 YeCYvTpYvadjYeTe = Matmul(Ye,CYvTpYvadjYeTe) 
 YeCYvTpYvCYvTpYv = Matmul(Ye,CYvTpYvCYvTpYv) 
 YeCYvTpYvCYvTpTv = Matmul(Ye,CYvTpYvCYvTpTv) 
 YeCYvTpTvadjYeYe = Matmul(Ye,CYvTpTvadjYeYe) 
 YeCYvTpTvCYvTpYv = Matmul(Ye,CYvTpTvCYvTpYv) 
 YeCYvkap1adjkap1TpYv = Matmul(Ye,CYvkap1adjkap1TpYv) 
 YeCYvkap2adjkap2TpYv = Matmul(Ye,CYvkap2adjkap2TpYv) 
 YeCYvkap3adjkap3TpYv = Matmul(Ye,CYvkap3adjkap3TpYv) 
 YeCYvTk1adjkap1TpYv = Matmul(Ye,CYvTk1adjkap1TpYv) 
 YeCYvTk2adjkap2TpYv = Matmul(Ye,CYvTk2adjkap2TpYv) 
 YeCYvTk3adjkap3TpYv = Matmul(Ye,CYvTk3adjkap3TpYv) 
 Yumq2adjYdYdadjYu = Matmul(Yu,mq2adjYdYdadjYu) 
 Yumq2adjYuYuadjYu = Matmul(Yu,mq2adjYuYuadjYu) 
 YuadjYdmd2YdadjYu = Matmul(Yu,adjYdmd2YdadjYu) 
Forall(i2=1:3)  YuadjYdmd2YdadjYu(i2,i2) =  Real(YuadjYdmd2YdadjYu(i2,i2),dp) 
 YuadjYdYdmq2adjYu = Matmul(Yu,adjYdYdmq2adjYu) 
 YuadjYdYdadjYdYd = Matmul(Yu,adjYdYdadjYdYd) 
 YuadjYdYdadjYdTd = Matmul(Yu,adjYdYdadjYdTd) 
 YuadjYdYdadjYumu2 = Matmul(Yu,adjYdYdadjYumu2) 
 YuadjYdYdadjYuYu = Matmul(Yu,adjYdYdadjYuYu) 
 YuadjYdYdadjYuTu = Matmul(Yu,adjYdYdadjYuTu) 
 YuadjYdTdadjYdYd = Matmul(Yu,adjYdTdadjYdYd) 
 YuadjYdTdadjYuYu = Matmul(Yu,adjYdTdadjYuYu) 
 YuadjYumu2YuadjYu = Matmul(Yu,adjYumu2YuadjYu) 
Forall(i2=1:3)  YuadjYumu2YuadjYu(i2,i2) =  Real(YuadjYumu2YuadjYu(i2,i2),dp) 
 YuadjYuYumq2adjYu = Matmul(Yu,adjYuYumq2adjYu) 
 YuadjYuYuadjYumu2 = Matmul(Yu,adjYuYuadjYumu2) 
 YuadjYuYuadjYuYu = Matmul(Yu,adjYuYuadjYuYu) 
 YuadjYuYuadjYuTu = Matmul(Yu,adjYuYuadjYuTu) 
 YuadjYuTuadjYuYu = Matmul(Yu,adjYuTuadjYuYu) 
 YvadjYvYvadjYvvL = Matmul(Yv,adjYvYvadjYvvL) 
 YvadjYvYvadjYvYv = Matmul(Yv,adjYvYvadjYvYv) 
 YvadjYvYvadjYvTv = Matmul(Yv,adjYvYvadjYvTv) 
 YvadjYvCml2YvadjYv = Matmul(Yv,adjYvCml2YvadjYv) 
Forall(i2=1:3)  YvadjYvCml2YvadjYv(i2,i2) =  Real(YvadjYvCml2YvadjYv(i2,i2),dp) 
 YvadjYvCml2TpYeCYe = Matmul(Yv,adjYvCml2TpYeCYe) 
 YvadjYvTvadjYvYv = Matmul(Yv,adjYvTvadjYvYv) 
 YvadjYvTpYeCme2CYe = Matmul(Yv,adjYvTpYeCme2CYe) 
 YvadjYvTpYeCYeYv = Matmul(Yv,adjYvTpYeCYeYv) 
 YvadjYvTpYeCYeCml2 = Matmul(Yv,adjYvTpYeCYeCml2) 
 YvadjYvTpYeCYeTv = Matmul(Yv,adjYvTpYeCYeTv) 
 YvadjYvTpTeCYeYv = Matmul(Yv,adjYvTpTeCYeYv) 
 YvCmv2adjYvYvadjYv = Matmul(Yv,Cmv2adjYvYvadjYv) 
 YvCmv2adjYvTpYeCYe = Matmul(Yv,Cmv2adjYvTpYeCYe) 
 YvCmv2Ckap1kap1adjYv = Matmul(Yv,Cmv2Ckap1kap1adjYv) 
 YvCmv2Ckap2kap2adjYv = Matmul(Yv,Cmv2Ckap2kap2adjYv) 
 YvCmv2Ckap3kap3adjYv = Matmul(Yv,Cmv2Ckap3kap3adjYv) 
 YvCkap1TpYvCYvTpTk1 = Matmul(Yv,Ckap1TpYvCYvTpTk1) 
 YvCkap1Tpkap1adjYvvL = Matmul(Yv,Ckap1Tpkap1adjYvvL) 
 YvCkap1Tpkap1adjYvYv = Matmul(Yv,Ckap1Tpkap1adjYvYv) 
 YvCkap1Tpkap1Cmv2adjYv = Matmul(Yv,Ckap1Tpkap1Cmv2adjYv) 
 YvCkap1TpTk1adjYvYv = Matmul(Yv,Ckap1TpTk1adjYvYv) 
 YvCkap1kap1adjYvCml2 = Matmul(Yv,Ckap1kap1adjYvCml2) 
 YvCkap2TpYvCYvTpTk2 = Matmul(Yv,Ckap2TpYvCYvTpTk2) 
 YvCkap2Tpkap2adjYvvL = Matmul(Yv,Ckap2Tpkap2adjYvvL) 
 YvCkap2Tpkap2adjYvYv = Matmul(Yv,Ckap2Tpkap2adjYvYv) 
 YvCkap2Tpkap2Cmv2adjYv = Matmul(Yv,Ckap2Tpkap2Cmv2adjYv) 
 YvCkap2TpTk2adjYvYv = Matmul(Yv,Ckap2TpTk2adjYvYv) 
 YvCkap2kap2adjYvCml2 = Matmul(Yv,Ckap2kap2adjYvCml2) 
 YvCkap3TpYvCYvTpTk3 = Matmul(Yv,Ckap3TpYvCYvTpTk3) 
 YvCkap3Tpkap3adjYvvL = Matmul(Yv,Ckap3Tpkap3adjYvvL) 
 YvCkap3Tpkap3adjYvYv = Matmul(Yv,Ckap3Tpkap3adjYvYv) 
 YvCkap3Tpkap3Cmv2adjYv = Matmul(Yv,Ckap3Tpkap3Cmv2adjYv) 
 YvCkap3TpTk3adjYvYv = Matmul(Yv,Ckap3TpTk3adjYvYv) 
 YvCkap3kap3adjYvCml2 = Matmul(Yv,Ckap3kap3adjYvCml2) 
 adjYdmd2YdadjYdYd = Matmul(adjYd,md2YdadjYdYd) 
 adjYdYdmq2adjYdYd = Matmul(adjYd,Ydmq2adjYdYd) 
Forall(i2=1:3)  adjYdYdmq2adjYdYd(i2,i2) =  Real(adjYdYdmq2adjYdYd(i2,i2),dp) 
 adjYdYdadjYdmd2Yd = Matmul(adjYd,YdadjYdmd2Yd) 
 adjYdYdadjYdYdmq2 = Matmul(adjYd,YdadjYdYdmq2) 
 adjYeme2YeadjYeYe = Matmul(adjYe,me2YeadjYeYe) 
 adjYeme2YeCYvlam = Matmul(adjYe,me2YeCYvlam) 
 adjYeYeml2adjYeYe = Matmul(adjYe,Yeml2adjYeYe) 
Forall(i2=1:3)  adjYeYeml2adjYeYe(i2,i2) =  Real(adjYeYeml2adjYeYe(i2,i2),dp) 
 adjYeYeml2CYvlam = Matmul(adjYe,Yeml2CYvlam) 
 adjYeYeadjYeme2Ye = Matmul(adjYe,YeadjYeme2Ye) 
 adjYeYeadjYeYeml2 = Matmul(adjYe,YeadjYeYeml2) 
 adjYeYeCYvmv2lam = Matmul(adjYe,YeCYvmv2lam) 
 adjYumu2YuadjYuYu = Matmul(adjYu,mu2YuadjYuYu) 
 adjYuYumq2adjYuYu = Matmul(adjYu,Yumq2adjYuYu) 
Forall(i2=1:3)  adjYuYumq2adjYuYu(i2,i2) =  Real(adjYuYumq2adjYuYu(i2,i2),dp) 
 adjYuYuadjYumu2Yu = Matmul(adjYu,YuadjYumu2Yu) 
 adjYuYuadjYuYumq2 = Matmul(adjYu,YuadjYuYumq2) 
 adjYvYvadjkap1mv2kap2 = Matmul(adjYv,Yvadjkap1mv2kap2) 
 adjYvYvadjkap1mv2kap3 = Matmul(adjYv,Yvadjkap1mv2kap3) 
 adjYvYvadjkap2mv2kap1 = Matmul(adjYv,Yvadjkap2mv2kap1) 
 adjYvYvadjkap2mv2kap3 = Matmul(adjYv,Yvadjkap2mv2kap3) 
 adjYvYvadjkap3mv2kap1 = Matmul(adjYv,Yvadjkap3mv2kap1) 
 adjYvYvadjkap3mv2kap2 = Matmul(adjYv,Yvadjkap3mv2kap2) 
 adjYvYvCmv2Ckap1kap2 = Matmul(adjYv,YvCmv2Ckap1kap2) 
 adjYvYvCmv2Ckap1kap3 = Matmul(adjYv,YvCmv2Ckap1kap3) 
 adjYvYvCmv2Ckap2kap1 = Matmul(adjYv,YvCmv2Ckap2kap1) 
 adjYvYvCmv2Ckap2kap3 = Matmul(adjYv,YvCmv2Ckap2kap3) 
 adjYvYvCmv2Ckap3kap1 = Matmul(adjYv,YvCmv2Ckap3kap1) 
 adjYvYvCmv2Ckap3kap2 = Matmul(adjYv,YvCmv2Ckap3kap2) 
 adjYvCml2YvCkap1kap2 = Matmul(adjYv,Cml2YvCkap1kap2) 
 adjYvCml2YvCkap1kap3 = Matmul(adjYv,Cml2YvCkap1kap3) 
 adjYvCml2YvCkap2kap1 = Matmul(adjYv,Cml2YvCkap2kap1) 
 adjYvCml2YvCkap2kap3 = Matmul(adjYv,Cml2YvCkap2kap3) 
 adjYvCml2YvCkap3kap1 = Matmul(adjYv,Cml2YvCkap3kap1) 
 adjYvCml2YvCkap3kap2 = Matmul(adjYv,Cml2YvCkap3kap2) 
 CYvmv2TpYvCYvlam = Matmul(Conjg(Yv),mv2TpYvCYvlam) 
 CYvmv2TpYvCYvTpYv = Matmul(Conjg(Yv),mv2TpYvCYvTpYv) 
 CYvmv2kap1adjkap1TpYv = Matmul(Conjg(Yv),mv2kap1adjkap1TpYv) 
 CYvmv2kap2adjkap2TpYv = Matmul(Conjg(Yv),mv2kap2adjkap2TpYv) 
 CYvmv2kap3adjkap3TpYv = Matmul(Conjg(Yv),mv2kap3adjkap3TpYv) 
 CYvTpYvml2CYvlam = Matmul(Conjg(Yv),TpYvml2CYvlam) 
 CYvTpYvml2CYvTpYv = Matmul(Conjg(Yv),TpYvml2CYvTpYv) 
Forall(i2=1:3)  CYvTpYvml2CYvTpYv(i2,i2) =  Real(CYvTpYvml2CYvTpYv(i2,i2),dp) 
 CYvTpYvCYvmv2lam = Matmul(Conjg(Yv),TpYvCYvmv2lam) 
 CYvTpYvCYvmv2TpYv = Matmul(Conjg(Yv),TpYvCYvmv2TpYv) 
 CYvTpYvCYvTpYvml2 = Matmul(Conjg(Yv),TpYvCYvTpYvml2) 
 CYvkap1adjkap1TpYvml2 = Matmul(Conjg(Yv),kap1adjkap1TpYvml2) 
 CYvkap2adjkap2TpYvml2 = Matmul(Conjg(Yv),kap2adjkap2TpYvml2) 
 CYvkap3adjkap3TpYvml2 = Matmul(Conjg(Yv),kap3adjkap3TpYvml2) 
 TdadjYdYdadjYdYd = Matmul(Td,adjYdYdadjYdYd) 
 TdadjYuYuadjYdYd = Matmul(Td,adjYuYuadjYdYd) 
 TdadjYuYuadjYuYu = Matmul(Td,adjYuYuadjYuYu) 
 TeadjYeYeadjYeYe = Matmul(Te,adjYeYeadjYeYe) 
 TeCYvTpYvadjYeYe = Matmul(Te,CYvTpYvadjYeYe) 
 TeCYvTpYvCYvTpYv = Matmul(Te,CYvTpYvCYvTpYv) 
 TuadjYdYdadjYdYd = Matmul(Tu,adjYdYdadjYdYd) 
 TuadjYdYdadjYuYu = Matmul(Tu,adjYdYdadjYuYu) 
 TuadjYuYuadjYuYu = Matmul(Tu,adjYuYuadjYuYu) 
 TvadjYvYvadjYvYv = Matmul(Tv,adjYvYvadjYvYv) 
 TvadjYvTpYeCYeYv = Matmul(Tv,adjYvTpYeCYeYv) 
 TpYeCYeTpYeCYevL = Matmul(Transpose(Ye),CYeTpYeCYevL) 
 TpYeCYeTpYeCYeYv = Matmul(Transpose(Ye),CYeTpYeCYeYv) 
 TpYeCYeTpYeCYeTv = Matmul(Transpose(Ye),CYeTpYeCYeTv) 
 TpYeCYeTpTeCYeYv = Matmul(Transpose(Ye),CYeTpTeCYeYv) 
 TpYvml2adjYeYeCYv = Matmul(Transpose(Yv),ml2adjYeYeCYv) 
 TpYvml2CYvTpYvCYv = Matmul(Transpose(Yv),ml2CYvTpYvCYv) 
 TpYvadjYeme2YeCYv = Matmul(Transpose(Yv),adjYeme2YeCYv) 
Forall(i2=1:3)  TpYvadjYeme2YeCYv(i2,i2) =  Real(TpYvadjYeme2YeCYv(i2,i2),dp) 
 TpYvadjYeYeml2CYv = Matmul(Transpose(Yv),adjYeYeml2CYv) 
 TpYvadjYeYeCYvmv2 = Matmul(Transpose(Yv),adjYeYeCYvmv2) 
 TpYvadjYeYeCYvvR = Matmul(Transpose(Yv),adjYeYeCYvvR) 
 TpYvadjYeYeCYvlam = Matmul(Transpose(Yv),adjYeYeCYvlam) 
 TpYvadjYeTeCYvlam = Matmul(Transpose(Yv),adjYeTeCYvlam) 
 TpYvCYvmv2TpYvCYv = Matmul(Transpose(Yv),CYvmv2TpYvCYv) 
Forall(i2=1:3)  TpYvCYvmv2TpYvCYv(i2,i2) =  Real(TpYvCYvmv2TpYvCYv(i2,i2),dp) 
 TpYvCYvTpYvml2CYv = Matmul(Transpose(Yv),CYvTpYvml2CYv) 
 TpYvCYvTpYvCYvmv2 = Matmul(Transpose(Yv),CYvTpYvCYvmv2) 
 TpYvCYvTpYvCYvvR = Matmul(Transpose(Yv),CYvTpYvCYvvR) 
 TpYvCYvTpYvCYvlam = Matmul(Transpose(Yv),CYvTpYvCYvlam) 
 TpYvCYvTpYvCYvTlam = Matmul(Transpose(Yv),CYvTpYvCYvTlam) 
 TpYvCYvTpTvCYvlam = Matmul(Transpose(Yv),CYvTpTvCYvlam) 
 TpYvCYvkap1adjkap1lam = Matmul(Transpose(Yv),CYvkap1adjkap1lam) 
 TpYvCYvkap2adjkap2lam = Matmul(Transpose(Yv),CYvkap2adjkap2lam) 
 TpYvCYvkap3adjkap3lam = Matmul(Transpose(Yv),CYvkap3adjkap3lam) 
 TpYvCYvTk1adjkap1lam = Matmul(Transpose(Yv),CYvTk1adjkap1lam) 
 TpYvCYvTk2adjkap2lam = Matmul(Transpose(Yv),CYvTk2adjkap2lam) 
 TpYvCYvTk3adjkap3lam = Matmul(Transpose(Yv),CYvTk3adjkap3lam) 
 TpTeCYeTpYeCYeYv = Matmul(Transpose(Te),CYeTpYeCYeYv) 
 TpTvadjYeYeCYvlam = Matmul(Transpose(Tv),adjYeYeCYvlam) 
 TpTvCYvTpYvCYvlam = Matmul(Transpose(Tv),CYvTpYvCYvlam) 
 kap1adjYvYvadjkap1mv2 = Matmul(kap1,adjYvYvadjkap1mv2) 
 kap1adjYvYvCmv2adjkap1 = Matmul(kap1,adjYvYvCmv2adjkap1) 
 kap1Cmv2adjYvYvadjkap1 = Matmul(kap1,Cmv2adjYvYvadjkap1) 
 kap1Ckap1Tpkap1Cmv2adjkap1 = Matmul(kap1,Ckap1Tpkap1Cmv2adjkap1) 
 kap1Ckap2Tpkap2Cmv2adjkap1 = Matmul(kap1,Ckap2Tpkap2Cmv2adjkap1) 
 kap1Ckap3Tpkap3Cmv2adjkap1 = Matmul(kap1,Ckap3Tpkap3Cmv2adjkap1) 
 kap2adjYvYvadjkap2mv2 = Matmul(kap2,adjYvYvadjkap2mv2) 
 kap2adjYvYvCmv2adjkap2 = Matmul(kap2,adjYvYvCmv2adjkap2) 
 kap2Cmv2adjYvYvadjkap2 = Matmul(kap2,Cmv2adjYvYvadjkap2) 
 kap2Ckap1Tpkap1Cmv2adjkap2 = Matmul(kap2,Ckap1Tpkap1Cmv2adjkap2) 
 kap2Ckap2Tpkap2Cmv2adjkap2 = Matmul(kap2,Ckap2Tpkap2Cmv2adjkap2) 
 kap2Ckap3Tpkap3Cmv2adjkap2 = Matmul(kap2,Ckap3Tpkap3Cmv2adjkap2) 
 kap3adjYvYvadjkap3mv2 = Matmul(kap3,adjYvYvadjkap3mv2) 
 kap3adjYvYvCmv2adjkap3 = Matmul(kap3,adjYvYvCmv2adjkap3) 
 kap3Cmv2adjYvYvadjkap3 = Matmul(kap3,Cmv2adjYvYvadjkap3) 
 kap3Ckap1Tpkap1Cmv2adjkap3 = Matmul(kap3,Ckap1Tpkap1Cmv2adjkap3) 
 kap3Ckap2Tpkap2Cmv2adjkap3 = Matmul(kap3,Ckap2Tpkap2Cmv2adjkap3) 
 kap3Ckap3Tpkap3Cmv2adjkap3 = Matmul(kap3,Ckap3Tpkap3Cmv2adjkap3) 
 Tk1adjYvYvCkap1lam = Matmul(Tk1,adjYvYvCkap1lam) 
 Tk2adjYvYvCkap2lam = Matmul(Tk2,adjYvYvCkap2lam) 
 Tk3adjYvYvCkap3lam = Matmul(Tk3,adjYvYvCkap3lam) 
 TradjTk1Tk1 = cTrace(adjTk1Tk1) 
 TradjTk2Tk2 = cTrace(adjTk2Tk2) 
 TradjTk3Tk3 = cTrace(adjTk3Tk3) 
 TrCTdTpYd = cTrace(CTdTpYd) 
 TrCTeTpYe = cTrace(CTeTpYe) 
 TrCTuTpYu = cTrace(CTuTpYu) 
 TrCTvTpYv = cTrace(CTvTpYv) 
 TrCkap1Tpkap1 = cTrace(Ckap1Tpkap1) 
 TrCkap1Tpkap2 = cTrace(Ckap1Tpkap2) 
 TrCkap1Tpkap3 = cTrace(Ckap1Tpkap3) 
 TrCkap1kap1 = cTrace(Ckap1kap1) 
 TrCkap1kap2 = cTrace(Ckap1kap2) 
 TrCkap1kap3 = cTrace(Ckap1kap3) 
 TrCkap1Tk1 = cTrace(Ckap1Tk1) 
 TrCkap1Tk2 = cTrace(Ckap1Tk2) 
 TrCkap1Tk3 = cTrace(Ckap1Tk3) 
 TrCkap2Tpkap1 = cTrace(Ckap2Tpkap1) 
 TrCkap2Tpkap2 = cTrace(Ckap2Tpkap2) 
 TrCkap2Tpkap3 = cTrace(Ckap2Tpkap3) 
 TrCkap2kap1 = cTrace(Ckap2kap1) 
 TrCkap2kap2 = cTrace(Ckap2kap2) 
 TrCkap2kap3 = cTrace(Ckap2kap3) 
 TrCkap2Tk1 = cTrace(Ckap2Tk1) 
 TrCkap2Tk2 = cTrace(Ckap2Tk2) 
 TrCkap2Tk3 = cTrace(Ckap2Tk3) 
 TrCkap3Tpkap1 = cTrace(Ckap3Tpkap1) 
 TrCkap3Tpkap2 = cTrace(Ckap3Tpkap2) 
 TrCkap3Tpkap3 = cTrace(Ckap3Tpkap3) 
 TrCkap3kap1 = cTrace(Ckap3kap1) 
 TrCkap3kap2 = cTrace(Ckap3kap2) 
 TrCkap3kap3 = cTrace(Ckap3kap3) 
 TrCkap3Tk1 = cTrace(Ckap3Tk1) 
 TrCkap3Tk2 = cTrace(Ckap3Tk2) 
 TrCkap3Tk3 = cTrace(Ckap3Tk3) 
 TrCTk1Tpkap1 = cTrace(CTk1Tpkap1) 
 TrCTk1Tpkap2 = cTrace(CTk1Tpkap2) 
 TrCTk1Tpkap3 = cTrace(CTk1Tpkap3) 
 TrCTk2Tpkap1 = cTrace(CTk2Tpkap1) 
 TrCTk2Tpkap2 = cTrace(CTk2Tpkap2) 
 TrCTk2Tpkap3 = cTrace(CTk2Tpkap3) 
 TrCTk3Tpkap1 = cTrace(CTk3Tpkap1) 
 TrCTk3Tpkap2 = cTrace(CTk3Tpkap2) 
 TrCTk3Tpkap3 = cTrace(CTk3Tpkap3) 
 Trml2YvadjYv = cTrace(ml2YvadjYv) 
 Trmv2kap1adjkap1 = cTrace(mv2kap1adjkap1) 
 Trmv2kap1adjkap2 = cTrace(mv2kap1adjkap2) 
 Trmv2kap1adjkap3 = cTrace(mv2kap1adjkap3) 
 Trmv2kap2adjkap1 = cTrace(mv2kap2adjkap1) 
 Trmv2kap2adjkap2 = cTrace(mv2kap2adjkap2) 
 Trmv2kap2adjkap3 = cTrace(mv2kap2adjkap3) 
 Trmv2kap3adjkap1 = cTrace(mv2kap3adjkap1) 
 Trmv2kap3adjkap2 = cTrace(mv2kap3adjkap2) 
 Trmv2kap3adjkap3 = cTrace(mv2kap3adjkap3) 
 TrYdadjYdCmd2 = cTrace(YdadjYdCmd2) 
 TrYdCmq2adjYd = cTrace(YdCmq2adjYd) 
 TrYeadjYeCme2 = cTrace(YeadjYeCme2) 
 TrYeCml2adjYe = cTrace(YeCml2adjYe) 
 TrYuadjYuCmu2 = cTrace(YuadjYuCmu2) 
 TrYuCmq2adjYu = cTrace(YuCmq2adjYu) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdTdadjYd = cTrace(YdadjYdTdadjYd) 
 TrYdadjYdTdadjTd = cTrace(YdadjYdTdadjTd) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYdadjYuTuadjYd = cTrace(YdadjYuTuadjYd) 
 TrYdadjYuTuadjTd = cTrace(YdadjYuTuadjTd) 
 TrYdadjTdTdadjYd = cTrace(YdadjTdTdadjYd) 
 TrYdadjTuTuadjYd = cTrace(YdadjTuTuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYeTeadjYe = cTrace(YeadjYeTeadjYe) 
 TrYeadjYeTeadjTe = cTrace(YeadjYeTeadjTe) 
 TrYeadjTeTeadjYe = cTrace(YeadjTeTeadjYe) 
 TrYeCTvTpTvadjYe = cTrace(YeCTvTpTvadjYe) 
 TrYuadjYdTdadjYu = cTrace(YuadjYdTdadjYu) 
 TrYuadjYdTdadjTu = cTrace(YuadjYdTdadjTu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYuadjYuTuadjYu = cTrace(YuadjYuTuadjYu) 
 TrYuadjYuTuadjTu = cTrace(YuadjYuTuadjTu) 
 TrYuadjTdTdadjYu = cTrace(YuadjTdTdadjYu) 
 TrYuadjTuTuadjYu = cTrace(YuadjTuTuadjYu) 
 TrYvadjYvYvadjYv = cTrace(YvadjYvYvadjYv) 
 TrYvadjYvTvadjYv = cTrace(YvadjYvTvadjYv) 
 TrYvadjYvTvadjTv = cTrace(YvadjYvTvadjTv) 
 TrYvadjYvTpYeCYe = cTrace(YvadjYvTpYeCYe) 
 TrYvadjYvTpTeCTe = cTrace(YvadjYvTpTeCTe) 
 TrYvadjTvTvadjYv = cTrace(YvadjTvTvadjYv) 
 TrYvCkap1Tpkap1adjYv = cTrace(YvCkap1Tpkap1adjYv) 
 TrYvCkap1TpTk1adjYv = cTrace(YvCkap1TpTk1adjYv) 
 TrYvCkap1kap1adjYv = cTrace(YvCkap1kap1adjYv) 
 TrYvCkap1Tk1adjTv = cTrace(YvCkap1Tk1adjTv) 
 TrYvCkap2Tpkap2adjYv = cTrace(YvCkap2Tpkap2adjYv) 
 TrYvCkap2TpTk2adjYv = cTrace(YvCkap2TpTk2adjYv) 
 TrYvCkap2kap2adjYv = cTrace(YvCkap2kap2adjYv) 
 TrYvCkap2Tk2adjTv = cTrace(YvCkap2Tk2adjTv) 
 TrYvCkap3Tpkap3adjYv = cTrace(YvCkap3Tpkap3adjYv) 
 TrYvCkap3TpTk3adjYv = cTrace(YvCkap3TpTk3adjYv) 
 TrYvCkap3kap3adjYv = cTrace(YvCkap3kap3adjYv) 
 TrYvCkap3Tk3adjTv = cTrace(YvCkap3Tk3adjTv) 
 TrYvCTk1TpTk1adjYv = cTrace(YvCTk1TpTk1adjYv) 
 TrYvCTk2TpTk2adjYv = cTrace(YvCTk2TpTk2adjYv) 
 TrYvCTk3TpTk3adjYv = cTrace(YvCTk3TpTk3adjYv) 
 TradjYeTeCYvTpYv = cTrace(adjYeTeCYvTpYv) 
 TradjYeTeCTvTpYv = cTrace(adjYeTeCTvTpYv) 
 TradjYvTvadjkap1kap1 = cTrace(adjYvTvadjkap1kap1) 
 TradjYvTvadjkap2kap2 = cTrace(adjYvTvadjkap2kap2) 
 TradjYvTvadjkap3kap3 = cTrace(adjYvTvadjkap3kap3) 
 TradjYvTpYeCYeTv = cTrace(adjYvTpYeCYeTv) 
 TradjYvTpYeCTeTv = cTrace(adjYvTpYeCTeTv) 
 Trmd2YdadjYdYdadjYd = cTrace(md2YdadjYdYdadjYd) 
 Trmd2YdadjYuYuadjYd = cTrace(md2YdadjYuYuadjYd) 
 Trme2YeadjYeYeadjYe = cTrace(me2YeadjYeYeadjYe) 
 Trml2adjYeYeadjYeYe = cTrace(ml2adjYeYeadjYeYe) 
 Trmq2adjYdYdadjYdYd = cTrace(mq2adjYdYdadjYdYd) 
 Trmq2adjYdYdadjYuYu = cTrace(mq2adjYdYdadjYuYu) 
 Trmq2adjYuYuadjYdYd = cTrace(mq2adjYuYuadjYdYd) 
 Trmq2adjYuYuadjYuYu = cTrace(mq2adjYuYuadjYuYu) 
 Trmu2YuadjYdYdadjYu = cTrace(mu2YuadjYdYdadjYu) 
 Trmu2YuadjYuYuadjYu = cTrace(mu2YuadjYuYuadjYu) 
 Trmv2kap1adjYvYvadjkap1 = cTrace(mv2kap1adjYvYvadjkap1) 
 Trmv2kap2adjYvYvadjkap2 = cTrace(mv2kap2adjYvYvadjkap2) 
 Trmv2kap3adjYvYvadjkap3 = cTrace(mv2kap3adjYvYvadjkap3) 
 TrYvadjYvCml2YvadjYv = cTrace(YvadjYvCml2YvadjYv) 
 TrYvadjYvCml2TpYeCYe = cTrace(YvadjYvCml2TpYeCYe) 
 TrYvadjYvTpYeCme2CYe = cTrace(YvadjYvTpYeCme2CYe) 
 TrYvadjYvTpYeCYeCml2 = cTrace(YvadjYvTpYeCYeCml2) 
 TrYvCmv2adjYvYvadjYv = cTrace(YvCmv2adjYvYvadjYv) 
 TrYvCmv2adjYvTpYeCYe = cTrace(YvCmv2adjYvTpYeCYe) 
 TrYvCmv2Ckap1kap1adjYv = cTrace(YvCmv2Ckap1kap1adjYv) 
 TrYvCmv2Ckap2kap2adjYv = cTrace(YvCmv2Ckap2kap2adjYv) 
 TrYvCmv2Ckap3kap3adjYv = cTrace(YvCmv2Ckap3kap3adjYv) 
 TrYvCkap1Tpkap1Cmv2adjYv = cTrace(YvCkap1Tpkap1Cmv2adjYv) 
 TrYvCkap1kap1adjYvCml2 = cTrace(YvCkap1kap1adjYvCml2) 
 TrYvCkap2Tpkap2Cmv2adjYv = cTrace(YvCkap2Tpkap2Cmv2adjYv) 
 TrYvCkap2kap2adjYvCml2 = cTrace(YvCkap2kap2adjYvCml2) 
 TrYvCkap3Tpkap3Cmv2adjYv = cTrace(YvCkap3Tpkap3Cmv2adjYv) 
 TrYvCkap3kap3adjYvCml2 = cTrace(YvCkap3kap3adjYvCml2) 
 SPmlHd2xxadjYeYeCYvlam = DOT_PRODUCT(mlHd2,adjYeYeCYvlam) 
 SPvLxxadjYeYeCYvlam = DOT_PRODUCT(vL,adjYeYeCYvlam) 
 SPvRxxadjYvYvClam = DOT_PRODUCT(vR,adjYvYvClam) 
 SPlamxxCTlam = DOT_PRODUCT(Conjg(lam),Conjg(Tlam)) 
 SPlamxxadjkap1lam = DOT_PRODUCT(Conjg(lam),adjkap1lam) 
 SPlamxxadjkap2lam = DOT_PRODUCT(Conjg(lam),adjkap2lam) 
 SPlamxxadjkap3lam = DOT_PRODUCT(Conjg(lam),adjkap3lam) 
 SPlamxxadjYvYvClam = DOT_PRODUCT(Conjg(lam),adjYvYvClam) 
 SPlamxxadjTvTvClam = DOT_PRODUCT(Conjg(lam),adjTvTvClam) 
 SPlamxxadjYvYvadjYvmlHd2 = DOT_PRODUCT(Conjg(lam),adjYvYvadjYvmlHd2) 
 SPlamxxadjYvYvadjYvvL = DOT_PRODUCT(Conjg(lam),adjYvYvadjYvvL) 
 SPlamxxadjYvYvCmv2Clam = DOT_PRODUCT(Conjg(lam),adjYvYvCmv2Clam) 
 SPlamxxCmv2adjYvYvClam = DOT_PRODUCT(Conjg(lam),Cmv2adjYvYvClam) 
 SPClamxxTpYvmlHd2 = DOT_PRODUCT(Conjg(Conjg(lam)),TpYvmlHd2) 
 SPClamxxTpTvCYvlam = DOT_PRODUCT(Conjg(Conjg(lam)),TpTvCYvlam) 
 SPClamxxTpYvml2CYvlam = DOT_PRODUCT(Conjg(Conjg(lam)),TpYvml2CYvlam) 
 SPadjkap1lamxxTpkap1adjYvvL = DOT_PRODUCT(Conjg(adjkap1lam),Tpkap1adjYvvL) 
 SPadjkap2lamxxTpkap2adjYvvL = DOT_PRODUCT(Conjg(adjkap2lam),Tpkap2adjYvvL) 
 SPadjkap3lamxxTpkap3adjYvvL = DOT_PRODUCT(Conjg(adjkap3lam),Tpkap3adjYvvL) 
 SPTpkap1Clamxxadjkap1lam = DOT_PRODUCT(Conjg(Tpkap1Clam),adjkap1lam) 
 SPTpkap2Clamxxadjkap2lam = DOT_PRODUCT(Conjg(Tpkap2Clam),adjkap2lam) 
 SPTpkap3Clamxxadjkap3lam = DOT_PRODUCT(Conjg(Tpkap3Clam),adjkap3lam) 
 SPTpTk1ClamxxadjTk1lam = DOT_PRODUCT(Conjg(TpTk1Clam),adjTk1lam) 
 SPTpTk2ClamxxadjTk2lam = DOT_PRODUCT(Conjg(TpTk2Clam),adjTk2lam) 
 SPTpTk3ClamxxadjTk3lam = DOT_PRODUCT(Conjg(TpTk3Clam),adjTk3lam) 
 SPadjYvYvClamxxTlam = DOT_PRODUCT(Conjg(adjYvYvClam),Tlam) 
 SPadjYvYvCTlamxxTlam = DOT_PRODUCT(Conjg(adjYvYvCTlam),Tlam) 
 SPadjTvYvClamxxTlam = DOT_PRODUCT(Conjg(adjTvYvClam),Tlam) 
 SPTpTvCYvlamxxCTlam = DOT_PRODUCT(Conjg(TpTvCYvlam),Conjg(Tlam)) 
 g1p4 =g1**4 
 g1p5 =g1**5 
 g2p4 =g2**4 
 g2p5 =g2**5 
 g3p4 =g3**4 
 g3p5 =g3**5 
 Xip2 =Xi**2 
 SPlamxxClamp2 =SPlamxxClam**2 
Do i1=1,3
  Do i2=1,3
Dylami1CTlami2(i1,i2) = Conjg(Tlam(i2))*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1YvClami2(i1,i2) = CYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTlami1YvCTlami2(i1,i2) = CYvTlam(i1)*YvCTlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvTlami1YvClami2(i1,i2) = CTvTlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvlami1TvClami2(i1,i2) = CTvlam(i1)*TvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1TvCTlami2(i1,i2) = CYvlam(i1)*TvCTlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Clami2(i1,i2) = kap1Clam(i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Clami2(i1,i2) = kap1Clam(i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Clami2(i1,i2) = kap1Clam(i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Clami2(i1,i2) = kap2Clam(i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Clami2(i1,i2) = kap2Clam(i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Clami2(i1,i2) = kap2Clam(i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Clami2(i1,i2) = kap3Clam(i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Clami2(i1,i2) = kap3Clam(i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Clami2(i1,i2) = kap3Clam(i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1Clami1CTk1lami2(i1,i2) = CTk1lam(i2)*Tk1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2Clami1CTk2lami2(i1,i2) = CTk2lam(i2)*Tk2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3Clami1CTk3lami2(i1,i2) = CTk3lam(i2)*Tk3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvlami1YvClami2(i1,i2) = ml2CYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1Clami1Ckap1lami2(i1,i2) = Ckap1lam(i2)*mv2kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2Clami1Ckap2lami2(i1,i2) = Ckap2lam(i2)*mv2kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3Clami1Ckap3lami2(i1,i2) = Ckap3lam(i2)*mv2kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvlami1YvClami2(i1,i2) = YeCYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvlami1TvClami2(i1,i2) = TvClam(i2)*YeCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvTlami1YvClami2(i1,i2) = YeCYvTlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1YvadjYvmlHd2i2(i1,i2) = CYvlam(i1)*YvadjYvmlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap1i11adjYvTv1i2(i1,i2) = adjYvTv(1,i2)*Yvadjkap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap1kap1i21(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi11Yvadjkap1kap1i21(i1,i2) = ml2CYv(i1,1)*Yvadjkap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi11Yvadjkap1kap1i21(i1,i2) = TeCYv(i1,1)*Yvadjkap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap2i11adjYvTv2i2(i1,i2) = adjYvTv(2,i2)*Yvadjkap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap1kap2i21(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi12Yvadjkap1kap2i21(i1,i2) = ml2CYv(i1,2)*Yvadjkap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi12Yvadjkap1kap2i21(i1,i2) = TeCYv(i1,2)*Yvadjkap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap3i11adjYvTv3i2(i1,i2) = adjYvTv(3,i2)*Yvadjkap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap1kap3i21(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi13Yvadjkap1kap3i21(i1,i2) = ml2CYv(i1,3)*Yvadjkap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi13Yvadjkap1kap3i21(i1,i2) = TeCYv(i1,3)*Yvadjkap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1Tk1i11lami2(i1,i2) = Yvadjkap1Tk1(i1,1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1Tk2i11lami2(i1,i2) = Yvadjkap1Tk2(i1,1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1Tk3i11lami2(i1,i2) = Yvadjkap1Tk3(i1,1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap1i12adjYvTv1i2(i1,i2) = adjYvTv(1,i2)*Yvadjkap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap2kap1i22(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi11Yvadjkap2kap1i22(i1,i2) = ml2CYv(i1,1)*Yvadjkap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi11Yvadjkap2kap1i22(i1,i2) = TeCYv(i1,1)*Yvadjkap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap2i12adjYvTv2i2(i1,i2) = adjYvTv(2,i2)*Yvadjkap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap2kap2i22(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi12Yvadjkap2kap2i22(i1,i2) = ml2CYv(i1,2)*Yvadjkap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi12Yvadjkap2kap2i22(i1,i2) = TeCYv(i1,2)*Yvadjkap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap3i12adjYvTv3i2(i1,i2) = adjYvTv(3,i2)*Yvadjkap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap2kap3i22(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi13Yvadjkap2kap3i22(i1,i2) = ml2CYv(i1,3)*Yvadjkap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi13Yvadjkap2kap3i22(i1,i2) = TeCYv(i1,3)*Yvadjkap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2Tk1i12lami2(i1,i2) = Yvadjkap2Tk1(i1,2)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2Tk2i12lami2(i1,i2) = Yvadjkap2Tk2(i1,2)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2Tk3i12lami2(i1,i2) = Yvadjkap2Tk3(i1,2)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap1i13adjYvTv1i2(i1,i2) = adjYvTv(1,i2)*Yvadjkap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap3kap1i23(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi11Yvadjkap3kap1i23(i1,i2) = ml2CYv(i1,1)*Yvadjkap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi11Yvadjkap3kap1i23(i1,i2) = TeCYv(i1,1)*Yvadjkap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap2i13adjYvTv2i2(i1,i2) = adjYvTv(2,i2)*Yvadjkap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap3kap2i23(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi12Yvadjkap3kap2i23(i1,i2) = ml2CYv(i1,2)*Yvadjkap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi12Yvadjkap3kap2i23(i1,i2) = TeCYv(i1,2)*Yvadjkap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap3i13adjYvTv3i2(i1,i2) = adjYvTv(3,i2)*Yvadjkap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap3kap3i23(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dyml2CYvi13Yvadjkap3kap3i23(i1,i2) = ml2CYv(i1,3)*Yvadjkap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvi13Yvadjkap3kap3i23(i1,i2) = TeCYv(i1,3)*Yvadjkap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3Tk1i13lami2(i1,i2) = Yvadjkap3Tk1(i1,3)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3Tk2i13lami2(i1,i2) = Yvadjkap3Tk2(i1,3)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3Tk3i13lami2(i1,i2) = Yvadjkap3Tk3(i1,3)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1YvCmv2Clami2(i1,i2) = CYvlam(i1)*YvCmv2Clam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1lami1Tk1Clami2(i1,i2) = Tk1Clam(i2)*YvCkap1lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11YvCkap1Tk1i21(i1,i2) = Conjg(Tv(i1,1))*YvCkap1Tk1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12YvCkap1Tk2i21(i1,i2) = Conjg(Tv(i1,2))*YvCkap1Tk2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13YvCkap1Tk3i21(i1,i2) = Conjg(Tv(i1,3))*YvCkap1Tk3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2lami1Tk2Clami2(i1,i2) = Tk2Clam(i2)*YvCkap2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11YvCkap2Tk1i22(i1,i2) = Conjg(Tv(i1,1))*YvCkap2Tk1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12YvCkap2Tk2i22(i1,i2) = Conjg(Tv(i1,2))*YvCkap2Tk2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13YvCkap2Tk3i22(i1,i2) = Conjg(Tv(i1,3))*YvCkap2Tk3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3lami1Tk3Clami2(i1,i2) = Tk3Clam(i2)*YvCkap3lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11YvCkap3Tk1i23(i1,i2) = Conjg(Tv(i1,1))*YvCkap3Tk1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12YvCkap3Tk2i23(i1,i2) = Conjg(Tv(i1,2))*YvCkap3Tk2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13YvCkap3Tk3i23(i1,i2) = Conjg(Tv(i1,3))*YvCkap3Tk3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvYvClami2(i1,i2) = adjYvYvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2lami1adjYvYvClami2(i1,i2) = adjYvYvClam(i2)*mv2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvTvCTlami2(i1,i2) = adjYvTvCTlam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjTvTvClami2(i1,i2) = adjTvTvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvlami1Cml2YvClami2(i1,i2) = Cml2YvClam(i2)*CYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvmv2lami1YvClami2(i1,i2) = CYvmv2lam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap1Ckap1i11mlHd2i2(i1,i2) = CYvTpkap1Ckap1(i1,1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap1Ckap2i11mlHd2i2(i1,i2) = CYvTpkap1Ckap2(i1,1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap1Ckap3i11mlHd2i2(i1,i2) = CYvTpkap1Ckap3(i1,1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap2Ckap1i12mlHd2i2(i1,i2) = CYvTpkap2Ckap1(i1,2)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap2Ckap2i12mlHd2i2(i1,i2) = CYvTpkap2Ckap2(i1,2)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap2Ckap3i12mlHd2i2(i1,i2) = CYvTpkap2Ckap3(i1,2)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap3Ckap1i13mlHd2i2(i1,i2) = CYvTpkap3Ckap1(i1,3)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap3Ckap2i13mlHd2i2(i1,i2) = CYvTpkap3Ckap2(i1,3)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpkap3Ckap3i13mlHd2i2(i1,i2) = CYvTpkap3Ckap3(i1,3)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Clami1Ckap1mv2lami2(i1,i2) = Ckap1mv2lam(i2)*kap1Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i11Ckap1kap1Ckap1i21(i1,i2) = Ckap1kap1Ckap1(i2,1)*mv2kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i12Ckap1kap1Ckap2i21(i1,i2) = Ckap1kap1Ckap2(i2,1)*mv2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i13Ckap1kap1Ckap3i21(i1,i2) = Ckap1kap1Ckap3(i2,1)*mv2kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i11Ckap1kap2Ckap1i22(i1,i2) = Ckap1kap2Ckap1(i2,2)*mv2kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i12Ckap1kap2Ckap2i22(i1,i2) = Ckap1kap2Ckap2(i2,2)*mv2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i13Ckap1kap2Ckap3i22(i1,i2) = Ckap1kap2Ckap3(i2,2)*mv2kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i11Ckap1kap3Ckap1i23(i1,i2) = Ckap1kap3Ckap1(i2,3)*mv2kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i12Ckap1kap3Ckap2i23(i1,i2) = Ckap1kap3Ckap2(i2,3)*mv2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap1i13Ckap1kap3Ckap3i23(i1,i2) = Ckap1kap3Ckap3(i2,3)*mv2kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Clami1Ckap2mv2lami2(i1,i2) = Ckap2mv2lam(i2)*kap2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i11Ckap2kap1Ckap1i21(i1,i2) = Ckap2kap1Ckap1(i2,1)*mv2kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i12Ckap2kap1Ckap2i21(i1,i2) = Ckap2kap1Ckap2(i2,1)*mv2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i13Ckap2kap1Ckap3i21(i1,i2) = Ckap2kap1Ckap3(i2,1)*mv2kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i11Ckap2kap2Ckap1i22(i1,i2) = Ckap2kap2Ckap1(i2,2)*mv2kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i12Ckap2kap2Ckap2i22(i1,i2) = Ckap2kap2Ckap2(i2,2)*mv2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i13Ckap2kap2Ckap3i22(i1,i2) = Ckap2kap2Ckap3(i2,2)*mv2kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i11Ckap2kap3Ckap1i23(i1,i2) = Ckap2kap3Ckap1(i2,3)*mv2kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i12Ckap2kap3Ckap2i23(i1,i2) = Ckap2kap3Ckap2(i2,3)*mv2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap2i13Ckap2kap3Ckap3i23(i1,i2) = Ckap2kap3Ckap3(i2,3)*mv2kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Clami1Ckap3mv2lami2(i1,i2) = Ckap3mv2lam(i2)*kap3Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i11Ckap3kap1Ckap1i21(i1,i2) = Ckap3kap1Ckap1(i2,1)*mv2kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i12Ckap3kap1Ckap2i21(i1,i2) = Ckap3kap1Ckap2(i2,1)*mv2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i13Ckap3kap1Ckap3i21(i1,i2) = Ckap3kap1Ckap3(i2,1)*mv2kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i11Ckap3kap2Ckap1i22(i1,i2) = Ckap3kap2Ckap1(i2,2)*mv2kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i12Ckap3kap2Ckap2i22(i1,i2) = Ckap3kap2Ckap2(i2,2)*mv2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i13Ckap3kap2Ckap3i22(i1,i2) = Ckap3kap2Ckap3(i2,2)*mv2kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i11Ckap3kap3Ckap1i23(i1,i2) = Ckap3kap3Ckap1(i2,3)*mv2kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i12Ckap3kap3Ckap2i23(i1,i2) = Ckap3kap3Ckap2(i2,3)*mv2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2kap3i13Ckap3kap3Ckap3i23(i1,i2) = Ckap3kap3Ckap3(i2,3)*mv2kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTeCYvlami1YvClami2(i1,i2) = TeCYvlam(i1)*YvClam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap1kap1i11TpYvCYvi21(i1,i2) = TpYvCYv(i2,1)*TvCkap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi11TvCkap1kap1i21(i1,i2) = TvCkap1kap1(i2,1)*YeCYv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap1kap2i11TpYvCYvi22(i1,i2) = TpYvCYv(i2,2)*TvCkap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi12TvCkap1kap2i21(i1,i2) = TvCkap1kap2(i2,1)*YeCYv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap1kap3i11TpYvCYvi23(i1,i2) = TpYvCYv(i2,3)*TvCkap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi13TvCkap1kap3i21(i1,i2) = TvCkap1kap3(i2,1)*YeCYv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap2kap1i12TpYvCYvi21(i1,i2) = TpYvCYv(i2,1)*TvCkap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi11TvCkap2kap1i22(i1,i2) = TvCkap2kap1(i2,2)*YeCYv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap2kap2i12TpYvCYvi22(i1,i2) = TpYvCYv(i2,2)*TvCkap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi12TvCkap2kap2i22(i1,i2) = TvCkap2kap2(i2,2)*YeCYv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap2kap3i12TpYvCYvi23(i1,i2) = TpYvCYv(i2,3)*TvCkap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi13TvCkap2kap3i22(i1,i2) = TvCkap2kap3(i2,2)*YeCYv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap3kap1i13TpYvCYvi21(i1,i2) = TpYvCYv(i2,1)*TvCkap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi11TvCkap3kap1i23(i1,i2) = TvCkap3kap1(i2,3)*YeCYv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap3kap2i13TpYvCYvi22(i1,i2) = TpYvCYv(i2,2)*TvCkap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi12TvCkap3kap2i23(i1,i2) = TvCkap3kap2(i2,3)*YeCYv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvCkap3kap3i13TpYvCYvi23(i1,i2) = TpYvCYv(i2,3)*TvCkap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYeCYvi13TvCkap3kap3i23(i1,i2) = TvCkap3kap3(i2,3)*YeCYv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvlami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvlami1adjYvmlHd2i2(i1,i2) = adjYvmlHd2(i2)*TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvlami1Cmv2Clami2(i1,i2) = Cmv2Clam(i2)*TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1TpYvCYvlami2(i1,i2) = TpYvCYvlam(i2)*YvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvClami1TpYvCYvlami2(i1,i2) = TpYvCYvlam(i2)*TvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTlami1CTlami2(i1,i2) = Conjg(Tlam(i2))*TpYvCYvTlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1TpYvCYvTlami2(i1,i2) = TpYvCYvTlam(i2)*YvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap11i1Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*TpYvCYvTpkap1(1,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11TpYvCYvTpkap11i2(i1,i2) = TpYvCYvTpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap12i1Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*TpYvCYvTpkap1(2,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11TpYvCYvTpkap12i2(i1,i2) = TpYvCYvTpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap13i1Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*TpYvCYvTpkap1(3,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11TpYvCYvTpkap13i2(i1,i2) = TpYvCYvTpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap21i1Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*TpYvCYvTpkap2(1,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12TpYvCYvTpkap21i2(i1,i2) = TpYvCYvTpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap22i1Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*TpYvCYvTpkap2(2,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12TpYvCYvTpkap22i2(i1,i2) = TpYvCYvTpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap23i1Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*TpYvCYvTpkap2(3,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12TpYvCYvTpkap23i2(i1,i2) = TpYvCYvTpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap31i1Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*TpYvCYvTpkap3(1,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13TpYvCYvTpkap31i2(i1,i2) = TpYvCYvTpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap32i1Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*TpYvCYvTpkap3(2,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13TpYvCYvTpkap32i2(i1,i2) = TpYvCYvTpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvTpkap33i1Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*TpYvCYvTpkap3(3,i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13TpYvCYvTpkap33i2(i1,i2) = TpYvCYvTpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCTvTlami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvCTvTlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpTvCYvlami1CTlami2(i1,i2) = Conjg(Tlam(i2))*TpTvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvClami1TpTvCYvlami2(i1,i2) = TpTvCYvlam(i2)*YvClam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpTvCTvlami1Clami2(i1,i2) = Conjg(lam(i2))*TpTvCTvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap1Ckap1TpTv1i1lami2(i1,i2) = Tpkap1Ckap1TpTv(1,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap1Ckap2TpTv2i1lami2(i1,i2) = Tpkap1Ckap2TpTv(2,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap1Ckap3TpTv3i1lami2(i1,i2) = Tpkap1Ckap3TpTv(3,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap2Ckap1TpTv1i1lami2(i1,i2) = Tpkap2Ckap1TpTv(1,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap2Ckap2TpTv2i1lami2(i1,i2) = Tpkap2Ckap2TpTv(2,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap2Ckap3TpTv3i1lami2(i1,i2) = Tpkap2Ckap3TpTv(3,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap3Ckap1TpTv1i1lami2(i1,i2) = Tpkap3Ckap1TpTv(1,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap3Ckap2TpTv2i1lami2(i1,i2) = Tpkap3Ckap2TpTv(2,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpkap3Ckap3TpTv3i1lami2(i1,i2) = Tpkap3Ckap3TpTv(3,i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1adjYvTvi21(i1,i2) = kap1adjYvTv(i2,1)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1adjYvTvi22(i1,i2) = kap1adjYvTv(i2,2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1adjYvTvi23(i1,i2) = kap1adjYvTv(i2,3)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1adjkap1Tpkap11i2(i1,i2) = kap1adjkap1Tpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1adjkap1Tpkap12i2(i1,i2) = kap1adjkap1Tpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1adjkap1Tpkap13i2(i1,i2) = kap1adjkap1Tpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap1adjkap1Tpkap21i2(i1,i2) = kap1adjkap1Tpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap1adjkap1Tpkap22i2(i1,i2) = kap1adjkap1Tpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap1adjkap1Tpkap23i2(i1,i2) = kap1adjkap1Tpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap1adjkap1Tpkap31i2(i1,i2) = kap1adjkap1Tpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap1adjkap1Tpkap32i2(i1,i2) = kap1adjkap1Tpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap1adjkap1Tpkap33i2(i1,i2) = kap1adjkap1Tpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Clami1Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Clami1Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Clami1Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap1i11Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap1i11Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap2i11Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap2i11Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap3i11Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap1kap3i11Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Ckap1Tk1i21(i1,i2) = kap1Ckap1Tk1(i2,1)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Ckap1Tk2i21(i1,i2) = kap1Ckap1Tk2(i2,1)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Ckap1Tk3i21(i1,i2) = kap1Ckap1Tk3(i2,1)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap1i12Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap1i12Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap2i12Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap2i12Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap3i12Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap2kap3i12Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Ckap2Tk1i22(i1,i2) = kap1Ckap2Tk1(i2,2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Ckap2Tk2i22(i1,i2) = kap1Ckap2Tk2(i2,2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Ckap2Tk3i22(i1,i2) = kap1Ckap2Tk3(i2,2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap1i13Cmv2Ckap1i21(i1,i2) = Cmv2Ckap1(i2,1)*kap1Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap1i13Ckap1mv2i21(i1,i2) = Ckap1mv2(i2,1)*kap1Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap2i13Cmv2Ckap2i21(i1,i2) = Cmv2Ckap2(i2,1)*kap1Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap2i13Ckap2mv2i21(i1,i2) = Ckap2mv2(i2,1)*kap1Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap3i13Cmv2Ckap3i21(i1,i2) = Cmv2Ckap3(i2,1)*kap1Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Ckap3kap3i13Ckap3mv2i21(i1,i2) = Ckap3mv2(i2,1)*kap1Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap1Ckap3Tk1i23(i1,i2) = kap1Ckap3Tk1(i2,3)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap1Ckap3Tk2i23(i1,i2) = kap1Ckap3Tk2(i2,3)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap1Ckap3Tk3i23(i1,i2) = kap1Ckap3Tk3(i2,3)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2adjYvTvi21(i1,i2) = kap2adjYvTv(i2,1)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2adjYvTvi22(i1,i2) = kap2adjYvTv(i2,2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2adjYvTvi23(i1,i2) = kap2adjYvTv(i2,3)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap2adjkap2Tpkap11i2(i1,i2) = kap2adjkap2Tpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap2adjkap2Tpkap12i2(i1,i2) = kap2adjkap2Tpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap2adjkap2Tpkap13i2(i1,i2) = kap2adjkap2Tpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2adjkap2Tpkap21i2(i1,i2) = kap2adjkap2Tpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2adjkap2Tpkap22i2(i1,i2) = kap2adjkap2Tpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2adjkap2Tpkap23i2(i1,i2) = kap2adjkap2Tpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap2adjkap2Tpkap31i2(i1,i2) = kap2adjkap2Tpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap2adjkap2Tpkap32i2(i1,i2) = kap2adjkap2Tpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap2adjkap2Tpkap33i2(i1,i2) = kap2adjkap2Tpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Clami1Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Clami1Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Clami1Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap1i11Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap1i11Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap2i11Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap2i11Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap3i11Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap1kap3i11Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Ckap1Tk1i21(i1,i2) = kap2Ckap1Tk1(i2,1)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Ckap1Tk2i21(i1,i2) = kap2Ckap1Tk2(i2,1)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Ckap1Tk3i21(i1,i2) = kap2Ckap1Tk3(i2,1)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap1i12Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap1i12Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap2i12Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap2i12Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap3i12Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap2kap3i12Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Ckap2Tk1i22(i1,i2) = kap2Ckap2Tk1(i2,2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Ckap2Tk2i22(i1,i2) = kap2Ckap2Tk2(i2,2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Ckap2Tk3i22(i1,i2) = kap2Ckap2Tk3(i2,2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap1i13Cmv2Ckap1i22(i1,i2) = Cmv2Ckap1(i2,2)*kap2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap1i13Ckap1mv2i22(i1,i2) = Ckap1mv2(i2,2)*kap2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap2i13Cmv2Ckap2i22(i1,i2) = Cmv2Ckap2(i2,2)*kap2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap2i13Ckap2mv2i22(i1,i2) = Ckap2mv2(i2,2)*kap2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap3i13Cmv2Ckap3i22(i1,i2) = Cmv2Ckap3(i2,2)*kap2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Ckap3kap3i13Ckap3mv2i22(i1,i2) = Ckap3mv2(i2,2)*kap2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap2Ckap3Tk1i23(i1,i2) = kap2Ckap3Tk1(i2,3)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap2Ckap3Tk2i23(i1,i2) = kap2Ckap3Tk2(i2,3)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap2Ckap3Tk3i23(i1,i2) = kap2Ckap3Tk3(i2,3)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3adjYvTvi21(i1,i2) = kap3adjYvTv(i2,1)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3adjYvTvi22(i1,i2) = kap3adjYvTv(i2,2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3adjYvTvi23(i1,i2) = kap3adjYvTv(i2,3)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11kap3adjkap3Tpkap11i2(i1,i2) = kap3adjkap3Tpkap1(1,i2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11kap3adjkap3Tpkap12i2(i1,i2) = kap3adjkap3Tpkap1(2,i2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11kap3adjkap3Tpkap13i2(i1,i2) = kap3adjkap3Tpkap1(3,i2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12kap3adjkap3Tpkap21i2(i1,i2) = kap3adjkap3Tpkap2(1,i2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12kap3adjkap3Tpkap22i2(i1,i2) = kap3adjkap3Tpkap2(2,i2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12kap3adjkap3Tpkap23i2(i1,i2) = kap3adjkap3Tpkap2(3,i2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3adjkap3Tpkap31i2(i1,i2) = kap3adjkap3Tpkap3(1,i2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3adjkap3Tpkap32i2(i1,i2) = kap3adjkap3Tpkap3(2,i2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3adjkap3Tpkap33i2(i1,i2) = kap3adjkap3Tpkap3(3,i2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Clami1Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Clami1Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Clami1Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Clam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap1i11Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap1i11Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap2i11Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap2i11Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap3i11Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap1kap3i11Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Ckap1Tk1i21(i1,i2) = kap3Ckap1Tk1(i2,1)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Ckap1Tk2i21(i1,i2) = kap3Ckap1Tk2(i2,1)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Ckap1Tk3i21(i1,i2) = kap3Ckap1Tk3(i2,1)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap1i12Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap1i12Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap2i12Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap2i12Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap3i12Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap2kap3i12Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Ckap2Tk1i22(i1,i2) = kap3Ckap2Tk1(i2,2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Ckap2Tk2i22(i1,i2) = kap3Ckap2Tk2(i2,2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Ckap2Tk3i22(i1,i2) = kap3Ckap2Tk3(i2,2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap1i13Cmv2Ckap1i23(i1,i2) = Cmv2Ckap1(i2,3)*kap3Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap1i13Ckap1mv2i23(i1,i2) = Ckap1mv2(i2,3)*kap3Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap2i13Cmv2Ckap2i23(i1,i2) = Cmv2Ckap2(i2,3)*kap3Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap2i13Ckap2mv2i23(i1,i2) = Ckap2mv2(i2,3)*kap3Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap3i13Cmv2Ckap3i23(i1,i2) = Cmv2Ckap3(i2,3)*kap3Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Ckap3kap3i13Ckap3mv2i23(i1,i2) = Ckap3mv2(i2,3)*kap3Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13kap3Ckap3Tk1i23(i1,i2) = kap3Ckap3Tk1(i2,3)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13kap3Ckap3Tk2i23(i1,i2) = kap3Ckap3Tk2(i2,3)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13kap3Ckap3Tk3i23(i1,i2) = kap3Ckap3Tk3(i2,3)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11Tk1Ckap1kap1i21(i1,i2) = Tk1Ckap1kap1(i2,1)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12Tk1Ckap1kap2i21(i1,i2) = Tk1Ckap1kap2(i2,1)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13Tk1Ckap1kap3i21(i1,i2) = Tk1Ckap1kap3(i2,1)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11Tk1Ckap2kap1i22(i1,i2) = Tk1Ckap2kap1(i2,2)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12Tk1Ckap2kap2i22(i1,i2) = Tk1Ckap2kap2(i2,2)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13Tk1Ckap2kap3i22(i1,i2) = Tk1Ckap2kap3(i2,2)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i11Tk1Ckap3kap1i23(i1,i2) = Tk1Ckap3kap1(i2,3)*YvCkap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i12Tk1Ckap3kap2i23(i1,i2) = Tk1Ckap3kap2(i2,3)*YvCkap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1i13Tk1Ckap3kap3i23(i1,i2) = Tk1Ckap3kap3(i2,3)*YvCkap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11Tk2Ckap1kap1i21(i1,i2) = Tk2Ckap1kap1(i2,1)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12Tk2Ckap1kap2i21(i1,i2) = Tk2Ckap1kap2(i2,1)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13Tk2Ckap1kap3i21(i1,i2) = Tk2Ckap1kap3(i2,1)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11Tk2Ckap2kap1i22(i1,i2) = Tk2Ckap2kap1(i2,2)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12Tk2Ckap2kap2i22(i1,i2) = Tk2Ckap2kap2(i2,2)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13Tk2Ckap2kap3i22(i1,i2) = Tk2Ckap2kap3(i2,2)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i11Tk2Ckap3kap1i23(i1,i2) = Tk2Ckap3kap1(i2,3)*YvCkap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i12Tk2Ckap3kap2i23(i1,i2) = Tk2Ckap3kap2(i2,3)*YvCkap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2i13Tk2Ckap3kap3i23(i1,i2) = Tk2Ckap3kap3(i2,3)*YvCkap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11Tk3Ckap1kap1i21(i1,i2) = Tk3Ckap1kap1(i2,1)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12Tk3Ckap1kap2i21(i1,i2) = Tk3Ckap1kap2(i2,1)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13Tk3Ckap1kap3i21(i1,i2) = Tk3Ckap1kap3(i2,1)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11Tk3Ckap2kap1i22(i1,i2) = Tk3Ckap2kap1(i2,2)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12Tk3Ckap2kap2i22(i1,i2) = Tk3Ckap2kap2(i2,2)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13Tk3Ckap2kap3i22(i1,i2) = Tk3Ckap2kap3(i2,2)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i11Tk3Ckap3kap1i23(i1,i2) = Tk3Ckap3kap1(i2,3)*YvCkap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i12Tk3Ckap3kap2i23(i1,i2) = Tk3Ckap3kap2(i2,3)*YvCkap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3i13Tk3Ckap3kap3i23(i1,i2) = Tk3Ckap3kap3(i2,3)*YvCkap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dymv2TpYvCYvlami1Clami2(i1,i2) = Conjg(lam(i2))*mv2TpYvCYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjYvYvClami1lami2(i1,i2) = YvadjYvYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjYvTvClami1lami2(i1,i2) = YvadjYvTvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap1mv2kap1i21(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap1mv2kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap1mv2kap2i21(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap1mv2kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap1mv2kap3i21(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap1mv2kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap2mv2kap1i22(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap2mv2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap2mv2kap2i22(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap2mv2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap2mv2kap3i22(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap2mv2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Yvadjkap3mv2kap1i23(i1,i2) = Conjg(Yv(i1,1))*Yvadjkap3mv2kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Yvadjkap3mv2kap2i23(i1,i2) = Conjg(Yv(i1,2))*Yvadjkap3mv2kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Yvadjkap3mv2kap3i23(i1,i2) = Conjg(Yv(i1,3))*Yvadjkap3mv2kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11YvCmv2Ckap1kap1i21(i1,i2) = Conjg(Yv(i1,1))*YvCmv2Ckap1kap1(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12YvCmv2Ckap1kap2i21(i1,i2) = Conjg(Yv(i1,2))*YvCmv2Ckap1kap2(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13YvCmv2Ckap1kap3i21(i1,i2) = Conjg(Yv(i1,3))*YvCmv2Ckap1kap3(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11YvCmv2Ckap2kap1i22(i1,i2) = Conjg(Yv(i1,1))*YvCmv2Ckap2kap1(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12YvCmv2Ckap2kap2i22(i1,i2) = Conjg(Yv(i1,2))*YvCmv2Ckap2kap2(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13YvCmv2Ckap2kap3i22(i1,i2) = Conjg(Yv(i1,3))*YvCmv2Ckap2kap3(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11YvCmv2Ckap3kap1i23(i1,i2) = Conjg(Yv(i1,1))*YvCmv2Ckap3kap1(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12YvCmv2Ckap3kap2i23(i1,i2) = Conjg(Yv(i1,2))*YvCmv2Ckap3kap2(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13YvCmv2Ckap3kap3i23(i1,i2) = Conjg(Yv(i1,3))*YvCmv2Ckap3kap3(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap1Tpkap1Clami1lami2(i1,i2) = YvCkap1Tpkap1Clam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap2Tpkap2Clami1lami2(i1,i2) = YvCkap2Tpkap2Clam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvCkap3Tpkap3Clami1lami2(i1,i2) = YvCkap3Tpkap3Clam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyadjYeYeCYvlami1mlHd2i2(i1,i2) = adjYeYeCYvlam(i1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvYvadjYvmlHd2i2(i1,i2) = adjYvYvadjYvmlHd2(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvYvCmv2Clami2(i1,i2) = adjYvYvCmv2Clam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1adjYvCml2YvClami2(i1,i2) = adjYvCml2YvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dylami1Cmv2adjYvYvClami2(i1,i2) = Cmv2adjYvYvClam(i2)*lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvTpYvCYvlami1mlHd2i2(i1,i2) = CYvTpYvCYvlam(i1)*mlHd2(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvadjYvYvClami1lami2(i1,i2) = TvadjYvYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYeCYeYvClami1lami2(i1,i2) = TpYeCYeYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYeCYeTvClami1lami2(i1,i2) = TpYeCYeTvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvml2CYvlami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvml2CYvlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYvCYvmv2lami1Clami2(i1,i2) = Conjg(lam(i2))*TpYvCYvmv2lam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpTeCYeYvClami1lami2(i1,i2) = TpTeCYeYvClam(i1)*lam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap1kap1i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap1kap2i11Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap1kap3i11Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap2kap1i12Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap2kap2i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap2kap3i12Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap3kap1i13Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*kap1Cmv2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap3kap2i13Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*kap1Cmv2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1Cmv2Ckap3kap3i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*kap1Cmv2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap1kap1i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap1kap2i11Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap1kap3i11Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap2kap1i12Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap2kap2i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap2kap3i12Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap3kap1i13Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*kap2Cmv2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap3kap2i13Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*kap2Cmv2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2Cmv2Ckap3kap3i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*kap2Cmv2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap1kap1i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Ckap1kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap1kap2i11Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Ckap1kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap1kap3i11Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Ckap1kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap2kap1i12Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Ckap2kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap2kap2i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Ckap2kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap2kap3i12Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Ckap2kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap3kap1i13Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*kap3Cmv2Ckap3kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap3kap2i13Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*kap3Cmv2Ckap3kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3Cmv2Ckap3kap3i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*kap3Cmv2Ckap3kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1i22(i1,i2) = kap1Ckap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1i23(i1,i2) = kap1Ckap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1i21(i1,i2) = kap2Ckap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1i23(i1,i2) = kap2Ckap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1i21(i1,i2) = kap3Ckap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1i22(i1,i2) = kap3Ckap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1Ckap1lami2(i1,i2) = kap1Ckap1lam(i2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2Ckap1lami2(i1,i2) = kap2Ckap1lam(i2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3Ckap1lami2(i1,i2) = kap3Ckap1lam(i2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap1adjYvYvCkap1i21(i1,i2) = kap1adjYvYvCkap1(i2,1)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap2adjYvYvCkap1i22(i1,i2) = kap2adjYvYvCkap1(i2,2)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi11kap3adjYvYvCkap1i23(i1,i2) = kap3adjYvYvCkap1(i2,3)*Tv(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2i22(i1,i2) = kap1Ckap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2i23(i1,i2) = kap1Ckap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2i21(i1,i2) = kap2Ckap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2i23(i1,i2) = kap2Ckap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2i21(i1,i2) = kap3Ckap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2i22(i1,i2) = kap3Ckap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1Ckap2lami2(i1,i2) = kap1Ckap2lam(i2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2Ckap2lami2(i1,i2) = kap2Ckap2lam(i2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3Ckap2lami2(i1,i2) = kap3Ckap2lam(i2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap1adjYvYvCkap2i21(i1,i2) = kap1adjYvYvCkap2(i2,1)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap2adjYvYvCkap2i22(i1,i2) = kap2adjYvYvCkap2(i2,2)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi12kap3adjYvYvCkap2i23(i1,i2) = kap3adjYvYvCkap2(i2,3)*Tv(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3i22(i1,i2) = kap1Ckap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3i23(i1,i2) = kap1Ckap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3i21(i1,i2) = kap2Ckap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3i23(i1,i2) = kap2Ckap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3i21(i1,i2) = kap3Ckap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3i22(i1,i2) = kap3Ckap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1Ckap3lami2(i1,i2) = kap1Ckap3lam(i2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2Ckap3lami2(i1,i2) = kap2Ckap3lam(i2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3Ckap3lami2(i1,i2) = kap3Ckap3lam(i2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap1adjYvYvCkap3i21(i1,i2) = kap1adjYvYvCkap3(i2,1)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap2adjYvYvCkap3i22(i1,i2) = kap2adjYvYvCkap3(i2,2)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTvi13kap3adjYvYvCkap3i23(i1,i2) = kap3adjYvYvCkap3(i2,3)*Tv(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Tvi21(i1,i2) = Conjg(Yv(i1,1))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Tvi21(i1,i2) = Conjg(Yv(i1,2))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Tvi21(i1,i2) = Conjg(Yv(i1,3))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11Tvi21(i1,i2) = Conjg(Tv(i1,1))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12Tvi21(i1,i2) = Conjg(Tv(i1,2))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13Tvi21(i1,i2) = Conjg(Tv(i1,3))*Tv(i2,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Tvi22(i1,i2) = Conjg(Yv(i1,1))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Tvi22(i1,i2) = Conjg(Yv(i1,2))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Tvi22(i1,i2) = Conjg(Yv(i1,3))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11Tvi22(i1,i2) = Conjg(Tv(i1,1))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12Tvi22(i1,i2) = Conjg(Tv(i1,2))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13Tvi22(i1,i2) = Conjg(Tv(i1,3))*Tv(i2,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi11Tvi23(i1,i2) = Conjg(Yv(i1,1))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi12Tvi23(i1,i2) = Conjg(Yv(i1,2))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCYvi13Tvi23(i1,i2) = Conjg(Yv(i1,3))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi11Tvi23(i1,i2) = Conjg(Tv(i1,1))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi12Tvi23(i1,i2) = Conjg(Tv(i1,2))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyCTvi13Tvi23(i1,i2) = Conjg(Tv(i1,3))*Tv(i2,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1Clami2(i1,i2) = Conjg(lam(i2))*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1adjYvYvCTlami2(i1,i2) = adjYvYvCTlam(i2)*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTlami1adjTvYvClami2(i1,i2) = adjTvYvClam(i2)*Tlam(i1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap1i11Tlami2(i1,i2) = Yvadjkap1kap1(i1,1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap2i11Tlami2(i1,i2) = Yvadjkap1kap2(i1,1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap1kap3i11Tlami2(i1,i2) = Yvadjkap1kap3(i1,1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap1i12Tlami2(i1,i2) = Yvadjkap2kap1(i1,2)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap2i12Tlami2(i1,i2) = Yvadjkap2kap2(i1,2)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap2kap3i12Tlami2(i1,i2) = Yvadjkap2kap3(i1,2)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap1i13Tlami2(i1,i2) = Yvadjkap3kap1(i1,3)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap2i13Tlami2(i1,i2) = Yvadjkap3kap2(i1,3)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjkap3kap3i13Tlami2(i1,i2) = Yvadjkap3kap3(i1,3)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyYvadjYvYvClami1Tlami2(i1,i2) = YvadjYvYvClam(i1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTpYeCYeYvClami1Tlami2(i1,i2) = TpYeCYeYvClam(i1)*Tlam(i2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i11Ckap1TpYvCYvi21(i1,i2) = Ckap1TpYvCYv(i2,1)*kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i11adjYvCml2Yvadjkap11i2(i1,i2) = adjYvCml2Yvadjkap1(1,i2)*kap1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i12Ckap2TpYvCYvi21(i1,i2) = Ckap2TpYvCYv(i2,1)*kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i12adjYvCml2Yvadjkap21i2(i1,i2) = adjYvCml2Yvadjkap2(1,i2)*kap1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i13Ckap3TpYvCYvi21(i1,i2) = Ckap3TpYvCYv(i2,1)*kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap1i13adjYvCml2Yvadjkap31i2(i1,i2) = adjYvCml2Yvadjkap3(1,i2)*kap1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i11Ckap1TpYvCYvi22(i1,i2) = Ckap1TpYvCYv(i2,2)*kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i11adjYvCml2Yvadjkap12i2(i1,i2) = adjYvCml2Yvadjkap1(2,i2)*kap2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i12Ckap2TpYvCYvi22(i1,i2) = Ckap2TpYvCYv(i2,2)*kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i12adjYvCml2Yvadjkap22i2(i1,i2) = adjYvCml2Yvadjkap2(2,i2)*kap2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i13Ckap3TpYvCYvi22(i1,i2) = Ckap3TpYvCYv(i2,2)*kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap2i13adjYvCml2Yvadjkap32i2(i1,i2) = adjYvCml2Yvadjkap3(2,i2)*kap2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i11Ckap1TpYvCYvi23(i1,i2) = Ckap1TpYvCYv(i2,3)*kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i11adjYvCml2Yvadjkap13i2(i1,i2) = adjYvCml2Yvadjkap1(3,i2)*kap3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i12Ckap2TpYvCYvi23(i1,i2) = Ckap2TpYvCYv(i2,3)*kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i12adjYvCml2Yvadjkap23i2(i1,i2) = adjYvCml2Yvadjkap2(3,i2)*kap3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i13Ckap3TpYvCYvi23(i1,i2) = Ckap3TpYvCYv(i2,3)*kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
Dykap3i13adjYvCml2Yvadjkap33i2(i1,i2) = adjYvCml2Yvadjkap3(3,i2)*kap3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11CTk1i21(i1,i2) = Conjg(Tk1(i2,1))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11CTk1i22(i1,i2) = Conjg(Tk1(i2,2))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11CTk1i23(i1,i2) = Conjg(Tk1(i2,3))*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1lami2(i1,i2) = Ckap1lam(i2)*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i11Ckap1TpYvCTvi21(i1,i2) = Ckap1TpYvCTv(i2,1)*Tk1(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12CTk1i21(i1,i2) = Conjg(Tk1(i2,1))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12CTk1i22(i1,i2) = Conjg(Tk1(i2,2))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12CTk1i23(i1,i2) = Conjg(Tk1(i2,3))*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2lami2(i1,i2) = Ckap2lam(i2)*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i12Ckap2TpYvCTvi21(i1,i2) = Ckap2TpYvCTv(i2,1)*Tk1(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13CTk1i21(i1,i2) = Conjg(Tk1(i2,1))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13CTk1i22(i1,i2) = Conjg(Tk1(i2,2))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13CTk1i23(i1,i2) = Conjg(Tk1(i2,3))*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3lami2(i1,i2) = Ckap3lam(i2)*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk1i13Ckap3TpYvCTvi21(i1,i2) = Ckap3TpYvCTv(i2,1)*Tk1(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11CTk2i21(i1,i2) = Conjg(Tk2(i2,1))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11CTk2i22(i1,i2) = Conjg(Tk2(i2,2))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11CTk2i23(i1,i2) = Conjg(Tk2(i2,3))*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1lami2(i1,i2) = Ckap1lam(i2)*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i11Ckap1TpYvCTvi22(i1,i2) = Ckap1TpYvCTv(i2,2)*Tk2(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12CTk2i21(i1,i2) = Conjg(Tk2(i2,1))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12CTk2i22(i1,i2) = Conjg(Tk2(i2,2))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12CTk2i23(i1,i2) = Conjg(Tk2(i2,3))*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2lami2(i1,i2) = Ckap2lam(i2)*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i12Ckap2TpYvCTvi22(i1,i2) = Ckap2TpYvCTv(i2,2)*Tk2(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13CTk2i21(i1,i2) = Conjg(Tk2(i2,1))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13CTk2i22(i1,i2) = Conjg(Tk2(i2,2))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13CTk2i23(i1,i2) = Conjg(Tk2(i2,3))*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3lami2(i1,i2) = Ckap3lam(i2)*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk2i13Ckap3TpYvCTvi22(i1,i2) = Ckap3TpYvCTv(i2,2)*Tk2(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1i21(i1,i2) = Conjg(kap1(i2,1))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1i22(i1,i2) = Conjg(kap1(i2,2))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1i23(i1,i2) = Conjg(kap1(i2,3))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11CTk3i21(i1,i2) = Conjg(Tk3(i2,1))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11CTk3i22(i1,i2) = Conjg(Tk3(i2,2))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11CTk3i23(i1,i2) = Conjg(Tk3(i2,3))*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1lami2(i1,i2) = Ckap1lam(i2)*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i11Ckap1TpYvCTvi23(i1,i2) = Ckap1TpYvCTv(i2,3)*Tk3(i1,1) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2i21(i1,i2) = Conjg(kap2(i2,1))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2i22(i1,i2) = Conjg(kap2(i2,2))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2i23(i1,i2) = Conjg(kap2(i2,3))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12CTk3i21(i1,i2) = Conjg(Tk3(i2,1))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12CTk3i22(i1,i2) = Conjg(Tk3(i2,2))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12CTk3i23(i1,i2) = Conjg(Tk3(i2,3))*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2lami2(i1,i2) = Ckap2lam(i2)*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i12Ckap2TpYvCTvi23(i1,i2) = Ckap2TpYvCTv(i2,3)*Tk3(i1,2) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3i21(i1,i2) = Conjg(kap3(i2,1))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3i22(i1,i2) = Conjg(kap3(i2,2))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3i23(i1,i2) = Conjg(kap3(i2,3))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13CTk3i21(i1,i2) = Conjg(Tk3(i2,1))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13CTk3i22(i1,i2) = Conjg(Tk3(i2,2))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13CTk3i23(i1,i2) = Conjg(Tk3(i2,3))*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3lami2(i1,i2) = Ckap3lam(i2)*Tk3(i1,3) 
  End Do 
End Do 
Do i1=1,3
  Do i2=1,3
DyTk3i13Ckap3TpYvCTvi23(i1,i2) = Ckap3TpYvCTv(i2,3)*Tk3(i1,3) 
  End Do 
End Do 
End If 
 
 
Tr1(1) = g1*sqrt3ov5*(-1._dp*(mHd2) + mHu2 + Trmd2 + Trme2 - Trml2 + Trmq2 -          & 
&  2._dp*(Trmu2))

If (TwoLoopRGE) Then 
Tr2U1(1, 1) = (g1p2*(3._dp*(mHd2) + 3._dp*(mHu2) + 2._dp*(Trmd2) + 6._dp*(Trme2)      & 
&  + 3._dp*(Trml2) + Trmq2 + 8._dp*(Trmu2)))/10._dp

Tr3(1) = (g1*ooSqrt15*(-9*g1p2*mHd2 - 45*g2p2*mHd2 + 9*g1p2*mHu2 + 45*g2p2*mHu2 -     & 
&  30._dp*(SPClamxxTpYvmlHd2) + 30*(mHd2 - mHu2)*SPlamxxClam + 4*g1p2*Trmd2 +            & 
&  80*g3p2*Trmd2 + 36*g1p2*Trme2 - 9*g1p2*Trml2 - 45*g2p2*Trml2 + 30._dp*(Trml2YvadjYv)  & 
&  + g1p2*Trmq2 + 45*g2p2*Trmq2 + 80*g3p2*Trmq2 - 32*g1p2*Trmu2 - 160*g3p2*Trmu2 +       & 
&  90*mHd2*TrYdadjYd - 60._dp*(TrYdadjYdCmd2) - 30._dp*(TrYdCmq2adjYd) + 30*mHd2*TrYeadjYe -& 
&  60._dp*(TrYeadjYeCme2) + 30._dp*(TrYeCml2adjYe) - 90*mHu2*TrYuadjYu + 120._dp*(TrYuadjYuCmu2)& 
&  - 30._dp*(TrYuCmq2adjYu) - 30*mHu2*TrYvadjYv))/20._dp

Tr2(2) = (mHd2 + mHu2 + Trml2 + 3._dp*(Trmq2))/2._dp

Tr2(3) = (Trmd2 + 2._dp*(Trmq2) + Trmu2)/2._dp

End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = 33._dp*(g1p3)/5._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (199._dp*(g1p5) + 5*g1p3*(27._dp*(g2p2) + 88._dp*(g3p2) - 6._dp*(SPlamxxClam) -       & 
&  14._dp*(TrYdadjYd) - 18._dp*(TrYeadjYe) - 26._dp*(TrYuadjYu) - 6._dp*(TrYvadjYv)))/25._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = g2p3

 
 
If (TwoLoopRGE) Then 
betag22 = (9*g1p2*g2p3)/5._dp + 25._dp*(g2p5) + 2*g2p3*(12._dp*(g3p2) - SPlamxxClam -           & 
&  3._dp*(TrYdadjYd) - TrYeadjYe - 3._dp*(TrYuadjYu) - TrYvadjYv)

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = -3._dp*(g3p3)

 
 
If (TwoLoopRGE) Then 
betag32 = (11*g1p2*g3p3)/5._dp + 9*g2p2*g3p3 + 14._dp*(g3p5) - 4*g3p3*TrYdadjYd -               & 
&  4*g3p3*TrYuadjYu

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = (-7._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)/3._dp +               & 
&  SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd) + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = ((287._dp*(g1p4) + 2*g1p2*(45._dp*(g2p2) + 40._dp*(g3p2) - 18._dp*(TrYdadjYd) +       & 
&  54._dp*(TrYeadjYe)) + 5*(135._dp*(g2p4) + 2*(72*g2p2*g3p2 - 16._dp*(g3p4) -           & 
&  9*(2._dp*(SPlamxxadjYvYvClam) + 3._dp*(SPlamxxClamp2) + 2._dp*(SPTpkap1Clamxxadjkap1lam) +& 
&  2._dp*(SPTpkap2Clamxxadjkap2lam) + 2._dp*(SPTpkap3Clamxxadjkap3lam) - 16*g3p2*TrYdadjYd +& 
&  9._dp*(TrYdadjYdYdadjYd) + 3._dp*(TrYdadjYuYuadjYd) + 3._dp*(TrYeadjYeYeadjYe) +      & 
&  3*SPlamxxClam*TrYuadjYu + SPlamxxClam*TrYvadjYv + TrYvadjYvTpYeCYe))))*Yd)/90._dp +   & 
&  (4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 3*(SPlamxxClam + 3._dp*(TrYdadjYd) +             & 
&  TrYeadjYe))*YdadjYdYd - 4._dp*(YdadjYdYdadjYdYd) + (4*g1p2*YdadjYuYu)/5._dp -         & 
&  SPlamxxClam*YdadjYuYu - 3*TrYuadjYu*YdadjYuYu - TrYvadjYv*YdadjYuYu - 2._dp*(YdadjYuYuadjYdYd) -& 
&  2._dp*(YdadjYuYuadjYuYu)

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
Call Chop(DYd) 

!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (-9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + SPlamxxClam + 3._dp*(TrYdadjYd)      & 
&  + TrYeadjYe)*Ye + 3._dp*(YeadjYeYe) + YeCYvTpYv

 
 
If (TwoLoopRGE) Then 
betaYe2 = -2._dp*(DyYeCYvlami1YvClami2) + ((135._dp*(g1p4) + 2*g1p2*(9._dp*(g2p2) -             & 
&  2._dp*(TrYdadjYd) + 6._dp*(TrYeadjYe)) + 5*(15._dp*(g2p4) - 2*(2._dp*(SPlamxxadjYvYvClam) +& 
&  3._dp*(SPlamxxClamp2) + 2._dp*(SPTpkap1Clamxxadjkap1lam) + 2._dp*(SPTpkap2Clamxxadjkap2lam) +& 
&  2._dp*(SPTpkap3Clamxxadjkap3lam) - 16*g3p2*TrYdadjYd + 9._dp*(TrYdadjYdYdadjYd) +     & 
&  3._dp*(TrYdadjYuYuadjYd) + 3._dp*(TrYeadjYeYeadjYe) + 3*SPlamxxClam*TrYuadjYu +       & 
&  SPlamxxClam*TrYvadjYv + TrYvadjYvTpYeCYe)))*Ye)/10._dp + 6*g2p2*YeadjYeYe -           & 
&  3*SPlamxxClam*YeadjYeYe - 9*TrYdadjYd*YeadjYeYe - 3*TrYeadjYe*YeadjYeYe -             & 
&  4._dp*(YeadjYeYeadjYeYe) - 2._dp*(YeCYvkap1adjkap1TpYv) - 2._dp*(YeCYvkap2adjkap2TpYv) -& 
&  2._dp*(YeCYvkap3adjkap3TpYv) - SPlamxxClam*YeCYvTpYv - 3*TrYuadjYu*YeCYvTpYv -        & 
&  TrYvadjYv*YeCYvTpYv - 2._dp*(YeCYvTpYvadjYeYe) - 2._dp*(YeCYvTpYvCYvTpYv)

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
Call Chop(DYe) 

!-------------------- 
! lam 
!-------------------- 
 
Do i1 = 1,3
betalam1(i1) = 3*TpYvCYvlam(i1) + 2*kap1Ckap1(i1,1)*lam(1) + 2*kap2Ckap1(i1,2)*lam(1) +              & 
&  2*kap3Ckap1(i1,3)*lam(1) + 2*kap1Ckap2(i1,1)*lam(2) + 2*kap2Ckap2(i1,2)*lam(2) +      & 
&  2*kap3Ckap2(i1,3)*lam(2) + 2*kap1Ckap3(i1,1)*lam(3) + 2*kap2Ckap3(i1,2)*lam(3) +      & 
&  2*kap3Ckap3(i1,3)*lam(3) - (3*g1p2*lam(i1))/5._dp - 3*g2p2*lam(i1) + 4*SPlamxxClam*lam(i1) +& 
&  3*TrYdadjYd*lam(i1) + TrYeadjYe*lam(i1) + 3*TrYuadjYu*lam(i1) + TrYvadjYv*lam(i1)
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
betalam2(i1) = -8*adjkap1lam(1)*kap1adjkap1Tpkap1(1,i1) - 8*adjkap2lam(1)*kap1adjkap1Tpkap1(2,       & 
& i1) - 8*adjkap3lam(1)*kap1adjkap1Tpkap1(3,i1) - 8*adjkap1lam(2)*kap1adjkap1Tpkap2(1,   & 
& i1) - 8*adjkap2lam(2)*kap1adjkap1Tpkap2(2,i1) - 8*adjkap3lam(2)*kap1adjkap1Tpkap2(3,   & 
& i1) - 8*adjkap1lam(3)*kap1adjkap1Tpkap3(1,i1) - 8*adjkap2lam(3)*kap1adjkap1Tpkap3(2,   & 
& i1) - 8*adjkap3lam(3)*kap1adjkap1Tpkap3(3,i1) - 8*adjkap1lam(1)*kap2adjkap2Tpkap1(1,   & 
& i1) - 8*adjkap2lam(1)*kap2adjkap2Tpkap1(2,i1) - 8*adjkap3lam(1)*kap2adjkap2Tpkap1(3,   & 
& i1) - 8*adjkap1lam(2)*kap2adjkap2Tpkap2(1,i1) - 8*adjkap2lam(2)*kap2adjkap2Tpkap2(2,   & 
& i1) - 8*adjkap3lam(2)*kap2adjkap2Tpkap2(3,i1) - 8*adjkap1lam(3)*kap2adjkap2Tpkap3(1,   & 
& i1) - 8*adjkap2lam(3)*kap2adjkap2Tpkap3(2,i1) - 8*adjkap3lam(3)*kap2adjkap2Tpkap3(3,   & 
& i1) - 8*adjkap1lam(1)*kap3adjkap3Tpkap1(1,i1) - 8*adjkap2lam(1)*kap3adjkap3Tpkap1(2,   & 
& i1) - 8*adjkap3lam(1)*kap3adjkap3Tpkap1(3,i1) - 8*adjkap1lam(2)*kap3adjkap3Tpkap2(1,   & 
& i1) - 8*adjkap2lam(2)*kap3adjkap3Tpkap2(2,i1) - 8*adjkap3lam(2)*kap3adjkap3Tpkap2(3,   & 
& i1) - 8*adjkap1lam(3)*kap3adjkap3Tpkap3(1,i1) - 8*adjkap2lam(3)*kap3adjkap3Tpkap3(2,   & 
& i1) - 8*adjkap3lam(3)*kap3adjkap3Tpkap3(3,i1) - TpYvadjYeYeCYvlam(i1) - 2*TpYvCYvkap1adjkap1lam(i1) -& 
&  2*TpYvCYvkap2adjkap2lam(i1) - 2*TpYvCYvkap3adjkap3lam(i1) + (6*g1p2*TpYvCYvlam(i1))/5._dp +& 
&  6*g2p2*TpYvCYvlam(i1) - 7*SPlamxxClam*TpYvCYvlam(i1) - 9*TrYuadjYu*TpYvCYvlam(i1) -   & 
&  3*TrYvadjYv*TpYvCYvlam(i1) - 8*adjkap1lam(1)*TpYvCYvTpkap1(1,i1) - 8*adjkap2lam(1)*TpYvCYvTpkap1(2,& 
& i1) - 8*adjkap3lam(1)*TpYvCYvTpkap1(3,i1) - 8*adjkap1lam(2)*TpYvCYvTpkap2(1,           & 
& i1) - 8*adjkap2lam(2)*TpYvCYvTpkap2(2,i1) - 8*adjkap3lam(2)*TpYvCYvTpkap2(3,           & 
& i1) - 8*adjkap1lam(3)*TpYvCYvTpkap3(1,i1) - 8*adjkap2lam(3)*TpYvCYvTpkap3(2,           & 
& i1) - 8*adjkap3lam(3)*TpYvCYvTpkap3(3,i1) - 4*TpYvCYvTpYvCYvlam(i1) - 8*adjkap1lam(1)*kap1Clam(i1)*lam(1)  
betalam2(i1) =  betalam2(i1)- 8*adjkap1lam(2)*kap2Clam(i1)*lam(1) - 8*adjkap1lam(3)*kap3Clam(i1)*lam(1) -           & 
&  8*adjkap2lam(1)*kap1Clam(i1)*lam(2) - 8*adjkap2lam(2)*kap2Clam(i1)*lam(2) -           & 
&  8*adjkap2lam(3)*kap3Clam(i1)*lam(2) - 8*adjkap3lam(1)*kap1Clam(i1)*lam(3) -           & 
&  8*adjkap3lam(2)*kap2Clam(i1)*lam(3) - 8*adjkap3lam(3)*kap3Clam(i1)*lam(3) +           & 
&  (207*g1p4*lam(i1))/50._dp + (9*g1p2*g2p2*lam(i1))/5._dp + (15*g2p4*lam(i1))/2._dp -   & 
&  10*SPlamxxadjYvYvClam*lam(i1) + (6*g1p2*SPlamxxClam*lam(i1))/5._dp + 6*g2p2*SPlamxxClam*lam(i1) -& 
&  10*SPlamxxClamp2*lam(i1) - 4*SPTpkap1Clamxxadjkap1lam*lam(i1) - 4*SPTpkap2Clamxxadjkap2lam*lam(i1) -& 
&  4*SPTpkap3Clamxxadjkap3lam*lam(i1) - (2*g1p2*TrYdadjYd*lam(i1))/5._dp +               & 
&  16*g3p2*TrYdadjYd*lam(i1) - 9*SPlamxxClam*TrYdadjYd*lam(i1) - 9*TrYdadjYdYdadjYd*lam(i1) -& 
&  6*TrYdadjYuYuadjYd*lam(i1) + (6*g1p2*TrYeadjYe*lam(i1))/5._dp - 3*SPlamxxClam*TrYeadjYe*lam(i1) -& 
&  3*TrYeadjYeYeadjYe*lam(i1) + (4*g1p2*TrYuadjYu*lam(i1))/5._dp + 16*g3p2*TrYuadjYu*lam(i1) -& 
&  9*SPlamxxClam*TrYuadjYu*lam(i1) - 9*TrYuadjYuYuadjYu*lam(i1) - 3*SPlamxxClam*TrYvadjYv*lam(i1) -& 
&  2*TrYvadjYvTpYeCYe*lam(i1) - 3*TrYvadjYvYvadjYv*lam(i1) - 2*TrYvCkap1Tpkap1adjYv*lam(i1) -& 
&  2*TrYvCkap2Tpkap2adjYv*lam(i1) - 2*TrYvCkap3Tpkap3adjYv*lam(i1)
End Do

 
Dlam = oo16pi2*( betalam1 + oo16pi2 * betalam2 ) 

 
Else 
Dlam = oo16pi2* betalam1 
End If 
 
 
Call Chop(Dlam) 

!-------------------- 
! Yv 
!-------------------- 
 
betaYv1  = 3._dp*(DyYvClami1lami2) + TpYeCYeYv - (3*g1p2*Yv)/5._dp - 3*g2p2*Yv +      & 
&  SPlamxxClam*Yv + 3*TrYuadjYu*Yv + TrYvadjYv*Yv + 3._dp*(YvadjYvYv) + 2._dp*(YvCkap1Tpkap1)& 
&  + 2._dp*(YvCkap2Tpkap2) + 2._dp*(YvCkap3Tpkap3)

 
 
If (TwoLoopRGE) Then 
betaYv2 = DyTpYeCYeYvClami1lami2 - 4._dp*(DyYvadjYvYvClami1lami2) - 8._dp*(DyYvCkap1i11kap1adjkap1Tpkap11i2) -& 
&  8._dp*(DyYvCkap1i11kap2adjkap2Tpkap11i2) - 8._dp*(DyYvCkap1i11kap3adjkap3Tpkap11i2) - & 
&  8._dp*(DyYvCkap1i11TpYvCYvTpkap11i2) - 8._dp*(DyYvCkap1i12kap1adjkap1Tpkap21i2) -     & 
&  8._dp*(DyYvCkap1i12kap2adjkap2Tpkap21i2) - 8._dp*(DyYvCkap1i12kap3adjkap3Tpkap21i2) - & 
&  8._dp*(DyYvCkap1i12TpYvCYvTpkap21i2) - 8._dp*(DyYvCkap1i13kap1adjkap1Tpkap31i2) -     & 
&  8._dp*(DyYvCkap1i13kap2adjkap2Tpkap31i2) - 8._dp*(DyYvCkap1i13kap3adjkap3Tpkap31i2) - & 
&  8._dp*(DyYvCkap1i13TpYvCYvTpkap31i2) - 2._dp*(DyYvCkap1Tpkap1Clami1lami2) -           & 
&  8._dp*(DyYvCkap2i11kap1adjkap1Tpkap12i2) - 8._dp*(DyYvCkap2i11kap2adjkap2Tpkap12i2) - & 
&  8._dp*(DyYvCkap2i11kap3adjkap3Tpkap12i2) - 8._dp*(DyYvCkap2i11TpYvCYvTpkap12i2) -     & 
&  8._dp*(DyYvCkap2i12kap1adjkap1Tpkap22i2) - 8._dp*(DyYvCkap2i12kap2adjkap2Tpkap22i2) - & 
&  8._dp*(DyYvCkap2i12kap3adjkap3Tpkap22i2) - 8._dp*(DyYvCkap2i12TpYvCYvTpkap22i2) -     & 
&  8._dp*(DyYvCkap2i13kap1adjkap1Tpkap32i2) - 8._dp*(DyYvCkap2i13kap2adjkap2Tpkap32i2) - & 
&  8._dp*(DyYvCkap2i13kap3adjkap3Tpkap32i2) - 8._dp*(DyYvCkap2i13TpYvCYvTpkap32i2) -     & 
&  2._dp*(DyYvCkap2Tpkap2Clami1lami2) - 8._dp*(DyYvCkap3i11kap1adjkap1Tpkap13i2) -       & 
&  8._dp*(DyYvCkap3i11kap2adjkap2Tpkap13i2) - 8._dp*(DyYvCkap3i11kap3adjkap3Tpkap13i2) - & 
&  8._dp*(DyYvCkap3i11TpYvCYvTpkap13i2) - 8._dp*(DyYvCkap3i12kap1adjkap1Tpkap23i2) -     & 
&  8._dp*(DyYvCkap3i12kap2adjkap2Tpkap23i2) - 8._dp*(DyYvCkap3i12kap3adjkap3Tpkap23i2) - & 
&  8._dp*(DyYvCkap3i12TpYvCYvTpkap23i2) - 8._dp*(DyYvCkap3i13kap1adjkap1Tpkap33i2) -     & 
&  8._dp*(DyYvCkap3i13kap2adjkap2Tpkap33i2) - 8._dp*(DyYvCkap3i13kap3adjkap3Tpkap33i2) - & 
&  8._dp*(DyYvCkap3i13TpYvCYvTpkap33i2) - 2._dp*(DyYvCkap3Tpkap3Clami1lami2) -           & 
&  4._dp*(DyYvClami1TpYvCYvlami2) + (6*DyYvClami1lami2*g1p2)/5._dp + 6*DyYvClami1lami2*g2p2  
betaYv2 =  betaYv2- 7*DyYvClami1lami2*SPlamxxClam - 2._dp*(TpYeCYeTpYeCYeYv) + (6*g1p2*TpYeCYeYv)/5._dp - & 
&  SPlamxxClam*TpYeCYeYv - 6*DyYvClami1lami2*TrYdadjYd - 3*TpYeCYeYv*TrYdadjYd -         & 
&  2*DyYvClami1lami2*TrYeadjYe - TpYeCYeYv*TrYeadjYe - 9*DyYvClami1lami2*TrYuadjYu -     & 
&  3*DyYvClami1lami2*TrYvadjYv + (207*g1p4*Yv)/50._dp + (9*g1p2*g2p2*Yv)/5._dp +         & 
&  (15*g2p4*Yv)/2._dp - 6*SPlamxxadjYvYvClam*Yv - 3*SPlamxxClamp2*Yv - 2*SPTpkap1Clamxxadjkap1lam*Yv -& 
&  2*SPTpkap2Clamxxadjkap2lam*Yv - 2*SPTpkap3Clamxxadjkap3lam*Yv - 3*SPlamxxClam*TrYdadjYd*Yv -& 
&  3*TrYdadjYuYuadjYd*Yv - SPlamxxClam*TrYeadjYe*Yv + (4*g1p2*TrYuadjYu*Yv)/5._dp +      & 
&  16*g3p2*TrYuadjYu*Yv - 9*TrYuadjYuYuadjYu*Yv - TrYvadjYvTpYeCYe*Yv - 3*TrYvadjYvYvadjYv*Yv -& 
&  2*TrYvCkap1Tpkap1adjYv*Yv - 2*TrYvCkap2Tpkap2adjYv*Yv - 2*TrYvCkap3Tpkap3adjYv*Yv -   & 
&  2._dp*(YvadjYvTpYeCYeYv) + (6*g1p2*YvadjYvYv)/5._dp + 6*g2p2*YvadjYvYv -              & 
&  3*SPlamxxClam*YvadjYvYv - 9*TrYuadjYu*YvadjYvYv - 3*TrYvadjYv*YvadjYvYv -             & 
&  4._dp*(YvadjYvYvadjYvYv) - 2._dp*(YvCkap1Tpkap1adjYvYv) - 2._dp*(YvCkap2Tpkap2adjYvYv) -& 
&  2._dp*(YvCkap3Tpkap3adjYvYv) - 8*DyYvCkap1i11kap1Clami2*lam(1) - 8*DyYvCkap1i12kap2Clami2*lam(1) -& 
&  8*DyYvCkap1i13kap3Clami2*lam(1) - 8*DyYvCkap2i11kap1Clami2*lam(2) - 8*DyYvCkap2i12kap2Clami2*lam(2) -& 
&  8*DyYvCkap2i13kap3Clami2*lam(2) - 8*DyYvCkap3i11kap1Clami2*lam(3) - 8*DyYvCkap3i12kap2Clami2*lam(3) -& 
&  8*DyYvCkap3i13kap3Clami2*lam(3)

 
DYv = oo16pi2*( betaYv1 + oo16pi2 * betaYv2 ) 

 
Else 
DYv = oo16pi2* betaYv1 
End If 
 
 
Call Chop(DYv) 

!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = (-13._dp*(g1p2)/15._dp - 3._dp*(g2p2) - 16._dp*(g3p2)/3._dp +              & 
&  SPlamxxClam + 3._dp*(TrYuadjYu) + TrYvadjYv)*Yu + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = ((2743._dp*(g1p4) + 5*(g1p2*(90._dp*(g2p2) + 272._dp*(g3p2) + 72._dp*(TrYuadjYu)) +   & 
&  5*(135._dp*(g2p4) + 2*(72*g2p2*g3p2 - 16._dp*(g3p4) - 9*(6._dp*(SPlamxxadjYvYvClam) + & 
&  3._dp*(SPlamxxClamp2) + 2._dp*(SPTpkap1Clamxxadjkap1lam) + 2._dp*(SPTpkap2Clamxxadjkap2lam) +& 
&  2._dp*(SPTpkap3Clamxxadjkap3lam) + 3*SPlamxxClam*TrYdadjYd + 3._dp*(TrYdadjYuYuadjYd) +& 
&  SPlamxxClam*TrYeadjYe - 16*g3p2*TrYuadjYu + 9._dp*(TrYuadjYuYuadjYu) + TrYvadjYvTpYeCYe +& 
&  3._dp*(TrYvadjYvYvadjYv) + 2._dp*(TrYvCkap1Tpkap1adjYv) + 2._dp*(TrYvCkap2Tpkap2adjYv) +& 
&  2._dp*(TrYvCkap3Tpkap3adjYv))))))*Yu)/450._dp + ((2._dp*(g1p2) - 5*(SPlamxxClam +     & 
&  3._dp*(TrYdadjYd) + TrYeadjYe))*YuadjYdYd)/5._dp - 2._dp*(YuadjYdYdadjYdYd) -         & 
&  2._dp*(YuadjYdYdadjYuYu) + (2*g1p2*YuadjYuYu)/5._dp + 6*g2p2*YuadjYuYu -              & 
&  3*SPlamxxClam*YuadjYuYu - 9*TrYuadjYu*YuadjYuYu - 3*TrYvadjYv*YuadjYuYu -             & 
&  4._dp*(YuadjYuYuadjYuYu)

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
Call Chop(DYu) 

!-------------------- 
! kap 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betakap1(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap1(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap1(i1,i2,3) = 0
End Do
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betakap2(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap2(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betakap2(i1,i2,3) = 0
End Do
End Do

 
Dkap = oo16pi2*( betakap1 + oo16pi2 * betakap2 ) 

 
Else 
Dkap = oo16pi2* betakap1 
End If 
 
 
!-------------------- 
! Td 
!-------------------- 
 
betaTd1  = 5._dp*(TdadjYdYd) + TdadjYuYu + (14*g1p2*M1*Yd)/15._dp + (32*g3p2*M3*Yd)   & 
& /3._dp + 6*g2p2*M2*Yd + 2*SPClamxxTlam*Yd + 6*TradjYdTd*Yd + 2*TradjYeTe*Yd +          & 
&  4._dp*(YdadjYdTd) + 2._dp*(YdadjYuTu) - (7*g1p2*Td)/15._dp - 3*g2p2*Td -              & 
&  (16*g3p2*Td)/3._dp + SPlamxxClam*Td + 3*TrYdadjYd*Td + TrYeadjYe*Td

 
 
If (TwoLoopRGE) Then 
betaTd2 = -6._dp*(TdadjYdYdadjYdYd) + (4*g1p2*TdadjYuYu)/5._dp - SPlamxxClam*TdadjYuYu -        & 
&  4._dp*(TdadjYuYuadjYdYd) - 2._dp*(TdadjYuYuadjYuYu) + TdadjYdYd*(6._dp*(g1p2)/5._dp + & 
&  12._dp*(g2p2) - 5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)) - 3*TdadjYuYu*TrYuadjYu -& 
&  TdadjYuYu*TrYvadjYv - (574*g1p4*M1*Yd)/45._dp - 2*g1p2*g2p2*M1*Yd - (16*g1p2*g3p2*M1*Yd)/9._dp -& 
&  (16*g1p2*g3p2*M3*Yd)/9._dp - 16*g2p2*g3p2*M3*Yd + (64*g3p4*M3*Yd)/9._dp -             & 
&  2*g1p2*g2p2*M2*Yd - 30*g2p4*M2*Yd - 16*g2p2*g3p2*M2*Yd - 4*SPadjYvYvClamxxTlam*Yd -   & 
&  4*SPClamxxTpTvCYvlam*Yd - 12*SPClamxxTlam*SPlamxxClam*Yd - (4*g1p2*TradjYdTd*Yd)/5._dp +& 
&  32*g3p2*TradjYdTd*Yd + (12*g1p2*TradjYeTe*Yd)/5._dp - 2*TradjYeTeCYvTpYv*Yd -         & 
&  6*SPlamxxClam*TradjYuTu*Yd - 2*TradjYvTpYeCYeTv*Yd - 2*SPlamxxClam*TradjYvTv*Yd +     & 
&  (4*g1p2*M1*TrYdadjYd*Yd)/5._dp - 32*g3p2*M3*TrYdadjYd*Yd - 36*TrYdadjYdTdadjYd*Yd -   & 
&  6*TrYdadjYuTuadjYd*Yd - (12*g1p2*M1*TrYeadjYe*Yd)/5._dp - 12*TrYeadjYeTeadjYe*Yd -    & 
&  6*TrYuadjYdTdadjYu*Yd - 6*SPClamxxTlam*TrYuadjYu*Yd - 2*SPClamxxTlam*TrYvadjYv*Yd +   & 
&  (6*g1p2*YdadjYdTd)/5._dp + 6*g2p2*YdadjYdTd - 4*SPlamxxClam*YdadjYdTd -               & 
&  12*TrYdadjYd*YdadjYdTd - 4*TrYeadjYe*YdadjYdTd - 8._dp*(YdadjYdTdadjYdYd) -           & 
&  (8*g1p2*M1*YdadjYdYd)/5._dp - 12*g2p2*M2*YdadjYdYd - 6*SPClamxxTlam*YdadjYdYd -       & 
&  18*TradjYdTd*YdadjYdYd - 6*TradjYeTe*YdadjYdYd - 6._dp*(YdadjYdYdadjYdTd) +           & 
&  (8*g1p2*YdadjYuTu)/5._dp - 2*SPlamxxClam*YdadjYuTu - 6*TrYuadjYu*YdadjYuTu -          & 
&  2*TrYvadjYv*YdadjYuTu - 4._dp*(YdadjYuTuadjYdYd) - 4._dp*(YdadjYuTuadjYuYu) -         & 
&  (8*g1p2*M1*YdadjYuYu)/5._dp - 2*SPClamxxTlam*YdadjYuYu - 6*TradjYuTu*YdadjYuYu -      & 
&  2*TradjYvTv*YdadjYuYu - 2._dp*(YdadjYuYuadjYdTd) - 4._dp*(YdadjYuYuadjYuTu) +         & 
&  (287*g1p4*Td)/90._dp + g1p2*g2p2*Td + (15*g2p4*Td)/2._dp + (8*g1p2*g3p2*Td)/9._dp  
betaTd2 =  betaTd2+ 8*g2p2*g3p2*Td - (16*g3p4*Td)/9._dp - 2*SPlamxxadjYvYvClam*Td - 3*SPlamxxClamp2*Td -  & 
&  (2*g1p2*TrYdadjYd*Td)/5._dp + 16*g3p2*TrYdadjYd*Td - 9*TrYdadjYdYdadjYd*Td -          & 
&  3*TrYdadjYuYuadjYd*Td + (6*g1p2*TrYeadjYe*Td)/5._dp - 3*TrYeadjYeYeadjYe*Td -         & 
&  3*SPlamxxClam*TrYuadjYu*Td - SPlamxxClam*TrYvadjYv*Td - TrYvadjYvTpYeCYe*Td -         & 
&  2*Conjg(lam(1))*Td*Tpkap1Ckap1lam(1) - 4*Yd*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) -        & 
&  2*Conjg(lam(1))*Td*Tpkap1Ckap2lam(2) - 4*Yd*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) -        & 
&  2*Conjg(lam(1))*Td*Tpkap1Ckap3lam(3) - 4*Yd*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) -        & 
&  2*Conjg(lam(2))*Td*Tpkap2Ckap1lam(1) - 4*Yd*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) -        & 
&  2*Conjg(lam(2))*Td*Tpkap2Ckap2lam(2) - 4*Yd*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) -        & 
&  2*Conjg(lam(2))*Td*Tpkap2Ckap3lam(3) - 4*Yd*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) -        & 
&  2*Conjg(lam(3))*Td*Tpkap3Ckap1lam(1) - 4*Yd*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) -        & 
&  2*Conjg(lam(3))*Td*Tpkap3Ckap2lam(2) - 4*Yd*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) -        & 
&  2*Conjg(lam(3))*Td*Tpkap3Ckap3lam(3) - 4*Yd*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) -        & 
&  4*Yd*Conjg(lam(1))*TpTk1Ckap1lam(1) - 4*Yd*Conjg(lam(1))*TpTk1Ckap2lam(2) -           & 
&  4*Yd*Conjg(lam(1))*TpTk1Ckap3lam(3) - 4*Yd*Conjg(lam(2))*TpTk2Ckap1lam(1) -           & 
&  4*Yd*Conjg(lam(2))*TpTk2Ckap2lam(2) - 4*Yd*Conjg(lam(2))*TpTk2Ckap3lam(3) -           & 
&  4*Yd*Conjg(lam(3))*TpTk3Ckap1lam(1) - 4*Yd*Conjg(lam(3))*TpTk3Ckap2lam(2) -           & 
&  4*Yd*Conjg(lam(3))*TpTk3Ckap3lam(3)

 
DTd = oo16pi2*( betaTd1 + oo16pi2 * betaTd2 ) 

 
Else 
DTd = oo16pi2* betaTd1 
End If 
 
 
Call Chop(DTd) 

!-------------------- 
! Te 
!-------------------- 
 
betaTe1  = 5._dp*(TeadjYeYe) + TeCYvTpYv + (18*g1p2*M1*Ye)/5._dp + 6*g2p2*M2*Ye +     & 
&  2*SPClamxxTlam*Ye + 6*TradjYdTd*Ye + 2*TradjYeTe*Ye + 4._dp*(YeadjYeTe)               & 
&  + 2._dp*(YeCYvTpTv) - (9*g1p2*Te)/5._dp - 3*g2p2*Te + SPlamxxClam*Te + 3*TrYdadjYd*Te +& 
&  TrYeadjYe*Te

 
 
If (TwoLoopRGE) Then 
betaTe2 = -2._dp*(DyTeCYvi11Yvadjkap1kap1i21) - 2._dp*(DyTeCYvi11Yvadjkap2kap1i22) -            & 
&  2._dp*(DyTeCYvi11Yvadjkap3kap1i23) - 2._dp*(DyTeCYvi12Yvadjkap1kap2i21) -             & 
&  2._dp*(DyTeCYvi12Yvadjkap2kap2i22) - 2._dp*(DyTeCYvi12Yvadjkap3kap2i23) -             & 
&  2._dp*(DyTeCYvi13Yvadjkap1kap3i21) - 2._dp*(DyTeCYvi13Yvadjkap2kap3i22) -             & 
&  2._dp*(DyTeCYvi13Yvadjkap3kap3i23) - 2._dp*(DyTeCYvlami1YvClami2) - 4._dp*(DyYeCYvi11TvCkap1kap1i21) -& 
&  4._dp*(DyYeCYvi11TvCkap2kap1i22) - 4._dp*(DyYeCYvi11TvCkap3kap1i23) - 4._dp*(DyYeCYvi12TvCkap1kap2i21) -& 
&  4._dp*(DyYeCYvi12TvCkap2kap2i22) - 4._dp*(DyYeCYvi12TvCkap3kap2i23) - 4._dp*(DyYeCYvi13TvCkap1kap3i21) -& 
&  4._dp*(DyYeCYvi13TvCkap2kap3i22) - 4._dp*(DyYeCYvi13TvCkap3kap3i23) - 4._dp*(DyYeCYvlami1TvClami2) -& 
&  4._dp*(DyYeCYvTlami1YvClami2) - (6*g1p2*TeadjYeYe)/5._dp + 12*g2p2*TeadjYeYe -        & 
&  5*SPlamxxClam*TeadjYeYe - 6._dp*(TeadjYeYeadjYeYe) - SPlamxxClam*TeCYvTpYv -          & 
&  4._dp*(TeCYvTpYvadjYeYe) - 2._dp*(TeCYvTpYvCYvTpYv) - 15*TeadjYeYe*TrYdadjYd -        & 
&  5*TeadjYeYe*TrYeadjYe - 3*TeCYvTpYv*TrYuadjYu - TeCYvTpYv*TrYvadjYv - 54*g1p4*M1*Ye - & 
&  (18*g1p2*g2p2*M1*Ye)/5._dp - (18*g1p2*g2p2*M2*Ye)/5._dp - 30*g2p4*M2*Ye -             & 
&  4*SPadjYvYvClamxxTlam*Ye - 4*SPClamxxTpTvCYvlam*Ye - 12*SPClamxxTlam*SPlamxxClam*Ye - & 
&  (4*g1p2*TradjYdTd*Ye)/5._dp + 32*g3p2*TradjYdTd*Ye + (12*g1p2*TradjYeTe*Ye)/5._dp -   & 
&  2*TradjYeTeCYvTpYv*Ye - 6*SPlamxxClam*TradjYuTu*Ye - 2*TradjYvTpYeCYeTv*Ye -          & 
&  2*SPlamxxClam*TradjYvTv*Ye + (4*g1p2*M1*TrYdadjYd*Ye)/5._dp - 32*g3p2*M3*TrYdadjYd*Ye -& 
&  36*TrYdadjYdTdadjYd*Ye - 6*TrYdadjYuTuadjYd*Ye - (12*g1p2*M1*TrYeadjYe*Ye)/5._dp -    & 
&  12*TrYeadjYeTeadjYe*Ye - 6*TrYuadjYdTdadjYu*Ye - 6*SPClamxxTlam*TrYuadjYu*Ye -        & 
&  2*SPClamxxTlam*TrYvadjYv*Ye + (6*g1p2*YeadjYeTe)/5._dp + 6*g2p2*YeadjYeTe -           & 
&  4*SPlamxxClam*YeadjYeTe - 12*TrYdadjYd*YeadjYeTe - 4*TrYeadjYe*YeadjYeTe  
betaTe2 =  betaTe2- 8._dp*(YeadjYeTeadjYeYe) - 12*g2p2*M2*YeadjYeYe - 6*SPClamxxTlam*YeadjYeYe -          & 
&  18*TradjYdTd*YeadjYeYe - 6*TradjYeTe*YeadjYeYe - 6._dp*(YeadjYeYeadjYeTe) -           & 
&  4._dp*(YeCYvTk1adjkap1TpYv) - 4._dp*(YeCYvTk2adjkap2TpYv) - 4._dp*(YeCYvTk3adjkap3TpYv) -& 
&  2*SPlamxxClam*YeCYvTpTv - 6*TrYuadjYu*YeCYvTpTv - 2*TrYvadjYv*YeCYvTpTv -             & 
&  4._dp*(YeCYvTpTvadjYeYe) - 4._dp*(YeCYvTpTvCYvTpYv) - 2*SPClamxxTlam*YeCYvTpYv -      & 
&  6*TradjYuTu*YeCYvTpYv - 2*TradjYvTv*YeCYvTpYv - 2._dp*(YeCYvTpYvadjYeTe) -            & 
&  4._dp*(YeCYvTpYvCYvTpTv) + (27*g1p4*Te)/2._dp + (9*g1p2*g2p2*Te)/5._dp +              & 
&  (15*g2p4*Te)/2._dp - 2*SPlamxxadjYvYvClam*Te - 3*SPlamxxClamp2*Te - (2*g1p2*TrYdadjYd*Te)/5._dp +& 
&  16*g3p2*TrYdadjYd*Te - 9*TrYdadjYdYdadjYd*Te - 3*TrYdadjYuYuadjYd*Te + (6*g1p2*TrYeadjYe*Te)/5._dp -& 
&  3*TrYeadjYeYeadjYe*Te - 3*SPlamxxClam*TrYuadjYu*Te - SPlamxxClam*TrYvadjYv*Te -       & 
&  TrYvadjYvTpYeCYe*Te - 2*Conjg(lam(1))*Te*Tpkap1Ckap1lam(1) - 4*Ye*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) -& 
&  2*Conjg(lam(1))*Te*Tpkap1Ckap2lam(2) - 4*Ye*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) -        & 
&  2*Conjg(lam(1))*Te*Tpkap1Ckap3lam(3) - 4*Ye*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) -        & 
&  2*Conjg(lam(2))*Te*Tpkap2Ckap1lam(1) - 4*Ye*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) -        & 
&  2*Conjg(lam(2))*Te*Tpkap2Ckap2lam(2) - 4*Ye*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) -        & 
&  2*Conjg(lam(2))*Te*Tpkap2Ckap3lam(3) - 4*Ye*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) -        & 
&  2*Conjg(lam(3))*Te*Tpkap3Ckap1lam(1) - 4*Ye*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) -        & 
&  2*Conjg(lam(3))*Te*Tpkap3Ckap2lam(2) - 4*Ye*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) -        & 
&  2*Conjg(lam(3))*Te*Tpkap3Ckap3lam(3) - 4*Ye*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) -        & 
&  4*Ye*Conjg(lam(1))*TpTk1Ckap1lam(1) - 4*Ye*Conjg(lam(1))*TpTk1Ckap2lam(2) -           & 
&  4*Ye*Conjg(lam(1))*TpTk1Ckap3lam(3) - 4*Ye*Conjg(lam(2))*TpTk2Ckap1lam(1)  
betaTe2 =  betaTe2- 4*Ye*Conjg(lam(2))*TpTk2Ckap2lam(2) - 4*Ye*Conjg(lam(2))*TpTk2Ckap3lam(3) -           & 
&  4*Ye*Conjg(lam(3))*TpTk3Ckap1lam(1) - 4*Ye*Conjg(lam(3))*TpTk3Ckap2lam(2) -           & 
&  4*Ye*Conjg(lam(3))*TpTk3Ckap3lam(3)

 
DTe = oo16pi2*( betaTe1 + oo16pi2 * betaTe2 ) 

 
Else 
DTe = oo16pi2* betaTe1 
End If 
 
 
Call Chop(DTe) 

!-------------------- 
! Tlam 
!-------------------- 
 
Do i1 = 1,3
betaTlam1(i1) = 4*Tk1adjkap1lam(i1) + 4*Tk2adjkap2lam(i1) + 4*Tk3adjkap3lam(i1) + 5*TpTvCYvlam(i1) +  & 
&  4*TpYvCYvTlam(i1) + (6*g1p2*M1*lam(i1))/5._dp + 6*g2p2*M2*lam(i1) + 6*SPClamxxTlam*lam(i1) +& 
&  6*TradjYdTd*lam(i1) + 2*TradjYeTe*lam(i1) + 6*TradjYuTu*lam(i1) + 2*TradjYvTv*lam(i1) +& 
&  2*kap1Ckap1(i1,1)*Tlam(1) + 2*kap2Ckap1(i1,2)*Tlam(1) + 2*kap3Ckap1(i1,               & 
& 3)*Tlam(1) + 2*kap1Ckap2(i1,1)*Tlam(2) + 2*kap2Ckap2(i1,2)*Tlam(2) + 2*kap3Ckap2(i1,   & 
& 3)*Tlam(2) + 2*kap1Ckap3(i1,1)*Tlam(3) + 2*kap2Ckap3(i1,2)*Tlam(3) + 2*kap3Ckap3(i1,   & 
& 3)*Tlam(3) - (3*g1p2*Tlam(i1))/5._dp - 3*g2p2*Tlam(i1) + 6*SPlamxxClam*Tlam(i1) +      & 
&  3*TrYdadjYd*Tlam(i1) + TrYeadjYe*Tlam(i1) + 3*TrYuadjYu*Tlam(i1) + TrYvadjYv*Tlam(i1)
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
betaTlam2(i1) = -16*TradjYvTv*kap1Ckap1lam(i1) - 16*TrCkap1Tk1*kap1Ckap1lam(i1) - 16*Ckap2Tk1(1,      & 
& 2)*kap1Ckap1lam(i1) - 16*Ckap3Tk1(1,3)*kap1Ckap1lam(i1) - 16*TrCkap1Tk2*kap1Ckap2lam(i1) -& 
&  16*adjYvTv(1,2)*kap1Ckap2lam(i1) - 16*Ckap2Tk2(1,2)*kap1Ckap2lam(i1) - 16*Ckap3Tk2(1, & 
& 3)*kap1Ckap2lam(i1) - 16*TrCkap1Tk3*kap1Ckap3lam(i1) - 16*adjYvTv(1,3)*kap1Ckap3lam(i1) -& 
&  16*Ckap2Tk3(1,2)*kap1Ckap3lam(i1) - 16*Ckap3Tk3(1,3)*kap1Ckap3lam(i1) -               & 
&  16*TrCkap2Tk1*kap2Ckap1lam(i1) - 16*adjYvTv(2,1)*kap2Ckap1lam(i1) - 16*Ckap1Tk1(2,    & 
& 1)*kap2Ckap1lam(i1) - 16*Ckap3Tk1(2,3)*kap2Ckap1lam(i1) - 16*TradjYvTv*kap2Ckap2lam(i1) -& 
&  16*TrCkap2Tk2*kap2Ckap2lam(i1) - 16*Ckap1Tk2(2,1)*kap2Ckap2lam(i1) - 16*Ckap3Tk2(2,   & 
& 3)*kap2Ckap2lam(i1) - 16*TrCkap2Tk3*kap2Ckap3lam(i1) - 16*adjYvTv(2,3)*kap2Ckap3lam(i1) -& 
&  16*Ckap1Tk3(2,1)*kap2Ckap3lam(i1) - 16*Ckap3Tk3(2,3)*kap2Ckap3lam(i1) -               & 
&  16*TrCkap3Tk1*kap3Ckap1lam(i1) - 16*adjYvTv(3,1)*kap3Ckap1lam(i1) - 16*Ckap1Tk1(3,    & 
& 1)*kap3Ckap1lam(i1) - 16*Ckap2Tk1(3,2)*kap3Ckap1lam(i1) - 16*TrCkap3Tk2*kap3Ckap2lam(i1) -& 
&  16*adjYvTv(3,2)*kap3Ckap2lam(i1) - 16*Ckap1Tk2(3,1)*kap3Ckap2lam(i1) - 16*Ckap2Tk2(3, & 
& 2)*kap3Ckap2lam(i1) - 16*TradjYvTv*kap3Ckap3lam(i1) - 16*TrCkap3Tk3*kap3Ckap3lam(i1) - & 
&  16*Ckap1Tk3(3,1)*kap3Ckap3lam(i1) - 16*Ckap2Tk3(3,2)*kap3Ckap3lam(i1) -               & 
&  16*Tk1adjYvYvCkap1lam(i1) - 16*SPlamxxadjkap1lam*Tk1Clam(i1) - 16*Tk2adjYvYvCkap2lam(i1) -& 
&  16*SPlamxxadjkap2lam*Tk2Clam(i1) - 16*Tk3adjYvYvCkap3lam(i1) - 16*SPlamxxadjkap3lam*Tk3Clam(i1) -& 
&  2*adjYvTv(1,i1)*Tpkap1Ckap1lam(1) - 16*Ckap1lam(1)*Tpkap1Ckap1TpTk1(1,i1) -           & 
&  16*Ckap2lam(1)*Tpkap1Ckap1TpTk2(1,i1) - 16*Ckap3lam(1)*Tpkap1Ckap1TpTk3(1,            & 
& i1) - 2*adjYvTv(1,i1)*Tpkap1Ckap2lam(2) - 16*Ckap1lam(1)*Tpkap1Ckap2TpTk1(2,           & 
& i1) - 16*Ckap2lam(1)*Tpkap1Ckap2TpTk2(2,i1) - 16*Ckap3lam(1)*Tpkap1Ckap2TpTk3(2,       & 
& i1) - 2*adjYvTv(1,i1)*Tpkap1Ckap3lam(3) - 16*Ckap1lam(1)*Tpkap1Ckap3TpTk1(3,           & 
& i1) - 16*Ckap2lam(1)*Tpkap1Ckap3TpTk2(3,i1) - 16*Ckap3lam(1)*Tpkap1Ckap3TpTk3(3,       & 
& i1) - 2*adjYvTv(2,i1)*Tpkap2Ckap1lam(1) - 16*Ckap1lam(2)*Tpkap2Ckap1TpTk1(1,           & 
& i1) - 16*Ckap2lam(2)*Tpkap2Ckap1TpTk2(1,i1) - 16*Ckap3lam(2)*Tpkap2Ckap1TpTk3(1,       & 
& i1) - 2*adjYvTv(2,i1)*Tpkap2Ckap2lam(2) - 16*Ckap1lam(2)*Tpkap2Ckap2TpTk1(2,           & 
& i1) - 16*Ckap2lam(2)*Tpkap2Ckap2TpTk2(2,i1) - 16*Ckap3lam(2)*Tpkap2Ckap2TpTk3(2,       & 
& i1) - 2*adjYvTv(2,i1)*Tpkap2Ckap3lam(3) - 16*Ckap1lam(2)*Tpkap2Ckap3TpTk1(3,           & 
& i1) - 16*Ckap2lam(2)*Tpkap2Ckap3TpTk2(3,i1) - 16*Ckap3lam(2)*Tpkap2Ckap3TpTk3(3,       & 
& i1) - 2*adjYvTv(3,i1)*Tpkap3Ckap1lam(1) - 16*Ckap1lam(3)*Tpkap3Ckap1TpTk1(1,           & 
& i1) - 16*Ckap2lam(3)*Tpkap3Ckap1TpTk2(1,i1) - 16*Ckap3lam(3)*Tpkap3Ckap1TpTk3(1,       & 
& i1) - 2*adjYvTv(3,i1)*Tpkap3Ckap2lam(2) - 16*Ckap1lam(3)*Tpkap3Ckap2TpTk1(2,           & 
& i1) - 16*Ckap2lam(3)*Tpkap3Ckap2TpTk2(2,i1) - 16*Ckap3lam(3)*Tpkap3Ckap2TpTk3(2,       & 
& i1) - 2*adjYvTv(3,i1)*Tpkap3Ckap3lam(3) - 16*Ckap1lam(3)*Tpkap3Ckap3TpTk1(3,           & 
& i1) - 16*Ckap2lam(3)*Tpkap3Ckap3TpTk2(3,i1) - 16*Ckap3lam(3)*Tpkap3Ckap3TpTk3(3,       & 
& i1) - 3*TpTvadjYeYeCYvlam(i1) + (12*g1p2*TpTvCYvlam(i1))/5._dp + 12*g2p2*TpTvCYvlam(i1)  
betaTlam2(i1) =  betaTlam2(i1)- 11*SPlamxxClam*TpTvCYvlam(i1) - 15*TrYuadjYu*TpTvCYvlam(i1) - 5*TrYvadjYv*TpTvCYvlam(i1) -& 
&  6*TpTvCYvTpYvCYvlam(i1) - 2*TpYvadjYeTeCYvlam(i1) - 4*Tpkap1adjkap1Tlam(1)*TpYvCYv(i1,& 
& 1) - 4*Tpkap1adjkap2Tlam(2)*TpYvCYv(i1,1) - 4*Tpkap1adjkap3Tlam(3)*TpYvCYv(i1,         & 
& 1) - 4*Tpkap2adjkap1Tlam(1)*TpYvCYv(i1,2) - 4*Tpkap2adjkap2Tlam(2)*TpYvCYv(i1,         & 
& 2) - 4*Tpkap2adjkap3Tlam(3)*TpYvCYv(i1,2) - 4*Tpkap3adjkap1Tlam(1)*TpYvCYv(i1,         & 
& 3) - 4*Tpkap3adjkap2Tlam(2)*TpYvCYv(i1,3) - 4*Tpkap3adjkap3Tlam(3)*TpYvCYv(i1,         & 
& 3) - (12*g1p2*M1*TpYvCYvlam(i1))/5._dp - 12*g2p2*M2*TpYvCYvlam(i1) - 12*SPClamxxTlam*TpYvCYvlam(i1) -& 
&  18*TradjYuTu*TpYvCYvlam(i1) - 6*TradjYvTv*TpYvCYvlam(i1) - 4*TpYvCYvTk1adjkap1lam(i1) -& 
&  4*TpYvCYvTk2adjkap2lam(i1) - 4*TpYvCYvTk3adjkap3lam(i1) + (6*g1p2*TpYvCYvTlam(i1))/5._dp +& 
&  6*g2p2*TpYvCYvTlam(i1) - 12*SPlamxxClam*TpYvCYvTlam(i1) - 12*TrYuadjYu*TpYvCYvTlam(i1) -& 
&  4*TrYvadjYv*TpYvCYvTlam(i1) - 8*TpYvCYvTpTvCYvlam(i1) - 6*TpYvCYvTpYvCYvTlam(i1) -    & 
&  (414*g1p4*M1*lam(i1))/25._dp - (18*g1p2*g2p2*M1*lam(i1))/5._dp - (18*g1p2*g2p2*M2*lam(i1))/5._dp -& 
&  30*g2p4*M2*lam(i1) - 18*SPadjYvYvClamxxTlam*lam(i1) + (6*g1p2*SPClamxxTlam*lam(i1))/5._dp +& 
&  6*g2p2*SPClamxxTlam*lam(i1) - 20*SPClamxxTpTvCYvlam*lam(i1) - (12*g1p2*M1*SPlamxxClam*lam(i1))/5._dp -& 
&  12*g2p2*M2*SPlamxxClam*lam(i1) - 36*SPClamxxTlam*SPlamxxClam*lam(i1) - (4*g1p2*TradjYdTd*lam(i1))/5._dp +& 
&  32*g3p2*TradjYdTd*lam(i1) - 18*SPlamxxClam*TradjYdTd*lam(i1) + (12*g1p2*TradjYeTe*lam(i1))/5._dp -& 
&  6*SPlamxxClam*TradjYeTe*lam(i1) - 4*TradjYeTeCYvTpYv*lam(i1) + (8*g1p2*TradjYuTu*lam(i1))/5._dp +& 
&  32*g3p2*TradjYuTu*lam(i1) - 18*SPlamxxClam*TradjYuTu*lam(i1) - 4*TradjYvTpYeCYeTv*lam(i1) -& 
&  6*SPlamxxClam*TradjYvTv*lam(i1) - 4*TradjYvTvadjkap1kap1*lam(i1) - 4*TradjYvTvadjkap2kap2*lam(i1) -& 
&  4*TradjYvTvadjkap3kap3*lam(i1) + (4*g1p2*M1*TrYdadjYd*lam(i1))/5._dp - 32*g3p2*M3*TrYdadjYd*lam(i1) -& 
&  12*SPClamxxTlam*TrYdadjYd*lam(i1) - 36*TrYdadjYdTdadjYd*lam(i1) - 12*TrYdadjYuTuadjYd*lam(i1)  
betaTlam2(i1) =  betaTlam2(i1)- (12*g1p2*M1*TrYeadjYe*lam(i1))/5._dp - 4*SPClamxxTlam*TrYeadjYe*lam(i1) -             & 
&  12*TrYeadjYeTeadjYe*lam(i1) - 12*TrYuadjYdTdadjYu*lam(i1) - (8*g1p2*M1*TrYuadjYu*lam(i1))/5._dp -& 
&  32*g3p2*M3*TrYuadjYu*lam(i1) - 12*SPClamxxTlam*TrYuadjYu*lam(i1) - 36*TrYuadjYuTuadjYu*lam(i1) -& 
&  4*SPClamxxTlam*TrYvadjYv*lam(i1) - 12*TrYvadjYvTvadjYv*lam(i1) - 4*TrYvCkap1TpTk1adjYv*lam(i1) -& 
&  4*TrYvCkap2TpTk2adjYv*lam(i1) - 4*TrYvCkap3TpTk3adjYv*lam(i1) - 4*adjYvTvadjkap1kap2(2,& 
& 1)*lam(i1) - 4*adjYvTvadjkap1kap3(3,1)*lam(i1) - 4*adjYvTvadjkap2kap1(1,               & 
& 2)*lam(i1) - 4*adjYvTvadjkap2kap3(3,2)*lam(i1) - 4*adjYvTvadjkap3kap1(1,               & 
& 3)*lam(i1) - 4*adjYvTvadjkap3kap2(2,3)*lam(i1) - 8*Conjg(lam(1))*Tpkap1Ckap1Tlam(1)*lam(i1) -& 
&  8*Conjg(lam(1))*Tpkap1Ckap2Tlam(2)*lam(i1) - 8*Conjg(lam(1))*Tpkap1Ckap3Tlam(3)*lam(i1) -& 
&  8*Conjg(lam(2))*Tpkap2Ckap1Tlam(1)*lam(i1) - 8*Conjg(lam(2))*Tpkap2Ckap2Tlam(2)*lam(i1) -& 
&  8*Conjg(lam(2))*Tpkap2Ckap3Tlam(3)*lam(i1) - 8*Conjg(lam(3))*Tpkap3Ckap1Tlam(1)*lam(i1) -& 
&  8*Conjg(lam(3))*Tpkap3Ckap2Tlam(2)*lam(i1) - 8*Conjg(lam(3))*Tpkap3Ckap3Tlam(3)*lam(i1) -& 
&  8*Conjg(lam(1))*TpTk1Ckap1lam(1)*lam(i1) - 8*Conjg(lam(1))*TpTk1Ckap2lam(2)*lam(i1) - & 
&  8*Conjg(lam(1))*TpTk1Ckap3lam(3)*lam(i1) - 8*Conjg(lam(2))*TpTk2Ckap1lam(1)*lam(i1) - & 
&  8*Conjg(lam(2))*TpTk2Ckap2lam(2)*lam(i1) - 8*Conjg(lam(2))*TpTk2Ckap3lam(3)*lam(i1) - & 
&  8*Conjg(lam(3))*TpTk3Ckap1lam(1)*lam(i1) - 8*Conjg(lam(3))*TpTk3Ckap2lam(2)*lam(i1) - & 
&  8*Conjg(lam(3))*TpTk3Ckap3lam(3)*lam(i1) - 8*kap1adjYvYvCkap1(i1,1)*Tlam(1) -         & 
&  8*TrCkap1kap1*kap1Ckap1(i1,1)*Tlam(1) - 8*Ckap2kap1(1,2)*kap1Ckap1(i1,1)*Tlam(1) -    & 
&  8*Ckap3kap1(1,3)*kap1Ckap1(i1,1)*Tlam(1) - 8*TrCkap1kap2*kap1Ckap1(i1,2)*Tlam(1) -    & 
&  8*Ckap2kap2(1,2)*kap1Ckap1(i1,2)*Tlam(1) - 8*Ckap3kap2(1,3)*kap1Ckap1(i1,             & 
& 2)*Tlam(1) - 8*TrCkap1kap3*kap1Ckap1(i1,3)*Tlam(1) - 8*Ckap2kap3(1,2)*kap1Ckap1(i1,    & 
& 3)*Tlam(1) - 8*Ckap3kap3(1,3)*kap1Ckap1(i1,3)*Tlam(1) - 24*Conjg(lam(1))*kap1Ckap1lam(i1)*Tlam(1)  
betaTlam2(i1) =  betaTlam2(i1)- 8*kap2adjYvYvCkap1(i1,2)*Tlam(1) - 8*TrCkap2kap1*kap2Ckap1(i1,1)*Tlam(1) -            & 
&  8*Ckap1kap1(2,1)*kap2Ckap1(i1,1)*Tlam(1) - 8*Ckap3kap1(2,3)*kap2Ckap1(i1,             & 
& 1)*Tlam(1) - 8*TrCkap2kap2*kap2Ckap1(i1,2)*Tlam(1) - 8*Ckap1kap2(2,1)*kap2Ckap1(i1,    & 
& 2)*Tlam(1) - 8*Ckap3kap2(2,3)*kap2Ckap1(i1,2)*Tlam(1) - 8*TrCkap2kap3*kap2Ckap1(i1,    & 
& 3)*Tlam(1) - 8*Ckap1kap3(2,1)*kap2Ckap1(i1,3)*Tlam(1) - 8*Ckap3kap3(2,3)*kap2Ckap1(i1, & 
& 3)*Tlam(1) - 24*Conjg(lam(2))*kap2Ckap1lam(i1)*Tlam(1) - 8*kap3adjYvYvCkap1(i1,        & 
& 3)*Tlam(1) - 8*TrCkap3kap1*kap3Ckap1(i1,1)*Tlam(1) - 8*Ckap1kap1(3,1)*kap3Ckap1(i1,    & 
& 1)*Tlam(1) - 8*Ckap2kap1(3,2)*kap3Ckap1(i1,1)*Tlam(1) - 8*TrCkap3kap2*kap3Ckap1(i1,    & 
& 2)*Tlam(1) - 8*Ckap1kap2(3,1)*kap3Ckap1(i1,2)*Tlam(1) - 8*Ckap2kap2(3,2)*kap3Ckap1(i1, & 
& 2)*Tlam(1) - 8*TrCkap3kap3*kap3Ckap1(i1,3)*Tlam(1) - 8*Ckap1kap3(3,1)*kap3Ckap1(i1,    & 
& 3)*Tlam(1) - 8*Ckap2kap3(3,2)*kap3Ckap1(i1,3)*Tlam(1) - 24*Conjg(lam(3))*kap3Ckap1lam(i1)*Tlam(1) -& 
&  8*kap1adjYvYvCkap2(i1,1)*Tlam(2) - 8*TrCkap1kap1*kap1Ckap2(i1,1)*Tlam(2) -            & 
&  8*Ckap2kap1(1,2)*kap1Ckap2(i1,1)*Tlam(2) - 8*Ckap3kap1(1,3)*kap1Ckap2(i1,             & 
& 1)*Tlam(2) - 8*TrCkap1kap2*kap1Ckap2(i1,2)*Tlam(2) - 8*Ckap2kap2(1,2)*kap1Ckap2(i1,    & 
& 2)*Tlam(2) - 8*Ckap3kap2(1,3)*kap1Ckap2(i1,2)*Tlam(2) - 8*TrCkap1kap3*kap1Ckap2(i1,    & 
& 3)*Tlam(2) - 8*Ckap2kap3(1,2)*kap1Ckap2(i1,3)*Tlam(2) - 8*Ckap3kap3(1,3)*kap1Ckap2(i1, & 
& 3)*Tlam(2) - 24*Conjg(lam(1))*kap1Ckap2lam(i1)*Tlam(2) - 8*kap2adjYvYvCkap2(i1,        & 
& 2)*Tlam(2) - 8*TrCkap2kap1*kap2Ckap2(i1,1)*Tlam(2) - 8*Ckap1kap1(2,1)*kap2Ckap2(i1,    & 
& 1)*Tlam(2) - 8*Ckap3kap1(2,3)*kap2Ckap2(i1,1)*Tlam(2) - 8*TrCkap2kap2*kap2Ckap2(i1,    & 
& 2)*Tlam(2) - 8*Ckap1kap2(2,1)*kap2Ckap2(i1,2)*Tlam(2) - 8*Ckap3kap2(2,3)*kap2Ckap2(i1, & 
& 2)*Tlam(2) - 8*TrCkap2kap3*kap2Ckap2(i1,3)*Tlam(2) - 8*Ckap1kap3(2,1)*kap2Ckap2(i1,    & 
& 3)*Tlam(2) - 8*Ckap3kap3(2,3)*kap2Ckap2(i1,3)*Tlam(2) - 24*Conjg(lam(2))*kap2Ckap2lam(i1)*Tlam(2)  
betaTlam2(i1) =  betaTlam2(i1)- 8*kap3adjYvYvCkap2(i1,3)*Tlam(2) - 8*TrCkap3kap1*kap3Ckap2(i1,1)*Tlam(2) -            & 
&  8*Ckap1kap1(3,1)*kap3Ckap2(i1,1)*Tlam(2) - 8*Ckap2kap1(3,2)*kap3Ckap2(i1,             & 
& 1)*Tlam(2) - 8*TrCkap3kap2*kap3Ckap2(i1,2)*Tlam(2) - 8*Ckap1kap2(3,1)*kap3Ckap2(i1,    & 
& 2)*Tlam(2) - 8*Ckap2kap2(3,2)*kap3Ckap2(i1,2)*Tlam(2) - 8*TrCkap3kap3*kap3Ckap2(i1,    & 
& 3)*Tlam(2) - 8*Ckap1kap3(3,1)*kap3Ckap2(i1,3)*Tlam(2) - 8*Ckap2kap3(3,2)*kap3Ckap2(i1, & 
& 3)*Tlam(2) - 24*Conjg(lam(3))*kap3Ckap2lam(i1)*Tlam(2) - 8*kap1adjYvYvCkap3(i1,        & 
& 1)*Tlam(3) - 8*TrCkap1kap1*kap1Ckap3(i1,1)*Tlam(3) - 8*Ckap2kap1(1,2)*kap1Ckap3(i1,    & 
& 1)*Tlam(3) - 8*Ckap3kap1(1,3)*kap1Ckap3(i1,1)*Tlam(3) - 8*TrCkap1kap2*kap1Ckap3(i1,    & 
& 2)*Tlam(3) - 8*Ckap2kap2(1,2)*kap1Ckap3(i1,2)*Tlam(3) - 8*Ckap3kap2(1,3)*kap1Ckap3(i1, & 
& 2)*Tlam(3) - 8*TrCkap1kap3*kap1Ckap3(i1,3)*Tlam(3) - 8*Ckap2kap3(1,2)*kap1Ckap3(i1,    & 
& 3)*Tlam(3) - 8*Ckap3kap3(1,3)*kap1Ckap3(i1,3)*Tlam(3) - 24*Conjg(lam(1))*kap1Ckap3lam(i1)*Tlam(3) -& 
&  8*kap2adjYvYvCkap3(i1,2)*Tlam(3) - 8*TrCkap2kap1*kap2Ckap3(i1,1)*Tlam(3) -            & 
&  8*Ckap1kap1(2,1)*kap2Ckap3(i1,1)*Tlam(3) - 8*Ckap3kap1(2,3)*kap2Ckap3(i1,             & 
& 1)*Tlam(3) - 8*TrCkap2kap2*kap2Ckap3(i1,2)*Tlam(3) - 8*Ckap1kap2(2,1)*kap2Ckap3(i1,    & 
& 2)*Tlam(3) - 8*Ckap3kap2(2,3)*kap2Ckap3(i1,2)*Tlam(3) - 8*TrCkap2kap3*kap2Ckap3(i1,    & 
& 3)*Tlam(3) - 8*Ckap1kap3(2,1)*kap2Ckap3(i1,3)*Tlam(3) - 8*Ckap3kap3(2,3)*kap2Ckap3(i1, & 
& 3)*Tlam(3) - 24*Conjg(lam(2))*kap2Ckap3lam(i1)*Tlam(3) - 8*kap3adjYvYvCkap3(i1,        & 
& 3)*Tlam(3) - 8*TrCkap3kap1*kap3Ckap3(i1,1)*Tlam(3) - 8*Ckap1kap1(3,1)*kap3Ckap3(i1,    & 
& 1)*Tlam(3) - 8*Ckap2kap1(3,2)*kap3Ckap3(i1,1)*Tlam(3) - 8*TrCkap3kap2*kap3Ckap3(i1,    & 
& 2)*Tlam(3) - 8*Ckap1kap2(3,1)*kap3Ckap3(i1,2)*Tlam(3) - 8*Ckap2kap2(3,2)*kap3Ckap3(i1, & 
& 2)*Tlam(3) - 8*TrCkap3kap3*kap3Ckap3(i1,3)*Tlam(3) - 8*Ckap1kap3(3,1)*kap3Ckap3(i1,    & 
& 3)*Tlam(3) - 8*Ckap2kap3(3,2)*kap3Ckap3(i1,3)*Tlam(3) - 24*Conjg(lam(3))*kap3Ckap3lam(i1)*Tlam(3)  
betaTlam2(i1) =  betaTlam2(i1)+ (207*g1p4*Tlam(i1))/50._dp + (9*g1p2*g2p2*Tlam(i1))/5._dp + (15*g2p4*Tlam(i1))/2._dp -& 
&  12*SPlamxxadjYvYvClam*Tlam(i1) + (12*g1p2*SPlamxxClam*Tlam(i1))/5._dp +               & 
&  12*g2p2*SPlamxxClam*Tlam(i1) - 14*SPlamxxClamp2*Tlam(i1) - (2*g1p2*TrYdadjYd*Tlam(i1))/5._dp +& 
&  16*g3p2*TrYdadjYd*Tlam(i1) - 15*SPlamxxClam*TrYdadjYd*Tlam(i1) - 9*TrYdadjYdYdadjYd*Tlam(i1) -& 
&  6*TrYdadjYuYuadjYd*Tlam(i1) + (6*g1p2*TrYeadjYe*Tlam(i1))/5._dp - 5*SPlamxxClam*TrYeadjYe*Tlam(i1) -& 
&  3*TrYeadjYeYeadjYe*Tlam(i1) + (4*g1p2*TrYuadjYu*Tlam(i1))/5._dp + 16*g3p2*TrYuadjYu*Tlam(i1) -& 
&  15*SPlamxxClam*TrYuadjYu*Tlam(i1) - 9*TrYuadjYuYuadjYu*Tlam(i1) - 5*SPlamxxClam*TrYvadjYv*Tlam(i1) -& 
&  2*TrYvadjYvTpYeCYe*Tlam(i1) - 3*TrYvadjYvYvadjYv*Tlam(i1) - 2*TrYvCkap1kap1adjYv*Tlam(i1) -& 
&  2*TrYvCkap2kap2adjYv*Tlam(i1) - 2*TrYvCkap3kap3adjYv*Tlam(i1) - 2*adjYvYvCkap1kap2(2, & 
& 1)*Tlam(i1) - 2*adjYvYvCkap1kap3(3,1)*Tlam(i1) - 2*adjYvYvCkap2kap1(1,2)*Tlam(i1) -    & 
&  2*adjYvYvCkap2kap3(3,2)*Tlam(i1) - 2*adjYvYvCkap3kap1(1,3)*Tlam(i1) - 2*adjYvYvCkap3kap2(2,& 
& 3)*Tlam(i1) - 4*Conjg(lam(1))*Tpkap1Ckap1lam(1)*Tlam(i1) - 4*Conjg(lam(1))*Tpkap1Ckap2lam(2)*Tlam(i1) -& 
&  4*Conjg(lam(1))*Tpkap1Ckap3lam(3)*Tlam(i1) - 4*Conjg(lam(2))*Tpkap2Ckap1lam(1)*Tlam(i1) -& 
&  4*Conjg(lam(2))*Tpkap2Ckap2lam(2)*Tlam(i1) - 4*Conjg(lam(2))*Tpkap2Ckap3lam(3)*Tlam(i1) -& 
&  4*Conjg(lam(3))*Tpkap3Ckap1lam(1)*Tlam(i1) - 4*Conjg(lam(3))*Tpkap3Ckap2lam(2)*Tlam(i1) -& 
&  4*Conjg(lam(3))*Tpkap3Ckap3lam(3)*Tlam(i1)
End Do

 
DTlam = oo16pi2*( betaTlam1 + oo16pi2 * betaTlam2 ) 

 
Else 
DTlam = oo16pi2* betaTlam1 
End If 
 
 
Call Chop(DTlam) 

!-------------------- 
! Tv 
!-------------------- 
 
betaTv1  = 4._dp*(DyTvClami1lami2) + 2._dp*(DyTvi11kap1Ckap1i21) + 2._dp*(DyTvi11kap2Ckap1i22)& 
&  + 2._dp*(DyTvi11kap3Ckap1i23) + 2._dp*(DyTvi12kap1Ckap2i21) + 2._dp*(DyTvi12kap2Ckap2i22)& 
&  + 2._dp*(DyTvi12kap3Ckap2i23) + 2._dp*(DyTvi13kap1Ckap3i21) + 2._dp*(DyTvi13kap2Ckap3i22)& 
&  + 2._dp*(DyTvi13kap3Ckap3i23) + 5._dp*(DyYvClami1Tlami2) + 2._dp*(TpTeCYeYv)          & 
&  + TpYeCYeTv + 4._dp*(TvadjYvYv) + (6*g1p2*M1*Yv)/5._dp + 6*g2p2*M2*Yv +               & 
&  2*SPClamxxTlam*Yv + 6*TradjYuTu*Yv + 2*TradjYvTv*Yv + 5._dp*(YvadjYvTv)               & 
&  + 4._dp*(YvCkap1TpTk1) + 4._dp*(YvCkap2TpTk2) + 4._dp*(YvCkap3TpTk3) - (3*g1p2*Tv)    & 
& /5._dp - 3*g2p2*Tv + SPlamxxClam*Tv + 3*TrYuadjYu*Tv + TrYvadjYv*Tv

 
 
If (TwoLoopRGE) Then 
betaTv2 = 2._dp*(DyTpTeCYeYvClami1lami2) + 2._dp*(DyTpYeCYeTvClami1lami2) + DyTpYeCYeYvClami1Tlami2 -& 
&  6._dp*(DyTvadjYvYvClami1lami2) - 4._dp*(DyTvCkap1kap1i11TpYvCYvi21) - 4._dp*(DyTvCkap1kap2i11TpYvCYvi22) -& 
&  4._dp*(DyTvCkap1kap3i11TpYvCYvi23) - 4._dp*(DyTvCkap2kap1i12TpYvCYvi21) -             & 
&  4._dp*(DyTvCkap2kap2i12TpYvCYvi22) - 4._dp*(DyTvCkap2kap3i12TpYvCYvi23) -             & 
&  4._dp*(DyTvCkap3kap1i13TpYvCYvi21) - 4._dp*(DyTvCkap3kap2i13TpYvCYvi22) -             & 
&  4._dp*(DyTvCkap3kap3i13TpYvCYvi23) - 6._dp*(DyTvClami1TpYvCYvlami2) - 8._dp*(DyTvi11kap1adjYvYvCkap1i21) -& 
&  8._dp*(DyTvi11kap2adjYvYvCkap1i22) - 8._dp*(DyTvi11kap3adjYvYvCkap1i23) -             & 
&  8._dp*(DyTvi12kap1adjYvYvCkap2i21) - 8._dp*(DyTvi12kap2adjYvYvCkap2i22) -             & 
&  8._dp*(DyTvi12kap3adjYvYvCkap2i23) - 8._dp*(DyTvi13kap1adjYvYvCkap3i21) -             & 
&  8._dp*(DyTvi13kap2adjYvYvCkap3i22) - 8._dp*(DyTvi13kap3adjYvYvCkap3i23) -             & 
&  2._dp*(DyYvadjkap1kap1i11adjYvTv1i2) - 2._dp*(DyYvadjkap1kap2i11adjYvTv2i2) -         & 
&  2._dp*(DyYvadjkap1kap3i11adjYvTv3i2) - 2._dp*(DyYvadjkap2kap1i12adjYvTv1i2) -         & 
&  2._dp*(DyYvadjkap2kap2i12adjYvTv2i2) - 2._dp*(DyYvadjkap2kap3i12adjYvTv3i2) -         & 
&  2._dp*(DyYvadjkap3kap1i13adjYvTv1i2) - 2._dp*(DyYvadjkap3kap2i13adjYvTv2i2) -         & 
&  2._dp*(DyYvadjkap3kap3i13adjYvTv3i2) - 8._dp*(DyYvadjYvTvClami1lami2) -               & 
&  6._dp*(DyYvadjYvYvClami1Tlami2) - 16._dp*(DyYvCkap1i11kap1adjYvTvi21) -               & 
&  16._dp*(DyYvCkap1i11kap1Ckap1Tk1i21) - 16._dp*(DyYvCkap1i11kap1Ckap2Tk1i22) -         & 
&  16._dp*(DyYvCkap1i11kap1Ckap3Tk1i23) - 16._dp*(DyYvCkap1i11Tk1Ckap1kap1i21) -         & 
&  16._dp*(DyYvCkap1i11Tk1Ckap2kap1i22) - 16._dp*(DyYvCkap1i11Tk1Ckap3kap1i23) -         & 
&  16._dp*(DyYvCkap1i12kap2adjYvTvi21) - 16._dp*(DyYvCkap1i12kap2Ckap1Tk1i21) -          & 
&  16._dp*(DyYvCkap1i12kap2Ckap2Tk1i22) - 16._dp*(DyYvCkap1i12kap2Ckap3Tk1i23)  
betaTv2 =  betaTv2- 16._dp*(DyYvCkap1i12Tk1Ckap1kap2i21) - 16._dp*(DyYvCkap1i12Tk1Ckap2kap2i22) -         & 
&  16._dp*(DyYvCkap1i12Tk1Ckap3kap2i23) - 16._dp*(DyYvCkap1i13kap3adjYvTvi21) -          & 
&  16._dp*(DyYvCkap1i13kap3Ckap1Tk1i21) - 16._dp*(DyYvCkap1i13kap3Ckap2Tk1i22) -         & 
&  16._dp*(DyYvCkap1i13kap3Ckap3Tk1i23) - 16._dp*(DyYvCkap1i13Tk1Ckap1kap3i21) -         & 
&  16._dp*(DyYvCkap1i13Tk1Ckap2kap3i22) - 16._dp*(DyYvCkap1i13Tk1Ckap3kap3i23) -         & 
&  16._dp*(DyYvCkap1lami1Tk1Clami2) - 16._dp*(DyYvCkap2i11kap1adjYvTvi22) -              & 
&  16._dp*(DyYvCkap2i11kap1Ckap1Tk2i21) - 16._dp*(DyYvCkap2i11kap1Ckap2Tk2i22) -         & 
&  16._dp*(DyYvCkap2i11kap1Ckap3Tk2i23) - 16._dp*(DyYvCkap2i11Tk2Ckap1kap1i21) -         & 
&  16._dp*(DyYvCkap2i11Tk2Ckap2kap1i22) - 16._dp*(DyYvCkap2i11Tk2Ckap3kap1i23) -         & 
&  16._dp*(DyYvCkap2i12kap2adjYvTvi22) - 16._dp*(DyYvCkap2i12kap2Ckap1Tk2i21) -          & 
&  16._dp*(DyYvCkap2i12kap2Ckap2Tk2i22) - 16._dp*(DyYvCkap2i12kap2Ckap3Tk2i23) -         & 
&  16._dp*(DyYvCkap2i12Tk2Ckap1kap2i21) - 16._dp*(DyYvCkap2i12Tk2Ckap2kap2i22) -         & 
&  16._dp*(DyYvCkap2i12Tk2Ckap3kap2i23) - 16._dp*(DyYvCkap2i13kap3adjYvTvi22) -          & 
&  16._dp*(DyYvCkap2i13kap3Ckap1Tk2i21) - 16._dp*(DyYvCkap2i13kap3Ckap2Tk2i22) -         & 
&  16._dp*(DyYvCkap2i13kap3Ckap3Tk2i23) - 16._dp*(DyYvCkap2i13Tk2Ckap1kap3i21) -         & 
&  16._dp*(DyYvCkap2i13Tk2Ckap2kap3i22) - 16._dp*(DyYvCkap2i13Tk2Ckap3kap3i23) -         & 
&  16._dp*(DyYvCkap2lami1Tk2Clami2) - 16._dp*(DyYvCkap3i11kap1adjYvTvi23) -              & 
&  16._dp*(DyYvCkap3i11kap1Ckap1Tk3i21) - 16._dp*(DyYvCkap3i11kap1Ckap2Tk3i22) -         & 
&  16._dp*(DyYvCkap3i11kap1Ckap3Tk3i23) - 16._dp*(DyYvCkap3i11Tk3Ckap1kap1i21) -         & 
&  16._dp*(DyYvCkap3i11Tk3Ckap2kap1i22) - 16._dp*(DyYvCkap3i11Tk3Ckap3kap1i23) -         & 
&  16._dp*(DyYvCkap3i12kap2adjYvTvi23) - 16._dp*(DyYvCkap3i12kap2Ckap1Tk3i21)  
betaTv2 =  betaTv2- 16._dp*(DyYvCkap3i12kap2Ckap2Tk3i22) - 16._dp*(DyYvCkap3i12kap2Ckap3Tk3i23) -         & 
&  16._dp*(DyYvCkap3i12Tk3Ckap1kap2i21) - 16._dp*(DyYvCkap3i12Tk3Ckap2kap2i22) -         & 
&  16._dp*(DyYvCkap3i12Tk3Ckap3kap2i23) - 16._dp*(DyYvCkap3i13kap3adjYvTvi23) -          & 
&  16._dp*(DyYvCkap3i13kap3Ckap1Tk3i21) - 16._dp*(DyYvCkap3i13kap3Ckap2Tk3i22) -         & 
&  16._dp*(DyYvCkap3i13kap3Ckap3Tk3i23) - 16._dp*(DyYvCkap3i13Tk3Ckap1kap3i21) -         & 
&  16._dp*(DyYvCkap3i13Tk3Ckap2kap3i22) - 16._dp*(DyYvCkap3i13Tk3Ckap3kap3i23) -         & 
&  16._dp*(DyYvCkap3lami1Tk3Clami2) - 6._dp*(DyYvClami1TpTvCYvlami2) - 8._dp*(DyYvClami1TpYvCYvTlami2) +& 
&  (6*DyTvClami1lami2*g1p2)/5._dp + (12*DyYvClami1Tlami2*g1p2)/5._dp + 6*DyTvClami1lami2*g2p2 +& 
&  12*DyYvClami1Tlami2*g2p2 - (12*DyYvClami1lami2*g1p2*M1)/5._dp - 12*DyYvClami1lami2*g2p2*M2 -& 
&  14*DyYvClami1lami2*SPClamxxTlam - 10*DyTvClami1lami2*SPlamxxClam - 11*DyYvClami1Tlami2*SPlamxxClam -& 
&  4._dp*(TpTeCYeTpYeCYeYv) + (12*g1p2*TpTeCYeYv)/5._dp - 2*SPlamxxClam*TpTeCYeYv -      & 
&  4._dp*(TpYeCYeTpTeCYeYv) - 2._dp*(TpYeCYeTpYeCYeTv) + (6*g1p2*TpYeCYeTv)/5._dp -      & 
&  SPlamxxClam*TpYeCYeTv - (12*g1p2*M1*TpYeCYeYv)/5._dp - 2*SPClamxxTlam*TpYeCYeYv -     & 
&  12*DyYvClami1lami2*TradjYdTd - 6*TpYeCYeYv*TradjYdTd - 4*DyYvClami1lami2*TradjYeTe -  & 
&  2*TpYeCYeYv*TradjYeTe - 18*DyYvClami1lami2*TradjYuTu - 6*DyYvClami1lami2*TradjYvTv -  & 
&  8*DyTvi11kap1Ckap1i21*TrCkap1kap1 - 8*DyTvi12kap1Ckap2i21*TrCkap1kap1 -               & 
&  8*DyTvi13kap1Ckap3i21*TrCkap1kap1 - 8*DyTvi11kap1Ckap1i22*TrCkap1kap2 -               & 
&  8*DyTvi12kap1Ckap2i22*TrCkap1kap2 - 8*DyTvi13kap1Ckap3i22*TrCkap1kap2 -               & 
&  8*DyTvi11kap1Ckap1i23*TrCkap1kap3 - 8*DyTvi12kap1Ckap2i23*TrCkap1kap3 -               & 
&  8*DyTvi13kap1Ckap3i23*TrCkap1kap3 - 8*DyTvi11kap2Ckap1i21*TrCkap2kap1 -               & 
&  8*DyTvi12kap2Ckap2i21*TrCkap2kap1 - 8*DyTvi13kap2Ckap3i21*TrCkap2kap1  
betaTv2 =  betaTv2- 8*DyTvi11kap2Ckap1i22*TrCkap2kap2 - 8*DyTvi12kap2Ckap2i22*TrCkap2kap2 -               & 
&  8*DyTvi13kap2Ckap3i22*TrCkap2kap2 - 8*DyTvi11kap2Ckap1i23*TrCkap2kap3 -               & 
&  8*DyTvi12kap2Ckap2i23*TrCkap2kap3 - 8*DyTvi13kap2Ckap3i23*TrCkap2kap3 -               & 
&  8*DyTvi11kap3Ckap1i21*TrCkap3kap1 - 8*DyTvi12kap3Ckap2i21*TrCkap3kap1 -               & 
&  8*DyTvi13kap3Ckap3i21*TrCkap3kap1 - 8*DyTvi11kap3Ckap1i22*TrCkap3kap2 -               & 
&  8*DyTvi12kap3Ckap2i22*TrCkap3kap2 - 8*DyTvi13kap3Ckap3i22*TrCkap3kap2 -               & 
&  8*DyTvi11kap3Ckap1i23*TrCkap3kap3 - 8*DyTvi12kap3Ckap2i23*TrCkap3kap3 -               & 
&  8*DyTvi13kap3Ckap3i23*TrCkap3kap3 - 6*DyTvClami1lami2*TrYdadjYd - 12*DyYvClami1Tlami2*TrYdadjYd -& 
&  6*TpTeCYeYv*TrYdadjYd - 3*TpYeCYeTv*TrYdadjYd - 2*DyTvClami1lami2*TrYeadjYe -         & 
&  4*DyYvClami1Tlami2*TrYeadjYe - 2*TpTeCYeYv*TrYeadjYe - TpYeCYeTv*TrYeadjYe -          & 
&  12*DyTvClami1lami2*TrYuadjYu - 15*DyYvClami1Tlami2*TrYuadjYu - 4*DyTvClami1lami2*TrYvadjYv -& 
&  5*DyYvClami1Tlami2*TrYvadjYv - 2._dp*(TvadjYvTpYeCYeYv) + (6*g1p2*TvadjYvYv)/5._dp +  & 
&  6*g2p2*TvadjYvYv - 4*SPlamxxClam*TvadjYvYv - 12*TrYuadjYu*TvadjYvYv - 4*TrYvadjYv*TvadjYvYv -& 
&  6._dp*(TvadjYvYvadjYvYv) - (414*g1p4*M1*Yv)/25._dp - (18*g1p2*g2p2*M1*Yv)/5._dp -     & 
&  (18*g1p2*g2p2*M2*Yv)/5._dp - 30*g2p4*M2*Yv - 12*SPadjYvYvClamxxTlam*Yv -              & 
&  12*SPClamxxTpTvCYvlam*Yv - 12*SPClamxxTlam*SPlamxxClam*Yv - 6*SPlamxxClam*TradjYdTd*Yv -& 
&  2*SPlamxxClam*TradjYeTe*Yv - 2*TradjYeTeCYvTpYv*Yv + (8*g1p2*TradjYuTu*Yv)/5._dp +    & 
&  32*g3p2*TradjYuTu*Yv - 2*TradjYvTpYeCYeTv*Yv - 4*TradjYvTvadjkap1kap1*Yv -            & 
&  4*TradjYvTvadjkap2kap2*Yv - 4*TradjYvTvadjkap3kap3*Yv - 6*SPClamxxTlam*TrYdadjYd*Yv - & 
&  6*TrYdadjYuTuadjYd*Yv - 2*SPClamxxTlam*TrYeadjYe*Yv - 6*TrYuadjYdTdadjYu*Yv -         & 
&  (8*g1p2*M1*TrYuadjYu*Yv)/5._dp - 32*g3p2*M3*TrYuadjYu*Yv - 36*TrYuadjYuTuadjYu*Yv  
betaTv2 =  betaTv2- 12*TrYvadjYvTvadjYv*Yv - 4*TrYvCkap1TpTk1adjYv*Yv - 4*TrYvCkap2TpTk2adjYv*Yv -        & 
&  4*TrYvCkap3TpTk3adjYv*Yv - 4._dp*(YvadjYvTpTeCYeYv) - 4._dp*(YvadjYvTpYeCYeTv) +      & 
&  (12*g1p2*YvadjYvTv)/5._dp + 12*g2p2*YvadjYvTv - 5*SPlamxxClam*YvadjYvTv -             & 
&  15*TrYuadjYu*YvadjYvTv - 5*TrYvadjYv*YvadjYvTv - 8._dp*(YvadjYvTvadjYvYv) -           & 
&  (12*g1p2*M1*YvadjYvYv)/5._dp - 12*g2p2*M2*YvadjYvYv - 6*SPClamxxTlam*YvadjYvYv -      & 
&  18*TradjYuTu*YvadjYvYv - 6*TradjYvTv*YvadjYvYv - 6._dp*(YvadjYvYvadjYvTv) -           & 
&  4._dp*(YvCkap1TpTk1adjYvYv) - 16._dp*(YvCkap1TpYvCYvTpTk1) - 4._dp*(YvCkap2TpTk2adjYvYv) -& 
&  16._dp*(YvCkap2TpYvCYvTpTk2) - 4._dp*(YvCkap3TpTk3adjYvYv) - 16._dp*(YvCkap3TpYvCYvTpTk3) -& 
&  4*Yv*adjYvTvadjkap1kap2(2,1) - 4*Yv*adjYvTvadjkap1kap3(3,1) - 4*Yv*adjYvTvadjkap2kap1(1,& 
& 2) - 4*Yv*adjYvTvadjkap2kap3(3,2) - 4*Yv*adjYvTvadjkap3kap1(1,3) - 4*Yv*adjYvTvadjkap3kap2(2,& 
& 3) - 8*DyTvi11kap2Ckap1i21*Ckap1kap1(2,1) - 8*DyTvi12kap2Ckap2i21*Ckap1kap1(2,         & 
& 1) - 8*DyTvi13kap2Ckap3i21*Ckap1kap1(2,1) - 8*DyTvi11kap3Ckap1i21*Ckap1kap1(3,         & 
& 1) - 8*DyTvi12kap3Ckap2i21*Ckap1kap1(3,1) - 8*DyTvi13kap3Ckap3i21*Ckap1kap1(3,         & 
& 1) - 8*DyTvi11kap2Ckap1i22*Ckap1kap2(2,1) - 8*DyTvi12kap2Ckap2i22*Ckap1kap2(2,         & 
& 1) - 8*DyTvi13kap2Ckap3i22*Ckap1kap2(2,1) - 8*DyTvi11kap3Ckap1i22*Ckap1kap2(3,         & 
& 1) - 8*DyTvi12kap3Ckap2i22*Ckap1kap2(3,1) - 8*DyTvi13kap3Ckap3i22*Ckap1kap2(3,         & 
& 1) - 8*DyTvi11kap2Ckap1i23*Ckap1kap3(2,1) - 8*DyTvi12kap2Ckap2i23*Ckap1kap3(2,         & 
& 1) - 8*DyTvi13kap2Ckap3i23*Ckap1kap3(2,1) - 8*DyTvi11kap3Ckap1i23*Ckap1kap3(3,         & 
& 1) - 8*DyTvi12kap3Ckap2i23*Ckap1kap3(3,1) - 8*DyTvi13kap3Ckap3i23*Ckap1kap3(3,         & 
& 1) - 8*DyTvi11kap1Ckap1i21*Ckap2kap1(1,2) - 8*DyTvi12kap1Ckap2i21*Ckap2kap1(1,         & 
& 2) - 8*DyTvi13kap1Ckap3i21*Ckap2kap1(1,2) - 8*DyTvi11kap3Ckap1i21*Ckap2kap1(3,         & 
& 2) - 8*DyTvi12kap3Ckap2i21*Ckap2kap1(3,2) - 8*DyTvi13kap3Ckap3i21*Ckap2kap1(3,         & 
& 2) - 8*DyTvi11kap1Ckap1i22*Ckap2kap2(1,2) - 8*DyTvi12kap1Ckap2i22*Ckap2kap2(1,         & 
& 2) - 8*DyTvi13kap1Ckap3i22*Ckap2kap2(1,2) - 8*DyTvi11kap3Ckap1i22*Ckap2kap2(3,         & 
& 2) - 8*DyTvi12kap3Ckap2i22*Ckap2kap2(3,2) - 8*DyTvi13kap3Ckap3i22*Ckap2kap2(3,         & 
& 2) - 8*DyTvi11kap1Ckap1i23*Ckap2kap3(1,2) - 8*DyTvi12kap1Ckap2i23*Ckap2kap3(1,         & 
& 2) - 8*DyTvi13kap1Ckap3i23*Ckap2kap3(1,2) - 8*DyTvi11kap3Ckap1i23*Ckap2kap3(3,         & 
& 2) - 8*DyTvi12kap3Ckap2i23*Ckap2kap3(3,2) - 8*DyTvi13kap3Ckap3i23*Ckap2kap3(3,         & 
& 2) - 8*DyTvi11kap1Ckap1i21*Ckap3kap1(1,3) - 8*DyTvi12kap1Ckap2i21*Ckap3kap1(1,         & 
& 3) - 8*DyTvi13kap1Ckap3i21*Ckap3kap1(1,3) - 8*DyTvi11kap2Ckap1i21*Ckap3kap1(2,         & 
& 3) - 8*DyTvi12kap2Ckap2i21*Ckap3kap1(2,3) - 8*DyTvi13kap2Ckap3i21*Ckap3kap1(2,         & 
& 3) - 8*DyTvi11kap1Ckap1i22*Ckap3kap2(1,3) - 8*DyTvi12kap1Ckap2i22*Ckap3kap2(1,         & 
& 3) - 8*DyTvi13kap1Ckap3i22*Ckap3kap2(1,3) - 8*DyTvi11kap2Ckap1i22*Ckap3kap2(2,         & 
& 3) - 8*DyTvi12kap2Ckap2i22*Ckap3kap2(2,3) - 8*DyTvi13kap2Ckap3i22*Ckap3kap2(2,         & 
& 3) - 8*DyTvi11kap1Ckap1i23*Ckap3kap3(1,3) - 8*DyTvi12kap1Ckap2i23*Ckap3kap3(1,         & 
& 3) - 8*DyTvi13kap1Ckap3i23*Ckap3kap3(1,3) - 8*DyTvi11kap2Ckap1i23*Ckap3kap3(2,         & 
& 3) - 8*DyTvi12kap2Ckap2i23*Ckap3kap3(2,3) - 8*DyTvi13kap2Ckap3i23*Ckap3kap3(2,         & 
& 3) - 4*DyTpkap1Ckap1TpTv1i1lami2*Conjg(lam(1)) - 4*DyTpkap1Ckap2TpTv2i1lami2*Conjg(lam(1))  
betaTv2 =  betaTv2- 4*DyTpkap1Ckap3TpTv3i1lami2*Conjg(lam(1)) - 8*DyTvi11kap1Ckap1lami2*Conjg(lam(1)) -   & 
&  8*DyTvi12kap1Ckap2lami2*Conjg(lam(1)) - 8*DyTvi13kap1Ckap3lami2*Conjg(lam(1)) -       & 
&  2*DyYvadjkap1kap1i11Tlami2*Conjg(lam(1)) - 4*DyYvadjkap1Tk1i11lami2*Conjg(lam(1)) -   & 
&  2*DyYvadjkap2kap1i12Tlami2*Conjg(lam(1)) - 4*DyYvadjkap2Tk1i12lami2*Conjg(lam(1)) -   & 
&  2*DyYvadjkap3kap1i13Tlami2*Conjg(lam(1)) - 4*DyYvadjkap3Tk1i13lami2*Conjg(lam(1)) -   & 
&  4*DyTpkap2Ckap1TpTv1i1lami2*Conjg(lam(2)) - 4*DyTpkap2Ckap2TpTv2i1lami2*Conjg(lam(2)) -& 
&  4*DyTpkap2Ckap3TpTv3i1lami2*Conjg(lam(2)) - 8*DyTvi11kap2Ckap1lami2*Conjg(lam(2)) -   & 
&  8*DyTvi12kap2Ckap2lami2*Conjg(lam(2)) - 8*DyTvi13kap2Ckap3lami2*Conjg(lam(2)) -       & 
&  2*DyYvadjkap1kap2i11Tlami2*Conjg(lam(2)) - 4*DyYvadjkap1Tk2i11lami2*Conjg(lam(2)) -   & 
&  2*DyYvadjkap2kap2i12Tlami2*Conjg(lam(2)) - 4*DyYvadjkap2Tk2i12lami2*Conjg(lam(2)) -   & 
&  2*DyYvadjkap3kap2i13Tlami2*Conjg(lam(2)) - 4*DyYvadjkap3Tk2i13lami2*Conjg(lam(2)) -   & 
&  4*DyTpkap3Ckap1TpTv1i1lami2*Conjg(lam(3)) - 4*DyTpkap3Ckap2TpTv2i1lami2*Conjg(lam(3)) -& 
&  4*DyTpkap3Ckap3TpTv3i1lami2*Conjg(lam(3)) - 8*DyTvi11kap3Ckap1lami2*Conjg(lam(3)) -   & 
&  8*DyTvi12kap3Ckap2lami2*Conjg(lam(3)) - 8*DyTvi13kap3Ckap3lami2*Conjg(lam(3)) -       & 
&  2*DyYvadjkap1kap3i11Tlami2*Conjg(lam(3)) - 4*DyYvadjkap1Tk3i11lami2*Conjg(lam(3)) -   & 
&  2*DyYvadjkap2kap3i12Tlami2*Conjg(lam(3)) - 4*DyYvadjkap2Tk3i12lami2*Conjg(lam(3)) -   & 
&  2*DyYvadjkap3kap3i13Tlami2*Conjg(lam(3)) - 4*DyYvadjkap3Tk3i13lami2*Conjg(lam(3)) +   & 
&  (207*g1p4*Tv)/50._dp + (9*g1p2*g2p2*Tv)/5._dp + (15*g2p4*Tv)/2._dp - 6*SPlamxxadjYvYvClam*Tv -& 
&  3*SPlamxxClamp2*Tv - 3*SPlamxxClam*TrYdadjYd*Tv - 3*TrYdadjYuYuadjYd*Tv -             & 
&  SPlamxxClam*TrYeadjYe*Tv + (4*g1p2*TrYuadjYu*Tv)/5._dp + 16*g3p2*TrYuadjYu*Tv -       & 
&  9*TrYuadjYuYuadjYu*Tv - TrYvadjYvTpYeCYe*Tv - 3*TrYvadjYvYvadjYv*Tv - 2*TrYvCkap1kap1adjYv*Tv  
betaTv2 =  betaTv2- 2*TrYvCkap2kap2adjYv*Tv - 2*TrYvCkap3kap3adjYv*Tv - 2*adjYvYvCkap1kap2(2,             & 
& 1)*Tv - 2*adjYvYvCkap1kap3(3,1)*Tv - 2*adjYvYvCkap2kap1(1,2)*Tv - 2*adjYvYvCkap2kap3(3,& 
& 2)*Tv - 2*adjYvYvCkap3kap1(1,3)*Tv - 2*adjYvYvCkap3kap2(2,3)*Tv - 2*Conjg(lam(1))*Tv*Tpkap1Ckap1lam(1) -& 
&  4*Yv*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) - 2*Conjg(lam(1))*Tv*Tpkap1Ckap2lam(2) -        & 
&  4*Yv*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) - 2*Conjg(lam(1))*Tv*Tpkap1Ckap3lam(3) -        & 
&  4*Yv*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) - 2*Conjg(lam(2))*Tv*Tpkap2Ckap1lam(1) -        & 
&  4*Yv*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) - 2*Conjg(lam(2))*Tv*Tpkap2Ckap2lam(2) -        & 
&  4*Yv*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) - 2*Conjg(lam(2))*Tv*Tpkap2Ckap3lam(3) -        & 
&  4*Yv*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) - 2*Conjg(lam(3))*Tv*Tpkap3Ckap1lam(1) -        & 
&  4*Yv*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) - 2*Conjg(lam(3))*Tv*Tpkap3Ckap2lam(2) -        & 
&  4*Yv*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) - 2*Conjg(lam(3))*Tv*Tpkap3Ckap3lam(3) -        & 
&  4*Yv*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) - 4*Yv*Conjg(lam(1))*TpTk1Ckap1lam(1) -         & 
&  4*Yv*Conjg(lam(1))*TpTk1Ckap2lam(2) - 4*Yv*Conjg(lam(1))*TpTk1Ckap3lam(3) -           & 
&  4*Yv*Conjg(lam(2))*TpTk2Ckap1lam(1) - 4*Yv*Conjg(lam(2))*TpTk2Ckap2lam(2) -           & 
&  4*Yv*Conjg(lam(2))*TpTk2Ckap3lam(3) - 4*Yv*Conjg(lam(3))*TpTk3Ckap1lam(1) -           & 
&  4*Yv*Conjg(lam(3))*TpTk3Ckap2lam(2) - 4*Yv*Conjg(lam(3))*TpTk3Ckap3lam(3) -           & 
&  16*DyYvCkap1i11kap1Clami2*Tlam(1) - 16*DyYvCkap1i12kap2Clami2*Tlam(1) -               & 
&  16*DyYvCkap1i13kap3Clami2*Tlam(1) - 16*DyYvCkap2i11kap1Clami2*Tlam(2) -               & 
&  16*DyYvCkap2i12kap2Clami2*Tlam(2) - 16*DyYvCkap2i13kap3Clami2*Tlam(2) -               & 
&  16*DyYvCkap3i11kap1Clami2*Tlam(3) - 16*DyYvCkap3i12kap2Clami2*Tlam(3) -               & 
&  16*DyYvCkap3i13kap3Clami2*Tlam(3)

 
DTv = oo16pi2*( betaTv1 + oo16pi2 * betaTv2 ) 

 
Else 
DTv = oo16pi2* betaTv1 
End If 
 
 
Call Chop(DTv) 

!-------------------- 
! Tu 
!-------------------- 
 
betaTu1  = TuadjYdYd + 5._dp*(TuadjYuYu) + (26*g1p2*M1*Yu)/15._dp + (32*g3p2*M3*Yu)   & 
& /3._dp + 6*g2p2*M2*Yu + 2*SPClamxxTlam*Yu + 6*TradjYuTu*Yu + 2*TradjYvTv*Yu +          & 
&  2._dp*(YuadjYdTd) + 4._dp*(YuadjYuTu) - (13*g1p2*Tu)/15._dp - 3*g2p2*Tu -             & 
&  (16*g3p2*Tu)/3._dp + SPlamxxClam*Tu + 3*TrYuadjYu*Tu + TrYvadjYv*Tu

 
 
If (TwoLoopRGE) Then 
betaTu2 = ((2._dp*(g1p2) - 5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe))*TuadjYdYd)/5._dp -  & 
&  2._dp*(TuadjYdYdadjYdYd) - 4._dp*(TuadjYdYdadjYuYu) + 12*g2p2*TuadjYuYu -             & 
&  5*SPlamxxClam*TuadjYuYu - 15*TrYuadjYu*TuadjYuYu - 5*TrYvadjYv*TuadjYuYu -            & 
&  6._dp*(TuadjYuYuadjYuYu) - (5486*g1p4*M1*Yu)/225._dp - 2*g1p2*g2p2*M1*Yu -            & 
&  (272*g1p2*g3p2*M1*Yu)/45._dp - (272*g1p2*g3p2*M3*Yu)/45._dp - 16*g2p2*g3p2*M3*Yu +    & 
&  (64*g3p4*M3*Yu)/9._dp - 2*g1p2*g2p2*M2*Yu - 30*g2p4*M2*Yu - 16*g2p2*g3p2*M2*Yu -      & 
&  12*SPadjYvYvClamxxTlam*Yu - 12*SPClamxxTpTvCYvlam*Yu - 12*SPClamxxTlam*SPlamxxClam*Yu -& 
&  6*SPlamxxClam*TradjYdTd*Yu - 2*SPlamxxClam*TradjYeTe*Yu - 2*TradjYeTeCYvTpYv*Yu +     & 
&  (8*g1p2*TradjYuTu*Yu)/5._dp + 32*g3p2*TradjYuTu*Yu - 2*TradjYvTpYeCYeTv*Yu -          & 
&  4*TradjYvTvadjkap1kap1*Yu - 4*TradjYvTvadjkap2kap2*Yu - 4*TradjYvTvadjkap3kap3*Yu -   & 
&  6*SPClamxxTlam*TrYdadjYd*Yu - 6*TrYdadjYuTuadjYd*Yu - 2*SPClamxxTlam*TrYeadjYe*Yu -   & 
&  6*TrYuadjYdTdadjYu*Yu - (8*g1p2*M1*TrYuadjYu*Yu)/5._dp - 32*g3p2*M3*TrYuadjYu*Yu -    & 
&  36*TrYuadjYuTuadjYu*Yu - 12*TrYvadjYvTvadjYv*Yu - 4*TrYvCkap1TpTk1adjYv*Yu -          & 
&  4*TrYvCkap2TpTk2adjYv*Yu - 4*TrYvCkap3TpTk3adjYv*Yu + (4*g1p2*YuadjYdTd)/5._dp -      & 
&  2*SPlamxxClam*YuadjYdTd - 6*TrYdadjYd*YuadjYdTd - 2*TrYeadjYe*YuadjYdTd -             & 
&  4._dp*(YuadjYdTdadjYdYd) - 4._dp*(YuadjYdTdadjYuYu) - (4*g1p2*M1*YuadjYdYd)/5._dp -   & 
&  2*SPClamxxTlam*YuadjYdYd - 6*TradjYdTd*YuadjYdYd - 2*TradjYeTe*YuadjYdYd -            & 
&  4._dp*(YuadjYdYdadjYdTd) - 2._dp*(YuadjYdYdadjYuTu) + (6*g1p2*YuadjYuTu)/5._dp +      & 
&  6*g2p2*YuadjYuTu - 4*SPlamxxClam*YuadjYuTu - 12*TrYuadjYu*YuadjYuTu - 4*TrYvadjYv*YuadjYuTu -& 
&  8._dp*(YuadjYuTuadjYuYu) - (4*g1p2*M1*YuadjYuYu)/5._dp - 12*g2p2*M2*YuadjYuYu -       & 
&  6*SPClamxxTlam*YuadjYuYu - 18*TradjYuTu*YuadjYuYu - 6*TradjYvTv*YuadjYuYu  
betaTu2 =  betaTu2- 6._dp*(YuadjYuYuadjYuTu) - 4*Yu*adjYvTvadjkap1kap2(2,1) - 4*Yu*adjYvTvadjkap1kap3(3,  & 
& 1) - 4*Yu*adjYvTvadjkap2kap1(1,2) - 4*Yu*adjYvTvadjkap2kap3(3,2) - 4*Yu*adjYvTvadjkap3kap1(1,& 
& 3) - 4*Yu*adjYvTvadjkap3kap2(2,3) + (2743*g1p4*Tu)/450._dp + g1p2*g2p2*Tu +            & 
&  (15*g2p4*Tu)/2._dp + (136*g1p2*g3p2*Tu)/45._dp + 8*g2p2*g3p2*Tu - (16*g3p4*Tu)/9._dp -& 
&  6*SPlamxxadjYvYvClam*Tu - 3*SPlamxxClamp2*Tu - 3*SPlamxxClam*TrYdadjYd*Tu -           & 
&  3*TrYdadjYuYuadjYd*Tu - SPlamxxClam*TrYeadjYe*Tu + (4*g1p2*TrYuadjYu*Tu)/5._dp +      & 
&  16*g3p2*TrYuadjYu*Tu - 9*TrYuadjYuYuadjYu*Tu - TrYvadjYvTpYeCYe*Tu - 3*TrYvadjYvYvadjYv*Tu -& 
&  2*TrYvCkap1kap1adjYv*Tu - 2*TrYvCkap2kap2adjYv*Tu - 2*TrYvCkap3kap3adjYv*Tu -         & 
&  2*adjYvYvCkap1kap2(2,1)*Tu - 2*adjYvYvCkap1kap3(3,1)*Tu - 2*adjYvYvCkap2kap1(1,       & 
& 2)*Tu - 2*adjYvYvCkap2kap3(3,2)*Tu - 2*adjYvYvCkap3kap1(1,3)*Tu - 2*adjYvYvCkap3kap2(2,& 
& 3)*Tu - 2*Conjg(lam(1))*Tu*Tpkap1Ckap1lam(1) - 4*Yu*Conjg(lam(1))*Tpkap1Ckap1Tlam(1) - & 
&  2*Conjg(lam(1))*Tu*Tpkap1Ckap2lam(2) - 4*Yu*Conjg(lam(1))*Tpkap1Ckap2Tlam(2) -        & 
&  2*Conjg(lam(1))*Tu*Tpkap1Ckap3lam(3) - 4*Yu*Conjg(lam(1))*Tpkap1Ckap3Tlam(3) -        & 
&  2*Conjg(lam(2))*Tu*Tpkap2Ckap1lam(1) - 4*Yu*Conjg(lam(2))*Tpkap2Ckap1Tlam(1) -        & 
&  2*Conjg(lam(2))*Tu*Tpkap2Ckap2lam(2) - 4*Yu*Conjg(lam(2))*Tpkap2Ckap2Tlam(2) -        & 
&  2*Conjg(lam(2))*Tu*Tpkap2Ckap3lam(3) - 4*Yu*Conjg(lam(2))*Tpkap2Ckap3Tlam(3) -        & 
&  2*Conjg(lam(3))*Tu*Tpkap3Ckap1lam(1) - 4*Yu*Conjg(lam(3))*Tpkap3Ckap1Tlam(1) -        & 
&  2*Conjg(lam(3))*Tu*Tpkap3Ckap2lam(2) - 4*Yu*Conjg(lam(3))*Tpkap3Ckap2Tlam(2) -        & 
&  2*Conjg(lam(3))*Tu*Tpkap3Ckap3lam(3) - 4*Yu*Conjg(lam(3))*Tpkap3Ckap3Tlam(3) -        & 
&  4*Yu*Conjg(lam(1))*TpTk1Ckap1lam(1) - 4*Yu*Conjg(lam(1))*TpTk1Ckap2lam(2) -           & 
&  4*Yu*Conjg(lam(1))*TpTk1Ckap3lam(3) - 4*Yu*Conjg(lam(2))*TpTk2Ckap1lam(1)  
betaTu2 =  betaTu2- 4*Yu*Conjg(lam(2))*TpTk2Ckap2lam(2) - 4*Yu*Conjg(lam(2))*TpTk2Ckap3lam(3) -           & 
&  4*Yu*Conjg(lam(3))*TpTk3Ckap1lam(1) - 4*Yu*Conjg(lam(3))*TpTk3Ckap2lam(2) -           & 
&  4*Yu*Conjg(lam(3))*TpTk3Ckap3lam(3)

 
DTu = oo16pi2*( betaTu1 + oo16pi2 * betaTu2 ) 

 
Else 
DTu = oo16pi2* betaTu1 
End If 
 
 
Call Chop(DTu) 

!-------------------- 
! Tk 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaTk1(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk1(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk1(i1,i2,3) = 0
End Do
End Do

 
 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaTk2(i1,i2,1) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk2(i1,i2,2) = 0
End Do
End Do
Do i1 = 1,3
Do i2 = 1,3
betaTk2(i1,i2,3) = 0
End Do
End Do

 
DTk = oo16pi2*( betaTk1 + oo16pi2 * betaTk2 ) 

 
Else 
DTk = oo16pi2* betaTk1 
End If 
 
 
!-------------------- 
! mq2 
!-------------------- 
 
betamq21  = 2._dp*(adjTdTd) + 2._dp*(adjTuTu) + 2._dp*(adjYdmd2Yd) + adjYdYdmq2 +     & 
&  2._dp*(adjYumu2Yu) + adjYuYumq2 - (2*AbsM1*g1p2*id3R)/15._dp - 6*AbsM2*g2p2*id3R -    & 
&  (32*AbsM3*g3p2*id3R)/3._dp + 2*adjYdYd*mHd2 + 2*adjYuYu*mHu2 + mq2adjYdYd +           & 
&  mq2adjYuYu + g1*id3R*ooSqrt15*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamq22 = -4._dp*(adjTdTdadjYdYd) - 4._dp*(adjTdYdadjYdTd) - 4._dp*(adjTuTuadjYuYu) -           & 
&  4._dp*(adjTuYuadjYuTu) - 4._dp*(adjYdmd2YdadjYdYd) - 4._dp*(adjYdTdadjTdYd) -         & 
&  4._dp*(adjYdYdadjTdTd) - 4._dp*(adjYdYdadjYdmd2Yd) - 2._dp*(adjYdYdadjYdYdmq2) -      & 
&  4._dp*(adjYdYdmq2adjYdYd) - 4._dp*(adjYumu2YuadjYuYu) - 4._dp*(adjYuTuadjTuYu) -      & 
&  4._dp*(adjYuYuadjTuTu) - 4._dp*(adjYuYuadjYumu2Yu) - 2._dp*(adjYuYuadjYuYumq2) -      & 
&  4._dp*(adjYuYumq2adjYuYu) + (8*adjTuTu*g1p2)/5._dp + (4*adjYdmd2Yd*g1p2)/5._dp +      & 
&  (8*AbsM1*adjYdYd*g1p2)/5._dp + (2*adjYdYdmq2*g1p2)/5._dp + (8*adjYumu2Yu*g1p2)/5._dp +& 
&  (16*AbsM1*adjYuYu*g1p2)/5._dp + (4*adjYuYumq2*g1p2)/5._dp + (199*AbsM1*g1p4*id3R)/75._dp +& 
&  (2*AbsM1*g1p2*g2p2*id3R)/5._dp + (2*AbsM2*g1p2*g2p2*id3R)/5._dp + 33*AbsM2*g2p4*id3R +& 
&  (32*AbsM1*g1p2*g3p2*id3R)/45._dp + (32*AbsM3*g1p2*g3p2*id3R)/45._dp + 32*AbsM2*g2p2*g3p2*id3R +& 
&  32*AbsM3*g2p2*g3p2*id3R - (128*AbsM3*g3p4*id3R)/3._dp - (4*adjTdYd*g1p2*M1)/5._dp -   & 
&  (8*adjTuYu*g1p2*M1)/5._dp - 8*adjYdYdadjYdYd*mHd2 + (4*adjYdYd*g1p2*mHd2)/5._dp -     & 
&  8*adjYuYuadjYuYu*mHu2 + (8*adjYuYu*g1p2*mHu2)/5._dp + (2*g1p2*mq2adjYdYd)/5._dp -     & 
&  2._dp*(mq2adjYdYdadjYdYd) + (4*g1p2*mq2adjYuYu)/5._dp - 2._dp*(mq2adjYuYuadjYuYu) -   & 
&  2*adjTdYd*SPClamxxTlam - 2*adjTuYu*SPClamxxTlam + 2*adjYdYd*SPlamxxadjYvmlHd2 +       & 
&  2*adjYuYu*SPlamxxadjYvmlHd2 - 2*adjTuTu*SPlamxxClam - 2*adjYdmd2Yd*SPlamxxClam -      & 
&  adjYdYdmq2*SPlamxxClam - 2*adjYumu2Yu*SPlamxxClam - adjYuYumq2*SPlamxxClam -          & 
&  4*adjYdYd*mHd2*SPlamxxClam - 2*adjYuYu*mHd2*SPlamxxClam - 2*adjYdYd*mHu2*SPlamxxClam -& 
&  4*adjYuYu*mHu2*SPlamxxClam - mq2adjYdYd*SPlamxxClam - mq2adjYuYu*SPlamxxClam -        & 
&  2*adjYdYd*SPlamxxCmv2Clam - 2*adjYuYu*SPlamxxCmv2Clam - 2*adjYdTd*SPlamxxCTlam -      & 
&  2*adjYuTu*SPlamxxCTlam - 2*adjYdYd*SPTlamxxCTlam - 2*adjYuYu*SPTlamxxCTlam  
betamq22 =  betamq22- 6*adjTdYd*TradjYdTd - 2*adjTdYd*TradjYeTe - 6*adjTuYu*TradjYuTu - 2*adjTuYu*TradjYvTv -& 
&  6*adjYdYd*TrCTdTpTd - 6*adjYdTd*TrCTdTpYd - 2*adjYdYd*TrCTeTpTe - 2*adjYdTd*TrCTeTpYe -& 
&  6*adjYuYu*TrCTuTpTu - 6*adjYuTu*TrCTuTpYu - 2*adjYuYu*TrCTvTpTv - 2*adjYuTu*TrCTvTpYv -& 
&  6*adjYdYd*Trmd2YdadjYd - 2*adjYdYd*Trme2YeadjYe - 2*adjYdYd*Trml2adjYeYe -            & 
&  6*adjYdYd*Trmq2adjYdYd - 6*adjYuYu*Trmq2adjYuYu - 6*adjYuYu*Trmu2YuadjYu -            & 
&  6*adjYdmd2Yd*TrYdadjYd - 3*adjYdYdmq2*TrYdadjYd - 12*adjYdYd*mHd2*TrYdadjYd -         & 
&  3*mq2adjYdYd*TrYdadjYd - 2*adjYdmd2Yd*TrYeadjYe - adjYdYdmq2*TrYeadjYe -              & 
&  4*adjYdYd*mHd2*TrYeadjYe - mq2adjYdYd*TrYeadjYe + (2*adjTdTd*(2._dp*(g1p2) -          & 
&  5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)))/5._dp - 6*adjTuTu*TrYuadjYu -       & 
&  6*adjYumu2Yu*TrYuadjYu - 3*adjYuYumq2*TrYuadjYu - 12*adjYuYu*mHu2*TrYuadjYu -         & 
&  3*mq2adjYuYu*TrYuadjYu - 2*adjTuTu*TrYvadjYv - 2*adjYumu2Yu*TrYvadjYv -               & 
&  adjYuYumq2*TrYvadjYv - 4*adjYuYu*mHu2*TrYvadjYv - mq2adjYuYu*TrYvadjYv -              & 
&  2*adjYuYu*TrYvadjYvCml2 - 2*adjYuYu*TrYvCmv2adjYv - (4*adjYdTd*g1p2*Conjg(M1))/5._dp -& 
&  (8*adjYuTu*g1p2*Conjg(M1))/5._dp + (16*g1p2*g3p2*id3R*M3*Conjg(M1))/45._dp +          & 
&  (g1p2*g2p2*id3R*M2*Conjg(M1))/5._dp + (16*g1p2*g3p2*id3R*M1*Conjg(M3))/45._dp +       & 
&  16*g2p2*g3p2*id3R*M2*Conjg(M3) + (g1p2*g2p2*id3R*M1*Conjg(M2))/5._dp + 16*g2p2*g3p2*id3R*M3*Conjg(M2) +& 
&  6*g2p4*id3R*Tr2(2) + (32*g3p4*id3R*Tr2(3))/3._dp + (2*g1p2*id3R*Tr2U1(1,              & 
& 1))/15._dp + 4*g1*id3R*ooSqrt15*Tr3(1)

 
Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 ) 

 
Else 
Dmq2 = oo16pi2* betamq21 
End If 
 
 
Call Chop(Dmq2) 

Forall(i1=1:3) Dmq2(i1,i1) =  Real(Dmq2(i1,i1),dp) 
Dmq2 = 0.5_dp*(Dmq2+ Conjg(Transpose(Dmq2)) ) 
!-------------------- 
! ml2 
!-------------------- 
 
betaml21  = 2._dp*(adjTeTe) + 2._dp*(adjYeme2Ye) + adjYeYeml2 + 2._dp*(CTvTpTv)       & 
&  + 2._dp*(CYvmv2TpYv) + CYvTpYvml2 - DyCYvlami1mlHd2i2 - (6*AbsM1*g1p2*id3R)           & 
& /5._dp - 6*AbsM2*g2p2*id3R + 2*adjYeYe*mHd2 + 2*CYvTpYv*mHu2 + ml2adjYeYe +            & 
&  ml2CYvTpYv - g1*id3R*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betaml22 = -4._dp*(adjTeTeadjYeYe) - 4._dp*(adjTeYeadjYeTe) - 4._dp*(adjYeme2YeadjYeYe) -        & 
&  4._dp*(adjYeTeadjTeYe) - 4._dp*(adjYeYeadjTeTe) - 4._dp*(adjYeYeadjYeme2Ye) -         & 
&  2._dp*(adjYeYeadjYeYeml2) - 4._dp*(adjYeYeml2adjYeYe) - 4._dp*(CTvTpTvCYvTpYv) -      & 
&  4._dp*(CTvTpYvCYvTpTv) - 2._dp*(CYvkap1adjkap1TpYvml2) - 2._dp*(CYvkap2adjkap2TpYvml2) -& 
&  2._dp*(CYvkap3adjkap3TpYvml2) - 4._dp*(CYvmv2kap1adjkap1TpYv) - 4._dp*(CYvmv2kap2adjkap2TpYv) -& 
&  4._dp*(CYvmv2kap3adjkap3TpYv) - 4._dp*(CYvmv2TpYvCYvTpYv) - 4._dp*(CYvTk1adjTk1TpYv) -& 
&  4._dp*(CYvTk2adjTk2TpYv) - 4._dp*(CYvTk3adjTk3TpYv) - 4._dp*(CYvTpTvCTvTpYv) -        & 
&  4._dp*(CYvTpYvCTvTpTv) - 4._dp*(CYvTpYvCYvmv2TpYv) - 2._dp*(CYvTpYvCYvTpYvml2) -      & 
&  4._dp*(CYvTpYvml2CYvTpYv) - DyadjYeYeCYvlami1mlHd2i2 - 4._dp*(DyCTvi11YvCkap1Tk1i21) -& 
&  4._dp*(DyCTvi11YvCkap2Tk1i22) - 4._dp*(DyCTvi11YvCkap3Tk1i23) - 4._dp*(DyCTvi12YvCkap1Tk2i21) -& 
&  4._dp*(DyCTvi12YvCkap2Tk2i22) - 4._dp*(DyCTvi12YvCkap3Tk2i23) - 4._dp*(DyCTvi13YvCkap1Tk3i21) -& 
&  4._dp*(DyCTvi13YvCkap2Tk3i22) - 4._dp*(DyCTvi13YvCkap3Tk3i23) - 4._dp*(DyCTvlami1TvClami2) -& 
&  4._dp*(DyCTvTlami1YvClami2) - 8._dp*(DyCYvi11Yvadjkap1mv2kap1i21) - 8._dp*(DyCYvi11Yvadjkap2mv2kap1i22) -& 
&  8._dp*(DyCYvi11Yvadjkap3mv2kap1i23) - 4._dp*(DyCYvi11YvCmv2Ckap1kap1i21) -            & 
&  4._dp*(DyCYvi11YvCmv2Ckap2kap1i22) - 4._dp*(DyCYvi11YvCmv2Ckap3kap1i23) -             & 
&  8._dp*(DyCYvi12Yvadjkap1mv2kap2i21) - 8._dp*(DyCYvi12Yvadjkap2mv2kap2i22) -           & 
&  8._dp*(DyCYvi12Yvadjkap3mv2kap2i23) - 4._dp*(DyCYvi12YvCmv2Ckap1kap2i21) -            & 
&  4._dp*(DyCYvi12YvCmv2Ckap2kap2i22) - 4._dp*(DyCYvi12YvCmv2Ckap3kap2i23) -             & 
&  8._dp*(DyCYvi13Yvadjkap1mv2kap3i21) - 8._dp*(DyCYvi13Yvadjkap2mv2kap3i22) -           & 
&  8._dp*(DyCYvi13Yvadjkap3mv2kap3i23) - 4._dp*(DyCYvi13YvCmv2Ckap1kap3i21) -            & 
&  4._dp*(DyCYvi13YvCmv2Ckap2kap3i22) - 4._dp*(DyCYvi13YvCmv2Ckap3kap3i23)  
betaml22 =  betaml22- 2._dp*(DyCYvlami1Cml2YvClami2) - 4._dp*(DyCYvlami1TvCTlami2) + 4._dp*(DyCYvlami1YvadjYvmlHd2i2) -& 
&  4._dp*(DyCYvlami1YvCmv2Clami2) - 4._dp*(DyCYvmv2lami1YvClami2) - 4._dp*(DyCYvTlami1YvCTlami2) +& 
&  2._dp*(DyCYvTpYvCYvlami1mlHd2i2) - 2._dp*(Dyml2CYvi11Yvadjkap1kap1i21) -              & 
&  2._dp*(Dyml2CYvi11Yvadjkap2kap1i22) - 2._dp*(Dyml2CYvi11Yvadjkap3kap1i23) -           & 
&  2._dp*(Dyml2CYvi12Yvadjkap1kap2i21) - 2._dp*(Dyml2CYvi12Yvadjkap2kap2i22) -           & 
&  2._dp*(Dyml2CYvi12Yvadjkap3kap2i23) - 2._dp*(Dyml2CYvi13Yvadjkap1kap3i21) -           & 
&  2._dp*(Dyml2CYvi13Yvadjkap2kap3i22) - 2._dp*(Dyml2CYvi13Yvadjkap3kap3i23) -           & 
&  2._dp*(Dyml2CYvlami1YvClami2) + (12*adjYeme2Ye*g1p2)/5._dp + (24*AbsM1*adjYeYe*g1p2)/5._dp +& 
&  (6*adjYeYeml2*g1p2)/5._dp + (621*AbsM1*g1p4*id3R)/25._dp + (18*AbsM1*g1p2*g2p2*id3R)/5._dp +& 
&  (18*AbsM2*g1p2*g2p2*id3R)/5._dp + 33*AbsM2*g2p4*id3R - (12*adjTeYe*g1p2*M1)/5._dp -   & 
&  8*adjYeYeadjYeYe*mHd2 - 4*DyCYvlami1YvClami2*mHd2 + (12*adjYeYe*g1p2*mHd2)/5._dp -    & 
&  8*CYvTpYvCYvTpYv*mHu2 - 4*DyCYvi11Yvadjkap1kap1i21*mHu2 - 4*DyCYvi11Yvadjkap2kap1i22*mHu2 -& 
&  4*DyCYvi11Yvadjkap3kap1i23*mHu2 - 4*DyCYvi12Yvadjkap1kap2i21*mHu2 - 4*DyCYvi12Yvadjkap2kap2i22*mHu2 -& 
&  4*DyCYvi12Yvadjkap3kap2i23*mHu2 - 4*DyCYvi13Yvadjkap1kap3i21*mHu2 - 4*DyCYvi13Yvadjkap2kap3i22*mHu2 -& 
&  4*DyCYvi13Yvadjkap3kap3i23*mHu2 - 8*DyCYvlami1YvClami2*mHu2 + (6*g1p2*ml2adjYeYe)/5._dp -& 
&  2._dp*(ml2adjYeYeadjYeYe) - 2._dp*(ml2CYvTpYvCYvTpYv) - 2*adjTeYe*SPClamxxTlam -      & 
&  2*CTvTpYv*SPClamxxTlam + 2*adjYeYe*SPlamxxadjYvmlHd2 + 2*CYvTpYv*SPlamxxadjYvmlHd2 -  & 
&  2*adjYeme2Ye*SPlamxxClam - adjYeYeml2*SPlamxxClam - 2*CTvTpTv*SPlamxxClam -           & 
&  2*CYvmv2TpYv*SPlamxxClam - CYvTpYvml2*SPlamxxClam + 3*DyCYvlami1mlHd2i2*SPlamxxClam - & 
&  4*adjYeYe*mHd2*SPlamxxClam - 2*CYvTpYv*mHd2*SPlamxxClam - 2*adjYeYe*mHu2*SPlamxxClam -& 
&  4*CYvTpYv*mHu2*SPlamxxClam - ml2adjYeYe*SPlamxxClam - ml2CYvTpYv*SPlamxxClam  
betaml22 =  betaml22- 2*adjYeYe*SPlamxxCmv2Clam - 2*CYvTpYv*SPlamxxCmv2Clam - 2*adjYeTe*SPlamxxCTlam -      & 
&  2*CYvTpTv*SPlamxxCTlam - 2*adjYeYe*SPTlamxxCTlam - 2*CYvTpYv*SPTlamxxCTlam -          & 
&  6*adjTeYe*TradjYdTd - 2*adjTeYe*TradjYeTe - 6*CTvTpYv*TradjYuTu - 2*CTvTpYv*TradjYvTv -& 
&  4*DyCTvi11Tvi21*TrCkap1Tpkap1 - 4*DyCTvi12Tvi21*TrCkap1Tpkap2 - 4*DyCTvi13Tvi21*TrCkap1Tpkap3 -& 
&  4*DyCTvi11Tvi22*TrCkap2Tpkap1 - 4*DyCTvi12Tvi22*TrCkap2Tpkap2 - 4*DyCTvi13Tvi22*TrCkap2Tpkap3 -& 
&  4*DyCTvi11Tvi23*TrCkap3Tpkap1 - 4*DyCTvi12Tvi23*TrCkap3Tpkap2 - 4*DyCTvi13Tvi23*TrCkap3Tpkap3 -& 
&  6*adjYeYe*TrCTdTpTd - 6*adjYeTe*TrCTdTpYd - 2*adjYeYe*TrCTeTpTe - 2*adjYeTe*TrCTeTpYe -& 
&  4*DyCYvi11Tvi21*TrCTk1Tpkap1 - 4*DyCYvi12Tvi21*TrCTk1Tpkap2 - 4*DyCYvi13Tvi21*TrCTk1Tpkap3 -& 
&  4*DyCYvi11Tvi22*TrCTk2Tpkap1 - 4*DyCYvi12Tvi22*TrCTk2Tpkap2 - 4*DyCYvi13Tvi22*TrCTk2Tpkap3 -& 
&  4*DyCYvi11Tvi23*TrCTk3Tpkap1 - 4*DyCYvi12Tvi23*TrCTk3Tpkap2 - 4*DyCYvi13Tvi23*TrCTk3Tpkap3 -& 
&  6*CYvTpYv*TrCTuTpTu - 6*CYvTpTv*TrCTuTpYu - 2*CYvTpYv*TrCTvTpTv - 2*CYvTpTv*TrCTvTpYv -& 
&  6*adjYeYe*Trmd2YdadjYd - 2*adjYeYe*Trme2YeadjYe - 2*adjYeYe*Trml2adjYeYe -            & 
&  6*adjYeYe*Trmq2adjYdYd - 6*CYvTpYv*Trmq2adjYuYu - 6*CYvTpYv*Trmu2YuadjYu -            & 
&  6*adjYeme2Ye*TrYdadjYd - 3*adjYeYeml2*TrYdadjYd - 12*adjYeYe*mHd2*TrYdadjYd -         & 
&  3*ml2adjYeYe*TrYdadjYd - 2*adjYeme2Ye*TrYeadjYe - adjYeYeml2*TrYeadjYe -              & 
&  4*adjYeYe*mHd2*TrYeadjYe - ml2adjYeYe*TrYeadjYe + (2*adjTeTe*(6._dp*(g1p2) -          & 
&  5*(SPlamxxClam + 3._dp*(TrYdadjYd) + TrYeadjYe)))/5._dp - 6*CTvTpTv*TrYuadjYu -       & 
&  6*CYvmv2TpYv*TrYuadjYu - 3*CYvTpYvml2*TrYuadjYu + 3*DyCYvlami1mlHd2i2*TrYuadjYu -     & 
&  12*CYvTpYv*mHu2*TrYuadjYu - 3*ml2CYvTpYv*TrYuadjYu - 2*CTvTpTv*TrYvadjYv -            & 
&  2*CYvmv2TpYv*TrYvadjYv - CYvTpYvml2*TrYvadjYv + DyCYvlami1mlHd2i2*TrYvadjYv -         & 
&  4*CYvTpYv*mHu2*TrYvadjYv - ml2CYvTpYv*TrYvadjYv - 2*CYvTpYv*TrYvadjYvCml2  
betaml22 =  betaml22- 2*CYvTpYv*TrYvCmv2adjYv - (12*adjYeTe*g1p2*Conjg(M1))/5._dp + (9*g1p2*g2p2*id3R*M2*Conjg(M1))/5._dp +& 
&  (9*g1p2*g2p2*id3R*M1*Conjg(M2))/5._dp + 6*g2p4*id3R*Tr2(2) + (6*g1p2*id3R*Tr2U1(1,    & 
& 1))/5._dp - 4*g1*id3R*sqrt3ov5*Tr3(1) + 2*DyCYvTpkap1Ckap1i11mlHd2i2*lam(1) +          & 
&  2*DyCYvTpkap2Ckap1i12mlHd2i2*lam(1) + 2*DyCYvTpkap3Ckap1i13mlHd2i2*lam(1) +           & 
&  2*DyCYvTpkap1Ckap2i11mlHd2i2*lam(2) + 2*DyCYvTpkap2Ckap2i12mlHd2i2*lam(2) +           & 
&  2*DyCYvTpkap3Ckap2i13mlHd2i2*lam(2) + 2*DyCYvTpkap1Ckap3i11mlHd2i2*lam(3) +           & 
&  2*DyCYvTpkap2Ckap3i12mlHd2i2*lam(3) + 2*DyCYvTpkap3Ckap3i13mlHd2i2*lam(3)

 
Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 ) 

 
Else 
Dml2 = oo16pi2* betaml21 
End If 
 
 
Call Chop(Dml2) 

Forall(i1=1:3) Dml2(i1,i1) =  Real(Dml2(i1,i1),dp) 
Dml2 = 0.5_dp*(Dml2+ Conjg(Transpose(Dml2)) ) 
!-------------------- 
! mHd2 
!-------------------- 
 
betamHd21  = (-6*AbsM1*g1p2)/5._dp - 6*AbsM2*g2p2 - SPlamxxadjYvmlHd2 +               & 
&  2*mHd2*SPlamxxClam + 2*mHu2*SPlamxxClam + 2._dp*(SPlamxxCmv2Clam) + 2._dp*(SPTlamxxCTlam)& 
&  + 6._dp*(TrCTdTpTd) + 2._dp*(TrCTeTpTe) + 6._dp*(Trmd2YdadjYd) + 2._dp*(Trme2YeadjYe) & 
&  + 2._dp*(Trml2adjYeYe) + 6._dp*(Trmq2adjYdYd) + 6*mHd2*TrYdadjYd + 2*mHd2*TrYeadjYe - & 
&  g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHd22 = (621*AbsM1*g1p4)/25._dp + (18*AbsM1*g1p2*g2p2)/5._dp + (18*AbsM2*g1p2*g2p2)/5._dp +   & 
&  33*AbsM2*g2p4 - 4._dp*(SPadjTvYvClamxxTlam) - 4._dp*(SPadjYvYvCTlamxxTlam) -          & 
&  4._dp*(SPClamxxTpYvml2CYvlam) - 4._dp*(SPlamxxadjTvTvClam) + 2._dp*(SPlamxxadjYvYvadjYvmlHd2) -& 
&  4*mHd2*SPlamxxadjYvYvClam - 8*mHu2*SPlamxxadjYvYvClam - 4._dp*(SPlamxxadjYvYvCmv2Clam) +& 
&  9*SPlamxxadjYvmlHd2*SPlamxxClam - 12*mHd2*SPlamxxClamp2 - 12*mHu2*SPlamxxClamp2 -     & 
&  4._dp*(SPlamxxCmv2adjYvYvClam) - 12*SPlamxxClam*SPlamxxCmv2Clam - 12*SPClamxxTlam*SPlamxxCTlam +& 
&  SPmlHd2xxadjYeYeCYvlam - 12*SPlamxxClam*SPTlamxxCTlam - 4._dp*(SPTpTk1ClamxxadjTk1lam) -& 
&  4._dp*(SPTpTk2ClamxxadjTk2lam) - 4._dp*(SPTpTk3ClamxxadjTk3lam) - 4._dp*(SPTpTvCYvlamxxCTlam) -& 
&  2._dp*(TradjYeTeCTvTpYv) - 6*SPlamxxCTlam*TradjYuTu - 2._dp*(TradjYvTpYeCTeTv) -      & 
&  2*SPlamxxCTlam*TradjYvTv - (4*g1p2*TrCTdTpTd)/5._dp + 32*g3p2*TrCTdTpTd +             & 
&  (4*g1p2*M1*TrCTdTpYd)/5._dp - 32*g3p2*M3*TrCTdTpYd + (12*g1p2*TrCTeTpTe)/5._dp -      & 
&  (12*g1p2*M1*TrCTeTpYe)/5._dp - 6*SPlamxxClam*TrCTuTpTu - 6*SPClamxxTlam*TrCTuTpYu -   & 
&  2*SPlamxxClam*TrCTvTpTv - 2*SPClamxxTlam*TrCTvTpYv - (4*g1p2*Trmd2YdadjYd)/5._dp +    & 
&  32*g3p2*Trmd2YdadjYd - 36._dp*(Trmd2YdadjYdYdadjYd) - 6._dp*(Trmd2YdadjYuYuadjYd) +   & 
&  (12*g1p2*Trme2YeadjYe)/5._dp - 12._dp*(Trme2YeadjYeYeadjYe) + (12*g1p2*Trml2adjYeYe)/5._dp -& 
&  12._dp*(Trml2adjYeYeadjYeYe) - (4*g1p2*Trmq2adjYdYd)/5._dp + 32*g3p2*Trmq2adjYdYd -   & 
&  36._dp*(Trmq2adjYdYdadjYdYd) - 6._dp*(Trmq2adjYdYdadjYuYu) - 6*SPlamxxClam*Trmq2adjYuYu -& 
&  6._dp*(Trmq2adjYuYuadjYdYd) - 6._dp*(Trmu2YuadjYdYdadjYu) - 6*SPlamxxClam*Trmu2YuadjYu -& 
&  36._dp*(TrYdadjTdTdadjYd) - 6._dp*(TrYdadjTuTuadjYd) - (8*AbsM1*g1p2*TrYdadjYd)/5._dp +& 
&  64*AbsM3*g3p2*TrYdadjYd - (4*g1p2*mHd2*TrYdadjYd)/5._dp + 32*g3p2*mHd2*TrYdadjYd -    & 
&  36._dp*(TrYdadjYdTdadjTd) - 36*mHd2*TrYdadjYdYdadjYd - 6._dp*(TrYdadjYuTuadjTd)  
betamHd22 =  betamHd22- 6*mHd2*TrYdadjYuYuadjYd - 6*mHu2*TrYdadjYuYuadjYd - 12._dp*(TrYeadjTeTeadjYe) +       & 
&  (24*AbsM1*g1p2*TrYeadjYe)/5._dp + (12*g1p2*mHd2*TrYeadjYe)/5._dp - 12._dp*(TrYeadjYeTeadjTe) -& 
&  12*mHd2*TrYeadjYeYeadjYe - 2._dp*(TrYeCTvTpTvadjYe) - 6._dp*(TrYuadjTdTdadjYu) -      & 
&  6._dp*(TrYuadjYdTdadjTu) + 3*SPlamxxadjYvmlHd2*TrYuadjYu - 6*mHd2*SPlamxxClam*TrYuadjYu -& 
&  12*mHu2*SPlamxxClam*TrYuadjYu - 6*SPlamxxCmv2Clam*TrYuadjYu - 6*SPTlamxxCTlam*TrYuadjYu +& 
&  SPlamxxadjYvmlHd2*TrYvadjYv - 2*mHd2*SPlamxxClam*TrYvadjYv - 4*mHu2*SPlamxxClam*TrYvadjYv -& 
&  2*SPlamxxCmv2Clam*TrYvadjYv - 2*SPTlamxxCTlam*TrYvadjYv - 2*SPlamxxClam*TrYvadjYvCml2 -& 
&  2._dp*(TrYvadjYvCml2TpYeCYe) - 2._dp*(TrYvadjYvTpTeCTe) - 2._dp*(TrYvadjYvTpYeCme2CYe) -& 
&  2*mHd2*TrYvadjYvTpYeCYe - 2*mHu2*TrYvadjYvTpYeCYe - 2._dp*(TrYvadjYvTpYeCYeCml2) -    & 
&  2*SPlamxxClam*TrYvCmv2adjYv - 2._dp*(TrYvCmv2adjYvTpYeCYe) + (g1p2*(9*g2p2*M2 +       & 
&  4*(TradjYdTd - 3._dp*(TradjYeTe)))*Conjg(M1))/5._dp - 32*g3p2*TradjYdTd*Conjg(M3) +   & 
&  (9*g1p2*g2p2*M1*Conjg(M2))/5._dp + 2*adjYvmlHd2(1)*Tpkap1adjkap1lam(1) -              & 
&  4*Conjg(lam(1))*Tpkap1adjkap1mv2lam(1) + 2*adjYvmlHd2(1)*Tpkap1adjkap2lam(2) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap2mv2lam(2) + 2*adjYvmlHd2(1)*Tpkap1adjkap3lam(3) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap3mv2lam(3) - 2*mHd2*Conjg(lam(1))*Tpkap1Ckap1lam(1) -     & 
&  4*mHu2*Conjg(lam(1))*Tpkap1Ckap1lam(1) - 2*mHd2*Conjg(lam(1))*Tpkap1Ckap2lam(2) -     & 
&  4*mHu2*Conjg(lam(1))*Tpkap1Ckap2lam(2) - 2*mHd2*Conjg(lam(1))*Tpkap1Ckap3lam(3) -     & 
&  4*mHu2*Conjg(lam(1))*Tpkap1Ckap3lam(3) + 2*adjYvmlHd2(2)*Tpkap2adjkap1lam(1) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap1mv2lam(1) + 2*adjYvmlHd2(2)*Tpkap2adjkap2lam(2) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap2mv2lam(2) + 2*adjYvmlHd2(2)*Tpkap2adjkap3lam(3) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap3mv2lam(3) - 2*mHd2*Conjg(lam(2))*Tpkap2Ckap1lam(1)  
betamHd22 =  betamHd22- 4*mHu2*Conjg(lam(2))*Tpkap2Ckap1lam(1) - 2*mHd2*Conjg(lam(2))*Tpkap2Ckap2lam(2) -     & 
&  4*mHu2*Conjg(lam(2))*Tpkap2Ckap2lam(2) - 2*mHd2*Conjg(lam(2))*Tpkap2Ckap3lam(3) -     & 
&  4*mHu2*Conjg(lam(2))*Tpkap2Ckap3lam(3) + 2*adjYvmlHd2(3)*Tpkap3adjkap1lam(1) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap1mv2lam(1) + 2*adjYvmlHd2(3)*Tpkap3adjkap2lam(2) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap2mv2lam(2) + 2*adjYvmlHd2(3)*Tpkap3adjkap3lam(3) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap3mv2lam(3) - 2*mHd2*Conjg(lam(3))*Tpkap3Ckap1lam(1) -     & 
&  4*mHu2*Conjg(lam(3))*Tpkap3Ckap1lam(1) - 2*mHd2*Conjg(lam(3))*Tpkap3Ckap2lam(2) -     & 
&  4*mHu2*Conjg(lam(3))*Tpkap3Ckap2lam(2) - 2*mHd2*Conjg(lam(3))*Tpkap3Ckap3lam(3) -     & 
&  4*mHu2*Conjg(lam(3))*Tpkap3Ckap3lam(3) - 4*Conjg(Tlam(1))*TpTk1adjkap1lam(1) -        & 
&  4*Conjg(Tlam(1))*TpTk1adjkap2lam(2) - 4*Conjg(Tlam(1))*TpTk1adjkap3lam(3) -           & 
&  4*Conjg(Tlam(2))*TpTk2adjkap1lam(1) - 4*Conjg(Tlam(2))*TpTk2adjkap2lam(2) -           & 
&  4*Conjg(Tlam(2))*TpTk2adjkap3lam(3) - 4*Conjg(Tlam(3))*TpTk3adjkap1lam(1) -           & 
&  4*Conjg(Tlam(3))*TpTk3adjkap2lam(2) - 4*Conjg(Tlam(3))*TpTk3adjkap3lam(3) +           & 
&  6*g2p4*Tr2(2) + (6*g1p2*Tr2U1(1,1))/5._dp - 4*g1*sqrt3ov5*Tr3(1) - 2*mHd2*adjkap1kap1Clam(1)*lam(1) -& 
&  2*mHd2*adjkap1kap2Clam(2)*lam(1) - 2*mHd2*adjkap1kap3Clam(3)*lam(1) - 4*adjkap1Tpkap1Cmv2Clam(1)*lam(1) -& 
&  4*adjkap1Tpkap2Cmv2Clam(2)*lam(1) - 4*adjkap1Tpkap3Cmv2Clam(3)*lam(1) -               & 
&  8*Trmv2kap1adjkap1*Conjg(lam(1))*lam(1) - 8*Trmv2kap2adjkap1*Conjg(lam(2))*lam(1) -   & 
&  8*Trmv2kap3adjkap1*Conjg(lam(3))*lam(1) - 2*mHd2*adjkap2kap1Clam(1)*lam(2) -          & 
&  2*mHd2*adjkap2kap2Clam(2)*lam(2) - 2*mHd2*adjkap2kap3Clam(3)*lam(2) - 4*adjkap2Tpkap1Cmv2Clam(1)*lam(2) -& 
&  4*adjkap2Tpkap2Cmv2Clam(2)*lam(2) - 4*adjkap2Tpkap3Cmv2Clam(3)*lam(2) -               & 
&  8*Trmv2kap1adjkap2*Conjg(lam(1))*lam(2) - 8*Trmv2kap2adjkap2*Conjg(lam(2))*lam(2)  
betamHd22 =  betamHd22- 8*Trmv2kap3adjkap2*Conjg(lam(3))*lam(2) - 2*mHd2*adjkap3kap1Clam(1)*lam(3) -          & 
&  2*mHd2*adjkap3kap2Clam(2)*lam(3) - 2*mHd2*adjkap3kap3Clam(3)*lam(3) - 4*adjkap3Tpkap1Cmv2Clam(1)*lam(3) -& 
&  4*adjkap3Tpkap2Cmv2Clam(2)*lam(3) - 4*adjkap3Tpkap3Cmv2Clam(3)*lam(3) -               & 
&  8*Trmv2kap1adjkap3*Conjg(lam(1))*lam(3) - 8*Trmv2kap2adjkap3*Conjg(lam(2))*lam(3) -   & 
&  8*Trmv2kap3adjkap3*Conjg(lam(3))*lam(3) - 4*TrCTk1Tpkap1*Conjg(lam(1))*Tlam(1) -      & 
&  4*TrCTk1Tpkap2*Conjg(lam(2))*Tlam(1) - 4*TrCTk1Tpkap3*Conjg(lam(3))*Tlam(1) -         & 
&  4*TrCkap1Tpkap1*Conjg(Tlam(1))*Tlam(1) - 4*TrCkap1Tpkap2*Conjg(Tlam(2))*Tlam(1) -     & 
&  4*TrCkap1Tpkap3*Conjg(Tlam(3))*Tlam(1) - 4*TrCTk2Tpkap1*Conjg(lam(1))*Tlam(2) -       & 
&  4*TrCTk2Tpkap2*Conjg(lam(2))*Tlam(2) - 4*TrCTk2Tpkap3*Conjg(lam(3))*Tlam(2) -         & 
&  4*TrCkap2Tpkap1*Conjg(Tlam(1))*Tlam(2) - 4*TrCkap2Tpkap2*Conjg(Tlam(2))*Tlam(2) -     & 
&  4*TrCkap2Tpkap3*Conjg(Tlam(3))*Tlam(2) - 4*TrCTk3Tpkap1*Conjg(lam(1))*Tlam(3) -       & 
&  4*TrCTk3Tpkap2*Conjg(lam(2))*Tlam(3) - 4*TrCTk3Tpkap3*Conjg(lam(3))*Tlam(3) -         & 
&  4*TrCkap3Tpkap1*Conjg(Tlam(1))*Tlam(3) - 4*TrCkap3Tpkap2*Conjg(Tlam(2))*Tlam(3) -     & 
&  4*TrCkap3Tpkap3*Conjg(Tlam(3))*Tlam(3)

 
DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 ) 

 
Else 
DmHd2 = oo16pi2* betamHd21 
End If 
 
 
!-------------------- 
! mHu2 
!-------------------- 
 
betamHu21  = (-6*AbsM1*g1p2)/5._dp - 6*AbsM2*g2p2 - 2._dp*(SPlamxxadjYvmlHd2)         & 
&  + 2*mHd2*SPlamxxClam + 2*mHu2*SPlamxxClam + 2._dp*(SPlamxxCmv2Clam) + 2._dp*(SPTlamxxCTlam)& 
&  + 6._dp*(TrCTuTpTu) + 2._dp*(TrCTvTpTv) + 6._dp*(Trmq2adjYuYu) + 6._dp*(Trmu2YuadjYu) & 
&  + 6*mHu2*TrYuadjYu + 2*mHu2*TrYvadjYv + 2._dp*(TrYvadjYvCml2) + 2._dp*(TrYvCmv2adjYv) & 
&  + g1*sqrt3ov5*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHu22 = (621*AbsM1*g1p4)/25._dp + (18*AbsM1*g1p2*g2p2)/5._dp + (18*AbsM2*g1p2*g2p2)/5._dp +   & 
&  33*AbsM2*g2p4 - 12._dp*(SPadjTvYvClamxxTlam) - 12._dp*(SPadjYvYvCTlamxxTlam) -        & 
&  12._dp*(SPClamxxTpYvml2CYvlam) - 12._dp*(SPlamxxadjTvTvClam) + 12._dp*(SPlamxxadjYvYvadjYvmlHd2) -& 
&  12*mHd2*SPlamxxadjYvYvClam - 24*mHu2*SPlamxxadjYvYvClam - 12._dp*(SPlamxxadjYvYvCmv2Clam) +& 
&  12*SPlamxxadjYvmlHd2*SPlamxxClam - 12*mHd2*SPlamxxClamp2 - 12*mHu2*SPlamxxClamp2 -    & 
&  12._dp*(SPlamxxCmv2adjYvYvClam) - 12*SPlamxxClam*SPlamxxCmv2Clam - 12*SPClamxxTlam*SPlamxxCTlam -& 
&  12*SPlamxxClam*SPTlamxxCTlam - 4._dp*(SPTpTk1ClamxxadjTk1lam) - 4._dp*(SPTpTk2ClamxxadjTk2lam) -& 
&  4._dp*(SPTpTk3ClamxxadjTk3lam) - 12._dp*(SPTpTvCYvlamxxCTlam) - 6*SPlamxxCTlam*TradjYdTd -& 
&  2*SPlamxxCTlam*TradjYeTe - 2._dp*(TradjYeTeCTvTpYv) - 2._dp*(TradjYvTpYeCTeTv) -      & 
&  6*SPlamxxClam*TrCTdTpTd - 6*SPClamxxTlam*TrCTdTpYd - 2*SPlamxxClam*TrCTeTpTe -        & 
&  2*SPClamxxTlam*TrCTeTpYe - 4*TradjYvTv*TrCTk1Tpkap1 - 4*TradjYvTv*TrCTk2Tpkap2 -      & 
&  4*TradjYvTv*TrCTk3Tpkap3 + (8*g1p2*TrCTuTpTu)/5._dp + 32*g3p2*TrCTuTpTu -             & 
&  (8*g1p2*M1*TrCTuTpYu)/5._dp - 32*g3p2*M3*TrCTuTpYu - 4*TrCkap1Tpkap1*TrCTvTpTv -      & 
&  4*TrCkap2Tpkap2*TrCTvTpTv - 4*TrCkap3Tpkap3*TrCTvTpTv - 6*SPlamxxClam*Trmd2YdadjYd -  & 
&  6._dp*(Trmd2YdadjYuYuadjYd) - 2*SPlamxxClam*Trme2YeadjYe - 2*SPlamxxClam*Trml2adjYeYe -& 
&  6*SPlamxxClam*Trmq2adjYdYd - 6._dp*(Trmq2adjYdYdadjYuYu) + (8*g1p2*Trmq2adjYuYu)/5._dp +& 
&  32*g3p2*Trmq2adjYuYu - 6._dp*(Trmq2adjYuYuadjYdYd) - 36._dp*(Trmq2adjYuYuadjYuYu) -   & 
&  6._dp*(Trmu2YuadjYdYdadjYu) + (8*g1p2*Trmu2YuadjYu)/5._dp + 32*g3p2*Trmu2YuadjYu -    & 
&  36._dp*(Trmu2YuadjYuYuadjYu) - 8._dp*(Trmv2kap1adjYvYvadjkap1) - 8._dp*(Trmv2kap2adjYvYvadjkap2) -& 
&  8._dp*(Trmv2kap3adjYvYvadjkap3) - 6._dp*(TrYdadjTuTuadjYd) + 6*SPlamxxadjYvmlHd2*TrYdadjYd -& 
&  12*mHd2*SPlamxxClam*TrYdadjYd - 6*mHu2*SPlamxxClam*TrYdadjYd - 6*SPlamxxCmv2Clam*TrYdadjYd  
betamHu22 =  betamHu22- 6*SPTlamxxCTlam*TrYdadjYd - 6._dp*(TrYdadjYuTuadjTd) - 6*mHd2*TrYdadjYuYuadjYd -      & 
&  6*mHu2*TrYdadjYuYuadjYd + 2*SPlamxxadjYvmlHd2*TrYeadjYe - 4*mHd2*SPlamxxClam*TrYeadjYe -& 
&  2*mHu2*SPlamxxClam*TrYeadjYe - 2*SPlamxxCmv2Clam*TrYeadjYe - 2*SPTlamxxCTlam*TrYeadjYe -& 
&  2._dp*(TrYeCTvTpTvadjYe) - 6._dp*(TrYuadjTdTdadjYu) - 36._dp*(TrYuadjTuTuadjYu) -     & 
&  6._dp*(TrYuadjYdTdadjTu) + (16*AbsM1*g1p2*TrYuadjYu)/5._dp + 64*AbsM3*g3p2*TrYuadjYu +& 
&  (8*g1p2*mHu2*TrYuadjYu)/5._dp + 32*g3p2*mHu2*TrYuadjYu - 36._dp*(TrYuadjYuTuadjTu) -  & 
&  36*mHu2*TrYuadjYuYuadjYu - 12._dp*(TrYvadjTvTvadjYv) - 2._dp*(TrYvadjYvCml2TpYeCYe) - & 
&  12._dp*(TrYvadjYvCml2YvadjYv) - 2._dp*(TrYvadjYvTpTeCTe) - 2._dp*(TrYvadjYvTpYeCme2CYe) -& 
&  2*mHd2*TrYvadjYvTpYeCYe - 2*mHu2*TrYvadjYvTpYeCYe - 2._dp*(TrYvadjYvTpYeCYeCml2) -    & 
&  12._dp*(TrYvadjYvTvadjTv) - 12*mHu2*TrYvadjYvYvadjYv - 2*mHu2*TrYvCkap1kap1adjYv -    & 
&  4._dp*(TrYvCkap1kap1adjYvCml2) - 4._dp*(TrYvCkap1Tk1adjTv) - 2*mHu2*TrYvCkap1Tpkap1adjYv -& 
&  4._dp*(TrYvCkap1Tpkap1Cmv2adjYv) - 2*mHu2*TrYvCkap2kap2adjYv - 4._dp*(TrYvCkap2kap2adjYvCml2) -& 
&  4._dp*(TrYvCkap2Tk2adjTv) - 2*mHu2*TrYvCkap2Tpkap2adjYv - 4._dp*(TrYvCkap2Tpkap2Cmv2adjYv) -& 
&  2*mHu2*TrYvCkap3kap3adjYv - 4._dp*(TrYvCkap3kap3adjYvCml2) - 4._dp*(TrYvCkap3Tk3adjTv) -& 
&  2*mHu2*TrYvCkap3Tpkap3adjYv - 4._dp*(TrYvCkap3Tpkap3Cmv2adjYv) - 2._dp*(TrYvCmv2adjYvTpYeCYe) -& 
&  12._dp*(TrYvCmv2adjYvYvadjYv) - 4._dp*(TrYvCmv2Ckap1kap1adjYv) - 4._dp*(TrYvCmv2Ckap2kap2adjYv) -& 
&  4._dp*(TrYvCmv2Ckap3kap3adjYv) - 4._dp*(TrYvCTk1TpTk1adjYv) - 4._dp*(TrYvCTk2TpTk2adjYv) -& 
&  4._dp*(TrYvCTk3TpTk3adjYv) - 4*adjTvYvCkap1Tk2(2,1) - 4*adjTvYvCkap1Tk3(3,            & 
& 1) - 4*adjTvYvCkap2Tk1(1,2) - 4*adjTvYvCkap2Tk3(3,2) - 4*adjTvYvCkap3Tk1(1,            & 
& 3) - 4*adjTvYvCkap3Tk2(2,3) - 4*adjYvCml2YvCkap1kap2(2,1) - 4*adjYvCml2YvCkap1kap3(3,  & 
& 1) - 4*adjYvCml2YvCkap2kap1(1,2) - 4*adjYvCml2YvCkap2kap3(3,2) - 4*adjYvCml2YvCkap3kap1(1,& 
& 3) - 4*adjYvCml2YvCkap3kap2(2,3) - 4*TrCTk2Tpkap1*adjYvTv(1,2) - 4*TrCTk3Tpkap1*adjYvTv(1,& 
& 3) - 4*TrCTk1Tpkap2*adjYvTv(2,1) - 4*TrCTk3Tpkap2*adjYvTv(2,3) - 4*TrCTk1Tpkap3*adjYvTv(3,& 
& 1) - 4*TrCTk2Tpkap3*adjYvTv(3,2) - 8*adjYvYvadjkap1mv2kap2(2,1) - 8*adjYvYvadjkap1mv2kap3(3,& 
& 1) - 8*adjYvYvadjkap2mv2kap1(1,2) - 8*adjYvYvadjkap2mv2kap3(3,2) - 8*adjYvYvadjkap3mv2kap1(1,& 
& 3) - 8*adjYvYvadjkap3mv2kap2(2,3) - 2*mHu2*adjYvYvCkap1kap2(2,1) - 2*mHu2*adjYvYvCkap1kap3(3,& 
& 1) - 2*mHu2*adjYvYvCkap2kap1(1,2) - 2*mHu2*adjYvYvCkap2kap3(3,2) - 2*mHu2*adjYvYvCkap3kap1(1,& 
& 3) - 2*mHu2*adjYvYvCkap3kap2(2,3) - 4*adjYvYvCmv2Ckap1kap2(2,1) - 4*adjYvYvCmv2Ckap1kap3(3,& 
& 1) - 4*adjYvYvCmv2Ckap2kap1(1,2) - 4*adjYvYvCmv2Ckap2kap3(3,2) - 4*adjYvYvCmv2Ckap3kap1(1,& 
& 3) - 4*adjYvYvCmv2Ckap3kap2(2,3) + (g1p2*(9*g2p2*M2 - 8._dp*(TradjYuTu))*Conjg(M1))/5._dp  
betamHu22 =  betamHu22- 32*g3p2*TradjYuTu*Conjg(M3) + (9*g1p2*g2p2*M1*Conjg(M2))/5._dp + 4*adjYvmlHd2(1)*Tpkap1adjkap1lam(1) -& 
&  4*Conjg(lam(1))*Tpkap1adjkap1mv2lam(1) + 4*adjYvmlHd2(1)*Tpkap1adjkap2lam(2) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap2mv2lam(2) + 4*adjYvmlHd2(1)*Tpkap1adjkap3lam(3) -        & 
&  4*Conjg(lam(1))*Tpkap1adjkap3mv2lam(3) - 4*mHd2*Conjg(lam(1))*Tpkap1Ckap1lam(1) -     & 
&  2*mHu2*Conjg(lam(1))*Tpkap1Ckap1lam(1) - 4*mHd2*Conjg(lam(1))*Tpkap1Ckap2lam(2) -     & 
&  2*mHu2*Conjg(lam(1))*Tpkap1Ckap2lam(2) - 4*mHd2*Conjg(lam(1))*Tpkap1Ckap3lam(3) -     & 
&  2*mHu2*Conjg(lam(1))*Tpkap1Ckap3lam(3) + 4*adjYvmlHd2(2)*Tpkap2adjkap1lam(1) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap1mv2lam(1) + 4*adjYvmlHd2(2)*Tpkap2adjkap2lam(2) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap2mv2lam(2) + 4*adjYvmlHd2(2)*Tpkap2adjkap3lam(3) -        & 
&  4*Conjg(lam(2))*Tpkap2adjkap3mv2lam(3) - 4*mHd2*Conjg(lam(2))*Tpkap2Ckap1lam(1) -     & 
&  2*mHu2*Conjg(lam(2))*Tpkap2Ckap1lam(1) - 4*mHd2*Conjg(lam(2))*Tpkap2Ckap2lam(2) -     & 
&  2*mHu2*Conjg(lam(2))*Tpkap2Ckap2lam(2) - 4*mHd2*Conjg(lam(2))*Tpkap2Ckap3lam(3) -     & 
&  2*mHu2*Conjg(lam(2))*Tpkap2Ckap3lam(3) + 4*adjYvmlHd2(3)*Tpkap3adjkap1lam(1) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap1mv2lam(1) + 4*adjYvmlHd2(3)*Tpkap3adjkap2lam(2) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap2mv2lam(2) + 4*adjYvmlHd2(3)*Tpkap3adjkap3lam(3) -        & 
&  4*Conjg(lam(3))*Tpkap3adjkap3mv2lam(3) - 4*mHd2*Conjg(lam(3))*Tpkap3Ckap1lam(1) -     & 
&  2*mHu2*Conjg(lam(3))*Tpkap3Ckap1lam(1) - 4*mHd2*Conjg(lam(3))*Tpkap3Ckap2lam(2) -     & 
&  2*mHu2*Conjg(lam(3))*Tpkap3Ckap2lam(2) - 4*mHd2*Conjg(lam(3))*Tpkap3Ckap3lam(3) -     & 
&  2*mHu2*Conjg(lam(3))*Tpkap3Ckap3lam(3) - 4*Conjg(Tlam(1))*TpTk1adjkap1lam(1) -        & 
&  4*Conjg(Tlam(1))*TpTk1adjkap2lam(2) - 4*Conjg(Tlam(1))*TpTk1adjkap3lam(3) -           & 
&  4*Conjg(Tlam(2))*TpTk2adjkap1lam(1) - 4*Conjg(Tlam(2))*TpTk2adjkap2lam(2)  
betamHu22 =  betamHu22- 4*Conjg(Tlam(2))*TpTk2adjkap3lam(3) - 4*Conjg(Tlam(3))*TpTk3adjkap1lam(1) -           & 
&  4*Conjg(Tlam(3))*TpTk3adjkap2lam(2) - 4*Conjg(Tlam(3))*TpTk3adjkap3lam(3) -           & 
&  4*TrCkap1Tpkap2*TpTvCTv(1,2) - 4*TrCkap1Tpkap3*TpTvCTv(1,3) - 4*TrCkap2Tpkap1*TpTvCTv(2,& 
& 1) - 4*TrCkap2Tpkap3*TpTvCTv(2,3) - 4*TrCkap3Tpkap1*TpTvCTv(3,1) - 4*TrCkap3Tpkap2*TpTvCTv(3,& 
& 2) + 6*g2p4*Tr2(2) + (6*g1p2*Tr2U1(1,1))/5._dp + 4*g1*sqrt3ov5*Tr3(1) - 2*mHu2*adjkap1kap1Clam(1)*lam(1) -& 
&  2*mHu2*adjkap1kap2Clam(2)*lam(1) - 2*mHu2*adjkap1kap3Clam(3)*lam(1) - 4*adjkap1Tpkap1Cmv2Clam(1)*lam(1) -& 
&  4*adjkap1Tpkap2Cmv2Clam(2)*lam(1) - 4*adjkap1Tpkap3Cmv2Clam(3)*lam(1) -               & 
&  8*Trmv2kap1adjkap1*Conjg(lam(1))*lam(1) - 8*Trmv2kap2adjkap1*Conjg(lam(2))*lam(1) -   & 
&  8*Trmv2kap3adjkap1*Conjg(lam(3))*lam(1) - 2*mHu2*adjkap2kap1Clam(1)*lam(2) -          & 
&  2*mHu2*adjkap2kap2Clam(2)*lam(2) - 2*mHu2*adjkap2kap3Clam(3)*lam(2) - 4*adjkap2Tpkap1Cmv2Clam(1)*lam(2) -& 
&  4*adjkap2Tpkap2Cmv2Clam(2)*lam(2) - 4*adjkap2Tpkap3Cmv2Clam(3)*lam(2) -               & 
&  8*Trmv2kap1adjkap2*Conjg(lam(1))*lam(2) - 8*Trmv2kap2adjkap2*Conjg(lam(2))*lam(2) -   & 
&  8*Trmv2kap3adjkap2*Conjg(lam(3))*lam(2) - 2*mHu2*adjkap3kap1Clam(1)*lam(3) -          & 
&  2*mHu2*adjkap3kap2Clam(2)*lam(3) - 2*mHu2*adjkap3kap3Clam(3)*lam(3) - 4*adjkap3Tpkap1Cmv2Clam(1)*lam(3) -& 
&  4*adjkap3Tpkap2Cmv2Clam(2)*lam(3) - 4*adjkap3Tpkap3Cmv2Clam(3)*lam(3) -               & 
&  8*Trmv2kap1adjkap3*Conjg(lam(1))*lam(3) - 8*Trmv2kap2adjkap3*Conjg(lam(2))*lam(3) -   & 
&  8*Trmv2kap3adjkap3*Conjg(lam(3))*lam(3) - 4*TrCTk1Tpkap1*Conjg(lam(1))*Tlam(1) -      & 
&  4*TrCTk1Tpkap2*Conjg(lam(2))*Tlam(1) - 4*TrCTk1Tpkap3*Conjg(lam(3))*Tlam(1) -         & 
&  4*TrCkap1Tpkap1*Conjg(Tlam(1))*Tlam(1) - 4*TrCkap1Tpkap2*Conjg(Tlam(2))*Tlam(1) -     & 
&  4*TrCkap1Tpkap3*Conjg(Tlam(3))*Tlam(1) - 4*TrCTk2Tpkap1*Conjg(lam(1))*Tlam(2) -       & 
&  4*TrCTk2Tpkap2*Conjg(lam(2))*Tlam(2) - 4*TrCTk2Tpkap3*Conjg(lam(3))*Tlam(2)  
betamHu22 =  betamHu22- 4*TrCkap2Tpkap1*Conjg(Tlam(1))*Tlam(2) - 4*TrCkap2Tpkap2*Conjg(Tlam(2))*Tlam(2) -     & 
&  4*TrCkap2Tpkap3*Conjg(Tlam(3))*Tlam(2) - 4*TrCTk3Tpkap1*Conjg(lam(1))*Tlam(3) -       & 
&  4*TrCTk3Tpkap2*Conjg(lam(2))*Tlam(3) - 4*TrCTk3Tpkap3*Conjg(lam(3))*Tlam(3) -         & 
&  4*TrCkap3Tpkap1*Conjg(Tlam(1))*Tlam(3) - 4*TrCkap3Tpkap2*Conjg(Tlam(2))*Tlam(3) -     & 
&  4*TrCkap3Tpkap3*Conjg(Tlam(3))*Tlam(3)

 
DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 ) 

 
Else 
DmHu2 = oo16pi2* betamHu21 
End If 
 
 
!-------------------- 
! md2 
!-------------------- 
 
betamd21  = 2*(md2YdadjYd + 2._dp*(TdadjTd) + 2*mHd2*YdadjYd + YdadjYdmd2 +           & 
&  2._dp*(Ydmq2adjYd)) - (2*id3R*(4*AbsM1*g1p2 + 80*AbsM3*g3p2 - 15*g1*ooSqrt15*Tr1(1)))/15._dp

 
 
If (TwoLoopRGE) Then 
betamd22 = md2YdadjYd*(2._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 2*(SPlamxxClam + 3._dp*(TrYdadjYd) +  & 
&  TrYeadjYe)) - (2*(5._dp*(md2YdadjYdYdadjYd) + 5._dp*(md2YdadjYuYuadjYd) -             & 
&  2*g1p2*TdadjTd - 30*g2p2*TdadjTd + 10*SPlamxxClam*TdadjTd + 10._dp*(TdadjTdYdadjYd) + & 
&  10._dp*(TdadjTuYuadjYd) + 10*SPlamxxCTlam*TdadjYd + 10._dp*(TdadjYdYdadjTd) +         & 
&  10._dp*(TdadjYuYuadjTd) + 30*TdadjYd*TrCTdTpYd + 10*TdadjYd*TrCTeTpYe +               & 
&  30*TdadjTd*TrYdadjYd + 10*TdadjTd*TrYeadjYe + 2*g1p2*M1*YdadjTd + 30*g2p2*M2*YdadjTd +& 
&  10*SPClamxxTlam*YdadjTd + 30*TradjYdTd*YdadjTd + 10*TradjYeTe*YdadjTd +               & 
&  10._dp*(YdadjTdTdadjYd) + 10._dp*(YdadjTuTuadjYd) - 4*AbsM1*g1p2*YdadjYd -            & 
&  60*AbsM2*g2p2*YdadjYd - 2*g1p2*mHd2*YdadjYd - 30*g2p2*mHd2*YdadjYd - 10*SPlamxxadjYvmlHd2*YdadjYd +& 
&  20*mHd2*SPlamxxClam*YdadjYd + 10*mHu2*SPlamxxClam*YdadjYd + 10*SPlamxxCmv2Clam*YdadjYd +& 
&  10*SPTlamxxCTlam*YdadjYd + 30*TrCTdTpTd*YdadjYd + 10*TrCTeTpTe*YdadjYd +              & 
&  30*Trmd2YdadjYd*YdadjYd + 10*Trme2YeadjYe*YdadjYd + 10*Trml2adjYeYe*YdadjYd +         & 
&  30*Trmq2adjYdYd*YdadjYd + 60*mHd2*TrYdadjYd*YdadjYd + 20*mHd2*TrYeadjYe*YdadjYd -     & 
&  g1p2*YdadjYdmd2 - 15*g2p2*YdadjYdmd2 + 5*SPlamxxClam*YdadjYdmd2 + 15*TrYdadjYd*YdadjYdmd2 +& 
&  5*TrYeadjYe*YdadjYdmd2 + 10._dp*(YdadjYdmd2YdadjYd) + 10._dp*(YdadjYdTdadjTd) +       & 
&  20*mHd2*YdadjYdYdadjYd + 5._dp*(YdadjYdYdadjYdmd2) + 10._dp*(YdadjYdYdmq2adjYd) +     & 
&  10._dp*(YdadjYumu2YuadjYd) + 10._dp*(YdadjYuTuadjTd) + 10*mHd2*YdadjYuYuadjYd +       & 
&  10*mHu2*YdadjYuYuadjYd + 5._dp*(YdadjYuYuadjYdmd2) + 10._dp*(YdadjYuYumq2adjYd) -     & 
&  2*g1p2*Ydmq2adjYd - 30*g2p2*Ydmq2adjYd + 10*SPlamxxClam*Ydmq2adjYd + 30*TrYdadjYd*Ydmq2adjYd +& 
&  10*TrYeadjYe*Ydmq2adjYd + 10._dp*(Ydmq2adjYdYdadjYd) + 10._dp*(Ydmq2adjYuYuadjYd) +   & 
&  2*g1p2*TdadjYd*Conjg(M1) + 30*g2p2*TdadjYd*Conjg(M2)))/5._dp + (8*id3R*(303*AbsM1*g1p4 +& 
&  80*AbsM1*g1p2*g3p2 + 80*AbsM3*g1p2*g3p2 - 1200*AbsM3*g3p4 + 40*g1p2*g3p2*M3*Conjg(M1) +& 
&  40*g1p2*g3p2*M1*Conjg(M3) + 300*g3p4*Tr2(3) + 15*g1p2*Tr2U1(1,1) + 225*g1*ooSqrt15*Tr3(1)))/225._dp

 
Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 ) 

 
Else 
Dmd2 = oo16pi2* betamd21 
End If 
 
 
Call Chop(Dmd2) 

Forall(i1=1:3) Dmd2(i1,i1) =  Real(Dmd2(i1,i1),dp) 
Dmd2 = 0.5_dp*(Dmd2+ Conjg(Transpose(Dmd2)) ) 
!-------------------- 
! mu2 
!-------------------- 
 
betamu21  = 2*(mu2YuadjYu + 2._dp*(TuadjTu) + 2*mHu2*YuadjYu + YuadjYumu2 +           & 
&  2._dp*(Yumq2adjYu)) - (4*id3R*(8*AbsM1*g1p2 + 40*AbsM3*g3p2 + 15*g1*ooSqrt15*Tr1(1)))/15._dp

 
 
If (TwoLoopRGE) Then 
betamu22 = (-2*(5._dp*(mu2YuadjYdYdadjYu) + 5._dp*(mu2YuadjYuYuadjYu) + mu2YuadjYu*(g1p2 +       & 
&  5*(-3._dp*(g2p2) + SPlamxxClam + 3._dp*(TrYuadjYu) + TrYvadjYv)) + 10._dp*(TuadjTdYdadjYu) +& 
&  2*g1p2*TuadjTu - 30*g2p2*TuadjTu + 10*SPlamxxClam*TuadjTu + 30*TrYuadjYu*TuadjTu +    & 
&  10*TrYvadjYv*TuadjTu + 10._dp*(TuadjTuYuadjYu) + 10._dp*(TuadjYdYdadjTu) +            & 
&  10*SPlamxxCTlam*TuadjYu + 30*TrCTuTpYu*TuadjYu + 10*TrCTvTpYv*TuadjYu +               & 
&  10._dp*(TuadjYuYuadjTu) + 10._dp*(YuadjTdTdadjYu) - 2*g1p2*M1*YuadjTu +               & 
&  30*g2p2*M2*YuadjTu + 10*SPClamxxTlam*YuadjTu + 30*TradjYuTu*YuadjTu + 10*TradjYvTv*YuadjTu +& 
&  10._dp*(YuadjTuTuadjYu) + 10._dp*(YuadjYdmd2YdadjYu) + 10._dp*(YuadjYdTdadjTu) +      & 
&  10*mHd2*YuadjYdYdadjYu + 10*mHu2*YuadjYdYdadjYu + 5._dp*(YuadjYdYdadjYumu2) +         & 
&  10._dp*(YuadjYdYdmq2adjYu) + 4*AbsM1*g1p2*YuadjYu - 60*AbsM2*g2p2*YuadjYu +           & 
&  2*g1p2*mHu2*YuadjYu - 30*g2p2*mHu2*YuadjYu - 10*SPlamxxadjYvmlHd2*YuadjYu +           & 
&  10*mHd2*SPlamxxClam*YuadjYu + 20*mHu2*SPlamxxClam*YuadjYu + 10*SPlamxxCmv2Clam*YuadjYu +& 
&  10*SPTlamxxCTlam*YuadjYu + 30*TrCTuTpTu*YuadjYu + 10*TrCTvTpTv*YuadjYu +              & 
&  30*Trmq2adjYuYu*YuadjYu + 30*Trmu2YuadjYu*YuadjYu + 60*mHu2*TrYuadjYu*YuadjYu +       & 
&  20*mHu2*TrYvadjYv*YuadjYu + 10*TrYvadjYvCml2*YuadjYu + 10*TrYvCmv2adjYv*YuadjYu +     & 
&  g1p2*YuadjYumu2 - 15*g2p2*YuadjYumu2 + 5*SPlamxxClam*YuadjYumu2 + 15*TrYuadjYu*YuadjYumu2 +& 
&  5*TrYvadjYv*YuadjYumu2 + 10._dp*(YuadjYumu2YuadjYu) + 10._dp*(YuadjYuTuadjTu) +       & 
&  20*mHu2*YuadjYuYuadjYu + 5._dp*(YuadjYuYuadjYumu2) + 10._dp*(YuadjYuYumq2adjYu) +     & 
&  10._dp*(Yumq2adjYdYdadjYu) + 2*g1p2*Yumq2adjYu - 30*g2p2*Yumq2adjYu + 10*SPlamxxClam*Yumq2adjYu +& 
&  30*TrYuadjYu*Yumq2adjYu + 10*TrYvadjYv*Yumq2adjYu + 10._dp*(Yumq2adjYuYuadjYu) -      & 
&  2*g1p2*TuadjYu*Conjg(M1) + 30*g2p2*TuadjYu*Conjg(M2)))/5._dp + (16*id3R*(642*AbsM1*g1p4 +& 
&  160*AbsM1*g1p2*g3p2 + 160*AbsM3*g1p2*g3p2 - 600*AbsM3*g3p4 + 80*g1p2*g3p2*M3*Conjg(M1) +& 
&  80*g1p2*g3p2*M1*Conjg(M3) + 150*g3p4*Tr2(3) + 30*g1p2*Tr2U1(1,1) - 225*g1*ooSqrt15*Tr3(1)))/225._dp

 
Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 ) 

 
Else 
Dmu2 = oo16pi2* betamu21 
End If 
 
 
Call Chop(Dmu2) 

Forall(i1=1:3) Dmu2(i1,i1) =  Real(Dmu2(i1,i1),dp) 
Dmu2 = 0.5_dp*(Dmu2+ Conjg(Transpose(Dmu2)) ) 
!-------------------- 
! me2 
!-------------------- 
 
betame21  = 2*(me2YeadjYe + 2._dp*(TeadjTe) + 2*mHd2*YeadjYe + YeadjYeme2 +           & 
&  2._dp*(Yeml2adjYe)) + id3R*((-24*AbsM1*g1p2)/5._dp + 2*g1*sqrt3ov5*Tr1(1))

 
 
If (TwoLoopRGE) Then 
betame22 = (-2*(5._dp*(me2YeadjYeYeadjYe) + 5._dp*(me2YeCYvTpYvadjYe) + 6*g1p2*TeadjTe -         & 
&  30*g2p2*TeadjTe + 10*SPlamxxClam*TeadjTe + 10._dp*(TeadjTeYeadjYe) + 10*SPlamxxCTlam*TeadjYe +& 
&  10._dp*(TeadjYeYeadjTe) + 10._dp*(TeCTvTpYvadjYe) + 10._dp*(TeCYvTpYvadjTe) +         & 
&  30*TeadjYe*TrCTdTpYd + 10*TeadjYe*TrCTeTpYe + 30*TeadjTe*TrYdadjYd + 10*TeadjTe*TrYeadjYe +& 
&  me2YeadjYe*(3._dp*(g1p2) + 5*(-3._dp*(g2p2) + SPlamxxClam + 3._dp*(TrYdadjYd) +       & 
&  TrYeadjYe)) - 6*g1p2*M1*YeadjTe + 30*g2p2*M2*YeadjTe + 10*SPClamxxTlam*YeadjTe +      & 
&  30*TradjYdTd*YeadjTe + 10*TradjYeTe*YeadjTe + 10._dp*(YeadjTeTeadjYe) +               & 
&  12*AbsM1*g1p2*YeadjYe - 60*AbsM2*g2p2*YeadjYe + 6*g1p2*mHd2*YeadjYe - 30*g2p2*mHd2*YeadjYe -& 
&  10*SPlamxxadjYvmlHd2*YeadjYe + 20*mHd2*SPlamxxClam*YeadjYe + 10*mHu2*SPlamxxClam*YeadjYe +& 
&  10*SPlamxxCmv2Clam*YeadjYe + 10*SPTlamxxCTlam*YeadjYe + 30*TrCTdTpTd*YeadjYe +        & 
&  10*TrCTeTpTe*YeadjYe + 30*Trmd2YdadjYd*YeadjYe + 10*Trme2YeadjYe*YeadjYe +            & 
&  10*Trml2adjYeYe*YeadjYe + 30*Trmq2adjYdYd*YeadjYe + 60*mHd2*TrYdadjYd*YeadjYe +       & 
&  20*mHd2*TrYeadjYe*YeadjYe + 3*g1p2*YeadjYeme2 - 15*g2p2*YeadjYeme2 + 5*SPlamxxClam*YeadjYeme2 +& 
&  15*TrYdadjYd*YeadjYeme2 + 5*TrYeadjYe*YeadjYeme2 + 10._dp*(YeadjYeme2YeadjYe) +       & 
&  10._dp*(YeadjYeTeadjTe) + 20*mHd2*YeadjYeYeadjYe + 5._dp*(YeadjYeYeadjYeme2) +        & 
&  10._dp*(YeadjYeYeml2adjYe) + 10._dp*(YeCTvTpTvadjYe) + 10._dp*(YeCYvmv2TpYvadjYe) +   & 
&  10._dp*(YeCYvTpTvadjTe) + 10*mHd2*YeCYvTpYvadjYe + 10*mHu2*YeCYvTpYvadjYe +           & 
&  5._dp*(YeCYvTpYvadjYeme2) + 10._dp*(YeCYvTpYvml2adjYe) + 6*g1p2*Yeml2adjYe -          & 
&  30*g2p2*Yeml2adjYe + 10*SPlamxxClam*Yeml2adjYe + 30*TrYdadjYd*Yeml2adjYe +            & 
&  10*TrYeadjYe*Yeml2adjYe + 10._dp*(Yeml2adjYeYeadjYe) + 10._dp*(Yeml2CYvTpYvadjYe) -   & 
&  6*g1p2*TeadjYe*Conjg(M1) + 30*g2p2*TeadjYe*Conjg(M2)))/5._dp + (8*id3R*(351*AbsM1*g1p4 +& 
&  15*g1p2*Tr2U1(1,1) + 25*g1*sqrt3ov5*Tr3(1)))/25._dp

 
Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 ) 

 
Else 
Dme2 = oo16pi2* betame21 
End If 
 
 
Call Chop(Dme2) 

Forall(i1=1:3) Dme2(i1,i1) =  Real(Dme2(i1,i1),dp) 
Dme2 = 0.5_dp*(Dme2+ Conjg(Transpose(Dme2)) ) 
!-------------------- 
! mv2 
!-------------------- 
 
betamv21  = -4._dp*(Dylami1adjYvmlHd2i2) + 2._dp*(Dylami1Cmv2Clami2) + 2._dp*(Dymv2lami1Clami2)& 
&  + 4._dp*(DyTlami1CTlami2) + 2._dp*(kap1adjkap1mv2) + 8._dp*(kap1Cmv2adjkap1)          & 
&  + 2._dp*(kap2adjkap2mv2) + 8._dp*(kap2Cmv2adjkap2) + 2._dp*(kap3adjkap3mv2)           & 
&  + 8._dp*(kap3Cmv2adjkap3) + 4*Dylami1Clami2*mHd2 + 4*Dylami1Clami2*mHu2 +             & 
&  2._dp*(mv2kap1adjkap1) + 2._dp*(mv2kap2adjkap2) + 2._dp*(mv2kap3adjkap3)              & 
&  + 2._dp*(mv2TpYvCYv) + 4._dp*(Tk1adjTk1) + 4._dp*(Tk2adjTk2) + 4._dp*(Tk3adjTk3)      & 
&  + 4._dp*(TpTvCTv) + 4*mHu2*TpYvCYv + 2._dp*(TpYvCYvmv2) + 4._dp*(TpYvml2CYv)

 
 
If (TwoLoopRGE) Then 
betamv22 = -16._dp*(Dykap1Ckap1kap1i11Ckap1mv2i21) - 8._dp*(Dykap1Ckap1kap1i11Cmv2Ckap1i21) -    & 
&  16._dp*(Dykap1Ckap1kap2i11Ckap2mv2i21) - 8._dp*(Dykap1Ckap1kap2i11Cmv2Ckap2i21) -     & 
&  16._dp*(Dykap1Ckap1kap3i11Ckap3mv2i21) - 8._dp*(Dykap1Ckap1kap3i11Cmv2Ckap3i21) -     & 
&  16._dp*(Dykap1Ckap2kap1i12Ckap1mv2i21) - 8._dp*(Dykap1Ckap2kap1i12Cmv2Ckap1i21) -     & 
&  16._dp*(Dykap1Ckap2kap2i12Ckap2mv2i21) - 8._dp*(Dykap1Ckap2kap2i12Cmv2Ckap2i21) -     & 
&  16._dp*(Dykap1Ckap2kap3i12Ckap3mv2i21) - 8._dp*(Dykap1Ckap2kap3i12Cmv2Ckap3i21) -     & 
&  16._dp*(Dykap1Ckap3kap1i13Ckap1mv2i21) - 8._dp*(Dykap1Ckap3kap1i13Cmv2Ckap1i21) -     & 
&  16._dp*(Dykap1Ckap3kap2i13Ckap2mv2i21) - 8._dp*(Dykap1Ckap3kap2i13Cmv2Ckap2i21) -     & 
&  16._dp*(Dykap1Ckap3kap3i13Ckap3mv2i21) - 8._dp*(Dykap1Ckap3kap3i13Cmv2Ckap3i21) -     & 
&  16._dp*(Dykap1Clami1Ckap1mv2lami2) - 16._dp*(Dykap1Cmv2Ckap1kap1i11Ckap1i21) -        & 
&  16._dp*(Dykap1Cmv2Ckap1kap2i11Ckap2i21) - 16._dp*(Dykap1Cmv2Ckap1kap3i11Ckap3i21) -   & 
&  16._dp*(Dykap1Cmv2Ckap2kap1i12Ckap1i21) - 16._dp*(Dykap1Cmv2Ckap2kap2i12Ckap2i21) -   & 
&  16._dp*(Dykap1Cmv2Ckap2kap3i12Ckap3i21) - 16._dp*(Dykap1Cmv2Ckap3kap1i13Ckap1i21) -   & 
&  16._dp*(Dykap1Cmv2Ckap3kap2i13Ckap2i21) - 16._dp*(Dykap1Cmv2Ckap3kap3i13Ckap3i21) -   & 
&  16._dp*(Dykap1i11adjYvCml2Yvadjkap11i2) - 16._dp*(Dykap1i12adjYvCml2Yvadjkap21i2) -   & 
&  16._dp*(Dykap1i13adjYvCml2Yvadjkap31i2) - 16._dp*(Dykap2Ckap1kap1i11Ckap1mv2i22) -    & 
&  8._dp*(Dykap2Ckap1kap1i11Cmv2Ckap1i22) - 16._dp*(Dykap2Ckap1kap2i11Ckap2mv2i22) -     & 
&  8._dp*(Dykap2Ckap1kap2i11Cmv2Ckap2i22) - 16._dp*(Dykap2Ckap1kap3i11Ckap3mv2i22) -     & 
&  8._dp*(Dykap2Ckap1kap3i11Cmv2Ckap3i22) - 16._dp*(Dykap2Ckap2kap1i12Ckap1mv2i22) -     & 
&  8._dp*(Dykap2Ckap2kap1i12Cmv2Ckap1i22) - 16._dp*(Dykap2Ckap2kap2i12Ckap2mv2i22) -     & 
&  8._dp*(Dykap2Ckap2kap2i12Cmv2Ckap2i22) - 16._dp*(Dykap2Ckap2kap3i12Ckap3mv2i22)  
betamv22 =  betamv22- 8._dp*(Dykap2Ckap2kap3i12Cmv2Ckap3i22) - 16._dp*(Dykap2Ckap3kap1i13Ckap1mv2i22) -     & 
&  8._dp*(Dykap2Ckap3kap1i13Cmv2Ckap1i22) - 16._dp*(Dykap2Ckap3kap2i13Ckap2mv2i22) -     & 
&  8._dp*(Dykap2Ckap3kap2i13Cmv2Ckap2i22) - 16._dp*(Dykap2Ckap3kap3i13Ckap3mv2i22) -     & 
&  8._dp*(Dykap2Ckap3kap3i13Cmv2Ckap3i22) - 16._dp*(Dykap2Clami1Ckap2mv2lami2) -         & 
&  16._dp*(Dykap2Cmv2Ckap1kap1i11Ckap1i22) - 16._dp*(Dykap2Cmv2Ckap1kap2i11Ckap2i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap1kap3i11Ckap3i22) - 16._dp*(Dykap2Cmv2Ckap2kap1i12Ckap1i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap2kap2i12Ckap2i22) - 16._dp*(Dykap2Cmv2Ckap2kap3i12Ckap3i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap3kap1i13Ckap1i22) - 16._dp*(Dykap2Cmv2Ckap3kap2i13Ckap2i22) -   & 
&  16._dp*(Dykap2Cmv2Ckap3kap3i13Ckap3i22) - 16._dp*(Dykap2i11adjYvCml2Yvadjkap12i2) -   & 
&  16._dp*(Dykap2i12adjYvCml2Yvadjkap22i2) - 16._dp*(Dykap2i13adjYvCml2Yvadjkap32i2) -   & 
&  16._dp*(Dykap3Ckap1kap1i11Ckap1mv2i23) - 8._dp*(Dykap3Ckap1kap1i11Cmv2Ckap1i23) -     & 
&  16._dp*(Dykap3Ckap1kap2i11Ckap2mv2i23) - 8._dp*(Dykap3Ckap1kap2i11Cmv2Ckap2i23) -     & 
&  16._dp*(Dykap3Ckap1kap3i11Ckap3mv2i23) - 8._dp*(Dykap3Ckap1kap3i11Cmv2Ckap3i23) -     & 
&  16._dp*(Dykap3Ckap2kap1i12Ckap1mv2i23) - 8._dp*(Dykap3Ckap2kap1i12Cmv2Ckap1i23) -     & 
&  16._dp*(Dykap3Ckap2kap2i12Ckap2mv2i23) - 8._dp*(Dykap3Ckap2kap2i12Cmv2Ckap2i23) -     & 
&  16._dp*(Dykap3Ckap2kap3i12Ckap3mv2i23) - 8._dp*(Dykap3Ckap2kap3i12Cmv2Ckap3i23) -     & 
&  16._dp*(Dykap3Ckap3kap1i13Ckap1mv2i23) - 8._dp*(Dykap3Ckap3kap1i13Cmv2Ckap1i23) -     & 
&  16._dp*(Dykap3Ckap3kap2i13Ckap2mv2i23) - 8._dp*(Dykap3Ckap3kap2i13Cmv2Ckap2i23) -     & 
&  16._dp*(Dykap3Ckap3kap3i13Ckap3mv2i23) - 8._dp*(Dykap3Ckap3kap3i13Cmv2Ckap3i23) -     & 
&  16._dp*(Dykap3Clami1Ckap3mv2lami2) - 16._dp*(Dykap3Cmv2Ckap1kap1i11Ckap1i23) -        & 
&  16._dp*(Dykap3Cmv2Ckap1kap2i11Ckap2i23) - 16._dp*(Dykap3Cmv2Ckap1kap3i11Ckap3i23)  
betamv22 =  betamv22- 16._dp*(Dykap3Cmv2Ckap2kap1i12Ckap1i23) - 16._dp*(Dykap3Cmv2Ckap2kap2i12Ckap2i23) -   & 
&  16._dp*(Dykap3Cmv2Ckap2kap3i12Ckap3i23) - 16._dp*(Dykap3Cmv2Ckap3kap1i13Ckap1i23) -   & 
&  16._dp*(Dykap3Cmv2Ckap3kap2i13Ckap2i23) - 16._dp*(Dykap3Cmv2Ckap3kap3i13Ckap3i23) -   & 
&  16._dp*(Dykap3i11adjYvCml2Yvadjkap13i2) - 16._dp*(Dykap3i12adjYvCml2Yvadjkap23i2) -   & 
&  16._dp*(Dykap3i13adjYvCml2Yvadjkap33i2) - 4._dp*(Dylami1adjTvTvClami2) -              & 
&  4._dp*(Dylami1adjYvCml2YvClami2) - 4._dp*(Dylami1adjYvTvCTlami2) + 4._dp*(Dylami1adjYvYvadjYvmlHd2i2) -& 
&  4._dp*(Dylami1adjYvYvCmv2Clami2) - 2._dp*(Dylami1Cmv2adjYvYvClami2) - 8._dp*(Dymv2kap1Clami1Ckap1lami2) -& 
&  8._dp*(Dymv2kap1i11Ckap1kap1Ckap1i21) - 8._dp*(Dymv2kap1i11Ckap1kap2Ckap1i22) -       & 
&  8._dp*(Dymv2kap1i11Ckap1kap3Ckap1i23) - 8._dp*(Dymv2kap1i12Ckap1kap1Ckap2i21) -       & 
&  8._dp*(Dymv2kap1i12Ckap1kap2Ckap2i22) - 8._dp*(Dymv2kap1i12Ckap1kap3Ckap2i23) -       & 
&  8._dp*(Dymv2kap1i13Ckap1kap1Ckap3i21) - 8._dp*(Dymv2kap1i13Ckap1kap2Ckap3i22) -       & 
&  8._dp*(Dymv2kap1i13Ckap1kap3Ckap3i23) - 8._dp*(Dymv2kap2Clami1Ckap2lami2) -           & 
&  8._dp*(Dymv2kap2i11Ckap2kap1Ckap1i21) - 8._dp*(Dymv2kap2i11Ckap2kap2Ckap1i22) -       & 
&  8._dp*(Dymv2kap2i11Ckap2kap3Ckap1i23) - 8._dp*(Dymv2kap2i12Ckap2kap1Ckap2i21) -       & 
&  8._dp*(Dymv2kap2i12Ckap2kap2Ckap2i22) - 8._dp*(Dymv2kap2i12Ckap2kap3Ckap2i23) -       & 
&  8._dp*(Dymv2kap2i13Ckap2kap1Ckap3i21) - 8._dp*(Dymv2kap2i13Ckap2kap2Ckap3i22) -       & 
&  8._dp*(Dymv2kap2i13Ckap2kap3Ckap3i23) - 8._dp*(Dymv2kap3Clami1Ckap3lami2) -           & 
&  8._dp*(Dymv2kap3i11Ckap3kap1Ckap1i21) - 8._dp*(Dymv2kap3i11Ckap3kap2Ckap1i22) -       & 
&  8._dp*(Dymv2kap3i11Ckap3kap3Ckap1i23) - 8._dp*(Dymv2kap3i12Ckap3kap1Ckap2i21) -       & 
&  8._dp*(Dymv2kap3i12Ckap3kap2Ckap2i22) - 8._dp*(Dymv2kap3i12Ckap3kap3Ckap2i23) -       & 
&  8._dp*(Dymv2kap3i13Ckap3kap1Ckap3i21) - 8._dp*(Dymv2kap3i13Ckap3kap2Ckap3i22)  
betamv22 =  betamv22- 8._dp*(Dymv2kap3i13Ckap3kap3Ckap3i23) - 2._dp*(Dymv2lami1adjYvYvClami2) -             & 
&  2._dp*(Dymv2TpYvCYvlami1Clami2) - 16._dp*(DyTk1Clami1CTk1lami2) - 16._dp*(DyTk1i11Ckap1TpYvCTvi21) -& 
&  16._dp*(DyTk1i12Ckap2TpYvCTvi21) - 16._dp*(DyTk1i13Ckap3TpYvCTvi21) - 16._dp*(DyTk2Clami1CTk2lami2) -& 
&  16._dp*(DyTk2i11Ckap1TpYvCTvi22) - 16._dp*(DyTk2i12Ckap2TpYvCTvi22) - 16._dp*(DyTk2i13Ckap3TpYvCTvi22) -& 
&  16._dp*(DyTk3Clami1CTk3lami2) - 16._dp*(DyTk3i11Ckap1TpYvCTvi23) - 16._dp*(DyTk3i12Ckap2TpYvCTvi23) -& 
&  16._dp*(DyTk3i13Ckap3TpYvCTvi23) - 4._dp*(DyTlami1adjTvYvClami2) - 4._dp*(DyTlami1adjYvYvCTlami2) -& 
&  4._dp*(DyTpTvCTvlami1Clami2) - 4._dp*(DyTpTvCYvlami1CTlami2) - 4._dp*(DyTpYvCTvTlami1Clami2) +& 
&  4._dp*(DyTpYvCYvlami1adjYvmlHd2i2) - 2._dp*(DyTpYvCYvlami1Cmv2Clami2) -               & 
&  4._dp*(DyTpYvCYvmv2lami1Clami2) - 4._dp*(DyTpYvCYvTlami1CTlami2) - 16._dp*(DyTpYvCYvTpkap11i1Ckap1mv2i21) -& 
&  16._dp*(DyTpYvCYvTpkap12i1Ckap2mv2i21) - 16._dp*(DyTpYvCYvTpkap13i1Ckap3mv2i21) -     & 
&  16._dp*(DyTpYvCYvTpkap21i1Ckap1mv2i22) - 16._dp*(DyTpYvCYvTpkap22i1Ckap2mv2i22) -     & 
&  16._dp*(DyTpYvCYvTpkap23i1Ckap3mv2i22) - 16._dp*(DyTpYvCYvTpkap31i1Ckap1mv2i23) -     & 
&  16._dp*(DyTpYvCYvTpkap32i1Ckap2mv2i23) - 16._dp*(DyTpYvCYvTpkap33i1Ckap3mv2i23) -     & 
&  4._dp*(DyTpYvml2CYvlami1Clami2) - (12*Dylami1adjYvmlHd2i2*g1p2)/5._dp +               & 
&  (24*AbsM1*Dylami1Clami2*g1p2)/5._dp + (6*Dylami1Cmv2Clami2*g1p2)/5._dp +              & 
&  (6*Dymv2lami1Clami2*g1p2)/5._dp + (12*DyTlami1CTlami2*g1p2)/5._dp - 12*Dylami1adjYvmlHd2i2*g2p2 +& 
&  24*AbsM2*Dylami1Clami2*g2p2 + 6*Dylami1Cmv2Clami2*g2p2 + 6*Dymv2lami1Clami2*g2p2 +    & 
&  12*DyTlami1CTlami2*g2p2 - 8._dp*(kap1adjYvYvadjkap1mv2) - 16._dp*(kap1adjYvYvCmv2adjkap1) -& 
&  16._dp*(kap1Ckap1Tpkap1Cmv2adjkap1) - 16._dp*(kap1Ckap2Tpkap2Cmv2adjkap1) -           & 
&  16._dp*(kap1Ckap3Tpkap3Cmv2adjkap1) - 16._dp*(kap1Cmv2adjYvYvadjkap1) -               & 
&  8._dp*(kap2adjYvYvadjkap2mv2) - 16._dp*(kap2adjYvYvCmv2adjkap2) - 16._dp*(kap2Ckap1Tpkap1Cmv2adjkap2)  
betamv22 =  betamv22- 16._dp*(kap2Ckap2Tpkap2Cmv2adjkap2) - 16._dp*(kap2Ckap3Tpkap3Cmv2adjkap2) -           & 
&  16._dp*(kap2Cmv2adjYvYvadjkap2) - 8._dp*(kap3adjYvYvadjkap3mv2) - 16._dp*(kap3adjYvYvCmv2adjkap3) -& 
&  16._dp*(kap3Ckap1Tpkap1Cmv2adjkap3) - 16._dp*(kap3Ckap2Tpkap2Cmv2adjkap3) -           & 
&  16._dp*(kap3Ckap3Tpkap3Cmv2adjkap3) - 16._dp*(kap3Cmv2adjYvYvadjkap3) -               & 
&  (12*Dylami1CTlami2*g1p2*M1)/5._dp - 12*Dylami1CTlami2*g2p2*M2 - 4*Dylami1adjYvYvClami2*mHd2 -& 
&  4*DyTpYvCYvlami1Clami2*mHd2 + (12*Dylami1Clami2*g1p2*mHd2)/5._dp + 12*Dylami1Clami2*g2p2*mHd2 -& 
&  16*Dykap1i11Ckap1TpYvCYvi21*mHu2 - 16*Dykap1i12Ckap2TpYvCYvi21*mHu2 - 16*Dykap1i13Ckap3TpYvCYvi21*mHu2 -& 
&  16*Dykap2i11Ckap1TpYvCYvi22*mHu2 - 16*Dykap2i12Ckap2TpYvCYvi22*mHu2 - 16*Dykap2i13Ckap3TpYvCYvi22*mHu2 -& 
&  16*Dykap3i11Ckap1TpYvCYvi23*mHu2 - 16*Dykap3i12Ckap2TpYvCYvi23*mHu2 - 16*Dykap3i13Ckap3TpYvCYvi23*mHu2 -& 
&  8*Dylami1adjYvYvClami2*mHu2 - 8*DyTpYvCYvlami1Clami2*mHu2 + (12*Dylami1Clami2*g1p2*mHu2)/5._dp +& 
&  12*Dylami1Clami2*g2p2*mHu2 - 8._dp*(mv2kap1adjYvYvadjkap1) - 8._dp*(mv2kap2adjYvYvadjkap2) -& 
&  8._dp*(mv2kap3adjYvYvadjkap3) - 2._dp*(mv2TpYvadjYeYeCYv) + (6*g1p2*mv2TpYvCYv)/5._dp +& 
&  6*g2p2*mv2TpYvCYv - 2._dp*(mv2TpYvCYvTpYvCYv) - 8*Dylami1CTlami2*SPClamxxTlam +       & 
&  8*Dylami1Clami2*SPlamxxadjYvmlHd2 + 8*Dylami1adjYvmlHd2i2*SPlamxxClam -               & 
&  4*Dylami1Cmv2Clami2*SPlamxxClam - 4*Dymv2lami1Clami2*SPlamxxClam - 8*DyTlami1CTlami2*SPlamxxClam -& 
&  16*Dylami1Clami2*mHd2*SPlamxxClam - 16*Dylami1Clami2*mHu2*SPlamxxClam -               & 
&  2*mv2TpYvCYv*SPlamxxClam - 8*Dylami1Clami2*SPlamxxCmv2Clam - 8*DyTlami1Clami2*SPlamxxCTlam -& 
&  8*Dylami1Clami2*SPTlamxxCTlam - 16._dp*(Tk1adjYvYvadjTk1) - 16._dp*(Tk2adjYvYvadjTk2) -& 
&  16._dp*(Tk3adjYvYvadjTk3) - 4._dp*(TpTvadjTeYeCYv) - 4._dp*(TpTvadjYeYeCTv) +         & 
&  (12*g1p2*TpTvCTv)/5._dp + 12*g2p2*TpTvCTv - 4*SPlamxxClam*TpTvCTv - 4._dp*(TpTvCTvTpYvCYv) -& 
&  4*SPlamxxCTlam*TpTvCYv - 4._dp*(TpTvCYvTpYvCTv) - 4._dp*(TpYvadjTeTeCYv)  
betamv22 =  betamv22- 4._dp*(TpYvadjYeme2YeCYv) - 4._dp*(TpYvadjYeTeCTv) - 4*mHd2*TpYvadjYeYeCYv -          & 
&  4*mHu2*TpYvadjYeYeCYv - 2._dp*(TpYvadjYeYeCYvmv2) - 4._dp*(TpYvadjYeYeml2CYv) -       & 
&  (12*g1p2*M1*TpYvCTv)/5._dp - 12*g2p2*M2*TpYvCTv - 4*SPClamxxTlam*TpYvCTv -            & 
&  4._dp*(TpYvCTvTpTvCYv) + (24*AbsM1*g1p2*TpYvCYv)/5._dp + 24*AbsM2*g2p2*TpYvCYv +      & 
&  (12*g1p2*mHu2*TpYvCYv)/5._dp + 12*g2p2*mHu2*TpYvCYv + 4*SPlamxxadjYvmlHd2*TpYvCYv -   & 
&  4*mHd2*SPlamxxClam*TpYvCYv - 8*mHu2*SPlamxxClam*TpYvCYv - 4*SPlamxxCmv2Clam*TpYvCYv - & 
&  4*SPTlamxxCTlam*TpYvCYv + (6*g1p2*TpYvCYvmv2)/5._dp + 6*g2p2*TpYvCYvmv2 -             & 
&  2*SPlamxxClam*TpYvCYvmv2 - 4._dp*(TpYvCYvmv2TpYvCYv) - 4._dp*(TpYvCYvTpTvCTv) -       & 
&  8*mHu2*TpYvCYvTpYvCYv - 2._dp*(TpYvCYvTpYvCYvmv2) - 4._dp*(TpYvCYvTpYvml2CYv) -       & 
&  4._dp*(TpYvml2adjYeYeCYv) + (12*g1p2*TpYvml2CYv)/5._dp + 12*g2p2*TpYvml2CYv -         & 
&  4*SPlamxxClam*TpYvml2CYv - 4._dp*(TpYvml2CYvTpYvCYv) - 16*kap1adjkap1*TradjTk1Tk1 -   & 
&  16*kap2adjkap2*TradjTk1Tk1 - 16*kap3adjkap3*TradjTk1Tk1 - 16*kap1adjkap1*TradjTk2Tk2 -& 
&  16*kap2adjkap2*TradjTk2Tk2 - 16*kap3adjkap3*TradjTk2Tk2 - 16*kap1adjkap1*TradjTk3Tk3 -& 
&  16*kap2adjkap2*TradjTk3Tk3 - 16*kap3adjkap3*TradjTk3Tk3 - 12*Dylami1CTlami2*TradjYdTd -& 
&  4*Dylami1CTlami2*TradjYeTe - 12*Dylami1CTlami2*TradjYuTu - 12*TpYvCTv*TradjYuTu -     & 
&  4*Dylami1CTlami2*TradjYvTv - 16*kap1adjTk1*TradjYvTv - 16*kap2adjTk2*TradjYvTv -      & 
&  16*kap3adjTk3*TradjYvTv - 4*TpYvCTv*TradjYvTv - 16*kap1adjTk1*TrCkap1Tk1 -            & 
&  16*kap1adjTk2*TrCkap1Tk2 - 16*kap1adjTk3*TrCkap1Tk3 - 16*DyTk1i11CTk1i21*TrCkap1Tpkap1 -& 
&  16*DyTk2i11CTk2i21*TrCkap1Tpkap1 - 16*DyTk3i11CTk3i21*TrCkap1Tpkap1 - 16*DyTk1i11CTk1i22*TrCkap1Tpkap2 -& 
&  16*DyTk2i11CTk2i22*TrCkap1Tpkap2 - 16*DyTk3i11CTk3i22*TrCkap1Tpkap2 - 16*DyTk1i11CTk1i23*TrCkap1Tpkap3 -& 
&  16*DyTk2i11CTk2i23*TrCkap1Tpkap3 - 16*DyTk3i11CTk3i23*TrCkap1Tpkap3 - 16*kap2adjTk1*TrCkap2Tk1  
betamv22 =  betamv22- 16*kap2adjTk2*TrCkap2Tk2 - 16*kap2adjTk3*TrCkap2Tk3 - 16*DyTk1i12CTk1i21*TrCkap2Tpkap1 -& 
&  16*DyTk2i12CTk2i21*TrCkap2Tpkap1 - 16*DyTk3i12CTk3i21*TrCkap2Tpkap1 - 16*DyTk1i12CTk1i22*TrCkap2Tpkap2 -& 
&  16*DyTk2i12CTk2i22*TrCkap2Tpkap2 - 16*DyTk3i12CTk3i22*TrCkap2Tpkap2 - 16*DyTk1i12CTk1i23*TrCkap2Tpkap3 -& 
&  16*DyTk2i12CTk2i23*TrCkap2Tpkap3 - 16*DyTk3i12CTk3i23*TrCkap2Tpkap3 - 16*kap3adjTk1*TrCkap3Tk1 -& 
&  16*kap3adjTk2*TrCkap3Tk2 - 16*kap3adjTk3*TrCkap3Tk3 - 16*DyTk1i13CTk1i21*TrCkap3Tpkap1 -& 
&  16*DyTk2i13CTk2i21*TrCkap3Tpkap1 - 16*DyTk3i13CTk3i21*TrCkap3Tpkap1 - 16*DyTk1i13CTk1i22*TrCkap3Tpkap2 -& 
&  16*DyTk2i13CTk2i22*TrCkap3Tpkap2 - 16*DyTk3i13CTk3i22*TrCkap3Tpkap2 - 16*DyTk1i13CTk1i23*TrCkap3Tpkap3 -& 
&  16*DyTk2i13CTk2i23*TrCkap3Tpkap3 - 16*DyTk3i13CTk3i23*TrCkap3Tpkap3 - 12*Dylami1Clami2*TrCTdTpTd -& 
&  12*DyTlami1Clami2*TrCTdTpYd - 4*Dylami1Clami2*TrCTeTpTe - 4*DyTlami1Clami2*TrCTeTpYe -& 
&  16*DyTk1i11Ckap1i21*TrCTk1Tpkap1 - 16*DyTk1i12Ckap2i21*TrCTk1Tpkap1 - 16*DyTk1i13Ckap3i21*TrCTk1Tpkap1 -& 
&  16*DyTk1i11Ckap1i22*TrCTk1Tpkap2 - 16*DyTk1i12Ckap2i22*TrCTk1Tpkap2 - 16*DyTk1i13Ckap3i22*TrCTk1Tpkap2 -& 
&  16*DyTk1i11Ckap1i23*TrCTk1Tpkap3 - 16*DyTk1i12Ckap2i23*TrCTk1Tpkap3 - 16*DyTk1i13Ckap3i23*TrCTk1Tpkap3 -& 
&  16*DyTk2i11Ckap1i21*TrCTk2Tpkap1 - 16*DyTk2i12Ckap2i21*TrCTk2Tpkap1 - 16*DyTk2i13Ckap3i21*TrCTk2Tpkap1 -& 
&  16*DyTk2i11Ckap1i22*TrCTk2Tpkap2 - 16*DyTk2i12Ckap2i22*TrCTk2Tpkap2 - 16*DyTk2i13Ckap3i22*TrCTk2Tpkap2 -& 
&  16*DyTk2i11Ckap1i23*TrCTk2Tpkap3 - 16*DyTk2i12Ckap2i23*TrCTk2Tpkap3 - 16*DyTk2i13Ckap3i23*TrCTk2Tpkap3 -& 
&  16*DyTk3i11Ckap1i21*TrCTk3Tpkap1 - 16*DyTk3i12Ckap2i21*TrCTk3Tpkap1 - 16*DyTk3i13Ckap3i21*TrCTk3Tpkap1 -& 
&  16*DyTk3i11Ckap1i22*TrCTk3Tpkap2 - 16*DyTk3i12Ckap2i22*TrCTk3Tpkap2 - 16*DyTk3i13Ckap3i22*TrCTk3Tpkap2 -& 
&  16*DyTk3i11Ckap1i23*TrCTk3Tpkap3 - 16*DyTk3i12Ckap2i23*TrCTk3Tpkap3 - 16*DyTk3i13Ckap3i23*TrCTk3Tpkap3 -& 
&  12*Dylami1Clami2*TrCTuTpTu - 12*TpYvCYv*TrCTuTpTu - 12*DyTlami1Clami2*TrCTuTpYu -     & 
&  12*TpTvCYv*TrCTuTpYu - 4*Dylami1Clami2*TrCTvTpTv - 16*kap1adjkap1*TrCTvTpTv -         & 
&  16*kap2adjkap2*TrCTvTpTv - 16*kap3adjkap3*TrCTvTpTv - 4*TpYvCYv*TrCTvTpTv  
betamv22 =  betamv22- 4*DyTlami1Clami2*TrCTvTpYv - 4*TpTvCYv*TrCTvTpYv - 12*Dylami1Clami2*Trmd2YdadjYd -    & 
&  4*Dylami1Clami2*Trme2YeadjYe - 4*Dylami1Clami2*Trml2adjYeYe - 12*Dylami1Clami2*Trmq2adjYdYd -& 
&  12*Dylami1Clami2*Trmq2adjYuYu - 12*TpYvCYv*Trmq2adjYuYu - 12*Dylami1Clami2*Trmu2YuadjYu -& 
&  12*TpYvCYv*Trmu2YuadjYu - 32*kap1adjkap1*Trmv2kap1adjkap1 - 32*kap2adjkap1*Trmv2kap1adjkap2 -& 
&  32*kap3adjkap1*Trmv2kap1adjkap3 - 32*kap1adjkap2*Trmv2kap2adjkap1 - 32*kap2adjkap2*Trmv2kap2adjkap2 -& 
&  32*kap3adjkap2*Trmv2kap2adjkap3 - 32*kap1adjkap3*Trmv2kap3adjkap1 - 32*kap2adjkap3*Trmv2kap3adjkap2 -& 
&  32*kap3adjkap3*Trmv2kap3adjkap3 + 12*Dylami1adjYvmlHd2i2*TrYdadjYd - 6*Dylami1Cmv2Clami2*TrYdadjYd -& 
&  6*Dymv2lami1Clami2*TrYdadjYd - 12*DyTlami1CTlami2*TrYdadjYd - 24*Dylami1Clami2*mHd2*TrYdadjYd -& 
&  12*Dylami1Clami2*mHu2*TrYdadjYd + 4*Dylami1adjYvmlHd2i2*TrYeadjYe - 2*Dylami1Cmv2Clami2*TrYeadjYe -& 
&  2*Dymv2lami1Clami2*TrYeadjYe - 4*DyTlami1CTlami2*TrYeadjYe - 8*Dylami1Clami2*mHd2*TrYeadjYe -& 
&  4*Dylami1Clami2*mHu2*TrYeadjYe + 12*Dylami1adjYvmlHd2i2*TrYuadjYu - 6*Dylami1Cmv2Clami2*TrYuadjYu -& 
&  6*Dymv2lami1Clami2*TrYuadjYu - 12*DyTlami1CTlami2*TrYuadjYu - 12*Dylami1Clami2*mHd2*TrYuadjYu -& 
&  24*Dylami1Clami2*mHu2*TrYuadjYu - 6*mv2TpYvCYv*TrYuadjYu - 12*TpTvCTv*TrYuadjYu -     & 
&  24*mHu2*TpYvCYv*TrYuadjYu - 6*TpYvCYvmv2*TrYuadjYu - 12*TpYvml2CYv*TrYuadjYu +        & 
&  4*Dylami1adjYvmlHd2i2*TrYvadjYv - 2*Dylami1Cmv2Clami2*TrYvadjYv - 2*Dymv2lami1Clami2*TrYvadjYv -& 
&  4*DyTlami1CTlami2*TrYvadjYv - 4*Dylami1Clami2*mHd2*TrYvadjYv - 8*Dylami1Clami2*mHu2*TrYvadjYv -& 
&  2*mv2TpYvCYv*TrYvadjYv - 4*TpTvCTv*TrYvadjYv - 8*mHu2*TpYvCYv*TrYvadjYv -             & 
&  2*TpYvCYvmv2*TrYvadjYv - 4*TpYvml2CYv*TrYvadjYv - 4*Dylami1Clami2*TrYvadjYvCml2 -     & 
&  4*TpYvCYv*TrYvadjYvCml2 - 4*Dylami1Clami2*TrYvCmv2adjYv - 4*TpYvCYv*TrYvCmv2adjYv -   & 
&  32*kap2adjkap1*adjkap1mv2kap1(2,1) - 32*kap3adjkap1*adjkap1mv2kap1(3,1) -             & 
&  32*kap2adjkap2*adjkap1mv2kap2(2,1) - 32*kap3adjkap2*adjkap1mv2kap2(3,1)  
betamv22 =  betamv22- 32*kap2adjkap3*adjkap1mv2kap3(2,1) - 32*kap3adjkap3*adjkap1mv2kap3(3,1) -             & 
&  32*kap1adjkap1*adjkap2mv2kap1(1,2) - 32*kap3adjkap1*adjkap2mv2kap1(3,2) -             & 
&  32*kap1adjkap2*adjkap2mv2kap2(1,2) - 32*kap3adjkap2*adjkap2mv2kap2(3,2) -             & 
&  32*kap1adjkap3*adjkap2mv2kap3(1,2) - 32*kap3adjkap3*adjkap2mv2kap3(3,2) -             & 
&  32*kap1adjkap1*adjkap3mv2kap1(1,3) - 32*kap2adjkap1*adjkap3mv2kap1(2,3) -             & 
&  32*kap1adjkap2*adjkap3mv2kap2(1,3) - 32*kap2adjkap2*adjkap3mv2kap2(2,3) -             & 
&  32*kap1adjkap3*adjkap3mv2kap3(1,3) - 32*kap2adjkap3*adjkap3mv2kap3(2,3) -             & 
&  16*kap1adjTk2*adjYvTv(1,2) - 16*kap1adjTk3*adjYvTv(1,3) - 16*kap2adjTk1*adjYvTv(2,    & 
& 1) - 16*kap2adjTk3*adjYvTv(2,3) - 16*kap3adjTk1*adjYvTv(3,1) - 16*kap3adjTk2*adjYvTv(3,& 
& 2) - 16*kap2adjTk1*Ckap1Tk1(2,1) - 16*kap3adjTk1*Ckap1Tk1(3,1) - 16*kap2adjTk2*Ckap1Tk2(2,& 
& 1) - 16*kap3adjTk2*Ckap1Tk2(3,1) - 16*kap2adjTk3*Ckap1Tk3(2,1) - 16*kap3adjTk3*Ckap1Tk3(3,& 
& 1) - 16*kap1adjTk1*Ckap2Tk1(1,2) - 16*kap3adjTk1*Ckap2Tk1(3,2) - 16*kap1adjTk2*Ckap2Tk2(1,& 
& 2) - 16*kap3adjTk2*Ckap2Tk2(3,2) - 16*kap1adjTk3*Ckap2Tk3(1,2) - 16*kap3adjTk3*Ckap2Tk3(3,& 
& 2) - 16*kap1adjTk1*Ckap3Tk1(1,3) - 16*kap2adjTk1*Ckap3Tk1(2,3) - 16*kap1adjTk2*Ckap3Tk2(1,& 
& 3) - 16*kap2adjTk2*Ckap3Tk2(2,3) - 16*kap1adjTk3*Ckap3Tk3(1,3) - 16*kap2adjTk3*Ckap3Tk3(2,& 
& 3) - (12*DyTlami1Clami2*g1p2*Conjg(M1))/5._dp - (12*g1p2*TpTvCYv*Conjg(M1))/5._dp -    & 
&  12*DyTlami1Clami2*g2p2*Conjg(M2) - 12*g2p2*TpTvCYv*Conjg(M2) - 16*DyTk1i11Ckap1lami2*Conjg(Tlam(1)) -& 
&  16*DyTk1i12Ckap2lami2*Conjg(Tlam(1)) - 16*DyTk1i13Ckap3lami2*Conjg(Tlam(1)) -         & 
&  16*DyTk2i11Ckap1lami2*Conjg(Tlam(2)) - 16*DyTk2i12Ckap2lami2*Conjg(Tlam(2)) -         & 
&  16*DyTk2i13Ckap3lami2*Conjg(Tlam(2)) - 16*DyTk3i11Ckap1lami2*Conjg(Tlam(3)) -         & 
&  16*DyTk3i12Ckap2lami2*Conjg(Tlam(3)) - 16*DyTk3i13Ckap3lami2*Conjg(Tlam(3))  
betamv22 =  betamv22- 16*kap2adjkap1*Tk1adjTk1(1,2) - 16*kap3adjkap1*Tk1adjTk1(1,3) - 16*kap1adjkap2*Tk1adjTk1(2,& 
& 1) - 16*kap3adjkap2*Tk1adjTk1(2,3) - 16*kap1adjkap3*Tk1adjTk1(3,1) - 16*kap2adjkap3*Tk1adjTk1(3,& 
& 2) - 16*kap2adjkap1*Tk2adjTk2(1,2) - 16*kap3adjkap1*Tk2adjTk2(1,3) - 16*kap1adjkap2*Tk2adjTk2(2,& 
& 1) - 16*kap3adjkap2*Tk2adjTk2(2,3) - 16*kap1adjkap3*Tk2adjTk2(3,1) - 16*kap2adjkap3*Tk2adjTk2(3,& 
& 2) - 16*kap2adjkap1*Tk3adjTk3(1,2) - 16*kap3adjkap1*Tk3adjTk3(1,3) - 16*kap1adjkap2*Tk3adjTk3(2,& 
& 1) - 16*kap3adjkap2*Tk3adjTk3(2,3) - 16*kap1adjkap3*Tk3adjTk3(3,1) - 16*kap2adjkap3*Tk3adjTk3(3,& 
& 2) - 16*kap2adjkap1*TpTvCTv(1,2) - 16*kap3adjkap1*TpTvCTv(1,3) - 16*kap1adjkap2*TpTvCTv(2,& 
& 1) - 16*kap3adjkap2*TpTvCTv(2,3) - 16*kap1adjkap3*TpTvCTv(3,1) - 16*kap2adjkap3*TpTvCTv(3,& 
& 2) - 16*Dykap1Clami1Ckap1mv2i21*lam(1) - 8*Dykap1Clami1Cmv2Ckap1i21*lam(1) -           & 
&  16*Dykap1Cmv2Clami1Ckap1i21*lam(1) - 16*Dykap2Clami1Ckap1mv2i22*lam(1) -              & 
&  8*Dykap2Clami1Cmv2Ckap1i22*lam(1) - 16*Dykap2Cmv2Clami1Ckap1i22*lam(1) -              & 
&  16*Dykap3Clami1Ckap1mv2i23*lam(1) - 8*Dykap3Clami1Cmv2Ckap1i23*lam(1) -               & 
&  16*Dykap3Cmv2Clami1Ckap1i23*lam(1) + 16*kap1adjkap1*adjYvmlHd2(1)*lam(1) +            & 
&  16*kap2adjkap1*adjYvmlHd2(2)*lam(1) + 16*kap3adjkap1*adjYvmlHd2(3)*lam(1) -           & 
&  16*kap1adjkap1*mHd2*Conjg(lam(1))*lam(1) - 16*kap1adjkap1*mHu2*Conjg(lam(1))*lam(1) - & 
&  16*kap2adjkap1*mHd2*Conjg(lam(2))*lam(1) - 16*kap2adjkap1*mHu2*Conjg(lam(2))*lam(1) - & 
&  16*kap3adjkap1*mHd2*Conjg(lam(3))*lam(1) - 16*kap3adjkap1*mHu2*Conjg(lam(3))*lam(1) - & 
&  16*Dykap1Clami1Ckap2mv2i21*lam(2) - 8*Dykap1Clami1Cmv2Ckap2i21*lam(2) -               & 
&  16*Dykap1Cmv2Clami1Ckap2i21*lam(2) - 16*Dykap2Clami1Ckap2mv2i22*lam(2) -              & 
&  8*Dykap2Clami1Cmv2Ckap2i22*lam(2) - 16*Dykap2Cmv2Clami1Ckap2i22*lam(2) -              & 
&  16*Dykap3Clami1Ckap2mv2i23*lam(2) - 8*Dykap3Clami1Cmv2Ckap2i23*lam(2)  
betamv22 =  betamv22- 16*Dykap3Cmv2Clami1Ckap2i23*lam(2) + 16*kap1adjkap2*adjYvmlHd2(1)*lam(2) +            & 
&  16*kap2adjkap2*adjYvmlHd2(2)*lam(2) + 16*kap3adjkap2*adjYvmlHd2(3)*lam(2) -           & 
&  16*kap1adjkap2*mHd2*Conjg(lam(1))*lam(2) - 16*kap1adjkap2*mHu2*Conjg(lam(1))*lam(2) - & 
&  16*kap2adjkap2*mHd2*Conjg(lam(2))*lam(2) - 16*kap2adjkap2*mHu2*Conjg(lam(2))*lam(2) - & 
&  16*kap3adjkap2*mHd2*Conjg(lam(3))*lam(2) - 16*kap3adjkap2*mHu2*Conjg(lam(3))*lam(2) - & 
&  16*Dykap1Clami1Ckap3mv2i21*lam(3) - 8*Dykap1Clami1Cmv2Ckap3i21*lam(3) -               & 
&  16*Dykap1Cmv2Clami1Ckap3i21*lam(3) - 16*Dykap2Clami1Ckap3mv2i22*lam(3) -              & 
&  8*Dykap2Clami1Cmv2Ckap3i22*lam(3) - 16*Dykap2Cmv2Clami1Ckap3i22*lam(3) -              & 
&  16*Dykap3Clami1Ckap3mv2i23*lam(3) - 8*Dykap3Clami1Cmv2Ckap3i23*lam(3) -               & 
&  16*Dykap3Cmv2Clami1Ckap3i23*lam(3) + 16*kap1adjkap3*adjYvmlHd2(1)*lam(3) +            & 
&  16*kap2adjkap3*adjYvmlHd2(2)*lam(3) + 16*kap3adjkap3*adjYvmlHd2(3)*lam(3) -           & 
&  16*kap1adjkap3*mHd2*Conjg(lam(1))*lam(3) - 16*kap1adjkap3*mHu2*Conjg(lam(1))*lam(3) - & 
&  16*kap2adjkap3*mHd2*Conjg(lam(2))*lam(3) - 16*kap2adjkap3*mHu2*Conjg(lam(2))*lam(3) - & 
&  16*kap3adjkap3*mHd2*Conjg(lam(3))*lam(3) - 16*kap3adjkap3*mHu2*Conjg(lam(3))*lam(3) - & 
&  16*kap1adjTk1*Conjg(lam(1))*Tlam(1) - 16*kap2adjTk1*Conjg(lam(2))*Tlam(1) -           & 
&  16*kap3adjTk1*Conjg(lam(3))*Tlam(1) - 16*kap1adjkap1*Conjg(Tlam(1))*Tlam(1) -         & 
&  16*kap2adjkap1*Conjg(Tlam(2))*Tlam(1) - 16*kap3adjkap1*Conjg(Tlam(3))*Tlam(1) -       & 
&  16*kap1adjTk2*Conjg(lam(1))*Tlam(2) - 16*kap2adjTk2*Conjg(lam(2))*Tlam(2) -           & 
&  16*kap3adjTk2*Conjg(lam(3))*Tlam(2) - 16*kap1adjkap2*Conjg(Tlam(1))*Tlam(2) -         & 
&  16*kap2adjkap2*Conjg(Tlam(2))*Tlam(2) - 16*kap3adjkap2*Conjg(Tlam(3))*Tlam(2) -       & 
&  16*kap1adjTk3*Conjg(lam(1))*Tlam(3) - 16*kap2adjTk3*Conjg(lam(2))*Tlam(3)  
betamv22 =  betamv22- 16*kap3adjTk3*Conjg(lam(3))*Tlam(3) - 16*kap1adjkap3*Conjg(Tlam(1))*Tlam(3) -         & 
&  16*kap2adjkap3*Conjg(Tlam(2))*Tlam(3) - 16*kap3adjkap3*Conjg(Tlam(3))*Tlam(3)

 
Dmv2 = oo16pi2*( betamv21 + oo16pi2 * betamv22 ) 

 
Else 
Dmv2 = oo16pi2* betamv21 
End If 
 
 
Call Chop(Dmv2) 

Forall(i1=1:3) Dmv2(i1,i1) =  Real(Dmv2(i1,i1),dp) 
Dmv2 = 0.5_dp*(Dmv2+ Conjg(Transpose(Dmv2)) ) 
!-------------------- 
! mlHd2 
!-------------------- 
 
betamlHd21  = -2._dp*(CTvTlam) - 2._dp*(CYvmv2lam) - CYvlam*(mHd2 + 2._dp*(mHu2))     & 
&  - ml2CYvlam

 
If (TwoLoopRGE) Then 
!rruiz
i1 = 1
betamlHd22 = -2._dp*(adjTeTeCYvlam) - 2._dp*(adjTeYeCYvTlam) - 2._dp*(adjYeme2YeCYvlam) -          & 
&  2._dp*(adjYeTeCTvlam) - 2._dp*(adjYeYeCTvTlam) - 2._dp*(adjYeYeCYvmv2lam) -           & 
&  2._dp*(adjYeYeml2CYvlam) + 4._dp*(CTvTpTvCYvlam) + 4._dp*(CTvTpYvCYvTlam) +           & 
&  4._dp*(CYvmv2TpYvCYvlam) + 4._dp*(CYvTk1adjTk1lam) + 4._dp*(CYvTk2adjTk2lam) +        & 
&  4._dp*(CYvTk3adjTk3lam) + 4._dp*(CYvTpTvCTvlam) + 4._dp*(CYvTpYvCTvTlam) +            & 
&  4._dp*(CYvTpYvCYvmv2lam) + 4._dp*(CYvTpYvml2CYvlam) - 3*adjYeYeCYvlam*mHd2 +          & 
&  2*CYvTpYvCYvlam*mHd2 - 2*adjYeYeCYvlam*mHu2 + 8*CYvTpYvCYvlam*mHu2 - ml2adjYeYeCYvlam +& 
&  2._dp*(ml2CYvTpYvCYvlam) + 6*CTvlam*SPClamxxTlam - 6*CYvlam*SPlamxxadjYvmlHd2 +       & 
&  6*CTvTlam*SPlamxxClam + 6*CYvmv2lam*SPlamxxClam + 9*CYvlam*mHd2*SPlamxxClam +         & 
&  12*CYvlam*mHu2*SPlamxxClam + 3*ml2CYvlam*SPlamxxClam + 6*CYvlam*SPlamxxCmv2Clam +     & 
&  6*CYvTlam*SPlamxxCTlam + 6*CYvlam*SPTlamxxCTlam + 6*CTvlam*TradjYuTu + 2*CTvlam*TradjYvTv +& 
&  6*CYvlam*TrCTuTpTu + 6*CYvTlam*TrCTuTpYu + 2*CYvlam*TrCTvTpTv + 2*CYvTlam*TrCTvTpYv + & 
&  6*CYvlam*Trmq2adjYuYu + 6*CYvlam*Trmu2YuadjYu + 6*CTvTlam*TrYuadjYu + 6*CYvmv2lam*TrYuadjYu +& 
&  3*CYvlam*mHd2*TrYuadjYu + 12*CYvlam*mHu2*TrYuadjYu + 3*ml2CYvlam*TrYuadjYu +          & 
&  2*CTvTlam*TrYvadjYv + 2*CYvmv2lam*TrYvadjYv + CYvlam*mHd2*TrYvadjYv + 4*CYvlam*mHu2*TrYvadjYv +& 
&  ml2CYvlam*TrYvadjYv + 2*CYvlam*TrYvadjYvCml2 + 2*CYvlam*TrYvCmv2adjYv +               & 
&  4*Conjg(Yv(i1,1))*Tpkap1adjkap1mv2lam(1) + 4*Conjg(Yv(i1,1))*Tpkap1adjkap2mv2lam(2) + & 
&  4*Conjg(Yv(i1,1))*Tpkap1adjkap3mv2lam(3) + 4*mHu2*Conjg(Yv(i1,1))*Tpkap1Ckap1lam(1) + & 
&  2*ml2CYv(i1,1)*Tpkap1Ckap1lam(1) + 4*mHu2*Conjg(Yv(i1,1))*Tpkap1Ckap2lam(2) +         & 
&  2*ml2CYv(i1,1)*Tpkap1Ckap2lam(2) + 4*mHu2*Conjg(Yv(i1,1))*Tpkap1Ckap3lam(3) +         & 
&  2*ml2CYv(i1,1)*Tpkap1Ckap3lam(3) + 4*Conjg(Yv(i1,2))*Tpkap2adjkap1mv2lam(1)  
i1 = 2
betamlHd22 =  betamlHd22+ 4*Conjg(Yv(i1,2))*Tpkap2adjkap2mv2lam(2) + 4*Conjg(Yv(i1,2))*Tpkap2adjkap3mv2lam(3) + & 
&  4*mHu2*Conjg(Yv(i1,2))*Tpkap2Ckap1lam(1) + 2*ml2CYv(i1,2)*Tpkap2Ckap1lam(1) +         & 
&  4*mHu2*Conjg(Yv(i1,2))*Tpkap2Ckap2lam(2) + 2*ml2CYv(i1,2)*Tpkap2Ckap2lam(2) +         & 
&  4*mHu2*Conjg(Yv(i1,2))*Tpkap2Ckap3lam(3) + 2*ml2CYv(i1,2)*Tpkap2Ckap3lam(3) +         & 
&  4*Conjg(Yv(i1,3))*Tpkap3adjkap1mv2lam(1) + 4*Conjg(Yv(i1,3))*Tpkap3adjkap2mv2lam(2) + & 
&  4*Conjg(Yv(i1,3))*Tpkap3adjkap3mv2lam(3) + 4*mHu2*Conjg(Yv(i1,3))*Tpkap3Ckap1lam(1) + & 
&  2*ml2CYv(i1,3)*Tpkap3Ckap1lam(1) + 4*mHu2*Conjg(Yv(i1,3))*Tpkap3Ckap2lam(2) +         & 
&  2*ml2CYv(i1,3)*Tpkap3Ckap2lam(2) + 4*mHu2*Conjg(Yv(i1,3))*Tpkap3Ckap3lam(3) +         & 
&  2*ml2CYv(i1,3)*Tpkap3Ckap3lam(3) + 4*Conjg(Tv(i1,1))*TpTk1adjkap1lam(1) +             & 
&  4*Conjg(Tv(i1,1))*TpTk1adjkap2lam(2) + 4*Conjg(Tv(i1,1))*TpTk1adjkap3lam(3) +         & 
&  4*Conjg(Tv(i1,2))*TpTk2adjkap1lam(1) + 4*Conjg(Tv(i1,2))*TpTk2adjkap2lam(2) +         & 
&  4*Conjg(Tv(i1,2))*TpTk2adjkap3lam(3) + 4*Conjg(Tv(i1,3))*TpTk3adjkap1lam(1) +         & 
&  4*Conjg(Tv(i1,3))*TpTk3adjkap2lam(2) + 4*Conjg(Tv(i1,3))*TpTk3adjkap3lam(3) +         & 
&  8*Trmv2kap1adjkap1*Conjg(Yv(i1,1))*lam(1) + 8*Trmv2kap2adjkap1*Conjg(Yv(i1,           & 
& 2))*lam(1) + 8*Trmv2kap3adjkap1*Conjg(Yv(i1,3))*lam(1) + 4*CYvmv2kap1Ckap1(i1,         & 
& 1)*lam(1) + 4*CYvmv2kap2Ckap1(i1,2)*lam(1) + 4*CYvmv2kap3Ckap1(i1,3)*lam(1) +          & 
&  2*mHd2*CYvTpkap1Ckap1(i1,1)*lam(1) + 2*mHd2*CYvTpkap2Ckap1(i1,2)*lam(1) +             & 
&  2*mHd2*CYvTpkap3Ckap1(i1,3)*lam(1) + 8*Trmv2kap1adjkap2*Conjg(Yv(i1,1))*lam(2) +      & 
&  8*Trmv2kap2adjkap2*Conjg(Yv(i1,2))*lam(2) + 8*Trmv2kap3adjkap2*Conjg(Yv(i1,           & 
& 3))*lam(2) + 4*CYvmv2kap1Ckap2(i1,1)*lam(2) + 4*CYvmv2kap2Ckap2(i1,2)*lam(2) +         & 
&  4*CYvmv2kap3Ckap2(i1,3)*lam(2) + 2*mHd2*CYvTpkap1Ckap2(i1,1)*lam(2) + 2*mHd2*CYvTpkap2Ckap2(i1,& 
& 2)*lam(2) + 2*mHd2*CYvTpkap3Ckap2(i1,3)*lam(2) + 8*Trmv2kap1adjkap3*Conjg(Yv(i1,       & 
& 1))*lam(3) + 8*Trmv2kap2adjkap3*Conjg(Yv(i1,2))*lam(3) + 8*Trmv2kap3adjkap3*Conjg(Yv(i1,& 
& 3))*lam(3) + 4*CYvmv2kap1Ckap3(i1,1)*lam(3) + 4*CYvmv2kap2Ckap3(i1,2)*lam(3)  
i1 = 3
betamlHd22 =  betamlHd22+ 4*CYvmv2kap3Ckap3(i1,3)*lam(3) + 2*mHd2*CYvTpkap1Ckap3(i1,1)*lam(3) + 2*mHd2*CYvTpkap2Ckap3(i1,& 
& 2)*lam(3) + 2*mHd2*CYvTpkap3Ckap3(i1,3)*lam(3) + 4*TrCTk1Tpkap1*Conjg(Yv(i1,           & 
& 1))*Tlam(1) + 4*TrCTk1Tpkap2*Conjg(Yv(i1,2))*Tlam(1) + 4*TrCTk1Tpkap3*Conjg(Yv(i1,     & 
& 3))*Tlam(1) + 4*TrCkap1Tpkap1*Conjg(Tv(i1,1))*Tlam(1) + 4*TrCkap1Tpkap2*Conjg(Tv(i1,   & 
& 2))*Tlam(1) + 4*TrCkap1Tpkap3*Conjg(Tv(i1,3))*Tlam(1) + 4*TrCTk2Tpkap1*Conjg(Yv(i1,    & 
& 1))*Tlam(2) + 4*TrCTk2Tpkap2*Conjg(Yv(i1,2))*Tlam(2) + 4*TrCTk2Tpkap3*Conjg(Yv(i1,     & 
& 3))*Tlam(2) + 4*TrCkap2Tpkap1*Conjg(Tv(i1,1))*Tlam(2) + 4*TrCkap2Tpkap2*Conjg(Tv(i1,   & 
& 2))*Tlam(2) + 4*TrCkap2Tpkap3*Conjg(Tv(i1,3))*Tlam(2) + 4*TrCTk3Tpkap1*Conjg(Yv(i1,    & 
& 1))*Tlam(3) + 4*TrCTk3Tpkap2*Conjg(Yv(i1,2))*Tlam(3) + 4*TrCTk3Tpkap3*Conjg(Yv(i1,     & 
& 3))*Tlam(3) + 4*TrCkap3Tpkap1*Conjg(Tv(i1,1))*Tlam(3) + 4*TrCkap3Tpkap2*Conjg(Tv(i1,   & 
& 2))*Tlam(3) + 4*TrCkap3Tpkap3*Conjg(Tv(i1,3))*Tlam(3)

 
DmlHd2 = oo16pi2*( betamlHd21 + oo16pi2 * betamlHd22 ) 

 
Else 
DmlHd2 = oo16pi2* betamlHd21 
End If 
 
 
!-------------------- 
! M1 
!-------------------- 
 
betaM11  = (66*g1p2*M1)/5._dp

 
 
If (TwoLoopRGE) Then 
betaM12 = (2*(398*g1p4*M1 + 5*g1p2*(27*g2p2*(M1 + M2) + 2*(44*g3p2*(M1 + M3) + 3._dp*(SPClamxxTlam) -& 
&  3*M1*SPlamxxClam + 7._dp*(TradjYdTd) + 9._dp*(TradjYeTe) + 13._dp*(TradjYuTu) +       & 
&  3._dp*(TradjYvTv) - 7*M1*TrYdadjYd - 9*M1*TrYeadjYe - 13*M1*TrYuadjYu -               & 
&  3*M1*TrYvadjYv))))/25._dp

 
DM1 = oo16pi2*( betaM11 + oo16pi2 * betaM12 ) 

 
Else 
DM1 = oo16pi2* betaM11 
End If 
 
 
Call Chop(DM1) 

!-------------------- 
! M2 
!-------------------- 
 
betaM21  = 2*g2p2*M2

 
 
If (TwoLoopRGE) Then 
betaM22 = (18*g1p2*g2p2*(M1 + M2))/5._dp + 4*(25*g2p4*M2 + g2p2*(12*g3p2*(M3 + M2) +            & 
&  SPClamxxTlam - M2*SPlamxxClam + 3._dp*(TradjYdTd) + TradjYeTe + 3._dp*(TradjYuTu) +   & 
&  TradjYvTv - 3*M2*TrYdadjYd - M2*TrYeadjYe - 3*M2*TrYuadjYu - M2*TrYvadjYv))

 
DM2 = oo16pi2*( betaM21 + oo16pi2 * betaM22 ) 

 
Else 
DM2 = oo16pi2* betaM21 
End If 
 
 
Call Chop(DM2) 

!-------------------- 
! M3 
!-------------------- 
 
betaM31  = -6*g3p2*M3

 
 
If (TwoLoopRGE) Then 
betaM32 = (22*g1p2*g3p2*(M1 + M3))/5._dp + 2*(9*g2p2*g3p2*(M3 + M2) + 4*(7*g3p4*M3 +            & 
&  g3p2*(TradjYdTd + TradjYuTu - M3*(TrYdadjYd + TrYuadjYu))))

 
DM3 = oo16pi2*( betaM31 + oo16pi2 * betaM32 ) 

 
Else 
DM3 = oo16pi2* betaM31 
End If 
 
 
Call Chop(DM3) 

!-------------------- 
! vd 
!-------------------- 
 
betavd1  = SPlamxxadjYvvL + (vd*(3*g1p2*(1 + Xi) + 5*(-4*(SPlamxxClam +               & 
&  3._dp*(TrYdadjYd) + TrYeadjYe) + 3*g2p2*(1 + Xi))))/20._dp

 
 
If (TwoLoopRGE) Then 
betavd2 = -2._dp*(SPadjkap1lamxxTpkap1adjYvvL) - 2._dp*(SPadjkap2lamxxTpkap2adjYvvL) -          & 
&  2._dp*(SPadjkap3lamxxTpkap3adjYvvL) - 2._dp*(SPlamxxadjYvYvadjYvvL) - 3*SPlamxxadjYvvL*SPlamxxClam +& 
&  SPvLxxadjYeYeCYvlam - 3*SPlamxxadjYvvL*TrYuadjYu - SPlamxxadjYvvL*TrYvadjYv -         & 
&  (207*g1p4*vd)/200._dp - (9*g1p2*g2p2*vd)/20._dp - 3*g2p4*vd + 2*SPlamxxadjYvYvClam*vd +& 
&  3*SPlamxxClamp2*vd + 2*SPTpkap1Clamxxadjkap1lam*vd + 2*SPTpkap2Clamxxadjkap2lam*vd +  & 
&  2*SPTpkap3Clamxxadjkap3lam*vd + (2*g1p2*TrYdadjYd*vd)/5._dp - 16*g3p2*TrYdadjYd*vd +  & 
&  9*TrYdadjYdYdadjYd*vd + 3*TrYdadjYuYuadjYd*vd - (6*g1p2*TrYeadjYe*vd)/5._dp +         & 
&  3*TrYeadjYeYeadjYe*vd + 3*SPlamxxClam*TrYuadjYu*vd + SPlamxxClam*TrYvadjYv*vd +       & 
&  TrYvadjYvTpYeCYe*vd - (9*g1p4*vd*Xi)/400._dp - (9*g1p2*g2p2*vd*Xi)/40._dp +           & 
&  (35*g2p4*vd*Xi)/16._dp - (3*g1p2*SPlamxxClam*vd*Xi)/10._dp - (3*g2p2*SPlamxxClam*vd*Xi)/2._dp -& 
&  (9*g1p2*TrYdadjYd*vd*Xi)/10._dp - (9*g2p2*TrYdadjYd*vd*Xi)/2._dp - (3*g1p2*TrYeadjYe*vd*Xi)/10._dp -& 
&  (3*g2p2*TrYeadjYe*vd*Xi)/2._dp + (9*g1p4*vd*Xip2)/400._dp + (9*g1p2*g2p2*vd*Xip2)/40._dp -& 
&  (9*g2p4*vd*Xip2)/16._dp

 
Dvd = oo16pi2*( betavd1 + oo16pi2 * betavd2 ) 

 
Else 
Dvd = oo16pi2* betavd1 
End If 
 
 
!-------------------- 
! vu 
!-------------------- 
 
betavu1  = (vu*(3*g1p2*(1 + Xi) + 5*(-4*(SPlamxxClam + 3._dp*(TrYuadjYu)              & 
&  + TrYvadjYv) + 3*g2p2*(1 + Xi))))/20._dp

 
 
If (TwoLoopRGE) Then 
betavu2 = -(vu*(9*g1p4*(46 + Xi - Xip2) + 5*(2*g1p2*(4*(3*(SPlamxxClam + TrYvadjYv)*Xi +        & 
&  TrYuadjYu*(8 + 9._dp*(Xi))) + 9*g2p2*(2 + Xi - Xip2)) + 5*(-8*(12._dp*(SPlamxxadjYvYvClam) +& 
&  6._dp*(SPlamxxClamp2) + 4._dp*(SPTpkap1Clamxxadjkap1lam) + 4._dp*(SPTpkap2Clamxxadjkap2lam) +& 
&  4._dp*(SPTpkap3Clamxxadjkap3lam) + 6*SPlamxxClam*TrYdadjYd + 6._dp*(TrYdadjYuYuadjYd) +& 
&  2*SPlamxxClam*TrYeadjYe - 32*g3p2*TrYuadjYu + 18._dp*(TrYuadjYuYuadjYu) +             & 
&  2._dp*(TrYvadjYvTpYeCYe) + 6._dp*(TrYvadjYvYvadjYv) + 4._dp*(TrYvCkap1Tpkap1adjYv) +  & 
&  4._dp*(TrYvCkap2Tpkap2adjYv) + 4._dp*(TrYvCkap3Tpkap3adjYv) - 3*g2p2*SPlamxxClam*Xi - & 
&  9*g2p2*TrYuadjYu*Xi - 3*g2p2*TrYvadjYv*Xi) + g2p4*(48 - 35._dp*(Xi) + 9._dp*(Xip2))))))/400._dp

 
Dvu = oo16pi2*( betavu1 + oo16pi2 * betavu2 ) 

 
Else 
Dvu = oo16pi2* betavu1 
End If 
 
 
!-------------------- 
! vL 
!-------------------- 
 
betavL1  = -1._dp*(TpYeCYevL) + (3*(g1p2 + 5._dp*(g2p2))*vL*(1 + Xi) + 20*(-          & 
& 1._dp*(YvadjYvvL) + vd*YvClam))/20._dp

 
 
If (TwoLoopRGE) Then 
betavL2 = 2._dp*(TpYeCYeTpYeCYevL) - (6*g1p2*TpYeCYevL)/5._dp + SPlamxxClam*TpYeCYevL +         & 
&  3*TpYeCYevL*TrYdadjYd + TpYeCYevL*TrYeadjYe + TpYeCYeYvClam*vd - (207*g1p4*vL)/200._dp -& 
&  (9*g1p2*g2p2*vL)/20._dp - 3*g2p4*vL - (3*adjYeYevL*(g1p2 + 5._dp*(g2p2))*Xi)/20._dp - & 
&  (3*CYvTpYvvL*(g1p2 + 5._dp*(g2p2))*Xi)/20._dp - (3*g1p2*TpYeCYevL*Xi)/20._dp -        & 
&  (3*g2p2*TpYeCYevL*Xi)/4._dp - (9*g1p4*vL*Xi)/400._dp - (9*g1p2*g2p2*vL*Xi)/40._dp +   & 
&  (35*g2p4*vL*Xi)/16._dp + (9*g1p4*vL*Xip2)/400._dp + (9*g1p2*g2p2*vL*Xip2)/40._dp -    & 
&  (9*g2p4*vL*Xip2)/16._dp + SPlamxxClam*YvadjYvvL + 3*TrYuadjYu*YvadjYvvL +             & 
&  TrYvadjYv*YvadjYvvL - (3*g1p2*Xi*YvadjYvvL)/20._dp - (3*g2p2*Xi*YvadjYvvL)/4._dp +    & 
&  2._dp*(YvadjYvYvadjYvvL) - 2*vd*YvadjYvYvClam + 2._dp*(YvCkap1Tpkap1adjYvvL) -        & 
&  2*vd*YvCkap1Tpkap1Clam + 2._dp*(YvCkap2Tpkap2adjYvvL) - 2*vd*YvCkap2Tpkap2Clam +      & 
&  2._dp*(YvCkap3Tpkap3adjYvvL) - 2*vd*YvCkap3Tpkap3Clam + 2*SPlamxxadjYvvL*YvClam -     & 
&  3*SPlamxxClam*vd*YvClam - 3*TrYuadjYu*vd*YvClam - TrYvadjYv*vd*YvClam

 
DvL = oo16pi2*( betavL1 + oo16pi2 * betavL2 ) 

 
Else 
DvL = oo16pi2* betavL1 
End If 
 
 
!-------------------- 
! vR 
!-------------------- 
 
betavR1  = -2*(kap1Ckap1vR + kap2Ckap2vR + kap3Ckap3vR + TpYvCYvvR + SPvRxxClam*lam)

 
 
If (TwoLoopRGE) Then 
betavR2 = 2._dp*(TpYvadjYeYeCYvvR) + 2*SPvRxxClam*TpYvCYvlam + 2._dp*(TpYvCYvTpYvCYvvR) -       & 
&  (6*g1p2*TpYvCYvvR)/5._dp - 6*g2p2*TpYvCYvvR + 2*SPlamxxClam*TpYvCYvvR +               & 
&  6*TpYvCYvvR*TrYuadjYu + 2*TpYvCYvvR*TrYvadjYv + 2*SPvRxxadjYvYvClam*lam -             & 
&  (6*g1p2*SPvRxxClam*lam)/5._dp - 6*g2p2*SPvRxxClam*lam + 4*SPlamxxClam*SPvRxxClam*lam +& 
&  6*SPvRxxClam*TrYdadjYd*lam + 2*SPvRxxClam*TrYeadjYe*lam + 6*SPvRxxClam*TrYuadjYu*lam +& 
&  2*SPvRxxClam*TrYvadjYv*lam + 8*adjkap1vR(1)*kap1adjkap1Tpkap1(1,i1) + 8*adjkap2vR(1)*kap1adjkap1Tpkap1(2,& 
& i1) + 8*adjkap3vR(1)*kap1adjkap1Tpkap1(3,i1) + 8*adjkap1vR(2)*kap1adjkap1Tpkap2(1,     & 
& i1) + 8*adjkap2vR(2)*kap1adjkap1Tpkap2(2,i1) + 8*adjkap3vR(2)*kap1adjkap1Tpkap2(3,     & 
& i1) + 8*adjkap1vR(3)*kap1adjkap1Tpkap3(1,i1) + 8*adjkap2vR(3)*kap1adjkap1Tpkap3(2,     & 
& i1) + 8*adjkap3vR(3)*kap1adjkap1Tpkap3(3,i1) + 8*adjkap1vR(1)*kap2adjkap2Tpkap1(1,     & 
& i1) + 8*adjkap2vR(1)*kap2adjkap2Tpkap1(2,i1) + 8*adjkap3vR(1)*kap2adjkap2Tpkap1(3,     & 
& i1) + 8*adjkap1vR(2)*kap2adjkap2Tpkap2(1,i1) + 8*adjkap2vR(2)*kap2adjkap2Tpkap2(2,     & 
& i1) + 8*adjkap3vR(2)*kap2adjkap2Tpkap2(3,i1) + 8*adjkap1vR(3)*kap2adjkap2Tpkap3(1,     & 
& i1) + 8*adjkap2vR(3)*kap2adjkap2Tpkap3(2,i1) + 8*adjkap3vR(3)*kap2adjkap2Tpkap3(3,     & 
& i1) + 8*adjkap1vR(1)*kap3adjkap3Tpkap1(1,i1) + 8*adjkap2vR(1)*kap3adjkap3Tpkap1(2,     & 
& i1) + 8*adjkap3vR(1)*kap3adjkap3Tpkap1(3,i1) + 8*adjkap1vR(2)*kap3adjkap3Tpkap2(1,     & 
& i1) + 8*adjkap2vR(2)*kap3adjkap3Tpkap2(2,i1) + 8*adjkap3vR(2)*kap3adjkap3Tpkap2(3,     & 
& i1) + 8*adjkap1vR(3)*kap3adjkap3Tpkap3(1,i1) + 8*adjkap2vR(3)*kap3adjkap3Tpkap3(2,     & 
& i1) + 8*adjkap3vR(3)*kap3adjkap3Tpkap3(3,i1) + 8*adjkap1vR(1)*TpYvCYvTpkap1(1,         & 
& i1) + 8*adjkap2vR(1)*TpYvCYvTpkap1(2,i1) + 8*adjkap3vR(1)*TpYvCYvTpkap1(3,             & 
& i1) + 8*adjkap1vR(2)*TpYvCYvTpkap2(1,i1) + 8*adjkap2vR(2)*TpYvCYvTpkap2(2,             & 
& i1) + 8*adjkap3vR(2)*TpYvCYvTpkap2(3,i1) + 8*adjkap1vR(3)*TpYvCYvTpkap3(1,             & 
& i1) + 8*adjkap2vR(3)*TpYvCYvTpkap3(2,i1) + 8*adjkap3vR(3)*TpYvCYvTpkap3(3,             & 
& i1) + 8*kap1Clam*adjkap1vR(1)*lam(1) + 8*kap2Clam*adjkap1vR(2)*lam(1) + 8*kap3Clam*adjkap1vR(3)*lam(1)  
betavR2 =  betavR2+ 8*kap1Clam*adjkap2vR(1)*lam(2) + 8*kap2Clam*adjkap2vR(2)*lam(2) + 8*kap3Clam*adjkap2vR(3)*lam(2) +& 
&  8*kap1Clam*adjkap3vR(1)*lam(3) + 8*kap2Clam*adjkap3vR(2)*lam(3) + 8*kap3Clam*adjkap3vR(3)*lam(3)

 
DvR = oo16pi2*( betavR1 + oo16pi2 * betavR2 ) 

 
Else 
DvR = oo16pi2* betavR1 
End If 
 
 
Call ParametersToG394(Dg1,Dg2,Dg3,DYd,DYe,Dlam,DYv,DYu,Dkap,DTd,DTe,DTlam,            & 
& DTv,DTu,DTk,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,Dmv2,DmlHd2,DM1,DM2,DM3,              & 
& Dvd,Dvu,DvL,DvR,f)

Iname = Iname - 1 
 
End Subroutine rge394  

End Module RGEs_munuSSM3G 
 
