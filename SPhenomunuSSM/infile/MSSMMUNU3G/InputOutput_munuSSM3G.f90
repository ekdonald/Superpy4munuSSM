! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.5.9b3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 14:40 on 26.8.2015   
! ----------------------------------------------------------------------  
 
 
Module InputOutput_munuSSM3G 
 
Use Control 
!Use Experiment 
Use Model_Data_munuSSM3G 
Use LoopFunctions 
Use SugraRuns_munuSSM3G 
 
Use EffPotFunctions 
Logical,Save::LesHouches_Format
Character(len=8),Save,Private::versionSARAH="4.5.9b3"
Integer,Private::i_cpv=0
Integer,Save,Private::in_kont(2)
Logical,Save::Add_Rparity= .False. 
Logical,Save::Write_HiggsBounds= .False. 
Character(len=40),Private::sp_info

Logical,Private::l_RP_Pythia= .False. 
Logical,Save,Private::Use_Flavour_States= .False. 
Real(dp),Save,Private::BrMin=1.e-4_dp 
Real(dp),Save,Private::SigMin=1.e-4_dp 
Character(len=60)::inputFileName,outputFileName 
Contains 
 
Subroutine LesHouches_Input(kont, Ecms, Pm, Pp, l_ISR, Fgmsb) 
 
Implicit None 
Integer, Intent(out) :: kont
Real(dp), Intent(out) :: Fgmsb, Ecms(:), Pm(:), Pp(:)
Logical, Intent(out) :: l_ISR(:)
Character(len=80) :: read_line
Integer :: i_mod=-1, i_sm=-1, i_par=-1, set_mod_par(25)=-1 &
& , i1, p_max, p_act, i_sp, i_model=-1, i_particles=-1
Real(dp) :: wert, Abs_Mu2, cosb2, cos2b, sinb2, RG0(3,3) &
 & , mat_D(3,3), R2(2,2), s12,s13,s23,c12,c13,c23
Logical :: check, calc_ferm, check_alpha(2)
Complex(dp) :: lam_vS
Logical, Save :: l_open = .True. 
 
Iname = Iname + 1
NameOfUnit(Iname) = "LesHouches_Input" 

check_alpha = .False. ! used to check consistency of alpha(mZ) calculation
in_kont = 0

Call InitializeStandardModel 
Call InitializeLoopFunctions 
 
i_mod = -1
i_sm = -1
i_par = -1
set_mod_par = -1 

ErrorLevel = -1
GenerationMixing=.False.
If (l_open) Then
   Open(ErrCan,file="Messages.out",status="unknown")
   Open(11,file="SPheno.out",status="unknown")
   l_open = .False.
End If 

Call Set_All_Parameters_0()

lam_vs = 0._dp
sp_info = " "
HighScaleModel="SARAH_Generated_Model" 
TwoLoopRGE = .True.
Fgmsb = 1.e12_dp
m32 = 1.e20_dp 
 
kont = 0
Open(99,file=inputFileName,status="old",err=200)
 
Do 
  Read(99,"(a80)",End=200,err=200) read_line 
  If (read_line(1:1).Eq."#") Cycle 
  If (read_line.Eq." ") Cycle 
  Call PutUpperCase(read_line) 
  If (read_line(1:5).Eq."BLOCK") Then  
    If (read_line(7:12).Eq."MODSEL") Then  
      kont = 0  
     Call Read_MODSEL(99,i_particles,i_model,i_cpv,kont)  
 CKMcomplex = CKM 
 If (i_cpv.Eq.0) Then 
 s12=lam_wolf 
 s23=s12**2*A_wolf 
 s13=s23*lam_wolf*Sqrt(eta_wolf**2+rho_wolf**2) 
 c12=Sqrt(1._dp-s12*s12) 
 c23=Sqrt(1._dp-s23*s23) 
 c13=Sqrt(1._dp-s13*s13) 
 CKM(1,1)=c12*c13 
 CKM(1,2)=s12*c13 
 CKM(1,3)=s13 
 CKM(2,1)=-s12*c23-c12*s23*s13 
 CKM(2,2)=c12*c23-s12*s23*s13 
 CKM(2,3)=s23*c13 
 CKM(3,1)=s12*s23-c12*c23*s13 
 CKM(3,2)=-c12*s23-s12*c23*s13 
 CKM(3,3)=c23*c13 
 End If 
    Else If (read_line(7:14).Eq."SMINPUTS") Then  
     Call Read_SMinput(99)  
    Else If (read_line(7:12).Eq."VCKMIN") Then  
     Call Read_CKM(99,i_cpv)  
    Else If (read_line(7:12).Eq."FCONST") Then  
     Call Read_FCONST(99)  
    Else If (read_line(7:11).Eq."FMASS") Then  
     Call Read_FMASS(99)  
    Else If (read_line(7:11).Eq."FLIFE") Then  
     Call Read_FLIFE(99)  
    Else If (read_line(7:17).Eq."SPHENOINPUT") Then  
     Call Read_SPhenoInput(99)  
    Else If (read_line(7:12).Eq."MINPAR") Then  
     Call Read_MINPAR(99,0,i_model,set_mod_par,kont)  
    Else If (read_line(7:14).Eq."IMMINPAR") Then  
       If (i_cpv.Lt.2) Then 
       Call Warn_CPV(i_cpv,"IMMINPAR") 
       End If 
    Call Read_MINPAR(99,1,i_model,set_mod_par,kont)  
    Else If (read_line(7:12).Eq."EXTPAR") Then  
     Call Read_EXTPAR(99,0,i_model,set_mod_par,kont)  
    Else If (read_line(7:14).Eq."IMEXTPAR") Then  
       If (i_cpv.Lt.2) Then 
       Call Warn_CPV(i_cpv,"IMEXTPAR") 
       End If 
    Call Read_EXTPAR(99,1,i_model,set_mod_par,kont)  
   Else If (read_line(7:10).Eq."YDIN") Then 
InputValueforYd= .True. 
    Call ReadMatrixC(99,3,3,YdIN,0, "YdIN",kont)

 
   Else If (read_line(7:12).Eq."IMYDIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMYd") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,YdIN,1, "YdIN",kont)

 
   Else If (read_line(7:10).Eq."YEIN") Then 
InputValueforYe= .True. 
    Call ReadMatrixC(99,3,3,YeIN,0, "YeIN",kont)

 
   Else If (read_line(7:12).Eq."IMYEIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMYe") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,YeIN,1, "YeIN",kont)

 
   Else If (read_line(7:16).Eq."NMSSMRUNIN") Then 
InputValueforlam= .True. 
    Call ReadVectorC(99,3,lamIN,0, "lamIN",kont)

   Else If (read_line(7:11).Eq."LAMIN") Then 
    InputValueforlam= .True. 
    Call ReadVectorC(99,3,lamIN,0, "lamIN",kont)

   Else If (read_line(7:13).Eq."IMLAMIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMlam") 
       Cycle 
     End If 
    Call ReadVectorC(99,3,lamIN,1, "lamIN",kont)

 
   Else If (read_line(7:10).Eq."YVIN") Then 
InputValueforYv= .True. 
    Call ReadMatrixC(99,3,3,YvIN,0, "YvIN",kont)

 
   Else If (read_line(7:12).Eq."IMYVIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMYv") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,YvIN,1, "YvIN",kont)

 
   Else If (read_line(7:10).Eq."YUIN") Then 
InputValueforYu= .True. 
    Call ReadMatrixC(99,3,3,YuIN,0, "YuIN",kont)

 
   Else If (read_line(7:12).Eq."IMYUIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMYu") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,YuIN,1, "YuIN",kont)

 
   Else If (read_line(7:16).Eq."NMSSMRUNIN") Then 
InputValueforkap= .True. 
    Call ReadTensorC(99,3,3,3,kapIN,0, "kapIN",kont)

   Else If (read_line(7:11).Eq."KAPIN") Then 
    InputValueforkap= .True. 
    Call ReadTensorC(99,3,3,3,kapIN,0, "kapIN",kont)
 
   Else If (read_line(7:13).Eq."IMKAPIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMkap") 
       Cycle 
     End If 
    Call ReadTensorC(99,3,3,3,kapIN,1, "kapIN",kont)

 
   Else If (read_line(7:10).Eq."TDIN") Then 
InputValueforTd= .True. 
    Call ReadMatrixC(99,3,3,TdIN,0, "TdIN",kont)

 
   Else If (read_line(7:12).Eq."IMTDIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMTd") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,TdIN,1, "TdIN",kont)

 
   Else If (read_line(7:10).Eq."TEIN") Then 
InputValueforTe= .True. 
    Call ReadMatrixC(99,3,3,TeIN,0, "TeIN",kont)

 
   Else If (read_line(7:12).Eq."IMTEIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMTe") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,TeIN,1, "TeIN",kont)

   Else If (read_line(7:12).Eq."TLAMIN") Then 
    InputValueforTlam= .True. 
    Call ReadVectorC(99,3,TlamIN,0, "TlamIN",kont)
 
   Else If (read_line(7:16).Eq."NMSSMRUNIN") Then 
InputValueforTlam= .True. 
    Call ReadVectorC(99,3,TlamIN,0, "TlamIN",kont)

 
   Else If (read_line(7:14).Eq."IMTLAMIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMTlam") 
       Cycle 
     End If 
    Call ReadVectorC(99,3,TlamIN,1, "TlamIN",kont)

 
   Else If (read_line(7:10).Eq."TVIN") Then 
InputValueforTv= .True. 
    Call ReadMatrixC(99,3,3,TvIN,0, "TvIN",kont)

 
   Else If (read_line(7:12).Eq."IMTVIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMTv") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,TvIN,1, "TvIN",kont)

 
   Else If (read_line(7:10).Eq."TUIN") Then 
InputValueforTu= .True. 
    Call ReadMatrixC(99,3,3,TuIN,0, "TuIN",kont)

 
   Else If (read_line(7:12).Eq."IMTUIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMTu") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,TuIN,1, "TuIN",kont)

 
   Else If (read_line(7:16).Eq."NMSSMRUNIN") Then 
InputValueforTk= .True. 
    Call ReadTensorC(99,3,3,3,TkIN,0, "TkIN",kont)

   Else If (read_line(7:10).Eq."TKIN") Then 
    InputValueforTk= .True.
    Call ReadTensorC(99,3,3,3,TkIN,0, "TkIN",kont) 

   Else If (read_line(7:12).Eq."IMTKIN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMTk") 
       Cycle 
     End If 
    Call ReadTensorC(99,3,3,3,TkIN,1, "TkIN",kont)

 
   Else If (read_line(7:12).Eq."MSQ2IN") Then 
InputValueformq2= .True. 
    Call ReadMatrixC(99,3,3,mq2IN,0, "mq2IN",kont)

 
   Else If (read_line(7:13).Eq."IMMQ2IN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMmq2") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,mq2IN,1, "mq2IN",kont)

 
   Else If (read_line(7:12).Eq."MSL2IN") Then 
InputValueforml2= .True. 
    Call ReadMatrixC(99,3,3,ml2IN,0, "ml2IN",kont)

 
   Else If (read_line(7:13).Eq."IMML2IN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMml2") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,ml2IN,1, "ml2IN",kont)

 
   Else If (read_line(7:12).Eq."MSD2IN") Then 
InputValueformd2= .True. 
    Call ReadMatrixC(99,3,3,md2IN,0, "md2IN",kont)

 
   Else If (read_line(7:13).Eq."IMMD2IN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMmd2") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,md2IN,1, "md2IN",kont)

 
   Else If (read_line(7:12).Eq."MSU2IN") Then 
InputValueformu2= .True. 
    Call ReadMatrixC(99,3,3,mu2IN,0, "mu2IN",kont)

 
   Else If (read_line(7:13).Eq."IMMU2IN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMmu2") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,mu2IN,1, "mu2IN",kont)

 
   Else If (read_line(7:12).Eq."MSE2IN") Then 
InputValueforme2= .True. 
    Call ReadMatrixC(99,3,3,me2IN,0, "me2IN",kont)

 
   Else If (read_line(7:13).Eq."IMME2IN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMme2") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,me2IN,1, "me2IN",kont)

 
   Else If (read_line(7:11).Eq."MV2IN") Then 
InputValueformv2= .True. 
    Call ReadMatrixC(99,3,3,mv2IN,0, "mv2IN",kont)

 
   Else If (read_line(7:13).Eq."IMMV2IN") Then 
     If (i_cpv.Lt.2) Then  
       Call Warn_CPV(i_cpv,"IMmv2") 
       Cycle 
     End If 
    Call ReadMatrixC(99,3,3,mv2IN,1, "mv2IN",kont)

 
   Else If (read_line(7:15).Eq."RVM2LH1IN") Then 
InputValueformlHd2= .True. 
    Call ReadVectorR(99,3,mlHd2IN, "mlHd2IN",kont)

 
   Else If (read_line(7:16).Eq."RIGHTVEVIN") Then 
InputValueforvR= .True. 
    Call ReadVectorR(99,3,vRIN, "vRIN",kont)

 
   Else If (read_line(7:15).Eq."RVSNVEVIN") Then 
InputValueforvL= .True. 
    Call ReadVectorR(99,3,vLIN, "vLIN",kont)

 
   Else If (read_line(7:13).Eq."GAUGEIN") Then 
    Call Read_GAUGEIN(99,0,i_model,set_mod_par,kont) 
 
   Else If (read_line(7:13).Eq."MSOFTIN") Then 
    Call Read_MSOFTIN(99,0,i_model,set_mod_par,kont) 
 
   Else If (read_line(7:12).Eq."HMIXIN") Then 
    Call Read_HMIXIN(99,0,i_model,set_mod_par,kont) 
 
   Else If (read_line(7:14).Eq."PHASESIN") Then 
    Call Read_PHASESIN(99,0,i_model,set_mod_par,kont) 
 
End if 
End If 
End Do 
200 Close(99) 
gmZ = gamZ * mZ
gmZ2 = gmZ**2
mW2 = mZ2 * (0.5_dp + Sqrt(0.25_dp-Alpha_Mz*pi / (sqrt2*G_F*mZ2))) / 0.987_dp
mW = Sqrt(mW2) 
mW_SM = MW 
gamW = 2.06_dp 
gamW2 = gamW**2
gmW = gamW * mW
gmW2 = gmW**2
Alpha_mZ = Alpha_MSbar(mZ, mW)
If (calc_ferm) Call CalculateRunningMasses(mf_l,mf_d,mf_u&
&,Q_light_quarks,alpha_mZ,alphas_mZ,mZ&
&,mf_l_mZ,mf_d_mZ,mf_u_mZ,kont)


Iname=Iname-1
Contains
Subroutine Read_MINPAR(io,i_c,i_model,set_mod_par,kont) 
Implicit None 
Integer,Intent(in)::io,i_c,i_model 
Integer,Intent(inout)::kont,set_mod_par(:) 
Integer::i_par 
Real(dp)::wert 
Character(len=80)::read_line 
Do 
Read(io,*,End=200) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line 
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_par,wert!,read_line 
If (i_par.Eq.1) Then 
m0= wert 
Else If (i_par.Eq.2) Then 
If (i_c.Eq.0) m12= Cmplx(wert,Aimag(m12),dp) 
If (i_c.Eq.1) m12= Cmplx(Real(m12,dp),wert,dp) 
Else If (i_par.Eq.3) Then 
TanBeta= wert 
Else If (i_par.Eq.5) Then 
If (i_c.Eq.0) Azero= Cmplx(wert,Aimag(Azero),dp) 
If (i_c.Eq.1) Azero= Cmplx(Real(Azero,dp),wert,dp) 
Else
Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block MINPAR ",i_par
If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMMINPAR ",i_par
If (i_c.Eq.0) Write(*,*) "Unknown entry for Block MINPAR ",i_par
If (i_c.Eq.1) Write(*,*) "Unknown entry for Block IMMINPAR ",i_par
Call AddError(304)
If (ErrorLevel.Eq.2) Call TerminateProgram
End If
End Do! i_par
200 Return
End Subroutine Read_MINPAR 
 
 
Subroutine Read_EXTPAR(io,i_c,i_model,set_mod_par,kont) 
Implicit None 
Integer,Intent(in)::io,i_c,i_model 
Integer,Intent(inout)::kont,set_mod_par(:) 
Integer::i_par 
Real(dp)::wert 
Character(len=80)::read_line 
Do 
Read(io,*,End=200) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line 
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_par,wert!,read_line 
If (i_par.Eq.61) Then 
If (i_c.Eq.0) LambdaInput= Cmplx(wert,Aimag(LambdaInput),dp) 
If (i_c.Eq.1) LambdaInput= Cmplx(Real(LambdaInput,dp),wert,dp) 
Else If (i_par.Eq.62) Then 
If (i_c.Eq.0) KappaInput= Cmplx(wert,Aimag(KappaInput),dp) 
If (i_c.Eq.1) KappaInput= Cmplx(Real(KappaInput,dp),wert,dp) 
Else If (i_par.Eq.63) Then 
If (i_c.Eq.0) ALambdaInput= Cmplx(wert,Aimag(ALambdaInput),dp) 
If (i_c.Eq.1) ALambdaInput= Cmplx(Real(ALambdaInput,dp),wert,dp) 
Else If (i_par.Eq.64) Then 
If (i_c.Eq.0) AKappaInput= Cmplx(wert,Aimag(AKappaInput),dp) 
If (i_c.Eq.1) AKappaInput= Cmplx(Real(AKappaInput,dp),wert,dp) 
Else If (i_par.Eq.65) Then 
If (i_c.Eq.0) vR1Input= Cmplx(wert,Aimag(vR1Input),dp) 
If (i_c.Eq.1) vR1Input= Cmplx(Real(vR1Input,dp),wert,dp) 
Else If (i_par.Eq.66) Then 
If (i_c.Eq.0) vR2Input= Cmplx(wert,Aimag(vR2Input),dp) 
If (i_c.Eq.1) vR2Input= Cmplx(Real(vR2Input,dp),wert,dp) 
Else If (i_par.Eq.67) Then 
If (i_c.Eq.0) vR3Input= Cmplx(wert,Aimag(vR3Input),dp) 
If (i_c.Eq.1) vR3Input= Cmplx(Real(vR3Input,dp),wert,dp) 
Else If (i_par.Eq.200) Then 
If (i_c.Eq.0) vL1Input= Cmplx(wert,Aimag(vL1Input),dp) 
If (i_c.Eq.1) vL1Input= Cmplx(Real(vL1Input,dp),wert,dp) 
Else If (i_par.Eq.201) Then 
If (i_c.Eq.0) vL2Input= Cmplx(wert,Aimag(vL2Input),dp) 
If (i_c.Eq.1) vL2Input= Cmplx(Real(vL2Input,dp),wert,dp) 
Else If (i_par.Eq.202) Then 
If (i_c.Eq.0) vL3Input= Cmplx(wert,Aimag(vL3Input),dp) 
If (i_c.Eq.1) vL3Input= Cmplx(Real(vL3Input,dp),wert,dp) 
Else
Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block EXTPAR ",i_par
If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMEXTPAR ",i_par
If (i_c.Eq.0) Write(*,*) "Unknown entry for Block EXTPAR ",i_par
If (i_c.Eq.1) Write(*,*) "Unknown entry for Block IMEXTPAR ",i_par
Call AddError(304)
If (ErrorLevel.Eq.2) Call TerminateProgram
End If
End Do! i_par
200 Return
End Subroutine Read_EXTPAR 
 
 
 Subroutine Read_MODSEL(io, i_particles, i_model, i_cpv, kont)
  Implicit None
   Integer, Intent(in) :: io
   Integer, Intent(out) :: i_particles, i_model, i_cpv
   Integer, Intent(inout) :: kont
    Real(dp) :: r_mod

   Integer :: i_mod, i_test, i_rp
   Character(len=80) :: read_line

   i_cpv = 0
   i_rp = 0

    Do 
     Read(io,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_test,r_mod ! ,read_line
     If (i_test.Ne.12) Then
      Backspace(io)
      Read(io,*) i_test,i_mod ! ,read_line
     End If

!      Read(io,*) i_test,i_mod,read_line
     If (i_test.Eq.1) Then
      i_particles = i_test
      i_model = i_mod
      If ((i_mod.Lt.0).Or.(i_mod.Gt.99)) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "MSSM, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      Else If (i_mod.Eq.1) Then
       HighScaleModel = "GUT"
      Else If (i_mod.Eq.0) Then
       HighScaleModel = "LOW"
      End If

     Else If (i_test.Eq.2) Then
      BoundaryCondition = i_mod
     Else If (i_test.Eq.4) Then
      If (i_mod.Eq.1) Then
       i_rp = 1

      Else If (i_mod.Ne.0) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "Unknown entry for Block MODSEL ",i_test,i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.5) Then
      i_cpv = i_mod
      If ((i_mod.Lt.0).Or.(i_mod.Gt.2)) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "Unknown entry for Block MODSEL ",i_test,i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.6) Then
      If (i_mod.Eq.0) Then
       GenerationMixing = .False.
      Else If ((i_mod.Ge.1).And.(i_mod.Le.3)) Then
       GenerationMixing = .True.
      Else
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "GenerationMixing, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

    Else If (i_test.Eq.12) Then
      Call SetRGEScale(r_mod**2)  ! set Q_EWSB

     End If
    End Do ! i_mod

  End Subroutine Read_MODSEL

  Subroutine Read_SMinput(io)
  Implicit None
   Integer, Intent(in) :: io
   
   Integer :: i_sm
   Real(dp) :: wert
   Character(len=80) :: read_line

    Do 
     Read(io,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_sm,wert,read_line

     Select Case(i_sm)
     Case(1)
      check_alpha(1) = .True.
      MZ_input = .True.
      Alpha_MZ_MS = 1._dp / wert

     Case(2) ! G_F
      G_F = wert

     Case(3) ! alpha_s(m_Z)
      alphaS_mZ = wert

     Case(4) ! m_Z
      mZ = wert
      mZ2 = mZ**2
      calc_ferm = .True.

     Case(5) ! m_b(m_b)^MSbar
      mf_d(3) = wert
      mf_d2(3) = mf_d(3)**2
      calc_ferm = .True.

     Case(6) ! m_t^pole
      mf_u(3) = wert
      mf_u2(3) = mf_u(3)**2

     Case(7) ! m_tau^pole
      mf_l(3) = wert
      mf_l2(3) = mf_l(3)**2
      calc_ferm = .True.

     Case(8) ! m_nu_3, input is in GeV
      Mf_nu(3) = wert

     Case(11) ! electron mass
      mf_l(1) = wert
      mf_l2(1) = wert**2
      calc_ferm = .True.

     Case(12) ! m_nu_1, input is in GeV
      Mf_nu(1) = wert 

     Case(13) ! muon mass
      mf_l(2) = wert
      mf_l2(2) = wert**2
      calc_ferm = .True.

     Case(14) ! m_nu_2, input is in eV, transform to GeV
      Mf_nu(2) = wert 

     Case(21) ! d-quark mass at 2 GeV
      mf_d(1) = wert
      mf_d2(1) = wert**2
      calc_ferm = .True.

     Case(22) ! u-quark mass at 2 GeV
      mf_u(1) = wert
      mf_u2(1) = wert**2
      calc_ferm = .True.

     Case(23) ! s-quark mass at 2 GeV
      mf_d(2) = wert
      mf_d2(2) = wert**2
      calc_ferm = .True.

     Case(24) ! c-quark mass at Q=m_c
      mf_u(2) = wert
      mf_u2(2) = wert**2
      calc_ferm = .True.

     Case Default
      If (output_screen) &
           & Write(*,*) "Ignoring unknown entry for Block SMINPUTS ",i_sm 
      Write(ErrCan,*) "Ignoring unknown entry for Block SMINPUTS ",i_sm 
     End Select

    End Do ! i_sm

  End Subroutine Read_SMinput

 Subroutine Read_CKM(io, i_cpv)
  Implicit None
   Integer, Intent(in) :: io, i_cpv

   Real(dp) :: s12, s13, s23, c12, c13, c23, phase
    
    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)     
     Case(1)
      lam_wolf = wert
     Case(2)
      A_wolf = wert
     Case(3)
      rho_wolf = wert
     Case(4)
      eta_wolf = wert
     Case default
     End Select

    End Do

 200   s12 = lam_wolf
    s23 = s12**2 * A_wolf
    s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2)
    If (i_cpv.Eq.0) Then
     Write(ErrCan,*) "Warning: CP violation is switched of, ignoring CKM phase."
     phase = 0._dp
    Else
     phase = Atan(eta_wolf/rho_wolf)
    End If


    c12 = Sqrt(1._dp-s12*s12)
    c23 = Sqrt(1._dp-s23*s23)
    c13 = Sqrt(1._dp-s13*s13)

    CKM(1,1) = c12 * c13
    CKM(1,2) = s12 * c13
    CKM(2,3) = s23 * c13
    CKM(3,3) = c23 * c13
    If (phase.Ne.0._dp) Then
     CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
     CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    Else
     CKM(1,3) = s13
     CKM(2,1) = -s12*c23 -c12*s23*s13
     CKM(2,2) = c12*c23 -s12*s23*s13
     CKM(3,1) = s12*s23 -c12*c23*s13
     CKM(3,2) = -c12*s23 - s12*c23*s13
    End If

  End Subroutine Read_CKM

 Subroutine Read_SPINFO(io, kont)
  Implicit None
   Integer, Intent(in) :: io
   Integer, Intent(inout) :: kont

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line

     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_sp, read_line

     If (i_sp.Eq.1) Then
      sp_info = Trim(read_line)//" "//Trim(sp_info)
     Else If (i_sp.Eq.2) Then
      sp_info = Trim(sp_info)//" "//Trim(read_line)
     Else If (i_sp.Eq.4) Then ! there is some inconsistency, exit
      kont = -306
      Call AddError(306)
      Iname = Iname - 1
      Return
     Else
      Write(ErrCan,*) "Unknown entry in BLOCK SPINFO, is ignored:"
      Write(ErrCan,*) i_sp, read_line
     End If
    End Do

   200 Return

  End Subroutine Read_SPINFO

  Subroutine Read_SPhenoInput(io)
  Implicit None
   Integer, Intent(in) :: io

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    ! This initialization is necessary for the arrar of production infos
    p_max = Size(Ecms)
    p_act = 0
    Ecms = 0._dp
    Pm = 0._dp
    Pp = 0._dp
    l_ISR = .False.
    Do 
     Read(io,*,End=200,err=200) read_line
!     Write(*,*) trim(read_line)
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*,End=200) i_par,wert,read_line
!     write(*,*) i_par,wert,trim(read_line)
     Select Case(i_par)
     Case(1)
      ErrorLevel = Int(wert)

     Case(2)
      If (Int(wert).Ne.0) Then
       SPA_convention = .True.
       Call SetRGEScale(1.e3_dp**2)
      End If

     Case(3)
      If (Int(wert).Ne.0) Then 
       External_Spectrum = .True.
       External_Higgs = .True.
      End If

     Case(4)
      If (Int(wert).Ne.0) Use_Flavour_States = .True.

     Case(5)
      If (Int(wert).Ne.0) FermionMassResummation = .False.
      
     Case(6)
!       RXiNew = wert      
       RXiNew = 1d0 
     Case(7)
       If (wert.eq.1) then
         CalculateTwoLoopHiggsMasses=.False.
       Else
         CalculateTwoLoopHiggsMasses=.True.
       End if

     Case(8)
        SELECT CASE (int(WERT))
        CASE ( 1 )
           PurelyNumericalEffPot = .true.
           CalculateMSSM2Loop = .false.
           TwoLoopMethod=1
        CASE ( 2 )
           PurelyNumericalEffPot = .false.
           CalculateMSSM2Loop = .false.
           TwoLoopMethod=2
        CASE ( 3 )
           CalculateMSSM2Loop = .false.
           TwoLoopMethod=3
        CASE ( 8 )
           CalculateMSSM2Loop = .True.
           TwoLoopMethod=8
        CASE ( 9 )
           CalculateMSSM2Loop = .True.
           TwoLoopMethod=9
        CASE DEFAULT 
           Write(*,*) "Unknown option for two-loop mass calculation"
           CalculateTwoLoopHiggsMasses=.False.
        END SELECT


 
     Case(9)
       If (wert.Ne.0) Then
        GaugelessLimit=.true.
       Else
        GaugelessLimit=.false.
       End If

     Case(400)
       hstep_pn = wert
     Case(401)
       hstep_sa = wert

     Case(10)
       If (wert.Ne.1) Then
        TwoLoopSafeMode=.false.
       Else
        TwoLoopSafeMode=.true.
       End If

     Case(11)  ! whether to calculate  branching ratios or not
      If (Int(wert).Eq.1) L_BR = .True.
      If (Int(wert).Eq.0) L_BR = .False.

     Case(12) ! minimal value such that a branching ratio is written out
      Call SetWriteMinBR(wert)

     Case(13) ! minimal value such that a branching ratio is written out
      If (wert.Eq.0) Then
           Enable3BDecaysF = .False.
           Enable3BDecaysS = .False.        
      Elseif (wert.Eq.1) Then
           Enable3BDecaysF = .True.
           Enable3BDecaysS = .False.        
      Elseif (wert.Eq.2) Then
           Enable3BDecaysS = .True.
           Enable3BDecaysF = .False.        
      Elseif (wert.Eq.3) Then
           Enable3BDecaysF = .True.
           Enable3BDecaysS = .True.        
      Else 
          Write(*,*) "Unknown option for flag 13 (three-body decays): ",wert
      End if


     Case(14) ! run SUSY couplings to scale of decaying particle
      If (wert.Eq.0) RunningCouplingsDecays = .False.    

     Case(15) ! run SUSY couplings to scale of decaying particle
      MinWidth = wert    

!      Case(21)  ! whether to calculate cross sections or not
!       If (Int(wert).Eq.1) L_CS = .True.
!       If (Int(wert).Eq.0) L_CS = .False.
! 
!      Case(22) ! cms energy
!       p_act = p_act + 1
!       ! this test is necessary to avoid a memory violation
!       If (p_act.Le.p_max) Then
!        Ecms(p_act) = wert
!       Else
!        If (output_screen) &
!            & Write(*,*) "The number of required points for the calculation"// &
!            &  " of cross sections exceeds",p_max
!        If (output_screen) &
!            & Write(*,*) "Ignoring this information"
!        If (output_screen) &
!      &  Write(*,*) "Please enlarge the corresponding arrays in the main program."
!        Write(ErrCan,*) "The number of required points for the calculation"// &
!                &   " of cross sections exceeds",p_max
!        Write(ErrCan,*) "Ignoring this information"
!        Write(ErrCan,*) &
!          &"Please enlarge the corresponding arrays in the main program."
!       End If

!      Case (23) ! polarisation of incoming e- beam
!       If (Abs(wert).Gt.1._dp) Then
!        If (output_screen) Write(*,*) &
!            & "e- beam polarisation has to between -1 and 1 and not",wert
!        If (output_screen) &
!            & Write(*,*) "using now unpolarised e- beam"
!        Write(ErrCan,*) &
!           & "e- beam polarisation has to between -1 and 1 and not",wert
!        Write(ErrCan,*) "using now unpolarised e- beam"
!        If (p_act.Le.p_max) Pm(p_act) = 0
!       Else
!        If (p_act.Le.p_max) Pm(p_act) = wert
!       End If
! 
!      Case (24) ! polarisation of incoming e+ beam
!       If (Abs(wert).Gt.1._dp) Then
!        If (output_screen) Write(*,*) &
!            & "e+ beam polarisation has to between -1 and 1 and not",wert
!        If (output_screen) &
!            & Write(*,*) "using now unpolarised e+ beam"
!        Write(ErrCan,*) &
!           & "e+ beam polarisation has to between -1 and 1 and not",wert
!        Write(ErrCan,*) "using now unpolarised e+ beam"
!        If (p_act.Le.p_max) Pp(p_act) = 0
!       Else
!        If (p_act.Le.p_max) Pp(p_act) = wert
!       End If

!      Case(25)
!       If ((wert.Eq.1._dp).And.(p_act.Le.p_max)) L_ISR(p_act) = .True.
! 
!      Case(26) ! minimal value such that a cross section is written out
!       Call SetWriteMinSig(wert)

     Case(31) ! setting a fixed GUT scale if wert > 0
      If (wert.Gt.0._dp) Call SetGUTScale(wert)
      Call SetGUTScale(wert)

     Case(32) ! requires strict unification
      If (Int(wert).Ne.0) check = SetStrictUnification(.True.)

     Case(33) ! setting a fixed renormalization scale if wert > 0
      If (wert.Gt.0._dp) Call SetRGEScale(wert**2)

     Case(34) ! precision of mass calculation
      delta_mass = wert

     Case(35) ! maximal number of iterations
      n_run = Int(wert)

     Case(36) ! minimal number of iterations
      MinimalNumberIterations = Int(wert)

!      Case(36) ! write out debug information
!       If (wert.Eq.0) Then
!        WriteOut = .False.
!       Else
!        WriteOut = .True.
!       End If

     Case(37) ! if =1 -> CKM thourgh V_u, if =2 CKM through V_d 
      If ((wert.Eq.1._dp).Or.(wert.Eq.2._dp)) i1 =  SetYukawaScheme(Int(wert))

     Case(38) ! set looplevel of RGEs
      If (wert.Ne.2._dp) Then
       TwoLoopRGE=.False.
      Else
       TwoLoopRGE=.True.
      End If

     Case(39) ! write out debug information
      If (wert.Eq.0) Then
       WriteSLHA1 = .False.
      Else
       WriteSLHA1 = .True.
      End If


     Case(40) ! alpha(0)
      check_alpha(2) = .True.
      Alpha = 1._dp / wert

     Case(41) ! Z-boson width
      gamZ = wert

     Case(42) ! W-boson width
      gamW = wert

     Case(50) ! W-boson width
      If (wert.Ne.1._dp) Then
       RotateNegativeFermionMasses=.False.
      Else
       RotateNegativeFermionMasses=.True.
      End If

     Case(51)
      If (wert.Ne.1._dp) Then
       SwitchToSCKM=.False.
      Else
       SwitchToSCKM=.True.
      End If

     Case(52)
      If (wert.Ne.1._dp) Then
       IgnoreNegativeMasses=.False.
      Else
       IgnoreNegativeMasses=.True.
      End If

     Case(53)
      If (wert.Ne.1._dp) Then
       IgnoreNegativeMassesMZ=.False.
      Else
       IgnoreNegativeMassesMZ=.True.
      End If

     Case(54)
      If (wert.Ne.1._dp) Then
       WriteOutputForNonConvergence=.False.
      Else
       WriteOutputForNonConvergence=.True.
      End If

     Case(55)
      If (wert.Ne.0._dp) Then
       CalculateOneLoopMasses=.True.
      Else
       CalculateOneLoopMasses=.False.
      End If

!      Case(56)
!       If (wert.Ne.0._dp) Then
!        CalculateTwoLoopHiggsMasses=.True.
!       Else
!        CalculateTwoLoopHiggsMasses=.False.
!       End If

     Case(57)
      If (wert.Ne.0._dp) Then
       CalculateLowEnergy=.True.
      Else
       CalculateLowEnergy=.False.
      End If

     Case(58)
      If (wert.Ne.0._dp) Then
        IncludeDeltaVB=.True.
        If (wert.Ne.2._dp) Then
         IncludeBSMdeltaVB=.True.
        Else
         IncludeBSMdeltaVB=.False.
        End If
      Else
       IncludeDeltaVB=.False.
      End If

     Case(60)
      If (wert.Ne.0._dp) Then
       KineticMixing=.True.
      Else
       KineticMixing=.False.
      End If

     Case(61)
      If (wert.Ne.0._dp) Then
       SMrunningLowScaleInput=.True.
      Else
       SMrunningLowScaleInput=.False.
      End If

     Case(62)
      If (wert.Ne.0._dp) Then
       RunningSUSYparametersLowEnergy=.True.
      Else
       RunningSUSYparametersLowEnergy=.False.
      End If

     Case(63)
      If (wert.Ne.0._dp) Then
       RunningSMparametersLowEnergy=.True.
      Else
       RunningSMparametersLowEnergy=.False.
      End If

     Case(64)
      If (wert.Ne.0._dp) Then
       WriteParametersAtQ=.True.
      Else
       WriteParametersAtQ=.False.
      End If

     Case(70)
      If (wert.Ne.0._dp) Then
       SUSYrunningFromMZ=.True.
      Else
       SUSYrunningFromMZ=.False.
      End If

     Case(65)
      If (wert.gt.0) SolutionTadpoleNr = wert 


     Case(75) ! Writes the parameter file for WHIZARD
      If (wert.Eq.1) Write_WHIZARD = .True.     

     Case(76) ! Writes input files for HiggsBounfs
      If (wert.Eq.1) Write_HiggsBounds = .True.      


     Case(80) ! exit for sure with non-zero value if a problem occurs
      If (wert.Eq.1) Non_Zero_Exit = .True.      

     Case(86) ! width to be counted as inivisble in HiggsBounds output
      WidthToBeInvisible = wert   

     Case(88) ! maximal mass allowed in loops
      MaxMassLoop = wert**2
   
     Case(89) ! maximal mass counted as numerical zero
      MaxMassNumericalZero = wert

     Case(95) ! Force mass matrices at 1-loop to be real
      If (wert.Eq.1) ForceRealMatrices  = .True.


!      Case(90) ! add R-parity at low energies
!       If (wert.Eq.1) Add_Rparity = .True.      
! 
!      Case(91) ! fit RP parameters such, that neutrino data are o.k.
!       If (wert.Eq.1) l_fit_RP_parameters = .True.      
! 
!      Case(92) ! for Pythia input
!       If (wert.Eq.1) l_RP_Pythia = .True.      
! 
!      Case(97) ! for Pythia input
!       If (wert.Eq.1) PrintPartialContributions = .True.     


     Case(510)
      If (wert.Ne.1._dp) Then
       WriteTreeLevelTadpoleSolutions=.False.
      Else
       WriteTreeLevelTadpoleSolutions=.True.
      End If

     Case(515)
      If (wert.Ne.0._dp) Then
       WriteGUTvalues=.True.
      Else
       WriteGUTvalues=.False.
      End If

     Case(520)
      If (wert.Ne.1._dp) Then
       WriteEffHiggsCouplingRatios=.False.
      Else
       WriteEffHiggsCouplingRatios=.True.
      End If


     Case(525)
      If (wert.Ne.1._dp) Then
       WriteHiggsDiphotonLoopContributions=.False.
      Else
       WriteHiggsDiphotonLoopContributions=.True.
      End If

     Case(530)
      If (wert.Ne.1._dp) Then
       WriteTreeLevelTadpoleParameters=.False.
      Else
       WriteTreeLevelTadpoleParameters=.True.
      End If

     Case(550)
      If (wert.Ne.1._dp) Then
       CalcFT=.False.
      Else
       CalcFT=.True.
      End If

     Case(551)
      If (wert.Ne.1._dp) Then
       OneLoopFT=.False.
      Else
       OneLoopFT=.True.
      End If



     Case(990)
      If (wert.Ne.1._dp) Then
       MakeQTEST=.False.
      Else
       MakeQTEST=.True.
      End If
      
!!-----------------------------------------------------------------
!! To be added in order to resume jobs from where it stopped
!! whenn postprocessing scan data
!     Case(995) !! Added by donald on Jan. 05 2018
!      If (int(wert) .Ne. 1) Then
!       ResumeJob=.False.
!      Else
!       ResumeJob=.True.
!      End If
!!-----------------------------------------------------------------

     Case(999)
      If (wert.Ne.1._dp) Then
       PrintDebugInformation=.False.
      Else
       PrintDebugInformation=.True.
      End If
 


    Case Default
      If (output_screen) Write(*,*) &
           & "Problem while reading SPhenoInput, ignoring unknown entry" &
           & ,i_par,wert
      Write(Errcan,*) &
          & "Problem while reading  SPhenoInput, ignoring unknown entry" &
               & ,i_par,wert
     End Select ! i_par

    End Do  ! i_par 

   200 Return

  End Subroutine Read_SPhenoInput
End Subroutine LesHouches_Input 
 
 
 
Subroutine LesHouches_Out(io_L,io,kont,M_GUT,ae,amu,atau,EDMe,EDMmu,EDMtau,           & 
& dRho,BrBsGamma,ratioBsGamma,BrDmunu,ratioDmunu,BrDsmunu,ratioDsmunu,BrDstaunu,         & 
& ratioDstaunu,BrBmunu,ratioBmunu,BrBtaunu,ratioBtaunu,BrKmunu,ratioKmunu,               & 
& RK,RKSM,muEgamma,tauEgamma,tauMuGamma,CRmuEAl,CRmuETi,CRmuESr,CRmuESb,CRmuEAu,         & 
& CRmuEPb,BRmuTo3e,BRtauTo3e,BRtauTo3mu,BRtauToemumu,BRtauTomuee,BRtauToemumu2,          & 
& BRtauTomuee2,BrZtoMuE,BrZtoTauE,BrZtoTauMu,BrhtoMuE,BrhtoTauE,BrhtoTauMu,              & 
& DeltaMBs,ratioDeltaMBs,DeltaMBq,ratioDeltaMBq,BrTautoEPi,BrTautoEEta,BrTautoEEtap,     & 
& BrTautoMuPi,BrTautoMuEta,BrTautoMuEtap,BrB0dEE,ratioB0dEE,BrB0sEE,ratioB0sEE,          & 
& BrB0dMuMu,ratioB0dMuMu,BrB0sMuMu,ratioB0sMuMu,BrB0dTauTau,ratioB0dTauTau,              & 
& BrB0sTauTau,ratioB0sTauTau,BrBtoSEE,ratioBtoSEE,BrBtoSMuMu,ratioBtoSMuMu,              & 
& BrBtoKmumu,ratioBtoKmumu,BrBtoSnunu,ratioBtoSnunu,BrBtoDnunu,ratioBtoDnunu,            & 
& BrKptoPipnunu,ratioKptoPipnunu,BrKltoPinunu,ratioKltoPinunu,DelMK,ratioDelMK,          & 
& epsK,ratioepsK,GenerationMixing,f_name)

Implicit None 
Integer, Intent(in) :: io_L, io, kont
Real(dp),Intent(in) :: ae,amu,atau,EDMe,EDMmu,EDMtau,dRho,BrBsGamma,ratioBsGamma,BrDmunu,ratioDmunu,         & 
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

Real(dp), Intent(in) :: M_GUT
Character(len=8)::Datum 
Character(len=10)::Zeit 
Logical,Intent(in)::GenerationMixing 
Logical,Save::l_open= .True. 
Integer,Parameter::n_max=500 
Integer :: i1, i2 
Logical :: WriteNextBlock 
Character(len=30),Dimension(n_max)::Fnames,Lnames 
Character(len=*),Intent(in),Optional::f_name
Real(dp) :: Q, MassLSP(2), facPP, facGG, facPZ 
Complex(dp) :: CKM_Q(3,3), ZD_ckm(6,6),ZU_ckm(6,6), &  
 & mq2_ckm(3,3),  mu2_ckm(3,3),  md2_ckm(3,3),  Tu_ckm(3,3), &  
 &Td_ckm(3,3),  Yd_ckm(3,3),  Yu_ckm(3,3) 
Integer :: CurrentPDG2(2), CurrentPDG3(3), PDGlsp(2) 
Integer::ierr,i_errors(1100),gt1,gt2,gt3,icount
Complex(dp) :: PDGSd(6),PDGSu(6),PDGhh(8),PDGAh(8),PDGHpm(8),PDGVP,PDGVZ,PDGVG,PDGVWm,               & 
& PDGgP,PDGgWm,PDGgWmC,PDGgZ,PDGgG,PDGFd(3),PDGFu(3),PDGGlu,PDGChi(10),PDGCha(5)

Character(len=30),Dimension(6):: NameParticleSd
Character(len=30),Dimension(6):: NameParticleSu
Character(len=30),Dimension(8):: NameParticlehh
Character(len=30),Dimension(8):: NameParticleAh
Character(len=30),Dimension(8):: NameParticleHpm
Character(len=30) :: NameParticleVP
Character(len=30) :: NameParticleVZ
Character(len=30) :: NameParticleVG
Character(len=30) :: NameParticleVWm
Character(len=30) :: NameParticlegP
Character(len=30) :: NameParticlegWm
Character(len=30) :: NameParticlegWmC
Character(len=30) :: NameParticlegZ
Character(len=30) :: NameParticlegG
Character(len=30),Dimension(3):: NameParticleFd
Character(len=30),Dimension(3):: NameParticleFu
Character(len=30) :: NameParticleGlu
Character(len=30),Dimension(10):: NameParticleChi
Character(len=30),Dimension(5):: NameParticleCha
Complex(dp) :: Zbottom(2,2), Ztop(2,2), Ztau(2,2) 

 
 
 ! ----------- Set names and PDGs -------- 
 
PDGSd(1)=1000001
NameParticleSd(1)="Sd_1"
PDGSd(2)=1000003
NameParticleSd(2)="Sd_2"
PDGSd(3)=1000005
NameParticleSd(3)="Sd_3"
PDGSd(4)=2000001
NameParticleSd(4)="Sd_4"
PDGSd(5)=2000003
NameParticleSd(5)="Sd_5"
PDGSd(6)=2000005
NameParticleSd(6)="Sd_6"
PDGSu(1)=1000002
NameParticleSu(1)="Su_1"
PDGSu(2)=1000004
NameParticleSu(2)="Su_2"
PDGSu(3)=1000006
NameParticleSu(3)="Su_3"
PDGSu(4)=2000002
NameParticleSu(4)="Su_4"
PDGSu(5)=2000004
NameParticleSu(5)="Su_5"
PDGSu(6)=2000006
NameParticleSu(6)="Su_6"
PDGhh(1)=25
NameParticlehh(1)="hh_1"
PDGhh(2)=35
NameParticlehh(2)="hh_2"
PDGhh(3)=1000012
NameParticlehh(3)="hh_3"
PDGhh(4)=1000014
NameParticlehh(4)="hh_4"
PDGhh(5)=1000016
NameParticlehh(5)="hh_5"
PDGhh(6)=2000012
NameParticlehh(6)="hh_6"
PDGhh(7)=2000014
NameParticlehh(7)="hh_7"
PDGhh(8)=2000016
NameParticlehh(8)="hh_8"
PDGAh(1)=0
NameParticleAh(1)="Ah_1"
PDGAh(2)=36
NameParticleAh(2)="Ah_2"
PDGAh(3)=1000017
NameParticleAh(3)="Ah_3"
PDGAh(4)=1000018
NameParticleAh(4)="Ah_4"
PDGAh(5)=1000019
NameParticleAh(5)="Ah_5"
PDGAh(6)=2000018
NameParticleAh(6)="Ah_6"
PDGAh(7)=2000019
NameParticleAh(7)="Ah_7"
PDGAh(8)=2000020
NameParticleAh(8)="Ah_8"
PDGHpm(1)=0
NameParticleHpm(1)="Hpm_1"
PDGHpm(2)=37
NameParticleHpm(2)="Hpm_2"
PDGHpm(3)=1000011
NameParticleHpm(3)="Hpm_3"
PDGHpm(4)=2000011
NameParticleHpm(4)="Hpm_4"
PDGHpm(5)=1000013
NameParticleHpm(5)="Hpm_5"
PDGHpm(6)=2000013
NameParticleHpm(6)="Hpm_6"
PDGHpm(7)=1000015
NameParticleHpm(7)="Hpm_7"
PDGHpm(8)=2000015
NameParticleHpm(8)="Hpm_8"
PDGVP=22
NameParticleVP="VP"
PDGVZ=23
NameParticleVZ="VZ"
PDGVG=21
NameParticleVG="VG"
PDGVWm=-24
NameParticleVWm="VWm"
PDGgP=0
NameParticlegP="gP"
PDGgWm=0
NameParticlegWm="gWm"
PDGgWmC=0
NameParticlegWmC="gWmC"
PDGgZ=0
NameParticlegZ="gZ"
PDGgG=0
NameParticlegG="gG"
PDGFd(1)=1
NameParticleFd(1)="Fd_1"
PDGFd(2)=3
NameParticleFd(2)="Fd_2"
PDGFd(3)=5
NameParticleFd(3)="Fd_3"
PDGFu(1)=2
NameParticleFu(1)="Fu_1"
PDGFu(2)=4
NameParticleFu(2)="Fu_2"
PDGFu(3)=6
NameParticleFu(3)="Fu_3"
PDGGlu=1000021
NameParticleGlu="Glu"
PDGChi(1)=12
NameParticleChi(1)="Chi_1"
PDGChi(2)=14
NameParticleChi(2)="Chi_2"
PDGChi(3)=16
NameParticleChi(3)="Chi_3"
PDGChi(4)=1000022
NameParticleChi(4)="Chi_4"
PDGChi(5)=1000023
NameParticleChi(5)="Chi_5"
PDGChi(6)=1000025
NameParticleChi(6)="Chi_6"
PDGChi(7)=1000039
NameParticleChi(7)="Chi_7"
PDGChi(8)=1000045
NameParticleChi(8)="Chi_8"
PDGChi(9)=1000055
NameParticleChi(9)="Chi_9"
PDGChi(10)=1000065
NameParticleChi(10)="Chi_10"
PDGCha(1)=11
NameParticleCha(1)="Cha_1"
PDGCha(2)=13
NameParticleCha(2)="Cha_2"
PDGCha(3)=15
NameParticleCha(3)="Cha_3"
PDGCha(4)=-1000024
NameParticleCha(4)="Cha_4"
PDGCha(5)=-1000037
NameParticleCha(5)="Cha_5"

 
 
 ! ----------- Use SLHA 1 conventions if demanded -------- 
 
If(WriteSLHA1) Write(*,*) "SLHA 1 output for given model not possible" 
Q=Sqrt(GetRenormalizationScale())
Call Date_and_time(datum,zeit)

!If (l_open) Then
!If (Present(f_name)) Then
!Open(io_L,file=Trim(f_name),status="unknown")
!Else
Open(io_L,file=outputFileName,status="unknown")
!End If
l_open= .False.
!End If

If (.Not.RotateNegativeFermionMasses) Then 
Do i1=1,10
  If (MaxVal(Abs(AImag(UV(i1,:)))).gt.MaxVal(Abs(Real(UV(i1,:),dp)))) Then 
MChi(i1) = - MChi(i1) 
UV(i1,:) = -(0._dp,1._dp)*UV(i1,:) 
  End If 
End Do 
End If 
 
Write(io_L,100) "# SUSY Les Houches Accord 2 - munuSSM Spectrum + Decays + Flavor Observables"
Write(io_L,100) "# SPheno module generated by SARAH" 
Write(io_L,100) "# ----------------------------------------------------------------------" 
Write(io_L,100) "# SPheno "//version 
Write(io_L,100) "#   W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
Write(io_L,100) "#   W. Porod, F.Staub, Comput.Phys.Commun.183 (2012) 2458-2469, arXiv:1104.1573"
Write(io_L,100) "# SARAH: "//versionSARAH 
Write(io_L,100) "#   F. Staub; arXiv:0806.0538 (online manual)"
Write(io_L,100) "#   F. Staub; Comput. Phys. Commun. 181 (2010) 1077-1086; arXiv:0909.2863"
Write(io_L,100) "#   F. Staub; Comput. Phys. Commun. 182 (2011)  808-833; arXiv:1002.0840"
Write(io_L,100) "#   F. Staub; Comput. Phys. Commun. 184 (2013)  1792-1809; arXiv:1207.0906"
Write(io_L,100) "#   F. Staub; Comput. Phys. Commun. 185 (2014)  1773-1790; arXiv:1309.7223 "
Write(io_L,100) "# Including the calculation of flavor observables based on the FlavorKit "
Write(io_L,100) "#   W. Porod, F. Staub, A. Vicente; Eur.Phys.J. C74 (2014) 8, 2992; arXiv:1405.1434 "
Write(io_L,100) "# Two-loop masss corrections to Higgs fields based on "
Write(io_L,100) "#   M. D. Goodsell, K. Nickel, F. Staub; arXiv:1411.0675 "
Write(io_L,100) "#   M. D. Goodsell, K. Nickel, F. Staub; arXiv:1503.03098"
Write(io_L,100) "#  "
Write(io_L,100) "# in case of problems send email to florian.staub@cern.ch and goodsell@lpthe.jussieu.fr "
Write(io_L,100) "# ----------------------------------------------------------------------" 
Write(io_L,100) "# Created: "//Datum(7:8)//"."//Datum(5:6)//"."//Datum(1:4)&
&//",  "//Zeit(1:2)//":"//Zeit(3:4)
Write(io_L,100) "Block SPINFO         # Program information"
Write(io_L,100) "     1   SPhenoSARAH      # spectrum calculator"
Write(io_L,100) "     2   "//version//"    # version number of SPheno"
Write(io_L,100) "     9   "//versionSARAH//"    # version number of SARAH"
Call GetError(i_errors)
If ((i_errors(1)+i_errors(3)+i_errors(5)+i_errors(7)+i_errors(8)& 
&+i_errors(10)+i_errors(12)+Sum(i_errors(14:19))).Gt.0)&
& Write(io_L,100)&
& "     3               # potential numerical problem, check file Messages.out"
Write(io_L,100) "Block MODSEL  # Input parameters"
If (HighScaleModel.Eq."LOW") Then 
 Write(io_L,110)  1,0, " SUSY Scale input"
Else  
 Write(io_L,110)  1, 1, " GUT scale input"
End If  
 Write(io_L,110) 2, BoundaryCondition, " Boundary conditions "
If (i_cpv.Gt.0) Write(io_L,110) 5,i_cpv," switching on CP violation"
If (GenerationMixing) Write(io_L,110) &
&     6,1, " switching on flavour violation" 
Write(io_L,100) "Block MINPAR  # Input parameters"
Write(io_L,101) 1, Real(m0,dp) ,"# m0"
Write(io_L,101) 2, Real(m12,dp) ,"# m12"
Write(io_L,101) 3, Real(TanBeta,dp) ,"# TanBeta"
Write(io_L,101) 5, Real(Azero,dp) ,"# Azero"
WriteNextBlock = .False. 
If (Abs(Aimag(m12)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(Azero)).gt.0._dp) WriteNextBlock = .True. 
If(WriteNextBlock) Then 
Write(io_L,100) "Block IMMINPAR  # Input parameters"
If (Abs(Aimag(m12)).gt.0._dp) Then 
Write(io_L,101) 2, Aimag(m12) ,"# m12"
End if 
If (Abs(Aimag(Azero)).gt.0._dp) Then 
Write(io_L,101) 5, Aimag(Azero) ,"# Azero"
End if 
End if 
WriteNextBlock = .false. 
If (Abs(Real(LambdaInput,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(KappaInput,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(ALambdaInput,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(AKappaInput,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(vR1Input,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(vR2Input,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(vR3Input,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(vL1Input,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(vL2Input,dp)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Real(vL3Input,dp)).gt.0._dp) WriteNextBlock = .True. 
If(WriteNextBlock) Then 
Write(io_L,100) "Block EXTPAR  # Input parameters"
If (Abs(Real(LambdaInput,dp)).gt.0._dp) Then 
Write(io_L,101) 61, Real(LambdaInput,dp) ,"# LambdaInput"
End if 
If (Abs(Real(KappaInput,dp)).gt.0._dp) Then 
Write(io_L,101) 62, Real(KappaInput,dp) ,"# KappaInput"
End if 
If (Abs(Real(ALambdaInput,dp)).gt.0._dp) Then 
Write(io_L,101) 63, Real(ALambdaInput,dp) ,"# ALambdaInput"
End if 
If (Abs(Real(AKappaInput,dp)).gt.0._dp) Then 
Write(io_L,101) 64, Real(AKappaInput,dp) ,"# AKappaInput"
End if 
If (Abs(Real(vR1Input,dp)).gt.0._dp) Then 
Write(io_L,101) 65, Real(vR1Input,dp) ,"# vR1Input"
End if 
If (Abs(Real(vR2Input,dp)).gt.0._dp) Then 
Write(io_L,101) 66, Real(vR2Input,dp) ,"# vR2Input"
End if 
If (Abs(Real(vR3Input,dp)).gt.0._dp) Then 
Write(io_L,101) 67, Real(vR3Input,dp) ,"# vR3Input"
End if 
If (Abs(Real(vL1Input,dp)).gt.0._dp) Then 
Write(io_L,101) 200, Real(vL1Input,dp) ,"# vL1Input"
End if 
If (Abs(Real(vL2Input,dp)).gt.0._dp) Then 
Write(io_L,101) 201, Real(vL2Input,dp) ,"# vL2Input"
End if 
If (Abs(Real(vL3Input,dp)).gt.0._dp) Then 
Write(io_L,101) 202, Real(vL3Input,dp) ,"# vL3Input"
End if 
End if 
WriteNextBlock = .false. 
If (Abs(Aimag(LambdaInput)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(KappaInput)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(ALambdaInput)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(AKappaInput)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(vR1Input)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(vR2Input)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(vR3Input)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(vL1Input)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(vL2Input)).gt.0._dp) WriteNextBlock = .True. 
If (Abs(Aimag(vL3Input)).gt.0._dp) WriteNextBlock = .True. 
If(WriteNextBlock) Then 
Write(io_L,100) "Block IMEXTPAR  # Input parameters"
If (Abs(Aimag(LambdaInput)).gt.0._dp) Then 
Write(io_L,101) 61, Aimag(LambdaInput) ,"# LambdaInput"
End if 
If (Abs(Aimag(KappaInput)).gt.0._dp) Then 
Write(io_L,101) 62, Aimag(KappaInput) ,"# KappaInput"
End if 
If (Abs(Aimag(ALambdaInput)).gt.0._dp) Then 
Write(io_L,101) 63, Aimag(ALambdaInput) ,"# ALambdaInput"
End if 
If (Abs(Aimag(AKappaInput)).gt.0._dp) Then 
Write(io_L,101) 64, Aimag(AKappaInput) ,"# AKappaInput"
End if 
If (Abs(Aimag(vR1Input)).gt.0._dp) Then 
Write(io_L,101) 65, Aimag(vR1Input) ,"# vR1Input"
End if 
If (Abs(Aimag(vR2Input)).gt.0._dp) Then 
Write(io_L,101) 66, Aimag(vR2Input) ,"# vR2Input"
End if 
If (Abs(Aimag(vR3Input)).gt.0._dp) Then 
Write(io_L,101) 67, Aimag(vR3Input) ,"# vR3Input"
End if 
If (Abs(Aimag(vL1Input)).gt.0._dp) Then 
Write(io_L,101) 200, Aimag(vL1Input) ,"# vL1Input"
End if 
If (Abs(Aimag(vL2Input)).gt.0._dp) Then 
Write(io_L,101) 201, Aimag(vL2Input) ,"# vL2Input"
End if 
If (Abs(Aimag(vL3Input)).gt.0._dp) Then 
Write(io_L,101) 202, Aimag(vL3Input) ,"# vL3Input"
End if 
End if 
Write(io_L,106) "Block gaugeGUT Q=",m_GUT,"# (GUT scale)" 
Write(io_L,104) 1,g1GUT, "# g1(Q)^DRbar" 
Write(io_L,104) 2,g2GUT, "# g2(Q)^DRbar" 
Write(io_L,104) 3,g3GUT, "# g3(Q)^DRbar" 
Write(io_L,100) "Block SMINPUTS  # SM parameters"
Write(io_L,102) 1,1._dp/alpha_MSbar(mZ,mW),"# alpha_em^-1(MZ)^MSbar"
Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
Write(io_L,102) 4,mZ,"# m_Z(pole)"
Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"
 
If (SwitchToSCKM) Then
Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
Write(io_L,102) 1,lam_wolf,"# lambda"
Write(io_L,102) 2,A_wolf,"# A"
Write(io_L,102) 3,rho_wolf,"# rho bar"
Write(io_L,102) 4,eta_wolf,"# eta bar"
ZU_ckm = ZU(1:6,1:6) 
ZD_ckm = ZD(1:6,1:6) 
Call Switch_to_superCKM(Yd(1:3,1:3),Yu(1:3,1:3),Td(1:3,1:3),Tu(1:3,1:3),md2(1:3,1:3),mq2(1:3,1:3),mu2(1:3,1:3) &
&,Td_ckm,Tu_ckm,md2_ckm,mq2_ckm,mu2_ckm, .False.&
&,ZD_ckm,ZU_ckm,ZD,ZU,CKM_Q,Yd_ckm,Yu_ckm)
Yd(1:3,1:3)=Yd_ckm 
Yu(1:3,1:3)=Yu_ckm 
Td(1:3,1:3)=Td_ckm 
Tu(1:3,1:3)=Tu_ckm 
md2(1:3,1:3)=md2_ckm 
mu2(1:3,1:3)=mu2_ckm 
mq2(1:3,1:3)=mq2_ckm 
ZU(1:6,1:6)=ZU_ckm 
ZD(1:6,1:6)=ZD_ckm 
Write(io_L,106) "Block VCKM Q=",Q,"# Re(CKM) at the SUSY Scale" 
Do i1=1,3
Do i2=1,3
Write(io_L,107) i1,i2,Real(CKM_Q(i1,i2),dp),"# Re(V_"//Bu(i1)//Bu(i2)//")"
End Do
End Do
If (Maxval(Abs(Aimag(CKM_Q))).Gt.0._dp) Then
Write(io_L,106) "Block IMVCKM Q=",Q,"# Im(CKM) at the SUSY Scale" 
Do i1=1,3
Do i2=1,3
Write(io_L,107) i1,i2,Aimag(CKM_Q(i1,i2)),"# Im(V_"//Bu(i1)//Bu(i2)//")"
End Do
End Do
End If
End If 
 
 
WriteNextBlock = .false. 
Write(io_L,106) "Block GAUGE Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Real(g1,dp), "# g1" 
Write(io_L,104) 2,Real(g2,dp), "# g2" 
Write(io_L,104) 3,Real(g3,dp), "# g3" 
If(WriteNextBlock) Then 
Write(io_L,106) "Block IMGAUGE Q=",Q,"# (SUSY Scale)" 
End if 
WriteNextBlock = .false. 
Write(io_L,106) "Block MSOFT Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 21,Real(mHd2,dp), "# mHd2" 
Write(io_L,104) 22,Real(mHu2,dp), "# mHu2" 
Write(io_L,104) 1,Real(M1,dp), "# M1" 
If (Abs(Aimag(M1)).gt.0._dp) WriteNextBlock = .True. 
Write(io_L,104) 2,Real(M2,dp), "# M2" 
If (Abs(Aimag(M2)).gt.0._dp) WriteNextBlock = .True. 
Write(io_L,104) 3,Real(M3,dp), "# M3" 
If (Abs(Aimag(M3)).gt.0._dp) WriteNextBlock = .True. 
If(WriteNextBlock) Then 
Write(io_L,106) "Block IMMSOFT Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Aimag(M1), "# M1" 
Write(io_L,104) 2,Aimag(M2), "# M2" 
Write(io_L,104) 3,Aimag(M3), "# M3" 
End if 
WriteNextBlock = .false. 
Write(io_L,106) "Block HMIX Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 102,Real(vd,dp), "# vd" 
Write(io_L,104) 103,Real(vu,dp), "# vu" 
If(WriteNextBlock) Then 
Write(io_L,106) "Block IMHMIX Q=",Q,"# (SUSY Scale)" 
End if 
WriteNextBlock = .false. 
Write(io_L,106) "Block PHASES Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Real(pG,dp), "# pG" 
If (Abs(Aimag(pG)).gt.0._dp) WriteNextBlock = .True. 
If(WriteNextBlock) Then 
Write(io_L,106) "Block IMPHASES Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Aimag(pG), "# pG" 
End if 
If (WriteTreeLevelTadpoleParameters) Then 
If (HighScaleModel.Eq."LOW") Then 
WriteNextBlock = .false. 
Write(io_L,106) "Block TREEMSOFT Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 21,Real(mHd2Tree,dp), "# mHd2" 
Write(io_L,104) 22,Real(mHu2Tree,dp), "# mHu2" 
If(WriteNextBlock) Then 
Write(io_L,106) "Block TREEIMMSOFT Q=",Q,"# (SUSY Scale)" 
End if 
WriteNextBlock = .false. 
Write(io_L,106) "Block LOOPMSOFT Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 21,Real(mHd21L,dp), "# mHd2" 
Write(io_L,104) 22,Real(mHu21L,dp), "# mHu2" 
If(WriteNextBlock) Then 
Write(io_L,106) "Block LOOPIMMSOFT Q=",Q,"# (SUSY Scale)" 
End if 
Else 
WriteNextBlock = .false. 
Write(io_L,106) "Block TREEMSOFT Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 21,Real(mHd2Tree,dp), "# mHd2" 
Write(io_L,104) 22,Real(mHu2Tree,dp), "# mHu2" 
If(WriteNextBlock) Then 
Write(io_L,106) "Block TREEIMMSOFT Q=",Q,"# (SUSY Scale)" 
End if 
WriteNextBlock = .false. 
Write(io_L,106) "Block LOOPMSOFT Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 21,Real(mHd21L,dp), "# mHd2" 
Write(io_L,104) 22,Real(mHu21L,dp), "# mHu2" 
If(WriteNextBlock) Then 
Write(io_L,106) "Block LOOPIMMSOFT Q=",Q,"# (SUSY Scale)" 
End if 
End if 
End if 
Write(io_L,106) "Block Yd Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Yd(1,1),dp), "# Real(Yd(1,1),dp)" 
Write(io_L,107)1,2,Real(Yd(1,2),dp), "# Real(Yd(1,2),dp)" 
Write(io_L,107)1,3,Real(Yd(1,3),dp), "# Real(Yd(1,3),dp)" 
Write(io_L,107)2,1,Real(Yd(2,1),dp), "# Real(Yd(2,1),dp)" 
Write(io_L,107)2,2,Real(Yd(2,2),dp), "# Real(Yd(2,2),dp)" 
Write(io_L,107)2,3,Real(Yd(2,3),dp), "# Real(Yd(2,3),dp)" 
Write(io_L,107)3,1,Real(Yd(3,1),dp), "# Real(Yd(3,1),dp)" 
Write(io_L,107)3,2,Real(Yd(3,2),dp), "# Real(Yd(3,2),dp)" 
Write(io_L,107)3,3,Real(Yd(3,3),dp), "# Real(Yd(3,3),dp)" 
If (MaxVal(Abs(AImag(Yd))).gt.0._dp) Then 
Write(io_L,106) "Block IMYd Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Yd(1,1)), "# Aimag(Yd(1,1))" 
Write(io_L,107)1,2,Aimag(Yd(1,2)), "# Aimag(Yd(1,2))" 
Write(io_L,107)1,3,Aimag(Yd(1,3)), "# Aimag(Yd(1,3))" 
Write(io_L,107)2,1,Aimag(Yd(2,1)), "# Aimag(Yd(2,1))" 
Write(io_L,107)2,2,Aimag(Yd(2,2)), "# Aimag(Yd(2,2))" 
Write(io_L,107)2,3,Aimag(Yd(2,3)), "# Aimag(Yd(2,3))" 
Write(io_L,107)3,1,Aimag(Yd(3,1)), "# Aimag(Yd(3,1))" 
Write(io_L,107)3,2,Aimag(Yd(3,2)), "# Aimag(Yd(3,2))" 
Write(io_L,107)3,3,Aimag(Yd(3,3)), "# Aimag(Yd(3,3))" 
End If 

Write(io_L,106) "Block Ye Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Ye(1,1),dp), "# Real(Ye(1,1),dp)" 
Write(io_L,107)1,2,Real(Ye(1,2),dp), "# Real(Ye(1,2),dp)" 
Write(io_L,107)1,3,Real(Ye(1,3),dp), "# Real(Ye(1,3),dp)" 
Write(io_L,107)2,1,Real(Ye(2,1),dp), "# Real(Ye(2,1),dp)" 
Write(io_L,107)2,2,Real(Ye(2,2),dp), "# Real(Ye(2,2),dp)" 
Write(io_L,107)2,3,Real(Ye(2,3),dp), "# Real(Ye(2,3),dp)" 
Write(io_L,107)3,1,Real(Ye(3,1),dp), "# Real(Ye(3,1),dp)" 
Write(io_L,107)3,2,Real(Ye(3,2),dp), "# Real(Ye(3,2),dp)" 
Write(io_L,107)3,3,Real(Ye(3,3),dp), "# Real(Ye(3,3),dp)" 
If (MaxVal(Abs(AImag(Ye))).gt.0._dp) Then 
Write(io_L,106) "Block IMYe Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Ye(1,1)), "# Aimag(Ye(1,1))" 
Write(io_L,107)1,2,Aimag(Ye(1,2)), "# Aimag(Ye(1,2))" 
Write(io_L,107)1,3,Aimag(Ye(1,3)), "# Aimag(Ye(1,3))" 
Write(io_L,107)2,1,Aimag(Ye(2,1)), "# Aimag(Ye(2,1))" 
Write(io_L,107)2,2,Aimag(Ye(2,2)), "# Aimag(Ye(2,2))" 
Write(io_L,107)2,3,Aimag(Ye(2,3)), "# Aimag(Ye(2,3))" 
Write(io_L,107)3,1,Aimag(Ye(3,1)), "# Aimag(Ye(3,1))" 
Write(io_L,107)3,2,Aimag(Ye(3,2)), "# Aimag(Ye(3,2))" 
Write(io_L,107)3,3,Aimag(Ye(3,3)), "# Aimag(Ye(3,3))" 
End If 

Write(io_L,106) "Block {NMSSMRUN, 1} Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Real(lam(1),dp), "# Real(lam(1) ,dp)" 
Write(io_L,104) 2,Real(lam(2),dp), "# Real(lam(2) ,dp)" 
Write(io_L,104) 3,Real(lam(3),dp), "# Real(lam(3) ,dp)" 
If (MaxVal(Abs(AImag(lam))).gt.0._dp) Then 
Write(io_L,106) "Block IM{NMSSMRUN, 1} Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Aimag(lam(1)), "# Aimag(lam(1) )" 
Write(io_L,104) 2,Aimag(lam(2)), "# Aimag(lam(2) )" 
Write(io_L,104) 3,Aimag(lam(3)), "# Aimag(lam(3) )" 
End If 

Write(io_L,106) "Block Yv Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Yv(1,1),dp), "# Real(Yv(1,1),dp)" 
Write(io_L,107)1,2,Real(Yv(1,2),dp), "# Real(Yv(1,2),dp)" 
Write(io_L,107)1,3,Real(Yv(1,3),dp), "# Real(Yv(1,3),dp)" 
Write(io_L,107)2,1,Real(Yv(2,1),dp), "# Real(Yv(2,1),dp)" 
Write(io_L,107)2,2,Real(Yv(2,2),dp), "# Real(Yv(2,2),dp)" 
Write(io_L,107)2,3,Real(Yv(2,3),dp), "# Real(Yv(2,3),dp)" 
Write(io_L,107)3,1,Real(Yv(3,1),dp), "# Real(Yv(3,1),dp)" 
Write(io_L,107)3,2,Real(Yv(3,2),dp), "# Real(Yv(3,2),dp)" 
Write(io_L,107)3,3,Real(Yv(3,3),dp), "# Real(Yv(3,3),dp)" 
If (MaxVal(Abs(AImag(Yv))).gt.0._dp) Then 
Write(io_L,106) "Block IMYv Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Yv(1,1)), "# Aimag(Yv(1,1))" 
Write(io_L,107)1,2,Aimag(Yv(1,2)), "# Aimag(Yv(1,2))" 
Write(io_L,107)1,3,Aimag(Yv(1,3)), "# Aimag(Yv(1,3))" 
Write(io_L,107)2,1,Aimag(Yv(2,1)), "# Aimag(Yv(2,1))" 
Write(io_L,107)2,2,Aimag(Yv(2,2)), "# Aimag(Yv(2,2))" 
Write(io_L,107)2,3,Aimag(Yv(2,3)), "# Aimag(Yv(2,3))" 
Write(io_L,107)3,1,Aimag(Yv(3,1)), "# Aimag(Yv(3,1))" 
Write(io_L,107)3,2,Aimag(Yv(3,2)), "# Aimag(Yv(3,2))" 
Write(io_L,107)3,3,Aimag(Yv(3,3)), "# Aimag(Yv(3,3))" 
End If 

Write(io_L,106) "Block Yu Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Yu(1,1),dp), "# Real(Yu(1,1),dp)" 
Write(io_L,107)1,2,Real(Yu(1,2),dp), "# Real(Yu(1,2),dp)" 
Write(io_L,107)1,3,Real(Yu(1,3),dp), "# Real(Yu(1,3),dp)" 
Write(io_L,107)2,1,Real(Yu(2,1),dp), "# Real(Yu(2,1),dp)" 
Write(io_L,107)2,2,Real(Yu(2,2),dp), "# Real(Yu(2,2),dp)" 
Write(io_L,107)2,3,Real(Yu(2,3),dp), "# Real(Yu(2,3),dp)" 
Write(io_L,107)3,1,Real(Yu(3,1),dp), "# Real(Yu(3,1),dp)" 
Write(io_L,107)3,2,Real(Yu(3,2),dp), "# Real(Yu(3,2),dp)" 
Write(io_L,107)3,3,Real(Yu(3,3),dp), "# Real(Yu(3,3),dp)" 
If (MaxVal(Abs(AImag(Yu))).gt.0._dp) Then 
Write(io_L,106) "Block IMYu Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Yu(1,1)), "# Aimag(Yu(1,1))" 
Write(io_L,107)1,2,Aimag(Yu(1,2)), "# Aimag(Yu(1,2))" 
Write(io_L,107)1,3,Aimag(Yu(1,3)), "# Aimag(Yu(1,3))" 
Write(io_L,107)2,1,Aimag(Yu(2,1)), "# Aimag(Yu(2,1))" 
Write(io_L,107)2,2,Aimag(Yu(2,2)), "# Aimag(Yu(2,2))" 
Write(io_L,107)2,3,Aimag(Yu(2,3)), "# Aimag(Yu(2,3))" 
Write(io_L,107)3,1,Aimag(Yu(3,1)), "# Aimag(Yu(3,1))" 
Write(io_L,107)3,2,Aimag(Yu(3,2)), "# Aimag(Yu(3,2))" 
Write(io_L,107)3,3,Aimag(Yu(3,3)), "# Aimag(Yu(3,3))" 
End If 

Write(io_L,106) "Block {NMSSMRUN, 2} Q=",Q,"# (SUSY Scale)" 
Write(io_L,127) 1,1,1,Real(kap(1,1,1),dp), "# Real(kap(1,1,1),dp)" 
Write(io_L,127) 1,1,2,Real(kap(1,1,2),dp), "# Real(kap(1,1,2),dp)" 
Write(io_L,127) 1,1,3,Real(kap(1,1,3),dp), "# Real(kap(1,1,3),dp)" 
Write(io_L,127) 1,2,1,Real(kap(1,2,1),dp), "# Real(kap(1,2,1),dp)" 
Write(io_L,127) 1,2,2,Real(kap(1,2,2),dp), "# Real(kap(1,2,2),dp)" 
Write(io_L,127) 1,2,3,Real(kap(1,2,3),dp), "# Real(kap(1,2,3),dp)" 
Write(io_L,127) 1,3,1,Real(kap(1,3,1),dp), "# Real(kap(1,3,1),dp)" 
Write(io_L,127) 1,3,2,Real(kap(1,3,2),dp), "# Real(kap(1,3,2),dp)" 
Write(io_L,127) 1,3,3,Real(kap(1,3,3),dp), "# Real(kap(1,3,3),dp)" 
Write(io_L,127) 2,1,1,Real(kap(2,1,1),dp), "# Real(kap(2,1,1),dp)" 
Write(io_L,127) 2,1,2,Real(kap(2,1,2),dp), "# Real(kap(2,1,2),dp)" 
Write(io_L,127) 2,1,3,Real(kap(2,1,3),dp), "# Real(kap(2,1,3),dp)" 
Write(io_L,127) 2,2,1,Real(kap(2,2,1),dp), "# Real(kap(2,2,1),dp)" 
Write(io_L,127) 2,2,2,Real(kap(2,2,2),dp), "# Real(kap(2,2,2),dp)" 
Write(io_L,127) 2,2,3,Real(kap(2,2,3),dp), "# Real(kap(2,2,3),dp)" 
Write(io_L,127) 2,3,1,Real(kap(2,3,1),dp), "# Real(kap(2,3,1),dp)" 
Write(io_L,127) 2,3,2,Real(kap(2,3,2),dp), "# Real(kap(2,3,2),dp)" 
Write(io_L,127) 2,3,3,Real(kap(2,3,3),dp), "# Real(kap(2,3,3),dp)" 
Write(io_L,127) 3,1,1,Real(kap(3,1,1),dp), "# Real(kap(3,1,1),dp)" 
Write(io_L,127) 3,1,2,Real(kap(3,1,2),dp), "# Real(kap(3,1,2),dp)" 
Write(io_L,127) 3,1,3,Real(kap(3,1,3),dp), "# Real(kap(3,1,3),dp)" 
Write(io_L,127) 3,2,1,Real(kap(3,2,1),dp), "# Real(kap(3,2,1),dp)" 
Write(io_L,127) 3,2,2,Real(kap(3,2,2),dp), "# Real(kap(3,2,2),dp)" 
Write(io_L,127) 3,2,3,Real(kap(3,2,3),dp), "# Real(kap(3,2,3),dp)" 
Write(io_L,127) 3,3,1,Real(kap(3,3,1),dp), "# Real(kap(3,3,1),dp)" 
Write(io_L,127) 3,3,2,Real(kap(3,3,2),dp), "# Real(kap(3,3,2),dp)" 
Write(io_L,127) 3,3,3,Real(kap(3,3,3),dp), "# Real(kap(3,3,3),dp)" 
If (MaxVal(Abs(AImag(kap))).gt.0._dp) Then 
Write(io_L,106) "Block IM{NMSSMRUN, 2} Q=",Q,"# (SUSY Scale)" 
Write(io_L,127) 1,1,1,Aimag(kap(1,1,1)), "# Aimag(kap(1,1,1))" 
Write(io_L,127) 1,1,2,Aimag(kap(1,1,2)), "# Aimag(kap(1,1,2))" 
Write(io_L,127) 1,1,3,Aimag(kap(1,1,3)), "# Aimag(kap(1,1,3))" 
Write(io_L,127) 1,2,1,Aimag(kap(1,2,1)), "# Aimag(kap(1,2,1))" 
Write(io_L,127) 1,2,2,Aimag(kap(1,2,2)), "# Aimag(kap(1,2,2))" 
Write(io_L,127) 1,2,3,Aimag(kap(1,2,3)), "# Aimag(kap(1,2,3))" 
Write(io_L,127) 1,3,1,Aimag(kap(1,3,1)), "# Aimag(kap(1,3,1))" 
Write(io_L,127) 1,3,2,Aimag(kap(1,3,2)), "# Aimag(kap(1,3,2))" 
Write(io_L,127) 1,3,3,Aimag(kap(1,3,3)), "# Aimag(kap(1,3,3))" 
Write(io_L,127) 2,1,1,Aimag(kap(2,1,1)), "# Aimag(kap(2,1,1))" 
Write(io_L,127) 2,1,2,Aimag(kap(2,1,2)), "# Aimag(kap(2,1,2))" 
Write(io_L,127) 2,1,3,Aimag(kap(2,1,3)), "# Aimag(kap(2,1,3))" 
Write(io_L,127) 2,2,1,Aimag(kap(2,2,1)), "# Aimag(kap(2,2,1))" 
Write(io_L,127) 2,2,2,Aimag(kap(2,2,2)), "# Aimag(kap(2,2,2))" 
Write(io_L,127) 2,2,3,Aimag(kap(2,2,3)), "# Aimag(kap(2,2,3))" 
Write(io_L,127) 2,3,1,Aimag(kap(2,3,1)), "# Aimag(kap(2,3,1))" 
Write(io_L,127) 2,3,2,Aimag(kap(2,3,2)), "# Aimag(kap(2,3,2))" 
Write(io_L,127) 2,3,3,Aimag(kap(2,3,3)), "# Aimag(kap(2,3,3))" 
Write(io_L,127) 3,1,1,Aimag(kap(3,1,1)), "# Aimag(kap(3,1,1))" 
Write(io_L,127) 3,1,2,Aimag(kap(3,1,2)), "# Aimag(kap(3,1,2))" 
Write(io_L,127) 3,1,3,Aimag(kap(3,1,3)), "# Aimag(kap(3,1,3))" 
Write(io_L,127) 3,2,1,Aimag(kap(3,2,1)), "# Aimag(kap(3,2,1))" 
Write(io_L,127) 3,2,2,Aimag(kap(3,2,2)), "# Aimag(kap(3,2,2))" 
Write(io_L,127) 3,2,3,Aimag(kap(3,2,3)), "# Aimag(kap(3,2,3))" 
Write(io_L,127) 3,3,1,Aimag(kap(3,3,1)), "# Aimag(kap(3,3,1))" 
Write(io_L,127) 3,3,2,Aimag(kap(3,3,2)), "# Aimag(kap(3,3,2))" 
Write(io_L,127) 3,3,3,Aimag(kap(3,3,3)), "# Aimag(kap(3,3,3))" 
End If 

Write(io_L,106) "Block Td Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Td(1,1),dp), "# Real(Td(1,1),dp)" 
Write(io_L,107)1,2,Real(Td(1,2),dp), "# Real(Td(1,2),dp)" 
Write(io_L,107)1,3,Real(Td(1,3),dp), "# Real(Td(1,3),dp)" 
Write(io_L,107)2,1,Real(Td(2,1),dp), "# Real(Td(2,1),dp)" 
Write(io_L,107)2,2,Real(Td(2,2),dp), "# Real(Td(2,2),dp)" 
Write(io_L,107)2,3,Real(Td(2,3),dp), "# Real(Td(2,3),dp)" 
Write(io_L,107)3,1,Real(Td(3,1),dp), "# Real(Td(3,1),dp)" 
Write(io_L,107)3,2,Real(Td(3,2),dp), "# Real(Td(3,2),dp)" 
Write(io_L,107)3,3,Real(Td(3,3),dp), "# Real(Td(3,3),dp)" 
If (MaxVal(Abs(AImag(Td))).gt.0._dp) Then 
Write(io_L,106) "Block IMTd Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Td(1,1)), "# Aimag(Td(1,1))" 
Write(io_L,107)1,2,Aimag(Td(1,2)), "# Aimag(Td(1,2))" 
Write(io_L,107)1,3,Aimag(Td(1,3)), "# Aimag(Td(1,3))" 
Write(io_L,107)2,1,Aimag(Td(2,1)), "# Aimag(Td(2,1))" 
Write(io_L,107)2,2,Aimag(Td(2,2)), "# Aimag(Td(2,2))" 
Write(io_L,107)2,3,Aimag(Td(2,3)), "# Aimag(Td(2,3))" 
Write(io_L,107)3,1,Aimag(Td(3,1)), "# Aimag(Td(3,1))" 
Write(io_L,107)3,2,Aimag(Td(3,2)), "# Aimag(Td(3,2))" 
Write(io_L,107)3,3,Aimag(Td(3,3)), "# Aimag(Td(3,3))" 
End If 

Write(io_L,106) "Block Te Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Te(1,1),dp), "# Real(Te(1,1),dp)" 
Write(io_L,107)1,2,Real(Te(1,2),dp), "# Real(Te(1,2),dp)" 
Write(io_L,107)1,3,Real(Te(1,3),dp), "# Real(Te(1,3),dp)" 
Write(io_L,107)2,1,Real(Te(2,1),dp), "# Real(Te(2,1),dp)" 
Write(io_L,107)2,2,Real(Te(2,2),dp), "# Real(Te(2,2),dp)" 
Write(io_L,107)2,3,Real(Te(2,3),dp), "# Real(Te(2,3),dp)" 
Write(io_L,107)3,1,Real(Te(3,1),dp), "# Real(Te(3,1),dp)" 
Write(io_L,107)3,2,Real(Te(3,2),dp), "# Real(Te(3,2),dp)" 
Write(io_L,107)3,3,Real(Te(3,3),dp), "# Real(Te(3,3),dp)" 
If (MaxVal(Abs(AImag(Te))).gt.0._dp) Then 
Write(io_L,106) "Block IMTe Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Te(1,1)), "# Aimag(Te(1,1))" 
Write(io_L,107)1,2,Aimag(Te(1,2)), "# Aimag(Te(1,2))" 
Write(io_L,107)1,3,Aimag(Te(1,3)), "# Aimag(Te(1,3))" 
Write(io_L,107)2,1,Aimag(Te(2,1)), "# Aimag(Te(2,1))" 
Write(io_L,107)2,2,Aimag(Te(2,2)), "# Aimag(Te(2,2))" 
Write(io_L,107)2,3,Aimag(Te(2,3)), "# Aimag(Te(2,3))" 
Write(io_L,107)3,1,Aimag(Te(3,1)), "# Aimag(Te(3,1))" 
Write(io_L,107)3,2,Aimag(Te(3,2)), "# Aimag(Te(3,2))" 
Write(io_L,107)3,3,Aimag(Te(3,3)), "# Aimag(Te(3,3))" 
End If 

Write(io_L,106) "Block {NMSSMRUN, 3} Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Real(Tlam(1),dp), "# Real(Tlam(1) ,dp)" 
Write(io_L,104) 2,Real(Tlam(2),dp), "# Real(Tlam(2) ,dp)" 
Write(io_L,104) 3,Real(Tlam(3),dp), "# Real(Tlam(3) ,dp)" 
If (MaxVal(Abs(AImag(Tlam))).gt.0._dp) Then 
Write(io_L,106) "Block IM{NMSSMRUN, 3} Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,Aimag(Tlam(1)), "# Aimag(Tlam(1) )" 
Write(io_L,104) 2,Aimag(Tlam(2)), "# Aimag(Tlam(2) )" 
Write(io_L,104) 3,Aimag(Tlam(3)), "# Aimag(Tlam(3) )" 
End If 

Write(io_L,106) "Block Tv Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Tv(1,1),dp), "# Real(Tv(1,1),dp)" 
Write(io_L,107)1,2,Real(Tv(1,2),dp), "# Real(Tv(1,2),dp)" 
Write(io_L,107)1,3,Real(Tv(1,3),dp), "# Real(Tv(1,3),dp)" 
Write(io_L,107)2,1,Real(Tv(2,1),dp), "# Real(Tv(2,1),dp)" 
Write(io_L,107)2,2,Real(Tv(2,2),dp), "# Real(Tv(2,2),dp)" 
Write(io_L,107)2,3,Real(Tv(2,3),dp), "# Real(Tv(2,3),dp)" 
Write(io_L,107)3,1,Real(Tv(3,1),dp), "# Real(Tv(3,1),dp)" 
Write(io_L,107)3,2,Real(Tv(3,2),dp), "# Real(Tv(3,2),dp)" 
Write(io_L,107)3,3,Real(Tv(3,3),dp), "# Real(Tv(3,3),dp)" 
If (MaxVal(Abs(AImag(Tv))).gt.0._dp) Then 
Write(io_L,106) "Block IMTv Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Tv(1,1)), "# Aimag(Tv(1,1))" 
Write(io_L,107)1,2,Aimag(Tv(1,2)), "# Aimag(Tv(1,2))" 
Write(io_L,107)1,3,Aimag(Tv(1,3)), "# Aimag(Tv(1,3))" 
Write(io_L,107)2,1,Aimag(Tv(2,1)), "# Aimag(Tv(2,1))" 
Write(io_L,107)2,2,Aimag(Tv(2,2)), "# Aimag(Tv(2,2))" 
Write(io_L,107)2,3,Aimag(Tv(2,3)), "# Aimag(Tv(2,3))" 
Write(io_L,107)3,1,Aimag(Tv(3,1)), "# Aimag(Tv(3,1))" 
Write(io_L,107)3,2,Aimag(Tv(3,2)), "# Aimag(Tv(3,2))" 
Write(io_L,107)3,3,Aimag(Tv(3,3)), "# Aimag(Tv(3,3))" 
End If 

Write(io_L,106) "Block Tu Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(Tu(1,1),dp), "# Real(Tu(1,1),dp)" 
Write(io_L,107)1,2,Real(Tu(1,2),dp), "# Real(Tu(1,2),dp)" 
Write(io_L,107)1,3,Real(Tu(1,3),dp), "# Real(Tu(1,3),dp)" 
Write(io_L,107)2,1,Real(Tu(2,1),dp), "# Real(Tu(2,1),dp)" 
Write(io_L,107)2,2,Real(Tu(2,2),dp), "# Real(Tu(2,2),dp)" 
Write(io_L,107)2,3,Real(Tu(2,3),dp), "# Real(Tu(2,3),dp)" 
Write(io_L,107)3,1,Real(Tu(3,1),dp), "# Real(Tu(3,1),dp)" 
Write(io_L,107)3,2,Real(Tu(3,2),dp), "# Real(Tu(3,2),dp)" 
Write(io_L,107)3,3,Real(Tu(3,3),dp), "# Real(Tu(3,3),dp)" 
If (MaxVal(Abs(AImag(Tu))).gt.0._dp) Then 
Write(io_L,106) "Block IMTu Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(Tu(1,1)), "# Aimag(Tu(1,1))" 
Write(io_L,107)1,2,Aimag(Tu(1,2)), "# Aimag(Tu(1,2))" 
Write(io_L,107)1,3,Aimag(Tu(1,3)), "# Aimag(Tu(1,3))" 
Write(io_L,107)2,1,Aimag(Tu(2,1)), "# Aimag(Tu(2,1))" 
Write(io_L,107)2,2,Aimag(Tu(2,2)), "# Aimag(Tu(2,2))" 
Write(io_L,107)2,3,Aimag(Tu(2,3)), "# Aimag(Tu(2,3))" 
Write(io_L,107)3,1,Aimag(Tu(3,1)), "# Aimag(Tu(3,1))" 
Write(io_L,107)3,2,Aimag(Tu(3,2)), "# Aimag(Tu(3,2))" 
Write(io_L,107)3,3,Aimag(Tu(3,3)), "# Aimag(Tu(3,3))" 
End If 

Write(io_L,106) "Block {NMSSMRUN, 4} Q=",Q,"# (SUSY Scale)" 
Write(io_L,127) 1,1,1,Real(Tk(1,1,1),dp), "# Real(Tk(1,1,1),dp)" 
Write(io_L,127) 1,1,2,Real(Tk(1,1,2),dp), "# Real(Tk(1,1,2),dp)" 
Write(io_L,127) 1,1,3,Real(Tk(1,1,3),dp), "# Real(Tk(1,1,3),dp)" 
Write(io_L,127) 1,2,1,Real(Tk(1,2,1),dp), "# Real(Tk(1,2,1),dp)" 
Write(io_L,127) 1,2,2,Real(Tk(1,2,2),dp), "# Real(Tk(1,2,2),dp)" 
Write(io_L,127) 1,2,3,Real(Tk(1,2,3),dp), "# Real(Tk(1,2,3),dp)" 
Write(io_L,127) 1,3,1,Real(Tk(1,3,1),dp), "# Real(Tk(1,3,1),dp)" 
Write(io_L,127) 1,3,2,Real(Tk(1,3,2),dp), "# Real(Tk(1,3,2),dp)" 
Write(io_L,127) 1,3,3,Real(Tk(1,3,3),dp), "# Real(Tk(1,3,3),dp)" 
Write(io_L,127) 2,1,1,Real(Tk(2,1,1),dp), "# Real(Tk(2,1,1),dp)" 
Write(io_L,127) 2,1,2,Real(Tk(2,1,2),dp), "# Real(Tk(2,1,2),dp)" 
Write(io_L,127) 2,1,3,Real(Tk(2,1,3),dp), "# Real(Tk(2,1,3),dp)" 
Write(io_L,127) 2,2,1,Real(Tk(2,2,1),dp), "# Real(Tk(2,2,1),dp)" 
Write(io_L,127) 2,2,2,Real(Tk(2,2,2),dp), "# Real(Tk(2,2,2),dp)" 
Write(io_L,127) 2,2,3,Real(Tk(2,2,3),dp), "# Real(Tk(2,2,3),dp)" 
Write(io_L,127) 2,3,1,Real(Tk(2,3,1),dp), "# Real(Tk(2,3,1),dp)" 
Write(io_L,127) 2,3,2,Real(Tk(2,3,2),dp), "# Real(Tk(2,3,2),dp)" 
Write(io_L,127) 2,3,3,Real(Tk(2,3,3),dp), "# Real(Tk(2,3,3),dp)" 
Write(io_L,127) 3,1,1,Real(Tk(3,1,1),dp), "# Real(Tk(3,1,1),dp)" 
Write(io_L,127) 3,1,2,Real(Tk(3,1,2),dp), "# Real(Tk(3,1,2),dp)" 
Write(io_L,127) 3,1,3,Real(Tk(3,1,3),dp), "# Real(Tk(3,1,3),dp)" 
Write(io_L,127) 3,2,1,Real(Tk(3,2,1),dp), "# Real(Tk(3,2,1),dp)" 
Write(io_L,127) 3,2,2,Real(Tk(3,2,2),dp), "# Real(Tk(3,2,2),dp)" 
Write(io_L,127) 3,2,3,Real(Tk(3,2,3),dp), "# Real(Tk(3,2,3),dp)" 
Write(io_L,127) 3,3,1,Real(Tk(3,3,1),dp), "# Real(Tk(3,3,1),dp)" 
Write(io_L,127) 3,3,2,Real(Tk(3,3,2),dp), "# Real(Tk(3,3,2),dp)" 
Write(io_L,127) 3,3,3,Real(Tk(3,3,3),dp), "# Real(Tk(3,3,3),dp)" 
If (MaxVal(Abs(AImag(Tk))).gt.0._dp) Then 
Write(io_L,106) "Block IM{NMSSMRUN, 4} Q=",Q,"# (SUSY Scale)" 
Write(io_L,127) 1,1,1,Aimag(Tk(1,1,1)), "# Aimag(Tk(1,1,1))" 
Write(io_L,127) 1,1,2,Aimag(Tk(1,1,2)), "# Aimag(Tk(1,1,2))" 
Write(io_L,127) 1,1,3,Aimag(Tk(1,1,3)), "# Aimag(Tk(1,1,3))" 
Write(io_L,127) 1,2,1,Aimag(Tk(1,2,1)), "# Aimag(Tk(1,2,1))" 
Write(io_L,127) 1,2,2,Aimag(Tk(1,2,2)), "# Aimag(Tk(1,2,2))" 
Write(io_L,127) 1,2,3,Aimag(Tk(1,2,3)), "# Aimag(Tk(1,2,3))" 
Write(io_L,127) 1,3,1,Aimag(Tk(1,3,1)), "# Aimag(Tk(1,3,1))" 
Write(io_L,127) 1,3,2,Aimag(Tk(1,3,2)), "# Aimag(Tk(1,3,2))" 
Write(io_L,127) 1,3,3,Aimag(Tk(1,3,3)), "# Aimag(Tk(1,3,3))" 
Write(io_L,127) 2,1,1,Aimag(Tk(2,1,1)), "# Aimag(Tk(2,1,1))" 
Write(io_L,127) 2,1,2,Aimag(Tk(2,1,2)), "# Aimag(Tk(2,1,2))" 
Write(io_L,127) 2,1,3,Aimag(Tk(2,1,3)), "# Aimag(Tk(2,1,3))" 
Write(io_L,127) 2,2,1,Aimag(Tk(2,2,1)), "# Aimag(Tk(2,2,1))" 
Write(io_L,127) 2,2,2,Aimag(Tk(2,2,2)), "# Aimag(Tk(2,2,2))" 
Write(io_L,127) 2,2,3,Aimag(Tk(2,2,3)), "# Aimag(Tk(2,2,3))" 
Write(io_L,127) 2,3,1,Aimag(Tk(2,3,1)), "# Aimag(Tk(2,3,1))" 
Write(io_L,127) 2,3,2,Aimag(Tk(2,3,2)), "# Aimag(Tk(2,3,2))" 
Write(io_L,127) 2,3,3,Aimag(Tk(2,3,3)), "# Aimag(Tk(2,3,3))" 
Write(io_L,127) 3,1,1,Aimag(Tk(3,1,1)), "# Aimag(Tk(3,1,1))" 
Write(io_L,127) 3,1,2,Aimag(Tk(3,1,2)), "# Aimag(Tk(3,1,2))" 
Write(io_L,127) 3,1,3,Aimag(Tk(3,1,3)), "# Aimag(Tk(3,1,3))" 
Write(io_L,127) 3,2,1,Aimag(Tk(3,2,1)), "# Aimag(Tk(3,2,1))" 
Write(io_L,127) 3,2,2,Aimag(Tk(3,2,2)), "# Aimag(Tk(3,2,2))" 
Write(io_L,127) 3,2,3,Aimag(Tk(3,2,3)), "# Aimag(Tk(3,2,3))" 
Write(io_L,127) 3,3,1,Aimag(Tk(3,3,1)), "# Aimag(Tk(3,3,1))" 
Write(io_L,127) 3,3,2,Aimag(Tk(3,3,2)), "# Aimag(Tk(3,3,2))" 
Write(io_L,127) 3,3,3,Aimag(Tk(3,3,3)), "# Aimag(Tk(3,3,3))" 
End If 

Write(io_L,106) "Block MSQ2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(mq2(1,1),dp), "# Real(mq2(1,1),dp)" 
Write(io_L,107)1,2,Real(mq2(1,2),dp), "# Real(mq2(1,2),dp)" 
Write(io_L,107)1,3,Real(mq2(1,3),dp), "# Real(mq2(1,3),dp)" 
Write(io_L,107)2,1,Real(mq2(2,1),dp), "# Real(mq2(2,1),dp)" 
Write(io_L,107)2,2,Real(mq2(2,2),dp), "# Real(mq2(2,2),dp)" 
Write(io_L,107)2,3,Real(mq2(2,3),dp), "# Real(mq2(2,3),dp)" 
Write(io_L,107)3,1,Real(mq2(3,1),dp), "# Real(mq2(3,1),dp)" 
Write(io_L,107)3,2,Real(mq2(3,2),dp), "# Real(mq2(3,2),dp)" 
Write(io_L,107)3,3,Real(mq2(3,3),dp), "# Real(mq2(3,3),dp)" 
If (MaxVal(Abs(AImag(mq2))).gt.0._dp) Then 
Write(io_L,106) "Block IMMSQ2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(mq2(1,1)), "# Aimag(mq2(1,1))" 
Write(io_L,107)1,2,Aimag(mq2(1,2)), "# Aimag(mq2(1,2))" 
Write(io_L,107)1,3,Aimag(mq2(1,3)), "# Aimag(mq2(1,3))" 
Write(io_L,107)2,1,Aimag(mq2(2,1)), "# Aimag(mq2(2,1))" 
Write(io_L,107)2,2,Aimag(mq2(2,2)), "# Aimag(mq2(2,2))" 
Write(io_L,107)2,3,Aimag(mq2(2,3)), "# Aimag(mq2(2,3))" 
Write(io_L,107)3,1,Aimag(mq2(3,1)), "# Aimag(mq2(3,1))" 
Write(io_L,107)3,2,Aimag(mq2(3,2)), "# Aimag(mq2(3,2))" 
Write(io_L,107)3,3,Aimag(mq2(3,3)), "# Aimag(mq2(3,3))" 
End If 

Write(io_L,106) "Block MSL2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(ml2(1,1),dp), "# Real(ml2(1,1),dp)" 
Write(io_L,107)1,2,Real(ml2(1,2),dp), "# Real(ml2(1,2),dp)" 
Write(io_L,107)1,3,Real(ml2(1,3),dp), "# Real(ml2(1,3),dp)" 
Write(io_L,107)2,1,Real(ml2(2,1),dp), "# Real(ml2(2,1),dp)" 
Write(io_L,107)2,2,Real(ml2(2,2),dp), "# Real(ml2(2,2),dp)" 
Write(io_L,107)2,3,Real(ml2(2,3),dp), "# Real(ml2(2,3),dp)" 
Write(io_L,107)3,1,Real(ml2(3,1),dp), "# Real(ml2(3,1),dp)" 
Write(io_L,107)3,2,Real(ml2(3,2),dp), "# Real(ml2(3,2),dp)" 
Write(io_L,107)3,3,Real(ml2(3,3),dp), "# Real(ml2(3,3),dp)" 
If (MaxVal(Abs(AImag(ml2))).gt.0._dp) Then 
Write(io_L,106) "Block IMMSL2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(ml2(1,1)), "# Aimag(ml2(1,1))" 
Write(io_L,107)1,2,Aimag(ml2(1,2)), "# Aimag(ml2(1,2))" 
Write(io_L,107)1,3,Aimag(ml2(1,3)), "# Aimag(ml2(1,3))" 
Write(io_L,107)2,1,Aimag(ml2(2,1)), "# Aimag(ml2(2,1))" 
Write(io_L,107)2,2,Aimag(ml2(2,2)), "# Aimag(ml2(2,2))" 
Write(io_L,107)2,3,Aimag(ml2(2,3)), "# Aimag(ml2(2,3))" 
Write(io_L,107)3,1,Aimag(ml2(3,1)), "# Aimag(ml2(3,1))" 
Write(io_L,107)3,2,Aimag(ml2(3,2)), "# Aimag(ml2(3,2))" 
Write(io_L,107)3,3,Aimag(ml2(3,3)), "# Aimag(ml2(3,3))" 
End If 

Write(io_L,106) "Block MSD2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(md2(1,1),dp), "# Real(md2(1,1),dp)" 
Write(io_L,107)1,2,Real(md2(1,2),dp), "# Real(md2(1,2),dp)" 
Write(io_L,107)1,3,Real(md2(1,3),dp), "# Real(md2(1,3),dp)" 
Write(io_L,107)2,1,Real(md2(2,1),dp), "# Real(md2(2,1),dp)" 
Write(io_L,107)2,2,Real(md2(2,2),dp), "# Real(md2(2,2),dp)" 
Write(io_L,107)2,3,Real(md2(2,3),dp), "# Real(md2(2,3),dp)" 
Write(io_L,107)3,1,Real(md2(3,1),dp), "# Real(md2(3,1),dp)" 
Write(io_L,107)3,2,Real(md2(3,2),dp), "# Real(md2(3,2),dp)" 
Write(io_L,107)3,3,Real(md2(3,3),dp), "# Real(md2(3,3),dp)" 
If (MaxVal(Abs(AImag(md2))).gt.0._dp) Then 
Write(io_L,106) "Block IMMSD2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(md2(1,1)), "# Aimag(md2(1,1))" 
Write(io_L,107)1,2,Aimag(md2(1,2)), "# Aimag(md2(1,2))" 
Write(io_L,107)1,3,Aimag(md2(1,3)), "# Aimag(md2(1,3))" 
Write(io_L,107)2,1,Aimag(md2(2,1)), "# Aimag(md2(2,1))" 
Write(io_L,107)2,2,Aimag(md2(2,2)), "# Aimag(md2(2,2))" 
Write(io_L,107)2,3,Aimag(md2(2,3)), "# Aimag(md2(2,3))" 
Write(io_L,107)3,1,Aimag(md2(3,1)), "# Aimag(md2(3,1))" 
Write(io_L,107)3,2,Aimag(md2(3,2)), "# Aimag(md2(3,2))" 
Write(io_L,107)3,3,Aimag(md2(3,3)), "# Aimag(md2(3,3))" 
End If 

Write(io_L,106) "Block MSU2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(mu2(1,1),dp), "# Real(mu2(1,1),dp)" 
Write(io_L,107)1,2,Real(mu2(1,2),dp), "# Real(mu2(1,2),dp)" 
Write(io_L,107)1,3,Real(mu2(1,3),dp), "# Real(mu2(1,3),dp)" 
Write(io_L,107)2,1,Real(mu2(2,1),dp), "# Real(mu2(2,1),dp)" 
Write(io_L,107)2,2,Real(mu2(2,2),dp), "# Real(mu2(2,2),dp)" 
Write(io_L,107)2,3,Real(mu2(2,3),dp), "# Real(mu2(2,3),dp)" 
Write(io_L,107)3,1,Real(mu2(3,1),dp), "# Real(mu2(3,1),dp)" 
Write(io_L,107)3,2,Real(mu2(3,2),dp), "# Real(mu2(3,2),dp)" 
Write(io_L,107)3,3,Real(mu2(3,3),dp), "# Real(mu2(3,3),dp)" 
If (MaxVal(Abs(AImag(mu2))).gt.0._dp) Then 
Write(io_L,106) "Block IMMSU2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(mu2(1,1)), "# Aimag(mu2(1,1))" 
Write(io_L,107)1,2,Aimag(mu2(1,2)), "# Aimag(mu2(1,2))" 
Write(io_L,107)1,3,Aimag(mu2(1,3)), "# Aimag(mu2(1,3))" 
Write(io_L,107)2,1,Aimag(mu2(2,1)), "# Aimag(mu2(2,1))" 
Write(io_L,107)2,2,Aimag(mu2(2,2)), "# Aimag(mu2(2,2))" 
Write(io_L,107)2,3,Aimag(mu2(2,3)), "# Aimag(mu2(2,3))" 
Write(io_L,107)3,1,Aimag(mu2(3,1)), "# Aimag(mu2(3,1))" 
Write(io_L,107)3,2,Aimag(mu2(3,2)), "# Aimag(mu2(3,2))" 
Write(io_L,107)3,3,Aimag(mu2(3,3)), "# Aimag(mu2(3,3))" 
End If 

Write(io_L,106) "Block MSE2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(me2(1,1),dp), "# Real(me2(1,1),dp)" 
Write(io_L,107)1,2,Real(me2(1,2),dp), "# Real(me2(1,2),dp)" 
Write(io_L,107)1,3,Real(me2(1,3),dp), "# Real(me2(1,3),dp)" 
Write(io_L,107)2,1,Real(me2(2,1),dp), "# Real(me2(2,1),dp)" 
Write(io_L,107)2,2,Real(me2(2,2),dp), "# Real(me2(2,2),dp)" 
Write(io_L,107)2,3,Real(me2(2,3),dp), "# Real(me2(2,3),dp)" 
Write(io_L,107)3,1,Real(me2(3,1),dp), "# Real(me2(3,1),dp)" 
Write(io_L,107)3,2,Real(me2(3,2),dp), "# Real(me2(3,2),dp)" 
Write(io_L,107)3,3,Real(me2(3,3),dp), "# Real(me2(3,3),dp)" 
If (MaxVal(Abs(AImag(me2))).gt.0._dp) Then 
Write(io_L,106) "Block IMMSE2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(me2(1,1)), "# Aimag(me2(1,1))" 
Write(io_L,107)1,2,Aimag(me2(1,2)), "# Aimag(me2(1,2))" 
Write(io_L,107)1,3,Aimag(me2(1,3)), "# Aimag(me2(1,3))" 
Write(io_L,107)2,1,Aimag(me2(2,1)), "# Aimag(me2(2,1))" 
Write(io_L,107)2,2,Aimag(me2(2,2)), "# Aimag(me2(2,2))" 
Write(io_L,107)2,3,Aimag(me2(2,3)), "# Aimag(me2(2,3))" 
Write(io_L,107)3,1,Aimag(me2(3,1)), "# Aimag(me2(3,1))" 
Write(io_L,107)3,2,Aimag(me2(3,2)), "# Aimag(me2(3,2))" 
Write(io_L,107)3,3,Aimag(me2(3,3)), "# Aimag(me2(3,3))" 
End If 

Write(io_L,106) "Block mv2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Real(mv2(1,1),dp), "# Real(mv2(1,1),dp)" 
Write(io_L,107)1,2,Real(mv2(1,2),dp), "# Real(mv2(1,2),dp)" 
Write(io_L,107)1,3,Real(mv2(1,3),dp), "# Real(mv2(1,3),dp)" 
Write(io_L,107)2,1,Real(mv2(2,1),dp), "# Real(mv2(2,1),dp)" 
Write(io_L,107)2,2,Real(mv2(2,2),dp), "# Real(mv2(2,2),dp)" 
Write(io_L,107)2,3,Real(mv2(2,3),dp), "# Real(mv2(2,3),dp)" 
Write(io_L,107)3,1,Real(mv2(3,1),dp), "# Real(mv2(3,1),dp)" 
Write(io_L,107)3,2,Real(mv2(3,2),dp), "# Real(mv2(3,2),dp)" 
Write(io_L,107)3,3,Real(mv2(3,3),dp), "# Real(mv2(3,3),dp)" 
If (MaxVal(Abs(AImag(mv2))).gt.0._dp) Then 
Write(io_L,106) "Block IMmv2 Q=",Q,"# (SUSY Scale)" 
Write(io_L,107)1,1,Aimag(mv2(1,1)), "# Aimag(mv2(1,1))" 
Write(io_L,107)1,2,Aimag(mv2(1,2)), "# Aimag(mv2(1,2))" 
Write(io_L,107)1,3,Aimag(mv2(1,3)), "# Aimag(mv2(1,3))" 
Write(io_L,107)2,1,Aimag(mv2(2,1)), "# Aimag(mv2(2,1))" 
Write(io_L,107)2,2,Aimag(mv2(2,2)), "# Aimag(mv2(2,2))" 
Write(io_L,107)2,3,Aimag(mv2(2,3)), "# Aimag(mv2(2,3))" 
Write(io_L,107)3,1,Aimag(mv2(3,1)), "# Aimag(mv2(3,1))" 
Write(io_L,107)3,2,Aimag(mv2(3,2)), "# Aimag(mv2(3,2))" 
Write(io_L,107)3,3,Aimag(mv2(3,3)), "# Aimag(mv2(3,3))" 
End If 

Write(io_L,106) "Block RVM2LH1 Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,mlHd2(1), "# mlHd2(1) " 
Write(io_L,104) 2,mlHd2(2), "# mlHd2(2) " 
Write(io_L,104) 3,mlHd2(3), "# mlHd2(3) " 
Write(io_L,106) "Block RIGHTVEV Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,vR(1), "# vR(1) " 
Write(io_L,104) 2,vR(2), "# vR(2) " 
Write(io_L,104) 3,vR(3), "# vR(3) " 
Write(io_L,106) "Block RVSNVEV Q=",Q,"# (SUSY Scale)" 
Write(io_L,104) 1,vL(1), "# vL(1) " 
Write(io_L,104) 2,vL(2), "# vL(2) " 
Write(io_L,104) 3,vL(3), "# vL(3) " 
If (WriteGUTvalues) Then 
Write(io_L,106) "Block GAUGEGUT Q=",M_GUT,"# (GUT scale)" 
Write(io_L,104) 1,Real(g1GUT,dp), "# g1" 
Write(io_L,104) 2,Real(g2GUT,dp), "# g2" 
Write(io_L,104) 3,Real(g3GUT,dp), "# g3" 
Write(io_L,106) "Block MSOFTGUT Q=",M_GUT,"# (GUT scale)" 
Write(io_L,104) 21,Real(mHd2GUT,dp), "# mHd2" 
Write(io_L,104) 22,Real(mHu2GUT,dp), "# mHu2" 
Write(io_L,104) 1,Real(M1GUT,dp), "# M1" 
Write(io_L,104) 2,Real(M2GUT,dp), "# M2" 
Write(io_L,104) 3,Real(M3GUT,dp), "# M3" 
Write(io_L,106) "Block YdGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(YdGUT(1,1),dp), "# Real(YdGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(YdGUT(1,2),dp), "# Real(YdGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(YdGUT(1,3),dp), "# Real(YdGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(YdGUT(2,1),dp), "# Real(YdGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(YdGUT(2,2),dp), "# Real(YdGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(YdGUT(2,3),dp), "# Real(YdGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(YdGUT(3,1),dp), "# Real(YdGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(YdGUT(3,2),dp), "# Real(YdGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(YdGUT(3,3),dp), "# Real(YdGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(YdGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMYdGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(YdGUT(1,1)), "# Aimag(YdGUT(1,1))" 
Write(io_L,107)1,2,Aimag(YdGUT(1,2)), "# Aimag(YdGUT(1,2))" 
Write(io_L,107)1,3,Aimag(YdGUT(1,3)), "# Aimag(YdGUT(1,3))" 
Write(io_L,107)2,1,Aimag(YdGUT(2,1)), "# Aimag(YdGUT(2,1))" 
Write(io_L,107)2,2,Aimag(YdGUT(2,2)), "# Aimag(YdGUT(2,2))" 
Write(io_L,107)2,3,Aimag(YdGUT(2,3)), "# Aimag(YdGUT(2,3))" 
Write(io_L,107)3,1,Aimag(YdGUT(3,1)), "# Aimag(YdGUT(3,1))" 
Write(io_L,107)3,2,Aimag(YdGUT(3,2)), "# Aimag(YdGUT(3,2))" 
Write(io_L,107)3,3,Aimag(YdGUT(3,3)), "# Aimag(YdGUT(3,3))" 
End If 

Write(io_L,106) "Block YeGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(YeGUT(1,1),dp), "# Real(YeGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(YeGUT(1,2),dp), "# Real(YeGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(YeGUT(1,3),dp), "# Real(YeGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(YeGUT(2,1),dp), "# Real(YeGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(YeGUT(2,2),dp), "# Real(YeGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(YeGUT(2,3),dp), "# Real(YeGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(YeGUT(3,1),dp), "# Real(YeGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(YeGUT(3,2),dp), "# Real(YeGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(YeGUT(3,3),dp), "# Real(YeGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(YeGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMYeGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(YeGUT(1,1)), "# Aimag(YeGUT(1,1))" 
Write(io_L,107)1,2,Aimag(YeGUT(1,2)), "# Aimag(YeGUT(1,2))" 
Write(io_L,107)1,3,Aimag(YeGUT(1,3)), "# Aimag(YeGUT(1,3))" 
Write(io_L,107)2,1,Aimag(YeGUT(2,1)), "# Aimag(YeGUT(2,1))" 
Write(io_L,107)2,2,Aimag(YeGUT(2,2)), "# Aimag(YeGUT(2,2))" 
Write(io_L,107)2,3,Aimag(YeGUT(2,3)), "# Aimag(YeGUT(2,3))" 
Write(io_L,107)3,1,Aimag(YeGUT(3,1)), "# Aimag(YeGUT(3,1))" 
Write(io_L,107)3,2,Aimag(YeGUT(3,2)), "# Aimag(YeGUT(3,2))" 
Write(io_L,107)3,3,Aimag(YeGUT(3,3)), "# Aimag(YeGUT(3,3))" 
End If 

Write(io_L,106) "Block lamGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,104) 1,Real(lamGUT(1),dp), "# Real(lamGUT(1) ,dp)" 
Write(io_L,104) 2,Real(lamGUT(2),dp), "# Real(lamGUT(2) ,dp)" 
Write(io_L,104) 3,Real(lamGUT(3),dp), "# Real(lamGUT(3) ,dp)" 
If (MaxVal(Abs(AImag(lamGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMlamGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,104) 1,Aimag(lamGUT(1)), "# Aimag(lamGUT(1) )" 
Write(io_L,104) 2,Aimag(lamGUT(2)), "# Aimag(lamGUT(2) )" 
Write(io_L,104) 3,Aimag(lamGUT(3)), "# Aimag(lamGUT(3) )" 
End If 

Write(io_L,106) "Block YvGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(YvGUT(1,1),dp), "# Real(YvGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(YvGUT(1,2),dp), "# Real(YvGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(YvGUT(1,3),dp), "# Real(YvGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(YvGUT(2,1),dp), "# Real(YvGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(YvGUT(2,2),dp), "# Real(YvGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(YvGUT(2,3),dp), "# Real(YvGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(YvGUT(3,1),dp), "# Real(YvGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(YvGUT(3,2),dp), "# Real(YvGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(YvGUT(3,3),dp), "# Real(YvGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(YvGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMYvGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(YvGUT(1,1)), "# Aimag(YvGUT(1,1))" 
Write(io_L,107)1,2,Aimag(YvGUT(1,2)), "# Aimag(YvGUT(1,2))" 
Write(io_L,107)1,3,Aimag(YvGUT(1,3)), "# Aimag(YvGUT(1,3))" 
Write(io_L,107)2,1,Aimag(YvGUT(2,1)), "# Aimag(YvGUT(2,1))" 
Write(io_L,107)2,2,Aimag(YvGUT(2,2)), "# Aimag(YvGUT(2,2))" 
Write(io_L,107)2,3,Aimag(YvGUT(2,3)), "# Aimag(YvGUT(2,3))" 
Write(io_L,107)3,1,Aimag(YvGUT(3,1)), "# Aimag(YvGUT(3,1))" 
Write(io_L,107)3,2,Aimag(YvGUT(3,2)), "# Aimag(YvGUT(3,2))" 
Write(io_L,107)3,3,Aimag(YvGUT(3,3)), "# Aimag(YvGUT(3,3))" 
End If 

Write(io_L,106) "Block YuGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(YuGUT(1,1),dp), "# Real(YuGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(YuGUT(1,2),dp), "# Real(YuGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(YuGUT(1,3),dp), "# Real(YuGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(YuGUT(2,1),dp), "# Real(YuGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(YuGUT(2,2),dp), "# Real(YuGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(YuGUT(2,3),dp), "# Real(YuGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(YuGUT(3,1),dp), "# Real(YuGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(YuGUT(3,2),dp), "# Real(YuGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(YuGUT(3,3),dp), "# Real(YuGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(YuGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMYuGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(YuGUT(1,1)), "# Aimag(YuGUT(1,1))" 
Write(io_L,107)1,2,Aimag(YuGUT(1,2)), "# Aimag(YuGUT(1,2))" 
Write(io_L,107)1,3,Aimag(YuGUT(1,3)), "# Aimag(YuGUT(1,3))" 
Write(io_L,107)2,1,Aimag(YuGUT(2,1)), "# Aimag(YuGUT(2,1))" 
Write(io_L,107)2,2,Aimag(YuGUT(2,2)), "# Aimag(YuGUT(2,2))" 
Write(io_L,107)2,3,Aimag(YuGUT(2,3)), "# Aimag(YuGUT(2,3))" 
Write(io_L,107)3,1,Aimag(YuGUT(3,1)), "# Aimag(YuGUT(3,1))" 
Write(io_L,107)3,2,Aimag(YuGUT(3,2)), "# Aimag(YuGUT(3,2))" 
Write(io_L,107)3,3,Aimag(YuGUT(3,3)), "# Aimag(YuGUT(3,3))" 
End If 

Write(io_L,106) "Block kapGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,127) 1,1,1,Real(kapGUT(1,1,1),dp), "# Real(kapGUT(1,1,1),dp)" 
Write(io_L,127) 1,1,2,Real(kapGUT(1,1,2),dp), "# Real(kapGUT(1,1,2),dp)" 
Write(io_L,127) 1,1,3,Real(kapGUT(1,1,3),dp), "# Real(kapGUT(1,1,3),dp)" 
Write(io_L,127) 1,2,1,Real(kapGUT(1,2,1),dp), "# Real(kapGUT(1,2,1),dp)" 
Write(io_L,127) 1,2,2,Real(kapGUT(1,2,2),dp), "# Real(kapGUT(1,2,2),dp)" 
Write(io_L,127) 1,2,3,Real(kapGUT(1,2,3),dp), "# Real(kapGUT(1,2,3),dp)" 
Write(io_L,127) 1,3,1,Real(kapGUT(1,3,1),dp), "# Real(kapGUT(1,3,1),dp)" 
Write(io_L,127) 1,3,2,Real(kapGUT(1,3,2),dp), "# Real(kapGUT(1,3,2),dp)" 
Write(io_L,127) 1,3,3,Real(kapGUT(1,3,3),dp), "# Real(kapGUT(1,3,3),dp)" 
Write(io_L,127) 2,1,1,Real(kapGUT(2,1,1),dp), "# Real(kapGUT(2,1,1),dp)" 
Write(io_L,127) 2,1,2,Real(kapGUT(2,1,2),dp), "# Real(kapGUT(2,1,2),dp)" 
Write(io_L,127) 2,1,3,Real(kapGUT(2,1,3),dp), "# Real(kapGUT(2,1,3),dp)" 
Write(io_L,127) 2,2,1,Real(kapGUT(2,2,1),dp), "# Real(kapGUT(2,2,1),dp)" 
Write(io_L,127) 2,2,2,Real(kapGUT(2,2,2),dp), "# Real(kapGUT(2,2,2),dp)" 
Write(io_L,127) 2,2,3,Real(kapGUT(2,2,3),dp), "# Real(kapGUT(2,2,3),dp)" 
Write(io_L,127) 2,3,1,Real(kapGUT(2,3,1),dp), "# Real(kapGUT(2,3,1),dp)" 
Write(io_L,127) 2,3,2,Real(kapGUT(2,3,2),dp), "# Real(kapGUT(2,3,2),dp)" 
Write(io_L,127) 2,3,3,Real(kapGUT(2,3,3),dp), "# Real(kapGUT(2,3,3),dp)" 
Write(io_L,127) 3,1,1,Real(kapGUT(3,1,1),dp), "# Real(kapGUT(3,1,1),dp)" 
Write(io_L,127) 3,1,2,Real(kapGUT(3,1,2),dp), "# Real(kapGUT(3,1,2),dp)" 
Write(io_L,127) 3,1,3,Real(kapGUT(3,1,3),dp), "# Real(kapGUT(3,1,3),dp)" 
Write(io_L,127) 3,2,1,Real(kapGUT(3,2,1),dp), "# Real(kapGUT(3,2,1),dp)" 
Write(io_L,127) 3,2,2,Real(kapGUT(3,2,2),dp), "# Real(kapGUT(3,2,2),dp)" 
Write(io_L,127) 3,2,3,Real(kapGUT(3,2,3),dp), "# Real(kapGUT(3,2,3),dp)" 
Write(io_L,127) 3,3,1,Real(kapGUT(3,3,1),dp), "# Real(kapGUT(3,3,1),dp)" 
Write(io_L,127) 3,3,2,Real(kapGUT(3,3,2),dp), "# Real(kapGUT(3,3,2),dp)" 
Write(io_L,127) 3,3,3,Real(kapGUT(3,3,3),dp), "# Real(kapGUT(3,3,3),dp)" 
If (MaxVal(Abs(AImag(kapGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMkapGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,127) 1,1,1,Aimag(kapGUT(1,1,1)), "# Aimag(kapGUT(1,1,1))" 
Write(io_L,127) 1,1,2,Aimag(kapGUT(1,1,2)), "# Aimag(kapGUT(1,1,2))" 
Write(io_L,127) 1,1,3,Aimag(kapGUT(1,1,3)), "# Aimag(kapGUT(1,1,3))" 
Write(io_L,127) 1,2,1,Aimag(kapGUT(1,2,1)), "# Aimag(kapGUT(1,2,1))" 
Write(io_L,127) 1,2,2,Aimag(kapGUT(1,2,2)), "# Aimag(kapGUT(1,2,2))" 
Write(io_L,127) 1,2,3,Aimag(kapGUT(1,2,3)), "# Aimag(kapGUT(1,2,3))" 
Write(io_L,127) 1,3,1,Aimag(kapGUT(1,3,1)), "# Aimag(kapGUT(1,3,1))" 
Write(io_L,127) 1,3,2,Aimag(kapGUT(1,3,2)), "# Aimag(kapGUT(1,3,2))" 
Write(io_L,127) 1,3,3,Aimag(kapGUT(1,3,3)), "# Aimag(kapGUT(1,3,3))" 
Write(io_L,127) 2,1,1,Aimag(kapGUT(2,1,1)), "# Aimag(kapGUT(2,1,1))" 
Write(io_L,127) 2,1,2,Aimag(kapGUT(2,1,2)), "# Aimag(kapGUT(2,1,2))" 
Write(io_L,127) 2,1,3,Aimag(kapGUT(2,1,3)), "# Aimag(kapGUT(2,1,3))" 
Write(io_L,127) 2,2,1,Aimag(kapGUT(2,2,1)), "# Aimag(kapGUT(2,2,1))" 
Write(io_L,127) 2,2,2,Aimag(kapGUT(2,2,2)), "# Aimag(kapGUT(2,2,2))" 
Write(io_L,127) 2,2,3,Aimag(kapGUT(2,2,3)), "# Aimag(kapGUT(2,2,3))" 
Write(io_L,127) 2,3,1,Aimag(kapGUT(2,3,1)), "# Aimag(kapGUT(2,3,1))" 
Write(io_L,127) 2,3,2,Aimag(kapGUT(2,3,2)), "# Aimag(kapGUT(2,3,2))" 
Write(io_L,127) 2,3,3,Aimag(kapGUT(2,3,3)), "# Aimag(kapGUT(2,3,3))" 
Write(io_L,127) 3,1,1,Aimag(kapGUT(3,1,1)), "# Aimag(kapGUT(3,1,1))" 
Write(io_L,127) 3,1,2,Aimag(kapGUT(3,1,2)), "# Aimag(kapGUT(3,1,2))" 
Write(io_L,127) 3,1,3,Aimag(kapGUT(3,1,3)), "# Aimag(kapGUT(3,1,3))" 
Write(io_L,127) 3,2,1,Aimag(kapGUT(3,2,1)), "# Aimag(kapGUT(3,2,1))" 
Write(io_L,127) 3,2,2,Aimag(kapGUT(3,2,2)), "# Aimag(kapGUT(3,2,2))" 
Write(io_L,127) 3,2,3,Aimag(kapGUT(3,2,3)), "# Aimag(kapGUT(3,2,3))" 
Write(io_L,127) 3,3,1,Aimag(kapGUT(3,3,1)), "# Aimag(kapGUT(3,3,1))" 
Write(io_L,127) 3,3,2,Aimag(kapGUT(3,3,2)), "# Aimag(kapGUT(3,3,2))" 
Write(io_L,127) 3,3,3,Aimag(kapGUT(3,3,3)), "# Aimag(kapGUT(3,3,3))" 
End If 

Write(io_L,106) "Block TdGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(TdGUT(1,1),dp), "# Real(TdGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(TdGUT(1,2),dp), "# Real(TdGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(TdGUT(1,3),dp), "# Real(TdGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(TdGUT(2,1),dp), "# Real(TdGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(TdGUT(2,2),dp), "# Real(TdGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(TdGUT(2,3),dp), "# Real(TdGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(TdGUT(3,1),dp), "# Real(TdGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(TdGUT(3,2),dp), "# Real(TdGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(TdGUT(3,3),dp), "# Real(TdGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(TdGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMTdGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(TdGUT(1,1)), "# Aimag(TdGUT(1,1))" 
Write(io_L,107)1,2,Aimag(TdGUT(1,2)), "# Aimag(TdGUT(1,2))" 
Write(io_L,107)1,3,Aimag(TdGUT(1,3)), "# Aimag(TdGUT(1,3))" 
Write(io_L,107)2,1,Aimag(TdGUT(2,1)), "# Aimag(TdGUT(2,1))" 
Write(io_L,107)2,2,Aimag(TdGUT(2,2)), "# Aimag(TdGUT(2,2))" 
Write(io_L,107)2,3,Aimag(TdGUT(2,3)), "# Aimag(TdGUT(2,3))" 
Write(io_L,107)3,1,Aimag(TdGUT(3,1)), "# Aimag(TdGUT(3,1))" 
Write(io_L,107)3,2,Aimag(TdGUT(3,2)), "# Aimag(TdGUT(3,2))" 
Write(io_L,107)3,3,Aimag(TdGUT(3,3)), "# Aimag(TdGUT(3,3))" 
End If 

Write(io_L,106) "Block TeGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(TeGUT(1,1),dp), "# Real(TeGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(TeGUT(1,2),dp), "# Real(TeGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(TeGUT(1,3),dp), "# Real(TeGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(TeGUT(2,1),dp), "# Real(TeGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(TeGUT(2,2),dp), "# Real(TeGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(TeGUT(2,3),dp), "# Real(TeGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(TeGUT(3,1),dp), "# Real(TeGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(TeGUT(3,2),dp), "# Real(TeGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(TeGUT(3,3),dp), "# Real(TeGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(TeGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMTeGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(TeGUT(1,1)), "# Aimag(TeGUT(1,1))" 
Write(io_L,107)1,2,Aimag(TeGUT(1,2)), "# Aimag(TeGUT(1,2))" 
Write(io_L,107)1,3,Aimag(TeGUT(1,3)), "# Aimag(TeGUT(1,3))" 
Write(io_L,107)2,1,Aimag(TeGUT(2,1)), "# Aimag(TeGUT(2,1))" 
Write(io_L,107)2,2,Aimag(TeGUT(2,2)), "# Aimag(TeGUT(2,2))" 
Write(io_L,107)2,3,Aimag(TeGUT(2,3)), "# Aimag(TeGUT(2,3))" 
Write(io_L,107)3,1,Aimag(TeGUT(3,1)), "# Aimag(TeGUT(3,1))" 
Write(io_L,107)3,2,Aimag(TeGUT(3,2)), "# Aimag(TeGUT(3,2))" 
Write(io_L,107)3,3,Aimag(TeGUT(3,3)), "# Aimag(TeGUT(3,3))" 
End If 

Write(io_L,106) "Block TlamGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,104) 1,Real(TlamGUT(1),dp), "# Real(TlamGUT(1) ,dp)" 
Write(io_L,104) 2,Real(TlamGUT(2),dp), "# Real(TlamGUT(2) ,dp)" 
Write(io_L,104) 3,Real(TlamGUT(3),dp), "# Real(TlamGUT(3) ,dp)" 
If (MaxVal(Abs(AImag(TlamGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMTlamGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,104) 1,Aimag(TlamGUT(1)), "# Aimag(TlamGUT(1) )" 
Write(io_L,104) 2,Aimag(TlamGUT(2)), "# Aimag(TlamGUT(2) )" 
Write(io_L,104) 3,Aimag(TlamGUT(3)), "# Aimag(TlamGUT(3) )" 
End If 

Write(io_L,106) "Block TvGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(TvGUT(1,1),dp), "# Real(TvGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(TvGUT(1,2),dp), "# Real(TvGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(TvGUT(1,3),dp), "# Real(TvGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(TvGUT(2,1),dp), "# Real(TvGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(TvGUT(2,2),dp), "# Real(TvGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(TvGUT(2,3),dp), "# Real(TvGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(TvGUT(3,1),dp), "# Real(TvGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(TvGUT(3,2),dp), "# Real(TvGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(TvGUT(3,3),dp), "# Real(TvGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(TvGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMTvGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(TvGUT(1,1)), "# Aimag(TvGUT(1,1))" 
Write(io_L,107)1,2,Aimag(TvGUT(1,2)), "# Aimag(TvGUT(1,2))" 
Write(io_L,107)1,3,Aimag(TvGUT(1,3)), "# Aimag(TvGUT(1,3))" 
Write(io_L,107)2,1,Aimag(TvGUT(2,1)), "# Aimag(TvGUT(2,1))" 
Write(io_L,107)2,2,Aimag(TvGUT(2,2)), "# Aimag(TvGUT(2,2))" 
Write(io_L,107)2,3,Aimag(TvGUT(2,3)), "# Aimag(TvGUT(2,3))" 
Write(io_L,107)3,1,Aimag(TvGUT(3,1)), "# Aimag(TvGUT(3,1))" 
Write(io_L,107)3,2,Aimag(TvGUT(3,2)), "# Aimag(TvGUT(3,2))" 
Write(io_L,107)3,3,Aimag(TvGUT(3,3)), "# Aimag(TvGUT(3,3))" 
End If 

Write(io_L,106) "Block TuGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(TuGUT(1,1),dp), "# Real(TuGUT(1,1),dp)" 
Write(io_L,107)1,2,Real(TuGUT(1,2),dp), "# Real(TuGUT(1,2),dp)" 
Write(io_L,107)1,3,Real(TuGUT(1,3),dp), "# Real(TuGUT(1,3),dp)" 
Write(io_L,107)2,1,Real(TuGUT(2,1),dp), "# Real(TuGUT(2,1),dp)" 
Write(io_L,107)2,2,Real(TuGUT(2,2),dp), "# Real(TuGUT(2,2),dp)" 
Write(io_L,107)2,3,Real(TuGUT(2,3),dp), "# Real(TuGUT(2,3),dp)" 
Write(io_L,107)3,1,Real(TuGUT(3,1),dp), "# Real(TuGUT(3,1),dp)" 
Write(io_L,107)3,2,Real(TuGUT(3,2),dp), "# Real(TuGUT(3,2),dp)" 
Write(io_L,107)3,3,Real(TuGUT(3,3),dp), "# Real(TuGUT(3,3),dp)" 
If (MaxVal(Abs(AImag(TuGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMTuGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(TuGUT(1,1)), "# Aimag(TuGUT(1,1))" 
Write(io_L,107)1,2,Aimag(TuGUT(1,2)), "# Aimag(TuGUT(1,2))" 
Write(io_L,107)1,3,Aimag(TuGUT(1,3)), "# Aimag(TuGUT(1,3))" 
Write(io_L,107)2,1,Aimag(TuGUT(2,1)), "# Aimag(TuGUT(2,1))" 
Write(io_L,107)2,2,Aimag(TuGUT(2,2)), "# Aimag(TuGUT(2,2))" 
Write(io_L,107)2,3,Aimag(TuGUT(2,3)), "# Aimag(TuGUT(2,3))" 
Write(io_L,107)3,1,Aimag(TuGUT(3,1)), "# Aimag(TuGUT(3,1))" 
Write(io_L,107)3,2,Aimag(TuGUT(3,2)), "# Aimag(TuGUT(3,2))" 
Write(io_L,107)3,3,Aimag(TuGUT(3,3)), "# Aimag(TuGUT(3,3))" 
End If 

Write(io_L,106) "Block TkGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,127) 1,1,1,Real(TkGUT(1,1,1),dp), "# Real(TkGUT(1,1,1),dp)" 
Write(io_L,127) 1,1,2,Real(TkGUT(1,1,2),dp), "# Real(TkGUT(1,1,2),dp)" 
Write(io_L,127) 1,1,3,Real(TkGUT(1,1,3),dp), "# Real(TkGUT(1,1,3),dp)" 
Write(io_L,127) 1,2,1,Real(TkGUT(1,2,1),dp), "# Real(TkGUT(1,2,1),dp)" 
Write(io_L,127) 1,2,2,Real(TkGUT(1,2,2),dp), "# Real(TkGUT(1,2,2),dp)" 
Write(io_L,127) 1,2,3,Real(TkGUT(1,2,3),dp), "# Real(TkGUT(1,2,3),dp)" 
Write(io_L,127) 1,3,1,Real(TkGUT(1,3,1),dp), "# Real(TkGUT(1,3,1),dp)" 
Write(io_L,127) 1,3,2,Real(TkGUT(1,3,2),dp), "# Real(TkGUT(1,3,2),dp)" 
Write(io_L,127) 1,3,3,Real(TkGUT(1,3,3),dp), "# Real(TkGUT(1,3,3),dp)" 
Write(io_L,127) 2,1,1,Real(TkGUT(2,1,1),dp), "# Real(TkGUT(2,1,1),dp)" 
Write(io_L,127) 2,1,2,Real(TkGUT(2,1,2),dp), "# Real(TkGUT(2,1,2),dp)" 
Write(io_L,127) 2,1,3,Real(TkGUT(2,1,3),dp), "# Real(TkGUT(2,1,3),dp)" 
Write(io_L,127) 2,2,1,Real(TkGUT(2,2,1),dp), "# Real(TkGUT(2,2,1),dp)" 
Write(io_L,127) 2,2,2,Real(TkGUT(2,2,2),dp), "# Real(TkGUT(2,2,2),dp)" 
Write(io_L,127) 2,2,3,Real(TkGUT(2,2,3),dp), "# Real(TkGUT(2,2,3),dp)" 
Write(io_L,127) 2,3,1,Real(TkGUT(2,3,1),dp), "# Real(TkGUT(2,3,1),dp)" 
Write(io_L,127) 2,3,2,Real(TkGUT(2,3,2),dp), "# Real(TkGUT(2,3,2),dp)" 
Write(io_L,127) 2,3,3,Real(TkGUT(2,3,3),dp), "# Real(TkGUT(2,3,3),dp)" 
Write(io_L,127) 3,1,1,Real(TkGUT(3,1,1),dp), "# Real(TkGUT(3,1,1),dp)" 
Write(io_L,127) 3,1,2,Real(TkGUT(3,1,2),dp), "# Real(TkGUT(3,1,2),dp)" 
Write(io_L,127) 3,1,3,Real(TkGUT(3,1,3),dp), "# Real(TkGUT(3,1,3),dp)" 
Write(io_L,127) 3,2,1,Real(TkGUT(3,2,1),dp), "# Real(TkGUT(3,2,1),dp)" 
Write(io_L,127) 3,2,2,Real(TkGUT(3,2,2),dp), "# Real(TkGUT(3,2,2),dp)" 
Write(io_L,127) 3,2,3,Real(TkGUT(3,2,3),dp), "# Real(TkGUT(3,2,3),dp)" 
Write(io_L,127) 3,3,1,Real(TkGUT(3,3,1),dp), "# Real(TkGUT(3,3,1),dp)" 
Write(io_L,127) 3,3,2,Real(TkGUT(3,3,2),dp), "# Real(TkGUT(3,3,2),dp)" 
Write(io_L,127) 3,3,3,Real(TkGUT(3,3,3),dp), "# Real(TkGUT(3,3,3),dp)" 
If (MaxVal(Abs(AImag(TkGUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMTkGUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,127) 1,1,1,Aimag(TkGUT(1,1,1)), "# Aimag(TkGUT(1,1,1))" 
Write(io_L,127) 1,1,2,Aimag(TkGUT(1,1,2)), "# Aimag(TkGUT(1,1,2))" 
Write(io_L,127) 1,1,3,Aimag(TkGUT(1,1,3)), "# Aimag(TkGUT(1,1,3))" 
Write(io_L,127) 1,2,1,Aimag(TkGUT(1,2,1)), "# Aimag(TkGUT(1,2,1))" 
Write(io_L,127) 1,2,2,Aimag(TkGUT(1,2,2)), "# Aimag(TkGUT(1,2,2))" 
Write(io_L,127) 1,2,3,Aimag(TkGUT(1,2,3)), "# Aimag(TkGUT(1,2,3))" 
Write(io_L,127) 1,3,1,Aimag(TkGUT(1,3,1)), "# Aimag(TkGUT(1,3,1))" 
Write(io_L,127) 1,3,2,Aimag(TkGUT(1,3,2)), "# Aimag(TkGUT(1,3,2))" 
Write(io_L,127) 1,3,3,Aimag(TkGUT(1,3,3)), "# Aimag(TkGUT(1,3,3))" 
Write(io_L,127) 2,1,1,Aimag(TkGUT(2,1,1)), "# Aimag(TkGUT(2,1,1))" 
Write(io_L,127) 2,1,2,Aimag(TkGUT(2,1,2)), "# Aimag(TkGUT(2,1,2))" 
Write(io_L,127) 2,1,3,Aimag(TkGUT(2,1,3)), "# Aimag(TkGUT(2,1,3))" 
Write(io_L,127) 2,2,1,Aimag(TkGUT(2,2,1)), "# Aimag(TkGUT(2,2,1))" 
Write(io_L,127) 2,2,2,Aimag(TkGUT(2,2,2)), "# Aimag(TkGUT(2,2,2))" 
Write(io_L,127) 2,2,3,Aimag(TkGUT(2,2,3)), "# Aimag(TkGUT(2,2,3))" 
Write(io_L,127) 2,3,1,Aimag(TkGUT(2,3,1)), "# Aimag(TkGUT(2,3,1))" 
Write(io_L,127) 2,3,2,Aimag(TkGUT(2,3,2)), "# Aimag(TkGUT(2,3,2))" 
Write(io_L,127) 2,3,3,Aimag(TkGUT(2,3,3)), "# Aimag(TkGUT(2,3,3))" 
Write(io_L,127) 3,1,1,Aimag(TkGUT(3,1,1)), "# Aimag(TkGUT(3,1,1))" 
Write(io_L,127) 3,1,2,Aimag(TkGUT(3,1,2)), "# Aimag(TkGUT(3,1,2))" 
Write(io_L,127) 3,1,3,Aimag(TkGUT(3,1,3)), "# Aimag(TkGUT(3,1,3))" 
Write(io_L,127) 3,2,1,Aimag(TkGUT(3,2,1)), "# Aimag(TkGUT(3,2,1))" 
Write(io_L,127) 3,2,2,Aimag(TkGUT(3,2,2)), "# Aimag(TkGUT(3,2,2))" 
Write(io_L,127) 3,2,3,Aimag(TkGUT(3,2,3)), "# Aimag(TkGUT(3,2,3))" 
Write(io_L,127) 3,3,1,Aimag(TkGUT(3,3,1)), "# Aimag(TkGUT(3,3,1))" 
Write(io_L,127) 3,3,2,Aimag(TkGUT(3,3,2)), "# Aimag(TkGUT(3,3,2))" 
Write(io_L,127) 3,3,3,Aimag(TkGUT(3,3,3)), "# Aimag(TkGUT(3,3,3))" 
End If 

Write(io_L,106) "Block mq2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(mq2GUT(1,1),dp), "# Real(mq2GUT(1,1),dp)" 
Write(io_L,107)1,2,Real(mq2GUT(1,2),dp), "# Real(mq2GUT(1,2),dp)" 
Write(io_L,107)1,3,Real(mq2GUT(1,3),dp), "# Real(mq2GUT(1,3),dp)" 
Write(io_L,107)2,1,Real(mq2GUT(2,1),dp), "# Real(mq2GUT(2,1),dp)" 
Write(io_L,107)2,2,Real(mq2GUT(2,2),dp), "# Real(mq2GUT(2,2),dp)" 
Write(io_L,107)2,3,Real(mq2GUT(2,3),dp), "# Real(mq2GUT(2,3),dp)" 
Write(io_L,107)3,1,Real(mq2GUT(3,1),dp), "# Real(mq2GUT(3,1),dp)" 
Write(io_L,107)3,2,Real(mq2GUT(3,2),dp), "# Real(mq2GUT(3,2),dp)" 
Write(io_L,107)3,3,Real(mq2GUT(3,3),dp), "# Real(mq2GUT(3,3),dp)" 
If (MaxVal(Abs(AImag(mq2GUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMmq2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(mq2GUT(1,1)), "# Aimag(mq2GUT(1,1))" 
Write(io_L,107)1,2,Aimag(mq2GUT(1,2)), "# Aimag(mq2GUT(1,2))" 
Write(io_L,107)1,3,Aimag(mq2GUT(1,3)), "# Aimag(mq2GUT(1,3))" 
Write(io_L,107)2,1,Aimag(mq2GUT(2,1)), "# Aimag(mq2GUT(2,1))" 
Write(io_L,107)2,2,Aimag(mq2GUT(2,2)), "# Aimag(mq2GUT(2,2))" 
Write(io_L,107)2,3,Aimag(mq2GUT(2,3)), "# Aimag(mq2GUT(2,3))" 
Write(io_L,107)3,1,Aimag(mq2GUT(3,1)), "# Aimag(mq2GUT(3,1))" 
Write(io_L,107)3,2,Aimag(mq2GUT(3,2)), "# Aimag(mq2GUT(3,2))" 
Write(io_L,107)3,3,Aimag(mq2GUT(3,3)), "# Aimag(mq2GUT(3,3))" 
End If 

Write(io_L,106) "Block ml2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(ml2GUT(1,1),dp), "# Real(ml2GUT(1,1),dp)" 
Write(io_L,107)1,2,Real(ml2GUT(1,2),dp), "# Real(ml2GUT(1,2),dp)" 
Write(io_L,107)1,3,Real(ml2GUT(1,3),dp), "# Real(ml2GUT(1,3),dp)" 
Write(io_L,107)2,1,Real(ml2GUT(2,1),dp), "# Real(ml2GUT(2,1),dp)" 
Write(io_L,107)2,2,Real(ml2GUT(2,2),dp), "# Real(ml2GUT(2,2),dp)" 
Write(io_L,107)2,3,Real(ml2GUT(2,3),dp), "# Real(ml2GUT(2,3),dp)" 
Write(io_L,107)3,1,Real(ml2GUT(3,1),dp), "# Real(ml2GUT(3,1),dp)" 
Write(io_L,107)3,2,Real(ml2GUT(3,2),dp), "# Real(ml2GUT(3,2),dp)" 
Write(io_L,107)3,3,Real(ml2GUT(3,3),dp), "# Real(ml2GUT(3,3),dp)" 
If (MaxVal(Abs(AImag(ml2GUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMml2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(ml2GUT(1,1)), "# Aimag(ml2GUT(1,1))" 
Write(io_L,107)1,2,Aimag(ml2GUT(1,2)), "# Aimag(ml2GUT(1,2))" 
Write(io_L,107)1,3,Aimag(ml2GUT(1,3)), "# Aimag(ml2GUT(1,3))" 
Write(io_L,107)2,1,Aimag(ml2GUT(2,1)), "# Aimag(ml2GUT(2,1))" 
Write(io_L,107)2,2,Aimag(ml2GUT(2,2)), "# Aimag(ml2GUT(2,2))" 
Write(io_L,107)2,3,Aimag(ml2GUT(2,3)), "# Aimag(ml2GUT(2,3))" 
Write(io_L,107)3,1,Aimag(ml2GUT(3,1)), "# Aimag(ml2GUT(3,1))" 
Write(io_L,107)3,2,Aimag(ml2GUT(3,2)), "# Aimag(ml2GUT(3,2))" 
Write(io_L,107)3,3,Aimag(ml2GUT(3,3)), "# Aimag(ml2GUT(3,3))" 
End If 

Write(io_L,106) "Block md2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(md2GUT(1,1),dp), "# Real(md2GUT(1,1),dp)" 
Write(io_L,107)1,2,Real(md2GUT(1,2),dp), "# Real(md2GUT(1,2),dp)" 
Write(io_L,107)1,3,Real(md2GUT(1,3),dp), "# Real(md2GUT(1,3),dp)" 
Write(io_L,107)2,1,Real(md2GUT(2,1),dp), "# Real(md2GUT(2,1),dp)" 
Write(io_L,107)2,2,Real(md2GUT(2,2),dp), "# Real(md2GUT(2,2),dp)" 
Write(io_L,107)2,3,Real(md2GUT(2,3),dp), "# Real(md2GUT(2,3),dp)" 
Write(io_L,107)3,1,Real(md2GUT(3,1),dp), "# Real(md2GUT(3,1),dp)" 
Write(io_L,107)3,2,Real(md2GUT(3,2),dp), "# Real(md2GUT(3,2),dp)" 
Write(io_L,107)3,3,Real(md2GUT(3,3),dp), "# Real(md2GUT(3,3),dp)" 
If (MaxVal(Abs(AImag(md2GUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMmd2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(md2GUT(1,1)), "# Aimag(md2GUT(1,1))" 
Write(io_L,107)1,2,Aimag(md2GUT(1,2)), "# Aimag(md2GUT(1,2))" 
Write(io_L,107)1,3,Aimag(md2GUT(1,3)), "# Aimag(md2GUT(1,3))" 
Write(io_L,107)2,1,Aimag(md2GUT(2,1)), "# Aimag(md2GUT(2,1))" 
Write(io_L,107)2,2,Aimag(md2GUT(2,2)), "# Aimag(md2GUT(2,2))" 
Write(io_L,107)2,3,Aimag(md2GUT(2,3)), "# Aimag(md2GUT(2,3))" 
Write(io_L,107)3,1,Aimag(md2GUT(3,1)), "# Aimag(md2GUT(3,1))" 
Write(io_L,107)3,2,Aimag(md2GUT(3,2)), "# Aimag(md2GUT(3,2))" 
Write(io_L,107)3,3,Aimag(md2GUT(3,3)), "# Aimag(md2GUT(3,3))" 
End If 

Write(io_L,106) "Block mu2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(mu2GUT(1,1),dp), "# Real(mu2GUT(1,1),dp)" 
Write(io_L,107)1,2,Real(mu2GUT(1,2),dp), "# Real(mu2GUT(1,2),dp)" 
Write(io_L,107)1,3,Real(mu2GUT(1,3),dp), "# Real(mu2GUT(1,3),dp)" 
Write(io_L,107)2,1,Real(mu2GUT(2,1),dp), "# Real(mu2GUT(2,1),dp)" 
Write(io_L,107)2,2,Real(mu2GUT(2,2),dp), "# Real(mu2GUT(2,2),dp)" 
Write(io_L,107)2,3,Real(mu2GUT(2,3),dp), "# Real(mu2GUT(2,3),dp)" 
Write(io_L,107)3,1,Real(mu2GUT(3,1),dp), "# Real(mu2GUT(3,1),dp)" 
Write(io_L,107)3,2,Real(mu2GUT(3,2),dp), "# Real(mu2GUT(3,2),dp)" 
Write(io_L,107)3,3,Real(mu2GUT(3,3),dp), "# Real(mu2GUT(3,3),dp)" 
If (MaxVal(Abs(AImag(mu2GUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMmu2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(mu2GUT(1,1)), "# Aimag(mu2GUT(1,1))" 
Write(io_L,107)1,2,Aimag(mu2GUT(1,2)), "# Aimag(mu2GUT(1,2))" 
Write(io_L,107)1,3,Aimag(mu2GUT(1,3)), "# Aimag(mu2GUT(1,3))" 
Write(io_L,107)2,1,Aimag(mu2GUT(2,1)), "# Aimag(mu2GUT(2,1))" 
Write(io_L,107)2,2,Aimag(mu2GUT(2,2)), "# Aimag(mu2GUT(2,2))" 
Write(io_L,107)2,3,Aimag(mu2GUT(2,3)), "# Aimag(mu2GUT(2,3))" 
Write(io_L,107)3,1,Aimag(mu2GUT(3,1)), "# Aimag(mu2GUT(3,1))" 
Write(io_L,107)3,2,Aimag(mu2GUT(3,2)), "# Aimag(mu2GUT(3,2))" 
Write(io_L,107)3,3,Aimag(mu2GUT(3,3)), "# Aimag(mu2GUT(3,3))" 
End If 

Write(io_L,106) "Block me2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(me2GUT(1,1),dp), "# Real(me2GUT(1,1),dp)" 
Write(io_L,107)1,2,Real(me2GUT(1,2),dp), "# Real(me2GUT(1,2),dp)" 
Write(io_L,107)1,3,Real(me2GUT(1,3),dp), "# Real(me2GUT(1,3),dp)" 
Write(io_L,107)2,1,Real(me2GUT(2,1),dp), "# Real(me2GUT(2,1),dp)" 
Write(io_L,107)2,2,Real(me2GUT(2,2),dp), "# Real(me2GUT(2,2),dp)" 
Write(io_L,107)2,3,Real(me2GUT(2,3),dp), "# Real(me2GUT(2,3),dp)" 
Write(io_L,107)3,1,Real(me2GUT(3,1),dp), "# Real(me2GUT(3,1),dp)" 
Write(io_L,107)3,2,Real(me2GUT(3,2),dp), "# Real(me2GUT(3,2),dp)" 
Write(io_L,107)3,3,Real(me2GUT(3,3),dp), "# Real(me2GUT(3,3),dp)" 
If (MaxVal(Abs(AImag(me2GUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMme2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(me2GUT(1,1)), "# Aimag(me2GUT(1,1))" 
Write(io_L,107)1,2,Aimag(me2GUT(1,2)), "# Aimag(me2GUT(1,2))" 
Write(io_L,107)1,3,Aimag(me2GUT(1,3)), "# Aimag(me2GUT(1,3))" 
Write(io_L,107)2,1,Aimag(me2GUT(2,1)), "# Aimag(me2GUT(2,1))" 
Write(io_L,107)2,2,Aimag(me2GUT(2,2)), "# Aimag(me2GUT(2,2))" 
Write(io_L,107)2,3,Aimag(me2GUT(2,3)), "# Aimag(me2GUT(2,3))" 
Write(io_L,107)3,1,Aimag(me2GUT(3,1)), "# Aimag(me2GUT(3,1))" 
Write(io_L,107)3,2,Aimag(me2GUT(3,2)), "# Aimag(me2GUT(3,2))" 
Write(io_L,107)3,3,Aimag(me2GUT(3,3)), "# Aimag(me2GUT(3,3))" 
End If 

Write(io_L,106) "Block mv2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Real(mv2GUT(1,1),dp), "# Real(mv2GUT(1,1),dp)" 
Write(io_L,107)1,2,Real(mv2GUT(1,2),dp), "# Real(mv2GUT(1,2),dp)" 
Write(io_L,107)1,3,Real(mv2GUT(1,3),dp), "# Real(mv2GUT(1,3),dp)" 
Write(io_L,107)2,1,Real(mv2GUT(2,1),dp), "# Real(mv2GUT(2,1),dp)" 
Write(io_L,107)2,2,Real(mv2GUT(2,2),dp), "# Real(mv2GUT(2,2),dp)" 
Write(io_L,107)2,3,Real(mv2GUT(2,3),dp), "# Real(mv2GUT(2,3),dp)" 
Write(io_L,107)3,1,Real(mv2GUT(3,1),dp), "# Real(mv2GUT(3,1),dp)" 
Write(io_L,107)3,2,Real(mv2GUT(3,2),dp), "# Real(mv2GUT(3,2),dp)" 
Write(io_L,107)3,3,Real(mv2GUT(3,3),dp), "# Real(mv2GUT(3,3),dp)" 
If (MaxVal(Abs(AImag(mv2GUT))).gt.0._dp) Then 
Write(io_L,106) "Block IMmv2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,107)1,1,Aimag(mv2GUT(1,1)), "# Aimag(mv2GUT(1,1))" 
Write(io_L,107)1,2,Aimag(mv2GUT(1,2)), "# Aimag(mv2GUT(1,2))" 
Write(io_L,107)1,3,Aimag(mv2GUT(1,3)), "# Aimag(mv2GUT(1,3))" 
Write(io_L,107)2,1,Aimag(mv2GUT(2,1)), "# Aimag(mv2GUT(2,1))" 
Write(io_L,107)2,2,Aimag(mv2GUT(2,2)), "# Aimag(mv2GUT(2,2))" 
Write(io_L,107)2,3,Aimag(mv2GUT(2,3)), "# Aimag(mv2GUT(2,3))" 
Write(io_L,107)3,1,Aimag(mv2GUT(3,1)), "# Aimag(mv2GUT(3,1))" 
Write(io_L,107)3,2,Aimag(mv2GUT(3,2)), "# Aimag(mv2GUT(3,2))" 
Write(io_L,107)3,3,Aimag(mv2GUT(3,3)), "# Aimag(mv2GUT(3,3))" 
End If 

Write(io_L,106) "Block mlHd2GUT Q=",M_GUT,"# (GUT Scale)" 
Write(io_L,104) 1,mlHd2GUT(1), "# mlHd2GUT(1) " 
Write(io_L,104) 2,mlHd2GUT(2), "# mlHd2GUT(2) " 
Write(io_L,104) 3,mlHd2GUT(3), "# mlHd2GUT(3) " 
End if 
 
MassLSP = 100000._dp 
Write(io_L,100) "Block MASS  # Mass spectrum"
Write(io_L,100) "#   PDG code      mass          particle" 
 Write(io_L,102) INT(Abs(PDGSd(1))),MSd(1),"# "//Trim(NameParticleSd(1))// "" 
 Write(io_L,102) INT(Abs(PDGSd(2))),MSd(2),"# "//Trim(NameParticleSd(2))// "" 
 Write(io_L,102) INT(Abs(PDGSd(3))),MSd(3),"# "//Trim(NameParticleSd(3))// "" 
 Write(io_L,102) INT(Abs(PDGSd(4))),MSd(4),"# "//Trim(NameParticleSd(4))// "" 
 Write(io_L,102) INT(Abs(PDGSd(5))),MSd(5),"# "//Trim(NameParticleSd(5))// "" 
 Write(io_L,102) INT(Abs(PDGSd(6))),MSd(6),"# "//Trim(NameParticleSd(6))// "" 
 Write(io_L,102) INT(Abs(PDGSu(1))),MSu(1),"# "//Trim(NameParticleSu(1))// "" 
 Write(io_L,102) INT(Abs(PDGSu(2))),MSu(2),"# "//Trim(NameParticleSu(2))// "" 
 Write(io_L,102) INT(Abs(PDGSu(3))),MSu(3),"# "//Trim(NameParticleSu(3))// "" 
 Write(io_L,102) INT(Abs(PDGSu(4))),MSu(4),"# "//Trim(NameParticleSu(4))// "" 
 Write(io_L,102) INT(Abs(PDGSu(5))),MSu(5),"# "//Trim(NameParticleSu(5))// "" 
 Write(io_L,102) INT(Abs(PDGSu(6))),MSu(6),"# "//Trim(NameParticleSu(6))// "" 
 Write(io_L,102) INT(Abs(PDGhh(1))),Mhh(1),"# "//Trim(NameParticlehh(1))// "" 
 Write(io_L,102) INT(Abs(PDGhh(2))),Mhh(2),"# "//Trim(NameParticlehh(2))// "" 
 Write(io_L,102) INT(Abs(PDGhh(3))),Mhh(3),"# "//Trim(NameParticlehh(3))// "" 
 Write(io_L,102) INT(Abs(PDGhh(4))),Mhh(4),"# "//Trim(NameParticlehh(4))// "" 
 Write(io_L,102) INT(Abs(PDGhh(5))),Mhh(5),"# "//Trim(NameParticlehh(5))// "" 
 Write(io_L,102) INT(Abs(PDGhh(6))),Mhh(6),"# "//Trim(NameParticlehh(6))// "" 
 Write(io_L,102) INT(Abs(PDGhh(7))),Mhh(7),"# "//Trim(NameParticlehh(7))// "" 
 Write(io_L,102) INT(Abs(PDGhh(8))),Mhh(8),"# "//Trim(NameParticlehh(8))// "" 
 Write(io_L,102) INT(Abs(PDGAh(2))),MAh(2),"# "//Trim(NameParticleAh(2))// "" 
 Write(io_L,102) INT(Abs(PDGAh(3))),MAh(3),"# "//Trim(NameParticleAh(3))// "" 
 Write(io_L,102) INT(Abs(PDGAh(4))),MAh(4),"# "//Trim(NameParticleAh(4))// "" 
 Write(io_L,102) INT(Abs(PDGAh(5))),MAh(5),"# "//Trim(NameParticleAh(5))// "" 
 Write(io_L,102) INT(Abs(PDGAh(6))),MAh(6),"# "//Trim(NameParticleAh(6))// "" 
 Write(io_L,102) INT(Abs(PDGAh(7))),MAh(7),"# "//Trim(NameParticleAh(7))// "" 
 Write(io_L,102) INT(Abs(PDGAh(8))),MAh(8),"# "//Trim(NameParticleAh(8))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(2))),MHpm(2),"# "//Trim(NameParticleHpm(2))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(3))),MHpm(3),"# "//Trim(NameParticleHpm(3))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(4))),MHpm(4),"# "//Trim(NameParticleHpm(4))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(5))),MHpm(5),"# "//Trim(NameParticleHpm(5))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(6))),MHpm(6),"# "//Trim(NameParticleHpm(6))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(7))),MHpm(7),"# "//Trim(NameParticleHpm(7))// "" 
 Write(io_L,102) INT(Abs(PDGHpm(8))),MHpm(8),"# "//Trim(NameParticleHpm(8))// "" 
 Write(io_L,102) 23,MVZ,"# VZ" 
 Write(io_L,102) 24,MVWm,"# VWm" 
 Write(io_L,102) INT(Abs(PDGFd(1))),MFd(1),"# "//Trim(NameParticleFd(1))// "" 
 Write(io_L,102) INT(Abs(PDGFd(2))),MFd(2),"# "//Trim(NameParticleFd(2))// "" 
 Write(io_L,102) INT(Abs(PDGFd(3))),MFd(3),"# "//Trim(NameParticleFd(3))// "" 
 Write(io_L,102) INT(Abs(PDGFu(1))),MFu(1),"# "//Trim(NameParticleFu(1))// "" 
 Write(io_L,102) INT(Abs(PDGFu(2))),MFu(2),"# "//Trim(NameParticleFu(2))// "" 
 Write(io_L,102) INT(Abs(PDGFu(3))),MFu(3),"# "//Trim(NameParticleFu(3))// "" 
 Write(io_L,102) 1000021,MGlu,"# Glu" 
 Write(io_L,102) INT(Abs(PDGChi(1))),MChi(1),"# "//Trim(NameParticleChi(1))// "" 
 Write(io_L,102) INT(Abs(PDGChi(2))),MChi(2),"# "//Trim(NameParticleChi(2))// "" 
 Write(io_L,102) INT(Abs(PDGChi(3))),MChi(3),"# "//Trim(NameParticleChi(3))// "" 
 Write(io_L,102) INT(Abs(PDGChi(4))),MChi(4),"# "//Trim(NameParticleChi(4))// "" 
 Write(io_L,102) INT(Abs(PDGChi(5))),MChi(5),"# "//Trim(NameParticleChi(5))// "" 
 Write(io_L,102) INT(Abs(PDGChi(6))),MChi(6),"# "//Trim(NameParticleChi(6))// "" 
 Write(io_L,102) INT(Abs(PDGChi(7))),MChi(7),"# "//Trim(NameParticleChi(7))// "" 
 Write(io_L,102) INT(Abs(PDGChi(8))),MChi(8),"# "//Trim(NameParticleChi(8))// "" 
 Write(io_L,102) INT(Abs(PDGChi(9))),MChi(9),"# "//Trim(NameParticleChi(9))// "" 
 Write(io_L,102) INT(Abs(PDGChi(10))),MChi(10),"# "//Trim(NameParticleChi(10))// "" 
 Write(io_L,102) INT(Abs(PDGCha(1))),MCha(1),"# "//Trim(NameParticleCha(1))// "" 
 Write(io_L,102) INT(Abs(PDGCha(2))),MCha(2),"# "//Trim(NameParticleCha(2))// "" 
 Write(io_L,102) INT(Abs(PDGCha(3))),MCha(3),"# "//Trim(NameParticleCha(3))// "" 
 Write(io_L,102) INT(Abs(PDGCha(4))),MCha(4),"# "//Trim(NameParticleCha(4))// "" 
 Write(io_L,102) INT(Abs(PDGCha(5))),MCha(5),"# "//Trim(NameParticleCha(5))// "" 

 
Write(io_L,106) "Block DSQMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZD(1,1),dp), "# Real(ZD(1,1),dp)" 
Write(io_L,107)1,2,Real(ZD(1,2),dp), "# Real(ZD(1,2),dp)" 
Write(io_L,107)1,3,Real(ZD(1,3),dp), "# Real(ZD(1,3),dp)" 
Write(io_L,107)1,4,Real(ZD(1,4),dp), "# Real(ZD(1,4),dp)" 
Write(io_L,107)1,5,Real(ZD(1,5),dp), "# Real(ZD(1,5),dp)" 
Write(io_L,107)1,6,Real(ZD(1,6),dp), "# Real(ZD(1,6),dp)" 
Write(io_L,107)2,1,Real(ZD(2,1),dp), "# Real(ZD(2,1),dp)" 
Write(io_L,107)2,2,Real(ZD(2,2),dp), "# Real(ZD(2,2),dp)" 
Write(io_L,107)2,3,Real(ZD(2,3),dp), "# Real(ZD(2,3),dp)" 
Write(io_L,107)2,4,Real(ZD(2,4),dp), "# Real(ZD(2,4),dp)" 
Write(io_L,107)2,5,Real(ZD(2,5),dp), "# Real(ZD(2,5),dp)" 
Write(io_L,107)2,6,Real(ZD(2,6),dp), "# Real(ZD(2,6),dp)" 
Write(io_L,107)3,1,Real(ZD(3,1),dp), "# Real(ZD(3,1),dp)" 
Write(io_L,107)3,2,Real(ZD(3,2),dp), "# Real(ZD(3,2),dp)" 
Write(io_L,107)3,3,Real(ZD(3,3),dp), "# Real(ZD(3,3),dp)" 
Write(io_L,107)3,4,Real(ZD(3,4),dp), "# Real(ZD(3,4),dp)" 
Write(io_L,107)3,5,Real(ZD(3,5),dp), "# Real(ZD(3,5),dp)" 
Write(io_L,107)3,6,Real(ZD(3,6),dp), "# Real(ZD(3,6),dp)" 
Write(io_L,107)4,1,Real(ZD(4,1),dp), "# Real(ZD(4,1),dp)" 
Write(io_L,107)4,2,Real(ZD(4,2),dp), "# Real(ZD(4,2),dp)" 
Write(io_L,107)4,3,Real(ZD(4,3),dp), "# Real(ZD(4,3),dp)" 
Write(io_L,107)4,4,Real(ZD(4,4),dp), "# Real(ZD(4,4),dp)" 
Write(io_L,107)4,5,Real(ZD(4,5),dp), "# Real(ZD(4,5),dp)" 
Write(io_L,107)4,6,Real(ZD(4,6),dp), "# Real(ZD(4,6),dp)" 
Write(io_L,107)5,1,Real(ZD(5,1),dp), "# Real(ZD(5,1),dp)" 
Write(io_L,107)5,2,Real(ZD(5,2),dp), "# Real(ZD(5,2),dp)" 
Write(io_L,107)5,3,Real(ZD(5,3),dp), "# Real(ZD(5,3),dp)" 
Write(io_L,107)5,4,Real(ZD(5,4),dp), "# Real(ZD(5,4),dp)" 
Write(io_L,107)5,5,Real(ZD(5,5),dp), "# Real(ZD(5,5),dp)" 
Write(io_L,107)5,6,Real(ZD(5,6),dp), "# Real(ZD(5,6),dp)" 
Write(io_L,107)6,1,Real(ZD(6,1),dp), "# Real(ZD(6,1),dp)" 
Write(io_L,107)6,2,Real(ZD(6,2),dp), "# Real(ZD(6,2),dp)" 
Write(io_L,107)6,3,Real(ZD(6,3),dp), "# Real(ZD(6,3),dp)" 
Write(io_L,107)6,4,Real(ZD(6,4),dp), "# Real(ZD(6,4),dp)" 
Write(io_L,107)6,5,Real(ZD(6,5),dp), "# Real(ZD(6,5),dp)" 
Write(io_L,107)6,6,Real(ZD(6,6),dp), "# Real(ZD(6,6),dp)" 
If (MaxVal(Abs(AImag(ZD))).gt.0._dp) Then 
Write(io_L,106) "Block IMDSQMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZD(1,1)), "# Aimag(ZD(1,1))" 
Write(io_L,107)1,2,Aimag(ZD(1,2)), "# Aimag(ZD(1,2))" 
Write(io_L,107)1,3,Aimag(ZD(1,3)), "# Aimag(ZD(1,3))" 
Write(io_L,107)1,4,Aimag(ZD(1,4)), "# Aimag(ZD(1,4))" 
Write(io_L,107)1,5,Aimag(ZD(1,5)), "# Aimag(ZD(1,5))" 
Write(io_L,107)1,6,Aimag(ZD(1,6)), "# Aimag(ZD(1,6))" 
Write(io_L,107)2,1,Aimag(ZD(2,1)), "# Aimag(ZD(2,1))" 
Write(io_L,107)2,2,Aimag(ZD(2,2)), "# Aimag(ZD(2,2))" 
Write(io_L,107)2,3,Aimag(ZD(2,3)), "# Aimag(ZD(2,3))" 
Write(io_L,107)2,4,Aimag(ZD(2,4)), "# Aimag(ZD(2,4))" 
Write(io_L,107)2,5,Aimag(ZD(2,5)), "# Aimag(ZD(2,5))" 
Write(io_L,107)2,6,Aimag(ZD(2,6)), "# Aimag(ZD(2,6))" 
Write(io_L,107)3,1,Aimag(ZD(3,1)), "# Aimag(ZD(3,1))" 
Write(io_L,107)3,2,Aimag(ZD(3,2)), "# Aimag(ZD(3,2))" 
Write(io_L,107)3,3,Aimag(ZD(3,3)), "# Aimag(ZD(3,3))" 
Write(io_L,107)3,4,Aimag(ZD(3,4)), "# Aimag(ZD(3,4))" 
Write(io_L,107)3,5,Aimag(ZD(3,5)), "# Aimag(ZD(3,5))" 
Write(io_L,107)3,6,Aimag(ZD(3,6)), "# Aimag(ZD(3,6))" 
Write(io_L,107)4,1,Aimag(ZD(4,1)), "# Aimag(ZD(4,1))" 
Write(io_L,107)4,2,Aimag(ZD(4,2)), "# Aimag(ZD(4,2))" 
Write(io_L,107)4,3,Aimag(ZD(4,3)), "# Aimag(ZD(4,3))" 
Write(io_L,107)4,4,Aimag(ZD(4,4)), "# Aimag(ZD(4,4))" 
Write(io_L,107)4,5,Aimag(ZD(4,5)), "# Aimag(ZD(4,5))" 
Write(io_L,107)4,6,Aimag(ZD(4,6)), "# Aimag(ZD(4,6))" 
Write(io_L,107)5,1,Aimag(ZD(5,1)), "# Aimag(ZD(5,1))" 
Write(io_L,107)5,2,Aimag(ZD(5,2)), "# Aimag(ZD(5,2))" 
Write(io_L,107)5,3,Aimag(ZD(5,3)), "# Aimag(ZD(5,3))" 
Write(io_L,107)5,4,Aimag(ZD(5,4)), "# Aimag(ZD(5,4))" 
Write(io_L,107)5,5,Aimag(ZD(5,5)), "# Aimag(ZD(5,5))" 
Write(io_L,107)5,6,Aimag(ZD(5,6)), "# Aimag(ZD(5,6))" 
Write(io_L,107)6,1,Aimag(ZD(6,1)), "# Aimag(ZD(6,1))" 
Write(io_L,107)6,2,Aimag(ZD(6,2)), "# Aimag(ZD(6,2))" 
Write(io_L,107)6,3,Aimag(ZD(6,3)), "# Aimag(ZD(6,3))" 
Write(io_L,107)6,4,Aimag(ZD(6,4)), "# Aimag(ZD(6,4))" 
Write(io_L,107)6,5,Aimag(ZD(6,5)), "# Aimag(ZD(6,5))" 
Write(io_L,107)6,6,Aimag(ZD(6,6)), "# Aimag(ZD(6,6))" 
End If 

Write(io_L,106) "Block USQMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZU(1,1),dp), "# Real(ZU(1,1),dp)" 
Write(io_L,107)1,2,Real(ZU(1,2),dp), "# Real(ZU(1,2),dp)" 
Write(io_L,107)1,3,Real(ZU(1,3),dp), "# Real(ZU(1,3),dp)" 
Write(io_L,107)1,4,Real(ZU(1,4),dp), "# Real(ZU(1,4),dp)" 
Write(io_L,107)1,5,Real(ZU(1,5),dp), "# Real(ZU(1,5),dp)" 
Write(io_L,107)1,6,Real(ZU(1,6),dp), "# Real(ZU(1,6),dp)" 
Write(io_L,107)2,1,Real(ZU(2,1),dp), "# Real(ZU(2,1),dp)" 
Write(io_L,107)2,2,Real(ZU(2,2),dp), "# Real(ZU(2,2),dp)" 
Write(io_L,107)2,3,Real(ZU(2,3),dp), "# Real(ZU(2,3),dp)" 
Write(io_L,107)2,4,Real(ZU(2,4),dp), "# Real(ZU(2,4),dp)" 
Write(io_L,107)2,5,Real(ZU(2,5),dp), "# Real(ZU(2,5),dp)" 
Write(io_L,107)2,6,Real(ZU(2,6),dp), "# Real(ZU(2,6),dp)" 
Write(io_L,107)3,1,Real(ZU(3,1),dp), "# Real(ZU(3,1),dp)" 
Write(io_L,107)3,2,Real(ZU(3,2),dp), "# Real(ZU(3,2),dp)" 
Write(io_L,107)3,3,Real(ZU(3,3),dp), "# Real(ZU(3,3),dp)" 
Write(io_L,107)3,4,Real(ZU(3,4),dp), "# Real(ZU(3,4),dp)" 
Write(io_L,107)3,5,Real(ZU(3,5),dp), "# Real(ZU(3,5),dp)" 
Write(io_L,107)3,6,Real(ZU(3,6),dp), "# Real(ZU(3,6),dp)" 
Write(io_L,107)4,1,Real(ZU(4,1),dp), "# Real(ZU(4,1),dp)" 
Write(io_L,107)4,2,Real(ZU(4,2),dp), "# Real(ZU(4,2),dp)" 
Write(io_L,107)4,3,Real(ZU(4,3),dp), "# Real(ZU(4,3),dp)" 
Write(io_L,107)4,4,Real(ZU(4,4),dp), "# Real(ZU(4,4),dp)" 
Write(io_L,107)4,5,Real(ZU(4,5),dp), "# Real(ZU(4,5),dp)" 
Write(io_L,107)4,6,Real(ZU(4,6),dp), "# Real(ZU(4,6),dp)" 
Write(io_L,107)5,1,Real(ZU(5,1),dp), "# Real(ZU(5,1),dp)" 
Write(io_L,107)5,2,Real(ZU(5,2),dp), "# Real(ZU(5,2),dp)" 
Write(io_L,107)5,3,Real(ZU(5,3),dp), "# Real(ZU(5,3),dp)" 
Write(io_L,107)5,4,Real(ZU(5,4),dp), "# Real(ZU(5,4),dp)" 
Write(io_L,107)5,5,Real(ZU(5,5),dp), "# Real(ZU(5,5),dp)" 
Write(io_L,107)5,6,Real(ZU(5,6),dp), "# Real(ZU(5,6),dp)" 
Write(io_L,107)6,1,Real(ZU(6,1),dp), "# Real(ZU(6,1),dp)" 
Write(io_L,107)6,2,Real(ZU(6,2),dp), "# Real(ZU(6,2),dp)" 
Write(io_L,107)6,3,Real(ZU(6,3),dp), "# Real(ZU(6,3),dp)" 
Write(io_L,107)6,4,Real(ZU(6,4),dp), "# Real(ZU(6,4),dp)" 
Write(io_L,107)6,5,Real(ZU(6,5),dp), "# Real(ZU(6,5),dp)" 
Write(io_L,107)6,6,Real(ZU(6,6),dp), "# Real(ZU(6,6),dp)" 
If (MaxVal(Abs(AImag(ZU))).gt.0._dp) Then 
Write(io_L,106) "Block IMUSQMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZU(1,1)), "# Aimag(ZU(1,1))" 
Write(io_L,107)1,2,Aimag(ZU(1,2)), "# Aimag(ZU(1,2))" 
Write(io_L,107)1,3,Aimag(ZU(1,3)), "# Aimag(ZU(1,3))" 
Write(io_L,107)1,4,Aimag(ZU(1,4)), "# Aimag(ZU(1,4))" 
Write(io_L,107)1,5,Aimag(ZU(1,5)), "# Aimag(ZU(1,5))" 
Write(io_L,107)1,6,Aimag(ZU(1,6)), "# Aimag(ZU(1,6))" 
Write(io_L,107)2,1,Aimag(ZU(2,1)), "# Aimag(ZU(2,1))" 
Write(io_L,107)2,2,Aimag(ZU(2,2)), "# Aimag(ZU(2,2))" 
Write(io_L,107)2,3,Aimag(ZU(2,3)), "# Aimag(ZU(2,3))" 
Write(io_L,107)2,4,Aimag(ZU(2,4)), "# Aimag(ZU(2,4))" 
Write(io_L,107)2,5,Aimag(ZU(2,5)), "# Aimag(ZU(2,5))" 
Write(io_L,107)2,6,Aimag(ZU(2,6)), "# Aimag(ZU(2,6))" 
Write(io_L,107)3,1,Aimag(ZU(3,1)), "# Aimag(ZU(3,1))" 
Write(io_L,107)3,2,Aimag(ZU(3,2)), "# Aimag(ZU(3,2))" 
Write(io_L,107)3,3,Aimag(ZU(3,3)), "# Aimag(ZU(3,3))" 
Write(io_L,107)3,4,Aimag(ZU(3,4)), "# Aimag(ZU(3,4))" 
Write(io_L,107)3,5,Aimag(ZU(3,5)), "# Aimag(ZU(3,5))" 
Write(io_L,107)3,6,Aimag(ZU(3,6)), "# Aimag(ZU(3,6))" 
Write(io_L,107)4,1,Aimag(ZU(4,1)), "# Aimag(ZU(4,1))" 
Write(io_L,107)4,2,Aimag(ZU(4,2)), "# Aimag(ZU(4,2))" 
Write(io_L,107)4,3,Aimag(ZU(4,3)), "# Aimag(ZU(4,3))" 
Write(io_L,107)4,4,Aimag(ZU(4,4)), "# Aimag(ZU(4,4))" 
Write(io_L,107)4,5,Aimag(ZU(4,5)), "# Aimag(ZU(4,5))" 
Write(io_L,107)4,6,Aimag(ZU(4,6)), "# Aimag(ZU(4,6))" 
Write(io_L,107)5,1,Aimag(ZU(5,1)), "# Aimag(ZU(5,1))" 
Write(io_L,107)5,2,Aimag(ZU(5,2)), "# Aimag(ZU(5,2))" 
Write(io_L,107)5,3,Aimag(ZU(5,3)), "# Aimag(ZU(5,3))" 
Write(io_L,107)5,4,Aimag(ZU(5,4)), "# Aimag(ZU(5,4))" 
Write(io_L,107)5,5,Aimag(ZU(5,5)), "# Aimag(ZU(5,5))" 
Write(io_L,107)5,6,Aimag(ZU(5,6)), "# Aimag(ZU(5,6))" 
Write(io_L,107)6,1,Aimag(ZU(6,1)), "# Aimag(ZU(6,1))" 
Write(io_L,107)6,2,Aimag(ZU(6,2)), "# Aimag(ZU(6,2))" 
Write(io_L,107)6,3,Aimag(ZU(6,3)), "# Aimag(ZU(6,3))" 
Write(io_L,107)6,4,Aimag(ZU(6,4)), "# Aimag(ZU(6,4))" 
Write(io_L,107)6,5,Aimag(ZU(6,5)), "# Aimag(ZU(6,5))" 
Write(io_L,107)6,6,Aimag(ZU(6,6)), "# Aimag(ZU(6,6))" 
End If 

Write(io_L,106) "Block SCALARMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,ZH(1,1), "# ZH(1,1)" 
Write(io_L,107)1,2,ZH(1,2), "# ZH(1,2)" 
Write(io_L,107)1,3,ZH(1,3), "# ZH(1,3)" 
Write(io_L,107)1,4,ZH(1,4), "# ZH(1,4)" 
Write(io_L,107)1,5,ZH(1,5), "# ZH(1,5)" 
Write(io_L,107)1,6,ZH(1,6), "# ZH(1,6)" 
Write(io_L,107)1,7,ZH(1,7), "# ZH(1,7)" 
Write(io_L,107)1,8,ZH(1,8), "# ZH(1,8)" 
Write(io_L,107)2,1,ZH(2,1), "# ZH(2,1)" 
Write(io_L,107)2,2,ZH(2,2), "# ZH(2,2)" 
Write(io_L,107)2,3,ZH(2,3), "# ZH(2,3)" 
Write(io_L,107)2,4,ZH(2,4), "# ZH(2,4)" 
Write(io_L,107)2,5,ZH(2,5), "# ZH(2,5)" 
Write(io_L,107)2,6,ZH(2,6), "# ZH(2,6)" 
Write(io_L,107)2,7,ZH(2,7), "# ZH(2,7)" 
Write(io_L,107)2,8,ZH(2,8), "# ZH(2,8)" 
Write(io_L,107)3,1,ZH(3,1), "# ZH(3,1)" 
Write(io_L,107)3,2,ZH(3,2), "# ZH(3,2)" 
Write(io_L,107)3,3,ZH(3,3), "# ZH(3,3)" 
Write(io_L,107)3,4,ZH(3,4), "# ZH(3,4)" 
Write(io_L,107)3,5,ZH(3,5), "# ZH(3,5)" 
Write(io_L,107)3,6,ZH(3,6), "# ZH(3,6)" 
Write(io_L,107)3,7,ZH(3,7), "# ZH(3,7)" 
Write(io_L,107)3,8,ZH(3,8), "# ZH(3,8)" 
Write(io_L,107)4,1,ZH(4,1), "# ZH(4,1)" 
Write(io_L,107)4,2,ZH(4,2), "# ZH(4,2)" 
Write(io_L,107)4,3,ZH(4,3), "# ZH(4,3)" 
Write(io_L,107)4,4,ZH(4,4), "# ZH(4,4)" 
Write(io_L,107)4,5,ZH(4,5), "# ZH(4,5)" 
Write(io_L,107)4,6,ZH(4,6), "# ZH(4,6)" 
Write(io_L,107)4,7,ZH(4,7), "# ZH(4,7)" 
Write(io_L,107)4,8,ZH(4,8), "# ZH(4,8)" 
Write(io_L,107)5,1,ZH(5,1), "# ZH(5,1)" 
Write(io_L,107)5,2,ZH(5,2), "# ZH(5,2)" 
Write(io_L,107)5,3,ZH(5,3), "# ZH(5,3)" 
Write(io_L,107)5,4,ZH(5,4), "# ZH(5,4)" 
Write(io_L,107)5,5,ZH(5,5), "# ZH(5,5)" 
Write(io_L,107)5,6,ZH(5,6), "# ZH(5,6)" 
Write(io_L,107)5,7,ZH(5,7), "# ZH(5,7)" 
Write(io_L,107)5,8,ZH(5,8), "# ZH(5,8)" 
Write(io_L,107)6,1,ZH(6,1), "# ZH(6,1)" 
Write(io_L,107)6,2,ZH(6,2), "# ZH(6,2)" 
Write(io_L,107)6,3,ZH(6,3), "# ZH(6,3)" 
Write(io_L,107)6,4,ZH(6,4), "# ZH(6,4)" 
Write(io_L,107)6,5,ZH(6,5), "# ZH(6,5)" 
Write(io_L,107)6,6,ZH(6,6), "# ZH(6,6)" 
Write(io_L,107)6,7,ZH(6,7), "# ZH(6,7)" 
Write(io_L,107)6,8,ZH(6,8), "# ZH(6,8)" 
Write(io_L,107)7,1,ZH(7,1), "# ZH(7,1)" 
Write(io_L,107)7,2,ZH(7,2), "# ZH(7,2)" 
Write(io_L,107)7,3,ZH(7,3), "# ZH(7,3)" 
Write(io_L,107)7,4,ZH(7,4), "# ZH(7,4)" 
Write(io_L,107)7,5,ZH(7,5), "# ZH(7,5)" 
Write(io_L,107)7,6,ZH(7,6), "# ZH(7,6)" 
Write(io_L,107)7,7,ZH(7,7), "# ZH(7,7)" 
Write(io_L,107)7,8,ZH(7,8), "# ZH(7,8)" 
Write(io_L,107)8,1,ZH(8,1), "# ZH(8,1)" 
Write(io_L,107)8,2,ZH(8,2), "# ZH(8,2)" 
Write(io_L,107)8,3,ZH(8,3), "# ZH(8,3)" 
Write(io_L,107)8,4,ZH(8,4), "# ZH(8,4)" 
Write(io_L,107)8,5,ZH(8,5), "# ZH(8,5)" 
Write(io_L,107)8,6,ZH(8,6), "# ZH(8,6)" 
Write(io_L,107)8,7,ZH(8,7), "# ZH(8,7)" 
Write(io_L,107)8,8,ZH(8,8), "# ZH(8,8)" 
Write(io_L,106) "Block PSEUDOSCALARMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,ZA(1,1), "# ZA(1,1)" 
Write(io_L,107)1,2,ZA(1,2), "# ZA(1,2)" 
Write(io_L,107)1,3,ZA(1,3), "# ZA(1,3)" 
Write(io_L,107)1,4,ZA(1,4), "# ZA(1,4)" 
Write(io_L,107)1,5,ZA(1,5), "# ZA(1,5)" 
Write(io_L,107)1,6,ZA(1,6), "# ZA(1,6)" 
Write(io_L,107)1,7,ZA(1,7), "# ZA(1,7)" 
Write(io_L,107)1,8,ZA(1,8), "# ZA(1,8)" 
Write(io_L,107)2,1,ZA(2,1), "# ZA(2,1)" 
Write(io_L,107)2,2,ZA(2,2), "# ZA(2,2)" 
Write(io_L,107)2,3,ZA(2,3), "# ZA(2,3)" 
Write(io_L,107)2,4,ZA(2,4), "# ZA(2,4)" 
Write(io_L,107)2,5,ZA(2,5), "# ZA(2,5)" 
Write(io_L,107)2,6,ZA(2,6), "# ZA(2,6)" 
Write(io_L,107)2,7,ZA(2,7), "# ZA(2,7)" 
Write(io_L,107)2,8,ZA(2,8), "# ZA(2,8)" 
Write(io_L,107)3,1,ZA(3,1), "# ZA(3,1)" 
Write(io_L,107)3,2,ZA(3,2), "# ZA(3,2)" 
Write(io_L,107)3,3,ZA(3,3), "# ZA(3,3)" 
Write(io_L,107)3,4,ZA(3,4), "# ZA(3,4)" 
Write(io_L,107)3,5,ZA(3,5), "# ZA(3,5)" 
Write(io_L,107)3,6,ZA(3,6), "# ZA(3,6)" 
Write(io_L,107)3,7,ZA(3,7), "# ZA(3,7)" 
Write(io_L,107)3,8,ZA(3,8), "# ZA(3,8)" 
Write(io_L,107)4,1,ZA(4,1), "# ZA(4,1)" 
Write(io_L,107)4,2,ZA(4,2), "# ZA(4,2)" 
Write(io_L,107)4,3,ZA(4,3), "# ZA(4,3)" 
Write(io_L,107)4,4,ZA(4,4), "# ZA(4,4)" 
Write(io_L,107)4,5,ZA(4,5), "# ZA(4,5)" 
Write(io_L,107)4,6,ZA(4,6), "# ZA(4,6)" 
Write(io_L,107)4,7,ZA(4,7), "# ZA(4,7)" 
Write(io_L,107)4,8,ZA(4,8), "# ZA(4,8)" 
Write(io_L,107)5,1,ZA(5,1), "# ZA(5,1)" 
Write(io_L,107)5,2,ZA(5,2), "# ZA(5,2)" 
Write(io_L,107)5,3,ZA(5,3), "# ZA(5,3)" 
Write(io_L,107)5,4,ZA(5,4), "# ZA(5,4)" 
Write(io_L,107)5,5,ZA(5,5), "# ZA(5,5)" 
Write(io_L,107)5,6,ZA(5,6), "# ZA(5,6)" 
Write(io_L,107)5,7,ZA(5,7), "# ZA(5,7)" 
Write(io_L,107)5,8,ZA(5,8), "# ZA(5,8)" 
Write(io_L,107)6,1,ZA(6,1), "# ZA(6,1)" 
Write(io_L,107)6,2,ZA(6,2), "# ZA(6,2)" 
Write(io_L,107)6,3,ZA(6,3), "# ZA(6,3)" 
Write(io_L,107)6,4,ZA(6,4), "# ZA(6,4)" 
Write(io_L,107)6,5,ZA(6,5), "# ZA(6,5)" 
Write(io_L,107)6,6,ZA(6,6), "# ZA(6,6)" 
Write(io_L,107)6,7,ZA(6,7), "# ZA(6,7)" 
Write(io_L,107)6,8,ZA(6,8), "# ZA(6,8)" 
Write(io_L,107)7,1,ZA(7,1), "# ZA(7,1)" 
Write(io_L,107)7,2,ZA(7,2), "# ZA(7,2)" 
Write(io_L,107)7,3,ZA(7,3), "# ZA(7,3)" 
Write(io_L,107)7,4,ZA(7,4), "# ZA(7,4)" 
Write(io_L,107)7,5,ZA(7,5), "# ZA(7,5)" 
Write(io_L,107)7,6,ZA(7,6), "# ZA(7,6)" 
Write(io_L,107)7,7,ZA(7,7), "# ZA(7,7)" 
Write(io_L,107)7,8,ZA(7,8), "# ZA(7,8)" 
Write(io_L,107)8,1,ZA(8,1), "# ZA(8,1)" 
Write(io_L,107)8,2,ZA(8,2), "# ZA(8,2)" 
Write(io_L,107)8,3,ZA(8,3), "# ZA(8,3)" 
Write(io_L,107)8,4,ZA(8,4), "# ZA(8,4)" 
Write(io_L,107)8,5,ZA(8,5), "# ZA(8,5)" 
Write(io_L,107)8,6,ZA(8,6), "# ZA(8,6)" 
Write(io_L,107)8,7,ZA(8,7), "# ZA(8,7)" 
Write(io_L,107)8,8,ZA(8,8), "# ZA(8,8)" 
Write(io_L,106) "Block CHARGEMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,ZP(1,1), "# ZP(1,1)" 
Write(io_L,107)1,2,ZP(1,2), "# ZP(1,2)" 
Write(io_L,107)1,3,ZP(1,3), "# ZP(1,3)" 
Write(io_L,107)1,4,ZP(1,4), "# ZP(1,4)" 
Write(io_L,107)1,5,ZP(1,5), "# ZP(1,5)" 
Write(io_L,107)1,6,ZP(1,6), "# ZP(1,6)" 
Write(io_L,107)1,7,ZP(1,7), "# ZP(1,7)" 
Write(io_L,107)1,8,ZP(1,8), "# ZP(1,8)" 
Write(io_L,107)2,1,ZP(2,1), "# ZP(2,1)" 
Write(io_L,107)2,2,ZP(2,2), "# ZP(2,2)" 
Write(io_L,107)2,3,ZP(2,3), "# ZP(2,3)" 
Write(io_L,107)2,4,ZP(2,4), "# ZP(2,4)" 
Write(io_L,107)2,5,ZP(2,5), "# ZP(2,5)" 
Write(io_L,107)2,6,ZP(2,6), "# ZP(2,6)" 
Write(io_L,107)2,7,ZP(2,7), "# ZP(2,7)" 
Write(io_L,107)2,8,ZP(2,8), "# ZP(2,8)" 
Write(io_L,107)3,1,ZP(3,1), "# ZP(3,1)" 
Write(io_L,107)3,2,ZP(3,2), "# ZP(3,2)" 
Write(io_L,107)3,3,ZP(3,3), "# ZP(3,3)" 
Write(io_L,107)3,4,ZP(3,4), "# ZP(3,4)" 
Write(io_L,107)3,5,ZP(3,5), "# ZP(3,5)" 
Write(io_L,107)3,6,ZP(3,6), "# ZP(3,6)" 
Write(io_L,107)3,7,ZP(3,7), "# ZP(3,7)" 
Write(io_L,107)3,8,ZP(3,8), "# ZP(3,8)" 
Write(io_L,107)4,1,ZP(4,1), "# ZP(4,1)" 
Write(io_L,107)4,2,ZP(4,2), "# ZP(4,2)" 
Write(io_L,107)4,3,ZP(4,3), "# ZP(4,3)" 
Write(io_L,107)4,4,ZP(4,4), "# ZP(4,4)" 
Write(io_L,107)4,5,ZP(4,5), "# ZP(4,5)" 
Write(io_L,107)4,6,ZP(4,6), "# ZP(4,6)" 
Write(io_L,107)4,7,ZP(4,7), "# ZP(4,7)" 
Write(io_L,107)4,8,ZP(4,8), "# ZP(4,8)" 
Write(io_L,107)5,1,ZP(5,1), "# ZP(5,1)" 
Write(io_L,107)5,2,ZP(5,2), "# ZP(5,2)" 
Write(io_L,107)5,3,ZP(5,3), "# ZP(5,3)" 
Write(io_L,107)5,4,ZP(5,4), "# ZP(5,4)" 
Write(io_L,107)5,5,ZP(5,5), "# ZP(5,5)" 
Write(io_L,107)5,6,ZP(5,6), "# ZP(5,6)" 
Write(io_L,107)5,7,ZP(5,7), "# ZP(5,7)" 
Write(io_L,107)5,8,ZP(5,8), "# ZP(5,8)" 
Write(io_L,107)6,1,ZP(6,1), "# ZP(6,1)" 
Write(io_L,107)6,2,ZP(6,2), "# ZP(6,2)" 
Write(io_L,107)6,3,ZP(6,3), "# ZP(6,3)" 
Write(io_L,107)6,4,ZP(6,4), "# ZP(6,4)" 
Write(io_L,107)6,5,ZP(6,5), "# ZP(6,5)" 
Write(io_L,107)6,6,ZP(6,6), "# ZP(6,6)" 
Write(io_L,107)6,7,ZP(6,7), "# ZP(6,7)" 
Write(io_L,107)6,8,ZP(6,8), "# ZP(6,8)" 
Write(io_L,107)7,1,ZP(7,1), "# ZP(7,1)" 
Write(io_L,107)7,2,ZP(7,2), "# ZP(7,2)" 
Write(io_L,107)7,3,ZP(7,3), "# ZP(7,3)" 
Write(io_L,107)7,4,ZP(7,4), "# ZP(7,4)" 
Write(io_L,107)7,5,ZP(7,5), "# ZP(7,5)" 
Write(io_L,107)7,6,ZP(7,6), "# ZP(7,6)" 
Write(io_L,107)7,7,ZP(7,7), "# ZP(7,7)" 
Write(io_L,107)7,8,ZP(7,8), "# ZP(7,8)" 
Write(io_L,107)8,1,ZP(8,1), "# ZP(8,1)" 
Write(io_L,107)8,2,ZP(8,2), "# ZP(8,2)" 
Write(io_L,107)8,3,ZP(8,3), "# ZP(8,3)" 
Write(io_L,107)8,4,ZP(8,4), "# ZP(8,4)" 
Write(io_L,107)8,5,ZP(8,5), "# ZP(8,5)" 
Write(io_L,107)8,6,ZP(8,6), "# ZP(8,6)" 
Write(io_L,107)8,7,ZP(8,7), "# ZP(8,7)" 
Write(io_L,107)8,8,ZP(8,8), "# ZP(8,8)" 
Write(io_L,106) "Block UVMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(UV(1,1),dp), "# Real(UV(1,1),dp)" 
Write(io_L,107)1,2,Real(UV(1,2),dp), "# Real(UV(1,2),dp)" 
Write(io_L,107)1,3,Real(UV(1,3),dp), "# Real(UV(1,3),dp)" 
Write(io_L,107)1,4,Real(UV(1,4),dp), "# Real(UV(1,4),dp)" 
Write(io_L,107)1,5,Real(UV(1,5),dp), "# Real(UV(1,5),dp)" 
Write(io_L,107)1,6,Real(UV(1,6),dp), "# Real(UV(1,6),dp)" 
Write(io_L,107)1,7,Real(UV(1,7),dp), "# Real(UV(1,7),dp)" 
Write(io_L,107)1,8,Real(UV(1,8),dp), "# Real(UV(1,8),dp)" 
Write(io_L,107)1,9,Real(UV(1,9),dp), "# Real(UV(1,9),dp)" 
Write(io_L,107)1,10,Real(UV(1,10),dp), "# Real(UV(1,10),dp)" 
Write(io_L,107)2,1,Real(UV(2,1),dp), "# Real(UV(2,1),dp)" 
Write(io_L,107)2,2,Real(UV(2,2),dp), "# Real(UV(2,2),dp)" 
Write(io_L,107)2,3,Real(UV(2,3),dp), "# Real(UV(2,3),dp)" 
Write(io_L,107)2,4,Real(UV(2,4),dp), "# Real(UV(2,4),dp)" 
Write(io_L,107)2,5,Real(UV(2,5),dp), "# Real(UV(2,5),dp)" 
Write(io_L,107)2,6,Real(UV(2,6),dp), "# Real(UV(2,6),dp)" 
Write(io_L,107)2,7,Real(UV(2,7),dp), "# Real(UV(2,7),dp)" 
Write(io_L,107)2,8,Real(UV(2,8),dp), "# Real(UV(2,8),dp)" 
Write(io_L,107)2,9,Real(UV(2,9),dp), "# Real(UV(2,9),dp)" 
Write(io_L,107)2,10,Real(UV(2,10),dp), "# Real(UV(2,10),dp)" 
Write(io_L,107)3,1,Real(UV(3,1),dp), "# Real(UV(3,1),dp)" 
Write(io_L,107)3,2,Real(UV(3,2),dp), "# Real(UV(3,2),dp)" 
Write(io_L,107)3,3,Real(UV(3,3),dp), "# Real(UV(3,3),dp)" 
Write(io_L,107)3,4,Real(UV(3,4),dp), "# Real(UV(3,4),dp)" 
Write(io_L,107)3,5,Real(UV(3,5),dp), "# Real(UV(3,5),dp)" 
Write(io_L,107)3,6,Real(UV(3,6),dp), "# Real(UV(3,6),dp)" 
Write(io_L,107)3,7,Real(UV(3,7),dp), "# Real(UV(3,7),dp)" 
Write(io_L,107)3,8,Real(UV(3,8),dp), "# Real(UV(3,8),dp)" 
Write(io_L,107)3,9,Real(UV(3,9),dp), "# Real(UV(3,9),dp)" 
Write(io_L,107)3,10,Real(UV(3,10),dp), "# Real(UV(3,10),dp)" 
Write(io_L,107)4,1,Real(UV(4,1),dp), "# Real(UV(4,1),dp)" 
Write(io_L,107)4,2,Real(UV(4,2),dp), "# Real(UV(4,2),dp)" 
Write(io_L,107)4,3,Real(UV(4,3),dp), "# Real(UV(4,3),dp)" 
Write(io_L,107)4,4,Real(UV(4,4),dp), "# Real(UV(4,4),dp)" 
Write(io_L,107)4,5,Real(UV(4,5),dp), "# Real(UV(4,5),dp)" 
Write(io_L,107)4,6,Real(UV(4,6),dp), "# Real(UV(4,6),dp)" 
Write(io_L,107)4,7,Real(UV(4,7),dp), "# Real(UV(4,7),dp)" 
Write(io_L,107)4,8,Real(UV(4,8),dp), "# Real(UV(4,8),dp)" 
Write(io_L,107)4,9,Real(UV(4,9),dp), "# Real(UV(4,9),dp)" 
Write(io_L,107)4,10,Real(UV(4,10),dp), "# Real(UV(4,10),dp)" 
Write(io_L,107)5,1,Real(UV(5,1),dp), "# Real(UV(5,1),dp)" 
Write(io_L,107)5,2,Real(UV(5,2),dp), "# Real(UV(5,2),dp)" 
Write(io_L,107)5,3,Real(UV(5,3),dp), "# Real(UV(5,3),dp)" 
Write(io_L,107)5,4,Real(UV(5,4),dp), "# Real(UV(5,4),dp)" 
Write(io_L,107)5,5,Real(UV(5,5),dp), "# Real(UV(5,5),dp)" 
Write(io_L,107)5,6,Real(UV(5,6),dp), "# Real(UV(5,6),dp)" 
Write(io_L,107)5,7,Real(UV(5,7),dp), "# Real(UV(5,7),dp)" 
Write(io_L,107)5,8,Real(UV(5,8),dp), "# Real(UV(5,8),dp)" 
Write(io_L,107)5,9,Real(UV(5,9),dp), "# Real(UV(5,9),dp)" 
Write(io_L,107)5,10,Real(UV(5,10),dp), "# Real(UV(5,10),dp)" 
Write(io_L,107)6,1,Real(UV(6,1),dp), "# Real(UV(6,1),dp)" 
Write(io_L,107)6,2,Real(UV(6,2),dp), "# Real(UV(6,2),dp)" 
Write(io_L,107)6,3,Real(UV(6,3),dp), "# Real(UV(6,3),dp)" 
Write(io_L,107)6,4,Real(UV(6,4),dp), "# Real(UV(6,4),dp)" 
Write(io_L,107)6,5,Real(UV(6,5),dp), "# Real(UV(6,5),dp)" 
Write(io_L,107)6,6,Real(UV(6,6),dp), "# Real(UV(6,6),dp)" 
Write(io_L,107)6,7,Real(UV(6,7),dp), "# Real(UV(6,7),dp)" 
Write(io_L,107)6,8,Real(UV(6,8),dp), "# Real(UV(6,8),dp)" 
Write(io_L,107)6,9,Real(UV(6,9),dp), "# Real(UV(6,9),dp)" 
Write(io_L,107)6,10,Real(UV(6,10),dp), "# Real(UV(6,10),dp)" 
Write(io_L,107)7,1,Real(UV(7,1),dp), "# Real(UV(7,1),dp)" 
Write(io_L,107)7,2,Real(UV(7,2),dp), "# Real(UV(7,2),dp)" 
Write(io_L,107)7,3,Real(UV(7,3),dp), "# Real(UV(7,3),dp)" 
Write(io_L,107)7,4,Real(UV(7,4),dp), "# Real(UV(7,4),dp)" 
Write(io_L,107)7,5,Real(UV(7,5),dp), "# Real(UV(7,5),dp)" 
Write(io_L,107)7,6,Real(UV(7,6),dp), "# Real(UV(7,6),dp)" 
Write(io_L,107)7,7,Real(UV(7,7),dp), "# Real(UV(7,7),dp)" 
Write(io_L,107)7,8,Real(UV(7,8),dp), "# Real(UV(7,8),dp)" 
Write(io_L,107)7,9,Real(UV(7,9),dp), "# Real(UV(7,9),dp)" 
Write(io_L,107)7,10,Real(UV(7,10),dp), "# Real(UV(7,10),dp)" 
Write(io_L,107)8,1,Real(UV(8,1),dp), "# Real(UV(8,1),dp)" 
Write(io_L,107)8,2,Real(UV(8,2),dp), "# Real(UV(8,2),dp)" 
Write(io_L,107)8,3,Real(UV(8,3),dp), "# Real(UV(8,3),dp)" 
Write(io_L,107)8,4,Real(UV(8,4),dp), "# Real(UV(8,4),dp)" 
Write(io_L,107)8,5,Real(UV(8,5),dp), "# Real(UV(8,5),dp)" 
Write(io_L,107)8,6,Real(UV(8,6),dp), "# Real(UV(8,6),dp)" 
Write(io_L,107)8,7,Real(UV(8,7),dp), "# Real(UV(8,7),dp)" 
Write(io_L,107)8,8,Real(UV(8,8),dp), "# Real(UV(8,8),dp)" 
Write(io_L,107)8,9,Real(UV(8,9),dp), "# Real(UV(8,9),dp)" 
Write(io_L,107)8,10,Real(UV(8,10),dp), "# Real(UV(8,10),dp)" 
Write(io_L,107)9,1,Real(UV(9,1),dp), "# Real(UV(9,1),dp)" 
Write(io_L,107)9,2,Real(UV(9,2),dp), "# Real(UV(9,2),dp)" 
Write(io_L,107)9,3,Real(UV(9,3),dp), "# Real(UV(9,3),dp)" 
Write(io_L,107)9,4,Real(UV(9,4),dp), "# Real(UV(9,4),dp)" 
Write(io_L,107)9,5,Real(UV(9,5),dp), "# Real(UV(9,5),dp)" 
Write(io_L,107)9,6,Real(UV(9,6),dp), "# Real(UV(9,6),dp)" 
Write(io_L,107)9,7,Real(UV(9,7),dp), "# Real(UV(9,7),dp)" 
Write(io_L,107)9,8,Real(UV(9,8),dp), "# Real(UV(9,8),dp)" 
Write(io_L,107)9,9,Real(UV(9,9),dp), "# Real(UV(9,9),dp)" 
Write(io_L,107)9,10,Real(UV(9,10),dp), "# Real(UV(9,10),dp)" 
Write(io_L,107)10,1,Real(UV(10,1),dp), "# Real(UV(10,1),dp)" 
Write(io_L,107)10,2,Real(UV(10,2),dp), "# Real(UV(10,2),dp)" 
Write(io_L,107)10,3,Real(UV(10,3),dp), "# Real(UV(10,3),dp)" 
Write(io_L,107)10,4,Real(UV(10,4),dp), "# Real(UV(10,4),dp)" 
Write(io_L,107)10,5,Real(UV(10,5),dp), "# Real(UV(10,5),dp)" 
Write(io_L,107)10,6,Real(UV(10,6),dp), "# Real(UV(10,6),dp)" 
Write(io_L,107)10,7,Real(UV(10,7),dp), "# Real(UV(10,7),dp)" 
Write(io_L,107)10,8,Real(UV(10,8),dp), "# Real(UV(10,8),dp)" 
Write(io_L,107)10,9,Real(UV(10,9),dp), "# Real(UV(10,9),dp)" 
Write(io_L,107)10,10,Real(UV(10,10),dp), "# Real(UV(10,10),dp)" 
If (MaxVal(Abs(AImag(UV))).gt.0._dp) Then 
Write(io_L,106) "Block IMUVMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(UV(1,1)), "# Aimag(UV(1,1))" 
Write(io_L,107)1,2,Aimag(UV(1,2)), "# Aimag(UV(1,2))" 
Write(io_L,107)1,3,Aimag(UV(1,3)), "# Aimag(UV(1,3))" 
Write(io_L,107)1,4,Aimag(UV(1,4)), "# Aimag(UV(1,4))" 
Write(io_L,107)1,5,Aimag(UV(1,5)), "# Aimag(UV(1,5))" 
Write(io_L,107)1,6,Aimag(UV(1,6)), "# Aimag(UV(1,6))" 
Write(io_L,107)1,7,Aimag(UV(1,7)), "# Aimag(UV(1,7))" 
Write(io_L,107)1,8,Aimag(UV(1,8)), "# Aimag(UV(1,8))" 
Write(io_L,107)1,9,Aimag(UV(1,9)), "# Aimag(UV(1,9))" 
Write(io_L,107)1,10,Aimag(UV(1,10)), "# Aimag(UV(1,10))" 
Write(io_L,107)2,1,Aimag(UV(2,1)), "# Aimag(UV(2,1))" 
Write(io_L,107)2,2,Aimag(UV(2,2)), "# Aimag(UV(2,2))" 
Write(io_L,107)2,3,Aimag(UV(2,3)), "# Aimag(UV(2,3))" 
Write(io_L,107)2,4,Aimag(UV(2,4)), "# Aimag(UV(2,4))" 
Write(io_L,107)2,5,Aimag(UV(2,5)), "# Aimag(UV(2,5))" 
Write(io_L,107)2,6,Aimag(UV(2,6)), "# Aimag(UV(2,6))" 
Write(io_L,107)2,7,Aimag(UV(2,7)), "# Aimag(UV(2,7))" 
Write(io_L,107)2,8,Aimag(UV(2,8)), "# Aimag(UV(2,8))" 
Write(io_L,107)2,9,Aimag(UV(2,9)), "# Aimag(UV(2,9))" 
Write(io_L,107)2,10,Aimag(UV(2,10)), "# Aimag(UV(2,10))" 
Write(io_L,107)3,1,Aimag(UV(3,1)), "# Aimag(UV(3,1))" 
Write(io_L,107)3,2,Aimag(UV(3,2)), "# Aimag(UV(3,2))" 
Write(io_L,107)3,3,Aimag(UV(3,3)), "# Aimag(UV(3,3))" 
Write(io_L,107)3,4,Aimag(UV(3,4)), "# Aimag(UV(3,4))" 
Write(io_L,107)3,5,Aimag(UV(3,5)), "# Aimag(UV(3,5))" 
Write(io_L,107)3,6,Aimag(UV(3,6)), "# Aimag(UV(3,6))" 
Write(io_L,107)3,7,Aimag(UV(3,7)), "# Aimag(UV(3,7))" 
Write(io_L,107)3,8,Aimag(UV(3,8)), "# Aimag(UV(3,8))" 
Write(io_L,107)3,9,Aimag(UV(3,9)), "# Aimag(UV(3,9))" 
Write(io_L,107)3,10,Aimag(UV(3,10)), "# Aimag(UV(3,10))" 
Write(io_L,107)4,1,Aimag(UV(4,1)), "# Aimag(UV(4,1))" 
Write(io_L,107)4,2,Aimag(UV(4,2)), "# Aimag(UV(4,2))" 
Write(io_L,107)4,3,Aimag(UV(4,3)), "# Aimag(UV(4,3))" 
Write(io_L,107)4,4,Aimag(UV(4,4)), "# Aimag(UV(4,4))" 
Write(io_L,107)4,5,Aimag(UV(4,5)), "# Aimag(UV(4,5))" 
Write(io_L,107)4,6,Aimag(UV(4,6)), "# Aimag(UV(4,6))" 
Write(io_L,107)4,7,Aimag(UV(4,7)), "# Aimag(UV(4,7))" 
Write(io_L,107)4,8,Aimag(UV(4,8)), "# Aimag(UV(4,8))" 
Write(io_L,107)4,9,Aimag(UV(4,9)), "# Aimag(UV(4,9))" 
Write(io_L,107)4,10,Aimag(UV(4,10)), "# Aimag(UV(4,10))" 
Write(io_L,107)5,1,Aimag(UV(5,1)), "# Aimag(UV(5,1))" 
Write(io_L,107)5,2,Aimag(UV(5,2)), "# Aimag(UV(5,2))" 
Write(io_L,107)5,3,Aimag(UV(5,3)), "# Aimag(UV(5,3))" 
Write(io_L,107)5,4,Aimag(UV(5,4)), "# Aimag(UV(5,4))" 
Write(io_L,107)5,5,Aimag(UV(5,5)), "# Aimag(UV(5,5))" 
Write(io_L,107)5,6,Aimag(UV(5,6)), "# Aimag(UV(5,6))" 
Write(io_L,107)5,7,Aimag(UV(5,7)), "# Aimag(UV(5,7))" 
Write(io_L,107)5,8,Aimag(UV(5,8)), "# Aimag(UV(5,8))" 
Write(io_L,107)5,9,Aimag(UV(5,9)), "# Aimag(UV(5,9))" 
Write(io_L,107)5,10,Aimag(UV(5,10)), "# Aimag(UV(5,10))" 
Write(io_L,107)6,1,Aimag(UV(6,1)), "# Aimag(UV(6,1))" 
Write(io_L,107)6,2,Aimag(UV(6,2)), "# Aimag(UV(6,2))" 
Write(io_L,107)6,3,Aimag(UV(6,3)), "# Aimag(UV(6,3))" 
Write(io_L,107)6,4,Aimag(UV(6,4)), "# Aimag(UV(6,4))" 
Write(io_L,107)6,5,Aimag(UV(6,5)), "# Aimag(UV(6,5))" 
Write(io_L,107)6,6,Aimag(UV(6,6)), "# Aimag(UV(6,6))" 
Write(io_L,107)6,7,Aimag(UV(6,7)), "# Aimag(UV(6,7))" 
Write(io_L,107)6,8,Aimag(UV(6,8)), "# Aimag(UV(6,8))" 
Write(io_L,107)6,9,Aimag(UV(6,9)), "# Aimag(UV(6,9))" 
Write(io_L,107)6,10,Aimag(UV(6,10)), "# Aimag(UV(6,10))" 
Write(io_L,107)7,1,Aimag(UV(7,1)), "# Aimag(UV(7,1))" 
Write(io_L,107)7,2,Aimag(UV(7,2)), "# Aimag(UV(7,2))" 
Write(io_L,107)7,3,Aimag(UV(7,3)), "# Aimag(UV(7,3))" 
Write(io_L,107)7,4,Aimag(UV(7,4)), "# Aimag(UV(7,4))" 
Write(io_L,107)7,5,Aimag(UV(7,5)), "# Aimag(UV(7,5))" 
Write(io_L,107)7,6,Aimag(UV(7,6)), "# Aimag(UV(7,6))" 
Write(io_L,107)7,7,Aimag(UV(7,7)), "# Aimag(UV(7,7))" 
Write(io_L,107)7,8,Aimag(UV(7,8)), "# Aimag(UV(7,8))" 
Write(io_L,107)7,9,Aimag(UV(7,9)), "# Aimag(UV(7,9))" 
Write(io_L,107)7,10,Aimag(UV(7,10)), "# Aimag(UV(7,10))" 
Write(io_L,107)8,1,Aimag(UV(8,1)), "# Aimag(UV(8,1))" 
Write(io_L,107)8,2,Aimag(UV(8,2)), "# Aimag(UV(8,2))" 
Write(io_L,107)8,3,Aimag(UV(8,3)), "# Aimag(UV(8,3))" 
Write(io_L,107)8,4,Aimag(UV(8,4)), "# Aimag(UV(8,4))" 
Write(io_L,107)8,5,Aimag(UV(8,5)), "# Aimag(UV(8,5))" 
Write(io_L,107)8,6,Aimag(UV(8,6)), "# Aimag(UV(8,6))" 
Write(io_L,107)8,7,Aimag(UV(8,7)), "# Aimag(UV(8,7))" 
Write(io_L,107)8,8,Aimag(UV(8,8)), "# Aimag(UV(8,8))" 
Write(io_L,107)8,9,Aimag(UV(8,9)), "# Aimag(UV(8,9))" 
Write(io_L,107)8,10,Aimag(UV(8,10)), "# Aimag(UV(8,10))" 
Write(io_L,107)9,1,Aimag(UV(9,1)), "# Aimag(UV(9,1))" 
Write(io_L,107)9,2,Aimag(UV(9,2)), "# Aimag(UV(9,2))" 
Write(io_L,107)9,3,Aimag(UV(9,3)), "# Aimag(UV(9,3))" 
Write(io_L,107)9,4,Aimag(UV(9,4)), "# Aimag(UV(9,4))" 
Write(io_L,107)9,5,Aimag(UV(9,5)), "# Aimag(UV(9,5))" 
Write(io_L,107)9,6,Aimag(UV(9,6)), "# Aimag(UV(9,6))" 
Write(io_L,107)9,7,Aimag(UV(9,7)), "# Aimag(UV(9,7))" 
Write(io_L,107)9,8,Aimag(UV(9,8)), "# Aimag(UV(9,8))" 
Write(io_L,107)9,9,Aimag(UV(9,9)), "# Aimag(UV(9,9))" 
Write(io_L,107)9,10,Aimag(UV(9,10)), "# Aimag(UV(9,10))" 
Write(io_L,107)10,1,Aimag(UV(10,1)), "# Aimag(UV(10,1))" 
Write(io_L,107)10,2,Aimag(UV(10,2)), "# Aimag(UV(10,2))" 
Write(io_L,107)10,3,Aimag(UV(10,3)), "# Aimag(UV(10,3))" 
Write(io_L,107)10,4,Aimag(UV(10,4)), "# Aimag(UV(10,4))" 
Write(io_L,107)10,5,Aimag(UV(10,5)), "# Aimag(UV(10,5))" 
Write(io_L,107)10,6,Aimag(UV(10,6)), "# Aimag(UV(10,6))" 
Write(io_L,107)10,7,Aimag(UV(10,7)), "# Aimag(UV(10,7))" 
Write(io_L,107)10,8,Aimag(UV(10,8)), "# Aimag(UV(10,8))" 
Write(io_L,107)10,9,Aimag(UV(10,9)), "# Aimag(UV(10,9))" 
Write(io_L,107)10,10,Aimag(UV(10,10)), "# Aimag(UV(10,10))" 
End If 

Write(io_L,106) "Block UERMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZER(1,1),dp), "# Real(ZER(1,1),dp)" 
Write(io_L,107)1,2,Real(ZER(1,2),dp), "# Real(ZER(1,2),dp)" 
Write(io_L,107)1,3,Real(ZER(1,3),dp), "# Real(ZER(1,3),dp)" 
Write(io_L,107)1,4,Real(ZER(1,4),dp), "# Real(ZER(1,4),dp)" 
Write(io_L,107)1,5,Real(ZER(1,5),dp), "# Real(ZER(1,5),dp)" 
Write(io_L,107)2,1,Real(ZER(2,1),dp), "# Real(ZER(2,1),dp)" 
Write(io_L,107)2,2,Real(ZER(2,2),dp), "# Real(ZER(2,2),dp)" 
Write(io_L,107)2,3,Real(ZER(2,3),dp), "# Real(ZER(2,3),dp)" 
Write(io_L,107)2,4,Real(ZER(2,4),dp), "# Real(ZER(2,4),dp)" 
Write(io_L,107)2,5,Real(ZER(2,5),dp), "# Real(ZER(2,5),dp)" 
Write(io_L,107)3,1,Real(ZER(3,1),dp), "# Real(ZER(3,1),dp)" 
Write(io_L,107)3,2,Real(ZER(3,2),dp), "# Real(ZER(3,2),dp)" 
Write(io_L,107)3,3,Real(ZER(3,3),dp), "# Real(ZER(3,3),dp)" 
Write(io_L,107)3,4,Real(ZER(3,4),dp), "# Real(ZER(3,4),dp)" 
Write(io_L,107)3,5,Real(ZER(3,5),dp), "# Real(ZER(3,5),dp)" 
Write(io_L,107)4,1,Real(ZER(4,1),dp), "# Real(ZER(4,1),dp)" 
Write(io_L,107)4,2,Real(ZER(4,2),dp), "# Real(ZER(4,2),dp)" 
Write(io_L,107)4,3,Real(ZER(4,3),dp), "# Real(ZER(4,3),dp)" 
Write(io_L,107)4,4,Real(ZER(4,4),dp), "# Real(ZER(4,4),dp)" 
Write(io_L,107)4,5,Real(ZER(4,5),dp), "# Real(ZER(4,5),dp)" 
Write(io_L,107)5,1,Real(ZER(5,1),dp), "# Real(ZER(5,1),dp)" 
Write(io_L,107)5,2,Real(ZER(5,2),dp), "# Real(ZER(5,2),dp)" 
Write(io_L,107)5,3,Real(ZER(5,3),dp), "# Real(ZER(5,3),dp)" 
Write(io_L,107)5,4,Real(ZER(5,4),dp), "# Real(ZER(5,4),dp)" 
Write(io_L,107)5,5,Real(ZER(5,5),dp), "# Real(ZER(5,5),dp)" 
If (MaxVal(Abs(AImag(ZER))).gt.0._dp) Then 
Write(io_L,106) "Block IMUERMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZER(1,1)), "# Aimag(ZER(1,1))" 
Write(io_L,107)1,2,Aimag(ZER(1,2)), "# Aimag(ZER(1,2))" 
Write(io_L,107)1,3,Aimag(ZER(1,3)), "# Aimag(ZER(1,3))" 
Write(io_L,107)1,4,Aimag(ZER(1,4)), "# Aimag(ZER(1,4))" 
Write(io_L,107)1,5,Aimag(ZER(1,5)), "# Aimag(ZER(1,5))" 
Write(io_L,107)2,1,Aimag(ZER(2,1)), "# Aimag(ZER(2,1))" 
Write(io_L,107)2,2,Aimag(ZER(2,2)), "# Aimag(ZER(2,2))" 
Write(io_L,107)2,3,Aimag(ZER(2,3)), "# Aimag(ZER(2,3))" 
Write(io_L,107)2,4,Aimag(ZER(2,4)), "# Aimag(ZER(2,4))" 
Write(io_L,107)2,5,Aimag(ZER(2,5)), "# Aimag(ZER(2,5))" 
Write(io_L,107)3,1,Aimag(ZER(3,1)), "# Aimag(ZER(3,1))" 
Write(io_L,107)3,2,Aimag(ZER(3,2)), "# Aimag(ZER(3,2))" 
Write(io_L,107)3,3,Aimag(ZER(3,3)), "# Aimag(ZER(3,3))" 
Write(io_L,107)3,4,Aimag(ZER(3,4)), "# Aimag(ZER(3,4))" 
Write(io_L,107)3,5,Aimag(ZER(3,5)), "# Aimag(ZER(3,5))" 
Write(io_L,107)4,1,Aimag(ZER(4,1)), "# Aimag(ZER(4,1))" 
Write(io_L,107)4,2,Aimag(ZER(4,2)), "# Aimag(ZER(4,2))" 
Write(io_L,107)4,3,Aimag(ZER(4,3)), "# Aimag(ZER(4,3))" 
Write(io_L,107)4,4,Aimag(ZER(4,4)), "# Aimag(ZER(4,4))" 
Write(io_L,107)4,5,Aimag(ZER(4,5)), "# Aimag(ZER(4,5))" 
Write(io_L,107)5,1,Aimag(ZER(5,1)), "# Aimag(ZER(5,1))" 
Write(io_L,107)5,2,Aimag(ZER(5,2)), "# Aimag(ZER(5,2))" 
Write(io_L,107)5,3,Aimag(ZER(5,3)), "# Aimag(ZER(5,3))" 
Write(io_L,107)5,4,Aimag(ZER(5,4)), "# Aimag(ZER(5,4))" 
Write(io_L,107)5,5,Aimag(ZER(5,5)), "# Aimag(ZER(5,5))" 
End If 

Write(io_L,106) "Block UELMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZEL(1,1),dp), "# Real(ZEL(1,1),dp)" 
Write(io_L,107)1,2,Real(ZEL(1,2),dp), "# Real(ZEL(1,2),dp)" 
Write(io_L,107)1,3,Real(ZEL(1,3),dp), "# Real(ZEL(1,3),dp)" 
Write(io_L,107)1,4,Real(ZEL(1,4),dp), "# Real(ZEL(1,4),dp)" 
Write(io_L,107)1,5,Real(ZEL(1,5),dp), "# Real(ZEL(1,5),dp)" 
Write(io_L,107)2,1,Real(ZEL(2,1),dp), "# Real(ZEL(2,1),dp)" 
Write(io_L,107)2,2,Real(ZEL(2,2),dp), "# Real(ZEL(2,2),dp)" 
Write(io_L,107)2,3,Real(ZEL(2,3),dp), "# Real(ZEL(2,3),dp)" 
Write(io_L,107)2,4,Real(ZEL(2,4),dp), "# Real(ZEL(2,4),dp)" 
Write(io_L,107)2,5,Real(ZEL(2,5),dp), "# Real(ZEL(2,5),dp)" 
Write(io_L,107)3,1,Real(ZEL(3,1),dp), "# Real(ZEL(3,1),dp)" 
Write(io_L,107)3,2,Real(ZEL(3,2),dp), "# Real(ZEL(3,2),dp)" 
Write(io_L,107)3,3,Real(ZEL(3,3),dp), "# Real(ZEL(3,3),dp)" 
Write(io_L,107)3,4,Real(ZEL(3,4),dp), "# Real(ZEL(3,4),dp)" 
Write(io_L,107)3,5,Real(ZEL(3,5),dp), "# Real(ZEL(3,5),dp)" 
Write(io_L,107)4,1,Real(ZEL(4,1),dp), "# Real(ZEL(4,1),dp)" 
Write(io_L,107)4,2,Real(ZEL(4,2),dp), "# Real(ZEL(4,2),dp)" 
Write(io_L,107)4,3,Real(ZEL(4,3),dp), "# Real(ZEL(4,3),dp)" 
Write(io_L,107)4,4,Real(ZEL(4,4),dp), "# Real(ZEL(4,4),dp)" 
Write(io_L,107)4,5,Real(ZEL(4,5),dp), "# Real(ZEL(4,5),dp)" 
Write(io_L,107)5,1,Real(ZEL(5,1),dp), "# Real(ZEL(5,1),dp)" 
Write(io_L,107)5,2,Real(ZEL(5,2),dp), "# Real(ZEL(5,2),dp)" 
Write(io_L,107)5,3,Real(ZEL(5,3),dp), "# Real(ZEL(5,3),dp)" 
Write(io_L,107)5,4,Real(ZEL(5,4),dp), "# Real(ZEL(5,4),dp)" 
Write(io_L,107)5,5,Real(ZEL(5,5),dp), "# Real(ZEL(5,5),dp)" 
If (MaxVal(Abs(AImag(ZEL))).gt.0._dp) Then 
Write(io_L,106) "Block IMUELMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZEL(1,1)), "# Aimag(ZEL(1,1))" 
Write(io_L,107)1,2,Aimag(ZEL(1,2)), "# Aimag(ZEL(1,2))" 
Write(io_L,107)1,3,Aimag(ZEL(1,3)), "# Aimag(ZEL(1,3))" 
Write(io_L,107)1,4,Aimag(ZEL(1,4)), "# Aimag(ZEL(1,4))" 
Write(io_L,107)1,5,Aimag(ZEL(1,5)), "# Aimag(ZEL(1,5))" 
Write(io_L,107)2,1,Aimag(ZEL(2,1)), "# Aimag(ZEL(2,1))" 
Write(io_L,107)2,2,Aimag(ZEL(2,2)), "# Aimag(ZEL(2,2))" 
Write(io_L,107)2,3,Aimag(ZEL(2,3)), "# Aimag(ZEL(2,3))" 
Write(io_L,107)2,4,Aimag(ZEL(2,4)), "# Aimag(ZEL(2,4))" 
Write(io_L,107)2,5,Aimag(ZEL(2,5)), "# Aimag(ZEL(2,5))" 
Write(io_L,107)3,1,Aimag(ZEL(3,1)), "# Aimag(ZEL(3,1))" 
Write(io_L,107)3,2,Aimag(ZEL(3,2)), "# Aimag(ZEL(3,2))" 
Write(io_L,107)3,3,Aimag(ZEL(3,3)), "# Aimag(ZEL(3,3))" 
Write(io_L,107)3,4,Aimag(ZEL(3,4)), "# Aimag(ZEL(3,4))" 
Write(io_L,107)3,5,Aimag(ZEL(3,5)), "# Aimag(ZEL(3,5))" 
Write(io_L,107)4,1,Aimag(ZEL(4,1)), "# Aimag(ZEL(4,1))" 
Write(io_L,107)4,2,Aimag(ZEL(4,2)), "# Aimag(ZEL(4,2))" 
Write(io_L,107)4,3,Aimag(ZEL(4,3)), "# Aimag(ZEL(4,3))" 
Write(io_L,107)4,4,Aimag(ZEL(4,4)), "# Aimag(ZEL(4,4))" 
Write(io_L,107)4,5,Aimag(ZEL(4,5)), "# Aimag(ZEL(4,5))" 
Write(io_L,107)5,1,Aimag(ZEL(5,1)), "# Aimag(ZEL(5,1))" 
Write(io_L,107)5,2,Aimag(ZEL(5,2)), "# Aimag(ZEL(5,2))" 
Write(io_L,107)5,3,Aimag(ZEL(5,3)), "# Aimag(ZEL(5,3))" 
Write(io_L,107)5,4,Aimag(ZEL(5,4)), "# Aimag(ZEL(5,4))" 
Write(io_L,107)5,5,Aimag(ZEL(5,5)), "# Aimag(ZEL(5,5))" 
End If 

Write(io_L,106) "Block UDLMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZDL(1,1),dp), "# Real(ZDL(1,1),dp)" 
Write(io_L,107)1,2,Real(ZDL(1,2),dp), "# Real(ZDL(1,2),dp)" 
Write(io_L,107)1,3,Real(ZDL(1,3),dp), "# Real(ZDL(1,3),dp)" 
Write(io_L,107)2,1,Real(ZDL(2,1),dp), "# Real(ZDL(2,1),dp)" 
Write(io_L,107)2,2,Real(ZDL(2,2),dp), "# Real(ZDL(2,2),dp)" 
Write(io_L,107)2,3,Real(ZDL(2,3),dp), "# Real(ZDL(2,3),dp)" 
Write(io_L,107)3,1,Real(ZDL(3,1),dp), "# Real(ZDL(3,1),dp)" 
Write(io_L,107)3,2,Real(ZDL(3,2),dp), "# Real(ZDL(3,2),dp)" 
Write(io_L,107)3,3,Real(ZDL(3,3),dp), "# Real(ZDL(3,3),dp)" 
If (MaxVal(Abs(AImag(ZDL))).gt.0._dp) Then 
Write(io_L,106) "Block IMUDLMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZDL(1,1)), "# Aimag(ZDL(1,1))" 
Write(io_L,107)1,2,Aimag(ZDL(1,2)), "# Aimag(ZDL(1,2))" 
Write(io_L,107)1,3,Aimag(ZDL(1,3)), "# Aimag(ZDL(1,3))" 
Write(io_L,107)2,1,Aimag(ZDL(2,1)), "# Aimag(ZDL(2,1))" 
Write(io_L,107)2,2,Aimag(ZDL(2,2)), "# Aimag(ZDL(2,2))" 
Write(io_L,107)2,3,Aimag(ZDL(2,3)), "# Aimag(ZDL(2,3))" 
Write(io_L,107)3,1,Aimag(ZDL(3,1)), "# Aimag(ZDL(3,1))" 
Write(io_L,107)3,2,Aimag(ZDL(3,2)), "# Aimag(ZDL(3,2))" 
Write(io_L,107)3,3,Aimag(ZDL(3,3)), "# Aimag(ZDL(3,3))" 
End If 

Write(io_L,106) "Block UDRMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZDR(1,1),dp), "# Real(ZDR(1,1),dp)" 
Write(io_L,107)1,2,Real(ZDR(1,2),dp), "# Real(ZDR(1,2),dp)" 
Write(io_L,107)1,3,Real(ZDR(1,3),dp), "# Real(ZDR(1,3),dp)" 
Write(io_L,107)2,1,Real(ZDR(2,1),dp), "# Real(ZDR(2,1),dp)" 
Write(io_L,107)2,2,Real(ZDR(2,2),dp), "# Real(ZDR(2,2),dp)" 
Write(io_L,107)2,3,Real(ZDR(2,3),dp), "# Real(ZDR(2,3),dp)" 
Write(io_L,107)3,1,Real(ZDR(3,1),dp), "# Real(ZDR(3,1),dp)" 
Write(io_L,107)3,2,Real(ZDR(3,2),dp), "# Real(ZDR(3,2),dp)" 
Write(io_L,107)3,3,Real(ZDR(3,3),dp), "# Real(ZDR(3,3),dp)" 
If (MaxVal(Abs(AImag(ZDR))).gt.0._dp) Then 
Write(io_L,106) "Block IMUDRMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZDR(1,1)), "# Aimag(ZDR(1,1))" 
Write(io_L,107)1,2,Aimag(ZDR(1,2)), "# Aimag(ZDR(1,2))" 
Write(io_L,107)1,3,Aimag(ZDR(1,3)), "# Aimag(ZDR(1,3))" 
Write(io_L,107)2,1,Aimag(ZDR(2,1)), "# Aimag(ZDR(2,1))" 
Write(io_L,107)2,2,Aimag(ZDR(2,2)), "# Aimag(ZDR(2,2))" 
Write(io_L,107)2,3,Aimag(ZDR(2,3)), "# Aimag(ZDR(2,3))" 
Write(io_L,107)3,1,Aimag(ZDR(3,1)), "# Aimag(ZDR(3,1))" 
Write(io_L,107)3,2,Aimag(ZDR(3,2)), "# Aimag(ZDR(3,2))" 
Write(io_L,107)3,3,Aimag(ZDR(3,3)), "# Aimag(ZDR(3,3))" 
End If 

Write(io_L,106) "Block UULMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZUL(1,1),dp), "# Real(ZUL(1,1),dp)" 
Write(io_L,107)1,2,Real(ZUL(1,2),dp), "# Real(ZUL(1,2),dp)" 
Write(io_L,107)1,3,Real(ZUL(1,3),dp), "# Real(ZUL(1,3),dp)" 
Write(io_L,107)2,1,Real(ZUL(2,1),dp), "# Real(ZUL(2,1),dp)" 
Write(io_L,107)2,2,Real(ZUL(2,2),dp), "# Real(ZUL(2,2),dp)" 
Write(io_L,107)2,3,Real(ZUL(2,3),dp), "# Real(ZUL(2,3),dp)" 
Write(io_L,107)3,1,Real(ZUL(3,1),dp), "# Real(ZUL(3,1),dp)" 
Write(io_L,107)3,2,Real(ZUL(3,2),dp), "# Real(ZUL(3,2),dp)" 
Write(io_L,107)3,3,Real(ZUL(3,3),dp), "# Real(ZUL(3,3),dp)" 
If (MaxVal(Abs(AImag(ZUL))).gt.0._dp) Then 
Write(io_L,106) "Block IMUULMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZUL(1,1)), "# Aimag(ZUL(1,1))" 
Write(io_L,107)1,2,Aimag(ZUL(1,2)), "# Aimag(ZUL(1,2))" 
Write(io_L,107)1,3,Aimag(ZUL(1,3)), "# Aimag(ZUL(1,3))" 
Write(io_L,107)2,1,Aimag(ZUL(2,1)), "# Aimag(ZUL(2,1))" 
Write(io_L,107)2,2,Aimag(ZUL(2,2)), "# Aimag(ZUL(2,2))" 
Write(io_L,107)2,3,Aimag(ZUL(2,3)), "# Aimag(ZUL(2,3))" 
Write(io_L,107)3,1,Aimag(ZUL(3,1)), "# Aimag(ZUL(3,1))" 
Write(io_L,107)3,2,Aimag(ZUL(3,2)), "# Aimag(ZUL(3,2))" 
Write(io_L,107)3,3,Aimag(ZUL(3,3)), "# Aimag(ZUL(3,3))" 
End If 

Write(io_L,106) "Block UURMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Real(ZUR(1,1),dp), "# Real(ZUR(1,1),dp)" 
Write(io_L,107)1,2,Real(ZUR(1,2),dp), "# Real(ZUR(1,2),dp)" 
Write(io_L,107)1,3,Real(ZUR(1,3),dp), "# Real(ZUR(1,3),dp)" 
Write(io_L,107)2,1,Real(ZUR(2,1),dp), "# Real(ZUR(2,1),dp)" 
Write(io_L,107)2,2,Real(ZUR(2,2),dp), "# Real(ZUR(2,2),dp)" 
Write(io_L,107)2,3,Real(ZUR(2,3),dp), "# Real(ZUR(2,3),dp)" 
Write(io_L,107)3,1,Real(ZUR(3,1),dp), "# Real(ZUR(3,1),dp)" 
Write(io_L,107)3,2,Real(ZUR(3,2),dp), "# Real(ZUR(3,2),dp)" 
Write(io_L,107)3,3,Real(ZUR(3,3),dp), "# Real(ZUR(3,3),dp)" 
If (MaxVal(Abs(AImag(ZUR))).gt.0._dp) Then 
Write(io_L,106) "Block IMUURMIX Q=",Q,"# ()" 
Write(io_L,107)1,1,Aimag(ZUR(1,1)), "# Aimag(ZUR(1,1))" 
Write(io_L,107)1,2,Aimag(ZUR(1,2)), "# Aimag(ZUR(1,2))" 
Write(io_L,107)1,3,Aimag(ZUR(1,3)), "# Aimag(ZUR(1,3))" 
Write(io_L,107)2,1,Aimag(ZUR(2,1)), "# Aimag(ZUR(2,1))" 
Write(io_L,107)2,2,Aimag(ZUR(2,2)), "# Aimag(ZUR(2,2))" 
Write(io_L,107)2,3,Aimag(ZUR(2,3)), "# Aimag(ZUR(2,3))" 
Write(io_L,107)3,1,Aimag(ZUR(3,1)), "# Aimag(ZUR(3,1))" 
Write(io_L,107)3,2,Aimag(ZUR(3,2)), "# Aimag(ZUR(3,2))" 
Write(io_L,107)3,3,Aimag(ZUR(3,3)), "# Aimag(ZUR(3,3))" 
End If 

Write(io_L,100) "Block SPheno # SPheno internal parameters " 
Write(io_L,102) 1,1.*ErrorLevel,"# ErrorLevel"

If (SPA_convention) Then
Write(io_L,102) 2,1.,"# SPA_conventions"
Else
Write(io_L,102) 2,0.,"# SPA_conventions"
End if

If (L_BR) Then
Write(io_L,102) 11,1.,"# Branching ratios"
Else
Write(io_L,102) 11,0.,"# Branching ratios"
End if


If (Enable3BDecays) Then
Write(io_L,102) 13,1.,"# 3 Body decays"
Else
Write(io_L,102) 13,0.,"# 3 Body decays"
End if


Write(io_L,102) 31,m_GUT,"# GUT scale"

Write(io_L,102) 33,Q,"# Renormalization scale"

Write(io_L,102) 34,delta_mass,"# Precision"

Write(io_L,102) 35,1.*n_run,"# Iterations"


If (TwoLoopRGE) Then
Write(io_L,102) 38,2.,"# RGE level"
Else
Write(io_L,102) 38,1.,"# RGE level"
End if

Write(io_L,102) 40,Alpha,"# Alpha"

Write(io_L,102) 41,gamZ,"# Gamma_Z"

Write(io_L,102) 42,gamW,"# Gamma_W"

If (RotateNegativeFermionMasses) Then
Write(io_L,102) 50,1.,"# Rotate negative fermion masses"
Else
Write(io_L,102) 50,0.,"# Rotate negative fermion masses"
End if


If (SwitchToSCKM) Then
Write(io_L,102) 51,1.,"# Switch to SCKM matrix"
Else
Write(io_L,102) 51,0.,"# Switch to SCKM matrix"
End if


If (IgnoreNegativeMasses) Then
Write(io_L,102) 52,1.,"# Ignore negative masses"
Else
Write(io_L,102) 52,0.,"# Ignore negative masses"
End if


If (IgnoreNegativeMassesMZ) Then
Write(io_L,102) 53,1.,"# Ignore negative masses at MZ"
Else
Write(io_L,102) 53,0.,"# Ignore negative masses at MZ"
End if


If (CalculateOneLoopMasses) Then
Write(io_L,102) 55,1.,"# Calculate one loop masses"
Else
Write(io_L,102) 55,0.,"# Calculate one loop masses"
End if


If (CalculateTwoLoopHiggsMasses) Then
Write(io_L,102) 56,1.,"# Calculate two-loop Higgs masses"
Else
Write(io_L,102) 56,0.,"# Calculate two-loop Higgs masses"
End if

If (CalculateLowEnergy) Then
Write(io_L,102) 57,1.,"# Calculate low energy"
Else
Write(io_L,102) 57,0.,"# Calculate low energy"
End if

If (KineticMixing) Then
Write(io_L,102) 60,1.,"# Include kinetic mixing"
Else
Write(io_L,102) 60,0.,"# Include kinetic mixing"
End if

Write(io_L,102) 65,1.*SolutionTadpoleNr,"# Solution of tadpole equation"
 

SELECT CASE (TwoLoopMethod)
CASE ( 1 ) 
 Write(io_L,102) 8,1.,"# Two-Loop Method: purely numerical " 
CASE ( 2 ) 
 Write(io_L,102) 8,2.,"# Two-Loop Method: semi-analytical " 
CASE ( 3 ) 
 Write(io_L,102) 8,3.,"# Two-Loop Method: diagrammatic " 
CASE ( 8 ) 
 Write(io_L,102) 8,8.,"# Two-Loop Method: MSSM+ alpha_s alpha_t routines " 
CASE ( 9 ) 
 Write(io_L,102) 8,9.,"# Two-Loop Method: MSSM+ routines " 
CASE DEFAULT 
 Write(io_L,102) 8,0.,"# Two-Loop Method: no two loop calculation " 
END SELECT 
If (GaugelessLimit) Then 
 Write(io_L,102) 9,1.,"# Gauge-less limit" 
Else 
 Write(io_L,102) 9,0.,"# Gauge-less limit" 
End if 
Write(io_L,102) 400,hstep_pn,"# Step-size for purely-numerical methode for 2-loop calculation" 
Write(io_L,102) 401,hstep_sa,"# Step-size for semi-analytical methode for 2-loop calculation" 
Write(io_L,102) 410,err2L,"# indicative error in numerical derivation" 

 
If(Write_WHIZARD) Call WriteWHIZARD 
 
If(Write_HiggsBounds) Call WriteHiggsBounds
 

 
If (L_BR) Then 
Write(io_L,100) "Block HiggsLHC7 # Higgs production cross section at LHC7 [pb] " 
Do i1=1,8
CurrentPDG2(1) = Abs(PDGhh(i1)) 
If (CS_Higgs_LHC(1,i1,1).gt.0._dp) Write(io_L,119) 1, CurrentPDG2(1), CS_Higgs_LHC(1,i1,1), " # Gluon fusion " 
If (CS_Higgs_LHC(1,i1,2).gt.0._dp) Write(io_L,119) 2, CurrentPDG2(1), CS_Higgs_LHC(1,i1,2), " # Vector boson fusion " 
If (CS_Higgs_LHC(1,i1,3).gt.0._dp) Write(io_L,119) 3, CurrentPDG2(1), CS_Higgs_LHC(1,i1,3), " # W-H production " 
If (CS_Higgs_LHC(1,i1,4).gt.0._dp) Write(io_L,119) 4, CurrentPDG2(1), CS_Higgs_LHC(1,i1,4), " # Z-H production " 
If (CS_Higgs_LHC(1,i1,5).gt.0._dp) Write(io_L,119) 5, CurrentPDG2(1), CS_Higgs_LHC(1,i1,5), " # t-t-H production " 
End Do 
Write(io_L,100) "Block HiggsLHC8 # Higgs production cross section at LHC8 [pb] " 
Do i1=1,8
CurrentPDG2(1) = Abs(PDGhh(i1)) 
If (CS_Higgs_LHC(2,i1,1).gt.0._dp) Write(io_L,119) 1, CurrentPDG2(1), CS_Higgs_LHC(2,i1,1), " # Gluon fusion " 
If (CS_Higgs_LHC(2,i1,2).gt.0._dp) Write(io_L,119) 2, CurrentPDG2(1), CS_Higgs_LHC(2,i1,2), " # Vector boson fusion " 
If (CS_Higgs_LHC(2,i1,3).gt.0._dp) Write(io_L,119) 3, CurrentPDG2(1), CS_Higgs_LHC(2,i1,3), " # W-H production " 
If (CS_Higgs_LHC(2,i1,4).gt.0._dp) Write(io_L,119) 4, CurrentPDG2(1), CS_Higgs_LHC(2,i1,4), " # Z-H production " 
If (CS_Higgs_LHC(2,i1,5).gt.0._dp) Write(io_L,119) 5, CurrentPDG2(1), CS_Higgs_LHC(2,i1,5), " # t-t-H production " 
End Do 
If (WriteEffHiggsCouplingRatios) Then 
Write(io_L,100) "Block HiggsBoundsInputHiggsCouplingsFermions # " 
Write(io_L,1101) rHB_S_S_Fd(1,3),rHB_S_P_Fd(1,3), 3 ,25,5,5, " # h_1 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(1,2),rHB_S_P_Fd(1,2), 3 ,25,3,3, " # h_1 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(1,3),rHB_S_P_Fu(1,3), 3 ,25,6,6, " # h_1 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(1,2),rHB_S_P_Fu(1,2),3 ,25,4,4, " # h_1 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(1,3),rHB_S_P_Cha(1,3), 3 ,25,15,15, " # h_1 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(1,2),rHB_S_P_Cha(1,2), 3 ,25,13,13, " # h_1 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(2,3),rHB_S_P_Fd(2,3), 3 ,35,5,5, " # h_2 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(2,2),rHB_S_P_Fd(2,2), 3 ,35,3,3, " # h_2 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(2,3),rHB_S_P_Fu(2,3), 3 ,35,6,6, " # h_2 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(2,2),rHB_S_P_Fu(2,2),3 ,35,4,4, " # h_2 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(2,3),rHB_S_P_Cha(2,3), 3 ,35,15,15, " # h_2 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(2,2),rHB_S_P_Cha(2,2), 3 ,35,13,13, " # h_2 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(3,3),rHB_S_P_Fd(3,3), 3 ,1000012,5,5, " # h_3 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(3,2),rHB_S_P_Fd(3,2), 3 ,1000012,3,3, " # h_3 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(3,3),rHB_S_P_Fu(3,3), 3 ,1000012,6,6, " # h_3 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(3,2),rHB_S_P_Fu(3,2),3 ,1000012,4,4, " # h_3 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(3,3),rHB_S_P_Cha(3,3), 3 ,1000012,15,15, " # h_3 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(3,2),rHB_S_P_Cha(3,2), 3 ,1000012,13,13, " # h_3 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(4,3),rHB_S_P_Fd(4,3), 3 ,1000014,5,5, " # h_4 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(4,2),rHB_S_P_Fd(4,2), 3 ,1000014,3,3, " # h_4 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(4,3),rHB_S_P_Fu(4,3), 3 ,1000014,6,6, " # h_4 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(4,2),rHB_S_P_Fu(4,2),3 ,1000014,4,4, " # h_4 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(4,3),rHB_S_P_Cha(4,3), 3 ,1000014,15,15, " # h_4 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(4,2),rHB_S_P_Cha(4,2), 3 ,1000014,13,13, " # h_4 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(5,3),rHB_S_P_Fd(5,3), 3 ,1000016,5,5, " # h_5 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(5,2),rHB_S_P_Fd(5,2), 3 ,1000016,3,3, " # h_5 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(5,3),rHB_S_P_Fu(5,3), 3 ,1000016,6,6, " # h_5 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(5,2),rHB_S_P_Fu(5,2),3 ,1000016,4,4, " # h_5 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(5,3),rHB_S_P_Cha(5,3), 3 ,1000016,15,15, " # h_5 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(5,2),rHB_S_P_Cha(5,2), 3 ,1000016,13,13, " # h_5 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(6,3),rHB_S_P_Fd(6,3), 3 ,2000012,5,5, " # h_6 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(6,2),rHB_S_P_Fd(6,2), 3 ,2000012,3,3, " # h_6 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(6,3),rHB_S_P_Fu(6,3), 3 ,2000012,6,6, " # h_6 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(6,2),rHB_S_P_Fu(6,2),3 ,2000012,4,4, " # h_6 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(6,3),rHB_S_P_Cha(6,3), 3 ,2000012,15,15, " # h_6 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(6,2),rHB_S_P_Cha(6,2), 3 ,2000012,13,13, " # h_6 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(7,3),rHB_S_P_Fd(7,3), 3 ,2000014,5,5, " # h_7 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(7,2),rHB_S_P_Fd(7,2), 3 ,2000014,3,3, " # h_7 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(7,3),rHB_S_P_Fu(7,3), 3 ,2000014,6,6, " # h_7 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(7,2),rHB_S_P_Fu(7,2),3 ,2000014,4,4, " # h_7 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(7,3),rHB_S_P_Cha(7,3), 3 ,2000014,15,15, " # h_7 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(7,2),rHB_S_P_Cha(7,2), 3 ,2000014,13,13, " # h_7 mu mu coupling  " 
Write(io_L,1101) rHB_S_S_Fd(8,3),rHB_S_P_Fd(8,3), 3 ,2000016,5,5, " # h_8 b b coupling " 
Write(io_L,1101) rHB_S_S_Fd(8,2),rHB_S_P_Fd(8,2), 3 ,2000016,3,3, " # h_8 s s coupling " 
Write(io_L,1101) rHB_S_S_Fu(8,3),rHB_S_P_Fu(8,3), 3 ,2000016,6,6, " # h_8 t t coupling  " 
Write(io_L,1101) rHB_S_S_Fu(8,2),rHB_S_P_Fu(8,2),3 ,2000016,4,4, " # h_8 c c coupling " 
Write(io_L,1101) rHB_S_S_Cha(8,3),rHB_S_P_Cha(8,3), 3 ,2000016,15,15, " # h_8 tau tau coupling " 
Write(io_L,1101) rHB_S_S_Cha(8,2),rHB_S_P_Cha(8,2), 3 ,2000016,13,13, " # h_8 mu mu coupling  " 
Write(io_L,1101) rHB_P_S_Fd(2,3),rHB_P_P_Fd(2,3), 3 ,36,5,5, " # A_2 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(2,2),rHB_P_P_Fd(2,2), 3 ,36,3,3, " # A_2 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(2,3),rHB_P_P_Fu(2,3), 3 ,36,6,6, " # A_2 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(2,2),rHB_P_P_Fu(2,2), 3 ,36,4,4, " # A_2 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(2,3),rHB_P_P_Cha(2,3), 3 ,36,15,15, " # A_2 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(2,2),rHB_P_P_Cha(2,2), 3 ,36,13,13, " # A_2 mu mu coupling " 
Write(io_L,1101) rHB_P_S_Fd(3,3),rHB_P_P_Fd(3,3), 3 ,1000017,5,5, " # A_3 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(3,2),rHB_P_P_Fd(3,2), 3 ,1000017,3,3, " # A_3 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(3,3),rHB_P_P_Fu(3,3), 3 ,1000017,6,6, " # A_3 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(3,2),rHB_P_P_Fu(3,2), 3 ,1000017,4,4, " # A_3 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(3,3),rHB_P_P_Cha(3,3), 3 ,1000017,15,15, " # A_3 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(3,2),rHB_P_P_Cha(3,2), 3 ,1000017,13,13, " # A_3 mu mu coupling " 
Write(io_L,1101) rHB_P_S_Fd(4,3),rHB_P_P_Fd(4,3), 3 ,1000018,5,5, " # A_4 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(4,2),rHB_P_P_Fd(4,2), 3 ,1000018,3,3, " # A_4 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(4,3),rHB_P_P_Fu(4,3), 3 ,1000018,6,6, " # A_4 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(4,2),rHB_P_P_Fu(4,2), 3 ,1000018,4,4, " # A_4 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(4,3),rHB_P_P_Cha(4,3), 3 ,1000018,15,15, " # A_4 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(4,2),rHB_P_P_Cha(4,2), 3 ,1000018,13,13, " # A_4 mu mu coupling " 
Write(io_L,1101) rHB_P_S_Fd(5,3),rHB_P_P_Fd(5,3), 3 ,1000019,5,5, " # A_5 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(5,2),rHB_P_P_Fd(5,2), 3 ,1000019,3,3, " # A_5 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(5,3),rHB_P_P_Fu(5,3), 3 ,1000019,6,6, " # A_5 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(5,2),rHB_P_P_Fu(5,2), 3 ,1000019,4,4, " # A_5 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(5,3),rHB_P_P_Cha(5,3), 3 ,1000019,15,15, " # A_5 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(5,2),rHB_P_P_Cha(5,2), 3 ,1000019,13,13, " # A_5 mu mu coupling " 
Write(io_L,1101) rHB_P_S_Fd(6,3),rHB_P_P_Fd(6,3), 3 ,2000018,5,5, " # A_6 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(6,2),rHB_P_P_Fd(6,2), 3 ,2000018,3,3, " # A_6 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(6,3),rHB_P_P_Fu(6,3), 3 ,2000018,6,6, " # A_6 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(6,2),rHB_P_P_Fu(6,2), 3 ,2000018,4,4, " # A_6 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(6,3),rHB_P_P_Cha(6,3), 3 ,2000018,15,15, " # A_6 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(6,2),rHB_P_P_Cha(6,2), 3 ,2000018,13,13, " # A_6 mu mu coupling " 
Write(io_L,1101) rHB_P_S_Fd(7,3),rHB_P_P_Fd(7,3), 3 ,2000019,5,5, " # A_7 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(7,2),rHB_P_P_Fd(7,2), 3 ,2000019,3,3, " # A_7 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(7,3),rHB_P_P_Fu(7,3), 3 ,2000019,6,6, " # A_7 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(7,2),rHB_P_P_Fu(7,2), 3 ,2000019,4,4, " # A_7 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(7,3),rHB_P_P_Cha(7,3), 3 ,2000019,15,15, " # A_7 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(7,2),rHB_P_P_Cha(7,2), 3 ,2000019,13,13, " # A_7 mu mu coupling " 
Write(io_L,1101) rHB_P_S_Fd(8,3),rHB_P_P_Fd(8,3), 3 ,2000020,5,5, " # A_8 b b coupling " 
Write(io_L,1101) rHB_P_S_Fd(8,2),rHB_P_P_Fd(8,2), 3 ,2000020,3,3, " # A_8 s s coupling " 
Write(io_L,1101) rHB_P_S_Fu(8,3),rHB_P_P_Fu(8,3), 3 ,2000020,6,6, " # A_8 t t coupling "  
Write(io_L,1101) rHB_P_S_Fu(8,2),rHB_P_P_Fu(8,2), 3 ,2000020,4,4, " # A_8 c c coupling " 
Write(io_L,1101) rHB_P_S_Cha(8,3),rHB_P_P_Cha(8,3), 3 ,2000020,15,15, " # A_8 tau tau coupling " 
Write(io_L,1101) rHB_P_S_Cha(8,2),rHB_P_P_Cha(8,2), 3 ,2000020,13,13, " # A_8 mu mu coupling " 
Write(io_L,100) "Block HiggsBoundsInputHiggsCouplingsBosons # " 
Write(io_L,1102) rHB_S_VWm(1),3 ,25,24,24, " # h_1 W W coupling " 
Write(io_L,1102) rHB_S_VZ(1),3 ,25,23,23, " # h_1 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,25,23,22, " # h_1 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(1),dp),3 ,25,22,22, " # h_1 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(1),dp),3 ,25,21,21, " # h_1 g g coupling " 
Write(io_L,1103) 0._dp,4 ,25,21,21,23, " # h_1 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(2),3 ,35,24,24, " # h_2 W W coupling " 
Write(io_L,1102) rHB_S_VZ(2),3 ,35,23,23, " # h_2 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,35,23,22, " # h_2 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(2),dp),3 ,35,22,22, " # h_2 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(2),dp),3 ,35,21,21, " # h_2 g g coupling " 
Write(io_L,1103) 0._dp,4 ,35,21,21,23, " # h_2 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(3),3 ,1000012,24,24, " # h_3 W W coupling " 
Write(io_L,1102) rHB_S_VZ(3),3 ,1000012,23,23, " # h_3 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,1000012,23,22, " # h_3 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(3),dp),3 ,1000012,22,22, " # h_3 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(3),dp),3 ,1000012,21,21, " # h_3 g g coupling " 
Write(io_L,1103) 0._dp,4 ,1000012,21,21,23, " # h_3 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(4),3 ,1000014,24,24, " # h_4 W W coupling " 
Write(io_L,1102) rHB_S_VZ(4),3 ,1000014,23,23, " # h_4 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,1000014,23,22, " # h_4 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(4),dp),3 ,1000014,22,22, " # h_4 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(4),dp),3 ,1000014,21,21, " # h_4 g g coupling " 
Write(io_L,1103) 0._dp,4 ,1000014,21,21,23, " # h_4 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(5),3 ,1000016,24,24, " # h_5 W W coupling " 
Write(io_L,1102) rHB_S_VZ(5),3 ,1000016,23,23, " # h_5 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,1000016,23,22, " # h_5 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(5),dp),3 ,1000016,22,22, " # h_5 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(5),dp),3 ,1000016,21,21, " # h_5 g g coupling " 
Write(io_L,1103) 0._dp,4 ,1000016,21,21,23, " # h_5 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(6),3 ,2000012,24,24, " # h_6 W W coupling " 
Write(io_L,1102) rHB_S_VZ(6),3 ,2000012,23,23, " # h_6 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,2000012,23,22, " # h_6 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(6),dp),3 ,2000012,22,22, " # h_6 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(6),dp),3 ,2000012,21,21, " # h_6 g g coupling " 
Write(io_L,1103) 0._dp,4 ,2000012,21,21,23, " # h_6 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(7),3 ,2000014,24,24, " # h_7 W W coupling " 
Write(io_L,1102) rHB_S_VZ(7),3 ,2000014,23,23, " # h_7 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,2000014,23,22, " # h_7 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(7),dp),3 ,2000014,22,22, " # h_7 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(7),dp),3 ,2000014,21,21, " # h_7 g g coupling " 
Write(io_L,1103) 0._dp,4 ,2000014,21,21,23, " # h_7 g g Z coupling " 
Write(io_L,1102) rHB_S_VWm(8),3 ,2000016,24,24, " # h_8 W W coupling " 
Write(io_L,1102) rHB_S_VZ(8),3 ,2000016,23,23, " # h_8 Z Z coupling  " 
Write(io_L,1102) 0._dp ,3 ,2000016,23,22, " # h_8 Z gamma coupling " 
Write(io_L,1102) Real(ratioPP(8),dp),3 ,2000016,22,22, " # h_8 gamma gamma coupling " 
Write(io_L,1102) Real(ratioGG(8),dp),3 ,2000016,21,21, " # h_8 g g coupling " 
Write(io_L,1103) 0._dp,4 ,2000016,21,21,23, " # h_8 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(2),3 ,36,24,24, " # A_2 W W coupling " 
Write(io_L,1102) rHB_P_VZ(2),3 ,36,23,23, " # A_2 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,36,23,22, " # A_2 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(2),dp),3 ,36,22,22, " # A_2 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(2),dp),3 ,36,21,21, " # A_2 g g coupling " 
Write(io_L,1103) 0._dp,4 ,36,21,21,23, " # A_2 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(3),3 ,1000017,24,24, " # A_3 W W coupling " 
Write(io_L,1102) rHB_P_VZ(3),3 ,1000017,23,23, " # A_3 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,1000017,23,22, " # A_3 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(3),dp),3 ,1000017,22,22, " # A_3 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(3),dp),3 ,1000017,21,21, " # A_3 g g coupling " 
Write(io_L,1103) 0._dp,4 ,1000017,21,21,23, " # A_3 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(4),3 ,1000018,24,24, " # A_4 W W coupling " 
Write(io_L,1102) rHB_P_VZ(4),3 ,1000018,23,23, " # A_4 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,1000018,23,22, " # A_4 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(4),dp),3 ,1000018,22,22, " # A_4 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(4),dp),3 ,1000018,21,21, " # A_4 g g coupling " 
Write(io_L,1103) 0._dp,4 ,1000018,21,21,23, " # A_4 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(5),3 ,1000019,24,24, " # A_5 W W coupling " 
Write(io_L,1102) rHB_P_VZ(5),3 ,1000019,23,23, " # A_5 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,1000019,23,22, " # A_5 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(5),dp),3 ,1000019,22,22, " # A_5 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(5),dp),3 ,1000019,21,21, " # A_5 g g coupling " 
Write(io_L,1103) 0._dp,4 ,1000019,21,21,23, " # A_5 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(6),3 ,2000018,24,24, " # A_6 W W coupling " 
Write(io_L,1102) rHB_P_VZ(6),3 ,2000018,23,23, " # A_6 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,2000018,23,22, " # A_6 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(6),dp),3 ,2000018,22,22, " # A_6 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(6),dp),3 ,2000018,21,21, " # A_6 g g coupling " 
Write(io_L,1103) 0._dp,4 ,2000018,21,21,23, " # A_6 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(7),3 ,2000019,24,24, " # A_7 W W coupling " 
Write(io_L,1102) rHB_P_VZ(7),3 ,2000019,23,23, " # A_7 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,2000019,23,22, " # A_7 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(7),dp),3 ,2000019,22,22, " # A_7 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(7),dp),3 ,2000019,21,21, " # A_7 g g coupling " 
Write(io_L,1103) 0._dp,4 ,2000019,21,21,23, " # A_7 g g Z coupling " 
Write(io_L,1102) rHB_P_VWm(8),3 ,2000020,24,24, " # A_8 W W coupling " 
Write(io_L,1102) rHB_P_VZ(8),3 ,2000020,23,23, " # A_8 Z Z coupling " 
Write(io_L,1102) 0._dp ,3 ,2000020,23,22, " # A_8 Z gamma coupling " 
Write(io_L,1102) Real(ratioPPP(8),dp),3 ,2000020,22,22, " # A_8 gamma gamma coupling " 
Write(io_L,1102) Real(ratioPGG(8),dp),3 ,2000020,21,21, " # A_8 g g coupling " 
Write(io_L,1103) 0._dp,4 ,2000020,21,21,23, " # A_8 g g Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(1,1), dp),3 ,25,25,23, " # h_1 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,2), dp),3 ,25,35,23, " # h_1 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,3), dp),3 ,25,1000012,23, " # h_1 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,4), dp),3 ,25,1000014,23, " # h_1 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,5), dp),3 ,25,1000016,23, " # h_1 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,6), dp),3 ,25,2000012,23, " # h_1 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,7), dp),3 ,25,2000014,23, " # h_1 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(1,8), dp),3 ,25,2000016,23, " # h_1 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,1), dp),3 ,25,36,23, " # h_1 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,1), dp),3 ,25,1000017,23, " # h_1 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,1), dp),3 ,25,1000018,23, " # h_1 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,1), dp),3 ,25,1000019,23, " # h_1 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,1), dp),3 ,25,2000018,23, " # h_1 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,1), dp),3 ,25,2000019,23, " # h_1 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,1), dp),3 ,25,2000020,23, " # h_1 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(2,1), dp),3 ,35,25,23, " # h_2 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,2), dp),3 ,35,35,23, " # h_2 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,3), dp),3 ,35,1000012,23, " # h_2 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,4), dp),3 ,35,1000014,23, " # h_2 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,5), dp),3 ,35,1000016,23, " # h_2 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,6), dp),3 ,35,2000012,23, " # h_2 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,7), dp),3 ,35,2000014,23, " # h_2 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(2,8), dp),3 ,35,2000016,23, " # h_2 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,2), dp),3 ,35,36,23, " # h_2 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,2), dp),3 ,35,1000017,23, " # h_2 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,2), dp),3 ,35,1000018,23, " # h_2 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,2), dp),3 ,35,1000019,23, " # h_2 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,2), dp),3 ,35,2000018,23, " # h_2 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,2), dp),3 ,35,2000019,23, " # h_2 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,2), dp),3 ,35,2000020,23, " # h_2 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(3,1), dp),3 ,1000012,25,23, " # h_3 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,2), dp),3 ,1000012,35,23, " # h_3 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,3), dp),3 ,1000012,1000012,23, " # h_3 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,4), dp),3 ,1000012,1000014,23, " # h_3 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,5), dp),3 ,1000012,1000016,23, " # h_3 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,6), dp),3 ,1000012,2000012,23, " # h_3 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,7), dp),3 ,1000012,2000014,23, " # h_3 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(3,8), dp),3 ,1000012,2000016,23, " # h_3 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,3), dp),3 ,1000012,36,23, " # h_3 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,3), dp),3 ,1000012,1000017,23, " # h_3 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,3), dp),3 ,1000012,1000018,23, " # h_3 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,3), dp),3 ,1000012,1000019,23, " # h_3 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,3), dp),3 ,1000012,2000018,23, " # h_3 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,3), dp),3 ,1000012,2000019,23, " # h_3 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,3), dp),3 ,1000012,2000020,23, " # h_3 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(4,1), dp),3 ,1000014,25,23, " # h_4 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,2), dp),3 ,1000014,35,23, " # h_4 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,3), dp),3 ,1000014,1000012,23, " # h_4 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,4), dp),3 ,1000014,1000014,23, " # h_4 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,5), dp),3 ,1000014,1000016,23, " # h_4 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,6), dp),3 ,1000014,2000012,23, " # h_4 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,7), dp),3 ,1000014,2000014,23, " # h_4 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(4,8), dp),3 ,1000014,2000016,23, " # h_4 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,4), dp),3 ,1000014,36,23, " # h_4 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,4), dp),3 ,1000014,1000017,23, " # h_4 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,4), dp),3 ,1000014,1000018,23, " # h_4 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,4), dp),3 ,1000014,1000019,23, " # h_4 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,4), dp),3 ,1000014,2000018,23, " # h_4 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,4), dp),3 ,1000014,2000019,23, " # h_4 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,4), dp),3 ,1000014,2000020,23, " # h_4 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(5,1), dp),3 ,1000016,25,23, " # h_5 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,2), dp),3 ,1000016,35,23, " # h_5 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,3), dp),3 ,1000016,1000012,23, " # h_5 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,4), dp),3 ,1000016,1000014,23, " # h_5 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,5), dp),3 ,1000016,1000016,23, " # h_5 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,6), dp),3 ,1000016,2000012,23, " # h_5 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,7), dp),3 ,1000016,2000014,23, " # h_5 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(5,8), dp),3 ,1000016,2000016,23, " # h_5 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,5), dp),3 ,1000016,36,23, " # h_5 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,5), dp),3 ,1000016,1000017,23, " # h_5 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,5), dp),3 ,1000016,1000018,23, " # h_5 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,5), dp),3 ,1000016,1000019,23, " # h_5 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,5), dp),3 ,1000016,2000018,23, " # h_5 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,5), dp),3 ,1000016,2000019,23, " # h_5 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,5), dp),3 ,1000016,2000020,23, " # h_5 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(6,1), dp),3 ,2000012,25,23, " # h_6 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,2), dp),3 ,2000012,35,23, " # h_6 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,3), dp),3 ,2000012,1000012,23, " # h_6 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,4), dp),3 ,2000012,1000014,23, " # h_6 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,5), dp),3 ,2000012,1000016,23, " # h_6 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,6), dp),3 ,2000012,2000012,23, " # h_6 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,7), dp),3 ,2000012,2000014,23, " # h_6 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(6,8), dp),3 ,2000012,2000016,23, " # h_6 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,6), dp),3 ,2000012,36,23, " # h_6 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,6), dp),3 ,2000012,1000017,23, " # h_6 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,6), dp),3 ,2000012,1000018,23, " # h_6 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,6), dp),3 ,2000012,1000019,23, " # h_6 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,6), dp),3 ,2000012,2000018,23, " # h_6 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,6), dp),3 ,2000012,2000019,23, " # h_6 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,6), dp),3 ,2000012,2000020,23, " # h_6 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(7,1), dp),3 ,2000014,25,23, " # h_7 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,2), dp),3 ,2000014,35,23, " # h_7 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,3), dp),3 ,2000014,1000012,23, " # h_7 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,4), dp),3 ,2000014,1000014,23, " # h_7 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,5), dp),3 ,2000014,1000016,23, " # h_7 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,6), dp),3 ,2000014,2000012,23, " # h_7 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,7), dp),3 ,2000014,2000014,23, " # h_7 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(7,8), dp),3 ,2000014,2000016,23, " # h_7 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,7), dp),3 ,2000014,36,23, " # h_7 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,7), dp),3 ,2000014,1000017,23, " # h_7 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,7), dp),3 ,2000014,1000018,23, " # h_7 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,7), dp),3 ,2000014,1000019,23, " # h_7 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,7), dp),3 ,2000014,2000018,23, " # h_7 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,7), dp),3 ,2000014,2000019,23, " # h_7 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,7), dp),3 ,2000014,2000020,23, " # h_7 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_H_H_Z(8,1), dp),3 ,2000016,25,23, " # h_8 h_1 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,2), dp),3 ,2000016,35,23, " # h_8 h_2 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,3), dp),3 ,2000016,1000012,23, " # h_8 h_3 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,4), dp),3 ,2000016,1000014,23, " # h_8 h_4 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,5), dp),3 ,2000016,1000016,23, " # h_8 h_5 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,6), dp),3 ,2000016,2000012,23, " # h_8 h_6 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,7), dp),3 ,2000016,2000014,23, " # h_8 h_7 Z coupling  "
Write(io_L,1102) Real(CPL_H_H_Z(8,8), dp),3 ,2000016,2000016,23, " # h_8 h_8 Z coupling  "
Write(io_L,1102) Real(CPL_A_H_Z(2,8), dp),3 ,2000016,36,23, " # h_8 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(3,8), dp),3 ,2000016,1000017,23, " # h_8 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(4,8), dp),3 ,2000016,1000018,23, " # h_8 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(5,8), dp),3 ,2000016,1000019,23, " # h_8 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(6,8), dp),3 ,2000016,2000018,23, " # h_8 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(7,8), dp),3 ,2000016,2000019,23, " # h_8 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_H_Z(8,8), dp),3 ,2000016,2000020,23, " # h_8 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,2), dp),3 ,36,36,23, " # A_2 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,2), dp),3 ,36,1000017,23, " # A_2 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,2), dp),3 ,36,1000018,23, " # A_2 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,2), dp),3 ,36,1000019,23, " # A_2 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,2), dp),3 ,36,2000018,23, " # A_2 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,2), dp),3 ,36,2000019,23, " # A_2 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,2), dp),3 ,36,2000020,23, " # A_2 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,3), dp),3 ,1000017,36,23, " # A_3 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,3), dp),3 ,1000017,1000017,23, " # A_3 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,3), dp),3 ,1000017,1000018,23, " # A_3 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,3), dp),3 ,1000017,1000019,23, " # A_3 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,3), dp),3 ,1000017,2000018,23, " # A_3 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,3), dp),3 ,1000017,2000019,23, " # A_3 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,3), dp),3 ,1000017,2000020,23, " # A_3 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,4), dp),3 ,1000018,36,23, " # A_4 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,4), dp),3 ,1000018,1000017,23, " # A_4 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,4), dp),3 ,1000018,1000018,23, " # A_4 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,4), dp),3 ,1000018,1000019,23, " # A_4 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,4), dp),3 ,1000018,2000018,23, " # A_4 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,4), dp),3 ,1000018,2000019,23, " # A_4 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,4), dp),3 ,1000018,2000020,23, " # A_4 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,5), dp),3 ,1000019,36,23, " # A_5 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,5), dp),3 ,1000019,1000017,23, " # A_5 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,5), dp),3 ,1000019,1000018,23, " # A_5 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,5), dp),3 ,1000019,1000019,23, " # A_5 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,5), dp),3 ,1000019,2000018,23, " # A_5 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,5), dp),3 ,1000019,2000019,23, " # A_5 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,5), dp),3 ,1000019,2000020,23, " # A_5 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,6), dp),3 ,2000018,36,23, " # A_6 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,6), dp),3 ,2000018,1000017,23, " # A_6 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,6), dp),3 ,2000018,1000018,23, " # A_6 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,6), dp),3 ,2000018,1000019,23, " # A_6 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,6), dp),3 ,2000018,2000018,23, " # A_6 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,6), dp),3 ,2000018,2000019,23, " # A_6 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,6), dp),3 ,2000018,2000020,23, " # A_6 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,7), dp),3 ,2000019,36,23, " # A_7 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,7), dp),3 ,2000019,1000017,23, " # A_7 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,7), dp),3 ,2000019,1000018,23, " # A_7 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,7), dp),3 ,2000019,1000019,23, " # A_7 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,7), dp),3 ,2000019,2000018,23, " # A_7 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,7), dp),3 ,2000019,2000019,23, " # A_7 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,7), dp),3 ,2000019,2000020,23, " # A_7 A_8 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(2,8), dp),3 ,2000020,36,23, " # A_8 A_2 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(3,8), dp),3 ,2000020,1000017,23, " # A_8 A_3 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(4,8), dp),3 ,2000020,1000018,23, " # A_8 A_4 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(5,8), dp),3 ,2000020,1000019,23, " # A_8 A_5 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(6,8), dp),3 ,2000020,2000018,23, " # A_8 A_6 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(7,8), dp),3 ,2000020,2000019,23, " # A_8 A_7 Z coupling " 
Write(io_L,1102) Real(CPL_A_A_Z(8,8), dp),3 ,2000020,2000020,23, " # A_8 A_8 Z coupling " 
End If 

 
If (WriteHiggsDiphotonLoopContributions) Then 
Write(io_L,100) "Block HPPloops # Loop contributions to H-Photon-Photon coupling " 
Do i1=1,8
CurrentPDG2(1) = Abs(PDGhh(i1)) 
Do i2=1,1
CurrentPDG2(2) = Abs(PDGVWm) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopVWm(i1), " # h(",i1,")-VWm(",i2,")-loop " 
End do 
Do i2=1,6
CurrentPDG2(2) = Abs(PDGSd(i2)) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopSd(i2,i1), " # h(",i1,")-Sd(",i2,")-loop " 
End do 
Do i2=1,6
CurrentPDG2(2) = Abs(PDGSu(i2)) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopSu(i2,i1), " # h(",i1,")-Su(",i2,")-loop " 
End do 
Do i2=2,8
CurrentPDG2(2) = Abs(PDGHpm(i2)) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopHpm(i2,i1), " # h(",i1,")-Hpm(",i2,")-loop " 
End do 
Do i2=1,5
CurrentPDG2(2) = Abs(PDGCha(i2)) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopCha(i2,i1), " # h(",i1,")-Cha(",i2,")-loop " 
End do 
Do i2=1,3
CurrentPDG2(2) = Abs(PDGFd(i2)) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopFd(i2,i1), " # h(",i1,")-Fd(",i2,")-loop " 
End do 
Do i2=1,3
CurrentPDG2(2) = Abs(PDGFu(i2)) 
Write(io_L,122) CurrentPDG2(1), CurrentPDG2(2), HPPloopFu(i2,i1), " # h(",i1,")-Fu(",i2,")-loop " 
End do 
End Do 
End if 

 
Write(io_L,100) "Block EFFHIGGSCOUPLINGS # values of loop-induced couplings " 
facPP = Alpha*Sqrt(2._dp*G_F/sqrt(2._dp))/(2._dp*Pi) 
facGG = AlphaS_MZ*Sqrt(2._dp*G_F/sqrt(2._dp))/(Sqrt(2._dp)*2._dp*Pi)
facPZ = 0._dp 
Do i1=1,8
CurrentPDG3(1) = Abs(PDGhh(i1)) 
CurrentPDG3(2) = Abs(PDGVP) 
CurrentPDG3(3) = Abs(PDGVP) 
Write(io_L,121) CurrentPDG3(1), CurrentPDG3(2), CurrentPDG3(3), Abs(CoupHPP(i1))*facPP, " # H-Photon-Photon " 
CurrentPDG3(2) = Abs(PDGVG) 
CurrentPDG3(3) = Abs(PDGVG) 
Write(io_L,121) CurrentPDG3(1), CurrentPDG3(2), CurrentPDG3(3), Abs(CoupHGG(i1))*facGG, " # H-Gluon-Gluon " 
CurrentPDG3(2) = Abs(PDGVP) 
CurrentPDG3(3) = Abs(PDGVZ) 
Write(io_L,121) CurrentPDG3(1), CurrentPDG3(2), CurrentPDG3(3), 0._dp, " # H-Photon-Z (not yet calculated by SPheno) " 
End Do 
Do i1=2,8
CurrentPDG3(1) = Abs(PDGAh(i1)) 
CurrentPDG3(2) = Abs(PDGVP) 
CurrentPDG3(3) = Abs(PDGVP) 
Write(io_L,121) CurrentPDG3(1), CurrentPDG3(2), CurrentPDG3(3), Abs(CoupAPP(i1))*facPP, " # A-Photon-Photon " 
CurrentPDG3(2) = Abs(PDGVG) 
CurrentPDG3(3) = Abs(PDGVG) 
Write(io_L,121) CurrentPDG3(1), CurrentPDG3(2), CurrentPDG3(3), Abs(CoupAGG(i1))*facGG, " # A-Gluon-Gluon " 
CurrentPDG3(2) = Abs(PDGVP) 
CurrentPDG3(3) = Abs(PDGVZ) 
Write(io_L,121) CurrentPDG3(1), CurrentPDG3(2), CurrentPDG3(3), 0._dp, " # A-Photon-Z (not yet calculated by SPheno) " 
End Do 
End If 

 
Write(io_L,100) "Block SPhenoLowEnergy # low energy observables " 
Write(io_L,1010) 20,ae,  "# (g-2)_e" 
Write(io_L,1010) 21,amu,  "# (g-2)_mu" 
Write(io_L,1010) 22,atau,  "# (g-2)_tau" 
Write(io_L,1010) 23,EDMe,  "# EDM(e)" 
Write(io_L,1010) 24,EDMmu,  "# EDM(mu)" 
Write(io_L,1010) 25,EDMtau,  "# EDM(tau)" 
Write(io_L,1010) 39,dRho,  "# delta(rho)" 
Write(io_L,100) "Block FlavorKitQFV # quark flavor violating observables " 
Write(io_L,1010) 200,BrBsGamma,  "# BR(B->X_s gamma)" 
Write(io_L,1010) 201,ratioBsGamma,  "# BR(B->X_s gamma)/BR(B->X_s gamma)_SM" 
Write(io_L,1010) 300,BrDmunu,  "# BR(D->mu nu)" 
Write(io_L,1010) 301,ratioDmunu,  "# BR(D->mu nu)/BR(D->mu nu)_SM" 
Write(io_L,1010) 400,BrDsmunu,  "# BR(Ds->mu nu)" 
Write(io_L,1010) 401,ratioDsmunu,  "# BR(Ds->mu nu)/BR(Ds->mu nu)_SM" 
Write(io_L,1010) 402,BrDstaunu,  "# BR(Ds->tau nu)" 
Write(io_L,1010) 403,ratioDstaunu,  "# BR(Ds->tau nu)/BR(Ds->tau nu)_SM" 
Write(io_L,1010) 500,BrBmunu,  "# BR(B->mu nu)" 
Write(io_L,1010) 501,ratioBmunu,  "# BR(B->mu nu)/BR(B->mu nu)_SM" 
Write(io_L,1010) 502,BrBtaunu,  "# BR(B->tau nu)" 
Write(io_L,1010) 503,ratioBtaunu,  "# BR(B->tau nu)/BR(B->tau nu)_SM" 
Write(io_L,1010) 600,BrKmunu,  "# BR(K->mu nu)" 
Write(io_L,1010) 601,ratioKmunu,  "# BR(K->mu nu)/BR(K->mu nu)_SM" 
Write(io_L,1010) 602,RK,  "# R_K = BR(K->e nu)/(K->mu nu)" 
Write(io_L,1010) 603,RKSM,  "# R_K^SM = BR(K->e nu)_SM/(K->mu nu)_SM" 
Write(io_L,1010) 1900,DeltaMBs,  "# Delta(M_Bs)" 
Write(io_L,1010) 1901,ratioDeltaMBs,  "# Delta(M_Bs)/Delta(M_Bs)_SM" 
Write(io_L,1010) 1902,DeltaMBq,  "# Delta(M_Bd)" 
Write(io_L,1010) 1903,ratioDeltaMBq,  "# Delta(M_Bd)/Delta(M_Bd)_SM" 
Write(io_L,1010) 4000,BrB0dEE,  "# BR(B^0_d->e e)" 
Write(io_L,1010) 4001,ratioB0dEE,  "# BR(B^0_d->e e)/BR(B^0_d->e e)_SM" 
Write(io_L,1010) 4002,BrB0sEE,  "# BR(B^0_s->e e)" 
Write(io_L,1010) 4003,ratioB0sEE,  "# BR(B^0_s->e e)/BR(B^0_s->e e)_SM" 
Write(io_L,1010) 4004,BrB0dMuMu,  "# BR(B^0_d->mu mu)" 
Write(io_L,1010) 4005,ratioB0dMuMu,  "# BR(B^0_d->mu mu)/BR(B^0_d->mu mu)_SM" 
Write(io_L,1010) 4006,BrB0sMuMu,  "# BR(B^0_s->mu mu)" 
Write(io_L,1010) 4007,ratioB0sMuMu,  "# BR(B^0_s->mu mu)/BR(B^0_s->mu mu)_SM" 
Write(io_L,1010) 4008,BrB0dTauTau,  "# BR(B^0_d->tau tau)" 
Write(io_L,1010) 4009,ratioB0dTauTau,  "# BR(B^0_d->tau tau)/BR(B^0_d->tau tau)_SM" 
Write(io_L,1010) 4010,BrB0sTauTau,  "# BR(B^0_s->tau tau)" 
Write(io_L,1010) 4011,ratioB0sTauTau,  "# BR(B^0_s->tau tau)/BR(B^0_s->tau tau)_SM" 
Write(io_L,1010) 5000,BrBtoSEE,  "# BR(B-> s e e)" 
Write(io_L,1010) 5001,ratioBtoSEE,  "# BR(B-> s e e)/BR(B-> s e e)_SM" 
Write(io_L,1010) 5002,BrBtoSMuMu,  "# BR(B-> s mu mu)" 
Write(io_L,1010) 5003,ratioBtoSMuMu,  "# BR(B-> s mu mu)/BR(B-> s mu mu)_SM" 
Write(io_L,1010) 6000,BrBtoKmumu,  "# BR(B -> K mu mu)" 
Write(io_L,1010) 6001,ratioBtoKmumu,  "# BR(B -> K mu mu)/BR(B -> K mu mu)_SM" 
Write(io_L,1010) 7000,BrBtoSnunu,  "# BR(B->s nu nu)" 
Write(io_L,1010) 7001,ratioBtoSnunu,  "# BR(B->s nu nu)/BR(B->s nu nu)_SM" 
Write(io_L,1010) 7002,BrBtoDnunu,  "# BR(B->D nu nu)" 
Write(io_L,1010) 7003,ratioBtoDnunu,  "# BR(B->D nu nu)/BR(B->D nu nu)_SM" 
Write(io_L,1010) 8000,BrKptoPipnunu,  "# BR(K^+ -> pi^+ nu nu)" 
Write(io_L,1010) 8001,ratioKptoPipnunu,  "# BR(K^+ -> pi^+ nu nu)/BR(K^+ -> pi^+ nu nu)_SM" 
Write(io_L,1010) 8002,BrKltoPinunu,  "# BR(K_L -> pi^0 nu nu)" 
Write(io_L,1010) 8003,ratioKltoPinunu,  "# BR(K_L -> pi^0 nu nu)/BR(K_L -> pi^0 nu nu)_SM" 
Write(io_L,1010) 9100,DelMK,  "# Delta(M_K)" 
Write(io_L,1010) 9102,ratioDelMK,  "# Delta(M_K)/Delta(M_K)_SM" 
Write(io_L,1010) 9103,epsK,  "# epsilon_K" 
Write(io_L,1010) 9104,ratioepsK,  "# epsilon_K/epsilon_K^SM" 
Write(io_L,100) "Block FlavorKitLFV # lepton flavor violating observables " 
Write(io_L,1010) 701,muEgamma,  "# BR(mu->e gamma)" 
Write(io_L,1010) 702,tauEgamma,  "# BR(tau->e gamma)" 
Write(io_L,1010) 703,tauMuGamma,  "# BR(tau->mu gamma)" 
Write(io_L,1010) 800,CRmuEAl,  "# CR(mu-e, Al)" 
Write(io_L,1010) 801,CRmuETi,  "# CR(mu-e, Ti)" 
Write(io_L,1010) 802,CRmuESr,  "# CR(mu-e, Sr)" 
Write(io_L,1010) 803,CRmuESb,  "# CR(mu-e, Sb)" 
Write(io_L,1010) 804,CRmuEAu,  "# CR(mu-e, Au)" 
Write(io_L,1010) 805,CRmuEPb,  "# CR(mu-e, Pb)" 
Write(io_L,1010) 901,BRmuTo3e,  "# BR(mu->3e)" 
Write(io_L,1010) 902,BRtauTo3e,  "# BR(tau->3e)" 
Write(io_L,1010) 903,BRtauTo3mu,  "# BR(tau->3mu)" 
Write(io_L,1010) 904,BRtauToemumu,  "# BR(tau- -> e- mu+ mu-)" 
Write(io_L,1010) 905,BRtauTomuee,  "# BR(tau- -> mu- e+ e-)" 
Write(io_L,1010) 906,BRtauToemumu2,  "# BR(tau- -> e+ mu- mu-)" 
Write(io_L,1010) 907,BRtauTomuee2,  "# BR(tau- -> mu+ e- e-)" 
Write(io_L,1010) 1001,BrZtoMuE,  "# BR(Z->e mu)" 
Write(io_L,1010) 1002,BrZtoTauE,  "# BR(Z->e tau)" 
Write(io_L,1010) 1003,BrZtoTauMu,  "# BR(Z->mu tau)" 
Write(io_L,1010) 1101,BrhtoMuE,  "# BR(h->e mu)" 
Write(io_L,1010) 1102,BrhtoTauE,  "# BR(h->e tau)" 
Write(io_L,1010) 1103,BrhtoTauMu,  "# BR(h->mu tau)" 
Write(io_L,1010) 2001,BrTautoEPi,  "# BR(tau->e pi)" 
Write(io_L,1010) 2002,BrTautoEEta,  "# BR(tau->e eta)" 
Write(io_L,1010) 2003,BrTautoEEtap,  "# BR(tau->e eta')" 
Write(io_L,1010) 2004,BrTautoMuPi,  "# BR(tau->mu pi)" 
Write(io_L,1010) 2005,BrTautoMuEta,  "# BR(tau->mu eta)" 
Write(io_L,1010) 2006,BrTautoMuEtap,  "# BR(tau->mu eta')" 

 
Write(io_L,100) "Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q " 
Write(io_L,222) "     0305" , "4422" , "00", "0", Real(coeffC7sm,dp),  " # coeffC7sm"  
Write(io_L,222) "     0305" , "4422" , "00", "2", Real(coeffC7,dp),  " # coeffC7"  
Write(io_L,222) "     0305" , "4322" , "00", "2", Real(coeffC7p,dp),  " # coeffC7p"  
Write(io_L,222) "     0305" , "6421" , "00", "0", Real(coeffC8sm,dp),  " # coeffC8sm"  
Write(io_L,222) "     0305" , "6421" , "00", "2", Real(coeffC8,dp),  " # coeffC8"  
Write(io_L,222) "     0305" , "6321" , "00", "2", Real(coeffC8p,dp),  " # coeffC8p"  
Write(io_L,222) " 03051111" , "4133" , "00", "0", Real(coeffC9eeSM,dp),  " # coeffC9eeSM"  
Write(io_L,222) " 03051111" , "4133" , "00", "2", Real(coeffC9ee,dp),  " # coeffC9ee"  
Write(io_L,222) " 03051111" , "4233" , "00", "2", Real(coeffC9Pee,dp),  " # coeffC9Pee"  
Write(io_L,222) " 03051111" , "4137" , "00", "0", Real(coeffC10eeSM,dp),  " # coeffC10eeSM"  
Write(io_L,222) " 03051111" , "4137" , "00", "2", Real(coeffC10ee,dp),  " # coeffC10ee"  
Write(io_L,222) " 03051111" , "4237" , "00", "2", Real(coeffC10Pee,dp),  " # coeffC10Pee"  
Write(io_L,222) " 03051313" , "4133" , "00", "0", Real(coeffC9mumuSM,dp),  " # coeffC9mumuSM"  
Write(io_L,222) " 03051313" , "4133" , "00", "2", Real(coeffC9mumu,dp),  " # coeffC9mumu"  
Write(io_L,222) " 03051313" , "4233" , "00", "2", Real(coeffC9Pmumu,dp),  " # coeffC9Pmumu"  
Write(io_L,222) " 03051313" , "4137" , "00", "0", Real(coeffC10mumuSM,dp),  " # coeffC10mumuSM"  
Write(io_L,222) " 03051313" , "4137" , "00", "2", Real(coeffC10mumu,dp),  " # coeffC10mumu"  
Write(io_L,222) " 03051313" , "4237" , "00", "2", Real(coeffC10Pmumu,dp),  " # coeffC10Pmumu"  
Write(io_L,222) " 03051212" , "4137" , "00", "0", Real(coeffC11nu1nu1SM,dp),  " # coeffC11nu1nu1SM"  
Write(io_L,222) " 03051212" , "4137" , "00", "2", Real(coeffC11nu1nu1,dp),  " # coeffC11nu1nu1"  
Write(io_L,222) " 03051212" , "4237" , "00", "2", Real(coeffC11Pnu1nu1,dp),  " # coeffC11Pnu1nu1"  
Write(io_L,222) " 03051414" , "4137" , "00", "0", Real(coeffC11nu2nu2SM,dp),  " # coeffC11nu2nu2SM"  
Write(io_L,222) " 03051414" , "4137" , "00", "2", Real(coeffC11nu2nu2,dp),  " # coeffC11nu2nu2"  
Write(io_L,222) " 03051414" , "4237" , "00", "2", Real(coeffC11Pnu2nu2,dp),  " # coeffC11Pnu2nu2"  
Write(io_L,222) " 03051616" , "4137" , "00", "0", Real(coeffC11nu3nu3SM,dp),  " # coeffC11nu3nu3SM"  
Write(io_L,222) " 03051616" , "4137" , "00", "2", Real(coeffC11nu3nu3,dp),  " # coeffC11nu3nu3"  
Write(io_L,222) " 03051616" , "4237" , "00", "2", Real(coeffC11Pnu3nu3,dp),  " # coeffC11Pnu3nu3"  
Write(io_L,100) "Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q " 
Write(io_L,222) "     0305" , "4422" , "00", "0", Aimag(coeffC7sm),  " # coeffC7sm"  
Write(io_L,222) "     0305" , "4422" , "00", "2", Aimag(coeffC7),  " # coeffC7"  
Write(io_L,222) "     0305" , "4322" , "00", "2", Aimag(coeffC7p),  " # coeffC7p"  
Write(io_L,222) "     0305" , "6421" , "00", "0", Aimag(coeffC8sm),  " # coeffC8sm"  
Write(io_L,222) "     0305" , "6421" , "00", "2", Aimag(coeffC8),  " # coeffC8"  
Write(io_L,222) "     0305" , "6321" , "00", "2", Aimag(coeffC8p),  " # coeffC8p"  
Write(io_L,222) " 03051111" , "4133" , "00", "0", Aimag(coeffC9eeSM),  " # coeffC9eeSM"  
Write(io_L,222) " 03051111" , "4133" , "00", "2", Aimag(coeffC9ee),  " # coeffC9ee"  
Write(io_L,222) " 03051111" , "4233" , "00", "2", Aimag(coeffC9Pee),  " # coeffC9Pee"  
Write(io_L,222) " 03051111" , "4137" , "00", "0", Aimag(coeffC10eeSM),  " # coeffC10eeSM"  
Write(io_L,222) " 03051111" , "4137" , "00", "2", Aimag(coeffC10ee),  " # coeffC10ee"  
Write(io_L,222) " 03051111" , "4237" , "00", "2", Aimag(coeffC10Pee),  " # coeffC10Pee"  
Write(io_L,222) " 03051313" , "4133" , "00", "0", Aimag(coeffC9mumuSM),  " # coeffC9mumuSM"  
Write(io_L,222) " 03051313" , "4133" , "00", "2", Aimag(coeffC9mumu),  " # coeffC9mumu"  
Write(io_L,222) " 03051313" , "4233" , "00", "2", Aimag(coeffC9Pmumu),  " # coeffC9Pmumu"  
Write(io_L,222) " 03051313" , "4137" , "00", "0", Aimag(coeffC10mumuSM),  " # coeffC10mumuSM"  
Write(io_L,222) " 03051313" , "4137" , "00", "2", Aimag(coeffC10mumu),  " # coeffC10mumu"  
Write(io_L,222) " 03051313" , "4237" , "00", "2", Aimag(coeffC10Pmumu),  " # coeffC10Pmumu"  
Write(io_L,222) " 03051212" , "4137" , "00", "0", Aimag(coeffC11nu1nu1SM),  " # coeffC11nu1nu1SM"  
Write(io_L,222) " 03051212" , "4137" , "00", "2", Aimag(coeffC11nu1nu1),  " # coeffC11nu1nu1"  
Write(io_L,222) " 03051212" , "4237" , "00", "2", Aimag(coeffC11Pnu1nu1),  " # coeffC11Pnu1nu1"  
Write(io_L,222) " 03051414" , "4137" , "00", "0", Aimag(coeffC11nu2nu2SM),  " # coeffC11nu2nu2SM"  
Write(io_L,222) " 03051414" , "4137" , "00", "2", Aimag(coeffC11nu2nu2),  " # coeffC11nu2nu2"  
Write(io_L,222) " 03051414" , "4237" , "00", "2", Aimag(coeffC11Pnu2nu2),  " # coeffC11Pnu2nu2"  
Write(io_L,222) " 03051616" , "4137" , "00", "0", Aimag(coeffC11nu3nu3SM),  " # coeffC11nu3nu3SM"  
Write(io_L,222) " 03051616" , "4137" , "00", "2", Aimag(coeffC11nu3nu3),  " # coeffC11nu3nu3"  
Write(io_L,222) " 03051616" , "4237" , "00", "2", Aimag(coeffC11Pnu3nu3),  " # coeffC11Pnu3nu3"  

 
 !-------------------------------
!Sd
!-------------------------------
 
If(gTSd(1).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSd(1)),gTSd(1),Trim(NameParticleSd(1)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 5
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRSd(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSd(1,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,10
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGHpm(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleHpm(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
    Do gt3=1,6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSd(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(1))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
End if 
If(gTSd(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSd(2)),gTSd(2),Trim(NameParticleSd(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 5
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRSd(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSd(2,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,10
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGHpm(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleHpm(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
    Do gt3=1,6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSd(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(2))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
End if 
If(gTSd(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSd(3)),gTSd(3),Trim(NameParticleSd(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 5
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRSd(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSd(3,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,10
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGHpm(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleHpm(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
    Do gt3=1,6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSd(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(3))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
End if 
If(gTSd(4).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSd(4)),gTSd(4),Trim(NameParticleSd(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 5
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRSd(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSd(4,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,10
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGHpm(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleHpm(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
    Do gt3=1,6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSd(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(4))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
End if 
If(gTSd(5).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSd(5)),gTSd(5),Trim(NameParticleSd(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 5
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRSd(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSd(5,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,10
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGHpm(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleHpm(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
    Do gt3=1,6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSd(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(5))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
End if 
If(gTSd(6).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSd(6)),gTSd(6),Trim(NameParticleSd(6)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 5
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRSd(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSd(6,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,10
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGHpm(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleHpm(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
    Do gt3=1,6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGHpm(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSd(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSd(gt3) 
Write(io_L,202) BRSd(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSd(6))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
End if 

 
!-------------------------------
!Su
!-------------------------------
 
If(gTSu(1).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSu(1)),gTSu(1),Trim(NameParticleSu(1)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSu(1,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=1,5
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGHpm(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleHpm(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
  Do gt2=1,6
    Do gt3=1,5
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGSd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(1))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 
If(gTSu(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSu(2)),gTSu(2),Trim(NameParticleSu(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSu(2,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=1,5
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGHpm(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleHpm(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
  Do gt2=1,6
    Do gt3=1,5
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGSd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSu(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(2))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 
If(gTSu(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSu(3)),gTSu(3),Trim(NameParticleSu(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSu(3,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=1,5
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGHpm(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleHpm(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
  Do gt2=1,6
    Do gt3=1,5
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGSd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSu(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(3))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 
If(gTSu(4).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSu(4)),gTSu(4),Trim(NameParticleSu(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSu(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSu(4,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=1,5
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGHpm(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleHpm(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
  Do gt2=1,6
    Do gt3=1,5
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGSd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSu(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(4))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 
If(gTSu(5).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSu(5)),gTSu(5),Trim(NameParticleSu(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSu(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSu(5,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=1,5
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGHpm(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleHpm(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
  Do gt2=1,6
    Do gt3=1,5
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGSd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSu(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(5))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 
If(gTSu(6).gt.MinWidth) Then 
Write(io_L,200) INT(PDGSu(6)),gTSu(6),Trim(NameParticleSu(6)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 6
  Do gt2=2, 8
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 10
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
  Do gt2=1, 3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 8
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGSu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRSu(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
If (Maxval(BRSu(6,193:2880)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGAh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=gt2,10
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,10
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
    Do gt3=1,6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,6
  Do gt2=1,10
    Do gt3=1,5
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSd(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
    Do gt3=1,5
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,5
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=gt1,3
    Do gt3=1,6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGSu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleSu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=2,8
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGHpm(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleHpm(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSd(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,3
    Do gt3=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=2,8
  Do gt2=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGHpm(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,8
  Do gt2=1,3
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGhh(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
    Do gt3=1,6
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGGlu 
CurrentPDG3(3) = PDGSu(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt3))//" "//")"
End if 
icount = icount +1 
End Do 
 
  Do gt2=1,6
    Do gt3=1,5
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGGlu 
CurrentPDG3(2) = PDGSd(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSd(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
End Do 
 
End Do 
 
Do gt1=1,6
  Do gt2=1,5
    Do gt3=1,5
If (BRSu(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGSu(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRSu(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleSu(6))//" -> "//Trim(NameParticleSu(gt1))//" "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 

 
 !-------------------------------
!hh
!-------------------------------
 
If(gThh(1).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(1)),gThh(1),Trim(NameParticlehh(1)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(1) = BRinvH(1)+BRhh(1,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(1))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(2)),gThh(2),Trim(NameParticlehh(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(2) = BRinvH(2)+BRhh(2,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(2,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(2,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(2))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(3)),gThh(3),Trim(NameParticlehh(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(3) = BRinvH(3)+BRhh(3,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(3,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(3,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(3))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(4).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(4)),gThh(4),Trim(NameParticlehh(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(4) = BRinvH(4)+BRhh(4,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(4,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(4,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(4))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(5).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(5)),gThh(5),Trim(NameParticlehh(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(5) = BRinvH(5)+BRhh(5,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(5,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(5,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(5))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(6).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(6)),gThh(6),Trim(NameParticlehh(6)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(6) = BRinvH(6)+BRhh(6,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(6,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(6,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(6))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(7).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(7)),gThh(7),Trim(NameParticlehh(7)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(7) = BRinvH(7)+BRhh(7,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(7,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(7,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(7))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gThh(8).gt.MinWidth) Then 
Write(io_L,200) INT(PDGhh(8)),gThh(8),Trim(NameParticlehh(8)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGVWm 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleVWm)//"^* "//Trim(NameParticleVWm)//"_virt "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvH(8) = BRinvH(8)+BRhh(8,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(8,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRhh(8,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
If (BRhh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVZ 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRhh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticlehh(8))//" -> "//Trim(NameParticleVZ)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 

 
!-------------------------------
!Ah
!-------------------------------
 
If(gTAh(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(2)),gTAh(2),Trim(NameParticleAh(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(2) = BRinvA(2)+BRAh(2,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(2,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(2,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(2))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTAh(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(3)),gTAh(3),Trim(NameParticleAh(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(3) = BRinvA(3)+BRAh(3,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(3,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(3,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(3))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTAh(4).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(4)),gTAh(4),Trim(NameParticleAh(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(4) = BRinvA(4)+BRAh(4,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(4,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(4,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(4))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTAh(5).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(5)),gTAh(5),Trim(NameParticleAh(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(5) = BRinvA(5)+BRAh(5,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(5,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(5,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(5))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTAh(6).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(6)),gTAh(6),Trim(NameParticleAh(6)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(6) = BRinvA(6)+BRAh(6,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(6,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(6,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(6))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTAh(7).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(7)),gTAh(7),Trim(NameParticleAh(7)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(7) = BRinvA(7)+BRAh(7,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(7,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(7,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(7))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTAh(8).gt.MinWidth) Then 
Write(io_L,200) INT(PDGAh(8)),gTAh(8),Trim(NameParticleAh(8)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVP 
CurrentPDG2(2) = PDGVP 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleVP)//" "//Trim(NameParticleVP)//" "//")"
End if 
icount = icount +1 
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVG 
CurrentPDG2(2) = PDGVG 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleVG)//" "//Trim(NameParticleVG)//" "//")"
End if 
icount = icount +1 
Do gt1= 2, 8
  Do gt2= gt1, 8
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2=2, 8
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 5
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2= gt1, 10
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGChi(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//")"
If (((gTChi(gt1).le.WidthToBeInvisible).and.(gTChi(gt2).le.WidthToBeInvisible)).OR. & 
 & ((WidthToBeInvisible.le.-1._dp).And.(CurrentPDG2(1).eq.PDGLSP(1)).And.(CurrentPDG2(2).eq.PDGLSP(1)))) Then 
  BRinvA(8) = BRinvA(8)+BRAh(8,icount) 
End if 
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFu(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
  Do gt2= gt1, 8
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGHpm(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRAh(8,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRAh(8,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleHpm(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSd(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleSd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRAh(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRAh(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleAh(8))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 

 
 !-------------------------------
!Hpm
!-------------------------------
 
If(gTHpm(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(2)),gTHpm(2),Trim(NameParticleHpm(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(2))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gTHpm(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(3)),gTHpm(3),Trim(NameParticleHpm(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(3))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gTHpm(4).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(4)),gTHpm(4),Trim(NameParticleHpm(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(4))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gTHpm(5).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(5)),gTHpm(5),Trim(NameParticleHpm(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(5))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gTHpm(6).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(6)),gTHpm(6),Trim(NameParticleHpm(6)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(6))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gTHpm(7).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(7)),gTHpm(7),Trim(NameParticleHpm(7)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(7))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 
If(gTHpm(8).gt.MinWidth) Then 
Write(io_L,200) INT(PDGHpm(8)),gTHpm(8),Trim(NameParticleHpm(8)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 2, 8
  Do gt2=2, 8
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGAh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleAh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 5
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGCha(gt2) 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 3
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGFd(gt2) 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 2, 8
  Do gt2=1, 8
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 8
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGhh(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticlehh(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 2, 8
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGHpm(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleHpm(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 6
  Do gt2=1, 6
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGSu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleSu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (BRHpm(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGVWm 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRHpm(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleHpm(8))//" -> "//Trim(NameParticleVWm)//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
End if 

 
 !-------------------------------
!Glu
!-------------------------------
 
If(gTGlu.gt.MinWidth) Then 
Write(io_L,200) INT(PDGGlu),gTGlu,Trim(NameParticleGlu) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRGlu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRGlu(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRGlu(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRGlu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRGlu(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRGlu(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRGlu(1,37:261)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,10
If (BRGlu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRGlu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,5
If (BRGlu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRGlu(1,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleCha(gt3))//"^* "//")"
Write(io_L,202) BRGlu(1,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
    Do gt3=1,10
If (BRGlu(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRGlu(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleGlu)//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
End if 

 
 !-------------------------------
!Fu
!-------------------------------
 
If(gTFu(1).gt.MinWidth) Then 
Write(io_L,200) INT(PDGFu(1)),gTFu(1),Trim(NameParticleFu(1)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 3
  Do gt2=2, 8
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2=1, 6
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=2, 8
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 8
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
  Do gt2=1, 6
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 6
If (BRFu(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRFu(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTFu(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGFu(2)),gTFu(2),Trim(NameParticleFu(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 3
  Do gt2=2, 8
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2=1, 6
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=2, 8
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 8
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
  Do gt2=1, 6
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 6
If (BRFu(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRFu(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTFu(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGFu(3)),gTFu(3),Trim(NameParticleFu(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 3
  Do gt2=2, 8
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
  Do gt2=1, 6
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=2, 8
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 8
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
  Do gt2=1, 6
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGGlu 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleGlu)//" "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 6
If (BRFu(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRFu(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleFu(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 

 
 !-------------------------------
!Cha
!-------------------------------
 
If(gTCha(1).gt.MinWidth) Then 
Write(io_L,200) INT(PDGCha(1)),gTCha(1),Trim(NameParticleCha(1)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 8
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRCha(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRCha(1,197:1010)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,5
  Do gt2=1,5
    Do gt3=gt1,5
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRCha(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,10
    Do gt3=gt2,10
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRCha(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRCha(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRCha(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRCha(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRCha(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRCha(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTCha(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGCha(2)),gTCha(2),Trim(NameParticleCha(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 8
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRCha(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRCha(2,197:1010)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,5
  Do gt2=1,5
    Do gt3=gt1,5
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRCha(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,10
    Do gt3=gt2,10
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRCha(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRCha(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRCha(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRCha(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRCha(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRCha(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTCha(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGCha(3)),gTCha(3),Trim(NameParticleCha(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 8
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGHpm(gt2) 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVWm 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFu(gt1) 
CurrentPDG2(2) = PDGSd(gt2) 
Write(io_L,201) BRCha(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRCha(3,197:1010)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,5
  Do gt2=1,5
    Do gt3=gt1,5
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRCha(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,10
    Do gt3=gt2,10
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRCha(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRCha(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRCha(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRCha(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRCha(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRCha(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTCha(4).gt.MinWidth) Then 
Write(io_L,200) -INT(PDGCha(4)),gTCha(4),Trim(NameParticleCha(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 8
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRCha(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRCha(4,197:1010)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,5
  Do gt2=1,5
    Do gt3=gt1,5
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRCha(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,10
    Do gt3=gt2,10
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRCha(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRCha(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRCha(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRCha(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRCha(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRCha(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(4))//"^* -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTCha(5).gt.MinWidth) Then 
Write(io_L,200) -INT(PDGCha(5)),gTCha(5),Trim(NameParticleCha(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=1, 8
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGCha(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = -PDGFd(gt1) 
CurrentPDG2(2) = PDGSu(gt2) 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRCha(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRCha(5,197:1010)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,5
  Do gt2=1,5
    Do gt3=gt1,5
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGCha(gt2) 
CurrentPDG3(3) = -PDGCha(gt3) 
Write(io_L,202) BRCha(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleCha(gt2))//" "//Trim(NameParticleCha(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,10
    Do gt3=gt2,10
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRCha(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFd(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRCha(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGCha(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFu(gt3) 
Write(io_L,202) BRCha(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = -PDGFd(gt3) 
Write(io_L,202) BRCha(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleFd(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRCha(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = -PDGFd(gt1) 
CurrentPDG3(2) = PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRCha(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleCha(5))//"^* -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleFu(gt2))//" "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 

 
 !-------------------------------
!Chi
!-------------------------------
 
If(gTChi(1).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(1)),gTChi(1),Trim(NameParticleChi(1)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(1,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(1,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(1,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(1,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(1,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(1,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(1,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(1,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(1))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(2).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(2)),gTChi(2),Trim(NameParticleChi(2)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(2,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(2,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(2,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(2,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(2,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(2,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(2,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(2,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(2,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(2,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(2,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(2,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(2,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(2,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(2))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(3).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(3)),gTChi(3),Trim(NameParticleChi(3)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(3,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(3,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(3,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(3,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(3,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(3,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(3,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(3,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(3,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(3,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(3,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(3,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(3,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(3,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(3))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(4).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(4)),gTChi(4),Trim(NameParticleChi(4)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(4,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(4,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(4,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(4,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(4,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(4,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(4,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(4,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(4,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(4,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(4,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(4,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(4,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(4,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(4))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(5).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(5)),gTChi(5),Trim(NameParticleChi(5)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(5,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(5,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(5,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(5,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(5,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(5,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(5,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(5,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(5,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(5,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(5,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(5,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(5,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(5,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(5))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(6).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(6)),gTChi(6),Trim(NameParticleChi(6)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(6,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(6,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(6,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(6,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(6,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(6,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(6,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(6,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(6,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(6,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(6,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(6,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(6,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(6,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(6))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(7).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(7)),gTChi(7),Trim(NameParticleChi(7)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(7,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(7,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(7,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(7,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(7,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(7,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(7,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(7,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(7,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(7,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(7,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(7,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(7,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(7,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(7,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(7,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(7,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(7,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(7,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(7))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(8).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(8)),gTChi(8),Trim(NameParticleChi(8)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(8,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(8,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(8,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(8,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(8,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(8,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(8,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(8,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(8,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(8,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(8,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(8,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(8,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(8,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(8,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(8,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(8,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(8,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(8,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(8))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(9).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(9)),gTChi(9),Trim(NameParticleChi(9)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(9,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(9,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(9,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(9,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(9,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(9,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(9,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(9,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(9,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(9,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(9,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(9,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(9,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(9,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(9,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(9,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(9,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(9,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(9,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(9,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(9,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(9))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
If(gTChi(10).gt.MinWidth) Then 
Write(io_L,200) INT(PDGChi(10)),gTChi(10),Trim(NameParticleChi(10)) 
Write(io_L,100) "#    BR                NDA      ID1      ID2" 
icount = 1 
Do gt1= 1, 10
  Do gt2=2, 8
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGAh(gt2) 
Write(io_L,201) BRChi(10,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleAh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
  Do gt2=2, 8
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGHpm(gt2) 
Write(io_L,201) BRChi(10,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleHpm(gt2))//"^* "//")"
Write(io_L,201) BRChi(10,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleHpm(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 5
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGCha(gt1) 
CurrentPDG2(2) = -PDGVWm 
Write(io_L,201) BRChi(10,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleVWm)//"^* "//")"
Write(io_L,201) BRChi(10,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleVWm)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 10
  Do gt2=1, 8
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGhh(gt2) 
Write(io_L,201) BRChi(10,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticlehh(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 10
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGChi(gt1) 
CurrentPDG2(2) = PDGVZ 
Write(io_L,201) BRChi(10,icount),2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleVZ)//" "//")"
End if 
icount = icount +1 
  End Do 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFd(gt1) 
CurrentPDG2(2) = -PDGSd(gt2) 
Write(io_L,201) BRChi(10,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleSd(gt2))//"^* "//")"
Write(io_L,201) BRChi(10,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleFd(gt1))//"^* "//Trim(NameParticleSd(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1= 1, 3
  Do gt2=1, 6
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG2(1) = PDGFu(gt1) 
CurrentPDG2(2) = -PDGSu(gt2) 
Write(io_L,201) BRChi(10,icount)/2._dp,2,CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleSu(gt2))//"^* "//")"
Write(io_L,201) BRChi(10,icount)/2._dp,2,-CurrentPDG2, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleFu(gt1))//"^* "//Trim(NameParticleSu(gt2))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
If (Maxval(BRChi(10,237:1729)).Gt.BRmin) Then 
Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3" 
End If 
Do gt1=1,10
  Do gt2=1,5
    Do gt3=1,5
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGCha(gt2) 
CurrentPDG3(3) = PDGCha(gt3) 
Write(io_L,202) BRChi(10,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleCha(gt2))//"^* "//Trim(NameParticleCha(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=gt1,10
    Do gt3=gt2,10
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = PDGChi(gt2) 
CurrentPDG3(3) = PDGChi(gt3) 
Write(io_L,202) BRChi(10,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleChi(gt2))//" "//Trim(NameParticleChi(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFd(gt3) 
Write(io_L,202) BRChi(10,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFd(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,10
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGChi(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(10,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleChi(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,5
  Do gt2=1,3
    Do gt3=1,3
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGCha(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGFu(gt3) 
Write(io_L,202) BRChi(10,icount)/2._dp,3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleCha(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleFu(gt3))//" "//")"
Write(io_L,202) BRChi(10,icount)/2._dp,3,-CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleCha(gt1))//"^* "//Trim(NameParticleFd(gt2))//" "//Trim(NameParticleFu(gt3))//"^* "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFd(gt1) 
CurrentPDG3(2) = -PDGFd(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(10,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleFd(gt1))//" "//Trim(NameParticleFd(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
Do gt1=1,3
  Do gt2=1,3
If (BRChi(10,icount).Gt.BrMin) Then 
CurrentPDG3(1) = PDGFu(gt1) 
CurrentPDG3(2) = -PDGFu(gt2) 
CurrentPDG3(3) = PDGGlu 
Write(io_L,202) BRChi(10,icount),3,CurrentPDG3, & 
 & Trim(NameParticleChi(10))//" -> "//Trim(NameParticleFu(gt1))//" "//Trim(NameParticleFu(gt2))//"^* "//Trim(NameParticleGlu)//" "//")"
End if 
icount = icount +1 
  End Do 
End Do 
 
End if 
99 Format(1x,i5,3x,a) 
100 Format(a) 
101 Format(2x,i3,2x,1P,e16.8,2x,a) 
1010 Format(2x,i6,2x,1P,e16.8,2x,a) 
102 Format(1x,i9,3x,1P,e16.8,2x,a) 
103 Format(a13,1P,e16.8,2x,a) 
104 Format(i4,2x,1P,e16.8,2x,a) 
105 Format(1x,2i3,3x,1P,e16.8,3x,a) 
106 Format(a,1P,e16.8,2x,a) 
107 Format(2i3,3x,1P,e16.8,3x,a) 
127 Format(3i3,3x,1P,e16.8,3x,a) 
117 Format(i3,i8,3x,1P,e16.8,3x,a) 
118 Format(i3,i10,3x,1P,a) 
119 Format(i3,i10,3x,1P,3x,e16.8,a) 
120 Format(i3,i10,3x,1P,3x,e16.8,a,i2,a,i2,a) 
121 Format(i10,3x,i10,3x,i10,3x,e16.8,a) 
122 Format(i10,i10,3x,1P,3x,e16.8,a,i2,a,i2,a) 
108 Format(9x,1P,E16.8,0P,3x,a) 
109 Format(1x,3i3,3x,1P,e16.8,3x,a) 
110 Format(3x,2i3,3x,"# ",a) 
200 Format("DECAY",1x,I9,3x,1P,E16.8,0P,3x,"# ",a) 
201 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i10,1x),2x,"# BR(",a) 
202 Format(3x,1P,e16.8,0p,3x,I2,3x,3(i10,1x),2x,"# BR(",a) 
222 Format(1x,a8,1x,a4,3x,a2,3x,a1,3x,E16.8,3x,a) 
4711 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i10,1x),2x," # ",A)
4712 Format("XS 11 -11 ",F7.1," ",F5.2," ",F5.2," ",A)

5410 Format(a25,1p,e16.7) 
5411 Format(a25,1p,"(",e16.7,",",e16.7,")") 
1101 Format(1P,2x,e16.8,2x,e16.8,0P,5x,i4,5x,3i10,a) 
1102 Format(1P,2x,e16.8,0P,5x,i4,5x,3i10,a) 
1103 Format(1P,2x,e16.8,0P,5x,i4,5x,4i10,a) 
End Subroutine LesHouches_Out 
 
 
Subroutine WriteWHIZARD 
   Open(123,file="WHIZARD.par.munuSSM3G",status="unknown")
Write(123,*) "# Couplings and VEVs" 
 
Write(123,*) "" 
Write(123,*) "Yd11_r= ",Real(Yd(1,1),dp)
Write(123,*) "Yd12_r= ",Real(Yd(1,2),dp)
Write(123,*) "Yd13_r= ",Real(Yd(1,3),dp)
Write(123,*) "Yd21_r= ",Real(Yd(2,1),dp)
Write(123,*) "Yd22_r= ",Real(Yd(2,2),dp)
Write(123,*) "Yd23_r= ",Real(Yd(2,3),dp)
Write(123,*) "Yd31_r= ",Real(Yd(3,1),dp)
Write(123,*) "Yd32_r= ",Real(Yd(3,2),dp)
Write(123,*) "Yd33_r= ",Real(Yd(3,3),dp)
Write(123,*) "Yd11_i= ",AImag(Yd(1,1))
Write(123,*) "Yd12_i= ",AImag(Yd(1,2))
Write(123,*) "Yd13_i= ",AImag(Yd(1,3))
Write(123,*) "Yd21_i= ",AImag(Yd(2,1))
Write(123,*) "Yd22_i= ",AImag(Yd(2,2))
Write(123,*) "Yd23_i= ",AImag(Yd(2,3))
Write(123,*) "Yd31_i= ",AImag(Yd(3,1))
Write(123,*) "Yd32_i= ",AImag(Yd(3,2))
Write(123,*) "Yd33_i= ",AImag(Yd(3,3))
Write(123,*) "lam1_r= ",Real(lam(1),dp)
Write(123,*) "lam2_r= ",Real(lam(2),dp)
Write(123,*) "lam3_r= ",Real(lam(3),dp)
Write(123,*) "lam1_i= ",AImag(lam(1))
Write(123,*) "lam2_i= ",AImag(lam(2))
Write(123,*) "lam3_i= ",AImag(lam(3))
Write(123,*) "Yv11_r= ",Real(Yv(1,1),dp)
Write(123,*) "Yv12_r= ",Real(Yv(1,2),dp)
Write(123,*) "Yv13_r= ",Real(Yv(1,3),dp)
Write(123,*) "Yv21_r= ",Real(Yv(2,1),dp)
Write(123,*) "Yv22_r= ",Real(Yv(2,2),dp)
Write(123,*) "Yv23_r= ",Real(Yv(2,3),dp)
Write(123,*) "Yv31_r= ",Real(Yv(3,1),dp)
Write(123,*) "Yv32_r= ",Real(Yv(3,2),dp)
Write(123,*) "Yv33_r= ",Real(Yv(3,3),dp)
Write(123,*) "Yv11_i= ",AImag(Yv(1,1))
Write(123,*) "Yv12_i= ",AImag(Yv(1,2))
Write(123,*) "Yv13_i= ",AImag(Yv(1,3))
Write(123,*) "Yv21_i= ",AImag(Yv(2,1))
Write(123,*) "Yv22_i= ",AImag(Yv(2,2))
Write(123,*) "Yv23_i= ",AImag(Yv(2,3))
Write(123,*) "Yv31_i= ",AImag(Yv(3,1))
Write(123,*) "Yv32_i= ",AImag(Yv(3,2))
Write(123,*) "Yv33_i= ",AImag(Yv(3,3))
Write(123,*) "Yu11_r= ",Real(Yu(1,1),dp)
Write(123,*) "Yu12_r= ",Real(Yu(1,2),dp)
Write(123,*) "Yu13_r= ",Real(Yu(1,3),dp)
Write(123,*) "Yu21_r= ",Real(Yu(2,1),dp)
Write(123,*) "Yu22_r= ",Real(Yu(2,2),dp)
Write(123,*) "Yu23_r= ",Real(Yu(2,3),dp)
Write(123,*) "Yu31_r= ",Real(Yu(3,1),dp)
Write(123,*) "Yu32_r= ",Real(Yu(3,2),dp)
Write(123,*) "Yu33_r= ",Real(Yu(3,3),dp)
Write(123,*) "Yu11_i= ",AImag(Yu(1,1))
Write(123,*) "Yu12_i= ",AImag(Yu(1,2))
Write(123,*) "Yu13_i= ",AImag(Yu(1,3))
Write(123,*) "Yu21_i= ",AImag(Yu(2,1))
Write(123,*) "Yu22_i= ",AImag(Yu(2,2))
Write(123,*) "Yu23_i= ",AImag(Yu(2,3))
Write(123,*) "Yu31_i= ",AImag(Yu(3,1))
Write(123,*) "Yu32_i= ",AImag(Yu(3,2))
Write(123,*) "Yu33_i= ",AImag(Yu(3,3))
Write(123,*) "kap111_r= ",Real(kap(1,1,1),dp)
Write(123,*) "kap112_r= ",Real(kap(1,1,2),dp)
Write(123,*) "kap113_r= ",Real(kap(1,1,3),dp)
Write(123,*) "kap121_r= ",Real(kap(1,2,1),dp)
Write(123,*) "kap122_r= ",Real(kap(1,2,2),dp)
Write(123,*) "kap123_r= ",Real(kap(1,2,3),dp)
Write(123,*) "kap131_r= ",Real(kap(1,3,1),dp)
Write(123,*) "kap132_r= ",Real(kap(1,3,2),dp)
Write(123,*) "kap133_r= ",Real(kap(1,3,3),dp)
Write(123,*) "kap211_r= ",Real(kap(2,1,1),dp)
Write(123,*) "kap212_r= ",Real(kap(2,1,2),dp)
Write(123,*) "kap213_r= ",Real(kap(2,1,3),dp)
Write(123,*) "kap221_r= ",Real(kap(2,2,1),dp)
Write(123,*) "kap222_r= ",Real(kap(2,2,2),dp)
Write(123,*) "kap223_r= ",Real(kap(2,2,3),dp)
Write(123,*) "kap231_r= ",Real(kap(2,3,1),dp)
Write(123,*) "kap232_r= ",Real(kap(2,3,2),dp)
Write(123,*) "kap233_r= ",Real(kap(2,3,3),dp)
Write(123,*) "kap311_r= ",Real(kap(3,1,1),dp)
Write(123,*) "kap312_r= ",Real(kap(3,1,2),dp)
Write(123,*) "kap313_r= ",Real(kap(3,1,3),dp)
Write(123,*) "kap321_r= ",Real(kap(3,2,1),dp)
Write(123,*) "kap322_r= ",Real(kap(3,2,2),dp)
Write(123,*) "kap323_r= ",Real(kap(3,2,3),dp)
Write(123,*) "kap331_r= ",Real(kap(3,3,1),dp)
Write(123,*) "kap332_r= ",Real(kap(3,3,2),dp)
Write(123,*) "kap333_r= ",Real(kap(3,3,3),dp)
Write(123,*) "kap111_i= ",AImag(kap(1,1,1))
Write(123,*) "kap112_i= ",AImag(kap(1,1,2))
Write(123,*) "kap113_i= ",AImag(kap(1,1,3))
Write(123,*) "kap121_i= ",AImag(kap(1,2,1))
Write(123,*) "kap122_i= ",AImag(kap(1,2,2))
Write(123,*) "kap123_i= ",AImag(kap(1,2,3))
Write(123,*) "kap131_i= ",AImag(kap(1,3,1))
Write(123,*) "kap132_i= ",AImag(kap(1,3,2))
Write(123,*) "kap133_i= ",AImag(kap(1,3,3))
Write(123,*) "kap211_i= ",AImag(kap(2,1,1))
Write(123,*) "kap212_i= ",AImag(kap(2,1,2))
Write(123,*) "kap213_i= ",AImag(kap(2,1,3))
Write(123,*) "kap221_i= ",AImag(kap(2,2,1))
Write(123,*) "kap222_i= ",AImag(kap(2,2,2))
Write(123,*) "kap223_i= ",AImag(kap(2,2,3))
Write(123,*) "kap231_i= ",AImag(kap(2,3,1))
Write(123,*) "kap232_i= ",AImag(kap(2,3,2))
Write(123,*) "kap233_i= ",AImag(kap(2,3,3))
Write(123,*) "kap311_i= ",AImag(kap(3,1,1))
Write(123,*) "kap312_i= ",AImag(kap(3,1,2))
Write(123,*) "kap313_i= ",AImag(kap(3,1,3))
Write(123,*) "kap321_i= ",AImag(kap(3,2,1))
Write(123,*) "kap322_i= ",AImag(kap(3,2,2))
Write(123,*) "kap323_i= ",AImag(kap(3,2,3))
Write(123,*) "kap331_i= ",AImag(kap(3,3,1))
Write(123,*) "kap332_i= ",AImag(kap(3,3,2))
Write(123,*) "kap333_i= ",AImag(kap(3,3,3))
Write(123,*) "Td11_r= ",Real(Td(1,1),dp)
Write(123,*) "Td12_r= ",Real(Td(1,2),dp)
Write(123,*) "Td13_r= ",Real(Td(1,3),dp)
Write(123,*) "Td21_r= ",Real(Td(2,1),dp)
Write(123,*) "Td22_r= ",Real(Td(2,2),dp)
Write(123,*) "Td23_r= ",Real(Td(2,3),dp)
Write(123,*) "Td31_r= ",Real(Td(3,1),dp)
Write(123,*) "Td32_r= ",Real(Td(3,2),dp)
Write(123,*) "Td33_r= ",Real(Td(3,3),dp)
Write(123,*) "Td11_i= ",AImag(Td(1,1))
Write(123,*) "Td12_i= ",AImag(Td(1,2))
Write(123,*) "Td13_i= ",AImag(Td(1,3))
Write(123,*) "Td21_i= ",AImag(Td(2,1))
Write(123,*) "Td22_i= ",AImag(Td(2,2))
Write(123,*) "Td23_i= ",AImag(Td(2,3))
Write(123,*) "Td31_i= ",AImag(Td(3,1))
Write(123,*) "Td32_i= ",AImag(Td(3,2))
Write(123,*) "Td33_i= ",AImag(Td(3,3))
Write(123,*) "Te11_r= ",Real(Te(1,1),dp)
Write(123,*) "Te12_r= ",Real(Te(1,2),dp)
Write(123,*) "Te13_r= ",Real(Te(1,3),dp)
Write(123,*) "Te21_r= ",Real(Te(2,1),dp)
Write(123,*) "Te22_r= ",Real(Te(2,2),dp)
Write(123,*) "Te23_r= ",Real(Te(2,3),dp)
Write(123,*) "Te31_r= ",Real(Te(3,1),dp)
Write(123,*) "Te32_r= ",Real(Te(3,2),dp)
Write(123,*) "Te33_r= ",Real(Te(3,3),dp)
Write(123,*) "Te11_i= ",AImag(Te(1,1))
Write(123,*) "Te12_i= ",AImag(Te(1,2))
Write(123,*) "Te13_i= ",AImag(Te(1,3))
Write(123,*) "Te21_i= ",AImag(Te(2,1))
Write(123,*) "Te22_i= ",AImag(Te(2,2))
Write(123,*) "Te23_i= ",AImag(Te(2,3))
Write(123,*) "Te31_i= ",AImag(Te(3,1))
Write(123,*) "Te32_i= ",AImag(Te(3,2))
Write(123,*) "Te33_i= ",AImag(Te(3,3))
Write(123,*) "Tlam1_r= ",Real(Tlam(1),dp)
Write(123,*) "Tlam2_r= ",Real(Tlam(2),dp)
Write(123,*) "Tlam3_r= ",Real(Tlam(3),dp)
Write(123,*) "Tlam1_i= ",AImag(Tlam(1))
Write(123,*) "Tlam2_i= ",AImag(Tlam(2))
Write(123,*) "Tlam3_i= ",AImag(Tlam(3))
Write(123,*) "Tv11_r= ",Real(Tv(1,1),dp)
Write(123,*) "Tv12_r= ",Real(Tv(1,2),dp)
Write(123,*) "Tv13_r= ",Real(Tv(1,3),dp)
Write(123,*) "Tv21_r= ",Real(Tv(2,1),dp)
Write(123,*) "Tv22_r= ",Real(Tv(2,2),dp)
Write(123,*) "Tv23_r= ",Real(Tv(2,3),dp)
Write(123,*) "Tv31_r= ",Real(Tv(3,1),dp)
Write(123,*) "Tv32_r= ",Real(Tv(3,2),dp)
Write(123,*) "Tv33_r= ",Real(Tv(3,3),dp)
Write(123,*) "Tv11_i= ",AImag(Tv(1,1))
Write(123,*) "Tv12_i= ",AImag(Tv(1,2))
Write(123,*) "Tv13_i= ",AImag(Tv(1,3))
Write(123,*) "Tv21_i= ",AImag(Tv(2,1))
Write(123,*) "Tv22_i= ",AImag(Tv(2,2))
Write(123,*) "Tv23_i= ",AImag(Tv(2,3))
Write(123,*) "Tv31_i= ",AImag(Tv(3,1))
Write(123,*) "Tv32_i= ",AImag(Tv(3,2))
Write(123,*) "Tv33_i= ",AImag(Tv(3,3))
Write(123,*) "Tu11_r= ",Real(Tu(1,1),dp)
Write(123,*) "Tu12_r= ",Real(Tu(1,2),dp)
Write(123,*) "Tu13_r= ",Real(Tu(1,3),dp)
Write(123,*) "Tu21_r= ",Real(Tu(2,1),dp)
Write(123,*) "Tu22_r= ",Real(Tu(2,2),dp)
Write(123,*) "Tu23_r= ",Real(Tu(2,3),dp)
Write(123,*) "Tu31_r= ",Real(Tu(3,1),dp)
Write(123,*) "Tu32_r= ",Real(Tu(3,2),dp)
Write(123,*) "Tu33_r= ",Real(Tu(3,3),dp)
Write(123,*) "Tu11_i= ",AImag(Tu(1,1))
Write(123,*) "Tu12_i= ",AImag(Tu(1,2))
Write(123,*) "Tu13_i= ",AImag(Tu(1,3))
Write(123,*) "Tu21_i= ",AImag(Tu(2,1))
Write(123,*) "Tu22_i= ",AImag(Tu(2,2))
Write(123,*) "Tu23_i= ",AImag(Tu(2,3))
Write(123,*) "Tu31_i= ",AImag(Tu(3,1))
Write(123,*) "Tu32_i= ",AImag(Tu(3,2))
Write(123,*) "Tu33_i= ",AImag(Tu(3,3))
Write(123,*) "Tk111_r= ",Real(Tk(1,1,1),dp)
Write(123,*) "Tk112_r= ",Real(Tk(1,1,2),dp)
Write(123,*) "Tk113_r= ",Real(Tk(1,1,3),dp)
Write(123,*) "Tk121_r= ",Real(Tk(1,2,1),dp)
Write(123,*) "Tk122_r= ",Real(Tk(1,2,2),dp)
Write(123,*) "Tk123_r= ",Real(Tk(1,2,3),dp)
Write(123,*) "Tk131_r= ",Real(Tk(1,3,1),dp)
Write(123,*) "Tk132_r= ",Real(Tk(1,3,2),dp)
Write(123,*) "Tk133_r= ",Real(Tk(1,3,3),dp)
Write(123,*) "Tk211_r= ",Real(Tk(2,1,1),dp)
Write(123,*) "Tk212_r= ",Real(Tk(2,1,2),dp)
Write(123,*) "Tk213_r= ",Real(Tk(2,1,3),dp)
Write(123,*) "Tk221_r= ",Real(Tk(2,2,1),dp)
Write(123,*) "Tk222_r= ",Real(Tk(2,2,2),dp)
Write(123,*) "Tk223_r= ",Real(Tk(2,2,3),dp)
Write(123,*) "Tk231_r= ",Real(Tk(2,3,1),dp)
Write(123,*) "Tk232_r= ",Real(Tk(2,3,2),dp)
Write(123,*) "Tk233_r= ",Real(Tk(2,3,3),dp)
Write(123,*) "Tk311_r= ",Real(Tk(3,1,1),dp)
Write(123,*) "Tk312_r= ",Real(Tk(3,1,2),dp)
Write(123,*) "Tk313_r= ",Real(Tk(3,1,3),dp)
Write(123,*) "Tk321_r= ",Real(Tk(3,2,1),dp)
Write(123,*) "Tk322_r= ",Real(Tk(3,2,2),dp)
Write(123,*) "Tk323_r= ",Real(Tk(3,2,3),dp)
Write(123,*) "Tk331_r= ",Real(Tk(3,3,1),dp)
Write(123,*) "Tk332_r= ",Real(Tk(3,3,2),dp)
Write(123,*) "Tk333_r= ",Real(Tk(3,3,3),dp)
Write(123,*) "Tk111_i= ",AImag(Tk(1,1,1))
Write(123,*) "Tk112_i= ",AImag(Tk(1,1,2))
Write(123,*) "Tk113_i= ",AImag(Tk(1,1,3))
Write(123,*) "Tk121_i= ",AImag(Tk(1,2,1))
Write(123,*) "Tk122_i= ",AImag(Tk(1,2,2))
Write(123,*) "Tk123_i= ",AImag(Tk(1,2,3))
Write(123,*) "Tk131_i= ",AImag(Tk(1,3,1))
Write(123,*) "Tk132_i= ",AImag(Tk(1,3,2))
Write(123,*) "Tk133_i= ",AImag(Tk(1,3,3))
Write(123,*) "Tk211_i= ",AImag(Tk(2,1,1))
Write(123,*) "Tk212_i= ",AImag(Tk(2,1,2))
Write(123,*) "Tk213_i= ",AImag(Tk(2,1,3))
Write(123,*) "Tk221_i= ",AImag(Tk(2,2,1))
Write(123,*) "Tk222_i= ",AImag(Tk(2,2,2))
Write(123,*) "Tk223_i= ",AImag(Tk(2,2,3))
Write(123,*) "Tk231_i= ",AImag(Tk(2,3,1))
Write(123,*) "Tk232_i= ",AImag(Tk(2,3,2))
Write(123,*) "Tk233_i= ",AImag(Tk(2,3,3))
Write(123,*) "Tk311_i= ",AImag(Tk(3,1,1))
Write(123,*) "Tk312_i= ",AImag(Tk(3,1,2))
Write(123,*) "Tk313_i= ",AImag(Tk(3,1,3))
Write(123,*) "Tk321_i= ",AImag(Tk(3,2,1))
Write(123,*) "Tk322_i= ",AImag(Tk(3,2,2))
Write(123,*) "Tk323_i= ",AImag(Tk(3,2,3))
Write(123,*) "Tk331_i= ",AImag(Tk(3,3,1))
Write(123,*) "Tk332_i= ",AImag(Tk(3,3,2))
Write(123,*) "Tk333_i= ",AImag(Tk(3,3,3))
Write(123,*) "vR1= ",vR(1)
Write(123,*) "vR2= ",vR(2)
Write(123,*) "vR3= ",vR(3)
Write(123,*) "vL1= ",vL(1)
Write(123,*) "vL2= ",vL(2)
Write(123,*) "vL3= ",vL(3)
Write(123,*) "vd= ",vd
Write(123,*) "vu= ",vu
Write(123,*) "vR1= ",vR(1)
Write(123,*) "vR2= ",vR(2)
Write(123,*) "vR3= ",vR(3)
Write(123,*) "vL1= ",vL(1)
Write(123,*) "vL2= ",vL(2)
Write(123,*) "vL3= ",vL(3)
Write(123,*) "pG_r = ",Real(pG,dp) 
Write(123,*) "pG_i = ",AImag(pG) 
Write(123,*) "" 
Write(123,*) "" 

 
 
 Write(123,*) "# Dependent parameters " 
 
Write(123,*) "" 
Write(123,*) "" 
Write(123,*) "" 

 
 
 Write(123,*) "# Necessary MINPAR parameters " 
 
Write(123,*) "" 
Write(123,*) "" 
Write(123,*) "" 

 
 
 Write(123,*) "# Masses of particles" 
 
Write(123,*) "" 
Write(123,*) "Msd1= ", Abs(MSd(1)) 
Write(123,*) "Msd2= ", Abs(MSd(2)) 
Write(123,*) "Msd3= ", Abs(MSd(3)) 
Write(123,*) "Msd4= ", Abs(MSd(4)) 
Write(123,*) "Msd5= ", Abs(MSd(5)) 
Write(123,*) "Msd6= ", Abs(MSd(6)) 
Write(123,*) "Msu1= ", Abs(MSu(1)) 
Write(123,*) "Msu2= ", Abs(MSu(2)) 
Write(123,*) "Msu3= ", Abs(MSu(3)) 
Write(123,*) "Msu4= ", Abs(MSu(4)) 
Write(123,*) "Msu5= ", Abs(MSu(5)) 
Write(123,*) "Msu6= ", Abs(MSu(6)) 
Write(123,*) "Mh1= ", Abs(Mhh(1)) 
Write(123,*) "Mh2= ", Abs(Mhh(2)) 
Write(123,*) "Mh3= ", Abs(Mhh(3)) 
Write(123,*) "Mh4= ", Abs(Mhh(4)) 
Write(123,*) "Mh5= ", Abs(Mhh(5)) 
Write(123,*) "Mh6= ", Abs(Mhh(6)) 
Write(123,*) "Mh7= ", Abs(Mhh(7)) 
Write(123,*) "Mh8= ", Abs(Mhh(8)) 
Write(123,*) "MAh2= ", Abs(MAh(2)) 
Write(123,*) "MAh3= ", Abs(MAh(3)) 
Write(123,*) "MAh4= ", Abs(MAh(4)) 
Write(123,*) "MAh5= ", Abs(MAh(5)) 
Write(123,*) "MAh6= ", Abs(MAh(6)) 
Write(123,*) "MAh7= ", Abs(MAh(7)) 
Write(123,*) "MAh8= ", Abs(MAh(8)) 
Write(123,*) "MHm2= ", Abs(MHpm(2)) 
Write(123,*) "MHm3= ", Abs(MHpm(3)) 
Write(123,*) "MHm4= ", Abs(MHpm(4)) 
Write(123,*) "MHm5= ", Abs(MHpm(5)) 
Write(123,*) "MHm6= ", Abs(MHpm(6)) 
Write(123,*) "MHm7= ", Abs(MHpm(7)) 
Write(123,*) "MHm8= ", Abs(MHpm(8)) 
Write(123,*) "MZ= ", Abs(MVZ) 
Write(123,*) "Md1= ", Abs(MFd(1)) 
Write(123,*) "Md2= ", Abs(MFd(2)) 
Write(123,*) "Md3= ", Abs(MFd(3)) 
Write(123,*) "Mu1= ", Abs(MFu(1)) 
Write(123,*) "Mu2= ", Abs(MFu(2)) 
Write(123,*) "Mu3= ", Abs(MFu(3)) 
Write(123,*) "Mgo= ", Abs(MGlu) 
Write(123,*) "Mnu4= ", Abs(MChi(4)) 
Write(123,*) "Mnu5= ", Abs(MChi(5)) 
Write(123,*) "Mnu6= ", Abs(MChi(6)) 
Write(123,*) "Mnu7= ", Abs(MChi(7)) 
Write(123,*) "Mnu8= ", Abs(MChi(8)) 
Write(123,*) "Mnu9= ", Abs(MChi(9)) 
Write(123,*) "Mnu10= ", Abs(MChi(10)) 
Write(123,*) "Me1= ", Abs(MCha(1)) 
Write(123,*) "Me2= ", Abs(MCha(2)) 
Write(123,*) "Me3= ", Abs(MCha(3)) 
Write(123,*) "Me4= ", Abs(MCha(4)) 
Write(123,*) "Me5= ", Abs(MCha(5)) 
Write(123,*) "" 
Write(123,*) "" 

 
 
 Write(123,*) "# Widths of particles" 
 
Write(123,*) "" 
Write(123,*) "Wsd1 = ",gTSd(1)
Write(123,*) "Wsd2 = ",gTSd(2)
Write(123,*) "Wsd3 = ",gTSd(3)
Write(123,*) "Wsd4 = ",gTSd(4)
Write(123,*) "Wsd5 = ",gTSd(5)
Write(123,*) "Wsd6 = ",gTSd(6)
Write(123,*) "Wsu1 = ",gTSu(1)
Write(123,*) "Wsu2 = ",gTSu(2)
Write(123,*) "Wsu3 = ",gTSu(3)
Write(123,*) "Wsu4 = ",gTSu(4)
Write(123,*) "Wsu5 = ",gTSu(5)
Write(123,*) "Wsu6 = ",gTSu(6)
Write(123,*) "Wh1 = ",gThh(1)
Write(123,*) "Wh2 = ",gThh(2)
Write(123,*) "Wh3 = ",gThh(3)
Write(123,*) "Wh4 = ",gThh(4)
Write(123,*) "Wh5 = ",gThh(5)
Write(123,*) "Wh6 = ",gThh(6)
Write(123,*) "Wh7 = ",gThh(7)
Write(123,*) "Wh8 = ",gThh(8)
Write(123,*) "WAh2 = ",gTAh(2)
Write(123,*) "WAh3 = ",gTAh(3)
Write(123,*) "WAh4 = ",gTAh(4)
Write(123,*) "WAh5 = ",gTAh(5)
Write(123,*) "WAh6 = ",gTAh(6)
Write(123,*) "WAh7 = ",gTAh(7)
Write(123,*) "WAh8 = ",gTAh(8)
Write(123,*) "WHm2 = ",gTHpm(2)
Write(123,*) "WHm3 = ",gTHpm(3)
Write(123,*) "WHm4 = ",gTHpm(4)
Write(123,*) "WHm5 = ",gTHpm(5)
Write(123,*) "WHm6 = ",gTHpm(6)
Write(123,*) "WHm7 = ",gTHpm(7)
Write(123,*) "WHm8 = ",gTHpm(8)
Write(123,*) "Wgo = ",gTGlu
Write(123,*) "Wu3 = ",gTFu(3)
Write(123,*) "We4 = ",gTCha(4)
Write(123,*) "We5 = ",gTCha(5)
Write(123,*) "Wnu4 = ",gTChi(4)
Write(123,*) "Wnu5 = ",gTChi(5)
Write(123,*) "Wnu6 = ",gTChi(6)
Write(123,*) "Wnu7 = ",gTChi(7)
Write(123,*) "Wnu8 = ",gTChi(8)
Write(123,*) "Wnu9 = ",gTChi(9)
Write(123,*) "Wnu10 = ",gTChi(10)
Write(123,*) "" 
Write(123,*) "" 

 
 
 Write(123,*) "# Mixing matrices" 
 
Write(123,*) "" 
If (MSd(1).Gt.0._dp) Then 
Write(123,*) "ZD11_r = ", Real(ZD(1,1),dp)
Write(123,*) "ZD11_i = ", AImag(ZD(1,1))
Write(123,*) "ZD12_r = ", Real(ZD(1,2),dp)
Write(123,*) "ZD12_i = ", AImag(ZD(1,2))
Write(123,*) "ZD13_r = ", Real(ZD(1,3),dp)
Write(123,*) "ZD13_i = ", AImag(ZD(1,3))
Write(123,*) "ZD14_r = ", Real(ZD(1,4),dp)
Write(123,*) "ZD14_i = ", AImag(ZD(1,4))
Write(123,*) "ZD15_r = ", Real(ZD(1,5),dp)
Write(123,*) "ZD15_i = ", AImag(ZD(1,5))
Write(123,*) "ZD16_r = ", Real(ZD(1,6),dp)
Write(123,*) "ZD16_i = ", AImag(ZD(1,6))
Else 
Write(123,*) "ZD11_i = ", Real(ZD(1,1),dp)
Write(123,*) "ZD11_r = ", -AImag(ZD(1,1))
Write(123,*) "ZD12_i = ", Real(ZD(1,2),dp)
Write(123,*) "ZD12_r = ", -AImag(ZD(1,2))
Write(123,*) "ZD13_i = ", Real(ZD(1,3),dp)
Write(123,*) "ZD13_r = ", -AImag(ZD(1,3))
Write(123,*) "ZD14_i = ", Real(ZD(1,4),dp)
Write(123,*) "ZD14_r = ", -AImag(ZD(1,4))
Write(123,*) "ZD15_i = ", Real(ZD(1,5),dp)
Write(123,*) "ZD15_r = ", -AImag(ZD(1,5))
Write(123,*) "ZD16_i = ", Real(ZD(1,6),dp)
Write(123,*) "ZD16_r = ", -AImag(ZD(1,6))
End if 
If (MSd(2).Gt.0._dp) Then 
Write(123,*) "ZD21_r = ", Real(ZD(2,1),dp)
Write(123,*) "ZD21_i = ", AImag(ZD(2,1))
Write(123,*) "ZD22_r = ", Real(ZD(2,2),dp)
Write(123,*) "ZD22_i = ", AImag(ZD(2,2))
Write(123,*) "ZD23_r = ", Real(ZD(2,3),dp)
Write(123,*) "ZD23_i = ", AImag(ZD(2,3))
Write(123,*) "ZD24_r = ", Real(ZD(2,4),dp)
Write(123,*) "ZD24_i = ", AImag(ZD(2,4))
Write(123,*) "ZD25_r = ", Real(ZD(2,5),dp)
Write(123,*) "ZD25_i = ", AImag(ZD(2,5))
Write(123,*) "ZD26_r = ", Real(ZD(2,6),dp)
Write(123,*) "ZD26_i = ", AImag(ZD(2,6))
Else 
Write(123,*) "ZD21_i = ", Real(ZD(2,1),dp)
Write(123,*) "ZD21_r = ", -AImag(ZD(2,1))
Write(123,*) "ZD22_i = ", Real(ZD(2,2),dp)
Write(123,*) "ZD22_r = ", -AImag(ZD(2,2))
Write(123,*) "ZD23_i = ", Real(ZD(2,3),dp)
Write(123,*) "ZD23_r = ", -AImag(ZD(2,3))
Write(123,*) "ZD24_i = ", Real(ZD(2,4),dp)
Write(123,*) "ZD24_r = ", -AImag(ZD(2,4))
Write(123,*) "ZD25_i = ", Real(ZD(2,5),dp)
Write(123,*) "ZD25_r = ", -AImag(ZD(2,5))
Write(123,*) "ZD26_i = ", Real(ZD(2,6),dp)
Write(123,*) "ZD26_r = ", -AImag(ZD(2,6))
End if 
If (MSd(3).Gt.0._dp) Then 
Write(123,*) "ZD31_r = ", Real(ZD(3,1),dp)
Write(123,*) "ZD31_i = ", AImag(ZD(3,1))
Write(123,*) "ZD32_r = ", Real(ZD(3,2),dp)
Write(123,*) "ZD32_i = ", AImag(ZD(3,2))
Write(123,*) "ZD33_r = ", Real(ZD(3,3),dp)
Write(123,*) "ZD33_i = ", AImag(ZD(3,3))
Write(123,*) "ZD34_r = ", Real(ZD(3,4),dp)
Write(123,*) "ZD34_i = ", AImag(ZD(3,4))
Write(123,*) "ZD35_r = ", Real(ZD(3,5),dp)
Write(123,*) "ZD35_i = ", AImag(ZD(3,5))
Write(123,*) "ZD36_r = ", Real(ZD(3,6),dp)
Write(123,*) "ZD36_i = ", AImag(ZD(3,6))
Else 
Write(123,*) "ZD31_i = ", Real(ZD(3,1),dp)
Write(123,*) "ZD31_r = ", -AImag(ZD(3,1))
Write(123,*) "ZD32_i = ", Real(ZD(3,2),dp)
Write(123,*) "ZD32_r = ", -AImag(ZD(3,2))
Write(123,*) "ZD33_i = ", Real(ZD(3,3),dp)
Write(123,*) "ZD33_r = ", -AImag(ZD(3,3))
Write(123,*) "ZD34_i = ", Real(ZD(3,4),dp)
Write(123,*) "ZD34_r = ", -AImag(ZD(3,4))
Write(123,*) "ZD35_i = ", Real(ZD(3,5),dp)
Write(123,*) "ZD35_r = ", -AImag(ZD(3,5))
Write(123,*) "ZD36_i = ", Real(ZD(3,6),dp)
Write(123,*) "ZD36_r = ", -AImag(ZD(3,6))
End if 
If (MSd(4).Gt.0._dp) Then 
Write(123,*) "ZD41_r = ", Real(ZD(4,1),dp)
Write(123,*) "ZD41_i = ", AImag(ZD(4,1))
Write(123,*) "ZD42_r = ", Real(ZD(4,2),dp)
Write(123,*) "ZD42_i = ", AImag(ZD(4,2))
Write(123,*) "ZD43_r = ", Real(ZD(4,3),dp)
Write(123,*) "ZD43_i = ", AImag(ZD(4,3))
Write(123,*) "ZD44_r = ", Real(ZD(4,4),dp)
Write(123,*) "ZD44_i = ", AImag(ZD(4,4))
Write(123,*) "ZD45_r = ", Real(ZD(4,5),dp)
Write(123,*) "ZD45_i = ", AImag(ZD(4,5))
Write(123,*) "ZD46_r = ", Real(ZD(4,6),dp)
Write(123,*) "ZD46_i = ", AImag(ZD(4,6))
Else 
Write(123,*) "ZD41_i = ", Real(ZD(4,1),dp)
Write(123,*) "ZD41_r = ", -AImag(ZD(4,1))
Write(123,*) "ZD42_i = ", Real(ZD(4,2),dp)
Write(123,*) "ZD42_r = ", -AImag(ZD(4,2))
Write(123,*) "ZD43_i = ", Real(ZD(4,3),dp)
Write(123,*) "ZD43_r = ", -AImag(ZD(4,3))
Write(123,*) "ZD44_i = ", Real(ZD(4,4),dp)
Write(123,*) "ZD44_r = ", -AImag(ZD(4,4))
Write(123,*) "ZD45_i = ", Real(ZD(4,5),dp)
Write(123,*) "ZD45_r = ", -AImag(ZD(4,5))
Write(123,*) "ZD46_i = ", Real(ZD(4,6),dp)
Write(123,*) "ZD46_r = ", -AImag(ZD(4,6))
End if 
If (MSd(5).Gt.0._dp) Then 
Write(123,*) "ZD51_r = ", Real(ZD(5,1),dp)
Write(123,*) "ZD51_i = ", AImag(ZD(5,1))
Write(123,*) "ZD52_r = ", Real(ZD(5,2),dp)
Write(123,*) "ZD52_i = ", AImag(ZD(5,2))
Write(123,*) "ZD53_r = ", Real(ZD(5,3),dp)
Write(123,*) "ZD53_i = ", AImag(ZD(5,3))
Write(123,*) "ZD54_r = ", Real(ZD(5,4),dp)
Write(123,*) "ZD54_i = ", AImag(ZD(5,4))
Write(123,*) "ZD55_r = ", Real(ZD(5,5),dp)
Write(123,*) "ZD55_i = ", AImag(ZD(5,5))
Write(123,*) "ZD56_r = ", Real(ZD(5,6),dp)
Write(123,*) "ZD56_i = ", AImag(ZD(5,6))
Else 
Write(123,*) "ZD51_i = ", Real(ZD(5,1),dp)
Write(123,*) "ZD51_r = ", -AImag(ZD(5,1))
Write(123,*) "ZD52_i = ", Real(ZD(5,2),dp)
Write(123,*) "ZD52_r = ", -AImag(ZD(5,2))
Write(123,*) "ZD53_i = ", Real(ZD(5,3),dp)
Write(123,*) "ZD53_r = ", -AImag(ZD(5,3))
Write(123,*) "ZD54_i = ", Real(ZD(5,4),dp)
Write(123,*) "ZD54_r = ", -AImag(ZD(5,4))
Write(123,*) "ZD55_i = ", Real(ZD(5,5),dp)
Write(123,*) "ZD55_r = ", -AImag(ZD(5,5))
Write(123,*) "ZD56_i = ", Real(ZD(5,6),dp)
Write(123,*) "ZD56_r = ", -AImag(ZD(5,6))
End if 
If (MSd(6).Gt.0._dp) Then 
Write(123,*) "ZD61_r = ", Real(ZD(6,1),dp)
Write(123,*) "ZD61_i = ", AImag(ZD(6,1))
Write(123,*) "ZD62_r = ", Real(ZD(6,2),dp)
Write(123,*) "ZD62_i = ", AImag(ZD(6,2))
Write(123,*) "ZD63_r = ", Real(ZD(6,3),dp)
Write(123,*) "ZD63_i = ", AImag(ZD(6,3))
Write(123,*) "ZD64_r = ", Real(ZD(6,4),dp)
Write(123,*) "ZD64_i = ", AImag(ZD(6,4))
Write(123,*) "ZD65_r = ", Real(ZD(6,5),dp)
Write(123,*) "ZD65_i = ", AImag(ZD(6,5))
Write(123,*) "ZD66_r = ", Real(ZD(6,6),dp)
Write(123,*) "ZD66_i = ", AImag(ZD(6,6))
Else 
Write(123,*) "ZD61_i = ", Real(ZD(6,1),dp)
Write(123,*) "ZD61_r = ", -AImag(ZD(6,1))
Write(123,*) "ZD62_i = ", Real(ZD(6,2),dp)
Write(123,*) "ZD62_r = ", -AImag(ZD(6,2))
Write(123,*) "ZD63_i = ", Real(ZD(6,3),dp)
Write(123,*) "ZD63_r = ", -AImag(ZD(6,3))
Write(123,*) "ZD64_i = ", Real(ZD(6,4),dp)
Write(123,*) "ZD64_r = ", -AImag(ZD(6,4))
Write(123,*) "ZD65_i = ", Real(ZD(6,5),dp)
Write(123,*) "ZD65_r = ", -AImag(ZD(6,5))
Write(123,*) "ZD66_i = ", Real(ZD(6,6),dp)
Write(123,*) "ZD66_r = ", -AImag(ZD(6,6))
End if 
If (MSu(1).Gt.0._dp) Then 
Write(123,*) "ZU11_r = ", Real(ZU(1,1),dp)
Write(123,*) "ZU11_i = ", AImag(ZU(1,1))
Write(123,*) "ZU12_r = ", Real(ZU(1,2),dp)
Write(123,*) "ZU12_i = ", AImag(ZU(1,2))
Write(123,*) "ZU13_r = ", Real(ZU(1,3),dp)
Write(123,*) "ZU13_i = ", AImag(ZU(1,3))
Write(123,*) "ZU14_r = ", Real(ZU(1,4),dp)
Write(123,*) "ZU14_i = ", AImag(ZU(1,4))
Write(123,*) "ZU15_r = ", Real(ZU(1,5),dp)
Write(123,*) "ZU15_i = ", AImag(ZU(1,5))
Write(123,*) "ZU16_r = ", Real(ZU(1,6),dp)
Write(123,*) "ZU16_i = ", AImag(ZU(1,6))
Else 
Write(123,*) "ZU11_i = ", Real(ZU(1,1),dp)
Write(123,*) "ZU11_r = ", -AImag(ZU(1,1))
Write(123,*) "ZU12_i = ", Real(ZU(1,2),dp)
Write(123,*) "ZU12_r = ", -AImag(ZU(1,2))
Write(123,*) "ZU13_i = ", Real(ZU(1,3),dp)
Write(123,*) "ZU13_r = ", -AImag(ZU(1,3))
Write(123,*) "ZU14_i = ", Real(ZU(1,4),dp)
Write(123,*) "ZU14_r = ", -AImag(ZU(1,4))
Write(123,*) "ZU15_i = ", Real(ZU(1,5),dp)
Write(123,*) "ZU15_r = ", -AImag(ZU(1,5))
Write(123,*) "ZU16_i = ", Real(ZU(1,6),dp)
Write(123,*) "ZU16_r = ", -AImag(ZU(1,6))
End if 
If (MSu(2).Gt.0._dp) Then 
Write(123,*) "ZU21_r = ", Real(ZU(2,1),dp)
Write(123,*) "ZU21_i = ", AImag(ZU(2,1))
Write(123,*) "ZU22_r = ", Real(ZU(2,2),dp)
Write(123,*) "ZU22_i = ", AImag(ZU(2,2))
Write(123,*) "ZU23_r = ", Real(ZU(2,3),dp)
Write(123,*) "ZU23_i = ", AImag(ZU(2,3))
Write(123,*) "ZU24_r = ", Real(ZU(2,4),dp)
Write(123,*) "ZU24_i = ", AImag(ZU(2,4))
Write(123,*) "ZU25_r = ", Real(ZU(2,5),dp)
Write(123,*) "ZU25_i = ", AImag(ZU(2,5))
Write(123,*) "ZU26_r = ", Real(ZU(2,6),dp)
Write(123,*) "ZU26_i = ", AImag(ZU(2,6))
Else 
Write(123,*) "ZU21_i = ", Real(ZU(2,1),dp)
Write(123,*) "ZU21_r = ", -AImag(ZU(2,1))
Write(123,*) "ZU22_i = ", Real(ZU(2,2),dp)
Write(123,*) "ZU22_r = ", -AImag(ZU(2,2))
Write(123,*) "ZU23_i = ", Real(ZU(2,3),dp)
Write(123,*) "ZU23_r = ", -AImag(ZU(2,3))
Write(123,*) "ZU24_i = ", Real(ZU(2,4),dp)
Write(123,*) "ZU24_r = ", -AImag(ZU(2,4))
Write(123,*) "ZU25_i = ", Real(ZU(2,5),dp)
Write(123,*) "ZU25_r = ", -AImag(ZU(2,5))
Write(123,*) "ZU26_i = ", Real(ZU(2,6),dp)
Write(123,*) "ZU26_r = ", -AImag(ZU(2,6))
End if 
If (MSu(3).Gt.0._dp) Then 
Write(123,*) "ZU31_r = ", Real(ZU(3,1),dp)
Write(123,*) "ZU31_i = ", AImag(ZU(3,1))
Write(123,*) "ZU32_r = ", Real(ZU(3,2),dp)
Write(123,*) "ZU32_i = ", AImag(ZU(3,2))
Write(123,*) "ZU33_r = ", Real(ZU(3,3),dp)
Write(123,*) "ZU33_i = ", AImag(ZU(3,3))
Write(123,*) "ZU34_r = ", Real(ZU(3,4),dp)
Write(123,*) "ZU34_i = ", AImag(ZU(3,4))
Write(123,*) "ZU35_r = ", Real(ZU(3,5),dp)
Write(123,*) "ZU35_i = ", AImag(ZU(3,5))
Write(123,*) "ZU36_r = ", Real(ZU(3,6),dp)
Write(123,*) "ZU36_i = ", AImag(ZU(3,6))
Else 
Write(123,*) "ZU31_i = ", Real(ZU(3,1),dp)
Write(123,*) "ZU31_r = ", -AImag(ZU(3,1))
Write(123,*) "ZU32_i = ", Real(ZU(3,2),dp)
Write(123,*) "ZU32_r = ", -AImag(ZU(3,2))
Write(123,*) "ZU33_i = ", Real(ZU(3,3),dp)
Write(123,*) "ZU33_r = ", -AImag(ZU(3,3))
Write(123,*) "ZU34_i = ", Real(ZU(3,4),dp)
Write(123,*) "ZU34_r = ", -AImag(ZU(3,4))
Write(123,*) "ZU35_i = ", Real(ZU(3,5),dp)
Write(123,*) "ZU35_r = ", -AImag(ZU(3,5))
Write(123,*) "ZU36_i = ", Real(ZU(3,6),dp)
Write(123,*) "ZU36_r = ", -AImag(ZU(3,6))
End if 
If (MSu(4).Gt.0._dp) Then 
Write(123,*) "ZU41_r = ", Real(ZU(4,1),dp)
Write(123,*) "ZU41_i = ", AImag(ZU(4,1))
Write(123,*) "ZU42_r = ", Real(ZU(4,2),dp)
Write(123,*) "ZU42_i = ", AImag(ZU(4,2))
Write(123,*) "ZU43_r = ", Real(ZU(4,3),dp)
Write(123,*) "ZU43_i = ", AImag(ZU(4,3))
Write(123,*) "ZU44_r = ", Real(ZU(4,4),dp)
Write(123,*) "ZU44_i = ", AImag(ZU(4,4))
Write(123,*) "ZU45_r = ", Real(ZU(4,5),dp)
Write(123,*) "ZU45_i = ", AImag(ZU(4,5))
Write(123,*) "ZU46_r = ", Real(ZU(4,6),dp)
Write(123,*) "ZU46_i = ", AImag(ZU(4,6))
Else 
Write(123,*) "ZU41_i = ", Real(ZU(4,1),dp)
Write(123,*) "ZU41_r = ", -AImag(ZU(4,1))
Write(123,*) "ZU42_i = ", Real(ZU(4,2),dp)
Write(123,*) "ZU42_r = ", -AImag(ZU(4,2))
Write(123,*) "ZU43_i = ", Real(ZU(4,3),dp)
Write(123,*) "ZU43_r = ", -AImag(ZU(4,3))
Write(123,*) "ZU44_i = ", Real(ZU(4,4),dp)
Write(123,*) "ZU44_r = ", -AImag(ZU(4,4))
Write(123,*) "ZU45_i = ", Real(ZU(4,5),dp)
Write(123,*) "ZU45_r = ", -AImag(ZU(4,5))
Write(123,*) "ZU46_i = ", Real(ZU(4,6),dp)
Write(123,*) "ZU46_r = ", -AImag(ZU(4,6))
End if 
If (MSu(5).Gt.0._dp) Then 
Write(123,*) "ZU51_r = ", Real(ZU(5,1),dp)
Write(123,*) "ZU51_i = ", AImag(ZU(5,1))
Write(123,*) "ZU52_r = ", Real(ZU(5,2),dp)
Write(123,*) "ZU52_i = ", AImag(ZU(5,2))
Write(123,*) "ZU53_r = ", Real(ZU(5,3),dp)
Write(123,*) "ZU53_i = ", AImag(ZU(5,3))
Write(123,*) "ZU54_r = ", Real(ZU(5,4),dp)
Write(123,*) "ZU54_i = ", AImag(ZU(5,4))
Write(123,*) "ZU55_r = ", Real(ZU(5,5),dp)
Write(123,*) "ZU55_i = ", AImag(ZU(5,5))
Write(123,*) "ZU56_r = ", Real(ZU(5,6),dp)
Write(123,*) "ZU56_i = ", AImag(ZU(5,6))
Else 
Write(123,*) "ZU51_i = ", Real(ZU(5,1),dp)
Write(123,*) "ZU51_r = ", -AImag(ZU(5,1))
Write(123,*) "ZU52_i = ", Real(ZU(5,2),dp)
Write(123,*) "ZU52_r = ", -AImag(ZU(5,2))
Write(123,*) "ZU53_i = ", Real(ZU(5,3),dp)
Write(123,*) "ZU53_r = ", -AImag(ZU(5,3))
Write(123,*) "ZU54_i = ", Real(ZU(5,4),dp)
Write(123,*) "ZU54_r = ", -AImag(ZU(5,4))
Write(123,*) "ZU55_i = ", Real(ZU(5,5),dp)
Write(123,*) "ZU55_r = ", -AImag(ZU(5,5))
Write(123,*) "ZU56_i = ", Real(ZU(5,6),dp)
Write(123,*) "ZU56_r = ", -AImag(ZU(5,6))
End if 
If (MSu(6).Gt.0._dp) Then 
Write(123,*) "ZU61_r = ", Real(ZU(6,1),dp)
Write(123,*) "ZU61_i = ", AImag(ZU(6,1))
Write(123,*) "ZU62_r = ", Real(ZU(6,2),dp)
Write(123,*) "ZU62_i = ", AImag(ZU(6,2))
Write(123,*) "ZU63_r = ", Real(ZU(6,3),dp)
Write(123,*) "ZU63_i = ", AImag(ZU(6,3))
Write(123,*) "ZU64_r = ", Real(ZU(6,4),dp)
Write(123,*) "ZU64_i = ", AImag(ZU(6,4))
Write(123,*) "ZU65_r = ", Real(ZU(6,5),dp)
Write(123,*) "ZU65_i = ", AImag(ZU(6,5))
Write(123,*) "ZU66_r = ", Real(ZU(6,6),dp)
Write(123,*) "ZU66_i = ", AImag(ZU(6,6))
Else 
Write(123,*) "ZU61_i = ", Real(ZU(6,1),dp)
Write(123,*) "ZU61_r = ", -AImag(ZU(6,1))
Write(123,*) "ZU62_i = ", Real(ZU(6,2),dp)
Write(123,*) "ZU62_r = ", -AImag(ZU(6,2))
Write(123,*) "ZU63_i = ", Real(ZU(6,3),dp)
Write(123,*) "ZU63_r = ", -AImag(ZU(6,3))
Write(123,*) "ZU64_i = ", Real(ZU(6,4),dp)
Write(123,*) "ZU64_r = ", -AImag(ZU(6,4))
Write(123,*) "ZU65_i = ", Real(ZU(6,5),dp)
Write(123,*) "ZU65_r = ", -AImag(ZU(6,5))
Write(123,*) "ZU66_i = ", Real(ZU(6,6),dp)
Write(123,*) "ZU66_r = ", -AImag(ZU(6,6))
End if 
Write(123,*) "ZH11= ",ZH(1,1)
Write(123,*) "ZH12= ",ZH(1,2)
Write(123,*) "ZH13= ",ZH(1,3)
Write(123,*) "ZH14= ",ZH(1,4)
Write(123,*) "ZH15= ",ZH(1,5)
Write(123,*) "ZH16= ",ZH(1,6)
Write(123,*) "ZH17= ",ZH(1,7)
Write(123,*) "ZH18= ",ZH(1,8)
Write(123,*) "ZH21= ",ZH(2,1)
Write(123,*) "ZH22= ",ZH(2,2)
Write(123,*) "ZH23= ",ZH(2,3)
Write(123,*) "ZH24= ",ZH(2,4)
Write(123,*) "ZH25= ",ZH(2,5)
Write(123,*) "ZH26= ",ZH(2,6)
Write(123,*) "ZH27= ",ZH(2,7)
Write(123,*) "ZH28= ",ZH(2,8)
Write(123,*) "ZH31= ",ZH(3,1)
Write(123,*) "ZH32= ",ZH(3,2)
Write(123,*) "ZH33= ",ZH(3,3)
Write(123,*) "ZH34= ",ZH(3,4)
Write(123,*) "ZH35= ",ZH(3,5)
Write(123,*) "ZH36= ",ZH(3,6)
Write(123,*) "ZH37= ",ZH(3,7)
Write(123,*) "ZH38= ",ZH(3,8)
Write(123,*) "ZH41= ",ZH(4,1)
Write(123,*) "ZH42= ",ZH(4,2)
Write(123,*) "ZH43= ",ZH(4,3)
Write(123,*) "ZH44= ",ZH(4,4)
Write(123,*) "ZH45= ",ZH(4,5)
Write(123,*) "ZH46= ",ZH(4,6)
Write(123,*) "ZH47= ",ZH(4,7)
Write(123,*) "ZH48= ",ZH(4,8)
Write(123,*) "ZH51= ",ZH(5,1)
Write(123,*) "ZH52= ",ZH(5,2)
Write(123,*) "ZH53= ",ZH(5,3)
Write(123,*) "ZH54= ",ZH(5,4)
Write(123,*) "ZH55= ",ZH(5,5)
Write(123,*) "ZH56= ",ZH(5,6)
Write(123,*) "ZH57= ",ZH(5,7)
Write(123,*) "ZH58= ",ZH(5,8)
Write(123,*) "ZH61= ",ZH(6,1)
Write(123,*) "ZH62= ",ZH(6,2)
Write(123,*) "ZH63= ",ZH(6,3)
Write(123,*) "ZH64= ",ZH(6,4)
Write(123,*) "ZH65= ",ZH(6,5)
Write(123,*) "ZH66= ",ZH(6,6)
Write(123,*) "ZH67= ",ZH(6,7)
Write(123,*) "ZH68= ",ZH(6,8)
Write(123,*) "ZH71= ",ZH(7,1)
Write(123,*) "ZH72= ",ZH(7,2)
Write(123,*) "ZH73= ",ZH(7,3)
Write(123,*) "ZH74= ",ZH(7,4)
Write(123,*) "ZH75= ",ZH(7,5)
Write(123,*) "ZH76= ",ZH(7,6)
Write(123,*) "ZH77= ",ZH(7,7)
Write(123,*) "ZH78= ",ZH(7,8)
Write(123,*) "ZH81= ",ZH(8,1)
Write(123,*) "ZH82= ",ZH(8,2)
Write(123,*) "ZH83= ",ZH(8,3)
Write(123,*) "ZH84= ",ZH(8,4)
Write(123,*) "ZH85= ",ZH(8,5)
Write(123,*) "ZH86= ",ZH(8,6)
Write(123,*) "ZH87= ",ZH(8,7)
Write(123,*) "ZH88= ",ZH(8,8)
Write(123,*) "ZA11= ",ZA(1,1)
Write(123,*) "ZA12= ",ZA(1,2)
Write(123,*) "ZA13= ",ZA(1,3)
Write(123,*) "ZA14= ",ZA(1,4)
Write(123,*) "ZA15= ",ZA(1,5)
Write(123,*) "ZA16= ",ZA(1,6)
Write(123,*) "ZA17= ",ZA(1,7)
Write(123,*) "ZA18= ",ZA(1,8)
Write(123,*) "ZA21= ",ZA(2,1)
Write(123,*) "ZA22= ",ZA(2,2)
Write(123,*) "ZA23= ",ZA(2,3)
Write(123,*) "ZA24= ",ZA(2,4)
Write(123,*) "ZA25= ",ZA(2,5)
Write(123,*) "ZA26= ",ZA(2,6)
Write(123,*) "ZA27= ",ZA(2,7)
Write(123,*) "ZA28= ",ZA(2,8)
Write(123,*) "ZA31= ",ZA(3,1)
Write(123,*) "ZA32= ",ZA(3,2)
Write(123,*) "ZA33= ",ZA(3,3)
Write(123,*) "ZA34= ",ZA(3,4)
Write(123,*) "ZA35= ",ZA(3,5)
Write(123,*) "ZA36= ",ZA(3,6)
Write(123,*) "ZA37= ",ZA(3,7)
Write(123,*) "ZA38= ",ZA(3,8)
Write(123,*) "ZA41= ",ZA(4,1)
Write(123,*) "ZA42= ",ZA(4,2)
Write(123,*) "ZA43= ",ZA(4,3)
Write(123,*) "ZA44= ",ZA(4,4)
Write(123,*) "ZA45= ",ZA(4,5)
Write(123,*) "ZA46= ",ZA(4,6)
Write(123,*) "ZA47= ",ZA(4,7)
Write(123,*) "ZA48= ",ZA(4,8)
Write(123,*) "ZA51= ",ZA(5,1)
Write(123,*) "ZA52= ",ZA(5,2)
Write(123,*) "ZA53= ",ZA(5,3)
Write(123,*) "ZA54= ",ZA(5,4)
Write(123,*) "ZA55= ",ZA(5,5)
Write(123,*) "ZA56= ",ZA(5,6)
Write(123,*) "ZA57= ",ZA(5,7)
Write(123,*) "ZA58= ",ZA(5,8)
Write(123,*) "ZA61= ",ZA(6,1)
Write(123,*) "ZA62= ",ZA(6,2)
Write(123,*) "ZA63= ",ZA(6,3)
Write(123,*) "ZA64= ",ZA(6,4)
Write(123,*) "ZA65= ",ZA(6,5)
Write(123,*) "ZA66= ",ZA(6,6)
Write(123,*) "ZA67= ",ZA(6,7)
Write(123,*) "ZA68= ",ZA(6,8)
Write(123,*) "ZA71= ",ZA(7,1)
Write(123,*) "ZA72= ",ZA(7,2)
Write(123,*) "ZA73= ",ZA(7,3)
Write(123,*) "ZA74= ",ZA(7,4)
Write(123,*) "ZA75= ",ZA(7,5)
Write(123,*) "ZA76= ",ZA(7,6)
Write(123,*) "ZA77= ",ZA(7,7)
Write(123,*) "ZA78= ",ZA(7,8)
Write(123,*) "ZA81= ",ZA(8,1)
Write(123,*) "ZA82= ",ZA(8,2)
Write(123,*) "ZA83= ",ZA(8,3)
Write(123,*) "ZA84= ",ZA(8,4)
Write(123,*) "ZA85= ",ZA(8,5)
Write(123,*) "ZA86= ",ZA(8,6)
Write(123,*) "ZA87= ",ZA(8,7)
Write(123,*) "ZA88= ",ZA(8,8)
Write(123,*) "ZP11= ",ZP(1,1)
Write(123,*) "ZP12= ",ZP(1,2)
Write(123,*) "ZP13= ",ZP(1,3)
Write(123,*) "ZP14= ",ZP(1,4)
Write(123,*) "ZP15= ",ZP(1,5)
Write(123,*) "ZP16= ",ZP(1,6)
Write(123,*) "ZP17= ",ZP(1,7)
Write(123,*) "ZP18= ",ZP(1,8)
Write(123,*) "ZP21= ",ZP(2,1)
Write(123,*) "ZP22= ",ZP(2,2)
Write(123,*) "ZP23= ",ZP(2,3)
Write(123,*) "ZP24= ",ZP(2,4)
Write(123,*) "ZP25= ",ZP(2,5)
Write(123,*) "ZP26= ",ZP(2,6)
Write(123,*) "ZP27= ",ZP(2,7)
Write(123,*) "ZP28= ",ZP(2,8)
Write(123,*) "ZP31= ",ZP(3,1)
Write(123,*) "ZP32= ",ZP(3,2)
Write(123,*) "ZP33= ",ZP(3,3)
Write(123,*) "ZP34= ",ZP(3,4)
Write(123,*) "ZP35= ",ZP(3,5)
Write(123,*) "ZP36= ",ZP(3,6)
Write(123,*) "ZP37= ",ZP(3,7)
Write(123,*) "ZP38= ",ZP(3,8)
Write(123,*) "ZP41= ",ZP(4,1)
Write(123,*) "ZP42= ",ZP(4,2)
Write(123,*) "ZP43= ",ZP(4,3)
Write(123,*) "ZP44= ",ZP(4,4)
Write(123,*) "ZP45= ",ZP(4,5)
Write(123,*) "ZP46= ",ZP(4,6)
Write(123,*) "ZP47= ",ZP(4,7)
Write(123,*) "ZP48= ",ZP(4,8)
Write(123,*) "ZP51= ",ZP(5,1)
Write(123,*) "ZP52= ",ZP(5,2)
Write(123,*) "ZP53= ",ZP(5,3)
Write(123,*) "ZP54= ",ZP(5,4)
Write(123,*) "ZP55= ",ZP(5,5)
Write(123,*) "ZP56= ",ZP(5,6)
Write(123,*) "ZP57= ",ZP(5,7)
Write(123,*) "ZP58= ",ZP(5,8)
Write(123,*) "ZP61= ",ZP(6,1)
Write(123,*) "ZP62= ",ZP(6,2)
Write(123,*) "ZP63= ",ZP(6,3)
Write(123,*) "ZP64= ",ZP(6,4)
Write(123,*) "ZP65= ",ZP(6,5)
Write(123,*) "ZP66= ",ZP(6,6)
Write(123,*) "ZP67= ",ZP(6,7)
Write(123,*) "ZP68= ",ZP(6,8)
Write(123,*) "ZP71= ",ZP(7,1)
Write(123,*) "ZP72= ",ZP(7,2)
Write(123,*) "ZP73= ",ZP(7,3)
Write(123,*) "ZP74= ",ZP(7,4)
Write(123,*) "ZP75= ",ZP(7,5)
Write(123,*) "ZP76= ",ZP(7,6)
Write(123,*) "ZP77= ",ZP(7,7)
Write(123,*) "ZP78= ",ZP(7,8)
Write(123,*) "ZP81= ",ZP(8,1)
Write(123,*) "ZP82= ",ZP(8,2)
Write(123,*) "ZP83= ",ZP(8,3)
Write(123,*) "ZP84= ",ZP(8,4)
Write(123,*) "ZP85= ",ZP(8,5)
Write(123,*) "ZP86= ",ZP(8,6)
Write(123,*) "ZP87= ",ZP(8,7)
Write(123,*) "ZP88= ",ZP(8,8)
If (MChi(1).Gt.0._dp) Then 
Write(123,*) "UV11_r = ", Real(UV(1,1),dp)
Write(123,*) "UV11_i = ", AImag(UV(1,1))
Write(123,*) "UV12_r = ", Real(UV(1,2),dp)
Write(123,*) "UV12_i = ", AImag(UV(1,2))
Write(123,*) "UV13_r = ", Real(UV(1,3),dp)
Write(123,*) "UV13_i = ", AImag(UV(1,3))
Write(123,*) "UV14_r = ", Real(UV(1,4),dp)
Write(123,*) "UV14_i = ", AImag(UV(1,4))
Write(123,*) "UV15_r = ", Real(UV(1,5),dp)
Write(123,*) "UV15_i = ", AImag(UV(1,5))
Write(123,*) "UV16_r = ", Real(UV(1,6),dp)
Write(123,*) "UV16_i = ", AImag(UV(1,6))
Write(123,*) "UV17_r = ", Real(UV(1,7),dp)
Write(123,*) "UV17_i = ", AImag(UV(1,7))
Write(123,*) "UV18_r = ", Real(UV(1,8),dp)
Write(123,*) "UV18_i = ", AImag(UV(1,8))
Write(123,*) "UV19_r = ", Real(UV(1,9),dp)
Write(123,*) "UV19_i = ", AImag(UV(1,9))
Write(123,*) "UV110_r = ", Real(UV(1,10),dp)
Write(123,*) "UV110_i = ", AImag(UV(1,10))
Else 
Write(123,*) "UV11_i = ", Real(UV(1,1),dp)
Write(123,*) "UV11_r = ", -AImag(UV(1,1))
Write(123,*) "UV12_i = ", Real(UV(1,2),dp)
Write(123,*) "UV12_r = ", -AImag(UV(1,2))
Write(123,*) "UV13_i = ", Real(UV(1,3),dp)
Write(123,*) "UV13_r = ", -AImag(UV(1,3))
Write(123,*) "UV14_i = ", Real(UV(1,4),dp)
Write(123,*) "UV14_r = ", -AImag(UV(1,4))
Write(123,*) "UV15_i = ", Real(UV(1,5),dp)
Write(123,*) "UV15_r = ", -AImag(UV(1,5))
Write(123,*) "UV16_i = ", Real(UV(1,6),dp)
Write(123,*) "UV16_r = ", -AImag(UV(1,6))
Write(123,*) "UV17_i = ", Real(UV(1,7),dp)
Write(123,*) "UV17_r = ", -AImag(UV(1,7))
Write(123,*) "UV18_i = ", Real(UV(1,8),dp)
Write(123,*) "UV18_r = ", -AImag(UV(1,8))
Write(123,*) "UV19_i = ", Real(UV(1,9),dp)
Write(123,*) "UV19_r = ", -AImag(UV(1,9))
Write(123,*) "UV110_i = ", Real(UV(1,10),dp)
Write(123,*) "UV110_r = ", -AImag(UV(1,10))
End if 
If (MChi(2).Gt.0._dp) Then 
Write(123,*) "UV21_r = ", Real(UV(2,1),dp)
Write(123,*) "UV21_i = ", AImag(UV(2,1))
Write(123,*) "UV22_r = ", Real(UV(2,2),dp)
Write(123,*) "UV22_i = ", AImag(UV(2,2))
Write(123,*) "UV23_r = ", Real(UV(2,3),dp)
Write(123,*) "UV23_i = ", AImag(UV(2,3))
Write(123,*) "UV24_r = ", Real(UV(2,4),dp)
Write(123,*) "UV24_i = ", AImag(UV(2,4))
Write(123,*) "UV25_r = ", Real(UV(2,5),dp)
Write(123,*) "UV25_i = ", AImag(UV(2,5))
Write(123,*) "UV26_r = ", Real(UV(2,6),dp)
Write(123,*) "UV26_i = ", AImag(UV(2,6))
Write(123,*) "UV27_r = ", Real(UV(2,7),dp)
Write(123,*) "UV27_i = ", AImag(UV(2,7))
Write(123,*) "UV28_r = ", Real(UV(2,8),dp)
Write(123,*) "UV28_i = ", AImag(UV(2,8))
Write(123,*) "UV29_r = ", Real(UV(2,9),dp)
Write(123,*) "UV29_i = ", AImag(UV(2,9))
Write(123,*) "UV210_r = ", Real(UV(2,10),dp)
Write(123,*) "UV210_i = ", AImag(UV(2,10))
Else 
Write(123,*) "UV21_i = ", Real(UV(2,1),dp)
Write(123,*) "UV21_r = ", -AImag(UV(2,1))
Write(123,*) "UV22_i = ", Real(UV(2,2),dp)
Write(123,*) "UV22_r = ", -AImag(UV(2,2))
Write(123,*) "UV23_i = ", Real(UV(2,3),dp)
Write(123,*) "UV23_r = ", -AImag(UV(2,3))
Write(123,*) "UV24_i = ", Real(UV(2,4),dp)
Write(123,*) "UV24_r = ", -AImag(UV(2,4))
Write(123,*) "UV25_i = ", Real(UV(2,5),dp)
Write(123,*) "UV25_r = ", -AImag(UV(2,5))
Write(123,*) "UV26_i = ", Real(UV(2,6),dp)
Write(123,*) "UV26_r = ", -AImag(UV(2,6))
Write(123,*) "UV27_i = ", Real(UV(2,7),dp)
Write(123,*) "UV27_r = ", -AImag(UV(2,7))
Write(123,*) "UV28_i = ", Real(UV(2,8),dp)
Write(123,*) "UV28_r = ", -AImag(UV(2,8))
Write(123,*) "UV29_i = ", Real(UV(2,9),dp)
Write(123,*) "UV29_r = ", -AImag(UV(2,9))
Write(123,*) "UV210_i = ", Real(UV(2,10),dp)
Write(123,*) "UV210_r = ", -AImag(UV(2,10))
End if 
If (MChi(3).Gt.0._dp) Then 
Write(123,*) "UV31_r = ", Real(UV(3,1),dp)
Write(123,*) "UV31_i = ", AImag(UV(3,1))
Write(123,*) "UV32_r = ", Real(UV(3,2),dp)
Write(123,*) "UV32_i = ", AImag(UV(3,2))
Write(123,*) "UV33_r = ", Real(UV(3,3),dp)
Write(123,*) "UV33_i = ", AImag(UV(3,3))
Write(123,*) "UV34_r = ", Real(UV(3,4),dp)
Write(123,*) "UV34_i = ", AImag(UV(3,4))
Write(123,*) "UV35_r = ", Real(UV(3,5),dp)
Write(123,*) "UV35_i = ", AImag(UV(3,5))
Write(123,*) "UV36_r = ", Real(UV(3,6),dp)
Write(123,*) "UV36_i = ", AImag(UV(3,6))
Write(123,*) "UV37_r = ", Real(UV(3,7),dp)
Write(123,*) "UV37_i = ", AImag(UV(3,7))
Write(123,*) "UV38_r = ", Real(UV(3,8),dp)
Write(123,*) "UV38_i = ", AImag(UV(3,8))
Write(123,*) "UV39_r = ", Real(UV(3,9),dp)
Write(123,*) "UV39_i = ", AImag(UV(3,9))
Write(123,*) "UV310_r = ", Real(UV(3,10),dp)
Write(123,*) "UV310_i = ", AImag(UV(3,10))
Else 
Write(123,*) "UV31_i = ", Real(UV(3,1),dp)
Write(123,*) "UV31_r = ", -AImag(UV(3,1))
Write(123,*) "UV32_i = ", Real(UV(3,2),dp)
Write(123,*) "UV32_r = ", -AImag(UV(3,2))
Write(123,*) "UV33_i = ", Real(UV(3,3),dp)
Write(123,*) "UV33_r = ", -AImag(UV(3,3))
Write(123,*) "UV34_i = ", Real(UV(3,4),dp)
Write(123,*) "UV34_r = ", -AImag(UV(3,4))
Write(123,*) "UV35_i = ", Real(UV(3,5),dp)
Write(123,*) "UV35_r = ", -AImag(UV(3,5))
Write(123,*) "UV36_i = ", Real(UV(3,6),dp)
Write(123,*) "UV36_r = ", -AImag(UV(3,6))
Write(123,*) "UV37_i = ", Real(UV(3,7),dp)
Write(123,*) "UV37_r = ", -AImag(UV(3,7))
Write(123,*) "UV38_i = ", Real(UV(3,8),dp)
Write(123,*) "UV38_r = ", -AImag(UV(3,8))
Write(123,*) "UV39_i = ", Real(UV(3,9),dp)
Write(123,*) "UV39_r = ", -AImag(UV(3,9))
Write(123,*) "UV310_i = ", Real(UV(3,10),dp)
Write(123,*) "UV310_r = ", -AImag(UV(3,10))
End if 
If (MChi(4).Gt.0._dp) Then 
Write(123,*) "UV41_r = ", Real(UV(4,1),dp)
Write(123,*) "UV41_i = ", AImag(UV(4,1))
Write(123,*) "UV42_r = ", Real(UV(4,2),dp)
Write(123,*) "UV42_i = ", AImag(UV(4,2))
Write(123,*) "UV43_r = ", Real(UV(4,3),dp)
Write(123,*) "UV43_i = ", AImag(UV(4,3))
Write(123,*) "UV44_r = ", Real(UV(4,4),dp)
Write(123,*) "UV44_i = ", AImag(UV(4,4))
Write(123,*) "UV45_r = ", Real(UV(4,5),dp)
Write(123,*) "UV45_i = ", AImag(UV(4,5))
Write(123,*) "UV46_r = ", Real(UV(4,6),dp)
Write(123,*) "UV46_i = ", AImag(UV(4,6))
Write(123,*) "UV47_r = ", Real(UV(4,7),dp)
Write(123,*) "UV47_i = ", AImag(UV(4,7))
Write(123,*) "UV48_r = ", Real(UV(4,8),dp)
Write(123,*) "UV48_i = ", AImag(UV(4,8))
Write(123,*) "UV49_r = ", Real(UV(4,9),dp)
Write(123,*) "UV49_i = ", AImag(UV(4,9))
Write(123,*) "UV410_r = ", Real(UV(4,10),dp)
Write(123,*) "UV410_i = ", AImag(UV(4,10))
Else 
Write(123,*) "UV41_i = ", Real(UV(4,1),dp)
Write(123,*) "UV41_r = ", -AImag(UV(4,1))
Write(123,*) "UV42_i = ", Real(UV(4,2),dp)
Write(123,*) "UV42_r = ", -AImag(UV(4,2))
Write(123,*) "UV43_i = ", Real(UV(4,3),dp)
Write(123,*) "UV43_r = ", -AImag(UV(4,3))
Write(123,*) "UV44_i = ", Real(UV(4,4),dp)
Write(123,*) "UV44_r = ", -AImag(UV(4,4))
Write(123,*) "UV45_i = ", Real(UV(4,5),dp)
Write(123,*) "UV45_r = ", -AImag(UV(4,5))
Write(123,*) "UV46_i = ", Real(UV(4,6),dp)
Write(123,*) "UV46_r = ", -AImag(UV(4,6))
Write(123,*) "UV47_i = ", Real(UV(4,7),dp)
Write(123,*) "UV47_r = ", -AImag(UV(4,7))
Write(123,*) "UV48_i = ", Real(UV(4,8),dp)
Write(123,*) "UV48_r = ", -AImag(UV(4,8))
Write(123,*) "UV49_i = ", Real(UV(4,9),dp)
Write(123,*) "UV49_r = ", -AImag(UV(4,9))
Write(123,*) "UV410_i = ", Real(UV(4,10),dp)
Write(123,*) "UV410_r = ", -AImag(UV(4,10))
End if 
If (MChi(5).Gt.0._dp) Then 
Write(123,*) "UV51_r = ", Real(UV(5,1),dp)
Write(123,*) "UV51_i = ", AImag(UV(5,1))
Write(123,*) "UV52_r = ", Real(UV(5,2),dp)
Write(123,*) "UV52_i = ", AImag(UV(5,2))
Write(123,*) "UV53_r = ", Real(UV(5,3),dp)
Write(123,*) "UV53_i = ", AImag(UV(5,3))
Write(123,*) "UV54_r = ", Real(UV(5,4),dp)
Write(123,*) "UV54_i = ", AImag(UV(5,4))
Write(123,*) "UV55_r = ", Real(UV(5,5),dp)
Write(123,*) "UV55_i = ", AImag(UV(5,5))
Write(123,*) "UV56_r = ", Real(UV(5,6),dp)
Write(123,*) "UV56_i = ", AImag(UV(5,6))
Write(123,*) "UV57_r = ", Real(UV(5,7),dp)
Write(123,*) "UV57_i = ", AImag(UV(5,7))
Write(123,*) "UV58_r = ", Real(UV(5,8),dp)
Write(123,*) "UV58_i = ", AImag(UV(5,8))
Write(123,*) "UV59_r = ", Real(UV(5,9),dp)
Write(123,*) "UV59_i = ", AImag(UV(5,9))
Write(123,*) "UV510_r = ", Real(UV(5,10),dp)
Write(123,*) "UV510_i = ", AImag(UV(5,10))
Else 
Write(123,*) "UV51_i = ", Real(UV(5,1),dp)
Write(123,*) "UV51_r = ", -AImag(UV(5,1))
Write(123,*) "UV52_i = ", Real(UV(5,2),dp)
Write(123,*) "UV52_r = ", -AImag(UV(5,2))
Write(123,*) "UV53_i = ", Real(UV(5,3),dp)
Write(123,*) "UV53_r = ", -AImag(UV(5,3))
Write(123,*) "UV54_i = ", Real(UV(5,4),dp)
Write(123,*) "UV54_r = ", -AImag(UV(5,4))
Write(123,*) "UV55_i = ", Real(UV(5,5),dp)
Write(123,*) "UV55_r = ", -AImag(UV(5,5))
Write(123,*) "UV56_i = ", Real(UV(5,6),dp)
Write(123,*) "UV56_r = ", -AImag(UV(5,6))
Write(123,*) "UV57_i = ", Real(UV(5,7),dp)
Write(123,*) "UV57_r = ", -AImag(UV(5,7))
Write(123,*) "UV58_i = ", Real(UV(5,8),dp)
Write(123,*) "UV58_r = ", -AImag(UV(5,8))
Write(123,*) "UV59_i = ", Real(UV(5,9),dp)
Write(123,*) "UV59_r = ", -AImag(UV(5,9))
Write(123,*) "UV510_i = ", Real(UV(5,10),dp)
Write(123,*) "UV510_r = ", -AImag(UV(5,10))
End if 
If (MChi(6).Gt.0._dp) Then 
Write(123,*) "UV61_r = ", Real(UV(6,1),dp)
Write(123,*) "UV61_i = ", AImag(UV(6,1))
Write(123,*) "UV62_r = ", Real(UV(6,2),dp)
Write(123,*) "UV62_i = ", AImag(UV(6,2))
Write(123,*) "UV63_r = ", Real(UV(6,3),dp)
Write(123,*) "UV63_i = ", AImag(UV(6,3))
Write(123,*) "UV64_r = ", Real(UV(6,4),dp)
Write(123,*) "UV64_i = ", AImag(UV(6,4))
Write(123,*) "UV65_r = ", Real(UV(6,5),dp)
Write(123,*) "UV65_i = ", AImag(UV(6,5))
Write(123,*) "UV66_r = ", Real(UV(6,6),dp)
Write(123,*) "UV66_i = ", AImag(UV(6,6))
Write(123,*) "UV67_r = ", Real(UV(6,7),dp)
Write(123,*) "UV67_i = ", AImag(UV(6,7))
Write(123,*) "UV68_r = ", Real(UV(6,8),dp)
Write(123,*) "UV68_i = ", AImag(UV(6,8))
Write(123,*) "UV69_r = ", Real(UV(6,9),dp)
Write(123,*) "UV69_i = ", AImag(UV(6,9))
Write(123,*) "UV610_r = ", Real(UV(6,10),dp)
Write(123,*) "UV610_i = ", AImag(UV(6,10))
Else 
Write(123,*) "UV61_i = ", Real(UV(6,1),dp)
Write(123,*) "UV61_r = ", -AImag(UV(6,1))
Write(123,*) "UV62_i = ", Real(UV(6,2),dp)
Write(123,*) "UV62_r = ", -AImag(UV(6,2))
Write(123,*) "UV63_i = ", Real(UV(6,3),dp)
Write(123,*) "UV63_r = ", -AImag(UV(6,3))
Write(123,*) "UV64_i = ", Real(UV(6,4),dp)
Write(123,*) "UV64_r = ", -AImag(UV(6,4))
Write(123,*) "UV65_i = ", Real(UV(6,5),dp)
Write(123,*) "UV65_r = ", -AImag(UV(6,5))
Write(123,*) "UV66_i = ", Real(UV(6,6),dp)
Write(123,*) "UV66_r = ", -AImag(UV(6,6))
Write(123,*) "UV67_i = ", Real(UV(6,7),dp)
Write(123,*) "UV67_r = ", -AImag(UV(6,7))
Write(123,*) "UV68_i = ", Real(UV(6,8),dp)
Write(123,*) "UV68_r = ", -AImag(UV(6,8))
Write(123,*) "UV69_i = ", Real(UV(6,9),dp)
Write(123,*) "UV69_r = ", -AImag(UV(6,9))
Write(123,*) "UV610_i = ", Real(UV(6,10),dp)
Write(123,*) "UV610_r = ", -AImag(UV(6,10))
End if 
If (MChi(7).Gt.0._dp) Then 
Write(123,*) "UV71_r = ", Real(UV(7,1),dp)
Write(123,*) "UV71_i = ", AImag(UV(7,1))
Write(123,*) "UV72_r = ", Real(UV(7,2),dp)
Write(123,*) "UV72_i = ", AImag(UV(7,2))
Write(123,*) "UV73_r = ", Real(UV(7,3),dp)
Write(123,*) "UV73_i = ", AImag(UV(7,3))
Write(123,*) "UV74_r = ", Real(UV(7,4),dp)
Write(123,*) "UV74_i = ", AImag(UV(7,4))
Write(123,*) "UV75_r = ", Real(UV(7,5),dp)
Write(123,*) "UV75_i = ", AImag(UV(7,5))
Write(123,*) "UV76_r = ", Real(UV(7,6),dp)
Write(123,*) "UV76_i = ", AImag(UV(7,6))
Write(123,*) "UV77_r = ", Real(UV(7,7),dp)
Write(123,*) "UV77_i = ", AImag(UV(7,7))
Write(123,*) "UV78_r = ", Real(UV(7,8),dp)
Write(123,*) "UV78_i = ", AImag(UV(7,8))
Write(123,*) "UV79_r = ", Real(UV(7,9),dp)
Write(123,*) "UV79_i = ", AImag(UV(7,9))
Write(123,*) "UV710_r = ", Real(UV(7,10),dp)
Write(123,*) "UV710_i = ", AImag(UV(7,10))
Else 
Write(123,*) "UV71_i = ", Real(UV(7,1),dp)
Write(123,*) "UV71_r = ", -AImag(UV(7,1))
Write(123,*) "UV72_i = ", Real(UV(7,2),dp)
Write(123,*) "UV72_r = ", -AImag(UV(7,2))
Write(123,*) "UV73_i = ", Real(UV(7,3),dp)
Write(123,*) "UV73_r = ", -AImag(UV(7,3))
Write(123,*) "UV74_i = ", Real(UV(7,4),dp)
Write(123,*) "UV74_r = ", -AImag(UV(7,4))
Write(123,*) "UV75_i = ", Real(UV(7,5),dp)
Write(123,*) "UV75_r = ", -AImag(UV(7,5))
Write(123,*) "UV76_i = ", Real(UV(7,6),dp)
Write(123,*) "UV76_r = ", -AImag(UV(7,6))
Write(123,*) "UV77_i = ", Real(UV(7,7),dp)
Write(123,*) "UV77_r = ", -AImag(UV(7,7))
Write(123,*) "UV78_i = ", Real(UV(7,8),dp)
Write(123,*) "UV78_r = ", -AImag(UV(7,8))
Write(123,*) "UV79_i = ", Real(UV(7,9),dp)
Write(123,*) "UV79_r = ", -AImag(UV(7,9))
Write(123,*) "UV710_i = ", Real(UV(7,10),dp)
Write(123,*) "UV710_r = ", -AImag(UV(7,10))
End if 
If (MChi(8).Gt.0._dp) Then 
Write(123,*) "UV81_r = ", Real(UV(8,1),dp)
Write(123,*) "UV81_i = ", AImag(UV(8,1))
Write(123,*) "UV82_r = ", Real(UV(8,2),dp)
Write(123,*) "UV82_i = ", AImag(UV(8,2))
Write(123,*) "UV83_r = ", Real(UV(8,3),dp)
Write(123,*) "UV83_i = ", AImag(UV(8,3))
Write(123,*) "UV84_r = ", Real(UV(8,4),dp)
Write(123,*) "UV84_i = ", AImag(UV(8,4))
Write(123,*) "UV85_r = ", Real(UV(8,5),dp)
Write(123,*) "UV85_i = ", AImag(UV(8,5))
Write(123,*) "UV86_r = ", Real(UV(8,6),dp)
Write(123,*) "UV86_i = ", AImag(UV(8,6))
Write(123,*) "UV87_r = ", Real(UV(8,7),dp)
Write(123,*) "UV87_i = ", AImag(UV(8,7))
Write(123,*) "UV88_r = ", Real(UV(8,8),dp)
Write(123,*) "UV88_i = ", AImag(UV(8,8))
Write(123,*) "UV89_r = ", Real(UV(8,9),dp)
Write(123,*) "UV89_i = ", AImag(UV(8,9))
Write(123,*) "UV810_r = ", Real(UV(8,10),dp)
Write(123,*) "UV810_i = ", AImag(UV(8,10))
Else 
Write(123,*) "UV81_i = ", Real(UV(8,1),dp)
Write(123,*) "UV81_r = ", -AImag(UV(8,1))
Write(123,*) "UV82_i = ", Real(UV(8,2),dp)
Write(123,*) "UV82_r = ", -AImag(UV(8,2))
Write(123,*) "UV83_i = ", Real(UV(8,3),dp)
Write(123,*) "UV83_r = ", -AImag(UV(8,3))
Write(123,*) "UV84_i = ", Real(UV(8,4),dp)
Write(123,*) "UV84_r = ", -AImag(UV(8,4))
Write(123,*) "UV85_i = ", Real(UV(8,5),dp)
Write(123,*) "UV85_r = ", -AImag(UV(8,5))
Write(123,*) "UV86_i = ", Real(UV(8,6),dp)
Write(123,*) "UV86_r = ", -AImag(UV(8,6))
Write(123,*) "UV87_i = ", Real(UV(8,7),dp)
Write(123,*) "UV87_r = ", -AImag(UV(8,7))
Write(123,*) "UV88_i = ", Real(UV(8,8),dp)
Write(123,*) "UV88_r = ", -AImag(UV(8,8))
Write(123,*) "UV89_i = ", Real(UV(8,9),dp)
Write(123,*) "UV89_r = ", -AImag(UV(8,9))
Write(123,*) "UV810_i = ", Real(UV(8,10),dp)
Write(123,*) "UV810_r = ", -AImag(UV(8,10))
End if 
If (MChi(9).Gt.0._dp) Then 
Write(123,*) "UV91_r = ", Real(UV(9,1),dp)
Write(123,*) "UV91_i = ", AImag(UV(9,1))
Write(123,*) "UV92_r = ", Real(UV(9,2),dp)
Write(123,*) "UV92_i = ", AImag(UV(9,2))
Write(123,*) "UV93_r = ", Real(UV(9,3),dp)
Write(123,*) "UV93_i = ", AImag(UV(9,3))
Write(123,*) "UV94_r = ", Real(UV(9,4),dp)
Write(123,*) "UV94_i = ", AImag(UV(9,4))
Write(123,*) "UV95_r = ", Real(UV(9,5),dp)
Write(123,*) "UV95_i = ", AImag(UV(9,5))
Write(123,*) "UV96_r = ", Real(UV(9,6),dp)
Write(123,*) "UV96_i = ", AImag(UV(9,6))
Write(123,*) "UV97_r = ", Real(UV(9,7),dp)
Write(123,*) "UV97_i = ", AImag(UV(9,7))
Write(123,*) "UV98_r = ", Real(UV(9,8),dp)
Write(123,*) "UV98_i = ", AImag(UV(9,8))
Write(123,*) "UV99_r = ", Real(UV(9,9),dp)
Write(123,*) "UV99_i = ", AImag(UV(9,9))
Write(123,*) "UV910_r = ", Real(UV(9,10),dp)
Write(123,*) "UV910_i = ", AImag(UV(9,10))
Else 
Write(123,*) "UV91_i = ", Real(UV(9,1),dp)
Write(123,*) "UV91_r = ", -AImag(UV(9,1))
Write(123,*) "UV92_i = ", Real(UV(9,2),dp)
Write(123,*) "UV92_r = ", -AImag(UV(9,2))
Write(123,*) "UV93_i = ", Real(UV(9,3),dp)
Write(123,*) "UV93_r = ", -AImag(UV(9,3))
Write(123,*) "UV94_i = ", Real(UV(9,4),dp)
Write(123,*) "UV94_r = ", -AImag(UV(9,4))
Write(123,*) "UV95_i = ", Real(UV(9,5),dp)
Write(123,*) "UV95_r = ", -AImag(UV(9,5))
Write(123,*) "UV96_i = ", Real(UV(9,6),dp)
Write(123,*) "UV96_r = ", -AImag(UV(9,6))
Write(123,*) "UV97_i = ", Real(UV(9,7),dp)
Write(123,*) "UV97_r = ", -AImag(UV(9,7))
Write(123,*) "UV98_i = ", Real(UV(9,8),dp)
Write(123,*) "UV98_r = ", -AImag(UV(9,8))
Write(123,*) "UV99_i = ", Real(UV(9,9),dp)
Write(123,*) "UV99_r = ", -AImag(UV(9,9))
Write(123,*) "UV910_i = ", Real(UV(9,10),dp)
Write(123,*) "UV910_r = ", -AImag(UV(9,10))
End if 
If (MChi(10).Gt.0._dp) Then 
Write(123,*) "UV101_r = ", Real(UV(10,1),dp)
Write(123,*) "UV101_i = ", AImag(UV(10,1))
Write(123,*) "UV102_r = ", Real(UV(10,2),dp)
Write(123,*) "UV102_i = ", AImag(UV(10,2))
Write(123,*) "UV103_r = ", Real(UV(10,3),dp)
Write(123,*) "UV103_i = ", AImag(UV(10,3))
Write(123,*) "UV104_r = ", Real(UV(10,4),dp)
Write(123,*) "UV104_i = ", AImag(UV(10,4))
Write(123,*) "UV105_r = ", Real(UV(10,5),dp)
Write(123,*) "UV105_i = ", AImag(UV(10,5))
Write(123,*) "UV106_r = ", Real(UV(10,6),dp)
Write(123,*) "UV106_i = ", AImag(UV(10,6))
Write(123,*) "UV107_r = ", Real(UV(10,7),dp)
Write(123,*) "UV107_i = ", AImag(UV(10,7))
Write(123,*) "UV108_r = ", Real(UV(10,8),dp)
Write(123,*) "UV108_i = ", AImag(UV(10,8))
Write(123,*) "UV109_r = ", Real(UV(10,9),dp)
Write(123,*) "UV109_i = ", AImag(UV(10,9))
Write(123,*) "UV1010_r = ", Real(UV(10,10),dp)
Write(123,*) "UV1010_i = ", AImag(UV(10,10))
Else 
Write(123,*) "UV101_i = ", Real(UV(10,1),dp)
Write(123,*) "UV101_r = ", -AImag(UV(10,1))
Write(123,*) "UV102_i = ", Real(UV(10,2),dp)
Write(123,*) "UV102_r = ", -AImag(UV(10,2))
Write(123,*) "UV103_i = ", Real(UV(10,3),dp)
Write(123,*) "UV103_r = ", -AImag(UV(10,3))
Write(123,*) "UV104_i = ", Real(UV(10,4),dp)
Write(123,*) "UV104_r = ", -AImag(UV(10,4))
Write(123,*) "UV105_i = ", Real(UV(10,5),dp)
Write(123,*) "UV105_r = ", -AImag(UV(10,5))
Write(123,*) "UV106_i = ", Real(UV(10,6),dp)
Write(123,*) "UV106_r = ", -AImag(UV(10,6))
Write(123,*) "UV107_i = ", Real(UV(10,7),dp)
Write(123,*) "UV107_r = ", -AImag(UV(10,7))
Write(123,*) "UV108_i = ", Real(UV(10,8),dp)
Write(123,*) "UV108_r = ", -AImag(UV(10,8))
Write(123,*) "UV109_i = ", Real(UV(10,9),dp)
Write(123,*) "UV109_r = ", -AImag(UV(10,9))
Write(123,*) "UV1010_i = ", Real(UV(10,10),dp)
Write(123,*) "UV1010_r = ", -AImag(UV(10,10))
End if 
If (MFd(1).Gt.0._dp) Then 
Write(123,*) "ZDL11_r = ", Real(ZDL(1,1),dp)
Write(123,*) "ZDL11_i = ", AImag(ZDL(1,1))
Write(123,*) "ZDL12_r = ", Real(ZDL(1,2),dp)
Write(123,*) "ZDL12_i = ", AImag(ZDL(1,2))
Write(123,*) "ZDL13_r = ", Real(ZDL(1,3),dp)
Write(123,*) "ZDL13_i = ", AImag(ZDL(1,3))
Else 
Write(123,*) "ZDL11_i = ", Real(ZDL(1,1),dp)
Write(123,*) "ZDL11_r = ", -AImag(ZDL(1,1))
Write(123,*) "ZDL12_i = ", Real(ZDL(1,2),dp)
Write(123,*) "ZDL12_r = ", -AImag(ZDL(1,2))
Write(123,*) "ZDL13_i = ", Real(ZDL(1,3),dp)
Write(123,*) "ZDL13_r = ", -AImag(ZDL(1,3))
End if 
If (MFd(2).Gt.0._dp) Then 
Write(123,*) "ZDL21_r = ", Real(ZDL(2,1),dp)
Write(123,*) "ZDL21_i = ", AImag(ZDL(2,1))
Write(123,*) "ZDL22_r = ", Real(ZDL(2,2),dp)
Write(123,*) "ZDL22_i = ", AImag(ZDL(2,2))
Write(123,*) "ZDL23_r = ", Real(ZDL(2,3),dp)
Write(123,*) "ZDL23_i = ", AImag(ZDL(2,3))
Else 
Write(123,*) "ZDL21_i = ", Real(ZDL(2,1),dp)
Write(123,*) "ZDL21_r = ", -AImag(ZDL(2,1))
Write(123,*) "ZDL22_i = ", Real(ZDL(2,2),dp)
Write(123,*) "ZDL22_r = ", -AImag(ZDL(2,2))
Write(123,*) "ZDL23_i = ", Real(ZDL(2,3),dp)
Write(123,*) "ZDL23_r = ", -AImag(ZDL(2,3))
End if 
If (MFd(3).Gt.0._dp) Then 
Write(123,*) "ZDL31_r = ", Real(ZDL(3,1),dp)
Write(123,*) "ZDL31_i = ", AImag(ZDL(3,1))
Write(123,*) "ZDL32_r = ", Real(ZDL(3,2),dp)
Write(123,*) "ZDL32_i = ", AImag(ZDL(3,2))
Write(123,*) "ZDL33_r = ", Real(ZDL(3,3),dp)
Write(123,*) "ZDL33_i = ", AImag(ZDL(3,3))
Else 
Write(123,*) "ZDL31_i = ", Real(ZDL(3,1),dp)
Write(123,*) "ZDL31_r = ", -AImag(ZDL(3,1))
Write(123,*) "ZDL32_i = ", Real(ZDL(3,2),dp)
Write(123,*) "ZDL32_r = ", -AImag(ZDL(3,2))
Write(123,*) "ZDL33_i = ", Real(ZDL(3,3),dp)
Write(123,*) "ZDL33_r = ", -AImag(ZDL(3,3))
End if 
If (MFd(1).Gt.0._dp) Then 
Write(123,*) "ZDR11_r = ", Real(ZDR(1,1),dp)
Write(123,*) "ZDR11_i = ", AImag(ZDR(1,1))
Write(123,*) "ZDR12_r = ", Real(ZDR(1,2),dp)
Write(123,*) "ZDR12_i = ", AImag(ZDR(1,2))
Write(123,*) "ZDR13_r = ", Real(ZDR(1,3),dp)
Write(123,*) "ZDR13_i = ", AImag(ZDR(1,3))
Else 
Write(123,*) "ZDR11_i = ", Real(ZDR(1,1),dp)
Write(123,*) "ZDR11_r = ", -AImag(ZDR(1,1))
Write(123,*) "ZDR12_i = ", Real(ZDR(1,2),dp)
Write(123,*) "ZDR12_r = ", -AImag(ZDR(1,2))
Write(123,*) "ZDR13_i = ", Real(ZDR(1,3),dp)
Write(123,*) "ZDR13_r = ", -AImag(ZDR(1,3))
End if 
If (MFd(2).Gt.0._dp) Then 
Write(123,*) "ZDR21_r = ", Real(ZDR(2,1),dp)
Write(123,*) "ZDR21_i = ", AImag(ZDR(2,1))
Write(123,*) "ZDR22_r = ", Real(ZDR(2,2),dp)
Write(123,*) "ZDR22_i = ", AImag(ZDR(2,2))
Write(123,*) "ZDR23_r = ", Real(ZDR(2,3),dp)
Write(123,*) "ZDR23_i = ", AImag(ZDR(2,3))
Else 
Write(123,*) "ZDR21_i = ", Real(ZDR(2,1),dp)
Write(123,*) "ZDR21_r = ", -AImag(ZDR(2,1))
Write(123,*) "ZDR22_i = ", Real(ZDR(2,2),dp)
Write(123,*) "ZDR22_r = ", -AImag(ZDR(2,2))
Write(123,*) "ZDR23_i = ", Real(ZDR(2,3),dp)
Write(123,*) "ZDR23_r = ", -AImag(ZDR(2,3))
End if 
If (MFd(3).Gt.0._dp) Then 
Write(123,*) "ZDR31_r = ", Real(ZDR(3,1),dp)
Write(123,*) "ZDR31_i = ", AImag(ZDR(3,1))
Write(123,*) "ZDR32_r = ", Real(ZDR(3,2),dp)
Write(123,*) "ZDR32_i = ", AImag(ZDR(3,2))
Write(123,*) "ZDR33_r = ", Real(ZDR(3,3),dp)
Write(123,*) "ZDR33_i = ", AImag(ZDR(3,3))
Else 
Write(123,*) "ZDR31_i = ", Real(ZDR(3,1),dp)
Write(123,*) "ZDR31_r = ", -AImag(ZDR(3,1))
Write(123,*) "ZDR32_i = ", Real(ZDR(3,2),dp)
Write(123,*) "ZDR32_r = ", -AImag(ZDR(3,2))
Write(123,*) "ZDR33_i = ", Real(ZDR(3,3),dp)
Write(123,*) "ZDR33_r = ", -AImag(ZDR(3,3))
End if 
If (MFu(1).Gt.0._dp) Then 
Write(123,*) "ZUL11_r = ", Real(ZUL(1,1),dp)
Write(123,*) "ZUL11_i = ", AImag(ZUL(1,1))
Write(123,*) "ZUL12_r = ", Real(ZUL(1,2),dp)
Write(123,*) "ZUL12_i = ", AImag(ZUL(1,2))
Write(123,*) "ZUL13_r = ", Real(ZUL(1,3),dp)
Write(123,*) "ZUL13_i = ", AImag(ZUL(1,3))
Else 
Write(123,*) "ZUL11_i = ", Real(ZUL(1,1),dp)
Write(123,*) "ZUL11_r = ", -AImag(ZUL(1,1))
Write(123,*) "ZUL12_i = ", Real(ZUL(1,2),dp)
Write(123,*) "ZUL12_r = ", -AImag(ZUL(1,2))
Write(123,*) "ZUL13_i = ", Real(ZUL(1,3),dp)
Write(123,*) "ZUL13_r = ", -AImag(ZUL(1,3))
End if 
If (MFu(2).Gt.0._dp) Then 
Write(123,*) "ZUL21_r = ", Real(ZUL(2,1),dp)
Write(123,*) "ZUL21_i = ", AImag(ZUL(2,1))
Write(123,*) "ZUL22_r = ", Real(ZUL(2,2),dp)
Write(123,*) "ZUL22_i = ", AImag(ZUL(2,2))
Write(123,*) "ZUL23_r = ", Real(ZUL(2,3),dp)
Write(123,*) "ZUL23_i = ", AImag(ZUL(2,3))
Else 
Write(123,*) "ZUL21_i = ", Real(ZUL(2,1),dp)
Write(123,*) "ZUL21_r = ", -AImag(ZUL(2,1))
Write(123,*) "ZUL22_i = ", Real(ZUL(2,2),dp)
Write(123,*) "ZUL22_r = ", -AImag(ZUL(2,2))
Write(123,*) "ZUL23_i = ", Real(ZUL(2,3),dp)
Write(123,*) "ZUL23_r = ", -AImag(ZUL(2,3))
End if 
If (MFu(3).Gt.0._dp) Then 
Write(123,*) "ZUL31_r = ", Real(ZUL(3,1),dp)
Write(123,*) "ZUL31_i = ", AImag(ZUL(3,1))
Write(123,*) "ZUL32_r = ", Real(ZUL(3,2),dp)
Write(123,*) "ZUL32_i = ", AImag(ZUL(3,2))
Write(123,*) "ZUL33_r = ", Real(ZUL(3,3),dp)
Write(123,*) "ZUL33_i = ", AImag(ZUL(3,3))
Else 
Write(123,*) "ZUL31_i = ", Real(ZUL(3,1),dp)
Write(123,*) "ZUL31_r = ", -AImag(ZUL(3,1))
Write(123,*) "ZUL32_i = ", Real(ZUL(3,2),dp)
Write(123,*) "ZUL32_r = ", -AImag(ZUL(3,2))
Write(123,*) "ZUL33_i = ", Real(ZUL(3,3),dp)
Write(123,*) "ZUL33_r = ", -AImag(ZUL(3,3))
End if 
If (MFu(1).Gt.0._dp) Then 
Write(123,*) "ZUR11_r = ", Real(ZUR(1,1),dp)
Write(123,*) "ZUR11_i = ", AImag(ZUR(1,1))
Write(123,*) "ZUR12_r = ", Real(ZUR(1,2),dp)
Write(123,*) "ZUR12_i = ", AImag(ZUR(1,2))
Write(123,*) "ZUR13_r = ", Real(ZUR(1,3),dp)
Write(123,*) "ZUR13_i = ", AImag(ZUR(1,3))
Else 
Write(123,*) "ZUR11_i = ", Real(ZUR(1,1),dp)
Write(123,*) "ZUR11_r = ", -AImag(ZUR(1,1))
Write(123,*) "ZUR12_i = ", Real(ZUR(1,2),dp)
Write(123,*) "ZUR12_r = ", -AImag(ZUR(1,2))
Write(123,*) "ZUR13_i = ", Real(ZUR(1,3),dp)
Write(123,*) "ZUR13_r = ", -AImag(ZUR(1,3))
End if 
If (MFu(2).Gt.0._dp) Then 
Write(123,*) "ZUR21_r = ", Real(ZUR(2,1),dp)
Write(123,*) "ZUR21_i = ", AImag(ZUR(2,1))
Write(123,*) "ZUR22_r = ", Real(ZUR(2,2),dp)
Write(123,*) "ZUR22_i = ", AImag(ZUR(2,2))
Write(123,*) "ZUR23_r = ", Real(ZUR(2,3),dp)
Write(123,*) "ZUR23_i = ", AImag(ZUR(2,3))
Else 
Write(123,*) "ZUR21_i = ", Real(ZUR(2,1),dp)
Write(123,*) "ZUR21_r = ", -AImag(ZUR(2,1))
Write(123,*) "ZUR22_i = ", Real(ZUR(2,2),dp)
Write(123,*) "ZUR22_r = ", -AImag(ZUR(2,2))
Write(123,*) "ZUR23_i = ", Real(ZUR(2,3),dp)
Write(123,*) "ZUR23_r = ", -AImag(ZUR(2,3))
End if 
If (MFu(3).Gt.0._dp) Then 
Write(123,*) "ZUR31_r = ", Real(ZUR(3,1),dp)
Write(123,*) "ZUR31_i = ", AImag(ZUR(3,1))
Write(123,*) "ZUR32_r = ", Real(ZUR(3,2),dp)
Write(123,*) "ZUR32_i = ", AImag(ZUR(3,2))
Write(123,*) "ZUR33_r = ", Real(ZUR(3,3),dp)
Write(123,*) "ZUR33_i = ", AImag(ZUR(3,3))
Else 
Write(123,*) "ZUR31_i = ", Real(ZUR(3,1),dp)
Write(123,*) "ZUR31_r = ", -AImag(ZUR(3,1))
Write(123,*) "ZUR32_i = ", Real(ZUR(3,2),dp)
Write(123,*) "ZUR32_r = ", -AImag(ZUR(3,2))
Write(123,*) "ZUR33_i = ", Real(ZUR(3,3),dp)
Write(123,*) "ZUR33_r = ", -AImag(ZUR(3,3))
End if 
Write(123,*) "ZER11_r = ", Real(ZER(1,1),dp)
Write(123,*) "ZER11_i = ", AImag(ZER(1,1))
Write(123,*) "ZER12_r = ", Real(ZER(1,2),dp)
Write(123,*) "ZER12_i = ", AImag(ZER(1,2))
Write(123,*) "ZER13_r = ", Real(ZER(1,3),dp)
Write(123,*) "ZER13_i = ", AImag(ZER(1,3))
Write(123,*) "ZER14_r = ", Real(ZER(1,4),dp)
Write(123,*) "ZER14_i = ", AImag(ZER(1,4))
Write(123,*) "ZER15_r = ", Real(ZER(1,5),dp)
Write(123,*) "ZER15_i = ", AImag(ZER(1,5))
Write(123,*) "ZER21_r = ", Real(ZER(2,1),dp)
Write(123,*) "ZER21_i = ", AImag(ZER(2,1))
Write(123,*) "ZER22_r = ", Real(ZER(2,2),dp)
Write(123,*) "ZER22_i = ", AImag(ZER(2,2))
Write(123,*) "ZER23_r = ", Real(ZER(2,3),dp)
Write(123,*) "ZER23_i = ", AImag(ZER(2,3))
Write(123,*) "ZER24_r = ", Real(ZER(2,4),dp)
Write(123,*) "ZER24_i = ", AImag(ZER(2,4))
Write(123,*) "ZER25_r = ", Real(ZER(2,5),dp)
Write(123,*) "ZER25_i = ", AImag(ZER(2,5))
Write(123,*) "ZER31_r = ", Real(ZER(3,1),dp)
Write(123,*) "ZER31_i = ", AImag(ZER(3,1))
Write(123,*) "ZER32_r = ", Real(ZER(3,2),dp)
Write(123,*) "ZER32_i = ", AImag(ZER(3,2))
Write(123,*) "ZER33_r = ", Real(ZER(3,3),dp)
Write(123,*) "ZER33_i = ", AImag(ZER(3,3))
Write(123,*) "ZER34_r = ", Real(ZER(3,4),dp)
Write(123,*) "ZER34_i = ", AImag(ZER(3,4))
Write(123,*) "ZER35_r = ", Real(ZER(3,5),dp)
Write(123,*) "ZER35_i = ", AImag(ZER(3,5))
Write(123,*) "ZER41_r = ", Real(ZER(4,1),dp)
Write(123,*) "ZER41_i = ", AImag(ZER(4,1))
Write(123,*) "ZER42_r = ", Real(ZER(4,2),dp)
Write(123,*) "ZER42_i = ", AImag(ZER(4,2))
Write(123,*) "ZER43_r = ", Real(ZER(4,3),dp)
Write(123,*) "ZER43_i = ", AImag(ZER(4,3))
Write(123,*) "ZER44_r = ", Real(ZER(4,4),dp)
Write(123,*) "ZER44_i = ", AImag(ZER(4,4))
Write(123,*) "ZER45_r = ", Real(ZER(4,5),dp)
Write(123,*) "ZER45_i = ", AImag(ZER(4,5))
Write(123,*) "ZER51_r = ", Real(ZER(5,1),dp)
Write(123,*) "ZER51_i = ", AImag(ZER(5,1))
Write(123,*) "ZER52_r = ", Real(ZER(5,2),dp)
Write(123,*) "ZER52_i = ", AImag(ZER(5,2))
Write(123,*) "ZER53_r = ", Real(ZER(5,3),dp)
Write(123,*) "ZER53_i = ", AImag(ZER(5,3))
Write(123,*) "ZER54_r = ", Real(ZER(5,4),dp)
Write(123,*) "ZER54_i = ", AImag(ZER(5,4))
Write(123,*) "ZER55_r = ", Real(ZER(5,5),dp)
Write(123,*) "ZER55_i = ", AImag(ZER(5,5))
Write(123,*) "ZEL11_r = ", Real(ZEL(1,1),dp)
Write(123,*) "ZEL11_i = ", AImag(ZEL(1,1))
Write(123,*) "ZEL12_r = ", Real(ZEL(1,2),dp)
Write(123,*) "ZEL12_i = ", AImag(ZEL(1,2))
Write(123,*) "ZEL13_r = ", Real(ZEL(1,3),dp)
Write(123,*) "ZEL13_i = ", AImag(ZEL(1,3))
Write(123,*) "ZEL14_r = ", Real(ZEL(1,4),dp)
Write(123,*) "ZEL14_i = ", AImag(ZEL(1,4))
Write(123,*) "ZEL15_r = ", Real(ZEL(1,5),dp)
Write(123,*) "ZEL15_i = ", AImag(ZEL(1,5))
Write(123,*) "ZEL21_r = ", Real(ZEL(2,1),dp)
Write(123,*) "ZEL21_i = ", AImag(ZEL(2,1))
Write(123,*) "ZEL22_r = ", Real(ZEL(2,2),dp)
Write(123,*) "ZEL22_i = ", AImag(ZEL(2,2))
Write(123,*) "ZEL23_r = ", Real(ZEL(2,3),dp)
Write(123,*) "ZEL23_i = ", AImag(ZEL(2,3))
Write(123,*) "ZEL24_r = ", Real(ZEL(2,4),dp)
Write(123,*) "ZEL24_i = ", AImag(ZEL(2,4))
Write(123,*) "ZEL25_r = ", Real(ZEL(2,5),dp)
Write(123,*) "ZEL25_i = ", AImag(ZEL(2,5))
Write(123,*) "ZEL31_r = ", Real(ZEL(3,1),dp)
Write(123,*) "ZEL31_i = ", AImag(ZEL(3,1))
Write(123,*) "ZEL32_r = ", Real(ZEL(3,2),dp)
Write(123,*) "ZEL32_i = ", AImag(ZEL(3,2))
Write(123,*) "ZEL33_r = ", Real(ZEL(3,3),dp)
Write(123,*) "ZEL33_i = ", AImag(ZEL(3,3))
Write(123,*) "ZEL34_r = ", Real(ZEL(3,4),dp)
Write(123,*) "ZEL34_i = ", AImag(ZEL(3,4))
Write(123,*) "ZEL35_r = ", Real(ZEL(3,5),dp)
Write(123,*) "ZEL35_i = ", AImag(ZEL(3,5))
Write(123,*) "ZEL41_r = ", Real(ZEL(4,1),dp)
Write(123,*) "ZEL41_i = ", AImag(ZEL(4,1))
Write(123,*) "ZEL42_r = ", Real(ZEL(4,2),dp)
Write(123,*) "ZEL42_i = ", AImag(ZEL(4,2))
Write(123,*) "ZEL43_r = ", Real(ZEL(4,3),dp)
Write(123,*) "ZEL43_i = ", AImag(ZEL(4,3))
Write(123,*) "ZEL44_r = ", Real(ZEL(4,4),dp)
Write(123,*) "ZEL44_i = ", AImag(ZEL(4,4))
Write(123,*) "ZEL45_r = ", Real(ZEL(4,5),dp)
Write(123,*) "ZEL45_i = ", AImag(ZEL(4,5))
Write(123,*) "ZEL51_r = ", Real(ZEL(5,1),dp)
Write(123,*) "ZEL51_i = ", AImag(ZEL(5,1))
Write(123,*) "ZEL52_r = ", Real(ZEL(5,2),dp)
Write(123,*) "ZEL52_i = ", AImag(ZEL(5,2))
Write(123,*) "ZEL53_r = ", Real(ZEL(5,3),dp)
Write(123,*) "ZEL53_i = ", AImag(ZEL(5,3))
Write(123,*) "ZEL54_r = ", Real(ZEL(5,4),dp)
Write(123,*) "ZEL54_i = ", AImag(ZEL(5,4))
Write(123,*) "ZEL55_r = ", Real(ZEL(5,5),dp)
Write(123,*) "ZEL55_i = ", AImag(ZEL(5,5))
    Close(123) 
End Subroutine WriteWHIZARD 

 
 Subroutine WriteHiggsBounds 
Implicit None 
Open(87,file="MH_GammaTot.dat",status="unknown") 
Open(88,file="MHplus_GammaTot.dat",status="unknown") 
Open(89,file="effC.dat",status="unknown") 
Open(90,file="BR_H_NP.dat",status="unknown") 
Open(91,file="BR_Hplus.dat",status="unknown") 
Open(92,file="BR_t.dat",status="unknown") 
Open(93,file="LEP_HpHm_CS_ratios.dat",status="unknown") 
Write(87,"(I1)",advance="No") 1 
Write(88,"(I1)",advance="No") 1 
Write(87,"(2e16.8)",advance="No") Mhh(1)
Write(87,"(2e16.8)",advance="No") Mhh(2)
Write(87,"(2e16.8)",advance="No") Mhh(3)
Write(87,"(2e16.8)",advance="No") Mhh(4)
Write(87,"(2e16.8)",advance="No") Mhh(5)
Write(87,"(2e16.8)",advance="No") Mhh(6)
Write(87,"(2e16.8)",advance="No") Mhh(7)
Write(87,"(2e16.8)",advance="No") Mhh(8)
Write(87,"(2e16.8)",advance="No") MAh(2)
Write(87,"(2e16.8)",advance="No") MAh(3)
Write(87,"(2e16.8)",advance="No") MAh(4)
Write(87,"(2e16.8)",advance="No") MAh(5)
Write(87,"(2e16.8)",advance="No") MAh(6)
Write(87,"(2e16.8)",advance="No") MAh(7)
Write(87,"(2e16.8)",advance="No") MAh(8)
Write(88,"(2e16.8)",advance="No") MHpm(2)
Write(88,"(2e16.8)",advance="No") MHpm(3)
Write(88,"(2e16.8)",advance="No") MHpm(4)
Write(88,"(2e16.8)",advance="No") MHpm(5)
Write(88,"(2e16.8)",advance="No") MHpm(6)
Write(88,"(2e16.8)",advance="No") MHpm(7)
Write(88,"(2e16.8)",advance="No") MHpm(8)
Write(87,"(2e16.8)",advance="No") gThh(1)
Write(87,"(2e16.8)",advance="No") gThh(2)
Write(87,"(2e16.8)",advance="No") gThh(3)
Write(87,"(2e16.8)",advance="No") gThh(4)
Write(87,"(2e16.8)",advance="No") gThh(5)
Write(87,"(2e16.8)",advance="No") gThh(6)
Write(87,"(2e16.8)",advance="No") gThh(7)
Write(87,"(2e16.8)",advance="No") gThh(8)
Write(87,"(2e16.8)",advance="No") gTAh(2)
Write(87,"(2e16.8)",advance="No") gTAh(3)
Write(87,"(2e16.8)",advance="No") gTAh(4)
Write(87,"(2e16.8)",advance="No") gTAh(5)
Write(87,"(2e16.8)",advance="No") gTAh(6)
Write(87,"(2e16.8)",advance="No") gTAh(7)
Write(87,"(2e16.8)",advance="No") gTAh(8)
Write(88,"(2e16.8)",advance="No") gTHpm(2)
Write(88,"(2e16.8)",advance="No") gTHpm(3)
Write(88,"(2e16.8)",advance="No") gTHpm(4)
Write(88,"(2e16.8)",advance="No") gTHpm(5)
Write(88,"(2e16.8)",advance="No") gTHpm(6)
Write(88,"(2e16.8)",advance="No") gTHpm(7)
Write(88,"(2e16.8)",advance="No") gTHpm(8)
Write(90,"(I1)",advance="No") 1 
Write(90,"(e16.8)",advance="No") BRinvH(1) 
Write(90,"(e16.8)",advance="No") BRinvH(2) 
Write(90,"(e16.8)",advance="No") BRinvH(3) 
Write(90,"(e16.8)",advance="No") BRinvH(4) 
Write(90,"(e16.8)",advance="No") BRinvH(5) 
Write(90,"(e16.8)",advance="No") BRinvH(6) 
Write(90,"(e16.8)",advance="No") BRinvH(7) 
Write(90,"(e16.8)",advance="No") BRinvH(8) 
Write(90,"(e16.8)",advance="No") BRinvA(2) 
Write(90,"(e16.8)",advance="No") BRinvA(3) 
Write(90,"(e16.8)",advance="No") BRinvA(4) 
Write(90,"(e16.8)",advance="No") BRinvA(5) 
Write(90,"(e16.8)",advance="No") BRinvA(6) 
Write(90,"(e16.8)",advance="No") BRinvA(7) 
Write(90,"(e16.8)",advance="No") BRinvA(8) 
Write(90,"(e16.8)",advance="No") BRHHH(1,2) 
Write(90,"(e16.8)",advance="No") BRHHH(1,3) 
Write(90,"(e16.8)",advance="No") BRHHH(1,4) 
Write(90,"(e16.8)",advance="No") BRHHH(1,5) 
Write(90,"(e16.8)",advance="No") BRHHH(1,6) 
Write(90,"(e16.8)",advance="No") BRHHH(1,7) 
Write(90,"(e16.8)",advance="No") BRHHH(1,8) 
Write(90,"(e16.8)",advance="No") BRHAA(1,2) 
Write(90,"(e16.8)",advance="No") BRHAA(1,3) 
Write(90,"(e16.8)",advance="No") BRHAA(1,4) 
Write(90,"(e16.8)",advance="No") BRHAA(1,5) 
Write(90,"(e16.8)",advance="No") BRHAA(1,6) 
Write(90,"(e16.8)",advance="No") BRHAA(1,7) 
Write(90,"(e16.8)",advance="No") BRHAA(1,8) 
Write(90,"(e16.8)",advance="No") BRHHH(2,1) 
Write(90,"(e16.8)",advance="No") BRHHH(2,3) 
Write(90,"(e16.8)",advance="No") BRHHH(2,4) 
Write(90,"(e16.8)",advance="No") BRHHH(2,5) 
Write(90,"(e16.8)",advance="No") BRHHH(2,6) 
Write(90,"(e16.8)",advance="No") BRHHH(2,7) 
Write(90,"(e16.8)",advance="No") BRHHH(2,8) 
Write(90,"(e16.8)",advance="No") BRHAA(2,2) 
Write(90,"(e16.8)",advance="No") BRHAA(2,3) 
Write(90,"(e16.8)",advance="No") BRHAA(2,4) 
Write(90,"(e16.8)",advance="No") BRHAA(2,5) 
Write(90,"(e16.8)",advance="No") BRHAA(2,6) 
Write(90,"(e16.8)",advance="No") BRHAA(2,7) 
Write(90,"(e16.8)",advance="No") BRHAA(2,8) 
Write(90,"(e16.8)",advance="No") BRHHH(3,1) 
Write(90,"(e16.8)",advance="No") BRHHH(3,2) 
Write(90,"(e16.8)",advance="No") BRHHH(3,4) 
Write(90,"(e16.8)",advance="No") BRHHH(3,5) 
Write(90,"(e16.8)",advance="No") BRHHH(3,6) 
Write(90,"(e16.8)",advance="No") BRHHH(3,7) 
Write(90,"(e16.8)",advance="No") BRHHH(3,8) 
Write(90,"(e16.8)",advance="No") BRHAA(3,2) 
Write(90,"(e16.8)",advance="No") BRHAA(3,3) 
Write(90,"(e16.8)",advance="No") BRHAA(3,4) 
Write(90,"(e16.8)",advance="No") BRHAA(3,5) 
Write(90,"(e16.8)",advance="No") BRHAA(3,6) 
Write(90,"(e16.8)",advance="No") BRHAA(3,7) 
Write(90,"(e16.8)",advance="No") BRHAA(3,8) 
Write(90,"(e16.8)",advance="No") BRHHH(4,1) 
Write(90,"(e16.8)",advance="No") BRHHH(4,2) 
Write(90,"(e16.8)",advance="No") BRHHH(4,3) 
Write(90,"(e16.8)",advance="No") BRHHH(4,5) 
Write(90,"(e16.8)",advance="No") BRHHH(4,6) 
Write(90,"(e16.8)",advance="No") BRHHH(4,7) 
Write(90,"(e16.8)",advance="No") BRHHH(4,8) 
Write(90,"(e16.8)",advance="No") BRHAA(4,2) 
Write(90,"(e16.8)",advance="No") BRHAA(4,3) 
Write(90,"(e16.8)",advance="No") BRHAA(4,4) 
Write(90,"(e16.8)",advance="No") BRHAA(4,5) 
Write(90,"(e16.8)",advance="No") BRHAA(4,6) 
Write(90,"(e16.8)",advance="No") BRHAA(4,7) 
Write(90,"(e16.8)",advance="No") BRHAA(4,8) 
Write(90,"(e16.8)",advance="No") BRHHH(5,1) 
Write(90,"(e16.8)",advance="No") BRHHH(5,2) 
Write(90,"(e16.8)",advance="No") BRHHH(5,3) 
Write(90,"(e16.8)",advance="No") BRHHH(5,4) 
Write(90,"(e16.8)",advance="No") BRHHH(5,6) 
Write(90,"(e16.8)",advance="No") BRHHH(5,7) 
Write(90,"(e16.8)",advance="No") BRHHH(5,8) 
Write(90,"(e16.8)",advance="No") BRHAA(5,2) 
Write(90,"(e16.8)",advance="No") BRHAA(5,3) 
Write(90,"(e16.8)",advance="No") BRHAA(5,4) 
Write(90,"(e16.8)",advance="No") BRHAA(5,5) 
Write(90,"(e16.8)",advance="No") BRHAA(5,6) 
Write(90,"(e16.8)",advance="No") BRHAA(5,7) 
Write(90,"(e16.8)",advance="No") BRHAA(5,8) 
Write(90,"(e16.8)",advance="No") BRHHH(6,1) 
Write(90,"(e16.8)",advance="No") BRHHH(6,2) 
Write(90,"(e16.8)",advance="No") BRHHH(6,3) 
Write(90,"(e16.8)",advance="No") BRHHH(6,4) 
Write(90,"(e16.8)",advance="No") BRHHH(6,5) 
Write(90,"(e16.8)",advance="No") BRHHH(6,7) 
Write(90,"(e16.8)",advance="No") BRHHH(6,8) 
Write(90,"(e16.8)",advance="No") BRHAA(6,2) 
Write(90,"(e16.8)",advance="No") BRHAA(6,3) 
Write(90,"(e16.8)",advance="No") BRHAA(6,4) 
Write(90,"(e16.8)",advance="No") BRHAA(6,5) 
Write(90,"(e16.8)",advance="No") BRHAA(6,6) 
Write(90,"(e16.8)",advance="No") BRHAA(6,7) 
Write(90,"(e16.8)",advance="No") BRHAA(6,8) 
Write(90,"(e16.8)",advance="No") BRHHH(7,1) 
Write(90,"(e16.8)",advance="No") BRHHH(7,2) 
Write(90,"(e16.8)",advance="No") BRHHH(7,3) 
Write(90,"(e16.8)",advance="No") BRHHH(7,4) 
Write(90,"(e16.8)",advance="No") BRHHH(7,5) 
Write(90,"(e16.8)",advance="No") BRHHH(7,6) 
Write(90,"(e16.8)",advance="No") BRHHH(7,8) 
Write(90,"(e16.8)",advance="No") BRHAA(7,2) 
Write(90,"(e16.8)",advance="No") BRHAA(7,3) 
Write(90,"(e16.8)",advance="No") BRHAA(7,4) 
Write(90,"(e16.8)",advance="No") BRHAA(7,5) 
Write(90,"(e16.8)",advance="No") BRHAA(7,6) 
Write(90,"(e16.8)",advance="No") BRHAA(7,7) 
Write(90,"(e16.8)",advance="No") BRHAA(7,8) 
Write(90,"(e16.8)",advance="No") BRHHH(8,1) 
Write(90,"(e16.8)",advance="No") BRHHH(8,2) 
Write(90,"(e16.8)",advance="No") BRHHH(8,3) 
Write(90,"(e16.8)",advance="No") BRHHH(8,4) 
Write(90,"(e16.8)",advance="No") BRHHH(8,5) 
Write(90,"(e16.8)",advance="No") BRHHH(8,6) 
Write(90,"(e16.8)",advance="No") BRHHH(8,7) 
Write(90,"(e16.8)",advance="No") BRHAA(8,2) 
Write(90,"(e16.8)",advance="No") BRHAA(8,3) 
Write(90,"(e16.8)",advance="No") BRHAA(8,4) 
Write(90,"(e16.8)",advance="No") BRHAA(8,5) 
Write(90,"(e16.8)",advance="No") BRHAA(8,6) 
Write(90,"(e16.8)",advance="No") BRHAA(8,7) 
Write(90,"(e16.8)",advance="No") BRHAA(8,8) 
Write(90,"(e16.8)",advance="No") BRAHH(2,1) 
Write(90,"(e16.8)",advance="No") BRAHH(2,2) 
Write(90,"(e16.8)",advance="No") BRAHH(2,3) 
Write(90,"(e16.8)",advance="No") BRAHH(2,4) 
Write(90,"(e16.8)",advance="No") BRAHH(2,5) 
Write(90,"(e16.8)",advance="No") BRAHH(2,6) 
Write(90,"(e16.8)",advance="No") BRAHH(2,7) 
Write(90,"(e16.8)",advance="No") BRAHH(2,8) 
Write(90,"(e16.8)",advance="No") BRAAA(2,3) 
Write(90,"(e16.8)",advance="No") BRAAA(2,4) 
Write(90,"(e16.8)",advance="No") BRAAA(2,5) 
Write(90,"(e16.8)",advance="No") BRAAA(2,6) 
Write(90,"(e16.8)",advance="No") BRAAA(2,7) 
Write(90,"(e16.8)",advance="No") BRAAA(2,8) 
Write(90,"(e16.8)",advance="No") BRAHH(3,1) 
Write(90,"(e16.8)",advance="No") BRAHH(3,2) 
Write(90,"(e16.8)",advance="No") BRAHH(3,3) 
Write(90,"(e16.8)",advance="No") BRAHH(3,4) 
Write(90,"(e16.8)",advance="No") BRAHH(3,5) 
Write(90,"(e16.8)",advance="No") BRAHH(3,6) 
Write(90,"(e16.8)",advance="No") BRAHH(3,7) 
Write(90,"(e16.8)",advance="No") BRAHH(3,8) 
Write(90,"(e16.8)",advance="No") BRAAA(3,2) 
Write(90,"(e16.8)",advance="No") BRAAA(3,4) 
Write(90,"(e16.8)",advance="No") BRAAA(3,5) 
Write(90,"(e16.8)",advance="No") BRAAA(3,6) 
Write(90,"(e16.8)",advance="No") BRAAA(3,7) 
Write(90,"(e16.8)",advance="No") BRAAA(3,8) 
Write(90,"(e16.8)",advance="No") BRAHH(4,1) 
Write(90,"(e16.8)",advance="No") BRAHH(4,2) 
Write(90,"(e16.8)",advance="No") BRAHH(4,3) 
Write(90,"(e16.8)",advance="No") BRAHH(4,4) 
Write(90,"(e16.8)",advance="No") BRAHH(4,5) 
Write(90,"(e16.8)",advance="No") BRAHH(4,6) 
Write(90,"(e16.8)",advance="No") BRAHH(4,7) 
Write(90,"(e16.8)",advance="No") BRAHH(4,8) 
Write(90,"(e16.8)",advance="No") BRAAA(4,2) 
Write(90,"(e16.8)",advance="No") BRAAA(4,3) 
Write(90,"(e16.8)",advance="No") BRAAA(4,5) 
Write(90,"(e16.8)",advance="No") BRAAA(4,6) 
Write(90,"(e16.8)",advance="No") BRAAA(4,7) 
Write(90,"(e16.8)",advance="No") BRAAA(4,8) 
Write(90,"(e16.8)",advance="No") BRAHH(5,1) 
Write(90,"(e16.8)",advance="No") BRAHH(5,2) 
Write(90,"(e16.8)",advance="No") BRAHH(5,3) 
Write(90,"(e16.8)",advance="No") BRAHH(5,4) 
Write(90,"(e16.8)",advance="No") BRAHH(5,5) 
Write(90,"(e16.8)",advance="No") BRAHH(5,6) 
Write(90,"(e16.8)",advance="No") BRAHH(5,7) 
Write(90,"(e16.8)",advance="No") BRAHH(5,8) 
Write(90,"(e16.8)",advance="No") BRAAA(5,2) 
Write(90,"(e16.8)",advance="No") BRAAA(5,3) 
Write(90,"(e16.8)",advance="No") BRAAA(5,4) 
Write(90,"(e16.8)",advance="No") BRAAA(5,6) 
Write(90,"(e16.8)",advance="No") BRAAA(5,7) 
Write(90,"(e16.8)",advance="No") BRAAA(5,8) 
Write(90,"(e16.8)",advance="No") BRAHH(6,1) 
Write(90,"(e16.8)",advance="No") BRAHH(6,2) 
Write(90,"(e16.8)",advance="No") BRAHH(6,3) 
Write(90,"(e16.8)",advance="No") BRAHH(6,4) 
Write(90,"(e16.8)",advance="No") BRAHH(6,5) 
Write(90,"(e16.8)",advance="No") BRAHH(6,6) 
Write(90,"(e16.8)",advance="No") BRAHH(6,7) 
Write(90,"(e16.8)",advance="No") BRAHH(6,8) 
Write(90,"(e16.8)",advance="No") BRAAA(6,2) 
Write(90,"(e16.8)",advance="No") BRAAA(6,3) 
Write(90,"(e16.8)",advance="No") BRAAA(6,4) 
Write(90,"(e16.8)",advance="No") BRAAA(6,5) 
Write(90,"(e16.8)",advance="No") BRAAA(6,7) 
Write(90,"(e16.8)",advance="No") BRAAA(6,8) 
Write(90,"(e16.8)",advance="No") BRAHH(7,1) 
Write(90,"(e16.8)",advance="No") BRAHH(7,2) 
Write(90,"(e16.8)",advance="No") BRAHH(7,3) 
Write(90,"(e16.8)",advance="No") BRAHH(7,4) 
Write(90,"(e16.8)",advance="No") BRAHH(7,5) 
Write(90,"(e16.8)",advance="No") BRAHH(7,6) 
Write(90,"(e16.8)",advance="No") BRAHH(7,7) 
Write(90,"(e16.8)",advance="No") BRAHH(7,8) 
Write(90,"(e16.8)",advance="No") BRAAA(7,2) 
Write(90,"(e16.8)",advance="No") BRAAA(7,3) 
Write(90,"(e16.8)",advance="No") BRAAA(7,4) 
Write(90,"(e16.8)",advance="No") BRAAA(7,5) 
Write(90,"(e16.8)",advance="No") BRAAA(7,6) 
Write(90,"(e16.8)",advance="No") BRAAA(7,8) 
Write(90,"(e16.8)",advance="No") BRAHH(8,1) 
Write(90,"(e16.8)",advance="No") BRAHH(8,2) 
Write(90,"(e16.8)",advance="No") BRAHH(8,3) 
Write(90,"(e16.8)",advance="No") BRAHH(8,4) 
Write(90,"(e16.8)",advance="No") BRAHH(8,5) 
Write(90,"(e16.8)",advance="No") BRAHH(8,6) 
Write(90,"(e16.8)",advance="No") BRAHH(8,7) 
Write(90,"(e16.8)",advance="No") BRAHH(8,8) 
Write(90,"(e16.8)",advance="No") BRAAA(8,2) 
Write(90,"(e16.8)",advance="No") BRAAA(8,3) 
Write(90,"(e16.8)",advance="No") BRAAA(8,4) 
Write(90,"(e16.8)",advance="No") BRAAA(8,5) 
Write(90,"(e16.8)",advance="No") BRAAA(8,6) 
Write(90,"(e16.8)",advance="No") BRAAA(8,7) 

 
 
Write(92,"(I1)",advance="No") 1
Write(92,"(e16.8)",advance="No") BR_tWb 
Write(92,"(e16.8)",advance="No") BR_tHb(2) 
Write(92,"(e16.8)",advance="No") BR_tHb(3) 
Write(92,"(e16.8)",advance="No") BR_tHb(4) 
Write(92,"(e16.8)",advance="No") BR_tHb(5) 
Write(92,"(e16.8)",advance="No") BR_tHb(6) 
Write(92,"(e16.8)",advance="No") BR_tHb(7) 
Write(92,"(e16.8)",advance="No") BR_tHb(8) 
Write(91,"(I1)",advance="No") 1 
Write(91,"(3e16.8)",advance="No") BR_Hcs(2) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(2) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(2) 
Write(91,"(3e16.8)",advance="No") BR_Hcs(3) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(3) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(3) 
Write(91,"(3e16.8)",advance="No") BR_Hcs(4) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(4) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(4) 
Write(91,"(3e16.8)",advance="No") BR_Hcs(5) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(5) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(5) 
Write(91,"(3e16.8)",advance="No") BR_Hcs(6) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(6) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(6) 
Write(91,"(3e16.8)",advance="No") BR_Hcs(7) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(7) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(7) 
Write(91,"(3e16.8)",advance="No") BR_Hcs(8) 
Write(91,"(3e16.8)",advance="No")BR_Hcb(8) 
Write(91,"(3e16.8)",advance="No")BR_Htaunu(8) 
Write(89,"(I1)",advance="No") 1 
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(1,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(2,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(3,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(4,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(5,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(6,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(7,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(8,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(2,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(3,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(4,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(5,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(6,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(7,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(8,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(1,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(2,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(3,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(4,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(5,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(6,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(7,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(8,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(2,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(3,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(4,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(5,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(6,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(7,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(8,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(1,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(2,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(3,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(4,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(5,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(6,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(7,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(8,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(2,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(3,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(4,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(5,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(6,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(7,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(8,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(1,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(2,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(3,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(4,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(5,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(6,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(7,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(8,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(2,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(3,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(4,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(5,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(6,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(7,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(8,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(1,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(2,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(3,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(4,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(5,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(6,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(7,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fd(8,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(2,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(3,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(4,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(5,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(6,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(7,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fd(8,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(1,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(2,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(3,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(4,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(5,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(6,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(7,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fd(8,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(2,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(3,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(4,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(5,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(6,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(7,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fd(8,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(1,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(2,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(3,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(4,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(5,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(6,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(7,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Fu(8,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(2,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(3,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(4,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(5,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(6,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(7,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Fu(8,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(1,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(2,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(3,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(4,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(5,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(6,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(7,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Fu(8,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(2,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(3,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(4,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(5,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(6,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(7,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Fu(8,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(1,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(2,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(3,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(4,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(5,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(6,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(7,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(8,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(2,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(3,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(4,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(5,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(6,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(7,2)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(8,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(1,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(2,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(3,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(4,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(5,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(6,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(7,2)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(8,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(2,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(3,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(4,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(5,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(6,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(7,2)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(8,2)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(1,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(2,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(3,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(4,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(5,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(6,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(7,3)
Write(89,"(3e16.8)",advance="No") rHB_S_S_Cha(8,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(2,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(3,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(4,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(5,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(6,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(7,3)
Write(89,"(3e16.8)",advance="No") rHB_P_S_Cha(8,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(1,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(2,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(3,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(4,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(5,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(6,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(7,3)
Write(89,"(3e16.8)",advance="No") rHB_S_P_Cha(8,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(2,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(3,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(4,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(5,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(6,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(7,3)
Write(89,"(3e16.8)",advance="No") rHB_P_P_Cha(8,3)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(1)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(2)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(3)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(4)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(5)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(6)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(7)
Write(89,"(3e16.8)",advance="No") rHB_S_VWm(8)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(2)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(3)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(4)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(5)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(6)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(7)
Write(89,"(3e16.8)",advance="No") rHB_P_VWm(8)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(1)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(2)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(3)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(4)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(5)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(6)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(7)
Write(89,"(3e16.8)",advance="No") rHB_S_VZ(8)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(2)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(3)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(4)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(5)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(6)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(7)
Write(89,"(3e16.8)",advance="No") rHB_P_VZ(8)
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(1),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(2),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(3),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(4),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(5),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(6),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(7),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioPP(8),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(2),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(3),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(4),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(5),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(6),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(7),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPPP(8),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(1),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(2),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(3),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(4),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(5),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(6),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(7),dp) 
Write(89,"(3e16.8)",advance="No")  Real(ratioGG(8),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(2),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(3),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(4),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(5),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(6),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(7),dp) 
Write(89,"(3e16.8)",advance="No") Real(ratioPGG(8),dp) 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(3e16.8)",advance="No") 0._dp 
Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(3,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(3,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(4,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(3,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(4,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(5,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(3,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(4,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(5,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(6,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(3,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(4,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(5,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(6,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(7,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(1,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(2,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(3,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(4,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(5,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(6,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(7,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_H_H_Z(8,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(2,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(2,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(3,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(3,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(3,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(4,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(4,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(4,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(4,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(5,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(5,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(5,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(5,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(5,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(6,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(6,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(6,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(6,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(6,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(6,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(7,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(7,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(7,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(7,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(7,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(7,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(7,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,1), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_H_Z(8,8), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,2), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,3), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,4), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,5), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,6), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,7), dp) 
 Write(89,"(e16.8)",advance="No") Real(CPL_A_A_Z(8,8), dp) 
 
 
 
Write(93,"(I1)",advance="No") 1 
Write(93,"(e16.8)",advance="No") 0.0000 
Write(93,"(e16.8)",advance="No") 0.0000 
Write(93,"(e16.8)",advance="No") 0.0000 
Write(93,"(e16.8)",advance="No") 0.0000 
Write(93,"(e16.8)",advance="No") 0.0000 
Write(93,"(e16.8)",advance="No") 0.0000 
Write(93,"(e16.8)",advance="No") 0.0000 
Close(87) 
Close(88) 
Close(90) 
Close(91) 
Close(92) 
Close(93) 
End Subroutine WriteHiggsBounds 
 
 
 Subroutine ReadMatrixC(io, nmax1, nmax2, mat, ic, mat_name, kont, fill)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax1, nmax2, io, ic
   Integer, Intent(in), Optional :: fill
   Complex(dp), Intent(inout) :: mat(nmax1, nmax2)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax1)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax2)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If

    If (ic.Eq.0) Then
     mat(i1,i2) = Cmplx(0._dp,Aimag(mat(i1,i2)),dp) + wert
     If (Present(fill).And.(i1.Ne.i2)) &
       &  mat(i2,i1) = Cmplx(0._dp, Aimag(mat(i2,i1)), dp) + wert
    Else If (ic.Eq.1) Then
     mat(i1,i2) = Real(mat(i1,i2),dp) + Cmplx(0._dp, wert, dp)
     !-------------------------------------------------------------
     ! if fill==1 -> matrix is hermitian
     ! if fill==2 -> matrix is complex symmetric
     !-------------------------------------------------------------
     If (Present(fill).And.(i1.Ne.i2)) Then
      If (fill.Eq.1) mat(i2,i1) = Real(mat(i2,i1),dp) - Cmplx(0._dp, wert, dp)
      If (fill.Eq.2) mat(i2,i1) = Real(mat(i2,i1),dp) + Cmplx(0._dp, wert, dp)
     End If
    End If

   End Do

   200 Return

  End Subroutine ReadMatrixC

 
  Subroutine ReadMatrixR(io, nmax1, nmax2, mat, mat_name, kont)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax1, nmax2, io
   Real(dp), Intent(inout) :: mat(nmax1, nmax2)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax1)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax2)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If

    mat(i1,i2) = wert

   End Do

   200 Return

  End Subroutine ReadMatrixR

  
  Subroutine ReadVectorC(io, nmax, vec, ic, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io, ic
   Complex(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If

    If (ic.Eq.0) vec(i1) = Cmplx(0._dp, Aimag(vec(i1)), dp) + wert
    If (ic.Eq.1) vec(i1) = Real(vec(i1),dp) + Cmplx(0._dp, wert, dp)

   End Do

   200 Return

  End Subroutine ReadVectorC


Subroutine ReadScalarC(io, vec, ic, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: io, ic
   Complex(dp), Intent(inout) :: vec
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) wert, read_line

!     If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
!      Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
!         & Trim(NameOfUnit(Iname))//", index i1=",i1
!      Iname = Iname - 2
!      kont = -305 
!      Call TerminateProgram()
!     End If

    If (ic.Eq.0) vec = Cmplx(0._dp, Aimag(vec), dp) + wert
    If (ic.Eq.1) vec = Real(vec,dp) + Cmplx(0._dp, wert, dp)

   End Do

   200 Return

  End Subroutine ReadScalarC

  
  Subroutine ReadVectorR(io, nmax, vec, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io
   Real(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If

    vec(i1) = wert

   End Do

   200 Return

  End Subroutine ReadVectorR

Subroutine ReadScalarR(io, vec, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: io
   Real(dp), Intent(inout) :: vec
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) wert, read_line

!     If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
!      Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
!         & Trim(NameOfUnit(Iname))//", index i1=",i1
!      Iname = Iname - 2
!      kont = -305 
!      Call TerminateProgram()
!     End If

    vec = wert

   End Do

   200 Return

  End Subroutine ReadScalarR

  
  Subroutine ReadTensorC(io, nmax1, nmax2, nmax3, mat, ic, mat_name, kont)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax1, nmax2, nmax3, io, ic
   Complex(dp), Intent(inout) :: mat(nmax1, nmax2, nmax3)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2, i3
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadTensorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, i3, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax1)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax2)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If
    If ((i3.Lt.1).Or.(i3.Gt.nmax3)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i3=",i3
     Iname = Iname - 2
     kont = -305 
     Call TerminateProgram()
    End If

    If (ic.Eq.0) mat(i1,i2,i3) = Cmplx(0._dp, Aimag(mat(i1,i2,i3)), dp) + wert
    If (ic.Eq.1) mat(i1,i2,i3) = mat(i1,i2,i3) + Cmplx(0._dp, wert, dp)

   End Do

   200 Return

  End Subroutine ReadTensorC

 Subroutine SetWriteMinBR(wert)
 !-------------------------------------------------------------------
 ! sets the minimal branching ratio (=wert) appearing in the output
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: wert
  BrMin = wert
 End Subroutine SetWriteMinBR


 Subroutine SetWriteMinSig(wert)
 !-------------------------------------------------------------------
 ! sets the minimal cross section (=wert) appearing in the output
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: wert
  SigMin = wert
 End Subroutine SetWriteMinSig

 Subroutine Warn_CPV(i_cpv, name)
  Implicit None 
   Integer, Intent(in) :: i_cpv
   Character(len=*), Intent(in) :: name
   If (i_cpv.Eq.0) Write(ErrCan,*) "CP violation is switched off"
   If (i_cpv.Eq.1) Write(ErrCan,*) "CP violation beyond CKM is switched off"
   Write(ErrCan,*) "Ignoring block "//Trim(name)
   If (ErrorLevel.Eq.2) Call TerminateProgram
  End Subroutine Warn_CPV

Subroutine PutUpperCase(name)
 Implicit None
  Character(len=80), Intent(inout) :: name
  Integer :: len=80, i1
  Do i1=1,len
   If (name(i1:i1).Eq."a") name(i1:i1) = "A"
   If (name(i1:i1).Eq."b") name(i1:i1) = "B"
   If (name(i1:i1).Eq."c") name(i1:i1) = "C"
   If (name(i1:i1).Eq."d") name(i1:i1) = "D"
   If (name(i1:i1).Eq."e") name(i1:i1) = "E"
   If (name(i1:i1).Eq."f") name(i1:i1) = "F"
   If (name(i1:i1).Eq."g") name(i1:i1) = "G"
   If (name(i1:i1).Eq."h") name(i1:i1) = "H"
   If (name(i1:i1).Eq."i") name(i1:i1) = "I"
   If (name(i1:i1).Eq."j") name(i1:i1) = "J"
   If (name(i1:i1).Eq."k") name(i1:i1) = "K"
   If (name(i1:i1).Eq."l") name(i1:i1) = "L"
   If (name(i1:i1).Eq."m") name(i1:i1) = "M"
   If (name(i1:i1).Eq."n") name(i1:i1) = "N"
   If (name(i1:i1).Eq."o") name(i1:i1) = "O"
   If (name(i1:i1).Eq."p") name(i1:i1) = "P"
   If (name(i1:i1).Eq."q") name(i1:i1) = "Q"
   If (name(i1:i1).Eq."r") name(i1:i1) = "R"
   If (name(i1:i1).Eq."s") name(i1:i1) = "S"
   If (name(i1:i1).Eq."t") name(i1:i1) = "T"
   If (name(i1:i1).Eq."u") name(i1:i1) = "U"
   If (name(i1:i1).Eq."v") name(i1:i1) = "V"
   If (name(i1:i1).Eq."w") name(i1:i1) = "W"
   If (name(i1:i1).Eq."x") name(i1:i1) = "X"
   If (name(i1:i1).Eq."y") name(i1:i1) = "Y"
   If (name(i1:i1).Eq."z") name(i1:i1) = "Z"
  End Do
 End Subroutine PutUpperCase

Subroutine Read_FLIFE(io) 
Implicit None 
Integer,Intent(in)::io 
Real(dp)::r_mod,wert 
Integer::i_mod,i_test,i_rp 
Character(len=80)::read_line 
Do 
Read(io,*) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_test,wert!,read_line 
   if (i_test.Eq.111) Then 
    tau_pi0 =wert 
  Else if (i_test.Eq.211) Then 
    tau_pip =wert 
  Else if (i_test.Eq.113) Then 
    tau_rho0 =wert 
  Else if (i_test.Eq.421) Then 
    tau_D0 =wert 
  Else if (i_test.Eq.411) Then 
    tau_Dp =wert 
  Else if (i_test.Eq.431) Then 
    tau_DSp =wert 
  Else if (i_test.Eq.433) Then 
    tau_DSsp =wert 
  Else if (i_test.Eq.221) Then 
    tau_eta =wert 
  Else if (i_test.Eq.331) Then 
    tau_etap =wert 
  Else if (i_test.Eq.223) Then 
    tau_omega =wert 
  Else if (i_test.Eq.333) Then 
    tau_phi =wert 
  Else if (i_test.Eq.130) Then 
    tau_KL0 =wert 
  Else if (i_test.Eq.310) Then 
    tau_KS0 =wert 
  Else if (i_test.Eq.311) Then 
    tau_K0 =wert 
  Else if (i_test.Eq.321) Then 
    tau_Kp =wert 
  Else if (i_test.Eq.511) Then 
    tau_B0d =wert 
  Else if (i_test.Eq.531) Then 
    tau_B0s =wert 
  Else if (i_test.Eq.521) Then 
    tau_Bp =wert 
  Else if (i_test.Eq.513) Then 
    tau_B0c =wert 
  Else if (i_test.Eq.523) Then 
    tau_Bpc =wert 
  Else if (i_test.Eq.541) Then 
    tau_Bcp =wert 
  Else if (i_test.Eq.313) Then 
    tau_K0c =wert 
  Else if (i_test.Eq.323) Then 
    tau_Kpc =wert 
  Else if (i_test.Eq.441) Then 
    tau_etac =wert 
  Else if (i_test.Eq.443) Then 
    tau_JPsi =wert 
  Else if (i_test.Eq.553) Then 
    tau_Ups =wert 
End If 
End Do! i_mod 
End Subroutine Read_FLIFE 


Subroutine Read_FMASS(io) 
Implicit None 
Integer,Intent(in)::io 
Real(dp)::r_mod,wert 
Integer::i_mod,i_test,i_rp 
Character(len=80)::read_line 
Do 
Read(io,*) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_test,wert!,read_line 
   if (i_test.Eq.111) Then 
    mass_pi0 =wert 
  Else if (i_test.Eq.211) Then 
    mass_pip =wert 
  Else if (i_test.Eq.113) Then 
    mass_rho0 =wert 
  Else if (i_test.Eq.421) Then 
    mass_D0 =wert 
  Else if (i_test.Eq.411) Then 
    mass_Dp =wert 
  Else if (i_test.Eq.431) Then 
    mass_DSp =wert 
  Else if (i_test.Eq.433) Then 
    mass_DSsp =wert 
  Else if (i_test.Eq.221) Then 
    mass_eta =wert 
  Else if (i_test.Eq.331) Then 
    mass_etap =wert 
  Else if (i_test.Eq.223) Then 
    mass_omega =wert 
  Else if (i_test.Eq.333) Then 
    mass_phi =wert 
  Else if (i_test.Eq.130) Then 
    mass_KL0 =wert 
  Else if (i_test.Eq.310) Then 
    mass_KS0 =wert 
  Else if (i_test.Eq.311) Then 
    mass_K0 =wert 
  Else if (i_test.Eq.321) Then 
    mass_Kp =wert 
  Else if (i_test.Eq.511) Then 
    mass_B0d =wert 
  Else if (i_test.Eq.531) Then 
    mass_B0s =wert 
  Else if (i_test.Eq.521) Then 
    mass_Bp =wert 
  Else if (i_test.Eq.513) Then 
    mass_B0c =wert 
  Else if (i_test.Eq.523) Then 
    mass_Bpc =wert 
  Else if (i_test.Eq.541) Then 
    mass_Bcp =wert 
  Else if (i_test.Eq.313) Then 
    mass_K0c =wert 
  Else if (i_test.Eq.323) Then 
    mass_Kpc =wert 
  Else if (i_test.Eq.441) Then 
    mass_etac =wert 
  Else if (i_test.Eq.443) Then 
    mass_JPSi =wert 
  Else if (i_test.Eq.553) Then 
    mass_Ups =wert 
End If 
End Do! i_mod 
End Subroutine Read_FMASS 


Subroutine Read_FCONST(io) 
Implicit None 
Integer,Intent(in)::io 
Real(dp)::r_mod,wert 
Integer::i_mod,i_test,i_rp 
Character(len=80)::read_line 
Do 
Read(io,*) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_test, i_mod, wert!,read_line 
    If (i_test.Eq.111) Then 
    If (i_mod.Eq.1) Then 
    f_pi_CONST =wert 
    End If 
   Else If (i_test.Eq.213) Then 
    If (i_mod.Eq.1) Then 
    f_rho_CONST =wert 
   Else If (i_mod.Eq.11) Then 
    h_rho_CONST =wert 
    End If 
   Else If (i_test.Eq.221) Then 
    If (i_mod.Eq.1) Then 
    f_eta_q_CONST =wert 
   Else If (i_mod.Eq.2) Then 
    f_eta_s_CONST =wert 
   Else If (i_mod.Eq.11) Then 
    h_eta_q_CONST =wert 
   Else If (i_mod.Eq.12) Then 
    h_eta_s_CONST =wert 
    End If 
   Else If (i_test.Eq.223) Then 
    If (i_mod.Eq.1) Then 
    f_omega_q_CONST =wert 
   Else If (i_mod.Eq.2) Then 
    f_omega_s_CONST =wert 
   Else If (i_mod.Eq.11) Then 
    h_omega_q_CONST =wert 
   Else If (i_mod.Eq.12) Then 
    h_omega_s_CONST =wert 
    End If 
   Else If (i_test.Eq.231) Then 
    If (i_mod.Eq.1) Then 
    f_etap_CONST =wert 
    End If 
   Else If (i_test.Eq.311) Then 
    If (i_mod.Eq.1) Then 
    f_k_CONST =wert 
    End If 
   Else If (i_test.Eq.321) Then 
    If (i_mod.Eq.1) Then 
    f_Kp_CONST =wert 
   Else If (i_mod.Eq.11) Then 
    h_k_CONST =wert 
    End If 
   Else If (i_test.Eq.411) Then 
    If (i_mod.Eq.1) Then 
    f_Dp_CONST =wert 
    End If 
   Else If (i_test.Eq.421) Then 
    If (i_mod.Eq.1) Then 
    f_D_CONST =wert 
    End If 
   Else If (i_test.Eq.431) Then 
    If (i_mod.Eq.1) Then 
    f_Dsp_CONST =wert 
    End If 
   Else If (i_test.Eq.511) Then 
    If (i_mod.Eq.1) Then 
    f_B0d_CONST =wert 
    End If 
   Else If (i_test.Eq.521) Then 
    If (i_mod.Eq.1) Then 
    f_Bp_CONST =wert 
    End If 
   Else If (i_test.Eq.531) Then 
    If (i_mod.Eq.1) Then 
    f_B0s_CONST =wert 
    End If 
End If 
End Do! i_mod 
End Subroutine Read_FCONST 


Subroutine Read_GAUGEIN(io,i_c,i_model,set_mod_par,kont) 
Implicit None 
Integer,Intent(in)::io,i_c,i_model 
Integer,Intent(inout)::kont,set_mod_par(:) 
Integer::i_par 
Real(dp)::wert 
Character(len=80)::read_line 
Do 
Read(io,*,End=200) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line 
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_par,wert!,read_line 
If (i_par.Eq.1) Then 
g1IN= wert 
InputValueforg1= .True. 
Else If (i_par.Eq.2) Then 
g2IN= wert 
InputValueforg2= .True. 
Else If (i_par.Eq.3) Then 
g3IN= wert 
InputValueforg3= .True. 
Else
Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block GAUGEIN ",i_par
If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMGAUGEIN ",i_par
If (i_c.Eq.0) Write(*,*) "Unknown entry for Block GAUGEIN ",i_par
If (i_c.Eq.1) Write(*,*) "Unknown entry for Block IMGAUGEIN ",i_par
Call AddError(304)
If (ErrorLevel.Eq.2) Call TerminateProgram
End If
End Do! i_par
200 Return
End Subroutine Read_GAUGEIN 
 
 
Subroutine Read_MSOFTIN(io,i_c,i_model,set_mod_par,kont) 
Implicit None 
Integer,Intent(in)::io,i_c,i_model 
Integer,Intent(inout)::kont,set_mod_par(:) 
Integer::i_par 
Real(dp)::wert 
Character(len=80)::read_line 
Do 
Read(io,*,End=200) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line 
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_par,wert!,read_line 
If (i_par.Eq.21) Then 
mHd2IN= wert 
InputValueformHd2= .True. 
Else If (i_par.Eq.22) Then 
mHu2IN= wert 
InputValueformHu2= .True. 
Else If (i_par.Eq.1) Then 
If (i_c.Eq.0) M1IN= Cmplx(wert,Aimag(M1IN),dp) 
If (i_c.Eq.1) M1IN= Cmplx(Real(M1IN,dp),wert,dp) 
InputValueforM1= .True. 
Else If (i_par.Eq.2) Then 
If (i_c.Eq.0) M2IN= Cmplx(wert,Aimag(M2IN),dp) 
If (i_c.Eq.1) M2IN= Cmplx(Real(M2IN,dp),wert,dp) 
InputValueforM2= .True. 
Else If (i_par.Eq.3) Then 
If (i_c.Eq.0) M3IN= Cmplx(wert,Aimag(M3IN),dp) 
If (i_c.Eq.1) M3IN= Cmplx(Real(M3IN,dp),wert,dp) 
InputValueforM3= .True. 
Else
Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block MSOFTIN ",i_par
If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMMSOFTIN ",i_par
If (i_c.Eq.0) Write(*,*) "Unknown entry for Block MSOFTIN ",i_par
If (i_c.Eq.1) Write(*,*) "Unknown entry for Block IMMSOFTIN ",i_par
Call AddError(304)
If (ErrorLevel.Eq.2) Call TerminateProgram
End If
End Do! i_par
200 Return
End Subroutine Read_MSOFTIN 
 
 
Subroutine Read_HMIXIN(io,i_c,i_model,set_mod_par,kont) 
Implicit None 
Integer,Intent(in)::io,i_c,i_model 
Integer,Intent(inout)::kont,set_mod_par(:) 
Integer::i_par 
Real(dp)::wert 
Character(len=80)::read_line 
Do 
Read(io,*,End=200) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line 
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_par,wert!,read_line 
If (i_par.Eq.102) Then 
vdIN= wert 
Else If (i_par.Eq.103) Then 
vuIN= wert 
Else
Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block HMIXIN ",i_par
If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMHMIXIN ",i_par
If (i_c.Eq.0) Write(*,*) "Unknown entry for Block HMIXIN ",i_par
If (i_c.Eq.1) Write(*,*) "Unknown entry for Block IMHMIXIN ",i_par
Call AddError(304)
If (ErrorLevel.Eq.2) Call TerminateProgram
End If
End Do! i_par
200 Return
End Subroutine Read_HMIXIN 
 
 
Subroutine Read_PHASESIN(io,i_c,i_model,set_mod_par,kont) 
Implicit None 
Integer,Intent(in)::io,i_c,i_model 
Integer,Intent(inout)::kont,set_mod_par(:) 
Integer::i_par 
Real(dp)::wert 
Character(len=80)::read_line 
Do 
Read(io,*,End=200) read_line 
If (read_line(1:1).Eq."#") Cycle! this loop 
Backspace(io)! resetting to the beginning of the line 
If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b")) Exit! this loop 
Read(io,*) i_par,wert!,read_line 
End Do! i_par
200 Return
End Subroutine Read_PHASESIN 
 
 
Subroutine Switch_to_superCKM(Y_d, Y_u, Ad_in, Au_in, MD_in, MQ_in, MU_in &
                      &, Ad_out, Au_out, MD_out, MQ_out, MU_out, tr        &
                      &, RSd_in, RSu_in, RSd_out, RSu_out, CKM_out, Yd_out, Yu_out )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the electroweak basis to the super CKM basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Au_in, Ad_in, MD_in &
        & , MQ_in, MU_in
  Complex(dp), Optional, Intent(in), Dimension(6,6) :: RSu_in, RSd_in
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Au_out, Ad_out, MD_out, MQ_out &
        & , MU_out, Yd_out, Yu_out
  Complex(dp), Optional, Intent(out), Dimension(6,6) :: RSu_out, RSd_out
  Complex(dp), Optional, Intent(out) :: CKM_out(3,3)

  Complex(dp), Dimension(3,3) :: uU_L, uU_R, uD_L, uD_R, CKM_Q
  Complex(dp) :: rot(6,6), Ephi

  Real(dp) :: mf(3), s12, s23, aR, aI, s13, c13
  Integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------

   Call FermionMass(Y_u, 1._dp, mf, uU_L, uU_R, ierr)
   Call FermionMass(Y_d, 1._dp, mf, uD_L, uD_R, ierr)
   Yu_out = MatMul(MatMul(conjg(uU_L),Y_u),Transpose(conjg(uU_R)))
   Yd_out = MatMul(MatMul(conjg(uD_L),Y_d),Transpose(conjg(uD_R)))

  !---------------------------------------------------------
  ! CKM matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  CKM_Q =  Matmul(uU_R, Transpose(Conjg(ud_R)) )
  uD_L(1,:) = uD_L(1,:) / Conjg(CKM_Q(1,1)) * Abs(CKM_Q(1,1))
  uD_L(2,:) = uD_L(2,:) / Conjg(CKM_Q(1,2)) * Abs(CKM_Q(1,2))
  uU_L(2,:) = uU_L(2,:) / CKM_Q(2,3) * Abs(CKM_Q(2,3))
  uU_L(3,:) = uU_L(3,:) / CKM_Q(3,3) * Abs(CKM_Q(3,3))
  !-------------------------------------------------------------------
  ! also the right quark must be multiplied with the conjugate phase
  ! as otherwise the masses get complex
  !-------------------------------------------------------------------
  uD_R(1,:) = uD_R(1,:) / CKM_Q(1,1) * Abs(CKM_Q(1,1))
  uD_R(2,:) = uD_R(2,:) / CKM_Q(1,2) * Abs(CKM_Q(1,2))
  uU_R(2,:) = uU_R(2,:) * Abs(CKM_Q(2,3)) / Conjg(CKM_Q(2,3))
  uU_R(3,:) = uU_R(3,:) * Abs(CKM_Q(3,3)) / Conjg(CKM_Q(3,3))
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  !--------------------------------------------------------------
  ! one more freedom left
  !--------------------------------------------------------------
  s13 = Abs(CKM_Q(1,3))
  c13 = sqrt(1._dp - s13**2)
  s23 = Abs(CKM_Q(2,3))/c13
  s12 = Abs(CKM_Q(1,2))/c13

  aR = Real(CKM_Q(2,2),dp) + s12 * s23 * Real(CKM_Q(1,3),dp)
  aI =  s12 * s23 * Aimag(CKM_Q(1,3)) - Aimag(CKM_Q(2,2))
  Ephi = Cmplx(aR/Sqrt(aR**2+aI**2),aI/Sqrt(aR**2+aI**2),dp)

  uU_L(2:3,:) = Ephi * uU_L(2:3,:)
  uD_L(3,:) = Ephi * uD_L(3,:)
  Ephi = Conjg(Ephi)
  uU_R(2:3,:) = Ephi * uU_R(2:3,:)
  uD_R(3,:) = Ephi * uD_R(3,:)


  CKM_Q =  Matmul(uU_R, Transpose(Conjg(ud_R)) )

  If (Present(CKM_out)) CKM_out = CKM_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super CKM basis
  !-------------------------------------------------------------------

   Au_out = Matmul( Matmul(Conjg(uU_L), Au_in), Conjg(Transpose(uU_R)))

   Ad_out = Matmul( Matmul(Conjg(uD_L), Ad_in), Conjg(Transpose(uD_R)))


  MD_out = Matmul( Matmul( Conjg(uD_L), Transpose(MD_in)), Transpose(uD_L))
  MU_out = Matmul( Matmul( Conjg(uU_L), Transpose( MU_in)), Transpose(uU_L))
  MQ_out = Matmul( Matmul( uD_R, MQ_in), Transpose(Conjg(uD_R)) )

   If (Present(RSu_in).And.Present(RSu_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uU_L)
    rot(4:6,4:6) = uU_R
    RSu_out = Matmul(RSu_in,Transpose(rot))
   End If
   If (Present(RSd_in).And.Present(RSd_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uD_L)
    rot(4:6,4:6) = uD_R
    RSd_out = Matmul(RSd_in,Transpose(rot))
   End If

 End Subroutine Switch_to_superCKM




 Subroutine Switch_to_superPMNS(Y_l, uN_L, Al_in, ME_in, ML_in &
                      &, Al_out, ME_out, ML_out, tr            &
                      &, RSl_in, RSn_in, RSl_out, RSn_out, PMNS_out, Yl )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the electroweak basis to the super PMNS basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in), Dimension(3,3) :: Y_l, uN_L, Al_in, ME_in, ML_in
  Complex(dp), Optional, Intent(in) :: RSl_in(6,6), RSn_in(3,3)
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Al_out, ME_out, ML_out
  Complex(dp), Optional, Intent(out), Dimension(6,6) :: RSl_out(6,6)
  Complex(dp), Optional, Intent(out) :: PMNS_out(3,3), RSn_out(3,3)
  Complex(dp), Optional, Intent(out) :: Yl(3,3)

  Complex(dp), Dimension(3,3) :: uL_L, uL_R, PMNS_Q
  Complex(dp) :: rot(6,6)

  Real(dp) :: mf(3)
  Integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------



   If (tr) Then
   Call FermionMass(Transpose(Y_l), 1._dp, mf, uL_L, uL_R, ierr)
  Else
   Call FermionMass(Y_l, 1._dp, mf, uL_L, uL_R, ierr)
  End If

  If (Present(Yl)) Then 
    Yl = 0._dp
    Yl(1,1) = sqrt2 * mf(1)
    Yl(2,2) = sqrt2 * mf(2)
    Yl(3,3) = sqrt2 * mf(3)
  End if

  !---------------------------------------------------------
  ! PMNS matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  PMNS_Q =  Matmul(uL_L, uN_L)

  If (Present(PMNS_out)) PMNS_out = PMNS_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super PMNS basis
  !-------------------------------------------------------------------
  If (tr) Then
   Al_out = Matmul( Matmul(uL_R, Al_in), Transpose(Conjg(uL_L)))

   ME_out = Matmul( Matmul( uL_R, ME_in), Transpose(Conjg(uL_R)))
   ML_out = Matmul( Matmul( Transpose(uL_L), ML_in), Conjg(uL_L) )

  Else
   Al_out = Matmul( Matmul(Conjg(uL_L), Al_in), Transpose(uL_R))

   ME_out = Matmul( Matmul( Conjg(uL_R), ME_in), Transpose(uL_R))
   ML_out = Matmul( Matmul( uL_L, ML_in), Transpose(Conjg(uL_L)) )

  End If
  !------------------------------------------------------------------
  ! to avoid numerical problems ensure that matrices are hermitian
  !-----------------------------------------------------------------
  ME_out = 0.5_dp * ( ME_out + Conjg(Transpose(ME_out)) )
  ML_out = 0.5_dp * ( ML_out + Conjg(Transpose(ML_out)) )

   If (Present(RSn_in).And.Present(RSn_out)) Then
    RSn_out = Matmul(RSn_in, Conjg(uN_L) )
   End If
   If (Present(RSl_in).And.Present(RSl_out)) Then
    rot = 0._dp
    rot(1:3,1:3) = Transpose(uL_L)
    rot(4:6,4:6) = Transpose(Conjg(uL_R))
    RSl_out = Matmul(RSl_in, rot)
   End If

 End Subroutine Switch_to_superPMNS

Subroutine SLHA1converter(MSd,MSd2, MSu,MSu2, MSe, MSe2, MSv, MSv2, &
   & ZD,ZU,ZE,ZV,Ztop,Zbottom,Ztau,  &
   & PDGd, PDGu, PDGe, PDGv, NamesD, NamesU, NamesE, NamesV)
Implicit None
Real(dp), Intent(inout) :: MSd(6), MSu(6), MSe(6), MSv(3), MSd2(6), MSu2(6), MSe2(6), MSv2(3)
Character(len=30),Dimension(6), Intent(inout) :: NamesD, NamesU, NamesE
Character(len=30),Dimension(3), Intent(inout) :: NamesV
Complex(dp), Intent(inout) :: PDGd(6), PDGu(6), PDGe(6), PDGv(3)
Complex(dp), Intent(in) :: ZU(6,6), ZD(6,6), ZE(6,6), ZV(3,3)
Complex(dp), Intent(out) :: Ztop(2,2), Zbottom(2,2), Ztau(2,2)
Real(dp) :: MSdt(6), MSut(6), MSet(6), MSvt(3)
Character(len=30),Dimension(6) :: NamesDt, NamesUt, NamesEt
Character(len=30),Dimension(3) :: NamesVt
Complex(dp) :: PDGdt(6), PDGut(6), PDGet(6), PDGvt(3)


Integer :: i1, i2, i3, ii, jj, i_zaehl


! Down-Squark

Call CheckMixingMatrix6(ZD, Zbottom,MSd,MSdt, PDGd, PDGdt, NamesD, NamesDt, 1)
!GammaD = GammaDt
!BrD = BrDt
PDGd = PDGdt
NamesD = NamesDt
! MSd = MSdt
! MSd2 = MSdt**2
   

! Up-Squark

Call CheckMixingMatrix6(ZU, Ztop,MSu,MSut, PDGu, PDGut, NamesU, NamesUt,2)
!GammaU = GammaUt
!BrU = BrUt
PDGu = PDGut
NamesU = NamesUt
! MSu = MSut
! MSu2 = MSut**2


! Selectron

Call CheckMixingMatrix6(ZE, Ztau,MSe,MSet, PDGe, PDGet, NamesE, NamesEt,3)
!GammaE = GammaEt
PDGe = PDGet
NamesE = NamesEt
!BrE = BrEt
! MSe = MSet
! MSe2 = MSet**2

! Sneutrino

Call CheckMixingMatrix3(ZV, MSv,MSvt, PDGv, PDGvt, NamesV, NamesVt)
!GammaV = GammaVt
PDGv = PDGvt
NamesV = NamesVt
!BrV = BrVt
! MSv = MSvt
! MSv2 = MSvt**2

Contains 

Subroutine CheckMixingMatrix6(Z,Z_out,m_in,m_out, PDG_in, PDG_out, Names_in, Names_out,particle)
Implicit None
Complex(dp), Intent(in) :: Z(6,6)
Complex(dp), Intent(out) :: Z_out(2,2)
Character(len=30),Dimension(6), Intent(in) :: Names_in
Character(len=30),Dimension(6), Intent(out) :: Names_out
Character(len=30) :: Names_temp, Names_save(6)
Complex(dp), Intent(in) :: PDG_in(6)
Complex(dp), Intent(out) :: PDG_out(6)
Complex(dp) :: PDG_temp, PDG_save(6)
Real(dp), Intent(in) :: m_in(6)
Real(dp) :: mat6R(6,6), mtemp, Maxcont, Z_temp(2)
Integer, Intent(in) ::particle
Real(dp), Intent(out) :: m_out(6)
Integer :: i1, ii, jj, PDG3, PDG6

      Select Case(particle)
       Case(1)  ! d-squark
          Names_save(1) = "SdL"
          Names_save(2) = "SsL"
          Names_save(3) = "Sb1" 
          Names_save(4) = "SdR"
          Names_save(5) = "SsR"
          Names_save(6) = "Sb2" 
          PDG_save(1) = 1000001
          PDG_save(2) = 1000003
          PDG_save(3) = 1000005
          PDG_save(4) = 2000001
          PDG_save(5) = 2000003
          PDG_save(6) = 2000005
       Case(2)  ! u-squark
          Names_save(1) = "SuL"
          Names_save(2) = "ScL"
          Names_save(3) = "St1" 
          Names_save(4) = "SuR"
          Names_save(5) = "ScR"
          Names_save(6) = "St2"
          PDG_save(1) = 1000002
          PDG_save(2) = 1000004
          PDG_save(3) = 1000006
          PDG_save(4) = 2000002
          PDG_save(5) = 2000004
          PDG_save(6) = 2000006 
       Case(3)  ! selectron
          Names_save(1) = "SeL"
          Names_save(2) = "SmuL"
          Names_save(3) = "Stau1" 
          Names_save(4) = "SeR"
          Names_save(5) = "SmuR"
          Names_save(6) = "Stau2"
          PDG_save(1) = 1000011
          PDG_save(2) = 1000013
          PDG_save(3) = 1000015
          PDG_save(4) = 2000011
          PDG_save(5) = 2000013
          PDG_save(6) = 2000015 
      End Select


     mat6R = Abs(Z)
     Do i1=1,6
      jj = Maxloc(mat6R(i1,:),1)
      m_out(jj) = m_in(i1)
      PDG_out(i1) = PDG_save(jj)
      Names_out(i1) = Names_save(jj)
      If (jj.eq.3) Then
         Z_out(1,1) = MaxVal(mat6R(i1,:))
         PDG3 = i1
      End if
      If (jj.eq.6) Then
          Z_out(2,2) = MaxVal(mat6R(i1,:)) 
          PDG6 = i1
      End if
!       mat6R(ii,:) = 0._dp
!       mat6R(:,jj) = 0._dp
     End Do
      Z_out(1,2) = sqrt(1._dp - Z_out(1,1)**2)
      Z_out(2,1) = -Z_out(1,2)
     If (M_out(3).gt.M_out(6)) Then

      Names_out(PDG3) = Names_save(6)
      Names_out(PDG6) = Names_save(3) 

      PDG_out(PDG3) = PDG_save(6)
      PDG_out(PDG6) = PDG_save(3) 

      Z_temp = Z_out(1,:)
      Z_out(1,:) = Z_out(2,:)
      Z_out(2,:) = Z_temp
     End if






End Subroutine CheckMixingMatrix6


Subroutine CheckMixingMatrix3(Z,m_in,m_out, PDG_in, PDG_out, Names_in, Names_out)
Implicit None
Complex(dp), Intent(in) :: Z(3,3)
Real(dp), Intent(in) :: m_in(3)
Character(len=30),Dimension(3), Intent(in) :: Names_in
Character(len=30),Dimension(3), Intent(out) :: Names_out
Character(len=30),Dimension(3) :: Names_save
Complex(dp), Intent(in) :: PDG_in(3)
Complex(dp), Intent(out) :: PDG_out(3)
Complex(dp) :: PDG_save(3)
Real(dp), Intent(out) :: m_out(3)
Real(dp) :: mat6R(3,3), mtemp, Maxcont
Integer :: i1, ii, jj

          PDG_save(1) = 1000012
          PDG_save(2) = 1000014
          PDG_save(3) = 1000016
          Names_save(1) = "Snu_e"
          Names_save(2) = "Snu_mu"
          Names_save(3) = "Snu_tau" 

     mat6R = Abs(Z)
     Do i1=1,3
      jj = Maxloc(mat6R(i1,:),1)
      m_out(jj) = m_in(i1)
      PDG_out(i1) = PDG_save(jj)
      Names_out(i1) = Names_save(jj)
!       mat6R(ii,:) = 0._dp
!       mat6R(:,jj) = 0._dp
     End Do

End Subroutine CheckMixingMatrix3


End Subroutine SLHA1converter
End Module InputOutput_munuSSM3G 
 
