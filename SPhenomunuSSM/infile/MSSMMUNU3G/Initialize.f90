Subroutine Initialize(LFlags, NuisP) 

Use parameters

Use Control 
!Use Experiment 
Use Model_Data_munuSSM3G 
Use LoopFunctions 
Use SugraRuns_munuSSM3G 
 
Use EffPotFunctions 
 
Implicit None 

Type(Flags) :: LFlags 
Type(Nuisance_params):: NuisP

Integer :: kont
!Real(dp), Intent(out) :: Fgmsb, Ecms(:), Pm(:), Pp(:)
!Logical, Intent(out) :: l_ISR(:)
Character(len=80) :: read_line
Integer :: i_mod=-1, i_sm=-1, i_par=-1, set_mod_par(25)=-1 &
& , i1, p_max, p_act, i_sp, i_model=-1, i_particles=-1, wert, i_cpv
Real(dp) :: Abs_Mu2, cosb2, cos2b, sinb2, RG0(3,3) &
 & , mat_D(3,3), R2(2,2), s12,s13,s23,c12,c13,c23
Logical :: check, calc_ferm, check_alpha(2)
Complex(dp) :: lam_vS
Logical, Save :: l_open = .True. 
 

    check_alpha = .False. ! used to check consistency of alpha(mZ) calculation
    !in_kont = 0

    Call InitializeStandardModel 
    Call InitializeLoopFunctions

    i_mod = -1
    i_sm = -1
    i_par = -1
    set_mod_par = -1 

    ErrorLevel = -1
    GenerationMixing=.False.

    Call Set_All_Parameters_0()

    lam_vs = 0._dp
    !sp_info = " "
    HighScaleModel="SARAH_Generated_Model" 
    TwoLoopRGE = .True.
    !Fgmsb = 1.e12_dp
    m32 = 1.e20_dp 
 
    kont = 0

    !Assumed Gneration mixing
    RXiNew = 1d0  

    CKMcomplex = CKM 
    i_cpv = 0
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

    HighScaleModel = "LOW"
    BoundaryCondition = 1
    GenerationMixing = .True.
    Call SetRGEScale(1d3**2)  ! set Q_EWSB

    alphaS_mZ = NuisP%alphasmz

    mf_d(3) = NuisP%mbmb
    mf_d2(3) = mf_d(3)**2
    calc_ferm = .True.
 
    mf_u(3) = NuisP%mtop
    mf_u2(3) = mf_u(3)**2



    CalculateTwoLoopHiggsMasses=.False.

    wert = 3
    SELECT CASE (WERT)
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

    !fixed     
    GaugelessLimit=.true.

    !           Case(10)
    !   If (wert.Ne.1) Then
    TwoLoopSafeMode=.false.
    !   Else
    !    TwoLoopSafeMode=.true.
    !   End If

    !fixed ??? May be we need for the Higgs
!    L_BR = .False.
    L_BR = .True.

    !BrMin = 1d-4
    !Call SetWriteMinBR(1d-4)

    !
    Enable3BDecaysF = .False.
    Enable3BDecaysS = .False. 

    !Fixed
    RunningCouplingsDecays = .False. !.True.

    !Fixed
    MinWidth = 1d-30
 
    !Fixed
    Call SetGUTScale(-1d0)

    !Fixed
    delta_mass = 1d-3 !1d-4

    !Fixed
    n_run = 40
    MinimalNumberIterations = 0

    !Fixed
    i1 = SetYukawaScheme(1)

    !Fixed
    TwoLoopRGE=.True.

    RotateNegativeFermionMasses=.True.

    SwitchToSCKM=.False.

    IgnoreNegativeMasses=.False.

    CalculateOneLoopMasses=.True.

    CalculateLowEnergy=.True.

    SolutionTadpoleNr = 1 

    SUSYrunningFromMZ=.False.

    Write_WHIZARD =.False.     

    !Write_HiggsBounds = .True.   

    WidthToBeInvisible = 0

    WriteTreeLevelTadpoleSolutions=.False.

    WriteGUTvalues=.False.

    WriteEffHiggsCouplingRatios=.True.

    WriteHiggsDiphotonLoopContributions=.False.

    WriteTreeLevelTadpoleParameters=.True.

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

End Subroutine Initialize
