module UnsteadyAero
   
   use NWTC_Library   
   use UnsteadyAero_Types
   use AirfoilInfo
   
   implicit none 

private


   public :: UA_Init
   public :: UA_UpdateDiscState
   public :: UA_UpdateStates
   public :: UA_CalcOutput


contains

!==============================================================================   
subroutine GetSteadyOutputs(AFInfo, AOA, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
! Called by : UA_CalcOutput
! Calls  to : CubicSplineInterpM   
!..............................................................................
   type(AFInfoType), intent(in   ) :: AFInfo                ! Airfoil info structure 
   real(ReKi),       intent(in   ) :: AOA                   ! Angle of attack (rad)
   real(ReKi),       intent(  out) :: Cl                    ! Coefficient of lift (-)
   real(ReKi),       intent(  out) :: Cd                    ! Coefficient of drag (-)
   real(ReKi),       intent(  out) :: Cm                    ! Pitch moment coefficient (-)
   real(ReKi),       intent(  out) :: Cd0                   ! Minimum Cd value (-)
   integer(IntKi),   intent(  out) :: ErrStat               ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg                ! Error message if ErrStat /= ErrID_None
   
   real                            :: IntAFCoefs(4)         ! The interpolated airfoil coefficients.
   integer                         :: s1                    ! Number of columns in the AFInfo structure
   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
      ! TODO: Extend this to use the UnsteadyAero module to determine the Cl, Cd, Cm info, as needed.  We may need to be tracking whether this call is
      ! part of an UpdateStates action or a CalcOutput action, because the calculation chain differs for the two.
      
      
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   s1 = size(AFInfo%Table(1)%Coefs,2)
   if (s1 < 3) then
      ErrMsg  = 'The Airfoil info table must contains columns for lift, drag, and pitching moment'
      ErrStat = ErrID_Fatal
      return
   end if
   Cd0 =   AFInfo%Table(1)%UA_BL%Cd0
   IntAFCoefs(1:s1) = CubicSplineInterpM( 1.0*real( AOA ) &
                                             , AFInfo%Table(1)%Alpha &
                                             , AFInfo%Table(1)%Coefs &
                                             , AFInfo%Table(1)%SplineCoefs &
                                             , ErrStat, ErrMsg )
   if (ErrStat > ErrID_None) return
      
   Cl = IntAFCoefs(1)
   Cd = IntAFCoefs(2)
   Cm = IntAFCoefs(3)
   
end subroutine GetSteadyOutputs
!==============================================================================

!==============================================================================
real(ReKi) function Get_ds( U, c, dt )
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   real(ReKi), intent(in   ) :: U       ! air velocity magnitude relative to the airfoil (m/s)
   real(ReKi), intent(in   ) :: c       ! cord length (m)
   real(DbKi), intent(in   ) :: dt      ! time step length (s)
   
   Get_ds = 2.0_ReKi*U*dt/c
   
end function Get_ds
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Kupper( dt, x, x_minus1 )
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   real(DbKi), intent(in   ) :: dt          ! time step length (s)
   real(ReKi), intent(in   ) :: x           ! quantity at current timestep
   real(ReKi), intent(in   ) :: x_minus1    ! quantity at previous timestep
   
      ! Implements eqn 1.18b or 1.19b
   Get_Kupper = (x - x_minus1) / dt
   
end function Get_Kupper
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Beta_M( M ) 
! Compute the Prandlt-Glauert compressibility correction factor
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   real(ReKi), intent(in   ) :: M   ! Mach number (-)
   
   Get_Beta_M = sqrt(1 - M**2)
   
end function Get_Beta_M
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Beta_M_Sqrd( M ) 
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: M   ! Mach number (-)
   
   Get_Beta_M_Sqrd = 1 - M**2
   
end function Get_Beta_M_Sqrd
!==============================================================================   

!==============================================================================
real(ReKi) function Get_k_( M, beta_M, C_nalpha          )
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: M               ! Mach number
   real(ReKi), intent(in   ) :: beta_M          ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2)
   real(ReKi), intent(in   ) :: C_nalpha        ! slope of the 2D pitching moment curve 
   
      ! Implements equation 1.11a/b
   Get_k_   = 1.0_ReKi / ( (1.0_ReKi - M) + C_nalpha * M**2 * beta_M *  0.413_ReKi  )
   
end function Get_k_
!==============================================================================   

!==============================================================================
real(ReKi) function Get_ExpEqn( dt, T, Y_minus1, x, x_minus1 )
! Compute the deficiency function
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: dt                              ! numerator of the exponent
   real(ReKi), intent(in   ) :: T                               ! denominator of the exponent
   real(ReKi), intent(in   ) :: Y_minus1                        ! previous value of the computed decay function
   real(ReKi), intent(in   ) :: x                               ! current value of x
   real(ReKi), intent(in   ) :: x_minus1                        ! previous value of x
   
   Get_ExpEqn = Y_minus1*exp(-dt/T) + (x - x_minus1)*exp(-0.5_ReKi*dt/T)

end function Get_ExpEqn
!==============================================================================   

!==============================================================================  
real(ReKi) function Get_Pitchrate( dt, alpha_cur, alpha_minus1, U, c )
! This is the non-dimensional pitching rate, given as q in the documentation: equation 1.8  
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(DbKi), intent(in   ) :: dt                              ! time step length (s)
   real(ReKi), intent(in   ) :: alpha_cur                       ! angle of attack, current timestep  (radians)
   real(ReKi), intent(in   ) :: alpha_minus1                    ! angle of attack, previous timestep  (radians)
   real(ReKi), intent(in   ) :: U                               ! air velocity magnitude relative to the airfoil (m/s)
   real(ReKi), intent(in   ) :: c                               ! chord length (m)

!IF (EqualRealNos(U,0.0_ReKi)) then
!   call ProgWarn('GetPitchRate: division by zero.')
!   Get_Pitchrate = HUGE(Get_Pitchrate)
!   RETURN
!END IF

      ! Implements equation 1.8
   Get_Pitchrate = ( alpha_cur - alpha_minus1 ) * c / (U*dt)  
   
end function Get_Pitchrate
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cn_nc( M, T, K, Kprime )
!  Compute either Cn_alpha_nc or Cn_q_nc
!  Eqn 1.18a or 1.19a
!  Called by : ComputeKelvinChain
!  Calls  to : NONE
!..............................................................................

   real(ReKi)                :: M              ! mach number
   real(ReKi), intent(in   ) :: T              ! mach-dependent time constant related to alpha and non-circulatory terms (???)
   real(ReKi)                :: K              ! Either K_alpha or K_q, backward finite difference of either alpha or q
   real(ReKi)                :: Kprime         ! Either Kprime_alpha or Kprime_q, deficiency functions
   
      ! Implements equation 1.18a or 1.19a  
   Get_Cn_nc =  T * ( K - Kprime ) / M
   
end function Get_Cn_nc 
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cm_q_circ( C_nalpha, beta_M, q_cur, K3prime_q, c, U )
! Implements eqn 1.21
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: C_nalpha                       ! slope of the 2D pitching moment curve     
   real(ReKi), intent(in   ) :: beta_M                         ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2) 
   real(ReKi), intent(in   ) :: q_cur                          ! non-dimensional pitching rate, current timestep   
   real(ReKi), intent(in   ) :: K3prime_q                      ! deficiency function
   real(ReKi), intent(in   ) :: c                              ! chord length (m)   
   real(ReKi), intent(in   ) :: U                              ! air velocity magnitude relative to the airfoil (m/s)

   Get_Cm_q_circ = -C_nalpha*(q_cur -K3prime_q)*c/(16.0*beta_M*U)  ! TODO: Check units
   
end function Get_Cm_q_circ
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cm_q_nc( M, k_mq, T_I, Kq, Kprimeprime_q )
! Called by : ComputeKelvinChain
! Calls  to : NONE   
!..............................................................................

   real(ReKi), intent(in   ) :: M                     ! Mach number
   real(ReKi), intent(in   ) :: k_mq                  ! airfoil parameter                     
   real(ReKi), intent(in   ) :: T_I                   ! time constant        
   real(ReKi), intent(in   ) :: Kq                    ! backward finite difference of q    
   real(ReKi), intent(in   ) :: Kprimeprime_q         ! deficiency function             
    
   ! Implements eqn 1.25a
   Get_Cm_q_nc = -7*(k_mq**2)*T_I*(Kq-Kprimeprime_q)/(12.0*M)
   
end function Get_Cm_q_nc
!==============================================================================   

!==============================================================================
real(ReKi) function Get_alpha_e( alpha, alpha0, X1, X2 )
! Called by : ComputeKelvinChain
! Calls  to : NONE   
!..............................................................................
   
   real(ReKi), intent(in   ) :: alpha          ! angle of attack (radians)
   real(ReKi), intent(in   ) :: alpha0         ! zero lift angle of attack (radians)
   real(ReKi), intent(in   ) :: X1             ! deficiency function 
   real(ReKi), intent(in   ) :: X2             ! deficiency function 
   
   ! Implements equation 1.14
   Get_Alpha_e   = (alpha - alpha0) - X1 - X2
   
end function Get_Alpha_e
!==============================================================================   

   
!==============================================================================
! TE Flow Separation Equations                                                !
!==============================================================================

!==============================================================================
real(ReKi) function Get_f( alpha, alpha0, alpha1, alpha2, S1, S2, S3, S4 )
! Compute either fprime or fprimeprime
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: alpha         ! angle of attack (radians)
   real(ReKi), intent(in   ) :: alpha0        ! zero lift angle of attack (radians)
   real(ReKi), intent(in   ) :: alpha1        ! angle of attack at f = 0.7, approximately the stall angle; for alpha >= alpha0 (radians)
   real(ReKi), intent(in   ) :: alpha2        ! angle of attack at f = 0.7, approximately the stall angle; for alpha < alpha0 (radians)
   real(ReKi), intent(in   ) :: S1            ! constant in the f-curve best-fit, alpha >= alpha0  
   real(ReKi), intent(in   ) :: S2            ! constant in the f-curve best-fit, alpha >= alpha0  
   real(ReKi), intent(in   ) :: S3            ! constant in the f-curve best-fit, alpha < alpha0    
   real(ReKi), intent(in   ) :: S4            ! constant in the f-curve best-fit, alpha < alpha0   

      
      ! Implements Equation 1.30
   if (alpha > alpha1) then
      Get_f = 0.04_ReKi + 0.66_ReKi*exp((alpha1-alpha)/S2)
   else if (alpha < alpha2) then
      Get_f = 0.04_ReKi + 0.66_ReKi*exp((alpha-alpha2)/S4)
   else if (alpha>= alpha0) then
      Get_f = 1 - 0.3_ReKi*exp((alpha-alpha1)/S1)
   else ! alpha < alpha0
      Get_f = 1 - 0.3_ReKi*exp((alpha2-alpha)/S3)
   end if
   
end function Get_f
!==============================================================================   

!==============================================================================   
real(ReKi) function Get_Cn_FS( Cn_alpha_q_nc, Cn_alpha_q_circ, fprimeprime, DSMod )
!  
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

    
   real(ReKi), intent(in   ) :: Cn_alpha_q_nc       ! non-circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ     ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime         ! lagged version of fprime accounting for unsteady boundary layer response
   integer   , intent(in   ) :: DSMod               ! model for the dynamic stall equations.  [1 = Leishman/Beddoes, 2 = Gonzalez, 3 = Minnema]
   
   if ( DSMod == 1 ) then
      ! use equation 1.35
      Get_Cn_FS   = Cn_alpha_q_nc + Cn_alpha_q_circ * ( (1.0_ReKi + sqrt(fprimeprime) ) / 2.0_ReKi )**2
   else if ( DSMod == 2 ) then
      ! use equation 1.36
      Get_Cn_FS   = Cn_alpha_q_nc + Cn_alpha_q_circ *  ( (1.0_ReKi + 2.0_ReKi*sqrt(fprimeprime) ) / 3.0_ReKi  )**2
   else
      ! implementation error!  Cannot have DSMod other than 0,1,2
   end if


end function Get_Cn_FS
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cc_FS( eta_e, alpha_e, alpha0, Cn_alpha_q_circ, fprimeprime, DSMod )
!  
! Called by : NONE
! Calls  to : NONE
!..............................................................................
   
   real(ReKi), intent(in   ) :: eta_e
   real(ReKi), intent(in   ) :: alpha_e                    ! effective angle of attack at 3/4 chord  TODO: verify 3/4 and not 1/4
   real(ReKi), intent(in   ) :: alpha0                     ! zero lift angle of attack (radians)
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ            ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime                ! lagged version of fprime accounting for unsteady boundary layer response
   integer   , intent(in   ) :: DSMod                      ! model for the dynamic stall equations.  [1 = Leishman/Beddoes, 2 = Gonzalez, 3 = Minnema]
   
   if ( DSMod == 1 ) then
      ! use equation 1.37
      Get_Cc_FS   = Cn_alpha_q_circ * sqrt(fprimeprime) * tan(alpha_e + alpha0) * eta_e
   else if ( DSMod == 2 ) then
      ! use equation 1.38
      Get_Cc_FS   = Cn_alpha_q_circ * (sqrt(fprimeprime) - 0.2_ReKi) * tan(alpha_e + alpha0) * eta_e
   else
      ! implementation error!  Cannot have DSMod other than 0,1,2
   end if

end function Get_Cc_FS
!==============================================================================   


!==============================================================================
! Dynamic Stall Equations                                                     !
!==============================================================================

!==============================================================================
real(ReKi) function Get_C_V( Cn_alpha_q_circ, fprimeprime, DSMod )
! 
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: Cn_alpha_q_circ               ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime                   ! lagged version of fprime accounting for unsteady boundary layer response
   integer   , intent(in   ) :: DSMod                         ! model for the dynamic stall equations.  [1 = Leishman/Beddoes, 2 = Gonzalez, 3 = Minnema]
   
    if ( DSMod == 1 ) then
      ! use equation 1.47
      Get_C_V   = Cn_alpha_q_circ *  ( 1.0_ReKi - ( 0.5_ReKi + 0.5_ReKi*sqrt(fprimeprime) ) )**2
   else if ( DSMod == 2 ) then
      ! use equation 1.48
      Get_C_V   = Cn_alpha_q_circ * ( 1.0_ReKi - (1.0_ReKi + 2.0_ReKi*sqrt(fprimeprime) )/3.0_ReKi  )**2
   else
      ! implementation error!  Cannot have DSMod other than 0,1,2
   end if
   

end function Get_C_V
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cn_v ( ds, T_V, Cn_v_minus1, C_V, C_V_minus1, tau_V, T_VL, Kalpha, alpha, alpha0 ) ! do I pass VRTX flag or tau_V?
! Implements Equation 1.45 or 1.49 
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   
   real(ReKi), intent(in   ) :: ds                     ! non-dimensionalized distance parameter
   real(ReKi), intent(in   ) :: T_V                    ! backwards finite difference of the non-dimensionalized distance parameter
   real(ReKi), intent(in   ) :: Cn_v_minus1            ! normal force coefficient due to the presence of LE vortex, previous time step (-)
   real(ReKi), intent(in   ) :: C_V                    ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi), intent(in   ) :: C_V_minus1             ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex, previous time step (-)
   real(ReKi), intent(in   ) :: tau_V                  ! time variable, tracking the travel of the LE vortex over the airfoil suction surface (-)
   real(ReKi), intent(in   ) :: T_VL                   ! time variable associated with the vortex advection process; it represents the non-dimensional time in semi-chords needed for a vortex to travel from leading edge to trailing edge
   real(ReKi), intent(in   ) :: Kalpha                 ! backwards finite difference of alpha
   real(ReKi), intent(in   ) :: alpha                  ! angle of attack (rad)
   real(ReKi), intent(in   ) :: alpha0                 ! zero lift angle of attack (rad)
  
   real(ReKi)                :: factor
   
   factor = (alpha - alpha0) * Kalpha
   
   if (tau_V > T_VL .AND. factor >= 0.0_ReKi ) then  
      Get_Cn_v = Cn_v_minus1*exp(-ds/T_V)   ! Eqn 1.49
   else      
      Get_Cn_v = Get_ExpEqn( ds, T_V, Cn_v_minus1, C_V, C_V_minus1 )   ! Eqn 1.45
   end if
   
end function Get_Cn_v
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cm( Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, Cn_alpha_q_circ, fprimeprime, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, tau_v )
! Leishman's pitching moment coefficient about the 1/4 chord
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   real(ReKi), intent(in   ) :: Cm0                           ! 2D pitching moment coefficient at zero lift, positive if nose is up
   real(ReKi), intent(in   ) :: k0                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi), intent(in   ) :: k1                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi), intent(in   ) :: k2                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi), intent(in   ) :: k3                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi), intent(in   ) :: T_VL                          ! time variable associated with the vortex advection process; it represents the non-dimensional time in semi-chords needed for a vortex to travel from leading edge to trailing edge
   real(ReKi), intent(in   ) :: x_cp_bar                      ! airfoil parameter for calulating x_cp_v
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ               ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime                   ! lagged version of fprime accounting for unsteady boundary layer response
   real(ReKi), intent(in   ) :: Cm_q_circ                     ! circulatory component of the pitching moment coefficient response to a step change in q
   real(ReKi), intent(in   ) :: Cn_alpha_nc                   ! non-circulatory component of the normal force coefficient response to step change in alpha
   real(ReKi), intent(in   ) :: Cm_q_nc                       ! non-circulatory component of the moment coefficient response to step change in q
   real(ReKi), intent(in   ) :: Cn_v                          ! normal force coefficient due to the presence of LE vortex
   real(ReKi), intent(in   ) :: tau_v                         ! time variable, tracking the travel of the LE vortex over the airfoil suction surface (-)
   
      ! local variables
   real(ReKi)                :: x_cp_hat                      ! center-of-pressure distance from LE in chord fraction
   real(ReKi)                :: Cm_v                          ! pitching moment coefficient due to the presence of LE vortex
   
      ! Eqn 1.40
   x_cp_hat = k0 + k1*(1-fprimeprime) + k2*sin(pi*fprimeprime**k3)
   
      ! Eqn 1.154
   Cm_v     = -x_cp_bar*( 1-cos( pi*tau_v/T_VL ) )*Cn_v
   
      ! Eqn 1.55 = Cm0 - 1.13 + 1.21 - 1.18a/4 + 1.25a + 1.54  
   Get_Cm   = Cm0 - Cn_alpha_q_circ*(x_cp_hat - 0.25_ReKi) + Cm_q_circ - Cn_alpha_nc / 4.0_ReKi + Cm_q_nc + Cm_v
   
   
end function Get_Cm
!==============================================================================   


!==============================================================================                              
subroutine ComputeKelvinChain( i, j, u, p, xd, OtherState, AFInfo, Cn_prime, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                               St_sh, Kalpha, alpha_e, alpha0, dalpha0, eta_e, Kq, q_cur, X1, X2, &
                               Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, fprime, &
                               Cn_alpha_q_circ, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, errStat, errMsg )
! 
! Called by : DRIVER
! Calls  to : Get_Beta_M_Sqrd, Get_Beta_M, AFI_GetAirfoilParams, Get_ds, Get_Pitchrate, Get_k_, Get_Kupper, Get_ExpEqn
!             Get_Cn_nc, Get_alpha_e, Get_Cm_q_circ, Get_Cm_q_nc, Get_f, Get_Cn_FS, Get_C_V, Get_Cn_v
!...............................................................................

      
   integer(IntKi),                         intent(in   ) :: i,j               ! Blade node indices
   type(UA_InputType),                     intent(in   ) :: u                 ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   type(UA_ParameterType),                 intent(in   ) :: p                 ! Parameters   TODO: this should only be in, but is needed because of a copy of AFInfo in a called subroutine. GJH 1/5/2015
   type(UA_DiscreteStateType),             intent(in   ) :: xd                ! Input: Discrete states at t;
   type(UA_OtherStateType),                intent(in   ) :: OtherState        ! Other/optimization states
   type(AFInfoType),                       intent(in   ) :: AFInfo            ! The airfoil parameter data
   real(ReKi),                             intent(  out) :: Cn_prime          !
   real(ReKi),                             intent(  out) :: Cn1               ! critical value of Cn_prime at LE separation for alpha >= alpha0
   real(ReKi),                             intent(  out) :: Cn2               ! critical value of Cn_prime at LE separation for alpha < alpha0
   real(ReKi),                             intent(  out) :: Cd0               !
   real(ReKi),                             intent(  out) :: Cm0               ! 2D pitching moment coefficient at zero lift, positive if nose is up
   real(ReKi),                             intent(  out) :: k0                ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),                             intent(  out) :: k1                ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),                             intent(  out) :: k2                ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),                             intent(  out) :: k3                ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),                             intent(  out) :: T_VL              ! time variable associated with the vortex advection process; it represents the non-dimensional time in semi-chords needed for a vortex to travel from leading edge to trailing edge
   real(ReKi),                             intent(  out) :: x_cp_bar          ! airfoil parameter for calulating x_cp_v
   real(ReKi),                             intent(  out) :: St_sh             !
   real(ReKi),                             intent(  out) :: Kalpha            ! backwards finite difference of alpha
   real(ReKi),                             intent(  out) :: alpha_e           ! effective angle of attack at 3/4 chord  TODO: verify 3/4 and not 1/4                          
   real(ReKi),                             intent(  out) :: alpha0            ! zero lift angle of attack (radians)
   real(ReKi),                             intent(  out) :: dalpha0           !
   real(ReKi),                             intent(  out) :: eta_e             !
   real(ReKi),                             intent(  out) :: Kq                !
   real(ReKi),                             intent(  out) :: q_cur             !
   real(ReKi),                             intent(  out) :: X1                !
   real(ReKi),                             intent(  out) :: X2                !
   real(ReKi),                             intent(  out) :: Kprime_alpha      !
   real(ReKi),                             intent(  out) :: Kprime_q          !
   real(ReKi),                             intent(  out) :: Dp                !
   real(ReKi),                             intent(  out) :: Cn_pot            !
   real(ReKi),                             intent(  out) :: Cc_pot            !
   real(ReKi),                             intent(  out) :: Cn_alpha_q_circ   !
   real(ReKi),                             intent(  out) :: Cm_q_circ         !
   real(ReKi),                             intent(  out) :: Cn_alpha_nc       ! non-circulatory component of the normal force coefficient response to step change in alpha
   real(ReKi),                             intent(  out) :: Cm_q_nc           ! non-circulatory component of the moment coefficient response to step change in q
   real(ReKi),                             intent(  out) :: fprimeprime       !
   real(ReKi),                             intent(  out) :: Df                !
   real(ReKi),                             intent(  out) :: fprime            !
   real(ReKi),                             intent(  out) :: Cn_v              ! normal force coefficient due to the presence of LE vortex
   real(ReKi),                             intent(  out) :: C_V               ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi),                             intent(  out) :: Cn_FS             !
                                                                              !
   integer(IntKi),                         intent(  out) :: errStat           ! Error status of the operation
   character(*),                           intent(  out) :: errMsg            ! Error message if ErrStat /= ErrID_None


   real(ReKi)                :: ds                                            ! non-dimensionalized distance parameter
   real(ReKi)                :: M                                             ! Mach number (-)
   real(ReKi)                :: beta_M                                        ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2)
   real(ReKi)                :: beta_M_Sqrd                                   ! square of the Prandtl-Glauert compressibility correction factor,  (1-M**2)
   real(ReKi)                :: Cn_alpha_q_nc                                 ! non-circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi)                :: k_alpha                                       !
   real(ReKi)                :: T_alpha                                       !
   real(ReKi)                :: T_q                                           !
   real(ReKi)                :: T_f                                           !
   real(ReKi)                :: T_V                                           ! backwards finite difference of the non-dimensionalized distance parameter
   real(ReKi)                :: T_I                                           !
   real(ReKi)                :: C_nalpha_circ                                 ! slope of the circulatory normal force coefficient vs alpha curve
   real(ReKi)                :: T_f0                                          ! initial value of T_f, airfoil specific, used to compute D_f and fprimeprime
   real(ReKi)                :: T_V0                                          ! initial value of T_V, airfoil specific, time parameter associated with the vortex lift decay process, used in Cn_v
   real(ReKi)                :: T_p                                           ! boundary-layer, leading edge pressure gradient time parameter; used in D_p; airfoil specific
   real(ReKi)                :: b1                                            ! airfoil constant derived from experimental results, usually 0.14
   real(ReKi)                :: b2                                            ! airfoil constant derived from experimental results, usually 0.53
   real(ReKi)                :: b5                                            ! airfoil constant derived from experimental results, usually 5.0
   real(ReKi)                :: A1                                            ! airfoil constant derived from experimental results, usually 0.3
   real(ReKi)                :: A2                                            ! airfoil constant derived from experimental results, usually 0.7
   real(ReKi)                :: A5                                            ! airfoil constant derived from experimental results, usually 1.0
   real(ReKi)                :: C_nalpha                                      !
   real(ReKi)                :: alpha1                                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha >= alpha0 (radians)
   real(ReKi)                :: alpha2                                        ! angle of attack at f = 0.7, approximately the stall angle; for alpha < alpha0 (radians)
   real(ReKi)                :: S1                                            ! constant in the f-curve best-fit, alpha >= alpha0
   real(ReKi)                :: S2                                            ! constant in the f-curve best-fit, alpha >= alpha0
   real(ReKi)                :: S3                                            ! constant in the f-curve best-fit, alpha <  alpha0
   real(ReKi)                :: S4                                            ! constant in the f-curve best-fit, alpha <  alpha0
   real(ReKi)                :: Cn_q_nc                                       !
   real(ReKi)                :: alpha_f                                       !
   real(ReKi)                :: k_q                                           !
   real(ReKi)                :: Kalpha_minus1                                 !
   real(ReKi)                :: Kq_minus1                                     !
   real(ReKi)                :: alpha_minus1                                  !
   real(ReKi)                :: alpha_minus2                                  !
   real(ReKi)                :: q_minus1                                      !
   real(ReKi)                :: q_minus2                                      !
   real(ReKi)                :: fprime_minus1                                 !
   real(ReKi)                :: Cn_pot_minus1                                 !
   real(ReKi)                :: k1_hat                                        !
   real(ReKi)                :: K3prime_q                                     !
   real(ReKi)                :: k_mq                                          !
   real(ReKi)                :: Kprimeprime_q                                 !
   
   

   M           = u%U / p%a_s
   beta_M_Sqrd = Get_Beta_M_Sqrd(M)
   beta_M      = Get_Beta_M(M)
   T_I         = p%c(i,j) / p%a_s                                                  ! Eqn 1.11c
   ds          = Get_ds       ( u%U, p%c(i,j), p%dt )                              ! Eqn 1.5b
   
   ! Lookup values using Airfoil Info module
   call AFI_GetAirfoilParams( AFInfo, M, u%Re, u%alpha, alpha0, alpha1, alpha2, eta_e, C_nalpha, C_nalpha_circ, &
                              T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, b5, A1, A2, A5, S1, S2, S3, S4, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, k1_hat, x_cp_bar, errMsg, errStat )           
   
   if (OtherState%FirstPass) then
      alpha_minus1 = u%alpha
      alpha_minus2 = u%alpha
   else
      alpha_minus1 = xd%alpha_minus1(i,j)
      alpha_minus2 = xd%alpha_minus2(i,j)
   end if
   
   dalpha0  = u%alpha - alpha0
   q_cur    = Get_Pitchrate( p%dt, u%alpha, alpha_minus1, u%U, p%c(i,j)  )  ! Eqn 1.8
   
   if (OtherState%FirstPass) then
      q_minus1 = q_cur   
      q_minus2 = q_cur 
   else
      q_minus1 = xd%q_minus1(i,j)   
      q_minus2 = xd%q_minus2(i,j) 
   end if
   
   k_alpha  = Get_k_       ( M, beta_M, C_nalpha/2.0_ReKi )                 ! Eqn 1.11a
   k_q      = Get_k_       ( M, beta_M, C_nalpha          )                 ! Eqn 1.11b
   T_alpha  = T_I * k_alpha                                                 ! Eqn 1.10a
   T_q      = T_I * k_q                                                     ! Eqn 1.10b
      
      ! These quantities are needed for the update state calculations, but are then updated themselves based on the logic which follows
   T_f           = T_f0 / OtherState%sigma1(i,j)       ! Eqn 1.34
   T_V           = T_V0 / OtherState%sigma3(i,j)       ! Eqn 1.46
   
      ! Compute Kalpha using Eqn 1.18b and Kalpha_minus1 using Eqn 1.18b but with time-shifted alphas
   Kalpha          = Get_Kupper( p%dt, u%alpha     , alpha_minus1 )
   Kalpha_minus1   = Get_Kupper( p%dt, alpha_minus1, alpha_minus2 ) 
   
      ! Compute Kq using Eqn 1.19b and Kq_minus1 using Eqn 1.19b but with time-shifted alphas
   Kq              = Get_Kupper( p%dt, q_cur   , q_minus1     )
   Kq_minus1       = Get_Kupper( p%dt, q_minus1, q_minus2     )
   
      ! Compute Kprime_alpha using Eqn 1.18c
   Kprime_alpha = Get_ExpEqn( real(p%dt), T_alpha, xd%Kprime_alpha_minus1(i,j), Kalpha, Kalpha_minus1 )
   
      ! Compute Kprime_q using Eqn 1.19c    
   Kprime_q     = Get_ExpEqn( real(p%dt), T_q    , xd%Kprime_q_minus1(i,j)    ,  Kq   , Kq_minus1     )
   
      ! Compute Cn_alpha_nc using Eqn 1.18a
      ! Depends on T_alpha, M, alpha (1.18b), Kprime_alpha (1.18c), xd%alpha_minus1, xd%Kprime_alpha_minus1
   Cn_alpha_nc     = Get_Cn_nc( M, 4.0_ReKi*T_alpha, Kalpha, Kprime_alpha )
   
      ! Compute Cn_alpha_nc using Eqn 1.19a
   Cn_q_nc         = Get_Cn_nc( M,          T_q    , Kq    , Kprime_q     )
   
      ! Compute Cn_alpha_q_nc using Eqn 1.17
      ! Depends on Cn_alpha_nc (1.18a), Cn_q_nc (1.19a)
   Cn_alpha_q_nc   = Cn_alpha_nc + Cn_q_nc  
   
      ! Compute X1 using Eqn 1.15a   
   X1           = Get_ExpEqn( ds*beta_M_Sqrd*b1, 1.0_ReKi, xd%X1_minus1(i,j), A1*(u%alpha - alpha_minus1), 0.0_ReKi )
   
      ! Compute X2 using Eqn 1.15b
   X2           = Get_ExpEqn( ds*beta_M_Sqrd*b2, 1.0_ReKi, xd%X2_minus1(i,j), A2*(u%alpha - alpha_minus1), 0.0_ReKi )
   
      ! Compute alpha_e using Eqn 1.14
   alpha_e         = Get_alpha_e( u%alpha, alpha0, X1, X2 )
   
      ! Compute Cn_alpha_q_circ using Eqn 1.13
   Cn_alpha_q_circ = C_nalpha_circ * alpha_e
   
      ! Compute K3prime_q using Eqn 1.22
   K3prime_q       = Get_ExpEqn( b5*beta_M_Sqrd*ds, 1.0_ReKi, xd%K3prime_q_minus1(i,j),  A5*Kq*real(p%dt), A5*Kq_minus1*real(p%dt) )
   
      ! Compute Cm_q_circ using Eqn 1.21
   Cm_q_circ       = Get_Cm_q_circ( C_nalpha, beta_M, q_cur, K3prime_q, p%c(i,j), u%U )
   
      ! Compute Cn_pot using eqn 1.20
      ! Depends on Cn_alpha_q_circ (1.13) , Cn_alpha_q_nc (1.17)
   Cn_pot          = Cn_alpha_q_circ + Cn_alpha_q_nc
   
      ! Eqn 1.25b
   k_mq            = 7.0/(15.0*(1.0-M)+1.5*C_nalpha*A5*b5*beta_M*M**2)
   
      ! Eqn 1.25d
   Kprimeprime_q   = Get_ExpEqn( real(p%dt), k_mq**2*T_I   , xd%Kprimeprime_q_minus1(i,j)    ,  Kq   , Kq_minus1     )
   
      ! Compute Cm_q_nc using eqn 1.25a
   Cm_q_nc         = Get_Cm_q_nc( M, k_mq, T_I, Kq, Kprimeprime_q )
   
   
      ! Compute Cc_pot using eqn 1.28
   Cc_pot          = Cn_pot*tan(alpha_e+alpha0)
   
   
   if (OtherState%FirstPass) then
      Cn_pot_minus1 = Cn_pot
   else
      Cn_pot_minus1 = xd%Cn_pot_minus1(i,j)
   end if
   
      ! Compute Dp using Eqn 1.32b
   Dp            = Get_ExpEqn( ds, T_p, xd%Dp_minus1(i,j), Cn_pot, Cn_pot_minus1 )
   
      ! Compute Cn_prime using Eqn 1.32a
      ! Depends on Cn_pot (1.20) and Dp (1.32b)
   Cn_prime      = Cn_Pot - Dp
   
      ! Compute alpha_f using Eqn 1.31
   alpha_f       = Cn_prime / C_nalpha + alpha0
   
      ! Compute fprime using Eqn 1.30 and Eqn 1.31
   fprime        = Get_f( alpha_f, alpha0, alpha1, alpha2, S1, S2, S3, S4)
   if (OtherState%FirstPass) then
      fprime_minus1 = fprime
   else
      fprime_minus1 = xd%fprime_minus1(i,j)
   end if
   
      ! Compute Df using Eqn 1.33b   
   Df            = Get_ExpEqn( ds, T_f, xd%Df_minus1(i,j), fprime, fprime_minus1 )
   
      ! Compute fprimeprime using Eqn 1.33a
   fprimeprime   = fprime - Df
   
      ! Compute Cn_FS using Eqn 1.35 or 1.36 depending on option selected
   Cn_FS         = Get_Cn_FS( Cn_alpha_q_nc, Cn_alpha_q_circ, fprimeprime, p%DSMod )
   
      ! Compute Cc_FS usign Eqn 1.37 or 1.38, but this is not used in the final value of Cc and does not contribute to any states!  
      ! So we will leave it commented out. GJH 2/23/2015
   !Cc_FS         = Get_Cc_FS( eta_e, alpha_e, alpha0, Cn_alpha_q_circ, fprimeprime, p%DSMod )
   
      ! Compute C_V using Eqn 1.47 or 1.48 depending on option selected
   C_V           = Get_C_V  ( Cn_alpha_q_circ, fprimeprime, p%DSMod )
   
      ! Compute Cn_v using either Eqn 1.45 or 1.49 depending on operating conditions
   Cn_v          = Get_Cn_v ( ds, T_V, xd%Cn_v_minus1(i,j), C_V, xd%C_V_minus1(i,j), xd%tau_V(i,j), T_VL, Kalpha, u%alpha, alpha0 ) ! do I pass VRTX flag or tau_V?
   
      ! Finally, compute Cn using Eqn 1.50
   !Cn            = Cn_FS + Cn_v

 end subroutine ComputeKelvinChain
!==============================================================================   
                          

!==============================================================================
! Framework Routines                                                          !
!==============================================================================                               
      

!==============================================================================
subroutine UA_SetParameters( dt, InitInp, p, errStat, errMsg )
! 
! Called by : UA_Init
! Calls  to : NONE
!..............................................................................
   
   real(DbKi),                             intent(inout)  :: dt          ! time step length (s)
   type(UA_InitInputType),       intent(inout)  :: InitInp     ! input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UA_ParameterType),       intent(inout)  :: p           ! parameters
   integer(IntKi),                         intent(  out)  :: errStat     ! error status of the operation
   character(*),                           intent(  out)  :: errMsg      ! error message if ErrStat /= ErrID_None

      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""
   p%dt         = dt
   
   allocate(p%c(InitInp%nNodesPerBlade,InitInp%numBlades), stat = errStat)
   if (errStat > ErrID_None) then
      ! Set errmessage and return
      return
   end if

   p%c          = InitInp%c        
   p%numBlades  = InitInp%numBlades
   p%nNodesPerBlade  = InitInp%nNodesPerBlade
   p%DSMod      = InitInp%DSMod    
   p%a_s        = InitInp%a_s       
   
   
   
end subroutine UA_SetParameters
!==============================================================================   

!==============================================================================
subroutine UA_InitStates( p, xd, OtherState, errStat, errMsg )  
! Called by : UA_Init
! Calls  to : NONE
!..............................................................................

   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Initial other/optimization states
   integer(IntKi),               intent(  out)  :: errStat     ! Error status of the operation
   character(*),                 intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None
                                                                         
   
   errMsg  = ""
   errStat = ErrID_None
   
      ! allocate all the state arrays
   allocate(xd%alpha_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%alpha_minus2(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%q_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%q_minus2(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%X1_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%X2_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Kprime_alpha_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Kprime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Kprimeprime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%K3prime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Dp_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Cn_pot_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%fprimeprime_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Df_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%fprime_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%tau_V(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Cn_v_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%C_V_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(OtherState%sigma1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(OtherState%sigma3(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(OtherState%TESF(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(OtherState%LESF(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(OtherState%VRTX(p%nNodesPerBlade,p%numBlades), stat = errStat)
   
   OtherState%sigma1    = 1.0_ReKi
   OtherState%sigma3    = 1.0_ReKi
   OtherState%TESF      = .FALSE.  
   OtherState%LESF      = .FALSE.   
   OtherState%VRTX      = .FALSE.  
   OtherState%FirstPass = .true.
   
   xd%alpha_minus1         = 0.0_ReKi
   xd%alpha_minus2         = 0.0_ReKi
   xd%q_minus1             = 0.0_ReKi
   xd%q_minus2             = 0.0_ReKi
   xd%X1_minus1            = 0.0_ReKi
   xd%X2_minus1            = 0.0_ReKi
   xd%Kprime_alpha_minus1  = 0.0_ReKi
   xd%Kprime_q_minus1      = 0.0_ReKi
   xd%Dp_minus1            = 0.0_ReKi
   xd%Cn_pot_minus1        = 0.0_ReKi
   xd%K3prime_q_minus1     = 0.0_ReKi
   xd%Kprimeprime_q_minus1 = 0.0_ReKi  
   xd%fprimeprime_minus1   = 0.0_ReKi
   xd%Df_minus1            = 0.0_ReKi
   xd%fprime_minus1        = 0.0_ReKi
   xd%tau_V                = 0.0_ReKi 
   xd%Cn_v_minus1          = 0.0_ReKi
   xd%C_V_minus1           = 0.0_ReKi

end subroutine UA_InitStates 
!==============================================================================   

!==============================================================================
subroutine UA_Init( InitInp, u, p, xd, OtherState, y, Interval, &
                              InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
! Called by : DRIVER
! Calls  to : NWTC_Init, UA_SetParameters, UA_InitStates
!..............................................................................

   type(UA_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UA_InputType),           intent(in   )  :: u           ! An initial guess for the input; input mesh must be defined
   type(UA_ParameterType),       intent(  out)  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(  out)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(  out)  :: OtherState  ! Initial other/optimization states
   type(UA_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                         !   only the output mesh is initialized)
   real(DbKi),                   intent(inout)  :: interval    ! Coupling interval in seconds: the rate that
                                                                         !   (1) BEMT_UpdateStates() is called in loose coupling &
                                                                         !   (2) BEMT_UpdateDiscState() is called in tight coupling.
                                                                         !   Input is the suggested time from the glue code;
                                                                         !   Output is the actual coupling interval that will be used
                                                                         !   by the glue code.
   type(UA_InitOutputType),      intent(  out)  :: InitOut     ! Output for initialization routine
   integer(IntKi),               intent(  out)  :: errStat     ! Error status of the operation
   character(*),                 intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                       :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                               :: errStat2    ! temporary Error status of the operation
   integer(IntKi)                               :: i,j, iNode, iOffset
   character(64)                                :: chanPrefix
   
      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""


      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
   !call DispNVD( UA_Ver )
   
      ! Allocate and set parameter data structure using initialization data
   call UA_SetParameters( interval, InitInp, p, errStat, errMsg )
   if (errStat >= AbortErrLev) return
   
      ! initialize the discrete states
   call UA_InitStates( p, xd, OtherState, errStat, errMsg )     ! initialize the continuous states
   if (errStat >= AbortErrLev) return    
   
      ! Allocate and set the InitOut data
   p%NumOuts = 6
   p%Delim   = ''
   p%OutFmt  = 'ES11.4e2'
   p%OutSFmt = 'A11'
   allocate(InitOut%WriteOutputHdr(p%NumOuts*p%numBlades*p%nNodesPerBlade))
   
   allocate(InitOut%WriteOutputUnt(p%NumOuts*p%numBlades*p%nNodesPerBlade))
   
   iNode = 0
   do j = 1,p%numBlades
      do i = 1,p%nNodesPerBlade
         
         iOffset = (i-1)*p%NumOuts + (j-1)*p%nNodesPerBlade*p%NumOuts 
         chanPrefix = "B"//trim(num2lstr(j))//"N"//trim(num2lstr(i))  
         InitOut%WriteOutputHdr(iOffset+1) =trim(chanPrefix)//'AOA'
         InitOut%WriteOutputHdr(iOffset+2) =trim(chanPrefix)//'Cn' 
         InitOut%WriteOutputHdr(iOffset+3) =trim(chanPrefix)//'Cc' 
         InitOut%WriteOutputHdr(iOffset+4) =trim(chanPrefix)//'Cm' 
         InitOut%WriteOutputHdr(iOffset+5) =trim(chanPrefix)//'Cl' 
         InitOut%WriteOutputHdr(iOffset+6) =trim(chanPrefix)//'Cd' 
         InitOut%WriteOutputUnt(iOffset+1)   ='(deg)  '
         InitOut%WriteOutputUnt(iOffset+2) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+3) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+4) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+5) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+6) ='(-)   '
      end do
   end do
   
   allocate(y%WriteOutput(p%NumOuts*p%numBlades*p%nNodesPerBlade))
   
end subroutine UA_Init
!==============================================================================     
                              
                              
!============================================================================== 
subroutine UA_UpdateDiscState( i, j, u, p, xd, OtherState, AFInfo, ErrStat, ErrMsg )   
! Tight coupling routine for updating discrete states
!..............................................................................
   
   integer   ,                   intent(in   )  :: i           ! node index within a blade
   integer   ,                   intent(in   )  :: j           ! blade index    
   type(UA_InputType),           intent(in   )  :: u           ! Inputs at Time                       
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters                                 
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Input: Discrete states at Time; 
                                                                           ! Output: Discrete states at Time + Interval
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Other/optimization states  
   type(AFInfoType),             intent(in   )  :: AFInfo      ! The airfoil parameter data
   integer(IntKi),               intent(  out)  :: errStat     ! Error status of the operation
   character(*),                 intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

         ! Local Variables
  
   real(ReKi)                                   :: Cn_prime
   real(ReKi)                                   :: Cn1               ! critical value of Cn_prime at LE separation for alpha >= alpha0    
   real(ReKi)                                   :: Cn2               ! critical value of Cn_prime at LE separation for alpha < alpha0
   real(ReKi)                                   :: Cd0
   real(ReKi)                                   :: Cm0
   real(ReKi)                                   :: k0
   real(ReKi)                                   :: k1
   real(ReKi)                                   :: k2
   real(ReKi)                                   :: k3
   real(ReKi)                                   :: T_VL
   real(ReKi)                                   :: x_cp_bar  
   real(ReKi)                                   :: Kalpha            ! backwards finite difference of alpha
   real(ReKi)                                   :: Kq
   real(ReKi)                                   :: q_cur
   real(ReKi)                                   :: X1
   real(ReKi)                                   :: X2
   real(ReKi)                                   :: Kprime_alpha
   real(ReKi)                                   :: Kprime_q
   real(ReKi)                                   :: Dp
   real(ReKi)                                   :: Cn_pot
   real(ReKi)                                   :: Cc_pot
   real(ReKi)                                   :: Cn_alpha_q_circ
   real(ReKi)                                   :: Cm_q_circ
   real(ReKi)                                   :: Cn_alpha_nc
   real(ReKi)                                   :: Cm_q_nc
   real(ReKi)                                   :: fprimeprime
   real(ReKi)                                   :: Df
   real(ReKi)                                   :: fprime
   real(ReKi)                                   :: Cn_v              ! normal force coefficient due to the presence of LE vortex
   real(ReKi)                                   :: C_V               ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi)                                   :: Cn_FS
   real(ReKi)                                   :: alpha_e
   real(ReKi)                                   :: alpha0            ! zero lift angle of attack (radians)
   real(ReKi)                                   :: dalpha0
   real(ReKi)                                   :: eta_e
   real(ReKi)                                   :: Kafactor
   real(ReKi)                                   :: St_sh
   real(ReKi)                                   :: T_sh
            
      ! Initialize ErrStat
         
   errStat = ErrID_None         
   errMsg  = ""               

      call ComputeKelvinChain(i, j, u, p, xd, OtherState, AFInfo, Cn_prime, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                                     St_sh, Kalpha, alpha_e, alpha0, dalpha0, eta_e, Kq, q_cur, X1, X2, &
                                     Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, fprime, &
                                     Cn_alpha_q_circ, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, errStat, errMsg )
      
      
      T_sh = 2.0_ReKi*(1.0_ReKi-fprimeprime) / St_sh
      
      !---------------------------------------------------------
      ! Update the OtherStates
      !---------------------------------------------------------
   
      if ( (Cn_prime > Cn1) .or. ( Cn_prime < Cn2 ) ) then  ! assumption is that Cn2 <  0.0 and Cn1 > 0
         OtherState%LESF(i,j) = .true.  ! LE separation can occur
      else
         OtherState%LESF(i,j) = .false.
      end if
                                 
      if (fprimeprime < xd%fprimeprime_minus1(i,j)) then
         OtherState%TESF(i,j) = .true.
      else
         OtherState%TESF(i,j) = .false.
      end if
   
      ! Process VRTX-related quantities
      !!!!!!!!!!!!!!!!!!!!!
      !! NEW CODE 2/19/2015
      
      
      if ( (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) > 0.0_ReKi) ) then
         OtherState%VRTX(i,j) = .true.
      else
         OtherState%VRTX(i,j) = .false.
         
      end if
      
      if ( xd%tau_V(i,j) >= (T_VL + T_sh) ) then
         xd%tau_V(i,j)     = 0.0_ReKi
      end if
      
      OtherState%sigma1(i,j) = 1.0_ReKi
      Kafactor             = Kalpha*dalpha0
      
      if ( OtherState%TESF(i,j) ) then
         if (Kafactor < 0.0_ReKi) then
            OtherState%sigma1(i,j) = 2.0_ReKi
         else if (.not. OtherState%LESF(i,j) ) then
            OtherState%sigma1(i,j) = 1.0_ReKi
         else if (xd%fprimeprime_minus1(i,j) <= 0.7_ReKi) then
            OtherState%sigma1(i,j) = 2.0_ReKi
         else
            OtherState%sigma1(i,j) = 1.75_ReKi
         end if
      else
         if (.not. OtherState%LESF(i,j) ) then
            OtherState%sigma1(i,j) = 0.5_ReKi
         end if
         
         if ( OtherState%VRTX(i,j) .and. (xd%tau_V(i,j) <= T_VL) ) then  ! TODO: Verify 2/20/2015 GJH
            OtherState%sigma1(i,j) = 0.25_ReKi
         end if
         if (Kafactor > 0.0_ReKi) then
            OtherState%sigma1(i,j) = 0.75_ReKi
         end if
      end if
      
      OtherState%sigma3(i,j) = 1.0_ReKi
      
      if ( (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) >= T_VL) ) then
         OtherState%sigma3(i,j) = 3.0_ReKi
         if (.not. OtherState%TESF(i,j)) then
            OtherState%sigma3(i,j) = 4.0_ReKi
            if ( OtherState%VRTX(i,j) .and. (xd%tau_V(i,j) <= T_VL) ) then
               if (Kafactor < 0.0_ReKi) then
                  OtherState%sigma3(i,j) = 2.0_ReKi
               else
                  OtherState%sigma3(i,j) = 1.0_ReKi
               end if
            end if
         end if
      else if (Kafactor < 0 ) then
         OtherState%sigma3(i,j) = 4.0_ReKi
      end if
      
      
      if ((.not. OtherState%TESF(i,j)) .and. (Kq*dalpha0 < 0.0_ReKi)) then
         OtherState%sigma3(i,j) = 1.0_ReKi
      end if
   
      
     
      !---------------------------------------------------------
      ! Update the Discrete States, xd
      !---------------------------------------------------------
      if (OtherState%FirstPass) then
         xd%alpha_minus2(i,j)      = u%alpha
         xd%q_minus2(i,j)          = q_cur
      else
         xd%alpha_minus2(i,j)      = xd%alpha_minus1(i,j)
         xd%q_minus2(i,j)          = xd%q_minus1(i,j)
      end if
           
      xd%alpha_minus1(i,j)        = u%alpha
      xd%q_minus1(i,j)            = q_cur
      xd%X1_minus1(i,j)           = X1
      xd%X2_minus1(i,j)           = X2
      xd%Kprime_alpha_minus1(i,j) = Kprime_alpha
      xd%Kprime_q_minus1(i,j)     = Kprime_q
      xd%Dp_minus1(i,j)           = Dp
      xd%Cn_pot_minus1(i,j)       = Cn_pot
      xd%fprimeprime_minus1(i,j)  = fprimeprime
      xd%Df_minus1(i,j)           = Df
      xd%fprime_minus1(i,j)       = fprime
      xd%Cn_v_minus1(i,j)         = Cn_v
      xd%C_V_minus1(i,j)          = C_V
      OtherState%FirstPass      = .false.
   
      if ( xd%tau_V(i,j) > 0.0 .or. OtherState%LESF(i,j) ) then                !! TODO:  Verify this condition 2/20/2015 GJH
         xd%tau_V(i,j)          = xd%tau_V(i,j) + 2.0_ReKi*p%dt*u%U / p%c(i,j)  
      end if
   
end subroutine UA_UpdateDiscState
!==============================================================================   

!============================================================================== 
subroutine UA_UpdateStates( i, j, u, p, xd, OtherState, AFInfo, errStat, errMsg )
! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
! Continuous, constraint, and discrete states are updated to values at t + Interval.
!..............................................................................

   integer(IntKi),                intent(in   ) :: i               ! node index within a blade
   integer(IntKi),                intent(in   ) :: j               ! blade index 
   type(UA_InputType),            intent(inout) :: u               ! Input at current timestep, t
   type(UA_ParameterType),        intent(in   ) :: p               ! Parameters
   type(UA_DiscreteStateType),    intent(inout) :: xd              ! Input: Discrete states at t;
                                                                   !   Output: Discrete states at t + Interval
   type(UA_OtherStateType),       intent(inout) :: OtherState      ! Other/optimization states
   type(AFInfoType),              intent(in   ) :: AFInfo          ! The airfoil parameter data
   integer(IntKi),                intent(  out) :: errStat         ! Error status of the operation
   character(*),                  intent(  out) :: errMsg          ! Error message if ErrStat /= ErrID_None

      ! Local variables  
   
   
   integer(IntKi)                               :: errStat2        ! Error status of the operation (secondary error)
   character(len(ErrMsg))                       :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   
      ! Initialize variables

   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
         
   
   OtherState%iBladeNode = i; 
   OtherState%iBlade = j

      !BJJ: this seems to be the root cause of all sorts of numerical problems....
   IF (EqualRealNos(u%u, 0.0_ReKi) ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'UA_UpdateStates: U (air velocity magnitude relative to the airfoil) is zero.'
      RETURN
   END IF
   
   
      ! Update discrete states:
   call UA_UpdateDiscState( i, j, u, p, xd, OtherState, AFInfo, ErrStat, ErrMsg )
      
end subroutine UA_UpdateStates
!==============================================================================   

!============================================================================== 
subroutine UA_CalcOutput( u, p, xd, OtherState, AFInfo, y, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..............................................................................
   
   type(UA_InputType),           intent(inout)  :: u           ! Inputs at Time
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   )  :: OtherState  ! Other/optimization states
   type(AFInfoType),             intent(in   )  :: AFInfo      ! The airfoil parameter data
   type(UA_OutputType),          intent(inout)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                     !   nectivity information does not have to be recalculated)
   integer(IntKi),               intent(  out)  :: errStat     ! Error status of the operation
   character(*),                 intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

   
   integer(IntKi)                                         :: errStat2        ! Error status of the operation (secondary error)
   character(len(errMsg))                                 :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   real(ReKi)                                             :: Cn_prime
   real(ReKi)                                             :: Cn1             ! critical value of Cn_prime at LE separation for alpha >= alpha0    
   real(ReKi)                                             :: Cn2             ! critical value of Cn_prime at LE separation for alpha < alpha0
   real(ReKi)                                             :: Cd0
   real(ReKi)                                             :: Cm0
   real(ReKi)                                             :: k0
   real(ReKi)                                             :: k1
   real(ReKi)                                             :: k2
   real(ReKi)                                             :: k3
   real(ReKi)                                             :: T_VL
   real(ReKi)                                             :: x_cp_bar 
   real(ReKi)                                             :: Kalpha             ! backwards finite difference of alpha
   real(ReKi)                                             :: Kq
   real(ReKi)                                             :: q_cur
   real(ReKi)                                             :: X1
   real(ReKi)                                             :: X2
   real(ReKi)                                             :: Kprime_alpha
   real(ReKi)                                             :: Kprime_q
   real(ReKi)                                             :: Dp
   real(ReKi)                                             :: Cn_pot
   real(ReKi)                                             :: Cc_pot
   real(ReKi)                                             :: Cn_alpha_q_circ
   real(ReKi)                                             :: Cm_q_circ
   real(ReKi)                                             :: Cn_alpha_nc
   real(ReKi)                                             :: Cm_q_nc
   real(ReKi)                                             :: fprimeprime
   real(ReKi)                                             :: Df
   real(ReKi)                                             :: fprime
   real(ReKi)                                             :: Cn_v                ! normal force coefficient due to the presence of LE vortex
   real(ReKi)                                             :: C_V                 ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi)                                             :: Cn_FS  
   real(ReKi)                                             :: alpha_e
   real(ReKi)                                             :: alpha0              ! zero lift angle of attack (radians)
   real(ReKi)                                             :: dalpha0
   real(ReKi)                                             :: eta_e
   real(ReKi)                                             :: St_sh
   real(ReKi)                                             :: Cl, Cd
   integer                                                :: iOffset
   
   !BJJ: what are OtherState%iBladeNode, OtherState%iBlade here? I don't see them set in the routine that calls UA_CalcOutput
   
   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
   
  
   iOffset = (OtherState%iBladeNode-1)*p%NumOuts + (OtherState%iBlade-1)*p%nNodesPerBlade*p%NumOuts 
   
   if (OtherState%FirstPass) then
      
      
      call GetSteadyOutputs(AFInfo, u%alpha, y%Cl, y%Cd, y%Cm, Cd0, errStat, errMsg)
      
      y%Cn = y%Cl*cos(u%alpha) + (y%Cd-Cd0)*sin(u%alpha)
      y%Cc = y%Cl*sin(u%alpha) - (y%Cd-Cd0)*cos(u%alpha)
      y%WriteOutput(iOffset+1) = u%alpha*180.0/pi
      y%WriteOutput(iOffset+2) = y%Cn
      y%WriteOutput(iOffset+3) = y%Cc
      y%WriteOutput(iOffset+4) = y%Cm
      y%WriteOutput(iOffset+5) = y%Cl
      y%WriteOutput(iOffset+6) = y%Cd
      
   else
      
            
      call ComputeKelvinChain( OtherState%iBladeNode, OtherState%iBlade, u, p, xd, OtherState, AFInfo, Cn_prime, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                                 St_sh, Kalpha, alpha_e, alpha0, dalpha0, eta_e, Kq, q_cur, X1, X2, &
                                 Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, fprime, &
                                 Cn_alpha_q_circ, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, errStat, errMsg )
         ! Eqn 1.50
      y%Cn= Cn_FS + Cn_v
            
         ! Eqn 1.52
      y%Cc = eta_e*Cc_pot*sqrt(fprimeprime) + Cn_v*tan(alpha_e)*(1-xd%tau_v(OtherState%iBladeNode, OtherState%iBlade))
            
         ! Eqn 1.2
      y%Cl = y%Cn*cos(u%alpha) + y%Cc*sin(u%alpha)
      y%Cd = y%Cn*sin(u%alpha) - y%Cc*cos(u%alpha) + Cd0
            
         ! Eqn 1.55
      y%Cm = Get_Cm( Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, Cn_alpha_q_circ, fprimeprime, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, xd%tau_v(OtherState%iBladeNode, OtherState%iBlade) )
            
      y%WriteOutput(iOffset+1)   = u%alpha*180.0/pi  ! Angle of attack in degrees
      y%WriteOutput(iOffset+2) = y%Cn
      y%WriteOutput(iOffset+3) = y%Cc
      y%WriteOutput(iOffset+4) = y%Cm
      y%WriteOutput(iOffset+5) = y%Cl
      y%WriteOutput(iOffset+6) = y%Cd
      
   end if
   
end subroutine UA_CalcOutput
!==============================================================================   

 !integer function Test_Cn()
 !
 !  ! Driver constants /  parameters
 !  real(DbKi)  :: dt                = 1.0_DbKi
 !  real(ReKi)  :: a_s               = 1000.0_ReKi  ! speed of sound
 !  
 !  ! Inputs
 !  real(ReKi)  :: alpha_cur
 !  real(ReKi)  :: q_cur
 !  real(ReKi)  :: M  
 !  real(ReKi)  :: U
 !  
 !  
 !  ! States
 !  real(ReKi)  :: alpha_minus1
 !  real(ReKi)  :: alpha_minus2
 !  real(ReKi)  :: q_minus1
 !  real(ReKi)  :: q_minus2
 !  real(ReKi)  :: X_1_minus1
 !  real(ReKi)  :: X_2_minus1
 !  real(ReKi)  :: Kprime_alpha_minus1
 !  real(ReKi)  :: Kprime_q_minus1
 !  real(ReKi)  :: Dp_minus1
 !  real(ReKi)  :: Cn_pot_minus1
 !  real(ReKi)  :: Df_minus1
 !  real(ReKi)  :: fprime_minus1
 !  
 !  
 !  
 !  ! Parameters
 !  real(ReKi)  :: c   
 !  
 !  real(ReKi)  :: T_p           ! empirical, specific to a given Airfoil ??  
 !  real(ReKi)  :: T_f           ! function of M, Re, and specific for a given Airfoil
 !  real(ReKi)  :: C_nalpha_circ ! function of s, M, and airfoil
 !  real(ReKi)  :: b1
 !  real(ReKi)  :: b2
 !  real(ReKi)  :: A1
 !  real(ReKi)  :: A2
 !  real(ReKi)  :: C_nalpha
 !  real(ReKi)  :: alpha1
 !  real(ReKi)  :: alpha0
 !  real(ReKi)  :: S1
 !  real(ReKi)  :: S2
 !  real(ReKi)  :: beta_M_sqrd
 !  real(ReKi)  :: T_I
 !  integer     :: DSMod
 !
 !  
 !  
 !  ! Local derived quantities
 !  
 !  real(ReKi)  :: ds
 !  real(ReKi)  :: T_alpha       
 !  real(ReKi)  :: T_q           
 !  real(ReKi)  :: dalpha
 !  real(ReKi)  :: Cn
 !  
 !  real(ReKi)  :: T_alphaVer
 !  
 !  ! Set Parameters
 !  c             = 1.0_ReKi
 !  T_p           = 1.0_ReKi
 !  T_f           = 1.0_ReKi
 !  C_nalpha_circ = 1.0_ReKi
 !  b1            = 1.0_ReKi
 !  b2            = 1.0_ReKi
 !  A1            = 1.0_ReKi
 !  A2            = 1.0_ReKi
 !  C_nalpha      = 1.0_ReKi
 !  alpha1        = 1.0_ReKi
 !  alpha0        = 1.0_ReKi
 !  S1            = 1.0_ReKi
 !  S2            = 1.0_ReKi
 !  DSMod         = 1
 !  
 !  ! Initialize states
 !  alpha_minus1    = 0.0_ReKi
 !  alpha_minus2  = 0.0_ReKi
 !  q_minus1        = 0.0_ReKi
 !  q_minus2      = 0.0_ReKi
 !  X_1_minus1      = 0.0_ReKi
 !  X_2_minus1      = 0.0_ReKi
 !  Kprime_alpha_minus1 = 0.0_ReKi
 !  Kprime_q_minus1 = 0.0_ReKi
 !  Dp_minus1       = 0.0_ReKi
 !  Cn_pot_minus1   = 0.0_ReKi
 !  Df_minus1       = 0.0_ReKi
 !  fprime_minus1   = 0.0_ReKi
 !  
 !  
 !  ! set inputs
 !  alpha_cur     = 0.1_ReKi
 !  
 !  
 !  U             = 10.0_ReKi
 !  M             = U / a_s
 !  T_I           = c / a_s
 !  beta_M_sqrd   = Get_Beta_M_Sqrd(M)
 !  q_cur         = Get_Pitchrate( alpha_cur, alpha_minus1, U, c )
 !  dalpha   = alpha_cur - alpha_minus1
 !  ds      = Get_ds(U, c, dt)
 !  !T_alpha = Get_T_Alpha( M, a_s, c, C_nalpha )
 !  T_alphaVer = 1/ ( (1-M) + .5*M**2*(1-M**2)*.413)* c / a_s
 !  T_q     = Get_T_q( M, beta_M_sqrd, T_I, a_s, c, C_nalpha) 
 !  
 !  !Cn = Get_Cn( dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
 !  !            q_minus1, q_minus2, T_alpha, T_q, T_p, T_f, X_1_minus1, X_2_minus1, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, Df_minus1, fprime_minus1, b1, b2, A1, A2, C_nalpha, &
 !  !            alpha1, S1, S2, DSMod )
 !              !dt, ds, U, M, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
 !              !              q_minus1, q_minus2, T_alpha, T_q, T_p, T_f, X_1_minus1, X_2_minus1, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, Df_minus1, fprime_minus1, b1, b2, A1, A2, C_nalpha, &
 !              !              alpha1, S1, S2, DSMod )
 !  
 !  Test_Cn = 1
 !  
 !end function Test_Cn
 !  
   
   
   
   
end module UnsteadyAero