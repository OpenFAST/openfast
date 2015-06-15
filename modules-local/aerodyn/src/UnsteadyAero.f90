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
   IntAFCoefs(1:s1) = CubicSplineInterpM( 1.0*real( AOA*180.0/PI ) &
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
real(ReKi) function Get_k_( M, beta_M, C_nalpha, A1, A2, b1, b2          )
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: M               ! Mach number
   real(ReKi), intent(in   ) :: beta_M          ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2)
   real(ReKi), intent(in   ) :: C_nalpha        ! slope of the 2D pitching moment curve 
   real(ReKi)                :: A1              ! airfoil constant derived from experimental results, usually 0.3
   real(ReKi)                :: A2              ! airfoil constant derived from experimental results, usually 0.7
   real(ReKi)                :: b1              ! airfoil constant derived from experimental results, usually 0.14
   real(ReKi)                :: b2              ! airfoil constant derived from experimental results, usually 0.53

      ! Implements equation 1.11a/b
   Get_k_   = 1.0_ReKi / ( (1.0_ReKi - M) + C_nalpha * M**2 * beta_M * (A1*b1 + A2*b2) / 2.0  )
   !Get_k_   = 1.0_ReKi / ( (1.0_ReKi - M) + C_nalpha * M**2 * beta_M *  0.413_ReKi  )
   
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
real(ReKi) function Get_Cn_q_circ(q_cur, C_nalpha_circ, X3, X4)
! Called by : ComputeKelvinChain
! Calls  to : NONE   
!..............................................................................
   real(ReKi), intent(in   )     :: q_cur
   real(ReKi), intent(in   )     :: C_nalpha_circ           ! slope of the circulatory normal force coefficient vs alpha curve
   real(ReKi), intent(in   )     :: X3                      ! Exponential decay function associated with Cn_q_circ
   real(ReKi), intent(in   )     :: X4                      ! Exponential decay function associated with Cn_q_circ

      ! Implements eqn 1.16a
   Get_Cn_q_circ = C_nalpha_circ*q_cur/2.0 - X3 - X4
   
end function Get_Cn_q_circ
!==============================================================================


!==============================================================================
real(ReKi) function Get_Cm_q_nc( UAMod, M, k_mq, T_I, Cn_q_nc, k_alpha, Kq, Kprimeprime_q )
! Called by : ComputeKelvinChain
! Calls  to : NONE   
!..............................................................................

   integer,    intent(in   ) :: UAMod                 ! UA model
   real(ReKi), intent(in   ) :: M                     ! Mach number
   real(ReKi), intent(in   ) :: k_mq                  ! airfoil parameter                     
   real(ReKi), intent(in   ) :: T_I                   ! time constant        
   real(ReKi), intent(in   ) :: Cn_q_nc               !    
   real(ReKi), intent(in   ) :: k_alpha               !     
   real(ReKi), intent(in   ) :: Kq                    ! backward finite difference of q    
   real(ReKi), intent(in   ) :: Kprimeprime_q         ! deficiency function             
    
   if ( UAMod == 3 ) then
      ! Implements eqn 1.27
      Get_Cm_q_nc =  -Cn_q_nc/4.0 - (k_alpha**2)*T_I*(Kq-Kprimeprime_q)/(3.0*M)
   else  
         ! Implements eqn 1.25a
      Get_Cm_q_nc = -7*(k_mq**2)*T_I*(Kq-Kprimeprime_q)/(12.0*M)
      
   end if
   
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
real(ReKi) function Get_f_from_Lookup( UAMod, Re, alpha, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
! Compute either fprime or fprimeprime using an analytical equation (and eventually a table lookup)
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   integer,          intent(in   ) :: UAMod
   real(ReKi),       intent(in   ) :: Re            ! Reynolds number
   real(ReKi),       intent(in   ) :: alpha         ! angle of attack (radians)
   real(ReKi),       intent(in   ) :: alpha0
   real(ReKi),       intent(in   ) :: C_nalpha
   type(AFInfoType), intent(in   ) :: AFInfo        ! The airfoil parameter data
   integer(IntKi),   intent(  out) :: ErrStat               ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg                ! Error message if ErrStat /= ErrID_None
   
   !real                            :: IntAFCoefs(4)         ! The interpolated airfoil coefficients.
   real                            :: Cn, Cl, Cd, Cm, Cd0
   !integer                         :: s1                    ! Number of columns in the AFInfo structure
   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
      ! TODO: Extend this to use the AFInfo tables to directly look up f based on either Cn or Cc (or Cm?).
      
      
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   !s1 = size(AFInfo%Table(1)%Coefs,2)
   !if (s1 < 4) then
   !   ErrMsg  = 'The Airfoil info table must contains columns for lift, drag, and pitching moment, and f'
   !   ErrStat = ErrID_Fatal
   !   return
   !end if
   !
   !IntAFCoefs(1:s1) = CubicSplineInterpM( 1.0*real( alpha*180.0/PI ) &
   !                                          , AFInfo%Table(1)%Alpha &
   !                                          , AFInfo%Table(1)%Coefs &
   !                                          , AFInfo%Table(1)%SplineCoefs &
   !                                          , ErrStat, ErrMsg )
   !if (ErrStat > ErrID_None) return
  
   !Get_f_from_Lookup = IntAFCoefs(4)   
   
   if (abs(alpha-alpha0) < .01) then
      Get_f_from_Lookup = 1.0_ReKi
      return
   end if
   
   call GetSteadyOutputs(AFInfo, alpha, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
      if (ErrStat > ErrID_None) return
   
   Cn =  Cl*cos(alpha) + (Cd-Cd0)*sin(alpha)
   
   if (UAMod == 2) then
      Get_f_from_Lookup = ((3*sqrt(Cn/(C_nalpha*(alpha-alpha0)))-1)/2.0)**2
   else
      Get_f_from_Lookup = ( 2 * sqrt( Cn / ( C_nalpha*( alpha-alpha0 ) ) ) - 1 ) **2 
   end if
      
   if ( Get_f_from_Lookup > 1.0 ) then
      Get_f_from_Lookup = 1.0_ReKi
   end if
   
end function Get_f_from_Lookup      


!==============================================================================
real(ReKi) function Get_f_c_from_Lookup( Re, alpha, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
! Compute either fprime or fprimeprime using an analytical equation (and eventually a table lookup)
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi),       intent(in   ) :: Re            ! Reynolds number
   real(ReKi),       intent(in   ) :: alpha         ! angle of attack (radians)
   real(ReKi),       intent(in   ) :: alpha0
   real(ReKi),       intent(in   ) :: C_nalpha
   type(AFInfoType), intent(in   ) :: AFInfo        ! The airfoil parameter data
   integer(IntKi),   intent(  out) :: ErrStat               ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg                ! Error message if ErrStat /= ErrID_None
   
   !real                            :: IntAFCoefs(4)         ! The interpolated airfoil coefficients.
   real                            :: Cc, Cl, Cd, Cm, Cd0
   !integer                         :: s1                    ! Number of columns in the AFInfo structure
   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
      ! TODO: Extend this to use the AFInfo tables to directly look up f based on either Cn or Cc (or Cm?).
      
      
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   !s1 = size(AFInfo%Table(1)%Coefs,2)
   !if (s1 < 4) then
   !   ErrMsg  = 'The Airfoil info table must contains columns for lift, drag, and pitching moment, and f'
   !   ErrStat = ErrID_Fatal
   !   return
   !end if
   !
   !IntAFCoefs(1:s1) = CubicSplineInterpM( 1.0*real( alpha*180.0/PI ) &
   !                                          , AFInfo%Table(1)%Alpha &
   !                                          , AFInfo%Table(1)%Coefs &
   !                                          , AFInfo%Table(1)%SplineCoefs &
   !                                          , ErrStat, ErrMsg )
   !if (ErrStat > ErrID_None) return
  
   !Get_f_from_Lookup = IntAFCoefs(4)   
   
   call GetSteadyOutputs(AFInfo, alpha, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
      if (ErrStat > ErrID_None) return
   
   Cc =  Cl*sin(alpha) - (Cd-Cd0)*cos(alpha)
      ! Apply an offset of 0.2 to fix cases where f_c should be negative, but we are using **2 so can only return positive values
      ! Note: because we include this offset, it must be accounted for in the final value of Cc, eqn 1.38.  This will be applied
      ! For both UA_Mod = 1,2, and 3 when using Flookup = T
   Get_f_c_from_Lookup = (  Cc / ( C_nalpha*( alpha-alpha0 )*tan(alpha) )  + 0.2 ) **2 
   
   !Get_f_c_from_Lookup = (  Cc / ( C_nalpha*( alpha-alpha0 )*tan(alpha) )  ) **2 
   !if ( Cc < 0.0_ReKi ) then
   !   Get_f_c_from_Lookup = -Get_f_c_from_Lookup
   !end if
   
   if ( Get_f_c_from_Lookup > 1.0 ) then
      Get_f_c_from_Lookup = 1.0_ReKi
   end if
   
end function Get_f_c_from_Lookup      

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
real(ReKi) function Get_Cn_FS( Cn_alpha_q_nc, Cn_alpha_q_circ, fprimeprime, UAMod )
!  
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

    
   real(ReKi), intent(in   ) :: Cn_alpha_q_nc       ! non-circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ     ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime         ! lagged version of fprime accounting for unsteady boundary layer response
   integer   , intent(in   ) :: UAMod               ! model for the dynamic stall equations.  [1 = Leishman/Beddoes, 2 = Gonzalez, 3 = Minnema]
   
   if ( UAMod == 2 ) then
      ! use equation 1.36
      Get_Cn_FS   = Cn_alpha_q_nc + Cn_alpha_q_circ *  ( (1.0_ReKi + 2.0_ReKi*sqrt(fprimeprime) ) / 3.0_ReKi  )**2
   else
      ! use equation 1.35
      Get_Cn_FS   = Cn_alpha_q_nc + Cn_alpha_q_circ * ( (1.0_ReKi + sqrt(fprimeprime) ) / 2.0_ReKi )**2
   
   end if


end function Get_Cn_FS
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cc_FS( eta_e, alpha_e, alpha0, Cn_alpha_q_circ, fprimeprime, UAMod )
!  
! Called by : NONE
! Calls  to : NONE
!..............................................................................
   
   real(ReKi), intent(in   ) :: eta_e
   real(ReKi), intent(in   ) :: alpha_e                    ! effective angle of attack at 3/4 chord  TODO: verify 3/4 and not 1/4
   real(ReKi), intent(in   ) :: alpha0                     ! zero lift angle of attack (radians)
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ            ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime                ! lagged version of fprime accounting for unsteady boundary layer response
   integer   , intent(in   ) :: UAMod                      ! model for the dynamic stall equations.  [1 = Leishman/Beddoes, 2 = Gonzalez, 3 = Minnema]
   
   if ( UAMod == 1 ) then
      ! use equation 1.37
      Get_Cc_FS   = Cn_alpha_q_circ * sqrt(fprimeprime) * tan(alpha_e + alpha0) * eta_e
   else if ( UAMod == 2 ) then
      ! use equation 1.38
      Get_Cc_FS   = Cn_alpha_q_circ * (sqrt(fprimeprime) - 0.2_ReKi) * tan(alpha_e + alpha0) * eta_e
   else
      ! implementation error!  Cannot have UAMod other than 0,1,2
   end if

end function Get_Cc_FS
!==============================================================================   


!==============================================================================
! Dynamic Stall Equations                                                     !
!==============================================================================

!==============================================================================
real(ReKi) function Get_C_V( Cn_alpha_q_circ, fprimeprime, UAMod )
! 
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: Cn_alpha_q_circ               ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi), intent(in   ) :: fprimeprime                   ! lagged version of fprime accounting for unsteady boundary layer response
   integer   , intent(in   ) :: UAMod                         ! model for the dynamic stall equations.  [1 = Leishman/Beddoes, 2 = Gonzalez, 3 = Minnema]
   
   if ( UAMod == 2 ) then
      ! use equation 1.48
      Get_C_V   = Cn_alpha_q_circ * ( 1.0_ReKi - (1.0_ReKi + 2.0_ReKi*sqrt(fprimeprime) )/3.0_ReKi  )**2
   else
      ! use equation 1.47
      Get_C_V   = Cn_alpha_q_circ *  ( 1.0_ReKi - ( 0.5_ReKi + 0.5_ReKi*sqrt(fprimeprime) ) )**2
   
   end if
   

end function Get_C_V
!==============================================================================   

!==============================================================================
real(ReKi) function Get_Cn_v ( ds, T_V, Cn_v_minus1, C_V, C_V_minus1, tau_V, T_VL, T_V0, Kalpha, alpha, alpha0 ) ! do I pass VRTX flag or tau_V?
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
   real(ReKi), intent(in   ) :: T_V0                   !
   real(ReKi), intent(in   ) :: Kalpha                 ! backwards finite difference of alpha
   real(ReKi), intent(in   ) :: alpha                  ! angle of attack (rad)
   real(ReKi), intent(in   ) :: alpha0                 ! zero lift angle of attack (rad)
  
   real(ReKi)                :: factor
   
   factor = (alpha - alpha0) * Kalpha
   
   if (tau_V > T_VL .AND. factor >= 0.0_ReKi ) then  
      !TODO: Need to examine eqn 1.49.  The manual states T_V0 / sigma3 and forces sigma3 = 2.0 ! GJH 5/28/2015
      Get_Cn_v = Cn_v_minus1*exp(-ds/T_V)   ! Eqn 1.49
      !Get_Cn_v = Cn_v_minus1*exp(-2.0*ds/T_V0)   ! Eqn 1.49  alternative using T_V0 and forcing sigma3 for this eqn = 2.0
   else      
      Get_Cn_v = Get_ExpEqn( ds, T_V, Cn_v_minus1, C_V, C_V_minus1 )   ! Eqn 1.45
   end if
   
end function Get_Cn_v
!==============================================================================   


!==============================================================================
real(ReKi) function Get_Cm_FS( Cm0, k0, k1, k2, k3, Cn_alpha_q_circ, fprimeprime, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Dalphaf, alpha_f, Cn_FS, fprimeprime_m, AFInfo, UAMod, Cm_alpha_nc, ErrStat, ErrMsg )
! Leishman's pitching moment coefficient about the 1/4 chord
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   
   real(ReKi),       intent(in   ) :: Cm0                           ! 2D pitching moment coefficient at zero lift, positive if nose is up
   real(ReKi),       intent(in   ) :: k0                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(in   ) :: k1                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(in   ) :: k2                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(in   ) :: k3                            ! airfoil parameter in the x_cp_hat curve best-fit
   real(ReKi),       intent(in   ) :: Cn_alpha_q_circ               ! circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi),       intent(in   ) :: fprimeprime                   ! lagged version of fprime accounting for unsteady boundary layer response
   real(ReKi),       intent(in   ) :: Cm_q_circ                     ! circulatory component of the pitching moment coefficient response to a step change in q
   real(ReKi),       intent(in   ) :: Cn_alpha_nc                   ! non-circulatory component of the normal force coefficient response to step change in alpha
   real(ReKi),       intent(in   ) :: Cm_q_nc                       ! non-circulatory component of the moment coefficient response to step change in q
   real(ReKi),       intent(in   ) :: Dalphaf
   real(ReKi),       intent(in   ) :: alpha_f
   real(ReKi),       intent(in   ) :: Cn_FS
   real(ReKi),       intent(in   ) :: fprimeprime_m
   type(AFInfoType), intent(in   ) :: AFInfo
   integer,          intent(in   ) :: UAMod                         ! UA model
   real(ReKi),       intent(  out) :: Cm_alpha_nc
   integer(IntKi),   intent(  out) :: ErrStat           ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg            ! Error message if ErrStat /= ErrID_None
  
      ! local variables
   real(ReKi)                :: x_cp_hat                      ! center-of-pressure distance from LE in chord fraction
   real(ReKi)                :: Cm_common                     !
   real(ReKi)                :: alpha_prime_f
   real(ReKi)                :: Cl
   real(ReKi)                :: Cd
   real(ReKi)                :: Cm
   real(ReKi)                :: Cd0

   
      ! Eqn 1.21 + 1.23 + 1.25a
   Cm_alpha_nc = - Cn_alpha_nc / 4.0_ReKi 
   Cm_common = Cm_q_circ + Cm_alpha_nc + Cm_q_nc
   
   if ( UAMod == 1 ) then
         ! Eqn 1.40
      x_cp_hat = k0 + k1*(1-fprimeprime) + k2*sin(pi*fprimeprime**k3) + Cm_common
         ! Eqn 1.39
      Get_Cm_FS  = Cm0 - Cn_alpha_q_circ*(x_cp_hat - 0.25_ReKi)
      
   elseif ( UAMod == 3 ) then
      
         ! Eqn 1.41a
      alpha_prime_f = alpha_f - Dalphaf
         ! Look up Cm using alpha_prime_f
      call GetSteadyOutputs(AFInfo, alpha_prime_f, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
      Get_Cm_FS = Cm + Cm_common
      
   else ! UAMod == 2
         ! TODO: Incomplete because we are not computing fprimeprime_m yet. GJH 5/21/2015
      Get_Cm_FS = Cm0 + Cn_FS*fprimeprime_m + Cm_common
   end if
   
end function Get_Cm_FS



!==============================================================================
real(ReKi) function Get_Cm( Cm_FS, T_VL, x_cp_bar, Cn_v, tau_v )
! Leishman's pitching moment coefficient about the 1/4 chord
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   real(ReKi), intent(in   ) :: Cm_FS                         ! 
   real(ReKi), intent(in   ) :: T_VL                          ! time variable associated with the vortex advection process; it represents the non-dimensional time in semi-chords needed for a vortex to travel from leading edge to trailing edge
   real(ReKi), intent(in   ) :: x_cp_bar                      ! airfoil parameter for calulating x_cp_v
   real(ReKi), intent(in   ) :: Cn_v                          ! normal force coefficient due to the presence of LE vortex
   real(ReKi), intent(in   ) :: tau_v                         ! time variable, tracking the travel of the LE vortex over the airfoil suction surface (-)
   
      ! local variables
   
   real(ReKi)                :: Cm_v                          ! pitching moment coefficient due to the presence of LE vortex
   
   
      ! Eqn 1.54
   Cm_v     = -x_cp_bar*( 1-cos( pi*tau_v/T_VL ) )*Cn_v
   
      ! Eqn 1.55 = (1.39 or 1.42 or 1.43) +  1.54  
   Get_Cm   = Cm_FS + Cm_v
   
   
end function Get_Cm
!==============================================================================   


!==============================================================================                              
subroutine ComputeKelvinChain( i, j, u, p, xd, OtherState, AFInfo, Cn_prime, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                               St_sh, Kalpha, alpha_e, alpha0, dalpha0, alpha_f, eta_e, Kq, q_cur, X1, X2, X3, X4, &
                               Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, Df_c, Dalphaf, fprime, fprime_c, fprimeprime_c, fprimeprime_m, &
                               Cn_alpha_q_circ, Cn_q_circ, Cn_q_nc, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, T_f, T_V, ErrStat, ErrMsg )
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
   real(ReKi),                             intent(  out) :: alpha_f           !
   real(ReKi),                             intent(  out) :: eta_e             !
   real(ReKi),                             intent(  out) :: Kq                !
   real(ReKi),                             intent(  out) :: q_cur             !
   real(ReKi),                             intent(  out) :: X1                !
   real(ReKi),                             intent(  out) :: X2                !
   real(ReKi),                             intent(  out) :: X3                !
   real(ReKi),                             intent(  out) :: X4                !
   real(ReKi),                             intent(  out) :: Kprime_alpha      !
   real(ReKi),                             intent(  out) :: Kprime_q          !
   real(ReKi),                             intent(  out) :: Dp                !
   real(ReKi),                             intent(  out) :: Cn_pot            !
   real(ReKi),                             intent(  out) :: Cc_pot            !
   real(ReKi),                             intent(  out) :: Cn_alpha_q_circ   !
   real(ReKi),                             intent(  out) :: Cm_q_circ         !
   real(ReKi),                             intent(  out) :: Cn_alpha_nc       ! non-circulatory component of the normal force coefficient response to step change in alpha
   real(ReKi),                             intent(  out) :: Cn_q_circ
   real(ReKi),                             intent(  out) :: Cn_q_nc
   real(ReKi),                             intent(  out) :: Cm_q_nc           ! non-circulatory component of the moment coefficient response to step change in q
   real(ReKi),                             intent(  out) :: fprimeprime       !
   real(ReKi),                             intent(  out) :: Df                !
   real(ReKi),                             intent(  out) :: Df_c              !
   real(ReKi),                             intent(  out) :: Dalphaf           !
   real(ReKi),                             intent(  out) :: fprime            !
   real(ReKi),                             intent(  out) :: fprime_c          !
   real(ReKi),                             intent(  out) :: fprimeprime_c     !
   real(ReKi),                             intent(  out) :: fprimeprime_m     !
   real(ReKi),                             intent(  out) :: Cn_v              ! normal force coefficient due to the presence of LE vortex
   real(ReKi),                             intent(  out) :: C_V               ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi),                             intent(  out) :: Cn_FS             !
   real(ReKi),                             intent(  out) :: T_f                                           !
   real(ReKi),                             intent(  out) :: T_V                                           ! backwards finite difference of the non-dimensionalized distance parameter
                                                                              !
   integer(IntKi),                         intent(  out) :: ErrStat           ! Error status of the operation
   character(*),                           intent(  out) :: ErrMsg            ! Error message if ErrStat /= ErrID_None


   real(ReKi)                :: ds                                            ! non-dimensionalized distance parameter
   real(ReKi)                :: M                                             ! Mach number (-)
   real(ReKi)                :: beta_M                                        ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2)
   real(ReKi)                :: beta_M_Sqrd                                   ! square of the Prandtl-Glauert compressibility correction factor,  (1-M**2)
   real(ReKi)                :: Cn_alpha_q_nc                                 ! non-circulatory component of normal force coefficient response to step change in alpha and q
   real(ReKi)                :: k_alpha                                       !
   real(ReKi)                :: T_alpha                                       !
   real(ReKi)                :: T_q                                           !
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
   real(ReKi)                :: fprime_c_minus1
   real(ReKi)                :: alphaf_minus1
   
   

   M           = u%U / p%a_s
   beta_M_Sqrd = Get_Beta_M_Sqrd(M)
   beta_M      = Get_Beta_M(M)
   T_I         = p%c(i,j) / p%a_s                                                  ! Eqn 1.11c
   ds          = Get_ds       ( u%U, p%c(i,j), p%dt )                              ! Eqn 1.5b
   
   ! Lookup values using Airfoil Info module
   call AFI_GetAirfoilParams( AFInfo, M, u%Re, u%alpha, alpha0, alpha1, alpha2, eta_e, C_nalpha, C_nalpha_circ, &
                              T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, b5, A1, A2, A5, S1, S2, S3, S4, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, k1_hat, x_cp_bar, ErrMsg, ErrStat )           
   
      ! Override eta_e if we are using Flookup
   if ( p%Flookup ) then
      eta_e = 1.0
   end if
   
   if (OtherState%FirstPass(i,j)) then
      alpha_minus1 = u%alpha
      alpha_minus2 = u%alpha
   else
      alpha_minus1 = xd%alpha_minus1(i,j)
      alpha_minus2 = xd%alpha_minus2(i,j)
   end if
   
   dalpha0  = u%alpha - alpha0
   q_cur    = Get_Pitchrate( p%dt, u%alpha, alpha_minus1, u%U, p%c(i,j)  )  ! Eqn 1.8
   
   if (OtherState%FirstPass(i,j)) then
      q_minus1 = q_cur   
      q_minus2 = q_cur 
   else
      q_minus1 = xd%q_minus1(i,j)   
      q_minus2 = xd%q_minus2(i,j) 
   end if
   
   k_alpha  = Get_k_       ( M, beta_M, C_nalpha/2.0_ReKi, A1, A2, b1, b2 )                 ! Eqn 1.11a
   k_q      = Get_k_       ( M, beta_M, C_nalpha         , A1, A2, b1, b2 )                 ! Eqn 1.11b
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
   if ( p%UAMod == 2 ) then
         ! Compute X3 and X4 using Eqn 1.16b  and then add Cn_q_circ to the previously computed Cn_alpha_q_circ
      X3              = Get_ExpEqn( ds*beta_M_Sqrd*b1, 1.0_ReKi, xd%X3_minus1(i,j), A1*(q_cur - q_minus1), 0.0_ReKi )
      X4              = Get_ExpEqn( ds*beta_M_Sqrd*b2, 1.0_ReKi, xd%X4_minus1(i,j), A2*(q_cur - q_minus1), 0.0_ReKi )
      Cn_q_circ       = Get_Cn_q_circ( q_cur, C_nalpha_circ, X3, X4  )
      Cn_alpha_q_circ = Cn_alpha_q_circ + Cn_q_circ
   else
      Cn_q_circ       = 0.0
   end if
   
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
   Cm_q_nc         = Get_Cm_q_nc( p%UAMod, M, k_mq, T_I, Cn_q_nc, k_alpha, Kq, Kprimeprime_q )
   
   
      ! Compute Cc_pot using eqn 1.28
   Cc_pot          = Cn_alpha_q_circ*tan(alpha_e+alpha0)
   
   
   if (OtherState%FirstPass(i,j)) then
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
   if (p%flookup) then
      fprime        = Get_f_from_Lookup( p%UAMod, u%Re, alpha_f, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
   else   
      fprime        = Get_f( alpha_f, alpha0, alpha1, alpha2, S1, S2, S3, S4)
   end if
   
   if (OtherState%FirstPass(i,j)) then
      fprime_minus1 = fprime
   else
      fprime_minus1 = xd%fprime_minus1(i,j)
   end if
   
      ! Compute Df using Eqn 1.33b   
   Df            = Get_ExpEqn( ds, T_f, xd%Df_minus1(i,j), fprime, fprime_minus1 )
   
      ! Compute fprimeprime using Eqn 1.33a
   fprimeprime   = fprime - Df
   
   if (p%Flookup) then
         ! Compute fprime using Eqn 1.30 and Eqn 1.31
      !if (p%UAMod == 1) then
      !   fprime_c   = Get_f_from_Lookup( p%UAMod, u%Re, alpha_f, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
      !else   
         fprime_c   = Get_f_c_from_Lookup( u%Re, alpha_f, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
      !end if
      
   
      if (OtherState%FirstPass(i,j)) then
         fprime_c_minus1 = fprime_c
         alphaf_minus1   = alpha_f
      else
         fprime_c_minus1 = xd%fprime_c_minus1(i,j)
         alphaf_minus1   = xd%alphaf_minus1(i,j)
      end if
   
         ! Compute Df using Eqn 1.33b   
      Df_c            = Get_ExpEqn( ds, T_f, xd%Df_c_minus1(i,j), fprime_c, fprime_c_minus1 )
   
         ! Compute fprimeprime using Eqn 1.33a
      fprimeprime_c   = fprime_c - Df_c
   else
      fprimeprime_c   = fprimeprime
   end if
   
      ! Compute Cn_FS using Eqn 1.35 or 1.36 depending on option selected
   Cn_FS         = Get_Cn_FS( Cn_alpha_q_nc, Cn_alpha_q_circ, fprimeprime, p%UAMod )
   
   if ( p%UAMod == 3 ) then
         ! Eqn 1.41b
      Dalphaf    = Get_ExpEqn( ds, T_f, xd%Dalphaf_minus1(i,j), alpha_f, alphaf_minus1 )
   else
      Dalphaf    = 0.0_ReKi
   end if
   
   
      ! Compute Cc_FS usign Eqn 1.37 or 1.38, but this is not used in the final value of Cc and does not contribute to any states!  
      ! So we will leave it commented out. GJH 2/23/2015
   !Cc_FS         = Get_Cc_FS( eta_e, alpha_e, alpha0, Cn_alpha_q_circ, fprimeprime, p%UAMod )
   
      ! Compute C_V using Eqn 1.47 or 1.48 depending on option selected
   C_V           = Get_C_V  ( Cn_alpha_q_circ, fprimeprime, p%UAMod )
   
      ! Compute Cn_v using either Eqn 1.45 or 1.49 depending on operating conditions
   Cn_v          = Get_Cn_v ( ds, T_V, xd%Cn_v_minus1(i,j), C_V, xd%C_V_minus1(i,j), xd%tau_V(i,j), T_VL, T_V0, Kalpha, u%alpha, alpha0 ) ! do I pass VRTX flag or tau_V?
   
   if (OtherState%FirstPass(i,j)) then
      Cn_v = 0.0_ReKi
   end if
   
      ! Finally, compute Cn using Eqn 1.50
   !Cn            = Cn_FS + Cn_v

 end subroutine ComputeKelvinChain
!==============================================================================   
                          

!==============================================================================
! Framework Routines                                                          !
!==============================================================================                               
      

!==============================================================================
subroutine UA_SetParameters( dt, InitInp, p, ErrStat, ErrMsg )
! 
! Called by : UA_Init
! Calls  to : NONE
!..............................................................................
   
   real(DbKi),                             intent(inout)  :: dt          ! time step length (s)
   type(UA_InitInputType),       intent(inout)  :: InitInp     ! input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UA_ParameterType),       intent(inout)  :: p           ! parameters
   integer(IntKi),                         intent(  out)  :: ErrStat     ! error status of the operation
   character(*),                           intent(  out)  :: ErrMsg      ! error message if ErrStat /= ErrID_None

      ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""
   p%dt         = dt
   
   allocate(p%c(InitInp%nNodesPerBlade,InitInp%numBlades), stat = ErrStat)
   if (ErrStat > ErrID_None) then
      ! Set errmessage and return
      return
   end if

   p%c          = InitInp%c        
   p%numBlades  = InitInp%numBlades
   p%nNodesPerBlade  = InitInp%nNodesPerBlade
   p%UAMod      = InitInp%UAMod    
   p%a_s        = InitInp%a_s 
   p%Flookup    = InitInp%Flookup
   
   
   
end subroutine UA_SetParameters
!==============================================================================   

!==============================================================================
subroutine UA_InitStates( p, xd, OtherState, ErrStat, ErrMsg )  
! Called by : UA_Init
! Calls  to : NONE
!..............................................................................

   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Initial other/optimization states
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
                                                                         
   
   ErrMsg  = ""
   ErrStat = ErrID_None
   
      ! allocate all the state arrays
   allocate(xd%alpha_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%alpha_minus2(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%q_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%q_minus2(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%X1_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%X2_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%X3_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%X4_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Kprime_alpha_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Kprime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Kprimeprime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%K3prime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Dp_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Cn_pot_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%fprimeprime_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%fprimeprime_c_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Df_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Df_c_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Dalphaf_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%alphaf_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%fprime_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%fprime_c_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%tau_V(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%Cn_v_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(xd%C_V_minus1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(OtherState%sigma1(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(OtherState%sigma3(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(OtherState%TESF(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(OtherState%LESF(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(OtherState%VRTX(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
   allocate(OtherState%FirstPass(p%nNodesPerBlade,p%numBlades), stat = ErrStat)
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
   xd%X3_minus1            = 0.0_ReKi
   xd%X4_minus1            = 0.0_ReKi
   xd%Kprime_alpha_minus1  = 0.0_ReKi
   xd%Kprime_q_minus1      = 0.0_ReKi
   xd%Dp_minus1            = 0.0_ReKi
   xd%Cn_pot_minus1        = 0.0_ReKi
   xd%K3prime_q_minus1     = 0.0_ReKi
   xd%Kprimeprime_q_minus1 = 0.0_ReKi  
   xd%fprimeprime_minus1   = 0.0_ReKi
   xd%fprimeprime_c_minus1 = 0.0_ReKi
   xd%Df_minus1            = 0.0_ReKi
   xd%Df_c_minus1          = 0.0_ReKi
   xd%Dalphaf_minus1       = 0.0_ReKi
   xd%alphaf_minus1        = 0.0_ReKi
   xd%fprime_minus1        = 0.0_ReKi
   xd%fprime_c_minus1      = 0.0_ReKi
   xd%tau_V                = 0.0_ReKi 
   xd%Cn_v_minus1          = 0.0_ReKi
   xd%C_V_minus1           = 0.0_ReKi  ! This probably should not be set to 0.0, but should be set 

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
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(ErrMsg))                       :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                               :: errStat2    ! temporary Error status of the operation
   integer(IntKi)                               :: i,j, iNode, iOffset
   character(64)                                :: chanPrefix
   
      ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
   !call DispNVD( UA_Ver )
   
      ! Allocate and set parameter data structure using initialization data
   call UA_SetParameters( interval, InitInp, p, ErrStat, ErrMsg )
   if (ErrStat >= AbortErrLev) return
   
      ! initialize the discrete states
   call UA_InitStates( p, xd, OtherState, ErrStat, ErrMsg )     ! initialize the continuous states
   if (ErrStat >= AbortErrLev) return    
   
      ! Allocate and set the InitOut data
   p%NumOuts = 30
   p%Delim   = ''
   p%OutFmt  = 'ES15.4e2'
   p%OutSFmt = 'A15'
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
         InitOut%WriteOutputHdr(iOffset+7)  =trim(chanPrefix)//'Cn_alf_circ'
         InitOut%WriteOutputHdr(iOffset+8)  =trim(chanPrefix)//'Cn_alf_nc'
         InitOut%WriteOutputHdr(iOffset+9)  =trim(chanPrefix)//'Cn_q_circ'
         InitOut%WriteOutputHdr(iOffset+10) =trim(chanPrefix)//'Cn_q_nc'
         InitOut%WriteOutputHdr(iOffset+11) =trim(chanPrefix)//'Cm_alf_circ'
         InitOut%WriteOutputHdr(iOffset+12) =trim(chanPrefix)//'Cm_alf_nc'
         InitOut%WriteOutputHdr(iOffset+13) =trim(chanPrefix)//'Cm_q_circ'
         InitOut%WriteOutputHdr(iOffset+14) =trim(chanPrefix)//'Cm_q_nc'
         InitOut%WriteOutputHdr(iOffset+15) =trim(chanPrefix)//'Alpha_e'
         InitOut%WriteOutputHdr(iOffset+16) =trim(chanPrefix)//"f'"
         InitOut%WriteOutputHdr(iOffset+17) =trim(chanPrefix)//"f'_c"
         InitOut%WriteOutputHdr(iOffset+18) =trim(chanPrefix)//"f''"
         InitOut%WriteOutputHdr(iOffset+19) =trim(chanPrefix)//"f''_c"
         InitOut%WriteOutputHdr(iOffset+20) =trim(chanPrefix)//"f''_m"
         InitOut%WriteOutputHdr(iOffset+21) =trim(chanPrefix)//'T_f'
         InitOut%WriteOutputHdr(iOffset+22) =trim(chanPrefix)//'T_V'
         InitOut%WriteOutputHdr(iOffset+23) =trim(chanPrefix)//'LESF'
         InitOut%WriteOutputHdr(iOffset+24) =trim(chanPrefix)//'TESF'
         InitOut%WriteOutputHdr(iOffset+25) =trim(chanPrefix)//'VRTX'
         InitOut%WriteOutputHdr(iOffset+26) =trim(chanPrefix)//'tau_v'
         InitOut%WriteOutputHdr(iOffset+27) = trim(chanPrefix)//'Cn_v'         
         InitOut%WriteOutputHdr(iOffset+28) = trim(chanPrefix)//'C_V '         
         InitOut%WriteOutputHdr(iOffset+29) = trim(chanPrefix)//'Cn_FS' 
         InitOut%WriteOutputHdr(iOffset+30) = trim(chanPrefix)//'Cm_FS'
         InitOut%WriteOutputUnt(iOffset+1)   ='(deg)  '
         InitOut%WriteOutputUnt(iOffset+2) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+3) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+4) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+5) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+6) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+7) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+8) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+9) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+10) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+11) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+12) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+13) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+14) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+15) ='deg'
         InitOut%WriteOutputUnt(iOffset+16) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+17) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+18) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+19) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+20) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+21) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+22) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+23) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+24) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+25) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+26) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+27) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+28) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+29) ='(-)   '
         InitOut%WriteOutputUnt(iOffset+30) ='(-)   '
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
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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
   real(ReKi)                                   :: X3
   real(ReKi)                                   :: X4
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
   real(ReKi)                                   :: fprimeprime_c
   real(ReKi)                                   :: fprimeprime_m
   real(ReKi)                                   :: Df
   real(ReKi)                                   :: Df_c
   real(ReKi)                                   :: Dalphaf
   real(ReKi)                                   :: fprime
   real(ReKi)                                   :: fprime_c
   real(ReKi)                                   :: Cn_v              ! normal force coefficient due to the presence of LE vortex
   real(ReKi)                                   :: C_V               ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi)                                   :: Cn_FS
   real(ReKi)                                   :: alpha_e
   real(ReKi)                                   :: alpha0            ! zero lift angle of attack (radians)
   real(ReKi)                                   :: dalpha0
   real(ReKi)                                   :: alpha_f
   real(ReKi)                                   :: eta_e
   real(ReKi)                                   :: Kafactor
   real(ReKi)                                   :: Cn_q_circ
   real(ReKi)                                   :: Cn_q_nc
   real(ReKi)                                   :: St_sh
   real(ReKi)                                   :: T_sh, T_f, T_V
            
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               

      call ComputeKelvinChain(i, j, u, p, xd, OtherState, AFInfo, Cn_prime, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                                     St_sh, Kalpha, alpha_e, alpha0, dalpha0, alpha_f, eta_e, Kq, q_cur, X1, X2, X3, X4, &
                                     Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, Df_c,Dalphaf, fprime, fprime_c, fprimeprime_c, fprimeprime_m, &
                                     Cn_alpha_q_circ, Cn_q_circ, Cn_q_nc, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, T_f, T_V, ErrStat, ErrMsg )
      
      
      T_sh = 2.0_ReKi*(1.0_ReKi-fprimeprime) / St_sh
      
      !---------------------------------------------------------
      ! Update the OtherStates
      !---------------------------------------------------------
   
      if ( (Cn_prime > Cn1) .or. ( Cn_prime < Cn2 ) ) then  ! assumption is that Cn2 <  0.0 and Cn1 > 0
         OtherState%LESF(i,j) = .true.  ! LE separation can occur
      else
         OtherState%LESF(i,j) = .false.
      end if
                                 
      if ( fprimeprime < xd%fprimeprime_minus1(i,j)) then
         OtherState%TESF(i,j) = .true.  ! Separation point is moving towards the Leading Edge
      else
         OtherState%TESF(i,j) = .false. ! Separation point is moving towards the Trailing Edge
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
     ! if ( fprimeprime < 0.7 .AND.  fprimeprime > 0.3 ) then 
         if ( OtherState%TESF(i,j) ) then  ! Separating flow
            if (Kafactor < 0.0_ReKi) then
               OtherState%sigma1(i,j) = 2.0_ReKi  ! This must be the first check
            else if (.not. OtherState%LESF(i,j) ) then
               OtherState%sigma1(i,j) = 1.0_ReKi !.4_ReKi    ! Leading edge separation has not occurred
            else if (xd%fprimeprime_minus1(i,j) <= 0.7_ReKi) then ! For this else, LESF = True
               OtherState%sigma1(i,j) = 2.0_ReKi !1.0_ReKi 
            else
               OtherState%sigma1(i,j) = 1.75_ReKi
            end if
         else ! Reattaching flow
            if (.not. OtherState%LESF(i,j) ) then
               OtherState%sigma1(i,j) = 0.5_ReKi
            end if
         
            if ( OtherState%VRTX(i,j) .and. (xd%tau_V(i,j) <= T_VL) ) then  ! Still shedding a vortex?
               OtherState%sigma1(i,j) = 0.25_ReKi
            end if
            if (Kafactor > 0.0_ReKi) then
               OtherState%sigma1(i,j) = 0.75_ReKi
            end if
         end if
      
      OtherState%sigma3(i,j) = 1.0_ReKi
      
      !if ( (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) >= T_VL) ) then
      !   OtherState%sigma3(i,j) = 3.0_ReKi
      !   if (.not. OtherState%TESF(i,j)) then
      !      OtherState%sigma3(i,j) = 4.0_ReKi
      !      if ( OtherState%VRTX(i,j) .and. (xd%tau_V(i,j) <= T_VL) ) then
      !         if (Kafactor < 0.0_ReKi) then
      !            OtherState%sigma3(i,j) = 2.0_ReKi
      !         else
      !            OtherState%sigma3(i,j) = 1.0_ReKi
      !         end if
      !      end if
      !   end if
      !else if (Kafactor < 0 ) then 
      !   OtherState%sigma3(i,j) = 4.0_ReKi
      !end if
      
         ! We are testing this heirarchical logic instead of the above block 5/29/2015
      if ( (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) >= T_VL) ) then
         OtherState%sigma3(i,j) =  3.0_ReKi
      else if (.not. OtherState%TESF(i,j)) then
         OtherState%sigma3(i,j) =  4.0_ReKi
      else if ( OtherState%VRTX(i,j) .and. (xd%tau_V(i,j) <= T_VL) ) then
         if (Kafactor < 0.0_ReKi) then
            OtherState%sigma3(i,j) =  2.0_ReKi
         else
            OtherState%sigma3(i,j) = 1.0_ReKi
          end if           
      else if (Kafactor < 0 ) then 
         OtherState%sigma3(i,j) =  4.0_ReKi
      end if
      
      
      if ((.not. OtherState%TESF(i,j)) .and. (Kq*dalpha0 < 0.0_ReKi)) then
         OtherState%sigma3(i,j) = 1.0_ReKi
      end if
   
      
     
      !---------------------------------------------------------
      ! Update the Discrete States, xd
      !---------------------------------------------------------
      if (OtherState%FirstPass(i,j)) then
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
      xd%X3_minus1(i,j)           = X3
      xd%X4_minus1(i,j)           = X4
      xd%Kprime_alpha_minus1(i,j) = Kprime_alpha
      xd%Kprime_q_minus1(i,j)     = Kprime_q
      xd%Dp_minus1(i,j)           = Dp
      xd%Cn_pot_minus1(i,j)       = Cn_pot
      xd%fprimeprime_minus1(i,j)  = fprimeprime
      xd%Df_minus1(i,j)           = Df
      if (p%Flookup) then
         xd%Df_c_minus1(i,j)           = Df_c
         xd%fprimeprime_c_minus1(i,j)  = fprimeprime_c
         xd%fprime_c_minus1(i,j)       = fprime_c
      end if
      
      xd%Dalphaf_minus1(i,j)      = Dalphaf
      xd%fprime_minus1(i,j)       = fprime
      xd%alphaf_minus1(i,j)       = alpha_f
      xd%Cn_v_minus1(i,j)         = Cn_v
      xd%C_V_minus1(i,j)          = C_V
      OtherState%FirstPass(i,j)        = .false.
   
      if ( xd%tau_V(i,j) > 0.0 .or. OtherState%LESF(i,j) ) then                !! TODO:  Verify this condition 2/20/2015 GJH
         xd%tau_V(i,j)          = xd%tau_V(i,j) + 2.0_ReKi*p%dt*u%U / p%c(i,j)  
      end if
   
end subroutine UA_UpdateDiscState
!==============================================================================   

!============================================================================== 
subroutine UA_UpdateStates( i, j, u, p, xd, OtherState, AFInfo, ErrStat, ErrMsg )
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
   integer(IntKi),                intent(  out) :: ErrStat         ! Error status of the operation
   character(*),                  intent(  out) :: ErrMsg          ! Error message if ErrStat /= ErrID_None

      ! Local variables  
   
   
   integer(IntKi)                               :: errStat2        ! Error status of the operation (secondary error)
   character(len(ErrMsg))                       :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   
      ! Initialize variables

   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
         
   
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
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   integer(IntKi)                                         :: errStat2        ! Error status of the operation (secondary error)
   character(len(ErrMsg))                                 :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
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
   real(ReKi)                                             :: X3
   real(ReKi)                                             :: X4
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
   real(ReKi)                                             :: fprimeprime_c
   real(ReKi)                                             :: fprimeprime_m
   real(ReKi)                                             :: Df
   real(ReKi)                                             :: Df_c
   real(ReKi)                                             :: Dalphaf
   real(ReKi)                                             :: fprime
   real(ReKi)                                             :: fprime_c
   real(ReKi)                                             :: Cn_v                ! normal force coefficient due to the presence of LE vortex
   real(ReKi)                                             :: C_V                 ! contribution to the normal force coefficient due to accumulated vorticity in the LE vortex
   real(ReKi)                                             :: Cn_FS, Cm_FS  
   real(ReKi)                                             :: alpha_e
   real(ReKi)                                             :: alpha0              ! zero lift angle of attack (radians)
   real(ReKi)                                             :: dalpha0
   real(ReKi)                                             :: alpha_f
   real(ReKi)                                             :: eta_e
   real(ReKi)                                             :: St_sh, Kafactor, f_c_offset
   real(ReKi)                                             :: Cl_static, Cd_static, Cm_static, Cn_static, Cm_alpha_nc, Cn_q_circ, Cn_q_nc, T_f, T_V
   real(ReKi)                                             :: M, alpha1, alpha2, C_nalpha, C_nalpha_circ, T_f0, T_V0, T_p, b1,b2,b5,A1,A2,A5,S1,S2,S3,S4,k1_hat,f, k2_hat
   integer                                                :: iOffset
   
   !BJJ: what are OtherState%iBladeNode, OtherState%iBlade here? I don't see them set in the routine that calls UA_CalcOutput
   
   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
   
  
   iOffset = (OtherState%iBladeNode-1)*p%NumOuts + (OtherState%iBlade-1)*p%nNodesPerBlade*p%NumOuts 
   
   if (OtherState%FirstPass(OtherState%iBladeNode, OtherState%iBlade)) then
      
      
      call GetSteadyOutputs(AFInfo, u%alpha, y%Cl, y%Cd, y%Cm, Cd0, ErrStat, ErrMsg)
      
      y%Cn = y%Cl*cos(u%alpha) + (y%Cd-Cd0)*sin(u%alpha)
      y%Cc = y%Cl*sin(u%alpha) - (y%Cd-Cd0)*cos(u%alpha)
            
   else
      M           = u%U / p%a_s
      call AFI_GetAirfoilParams( AFInfo, M, u%Re, u%alpha, alpha0, alpha1, alpha2, eta_e, C_nalpha, C_nalpha_circ, &
                              T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, b5, A1, A2, A5, S1, S2, S3, S4, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, k1_hat, x_cp_bar, ErrMsg, ErrStat )           
   
      
      call ComputeKelvinChain( OtherState%iBladeNode, OtherState%iBlade, u, p, xd, OtherState, AFInfo, Cn_prime, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                                 St_sh, Kalpha, alpha_e, alpha0, dalpha0, alpha_f, eta_e, Kq, q_cur, X1, X2, X3, X4, &
                                 Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, Df_c, Dalphaf, fprime, fprime_c, fprimeprime_c, fprimeprime_m, &
                                 Cn_alpha_q_circ, Cn_q_circ, Cn_q_nc, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, T_f, T_V, ErrStat, ErrMsg )
      
            
         ! Eqn 1.50(a or b depending on UAMod)
      if (xd%tau_v(OtherState%iBladeNode, OtherState%iBlade) > 0.0) then
         y%Cn= Cn_FS + Cn_v
      else
         y%Cn= Cn_FS
      end if
      
         ! This offset is to account for a correction for negative values of f_c which we cannot return due to a squaring of the quantity when using a lookup.
      if (p%Flookup) then
         f_c_offset = 0.2_ReKi
      else
         f_c_offset = 0.0_ReKi
      end if
      
      if ( p%UAMod == 3 .OR. p%UAMod == 4) then
            ! Eqn 1.52a   TODO: NOTE:  This is what is used in AD v14, so we may need a way to use this when we want to match v14 as much as possible! GJH 5/21/2015
        ! TODO: Look at eliminating this vortex related contribution under reversing AOA and perhaps when VRTX = F
         Kafactor             = Kalpha*dalpha0
         if  ( Kafactor < 0.0_ReKi ) then ! .AND.  .NOT. OtherState%VRTX(OtherState%iBladeNode, OtherState%iBlade )  )then
            y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset)
         else         
            y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset) + Cn_v*tan(alpha_e)*(1-xd%tau_v(OtherState%iBladeNode, OtherState%iBlade))
         end if
         
      elseif ( p%UAMod == 2 ) then
            ! Eqn 1.52b
         y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset)
      else
         if ( Cn_prime <= Cn1 ) then
            y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset) !*sin(alpha_e + alpha0)
         else
            if ( p%flookup ) then 
              ! if (p%UAMod == 1) then
              !    f      = Get_f_from_Lookup( p%UAMod, u%Re, u%alpha, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
              ! else   
               ! TODO: Need to understand the effect of the offset on this f in the following equations and fprimeprime_c.
                  f      = Get_f_c_from_Lookup( u%Re, u%alpha, alpha0, C_nalpha, AFInfo, ErrStat, ErrMsg)
              ! endif
               
            else
               f      = Get_f( u%alpha, alpha0, alpha1, alpha2, S1, S2, S3, S4 )
            end if
            
            k2_hat = 2*(Cn_prime-Cn1) + fprimeprime_c - f
            y%Cc   = k1_hat + Cc_pot*sqrt(fprimeprime_c)*fprimeprime_c**k2_hat !*sin(alpha_e + alpha0)
            
         end if
         
      end if
      
         ! Eqn 1.2
      y%Cl = y%Cn*cos(u%alpha) + y%Cc*sin(u%alpha)
      y%Cd = y%Cn*sin(u%alpha) - y%Cc*cos(u%alpha) + Cd0
            
         ! Eqn 1.55
         ! Compute Cn_FS using Eqn 1.35 or 1.36 depending on option selected
      if ( p%UAMod == 2 ) then
         call GetSteadyOutputs(AFInfo, alpha_f, Cl_static, Cd_static, Cm_static, Cd0, ErrStat, ErrMsg)
         Cn_static = Cl_static*cos(alpha_f) + (Cd_static-Cd0)*sin(alpha_f)
         ! TODO: What about when Cn = 0  GJH 5/22/2015
         fprimeprime_m = (Cm_static - Cm0) / Cn_static 
      else
         fprimeprime_m = 0.0
      end if
      
      Cm_FS = Get_Cm_FS( Cm0, k0, k1, k2, k3, Cn_alpha_q_circ, fprimeprime, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Dalphaf, alpha_f, Cn_FS, fprimeprime_m, AFInfo, p%UAMod, Cm_alpha_nc, ErrStat, ErrMsg )
      
      y%Cm = Get_Cm( Cm_FS, T_VL, x_cp_bar, Cn_v, xd%tau_v(OtherState%iBladeNode, OtherState%iBlade) )
                  
   end if
   
   if (allocated(y%WriteOutput)) then  !bjj: because BEMT uses local variables for UA output, y%WriteOutput is not necessarially allocated. Need to figure out a better solution.
      y%WriteOutput(iOffset+1) = u%alpha*180.0/pi
      y%WriteOutput(iOffset+2) = y%Cn
      y%WriteOutput(iOffset+3) = y%Cc
      y%WriteOutput(iOffset+4) = y%Cm
      y%WriteOutput(iOffset+5) = y%Cl
      y%WriteOutput(iOffset+6) = y%Cd
      y%WriteOutput(iOffset+7) = C_nalpha_circ*alpha_e  ! = Cn_alpha_circ
      y%WriteOutput(iOffset+8) = Cn_alpha_nc
      y%WriteOutput(iOffset+9) = Cn_q_circ
      y%WriteOutput(iOffset+10) = Cn_q_nc
      y%WriteOutput(iOffset+11) = 0.0  ! = Cm_alpha_circ  NOTE: we do not use this quantity
      y%WriteOutput(iOffset+12) = Cm_alpha_nc
      y%WriteOutput(iOffset+13) = Cm_q_circ
      y%WriteOutput(iOffset+14) = Cm_q_nc
      y%WriteOutput(iOffset+15) = alpha_e*180/pi
      y%WriteOutput(iOffset+16) = fprime
      y%WriteOutput(iOffset+17) = fprime_c
      y%WriteOutput(iOffset+18) = fprimeprime
      y%WriteOutput(iOffset+19) = fprimeprime_c
      y%WriteOutput(iOffset+20) = fprimeprime_m
      y%WriteOutput(iOffset+21) = T_f
      y%WriteOutput(iOffset+22) = T_V
      if ( OtherState%LESF(OtherState%iBladeNode, OtherState%iBlade) ) then
         y%WriteOutput(iOffset+23) = 1.0_ReKi
      else
         y%WriteOutput(iOffset+23) = 0.0_ReKi
      end if
      
      if ( OtherState%TESF(OtherState%iBladeNode, OtherState%iBlade) ) then
         y%WriteOutput(iOffset+24) = 1.0_ReKi
      else
         y%WriteOutput(iOffset+24) = 0.0_ReKi
      end if
      if ( OtherState%VRTX(OtherState%iBladeNode, OtherState%iBlade) ) then
         y%WriteOutput(iOffset+25) = 1.0_ReKi
      else
         y%WriteOutput(iOffset+25) = 0.0_ReKi
      end if
      y%WriteOutput(iOffset+26) = xd%tau_v(OtherState%iBladeNode, OtherState%iBlade)
      y%WriteOutput(iOffset+27) = Cn_v         
      y%WriteOutput(iOffset+28) = C_V          
      y%WriteOutput(iOffset+29) = Cn_FS 
      y%WriteOutput(iOffset+30) = Cm_FS
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
 !  integer     :: UAMod
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
 !  UAMod         = 1
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
 !  !            alpha1, S1, S2, UAMod )
 !              !dt, ds, U, M, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
 !              !              q_minus1, q_minus2, T_alpha, T_q, T_p, T_f, X_1_minus1, X_2_minus1, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, Df_minus1, fprime_minus1, b1, b2, A1, A2, C_nalpha, &
 !              !              alpha1, S1, S2, UAMod )
 !  
 !  Test_Cn = 1
 !  
 !end function Test_Cn
 !  
   
   
   
   
end module UnsteadyAero