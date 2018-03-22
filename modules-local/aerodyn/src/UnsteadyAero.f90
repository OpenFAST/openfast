!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
module UnsteadyAero

   ! This module uses equations defined in the document "The Unsteady Aerodynamics Module for FAST 8" by Rick Damiani and Greg Hayman, 28-Feb-2017

use NWTC_Library   
   use UnsteadyAero_Types
   use AirfoilInfo
   
   implicit none 

private
   type(ProgDesc), parameter  :: UA_Ver = ProgDesc( 'UnsteadyAero', '', '' )

   public :: UA_Init
   public :: UA_UpdateDiscOtherState
   public :: UA_UpdateStates
   public :: UA_CalcOutput


   contains
   
! **************************************************
FUNCTION SAT( X, VAL, SLOPE )
 !  AOA saturation function 02/15/98
 ! **************************************************

IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi)                 :: SAT
REAL(ReKi),INTENT(IN)       :: SLOPE
REAL(ReKi),INTENT(IN)       :: VAL
REAL(ReKi),INTENT(IN)       :: X


IF ( ABS(X) <= VAL )  THEN
    SAT = X
ELSEIF ( X > VAL)  THEN
    SAT = SLOPE * X + VAL * ( 1. - SLOPE )
ELSE
    SAT = SLOPE * X - VAL * ( 1. - SLOPE )
ENDIF


RETURN
END FUNCTION SAT
   
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
   
   real(ReKi)                      :: IntAFCoefs(4)         ! The interpolated airfoil coefficients.
   integer                         :: s1                    ! Number of columns in the AFInfo structure
   real(ReKi)                      :: Alpha                 ! AOA in range [-pi,pi]

   
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   IntAFCoefs = 0.0_ReKi ! initialize in case we only don't have 4 columns in the airfoil data (i.e., so cm is zero if not in the file)
   
      
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   s1 = size(AFInfo%Table(1)%Coefs,2)
   !if (s1 < 3) then
   !   ErrMsg  = 'The Airfoil info table must contains columns for lift, drag, and pitching moment'
   !   ErrStat = ErrID_Fatal
   !   return
   !end if
   Cd0 =   AFInfo%Table(1)%UA_BL%Cd0
   
   Alpha = AOA
   call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
   IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha &
                                             , AFInfo%Table(1)%Alpha &
                                             , AFInfo%Table(1)%Coefs &
                                             , AFInfo%Table(1)%SplineCoefs &
                                             , ErrStat, ErrMsg )
   if (ErrStat >= AbortErrLev) return
      
   Cl = IntAFCoefs(1)
   Cd = IntAFCoefs(2)
   Cm = IntAFCoefs(3)
   
end subroutine GetSteadyOutputs
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
! TE Flow Separation Equations                                                !
!==============================================================================

!==============================================================================
real(ReKi) function Get_f_from_Lookup( UAMod, Re, alpha, alpha0, C_nalpha_circ, AFInfo, ErrStat, ErrMsg)
! Compute either fprime or fprimeprime using an analytical equation (and eventually a table lookup)
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   integer,          intent(in   ) :: UAMod
   real(ReKi),       intent(in   ) :: Re            ! Reynolds number
   real(ReKi),       intent(in   ) :: alpha         ! angle of attack (radians)
   real(ReKi),       intent(in   ) :: alpha0
   real(ReKi),       intent(in   ) :: C_nalpha_circ
   type(AFInfoType), intent(in   ) :: AFInfo        ! The airfoil parameter data
   integer(IntKi),   intent(  out) :: ErrStat               ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg                ! Error message if ErrStat /= ErrID_None
   
   !real                            :: IntAFCoefs(4)         ! The interpolated airfoil coefficients.
   real(ReKi)                       :: Cn, Cl, Cd, Cm, Cd0, tmpRoot, denom
   !integer                         :: s1                    ! Number of columns in the AFInfo structure
   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
    
   
   
   
   call GetSteadyOutputs(AFInfo, alpha, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
      if (ErrStat >= AbortErrLev ) return
   
   Cn =  Cl*cos(alpha) + (Cd-Cd0)*sin(alpha)
   denom = (C_nalpha_circ*(alpha-alpha0))
   
   
! #ifndef CHECK_DENOM_DIFF
   if (abs(denom) < .01) then
!#else
!   if (EqualRealNos(denom,0.0_ReKi)) then
!#endif   

       Get_f_from_Lookup = 1.0_ReKi
       return
   end if
   
   tmpRoot  = Cn/denom
   
   if (tmpRoot < 0.0_ReKi) then
      Get_f_from_Lookup = 1.0_ReKi
      !TODO: Should tmpRoot = 0.0 instead so that we can still solve the equations below instead of returning 2/28/2017 GJH
      return
   end if
   
   if (UAMod == 2) then
      Get_f_from_Lookup = ((3*sqrt(tmpRoot)-1)/2.0)**2
   else
      Get_f_from_Lookup = ( 2 * sqrt( tmpRoot ) - 1 ) **2 
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
   real(ReKi)                      :: Cc, Cl, Cd, Cm, Cd0, denom
   !integer                         :: s1                    ! Number of columns in the AFInfo structure
   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
  
   
   call GetSteadyOutputs(AFInfo, alpha, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
      if (ErrStat >= AbortErrLev) return
   
   denom = ( alpha-alpha0 )*tan(alpha)  !NOTE,NOTE: On 8/27/15 GJH Added back tan(alpha) because results for Fy did not match steady state without it. ! NOTE: We removed the tan(alpha) term from the equation and repace with another (alpha-alpha0) term, per Rick's suggestion 8/13/2015.    *tan(alpha)
   
   !denom = ( alpha-alpha0 )**2  !testing again On 9/16/15
   if (abs(denom) < .015 ) then
      Get_f_c_from_Lookup = 1.44_ReKi
      return
   end if
   denom = C_nalpha*denom
   Cc =  Cl*sin(alpha) - (Cd-Cd0)*cos(alpha)
      ! Apply an offset of 0.2 to fix cases where f_c should be negative, but we are using **2 so can only return positive values
      ! Note: because we include this offset, it must be accounted for in the final value of Cc, eqn 1.40.  This will be applied
      ! For both UA_Mod = 1,2, and 3 when using Flookup = T
   Get_f_c_from_Lookup = (  Cc / ( denom )  + 0.2 ) **2 
   
   
   !Get_f_c_from_Lookup = (  Cc / ( C_nalpha*( alpha-alpha0 )*tan(alpha) )  ) **2 
   !if ( Cc < 0.0_ReKi ) then
   !   Get_f_c_from_Lookup = -Get_f_c_from_Lookup
   !end if
   
   if ( Get_f_c_from_Lookup > 1.44 ) then
      Get_f_c_from_Lookup = 1.44_ReKi
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
subroutine ComputeKelvinChain( i, j, u, p, xd, OtherState, misc, AFInfo, Cn_prime, Cn_prime_diff, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                               St_sh, Kalpha_f, Kq_f, alpha_filt_cur, alpha_e, alpha0, dalpha0, alpha_f, eta_e, Kq, q_cur, q_f_cur, X1, X2, X3, X4, &
                               Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, Df_c, Dalphaf, fprime, fprime_c, fprimeprime_c, fprimeprime_m, &
                               Cn_alpha_q_nc, Cn_alpha_q_circ, Cn_q_circ, Cn_q_nc, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, T_f, T_V, ErrStat, ErrMsg )
! 
! Called by : DRIVER
! Calls  to : Get_Beta_M_Sqrd, Get_Beta_M, AFI_GetAirfoilParams, Get_ds, Get_Pitchrate, Get_k_, Get_Kupper, Get_ExpEqn
!             Get_Cn_nc, Get_alpha_e, Get_Cm_q_circ, Get_Cm_q_nc, Get_f, Get_Cn_FS, Get_C_V, Get_Cn_v
!...............................................................................

      
   integer(IntKi),                         intent(in   ) :: i,j               ! Blade node indices
   type(UA_InputType),                     intent(in   ) :: u                 ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   type(UA_ParameterType),                 intent(in   ) :: p                 ! Parameters   
   type(UA_DiscreteStateType),             intent(in   ) :: xd                ! Input: Discrete states at t;
   type(UA_OtherStateType),                intent(in   ) :: OtherState        ! Other states at t
   type(UA_MiscVarType),                   intent(inout) :: misc              ! Misc/optimization variables
   type(AFInfoType),                       intent(in   ) :: AFInfo            ! The airfoil parameter data
   real(ReKi),                             intent(  out) :: Cn_prime          !
   real(ReKi),                             intent(  out) :: Cn_prime_diff
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
   real(ReKi),                             intent(  out) :: Kalpha_f          ! filtered  backwards finite difference of alpha
   real(ReKi),                             intent(  out) :: Kq_f              ! filtered  backwards finite difference of q
   real(ReKi),                             intent(  out) :: alpha_filt_cur    ! filtered angle of attack                          
   real(ReKi),                             intent(  out) :: alpha_e           ! effective angle of attack at 3/4 chord  TODO: verify 3/4 and not 1/4                          
   real(ReKi),                             intent(  out) :: alpha0            ! zero lift angle of attack (radians)
   real(ReKi),                             intent(  out) :: dalpha0           !
   real(ReKi),                             intent(  out) :: alpha_f           !
   real(ReKi),                             intent(  out) :: eta_e             !
   real(ReKi),                             intent(  out) :: Kq                !
   real(ReKi),                             intent(  out) :: q_cur             !
   real(ReKi),                             intent(  out) :: q_f_cur             !
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
   real(ReKi),                             intent(  out) :: Cn_alpha_q_nc     ! non-circulatory component of normal force coefficient response to step change in alpha and q
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
   real(ReKi)                :: Kalpha                                        !
   real(ReKi)                :: Kalpha_f_minus1                               !
   real(ReKi)                :: Kq_f_minus1                                   !
   real(ReKi)                :: alpha_minus1                                  !
   real(ReKi)                :: alpha_filt_minus1                                  !
   real(ReKi)                :: q_minus1                                      !
   real(ReKi)                :: q_f_minus1
   real(ReKi)                :: fprime_minus1                                 !
   real(ReKi)                :: Cn_pot_minus1                                 !
   real(ReKi)                :: k1_hat                                        !
   real(ReKi)                :: K3prime_q                                     !
   real(ReKi)                :: k_mq                                          !
   real(ReKi)                :: Kprimeprime_q                                 !
   real(ReKi)                :: fprime_c_minus1
   real(ReKi)                :: alphaf_minus1
   real(ReKi)                :: LowPassConst
   real(ReKi)                :: filtCutOff
   real(ReKi)                :: factor

   integer(IntKi)            :: ErrStat2
   character(ErrMsgLen)      :: ErrMsg2
   character(*), parameter   :: RoutineName = 'ComputeKelvinChain'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   fprimeprime_m = 0

   M           = u%U / p%a_s
   call UA_CheckMachNumber(M, misc%FirstWarn_M, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
      
   beta_M_Sqrd = 1 - M**2
   beta_M      = sqrt(1 - M**2) 
   T_I         = p%c(i,j) / p%a_s                                       ! Eqn 1.11c
   ds          = 2.0_ReKi*u%U*p%dt/p%c(i,j)                             ! Eqn 1.5b
   
   ! Lookup values using Airfoil Info module
   call AFI_GetAirfoilParams( AFInfo, M, u%Re, alpha0, alpha1, alpha2, eta_e, C_nalpha, C_nalpha_circ, &
                              T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, b5, A1, A2, A5, S1, S2, S3, S4, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, k1_hat, x_cp_bar, filtCutOff, ErrMsg2, ErrStat2 )           
      ! AFI_GetAirfoilParams doesn't return error, so I will not check
   
      ! Override eta_e if we are using Flookup
   if ( p%Flookup ) then
      eta_e = 1.0
   end if
   
   if (OtherState%FirstPass(i,j)) then
      alpha_minus1 = u%alpha      
      alpha_filt_minus1 = u%alpha     
   else
      alpha_minus1 = xd%alpha_minus1(i,j)     
      alpha_filt_minus1 = xd%alpha_filt_minus1(i,j)   
   end if
   
   ! Low Pass filtering of alpha, q, and Kq in order to deal with numerical issues that arise for very small values of dt
   ! See equations 1.7 and 1.8 of the manual
   
   LowPassConst  =  exp(-2.0_ReKi*PI*p%dt*filtCutOff)
   
   alpha_filt_cur = LowPassConst*alpha_filt_minus1 + (1.0_ReKi-LowPassConst)*u%alpha
   
   
   dalpha0  = alpha_filt_cur - alpha0
   
    
      ! Compute Kalpha using Eqn 1.7
  
   Kalpha        = ( alpha_filt_cur - alpha_filt_minus1 ) / p%dt
   
   if (OtherState%FirstPass(i,j)) then
      Kalpha_f_minus1 = 0.0_ReKi
   else
      Kalpha_f_minus1 = xd%Kalpha_f_minus1(i,j)
   end if
   
   q_cur         = Kalpha   * p%c(i,j) / u%U   ! Kalpha here uses the low-pass filtered alphas
   
       
   if (OtherState%FirstPass(i,j)) then
      q_minus1   = q_cur 
      q_f_minus1 = q_cur    
   else    
      q_minus1   = xd%q_minus1(i,j) 
      q_f_minus1 = xd%q_f_minus1(i,j) 
   end if  
   
   q_f_cur       = LowPassConst*q_f_minus1 + (1.0_ReKi-LowPassConst)*q_cur
      
      
#ifdef TEST_THEORY
      ! Change found in ADv14
   q_f_cur = SAT( q_f_cur, 0.03_ReKi, 0.1_ReKi )
#endif   
   
   Kalpha_f      = q_f_cur * u%U / p%c(i,j)  ! Kalpha_f is using the low-pass filtered q

   k_alpha  = 1.0_ReKi / ( (1.0_ReKi - M) + (C_nalpha/2.0_ReKi) * M**2 * beta_M * (A1*b1 + A2*b2) )                 ! Eqn 1.11a
   k_q      = 1.0_ReKi / ( (1.0_ReKi - M) + C_nalpha            * M**2 * beta_M * (A1*b1 + A2*b2) )                 ! Eqn 1.11b
   T_alpha  = T_I * k_alpha * 0.75                                                 ! Eqn 1.10a
   T_q      = T_I * k_q * 0.75                                                     ! Eqn 1.10b
      
      ! These quantities are needed for the update state calculations, but are then updated themselves based on the logic which follows
   T_f           = T_f0 / OtherState%sigma1(i,j)       ! Eqn 1.37
   T_V           = T_V0 / OtherState%sigma3(i,j)       ! Eqn 1.48
   
    
   ! Compute Kq  using Eqn 1.9  with time-shifted q s
#ifdef TEST_THEORY
   Kq            = ( q_f_cur  - q_f_minus1 ) / p%dt
#else
   
   Kq            = ( q_cur  - q_minus1 ) / p%dt
#endif

   if (OtherState%FirstPass(i,j)) then
      Kq_f_minus1 = 0.0_ReKi
   else
      Kq_f_minus1 = xd%Kq_f_minus1(i,j)
   end if
   
   Kq_f          =  LowPassConst*Kq_f_minus1 + (1.0_ReKi-LowPassConst)*Kq
   
  
      ! Compute Kprime_alpha using Eqn 1.18b
   Kprime_alpha  = Get_ExpEqn( real(p%dt,ReKi), T_alpha, xd%Kprime_alpha_minus1(i,j), Kalpha_f, Kalpha_f_minus1 )
   
      ! Compute Kprime_q using Eqn 1.19b    
   Kprime_q      = Get_ExpEqn( real(p%dt,ReKi), T_q    , xd%Kprime_q_minus1(i,j)    ,  Kq_f   , Kq_f_minus1     )
   
      ! Compute Cn_alpha_nc using Eqn 1.18a 
   Cn_alpha_nc   = 4.0_ReKi*T_alpha * ( Kalpha_f - Kprime_alpha ) / M
   
      ! Compute Cn_alpha_nc using Eqn 1.19a  
   Cn_q_nc       = -1.0_ReKi*T_q * ( Kq_f - Kprime_q ) / M
   
      ! Compute Cn_alpha_q_nc using Eqn 1.17
   Cn_alpha_q_nc = Cn_alpha_nc + Cn_q_nc  
   
      ! Compute X1 using Eqn 1.15a   
   X1            = Get_ExpEqn( ds*beta_M_Sqrd*b1, 1.0_ReKi, xd%X1_minus1(i,j), A1*(alpha_filt_cur - alpha_filt_minus1), 0.0_ReKi )
   
      ! Compute X2 using Eqn 1.15b
   X2            = Get_ExpEqn( ds*beta_M_Sqrd*b2, 1.0_ReKi, xd%X2_minus1(i,j), A2*(alpha_filt_cur - alpha_filt_minus1), 0.0_ReKi )
   
      ! Compute alpha_e using Eqn 1.14
   alpha_e       = (alpha_filt_cur - alpha0) - X1 - X2  
   
      ! Compute Cn_alpha_q_circ using Eqn 1.13  
   Cn_alpha_q_circ = C_nalpha_circ * alpha_e
   
   if ( p%UAMod == 2 ) then
         ! Compute X3 and X4 using Eqn 1.16a  and then add Cn_q_circ (Eqn 1.16) to the previously computed Cn_alpha_q_circ
      X3              = Get_ExpEqn( ds*beta_M_Sqrd*b1, 1.0_ReKi, xd%X3_minus1(i,j), A1*(q_f_cur - q_f_minus1), 0.0_ReKi )
      X4              = Get_ExpEqn( ds*beta_M_Sqrd*b2, 1.0_ReKi, xd%X4_minus1(i,j), A2*(q_f_cur - q_f_minus1), 0.0_ReKi )
      Cn_q_circ       = C_nalpha_circ*q_f_cur/2.0 - X3 - X4  
   else
      Cn_q_circ       = 0.0
   end if
   
      ! Compute K3prime_q using Eqn 1.26
   K3prime_q       = Get_ExpEqn( b5*beta_M_Sqrd*ds, 1.0_ReKi, xd%K3prime_q_minus1(i,j),  A5*(q_f_cur - q_f_minus1), 0.0_ReKi )
   
      ! Compute Cm_q_circ using Eqn 1.25
   Cm_q_circ       = -C_nalpha*(q_f_cur -K3prime_q)*p%c(i,j)/(16.0*beta_M*u%U)
   
      ! Compute Cn_pot using eqn 1.20a
   Cn_pot          = Cn_alpha_q_circ + Cn_alpha_q_nc
   
      ! Eqn 1.29b
   k_mq            = 7.0/(15.0*(1.0-M)+1.5*C_nalpha*A5*b5*beta_M*M**2)
   
      ! Eqn 1.29c
   Kprimeprime_q   = Get_ExpEqn( real(p%dt,ReKi), k_mq**2*T_I   , xd%Kprimeprime_q_minus1(i,j)    ,  Kq_f   , Kq_f_minus1     )
   
      ! Compute Cm_q_nc 
   if ( p%UAMod == 3 ) then
      ! Implements eqn 1.31
      Cm_q_nc =  -Cn_q_nc/4.0 - (k_alpha**2)*T_I*(Kq_f-Kprimeprime_q)/(3.0*M)
   else  
         ! Implements eqn 1.29a
      Cm_q_nc = -7*(k_mq**2)*T_I*(Kq_f-Kprimeprime_q)/(12.0*M)
      
   end if
   
      ! Compute Cc_pot using eqn 1.21
   Cc_pot          = Cn_alpha_q_circ*tan(alpha_e+alpha0)
   
   if (OtherState%FirstPass(i,j)) then
      Cn_pot_minus1 = Cn_pot
   else
      Cn_pot_minus1 = xd%Cn_pot_minus1(i,j)
   end if
   
      ! Compute Dp using Eqn 1.35b
   Dp            = Get_ExpEqn( ds, T_p, xd%Dp_minus1(i,j), Cn_pot, Cn_pot_minus1 )
   
      ! Compute Cn_prime using Eqn 1.35a
   Cn_prime      = Cn_Pot - Dp
   
   if (OtherState%FirstPass(i,j)) then
      Cn_prime_diff = 0.0_ReKi
   else
      Cn_prime_diff = Cn_prime - xd%Cn_prime_minus1(i,j)
   end if
   

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! This code is taken from ADv14 but doesn't reflect the original intent of the UA theory document
#ifdef TEST_THEORY
   IF ( alpha_filt_cur * Cn_prime_diff < 0. ) THEN
      T_f   = T_f0*1.5
   ELSE
      T_f   = T_f0
   ENDIF
#endif   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
      ! Compute alpha_f using Eqn 1.34
   alpha_f       = Cn_prime / C_nalpha_circ + alpha0
   
      ! Compute fprime using Eqn 1.32 and Eqn 1.33
   if (p%flookup) then
      fprime        = Get_f_from_Lookup( p%UAMod, u%Re, alpha_f, alpha0, C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   else   
      fprime        = Get_f( alpha_f, alpha0, alpha1, alpha2, S1, S2, S3, S4)
   end if
   
   if (OtherState%FirstPass(i,j)) then
      fprime_minus1 = fprime
   else
      fprime_minus1 = xd%fprime_minus1(i,j)
   end if
   
      ! Compute Df using Eqn 1.36b   
   Df            = Get_ExpEqn( ds, T_f, xd%Df_minus1(i,j), fprime, fprime_minus1 )
   
      ! Compute fprimeprime using Eqn 1.36a
   fprimeprime   = fprime - Df
   
   if (p%Flookup) then
         ! Compute fprime using Eqn 1.32 and Eqn 1.33
      fprime_c   = Get_f_c_from_Lookup( u%Re, alpha_f, alpha0, C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
     
      if (OtherState%FirstPass(i,j)) then
         fprime_c_minus1 = fprime_c       
      else
         fprime_c_minus1 = xd%fprime_c_minus1(i,j)        
      end if
   
         ! Compute Df using Eqn 1.36b   
      Df_c            = Get_ExpEqn( ds, T_f, xd%Df_c_minus1(i,j), fprime_c, fprime_c_minus1 )
   
         ! Compute fprimeprime using Eqn 1.36a
      fprimeprime_c   = fprime_c - Df_c
   else
      fprimeprime_c   = fprimeprime
   end if
   
    
   if ( p%UAMod == 2 ) then
         ! use equation 1.39
      Cn_FS   = Cn_alpha_q_nc + Cn_alpha_q_circ *  ( (1.0_ReKi + 2.0_ReKi*sqrt(fprimeprime) ) / 3.0_ReKi  )**2 + Cn_q_circ
   else
         ! use equation 1.38
      Cn_FS   = Cn_alpha_q_nc + Cn_alpha_q_circ *  ( (1.0_ReKi + sqrt(fprimeprime)          ) / 2.0_ReKi )**2
   
   end if
   
      
   if ( p%UAMod == 3 ) then
      if (OtherState%FirstPass(i,j)) then     
         alphaf_minus1   = alpha_f
      else
         alphaf_minus1   = xd%alphaf_minus1(i,j)
      end if
         ! Eqn 1.43b
      Dalphaf    = Get_ExpEqn( ds, 0.1_ReKi*T_f, xd%Dalphaf_minus1(i,j), alpha_f, alphaf_minus1 )
   else
      Dalphaf    = 0.0_ReKi
   end if
   
      ! Compute C_V using Eqn 1.49 or 1.50 depending on option selected
   if ( p%UAMod == 2 ) then
      ! use equation 1.50
      C_V   = Cn_alpha_q_circ * ( 1.0_ReKi - ((1.0_ReKi + 2.0_ReKi*sqrt(fprimeprime) )/3.0_ReKi)**2  )
   else
      ! use equation 1.49
      C_V   = Cn_alpha_q_circ *  ( 1.0_ReKi - ( 0.5_ReKi + 0.5_ReKi*sqrt(fprimeprime) )**2 )
   end if
   
      ! Compute Cn_v using either Eqn 1.47 or 1.52 depending on operating conditions

   factor = (alpha_filt_cur - alpha0) * Kalpha_f
   
   if (xd%tau_V(i,j) > T_VL .AND. (factor > 0)) then 
         ! The assertion is the T_V will always equal T_V0/2 when this condition is satisfied
      Cn_v = xd%Cn_v_minus1(i,j)*exp(-ds/T_V)   ! Eqn 1.52    
   else      
      Cn_v = Get_ExpEqn( ds, T_V, xd%Cn_v_minus1(i,j), C_V, xd%C_V_minus1(i,j) )   ! Eqn 1.47
   end if
   
   if ( Cn_v < 0.0_ReKi ) then
      Cn_v = 0.0_ReKi
   end if
   
   if (OtherState%FirstPass(i,j)) then
      Cn_v = 0.0_ReKi
   end if

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

   integer(IntKi)            :: ErrStat2
   character(*), parameter   :: RoutineName = 'UA_SetParameters'
   
   
   
      ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""
   p%dt         = dt
   
   allocate(p%c(InitInp%nNodesPerBlade,InitInp%numBlades), stat = ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating p%c.', ErrStat, ErrMsg, RoutineName )
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
subroutine UA_InitStates_Misc( p, xd, OtherState, m, ErrStat, ErrMsg )  
! Called by : UA_Init
! Calls  to : NONE
!..............................................................................

   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Initial other states
   type(UA_MiscVarType),         intent(inout)  :: m           ! Initial misc/optimization variables
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer(IntKi)            :: ErrStat2
   character(ErrMsgLen)      :: ErrMsg2
   character(*), parameter   :: RoutineName = 'UA_InitStates'
         
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   
      ! allocate all the state arrays
   call AllocAry( xd%alpha_minus1,        p%nNodesPerBlade,p%numBlades, 'xd%alpha_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%alpha_filt_minus1,      p%nNodesPerBlade,p%numBlades, 'xd%alpha_filt_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%q_minus1,            p%nNodesPerBlade,p%numBlades, 'xd%q_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%q_f_minus1,          p%nNodesPerBlade,p%numBlades, 'xd%q_f_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%Kq_f_minus1,         p%nNodesPerBlade,p%numBlades, 'xd%Kq_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%Kalpha_f_minus1,     p%nNodesPerBlade,p%numBlades, 'xd%Kalpha_f_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%X1_minus1,           p%nNodesPerBlade,p%numBlades, 'xd%X1_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%X2_minus1,           p%nNodesPerBlade,p%numBlades, 'xd%X2_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%X3_minus1,           p%nNodesPerBlade,p%numBlades, 'xd%X3_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%X4_minus1,           p%nNodesPerBlade,p%numBlades, 'xd%X4_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%Kprime_alpha_minus1, p%nNodesPerBlade,p%numBlades, 'xd%Kprime_alpha_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%Kprime_q_minus1,     p%nNodesPerBlade,p%numBlades, 'xd%Kprime_q_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%Kprimeprime_q_minus1,p%nNodesPerBlade,p%numBlades, 'xd%Kprimeprime_q_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%K3prime_q_minus1,    p%nNodesPerBlade,p%numBlades, 'xd%K3prime_q_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry( xd%Dp_minus1,           p%nNodesPerBlade,p%numBlades, 'xd%Dp_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   call AllocAry(xd%Cn_prime_minus1     ,p%nNodesPerBlade,p%numBlades,'xd%Cn_prime_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%Cn_pot_minus1       ,p%nNodesPerBlade,p%numBlades,'xd%Cn_pot_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%fprimeprime_minus1  ,p%nNodesPerBlade,p%numBlades,'xd%fprimeprime_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%fprimeprime_c_minus1,p%nNodesPerBlade,p%numBlades,'xd%fprimeprime_c_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%Df_minus1           ,p%nNodesPerBlade,p%numBlades,'xd%Df_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%Df_c_minus1         ,p%nNodesPerBlade,p%numBlades,'xd%Df_c_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%Dalphaf_minus1      ,p%nNodesPerBlade,p%numBlades,'xd%Dalphaf_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%alphaf_minus1       ,p%nNodesPerBlade,p%numBlades,'xd%alphaf_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%fprime_minus1       ,p%nNodesPerBlade,p%numBlades,'xd%fprime_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%fprime_c_minus1     ,p%nNodesPerBlade,p%numBlades,'xd%fprime_c_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%tau_V               ,p%nNodesPerBlade,p%numBlades,'xd%tau_V',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%tau_V_minus1        ,p%nNodesPerBlade,p%numBlades,'xd%tau_V_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%Cn_v_minus1         ,p%nNodesPerBlade,p%numBlades,'xd%Cn_v_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(xd%C_V_minus1          ,p%nNodesPerBlade,p%numBlades,'xd%C_V_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   call AllocAry(OtherState%FirstPass,p%nNodesPerBlade,p%numBlades,'OtherState%FirstPass',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(OtherState%sigma1   ,p%nNodesPerBlade,p%numBlades,'OtherState%sigma1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(OtherState%sigma3   ,p%nNodesPerBlade,p%numBlades,'OtherState%sigma3',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

#ifdef UA_OUTS
   call AllocAry(m%TESF     ,p%nNodesPerBlade,p%numBlades,'m%TESF',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(m%LESF     ,p%nNodesPerBlade,p%numBlades,'m%LESF',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(m%VRTX     ,p%nNodesPerBlade,p%numBlades,'m%VRTX',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
#endif  
   
   if (ErrStat >= AbortErrLev) return
   
   m%FirstWarn_M = .true.   
   
   OtherState%sigma1    = 1.0_ReKi
   OtherState%sigma3    = 1.0_ReKi
   
#ifdef UA_OUTS
   m%TESF      = .FALSE.  
   m%LESF      = .FALSE.   
   m%VRTX      = .FALSE. 
#endif   
   
   OtherState%FirstPass = .true.
   
   
   xd%Cn_prime_minus1      = 0.0_ReKi
   xd%alpha_minus1         = 0.0_ReKi
   xd%alpha_filt_minus1    = 0.0_ReKi
   xd%q_minus1             = 0.0_ReKi
   xd%q_f_minus1           = 0.0_ReKi
   xd%Kq_f_minus1          = 0.0_ReKi
   xd%Kalpha_f_minus1      = 0.0_ReKi
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
   xd%tau_V_minus1         = 0.0_ReKi 
   xd%Cn_v_minus1          = 0.0_ReKi
   xd%C_V_minus1           = 0.0_ReKi  ! This probably should not be set to 0.0, but should be set 

   
end subroutine UA_InitStates_Misc 
!==============================================================================   

!==============================================================================
subroutine UA_Init( InitInp, u, p, xd, OtherState, y,  m, Interval, &
                              InitOut,ErrStat, ErrMsg )
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
   type(UA_OtherStateType),      intent(  out)  :: OtherState  ! Initial other states
   type(UA_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   type(UA_MiscVarType),         intent(  out)  :: m           ! Initial misc/optimization variables
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
   character(ErrMsgLen)                         :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                               :: errStat2    ! temporary Error status of the operation
   character(*), parameter                      :: RoutineName = 'UA_Init'
   integer(IntKi)                               :: i,j, iNode, iOffset
   character(64)                                :: chanPrefix
   
      ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
   call DispNVD( UA_Ver )
   
      ! Allocate and set parameter data structure using initialization data
   call UA_SetParameters( interval, InitInp, p, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return    
   
      ! initialize the discrete states, other states, and misc variables
   call UA_InitStates_Misc( p, xd, OtherState, m, ErrStat2, ErrMsg2 )     ! initialize the continuous states
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return    
      
! TODO: Wrap all of this output handling in a  #ifdef block to avoid allocating arrays when coupled to FAST
#ifdef UA_OUTS   
      ! Allocate and set the InitOut data
   p%NumOuts = 29
   
   
   allocate(InitOut%WriteOutputHdr(p%NumOuts*p%numBlades*p%nNodesPerBlade),STAT=ErrStat2)
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'Error allocating WriteOutputHdr.',ErrStat,ErrMsg,RoutineName)
   allocate(InitOut%WriteOutputUnt(p%NumOuts*p%numBlades*p%nNodesPerBlade),STAT=ErrStat2)
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'Error allocating WriteOutputUnt.',ErrStat,ErrMsg,RoutineName)
   allocate(y%WriteOutput(p%NumOuts*p%numBlades*p%nNodesPerBlade),STAT=ErrStat2)
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'Error allocating y%WriteOutput.',ErrStat,ErrMsg,RoutineName)
   if (ErrStat >= AbortErrLev) return
         
   
   iNode = 0
   do j = 1,p%numBlades
      do i = 1,p%nNodesPerBlade
         
         iOffset = (i-1)*p%NumOuts + (j-1)*p%nNodesPerBlade*p%NumOuts 
         
         
         chanPrefix = "B"//trim(num2lstr(j))//"N"//trim(num2lstr(i))  
         InitOut%WriteOutputHdr(iOffset+ 1)  = 'ALPHA'//chanPrefix      
         InitOut%WriteOutputHdr(iOffset+ 2)  = 'VREL'//chanPrefix  
         InitOut%WriteOutputHdr(iOffset+ 3)  = 'CN'//chanPrefix      
         InitOut%WriteOutputHdr(iOffset+ 4)  = 'CC'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+ 5)  = 'CL'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+ 6)  = 'CD'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+ 7)  = 'CM'//chanPrefix     
         InitOut%WriteOutputHdr(iOffset+ 8)  = 'CNCP'//chanPrefix    
         InitOut%WriteOutputHdr(iOffset+ 9)  = 'CNIQ'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+10)  = 'CNPOT'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+11)  = 'DPP'//chanPrefix       
         InitOut%WriteOutputHdr(iOffset+12)  = 'CNP'//chanPrefix        
         InitOut%WriteOutputHdr(iOffset+13)  = 'FSP'//chanPrefix       
         InitOut%WriteOutputHdr(iOffset+14)  = 'DF'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+15)  = 'CNV'//chanPrefix         
         InitOut%WriteOutputHdr(iOffset+16)  = 'TAUV'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+17)  = 'LESF'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+18)  = 'TESF'//chanPrefix 
         InitOut%WriteOutputHdr(iOffset+19)  = 'VRTX'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+20)  = 'CVN'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+21)  = 'CMI'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+22)  = 'CMQ'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+23)  = 'CMV'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+24)  = 'AFEP'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+25)  = 'DFAF'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+26)  = 'PMC'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+27)  = 'T_f'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+28)  = 'T_V'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+29)  = 'dS'//chanPrefix
         InitOut%WriteOutputUnt(iOffset+1)  ='(deg)  '                                                
         InitOut%WriteOutputUnt(iOffset+2)  ='(m/s)   '                                                  
         InitOut%WriteOutputUnt(iOffset+3)  ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+4)  ='(-)   '                                                   
         InitOut%WriteOutputUnt(iOffset+5)  ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+6)  ='(-)   '                                                   
         InitOut%WriteOutputUnt(iOffset+7)  ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+8)  ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+9)  ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+10) ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+11) ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+12) ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+13) ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+14) ='(-)   '                                                    
         InitOut%WriteOutputUnt(iOffset+15) ='(-)   '                                                    
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
         !InitOut%WriteOutputHdr(iOffset25+1) =trim(chanPrefix)//'AOA'
         !InitOut%WriteOutputHdr(iOffset26+2) =trim(chanPrefix)//'Cn' 
         !InitOut%WriteOutputHdr(iOffset27+3) =trim(chanPrefix)//'Cc' 
         !InitOut%WriteOutputHdr(iOffset+4) =trim(chanPrefix)//'Cm' 
         !InitOut%WriteOutputHdr(iOffset+5) =trim(chanPrefix)//'Cl' 
         !InitOut%WriteOutputHdr(iOffset+6) =trim(chanPrefix)//'Cd' 
         !InitOut%WriteOutputHdr(iOffset+7)  =trim(chanPrefix)//'Cn_alf_circ'
         !InitOut%WriteOutputHdr(iOffset+8)  =trim(chanPrefix)//'Cn_alf_nc'
         !InitOut%WriteOutputHdr(iOffset+9)  =trim(chanPrefix)//'Cn_q_circ'
         !InitOut%WriteOutputHdr(iOffset+10) =trim(chanPrefix)//'Cn_q_nc'
         !InitOut%WriteOutputHdr(iOffset+11) =trim(chanPrefix)//'Cm_alf_circ'
         !InitOut%WriteOutputHdr(iOffset+12) =trim(chanPrefix)//'Cm_alf_nc'
         !InitOut%WriteOutputHdr(iOffset+13) =trim(chanPrefix)//'Cm_q_circ'
         !InitOut%WriteOutputHdr(iOffset+14) =trim(chanPrefix)//'Cm_q_nc'
         !InitOut%WriteOutputHdr(iOffset+15) =trim(chanPrefix)//'Alpha_e'
         !InitOut%WriteOutputHdr(iOffset+16) =trim(chanPrefix)//"f'"
         !InitOut%WriteOutputHdr(iOffset+17) =trim(chanPrefix)//"f'_c"
         !InitOut%WriteOutputHdr(iOffset+18) =trim(chanPrefix)//"f''"
         !InitOut%WriteOutputHdr(iOffset+19) =trim(chanPrefix)//"f''_c"
         !InitOut%WriteOutputHdr(iOffset+20) =trim(chanPrefix)//"f''_m"
         !InitOut%WriteOutputHdr(iOffset+21) =trim(chanPrefix)//'T_f'
         !InitOut%WriteOutputHdr(iOffset+22) =trim(chanPrefix)//'T_V'
         !InitOut%WriteOutputHdr(iOffset+23) =trim(chanPrefix)//'LESF'
         !InitOut%WriteOutputHdr(iOffset+24) =trim(chanPrefix)//'TESF'
         !InitOut%WriteOutputHdr(iOffset+25) =trim(chanPrefix)//'VRTX'
         !InitOut%WriteOutputHdr(iOffset+26) =trim(chanPrefix)//'tau_v'
         !InitOut%WriteOutputHdr(iOffset+27) = trim(chanPrefix)//'Cn_v'         
         !InitOut%WriteOutputHdr(iOffset+28) = trim(chanPrefix)//'C_V '         
         !InitOut%WriteOutputHdr(iOffset+29) = trim(chanPrefix)//'Cn_FS' 
         !InitOut%WriteOutputHdr(iOffset+30) = trim(chanPrefix)//'Cm_FS'
         !InitOut%WriteOutputUnt(iOffset+1)   ='(deg)  '
         !InitOut%WriteOutputUnt(iOffset+2) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+3) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+4) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+5) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+6) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+7) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+8) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+9) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+10) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+11) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+12) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+13) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+14) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+15) ='deg'
         !InitOut%WriteOutputUnt(iOffset+16) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+17) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+18) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+19) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+20) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+21) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+22) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+23) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+24) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+25) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+26) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+27) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+28) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+29) ='(-)   '
         !InitOut%WriteOutputUnt(iOffset+30) ='(-)   '
      end do
   end do
#else
   p%NumOuts = 0
#endif   
   
end subroutine UA_Init
!==============================================================================     
                              
                              
!============================================================================== 
subroutine UA_UpdateDiscOtherState( i, j, u, p, xd, OtherState, AFInfo, m, ErrStat, ErrMsg )   
! Routine for updating discrete states and other states (note it breaks the framework)
!..............................................................................
   
   integer   ,                   intent(in   )  :: i           ! node index within a blade
   integer   ,                   intent(in   )  :: j           ! blade index    
   type(UA_InputType),           intent(in   )  :: u           ! Inputs at Time                       
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters                                 
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Input: Discrete states at Time; 
                                                               ! Output: Discrete states at Time + Interval
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Other states  
   type(UA_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
   type(AFInfoType),             intent(in   )  :: AFInfo      ! The airfoil parameter data
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local Variables
  
   LOGICAL                                     :: LESF      ! logical flag indicating if leading edge separation is possible [-]
   LOGICAL                                     :: VRTX      ! logical flag indicating if a vortex is being processed [-]
   LOGICAL                                     :: TESF      ! logical flag indicating if trailing edge separation is possible [-]
      
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
   real(ReKi)                                   :: Kalpha_f          ! filtered backwards finite difference of alpha
   real(ReKi)                                   :: Kq_f              ! filtered backwards finite difference of q
   real(ReKi)                                   :: Kq
   real(ReKi)                                   :: q_cur
   real(ReKi)                                   :: q_f_cur
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
   real(ReKi)                                   :: alpha_filt_cur
   real(ReKi)                                   :: alpha0            ! zero lift angle of attack (radians)
   real(ReKi)                                   :: dalpha0
   real(ReKi)                                   :: alpha_f
   real(ReKi)                                   :: eta_e
   real(ReKi)                                   :: Kafactor
   real(ReKi)                                   :: Cn_q_circ
   real(ReKi)                                   :: Cn_q_nc
   real(ReKi)                                   :: St_sh
   real(ReKi)                                   :: T_sh, T_f, T_V
   real(ReKi)                                   :: Cn_alpha_q_nc
   real(ReKi)                                   :: Cn_prime_diff
      
   character(ErrMsgLen)                         :: errMsg2
   integer(IntKi)                               :: errStat2
   character(*), parameter                      :: RoutineName = 'UA_UpdateDiscOtherState'
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               

!bjj: use sigma1 and sigma3 states values from t      
   
      call ComputeKelvinChain(i, j, u, p, xd, OtherState, m, AFInfo, Cn_prime, Cn_prime_diff, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                                     St_sh, Kalpha_f, Kq_f, alpha_filt_cur, alpha_e, alpha0, dalpha0, alpha_f, eta_e, Kq, q_cur, q_f_cur, X1, X2, X3, X4, &
                                     Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, Df_c,Dalphaf, fprime, fprime_c, fprimeprime_c, fprimeprime_m, &
                                     Cn_alpha_q_nc, Cn_alpha_q_circ, Cn_q_circ, Cn_q_nc, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, T_f, T_V, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat>= AbortErrLev) return
            
      
      T_sh = 2.0_ReKi*(1.0_ReKi-fprimeprime) / St_sh
      
      !---------------------------------------------------------
      ! Update the OtherStates
      !---------------------------------------------------------
   
      LESF = (Cn_prime > Cn1) .or. ( Cn_prime < Cn2 ) ! LE separation can occur when this is .true.; assumption is that Cn2 <  0.0 and Cn1 > 0
      TESF = fprimeprime < xd%fprimeprime_minus1(i,j) ! Separation point is moving towards the Leading Edge when .true.; otherwise separation point is moving toward trailing edge   
            
      ! Process VRTX-related quantities
      !!!!!!!!!!!!!!!!!!!!!
      !! NEW CODE 2/19/2015
      VRTX = (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) > 0.0_ReKi)
      
      
      if (( xd%tau_V(i,j) >= (T_VL + T_sh) ) .and. TESF)  then !.and. (TESF) RRD added 
         xd%tau_V(i,j)     = 0.0_ReKi
        ! LESF=.FALSE. !also added this
      end if
      
!bjj: update sigma1 to value at t + dt      
      OtherState%sigma1(i,j) = 1.0_ReKi
      Kafactor      = Kalpha_f*dalpha0
     ! if ( fprimeprime < 0.7 .AND.  fprimeprime > 0.3 ) then 
         if ( TESF ) then  ! Separating flow
            if (Kafactor < 0.0_ReKi) then
               OtherState%sigma1(i,j) = 2.0_ReKi  ! This must be the first check
            else if (.not. LESF ) then
               OtherState%sigma1(i,j) = 1.0_ReKi !.4_ReKi    ! Leading edge separation has not occurred
            else if (xd%fprimeprime_minus1(i,j) <= 0.7_ReKi) then ! For this else, LESF = True
               OtherState%sigma1(i,j) = 2.0_ReKi !1.0_ReKi 
            else
               OtherState%sigma1(i,j) = 1.75_ReKi
            end if
         else ! Reattaching flow
        
            if (.not. LESF ) then
               OtherState%sigma1(i,j) = 0.5_ReKi
            end if

            if ( VRTX .and. (xd%tau_V(i,j) <= T_VL) ) then  ! Still shedding a vortex?
               OtherState%sigma1(i,j) = 0.25_ReKi
            end if
            if (Kafactor > 0.0_ReKi) then
               OtherState%sigma1(i,j) = 0.75_ReKi
            end if
            
         end if
      !!!if (.not. LESF ) then  !RRD: trying to emulate the old AD14 SEPAR.f90 with SHIFT=NOT(LESF), go back to original commented!!! out above when done
      !!!        OtherState%sigma1(i,j) = 0.667_ReKi
      !!!end if
      
!bjj: update sigma3 to value at t + dt      
         
      OtherState%sigma3(i,j) = 1.0_ReKi
      
      if ( (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) >= T_VL) ) then
         OtherState%sigma3(i,j) = 3.0_ReKi
         if (.not. TESF) then
            OtherState%sigma3(i,j) = 4.0_ReKi
            if ( VRTX .and. (xd%tau_V(i,j) <= T_VL) ) then
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
      
      !   ! We are testing this heirarchical logic instead of the above block 5/29/2015
      !if ( (xd%tau_V(i,j) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i,j) >= T_VL) ) then
      !   OtherState%sigma3(i,j) =  3.0_ReKi
      !else if (.not. TESF) then
      !   OtherState%sigma3(i,j) =  4.0_ReKi
      !else if ( VRTX .and. (xd%tau_V(i,j) <= T_VL) ) then
      !   if (Kafactor < 0.0_ReKi) then
      !      OtherState%sigma3(i,j) =  2.0_ReKi
      !   else
      !      OtherState%sigma3(i,j) = 1.0_ReKi
      !    end if           
      !else if (Kafactor < 0 ) then 
      !   OtherState%sigma3(i,j) =  4.0_ReKi
      !end if
      
      !!!if ( (.not. VRTX) .OR. (.not. TESF) ) then !RRD: trying to emulate the old AD14 SEPAR.f90 with SHIFT=NOT(TESF), go back to original commented!!! out above when done
      !!!       OtherState%sigma3(i,j) =  2_ReKi
      !!!   else
      !!!      OtherState%sigma3(i,j) = 1.0_ReKi
      !!!endif
      
      if ((.not. TESF) .and. (Kq_f*dalpha0 < 0.0_ReKi)) then
         OtherState%sigma3(i,j) = 1.0_ReKi
      end if
   
      
     
      !---------------------------------------------------------
      ! Update the Discrete States, xd
      !---------------------------------------------------------
     
           
      xd%alpha_minus1(i,j)        = u%alpha
      xd%alpha_filt_minus1(i,j)   = alpha_filt_cur
      xd%Kalpha_f_minus1(i,j)     = Kalpha_f
      xd%Kq_f_minus1(i,j)         = Kq_f     
      xd%q_f_minus1(i,j)          = q_f_cur    
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
      OtherState%FirstPass(i,j)   = .false.
   
      if ( xd%tau_V(i,j) > 0.0 .or. LESF ) then                !! TODO:  Verify this condition 2/20/2015 GJH
         xd%tau_V(i,j)          = xd%tau_V(i,j) + 2.0_ReKi*p%dt*u%U / p%c(i,j)  
      end if
   
      
#ifdef UA_OUTS
   m%TESF(i,j) = TESF  
   m%LESF(i,j) = LESF   
   m%VRTX(i,j) = VRTX 
#endif
      
      
end subroutine UA_UpdateDiscOtherState
!==============================================================================   

!============================================================================== 

!============================================================================== 
subroutine UA_UpdateStates( i, j, u, p, xd, OtherState, AFInfo, m, ErrStat, ErrMsg )
! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
! Continuous, constraint, discrete, and other states are updated to values at t + Interval.
!..............................................................................

   integer(IntKi),                intent(in   ) :: i               ! node index within a blade
   integer(IntKi),                intent(in   ) :: j               ! blade index 
   type(UA_InputType),            intent(in   ) :: u               ! Input at current timestep, t
   type(UA_ParameterType),        intent(in   ) :: p               ! Parameters
   type(UA_DiscreteStateType),    intent(inout) :: xd              ! Input: Discrete states at t;
                                                                   !   Output: Discrete states at t + Interval
   type(UA_OtherStateType),       intent(inout) :: OtherState      ! Input: Other states at t;
                                                                   !   Output: Other states at t + Interval
   type(UA_MiscVarType),          intent(inout) :: m               ! Misc/optimization variables
   type(AFInfoType),              intent(in   ) :: AFInfo          ! The airfoil parameter data
   integer(IntKi),                intent(  out) :: ErrStat         ! Error status of the operation
   character(*),                  intent(  out) :: ErrMsg          ! Error message if ErrStat /= ErrID_None

      ! Local variables  
      
   character(ErrMsgLen)                         :: errMsg2
   integer(IntKi)                               :: errStat2
   character(*), parameter                      :: RoutineName = 'UA_UpdateStates'

      ! Initialize variables

   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
         
   
   m%iBladeNode = i; 
   m%iBlade = j

      !BJJ: this seems to be the root cause of all sorts of numerical problems....
   IF (EqualRealNos(u%u, 0.0_ReKi) ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = 'UA_UpdateStates: U (air velocity magnitude relative to the airfoil) is zero.'
      RETURN
   END IF
   
   
      ! Update discrete states:
#ifdef DEBUG_v14
   call UA_UpdateDiscOtherState2( i, j, u, p, xd, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
#else
   call UA_UpdateDiscOtherState( i, j, u, p, xd, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
#endif
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   

end subroutine UA_UpdateStates
!==============================================================================   

!============================================================================== 
subroutine UA_CalcOutput( u, p, xd, OtherState, AFInfo, y, misc, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..............................................................................
   
   type(UA_InputType),           intent(in   )  :: u           ! Inputs at Time
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   )  :: OtherState  ! Other states at Time
   type(AFInfoType),             intent(in   )  :: AFInfo      ! The airfoil parameter data
   type(UA_OutputType),          intent(inout)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   type(UA_MiscVarType),         intent(inout)  :: misc        ! Misc/optimization variables
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   integer(IntKi)                                         :: errStat2        ! Error status of the operation (secondary error)
   character(ErrMsgLen)                                   :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   character(*), parameter                                :: RoutineName = 'UA_CalcOutput'
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
   real(ReKi)                                             :: Kalpha_f             ! backwards finite difference of alpha
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
   real(ReKi)                                             :: alpha_filt_cur
   real(ReKi)                                             :: alpha0              ! zero lift angle of attack (radians)
   real(ReKi)                                             :: dalpha0
   real(ReKi)                                             :: alpha_f
   real(ReKi)                                             :: eta_e
   real(ReKi)                                             :: St_sh, f_c_offset
   real(ReKi)                                             :: Cl_temp, Cd_temp, Cm_temp, Cn_temp, Cm_alpha_nc, Cn_q_circ, Cn_q_nc, T_f, T_V
   real(ReKi)                                             :: M, alpha1, alpha2, C_nalpha, C_nalpha_circ, T_f0, T_V0, T_p, b1,b2,b5,A1,A2,A5,S1,S2,S3,S4,k1_hat,f, k2_hat
   integer                                                :: iOffset
   real(ReKi)                                             :: Cn_alpha_q_nc, Cm_v, Cm_Lookup, alpha_prime_f, Cn_prime_diff, Kq_f,q_f_cur
   real(ReKi)                                             :: filtCutOff                    ! airfoil parameter for the low-pass cut-off frequency for pitching rate and accelerations (Hz)  
   real(ReKi)                                             :: x_cp_hat                      ! center-of-pressure distance from LE in chord fraction
   real(ReKi)                                             :: Cm_common                     ! 

   
   !BJJ: what are misc%iBladeNode, misc%iBlade here? I don't see them set in the routine that calls UA_CalcOutput
   
   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
   
  
  
   
   if (OtherState%FirstPass(misc%iBladeNode, misc%iBlade)) then
      
      
      call GetSteadyOutputs(AFInfo, u%alpha, y%Cl, y%Cd, y%Cm, Cd0, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
      
      
      y%Cn = y%Cl*cos(u%alpha) + (y%Cd-Cd0)*sin(u%alpha)
      y%Cc = y%Cl*sin(u%alpha) - (y%Cd-Cd0)*cos(u%alpha)
      Cm_v              = 0.0_ReKi
      Cn_alpha_q_circ   = 0.0_ReKi
      Cn_alpha_q_nc     = 0.0_ReKi
      Cn_pot            = 0.0_ReKi
      Dp                = 0.0_ReKi
      Cn_prime          = 0.0_ReKi
      fprime            = 0.0_ReKi
      Df                = 0.0_ReKi
      Cn_V              = 0.0_ReKi
      C_V               = 0.0_ReKi
      y%Cc              = 0.0_ReKi
      y%Cm              = 0.0_ReKi
      Cm_alpha_nc       = 0.0_ReKi
      Cm_q_circ         = 0.0_ReKi
      Cn_prime_diff     = 0.0_ReKi
      Cm_q_nc           = 0.0_ReKi
      alpha_prime_f     = 0.0_ReKi
      Dalphaf           = 0.0_ReKi
      Cm_temp         = 0.0_ReKi
      T_f               = 0.0_ReKi
      T_V               = 0.0_ReKi
      alpha_filt_cur    = 0.0_ReKi
      Cm_alpha_nc       = 0.0_ReKi
   else
      M           = u%U / p%a_s
      
      call UA_CheckMachNumber(M, misc%FirstWarn_M, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
      
      call AFI_GetAirfoilParams( AFInfo, M, u%Re, alpha0, alpha1, alpha2, eta_e, C_nalpha, C_nalpha_circ, &
                              T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, b5, A1, A2, A5, S1, S2, S3, S4, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, k1_hat, x_cp_bar, filtCutOff, ErrMsg2, ErrStat2 )           
            !AFI_GetAirfoilParams doens't return error, so I won't check. ! call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      
      call ComputeKelvinChain( misc%iBladeNode, misc%iBlade, u, p, xd, OtherState, misc, AFInfo, Cn_prime, Cn_prime_diff, Cn1, Cn2, Cd0, Cm0, k0, k1, k2, k3, T_VL, x_cp_bar, &
                                 St_sh, Kalpha_f, Kq_f, alpha_filt_cur, alpha_e, alpha0, dalpha0, alpha_f, eta_e, Kq, q_cur, q_f_cur, X1, X2, X3, X4, &
                                 Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, Df_c, Dalphaf, fprime, fprime_c, fprimeprime_c, fprimeprime_m, &
                                 Cn_alpha_q_nc, Cn_alpha_q_circ, Cn_q_circ, Cn_q_nc, Cm_q_circ, Cn_alpha_nc, Cm_q_nc, Cn_v, C_V, Cn_FS, T_f, T_V, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
            
         ! Eqn 1.53 or 1.53b depending on UAMod
      if (xd%tau_v(misc%iBladeNode, misc%iBlade) > 0.0) then
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
      
      if ( p%UAMod == 3 ) then
            ! Eqn 1.55a   
#ifdef TEST_THEORY
            y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset) + Cn_v*tan(alpha_e)*(1-xd%tau_v(misc%iBladeNode, misc%iBlade)/(T_VL))
#else            
            y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset) + Cn_v*    alpha_e *(1-xd%tau_v(misc%iBladeNode, misc%iBlade)/(T_VL)) !testing without tan 9/16/15
#endif            
        
         
      elseif ( p%UAMod == 2 ) then
            ! Eqn 1.55b
         y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset)
      else
         ! TODO: This is UAMod = 1 and we still need to verify these equations!!!  GJH Dec 20 2015
         if ( Cn_prime <= Cn1 ) then
            y%Cc = eta_e*Cc_pot*(sqrt(fprimeprime_c) - f_c_offset)*sin(alpha_e + alpha0)
         else
            if ( p%flookup ) then 
              ! if (p%UAMod == 1) then
              !    f      = Get_f_from_Lookup( p%UAMod, u%Re, u%alpha, alpha0, C_nalpha_circ, AFInfo, ErrStat, ErrMsg)
              ! else   
               ! TODO: Need to understand the effect of the offset on this f in the following equations and fprimeprime_c.
                  f      = Get_f_c_from_Lookup( u%Re, alpha_filt_cur, alpha0, C_nalpha_circ, AFInfo, ErrStat, ErrMsg)
              ! endif
               
            else
               f      = Get_f( alpha_filt_cur, alpha0, alpha1, alpha2, S1, S2, S3, S4 )
            end if
            
            k2_hat = 2*(Cn_prime-Cn1) + fprimeprime_c - f

            y%Cc   = k1_hat + Cc_pot*sqrt(fprimeprime_c)*fprimeprime_c**k2_hat*sin(alpha_e + alpha0)
            
         end if
         
      end if
      
         ! Eqn 1.2
      y%Cl = y%Cn*cos(u%alpha) + y%Cc*sin(u%alpha)
      y%Cd = y%Cn*sin(u%alpha) - y%Cc*cos(u%alpha) + Cd0
      
         ! Check for Cm column in AFInfo data        
      s1 = size(AFInfo%Table(1)%Coefs,2)
      if (s1 < 3) then      
         y%Cm = 0.0_ReKi
         Cm_alpha_nc = 0.0_ReKi        
         Cm_temp     = 0.0_ReKi
         Cm_v        = 0.0_ReKi
         alpha_prime_f = 0.0_ReKi
     
      else
            
         if ( p%UAMod == 2 ) then
            if (abs(alpha_f-alpha0) < .01) then
               if (alpha_f < alpha0) then
                  alpha_f = alpha_f - .01
               else
                  alpha_f = alpha_f + .01
               end if
            end if
         
            call GetSteadyOutputs(AFInfo, alpha_f, Cl_temp, Cd_temp, Cm_temp, Cd0, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            Cn_temp = Cl_temp*cos(alpha_f) + (Cd_temp-Cd0)*sin(alpha_f)
            
            ! TODO: In the future we are going to look at a delayed version of fprime_m and that is why we are using
            ! the variable name fprimeprime_m even though we have not introduced the delay, yet.  GJH 4/12/2016
            
!#ifndef CHECK_DENOM_DIFF
            if (abs(Cn_temp) < 0.01 ) then
!#else
!            if (EqualRealNos(Cn_temp,0.0_ReKi)) then
!#endif   
               fprimeprime_m = 1.0
            else            
               fprimeprime_m = (Cm_temp - Cm0) / Cn_temp 
            end if
         
         else
            fprimeprime_m = 0.0
         end if
      
            
         Cm_alpha_nc = - Cn_alpha_nc / 4.0_ReKi  ! Eqn 1.27
         Cm_common   = Cm_q_circ + Cm_alpha_nc + Cm_q_nc  ! Eqn 1.25 + 1.27 + 1.29a/1.31
         Cm_temp     = 0.0_ReKi
   
         if ( p%UAMod == 1 ) then
               ! Eqn 1.42
            x_cp_hat = k0 + k1*(1-fprimeprime) + k2*sin(pi*fprimeprime**k3) 
               ! Eqn 1.41
            Cm_FS  = Cm0 - Cn_alpha_q_circ*(x_cp_hat - 0.25_ReKi) + Cm_common

         elseif ( p%UAMod == 3 ) then
      
               ! Eqn 1.43a
            alpha_prime_f = alpha_f - Dalphaf
               ! Look up Cm using alpha_prime_f
            call GetSteadyOutputs(AFInfo, alpha_prime_f, Cl_temp, Cd_temp, Cm_temp, Cd0, ErrStat, ErrMsg)
               ! Eqn 1.44
            Cm_FS = Cm_temp + Cm_common
      
         else ! UAMod == 2
               ! Eqn 1.45   : NOTE: This differs from manual in that Cm0 is added because it is subtracted in our version of fprimeprime_m
            Cm_FS = Cm0 + Cn_FS*fprimeprime_m + Cm_common
            alpha_prime_f = 0.0_ReKi
         end if   
      
            ! Eqn 1.57
         Cm_v     = -x_cp_bar*( 1-cos( pi*xd%tau_v(misc%iBladeNode, misc%iBlade)/T_VL ) )*Cn_v
   
            ! Eqn 1.58 - 1.60 
         y%Cm   = Cm_FS + Cm_v 
      end if
   end if
   
#ifdef UA_OUTS
   iOffset = (misc%iBladeNode-1)*p%NumOuts + (misc%iBlade-1)*p%nNodesPerBlade*p%NumOuts 
   if (allocated(y%WriteOutput)) then  !bjj: because BEMT uses local variables for UA output, y%WriteOutput is not necessarially allocated. Need to figure out a better solution.
      y%WriteOutput(iOffset+ 1)    = alpha_filt_cur*180.0/pi                                                    
      y%WriteOutput(iOffset+ 2)    = u%U                                                                 
      y%WriteOutput(iOffset+ 3)    = y%Cn                                                                
      y%WriteOutput(iOffset+ 4)    = y%Cc                                                                
      y%WriteOutput(iOffset+ 5)    = y%Cl                                                                
      y%WriteOutput(iOffset+ 6)    = y%Cd                                                                
      y%WriteOutput(iOffset+ 7)    = y%Cm                                                                
      y%WriteOutput(iOffset+ 8)    = Cn_alpha_q_circ               ! CNCP in ADv14                                      
      y%WriteOutput(iOffset+ 9)    = Cn_alpha_q_nc                 ! CNIQ in ADv14                                      
      y%WriteOutput(iOffset+10)    = Cn_pot                                                              
      y%WriteOutput(iOffset+11)    = Dp                                                                  
      y%WriteOutput(iOffset+12)    = Cn_prime                                                            
      y%WriteOutput(iOffset+13)    = fprime                                                              
      y%WriteOutput(iOffset+14)    = Df                                                                  
      y%WriteOutput(iOffset+15)    = Cn_V                                                                
      y%WriteOutput(iOffset+16)    = xd%tau_v(misc%iBladeNode, misc%iBlade)                  
                                                                                                         
      if ( misc%LESF(misc%iBladeNode, misc%iBlade) ) then  ! BEDSEP in v14 !bjj: note that this output is calculated in UpdateDiscState, so it may be out of sync here.
         y%WriteOutput(iOffset+17) = 1.0_ReKi                                                            
      else                                                                                               
         y%WriteOutput(iOffset+17) = 0.0_ReKi                                                            
      end if                                                                                             
                                                                                                         
      if ( misc%TESF(misc%iBladeNode, misc%iBlade) ) then  !SHIFT in v14   !bjj: note that this output is calculated in UpdateDiscState, so it may be out of sync here.
         y%WriteOutput(iOffset+18) = 1.0_ReKi                                                            
      else                                                                                               
         y%WriteOutput(iOffset+18) = 0.0_ReKi                                                            
      end if                                                                                             
      if ( misc%VRTX(misc%iBladeNode, misc%iBlade) ) then  ! VOR in v14    !bjj: note that this output is calculated in UpdateDiscState, so it may be out of sync here.
         y%WriteOutput(iOffset+19) = 1.0_ReKi
      else
         y%WriteOutput(iOffset+19) = 0.0_ReKi
      end if
      y%WriteOutput(iOffset+20)    = C_V  
      y%WriteOutput(iOffset+21)    = Cm_alpha_nc
      y%WriteOutput(iOffset+22)    = Cm_q_nc 
      y%WriteOutput(iOffset+23)    = Cm_v
      y%WriteOutput(iOffset+24)    = alpha_prime_f
      y%WriteOutput(iOffset+25)    = Dalphaf
      y%WriteOutput(iOffset+26)    = Cm_temp
      y%WriteOutput(iOffset+27)    = T_f
      y%WriteOutput(iOffset+28)    = T_V
      y%WriteOutput(iOffset+29)    = 2.0_ReKi*u%U*p%dt/p%c(misc%iBladeNode, misc%iBlade)   
   
      !y%WriteOutput(iOffset+1) = u%alpha*180.0/pi
      !y%WriteOutput(iOffset+2) = y%Cn
      !y%WriteOutput(iOffset+3) = y%Cc
      !y%WriteOutput(iOffset+4) = y%Cm
      !y%WriteOutput(iOffset+5) = y%Cl
      !y%WriteOutput(iOffset+6) = y%Cd
      !y%WriteOutput(iOffset+7) = C_nalpha_circ*alpha_e  ! = Cn_alpha_circ
      !y%WriteOutput(iOffset+8) = Cn_alpha_nc
      !y%WriteOutput(iOffset+9) = Cn_q_circ
      !y%WriteOutput(iOffset+10) = Cn_q_nc
      !y%WriteOutput(iOffset+11) = 0.0  ! = Cm_alpha_circ  NOTE: we do not use this quantity
      !y%WriteOutput(iOffset+12) = Cm_alpha_nc
      !y%WriteOutput(iOffset+13) = Cm_q_circ
      !y%WriteOutput(iOffset+14) = Cm_q_nc
      !y%WriteOutput(iOffset+15) = alpha_e*180/pi
      !y%WriteOutput(iOffset+16) = fprime
      !y%WriteOutput(iOffset+17) = fprime_c
      !y%WriteOutput(iOffset+18) = fprimeprime
      !y%WriteOutput(iOffset+19) = fprimeprime_c
      !y%WriteOutput(iOffset+20) = fprimeprime_m
      !y%WriteOutput(iOffset+21) = T_f
      !y%WriteOutput(iOffset+22) = T_V
      !if ( OtherState%LESF(misc%iBladeNode, misc%iBlade) ) then
      !   y%WriteOutput(iOffset+23) = 1.0_ReKi
      !else
      !   y%WriteOutput(iOffset+23) = 0.0_ReKi
      !end if
      !
      !if ( misc%TESF(misc%iBladeNode, misc%iBlade) ) then
      !   y%WriteOutput(iOffset+24) = 1.0_ReKi
      !else
      !   y%WriteOutput(iOffset+24) = 0.0_ReKi
      !end if
      !if ( misc%VRTX(misc%iBladeNode, misc%iBlade) ) then
      !   y%WriteOutput(iOffset+25) = 1.0_ReKi
      !else
      !   y%WriteOutput(iOffset+25) = 0.0_ReKi
      !end if
      !y%WriteOutput(iOffset+26) = xd%tau_v(misc%iBladeNode, misc%iBlade)
      !y%WriteOutput(iOffset+27) = Cn_v         
      !y%WriteOutput(iOffset+28) = C_V          
      !y%WriteOutput(iOffset+29) = Cn_FS 
      !y%WriteOutput(iOffset+30) = Cm_FS
!      y%WriteOutput(iOffset+31) = u%U
   end if
#endif

!   if (allocated(y%WriteOutput)) then  !bjj: because BEMT uses local variables for UA output, y%WriteOutput is not necessarially allocated. Need to figure out a better solution.
!      y%WriteOutput(iOffset+1) = u%alpha*180.0/pi
!      y%WriteOutput(iOffset+2) = y%Cn
!      y%WriteOutput(iOffset+3) = y%Cc
!      y%WriteOutput(iOffset+4) = y%Cm
!      y%WriteOutput(iOffset+5) = y%Cl
!      y%WriteOutput(iOffset+6) = y%Cd
!      y%WriteOutput(iOffset+7) = C_nalpha_circ*alpha_e  ! = Cn_alpha_circ
!      y%WriteOutput(iOffset+8) = Cn_alpha_nc
!      y%WriteOutput(iOffset+9) = Cn_q_circ
!      y%WriteOutput(iOffset+10) = Cn_q_nc
!      y%WriteOutput(iOffset+11) = 0.0  ! = Cm_alpha_circ  NOTE: we do not use this quantity
!      y%WriteOutput(iOffset+12) = Cm_alpha_nc
!      y%WriteOutput(iOffset+13) = Cm_q_circ
!      y%WriteOutput(iOffset+14) = Cm_q_nc
!      y%WriteOutput(iOffset+15) = alpha_e*180/pi
!      y%WriteOutput(iOffset+16) = fprime
!      y%WriteOutput(iOffset+17) = fprime_c
!      y%WriteOutput(iOffset+18) = fprimeprime
!      y%WriteOutput(iOffset+19) = fprimeprime_c
!      y%WriteOutput(iOffset+20) = fprimeprime_m
!      y%WriteOutput(iOffset+21) = T_f
!      y%WriteOutput(iOffset+22) = T_V
!      if ( misc%LESF(misc%iBladeNode, misc%iBlade) ) then
!         y%WriteOutput(iOffset+23) = 1.0_ReKi
!      else
!         y%WriteOutput(iOffset+23) = 0.0_ReKi
!      end if
!      
!      if ( misc%TESF(misc%iBladeNode, misc%iBlade) ) then
!         y%WriteOutput(iOffset+24) = 1.0_ReKi
!      else
!         y%WriteOutput(iOffset+24) = 0.0_ReKi
!      end if
!      if ( misc%VRTX(misc%iBladeNode, misc%iBlade) ) then
!         y%WriteOutput(iOffset+25) = 1.0_ReKi
!      else
!         y%WriteOutput(iOffset+25) = 0.0_ReKi
!      end if
!      y%WriteOutput(iOffset+26) = xd%tau_v(misc%iBladeNode, misc%iBlade)
!      y%WriteOutput(iOffset+27) = Cn_v         
!      y%WriteOutput(iOffset+28) = C_V          
!      y%WriteOutput(iOffset+29) = Cn_FS 
!      y%WriteOutput(iOffset+30) = Cm_FS
!!      y%WriteOutput(iOffset+31) = u%U
!   end if
   
end subroutine UA_CalcOutput
!==============================================================================   
!> This subroutine checks that the Mach number is valid. If M > 0.3, the theory 
!! is invalid. If M > 1, numerical issues result in the code.
subroutine UA_CheckMachNumber(M, FirstWarn_M, ErrStat, ErrMsg )

   real(ReKi),                   intent(in   )  :: M           !< Mach number (-)
   logical,                      intent(inout)  :: FirstWarn_M !< is this the first warning about Mach number?
   integer(IntKi),               intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   if (abs(M) > 0.3_ReKi) then
      if (abs(M) >= 1.0_ReKi) then
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Mach number exceeds 1.0. Equations cannot be evaluated.'
      else         
         if ( FirstWarn_M ) then
            ErrStat = ErrID_Warn
            ErrMsg  = 'Mach number exceeds 0.3. Theory is invalid. This warning will not be repeated though the condition may persist.'
            FirstWarn_M = .false.
         else
            ErrStat = ErrID_None
            ErrMsg = "" 
         end if         
      end if      
   else
      ErrStat = ErrID_None
      ErrMsg = ""
   end if
   
   
end subroutine UA_CheckMachNumber


 
   
   
end module UnsteadyAero
