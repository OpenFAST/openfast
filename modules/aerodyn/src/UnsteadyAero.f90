!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
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
module UnsteadyAero

   ! This module uses equations defined in the document "The Unsteady Aerodynamics Module for FAST 8" by Rick Damiani and Greg Hayman, 28-Feb-2017

   use NWTC_Library
   use UnsteadyAero_Types
   use AirfoilInfo
   
   implicit none 

private


   public :: UA_Init
   public :: UA_UpdateStates
   public :: UA_CalcOutput
   public :: UA_CalcContStateDeriv
   public :: UA_End
   public :: UA_WriteOutputToFile

   public :: UA_ReInit
   public :: UA_InitStates_AllNodes ! used for AD linearization initialization

   integer(intki), parameter         :: UA_Baseline      = 1   ! UAMod = 1 [Baseline model (Original)]
   integer(intki), parameter         :: UA_Gonzalez      = 2   ! UAMod = 2 [Gonzalez's variant (changes in Cn,Cc,Cm)]
   integer(intki), parameter         :: UA_MinnemaPierce = 3   ! UAMod = 3 [Minnema/Pierce variant (changes in Cc and Cm)]
   integer(intki), parameter, public :: UA_HGM           = 4   ! UAMod = 4 [continuous variant of HGM (Hansen) model]
   integer(intki), parameter, public :: UA_OYE           = 5   ! UAMod = 5 [continuous Oye model]
   
   real(ReKi),     parameter :: Gonzalez_factor = 0.2_ReKi   ! this factor, proposed by Gonzalez (for "all" models) is used to modify Cc to account for negative values seen at f=0 (see Eqn 1.40)

   real(ReKi),     parameter, public :: UA_u_min = 0.01_ReKi   ! m/s; used to provide a minimum value so UA equations don't blow up (this should be much lower than range where UA is turned off)

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
!> Compute the deficiency function:
!! \f$ Y_n = Y_{n-1} \exp \left(-\frac{\Delta t}{T}\right)+\left(x_n - x_{n-1}\right)\exp\left(-\frac{\Delta t}{2T}\right)\f$
real(ReKi) function Get_ExpEqn( dt, T, Y_minus1, x, x_minus1 )
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi), intent(in   ) :: dt                              !< numerator of the exponent, \f$\Delta t \f$
   real(ReKi), intent(in   ) :: T                               !< denominator of the exponent (time constant), \f$T\f$
   real(ReKi), intent(in   ) :: Y_minus1                        !< previous value of the computed decay function, \f$ Y_{n-1} \f$
   real(ReKi), intent(in   ) :: x                               !< current value of x, \f$x_n\f$
   real(ReKi), intent(in   ) :: x_minus1                        !< previous value of x, \f$x_{n-1}\f$
   
   real(ReKi)                :: tmp
   
   tmp = -dt/T ! tmp should always be negative... should we check this, or are there some physics that make this check unnecessary?
   
   Get_ExpEqn = Y_minus1*exp(tmp) + (x - x_minus1)*exp(0.5_ReKi*tmp)

end function Get_ExpEqn
!==============================================================================     

   
!==============================================================================
! TE Flow Separation Equations                                                !
!==============================================================================

!==============================================================================
subroutine Get_f_from_Lookup( UAMod, Re, UserProp, alpha_in, alpha0, C_nalpha_circ, AFInfo, ErrStat, ErrMsg, f, cn_fs)
! Compute either fprime or fprimeprime using an analytical equation (and eventually a table lookup)
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   integer,                 intent(in   ) :: UAMod
   real(ReKi),              intent(in   ) :: Re                    ! Reynolds number
   real(ReKi),              intent(in   ) :: UserProp              ! User property for interpolating AFI
   real(ReKi),              intent(in   ) :: alpha_in              ! angle of attack (radians)
   real(ReKi),              intent(in   ) :: alpha0
   real(ReKi),              intent(in   ) :: C_nalpha_circ
   type(AFI_ParameterType), intent(in   ) :: AFInfo                ! The airfoil parameter data
   integer(IntKi),          intent(  out) :: ErrStat               ! Error status of the operation
   character(*),            intent(  out) :: ErrMsg                ! Error message if ErrStat /= ErrID_None
   real(ReKi),optional,     intent(  out) :: f
   real(ReKi),optional,     intent(  out) :: cn_fs
   
   real(ReKi)                             :: f_st
   
   real(ReKi)                             :: Cn, tmpRoot, denom
   type(AFI_OutputType)                   :: AFI_Interp


   real(ReKi)                       :: alpha         ! angle of attack (radians)
   real(ReKi)                       :: alpha_minus_alpha0

   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
    
   
   ! this routine is solving Eqn 1.32 for f:
   ! cn = C_nalpha_circ * (alpha-alpha0) * ((1+sqrt(f))/2)**2
   
   
   !bjj: if cn = 0 or c_nalpha_circ = 0 or alpha=alpha0, f has infinitely many solutions to this equation
   
      ! ensure that these angles are in appropriate ranges
   alpha  = alpha_in
   
   call MPi2Pi(alpha)

   alpha_minus_alpha0 = alpha - alpha0
   call MPi2Pi(alpha_minus_alpha0)

   call AFI_ComputeAirfoilCoefs( alpha, Re, UserProp, AFInfo, AFI_interp, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev ) return
   
   Cn =  AFI_interp%Cl*cos(alpha) + (AFI_interp%Cd-AFI_interp%Cd0)*sin(alpha)
   
   
   if (EqualRealNos( real(c_nalpha_circ,SiKi), 0.0_SiKi )) then
      tmpRoot = 0.0_ReKi
   else if (EqualRealNos( real(alpha,SiKi), real(alpha0,SiKi) )) then
      tmpRoot = 0.0_ReKi
   else
      
   
      if (EqualRealNos( real(cn,SiKi), 0.0_SiKi )) then
         tmpRoot = 0.0_ReKi
      else 
      
         denom = (C_nalpha_circ*alpha_minus_alpha0)
      
         !    bjj: if tmpRoot=cn/(C_nalpha_circ*(alpha-alpha0)) is negative, this whole equation is bogus....
         tmpRoot  = Cn/denom
   
         if (tmpRoot < 0.0_ReKi) then
            tmpRoot = 0.0_ReKi
         end if
                        
      end if
            
   end if
   
   if (UAMod == UA_Gonzalez) then
      f_st = ((3.0_ReKi * sqrt(tmpRoot)-1.0_ReKi)/ 2.0_ReKi )**2
   else
      f_st = ( 2.0_ReKi * sqrt( tmpRoot ) - 1.0_ReKi )**2 
   end if
   
   
   if ( f_st > 1.0 ) then
      f_st = 1.0_ReKi
   end if
   
   if (present(f)) then
      f = f_st
   end if

   if (present(cn_fs)) then
      if (equalRealNos(f_st,1.0_ReKi)) then
         cn_fs = 0.0_ReKi
      else
         cn_fs = (Cn - C_nalpha_circ * alpha_minus_alpha0 * f_st) / (1.0_ReKi-f_st); ! modification by Envision Energy
      end if   
   end if
   
end subroutine Get_f_from_Lookup      
!==============================================================================

!==============================================================================
real(ReKi) function Get_f_c_from_Lookup( UAMod, Re, UserProp, alpha_in, alpha0_in, c_nalpha_circ, eta_e, AFInfo, ErrStat, ErrMsg)
! Compute either fprime or fprimeprime using an analytical equation (and eventually a table lookup)
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................
   integer,          intent(in   ) :: UAMod
   real(ReKi),       intent(in   ) :: Re            ! Reynolds number
   real(ReKi),       intent(in   ) :: UserProp      ! User property for 2D AFI interpolation
   real(ReKi),       intent(in   ) :: alpha_in      ! angle of attack (radians)
   real(ReKi),       intent(in   ) :: alpha0_in
   real(ReKi),       intent(in   ) :: c_nalpha_circ
   real(ReKi),       intent(in   ) :: eta_e
   type(AFI_ParameterType), intent(in   ) :: AFInfo        ! The airfoil parameter data   
   integer(IntKi),   intent(  out) :: ErrStat       ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
   
   
   real(ReKi), parameter           :: fc_limit = (1.0_ReKi + Gonzalez_factor)**2    ! normally, fc is limited by 1, but we're limiting (sqrt(fc)-Gonzalez_factor) to 1, so fc is limited to 1.44 instead (when Gonzalez_factor is 0.2)
   real(ReKi)                      :: Cc, denom
   type(AFI_OutputType)            :: AFI_Interp

   real(ReKi)                      :: alpha         ! angle of attack (radians)
   real(ReKi)                      :: alpha0

   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
      ! ensure that these angles are in appropriate ranges
   alpha  = alpha_in
   alpha0 = alpha0_in
   
   call MPi2Pi(alpha)
   call MPi2Pi(alpha0)

      ! in cases where denom is zero, Get_f_c_from_Lookup = min(fc_limit, inf)
   if (EqualRealNos(real(alpha,SiKi), 0.0_SiKi)) then
      
      Get_f_c_from_Lookup = fc_limit
      
   elseif (EqualRealNos(real(alpha,SiKi), real(alpha0,SiKi))) then
      
      Get_f_c_from_Lookup = fc_limit
      
   else if (EqualRealNos( real(c_nalpha_circ,SiKi), 0.0_SiKi )) then
      
      Get_f_c_from_Lookup = fc_limit

   else
         
      call AFI_ComputeAirfoilCoefs( alpha, Re, UserProp,  AFInfo, AFI_interp, ErrStat, ErrMsg)
         if (ErrStat >= AbortErrLev) return
   
      Cc =  AFI_interp%Cl*sin(alpha) - (AFI_interp%Cd-AFI_interp%Cd0)*cos(alpha)
   

      if (UAMod == UA_Gonzalez) then
         denom = eta_e*c_nalpha_circ*( alpha-alpha0 )*(alpha)    !NOTE: Added back (alpha) because idling cases with alpha 90-degrees show problems with tan(alpha), the code should match steady state if the formulation in the calculation of Cc is in agreement with this formulation
      else
         denom = eta_e*c_nalpha_circ*( alpha-alpha0 )*tan(alpha)
      endif
      Get_f_c_from_Lookup =  min(fc_limit, (  Cc / denom  + Gonzalez_factor ) **2 )
   end if
      ! Apply an offset of Gonzalez_factor = 0.2 to fix cases where f_c should be negative, but we are using **2 so can only return positive values
      ! Note: because we include this offset, it must be accounted for in the final value of Cc, eqn 1.40.  This will be applied
      ! For both UA_Mod = 1,2, and 3 when using Flookup = T
   
   
end function Get_f_c_from_Lookup      

!==============================================================================
real(ReKi) function Get_f( alpha, alpha0, alpha1, alpha2, S1, S2, S3, S4 )
! Compute fprime, implements Eqn 1.33
! Called by : ComputeKelvinChain
!..............................................................................

   real(ReKi), intent(in   ) :: alpha         ! angle of attack (radians)
   real(ReKi), intent(in   ) :: alpha0        ! zero lift angle of attack (radians)
   real(ReKi), intent(in   ) :: alpha1        ! angle of attack at f = 0.7, approximately the stall angle; for alpha >= alpha0 (radians)
   real(ReKi), intent(in   ) :: alpha2        ! angle of attack at f = 0.7, approximately the stall angle; for alpha < alpha0 (radians)
   real(ReKi), intent(in   ) :: S1            ! constant in the f-curve best-fit, alpha >= alpha0  
   real(ReKi), intent(in   ) :: S2            ! constant in the f-curve best-fit, alpha >= alpha0  
   real(ReKi), intent(in   ) :: S3            ! constant in the f-curve best-fit, alpha < alpha0    
   real(ReKi), intent(in   ) :: S4            ! constant in the f-curve best-fit, alpha < alpha0   

      
      ! Implements Equation 1.33
   if (alpha > alpha1) then
      Get_f = 0.04_ReKi + 0.66_ReKi*exp((alpha1-alpha)/S2)
   else if (alpha < alpha2) then
      Get_f = 0.04_ReKi + 0.66_ReKi*exp((alpha-alpha2)/S4)
   else if (alpha>= alpha0) then !alpha0<=alpha<=alpha1
      Get_f = 1.0_ReKi - 0.3_ReKi*exp((alpha-alpha1)/S1)
   else ! alpha2<= alpha < alpha0
      Get_f = 1.0_ReKi - 0.3_ReKi*exp((alpha2-alpha)/S3)
   end if
   
end function Get_f
!==============================================================================                              
subroutine ComputeKelvinChain( i, j, u, p, xd, OtherState, misc, AFInfo, KC, BL_p, ErrStat, ErrMsg )
! 
! Called by : DRIVER
! Calls  to : Get_Beta_M_Sqrd, Get_Beta_M, AFI_ComputeUACoefs, Get_ds, Get_Pitchrate, Get_k_, Get_Kupper, Get_ExpEqn
!             Get_Cn_nc, Get_alpha_e, Get_Cm_q_circ, Get_Cm_q_nc, Get_f, Get_Cn_FS, Get_C_V, Get_Cn_v
!...............................................................................

      
   integer(IntKi),                         intent(in   ) :: i,j               ! Blade node indices
   type(UA_InputType),                     intent(in   ) :: u                 ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   type(UA_ParameterType),                 intent(in   ) :: p                 ! Parameters   
   type(UA_DiscreteStateType),             intent(in   ) :: xd                ! Input: Discrete states at t;
   type(UA_OtherStateType),                intent(in   ) :: OtherState        ! Other states at t
   type(UA_MiscVarType),                   intent(inout) :: misc              ! Misc/optimization variables
   type(AFI_ParameterType),                intent(in   ) :: AFInfo            ! The airfoil parameter data
   
   
   type(AFI_UA_BL_Type),                   intent(  out) :: BL_p

   integer(IntKi),                         intent(  out) :: ErrStat           ! Error status of the operation
   character(*),                           intent(  out) :: ErrMsg            ! Error message if ErrStat /= ErrID_None

            
   type(UA_KelvinChainType),               intent(  out) :: KC

   real(ReKi)                :: M                                             ! Mach number (-)
   real(ReKi)                :: beta_M                                        ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2)
   real(ReKi)                :: beta_M_Sqrd                                   ! square of the Prandtl-Glauert compressibility correction factor,  (1-M**2)
                 
   real(ReKi)                :: Cn_temp
   real(ReKi)                :: Cn_fs_temp

   type(AFI_OutputType)      :: AFI_interp                 !  Cl, Cd, Cm, Cpmin

   real(ReKi)                :: T_I                                           !
   real(ReKi)                :: Kalpha                                        !
   real(ReKi)                :: Kalpha_f_minus1                               !
   real(ReKi)                :: Kq_f_minus1                                   !
   real(ReKi)                :: alpha_minus1                                  !
   real(ReKi)                :: alpha_filt_minus1                                  !
   real(ReKi)                :: q_minus1                                      !
   real(ReKi)                :: q_f_minus1
   real(ReKi)                :: Cn_pot_minus1                                 !
   real(ReKi)                :: K3prime_q                                     !
   real(ReKi)                :: k_mq                                          !
   real(ReKi)                :: Kprimeprime_q                                 !
   real(ReKi)                :: dynamicFilterCutoffHz                         ! find frequency based on reduced frequency of k = BL_p%filtCutOff
   real(ReKi)                :: LowPassConst
   

#ifdef TEST_THEORY
   real(ReKi)                                            :: Cn_prime_diff
#endif   
   
   integer(IntKi)            :: ErrStat2
   character(ErrMsgLen)      :: ErrMsg2
   character(*), parameter   :: RoutineName = 'ComputeKelvinChain'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   KC%fprimeprime_m = 0

   M           = u%U / p%a_s
   call UA_CheckMachNumber(M, misc%FirstWarn_M, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
      
   beta_M_Sqrd = 1.0_ReKi - M**2
   beta_M      = sqrt(beta_M_Sqrd) 
   
   ! Lookup values using Airfoil Info module
   call AFI_ComputeUACoefs( AFInfo, u%Re, u%UserProp, BL_p, ErrMsg2, ErrStat2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
   KC%C_nalpha_circ  =  BL_p%C_nalpha / beta_M
      
      ! Override eta_e if we are using Flookup
   if ( p%Flookup ) then
      BL_p%eta_e = 1.0
   end if
   
   ! Kelvin chain
   KC%ds       = 2.0_ReKi*u%U*p%dt/p%c(i,j)                             ! Eqn 1.5b
   
   if (OtherState%FirstPass(i,j)) then
      alpha_minus1 = u%alpha      
      alpha_filt_minus1 = u%alpha     
   else
      alpha_minus1 = xd%alpha_minus1(i,j)     
      alpha_filt_minus1 = xd%alpha_filt_minus1(i,j)   
   end if
   
   ! Low Pass filtering of alpha, q, and Kq in order to deal with numerical issues that arise for very small values of dt
   ! See equations 1.7 and 1.8 of the manual
   ! This filter is a Simple Infinite Impulse Response Filter
   ! See https://en.wikipedia.org/wiki/Low-pass_filter#Simple_infinite_impulse_response_filter
   
   dynamicFilterCutoffHz = max( 1.0_ReKi, u%U ) * BL_p%filtCutOff / PI / p%C(i,j)
   LowPassConst  =  exp(-2.0_ReKi*PI*p%dt*dynamicFilterCutoffHz) ! from Eqn 1.8 [7]
   
   KC%alpha_filt_cur = LowPassConst*alpha_filt_minus1 + (1.0_ReKi-LowPassConst)*u%alpha ! from eq 1.8 [1: typo in documentation, though]
   
   
   KC%dalpha0  = KC%alpha_filt_cur - BL_p%alpha0
   
    
      ! Compute Kalpha using Eqn 1.7
  
   Kalpha        = ( KC%alpha_filt_cur - alpha_filt_minus1 ) / p%dt ! Eqn 1.7, using filtered values of alpha
   
   if (OtherState%FirstPass(i,j)) then
      Kalpha_f_minus1 = 0.0_ReKi
   else
      Kalpha_f_minus1 = xd%Kalpha_f_minus1(i,j)
   end if
   
   KC%q_cur         = Kalpha   * p%c(i,j) / u%U   ! Kalpha here uses the low-pass filtered alphas (Eqn 1.8 [2])
   
       
   if (OtherState%FirstPass(i,j)) then
      q_minus1   = KC%q_cur 
      q_f_minus1 = KC%q_cur    
   else    
      q_minus1   = xd%q_minus1(i,j) 
      q_f_minus1 = xd%q_f_minus1(i,j) 
   end if  
   
   KC%q_f_cur       = LowPassConst*q_f_minus1 + (1.0_ReKi-LowPassConst)*KC%q_cur ! (Eqn 1.8 [3])
      
      
#ifdef TEST_THEORY
      ! Change found in ADv14
   KC%q_f_cur = SAT( KC%q_f_cur, 0.03_ReKi, 0.1_ReKi )
#endif   
   
   KC%Kalpha_f      = KC%q_f_cur * u%U / p%c(i,j)  ! Kalpha_f is using the low-pass filtered q (Eqn. 1.8 [4])

   
   ! Compute Kq  using Eqn 1.8  with time-shifted q s
#ifdef TEST_THEORY
   KC%Kq            = ( KC%q_f_cur  - q_f_minus1 ) / p%dt !bjj: jmj questions if this should be the way it's implemented
#else
   
   KC%Kq            = ( KC%q_cur  - q_minus1 ) / p%dt ! Eqn. 1.8 [5]
#endif

   if (OtherState%FirstPass(i,j)) then
      Kq_f_minus1 = 0.0_ReKi
   else
      Kq_f_minus1 = xd%Kq_f_minus1(i,j)
   end if

   KC%Kq_f          =  LowPassConst*Kq_f_minus1 + (1.0_ReKi-LowPassConst)*KC%Kq !Eqn. 1.8 [6]
   
   
   !bjj: todo: should we check the denominator to make sure it doesn't go to 0?

   KC%k_alpha = 1.0_ReKi / ( (1.0_ReKi - M) + (BL_p%C_nalpha/2.0_ReKi) * M**2 * beta_M * (BL_p%A1*BL_p%b1 + BL_p%A2*BL_p%b2) )  ! Eqn 1.11a
   KC%k_q     = 1.0_ReKi / ( (1.0_ReKi - M) +  BL_p%C_nalpha           * M**2 * beta_M * (BL_p%A1*BL_p%b1 + BL_p%A2*BL_p%b2) )  ! Eqn 1.11b   
   T_I        = p%c(i,j) / p%a_s                                                                                                ! Eqn 1.11c
   
   KC%T_alpha  = T_I * KC%k_alpha * 0.75                                                                                      ! Eqn 1.10a
   KC%T_q      = T_I * KC%k_q * 0.75                                                                                          ! Eqn 1.10b
      
   KC%T_f           = BL_p%T_f0 / OtherState%sigma1(i,j)                                                                      ! Eqn 1.37          
   KC%T_fc          = BL_p%T_f0 / OtherState%sigma1c(i,j)      ! NOTE: Added equations for time constants of fc (for Cc) and fm (for Cm) with UAMod=2
   KC%T_fm          = BL_p%T_f0 / OtherState%sigma1m(i,j)
  
   KC%Kprime_alpha  = Get_ExpEqn( real(p%dt,ReKi), KC%T_alpha, xd%Kprime_alpha_minus1(i,j), KC%Kalpha_f, Kalpha_f_minus1 )    ! Eqn 1.18b
   KC%Cn_alpha_nc   = 4.0_ReKi*KC%T_alpha * ( KC%Kalpha_f - KC%Kprime_alpha ) / M                                             ! Eqn 1.18a
   
   KC%Kprime_q      = Get_ExpEqn( real(p%dt,ReKi), KC%T_q    , xd%Kprime_q_minus1(i,j)    ,  KC%Kq_f   , Kq_f_minus1     )    ! Eqn 1.19b 
   KC%Cn_q_nc       = -1.0_ReKi*KC%T_q * ( KC%Kq_f - KC%Kprime_q ) / M                                                        ! Eqn 1.19a
         
   KC%Cn_alpha_q_nc = KC%Cn_alpha_nc + KC%Cn_q_nc                                                                             ! Eqn 1.17
   
if (p%ShedEffect) then
   KC%X1            = Get_ExpEqn( KC%ds*beta_M_Sqrd*BL_p%b1, 1.0_ReKi, xd%X1_minus1(i,j), BL_p%A1*(KC%alpha_filt_cur - alpha_filt_minus1), 0.0_ReKi ) ! Eqn 1.15a
   KC%X2            = Get_ExpEqn( KC%ds*beta_M_Sqrd*BL_p%b2, 1.0_ReKi, xd%X2_minus1(i,j), BL_p%A2*(KC%alpha_filt_cur - alpha_filt_minus1), 0.0_ReKi ) ! Eqn 1.15b
else
   KC%X1  = 0.0_ReKi  ! u%alpha (and alpha_filt_cur) contains shed vorticity effect already 
   KC%X2  = 0.0_ReKi  ! so that alpha_e = u%alpha-alpha0 directly
endif
   
   KC%alpha_e       = (KC%alpha_filt_cur - BL_p%alpha0) - KC%X1 - KC%X2                                                       ! Eqn 1.14
   
   KC%Cn_alpha_q_circ = KC%C_nalpha_circ * KC%alpha_e                                                                         ! Eqn 1.13
   
   if ( p%UAMod == UA_Gonzalez ) then
         ! Compute X3 and X4 using Eqn 1.16a  and then add Cn_q_circ (Eqn 1.16) to the previously computed Cn_alpha_q_circ
if (p%ShedEffect) then
      KC%X3              = Get_ExpEqn( KC%ds*beta_M_Sqrd*BL_p%b1, 1.0_ReKi, xd%X3_minus1(i,j), BL_p%A1*(KC%q_f_cur - q_f_minus1), 0.0_ReKi ) ! Eqn 1.16a [1]
      KC%X4              = Get_ExpEqn( KC%ds*beta_M_Sqrd*BL_p%b2, 1.0_ReKi, xd%X4_minus1(i,j), BL_p%A2*(KC%q_f_cur - q_f_minus1), 0.0_ReKi ) ! Eqn 1.16a [2]
else
      KC%X3 = 0.0_ReKi ! Similar to X1 and X2, we assumed that this effect is already included
      KC%X4 = 0.0_ReKi
endif
      
      KC%Cn_q_circ       = KC%C_nalpha_circ*KC%q_f_cur/2.0 - KC%X3 - KC%X4                                                    ! Eqn 1.16

   else ! these aren't used (they are possibly output to UA output file (when UA_OUTS defined) file, though)
      KC%X3              = 0.0_ReKi
      KC%X4              = 0.0_ReKi
      KC%Cn_q_circ       = 0.0_ReKi
   end if
   
   K3prime_q       = Get_ExpEqn( BL_p%b5*beta_M_Sqrd*KC%ds, 1.0_ReKi, xd%K3prime_q_minus1(i,j),  BL_p%A5*(KC%q_f_cur - q_f_minus1), 0.0_ReKi )  ! Eqn 1.26
   KC%Cm_q_circ    = -BL_p%C_nalpha*(KC%q_f_cur - K3prime_q)*p%c(i,j)/(16.0_ReKi*beta_M*u%U)                                  ! Eqn 1.25
   
   KC%Cn_pot       = KC%Cn_alpha_q_circ + KC%Cn_alpha_q_nc                                                                    ! Eqn 1.20 [2a]
   
   k_mq            = 7.0_ReKi / (15.0_ReKi*(1.0_ReKi-M) + 1.5_ReKi * BL_p%C_nalpha * BL_p%A5 * BL_p%b5 * beta_M * M**2)       ! Eqn 1.29 [2]      ! CHECK THAT DENOM ISN'T ZERO!
   Kprimeprime_q   = Get_ExpEqn( real(p%dt,ReKi), k_mq**2*T_I , xd%Kprimeprime_q_minus1(i,j) ,  KC%Kq_f , Kq_f_minus1  )      ! Eqn 1.29 [3]
   
      ! Compute Cm_q_nc 
   if ( p%UAMod == UA_MinnemaPierce ) then
      KC%Cm_q_nc =  -1.0_ReKi * KC%Cn_q_nc / 4.0_ReKi - (KC%k_alpha**2) * T_I * (KC%Kq_f - Kprimeprime_q) / (3.0_ReKi*M)      ! Eqn 1.31
   else  
      KC%Cm_q_nc = -7.0_ReKi * (k_mq**2) * T_I * (KC%Kq_f - Kprimeprime_q) / (12.0_ReKi*M)                                    ! Eqn 1.29 [1]       
   end if
   
   if ( p%UAMod == UA_Gonzalez ) then
      KC%Cc_pot = KC%C_nalpha_circ * KC%alpha_e * u%alpha  !Added this equation with (u%alpha) instead of tan(alpha_e+alpha0). First, tangent gives problems in idling conditions at angles of attack of 90 degrees. Second, the angle there is a physical concept according to the original BL model, and u%alpha could be more suitable 
   else   
   !!! THIS IS A PROBLEM IF KC%alpha_e+BL_p%alpha0 ARE NEAR +/-PI/2
      KC%Cc_pot = KC%Cn_alpha_q_circ * tan(KC%alpha_e+BL_p%alpha0)                                                           ! Eqn 1.21 with cn_pot_circ=KC%Cn_alpha_q_circ as from Eqn 1.20 [3]  
   endif
   
   if (OtherState%FirstPass(i,j)) then
      Cn_pot_minus1 = KC%Cn_pot
   else
      Cn_pot_minus1 = xd%Cn_pot_minus1(i,j)
   end if
   
   KC%Dp            = Get_ExpEqn( KC%ds, BL_p%T_p, xd%Dp_minus1(i,j), KC%Cn_pot, Cn_pot_minus1 )                              ! Eqn 1.35b
   KC%Cn_prime      = KC%Cn_Pot - KC%Dp                                                                                       ! Eqn 1.35a
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! This code is taken from ADv14 but doesn't reflect the original intent of the UA theory document
#ifdef TEST_THEORY
   if (OtherState%FirstPass(i,j)) then
      Cn_prime_diff = 0.0_ReKi
   else
      Cn_prime_diff = KC%Cn_prime - xd%Cn_prime_minus1(i,j)
   end if
      
IF ( p%UAMod /= UA_Gonzalez ) THEN
   IF ( KC%alpha_filt_cur * Cn_prime_diff < 0. ) THEN

      KC%T_f   = BL_p%T_f0*1.5
   ELSE

      KC%T_f   = BL_p%T_f0
   ENDIF
ENDIF
#endif   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   KC%alpha_f       = KC%Cn_prime / KC%C_nalpha_circ + BL_p%alpha0                                                            ! Eqn 1.34
   
   if (p%flookup) then
      call Get_f_from_Lookup( p%UAMod, u%Re, u%UserProp, KC%alpha_f, BL_p%alpha0, KC%C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2, f=KC%fprime)     ! Solve Eqn 1.32a for f (=KC%fprime) when alpha is replaced with alpha_f (see issue when KC%C_nalpha_circ is 0) 
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   else   
      KC%fprime = Get_f( KC%alpha_f, BL_p%alpha0, BL_p%alpha1, BL_p%alpha2, BL_p%S1, BL_p%S2, BL_p%S3, BL_p%S4)               ! Eqn 1.33
   end if
   
   if (OtherState%FirstPass(i,j)) then
      KC%Df = 0.0_ReKi
   else
      KC%Df = Get_ExpEqn( KC%ds, KC%T_f, xd%Df_minus1(i,j), KC%fprime, xd%fprime_minus1(i,j) )                                ! Eqn 1.36b
   end if
      
   KC%fprimeprime   = KC%fprime - KC%Df                                                                                       ! Eqn 1.36a
   
   if (p%Flookup) then
         ! Compute fprime using Eqn 1.32 and Eqn 1.33
      KC%fprime_c   = Get_f_c_from_Lookup( p%UAMod, u%Re, u%UserProp, KC%alpha_f, BL_p%alpha0, KC%C_nalpha_circ, BL_p%eta_e, AFInfo, ErrStat2, ErrMsg2) ! Solve Eqn 1.32b for f when alpha is replaced with alpha_f
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return

      if ( p%UAMod == UA_Gonzalez ) then   !Added this part of the code to obtain fm
         call AFI_ComputeAirfoilCoefs( KC%alpha_f, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat2, ErrMsg2)
           call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
           if (ErrStat >= AbortErrLev) return

         Cn_temp = AFI_interp%Cl*cos(KC%alpha_f) + (AFI_interp%Cd - AFI_interp%Cd0)*sin(KC%alpha_f)
         if (abs(Cn_temp) < 0.01_ReKi ) then
            KC%fprime_m = 0.0_ReKi
         else
            KC%fprime_m = (AFI_interp%Cm - AFI_interp%Cm0) / Cn_temp 
         end if
      else
         KC%fprime_m = 0.0_ReKi
      endif

     
      if (OtherState%FirstPass(i,j)) then
         KC%Df_c = 0.0_ReKi
         KC%Df_m = 0.0_ReKi
      else
         KC%Df_c = Get_ExpEqn( KC%ds, KC%T_fc, xd%Df_c_minus1(i,j), KC%fprime_c, xd%fprime_c_minus1(i,j)  )
         KC%Df_m = Get_ExpEqn( KC%ds, KC%T_fm, xd%Df_m_minus1(i,j), KC%fprime_m, xd%fprime_m_minus1(i,j)  )  ! used in UAMod=UA_Gonzalez only
      end if
   
         ! Compute Df using Eqn 1.36b   
   
         ! Compute fprimeprime using Eqn 1.36a
      KC%fprimeprime_c   = KC%fprime_c - KC%Df_c
      
      IF ( p%UAMod == UA_Gonzalez ) THEN
         KC%fprimeprime_m   = KC%fprime_m - KC%Df_m
      END IF
   else
      KC%fprime_c = KC%fprime
      KC%Df_c = KC%Df
      KC%fprimeprime_c = KC%fprimeprime

         ! variables used for UAMod=UA_Gonzalez
      KC%fprime_m = 0.0_ReKi
      KC%Df_m     = 0.0_ReKi
      KC%fprimeprime_m = KC%fprimeprime
   end if
   
   
   if ( p%UAMod == UA_Gonzalez ) then
      KC%Cn_FS   = KC%Cn_alpha_q_nc + KC%Cn_q_circ + KC%Cn_alpha_q_circ *  ( (1.0_ReKi + 2.0_ReKi*sqrt(KC%fprimeprime) ) / 3.0_ReKi )**2     ! Eqn 1.39
   else
   ! change proposed by Pariya:
     ! KC%Cn_FS   = KC%Cn_alpha_q_nc                + KC%Cn_alpha_q_circ *  ( (1.0_ReKi +          sqrt(KC%fprimeprime) ) / 2.0_ReKi )**2     ! Eqn 1.38
      call Get_f_from_Lookup( p%UAMod, u%Re, u%UserProp, KC%alpha_e+BL_p%alpha0, BL_p%alpha0, KC%C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2, cn_fs=Cn_fs_temp)
      KC%Cn_FS  = KC%Cn_alpha_q_nc + KC%C_nalpha_circ * KC%alpha_e*KC%fprimeprime  + Cn_fs_temp*(1-KC%fprimeprime)
   end if
      
   if ( p%UAMod == UA_MinnemaPierce ) then
      if (OtherState%FirstPass(i,j)) then     
         KC%Dalphaf    = 0.0_ReKi
      else
         KC%Dalphaf    = Get_ExpEqn( KC%ds, 0.1_ReKi*KC%T_f, xd%Dalphaf_minus1(i,j), KC%alpha_f, xd%alphaf_minus1(i,j) )         ! Eqn 1.43
      end if
   else
      KC%Dalphaf    = 0.0_ReKi
   end if
   
   if ( p%UAMod == UA_Gonzalez ) then
      KC%C_V   = KC%Cn_alpha_q_circ * ( 1.0_ReKi - ((1.0_ReKi + 2.0_ReKi*sqrt(KC%fprimeprime) )/3.0_ReKi)**2  )               ! Eqn. 1.50
   else
      KC%C_V   = KC%Cn_alpha_q_circ *  ( 1.0_ReKi - ( 0.5_ReKi + 0.5_ReKi*sqrt(KC%fprimeprime) )**2 )                         ! Eqn. 1.49
   end if

   KC%T_V = BL_p%T_V0 / OtherState%sigma3(i,j)                                                                                ! Eqn 1.48
   
   if (OtherState%FirstPass(i,j)) then
      KC%Cn_v = 0.0_ReKi
   else
      if (xd%tau_V(i,j) > BL_p%T_VL .AND. KC%Kalpha_f * KC%dalpha0 > 0 ) then ! .AND. (.not. LESF)
         ! We no longer require that T_V will always equal T_V0/2 when this condition is satisfied as was the case in AD v13 GJH 7/20/2017
         ! If we fall into this condition, we need to require we stay here until the current vortex is shed (i.e., tauV is reset to zero)
         if ( p%UAMod == UA_Gonzalez ) then   !Added this equation from the formulation used in UAMod=UA_Gonzalez
            KC%Cn_v = xd%Cn_v_minus1(i,j)*exp(-2.0_ReKi*KC%ds/KC%T_V)
         else    
            KC%Cn_v = xd%Cn_v_minus1(i,j)*exp(-KC%ds/KC%T_V)                                                                  ! Eqn 1.52
         end if
      else      
         KC%Cn_v = Get_ExpEqn( KC%ds, KC%T_V, xd%Cn_v_minus1(i,j), KC%C_V, xd%C_V_minus1(i,j) )                               ! Eqn 1.47
      end if
   
      if ( KC%Cn_v < 0.0_ReKi ) then
         KC%Cn_v = 0.0_ReKi
      end if      
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
   
   real(DbKi),                   intent(in   )  :: dt          ! time step length (s)
   type(UA_InitInputType),       intent(inout)  :: InitInp     ! input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UA_ParameterType),       intent(inout)  :: p           ! parameters
   integer(IntKi),               intent(  out)  :: ErrStat     ! error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! error message if ErrStat /= ErrID_None

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
   
   p%c          = InitInp%c         ! this can't be 0
   p%numBlades  = InitInp%numBlades
   p%nNodesPerBlade  = InitInp%nNodesPerBlade
   p%UAMod      = InitInp%UAMod    
   p%a_s        = InitInp%a_s ! this can't be 0
   p%Flookup    = InitInp%Flookup
   p%ShedEffect = InitInp%ShedEffect
   
   if (p%UAMod==UA_HGM) then
      p%lin_nx = p%numBlades*p%nNodesPerBlade*4
   else
      p%lin_nx = 0
   end if
   
end subroutine UA_SetParameters
!==============================================================================   

!==============================================================================
subroutine UA_InitStates_Misc( p, AFInfo, AFIndx, x, xd, OtherState, m, ErrStat, ErrMsg )  
! Called by : UA_Init
! Calls  to : NONE
!..............................................................................

   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(AFI_ParameterType),      intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   integer(IntKi),               intent(in   )  :: AFIndx(:,:)
   type(UA_ContinuousStateType), intent(inout)  :: x           ! Initial continuous states
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
   if (p%UAMod == UA_HGM) then
   
      allocate( x%element( p%nNodesPerBlade, p%numBlades ), stat=ErrStat2 )
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,"Cannot allocate x%x.",ErrStat,ErrMsg,RoutineName)

      allocate( OtherState%n(p%nNodesPerBlade, p%numBlades), stat=ErrStat2)
         if (ErrStat2 /= 0 ) call SetErrStat( ErrID_Fatal, " Error allocating OtherState%n.", ErrStat, ErrMsg, RoutineName)
      
      
   else
      call AllocAry( xd%alpha_minus1,        p%nNodesPerBlade,p%numBlades, 'xd%alpha_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry( xd%alpha_filt_minus1,   p%nNodesPerBlade,p%numBlades, 'xd%alpha_filt_minus1', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
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
      call AllocAry(xd%fprimeprime_m_minus1,p%nNodesPerBlade,p%numBlades,'xd%fprimeprime_m_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%Df_minus1           ,p%nNodesPerBlade,p%numBlades,'xd%Df_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%Df_c_minus1         ,p%nNodesPerBlade,p%numBlades,'xd%Df_c_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%Df_m_minus1         ,p%nNodesPerBlade,p%numBlades,'xd%Df_m_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%Dalphaf_minus1      ,p%nNodesPerBlade,p%numBlades,'xd%Dalphaf_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%alphaf_minus1       ,p%nNodesPerBlade,p%numBlades,'xd%alphaf_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%fprime_minus1       ,p%nNodesPerBlade,p%numBlades,'xd%fprime_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%fprime_c_minus1     ,p%nNodesPerBlade,p%numBlades,'xd%fprime_c_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%fprime_m_minus1     ,p%nNodesPerBlade,p%numBlades,'xd%fprime_m_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%tau_V               ,p%nNodesPerBlade,p%numBlades,'xd%tau_V',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%tau_V_minus1        ,p%nNodesPerBlade,p%numBlades,'xd%tau_V_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%Cn_v_minus1         ,p%nNodesPerBlade,p%numBlades,'xd%Cn_v_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(xd%C_V_minus1          ,p%nNodesPerBlade,p%numBlades,'xd%C_V_minus1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      call AllocAry(OtherState%sigma1   ,p%nNodesPerBlade,p%numBlades,'OtherState%sigma1',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(OtherState%sigma1c   ,p%nNodesPerBlade,p%numBlades,'OtherState%sigma1c',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(OtherState%sigma1m   ,p%nNodesPerBlade,p%numBlades,'OtherState%sigma1m',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(OtherState%sigma3   ,p%nNodesPerBlade,p%numBlades,'OtherState%sigma3',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      
#     ifdef UA_OUTS
         call AllocAry(m%TESF     ,p%nNodesPerBlade,p%numBlades,'m%TESF',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call AllocAry(m%LESF     ,p%nNodesPerBlade,p%numBlades,'m%LESF',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call AllocAry(m%VRTX     ,p%nNodesPerBlade,p%numBlades,'m%VRTX',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call AllocAry(m%T_Sh     ,p%nNodesPerBlade,p%numBlades,'m%T_Sh',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
#     endif
   end if
   
   call AllocAry(m%weight     ,p%nNodesPerBlade,p%numBlades,'m%weight',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   call AllocAry(OtherState%FirstPass,p%nNodesPerBlade,p%numBlades,'OtherState%FirstPass',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(OtherState%UA_off_forGood,p%nNodesPerBlade,p%numBlades,'OtherState%UA_off_forGood',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   if (ErrStat >= AbortErrLev) return

   
   call UA_ReInit( p, AFInfo, AFIndx, x, xd, OtherState, m, ErrStat2,ErrMsg2 )   ! initializes values of states and misc vars
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
end subroutine UA_InitStates_Misc 
!==============================================================================   
subroutine UA_ReInit( p, AFInfo, AFIndx, x, xd, OtherState, m, ErrStat, ErrMsg )  
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(AFI_ParameterType),      intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   integer(IntKi),               intent(in   )  :: AFIndx(:,:)
   type(UA_ContinuousStateType), intent(inout)  :: x           ! Initial continuous states
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Initial other states
   type(UA_MiscVarType),         intent(inout)  :: m           ! Initial misc/optimization variables

   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'UA_ReInit'
   integer                                      :: i
   integer                                      :: j
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   m%FirstWarn_M = .true.
   m%FirstWarn_UA = .true.
   m%FirstWarn_UA_off = .true.
   m%weight = 1.0
   
   OtherState%FirstPass = .true.
   
   OtherState%UA_off_forGood = .false. ! flag that determines if UA parameters are invalid and should be turned off for the whole simulation
   do j=1,size(OtherState%UA_off_forGood,2)
      do i=1,size(OtherState%UA_off_forGood,1)
      
         call UA_TurnOff_param(p, AFInfo(AFIndx(i,j)), ErrStat2, ErrMsg2)
         if (ErrStat2 > ErrID_None) then
            call WrScr( 'Warning: Turning off Unsteady Aerodynamics because '//trim(ErrMsg2)//' (node '//trim(num2lstr(i))//', blade '//trim(num2lstr(j))//')' )
            OtherState%UA_off_forGood(i,j) = .true.
            m%weight(i,j) = 0.0_ReKi
         end if
         
      end do
   end do
   
   
   if ( p%UAMod == UA_HGM ) then
   
      OtherState%n   = -1  ! we haven't updated OtherState%xdot, yet
      
      do j=1,size(x%element,2)
         do i=1,size(x%element,1)
            x%element(i,j)%x = 0.0_ReKi
         end do
      end do

      do i = 1, size(OtherState%xdot)
         call UA_CopyContState( x, OtherState%xdot(i), MESH_UPDATECOPY, ErrStat2, ErrMsg2) ! there are no meshes, so the control code is irrelevant
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end do
   
   else
      OtherState%sigma1    = 1.0_ReKi
      OtherState%sigma1c   = 1.0_ReKi
      OtherState%sigma1m   = 1.0_ReKi
      OtherState%sigma3    = 1.0_ReKi

#     ifdef UA_OUTS
         m%TESF      = .FALSE.
         m%LESF      = .FALSE.
         m%VRTX      = .FALSE.
         m%T_sh      = 0.0_ReKi
#     endif
   
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
      xd%fprimeprime_m_minus1 = 0.0_ReKi
      xd%Df_minus1            = 0.0_ReKi
      xd%Df_c_minus1          = 0.0_ReKi
      xd%Df_m_minus1          = 0.0_ReKi
      xd%Dalphaf_minus1       = 0.0_ReKi
      xd%alphaf_minus1        = 0.0_ReKi
      xd%fprime_minus1        = 0.0_ReKi
      xd%fprime_c_minus1      = 0.0_ReKi
      xd%fprime_m_minus1      = 0.0_ReKi
      xd%tau_V                = 0.0_ReKi 
      xd%tau_V_minus1         = 0.0_ReKi 
      xd%Cn_v_minus1          = 0.0_ReKi
      xd%C_V_minus1           = 0.0_ReKi  ! This probably should not be set to 0.0, but should be set 
   end if

end subroutine UA_ReInit
!==============================================================================
subroutine UA_Init( InitInp, u, p, x, xd, OtherState, y,  m, Interval, &
                    AFInfo, AFIndx, InitOut,ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
! Called by : DRIVER
! Calls  to : NWTC_Init, UA_SetParameters, UA_InitStates
!..............................................................................

   type(UA_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UA_InputType),           intent(in   )  :: u           ! An initial guess for the input; input mesh must be defined
   type(UA_ParameterType),       intent(  out)  :: p           ! Parameters
   type(UA_ContinuousStateType), intent(  out)  :: x           ! Initial continuous states
   type(UA_DiscreteStateType),   intent(  out)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(  out)  :: OtherState  ! Initial other states
   type(UA_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   type(UA_MiscVarType),         intent(  out)  :: m           ! Initial misc/optimization variables
   real(DbKi),                   intent(in   )  :: interval    ! Coupling interval in seconds: the rate that
                                                               !   (1) BEMT_UpdateStates() is called in loose coupling &
                                                               !   (2) BEMT_UpdateDiscState() is called in tight coupling.
                                                               !   Input is the suggested time from the glue code;
                                                               !   Output is the actual coupling interval that will be used
                                                               !   by the glue code.
   type(AFI_ParameterType),      intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   integer(IntKi),               intent(in   )  :: AFIndx(:,:)
   type(UA_InitOutputType),      intent(  out)  :: InitOut     ! Output for initialization routine
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(ErrMsgLen)                         :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                               :: errStat2    ! temporary Error status of the operation
   character(*), parameter                      :: RoutineName = 'UA_Init'
   
#ifdef UA_OUTS
   integer(IntKi)                               :: i,j, iNode, iOffset
   character(64)                                :: chanPrefix
#endif   
   
      ! Initialize variables for this routine
   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

   call UA_ValidateInput(InitInp, AFInfo, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   
      ! Allocate and set parameter data structure using initialization data
   call UA_SetParameters( interval, InitInp, p, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return    
   
      ! initialize the discrete states, other states, and misc variables
   call UA_InitStates_Misc( p, AFInfo, AFIndx, x, xd, OtherState, m, ErrStat2, ErrMsg2 )     ! initialize the continuous states
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return    
      
   
#ifdef UA_OUTS

      ! Allocate and set the InitOut data
   if (p%UAMod == UA_HGM) then
      p%NumOuts = 19
   else
      p%NumOuts = 45
   end if
      
   allocate(InitOut%WriteOutputHdr(p%NumOuts*p%numBlades*p%nNodesPerBlade),STAT=ErrStat2)
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'Error allocating WriteOutputHdr.',ErrStat,ErrMsg,RoutineName)
   allocate(InitOut%WriteOutputUnt(p%NumOuts*p%numBlades*p%nNodesPerBlade),STAT=ErrStat2)
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'Error allocating WriteOutputUnt.',ErrStat,ErrMsg,RoutineName)
   allocate(y%WriteOutput(p%NumOuts*p%numBlades*p%nNodesPerBlade),STAT=ErrStat2)
      if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'Error allocating y%WriteOutput.',ErrStat,ErrMsg,RoutineName)
   if (ErrStat >= AbortErrLev) return
   y%WriteOutput = 0.0_ReKi
         
   iNode = 0
   do j = 1,p%numBlades
      do i = 1,p%nNodesPerBlade
         
         iOffset = (i-1)*p%NumOuts + (j-1)*p%nNodesPerBlade*p%NumOuts 
                  
         chanPrefix = "B"//trim(num2lstr(j))//"N"//trim(num2lstr(i))
         
         InitOut%WriteOutputHdr(iOffset+ 1)  = trim(chanPrefix)//'ALPHA'
         InitOut%WriteOutputHdr(iOffset+ 2)  = trim(chanPrefix)//'VREL'
         InitOut%WriteOutputHdr(iOffset+ 3)  = trim(chanPrefix)//'Cn'
         InitOut%WriteOutputHdr(iOffset+ 4)  = trim(chanPrefix)//'Cc'
         InitOut%WriteOutputHdr(iOffset+ 5)  = trim(chanPrefix)//'Cl'
         InitOut%WriteOutputHdr(iOffset+ 6)  = trim(chanPrefix)//'Cd'
         InitOut%WriteOutputHdr(iOffset+ 7)  = trim(chanPrefix)//'Cm'
         
         InitOut%WriteOutputUnt(iOffset+ 1)  ='(deg)'
         InitOut%WriteOutputUnt(iOffset+ 2)  ='(m/s)'
         InitOut%WriteOutputUnt(iOffset+ 3)  ='(-)'
         InitOut%WriteOutputUnt(iOffset+ 4)  ='(-)'
         InitOut%WriteOutputUnt(iOffset+ 5)  ='(-)'
         InitOut%WriteOutputUnt(iOffset+ 6)  ='(-)'
         InitOut%WriteOutputUnt(iOffset+ 7)  ='(-)'
         
         if (p%UAmod == UA_HGM) then
            
            InitOut%WriteOutputHdr(iOffset+ 8)  = trim(chanPrefix)//'omega'
            InitOut%WriteOutputHdr(iOffset+ 9)  = trim(chanPrefix)//'alphaE'
            InitOut%WriteOutputHdr(iOffset+10)  = trim(chanPrefix)//'Tu'
            InitOut%WriteOutputHdr(iOffset+11)  = trim(chanPrefix)//'alpha_34'
            InitOut%WriteOutputHdr(iOffset+12)  = trim(chanPrefix)//'cl_fs'
            InitOut%WriteOutputHdr(iOffset+13)  = trim(chanPrefix)//'fs_aE'
         
            InitOut%WriteOutputHdr(iOffset+14)  = trim(chanPrefix)//'x1'
            InitOut%WriteOutputHdr(iOffset+15)  = trim(chanPrefix)//'x2'
            InitOut%WriteOutputHdr(iOffset+16)  = trim(chanPrefix)//'x3'
            InitOut%WriteOutputHdr(iOffset+17)  = trim(chanPrefix)//'x4'
            InitOut%WriteOutputHdr(iOffset+18)  = trim(chanPrefix)//'k'
            InitOut%WriteOutputHdr(iOffset+19)  = trim(chanPrefix)//'weight'
            

            InitOut%WriteOutputUnt(iOffset+ 8)  = '(deg/sec)'
            InitOut%WriteOutputUnt(iOffset+ 9)  = '(deg)'
            InitOut%WriteOutputUnt(iOffset+10)  = '(s)'
            InitOut%WriteOutputUnt(iOffset+11)  = '(deg)'
            InitOut%WriteOutputUnt(iOffset+12)  = '(-)'
            InitOut%WriteOutputUnt(iOffset+13)  = '(-)'
         
            InitOut%WriteOutputUnt(iOffset+14)  = '(rad)'
            InitOut%WriteOutputUnt(iOffset+15)  = '(rad)'
            InitOut%WriteOutputUnt(iOffset+16)  = '(-)'
            InitOut%WriteOutputUnt(iOffset+17)  = '(-)'
            InitOut%WriteOutputUnt(iOffset+18)  = '(-)'
            InitOut%WriteOutputUnt(iOffset+19)  = '(-)'

         else
            InitOut%WriteOutputHdr(iOffset+ 8)  = trim(chanPrefix)//'Cn_aq_circ'
            InitOut%WriteOutputHdr(iOffset+ 9)  = trim(chanPrefix)//'Cn_aq_nc'
            InitOut%WriteOutputHdr(iOffset+10)  = trim(chanPrefix)//'Cn_pot'
            InitOut%WriteOutputHdr(iOffset+11)  = trim(chanPrefix)//'Dp'
            InitOut%WriteOutputHdr(iOffset+12)  = trim(chanPrefix)//'Cn_prime'
            InitOut%WriteOutputHdr(iOffset+13)  = trim(chanPrefix)//'fprime'
            InitOut%WriteOutputHdr(iOffset+14)  = trim(chanPrefix)//'Df'
            InitOut%WriteOutputHdr(iOffset+15)  = trim(chanPrefix)//'Cn_v'
            InitOut%WriteOutputHdr(iOffset+16)  = trim(chanPrefix)//'Tau_V'
            InitOut%WriteOutputHdr(iOffset+17)  = trim(chanPrefix)//'LESF'
            InitOut%WriteOutputHdr(iOffset+18)  = trim(chanPrefix)//'TESF' 
            InitOut%WriteOutputHdr(iOffset+19)  = trim(chanPrefix)//'VRTX'
            InitOut%WriteOutputHdr(iOffset+20)  = trim(chanPrefix)//'C_v'
            InitOut%WriteOutputHdr(iOffset+21)  = trim(chanPrefix)//'Cm_a_nc'
            InitOut%WriteOutputHdr(iOffset+22)  = trim(chanPrefix)//'Cm_q_nc'
            InitOut%WriteOutputHdr(iOffset+23)  = trim(chanPrefix)//'Cm_v'
            InitOut%WriteOutputHdr(iOffset+24)  = trim(chanPrefix)//'alpha_p_f'
            InitOut%WriteOutputHdr(iOffset+25)  = trim(chanPrefix)//'Dalphaf'
            InitOut%WriteOutputHdr(iOffset+26)  = trim(chanPrefix)//'PMC'
            InitOut%WriteOutputHdr(iOffset+27)  = trim(chanPrefix)//'T_f'
            InitOut%WriteOutputHdr(iOffset+28)  = trim(chanPrefix)//'T_V'
            InitOut%WriteOutputHdr(iOffset+29)  = trim(chanPrefix)//'dS'
            InitOut%WriteOutputHdr(iOffset+30)  = trim(chanPrefix)//'T_alpha'
            InitOut%WriteOutputHdr(iOffset+31)  = trim(chanPrefix)//'T_q'
            InitOut%WriteOutputHdr(iOffset+32)  = trim(chanPrefix)//'k_alpha'
            InitOut%WriteOutputHdr(iOffset+33)  = trim(chanPrefix)//'k_q'
            InitOut%WriteOutputHdr(iOffset+34)  = trim(chanPrefix)//'alpha_e'
            InitOut%WriteOutputHdr(iOffset+35)  = trim(chanPrefix)//'X1'
            InitOut%WriteOutputHdr(iOffset+36)  = trim(chanPrefix)//'X2'
            InitOut%WriteOutputHdr(iOffset+37)  = trim(chanPrefix)//'cn_q_nc'
            InitOut%WriteOutputHdr(iOffset+38)  = trim(chanPrefix)//'alpha_f'
            InitOut%WriteOutputHdr(iOffset+39)  = trim(chanPrefix)//'fprimeprime'
            InitOut%WriteOutputHdr(iOffset+40)  = trim(chanPrefix)//'sigma1'
            InitOut%WriteOutputHdr(iOffset+41)  = trim(chanPrefix)//'sigma3'
            InitOut%WriteOutputHdr(iOffset+42)  = trim(chanPrefix)//'T_sh'
            InitOut%WriteOutputHdr(iOffset+43)  = trim(chanPrefix)//'k'
            InitOut%WriteOutputHdr(iOffset+44)  = trim(chanPrefix)//'ALPHA_filt'
            InitOut%WriteOutputHdr(iOffset+44)  = trim(chanPrefix)//'weight'

            
            InitOut%WriteOutputUnt(iOffset+8)  ='(-)'
            InitOut%WriteOutputUnt(iOffset+9)  ='(-)'
            InitOut%WriteOutputUnt(iOffset+10) ='(-)'
            InitOut%WriteOutputUnt(iOffset+11) ='(-)'
            InitOut%WriteOutputUnt(iOffset+12) ='(-)'
            InitOut%WriteOutputUnt(iOffset+13) ='(-)'
            InitOut%WriteOutputUnt(iOffset+14) ='(-)'
            InitOut%WriteOutputUnt(iOffset+15) ='(-)'
            InitOut%WriteOutputUnt(iOffset+16) ='(-)'
            InitOut%WriteOutputUnt(iOffset+17) ='(-)'
            InitOut%WriteOutputUnt(iOffset+18) ='(-)'
            InitOut%WriteOutputUnt(iOffset+19) ='(-)'
            InitOut%WriteOutputUnt(iOffset+20) ='(-)'
            InitOut%WriteOutputUnt(iOffset+21) ='(-)'
            InitOut%WriteOutputUnt(iOffset+22) ='(-)'
            InitOut%WriteOutputUnt(iOffset+23) ='(-)'
            InitOut%WriteOutputUnt(iOffset+24) ='(rad)'
            InitOut%WriteOutputUnt(iOffset+25) ='(-)'
            InitOut%WriteOutputUnt(iOffset+26) ='(-)'
            InitOut%WriteOutputUnt(iOffset+27) ='(-)'
            InitOut%WriteOutputUnt(iOffset+28) ='(-)'
            InitOut%WriteOutputUnt(iOffset+29) ='(-)'
            InitOut%WriteOutputUnt(iOffset+30) ='(-)'
            InitOut%WriteOutputUnt(iOffset+31) ='(-)'
            InitOut%WriteOutputUnt(iOffset+32) ='(-)'
            InitOut%WriteOutputUnt(iOffset+33) ='(-)'
            InitOut%WriteOutputUnt(iOffset+34) ='(rad)'
            InitOut%WriteOutputUnt(iOffset+35) ='(-)'
            InitOut%WriteOutputUnt(iOffset+36) ='(-)'
            InitOut%WriteOutputUnt(iOffset+37) ='(-)'
            InitOut%WriteOutputUnt(iOffset+38) ='(rad)'
            InitOut%WriteOutputUnt(iOffset+39) ='(-)'
            InitOut%WriteOutputUnt(iOffset+40) ='(-)'
            InitOut%WriteOutputUnt(iOffset+41) ='(-)'
            InitOut%WriteOutputUnt(iOffset+42) ='(-)'
            InitOut%WriteOutputUnt(iOffset+43) ='(-)'
            InitOut%WriteOutputUnt(iOffset+44)  ='(deg)'
            InitOut%WriteOutputUnt(iOffset+45) ='(-)'
         end if
         
      end do
   end do
   
   p%OutSFmt = 'A19'
   p%OutFmt  = 'ES19.5e3'
   p%Delim   =''
   
   if (p%NumOuts > 0) then
      CALL GetNewUnit( p%unOutFile, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

      CALL OpenFOutFile ( p%unOutFile, trim(InitInp%OutRootName)//'.out', ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return

         ! Heading:
      WRITE (p%unOutFile,'(/,A/)')  'These predictions were generated by UnsteadyAero on '//CurDate()//' at '//CurTime()//'.'
      WRITE (p%unOutFile,'(/,A/)', IOSTAT=ErrStat)  ' '
         
      ! Write the names of the output parameters:
      WRITE (p%unOutFile,'('//p%OutSFmt//')', ADVANCE='no' ) 'Time'
      do i=1,size(InitOut%WriteOutputHdr)
         WRITE (p%unOutFile,'(:,A,'//trim( p%OutSFmt )//')', ADVANCE='no' )  p%Delim, trim(InitOut%WriteOutputHdr(i))
      end do  
      WRITE (p%unOutFile,'()', IOSTAT=ErrStat2)          ! write the line return
   
      WRITE (p%unOutFile,'('//p%OutSFmt//')', ADVANCE='no' ) '(s)'
      do i=1,size(InitOut%WriteOutputUnt)
         WRITE (p%unOutFile,'(:,A,'//trim( p%OutSFmt )//')', ADVANCE='no' )  p%Delim, trim(InitOut%WriteOutputUnt(i))
      end do
      WRITE (p%unOutFile,'()', IOSTAT=ErrStat2)          ! write the line return
   end if
         

#else
   p%NumOuts = 0
   p%unOutFile = -1
   !.....................................
   ! add the following two lines only to avoid compiler warnings about uninitialized variables when not building the UA driver:
   y%cm = 0.0_ReKi 
   InitOut%Version = ProgDesc( 'Unsteady Aero', '', '' )
   !.....................................
   
#endif
   
end subroutine UA_Init
!==============================================================================     
subroutine UA_ValidateInput(InitInp, AFInfo, ErrStat, ErrMsg)
   type(UA_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(AFI_ParameterType),      intent(in   )  :: AFInfo(:)   !< The airfoil parameter data
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                               :: ErrStat2     ! Error status of the operation
   character(ErrMsgLen)                         :: ErrMsg2      ! Error message if ErrStat /= ErrID_None
   character(*), parameter                      :: RoutineName = 'UA_ValidateInput'

   integer(IntKi)                               :: i             ! loop counter
   
   ErrStat = ErrID_None
   ErrMsg  = ""

   !>>> remove after this feature gets tested better:
   if (InitInp%UAMod == UA_HGM ) then
      call SetErrStat( ErrID_Warn, "UAMod 4 (continuous HGM model) is in beta for this version of OpenFAST.", ErrStat, ErrMsg, RoutineName )
   end if
   !<<<

   if (InitInp%UAMod < UA_Gonzalez .or. InitInp%UAMod > UA_HGM ) call SetErrStat( ErrID_Fatal, &
      "In this version, UAMod must be 2 (Gonzalez's variant), 3 (Minnema/Pierce variant), or 4 (continuous HGM model).", ErrStat, ErrMsg, RoutineName )  ! NOTE: for later-  1 (baseline/original) 
      
   if (.not. InitInp%FLookUp ) call SetErrStat( ErrID_Fatal, 'FLookUp must be TRUE for this version.', ErrStat, ErrMsg, RoutineName )
   
   if (InitInp%a_s <= 0.0) call SetErrStat ( ErrID_Fatal, 'The speed of sound (SpdSound) must be greater than zero.', ErrStat, ErrMsg, RoutineName )

   ! check that the airfoils have appropriate data for UA
   do i=1,size(AFInfo,1)
      call UA_ValidateAFI(InitInp%UAMod, AFInfo(i), ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end do

   
   
end subroutine UA_ValidateInput
!==============================================================================     
subroutine UA_ValidateAFI(UAMod, AFInfo, ErrStat, ErrMsg)
   integer(IntKi),               intent(in   )  :: UAMod       ! which UA model we are using
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                               :: j
   character(*), parameter                      :: RoutineName = 'UA_ValidateAFI'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if (.not. allocated(AFInfo%Table)) then
      call SetErrStat(ErrID_Fatal, 'Airfoil table not allocated in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
   else

      do j=1, AFInfo%NumTabs
         if ( .not. AFInfo%Table(j)%InclUAdata ) then
            call SetErrStat(ErrID_Fatal, 'Airfoil file "'//trim(AFInfo%FileName)//'", table #'//trim(num2lstr(j))// &
                            ' does not contain parameters for UA data.', ErrStat, ErrMsg, RoutineName )
         else
            ! parameters used only for UAMod/=UA_HGM)
            if (UAMod /= UA_HGM) then
               if ( EqualRealNos(AFInfo%Table(j)%UA_BL%St_sh, 0.0_ReKi) ) then
                  call SetErrStat(ErrID_Fatal, 'UA St_sh parameter must not be 0 in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
               end if

               if ( AFInfo%Table(j)%UA_BL%alpha1 > pi .or. AFInfo%Table(j)%UA_BL%alpha1 < -pi ) then
                  call SetErrStat(ErrID_Fatal, 'UA alpha1 parameter must be between -180 and 180 degrees in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
               end if

               if ( AFInfo%Table(j)%UA_BL%alpha2 > pi .or. AFInfo%Table(j)%UA_BL%alpha2 < -pi ) then
                  call SetErrStat(ErrID_Fatal, 'UA alpha2 parameter must be between -180 and 180 degrees in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
               end if
               
               if ( AFInfo%Table(j)%UA_BL%filtCutOff < 0.0_ReKi ) then
                  call SetErrStat(ErrID_Fatal, 'UA filtCutOff parameter must be greater than 0 in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
               end if

            end if

               ! variables used in all UA models:
            if ( AFInfo%Table(j)%UA_BL%alpha0 > pi .or. AFInfo%Table(j)%UA_BL%alpha0 < -pi ) then
               call SetErrStat(ErrID_Fatal, 'UA alpha0 parameter must be between -180 and 180 degrees in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
            end if

            if ( AFInfo%Table(j)%UA_BL%T_f0 < 0.0_ReKi ) then
               call SetErrStat(ErrID_Fatal, 'UA T_f0 parameter must be greater than 0 in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
            end if
            
            if ( AFInfo%Table(j)%UA_BL%T_p < 0.0_ReKi ) then
               call SetErrStat(ErrID_Fatal, 'UA T_p parameter must be greater than 0 in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
            end if
            
         end if
            
         if ( AFInfo%Table(j)%UA_BL%UACutout < 0.0_ReKi ) then
            call SetErrStat(ErrID_Fatal, 'UA UACutout parameter must not be negative in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
         end if
            
            

      end do

      if (ErrStat >= AbortErrLev) return

      if (UAMod /= UA_HGM) then
         ! check interpolated values:
      
         do j=2, AFInfo%NumTabs
            if ( sign( 1.0_ReKi, AFInfo%Table(j)%UA_BL%St_sh) /= &
                 sign( 1.0_ReKi, AFInfo%Table(1)%UA_BL%St_sh) ) then
               call SetErrStat(ErrID_Fatal, 'UA St_sh parameter (interpolated value) must not be 0 in "'//trim(AFInfo%FileName)//'".', ErrStat, ErrMsg, RoutineName )
               exit
            end if
         end do
      end if
      
   end if
   

end subroutine UA_ValidateAFI
!==============================================================================
!> This routine checks if the UA parameters indicate that UA should not be used. (i.e., if C_nalpha = 0)
!! This should be called at initialization.
subroutine UA_TurnOff_param(p, AFInfo, ErrStat, ErrMsg)
   type(UA_ParameterType),       intent(in   )  :: p           ! The UA parameter data
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                               :: j

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (p%UAMod == UA_HGM) then
         ! unsteady aerodynamics will be turned off
      do j=1, AFInfo%NumTabs
         if ( EqualRealNos(AFInfo%Table(j)%UA_BL%C_lalpha, 0.0_ReKi) ) then
            ErrStat = ErrID_Fatal
            ErrMsg  = 'C_lalpha is 0.'
            return
         end if
      end do
      
         ! now check about interpolated values:
      do j=2, AFInfo%NumTabs
         if ( sign( 1.0_ReKi, AFInfo%Table(j)%UA_BL%C_lalpha) /= &
              sign( 1.0_ReKi, AFInfo%Table(1)%UA_BL%C_lalpha) ) then
            ErrStat = ErrID_Fatal
            ErrMsg  = 'C_lalpha (interpolated value) could be 0.'
            return
         end if
      end do
      
   else
         ! unsteady aerodynamics will be turned off
      do j=1, AFInfo%NumTabs
         if ( EqualRealNos(AFInfo%Table(j)%UA_BL%C_nalpha, 0.0_ReKi) ) then
            ErrStat = ErrID_Fatal
            ErrMsg  = 'C_nalpha is 0.'
            return
         end if
      end do

         ! now check about interpolated values:
      do j=2, AFInfo%NumTabs
         if ( sign( 1.0_ReKi, AFInfo%Table(j)%UA_BL%C_nalpha) /= &
              sign( 1.0_ReKi, AFInfo%Table(1)%UA_BL%C_nalpha) ) then
            ErrStat = ErrID_Fatal
            ErrMsg  = 'C_nalpha (interpolated value) could be 0.'
            return
         end if
      end do
      
   end if
   

end subroutine UA_TurnOff_param
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
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local Variables
  
   LOGICAL                                     :: LESF      ! logical flag indicating if leading edge separation is possible [-]
   LOGICAL                                     :: VRTX      ! logical flag indicating if a vortex is being processed [-]
   LOGICAL                                     :: TESF      ! logical flag indicating if trailing edge separation is possible [-]
      
   real(ReKi)                                   :: Kafactor
   real(ReKi)                                   :: T_sh
      
   type(AFI_UA_BL_Type)                         :: BL_p  ! airfoil UA parameters retrieved in Kelvin Chain 
   type(UA_KelvinChainType)                     :: KC    ! values computed in Kelvin Chain
   
   
   character(ErrMsgLen)                         :: errMsg2
   integer(IntKi)                               :: errStat2
   character(*), parameter                      :: RoutineName = 'UA_UpdateDiscOtherState'
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               

!bjj: use sigma1 and sigma3 states values from t      
   
      call ComputeKelvinChain(i, j, u, p, xd, OtherState, m, AFInfo, KC, BL_p, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat>= AbortErrLev) return
            
      
      T_sh = 2.0_ReKi*(1.0_ReKi-KC%fprimeprime) / BL_p%St_sh                      ! Eq 1.54
      
      !---------------------------------------------------------
      ! Update the OtherStates
      !---------------------------------------------------------
   
      LESF = (KC%Cn_prime > BL_p%Cn1) .or. ( KC%Cn_prime < BL_p%Cn2 ) ! LE separation can occur when this is .true.; assumption is that Cn2 <  0.0 and Cn1 > 0
      if (p%UAMod == UA_Gonzalez) then   !Added specific logic for UAMod=UA_Gonzalez
         TESF = ABS( KC%Cn_prime ) > ABS( xd%Cn_prime_minus1(i,j) )
      else
         TESF = KC%fprimeprime < xd%fprimeprime_minus1(i,j) ! Separation point is moving towards the Leading Edge when .true.; otherwise separation point is moving toward trailing edge   
      end if
   
       
      VRTX =  (xd%tau_V(i,j) <= BL_p%T_VL) .and. (xd%tau_V(i,j) > 0.0_ReKi) ! Is the vortex over the chord? 
      
      
   !---------------------------------------------------------
   ! Update the OtherStates
   !---------------------------------------------------------

!bjj: update sigma1 to value at t + dt 
   
      ! Set sigma1 to the default value.  T_f = T_f0 / sigma1
   OtherState%sigma1( i,j) = 1.0_ReKi
   OtherState%sigma1c(i,j) = 1.0_ReKi
   OtherState%sigma1m(i,j) = 1.0_ReKi

   Kafactor      = KC%Kalpha_f*KC%dalpha0  ! indicates if the airfoil is moving towards alpha0 [ Kafactor < 0.0 ] or away [ Kafactor > 0.0]
   
if (p%UAMod == UA_Gonzalez) then   !Added modifiers for Tfn, Tfc and Tfm depending on the aerodynamic state, this set of values can be optimized for different airfoils or wind conditions
      if ( TESF .AND. .NOT.LESF .AND. .NOT.VRTX ) then
         if ( KC%fprimeprime>0.7 ) then
            OtherState%sigma1( i,j)=1.0!0.2
            OtherState%sigma1c(i,j)=1.0!0.2
            OtherState%sigma1m(i,j)=1.0!0.1
         else !if ((KC%fprimeprime<=0.7).AND.(TESF).AND.(.NOT.LESF).AND. NOT.VRTX ) then
            OtherState%sigma1( i,j)=1.0!0.5
            OtherState%sigma1c(i,j)=1.0!0.8
            OtherState%sigma1m(i,j)=1.0!0.3
         endif
      elseif ((xd%tau_V(i,j) < BL_p%T_VL).AND.(xd%tau_V(i,j)>0.001).AND.(TESF)) then
         OtherState%sigma1( i,j)=2.0!4.0
         OtherState%sigma1c(i,j)=1.0!1.0
         OtherState%sigma1m(i,j)=1.0!0.2
      elseif ( LESF ) then
         if (TESF) then
            OtherState%sigma1( i,j)=2.0!4.0
            OtherState%sigma1c(i,j)=1.0!1.0
            OtherState%sigma1m(i,j)=1.0!0.3
         else !if ((LESF).AND.(.NOT.TESF)) then
            OtherState%sigma1( i,j)=1.0!0.2
            OtherState%sigma1c(i,j)=1.0!0.2
            OtherState%sigma1m(i,j)=1.0!0.2
         endif
      elseif ( .NOT.TESF ) then
         if (KC%fprimeprime<=0.7) then
            OtherState%sigma1( i,j)=0.5!0.5
            OtherState%sigma1c(i,j)=1.0!0.5
            OtherState%sigma1m(i,j)=1.0!2.0
          else !if ((KC%fprimeprime>0.7).AND.(.NOT.TESF)) then
            OtherState%sigma1( i,j)=0.5!4.0
            OtherState%sigma1c(i,j)=1.0!0.4
            OtherState%sigma1m(i,j)=1.0!2.0
         endif
      else
         OtherState%sigma1( i,j)=1.0_ReKi
         OtherState%sigma1c(i,j)=1.0_ReKi
         OtherState%sigma1m(i,j)=1.0_ReKi
      endif
else
   if ( TESF ) then  ! Flow is continuing or starting to separate
      if (Kafactor < 0.0_ReKi) then
            ! We are moving towards alpha0
         OtherState%sigma1(i,j) = 2.0_ReKi  ! This must be the first check
      else if (.not. LESF ) then
         OtherState%sigma1(i,j) = 1.0_ReKi                  ! Leading edge separation has not occurred and we are moving away from alpha0
      else if (xd%fprimeprime_minus1(i,j) <= 0.7_ReKi) then ! For this else, LESF = True and we are moving away from alpha0
         OtherState%sigma1(i,j) = 2.0_ReKi 
      else
         OtherState%sigma1(i,j) = 1.75_ReKi
      end if
   else ! Reattaching flow 
        
      if (.not. LESF ) then
         OtherState%sigma1(i,j) = 0.5_ReKi
      end if

      if ( VRTX ) then  ! Still shedding a vortex, i.e., the current vortex is still over the chord?
         OtherState%sigma1(i,j) = 0.25_ReKi
      end if
      
      if (Kafactor > 0.0_ReKi) then
         OtherState%sigma1(i,j) = 0.75_ReKi
      end if
            
   end if

end if
      
      

!bjj: update sigma3 to value at t + dt   
   
      ! This is the default value for sigma3 which effects T_V = T_V0 / sigma3
   OtherState%sigma3(i,j) = 1.0_ReKi
   
if (p%UAMod /= UA_Gonzalez) then  !this is not applied for UAMod=UZ_Gonzalez, Tv has always the same value
      ! Identify where the vortex is located relative to the chord 
   
      ! 1) Is the vortex past the trailing edge, but less than 2 chords?
   if ( (xd%tau_V(i,j) <= 2.0_ReKi*BL_p%T_VL) .and. (xd%tau_V(i,j) >= BL_p%T_VL) ) then
         ! We want to diminish the effects of this vortex on Cn by setting a high value of sigma3
      OtherState%sigma3(i,j) = 3.0_ReKi
      if (.not. TESF) then
            ! If we are reattaching the flow, then we when to diminish the current vortex's effects on Cn further

         OtherState%sigma3(i,j) = 4.0_ReKi
         
      end if
      
      ! 2) Is the vortex over the chord
   else if ( VRTX ) then
      if (Kafactor < 0.0_ReKi) then
            ! If we are moving towards alpha0, then we want to reduce the contribution of the vortex to Cn
         OtherState%sigma3(i,j) = 2.0_ReKi
      else
         OtherState%sigma3(i,j) = 1.0_ReKi
      end if
   
      ! 3) Is the vortex ((past 2 chords or at the leading edge), and is the airfoil moving away from the stall region (towards alpha0).
      ! NOTE: We could also end up here if tau_V = 0.0
   else if (Kafactor < 0 ) then 
         ! In this case, we want to diminish the effects of this vortex on Cn by setting a high value of sigma3
      OtherState%sigma3(i,j) = 4.0_ReKi
   end if
      

     
      ! Finally, we will override all the previous values of sigma1 if we are reattaching flow and the rate of change of the of the angle of attack is slowing down
      ! In this case we want to enhance the contribute of the vortex and set sigma3 = 1.0
   if ((.not. TESF) .and. (KC%Kq_f*KC%dalpha0 < 0.0_ReKi)) then
      OtherState%sigma3(i,j) = 1.0_ReKi
   end if
endif  
   
      
#ifdef TEST_UA_SIGMA
   OtherState%sigma3(i,j) = 1.0_ReKi
   OtherState%sigma1(i,j) = 1.0_ReKi
#endif
   
      !---------------------------------------------------------
      ! Update the Discrete States, xd
      !---------------------------------------------------------
     
           
      xd%alpha_minus1(i,j)        = u%alpha
      xd%alpha_filt_minus1(i,j)   = KC%alpha_filt_cur
      xd%Kalpha_f_minus1(i,j)     = KC%Kalpha_f
      xd%Kq_f_minus1(i,j)         = KC%Kq_f     
      xd%q_f_minus1(i,j)          = KC%q_f_cur    
      xd%q_minus1(i,j)            = KC%q_cur
      xd%X1_minus1(i,j)           = KC%X1
      xd%X2_minus1(i,j)           = KC%X2
      xd%X3_minus1(i,j)           = KC%X3
      xd%X4_minus1(i,j)           = KC%X4
      xd%Kprime_alpha_minus1(i,j) = KC%Kprime_alpha
      xd%Kprime_q_minus1(i,j)     = KC%Kprime_q
      xd%Dp_minus1(i,j)           = KC%Dp
      xd%Cn_pot_minus1(i,j)       = KC%Cn_pot
      xd%Cn_prime_minus1(i,j)     = KC%Cn_prime
      xd%fprimeprime_minus1(i,j)  = KC%fprimeprime
      xd%Df_minus1(i,j)           = KC%Df
      if (p%Flookup) then
         xd%Df_c_minus1(i,j)           = KC%Df_c
         xd%Df_m_minus1(i,j)           = KC%Df_m
         xd%fprimeprime_c_minus1(i,j)  = KC%fprimeprime_c
         xd%fprimeprime_m_minus1(i,j)  = KC%fprimeprime_m
         xd%fprime_c_minus1(i,j)       = KC%fprime_c
         xd%fprime_m_minus1(i,j)       = KC%fprime_m
      end if
      
      xd%Dalphaf_minus1(i,j)      = KC%Dalphaf
      xd%fprime_minus1(i,j)       = KC%fprime
      xd%alphaf_minus1(i,j)       = KC%alpha_f
      xd%Cn_v_minus1(i,j)         = KC%Cn_v
      xd%C_V_minus1(i,j)          = KC%C_V
 
         ! If we are currently tracking a vortex, or we are in the stall region, increment tau_V
      if (p%UAMod == UA_Gonzalez) then    !Added specific logic for UAMod=UA_Gonzalez
         if ( (.NOT.LESF .AND. .NOT.VRTX) .OR. & 
              (.NOT.TESF .AND. xd%tau_V(i,j)<0.0001_ReKi) .OR.& 
              (.NOT.TESF .AND. xd%tau_V(i,j) + KC%ds > 2.*BL_p%T_VL) ) then 
            xd%tau_V(i,j)=0.0 
         else
            xd%tau_V(i,j) = xd%tau_V(i,j) + KC%ds
            if (( xd%tau_V(i,j) >= (BL_p%T_VL + T_sh) ) .and. TESF)  then !.and. (TESF) RRD added  
               xd%tau_V(i,j) = xd%tau_V(i,j)-(BL_p%T_VL + T_sh) 
               ! LESF=.FALSE. !also added this 
            end if 
         endif 
         
      else
         if ( xd%tau_V(i,j) > 0.0 .or. LESF ) then
            xd%tau_V(i,j) = xd%tau_V(i,j) + KC%ds     ! Eqn 1.51
         end if

            ! If we a have been tracking a vortex and 1) it is now past the chord [T_VL] and 2) we have gone beyond the next shedding period [T_sh], and 
            !   3) we are continuing the flow serparation, we will shed the existing vortex so that we can create a new one at the leading edge
         if (( xd%tau_V(i,j) >= (BL_p%T_VL + T_sh) ) .and. TESF)  then !.and. (TESF) RRD added
            xd%tau_V(i,j) = 0.0_ReKi
         end if
      end if
      
#ifdef UA_OUTS
   m%TESF(i,j) = TESF  
   m%LESF(i,j) = LESF   
   m%VRTX(i,j) = VRTX 
   m%T_sh(i,j) = T_sh
#endif
      
      
end subroutine UA_UpdateDiscOtherState
!==============================================================================   

!============================================================================== 
subroutine UA_UpdateStates( i, j, t, n, u, uTimes, p, x, xd, OtherState, AFInfo, m, ErrStat, ErrMsg )
! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
! Continuous, constraint, discrete, and other states are updated to values at t + Interval.
!..............................................................................

   integer(IntKi),                intent(in   ) :: i               ! node index within a blade
   integer(IntKi),                intent(in   ) :: j               ! blade index 
   REAL(DbKi),                    INTENT(IN   ) :: t               ! Current simulation time in seconds
   integer(IntKi),                intent(in   ) :: n               !< Current simulation time step n = 0,1,...
   type(UA_InputType),            intent(in   ) :: u(:)            ! Input at current timestep, t and t+dt
   real(DbKi),                    intent(in   ) :: utimes(:)       !< Times associated with u(:), in seconds
   type(UA_ParameterType),        intent(in   ) :: p               ! Parameters
   type(UA_ContinuousStateType),  intent(inout) :: x               ! Input: Continuous states at t;
                                                                   !   Output: Continuous states at t + Interval
   type(UA_DiscreteStateType),    intent(inout) :: xd              ! Input: Discrete states at t;
                                                                   !   Output: Discrete states at t + Interval
   type(UA_OtherStateType),       intent(inout) :: OtherState      ! Input: Other states at t;
                                                                   !   Output: Other states at t + Interval
   type(UA_MiscVarType),          intent(inout) :: m               ! Misc/optimization variables
   type(AFI_ParameterType),       intent(in   ) :: AFInfo          ! The airfoil parameter data
   integer(IntKi),                intent(  out) :: ErrStat         ! Error status of the operation
   character(*),                  intent(  out) :: ErrMsg          ! Error message if ErrStat /= ErrID_None

      ! Local variables  
      
   character(ErrMsgLen)                         :: errMsg2
   integer(IntKi)                               :: errStat2
   character(*), parameter                      :: RoutineName = 'UA_UpdateStates'
   type(UA_InputType)                           :: u_interp_raw    ! Input at current timestep, t and t+dt
   type(UA_InputType)                           :: u_interp        ! Input at current timestep, t and t+dt

      ! Initialize variables

   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
         
   !BJJ: u%u == 0 seems to be the root cause of all sorts of numerical problems....

      
   if (OtherState%UA_off_forGood(i,j)) return   ! we don't have any states to update here
   
   if (p%UAMod == UA_HGM) then
      
         ! initialize states to steady-state values:
      if (OtherState%FirstPass(i,j)) then
         CALL UA_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         call HGM_Steady( i, j, u_interp, p, x%element(i,j), AFInfo, ErrStat2, ErrMsg2 )
      end if

      
      call UA_ABM4( i, j, t, n, u, utimes, p, x, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      x%element(i,j)%x(4) = max( min( x%element(i,j)%x(4), 1.0_R8Ki ), 0.0_R8Ki )
      
         ! these are angles that should not get too large, so I am fixing them here (should turn off UA if this exceeds reasonable numbers)
      if (abs(x%element(i,j)%x(1)) > pi .or. abs(x%element(i,j)%x(2)) > pi) then
         if (m%FirstWarn_UA) then
            call SetErrStat(ErrID_Severe, "Divergent states in UA HGM model", ErrStat, ErrMsg, RoutineName )
            m%FirstWarn_UA = .false.
         end if
         
         call Mpi2pi(x%element(i,j)%x(1))
         call Mpi2pi(x%element(i,j)%x(2))
      end if
      
   else
      if (n<=0) return ! previous logic (before adding UA_HGM required n > 0 before UA_UpdateStates was called)
      
      CALL UA_Input_ExtrapInterp( u, utimes, u_interp_raw, t, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! make sure that u%u is not zero (this previously turned off UA for the entire simulation. 
         ! Now, we keep it on, but we don't want the math to blow up when we divide by u%u)
      call UA_fixInputs(u_interp_raw, u_interp, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   

         ! Update discrete states:
#     ifdef DEBUG_v14
         call UA_UpdateDiscOtherState2( i, j, u_interp, p, xd, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
#     else
         call UA_UpdateDiscOtherState( i, j, u_interp, p, xd, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
#     endif
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   end if

   OtherState%FirstPass(i,j)   = .false.
   
end subroutine UA_UpdateStates
!==============================================================================
!!----------------------------------------------------------------------------------------------------------------------------------
!> routine to initialize the states based on inputs at t=0
!! used to obtain initial values in linearization so that they don't change with each call to calcOutput (or other routines)
subroutine UA_InitStates_AllNodes( u, p, x, OtherState, AFInfo, AFIndx )
   type(UA_InputType),           intent(in   ) :: u(:,:)     !< Inputs at t
   type(UA_ParameterType),       intent(in   ) :: p          !< Parameters
   type(UA_ContinuousStateType), intent(inout) :: x          !< Input: Continuous states at t;
   type(UA_OtherStateType),      intent(inout) :: OtherState !< Other/logical states at t on input; at t+dt on output
   type(AFI_ParameterType),      intent(in   ) :: AFInfo(:)  !< The airfoil parameter data
   INTEGER(IntKi),               intent(in)    :: AFIndx(:,:)
   
   INTEGER(IntKi)                              :: i          !< blade node counter
   INTEGER(IntKi)                              :: j          !< blade counter
               
   INTEGER(IntKi)                              :: ErrStat2
   CHARACTER(ErrMsgLen)                        :: ErrMsg2
   
   
      !...............................................................................................................................
      !  compute UA states at t=0 (with known inputs)
      !...............................................................................................................................
      if (p%UAMod == UA_HGM) then
      
         do j = 1,size(OtherState%UA_off_forGood,2) ! blades
            do i = 1,size(OtherState%UA_off_forGood,1) ! nodes

               ! We only update the UnsteadyAero states if we have unsteady aero turned on for this node
               if ( .not. OtherState%UA_off_forGood(i,j) .and. OtherState%FirstPass(i,j) ) then
               
                  ! initialize states to steady-state values:
                  call HGM_Steady( i, j, u(i,j), p, x%element(i,j), AFInfo(AFIndx(i,j)), ErrStat2, ErrMsg2 )
                     !call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                     
                  OtherState%FirstPass(i,j) = .false.
               end if

            end do
         end do
         
      end if

end subroutine UA_InitStates_AllNodes
!==============================================================================
SUBROUTINE HGM_Steady( i, j, u, p, x, AFInfo, ErrStat, ErrMsg )
! Routine to initialize the continuous states of the HGM model
!..................................................................................................................................

   INTEGER(IntKi),                      INTENT(IN   )  :: i           !< blade node counter
   INTEGER(IntKi),                      INTENT(IN   )  :: j           !< blade counter
   TYPE(UA_InputType),                  INTENT(IN   )  :: u           ! Inputs at t
   TYPE(UA_ParameterType),              INTENT(IN   )  :: p           ! Parameters
   TYPE(UA_ElementContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   type(AFI_ParameterType),             intent(in   )  :: AFInfo      ! The airfoil parameter data
   INTEGER(IntKi),                      INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                        INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables  
      
   type(AFI_UA_BL_Type)                         :: BL_p        ! potentially interpolated UA parameters
   type(AFI_OutputType)                         :: AFI_Interp
   character(ErrMsgLen)                         :: errMsg2
   integer(IntKi)                               :: errStat2
   character(*), parameter                      :: RoutineName = 'HGM_Steady'
   
   real(ReKi)                                   :: Tu
   real(ReKi)                                   :: alphaE
   real(ReKi)                                   :: alphaF
   real(ReKi)                                   :: alpha_34
  

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Lookup values using Airfoil Info module
   call AFI_ComputeUACoefs( AFInfo, u%Re, u%UserProp, BL_p, ErrMsg2, ErrStat2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
   call Get_HGM_constants(i, j, p, u, x, BL_p, Tu, alpha_34, alphaE) ! compute Tu, alpha_34, and alphaE
    
      
   ! States
   !x1: Downwash memory term 1 (rad)
   !x2: Downwash memory term 2 (rad)
   !x3: Clp', Lift coefficient with a time lag to the attached lift coeff
   !x4: f'' , Final separation point function


   ! Steady states
   x%x(1)     = BL_p%A1 * alpha_34
   x%x(2)     = BL_p%A2 * alpha_34
    
   alphaE   = alpha_34                                                    ! Eq. 12 (after substitute of x1 and x2 initializations)
   x%x(3)   = BL_p%c_lalpha * (alphaE-BL_p%alpha0)
    
      ! calculate x%x(4) = fs_aF = f_st(alphaF):
   alphaF  = x%x(3)/BL_p%c_lalpha + BL_p%alpha0                           ! p. 13
    
   call AFI_ComputeAirfoilCoefs( alphaF, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat, ErrMsg)
   x%x(4) = AFI_interp%f_st
   
end subroutine HGM_Steady
!----------------------------------------------------------------------------------------------------------------------------------
subroutine UA_CalcContStateDeriv( i, j, t, u_in, p, x, OtherState, AFInfo, m, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

   INTEGER(IntKi),                      INTENT(IN   )  :: i           !< blade node counter
   INTEGER(IntKi),                      INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                          INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(UA_InputType),                  INTENT(IN   )  :: u_in        ! Inputs at t
   TYPE(UA_ParameterType),              INTENT(IN   )  :: p           ! Parameters
   TYPE(UA_ElementContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(UA_OtherStateType),             INTENT(IN   )  :: OtherState  ! Other states at t
   type(UA_MiscVarType),                intent(inout)  :: m           ! Misc/optimization variables
   type(AFI_ParameterType),             intent(in   )  :: AFInfo      ! The airfoil parameter data
   TYPE(UA_ElementContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
   INTEGER(IntKi),                      INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                        INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables  
      
   type(AFI_UA_BL_Type)                         :: BL_p        ! potentially interpolated UA parameters
   type(AFI_OutputType)                         :: AFI_Interp
   character(ErrMsgLen)                         :: errMsg2
   integer(IntKi)                               :: errStat2
   character(*), parameter                      :: RoutineName = 'UA_CalcContStateDeriv'
   
   real(ReKi)                                   :: Tu
   real(ReKi)                                   :: alphaE
   real(ReKi)                                   :: alphaF
   real(ReKi)                                   :: Clp
   real(R8Ki)                                   :: x4
   real(ReKi)                                   :: alpha_34
   real(ReKi), parameter                        :: U_dot = 0.0_ReKi ! at some point we may add this term
   TYPE(UA_InputType)                           :: u        ! Inputs at t
  

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (OtherState%UA_off_forGood(i,j)) then
      dxdt%x = 0.0_R8Ki
      return
   end if
   
      ! make sure that u%u is not zero (this previously turned off UA for the entire simulation. 
      ! Now, we keep it on, but we don't want the math to blow up when we divide by u%u)
   call UA_fixInputs(u_in, u, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   
   
   ! Lookup values using Airfoil Info module
   call AFI_ComputeUACoefs( AFInfo, u%Re, u%UserProp, BL_p, ErrMsg2, ErrStat2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
      
   call Get_HGM_constants(i, j, p, u, x, BL_p, Tu, alpha_34, alphaE) ! compute Tu, alpha_34, and alphaE
    
   Clp = BL_p%c_lalpha * (alphaE - BL_p%alpha0) + pi * Tu * u%omega   ! Eq. 13
   
      ! fix definitions of T_f0 and T_p (based on email from Emmanuel 12-28-20 regarding HAWC2 default values)
   BL_p%T_f0 = BL_p%T_f0 * 2.0_ReKi * Tu
   BL_p%T_p  = BL_p%T_p  * Tu
      

      ! calculate fs_aF (stored in AFI_interp%f_st):
    
   !note: BL_p%c_lalpha cannot be zero. UA is turned off at initialization if this occurs.
   alphaF  = x%x(3)/BL_p%c_lalpha + BL_p%alpha0                           ! p. 13
   call AFI_ComputeAirfoilCoefs( alphaF, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
    ! States
    !x1: Downwash memory term 1 (rad)
    !x2: Downwash memory term 2 (rad)
    !x3: Clp', Lift coefficient with a time lag to the attached lift coeff
    !x4: f'' , Final separation point function
    
       ! Constraining x4 between 0 and 1 increases numerical stability (should be done elsewhere, but we'll double check here in case there were perturbations on the state value)
    x4 = max( min( x%x(4), 1.0_R8Ki ), 0.0_R8Ki )
    
if (p%ShedEffect) then
    dxdt%x(1) = -1.0_R8Ki / Tu * (BL_p%b1 + p%c(i,j) * U_dot/(2*u%u**2)) * x%x(1) + BL_p%b1 * BL_p%A1 / Tu * alpha_34
    dxdt%x(2) = -1.0_R8Ki / Tu * (BL_p%b2 + p%c(i,j) * U_dot/(2*u%u**2)) * x%x(2) + BL_p%b2 * BL_p%A2 / Tu * alpha_34
else
    dxdt%x(1) = 0.0_ReKi
    dxdt%x(2) = 0.0_ReKi
endif
    dxdt%x(3) = -1.0_R8Ki / BL_p%T_p                                     * x%x(3) +         1.0_ReKi / BL_p%T_p  * Clp
    dxdt%x(4) = -1.0_R8Ki / BL_p%T_f0                                    *    x4  +         1.0_ReKi / BL_p%T_f0 * AFI_interp%f_st

END SUBROUTINE UA_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Get_HGM_constants(i, j, p, u, x, BL_p, Tu, alpha_34, alphaE)
   INTEGER(IntKi),                      INTENT(IN   )  :: i           !< blade node counter
   INTEGER(IntKi),                      INTENT(IN   )  :: j           !< blade counter
   TYPE(UA_InputType),                  INTENT(IN   )  :: u           ! Inputs at t
   TYPE(UA_ParameterType),              INTENT(IN   )  :: p           ! Parameters
   TYPE(UA_ElementContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(AFI_UA_BL_Type),                INTENT(IN   )  :: BL_p        ! potentially interpolated UA parameters
   
   REAL(ReKi),                          INTENT(  OUT)  :: alpha_34
   REAL(ReKi),                          INTENT(  OUT)  :: Tu
   REAL(ReKi),                          INTENT(  OUT)  :: alphaE

      ! Local variables  
   real(ReKi)                                   :: vx_34

   
    ! Variables derived from inputs
    !u%u = U_ac = TwoNorm(u%v_ac)                                                 ! page 4 definitions
    
    Tu     = p%c(i,j) / (2.0_ReKi* max(u%u, UA_u_min))                            ! Eq. 23
    Tu     = min(Tu, 50.0_ReKi)   ! ensure the time constant doesn't exceed 50 s.
    Tu     = max(Tu,  0.001_ReKi) ! ensure the time constant doesn't get too small, either.

    vx_34 = u%v_ac(1) - u%omega * 0.5_ReKi*p%c(i,j)                        ! Eq. 1
    alpha_34 = atan2(vx_34, u%v_ac(2) )                                    ! page 5 definitions
    
    ! Variables derived from states
if (p%ShedEffect) then
    alphaE  = alpha_34*(1.0_ReKi - BL_p%A1 - BL_p%A2) + x%x(1) + x%x(2)    ! Eq. 12
else
    alphaE  = alpha_34
endif
    call MPi2Pi(alphaE)

END SUBROUTINE Get_HGM_constants
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE UA_RK4( i, j, t, n, u, utimes, p, x, OtherState, AFInfo, m, ErrStat, ErrMsg )
!..................................................................................................................................

   integer(IntKi),                  INTENT(IN   )  :: i           !< blade node counter
   integer(IntKi),                  INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                      INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   )  :: n           !< time step number
   TYPE(UA_InputType),              INTENT(IN   )  :: u(:)        !< Inputs at utimes
   REAL(DbKi),                      INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(UA_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(UA_ContinuousStateType),    INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(UA_OtherStateType),         INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(AFI_ParameterType),         INTENT(IN   )  :: AFInfo      ! The airfoil parameter data
   TYPE(UA_MiscVarType),            INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
         
   TYPE(UA_ElementContinuousStateType)             :: k1          ! RK4 constant; see above
   TYPE(UA_ElementContinuousStateType)             :: k2          ! RK4 constant; see above 
   TYPE(UA_ElementContinuousStateType)             :: k3          ! RK4 constant; see above 
   TYPE(UA_ElementContinuousStateType)             :: k4          ! RK4 constant; see above 
   TYPE(UA_ElementContinuousStateType)             :: x_tmp       ! Holds temporary modification to x
   TYPE(UA_InputType)                              :: u_interp    ! interpolated value of inputs
   
   REAL(DbKi)                                      :: TPlusHalfDt
   REAL(DbKi)                                      :: TPlusDt
   INTEGER(IntKi)                                  :: ErrStat2    ! local error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! local error message (ErrMsg)
   CHARACTER(*), PARAMETER                         :: RoutineName = 'UA_RK4'
      
      
   !NOTE: the error handling here assumes that we do not have any allocatable data in the inputs (u_interp) to be concerned with.
   !      Also, We assume that if there is going to be an error in UA_CalcContStateDeriv, it will happen only on the first call 
   !      to the routine.
   
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! interpolate u to find u_interp = u(t)
      CALL UA_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      x_tmp     = x%element(i,j)
      
      ! find xdot at t
      CALL UA_CalcContStateDeriv( i, j, t, u_interp, p, x_tmp, OtherState, AFInfo, m, k1, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k1%x     = p%dt * k1%x
  
      x_tmp%x  = x%element(i,j)%x + 0.5 * k1%x

      ! interpolate u to find u_interp = u(t + dt/2)
      TPlusHalfDt = t + 0.5_DbKi*p%dt
      CALL UA_Input_ExtrapInterp(u, utimes, u_interp, TPlusHalfDt, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt/2
      CALL UA_CalcContStateDeriv( i, j, TPlusHalfDt, u_interp, p, x_tmp, OtherState, AFInfo, m, k2, ErrStat2, ErrMsg2 )

      k2%x     = p%dt * k2%x

      x_tmp%x  = x%element(i,j)%x + 0.5 * k2%x

      ! find xdot at t + dt/2 (note x_tmp has changed)
      CALL UA_CalcContStateDeriv( i, j, TPlusHalfDt, u_interp, p, x_tmp, OtherState, AFInfo, m, k3, ErrStat2, ErrMsg2 )

      k3%x     = p%dt * k3%x

      x_tmp%x  = x%element(i,j)%x + k3%x

      ! interpolate u to find u_interp = u(t + dt)
      TPlusDt = t + p%dt
      CALL UA_Input_ExtrapInterp(u, utimes, u_interp, TPlusDt, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt
      CALL UA_CalcContStateDeriv( i, j, TPlusDt, u_interp, p, x_tmp, OtherState, AFInfo, m, k4, ErrStat2, ErrMsg2 )

      k4%x     = p%dt * k4%x

      x%element(i,j)%x = x%element(i,j)%x + ( k1%x + 2. * k2%x + 2. * k3%x + k4%x ) / 6.

END SUBROUTINE UA_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE UA_AB4( i, j, t, n, u, utimes, p, x, OtherState, AFInfo, m, ErrStat, ErrMsg )
!..................................................................................................................................

   integer(IntKi),               INTENT(IN   )  :: i           !< blade node counter
   integer(IntKi),               INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),               INTENT(IN   )  :: n           !< time step number
   TYPE(UA_InputType),           INTENT(IN   )  :: u(:)        !< Inputs at utimes
   REAL(DbKi),                   INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(UA_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(UA_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(UA_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(AFI_ParameterType),      INTENT(IN   )  :: AFInfo      ! The airfoil parameter data
   TYPE(UA_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   ! local variables
   TYPE(UA_InputType)                           :: u_interp
   TYPE(UA_ElementContinuousStateType)          :: x_tmp
   TYPE(UA_ElementContinuousStateType)          :: xdot
         
   INTEGER(IntKi)                               :: ErrStat2    ! local error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2     ! local error message (ErrMsg)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'UA_AB4'


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      if (OtherState%n(i,j) < n) then

         OtherState%n(i,j) = n
            
         OtherState%xdot(4)%element(i,j) = OtherState%xdot(3)%element(i,j)
         OtherState%xdot(3)%element(i,j) = OtherState%xdot(2)%element(i,j)
         OtherState%xdot(2)%element(i,j) = OtherState%xdot(1)%element(i,j)

      elseif (OtherState%n(i,j) > n) then

         CALL SetErrStat(ErrID_Fatal,'Backing up in time is not supported with a multistep method.',ErrStat,ErrMsg,RoutineName)
         RETURN

      endif
      

      ! need xdot at t, get inputs at t
      CALL UA_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
      x_tmp     = x%element(i,j)
      
      CALL UA_CalcContStateDeriv( i, j, t, u_interp, p, x_tmp, OtherState, AFInfo, m, xdot, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
         
        ! make sure OtherState%xdot( 1 ) is most up to date:
      OtherState%xdot( 1 )%element(i,j) = xdot

      if (n <= 2) then

         CALL UA_RK4(i, j, t, n, u, utimes, p, x, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

      else
         x%element(i,j)%x =  x%element(i,j)%x + p%DT/24. * ( 55.*OtherState%xdot(1)%element(i,j)%x - 59.*OtherState%xdot(2)%element(i,j)%x   &
                                                           + 37.*OtherState%xdot(3)%element(i,j)%x -  9.*OtherState%xdot(4)%element(i,j)%x )


      endif

      
END SUBROUTINE UA_AB4
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   Adams-Bashforth Predictor: \n
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector: \n
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE UA_ABM4( i, j, t, n, u, utimes, p, x, OtherState, AFInfo, m, ErrStat, ErrMsg )
!..................................................................................................................................

   integer(IntKi),               INTENT(IN   )  :: i           !< blade node counter
   integer(IntKi),               INTENT(IN   )  :: j           !< blade counter
   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),               INTENT(IN   )  :: n           !< time step number
   TYPE(UA_InputType),           INTENT(IN   )  :: u(:)        !< Inputs at utimes
   REAL(DbKi),                   INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(UA_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(UA_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(UA_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(AFI_ParameterType),      INTENT(IN   )  :: AFInfo      ! The airfoil parameter data
   TYPE(UA_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables

   TYPE(UA_InputType)                           :: u_interp    ! Inputs at t
   TYPE(UA_ElementContinuousStateType)          :: x_in        ! Continuous states at t
   TYPE(UA_ElementContinuousStateType)          :: xdot_pred   ! Derivative of continuous states at t

   INTEGER(IntKi)                               :: ErrStat2    ! local error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2     ! local error message (ErrMsg)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'UA_ABM4'
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! save copy of x(t):
      x_in     = x%element(i,j)
      
         ! predict: (note that we are overwritting x%element(i,j) here):
      CALL UA_AB4( i, j, t, n, u, utimes, p, x, OtherState, AFInfo, m, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      if (n > 2_IntKi) then
         
            ! correct:
         
         CALL UA_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         CALL UA_CalcContStateDeriv( i, j, t + p%dt, u_interp, p, x%element(i,j), OtherState, AFInfo, m, xdot_pred, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

         
         x%element(i,j)%x = x_in%x     + p%DT/24. * ( 9. * xdot_pred%x     + 19. * OtherState%xdot(1)%element(i,j)%x &
                                                                           -  5. * OtherState%xdot(2)%element(i,j)%x &
                                                                           +  1. * OtherState%xdot(3)%element(i,j)%x )

      endif
      
END SUBROUTINE UA_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------



!============================================================================== 
subroutine UA_CalcOutput( i, j, u_in, p, x, xd, OtherState, AFInfo, y, misc, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..............................................................................
   
   integer(IntKi),               intent(in   )  :: i           ! node index within a blade
   integer(IntKi),               intent(in   )  :: j           ! blade index 
   type(UA_InputType),           intent(in   )  :: u_in        ! Inputs at Time
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_ContinuousStateType), intent(in   )  :: x           ! Continuous states at Time
   type(UA_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   )  :: OtherState  ! Other states at Time
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
   type(UA_OutputType),          intent(inout)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   type(UA_MiscVarType),         intent(inout)  :: misc        ! Misc/optimization variables
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   integer(IntKi)                               :: errStat2        ! Error status of the operation (secondary error)
   character(ErrMsgLen)                         :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   character(*), parameter                      :: RoutineName = 'UA_CalcOutput'
   
   type(AFI_UA_BL_Type)                         :: BL_p        ! airfoil values computed in Kelvin Chain
   type(UA_KelvinChainType)                     :: KC          ! values computed in Kelvin Chain
   
   type(UA_InputType)                           :: u           ! Inputs at Time (with u%u set to appropriate value)
   
   
   real(ReKi)                                   :: Cm_FS
   real(ReKi)                                   :: Cc_FS
   real(ReKi)                                   :: Cm_alpha_nc
   real(ReKi)                                   :: M, f, k2_hat
   real(ReKi)                                   :: Cm_v, alpha_prime_f
   real(ReKi)                                   :: x_cp_hat                      ! center-of-pressure distance from LE in chord fraction
   real(ReKi)                                   :: Cm_common                     ! 
   real(ReKi)                                   :: k ! reduced frequency
   
   ! for UA_HGM
   real(ReKi)                                   :: alphaE
   real(ReKi)                                   :: Tu
   real(ReKi)                                   :: alpha_34
   real(ReKi)                                   :: fs_aE
   real(ReKi)                                   :: cl_fs
   real(ReKi)                                   :: x4
   real(ReKi)                                   :: delta_c_df_primeprime
   real(ReKi), parameter                        :: delta_c_mf_primeprime = 0.0_ReKi
   TYPE(UA_ElementContinuousStateType)          :: x_in        ! Continuous states at t
   

   type(AFI_OutputType)                         :: AFI_interp
   
#ifdef UA_OUTS
   integer                                      :: iOffset
#endif   
   
   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
   
   Cm_alpha_nc = 0.0_ReKi

   AFI_interp%Cm = 0.0_ReKi ! value will be output if not computed below
   alpha_prime_f = 0.0_ReKi ! value will be output if not computed below
   
      ! make sure that u%u is not zero (this previously turned off UA for the entire simulation. 
      ! Now, we keep it on, but we don't want the math to blow up when we divide by u%u)
   call UA_fixInputs(u_in, u, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   k = abs(u%omega * p%c(i, j) / (2.0_ReKi* u%u))
   
   if (  OtherState%UA_off_forGood(i,j) .or. (OtherState%FirstPass(i, j) .and. p%UAMod /= UA_HGM) ) then ! note: if u%U isn't zero because we've called UA_fixInputs
        
      misc%weight(i,j) = 0.0
      
      call AFI_ComputeAirfoilCoefs( u%alpha, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
      y%Cl = AFI_interp%Cl
      y%Cd = AFI_interp%Cd
      y%Cm = AFI_interp%Cm
      
      
      y%Cn = y%Cl*cos(u%alpha) + (y%Cd-AFI_interp%Cd0)*sin(u%alpha)
      y%Cc = y%Cl*sin(u%alpha) - (y%Cd-AFI_interp%Cd0)*cos(u%alpha)
      
         Cm_v              = 0.0_ReKi
      KC%Cn_alpha_q_circ   = 0.0_ReKi
      KC%Cn_alpha_q_nc     = 0.0_ReKi
      KC%Cn_pot            = 0.0_ReKi
      KC%Dp                = 0.0_ReKi
      KC%Cn_prime          = 0.0_ReKi
      KC%fprime            = 0.0_ReKi
      KC%Df                = 0.0_ReKi
      KC%Cn_V              = 0.0_ReKi
      KC%C_V               = 0.0_ReKi
      KC%Cm_q_circ         = 0.0_ReKi
      KC%Cm_q_nc           = 0.0_ReKi
      KC%Dalphaf           = 0.0_ReKi
      KC%T_f               = 0.0_ReKi
      KC%T_V               = 0.0_ReKi
      KC%alpha_filt_cur    = u%alpha
      KC%ds                = 2.0_ReKi*u%U*p%dt/p%c(i, j)
      
      alphaE   = 0.0
      Tu       = 0.0
      alpha_34 = 0.0
      cl_fs    = 0.0
      fs_aE    = 0.0
      
   elseif (p%UAMod == UA_HGM) then
      x_in = x%element(i,j)
      
      if (OtherState%FirstPass(i,j)) then
         call HGM_Steady( i, j, u, p, x_in, AFInfo, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end if
      
      call AFI_ComputeUACoefs( AFInfo, u%Re, u%UserProp, BL_p, ErrMsg2, ErrStat2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
      
      call Get_HGM_constants(i, j, p, u, x_in, BL_p, Tu, alpha_34, alphaE) ! compute Tu, alpha_34, and alphaE
   
      call AFI_ComputeAirfoilCoefs( alphaE, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
      
         ! calculate fs_aE and cl_fs:
      cl_fs = AFI_interp%cl_fs
      fs_aE = AFI_interp%f_st
    
       ! Constraining x4 between 0 and 1 increases numerical stability (should be done elsewhere, but we'll double check here in case there were perturbations on the state value)
      x4 = max( min( x_in%x(4), 1.0_R8Ki ), 0.0_R8Ki )
    
      delta_c_df_primeprime = 0.5_ReKi * (sqrt(fs_aE) - sqrt(x4)) - 0.25_ReKi * (fs_aE - x4)
      
! bjj: do we need to check that u%alpha is between -pi and + pi?
      ! y%Cl = AFI_interp%Cl < TODO consider using this in front of x4 for "true" Cl
      y%Cl = x4 * (alphaE - BL_p%alpha0) * BL_p%c_lalpha  + (1.0_ReKi - x4) * cl_fs  + pi * Tu * u%omega       ! Eq. 78
      y%Cd = AFI_interp%Cd + (u%alpha - alphaE) * y%Cl + (AFI_interp%Cd - BL_p%Cd0) * delta_c_df_primeprime    ! Eq. 79
      
      if (AFInfo%ColCm == 0) then ! we don't have a cm column, so make everything 0
         y%Cm          = 0.0_ReKi
      else
         y%Cm = AFI_interp%Cm + y%Cl * delta_c_mf_primeprime - piBy2 * Tu * u%omega                            ! Eq. 80
      end if

      y%Cn = y%Cl*cos(u%alpha) + y%Cd*sin(u%alpha)
      y%Cc = y%Cl*sin(u%alpha) - y%Cd*cos(u%alpha)
      
         ! now check if we should have turned off UA, and modify outputs accordingly (with linear combination of steady outputs)
      call UA_BlendSteady(k, u, p, AFInfo, y, misc%FirstWarn_UA_off, misc%weight(i,j), ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   else
      
      M           = u%U / p%a_s
      
      call UA_CheckMachNumber(M, misc%FirstWarn_M, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
            
      call ComputeKelvinChain( i, j, u, p, xd, OtherState, misc, AFInfo, KC, BL_p, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      

      !.............................
      ! Cn
      !.............................
         ! Eqn 1.53 or 1.53b depending on UAMod
      if (xd%tau_v(i, j) > 0.0) then
         y%Cn= KC%Cn_FS + KC%Cn_v                                                                                                                  ! Eqn 1.53
      else
         y%Cn= KC%Cn_FS
      end if

      !.............................
      ! Cc
      !.............................
      
      
      if ( p%flookup ) then 
      !if ( p%UAMod == UA_Gonzalez ) then
         Cc_FS = BL_p%eta_e*KC%Cc_pot *(sqrt(KC%fprimeprime_c) - Gonzalez_factor)                                                                  ! Eqn 1.40  
      else
         Cc_FS = BL_p%eta_e*KC%Cc_pot * sqrt(KC%fprimeprime_c)                                                                                     ! Eqn 1.40 without Gonzalez's modification for negative Cc
      end if
      
            
      if ( p%UAMod == UA_MinnemaPierce ) then
#ifdef TEST_THEORY
            y%Cc = Cc_FS + KC%Cn_v*tan(KC%alpha_e)*(1-xd%tau_v(i, j)/(BL_p%T_VL))                                          ! Eqn 1.55 with Eqn. 1.40
#else            
            y%Cc = Cc_FS + KC%Cn_v*    KC%alpha_e *(1.0_ReKi-xd%tau_v(i, j)/(BL_p%T_VL))                                   ! Eqn 1.55 with approximation of tan(KC%alpha_e)=KC%alpha_e and Eqn. 1.40 substitution
#endif            
                 
      elseif ( p%UAMod == UA_Gonzalez ) then
         y%Cc = Cc_FS                                                                                                                              ! Eqn. 1.55b
      else
         ! TODO: This is UAMod = 1 and we still need to verify these equations!!!  GJH Dec 20 2015
         if ( KC%Cn_prime <= BL_p%Cn1 ) then
            !y%Cc = BL_p%eta_e*KC%Cc_pot*(sqrt(KC%fprimeprime_c) - f_c_offset) !*sin(KC%alpha_e + BL_p%alpha0)
            y%Cc = Cc_FS * sin(KC%alpha_e + BL_p%alpha0)                                                                                           ! Eqn 1.56a [1]
         else
            if ( p%flookup ) then 
              ! if (p%UAMod == UA_Baseline) then
              !    call Get_f_from_Lookup( p%UAMod, u%Re, u%UserProp, KC%alpha_filt_cur, BL_p%alpha0, KC%C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2, f=f)
              ! else   
               ! TODO: Need to understand the effect of the offset on this f in the following equations and fprimeprime_c.
                  !bjj: fprimeprime_c is computed with the Gonzalez offset in the Kelvin Chain, so it's not just f that could be a problem.
                  f      = Get_f_c_from_Lookup( p%UAMod, u%Re, u%UserProp, KC%alpha_filt_cur, BL_p%alpha0, KC%C_nalpha_circ, BL_p%eta_e, AFInfo, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
              ! endif
               
            else
               f      = Get_f( KC%alpha_filt_cur, BL_p%alpha0, BL_p%alpha1, BL_p%alpha2, BL_p%S1, BL_p%S2, BL_p%S3, BL_p%S4 )
            end if
            
            k2_hat = 2.0_ReKi*(KC%Cn_prime - BL_p%Cn1) + KC%fprimeprime_c - f                                                                      ! Eqn 1.56b
            ! TODO: why did we comment out the sin() term below?  GJH 2/28/2017
            y%Cc   = BL_p%k1_hat + KC%Cc_pot*sqrt(KC%fprimeprime_c)*KC%fprimeprime_c**k2_hat *sin(KC%alpha_e + BL_p%alpha0)                        ! Eqn 1.56a [2]
            
         end if
         
      end if
      
      !.............................
      ! convert cc and cn to cl and cd
      !.............................
            
      y%Cl = y%Cn*cos(u%alpha) + y%Cc*sin(u%alpha)                                                                                                 ! Eqn 1.2a
      y%Cd = y%Cn*sin(u%alpha) - y%Cc*cos(u%alpha) + BL_p%Cd0                                                                                      ! Eqn 1.2b 

         ! Make Cn and CC consistent with the added contribution of Cd0 in Cd:
      y%Cn = y%Cl*cos(u%alpha) + y%Cd*sin(u%alpha)    !Added the contribution of Cd0 in Cn and Cc
      y%Cc = y%Cl*sin(u%alpha) - y%Cd*cos(u%alpha)
      
      !.............................
      ! convert cm
      !.............................
      
      alpha_prime_f = 0.0_ReKi ! initialize for output purposes
      
         ! Check for Cm column in AFInfo data
      if (AFInfo%ColCm == 0) then ! we don't have a cm column, so make everything 0
         
         y%Cm          = 0.0_ReKi
         Cm_v          = 0.0_ReKi
     
      else
            
         if ( p%UAMod == UA_Gonzalez ) then !bjj: what is this doing?
            
            if (abs(KC%alpha_f - BL_p%alpha0) < .01) then
               if (KC%alpha_f < BL_p%alpha0) then
                  KC%alpha_f = KC%alpha_f - .01
               else
                  KC%alpha_f = KC%alpha_f + .01
               end if
            end if       ! Removed part of the code related to fprimeprime_m, which has been added in other part of the code        
         
         end if
      
         !...........
         ! compute Cm_FS
         !...........
         
         Cm_alpha_nc = - KC%Cn_alpha_nc / 4.0_ReKi                                                                                                 ! Eqn 1.27
         Cm_common   = KC%Cm_q_circ + Cm_alpha_nc + KC%Cm_q_nc                                                                                     ! second parts of Eqn 1.41, Eqn 1.44, and Eqn 1.45
   
         if ( p%UAMod == UA_Baseline ) then
            
            x_cp_hat = BL_p%k0 + BL_p%k1*(1.0_ReKi-KC%fprimeprime) + BL_p%k2*sin(pi*KC%fprimeprime**BL_p%k3)                                       ! Eqn 1.42
            Cm_FS  = BL_p%Cm0 - KC%Cn_alpha_q_circ*(x_cp_hat - 0.25_ReKi) + Cm_common                                                              ! Eqn 1.41

         elseif ( p%UAMod == UA_MinnemaPierce ) then
      
               ! Look up Cm using alpha_prime_f
            alpha_prime_f = KC%alpha_f - KC%Dalphaf                                                                                                ! Eqn 1.43a

            call AFI_ComputeAirfoilCoefs( alpha_prime_f, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)            
            Cm_FS = AFI_interp%Cm + Cm_common                                                                                                      ! Eqn 1.44
         else ! UAMod == UA_Gonzalez
            Cm_FS = BL_p%Cm0 + KC%Cn_FS*KC%fprimeprime_m + Cm_common                                                                               ! Eqn 1.45

         end if   
      
         Cm_v     = -BL_p%x_cp_bar*( 1.0_ReKi - cos( pi*xd%tau_v(i, j)/BL_p%T_VL ) )*KC%Cn_v                               ! Eqn 1.57
         if (p%UAMod == UA_Gonzalez .and. xd%tau_v(i, j) <= 0.0 ) then !Added specific logic for UAMod=UA_Gonzalez
            y%Cm   = Cm_FS
         else
            y%Cm   = Cm_FS + Cm_v                                                                                                                  ! Eqn 1.58, Eqn 1.59, and Eqn 1.60
         end if

      end if
      
         ! now check if we should have turned off UA, and modify outputs accordingly (with linear combination of steady outputs)
      call UA_BlendSteady(k, u, p, AFInfo, y, misc%FirstWarn_UA_off, misc%weight(i,j), ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   end if
   
#ifdef UA_OUTS
   iOffset = (i-1)*p%NumOuts + (j-1)*p%nNodesPerBlade*p%NumOuts
   if (allocated(y%WriteOutput)) then  !bjj: because BEMT uses local variables for UA output, y%WriteOutput is not necessarially allocated. Need to figure out a better solution.
   
      y%WriteOutput(iOffset+ 1)    = u%alpha*R2D
      y%WriteOutput(iOffset+ 2)    = u%U
      y%WriteOutput(iOffset+ 3)    = y%Cn
      y%WriteOutput(iOffset+ 4)    = y%Cc
      y%WriteOutput(iOffset+ 5)    = y%Cl
      y%WriteOutput(iOffset+ 6)    = y%Cd
      y%WriteOutput(iOffset+ 7)    = y%Cm
   
      if (p%UAMod == UA_HGM) then

         y%WriteOutput(iOffset+ 8)    = u%omega*R2D
         y%WriteOutput(iOffset+ 9)    = alphaE*R2D
         y%WriteOutput(iOffset+10)    = Tu
         y%WriteOutput(iOffset+11)    = alpha_34*R2D
         y%WriteOutput(iOffset+12)    = cl_fs
         y%WriteOutput(iOffset+13)    = fs_aE
         
         y%WriteOutput(iOffset+14)    = x%element(i,j)%x(1)
         y%WriteOutput(iOffset+15)    = x%element(i,j)%x(2)
         y%WriteOutput(iOffset+16)    = x%element(i,j)%x(3)
         y%WriteOutput(iOffset+17)    = x%element(i,j)%x(4)
         y%WriteOutput(iOffset+18)    = k
         y%WriteOutput(iOffset+19)    = misc%weight(i,j)

      else
         y%WriteOutput(iOffset+ 8)    = KC%Cn_alpha_q_circ               ! CNCP in ADv14
         y%WriteOutput(iOffset+ 9)    = KC%Cn_alpha_q_nc                 ! CNIQ in ADv14
         y%WriteOutput(iOffset+10)    = KC%Cn_pot
         y%WriteOutput(iOffset+11)    = KC%Dp
         y%WriteOutput(iOffset+12)    = KC%Cn_prime
         y%WriteOutput(iOffset+13)    = KC%fprime
         y%WriteOutput(iOffset+14)    = KC%Df
         y%WriteOutput(iOffset+15)    = KC%Cn_V
         y%WriteOutput(iOffset+16)    = xd%tau_v(i, j)

         if ( misc%LESF(i, j) ) then  ! BEDSEP in v14 !bjj: note that this output is calculated in UpdateDiscState, so it may be out of sync here.
            y%WriteOutput(iOffset+17) = 1.0_ReKi
         else
            y%WriteOutput(iOffset+17) = 0.0_ReKi
         end if

         if ( misc%TESF(i, j) ) then  !SHIFT in v14   !bjj: note that this output is calculated in UpdateDiscState, so it may be out of sync here.
            y%WriteOutput(iOffset+18) = 1.0_ReKi
         else
            y%WriteOutput(iOffset+18) = 0.0_ReKi
         end if
         if ( misc%VRTX(i, j) ) then  ! VOR in v14    !bjj: note that this output is calculated in UpdateDiscState, so it may be out of sync here.
            y%WriteOutput(iOffset+19) = 1.0_ReKi
         else
            y%WriteOutput(iOffset+19) = 0.0_ReKi
         end if
         y%WriteOutput(iOffset+20)    = KC%C_V
         y%WriteOutput(iOffset+21)    = Cm_alpha_nc
         y%WriteOutput(iOffset+22)    = KC%Cm_q_nc
         y%WriteOutput(iOffset+23)    = Cm_v
         y%WriteOutput(iOffset+24)    = alpha_prime_f
         y%WriteOutput(iOffset+25)    = KC%Dalphaf
         y%WriteOutput(iOffset+26)    = AFI_interp%Cm
         y%WriteOutput(iOffset+27)    = KC%T_f
         y%WriteOutput(iOffset+28)    = KC%T_V
         y%WriteOutput(iOffset+29)    = KC%ds  ! Eqn 1.51 (tau_v) 
         y%WriteOutput(iOffset+30)    = KC%T_alpha
         y%WriteOutput(iOffset+31)    = KC%T_q
         y%WriteOutput(iOffset+32)    = KC%k_alpha
         y%WriteOutput(iOffset+33)    = KC%k_q
         y%WriteOutput(iOffset+34)    = KC%alpha_e
         y%WriteOutput(iOffset+35)    = KC%X1
         y%WriteOutput(iOffset+36)    = KC%X2
         y%WriteOutput(iOffset+37)    = KC%cn_q_nc
         y%WriteOutput(iOffset+38)    = KC%alpha_f
         y%WriteOutput(iOffset+39)    = KC%fprimeprime
         y%WriteOutput(iOffset+40)    = OtherState%sigma1(i, j)
         y%WriteOutput(iOffset+41)    = OtherState%sigma3(i, j)
         y%WriteOutput(iOffset+42)    = misc%T_sh(i, j)
         y%WriteOutput(iOffset+43)    = k
         y%WriteOutput(iOffset+44)    = KC%alpha_filt_cur*R2D
         y%WriteOutput(iOffset+45)    = misc%weight(i,j)
      end if
   end if
#endif
   
end subroutine UA_CalcOutput
!==============================================================================   
subroutine UA_WriteOutputToFile(t, p, y)
   real(DbKi),                   intent(in   )  :: t           ! current time (s)
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_OutputType),          intent(in   )  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
!   type(UA_ContinuousStateType), intent(in   )  :: x           ! Continuous states at Time
!   type(UA_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
!   type(UA_OtherStateType),      intent(in   )  :: OtherState  ! Other states at Time
!   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
!   type(UA_MiscVarType),         intent(inout)  :: misc        ! Misc/optimization variables
!   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
!   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   integer                                   :: k
      
      ! Generate file outputs
#ifdef UA_OUTS
   if (p%unOutFile > 0 .and. allocated(y%WriteOutput)) then
   
      write (p%unOutFile,"(F19.6)",ADVANCE='no')  t
      do k=1,size(y%WriteOutput)
         WRITE (p%unOutFile,'(:,A,'//trim( p%OutFmt )//')', ADVANCE='no' )  p%Delim, y%WriteOutput(k)
      end do  
      WRITE (p%unOutFile,'()')          ! write the line return

   end if
#endif

end subroutine UA_WriteOutputToFile
!==============================================================================
subroutine UA_End(p)
   type(UA_ParameterType),       intent(inout)  :: p           ! Parameters
!   type(UA_ContinuousStateType), intent(in   )  :: x           ! Continuous states at Time
!   type(UA_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
!   type(UA_OtherStateType),      intent(in   )  :: OtherState  ! Other states at Time
!   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
!   type(UA_OutputType),          intent(inout)  :: y           ! Outputs computed at Time (Input only so that mesh con-
!                                                               !   nectivity information does not have to be recalculated)
!   type(UA_MiscVarType),         intent(inout)  :: misc        ! Misc/optimization variables
!   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
!   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   if (p%NumOuts > 0 .and. p%UnOutFile > 0) CLOSE(p%UnOutFile)
   p%unOutFile = -1
end subroutine UA_End
!==============================================================================   
!>This subroutine blends the steady outputs with the unsteady-outputs so that
!! UA can turn back on if the angle of attack goes back into a reasonable range.
subroutine UA_BlendSteady(k, u, p, AFInfo, y, FirstWarn_UA_off, weight, ErrStat, ErrMsg)
   REAL(ReKi),                   intent(in   )  :: k           ! reduced frequency
   type(UA_InputType),           intent(in   )  :: u           ! "Fixed" Inputs at Time
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(AFI_ParameterType),      intent(in   )  :: AFInfo      ! The airfoil parameter data
   type(UA_OutputType),          intent(inout)  :: y           ! Outputs computed at Time (Input only so that mesh con-
   LOGICAL,                      intent(inout)  :: FirstWarn_UA_off      ! flag to determine if warning message should be displayed
   REAL(ReKi),                   intent(  out)  :: weight      ! scaling weight for UA vs steady outputs
   integer(IntKi),               intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   type(AFI_OutputType)                         :: AFI_steady
   REAL(ReKi)                                   :: W1,W2,W3      ! Weights for turning off UA temporarily
   REAL(ReKi)                                   :: AFI_steady_Cn ! Cn from steady coefficients
   REAL(ReKi)                                   :: AFI_steady_Cc ! Cc from steady coefficients
   TYPE(AFI_UA_BL_Type)                         :: UA_BL         ! The tables of Leishman-Beddoes unsteady-aero data for given Re and control setting [-]
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   CHARACTER(*), PARAMETER                      :: RoutineName = 'UA_BlendSteady'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   weight = 1.0_ReKi ! default in case of error
   
      ! Determine what the cutout angle of attack is
   call AFI_ComputeUACoefs( AFInfo, u%Re, u%UserProp, UA_BL, ErrMsg, ErrStat )
   if (ErrStat >= AbortErrLev) return ! would only have error if there is a memory problem
      
      ! put alpha in [-pi,pi] before checking its value
   
   W1 = 1.0_ReKi - BlendCosine( abs(u%alpha), UA_BL%UACutout_blend, UA_BL%UACutout ) ! shut off when AoA reaches UACutout, but blend it off 5 degrees before (5 degrees is set in AFI to avoid that math each time)
   W2 =            BlendCosine( abs(u%U),                 UA_u_min, 1.0_ReKi       ) ! turn off UA when inflow velocity is 0 m/s, but start blend it off 1 m/s before (make sure it is greater than u_min)
  !W3 =            BlendCosine( k,   0.0_ReKi, 0.02_ReKi       ) ! turn off UA when reduced frequency is 0, but start blend it off at k=0.02 (this is a quasi-static state)
   W3 = 1.0_ReKi
   weight = W1*W2*W3

   if (weight < 1.0_ReKi) then
   
      if (FirstWarn_UA_off) then
         CALL SetErrStat(ErrID_Warn,"Temporarily turning off UA due to high angle of attack or low relative velocity. This warning will not be repeated though the condition may persist.", ErrStat, ErrMsg, RoutineName)
         FirstWarn_UA_off = .false.
      end if
   
      ! calculate the steady coefficients
      call AFI_ComputeAirfoilCoefs( u%alpha, u%Re, u%UserProp, AFInfo, AFI_steady, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      AFI_steady_Cn = AFI_steady%Cl*cos(u%alpha) + (AFI_steady%Cd-AFI_steady%Cd0)*sin(u%alpha)
      AFI_steady_Cc = AFI_steady%Cl*sin(u%alpha) - (AFI_steady%Cd-AFI_steady%Cd0)*cos(u%alpha)
         
      y%Cn = y%Cn*weight + (1.0_ReKi - weight)*AFI_steady_Cn
      y%Cc = y%Cc*weight + (1.0_ReKi - weight)*AFI_steady_Cc
      y%Cm = y%Cm*weight + (1.0_ReKi - weight)*AFI_steady%Cm
      y%Cl = y%Cl*weight + (1.0_ReKi - weight)*AFI_steady%Cl
      y%Cd = y%Cd*weight + (1.0_ReKi - weight)*AFI_steady%Cd
         
      
   end if


end subroutine UA_BlendSteady
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


subroutine UA_fixInputs(u_in, u, errStat, errMsg)
   type(UA_InputType),           intent(in   )  :: u_in        ! Inputs at Time
   type(UA_InputType),           intent(inout)  :: u           ! Inputs at Time

   integer(IntKi)                               :: errStat        ! Error status of the operation (secondary error)
   character(ErrMsgLen)                         :: errMsg         ! Error message if ErrStat2 /= ErrID_None
   
      ! make sure that u%u is not zero (this previously turned off UA for the entire simulation. 
      ! Now, we keep it on, but we don't want the math to blow up when we divide by u%u)
   call UA_CopyInput(u_in, u, MESH_UPDATECOPY, errStat, errMsg)
   
   call mPi2Pi(u%alpha) ! make sure alpha is in a good range
   
   if (abs(u%u) < UA_u_min) then
      u%u = sign(UA_u_min, u%u)
      
      u%v_ac(1) = sin(u%alpha)*u%U
      u%v_ac(2) = cos(u%alpha)*u%U
   end if
   
end subroutine UA_fixInputs
   
   
end module UnsteadyAero
