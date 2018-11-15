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
   type(ProgDesc), parameter  :: UA_Ver = ProgDesc( 'UnsteadyAero', '', '' )

   public :: UA_Init
   public :: UA_UpdateDiscOtherState
   public :: UA_UpdateStates
   public :: UA_CalcOutput

   public :: UA_ReInit

   integer(intki), parameter :: UA_Baseline      = 1   ! UAMod = 1 [Baseline model (Original)]
   integer(intki), parameter :: UA_Gonzalez      = 2   ! UAMod = 2 [Gonzalez's variant (changes in Cn,Cc,Cm)]
   integer(intki), parameter :: UA_MinemmaPierce = 3   ! UAMod = 3 [Minemma/Pierce variant (changes in Cc and Cm)]
   
   real(ReKi),     parameter :: Gonzales_factor = 0.2_ReKi   ! this factor, proposed by Gonzales (for "all" models) is used to modify Cc to account for negative values seen at f=0 (see Eqn 1.40)
   
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
    
   
   ! this routine is solving Eqn 1.32 for f:
   ! cn = C_nalpha_circ * (alpha-alpha0) * ((1+sqrt(f))/2)**2
   
   
   !bjj: if cn = 0 or c_nalpha_circ = 0 or alpha=alpha0, f has infinitely many solutions to this equation
   
   if (EqualRealNos( real(c_nalpha_circ,SiKi), 0.0_SiKi )) then
      tmpRoot = 0.0_ReKi
   else if (EqualRealNos( real(alpha,SiKi), real(alpha0,SiKi) )) then
      tmpRoot = 0.0_ReKi
   else
      
      call GetSteadyOutputs(AFInfo, alpha, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
         if (ErrStat >= AbortErrLev ) return
   
      Cn =  Cl*cos(alpha) + (Cd-Cd0)*sin(alpha)
   
      if (EqualRealNos( real(cn,SiKi), 0.0_SiKi )) then
         tmpRoot = 0.0_ReKi
      else 
      
         denom = (C_nalpha_circ*(alpha-alpha0))
      
         !    bjj: if tmpRoot=cn/(C_nalpha_circ*(alpha-alpha0)) is negative, this whole equation is bogus....
         tmpRoot  = Cn/denom
   
         if (tmpRoot < 0.0_ReKi) then
            tmpRoot = 0.0_ReKi
         end if
                        
      end if
            
   end if
   
   if (UAMod == UA_Gonzalez) then
      Get_f_from_Lookup = ((3*sqrt(tmpRoot)-1)/2.0)**2
   else
      Get_f_from_Lookup = ( 2 * sqrt( tmpRoot ) - 1 ) **2 
   end if
   
   
   if ( Get_f_from_Lookup > 1.0 ) then
      Get_f_from_Lookup = 1.0_ReKi
   end if
   
   
end function Get_f_from_Lookup      


!==============================================================================
real(ReKi) function Get_f_c_from_Lookup( Re, alpha, alpha0, c_nalpha_circ, eta_e, AFInfo, ErrStat, ErrMsg)
! Compute either fprime or fprimeprime using an analytical equation (and eventually a table lookup)
! Called by : ComputeKelvinChain
! Calls  to : NONE
!..............................................................................

   real(ReKi),       intent(in   ) :: Re            ! Reynolds number
   real(ReKi),       intent(in   ) :: alpha         ! angle of attack (radians)
   real(ReKi),       intent(in   ) :: alpha0
   real(ReKi),       intent(in   ) :: c_nalpha_circ
   real(ReKi),       intent(in   ) :: eta_e
   type(AFInfoType), intent(in   ) :: AFInfo        ! The airfoil parameter data   
   integer(IntKi),   intent(  out) :: ErrStat       ! Error status of the operation
   character(*),     intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
   
   
   real(ReKi), parameter           :: fc_limit = (1.0_ReKi + Gonzales_factor)**2    ! normally, fc is limited by 1, but we're limiting (sqrt(fc)-Gonzales_factor) to 1
   real(ReKi)                      :: Cc, Cl, Cd, Cm, Cd0, denom
   ErrStat = ErrID_None
   ErrMsg  = ''
      ! NOTE:  This subroutine call cannot live in Blade Element because BE module calls UnsteadyAero module.
   
  
   if (EqualRealNos(real(alpha,SiKi), 0.0_SiKi)) then
      
      Get_f_c_from_Lookup = 0.0_ReKi
      
   elseif (EqualRealNos(real(alpha,SiKi), real(alpha0,SiKi))) then
      
      Get_f_c_from_Lookup = 0.0_ReKi
      
   else if (EqualRealNos( real(c_nalpha_circ,SiKi), 0.0_SiKi )) then
      
      Get_f_c_from_Lookup = 0.0_ReKi

   else
         
      call GetSteadyOutputs(AFInfo, alpha, Cl, Cd, Cm, Cd0, ErrStat, ErrMsg)
         if (ErrStat >= AbortErrLev) return
      Cc =  Cl*sin(alpha) - (Cd-Cd0)*cos(alpha)
   
      denom = eta_e*c_nalpha_circ*( alpha-alpha0 )*tan(alpha)
      Get_f_c_from_Lookup =  min(fc_limit, (  Cc / denom  + Gonzales_factor ) **2 )
                  
   end if
      ! Apply an offset of 0.2 to fix cases where f_c should be negative, but we are using **2 so can only return positive values
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

 





!==============================================================================                              
subroutine ComputeKelvinChain( i, j, u, p, xd, OtherState, misc, AFInfo, KC, BL_p, ErrStat, ErrMsg )
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
   
   
   type(AFI_UA_BL_Type),                   intent(  out) :: BL_p

   integer(IntKi),                         intent(  out) :: ErrStat           ! Error status of the operation
   character(*),                           intent(  out) :: ErrMsg            ! Error message if ErrStat /= ErrID_None

            
   type(UA_KelvinChainType),               intent(  out) :: KC

   real(ReKi)                :: ds                                            ! non-dimensionalized distance parameter
   real(ReKi)                :: M                                             ! Mach number (-)
   real(ReKi)                :: beta_M                                        ! Prandtl-Glauert compressibility correction factor,  sqrt(1-M**2)
   real(ReKi)                :: beta_M_Sqrd                                   ! square of the Prandtl-Glauert compressibility correction factor,  (1-M**2)
                 
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
      
   beta_M_Sqrd = 1 - M**2
   beta_M      = sqrt(beta_M_Sqrd) 
   
   ! Lookup values using Airfoil Info module
   call AFI_GetAirfoilParams( AFInfo, M, u%Re, BL_p, KC%C_nalpha_circ, ErrMsg2, ErrStat2 )           
      ! AFI_GetAirfoilParams doesn't return error, so I will not check
   
   !bjj: shouldn't this happen in AFI_GetAirfoilParams instead?
      ! Override eta_e if we are using Flookup
   if ( p%Flookup ) then
      BL_p%eta_e = 1.0
   end if
   
   ! Kelvin chain
   ds          = 2.0_ReKi*u%U*p%dt/p%c(i,j)                             ! Eqn 1.5b
   
   if (OtherState%FirstPass(i,j)) then
      alpha_minus1 = u%alpha      
      alpha_filt_minus1 = u%alpha     
   else
      alpha_minus1 = xd%alpha_minus1(i,j)     
      alpha_filt_minus1 = xd%alpha_filt_minus1(i,j)   
   end if
   
   ! Low Pass filtering of alpha, q, and Kq in order to deal with numerical issues that arise for very small values of dt
   ! See equations 1.7 and 1.8 of the manual
   
   LowPassConst  =  exp(-2.0_ReKi*PI*p%dt*BL_p%filtCutOff) ! from Eqn 1.8 [7]
   
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
  
   KC%Kprime_alpha  = Get_ExpEqn( real(p%dt,ReKi), KC%T_alpha, xd%Kprime_alpha_minus1(i,j), KC%Kalpha_f, Kalpha_f_minus1 )    ! Eqn 1.18b
   KC%Cn_alpha_nc   = 4.0_ReKi*KC%T_alpha * ( KC%Kalpha_f - KC%Kprime_alpha ) / M                                             ! Eqn 1.18a
   
   KC%Kprime_q      = Get_ExpEqn( real(p%dt,ReKi), KC%T_q    , xd%Kprime_q_minus1(i,j)    ,  KC%Kq_f   , Kq_f_minus1     )    ! Eqn 1.19b 
   KC%Cn_q_nc       = -1.0_ReKi*KC%T_q * ( KC%Kq_f - KC%Kprime_q ) / M                                                        ! Eqn 1.19a
         
   KC%Cn_alpha_q_nc = KC%Cn_alpha_nc + KC%Cn_q_nc                                                                             ! Eqn 1.17
   
   KC%X1            = Get_ExpEqn( ds*beta_M_Sqrd*BL_p%b1, 1.0_ReKi, xd%X1_minus1(i,j), BL_p%A1*(KC%alpha_filt_cur - alpha_filt_minus1), 0.0_ReKi ) ! Eqn 1.15a
   KC%X2            = Get_ExpEqn( ds*beta_M_Sqrd*BL_p%b2, 1.0_ReKi, xd%X2_minus1(i,j), BL_p%A2*(KC%alpha_filt_cur - alpha_filt_minus1), 0.0_ReKi ) ! Eqn 1.15b
   
   KC%alpha_e       = (KC%alpha_filt_cur - BL_p%alpha0) - KC%X1 - KC%X2                                                       ! Eqn 1.14
   
   KC%Cn_alpha_q_circ = KC%C_nalpha_circ * KC%alpha_e                                                                         ! Eqn 1.13
   
   if ( p%UAMod == UA_Gonzalez ) then
         ! Compute X3 and X4 using Eqn 1.16a  and then add Cn_q_circ (Eqn 1.16) to the previously computed Cn_alpha_q_circ
      KC%X3              = Get_ExpEqn( ds*beta_M_Sqrd*BL_p%b1, 1.0_ReKi, xd%X3_minus1(i,j), BL_p%A1*(KC%q_f_cur - q_f_minus1), 0.0_ReKi ) ! Eqn 1.16a [1]
      KC%X4              = Get_ExpEqn( ds*beta_M_Sqrd*BL_p%b2, 1.0_ReKi, xd%X4_minus1(i,j), BL_p%A2*(KC%q_f_cur - q_f_minus1), 0.0_ReKi ) ! Eqn 1.16a [2]
      
      KC%Cn_q_circ       = KC%C_nalpha_circ*KC%q_f_cur/2.0 - KC%X3 - KC%X4                                                    ! Eqn 1.16
      ! TODO: Why does the Cn_q_circ appear in the following equation GJH 2/28/2017
      KC%Cn_alpha_q_circ = KC%Cn_alpha_q_circ + KC%Cn_q_circ                                                                  ! add to Eqn 1.13
   else ! these aren't used (they are possibly output to UA_OUT file, though)
      KC%X3              = 0.0_ReKi
      KC%X4              = 0.0_ReKi
      KC%Cn_q_circ       = 0.0_ReKi
   end if
   
   K3prime_q       = Get_ExpEqn( BL_p%b5*beta_M_Sqrd*ds, 1.0_ReKi, xd%K3prime_q_minus1(i,j),  BL_p%A5*(KC%q_f_cur - q_f_minus1), 0.0_ReKi )  ! Eqn 1.26
   KC%Cm_q_circ    = -BL_p%C_nalpha*(KC%q_f_cur - K3prime_q)*p%c(i,j)/(16.0_ReKi*beta_M*u%U)                                  ! Eqn 1.25
   
   KC%Cn_pot       = KC%Cn_alpha_q_circ + KC%Cn_alpha_q_nc                                                                    ! Eqn 1.20 [2a]
   
   k_mq            = 7.0_ReKi / (15.0_ReKi*(1.0_ReKi-M) + 1.5_ReKi * BL_p%C_nalpha * BL_p%A5 * BL_p%b5 * beta_M * M**2)       ! Eqn 1.29 [2]      ! CHECK THAT DENOM ISN'T ZERO!
   Kprimeprime_q   = Get_ExpEqn( real(p%dt,ReKi), k_mq**2*T_I , xd%Kprimeprime_q_minus1(i,j) ,  KC%Kq_f , Kq_f_minus1  )      ! Eqn 1.29 [3]
   
      ! Compute Cm_q_nc 
   if ( p%UAMod == UA_MinemmaPierce ) then
      KC%Cm_q_nc =  -1.0_ReKi * KC%Cn_q_nc / 4.0_ReKi - (KC%k_alpha**2) * T_I * (KC%Kq_f - Kprimeprime_q) / (3.0_ReKi*M)      ! Eqn 1.31
   else  
      KC%Cm_q_nc = -7.0_ReKi * (k_mq**2) * T_I * (KC%Kq_f - Kprimeprime_q) / (12.0_ReKi*M)                                    ! Eqn 1.29 [1]       
   end if
   
   !!! THIS IS A PROBLEM IF KC%alpha_e+BL_p%alpha0 ARE NEAR +/-PI/2
   KC%Cc_pot          = KC%Cn_alpha_q_circ * tan(KC%alpha_e+BL_p%alpha0)                                                      ! Eqn 1.21 with cn_pot_circ=KC%Cn_alpha_q_circ as from Eqn 1.20 [3]  
   
   if (OtherState%FirstPass(i,j)) then
      Cn_pot_minus1 = KC%Cn_pot
   else
      Cn_pot_minus1 = xd%Cn_pot_minus1(i,j)
   end if
   
   KC%Dp            = Get_ExpEqn( ds, BL_p%T_p, xd%Dp_minus1(i,j), KC%Cn_pot, Cn_pot_minus1 )                                 ! Eqn 1.35b   
   KC%Cn_prime      = KC%Cn_Pot - KC%Dp                                                                                       ! Eqn 1.35a
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! This code is taken from ADv14 but doesn't reflect the original intent of the UA theory document
#ifdef TEST_THEORY
   if (OtherState%FirstPass(i,j)) then
      Cn_prime_diff = 0.0_ReKi
   else
      Cn_prime_diff = KC%Cn_prime - xd%Cn_prime_minus1(i,j)
   end if
      
   IF ( KC%alpha_filt_cur * Cn_prime_diff < 0. ) THEN

      KC%T_f   = BL_p%T_f0*1.5
   ELSE

      KC%T_f   = BL_p%T_f0
   ENDIF
#endif   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   KC%alpha_f       = KC%Cn_prime / KC%C_nalpha_circ + BL_p%alpha0                                                            ! Eqn 1.34
   
   if (p%flookup) then
      KC%fprime = Get_f_from_Lookup( p%UAMod, u%Re, KC%alpha_f, BL_p%alpha0, KC%C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2)     ! Solve Eqn 1.32a for f when alpha is replaced with alpha_f (see issue when KC%C_nalpha_circ is 0) 
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   else   
      KC%fprime = Get_f( KC%alpha_f, BL_p%alpha0, BL_p%alpha1, BL_p%alpha2, BL_p%S1, BL_p%S2, BL_p%S3, BL_p%S4)               ! Eqn 1.33
   end if
   
   if (OtherState%FirstPass(i,j)) then
      KC%Df = 0.0_ReKi
   else
      KC%Df = Get_ExpEqn( ds, KC%T_f, xd%Df_minus1(i,j), KC%fprime, xd%fprime_minus1(i,j) )                                   ! Eqn 1.36b
   end if
      
   KC%fprimeprime   = KC%fprime - KC%Df                                                                                       ! Eqn 1.36a
   
   if (p%Flookup) then
         ! Compute fprime using Eqn 1.32 and Eqn 1.33
      KC%fprime_c   = Get_f_c_from_Lookup( u%Re, KC%alpha_f, BL_p%alpha0, KC%C_nalpha_circ, BL_p%eta_e, AFInfo, ErrStat2, ErrMsg2) ! Solve Eqn 1.32b for f when alpha is replaced with alpha_f
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
     
      if (OtherState%FirstPass(i,j)) then
         KC%Df_c = 0.0_ReKi
      else
         KC%Df_c = Get_ExpEqn( ds, KC%T_f, xd%Df_c_minus1(i,j), KC%fprime_c, xd%fprime_c_minus1(i,j)  )
      end if
   
         ! Compute Df using Eqn 1.36b   
   
         ! Compute fprimeprime using Eqn 1.36a
      KC%fprimeprime_c   = KC%fprime_c - KC%Df_c
   else
      KC%fprime_c = KC%fprime
      KC%Df_c = KC%Df
      KC%fprimeprime_c = KC%fprimeprime
   end if
   
    
   if ( p%UAMod == UA_Gonzalez ) then
      KC%Cn_FS   = KC%Cn_alpha_q_nc + KC%Cn_alpha_q_circ *  ( (1.0_ReKi + 2.0_ReKi*sqrt(KC%fprimeprime) ) / 3.0_ReKi )**2     ! Eqn 1.39 [bjj: note that KC%Cn_alpha_q_circ doesn't match the equation 1.39 for the Gonzales model]
   else
      KC%Cn_FS   = KC%Cn_alpha_q_nc + KC%Cn_alpha_q_circ *  ( (1.0_ReKi +          sqrt(KC%fprimeprime) ) / 2.0_ReKi )**2     ! Eqn 1.38   
   end if
   
      
   if ( p%UAMod == UA_MinemmaPierce ) then
      if (OtherState%FirstPass(i,j)) then     
         KC%Dalphaf    = 0.0_ReKi
      else
         KC%Dalphaf    = Get_ExpEqn( ds, 0.1_ReKi*KC%T_f, xd%Dalphaf_minus1(i,j), KC%alpha_f, xd%alphaf_minus1(i,j) )         ! Eqn 1.43
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
      if (xd%tau_V(i,j) > BL_p%T_VL .AND. ( (KC%alpha_filt_cur - BL_p%alpha0) * KC%Kalpha_f) > 0) then 
            ! The assertion is the T_V will always equal BL_p%T_V0/2 when this condition is satisfied
         KC%Cn_v = xd%Cn_v_minus1(i,j)*exp(-ds/KC%T_V)                                                                        ! Eqn 1.52    
      else      
         KC%Cn_v = Get_ExpEqn( ds, KC%T_V, xd%Cn_v_minus1(i,j), KC%C_V, xd%C_V_minus1(i,j) )                                  ! Eqn 1.47
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
   
   p%c          = InitInp%c         ! this can't be 0
   p%numBlades  = InitInp%numBlades
   p%nNodesPerBlade  = InitInp%nNodesPerBlade
   p%UAMod      = InitInp%UAMod    
   p%a_s        = InitInp%a_s ! this can't be 0
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
   call AllocAry(m%T_Sh     ,p%nNodesPerBlade,p%numBlades,'m%T_Sh',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
#endif  
   
   if (ErrStat >= AbortErrLev) return
   
   
   call UA_ReInit( p, xd, OtherState, m )   ! initializes values of states and misc vars

   
end subroutine UA_InitStates_Misc 
!==============================================================================   
subroutine UA_ReInit( p, xd, OtherState, m )  
   type(UA_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UA_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(UA_OtherStateType),      intent(inout)  :: OtherState  ! Initial other states
   type(UA_MiscVarType),         intent(inout)  :: m           ! Initial misc/optimization variables

   m%FirstWarn_M = .true.   
   
   OtherState%sigma1    = 1.0_ReKi
   OtherState%sigma3    = 1.0_ReKi
   
#ifdef UA_OUTS
   m%TESF      = .FALSE.  
   m%LESF      = .FALSE.   
   m%VRTX      = .FALSE. 
   m%T_sh      = 0.0_ReKi 
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

end subroutine UA_ReInit
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
   
#ifdef UA_OUTS   
   integer(IntKi)                               :: i,j, iNode, iOffset
   character(64)                                :: chanPrefix
#endif   
   
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
      
#ifdef UA_OUTS   

      ! Allocate and set the InitOut data
   p%NumOuts = 42
      
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
         InitOut%WriteOutputHdr(iOffset+ 1)  = 'ALPHA_filt'//chanPrefix      
         InitOut%WriteOutputHdr(iOffset+ 2)  = 'VREL'//chanPrefix  
         InitOut%WriteOutputHdr(iOffset+ 3)  = 'Cn'//chanPrefix      
         InitOut%WriteOutputHdr(iOffset+ 4)  = 'Cc'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+ 5)  = 'Cl'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+ 6)  = 'Cd'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+ 7)  = 'Cm'//chanPrefix     
         InitOut%WriteOutputHdr(iOffset+ 8)  = 'Cn_aq_circ'//chanPrefix    
         InitOut%WriteOutputHdr(iOffset+ 9)  = 'Cn_aq_nc'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+10)  = 'Cn_pot'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+11)  = 'Dp'//chanPrefix       
         InitOut%WriteOutputHdr(iOffset+12)  = 'Cn_prime'//chanPrefix        
         InitOut%WriteOutputHdr(iOffset+13)  = 'fprime'//chanPrefix       
         InitOut%WriteOutputHdr(iOffset+14)  = 'Df'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+15)  = 'Cn_v'//chanPrefix         
         InitOut%WriteOutputHdr(iOffset+16)  = 'Tau_V'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+17)  = 'LESF'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+18)  = 'TESF'//chanPrefix 
         InitOut%WriteOutputHdr(iOffset+19)  = 'VRTX'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+20)  = 'C_v'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+21)  = 'Cm_a_nc'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+22)  = 'Cm_q_nc'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+23)  = 'Cm_v'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+24)  = 'alpha_p_f'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+25)  = 'Dalphaf'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+26)  = 'PMC'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+27)  = 'T_f'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+28)  = 'T_V'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+29)  = 'dS'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+30)  = 'T_alpha'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+31)  = 'T_q'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+32)  = 'k_alpha'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+33)  = 'k_q'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+34)  = 'alpha_e'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+35)  = 'X1'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+36)  = 'X2'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+37)  = 'cn_q_nc'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+38)  = 'alpha_f'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+39)  = 'fprimeprime'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+40)  = 'sigma1'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+41)  = 'sigma3'//chanPrefix
         InitOut%WriteOutputHdr(iOffset+42)  = 'T_sh'//chanPrefix                  
                           
         
         InitOut%WriteOutputUnt(iOffset+1)  ='(deg)'                                                
         InitOut%WriteOutputUnt(iOffset+2)  ='(m/s)'                                                  
         InitOut%WriteOutputUnt(iOffset+3)  ='(-)'                                                    
         InitOut%WriteOutputUnt(iOffset+4)  ='(-)'                                                   
         InitOut%WriteOutputUnt(iOffset+5)  ='(-)'                                                    
         InitOut%WriteOutputUnt(iOffset+6)  ='(-)'                                                   
         InitOut%WriteOutputUnt(iOffset+7)  ='(-)'                                                    
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
         InitOut%WriteOutputUnt(iOffset+24) ='(-)'
         InitOut%WriteOutputUnt(iOffset+25) ='(-)'
         InitOut%WriteOutputUnt(iOffset+26) ='(-)'
         InitOut%WriteOutputUnt(iOffset+27) ='(-)'
         InitOut%WriteOutputUnt(iOffset+28) ='(-)'
         InitOut%WriteOutputUnt(iOffset+29) ='(-)'
         InitOut%WriteOutputUnt(iOffset+30) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+31) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+32) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+33) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+34) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+35) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+36) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+37) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+38) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+39) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+40) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+41) ='(-)' 
         InitOut%WriteOutputUnt(iOffset+42) ='(-)' 
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
      TESF = KC%fprimeprime < xd%fprimeprime_minus1(i,j) ! Separation point is moving towards the Leading Edge when .true.; otherwise separation point is moving toward trailing edge   
            
      ! Process VRTX-related quantities
      !!!!!!!!!!!!!!!!!!!!!
      !! NEW CODE 2/19/2015
      VRTX = (xd%tau_V(i,j) <= 2.0_ReKi*BL_p%T_VL) .and. (xd%tau_V(i,j) > 0.0_ReKi)
      
      
      if (( xd%tau_V(i,j) >= (BL_p%T_VL + T_sh) ) .and. TESF)  then !.and. (TESF) RRD added 
         xd%tau_V(i,j)     = 0.0_ReKi
        ! LESF=.FALSE. !also added this
      end if
      
!bjj: update sigma1 to value at t + dt      
      OtherState%sigma1(i,j) = 1.0_ReKi
      Kafactor      = KC%Kalpha_f*KC%dalpha0
     ! if ( KC%fprimeprime < 0.7 .AND.  KC%fprimeprime > 0.3 ) then 
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

            if ( VRTX .and. (xd%tau_V(i,j) <= BL_p%T_VL) ) then  ! Still shedding a vortex?
               OtherState%sigma1(i,j) = 0.25_ReKi
            end if
            if (Kafactor > 0.0_ReKi) then
               OtherState%sigma1(i,j) = 0.75_ReKi
            end if
            
         end if
      
!bjj: update sigma3 to value at t + dt      
         
      OtherState%sigma3(i,j) = 1.0_ReKi
      
      if ( (xd%tau_V(i,j) <= 2.0_ReKi*BL_p%T_VL) .and. (xd%tau_V(i,j) >= BL_p%T_VL) ) then
         OtherState%sigma3(i,j) = 3.0_ReKi
         if (.not. TESF) then
            OtherState%sigma3(i,j) = 4.0_ReKi
            if ( VRTX .and. (xd%tau_V(i,j) <= BL_p%T_VL) ) then
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
      
      if ((.not. TESF) .and. (KC%Kq_f*KC%dalpha0 < 0.0_ReKi)) then
         OtherState%sigma3(i,j) = 1.0_ReKi
      end if
   
      
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
      xd%fprimeprime_minus1(i,j)  = KC%fprimeprime
      xd%Df_minus1(i,j)           = KC%Df
      if (p%Flookup) then
         xd%Df_c_minus1(i,j)           = KC%Df_c
         xd%fprimeprime_c_minus1(i,j)  = KC%fprimeprime_c
         xd%fprime_c_minus1(i,j)       = KC%fprime_c
      end if
      
      xd%Dalphaf_minus1(i,j)      = KC%Dalphaf
      xd%fprime_minus1(i,j)       = KC%fprime
      xd%alphaf_minus1(i,j)       = KC%alpha_f
      xd%Cn_v_minus1(i,j)         = KC%Cn_v
      xd%C_V_minus1(i,j)          = KC%C_V
      OtherState%FirstPass(i,j)   = .false.
   
      if ( xd%tau_V(i,j) > 0.0 .or. LESF ) then                !! TODO:  Verify this condition 2/20/2015 GJH
         xd%tau_V(i,j)          = xd%tau_V(i,j) + 2.0_ReKi*p%dt*u%U / p%c(i,j)     ! Eqn 1.51
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
   
   type(AFI_UA_BL_Type)                                   :: BL_p  ! airfoil values computed in Kelvin Chain
   type(UA_KelvinChainType)                               :: KC    ! values computed in Kelvin Chain
   
   real(ReKi)                                             :: Cm_FS  
   real(ReKi)                                             :: Cc_FS
   real(ReKi)                                             :: Cl_temp, Cd_temp, Cm_temp, Cn_temp, Cm_alpha_nc
   real(ReKi)                                             :: M, f, k2_hat
   real(ReKi)                                             :: Cm_v, alpha_prime_f
   real(ReKi)                                             :: x_cp_hat                      ! center-of-pressure distance from LE in chord fraction
   real(ReKi)                                             :: Cm_common                     ! 

   integer(intKi)                                         :: ncols
   
#ifdef UA_OUTS
   integer                                                :: iOffset
#endif   
   
   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""
   
   Cm_alpha_nc = 0.0_ReKi
   Cm_temp     = 0.0_ReKi
   
   
   if (OtherState%FirstPass(misc%iBladeNode, misc%iBlade) .or. EqualRealNos(u%U, 0.0_ReKi) ) then ! note: if u%U is = in UpdateStates, BEMT shuts off UA; however, it could still be called with u%U=0 here
            
      call GetSteadyOutputs(AFInfo, u%alpha, y%Cl, y%Cd, y%Cm, BL_p%Cd0, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
      
      
      y%Cn = y%Cl*cos(u%alpha) + (y%Cd-BL_p%Cd0)*sin(u%alpha)
      y%Cc = y%Cl*sin(u%alpha) - (y%Cd-BL_p%Cd0)*cos(u%alpha)
      
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
         alpha_prime_f     = 0.0_ReKi
      KC%Dalphaf           = 0.0_ReKi
      KC%T_f               = 0.0_ReKi
      KC%T_V               = 0.0_ReKi
      KC%alpha_filt_cur    = 0.0_ReKi
      
   else
      
      M           = u%U / p%a_s
      
      call UA_CheckMachNumber(M, misc%FirstWarn_M, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
            
      call ComputeKelvinChain( misc%iBladeNode, misc%iBlade, u, p, xd, OtherState, misc, AFInfo, KC, BL_p, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      

      !.............................
      ! Cn
      !.............................
         ! Eqn 1.53 or 1.53b depending on UAMod
      if (xd%tau_v(misc%iBladeNode, misc%iBlade) > 0.0) then
         y%Cn= KC%Cn_FS + KC%Cn_v                                                                                                                  ! Eqn 1.53
      else
         y%Cn= KC%Cn_FS
      end if

      !.............................
      ! Cc
      !.............................
      
      
      if ( p%flookup ) then 
      !if ( p%UAMod == UA_Gonzalez ) then
         Cc_FS = BL_p%eta_e*KC%Cc_pot *(sqrt(KC%fprimeprime_c) - Gonzales_factor)                                                                  ! Eqn 1.40  
      else
         Cc_FS = BL_p%eta_e*KC%Cc_pot * sqrt(KC%fprimeprime_c)                                                                                     ! Eqn 1.40 without Gonzalez's modification for negative Cc
      end if
      
            
      if ( p%UAMod == UA_MinemmaPierce ) then
#ifdef TEST_THEORY
            y%Cc = Cc_FS + KC%Cn_v*tan(KC%alpha_e)*(1-xd%tau_v(misc%iBladeNode, misc%iBlade)/(BL_p%T_VL))                                          ! Eqn 1.55 with Eqn. 1.40
#else            
            y%Cc = Cc_FS + KC%Cn_v*    KC%alpha_e *(1.0_ReKi-xd%tau_v(misc%iBladeNode, misc%iBlade)/(BL_p%T_VL))                                   ! Eqn 1.55 with approximation of tan(KC%alpha_e)=KC%alpha_e and Eqn. 1.40 substitution
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
              !    f      = Get_f_from_Lookup( p%UAMod, u%Re, KC%alpha_filt_cur, BL_p%alpha0, KC%C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2)
              ! else   
               ! TODO: Need to understand the effect of the offset on this f in the following equations and fprimeprime_c.
                  !bjj: fprimeprime_c is computed with the Gonzalez offset in the Kelvin Chain, so it's not just f that could be a problem.
                  f      = Get_f_c_from_Lookup( u%Re, KC%alpha_filt_cur, BL_p%alpha0, KC%C_nalpha_circ, BL_p%eta_e, AFInfo, ErrStat2, ErrMsg2)
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

      
      !.............................
      ! convert cm
      !.............................
      
         ! Check for Cm column in AFInfo data        
      ncols = size(AFInfo%Table(1)%Coefs,2)
      if (ncols < 3) then ! we don't have a cm column, so make everything 0
         
         y%Cm          = 0.0_ReKi
         Cm_v          = 0.0_ReKi
         alpha_prime_f = 0.0_ReKi
     
      else
            
         if ( p%UAMod == UA_Gonzalez ) then !bjj: what is this doing?
            
            if (abs(KC%alpha_f - BL_p%alpha0) < .01) then
               if (KC%alpha_f < BL_p%alpha0) then
                  KC%alpha_f = KC%alpha_f - .01
               else
                  KC%alpha_f = KC%alpha_f + .01
               end if
            end if
         
            call GetSteadyOutputs(AFInfo, KC%alpha_f, Cl_temp, Cd_temp, Cm_temp, BL_p%Cd0, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            Cn_temp = Cl_temp*cos(KC%alpha_f) + (Cd_temp-BL_p%Cd0)*sin(KC%alpha_f)
            
            ! TODO: In the future we are going to look at a delayed version of fprime_m and that is why we are using
            ! the variable name fprimeprime_m even though we have not introduced the delay, yet.  GJH 4/12/2016
            
!#ifndef CHECK_DENOM_DIFF
            if (abs(Cn_temp) < 0.01 ) then
!#else
!            if (EqualRealNos(Cn_temp,0.0_ReKi)) then
!#endif   
               KC%fprimeprime_m = 1.0
            else            
               KC%fprimeprime_m = (Cm_temp - BL_p%Cm0) / Cn_temp 
            end if
            
            Cm_temp = 0.0_ReKi !reset this to zero, only for output purposes

         
         else
            KC%fprimeprime_m = 0.0
         end if
      
         !...........
         ! compute Cm_FS
         !...........
         
         Cm_alpha_nc = - KC%Cn_alpha_nc / 4.0_ReKi                                                                                                 ! Eqn 1.27
         Cm_common   = KC%Cm_q_circ + Cm_alpha_nc + KC%Cm_q_nc                                                                                     ! second parts of Eqn 1.41, Eqn 1.44, and Eqn 1.45
   
         if ( p%UAMod == UA_Baseline ) then
            
            x_cp_hat = BL_p%k0 + BL_p%k1*(1.0_ReKi-KC%fprimeprime) + BL_p%k2*sin(pi*KC%fprimeprime**BL_p%k3)                                       ! Eqn 1.42
            Cm_FS  = BL_p%Cm0 - KC%Cn_alpha_q_circ*(x_cp_hat - 0.25_ReKi) + Cm_common                                                              ! Eqn 1.41

         elseif ( p%UAMod == UA_MinemmaPierce ) then
      
               ! Look up Cm using alpha_prime_f
            alpha_prime_f = KC%alpha_f - KC%Dalphaf                                                                                                ! Eqn 1.43a
            call GetSteadyOutputs(AFInfo, alpha_prime_f, Cl_temp, Cd_temp, Cm_temp, BL_p%Cd0, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)            
            Cm_FS = Cm_temp + Cm_common                                                                                                            ! Eqn 1.44
      
         else ! UAMod == UA_Gonzalez
            
            Cm_FS = BL_p%Cm0 + KC%Cn_FS*KC%fprimeprime_m + Cm_common                                                                               ! Eqn 1.45 (NOTE: This differs from manual in that Cm0 is added because it is subtracted in our version of fprimeprime_m
            alpha_prime_f = 0.0_ReKi
         end if   
      
         Cm_v     = -BL_p%x_cp_bar*( 1.0_ReKi - cos( pi*xd%tau_v(misc%iBladeNode, misc%iBlade)/BL_p%T_VL ) )*KC%Cn_v                               ! Eqn 1.57   
         y%Cm   = Cm_FS + Cm_v                                                                                                                     ! Eqn 1.58, Eqn 1.59, and Eqn 1.60
         
      end if
   end if
   
#ifdef UA_OUTS
   iOffset = (misc%iBladeNode-1)*p%NumOuts + (misc%iBlade-1)*p%nNodesPerBlade*p%NumOuts 
   if (allocated(y%WriteOutput)) then  !bjj: because BEMT uses local variables for UA output, y%WriteOutput is not necessarially allocated. Need to figure out a better solution.
      y%WriteOutput(iOffset+ 1)    = KC%alpha_filt_cur*R2D                                                    
      y%WriteOutput(iOffset+ 2)    = u%U                                                                 
      y%WriteOutput(iOffset+ 3)    = y%Cn                                                                
      y%WriteOutput(iOffset+ 4)    = y%Cc                                                                
      y%WriteOutput(iOffset+ 5)    = y%Cl                                                                
      y%WriteOutput(iOffset+ 6)    = y%Cd                                                                
      y%WriteOutput(iOffset+ 7)    = y%Cm                                                                
      y%WriteOutput(iOffset+ 8)    = KC%Cn_alpha_q_circ               ! CNCP in ADv14                                      
      y%WriteOutput(iOffset+ 9)    = KC%Cn_alpha_q_nc                 ! CNIQ in ADv14                                      
      y%WriteOutput(iOffset+10)    = KC%Cn_pot                                                              
      y%WriteOutput(iOffset+11)    = KC%Dp                                                                  
      y%WriteOutput(iOffset+12)    = KC%Cn_prime                                                            
      y%WriteOutput(iOffset+13)    = KC%fprime                                                              
      y%WriteOutput(iOffset+14)    = KC%Df                                                                  
      y%WriteOutput(iOffset+15)    = KC%Cn_V                                                                
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
      y%WriteOutput(iOffset+20)    = KC%C_V  
      y%WriteOutput(iOffset+21)    = Cm_alpha_nc
      y%WriteOutput(iOffset+22)    = KC%Cm_q_nc 
      y%WriteOutput(iOffset+23)    = Cm_v
      y%WriteOutput(iOffset+24)    = alpha_prime_f
      y%WriteOutput(iOffset+25)    = KC%Dalphaf
      y%WriteOutput(iOffset+26)    = Cm_temp
      y%WriteOutput(iOffset+27)    = KC%T_f
      y%WriteOutput(iOffset+28)    = KC%T_V
      y%WriteOutput(iOffset+29)    = 2.0_ReKi*u%U*p%dt/p%c(misc%iBladeNode, misc%iBlade)  ! Eqn 1.51 (tau_v) 
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
      y%WriteOutput(iOffset+40)    = OtherState%sigma1(misc%iBladeNode, misc%iBlade)
      y%WriteOutput(iOffset+41)    = OtherState%sigma3(misc%iBladeNode, misc%iBlade)
      y%WriteOutput(iOffset+42)    = misc%T_sh(misc%iBladeNode, misc%iBlade)
   end if
#endif
   
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
