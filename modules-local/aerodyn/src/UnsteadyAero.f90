module UnsteadyAero
use NWTC_Library   
use UnsteadyAero_Types
use AirfoilInfo
implicit none 

private

!INTEGER, PARAMETER              :: R8Ki     = SELECTED_REAL_KIND( 14, 300 )     ! Kind for eight-byte floating-point numbers
!INTEGER, PARAMETER              :: SiKi     = SELECTED_REAL_KIND(  6,  30 )     ! Kind for four-byte, floating-point numbers
!INTEGER, PARAMETER              :: ReKi     = SiKi                              ! Default kind for floating-point numbers
!INTEGER, PARAMETER              :: DbKi     = R8Ki                              ! Default kind for double floating-point numbers


public :: UnsteadyAero_Init
public :: UnsteadyAero_UpdateDiscState
public :: UnsteadyAero_UpdateStates
public :: UnsteadyAero_CalcOutput


   contains

   
subroutine GetSteadyOutputs(AFInfo, AOA, Re, Cl, Cd, ErrStat, ErrMsg)
   
      type(AFInfoType), intent(in   ) :: AFInfo
      real(ReKi),       intent(in   ) :: AOA
      real(ReKi),       intent(in   ) :: Re
      real(ReKi),       intent(  out) :: Cl
      real(ReKi),       intent(  out) :: Cd
      integer(IntKi),   intent(  out) :: ErrStat     ! Error status of the operation
      character(*),     intent(  out) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
      real                            :: IntAFCoefs(4)                ! The interpolated airfoil coefficients.
      integer                         :: s1
      ErrStat = ErrID_None
      ErrMsg  = ''
   
      ! TODO: Extend this to use the UnsteadyAero module to determine the Cl, Cd, Cm info, as needed.  We may need to be tracking whether this call is
      ! part of an UpdateStates action or a CalcOutput action, because the calculation chain differs for the two.
      
      
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
      s1 = size(AFInfo%Table(1)%Coefs,2)
   
      IntAFCoefs(1:s1) = CubicSplineInterpM( 1.0*real( AOA ) &
                                              , AFInfo%Table(1)%Alpha &
                                              , AFInfo%Table(1)%Coefs &
                                              , AFInfo%Table(1)%SplineCoefs &
                                              , ErrStat, ErrMsg )
   
      Cl = IntAFCoefs(1)
      Cd = IntAFCoefs(2)
   
   end subroutine GetSteadyOutputs
real(ReKi) function Get_ds( U, c, dt )
   ! Called by : Get_Alpha_e, UnsteadyAero_UpdateStates
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: U
   real(ReKi), intent(in   ) :: c
   real(DbKi), intent(in   ) :: dt
   Get_ds = 2.0_ReKi*U*dt/c
end function Get_ds

real(ReKi) function Get_Kupper( dt, x, x_minus1 )
   real(DbKi), intent(in   ) :: dt
   real(ReKi), intent(in   ) :: x
   real(ReKi), intent(in   ) :: x_minus1
   
   Get_Kupper = (x - x_minus1) / dt
   
end function Get_Kupper

real(ReKi) function Get_X(ds, betaSqrd, b, A, alpha, alpha_minus1, X_minus1 )
   real(ReKi), intent(in   ) :: ds
   real(ReKi), intent(in   ) :: betaSqrd
   real(ReKi), intent(in   ) :: b
   real(ReKi), intent(in   ) :: A
   real(ReKi), intent(in   ) :: alpha
   real(ReKi), intent(in   ) :: alpha_minus1
   real(ReKi), intent(in   ) :: X_minus1
   
   Get_X = X_minus1 * exp(-b*betaSqrd*ds) + A*exp(-0.5_ReKi*b*betaSqrd*ds)*(alpha - alpha_minus1)

end function Get_X

real(ReKi) function Get_Beta_M( M ) 
   ! Called by : NONE
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: M
   Get_Beta_M = sqrt(1 - M**2)
end function Get_Beta_M

real(ReKi) function Get_Beta_M_Sqrd( M ) 
   ! Called by :  UnsteadyAero_UpdateStates
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: M
   Get_Beta_M_Sqrd = 1 - M**2
end function Get_Beta_M_Sqrd

real(ReKi) function Get_k_( M, beta, C_nalpha          )
   ! Called by : UnsteadyAero_UpdateStates
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: M                               ! Mach number
   real(ReKi), intent(in   ) :: beta
   real(ReKi), intent(in   ) :: C_nalpha                        !  
   
   ! Implements equation 1.11
   Get_k_   = 1.0_ReKi / ( (1.0_ReKi - M) + C_nalpha * M**2 * beta *  0.413_ReKi  )
   
end function Get_k_

real(ReKi) function Get_ExpEqn( dt, T, Y_minus1, x, x_minus1 )
   real(ReKi), intent(in   ) :: dt
   real(ReKi), intent(in   ) :: T
   real(ReKi), intent(in   ) :: Y_minus1
   real(ReKi), intent(in   ) :: x
   real(ReKi), intent(in   ) :: x_minus1
   
   Get_ExpEqn = Y_minus1*exp(-dt/T) + (x - x_minus1)*exp(-0.5_ReKi*dt/T)

end function Get_ExpEqn


real(ReKi) function Get_T_Alpha( M, a_s, c, C_nalpha, k_alpha ) 
   ! Called by : UnsteadyAero_UpdateStates
   ! Calls  to : NONE

   real(ReKi), intent(in   ) :: M                               ! Mach number
   real(ReKi), intent(in   ) :: a_s                             ! speed of sound (m/s)
   real(ReKi), intent(in   ) :: c                               ! chord length (m)
   real(ReKi), intent(in   ) :: C_nalpha                        ! 
   real(ReKi), intent(in   ) :: k_alpha
   

   real(ReKi)                :: T_I
   
   ! Implements equation 1.10a
   
   T_I           = c / a_s
   Get_T_Alpha   = 0.75*k_alpha * T_I
   
end function Get_T_Alpha   


real(ReKi) function Get_T_q(M, beta_M_sqrd, T_I, a_s, c, C_nalpha) 
   ! Called by : UnsteadyAero_UpdateStates
   ! Calls  to : NONE

   real(ReKi), intent(in   ) :: M                               ! Mach number
   real(ReKi), intent(in   ) :: beta_M_sqrd
   real(ReKi), intent(in   ) :: a_s                             ! speed of sound (m/s)
   real(ReKi), intent(in   ) :: c                               ! chord length (m)
   real(ReKi), intent(in   ) :: C_nalpha                        ! 
   
   
   real(ReKi)                :: k_q
   real(ReKi)                :: T_I
   
   ! Implements equation 1.10b
   
   
   k_q           = 1.0_ReKi / ( (1.0_ReKi - M) + C_nalpha * M**2 * beta_M_sqrd *  0.413_ReKi  )
   T_I           = c / a_s
   Get_T_q       = k_q * T_I
   
end function Get_T_q   


   
real(ReKi) function Get_Pitchrate( dt, alpha_cur, alpha_minus1, U, c )
   ! this is the non-dimensional pitching rate, given as q in the documentation: equation 1.8
   !
   ! Called by : UnsteadyAero_UpdateStates
   ! Calls  to : NONE
   real(DbKi), intent(in   ) :: dt
   real(ReKi), intent(in   ) :: alpha_cur                         ! angle of attack, current timestep  (radians)
   real(ReKi), intent(in   ) :: alpha_minus1                        ! angle of attack, previous timestep  (radians)
   real(ReKi), intent(in   ) :: U                               ! air velocity magnitude relative to the airfoil (m/s)
   real(ReKi), intent(in   ) :: c                               ! chord length (m)

      ! Implements equation 1.8
   Get_Pitchrate = ( alpha_cur - alpha_minus1 ) * c / (U*dt)  
   
end function Get_Pitchrate

real(ReKi) function Get_Cn_nc( M, T, K, Kprime )

   ! Called by : UnsteadyAero_UpdateStates
   ! Calls  to : NONE
   real(ReKi)                :: M
   real(ReKi), intent(in   ) :: T                               ! mach-dependent time constant related to alpha and non-circulatory terms (???)
   real(ReKi)                :: K
   real(ReKi)                :: Kprime
   
   ! Implements equation 1.xx
   
   Get_Cn_nc =  T * ( K - Kprime ) / M
   
end function Get_Cn_nc 

real(ReKi) function Get_Cn_q_nc( T_q, q_cur, q_minus1, q_minus2, Kprima_q_minus1, dt )
   ! Called by : Get_Cn_Alpha_q_nc
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: T_q                             ! mach-dependent time constant related to alpha and non-circulatory terms (???)
   real(ReKi), intent(in   ) :: q_cur                           ! non-dimensional pitching rate, current timestep   (radians)
   real(ReKi), intent(in   ) :: q_minus1                          ! non-dimensional pitching rate, previous timestep  (radians)
   real(ReKi), intent(in   ) :: q_minus2                        ! non-dimensional pitching rate, two timesteps ago  (radians)
   real(ReKi), intent(in   ) :: Kprima_q_minus1
   real(DbKi), intent(in   ) :: dt                              ! time between simulation steps [note, the same dt is assumed between all q values] (sec)
   
   real(ReKi)                :: K_q
   real(ReKi)                :: K_q_minus1
   real(ReKi)                :: Kprime_q
   
   ! Implements equation 1.19
   K_q         = ( q_cur  - q_minus1   ) / dt
   K_q_minus1    = ( q_minus1 - q_minus2 ) / dt
   Kprime_q    = Kprima_q_minus1 * exp( -dt/T_q ) + ( K_q - K_q_minus1 ) *  exp( -dt / ( 2.0_ReKi*T_q ) )
   Get_Cn_q_nc = 4 * T_q * ( K_q - Kprime_q )
   
end function Get_Cn_q_nc 

!real(ReKi) function Get_Cn_Alpha_q_nc( T_alpha, alpha_cur, alpha_minus1, alpha_minus2, Kprime_alpha_minus1, Kprime_q_minus1, T_q, q_cur, q_minus1, q_minus2, dt )
!   ! Called by : Get_Cn_pot, Get_Cn
!   ! Calls  to : Get_Cn_Alpha_nc, Get_Cn_q_nc
!   real(ReKi), intent(in   ) :: T_alpha                         ! mach-dependent time constant related to alpha and non-circulatory terms (???)
!   real(ReKi), intent(in   ) :: alpha_cur                         ! angle of attack, current timestep   (radians)
!   real(ReKi), intent(in   ) :: alpha_minus1                        ! angle of attack, previous timestep  (radians)
!   real(ReKi), intent(in   ) :: alpha_minus2                      ! angle of attack, two timesteps ago  (radians)
!   real(ReKi), intent(in   ) :: Kprime_alpha_minus1
!   real(ReKi), intent(in   ) :: Kprime_q_minus1
!   real(ReKi), intent(in   ) :: T_q                             ! mach-dependent time constant related to alpha and non-circulatory terms (???)
!   real(ReKi), intent(in   ) :: q_cur                           ! non-dimensional pitching rate, current timestep   (radians)
!   real(ReKi), intent(in   ) :: q_minus1                          ! non-dimensional pitching rate, previous timestep  (radians)
!   real(ReKi), intent(in   ) :: q_minus2                        ! non-dimensional pitching rate, two timesteps ago  (radians)
!   real(DbKi), intent(in   ) :: dt                              ! time between simulation steps [note, the same dt is assumed between all q values] (sec)
!   
!   ! Implements equation 1.17
!   Get_Cn_Alpha_q_nc = Get_Cn_Alpha_nc( T_alpha, alpha_cur, alpha_minus1, alpha_minus2, Kprime_alpha_minus1, dt ) + Get_Cn_q_nc( T_q, q_cur, q_minus1, q_minus2, Kprime_q_minus1, dt )
!
!end function Get_Cn_Alpha_q_nc



real(ReKi) function Get_alpha_e( alpha, alpha0, X1, X2 )
   ! Called by : Get_Cn_alpha_q_circ
   ! Calls  to : NONE
   ! Implements equation 1.14
   real(ReKi), intent(in   ) :: alpha
   real(ReKi), intent(in   ) :: alpha0
   real(ReKi), intent(in   ) :: X1
   real(ReKi), intent(in   ) :: X2
   
   Get_Alpha_e   = (alpha - alpha0) - X1 - X2
   
end function Get_Alpha_e
   
!real(ReKi) function Get_Cn_alpha_q_circ( dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha0, X_1_minus1, X_2_minus1, b1, b2, A1, A2 )
!   ! Called by : Get_Cn_pot
!   ! Calls  to : Get_Alpha_e
!   ! Implements equation 1.13
!   real(DbKi), intent(in   ) :: dt
!   real(ReKi), intent(in   ) :: ds
!   real(ReKi), intent(in   ) :: U
!   real(ReKi), intent(in   ) :: M
!   real(ReKi), intent(in   ) :: beta_M_sqrd
!   real(ReKi), intent(in   ) :: c
!   real(ReKi), intent(in   ) :: C_nalpha_circ
!   real(ReKi), intent(in   ) :: alpha_cur
!   real(ReKi), intent(in   ) :: alpha_minus1
!   real(ReKi), intent(in   ) :: alpha0
!   real(ReKi), intent(in   ) :: X_1_minus1
!   real(ReKi), intent(in   ) :: X_2_minus1
!   real(ReKi), intent(in   ) :: b1
!   real(ReKi), intent(in   ) :: b2
!   real(ReKi), intent(in   ) :: A1
!   real(ReKi), intent(in   ) :: A2
!   
!   Get_Cn_alpha_q_circ = C_nalpha_circ * Get_Alpha_e( dt, ds, U, M, beta_M_sqrd, c, alpha_cur, alpha_minus1, alpha0, X_1_minus1, X_2_minus1, b1, b2, A1, A2 )
!                                                    
!end function Get_Cn_alpha_q_circ


!real(ReKi) function Get_Cn_pot( dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, &
!                               q_cur, q_minus1, q_minus2, T_alpha, T_q, X_1_minus1, X_2_minus1, Cn_alpha_q_nc, Kprime_alpha_minus1, Kprime_q_minus1, b1, b2, A1, A2 )
!   ! Called by : Get_Cn_prime
!   ! Calls  to : Get_Cn_alpha_q_circ
!   ! Implements equation 1.20
!   real(DbKi), intent(in   ) :: dt
!   real(ReKi), intent(in   ) :: ds
!   real(ReKi), intent(in   ) :: U
!   real(ReKi), intent(in   ) :: M
!   real(ReKi), intent(in   ) :: beta_M_sqrd
!   real(ReKi), intent(in   ) :: c
!   real(ReKi), intent(in   ) :: C_nalpha_circ
!   real(ReKi), intent(in   ) :: alpha_cur
!   real(ReKi), intent(in   ) :: alpha_minus1
!   real(ReKi), intent(in   ) :: alpha_minus2
!   real(ReKi), intent(in   ) :: alpha0
!   real(ReKi), intent(in   ) :: q_cur
!   real(ReKi), intent(in   ) :: q_minus1
!   real(ReKi), intent(in   ) :: q_minus2
!   real(ReKi), intent(in   ) :: T_alpha
!   real(ReKi), intent(in   ) :: T_q
!   real(ReKi), intent(in   ) :: X_1_minus1
!   real(ReKi), intent(in   ) :: X_2_minus1
!   real(ReKi), intent(in   ) :: Cn_alpha_q_nc
!   real(ReKi), intent(in   ) :: Kprime_alpha_minus1
!   real(ReKi), intent(in   ) :: Kprime_q_minus1
!   real(ReKi), intent(in   ) :: b1
!   real(ReKi), intent(in   ) :: b2
!   real(ReKi), intent(in   ) :: A1
!   real(ReKi), intent(in   ) :: A2
!   !          = Eqn 1.13 + Eqn 1.17
!   Get_Cn_pot = Get_Cn_alpha_q_circ(dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha0, X_1_minus1, X_2_minus1, b1, b2, A1, A2) + &
!                Cn_alpha_q_nc
!   
!end function Get_Cn_pot
   
!==============================================================================
! TE Flow Separation Equations                                                !
!==============================================================================

real(ReKi) function Get_f(alpha, alpha1, S1, S2)
   ! Called by : Get_f_prime
   ! Calls  to : NONE
   ! Implements Equation 1.30
   real(ReKi), intent(in   ) :: alpha
   real(ReKi), intent(in   ) :: alpha1
   real(ReKi), intent(in   ) :: S1
   real(ReKi), intent(in   ) :: S2

   real(ReKi)                :: dalpha
   
   dalpha = abs(alpha) - alpha1
   
   if (dalpha < 0.0_ReKi) then
      Get_f = 1 - 0.3_ReKi*exp(dalpha/S1)
   else
      Get_f = 0.04_ReKi + 0.66_ReKi*exp(-dalpha/S2)
   end if
   
end function Get_f

real(ReKi) function Get_Dp(ds, T_p, Dp_minus1, Cn_pot, Cn_pot_minus1)
   ! Called by : Get_Cn_prime
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: ds
   real(ReKi), intent(in   ) :: T_p
   real(ReKi), intent(in   ) :: Dp_minus1
   real(ReKi), intent(in   ) :: Cn_pot
   real(ReKi), intent(in   ) :: Cn_pot_minus1
   
   ! Implements Equation 1.32b
   Get_Dp = Dp_minus1 * exp(-ds/T_p) + (Cn_pot - Cn_pot_minus1) * exp(-0.5_ReKi*ds*T_p)

end function Get_Dp



!real(ReKi) function Get_Cn_prime( dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
!                                 q_minus1, q_minus2, T_alpha, T_q, T_p, X_1_minus1, X_2_minus1, Cn_alpha_q_nc, Kprime_alpha_minus1, Kprime_q_minus1,Dp_minus1, Cn_pot_minus1, b1, &
!                                 b2, A1, A2 )
!   ! Called by : Get_Cn, UnsteadyAero_UpdateStates
!   ! Calls  to : Get_Cn_pot, Get_Dp
!   ! Implements Equation 1.32a
!   real(DbKi), intent(in   ) :: dt
!   real(ReKi), intent(in   ) :: ds
!   real(ReKi), intent(in   ) :: U
!   real(ReKi), intent(in   ) :: M
!   real(ReKi), intent(in   ) :: beta_M_sqrd
!   real(ReKi), intent(in   ) :: c
!   real(ReKi), intent(in   ) :: C_nalpha_circ
!   real(ReKi), intent(in   ) :: alpha_cur
!   real(ReKi), intent(in   ) :: alpha_minus1
!   real(ReKi), intent(in   ) :: alpha_minus2
!   real(ReKi), intent(in   ) :: alpha0
!   real(ReKi), intent(in   ) :: q_cur
!   real(ReKi), intent(in   ) :: q_minus1
!   real(ReKi), intent(in   ) :: q_minus2
!   real(ReKi), intent(in   ) :: T_alpha
!   real(ReKi), intent(in   ) :: T_q
!   real(ReKi), intent(in   ) :: T_p
!   real(ReKi), intent(in   ) :: X_1_minus1
!   real(ReKi), intent(in   ) :: X_2_minus1
!   real(ReKi), intent(in   ) :: Cn_alpha_q_nc
!   real(ReKi), intent(in   ) :: Kprime_alpha_minus1
!   real(ReKi), intent(in   ) :: Kprime_q_minus1
!   real(ReKi), intent(in   ) :: Dp_minus1
!   real(ReKi), intent(in   ) :: Cn_pot_minus1
!   real(ReKi), intent(in   ) :: b1
!   real(ReKi), intent(in   ) :: b2
!   real(ReKi), intent(in   ) :: A1
!   real(ReKi), intent(in   ) :: A2
!   
!   real(ReKi)                :: Cn_pot
!   real(ReKi)                :: Dp
!   
!   
!   ! Implements equation 1.32a
!   Cn_pot = Get_Cn_pot( dt, ds,U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
!                        q_minus1, q_minus2, T_alpha, T_q, X_1_minus1, X_2_minus1, Cn_alpha_q_nc, Kprime_alpha_minus1, Kprime_q_minus1, b1, b2, A1, A2 ) 
!   Dp     = Get_Dp    ( ds, T_p, Dp_minus1, Cn_pot, Cn_pot_minus1 )
!   Get_Cn_prime = Cn_pot -  Dp
!
!end function Get_Cn_prime

real(ReKi) function Get_alpha_f( Cn_prime, C_nalpha, alpha0 )
   ! Called by : Get_f_prime
   ! Calls  to : NONE
   ! Implements Equation 1.31
   real(ReKi), intent(in   ) :: Cn_prime
   real(ReKi), intent(in   ) :: C_nalpha
   real(ReKi), intent(in   ) :: alpha0
   
   Get_alpha_f = Cn_prime / C_nalpha + alpha0
   
end function Get_alpha_f

real(ReKi) function Get_f_prime( dt, ds, U, M, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, &
                                q_cur, q_minus1, q_minus2, T_alpha, T_q, T_p, X_1_minus1, X_2_minus1, Cn_prime, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, b1, b2, A1, &
                                A2, C_nalpha, alpha1, S1, S2 )
   ! Called by : Get_f_primeprime
   ! Calls  to : Get_Cn_prime, Get_alpha_f, Get_f
   real(DbKi), intent(in   ) :: dt
   real(ReKi), intent(in   ) :: ds
   real(ReKi), intent(in   ) :: U
   real(ReKi), intent(in   ) :: M
   real(ReKi), intent(in   ) :: c
   real(ReKi), intent(in   ) :: C_nalpha_circ
   real(ReKi), intent(in   ) :: alpha_cur
   real(ReKi), intent(in   ) :: alpha_minus1
   real(ReKi), intent(in   ) :: alpha_minus2
   real(ReKi), intent(in   ) :: alpha0
   real(ReKi), intent(in   ) :: q_cur
   real(ReKi), intent(in   ) :: q_minus1
   real(ReKi), intent(in   ) :: q_minus2
   real(ReKi), intent(in   ) :: T_alpha
   real(ReKi), intent(in   ) :: T_q
   real(ReKi), intent(in   ) :: T_p
   real(ReKi), intent(in   ) :: X_1_minus1
   real(ReKi), intent(in   ) :: X_2_minus1
   real(ReKi), intent(in   ) :: Cn_prime
   real(ReKi), intent(in   ) :: Kprime_alpha_minus1
   real(ReKi), intent(in   ) :: Kprime_q_minus1
   real(ReKi), intent(in   ) :: Dp_minus1
   real(ReKi), intent(in   ) :: Cn_pot_minus1
   real(ReKi), intent(in   ) :: b1
   real(ReKi), intent(in   ) :: b2
   real(ReKi), intent(in   ) :: A1
   real(ReKi), intent(in   ) :: A2
   real(ReKi), intent(in   ) :: C_nalpha
   real(ReKi), intent(in   ) :: alpha1
   real(ReKi), intent(in   ) :: S1
   real(ReKi), intent(in   ) :: S2
      
   real(ReKi)                :: alpha_f
   
   
   alpha_f  = Get_alpha_f( Cn_prime, C_nalpha, alpha0 )

   ! Use equation 1.30 to obtain f_prime, substituting alpha_f for alpha in the call to Get_f()
   Get_f_prime = Get_f(alpha_f, alpha1, S1, S2)

end function Get_f_prime


real(ReKi) function Get_Df(ds, T_f, Df_minus1, fprime, fprime_minus1)
   ! Called by : Get_f_primeprime
   ! Calls  to : NONE
   real(ReKi), intent(in   ) :: ds
   real(ReKi), intent(in   ) :: T_f
   real(ReKi), intent(in   ) :: Df_minus1
   real(ReKi), intent(in   ) :: fprime
   real(ReKi), intent(in   ) :: fprime_minus1
   ! Implements Equation 1.33b

   real(ReKi)                :: factor
   factor = -ds/T_f
   Get_Df = Df_minus1 * exp(factor) + (fprime - fprime_minus1)*exp(0.5_ReKi*factor)
   
end function Get_Df

real(ReKi) function Get_f_primeprime( dt, ds, U, M, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
                        q_minus1, q_minus2, T_alpha, T_q, T_p, T_f, X_1_minus1, X_2_minus1, Cn_prime, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, Df_minus1, fprime_minus1, b1, b2, A1, A2, C_nalpha, &
                        alpha1, S1, S2 )
   ! Called by : Get_Cv, UnsteadyAero_UpdateStates
   ! Calls  to : Get_f_prime, Get_Df
   ! Implements Equation 1.33a

   real(DbKi), intent(in   ) :: dt
   real(ReKi), intent(in   ) :: ds
   real(ReKi), intent(in   ) :: U
   real(ReKi), intent(in   ) :: M
   real(ReKi), intent(in   ) :: c
   real(ReKi), intent(in   ) :: C_nalpha_circ
   real(ReKi), intent(in   ) :: alpha_cur
   real(ReKi), intent(in   ) :: alpha_minus1
   real(ReKi), intent(in   ) :: alpha_minus2
   real(ReKi), intent(in   ) :: alpha0
   real(ReKi), intent(in   ) :: q_cur
   real(ReKi), intent(in   ) :: q_minus1
   real(ReKi), intent(in   ) :: q_minus2
   real(ReKi), intent(in   ) :: T_alpha
   real(ReKi), intent(in   ) :: T_q
   real(ReKi), intent(in   ) :: T_p
   real(ReKi), intent(in   ) :: T_f
   real(ReKi), intent(in   ) :: X_1_minus1
   real(ReKi), intent(in   ) :: X_2_minus1
   real(ReKi), intent(in   ) :: Cn_prime
   real(ReKi), intent(in   ) :: Kprime_alpha_minus1
   real(ReKi), intent(in   ) :: Kprime_q_minus1
   real(ReKi), intent(in   ) :: Dp_minus1
   real(ReKi), intent(in   ) :: Cn_pot_minus1
   real(ReKi), intent(in   ) :: Df_minus1
   real(ReKi), intent(in   ) :: fprime_minus1
   real(ReKi), intent(in   ) :: b1
   real(ReKi), intent(in   ) :: b2
   real(ReKi), intent(in   ) :: A1
   real(ReKi), intent(in   ) :: A2
   real(ReKi), intent(in   ) :: C_nalpha
   real(ReKi), intent(in   ) :: alpha1
   real(ReKi), intent(in   ) :: S1
   real(ReKi), intent(in   ) :: S2
   
   real(ReKi)                :: fprime
   fprime = Get_f_prime(dt, ds, U, M, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
                        q_minus1, q_minus2, T_alpha, T_q, T_p, X_1_minus1, X_2_minus1, Cn_prime, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, b1, b2, A1, A2, C_nalpha, &
                        alpha1, S1, S2)
   Get_f_primeprime =  fprime - Get_Df(ds, T_f, Df_minus1, fprime, fprime_minus1)

end function Get_f_primeprime

real(ReKi) function Get_Cn_FS( Cn_alpha_q_nc, Cn_alpha_q_circ, fprimeprime, DSMod )


   ! Called by : 
   ! Calls  to : 

   
   real(ReKi), intent(in   ) :: Cn_alpha_q_nc
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ
   real(ReKi), intent(in   ) :: fprimeprime
   integer   , intent(in   ) :: DSMod
   
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

real(ReKi) function Get_Cc_FS( eta_e, alpha_e, alpha0, Cn_alpha_q_circ, fprimeprime, DSMod )


   ! Called by : 
   ! Calls  to : 

   
   real(ReKi), intent(in   ) :: eta_e
   real(ReKi), intent(in   ) :: alpha_e
   real(ReKi), intent(in   ) :: alpha0
   real(ReKi), intent(in   ) :: Cn_alpha_q_circ
   real(ReKi), intent(in   ) :: fprimeprime
   integer   , intent(in   ) :: DSMod
   
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
! Dynamic Stall Equations                                                     !
!==============================================================================

real(ReKi) function Get_C_V( Cn_alpha_q_circ, fprimeprime, DSMod )
   ! Called by : Get_Cn
   ! Calls  to : Get_f_primeprime

   real(ReKi), intent(in   ) :: Cn_alpha_q_circ
   real(ReKi), intent(in   ) :: fprimeprime
   integer   , intent(in   ) :: DSMod
   
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

!real(ReKi) function Get_Cn_v(Cv, Cv_minus1, Cn_v_minus1, ds, T_V, doAccel)
real(ReKi) function Get_Cn_v ( ds, T_V, Cn_v_minus1, C_V, C_V_minus1, tau_V, T_VL, Kalpha, alpha, alpha0 ) ! do I pass VRTX flag or tau_V?
   ! Called by : Get_Cn
   ! Calls  to : NONE
   ! Implements Equation 1.45 or 1.49
   real(ReKi), intent(in   ) :: ds
   real(ReKi), intent(in   ) :: T_V
   real(ReKi), intent(in   ) :: Cn_v_minus1
   real(ReKi), intent(in   ) :: C_V
   real(ReKi), intent(in   ) :: C_V_minus1
   real(ReKi), intent(in   ) :: tau_V
   real(ReKi), intent(in   ) :: T_VL
   real(ReKi), intent(in   ) :: Kalpha
   real(ReKi), intent(in   ) :: alpha
   real(ReKi), intent(in   ) :: alpha0

   
   real(ReKi)                :: factor
   
   factor = (alpha - alpha0) * Kalpha
   
   if (tau_V > T_VL .AND. factor >= 0.0_ReKi ) then  
      Get_Cn_v = Cn_v_minus1*exp(-ds/T_V)   ! Eqn 1.49
   else      
      Get_Cn_v = Get_ExpEqn( ds, T_V, Cn_v_minus1, C_V, C_V_minus1 )   ! Eqn 1.45
   end if
   
end function Get_Cn_v


   
!real(ReKi) function Get_Cn( dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
!                             q_minus1, q_minus2, T_alpha, T_q, T_p, T_f, T_V, X_1_minus1, X_2_minus1, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, Cn_pot_minus1, Df_minus1, fprime_minus1, Cv_minus1, Cn_v_minus1, b1, b2, A1, A2, C_nalpha, &
!                             alpha1, S1, S2, DSMod, doAccel )
!   ! Called by : NONE
!   ! Calls  to : Get_Cv, Get_Cn_Alpha_q_nc, Get_Cn_v
!
!   real(DbKi), intent(in   ) :: dt
!   real(ReKi), intent(in   ) :: ds
!   real(ReKi), intent(in   ) :: U
!   real(ReKi), intent(in   ) :: M
!   real(ReKi), intent(in   ) :: beta_M_sqrd
!   real(ReKi), intent(in   ) :: c
!   real(ReKi), intent(in   ) :: C_nalpha_circ
!   real(ReKi), intent(in   ) :: alpha_cur
!   real(ReKi), intent(in   ) :: alpha_minus1
!   real(ReKi), intent(in   ) :: alpha_minus2
!   real(ReKi), intent(in   ) :: alpha0
!   real(ReKi), intent(in   ) :: q_cur
!   real(ReKi), intent(in   ) :: q_minus1
!   real(ReKi), intent(in   ) :: q_minus2
!   real(ReKi), intent(in   ) :: T_alpha
!   real(ReKi), intent(in   ) :: T_q
!   real(ReKi), intent(in   ) :: T_p
!   real(ReKi), intent(in   ) :: T_f
!   real(ReKi), intent(in   ) :: T_V
!   real(ReKi), intent(in   ) :: X_1_minus1
!   real(ReKi), intent(in   ) :: X_2_minus1
!   real(ReKi), intent(in   ) :: Kprime_alpha_minus1
!   real(ReKi), intent(in   ) :: Kprime_q_minus1
!   real(ReKi), intent(in   ) :: Dp_minus1
!   real(ReKi), intent(in   ) :: Cn_pot_minus1
!   real(ReKi), intent(in   ) :: Df_minus1
!   real(ReKi), intent(in   ) :: fprime_minus1
!   real(ReKi), intent(in   ) :: Cv_minus1
!   real(ReKi), intent(in   ) :: Cn_v_minus1
!
!   real(ReKi), intent(in   ) :: b1
!   real(ReKi), intent(in   ) :: b2
!   real(ReKi), intent(in   ) :: A1
!   real(ReKi), intent(in   ) :: A2
!   real(ReKi), intent(in   ) :: C_nalpha
!   real(ReKi), intent(in   ) :: alpha1
!   real(ReKi), intent(in   ) :: S1
!   real(ReKi), intent(in   ) :: S2
!   integer   , intent(in   ) :: DSMod
!   logical   , intent(in   ) :: doAccel
!   
!   
!   real(ReKi)                :: Cv
!   real(ReKi)                :: Cn_alpha_q_nc
!   real(ReKi)                :: Cn_prime
!   
!   ! Implements equation 1.50
!   Cn_alpha_q_nc = Get_Cn_Alpha_q_nc( T_alpha, alpha_cur, alpha_minus1, alpha_minus2, Kprime_alpha_minus1, Kprime_q_minus1, T_q, q_cur, q_minus1, q_minus2, dt ) 
!   Cn_prime      = Get_Cn_prime( dt, ds, U, M, beta_M_sqrd, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
!                                 q_minus1, q_minus2, T_alpha, T_q, T_p, X_1_minus1, X_2_minus1, Cn_alpha_q_nc, Kprime_alpha_minus1, Kprime_q_minus1,Dp_minus1, Cn_pot_minus1, b1, &
!                                 b2, A1, A2 )
!   Cv            = Get_Cv( dt, ds, U, M, c, C_nalpha_circ, alpha_cur, alpha_minus1, alpha_minus2, alpha0, q_cur, &
!                             q_minus1, q_minus2, T_alpha, T_q, T_p, T_f, X_1_minus1, X_2_minus1, Cn_prime, Cn_alpha_q_nc, Kprime_alpha_minus1, Kprime_q_minus1, Dp_minus1, &
!                             Cn_pot_minus1, Df_minus1, fprime_minus1, b1, b2, A1, A2, C_nalpha, &
!                             alpha1, S1, S2, DSMod  ) 
!   Get_Cn        = C_n_FS + Get_Cn_v( Cv, Cv_minus1, Cn_v_minus1, ds, T_V, doAccel )
!   Get_Cn        = Cv + Cn_alpha_q_nc +  Get_Cn_v( Cv, Cv_minus1, Cn_v_minus1, ds, T_V, doAccel )
!   
!   
!end function Get_Cn


subroutine UnsteadyAero_SetParameters( dt, InitInp, p, errStat, errMsg )
   real(DbKi),                             intent(inout)  :: dt
   type(UnsteadyAero_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UnsteadyAero_ParameterType),       intent(inout)  :: p           ! Parameters
   integer(IntKi),                         intent(  out)  :: errStat     ! Error status of the operation
   character(*),                           intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

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
   p%AFI_Params = InitInp%AFI_Params
   allocate(p%AFIndx(InitInp%nNodesPerBlade,InitInp%numBlades), stat = errStat)
   if (errStat > ErrID_None) then
      ! Set errmessage and return
      deallocate(p%c)
      return
   end if
   p%AFIndx     = InitInp%AFIndx
   
end subroutine UnsteadyAero_SetParameters



subroutine UnsteadyAero_InitStates( p, xd, OtherState, errStat, errMsg )
   type(UnsteadyAero_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UnsteadyAero_DiscreteStateType),   intent(inout)  :: xd          ! Initial discrete states
   type(UnsteadyAero_OtherStateType),      intent(inout)  :: OtherState  ! Initial other/optimization states
   integer(IntKi),                         intent(  out)  :: errStat     ! Error status of the operation
   character(*),                           intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

   integer                   :: i, j
   real(ReKi)                :: C_nalpha_circ
   real(ReKi)                :: alpha0       
   real(ReKi)                :: alpha_e
   real(ReKi)                :: T_f0         
   real(ReKi)                :: T_V0  
   real(ReKi)                :: T_VL
   real(ReKi)                :: T_p                   
   real(ReKi)                :: b1           
   real(ReKi)                :: b2           
   real(ReKi)                :: A1           
   real(ReKi)                :: A2           
   real(ReKi)                :: C_nalpha              
   real(ReKi)                :: alpha1 
   real(ReKi)                :: eta_e
   real(ReKi)                :: S1           
   real(ReKi)                :: S2  
   real(ReKi)                :: Cn1
   real(ReKi)                :: Cn2
   real(ReKi)                :: St_sh
   
   ! allocate all the state arrays
   allocate(xd%alpha_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%alpha_minus2(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%q_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%q_minus2(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%X1_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%X2_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Kprime_alpha_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Kprime_q_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Dp_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%Cn_pot_minus1(p%nNodesPerBlade,p%numBlades), stat = errStat)
   allocate(xd%T_f(p%nNodesPerBlade,p%numBlades), stat = errStat)
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
   
   OtherState%sigma1 = 1.0_ReKi
   OtherState%sigma3 = 1.0_ReKi
   OtherState%TESF   = .FALSE.  
   OtherState%LESF   = .FALSE.   
   OtherState%VRTX    = .FALSE.  
   OtherState%FirstPass = .true.
   
   xd%alpha_minus1         = 0.0_ReKi
   xd%alpha_minus2       = 0.0_ReKi
   xd%q_minus1             = 0.0_ReKi
   xd%q_minus2           = 0.0_ReKi
   xd%X1_minus1            = 0.0_ReKi
   xd%X2_minus1            = 0.0_ReKi
   xd%Kprime_alpha_minus1  = 0.0_ReKi
   xd%Kprime_q_minus1      = 0.0_ReKi
   xd%Dp_minus1            = 0.0_ReKi
   xd%Cn_pot_minus1        = 0.0_ReKi
   


   

   do j = 1,p%numBlades
      do i = 1,p%nNodesPerBlade
            ! Lookup values using Airfoil Info module
            ! Use reference values for M = 0.75 (million), Re = 1.0, and alpha = 0.0
         call AFI_GetAirfoilParams( p%AFI_Params%AFInfo(p%AFIndx(i,j)), 0.75, 1.0, 0.0, alpha0, alpha1, eta_e, C_nalpha, &
                                 C_nalpha_circ, T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, A1, A2, S1, S2, Cn1, Cn2, errMsg, errStat )           
         xd%T_f(i,j)          = T_f0
      end do
   end do       
   xd%fprimeprime_minus1   = 0.0_ReKi
   xd%Df_minus1            = 0.0_ReKi
   xd%fprime_minus1        = 0.0_ReKi
   xd%tau_V              = 0.0_ReKi 
   xd%Cn_v_minus1          = 0.0_ReKi
   xd%C_V_minus1           = 0.0_ReKi

end subroutine UnsteadyAero_InitStates


!----------------------------------------------------------------------------------------------------------------------------------
subroutine UnsteadyAero_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   type(UnsteadyAero_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in BEMT_SetParameters()
   type(UnsteadyAero_InputType),           intent(in   )  :: u           ! An initial guess for the input; input mesh must be defined
   type(UnsteadyAero_ParameterType),       intent(  out)  :: p           ! Parameters
   type(UnsteadyAero_ContinuousStateType), intent(  out)  :: x           ! Initial continuous states
   type(UnsteadyAero_DiscreteStateType),   intent(  out)  :: xd          ! Initial discrete states
   type(UnsteadyAero_ConstraintStateType), intent(  out)  :: z           ! Initial guess of the constraint states
   type(UnsteadyAero_OtherStateType),      intent(  out)  :: OtherState  ! Initial other/optimization states
   type(UnsteadyAero_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                         !   only the output mesh is initialized)
   real(DbKi),                             intent(inout)  :: interval    ! Coupling interval in seconds: the rate that
                                                                         !   (1) BEMT_UpdateStates() is called in loose coupling &
                                                                         !   (2) BEMT_UpdateDiscState() is called in tight coupling.
                                                                         !   Input is the suggested time from the glue code;
                                                                         !   Output is the actual coupling interval that will be used
                                                                         !   by the glue code.
   type(UnsteadyAero_InitOutputType),      intent(  out)  :: InitOut     ! Output for initialization routine
   integer(IntKi),                         intent(  out)  :: errStat     ! Error status of the operation
   character(*),                           intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                         :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! temporary Error status of the operation
   integer(IntKi)                                 :: i,j, iNode
   character(64)                                :: chanPrefix
      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""


      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
   !call DispNVD( UnsteadyAero_Ver )
   
      ! Allocate and set parameter data structure using initialization data
   call UnsteadyAero_SetParameters( interval, InitInp, p, errStat, errMsg )
   if (errStat >= AbortErrLev) return

   
      ! initialize the discrete states
   call UnsteadyAero_InitStates( p, xd, OtherState, errStat, errMsg )     ! initialize the continuous states
   if (errStat >= AbortErrLev) return
      
   
      ! Allocate and set the InitOut data
   p%NumOuts = p%numBlades*p%nNodesPerBlade*6
   p%Delim   = ''
   p%OutFmt  = 'ES11.4e2'
   p%OutSFmt = 'A11'
   allocate(InitOut%WriteOutputHdr(p%NumOuts))
   
   allocate(InitOut%WriteOutputUnt(p%NumOuts))
   
   iNode = 0
   do i = 1,p%numBlades
      do j = 1,p%nNodesPerBlade
         iNode = iNode + 1
         chanPrefix = "B"//trim(num2lstr(i))//"N"//trim(num2lstr(j))  
         InitOut%WriteOutputHdr(iNode) =trim(chanPrefix)//'AOA'
         InitOut%WriteOutputHdr(iNode+1) =trim(chanPrefix)//'Cn' 
         InitOut%WriteOutputHdr(iNode+2) =trim(chanPrefix)//'Cc' 
         InitOut%WriteOutputHdr(iNode+3) =trim(chanPrefix)//'Cm' 
         InitOut%WriteOutputHdr(iNode+4) =trim(chanPrefix)//'Cl' 
         InitOut%WriteOutputHdr(iNode+5) =trim(chanPrefix)//'Cd' 
         InitOut%WriteOutputUnt(iNode)   ='(deg)  '
         InitOut%WriteOutputUnt(iNode+1) ='(-)   '
         InitOut%WriteOutputUnt(iNode+2) ='(-)   '
         InitOut%WriteOutputUnt(iNode+3) ='(-)   '
         InitOut%WriteOutputUnt(iNode+4) ='(-)   '
         InitOut%WriteOutputUnt(iNode+5) ='(-)   '
      end do
   end do
   
      ! Allocate outputs
   allocate(y%Cn(p%nNodesPerBlade, p%numBlades))
   allocate(y%Cc(p%nNodesPerBlade, p%numBlades))
   allocate(y%Cm(p%nNodesPerBlade, p%numBlades))
   allocate(y%Cl(p%nNodesPerBlade, p%numBlades))
   allocate(y%Cd(p%nNodesPerBlade, p%numBlades))
   allocate(y%WriteOutput(p%NumOuts))
   
end subroutine UnsteadyAero_Init
   
subroutine ComputeKelvinChain( i, j, u, p, xd, OtherState, Cn_prime, Cn1, Cn2, T_VL, St_sh, Kalpha, alpha_e, alpha0, dalpha0, eta_e, Kq, q_cur, X1, X2, &
                                 Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, fprime, Cn_v, C_V, Cn_FS, errStat, errMsg )
                              
! 
!
! 
!..................................................................................................................................

      
   integer(IntKi),                         intent(in   ) :: i,j          ! Blade node index
   type(UnsteadyAero_InputType),           intent(in   ) :: u          ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   type(UnsteadyAero_ParameterType),       intent(in   ) :: p          ! Parameters   TODO: this should only be in, but is needed because of a copy of AFInfo in a called subroutine. GJH 1/5/2015
   type(UnsteadyAero_DiscreteStateType),   intent(in   ) :: xd         ! Input: Discrete states at t;
   type(UnsteadyAero_OtherStateType),      intent(in   ) :: OtherState ! Other/optimization states
   real(ReKi),                             intent(  out) :: Cn_prime
   real(ReKi),                             intent(  out) :: Cn1
   real(ReKi),                             intent(  out) :: Cn2
   real(ReKi),                             intent(  out) :: T_VL
   real(ReKi),                             intent(  out) :: St_sh
   real(ReKi),                             intent(  out) :: Kalpha
   real(ReKi),                             intent(  out) :: alpha_e
   real(ReKi),                             intent(  out) :: alpha0
   real(ReKi),                             intent(  out) :: dalpha0
   real(ReKi),                             intent(  out) :: eta_e
   real(ReKi),                             intent(  out) :: Kq
   real(ReKi),                             intent(  out) :: q_cur
   real(ReKi),                             intent(  out) :: X1
   real(ReKi),                             intent(  out) :: X2
   real(ReKi),                             intent(  out) :: Kprime_alpha
   real(ReKi),                             intent(  out) :: Kprime_q
   real(ReKi),                             intent(  out) :: Dp
   real(ReKi),                             intent(  out) :: Cn_pot
   real(ReKi),                             intent(  out) :: Cc_pot
   real(ReKi),                             intent(  out) :: fprimeprime
   real(ReKi),                             intent(  out) :: Df
   real(ReKi),                             intent(  out) :: fprime
   real(ReKi),                             intent(  out) :: Cn_v
   real(ReKi),                             intent(  out) :: C_V
   real(ReKi),                             intent(  out) :: Cn_FS
      
   integer(IntKi),                         intent(  out) :: errStat    ! Error status of the operation
   character(*),                           intent(  out) :: errMsg     ! Error message if ErrStat /= ErrID_None


   real(ReKi)                :: ds
   real(ReKi)                :: M
   real(ReKi)                :: beta
   real(ReKi)                :: betaSqrd
   real(ReKi)                :: Cn_alpha_q_nc
   real(ReKi)                :: k_alpha
   real(ReKi)                :: T_alpha
   real(ReKi)                :: T_q
   real(ReKi)                :: T_f
   real(ReKi)                :: T_V
   real(ReKi)                :: T_I 
   real(ReKi)                :: C_nalpha_circ
   real(ReKi)                :: T_f0         
   real(ReKi)                :: T_V0         
   real(ReKi)                :: T_p                   
   real(ReKi)                :: b1           
   real(ReKi)                :: b2           
   real(ReKi)                :: A1           
   real(ReKi)                :: A2           
   real(ReKi)                :: C_nalpha              
   real(ReKi)                :: alpha1          
   real(ReKi)                :: S1           
   real(ReKi)                :: S2   
   real(ReKi)                :: Cn_alpha_nc
   real(ReKi)                :: Cn_alpha_q_circ
   real(ReKi)                :: Cn_q_nc
   real(ReKi)                :: alpha_f
   real(ReKi)                :: k_q
   real(ReKi)                :: Kalpha_minus1
   real(ReKi)                :: Kq_minus1
   real(ReKi)                :: alpha_minus1
   real(ReKi)                :: alpha_minus2
   real(ReKi)                :: q_minus1
   real(ReKi)                :: q_minus2
   real(ReKi)                :: fprime_minus1
   real(ReKi)                :: Cn_pot_minus1
   
   ! Called by : DRIVER
   ! Calls  to : Get_Beta_M_Sqrd, AFI_GetAirfoilParams, Get_ds, Get_Pitchrate, Get_T_Alpha, Get_T_q, Get_Cn_prime, Get_f_primeprime, Get_k_Alpha



   M        = u%U(i,j) / p%a_s
   betaSqrd = Get_Beta_M_Sqrd(M)
   beta     = Get_Beta_M(M)
   T_I      = p%c(i,j) / p%a_s                                                  ! Eqn 1.11c
   ds       = Get_ds       ( u%U(i,j), p%c(i,j), p%dt )
   
   ! Lookup values using Airfoil Info module
   call AFI_GetAirfoilParams( p%AFI_Params%AFInfo(p%AFIndx(i,j)), M, u%Re(i,j), u%alpha(i,j), alpha0, alpha1, eta_e, C_nalpha, C_nalpha_circ, T_f0, T_V0, T_p, T_VL, St_sh, b1, b2, A1, A2, S1, S2, Cn1, Cn2, errMsg, errStat )           
   
   if (OtherState%FirstPass) then
      alpha_minus1 = u%alpha(i,j) 
      alpha_minus2 = u%alpha(i,j) 
   else
      alpha_minus1 = xd%alpha_minus1(i,j)
      alpha_minus2 = xd%alpha_minus2(i,j)
   end if
   
   dalpha0  = u%alpha(i,j) - alpha0
   q_cur    = Get_Pitchrate( p%dt, u%alpha(i,j), alpha_minus1, u%U(i,j), p%c(i,j)  )  ! Eqn 1.8
   
   if (OtherState%FirstPass) then
      q_minus1 = q_cur   
      q_minus2 = q_cur 
   else
      q_minus1 = xd%q_minus1(i,j)   
      q_minus2 = xd%q_minus2(i,j) 
   end if
   
   k_alpha  = Get_k_       ( M, beta, C_nalpha/2.0_ReKi )                 ! Eqn 1.11a
   k_q      = Get_k_       ( M, beta, C_nalpha          )                 ! Eqn 1.11b
   T_alpha  = T_I * k_alpha                                                   ! Eqn 1.10a
   T_q      = T_I * k_q                                                       ! Eqn 1.10b
      
   ! These quantities are needed for the update state calculations, but are then updated themselves based on the logic which follows
   T_f           = T_f0 / OtherState%sigma1(i,j)
   T_V           = T_V0 / OtherState%sigma3(i,j)
   
   ! Compute Kalpha using Eqn 1.18b and Kalpha_minus1 using Eqn 1.18b but with time-shifted alphas
   Kalpha          = Get_Kupper( p%dt, u%alpha      (i,j), alpha_minus1   )
   Kalpha_minus1   = Get_Kupper( p%dt, alpha_minus1, alpha_minus2 ) 
   
   ! Compute Kq using Eqn 1.19b and Kq_minus1 using Eqn 1.19b but with time-shifted alphas
   Kq              = Get_Kupper( p%dt, q_cur           , q_minus1       )
   Kq_minus1       = Get_Kupper( p%dt, q_minus1    , q_minus2     )
   
   ! Compute Kprime_alpha using Eqn 1.18c
   !if (OtherState%FirstPass) then
   !   Kprime_alpha  = 0.0_ReKi
   !else     
      Kprime_alpha = Get_ExpEqn( real(p%dt), T_alpha, xd%Kprime_alpha_minus1(i,j), Kalpha, Kalpha_minus1 )
   !end if
   ! Compute Kprime_q using Eqn 1.19c
   !if (OtherState%FirstPass) then
   !   Kprime_q  = 0.0_ReKi
   !else    
      Kprime_q     = Get_ExpEqn( real(p%dt), T_q    , xd%Kprime_q_minus1(i,j)    ,  Kq   , Kq_minus1     )
   !end if
   
   ! Compute Cn_alpha_nc using Eqn 1.18a
   ! Depends on T_alpha, M, alpha (1.18b), Kprime_alpha (1.18c), xd%alpha_minus1, xd%Kprime_alpha_minus1
   Cn_alpha_nc     = Get_Cn_nc( M, 4.0_ReKi*T_alpha, Kalpha, Kprime_alpha )
   
   ! Compute Cn_alpha_nc using Eqn 1.19a
   Cn_q_nc         = Get_Cn_nc( M,          T_q    , Kq    , Kprime_q     )
   
   ! Compute Cn_alpha_q_nc using Eqn 1.17
   ! Depends on Cn_alpha_nc (1.18a), Cn_q_nc (1.19a)
   Cn_alpha_q_nc   = Cn_alpha_nc + Cn_q_nc  
   
   ! Compute X1 using Eqn 1.15a
   !if (OtherState%FirstPass) then
   !   X1  = 0.0_ReKi
   !else    
      X1           = Get_X(ds, betaSqrd, b1, A1, u%alpha(i,j), alpha_minus1, xd%X1_minus1(i,j) )
   !end if
   ! Compute X2 using Eqn 1.15b
   !if (OtherState%FirstPass) then
   !   X2  = 0.0_ReKi
   !else
      X2           = Get_X(ds, betaSqrd, b2, A2, u%alpha(i,j), alpha_minus1, xd%X2_minus1(i,j) )
   !end if
   
   ! Compute alpha_e using Eqn 1.14
   alpha_e         = Get_alpha_e( u%alpha(i,j), alpha0, X1, X2 )
   
   ! Compute Cn_alpha_q_circ using Eqn 1.13
   Cn_alpha_q_circ = C_nalpha_circ * alpha_e
   
   
   ! Compute Cn_pot using eqn 1.20
   ! Depends on Cn_alpha_q_circ (1.13) , Cn_alpha_q_nc (1.17)
   Cn_pot          = Cn_alpha_q_circ + Cn_alpha_q_nc
   
   ! Compute Cc_pot using eqn 1.28
   Cc_pot          = Cn_pot*tan(alpha_e+alpha0)
   
   
   if (OtherState%FirstPass) then
      Cn_pot_minus1 = Cn_pot
   else
      Cn_pot_minus1 = xd%Cn_pot_minus1(i,j)
   end if
   
   ! Compute Dp using Eqn 1.32b
   ! Depends on xd%Dp_minus1, Tp, ds, xd%Cn_pot_minus1, Cn_pot (1.20)
   !if (OtherState%FirstPass) then
   !   Dp  = 0.0_ReKi
   !else
      Dp            = Get_ExpEqn( ds, T_p, xd%Dp_minus1(i,j), Cn_pot, Cn_pot_minus1 )
   !end if
   
   ! Compute Cn_prime using Eqn 1.32a
   ! Depends on Cn_pot (1.20) and Dp (1.32b)
   Cn_prime      = Cn_Pot - Dp
   
   ! Compute alpha_f using Eqn 1.31
   alpha_f       = Cn_prime / C_nalpha + alpha0
   
   ! Compute fprime using Eqn 1.30 and Eqn 1.31
   
   fprime        = Get_f( alpha_f, alpha1, S1, S2)
   if (OtherState%FirstPass) then
      fprime_minus1 = fprime
   else
      fprime_minus1 = xd%fprime_minus1(i,j)
      end if
   ! Compute Df using Eqn 1.33b
   !if (OtherState%FirstPass) then
   !   Df  = 0.0_ReKi
   !else
      Df            = Get_ExpEqn( ds, T_f, xd%Df_minus1(i,j), fprime, fprime_minus1 )
   !end if
   
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
   Cn_v          = Get_Cn_v ( ds, T_V, xd%Cn_v_minus1(i,j), C_V, xd%C_V_minus1(i,j), xd%tau_V(i,j), T_VL, Kalpha, u%alpha(i,j), alpha0 ) ! do I pass VRTX flag or tau_V?
   
   ! Finally, compute Cn using Eqn 1.50
   !Cn            = Cn_FS + Cn_v


end subroutine ComputeKelvinChain
   
!----------------------------------------------------------------------------------------------------------------------------------
subroutine UnsteadyAero_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )   
! Tight coupling routine for updating discrete states
!..................................................................................................................................
   
   real(DbKi),                             intent(in   )  :: Time        ! Current simulation time in seconds 
   integer(IntKi),                         intent(in   )  :: n           ! Current step of the simulation: t = n*Interval
   type(UnsteadyAero_InputType),           intent(in   )  :: u           ! Inputs at Time                       
   type(UnsteadyAero_ParameterType),       intent(in   )  :: p           ! Parameters                                 
   type(UnsteadyAero_ContinuousStateType), intent(in   )  :: x           ! Continuous states at Time
   type(UnsteadyAero_DiscreteStateType),   intent(inout)  :: xd          ! Input: Discrete states at Time; 
                                                                           ! Output: Discrete states at Time + Interval
   type(UnsteadyAero_ConstraintStateType), intent(in   )  :: z           ! Constraint states at Time
   type(UnsteadyAero_OtherStateType),      intent(inout)  :: OtherState  ! Other/optimization states           
   integer(IntKi),                         intent(  out)  :: errStat     ! Error status of the operation
   character(*),                           intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

         ! Local Variables
   integer                                                :: i, j
   real(ReKi)                                             :: Cn_prime
   real(ReKi)                                             :: Cn1
   real(ReKi)                                             :: Cn2
   real(ReKi)                                             :: T_VL   
   real(ReKi)                                             :: Kalpha
   real(ReKi)                                             :: Kq
   real(ReKi)                                             :: q_cur
   real(ReKi)                                             :: X1
   real(ReKi)                                             :: X2
   real(ReKi)                                             :: Kprime_alpha
   real(ReKi)                                             :: Kprime_q
   real(ReKi)                                             :: Dp
   real(ReKi)                                             :: Cn_pot
   real(ReKi)                                             :: Cc_pot
   real(ReKi)                                             :: fprimeprime
   real(ReKi)                                             :: Df
   real(ReKi)                                             :: fprime
   real(ReKi)                                             :: Cn_v
   real(ReKi)                                             :: C_V 
   real(ReKi)                                             :: Cn_FS
   real(ReKi)                                             :: alpha_e
   real(ReKi)                                             :: alpha0
   real(ReKi)                                             :: dalpha0
   real(ReKi)                                             :: eta_e
   real(ReKi)                                             :: Kafactor
   real(ReKi)                                             :: St_sh
   real(ReKi)                                             :: T_sh
            
      ! Initialize ErrStat
         
   errStat = ErrID_None         
   errMsg  = ""               

      ! We need to loop over all blade stations and update the states for each station
do j = 1,p%numBlades
   do i = 1,p%nNodesPerBlade
      call ComputeKelvinChain(i, j,  u, p, xd, OtherState, Cn_prime, Cn1, Cn2, T_VL, St_sh, Kalpha, alpha_e, alpha0, dalpha0, eta_e, Kq, q_cur, X1, X2, &
                                 Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, fprime, Cn_v, C_V, Cn_FS, errStat, errMsg )
      
      
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
      
      if ( xd%tau_V(i,j) >= (1.0_ReKi + T_sh/ T_VL) ) then
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
   
      
      
      
      
      !
      !
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! Process an existing vortex
      !xd%tau_V(i)        = xd%tau_V(i) + 2.0_ReKi*p%dt*u%U(i) / p%c(i)  !!
      !      
      !if ( OtherState%LESF(i) ) then
      !   ! We are in a condition where a new vortex can be produced or we are processing a current vortex
      !   if ( .not. (OtherState%VRTX(i) ) ) then
      !      ! Produce a new vortex
      !      OtherState%VRTX(i) = .true.
      !      xd%tau_V(i)        = 0.0_ReKi
      !   else
      !      ! Process an existing vortex
      !     !! xd%tau_V(i)        = xd%tau_V(i) + 2.0_ReKi*p%dt*u%U(i) / p%c(i)
      !      ! Do we need to generate a new vortex?
      !      if ( xd%tau_V(i) >= (1.0_ReKi + p%T_sh/ T_VL) ) then
      !         xd%tau_V(i)     = 0.0_ReKi
      !      elseif ( xd%tau_V(i) > 2.0_ReKi*T_VL ) then
      !         OtherState%VRTX(i) = .false.
      !      end if
      !   
      !   end if   
      !elseif ( OtherState%VRTX(i) ) then
      !   ! We are in a condition where we are only processing an existing vortex
      !   !!xd%tau_V(i)        = xd%tau_V(i) + 2.0_ReKi*p%dt*u%U(i) / p%c(i)
      !   ! Did the current vortex leave the near wake?
      !   if  ( xd%tau_V(i) > 2.0_ReKi*T_VL ) then
      !     !! xd%tau_V(i)        = 0.0_ReKi
      !      OtherState%VRTX(i) = .false.
      !   end if
      !
      !end if
      !
      !!!!!!!!!!!!!!!!
      !!if ( OtherState%VRTX(i) ) then
      !!   ! We are in a condition where we are only processing an existing vortex
      !!   xd%tau_V(i)        = xd%tau_V(i) + 2.0_ReKi*p%dt*u%U(i) / p%c(i)
      !!   ! Did the current vortex leave the near wake?
      !!   if  ( xd%tau_V(i) > 2.0_ReKi*T_VL ) then
      !!      OtherState%VRTX(i) = .false.
      !!      xd%tau_V(i)        = 0.0_ReKi
      !!   end if
      !!   
      !!   if xd%tau_V(i) >= (1.0_ReKi + p%T_sh/ T_VL) ) then
      !!      xd%tau_V(i)        = 0.0_ReKi
      !!   end if
      !!end if
      !!
      !!   if (  OtherState%LESF(i) ) then
      !!      OtherState%VRTX(i) = .true.
      !!   else
      !!      OtherState%VRTX(i) = .false.
      !!   end if
      !!   
      !!   end if
      !!   
      !!   ! Do we need to generate a new vortex?
      !!   if ( xd%tau_V(i) >= (1.0_ReKi + p%T_sh/ T_VL) ) then
      !!      xd%tau_V(i)     = 0.0_ReKi
      !!   end if   
      !!         
      !!end if
      !
      !
      !
      !!!!!!!!!!!!!!!
      !
      !if (OtherState%TESF(i)) then
      !   if (Kalpha < 0.0_ReKi) then 
      !      OtherState%sigma1(i) = 2.0_ReKi
      !   else if (.not. OtherState%LESF(i)) then
      !      OtherState%sigma1(i) = 1.0_ReKi
      !   else if (fprimeprime <= 0.7_ReKi) then
      !      OtherState%sigma1(i) = 2.0_ReKi
      !   else
      !      OtherState%sigma1(i) = 1.75_ReKi
      !   end if
      !else
      !   if (.not. OtherState%LESF(i)) then
      !      OtherState%sigma1(i) = 0.5_ReKi
      !   else if (OtherState%VRTX(i) .and. ( (xd%tau_V(i) <= T_VL) .and. (xd%tau_V(i) >= 0.0_ReKi) ) ) then
      !      OtherState%sigma1(i) = 0.25_ReKi
      !   else if (Kalpha > 0.0_ReKi) then
      !      OtherState%sigma1(i) = 0.75_ReKi
      !   end if
      !end if   
      !
      !if ( (xd%tau_V(i) <= 2.0_ReKi*T_VL) .and. (xd%tau_V(i) >= T_VL) ) then
      !   OtherState%sigma3(i) = 3.0_ReKi
      !else if (.not. OtherState%TESF(i)) then
      !   OtherState%sigma3(i) = 4.0_ReKi
      !else if ( (xd%tau_V(i) <= T_VL) .and. (xd%tau_V(i) >= 0.0_ReKi) ) then
      !   if (Kalpha < 0.0_ReKi) then
      !      OtherState%sigma3(i) = 2.0_ReKi
      !   else
      !      OtherState%sigma3(i) = 1.0_ReKi
      !   end if
      !else if (Kalpha < 0.0_ReKi) then
      !   if ( OtherState%TESF(i) ) then !!
      !      
      !!   OtherState%sigma3(i) = 4.0_ReKi  !!!
      !!else if ( .not. OtherState%TESF(i) ) !!!.and. Kalpha < 0.0 ) then !!!
      !!    OtherState%sigma3(i) = 1.0_ReKi!!!!!
      !!end if  !!!
      !
      !   if (OtherState%TESF(i)) then
      !      OtherState%sigma3(i) = 4.0_ReKi
      !   else
      !      OtherState%sigma3(i) = 1.0_ReKi
      !   end if
      !
      !end if
      !
     
      !---------------------------------------------------------
      ! Update the Discrete States, xd
      !---------------------------------------------------------
      if (OtherState%FirstPass) then
         xd%alpha_minus2(i,j)      = u%alpha(i,j)
         xd%q_minus2(i,j)          = q_cur
      else
         xd%alpha_minus2(i,j)      = xd%alpha_minus1(i,j)
         xd%q_minus2(i,j)          = xd%q_minus1(i,j)
      end if
           
      xd%alpha_minus1(i,j)        = u%alpha(i,j)   
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
         xd%tau_V(i,j)          = xd%tau_V(i,j) + 2.0_ReKi*p%dt*u%U(i,j) / p%c(i,j)  
      end if
   
   end do  
end do

end subroutine UnsteadyAero_UpdateDiscState


!----------------------------------------------------------------------------------------------------------------------------------
subroutine UnsteadyAero_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, errStat, errMsg )
! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
! Continuous, constraint, and discrete states are updated to values at t + Interval.
!..................................................................................................................................

   real(DbKi),                              intent(in   ) :: t               ! Current simulation time in seconds
   integer(IntKi),                          intent(in   ) :: n               ! Current step of the simulation: t = n*Interval
   type(UnsteadyAero_InputType),            intent(inout) :: Inputs(:)       ! Inputs at InputTimes
   real(DbKi),                              intent(in   ) :: InputTimes(:)   ! Times in seconds associated with Inputs
   type(UnsteadyAero_ParameterType),        intent(in   ) :: p               ! Parameters
   type(UnsteadyAero_ContinuousStateType),  intent(inout) :: x               ! Input: Continuous states at t;
                                                                             !   Output: Continuous states at t + Interval
   type(UnsteadyAero_DiscreteStateType),    intent(inout) :: xd              ! Input: Discrete states at t;
                                                                             !   Output: Discrete states at t + Interval
   type(UnsteadyAero_ConstraintStateType),  intent(inout) :: z               ! Input: Constraint states at t;
                                                                             !   Output: Constraint states at t + Interval
   type(UnsteadyAero_OtherStateType),       intent(inout) :: OtherState      ! Other/optimization states
   integer(IntKi),                          intent(  out) :: errStat         ! Error status of the operation
   character(*),                            intent(  out) :: errMsg          ! Error message if ErrStat /= ErrID_None

      ! Local variables  
   type(UnsteadyAero_DiscreteStateType)                   :: xd_t            ! Discrete states at t (copy)
   type(UnsteadyAero_InputType)                           :: u               ! Instantaneous inputs
   integer(IntKi)                                         :: errStat2        ! Error status of the operation (secondary error)
   character(len(ErrMsg))                                 :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   
      ! Initialize variables

   errStat   = ErrID_None           ! no error has occurred
   errMsg    = ""
            
   allocate(u%U(p%nNodesPerBlade, p%numBlades))   
   allocate(u%Re(p%nNodesPerBlade, p%numBlades))   
   allocate(u%alpha(p%nNodesPerBlade, p%numBlades))  
   
   ! This subroutine contains an example of how the states could be updated. Developers will
   ! want to adjust the logic as necessary for their own situations.
 
   
      ! Get the inputs at time t, based on the array of values sent by the glue code:
      
   call UnsteadyAero_Input_ExtrapInterp( Inputs, InputTimes, u, t, ErrStat, ErrMsg )  
   if ( ErrStat >= AbortErrLev ) return
       
      
      ! Update discrete states:
      !   Note that xd [discrete state] is changed in UnsteadyAero_UpdateDiscState() so xd will now contain values at t+Interval
      !   We'll first make a copy that contains xd at time t, which will be used in computing the constraint states
   !TODO:  Check if this is really necessary GJH 2/18/2015
   call UnsteadyAero_CopyDiscState( xd, xd_t, MESH_NEWCOPY, ErrStat, ErrMsg )
   if ( ErrStat >= AbortErrLev ) return

   call UnsteadyAero_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      
      
      ! Destroy local variables before returning    
   call UnsteadyAero_DestroyDiscState(   xd_t,       ErrStat2, ErrMsg2) 
      
end subroutine UnsteadyAero_UpdateStates

subroutine UnsteadyAero_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
   real(DbKi),                             intent(in   )  :: Time        ! Current simulation time in seconds
   type(UnsteadyAero_InputType),           intent(inout)  :: u           ! Inputs at Time
   type(UnsteadyAero_ParameterType),       intent(in   )  :: p           ! Parameters
   type(UnsteadyAero_ContinuousStateType), intent(in   )  :: x           ! Continuous states at Time
   type(UnsteadyAero_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at Time
   type(UnsteadyAero_ConstraintStateType), intent(in   )  :: z           ! Constraint states at Time
   type(UnsteadyAero_OtherStateType),      intent(inout)  :: OtherState  ! Other/optimization states
   type(UnsteadyAero_OutputType),          intent(inout)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                     !   nectivity information does not have to be recalculated)
   integer(IntKi),                         intent(  out)  :: errStat     ! Error status of the operation
   character(*),                           intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

   integer                                                :: i,j           ! Generic counters     
   integer(IntKi)                                         :: errStat2        ! Error status of the operation (secondary error)
   character(len(errMsg))                                 :: errMsg2         ! Error message if ErrStat2 /= ErrID_None
   real(ReKi)                                             :: Cn_prime
   real(ReKi)                                             :: Cn1
   real(ReKi)                                             :: Cn2
   real(ReKi)                                             :: T_VL   
   real(ReKi)                                             :: Kalpha
   real(ReKi)                                             :: Kq
   real(ReKi)                                             :: q_cur
   real(ReKi)                                             :: X1
   real(ReKi)                                             :: X2
   real(ReKi)                                             :: Kprime_alpha
   real(ReKi)                                             :: Kprime_q
   real(ReKi)                                             :: Dp
   real(ReKi)                                             :: Cn_pot
   real(ReKi)                                             :: Cc_pot
   real(ReKi)                                             :: fprimeprime
   real(ReKi)                                             :: Df
   real(ReKi)                                             :: fprime
   real(ReKi)                                             :: Cn_v
   real(ReKi)                                             :: C_V 
   real(ReKi)                                             :: Cn_FS  
   real(ReKi)                                             :: alpha_e
   real(ReKi)                                             :: alpha0
   real(ReKi)                                             :: dalpha0
   real(ReKi)                                             :: eta_e
   real(ReKi)                                             :: St_sh
   real(ReKi)                                             :: Cl, Cd, Cd0
   integer                                                :: iNode
   
   iNode = 0
      
            ! We need to loop over all blade stations and update the states for each station
   if (OtherState%FirstPass) then
      do j = 1,p%numBlades
         do i = 1,p%nNodesPerBlade
            iNode = iNode + 1
            call GetSteadyOutputs(p%AFI_Params%AFInfo(p%AFIndx(i,j)), u%alpha(i,j), u%Re(i,j), y%Cl(i,j), y%Cd(i,j), errStat, errMsg)
            Cd0 = 0.0
            y%Cn(i,j) = y%Cl(i,j)*cos(u%alpha(i,j)) + (y%Cd(i,j)-Cd0)*sin(u%alpha(i,j))
            y%Cc(i,j) = y%Cl(i,j)*sin(u%alpha(i,j)) - (y%Cd(i,j)-Cd0)*cos(u%alpha(i,j))
            y%WriteOutput(iNode) = u%alpha(i,j)*180.0/pi
            y%WriteOutput(iNode+1) = y%Cn(i,j)
            y%WriteOutput(iNode+2) = y%Cc(i,j)
            y%WriteOutput(iNode+3) = y%Cm(i,j)
            y%WriteOutput(iNode+4) = y%Cl(i,j)
            y%WriteOutput(iNode+5) = y%Cd(i,j)
         end do
      end do
      
   else
      
      do j = 1,p%numBlades
         do i = 1,p%nNodesPerBlade
            
            iNode = iNode + 1
            Cd0 = 0.0
            
            call ComputeKelvinChain(i, j, u, p, xd, OtherState, Cn_prime, Cn1, Cn2, T_VL, St_sh, Kalpha, alpha_e, alpha0, dalpha0, eta_e, Kq, q_cur, X1, X2, &
                                    Kprime_alpha, Kprime_q, Dp, Cn_pot, Cc_pot, fprimeprime, Df, fprime, Cn_v, C_V, Cn_FS, errStat, errMsg )
            
            y%Cn(i,j) = Cn_FS + Cn_v
            
            y%Cc(i,j) = eta_e*Cc_pot*sqrt(fprimeprime) + Cn_v*tan(alpha_e)*(1-xd%tau_v(i,j))
            
            y%Cl(i,j) = y%Cn(i,j)*cos(u%alpha(i,j)) + y%Cc(i,j)*sin(u%alpha(i,j))
            y%Cd(i,j) = y%Cn(i,j)*sin(u%alpha(i,j)) - y%Cc(i,j)*cos(u%alpha(i,j)) + Cd0
            
            y%WriteOutput(iNode) = u%alpha(i,j)*180.0/pi
            y%WriteOutput(iNode+1) = y%Cn(i,j)
            y%WriteOutput(iNode+2) = y%Cc(i,j)
            y%WriteOutput(iNode+3) = y%Cm(i,j)
            y%WriteOutput(iNode+4) = y%Cl(i,j)
            y%WriteOutput(iNode+5) = y%Cd(i,j)
            
         end do
      end do
      
   
      
   end if
   
end subroutine UnsteadyAero_CalcOutput


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