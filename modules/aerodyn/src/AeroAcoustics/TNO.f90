MODULE TNO


   use NWTC_Library  ! ReKi, DBKi, R8Ki
   use NWTC_SLATEC   ! slatec_qk61 -- which is all that is in that library right now.

   implicit none
   PUBLIC :: SPL_integrate

   INTEGER,       PARAMETER :: TNOKi = ReKi

   REAL (TNOKi),  PARAMETER :: Cnuk = 5.5
   REAL (TNOKi),  PARAMETER :: kappa = 0.41
   REAL (TNOKi),  PARAMETER :: Cmu = 0.09
!   INTEGER(IntKi),PARAMETER :: limit = 5000

   !TNO variables
   REAL (TNOKi) :: Omega_TNO ! NOTE: not a constant and used by function f_int1 and f_int2

   !atmosphere variables
   REAL (TNOKi)  :: nu
   REAL (TNOKi)  :: co
   REAL (TNOKi)  :: rho

   ! Wavenumber variables
   REAL (TNOKi)  :: k
   REAL (TNOKi)  :: k1
   REAL (TNOKi)  :: k3

   ! Blade params
   REAL (TNOKi)  :: d99(2)
   REAL (TNOKi)  :: Cf(2)
   REAL (TNOKi)  :: edgevel(2)

   ! Airfoil
    REAL(TNOKi)  :: Mach_TNO
    LOGICAL      :: ISSUCTION_TNO


contains

!> Solve the spl generated at this location and frequency
function SPL_integrate(Omega,limits,ISSUCTION,   &
               Mach,SpdSound,AirDens,KinVisc,      &
               Cfall,d99all,EdgeVelAll) result(integrand)
   real(ReKi), intent(in   ) :: Omega              !< frequency
   real(ReKi), intent(in   ) :: limits(2)          !< integration limits
   logical,    intent(in   ) :: ISSUCTION          !< Is it the suction edge
   real(ReKi), intent(in   ) :: Mach               !< Mach number
   real(ReKi), intent(in   ) :: SpdSound           !< Speed of sound
   real(ReKi), intent(in   ) :: AirDens            !< Air density
   real(ReKi), intent(in   ) :: KinVisc            !< Kinetic air viscosity
   real(ReKi), intent(in   ) :: Cfall(2)           !< Skin friction coefficient   (-)
   real(ReKi), intent(in   ) :: d99all(2)          !< 
   real(ReKi), intent(in   ) :: EdgeVelAll(2)      !< 
   real(ReKi)                :: integrand          !< integrand result
   
   real(TNOKi)               :: answer             !< value returned from qk61, NOTE the typing

   ! local variables that are ignored
   real(TNOKi) :: abserr,resabs,resasc             !< accuracy estimates and residuals. Currently ignored

   ! Set module values from input
   ISSUCTION_TNO  = ISSUCTION
   Omega_TNO      = real(Omega,TNOKi)
   ! Mach number of segment
   Mach_TNO       = real(Mach,TNOKi)
   ! Atmospheric values
   co       = real(SpdSound,  TNOKi)
   rho      = real(AirDens,   TNOKi)
   nu       = real(KinVisc,   TNOKi)
   ! Blade node values
   Cf       = real(Cfall,     TNOKi)
   d99      = real(d99all,    TNOKi)
   edgevel  = real(ABS(EdgeVelAll),TNOKi)

   call slatec_qk61(f_int2,limits(1),limits(2),answer,abserr,resabs,resasc)
   integrand = real( answer, ReKi )

end function SPL_integrate



FUNCTION f_int1(x2)
   REAL(TNOKi):: alpha
   REAL(TNOKi):: alpha_gauss
   REAL(TNOKi):: Cfin
   REAL(TNOKi):: delta 
   REAL(TNOKi):: dudx
   REAL(TNOKi):: ke
   REAL(TNOKi):: k1_hat
   REAL(TNOKi):: k3_hat
   REAL(TNOKi):: kT
   REAL(TNOKi):: L
   REAL(TNOKi):: Nut
   REAL(TNOKi):: phi22
   REAL(TNOKi):: phim
   REAL(TNOKi):: ums
   REAL(TNOKi):: u_star
   REAL(TNOKi):: U
   REAL(TNOKi):: Uc
   REAL(TNOKi):: Uo
   REAL(TNOKi):: W
   REAL(TNOKi), intent(in) :: x2
   REAL(TNOKi):: f_int1
   
   ! changed and being multiplied with edge velocity taken from xfoil output
   ! Uo=Mach_TNO*co ISSUCTION_TNO use edgevel(1) 
   
   !constants from xfoil
   if (ISSUCTION_TNO) then
      alpha = 0.45 ! = 0.3 pressure, = 0.45 suction
      Cfin = Cf(1)
      delta = d99(1)
      Uo=Mach_TNO*co*edgevel(1)
   else
      alpha = 0.30
      Cfin = Cf(2)
      delta = d99(2)
      Uo=Mach_TNO*co*edgevel(2)
   endif
   if (Cfin .le. 0.) then
      write(*,*) 'Cf is less than zero, Cf = ',Cfin
      stop
   endif
   u_star = Uo*sqrt(Cfin/2.)
   
   L = 0.085*delta*tanh(kappa*x2/(0.085*delta))
   
   if (x2 .gt. delta)then
      U = Uo
      dudx = 0.
      f_int1 = 0.
      RETURN
   else
      W = 1.-cos(pi*x2/delta);
      U = u_star*(1./kappa*log(u_star*x2/nu) +Cnuk+ (Uo/u_star-1./kappa*log(u_star*delta/nu)-Cnuk)*0.5*W)
      dudx = u_star*(1./(kappa*x2)+(Uo/u_star-1./kappa*log(u_star*delta/nu)-Cnuk)* &
             0.5*(pi/delta)*sin(pi*x2/delta))
   endif
          
   ke=sqrt(pi)/L*0.4213560764 !gamma(5./6.)/gamma(1./3.)
   k1_hat = k1/ke
   k3_hat = k3/ke
   
   Nut = (L*kappa)**2.*abs(dudx)
   kT = sqrt((Nut*dudx)**2./Cmu)
   ums = alpha*kT
   
   Uc = 0.7*U
   alpha_gauss = 0.05*Uc/L
   
   phim = 1./(alpha_gauss*sqrt(pi))*exp(-((Omega_TNO-Uc*k1)/alpha_gauss)**2.)
   phi22 = 4./9./pi*1/ke**2.*(k1_hat**2.+k3_hat**2.)/(1.+k1_hat**2.+k3_hat**2.)**(7./3.)
   f_int1 = L*ums*(dudx)**2*phi22*phim*exp(-2*abs(k)*x2)
   
   RETURN
END FUNCTION f_int1


FUNCTION f_int2(k1)  ! changed name from 'int2' to avoid conflicts with intrinsic of same name
   REAL (TNOKi), intent(in)  :: k1
   REAL (TNOKi) :: f_int2
   f_int2 = Omega_TNO/co/k1*Pressure(k1)
   RETURN 
END FUNCTION f_int2


FUNCTION Pressure(k1_in)
    ! Variables
   REAL(TNOKi)  :: a,b,answer
   REAL(TNOKi)  :: omega
   REAL(TNOKi)  :: abserr,resabs,resasc
   REAL(TNOKi)  :: k1_in
   real(TNOKi)  :: Pressure

   ! Set variables used in f_int1
   k1 = k1_in

    a = 0.0_TNOKi !1e-4*d99(1)
    IF (ISSUCTION_TNO)THEN
        b = d99(1)
    ELSE
        b = d99(2)
    ENDIF

    k3 = 0.
    k= sqrt(k1**2+k3**2)

    CALL slatec_qk61(f_int1,a,b,answer,abserr,resabs,resasc)
               
    Pressure = 4.*rho**2*k1**2./(k1**2.+k3**2.)*answer
               
   RETURN
END FUNCTION Pressure


END MODULE TNO
