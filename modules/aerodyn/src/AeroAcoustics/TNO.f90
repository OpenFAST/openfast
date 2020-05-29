MODULE TNO


   use NWTC_Library

   implicit none
   PUBLIC:: Pressure, f_int1, f_int2

   INTEGER, PARAMETER :: TNOKi = ReKi

   REAL (ReKi),PARAMETER :: Cnuk = 5.5
   REAL (ReKi),PARAMETER :: kappa = 0.41
   REAL (ReKi),PARAMETER :: Cmu = 0.09
   INTEGER(IntKi),PARAMETER :: limit = 5000

   !TNO variables
   REAL (ReKi) :: Omega_TNO ! NOTE: not a constant and used by function f_int1 and f_int2

   !atmosphere variables
   REAL (ReKi)  :: nu
   REAL (ReKi)  :: co
   REAL (ReKi)  :: rho

   ! Wavenumber variables
   REAL (ReKi)  :: k
   REAL (ReKi)  :: k1
   REAL (ReKi)  :: k3

   ! Blade params
   REAL (ReKi)  :: d99(2)
   REAL (ReKi)  :: Cf(2)
   REAL (ReKi)  :: edgevel(2)

   ! Airfoil
    REAL(ReKi)  :: Mach
    LOGICAL     :: ISSUCTION

   interface solve_qk61
      module procedure wrap_qk61
      module procedure wrap_dqk61
   end interface
      


contains

! Single precision wrapper for qk61 routine
subroutine wrap_qk61(Omega_TNO_in,f,a,b,answer,abserr,resabs,resasc)
   real(SiKi), intent(in   ) :: Omega_TNO_in
   real(SiKi), intent(in   ) :: a,b
   real(SiKi), intent(  out) :: answer
   real(SiKi), intent(in   ) :: abserr,resabs,resasc
   real(SiKi), external :: f
   Omega_TNO = real(Omega_TNO_in, TNOKi)
   call qk61(f,a,b,answer,abserr,resabs,resasc)
end subroutine wrap_qk61

! double precision wrapper for dqk61 routine
subroutine wrap_dqk61(Omega_TNO_in,f,a,b,answer,abserr,resabs,resasc)
   real(R8Ki), intent(in   ) :: Omega_TNO_in
   real(R8Ki), intent(in   ) :: a,b
   real(R8Ki), intent(  out) :: answer
   real(R8Ki), intent(in   ) :: abserr,resabs,resasc
   real(R8Ki), external :: f
   Omega_TNO = real(Omega_TNO_in, TNOKi)
   call dqk61(f,a,b,answer,abserr,resabs,resasc)
end subroutine wrap_dqk61


FUNCTION f_int1(x2)
   REAL(ReKi):: alpha
   REAL(ReKi):: alpha_gauss
   REAL(ReKi):: Cfin
   REAL(ReKi):: delta 
   REAL(ReKi):: dudx
   REAL(ReKi):: f_int1
   REAL(ReKi):: ke
   REAL(ReKi):: k1_hat
   REAL(ReKi):: k3_hat
   REAL(ReKi):: kT
   REAL(ReKi):: L
   REAL(ReKi):: Nut
   REAL(ReKi):: phi22
   REAL(ReKi):: phim
   REAL(ReKi):: ums
   REAL(ReKi):: u_star
   REAL(ReKi):: U
   REAL(ReKi):: Uc
   REAL(ReKi):: Uo
   REAL(ReKi):: W
   REAL(ReKi), intent(in) :: x2
   
   ! changed and being multiplied with edge velocity taken from xfoil output
   ! Uo=Mach*co issuction use edgevel(1) 
   
   !constants from xfoil
   if (ISSUCTION) then
      alpha = 0.45 ! = 0.3 pressure, = 0.45 suction
      Cfin = Cf(1)
      delta = d99(1)
      Uo=Mach*co*edgevel(1)
   else
      alpha = 0.30
      Cfin = Cf(2)
      delta = d99(2)
      Uo=Mach*co*edgevel(2)
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
   REAL (ReKi), intent(in)  :: k1
   REAL (ReKi) :: f_int2
   f_int2 = Omega_TNO/co/k1*Pressure(k1)
   RETURN 
END FUNCTION f_int2


FUNCTION Pressure(k1_in)
    ! Variables
   REAL(ReKi)  :: a,b,answer
   REAL(ReKi)  :: omega
   REAL(ReKi)  :: abserr,resabs,resasc
   REAL(ReKi)  :: k1_in
   real(ReKi)  :: Pressure

   ! Set variables used in f_int1
   k1 = k1_in

    a = 0.0_ReKi !1e-4*d99(1)
    IF (ISSUCTION)THEN
        b = d99(1)
    ELSE
        b = d99(2)
    ENDIF

    k3 = 0.
    k= sqrt(k1**2+k3**2)

    CALL qk61(f_int1,a,b,answer,abserr,resabs,resasc)
               
    Pressure = 4.*rho**2*k1**2./(k1**2.+k3**2.)*answer
               
   RETURN
END FUNCTION Pressure


END MODULE TNO
