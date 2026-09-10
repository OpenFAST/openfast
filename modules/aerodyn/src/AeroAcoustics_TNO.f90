MODULE TNO


   use NWTC_Library  ! ReKi, DBKi, R8Ki
   !use NWTC_SLATEC   ! slatec_qk61: bjj: this has been replaced with a local qk61_ctx that passes a context (TNO_ContextType) to the integrand. See comments below.

   implicit none
   PRIVATE
   PUBLIC :: SPL_integrate

   INTEGER,       PARAMETER :: TNOKi = ReKi

   REAL (TNOKi),  PARAMETER :: Cnuk = 5.5
   REAL (TNOKi),  PARAMETER :: kappa = 0.41
   REAL (TNOKi),  PARAMETER :: Cmu = 0.09
!   INTEGER(IntKi),PARAMETER :: limit = 5000

   !> Everything the TNO integrands need in order to be evaluated.
   !!
   !! This type replaces what used to be a set of module-level (and therefore implicitly SAVEd) variables. The reason
   !! those existed is that the SLATEC integrator qk61 declares its integrand as "external f" and calls it as f(x), so
   !! there was nowhere to thread the flow conditions through. Passing this context as an explicit argument instead
   !! makes the module reentrant: every activation of SPL_integrate owns its own context, so concurrent calls (e.g.
   !! from an OpenMP region, or two rotors evaluated in parallel) can no longer overwrite each other's state
   !! mid-integration.
   type :: TNO_ContextType
      ! frequency
      real(TNOKi) :: Omega     = 0.0_TNOKi     !< radian frequency
      ! atmosphere
      real(TNOKi) :: nu        = 0.0_TNOKi     !< kinematic viscosity
      real(TNOKi) :: co        = 0.0_TNOKi     !< speed of sound
      real(TNOKi) :: rho       = 0.0_TNOKi     !< air density
      ! airfoil
      real(TNOKi) :: Mach      = 0.0_TNOKi     !< Mach number of this blade node
      logical     :: IsSuction = .false.       !< .true. for the suction side, .false. for the pressure side
      ! blade-node boundary-layer properties; index 1 = suction side, index 2 = pressure side
      real(TNOKi) :: d99(2)     = 0.0_TNOKi    !< boundary layer thickness
      real(TNOKi) :: Cf(2)      = 0.0_TNOKi    !< skin friction coefficient
      real(TNOKi) :: edgevel(2) = 0.0_TNOKi    !< edge velocity ratio
      ! Wavenumbers. Unlike the fields above (which are fixed for the whole integration), these are updated by
      ! Pressure() at each point of the OUTER integration over k1 and then read by f_int1 during the INNER
      ! integration over x2. That mid-integration coupling is why the context is passed as intent(inout).
      real(TNOKi) :: k1        = 0.0_TNOKi
      real(TNOKi) :: k3        = 0.0_TNOKi
      real(TNOKi) :: k         = 0.0_TNOKi
   end type TNO_ContextType

   abstract interface
      !> Integrand accepted by qk61_ctx: a function of the integration variable plus a context.
      real(TNOKi) function TNO_Integrand(x, ctx)
         import :: TNOKi, TNO_ContextType
         real(TNOKi),           intent(in   ) :: x
         type(TNO_ContextType), intent(inout) :: ctx
      end function TNO_Integrand
   end interface

   ! ------------------------------------------------------------------------------------------------------------------
   ! Data for the 61-point Gauss-Kronrod rule, transcribed from modules/nwtc-library/src/NetLib/slatec/dqk61.f
   ! (SLATEC/QUADPACK; weights and abscissae evaluated with 80-decimal-digit arithmetic by L. W. Fullerton, Bell Labs,
   ! Nov. 1981). The abscissae and weights are given for the interval (-1,1); by symmetry only the non-negative
   ! abscissae and their corresponding weights are stored.
   !
   ! The literals below carry dqk61's full precision and are converted to TNOKi, so a single-precision build rounds to
   ! exactly the values hard-coded in qk61.f and a double-precision build matches dqk61.f. That is what lets this
   ! routine reproduce the previous results bit-for-bit in both precisions.
   !
   !   xgk - abscissae of the 61-point Kronrod rule.
   !         xgk(2), xgk(4), ... are the abscissae of the 30-point Gauss rule.
   !         xgk(1), xgk(3), ... are the optimally added abscissae.
   !   wgk - weights of the 61-point Kronrod rule.
   !   wg  - weights of the 30-point Gauss rule.
   ! ------------------------------------------------------------------------------------------------------------------
   real(TNOKi), parameter :: xgk(31) = (/ &
        0.999484410050490637571325895705811_TNOKi, 0.996893484074649540271630050918695_TNOKi, &
        0.991630996870404594858628366109486_TNOKi, 0.983668123279747209970032581605663_TNOKi, &
        0.973116322501126268374693868423707_TNOKi, 0.960021864968307512216871025581798_TNOKi, &
        0.944374444748559979415831324037439_TNOKi, 0.926200047429274325879324277080474_TNOKi, &
        0.905573307699907798546522558925958_TNOKi, 0.882560535792052681543116462530226_TNOKi, &
        0.857205233546061098958658510658944_TNOKi, 0.829565762382768397442898119732502_TNOKi, &
        0.799727835821839083013668942322683_TNOKi, 0.767777432104826194917977340974503_TNOKi, &
        0.733790062453226804726171131369528_TNOKi, 0.697850494793315796932292388026640_TNOKi, &
        0.660061064126626961370053668149271_TNOKi, 0.620526182989242861140477556431189_TNOKi, &
        0.579345235826361691756024932172540_TNOKi, 0.536624148142019899264169793311073_TNOKi, &
        0.492480467861778574993693061207709_TNOKi, 0.447033769538089176780609900322854_TNOKi, &
        0.400401254830394392535476211542661_TNOKi, 0.352704725530878113471037207089374_TNOKi, &
        0.304073202273625077372677107199257_TNOKi, 0.254636926167889846439805129817805_TNOKi, &
        0.204525116682309891438957671002025_TNOKi, 0.153869913608583546963794672743256_TNOKi, &
        0.102806937966737030147096751318001_TNOKi, 0.051471842555317695833025213166723_TNOKi, &
        0.000000000000000000000000000000000_TNOKi /)

   real(TNOKi), parameter :: wgk(31) = (/ &
        0.001389013698677007624551591226760_TNOKi, 0.003890461127099884051267201844516_TNOKi, &
        0.006630703915931292173319826369750_TNOKi, 0.009273279659517763428441146892024_TNOKi, &
        0.011823015253496341742232898853251_TNOKi, 0.014369729507045804812451432443580_TNOKi, &
        0.016920889189053272627572289420322_TNOKi, 0.019414141193942381173408951050128_TNOKi, &
        0.021828035821609192297167485738339_TNOKi, 0.024191162078080601365686370725232_TNOKi, &
        0.026509954882333101610601709335075_TNOKi, 0.028754048765041292843978785354334_TNOKi, &
        0.030907257562387762472884252943092_TNOKi, 0.032981447057483726031814191016854_TNOKi, &
        0.034979338028060024137499670731468_TNOKi, 0.036882364651821229223911065617136_TNOKi, &
        0.038678945624727592950348651532281_TNOKi, 0.040374538951535959111995279752468_TNOKi, &
        0.041969810215164246147147541285970_TNOKi, 0.043452539701356069316831728117073_TNOKi, &
        0.044814800133162663192355551616723_TNOKi, 0.046059238271006988116271735559374_TNOKi, &
        0.047185546569299153945261478181099_TNOKi, 0.048185861757087129140779492298305_TNOKi, &
        0.049055434555029778887528165367238_TNOKi, 0.049795683427074206357811569379942_TNOKi, &
        0.050405921402782346840893085653585_TNOKi, 0.050881795898749606492297473049805_TNOKi, &
        0.051221547849258772170656282604944_TNOKi, 0.051426128537459025933862879215781_TNOKi, &
        0.051494729429451567558340433647099_TNOKi /)

   real(TNOKi), parameter :: wg(15) = (/ &
        0.007968192496166605615465883474674_TNOKi, 0.018466468311090959142302131912047_TNOKi, &
        0.028784707883323369349719179611292_TNOKi, 0.038799192569627049596801936446348_TNOKi, &
        0.048402672830594052902938140422808_TNOKi, 0.057493156217619066481721689402056_TNOKi, &
        0.065974229882180495128128515115962_TNOKi, 0.073755974737705206268243850022191_TNOKi, &
        0.080755895229420215354694938460530_TNOKi, 0.086899787201082979802387530715126_TNOKi, &
        0.092122522237786128717632707087619_TNOKi, 0.096368737174644259639468626351810_TNOKi, &
        0.099593420586795267062780282103569_TNOKi, 0.101762389748405504596428952168554_TNOKi, &
        0.102852652893558840341285636705415_TNOKi /)


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
   
   real(TNOKi)               :: answer             !< value returned from qk61_ctx, NOTE the typing

   !> All state needed by the integrands. This is a local variable (not module data), which is what makes this
   !! routine safe to call from more than one thread at a time.
   type(TNO_ContextType)     :: ctx

   ! local variables that are ignored
   real(TNOKi) :: abserr,resabs,resasc             !< accuracy estimates and residuals. Currently ignored

   ! Set context values from input
   ctx%IsSuction = ISSUCTION
   ctx%Omega     = real(Omega,TNOKi)
   
   ! Mach number of segment
   ctx%Mach      = real(Mach,TNOKi)
   
   ! Atmospheric values
   ctx%co        = real(SpdSound,  TNOKi)
   ctx%rho       = real(AirDens,   TNOKi)
   ctx%nu        = real(KinVisc,   TNOKi)
   
   ! Blade node values
   ctx%Cf        = real(Cfall,     TNOKi)
   ctx%d99       = real(d99all,    TNOKi)
   ctx%edgevel   = real(ABS(EdgeVelAll),TNOKi)

   call qk61_ctx(f_int2,limits(1),limits(2),ctx,answer,abserr,resabs,resasc)
   integrand = real( answer, ReKi )

end function SPL_integrate


!==================================================================================================================================
!> Integrate FUNC over (a,b) using the 61-point Gauss-Kronrod rule, passing CTX through to the integrand.
!!
!! This is a direct port of the SLATEC routine qk61 (see modules/nwtc-library/src/NetLib/slatec/qk61.f), changed only
!! so that the integrand takes a context argument. The accumulation order is deliberately identical to the original so
!! that the computed RESULT is bit-for-bit the same as before this port.
!!
!! Two things the original could not offer:
!!   1. an explicit interface for the integrand (qk61 used "external f", an implicit interface), and
!!   2. a way to give the integrand its parameters other than global data.
!!
!! Declared RECURSIVE because the TNO model nests two integrations: the outer integral over k1 calls f_int2, which
!! calls Pressure, which integrates over x2. The original relied on per-file compiler flags (-frecursive /
!! -assume recursion, see modules/nwtc-library/CMakeLists.txt) to make qk61's locals automatic; stating RECURSIVE in
!! the source expresses that requirement in the code instead of in the build system.
recursive subroutine qk61_ctx(func, a, b, ctx, result, abserr, resabs, resasc)
   procedure(TNO_Integrand)             :: func     !< integrand, evaluated as func(x, ctx)
   real(TNOKi),           intent(in   ) :: a        !< lower limit of integration
   real(TNOKi),           intent(in   ) :: b        !< upper limit of integration
   type(TNO_ContextType), intent(inout) :: ctx      !< context handed to the integrand
   real(TNOKi),           intent(  out) :: result   !< approximation to the integral, from the 61-point Kronrod rule
   real(TNOKi),           intent(  out) :: abserr   !< estimate of the modulus of the absolute error
   real(TNOKi),           intent(  out) :: resabs   !< approximation to the integral of abs(func)
   real(TNOKi),           intent(  out) :: resasc   !< approximation to the integral of abs(func-i/(b-a))

   ! Local variables (names kept from qk61 to keep this reviewable against the original)
   real(TNOKi) :: absc                     ! abscissa
   real(TNOKi) :: centr                    ! mid point of the interval
   real(TNOKi) :: dhlgth                   ! abs(hlgth)
   real(TNOKi) :: epmach                   ! the largest relative spacing
   real(TNOKi) :: fc                       ! function value at the mid point
   real(TNOKi) :: fsum
   real(TNOKi) :: fval1, fval2             ! function values
   real(TNOKi) :: fv1(30), fv2(30)
   real(TNOKi) :: hlgth                    ! half-length of the interval
   real(TNOKi) :: resg                     ! result of the 30-point Gauss rule
   real(TNOKi) :: resk                     ! result of the 61-point Kronrod rule
   real(TNOKi) :: reskh                    ! approximation to the mean value of func over (a,b), i.e. to i/(b-a)
   real(TNOKi) :: uflow                    ! the smallest positive magnitude
   integer     :: j, jtw, jtwm1

   ! r1mach(4)/d1mach(4) is the largest relative spacing and r1mach(1)/d1mach(1) the smallest positive magnitude;
   ! EPSILON and TINY are the standard intrinsics for exactly those quantities, so the SLATEC *1mach dependency is
   ! not needed here.
   epmach = epsilon(1.0_TNOKi)
   uflow  = tiny(1.0_TNOKi)

   centr  = 0.5_TNOKi*(b+a)
   hlgth  = 0.5_TNOKi*(b-a)
   dhlgth = abs(hlgth)

   ! Compute the 61-point Kronrod approximation to the integral, and estimate the absolute error.
   resg   = 0.0_TNOKi
   fc     = func(centr, ctx)
   resk   = wgk(31)*fc
   resabs = abs(resk)

   do j = 1,15
      jtw       = j*2
      absc      = hlgth*xgk(jtw)
      fval1     = func(centr-absc, ctx)
      fval2     = func(centr+absc, ctx)
      fv1(jtw)  = fval1
      fv2(jtw)  = fval2
      fsum      = fval1+fval2
      resg      = resg+wg(j)*fsum
      resk      = resk+wgk(jtw)*fsum
      resabs    = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   end do

   do j = 1,15
      jtwm1       = j*2-1
      absc        = hlgth*xgk(jtwm1)
      fval1       = func(centr-absc, ctx)
      fval2       = func(centr+absc, ctx)
      fv1(jtwm1)  = fval1
      fv2(jtwm1)  = fval2
      fsum        = fval1+fval2
      resk        = resk+wgk(jtwm1)*fsum
      resabs      = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   end do

   reskh  = resk*0.5_TNOKi
   resasc = wgk(31)*abs(fc-reskh)
   do j = 1,30
      resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   end do

   result = resk*hlgth
   resabs = resabs*dhlgth
   resasc = resasc*dhlgth
   abserr = abs((resk-resg)*hlgth)
   if (resasc /= 0.0_TNOKi .and. abserr /= 0.0_TNOKi) &
      abserr = resasc*min(1.0_TNOKi,(200.0_TNOKi*abserr/resasc)**1.5_TNOKi)
   if (resabs > uflow/(50.0_TNOKi*epmach)) &
      abserr = max((epmach*50.0_TNOKi)*resabs,abserr)

end subroutine qk61_ctx
!==================================================================================================================================


FUNCTION f_int1(x2, ctx) result(f)
   REAL(TNOKi), intent(in) :: x2
   type(TNO_ContextType), intent(inout) :: ctx
   REAL(TNOKi):: f

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
   
   ! changed and being multiplied with edge velocity taken from xfoil output
   ! Uo=ctx%Mach*ctx%co ctx%IsSuction use ctx%edgevel(1) 
   
   !constants from xfoil
   if (ctx%IsSuction) then
      alpha = 0.45 ! = 0.3 pressure, = 0.45 suction
      Cfin = ctx%Cf(1)
      delta = ctx%d99(1)
      Uo=ctx%Mach*ctx%co*ctx%edgevel(1)
   else
      alpha = 0.30
      Cfin = ctx%Cf(2)
      delta = ctx%d99(2)
      Uo=ctx%Mach*ctx%co*ctx%edgevel(2)
   endif
   ! bjj: Bail out (contributing nothing to the integral) instead of producing NaN/Inf or killing the program.
   !  - Cf <= 0 used to execute a bare "stop", which aborts without OpenFAST's error handling and without closing
   !    output files. TBLTE_TNO already skips the side whose Cf is non-positive, so this is a backstop.
   !  - delta (d99) can be zero for unconverged entries in the pre-tabulated boundary-layer files. That makes L = 0/0
   !    and pi*x2/delta infinite below.
   !  - u_star is zero when the edge velocity ratio or the Mach number is zero, and log(u_star*x2/nu) is then -Inf,
   !    so U evaluates to 0*(-Inf) = NaN.
   !  - L is zero at x2 = 0 (the lower integration limit), which makes ke = sqrt(pi)/L infinite.
   if (Cfin .le. 0. .or. delta .le. 0. .or. x2 .le. 0.) then
      f = 0.
      RETURN
   endif
   
   u_star = Uo*sqrt(Cfin/2.)
   
   if (u_star .le. 0.) then
      f = 0.
      RETURN
   endif
   
   L = 0.085*delta*tanh(kappa*x2/(0.085*delta))
   
   if (L .le. 0.) then
      f = 0.
      RETURN
   endif
   
   if (x2 .gt. delta)then
      U = Uo
      dudx = 0.
      f = 0.
      RETURN
   else
      W = 1.-cos(pi*x2/delta);
      U = u_star*(1./kappa*log(u_star*x2/ctx%nu) +Cnuk+ (Uo/u_star-1./kappa*log(u_star*delta/ctx%nu)-Cnuk)*0.5*W)
      dudx = u_star*(1./(kappa*x2)+(Uo/u_star-1./kappa*log(u_star*delta/ctx%nu)-Cnuk)* &
             0.5*(pi/delta)*sin(pi*x2/delta))
   endif
          
   ke=sqrt(pi)/L*0.4213560764 !gamma(5./6.)/gamma(1./3.)
   k1_hat = ctx%k1/ke
   k3_hat = ctx%k3/ke
   
   Nut = (L*kappa)**2.*abs(dudx)
   kT = sqrt((Nut*dudx)**2./Cmu)
   ums = alpha*kT
   
   Uc = 0.7*U
   alpha_gauss = 0.05*Uc/L
   
   phim = 1./(alpha_gauss*sqrt(pi))*exp(-((ctx%Omega-Uc*ctx%k1)/alpha_gauss)**2.)
   phi22 = 4./9./pi*1/ke**2.*(k1_hat**2.+k3_hat**2.)/(1.+k1_hat**2.+k3_hat**2.)**(7./3.)
   f = L*ums*(dudx)**2*phi22*phim*exp(-2*abs(ctx%k)*x2)
   
   RETURN
END FUNCTION f_int1


FUNCTION f_int2(k1_in, ctx) result(f)  ! changed name from 'int2' to avoid conflicts with intrinsic of same name
   REAL (TNOKi), intent(in)  :: k1_in
   type(TNO_ContextType), intent(inout) :: ctx
   REAL (TNOKi) :: f
   
   ! bjj: The lower integration limit passed to qk61_ctx is exactly zero. The 61-point Gauss-Kronrod rule does not
   ! evaluate the integrand at the interval end points, so k1_in is not zero in practice, but guard the 1/k1_in here
   ! (and the k1**2/(k1**2+k3**2) = 0/0 in Pressure) so this does not depend on the internals of the quadrature rule.
   if (k1_in .le. 0.0_TNOKi) then
      f = 0.0_TNOKi
      RETURN
   endif
   
   f = ctx%Omega/ctx%co/k1_in*Pressure(k1_in, ctx)
   RETURN 
END FUNCTION f_int2


FUNCTION Pressure(k1_in, ctx) result(P)
    ! Variables
   REAL(TNOKi), intent(in)  :: k1_in
   type(TNO_ContextType), intent(inout) :: ctx
   real(TNOKi)  :: P

   REAL(TNOKi)  :: a,b,answer
   REAL(TNOKi)  :: abserr,resabs,resasc

   ! Set wavenumbers used in f_int1
   ctx%k1 = k1_in

    a = 0.0_TNOKi !1e-4*ctx%d99(1)
    IF (ctx%IsSuction)THEN
        b = ctx%d99(1)
    ELSE
        b = ctx%d99(2)
    ENDIF

    ctx%k3 = 0.
    ctx%k  = sqrt(ctx%k1**2+ctx%k3**2)

    CALL qk61_ctx(f_int1,a,b,ctx,answer,abserr,resabs,resasc)
               
    P = 4.0_TNOKi*ctx%rho**2 * ctx%k1**2 / (ctx%k1**2 + ctx%k3**2)*answer
               
   RETURN
END FUNCTION Pressure


END MODULE TNO
