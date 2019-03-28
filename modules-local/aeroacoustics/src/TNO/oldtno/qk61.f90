subroutine qk61 ( f, a, b, result, abserr, resabs, resasc ) 

!*****************************************************************************80
!
!! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!                    result is computed by applying the 61-point
!                    Kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point Gauss rule (resg).
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 4 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 4 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 30-point Gauss rule
!           resk   - result of the 61-point Kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absc
  real ( kind = 4 ) abserr
  real ( kind = 4 ) b
  real ( kind = 4 ) centr
  real ( kind = 4 ) dhlgth
  real ( kind = 4 ), external :: f
  real ( kind = 4 ) fc
  real ( kind = 4 ) fsum
  real ( kind = 4 ) fval1
  real ( kind = 4 ) fval2
  real ( kind = 4 ) fv1(30)
  real ( kind = 4 ) fv2(30)
  real ( kind = 4 ) hlgth
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtw
  integer ( kind = 4 ) jtwm1
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) resg
  real ( kind = 4 ) resk
  real ( kind = 4 ) reskh
  real ( kind = 4 ) result
  real ( kind = 4 ) wg(15)
  real ( kind = 4 ) wgk(31)
  real ( kind = 4 ) xgk(31)
!
!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.
!
!           xgk   - abscissae of the 61-point Kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   Gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point Gauss rule
!
!           wgk   - weights of the 61-point Kronrod rule
!
!           wg    - weigths of the 30-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
     xgk(9),xgk(10)/ &
       9.994844100504906E-01,     9.968934840746495E-01, &
       9.916309968704046E-01,     9.836681232797472E-01, &
       9.731163225011263E-01,     9.600218649683075E-01, &
       9.443744447485600E-01,     9.262000474292743E-01, &
       9.055733076999078E-01,     8.825605357920527E-01/
  data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
    xgk(18),xgk(19),xgk(20)/ &
       8.572052335460611E-01,     8.295657623827684E-01, &
       7.997278358218391E-01,     7.677774321048262E-01, &
       7.337900624532268E-01,     6.978504947933158E-01, &
       6.600610641266270E-01,     6.205261829892429E-01, &
       5.793452358263617E-01,     5.366241481420199E-01/
  data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
    xgk(28),xgk(29),xgk(30),xgk(31)/ &
       4.924804678617786E-01,     4.470337695380892E-01, &
       4.004012548303944E-01,     3.527047255308781E-01, &
       3.040732022736251E-01,     2.546369261678898E-01, &
       2.045251166823099E-01,     1.538699136085835E-01, &
       1.028069379667370E-01,     5.147184255531770E-02, &
       0.0E+00                   /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10)/ &
       1.389013698677008E-03,     3.890461127099884E-03, &
       6.630703915931292E-03,     9.273279659517763E-03, &
       1.182301525349634E-02,     1.436972950704580E-02, &
       1.692088918905327E-02,     1.941414119394238E-02, &
       2.182803582160919E-02,     2.419116207808060E-02/
  data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
    wgk(18),wgk(19),wgk(20)/ &
       2.650995488233310E-02,     2.875404876504129E-02, &
       3.090725756238776E-02,     3.298144705748373E-02, &
       3.497933802806002E-02,     3.688236465182123E-02, &
       3.867894562472759E-02,     4.037453895153596E-02, &
       4.196981021516425E-02,     4.345253970135607E-02/
  data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
    wgk(28),wgk(29),wgk(30),wgk(31)/ &
       4.481480013316266E-02,     4.605923827100699E-02, &
       4.718554656929915E-02,     4.818586175708713E-02, &
       4.905543455502978E-02,     4.979568342707421E-02, &
       5.040592140278235E-02,     5.088179589874961E-02, &
       5.122154784925877E-02,     5.142612853745903E-02, &
       5.149472942945157E-02/
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       7.968192496166606E-03,     1.846646831109096E-02, &
       2.878470788332337E-02,     3.879919256962705E-02, &
       4.840267283059405E-02,     5.749315621761907E-02, &
       6.597422988218050E-02,     7.375597473770521E-02/
  data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
       8.075589522942022E-02,     8.689978720108298E-02, &
       9.212252223778613E-02,     9.636873717464426E-02, &
       9.959342058679527E-02,     1.017623897484055E-01, &
       1.028526528935588E-01/

  centr = 5.0E-01*(b+a)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!	print*, 'ab',a,b
!
!  Compute the 61-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!

  resg = 0.0E+00
  fc = f(centr)
  resk = wgk(31)*fc
  resabs = abs(resk)

!        print*, 'SLEEPING '
!	CALL SLEEP(50)
  do j = 1, 15
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 15
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0E-01
  resasc = wgk(31)*abs(fc-reskh)

  do j = 1, 30
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00 .and. abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / (5.0E+01* epsilon ( resabs ) )) then
    abserr = max ( ( epsilon ( resabs ) *5.0E+01)*resabs, abserr )
  end if
!	print*,'result',result
  return
end
