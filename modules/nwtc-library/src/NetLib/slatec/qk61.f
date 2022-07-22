      subroutine qk61(f,a,b,result,abserr,resabs,resasc)
c***begin prologue  qk61
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  61-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f over (a,b) with error
c                           estimate
c                       j = integral of dabs(f) over (a,b)
c***description
c
c        integration rule
c        standard fortran subroutine
c        real version
c
c
c        parameters
c         on entry
c           f      - real
c                    function subprogram defining the integrand
c                    function f(x). the actual name for f needs to be
c                    declared e x t e r n a l in the calling program.
c
c           a      - real
c                    lower limit of integration
c
c           b      - real
c                    upper limit of integration
c
c         on return
c           result - real
c                    approximation to the integral i
c                    result is computed by applying the 61-point
c                    kronrod rule (resk) obtained by optimal addition of
c                    abscissae to the 30-point gauss rule (resg).
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed dabs(i-result)
c
c           resabs - real
c                    approximation to the integral j
c
c           resasc - real
c                    approximation to the integral of dabs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk61
c
      real a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,
     *  fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,r1mach,uflow,
     *  wg,wgk,xgk
      integer j,jtw,jtwm1
      external f
c
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
c
c           the abscissae and weights are given for the
c           interval (-1,1). because of symmetry only the positive
c           abscissae and their corresponding weights are given.
c
c           xgk   - abscissae of the 61-point kronrod rule
c                   xgk(2), xgk(4)  ... abscissae of the 30-point
c                   gauss rule
c                   xgk(1), xgk(3)  ... optimally added abscissae
c                   to the 30-point gauss rule
c
c           wgk   - weights of the 61-point kronrod rule
c
c           wg    - weigths of the 30-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),
     *   xgk(9),xgk(10)/
     *     0.9994844100504906e+00,     0.9968934840746495e+00,
     *     0.9916309968704046e+00,     0.9836681232797472e+00,
     *     0.9731163225011263e+00,     0.9600218649683075e+00,
     *     0.9443744447485600e+00,     0.9262000474292743e+00,
     *     0.9055733076999078e+00,     0.8825605357920527e+00/
      data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),
     *  xgk(17),xgk(18),xgk(19),xgk(20)/
     *     0.8572052335460611e+00,     0.8295657623827684e+00,
     *     0.7997278358218391e+00,     0.7677774321048262e+00,
     *     0.7337900624532268e+00,     0.6978504947933158e+00,
     *     0.6600610641266270e+00,     0.6205261829892429e+00,
     *     0.5793452358263617e+00,     0.5366241481420199e+00/
      data xgk(21),xgk(22),xgk(23),xgk(24),
     *  xgk(25),xgk(26),xgk(27),xgk(28),xgk(29),xgk(30),xgk(31)/
     *     0.4924804678617786e+00,     0.4470337695380892e+00,
     *     0.4004012548303944e+00,     0.3527047255308781e+00,
     *     0.3040732022736251e+00,     0.2546369261678898e+00,
     *     0.2045251166823099e+00,     0.1538699136085835e+00,
     *     0.1028069379667370e+00,     0.5147184255531770e-01,
     *     0.0e+00                   /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),
     *  wgk(9),wgk(10)/
     *     0.1389013698677008e-02,     0.3890461127099884e-02,
     *     0.6630703915931292e-02,     0.9273279659517763e-02,
     *     0.1182301525349634e-01,     0.1436972950704580e-01,
     *     0.1692088918905327e-01,     0.1941414119394238e-01,
     *     0.2182803582160919e-01,     0.2419116207808060e-01/
      data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),
     *  wgk(17),wgk(18),wgk(19),wgk(20)/
     *     0.2650995488233310e-01,     0.2875404876504129e-01,
     *     0.3090725756238776e-01,     0.3298144705748373e-01,
     *     0.3497933802806002e-01,     0.3688236465182123e-01,
     *     0.3867894562472759e-01,     0.4037453895153596e-01,
     *     0.4196981021516425e-01,     0.4345253970135607e-01/
      data wgk(21),wgk(22),wgk(23),wgk(24),
     *  wgk(25),wgk(26),wgk(27),wgk(28),wgk(29),wgk(30),wgk(31)/
     *     0.4481480013316266e-01,     0.4605923827100699e-01,
     *     0.4718554656929915e-01,     0.4818586175708713e-01,
     *     0.4905543455502978e-01,     0.4979568342707421e-01,
     *     0.5040592140278235e-01,     0.5088179589874961e-01,
     *     0.5122154784925877e-01,     0.5142612853745903e-01,
     *     0.5149472942945157e-01/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/
     *     0.7968192496166606e-02,     0.1846646831109096e-01,
     *     0.2878470788332337e-01,     0.3879919256962705e-01,
     *     0.4840267283059405e-01,     0.5749315621761907e-01,
     *     0.6597422988218050e-01,     0.7375597473770521e-01/
      data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/
     *     0.8075589522942022e-01,     0.8689978720108298e-01,
     *     0.9212252223778613e-01,     0.9636873717464426e-01,
     *     0.9959342058679527e-01,     0.1017623897484055e+00,
     *     0.1028526528935588e+00/
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 30-point gauss rule
c           resk   - result of the 61-point kronrod rule
c           reskh  - approximation to the mean value of f
c                    over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk61
      epmach = r1mach(4)
      uflow = r1mach(1)
c
      centr = 0.5e+00*(b+a)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 61-point kronrod approximation to the
c           integral, and estimate the absolute error.
c
      resg = 0.0e+00
      fc = f(centr)
      resk = wgk(31)*fc
      resabs = abs(resk)
      do 10 j=1,15
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
   10 continue
      do 15 j=1,15
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  15    continue
      reskh = resk*0.5e+00
      resasc = wgk(31)*abs(fc-reskh)
      do 20 j=1,30
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1
     *  ((epmach*0.5e+02)*resabs,abserr)
      return
      end
