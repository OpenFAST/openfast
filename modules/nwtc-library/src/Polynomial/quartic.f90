!+
MODULE PolynomialRoots
! ---------------------------------------------------------------------------
! PURPOSE - Solve for the roots of a polynomial equation with real
!   coefficients, up to quartic order. Retrun a code indicating the nature
!   of the roots found.

! AUTHORS - Alfred H. Morris, Naval Surface Weapons Center, Dahlgren,VA
!           William L. Davis, Naval Surface Weapons Center, Dahlgren,VA
!           Alan Miller,  CSIRO Mathematical & Information Sciences
!                         CLAYTON, VICTORIA, AUSTRALIA 3169
!                         http://www.mel.dms.csiro.au/~alan
!           Ralph L. Carmichael, Public Domain Aeronautical Software
!                         http://www.pdas.com
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!    ??    1.0 AHM&WLD Original coding                     
! 27Feb97  1.1   AM    Converted to be compatible with ELF90
! 12Jul98  1.2   RLC   Module format; numerous style changes
!  4Jan99  1.3   RLC   Made the tests for zero constant term exactly zero


  IMPLICIT NONE

  INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0_4), DP=KIND(1.0_8)
  REAL(DP),PARAMETER,PRIVATE:: ZERO=0.0D0, FOURTH=0.25D0, HALF=0.5D0
  REAL(DP),PARAMETER,PRIVATE:: ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, FOUR=4.0D0
  COMPLEX(DP),PARAMETER,PRIVATE:: CZERO=(0.D0,0.D0)

  REAL(DP),PARAMETER,PRIVATE:: EPS=EPSILON(ONE)

  CHARACTER(LEN=*),PARAMETER,PUBLIC:: POLYROOTS_VERSION= "1.3 (4 Jan 1999)"
  INTEGER,PRIVATE:: outputCode
!    =0 degenerate equation
!    =1 one real root
!    =21 two identical real roots
!    =22 two distinct real roots
!    =23 two complex roots
!    =31 multiple real roots
!    =32 one real and two complex roots
!    =33 three distinct real roots
!    =41
!    =42 two real and two complex roots
!    =43
!    =44 four complex roots

  PRIVATE:: CubeRoot
  PUBLIC:: LinearRoot
  PRIVATE:: OneLargeTwoSmall
  PUBLIC:: QuadraticRoots
  PUBLIC:: CubicRoots
  PUBLIC:: QuarticRoots
  PUBLIC:: SolvePolynomial
!----------------------------------------------------------------------------

  INTERFACE Swap
    MODULE PROCEDURE SwapDouble 
    MODULE PROCEDURE SwapSingle
  END INTERFACE

CONTAINS

!+
FUNCTION CubeRoot(x) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the Cube Root of a REAL(DP) number. If the argument is
!   negative, then the cube root is also negative.

  REAL(DP),INTENT(IN) :: x
  REAL(DP):: f
!----------------------------------------------------------------------------
  IF (x < ZERO) THEN
    f=-EXP(LOG(-x)/THREE)
  ELSE IF (x > ZERO) THEN
    f=EXP(LOG(x)/THREE)
  ELSE
    f=ZERO
  END IF
  RETURN
END Function CubeRoot   ! ---------------------------------------------------

!+
SUBROUTINE LinearRoot(a, z)
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z 
!     AND STORES THE RESULTS IN Z. It is assumed that a(2) is non-zero.
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  REAL(DP),INTENT(OUT):: z
!----------------------------------------------------------------------------
  IF (a(2)==0.0) THEN
    z=ZERO
  ELSE
    z=-a(1)/a(2)
  END IF
  RETURN
END Subroutine LinearRoot   ! -----------------------------------------------

!+
SUBROUTINE OneLargeTwoSmall(a1,a2,a4,w, z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the roots of a cubic when one root, w, is known to be
!   much larger in magnitude than the other two

  REAL(DP),INTENT(IN):: a1,a2,a4
  REAL(DP),INTENT(IN):: w
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


  REAL(DP),DIMENSION(3):: aq
!----------------------------------------------------------------------------
  aq(1)=a1
  aq(2)=a2+a1/w
  aq(3)=-a4*w
  CALL QuadraticRoots(aq, z)
  z(3)=CMPLX(w,ZERO,DP)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3)=z(2)
  z(2)=z(1)
  z(1)=CMPLX(w,ZERO,DP)
  RETURN
END Subroutine OneLargeTwoSmall   ! -----------------------------------------

!+
SUBROUTINE QuadraticRoots(a, z)
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.

  REAL(DP),INTENT(IN),DIMENSION(:):: a
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


  REAL(DP):: d, r, w, x, y
!----------------------------------------------------------------------------
  IF(a(1)==0.0) THEN     ! EPS is a global module constant (private)
    z(1) = CZERO               ! one root is obviously zero
    z(2) = CMPLX(-a(2)/a(3), ZERO,DP)    ! remainder is a linear eq.
    outputCode=21   ! two identical real roots
    RETURN
  END IF

  d = a(2)*a(2) - FOUR*a(1)*a(3)             ! the discriminant
  IF (ABS(d) <= TWO*eps*a(2)*a(2)) THEN
    z(1) = CMPLX(-HALF*a(2)/a(3), ZERO, DP) ! discriminant is tiny
    z(2) = z(1)
    outputCode=22  ! two distinct real roots
    RETURN
  END IF

  r = SQRT(ABS(d))
  IF (d < ZERO) THEN
    x = -HALF*a(2)/a(3)        ! negative discriminant => roots are complex   
    y = ABS(HALF*r/a(3))
    z(1) = CMPLX(x, y, DP)
    z(2) = CMPLX(x,-y, DP)   ! its conjugate
    outputCode=23                        !  COMPLEX ROOTS
    RETURN
  END IF

  IF (a(2) /= ZERO) THEN              ! see Numerical Recipes, sec. 5.5
    w = -(a(2) + SIGN(r,a(2)))
    z(1) = CMPLX(TWO*a(1)/w,  ZERO, DP)
    z(2) = CMPLX(HALF*w/a(3), ZERO, DP)
    outputCode=22           ! two real roots
    RETURN
  END IF

  x = ABS(HALF*r/a(3))   ! a(2)=0 if you get here
  z(1) = CMPLX( x, ZERO, DP)
  z(2) = CMPLX(-x, ZERO, DP)
  outputCode=22
  RETURN
END Subroutine QuadraticRoots   ! -------------------------------------------

!+
SUBROUTINE CubicRoots(a, z)
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z

  REAL(DP),PARAMETER:: RT3=1.7320508075689D0    ! (Sqrt(3)
  REAL (DP) :: aq(3), arg, c, cf, d, p, p1, q, q1
  REAL(DP):: r, ra, rb, rq, rt
  REAL(DP):: r1, s, sf, sq, sum, t, tol, t1, w
  REAL(DP):: w1, w2, x, x1, x2, x3, y, y1, y2, y3

! NOTE -   It is assumed that a(4) is non-zero. No test is made here.
!----------------------------------------------------------------------------
  IF (a(1)==0.0) THEN
    z(1) = CZERO  ! one root is obviously zero
    CALL QuadraticRoots(a(2:4), z(2:3))   ! remaining 2 roots here
    RETURN
  END IF

  p = a(3)/(THREE*a(4))
  q = a(2)/a(4)
  r = a(1)/a(4)
  tol = FOUR*EPS

  c = ZERO
  t = a(2) - p*a(3)
  IF (ABS(t) > tol*ABS(a(2))) c = t/a(4)

  t = TWO*p*p - q
  IF (ABS(t) <= tol*ABS(q)) t = ZERO
  d = r + p*t
  IF (ABS(d) <= tol*ABS(r)) GO TO 110

!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

  s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
  p1 = a(3)/(THREE*s)
  q1 = a(2)/s
  r1 = a(1)/s

  t1 = q - 2.25D0*p*p
  IF (ABS(t1) <= tol*ABS(q)) t1 = ZERO
  w = FOURTH*r1*r1
  w1 = HALF*p1*r1*t
  w2 = q1*q1*t1/27.0D0

  IF (w1 >= ZERO) THEN
    w = w + w1
    sq = w + w2
  ELSE IF (w2 < ZERO) THEN
    sq = w + (w1 + w2)
  ELSE
    w = w + w2
    sq = w + w1
  END IF

  IF (ABS(sq) <= tol*w) sq = ZERO
  rq = ABS(s/a(4))*SQRT(ABS(sq))
  IF (sq >= ZERO) GO TO 40

!                   ALL ROOTS ARE REAL

  arg = ATAN2(rq, -HALF*d)
  cf = COS(arg/THREE)
  sf = SIN(arg/THREE)
  rt = SQRT(-c/THREE)
  y1 = TWO*rt*cf
  y2 = -rt*(cf + rt3*sf)
  y3 = -(d/y1)/y2

  x1 = y1 - p
  x2 = y2 - p
  x3 = y3 - p

  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)
  IF (ABS(x2) > ABS(x3)) CALL Swap(x2,x3)
  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)

  w = x3

  IF (ABS(x2) < 0.1D0*ABS(x3)) GO TO 70
  IF (ABS(x1) < 0.1D0*ABS(x2)) x1 = - (r/x3)/x2
  z(1) = CMPLX(x1, ZERO,DP)
  z(2) = CMPLX(x2, ZERO,DP)
  z(3) = CMPLX(x3, ZERO,DP)
  RETURN

!                  REAL AND COMPLEX ROOTS

40 ra =CubeRoot(-HALF*d - SIGN(rq,d))
  rb = -c/(THREE*ra)
  t = ra + rb
  w = -p
  x = -p
  IF (ABS(t) <= tol*ABS(ra)) GO TO 41
  w = t - p
  x = -HALF*t - p
  IF (ABS(x) <= tol*ABS(p)) x = ZERO
  41 t = ABS(ra - rb)
  y = HALF*rt3*t
  
  IF (t <= tol*ABS(ra)) GO TO 60
  IF (ABS(x) < ABS(y)) GO TO 50
  s = ABS(x)
  t = y/x
  GO TO 51
50 s = ABS(y)
  t = x/y
51 IF (s < 0.1D0*ABS(w)) GO TO 70
  w1 = w/s
  sum = ONE + t*t
  IF (w1*w1 < 0.01D0*sum) w = - ((r/sum)/s)/s
  z(1) = CMPLX(w, ZERO,DP)
  z(2) = CMPLX(x, y,DP)
  z(3) = CMPLX(x,-y,DP)
  RETURN

!               AT LEAST TWO ROOTS ARE EQUAL

60 IF (ABS(x) < ABS(w)) GO TO 61
  IF (ABS(w) < 0.1D0*ABS(x)) w = - (r/x)/x
  z(1) = CMPLX(w, ZERO,DP)
  z(2) = CMPLX(x, ZERO,DP)
  z(3) = z(2)
  RETURN
  61 IF (ABS(x) < 0.1D0*ABS(w)) GO TO 70
  z(1) = CMPLX(x, ZERO,DP)
  z(2) = z(1)
  z(3) = CMPLX(w, ZERO,DP)
  RETURN

!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

70 aq(1) = a(1)
  aq(2) = a(2) + a(1)/w
  aq(3) = -a(4)*w
  CALL QuadraticRoots(aq, z)
  z(3) = CMPLX(w, ZERO,DP)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3) = z(2)
  z(2) = z(1)
  z(1) = CMPLX(w, ZERO,DP)
  RETURN
!-----------------------------------------------------------------------


!                   CASE WHEN D = 0

110 z(1) = CMPLX(-p, ZERO,DP)
  w = SQRT(ABS(c))
  IF (c < ZERO) GO TO 120
  z(2) = CMPLX(-p, w,DP)
  z(3) = CMPLX(-p,-w,DP)
  RETURN

120 IF (p /= ZERO) GO TO 130
  z(2) = CMPLX(w, ZERO,DP)
  z(3) = CMPLX(-w, ZERO,DP)
  RETURN

130 x = -(p + SIGN(w,p))
  z(3) = CMPLX(x, ZERO,DP)
  t = THREE*a(1)/(a(3)*x)
  IF (ABS(p) > ABS(t)) GO TO 131
  z(2) = CMPLX(t, ZERO,DP)
  RETURN
131 z(2) = z(1)
  z(1) = CMPLX(t, ZERO,DP)
  RETURN
END Subroutine CubicRoots   ! -----------------------------------------------


!+
SUBROUTINE QuarticRoots(a,z)
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!               A(1) + A(2)*Z + ... + A(5)*Z**4

  REAL(DP), INTENT(IN)     :: a(:)
  COMPLEX(DP), INTENT(OUT) :: z(:)

  COMPLEX(DP) :: w
  REAL(DP):: b,b2, c, d, e, h, p, q, r, t
  REAL(DP),DIMENSION(4):: temp
  REAL(DP):: u, v, v1, v2, x, x1, x2, x3, y


! NOTE - It is assumed that a(5) is non-zero. No test is made here

!----------------------------------------------------------------------------

  IF (a(1)==0.0) THEN
    z(1) = CZERO    !  one root is obviously zero
    CALL CubicRoots(a(2:), z(2:))
    RETURN
  END IF


  b = a(4)/(FOUR*a(5))
  c = a(3)/a(5)
  d = a(2)/a(5)
  e = a(1)/a(5)
  b2 = b*b

  p = HALF*(c - 6.0D0*b2)
  q = d - TWO*b*(c - FOUR*b2)
  r = b2*(c - THREE*b2) - b*d + e

! SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
! NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
! THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE
!     Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
! WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
! THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.

  temp(1) = -q*q/64.0D0
  temp(2) = 0.25D0*(p*p - r)
  temp(3) =  p
  temp(4) = ONE
  CALL CubicRoots(temp,z)
  IF (AIMAG(z(2)) /= ZERO) GO TO 60

!         THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
!         REORDER THE ROOTS IN INCREASING ORDER

  x1 = REAL(z(1),8)
  x2 = REAL(z(2),8)
  x3 = REAL(z(3),8)
  IF (x1 > x2) CALL Swap(x1,x2)
  IF (x2 > x3) CALL Swap(x2,x3)
  IF (x1 > x2) CALL Swap(x1,x2)

  u = ZERO
  IF (x3 > ZERO) u = SQRT(x3)
  IF (x2 <= ZERO) GO TO 41
  IF (x1 >= ZERO) GO TO 30
  IF (ABS(x1) > x2) GO TO 40
  x1 = ZERO

30 x1 = SQRT(x1)
  x2 = SQRT(x2)
  IF (q > ZERO) x1 = -x1
  temp(1) = (( x1 + x2) + u) - b
  temp(2) = ((-x1 - x2) + u) - b
  temp(3) = (( x1 - x2) - u) - b
  temp(4) = ((-x1 + x2) - u) - b
  CALL SelectSort(temp)
  IF (ABS(temp(1)) >= 0.1D0*ABS(temp(4))) GO TO 31
  t = temp(2)*temp(3)*temp(4)
  IF (t /= ZERO) temp(1) = e/t
31 z(1) = CMPLX(temp(1), ZERO,DP)
  z(2) = CMPLX(temp(2), ZERO,DP)
  z(3) = CMPLX(temp(3), ZERO,DP)
  z(4) = CMPLX(temp(4), ZERO,DP)
  RETURN

40 v1 = SQRT(ABS(x1))
v2 = ZERO
GO TO 50
41 v1 = SQRT(ABS(x1))
v2 = SQRT(ABS(x2))
IF (q < ZERO) u = -u

50 x = -u - b
y = v1 - v2
z(1) = CMPLX(x, y,DP)
z(2) = CMPLX(x,-y,DP)
x =  u - b
y = v1 + v2
z(3) = CMPLX(x, y,DP)
z(4) = CMPLX(x,-y,DP)
RETURN

!                THE RESOLVENT CUBIC HAS COMPLEX ROOTS

60 t = REAL(z(1),8)
x = ZERO
IF (t < ZERO) THEN
  GO TO 61
ELSE IF (t == ZERO) THEN
  GO TO 70
ELSE
  GO TO 62
END IF
61 h = ABS(REAL(z(2),8)) + ABS(AIMAG(z(2)))
IF (ABS(t) <= h) GO TO 70
GO TO 80
62 x = SQRT(t)
IF (q > ZERO) x = -x

70 w = SQRT(z(2))
  u = TWO*REAL(w,8)
  v = TWO*ABS(AIMAG(w))
  t =  x - b
  x1 = t + u
  x2 = t - u
  IF (ABS(x1) <= ABS(x2)) GO TO 71
  t = x1
  x1 = x2
  x2 = t
71 u = -x - b
  h = u*u + v*v
  IF (x1*x1 < 0.01D0*MIN(x2*x2,h)) x1 = e/(x2*h)
  z(1) = CMPLX(x1, ZERO,DP)
  z(2) = CMPLX(x2, ZERO,DP)
  z(3) = CMPLX(u, v,DP)
  z(4) = CMPLX(u,-v,DP)
  RETURN

80 v = SQRT(ABS(t))
  z(1) = CMPLX(-b, v,DP)
  z(2) = CMPLX(-b,-v,DP)
  z(3) = z(1)
  z(4) = z(2)
  RETURN

END Subroutine QuarticRoots

!+
SUBROUTINE SelectSort(a)
! ---------------------------------------------------------------------------
! PURPOSE - Reorder the elements of in increasing order.
  REAL(DP),INTENT(IN OUT),DIMENSION(:):: a

  INTEGER:: j
  INTEGER,DIMENSION(1):: k
! NOTE - This is a n**2 method. It should only be used for small arrays. <25
!----------------------------------------------------------------------------
  DO j=1,SIZE(a)-1
    k=MINLOC(a(j:))
    IF (j /= k(1)) CALL Swap(a(k(1)),a(j))
  END DO
  RETURN
END Subroutine SelectSort   ! -----------------------------------------------

!+
SUBROUTINE SolvePolynomial(quarticCoeff, cubicCoeff, quadraticCoeff, &
  linearCoeff, constantCoeff, code, root1,root2,root3,root4)
! ---------------------------------------------------------------------------
  REAL(DP),INTENT(IN):: quarticCoeff
  REAL(DP),INTENT(IN):: cubicCoeff, quadraticCoeff
  REAL(DP),INTENT(IN):: linearCoeff, constantCoeff
  INTEGER,INTENT(OUT):: code
  COMPLEX(DP),INTENT(OUT):: root1,root2,root3,root4
  REAL(DP),DIMENSION(5):: a
  COMPLEX(DP),DIMENSION(5):: z
!----------------------------------------------------------------------------
  a(1)=constantCoeff
  a(2)=linearCoeff
  a(3)=quadraticCoeff
  a(4)=cubicCoeff
  a(5)=quarticCoeff

  IF (quarticCoeff /= ZERO) THEN
    CALL QuarticRoots(a,z)  
  ELSE IF (cubicCoeff /= ZERO) THEN
    CALL CubicRoots(a,z)
  ELSE IF (quadraticCoeff /= ZERO) THEN
    CALL QuadraticRoots(a,z)
  ELSE IF (linearCoeff /= ZERO) THEN
    z(1)=CMPLX(-constantCoeff/linearCoeff, 0, DP)
    outputCode=1
  ELSE
    outputCode=0    !  { no roots }
  END IF

  code=outputCode
  IF (outputCode > 0) root1=z(1)
  IF (outputCode > 1) root2=z(2)
  IF (outputCode > 23) root3=z(3)
  IF (outputCode > 99) root4=z(4)
  RETURN
END Subroutine SolvePolynomial   ! ------------------------------------------

!+
SUBROUTINE SwapDouble(a,b)
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(DP),INTENT(IN OUT):: a,b
  REAL(DP):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine SwapDouble   ! -----------------------------------------------

!+
SUBROUTINE SwapSingle(a,b)
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(SP),INTENT(IN OUT):: a,b
  REAL(SP):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine SwapSingle   ! -----------------------------------------------


END Module PolynomialRoots   ! ==============================================


