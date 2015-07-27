!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
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
!**********************************************************************************************************************************
MODULE BeamDyn_Subs

   USE BeamDyn_Types
   !USE NWTC_LAPACK

   IMPLICIT NONE


CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateGLL(N, x, w, ErrStat, ErrMsg)
!
! This subroutine determines the (N+1) Gauss-Lobatto-Legendre points x and weights w
!
! For details, see
! @book{Deville-etal:2002,
!  author =    {M. O. Deville and P. F. Fischer and E. H. Mund},
!  title =     {High-Order Methods for Incompressible Fluid Flow},
!  publisher = {Cambridge University Press},
!  address = {Cambridge},
!  year =      2002
!}
!
!..................................................................................................................................

   ! input variables

   INTEGER(IntKi), INTENT(IN   ):: N           ! Order of spectral element
   REAL(ReKi),     INTENT(  OUT):: x(:)      ! location of GLL nodes
   REAL(ReKi),     INTENT(  OUT):: w(:)      ! quadrature weights at GLL nodes
   INTEGER(IntKi), INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),   INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                   :: tol       ! tolerance for newton-raphson solve
   INTEGER(IntKi)               :: maxit     ! maximum allowable iterations in newton-raphson solve
   REAL(ReKi)                   :: x_it      ! current NR-iteration value
   REAL(ReKi)                   :: x_old     ! last NR-iteration value
   REAL(ReKi)                   :: dleg(N+1)   ! legendre polynomial
   INTEGER(IntKi)               :: N1        ! N+1
   INTEGER(IntKi)               :: i         ! do-loop counter
   INTEGER(IntKi)               :: j         ! do-loop counter
   INTEGER(IntKi)               :: k         ! do-loop counter


   ErrStat = ErrID_None
   ErrMsg  = ""

   tol = 1.0D-15
   N1 = N+1
   maxit = 1.0D+03

   ! enter known endpoints  [-1.0, 1.0]
   x(1) = -1.0D+00
   x(N1) = 1.0D+00

   DO i = 1, N1
      x_it = -COS(pi * FLOAT(i-1) / N) ! initial guess - chebyshev points
      DO j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0
         dleg(2) = x_it
         DO k = 2,N
            dleg(k+1) = (  (2.0*REAL(k,DbKi) - 1.0) * dleg(k) * x_it &
                            - (REAL(k,DbKi)-1.0)*dleg(k-1) ) / REAL(k,DbKi)
         ENDDO
         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (REAL(N1,DbKi) * dleg(N1) )
         IF (ABS(x_it - x_old) .lt. tol) THEN
            EXIT
         ENDIF
      ENDDO

      x(i) = x_it
      w(i) = 2.0D0 / (REAL(N * N1, DbKi) * dleg(N1)**2 )

   ENDDO

END SUBROUTINE BD_GenerateGLL
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION BD_Tilde(vect)

   REAL(ReKi),INTENT(IN):: vect(3)
   REAL(ReKi):: BD_Tilde(3,3)

   BD_Tilde = 0.0D0

   BD_Tilde(1,2) = -vect(3)
   BD_Tilde(1,3) = vect(2)
   BD_Tilde(2,1) = vect(3)
   BD_Tilde(2,3) = -vect(1)
   BD_Tilde(3,1) = -vect(2)
   BD_Tilde(3,2) = vect(1)

END FUNCTION BD_Tilde
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvMatrixR(cc,Rr,ErrStat,ErrMsg)
!--------------------------------------------------
! This subroutine computes the rotation tensor (RT)
! given Wiener-Milenkovic rotation parameters
!--------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: cc(:)
   REAL(ReKi),    INTENT(  OUT):: Rr(:,:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: c0
   REAL(ReKi)                  :: c1
   REAL(ReKi)                  :: c2
   REAL(ReKi)                  :: c3
   REAL(ReKi)                  :: tr0

   ErrStat = ErrID_None
   ErrMsg  = ""

   Rr = 0.0D0
   c1 = cc(1)/4.0D0
   c2 = cc(2)/4.0D0
   c3 = cc(3)/4.0D0
   c0 = 0.5D0*(1.0D0-c1*c1-c2*c2-c3*c3)
   tr0 = 1.0D0 - c0
   tr0 = 2.0D0/(tr0*tr0)

   Rr(1,1) = tr0*(c1*c1 + c0*c0) - 1.0D0
   Rr(2,1) = tr0*(c1*c2 + c0*c3)
   Rr(3,1) = tr0*(c1*c3 - c0*c2)
   Rr(1,2) = tr0*(c1*c2 - c0*c3)
   Rr(2,2) = tr0*(c2*c2 + c0*c0) - 1.0D0
   Rr(3,2) = tr0*(c2*c3 + c0*c1)
   Rr(1,3) = tr0*(c1*c3 + c0*c2)
   Rr(2,3) = tr0*(c2*c3 - c0*c1)
   Rr(3,3) = tr0*(c3*c3 + c0*c0) - 1.0D0

END SUBROUTINE BD_CrvMatrixR
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvMatrixH(cc,Hh)

   REAL(ReKi),INTENT(IN)::cc(:)
   REAL(ReKi),INTENT(OUT)::Hh(:,:)

   REAL(ReKi):: cf1,cf2,cf3,cq,ocq,aa,cb0,cb1,cb2,cb3

   cf1 = cc(1)/4.0D0
   cf2 = cc(2)/4.0D0
   cf3 = cc(3)/4.0D0
   cq = cf1 * cf1 + cf2 * cf2 + cf3 * cf3
   ocq = 1.0D0 + cq
   aa = 2.0D0 * ocq * ocq
   cb0 = 2.0D0 * (1.0D0 - cq) / aa
   cb1 = cc(1)/aa
   cb2 = cc(2)/aa
   cb3 = cc(3)/aa

   Hh = 0.0D0

   Hh(1,1) = cb1 * cf1 + cb0
   Hh(2,1) = cb2 * cf1 + cb3
   Hh(3,1) = cb3 * cf1 - cb2
   Hh(1,2) = cb1 * cf2 - cb3
   Hh(2,2) = cb2 * cf2 + cb0
   Hh(3,2) = cb3 * cf2 + cb1
   Hh(1,3) = cb1 * cf3 + cb2
   Hh(2,3) = cb2 * cf3 - cb1
   Hh(3,3) = cb3 * cf3 + cb0

END SUBROUTINE BD_CrvMatrixH
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvCompose( rr, pp, qq, flag, ErrStat, ErrMsg)

!************************************************************************************************************
!   This subroutine composes two Wiener-Milenkovic parameters pp and qq to find the resulting parameter rr
!   This method is detailed in the paper: Bauchau, O.A., 2008, "Interpolation of finite rotations in flexible
!   multi-body dynamics simulations", IMechE, Equation (9).
!   flag = 0: R(rr) = R    (pp) R    (qq)
!   flag = 1: R(rr) = R(T) (pp) R    (qq)
!   flag = 2: R(rr) = R    (pp) R(T) (qq)
!   flag = 3: R(rr) = R(T) (pp) R(T) (qq)
!************************************************************************************************************

   REAL(ReKi),    INTENT(IN   ):: pp(:)     ! Input rotation 1
   REAL(ReKi),    INTENT(IN   ):: qq(:)     ! Input rotation 2
   INTEGER(IntKi),INTENT(IN   ):: flag      ! Option flag
   REAL(ReKi),    INTENT(  OUT):: rr(:)     ! Composed rotation
   INTEGER(IntKi),INTENT(  OUT):: ErrStat   ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg    ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: pp0
   REAL(ReKi)                  :: pp1
   REAL(ReKi)                  :: pp2
   REAL(ReKi)                  :: pp3
   REAL(ReKi)                  :: qq0
   REAL(ReKi)                  :: qq1
   REAL(ReKi)                  :: qq2
   REAL(ReKi)                  :: qq3
   REAL(ReKi)                  :: tr1
   REAL(ReKi)                  :: tr2
   REAL(ReKi)                  :: dd1
   REAL(ReKi)                  :: dd2

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF(flag==1 .OR. flag==3) THEN
       pp1 = -pp(1)
       pp2 = -pp(2)
       pp3 = -pp(3)
   ELSE
       pp1 = pp(1)
       pp2 = pp(2)
       pp3 = pp(3)
   ENDIF
   pp0 = 2.0D0 - (pp1 * pp1 + pp2 * pp2 + pp3 * pp3) / 8.0D0

   IF(flag==2 .OR. flag==3) THEN
       qq1 = -qq(1)
       qq2 = -qq(2)
       qq3 = -qq(3)
   ELSE
       qq1 = qq(1)
       qq2 = qq(2)
       qq3 = qq(3)
   ENDIF
   qq0 = 2.0D0 - (qq1 * qq1 + qq2 * qq2 + qq3 * qq3)/8.0D0

   tr1 = (4.0D0 - pp0) * (4.0D0 - qq0)
   tr2 = pp0 * qq0 - pp1 * qq1 - pp2 * qq2 - pp3 * qq3
   dd1 = tr1 + tr2
   dd2 = tr1 - tr2

   IF(dd1>dd2) THEN
       tr1 = 4.0D0 / dd1
   ELSE
       tr1 = -4.0D0 / dd2
   ENDIF

   rr(1) = tr1 * (pp1 * qq0 + pp0 * qq1 - pp3 * qq2 + pp2 * qq3)
   rr(2) = tr1 * (pp2 * qq0 + pp3 * qq1 + pp0 * qq2 - pp1 * qq3)
   rr(3) = tr1 * (pp3 * qq0 - pp2 * qq1 + pp1 * qq2 + pp0 * qq3)


END SUBROUTINE BD_CrvCompose
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvExtractCrv(Rr,cc,ErrStat,ErrMsg)
!--------------------------------------------------
! This subroutine computes the CRV parameters given
! the rotation matrix
!--------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: Rr(:,:)       ! Rotation Matrix
   REAL(ReKi),    INTENT(  OUT):: cc(:)         ! Crv paramteres
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   !Local variables
   REAL(ReKi)                  :: pivot
   REAL(ReKi)                  :: sm0
   REAL(ReKi)                  :: sm1
   REAL(ReKi)                  :: sm2
   REAL(ReKi)                  :: sm3
   REAL(ReKi)                  :: em
   REAL(ReKi)                  :: temp
   INTEGER(IntKi)              :: ipivot

   ErrStat = ErrID_None
   ErrMsg  = ""

   cc = 0.0D0
   ipivot = 0
   pivot = Rr(1,1) + Rr(2,2) + Rr(3,3)

   IF(Rr(1,1) .GT. pivot) THEN
       pivot = Rr(1,1)
       ipivot = 1
   ENDIF
   IF(Rr(2,2) .GT. pivot) THEN
       pivot = Rr(2,2)
       ipivot = 2
   ENDIF
   IF(Rr(3,3) .GT. pivot) THEN
       pivot = Rr(3,3)
       ipivot = 3
   ENDIF

   IF(ipivot .EQ. 0) THEN
       sm0 = 1.0D0 + Rr(1,1) + Rr(2,2) + Rr(3,3)
       sm1 = Rr(3,2) - Rr(2,3)
       sm2 = Rr(1,3) - Rr(3,1)
       sm3 = Rr(2,1) - Rr(1,2)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm0))
       ELSE
           temp = ABS(2.0D0*SQRT(sm0))
       ENDIF
       em = sm0 + temp
   ELSEIF(ipivot .EQ. 1) THEN
       sm0 = Rr(3,2) - Rr(2,3)
       sm1 = 1.0D0 + Rr(1,1) - Rr(2,2) - Rr(3,3)
       sm2 = Rr(1,2) + Rr(2,1)
       sm3 = Rr(1,3) + Rr(3,1)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm1))
       ELSE
           temp = ABS(2.0D0*SQRT(sm1))
       ENDIF
       em = sm0 + temp
   ELSEIF(ipivot .EQ. 2) THEN
       sm0 = Rr(1,3) - Rr(3,1)
       sm1 = Rr(1,2) + Rr(2,1)
       sm2 = 1.0D0 - Rr(1,1) + Rr(2,2) - Rr(3,3)
       sm3 = Rr(2,3) + Rr(3,2)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm2))
       ELSE
           temp = ABS(2.0D0*SQRT(sm2))
       ENDIF
       em = sm0 + temp
   ELSE
       sm0 = Rr(2,1) - Rr(1,2)
       sm1 = Rr(1,3) + Rr(3,1)
       sm2 = Rr(2,3) + Rr(3,2)
       sm3 = 1.0D0 - Rr(1,1) - Rr(2,2) + Rr(3,3)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm3))
       ELSE
           temp = ABS(2.0D0*SQRT(sm3))
       ENDIF
       em = sm0 + temp
   ENDIF

   em = 4.0D0/em
   cc(1) = em*sm1
   cc(2) = em*sm2
   cc(3) = em*sm3

END SUBROUTINE BD_CrvExtractCrv
!-----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_GaussPointWeight(n, x, w, ErrStat, ErrMsg)
   !-------------------------------------------------------------------------------
   ! This subroutine generates n-point gauss-legendre quadrature points and weights
   !-------------------------------------------------------------------------------

   INTEGER(IntKi),INTENT(IN   ):: n       ! Number of Gauss point
   REAL(ReKi),    INTENT(  OUT):: x(:)   ! Gauss point location
   REAL(ReKi),    INTENT(  OUT):: w(:)   ! Gauss point weight
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: x1
   REAL(ReKi)                  :: x2
   REAL(ReKi),        PARAMETER:: eps = 3.d-07
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   REAL(ReKi)                  :: p1
   REAL(ReKi)                  :: p2
   REAL(ReKi)                  :: p3
   REAL(ReKi)                  :: pp
   REAL(ReKi)                  :: xl
   REAL(ReKi)                  :: xm
   REAL(ReKi)                  :: z
   REAL(ReKi)                  :: z1

   ErrStat = ErrID_None
   ErrMsg  = ""

   m=(n+1)/2

   x1 = -1.d0
   x2 = +1.d0

   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   DO 12 i=1,m
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
1      CONTINUE
       p1=1.d0
       p2=0.d0
       DO 11 j=1,n
           p3=p2
           p2=p1
           p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11         CONTINUE
           pp=n*(z*p1-p2)/(z*z-1.d0)
           z1=z
           z=z1-p1/pp
           IF(ABS(z-z1).GT.eps) GOTO 1
           x(i)=xm-xl*z
           x(n+1-i)=xm+xl*z
           w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
           w(n+1-i)=w(i)
12         CONTINUE
   END SUBROUTINE BD_GaussPointWeight
!  (c) copr. 1986-92 numerical recipes software +k$<,(5cl.
!-----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_MotionTensor(RotTen,Pos,MotTen,flag)

   REAL(ReKi),     INTENT(IN   ):: RotTen(:,:)
   REAL(ReKi),     INTENT(IN   ):: Pos(:)
   REAL(ReKi),     INTENT(  OUT):: MotTen(:,:)
   INTEGER(IntKi), INTENT(IN   ):: flag            ! 0: Motion Tensor;
                                                   ! 1: Inverse of Motion Tensor

   MotTen(:,:) = 0.0D0
   IF (flag .EQ. 0) THEN
       MotTen(1:3,1:3) = RotTen(1:3,1:3)
       MotTen(4:6,4:6) = RotTen(1:3,1:3)
       MotTen(1:3,4:6) = MATMUL(BD_Tilde(Pos),RotTen)
   ELSEIF(flag .EQ. 1) THEN
       MotTen(1:3,1:3) = TRANSPOSE(RotTen(1:3,1:3))
       MotTen(4:6,4:6) = TRANSPOSE(RotTen(1:3,1:3))
       MotTen(1:3,4:6) = TRANSPOSE(MATMUL(BD_Tilde(Pos),RotTen))
   ENDIF

   END SUBROUTINE BD_MotionTensor
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE BeamDyn_Subs
