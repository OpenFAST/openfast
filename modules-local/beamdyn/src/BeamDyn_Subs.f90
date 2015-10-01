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
  
   INTEGER, PARAMETER :: BDKi = ReKi
   
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
   REAL(BDKi),     INTENT(  OUT):: x(:)        ! location of GLL nodes
   REAL(BDKi),     INTENT(  OUT):: w(:)        ! quadrature weights at GLL nodes
   INTEGER(IntKi), INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),   INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(BDKi)      , PARAMETER  :: tol   = 10.0_BDKi*EPSILON(tol) / 2.0_BDKi   ! tolerance for newton-raphson solve (ignores 1 significant digit)
   INTEGER(IntKi)  , PARAMETER  :: maxit = 1000                                ! maximum allowable iterations in newton-raphson solve

   REAL(BDKi)                   :: x_it        ! current NR-iteration value
   REAL(BDKi)                   :: x_old       ! last NR-iteration value
   REAL(BDKi)                   :: dleg(N+1)   ! legendre polynomial
   INTEGER(IntKi)               :: N1          ! N+1
   INTEGER(IntKi)               :: i           ! do-loop counter
   INTEGER(IntKi)               :: j           ! do-loop counter
   INTEGER(IntKi)               :: k           ! do-loop counter


   ErrStat = ErrID_None
   ErrMsg  = ""

   N1 = N+1

   ! enter known endpoints  [-1.0, 1.0]
   x(1) = -1.0_BDKi
   x(N1) = 1.0_BDKi

   DO i = 1, N1
      
      x_it = -COS(pi_D * REAL(i-1,BDKi) / N) ! initial guess - chebyshev points
      DO j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0_BDKi
         dleg(2) = x_it
         DO k = 2,N
            dleg(k+1) = (  REAL(2*k - 1,BDKi) * dleg(k) * x_it &
                         - REAL(  k - 1,BDKi) * dleg(k-1) ) / REAL(k,BDKi)
         ENDDO
         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (REAL(N1,BDKi) * dleg(N1) )
         IF (ABS(x_it - x_old) .lt. tol) THEN
            EXIT
         ENDIF
      ENDDO

      x(i) = x_it
      w(i) = 2.0_BDKi / (REAL(N * N1, BDKi) * dleg(N1)**2 )

   ENDDO

END SUBROUTINE BD_GenerateGLL
!-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION BD_Tilde(vect)
! this function returns the skew-symmetric matrix formed by the values of vect

   REAL(BDKi),INTENT(IN):: vect(3)   
   REAL(BDKi)           :: BD_Tilde(3,3)

   BD_Tilde(1,1) =  0.0_BDKi
   BD_Tilde(2,1) =  vect(3)
   BD_Tilde(3,1) = -vect(2)

   BD_Tilde(1,2) = -vect(3)
   BD_Tilde(2,2) =  0.0_BDKi
   BD_Tilde(3,2) =  vect(1)

   BD_Tilde(1,3) =  vect(2)   
   BD_Tilde(2,3) = -vect(1)
   BD_Tilde(3,3) =  0.0_BDKi

END FUNCTION BD_Tilde
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvMatrixR(cc,Rr,ErrStat,ErrMsg)
!--------------------------------------------------
! This subroutine computes the rotation tensor (RT)
! given Wiener-Milenkovic rotation parameters
!--------------------------------------------------
   REAL(BdKi),    INTENT(IN   ):: cc(:)
   REAL(BdKi),    INTENT(  OUT):: Rr(:,:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(BDKi)                  :: c0
   REAL(BDKi)                  :: c1
   REAL(BDKi)                  :: c2
   REAL(BDKi)                  :: c3
   REAL(BDKi)                  :: tr0

   ErrStat = ErrID_None
   ErrMsg  = ""

   c1  = cc(1)/4.0_BDKi
   c2  = cc(2)/4.0_BDKi
   c3  = cc(3)/4.0_BDKi
   c0  = 0.5_BDKi * (1.0_BDKi-c1*c1-c2*c2-c3*c3)
   tr0 = 1.0_BDKi - c0
   tr0 = 2.0_BDKi/(tr0*tr0)

   Rr(1,1) = tr0*(c1*c1 + c0*c0) - 1.0_BDKi
   Rr(2,1) = tr0*(c1*c2 + c0*c3)
   Rr(3,1) = tr0*(c1*c3 - c0*c2)
   
   Rr(1,2) = tr0*(c1*c2 - c0*c3)
   Rr(2,2) = tr0*(c2*c2 + c0*c0) - 1.0_BDKi
   Rr(3,2) = tr0*(c2*c3 + c0*c1)
   
   Rr(1,3) = tr0*(c1*c3 + c0*c2)
   Rr(2,3) = tr0*(c2*c3 - c0*c1)
   Rr(3,3) = tr0*(c3*c3 + c0*c0) - 1.0_BDKi

END SUBROUTINE BD_CrvMatrixR
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvMatrixH(cc,Hh)

   REAL(BDKi),INTENT(IN) ::cc(:)
   REAL(BDKi),INTENT(OUT)::Hh(:,:)

   REAL(BDKi):: cf1,cf2,cf3,cq,ocq,aa,cb0,cb1,cb2,cb3

   cf1 = cc(1)/4.0_BDKi
   cf2 = cc(2)/4.0_BDKi
   cf3 = cc(3)/4.0_BDKi
   cq  = cf1 * cf1 + cf2 * cf2 + cf3 * cf3
   ocq = 1.0_BDKi + cq
   aa  = 2.0_BDKi * ocq * ocq
   cb0 = 2.0_BDKi * (1.0_BDKi - cq) / aa
   cb1 = cc(1)/aa
   cb2 = cc(2)/aa
   cb3 = cc(3)/aa

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

   REAL(BDKi),    INTENT(IN   ):: pp(:)     ! Input rotation 1
   REAL(BDKi),    INTENT(IN   ):: qq(:)     ! Input rotation 2
   INTEGER       ,INTENT(IN   ):: flag      ! Option flag
   REAL(BDKi),    INTENT(  OUT):: rr(:)     ! Composed rotation
   INTEGER(IntKi),INTENT(  OUT):: ErrStat   ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg    ! Error message if ErrStat /= ErrID_None

   REAL(BDKi)                  :: pp0
   REAL(BDKi)                  :: pp1
   REAL(BDKi)                  :: pp2
   REAL(BDKi)                  :: pp3
   REAL(BDKi)                  :: qq0
   REAL(BDKi)                  :: qq1
   REAL(BDKi)                  :: qq2
   REAL(BDKi)                  :: qq3
   REAL(BDKi)                  :: tr1
   REAL(BDKi)                  :: tr2
   REAL(BDKi)                  :: dd1
   REAL(BDKi)                  :: dd2

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
   pp0 = 2.0_BDKi - (pp1 * pp1 + pp2 * pp2 + pp3 * pp3) / 8.0_BDKi

   IF(flag==2 .OR. flag==3) THEN
       qq1 = -qq(1)
       qq2 = -qq(2)
       qq3 = -qq(3)
   ELSE
       qq1 = qq(1)
       qq2 = qq(2)
       qq3 = qq(3)
   ENDIF
   qq0 = 2.0_BDKi - (qq1 * qq1 + qq2 * qq2 + qq3 * qq3)/8.0_BDKi

   tr1 = (4.0_BDKi - pp0) * (4.0_BDKi - qq0)
   tr2 = pp0 * qq0 - pp1 * qq1 - pp2 * qq2 - pp3 * qq3
   dd1 = tr1 + tr2
   dd2 = tr1 - tr2

   IF(dd1>dd2) THEN
       tr1 = 4.0_BDKi / dd1
   ELSE
       tr1 = -4.0_BDKi / dd2
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

   REAL(BDKi),    INTENT(IN   ):: Rr(:,:)       ! Rotation Matrix
   REAL(BDKi),    INTENT(  OUT):: cc(:)         ! Crv paramteres
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   !Local variables
   REAL(BDKi)                  :: pivot
   REAL(BDKi)                  :: sm0
   REAL(BDKi)                  :: sm1
   REAL(BDKi)                  :: sm2
   REAL(BDKi)                  :: sm3
   REAL(BDKi)                  :: em
   REAL(BDKi)                  :: temp

   ErrStat = ErrID_None
   ErrMsg  = ""

   pivot = Rr(1,1) + Rr(2,2) + Rr(3,3)

   IF (Rr(3,3) .GT. pivot) THEN
      !pivot = Rr(3,3)
      !ipivot = 3
      sm0  = Rr(2,1) - Rr(1,2)
      sm1  = Rr(1,3) + Rr(3,1)
      sm2  = Rr(2,3) + Rr(3,2)
      sm3  = 1.0_BDKi - Rr(1,1) - Rr(2,2) + Rr(3,3)
      temp = SIGN( 2.0_BDKi*SQRT(ABS(sm3)), sm0 ) 
   ELSEIF (Rr(2,2) .GT. pivot) THEN
      !pivot = Rr(2,2)
      !ipivot = 2
      sm0  = Rr(1,3) - Rr(3,1)
      sm1  = Rr(1,2) + Rr(2,1)
      sm2  = 1.0_BDKi - Rr(1,1) + Rr(2,2) - Rr(3,3)
      sm3  = Rr(2,3) + Rr(3,2)
      temp = SIGN( 2.0_BDKi*SQRT(ABS(sm2)), sm0 )       
   ELSEIF (Rr(1,1) .GT. pivot) THEN
      !pivot = Rr(1,1)
      !ipivot = 1
      sm0  = Rr(3,2) - Rr(2,3)
      sm1  = 1.0_BDKi + Rr(1,1) - Rr(2,2) - Rr(3,3)
      sm2  = Rr(1,2) + Rr(2,1)
      sm3  = Rr(1,3) + Rr(3,1)
      temp = SIGN( 2.0_BDKi*SQRT(ABS(sm1)), sm0 )
   ELSE
      !ipivot = 0
      sm0  = 1.0_BDKi + Rr(1,1) + Rr(2,2) + Rr(3,3)
      sm1  = Rr(3,2) - Rr(2,3)
      sm2  = Rr(1,3) - Rr(3,1)
      sm3  = Rr(2,1) - Rr(1,2)
      temp = SIGN( 2.0_BDKi*SQRT(ABS(sm0)), sm0 )       
   ENDIF

   em = sm0 + temp
   em = 4.0_BDKi/em
   cc(1) = em*sm1
   cc(2) = em*sm2
   cc(3) = em*sm3

END SUBROUTINE BD_CrvExtractCrv
!------------------------------------------------------------------------------
SUBROUTINE BD_GaussPointWeight(n, x, w, ErrStat, ErrMsg)
!---------------------------------------------------------------------------
! This subroutine generates n-point gauss-legendre quadrature points and weights
! 
! Subroutine is based on well-known formulas for generating Legendre polynomials, the 
! roots of which are the Gauss-Legendre quadrature points.  Also used are well-known
! expressions for the quadrature weights associated with the GL quadrature points.  
! The basic idea of the logic is to use the roots of the Chebyshev polynomial as
! an initial guess for the roots of the Legendre polynomial, and to then use Newton
! iteration to find the "exact" roots.
!-------------------------------------------------------------------------

   INTEGER(IntKi),INTENT(IN   ):: n       ! Number of Gauss point
   REAL(BDKi),    INTENT(  OUT):: x(:)    ! Gauss point location
   REAL(BDKi),    INTENT(  OUT):: w(:)    ! Gauss point weight
   INTEGER(IntKi),INTENT(  OUT):: ErrStat ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg  ! Error message if ErrStat /=

   ! local variables

   REAL(BDKi),  PARAMETER:: eps = 1000.0_BDKi*EPSILON(eps) / 2.0_BDKi   ! tolerance for newton-raphson solve (ignores 3 significant digits)

   
   INTEGER(IntKi) :: n_half 

   REAL(BDKi)     :: n_real  ! real version of loop index
   INTEGER(IntKi) :: i       ! loop index 
   INTEGER(IntKi) :: j       ! loop index
   REAL(BDKi)     :: j_real  ! real version of loop index

   INTEGER(IntKi) :: newton     ! newton-iteration index
   INTEGER(IntKi), parameter :: newton_max = 10 ! maximum number of newton iterations; should not need many!

   REAL(BDKi)     :: p1  ! legendre polynomial j-2
   REAL(BDKi)     :: p2  ! legendre polynomial j-1
   REAL(BDKi)     :: p3  ! legendre polynomial j

   REAL(BDKi)     :: dp3  ! derivative of legendre polynomial j

   if (n .lt. 2) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'BD_GaussPointWeight: invalid value of n.'
   else
      ErrStat = ErrID_None
      ErrMsg  = ''      
   endif

   n_real = REAL(n,BDKi)

   n_half = (n+1) / 2  ! integer for "half" the points: e.g., n=5 => n_half=3, n=4 => n_half = 2

   ! loop through each of the "half" n points
   do i = 1, n_half

      ! Intial guess for ith root; based on ith root of Chebyshev polynomial
      x(i) = - cos( REAL( 2*i - 1, BDKi)  * pi_D / (2.0_BDKi * n_real) )

      ! initialize newton iteration counter
      newton = 0

      p3 = 1. ! some intial big number;  logic will always take at least one newton iteration
      do while (abs(p3) .gt. eps .and. newton .le. newton_max) 

         newton = newton + 1  ! keep track of number of newton iterations

         ! build the legendre polynomial and its derivative evaluated at current root guess
         p1 = 1.0_BDKi  ! zeroth order legendre polynomial
         p2 = x(i)      ! first-order legendre polynomial
 
         do j = 2, n

            j_real = REAL(j, BDKi)

            ! recurrence relation for legendre polynomial
            p3 = ( REAL(2*j - 1, BDKi) * x(i) * p2 - REAL(j - 1, BDKi) * p1 ) / j_real

            ! derivative of legendre polynomial; needed for newton iteration
            dp3 = j_real * ( - x(i) * p3 + p2 ) / ( 1.0_BDKi - x(i)*x(i) )
    
            ! save values for next iteration 
            p1 = p2
            p2 = p3

         enddo
  
         ! newton iteration: xnew = xold - P_n(xold) / P_n'(xold)
         x(i) = x(i) - p3 / dp3

      enddo
     
      if (abs(x(i)) .lt. eps) x(i) = 0.0_BDKi

      ! fill in the other "side" of the array 
      x(n-i+1) = - x(i)

      ! calculate weight for this Gauss point location
      w(i) = 2.0_BDKi / ( ( 1.0_BDKi - x(i)*x(i) ) * dp3 * dp3 )
      w(n-i+1) = w(i)

   enddo

   END SUBROUTINE BD_GaussPointWeight
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_MotionTensor(RotTen,Pos,MotTen,flag)

   REAL(BDKi),     INTENT(IN   ):: RotTen(:,:)
   REAL(BDKi),     INTENT(IN   ):: Pos(:)
   REAL(BDKi),     INTENT(  OUT):: MotTen(:,:)
   INTEGER(IntKi), INTENT(IN   ):: flag            ! 0: Motion Tensor;
                                                   ! 1: Inverse of Motion Tensor

   MotTen = 0.0_BDKi
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
