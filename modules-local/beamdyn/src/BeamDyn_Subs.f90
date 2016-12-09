!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
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
  
   INTEGER, PARAMETER :: BDKi = R8Ki
   INTEGER, PARAMETER :: FLAG_R1R2   = 0 !<   BD_CrvCompose flag = 0: R(rr) = R    (pp) R    (qq)
   INTEGER, PARAMETER :: FLAG_R1TR2  = 1 !<   BD_CrvCompose flag = 1: R(rr) = R(T) (pp) R    (qq)
   INTEGER, PARAMETER :: FLAG_R1R2T  = 2 !<   BD_CrvCompose flag = 2: R(rr) = R    (pp) R(T) (qq)
   INTEGER, PARAMETER :: FLAG_R1TR2T = 3 !<   BD_CrvCompose flag = 3: R(rr) = R(T) (pp) R(T) (qq)

   INTEGER, PARAMETER :: GAUSS_QUADRATURE = 1 !<   p%quadrature method: gaussian quadrature
   INTEGER, PARAMETER :: TRAP_QUADRATURE  = 2 !<   p%quadrature method: trapeziodal quadrature
         
CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine determines the (N+1) Gauss-Lobatto-Legendre points x and weights w
!!
!! For details, see
!! \@book{Deville-etal:2002,
!!  author =    {M. O. Deville and P. F. Fischer and E. H. Mund},
!!  title =     {High-Order Methods for Incompressible Fluid Flow},
!!  publisher = {Cambridge University Press},
!!  address = {Cambridge},
!!  year =      2002
!!}
SUBROUTINE BD_GenerateGLL(N1, GLL_nodes, ErrStat, ErrMsg)

   ! input variables

   INTEGER(IntKi),          INTENT(IN   ) :: N1             !< 1 + the order of spectral element, or equivalently p%node_elem
   REAL(BDKi), ALLOCATABLE, INTENT(  OUT) :: GLL_nodes(:)   !< location of GLL nodes
   INTEGER(IntKi),          INTENT(  OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   REAL(BDKi)      , PARAMETER            :: tol   = 10.0_BDKi*EPSILON(tol) / 2.0_BDKi   ! tolerance for newton-raphson solve (ignores 1 significant digit)
   INTEGER(IntKi)  , PARAMETER            :: maxit = 1000                                ! maximum allowable iterations in newton-raphson solve
                                          
   REAL(BDKi)                             :: x_it        ! current NR-iteration value
   REAL(BDKi)                             :: x_old       ! last NR-iteration value
   REAL(BDKi)                             :: dleg(N1)    ! legendre polynomial
   INTEGER(IntKi)                         :: N           ! Order of spectral element 
   INTEGER(IntKi)                         :: i           ! do-loop counter
   INTEGER(IntKi)                         :: j           ! do-loop counter
   INTEGER(IntKi)                         :: k           ! do-loop counter

   INTEGER(IntKi)                         :: ErrStat2
   CHARACTER(ErrMsgLen)                   :: ErrMsg2
   character(*), parameter                :: RoutineName = 'BD_GenerateGLL'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   N = N1 - 1

   
   CALL AllocAry(GLL_nodes,N1,'GLL points array',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
  !CALL AllocAry(GLL_weights,N1,'GLL weight array',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return
   
   
   !bjj: these are getting overwritten later....
   ! enter known endpoints  [-1.0, 1.0] 
   GLL_nodes(1) = -1.0_BDKi
   GLL_nodes(N1) = 1.0_BDKi

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
         IF (ABS(x_it - x_old) .lt. tol) EXIT
      ENDDO

      GLL_nodes(i) = x_it
     !GLL_weights(i) = 2.0_BDKi / (REAL(N * N1, BDKi) * dleg(N1)**2 )

   ENDDO

END SUBROUTINE BD_GenerateGLL
!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the rotation tensor (RT=\f$ \phi \f$ )
!! given Wiener-Milenkovic rotation parameters, \f$ \vec{c} \f$, 
!! where \f$ \vec{c}=4\tan\left(\frac{\phi}{4}\right)\bar{n} \f$.
SUBROUTINE BD_CrvMatrixR(cc,Rr)

   REAL(BdKi),    INTENT(IN   ):: cc(3)         !< \f$ \vec{c}\f$
   REAL(BdKi),    INTENT(  OUT):: Rr(3,3)       !< returned rotation tensor (transpose of DCM matrix)

   REAL(BDKi)                  :: c0
   REAL(BDKi)                  :: c1
   REAL(BDKi)                  :: c2
   REAL(BDKi)                  :: c3
   REAL(BDKi)                  :: tr0


   !bjj: these are 1/4 the values described in the BeamDyn manual.
   c1  = cc(1)/4.0_BDKi
   c2  = cc(2)/4.0_BDKi
   c3  = cc(3)/4.0_BDKi
   c0  = 0.5_BDKi * (1.0_BDKi-c1*c1-c2*c2-c3*c3)
   
   tr0 = 1.0_BDKi - c0
   tr0 = 2.0_BDKi/(tr0*tr0)

   !bjj: this is the equation from the BD manual, however it doesn't seem to equate to what is calculated here...
   !> \f{equation}{
   !!  \tens{R} (\vec{c}) = \frac{1}{(4-c_0)^2}
   !! \begin{bmatrix}
   !! c_0^2 + c_1^2 - c_2^2 - c_3^2 & 2(c_1 c_2 - c_0 c_3) & 2(c_1 c_3 + c_0 c_2) \\
   !! 2(c_1 c_2 + c_0 c_3) & c_0^2 - c_1^2 + c_2^2 - c_3^2 & 2(c_2 c_3 - c_0 c_1) \\
   !! 2(c_1 c_3 - c_0 c_2)  & 2(c_2 c_3 + c_0 c_1) & c_0^2 - c_1^2 - c_2^2 + c_3^2 \\
   !! \end{bmatrix}
   !! }\f

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

   REAL(BDKi),INTENT(IN) ::cc(3)
   REAL(BDKi),INTENT(OUT)::Hh(3,3)

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
!>   This subroutine composes two Wiener-Milenkovic parameters pp and qq to find the resulting parameter rr
!!   This method is detailed in the paper: Bauchau, O.A., 2008, "Interpolation of finite rotations in flexible
!!   multi-body dynamics simulations", IMechE, Equation (9). 
!!   see http://www.dymoresolutions.com/publications/BauchauEppleHeo08.pdf \n
!!   FLAG_R1R2:   flag = 0: R(rr) = R    (pp) R    (qq) \n
!!   FLAG_R1TR2:  flag = 1: R(rr) = R(T) (pp) R    (qq) \n
!!   FLAG_R1R2T:  flag = 2: R(rr) = R    (pp) R(T) (qq) \n
!!   FLAG_R1TR2T: flag = 3: R(rr) = R(T) (pp) R(T) (qq)
SUBROUTINE BD_CrvCompose( rr, pp, qq, flag)

   REAL(BDKi),    INTENT(IN   ):: pp(3)     !< Input rotation 1
   REAL(BDKi),    INTENT(IN   ):: qq(3)     !< Input rotation 2
   INTEGER       ,INTENT(IN   ):: flag      !< Option flag
   REAL(BDKi),    INTENT(  OUT):: rr(3)     !< Composed rotation

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


   IF(flag==FLAG_R1TR2 .OR. flag==FLAG_R1TR2T) THEN ! "transpose" (negative) of first rotation parameter
       pp1 = -pp(1)
       pp2 = -pp(2)
       pp3 = -pp(3)
   ELSE
       pp1 = pp(1)
       pp2 = pp(2)
       pp3 = pp(3)
   ENDIF
   pp0 = 2.0_BDKi - (pp1 * pp1 + pp2 * pp2 + pp3 * pp3) / 8.0_BDKi

   IF(flag==FLAG_R1R2T .OR. flag==FLAG_R1TR2T) THEN ! "transpose" (negative) of second rotation parameter
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
!> This subroutine computes the CRV parameters given
!! the rotation matrix
SUBROUTINE BD_CrvExtractCrv(Rr,cc)

   REAL(BDKi),    INTENT(IN   ):: Rr(3,3)       !< Rotation Matrix
   REAL(BDKi),    INTENT(  OUT):: cc(3)         !< Crv paramteres

   !Local variables
   REAL(BDKi)                  :: pivot
   REAL(BDKi)                  :: sm0
   REAL(BDKi)                  :: sm1
   REAL(BDKi)                  :: sm2
   REAL(BDKi)                  :: sm3
   REAL(BDKi)                  :: em
   REAL(BDKi)                  :: temp


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
!> This subroutine generates n-point gauss-legendre quadrature points and weights
!! 
!! Subroutine is based on well-known formulas for generating Legendre polynomials, the 
!! roots of which are the Gauss-Legendre quadrature points.  Also used are well-known
!! expressions for the quadrature weights associated with the GL quadrature points.  
!! The basic idea of the logic is to use the roots of the Chebyshev polynomial as
!! an initial guess for the roots of the Legendre polynomial, and to then use Newton
!! iteration to find the "exact" roots.
SUBROUTINE BD_GaussPointWeight(n, x, w, ErrStat, ErrMsg)

   INTEGER(IntKi),INTENT(IN   ):: n       !< Number of Gauss point
   REAL(BDKi),    INTENT(  OUT):: x(n)    !< Gauss point location
   REAL(BDKi),    INTENT(  OUT):: w(n)    !< Gauss point weight
   INTEGER(IntKi),INTENT(  OUT):: ErrStat !< Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg  !< Error message if ErrStat /=

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
! This subroutine computes trapezoidal quadrature points and weights, p%GL and p%GLw
SUBROUTINE BD_TrapezoidalPointWeight(p, InputFileData)

   TYPE(BD_ParameterType),INTENT(INOUT):: p              !< BeamDyn parameters
   TYPE(BD_InputFile),    INTENT(IN   ):: InputFileData  !< BeamDyn input-file data
   
   ! local variables
   REAL(BDKi)                 :: denom ! denominator for quadrature weight computations
   
   INTEGER(IntKi)             :: indx, temp_id0, temp_id1  ! temporary index into GL array
   INTEGER(IntKi), parameter  :: id0 = 1
   INTEGER(IntKi)             :: id1, j
   
  
      ! compute the trapezoidal quadrature points, p%GL and scale to range [-1,1]:
      !  If there is refinement, this will add new points between the specified ones. If p%refine == 1, can skip this.
   p%GL(1) = InputFileData%InpBl%station_eta(id0)
   DO j = 2,p%ngp
      indx =  id0+(j-2_IntKi)/p%refine       ! forced rounding to integer???  --> (J-2)/p%refine may not be integer.
      p%GL(j) =  InputFileData%InpBl%station_eta(indx) + &
               ((InputFileData%InpBl%station_eta(indx+1) - InputFileData%InpBl%station_eta(indx))/p%refine) * (MOD(j-2,p%refine) + 1)
   ENDDO   
   p%GL = 2.0_BDKi*p%GL - 1.0_BDKi     ! rescale to range [-1,1]


      ! compute the trapezoidal quadrature weights, p%GLw:
   id1 = InputFileData%kp_member(1)             !adp: Why is this only checking the first member (size(kp_member) will be number of GL points)????
   temp_id0 = (id0 - 1)*p%refine + 1            ! Starting index in GL --> always going to be 1
   temp_id1 = (id1 - 1)*p%refine + 1            ! ending index in GL --> will be  size(p%GL)
   denom = p%GL(temp_id1) - p%GL(temp_id0)      ! This is the range of GL --> for single member, is always == 2

   p%GLw(1)     =  (p%GL(temp_id0+1) - p%GL(temp_id0    ))/denom   
   DO j=2,p%ngp-1      
      p%GLw(j)  =  (p%GL(temp_id0+j) - p%GL(temp_id0+j-2))/denom
   ENDDO
   p%GLw(p%ngp) =  (p%GL(temp_id1  ) - p%GL(temp_id1-1  ))/denom
   
   
      
END SUBROUTINE BD_TrapezoidalPointWeight
!-----------------------------------------------------------------------------------------------------------------------------------
!adp:  This routine appears to be useless.
!SUBROUTINE BD_MotionTensor(RotTen,Pos,MotTen,flag)
!
!   REAL(BDKi),     INTENT(IN   ):: RotTen(3,3)
!   REAL(BDKi),     INTENT(IN   ):: Pos(3)
!   REAL(BDKi),     INTENT(  OUT):: MotTen(6,6)
!   INTEGER(IntKi), INTENT(IN   ):: flag            !< 0: Motion Tensor;
!                                                   !! 1: Inverse of Motion Tensor
!
!   MotTen = 0.0_BDKi
!   IF (flag .EQ. 0) THEN
!       MotTen(1:3,1:3) = RotTen(1:3,1:3)
!       MotTen(4:6,4:6) = RotTen(1:3,1:3)
!       MotTen(1:3,4:6) = MATMUL(SkewSymMat(Pos),RotTen)
!   ELSEIF(flag .EQ. 1) THEN
!       MotTen(1:3,1:3) = TRANSPOSE(RotTen(1:3,1:3))
!       MotTen(4:6,4:6) = TRANSPOSE(RotTen(1:3,1:3))
!       MotTen(1:3,4:6) = TRANSPOSE(MATMUL(SkewSymMat(Pos),RotTen))
!   ENDIF
!
!END SUBROUTINE BD_MotionTensor
!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates y%BldMotion%TranslationDisp, y%BldMotion%Orientation, y%BldMotion%TranslationVel, and
!! y%BldMotion%RotationVel, which depend only on states and parameters.
SUBROUTINE Set_BldMotion_NoAcc(p, x, y)

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2
   REAL(BDKi)                                   :: cc(3)
   REAL(BDKi)                                   :: cc0(3)
   REAL(BDKi)                                   :: temp_cc(3)
   REAL(BDKi)                                   :: temp_R(3,3)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Set_BldMotion_NoAcc'
   

   DO i=1,p%elem_total
      DO j=1,p%node_elem
         temp_id = (i-1)*(p%node_elem-1)+j       ! The last node of the first element is used as the first node in the second element.
         temp_id2= (i-1)*p%node_elem+j
         
         temp_cc = MATMUL(p%GlbRot,x%q(1:3,temp_id))
         y%BldMotion%TranslationDisp(1,temp_id2) = temp_cc(2)
         y%BldMotion%TranslationDisp(2,temp_id2) = temp_cc(3)
         y%BldMotion%TranslationDisp(3,temp_id2) = temp_cc(1)
         
         cc = MATMUL(p%GlbRot,x%q(4:6,temp_id))
         cc0 = p%uuN0(4:6,j,i)
         cc0 = MATMUL(p%GlbRot,cc0)
         CALL BD_CrvCompose(temp_cc,p%Glb_crv,cc0,FLAG_R1R2) ! temp_cc = p%Glb_crv composed with cc0
         CALL BD_CrvCompose(cc0,cc,temp_cc,FLAG_R1R2) ! cc0 = cc composed with temp_cc
         temp_cc(1) = cc0(2)
         temp_cc(2) = cc0(3)
         temp_cc(3) = cc0(1)
         CALL BD_CrvMatrixR(temp_cc,temp_R)  ! returns temp_R (the transpose of the DCM orientation matrix)
         y%BldMotion%Orientation(1:3,1:3,temp_id2) = TRANSPOSE(temp_R(1:3,1:3))

         temp_id = (i-1)*(p%node_elem-1)+j
         temp_cc = MATMUL(p%GlbRot,x%dqdt(1:3,temp_id))
         y%BldMotion%TranslationVel(1,temp_id2) = temp_cc(2)
         y%BldMotion%TranslationVel(2,temp_id2) = temp_cc(3)
         y%BldMotion%TranslationVel(3,temp_id2) = temp_cc(1)
         
         temp_cc(:) = MATMUL(p%GlbRot,x%dqdt(4:6,temp_id))
         y%BldMotion%RotationVel(1,temp_id2) = temp_cc(2)
         y%BldMotion%RotationVel(2,temp_id2) = temp_cc(3)
         y%BldMotion%RotationVel(3,temp_id2) = temp_cc(1)
      ENDDO
   ENDDO
   
END SUBROUTINE Set_BldMotion_NoAcc
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE BeamDyn_Subs
