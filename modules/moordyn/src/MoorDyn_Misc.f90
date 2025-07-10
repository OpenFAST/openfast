!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020-2021 Alliance for Sustainable Energy, LLC
! Copyright (C) 2015-2019 Matthew Hall
!
!    This file is part of MoorDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http:!www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
MODULE MoorDyn_Misc

   USE MoorDyn_Types
   USE SeaSt_WaveField
   USE Current
   USE NWTC_Library
   USE NWTC_FFTPACK
   
   IMPLICIT NONE

   PRIVATE

   INTEGER(IntKi), PARAMETER            :: wordy = 0   ! verbosity level. >1 = more console output

   PUBLIC :: UnitVector
   PUBLIC :: ScaleVector
   PUBLIC :: GetCurvature
   PUBLIC :: GetOrientationAngles
   PUBLIC :: TransformKinematics
   PUBLIC :: TransformKinematicsA
   PUBLIC :: TransformKinematicsAtoB
   PUBLIC :: TranslateForce3to6DOF
   PUBLIC :: TranslateMass3to6DOF
   PUBLIC :: TranslateMass6to6DOF
   PUBLIC :: GetH
   PUBLIC :: RotateM6
   PUBLIC :: RotateM3
   PUBLIC :: CalcOrientation
   PUBLIC :: Inverse3by3
   PUBLIC :: LUsolve
   
   PUBLIC :: getInterpNums
   PUBLIC :: calculate4Dinterpolation
   PUBLIC :: calculate3Dinterpolation
   PUBLIC :: calculate2Dinterpolation
   
   PUBLIC :: getDepthFromBathymetry
   
   PUBLIC :: getWaterKin
   PUBLIC :: setupWaterKin

CONTAINS


   ! ::::::::::::::::::::::::::::::::: math convenience functions ::::::::::::::::::::::::::::::::::
   ! should add error checking if I keep these, but hopefully there are existing NWTCLib functions to replace them

   ! return unit vector (u) and in direction from r1 to r2 and distance between points
   !-----------------------------------------------------------------------
   SUBROUTINE UnitVector( r1, r2, u, Length )          ! note: order of parameters chagned in this function
      
      REAL(DbKi),       INTENT(IN   )  :: r1(:)
      REAL(DbKi),       INTENT(IN   )  :: r2(:)
      REAL(DbKi),       INTENT(  OUT)  :: u(:)
      REAL(DbKi),       INTENT(  OUT)  :: length

      u = r2 - r1
      length = TwoNorm(u)

      IF ( .NOT. EqualRealNos(length, 0.0_DbKi ) ) THEN
        u = u / Length
      END IF

   END SUBROUTINE UnitVector
   !-----------------------------------------------------------------------

   ! scale vector to desired length 
   !-----------------------------------------------------------------------
   SUBROUTINE ScaleVector( u_in, newlength, u_out )
      REAL(DbKi),       INTENT(IN   )  :: u_in(3)     ! input vector
      REAL(DbKi),       INTENT(IN   )  :: newlength   ! desired length of output vector
      REAL(DbKi),       INTENT(INOUT)  :: u_out(3)    ! output vector (hopefully can be the same as u_in without issue)

      REAL(DbKi)                       :: length_squared
      REAL(DbKi)                       :: scaler
      INTEGER(IntKi)                   :: J
      
      length_squared = 0.0;
      DO J=1,3
         length_squared = length_squared + u_in(J)*u_in(J)
      END DO
      
      IF (length_squared > 0) THEN
         scaler = newlength/sqrt(length_squared)
      ELSE                   ! if original vector is zero, return zero
         scaler = 0.0_DbKi
      END IF
      
      DO J=1,3
         u_out(J) = u_in(J) * scaler 
      END DO

   END SUBROUTINE ScaleVector
   !-----------------------------------------------------------------------


   ! convenience function to calculate curvature based on adjacent segments' direction vectors and their combined length
   function GetCurvature(length, q1, q2)
   
      Real(DbKi),   intent(in   ) :: length
      Real(DbKi),   intent(in   ) :: q1(3)
      Real(DbKi),   intent(in   ) :: q2(3)
      Real(DbKi)                  :: GetCurvature
      
      
      Real(DbKi)                  :: q1_dot_q2
      
      ! note "length" here is combined from both segments
      
      q1_dot_q2 = dot_product( q1, q2 )
      
      IF (q1_dot_q2 > 1.0) THEN           ! this is just a small numerical error, so set q1_dot_q2 to 1
         GetCurvature = 0.0_DbKi          ! this occurs when there's no curvature, so return zero curvature
         
      !ELSE IF (q1_dot_q2 < 0)   ! this is a bend of more than 90 degrees, too much, call an error!

      ELSE                                                        ! normal case
         GetCurvature = 4.0/length * sqrt(0.5*(1.0 - q1_dot_q2))    ! this is the normal curvature calculation
      END IF
      
      RETURN
   END function GetCurvature
   

   ! calculate orientation angles of a direction vector
   !-----------------------------------------------------------------------
   subroutine GetOrientationAngles(vec, phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat)
      Real(DbKi),   intent(in   ) :: vec(3) !p1(3),p2(3)
      Real(DbKi),   intent(  out) :: phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat(3)
            
      Real(DbKi)                  ::  vecLen, vecLen2D
 
      vecLen   = SQRT(Dot_Product(vec,vec))
      vecLen2D = SQRT(vec(1)**2+vec(2)**2)
      IF ( vecLen < 0.000001 ) THEN
         IF (wordy > 0) THEN
            print *, "ERROR in GetOrientationAngles in MoorDyn. Supplied vector is near zero" 
            print *, vec
         ENDIF
         k_hat = NaN ! 1.0/0.0
      ELSE
         k_hat = vec / vecLen 
         phi   = atan2(vecLen2D, vec(3))  ! incline angle   
      END IF
      IF ( phi < 0.000001) THEN
         beta = 0.0_ReKi
      ELSE
         beta = atan2(vec(2), vec(1))                    ! heading of incline     
      ENDIF
      sinPhi  = sin(phi)
      cosPhi  = cos(phi)  
      tanPhi  = tan(phi)     
      sinBeta = sin(beta)
      cosBeta = cos(beta)
            
   END subroutine GetOrientationAngles
   !-----------------------------------------------------------------------
   

   ! calculate position and velocity of point based on its position relative to moving 6DOF body
   !-----------------------------------------------------------------------
   SUBROUTINE TransformKinematics(rRelBody, r_in, TransMat, rd_in, r_out, rd_out)
      REAL(DbKi),       INTENT(IN   )  :: rRelBody(:)  ! coordinate of end A
      REAL(DbKi),       INTENT(IN   )  :: r_in(3)      ! Rod unit vector
      REAL(DbKi),       INTENT(IN   )  :: TransMat(3,3)! 
      REAL(DbKi),       INTENT(IN   )  :: rd_in(6)     ! 6DOF velecity vector about Rod end A, in global orientation frame
      REAL(DbKi),       INTENT(  OUT)  :: r_out(3)     ! coordinates of end B
      REAL(DbKi),       INTENT(  OUT)  :: rd_out(3)    ! velocity of end B

      REAL(DbKi)                       :: rRel(3)

      ! rd_in should be in global orientation frame
      ! note: it's okay if r_out and rd_out are 6-size. Only the first 3 will be written, and 4-6 will
      !       already be correct or can be assigned separately from r_in and rd_in (assuming orientation frames are identical)


      ! locations (unrotated reference frame) about platform reference point  (2021-01-05: just transposed TransMat, it was incorrect before)
      rRel(1) = TransMat(1,1)*rRelBody(1) + TransMat(1,2)*rRelBody(2) + TransMat(1,3)*rRelBody(3) ! x
      rRel(2) = TransMat(2,1)*rRelBody(1) + TransMat(2,2)*rRelBody(2) + TransMat(2,3)*rRelBody(3) ! y
      rRel(3) = TransMat(3,1)*rRelBody(1) + TransMat(3,2)*rRelBody(2) + TransMat(3,3)*rRelBody(3) ! z

      ! absolute locations
      r_out = rRel + r_in

      ! absolute velocities
      rd_out(1) =                   - rd_in(6)*rRel(2) + rd_in(5)*rRel(3) + rd_in(1) ! x   
      rd_out(2) =  rd_in(6)*rRel(1)                    - rd_in(4)*rRel(3) + rd_in(2) ! y
      rd_out(3) = -rd_in(5)*rRel(1) + rd_in(4)*rRel(2)                    + rd_in(3) ! z
      
      ! absolute accelerations
      rd_out(1) =                   - rd_in(6)*rRel(2) + rd_in(5)*rRel(3) + rd_in(1) ! x   
      rd_out(2) =  rd_in(6)*rRel(1)                    - rd_in(4)*rRel(3) + rd_in(2) ! y
      rd_out(3) = -rd_in(5)*rRel(1) + rd_in(4)*rRel(2)                    + rd_in(3) ! z



      !rRel = MATMUL(TransMat, rRelBody)        
      !H = getH(rRel)      
      !! absolute locations
      !r_out = rRel + r_in 
      !! absolute velocities
      !rd_out = MATMUL( H, rd_in(4:6)) + rd_in(1:3)  
      

   END SUBROUTINE TransformKinematics
   !-----------------------------------------------------------------------
   
   

   ! calculate position, velocity, and acceleration of point based on its position relative to moving 6DOF body
   !-----------------------------------------------------------------------
   SUBROUTINE TransformKinematicsA(rRelBody, r_in, TransMat, v_in, a_in, r_out, v_out, a_out)
      REAL(DbKi),       INTENT(IN   )  :: rRelBody(:)  ! relative location of point about reference point, in local/reference coordinate system
      REAL(DbKi),       INTENT(IN   )  :: r_in(3)      ! translation applied to reference point
      REAL(DbKi),       INTENT(IN   )  :: TransMat(3,3)! rotation matrix describing rotation about reference point
      REAL(DbKi),       INTENT(IN   )  :: v_in(6)      ! 6DOF velecity vector about ref point in global orientation frame
      REAL(DbKi),       INTENT(IN   )  :: a_in(6)      ! 6DOF acceleration vector
      REAL(DbKi),       INTENT(  OUT)  :: r_out(3)     ! coordinates of point of interest
      REAL(DbKi),       INTENT(  OUT)  :: v_out(3)     ! velocity of point
      REAL(DbKi),       INTENT(  OUT)  :: a_out(3)     ! acceleration of point

      REAL(DbKi)                       :: rRel(3)
!      REAL(DbKi)                       :: rRel2(3)

!      REAL(DbKi)                       :: r_out2(3)
!      REAL(DbKi)                       :: rd_out2(3)
      REAL(DbKi)                       :: H(3,3)

      ! rd_in should be in global orientation frame
      ! note: it's okay if r_out and rd_out are 6-size. Only the first 3 will be written, and 4-6 will
      !       already be correct or can be assigned separately from r_in and rd_in (assuming orientation frames are identical)


      ! locations about ref point in *unrotated* reference frame
      !rRel2(1) = TransMat(1,1)*rRelBody(1) + TransMat(2,1)*rRelBody(2) + TransMat(3,1)*rRelBody(3)  ! x
      !rRel2(2) = TransMat(1,2)*rRelBody(1) + TransMat(2,2)*rRelBody(2) + TransMat(3,2)*rRelBody(3)  ! y
      !rRel2(3) = TransMat(1,3)*rRelBody(1) + TransMat(2,3)*rRelBody(2) + TransMat(3,3)*rRelBody(3)  ! z
      
      rRel = MATMUL(TransMat, rRelBody)  
      
      H = getH(rRel)
      
      ! absolute locations
      r_out = rRel + r_in 

      ! absolute velocities
      !rd_out2(1) =                  - v_in(6)*rRel(2) + v_in(5)*rRel(3) + v_in(1) ! x   
      !rd_out2(2) =  v_in(6)*rRel(1)                   - v_in(4)*rRel(3) + v_in(2) ! y
      !rd_out2(3) = -v_in(5)*rRel(1) + v_in(4)*rRel(2)                   + v_in(3) ! z
      
      v_out = MATMUL( H, v_in(4:6)) + v_in(1:3)  
      
      ! absolute accelerations
      a_out = MATMUL( H, a_in(4:6)) + a_in(1:3)   ! << should add second order terms!
        

   END SUBROUTINE TransformKinematicsA
   !-----------------------------------------------------------------------
   
   ! calculate position and velocity of point along rod (distance L along direction u)
   !-----------------------------------------------------------------------
   SUBROUTINE TransformKinematicsAtoB(rA, u, L, rd_in, r_out, rd_out)
      REAL(DbKi),       INTENT(IN   )  :: rA(3)        ! coordinate of end A
      REAL(DbKi),       INTENT(IN   )  :: u(3)         ! Rod unit vector
      REAL(DbKi),       INTENT(IN   )  :: L            ! Rod length from end A to B
      REAL(DbKi),       INTENT(IN   )  :: rd_in(6)     ! 6DOF velecity vector about Rod end A, in global orientation frame
      REAL(DbKi),       INTENT(  OUT)  :: r_out(3)     ! coordinates of end B
      REAL(DbKi),       INTENT(  OUT)  :: rd_out(3)    ! velocity of end B

      REAL(DbKi)                       :: rRel(3)
      
      
      ! locations (unrotated reference frame)
      rRel = L*u              ! relative location of point B from point A
      r_out = rRel + rA           ! absolute location of point B
      
      ! absolute velocities
      rd_out(1) =                   - rd_in(6)*rRel(2) + rd_in(5)*rRel(3) + rd_in(1)  ! x   
      rd_out(2) =  rd_in(6)*rRel(1)                    - rd_in(4)*rRel(3) + rd_in(2)  ! y
      rd_out(3) = -rd_in(5)*rRel(1) + rd_in(4)*rRel(2)                    + rd_in(3)  ! z
      

   END SUBROUTINE TransformKinematicsAtoB
   !-----------------------------------------------------------------------
   
   !
   !-----------------------------------------------------------------------
   SUBROUTINE TranslateForce3to6DOF(dx, F, Fout)
      REAL(DbKi),       INTENT(IN   )  :: dx(3)       ! displacement vector from ref point to point of force (F) application
      REAL(DbKi),       INTENT(IN   )  :: F(3)        ! applied force
      REAL(DbKi),       INTENT(  OUT)  :: Fout(6)     ! resultant applied force and moment about ref point
        
      Fout(1:3) = F
      
      Fout(4:6) = CROSS_PRODUCT(dx, F)

   END SUBROUTINE TranslateForce3to6DOF
   !-----------------------------------------------------------------------
   
   
   !
   !-----------------------------------------------------------------------
   SUBROUTINE TranslateMass3to6DOF(dx, Min, Mout)
      REAL(DbKi),       INTENT(IN   )  :: dx(3)       ! displacement vector from ref point to point of mass matrix (Min)
      REAL(DbKi),       INTENT(IN   )  :: Min( 3,3)   ! original mass matrix (assumed at center of mass, or a point mass)
      REAL(DbKi),       INTENT(  OUT)  :: Mout(6,6)   ! resultant mass and inertia matrix about ref point
  
      REAL(DbKi)                       :: H(     3,3) ! "anti-symmetric tensor components" from Sadeghi and Incecik
!      REAL(DbKi)                       :: tempM( 3,3)
!      REAL(DbKi)                       :: tempM2(3,3)
!      REAL(DbKi)                       :: Htrans(3,3)
!      Integer(IntKi)                   :: I
   
      ! sub-matrix definitions are accordint to  | m    J |
      !                                          | J^T  I |
   
      H = getH(dx);

      ! mass matrix  [m'] = [m]
      Mout(1:3,1:3) = Min

      ! product of inertia matrix  [J'] = [m][H] + [J]
      Mout(1:3,4:6) = MATMUL(Min, H)
      Mout(4:6,1:3) = TRANSPOSE(Mout(1:3,4:6)) 

      !moment of inertia matrix  [I'] = [H][m][H]^T + [J]^T [H] + [H]^T [J] + [I]
      Mout(4:6,4:6) = MATMUL(MATMUL(H, Min), TRANSPOSE(H))

   END SUBROUTINE TranslateMass3to6DOF
   !-----------------------------------------------------------------------
   
   !
   !-----------------------------------------------------------------------
   SUBROUTINE TranslateMass6to6DOF(dx, Min, Mout)
      REAL(DbKi),       INTENT(IN   )  :: dx(3)       ! displacement vector from ref point to point of mass matrix (Min)
      REAL(DbKi),       INTENT(IN   )  :: Min( 6,6)   ! original mass matrix 
      REAL(DbKi),       INTENT(  OUT)  :: Mout(6,6)   ! resultant mass and inertia matrix about ref point
  
      REAL(DbKi)                       :: H(     3,3) ! "anti-symmetric tensor components" from Sadeghi and Incecik
         
      H = getH(dx);

      ! mass matrix  [m'] = [m]
      Mout(1:3,1:3) = Min(1:3,1:3)

      ! product of inertia matrix  [J'] = [m][H] + [J]
      Mout(1:3,4:6) = MATMUL(Min(1:3,1:3), H) + Min(1:3,4:6)
      Mout(4:6,1:3) = TRANSPOSE(Mout(1:3,4:6))      

      !moment of inertia matrix  [I'] = [H][m][H]^T + [J]^T [H] + [H]^T [J] + [I]
      Mout(4:6,4:6) = MATMUL(MATMUL(H, Min(1:3,1:3)), TRANSPOSE(H)) + MATMUL(Min(4:6,1:3),H) + MATMUL(TRANSPOSE(H),Min(1:3,4:6)) + Min(4:6,4:6)

   END SUBROUTINE TranslateMass6to6DOF
   !-----------------------------------------------------------------------
   
   ! produce alternator matrix
   !-----------------------------------------------------------------------
   FUNCTION GetH(r)
      Real(DbKi), INTENT(IN)      :: r(3)     ! inputted vector
      Real(DbKi)                  :: GetH(3,3) ! outputted matrix
      
      GetH(2,1) = -r(3)
      GetH(1,2) =  r(3)
      GetH(3,1) =  r(2)
      GetH(1,3) = -r(2)
      GetH(3,2) = -r(1)
      GetH(2,3) =  r(1)
      
      GetH(1,1) = 0.0_DbKi
      GetH(2,2) = 0.0_DbKi
      GetH(3,3) = 0.0_DbKi

   END FUNCTION GetH
   !-----------------------------------------------------------------------
   
   
   
   ! apply a rotation to a 6-by-6 mass/inertia tensor (see Sadeghi and Incecik 2005 for theory)
   !-----------------------------------------------------------------------
   FUNCTION RotateM6(Min, rotMat) result(outMat)
      
      Real(DbKi), INTENT(IN)      :: Min(6,6)     ! inputted matrix to be rotated
      Real(DbKi), INTENT(IN)      :: rotMat(3,3)  ! rotation matrix (DCM)
      Real(DbKi)                  :: outMat(6,6)  ! rotated matrix
      
      ! the process for each of the following is to 
      ! 1. copy out the relevant 3x3 matrix section,
      ! 2. rotate it, and
      ! 3. paste it into the output 6x6 matrix

      ! mass matrix
      outMat(1:3,1:3) = rotateM3(Min(1:3,1:3), rotMat)

      ! product of inertia matrix
      outMat(1:3,4:6) = rotateM3(Min(1:3,4:6), rotMat)      
      outMat(4:6,1:3) = TRANSPOSE(outMat(1:3,4:6))

      ! moment of inertia matrix
      outMat(4:6,4:6) = rotateM3(Min(4:6,4:6), rotMat)
   
   END FUNCTION RotateM6


   ! apply a rotation to a 3-by-3 mass matrix or any other second order tensor
   !-----------------------------------------------------------------------
   FUNCTION RotateM3(Min, rotMat) result(outMat)
      
      Real(DbKi), INTENT(IN)      :: Min(3,3)     ! inputted matrix to be rotated
      Real(DbKi), INTENT(IN)      :: rotMat(3,3)  ! rotation matrix (DCM)
      Real(DbKi)                  :: outMat(3,3)  ! rotated matrix
   
      ! overall operation is [m'] = [a]*[m]*[a]^T

      outMat = MATMUL( MATMUL(rotMat, Min), TRANSPOSE(rotMat) )

   END FUNCTION RotateM3
   
   
   
   
   
   ! calculates rotation matrix R to rotate from global axes to a member's local axes
   !-----------------------------------------------------------------------
   FUNCTION CalcOrientation(phi, beta, gamma) result(R)
      
      REAL(DbKi),      INTENT ( IN    )  :: phi     ! member incline angle
      REAL(DbKi),      INTENT ( IN    )  :: beta    ! member incline heading
      REAL(DbKi),      INTENT ( IN    )  :: gamma   ! member twist angle
      REAL(DbKi)                         :: R(3,3)  ! rotation matrix 
     
!      INTEGER(IntKi)                     :: errStat  
!      CHARACTER(100)                     :: errMsg   

      REAL(DbKi)                         :: s1, c1, s2, c2, s3, c3


      ! trig terms for Euler angles rotation based on beta, phi, and gamma
      s1 = sin(beta) 
      c1 = cos(beta)
      s2 = sin(phi) 
      c2 = cos(phi)
      s3 = sin(gamma) 
      c3 = cos(gamma)
      
      ! calculate rotation matrix based on Z1Y2Z3 Euler rotation sequence from https:!en.wikipedia.org/wiki/Euler_angles#Rotation_matrix

      R(1,1) = c1*c2*c3-s1*s3
      R(1,2) = -c3*s1-c1*c2*s3
      R(1,3) = c1*s2
      R(2,1) = c1*s3+c2*c3*s1
      R(2,2) = c1*c3-c2*s1*s3
      R(2,3) = s1*s2
      R(3,1) = -c3*s2
      R(3,2) = s2*s3
      R(3,3) = c2

      ! could also calculate unit normals p1 and p2 for rectangular cross sections   
      !p1 = matmul( R, [1,0,0] )               ! unit vector that is perpendicular to the 'beta' plane if gamma is zero
      !p2 = cross( q, p1 )                     ! unit vector orthogonal to both p1 and q
      
   END FUNCTION CalcOrientation
   

   !compute the inverse of a 3-by-3 matrix m
   !-----------------------------------------------------------------------
   SUBROUTINE Inverse3by3( Minv, M )
      Real(DbKi), INTENT(OUT)   :: Minv(3,3)  ! returned inverse matrix
      Real(DbKi), INTENT(IN)    :: M(3,3)     ! inputted matrix

      Real(DbKi)                :: det        ! the determinant
      Real(DbKi)                :: invdet     ! inverse of the determinant

      det = M(1, 1) * (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) - &
            M(1, 2) * (M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1)) + &
            M(1, 3) * (M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1));

      invdet = 1.0 / det   ! because multiplying is faster than dividing

      Minv(1, 1) = (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) * invdet
      Minv(1, 2) = (M(1, 3) * M(3, 2) - M(1, 2) * M(3, 3)) * invdet
      Minv(1, 3) = (M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2)) * invdet
      Minv(2, 1) = (M(2, 3) * M(3, 1) - M(2, 1) * M(3, 3)) * invdet
      Minv(2, 2) = (M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1)) * invdet
      Minv(2, 3) = (M(2, 1) * M(1, 3) - M(1, 1) * M(2, 3)) * invdet
      Minv(3, 1) = (M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2)) * invdet
      Minv(3, 2) = (M(3, 1) * M(1, 2) - M(1, 1) * M(3, 2)) * invdet
      Minv(3, 3) = (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)) * invdet

   END SUBROUTINE Inverse3by3
   !-----------------------------------------------------------------------


   ! One-function implementation of Crout LU Decomposition. Solves Ax=b for x
   SUBROUTINE LUsolve(n, A, LU, b, y, x)

      INTEGER(intKi),   INTENT(IN   )  :: n           ! size of matrices and vectors
      Real(DbKi),       INTENT(IN   )  :: A( n,n)     ! LHS matrix (e.g. mass matrix)
      Real(DbKi),       INTENT(INOUT)  :: LU(n,n)     ! stores LU matrix data
      Real(DbKi),       INTENT(IN   )  :: b(n)        ! RHS vector
      Real(DbKi),       INTENT(INOUT)  :: y(n)        ! temporary vector
      Real(DbKi),       INTENT(  OUT)  :: x(n)        ! LHS vector to solve for

      INTEGER(intKi)                   :: i,j,k,p
      Real(DbKi)                       :: sum

      DO k = 1,n
         DO i = k,n
         
            sum = 0.0_DbKi
            
            DO p=1,k-1   !for(int p=0; p<k; ++p)
               sum = sum + LU(i,p)*LU(p,k)
            END DO
            
            LU(i,k) = A(i,k) - sum
         END DO !I
         
         DO j=k+1,n  !for(int j=k+1;j<n;++j)
            
            sum = 0.0_DbKi

            DO p=1,k-1   !for(int p=0;p<k;++p)
               sum = sum + LU(k,p)*LU(p,j)
            END DO
            
            IF ( EqualRealNos( LU(k,k), 0.0_DbKi) ) THEN
               LU(k,j) = 0.0_DbKi   ! avoid divide by zero <<< numerator likely zero too. check IF this is safe <<<
            ELSE
               LU(k,j) = (A(k,j)-sum)/LU(k,k)
            END IF
         END DO !j
         
      END DO !K
      
      DO i=1,n
      
         sum = 0.0_DbKi
         
         DO k=1,i-1  !for(int k=0; k<i; ++k)
            sum = sum + LU(i,k)*y(k);
         END DO
         
         IF ( EqualRealNos( LU(i,i), 0.0_DbKi) ) THEN
            y(i) = 0.0_DbKi   ! avoid divide by zero <<< numerator likely zero too. check IF this is safe <<<
         ELSE
            y(i) = (b(i)-sum)/LU(i,i)
         END IF
      END DO
      
      DO j=1,n       ! this is actually for looping through i in reverse
         i = n+1-j   ! for(int i=n-1; i>=0; --i)
      
         sum = 0.0_DbKi
         
         DO k=i+1, n 
            sum = sum + LU(i,k)*x(k)
         END DO
         
         x(i) = (y(i)-sum) 
         
      END DO !j (actually decrementing i)
      
   END SUBROUTINE LUsolve


   
   ! :::::::::::::::::::::::::: interpolation subroutines :::::::::::::::::::::::::::::::
   
   
   SUBROUTINE getInterpNums(xlist, xin, istart, i, fout)
      
      Real(DbKi),    INTENT (IN   )            :: xlist(:)      ! list of x values
      Real(DbKi),    INTENT (IN   )            :: xin           ! x value to be interpolated
      Integer(IntKi),INTENT (IN   )            :: istart        ! first lower index to try
      Integer(IntKi),INTENT (  OUT)            :: i             ! lower index to interpolate from
      Real(DbKi),    INTENT (  OUT)            :: fout          ! fraction to return   such that y* = y[i] + fout*(y[i+1]-y[i])
      
      Integer(IntKi)                           :: i1
      Integer(IntKi)                           :: nx
      
      i1 = 1   ! Setting in declaration causes an implied save, which would never allow this routine to find anything at the start of the array.

      nx = SIZE(xlist)
      
      IF (xin <= xlist(1)) THEN                !  below lowest data point
         i = 1_IntKi
         fout = 0.0_DbKi
      
      ELSE IF (xlist(nx) <= xin) THEN          ! above highest data point
         i = nx
         fout = 0.0_DbKi
      
      ELSE                                     ! within the data range
     
         IF (xlist(min(istart,nx)) < xin) i1 = istart  ! if istart is below the actual value, start with it instead of starting at 1 to save time, but make sure it doesn't overstep the array
      
         DO i = i1, nx-1
            IF (xlist(i+1) > xin) THEN
               fout = (xin - xlist(i) )/( xlist(i+1) - xlist(i) )
               exit
            END IF
         END DO
      END IF
      
   END SUBROUTINE getInterpNums
   
   
   SUBROUTINE getInterpNumsSiKi(xlist, xin, istart, i, fout)
      
      Real(SiKi),    INTENT (IN   )            :: xlist(:)      ! list of x values
      Real(SiKi),    INTENT (IN   )            :: xin           ! x value to be interpolated
      Integer(IntKi),INTENT (IN   )            :: istart        ! first lower index to try
      Integer(IntKi),INTENT (  OUT)            :: i             ! lower index to interpolate from
      Real(SiKi),    INTENT (  OUT)            :: fout          ! fraction to return   such that y* = y[i] + fout*(y[i+1]-y[i])
      
      Integer(IntKi)                           :: i1
      Integer(IntKi)                           :: nx
      
      i1 = 1   ! Setting in declaration causes an implied save, which would never allow this routine to find anything at the start of the array.

      nx = SIZE(xlist)
      
      IF (xin <= xlist(1)) THEN                !  below lowest data point
         i = 1_IntKi
         fout = 0.0_SiKi
      
      ELSE IF (xlist(nx) <= xin) THEN          ! above highest data point
         i = nx
         fout = 0.0_SiKi
      
      ELSE                                     ! within the data range
     
         IF (xlist(min(istart,nx)) < xin) i1 = istart  ! if istart is below the actual value, start with it instead of starting at 1 to save time, but make sure it doesn't overstep the array
      
         DO i = i1, nx-1
            IF (xlist(i+1) > xin) THEN
               fout = (xin - xlist(i) )/( xlist(i+1) - xlist(i) )
               exit
            END IF
         END DO
      END IF
      
   END SUBROUTINE getInterpNumsSiKi
  
   SUBROUTINE calculate4Dinterpolation(f, ix0, iy0, iz0, it0, fx, fy, fz, ft, c)

      Real(SiKi),     INTENT (IN   )        :: f(:,:,:,:)                ! data array
      INTEGER(IntKi), INTENT (IN   )        :: ix0, iy0, iz0, it0        ! indices for interpolation
      Real(SiKi),     INTENT (IN   )        :: fx, fy, fz, ft            ! interpolation fractions
      Real(DbKi),     INTENT (  OUT)        :: c                         ! the output value
                                         
      INTEGER(IntKi)                        :: ix1, iy1, iz1, it1        ! second indices
      REAL(SiKi)                            :: c000, c001, c010, c011, c100, c101, c110, c111
      REAL(SiKi)                            :: c00, c01, c10, c11, c0, c1  
      
      ! handle end case conditions
      IF (fx == 0) THEN 
         ix1 = ix0
      ELSE  
         ix1 = min(ix0+1,size(f,4))    ! don't overstep bounds
      END IF
      
      IF (fy == 0) THEN
         iy1 = iy0
      ELSE
         iy1 = min(iy0+1,size(f,3))    ! don't overstep bounds
      END IF
      
      IF (fz == 0) THEN
         iz1 = iz0
      ELSE         
         iz1 = min(iz0+1,size(f,2))    ! don't overstep bounds
      END IF
      
      IF (ft == 0) THEN
         it1 = it0
      ELSE  
         it1 = min(it0+1,size(f,1))    ! don't overstep bounds
      END IF
      
      c000 = f(it0,iz0,iy0,ix0)*(1.0-ft) + f(it1,iz0,iy0,ix0)*ft
      c001 = f(it0,iz1,iy0,ix0)*(1.0-ft) + f(it1,iz1,iy0,ix0)*ft
      c010 = f(it0,iz0,iy1,ix0)*(1.0-ft) + f(it1,iz0,iy1,ix0)*ft
      c011 = f(it0,iz1,iy1,ix0)*(1.0-ft) + f(it1,iz1,iy1,ix0)*ft
      c100 = f(it0,iz0,iy0,ix1)*(1.0-ft) + f(it1,iz0,iy0,ix1)*ft
      c101 = f(it0,iz1,iy0,ix1)*(1.0-ft) + f(it1,iz1,iy0,ix1)*ft
      c110 = f(it0,iz0,iy1,ix1)*(1.0-ft) + f(it1,iz0,iy1,ix1)*ft
      c111 = f(it0,iz1,iy1,ix1)*(1.0-ft) + f(it1,iz1,iy1,ix1)*ft
      
      c00 = c000*(1.0-fx) + c100*fx
      c01 = c001*(1.0-fx) + c101*fx
      c10 = c010*(1.0-fx) + c110*fx
      c11 = c011*(1.0-fx) + c111*fx

      c0  = c00 *(1.0-fy) + c10 *fy
      c1  = c01 *(1.0-fy) + c11 *fy

      c   = c0  *(1.0-fz) + c1  *fz
            
   END SUBROUTINE


   SUBROUTINE calculate3Dinterpolation(f, ix0, iy0, iz0, fx, fy, fz, c)

         Real(SiKi),     INTENT (IN   )        :: f(:,:,:)                  ! data array
         INTEGER(IntKi), INTENT (IN   )        :: ix0, iy0, iz0             ! indices for interpolation
         Real(SiKi),     INTENT (IN   )        :: fx, fy, fz                ! interpolation fractions
         Real(DbKi),     INTENT (  OUT)        :: c                         ! the output value
         
         INTEGER(IntKi)                        :: ix1, iy1, iz1             ! second indices
         REAL(SiKi)                            :: c000, c001, c010, c011, c100, c101, c110, c111
         REAL(SiKi)                            :: c00, c01, c10, c11, c0, c1  
         
      ! note that "z" could also be "t" - dimension names are arbitrary

      ! handle end case conditions
      IF (fx == 0) THEN 
         ix1 = ix0
      ELSE  
         ix1 = min(ix0+1,size(f,3))    ! don't overstep bounds
      END IF
      
      IF (fy == 0) THEN
         iy1 = iy0
      ELSE
         iy1 = min(iy0+1,size(f,2))    ! don't overstep bounds
      END IF
      
      IF (fz == 0) THEN
         iz1 = iz0
      ELSE         
         iz1 = min(iz0+1,size(f,1))    ! don't overstep bounds
      END IF
      
      c000 = f(iz0,iy0,ix0)
      c001 = f(iz1,iy0,ix0)
      c010 = f(iz0,iy1,ix0)
      c011 = f(iz1,iy1,ix0)
      c100 = f(iz0,iy0,ix1)
      c101 = f(iz1,iy0,ix1)
      c110 = f(iz0,iy1,ix1)
      c111 = f(iz1,iy1,ix1)
      
      c00 = c000*(1.0-fx) + c100*fx
      c01 = c001*(1.0-fx) + c101*fx
      c10 = c010*(1.0-fx) + c110*fx
      c11 = c011*(1.0-fx) + c111*fx

      c0  = c00 *(1.0-fy) + c10 *fy
      c1  = c01 *(1.0-fy) + c11 *fy

      c   = c0  *(1.0-fz) + c1  *fz

   END SUBROUTINE

   SUBROUTINE calculate2Dinterpolation(f, ix0, iy0, fx, fy, c)
      REAL(DbKi),     INTENT (IN   )        :: f(:,:)                    ! data array
      INTEGER(IntKi), INTENT (IN   )        :: ix0, iy0                  ! indices for interpolation
      REAL(DbKi),     INTENT (IN   )        :: fx, fy                    ! interpolation fractions
      REAL(DbKi),     INTENT (  OUT)        :: c                         ! the output value
      
      INTEGER(IntKi)                        :: ix1, iy1                  ! second indices
      REAL(DbKi)                            :: c00, c01, c10, c11, c0, c1  

      ! handle end case conditions
      IF (fx == 0) THEN
         ix1 = ix0
      ELSE
         ix1 = min(ix0+1,size(f,2))  ! don't overstep bounds
      END IF
      IF (fy == 0) THEN
         iy1 = iy0
      ELSE
         iy1 = min(iy0+1,size(f,1))  ! don't overstep bounds
      END IF
      c00 = f(iy0, ix0)
      c01 = f(iy1, ix0)
      c10 = f(iy0, ix1)
      c11 = f(iy1, ix1)
      c0  = c00 *(1.0-fx) + c10 *fx
      c1  = c01 *(1.0-fx) + c11 *fx
      c   = c0  *(1.0-fy) + c1  *fy
   END SUBROUTINE calculate2Dinterpolation


   SUBROUTINE calculate1Dinterpolation(f, ix0, fx, c)
      REAL(DbKi),     INTENT (IN   )        :: f(:)                      ! data array
      INTEGER(IntKi), INTENT (IN   )        :: ix0                       ! indices for interpolation
      REAL(DbKi),     INTENT (IN   )        :: fx                        ! interpolation fractions
      REAL(DbKi),     INTENT (  OUT)        :: c                         ! the output value
      
      INTEGER(IntKi)                        :: ix1                       ! second index
      REAL(DbKi)                            :: c0, c1  

      ! handle end case conditions
      IF (fx == 0) THEN
         ix1 = ix0
      ELSE
         ix1 = min(ix0+1,size(f,1))  ! don't overstep bounds
      END IF

      c0 = f(ix0)
      c1 = f(ix1)
      c   = c0*(1.0-fx) + c1*fx
   END SUBROUTINE calculate1Dinterpolation




   ! :::::::::::::::::::::::::: bathymetry subroutines :::::::::::::::::::::::::::::::
   
   ! interpolates local seabed depth and normal vector
   SUBROUTINE getDepthFromBathymetry(BathymetryGrid, BathGrid_Xs, BathGrid_Ys, LineX, LineY, depth, nvec)

      REAL(DbKi),      INTENT(IN   )       :: BathymetryGrid(:,:) ! need colons or some sort of dimension setting
      REAL(DbKi),      INTENT(IN   )       :: BathGrid_Xs(:)
      REAL(DbKi),      INTENT(IN   )       :: BathGrid_Ys(:)
      REAL(DbKi),      INTENT(IN   )       :: LineX
      REAL(DbKi),      INTENT(IN   )       :: LineY 
      REAL(DbKi),      INTENT(  OUT)       :: depth      ! local seabed depth (positive down) [m]
      REAL(DbKi),      INTENT(  OUT)       :: nvec(3)       ! local seabed surface normal vector (positive out) 

      INTEGER(IntKi)                       :: ix0, iy0             ! indeces for interpolation   
      INTEGER(IntKi)                       :: ix1, iy1             ! second indices   
      Real(DbKi)                           :: fx, fy               ! interpolation fractions
      REAL(DbKi)                           :: c00, c01, c10, c11, cx0, cx1, c0y, c1y  ! temporary depth values
      Real(DbKi)                           :: dx, dy               ! x and y spacing of local grid panel [m]
      Real(DbKi)                           :: dc_dx, dc_dy         ! local slope
      Real(DbKi)                           :: tempVector(3)        ! normal vector before scaling to unit

      ! get interpolation indices and fractions for the relevant grid panel
      CALL getInterpNums(BathGrid_Xs, LineX, 1, ix0, fx)
      CALL getInterpNums(BathGrid_Ys, LineY, 1, iy0, fy)

      !CALL calculate2Dinterpolation(BathymetryGrid, ix, iy, fx, fy, depth)
      
      ! handle end case conditions
      IF (fx == 0) THEN
         ix1 = ix0
      ELSE
         ix1 = min(ix0+1,size(BathymetryGrid,2))  ! don't overstep bounds
      END IF
      IF (fy == 0) THEN
         iy1 = iy0
      ELSE
         iy1 = min(iy0+1,size(BathymetryGrid,1))  ! don't overstep bounds
      END IF
      
      ! get corner points of the panel
      c00 = BathymetryGrid(iy0, ix0)
      c01 = BathymetryGrid(iy1, ix0)
      c10 = BathymetryGrid(iy0, ix1)
      c11 = BathymetryGrid(iy1, ix1)
      
      ! get interpolated points and local value
      cx0    = c00 *(1.0-fx) + c10 *fx
      cx1    = c01 *(1.0-fx) + c11 *fx
      c0y    = c00 *(1.0-fy) + c01 *fy
      c1y    = c10 *(1.0-fy) + c11 *fy
      depth  = cx0 *(1.0-fy) + cx1 *fy
      
      ! get local slope
      dx = BathGrid_Xs(ix1) - BathGrid_Xs(ix0)
      dy = BathGrid_Ys(iy1) - BathGrid_Ys(iy0)
      IF ( dx > 0.0 ) THEN
         dc_dx = (c1y-c0y)/dx
      ELSE
         dc_dx = 0.0_DbKi   ! maybe this should raise an error
      END IF
      IF ( dy > 0.0 ) THEN
         dc_dy = (cx1-cx0)/dy
      ELSE
         dc_dy = 0.0_DbKi   ! maybe this should raise an error
      END IF
      
      tempVector(1) = dc_dx
      tempVector(2) = dc_dy
      tempVector(3) = 1.0_DbKi
      CALL ScaleVector( tempVector, 1.0_DbKi, nvec ) ! compute unit vector      

   END SUBROUTINE getDepthFromBathymetry
   
   
   ! :::::::::::::::::::::::::: wave and current subroutines :::::::::::::::::::::::::::::::
   
   
   ! master function to get wave/water kinematics at a given point -- called by each object 
   SUBROUTINE getWaterKin(p, m, x, y, z, t, U, Ud, zeta, PDyn, ErrStat, ErrMsg)
   
      ! This whole approach assuems that px, py, and pz are in increasing order.
      ! Wheeler stretching is built in for WaterKin 1 and 2.
   
      TYPE(MD_ParameterType),           INTENT (IN   )       :: p           ! MoorDyn parameters
      TYPE(MD_MiscVarType),             INTENT (INOUT)       :: m           ! MoorDyn misc data
      Real(DbKi),                       INTENT (IN   )       :: x           ! node position
      Real(DbKi),                       INTENT (IN   )       :: y           ! node position
      Real(DbKi),                       INTENT (IN   )       :: z           ! node position
      Real(DbKi),                       INTENT (IN   )       :: t           ! time
      Real(DbKi),                       INTENT (INOUT)       :: U(3)        ! fluid speed
      Real(DbKi),                       INTENT (INOUT)       :: Ud(3)       ! fluid acceleration 
      Real(DbKi),                       INTENT (INOUT)       :: zeta        ! water height above node 
      Real(DbKi),                       INTENT (INOUT)       :: PDyn        ! dynamic fluid pressure
      INTEGER(IntKi),                   INTENT(  OUT)        :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)        :: ErrMsg      ! Error message IF ErrStat /= ErrID_None

      ! outputs of WaveField_GetNodeWaveKin not in above list (single precision required by wave grid)
      Real(ReKi)                 :: xyz_sp(3) ! single precision
      Real(SiKi)                 :: WaveElev1 ! Unused by MD. per SeaSt_WaveField:  zeta = WaveElev  = WaveElev1 + WaveElev2
      Real(SiKi)                 :: WaveElev2 ! Unused by MD. per SeaSt_WaveField:  zeta = WaveElev  = WaveElev1 + WaveElev2
      Real(SiKi)                 :: zeta_sp  ! single precison
      Real(SiKi)                 :: U_sp(3)  ! single precision
      Real(SiKi)                 :: Ud_sp(3) ! single precision 
      Real(SiKi)                 :: FAMCF(3) ! Unused by MD
      Real(SiKi)                 :: PDyn_sp  ! single precision
      INTEGER(IntKi)             :: nodeInWater ! Unused by MD
      INTEGER(IntKi)             :: ErrStat2
      CHARACTER(ErrMsgLen)       :: ErrMsg2   
      
      ! interpolation variables for Waterkin = 1, 2
      INTEGER(IntKi)             :: ix, iy, iz, it        ! indices for interpolation      
      INTEGER(IntKi)             :: iz0, iz1              ! special indices for currrent interpolation  
      Real(SiKi)                 :: fx, fy, fz, ft        ! interpolation fractions
      Real(DbKi)                 :: zp                    ! zprime coordinate used for Wheeler stretching

      ! local 
      CHARACTER(120)             :: RoutineName = 'getWaterKin'   
   
      ErrStat = ErrID_None
      ErrMsg  = ""

      IF (p%WaterKin == 3 .AND. (.NOT. m%IC_gen)) THEN ! disable wavekin 3 during IC_gen, otherwise will never find steady state (becasue of waves)
         
         ! SeaState throws warning when queried location is out of bounds from the SeaState grid, so no need to handle here

         ! Pack all MD inputs to WaveGrid input data types (double to single)
         !   (only pos needed becasue time is double in wave field, all other are outputs that will be set by WaveField_GetNodeWaveKin)
         xyz_sp = REAL((/ x, y, z /),SiKi)

         ! for now we will force the node to be in the water (forceNodeInWater = True). Rods handle partial submergence seperately so they need to get information from SeaState 
         CALL WaveField_GetNodeWaveKin(p%WaveField, m%WaveField_m, t, xyz_sp, .TRUE., nodeInWater, WaveElev1, WaveElev2, zeta_sp, PDyn_sp, U_sp, Ud_sp, FAMCF, ErrStat2, ErrMsg2 ) ! outputs: nodeInWater, WaveElev1, WaveElev2, FAMCF all unused
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         ! Unpack all WaveGrid outputs to MD output data types (single to double)
         U = REAL(U_sp,DbKi)
         Ud = REAL(Ud_sp,DbKi)
         zeta = REAL(zeta_sp,DbKi)
         PDyn = REAL(PDyn_sp,DbKi)

      ELSEIF (p%WaterKin == 1 .OR. p%WaterKin == 2) THEN ! old or hybrid approach. SeaState contributions handeled in setupWaterKin, just proceed using old method

         ! If wave kinematics enabled, get interpolated values from grid
         IF (p%WaveKin > 0) THEN      

            ! find time interpolation indices and coefficients
            it = floor(t/ p%dtWave) + 1    ! add 1 because Fortran indexing starts at 1
            ft = (t - (it-1)*p%dtWave)/p%dtWave
            m%WaveTi = it                  
         
            ! find x-y interpolation indices and coefficients
            CALL getInterpNumsSiKi(p%pxWave   , REAL(x,SiKi),  1, ix, fx) ! wave grid
            CALL getInterpNumsSiKi(p%pyWave   , REAL(y,SiKi),  1, iy, fy) ! wave grid
            
            ! interpolate wave elevation
            CALL calculate3Dinterpolation(p%zeta, ix, iy, it, fx, fy, ft, zeta)   
            
            ! compute modified z coordinate to be used for interpolating velocities and accelerations with Wheeler stretching
            zp = ( z - zeta ) * p%WtrDpth/( p%WtrDpth + zeta )
            
            CALL getInterpNumsSiKi(p%pzWave   , REAL(zp,SiKi),  1, iz, fz) ! wave grid

            ! interpolate everything else
            CALL calculate4Dinterpolation(p%PDyn  , ix, iy, iz, it, fx, fy, fz, ft, PDyn)         
            CALL calculate4Dinterpolation(p%uxWave, ix, iy, iz, it, fx, fy, fz, ft, U(1)  )
            CALL calculate4Dinterpolation(p%uyWave, ix, iy, iz, it, fx, fy, fz, ft, U(2)  )
            CALL calculate4Dinterpolation(p%uzWave, ix, iy, iz, it, fx, fy, fz, ft, U(3)  )      
            CALL calculate4Dinterpolation(p%axWave, ix, iy, iz, it, fx, fy, fz, ft, Ud(1) )
            CALL calculate4Dinterpolation(p%ayWave, ix, iy, iz, it, fx, fy, fz, ft, Ud(2) )
            CALL calculate4Dinterpolation(p%azWave, ix, iy, iz, it, fx, fy, fz, ft, Ud(3) )
         
         ELSE ! set things to zero if wave kinematics not enabled
            U  = 0.0_DbKi
            Ud = 0.0_DbKi
            zeta = 0.0_DbKi
            PDyn = 0.0_DbKi

         ENDIF

         ! If current kinematics enabled, add interpolated current values from profile
         IF (p%Current > 0) THEN 
         
            CALL getInterpNumsSiKi(p%pzCurrent, REAL(z,SiKi), 1, iz0, fz)
                     
            IF (fz == 0) THEN  ! handle end case conditions
               iz1 = iz0
            ELSE
               iz1 = min(iz0+1,size(p%pzCurrent))  ! don't overstep bounds
            END IF
            
            ! Add the current velocities to the wave velocities (if any)
            U(1) = U(1) + (1.0-fz)*p%uxCurrent(iz0) + fz*p%uxCurrent(iz1)
            U(2) = U(2) + (1.0-fz)*p%uyCurrent(iz0) + fz*p%uyCurrent(iz1)
         END IF

      ELSEIF (p%WaterKin > 3) THEN 
         CALL SetErrStat(ErrID_Fatal, "Invalid p%WaterKin value found in getWaterKin", ErrStat, ErrMsg, RoutineName)
      
      ELSE ! set things to zero if Water Kinematics not enabled
         U  = 0.0_DbKi
         Ud = 0.0_DbKi
         zeta = 0.0_DbKi
         PDyn = 0.0_DbKi
      ENDIF

   END SUBROUTINE getWaterKin

   ! ----- process WaterKin input value, potentially reading wave inputs and generating wave field -----
   SUBROUTINE setupWaterKin(WaterKinString, p, Tmax, ErrStat, ErrMsg)

      CHARACTER(40),           INTENT(IN   )  :: WaterKinString      ! string describing water kinematics filename
      TYPE(MD_ParameterType),  INTENT(INOUT)  :: p                   ! Parameters
      REAL(ReKi),              INTENT(IN   )  :: Tmax
      INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat             ! Error status of the operation
      CHARACTER(*),            INTENT(  OUT)  :: ErrMsg              ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                   :: I, iIn, ix, iy, iz, numHdrLn
      INTEGER(IntKi)                   :: ntIn   ! number of time series inputs from file      
      INTEGER(IntKi)                   :: UnIn   ! unit number for coefficient input file
      INTEGER(IntKi)                   :: UnEcho       
      REAL(SiKi), ALLOCATABLE          :: pzCurrentTemp(:)   ! current depth increments read in from input file (positive-down at this stage)
      REAL(SiKi)                       :: uxCurrentTemp(100)
      REAL(SiKi)                       :: uyCurrentTemp(100)
      TYPE(Current_InitOutputType)     :: Current_InitOutput     ! Current outputs from SS subroutine Current_Init (for current mod = 2)
      
      CHARACTER(120)                   :: tmpString   
      CHARACTER(120)                   :: WaveKinFile   
      INTEGER(IntKi)                   :: UnElev  ! unit number for coefficient input file
      REAL(SiKi), ALLOCATABLE          :: WaveTimeIn(:)  ! temporarily holds wave input time series
      REAL(SiKi), ALLOCATABLE          :: WaveElevIn(:)
      REAL(SiKi), ALLOCATABLE          :: WaveElev0(:)   ! interpolated reference wave elevation time series
      REAL(SiKi)                       :: WaveDir
      REAL(SiKi)                       :: t, Frac
      CHARACTER(1024)                  :: FileName             ! Name of MoorDyn input file  
      CHARACTER(120)                   :: Line
      CHARACTER(4096)                  :: entries2  
      INTEGER(IntKi)                   :: coordtype
      LOGICAL                          :: dataBegin
   
      INTEGER(IntKi)                   :: NStepWave    ! 
      INTEGER(IntKi)                   :: NStepWave2   ! 
      REAL(SiKi)                       :: WaveTMax     ! max wave elevation time series duration after optimizing length for FFT
      REAL(SiKi)                       :: WaveDOmega   ! frequency step
      REAL(SiKi), ALLOCATABLE          :: SinWaveDir(:)                                      ! SIN( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction.
      REAL(SiKi), ALLOCATABLE          :: CosWaveDir(:)                                      ! COS( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction.
      LOGICAL                          :: WaveMultiDir = .FALSE. ! Indicates the waves are multidirectional -- set by WaveField pointer if enabled

      REAL(SiKi),  ALLOCATABLE         :: TmpFFTWaveElev(:)     ! Data for the FFT calculation
      TYPE(FFT_DataType)               :: FFT_Data              ! the instance of the FFT module we're using
      
      REAL(SiKi)                       :: tmpReal            ! A temporary real number
      COMPLEX(SiKi),ALLOCATABLE        :: tmpComplex(:)      ! A temporary array (0:NStepWave2-1) for FFT use. 
   
      REAL(SiKi)                       :: Omega                 ! Wave frequency (rad/s)
      COMPLEX(SiKi), PARAMETER         :: ImagNmbr = (0.0,1.0)  ! The imaginary number, SQRT(-1.0)
      COMPLEX(SiKi)                    :: ImagOmega             ! = ImagNmbr*Omega (rad/s)
      REAL(DbKi), ALLOCATABLE          :: WaveNmbr(:)           ! wave number for frequency array
      REAL(SiKi),  ALLOCATABLE         :: WaveDirArr(:)         ! Wave direction array. Each frequency has a unique direction of WaveNDir > 1 (degrees). 0's for WaveKin = 1 or if disabled in SeaState.
      REAL(SiKi), ALLOCATABLE          :: WaveElevC0(:,:)       ! Discrete Fourier transform of the instantaneous elevation of incident waves at the ref point (meters)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveElevC( :)         ! Discrete Fourier transform of the instantaneous elevation of incident waves at the ref point (meters)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveAccCHx(:)         ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveAccCHy(:)         ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveAccCV( :)         ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveDynPC( :)         ! Discrete Fourier transform of the instantaneous dynamic pressure                       of incident waves before applying stretching at the zi-coordinates for points (N/m^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveVelCHx(:)         ! Discrete Fourier transform of the instantaneous horizontal velocity  in x-direction    of incident waves before applying stretching at the zi-coordinates for points (m/s)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveVelCHy(:)         ! Discrete Fourier transform of the instantaneous horizontal velocity in y-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveVelCV( :)         ! Discrete Fourier transform of the instantaneous vertical   velocity      of incident waves before applying stretching at the zi-coordinates for points (m/s)
!      COMPLEX(SiKi)                    :: WGNC                  ! Discrete Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component (-)

      INTEGER(IntKi)                   :: ErrStat2
      CHARACTER(ErrMsgLen)             :: ErrMsg2
      CHARACTER(120)                   :: RoutineName = 'SetupWaveKin'   

      ErrStat2 = ErrID_None
      ErrMsg2  = ""

      IF (LEN_TRIM(WaterKinString) == 0) THEN
         ! If the input is empty (not provided), there are no water kinematics to be included
         p%WaveKin  = 0
         p%Current  = 0
         p%WaterKin = 0
         RETURN
         
      ELSE IF (SCAN(WaterKinString, "abcdfghijklmnopqrstuvwxyzABCDFGHIJKLMNOPQRSTUVWXYZ") == 0) THEN
         ! If the input has no letters, let's assume it's a number         
         IF (p%writeLog > 0) THEN
            WRITE(p%UnLog, '(A)'        ) "   ERROR WaveKin option does not currently support numeric entries. It must be a filename."
         ENDIF
         p%WaveKin  = 0
         p%Current  = 0
         p%WaterKin = 0
         CALL SetErrStat( ErrID_Fatal, "ERROR WaveKin option does not currently support numeric entries. It must be a filename.", ErrStat, ErrMsg, RoutineName); RETURN
      END IF

      tmpString = WaterKinString
      CALL Conv2UC(tmpString)

      ! Error check for SeaState rhoW and MoorDyn rhoW, SeaState gravity and MoorDyn gravity

      IF (tmpString == "SEASTATE") THEN ! if its the SeaState Keyword, then we are using the SeaState method (WaterKin = 3)

         ! Error check to make sure we are not running in FAST.Farm
         IF (p%nTurbines > 1) THEN ! if running in FAST.Farm, cannot use the SeaState coupling method because of multiple SeaState instances and only a single MoorDyn instance
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   MoorDyn fully coupled with SeaState is not compatible with FAST.Farm. Please use the WaterKin 1 or 2 methods."
            ENDIF
            CALL SetErrStat(ErrID_Fatal, "MoorDyn fully coupled with SeaState is not compatible with FAST.Farm. Please use the WaterKin 1 or 2 methods.", ErrStat, ErrMsg, RoutineName); RETURN
         END IF

         ! Error check to make sure the wave field pointer is not null
         IF (.NOT. ASSOCIATED(p%WaveField)) THEN
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   ERROR WaveField pointer is null. SeaState method requires SeaState to be enabled. Please check input files."
            ENDIF
            CALL SetErrStat(ErrID_Fatal, "WaveField pointer is null. SeaState method requires SeaState to be enabled. Please check input files.", ErrStat, ErrMsg, RoutineName); RETURN
         END IF

         ! Warning check to make sure SeaState and MoorDyn have the same water depth
         IF (p%WaveField%WtrDpth /= p%WtrDpth) THEN
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   WARNING SeaState water depth does not match MoorDyn water depth. Using SeaState values for water kinematics."
            ENDIF
            CALL SetErrStat(ErrID_Warn, "SeaState water depth does not match MoorDyn water depth. Using SeaState values for water kinematics.", ErrStat, ErrMsg, RoutineName)
         END IF

         ! Error check to make sure SeaState and MoorDyn have the same water density and gravity. This also checks that the pointer is valid
         IF (p%WaveField%RhoXg /= REAL((p%rhoW * p%g), SiKi)) THEN
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   ERROR SeaState (water density * gravity) ["//trim(num2lstr(p%WaveField%RhoXg))//"] does not match MoorDyn water (density * gravity) ["//trim(num2lstr(REAL((p%rhoW * p%g), SiKi)))//"]. The SeaState pointer may be corrupted."
            ENDIF
            CALL SetErrStat(ErrID_Fatal, "   ERROR SeaState (water density * gravity) ["//trim(num2lstr(p%WaveField%RhoXg))//"] does not match MoorDyn water (density * gravity) ["//trim(num2lstr(REAL((p%rhoW * p%g), SiKi)))//"]. The SeaState pointer may be corrupted.", ErrStat, ErrMsg, RoutineName)
         END IF

         ! Check for if SeaState grid does not match water depth
         IF (p%WaveField%GridParams%Z_Depth /= p%WtrDpth) THEN
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   INFO SeaState grid depth does not match MoorDyn water depth."
            ENDIF
         END IF

         ! turn on flags, for getWaterKin to know to just pull data from the WaveField pointer. Nothing more is needed for setup
         p%WaveKin  = 3
         p%Current  = 3
         p%WaterKin = 3

      ELSE ! we must be using the old or hybrid methods. Start setting up the user grids by reading the WaveKin File

         ! otherwise interpret the input as a file name to load the bathymetry lookup data from
         CALL WrScr( "   The waterKin input contains a filename so will load a water kinematics input file" )
         IF (p%writeLog > 0) THEN
            WRITE(p%UnLog, '(A)'        ) "   The waterKin input contains a filename so will load a water kinematics input file"
         ENDIF
         
         
         ! -------- load water kinematics input file -------------
         
         IF ( PathIsRelative( WaterKinString ) ) THEN   ! properly handle relative path <<<
            !CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
            FileName = TRIM(p%PriPath)//TRIM(WaterKinString)
         ELSE
            FileName = trim(WaterKinString)
         END IF
         


         UnEcho=-1
         CALL GetNewUnit( UnIn )   
         CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2); IF(Failed()) RETURN


         CALL ReadCom( UnIn, FileName, 'MoorDyn water kinematics input file header', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         CALL ReadCom( UnIn, FileName, 'MoorDyn water kinematics input file header', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         ! ----- waves -----
         CALL ReadCom( UnIn, FileName,                               'waves header', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         CALL ReadVar( UnIn, FileName, p%WaveKin  , 'WaveKinMod' ,  'WaveKinMod'   , ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         ! log the method being used
         IF (p%writeLog > 0) THEN
            IF (p%WaveKin == 2) THEN
               WRITE(p%UnLog, '(A)'        ) "    WaveKinMod = 2. Reading in the user provided wave grid and using SeaState for frequency calculation."
            ELSE IF (p%WaveKin == 1) THEN
               WRITE(p%UnLog, '(A)'        ) "    WaveKinMod = 1. Reading in the user provided wave grid and frequency information."
            ELSE IF (p%WaveKin == 0) THEN
               WRITE(p%UnLog, '(A)'        ) "    WaveKinMod = 0. No wave kinematics enabled."
            ELSE 
               WRITE(p%UnLog, '(A)'        ) "    Invalid value for WaveKinMod"
               Call SetErrStat(ErrID_Fatal, "Invalid value for WaveKinMod", ErrStat, ErrMsg, RoutineName); RETURN
            ENDIF
         ENDIF
         CALL ReadVar( UnIn, FileName, WaveKinFile, 'WaveKinFile',  'WaveKinFile'  , ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         CALL ReadVar( UnIn, FileName, p%dtWave   , 'dtWave', 'time step for waves', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         CALL ReadVar( UnIn, FileName, WaveDir    , 'WaveDir'    , 'wave direction', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         ! X grid points
         CALL ReadVar( UnIn, FileName, coordtype   , 'coordtype'   , '', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN        ! get the entry type
         READ(UnIn, '(A)', IOSTAT=ErrStat2) entries2; IF(ErrStat2 /= ErrID_None) ErrMsg2 = "There was an error reading in the grid description"; IF(Failed()) RETURN ! get entries as string to be processed
         CALL gridAxisCoords(coordtype, entries2, p%pxWave, p%nxWave); IF(Failed()) RETURN 
         ! Y grid points
         CALL ReadVar( UnIn, FileName, coordtype   , 'coordtype'   , '', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN        ! get the entry type
         READ(UnIn, '(A)', IOSTAT=ErrStat2) entries2; IF(ErrStat2 /= ErrID_None) ErrMsg2 = "There was an error reading in the grid description"; IF(Failed()) RETURN ! get entries as string to be processed
         CALL gridAxisCoords(coordtype, entries2, p%pyWave, p%nyWave); IF(Failed()) RETURN 
         ! Z grid points
         CALL ReadVar( UnIn, FileName, coordtype   , 'coordtype'   , '', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN        ! get the entry type
         READ(UnIn, '(A)', IOSTAT=ErrStat2) entries2; IF(ErrStat2 /= ErrID_None) ErrMsg2 = "There was an error reading in the grid description"; IF(Failed()) RETURN ! get entries as string to be processed
         CALL gridAxisCoords(coordtype, entries2, p%pzWave, p%nzWave); IF(Failed()) RETURN 
         ! ----- current -----
         CALL ReadCom( UnIn, FileName,                        'current header', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         CALL ReadVar( UnIn, FileName, p%Current,   'CurrentMod', 'CurrentMod', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN
         ! read in the current profile
         IF (p%Current == 2) THEN

            ! log the method being used
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "    CurrentMod = 2. Reading in the user provided current grid and using SeaState for profile calculation."
            ENDIF

            ! Z grid points
            CALL ReadVar( UnIn, FileName, coordtype   , 'coordtype'   , '', ErrStat2, ErrMsg2, UnEcho); IF(Failed()) RETURN         ! get the entry type
            READ(UnIn, '(A)', IOSTAT=ErrStat2) entries2; IF(ErrStat2 /= ErrID_None) ErrMsg2 = "There was an error reading in the grid description"; IF(Failed()) RETURN ! get entries as string to be processed
            CALL gridAxisCoords(coordtype, entries2, pzCurrentTemp, p%nzCurrent); IF(Failed()) RETURN  ! max size of 100 because gridAxisCoords has a 100 element temporary array used for processing entries2 
            uxCurrentTemp = 0.0 ! set these to zero to avoid unitialized values. This will be set later in this routine by SeaState
            uyCurrentTemp = 0.0 ! set these to zero to avoid unitialized values. This will be set later in this routine by SeaState

         ELSE IF (p%Current == 1) THEN 

            ! read two header lines. Old file format has 2 header lines, new file format has four (two info lines for CurrentMod = 2 and then the two headers).
            CALL ReadCom( UnIn, FileName,                'current profile header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
            CALL ReadCom( UnIn, FileName,                'current profile header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
            
            CALL AllocAry(pzCurrentTemp, 100, 'pzCurrentTemp', ErrStat2, ErrMsg2 ); IF(Failed()) RETURN ! allocate pzCurrentTemp if reading in user provided current table
            
            DO I=1,4 ! ignore lines until table is found, or exit after 3 lines. Old file format has 2 header lines, new file format has four (two info lines for CurrentMod = 2 and then the two headers).
               READ(UnIn, *, IOSTAT=ErrStat2) pzCurrentTemp(1), uxCurrentTemp(1), uyCurrentTemp(1)     ! try to read into a line to first elements in the array     
               IF (ErrStat2 == 0) THEN
                  EXIT      ! break out of the loop if successfully reads first table line
               END IF
            ENDDO

            IF (I == 4) THEN
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "   ERROR: Could not read the current profile table from the input file. Check the file format."
               ENDIF
               CALL SetErrStat(ErrID_Fatal, "Could not read the current profile table from the input file. Check the file format.", ErrStat, ErrMsg, RoutineName); RETURN
            ENDIF

            ! log the first line read
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "    CurrentMod = 1. Reading in the user provided current profile."
               IF (p%writeLog > 1) THEN
                  WRITE(p%UnLog, '(A)'        ) "     The first line read from the current table is:"
                  WRITE(p%UnLog, '(A)'        ) "     "//trim(num2lstr(pzCurrentTemp(1)))//"     "//trim(num2lstr(uxCurrentTemp(1)))//"     "//trim(num2lstr(uyCurrentTemp(1)))
               ENDIF
            ENDIF

            ! current profile table... (read until no more rows in table)
            DO I=2,100 ! start at second row because first is already read
               READ(UnIn, *, IOSTAT=ErrStat2) pzCurrentTemp(i), uxCurrentTemp(i), uyCurrentTemp(i)     ! read into a line      
               IF (ErrStat2 /= ErrID_None) THEN
                  p%nzCurrent = i-1 ! save number of valid current depth points in profile
                  EXIT      ! break out of the loop if it couldn't read the line (i.e. if at end of file)
               END IF
               IF (i == 100) THEN
                  p%nzCurrent = 100
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "    WARNING: MD can handle a maximum of 100 current profile points"
                  ENDIF      
                  CALL SetErrStat(ErrID_Warn, "MD can handle a maximum of 100 current profile points", ErrStat, ErrMsg, RoutineName)
                  EXIT
               END IF
            END DO

            ! Check that all z values are below the water line
            DO I=1,p%nzCurrent
               IF (pzCurrentTemp(I) > 0) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "    WARNING: MoorDyn current profile z values are above the water line."
                  ENDIF
                  CALL SetErrStat(ErrID_Warn, "Mooring current profile z values are above the water line.", ErrStat, ErrMsg, RoutineName)
               END IF
            END DO

            ! Check that order is valid (deepest to shallowest depths), flip array if necessary 
            IF (pzCurrentTemp(1) > pzCurrentTemp(p%nzCurrent)) THEN
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "    INFO: Current profile is decreasing depths. MoorDyn needs increasing depths. Flipping arrays"
               ENDIF
               DO I=1,p%nzCurrent/2
                  tmpReal = pzCurrentTemp(I)
                  pzCurrentTemp(I) = pzCurrentTemp(p%nzCurrent-I+1)
                  pzCurrentTemp(p%nzCurrent-I+1) = tmpReal
                  tmpReal = uxCurrentTemp(I)
                  uxCurrentTemp(I) = uxCurrentTemp(p%nzCurrent-I+1)
                  uxCurrentTemp(p%nzCurrent-I+1) = tmpReal
                  tmpReal = uyCurrentTemp(I)
                  uyCurrentTemp(I) = uyCurrentTemp(p%nzCurrent-I+1)
                  uyCurrentTemp(p%nzCurrent-I+1) = tmpReal
               END DO
            END IF

            ! Check for valid array (strictly increasing z values)
            DO I=1,p%nzCurrent-1
               IF (pzCurrentTemp(I) > pzCurrentTemp(I+1)) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   ERROR: Current profile z values are not strictly increasing. Check the file format."
                  ENDIF
                  CALL SetErrStat(ErrID_Fatal, "Current profile z values are not strictly increasing. Check the file format.", ErrStat, ErrMsg, RoutineName); RETURN
               END IF
            END DO

         ELSE IF (p%Current == 0) THEN
            ! log the method being used
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "    CurrentMod = 0. No currents will be simulated."
            ENDIF
         ELSE 
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "    Invalid value for CurrentMod"
            ENDIF
            CALL SetErrStat(ErrID_Fatal, "Invalid value for CurrentMod", ErrStat, ErrMsg, RoutineName); RETURN
         ENDIF ! if p%Current != 1 or 2, then no current
         

         CLOSE(UnIn)

         IF (p%Current > 2 .OR. p%WaveKin > 2) THEN
            CALL SetErrStat( ErrID_Fatal,'WaveKinMod and CurrentMod must be less than 3 if using user defined MoorDyn grid',ErrStat, ErrMsg, RoutineName); RETURN
            RETURN
         END IF

         ! Error check to make sure the wave field pointer is not null if CurrentMod = 2 or WaveKinMod = 2
         IF (p%WaveKin == 2 .OR. p%Current == 2) THEN
            IF (.NOT. ASSOCIATED(p%WaveField)) THEN
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "   ERROR WaveField pointer is null. Hybrid method requires SeaState to be enabled. Please check input files."
               ENDIF
               CALL SetErrStat(ErrID_Fatal, "WaveField pointer is null. Hybrid method requires SeaState to be enabled. Please check input files.", ErrStat, ErrMsg, RoutineName); RETURN
            END IF
         END IF

         ! set p%WaterKin, note that p%WaterKin is not used in this routine, but is necessary for getWaterKin
         IF (p%Current == p%WaveKin) THEN ! if they are the same, set WaveKin as that value
            p%WaterKin = p%WaveKin
         ELSEIF (p%Current == 0) THEN ! if one is zero, use the other for water kin
            p%WaterKin = p%WaveKin
         ELSEIF (p%WaveKin == 0) THEN ! if one is zero, use the other for water kin
            p%WaterKin = p%Current
         ELSEIF (p%WaveKin == 2 .AND. p%Current == 1 .AND. p%WaveField%Current_InitInput%CurrMod == 0) THEN ! if WaveKin = 2 and SeaState has no currents, allow users to provide custom current profile. WaterKin = 2 becasue of waves
            p%WaterKin = p%WaveKin
         ELSE
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   ERROR WaveKinMod and CurrentMod must be equal, one must be zero, or WaveKinMod must be 2 and CurrentMod must be 1 with SeaState currents disabled."
            ENDIF
            CALL SetErrStat( ErrID_Fatal,'WaveKinMod and CurrentMod must be equal, one must be zero, or WaveKinMod must be 2 and CurrentMod must be 1 with SeaState currents disabled',ErrStat, ErrMsg, RoutineName); RETURN
         ENDIF
            
         ! ------------------- start with wave kinematics -----------------------
         
         IF (p%WaveKin > 0) THEN 

            ! Check that all wave grid z values are below the water line, otherwise COSHNumOvrCOSHDen calcs will nan
            DO I=1,p%nzWave
               IF (p%pzWave(I) > 0) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "    ERROR: MoorDyn wave grid z values are above the water line. The wave grid cannot contain values greater than 0. This may be caused by interpolation over integer z-grid spacing."
                  ENDIF
                  CALL SetErrStat(ErrID_Fatal, "Mooring wave grid z values are above the water line. The wave grid cannot contain values greater than 0.", ErrStat, ErrMsg, RoutineName); RETURN
               END IF
            END DO
            
            IF (p%WaveKin == 2) THEN ! must be using the hybrid approach
      
               ! Warning check to make sure SeaState and MoorDyn have the same water depth
               IF (p%WaveField%WtrDpth /= p%WtrDpth) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   WARNING SeaState water depth does not match MoorDyn water depth. Using MoorDyn values for interpolating SeaState data to MoorDyn grid."
                  ENDIF
                  CALL SetErrStat(ErrID_Warn, "SeaState water depth does not match MoorDyn water depth. Using MoorDyn values for interpolating SeaState data to MoorDyn grid.", ErrStat, ErrMsg, RoutineName)
               END IF
      
               ! Error check to make sure SeaState and MoorDyn have the same water density and gravity. This also checks the pointer is valid.
               IF (p%WaveField%RhoXg /= REAL((p%rhoW * p%g), SiKi)) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   ERROR SeaState (water density * gravity) ["//trim(num2lstr(p%WaveField%RhoXg))//"] does not match MoorDyn water (density * gravity) ["//trim(num2lstr(REAL((p%rhoW * p%g), SiKi)))//"]. The SeaState pointer may be corrupted."
                  ENDIF
                  CALL SetErrStat(ErrID_Fatal, "   ERROR SeaState (water density * gravity) ["//trim(num2lstr(p%WaveField%RhoXg))//"] does not match MoorDyn water (density * gravity) ["//trim(num2lstr(REAL((p%rhoW * p%g), SiKi)))//"]. The SeaState pointer may be corrupted.", ErrStat, ErrMsg, RoutineName)
               END IF

               ! Warning check to make sure SeaState and MoorDyn have the same wave dir. For now, no wave spreading. This can be updated
               IF (p%WaveField%WaveDir /= WaveDir) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   WARNING SeaState WaveDir does not match MoorDyn WaveDir. Using MoorDyn values for interpolating SeaState data to MoorDyn grid."
                  ENDIF
                  CALL SetErrStat(ErrID_Warn, "SeaState WaveDir does not match MoorDyn WaveDir.Using MoorDyn values for interpolating SeaState data to MoorDyn grid.", ErrStat, ErrMsg, RoutineName)
               END IF

               ! Info check for if MoorDyn dtWave is non-zero. Users may set this accidentially, or could be left over in an input file
               IF (p%dtWave > 0) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   MoorDyn dtWave is ignored when using WaveKinMod = 2 becasue wave frequency information is supplied by SeaState"
                  ENDIF
               END IF

               ! check SS only is being run w/ first order waves, because MoorDyn is not compatiable with those(for now)
               IF (ALLOCATED(p%WaveField%WaveElev2)) THEN 
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   ERROR MoorDyn WaveKinMod2 does not support second order waves. Please use disable in SeaState."
                  ENDIF
                  CALL SetErrStat(ErrID_Fatal, "MoorDyn WaveKinMod2 does not support second order waves. Please use disable in SeaState.", ErrStat, ErrMsg, RoutineName); RETURN
               END IF

               ! set dtWave to WaveField dtWave (for interpolation use in getWaterKin)
               p%dtWave = p%WaveField%WaveTime(2)-p%WaveField%WaveTime(1) ! from Waves.f90 (line 1040),  DO I = 0,WaveField%NStepWave; WaveField%WaveTime(I) = I*REAL(InitInp%WaveDT,SiKi)

               ! Interpolations need the following from SeaState: WaveElevC0, NStepWave2, NStepWave, WaveDOmega. In p%WaveKin = 1, these values are extracted from the wave elevation time series
               ! Note: allocations not needed here because they are already allocated in SeaState
               WaveElevC0 = p%WaveField%WaveElevC0
               WaveDOmega = p%WaveField%WaveDOmega
               NStepWave2 = p%WaveField%NStepWave2
               NStepWave  = p%WaveField%NStepWave

               ! Pull some other things out of the WaveField pointer
               WaveMultiDir = p%WaveField%WaveMultiDir
               p%ntWave = NStepWave ! set ntWave to NStepWave

               ! Set wave spreading array if enabled in SeaState, otherwise set to zero
               If (WaveMultiDir) THEN
                  ! Note: allocations not needed here because they are already allocated in SeaState
                  WaveDirArr = p%WaveField%WaveDirArr
               ENDIF

            ELSEIF (p%WaveKin == 1) THEN ! must be a filepath therefore read wave elevations from timeseries

               ! NOTE: there is a decent ammount of code duplication (intentional for now) with what is in SeaState that eventually 
               ! we will want to synchronize with a standard library at some point

               ! --------------------- set from inputted wave elevation time series, grid approach -------------------
               CALL WrScr( '    WaveKinMod = 1. Reading wave elevation time series from file' )
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "    WaveKinMod = 1. Reading wave elevation time series from file"
               ENDIF
               
               IF ( LEN_TRIM( WaveKinFile ) == 0 )  THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   ERROR WaveKinFile must not be an empty string."
                  ENDIF
                  CALL SetErrStat( ErrID_Fatal,'WaveKinFile must not be an empty string.',ErrStat, ErrMsg, RoutineName); RETURN
                  RETURN
               END IF

               IF ( PathIsRelative( WaveKinFile ) ) THEN   ! properly handle relative path <<<
                  !CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
                  WaveKinFile = TRIM(p%PriPath)//TRIM(WaveKinFile)
               END IF
               
               ! note: following is adapted from MoorDyn_Driver
               
               CALL GetNewUnit( UnElev ) 
            
               CALL OpenFInpFile ( UnElev, WaveKinFile, ErrStat2, ErrMsg2 ); IF(Failed()) RETURN
            
               CALL WrScr( '    Reading wave elevation data from '//trim(WaveKinFile) )
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) '    Reading wave elevation data from '//trim(WaveKinFile)
               ENDIF
               
               ! Read through length of file to find its length
               i         = 0 ! start line counter
               numHdrLn  = 0 ! start header-line counter
               dataBegin = .FALSE. ! started reading the data section
               DO
                  READ(UnElev,'(A)',IOSTAT=ErrStat2) Line     !read into a line
                  IF (ErrStat2 /= 0) EXIT      ! break out of the loop if it couldn't read the line (i.e. if at end of file)
                  i = i+1
                  READ(Line,*,IOSTAT=ErrStat2) tmpReal
                  IF (ErrStat2/=0) THEN  ! Not a number
                     IF (dataBegin) THEN
                        IF (p%writeLog > 0) THEN
                           WRITE(p%UnLog, '(A)'        ) "   ERROR Non-data line detected in WaveKinFile past the header lines."
                        ENDIF
                        CALL SetErrStat( ErrID_Fatal,' Non-data line detected in WaveKinFile past the header lines.',ErrStat, ErrMsg, RoutineName); RETURN
                     END IF
                     numHdrLn = numHdrLn + 1
                  ELSE
                     dataBegin = .TRUE.
                  END IF
               END DO

               ! rewind to start of input file to re-read things now that we know how long it is
               REWIND(UnElev)      

               ntIn = i-numHdrLn     ! save number of lines of file
               

               ! allocate space for input wave elevation array (including time column)
               CALL AllocAry(WaveTimeIn,  ntIn, 'WaveTimeIn', ErrStat2, ErrMsg2 ); IF(Failed()) RETURN
               CALL AllocAry(WaveElevIn,  ntIn, 'WaveElevIn', ErrStat2, ErrMsg2 ); IF(Failed()) RETURN

               ! read the data in from the file
               DO i = 1, numHdrLn
                  READ(UnElev,'(A)',IOSTAT=ErrStat2) Line     ! skip header lines
               END DO
               
               DO i = 1, ntIn
                  READ (UnElev, *, IOSTAT=ErrStat2) WaveTimeIn(i), WaveElevIn(i) ! read wave elevation time series
                     
                  IF ( ErrStat2 /= 0 ) THEN
                     IF (p%writeLog > 0) THEN
                        WRITE(p%UnLog, '(A)'        ) "   ERROR reading WaveElev input file."
                     ENDIF
                     CALL SetErrStat( ErrID_Fatal,'Error reading WaveElev input file.',ErrStat, ErrMsg, RoutineName); RETURN
                  END IF 
               END DO  

               ! Close the inputs file 
               CLOSE ( UnElev ) 

               IF (WaveTimeIn(1) .NE. 0.0) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)') "ERROR MoorDyn WaveElev time series should start at t = 0 seconds."
                  ENDIF
                  CALL SetErrStat( ErrID_Fatal, ' MoorDyn WaveElev time series should start at t = 0 seconds.',ErrStat, ErrMsg, RoutineName); RETURN
               ENDIF
               
               CALL WrScr( "    Read "//trim(num2lstr(ntIn))//" time steps from input file." )
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "    Read "//trim(num2lstr(ntIn))//" time steps from input file."
               ENDIF

               ! IF (WaveTimeIn(ntIn) < TMax) THEN <<<< need to handle if time series is too short?
               
               ! specify stepping details 
               p%ntWave = CEILING(Tmax/p%dtWave)          ! number of wave time steps

               
               ! allocate space for processed reference wave elevation time series
               ALLOCATE ( WaveElev0( 0:p%ntWave ), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate array WaveElev0.'; IF(Failed0()) RETURN  ! this has an extra entry of zero in case it needs to be padded to be even
               WaveElev0 = 0.0_SiKi
               
               ! go through and interpolate (should replace with standard function)
               DO i = 1, p%ntWave       
                  t = p%dtWave*(i-1)
                  
                  ! interpolation routine to dtWave spacing
                  DO iIn = 1,ntIn-1      
                     IF (WaveTimeIn(iIn+1) > t) THEN   ! find the right two points to interpolate between (remember that the first column of PtfmMotIn is time)
                        frac = (t - WaveTimeIn(iIn) )/( WaveTimeIn(iIn+1) - WaveTimeIn(iIn) )  ! interpolation fraction (0-1) between two interpolation points
                        WaveElev0(i-1) = WaveElevIn(iIn) + frac*(WaveElevIn(iIn+1) - WaveElevIn(iIn))  ! get interpolated wave elevation
                        EXIT   ! break out of the loop for this time step once we've DOne its interpolation
                     END IF
                  END DO      
               END DO
               
               ! note: following is adapted from UserWaves.v90 UserWaveElevations_Init
               
               
               
               ! Set new value for NStepWave so that the FFT algorithms are efficient. We will use the values passed in rather than what is read from the file
               
               IF ( MOD(p%ntWave,2) == 1 )  p%ntWave = p%ntWave + 1              ! Set NStepWave to an even integer   
               NStepWave2 = MAX( p%ntWave/2, 1 )                             ! Make sure that NStepWave is an even product of small factors (PSF) that is
               NStepWave  = 2*PSF ( NStepWave2, 9 )                        !   greater or equal to WaveTMax/WaveDT to ensure that the FFT is efficient.
               NStepWave2 = NStepWave/2                                    ! Update the value of NStepWave2 based on the value needed for NStepWave.
               WaveTMax   = NStepWave*p%dtWave                               ! Update the value of WaveTMax   based on the value needed for NStepWave.
               WaveDOmega = TwoPi/TMax                                     ! Compute the frequency step for incident wave calculations.
               p%ntWave = NStepWave
               
               
               

               ! Allocate array to hold the wave elevations for calculation of FFT.
               ALLOCATE ( TmpFFTWaveElev( 0:NStepWave-1), STAT=ErrStat2 ); ErrMsg2 = 'Cannot allocate array TmpFFTWaveElev.'; IF (Failed0()) RETURN

               ! Allocate frequency array for the wave elevation information in frequency space
               ALLOCATE ( WaveElevC0(2, 0:NStepWave2), STAT=ErrStat2 ); ErrMsg2 = 'Cannot allocate array WaveElevC0.'; IF (Failed0()) RETURN
               

               ! Now check if all the allocations worked properly
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
               
               ! Set the values
               TmpFFTWaveElev  =  0.0_DbKi
               WaveElevC0(:,:) =  0.0_DbKi


               ! Copy values over
               DO I=0, MIN(SIZE(WaveElev0), NStepWave)-1
                  TmpFFTWaveElev(I) = WaveElev0(I)
               ENDDO

               ! Initialize the FFT
               CALL InitFFT ( NStepWave, FFT_Data, .FALSE., ErrStat2 ); ErrMsg2 = 'Error occurred while initializing the FFT.'; IF(Failed()) RETURN

               ! Apply the forward FFT on the wave elevation timeseries to get the real and imaginary parts of the frequency information.      
               CALL    ApplyFFT_f (  TmpFFTWaveElev(:), FFT_Data, ErrStat2 ); ErrMsg2 = 'Error occurred while applying the forwards FFT to TmpFFTWaveElev array.'; IF(Failed()) RETURN  ! Note that the TmpFFTWaveElev now contains the real and imaginary bits.

               ! Copy the resulting TmpFFTWaveElev(:) data over to the WaveElevC0 array
               DO I=1,NStepWave2-1
                  WaveElevC0     (1,I) = TmpFFTWaveElev(2*I-1) ! all the odd indicies in array (real part ?)
                  WaveElevC0     (2,I) = TmpFFTWaveElev(2*I) ! all the even indicies in array (imaginary part ?)
               ENDDO
               WaveElevC0(:,NStepWave2) = 0.0_SiKi

               CALL  ExitFFT(FFT_Data, ErrStat2); ErrMsg2 = 'Error occurred while cleaning up after the FFTs.'; IF(Failed()) RETURN

               IF (ALLOCATED( WaveElev0      )) DEALLOCATE( WaveElev0     , STAT=ErrStat2); ErrMsg2 = 'Cannot deallocate WaveElev0.'     ; IF (Failed0()) RETURN  
               IF (ALLOCATED( TmpFFTWaveElev )) DEALLOCATE( TmpFFTWaveElev, STAT=ErrStat2); ErrMsg2 = 'Cannot deallocate TmpFFTWaveElev.'; IF (Failed0()) RETURN

            ENDIF ! End getting frequency data either from time series or WaveField. The below needs to happen for both old and hybrid approach
            
            ! NOTE: there is a decent ammount of code duplication (intentional for now) with what is in SeaState that eventually 
            ! we will want to synchronize with a standard library at some point

            ! allocate all the wave kinematics FFT arrays  
            ALLOCATE( WaveNmbr  (0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveNmbr.'  ; IF (Failed0()) RETURN
            ALLOCATE( tmpComplex(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate tmpComplex.'; IF (Failed0()) RETURN
            ALLOCATE( WaveElevC (0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveElevC.' ; IF (Failed0()) RETURN
            ALLOCATE( WaveDynPC (0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveDynPC.' ; IF (Failed0()) RETURN
            ALLOCATE( WaveVelCHx(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveVelCHx.'; IF (Failed0()) RETURN
            ALLOCATE( WaveVelCHy(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveVelCHy.'; IF (Failed0()) RETURN
            ALLOCATE( WaveVelCV (0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveVelCV.' ; IF (Failed0()) RETURN
            ALLOCATE( WaveAccCHx(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveAccCHx.'; IF (Failed0()) RETURN
            ALLOCATE( WaveAccCHy(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveAccCHy.'; IF (Failed0()) RETURN
            ALLOCATE( WaveAccCV (0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate WaveAccCV .'; IF (Failed0()) RETURN
            
            ! allocate time series grid data arrays (now that we know the number of time steps coming from the IFFTs)
            CALL allocateKinematicsArrays() 
            
            ! allocate the wave direction arrays
            ALLOCATE ( CosWaveDir(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate CosWaveDir.'; IF (Failed0()) RETURN
            ALLOCATE ( SinWaveDir(0:NStepWave2), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate SinWaveDir.'; IF (Failed0()) RETURN       

            ! Set the CosWaveDir and SinWaveDir values. 
            IF (WaveMultiDir) THEN ! This is only possible with WaveKinMod = 2
               CosWaveDir=COS(D2R*WaveDirArr)
               SinWaveDir=SIN(D2R*WaveDirArr)
            ELSE
               CosWaveDir=COS(D2R*WaveDir)
               SinWaveDir=SIN(D2R*WaveDir)
            ENDIF
            
            ! get wave number array once
            DO I = 0, NStepWave2 
               WaveNmbr(i)   = WaveNumber ( REAL(I*WaveDOmega, R8Ki), p%g, p%WtrDpth )
               tmpComplex(I)    =  CMPLX(WaveElevC0(1,I), WaveElevC0(2,I)) ! 1 are real, 2 are imaginary 
            END DO    
            
            ! set up FFTer for doing IFFTs
            CALL InitFFT ( NStepWave, FFT_Data, .TRUE., ErrStat2 ); ErrMsg2 = 'Error occurred while initializing the FFT.'; IF(Failed()) RETURN

            ! Loop through all points where the incident wave kinematics will be computed      
            DO ix = 1,p%nxWave 
               DO iy = 1,p%nyWave
                  DO iz = 1,p%nzWave
                  
                     ! Compute the discrete Fourier transform of the incident wave kinematics
                     DO i = 0, NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms

                        Omega = i*WaveDOmega
                        ImagOmega = ImagNmbr*Omega

                        WaveElevC (i) = tmpComplex(i) * EXP( -ImagNmbr*WaveNmbr(i)*( p%pxWave(ix)*CosWaveDir(i) + p%pyWave(iy)*SinWaveDir(i) ))         ! Discrete Fourier transform of the instantaneous elevation of incident waves at the ref point (meters)                                                             
                        WaveDynPC (i) = p%rhoW*p%g* WaveElevC(i) * COSHNumOvrCOSHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) )                ! Discrete Fourier transform of the instantaneous dynamic pressure of incident waves before applying stretching at the zi-coordinates for points (N/m^2)
                        WaveVelCHx(i) =       Omega*WaveElevC(i) * COSHNumOvrSINHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) ) *CosWaveDir(i) ! Discrete Fourier transform of the instantaneous horizontal velocity     in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
                        WaveVelCHy(i) =       Omega*WaveElevC(i) * COSHNumOvrSINHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) ) *SinWaveDir(i) ! Discrete Fourier transform of the instantaneous horizontal velocity     in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)           
                        WaveVelCV (i) =   ImagOmega*WaveElevC(i) * SINHNumOvrSINHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) )                ! Discrete Fourier transform of the instantaneous vertical   velocity                    of incident waves before applying stretching at the zi-coordinates for points (m/s)
                        WaveAccCHx(i) =   ImagOmega*WaveVelCHx(i)                                                                                       ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
                        WaveAccCHy(i) =   ImagOmega*WaveVelCHy(i)                                                                                       ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
                        WaveAccCV (i) =   ImagOmega*WaveVelCV (i)                                                                                       ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
                     END DO  ! I, frequencies
                     
                     ! now IFFT all the wave kinematics except surface elevation and save it into the grid of data
                     CALL ApplyFFT_cx( p%PDyn  (:,iz,iy,ix), WaveDynPC , FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveDynPC.'; IF (Failed()) RETURN 
                     CALL ApplyFFT_cx( p%uxWave(:,iz,iy,ix), WaveVelCHx, FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveVelHx.'; IF (Failed()) RETURN
                     CALL ApplyFFT_cx( p%uyWave(:,iz,iy,ix), WaveVelCHy, FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveVelHy.'; IF (Failed()) RETURN
                     CALL ApplyFFT_cx( p%uzWave(:,iz,iy,ix), WaveVelCV , FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveVelV.' ; IF (Failed()) RETURN
                     CALL ApplyFFT_cx( p%axWave(:,iz,iy,ix), WaveAccCHx, FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveAccHx.'; IF (Failed()) RETURN
                     CALL ApplyFFT_cx( p%ayWave(:,iz,iy,ix), WaveAccCHy, FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveAccHy.'; IF (Failed()) RETURN
                     CALL ApplyFFT_cx( p%azWave(:,iz,iy,ix), WaveAccCV , FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveAccV.' ; IF (Failed()) RETURN

                  END DO ! iz
                  
                  ! IFFT wave elevation here because it's only at the surface
                  CALL ApplyFFT_cx( p%zeta(:,iy,ix) , WaveElevC , FFT_Data, ErrStat2 ); ErrMsg2 = 'Error IFFTing WaveElev.'; IF (Failed()) RETURN
               END DO ! iy
            END DO ! ix

            ! could also reproduce the wave elevation at 0,0,0 on a separate channel for verIFication...
            
            CALL  ExitFFT(FFT_Data, ErrStat2); ErrMsg2 = 'Error occurred while cleaning up after the IFFTs.'; IF(Failed()) RETURN
      
         END IF ! p%WaveKin > 0


         ! --------------------------------- now do currents --------------------------------
            
         IF (p%Current > 0) THEN
         
            IF (p%Current == 2) THEN

               ! Overwrite the SeaState current grid to use the MoorDyn grid
               p%WaveField%Current_InitInput%WaveKinGridzi = pzCurrentTemp
               p%WaveField%Current_InitInput%NGridPts = p%nzCurrent

               ! Calculate the current profile in the MD grid with SS inputs
               CALL Current_Init(p%WaveField%Current_InitInput, Current_InitOutput, ErrStat2, ErrMsg2); IF(Failed()) RETURN

               ! Check output current arrays are the right size
               IF (SIZE(Current_InitOutput%CurrVxi) /= p%nzCurrent) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   ERROR Current_Init output size does not match the MoorDyn grid size. "//trim(num2lstr(SIZE(Current_InitOutput%CurrVxi)))//" vs "//trim(num2lstr(REAL(p%nzCurrent, SiKi)))
                  ENDIF
                  CALL SetErrStat( ErrID_Fatal,'Current_Init output size does not match the MoorDyn grid size.',ErrStat, ErrMsg, RoutineName); RETURN
               ENDIF

               ! extract the currents from outputs
               DO I = 1, p%nzCurrent
                  uxCurrentTemp(I) = Current_InitOutput%CurrVxi(I)
                  uyCurrentTemp(I) = Current_InitOutput%CurrVyi(I)
               END DO

            ENDIF

            ! allocate current profile arrays to correct size
            CALL AllocAry( p%pzCurrent, p%nzCurrent, 'pzCurrent', ErrStat2, ErrMsg2 ); IF(Failed()) RETURN
            CALL AllocAry( p%uxCurrent, p%nzCurrent, 'uxCurrent', ErrStat2, ErrMsg2 ); IF(Failed()) RETURN
            CALL AllocAry( p%uyCurrent, p%nzCurrent, 'uyCurrent', ErrStat2, ErrMsg2 ); IF(Failed()) RETURN
            
            ! copy over data
            DO i = 1,p%nzCurrent
               p%pzCurrent(i) = pzCurrentTemp(i)
               p%uxCurrent(i) = uxCurrentTemp(i)
               p%uyCurrent(i) = uyCurrentTemp(i)
            END DO

         ENDIF ! p%Current >0

      ENDIF ! (tmpString == "SEASTATE")

      ! some status messages:
      IF (p%WaterKin == 0) THEN
         CALL WrScr("    No water kinematics set up")
         IF (p%writeLog > 0) THEN
            WRITE(p%UnLog, '(A)'        ) "    No water kinematics set up"
         ENDIF
      ELSEIF (p%WaterKin == 1) THEN
         CALL WrScr("    Water kinematics will be simulated using the old method (user provided grid and water kinematic data)")
         IF (p%writeLog > 0) THEN
            WRITE(p%UnLog, '(A)'        ) "    Water kinematics will be simulated using the old method (user provided grid and water kinematic data)"
         ENDIF
      ELSEIF (p%WaterKin == 2) THEN
         CALL WrScr("    Water kinematics will be simulated using the hybrid method (user provided grid with SeaState water kinematics data)")
         IF (p%writeLog > 0) THEN
            WRITE(p%UnLog, '(A)'        ) "    Water kinematics will be simulated using the hybrid method (user provided grid with SeaState water kinematics data)"
         ENDIF
      ELSEIF (p%WaterKin == 3) THEN
         CALL WrScr("    Water kinematics will be simulated using the SeaState method (SeaState provided grid and water kinematics data). These will be disabled during IC solve.")
         IF (p%writeLog > 0) THEN
            WRITE(p%UnLog, '(A)'        ) "    Water kinematics will be simulated using the SeaState method (SeaState provided grid and water kinematics data). These will be disabled during IC solve."
         ENDIF
      ELSE
         CALL SetErrStat( ErrID_Fatal,"Invalid value for WaterKin",ErrStat, ErrMsg, RoutineName); RETURN
      ENDIF

      ! ------------------------------ clean up and finished ---------------------------
      CALL cleanup()
      
      
   CONTAINS
   
   
      ! get grid axis coordinates, initialize/record in array, and return size
      SUBROUTINE gridAxisCoords(coordtype, entries, coordarray, n)
         
         INTEGER(IntKi),          INTENT(IN   )  :: coordtype
         CHARACTER(*),            INTENT(INOUT)  :: entries
         REAL(SiKi), ALLOCATABLE,  INTENT(INOUT)  :: coordarray(:)
         INTEGER(IntKi),          INTENT(  OUT)  :: n
            
         REAL(ReKi) :: tempArray (100)
         REAL(ReKi) :: dx
         INTEGER(IntKi)                   :: nEntries, I

         IF (len(trim(entries)) == len(entries)) THEN
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   WARNING in gridAxisCoords: Only "//trim(num2lstr(len(entries)))//" characters read from wave grid coordinates"
            ENDIF
            CALL SetErrStat(ErrID_Warn, "Only "//trim(num2lstr(len(entries)))//" characters read from wave grid coordinates", ErrStat, ErrMsg, RoutineName)
         END IF

         IF (entries(len(entries):len(entries)) == ',') THEN
            IF (p%writeLog > 0) THEN
               WRITE(p%UnLog, '(A)'        ) "   ERROR in gridAxisCoords: Last character of wave grid coordinate list cannot be comma"
            ENDIF
            CALL SetErrStat(ErrID_Fatal, "Last character of wave grid coordinate list cannot be comma", ErrStat, ErrMsg, RoutineName); RETURN
         ELSE
            ! get array of coordinate entries 
            CALL stringToArray(entries, nEntries, tempArray); IF(Failed()) RETURN 
            
            ! set number of coordinates
            IF (     coordtype==0) THEN   ! 0: not used - make one grid point at zero
               n = 1;
            ELSE IF (coordtype==1) THEN   ! 1: list values in ascending order
               n = nEntries
            ELSE IF (coordtype==2) THEN   ! 2: uniform specified by -xlim, xlim, num
               n = int(tempArray(3))

               ! Check that tmpArray range / n is not an integer value, warn in log file if so
               IF (MOD((tempArray(2)-tempArray(1)),REAL(n,ReKi)) == 0.0_ReKi) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   WARNING in gridAxisCoords: grid range / num grid points is an integer. This may cause interpolation issues. For integer grid spacing, add 1 to num grid points."
                  ENDIF
               ENDIF
               
            ELSE
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "   ERROR in gridAxisCoords: invalid coordinate type specified to gridAxisCoords. Check WaterKin input file format for missing lines"
               ENDIF
               CALL SetErrStat(ErrID_Fatal, "Invalid coordinate type. Check WaterKin input file format for missing lines", ErrStat, ErrMsg, RoutineName); RETURN
            END IF
            
            ! allocate coordinate array
            CALL AllocAry(coordarray, n, 'x,y, or z grid points' , ErrStat, ErrMsg)
            !ALLOCATE ( coordarray(n), STAT=ErrStat) 
            
            ! fill in coordinates
            IF (     coordtype==0) THEN
               coordarray(1) = 0.0_ReKi
            
            ELSE IF (coordtype==1) THEN
               coordarray(1:n) = tempArray(1:n)
            
            ELSE IF (coordtype==2) THEN  
               coordarray(1) = tempArray(1)
               coordarray(n) = tempArray(2)
               dx = (coordarray(n)-coordarray(1))/REAL(n-1)
               DO i=2,n
                  coordarray(i) = coordarray(i-1) + dx
               END DO
            
            ELSE
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "   ERROR in gridAxisCoords: invalid coordinate type specified to gridAxisCoords. Check WaterKin input file format for missing lines"
               ENDIF
               CALL SetErrStat(ErrID_Fatal, "Invalid coordinate type specified to gridAxisCoords. Check WaterKin input file format for missing lines", ErrStat, ErrMsg, RoutineName); RETURN
            END IF
            
            ! print *, "Set water grid coordinates to :"
            ! DO i=1,n
            !    print *, " ", coordarray(i)
            ! END DO
         END IF

      END SUBROUTINE gridAxisCoords
   
   
      ! Extract an array of numbers out of a string with comma-separated numbers (cannot use NWTC function ReadArr for now becasue len not known)
      SUBROUTINE stringToArray(instring, n, outarray)
   
         CHARACTER(*),          INTENT(INOUT)  :: instring
         INTEGER(IntKi),        INTENT(  OUT)  :: n
         REAL(ReKi),            INTENT(  OUT)  :: outarray(100)  ! array of output numbers (100 maximum)

         CHARACTER(40)                         :: tempstring         
         INTEGER                               :: pos1, pos2, i
    
         outarray = 0.0_ReKi
    
         n = 0
         pos1=1
    
         DO
            pos2 = INDEX(instring(pos1:), ",")  ! find index of next comma
            IF (pos2 == 0) THEN                 ! if there isn't another comma, read the last entry and call it done (this could be the only entry if no commas)
               n = n + 1
               READ(instring(pos1:), *, IOSTAT=ErrStat2) outarray(n)
               IF (ErrStat2 /= ErrID_None) THEN
                  IF (p%writeLog > 0) THEN
                     WRITE(p%UnLog, '(A)'        ) "   ERROR in stringToArray: invalid value in string"
                  ENDIF
                  CALL SetErrStat(ErrID_Fatal, "Invalid value in string", ErrStat, ErrMsg, RoutineName); RETURN
               ENDIF
               EXIT
            END IF
            n = n + 1
            IF (n > 100) THEN
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "   ERROR in stringToArray: cannot do more than 100 entries"
               ENDIF
               CALL SetErrStat(ErrID_Fatal, "Cannot do more than 100 entries", ErrStat, ErrMsg, RoutineName); RETURN
            END IF            
            READ(instring(pos1:pos1+pos2-2), *, IOSTAT=ErrStat2) outarray(n)
            IF (ErrStat2 /= ErrID_None) THEN
               IF (p%writeLog > 0) THEN
                  WRITE(p%UnLog, '(A)'        ) "   ERROR in stringToArray: invalid value in string"
               ENDIF
               CALL SetErrStat(ErrID_Fatal, "Invalid value in string", ErrStat, ErrMsg, RoutineName); RETURN
            ENDIF

            pos1 = pos2+pos1
         END DO
         
      END SUBROUTINE stringToArray
   
      
      ! allocate water kinematics arrays
      SUBROUTINE allocateKinematicsArrays()
        !  error check print *, "Error in Waves::makeGrid, a time or space array is size zero." << endl;

        ALLOCATE ( p%uxWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate uxWave'; IF (Failed0()) RETURN
        ALLOCATE ( p%uyWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate uyWave'; IF (Failed0()) RETURN
        ALLOCATE ( p%uzWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate uzWave'; IF (Failed0()) RETURN
        ALLOCATE ( p%axWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate axWave'; IF (Failed0()) RETURN
        ALLOCATE ( p%ayWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate ayWave'; IF (Failed0()) RETURN
        ALLOCATE ( p%azWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate azWave'; IF (Failed0()) RETURN
        ALLOCATE ( p%PDyn  ( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStat2); ErrMsg2 = 'Cannot allocate PDyn'  ; IF (Failed0()) RETURN
        ALLOCATE ( p%zeta  ( p%ntWave,p%nyWave,p%nxWave)         , STAT=ErrStat2); ErrMsg2 = 'Cannot allocate zeta'  ; IF (Failed0()) RETURN    ! 2D grid over x and y only
        
      END SUBROUTINE allocateKinematicsArrays

      
      ! compact way to set the right error status and check if an abort is needed (and do cleanup if so)
      LOGICAL FUNCTION Failed()
           CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetupWaterKin') 
           Failed =  ErrStat >= AbortErrLev
           IF (Failed) CALL CleanUp()
           CALL SetErrStat(ErrStat2, 'Cleanup was not succcessful', ErrStat, ErrMsg, 'SetupWaterKin') ! check that deallocation was successful
      END FUNCTION Failed

      ! check for failed where /= 0 is fatal
      LOGICAL FUNCTION Failed0()
         IF (ErrStat2 /= ErrID_None) THEN
            ErrStat2 = ErrID_Fatal
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         ENDIF
         Failed0 = ErrStat >= AbortErrLev
         IF (Failed0) CALL CleanUp()
         CALL SetErrStat(ErrStat2, 'Cleanup was not succcessful', ErrStat, ErrMsg, 'SetupWaterKin') ! check that deallocation was successful
      END FUNCTION Failed0
   

      SUBROUTINE CleanUp

         !IF (ALLOCATED( WaveElev   ))         DEALLOCATE( WaveElev,          STAT=ErrStat2)
         !IF (ALLOCATED( WaveTime   ))         DEALLOCATE( WaveTime,          STAT=ErrStat2)
         IF (ALLOCATED( TmpFFTWaveElev  ))    DEALLOCATE( TmpFFTWaveElev,    STAT=ErrStat2)
         IF (ALLOCATED( WaveElevC0      ))    DEALLOCATE( WaveElevC0,        STAT=ErrStat2)
         
         ! >>> missing some things <<<
         
         IF (ALLOCATED( WaveNmbr   )) DEALLOCATE( WaveNmbr   , STAT=ErrStat2)
         IF (ALLOCATED( tmpComplex )) DEALLOCATE( tmpComplex , STAT=ErrStat2)
         IF (ALLOCATED( WaveElevC  )) DEALLOCATE( WaveElevC  , STAT=ErrStat2)
         IF (ALLOCATED( WaveDynPC  )) DEALLOCATE( WaveDynPC  , STAT=ErrStat2)
         IF (ALLOCATED( WaveVelCHx )) DEALLOCATE( WaveVelCHx , STAT=ErrStat2)
         IF (ALLOCATED( WaveVelCHy )) DEALLOCATE( WaveVelCHy , STAT=ErrStat2)
         IF (ALLOCATED( WaveVelCV  )) DEALLOCATE( WaveVelCV  , STAT=ErrStat2)
         IF (ALLOCATED( WaveAccCHx )) DEALLOCATE( WaveAccCHx , STAT=ErrStat2)
         IF (ALLOCATED( WaveAccCHy )) DEALLOCATE( WaveAccCHy , STAT=ErrStat2)
         IF (ALLOCATED( WaveAccCV  )) DEALLOCATE( WaveAccCV  , STAT=ErrStat2)

      END SUBROUTINE CleanUp
      
      
            !=======================================================================
      FUNCTION WaveNumber ( Omega, g, h )


         ! This FUNCTION solves the finite depth dispersion relationship:
         !
         !                   k*tanh(k*h)=(Omega^2)/g
         !
         ! for k, the wavenumber (WaveNumber) given the frequency, Omega,
         ! gravitational constant, g, and water depth, h, as inputs.  A
         ! high order initial guess is used in conjunction with a quadratic
         ! Newton's method for the solution with seven significant digits
         ! accuracy using only one iteration pass.  The method is due to
         ! Professor J.N. Newman of M.I.T. as found in routine EIGVAL of
         ! the SWIM-MOTION-LINES (SML) software package in source file
         ! Solve.f of the SWIM module.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(DbKi), INTENT(IN )      :: g                                               ! Gravitational acceleration (m/s^2)
      REAL(DbKi), INTENT(IN )      :: h                                               ! Water depth (meters)
      REAL(DbKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(DbKi)                   :: WaveNumber                                      ! This function = wavenumber, k (1/m)


         ! Local Variables:

      REAL(DbKi)                   :: A                                               ! A temporary variable used in the solution.
      REAL(DbKi)                   :: B                                               ! A temporary variable used in the solution.
      REAL(DbKi)                   :: C                                               ! A temporary variable used in the solution.
      REAL(DbKi)                   :: C2                                              ! A temporary variable used in the solution.
      REAL(DbKi)                   :: CC                                              ! A temporary variable used in the solution.
      REAL(DbKi)                   :: E2                                              ! A temporary variable used in the solution.
      REAL(DbKi)                   :: X0                                              ! A temporary variable used in the solution.



         ! Compute the wavenumber, unless Omega is zero, in which case, return
         !   zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         WaveNumber = 0.0


      ELSE                       ! Omega > 0.0; solve for the wavenumber as usual.


         C  = Omega*Omega*h/REAL(g,DbKi)
         CC = C*C


         ! Find X0:

         IF ( C <= 2.0 )  THEN

            X0 = SQRT(C)*( 1.0 + C*( 0.169 + (0.031*C) ) )

         ELSE

            E2 = EXP(-2.0*C)

            X0 = C*( 1.0 + ( E2*( 2.0 - (12.0*E2) ) ) )

         END IF


         ! Find the WaveNumber:

         IF ( C <= 4.8 )  THEN

            C2 = CC - X0*X0
            A  = 1.0/( C - C2 )
            B  = A*( ( 0.5*LOG( ( X0 + C )/( X0 - C ) ) ) - X0 )

            WaveNumber = ( X0 - ( B*C2*( 1.0 + (A*B*C*X0) ) ) )/h

         ELSE

            WaveNumber = X0/h

         END IF


      END IF



      RETURN
      END FUNCTION WaveNumber
      
       !=======================================================================
      FUNCTION COSHNumOvrCOSHDen ( k, h, z )

      
         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/COSH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.

      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi)                   :: COSHNumOvrCOSHDen                               ! This function = COSH( k*( z + h ) )/COSH( k*h ) (-)
      REAL(DbKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(DbKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(DbKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF ( k*h  > 89.4_DbKi )  THEN   ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/COSH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrCOSHDen = REAL(EXP( k*z ) + EXP( -k*( z + 2.0_DbKi*h ) ))

      ELSE                       ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrCOSHDen =REAL( COSH( k*( z + h ) ),R8Ki)/COSH( k*h )

      END IF



      RETURN
      END FUNCTION COSHNumOvrCOSHDen
!=======================================================================
      FUNCTION COSHNumOvrSINHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi)                   :: COSHNumOvrSINHDen                               ! This function = COSH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(DbKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(DbKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(DbKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:


      IF (   k  < EPSILON(0.0_DbKi)  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, HUGE(k) is returned to approximate the known value of infinity.

         COSHNumOvrSINHDen = 1.0E20   ! HUGE( k )

      ELSEIF ( k*h  > 89.4_DbKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrSINHDen = EXP( k*z ) + EXP( -k*( z + 2*h ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrSINHDen = COSH( k*( z + h ) )/SINH( k*h )

      END IF



      RETURN
      END FUNCTION COSHNumOvrSINHDen
!=======================================================================
      FUNCTION COTH ( X )


         ! This FUNCTION computes the hyperbolic cotangent,
         ! COSH(X)/SINH(X).


      USE                             Precision


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(DbKi)                   :: COTH                                            ! This function = COSH( X )/SINH( X ) (-)
      REAL(DbKi), INTENT(IN )      :: X                                               ! The argument (-)



         ! Compute the hyperbolic cotangent:

      IF ( X == 0.0_DbKi )  THEN   ! When .TRUE., the formulation below is ill-conditioned; thus, HUGE(X) is returned to approximate the known value of infinity.

         COTH = HUGE( X )

      ELSE                    ! X /= 0.0; use the numerically-stable computation of COTH(X) by means of TANH(X).

         COTH = 1.0_DbKi/TANH( X ) ! = COSH( X )/SINH( X )

      END IF



      RETURN
      END FUNCTION COTH
    
      !=======================================================================
      FUNCTION SINHNumOvrSINHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    SINH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi)                   :: SINHNumOvrSINHDen                               ! This function = SINH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(DbKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(DbKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(DbKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF (     k   == 0.0_DbKi  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.

         SINHNumOvrSINHDen = 1.0

      ELSEIF ( k*h >  89.4_DbKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, SINH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) - EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         SINHNumOvrSINHDen = EXP( k*z ) - EXP( -k*( z + 2.0_DbKi*h ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         SINHNumOvrSINHDen = SINH( k*( z + h ) )/SINH( k*h )

      END IF



      RETURN
      END FUNCTION SINHNumOvrSINHDen

   END SUBROUTINE setupWaterKin

   ! Unused water kinematics file writing subroutines:
        
   ! ----- write wave grid spacing to output file -----
   SUBROUTINE WriteWaveGrid(p, ErrStat, ErrMsg)
   
      TYPE(MD_ParameterType),       INTENT(INOUT)  :: p     ! Parameters
      
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat            ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg             ! Error message if ErrStat /= ErrID_None
  
      INTEGER(IntKi)                   :: ErrStat2
  !    CHARACTER(ErrMsgLen)             :: ErrMsg2   
   
      CHARACTER(120)                   :: Frmt       
      INTEGER(IntKi)   :: UnOut    ! for outputing wave kinematics data
      INTEGER(IntKi)   :: I
          
       
       CALL GetNewUnit( UnOut)
    
       CALL OpenFOutFile ( UnOut, "waves.txt", ErrStat, ErrMsg )
       IF ( ErrStat > ErrID_None ) THEN
          ErrMsg = ' Error opening wave grid file: '//TRIM(ErrMsg)
          ErrStat = ErrID_Fatal
          RETURN
       END IF
    
       WRITE(UnOut, *, IOSTAT=ErrStat2)  TRIM( 'MoorDyn v2 wave/current kinematics grid file' )
       WRITE(UnOut, *, IOSTAT=ErrStat2)  TRIM( '---------------------------------------------' )
       WRITE(UnOut, *, IOSTAT=ErrStat2)  TRIM( 'The following 6 lines (4-9) specify the input type then the inputs for x, then, y, then z coordinates.' )
       
       WRITE(UnOut,*, IOSTAT=ErrStat2)  TRIM( '1  - X input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)' )
       Frmt = '('//TRIM(Int2LStr(5))//'(A1,e10.4))'      
       WRITE(UnOut,*, IOSTAT=ErrStat2)  ( " ", TRIM(Num2LStr(p%pxWave(I))), I=1,p%nxWave )
       
       WRITE(UnOut,*, IOSTAT=ErrStat2)  TRIM( '1  - Y input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)' )
       Frmt = '('//TRIM(Int2LStr(5))//'(A1,e10.4))'      
       WRITE(UnOut,*, IOSTAT=ErrStat2)  ( " ", TRIM(Num2LStr(p%pyWave(I))), I=1,p%nyWave )
       
       WRITE(UnOut,*, IOSTAT=ErrStat2)  TRIM( '1  - Z input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)' )
       Frmt = '('//TRIM(Int2LStr(8))//'(A1,e10.4))'      
       WRITE(UnOut,*, IOSTAT=ErrStat2)  ( " ", TRIM(Num2LStr(p%pzWave(I))), I=1,p%nzWave )
       
       CLOSE(UnOut, IOSTAT = ErrStat2 )
       IF ( ErrStat2 /= 0 ) THEN
          ErrStat = ErrID_Severe
          ErrMsg = 'Error closing wave grid file'
       END IF
       
     END SUBROUTINE WriteWaveGrid       
       
       
     ! ----- write wave kinematics grid data to output file -----
     SUBROUTINE WriteWaveData(p, ErrStat, ErrMsg)
     
      TYPE(MD_ParameterType),       INTENT(INOUT)  :: p     ! Parameters
      
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat            ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg             ! Error message if ErrStat /= ErrID_None
  
      INTEGER(IntKi)                   :: ErrStat2
  !    CHARACTER(ErrMsgLen)             :: ErrMsg2   
      
      INTEGER(IntKi)   :: UnOut    ! for outputing wave kinematics data
      INTEGER(IntKi)   :: I,J,K, l, Itemp
      
      CALL GetNewUnit( UnOut)
    
       CALL OpenFOutFile ( UnOut, "wave data.txt", ErrStat, ErrMsg )
       IF ( ErrStat > ErrID_None ) THEN
          ErrMsg = ' Error opening wave grid file: '//TRIM(ErrMsg)
          ErrStat = ErrID_Fatal
          RETURN
       END IF
       
       ! write channel labels
       
       
       ! time
       WRITE(UnOut,"(A10)", IOSTAT=ErrStat2, advance="no") "Time"
    
       DO J = 1,p%nyWave     !y
          DO K = 1,p%nxWave  !x 
             WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ze0", Num2Lstr(J+10*K)
          END DO
       END DO
       DO I=1,p%nzWave          !z
          DO J = 1,p%nyWave     !y
             DO K = 1,p%nxWave  !x 
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ux", Num2Lstr(I+10*J+100*K)
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " uy", Num2Lstr(I+10*J+100*K)
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " uz", Num2Lstr(I+10*J+100*K)
             END DO
          END DO
       END DO
       DO I=1,p%nzWave          !z
          DO J = 1,p%nyWave     !y
             DO K = 1,p%nxWave  !x 
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ax", Num2Lstr(I+10*J+100*K)
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ay", Num2Lstr(I+10*J+100*K)
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " az", Num2Lstr(I+10*J+100*K)
             END DO
          END DO
       END DO
       DO I=1,p%nzWave          !z
          DO J = 1,p%nyWave     !y
             DO K = 1,p%nxWave  !x 
                WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " pd", Num2Lstr(I+10*J+100*K)
             END DO
          END DO
       END DO
    
       ! end the line
       WRITE(UnOut, "(A1)", IOSTAT=ErrStat2) " "
       
       
       
       DO l=1, p%ntWave  ! loop through all time steps
       
          ! time
          WRITE(UnOut,"(F10.4)", IOSTAT=ErrStat2, advance="no") p%dtWave*(l-1)
          !WRITE(UnOut,"(F10.4)", IOSTAT=ErrStat2, advance="no") InitInp%WaveTime(l)
       
          ! wave elevation (all slices for now, to check)
          DO J = 1,p%nyWave     !y
             DO K = 1,p%nxWave  !x 
                Itemp = (J-1)*p%nxWave + K    ! index of actual node
                
                WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%zeta(l,J,K)
             END DO
          END DO
       
          ! wave velocities
          DO I=1,p%nzWave          !z
             DO J = 1,p%nyWave     !y
                DO K = 1,p%nxWave  !x 
                   Itemp = (I-1)*p%nxWave*p%nyWave + (J-1)*p%nxWave + K    ! index of actual node
                   
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%uxWave(l,I,J,K)
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%uyWave(l,I,J,K)
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%uzWave(l,I,J,K)
                END DO
             END DO
          END DO
          
          ! wave accelerations
          DO I=1,p%nzWave          !z
             DO J = 1,p%nyWave     !y
                DO K = 1,p%nxWave  !x 
                   Itemp = (I-1)*p%nxWave*p%nyWave + (J-1)*p%nxWave + K    ! index of actual node
                   
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%axWave(l,I,J,K)
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%ayWave(l,I,J,K)
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%azWave(l,I,J,K)
                END DO
             END DO
          END DO
       
          ! dynamic pressure
          DO I=1,p%nzWave          !z
             DO J = 1,p%nyWave     !y
                DO K = 1,p%nxWave  !x 
                   Itemp = (I-1)*p%nxWave*p%nyWave + (J-1)*p%nxWave + K    ! index of actual node
                   
                   WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%PDyn(l,I,J,K)
                END DO
             END DO
          END DO
       
          ! end the line
          WRITE(UnOut, "(A1)", IOSTAT=ErrStat2) " "
       
       
       END DO
       
       
       CLOSE(UnOut, IOSTAT = ErrStat )
       IF ( ErrStat /= 0 ) THEN
          ErrMsg = 'Error closing wave grid file'
       END IF
    
     END SUBROUTINE WriteWaveData    

END MODULE MoorDyn_Misc
