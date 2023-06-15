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

      if ( .NOT. EqualRealNos(length, 0.0_DbKi ) ) THEN
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
      
      if (length_squared > 0) then
         scaler = newlength/sqrt(length_squared)
      else                   ! if original vector is zero, return zero
         scaler = 0.0_DbKi
      end if
      
      DO J=1,3
         u_out(J) = u_in(J) * scaler 
      END DO

   END SUBROUTINE ScaleVector
   !-----------------------------------------------------------------------


   ! convenience function to calculate curvature based on adjacent segments' direction vectors and their combined length
   function GetCurvature(length, q1, q2)
   
      real(DbKi),   intent(in   ) :: length
      real(DbKi),   intent(in   ) :: q1(3)
      real(DbKi),   intent(in   ) :: q2(3)
      real(DbKi)                  :: GetCurvature
      
      
      real(DbKi)                  :: q1_dot_q2
      
      ! note "length" here is combined from both segments
      
      q1_dot_q2 = dot_product( q1, q2 )
      
      if (q1_dot_q2 > 1.0) then           ! this is just a small numerical error, so set q1_dot_q2 to 1
         GetCurvature = 0.0_DbKi          ! this occurs when there's no curvature, so return zero curvature
         
      !else if (q1_dot_q2 < 0)   ! this is a bend of more than 90 degrees, too much, call an error!

      else                                                        ! normal case
         GetCurvature = 4.0/length * sqrt(0.5*(1.0 - q1_dot_q2))    ! this is the normal curvature calculation
      end if
      
      return
   end function GetCurvature
   

   ! calculate orientation angles of a direction vector
   !-----------------------------------------------------------------------
   subroutine GetOrientationAngles(vec, phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat)
      real(DbKi),   intent(in   ) :: vec(3) !p1(3),p2(3)
      real(DbKi),   intent(  out) :: phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat(3)
            
      real(DbKi)                  ::  vecLen, vecLen2D
 
      vecLen   = SQRT(Dot_Product(vec,vec))
      vecLen2D = SQRT(vec(1)**2+vec(2)**2)
      if ( vecLen < 0.000001 ) then
         print *, "ERROR in GetOrientationAngles in MoorDyn. Supplied vector is near zero" 
         print *, vec
         k_hat = NaN ! 1.0/0.0
      else
         k_hat = vec / vecLen 
         phi   = atan2(vecLen2D, vec(3))  ! incline angle   
      end if
      if ( phi < 0.000001) then
         beta = 0.0_ReKi
      else
         beta = atan2(vec(2), vec(1))                    ! heading of incline     
      endif
      sinPhi  = sin(phi)
      cosPhi  = cos(phi)  
      tanPhi  = tan(phi)     
      sinBeta = sin(beta)
      cosBeta = cos(beta)
            
   end subroutine GetOrientationAngles
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
      !       already be correct or can be assigned seperately from r_in and rd_in (assuming orientation frames are identical)


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
      !       already be correct or can be assigned seperately from r_in and rd_in (assuming orientation frames are identical)


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
            
            if ( EqualRealNos( LU(k,k), 0.0_DbKi) ) then
               LU(k,j) = 0.0_DbKi   ! avoid divide by zero <<< numerator likely zero too. check if this is safe <<<
            else
               LU(k,j) = (A(k,j)-sum)/LU(k,k)
            end if
         END DO !j
         
      END DO !K
      
      DO i=1,n
      
         sum = 0.0_DbKi
         
         DO k=1,i-1  !for(int k=0; k<i; ++k)
            sum = sum + LU(i,k)*y(k);
         END DO
         
         if ( EqualRealNos( LU(i,i), 0.0_DbKi) ) then
            y(i) = 0.0_DbKi   ! avoid divide by zero <<< numerator likely zero too. check if this is safe <<<
         else
            y(i) = (b(i)-sum)/LU(i,i)
         end if
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
      
      if (xin <= xlist(1)) THEN                !  below lowest data point
         i = 1_IntKi
         fout = 0.0_DbKi
      
      else if (xlist(nx) <= xin) THEN          ! above highest data point
         i = nx
         fout = 0.0_DbKi
      
      else                                     ! within the data range
     
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
      
      if (xin <= xlist(1)) THEN                !  below lowest data point
         i = 1_IntKi
         fout = 0.0_SiKi
      
      else if (xlist(nx) <= xin) THEN          ! above highest data point
         i = nx
         fout = 0.0_SiKi
      
      else                                     ! within the data range
     
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
      if (fx == 0) then 
         ix1 = ix0
      else  
         ix1 = min(ix0+1,size(f,4))    ! don't overstep bounds
      end if
      
      if (fy == 0) then
         iy1 = iy0
      else
         iy1 = min(iy0+1,size(f,3))    ! don't overstep bounds
      end if
      
      if (fz == 0) then
         iz1 = iz0
      else         
         iz1 = min(iz0+1,size(f,2))    ! don't overstep bounds
      end if
      
      if (ft == 0) then
         it1 = it0
      else  
         it1 = min(it0+1,size(f,1))    ! don't overstep bounds
      end if
      
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
      if (fx == 0) then 
         ix1 = ix0
      else  
         ix1 = min(ix0+1,size(f,3))    ! don't overstep bounds
      end if
      
      if (fy == 0) then
         iy1 = iy0
      else
         iy1 = min(iy0+1,size(f,2))    ! don't overstep bounds
      end if
      
      if (fz == 0) then
         iz1 = iz0
      else         
         iz1 = min(iz0+1,size(f,1))    ! don't overstep bounds
      end if
      
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
      if ( dx > 0.0 ) then
         dc_dx = (c1y-c0y)/dx
      else
         dc_dx = 0.0_DbKi   ! maybe this should raise an error
      end if
      if ( dx > 0.0 ) then
         dc_dy = (cx1-cx0)/dy
      else
         dc_dy = 0.0_DbKi   ! maybe this should raise an error
      end if
      
      tempVector(1) = dc_dx
      tempVector(2) = dc_dy
      tempVector(3) = 1.0_DbKi
      CALL ScaleVector( tempVector, 1.0_DbKi, nvec ) ! compute unit vector      

   END SUBROUTINE getDepthFromBathymetry
   
   
   ! :::::::::::::::::::::::::: wave and current subroutines :::::::::::::::::::::::::::::::
   
   
   ! master function to get wave/water kinematics at a given point -- called by each object from grid-based data
   SUBROUTINE getWaterKin(p, x, y, z, t, tindex, U, Ud, zeta, PDyn)
   
      ! This whole approach assuems that px, py, and pz are in increasing order.
      ! Wheeler stretching is now built in.
   
      TYPE(MD_ParameterType),INTENT (IN   )       :: p       ! MoorDyn parameters (contains the wave info for now)
      Real(DbKi),            INTENT (IN   )       :: x
      Real(DbKi),            INTENT (IN   )       :: y
      Real(DbKi),            INTENT (IN   )       :: z
      Real(DbKi),            INTENT (IN   )       :: t
      INTEGER(IntKi),        INTENT (INOUT)       :: tindex  ! pass time index to try starting from, returns identified time index
      Real(DbKi),            INTENT (INOUT)       :: U(3)
      Real(DbKi),            INTENT (INOUT)       :: Ud(3)
      Real(DbKi),            INTENT (INOUT)       :: zeta
      Real(DbKi),            INTENT (INOUT)       :: PDyn


      INTEGER(IntKi)             :: ix, iy, iz, it        ! indices for interpolation      
      INTEGER(IntKi)             :: iz0, iz1              ! special indices for currrent interpolation  
!      INTEGER(IntKi)             :: N                     ! number of rod elements for convenience
      Real(SiKi)                 :: fx, fy, fz, ft        ! interpolation fractions
      Real(DbKi)                 :: zp                    ! zprime coordinate used for Wheeler stretching
   
   
      ! if wave kinematics enabled, get interpolated values from grid
      if (p%WaveKin > 0) then      
      
         ! find time interpolation indices and coefficients
         !CALL getInterpNums(p%tWave, t, tindex, it, ft)
         it = floor(t/ p%dtWave) + 1    ! add 1 because Fortran indexing starts at 1
         ft = (t - (it-1)*p%dtWave)/p%dtWave
         tindex = it                  
      
         ! find x-y interpolation indices and coefficients
         CALL getInterpNumsSiKi(p%pxWave   , REAL(x,SiKi),  1, ix, fx)
         CALL getInterpNumsSiKi(p%pyWave   , REAL(y,SiKi),  1, iy, fy)
         
         ! interpolate wave elevation
         CALL calculate3Dinterpolation(p%zeta, ix, iy, it, fx, fy, ft, zeta)   
         
         ! compute modified z coordinate to be used for interpolating velocities and accelerations with Wheeler stretching
         zp = ( z - zeta ) * p%WtrDpth/( p%WtrDpth + zeta )
         
         CALL getInterpNumsSiKi(p%pzWave   , REAL(zp,SiKi),  1, iz, fz)

         ! interpolate everything else
         CALL calculate4Dinterpolation(p%PDyn  , ix, iy, iz, it, fx, fy, fz, ft, PDyn)         
         CALL calculate4Dinterpolation(p%uxWave, ix, iy, iz, it, fx, fy, fz, ft, U(1)  )
         CALL calculate4Dinterpolation(p%uyWave, ix, iy, iz, it, fx, fy, fz, ft, U(2)  )
         CALL calculate4Dinterpolation(p%uzWave, ix, iy, iz, it, fx, fy, fz, ft, U(3)  )      
         CALL calculate4Dinterpolation(p%axWave, ix, iy, iz, it, fx, fy, fz, ft, Ud(1) )
         CALL calculate4Dinterpolation(p%ayWave, ix, iy, iz, it, fx, fy, fz, ft, Ud(2) )
         CALL calculate4Dinterpolation(p%azWave, ix, iy, iz, it, fx, fy, fz, ft, Ud(3) )
      else
         U  = 0.0_DbKi
         Ud = 0.0_DbKi
         zeta = 0.0_DbKi
         PDyn = 0.0_DbKi
      end if
      
      
      ! if current kinematics enabled, add interpolated current values from profile
      if (p%Current > 0) then      
      
         CALL getInterpNumsSiKi(p%pzCurrent, REAL(z,SiKi), 1, iz0, fz)
                  
         IF (fz == 0) THEN  ! handle end case conditions
            iz1 = iz0
         ELSE
            iz1 = min(iz0+1,size(p%pzCurrent))  ! don't overstep bounds
         END IF
         
         U(1) = U(1) + (1.0-fz)*p%uxCurrent(iz0) + fz*p%uxCurrent(iz1)
         U(2) = U(2) + (1.0-fz)*p%uyCurrent(iz0) + fz*p%uyCurrent(iz1)
      end if
      
   END SUBROUTINE getWaterKin


  !! ! unused routine with old code for taking wave kinematic grid inputs from HydroDyn
  !! SUBROUTINE CopyWaterKinFromHydroDyn(p, InitInp)
  !!
  !!    TYPE(MD_InitInputType),       INTENT(IN   )  :: InitInp     ! INTENT(INOUT) : Input data for initialization routine
  !!    TYPE(MD_ParameterType),       INTENT(  OUT)  :: p           ! INTENT( OUT) : Parameters
  !!    
  !!    INTEGER(IntKi)                               :: I, J, K, Itemp
  !!
  !!
  !!    ! ----------------------------- Arrays for wave kinematics -----------------------------
  !!    
  !!    
  !!!   :::::::::::::: BELOW WILL BE USED EVENTUALLY WHEN WAVE INFO IS AN INPUT ::::::::::::::::::
  !!!   ! The rAll array contains all nodes or reference points in the system 
  !!!   ! (x,y,z global coordinates for each) in the order of bodies, rods, points, internal line nodes.      
  !!!   
  !!!   ! count the number of nodes to use for passing wave kinematics
  !!!   J=0 
  !!!   ! Body reference point coordinates
  !!!   J = J + p%nBodies
  !!!   ! Rod node coordinates (including ends)
  !!!   DO l = 1, p%nRods
  !!!      J = J + (m%RodList(l)%N + 1)
  !!!   END DO
  !!!   ! Point reference point coordinates
  !!!   J = J + p%nConnects
  !!!   ! Line internal node coordinates
  !!!   DO l = 1, p%nLines
  !!!      J = J + (m%LineList(l)%N - 1)
  !!!   END DO
  !!!
  !!!   ! allocate all relevant arrays
  !!!   ! allocate state vector and temporary state vectors based on size just calculated
  !!!   ALLOCATE ( y%rAll(3,J), u%U(3,J), u%Ud(3,J), u%zeta(J), u%PDyn(J), STAT = ErrStat )
  !!!   IF ( ErrStat /= ErrID_None ) THEN
  !!!     ErrMsg  = ' Error allocating wave kinematics vectors.'
  !!!     RETURN
  !!!   END IF
  !!!
  !!!
  !!!   ! go through the nodes and fill in the data (this should maybe be turned into a global function)
  !!!   J=0 
  !!!   ! Body reference point coordinates
  !!!   DO I = 1, p%nBodies
  !!!      J = J + 1                     
  !!!      y%rAll(:,J) = m%BodyList(I)%r6(1:3)         
  !!!   END DO
  !!!   ! Rod node coordinates
  !!!   DO I = 1, p%nRods
  !!!      DO K = 0,m%RodList(I)%N  
  !!!         J = J + 1             
  !!!         y%rAll(:,J) = m%RodList(I)%r(:,K)
  !!!      END DO
  !!!   END DO
  !!!   ! Point reference point coordinates
  !!!   DO I = 1, p%nConnects
  !!!      J = J + 1
  !!!      y%rAll(:,J) = m%ConnectList(I)%r
  !!!   END DO      
  !!!   ! Line internal node coordinates
  !!!   DO I = 1, p%nLines
  !!!      DO K = 1,m%LineList(I)%N-1
  !!!         J = J + 1               
  !!!         y%rAll(:,J) = m%LineList(I)%r(:,K)
  !!!      END DO
  !!!   END DO        
  !!   ! :::::::::::::::: the above might be used eventually. For now, let's store wave info grids within this module :::::::::::::::::
  !!   
  !!   
  !!   ! ----- copy wave grid data over from HydroDyn (as was done in USFLOWT branch) -----
  !!   
  !!   ! get grid and time info (currently this is hard-coded to match what's in HydroDyn_Input
  !!   ! DO I=1,p%nzWave
  !!   !    p%pz(I) =  1.0 - 2.0**(p%nzWave-I)       !  -127,  -63,  -31,  -15,   -7,   -3,   -1,    0
  !!   ! END DO
  !!   ! DO J = 1,p%nyWave
  !!   !    p%py(J) = WaveGrid_y0 + WaveGrid_dy*(J-1)
  !!   ! END DO
  !!   ! DO K = 1,p%nxWave
  !!   !    p%px(K) = WaveGrid_x0 + WaveGrid_dx*(K-1)
  !!   ! END DO
  !!   ! 
  !!   ! p%tWave = InitInp%WaveTime
  !!   
  !!   DO I=1,p%nzWave
  !!      DO J = 1,p%nyWave
  !!         DO K = 1,p%nxWave 
  !!            Itemp = (I-1)*p%nxWave*p%nyWave + (J-1)*p%nxWave + K    ! index of actual node on 3D grid
  !!            
  !!            p%uxWave (:,I,J,K) = InitInp%WaveVel( :,Itemp,1)  ! note: indices are t, z, y, x
  !!            p%uyWave (:,I,J,K) = InitInp%WaveVel( :,Itemp,2)
  !!            p%uzWave (:,I,J,K) = InitInp%WaveVel( :,Itemp,3)
  !!            p%axWave (:,I,J,K) = InitInp%WaveAcc( :,Itemp,1)
  !!            p%ayWave (:,I,J,K) = InitInp%WaveAcc( :,Itemp,2)
  !!            p%azWave (:,I,J,K) = InitInp%WaveAcc( :,Itemp,3)
  !!            p%PDyn(   :,I,J,K) = InitInp%WavePDyn(:,Itemp)
  !!         END DO
  !!      END DO
  !!   END DO
  !!      
  !!   DO J = 1,p%nyWave
  !!      DO K = 1,p%nxWave
  !!         Itemp = (J-1)*p%nxWave + K    ! index of actual node on surface 2D grid   
  !!         p%zeta(:,J,K) = InitInp%WaveElev(:,Itemp)
  !!      END DO
  !!   END DO
  !!   
  !! END SUBROUTINE CopyWaterKinFromHydroDyn
   
     
   ! ----- write wave grid spacing to output file -----
   SUBROUTINE WriteWaveGrid(p, ErrStat, ErrMsg)
   
    TYPE(MD_ParameterType),       INTENT(INOUT)  :: p     ! Parameters
    
    INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat            ! Error status of the operation
    CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg             ! Error message if ErrStat /= ErrID_None

    INTEGER(IntKi)                   :: ErrStat2
!    CHARACTER(120)                   :: ErrMsg2   
 
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
!    CHARACTER(120)                   :: ErrMsg2   
    
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


   ! ----- process WaterKin input value, potentially reading wave inputs and generating wave field -----
   SUBROUTINE setupWaterKin(WaterKinString, p, Tmax, ErrStat, ErrMsg)

      CHARACTER(40),           INTENT(IN   )  :: WaterKinString      ! string describing water kinematics filename
      TYPE(MD_ParameterType),  INTENT(INOUT)  :: p                   ! Parameters
      REAL(ReKi),              INTENT(IN   )  :: Tmax
      INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat             ! Error status of the operation
      CHARACTER(*),            INTENT(  OUT)  :: ErrMsg              ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                   :: I, iIn, ix, iy, iz
      INTEGER(IntKi)                   :: ntIn   ! number of time series inputs from file      
      INTEGER(IntKi)                   :: UnIn   ! unit number for coefficient input file
      INTEGER(IntKi)                   :: UnEcho       
      REAL(SiKi)                       :: pzCurrentTemp(100)   ! current depth increments read in from input file (positive-down at this stage)
      REAL(SiKi)                       :: uxCurrentTemp(100)
      REAL(SiKi)                       :: uyCurrentTemp(100)
      
      CHARACTER(120)                   :: WaveKinFile   
      INTEGER(IntKi)                   :: UnElev  ! unit number for coefficient input file
      REAL(SiKi), ALLOCATABLE          :: WaveTimeIn(:)  ! temporarily holds wave input time series
      REAL(SiKi), ALLOCATABLE          :: WaveElevIn(:)
      REAL(SiKi), ALLOCATABLE          :: WaveElev0(:)   ! interpolated reference wave elevation time series
      REAL(SiKi)                       :: WaveDir
      REAL(SiKi)                       :: t, Frac
      CHARACTER(1024)                  :: FileName             ! Name of MoorDyn input file  
      CHARACTER(120)                   :: Line
!      CHARACTER(120)                   :: Line2  
      CHARACTER(120)                   :: entries2  
      INTEGER(IntKi)                   :: coordtype
   
      INTEGER(IntKi)                   :: NStepWave    ! 
      INTEGER(IntKi)                   :: NStepWave2   ! 
      REAL(SiKi)                       :: WaveTMax     ! max wave elevation time series duration after optimizing lenght for FFT
      REAL(SiKi)                       :: WaveDOmega   
      REAL(SiKi)                       :: SinWaveDir                                      ! SIN( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction.
      REAL(SiKi)                       :: CosWaveDir                                      ! COS( WaveDirArr(I) ) -- Each wave frequency has a unique wave direction.

      REAL(SiKi),  ALLOCATABLE         :: TmpFFTWaveElev(:)     ! Data for the FFT calculation
      TYPE(FFT_DataType)               :: FFT_Data              ! the instance of the FFT module we're using
      
      
      COMPLEX(SiKi),ALLOCATABLE        :: tmpComplex(:)      ! A temporary array (0:NStepWave2-1) for FFT use. 
   
      REAL(SiKi)                       :: Omega                 ! Wave frequency (rad/s)
      COMPLEX(SiKi), PARAMETER         :: ImagNmbr = (0.0,1.0)  ! The imaginary number, SQRT(-1.0)
      COMPLEX(SiKi)                    :: ImagOmega             ! = ImagNmbr*Omega (rad/s)
      REAL(DbKi), ALLOCATABLE          :: WaveNmbr(:)           ! wave number for frequency array
      REAL(SiKi), ALLOCATABLE          :: WaveElevC0(:,:)        ! Discrete Fourier transform of the instantaneous elevation of incident waves at the ref point (meters)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveElevC( :)        ! Discrete Fourier transform of the instantaneous elevation of incident waves at the ref point (meters)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveAccCHx(:)       ! Discrete Fourier transform of the instantaneous horizontal acceleration in x-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveAccCHy(:)       ! Discrete Fourier transform of the instantaneous horizontal acceleration in y-direction of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveAccCV( :)        ! Discrete Fourier transform of the instantaneous vertical   acceleration                of incident waves before applying stretching at the zi-coordinates for points (m/s^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveDynPC( :)        ! Discrete Fourier transform of the instantaneous dynamic pressure                       of incident waves before applying stretching at the zi-coordinates for points (N/m^2)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveVelCHx(:)       ! Discrete Fourier transform of the instantaneous horizontal velocity                    of incident waves before applying stretching at the zi-coordinates for points (m/s)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveVelCHy(:)       ! Discrete Fourier transform of the instantaneous horizontal velocity in x-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)
      COMPLEX(SiKi), ALLOCATABLE       :: WaveVelCV( :)        ! Discrete Fourier transform of the instantaneous vertical   velocity in y-direction     of incident waves before applying stretching at the zi-coordinates for points (m/s)
!      COMPLEX(SiKi)                    :: WGNC                  ! Discrete Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component (-)

      INTEGER(IntKi)                   :: ErrStatTmp
      INTEGER(IntKi)                   :: ErrStat2
      CHARACTER(120)                   :: ErrMsg2   
      CHARACTER(120)                   :: RoutineName = 'SetupWaveKin'   


      ErrStatTmp = ErrID_None  ! TODO: get rid of redundancy <<<
      ErrStat2 = ErrID_None
      ErrMsg2  = ""

      IF (LEN_TRIM(WaterKinString) == 0) THEN
         ! If the input is empty (not provided), there are no water kinematics to be included
         p%WaveKin = 0
         p%Current = 0
         return
         
      ELSE IF (SCAN(WaterKinString, "abcdfghijklmnopqrstuvwxyzABCDFGHIJKLMNOPQRSTUVWXYZ") == 0) THEN
         ! If the input has no letters, let's assume it's a number         
         print *, "ERROR WaveKin option does not currently support numeric entries. It must be a filename."
         p%WaveKin = 0
         p%Current = 0
         return
      END IF


      ! otherwise interpret the input as a file name to load the bathymetry lookup data from
      print *, "   The waterKin input contains letters so will load a water kinematics input file"
      
      
      ! -------- load water kinematics input file -------------
      
      IF ( PathIsRelative( WaterKinString ) ) THEN   ! properly handle relative path <<<
         !CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
         FileName = TRIM(p%PriPath)//TRIM(WaterKinString)
      ELSE
         FileName = trim(WaterKinString)
      END IF
      
      
      
      UnEcho=-1
      CALL GetNewUnit( UnIn )   
      CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2); if(Failed()) return


      CALL ReadCom( UnIn, FileName, 'MoorDyn water kinematics input file header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadCom( UnIn, FileName, 'MoorDyn water kinematics input file header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      ! ----- waves -----
      CALL ReadCom( UnIn, FileName,                               'waves header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadVar( UnIn, FileName, p%WaveKin  , 'WaveKinMod' ,  'WaveKinMod'   , ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadVar( UnIn, FileName, WaveKinFile, 'WaveKinFile',  'WaveKinFile'  , ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadVar( UnIn, FileName, p%dtWave   , 'dtWave', 'time step for waves', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadVar( UnIn, FileName, WaveDir    , 'WaveDir'    , 'wave direction', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      ! X grid points
      READ(UnIn,*, IOSTAT=ErrStat2)   coordtype         ! get the entry type
      READ(UnIn,'(A)', IOSTAT=ErrStat2)   entries2          ! get entries as string to be processed
      CALL gridAxisCoords(coordtype, entries2, p%pxWave, p%nxWave, ErrStat2, ErrMsg2)
      ! Y grid points
      READ(UnIn,*, IOSTAT=ErrStat2)   coordtype         ! get the entry type
      READ(UnIn,'(A)', IOSTAT=ErrStat2)   entries2          ! get entries as string to be processed
      CALL gridAxisCoords(coordtype, entries2, p%pyWave, p%nyWave, ErrStat2, ErrMsg2)
      ! Z grid points
      READ(UnIn,*, IOSTAT=ErrStat2)   coordtype         ! get the entry type
      READ(UnIn,'(A)', IOSTAT=ErrStat2)   entries2          ! get entries as string to be processed
      CALL gridAxisCoords(coordtype, entries2, p%pzWave, p%nzWave, ErrStat2, ErrMsg2)
      ! ----- current -----
      CALL ReadCom( UnIn, FileName,                        'current header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadVar( UnIn, FileName, p%Current,   'CurrentMod', 'CurrentMod', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadCom( UnIn, FileName,                'current profile header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      CALL ReadCom( UnIn, FileName,                'current profile header', ErrStat2, ErrMsg2, UnEcho); if(Failed()) return
      ! current profile table... (read through to end of file or ---)
      DO I=1,100
         READ(UnIn, *, IOSTAT=ErrStat2) pzCurrentTemp(i), uxCurrentTemp(i), uyCurrentTemp(i)     ! read into a line      
         if (ErrStat2 /= 0) then
            p%nzCurrent = i-1 ! save number of valid current depth points in profile
            EXIT      ! break out of the loop if it couldn't read the line (i.e. if at end of file)
         end if
         if (i == 100) then
            print*,"WARNING: MD can handle a maximum of 100 current profile points"
            exit
         end if
      END DO
      

      CLOSE(UnIn)

         
      ! ------------------- start with wave kinematics -----------------------

      ! WaveKin options: 0 - none or set externally during the sim (Waves object not needed unless there's current) [default]
      !                  1 - set externally for each node in each object (Waves object not needed unless there's current) (TBD)
      !                  2 - set from inputted wave elevation FFT, grid approach* (TBD)
      !                  3 - set from inputted wave elevation time series, grid approach* [supported]
      !                  4 - set from inputted wave elevation FFT, node approach (TBD)
      !                  5 - set from inputted wave elevation time series, node approach (TBD)
      !                  6 - set from inputted velocity, acceleration, and wave elevation grid data (TBD)**

      ! Current options: 0 - no currents or set externally (as part of WaveKin =0 or 1 approach) [default]
      !                  1 - read in steady current profile, grid approach (current_profile.txt)** [supported]
      !                  2 - read in dynamic current profile, grid approach (current_profile_dynamic.txt)** (TBD)
      !                  3 - read in steady current profile, node approach (current_profile.txt) (TBD)
      !                  4 - read in dynamic current profile, node approach (current_profile_dynamic.txt) (TBD)
       
      ! * the first call to any of these will attempt to load water_grid.txt to define the grid to put things on
      ! ** if a grid has already been set, these will interpolate onto it, otherwise they'll make a new grid based on their provided coordinates

      ! NOTE: lots of partial code is available from MD-C for supporting various wave kinematics input options
   
      ! WaveKin and Current compatibility check could go here in future
      
      
      ! --------------------- set from inputted wave elevation time series, grid approach -------------------
      if (p%WaveKin == 3) then

         print *, 'Setting up WaveKin 3 option: read wave elevation time series from file'

         IF ( LEN_TRIM( WaveKinFile ) == 0 )  THEN
            CALL SetErrStat( ErrID_Fatal,'WaveKinFile must not be an empty string.',ErrStat, ErrMsg, RoutineName); return
            RETURN
         END IF

         IF ( PathIsRelative( WaveKinFile ) ) THEN   ! properly handle relative path <<<
            !CALL GetPath( TRIM(InitInp%InputFile), TmpPath )
            WaveKinFile = TRIM(p%PriPath)//TRIM(WaveKinFile)
         END IF
         
         ! note: following is adapted from MoorDyn_Driver
         
         CALL GetNewUnit( UnElev ) 
      
         CALL OpenFInpFile ( UnElev, WaveKinFile, ErrStat2, ErrMsg2 ); if(Failed()) return
        
         print *, 'Reading wave elevation data from ', trim(WaveKinFile)
         
         ! Read through length of file to find its length
         i = 1  ! start counter
         DO
            READ(UnElev,'(A)',IOSTAT=ErrStat2) Line     !read into a line
            IF (ErrStat2 /= 0) EXIT      ! break out of the loop if it couldn't read the line (i.e. if at end of file)
            i = i+1
         END DO

         ! rewind to start of input file to re-read things now that we know how long it is
         REWIND(UnElev)      

         ntIn = i-3     ! save number of lines of file
         

         ! allocate space for input wave elevation array (including time column)
         CALL AllocAry(WaveTimeIn,  ntIn, 'WaveTimeIn', ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry(WaveElevIn,  ntIn, 'WaveElevIn', ErrStat2, ErrMsg2 ); if(Failed()) return

         ! read the data in from the file
         READ(UnElev,'(A)',IOSTAT=ErrStat2) Line     ! skip the first two lines as headers
         READ(UnElev,'(A)',IOSTAT=ErrStat2) Line     !
         
         DO i = 1, ntIn
            READ (UnElev, *, IOSTAT=ErrStat2) WaveTimeIn(i), WaveElevIn(i)
               
            IF ( ErrStat2 /= 0 ) THEN
              CALL SetErrStat( ErrID_Fatal,'Error reading WaveElev input file.',ErrStat, ErrMsg, RoutineName); return
            END IF 
         END DO  

         ! Close the inputs file 
         CLOSE ( UnElev ) 
         
         print *, "Read ", ntIn, " time steps from input file."

         ! if (WaveTimeIn(ntIn) < TMax) then <<<< need to handle if time series is too short?
           
         ! specify stepping details 
         p%ntWave = CEILING(Tmax/p%dtWave)          ! number of wave time steps

         
         ! allocate space for processed reference wave elevation time series
         ALLOCATE ( WaveElev0( 0:p%ntWave ), STAT=ErrStatTmp )  ! this has an extra entry of zero in case it needs to be padded to be even
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElev0.',ErrStat,ErrMsg,RoutineName)
         WaveElev0 = 0.0_SiKi
         
         ! go through and interpolate (should replace with standard function)
         DO i = 1, p%ntWave       
            t = p%dtWave*(i-1)
            
            ! interpolation routine 
            DO iIn = 1,ntIn-1      
               IF (WaveTimeIn(iIn+1) > t) THEN   ! find the right two points to interpolate between (remember that the first column of PtfmMotIn is time)
                  frac = (t - WaveTimeIn(iIn) )/( WaveTimeIn(iIn+1) - WaveTimeIn(iIn) )  ! interpolation fraction (0-1) between two interpolation points
                  WaveElev0(i-1) = WaveElevIn(iIn) + frac*(WaveElevIn(iIn+1) - WaveElevIn(iIn))  ! get interpolated wave elevation
                  EXIT   ! break out of the loop for this time step once we've done its interpolation
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
         ALLOCATE ( TmpFFTWaveElev( 0:NStepWave-1 ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpFFTWaveElev.',ErrStat,ErrMsg,RoutineName)

         ! Allocate frequency array for the wave elevation information in frequency space
         ALLOCATE ( WaveElevC0(2, 0:NStepWave2                ), STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElevC0.',ErrStat,ErrMsg,RoutineName)
         

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
         CALL InitFFT ( NStepWave, FFT_Data, .FALSE., ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName); if(Failed()) return

         ! Apply the forward FFT to get the real and imaginary parts of the frequency information.      
         CALL    ApplyFFT_f (  TmpFFTWaveElev(:), FFT_Data, ErrStatTmp )    ! Note that the TmpFFTWaveElev now contains the real and imaginary bits.
         CALL SetErrStat(ErrStatTmp,'Error occured while applying the forwards FFT to TmpFFTWaveElev array.',ErrStat,ErrMsg,RoutineName); if(Failed()) return

         ! Copy the resulting TmpFFTWaveElev(:) data over to the WaveElevC0 array
         DO I=1,NStepWave2-1
            WaveElevC0     (1,I) = TmpFFTWaveElev(2*I-1)
            WaveElevC0     (2,I) = TmpFFTWaveElev(2*I)
         ENDDO
         WaveElevC0(:,NStepWave2) = 0.0_SiKi

         CALL  ExitFFT(FFT_Data, ErrStatTmp)
         CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,RoutineName); if(Failed()) return


         IF (ALLOCATED( WaveElev0      )) DEALLOCATE( WaveElev0     , STAT=ErrStatTmp)
         IF (ALLOCATED( TmpFFTWaveElev )) DEALLOCATE( TmpFFTWaveElev, STAT=ErrStatTmp)


         
         ! note: following is a very streamlined adaptation from from Waves.v90 VariousWaves_Init
         
         ! allocate all the wave kinematics FFT arrays  
         ALLOCATE( WaveNmbr  (0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveNmbr.  ',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( tmpComplex(0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate tmpComplex.',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveElevC (0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveElevC .',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveDynPC (0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveDynPC .',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveVelCHx(0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveVelCHx.',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveVelCHy(0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveVelCHy.',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveVelCV (0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveVelCV .',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveAccCHx(0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveAccCHx.',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveAccCHy(0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveAccCHy.',ErrStat,ErrMsg,RoutineName)
         ALLOCATE( WaveAccCV (0:NStepWave2), STAT=ErrStatTmp); CALL SetErrStat(ErrStatTmp,'Cannot allocate WaveAccCV .',ErrStat,ErrMsg,RoutineName)
         
         ! allocate time series grid data arrays (now that we know the number of time steps coming from the IFFTs)
         CALL allocateKinematicsArrays() 
         
         
         ! Set the CosWaveDir and SinWaveDir values
         CosWaveDir=COS(D2R*WaveDir)
         SinWaveDir=SIN(D2R*WaveDir)
         
         ! get wave number array once
         DO I = 0, NStepWave2 
            WaveNmbr(i)   = WaveNumber ( REAL(I*WaveDOmega, R8Ki), p%g, p%WtrDpth )
            tmpComplex(I)    =  CMPLX(WaveElevC0(1,I), WaveElevC0(2,I))
         END DO    
         
         ! set up FFTer for doing IFFTs
         CALL InitFFT ( NStepWave, FFT_Data, .TRUE., ErrStatTmp )
         CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.', ErrStat, ErrMsg, routineName); if(Failed()) return

         ! Loop through all points where the incident wave kinematics will be computed      
         do ix = 1,p%nxWave 
            do iy = 1,p%nyWave
               do iz = 1,p%nzWave
                 
                  ! Compute the discrete Fourier transform of the incident wave kinematics
                  do i = 0, NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms

                     Omega = i*WaveDOmega
                     ImagOmega = ImagNmbr*Omega

                     WaveElevC (i) = tmpComplex(i) * EXP( -ImagNmbr*WaveNmbr(i)*( p%pxWave(ix)*CosWaveDir + p%pyWave(iy)*SinWaveDir ))                                                                 
                     WaveDynPC (i) = p%rhoW*p%g* WaveElevC(i) * COSHNumOvrCOSHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) )       
                     WaveVelCHx(i) =       Omega*WaveElevC(i) * COSHNumOvrSINHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) ) *CosWaveDir
                     WaveVelCHy(i) =       Omega*WaveElevC(i) * COSHNumOvrSINHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) ) *SinWaveDir             
                     WaveVelCV (i) =   ImagOmega*WaveElevC(i) * SINHNumOvrSINHDen( WaveNmbr(i), p%WtrDpth, REAL(p%pzWave(iz), R8Ki) )
                     WaveAccCHx(i) =   ImagOmega*WaveVelCHx(i)
                     WaveAccCHy(i) =   ImagOmega*WaveVelCHy(i)
                     WaveAccCV (i) =   ImagOmega*WaveVelCV (i)
                  end do  ! I, frequencies
                  
                  ! now IFFT all the wave kinematics except surface elevation and save it into the grid of data
                  CALL ApplyFFT_cx( p%PDyn  (:,iz,iy,ix), WaveDynPC , FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveDynP.', ErrStat,ErrMsg,RoutineName)
                  CALL ApplyFFT_cx( p%uxWave(:,iz,iy,ix), WaveVelCHx, FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveVelHx.',ErrStat,ErrMsg,RoutineName)
                  CALL ApplyFFT_cx( p%uyWave(:,iz,iy,ix), WaveVelCHy, FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveVelHy.',ErrStat,ErrMsg,RoutineName)
                  CALL ApplyFFT_cx( p%uzWave(:,iz,iy,ix), WaveVelCV , FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveVelV.', ErrStat,ErrMsg,RoutineName)
                  CALL ApplyFFT_cx( p%axWave(:,iz,iy,ix), WaveAccCHx, FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveAccHx.',ErrStat,ErrMsg,RoutineName)
                  CALL ApplyFFT_cx( p%ayWave(:,iz,iy,ix), WaveAccCHy, FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveAccHy.',ErrStat,ErrMsg,RoutineName)
                  CALL ApplyFFT_cx( p%azWave(:,iz,iy,ix), WaveAccCV , FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveAccV.', ErrStat,ErrMsg,RoutineName)

               end do ! iz
                 
               ! IFFT wave elevation here because it's only at the surface
               CALL ApplyFFT_cx( p%zeta(:,iy,ix) , WaveElevC , FFT_Data, ErrStatTmp ); CALL SetErrStat(ErrStatTmp,'Error IFFTing WaveElev.', ErrStat,ErrMsg,RoutineName)
            end do ! iy
         end do ! ix

         ! could also reproduce the wave elevation at 0,0,0 on a separate channel for verification...
         
         CALL  ExitFFT(FFT_Data, ErrStatTmp)
         CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the IFFTs.', ErrStat,ErrMsg,RoutineName); if(Failed()) return
      
      end if ! p%WaveKin == 3


      ! --------------------------------- now do currents --------------------------------
      if (p%Current == 1) then
      
         ! allocate current profile arrays to correct size
         CALL AllocAry( p%pzCurrent, p%nzCurrent, 'pzCurrent', ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry( p%uxCurrent, p%nzCurrent, 'uxCurrent', ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry( p%uyCurrent, p%nzCurrent, 'uyCurrent', ErrStat2, ErrMsg2 ); if(Failed()) return
         
         ! copy over data, flipping sign of depth values (to be positive-up) and reversing order
         do i = 1,p%nzCurrent
            p%pzCurrent(i) = -pzCurrentTemp(p%nzCurrent + 1 - i)  ! flip sign so depth is positive-up
            p%uxCurrent(i) =  uxCurrentTemp(p%nzCurrent + 1 - i) 
            p%uyCurrent(i) =  uyCurrentTemp(p%nzCurrent + 1 - i)
         end do
      
      end if  ! p%Current == 1


      ! ------------------------------ clean up and finished ---------------------------
      CALL cleanup()
      
      
   CONTAINS
   
   
      ! get grid axis coordinates, initialize/record in array, and return size
      SUBROUTINE gridAxisCoords(coordtype, entries, coordarray, n, ErrStat, ErrMsg)
         
         INTEGER(IntKi),          INTENT(IN   )  :: coordtype
         CHARACTER(*),            INTENT(INOUT)  :: entries
         REAL(SiKi), ALLOCATABLE,  INTENT(INOUT)  :: coordarray(:)
         INTEGER(IntKi),          INTENT(  OUT)  :: n
      
      
         INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat             ! Error status of the operation
         CHARACTER(*),            INTENT(  OUT)  :: ErrMsg              ! Error message if ErrStat /= ErrID_None
      
         REAL(ReKi) :: tempArray (100)
         REAL(ReKi) :: dx
         INTEGER(IntKi)                   :: nEntries, I
         
         ! get array of coordinate entries 
         CALL stringToArray(entries, nEntries, tempArray)
         
         ! set number of coordinates
         if (     coordtype==0) then   ! 0: not used - make one grid point at zero
            n = 1;
         else if (coordtype==1) then   ! 1: list values in ascending order
            n = nEntries
         else if (coordtype==2) then   ! 2: uniform specified by -xlim, xlim, num
            n = int(tempArray(3))
         else
            print *, "Error: invalid coordinate type specified to gridAxisCoords"
         end if
         
         ! allocate coordinate array
         CALL AllocAry(coordarray, n, 'x,y, or z grid points' , ErrStat, ErrMsg)
         !ALLOCATE ( coordarray(n), STAT=ErrStat) 
         
         ! fill in coordinates
         if (     coordtype==0) then
            coordarray(1) = 0.0_ReKi
         
         else if (coordtype==1) then
            coordarray(1:n) = tempArray(1:n)
         
         else if (coordtype==2) then  
            coordarray(1) = tempArray(1)
            coordarray(n) = tempArray(2)
            dx = (coordarray(n)-coordarray(0))/REAL(n-1)
            do i=2,n-1
               coordarray(i) = coordarray(1) + REAL(i)*dx
            end do
         
         else
            print *, "Error: invalid coordinate type specified to gridAxisCoords" 
         end if
         
         print *, "Set water grid coordinates to :"
         DO i=1,n
            print *, " ", coordarray(i)
         end do
         
      END SUBROUTINE gridAxisCoords
   
   
      ! Extract an array of numbers out of a string with comma-separated numbers (this could go in a more general location)
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
               READ(instring(pos1:), *) outarray(n)
               EXIT
            END IF
            n = n + 1
            if (n > 100) then
               print *, "ERROR - stringToArray cannot do more than 100 entries"
            end if            
            READ(instring(pos1:pos1+pos2-2), *) outarray(n)

            pos1 = pos2+pos1
         END DO
         
      END SUBROUTINE stringToArray
   
      
      ! allocate water kinematics arrays
      SUBROUTINE allocateKinematicsArrays()
        !  error check print *, "Error in Waves::makeGrid, a time or space array is size zero." << endl;

        ALLOCATE ( p%uxWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%uyWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%uzWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%axWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%ayWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%azWave( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%PDyn  ( p%ntWave,p%nzWave,p%nyWave,p%nxWave), STAT=ErrStatTmp)
        ALLOCATE ( p%zeta  ( p%ntWave,p%nyWave,p%nxWave), STAT = ErrStatTmp )    ! 2D grid over x and y only
        
      END SUBROUTINE allocateKinematicsArrays

      
      ! compact way to set the right error status and check if an abort is needed (and do cleanup if so)
      LOGICAL FUNCTION Failed()
           call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetupWaterKin') 
           Failed =  ErrStat >= AbortErrLev
           if (Failed) call CleanUp()
      END FUNCTION Failed
   

      SUBROUTINE CleanUp

         !IF (ALLOCATED( WaveElev   ))         DEALLOCATE( WaveElev,          STAT=ErrStatTmp)
         !IF (ALLOCATED( WaveTime   ))         DEALLOCATE( WaveTime,          STAT=ErrStatTmp)
         IF (ALLOCATED( TmpFFTWaveElev  ))    DEALLOCATE( TmpFFTWaveElev,    STAT=ErrStatTmp)
         IF (ALLOCATED( WaveElevC0      ))    DEALLOCATE( WaveElevC0,        STAT=ErrStatTmp)
         
         ! >>> missing some things <<<
         
         IF (ALLOCATED( WaveNmbr   )) DEALLOCATE( WaveNmbr   , STAT=ErrStatTmp)
         IF (ALLOCATED( tmpComplex )) DEALLOCATE( tmpComplex , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveElevC  )) DEALLOCATE( WaveElevC  , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveDynPC  )) DEALLOCATE( WaveDynPC  , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveVelCHx )) DEALLOCATE( WaveVelCHx , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveVelCHy )) DEALLOCATE( WaveVelCHy , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveVelCV  )) DEALLOCATE( WaveVelCV  , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveAccCHx )) DEALLOCATE( WaveAccCHx , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveAccCHy )) DEALLOCATE( WaveAccCHy , STAT=ErrStatTmp)
         IF (ALLOCATED( WaveAccCV  )) DEALLOCATE( WaveAccCV  , STAT=ErrStatTmp)

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
   
   



END MODULE MoorDyn_Misc
