   SUBROUTINE ComputeIniNodalCrv(EndP1,EndP2,MidP,phi,eta,cc)
   !-----------------------------------------------------
   ! This subroutine computes initial CRV parameters
   ! given geometry information
   !-----------------------------------------------------

   REAL(ReKi),INTENT(IN   ):: EndP1(:)    ! Position vector of End Point 1
   REAL(ReKi),INTENT(IN   ):: EndP2(:)    ! Position vector of End Point 2
   REAL(ReKi),INTENT(IN   ):: MidP(:)     ! Position vector of Mid Point
   REAL(ReKi),INTENT(IN   ):: eta         ! Nondimensional coordinate, [-1,1]
   REAL(ReKi),INTENT(IN   ):: phi         ! Initial twist angle
   REAL(ReKi),INTENT(  OUT):: cc(:)       ! Initial Crv Parameter

   REAL(ReKi):: hpx(3)                    ! Derivatives of second-order shape function 
   REAL(ReKi):: e1(3)                     ! Unit tangent vector
   REAL(ReKi):: e2(3)                     ! Unit normal vector
   REAL(ReKi):: e3(3)                     ! Unit e3 = e1 * e2, cross-product
   REAL(ReKi):: Rr(3,3)                   ! Initial rotation matrix
   REAL(ReKi):: temp
   REAL(ReKi):: temp2
   REAL(ReKi):: Delta
   INTEGER(IntKi):: i

   hpx = 0.0D0
   hpx(1) = 0.5D0*(2.0D0*eta - 1.0D0)
   hpx(2) = -2.0D0*eta
   hpx(3) = 0.5D0*(2.0D0*eta + 1.0D0)

   Rr = 0.0D0

   e1 = 0.0D0
   DO i=1,3
       e1(i) = hpx(1)*EndP1(i) + hpx(2)* MidP(i) + hpx(3)*EndP2(i) 
   ENDDO
   e1 = e1/Norm(e1)
   DO i=1,3
       Rr(i,1) = e1(i)
   ENDDO

   e2 = 0.0D0
   temp = phi*ACOS(-1.0D0)/180.0D0
   temp2 = ((e1(2)*COS(temp) + e1(3)*SIN(temp))/e1(1))
   Delta = SQRT(1.0D0 + temp2*temp2)
   e2(1) = -(e1(2)*COS(temp)+e1(3)*SIN(temp))/e1(1)
   e2(2) = COS(temp)
   e2(3) = SIN(temp)
   e2 = e2/Delta
   DO i=1,3
       Rr(i,2) = e2(i)
   ENDDO

   e3 = 0.0D0
   e3 = CrossProduct(e1,e2)
   DO i=1,3
       Rr(i,3) = e3(i)
   ENDDO

   CALL CrvExtractCrv(Rr,cc)

   END SUBROUTINE ComputeIniNodalCrv
