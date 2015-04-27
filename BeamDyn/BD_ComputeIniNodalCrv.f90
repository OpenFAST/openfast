   SUBROUTINE BD_ComputeIniNodalCrv(e1,phi,cc)
   !-----------------------------------------------------
   ! This subroutine computes initial CRV parameters
   ! given geometry information
   !-----------------------------------------------------

   REAL(ReKi),INTENT(IN   ):: e1(:)       ! Unit tangent vector
   REAL(ReKi),INTENT(IN   ):: phi         ! Initial twist angle
   REAL(ReKi),INTENT(  OUT):: cc(:)       ! Initial Crv Parameter

   REAL(ReKi):: e2(3)                     ! Unit normal vector
   REAL(ReKi):: e3(3)                     ! Unit e3 = e1 * e2, cross-product
   REAL(ReKi):: Rr(3,3)                   ! Initial rotation matrix
   REAL(ReKi):: temp
   REAL(ReKi):: temp2
   REAL(ReKi):: Delta
   INTEGER(IntKi):: i


   Rr = 0.0D0
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

   CALL BD_CrvExtractCrv(Rr,cc)

   END SUBROUTINE BD_ComputeIniNodalCrv
