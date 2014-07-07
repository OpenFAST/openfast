   SUBROUTINE ComputeIniNodalPosition(EndP1,EndP2,MidP,eta,PosiVec)

   REAL(ReKi),INTENT(IN):: EndP1(:)
   REAL(ReKi),INTENT(IN):: EndP2(:)
   REAL(ReKi),INTENT(IN):: MidP(:)
   REAL(ReKi),INTENT(IN):: eta
   REAL(ReKi),INTENT(OUT):: PosiVec(:)

   REAL(ReKi):: temp_a
   REAL(ReKi):: temp_b
   REAL(ReKi):: temp_c
!   REAL(ReKi):: hhx(3)
!   INTEGER(IntKi):: i

!   hhx = 0.0D0
!   hhx(1) = 0.5D0*(eta*eta - eta)
!   hhx(2) = 1.0D0 - eta*eta
!   hhx(3) = 0.5D0*(eta*eta + eta)
  
!   PosiVec = 0.0D0
!   DO i=1,3
!       PosiVec(i) = hhx(1)*EndP1(i) + hhx(2)*MidP(i) + hhx(3)*EndP2(i) 
!   ENDDO
   temp_v1(1:3) = MidP(1:3) - EndP1(1:3)
   temp_v2(1:3) = EndP2(1:3) - EndP1(1:3)
   temp_v3(1:3) = CrossProduct(temp_v1,temp_v2)
   CALL Norm(3,temp_v3,temp_area)
   temp_area = 0.5D0*temp_area
   PosiVec(:) = 0.0D0 
   IF(temp_area <= eps) THEN
       temp_a = (eta + 1.0D0)/2.0D0
       PosiVec(1:3) = (1.0D0 - temp_a)*EndP1(1:3) + temp_a*EndP2(1:3)
   ELSE
       CALL Norm(3,EndP1-MidP,temp_a)
       CALL Norm(3,MidP-EndP2,temp_b)
       CALL Norm(3,EndP2-EndP1,temp_c)
       temp_a = temp_a*temp_b*temp_c
       temp_v1(1:3) = EndP1(1:3) - MidP(1:3)
       temp_v2(1:3) = MidP(1:3) - EndP2(1:3)
       temp_v3(1:3) = CrossProduct(temp_v1,temp_v2)
       CALL Norm(3,temp_v3,temp_b)
       temp_radius = 0.5D0*temp_a/temp_b
       temp_v1(1:3) = EndP1(1:3) - MidP(1:3)
       temp_v2(1:3) = MidP(1:3) - EndP2(1:3)
       temp_v3(1:3) = EndP1(1:3) - EndP2(1:3)
       temp_v4 = CrossProduct(temp_v1,temp_v2)
       CALL Norm(3,temp_v4,temp1)
       temp1 = 2.0D0*temp1*temp1
       CALL Norm(3,temp_v2,temp2)
       temp2 = temp2*temp2
       temp_a = temp2*DOT_PRODUCT(temp_v1,temp_v3)/temp1
       CALL Norm(3,temp_v3,temp2)
       temp2 = temp2*temp2
       temp_b = temp2*DOT_PRODUCT(-temp_v1,temp_v2)/temp1
       CALL Norm(3,temp_v1,temp2)
       temp2 = temp2*temp2
       temp_c = temp2*DOT_PRODUCT(-temp_v3,-temp_v2)/temp1
       temp_center(1:3) = temp_a*EndP1(1:3) + temp_b*MidP(1:3) + temp_c*EndP2(1:3)
   
       temp_v1(1:3) = EndP1(1:3) - temp_center(1:3)
       CALL Norm(3,temp_v1,temp1)
       temp_v1(:) = temp_v1(:)/temp1
       temp_v2(1:3) = EndP2(1:3) - temp_center(1:3)
       CALL Norm(3,temp_v2,temp1)
       temp_v2(:) = temp_v2(:)/temp1
       temp_angle = ACOS(DOT_PRODUCT(temp_v1,temp_v2))
       temp_a = (eta + 1.0D0)*temp_angle/2.0D0
       temp_v3(:) = CrossProduct(temp_v1,temp_v2)
       CALL Norm(3,temp_v3,temp1)
       temp_v3(:) = temp_v3(:)/temp1
       PosiVec(:) = temp_center(:) + temp_radius*COS(temp_a)*temp_v1(:)&
                   & + temp_radius*SIN(temp_a)*CrossProduct(temp_v3,temp_v1)
   ENDIF


  
   END SUBROUTINE ComputeIniNodalPosition
