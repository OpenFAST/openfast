   SUBROUTINE ComputeIniNodalPosition(EndP1,EndP2,MidP,eta,PosiVec,e1)

   REAL(ReKi),INTENT(IN   ):: EndP1(:)
   REAL(ReKi),INTENT(IN   ):: EndP2(:)
   REAL(ReKi),INTENT(IN   ):: MidP(:)
   REAL(ReKi),INTENT(IN   ):: eta
   REAL(ReKi),INTENT(  OUT):: PosiVec(:)
   REAL(ReKi),INTENT(  OUT):: e1(:)

   REAL(ReKi)              :: temp_a
   REAL(ReKi)              :: temp_b
   REAL(ReKi)              :: temp_c
   REAL(ReKi)              :: temp1
   REAL(ReKi)              :: temp2
   REAL(ReKi)              :: temp_area
   REAL(ReKi)              :: temp_radius
   REAL(ReKi)              :: temp_angle
   REAL(ReKi)              :: temp_v1(3)
   REAL(ReKi)              :: temp_v2(3)
   REAL(ReKi)              :: temp_v3(3)
   REAL(ReKi)              :: temp_v4(3)
   REAL(ReKi)              :: temp_center(3)
   REAL(ReKi),PARAMETER    :: eps = 1.0D-10
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
   temp_v1(:) = MidP(:) - EndP1(:)
   temp_v2(:) = EndP2(:) - EndP1(:)
   temp_v3(:) = CrossProduct(temp_v1,temp_v2)
   temp_area = 0.5D0*Norm(temp_v3)
   PosiVec(:) = 0.0D0 
   e1(:) = 0.0D0
   IF(temp_area <= eps) THEN
       temp_a = (eta + 1.0D0)/2.0D0
       PosiVec(:) = (1.0D0 - temp_a)*EndP1(:) + temp_a*EndP2(:)
       e1(:) = EndP2(:) - EndP1(:)
       e1(:) = e1(:)/Norm(e1)
   ELSE
       temp_a = Norm(EndP1-MidP)
       temp_b = Norm(MidP-EndP2)
       temp_c = Norm(EndP2-EndP1)
       temp_a = temp_a*temp_b*temp_c
       temp_v1(:) = EndP1(:) - MidP(:)
       temp_v2(:) = MidP(:) - EndP2(:)
       temp_v3(:) = CrossProduct(temp_v1,temp_v2)
       temp_b = 2.0D0*Norm(temp_v3)
       temp_radius = temp_a/temp_b

       temp_v1(:) = EndP1(:) - MidP(:)
       temp_v2(:) = MidP(:) - EndP2(:)
       temp_v3(:) = EndP1(:) - EndP2(:)
       temp_v4(:) = CrossProduct(temp_v1,temp_v2)
       temp1 = 2.0D0*Norm(temp_v4)*Norm(temp_v4)
       temp2 = Norm(temp_v2)*Norm(temp_v2)
       temp_a = temp2*DOT_PRODUCT(temp_v1,temp_v3)/temp1
       temp2 = Norm(temp_v3)*Norm(temp_v3)
       temp_b = temp2*DOT_PRODUCT(-temp_v1,temp_v2)/temp1
       temp2 = Norm(temp_v1)*Norm(temp_v1)
       temp_c = temp2*DOT_PRODUCT(-temp_v3,-temp_v2)/temp1
       temp_center(:) = temp_a*EndP1(:) + temp_b*MidP(:) + temp_c*EndP2(:)
   
       temp_v1(:) = EndP1(:) - temp_center(:)
       temp_v1(:) = temp_v1(:)/Norm(temp_v1)
       temp_v2(:) = EndP2(:) - temp_center(:)
       temp_v2(:) = temp_v2(:)/Norm(temp_v2)
       temp_angle = ACOS(DOT_PRODUCT(temp_v1,temp_v2))
       temp_a = (eta + 1.0D0)*temp_angle/2.0D0
       temp_v3(:) = CrossProduct(temp_v1,temp_v2)
       temp_v3(:) = temp_v3(:)/Norm(temp_v3)
       PosiVec(:) = temp_center(:) + temp_radius*COS(temp_a)*temp_v1(:)&
                   & + temp_radius*SIN(temp_a)*CrossProduct(temp_v3,temp_v1)
       e1(:) = -temp_radius*SIN(temp_a)*temp_v1(:) + temp_radius*COS(temp_a)*CrossProduct(temp_v3,temp_v1)
       e1(:) = e1(:)/Norm(e1)
       WRITE(*,*) "e1: ",e1(:)
   ENDIF


  
   END SUBROUTINE ComputeIniNodalPosition
