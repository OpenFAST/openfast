   FUNCTION MemberArcLength(EndP1,MidP,EndP2)

   REAL(ReKi),    INTENT(IN):: EndP1(:)
   REAL(ReKi),    INTENT(IN):: MidP(:)
   REAL(ReKi),    INTENT(IN):: EndP2(:)
   REAL(ReKi)               :: MemberArcLength

   REAL(ReKi)               :: temp_a
   REAL(ReKi)               :: temp_b
   REAL(ReKi)               :: temp_c
   REAL(ReKi)               :: temp_radius
   REAL(ReKi)               :: temp_center
   REAL(ReKi)               :: temp_angle
   REAL(ReKi)               :: temp_v1(3)
   REAL(ReKi)               :: temp_v2(3)
   REAL(ReKi)               :: temp_v3(3)
   REAL(ReKi)               :: temp_v4(3)


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

   MemberArcLength = temp_angle*temp_radius

   END FUNCTION MemberArcLength
