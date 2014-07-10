   SUBROUTINE BldComputeMemberLength(member_total,kp_coordinate,member_length,total_length)

   INTEGER(IntKi),INTENT(IN   ):: member_total
   REAL(ReKi),    INTENT(IN   ):: kp_coordinate(:,:)

   REAL(ReKi),    INTENT(  OUT):: member_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: total_length

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: temp_id
   REAL(ReKi)                  :: temp_area
   REAL(ReKi)                  :: EndP1(3)
   REAL(ReKi)                  :: MidP(3)
   REAL(ReKi)                  :: EndP2(3)
   REAL(ReKi)                  :: temp_v1(3)
   REAL(ReKi)                  :: temp_v2(3)
   REAL(ReKi)                  :: temp_v3(3)
   REAL(ReKi),        PARAMETER:: eps = 1.0D-10

   temp_area  = 0.0D0
   EndP1(:)   = 0.0D0
   MidP(:)    = 0.0D0
   EndP2(:)   = 0.0D0
   temp_v1(:) = 0.0D0
   temp_v2(:) = 0.0D0
   temp_v3(:) = 0.0D0


   DO i=1,member_total
       temp_id = (i - 1)*2
       EndP1(1:3) = kp_coordinate(temp_id+1,1:3)
       MidP(1:3)  = kp_coordinate(temp_id+2,1:3)
       EndP2(1:3) = kp_coordinate(temp_id+3,1:3)
       temp_v1(:) = MidP(:) - EndP1(:)
       temp_v2(:) = EndP2(:) - EndP1(:)
       temp_v3(:) = CrossProduct(temp_v1,temp_v2)
       temp_area = 0.5D0*Norm(temp_v3)
       IF(temp_area <= eps) THEN
           member_length(i,1) = SQRT((EndP1(1)-EndP2(1))**2+(EndP1(2)-EndP2(2))**2+(EndP1(3)-EndP2(3))**2)
       ELSE
           member_length(i,1) = MemberArcLength(EndP1,MidP,EndP2)
       ENDIF
       total_length = total_length + member_length(i,1)
   ENDDO


   DO i=1,member_total
       member_length(i,2) = member_length(i,1)/total_length
   ENDDO


   END SUBROUTINE BldComputeMemberLength
