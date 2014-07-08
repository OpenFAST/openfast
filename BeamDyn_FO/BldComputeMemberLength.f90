   SUBROUTINE BldComputeMemberLength(member_total,kp_coordinate,member_length,total_length)

   INTEGER(IntKi),INTENT(IN   ):: member_total
   REAL(ReKi),    INTENT(IN   ):: kp_coordinate(:,:)

   REAL(ReKi),    INTENT(  OUT):: member_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: total_length

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: temp_id
   REAL(ReKi)                  :: temp_area
   REAL(ReKi)                  :: x1
   REAL(ReKi)                  :: x2
   REAL(ReKi)                  :: x3
   REAL(ReKi)                  :: y1
   REAL(ReKi)                  :: y2
   REAL(ReKi)                  :: y3
   REAL(ReKi)                  :: z1
   REAL(ReKi)                  :: z2
   REAL(ReKi)                  :: z3
   REAL(ReKi),        PARAMETER:: eps = 1.0D-10

   x1 = 0.0D0
   x2 = 0.0D0
   x3 = 0.0D0
   y1 = 0.0D0
   y2 = 0.0D0
   y3 = 0.0D0
   z1 = 0.0D0
   z2 = 0.0D0
   z3 = 0.0D0
   temp_area = 0.0D0


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
           member_length(i,1) = MemberArcLength(x1,x2,x3,y1,y2,y3,z1,z2,z3)
       ENDIF
       total_length = total_length + member_length(i,1)
   ENDDO


   DO i=1,member_total
       member_length(i,2) = member_length(i,1)/total_length
   ENDDO


   END SUBROUTINE BldComputeMemberLength
