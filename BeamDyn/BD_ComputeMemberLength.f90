   SUBROUTINE BD_ComputeMemberLength(member_total,kp_member,Coef,seg_length,member_length,total_length)

   INTEGER(IntKi),INTENT(IN   ):: member_total
   INTEGER(IntKi),INTENT(IN   ):: kp_member(:)
   REAL(ReKi),    INTENT(IN   ):: Coef(:,:,:)
   REAL(ReKi),    INTENT(  OUT):: seg_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: member_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: total_length

   REAL(ReKi)                  :: eta0
   REAL(ReKi)                  :: eta1
   REAL(ReKi)                  :: temp_pos0(3)
   REAL(ReKi)                  :: temp_pos1(3)
   REAL(ReKi)                  :: sample_step
   REAL(ReKi)                  :: temp
   INTEGER(IntKi)              :: sample_total
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: k
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id

   sample_total = 101
   sample_step = 1.0D0/(sample_total-1)

   temp_id = 0
   DO i=1,member_total
       DO m=1,kp_member(i)-1
           temp_id = temp_id + 1
           DO j=1,sample_total-1
               eta0 = (j-1)*sample_step
               eta1 = j*sample_step
               DO k=1,3
                   temp_pos0(k) = Coef(temp_id,1,k) + Coef(temp_id,2,k)*eta0 + Coef(temp_id,3,k)*eta0*eta0 + Coef(temp_id,4,k)*eta0*eta0*eta0
                   temp_pos1(k) = Coef(temp_id,1,k) + Coef(temp_id,2,k)*eta1 + Coef(temp_id,3,k)*eta1*eta1 + Coef(temp_id,4,k)*eta1*eta1*eta1
               ENDDO
               temp_pos1(:) = temp_pos1(:) - temp_pos0(:)
               seg_length(temp_id,1) = seg_length(temp_id,1) + SQRT(DOT_PRODUCT(temp_pos1,temp_pos1))
               member_length(i,1) = member_length(i,1) + SQRT(DOT_PRODUCT(temp_pos1,temp_pos1))
           ENDDO
       ENDDO
       total_length = total_length + member_length(i,1)
   ENDDO

   temp_id = 0
   DO i=1,member_total
       temp = 0.0D0
       DO j=1,kp_member(i)-1
           temp_id = temp_id + 1
           IF(j == 1) THEN 
               seg_length(temp_id,2) = 0.0D0
           ELSE
               seg_length(temp_id,2) = seg_length(temp_id-1,3)
           ENDIF
           temp = temp + seg_length(temp_id,1)
           seg_length(temp_id,3) = temp/member_length(i,1)     
       ENDDO    
   ENDDO

   DO i=1,member_total
       member_length(i,2) = member_length(i,1)/total_length
   ENDDO


   END SUBROUTINE BD_ComputeMemberLength
