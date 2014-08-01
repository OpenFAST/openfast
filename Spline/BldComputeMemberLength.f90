   SUBROUTINE BldComputeMemberLength(member_total,Coef,member_length,total_length)

   INTEGER(IntKi),INTENT(IN   ):: member_total
   REAL(ReKi),    INTENT(IN   ):: Coef(:,:,:)
   REAL(ReKi),    INTENT(  OUT):: member_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: total_length

   REAL(ReKi)                  :: sample_step
   REAL(ReKi)                  :: eta0
   REAL(ReKi)                  :: eta1
   REAL(ReKi)                  :: temp_pos0(3)
   REAL(ReKi)                  :: temp_pos1(3)
   REAL(ReKi)                  :: sample_step
   INTEGER(IntKi)              :: sample_total
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: k

   sample_total = 1001
   sample_step = 1.0D0/(sample_total-1)

   DO i=1,member_total
       DO j=1,sample_total-1
           eta0 = (j-1)*sample_step
           eta1 = j*sample_step
           DO k=1,3
               temp_pos0(k) = Coef(i,1,k) + Coef(i,2,k)*eta0 + Coef(i,3,k)*eta0*eta0 + Coef(i,4,k)*eta0*eta0*eta0
               temp_pos1(k) = Coef(i,1,k) + Coef(i,2,k)*eta1 + Coef(i,3,k)*eta1*eta1 + Coef(i,4,k)*eta1*eta1*eta1
           ENDDO
           temp_pos1(:) = temp_pos1(:) - temp_pos0(:)
           member_length(i,1) = member_length(i,1) + SQRT(DOT_PRODUCT(temp_pos1,temp_pos1))
       ENDDO
       total_length = total_length + member_length(i,1)
   ENDDO


   DO i=1,member_total
       member_length(i,2) = member_length(i,1)/total_length
   ENDDO


   END SUBROUTINE BldComputeMemberLength
