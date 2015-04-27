   SUBROUTINE BD_MotionTensor(RotTen,Pos,MotTen,flag)

   REAL(ReKi),     INTENT(IN   ):: RotTen(:,:) 
   REAL(ReKi),     INTENT(IN   ):: Pos(:) 
   REAL(ReKi),     INTENT(  OUT):: MotTen(:,:) 
   INTEGER(IntKi), INTENT(IN   ):: flag            ! 0: Motion Tensor; 
                                                   ! 1: Inverse of Motion Tensor
   
   MotTen(:,:) = 0.0D0
   IF (flag .EQ. 0) THEN
       MotTen(1:3,1:3) = RotTen(1:3,1:3)
       MotTen(4:6,4:6) = RotTen(1:3,1:3)
       MotTen(1:3,4:6) = MATMUL(Tilde(Pos),RotTen)
   ELSEIF(flag .EQ. 1) THEN
       MotTen(1:3,1:3) = TRANSPOSE(RotTen(1:3,1:3))
       MotTen(4:6,4:6) = TRANSPOSE(RotTen(1:3,1:3))
       MotTen(1:3,4:6) = TRANSPOSE(MATMUL(Tilde(Pos),RotTen))
   ENDIF

   END SUBROUTINE BD_MotionTensor 
