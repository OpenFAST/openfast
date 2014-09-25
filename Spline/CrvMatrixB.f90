   SUBROUTINE CrvMatrixB(cc,aa,Bb)

   REAL(ReKi),INTENT(IN   ):: cc(:)
   REAL(ReKi),INTENT(IN   ):: aa(:)
   REAL(ReKi),INTENT(  OUT):: Bb(:,:)

   REAL(ReKi)              :: c0
   REAL(ReKi)              :: temp0
   REAL(ReKi)              :: temp1
   REAL(ReKi)              :: temp2
   REAL(ReKi)              :: temp11
   REAL(ReKi)              :: temp22
   REAL(ReKi)              :: temp33(3,3)
  
   c0 = 2.0D0 - 0.125*DOT_PRODUCT(cc,cc)

   temp0 = 1.0D0/(2.0D0*(4.0D0-c0)*(4.0D0-c0))
   temp1 = 1.0D0/((4.0D0-c0)*(4.0D0-c0)*(4.0D0-c0))
   temp2 = 0.25D0*1.0D0/((4.0D0-c0)*(4.0D0-c0)*(4.0D0-c0))
   temp11= 2.0D0/((4.0D0-c0)*(4.0D0-c0))
   temp22= temp0

   Bb(:,:) = 0.0D0
   temp33(:,:) = OuterProduct(aa,cc)!MATMUL(aa,TRANSPOSE(cc))

   Bb(:,:) = temp0 * temp33(:,:) + &
            &temp1 * MATMUL(Tilde(cc),temp33) + &
            &temp2 * MATMUL(Tilde(cc),MATMUL(Tilde(cc),temp33)) - &
            &temp11 * Tilde(aa) - &
            &temp22 * 2.0D0 * MATMUL(Tilde(cc),Tilde(aa)) + &
            &temp22 * MATMUL(Tilde(aa),Tilde(cc))

   END SUBROUTINE CrvMatrixB
