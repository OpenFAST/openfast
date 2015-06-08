   SUBROUTINE AM2LinearizationMatrix(uuu,vvv,uuu0,vvv0,m00,mEta,rho,&
                                    &A2,A3,A4,A5,A6,A7)

   REAL(ReKi),INTENT(IN   ):: uuu(:)
   REAL(ReKi),INTENT(IN   ):: vvv(:)
   REAL(ReKi),INTENT(IN   ):: uuu0(:)
   REAL(ReKi),INTENT(IN   ):: vvv0(:)
   REAL(ReKi),INTENT(IN   ):: m00
   REAL(ReKi),INTENT(IN   ):: mEta(:)
   REAL(ReKi),INTENT(IN   ):: rho(:,:)
!   REAL(ReKi),INTENT(  OUT):: A1(:,:)
   REAL(ReKi),INTENT(  OUT):: A2(:,:)
   REAL(ReKi),INTENT(  OUT):: A3(:,:)
   REAL(ReKi),INTENT(  OUT):: A4(:,:)
   REAL(ReKi),INTENT(  OUT):: A5(:,:)
   REAL(ReKi),INTENT(  OUT):: A6(:,:)
   REAL(ReKi),INTENT(  OUT):: A7(:,:)

   REAL(ReKi)              :: temp_B(3,3)
   REAL(ReKi)              :: temp_B0(3,3)
   REAL(ReKi)              :: temp_H0(3,3)
   REAL(ReKi)              :: temp_H(3,3)
   REAL(ReKi)              :: vel(3)
   REAL(ReKi)              :: omega(3)
   REAL(ReKi)              :: vel0(3)
   REAL(ReKi)              :: omega0(3)
   INTEGER(IntKi)          :: i

   A6(:,:) = 0.0D0
   DO i=1,3
       A6(i,i) = 2.0D0
   ENDDO
   CALL CrvMatrixB(uuu(4:6),uuu(4:6),temp_B)
   CALL CrvMatrixB(uuu(4:6),uuu0(4:6),temp_B0)
   CALL CrvMatrixH(uuu0(4:6),temp_H0)
   CALL CrvMatrixH(uuu(4:6),temp_H)
   temp_B(:,:) = temp_B(:,:) + temp_H(:,:) - temp_B0(:,:) + temp_H0(:,:)
!   temp_B(:,:) = temp_B(:,:) - temp_B0(:,:) + temp_H0(:,:)
   A6(4:6,4:6) = temp_B(1:3,1:3)

   A7(:,:) = 0.0D0
   DO i=1,6
       A7(i,i) = 1.0D0
   ENDDO

!   A1(:,:) = 0.0D0
!   CALL MassMatrix(m00,mEta,rho,A1)

   vel(1:3) = vvv(1:3)
   omega(1:3) = vvv(4:6)
   vel0(1:3) = vvv0(1:3)
   omega0(1:3) = vvv0(4:6)

   A2(:,:) = 0.0D0
   A2(1:3,4:6) = MATMUL(Tilde(omega),TRANSPOSE(Tilde(mEta)))
   A2(4:6,4:6) = MATMUL(Tilde(vel),Tilde(mEta)) + &
                &MATMUL(rho,Tilde(omega)) - &
                &Tilde(MATMUL(rho,omega))

   A3(:,:) = 0.0D0
   A3(1:3,4:6) = MATMUL(Tilde(omega0),TRANSPOSE(Tilde(mEta)))
   A3(4:6,4:6) = MATMUL(Tilde(vel0),Tilde(mEta)) + &
                &MATMUL(rho,Tilde(omega0)) - &
                &Tilde(MATMUL(rho,omega0))

   A4(:,:) = 0.0D0
   A4(1:3,4:6) = TRANSPOSE(Tilde(MATMUL(Tilde(omega),mEta))) + &
                &MATMUL(Tilde(omega),TRANSPOSE(Tilde(mEta)))
   A4(4:6,4:6) = TRANSPOSE(Tilde(MATMUL(rho,omega))) + &
                &MATMUL(Tilde(omega),rho)
   
   A5(:,:) = 0.0D0
   A5(1:3,4:6) = MATMUL(Tilde(omega),MATMUL(Tilde(omega),TRANSPOSE(Tilde(mEta))))
   A5(4:6,4:6) = MATMUL(Tilde(omega),MATMUL(rho,Tilde(omega))) - &
                &MATMUL(Tilde(omega),Tilde(MATMUL(rho,omega)))


   END SUBROUTINE AM2LinearizationMatrix
