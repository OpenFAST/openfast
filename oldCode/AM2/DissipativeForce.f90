   SUBROUTINE DissipativeForce(beta,Stiff,vvv,vvp,E1,Fc,Fd,&
                              &Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)

   REAL(ReKi),INTENT(IN   ):: beta(:)
   REAL(ReKi),INTENT(IN   ):: Stiff(:,:)
   REAL(ReKi),INTENT(IN   ):: vvv(:)
   REAL(ReKi),INTENT(IN   ):: vvp(:)
   REAL(ReKi),INTENT(IN   ):: E1(:)
   REAL(ReKi),INTENT(INOUT):: Fc(:)
   REAL(ReKi),INTENT(INOUT):: Fd(:)
   REAL(ReKi),INTENT(  OUT):: Sd(:,:)
   REAL(ReKi),INTENT(  OUT):: Od(:,:)
   REAL(ReKi),INTENT(  OUT):: Pd(:,:)
   REAL(ReKi),INTENT(  OUT):: Qd(:,:)
   REAL(ReKi),INTENT(  OUT):: betaC(:,:)
   REAL(ReKi),INTENT(  OUT):: Gd(:,:)
   REAL(ReKi),INTENT(  OUT):: Xd(:,:)
   REAL(ReKi),INTENT(  OUT):: Yd(:,:)

   REAL(ReKi)              :: ome(3)
   REAL(ReKi)              :: eed(6)
   REAL(ReKi)              :: ffd(6)
   REAL(ReKi)              :: D11(3,3)
   REAL(ReKi)              :: D12(3,3)
   REAL(ReKi)              :: D21(3,3)
   REAL(ReKi)              :: D22(3,3)
   REAL(ReKi)              :: b11(3,3)
   REAL(ReKi)              :: b12(3,3)
   REAL(ReKi)              :: alpha(3,3)
   REAL(ReKi)              :: temp_b(6,6)
   INTEGER(IntKi)          :: i

   ome(1:3) = vvv(4:6)
   !---------------------------
   ! Compute strain rates
   !---------------------------
   eed(1:6) = vvp(1:6)
   eed(1:3) = eed(1:3) + MATMUL(Tilde(E1),ome)

   !---------------------------
   ! Compute damping matrix
   !---------------------------
   temp_b(:,:) = 0.0D0
   DO i=1,6
       temp_b(i,i) = beta(i)
   ENDDO
   betaC(:,:) = 0.0D0
   betaC(1:6,1:6) = MATMUL(temp_b,Stiff(:,:))
   D11(1:3,1:3) = betaC(1:3,1:3)
   D12(1:3,1:3) = betaC(1:3,4:6)
   D21(1:3,1:3) = betaC(4:6,1:3)
   D22(1:3,1:3) = betaC(4:6,4:6)

   !---------------------------
   ! Compute dissipative force
   !---------------------------
   ffd(1:6) = MATMUL(betaC,eed)
   Fc(1:6) = Fc(1:6) + ffd(1:6)
   Fd(4:6) = Fd(4:6) + MATMUL(Tilde(ffd(1:3)),E1)

   !----------------------------
   ! Compute stiffness matrix Sd
   !----------------------------
   Sd(:,:) = 0.0D0
   Sd(1:3,1:3) = MATMUL(D11,TRANSPOSE(Tilde(ome)))
   Sd(1:3,4:6) = MATMUL(D12,TRANSPOSE(Tilde(ome)))
   Sd(4:6,1:3) = MATMUL(D21,TRANSPOSE(Tilde(ome)))
   Sd(4:6,4:6) = MATMUL(D22,TRANSPOSE(Tilde(ome)))

   !----------------------------
   ! Compute stiffness matrix Pd
   !----------------------------
   Pd(:,:) = 0.0D0
   b11(1:3,1:3) = MATMUL(TRANSPOSE(Tilde(E1)),D11)
   b12(1:3,1:3) = MATMUL(TRANSPOSE(Tilde(E1)),D12)
   Pd(4:6,1:3) = Tilde(ffd(1:3)) + MATMUL(b11,TRANSPOSE(Tilde(ome)))
   Pd(4:6,1:3) = MATMUL(b12,TRANSPOSE(Tilde(ome)))

   !----------------------------
   ! Compute stiffness matrix Od
   !----------------------------
   Od(:,:) = 0.0D0
   alpha(1:3,1:3) = Tilde(vvp(1:3)) - MATMUL(Tilde(ome),Tilde(E1))
   Od(1:3,4:6) = MATMUL(D11,alpha) - Tilde(ffd(1:3))
   Od(4:6,4:6) = MATMUL(D21,alpha) - Tilde(ffd(4:6))

   !----------------------------
   ! Compute stiffness matrix Qd
   !----------------------------
   Qd(:,:) = 0.0D0
   Qd(4:6,4:6) = MATMUL(TRANSPOSE(Tilde(E1)),Od(1:3,4:6))

   !-----------------------------
   ! Compute gyroscopic matrix Gd
   !-----------------------------
   Gd(:,:) = 0.0D0
   Gd(1:3,4:6) = TRANSPOSE(b11)
   Gd(4:6,4:6) = TRANSPOSE(b12)

   !-----------------------------
   ! Compute gyroscopic matrix Xd
   !-----------------------------
   Xd(:,:) = 0.0D0
   Xd(4:6,4:6) = MATMUL(TRANSPOSE(Tilde(E1)),Gd(1:3,4:6))

   !-----------------------------
   ! Compute gyroscopic matrix Yd
   !-----------------------------
   Yd(:,:) = 0.0D0
   Yd(4:6,1:3) = b11
   Yd(4:6,4:6) = b12

   END SUBROUTINE DissipativeForce
