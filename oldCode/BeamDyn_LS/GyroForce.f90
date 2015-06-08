   SUBROUTINE GyroForce(mEta,rho,uuu,vvv,Fb)
!----------------------------------------------------------------------------------------
! This subroutine computes gyroscopic forces 
!----------------------------------------------------------------------------------------
   REAL(ReKi),INTENT(IN):: mEta(:) ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: rho(:,:) ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: uuu(:) ! Displacement(and rotation)  array at Gauss point
   REAL(ReKi),INTENT(IN):: vvv(:) ! Velocities at Gauss point (including linear and angular velocities)
   
   REAL(ReKi),INTENT(OUT):: Fb(:) ! Gyroscopic forces

   INTEGER(IntKi):: i,j
!   REAL(ReKi):: H(3,3),Hinv(3,3),HD(3,3)
   REAL(ReKi):: Bi(6,6),ome(3)
!   REAL(ReKi):: cc(3),cd(3)
   REAL(ReKi):: temp33(3,3),temp6(6)
   
   ome = 0.0D0
!   cc = 0.0D0
!   cd = 0.0D0
!   H = 0.0D0
!   Hinv = 0.0D0
!   HD = 0.0D0

   ome(:) = vvv(4:6)
!   cc(:) = uuu(4:6)
!   CALL CrvMatrixH(cc,H)
!   CALL CrvMatrixHinv(cc,Hinv)
!   cd = MATMUL(Hinv,ome)
!   CALL CrvMatrixHD(cc,cd,H,HD)

   Bi= 0.0D0
   temp33 = 0.0D0
!   temp33 = MATMUL(TRANSPOSE(Tilde(mEta)),H)
!   temp33 = MATMUL(Tilde(ome),temp33)
!   temp33 = temp33 + MANUAL(TRANSPOSE(Tilde(mEta)),HD)
   temp33 = MATMUL(Tilde(ome),TRANSPOSE(Tilde(mEta)))
   DO i=1,3
       DO j=1,3
           Bi(i,j+3) = temp33(i,j)
       ENDDO
   ENDDO
   temp33 = 0.0D0   
!   temp33 = MATMUL(rho,HD) + MATMUL(Tilde(ome),MATMUL(rho,H))
   temp33 = MATMUL(Tilde(ome),rho)
   DO i=1,3
       DO j=1,3
           Bi(i+3,j+3) = temp33(i,j)
       ENDDO
   ENDDO

   temp6 = 0.0D0
!   temp6(1:3) = vvv(1:3)
!   temp6(4:6) = cd(1:3)
   temp6(1:6) = vvv(1:6)  

   Fb = 0.0D0
   Fb = MATMUL(Bi,temp6)

   END SUBROUTINE GyroForce
