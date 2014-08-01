   SUBROUTINE MassMatrix(m00,mEta,rho,uuu,Mi)
!----------------------------------------------------------------------------------------
! This subroutine computes the mass matrix. 
!----------------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN):: m00 ! Mass density at Gauss point
   REAL(ReKi),INTENT(IN):: mEta(:) ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: rho(:,:) ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: uuu(:) ! Displacement(and rotation)  arrary at Gauss point
   REAL(ReKi),INTENT(OUT)::Mi(6,6) ! Mass matrix

   REAL(ReKi)::cc(3),H(3,3)
   INTEGER(IntKi)::i

   cc = 0.0D0
   H = 0.0D0

   cc(:) = uuu(4:6)
!   CALL CrvMatrixH(cc,H)
  
   !Mass Matrix
   Mi = 0.0D0
   DO i=1,3
       Mi(i,i) = m00
   ENDDO
   Mi(1:3,4:6) = TRANSPOSE(Tilde(mEta))
   Mi(4:6,1:3) = Tilde(mEta)
   Mi(4:6,4:6) = rho

   
   END SUBROUTINE MassMatrix
