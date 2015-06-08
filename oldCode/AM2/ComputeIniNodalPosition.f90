   SUBROUTINE ComputeIniNodalPosition(EndP1,EndP2,MidP,eta,PosiVec)

   REAL(ReKi),INTENT(IN):: EndP1(:)
   REAL(ReKi),INTENT(IN):: EndP2(:)
   REAL(ReKi),INTENT(IN):: MidP(:)
   REAL(ReKi),INTENT(IN):: eta
   REAL(ReKi),INTENT(OUT):: PosiVec(:)

   REAL(ReKi):: hhx(3)
   INTEGER(IntKi):: i

   hhx = 0.0D0
   hhx(1) = 0.5D0*(eta*eta - eta)
   hhx(2) = 1.0D0 - eta*eta
   hhx(3) = 0.5D0*(eta*eta + eta)
  
   PosiVec = 0.0D0
   DO i=1,3
       PosiVec(i) = hhx(1)*EndP1(i) + hhx(2)*MidP(i) + hhx(3)*EndP2(i) 
   ENDDO
   
   END SUBROUTINE ComputeIniNodalPosition
