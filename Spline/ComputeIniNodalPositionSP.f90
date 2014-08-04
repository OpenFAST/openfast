   SUBROUTINE ComputeIniNodalPositionSP(Coef,eta,PosiVec,e1)

   REAL(ReKi),INTENT(IN   ):: Coef(:,:)
   REAL(ReKi),INTENT(IN   ):: eta
   REAL(ReKi),INTENT(  OUT):: PosiVec(:)
   REAL(ReKi),INTENT(  OUT):: e1(:)

   REAL(ReKi)              :: temp
   INTEGER(IntKi)          :: i

   temp = (eta + 1.0D0)/2.0D0
   DO i=1,3
       PosiVec(i) = Coef(i,1) + Coef(i,2)*temp + Coef(i,3)*temp*temp + Coef(i,4)*temp*temp*temp 
       e1(i) = Coef(i,2) + 2.0D0*Coef(i,3)*temp + 3.0D0*Coef(i,4)*temp*temp
   ENDDO
   e1(:) = e1(:)/Norm(e1)


   END SUBROUTINE ComputeIniNodalPositionSP
