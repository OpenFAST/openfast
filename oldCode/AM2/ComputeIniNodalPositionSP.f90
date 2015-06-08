   SUBROUTINE ComputeIniNodalPositionSP(Coef,eta,PosiVec,e1,Twist_Angle)

   REAL(ReKi),INTENT(IN   ):: Coef(:,:)
   REAL(ReKi),INTENT(IN   ):: eta
   REAL(ReKi),INTENT(  OUT):: PosiVec(:)
   REAL(ReKi),INTENT(  OUT):: e1(:)
   REAL(ReKi),INTENT(  OUT):: Twist_Angle

   REAL(ReKi)              :: temp
   INTEGER(IntKi)          :: i

   temp = eta 
   PosiVec(:) = 0.0D0
   e1(:) = 0.0D0
   DO i=1,3
       PosiVec(i) = Coef(i,1) + Coef(i,2)*temp + Coef(i,3)*temp*temp + Coef(i,4)*temp*temp*temp 
!       WRITE(*,*) "Coef",i,Coef(i,1:4)
       e1(i) = Coef(i,2) + 2.0D0*Coef(i,3)*temp + 3.0D0*Coef(i,4)*temp*temp
   ENDDO
!   WRITE(*,*) "e1",e1(:)
   e1(:) = e1(:)/Norm(e1)
   Twist_Angle = 0.0D0
   Twist_Angle = Coef(4,1) + Coef(4,2)*temp + Coef(4,3)*temp*temp + Coef(4,4)*temp*temp*temp 
!   WRITE(*,*) "e1",e1(:)


   END SUBROUTINE ComputeIniNodalPositionSP
