   SUBROUTINE CrvMatrixHinv(cc,Hinv)

   REAL(ReKi),INTENT(IN):: cc(:)

   REAL(ReKi),INTENT(OUT):: Hinv(:,:)

   REAL(ReKi):: Htemp(3,3)
   REAL(ReKi):: c0, tr0

   c0 = 0.0D0
   tr0 = 0.0D0
   c0 = 2.0D0 - (cc(1)*cc(1) + cc(2)*cc(2) + cc(3)*cc(3))/8.0D0
   tr0 = (4.0D0 - c0)*(4.0D0 - c0)/4.0D0

   Htemp = 0.0D0
   CALL CrvMatrixh(cc,Htemp)

   Hinv = 0.0D0
   Hinv = tr0*TRANSPOSE(Htemp)


   END SUBROUTINE CrvMatrixHinv
