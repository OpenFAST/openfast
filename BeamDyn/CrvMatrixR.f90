   SUBROUTINE CrvMatrixR(cc,Rr) 

      REAL(ReKi),INTENT(IN)::cc(:)
      REAL(ReKi),INTENT(OUT)::Rr(:,:)

      INTEGER:: i, j
      REAL(ReKi):: c1,c2,c3,c0,tr0

      Rr = ZERO
      
      c1 = cc(1)/FOUR
      c2 = cc(2)/FOUR
      c3 = cc(3)/FOUR

      c0 = HALF*(ONE-c1*c1-c2*c2-c3*c3)
      tr0 = ONE - c0
      tr0 = TWO/(tr0*tr0)

      Rr(1,1) = tr0*(c1*c1 + c0*c0) - ONE 
      Rr(2,1) = tr0*(c1*c2 + c0*c3)
      Rr(3,1) = tr0*(c1*c3 - c0*c2)
      Rr(1,2) = tr0*(c1*c2 - c0*c3)
      Rr(2,2) = tr0*(c2*c2 + c0*c0) - ONE 
      Rr(3,2) = tr0*(c2*c3 + c0*c1)
      Rr(1,3) = tr0*(c1*c3 + c0*c2)
      Rr(2,3) = tr0*(c2*c3 - c0*c1)
      Rr(3,3) = tr0*(c3*c3 + c0*c0) - ONE

      END SUBROUTINE CrvMatrixR