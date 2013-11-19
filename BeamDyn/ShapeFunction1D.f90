   SUBROUTINE ShapeFunction1D(flag,rr,hhh,hhr)
   REAL(ReKi),INTENT(IN)::rr
   REAL(ReKi),INTENT(INOUT):: hhh(:),hhr(:)
   INTEGER(IntKi),INTENT(IN):: flag

   REAL(ReKi):: temp1,temp2,temp3,temp4

   temp1 = 0.1111111111111111D0 !1/9
   temp2 = 0.5625D0 !9/16
   temp3 = 1.6875 !27/16
   temp4 = 0.3333333333333333D0 !1/3

   SELECT CASE (flag)
       CASE (2)
           hhh(1) = 0.5D0 * (1.0D0 - rr)
           hhh(2) = 0.5D0 * (1.0D0 + rr)
           hhr(1) = -0.5D0
           hhr(2) = 0.5D0
       CASE (3)
           hhh(1) = -0.5D0 * rr * (1.0D0 - rr)
           hhh(2) = 1.0D0 - rr * rr
           hhh(3) = 0.5D0 * rr * (1.0D0 + rr)
           hhr(1) = -0.5D0 + rr 
           hhr(2) = -2.0D0 * rr
           hhr(3) = 0.5D0 + rr
       CASE(4)
           hhh(1) = temp2 * (rr * rr - temp1)*(1.0D0 - rr)
           hhh(2) = -temp3 * (1.0D0 - rr * rr)*(rr - temp4)
           hhh(3) = temp3 * (1.0D0 - rr * rr)*(rr + temp4)
           hhh(4) = temp2 * (rr * rr - temp1)*(1.0D0 + rr)
           hhr(1) = temp2 * (2.0D0 * rr - 3.0D0 * rr * rr + temp1)
           hhr(2) = -temp3 * (1.0D0 - 3.0D0 * rr * rr + temp4 * 2.0D0 * rr)
           hhr(3) = temp3 * (1.0D0 - 3.0D0 * rr * rr - temp4 * 2.0D0 * rr)
           hhr(4) = temp2 *(2.0D0 * rr + 3.0D0 * rr * rr - temp1)
   END SELECT
           


   END SUBROUTINE ShapeFunction1D
