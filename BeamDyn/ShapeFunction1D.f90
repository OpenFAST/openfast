   SUBROUTINE ShapeFunction1D(flag,rr,hhh,hhr)
   REAL(ReKi),INTENT(IN)::rr
   REAL(ReKi),INTENT(INOUT):: hhh(:),hhr(:)
   INTEGER(IntKi),INTENT(IN):: flag

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
           hhr(3) = 0.5D0 * rr
       CASE(4)
           hhh(1) = (9.0D0/16.0D0) * (rr * rr - 1.0D0/9.0D0)*(1.0D0 - rr)
           hhh(2) = -(27.0D0/16.0D0) * (1.0D0 - rr * rr)*(rr - 1.0D0/3.0D0)
           hhh(3) = (27.0D0/16.0D0) * (1.0D0 - rr * rr)*(rr + 1.0D0/3.0D0)
           hhh(4) = (9.0D0/16.0D0) * (rr * rr - 1.0D0/9.0D0)*(1.0D0 + rr)
           hhr(1) = (9.0D0/16.0D0) * (2.0D0 * rr - 3.0D0 * rr * rr + 1.0D0/9.0D0)
           hhr(2) = -(27.0D0/16.0D0) * (rr - 1.0D0/3.0D0 - rr*rr*rr + rr*rr/3.0D0)
           hhr(3) = (27.0D0/16.0D0) * (rr + 1.0D0/3.0D0 -rr*rr*rr - rr*rr/3.0D0)
           hhr(4) = (9.0D0/16.0D0) *(2.0D0 * rr + 3.0D0 * rr * rr - 1.0D0/9.0D0)
   END SELECT
           


   END SUBROUTINE ShapeFunction1D
