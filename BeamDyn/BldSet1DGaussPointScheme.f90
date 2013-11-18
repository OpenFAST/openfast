   SUBROUTINE BldSet1DGaussPointScheme(flag,Gpr,Gw)

   REAL(ReKi),INTENT(INOUT):: Gpr(:),Gw(:)
   INTEGER(IntKi),INTENT(IN):: flag

   SELECT CASE (flag)
       CASE (1)
           Gpr(1) = 0.0D0
           Gw(1) = 2.0D0
       CASE (2)
           Gpr(1) = - 0.577350269189626D0
           Gpr(2) = -Gpr(1)
           Gw(1) = 1.0D0
           Gw(2) = 1.0D0
       CASE (3)
           Gpr(1) = - 0.774596669241483D0
           Gpr(2) = 0.0D0
           Gpr(3) = -Gpr(1)
           Gw(1) = 0.555555555555556D0
           Gw(2) = 0.888888888888889D0
           Gw(3) = Gw(1)
   END SELECT
  

   END SUBROUTINE BldSet1DGaussPointScheme
