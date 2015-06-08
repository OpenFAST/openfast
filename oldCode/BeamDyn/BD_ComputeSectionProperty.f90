   SUBROUTINE BD_ComputeSectionProperty(flp,edg,ang,qflp,qedg,qcrs)

   REAL(ReKi),INTENT(IN   ) ::flp
   REAL(ReKi),INTENT(IN   ) ::edg
   REAL(ReKi),INTENT(IN   ) ::ang
   REAL(ReKi),INTENT(  OUT) ::qflp
   REAL(ReKi),INTENT(  OUT) ::qedg
   REAL(ReKi),INTENT(  OUT) ::qcrs

   REAL(ReKi)               ::temp

   temp = 0.0D0
   qflp = 0.0D0
   qedg = 0.0D0
   qcrs = 0.0D0  

   temp = SIN(2.0D0*ang*ACOS(-1.0D0)/180.0D0)
   qedg = 0.5D0*(flp+SQRT((1.0D0-temp*temp)*(flp-edg)*(flp-edg))+edg)
   qflp = 0.5D0*(flp-SQRT((1.0D0-temp*temp)*(flp-edg)*(flp-edg))+edg)
   qcrs = 0.5D0*temp*(edg-flp)

   END SUBROUTINE BD_ComputeSectionProperty
