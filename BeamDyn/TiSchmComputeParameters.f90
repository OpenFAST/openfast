   SUBROUTINE TiSchmComputeParameters(rhoinf, alfam,alfaf,gama,beta)

   REAL(ReKi),INTENT(IN)::rhoinf
   REAL(ReKi),INTENT(OUT)::alfam, alfaf, gama, beta

   REAL(ReKi)::tr0


   tr0 = rhoinf + 1.0D0
   alfam = (2.0D0 * rhoinf - 1.0D0) / tr0
   alfaf = rhoinf / tr0
   gama = 0.5D0 - alfam + alfaf
   beta = 0.25 * (1.0D0 - alfam + alfaf) * (1.0D0 - alfam + alfaf)
   

   END SUBROUTINE TiSchmComputeParameters
