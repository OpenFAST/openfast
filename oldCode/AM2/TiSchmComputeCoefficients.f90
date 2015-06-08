   SUBROUTINE TiSchmComputeCoefficients(rhoinf,deltat,coef)

   REAL(DbKi),INTENT(IN   ):: rhoinf
   REAL(DbKi),INTENT(IN   ):: deltat
   REAL(DbKi),INTENT(  OUT):: coef(:)
   
   REAL(DbKi)              :: tr0
   REAL(DbKi)              :: tr1
   REAL(DbKi)              :: tr2
   REAL(DbKi)              :: alfam
   REAL(DbKi)              :: alfaf
   REAL(DbKi)              :: gama
   REAL(DbKi)              :: beta
   REAL(DbKi)              :: oalfaM
   REAL(DbKi)              :: deltat2

   tr0 = rhoinf + 1.0D0
   alfam = (2.0D0 * rhoinf - 1.0D0) / tr0
   alfaf = rhoinf / tr0
   gama = 0.5D0 - alfam + alfaf
   beta = 0.25 * (1.0D0 - alfam + alfaf) * (1.0D0 - alfam + alfaf)

   deltat2 = deltat * deltat
   oalfaM = 1.0D0 - alfam
   tr0 = alfaf / oalfaM
   tr1 = alfam / oalfaM
   tr2 = (1.0D0 - alfaf) / oalfaM

   coef(1) = beta * tr0 * deltat2
   coef(2) = (0.5D0 - beta/oalfaM) * deltat2
   coef(3) = gama * tr0 * deltat
   coef(4) = (1.0D0 - gama / oalfaM) * deltat
   coef(5) = tr0
   coef(6) = -tr1
   coef(7) = gama * tr2 * deltat
   coef(8) = beta * tr2 * deltat2
   coef(9) = tr2 

   END SUBROUTINE TiSchmComputeCoefficients