   SUBROUTINE TiSchmComputeCoefficients(beta,gama,deltat,alfaM,alfaF,coef)

   REAL(ReKi),INTENT(IN)::beta, gama, deltat, alfaM, alfaF

   REAL(ReKi),INTENT(INOUT):: coef(:)

   REAL(ReKi)::deltat2, oalfaM, tr0, tr1, tr2

   deltat2 = deltat * deltat
   oalfaM = 1.0D0 - alfaM
   tr0 =  alfaF / oalfaM
   tr1 = alfaM / oalfaM
   tr2 = (1.0D0 - alfaF) / oalfaM

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
