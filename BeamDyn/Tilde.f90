   FUNCTION Tilde(vect)

   REAL(ReKi),INTENT(IN):: vect(3)
   REAL(ReKi):: Tilde(3,3)

      Tilde = 0.d0

      Tilde(1,2) = -vect(3)
      Tilde(1,3) = vect(2)
      Tilde(2,1) = vect(3)
      Tilde(2,3) = -vect(1)
      Tilde(3,1) = -vect(2)
      Tilde(3,2) = vect(1)

      END FUNCTION Tilde