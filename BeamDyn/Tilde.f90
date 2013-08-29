   SUBROUTINE Tilde(vect,tilded)

   DOUBLE PRECISION,INTENT(IN):: vect(3)
   DOUBLE PRECISION, INTENT(INOUT):: tilded(3,3)

   tilded = 0.d0

   tilded(1,2) = -vect(3)
   tilded(1,3) = vect(2) 
   tilded(2,1) = vect(3)
   tilded(2,3) = -vect(1)
   tilded(3,1) = -vect(2)
   tilded(3,2) = vect(1)

   END SUBROUTINE Tilde
