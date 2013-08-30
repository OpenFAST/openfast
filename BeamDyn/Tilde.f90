!   SUBROUTINE Tilde(vect,tilded)

!   REAL(ReKi),INTENT(IN):: vect(3)
!   REAL(ReKi), INTENT(INOUT):: tilded(3,3)

 !  tilded = 0.d0

!   tilded(1,2) = -vect(3)
!   tilded(1,3) = vect(2) 
!   tilded(2,1) = vect(3)
!   tilded(2,3) = -vect(1)
!   tilded(3,1) = -vect(2)
!   tilded(3,2) = vect(1)

!   END SUBROUTINE Tilde

   FUNCTION Tilde(vect)

   REAL(ReKi),INTENT(IN):: vect(3)
   REAL(ReKi):: Tilde(3,3)

   Tilde = 0.0D0

   Tilde(1,2) = -vect(3)
   Tilde(1,3) = vect(2)
   Tilde(2,1) = vect(3)
   Tilde(2,3) = -vect(1)
   Tilde(3,1) = -vect(2)
   Tilde(3,2) = vect(1)

   END FUNCTION Tilde
