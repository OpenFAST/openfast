!   SUBROUTINE BD_Tilde(vect,BD_Tilded)

!   REAL(ReKi),INTENT(IN):: vect(3)
!   REAL(ReKi), INTENT(INOUT):: BD_Tilded(3,3)

 !  BD_Tilded = 0.d0

!   BD_Tilded(1,2) = -vect(3)
!   BD_Tilded(1,3) = vect(2) 
!   BD_Tilded(2,1) = vect(3)
!   BD_Tilded(2,3) = -vect(1)
!   BD_Tilded(3,1) = -vect(2)
!   BD_Tilded(3,2) = vect(1)

!   END SUBROUTINE BD_Tilde

   FUNCTION BD_Tilde(vect)

   REAL(ReKi),INTENT(IN):: vect(3)
   REAL(ReKi):: BD_Tilde(3,3)

   BD_Tilde = 0.0D0

   BD_Tilde(1,2) = -vect(3)
   BD_Tilde(1,3) = vect(2)
   BD_Tilde(2,1) = vect(3)
   BD_Tilde(2,3) = -vect(1)
   BD_Tilde(3,1) = -vect(2)
   BD_Tilde(3,2) = vect(1)

   END FUNCTION BD_Tilde
