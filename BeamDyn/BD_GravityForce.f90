   SUBROUTINE BD_GravityForce(m00,mEta,grav,Fg)

   REAL(ReKi),INTENT(IN   ):: m00
   REAL(ReKi),INTENT(IN   ):: mEta(:)
   REAL(ReKi),INTENT(IN   ):: grav(:)
   REAL(ReKi),INTENT(  OUT):: Fg(:)

   Fg = 0.0D0
   Fg(1:3) = m00 * grav(1:3)
   Fg(4:6) = MATMUL(Tilde(mEta),grav)
   

   END SUBROUTINE BD_GravityForce
