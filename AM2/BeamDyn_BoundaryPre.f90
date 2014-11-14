   SUBROUTINE BeamDyn_BoundaryPre(x,u,t,counter,ErrStat,ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t

   REAL(DbKi),           INTENT(IN   )  :: t           ! Inputs at t

   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   INTEGER(IntKi),               INTENT(IN   )  :: counter     ! Error status of the operation
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                                   :: temp_R(3,3)
   REAL(ReKi)                                   :: temp_rot(3)

   REAL(ReKi)                                   :: temp_pp(3)
   REAL(ReKi)                                   :: temp_qq(3)
   REAL(ReKi)                                   :: pi

   pi = ACOS(-1.0D0)
   temp_pp(:) = 0.0D0
!   temp_pp(2) = 4.0D0*TAN((-pi*t*1.0D0/3.0D0+(2.0D0*pi*counter))/4.0D0)
   temp_qq(:) = 0.0D0
   x%q(4:6) = temp_pp(:)



   END SUBROUTINE BeamDyn_BoundaryPre
