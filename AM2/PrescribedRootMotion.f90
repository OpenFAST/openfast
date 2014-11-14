   SUBROUTINE PrescribedRootMotion(time,x,p)
   !-----------------------------------------------------
   ! This subroutine computes the motion of the blade
   !-----------------------------------------------------
!   TYPE(BD_InputType),INTENT(IN):: u
   REAL(DbKi),INTENT(IN):: time ! Time
   TYPE(BD_ContinuousStateType),INTENT(INOUT):: x ! States: q and first time derivative of q
   TYPE(BD_ParameterType),INTENT(IN):: p ! DOF per node

   INTEGER(IntKi):: i


   DO i=1,p%dof_node
       x%q(i) = 0.0D0
       x%dqdt(i) = 0.0D0
   ENDDO

!   DO i=1,3
!       x%q(i) = u%PointMesh%TranslationDisp(i,1)
!       x%q(i+3) = 0.0D0
!       x%dqdt(i) = u%PointMesh%TranslationVel(i,1)
!       x%dqdt(i+3) = 0.0D0
!   ENDDO

   x%q(5) = -4.0D0*TAN((3.1415926D0*time/3.0D0)/4.0D0)
   IF(ABS(x%q(5)) .GT. 4.0D0) THEN
       x%q(5) = -4.0D0*TAN((3.1415926D0*time/3.0D0+2.0D0*3.1415926D0)/4.0D0)
   ENDIF
   x%dqdt(5) = -3.1415926D0/3.0D0

   END SUBROUTINE PrescribedRootMotion
