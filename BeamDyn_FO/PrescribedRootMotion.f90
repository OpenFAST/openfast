   SUBROUTINE PrescribedRootMotion(u,x,p)

   TYPE(BDyn_InputType),INTENT(IN):: u
   TYPE(BDyn_ContinuousStateType),INTENT(INOUT):: x
   TYPE(BDyn_ParameterType),INTENT(IN):: p

   INTEGER(IntKi):: i


   DO i=1,p%dof_node
       x%q(i) = 0.0D0
       x%dqdt(i) = 0.0D0
   ENDDO

   DO i=1,3
       x%q(i) = u%PointMesh%TranslationDisp(i,1)
       x%q(i+3) = 0.0D0
       x%dqdt(i) = u%PointMesh%TranslationVel(i,1)
       x%dqdt(i+3) = 0.0D0
   ENDDO

   END SUBROUTINE PrescribedRootMotion
