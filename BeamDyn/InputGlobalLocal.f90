   SUBROUTINE InputGlobalLocal(p,u,flag)

   TYPE(BD_ParameterType), INTENT(IN   ):: p 
   TYPE(BD_InputType),     INTENT(INOUT):: u
   INTEGER(IntKi),         INTENT(IN   ):: flag            ! 0: Global to Blade; 
                                                           ! 1: Blade to Global
   REAL(ReKi)                           :: RotTen(3,3) 
   REAL(ReKi)                           :: Pos(3)
   REAL(ReKi)                           :: temp66(6,6)
   REAL(ReKi)                           :: temp6(6)   
   INTEGER(IntKi)                       :: i                             

   Pos(1:3)        = p%GlbPos(:)   
   RotTen(1:3,1:3) = p%GlbRot(:,:)
   IF (flag .EQ. 0) THEN
       ! Transform Root Motion from Global to Local (Blade) frame
       u%RootMotion%TranslationDisp(:,1) = MATMUL(TRANSPOSE(RotTen),u%RootMotion%TranslationDisp(:,1))
       u%RootMotion%Orientation(:,:,1) = MATMUL(u%RootMotion%Orientation(:,:,1),TRANSPOSE(RotTen))
       u%RootMotion%Orientation(:,:,1) = MATMUL(u%RootMotion%Orientation(:,:,1),RotTen)
       u%RootMotion%Orientation(:,:,1) = MATMUL(TRANSPOSE(RotTen),u%RootMotion%Orientation(:,:,1))
       CALL MotionTensor(RotTen,Pos,temp66,1)
       temp6(:) = 0.0D0
       temp6(1:3) = u%RootMotion%TranslationVel(1:3,1)
       temp6(4:6) = u%RootMotion%RotationVel(1:3,1)
       temp6(:) = MATMUL(temp66,temp6)
       u%RootMotion%TranslationVel(1:3,1) = temp6(1:3)
       u%RootMotion%RotationVel(1:3,1) = temp6(4:6)
       temp6(:) = 0.0D0
       temp6(1:3) = u%RootMotion%TranslationAcc(1:3,1)
       temp6(4:6) = u%RootMotion%RotationAcc(1:3,1)
       temp6(:) = MATMUL(temp66,temp6)
       u%RootMotion%TranslationAcc(1:3,1) = temp6(1:3)
       u%RootMotion%RotationAcc(1:3,1) = temp6(4:6)
       ! Transform Applied Forces from Global to Local (Blade) frame
       CALL MotionTensor(RotTen,Pos,temp66,0)
       DO i=1,p%node_total
           temp6(:) = 0.0D0
           temp6(1:3) = u%PointLoad%Force(1:3,i)
           temp6(4:6) = u%PointLoad%Moment(1:3,i)
           temp6(:) = MATMUL(TRANSPOSE(temp66),temp6)
           u%PointLoad%Force(1:3,i) = temp6(1:3)
           u%PointLoad%Moment(1:3,i) = temp6(4:6)
       ENDDO
       
       DO i=1,p%ngp * p%elem_total + 2
           temp6(:) = 0.0D0
           temp6(1:3) = u%DistrLoad%Force(1:3,i)
           temp6(4:6) = u%DistrLoad%Moment(1:3,i)
           temp6(:) = MATMUL(TRANSPOSE(temp66),temp6)
           u%DistrLoad%Force(1:3,i) = temp6(1:3)
           u%DistrLoad%Moment(1:3,i) = temp6(4:6)
       ENDDO
   ELSEIF(flag .EQ. 1) THEN
       
   ENDIF

   END SUBROUTINE InputGlobalLocal
