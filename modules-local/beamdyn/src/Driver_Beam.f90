!..................................................................................................................................
! LICENSING                                                                                                                         
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    Glue is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!  
!    ADD DESCRIPTION
!	
!    References:
!
!
!**********************************************************************************************************************************
PROGRAM MAIN

   USE BeamDyn

   USE BeamDyn_Subs  ! for crv extract routines
   USE BeamDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                     :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                    :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   REAL(DbKi)                         :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                         :: t_initial        ! time at initialization
   REAL(DbKi)                         :: t_final          ! time at simulation end 
   REAL(DbKi)                         :: t_global         ! global-loop time marker
   INTEGER(IntKi)                     :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                     :: n_t_global       ! global-loop time counter
   INTEGER(IntKi)                     :: BD_interp_order     ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InitInputType)           :: BD_InitInput
   TYPE(BD_ParameterType)           :: BD_Parameter
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousStateDeriv
   TYPE(BD_InitOutputType)          :: BD_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   TYPE(BD_OtherStateType)          :: BD_OtherState
   TYPE(BD_InputType),ALLOCATABLE  :: BD_Input(:)
   REAL(DbKi),        ALLOCATABLE  :: BD_InputTimes(:)
   TYPE(BD_OutputType),ALLOCATABLE  :: BD_Output(:)
   REAL(DbKi),ALLOCATABLE           :: BD_OutputTimes(:)

   CHARACTER(256)    :: DvrInputFile


   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   INTEGER(IntKi),PARAMETER:: QiTipDisp = 20
   INTEGER(IntKi),PARAMETER:: QiTipVel = 50
   INTEGER(IntKi),PARAMETER:: QiMidDisp = 21
   INTEGER(IntKi),PARAMETER:: QiMidForce = 22
   INTEGER(IntKi),PARAMETER:: QiMidAcc = 23
   INTEGER(IntKi),PARAMETER:: QiMidVel = 24
   INTEGER(IntKi),PARAMETER:: QiRootUnit = 30
   INTEGER(IntKi),PARAMETER:: QiReacUnit = 40
   INTEGER(IntKi),PARAMETER:: QiRootDispUnit = 50
   INTEGER(IntKi),PARAMETER:: QiRootVelUnit = 60
   INTEGER(IntKi),PARAMETER:: QiRootAccUnit = 70
   INTEGER(IntKi),PARAMETER:: QiInputUnit = 80

   REAL(ReKi):: temp_cc(3)
   REAL(ReKi):: temp_R(3,3)
   REAL(DbKi):: start, finish

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------
   CALL GET_COMMAND_ARGUMENT(1,DvrInputFile)
   CALL BD_ReadDvrFile(DvrInputFile,t_initial,t_final,dt_global,BD_InitInput,&
                          ErrStat,ErrMsg)
   BD_InitInput%RootName  = TRIM(BD_Initinput%InputFile)
   CALL BD_CrvMatrixR(BD_InitInput%GlbRot(:,1),temp_R,ErrStat,ErrMsg)
   BD_InitInput%GlbRot(1:3,1:3) = TRANSPOSE(temp_R(1:3,1:3))
   BD_InitInput%RootDisp(:) = 0.0D+00
   BD_InitInput%RootOri(:,:) = 0.0D0
   BD_InitInput%RootOri(1:3,1:3) = BD_InitInput%GlbRot(1:3,1:3)
   t_global = t_initial
   n_t_final = ((t_final - t_initial) / dt_global )

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BD_interp_order = 0

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1)) 
   ALLOCATE(BD_InputTimes(BD_interp_order + 1)) 
   ALLOCATE(BD_Output(BD_interp_order + 1)) 
   ALLOCATE(BD_OutputTimes(BD_interp_order + 1)) 

    OPEN(unit = QiTipDisp, file = 'Qi_Tip_Disp_Single.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiTipVel, file = 'Qi_Tip_Vel_Single.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiMidDisp, file = 'Qi_Mid_Disp.out', status = 'REPLACE',ACTION = 'WRITE')
!    OPEN(unit = QiMidForce, file = 'Qi_Mid_Force.out', status = 'REPLACE',ACTION = 'WRITE')
!    OPEN(unit = QiMidAcc, file = 'Qi_Mid_Acc_Single.out', status = 'REPLACE',ACTION = 'WRITE')
!    OPEN(unit = QiMidVel, file = 'Qi_Mid_Vel.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRootUnit,file = 'QiRoot_Single.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiReacUnit,file = 'QiReac_Single.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRootDispUnit,file = 'Qi_Root_Disp.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRootVelUnit,file = 'Qi_Root_Vel.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRootAccUnit,file = 'Qi_Root_Acc.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiInputUnit,file = 'BD_Input_File.dat', status = 'REPLACE',ACTION = 'WRITE')


   CALL BD_Init(BD_InitInput        &
                   , BD_Input(1)         &
                   , BD_Parameter        &
                   , BD_ContinuousState  &
                   , BD_DiscreteState    &
                   , BD_ConstraintState  &
                   , BD_OtherState       &
                   , BD_Output(1)        &
                   , dt_global             &
                   , BD_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )


   BD_InputTimes(1) = t_initial 
   BD_OutputTimes(1) = t_initial

   CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), BD_Parameter, BD_InitInput,ErrStat, ErrMsg)
   WRITE(QiInputUnit,9000) BD_InputTimes(1),&
                        &BD_Input(1)%RootMotion%Orientation(:,:,1),&
                        &BD_Input(1)%RootMotion%TranslationDisp(1:3,1),&
                        &BD_Input(1)%RootMotion%TranslationVel(1:3,1),&
                        &BD_Input(1)%RootMotion%RotationVel(1:3,1),&
                        &BD_Input(1)%RootMotion%TranslationAcc(1:3,1),&
                        &BD_Input(1)%RootMotion%RotationAcc(1:3,1)
CALL CPU_TIME(start)
   DO n_t_global = 0, n_t_final
WRITE(*,*) "Time Step: ", n_t_global
!IF(n_t_global == 0) STOP 
     CALL BD_CalcOutput( t_global, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                             BD_ConstraintState, &
                             BD_OtherState,  BD_Output(1), ErrStat, ErrMsg)
     IF(BD_Parameter%analysis_type .EQ. 2 .AND. n_t_global .EQ. 0) THEN
         CALL BD_InitAcc( t_global, BD_Input(1), BD_Parameter, &
               BD_ContinuousState,BD_OtherState,ErrStat,ErrMsg)
      WRITE(*,*) 'Initial Acc'
      WRITE(*,*) BD_OtherState%acc(:)
      WRITE(*,*) 'Initial Xcc'
      WRITE(*,*) BD_OtherState%xcc(:)
     ENDIF


CALL BD_CrvExtractCrv(TRANSPOSE(BD_OutPut(1)%BldMotion%Orientation(1:3,1:3,BD_Parameter%node_elem*BD_Parameter%elem_total)),temp_cc,&
                      ErrStat,ErrMsg)
      WRITE(QiTipDisp,6000) t_global,&
                           &BD_OutPut(1)%BldMotion%TranslationDisp(1:3,BD_Parameter%node_elem*BD_Parameter%elem_total),&
                           &temp_cc(1:3)
      WRITE(QiTipVel,6000) t_global,&
                           &BD_OutPut(1)%BldMotion%TranslationVel(1:3,BD_Parameter%node_elem*BD_Parameter%elem_total),&
                           &BD_OutPut(1)%BldMotion%RotationVel(1:3,BD_Parameter%node_elem*BD_Parameter%elem_total)
      WRITE(QiRootVelUnit,6000) t_global,&
                           &BD_OutPut(1)%BldMotion%TranslationVel(1:3,1),&
                           &BD_OutPut(1)%BldMotion%RotationVel(1:3,1)
      WRITE(QiRootAccUnit,6000) t_global,&
                           &BD_OutPut(1)%BldMotion%TranslationAcc(1:3,1),&
                           &BD_OutPut(1)%BldMotion%RotationAcc(1:3,1)
      WRITE(QiRootUnit,6000) t_global,&
                           &BD_OutPut(1)%BldForce%Force(1:3,1),&
                           &BD_OutPut(1)%BldForce%Moment(1:3,1)
      WRITE(QiReacUnit,6000) t_global,&
                           &BD_OutPut(1)%ReactionForce%Force(1:3,1),&
                           &BD_OutPut(1)%ReactionForce%Moment(1:3,1)
CALL BD_CrvExtractCrv(TRANSPOSE(BD_OutPut(1)%BldMotion%Orientation(1:3,1:3,1)),temp_cc,&
                      ErrStat,ErrMsg)
      WRITE(QiMidDisp,6000) t_global,&
                           &BD_OutPut(1)%BldMotion%TranslationDisp(1:3,1),&
                           &temp_cc(1:3)
!                           &BD_OutPut(1)%BldMotion%TranslationVel(1:3,1)
!      WRITE(QiMidAcc,6000) t_global,&
!                           &BD_OutPut(1)%BldMotion%TranslationAcc(1:3,BD_Parameter%node_elem),&
!                           &BD_OutPut(1)%BldMotion%RotationAcc(1:3,BD_Parameter%node_elem)

     IF(BD_Parameter%analysis_type .EQ. 1 .AND. n_t_global .EQ. 1) EXIT !RETURN

     CALL BD_InputSolve( t_global + dt_global, BD_Input(1), BD_Parameter, BD_InitInput, ErrStat, ErrMsg)
     WRITE(QiInputUnit,9000) t_global+dt_global,&
                        &BD_Input(1)%RootMotion%Orientation(:,:,1),&
                        &BD_Input(1)%RootMotion%TranslationDisp(1:3,1),&
                        &BD_Input(1)%RootMotion%TranslationVel(1:3,1),&
                        &BD_Input(1)%RootMotion%RotationVel(1:3,1),&
                        &BD_Input(1)%RootMotion%TranslationAcc(1:3,1),&
                        &BD_Input(1)%RootMotion%RotationAcc(1:3,1)

     CALL BD_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, &
                               BD_ContinuousState, &
                               BD_DiscreteState, BD_ConstraintState, &
                               BD_OtherState, ErrStat, ErrMsg )

!WRITE(*,*) 'x%q'
!WRITE(*,*) BD_ContinuousState%q
!WRITE(*,*) 'x%dqdt'
!WRITE(*,*) BD_ContinuousState%dqdt
!WRITE(*,*) 'OtherState%acc'
!WRITE(*,*) BD_OtherState%acc
!WRITE(*,*) 'OtherState%xcc'
!WRITE(*,*) BD_OtherState%xcc


      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

   ENDDO
CALL CPU_TIME(finish)
!print '("Time = ",f6.3," seconds.")',finish-start
WRITE(*,*) 'Start: ', start
WRITE(*,*) 'Finish: ', finish
WRITE(*,*) 'Time: ', finish-start

   

   CALL BD_End( BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                    BD_ConstraintState, BD_OtherState, BD_Output(1), ErrStat, ErrMsg )


   DEALLOCATE(BD_Input)
   DEALLOCATE(BD_InputTimes)
   DEALLOCATE(BD_Output)
   DEALLOCATE(BD_OutputTimes)

   6000 FORMAT (ES12.5,6ES21.12)
   9000 FORMAT (ES12.5,24ES21.12)
   CLOSE (QiTipDisp)
   CLOSE (QiTipVel)
   CLOSE (QiRootUnit)
   CLOSE (QiReacUnit)
   CLOSE (QiMidDisp)
!   CLOSE (QiMidForce)
!   CLOSE (QiMidAcc)
   CLOSE (QiRootDispUnit)
   CLOSE (QiRootVelUnit)
   CLOSE (QiRootAccUnit)
   CLOSE (QiInputUnit)

7000 FORMAT (ES12.5,9ES21.12)
!CLOSE (QiHUnit)

END PROGRAM MAIN



SUBROUTINE BD_InputSolve( t, u,  p, InitInput, ErrStat, ErrMsg)
 
   USE BeamDyn
   USE BeamDyn_Subs
   USE BeamDyn_Types

   REAL(DbKi),                     INTENT(IN   ):: t
   TYPE(BD_InputType),             INTENT(INOUT):: u
   TYPE(BD_ParameterType),         INTENT(IN   ):: p
   TYPE(BD_InitInputType),         INTENT(IN   ):: InitInput
   INTEGER(IntKi),                 INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)          :: i                ! do-loop counter
   REAL(ReKi)              :: temp_vec(3)
   REAL(ReKi)              :: temp_vec2(3)
   REAL(ReKi)              :: temp_rr(3)
   REAL(ReKi)              :: temp_R(3,3)
   REAL(ReKi)              :: temp_r0(3)
   REAL(ReKi)              :: temp_theta(3)

   ErrStat = ErrID_None
   ErrMsg  = ''

   temp_r0(:) = 0.0D0 
   temp_rr(:)     = 0.0D0
   temp_R(:,:)    = 0.0D0

   temp_r0(1) = p%GlbPos(2)
   temp_r0(2) = p%GlbPos(3)
   temp_r0(3) = p%GlbPos(1)

   temp_theta(1) = p%IniVelo(5)*t
   temp_theta(2) = p%IniVelo(6)*t
   temp_theta(3) = p%IniVelo(4)*t
   temp_vec(:) = 0.0D0
   temp_vec(:) = 4.0D0*TAN(temp_theta(:)/4.0D0)
   ! gather point forces and line forces

!------------------
!  Rotating beam
!------------------
   ! Point mesh: RootMotion 
   ! Calculate root displacements and rotations
   u%RootMotion%Orientation(:,:,:) = 0.0D0
   DO i=1,3
       u%RootMotion%Orientation(i,i,1) = 1.0D0
   ENDDO
   CALL BD_CrvMatrixR(temp_vec,u%RootMotion%Orientation(:,:,1),ErrStat,ErrMsg)
   temp_rr(:) = MATMUL(u%RootMotion%Orientation(:,:,1),temp_r0)
   u%RootMotion%Orientation(:,:,1) = TRANSPOSE(u%RootMotion%Orientation(:,:,1))
   u%RootMotion%TranslationDisp(:,:)  = 0.0D0
   u%RootMotion%TranslationDisp(:,1) = temp_rr(:) - temp_r0(:)
   ! END Calculate root displacements and rotations

   ! Calculate root translational and angular velocities
   u%RootMotion%RotationVel(:,:) = 0.0D0
   u%RootMotion%RotationVel(1,1) = p%IniVelo(5)
   u%RootMotion%RotationVel(2,1) = p%IniVelo(6)
   u%RootMotion%RotationVel(3,1) = p%IniVelo(1)
   u%RootMotion%TranslationVel(:,:) = 0.0D0
   u%RootMotion%TranslationVel(:,1) = MATMUL(BD_Tilde(u%RootMotion%RotationVel(:,1)),temp_rr)
   ! END Calculate root translational and angular velocities


   ! Calculate root translational and angular accelerations
   u%RootMotion%TranslationAcc(:,:) = 0.0D0
   u%RootMotion%RotationAcc(:,:) = 0.0D0
   u%RootMotion%TranslationAcc(:,1) = MATMUL(BD_Tilde(u%RootMotion%RotationVel(:,1)), &
               MATMUL(BD_Tilde(u%RootMotion%RotationVel(:,1)),temp_rr))
   ! END Calculate root translational and angular accelerations
!------------------
! End rotating beam
!------------------

   u%PointLoad%Force(:,:)  = 0.0D0
   u%PointLoad%Moment(:,:) = 0.0D0
!   u%PointLoad%Force(1,p%node_total)  = 5.0D+04
   u%PointLoad%Force(1:3,p%node_total)  = InitInput%TipLoad(1:3)
   u%PointLoad%Moment(1:3,p%node_total) = InitInput%TipLoad(4:6)
   
!   u%RootMotion%TranslationAcc(3,1) = 1.76991032448401212D-02
!   u%PointLoad%Force(1,p%node_total) = 1.0D+02*SIN(2.0*t)
!   u%PointLoad%Force(1,p%node_total) = 2.0
!   u%PointLoad%Force(3,p%node_total) = 1.0D+03*0.5*(1.0D0-COS(0.2*t))
!   u%PointLoad%Force(3,3) = 1.0D+00

   ! LINE2 mesh: DistrLoad
   u%DistrLoad%Force(:,:)  = 0.0D0
   u%DistrLoad%Moment(:,:) = 0.0D0

   IF(p%quadrature .EQ. 1) THEN
       DO i=1,p%ngp(1)*p%elem_total+2
           u%DistrLoad%Force(1:3,i) = InitInput%DistrLoad(1:3)
           u%DistrLoad%Moment(1:3,i)= InitInput%DistrLoad(4:6)
       ENDDO
   ELSEIF(p%quadrature .EQ. 2) THEN
       DO i=1,p%kp_total
           u%DistrLoad%Force(1:3,i) = InitInput%DistrLoad(1:3)
           u%DistrLoad%Moment(1:3,i)= InitInput%DistrLoad(4:6)
       ENDDO
   ENDIF

END SUBROUTINE BD_InputSolve

SUBROUTINE BD_ReadDvrFile(DvrInputFile,t_ini,t_f,dt,InitInputData,&
                          ErrStat,ErrMsg)
!------------------------------------------------------------------------------------
! This routine reads in the primary BeamDyn input file and places the values it reads
! in the InputFileData structure.
!   It opens an echo file if requested and returns the (still-open) echo file to the
!     calling routine.
!   It also returns the names of the BldFile, FurlFile, and TrwFile for further
!     reading of inputs.
!------------------------------------------------------------------------------------
   USE BeamDyn
   USE BeamDyn_Types
   USE NWTC_Library

   ! Passed variables
   CHARACTER(*),                 INTENT(IN   ) :: DvrInputFile
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg
   TYPE(BD_InitInputType),       INTENT(  OUT) :: InitInputData
   REAL(DbKi),                   INTENT(  OUT) :: t_ini
   REAL(DbKi),                   INTENT(  OUT) :: t_f
   REAL(DbKi),                   INTENT(  OUT) :: dt

   ! Local variables:
   INTEGER(IntKi)               :: UnIn                         ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
   character(*), parameter      :: RoutineName = 'BD_ReadDvrFile'
   INTEGER(IntKi)               :: UnEc
   
   CHARACTER(1024)              :: FTitle                       ! "File Title": the 2nd line of the input file, which contains a description of its contents

   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j

   ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   CALL GetNewUnit(UnIn,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL OpenFInpFile(UnIn,DvrInputFile,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
   !-------------------------- HEADER ---------------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'File Header: Module Version (line 1)',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL ReadStr(UnIn,DvrInputFile,FTitle,'FTitle','File Header: File Description (line 2)',ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   !---------------------- SIMULATION CONTROL --------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Simulation Control',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   t_ini = 0.0D0 
   t_f   = 0.0D0 
   dt    = 0.0D0 
   CALL ReadVar(UnIn,DvrInputFile,t_ini,'t_initial','Starting time of simulation',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar(UnIn,DvrInputFile,t_f,"t_final", "Ending time of simulation",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL ReadVar(UnIn,DvrInputFile,dt,"dt", "Time increment size",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !---------------------- GRAVITY PARAMETER --------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Gravity Parameter',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   InitInputData%gravity(:) = 0.0D0
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%gravity(1),"InitInputData%gravity(1)", "gravity vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%gravity(2),"InitInputData%gravity(2)", "gravity vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%gravity(3),"InitInputData%gravity(3)", "gravity vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call cleanup()
       return
   end if
   !---------------------- FRAME PARAMETER --------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Frame Parameter',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   InitInputData%GlbPos(:)   = 0.0D0
   InitInputData%GlbRot(:,:) = 0.0D0
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbPos(1),"InitInputData%GlbPos(1)", "position vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbPos(2),"InitInputData%GlbPos(2)", "position vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbPos(3),"InitInputData%GlbPos(3)", "position vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbRot(1,1),"InitInputData%GlbPos(1,1)", "rotation angle X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbRot(2,1),"InitInputData%GlbPos(2,1)", "rotation angle Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbRot(3,1),"InitInputData%GlbPos(3,1)", "rotation angle Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call cleanup()
       return
   end if
   !---------------------- INITIAL VELOCITY PARAMETER --------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Initial Velocity Parameter',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   InitInputData%RootVel(:)   = 0.0D0
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(1),"InitInputData%IniRootVel(1)", "velocity vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(2),"InitInputData%IniRootVel(2)", "velocity vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(3),"InitInputData%IniRootVel(3)", "velocity vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(4),"InitInputData%IniRootVel(4)", "angular velocity vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(5),"InitInputData%IniRootVel(5)", "angular velocity vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(6),"InitInputData%IniRootVel(6)", "angular velocity vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   if (ErrStat >= AbortErrLev) then
       call cleanup()
       return
   end if
  
   !---------------------- APPLIED FORCE --------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Applied Force',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   InitInputData%DistrLoad(:)   = 0.0D0
   InitInputData%TipLoad(:)     = 0.0D0
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%DistrLoad(1),"InitInputData%DistrLoad(1)", "Distributed load vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%DistrLoad(2),"InitInputData%DistrLoad(2)", "Distributed load vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%DistrLoad(3),"InitInputData%DistrLoad(3)", "Distributed load vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%DistrLoad(4),"InitInputData%DistrLoad(4)", "Distributed load vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%DistrLoad(5),"InitInputData%DistrLoad(5)", "Distributed load vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%DistrLoad(6),"InitInputData%DistrLoad(6)", "Distributed load vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%TipLoad(1),"InitInputData%TipLoad(1)", "Tip load vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%TipLoad(2),"InitInputData%TipLoad(2)", "Tip load vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%TipLoad(3),"InitInputData%TipLoad(3)", "Tip load vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%TipLoad(4),"InitInputData%TipLoad(4)", "Tip load vector X",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%TipLoad(5),"InitInputData%TipLoad(5)", "Tip load vector Y",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%TipLoad(6),"InitInputData%TipLoad(6)", "Tip load vector Z",ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   if (ErrStat >= AbortErrLev) then
       call cleanup()
       return
   end if
   !---------------------- BEAM SECTIONAL PARAMETER ----------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Primary input file',ErrStat2,ErrMsg2,UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadVar ( UnIn, DvrInputFile, InitInputData%InputFile, 'InputFile', 'Name of the primary input file', ErrStat2,ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call cleanup()
       return
   end if

   call cleanup()
   return
      
contains
   subroutine cleanup() 
      close(UnIn)
      return
   end subroutine cleanup         
END SUBROUTINE BD_ReadDvrFile
