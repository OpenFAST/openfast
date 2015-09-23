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
   INTEGER(IntKi)                     :: DvrOut 

   CHARACTER(256)    :: DvrInputFile
   CHARACTER(256)    :: RootName


   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   REAL(ReKi):: temp_R(3,3)
   REAL(DbKi):: start, finish

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------
   CALL GET_COMMAND_ARGUMENT(1,DvrInputFile)
   CALL GetRoot(DvrInputFile,RootName)
   CALL BD_ReadDvrFile(DvrInputFile,t_initial,t_final,dt_global,BD_InitInput,&
                          ErrStat,ErrMsg)
   BD_InitInput%RootName  = TRIM(BD_Initinput%InputFile)
   BD_InitInput%RootDisp(:) = 0.0D+00
   BD_InitInput%RootOri(:,:) = 0.0D0
   BD_InitInput%RootOri(1:3,1:3) = BD_InitInput%GlbRot(1:3,1:3)
   t_global = t_initial
   n_t_final = ((t_final - t_initial) / dt_global )

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BD_interp_order = 1

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1)) 
   ALLOCATE(BD_InputTimes(BD_interp_order + 1)) 
   ALLOCATE(BD_Output(BD_interp_order + 1)) 
   ALLOCATE(BD_OutputTimes(BD_interp_order + 1)) 

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
      CALL CheckError()

   CALL Dvr_InitializeOutputFile(DvrOut,BD_InitOutput,RootName,ErrStat,ErrMsg)
      CALL CheckError()

   BD_InputTimes(1) = t_initial
   BD_InputTimes(2) = t_initial 
   BD_OutputTimes(1) = t_initial
   BD_OutputTimes(2) = t_initial


   CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), BD_Parameter, BD_InitInput,ErrStat, ErrMsg)
   CALL BD_CopyInput(BD_Input(1), BD_Input(2), MESH_NEWCOPY, ErrStat, ErrMsg)
   CALL BD_CopyOutput(BD_Output(1), BD_Output(2), MESH_NEWCOPY, ErrStat, ErrMsg)
      CALL CheckError()

   CALL CPU_TIME(start)

   DO n_t_global = 0, n_t_final
     WRITE(*,*) "Time Step: ", n_t_global
     BD_InputTimes(2) = BD_InputTimes(1) 
     BD_InputTimes(1) = t_global + dt_global
     BD_OutputTimes(2) = BD_OutputTimes(1) 
     BD_OutputTimes(1) = t_global + dt_global
     CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), BD_Parameter, BD_InitInput, ErrStat, ErrMsg)
     CALL BD_InputSolve( BD_InputTimes(2), BD_Input(2), BD_Parameter, BD_InitInput, ErrStat, ErrMsg)
        CALL CheckError()

     CALL BD_CalcOutput( t_global, BD_Input(2), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                             BD_ConstraintState, &
                             BD_OtherState,  BD_Output(2), ErrStat, ErrMsg)
        CALL CheckError()

     CALL Dvr_WriteOutputLine(t_global,DvrOut,BD_Parameter%OutFmt,BD_Output(2),ErrStat,ErrMsg)
        CALL CheckError()

     IF(BD_Parameter%analysis_type .EQ. 1 .AND. n_t_global .EQ. 1) EXIT 

     CALL BD_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, &
                               BD_ContinuousState, &
                               BD_DiscreteState, BD_ConstraintState, &
                               BD_OtherState, ErrStat, ErrMsg )
        CALL CheckError()

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

   ENDDO

   CALL CPU_TIME(finish)
   
   WRITE(*,*) 'Start: ', start
   WRITE(*,*) 'Finish: ', finish
   WRITE(*,*) 'Time: ', finish-start

   CALL Dvr_End()

CONTAINS

   SUBROUTINE Dvr_End()

      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /=
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      character(*), parameter                       :: RoutineName = 'Dvr_End'

      IF(DvrOut >0) CLOSE(DvrOut)

      DO i=1,BD_interp_order + 1
          CALL BD_End( BD_Input(i), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                           BD_ConstraintState, BD_OtherState, BD_Output(i), ErrStat2, ErrMsg2 )
      ENDDO 
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      IF(ALLOCATED(BD_InputTimes )) DEALLOCATE(BD_InputTimes )
      IF(ALLOCATED(BD_OutputTimes)) DEALLOCATE(BD_OutputTimes)

      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'BeamDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
   END SUBROUTINE Dvr_End

   subroutine CheckError()
   
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            call Dvr_End()
         end if
      end if
         
   end subroutine CheckError

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

   u%PointLoad%Force(:,:)  = 0.0D0
   u%PointLoad%Moment(:,:) = 0.0D0
   u%PointLoad%Force(1:3,p%node_total)  = InitInput%TipLoad(1:3)
   u%PointLoad%Moment(1:3,p%node_total) = InitInput%TipLoad(4:6)
   
   ! LINE2 mesh: DistrLoad
   u%DistrLoad%Force(:,:)  = 0.0D0
   u%DistrLoad%Moment(:,:) = 0.0D0

   IF(p%quadrature .EQ. 1) THEN
       DO i=1,p%ngp*p%elem_total+2
           u%DistrLoad%Force(1:3,i) = InitInput%DistrLoad(1:3)
           u%DistrLoad%Moment(1:3,i)= InitInput%DistrLoad(4:6)
       ENDDO
   ELSEIF(p%quadrature .EQ. 2) THEN
       DO i=1,p%ngp
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
   UnEc = -1
   
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
   DO i=1,3
       CALL ReadAry(UnIn,DvrInputFile,InitInputData%GlbRot(i,:),3,"InitInputData%GlbPos",&
               "Global DCM",ErrStat2,ErrMsg2,UnEc)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ENDDO
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

SUBROUTINE Dvr_InitializeOutputFile(OutUnit,IntOutput,RootName,ErrStat,ErrMsg)

   USE BeamDyn
   USE BeamDyn_Types
   USE NWTC_Library

   INTEGER(IntKi),                    INTENT(  OUT):: OutUnit
   TYPE(BD_InitOutputType),           INTENT(IN   ):: IntOutput     ! Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),                      INTENT(IN   ):: RootName

   integer(IntKi)                            ::  i      
   integer(IntKi)                            :: numOuts
   INTEGER(IntKi)          :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)    :: ErrMsg2                      ! Temporary Error message
   character(*), parameter :: RoutineName = 'Dvr_InitializeOutputFile'

   CALL GetNewUnit(OutUnit,ErrStat2,ErrMsg2)
   CALL OpenFOutFile ( OutUnit, trim(RootName)//'.out', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   write (OutUnit,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using'//trim(GetNVD(IntOutput%Ver))
   write (OutUnit,'()' )    !print a blank line
   
   numOuts = size(IntOutput%WriteOutputHdr)
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................

   write (OutUnit,'()')
   write (OutUnit,'()')
   write (OutUnit,'()')

   call WrFileNR ( OutUnit, 'Time        ' )

   do i=1,NumOuts
      call WrFileNR ( OutUnit, tab//IntOutput%WriteOutputHdr(i) )
   end do ! i

   write (OutUnit,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................

   call WrFileNR ( OutUnit, ' (s)            ' )

   do i=1,NumOuts
      call WrFileNR ( Outunit, tab//IntOutput%WriteOutputUnt(i) )
   end do ! i

   write (OutUnit,'()')  


END SUBROUTINE Dvr_InitializeOutputFile


SUBROUTINE Dvr_WriteOutputLine(t,OutUnit, OutFmt, Output, errStat, errMsg)

   USE BeamDyn_Types
   USE NWTC_Library

   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   INTEGER(IntKi)         ,  intent(in   )   :: OutUnit              ! Status of error message
   CHARACTER(*)           ,  intent(in   )   :: OutFmt
!   real(ReKi)             ,  intent(in   )   :: output(:)            ! Rootname for the output file
   TYPE(BD_OutputType),      INTENT(IN   )   :: Output
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      
   ! Local variables.

   character(200)                   :: frmt                                      ! A string to hold a format specifier
   character(15)                    :: tmpStr                                    ! temporary string to print the time output as text

   integer :: numOuts
   
   errStat = ErrID_None
   errMsg  = ''

   numOuts = size(Output%WriteOutput,1)
   frmt = '"'//tab//'"'//trim(OutFmt)      ! format for array elements from individual modules
   
      ! time
   write( tmpStr, '(F15.6)' ) t
   call WrFileNR( OutUnit, tmpStr )
   call WrNumAryFileNR ( OutUnit, Output%WriteOutput,  frmt, errStat, errMsg )
   if ( errStat >= AbortErrLev ) return
   
     ! write a new line (advance to the next line)
   write (OutUnit,'()')
      
end subroutine Dvr_WriteOutputLine
