!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
!    This file is part of the NWTC Subroutine Library.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
PROGRAM BeamDyn_Driver_Program

   USE BeamDyn_driver_subs  ! all other modules inherited through this one
   USE VersionInfo

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                   :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                  :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   REAL(DbKi)                       :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                       :: t_global         ! global-loop time marker
   INTEGER(IntKi)                   :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                   :: n_t_global       ! global-loop time counter
   INTEGER(IntKi)                   :: n_t_vtk          ! global vtk step counter
   INTEGER(IntKi), parameter        :: BD_interp_order = 1  ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InitInputType)           :: BD_InitInput
   TYPE(BD_ParameterType)           :: BD_Parameter
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   TYPE(BD_InitOutputType)          :: BD_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   TYPE(BD_OtherStateType)          :: BD_OtherState
   TYPE(BD_MiscVarType)             :: BD_MiscVar
   TYPE(BD_InputType) ,ALLOCATABLE  :: BD_Input(:)
   REAL(DbKi),         ALLOCATABLE  :: BD_InputTimes(:)
   TYPE(BD_OutputType)              :: BD_Output
   INTEGER(IntKi)                   :: DvrOut 
   
   TYPE(BD_DriverInternalType)      :: DvrData

   ! local variables
   
   CHARACTER(256)                   :: DvrInputFile
   CHARACTER(256)                   :: RootName
   INTEGER(IntKi)                   :: j             ! counter for various loops
   INTEGER(IntKi)                   :: i             ! counter for various loops   
   INTEGER(IntKi)                   :: max_ld_step=8 ! maximum load steps for static runs.
   REAL(DbKi)                       :: TiLstPrn      ! The simulation time of the last print (to file) [(s)]
   REAL(ReKi)                       :: PrevClockTime ! Clock time at start of simulation in seconds [(s)]
   REAL(ReKi)                       :: UsrTime1      ! User CPU time for simulation initialization [(s)]
   REAL(ReKi)                       :: UsrTime2      ! User CPU time for simulation (without intialization) [(s)]
   INTEGER(IntKi) , DIMENSION(1:8)  :: StrtTime      ! Start time of simulation (including intialization) [-]
   INTEGER(IntKi) , DIMENSION(1:8)  :: SimStrtTime   ! Start time of simulation (after initialization) [-]
   CHARACTER(200)                   :: git_commit    ! String containing the current git commit hash

   TYPE(ProgDesc), PARAMETER        :: version   = ProgDesc( 'BeamDyn Driver', '', '' )  ! The version number of this program.
   

   ! -------------------------------------------------------------------------
   ! Initialization of library (especially for screen output)
   ! -------------------------------------------------------------------------  
   
   CALL DATE_AND_TIME ( Values=StrtTime )                 ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                             ! Initial time (this zeros the start time when used as a MATLAB function)
   UsrTime1 = MAX( 0.0_ReKi, UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned

   
   CALL NWTC_Init()
      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
   
   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------   
   
   CALL GET_COMMAND_ARGUMENT(1,DvrInputFile)
   CALL GetRoot(DvrInputFile,RootName)
   CALL BD_ReadDvrFile(DvrInputFile,dt_global,BD_InitInput,DvrData,ErrStat,ErrMsg)
      CALL CheckError()
      
      ! initialize the BD_InitInput values not in the driver input file
   BD_InitInput%RootName = TRIM(BD_Initinput%InputFile)
   BD_InitInput%RootName = TRIM(RootName)//'.BD'
   BD_InitInput%RootDisp = MATMUL(BD_InitInput%GlbPos(:),DvrData%RootRelInit) - BD_InitInput%GlbPos(:)
   BD_InitInput%DynamicSolve = DvrData%DynamicSolve      ! QuasiStatic options handled within the BD code.

   t_global = DvrData%t_initial
   n_t_final = ((DvrData%t_final - DvrData%t_initial) / dt_global )

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1)) 
   ALLOCATE(BD_InputTimes(BD_interp_order + 1)) 

   CALL BD_Init(BD_InitInput             &
                   , BD_Input(1)         &
                   , BD_Parameter        &
                   , BD_ContinuousState  &
                   , BD_DiscreteState    &
                   , BD_ConstraintState  &
                   , BD_OtherState       &
                   , BD_Output           &
                   , BD_MiscVar          &
                   , dt_global           &
                   , BD_InitOutput       &
                   , ErrStat             &
                   , ErrMsg )
      CALL CheckError()
   
      ! If the Quasi-Static solve is in use, rerun the initialization with loads at t=0 
      ! (HACK: set in the driver only because computing Jacobians with this option [as in FAST glue code] is problematic)
   BD_OtherState%RunQuasiStaticInit = BD_Parameter%analysis_type == BD_DYN_SSS_ANALYSIS


     ! Set the Initial root orientation
   BD_Input(1)%RootMotion%Orientation(1:3,1:3,1) = DvrData%RootRelInit
      
   call Init_RotationCenterMesh(DvrData, BD_InitInput, BD_Input(1)%RootMotion, ErrStat, ErrMsg)
      CALL CheckError()

   call CreateMultiPointMeshes(DvrData,BD_InitInput,BD_InitOutput,BD_Parameter, BD_Output, BD_Input(1), ErrStat, ErrMsg)   
   call Transfer_MultipointLoads(DvrData, BD_Output, BD_Input(1), ErrStat, ErrMsg)   
   
   CALL Dvr_InitializeOutputFile(DvrOut,BD_InitOutput,RootName,ErrStat,ErrMsg)
      CALL CheckError()
      
      
      ! initialize BD_Input and BD_InputTimes
   BD_InputTimes(1) = DvrData%t_initial
   CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), DvrData, ErrStat, ErrMsg)
   
   DO j = 2,BD_interp_order+1
         ! create new meshes
      CALL BD_CopyInput (BD_Input(1) , BD_Input(j) , MESH_NEWCOPY, ErrStat, ErrMsg)
         CALL CheckError()
         
         ! solve for inputs at previous time steps
      BD_InputTimes(j) = DvrData%t_initial - (j - 1) * dt_global
      CALL BD_InputSolve( BD_InputTimes(j), BD_Input(j), DvrData, ErrStat, ErrMsg)
         CALL CheckError()
   END DO
   
      ! Write VTK reference if requested (ref is (0,0,0)
   if (DvrData%WrVTK > 0) then
      call SetVTKvars()
      call MeshWrVTKreference( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Output%BldMotion,   trim(DvrData%VTK_OutFileRoot)//'_BldMotion', ErrStat, ErrMsg );  call CheckError()
      call MeshWrVTKreference( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Input(1)%PointLoad, trim(DvrData%VTK_OutFileRoot)//'_PointLoad', ErrStat, ErrMsg );  call CheckError()
      call MeshWrVTKreference( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Input(1)%DistrLoad, trim(DvrData%VTK_OutFileRoot)//'_DistrLoad', ErrStat, ErrMsg );  call CheckError()
   endif
      ! Write VTK reference if requested (ref is (0,0,0)
   if (DvrData%WrVTK == 2) then
      n_t_vtk = 0
      call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Output%BldMotion,   trim(DvrData%VTK_OutFileRoot)//'_BldMotion',  n_t_vtk, .true., ErrStat, ErrMsg, DvrData%VTK_tWidth )
      call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Input(1)%PointLoad, trim(DvrData%VTK_OutFileRoot)//'_PointLoad',  n_t_vtk, .true., ErrStat, ErrMsg, DvrData%VTK_tWidth )
      call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Input(1)%DistrLoad, trim(DvrData%VTK_OutFileRoot)//'_DistrLoad',  n_t_vtk, .true., ErrStat, ErrMsg, DvrData%VTK_tWidth )
      call CheckError()
   endif


      !.........................
      ! calculate outputs at t=0
      !.........................
   CALL SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, t_global, DvrData%t_final )

   CALL BD_CalcOutput( t_global, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                           BD_ConstraintState, BD_OtherState,  BD_Output, BD_MiscVar, ErrStat, ErrMsg)
      CALL CheckError()
   
   CALL Dvr_WriteOutputLine(t_global,DvrOut,BD_Parameter%OutFmt,BD_Output)
   
      !.........................
      ! time marching
      !.........................
     
   DO n_t_global = 0, n_t_final

      ! Shift "window" of BD_Input 
      DO j = BD_interp_order, 1, -1
         CALL BD_CopyInput (BD_Input(j),  BD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL CheckError()
         BD_InputTimes(j+1) = BD_InputTimes(j)
      END DO
      
      BD_InputTimes(1)  = t_global + dt_global
      CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), DvrData, ErrStat, ErrMsg)
         CALL CheckError()
      
                       
     IF(BD_Parameter%analysis_type .EQ. BD_STATIC_ANALYSIS .AND. n_t_global > max_ld_step) EXIT

      ! update states from n_t_global to n_t_global + 1
     CALL BD_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, &
                               BD_ContinuousState, &
                               BD_DiscreteState, BD_ConstraintState, &
                               BD_OtherState, BD_MiscVar, ErrStat, ErrMsg )
        CALL CheckError()

        
      ! advance time
     t_global = (n_t_global+1) * dt_global + DvrData%t_initial
           
      ! calculate outputs at n_t_global + 1
     CALL BD_CalcOutput( t_global, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                             BD_ConstraintState, BD_OtherState,  BD_Output, BD_MiscVar, ErrStat, ErrMsg)
        CALL CheckError()

     CALL Dvr_WriteOutputLine(t_global,DvrOut,BD_Parameter%OutFmt,BD_Output)

         ! Write VTK reference if requested (ref is (0,0,0)
      if (DvrData%WrVTK == 2) then
         if ( MOD( n_t_global, DvrData%n_VTKTime ) == 0 ) then
            n_t_vtk = n_t_vtk + 1
            call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Output%BldMotion,   trim(DvrData%VTK_OutFileRoot)//'_BldMotion', n_t_vtk, .true., ErrStat, ErrMsg, DvrData%VTK_tWidth )
            call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Input(1)%PointLoad, trim(DvrData%VTK_OutFileRoot)//'_PointLoad', n_t_vtk, .true., ErrStat, ErrMsg, DvrData%VTK_tWidth )
            call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), BD_Input(1)%DistrLoad, trim(DvrData%VTK_OutFileRoot)//'_DistrLoad', n_t_vtk, .true., ErrStat, ErrMsg, DvrData%VTK_tWidth )
            call CheckError()
         endif
      endif


     if ( MOD( n_t_global + 1, 100 ) == 0 ) call SimStatus( TiLstPrn, PrevClockTime, t_global, DvrData%t_final )
   ENDDO
      
   CALL RunTimes( StrtTime, UsrTime1, SimStrtTime, UsrTime2, t_global )

   
   CALL Dvr_End()

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Dvr_End()

      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /=
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      character(*), parameter                       :: RoutineName = 'Dvr_End'

      IF(DvrOut >0) CLOSE(DvrOut)

      IF ( ALLOCATED(BD_Input) ) THEN
         CALL BD_End( BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
               BD_ConstraintState, BD_OtherState, BD_Output, BD_MiscVar, ErrStat2, ErrMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
         DO i=2,BD_interp_order + 1
            CALL BD_DestroyInput( BD_Input(i), ErrStat2, ErrMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         ENDDO
         
         DEALLOCATE(BD_Input)
      END IF

      IF(ALLOCATED(BD_InputTimes )) DEALLOCATE(BD_InputTimes )
      if(allocated(DvrData%MultiPointLoad)) deallocate(DvrData%MultiPointLoad)

      
      
      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'BeamDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
   END SUBROUTINE Dvr_End
!----------------------------------------------------------------------------------------------------------------------------------
   subroutine CheckError()
   
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            call Dvr_End()
         end if
      end if
         
   end subroutine CheckError

   subroutine SetVTKvars()
      real(R8Ki)  :: TmpTime
      real(R8Ki)  :: TmpRate
      real(R8Ki)  :: TotalTime

      DvrData%VTK_OutFileRoot = trim(BD_InitInput%RootName)
      n_t_vtk = 0    ! first VTK output number

      ! convert frames-per-second to seconds per sample:
      TotalTime = DvrData%t_final - DvrData%t_initial
      if ( DvrData%VTK_fps == 0 ) then
         TmpTime = TotalTime + dt_global 
      else
         TmpTime = 1.0_R8Ki / DvrData%VTK_fps
      endif

      ! now save the number of time steps between VTK file output:
      if (TmpTime > TotalTime) then
         DvrData%n_VTKTime = HUGE(DvrData%n_VTKTime)
      else
         DvrData%n_VTKTime = NINT( TmpTime / dt_global )
         ! I'll warn if p%n_VTKTime*p%DT is not TmpTime
         IF (DvrData%WrVTK == 2) THEN
            TmpRate = DvrData%n_VTKTime*dt_global
            if (.not. EqualRealNos(TmpRate, TmpTime)) then
               call WrScr('1/VTK_fps is not an integer multiple of DT. FAST will output VTK information at '//&
                              trim(num2lstr(1.0_DbKi/TmpRate))//' fps, the closest rate possible.')
            end if
         end if
      end if

      DvrData%VTK_tWidth = CEILING( log10( real(n_t_final, ReKi) / DvrData%n_VTKTime ) ) + 1
   end subroutine SetVTKvars

END PROGRAM BeamDyn_Driver_Program
