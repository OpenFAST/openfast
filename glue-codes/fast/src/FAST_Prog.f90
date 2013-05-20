!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework. 
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
!    ElastoDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with ElastoDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
PROGRAM FAST
! This program models 2- or 3-bladed turbines of a standard configuration.
!.................................................................................................

   USE NWTC_Library
   USE FAST_Types

USE FAST_IO_Subs       

   USE ElastoDyn
   USE ElastoDyn_Types

   USE ServoDyn
   USE ServoDyn_Types

   USE AeroDyn
   USE AeroDyn_Types

   USE HydroDyn
   USE HydroDyn_Types
   
   
IMPLICIT  NONE


   ! Local variables:

   ! Data for the glue code:
TYPE(FAST_ParameterType)       :: p_FAST                                     ! Parameters for the glue code (bjj: made global for now)
TYPE(FAST_OutputType)          :: y_FAST                                     ! Output variables for the glue code 

   ! Data for the ElastoDyn module:
TYPE(ED_InitInputType)         :: InitInData_ED                              ! Initialization input data
TYPE(ED_InitOutputType)        :: InitOutData_ED                             ! Initialization output data
TYPE(ED_ContinuousStateType)   :: x_ED                                       ! Continuous states
TYPE(ED_DiscreteStateType)     :: xd_ED                                      ! Discrete states
TYPE(ED_ConstraintStateType)   :: z_ED                                       ! Constraint states
TYPE(ED_OtherStateType)        :: OtherSt_ED                                 ! Other/optimization states
TYPE(ED_ParameterType)         :: p_ED                                       ! Parameters
TYPE(ED_InputType)             :: u_ED                                       ! System inputs
TYPE(ED_OutputType)            :: y_ED                                       ! System outputs

   ! Data for the ServoDyn module:
TYPE(SrvD_InitInputType)       :: InitInData_SrvD                            ! Initialization input data
TYPE(SrvD_InitOutputType)      :: InitOutData_SrvD                           ! Initialization output data
TYPE(SrvD_ContinuousStateType) :: x_SrvD                                     ! Continuous states
TYPE(SrvD_DiscreteStateType)   :: xd_SrvD                                    ! Discrete states
TYPE(SrvD_ConstraintStateType) :: z_SrvD                                     ! Constraint states
TYPE(SrvD_OtherStateType)      :: OtherSt_SrvD                               ! Other/optimization states
TYPE(SrvD_ParameterType)       :: p_SrvD                                     ! Parameters
TYPE(SrvD_InputType)           :: u_SrvD                                     ! System inputs
TYPE(SrvD_OutputType)          :: y_SrvD                                     ! System outputs

   ! Data for the AeroDyn module:

   ! Data for InflowWind module:
REAL(ReKi)                     :: IfW_WriteOutput(3)                         ! Temporary hack for getting wind speeds from InflowWind
   
   ! Data for the HydroDyn module:

   ! Other/Misc variables
REAL(DbKi)                     :: TiLstPrn                                   ! The time of the last print
REAL(DbKi)                     :: ZTime                                      ! Current simulation time
REAL(DbKi)                     :: OutTime                                    ! Used to determine if output should be generated at this simulation time
REAL(ReKi)                     :: PrevClockTime                              ! Clock time at start of simulation in seconds 
REAL                           :: UsrTime1                                   ! User CPU time for simulation initialization

INTEGER(IntKi)                 :: J                                          ! generic loop counter
INTEGER                        :: StrtTime (8)                               ! Start time of simulation
INTEGER(IntKi)                 :: Step                                       ! Current simulation time step.
INTEGER(IntKi)                 :: ErrStat                                    ! Error status
CHARACTER(1024)                :: ErrMsg                                     ! Error message
REAL(DbKi), PARAMETER          :: t_initial = 0.0                            ! Initial time
REAL(DbKi)                     :: dt_global                                  ! we're limiting our simulation to lock-step time steps for now

INTEGER, PARAMETER :: ED_interp_order = 0
TYPE(ED_InputType)   :: ED_Input(ED_interp_order+1)
REAL(DbKi)           :: ED_InputTimes(ED_interp_order+1)
  

   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

      ! Get the current time 
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
   Step          = 0                                                    ! The first step counter

   AbortErrLev   = ErrID_Fatal                                          ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
   
      ! Initialize NWTC Library (open console, set pi constants)  
   CALL NWTC_Init( ProgNameIN=FAST_ver%Name, EchoLibVer=.FALSE. )       ! sets the pi constants, open console for output, etc...
   

      ! Open and read input files, initialize global parameters.
   CALL FAST_Init( p_FAST, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from FAST_Init: '//NewLine//ErrMsg )

   dt_global = p_FAST%dt

   
   ! We fill ED_InputTimes with negative times, but the ED_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(ED_Input)
   do j = 1, ED_interp_order + 1  
      ED_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      !Mod1_OutputTimes(i) = t_initial - (j - 1) * dt
   enddo   
   
      ! initialize ElastoDyn (must be done first)
   InitInData_ED%InputFile     = p_FAST%EDFile
   InitInData_ED%ADInputFile   = p_FAST%ADFile
   InitInData_ED%RootName      = p_FAST%OutFileRoot
   CALL ED_Init( InitInData_ED, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, dt_global, InitOutData_ED, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from ED_Init: '//NewLine//ErrMsg )

   IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in ElastoDyn must be the same as the value of DT in FAST.")
   
   
   CALL ED_CopyInput (ED_Input(1),  u_ED, 0, Errstat, ErrMsg)   !Bjj copying input here only to allocate the arrays, until we get a better solution...
   
      ! initialize ServoDyn
   IF ( p_FAST%CompServo ) THEN
      InitInData_SrvD%InputFile     = p_FAST%SrvDFile
      InitInData_SrvD%RootName      = p_FAST%OutFileRoot
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      CALL AllocAry(InitInData_SrvD%BlPitchInit, InitOutData_ED%NumBl, 'BlPitchInit', ErrStat, ErrMsg)
      CALL CheckError( ErrStat, ErrMsg )
      
      InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, dt_global, InitOutData_SrvD, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, 'Message from SrvD_Init: '//NewLine//ErrMsg )
      
      !IF ( InitOutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!
      
      IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in ServoDyn must be the same as the value of DT in FAST.")
      
      !CALL SrvD_CopyInput (SrvD_Input(1),  u_ED, 0, Errstat, ErrMsg)   !Bjj copying input here only to allocate the arrays, until we get a better solution...
      
      ! initialize y%ElecPwr and y%GenTq because they are one timestep different (used as input for the next step)
      y_SrvD%ElecPwr = 0.0
      y_SrvD%GenTrq   = 0.0
      
   END IF

   
      ! initialize AeroDyn
   IF ( p_FAST%CompAero ) THEN
   ! we need the air density (and wind speed) yet.... some strangeness still going on.
      CALL AeroInput(p_ED, p_FAST)            ! Read in the ADFile
   
         ! some weirdness that we probably won't need anymore....
      p_ED%AirDens   = AD_GetConstant('AirDensity', ErrStat)
      
   ELSE
      p_ED%AirDens = 0
      IfW_WriteOutput = 0.0
   END IF



   ! initialize HydroDyn
   IF ( p_FAST%CompHydro ) THEN
      
      HydroDyn_InitData%Gravity     = InitOutData_ED%Gravity      ! m/s^2   
      HydroDyn_InitData%FileName    = p_FAST%HDFile    
      HydroDyn_InitData%OutRootName = p_FAST%OutFileRoot
      HD_ConfigMarkers%Substructure%Position = (/0._ReKi, 0._ReKi, 0._ReKi/)   

      CALL HD_Init(HydroDyn_InitData, HD_ConfigMarkers, HD_AllMarkers, HydroDyn_data, ErrStat)
      CALL CheckError( ErrStat, 'Error initializing HydroDyn.' )


      ! get the mapping of the Hydro markers to the structural markers.  For now, we require that they
      ! be the same (otherwise, we must change RtHS() to interpolate the hydro loads to the correct structural markers.

         ! >>> BJJ this is a hack job for now.  We need to actually determine if the loads are tower (per-unit-length) or platform (lumped-sum)
      HD_TwrNodes = .FALSE.
   
      IF ( SIZE( HD_AllMarkers%Substructure ) ==  1        ) THEN !.AND. LUMPED-SUM LOAD  (otherwise, we can't have tower loads with only 1 element)

         DO J=1,3
            IF ( .NOT. EqualRealNos( HD_AllMarkers%Substructure(1)%Position(J), HD_ConfigMarkers%Substructure%Position(J) ) ) THEN
               CALL CheckError( ErrID_Fatal, ' ElastoDyn and HydroDyn must have the same substructure node.' )
            END IF
         END DO

      ELSEIF ( SIZE( HD_AllMarkers%Substructure ) ==  p_ED%TwrNodes )  THEN ! .AND. PER-UNIT-LENGTH LOADS

            ! We currently require that the tower nodes in FAST be the same as in HydroDyn
         DO J=1,p_ED%TwrNodes      
            IF ( .NOT. EqualRealNos( p_ED%HNodes(J) + p_ED%TwrRBHt - p_ED%TwrDraft, HD_AllMarkers%Substructure(J)%Position(3) ) ) THEN
               CALL CheckError( ErrID_Fatal, ' ElastoDyn and HydroDyn must have the same tower nodes.' )
            ELSEIF ( .NOT. EqualRealNos( 0.0_ReKi, HD_AllMarkers%Substructure(J)%Position(1) ) ) THEN
               CALL CheckError( ErrID_Fatal, ' HydroDyn tower markers must have X = 0 m.' )
            ELSEIF ( .NOT. EqualRealNos( 0.0_ReKi, HD_AllMarkers%Substructure(J)%Position(2) ) ) THEN
               CALL CheckError( ErrID_Fatal, ' HydroDyn tower markers must have Y = 0 m.' )
            ENDIF
         END DO !J: TwrNodes

         HD_TwrNodes = .TRUE.

      ELSE
         CALL CheckError( ErrID_Fatal, " Unable to discern HydroDyn's discretization." )
      ENDIF

      ! <<<<
      
   END IF   ! CompHydro
      

   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_SrvD, AD_Prog, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from FAST_InitOutput: '//NewLine//ErrMsg )
    
   
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................

   CALL ED_DestroyInitInput(  InitInData_ED, ErrStat, ErrMsg )
   CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat, ErrMsg )
   
   CALL SrvD_DestroyInitInput(  InitInData_SrvD, ErrStat, ErrMsg )
   CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat, ErrMsg )

   !...............................................................................................................................
   ! loose coupling
   !...............................................................................................................................

      ! Start simulation.  Initialize the simulation status.

   CALL WrScr1 ( '' )
!   CALL SimStatus ()  


!.................................................................
!BJJ: NOTE: there is currently a time shift in this algorithm that we need to fix, 
!  but I will wait until we get AeroDyn and InflowWind merged in this mix, then
!  use the glue code developed by M. Sprague in the Gasmi Paper Examples.
!.................................................................

      ! Loop through time.

   Step  = 0_IntKi
   ZTime = 0.0_DbKi
   DO 
      
            ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.
      
      CALL ED_Input_ExtrapInterp(ED_Input, ED_InputTimes, u_ED, ZTime + p_FAST%dt, ErrStat, ErrMsg)
      
         ! Shift "window" of the Mod1_Input and Mod1_Output
      
      do j = ED_interp_order, 1, -1
         Call ED_CopyInput (ED_Input(j),  ED_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(j+1),  0, Errstat, ErrMsg)
         ED_InputTimes(j+1) = ED_InputTimes(j)
         !Mod1_OutputTimes(j+1) = Mod1_OutputTimes(j)
      enddo
      
      CALL ED_CopyInput (u_ED,  ED_Input(1),  0, Errstat, ErrMsg)
      !Call Mod1_CopyOutput (y1,  Mod1_Output(1),  0, Errstat, ErrMsg)
      ED_InputTimes(1) = ZTime + p_FAST%dt

      !.....................................................
      ! Call predictor-corrector routine:
      !.....................................................

         ! ElastoDyn
      CALL ED_UpdateStates( ZTime, Step, ED_Input, ED_InputTimes, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, 'Message from ED_UpdateStates: '//NewLine//ErrMsg )

         ! ServoDyn
         
         ! AeroDyn
      
         ! HydroDyn
      IF ( p_FAST%CompHydro ) THEN
      END IF         

      
      
      
      !.....................................................
      ! Advance time:
      !.....................................................

      Step  = Step + 1
      ZTime = Step*p_FAST%DT


      !.....................................................
      ! Input-Output solve:
      !.....................................................      
      
      CALL ED_CalcOutput( ZTime, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from ED_CalcOutput: '//NewLine//ErrMsg  )

      IF ( p_FAST%CompAero ) THEN
         CALL AD_InputSolve( p_ED, x_ED, OtherSt_ED, u_ED, y_ED, ErrStat, ErrMsg )
         ADAeroLoads = AD_CalculateLoads( REAL(ZTime, ReKi), ADAeroMarkers, ADInterfaceComponents, ADIntrfaceOptions, ErrStat )
            CALL CheckError( ErrStat, ' Error calculating hydrodynamic loads in AeroDyn.'  )
         
            !InflowWind outputs
         IfW_WriteOutput = AD_GetUndisturbedWind( REAL(ZTime, ReKi), (/0.0_ReKi, 0.0_ReKi, p_ED%FASTHH /), ErrStat )            
            CALL CheckError( ErrStat, 'Message from IfW_CalcOutput: '//NewLine//ErrMsg  )
            
      END IF

      IF ( p_FAST%CompServo ) THEN
         CALL SrvD_InputSolve( p_FAST, u_SrvD, y_ED, IfW_WriteOutput, y_SrvD   )  !use the ServoDyn outputs from Step = Step-1 (bjj: need to think about this for predictor-corrector (make sure it doesn't get changed...))

         CALL SrvD_CalcOutput( ZTime, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from SrvD_CalcOutput: '//NewLine//ErrMsg  )
      END IF
            
      IF ( p_FAST%CompHydro ) THEN
         CALL HD_InputSolve( p_ED, x_ED, OtherSt_ED, u_ED, y_ED, ErrStat, ErrMsg )
         CALL HD_CalculateLoads( ZTime,  HD_AllMarkers,  HydroDyn_data, HD_AllLoads,  ErrStat )
            CALL CheckError( ErrStat, 'Error calculating hydrodynamic loads in HydroDyn.'  )
      END IF
      
         ! User Tower Loading
      IF ( p_FAST%CompUserTwrLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrTwr_CalcOutput()
      !   CALL UserTwrLd ( JNode, X, XD, t, p_FAST%DirRoot, y_UsrTwr%AddedMass(1:6,1:6,J), (/ y_UsrTwr%Force(:,J),y_UsrTwr%Moment(:,J) /) )
      END IF
      
         ! User Platform Loading
      IF ( p_FAST%CompUserPtfmLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrPtfm_CalcOutput()
      !   
      !   CALL UserPtfmLd ( x_ED%QT(1:6), x_%QDT(1:6), t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, (/ y_UsrPtfm%Force,y_UsrPtfm%Moment /) )
      !
      !      ! Ensure that the platform added mass matrix returned by UserPtfmLd, PtfmAM, is symmetric; Abort if necessary:
      !   IF ( .NOT. IsSymmetric( y_UsrPtfm%AddedMass ) ) THEN
      !      CALL CheckError ( ErrID_Fatal, ' The user-defined platform added mass matrix is unsymmetric.'// &
      !                        '  Make sure AddedMass returned by UserPtfmLd() is symmetric.'        )
      !   END IF
      !   
      END IF      
      
      
      CALL ED_InputSolve( p_FAST, p_ED,  ED_Input(1),   y_SrvD )
                   
                  
      !......................................................
      ! Check to see if we should output data this time step:
      !......................................................

      IF ( ZTime >= p_FAST%TStart )  THEN
         
            !bjj FIX THIS algorithm!!! this assumes dt_out is an integer multiple of dt; we will probably have to do some interpolation to get these outputs at the times we want them....
         OutTime = NINT( ZTime / p_FAST%DT_out ) * p_FAST%DT_out      
         IF ( EqualRealNos( ZTime, OutTime ) )  THEN 
                           
               ! Generate glue-code output file
            CALL WrOutputLine( ZTime, p_FAST, y_FAST, EDOutput=y_ED%WriteOutput, SrvDOutput=y_SrvD%WriteOutput, IfWOutput=IfW_WriteOutput, ErrStat=ErrStat, ErrMsg=ErrMsg )
            CALL CheckError( ErrStat, ErrMsg )
            
               ! Generate AeroDyn's element data if desired:
            CALL ElemOut()
            
         END IF
         
      ENDIF

      !.....................................................
      ! Display simulation status every SttsTime-seconds:
      !.....................................................

      IF ( ZTime - TiLstPrn >= p_FAST%SttsTime )  THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, ZTime, p_FAST%TMax )

      ENDIF


      ! If we've reached TMax, exit the DO loop:

      IF ( ZTime > p_FAST%TMax )  EXIT

   ENDDO        

   !...............................................................................................................................
   !  Write simulation times and stop
   !...............................................................................................................................

   CALL ExitThisProgram( Error=.FALSE. )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE ExitThisProgram( Error, ErrLev )
   ! This subroutine cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      LOGICAL,        INTENT(IN)           :: Error        ! flag to determine if this is an abort or normal stop
      INTEGER(IntKi), INTENT(IN), OPTIONAL :: ErrLev       ! Error level when Error == .TRUE. (required when Error is .TRUE.)
      
      !...............................................................................................................................
      ! Clean up modules (and write binary FAST output file), destroy any other variables
      !...............................................................................................................................
   
      CALL FAST_End( p_FAST, y_FAST, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL ED_End(   u_ED,   p_ED,   x_ED,   xd_ED,   z_ED,   OtherSt_ED,   y_ED,   ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL SrvD_End( u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL AeroDyn_End( ErrStat )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( 'Error ending AeroDyn' )

      CALL HD_Terminate( HydroDyn_data, ErrStat )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      
      ! in case we didn't get these destroyed earlier....
      
      CALL ED_DestroyInitInput(  InitInData_ED, ErrStat, ErrMsg )
      CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat, ErrMsg )
   
      CALL SrvD_DestroyInitInput(  InitInData_SrvD, ErrStat, ErrMsg )
      CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat, ErrMsg )
      
      
      !............................................................................................................................
      ! Set exit error code if there was an error;
      !............................................................................................................................

      IF (Error) CALL ProgAbort( ' Abort error level: '//TRIM(GetErrStr(ErrLev) ) )  !This assumes PRESENT(ErrID) is .TRUE.
      
      !............................................................................................................................
      !  Write simulation times and stop
      !............................................................................................................................

      CALL RunTimes( StrtTime, UsrTime1, ZTime )     
      
      CALL NormStop( )
      

   END SUBROUTINE ExitThisProgram    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)
      

      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         IF ( ErrID >= AbortErrLev ) CALL ExitThisProgram( Error=.TRUE., ErrLev=ErrID )
      END IF
            

   END SUBROUTINE CheckError       
END PROGRAM FAST
!=======================================================================
