!**********************************************************************************************************************************
!! This module is used to read and process the (undisturbed) inflow winds.  It must be initialized
!! using InflowWind_Init() with the name of the file, the file type, and possibly reference height and
!! width (depending on the type of wind file being used).  This module calls appropriate routines
!! in the wind modules so that the type of wind becomes seamless to the user.  InflowWind_End()
!! should be called when the program has finshed.
!!
!! Data are assumed to be in units of meters and seconds.  Z is measured from the ground (NOT the hub!).
!!
!!  7 Oct 2009    Initial Release with AeroDyn 13.00.00         B. Jonkman, NREL/NWTC
!! 14 Nov 2011    v1.00.01b-bjj                                 B. Jonkman
!!  1 Aug 2012    v1.01.00a-bjj                                 B. Jonkman
!! 10 Aug 2012    v1.01.00b-bjj                                 B. Jonkman
!!    Feb 2013    v2.00.00a-adp   conversion to Framework       A. Platt
!!    Sep 2015    v3.00.00a-adb   added separate input file     A. Platt
!
!..................................................................................................................................
! Files with this module:
!  InflowWind_Subs.f90
!  InflowWind.txt       -- InflowWind_Types will be auto-generated based on the descriptions found in this file.
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2017  Envision Energy
!
!    This file is part of InflowWind.
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
MODULE InflowWind


   USE                              NWTC_Library
   USE                              InflowWind_Types
   USE                              InflowWind_Subs
   USE                              InflowWind_IO_Types
   USE                              InflowWind_IO
   USE                              IfW_FlowField

   USE                              Lidar                      ! module for obtaining sensor data
   
   IMPLICIT NONE
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: IfW_Ver = ProgDesc( 'InflowWind', '', '' )
   integer,        parameter            :: NumExtendedInputs = 3 !: V, VShr, PropDir




      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: InflowWind_Init                                   !< Initialization routine
   PUBLIC :: InflowWind_CalcOutput                             !< Calculate the wind velocities
   PUBLIC :: InflowWind_End                                    !< Ending routine (includes clean up)

      ! These routines satisfy the framework, but do nothing at present.
   PUBLIC :: InflowWind_UpdateStates               !< Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   PUBLIC :: InflowWind_CalcConstrStateResidual    !< Tight coupling routine for returning the constraint state residual
   PUBLIC :: InflowWind_CalcContStateDeriv         !< Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: InflowWind_UpdateDiscState            !< Tight coupling routine for updating discrete states


      ! These routines compute Jacobians; only dYdu is defined.
   PUBLIC :: InflowWind_JacobianPInput             !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
   PUBLIC :: InflowWind_JacobianPContState         !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                                   !!   states(x)
   PUBLIC :: InflowWind_JacobianPDiscState         !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                                   !!   states(xd)
   PUBLIC :: InflowWind_JacobianPConstrState       !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                                   !!   states(z)
   PUBLIC :: InflowWind_GetOP                      !< Routine to pack the operating point values (for linearization) into arrays
   


CONTAINS
!====================================================================================================
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! Since this module acts as an interface to other modules, on some things are set before initiating
!! calls to the lower modules.
!----------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                           y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )

   ! Initialization data and guesses
   TYPE(InflowWind_InitInputType),        INTENT(IN   )  :: InitInp           !< Input data for initialization
   TYPE(InflowWind_InputType),            INTENT(  OUT)  :: InputGuess        !< An initial guess for the input; the input mesh must be defined
   TYPE(InflowWind_ParameterType),        INTENT(  OUT)  :: p                 !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(  OUT)  :: ContStates        !< Initial continuous states
   TYPE(InflowWind_DiscreteStateType),    INTENT(  OUT)  :: DiscStates        !< Initial discrete states
   TYPE(InflowWind_ConstraintStateType),  INTENT(  OUT)  :: ConstrStateGuess  !< Initial guess of the constraint states
   TYPE(InflowWind_OtherStateType),       INTENT(  OUT)  :: OtherStates       !< Initial other/optimization states
   TYPE(InflowWind_OutputType),           INTENT(  OUT)  :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   TYPE(InflowWind_MiscVarType),          INTENT(  OUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   REAL(DbKi),                            INTENT(IN   )  :: TimeInterval      !< Coupling time interval in seconds: InflowWind does not change this.
   TYPE(InflowWind_InitOutputType),       INTENT(  OUT)  :: InitOutData       !< Initial output data -- Names, units, and version info.
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   ! Local variables
   CHARACTER(*), PARAMETER                               :: RoutineName="InflowWind_Init"

   TYPE(InflowWind_InputFile)                            :: InputFileData        !< Data from input file

   Type(Steady_InitInputType)                            :: Steady_InitInput
   Type(Uniform_InitInputType)                           :: Uniform_InitInput
   Type(TurbSim_InitInputType)                           :: TurbSim_InitInput
   Type(Bladed_InitInputType)                            :: Bladed_InitInput
   Type(Bladed_InitOutputType)                           :: Bladed_InitOutput
   Type(HAWC_InitInputType)                              :: HAWC_InitInput
   Type(User_InitInputType)                              :: User_InitInput
   Type(Grid4D_InitInputType)                            :: Grid4D_InitInput
   Type(Points_InitInputType)                            :: Points_InitInput

   Type(WindFileDat)                                     :: FileDat
   

   ! TYPE(InflowWind_IO_InitInputType)                      :: FlowField_InitData     !< initialization info
   ! TYPE(InflowWind_IO_InitOutputType)                     :: FlowField_InitOutData  !< initialization output info

   TYPE(FileInfoType)                                    :: InFileInfo    !< The derived type for holding the full input file for parsing -- we may pass this in the future
   CHARACTER(1024)                                       :: PriPath

   INTEGER(IntKi)                                        :: I, j              !< Generic counter
   INTEGER(IntKi)                                        :: Lin_indx          !< Generic counter
   INTEGER(IntKi)                                        :: SumFileUnit       !< Unit number for the summary file
   CHARACTER(256)                                        :: SumFileName       !< Name of the summary file
   CHARACTER(256)                                        :: EchoFileName      !< Name of the summary file
   CHARACTER(1), PARAMETER                               :: UVW(3) = (/'U','V','W'/)
   CHARACTER(1), PARAMETER                               :: XYZ(3) = (/'X','Y','Z'/)
   INTEGER(IntKi)                                        :: TmpErrStat
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg         !< temporary error message

   !----------------------------------------------------------------------------------------------
   ! Initialize variables and check to see if this module has been initialized before.
   !----------------------------------------------------------------------------------------------

   ErrStat        =  ErrID_None
   ErrMsg         =  ""
   SumFileUnit    =  -1_IntKi ! set at beginning in case of error

   ! Set a few variables.

   p%DT   = TimeInterval             ! InflowWind does not require a specific time interval, so this is never changed.
   CALL NWTC_Init()
   CALL DispNVD( IfW_Ver )

   !----------------------------------------------------------------------------------------------
   ! Read the input file
   !----------------------------------------------------------------------------------------------

   ! Set the names of the files based on the inputfilename
   p%RootFileName  = InitInp%RootName
   IF (LEN_TRIM(p%RootFileName) == 0) CALL GetRoot( InitInp%InputFileName, p%RootFileName )
   EchoFileName  = TRIM(p%RootFileName)//".ech"
   SumFileName   = TRIM(p%RootFileName)//".sum"







   ! Parse all the InflowWind related input files and populate the *_InitDataType derived types
   CALL GetPath( InitInp%InputFileName, PriPath )

   IF ( InitInp%UseInputFile ) THEN
      CALL ProcessComFile( InitInp%InputFileName, InFileInfo, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF        
   ELSE
      CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, InFileInfo, MESH_NEWCOPY, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF          
   ENDIF

   CALL InflowWind_ParseInputFileInfo( InputFileData,  InFileInfo, PriPath, InitInp%InputFileName, EchoFileName, &
                                       InitInp%FixedWindFileRootName, InitInp%TurbineID, TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   ENDIF

   ! If wind is Grid4D from FAST.Farm, set input file values
   IF (InitInp%Use4Dext) then
      InputFileData%WindType = FDext_WindNumber      
      InputFileData%PropagationDir = 0.0_ReKi ! wind is in XYZ coordinates (already rotated if necessary), so don't rotate it again
      InputFileData%VFlowAngle = 0.0_ReKi
      InputFileData%VelInterpCubic = .false.
   END IF

   ! initialize sensor data:   
   p%lidar%NumBeam = InputFileData%NumBeam
   p%lidar%RotorApexOffsetPos = InputFileData%RotorApexOffsetPos
   p%lidar%SensorType = InputFileData%SensorType      
   p%lidar%LidRadialVel   = InputFileData%LidRadialVel
   p%lidar%NumPulseGate = InputFileData%NumPulseGate
   p%lidar%FocalDistanceX =  InputFileData%FocalDistanceX
   p%lidar%FocalDistanceY =  InputFileData%FocalDistanceY
   p%lidar%FocalDistanceZ =  InputFileData%FocalDistanceZ
   p%lidar%MeasurementInterval = InputFileData%MeasurementInterval
   p%lidar%PulseSpacing = InputFileData%PulseSpacing
   p%lidar%URefLid = InputFileData%URefLid
   p%lidar%ConsiderHubMotion = InputFileData%ConsiderHubMotion  
         
         
   CALL Lidar_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates,   &
                    y, m, TimeInterval, InitOutData, TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )      
   

      ! Validate the InflowWind input file information.

   CALL InflowWind_ValidateInput( InitInp, InputFileData, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      
   IF ( ErrStat>= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   ENDIF


   
      ! If a summary file was requested, open it.
   IF ( InputFileData%SumPrint ) THEN

         ! Open the summary file and write some preliminary info to it
      CALL InflowWind_OpenSumFile( SumFileUnit, SumFileName, IfW_Ver, InputFileData%WindType, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         ENDIF
   ELSE
      SumFileUnit =  -1_IntKi       ! So that we don't try to write to something.  Used as indicator in submodules.
   ENDIF

   ! Allocate the array for passing points
   CALL AllocAry( InputGuess%PositionXYZ, 3, InitInp%NumWindPoints, "Array of positions at which to find wind velocities", TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   InputGuess%PositionXYZ = 0.0_ReKi
   InputGuess%HubPosition = 0.0_ReKi
   CALL Eye(InputGuess%HubOrientation,TmpErrStat,TmpErrMsg)  

   ! Allocate the array for passing velocities out
   CALL AllocAry( y%VelocityUVW, 3, InitInp%NumWindPoints, "Array of wind velocities returned by InflowWind", TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF (ErrStat>= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   ENDIF
   y%VelocityUVW = 0.0_ReKi

   ! If requested, allocate the array for passing accelerations out
   IF ( InitInp%OutputAccel ) THEN
      CALL AllocAry( y%AccelUVW, 3, InitInp%NumWindPoints, "Array of wind accelerations returned by InflowWind", TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat>= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF
      y%AccelUVW = 0.0_ReKi
   ENDIF

   !----------------------------------------------------------------------------
   ! Set flow field input data based on wind type
   !----------------------------------------------------------------------------

   InitOutData%WindFileInfo%MWS = HUGE(InitOutData%WindFileInfo%MWS)

   select case(InputFileData%WindType)

   case (Steady_WindNumber)

      Steady_InitInput%HWindSpeed = InputFileData%Steady_HWindSpeed
      Steady_InitInput%RefHt = InputFileData%Steady_RefHt
      Steady_InitInput%PLExp = InputFileData%Steady_PLexp

      p%FlowField%FieldType = Uniform_FieldType
      call IfW_SteadyWind_Init(Steady_InitInput, SumFileUnit, p%FlowField%Uniform, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Uniform%RefHeight]

   case (Uniform_WindNumber)

      Uniform_InitInput%WindFileName = InputFileData%Uniform_FileName
      Uniform_InitInput%RefHt = InputFileData%Uniform_RefHt
      Uniform_InitInput%RefLength = InputFileData%Uniform_RefLength
      Uniform_InitInput%PropagationDir = InputFileData%PropagationDir
      Uniform_InitInput%UseInputFile = InitInp%WindType2UseInputFile
      Uniform_InitInput%PassedFileData = InitInp%WindType2Data

      p%FlowField%FieldType = Uniform_FieldType
      call IfW_UniformWind_Init(Uniform_InitInput, SumFileUnit, p%FlowField%Uniform, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Uniform%RefHeight]

   case (TSFF_WindNumber)

      TurbSim_InitInput%WindFileName = InputFileData%TSFF_FileName

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_TurbSim_Init(TurbSim_InitInput, SumFileUnit, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid3D%RefHeight]

   case (BladedFF_WindNumber)

      Bladed_InitInput%TurbineID = InitInp%TurbineID
      IF ( InitInp%FixedWindFileRootName ) THEN ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
         IF ( InitInp%TurbineID == 0 ) THEN     ! .TRUE. for the FAST.Farm low-resolution domain
            InputFileData%BladedFF_FileName = TRIM(InputFileData%BladedFF_FileName)//TRIM(PathSep)//'Low'
         ELSE                                   ! FAST.Farm high-resolution domain(s)
            InputFileData%BladedFF_FileName = TRIM(InputFileData%BladedFF_FileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
         ENDIF
      ENDIF
      Bladed_InitInput%WindType           = BladedFF_WindNumber
      Bladed_InitInput%WindFileName       = TRIM(InputFileData%BladedFF_FileName)//'.wnd'
      Bladed_InitInput%TowerFileExist     =  InputFileData%BladedFF_TowerFile
      Bladed_InitInput%NativeBladedFmt    = .false.

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_Bladed_Init(Bladed_InitInput, SumFileUnit, Bladed_InitOutput, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)   
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid3D%RefHeight]
      
   case (BladedFF_Shr_WindNumber)

      Bladed_InitInput%WindType           = BladedFF_Shr_WindNumber
      Bladed_InitInput%TurbineID          = InitInp%TurbineID
      Bladed_InitInput%WindFileName       = InputFileData%BladedFF_FileName
      Bladed_InitInput%TowerFileExist     = .false.
      Bladed_InitInput%NativeBladedFmt    = .true.

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_Bladed_Init(Bladed_InitInput, SumFileUnit, Bladed_InitOutput, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)   
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Overwrite the values of PropagationDir and VFlowAngle with values from the native Bladed file
      InputFileData%PropagationDir = Bladed_InitOutput%PropagationDir
      InputFileData%VFlowAngle     = Bladed_InitOutput%VFlowAngle 

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid3D%RefHeight]

   case (HAWC_WindNumber)

      HAWC_InitInput%WindFileName(1) = InputFileData%HAWC_FileName_u
      HAWC_InitInput%WindFileName(2) = InputFileData%HAWC_FileName_v
      HAWC_InitInput%WindFileName(3) = InputFileData%HAWC_FileName_w
      HAWC_InitInput%nx = InputFileData%HAWC_nx
      HAWC_InitInput%ny = InputFileData%HAWC_ny
      HAWC_InitInput%nz = InputFileData%HAWC_nz
      HAWC_InitInput%dx = InputFileData%HAWC_dx
      HAWC_InitInput%dy = InputFileData%HAWC_dy
      HAWC_InitInput%dz = InputFileData%HAWC_dz
      HAWC_InitInput%G3D%RefHt = InputFileData%FF%RefHt
      HAWC_InitInput%G3D%ScaleMethod = InputFileData%FF%ScaleMethod
      HAWC_InitInput%G3D%SF = InputFileData%FF%SF
      HAWC_InitInput%G3D%SigmaF = InputFileData%FF%SigmaF
      HAWC_InitInput%G3D%URef = InputFileData%FF%URef
      HAWC_InitInput%G3D%WindProfileType = InputFileData%FF%WindProfileType
      HAWC_InitInput%G3D%PLExp = InputFileData%FF%PLExp
      HAWC_InitInput%G3D%Z0 = InputFileData%FF%Z0
      HAWC_InitInput%G3D%XOffset = InputFileData%FF%XOffset

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_HAWC_Init(HAWC_InitInput, SumFileUnit, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid3D%RefHeight]

   case (User_WindNumber)

      p%FlowField%FieldType = User_FieldType
      call IfW_User_Init(User_InitInput, SumFileUnit, p%FlowField%User, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%User%RefHeight]

   case (FDext_WindNumber)

      p%FlowField%FieldType = Grid4D_FieldType
      call IfW_Grid4D_Init(InitInp%FDext, p%FlowField%Grid4D, TmpErrStat, TmpErrMsg)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid4D%RefHeight]

   case (Point_WindNumber)

      p%FlowField%FieldType = Point_FieldType
      Points_InitInput%NumWindPoints = InitInp%NumWindPoints
      call IfW_Points_Init(Points_InitInput, p%FlowField%Points, TmpErrStat, TmpErrMsg)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      endif

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = 0.0_ReKi

   case default  
      call SetErrStat(ErrID_Fatal, ' Undefined wind type.', ErrStat, ErrMsg, RoutineName)
      return
   end select

   !----------------------------------------------------------------------------
   ! Initialization by Field Type
   !----------------------------------------------------------------------------

   ! Reset flag indicating that acceleration field is valid
   p%FlowField%AccFieldValid = .false.

   ! Copy flag for enabling cubic velocity interpolation
   p%FlowField%VelInterpCubic = InputFileData%VelInterpCubic

   ! If cubic velocity interpolation requested and linearization is performed,
   ! display message that cubic interpolation is incompatible with linearization
   ! and will be disabled
   if (p%FlowField%VelInterpCubic .and. InitInp%Linearize) then
      call WrScr("InflowWind: Cubic interpolation of wind velocity is disabled for linearization")
      p%FlowField%VelInterpCubic = .false.
   end if

   ! Set box exceed flag and index
   p%FlowField%Grid3D%BoxExceedAllowF = InitInp%BoxExceedAllowF
   p%FlowField%Grid3D%BoxExceedAllowIdx = huge(1_IntKi)
   if (InitInp%BoxExceedAllowF .and. (InitInp%BoxExceedAllowIdx <= InitInp%NumWindPoints)) then
      p%FlowField%Grid3D%BoxExceedAllowIdx = InitInp%BoxExceedAllowIdx
   end if

   ! Select based on field type
   select case (p%FlowField%FieldType)

   case (Uniform_FieldType)

      if (InitInp%OutputAccel .or. p%FlowField%VelInterpCubic) then
         call IfW_UniformField_CalcAccel(p%FlowField%Uniform, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
         p%FlowField%AccFieldValid = .true.
      end if

   case (Grid3D_FieldType)

      ! Calculate acceleration
      if (InitInp%OutputAccel .or. p%FlowField%VelInterpCubic) then
         call IfW_Grid3DField_CalcAccel(p%FlowField%Grid3D, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
         p%FlowField%AccFieldValid = .true.
      end if

      ! Calculate field average if box is allowed to be exceeded
      if (p%FlowField%Grid3D%BoxExceedAllowF .and. p%FlowField%Grid3D%BoxExceedAllowIdx > 0) then
         call IfW_Grid3DField_CalcVelAvgProfile(p%FlowField%Grid3D, p%FlowField%AccFieldValid, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      end if

   case default

      if (InitInp%OutputAccel) then
         call SetErrStat(ErrID_Fatal, "Acceleration not implemented for field type "// &
                         num2LStr(p%FlowField%FieldType), ErrStat, ErrMsg, RoutineName)
         return
      end if
      if (p%FlowField%VelInterpCubic) then
         call WrScr(' Cubic velocity interpolation not implemented for WindType '// &
                     num2LStr(InputFileData%WindType))
         p%FlowField%VelInterpCubic = .false.
      end if

   end select

   !----------------------------------------------------------------------------
   ! Set the p and OtherStates for InflowWind using the input file information.
   ! (set this after initializing modules so that we can use PropagationDir 
   ! and VFlowAng from native-Bladed files
   !----------------------------------------------------------------------------

   CALL InflowWind_SetParameters( InitInp, InputFileData, p, m, TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat>= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   ENDIF
   
   ! Allocate arrays for the WriteOutput
   CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat>= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   ENDIF
   y%WriteOutput = 0.0_ReKi
   
   CALL AllocAry( InitOutData%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat>= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   ENDIF

   CALL AllocAry( InitOutData%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat>= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   ENDIF

   InitOutData%WriteOutputHdr = p%OutParam(1:p%NumOuts)%Name
   InitOutData%WriteOutputUnt = p%OutParam(1:p%NumOuts)%Units     

   !----------------------------------------------------------------------------
   ! Linearization
   !----------------------------------------------------------------------------
      
   ! allocate and fill variables for linearization
   if (InitInp%Linearize) then

      ! If field is uniform and there is any nonzero upflow, return error
      ! Math needs work before this can be implemented
      if (p%FlowField%FieldType == Uniform_FieldType) then
         if (any(p%FlowField%Uniform%AngleV /= 0.0_ReKi)) then
            call SetErrStat(ErrID_Fatal, 'Upflow in uniform wind files must be 0 for linearization analysis in InflowWind.', ErrStat, ErrMsg, RoutineName)
            call Cleanup()
            return
         end if
      end if

      ! also need to add InputGuess%HubOrientation to the u%Linear items
      CALL AllocAry(InitOutData%LinNames_u, InitInp%NumWindPoints*3 + size(InputGuess%HubPosition) + 3 + NumExtendedInputs, 'LinNames_u', TmpErrStat, TmpErrMsg) ! add hub position, orientation(3) + extended inputs
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry(InitOutData%RotFrame_u, InitInp%NumWindPoints*3 + size(InputGuess%HubPosition) + 3 + NumExtendedInputs, 'RotFrame_u', TmpErrStat, TmpErrMsg)
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry(InitOutData%IsLoad_u, InitInp%NumWindPoints*3 + size(InputGuess%HubPosition) + 3 + NumExtendedInputs, 'IsLoad_u', TmpErrStat, TmpErrMsg)
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry(InitOutData%LinNames_y, InitInp%NumWindPoints*3 + size(y%DiskVel) + p%NumOuts, 'LinNames_y', TmpErrStat, TmpErrMsg)
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry(InitOutData%RotFrame_y, InitInp%NumWindPoints*3 + size(y%DiskVel) + p%NumOuts, 'RotFrame_y', TmpErrStat, TmpErrMsg)
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ENDIF
      
      do i=1,InitInp%NumWindPoints
         do j=1,3
            InitOutData%LinNames_y((i-1)*3+j) = UVW(j)//'-component inflow velocity at node '//trim(num2lstr(i))//', m/s'
            InitOutData%LinNames_u((i-1)*3+j) = XYZ(j)//'-component position of node '//trim(num2lstr(i))//', m'
         end do
      end do

      ! hub position
      Lin_Indx = InitInp%NumWindPoints*3
      do j=1,3
         InitOutData%LinNames_y(Lin_Indx+j) = 'average '//UVW(j)//'-component rotor-disk velocity, m/s'
         InitOutData%LinNames_u(Lin_Indx+j) = XYZ(j)//'-component position of moving hub, m'
      end do
      Lin_Indx = Lin_Indx + 3
      
      ! hub orientation angles
      do j=1,3
         InitOutData%LinNames_u(Lin_Indx+j) = XYZ(j)//' orientation of moving hub, rad'
      end do
      Lin_Indx = Lin_Indx + 3
      
      InitOutData%LinNames_u(Lin_Indx + 1) = 'Extended input: horizontal wind speed (steady/uniform wind), m/s'
      InitOutData%LinNames_u(Lin_Indx + 2) = 'Extended input: vertical power-law shear exponent, -'
      InitOutData%LinNames_u(Lin_Indx + 3) = 'Extended input: propagation direction, rad'         
      
      do i=1,p%NumOuts
         InitOutData%LinNames_y(i+3*InitInp%NumWindPoints+size(y%DiskVel)) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
      end do

      ! IfW inputs and outputs are in the global, not rotating frame
      InitOutData%RotFrame_u = .false. 
      InitOutData%RotFrame_y = .false. 

      InitOutData%IsLoad_u = .false. ! IfW inputs for linearization are not loads
      
   end if
               
   ! Set the version information in InitOutData
   InitOutData%Ver = IfW_Ver

   CALL CleanUp()

CONTAINS

   SUBROUTINE CleanUp()

      ! add in stuff that we need to dispose of here
      CALL InflowWind_DestroyInputFile( InputFileData,  TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

      ! Ignore error messages from InFileInfo destruction
      call NWTC_Library_DestroyFileInfoType( InFileInfo, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

         ! Close the summary file if we were writing one
      IF ( SumFileUnit > 0 ) THEN
         CALL InflowWind_CloseSumFile( SumFileUnit, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      ENDIF

   END SUBROUTINE CleanUp
   
END SUBROUTINE InflowWind_Init


!====================================================================================================
!> This routine takes an input dataset of type InputType which contains a position array of dimensions 3*n. It then calculates
!! and returns the output dataset of type OutputType which contains a corresponding velocity array of dimensions 3*n. The input
!! array contains XYZ triplets for each position of interest (first index is X/Y/Z for values 1/2/3, second index is the point
!! number to evaluate). The returned values in the OutputData are similar with U/V/W for the first index of 1/2/3.
!!
!! _Coordinate transformation:_
!! The coordinates passed in are copied to the PositionXYZPrime array, then rotated by -(p%PropagationDir) (radians) so
!! that the wind direction lies along the X-axis (all wind files are given this way).  The submodules are then called with
!! these PositionXYZPrime coordinates.
!!
!! After the calculation by the submodule, the PositionXYZPrime coordinate array is deallocated.  The returned VelocityUVW
!! array is then rotated by p%PropagationDir so that it now corresponds the the global coordinate UVW values for wind
!! with that direction.
!----------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_CalcOutput( Time, InputData, p, &
                              ContStates, DiscStates, ConstrStates, &   ! Framework required states -- empty in this case.
                              OtherStates, OutputData, m, ErrStat, ErrMsg )

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_CalcOutput"

   REAL(DbKi),                               INTENT(IN   )  :: Time              !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: InputData         !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: ContStates        !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: DiscStates        !< Discrete states at Time
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: ConstrStates      !< Constraint states at Time
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherStates       !< Other/optimization states at Time
   TYPE(InflowWind_OutputType),              INTENT(INOUT)  :: OutputData        !< Outputs computed at Time (IN for mesh reasons and data allocation)
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)                                           :: i
      ! Temporary variables for error handling
   INTEGER(IntKi)                                           :: TmpErrStat
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg            ! temporary error message

   ErrStat  = ErrID_None
   ErrMsg   = ""

      ! Allocate the velocity array to get out
   IF ( .NOT. ALLOCATED(OutputData%VelocityUVW) ) THEN
      CALL AllocAry( OutputData%VelocityUVW, 3, SIZE(InputData%PositionXYZ,DIM=2), &
                  "Velocity array returned from IfW_CalcOutput", TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN
   ELSEIF ( SIZE(InputData%PositionXYZ,DIM=2) /= SIZE(OutputData%VelocityUVW,DIM=2) ) THEN
      CALL SetErrStat( ErrID_Fatal," Programming error: Position and Velocity arrays are not sized the same.",  &
            ErrStat, ErrMsg, RoutineName)
   ENDIF

   !-----------------------------
   ! Outputs: OutputData%VelocityUVW
   !-----------------------------
   CALL CalculateOutput( Time, InputData, p, ContStates, DiscStates, ConstrStates, & 
                     OtherStates, OutputData, m, .TRUE., TmpErrStat, TmpErrMsg )      
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )

   !-----------------------------
   ! Output: OutputData%DiskVel
   !-----------------------------
   CALL InflowWind_GetSpatialAverage( Time, InputData, p, ContStates, DiscStates, ConstrStates, &
                     OtherStates, m, OutputData%DiskVel, TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      
   !-----------------------------
   ! Output: OutputData%HubVel
   !-----------------------------
   CALL InflowWind_GetHubValues( Time, InputData, p, ContStates, DiscStates, ConstrStates, &
                     OtherStates, m, OutputData%HubVel, TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      
   !-----------------------------
   ! Outputs: OutputData%lidar%LidSpeed and OutputData%lidar%WtTrunc
   !-----------------------------
      
      ! return sensor values
   IF (p%lidar%SensorType /= SensorType_None) THEN
         
      CALL Lidar_CalcOutput(Time, InputData, p, &
                           ContStates, DiscStates, ConstrStates, OtherStates, &  
                           OutputData, m, TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
         
   END IF      
       
      
   !-----------------------------
   ! Outputs: OutputData%WriteOutput from m%AllOuts and OutputData%lidar%LidSpeed
   !-----------------------------

   CALL SetAllOuts( p, OutputData, m, TmpErrStat, TmpErrMsg ) 
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
   
      ! Map to the outputs
   DO I = 1,p%NumOuts  ! Loop through all selected output channels
      OutputData%WriteOutput(I) = p%OutParam(I)%SignM * m%AllOuts( p%OutParam(I)%Indx )
   ENDDO             ! I - All selected output channels


END SUBROUTINE InflowWind_CalcOutput



!====================================================================================================
!> Clean up the allocated variables and close all open files.  Reset the initialization flag so
!! that we have to reinitialize before calling the routines again.
SUBROUTINE InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                           y, m, ErrStat, ErrMsg )

   TYPE(InflowWind_InputType),               INTENT(INOUT)  :: InputData         !< Input data for initialization
   TYPE(InflowWind_ParameterType),           INTENT(INOUT)  :: p                 !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(INOUT)  :: ContStates        !< Continuous states
   TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: DiscStates        !< Discrete states
   TYPE(InflowWind_ConstraintStateType),     INTENT(INOUT)  :: ConstrStateGuess  !< Guess of the constraint states
   TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherStates       !< Other/optimization states
   TYPE(InflowWind_OutputType),              INTENT(INOUT)  :: y                 !< Output data
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   INTEGER( IntKi ),                         INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< error message

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_End"

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Reset the wind type so that the initialization routine must be called
   p%WindType      = Undef_WindNumber

   ! Destroy all inflow wind derived types
   CALL InflowWind_DestroyInput( InputData, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyParam( p, ErrStat, ErrMsg, DeallocatePointers=.true. )         
   CALL InflowWind_DestroyContState( ContStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyDiscState( DiscStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyConstrState( ConstrStateGuess, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyOtherState( OtherStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyOutput( y, ErrStat, ErrMsg )                     
   CALL InflowWind_DestroyMisc( m, ErrStat, ErrMsg )                     

END SUBROUTINE InflowWind_End


!====================================================================================================
! The following routines were added to satisfy the framework, but do nothing useful.
!====================================================================================================
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other 
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE InflowWind_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                            INTENT(IN   ) :: t               !< Current simulation time in seconds
   INTEGER(IntKi),                        INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
   TYPE(InflowWind_InputType),            INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes (output only for mesh record-keeping in ExtrapInterp routine)
   REAL(DbKi),                            INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
   TYPE(InflowWind_ParameterType),        INTENT(IN   ) :: p               !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!    Output: Continuous states at t + Interval
   TYPE(InflowWind_DiscreteStateType),    INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!    Output: Discrete states at t  + Interval
   TYPE(InflowWind_ConstraintStateType),  INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
   TYPE(InflowWind_OtherStateType),       INTENT(INOUT) :: OtherState      !< Other states: Other states at t;
                                                                           !!   Output: Other states at t + Interval
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT) :: m               !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                        INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   x%DummyContState     = 0.0_ReKi
   xd%DummyDiscState    = 0.0_ReKi
   z%DummyConstrState   = 0.0_ReKi
      
   RETURN


END SUBROUTINE InflowWind_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE InflowWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(InflowWind_ContinuousStateType),     INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Compute the first time derivatives of the continuous states here:

   dxdt%DummyContState = 0.0_ReKi


END SUBROUTINE InflowWind_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE InflowWind_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                           !! Output: Discrete states at Time + Interval
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Update discrete states here:

   ! StateData%DiscState =

END SUBROUTINE InflowWind_UpdateDiscState

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE InflowWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(InflowWind_ConstraintStateType),     INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                           !! the input values described above
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Solve for the constraint states here:

   z_residual%DummyConstrState = 0

END SUBROUTINE InflowWind_CalcConstrStateResidual

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in IfW_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE InflowWind_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) 
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) 
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) 
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z)  
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
 
      ! local variables: 
   INTEGER(IntKi)                                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                                           :: ErrMsg2            ! temporary error message
   CHARACTER(*), PARAMETER                                        :: RoutineName = 'InflowWind_JacobianPInput'
      
   
   REAL(R8Ki)                                                     :: local_dYdu(3,3+NumExtendedInputs)
   integer                                                        :: i, n
   integer                                                        :: i_start, i_end  ! indices for input/output start and end
   integer                                                        :: node, comp
   integer                                                        :: n_inputs
   integer                                                        :: n_outputs
   integer                                                        :: i_ExtendedInput_start
   integer                                                        :: i_WriteOutput
   
      
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdu ) ) THEN

      n_outputs = SIZE(u%PositionXYZ)+p%NumOuts + size(y%DiskVel)
      n_inputs  = SIZE(u%PositionXYZ)+size(u%HubPosition) + 3 + NumExtendedInputs ! need to add 3 for u%HubOrientation
      i_ExtendedInput_start = n_inputs - NumExtendedInputs + 1 ! index for extended inputs starts 2 from end (encompasses 3 values: V, VShr, PropDir)
      i_WriteOutput         = n_outputs - p%NumOuts ! index for where write outputs begin is i_WriteOutput + 1
      
      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
         
         ! outputs are all velocities at all positions plus the WriteOutput values
         !
      if (.not. ALLOCATED(dYdu)) then
         CALL AllocAry( dYdu, n_outputs, n_inputs, 'dYdu', ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
         
         
      SELECT CASE ( p%WindType )
      CASE (Steady_WindNumber, Uniform_WindNumber)

            ! note that we are including the propagation direction in the analytical derivative calculated
            ! inside IfW_UniformWind_JacobianPInput, so no need to transform input position vectors first
         
         dYdu = 0.0_R8Ki ! initialize all non-diagonal entries to zero (position of node effects the output of only that node) 
         
         n = SIZE(u%PositionXYZ,2)
            ! these are the positions used in the module coupling
         do i=1,n
            ! note that p%FlowField%RotToWind(1,1) = cos(p%PropagationDir) and p%FlowField%RotToWind(2,1) = sin(p%PropagationDir), which are the
            ! values we need to compute the jacobian.
!!!FIX ME with the propagation values!!!!         
            call IfW_UniformWind_JacobianPInput( p%FlowField%Uniform, t, u%PositionXYZ(:,i), p%FlowField%RotToWind(1,1), p%FlowField%RotToWind(2,1), local_dYdu )
            
            i_end  = 3*i
            i_start= i_end - 2
            
            dYdu(i_start:i_end,i_start:i_end) = local_dYdu(:,1:3)
            
            dYdu(i_start:i_end, i_ExtendedInput_start:) = local_dYdu(:,4:6) ! extended inputs
            
         end do
         
         
         ! see InflowWind_GetSpatialAverage():
         
         ! location of y%DiskAvg
         i_start = 3*n + 1
         i_end   = i_start + 2
         
         dYdu(i_start:i_end,:) = 0.0_R8Ki ! initialize because we're going to create averages
         
         do i=1,IfW_NumPtsAvg
            m%u_Avg%PositionXYZ(:,i) = matmul(u%HubOrientation,p%PositionAvg(:,i)) + u%HubPosition
!!!FIX ME with the propagation values!!!!         
            call IfW_UniformWind_JacobianPInput( p%FlowField%Uniform, t, m%u_Avg%PositionXYZ(:,i), p%FlowField%RotToWind(1,1), p%FlowField%RotToWind(2,1), local_dYdu )
         
            ! y%DiskAvg has the same index as u%HubPosition
            ! Also note that partial_(m%u_Avg%PositionXYZ) / partial_(u%HubPosition) is identity, so we can skip that part of the chain rule for these derivatives:
            dYdu(i_start:i_end,i_start:i_end)           = dYdu(i_start:i_end, i_start:i_end)          + local_dYdu(:,1:3)
            dYdu(i_start:i_end, i_ExtendedInput_start:) = dYdu(i_start:i_end, i_ExtendedInput_start:) + local_dYdu(:,4:6) ! extended inputs
         end do
         dYdu(i_start:i_end,i_start:i_end)           = dYdu(i_start:i_end, i_start:i_end)          / REAL(IfW_NumPtsAvg,R8Ki)
         dYdu(i_start:i_end,i_ExtendedInput_start:)  = dYdu(i_start:i_end, i_ExtendedInput_start:) / REAL(IfW_NumPtsAvg,R8Ki)
!FIX ME:
         ! need to calculate dXYZdHubOrient = partial_(m%u_Avg%PositionXYZ) / partial_(u%HubOrientation)
         !dYdu(i_start:i_end,(i_start+3):(i_end+3)) = matmul( dYdu(i_start:i_end,i_start:i_end), dXYZdHubOrient )


            ! these are the InflowWind WriteOutput velocities (and note that we may not have all of the components of each point) 
         ! they do not depend on the inputs, so the derivatives w.r.t. X, Y, Z are all zero
         do i=1, p%NumOuts
            node  = p%OutParamLinIndx(1,i) ! output node
            comp  = p%OutParamLinIndx(2,i) ! component of output node

            if (node > 0) then
!!!FIX ME with the propagation values!!!!         
               call IfW_UniformWind_JacobianPInput( p%FlowField%Uniform, t, p%WindViXYZ(:,node), p%FlowField%RotToWind(1,1), p%FlowField%RotToWind(2,1), local_dYdu )
            else
               local_dYdu = 0.0_R8Ki
               comp = 1
            end if
            
            dYdu(i_WriteOutput+i, i_ExtendedInput_start:) = p%OutParam(i)%SignM * local_dYdu( comp , 4:6)
         end do

      CASE DEFAULT

      END SELECT
                  
   END IF

   IF ( PRESENT( dXdu ) ) THEN
      if (allocated(dXdu)) deallocate(dXdu) 
   END IF

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu) 
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu) 
   END IF


END SUBROUTINE InflowWind_JacobianPInput
!..................................................................................................................................
!> Routine to compute the Jacobians of the output (Y) function with respect to the inputs (u). The partial 
!! derivative dY/du is returned. This submodule does not follow the modularization framework.
SUBROUTINE IfW_UniformWind_JacobianPInput(UF, t, Position, CosPropDir, SinPropDir, dYdu)
   USE IfW_FlowField, only : UniformField_InterpLinear, UniformField_InterpCubic

   TYPE(UniformFieldType),    INTENT(IN   )  :: UF                !< Uniform field derived type
   REAL(DbKi),                INTENT(IN   )  :: t                 !< Current simulation time in seconds
   REAL(ReKi),                INTENT(IN   )  :: Position(3)       !< XYZ Position at which to find velocity (operating point)
   REAL(ReKi),                INTENT(IN   )  :: CosPropDir        !< cosine of InflowWind propagation direction
   REAL(ReKi),                INTENT(IN   )  :: SinPropDir        !< sine of InflowWind propagation direction
   REAL(R8Ki),                INTENT(INOUT)  :: dYdu(3,6)         !< Partial derivatives of output functions (Y) with respect to the inputs (u)

   TYPE(UniformField_Interp)                 :: op                ! interpolated values of InterpParams
   REAL(R8Ki)                                :: RotatePosition(3) !< rotated position
   REAL(R8Ki)                                :: dVhdx             ! temporary value to hold partial v_h partial X   
   REAL(R8Ki)                                :: dVhdy             ! temporary value to hold partial v_h partial Y   
   REAL(R8Ki)                                :: dVhdz             ! temporary value to hold partial v_h partial Z   
   REAL(R8Ki)                                :: tmp_du            ! temporary value to hold calculations that are part of multiple components   
   REAL(R8Ki)                                :: tmp_dv            ! temporary value to hold calculations that are part of multiple components   
   REAL(R8Ki)                                :: dVhdPD            ! temporary value to hold partial v_h partial propagation direction
   REAL(R8Ki)                                :: dVhdV             ! temporary value to hold partial v_h partial V   
   REAL(R8Ki)                                :: Vh                ! temporary value to hold v_h    
   REAL(R8Ki)                                :: dVhdVShr          ! temporary value to hold partial v_h partial VShr   
   REAL(R8Ki)                                :: zr 

   if ( Position(3) < 0.0_ReKi .or. EqualRealNos(Position(3), 0.0_ReKi)) then
      dYdu = 0.0_R8Ki
      return
   end if      
      
   !-------------------------------------------------------------------------------------------------
   !> 1. Interpolate uniform field to get values at operating point
   !-------------------------------------------------------------------------------------------------

   op = UniformField_InterpLinear(UF, t)
   
   RotatePosition(1) = Position(1)*cosPropDir - Position(2)*sinPropDir
   RotatePosition(2) = Position(1)*sinPropDir + Position(2)*cosPropDir
   RotatePosition(3) = Position(3)
   
   !-------------------------------------------------------------------------------------------------
   !> 2. Calculate \f$ \frac{\partial Y_{Output \, Equations}}{\partial u_{inputs}} = \begin{bmatrix}
   !! \frac{\partial Vt_u}{\partial X} & \frac{\partial Vt_u}{\partial Y} & \frac{\partial Vt_u}{\partial Z} \\
   !! \frac{\partial Vt_v}{\partial X} & \frac{\partial Vt_v}{\partial Y} & \frac{\partial Vt_v}{\partial Z} \\
   !! \frac{\partial Vt_w}{\partial X} & \frac{\partial Vt_w}{\partial Y} & \frac{\partial Vt_w}{\partial Z} \\
   !! \end{bmatrix} \f$
   !-------------------------------------------------------------------------------------------------

   zr = RotatePosition(3)/UF%RefHeight
   tmp_du = op%VelH * op%ShrH / UF%RefLength * CosPropDir
   dVhdx  = tmp_du * op%SinAngleH
   dVhdy  = tmp_du * op%CosAngleH   
   dVhdz  = op%VelH * ( op%ShrV / UF%RefHeight * zr**(op%ShrV-1.0_R8Ki) + op%LinShrV/UF%RefLength)

   dVhdV = ( ( RotatePosition(3)/UF%RefHeight ) ** op%ShrV &                                             ! power-law wind shear
             + ( op%ShrH   * ( RotatePosition(2) * op%CosAngleH + RotatePosition(1) * op%SinAngleH ) &   ! horizontal linear shear
             +  op%LinShrV * ( RotatePosition(3) - UF%RefHeight ) )/UF%RefLength  )                      ! vertical linear shear   
   Vh = op%VelH * dVhdV + op%VelGust

   dVhdVShr = op%VelH * zr**op%ShrV * log(zr)
   dVhdPD   = op%VelH * op%ShrH / UF%RefLength * ( RotatePosition(1) * op%CosAngleH - RotatePosition(2) * op%SinAngleH )

   tmp_du =  CosPropDir*op%CosAngleH  - SinPropDir*op%SinAngleH
   tmp_dv = -SinPropDir*op%CosAngleH  - CosPropDir*op%SinAngleH

      !> \f$ \frac{\partial Vt_u}{\partial X} = \left[\cos(PropagationDir)\cos(Delta) - \sin(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \sin(Delta) \cos(PropagationDir) \f$
   dYdu(1,1) = tmp_du*dVhdx
      !> \f$ \frac{\partial Vt_v}{\partial X} = \left[-\sin(PropagationDir)\cos(Delta) - \cos(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \sin(Delta) \cos(PropagationDir) \f$
   dYdu(2,1) = tmp_dv*dVhdx
      !> \f$ \frac{\partial Vt_w}{\partial X} = 0 \f$
   dYdu(3,1) = 0.0_R8Ki

      !> \f$ \frac{\partial Vt_u}{\partial Y} = \left[\cos(PropagationDir)\cos(Delta) - \sin(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \cos(Delta) \cos(PropagationDir) \f$
   dYdu(1,2) = tmp_du*dVhdy
      !> \f$ \frac{\partial Vt_v}{\partial Y} = \left[-\sin(PropagationDir)\cos(Delta) - \cos(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \cos(Delta) \cos(PropagationDir) \f$
   dYdu(2,2) = tmp_dv*dVhdy
      !> \f$ \frac{\partial Vt_w}{\partial Y} = 0 \f$
   dYdu(3,2) = 0.0_R8Ki

      !> \f$ \frac{\partial Vt_u}{\partial Z} = \left[\cos(PropagationDir)\cos(Delta) - \sin(PropagationDir)\sin(Delta) \right]
      !! V \, \left[ \frac{V_{shr}}{Z_{ref}} \left( \frac{Z}{Z_{ref}} \right) ^ {V_{shr}-1} + \frac{V_{LinShr}}{RefWid} \right] \f$
   dYdu(1,3) = tmp_du*dVhdz 
      !> \f$ \frac{\partial Vt_v}{\partial Z} = \left[-\sin(PropagationDir)\cos(Delta) - \cos(PropagationDir)\sin(Delta) \right]
      !! V \, \left[ \frac{V_{shr}}{Z_{ref}} \left( \frac{Z}{Z_{ref}} \right) ^ {V_{shr}-1} + \frac{V_{LinShr}}{RefWid} \right] \f$      
   dYdu(2,3) = tmp_dv*dVhdz
      !> \f$ \frac{\partial Vt_w}{\partial Z} = 0 \f$
   dYdu(3,3) = 0.0_R8Ki

      ! \f$ \frac{\partial Vt_u}{\partial V} =  \f$
   dYdu(1,4) = tmp_du*dVhdV      
      ! \f$ \frac{\partial Vt_v}{\partial V} =  \f$
   dYdu(2,4) = tmp_dv*dVhdV
      !> \f$ \frac{\partial Vt_w}{\partial V} = 0 \f$
   dYdu(3,4) = 0.0_R8Ki

      ! \f$ \frac{\partial Vt_u}{\partial VShr} =  \f$
   dYdu(1,5) = tmp_du*dVhdVShr
      ! \f$ \frac{\partial Vt_v}{\partial VShr} =  \f$
   dYdu(2,5) = tmp_dv*dVhdVShr
      !> \f$ \frac{\partial Vt_w}{\partial VShr} = 0 \f$
   dYdu(3,5) = 0.0_R8Ki

      ! \f$ \frac{\partial Vt_u}{\partial PropDir} =  \f$
   dYdu(1,6) = tmp_dv*Vh + tmp_du*dVhdPD
      ! \f$ \frac{\partial Vt_v}{\partial PropDir} =  \f$
   dYdu(2,6) = -tmp_du*Vh + tmp_dv*dVhdPD
      !> \f$ \frac{\partial Vt_w}{\partial PropDir} = 0 \f$
   dYdu(3,6) = 0.0_R8Ki

END SUBROUTINE IfW_UniformWind_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE InflowWind_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdx.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions
                                                                                !!   (Y) with respect to the continuous
                                                                                !!   states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state
                                                                                !!   functions (X) with respect to
                                                                                !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state
                                                                                !!   functions (Xd) with respect to
                                                                                !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state
                                                                                !!   functions (Z) with respect to
                                                                                !!   the continuous states (x) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''



   IF ( PRESENT( dYdx ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:

      ! allocate and set dYdx

   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate and set dXdx

   END IF

   IF ( PRESENT( dXddx ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the continuous states (x) here:

      ! allocate and set dXddx

   END IF

   IF ( PRESENT( dZdx ) ) THEN


      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the continuous states (x) here:

      ! allocate and set dZdx

   END IF


END SUBROUTINE InflowWind_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
SUBROUTINE InflowWind_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdxd.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                                !!  (Y) with respect to the discrete
                                                                                !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                                !!   functions (X) with respect to the
                                                                                !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddxd(:,:) !< Partial derivatives of discrete state
                                                                                !!   functions (Xd) with respect to the
                                                                                !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                                !!   functions (Z) with respect to the
                                                                                !!   discrete states (xd) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdxd ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:

      ! allocate and set dYdxd

   END IF

   IF ( PRESENT( dXdxd ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:

      ! allocate and set dXdxd

   END IF

   IF ( PRESENT( dXddxd ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:

      ! allocate and set dXddxd

   END IF

   IF ( PRESENT( dZdxd ) ) THEN

      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:

      ! allocate and set dZdxd

   END IF


END SUBROUTINE InflowWind_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
SUBROUTINE InflowWind_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdz.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output
                                                                                !!  functions (Y) with respect to the
                                                                                !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous
                                                                                !!  state functions (X) with respect to
                                                                                !!  the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state
                                                                                !!  functions (Xd) with respect to the
                                                                                !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint
                                                                                !! state functions (Z) with respect to
                                                                                !!  the constraint states (z) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdz ) ) THEN

         ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z) here:

      ! allocate and set dYdz

   END IF

   IF ( PRESENT( dXdz ) ) THEN

         ! Calculate the partial derivative of the continuous state functions (X) with respect to the constraint states (z) here:

      ! allocate and set dXdz

   END IF

   IF ( PRESENT( dXddz ) ) THEN

         ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the constraint states (z) here:

      ! allocate and set dXddz

   END IF

   IF ( PRESENT( dZdz ) ) THEN

         ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z) here:

      ! allocate and set dZdz

   END IF


END SUBROUTINE InflowWind_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE InflowWind_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),           INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),       INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType), INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),   INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType), INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),      INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),          INTENT(IN   )           :: y          !< Output at operating point
   TYPE(InflowWind_MiscVarType),         INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states

   
   INTEGER(IntKi)                                    :: index, i, j
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'InflowWind_GetOP'
   

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( u_op ) ) THEN
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, size(u%PositionXYZ) + size(u%HubPosition) + 3 + NumExtendedInputs, 'u_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
         
      index = 0
      do i=1,size(u%PositionXYZ,2)
         do j=1,size(u%PositionXYZ,1)
            index = index + 1 !(i-1)*size(u%PositionXYZ,1)+j
            u_op(index) = u%PositionXYZ(j,i)
         end do            
      end do  
      
      do i=1,3
         index = index + 1
         u_op(index) = u%HubPosition(i)
      end do
      
      u_op((index+1):(index+3)) = EulerExtract(u%HubOrientation)
      index = index + 3
      
      call IfW_UniformWind_GetOP( p%FlowField%Uniform, t, p%FlowField%VelInterpCubic, u_op(index+1:index+2) )
      u_op(index + 3) = p%FlowField%PropagationDir
      
   END IF

   IF ( PRESENT( y_op ) ) THEN
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, size(u%PositionXYZ)+p%NumOuts+3, 'y_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
      index = 0
      do i=1,size(u%PositionXYZ,2)
         do j=1,size(u%PositionXYZ,1)
            index = index + 1 !(i-1)*size(u%PositionXYZ,1)+j
            y_op(index) = y%VelocityUVW(j,i)
         end do
      end do
         
      do j=1,size(y%DiskVel)
         index = index + 1
         y_op(index) = y%DiskVel(j)
      end do
      
      do i=1,p%NumOuts
         y_op(i+index) = y%WriteOutput( i )
      end do
         
   END IF

   IF ( PRESENT( x_op ) ) THEN

   END IF

   IF ( PRESENT( dx_op ) ) THEN

   END IF

   IF ( PRESENT( xd_op ) ) THEN

   END IF
   
   IF ( PRESENT( z_op ) ) THEN

   END IF

END SUBROUTINE InflowWind_GetOP

END MODULE InflowWind
