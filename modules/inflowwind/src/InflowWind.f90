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
!> InflowWind is used to read and process the (undisturbed) inflow winds.  It must be initialized
!! using InflowWind_Init() with the name of the file, the file type, and possibly reference height and
!! width (depending on the type of wind file being used).  This module calls appropriate routines
!! in the wind modules so that the type of wind becomes seamless to the user.  InflowWind_End()
!! should be called when the program has finshed.
!!
!! Data are assumed to be in units of meters and seconds.  Z is measured from the ground (NOT the hub!).
!!
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
   integer,        parameter            :: NumExtendedIO = 3 ! Number of extended inputs or outputs (same set): HWindSpeed, PlExp, PropDir

   ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: InflowWind_Init                                   !< Initialization routine
   PUBLIC :: InflowWind_CalcOutput                             !< Calculate the wind velocities
   PUBLIC :: InflowWind_End                                    !< Ending routine (includes clean up)

   ! These routines compute Jacobians; only dYdu is defined.
   PUBLIC :: InflowWind_JacobianPInput
   PUBLIC :: InflowWind_JacobianPContState
   PUBLIC :: InflowWind_JacobianPDiscState
   PUBLIC :: InflowWind_JacobianPConstrState
   PUBLIC :: InflowWind_GetOP

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

   IF ( InitInp%FilePassingMethod == 0_IntKi ) THEN    ! Normal calling with an input file

      CALL GetPath( InitInp%InputFileName, PriPath )
      CALL ProcessComFile( InitInp%InputFileName, InFileInfo, TmpErrStat, TmpErrMsg ); if (Failed()) return
      ! For diagnostic purposes, the following can be used to display the contents
      ! of the InFileInfo data structure.
      ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.

      CALL InflowWind_ParseInputFileInfo( InputFileData,  InFileInfo, PriPath, InitInp%InputFileName, EchoFileName, InitInp%FixedWindFileRootName, InitInp%TurbineID, TmpErrStat, TmpErrMsg ); if (Failed()) return

   ELSEIF ( InitInp%FilePassingMethod == 1_IntKi ) THEN        ! passing the FileInfoType structure

      CALL GetPath( InitInp%InputFileName, PriPath )           ! in case a summary file is written
      CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileInfo, InFileInfo, MESH_NEWCOPY, TmpErrStat, TmpErrMsg ); if (Failed()) return
      CALL InflowWind_ParseInputFileInfo( InputFileData,  InFileInfo, PriPath, InitInp%InputFileName, EchoFileName, InitInp%FixedWindFileRootName, InitInp%TurbineID, TmpErrStat, TmpErrMsg ); if (Failed()) return

   ELSEIF ( InitInp%FilePassingMethod == 2_IntKi ) THEN        ! passing the InputFileData structure

      CALL InflowWind_CopyInputFile( InitInp%PassedFileData, InputFileData, MESH_NEWCOPY, TmpErrStat, TmpErrMsg ); if (Failed()) return

   ENDIF

   ! If wind is Grid4D from FAST.Farm, set input file values
   IF (InitInp%Use4Dext) then
      InputFileData%WindType = FDext_WindNumber      
      InputFileData%PropagationDir = 0.0_ReKi ! wind is in XYZ coordinates (already rotated if necessary), so don't rotate it again
      InputFileData%VFlowAngle = 0.0_ReKi
      InputFileData%VelInterpCubic = .false.
   END IF


      ! Validate the InflowWind input file information.
   CALL InflowWind_ValidateInput( InitInp, InputFileData, TmpErrStat, TmpErrMsg ); if (Failed()) return

      
   ! Disable Lidar if not allowed (FAST.Farm doesn't allow this)
   if (.not. InitInp%LidarEnabled) then
      if (p%lidar%SensorType /= SensorType_None) then
         call WrScr('  WARNING: LiDAR cannot be used with this instance of InflowWind (not usable with FAST.Farm).')
         call WrScr('   --> Disabling LiDAR.')
         p%lidar%SensorType = SensorType_None
      end if
   endif

   ! initialize sensor data:
   CALL Lidar_Init( InitInp, InputFileData, InputGuess, p, y, m, TimeInterval, TmpErrStat, TmpErrMsg )
      if (Failed()) return

      ! If a summary file was requested, open it.
   IF ( InputFileData%SumPrint ) THEN

         ! Open the summary file and write some preliminary info to it
      CALL InflowWind_OpenSumFile( SumFileUnit, SumFileName, IfW_Ver, InputFileData%WindType, TmpErrStat, TmpErrMsg ); if (Failed()) return
   ELSE
      SumFileUnit =  -1_IntKi       ! So that we don't try to write to something.  Used as indicator in submodules.
   ENDIF

   ! Allocate the array for passing points
   CALL AllocAry( InputGuess%PositionXYZ, 3, InitInp%NumWindPoints, "Array of positions at which to find wind velocities", TmpErrStat, TmpErrMsg ); if (Failed()) return
   InputGuess%PositionXYZ = 0.0_ReKi
   InputGuess%HubPosition = InitInp%HubPosition
   CALL Eye(InputGuess%HubOrientation,TmpErrStat,TmpErrMsg); if (Failed()) return

   ! Allocate the array for passing velocities out
   CALL AllocAry( y%VelocityUVW, 3, InitInp%NumWindPoints, "Array of wind velocities returned by InflowWind", TmpErrStat, TmpErrMsg ); if (Failed()) return
   y%VelocityUVW = 0.0_ReKi

   ! If requested, allocate the array for passing accelerations out
   IF ( InitInp%OutputAccel ) THEN
      CALL AllocAry( y%AccelUVW, 3, InitInp%NumWindPoints, "Array of wind accelerations returned by InflowWind", TmpErrStat, TmpErrMsg ); if (Failed()) return
      y%AccelUVW = 0.0_ReKi
   ENDIF

   !----------------------------------------------------------------------------
   ! Set flow field input data based on wind type
   !----------------------------------------------------------------------------

   ! If flowfield is allocated, deallocate and allocate again to clear old data
   if (associated(p%FlowField)) deallocate(p%FlowField)
   allocate(p%FlowField)

   ! Associate initialization output to flow field
   InitOutData%FlowField => p%FlowField

   ! Initialize mean wind speed to a very large number
   InitOutData%WindFileInfo%MWS = HUGE(InitOutData%WindFileInfo%MWS)

   ! Switch based on the wind type specified in the input file
   select case(InputFileData%WindType)

   case (Steady_WindNumber)

      Steady_InitInput%HWindSpeed = InputFileData%Steady_HWindSpeed
      Steady_InitInput%RefHt = InputFileData%Steady_RefHt
      Steady_InitInput%PLExp = InputFileData%Steady_PLexp

      p%FlowField%FieldType = Uniform_FieldType
      call IfW_SteadyWind_Init(Steady_InitInput, SumFileUnit, p%FlowField%Uniform, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Uniform%RefHeight]

   case (Uniform_WindNumber)

      Uniform_InitInput%WindFileName   = InputFileData%Uniform_FileName
      Uniform_InitInput%RefHt          = InputFileData%Uniform_RefHt
      Uniform_InitInput%RefLength      = InputFileData%Uniform_RefLength
      Uniform_InitInput%PropagationDir = InputFileData%PropagationDir

      p%FlowField%FieldType = Uniform_FieldType
      call IfW_UniformWind_Init(Uniform_InitInput, SumFileUnit, p%FlowField%Uniform, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Uniform%RefHeight]

   case (TSFF_WindNumber)

      TurbSim_InitInput%WindFileName = InputFileData%TSFF_FileName

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_TurbSim_Init(TurbSim_InitInput, SumFileUnit, p%FlowField%Grid3D, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

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
      call IfW_Bladed_Init(Bladed_InitInput, SumFileUnit, Bladed_InitOutput, p%FlowField%Grid3D, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid3D%RefHeight]

   case (BladedFF_Shr_WindNumber)

      Bladed_InitInput%WindType           = BladedFF_Shr_WindNumber
      Bladed_InitInput%TurbineID          = InitInp%TurbineID
      Bladed_InitInput%WindFileName       = InputFileData%BladedFF_FileName
      Bladed_InitInput%TowerFileExist     = .false.
      Bladed_InitInput%NativeBladedFmt    = .true.

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_Bladed_Init(Bladed_InitInput, SumFileUnit, Bladed_InitOutput, p%FlowField%Grid3D, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

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
      call IfW_HAWC_Init(HAWC_InitInput, SumFileUnit, p%FlowField%Grid3D, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid3D%RefHeight]

   case (User_WindNumber)

      p%FlowField%FieldType = User_FieldType
      call IfW_User_Init(User_InitInput, SumFileUnit, p%FlowField%User, InitOutData%WindFileInfo, TmpErrStat, TmpErrMsg); if (Failed()) return

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%User%RefHeight]

   case (FDext_WindNumber)

      p%FlowField%FieldType = Grid4D_FieldType
      call IfW_Grid4D_Init(InitInp%FDext, p%FlowField%Grid4D, TmpErrStat, TmpErrMsg); if (Failed()) return

      ! Set reference position for wind rotation
      p%FlowField%RefPosition = [0.0_ReKi, 0.0_ReKi, p%FlowField%Grid4D%RefHeight]

   case (Point_WindNumber)

      p%FlowField%FieldType = Point_FieldType
      Points_InitInput%NumWindPoints = InitInp%NumWindPoints
      call IfW_Points_Init(Points_InitInput, p%FlowField%Points, TmpErrStat, TmpErrMsg); if (Failed()) return

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

   ! Select based on field type
   select case (p%FlowField%FieldType)

   case (Uniform_FieldType)

      if (InitInp%OutputAccel .or. p%FlowField%VelInterpCubic) then
         call IfW_UniformField_CalcAccel(p%FlowField%Uniform, TmpErrStat, TmpErrMsg); if (Failed()) return
         p%FlowField%AccFieldValid = .true.
      end if

   case (Grid3D_FieldType)

      ! Calculate acceleration
      if (InitInp%OutputAccel .or. p%FlowField%VelInterpCubic) then
         call IfW_Grid3DField_CalcAccel(p%FlowField%Grid3D, TmpErrStat, TmpErrMsg); if (Failed()) return
         p%FlowField%AccFieldValid = .true.
      end if

      ! If input requested points to exceed box or if lidar is enabled,
      ! set flag to allow box to be exceeded
      p%FlowField%Grid3D%BoxExceedAllow = &
         InitInp%BoxExceedAllow .or. (p%lidar%SensorType /= SensorType_None)

      ! Calculate field average if box is allowed to be exceeded
      if (p%FlowField%Grid3D%BoxExceedAllow) then
         call IfW_Grid3DField_CalcVelAvgProfile(p%FlowField%Grid3D, p%FlowField%AccFieldValid, TmpErrStat, TmpErrMsg); if (Failed()) return
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

   CALL InflowWind_SetParameters( InitInp, InputFileData, p, m, TmpErrStat, TmpErrMsg ); if (Failed()) return
   
   ! Allocate arrays for the WriteOutput
   CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', TmpErrStat, TmpErrMsg ); if (Failed()) return
   y%WriteOutput = 0.0_ReKi
   
   CALL AllocAry( InitOutData%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', TmpErrStat, TmpErrMsg ); if (Failed()) return
   CALL AllocAry( InitOutData%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', TmpErrStat, TmpErrMsg ); if (Failed()) return

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
      CALL AllocAry(InitOutData%LinNames_u, NumExtendedIO, 'LinNames_u', TmpErrStat, TmpErrMsg); if (Failed()) return
      CALL AllocAry(InitOutData%RotFrame_u, NumExtendedIO, 'RotFrame_u', TmpErrStat, TmpErrMsg); if (Failed()) return
      CALL AllocAry(InitOutData%IsLoad_u, NumExtendedIO, 'IsLoad_u', TmpErrStat, TmpErrMsg); if (Failed()) return
      CALL AllocAry(InitOutData%LinNames_y, NumExtendedIO + p%NumOuts, 'LinNames_y', TmpErrStat, TmpErrMsg); if (Failed()) return
      CALL AllocAry(InitOutData%RotFrame_y, NumExtendedIO + p%NumOuts, 'RotFrame_y', TmpErrStat, TmpErrMsg); if (Failed()) return
      
      ! Extended Inputs
      Lin_Indx = 0
      InitOutData%LinNames_u(Lin_Indx + 1) = 'Extended input: horizontal wind speed (steady/uniform wind), m/s'
      InitOutData%LinNames_u(Lin_Indx + 2) = 'Extended input: vertical power-law shear exponent, -'
      InitOutData%LinNames_u(Lin_Indx + 3) = 'Extended input: propagation direction, rad'

      ! Extended Outputs
      Lin_Indx = 0
      InitOutData%LinNames_y(Lin_Indx + 1) = 'Extended output: horizontal wind speed (steady/uniform wind), m/s'
      InitOutData%LinNames_y(Lin_Indx + 2) = 'Extended output: vertical power-law shear exponent, -'
      InitOutData%LinNames_y(Lin_Indx + 3) = 'Extended output: propagation direction, rad'

      ! Outputs
      do i=1,p%NumOuts
         InitOutData%LinNames_y(i+NumExtendedIO) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
      end do

      ! IfW inputs and outputs are in the global, not rotating frame
      InitOutData%RotFrame_u = .false.
      InitOutData%RotFrame_y = .false.
      InitOutData%IsLoad_u   = .false. ! IfW inputs for linearization are not loads

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
   
   
   logical function Failed()
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed
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
!! array is then rotated by p%PropagationDir so that it now corresponds to the global coordinate UVW values for wind
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
   CALL InflowWind_GetRotorSpatialAverage( Time, InputData, p, ContStates, DiscStates, ConstrStates, &
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
         
      CALL Lidar_CalcOutput(Time, InputData, p, OutputData, m, TmpErrStat, TmpErrMsg )
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

   ! Destroy all inflow wind derived types
   CALL InflowWind_DestroyInput( InputData, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyParam( p, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyContState( ContStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyDiscState( DiscStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyConstrState( ConstrStateGuess, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyOtherState( OtherStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyOutput( y, ErrStat, ErrMsg )                     
   CALL InflowWind_DestroyMisc( m, ErrStat, ErrMsg )                     

END SUBROUTINE InflowWind_End



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in IfW_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE InflowWind_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) 
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) 
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) 
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z)  
 
      ! local variables: 
   INTEGER(IntKi)                                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                                           :: ErrMsg2            ! temporary error message
   CHARACTER(*), PARAMETER                                        :: RoutineName = 'InflowWind_JacobianPInput'
   REAL(R8Ki)                                                     :: local_dYdu(3,NumExtendedIO)
   integer                                                        :: i,j, n
   integer                                                        :: i_start, i_end  ! indices for input/output start and end
   integer                                                        :: node, comp
      
      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdu ) ) THEN

      ! If dYdu is allocated, make sure it is the correct size
      if (allocated(dYdu)) then
         if (size(dYdu,1) /= NumExtendedIO + p%NumOuts)  deallocate (dYdu)
         if (size(dYdu,2) /= NumExtendedIO)              deallocate (dYdu)
      endif

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
      !  -  inputs are extended inputs only
      !  -  outputs are the extended outputs and the WriteOutput values
      if (.not. ALLOCATED(dYdu)) then
         CALL AllocAry( dYdu, NumExtendedIO + p%NumOuts, NumExtendedIO, 'dYdu', ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) return
      end if
         
         
      SELECT CASE ( p%FlowField%FieldType )
      CASE (Uniform_FieldType)
         dYdu = 0.0_R8Ki ! initialize all non-diagonal entries to zero (position of node effects the output of only that node) 

         ! Extended inputs to extended outputs (direct pass-through)
         do i=1,NumExtendedIO
            dYdu(i,i) = 1.0_R8Ki
         enddo

         ! WriteOutput velocities (note: may not have all of the components of each point) 
         do i=1, p%NumOuts
            node  = p%OutParamLinIndx(1,i) ! output node
            comp  = p%OutParamLinIndx(2,i) ! component of output node

            if (node > 0) then
               call IfW_UniformWind_JacobianPInput( p%FlowField%Uniform, t, p%WindViXYZ(:,node), p%FlowField%RotToWind(1,1), p%FlowField%RotToWind(2,1), local_dYdu )
            else
               local_dYdu = 0.0_R8Ki
               comp = 1
            end if
            dYdu(NumExtendedIO+i, 1:NumExtendedIO) = p%OutParam(i)%SignM * local_dYdu( comp , 1:NumExtendedIO)
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
   REAL(R8Ki),                INTENT(INOUT)  :: dYdu(3,NumExtendedIO) !< Partial derivatives of output functions (Y) with respect to the inputs (u)

   TYPE(UniformField_Interp)                 :: op                ! interpolated values of InterpParams
   REAL(R8Ki)                                :: RotatePosition(3) !< rotated position
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
   !! \frac{\partial Vt_u}{\partial V} & \frac{\partial Vt_u}{\partial VShr} & \frac{\partial Vt_u}{\partial PropDir} \\
   !! \frac{\partial Vt_v}{\partial V} & \frac{\partial Vt_v}{\partial VShr} & \frac{\partial Vt_v}{\partial PropDir} \\
   !! \frac{\partial Vt_w}{\partial V} & \frac{\partial Vt_w}{\partial VShr} & \frac{\partial Vt_w}{\partial PropDir} \\
   !! \end{bmatrix} \f$
   !-------------------------------------------------------------------------------------------------

   zr = RotatePosition(3)/UF%RefHeight

   dVhdV = ( ( RotatePosition(3)/UF%RefHeight ) ** op%ShrV &                                             ! power-law wind shear
             + ( op%ShrH   * ( RotatePosition(2) * op%CosAngleH + RotatePosition(1) * op%SinAngleH ) &   ! horizontal linear shear
             +  op%LinShrV * ( RotatePosition(3) - UF%RefHeight ) )/UF%RefLength  )                      ! vertical linear shear   
   Vh = op%VelH * dVhdV + op%VelGust

   dVhdVShr = op%VelH * zr**op%ShrV * log(zr)
   dVhdPD   = op%VelH * op%ShrH / UF%RefLength * ( RotatePosition(1) * op%CosAngleH - RotatePosition(2) * op%SinAngleH )

   tmp_du =  CosPropDir*op%CosAngleH  - SinPropDir*op%SinAngleH
   tmp_dv = -SinPropDir*op%CosAngleH  - CosPropDir*op%SinAngleH

      ! \f$ \frac{\partial Vt_u}{\partial V} =  \f$
   dYdu(1,1) = tmp_du*dVhdV      
      ! \f$ \frac{\partial Vt_v}{\partial V} =  \f$
   dYdu(2,1) = tmp_dv*dVhdV
      !> \f$ \frac{\partial Vt_w}{\partial V} = 0 \f$
   dYdu(3,1) = 0.0_R8Ki

      ! \f$ \frac{\partial Vt_u}{\partial VShr} =  \f$
   dYdu(1,2) = tmp_du*dVhdVShr
      ! \f$ \frac{\partial Vt_v}{\partial VShr} =  \f$
   dYdu(2,2) = tmp_dv*dVhdVShr
      !> \f$ \frac{\partial Vt_w}{\partial VShr} = 0 \f$
   dYdu(3,2) = 0.0_R8Ki

      ! \f$ \frac{\partial Vt_u}{\partial PropDir} =  \f$
   dYdu(1,3) = tmp_dv*Vh + tmp_du*dVhdPD
      ! \f$ \frac{\partial Vt_v}{\partial PropDir} =  \f$
   dYdu(2,3) = -tmp_du*Vh + tmp_dv*dVhdPD
      !> \f$ \frac{\partial Vt_w}{\partial PropDir} = 0 \f$
   dYdu(3,3) = 0.0_R8Ki

END SUBROUTINE IfW_UniformWind_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
!! Note: there are no states, so this routine is simply a placeholder to satisfy the framework and automate some glue code
SUBROUTINE InflowWind_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   return
!   IF ( PRESENT( dYdx ) ) THEN
!   END IF
!   IF ( PRESENT( dXdx ) ) THEN
!   END IF
!   IF ( PRESENT( dXddx ) ) THEN
!   END IF
!   IF ( PRESENT( dZdx ) ) THEN
!   END IF
END SUBROUTINE InflowWind_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
!! Note: there are no states, so this routine is simply a placeholder to satisfy the framework and automate some glue code
SUBROUTINE InflowWind_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddxd(:,:) !< Partial derivatives of discrete state
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   return

!   IF ( PRESENT( dYdxd ) ) THEN
!   END IF
!   IF ( PRESENT( dXdxd ) ) THEN
!   END IF
!   IF ( PRESENT( dXddxd ) ) THEN
!   END IF
!   IF ( PRESENT( dZdxd ) ) THEN
!   END IF
END SUBROUTINE InflowWind_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
!! Note: there are no states, so this routine is simply a placeholder to satisfy the framework and automate some glue code
SUBROUTINE InflowWind_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   return

!   IF ( PRESENT( dYdz ) ) THEN
!   END IF
!   IF ( PRESENT( dXdz ) ) THEN
!   END IF
!   IF ( PRESENT( dXddz ) ) THEN
!   END IF
!   IF ( PRESENT( dZdz ) ) THEN
!   END IF
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

   INTEGER(IntKi)                                    :: i
   real(ReKi)                                        :: tmp_op(NumExtendedIO)
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'InflowWind_GetOP'
   

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Since both u_op and y_op need this, calculate it up front
   if (present(u_op) .or. present(y_op)) then
      call IfW_UniformWind_GetOP( p%FlowField%Uniform, t, p%FlowField%VelInterpCubic, tmp_op )
      tmp_op(3) = p%FlowField%PropagationDir + tmp_op(3)  ! include the AngleH from Uniform Wind input files
   endif

   if ( PRESENT( u_op ) ) then
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, NumExtendedIO, 'u_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if

      u_op(1:NumExtendedIO) = tmp_op(1:NumExtendedIO)

   end if

   if ( PRESENT( y_op ) ) then
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, NumExtendedIO + p%NumOuts, 'y_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if

      y_op(1:NumExtendedIO) = tmp_op(1:NumExtendedIO)
      do i=1,p%NumOuts
         y_op(NumExtendedIO + i) = y%WriteOutput( i )
      end do
   end if

   return

!   IF ( PRESENT( x_op ) ) THEN
!   END IF
!   IF ( PRESENT( dx_op ) ) THEN
!   END IF
!   IF ( PRESENT( xd_op ) ) THEN
!   END IF
!   IF ( PRESENT( z_op ) ) THEN
!   END IF
END SUBROUTINE InflowWind_GetOP

END MODULE InflowWind
