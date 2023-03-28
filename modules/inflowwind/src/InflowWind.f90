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
module InflowWind

use NWTC_Library
use InflowWind_Types
use InflowWind_Subs
use InflowWind_IO_Types
use InflowWind_IO
use IfW_FlowField

use Lidar                      ! module for obtaining sensor data

implicit none
private

type(ProgDesc), parameter            :: IfW_Ver = ProgDesc('InflowWind', '', '')
integer, parameter            :: NumExtendedInputs = 3 !: V, VShr, PropDir

! ..... Public Subroutines ...................................................................................................

public :: InflowWind_Init                                   !< Initialization routine
public :: InflowWind_CalcOutput                             !< Calculate the wind velocities
public :: InflowWind_End                                    !< Ending routine (includes clean up)

! These routines satisfy the framework, but do nothing at present.
public :: InflowWind_UpdateStates               !< Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
public :: InflowWind_CalcConstrStateResidual    !< Tight coupling routine for returning the constraint state residual
public :: InflowWind_CalcContStateDeriv         !< Tight coupling routine for computing derivatives of continuous states
public :: InflowWind_UpdateDiscState            !< Tight coupling routine for updating discrete states

! These routines compute Jacobians; only dYdu is defined.
public :: InflowWind_JacobianPInput             !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
public :: InflowWind_JacobianPContState         !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                                   !!   states(x)
public :: InflowWind_JacobianPDiscState         !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                                   !!   states(xd)
public :: InflowWind_JacobianPConstrState       !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                                   !!   states(z)
public :: InflowWind_GetOP                      !< Routine to pack the operating point values (for linearization) into arrays

contains
!====================================================================================================
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! Since this module acts as an interface to other modules, on some things are set before initiating
!! calls to the lower modules.
!----------------------------------------------------------------------------------------------------
subroutine InflowWind_Init(InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                           y, m, TimeInterval, InitOutData, ErrStat, ErrMsg)

   ! Initialization data and guesses
   type(InflowWind_InitInputType), intent(IN)  :: InitInp           !< Input data for initialization
   type(InflowWind_InputType), intent(OUT)  :: InputGuess        !< An initial guess for the input; the input mesh must be defined
   type(InflowWind_ParameterType), intent(OUT)  :: p                 !< Parameters
   type(InflowWind_ContinuousStateType), intent(OUT)  :: ContStates        !< Initial continuous states
   type(InflowWind_DiscreteStateType), intent(OUT)  :: DiscStates        !< Initial discrete states
   type(InflowWind_ConstraintStateType), intent(OUT)  :: ConstrStateGuess  !< Initial guess of the constraint states
   type(InflowWind_OtherStateType), intent(OUT)  :: OtherStates       !< Initial other/optimization states
   type(InflowWind_OutputType), intent(OUT)  :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   type(InflowWind_MiscVarType), intent(OUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   real(DbKi), intent(IN)  :: TimeInterval      !< Coupling time interval in seconds: InflowWind does not change this.
   type(InflowWind_InitOutputType), intent(OUT)  :: InitOutData       !< Initial output data -- Names, units, and version info.
   integer(IntKi), intent(OUT)  :: ErrStat           !< Error status of the operation
   character(*), intent(OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   ! Local variables
   character(*), parameter                               :: RoutineName = "InflowWind_Init"

   type(InflowWind_InputFile)                            :: InputFileData        !< Data from input file

   type(Steady_InitInputType)                            :: Steady_InitInput
   type(Uniform_InitInputType)                           :: Uniform_InitInput
   type(TurbSim_InitInputType)                           :: TurbSim_InitInput
   type(Bladed_InitInputType)                            :: Bladed_InitInput
   type(Bladed_InitOutputType)                           :: Bladed_InitOutput
   type(HAWC_InitInputType)                              :: HAWC_InitInput
   type(User_InitInputType)                              :: User_InitInput
   type(Grid4D_InitInputType)                            :: Grid4D_InitInput
   type(Points_InitInputType)                            :: Points_InitInput

   type(WindFileDat)                                     :: FileDat

   ! TYPE(InflowWind_IO_InitInputType)                      :: FlowField_InitData     !< initialization info
   ! TYPE(InflowWind_IO_InitOutputType)                     :: FlowField_InitOutData  !< initialization output info

   type(FileInfoType)                                    :: InFileInfo    !< The derived type for holding the full input file for parsing -- we may pass this in the future
   character(1024)                                       :: PriPath

   integer(IntKi)                                        :: I, j              !< Generic counter
   integer(IntKi)                                        :: Lin_indx          !< Generic counter
   integer(IntKi)                                        :: SumFileUnit       !< Unit number for the summary file
   character(256)                                        :: SumFileName       !< Name of the summary file
   character(256)                                        :: EchoFileName      !< Name of the summary file
   character(1), parameter                               :: UVW(3) = (/'U', 'V', 'W'/)
   character(1), parameter                               :: XYZ(3) = (/'X', 'Y', 'Z'/)
   integer(IntKi)                                        :: TmpErrStat
   character(ErrMsgLen)                                  :: TmpErrMsg         !< temporary error message

   !----------------------------------------------------------------------------------------------
   ! Initialize variables and check to see if this module has been initialized before.
   !----------------------------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg = ""
   SumFileUnit = -1_IntKi ! set at beginning in case of error

   ! Set a few variables.

   p%DT = TimeInterval             ! InflowWind does not require a specific time interval, so this is never changed.
   call NWTC_Init()
   call DispNVD(IfW_Ver)

   !----------------------------------------------------------------------------------------------
   ! Read the input file
   !----------------------------------------------------------------------------------------------

   ! Set the names of the files based on the inputfilename
   p%RootFileName = InitInp%RootName
   if (LEN_TRIM(p%RootFileName) == 0) call GetRoot(InitInp%InputFileName, p%RootFileName)
   EchoFileName = TRIM(p%RootFileName)//".ech"
   SumFileName = TRIM(p%RootFileName)//".sum"

   ! these values (and others hard-coded in lidar_init) should be set in the input file, too
   InputFileData%SensorType = InitInp%lidar%SensorType
   InputFileData%NumPulseGate = InitInp%lidar%NumPulseGate
   InputFileData%RotorApexOffsetPos = InitInp%lidar%RotorApexOffsetPos
   InputFileData%LidRadialVel = InitInp%lidar%LidRadialVel

   ! Parse all the InflowWind related input files and populate the *_InitDataType derived types
   call GetPath(InitInp%InputFileName, PriPath)

   if (InitInp%UseInputFile) then
      call ProcessComFile(InitInp%InputFileName, InFileInfo, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   else
      call NWTC_Library_CopyFileInfoType(InitInp%PassedFileData, InFileInfo, MESH_NEWCOPY, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   end if

   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, InitInp%InputFileName, EchoFileName, &
                                      InitInp%FixedWindFileRootName, InitInp%TurbineID, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   ! If wind is Grid4D from FAST.Farm, set input file values
   if (InitInp%Use4Dext) then
      InputFileData%WindType = FDext_WindNumber
      InputFileData%PropagationDir = 0.0_ReKi ! wind is in XYZ coordinates (already rotated if necessary), so don't rotate it again
      InputFileData%VFlowAngle = 0.0_ReKi
      InputFileData%VelInterpCubic = .false.
   end if

   ! initialize sensor data
   call Lidar_Init(InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                   y, m, TimeInterval, InitOutData, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

   ! Validate the InflowWind input file information
   call InflowWind_ValidateInput(InitInp, InputFileData, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   ! If a summary file was requested, open it and write preliminary data.
   if (InputFileData%SumPrint) then
      call InflowWind_OpenSumFile(SumFileUnit, SumFileName, IfW_Ver, InputFileData%WindType, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   end if

   ! Allocate the array for passing points
   call AllocAry(InputGuess%PositionXYZ, 3, InitInp%NumWindPoints, "Array of positions at which to find wind velocities", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   InputGuess%PositionXYZ = 0.0_ReKi
   InputGuess%HubPosition = 0.0_ReKi
   call Eye(InputGuess%HubOrientation, TmpErrStat, TmpErrMsg)

   ! Allocate the array for passing velocities out
   call AllocAry(y%VelocityUVW, 3, InitInp%NumWindPoints, "Array of wind velocities returned by InflowWind", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   y%VelocityUVW = 0.0_ReKi

   ! If requested, allocate the array for passing accelerations out
   if (InitInp%OutputAccel) then
      call AllocAry(y%AccelUVW, 3, InitInp%NumWindPoints, "Array of wind accelerations returned by InflowWind", TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      y%AccelUVW = 0.0_ReKi
   end if

   !----------------------------------------------------------------------------
   ! Set flow field input data based on wind type
   !----------------------------------------------------------------------------

   InitOutData%WindFileInfo%MWS = HUGE(InitOutData%WindFileInfo%MWS)

   select case (InputFileData%WindType)

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
      end if

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
      end if

   case (TSFF_WindNumber)

      TurbSim_InitInput%WindFileName = InputFileData%TSFF_FileName

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_TurbSim_Init(TurbSim_InitInput, SumFileUnit, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

   case (BladedFF_WindNumber)

      Bladed_InitInput%TurbineID = InitInp%TurbineID
      if (InitInp%FixedWindFileRootName) then ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
         if (InitInp%TurbineID == 0) then     ! .TRUE. for the FAST.Farm low-resolution domain
            InputFileData%BladedFF_FileName = TRIM(InputFileData%BladedFF_FileName)//TRIM(PathSep)//'Low'
         else                                   ! FAST.Farm high-resolution domain(s)
            InputFileData%BladedFF_FileName = TRIM(InputFileData%BladedFF_FileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
         end if
      end if
      Bladed_InitInput%WindType = BladedFF_WindNumber
      Bladed_InitInput%WindFileName = TRIM(InputFileData%BladedFF_FileName)//'.wnd'
      Bladed_InitInput%TowerFileExist = InputFileData%BladedFF_TowerFile
      Bladed_InitInput%NativeBladedFmt = .false.

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_Bladed_Init(Bladed_InitInput, SumFileUnit, Bladed_InitOutput, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

   case (BladedFF_Shr_WindNumber)

      Bladed_InitInput%WindType = BladedFF_Shr_WindNumber
      Bladed_InitInput%TurbineID = InitInp%TurbineID
      Bladed_InitInput%WindFileName = InputFileData%BladedFF_FileName
      Bladed_InitInput%TowerFileExist = .false.
      Bladed_InitInput%NativeBladedFmt = .true.

      p%FlowField%FieldType = Grid3D_FieldType
      call IfW_Bladed_Init(Bladed_InitInput, SumFileUnit, Bladed_InitOutput, p%FlowField%Grid3D, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      ! Overwrite the values of PropagationDir and VFlowAngle with values from the native Bladed file
      InputFileData%PropagationDir = Bladed_InitOutput%PropagationDir
      InputFileData%VFlowAngle = Bladed_InitOutput%VFlowAngle

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
      end if

   case (User_WindNumber)

      p%FlowField%FieldType = User_FieldType
      call IfW_User_Init(User_InitInput, SumFileUnit, p%FlowField%User, FileDat, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

   case (FDext_WindNumber)

      p%FlowField%FieldType = Grid4D_FieldType
      call IfW_Grid4D_Init(InitInp%FDext, p%FlowField%Grid4D, FileDat, TmpErrStat, TmpErrMsg)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      ! case (Point_WindNumber)

      !    p%FlowField%FieldType = Point_FieldType

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

      ! Set box exceedence averaging grid
      p%FlowField%G3D%BoxExceedAllowF = InitInp%BoxExceedAllowF
      p%FlowField%G3D%BoxExceedAllowIdx = InitInp%BoxExceedAllowIdx
      if (p%FlowField%G3D%BoxExceedAllowF) then
         call GenMeanGridProfileTimeSeries(p%FlowField%G3D, TmpErrStat, TmpErrMsg)
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

   call InflowWind_SetParameters(InitInp, InputFileData, p, m, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   ! Allocate arrays for the WriteOutput
   call AllocAry(y%WriteOutput, p%NumOuts, 'WriteOutput', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   y%WriteOutput = 0.0_ReKi

   call AllocAry(InitOutData%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   call AllocAry(InitOutData%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

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
      call AllocAry(InitOutData%LinNames_u, InitInp%NumWindPoints*3 + size(InputGuess%HubPosition) + 3 + NumExtendedInputs, 'LinNames_u', TmpErrStat, TmpErrMsg) ! add hub position, orientation(3) + extended inputs
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      call AllocAry(InitOutData%RotFrame_u, InitInp%NumWindPoints*3 + size(InputGuess%HubPosition) + 3 + NumExtendedInputs, 'RotFrame_u', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      call AllocAry(InitOutData%IsLoad_u, InitInp%NumWindPoints*3 + size(InputGuess%HubPosition) + 3 + NumExtendedInputs, 'IsLoad_u', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      call AllocAry(InitOutData%LinNames_y, InitInp%NumWindPoints*3 + size(y%DiskVel) + p%NumOuts, 'LinNames_y', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      call AllocAry(InitOutData%RotFrame_y, InitInp%NumWindPoints*3 + size(y%DiskVel) + p%NumOuts, 'RotFrame_y', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      do i = 1, InitInp%NumWindPoints
         do j = 1, 3
            InitOutData%LinNames_y((i - 1)*3 + j) = UVW(j)//'-component inflow velocity at node '//trim(num2lstr(i))//', m/s'
            InitOutData%LinNames_u((i - 1)*3 + j) = XYZ(j)//'-component position of node '//trim(num2lstr(i))//', m'
         end do
      end do

      ! hub position
      Lin_Indx = InitInp%NumWindPoints*3
      do j = 1, 3
         InitOutData%LinNames_y(Lin_Indx + j) = 'average '//UVW(j)//'-component rotor-disk velocity, m/s'
         InitOutData%LinNames_u(Lin_Indx + j) = XYZ(j)//'-component position of moving hub, m'
      end do
      Lin_Indx = Lin_Indx + 3

      ! hub orientation angles
      do j = 1, 3
         InitOutData%LinNames_u(Lin_Indx + j) = XYZ(j)//' orientation of moving hub, rad'
      end do
      Lin_Indx = Lin_Indx + 3

      InitOutData%LinNames_u(Lin_Indx + 1) = 'Extended input: horizontal wind speed (steady/uniform wind), m/s'
      InitOutData%LinNames_u(Lin_Indx + 2) = 'Extended input: vertical power-law shear exponent, -'
      InitOutData%LinNames_u(Lin_Indx + 3) = 'Extended input: propagation direction, rad'

      do i = 1, p%NumOuts
         InitOutData%LinNames_y(i + 3*InitInp%NumWindPoints + size(y%DiskVel)) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
      end do

      ! IfW inputs and outputs are in the global, not rotating frame
      InitOutData%RotFrame_u = .false.
      InitOutData%RotFrame_y = .false.

      InitOutData%IsLoad_u = .false. ! IfW inputs for linearization are not loads

   end if

   ! Set the version information in InitOutData
   InitOutData%Ver = IfW_Ver

   call CleanUp()

contains

   subroutine CleanUp()

      ! add in stuff that we need to dispose of here
      call InflowWind_DestroyInputFile(InputFileData, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

      ! Ignore error messages from InFileInfo destruction
      call NWTC_Library_DestroyFileInfoType(InFileInfo, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

      ! Close the summary file if we were writing one
      if (SumFileUnit > 0) then
         call InflowWind_CloseSumFile(SumFileUnit, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      end if

   end subroutine CleanUp

end subroutine InflowWind_Init

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
subroutine InflowWind_CalcOutput(Time, InputData, p, &
                                 ContStates, DiscStates, ConstrStates, &   ! Framework required states -- empty in this case.
                                 OtherStates, OutputData, m, ErrStat, ErrMsg)

   character(*), parameter                     :: RoutineName = "InflowWind_CalcOutput"

   real(DbKi), intent(IN)  :: Time              !< Current simulation time in seconds
   type(InflowWind_InputType), intent(IN)  :: InputData         !< Inputs at Time
   type(InflowWind_ParameterType), intent(IN)  :: p                 !< Parameters
   type(InflowWind_ContinuousStateType), intent(IN)  :: ContStates        !< Continuous states at Time
   type(InflowWind_DiscreteStateType), intent(IN)  :: DiscStates        !< Discrete states at Time
   type(InflowWind_ConstraintStateType), intent(IN)  :: ConstrStates      !< Constraint states at Time
   type(InflowWind_OtherStateType), intent(IN)  :: OtherStates       !< Other/optimization states at Time
   type(InflowWind_OutputType), intent(INOUT)  :: OutputData        !< Outputs computed at Time (IN for mesh reasons and data allocation)
   type(InflowWind_MiscVarType), intent(INOUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   integer(IntKi), intent(OUT)  :: ErrStat           !< Error status of the operation
   character(*), intent(OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   ! Local variables
   integer(IntKi)                                           :: i
   ! Temporary variables for error handling
   integer(IntKi)                                           :: TmpErrStat
   character(ErrMsgLen)                                     :: TmpErrMsg            ! temporary error message

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Allocate the velocity array to get out
   if (.not. ALLOCATED(OutputData%VelocityUVW)) then
      call AllocAry(OutputData%VelocityUVW, 3, SIZE(InputData%PositionXYZ, DIM=2), &
                    "Velocity array returned from IfW_CalcOutput", TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   elseif (SIZE(InputData%PositionXYZ, DIM=2) /= SIZE(OutputData%VelocityUVW, DIM=2)) then
      call SetErrStat(ErrID_Fatal, " Programming error: Position and Velocity arrays are not sized the same.", &
                      ErrStat, ErrMsg, RoutineName)
   end if

   !-----------------------------
   ! Outputs: OutputData%VelocityUVW
   !-----------------------------
   call CalculateOutput(Time, InputData, p, ContStates, DiscStates, ConstrStates, &
                        OtherStates, OutputData, m, .true., TmpErrStat, TmpErrMsg)

   !-----------------------------
   ! Output: OutputData%DiskVel
   !-----------------------------
   call InflowWind_GetSpatialAverage(Time, InputData, p, ContStates, DiscStates, ConstrStates, &
                                     OtherStates, m, OutputData%DiskVel, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

   !-----------------------------
   ! Output: OutputData%HubVel
   !-----------------------------
   call InflowWind_GetHubValues(Time, InputData, p, ContStates, DiscStates, ConstrStates, &
                                OtherStates, m, OutputData%HubVel, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

   !-----------------------------
   ! Outputs: OutputData%lidar%LidSpeed and OutputData%lidar%WtTrunc
   !-----------------------------

   ! return sensor values
   if (p%lidar%SensorType /= SensorType_None) then

      call Lidar_CalcOutput(Time, InputData, p, &
                            ContStates, DiscStates, ConstrStates, OtherStates, &
                            OutputData, m, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

   end if

   !-----------------------------
   ! Outputs: OutputData%WriteOutput from m%AllOuts and OutputData%lidar%LidSpeed
   !-----------------------------

   call SetAllOuts(p, OutputData, m, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

   ! Map to the outputs
   do I = 1, p%NumOuts  ! Loop through all selected output channels
      OutputData%WriteOutput(I) = p%OutParam(I)%SignM*m%AllOuts(p%OutParam(I)%Indx)
   end do             ! I - All selected output channels

end subroutine InflowWind_CalcOutput

!====================================================================================================
!> Clean up the allocated variables and close all open files.  Reset the initialization flag so
!! that we have to reinitialize before calling the routines again.
subroutine InflowWind_End(InputData, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                          y, m, ErrStat, ErrMsg)

   type(InflowWind_InputType), intent(INOUT)  :: InputData         !< Input data for initialization
   type(InflowWind_ParameterType), intent(INOUT)  :: p                 !< Parameters
   type(InflowWind_ContinuousStateType), intent(INOUT)  :: ContStates        !< Continuous states
   type(InflowWind_DiscreteStateType), intent(INOUT)  :: DiscStates        !< Discrete states
   type(InflowWind_ConstraintStateType), intent(INOUT)  :: ConstrStateGuess  !< Guess of the constraint states
   type(InflowWind_OtherStateType), intent(INOUT)  :: OtherStates       !< Other/optimization states
   type(InflowWind_OutputType), intent(INOUT)  :: y                 !< Output data
   type(InflowWind_MiscVarType), intent(INOUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   integer(IntKi), intent(OUT)  :: ErrStat           !< error status
   character(*), intent(OUT)  :: ErrMsg            !< error message

   character(*), parameter                     :: RoutineName = "InflowWind_End"

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Reset the wind type so that the initialization routine must be called
   p%WindType = Undef_WindNumber

   ! Destroy all inflow wind derived types
   call InflowWind_DestroyInput(InputData, ErrStat, ErrMsg)
   call InflowWind_DestroyParam(p, ErrStat, ErrMsg, DeallocatePointers=.true.)
   call InflowWind_DestroyContState(ContStates, ErrStat, ErrMsg)
   call InflowWind_DestroyDiscState(DiscStates, ErrStat, ErrMsg)
   call InflowWind_DestroyConstrState(ConstrStateGuess, ErrStat, ErrMsg)
   call InflowWind_DestroyOtherState(OtherStates, ErrStat, ErrMsg)
   call InflowWind_DestroyOutput(y, ErrStat, ErrMsg)
   call InflowWind_DestroyMisc(m, ErrStat, ErrMsg)

end subroutine InflowWind_End

!====================================================================================================
! The following routines were added to satisfy the framework, but do nothing useful.
!====================================================================================================
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
subroutine InflowWind_UpdateStates(t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg)

   real(DbKi), intent(IN) :: t               !< Current simulation time in seconds
   integer(IntKi), intent(IN) :: n               !< Current step of the simulation: t = n*Interval
   type(InflowWind_InputType), intent(INOUT) :: Inputs(:)       !< Inputs at InputTimes (output only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi), intent(IN) :: InputTimes(:)   !< Times in seconds associated with Inputs
   type(InflowWind_ParameterType), intent(IN) :: p               !< Parameters
   type(InflowWind_ContinuousStateType), intent(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!    Output: Continuous states at t + Interval
   type(InflowWind_DiscreteStateType), intent(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!    Output: Discrete states at t  + Interval
   type(InflowWind_ConstraintStateType), intent(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
   type(InflowWind_OtherStateType), intent(INOUT) :: OtherState      !< Other states: Other states at t;
                                                                           !!   Output: Other states at t + Interval
   type(InflowWind_MiscVarType), intent(INOUT) :: m               !< Misc variables for optimization (not copied in glue code)
   integer(IntKi), intent(OUT) :: ErrStat         !< Error status of the operation
   character(*), intent(OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg = ""

   x%DummyContState = 0.0_ReKi
   xd%DummyDiscState = 0.0_ReKi
   z%DummyConstrState = 0.0_ReKi

   return

end subroutine InflowWind_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
subroutine InflowWind_CalcContStateDeriv(Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg)
!..................................................................................................................................

   real(DbKi), intent(IN)  :: Time        !< Current simulation time in seconds
   type(InflowWind_InputType), intent(IN)  :: u           !< Inputs at Time
   type(InflowWind_ParameterType), intent(IN)  :: p           !< Parameters
   type(InflowWind_ContinuousStateType), intent(IN)  :: x           !< Continuous states at Time
   type(InflowWind_DiscreteStateType), intent(IN)  :: xd          !< Discrete states at Time
   type(InflowWind_ConstraintStateType), intent(IN)  :: z           !< Constraint states at Time
   type(InflowWind_OtherStateType), intent(IN)  :: OtherState  !< Other states at Time
   type(InflowWind_MiscVarType), intent(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(InflowWind_ContinuousStateType), intent(OUT)  :: dxdt        !< Continuous state derivatives at Time
   integer(IntKi), intent(OUT)  :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Compute the first time derivatives of the continuous states here:

   dxdt%DummyContState = 0.0_ReKi

end subroutine InflowWind_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
subroutine InflowWind_UpdateDiscState(Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg)

   real(DbKi), intent(IN)  :: Time        !< Current simulation time in seconds
   type(InflowWind_InputType), intent(IN)  :: u           !< Inputs at Time
   type(InflowWind_ParameterType), intent(IN)  :: p           !< Parameters
   type(InflowWind_ContinuousStateType), intent(IN)  :: x           !< Continuous states at Time
   type(InflowWind_DiscreteStateType), intent(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                           !! Output: Discrete states at Time + Interval
   type(InflowWind_ConstraintStateType), intent(IN)  :: z           !< Constraint states at Time
   type(InflowWind_OtherStateType), intent(IN)  :: OtherState  !< Other states at Time
   type(InflowWind_MiscVarType), intent(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   integer(IntKi), intent(OUT)  :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Update discrete states here:

   ! StateData%DiscState =

end subroutine InflowWind_UpdateDiscState

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
subroutine InflowWind_CalcConstrStateResidual(Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg)

   real(DbKi), intent(IN)  :: Time        !< Current simulation time in seconds
   type(InflowWind_InputType), intent(IN)  :: u           !< Inputs at Time
   type(InflowWind_ParameterType), intent(IN)  :: p           !< Parameters
   type(InflowWind_ContinuousStateType), intent(IN)  :: x           !< Continuous states at Time
   type(InflowWind_DiscreteStateType), intent(IN)  :: xd          !< Discrete states at Time
   type(InflowWind_ConstraintStateType), intent(IN)  :: z           !< Constraint states at Time (possibly a guess)
   type(InflowWind_OtherStateType), intent(IN)  :: OtherState  !< Other states at Time
   type(InflowWind_MiscVarType), intent(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(InflowWind_ConstraintStateType), intent(OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                           !! the input values described above
   integer(IntKi), intent(OUT)  :: ErrStat     !< Error status of the operation
   character(*), intent(OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Solve for the constraint states here:

   z_residual%DummyConstrState = 0

end subroutine InflowWind_CalcConstrStateResidual

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in IfW_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
subroutine InflowWind_JacobianPInput(t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu)
!..................................................................................................................................

   real(DbKi), intent(IN)           :: t          !< Time in seconds at operating point
   type(InflowWind_InputType), intent(IN)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(InflowWind_ParameterType), intent(IN)           :: p          !< Parameters
   type(InflowWind_ContinuousStateType), intent(IN)           :: x          !< Continuous states at operating point
   type(InflowWind_DiscreteStateType), intent(IN)           :: xd         !< Discrete states at operating point
   type(InflowWind_ConstraintStateType), intent(IN)           :: z          !< Constraint states at operating point
   type(InflowWind_OtherStateType), intent(IN)           :: OtherState !< Other states at operating point
   type(InflowWind_OutputType), intent(IN)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                                !!   Output fields are not used by this routine, but type is
                                                                                !!   available here so that mesh parameter information (i.e.,
                                                                                !!   connectivity) does not have to be recalculated for dYdu.
   type(InflowWind_MiscVarType), intent(INOUT)           :: m          !< Misc/optimization variables
   integer(IntKi), intent(OUT)           :: ErrStat    !< Error status of the operation
   character(*), intent(OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional, intent(INOUT)           :: dYdu(:, :)  !< Partial derivatives of output functions (Y)
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional, intent(INOUT)           :: dXdu(:, :)  !< Partial derivatives of continuous state functions (X)
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional, intent(INOUT)           :: dXddu(:, :) !< Partial derivatives of discrete state functions (Xd)
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional, intent(INOUT)           :: dZdu(:, :)  !< Partial derivatives of constraint state functions (Z)
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]

   ! local variables:
   integer(IntKi)                                                 :: ErrStat2
   character(ErrMsgLen)                                           :: ErrMsg2            ! temporary error message
   character(*), parameter                                        :: RoutineName = 'InflowWind_JacobianPInput'

   real(R8Ki)                                                     :: local_dYdu(3, 3 + NumExtendedInputs)
   integer                                                        :: i, n
   integer                                                        :: i_start, i_end  ! indices for input/output start and end
   integer                                                        :: node, comp
   integer                                                        :: n_inputs
   integer                                                        :: n_outputs
   integer                                                        :: i_ExtendedInput_start
   integer                                                        :: i_WriteOutput

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg = ''

   if (PRESENT(dYdu)) then

      n_outputs = SIZE(u%PositionXYZ) + p%NumOuts + size(y%DiskVel)
      n_inputs = SIZE(u%PositionXYZ) + size(u%HubPosition) + 3 + NumExtendedInputs ! need to add 3 for u%HubOrientation
      i_ExtendedInput_start = n_inputs - NumExtendedInputs + 1 ! index for extended inputs starts 2 from end (encompasses 3 values: V, VShr, PropDir)
      i_WriteOutput = n_outputs - p%NumOuts ! index for where write outputs begin is i_WriteOutput + 1

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:

      ! outputs are all velocities at all positions plus the WriteOutput values
      !
      if (.not. ALLOCATED(dYdu)) then
         call AllocAry(dYdu, n_outputs, n_inputs, 'dYdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if

      select case (p%WindType)
      case (Steady_WindNumber, Uniform_WindNumber)

         ! note that we are including the propagation direction in the analytical derivative calculated
         ! inside IfW_UniformWind_JacobianPInput, so no need to transform input position vectors first

         dYdu = 0.0_R8Ki ! initialize all non-diagonal entries to zero (position of node effects the output of only that node)

         n = SIZE(u%PositionXYZ, 2)
         ! these are the positions used in the module coupling
         do i = 1, n
            ! note that p%FlowField%RotToWind(1,1) = cos(p%PropagationDir) and p%FlowField%RotToWind(2,1) = sin(p%PropagationDir), which are the
            ! values we need to compute the jacobian.
!!!FIX ME with the propagation values!!!!
            call IfW_UniformWind_JacobianPInput(p%FlowField%Uniform, t, u%PositionXYZ(:, i), p%FlowField%RotToWind(1, 1), p%FlowField%RotToWind(2, 1), local_dYdu)

            i_end = 3*i
            i_start = i_end - 2

            dYdu(i_start:i_end, i_start:i_end) = local_dYdu(:, 1:3)

            dYdu(i_start:i_end, i_ExtendedInput_start:) = local_dYdu(:, 4:6) ! extended inputs

         end do

         ! see InflowWind_GetSpatialAverage():

         ! location of y%DiskAvg
         i_start = 3*n + 1
         i_end = i_start + 2

         dYdu(i_start:i_end, :) = 0.0_R8Ki ! initialize because we're going to create averages

         do i = 1, IfW_NumPtsAvg
            m%u_Avg%PositionXYZ(:, i) = matmul(u%HubOrientation, p%PositionAvg(:, i)) + u%HubPosition
!!!FIX ME with the propagation values!!!!
            call IfW_UniformWind_JacobianPInput(p%FlowField%Uniform, t, m%u_Avg%PositionXYZ(:, i), p%FlowField%RotToWind(1, 1), p%FlowField%RotToWind(2, 1), local_dYdu)

            ! y%DiskAvg has the same index as u%HubPosition
            ! Also note that partial_(m%u_Avg%PositionXYZ) / partial_(u%HubPosition) is identity, so we can skip that part of the chain rule for these derivatives:
            dYdu(i_start:i_end, i_start:i_end) = dYdu(i_start:i_end, i_start:i_end) + local_dYdu(:, 1:3)
            dYdu(i_start:i_end, i_ExtendedInput_start:) = dYdu(i_start:i_end, i_ExtendedInput_start:) + local_dYdu(:, 4:6) ! extended inputs
         end do
         dYdu(i_start:i_end, i_start:i_end) = dYdu(i_start:i_end, i_start:i_end)/real(IfW_NumPtsAvg, R8Ki)
         dYdu(i_start:i_end, i_ExtendedInput_start:) = dYdu(i_start:i_end, i_ExtendedInput_start:)/real(IfW_NumPtsAvg, R8Ki)
!FIX ME:
         ! need to calculate dXYZdHubOrient = partial_(m%u_Avg%PositionXYZ) / partial_(u%HubOrientation)
         !dYdu(i_start:i_end,(i_start+3):(i_end+3)) = matmul( dYdu(i_start:i_end,i_start:i_end), dXYZdHubOrient )

         ! these are the InflowWind WriteOutput velocities (and note that we may not have all of the components of each point)
         ! they do not depend on the inputs, so the derivatives w.r.t. X, Y, Z are all zero
         do i = 1, p%NumOuts
            node = p%OutParamLinIndx(1, i) ! output node
            comp = p%OutParamLinIndx(2, i) ! component of output node

            if (node > 0) then
!!!FIX ME with the propagation values!!!!
               call IfW_UniformWind_JacobianPInput(p%FlowField%Uniform, t, p%WindViXYZ(:, node), p%FlowField%RotToWind(1, 1), p%FlowField%RotToWind(2, 1), local_dYdu)
            else
               local_dYdu = 0.0_R8Ki
               comp = 1
            end if

            dYdu(i_WriteOutput + i, i_ExtendedInput_start:) = p%OutParam(i)%SignM*local_dYdu(comp, 4:6)
         end do

      case DEFAULT

      end select

   end if

   if (PRESENT(dXdu)) then
      if (allocated(dXdu)) deallocate (dXdu)
   end if

   if (PRESENT(dXddu)) then
      if (allocated(dXddu)) deallocate (dXddu)
   end if

   if (PRESENT(dZdu)) then
      if (allocated(dZdu)) deallocate (dZdu)
   end if

end subroutine InflowWind_JacobianPInput
!..................................................................................................................................
!> Routine to compute the Jacobians of the output (Y) function with respect to the inputs (u). The partial
!! derivative dY/du is returned. This submodule does not follow the modularization framework.
subroutine IfW_UniformWind_JacobianPInput(UF, t, Position, CosPropDir, SinPropDir, dYdu)
   use IfW_FlowField, only: UniformField_InterpLinear, UniformField_InterpCubic

   type(UniformFieldType), intent(IN)  :: UF                !< Uniform field derived type
   real(DbKi), intent(IN)  :: t                 !< Current simulation time in seconds
   real(ReKi), intent(IN)  :: Position(3)       !< XYZ Position at which to find velocity (operating point)
   real(ReKi), intent(IN)  :: CosPropDir        !< cosine of InflowWind propagation direction
   real(ReKi), intent(IN)  :: SinPropDir        !< sine of InflowWind propagation direction
   real(R8Ki), intent(INOUT)  :: dYdu(3, 6)         !< Partial derivatives of output functions (Y) with respect to the inputs (u)

   type(UniformField_Interp)                 :: op                ! interpolated values of InterpParams
   real(R8Ki)                                :: RotatePosition(3) !< rotated position
   real(R8Ki)                                :: dVhdx             ! temporary value to hold partial v_h partial X
   real(R8Ki)                                :: dVhdy             ! temporary value to hold partial v_h partial Y
   real(R8Ki)                                :: dVhdz             ! temporary value to hold partial v_h partial Z
   real(R8Ki)                                :: tmp_du            ! temporary value to hold calculations that are part of multiple components
   real(R8Ki)                                :: tmp_dv            ! temporary value to hold calculations that are part of multiple components
   real(R8Ki)                                :: dVhdPD            ! temporary value to hold partial v_h partial propagation direction
   real(R8Ki)                                :: dVhdV             ! temporary value to hold partial v_h partial V
   real(R8Ki)                                :: Vh                ! temporary value to hold v_h
   real(R8Ki)                                :: dVhdVShr          ! temporary value to hold partial v_h partial VShr
   real(R8Ki)                                :: zr

   if (Position(3) < 0.0_ReKi .or. EqualRealNos(Position(3), 0.0_ReKi)) then
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
   !> 2. Calculate $ rac{\partial Y_{Output \, Equations}}{\partial u_{inputs}} = egin{bmatrix}
   !! rac{\partial Vt_u}{\partial X} & rac{\partial Vt_u}{\partial Y} & rac{\partial Vt_u}{\partial Z}    !! rac{\partial Vt_v}{\partial X} & rac{\partial Vt_v}{\partial Y} & rac{\partial Vt_v}{\partial Z}    !! rac{\partial Vt_w}{\partial X} & rac{\partial Vt_w}{\partial Y} & rac{\partial Vt_w}{\partial Z}    !! nd{bmatrix} $
   !-------------------------------------------------------------------------------------------------

   zr = RotatePosition(3)/UF%RefHeight
   tmp_du = op%VelH*op%ShrH/UF%RefLength*CosPropDir
   dVhdx = tmp_du*op%SinAngleH
   dVhdy = tmp_du*op%CosAngleH
   dVhdz = op%VelH*(op%ShrV/UF%RefHeight*zr**(op%ShrV - 1.0_R8Ki) + op%LinShrV/UF%RefLength)

   dVhdV = ((RotatePosition(3)/UF%RefHeight)**op%ShrV &                                             ! power-law wind shear
            + (op%ShrH*(RotatePosition(2)*op%CosAngleH + RotatePosition(1)*op%SinAngleH) &   ! horizontal linear shear
               + op%LinShrV*(RotatePosition(3) - UF%RefHeight))/UF%RefLength)                      ! vertical linear shear
   Vh = op%VelH*dVhdV + op%VelGust

   dVhdVShr = op%VelH*zr**op%ShrV*log(zr)
   dVhdPD = op%VelH*op%ShrH/UF%RefLength*(RotatePosition(1)*op%CosAngleH - RotatePosition(2)*op%SinAngleH)

   tmp_du = CosPropDir*op%CosAngleH - SinPropDir*op%SinAngleH
   tmp_dv = -SinPropDir*op%CosAngleH - CosPropDir*op%SinAngleH

   !> $ rac{\partial Vt_u}{\partial X} = \left[
