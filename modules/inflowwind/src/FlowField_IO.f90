!>  This module uses grid-field binary wind files to determine the wind inflow.
!!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!!  and that all units are specified in the metric system (using meters and seconds).
!!  Data is shifted by half the grid width to account for turbine yaw (so that data in the X
!!  direction actually starts at -1*p%FFYHWid meters).
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2022  National Renewable Energy Laboratory
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

module FlowField_IO

use NWTC_Library
use FlowField_IO_Types
use FlowField

implicit none
private
public :: FlowField_IO_Init

type(ProgDesc), parameter :: FlowField_IO_Ver = ProgDesc('FlowField_IO', '', '')

integer(IntKi), parameter :: ScaleMethod_None = 0, &           !< no scaling
                             ScaleMethod_Direct = 1, &         !< direct scaling factors
                             ScaleMethod_StdDev = 2            !< requested standard deviation

type :: TurbSimHeaderType
   integer(B2Ki)  :: FileID
   integer(B4Ki)  :: NZGrids, NYGrids, NTGrids, NSteps
   real(SiKi)     :: dz, dy, dt
   real(SiKi)     :: mws, ref_height, grid_base_height
   real(SiKi)     :: VslopeX, VoffsetX
   real(SiKi)     :: VslopeY, VoffsetY
   real(SiKi)     :: VslopeZ, VoffsetZ
   integer(B4Ki)  :: DescLen
end type

contains

!----------------------------------------------------------------------------------------------------
!> A subroutine to initialize the UserWind module. This routine will initialize the module.
!----------------------------------------------------------------------------------------------------
subroutine FlowField_IO_Init(InitInp, FF, InitOut, ErrStat, ErrMsg)

   type(FlowField_IO_InitInputType), intent(in)   :: InitInp           !< Input data for initialization
   type(FlowFieldType), intent(out)                :: FF                !< Flow field
   type(FlowField_IO_InitOutputType), intent(out) :: InitOut           !< Misc variables for optimization (not copied in glue code)
   integer(IntKi), intent(out)                     :: ErrStat           !< determines if an error has been encountered
   character(*), intent(out)                       :: ErrMsg            !< A message about the error.  See NWTC_Library info for ErrID_* levels.

   character(*), parameter                         :: RoutineName = 'FlowField_IO_Init'
   integer(IntKi)                                  :: TmpErrStat        ! Temp variable for the error status
   character(ErrMsgLen)                            :: TmpErrMsg         ! temporary error message

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Parameter initialization
   !----------------------------------------------------------------------------

   FF%PropagationDir = InitInp%PropagationDir*D2R
   FF%VFlowAngle = InitInp%VFlowAngle*D2R

   ! Shift propagation direction so it is between -pi and pi
   call MPi2Pi(FF%PropagationDir)

   !----------------------------------------------------------------------------
   ! Wind Type Initialization
   !----------------------------------------------------------------------------

   select case (InitInp%WindType)
   case (1) ! Steady
      FF%FieldType = Uniform_FieldType
      call SteadyWind_Init(InitInp%Steady, InitInp%SumFileUnit, FF%Uniform, InitOut%FileDat, TmpErrStat, TmpErrMsg)

   case (2) ! Uniform
      FF%FieldType = Uniform_FieldType
      call UniformWind_Init(InitInp%Uniform, InitInp%SumFileUnit, FF%Uniform, InitOut%FileDat, TmpErrStat, TmpErrMsg)

   case (3) ! Binary TurbSim FF
      FF%FieldType = Grid3D_FieldType
      call TurbSim_Init(InitInp%TurbSim, InitInp%SumFileUnit, FF%Grid3D, InitOut%FileDat, TmpErrStat, TmpErrMsg)

   case (4) ! Binary Bladed-Style FF
      FF%FieldType = Grid3D_FieldType
      ! call Read_Bladed_Binary(p, InitInp, InitOut, WindFileUnit, SumFileUnit, ErrStat, ErrMsg)
      TmpErrStat = ErrID_Fatal
      TmpErrMsg = "Binary Bladed Wind is not implemented"

   case (5) ! HAWC
      FF%FieldType = Grid3D_FieldType
      call HAWC_Init(InitInp%HAWC, InitInp%SumFileUnit, FF%Grid3D, InitOut%FileDat, TmpErrStat, TmpErrMsg)

   case (6) ! User Defined
      FF%FieldType = Grid3D_FieldType
      TmpErrStat = ErrID_Fatal
      TmpErrMsg = "User Wind is not implemented"

   case (7) ! Native Bladed FF
      FF%FieldType = Grid3D_FieldType
      ! call Read_Bladed_Native(p, InitInp, InitOut, WindFileUnit, SumFileUnit, ErrStat, ErrMsg)
      TmpErrStat = ErrID_Fatal
      TmpErrMsg = "Native Bladed Wind is not implemented"

   case (8) ! External Grid
      FF%FieldType = Grid4D_FieldType
      TmpErrStat = ErrID_Fatal
      TmpErrMsg = "Grid4D Wind is not implemented"

   case (9) ! External Point
      FF%FieldType = Point_FieldType
      TmpErrStat = ErrID_Fatal
      TmpErrMsg = "Points Wind is not implemented"

   case default ! Others
      TmpErrStat = ErrID_Fatal
      TmpErrMsg = "Unknown wind type"

   end select
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Field Type Initialization
   !----------------------------------------------------------------------------

   ! Reset flag indicating that acceleration field is valid
   FF%AccFieldValid = .false.

   select case (FF%FieldType)
   case (Uniform_FieldType)

      if (InitInp%OutputAccel .or. (InitInp%VelInterpCubic)) then
         call UniformField_CalcAccel(FF%Uniform, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
         FF%AccFieldValid = .true.
      end if
      FF%VelInterpCubic = InitInp%VelInterpCubic

   case (Grid3D_FieldType)
      if (InitInp%OutputAccel .or. (InitInp%VelInterpCubic)) then
         call Grid3DField_CalcAccel(FF%Grid3D, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
         FF%AccFieldValid = .true.
      end if
      FF%VelInterpCubic = InitInp%VelInterpCubic

   case default
      if (InitInp%OutputAccel) then
         call SetErrStat(ErrID_Fatal, "Acceleration not implemented for field type", &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if
      if (InitInp%VelInterpCubic) then
         CALL WrScr ( ' Cubic velocity interpolation not implemented for WindType '//&
                       num2LStr(InitInp%WindType)//', using linear interpolation' )
         FF%VelInterpCubic = .false.
      end if
   end select

   !----------------------------------------------------------------------------
   ! Setup the coordinate transforms for rotating the wind field.
   !----------------------------------------------------------------------------

   ! Create the rotation matrices
   ! rotate from XYZ to X'Y'Z' (wind aligned along X) coordinates
   ! Includes the wind upflow (inclination) angle (rotation about Y axis)
   FF%RotToWind(1, :) = [cos(-FF%VFlowAngle)*cos(-FF%PropagationDir), &
                         cos(-FF%VFlowAngle)*sin(-FF%PropagationDir), &
                         -sin(-FF%VFlowAngle)]
   FF%RotToWind(2, :) = [-sin(-FF%PropagationDir), &
                         cos(-FF%PropagationDir), &
                         0.0_ReKi]
   FF%RotToWind(3, :) = [sin(-FF%VFlowAngle)*cos(-FF%PropagationDir), &
                         sin(-FF%VFlowAngle)*sin(-FF%PropagationDir), &
                         cos(-FF%VFlowAngle)]

   ! Create the rotation matrices -- rotate from X'Y'Z' (wind aligned along X)
   ! to global XYZ coordinates: this is the same as a rotation about the
   ! (positive) upflow angle multiplied by a rotation about the (positive) wind direction:
   ! Global wind = R(p%PropagationDir) * R(p%VFlowAngle) * [local wind]
   ! local wind = R( -p%VFlowAngle) * R (-p%PropagationDir) [global wind]
   !            = R^T(p%VFlowAngle) * R^T(p%PropagationDir) [global wind]
   !            = (R(p%PropagationDir) * R(p%VFlowAngle))^T [global wind]
   FF%RotFromWind = transpose(FF%RotToWind)

   FF%RotateWindBox = .not. (EqualRealNos(FF%PropagationDir, 0.0_ReKi) .and. &
                             EqualRealNos(FF%VFlowAngle, 0.0_ReKi))

   FF%RefPosition = [0.0_ReKi, 0.0_ReKi, InitOut%FileDat%RefHt]

   !----------------------------------------------------------------------------
   ! Initialization Output
   !----------------------------------------------------------------------------

   InitOut%Ver = FlowField_IO_Ver

end subroutine

subroutine SteadyWind_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
   type(SteadyInitInputType), intent(in)  :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(UniformFieldType), intent(out)    :: UF
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'SteadyWind_Init'
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Set parameters from inititialization input
   UF%DataSize = 1
   UF%RefHeight = InitInp%RefHt
   UF%RefLength = 1.0_ReKi

   ! Allocate uniform wind data arrays
   call Alloc_UniformWindArrays(UF, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Set data values
   UF%Time = 0.0_ReKi
   UF%VelH = InitInp%HWindSpeed
   UF%VelV = 0.0_ReKi
   UF%VelGust = 0.0_ReKi
   UF%AngleH = 0.0_ReKi
   UF%AngleV = 0.0_ReKi
   UF%ShrH = 0.0_ReKi
   UF%ShrV = InitInp%PLExp
   UF%LinShrV = 0.0_ReKi

   !----------------------------------------------------------------------------
   ! Set wind file data
   !----------------------------------------------------------------------------

   FileDat%FileName = ""
   FileDat%WindType = 1
   FileDat%RefHt = UF%RefHeight
   FileDat%RefHt_Set = .false.
   FileDat%DT = 0.0_ReKi
   FileDat%NumTSteps = UF%DataSize
   FileDat%ConstantDT = .false.
   FileDat%TRange = 0.0_ReKi
   FileDat%TRange_Limited = .false.
   FileDat%YRange = 0.0_ReKi
   FileDat%YRange_Limited = .false.
   FileDat%ZRange = 0.0_ReKi
   FileDat%ZRange_Limited = .false.
   FileDat%BinaryFormat = 0_IntKi
   FileDat%IsBinary = .false.
   FileDat%TI = 0.0_ReKi
   FileDat%TI_listed = .false.
   FileDat%MWS = InitInp%HWindSpeed

   !----------------------------------------------------------------------------
   ! Write summary file if applicable
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A80)') 'Steady wind -- Constant wind profile for entire simulation. No windfile read in.'
      write (SumFileUnit, '(A40,G12.4)') '     Reference height:                  ', UF%RefHeight
      write (SumFileUnit, '(A40,G12.4)') '     Horizontal velocity:               ', UF%VelH(1)
      write (SumFileUnit, '(A40,G12.4)') '     Vertical sheer power law exponent: ', UF%ShrV(1)

      ! Get IO status for unit
      inquire (SumFileUnit, iostat=TmpErrStat)
      if (TmpErrStat /= 0_IntKi) then
         call SetErrStat(ErrID_Fatal, 'Error writing to summary file.', ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

end subroutine

subroutine UniformWind_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
   type(UniformInitInputType), intent(in)    :: InitInp
   integer(IntKi), intent(in)                :: SumFileUnit
   type(UniformFieldType), intent(out)       :: UF
   type(WindFileDat), intent(out)            :: FileDat
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter                   :: RoutineName = 'UniformWind_Init'
   integer(IntKi), parameter                 :: MaxNumCols = 9
   integer(IntKi), parameter                 :: MaxTries = 100
   integer(IntKi)                            :: NumCols
   integer(IntKi)                            :: LineNo, ILine
   integer(IntKi)                            :: i
   type(FileInfoType)                        :: WindFileInfo
   logical                                   :: WindFileConstantDT
   real(DbKi)                                :: WindFileDT
   real(ReKi)                                :: DirDiff
   real(ReKi)                                :: TmpData(MaxNumCols)  ! Temp variable for storing wind file data
   integer(IntKi)                            :: TmpErrStat
   character(ErrMsgLen)                      :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Set parameters from inititialization input
   UF%RefHeight = InitInp%RefHt
   UF%RefLength = InitInp%RefLength

   ! Read wind data from file or init input data
   if (InitInp%UseInputFile) then
      call ProcessComFile(InitInp%WindFileName, WindFileInfo, TmpErrStat, TmpErrMsg)
   else
      call NWTC_Library_CopyFileInfoType(InitInp%PassedFileData, WindFileInfo, MESH_NEWCOPY, TmpErrStat, TmpErrMsg)
   end if
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Get number of data lines in file
   UF%DataSize = WindFileInfo%NumLines

   ! Allocate uniform wind data arrays
   call Alloc_UniformWindArrays(UF, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Parse wind file data
   !----------------------------------------------------------------------------

   ! Initialize reading variables
   NumCols = MaxNumCols
   LineNo = 1

   ! Attempt to read 9 columns from line; if error, upflow column is not in file
   ! so set upflow=0 and reduce number of columns to read remaning data
   call ParseAry(WindFileInfo, LineNo, "HH file line", TmpData(1:NumCols), NumCols, TmpErrStat, TmpErrMsg)
   if (TmpErrStat /= 0) then
      CALL SetErrStat(ErrID_Info,' Could not read upflow column in uniform wind files. Assuming upflow is 0.', ErrStat, ErrMsg, RoutineName)
      UF%AngleV = 0
      NumCols = MaxNumCols - 1
   end if

   ! Parse the data and store it
   LineNo = 1
   do i = 1, UF%DataSize
      call ParseAry(WindFileInfo, LineNo, "HH file line", TmpData(1:NumCols), NumCols, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      UF%Time(i) = TmpData(1)
      UF%VelH(i) = TmpData(2)
      UF%AngleH(i) = TmpData(3)*D2R
      UF%VelV(i) = TmpData(4)
      UF%ShrH(i) = TmpData(5)
      UF%ShrV(i) = TmpData(6)
      UF%LinShrV(i) = TmpData(7)
      UF%VelGust(i) = TmpData(8)
      if (NumCols > 8) UF%AngleV(i) = TmpData(9)*D2R
   end do

   !----------------------------------------------------------------------------
   ! Ensure the wind direction isn't jumping more than 180 degrees between
   ! any 2 consecutive input times.
   ! (Avoids interpolation errors with modular arithemetic.)
   !----------------------------------------------------------------------------

   do i = 2, UF%DataSize
      ILine = 1
      do while (ILine < MaxTries)
         DirDiff = (UF%AngleH(i) - UF%AngleH(i - 1))
         if (abs(DirDiff) < Pi) exit
         UF%AngleH(i) = UF%AngleH(i) - SIGN(TwoPi, DirDiff)
         ILine = ILine + 1
      end do

      if (ILine >= MaxTries) then
         TmpErrMsg = ' Error calculating wind direction from uniform wind file. p%AngleH(' &
                     //TRIM(Num2LStr(i))//') = '//TRIM(Num2LStr(UF%AngleH(i)))//'; p%AngleH(' &
                     //TRIM(Num2LStr(i + 1))//') = '//TRIM(Num2LStr(UF%AngleH(i + 1)))
         call SetErrStat(ErrID_Fatal, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      end if
   end do

   !----------------------------------------------------------------------------
   ! Find out information on the timesteps and range
   !----------------------------------------------------------------------------

   ! Uniform timesteps
   if (UF%DataSize > 3) then
      WindFileConstantDT = .true.
      WindFileDT = UF%Time(2) - UF%Time(1)
      do I = 3, UF%DataSize
         if (.not. EqualRealNos((UF%Time(i) - UF%Time(i - 1)), real(WindFileDT, ReKi))) then
            WindFileConstantDT = .false.
            exit
         end if
      end do
   else
      ! Insufficient points to check, report that the timesteps are not uniform
      WindFileConstantDT = .false.
      WindFileDT = 0.0_ReKi
   end if

   !----------------------------------------------------------------------------
   ! Store wind file metadata
   !----------------------------------------------------------------------------

   FileDat%FileName = InitInp%WindFileName
   FileDat%WindType = 2
   FileDat%RefHt = UF%RefHeight
   FileDat%RefHt_Set = .false.
   FileDat%DT = WindFileDT
   FileDat%NumTSteps = UF%DataSize
   FileDat%ConstantDT = WindFileConstantDT
   FileDat%TRange = [UF%Time(1), UF%Time(UF%DataSize)]
   FileDat%TRange_Limited = .false.
   FileDat%YRange = (/0.0_ReKi, 0.0_ReKi/)
   FileDat%YRange_Limited = .false.
   FileDat%ZRange = (/0.0_ReKi, 0.0_ReKi/)
   FileDat%ZRange_Limited = .false.
   FileDat%BinaryFormat = 0_IntKi
   FileDat%IsBinary = .false.
   FileDat%TI = 0.0_ReKi
   FileDat%TI_listed = .false.

   if (UF%DataSize == 1) then
      FileDat%MWS = UF%VelH(1)
   else
      FileDat%MWS = 0.0_ReKi
      do i = 2, UF%DataSize
         FileDat%MWS = FileDat%MWS + 0.5_ReKi*(UF%VelH(i) + UF%VelH(i - 1))* &
                       (UF%Time(i) - UF%Time(i - 1))
      end do
      FileDat%MWS = FileDat%MWS/(UF%Time(UF%DataSize) - UF%Time(1))
   end if

   ! Check if the fist data point from the file is not along the X-axis while applying the windfield rotation
   if ((.not. EqualRealNos(UF%AngleH(1), 0.0_ReKi)) .and. &
       (.not. EqualRealNos(InitInp%PropagationDir, 0.0_ReKi))) then
      call SetErrStat(ErrID_Warn, ' Possible double rotation of wind field! Uniform wind file starts with a wind direction of '// &
                      TRIM(Num2LStr(UF%AngleH(1)*R2D))// &
                      ' degrees and the InflowWind input file specifies a PropagationDir of '// &
                      TRIM(Num2LStr(InitInp%PropagationDir*R2D))//' degrees.', &
                      ErrStat, ErrMsg, RoutineName)
   end if

   !----------------------------------------------------------------------------
   ! Print warnings and messages
   !----------------------------------------------------------------------------

   if (UF%Time(1) > 0.0) then
      TmpErrMsg = 'The uniform wind file : "'//TRIM(ADJUSTL(InitInp%WindFileName))// &
                  '" starts at a time '//'greater than zero. Interpolation errors may result.'
      call SetErrStat(ErrID_Warn, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   end if

   if (UF%DataSize == 1) then
      TmpErrMsg = ' Only 1 line in uniform wind file. Steady, horizontal wind speed at the hub height is '// &
                  TRIM(Num2LStr(UF%VelH(1)))//' m/s.'
      call SetErrStat(ErrID_Info, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   end if

   !----------------------------------------------------------------------------
   ! Write to the summary file
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A)') 'Uniform wind.  Module '//TRIM(FlowField_IO_Ver%Name)//' '//TRIM(FlowField_IO_Ver%Ver)
      write (SumFileUnit, '(A)') '     FileName:                    '//TRIM(InitInp%WindFileName)
      write (SumFileUnit, '(A34,G12.4)') '     Reference height (m):        ', UF%RefHeight
      write (SumFileUnit, '(A34,G12.4)') '     Reference length (m):        ', UF%RefLength
      write (SumFileUnit, '(A32,I8)') '     Number of data lines:        ', UF%DataSize
      write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
         TRIM(Num2LStr(UF%Time(1)))//' : '//TRIM(Num2LStr(UF%Time(UF%DataSize)))//' ]'

      ! Get IO status for unit
      inquire (SumFileUnit, iostat=TmpErrStat)
      if (TmpErrStat /= 0_IntKi) then
         call SetErrStat(ErrID_Fatal, 'Error writing to summary file.', ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

end subroutine

subroutine Alloc_UniformWindArrays(UF, ErrStat, ErrMsg)
   type(UniformFieldType), intent(inout)     :: UF
   integer(IntKi), intent(OUT)               :: ErrStat
   character(*), intent(OUT)                 :: ErrMsg

   character(*), parameter                   :: RoutineName = 'Alloc_UniformWindArrays'
   integer(IntKi)                            :: TmpErrStat
   character(ErrMsgLen)                      :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   if (.not. allocated(UF%Time)) then
      call AllocAry(UF%Time, UF%DataSize, 'Uniform wind time', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%VelH)) then
      call AllocAry(UF%VelH, UF%DataSize, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%AngleH)) then
      call AllocAry(UF%AngleH, UF%DataSize, 'Uniform wind direction', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%AngleV)) then
      call AllocAry(UF%AngleV, UF%DataSize, 'Uniform wind upflow angle', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%VelV)) then
      call AllocAry(UF%VelV, UF%DataSize, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%ShrH)) then
      call AllocAry(UF%ShrH, UF%DataSize, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%ShrV)) then
      call AllocAry(UF%ShrV, UF%DataSize, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%LinShrV)) then
      call AllocAry(UF%LinShrV, UF%DataSize, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

   if (.not. allocated(UF%VelGust)) then
      call AllocAry(UF%VelGust, UF%DataSize, 'Uniform gust velocity', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return
   end if

end subroutine

!> Read_TurbSim reads the binary TurbSim-format FF file (.bts).  It fills the FFData array with
!! velocity data for the grids and fills the Tower array with velocities at points on the tower
!! (if data exists).
subroutine TurbSim_Init(InitInp, SumFileUnit, GF, FileDat, ErrStat, ErrMsg)

   type(TurbSimInitInputType), intent(in)    :: InitInp
   integer(IntKi), intent(in)                :: SumFileUnit
   type(Grid3DFieldType), intent(out)          :: GF
   type(WindFileDat), intent(out)            :: FileDat
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter       :: RoutineName = "TurbSim_Init"
   integer(IntKi)                :: WindFileUnit
   integer(B2Ki), allocatable    :: VelRaw(:, :, :)   ! raw grid-field velocity data at one time step
   integer(B2Ki), allocatable    :: TwrRaw(:, :)      ! raw towrer velocity data at one time step
   character(:), allocatable     :: DescStr           ! description string contained in the file
   integer(IntKi)                :: IC                ! loop counter for wind components
   integer(IntKi)                :: IT                ! loop counter for time
   integer(IntKi)                :: NChar             ! number of characters in the description string
   real(SiKi)                    :: Vslope(3)         ! slope  for "un-normalizing" data
   real(SiKi)                    :: Voffset(3)        ! offset for "un-normalizing" data
   integer(IntKi)                :: TmpErrStat        ! temporary error status
   character(ErrMsgLen)          :: TmpErrMsg         ! temporary error message

   type(TurbSimHeaderType) :: header

   !----------------------------------------------------------------------------
   ! Initialize error variables
   !----------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Open the binary wind file and read header
   !----------------------------------------------------------------------------

   ! Get a unit number to use for the wind file
   call GetNewUnit(WindFileUnit, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Open binary file
   call OpenBInpFile(WindFileUnit, TRIM(InitInp%WindFileName), TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Read header from file
   read (WindFileUnit, IOSTAT=TmpErrStat) header
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading header of the FF binary file "' &
                      //TRIM(InitInp%WindFileName)//'."', ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Populate parameter data from header
   !----------------------------------------------------------------------------

   GF%WindFileFormat = header%FileID    ! file format identifier
   GF%Periodic = header%FileID == 8     ! 7 is used for non-periodic wind files; 8 is periodic wind
   GF%InterpTower = .false.             ! wind should not be interpolated at tower
   GF%AddMeanAfterInterp = .false.      ! do not add mean wind speed after interpolation

   GF%DTime = real(header%dt, ReKi)     ! grid spacing in time (dt), m/s
   GF%Rate = 1.0_ReKi/GF%DTime           ! Data rate (1/DTime), Hertz

   GF%NComp = 3                         ! TurbSim file file contains 3 wind components
   GF%NYGrids = header%NYGrids          ! the number of grid points laterally
   GF%NZGrids = header%NZGrids          ! the number of grid points vertically
   GF%NTGrids = header%NTGrids          ! the number of tower points
   GF%NSteps = header%NSteps          ! the number of time steps

   GF%InvDY = 1.0_ReKi/real(header%dy, ReKi)     ! 1/dy
   GF%YHWid = 0.5_ReKi*(GF%NYGrids - 1)/GF%InvDY   ! half the grid width (m)

   GF%InvDZ = 1.0_ReKi/real(header%dz, ReKi)     ! 1/dz
   GF%ZHWid = 0.5_ReKi*(GF%NZGrids - 1)/GF%InvDZ   ! half the grid height (m)

   GF%MeanWS = real(header%mws, ReKi)                 ! the mean wind speed at hub height (m/s)
   GF%InvMWS = 1.0_ReKi/GF%MeanWS                      ! inverse of mean wind speed

   GF%RefHeight = real(header%ref_height, ReKi)       ! height of the hub (m)
   GF%GridBase = real(header%grid_base_height, ReKi)  ! height of the bottom of the grid (m)

   if (GF%Periodic) then
      GF%InitXPosition = 0                       ! start at the hub
      GF%TotalTime = GF%NSteps*GF%DTime
   else
      GF%InitXPosition = GF%YHWid                 ! start half the grid width ahead of the turbine
      GF%TotalTime = (GF%NSteps - 1)*GF%DTime
   end if

   GF%WindProfileType = WindProfileType_None     ! unused for turbsim
   GF%PLExp = 0                                  ! unused for turbsim
   GF%Z0 = 0                                     ! unused for turbsim

   !----------------------------------------------------------------------------
   ! Binary scaling factors from header
   !----------------------------------------------------------------------------

   ! Wind component slope for scaling, REAL(4)
   Vslope = [header%VslopeX, header%VslopeY, header%VslopeZ]

   ! Wind component offset for scaling, REAL(4)
   Voffset = [header%VoffsetX, header%VoffsetY, header%VoffsetZ]

   !----------------------------------------------------------------------------------------------
   ! Read the description string: "Generated by TurbSim (vx.xx, dd-mmm-yyyy) on dd-mmm-yyyy at hh:mm:ss."
   !----------------------------------------------------------------------------------------------

   ! The number of characters in the description string, max 200, INT(4)
   NChar = header%DescLen

   ! Read description bytes, INT(1)
   allocate (character(NChar)::DescStr)
   read (WindFileUnit, IOSTAT=TmpErrStat) DescStr
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading description line in the FF binary file "' &
                      //TRIM(InitInp%WindFileName)//'."', ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Allocate arrays for the grid-field grid and tower if applicable
   !----------------------------------------------------------------------------

   ! Allocate storage for grid-field velocity data
   call AllocAry(GF%Vel, GF%NComp, GF%NYGrids, GF%NZGrids, GF%NSteps, &
                 'grid-field velocity data', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for raw grid-field velocity for each time step
   allocate (VelRaw(GF%NComp, GF%NYGrids, GF%NZGrids), stat=TmpErrStat)
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, "error allocating grid-field time step velocity data", &
                      ErrStat, ErrMsg, RoutineName)
   end if
   if (ErrStat >= AbortErrLev) return

   ! If tower grids specified
   if (GF%NTGrids > 0) then

      ! Allocate storage for tower velocity data
      call AllocAry(GF%VelTower, GF%NComp, GF%NTGrids, GF%NSteps, &
                    'tower wind velocity data.', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Allocate storage for raw tower data for each timestep
      allocate (TwrRaw(GF%NComp, GF%NTGrids), stat=TmpErrStat)
      if (TmpErrStat /= 0) then
         call SetErrStat(ErrID_Fatal, "error allocating tower time step velocity data", &
                         ErrStat, ErrMsg, RoutineName)
      end if
      if (ErrStat >= AbortErrLev) return
   end if

   !----------------------------------------------------------------------------
   ! Read the 16-bit raw velocity data and scale it to 32-bit real
   !----------------------------------------------------------------------------

   ! This could take a while, so we'll write a message indicating what's going on:
   call WrScr(NewLine//'   Reading a ' &
              //TRIM(Num2LStr(GF%NYGrids))//'x' &
              //TRIM(Num2LStr(GF%NZGrids))// &
              ' grid ('//TRIM(Num2LStr(GF%YHWid*2))//' m wide, '// &
              TRIM(Num2LStr(GF%GridBase))//' m to '// &
              TRIM(Num2LStr(GF%GridBase + GF%ZHWid*2))// &
              ' m above ground) with a characteristic wind speed of '// &
              TRIM(Num2LStr(GF%MeanWS))//' m/s. '//TRIM(DescStr))

   ! Loop through time steps
   do IT = 1, GF%NSteps

      ! Read grid-field raw wind data (normalized) comprised of 2-byte integers, INT(2)
      ! Indices are Velocity components, Y coordinates, Z coordinates
      read (WindFileUnit, IOSTAT=TmpErrStat) VelRaw
      if (TmpErrStat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading grid wind components in the FF binary file "'// &
                         TRIM(InitInp%WindFileName)//'."', ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Loop through wind components (U, V, W), calculate de-normalized velocity (m/s)
      do IC = 1, 3
         GF%Vel(IC, :, :, IT) = (real(VelRaw(IC, :, :), SiKi) - Voffset(IC))/VSlope(IC)
      end do !IC

      ! Read tower raw wind data (normalized) comprised of 2-byte integers, INT(2)
      ! Indices are Velocity components, Z coordinates
      if (GF%NTGrids > 0) then
         read (WindFileUnit, IOSTAT=TmpErrStat) TwrRaw
         if (TmpErrStat /= 0) then
            call SetErrStat(ErrID_Fatal, ' Error reading tower wind components in the FF binary file "'// &
                            TRIM(InitInp%WindFileName)//'."', ErrStat, ErrMsg, RoutineName)
            return
         end if

         ! Loop through wind components (U, V, W), calculate de-normalized velocity (m/s)
         do IC = 1, 3
            GF%VelTower(IC, :, IT) = (real(TwrRaw(IC, :), SiKi) - Voffset(IC))/VSlope(IC)
         end do
      end if
   end do

   !----------------------------------------------------------------------------
   ! Close the file
   !----------------------------------------------------------------------------

   close (WindFileUnit)

   if (GF%Periodic) then
      call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(GF%NSteps))//' time steps of '// &
                 TRIM(Num2LStr(GF%Rate))//'-Hz grid-field data (period of '// &
                 TRIM(Num2LStr(GF%DTime*(GF%NSteps)))//' seconds).')
   else
      call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(GF%NSteps))//' time steps of '// &
                 TRIM(Num2LStr(GF%Rate))//'-Hz grid-field data ('// &
                 TRIM(Num2LStr(GF%DTime*(GF%NSteps - 1)))//' seconds).')
   end if

   !----------------------------------------------------------------------------
   ! Store wind file metadata
   !----------------------------------------------------------------------------

   call PopulateWindFileDatFromGrid3DField(GF, InitInp%WindFileName, 3, GF%NTGrids > 0, FileDat)

   !----------------------------------------------------------------------------
   ! Write the summary file
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A)') 'TurbSim wind type.  Read by InflowWind sub-module ' &
         //TRIM(FlowField_IO_Ver%Name)//' '//TRIM(FlowField_IO_Ver%Ver)
      write (SumFileUnit, '(A)') TRIM(TmpErrMsg)
      write (SumFileUnit, '(5x,A)') 'FileName:                    '//TRIM(InitInp%WindFileName)
      write (SumFileUnit, '(5x,A29,I3)') 'Binary file format id:       ', GF%WindFileFormat
      write (SumFileUnit, '(5x,A29,G12.4)') 'Reference height (m):        ', GF%RefHeight
      write (SumFileUnit, '(5x,A29,G12.4)') 'Timestep (s):                ', GF%DTime
      write (SumFileUnit, '(5x,A29,I12)') 'Number of timesteps:         ', GF%NSteps
      write (SumFileUnit, '(5x,A29,G12.4)') 'Mean windspeed (m/s):        ', GF%MeanWS
      write (SumFileUnit, '(5x,A29,L1)') 'Windfile is periodic:        ', GF%Periodic
      write (SumFileUnit, '(5x,A29,L1)') 'Windfile includes tower:     ', GF%NTGrids > 0

      if (GF%Periodic) then
         write (SumFileUnit, '(5x,A)') 'Time range (s):              [ '// &
            TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(GF%TotalTime))//' ]'
      else  ! Shift the time range to compensate for the shifting of the wind grid
         write (SumFileUnit, '(5x,A)') 'Time range (s):              [ '// &
            TRIM(Num2LStr(-GF%InitXPosition*GF%InvMWS))//' : '// &
            TRIM(Num2LStr(GF%TotalTime - GF%InitXPosition*GF%InvMWS))//' ]'
      end if

      write (SumFileUnit, '(5x,A)') 'Y range (m):                 [ '// &
         TRIM(Num2LStr(-GF%YHWid))//' : '//TRIM(Num2LStr(GF%YHWid))//' ]'

      if (GF%NTGrids > 0) then
         write (SumFileUnit, '(5x,A)') 'Z range (m):                 [ '// &
            TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(GF%RefHeight + GF%ZHWid))//' ]'
      else
         write (SumFileUnit, '(5x,A)') 'Z range (m):                 [ '// &
            TRIM(Num2LStr(GF%RefHeight - GF%ZHWid))//' : '//TRIM(Num2LStr(GF%RefHeight + GF%ZHWid))//' ]'
      end if

      ! Get IO status for unit
      inquire (SumFileUnit, iostat=TmpErrStat)
      if (TmpErrStat /= 0_IntKi) then
         call SetErrStat(ErrID_Fatal, 'Error writing to summary file.', ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

end subroutine

subroutine HAWC_Init(InitInp, SumFileUnit, GF, FileDat, ErrStat, ErrMsg)

   type(HAWCInitInputType), intent(in)    :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(Grid3DFieldType), intent(out)       :: GF
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter       :: RoutineName = "Read_HAWC"
   integer(IntKi)                :: WindFileUnit
   real(SiKi), allocatable       :: VelRaw(:, :)   ! grid-field data for one timestep
   integer                       :: IC                ! Loop counter for the number of wind components
   integer                       :: IX, IY, IZ        ! Loop counters for the number of grid points in the X,Y,Z directions
   real(DbKi)                    :: vMean             ! average wind speeds over time at target position
   real(DbKi)                    :: vSum2             ! sum of wind speeds squared
   real(ReKi)                    :: ActualSigma       ! computed standard deviation
   real(ReKi)                    :: ScaleFactors(3)   ! scale factors
   integer(IntKi)                :: TmpErrStat        ! temporary error status
   character(ErrMsgLen)          :: TmpErrMsg

   !----------------------------------------------------------------------------
   ! Initialize variables
   !----------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg = ""

   GF%WindFileFormat = 0
   GF%Periodic = .true.
   GF%InterpTower = .true.
   GF%AddMeanAfterInterp = .true.

   GF%DTime = InitInp%dx/InitInp%URef
   GF%Rate = 1.0_ReKi/GF%DTime

   GF%NComp = 3
   GF%NYGrids = InitInp%ny
   GF%NZGrids = InitInp%nz
   GF%NTGrids = 0
   GF%NSteps = InitInp%nx

   GF%YHWid = 0.5_ReKi*InitInp%dy*(GF%NYGrids - 1)
   GF%InvDY = 1.0/InitInp%dy

   GF%ZHWid = 0.5_ReKi*InitInp%dz*(GF%NZGrids - 1)
   GF%InvDZ = 1.0_ReKi/InitInp%dz

   GF%MeanWS = InitInp%URef
   GF%InvMWS = 1.0_ReKi/GF%MeanWS

   GF%RefHeight = InitInp%RefHt
   GF%GridBase = GF%RefHeight - GF%ZHWid

   GF%InitXPosition = InitInp%XOffset
   GF%TotalTime = GF%NSteps*InitInp%dx/GF%MeanWS

   GF%WindProfileType = InitInp%WindProfileType
   GF%Z0 = InitInp%Z0
   GF%PLExp = InitInp%PLExp

   ScaleFactors = 0.0_ReKi

   !----------------------------------------------------------------------------
   ! Validate arguments
   !----------------------------------------------------------------------------

   if ((InitInp%ScaleMethod == ScaleMethod_Direct) .and. any(InitInp%SF(:GF%NComp) < 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'Turbulence scaling factors must not be negative.', ErrStat, ErrMsg, RoutineName)
   elseif ((InitInp%ScaleMethod == ScaleMethod_StdDev) .and. any(InitInp%SigmaF(:GF%NComp) < 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'Turbulence standard deviations must not be negative.', ErrStat, ErrMsg, RoutineName)
   else if ((InitInp%ScaleMethod /= 0) .and. (InitInp%ScaleMethod /= 1) .and. (InitInp%ScaleMethod /= 2)) then
      call SetErrStat(ErrID_Fatal, 'Turbulence scaling method must be 0, 1, or 2.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%RefHt < 0.0_ReKi .or. EqualRealNos(GF%RefHeight, 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'The grid reference height must be larger than 0.', ErrStat, ErrMsg, RoutineName)
   else if (GF%MeanWS < 0.0_ReKi) then
      call SetErrStat(ErrID_Fatal, 'The reference wind speed must not be negative.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%nx < 1) then
      call SetErrStat(ErrID_Fatal, 'Number of grid points in the X direction must be at least 1.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%ny < 1) then
      call SetErrStat(ErrID_Fatal, 'Number of grid points in the Y direction must be at least 1.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%nz < 1) then
      call SetErrStat(ErrID_Fatal, 'Number of grid points in the Z direction must be at least 1.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%dx < 0.0_ReKi .or. EqualRealNos(InitInp%dx, 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'The grid spacing in the X direction must be larger than 0.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%dy < 0.0_ReKi .or. EqualRealNos(InitInp%dy, 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'The grid spacing in the Y direction must be larger than 0.', ErrStat, ErrMsg, RoutineName)
   else if (InitInp%dz < 0.0_ReKi .or. EqualRealNos(InitInp%dz, 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'The grid spacing in the Z direction must be larger than 0.', ErrStat, ErrMsg, RoutineName)
   else if (GF%GridBase < 0.0_ReKi) then
      call SetErrStat(ErrID_Severe, 'WARNING: The bottom of the grid is located at a height of '// &
                      trim(Num2LStr(GF%GridBase))//' meters, which is below the ground.'// &
                      ' Winds below the ground will be set to 0.', ErrStat, ErrMsg, RoutineName)
   end if
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Allocate storage for grid-field velocity data
   !----------------------------------------------------------------------------

   ! Allocate storage for grid-field velocity data
   call AllocAry(GF%Vel, GF%NComp, GF%NYGrids, GF%NZGrids, GF%NSteps, &
                 'grid-field velocity data', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for raw grid-field velocity for each time step
   allocate (VelRaw(GF%NZGrids, GF%NYGrids), stat=TmpErrStat)
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, "error allocating grid-field time step velocity data", &
                      ErrStat, ErrMsg, RoutineName)
   end if
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Loop through files and read data
   !----------------------------------------------------------------------------

   ! Display message indicating that file is being read
   call WrScr(NewLine//'   Reading HAWC wind files with grids of '// &
              TRIM(Num2LStr(GF%NSteps))//' x '//TRIM(Num2LStr(GF%NYGrids))//' x '//TRIM(Num2LStr(GF%NZGrids))//' points'// &
              ' ('//TRIM(Num2LStr(GF%YHWid*2))//' m wide, '//TRIM(Num2LStr(GF%GridBase))//' m to '// &
              TRIM(Num2LStr(GF%GridBase + GF%ZHWid*2))// &
              ' m above ground) with a characteristic wind speed of '//TRIM(Num2LStr(GF%MeanWS))//' m/s. ')

   ! Get a unit number to use for the wind file
   call GetNewUnit(WindFileUnit, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Loop through wind components (X, Y, Z)
   do IC = 1, GF%NComp

      ! Open wind file for this component
      call OpenBInpFile(WindFileUnit, InitInp%WindFileName(IC), TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Loop through time steps
      do IX = 1, GF%NSteps

         ! Read file data for this timestep (Z(1:N),Y(N:1),C(1:3))
         read (WindFileUnit, IOSTAT=TmpErrStat) VelRaw
         if (TmpErrStat /= 0) then
            TmpErrMsg = ' Error reading binary data from "'//TRIM(InitInp%WindFileName(IC)) &
                        //'". I/O error '//TRIM(Num2LStr(TmpErrStat)) &
                        //' occurred at IX='//TRIM(Num2LStr(IX))//'.'
            close (WindFileUnit)
            call SetErrStat(ErrID_Fatal, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            return
         end if

         ! Reorganize raw data into grid-field array (reverse Y indices)
         do IZ = 1, GF%NZGrids
            ! Vel(NFFComp, NYGrids, NZGrids, NSteps)
            ! VelRaw(NZGrids, NYGrids)
            GF%Vel(IC, :, IZ, IX) = VelRaw(IZ, GF%NYGrids:1:-1)
         end do
      end do

      ! Close file
      close (WindFileUnit)

   end do

   !----------------------------------------------------------------------------
   ! Scale turbulence
   !----------------------------------------------------------------------------

   select case (InitInp%ScaleMethod)

   case (ScaleMethod_None) ! No scaling
      ScaleFactors = 1.0_ReKi

   case (ScaleMethod_Direct) ! Apply scale factors from file
      ScaleFactors = InitInp%SF

   case (ScaleMethod_StdDev) ! Use scale factor to get requested sigma

      ! find the center point of the grid (if we don't have an odd number
      ! of grid points, we'll pick the point closest to the center)
      iz = ishft(InitInp%nz + 1, -1) ! integer divide by 2
      iy = ishft(InitInp%ny + 1, -1) ! integer divide by 2

      ! Loop through components
      do ic = 1, GF%NComp

         ! Mean of component for all time
         vMean = sum(GF%Vel(ic, iy, iz, :))/InitInp%nx

         ! Sum of square of component value for all time
         vSum2 = dot_product(GF%Vel(ic, iy, iz, :), GF%Vel(ic, iy, iz, :))

         ! Standard deviation of component
         ActualSigma = real(sqrt(abs(vSum2/InitInp%nx - vMean**2)), kind(ActualSigma))

         ! If actual sigma is nearly zero,
         if (EqualRealNos(ActualSigma, 0.0_ReKi)) then
            ScaleFactors(ic) = 0.0_ReKi
            if (.not. EqualRealNos(InitInp%SigmaF(ic), 0.0_ReKi)) then
               call SetErrStat(ErrID_Fatal,"Computed standard deviation is zero; cannot scale to achieve target non-zero standard deviation.", ErrStat, ErrMsg, RoutineName )                  
            end if
         else
            ScaleFactors(ic) = InitInp%SigmaF(ic)/ActualSigma
         end if
      end do
   end select

   ! Loop through components and apply scale factors
   do ic = 1, GF%NComp
      GF%Vel(ic, :, :, :) = real(ScaleFactors(ic)*GF%Vel(ic, :, :, :), SiKi)
   end do

   !----------------------------------------------------------------------------
   ! Remove the U component mean wind speed
   !----------------------------------------------------------------------------

   ! If scaling method is not none, remove mean value of X component at each grid point
   if (InitInp%ScaleMethod /= ScaleMethod_None) then
      do iz = 1, GF%NZGrids
         do iy = 1, GF%NYGrids
            vMean = sum(GF%Vel(1, iy, iz, :))/GF%NSteps
            GF%Vel(1, iy, iz, :) = real(GF%Vel(1, iy, iz, :) - vMean, SiKi)
         end do
      end do
   end if

   !----------------------------------------------------------------------------
   ! Store wind file metadata
   !----------------------------------------------------------------------------

   call PopulateWindFileDatFromGrid3DField(GF, InitInp%WindFileName(1), 5, .false., FileDat)

   !----------------------------------------------------------------------------
   ! Write the summary file
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A)') 'HAWC wind type.  Read by InflowWind sub-module FlowField_IO'

      write (SumFileUnit, '(A34,G12.4)') '     Reference height (m):        ', GF%RefHeight
      write (SumFileUnit, '(A34,G12.4)') '     Timestep (s):                ', GF%DTime
      write (SumFileUnit, '(A34,I12)') '     Number of timesteps:         ', GF%NSteps
      write (SumFileUnit, '(A34,G12.4)') '     Mean windspeed (m/s):        ', GF%MeanWS
      write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
         TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(GF%TotalTime))//' ]'
      write (SumFileUnit, '(A)') '     X range (m):                 [ '// &
         TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(GF%TotalTime*GF%MeanWS))//' ]'
      write (SumFileUnit, '(A)') '     Y range (m):                 [ '// &
         TRIM(Num2LStr(-GF%YHWid))//' : '//TRIM(Num2LStr(GF%YHWid))//' ]'
      write (SumFileUnit, '(A)') '     Z range (m):                 [ '// &
         TRIM(Num2LStr(GF%GridBase))//' : '//TRIM(Num2LStr(GF%GridBase + GF%ZHWid*2.0))//' ]'

      write (SumFileUnit, '(A)') 'Scaling factors used:'
      write (SumFileUnit, '(A)') '  u           v           w       '
      write (SumFileUnit, '(A)') '----------  ----------  ----------'
      write (SumFileUnit, '(F10.3,2x,F10.3,2x,F10.3)') ScaleFactors
   end if

end subroutine

!> User_Init initializes a user defined wind field.
subroutine User_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)

   type(UserInitInputType), intent(in)    :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(UserFieldType), intent(out)       :: UF
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = "User_Init"

   ErrStat = ErrID_None
   ErrMsg = ""

   UF%Dummy = 1

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)') InitInp%Dummy
   end if

end subroutine

!> Grid4D_Init initializes a wind field defined by a 4D grid.
subroutine Grid4D_Init(InitInp, SumFileUnit, EGF, FileDat, ErrStat, ErrMsg)

   type(Grid4DInitInputType), intent(in) :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(Grid4DFieldType), intent(out)    :: EGF
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = "Grid4D_Init"
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Initialize field from inputs
   EGF%n = InitInp%n
   EGF%delta = InitInp%delta
   EGF%pZero = InitInp%pZero
   EGF%TimeStart = 0.0_ReKi

   ! uvw velocity components at x,y,z,t coordinates
   call AllocAry(EGF%Vel, 3, EGF%n(1), EGF%n(2), EGF%n(3), EGF%n(4), &
                 'External Grid Velocity', TmpErrStat, TmpErrMsg)
   call SetErrStat(ErrStat, ErrMsg, TmpErrStat, TmpErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Initialize velocities to zero
   EGF%Vel = 0.0_SiKi

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)') InitInp%n
   end if

end subroutine

subroutine PopulateWindFileDatFromGrid3DField(Grid3DField, FileName, WindType, HasTower, FileDat)

   type(Grid3DFieldType), intent(in)  :: Grid3DField
   character(*), intent(in)         :: FileName
   integer(IntKi), intent(in)       :: WindType
   logical, intent(in)              :: HasTower
   type(WindFileDat), intent(out)   :: FileDat

   FileDat%FileName = FileName
   FileDat%WindType = WindType
   FileDat%RefHt = Grid3DField%RefHeight
   FileDat%RefHt_Set = .true.
   FileDat%DT = Grid3DField%DTime
   FileDat%NumTSteps = Grid3DField%NSteps
   FileDat%ConstantDT = .true.

   if (Grid3DField%Periodic) then
      FileDat%TRange = [0.0_ReKi, Grid3DField%TotalTime]
      FileDat%TRange_Limited = .false.
   else  ! Shift the time range to compensate for the shifting of the wind grid
      FileDat%TRange = [0.0_ReKi, Grid3DField%TotalTime] - Grid3DField%InitXPosition*Grid3DField%InvMWS
      FileDat%TRange_Limited = .true.
   end if

   FileDat%YRange = [-Grid3DField%YHWid, Grid3DField%YHWid]
   FileDat%YRange_Limited = .true.   ! Hard boundaries enforced in y-direction

   ! If has tower data
   if (HasTower) then
      FileDat%ZRange = [0.0_Reki, Grid3DField%RefHeight + Grid3DField%ZHWid]
   else
      FileDat%ZRange = [Grid3DField%GridBase, Grid3DField%GridBase + Grid3DField%ZHWid*2.0]
   end if

   FileDat%ZRange_Limited = .true.
   FileDat%BinaryFormat = Grid3DField%WindFileFormat
   FileDat%IsBinary = .true.
   FileDat%MWS = Grid3DField%MeanWS

   FileDat%TI = 0.0_ReKi
   FileDat%TI_listed = .false.

end subroutine

! subroutine Read_Bladed_Native(p, InitInp, InitOut, WindFileUnit, SumFileUnit, ErrStat, ErrMsg)

!    ! Passed Variables:
!    type(IfW_Interp_ParameterType), intent(inout)   :: p        !< Parameters
!    type(IfW_Interp_InitInputType), intent(inout)   :: InitInp  !< Initialization inputs
!    type(IfW_Interp_InitOutputType), intent(out)    :: InitOut  !< Initialization outputs
!    integer(IntKi), intent(in)       :: WindFileUnit            !< unit number for the wind file
!    integer(IntKi), intent(in)       :: SumFileUnit             !< unit number for the summary file
!    integer(IntKi), intent(out)      :: ErrStat                 !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)        :: ErrMsg                  !< message about the error encountered

!    ! Local Parameters:
!    character(*), parameter          :: RoutineName = "Read_Bladed_Native"

!    ! Local Variables:
!    character(1024)                  :: WindFileName
!    logical                          :: NativeFormat      !< file is in native bladed format
!    logical                          :: TowerFileExist    !< tower file exists
!    character(1024)                  :: BinFileName
!    character(1024)                  :: PriPath
!    character(1028)                  :: SumFile           ! length is LEN(ParamData%WindFileName) + the 4-character extension.
!    character(1028)                  :: TwrFile           ! length is LEN(ParamData%WindFileName) + the 4-character extension.

!    integer(IntKi)                   :: ScaleMethod       ! Turbulence scaling method [0=none, 1=direct scaling, 2= calculate scaling factor based on a desired standard deviation]
!    real(ReKi)                       :: NatTI(3), TI      ! Turbulence scaling factor for each direction [ScaleMethod=1]
!    real(ReKi)                       :: SF(3)             ! Turbulence scaling factor for each direction [ScaleMethod=1]
!    real(ReKi)                       :: SigmaF(3)         ! Turbulence standard deviation to calculate scaling from in each direction [ScaleMethod=2]

!    real(ReKi)                       :: FFXDelt, FFYDelt, FFZDelt
!    real(ReKi)                       :: XOffset
!    real(ReKi)                       :: UBar
!    real(ReKi)                       :: ZCenter

!    logical                          :: CWise
!    logical                          :: LHR               ! Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
!    logical                          :: exists

!    integer(B2Ki)                    :: FileFormat
!    integer(IntKi)                   :: TmpErrStat        ! temporary error status
!    character(ErrMsgLen)             :: TmpErrMsg         ! temporary error message

!    ! Wind file is specified directly in input file
!    WindFileName = InitInp%WindFileNames(1)

!    ! No tower file is used
!    TowerFileExist = .false.

!    ! Read the native bladed summary file
!    call Read_Bladed_Native_Summary( &
!       WindFileName, p%MeanFFWS, p%RefHt, InitOut%TI, InitOut%PropagationDir, &
!       InitOut%VFlowAngle, BinFileName, p%PLExp, p%InitXPosition, ErrStat, ErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    ! Binary file will be relative to the path where the primary input file is located.
!    if (pathIsRelative(BinFileName)) then
!       call GetPath(WindFileName, PriPath)
!       BinFileName = TRIM(PriPath)//TRIM(BinFileName)
!    end if

!    ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
!    if (InitInp%FixedWindFileRootName) then
!       if (InitInp%TurbineID == 0) then
!          ! .TRUE. for the FAST.Farm low-resolution domain
!          BinFileName = TRIM(BinFileName)//TRIM(PathSep)//'Low'
!       else
!          ! FAST.Farm high-resolution domain(s)
!          BinFileName = TRIM(BinFileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
!       end if
!    end if

!    p%WindProfileType = WindProfileType_PL
!    p%Periodic = .true.

!    InitInp%ScaleMethod = ScaleMethod_StdDev
!    InitInp%SigmaF = InitOut%TI*p%MeanFFWS
!    InitInp%SF = InitInp%SigmaF

!    InitInp%RefHt = p%RefHt
!    InitInp%URef = p%MeanFFWS
!    InitInp%WindProfileType = WindProfileType_PL

!    CWise = .false.
!    ZCenter = p%RefHt
!    TI = 100.0_ReKi
!    UBar = 0.0_ReKi
!    LHR = .true.

! end subroutine

! !> Read_Bladed reads the bladed full-field wind file.
! subroutine Read_Bladed_Binary(p, InitInp, InitOut, WindFileUnit, SumFileUnit, ErrStat, ErrMsg)

!    ! Passed Variables:
!    type(IfW_Interp_ParameterType), intent(inout)   :: p        !< Parameters
!    type(IfW_Interp_InitInputType), intent(inout)   :: InitInp  !< Initialization inputs
!    type(IfW_Interp_InitOutputType), intent(out)    :: InitOut  !< Initialization outputs
!    integer(IntKi), intent(in)       :: WindFileUnit            !< unit number for the wind file
!    integer(IntKi), intent(in)       :: SumFileUnit             !< unit number for the summary file
!    integer(IntKi), intent(out)      :: ErrStat                 !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)        :: ErrMsg                  !< message about the error encountered

!    ! Local Parameters:
!    character(*), parameter          :: RoutineName = "Read_Bladed_Binary"
!    logical, parameter               :: NativeFormat = .false.

!    ! Local Variables:
!    character(1024)                  :: WindFileName
!    logical                          :: TowerFileExist    !< tower file exists
!    character(1024)                  :: BinFileName
!    character(1024)                  :: PriPath
!    character(1028)                  :: SumFile           ! length is LEN(ParamData%WindFileName) + the 4-character extension.
!    character(1028)                  :: TwrFile           ! length is LEN(ParamData%WindFileName) + the 4-character extension.

!    integer(IntKi)                   :: ScaleMethod       ! Turbulence scaling method [0=none, 1=direct scaling, 2= calculate scaling factor based on a desired standard deviation]

!    real(ReKi)                       :: TI(3)             ! Turbulence scaling factor for each direction [ScaleMethod=1]
!    real(ReKi)                       :: SF(3)             ! Turbulence scaling factor for each direction [ScaleMethod=1]
!    real(ReKi)                       :: SigmaF(3)         ! Turbulence standard deviation to calculate scaling from in each direction [ScaleMethod=2]

!    real(ReKi)                       :: XOffset
!    real(ReKi)                       :: UBar
!    real(ReKi)                       :: ZCenter

!    logical                          :: CWise
!    logical                          :: LHR               ! Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
!    logical                          :: exists

!    integer(IntKi)                   :: TmpErrStat        ! temporary error status
!    character(ErrMsgLen)             :: TmpErrMsg         ! temporary error message

!    !----------------------------------------------------------------------------
!    ! Initialize variables
!    !----------------------------------------------------------------------------

!    ErrStat = ErrID_None
!    ErrMsg = ""

!    InitOut%PropagationDir = 0.0_ReKi
!    InitOut%VFlowAngle = 0.0_ReKi
!    BinFileName = InitInp%WindFileNames(1)

!    ! Create summary and tower file paths from bin file path
!    call GetRoot(BinFileName, SumFile)
!    TwrFile = TRIM(SumFile)//'.twr'
!    SumFile = TRIM(SumFile)//'.sum'

!    ! Read the summary file to get necessary scaling information
!    call Read_Bladed_TurbSim_Summary(WindFileUnit, SumFile, CWise, ZCenter, TI, UBar, &
!                                     p%RefHt, p%Periodic, LHR, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
!    if (InitInp%FixedWindFileRootName) then
!       if (InitInp%TurbineID == 0) then
!          ! .TRUE. for the FAST.Farm low-resolution domain
!          WindFileName = TRIM(WindFileName)//TRIM(PathSep)//'Low'
!       else
!          ! FAST.Farm high-resolution domain(s)
!          WindFileName = TRIM(WindFileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
!       end if
!    end if

!    WindFileName = TRIM(WindFileName)//'.wnd'

!    InitOut%PropagationDir = 0.0_ReKi
!    InitOut%VFlowAngle = 0.0_ReKi
!    BinFileName = WindFileName

!    !----------------------------------------------------------------------------
!    ! Write the summary file
!    !----------------------------------------------------------------------------

!    if (SumFileUnit > 0) then
!       write (SumFileUnit, '(A)')
!       write (SumFileUnit, '(A)') 'Bladed-style wind type.  Read by InflowWind sub-module '// &
!          TRIM(IfW_Interp_Ver%Name)//' '//TRIM(IfW_Interp_Ver%Ver)
!       write (SumFileUnit, '(A)') TRIM(TmpErrMsg)
!       write (SumFileUnit, '(A)') '     FileName:                    '//TRIM(InitInp%WindFileNames(1))
!       write (SumFileUnit, '(A34,I3)') '     Binary file format id:       ', p%WindFileFormat
!       write (SumFileUnit, '(A34,G12.4)') '     Reference height (m):        ', p%RefHt
!       write (SumFileUnit, '(A34,G12.4)') '     Timestep (s):                ', p%FFDTime
!       write (SumFileUnit, '(A34,I12)') '     Number of timesteps:         ', p%NSteps
!       write (SumFileUnit, '(A34,G12.4)') '     Mean windspeed (m/s):        ', p%MeanFFWS
!       write (SumFileUnit, '(A)') '     Characteristic TI:            [ '// &
!          TRIM(Num2LStr(TI(1)))//', '//TRIM(Num2LStr(TI(2)))//', '//TRIM(Num2LStr(TI(3)))//' ] '
!       write (SumFileUnit, '(A34,L1)') '     Windfile is periodic:        ', p%Periodic
!       write (SumFileUnit, '(A34,L1)') '     Windfile includes tower:     ', p%NTGrids > 0

!       if (p%Periodic) then
!          write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
!             TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(p%TotalTime))//' ]'
!       else  ! Shift the time range to compensate for the shifting of the wind grid
!          write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
!             TRIM(Num2LStr(-p%InitXPosition*p%InvMFFWS))//' : '// &
!             TRIM(Num2LStr(p%TotalTime - p%InitXPosition*p%InvMFFWS))//' ]'
!       end if

!       write (SumFileUnit, '(A)') '     Y range (m):                 [ '// &
!          TRIM(Num2LStr(-p%FFYHWid))//' : '//TRIM(Num2LStr(p%FFYHWid))//' ]'

!       if (p%NTGrids > 0) then
!          write (SumFileUnit, '(A)') '     Z range (m):                 [ '// &
!             TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(p%RefHt + p%FFZHWid))//' ]'
!       else
!          write (SumFileUnit, '(A)') '     Z range (m):                 [ '// &
!             TRIM(Num2LStr(p%RefHt - p%FFZHWid))//' : '//TRIM(Num2LStr(p%RefHt + p%FFZHWid))//' ]'
!       end if

!    end if

! end subroutine

! subroutine Read_Bladed_Wind_File(p, WindFileUnit, BinFileName, NativeFormat, UBar, ZCenter, TI, ErrStat, ErrMsg)

!    ! Passed Variables:
!    type(IfW_Interp_ParameterType), intent(inout)  :: p         !< Parameters
!    integer(IntKi), intent(in)       :: WindFileUnit            !< unit number for the wind file
!    character(*), intent(in)         :: BinFileName
!    logical, intent(in)              :: NativeFormat
!    real(ReKi), intent(in)           :: UBar
!    real(ReKi), intent(in)           :: ZCenter
!    real(ReKi), intent(inout)        :: TI(3)
!    integer(IntKi), intent(out)      :: ErrStat                 !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)        :: ErrMsg                  !< message about the error encountered

!    ! Local Parameters:
!    character(*), parameter          :: RoutineName = "Read_Bladed_Wind_File"

!    ! Local Variables:
!    integer(B2Ki)                    :: FileFormat
!    real(ReKi)                       :: FFXDelt, FFYDelt, FFZDelt
!    real(ReKi)                       :: BinTI(3)          ! Turbulence scaling factor for each direction [ScaleMethod=1]
!    integer(IntKi)                   :: TmpErrStat        ! temporary error status
!    character(ErrMsgLen)             :: TmpErrMsg         ! temporary error message

!    !----------------------------------------------------------------------------
!    ! Open the binary file and read contents
!    !----------------------------------------------------------------------------

!    ! Open the wind file
!    call OpenBInpFile(WindFileUnit, TRIM(BinFileName), TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    ! Read the first binary integer from the file to get info on the type.
!    ! Cannot use library read routines since this is a 2-byte integer.
!    read (WindFileUnit, IOSTAT=TmpErrStat) FileFormat
!    if (TmpErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading first binary integer from file "' &
!                       //TRIM(BinFileName)//'."', ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Save file format in parameters struct
!    p%WindFileFormat = FileFormat

!    ! Read file headers
!    select case (p%WindFileFormat)
!    case (-1, -2, -3)
!       call Read_Bladed_FF_Header0(p, WindFileUnit, NativeFormat, FFXDelt, FFYDelt, FFZDelt, TmpErrStat, TmpErrMsg)
!    case (-99)
!       call Read_Bladed_FF_Header1(p, WindFileUnit, NativeFormat, BinTI, FFXDelt, FFYDelt, FFZDelt, TmpErrStat, TmpErrMsg)
!       where (BinTI > 0) TI = BinTI
!    case default
!       call SetErrStat(ErrID_Fatal, ' This is not a bladed-style binary wind file (binary format identifier: '// &
!                       TRIM(Num2LStr(p%WindFileFormat))//'.  This might be a TurbSim binary wind file.', &
!                       ErrStat, ErrMsg, RoutineName)
!    end select
!    if (ErrStat >= AbortErrLev) return

!    ! If not native format and mean wind speed is different from summary file
!    if ((.not. NativeFormat) .and. (ABS(UBar - p%MeanFFWS) > 0.1)) then
!       call SetErrStat(ErrID_Fatal, ' Error: Incompatible mean hub-height wind speeds in FF wind files. '// &
!                       '(Check that the .sum and .wnd files were generated together.)', ErrStat, ErrMsg, RoutineName)
!       close (WindFileUnit)
!       return
!    end if

!    ! Calculate the height of the bottom of the grid (meters), error if less than zero
!    p%GridBase = ZCenter - p%FFZHWid
!    if (p%GridBase < 0.0_ReKi) then
!       call SetErrStat(ErrID_Severe, 'WARNING: The bottom of the grid is located at a height of '// &
!                       TRIM(Num2LStr(p%GridBase))//' meters, which is below the ground.'// &
!                       ' Winds below the ground will be set to 0.', ErrStat, ErrMsg, RoutineName)
!    end if

!    ! Read binary grid data
!    call Read_Bladed_Grids(p, WindFileUnit, LHR, TI, CWise, NativeFormat, p, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    close (WindFileUnit)
!    if (ErrStat >= AbortErrLev) return

!    !----------------------------------------------------------------------------
!    ! Read tower data if requested
!    !----------------------------------------------------------------------------

!    ! If a tower file was specified
!    if (TowerFileExist) then

!       inquire (FILE=TRIM(TwrFile), EXIST=Exists)

!       ! Double check that the tower file exists and read it.  If it was requested but doesn't exist,
!       ! throw fatal error and exit.
!       if (Exists) then
!          call Read_Bladed_Tower()
!          call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!          if (ErrStat >= AbortErrLev) then
!             close (WindFileUnit)
!             return
!          end if
!       else
!          call SetErrStat(ErrID_Fatal, ' Tower file '//TRIM(TwrFile)//' specified for Bladed full-field '// &
!                          'wind files does not exist.', ErrStat, ErrMsg, RoutineName)
!          close (WindFileUnit)
!          return
!       end if
!    else
!       p%NTGrids = 0
!    end if

!    ! Report native turbulence intensities in percent
!    if (NativeFormat) TI = TI*100.0_ReKi

! end subroutine

! !> Read_Bladed_TurbSim_Summary reads the summary file generated by TurbSim
!    !! to get the normalizing parameters needed to read the bladed binary file
!    !! which was genered by TurbSim.
! subroutine Read_Bladed_TurbSim_Summary(WindFileUnit, FileName, CWise, ZCenter, TI, UBar, RefHt, Periodic, LHR, ErrStat, ErrMsg)

!    integer(IntKi), intent(IN)    :: WindFileUnit   !< unit number for the file to open
!    character(*), intent(IN)      :: FileName       !< name of the summary file
!    logical, intent(OUT)          :: CWise          !< rotation (for reading the order of the binary data)
!    real(ReKi), intent(OUT)       :: ZCenter        !< the height at the center of the grid
!    real(ReKi), intent(OUT)       :: TI(3)          !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
!    real(ReKi), intent(OUT)       :: UBar           !< mean (advection) wind speed
!    real(ReKi), intent(OUT)       :: RefHt          !< Reference height
!    logical, intent(OUT)          :: Periodic       !< rotation (for reading the order of the binary data)
!    logical, intent(OUT)          :: LHR            !< Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
!    integer(IntKi), intent(OUT)   :: ErrStat        !< returns 0 if no error encountered in the subroutine
!    character(*), intent(OUT)     :: ErrMsg         !< holds the error messages

!    character(*), parameter       :: RoutineName = "Read_Bladed_TurbSim_Summary"

!    real(ReKi)                    :: ZGOffset       ! The vertical offset of the turbine on rectangular grid (allows turbulence not centered on turbine hub)
!    integer(IntKi)                :: I, J           ! Counters
!    character(1024)               :: Line           ! temporary storage for reading a line from the file
!    integer(IntKi)                :: ContentSize    ! summary file size
!    character(:), allocatable     :: FileContents   ! summary file contents
!    integer(IntKi)                :: TmpErrStat     ! temporary error status
!    character(ErrMsgLen)          :: TmpErrMsg      ! temporary error message

!    !-------------------------------------------------------------------------
!    ! Initialize variables
!    !-------------------------------------------------------------------------

!    ErrStat = ErrID_None
!    ErrMsg = ''

!    !-------------------------------------------------------------------------
!    ! Read summary file contents
!    !-------------------------------------------------------------------------

!    call OpenBInpFile(WindFileUnit, FileName, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    inquire (WindFileUnit, size=ContentSize)
!    allocate (character(ContentSize) :: FileContents)
!    read (WindFileUnit, IOSTAT=TmpErrStat) FileContents
!    close (WindFileUnit)

!    if (TmpErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading TurbSim summary file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Convert file contents to uppercase
!    call Conv2UC(FileContents)

!    ! Find 'CLOCKWISE' in file
!    i = index(FileContents, 'CLOCKWISE')
!    if (i > 0) then

!       ! Find beginning of line where clockwise appears
!       j = index(FileContents(:i), NewLine, back=.true.) + 1

!       ! Get line with leading spaces removed
!       line = adjustl(FileContents(j:i))

!       ! Get value from beginning of line
!       select case (line(1:1))
!       case ("T", "Y")
!          CWise = .true.
!       case ("F", "N")
!          CWise = .false.
!       end select
!    else
!       CWise = .false.
!    end if

!    ! Find 'HUB HEIGHT' or 'ZHUB' in file
!    i = index(FileContents, 'HUB HEIGHT')
!    if (i == 0) i = index(FileContents, 'ZHUB')
!    if (i == 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading hub height from FF summary file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Find beginning of line where 'HUB HEIGHT' or 'ZHUB' occur
!    j = index(FileContents(:i), NewLine, back=.true.) + 1

!    ! Read value from line
!    read (FileContents(j:), *, IOSTAT=TmpErrStat) RefHt
!    if (TmpErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading hub height from FF summary file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Find the mean speed
!    i = index(FileContents, 'UBAR')
!    if (i == 0) then
!       call SetErrStat(ErrID_Fatal, ' Error find UBar in FF summary file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Get index of next equal sign
!    i = index(FileContents(i:), '=') + i + 1
!    read (FileContents(i:), *, IOSTAT=TmpErrStat) UBar
!    if (TmpErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading UBar from FF summary file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Read turbulence intensities on lines after mean speed
!    do j = 1, 3
!       i = INDEX(FileContents(i:), '=') + i + 1
!       read (FileContents(i:), *, IOSTAT=TmpErrStat) TI(j)
!       if (TmpErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(j))// &
!                          ') from FF summary file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if
!    end do

!    ! Find the grid "HEIGHT OFFSET" (optional, default to zero)
!    i = index(FileContents, "HEIGHT OFFSET")
!    if (i > 0) then
!       i = INDEX(FileContents(i:), '=') + i + 1
!       read (FileContents(i:), *, IOSTAT=TmpErrStat) ZGOffset
!       if (TmpErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading height offset from FF summary file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if
!    else
!       ZGOffset = 0.0_ReKi
!    end if

!    ! Calculate the height of the grid center
!    ZCenter = RefHt - ZGOffset

!    ! If "PERIODIC" is in summary, then file is periodic
!    Periodic = index(FileContents(i:), 'PERIODIC') > 0

!    ! If "BLADED LEFT-HAND RULE" is in file, then file uses LHR
!    LHR = index(FileContents(i:), 'BLADED LEFT-HAND RULE') > 0

! end subroutine Read_Bladed_TurbSim_Summary

! !> This subroutine reads the text summary file to get normalizing parameters, the location of the
! !! grid, and the direction the grid was written to the binary file
! subroutine Read_Bladed_Native_Summary(SumFileName, UBar, RefHt, TI, PropagationDir, &
!                                       VFlowAngle, BinFileName, PLExp, XOffset, ErrStat, ErrMsg)

!    character(*), intent(in)     :: SumFileName          !< Summary file name
!    real(ReKi), intent(out)       :: UBar
!    real(ReKi), intent(out)       :: RefHt
!    real(ReKi), intent(out)       :: TI(3)             ! turbulenc
!    real(ReKi), intent(out)       :: PropagationDir
!    real(ReKi), intent(out)       :: VFlowAngle
!    character(1024), intent(out)  :: BinFileName
!    real(ReKi), intent(out)       :: PLExp
!    real(ReKi), intent(out)       :: XOffset
!    integer(IntKi), intent(out)   :: ErrStat              !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)     :: ErrMsg               !< message about the error encountered

!    character(*), parameter       :: RoutineName = "Read_Bladed_Native_Summary"
!    integer(IntKi), parameter     :: UnEc = -1      ! echo file unit number (set to something else > 0 for debugging)

!    type(FileInfoType)            :: FileInfo       ! The derived type for holding the file information.
!    integer(IntKi)                :: CurLine        ! Current line to parse in FileInfo data structure
!    integer(IntKi)                :: TmpErrStat       ! temporary error status
!    character(ErrMsgLen)          :: TmpErrMsg        ! temporary error message

!    ! Initialize error variables
!    ErrStat = ErrID_None
!    ErrMsg = ''

!    !-------------------------------------------------------------------------
!    ! Open and read the summary file; store data in FileInfo structure.
!    !-------------------------------------------------------------------------

!    call ProcessComFile(SumFileName, FileInfo, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    !-------------------------------------------------------------------------
!    ! Process the lines stored in FileInfo
!    !-------------------------------------------------------------------------

!    CurLine = 1

!    call ParseVar(FileInfo, CurLine, 'UBAR', UBar, TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    call ParseVar(FileInfo, CurLine, 'REFHT', RefHt, TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    call ParseVar(FileInfo, CurLine, 'TI', TI(1), TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    call ParseVar(FileInfo, CurLine, 'TI_V', TI(2), TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    call ParseVar(FileInfo, CurLine, 'TI_W', TI(3), TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    call ParseVar(FileInfo, CurLine, 'WDIR', PropagationDir, TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    PropagationDir = R2D*PropagationDir

!    call ParseVar(FileInfo, CurLine, 'FLINC', VFlowAngle, TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    VFlowAngle = R2D*VFlowAngle ! convert to degrees

!    call ParseVar(FileInfo, CurLine, 'WINDF', BinFileName, TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    call ParseVar(FileInfo, CurLine, 'WSHEAR', PLExp, TmpErrStat, TmpErrMsg, UnEc)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

!    if (ErrStat >= AbortErrLev) return

!    call ParseVar(FileInfo, CurLine, 'XOffset', XOffset, TmpErrStat, TmpErrMsg, UnEc)
!    if (TmpErrStat /= ErrID_None) then
!       XOffset = 0.0_ReKi ! this will be the default if offset is not in the file
!    end if

! end subroutine Read_Bladed_Native_Summary

! subroutine Read_Bladed_FF_Header0(p, WindFileUnit, NativeFormat, FFXDelt, FFYDelt, FFZDelt, ErrStat, ErrMsg)

!    type(IfW_Interp_ParameterType), intent(inout)   :: p              !< Parameters
!    integer(IntKi), intent(IN)                      :: WindFileUnit   !< unit number for the file to open
!    logical, intent(in)                             :: NativeFormat
!    real(ReKi), intent(out)                         :: FFXDelt
!    real(ReKi), intent(out)                         :: FFYDelt
!    real(ReKi), intent(out)                         :: FFZDelt
!    integer(IntKi), intent(out)                     :: ErrStat        !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)                       :: ErrMsg         !< message about the error encountered

!    character(*), parameter                         :: RoutineName = "Read_Bladed_FF_Header0"

!    integer(IntKi)                                  :: TmpErrStat     ! temporary error status
!    character(ErrMsgLen)                            :: TmpErrMsg      ! temporary error message

!    type :: HeaderType
!       integer(B2Ki)  :: DeltaZ, DeltaY, DeltaX
!       integer(B2Ki)  :: NumTimeStepsHalf, MFFWS10
!       integer(B2Ki)  :: zLu, yLu, xLu, dummy, rnd
!       integer(B2Ki)  :: nz1000, ny1000
!    end type

!    type(HeaderType)  :: header
!    integer(B2Ki)     :: dummy(6)

!    ! Read header
!    read (WindFileUnit, IOSTAT=ErrStat) header
!    if (ErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading header 0 from binary FF file.', ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Number of full-field velocity components
!    p%NFFComp = -1*p%WindFileFormat

!    ! Grid Deltas
!    FFZDelt = 0.001_ReKi*real(header%DeltaZ, ReKi)
!    p%InvFFZD = 1.0_ReKi/FFZDelt
!    FFYDelt = 0.001_ReKi*real(header%DeltaY, ReKi)
!    p%InvFFYD = 1.0_ReKi/FFYDelt
!    FFXDelt = 0.001_ReKi*real(header%DeltaX, ReKi)

!    ! Number of time steps
!    p%NSteps = 2*header%NumTimeStepsHalf

!    ! Mean full field wind speed
!    if (.not. NativeFormat) p%MeanFFWS = 0.1*header%MFFWS10
!    p%InvMFFWS = 1.0/p%MeanFFWS
!    p%FFDTime = FFXDelt/p%MeanFFWS
!    p%FFRate = 1.0/p%FFDTime

!    ! Z and Y Grid size
!    p%NZGrids = header%nz1000/1000
!    p%FFZHWid = 0.5*FFZDelt*(p%NZGrids - 1)
!    p%NYGrids = header%ny1000/1000
!    p%FFYHWid = 0.5*FFYDelt*(p%NYGrids - 1)

!    ! Read dummy variables (zLv, yLv, xLv, zLw, yLw, xLw)
!    if (p%NFFComp == 3) then
!       read (WindFileUnit, IOSTAT=TmpErrStat) dummy
!       if (TmpErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading header 0 from binary FF file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if
!    end if

! end subroutine Read_Bladed_FF_Header0

! subroutine Read_Bladed_FF_Header1(p, WindFileUnit, NativeFormat, TI, FFXDelt, FFYDelt, FFZDelt, ErrStat, ErrMsg)

!    type(IfW_Interp_ParameterType), intent(inout)   :: p              !< Parameters
!    integer(IntKi), intent(IN)                      :: WindFileUnit   !< unit number for the file to open
!    logical, intent(in)                             :: NativeFormat
!    real(ReKi), intent(out)                        :: TI(3)          !< turbulence intensities
!    real(ReKi), intent(out)                         :: FFXDelt
!    real(ReKi), intent(out)                         :: FFYDelt
!    real(ReKi), intent(out)                         :: FFZDelt
!    integer(IntKi), intent(out)                     :: ErrStat        !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)                       :: ErrMsg         !< message about the error encountered

!    character(*), parameter                         :: RoutineName = "Read_Bladed_FF_Header1"

!    integer(B2Ki)                                   :: TurbType
!    integer(IntKi)                                  :: i

!    type :: Turb4Type
!       integer(B4Ki)  :: NComp
!       real(SiKi)  :: Latitude, RoughLen, RefHeight, TurbInt(3)
!    end type

!    type :: Turb78Type
!       integer(B4Ki)  :: HeaderSize, NComp
!    end type

!    type :: Turb7Type
!       real(SiKi)     :: CoherenceDecay, CoherenceScale
!    end type

!    type :: Turb8Type
!       real(SiKi)     :: dummy1(6)
!       integer(B4Ki)  :: dummy2(3)
!       real(SiKi)     :: dummy3(2)
!       integer(B4Ki)  :: dummy4(3)
!       real(SiKi)     :: dummy5(2)
!    end type

!    type :: Sub1Type
!       real(SiKi)     :: DeltaZ, DeltaY, DeltaX
!       integer(B4Ki)  :: NumTimeStepsHalf
!       real(SiKi)     :: MFFWS, zLu, yLu, xLu
!       integer(B4Ki)  :: dummy, rnd, nz, ny
!    end type

!    type :: Sub2Type
!       real(SiKi)     :: zLv, yLv, xLv, zLw, yLw, xLw
!    end type

!    type(Turb4Type)  :: Turb4
!    type(Turb78Type) :: Turb78
!    type(Sub1Type)   :: Sub1
!    type(Sub2Type)   :: Sub2
!    type(Turb7Type)  :: Turb7
!    type(Turb8Type)  :: Turb8

!    ! Read turbulence type (int16)
!    read (WindFileUnit, IOSTAT=ErrStat) TurbType
!    if (ErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading header 0 from binary FF file.', ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    ! Switch based on turbulence type
!    select case (TurbType)

!    case (1, 2)    ! 1-component Von Karman (1) or Kaimal (2)
!       p%NFFComp = 1

!    case (3, 5)    ! 3-component Von Karman (3) or IEC-2 Kaimal (5)
!       p%NFFComp = 3

!    case (4)       ! Improved Von Karman
!       read (WindFileUnit, IOSTAT=ErrStat) Turb4
!       if (ErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading header for TurbType=4 from binary FF file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if

!       ! Number of components (should be 3)
!       p%NFFComp = Turb4%NComp

!       ! Turbulence intensity
!       where (Turb4%TurbInt > 0) TI = Turb4%TurbInt

!    case (7, 8)       ! General Kaimal (7) or Mann model (8)
!       read (WindFileUnit, IOSTAT=ErrStat) Turb78
!       if (ErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading header for TurbType=7,8 from binary FF file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if

!       ! Number of components
!       p%NFFComp = Turb78%NComp

!    case default      ! Unknown turbulence type
!       call SetErrStat(ErrID_Warn, ' InflowWind does not recognize the full-field turbulence file type ='// &
!                       TRIM(Num2LStr(int(TurbType, IntKi)))//'.', ErrStat, ErrMsg, RoutineName)
!       if (ErrStat >= AbortErrLev) return
!    end select

!    ! Read header sub 1
!    read (WindFileUnit, IOSTAT=ErrStat) Sub1
!    if (ErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, ' Error reading header for Sub1 from binary FF file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    FFZDelt = Sub1%DeltaZ
!    p%InvFFZD = 1/FFZDelt
!    FFYDelt = Sub1%DeltaY
!    p%InvFFYD = 1/FFYDelt
!    FFXDelt = Sub1%DeltaX

!    if (.not. NativeFormat) p%MeanFFWS = Sub1%MFFWS
!    p%InvMFFWS = 1.0/p%MeanFFWS
!    p%FFDTime = FFXDelt/p%MeanFFWS
!    p%FFRate = 1.0/p%FFDTime

!    p%NZGrids = Sub1%nz
!    p%FFZHWid = 0.5*FFZDelt*(p%NZGrids - 1)
!    p%NYGrids = Sub1%ny
!    p%FFYHWid = 0.5*FFYDelt*(p%NYGrids - 1)

!    if (p%NFFComp == 3) then
!       read (WindFileUnit, IOSTAT=ErrStat) Sub2
!       if (ErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading header for Sub2 from binary FF file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if
!    end if

!    ! Read additional parameters based on turbulence type
!    select case (TurbType)
!    case (7)
!       read (WindFileUnit, IOSTAT=ErrStat) Turb7
!       if (ErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading header for Turb7 from binary FF file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if
!    case (8)
!       read (WindFileUnit, IOSTAT=ErrStat) Turb8
!       if (ErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading header for Turb8 from binary FF file.', &
!                          ErrStat, ErrMsg, RoutineName)
!          return
!       end if
!    end select

! end subroutine

! subroutine Read_Bladed_Grids(p, WindFileUnit, LHR, TI, NativeFormat, CWise, FFXDelt, FFYDelt, FFZDelt, ErrStat, ErrMsg)

!    type(IfW_Interp_ParameterType), intent(inout)   :: p              !< Parameters
!    integer(IntKi), intent(IN)                      :: WindFileUnit   !< unit number for the file to open
!    logical, intent(in)                             :: LHR
!    real(ReKi), intent(out)                         :: TI(3)          !< turbulence intensities
!    logical, intent(in)                             :: NativeFormat
!    logical, intent(in)                             :: Cwise
!    real(ReKi), intent(out)                         :: FFXDelt
!    real(ReKi), intent(out)                         :: FFYDelt
!    real(ReKi), intent(out)                         :: FFZDelt
!    integer(IntKi), intent(out)                     :: ErrStat        !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)                       :: ErrMsg         !< message about the error encountered

!    character(*), parameter                         :: RoutineName = "Read_Bladed_Grids"

!    real(ReKi)                                      :: FF_Scale(3)       !< used for "un-normalizing" the data
!    real(ReKi)                                      :: FF_Offset(3)      !< used for "un-normalizing" the data
!    integer(B2Ki), allocatable                      :: raw_ff(:, :, :)
!    integer(IntKi)                                  :: CFirst, CLast, CStep
!    integer(IntKi)                                  :: IC, IT
!    integer(IntKi)                                  :: TmpErrStat     ! temporary error status
!    character(ErrMsgLen)                            :: TmpErrMsg      ! temporary error message

!    ! Initialize error variables
!    ErrMsg = ""
!    ErrStat = ErrID_None

!    ! Calculate wind scale and offset
!    if (NativeFormat) then
!       FF_Scale = 0.001_ReKi
!       FF_Offset = 0.0_ReKi
!    else
!       FF_Scale = 0.001_ReKi*p%MeanFFWS*TI/100.0_ReKi
!       FF_Offset = (/p%MeanFFWS, 0.0_ReKi, 0.0_ReKi/)
!    end if

!    ! Bladed convention has positive V pointed along negative Y
!    if (LHR) FF_Scale(2) = -FF_Scale(2) ! left-hand rule

!    ! Print message to inform user that wind file is being read (may take a while).
!    call WrScr(NewLine//'   Reading a '//TRIM(Num2LStr(p%NYGrids))// &
!               'x'//TRIM(Num2LStr(p%NZGrids))// &
!               ' grid ('//TRIM(Num2LStr(p%FFYHWid*2))// &
!               ' m wide, '//TRIM(Num2LStr(p%GridBase))// &
!               ' m to '//TRIM(Num2LStr(p%GridBase + p%FFZHWid*2))// &
!               ' m above ground) with a characteristic wind speed of '// &
!               TRIM(Num2LStr(p%MeanFFWS))//' m/s. ')

!    !-------------------------------------------------------------------------
!    ! Allocate space for the FF array
!    !-------------------------------------------------------------------------

!    ! Add another step, just in case there is an odd number of steps.
!    p%NSteps = p%NSteps + 1

!    ! Allocate wind data
!    call Allocate_FF_Wind_Data(p%FFData, p%NSteps, p%NFFComp, &
!                               p%NYGrids, p%NZGrids, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    !-------------------------------------------------------------------------
!    ! Initialize the data and set column indexing to account for
!    ! direction of turbine rotation (CWise)
!    !-------------------------------------------------------------------------

!    ! Initialize entire array
!    p%FFData(:, :, :, :) = 0.0_SiKi     ! we may have only one component

!    if (CWise) then
!       CFirst = p%NYGrids
!       CLast = 1
!       CStep = -1
!    else
!       CFirst = 1
!       CLast = p%NYGrids
!       CStep = 1
!    end if

!    !----------------------------------------------------------------------------
!    ! Loop through all the time steps, reading the data and converting to m/s
!    !----------------------------------------------------------------------------

!    ! Allocate raw full-field array to hold grid data for single time step
!    allocate (raw_ff(p%NFFComp, p%NYGrids, p%NZGrids))

!    ! Loop through time steps
!    do IT = 1, p%NSteps

!       ! Read raw data (NFFComp,NYGrids,NZGrids)
!       read (WindFileUnit, IOStat=TmpErrStat) raw_ff

!       ! If there was an error
!       if (TmpErrStat /= 0) then

!          ! If there really were an even number of steps
!          if (IT == p%NSteps) then
!             p%NSteps = p%NSteps - 1
!             ErrStat = 0
!             exit
!          end if

!          ! Error reading file
!          call SetErrStat(ErrID_Fatal, ' Error reading binary data file. '// &
!                          ', it = '//TRIM(Num2LStr(it))// &
!                          ', nsteps = '//TRIM(Num2LStr(p%NSteps)), ErrStat, ErrMsg, RoutineName)
!          return
!       end if

!       ! Scale wind components and store in array
!       do IC = 1, p%NFFComp
!          p%FFData(IC, :, :, IT) = real(FF_Offset(IC) + FF_Scale(IC)*raw_ff(IC, CFirst:CLast:CStep, :), SiKi)
!       end do

!    end do

!    if (p%Periodic) then
!       call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(p%NSteps))//' time steps of '// &
!                  TRIM(Num2LStr(p%FFRate))//'-Hz full-field data (period of '// &
!                  TRIM(Num2LStr(p%FFDTime*p%NSteps))//' seconds).')

!    else
!       call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(p%NSteps))//' time steps of '// &
!                  TRIM(Num2LStr(p%FFRate))//'-Hz full-field data ('// &
!                  TRIM(Num2LStr(p%FFDTime*(p%NSteps - 1)))//' seconds).')
!    end if

! end subroutine Read_Bladed_Grids

! !> This subroutine reads the binary tower file that corresponds with the
! !! Bladed-style FF binary file. The FF grid must be read before this subroutine
! !! is called! (many checks are made to ensure the files belong together)
! subroutine Read_Bladed_Tower(p, ErrStat, ErrMsg)

!    type(IfW_Interp_ParameterType), intent(inout)   :: p              !< Parameters
!    integer(IntKi), intent(out)                     :: ErrStat        !< error status return value (0=no error; non-zero is error)
!    character(*), intent(out)                       :: ErrMsg         !< message about the error encountered

!    type :: HeaderType
!       real(SiKi)     :: dz, dx, Zmax
!       integer(B4Ki)  :: NumOutSteps, NumZ
!       real(SiKi)     :: UHub, TI(3)
!    end type

!    real(ReKi), parameter      :: TOL = 1E-4     ! tolerence for wind file comparisons
!    real(ReKi), parameter      :: FF_Offset(3) = (/1.0, 0.0, 0.0/)  ! used for "un-normalizing" the data

!    integer(IntKi)             :: IC, IT         ! loop counters
!    integer(B2Ki), allocatable :: raw_twr(:, :)   ! holds tower velocity for one timestep

!    type(HeaderType)           :: header

!    !-------------------------------------------------------------------------
!    ! Initialization
!    !-------------------------------------------------------------------------

!    ErrMsg = ''
!    ErrStat = ErrID_None

!    p%NTGrids = 0

!    if (p%NFFComp /= 3) then
!       call SetErrStat(ErrID_Fatal, ' Error: Tower binary files require 3 wind components.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    !-------------------------------------------------------------------------
!    ! Open the file
!    !-------------------------------------------------------------------------

!    call OpenBInpFile(WindFileUnit, TwrFile, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    !-------------------------------------------------------------------------
!    ! Read the header information and check that it's compatible with the FF Bladed-style binary
!    ! parameters already read.
!    !-------------------------------------------------------------------------

!    read (WindFileUnit, IOSTAT=TmpErrStat) header
!    if (TmpErrStat /= 0) then
!       call SetErrStat(ErrID_Fatal, &
!                       ' Error reading header of the binary tower file "' &
!                       //TRIM(TwrFile)//'."', ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    if (ABS(header%dz*p%InvFFZD - 1) > TOL) then
!       call SetErrStat(ErrID_Fatal, &
!                       ' Resolution in the FF binary file does not match the tower file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    if (ABS(header%dx*p%InvMFFWS/p%FFDTime - 1) > TOL) then
!       call SetErrStat(ErrID_Fatal, &
!                       ' Time resolution in the FF binary file does not match the tower file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    if (ABS(header%Zmax/p%GridBase - 1) > TOL) then
!       call SetErrStat(ErrID_Fatal, &
!                       ' Height in the FF binary file does not match the tower file "' &
!                       //TRIM(TwrFile)//'."', ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    if (header%NumOutSteps /= p%NSteps) then
!       call SetErrStat(ErrID_Fatal, &
!                       ' Number of time steps in the FF binary file does not match the tower file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       return
!    end if

!    p%NTGrids = header%NumZ

!    if (ABS(header%UHub*p%InvMFFWS - 1) > TOL) then
!       call SetErrStat(ErrID_Fatal, &
!                       ' Mean wind speed in the FF binary file does not match the tower file.', &
!                       ErrStat, ErrMsg, RoutineName)
!       p%NTGrids = 0
!       return
!    end if

!    ! If number of tower grids is zero, close file and return
!    if (p%NTGrids == 0) then
!       close (WindFileUnit)
!       return
!    end if

!    !-------------------------------------------------------------------------
!    ! Allocate arrays for the tower points
!    !-------------------------------------------------------------------------

!    call Allocate_FF_Tower_Data(p%Tower, p%NSteps, p%NFFComp, &
!                                p%NTGrids, TmpErrStat, TmpErrMsg)
!    call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
!    if (ErrStat >= AbortErrLev) return

!    !-------------------------------------------------------------------------
!    ! Read the 16-bit time-series data and scale it to 32-bit reals
!    !-------------------------------------------------------------------------

!    ! Loop through time.
!    do IT = 1, p%NSteps

!       ! Read normalized wind-component, INT(2) for this timestep
!       read (WindFileUnit, IOSTAT=TmpErrStat) raw_twr
!       if (TmpErrStat /= 0) then
!          call SetErrStat(ErrID_Fatal, ' Error reading binary tower data file. it = '//TRIM(Num2LStr(it))// &
!                          ', nsteps = '//TRIM(Num2LStr(p%NSteps)), ErrStat, ErrMsg, RoutineName)
!          p%NTGrids = 0
!          return
!       end if

!       ! Convert wind components to m/s
!       do IC = 1, p%NFFComp
!          p%Tower(IC, :, IT) = real(p%MeanFFWS*(FF_Offset(IC) + 0.00001*TI(IC)*raw_twr(IC, :)), SiKi)
!       end do
!    end do

!    !-------------------------------------------------------------------------
!    ! Close the file
!    !-------------------------------------------------------------------------

!    close (WindFileUnit)

!    TmpErrMsg = '   Processed '//TRIM(Num2LStr(p%NSteps))//' time steps of '// &
!                TRIM(Num2LStr(p%NTGrids))//'x1 tower data grids.'

!    call WrScr(NewLine//TRIM(TmpErrMsg))

! end subroutine Read_Bladed_Tower

end module FlowField_IO
