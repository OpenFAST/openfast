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

module InflowWind_IO

use NWTC_Library
use InflowWind_IO_Types
use IfW_FlowField

implicit none
private

public :: IfW_SteadyWind_Init, &
          IfW_UniformWind_Init, &
          IfW_TurbSim_Init, &
          IfW_Bladed_Init, &
          IfW_HAWC_Init, &
          IfW_User_Init, &
          IfW_Grid4D_Init, &
          IfW_Points_Init

public :: Uniform_WriteHH, &
          Grid3D_WriteBladed, &
          Grid3D_WriteHAWC, &
          Grid3D_WriteVTK

type(ProgDesc), parameter :: InflowWind_IO_Ver = ProgDesc('InflowWind_IO', '', '')

integer(IntKi), parameter :: ScaleMethod_None = 0, &           !< no scaling
                             ScaleMethod_Direct = 1, &         !< direct scaling factors
                             ScaleMethod_StdDev = 2            !< requested standard deviation

contains

subroutine IfW_Points_Init(InitInp, PF, ErrStat, ErrMsg)
   type(Points_InitInputType), intent(in) :: InitInp
   type(PointsFieldType), intent(out)     :: PF
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'IfW_Points_Init'
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! UVW components at points
   call AllocAry(PF%Vel, 3, InitInp%NumWindPoints, &
                 'Point Velocity Array', TmpErrStat, TmpErrMsg)
   call SetErrStat(ErrStat, ErrMsg, TmpErrStat, TmpErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

end subroutine

!> IfW_SteadyWind_Init initializes a Uniform field with with one set of values.
subroutine IfW_SteadyWind_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
   type(Steady_InitInputType), intent(in) :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(UniformFieldType), intent(out)    :: UF
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'IfW_SteadyWind_Init'
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Set parameters from inititialization input
   UF%DataSize = 1
   UF%RefHeight = InitInp%RefHt
   UF%RefLength = 1.0_ReKi

   ! Allocate uniform wind data arrays
   call UniformWind_AllocArrays(UF, TmpErrStat, TmpErrMsg)
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

!> IfW_UniformWind_Init initializes a Uniform field from file.
subroutine IfW_UniformWind_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)
   type(Uniform_InitInputType), intent(in)   :: InitInp
   integer(IntKi), intent(in)                :: SumFileUnit
   type(UniformFieldType), intent(out)       :: UF
   type(WindFileDat), intent(out)            :: FileDat
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter                   :: RoutineName = 'IfW_UniformWind_Init'
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
   call UniformWind_AllocArrays(UF, TmpErrStat, TmpErrMsg)
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
      call SetErrStat(ErrID_Info, ' Could not read upflow column in uniform wind files. Assuming upflow is 0.', &
                      ErrStat, ErrMsg, RoutineName)
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

   ! Check if the fist data point from the file is not along the X-axis while
   ! applying the windfield rotation
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
      write (SumFileUnit, '(A)') 'Uniform wind.  Module '//TRIM(InflowWind_IO_Ver%Name)//' '//TRIM(InflowWind_IO_Ver%Ver)
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

!> UniformWind_AllocArrays allocates the data arrays in the Uniform field.
subroutine UniformWind_AllocArrays(UF, ErrStat, ErrMsg)
   type(UniformFieldType), intent(inout)     :: UF
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter                   :: RoutineName = 'UniformWind_AllocArrays'
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

!> Uniform_WriteHH writes a Uniform field hub-height wind file.
subroutine Uniform_WriteHH(UF, FileRootName, unit, ErrStat, ErrMsg)

   type(UniformFieldType), intent(in)  :: UF             !< Parameter
   character(*), intent(in)            :: FileRootName   !< RootName for output files
   integer(IntKi), intent(in)          :: Unit           !< Indicates whether an error occurred (see NWTC_Library)
   integer(IntKi), intent(out)         :: ErrStat        !< Error status of the operation
   character(*), intent(out)           :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter             :: RoutineName = 'Uniform_WriteHH'
   integer(IntKi)                      :: i
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   call OpenFOutFile(unit, trim(FileRootName)//'.UniformWind.dat', ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   write (unit, "(A)") "#"
   write (unit, "(A)") '# Uniform Wind (deterministic) file for ENFAST generated by InflowWind'
   write (unit, "(A)") "#"
   write (unit, "(A)") "!    Time             Wind        Wind       Vertical    Horiz.    Pwr.Law     Lin.Vert.     Gust        Upflow"
   write (unit, "(A)") "!                     Speed       Dir        Speed       Shear     Vert.Shr    Shear         Speed       Angle"
   write (unit, "(A)") "!    (sec)            (m/s)       (Deg)      (m/s)                                           (m/s)       (deg)"

   do i = 1, UF%DataSize
      write (unit, "(F15.5,8(1x,F11.4))") UF%Time(i), UF%VelH(i), UF%AngleH(i)*R2D, UF%VelV(i), &
         UF%ShrH(i), UF%ShrV(i), UF%LinShrV(i), UF%VelGust(i), UF%AngleV(i)*R2D
   end do

   close (unit)

end subroutine Uniform_WriteHH

!> Read_TurbSim reads the binary TurbSim-format FF file (.bts).  It fills the FFData array with
!! velocity data for the grids and fills the Tower array with velocities at points on the tower
!! (if data exists).
subroutine IfW_TurbSim_Init(InitInp, SumFileUnit, G3D, FileDat, ErrStat, ErrMsg)

   type(TurbSim_InitInputType), intent(in)    :: InitInp
   integer(IntKi), intent(in)                :: SumFileUnit
   type(Grid3DFieldType), intent(out)        :: G3D
   type(WindFileDat), intent(out)            :: FileDat
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   character(*), parameter       :: RoutineName = "IfW_TurbSim_Init"
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

   type :: TurbSimHeaderType
      sequence
      integer(B2Ki)  :: FileID
      integer(B4Ki)  :: NZGrids, NYGrids, NTGrids, NSteps
      real(SiKi)     :: dz, dy, dt
      real(SiKi)     :: mws, ref_height, grid_base_height
      real(SiKi)     :: VslopeX, VoffsetX
      real(SiKi)     :: VslopeY, VoffsetY
      real(SiKi)     :: VslopeZ, VoffsetZ
      integer(B4Ki)  :: DescLen
   end type

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

   G3D%WindFileFormat = header%FileID     ! file format identifier
   G3D%Periodic = header%FileID == 8      ! 7 is used for non-periodic wind files; 8 is periodic wind
   G3D%InterpTower = .false.              ! wind should not be interpolated at tower
   G3D%AddMeanAfterInterp = .false.       ! do not add mean wind speed after interpolation

   G3D%DTime = real(header%dt, ReKi)      ! grid spacing in time (dt), m/s
   G3D%Rate = 1.0_ReKi/G3D%DTime          ! Data rate (1/DTime), Hertz

   G3D%NComp = 3                          ! TurbSim file file contains 3 wind components
   G3D%NYGrids = header%NYGrids           ! the number of grid points laterally
   G3D%NZGrids = header%NZGrids           ! the number of grid points vertically
   G3D%NTGrids = header%NTGrids           ! the number of tower points
   G3D%NSteps = header%NSteps             ! the number of time steps

   G3D%InvDY = 1.0_ReKi/real(header%dy, ReKi)     ! 1/dy
   G3D%YHWid = 0.5_ReKi*(G3D%NYGrids - 1)/G3D%InvDY   ! half the grid width (m)

   G3D%InvDZ = 1.0_ReKi/real(header%dz, ReKi)     ! 1/dz
   G3D%ZHWid = 0.5_ReKi*(G3D%NZGrids - 1)/G3D%InvDZ   ! half the grid height (m)

   G3D%MeanWS = real(header%mws, ReKi)                 ! the mean wind speed at hub height (m/s)
   G3D%InvMWS = 1.0_ReKi/G3D%MeanWS                      ! inverse of mean wind speed

   G3D%RefHeight = real(header%ref_height, ReKi)       ! height of the hub (m)
   G3D%GridBase = real(header%grid_base_height, ReKi)  ! height of the bottom of the grid (m)

   if (G3D%Periodic) then
      G3D%InitXPosition = 0                       ! start at the hub
      G3D%TotalTime = G3D%NSteps*G3D%DTime
   else
      G3D%InitXPosition = G3D%YHWid                 ! start half the grid width ahead of the turbine
      G3D%TotalTime = (G3D%NSteps - 1)*G3D%DTime
   end if

   G3D%WindProfileType = WindProfileType_None     ! unused for turbsim
   G3D%PLExp = 0                                  ! unused for turbsim
   G3D%Z0 = 0                                     ! unused for turbsim

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
   call AllocAry(G3D%Vel, G3D%NComp, G3D%NYGrids, G3D%NZGrids, G3D%NSteps, &
                 'grid-field velocity data', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for raw grid-field velocity for each time step
   allocate (VelRaw(G3D%NComp, G3D%NYGrids, G3D%NZGrids), stat=TmpErrStat)
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, "error allocating grid-field time step velocity data", &
                      ErrStat, ErrMsg, RoutineName)
   end if
   if (ErrStat >= AbortErrLev) return

   ! If tower grids specified
   if (G3D%NTGrids > 0) then

      ! Allocate storage for tower velocity data
      call AllocAry(G3D%VelTower, G3D%NComp, G3D%NTGrids, G3D%NSteps, &
                    'tower wind velocity data.', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Allocate storage for raw tower data for each timestep
      allocate (TwrRaw(G3D%NComp, G3D%NTGrids), stat=TmpErrStat)
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
              //TRIM(Num2LStr(G3D%NYGrids))//'x' &
              //TRIM(Num2LStr(G3D%NZGrids))// &
              ' grid ('//TRIM(Num2LStr(G3D%YHWid*2))//' m wide, '// &
              TRIM(Num2LStr(G3D%GridBase))//' m to '// &
              TRIM(Num2LStr(G3D%GridBase + G3D%ZHWid*2))// &
              ' m above ground) with a characteristic wind speed of '// &
              TRIM(Num2LStr(G3D%MeanWS))//' m/s. '//TRIM(DescStr))

   ! Loop through time steps
   do IT = 1, G3D%NSteps

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
         G3D%Vel(IC, :, :, IT) = (real(VelRaw(IC, :, :), SiKi) - Voffset(IC))/VSlope(IC)
      end do !IC

      ! Read tower raw wind data (normalized) comprised of 2-byte integers, INT(2)
      ! Indices are Velocity components, Z coordinates
      if (G3D%NTGrids > 0) then
         read (WindFileUnit, IOSTAT=TmpErrStat) TwrRaw
         if (TmpErrStat /= 0) then
            call SetErrStat(ErrID_Fatal, ' Error reading tower wind components in the FF binary file "'// &
                            TRIM(InitInp%WindFileName)//'."', ErrStat, ErrMsg, RoutineName)
            return
         end if

         ! Loop through wind components (U, V, W), calculate de-normalized velocity (m/s)
         do IC = 1, 3
            G3D%VelTower(IC, :, IT) = (real(TwrRaw(IC, :), SiKi) - Voffset(IC))/VSlope(IC)
         end do
      end if
   end do

   !----------------------------------------------------------------------------
   ! Close the file
   !----------------------------------------------------------------------------

   close (WindFileUnit)

   if (G3D%Periodic) then
      call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(G3D%NSteps))//' time steps of '// &
                 TRIM(Num2LStr(G3D%Rate))//'-Hz grid-field data (period of '// &
                 TRIM(Num2LStr(G3D%DTime*(G3D%NSteps)))//' seconds).')
   else
      call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(G3D%NSteps))//' time steps of '// &
                 TRIM(Num2LStr(G3D%Rate))//'-Hz grid-field data ('// &
                 TRIM(Num2LStr(G3D%DTime*(G3D%NSteps - 1)))//' seconds).')
   end if

   !----------------------------------------------------------------------------
   ! Store wind file metadata
   !----------------------------------------------------------------------------

   call Grid3D_PopulateWindFileDat(G3D, InitInp%WindFileName, 3, G3D%NTGrids > 0, FileDat)

   !----------------------------------------------------------------------------
   ! Write the summary file
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A)') 'TurbSim wind type.  Read by InflowWind sub-module ' &
         //TRIM(InflowWind_IO_Ver%Name)//' '//TRIM(InflowWind_IO_Ver%Ver)
      write (SumFileUnit, '(A)') TRIM(TmpErrMsg)
      write (SumFileUnit, '(5x,A)') 'FileName:                    '//TRIM(InitInp%WindFileName)
      write (SumFileUnit, '(5x,A29,I3)') 'Binary file format id:       ', G3D%WindFileFormat
      write (SumFileUnit, '(5x,A29,G12.4)') 'Reference height (m):        ', G3D%RefHeight
      write (SumFileUnit, '(5x,A29,G12.4)') 'Timestep (s):                ', G3D%DTime
      write (SumFileUnit, '(5x,A29,I12)') 'Number of timesteps:         ', G3D%NSteps
      write (SumFileUnit, '(5x,A29,G12.4)') 'Mean windspeed (m/s):        ', G3D%MeanWS
      write (SumFileUnit, '(5x,A29,L1)') 'Windfile is periodic:        ', G3D%Periodic
      write (SumFileUnit, '(5x,A29,L1)') 'Windfile includes tower:     ', G3D%NTGrids > 0

      if (G3D%Periodic) then
         write (SumFileUnit, '(5x,A)') 'Time range (s):              [ '// &
            TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(G3D%TotalTime))//' ]'
      else  ! Shift the time range to compensate for the shifting of the wind grid
         write (SumFileUnit, '(5x,A)') 'Time range (s):              [ '// &
            TRIM(Num2LStr(-G3D%InitXPosition*G3D%InvMWS))//' : '// &
            TRIM(Num2LStr(G3D%TotalTime - G3D%InitXPosition*G3D%InvMWS))//' ]'
      end if

      write (SumFileUnit, '(5x,A)') 'Y range (m):                 [ '// &
         TRIM(Num2LStr(-G3D%YHWid))//' : '//TRIM(Num2LStr(G3D%YHWid))//' ]'

      if (G3D%NTGrids > 0) then
         write (SumFileUnit, '(5x,A)') 'Z range (m):                 [ '// &
            TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(G3D%RefHeight + G3D%ZHWid))//' ]'
      else
         write (SumFileUnit, '(5x,A)') 'Z range (m):                 [ '// &
            TRIM(Num2LStr(G3D%RefHeight - G3D%ZHWid))//' : '//TRIM(Num2LStr(G3D%RefHeight + G3D%ZHWid))//' ]'
      end if

      if (G3D%BoxExceedAllowF) then
         write (SumFileUnit, '(A)') '     Wind grid exceedence allowed:  '// &
            'True      -- Only for points requested by OLAF free vortex wake, or LidarSim module'
         write (SumFileUnit, '(A)') '                                    '// &
            '             Out of bounds values are linearly interpolated to mean at Z loction for'
         write (SumFileUnit, '(A)') '                                    '// &
            '             given timestep and X,T value. Values above grid are held to top of wind'
         write (SumFileUnit, '(A)') '                                    '// &
            '             grid value'
      else
         write (SumFileUnit, '(A)') '     Wind grid exceedence allowed:  False'
      end if

      ! Get IO status for unit
      inquire (SumFileUnit, iostat=TmpErrStat)
      if (TmpErrStat /= 0_IntKi) then
         call SetErrStat(ErrID_Fatal, 'Error writing to summary file.', ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

end subroutine

subroutine IfW_HAWC_Init(InitInp, SumFileUnit, G3D, FileDat, ErrStat, ErrMsg)

   type(HAWC_InitInputType), intent(in)   :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(Grid3DFieldType), intent(out)     :: G3D
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter       :: RoutineName = "IfW_HAWC_Init"
   integer(IntKi)                :: WindFileUnit
   real(SiKi), allocatable       :: VelRaw(:, :)      ! grid-field data for one timestep
   integer                       :: IC                ! Loop counter for the number of wind components
   integer                       :: IX, IY, IZ        ! Loop counters for the number of grid points in the X,Y,Z directions
   real(DbKi)                    :: vMean             ! average wind speeds over time at target position
   real(ReKi)                    :: ScaleFactors(3)   ! scale factors
   integer(IntKi)                :: TmpErrStat        ! temporary error status
   character(ErrMsgLen)          :: TmpErrMsg

   !----------------------------------------------------------------------------
   ! Initialize variables
   !----------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg = ""

   G3D%WindFileFormat = 0
   G3D%Periodic = .true.
   G3D%InterpTower = .true.
   G3D%AddMeanAfterInterp = .true.

   G3D%DTime = InitInp%dx/InitInp%G3D%URef
   G3D%Rate = 1.0_ReKi/G3D%DTime

   G3D%NComp = 3
   G3D%NYGrids = InitInp%ny
   G3D%NZGrids = InitInp%nz
   G3D%NTGrids = 0
   G3D%NSteps = InitInp%nx

   G3D%YHWid = 0.5_ReKi*InitInp%dy*(G3D%NYGrids - 1)
   G3D%InvDY = 1.0/InitInp%dy

   G3D%ZHWid = 0.5_ReKi*InitInp%dz*(G3D%NZGrids - 1)
   G3D%InvDZ = 1.0_ReKi/InitInp%dz

   G3D%MeanWS = InitInp%G3D%URef
   G3D%InvMWS = 1.0_ReKi/G3D%MeanWS

   G3D%RefHeight = InitInp%G3D%RefHt
   G3D%GridBase = G3D%RefHeight - G3D%ZHWid

   G3D%InitXPosition = InitInp%G3D%XOffset
   G3D%TotalTime = G3D%NSteps*InitInp%dx/G3D%MeanWS

   G3D%WindProfileType = InitInp%G3D%WindProfileType
   G3D%Z0 = InitInp%G3D%Z0
   G3D%PLExp = InitInp%G3D%PLExp

   ScaleFactors = 0.0_ReKi

   !----------------------------------------------------------------------------
   ! Validate initialization iput
   !----------------------------------------------------------------------------

   call Grid3D_ValidateInput(InitInp%G3D, 3, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   if (InitInp%nx < 1) then
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
   end if
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Allocate storage for grid-field velocity data
   !----------------------------------------------------------------------------

   ! Allocate storage for grid-field velocity data
   call AllocAry(G3D%Vel, G3D%NComp, G3D%NYGrids, G3D%NZGrids, G3D%NSteps, &
                 'grid-field velocity data', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for raw grid-field velocity for each time step
   allocate (VelRaw(G3D%NZGrids, G3D%NYGrids), stat=TmpErrStat)
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
              TRIM(Num2LStr(G3D%NSteps))//' x '//TRIM(Num2LStr(G3D%NYGrids))//' x '//TRIM(Num2LStr(G3D%NZGrids))//' points'// &
              ' ('//TRIM(Num2LStr(G3D%YHWid*2))//' m wide, '//TRIM(Num2LStr(G3D%GridBase))//' m to '// &
              TRIM(Num2LStr(G3D%GridBase + G3D%ZHWid*2))// &
              ' m above ground) with a characteristic wind speed of '//TRIM(Num2LStr(G3D%MeanWS))//' m/s. ')

   ! Get a unit number to use for the wind file
   call GetNewUnit(WindFileUnit, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Loop through wind components (X, Y, Z)
   do IC = 1, G3D%NComp

      ! Open wind file for this component
      call OpenBInpFile(WindFileUnit, InitInp%WindFileName(IC), TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Loop through time steps
      do IX = 1, G3D%NSteps

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
         do IZ = 1, G3D%NZGrids
            G3D%Vel(IC, :, IZ, IX) = VelRaw(IZ, G3D%NYGrids:1:-1)
         end do
      end do

      ! Close file
      close (WindFileUnit)

   end do

   !----------------------------------------------------------------------------
   ! Scale turbulence to requested intensity
   !----------------------------------------------------------------------------

   call Grid3D_ScaleTurbulence(InitInp%G3D, G3D%Vel, ScaleFactors, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Remove the U component mean wind speed
   !----------------------------------------------------------------------------

   ! If scaling method is not none, remove mean value of X component at each grid point
   if (InitInp%G3D%ScaleMethod /= ScaleMethod_None) then
      do iz = 1, G3D%NZGrids
         do iy = 1, G3D%NYGrids
            vMean = sum(G3D%Vel(1, iy, iz, :))/G3D%NSteps
            G3D%Vel(1, iy, iz, :) = real(G3D%Vel(1, iy, iz, :) - vMean, SiKi)
         end do
      end do
   end if

   !----------------------------------------------------------------------------
   ! Store wind file metadata
   !----------------------------------------------------------------------------

   call Grid3D_PopulateWindFileDat(G3D, InitInp%WindFileName(1), 5, .false., FileDat)

   !----------------------------------------------------------------------------
   ! Write the summary file
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A)') 'HAWC wind type.  Read by InflowWind sub-module InflowWind_IO'

      write (SumFileUnit, '(A34,G12.4)') '     Reference height (m):        ', G3D%RefHeight
      write (SumFileUnit, '(A34,G12.4)') '     Timestep (s):                ', G3D%DTime
      write (SumFileUnit, '(A34,I12)') '     Number of timesteps:         ', G3D%NSteps
      write (SumFileUnit, '(A34,G12.4)') '     Mean windspeed (m/s):        ', G3D%MeanWS
      write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
         TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(G3D%TotalTime))//' ]'
      write (SumFileUnit, '(A)') '     X range (m):                 [ '// &
         TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(G3D%TotalTime*G3D%MeanWS))//' ]'
      write (SumFileUnit, '(A)') '     Y range (m):                 [ '// &
         TRIM(Num2LStr(-G3D%YHWid))//' : '//TRIM(Num2LStr(G3D%YHWid))//' ]'
      write (SumFileUnit, '(A)') '     Z range (m):                 [ '// &
         TRIM(Num2LStr(G3D%GridBase))//' : '//TRIM(Num2LStr(G3D%GridBase + G3D%ZHWid*2.0))//' ]'

      if (G3D%BoxExceedAllowF) then
         write (SumFileUnit, '(A)') '     Wind grid exceedence allowed:  '// &
            'True      -- Only for points requested by OLAF free vortex wake, or LidarSim module'
         write (SumFileUnit, '(A)') '                                    '// &
            '             Out of bounds values are linearly interpolated to mean at Z loction for'
         write (SumFileUnit, '(A)') '                                    '// &
            '             given timestep and X,T value. Values above grid are held to top of wind'
         write (SumFileUnit, '(A)') '                                    '// &
            '             grid value'
      else
         write (SumFileUnit, '(A)') '     Wind grid exceedence allowed:  False'
      end if

      write (SumFileUnit, '(A)') 'Scaling factors used:'
      write (SumFileUnit, '(A)') '  u           v           w       '
      write (SumFileUnit, '(A)') '----------  ----------  ----------'
      write (SumFileUnit, '(F10.3,2x,F10.3,2x,F10.3)') ScaleFactors
   end if

end subroutine

!> User_Init initializes a user defined wind field.
subroutine IfW_User_Init(InitInp, SumFileUnit, UF, FileDat, ErrStat, ErrMsg)

   type(User_InitInputType), intent(in)    :: InitInp
   integer(IntKi), intent(in)             :: SumFileUnit
   type(UserFieldType), intent(out)       :: UF
   type(WindFileDat), intent(out)         :: FileDat
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = "IfW_User_Init"

   ErrStat = ErrID_None
   ErrMsg = ""

   UF%RefHeight = 0.0_Reki

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)') UF%RefHeight
   end if

end subroutine

!> IfW_Grid4D_Init initializes a wind field defined by a 4D grid.
subroutine IfW_Grid4D_Init(InitInp, G4D, ErrStat, ErrMsg)

   type(Grid4D_InitInputType), intent(in) :: InitInp
   type(Grid4DFieldType), intent(out)     :: G4D
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = "IfW_Grid4D_Init"
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Initialize field from inputs
   G4D%n = InitInp%n
   G4D%delta = InitInp%delta
   G4D%pZero = InitInp%pZero
   G4D%TimeStart = 0.0_ReKi
   G4D%RefHeight = InitInp%pZero(3) + (InitInp%n(3)/2) * InitInp%delta(3)

   ! uvw velocity components at x,y,z,t coordinates
   call AllocAry(G4D%Vel, 3, G4D%n(1), G4D%n(2), G4D%n(3), G4D%n(4), &
                 'External Grid Velocity', TmpErrStat, TmpErrMsg)
   call SetErrStat(ErrStat, ErrMsg, TmpErrStat, TmpErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Initialize velocities to zero
   G4D%Vel = 0.0_SiKi

end subroutine

subroutine IfW_Bladed_Init(InitInp, SumFileUnit, InitOut, G3D, FileDat, ErrStat, ErrMsg)

   type(Bladed_InitInputType), intent(in)    :: InitInp  !< Initialization data passed to the module
   integer(IntKi), intent(in)                :: SumFileUnit
   type(Bladed_InitOutputType), intent(out)  :: InitOut  !< Initial output
   type(Grid3DFieldType), intent(out)        :: G3D      !< Parameters
   type(WindFileDat), intent(out)            :: FileDat
   integer(IntKi), intent(out)               :: ErrStat  !< determines if an error has been encountered
   character(*), intent(out)                 :: ErrMsg   !< Message about errors

   character(*), parameter    :: RoutineName = "IfW_Bladed_Init"
   real(ReKi)                 :: TI(3)             ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   type(Grid3D_InitInputType) :: G3D_InitInp       ! initialization input for grid 3d field
   real(ReKi)                 :: BinTI(3)          ! turbulence intensities of the wind components as defined in the FF binary file, not necessarially the actual TI
   real(ReKi)                 :: NatTI(3)          ! turbulence intensities of the wind components as defined in the native FF summary file
   real(ReKi)                 :: UBar
   real(ReKi)                 :: ZCenter
   real(ReKi)                 :: ScaleFactors(3)   ! turbulence scaling factors
   integer(IntKi)             :: UnitWind          ! Unit number for the InflowWind input file
   integer(B2Ki)              :: Dum_Int2
   integer(IntKi)             :: I
   logical                    :: CWise
   logical                    :: LHR               ! Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
   logical                    :: Exists
   character(1028)            :: SumFile           ! length is LEN(ParamData%WindFileName) + the 4-character extension.
   character(1028)            :: TwrFile           ! length is LEN(ParamData%WindFileName) + the 4-character extension.
   character(1024)            :: BinFileName
   character(1024)            :: PriPath
   character(ErrMsgLen)       :: TmpErrMsg         ! temporary error message
   integer(IntKi)             :: TmpErrStat        ! temporary error status

   ErrMsg = ''
   ErrStat = ErrID_None

   if (InitInp%NativeBladedFmt) then

      call Bladed_ReadNativeSummary(InitInp%WindFileName, G3D_InitInp%PLExp, G3D_InitInp%VLinShr, &
                                    G3D_InitInp%HLinShr, G3D_InitInp%RefLength, NatTI, G3D%MeanWS, &
                                    G3D%RefHeight, InitOut%PropagationDir, InitOut%VFlowAngle, &
                                    BinFileName, G3D_InitInp%XOffset, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      if (pathIsRelative(BinFileName)) then
         call GetPath(InitInp%WindFileName, PriPath)     ! Binary file will be relative to the path where the primary input file is located.
         BinFileName = TRIM(PriPath)//TRIM(BinFileName)
      end if

      if (InitInp%FixedWindFileRootName) then ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
         if (InitInp%TurbineID == 0) then     ! .TRUE. for the FAST.Farm low-resolution domain
            BinFileName = TRIM(BinFileName)//TRIM(PathSep)//'Low'
         else                                   ! FAST.Farm high-resolution domain(s)
            BinFileName = TRIM(BinFileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
         end if
      end if

      ! default values for Bladed Format
      CWise = .false.
      ZCenter = G3D%RefHeight
      G3D%Periodic = .true.

      G3D_InitInp%ScaleMethod = ScaleMethod_StdDev
      G3D_InitInp%SigmaF = NatTI*G3D%MeanWS
      G3D_InitInp%SF = G3D_InitInp%SigmaF

      ! it could also have logarithmic, but I'm going to leave that off for now
      G3D_InitInp%RefHt = G3D%RefHeight
      G3D_InitInp%URef = G3D%MeanWS
      G3D_InitInp%WindProfileType = WindProfileType_PL ! it could also have logarithmic, but I'm going to leave that off for now

      TI = 100.0_ReKi
      UBar = 0.0_ReKi
      LHR = .true.

   else
      InitOut%PropagationDir = 0.0_ReKi
      InitOut%VFlowAngle = 0.0_ReKi
      G3D%VLinShr = 0.0_ReKi
      G3D%HLinShr = 0.0_ReKi
      G3D%RefLength = 1.0_ReKi

      BinFileName = InitInp%WindFileName
   end if

   ! Get a unit number to use
   call GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Open the binary file, read its "header" (first 2-byte integer) to
   ! determine what format binary file it is, and close it.
   !----------------------------------------------------------------------------

   call OpenBInpFile(UnitWind, TRIM(BinFileName), TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Read the first binary integer from the file to get info on the type.
   ! Cannot use library read routines since this is a 2-byte integer.
   read (UnitWind, IOSTAT=TmpErrStat) Dum_Int2
   close (UnitWind)

   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading first binary integer from file "' &
                      //TRIM(BinFileName)//'."', ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Read the files to get the required FF data.
   !----------------------------------------------------------------------------

   ! Store the binary format information so the InflowWind code can use it.
   ! Also changes to IntKi from INT(2) to compare in the SELECT below
   G3D%WindFileFormat = Dum_Int2

   select case (G3D%WindFileFormat)

   case (-1, -2, -3, -99)                          ! Bladed-style binary format

      if (.not. InitInp%NativeBladedFmt) then

         !----------------------------------------------------------------------
         ! Create full-field summary file name from binary file root name.
         ! Also get tower file name.
         !----------------------------------------------------------------------

         call GetRoot(BinFileName, SumFile)

         TwrFile = TRIM(SumFile)//'.twr'
         SumFile = TRIM(SumFile)//'.sum'

         !----------------------------------------------------------------------
         ! Read the summary file to get necessary scaling information
         !----------------------------------------------------------------------

         call Bladed_ReadTurbSimSummary(UnitWind, TRIM(SumFile), CWise, ZCenter, TI, &
                                        UBar, G3D%RefHeight, G3D%Periodic, LHR, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            close (UnitWind)
            return
         end if

      end if

      !-------------------------------------------------------------------------
      ! Open the binary file and read its header
      !-------------------------------------------------------------------------

      call OpenBInpFile(UnitWind, TRIM(BinFileName), TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         close (UnitWind)
         return
      end if

      if (Dum_Int2 == -99) then                      ! Newer-style BLADED format
         call Bladed_ReadHeader1(UnitWind, BinTI, G3D, InitInp%NativeBladedFmt, TmpErrStat, TmpErrMsg)
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            close (UnitWind)
            return
         end if

         ! If the TIs are also in the binary file (BinTI > 0),
         ! use those numbers instead of ones from the summary file

         if (.not. InitInp%NativeBladedFmt) then
            do I = 1, G3D%NComp
               if (BinTI(I) > 0) TI(I) = BinTI(I)
            end do
         end if

      else
         call Bladed_ReadHeader0(UnitWind, G3D, InitInp%NativeBladedFmt, TmpErrStat, TmpErrMsg)     ! Older-style BLADED format
         call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            close (UnitWind)
            return
         end if

      end if

      !-------------------------------------------------------------------------
      ! Let's see if the summary and binary FF wind files go together before continuing.
      !-------------------------------------------------------------------------

      if (.not. InitInp%NativeBladedFmt) then
         if (ABS(UBar - G3D%MeanWS) > 0.1) then
            call SetErrStat(ErrID_Fatal, ' Error: Incompatible mean hub-height wind speeds in FF wind files. '// &
                            '(Check that the .sum and .wnd files were generated together.)', ErrStat, ErrMsg, RoutineName)
            close (UnitWind)
            return
         end if

      end if

      !-------------------------------------------------------------------------
      ! Calculate the height of the bottom of the grid
      !-------------------------------------------------------------------------

      G3D%GridBase = ZCenter - G3D%ZHWid         ! the location, in meters, of the bottom of the grid
      if (G3D%GridBase < 0.0_ReKi) then
         call SetErrStat(ErrID_Severe, 'WARNING: The bottom of the grid is located at a height of '// &
                         TRIM(Num2LStr(G3D%GridBase))//' meters, which is below the ground.'// &
                         ' Winds below the ground will be set to 0.', ErrStat, ErrMsg, RoutineName)
      end if

      !-------------------------------------------------------------------------
      ! Read the binary grids (converted tƒo m/s) and close the file
      !-------------------------------------------------------------------------

      call Bladed_ReadGrids(UnitWind, InitInp%NativeBladedFmt, CWise, LHR, TI, G3D, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

      close (UnitWind)
      if (InitInp%NativeBladedFmt) TI = NatTI*100.0_ReKi  ! report these TI for the native Bladed format in percent

      if (ErrStat >= AbortErrLev) return

      !-------------------------------------------------------------------------
      ! Read the tower points file
      !-------------------------------------------------------------------------

      if (InitInp%TowerFileExist .and. .not. InitInp%NativeBladedFmt) then      ! If we specified a tower file
         inquire (FILE=TRIM(TwrFile), EXIST=Exists)

         ! Double check that the tower file exists and read it.  If it was requested but doesn't exist,
         ! throw fatal error and exit.
         if (Exists) then
            call Bladed_ReadTower(UnitWind, G3D, TwrFile, TmpErrStat, TmpErrMsg)
            call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) then
               close (UnitWind)
               return
            end if
         else
            call SetErrStat(ErrID_Fatal, ' Tower file '//TRIM(TwrFile)//' specified for Bladed full-field '// &
                            'wind files does not exist.', ErrStat, ErrMsg, RoutineName)
            close (UnitWind)
            return
         end if
      else
         G3D%NTGrids = 0_IntKi
      end if

   case DEFAULT
      call SetErrStat(ErrID_Fatal, ' This is not a bladed-style binary wind file (binary format identifier: '// &
                      TRIM(Num2LStr(G3D%WindFileFormat))//'.  This might be a TurbSim binary wind file.', &
                      ErrStat, ErrMsg, RoutineName)
      return

   end select

   !----------------------------------------------------------------------------
   ! Post reading
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! If the wind file has zero-mean and unit standard deviation (native Bladed format), scale the data:
   !----------------------------------------------------------------------------

   G3D%AddMeanAfterInterp = .false.
   G3D%Z0 = G3D_InitInp%Z0
   G3D%PLExp = G3D_InitInp%PLExp
   G3D%VLinShr = G3D_InitInp%VLinShr
   G3D%HLinShr = G3D_InitInp%HLinShr
   G3D%RefLength = G3D_InitInp%RefLength

   if (InitInp%NativeBladedFmt) then

      G3D%InterpTower = .true.
      G3D%AddMeanAfterInterp = .true.
      G3D%WindProfileType = G3D_InitInp%WindProfileType

      ! Validate scaling data if we've got native-Bladed format
      call Grid3D_ValidateInput(G3D_InitInp, G3D%NComp, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! scale to requested TI (or use requested scale factors)
      call Grid3D_ScaleTurbulence(G3D_InitInp, G3D%Vel(:, :, :, 1:G3D%NSteps), ScaleFactors, TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   else
      G3D%InterpTower = .false.
      G3D%WindProfileType = WindProfileType_None
   end if

   if (G3D%Periodic) then
      G3D%InitXPosition = 0                ! start at the hub
      G3D%TotalTime = G3D%NSteps*G3D%DTime
   else
      G3D%InitXPosition = G3D%YHWid          ! start half the grid width ahead of the turbine
      G3D%TotalTime = (G3D%NSteps - 1)*G3D%DTime
   end if

   if (InitInp%NativeBladedFmt) then
      G3D%InitXPosition = G3D_InitInp%XOffset
   end if

   !----------------------------------------------------------------------------
   ! Store wind file metadata
   !----------------------------------------------------------------------------

   call Grid3D_PopulateWindFileDat(G3D, InitInp%WindFileName, InitInp%WindType, G3D%NTGrids > 0, FileDat)

   FileDat%TI = TI
   FileDat%TI_listed = .true.

   !----------------------------------------------------------------------------
   ! Write to the summary file
   !----------------------------------------------------------------------------

   if (SumFileUnit > 0) then
      write (SumFileUnit, '(A)')
      write (SumFileUnit, '(A)') 'Bladed-style wind type.  Read by InflowWind sub-module '// &
         TRIM(InflowWind_IO_Ver%Name)//' '//TRIM(InflowWind_IO_Ver%Ver)
      write (SumFileUnit, '(A)') TRIM(TmpErrMsg)
      write (SumFileUnit, '(A)') '     FileName:                    '//TRIM(InitInp%WindFileName)
      write (SumFileUnit, '(A34,I3)') '     Binary file format id:       ', G3D%WindFileFormat
      write (SumFileUnit, '(A34,G12.4)') '     Reference height (m):        ', G3D%RefHeight
      write (SumFileUnit, '(A34,G12.4)') '     Timestep (s):                ', G3D%DTime
      write (SumFileUnit, '(A34,I12)') '     Number of timesteps:         ', G3D%NSteps
      write (SumFileUnit, '(A34,G12.4)') '     Mean windspeed (m/s):        ', G3D%MeanWS
      write (SumFileUnit, '(A)') '     Characteristic TI:            [ '// &
         TRIM(Num2LStr(TI(1)))//', '//TRIM(Num2LStr(TI(2)))//', '//TRIM(Num2LStr(TI(3)))//' ] '
      write (SumFileUnit, '(A34,L1)') '     Windfile is periodic:        ', G3D%Periodic
      write (SumFileUnit, '(A34,L1)') '     Windfile includes tower:     ', G3D%NTGrids > 0

      if (G3D%Periodic) then
         write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
            TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(G3D%TotalTime))//' ]'
      else  ! Shift the time range to compensate for the shifting of the wind grid
         write (SumFileUnit, '(A)') '     Time range (s):              [ '// &
            TRIM(Num2LStr(-G3D%InitXPosition*G3D%InvMWS))//' : '// &
            TRIM(Num2LStr(G3D%TotalTime - G3D%InitXPosition*G3D%InvMWS))//' ]'
      end if

      write (SumFileUnit, '(A)') '     Y range (m):                 [ '// &
         TRIM(Num2LStr(-G3D%YHWid))//' : '//TRIM(Num2LStr(G3D%YHWid))//' ]'

      if (G3D%NTGrids > 0) then
         write (SumFileUnit, '(A)') '     Z range (m):                 [ '// &
            TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr(G3D%RefHeight + G3D%ZHWid))//' ]'
      else
         write (SumFileUnit, '(A)') '     Z range (m):                 [ '// &
            TRIM(Num2LStr(G3D%RefHeight - G3D%ZHWid))//' : '//TRIM(Num2LStr(G3D%RefHeight + G3D%ZHWid))//' ]'
      end if

      if (G3D%BoxExceedAllowF) then
         write (SumFileUnit, '(A)') '     Wind grid exceedence allowed:  '// &
            'True      -- Only for points requested by OLAF free vortex wake, or LidarSim module'
         write (SumFileUnit, '(A)') '                                    '// &
            '             Out of bounds values are linearly interpolated to mean at Z loction for'
         write (SumFileUnit, '(A)') '                                    '// &
            '             given timestep and X,T value. Values above grid are held to top of wind'
         write (SumFileUnit, '(A)') '                                    '// &
            '             grid value'
      else
         write (SumFileUnit, '(A)') '     Wind grid exceedence allowed:  '// &
            'False'
      end if

      ! We are assuming that if the last line was written ok, then all of them were.
      if (TmpErrStat /= 0_IntKi) then
         call SetErrStat(ErrID_Fatal, 'Error writing to summary file.', ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

end subroutine IfW_Bladed_Init

subroutine Bladed_ReadTurbSimSummary(UnitWind, FileName, CWise, ZCenter, TI, UBar, RefHt, Periodic, LHR, ErrStat, ErrMsg)

   integer(IntKi), intent(in)    :: UnitWind       !< unit number for the file to open
   character(*), intent(in)      :: FileName       !< name of the summary file
   logical, intent(out)          :: CWise          !< rotation (for reading the order of the binary data)
   real(ReKi), intent(out)       :: ZCenter        !< the height at the center of the grid
   real(ReKi), intent(out)       :: TI(3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   real(ReKi), intent(out)       :: UBar           !< mean (advection) wind speed
   real(ReKi), intent(out)       :: RefHt          !< Reference height
   logical, intent(out)          :: Periodic       !< rotation (for reading the order of the binary data)
   logical, intent(out)          :: LHR            !< Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
   integer(IntKi), intent(out)   :: ErrStat        !< returns 0 if no error encountered in the subroutine
   character(*), intent(out)     :: ErrMsg         !< holds the error messages

   character(*), parameter       :: RoutineName = "Bladed_ReadTurbSimSummary"
   integer(IntKi)                :: TmpErrStat     ! temporary error status
   character(ErrMsgLen)          :: TmpErrMsg      ! temporary error message
   real(ReKi)                    :: ZGOffset       ! The vertical offset of the turbine on rectangular grid (allows turbulence not centered on turbine hub)
   integer, parameter            :: NumStrings = 7 ! number of strings to be looking for in the file
   integer(IntKi)                :: FirstIndx      ! The first character of a line where data is located
   integer(IntKi)                :: I              ! A loop counter
   integer(IntKi)                :: LastIndx       ! The last  character of a line where data is located
   integer(IntKi)                :: LineCount      ! Number of lines that have been read in the file
   logical                       :: StrNeeded(NumStrings)   ! if the string has been found
   character(1024)               :: LINE           ! temporary storage for reading a line from the file

   !----------------------------------------------------------------------------------------------
   ! Initialize some variables
   !----------------------------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg = ''

   LineCount = 0
   StrNeeded(:) = .true.
   ZGOffset = 0.0
   RefHt = 0.0
   Periodic = .false.
   LHR = .false.
   CWise = .false.  ! default value, in case it is not in this file

   !----------------------------------------------------------------------------------------------
   ! Open summary file.
   !----------------------------------------------------------------------------------------------

   call OpenFInpFile(UnitWind, TRIM(FileName), TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------------------------
   ! Read the summary file.
   !----------------------------------------------------------------------------------------------

   ! Here are the strings we're looking for, in this order:
   ! 1) 'CLOCKWISE' (optional)
   ! 2) 'HUB HEIGHT'
   ! 3)     (unused; decided we didn't need to read data also stored in the binary file)
   ! 4) 'UBAR'
   ! 5) 'HEIGHT OFFSET' (optional)
   ! 6) 'PERIODIC' (optional)
   ! 7) 'BLADED LEFT-HAND RULE' (optional)

   do while ((ErrStat == ErrID_None) .and. StrNeeded(NumStrings))

      LineCount = LineCount + 1

      read (UnitWind, '(A)', IOSTAT=TmpErrStat) LINE
      if (TmpErrStat /= 0) then

         ! the "HEIGHT OFFSET", "PERIODIC", and "BLADED LEFT-HAND RULE" parameters are not necessary.  We'll assume they are zero/false if we didn't find it.
         ! We will also assume "CLOCKWISE" is false if we didn't find it.
         if (StrNeeded(2) .or. StrNeeded(4)) then
            call SetErrStat(ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'// &
                            TRIM(FileName)//'". Could not find all of the required parameters.', ErrStat, ErrMsg, RoutineName)
            return
         else
            exit
         end if

      end if

      call Conv2UC(LINE)

      if (StrNeeded(2)) then ! if "CLOCKWISE" (StrNeeded(1)) is in the file, we would have already read it. If not, it's not in this file.

         if (StrNeeded(1)) then

            !-------------------------------------------------------------------------------------------
            ! #1: Get the rotation direction, using the string "CLOCKWISE"
            !-------------------------------------------------------------------------------------------

            if (INDEX(LINE, 'CLOCKWISE') > 0) then

               read (LINE, *, IOSTAT=TmpErrStat) CWise          ! Look for True/False values

               if (TmpErrStat /= 0) then                         ! Look for Yes/No values instead

                  LINE = ADJUSTL(LINE)                      ! Remove leading spaces from input line

                  select case (LINE(1:1))
                  case ('Y')
                     CWise = .true.
                  case ('N')
                     CWise = .false.
                  case DEFAULT
                     call SetErrStat(ErrID_Fatal, ' Error reading rotation direction (CLOCKWISE) from FF summary file.', ErrStat, ErrMsg, RoutineName)
                     return
                  end select
                  cycle

               end if ! TmpErrStat /= 0
               StrNeeded(1) = .false.

            end if   ! INDEX for "CLOCKWISE"

         end if

         !-------------------------------------------------------------------------------------------
         ! #2: Get the hub height, using the strings "HUB HEIGHT" or "ZHUB"
         !-------------------------------------------------------------------------------------------

         if (INDEX(LINE, 'HUB HEIGHT') > 0 .or. INDEX(LINE, 'ZHUB') > 0) then

            read (LINE, *, IOSTAT=TmpErrStat) RefHt

            if (TmpErrStat /= 0) then
               call SetErrStat(ErrID_Fatal, ' Error reading hub height from FF summary file.', ErrStat, ErrMsg, RoutineName)
               return
            end if
            StrNeeded(2) = .false.

         end if !INDEX for "HUB HEIGHT" or "ZHUB"

         !      ELSEIF ( StrNeeded(3) ) THEN
         !
         !         !-------------------------------------------------------------------------------------------
         !         ! #3: Get the grid width (& height, if available), using the strings "GRID WIDTH" or "RDIAM"
         !         !    If GRID HEIGHT is specified, use it, too. -- THIS IS UNNECESSARY AS IT'S STORED IN THE BINARY FILE
         !         !-------------------------------------------------------------------------------------------

      elseif (StrNeeded(4)) then

         !-------------------------------------------------------------------------------------------
         ! #4: Get the mean wind speed "UBAR" and turbulence intensities from following lines for
         !     scaling Bladed-style FF binary files
         !-------------------------------------------------------------------------------------------

         if (INDEX(LINE, 'UBAR') > 0) then

            FirstIndx = INDEX(LINE, '=') + 1        ! Look for the equal siqn to find the number we're looking for

            read (LINE(FirstIndx:LEN(LINE)), *, IOSTAT=TmpErrStat) UBar

            if (TmpErrStat /= 0) then
               call SetErrStat(ErrID_Fatal, ' Error reading UBar binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, RoutineName)
               return
            end if

            do I = 1, 3

               LineCount = LineCount + 1

               read (UnitWind, '(A)', IOSTAT=TmpErrStat) LINE
               if (TmpErrStat /= 0) then
                  call SetErrStat(ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'//TRIM(FileName)// &
                                  '". Could not find all of the required parameters.', ErrStat, ErrMsg, RoutineName)
                  return
               end if

               FirstIndx = INDEX(LINE, '=') + 1     ! Read the number between the = and % signs
               LastIndx = INDEX(LINE, '%') - 1

               if (LastIndx <= FirstIndx) LastIndx = LEN(LINE)   ! If there's no % sign, read to the end of the line

               read (LINE(FirstIndx:LastIndx), *, IOSTAT=TmpErrStat) TI(I)
               if (TmpErrStat /= 0) then
                  call SetErrStat(ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(I))// &
                                  ') binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, RoutineName)
                  return
               end if

            end do !I

            StrNeeded(4) = .false.

         end if

      elseif (StrNeeded(5)) then

         !-------------------------------------------------------------------------------------------
         ! #5: Get the grid "HEIGHT OFFSET", if it exists (in TurbSim). Otherwise, assume it's zero
         !           ZGOffset = HH - GridBase - ParamData%FF%ZHWid
         !-------------------------------------------------------------------------------------------
         if (INDEX(LINE, 'HEIGHT OFFSET') > 0) then

            FirstIndx = INDEX(LINE, '=') + 1

            read (LINE(FirstIndx:LEN(LINE)), *, IOSTAT=TmpErrStat) ZGOffset

            if (TmpErrStat /= 0) then
               call SetErrStat(ErrID_Fatal, ' Error reading height offset from FF summary file.', ErrStat, ErrMsg, RoutineName)
               return
            end if

            StrNeeded(5) = .false.

         end if !INDEX for "HEIGHT OFFSET"

      else

         if (StrNeeded(6)) then

            !-------------------------------------------------------------------------------------------
            ! #6: Get the grid "PERIODIC", if it exists (in TurbSim). Otherwise, assume it's
            !        not a periodic file (would only show up if the HEIGHT OFFSET is in the file)
            !-------------------------------------------------------------------------------------------
            if (INDEX(LINE, 'PERIODIC') > 0) then

               Periodic = .true.
               StrNeeded(6) = .false.
               cycle
            end if !INDEX for "PERIODIC"
         end if

         if (StrNeeded(7)) then

            if (INDEX(LINE, 'BLADED LEFT-HAND RULE') > 0) then
               LHR = .true.
               StrNeeded(7) = .false.
            end if ! INDEX for "BLADED LEFT-HAND RULE"

         end if

      end if ! StrNeeded

   end do !WHILE

   !----------------------------------------------------------------------------
   ! Close the summary file
   !----------------------------------------------------------------------------

   close (UnitWind)

   !----------------------------------------------------------------------------
   ! Calculate the height of the grid center
   !----------------------------------------------------------------------------

   ZCenter = RefHt - ZGOffset

end subroutine Bladed_ReadTurbSimSummary

subroutine Bladed_ReadNativeSummary(FileName, PLExp, VLinShr, HLinShr, RefLength, TI, &
                                    UBar, RefHt, PropagationDir, VFlowAngle, BinFileName, &
                                    XOffset, ErrStat, ErrMsg)

   character(*), intent(in)         :: FileName       !< name of the summary file
   real(ReKi), intent(out)          :: PLExp          !< the power-law exponent for vertical wind shear
   real(ReKi), intent(out)          :: VLinShr        !< the linear shape for vertical wind shear
   real(ReKi), intent(out)          :: HLinShr        !< the linear shape for horizontal wind shear
   real(ReKi), intent(out)          :: RefLength      !< Reference (rotor) diameter
   real(ReKi), intent(out)          :: TI(3)          !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   real(ReKi), intent(out)          :: UBar           !< mean (advection) wind speed
   real(ReKi), intent(out)          :: RefHt          !< Reference height
   real(ReKi), intent(out)          :: PropagationDir !< propagation direction
   real(ReKi), intent(out)          :: VFlowAngle     !< vertical flow angle
   character(*), intent(out)        :: BinFileName    !< name of the binary file containing wind data
   real(ReKi), intent(out)          :: XOffset        !< distance offset for start of wind files
   integer(IntKi), intent(out)      :: ErrStat        !< returns 0 if no error encountered in the subroutine
   character(*), intent(out)        :: ErrMsg         !< holds the error messages

   character(*), parameter          :: RoutineName = "Bladed_ReadNativeSummary"
   integer(IntKi)                   :: ErrStat2       ! temporary error status
   character(ErrMsgLen)             :: ErrMsg2        ! temporary error message
   integer(IntKi), parameter        :: UnEc = -1      ! echo file unit number (set to something else > 0 for debugging)
   integer(IntKi)                   :: CurLine        ! Current line to parse in FileInfo data structure

   type(FileInfoType)               :: FileInfo       ! The derived type for holding the file information.

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Open and read the summary file; store data in FileInfo structure.
   !----------------------------------------------------------------------------

   call ProcessComFile(FileName, FileInfo, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   !----------------------------------------------------------------------------
   ! Process the lines stored in FileInfo
   !----------------------------------------------------------------------------

   CurLine = 1

   call ParseVar(FileInfo, CurLine, 'UBAR', UBar, ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call ParseVar(FileInfo, CurLine, 'REFHT', RefHt, ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call ParseVar(FileInfo, CurLine, 'TI', TI(1), ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call ParseVar(FileInfo, CurLine, 'TI_V', TI(2), ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call ParseVar(FileInfo, CurLine, 'TI_W', TI(3), ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call ParseVar(FileInfo, CurLine, 'WDIR', PropagationDir, ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   PropagationDir = R2D*PropagationDir

   call ParseVar(FileInfo, CurLine, 'FLINC', VFlowAngle, ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   VFlowAngle = R2D*VFlowAngle ! convert to degrees

   call ParseVar(FileInfo, CurLine, 'WINDF', BinFileName, ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call ParseVar(FileInfo, CurLine, 'WSHEAR', PLExp, ErrStat2, ErrMsg2, UnEc)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   call ParseVar(FileInfo, CurLine, 'VLINSHEAR', VLinShr, ErrStat2, ErrMsg2, UnEc)
   if (ErrStat2 /= ErrID_None) then
      VLinShr = 0.0_ReKi ! this will be the default if VLINSHEAR is not in the file
   end if

   call ParseVar(FileInfo, CurLine, 'HLINSHEAR', HLinShr, ErrStat2, ErrMsg2, UnEc)
   if (ErrStat2 /= ErrID_None) then
      HLinShr = 0.0_ReKi ! this will be the default if HLINSHEAR is not in the file
   end if

   call ParseVar(FileInfo, CurLine, 'REFLENGTH', RefLength, ErrStat2, ErrMsg2, UnEc)
   if (ErrStat2 /= ErrID_None) then
      RefLength = 0.0_ReKi ! this will be the default if RefLength is not in the file; it will cause an error if either of the linear shears are non-zero
   end if

   call ParseVar(FileInfo, CurLine, 'XOffset', XOffset, ErrStat2, ErrMsg2, UnEc)
   if (ErrStat2 /= ErrID_None) then
      XOffset = 0.0_ReKi ! this will be the default if offset is not in the file
   end if

   !----------------------------------------------------------------------------
   ! Clean FileInfo data structure (including pointers and allocatable arrays)
   !----------------------------------------------------------------------------

   call Cleanup()

contains

   subroutine Cleanup()
      call NWTC_Library_DestroyFileInfoType(FileInfo, ErrStat2, ErrMsg2)
   end subroutine Cleanup

end subroutine Bladed_ReadNativeSummary

!>   Reads the binary headers from the turbulence files of the old Bladed variety.  Note that
!!   because of the normalization, neither ParamData%FF%NZGrids or ParamData%FF%NYGrids are larger than 32 points.
!!   21-Sep-2009 - B. Jonkman, NREL/NWTC.
!!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
subroutine Bladed_ReadHeader0(WindFileUnit, G3D, NativeBladedFmt, ErrStat, ErrMsg)

   integer(IntKi), intent(in)             :: WindFileUnit      !< unit number of already-opened wind file
   type(Grid3DFieldType), intent(inout)   :: G3D               !< Parameters
   logical, intent(in)                    :: NativeBladedFmt   !< Whether this should ignore the advection speed in the binary file
   integer(IntKi), intent(out)            :: ErrStat           !< error status
   character(*), intent(out)              :: ErrMsg            !< error message

   type :: HeaderType
      sequence
      integer(B2Ki)  :: NComp
      integer(B2Ki)  :: DeltaZ, DeltaY, DeltaX
      integer(B2Ki)  :: NStepsHalf, MWS10
      integer(B2Ki)  :: zLu, yLu, xLu, dummy, rnd
      integer(B2Ki)  :: NZ1000, NY1000
   end type

   character(*), parameter    :: RoutineName = "Bladed_ReadHeader0"
   real(ReKi)                 :: FFXDelt
   real(ReKi)                 :: FFYDelt
   real(ReKi)                 :: FFZDelt

   type(HeaderType)           :: header
   integer(B2Ki)              :: dummy(6)
   integer(IntKi)             :: iostat     ! for checking the IOSTAT from a READ or Open statement

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Read the header (file has just been opened)
   !----------------------------------------------------------------------------

   ! Read header
   read (WindFileUnit, IOSTAT=iostat) header
   if (iostat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading header 0 from binary FF file.', &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   G3D%NComp = -1*header%NComp

   FFZDelt = 0.001*header%DeltaZ
   G3D%InvDZ = 1.0/FFZDelt

   FFYDelt = 0.001*header%DeltaY
   G3D%InvDY = 1.0/FFYDelt

   FFXDelt = 0.001*header%DeltaX
   G3D%NSteps = 2*header%NStepsHalf

   if (.not. NativeBladedFmt) G3D%MeanWS = 0.1*header%MWS10
   G3D%InvMWS = 1.0/G3D%MeanWS
   G3D%DTime = FFXDelt/G3D%MeanWS
   G3D%Rate = 1.0/G3D%DTime

   G3D%NZGrids = header%NZ1000/1000
   G3D%ZHWid = 0.5*FFZDelt*(G3D%NZGrids - 1)

   G3D%NYGrids = header%NY1000/1000
   G3D%YHWid = 0.5*FFYDelt*(G3D%NYGrids - 1)

   if (G3D%NComp == 3) then
      read (WindFileUnit, IOSTAT=iostat) dummy
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading header 0 from binary FF file.', &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

end subroutine Bladed_ReadHeader0

subroutine Bladed_ReadHeader1(UnitWind, TI, G3D, NativeBladedFmt, ErrStat, ErrMsg)

   integer(IntKi), intent(in)             :: UnitWind          !< unit number of already-opened wind file
   real(ReKi), intent(out)                :: TI(3)             !< turbulence intensity contained in file header
   type(Grid3DFieldType), intent(inout)   :: G3D               !< Parameters
   logical, intent(in)                    :: NativeBladedFmt   !< Whether this should ignore the advection speed in the binary file
   integer(IntKi), intent(out)            :: ErrStat           !< error status
   character(*), intent(out)              :: ErrMsg            !< error message

   type :: Turb4Type
      sequence
      integer(B4Ki)  :: NComp
      real(SiKi)     :: Latitude, RoughLen, RefHeight, TurbInt(3)
   end type

   type :: Turb78Type
      sequence
      integer(B4Ki)  :: HeaderSize, NComp
   end type

   type :: Turb7Type
      sequence
      real(SiKi)     :: CoherenceDecay, CoherenceScale
   end type

   type :: Turb8Type
      sequence
      real(SiKi)     :: dummy1(6)
      integer(B4Ki)  :: dummy2(3)
      real(SiKi)     :: dummy3(2)
      integer(B4Ki)  :: dummy4(3)
      real(SiKi)     :: dummy5(2)
   end type

   type :: Sub1Type
      sequence
      real(SiKi)     :: DeltaZ, DeltaY, DeltaX
      integer(B4Ki)  :: NStepsHalf
      real(SiKi)     :: MeanWS, zLu, yLu, xLu
      integer(B4Ki)  :: dummy, rnd, NZ, NY
   end type

   type :: Sub2Type
      sequence
      real(SiKi)     :: zLv, yLv, xLv, zLw, yLw, xLw
   end type

   character(*), parameter    :: RoutineName = "Bladed_ReadHeader1"

   type(Turb4Type)            :: Turb4
   type(Turb78Type)           :: Turb78
   type(Sub1Type)             :: Sub1
   type(Sub2Type)             :: Sub2
   type(Turb7Type)            :: Turb7
   type(Turb8Type)            :: Turb8

   real(ReKi)                 :: FFXDelt
   real(ReKi)                 :: FFYDelt
   real(ReKi)                 :: FFZDelt
   integer(B2Ki)              :: Dum_Int2
   integer(B2Ki)              :: TurbType
   integer(IntKi)             :: iostat

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Initialize turbulence intensities
   TI(:) = -1                                                                                !Initialize to -1 (not all models contain TI)

   !----------------------------------------------------------------------------
   ! File reading
   !----------------------------------------------------------------------------

   ! Read 2-byte integer. Can't use library routines for this.
   read (UnitWind, IOSTAT=iostat) Dum_Int2                                                 ! -99 (file ID)
   if (iostat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading integer from binary FF file.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Read 2-byte integer. Can't use library routines for this.
   read (UnitWind, IOSTAT=iostat) TurbType                                                 ! turbulence type
   if (iostat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading turbulence type from binary FF file.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   select case (TurbType)

   case (1, 2)          ! 1-component Von Karman (1) or Kaimal (2)
      G3D%NComp = 1

   case (3, 5)          ! 3-component Von Karman (3) or IEC-2 Kaimal (5)
      G3D%NComp = 3

   case (4)             ! improved Von Karman

      read (UnitWind, IOSTAT=iostat) Turb4
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading number of components from binary FF file.', &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if

      G3D%NComp = Turb4%NComp
      TI = Turb4%TurbInt

   case (7, 8)          ! General Kaimal (7) or  Mann model (8)

      read (UnitWind, IOSTAT=iostat) Turb78
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading number of header records from binary FF file.', &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if

      G3D%NComp = Turb78%NComp

   case DEFAULT

      call SetErrStat(ErrID_Warn, ' InflowWind does not recognize the full-field turbulence file type ='// &
                      TRIM(Num2LStr(int(TurbType)))//'.', ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

   end select

   read (UnitWind, IOSTAT=iostat) Sub1
   if (iostat /= 0) then
      call SetErrStat(ErrID_Fatal, ' Error reading header Sub1 from binary FF file.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   FFZDelt = Sub1%DeltaZ
   G3D%InvDZ = 1.0/FFZDelt

   FFYDelt = Sub1%DeltaY
   G3D%InvDY = 1.0/FFYDelt

   FFXDelt = Sub1%DeltaX

   G3D%NSteps = 2*Sub1%NStepsHalf

   if (.not. NativeBladedFmt) G3D%MeanWS = Sub1%MeanWS
   G3D%InvMWS = 1.0/G3D%MeanWS
   G3D%DTime = FFXDelt/G3D%MeanWS
   G3D%Rate = 1.0/G3D%DTime

   G3D%NZGrids = Sub1%NZ
   G3D%ZHWid = 0.5*FFZDelt*(G3D%NZGrids - 1)    ! half the vertical size of the grid

   G3D%NYGrids = Sub1%NY
   G3D%YHWid = 0.5*FFYDelt*(G3D%NYGrids - 1)

   ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw
   if (G3D%NComp == 3) then
      read (UnitWind, IOSTAT=iostat) Sub2
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
         return
      end if
   end if

   select case (TurbType)
   case (7)       ! General Kaimal model

      ! Unused variables: coherence decay constant and coherence scale parameter in m
      read (UnitWind, IOSTAT=iostat) Turb7
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading header for Turb7 from binary FF file.', &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if

   case (8)       ! Mann model

      ! Unused variables
      read (UnitWind, IOSTAT=iostat) Turb8
      if (iostat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading header for Turb8 from binary FF file.', &
                         ErrStat, ErrMsg, RoutineName)
         return
      end if

   end select !TurbType

end subroutine Bladed_ReadHeader1

subroutine Bladed_ReadGrids(UnitWind, NativeBladedFmt, CWise, LHR, TI, G3D, ErrStat, ErrMsg)

   integer(IntKi), intent(in)             :: UnitWind          !< unit number of already-opened wind file
   logical, intent(in)                    :: NativeBladedFmt   !< whether this data is in native Bladed format (scale to zero mean and unit standard deviation)
   logical, intent(in)                    :: CWise             !< clockwise flag (determines if y is increasing or decreasing in file)
   logical, intent(in)                    :: LHR               !< Left-hand rule for Bladed files (is the v component aligned along *negative* Y?)
   real(ReKi), intent(in)                 :: TI(3)             !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   type(Grid3DFieldType), intent(inout)   :: G3D               !< Parameters
   integer(IntKi), intent(out)            :: ErrStat           !< error status
   character(*), intent(out)              :: ErrMsg            !< error message

   character(*), parameter    :: RoutineName = "Bladed_ReadGrids"
   integer(IntKi)             :: TmpErrStat
   character(ErrMsgLen)       :: TmpErrMsg
   real(ReKi)                 :: FF_Scale(3)    !< used for "un-normalizing" the data
   real(ReKi)                 :: FF_Offset(3)   !< used for "un-normalizing" the data
   integer(B2Ki), allocatable :: raw_ff(:, :, :)
   integer(IntKi)             :: CFirst, CLast, CStep
   integer(IntKi)             :: IC, IR, IT

   ErrMsg = ""
   ErrStat = ErrID_None

   if (NativeBladedFmt) then
      FF_Scale = 0.001_ReKi
      FF_Offset = 0.0_ReKi
   else
      FF_Scale = 0.001_ReKi*G3D%MeanWS*TI/100.0_ReKi
      FF_Offset = (/G3D%MeanWS, 0.0_ReKi, 0.0_ReKi/)  ! used for "un-normalizing" the data
   end if

   ! Bladed convention has positive V pointed along negative Y
   if (LHR) then ! left-hand rule
      FF_Scale(2) = -FF_Scale(2)
   end if

   !----------------------------------------------------------------------------
   ! Generate an informative message
   !----------------------------------------------------------------------------
   ! This could take a while, so we'll write a message to tell users what's going on:

   call WrScr(NewLine//'   Reading a '//TRIM(Num2LStr(G3D%NYGrids))//'x'// &
              TRIM(Num2LStr(G3D%NZGrids))// &
              ' grid ('//TRIM(Num2LStr(G3D%YHWid*2))//' m wide, '// &
              TRIM(Num2LStr(G3D%GridBase))//' m to '// &
              TRIM(Num2LStr(G3D%GridBase + G3D%ZHWid*2))// &
              ' m above ground) with a characteristic wind speed of '// &
              TRIM(Num2LStr(G3D%MeanWS))//' m/s. ')

   !----------------------------------------------------------------------------
   ! Allocate space for the data array
   !----------------------------------------------------------------------------

   ! Add another step, just in case there is an odd number of steps.
   G3D%NSteps = G3D%NSteps + 1

   ! If velocity array is allocated and the size is wrong, deallocate
   if (allocated(G3D%Vel)) then
      if (SIZE(G3D%Vel, 1) /= G3D%NZGrids .or. &
          SIZE(G3D%Vel, 2) /= G3D%NYGrids .or. &
          SIZE(G3D%Vel, 3) /= G3D%NComp .or. &
          SIZE(G3D%Vel, 3) /= G3D%NSteps) then
         deallocate (G3D%Vel)
      end if
   end if

   ! If velocity array isn't allocated, allocate it with proper size
   if (.not. ALLOCATED(G3D%Vel)) then
      call AllocAry(G3D%Vel, G3D%NComp, G3D%NYGrids, G3D%NZGrids, G3D%NSteps, &
                    'Full-field wind data array.', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   end if

   ! Allocate space to read all data for a given time slice
   allocate (raw_ff(G3D%NComp, G3D%NYGrids, G3D%NZGrids), STAT=TmpErrStat)
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, 'Error allocating memory for '// &
                      TRIM(Num2LStr(G3D%NComp*G3D%NYGrids*G3D%NZGrids))// &
                      ' B2Ki in the raw data array.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Initialize the data and set column indexing to account for
   ! direction of turbine rotation (CWise)
   !----------------------------------------------------------------------------

   ! Initialize velocity components not in file to zero
   do IC = G3D%NComp + 1, 3
      G3D%Vel(IC, :, :, :) = 0.0_SiKi
   end do

   if (CWise) then
      CFirst = G3D%NYGrids
      CLast = 1
      CStep = -1
   else
      CFirst = 1
      CLast = G3D%NYGrids
      CStep = 1
   end if

   !----------------------------------------------------------------------------
   ! Loop through all the time steps, reading the data and converting to m/s
   !----------------------------------------------------------------------------

   ! Loop through time steps
   do IT = 1, G3D%NSteps

      ! Read raw data (NComp,NYGrids,NZGrids)
      read (UnitWind, IOStat=TmpErrStat) raw_ff

      ! If data was read successfully, transfer into velocity array, continue
      if (TmpErrStat == 0) then
         do IC = 1, G3D%NComp
            G3D%Vel(IC, :, :, IT) = real(FF_Offset(IC) + FF_Scale(IC)*raw_ff(IC, CFirst:CLast:CStep, :), SiKi)
         end do
         cycle
      end if

      ! If time iterator equals number of steps, then the number of steps were even
      ! so fix the number of steps and exit loop
      if (IT == G3D%NSteps) then
         G3D%NSteps = G3D%NSteps - 1
         ErrStat = ErrID_None
         exit
      end if

      ! Otherwise, an error occurred reading the file, return
      call SetErrStat(ErrID_Fatal, ' Error reading binary data file. '// &
                      'ic = '//TRIM(Num2LStr(ic))// &
                      ', ir = '//TRIM(Num2LStr(ir))// &
                      ', it = '//TRIM(Num2LStr(it))// &
                      ', nsteps = '//TRIM(Num2LStr(G3D%NSteps)), ErrStat, ErrMsg, RoutineName)
      return
   end do

   if (G3D%Periodic) then
      call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(G3D%NSteps))//' time steps of '// &
                 TRIM(Num2LStr(G3D%Rate))//'-Hz full-field data (period of '// &
                 TRIM(Num2LStr(G3D%DTime*G3D%NSteps))//' seconds).')

   else
      call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(G3D%NSteps))//' time steps of '// &
                 TRIM(Num2LStr(G3D%Rate))//'-Hz full-field data ('// &
                 TRIM(Num2LStr(G3D%DTime*(G3D%NSteps - 1)))//' seconds).')
   end if

end subroutine Bladed_ReadGrids

subroutine Bladed_ReadTower(UnitWind, G3D, TwrFileName, ErrStat, ErrMsg)

   integer(IntKi)                         :: UnitWind       !< unit number of wind file to be opened
   type(Grid3DFieldType), intent(inout)   :: G3D            !< Parameters
   character(*), intent(in)               :: TwrFileName
   integer(IntKi), intent(out)            :: ErrStat        !< error status return value (0=no error; non-zero is error)
   character(*), intent(out)              :: ErrMsg         !< a message for errors that occur

   type :: HeaderType
      sequence
      real(SiKi)     :: DZ, DX, Zmax
      integer(B4Ki)  :: NumOutSteps, NumZ
      real(SiKi)     :: UHub, TI(3)
   end type

   character(*), parameter    :: RoutineName = "Bladed_ReadTower"
   integer(IntKi)             :: TmpErrStat     ! IOSTAT value.
   character(ErrMsgLen)       :: TmpErrMsg
   real(ReKi), parameter      :: FF_Offset(3) = (/1.0, 0.0, 0.0/)  ! used for "un-normalizing" the data
   real(ReKi), parameter      :: TOL = 1E-4     ! tolerence for wind file comparisons
   integer(IntKi)             :: IC, IT         ! loop counters
   real(SiKi)                 :: TI(3)          ! scaling values for "un-normalizing the data" [approx. turbulence intensities of the wind components]
   integer(B2Ki), allocatable :: raw_twr(:, :)  ! holds tower velocity for one timestep
   type(HeaderType)           :: header

   ErrMsg = ''
   ErrStat = ErrID_None

   !----------------------------------------------------------------------------
   ! Initialization
   !----------------------------------------------------------------------------

   ! If number of wind components is not three, return with error
   if (G3D%NComp /= 3) then
      call SetErrStat(ErrID_Fatal, ' Error: Tower binary files require 3 wind components.', &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Initialize the number of tower grids to zero
   G3D%NTGrids = 0

   !----------------------------------------------------------------------------
   ! Open the file
   !----------------------------------------------------------------------------

   call OpenBInpFile(UnitWind, TRIM(TwrFileName), TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Read the header information and check that it's compatible with the
   ! FF Bladed-style binary parameters already read.
   !----------------------------------------------------------------------------

   read (UnitWind, IOSTAT=TmpErrStat) header
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, &
                      ' Error reading header of the binary tower file "' &
                      //TRIM(TwrFileName)//'."', ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! If delta Z in file doesn't match tower file
   if (ABS(header%DZ*G3D%InvDZ - 1) > TOL) then
      call SetErrStat(ErrID_Fatal, ' Z resolution in the FF binary file does not match the tower file.', &
                      ErrStat, ErrMsg, RoutineName)
   end if

   ! If X resolution doesn't match the tower file
   if (ABS(header%DX*G3D%InvMWS/G3D%DTime - 1) > TOL) then
      call SetErrStat(ErrID_Fatal, ' Time resolution in the FF binary file does not match the tower file.', &
                      ErrStat, ErrMsg, RoutineName)
   end if

   ! If the height doesn't match the tower file
   if (ABS(header%Zmax/G3D%GridBase - 1) > TOL) then
      call SetErrStat(ErrID_Fatal, ' Height in the FF binary file does not match the tower file "'//TRIM(TwrFileName)//'."', &
                      ErrStat, ErrMsg, RoutineName)
   end if

   ! Number of time steps doesn't match the tower file
   if (header%NumOutSteps /= G3D%NSteps) then
      call SetErrStat(ErrID_Fatal, ' Number of time steps in the FF binary file does not match the tower file.', &
                      ErrStat, ErrMsg, RoutineName)
   end if

   ! If mean wind speed doesn't match tower file
   if (ABS(header%UHub*G3D%InvMWS - 1) > TOL) then
      call SetErrStat(ErrID_Fatal, ' Mean wind speed in the FF binary file does not match the tower file.', &
                      ErrStat, ErrMsg, RoutineName)
   end if

   ! If any of the previous checks failed, or the number of tower grids is zero,
   ! close the wind file and return
   if (ErrStat >= AbortErrLev .or. header%NumZ == 0) then
      close (UnitWind)
      return
   end if

   ! Set number of tower grids from header
   G3D%NTGrids = header%NumZ

   ! Set turbulence intensity from header
   TI = header%TI

   !----------------------------------------------------------------------------
   ! Allocate arrays for the tower points
   !----------------------------------------------------------------------------

   ! If tower array is allocated and the wrong size, deallocate
   if (allocated(G3D%VelTower)) then
      if (size(G3D%VelTower, 1) /= G3D%NComp .or. &
          size(G3D%VelTower, 1) /= G3D%NTGrids .or. &
          size(G3D%VelTower, 1) /= G3D%NSteps) then
         deallocate (G3D%VelTower)
      end if
   end if

   ! If the tower array isn't allocated, allocate it
   if (.not. ALLOCATED(G3D%VelTower)) then
      call AllocAry(G3D%VelTower, G3D%NComp, G3D%NTGrids, G3D%NSteps, &
                    'Tower wind data array.', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   end if

   ! Allocate space to read all data for a given time slice
   allocate (raw_twr(G3D%NComp, G3D%NTGrids), STAT=TmpErrStat)
   if (TmpErrStat /= 0) then
      call SetErrStat(ErrID_Fatal, 'Error allocating memory for '// &
                      TRIM(Num2LStr(G3D%NComp*G3D%NYGrids*G3D%NZGrids))// &
                      ' B2Ki in the raw data array.', ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Read the 16-bit time-series data and scale it to 32-bit reals
   !----------------------------------------------------------------------------

   ! Loop through time.
   do IT = 1, G3D%NSteps

      ! Read wind compnents for this time slice. Can't use library read routines for this.
      read (UnitWind, IOSTAT=TmpErrStat) raw_twr       ! normalized wind-component, INT(2)
      if (TmpErrStat /= 0) then
         call SetErrStat(ErrID_Fatal, ' Error reading binary tower data file. it = '//TRIM(Num2LStr(it))// &
                         ', nsteps = '//TRIM(Num2LStr(G3D%NSteps)), ErrStat, ErrMsg, RoutineName)
         G3D%NTGrids = 0
         return
      end if

      ! Loop through wind components and populate array after scaling to m/s
      do IC = 1, G3D%NComp
         G3D%VelTower(IC, :, IT) = real(G3D%MeanWS*(FF_Offset(IC) + 0.00001*TI(IC)*raw_twr(IC, :)), SiKi)
      end do
   end do

   !----------------------------------------------------------------------------
   ! Close the file
   !----------------------------------------------------------------------------

   close (UnitWind)

   call WrScr(NewLine//'   Processed '//TRIM(Num2LStr(G3D%NSteps))//' time steps of '// &
              TRIM(Num2LStr(G3D%NTGrids))//'x1 tower data grids.')

end subroutine Bladed_ReadTower

subroutine Grid3D_PopulateWindFileDat(Grid3DField, FileName, WindType, HasTower, FileDat)

   type(Grid3DFieldType), intent(in)   :: Grid3DField
   character(*), intent(in)            :: FileName
   integer(IntKi), intent(in)          :: WindType
   logical, intent(in)                 :: HasTower
   type(WindFileDat), intent(out)      :: FileDat

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
   else
      ! Shift the time range to compensate for the shifting of the wind grid
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

subroutine Grid3D_AddMeanVelocity(InitInp, G3D)

   type(Grid3D_InitInputType), intent(in)    :: InitInp  !< Initialization input data passed to the module
   type(Grid3DFieldType), intent(inout)      :: G3D      !< Initialization input data passed to the module

   real(ReKi)        :: Z           ! height
   real(ReKi)        :: Y           ! distance from centre in horizontal direction
   real(ReKi)        :: U           ! mean wind speed
   integer(IntKi)    :: iz, iy      ! loop counter
   integer(IntKi)    :: centre_y    ! index of centre in y direction
   integer(IntKi)    :: first_positive_iz ! index of first height that is above ground

   ! Loop through grid elevations
   first_positive_iz = 1
   do iz = 1, G3D%NZGrids

      ! calculate height
      Z = G3D%GridBase + (iz - 1)/G3D%InvDZ
      if (Z <= 0.0_ReKi) then
         first_positive_iz = iz + 1
         cycle
      end if

      ! Calculate wind speed to add based on profile
      select case (G3D%WindProfileType)

      case (WindProfileType_PL)
         U = G3D%MeanWS*(Z/G3D%RefHeight)**G3D%PLExp ! [IEC 61400-1 6.3.1.2 (10)]

      case (WindProfileType_Log)
         if (.not. EqualRealNos(G3D%RefHeight, G3D%Z0) .and. Z > 0.0_ReKi) then
            U = G3D%MeanWS*(LOG(Z/G3D%Z0))/(LOG(G3D%RefHeight/G3D%Z0))
         else
            U = 0.0_ReKi
         end if

      case (WindProfileType_Constant)
         U = G3D%MeanWS

      case DEFAULT
         U = 0.0_ReKi

      end select

      ! Add vertical linear shear, if nonzero
      if (InitInp%VLinShr /= 0.0_ReKi) then
         U = U + G3D%MeanWS*InitInp%VLinShr*(Z - G3D%RefHeight)/G3D%RefLength
      end if

      ! Add velocity
      G3D%Vel(1, :, iz, :) = real(G3D%Vel(1, :, iz, :) + U, SiKi)

   end do

   ! Add horizontal linear shear, if nonzero (only to the points above the ground)
   if (InitInp%HLinShr /= 0.0_ReKi .and. first_positive_iz <= G3D%NZGrids) then

      ! find the center point of the grid (if we don't have an odd number of grid points, we'll pick the point closest to the center)
      centre_y = (G3D%NYGrids + 1)/2 ! integer division

      ! Loop through grid Y coordinates
      do iy = 1, G3D%NYGrids
         Y = (iy - centre_y)/G3D%InvDY
         U = G3D%MeanWS*InitInp%HLinShr*Y/G3D%RefLength
         G3D%Vel(1, iy, first_positive_iz:, :) = real(G3D%Vel(1, iy, first_positive_iz:, :) + U, SiKi)
      end do
   end if

end subroutine Grid3D_AddMeanVelocity

subroutine Grid3D_ScaleTurbulence(InitInp, Vel, ScaleFactors, ErrStat, ErrMsg)

   type(Grid3D_InitInputType), intent(in)   :: InitInp           !< Initialization input data passed to the module
   real(SiKi), intent(INOUT)               :: Vel(:, :, :, :)   !< full-field wind inflow data
   real(ReKi), intent(out)                 :: ScaleFactors(3)   !< scaling factors that were used
   integer(IntKi), intent(out)             :: ErrStat           !< determines if an error has been encountered
   character(*), intent(out)               :: ErrMsg            !< Message about errors

   character(*), parameter    :: RoutineName = 'Grid3D_ScaleTurbulence'
   real(DbKi)                 :: vMean(3)          ! average wind speeds over time at target position
   real(DbKi)                 :: vSum(3)           ! sum over time of wind speeds at target position
   real(DbKi)                 :: vSum2(3)          ! sum of wind speeds squared
   real(ReKi)                 :: ActualSigma(3)    ! computed standard deviation

   integer                    :: ic                ! Loop counter for wind component
   integer                    :: iy                ! Loop counter for y
   integer                    :: iz                ! Loop counter for z

   integer                    :: nc                ! number of FF wind components
   integer                    :: nt                ! size of x (or t) dimension of turbulence box
   integer                    :: ny                ! size of y dimension of turbulence box
   integer                    :: nz                ! size of z dimension of turbulence box

   ErrStat = ErrID_None
   ErrMsg = ""

   nc = size(Vel, 1)
   ny = size(Vel, 2)
   nz = size(Vel, 3)
   nt = size(Vel, 4)

   ! If scaling method is none, set factors to 1 and return (no scaling)
   if (InitInp%ScaleMethod == ScaleMethod_None) then
      ScaleFactors = 1.0_ReKi
      return
   end if

   !----------------------------------------------------------------------------
   ! Determine the scaling factors:
   !----------------------------------------------------------------------------

   ! Use the scaling factors specified in the input file
   if (InitInp%ScaleMethod == ScaleMethod_Direct) then
      ScaleFactors = InitInp%sf

   else ! compute the scaling factors to get requested sigma:

      ! find the center point of the grid (if we don't have an odd number of grid points, we'll pick the point closest to the center)
      iz = (nz + 1)/2 ! integer division
      iy = (ny + 1)/2 ! integer division

      ! compute the actual sigma at the point specified by (iy,iz). (This sigma should be close to 1.)
      vSum = sum(Vel(:, iy, iz, :), dim=2)
      vSum2 = sum(Vel(:, iy, iz, :)**2, dim=2)
      vMean = vSum/nt
      ActualSigma = real(SQRT(ABS((vSum2/nt) - vMean**2)), ReKi)

      ! check that the ActualSigma isn't 0
      ! InitOut%sf = InitInp%SigmaF / ActualSigma  ! factor = Target / actual
      do ic = 1, nc
         if (EqualRealNos(ActualSigma(ic), 0.0_ReKi)) then
            ScaleFactors(ic) = 0.0_ReKi
            if (.not. EqualRealNos(InitInp%SigmaF(ic), 0.0_ReKi)) then
               call SetErrStat(ErrID_Fatal, "Computed standard deviation is zero; cannot scale to achieve target non-zero standard deviation.", &
                               ErrStat, ErrMsg, RoutineName)
            end if
         else
            ScaleFactors(ic) = InitInp%SigmaF(ic)/ActualSigma(ic)
         end if
      end do

   end if

   !----------------------------------------------------------------------------
   ! scale the data using our scaling factors:
   !----------------------------------------------------------------------------

   do ic = 1, nc
      Vel(ic, :, :, :) = real(ScaleFactors(ic)*Vel(ic, :, :, :), SiKi)
   end do

end subroutine Grid3D_ScaleTurbulence

subroutine Grid3D_ValidateInput(InitInp, NComp, ErrStat, ErrMsg)

   type(Grid3D_InitInputType), intent(in)    :: InitInp     !< Initialization input data passed to the module
   integer(IntKi), intent(in)                :: NComp       !< number of full-field wind components (normally 3)

   character(*), parameter                   :: RoutineName = 'Grid3D_ValidateInput'
   integer(IntKi), intent(out)               :: ErrStat     !< determines if an error has been encountered
   character(*), intent(out)                 :: ErrMsg      !< Message about errors

   ErrStat = ErrID_None
   ErrMsg = ""

   if (InitInp%RefHt < 0.0_ReKi .or. EqualRealNos(InitInp%RefHt, 0.0_ReKi)) call SetErrStat(ErrID_Fatal, 'The grid reference height must be larger than 0.', ErrStat, ErrMsg, RoutineName)

   if (InitInp%ScaleMethod == ScaleMethod_Direct) then
      if (any(InitInp%sf < 0.0_ReKi)) call SetErrStat(ErrID_Fatal, 'Turbulence scaling factors must not be negative.', ErrStat, ErrMsg, RoutineName)
   elseif (InitInp%ScaleMethod == ScaleMethod_StdDev) then
      if (any(InitInp%sigmaf < 0.0_ReKi)) call SetErrStat(ErrID_Fatal, 'Turbulence standard deviations must not be negative.', ErrStat, ErrMsg, RoutineName)
   elseif (InitInp%ScaleMethod /= ScaleMethod_None) then
      call SetErrStat(ErrID_Fatal, 'Turbulence scaling method must be 0 (none), 1 (direct scaling factors), or 2 (target standard deviation).', ErrStat, ErrMsg, RoutineName)
   end if

   if (InitInp%WindProfileType == WindProfileType_Log) then
      if (InitInp%z0 < 0.0_ReKi .or. EqualRealNos(InitInp%z0, 0.0_ReKi)) &
         call SetErrStat(ErrID_Fatal, 'The surface roughness length, Z0, must be greater than zero', ErrStat, ErrMsg, RoutineName)
   elseif (InitInp%WindProfileType < WindProfileType_Constant .or. InitInp%WindProfileType > WindProfileType_PL) then
      call SetErrStat(ErrID_Fatal, 'The WindProfile type must be 0 (constant), 1 (logarithmic) or 2 (power law).', ErrStat, ErrMsg, RoutineName)
   end if

   if (InitInp%URef < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'The reference wind speed must not be negative.', ErrStat, ErrMsg, RoutineName)

   if (EqualRealNos(InitInp%RefLength, 0.0_ReKi) .or. InitInp%RefLength < 0.0_ReKi) then
      if (InitInp%VLinShr /= 0.0_ReKi .or. InitInp%HLinShr /= 0.0_ReKi) then
         call SetErrStat(ErrID_Fatal, 'The reference length must be a positive number when vertical or horizontal linear shear is used.', ErrStat, ErrMsg, RoutineName)
      end if
   end if

end subroutine

subroutine Grid3D_WriteBladed(G3D, FileRootName, unit, ErrStat, ErrMsg)

   type(Grid3DFieldType), intent(in)  :: G3D             !< Parameters
   character(*), intent(in)           :: FileRootName    !< Name of the file to write the output in
   integer(IntKi), intent(in)         :: Unit            !< Indicates whether an error occurred (see NWTC_Library)
   integer(IntKi), intent(out)        :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   character(*), intent(out)          :: ErrMsg          !< Error message associated with the ErrStat

   character(*), parameter       :: RoutineName = 'Grid3D_WriteBladed'
   real(SiKi), parameter         :: Tolerance = 0.0001   ! The largest difference between two numbers that are assumed to be equal
   integer(IntKi)                :: ic, it, iy, iz
   real(SiKi), allocatable       :: MeanVal(:, :)
   real(SiKi), allocatable       :: SigmaGrid(:, :)
   real(SiKi)                    :: TI(3)                !< array containing turbulence intensity (for scaling factors)
   real(SiKi)                    :: Sigma(3)             !< array containing standard deviations (for scaling factors)
   real(SiKi)                    :: Scl(3)               !< array containing scaling factors
   real(SiKi)                    :: Off(3)               !< array containing offsets
   real(SiKi)                    :: Tmp
   real(ReKi)                    :: MeanWS_nonZero       !< advection speed (mean wind speed at hub)
   real(ReKi)                    :: delta(3)
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   delta = [G3D%MeanWS*G3D%DTime, 1.0_ReKi/G3D%InvDY, 1.0_ReKi/G3D%InvDZ]

   call AllocAry(MeanVal, G3D%NYGrids, G3D%NZGrids, "MeanVal", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   call AllocAry(SigmaGrid, G3D%NYGrids, G3D%NZGrids, "SigmaGrid", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Loop through components
   do ic = 3, 1, -1

      ! mean values:
      MeanVal = sum(G3D%Vel(ic, :, :, :), dim=3)/G3D%NSteps

      ! standard deviations (with 1/N scaling factor):
      SigmaGrid = sum(G3D%Vel(ic, :, :, :)**2, dim=3)/G3D%NSteps
      SigmaGrid = sqrt(max(SigmaGrid - MeanVal**2, 0.0_SiKi))

      ! now get the average standard deviation for each component:
      Sigma(ic) = sum(SigmaGrid)/size(SigmaGrid) ! get the average sigma over the grid
      Sigma(ic) = max(100.0_SiKi*Tolerance, Sigma(ic)) ! make sure this scaling isn't too small

   end do

   ! We need to take into account the shear across the grid in the sigma calculations for scaling the data,
   ! and ensure that 32.767*sigma_u >= |V-UHub| so that we don't get values out of the range of our scaling values
   ! in this BLADED-style binary output.  Tmp is |V-UHub|
   ! Get the range of wind speed values for scaling in BLADED-format .wnd files
   Tmp = real(max(abs(maxval(G3D%Vel(:, :, 1, :)) - G3D%MeanWS), abs(minval(G3D%Vel(1, :, :, :)) - G3D%MeanWS)), SiKi)
   Sigma(1) = max(Sigma(1), 0.05_SiKi*Tmp)
   do ic = 2, 3
      ! put the abs() after the maxval() and minval() to avoid stack-overflow issues with large wind files
      Sigma(ic) = max(Sigma(ic), 0.05_SiKi*abs(maxVAL(G3D%Vel(ic, :, :, :))), 0.05_SiKi*abs(minval(G3D%Vel(ic, :, :, :))))
   end do

   ! Put normalizing factors into the summary file.  The user can use them to
   ! tell a simulation program how to rescale the data.
   if (abs(G3D%MeanWS) < 0.1_ReKi) then
      MeanWS_nonZero = sign(0.1_ReKi, G3D%MeanWS)
   else
      MeanWS_nonZero = G3D%MeanWS
   end if

   TI = real(Sigma/MeanWS_nonZero, SiKi)

   !----------------------------------------------------------------------------
   ! The summary file
   !----------------------------------------------------------------------------

   call OpenFOutFile(unit, trim(FileRootName)//'-Bladed.sum', ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! The string "TurbSim" needs to be in the 2nd line of the summary file if AeroDyn will read this.
   write (unit, "( / 'TurbSim - This summary file was generated by ', A, ' on ' , A , ' at ' , A , '.' / )") "NWTC_Library", CurDate(), CurTime()
   write (unit, '(/)')
   write (unit, '(/)')
   write (unit, '( L10,   2X, "Clockwise rotation when looking downwind?")') .false.
   write (unit, '( F10.3, 2X, "Hub height [m]")') G3D%RefHeight
   write (unit, '( F10.3, 2X, "Grid height [m]")') delta(3)*(G3D%NZGrids - 1)
   write (unit, '( F10.3, 2X, "Grid width [m]")') delta(2)*(G3D%NYGrids - 1)
   write (unit, '(/"BLADED-style binary scaling parameters:"/)')
   write (unit, '( 2X, "UBar  = ", F9.4, " m/s")') MeanWS_nonZero
   write (unit, '( 2X, "TI(u) = ", F9.4, " %")') 100.0*TI(1)
   write (unit, '( 2X, "TI(v) = ", F9.4, " %")') 100.0*TI(2)
   write (unit, '( 2X, "TI(w) = ", F9.4, " %")') 100.0*TI(3)
   write (unit, '(/)')
   write (unit, '( 2X, "Height offset = ", F9.4, " m" )') G3D%RefHeight - 0.5*delta(3)*(G3D%NZGrids - 1) - G3D%GridBase
   write (unit, '( 2X, "Grid Base     = ", F9.4, " m" )') G3D%GridBase
   if (G3D%Periodic) then
      write (unit, '()')
      write (unit, '( A)') 'Creating a PERIODIC output file.'
   end if
   write (unit, '( A)') 'Creating a BLADED LEFT-HAND RULE output file.'

   close (unit)

   !----------------------------------------------------------------------------
   ! The BINARY file
   !----------------------------------------------------------------------------

   call OpenBOutFile(unit, TRIM(FileRootName)//'-Bladed.wnd', ErrStat, ErrMsg)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   write (unit) int(-99, B2Ki)               ! -99 = New Bladed format
   write (unit) int(4, B2Ki)               ! 4 = improved von karman (not used, but needed for next 7 inputs)
   write (unit) int(size(G3D%Vel, 3), B4Ki)               ! size(FFWind,3) = 3 = number of wind components
   write (unit) real(45.0_SiKi, SiKi)               ! Latitude (degrees)   (informational, not used in FAST)
   write (unit) real(0.03_SiKi, SiKi)               ! Roughness length (m) (informational, not used in FAST)
   write (unit) real(G3D%RefHeight, SiKi)               ! Reference Height (m) (informational, not used in FAST)
   write (unit) real(100.0*TI(1), SiKi)               ! Longitudinal turbulence intensity (%)
   write (unit) real(100.0*TI(2), SiKi)               ! Lateral turbulence intensity (%)
   write (unit) real(100.0*TI(3), SiKi)               ! Vertical turbulence intensity (%)

   write (unit) real(delta(3), SiKi)               ! grid spacing in vertical direction, in m
   write (unit) real(delta(2), SiKi)               ! grid spacing in lateral direction, in m
   write (unit) real(delta(1), SiKi)               ! grid spacing in longitudinal direciton, in m
   write (unit) int(G3D%NSteps/2, B4Ki)               ! half the number of points in alongwind direction
   write (unit) real(MeanWS_nonZero, SiKi)               ! the mean wind speed in m/s
   write (unit) real(0, SiKi)               ! the vertical length scale of the longitudinal component in m
   write (unit) real(0, SiKi)               ! the lateral length scale of the longitudinal component in m
   write (unit) real(0, SiKi)               ! the longitudinal length scale of the longitudinal component in m
   write (unit) int(0, B4Ki)               ! an unused integer
   write (unit) int(0, B4Ki)               ! the random number seed
   write (unit) int(G3D%NZGrids, B4Ki)               ! the number of grid points vertically
   write (unit) int(G3D%NYGrids, B4Ki)               ! the number of grid points laterally
   write (unit) int(0, B4Ki)               ! the vertical length scale of the lateral component, not used
   write (unit) int(0, B4Ki)               ! the lateral length scale of the lateral component, not used
   write (unit) int(0, B4Ki)               ! the longitudinal length scale of the lateral component, not used
   write (unit) int(0, B4Ki)               ! the vertical length scale of the vertical component, not used
   write (unit) int(0, B4Ki)               ! the lateral length scale of the vertical component, not used
   write (unit) int(0, B4Ki)               ! the longitudinal length scale of the vertical component, not used

   ! Scaling value to convert wind speeds to 16-bit integers
   do ic = 1, 3
      if (.not. EqualRealNos(Sigma(ic), 0.0_SiKi)) then
         Scl(ic) = 1000.0/(Sigma(ic))
      else
         Scl(ic) = 1.0_SiKi
      end if
   end do

   ! Bladed convention is positive V is pointed along negative Y (IEC turbine coordinate)
   Scl(2) = -Scl(2)

   ! Offset value to convert wind speeds to 16-bit integers
   if (G3D%AddMeanAfterInterp) then ! Note that this will not take into account any shear!!!
      Off(1) = 0.0
   else
      Off(1) = real(G3D%MeanWS*Scl(1), SiKi)
   end if
   Off(2) = 0.0
   Off(3) = 0.0

   ! Scale velocity for 16-bit integers and write to file
   do it = 1, G3D%NSteps
      do iz = 1, G3D%NZGrids
         do iy = 1, G3D%NYGrids
            write (unit) NINT(G3D%Vel(:, iy, iz, it)*Scl - Off, B2Ki) ! scale to int16
         end do !IY
      end do !IZ
   end do !IT

   close (unit)

end subroutine Grid3D_WriteBladed

subroutine Grid3D_WriteVTK(G3D, FileRootName, unit, ErrStat, ErrMsg)

   type(Grid3DFieldType), intent(in)   :: G3D            !< Parameters
   character(*), intent(in)            :: FileRootName   !< RootName for output files
   integer(IntKi), intent(in)           :: unit           !< Error status of the operation
   integer(IntKi), intent(out)          :: ErrStat        !< Error status of the operation
   character(*), intent(out)            :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter                :: RoutineName = 'ConvertField_toVTK'
   character(1024)                        :: RootPathName
   character(1024)                        :: FileName
   integer                                :: i
   integer                                :: iy
   integer                                :: iz
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2

   call GetPath(FileRootName, RootPathName)
   RootPathName = trim(RootPathName)//PathSep//"vtk"
   call MkDir(trim(RootPathName))  ! make this directory if it doesn't already exist

   ! Loop through time steps
   do i = 1, G3D%NSteps

      ! Create the output vtk file with naming <WindFilePath>/vtk/DisYZ.t<i>.vtk
      FileName = trim(RootPathName)//PathSep//"DisYZ.t"//trim(num2lstr(i))//".vtp"

      ! see WrVTK_SP_header
      call OpenFOutFile(unit, TRIM(FileName), ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      write (unit, '(A)') '# vtk DataFile Version 3.0'
      write (unit, '(A)') "InflowWind YZ Slice at T= "//trim(num2lstr((i - 1)*G3D%DTime))//" s"
      write (unit, '(A)') 'ASCII'
      write (unit, '(A)') 'DATASET STRUCTURED_POINTS'

      ! Note: gridVals must be stored such that the left-most dimension is X
      ! and the right-most dimension is Z (see WrVTK_SP_vectors3D)
      write (unit, '(A,3(i5,1X))') 'DIMENSIONS ', 1, G3D%NYGrids, G3D%NZGrids
      write (unit, '(A,3(f10.2,1X))') 'ORIGIN ', G3D%InitXPosition, -G3D%YHWid, G3D%GridBase
      write (unit, '(A,3(f10.2,1X))') 'SPACING ', 0.0_ReKi, 1.0_ReKi/G3D%InvDY, 1.0_ReKi/G3D%InvDZ
      write (unit, '(A,i5)') 'POINT_DATA ', G3D%NYGrids*G3D%NZGrids
      write (unit, '(A)') 'VECTORS DisYZ float'

      do iz = 1, G3D%NZGrids
         do iy = 1, G3D%NYGrids
            write (unit, '(3(f10.2,1X))') G3D%Vel(:, iy, iz, i)
         end do
      end do

      close (unit)

   end do

end subroutine Grid3D_WriteVTK

subroutine Grid3D_WriteHAWC(G3D, FileRootName, unit, ErrStat, ErrMsg)

   character(*), intent(in)            :: FileRootName   !< Name of the file to write the output in
   type(Grid3DFieldType), intent(in)   :: G3D            !< Parameters
   integer(IntKi), intent(in)          :: unit           !< Error status of the operation
   integer(IntKi), intent(out)         :: ErrStat        !< Indicates whether an error occurred (see NWTC_Library)
   character(*), intent(out)           :: ErrMsg         !< Error message associated with the ErrStat

   character(*), parameter       :: RoutineName = 'Grid3D_WriteHAWC'
   character(*), parameter       :: Comp(3) = (/'u', 'v', 'w'/)
   real(ReKi)                    :: delta(3)
   integer(IntKi)                :: IC, IX, IY, IZ
   real(SiKi), allocatable       :: MeanVal(:)
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(1024)               :: RootWithoutPathName

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Array containing dx, dy, dz in meters
   delta = [G3D%MeanWS*G3D%DTime, 1.0_ReKi/G3D%InvDY, 1.0_ReKi/G3D%InvDZ]

   ! Allocate array to hold mean value by Z location
   call AllocAry(MeanVal, G3D%NZGrids, "MeanVal", ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Calculate mean value at each Z location
   MeanVal = sum(G3D%Vel(1, 1, :, :), dim=2)/G3D%NSteps

   !----------------------------------------------------------------------------
   ! Write summary file
   !----------------------------------------------------------------------------

   call OpenFOutFile(unit, trim(FileRootName)//'-HAWC.sum', ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   write (unit, '(A)') '; Wind file converted to HAWC format on '//CurDate()//' at '//CurTime()

   write (unit, '()')
   do IZ = G3D%NZGrids, 1, -1
      write (unit, '(A,I3,A,F15.5)') '; mean removed at z(', iz, ') = ', MeanVal(iz)
   end do

   write (unit, '(A)') 'turb_format 1 ;'

   write (unit, '()')
   write (unit, '(A)') 'begin mann;'

   ic = INDEX(FileRootName, '\', BACK=.true.)
   ic = MAX(ic, INDEX(FileRootName, '/', BACK=.true.))
   RootWithoutPathName = FileRootName((ic + 1):)

   write (unit, '(2x,A, T30, A, " ;")') 'filename_u', trim(RootWithoutPathName)//'-HAWC-u.bin'
   write (unit, '(2x,A, T30, A, " ;")') 'filename_v', trim(RootWithoutPathName)//'-HAWC-v.bin'
   write (unit, '(2x,A, T30, A, " ;")') 'filename_w', trim(RootWithoutPathName)//'-HAWC-w.bin'

   write (unit, '(2x,A, T30, I8, 1x, F15.5, " ;")') 'box_dim_u', G3D%NSteps, delta(1)
   write (unit, '(2x,A, T30, I8, 1x, F15.5, " ;")') 'box_dim_v', G3D%NYGrids, delta(2)
   write (unit, '(2x,A, T30, I8, 1x, F15.5, " ;")') 'box_dim_w', G3D%NZGrids, delta(3)

   write (unit, '(2x,A)') 'dont_scale 1;  converter did not rescale turbulence to unit standard deviation'
   write (unit, '(A)') 'end mann;'
   close (unit)

   !----------------------------------------------------------------------------
   ! Write the binary files for each component
   !----------------------------------------------------------------------------

   do IC = 1, G3D%NComp

      call OpenBOutFile(unit, trim(FileRootName)//'-HAWC-'//Comp(ic)//'.bin', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      do IX = 1, G3D%NSteps
         do IY = G3D%NYGrids, 1, -1
            write (unit, IOSTAT=ErrStat2) real(G3D%Vel(ic, iy, :, ix) - MeanVal, SiKi)
         end do
      end do

      close (unit)

      MeanVal = 0.0_SiKi

   end do

end subroutine Grid3D_WriteHAWC

end module InflowWind_IO
