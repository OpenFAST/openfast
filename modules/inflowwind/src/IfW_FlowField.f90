!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2022  National Renewable Energy Laboratory
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
!**********************************************************************************************************************************

module IfW_FlowField

use NWTC_Library
use IfW_FlowField_Types

implicit none

public IfW_FlowField_GetVelAcc
public IfW_UniformField_CalcAccel, IfW_Grid3DField_CalcAccel
public IfW_UniformWind_GetOP
public Grid3D_to_Uniform, Uniform_to_Grid3D

integer(IntKi), parameter  :: WindProfileType_None = -1     !< don't add wind profile; already included in input
integer(IntKi), parameter  :: WindProfileType_Constant = 0  !< constant wind
integer(IntKi), parameter  :: WindProfileType_Log = 1       !< logarithmic
integer(IntKi), parameter  :: WindProfileType_PL = 2        !< power law

real(ReKi), parameter      :: GridTol = 1.0E-3              ! Tolerance for determining if position is within grid

contains

!> IfW_FlowField_GetVelAcc gets the velocities (and accelerations) at the given point positions.
!! Accelerations are only calculated if the AccelUVW array is allocated.
subroutine IfW_FlowField_GetVelAcc(FF, IStart, Time, PositionXYZ, VelocityUVW, AccelUVW, ErrStat, ErrMsg)

   type(FlowFieldType), intent(in)           :: FF                !< FlowField data structure
   integer(IntKi), intent(in)                :: IStart            !< Start index for returning velocities for external field
   real(DbKi), intent(in)                    :: Time              !< Time to evaluate velocities/accelerations
   real(ReKi), intent(in)                    :: PositionXYZ(:, :) !< Array of positions to evaluate velocites/accelerations
   real(ReKi), intent(inout)                 :: VelocityUVW(:, :) !< Array of velocity outputs
   real(ReKi), allocatable, intent(inout)    :: AccelUVW(:, :)    !< Array of acceleration outputs
   integer(IntKi), intent(out)               :: ErrStat           !< Error status
   character(*), intent(out)                 :: ErrMsg            !< Error message

   character(*), parameter                   :: RoutineName = "IfW_FlowField_GetVelAcc"
   integer(IntKi)                            :: i
   integer(IntKi)                            :: NumPoints
   logical                                   :: OutputAccel, AddMeanAfterInterp
   real(ReKi), allocatable                   :: Position(:, :)
   integer(IntKi)                            :: Grid3D_AccelInterp
   integer(IntKi)                            :: TmpErrStat
   character(ErrMsgLen)                      :: TmpErrMsg

   ! Uniform Field
   type(UniformField_Interp)                 :: UFopVel, UFopAcc

   ! Grid3D Field
   real(ReKi)                                :: Xi(3)
   real(ReKi)                                :: VelCell(8, 3), AccCell(8, 3)
   logical                                   :: Is3D
   logical                                   :: GridExceedAllow   ! is this point allowed to exceed bounds of wind grid

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get number of points to evaluate
   NumPoints = size(PositionXYZ, dim=2)

   ! Determine if acceleration should be calculated and returned
   OutputAccel = allocated(AccelUVW)
   if (OutputAccel .and. .not. FF%AccFieldValid) then
      call SetErrStat(ErrID_Fatal, "Accel output requested, but accel field is not valid", &
                      ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Allocate position array
   call AllocAry(Position, 3, NumPoints, "Rotated position data", TmpErrStat, TmpErrMsg)
   if (TmpErrStat >= AbortErrLev) then
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   !  Wind speed coordinate transforms from Wind file coordinates to global
   !----------------------------------------------------------------------------

   ! The VelocityUVW array data that has been returned from the sub-modules is
   ! in the wind file (X'Y'Z') coordinates. These must be rotated to the global
   ! XYZ coordinates if PropagationDir or UpflowAngle is nonzero. Apply the
   ! coordinate transformation to the VelocityUVW(prime) coordinates
   ! (in wind X'Y'Z' coordinate frame) returned from the submodules to the XYZ
   ! coordinate frame, if PropagationDir is not zero. This is only a rotation
   ! of the returned wind field, so UVW contains the direction components of
   ! the wind at XYZ after translation from the U'V'W' wind velocity components
   ! in the X'Y'Z' (wind file) coordinate frame.
   ! NOTE: rotations are about the hub at [ 0 0 H ].

   ! Copy positions or transform based on wind box rotation
   if (FF%RotateWindBox) then
      do i = 1, NumPoints
         Position(:, i) = GetPrimePosition(PositionXYZ(:, i))
      end do
   else
      Position = PositionXYZ
   end if

   !----------------------------------------------------------------------------
   ! Get velocities/accelerations based on flow field type
   !----------------------------------------------------------------------------

   ! Switch based on flow type
   select case (FF%FieldType)
   case (Uniform_FieldType)

      !-------------------------------------------------------------------------
      ! Uniform Flow Field
      !-------------------------------------------------------------------------

      ! If cubic velocity interp requested, calculate operating point using
      ! cubic interpolation. Otherwise, use linear interpolation
      if (FF%VelInterpCubic) then
         UFopVel = UniformField_InterpCubic(FF%Uniform, Time)
      else
         UFopVel = UniformField_InterpLinear(FF%Uniform, Time)
      end if

      ! If velocity and acceleration output is requested
      if (OutputAccel) then

         ! If cubic interpolation was used for the velocity, use same for accel
         ! otherwise, calculate operating point via cubic interpolation
         if (FF%VelInterpCubic) then
            UFopAcc = UFopVel
         else
            UFopAcc = UniformField_InterpCubic(FF%Uniform, Time)
         end if

         ! Loop throuh points and calcualate velocity and acceleration
         do i = 1, NumPoints
            if (Position(3, i) > 0.0_ReKi) then
               VelocityUVW(:, i) = UniformField_GetVel(FF%Uniform, UFopVel, Position(:, i))
               AccelUVW(:, i) = UniformField_GetAcc(FF%Uniform, UFopAcc, Position(:, i))
            else
               VelocityUVW(:, i) = 0.0_ReKi
               AccelUVW(:, i) = 0.0_ReKi
            end if
         end do

      else  ! Otherwise, only velocity requested

         ! Loop throuh points and calcualate velocity
         do i = 1, NumPoints
            if (Position(3, i) > 0.0_ReKi) then
               VelocityUVW(:, i) = UniformField_GetVel(FF%Uniform, UFopVel, Position(:, i))
            else
               VelocityUVW(:, i) = 0.0_ReKi
            end if
         end do
      end if

   case (Grid3D_FieldType)

      !-------------------------------------------------------------------------
      ! Grid3D Flow Field
      !-------------------------------------------------------------------------

      ! Determine select case value for calculating acceleration and 
      ! interpolation so it doesn't have to be done in the loop (optimization)
      if (OutputAccel) then
         if (FF%VelInterpCubic) then
            Grid3D_AccelInterp = 1     ! Output accel, cubic interp
         else
            Grid3D_AccelInterp = 2     ! Output accel, linear interp
         end if
      else 
         if (FF%VelInterpCubic) then
            Grid3D_AccelInterp = 3     ! No accel, cubic interp
         else
            Grid3D_AccelInterp = 4     ! No accel, linear interp
         end if
      end if

      ! Store flag value since it doesn't change during loop
      AddMeanAfterInterp = FF%Grid3D%AddMeanAfterInterp

      ! Loop through points
      do i = 1, NumPoints

         ! If height < zero, set velocity/acceleration to zero, continue
         if (Position(3, i) <= 0.0_ReKi) then
            VelocityUVW(:, i) = 0.0_ReKi
            if (OutputAccel) AccelUVW(:, i) = 0.0_ReKi
            cycle
         end if

         ! Is this point allowed beyond the bounds of the wind box?
         GridExceedAllow = FF%Grid3D%BoxExceedAllowF .and. (i >= FF%Grid3D%BoxExceedAllowIdx)

         ! Calculate grid cells for interpolation, returns velocity and acceleration
         ! components at corners of grid cell containing time and position. Also
         ! returns interpolation values Xi.
         call Grid3DField_GetCell(FF%Grid3D, Time, Position(:, i), OutputAccel, GridExceedAllow, &
                                  VelCell, AccCell, Xi, Is3D, TmpErrStat, TmpErrMsg)
         if (TmpErrStat >= AbortErrLev) then
            call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            return
         end if

         ! Switch based on if acceleration is output and velocity interpolation is cubic
         select case(Grid3D_AccelInterp)
         case (1)    ! Cubic velocity and cubic acceleration
            call Grid3DField_GetVelAccCubic(FF%Grid3D%DTime, VelCell, AccCell, Xi, Is3D, &
                                            Velocity=VelocityUVW(:, i), Accel=AccelUVW(:, i))
         case (2)    ! Linear velocity and cubic acceleration
            VelocityUVW(:, i) = Grid3DField_GetVelLinear(VelCell, Xi, Is3D)
            call Grid3DField_GetVelAccCubic(FF%Grid3D%DTime, VelCell, AccCell, Xi, Is3D, &
                                            Accel=AccelUVW(:, i))
         case (3)    ! Cubic velocity and no acceleration
            call Grid3DField_GetVelAccCubic(FF%Grid3D%DTime, VelCell, AccCell, Xi, Is3D, &
                                            Velocity=VelocityUVW(:, i))
         case (4)    ! Linear velocity and no acceleration
            VelocityUVW(:, i) = Grid3DField_GetVelLinear(VelCell, Xi, Is3D)
         end select

         ! Add mean wind speed after interpolation if flag is set
         if (AddMeanAfterInterp) then
            VelocityUVW(1, i) = VelocityUVW(1, i) + &
                                CalculateMeanVelocity(FF%Grid3D, Position(3, i), Position(2, i))
         end if

      end do

   case (Grid4D_FieldType)

      !-------------------------------------------------------------------------
      ! Grid4D Flow Field
      !-------------------------------------------------------------------------

      ! If field is not allocated, return error
      if (.not. allocated(FF%Grid4D%Vel)) then
         call SetErrStat(ErrID_Fatal, "Grid4D Field not allocated", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Loop through points
      do i = 1, NumPoints

         ! If height less than or equal to zero, set velocity to zero
         if (Position(3, i) <= 0.0_ReKi) then
            VelocityUVW(:, i) = 0.0_ReKi
            cycle
         end if

         call Grid4DField_GetVel(FF%Grid4D, Time, Position(:, i), VelocityUVW(:, i), TmpErrStat, TmpErrMsg)
         if (TmpErrStat >= AbortErrLev) then
            call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            return
         end if
      end do

   case (Point_FieldType)

      !-------------------------------------------------------------------------
      ! Point Flow Field
      !-------------------------------------------------------------------------

      ! If points field is not allocated, return error
      if (.not. allocated(FF%Points%Vel)) then
         call SetErrStat(ErrID_Fatal, "Points Point Field not allocated", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Set velocities directly from velocity array
      VelocityUVW = FF%Points%Vel(:, IStart:IStart + NumPoints - 1)

   case (User_FieldType)

      !-------------------------------------------------------------------------
      ! User Flow Field
      !-------------------------------------------------------------------------

      call SetErrStat(ErrID_Fatal, "User Field not to be implemented", ErrStat, ErrMsg, RoutineName)
      return

   case default
      call SetErrStat(ErrID_Fatal, "Invalid FieldType "//trim(num2lstr(FF%FieldType)), ErrStat, ErrMsg, RoutineName)
      return
   end select

   !----------------------------------------------------------------------------
   ! If wind box was rotated, apply rotation to velocity/acceleration
   !----------------------------------------------------------------------------

   if (FF%RotateWindBox) then
      if (.not. OutputAccel) then
         do i = 1, NumPoints
            VelocityUVW(:, i) = matmul(FF%RotFromWind, VelocityUVW(:, i))
         end do
      else
         do i = 1, NumPoints
            VelocityUVW(:, i) = matmul(FF%RotFromWind, VelocityUVW(:, i))
            AccelUVW(:, i) = matmul(FF%RotFromWind, AccelUVW(:, i))
         end do
      end if
   end if

contains

   pure function GetPrimePosition(Pos) result(PrimePos)
      real(ReKi), dimension(3), intent(in)      :: Pos
      real(ReKi), dimension(3)                  :: PrimePos
      PrimePos = matmul(FF%RotToWind, (Pos - FF%RefPosition)) + FF%RefPosition
   end function

end subroutine

pure function UniformField_GetVel(UF, op, Position) result(Velocity)
   use, intrinsic :: ieee_arithmetic

   type(UniformFieldType), intent(in)        :: UF
   type(UniformField_Interp), intent(in)     :: op
   real(ReKi), dimension(3), intent(in)      :: Position
   real(ReKi), dimension(3)                  :: Velocity

   character(*), parameter                   :: RoutineName = "UniformField_GetVel"
   real(ReKi)                                :: V1
   real(ReKi)                                :: V1_rotate, VZ_rotate

   ! Calculate horizontal velocity if position is above ground
   V1 = op%VelH*((Position(3)/UF%RefHeight)**op%ShrV &                                 ! power-law wind shear
                 + (op%ShrH*(Position(2)*op%CosAngleH + Position(1)*op%SinAngleH) & ! horizontal linear shear
                    + op%LinShrV*(Position(3) - UF%RefHeight))/UF%RefLength) &         ! vertical linear shear
        + op%VelGust                                                                      ! gust speed

   ! Apply upflow angle
   V1_rotate = op%CosAngleV*V1 - op%SinAngleV*op%VelV
   VZ_rotate = op%SinAngleV*V1 + op%CosAngleV*op%VelV

   ! Apply wind direction
   Velocity = [V1_rotate*op%CosAngleH, -V1_rotate*op%SinAngleH, VZ_rotate]

end function

function UniformField_GetAcc(UF, op, Position) result(Accel)
   type(UniformFieldType), intent(in)     :: UF
   type(UniformField_Interp), intent(in)  :: op
   real(ReKi), dimension(3), intent(in)   :: Position
   real(ReKi), dimension(3)               :: Accel

   character(*), parameter                :: RoutineName = "UniformField_GetAcc"
   real(ReKi)                             :: C1, C2, C3, C4, C5

   C1 = (Position(3)/UF%RefHeight)**op%ShrV + &
        (op%LinShrV*(Position(3) - UF%RefHeight) + &
         op%ShrH*(Position(1)*op%SinAngleH + &
                  Position(2)*op%CosAngleH))/UF%RefLength

   C2 = op%CosAngleV*(op%VelGust + op%VelH*(C1))

   C3 = op%AngleVDot*op%SinAngleV*(op%VelGust + op%VelH*(C1))

   C4 = op%LinShrVDot*(Position(3) - UF%RefHeight) + &
        op%ShrHDot*(Position(1)*op%SinAngleH + Position(2)*op%CosAngleH) + &
        op%ShrH*(Position(1)*op%AngleHDot*op%CosAngleH - &
                 Position(2)*op%AngleHDot*op%SinAngleH)

   C5 = op%VelGustDot + op%VelHDot*C1 + &
        op%VelH*(op%ShrVDot*(Position(3)/UF%RefHeight)**op%ShrV* &
                 log(Position(3)/UF%RefHeight) + C4/UF%RefLength)

   Accel(1) = -op%AngleHDot*op%SinAngleH*(C2 - op%SinAngleV*op%VelV) + &
              op%CosAngleH*(-op%AngleVDot*op%CosAngleV*op%VelV - C3 - &
                            op%VelVDot*op%SinAngleV + op%CosAngleV*C5)

   Accel(2) = op%AngleHDot*op%CosAngleH*(-C2 + op%SinAngleV*op%VelV) + &
              op%SinAngleH*(op%AngleVDot*op%CosAngleV*op%VelV + C3 + &
                            op%VelVDot*op%SinAngleV - op%CosAngleV*C5)

   Accel(3) = op%AngleVDot*C2 - op%AngleVDot*op%SinAngleV*op%VelV + &
              op%VelVDot*op%CosAngleV + op%SinAngleV*C5

end function

pure function UniformField_InterpLinear(UF, Time) result(op)

   type(UniformFieldType), intent(in)  :: UF
   real(DbKi), intent(in)              :: Time
   type(UniformField_Interp)           :: op
   integer(IntKi)                      :: i
   real(ReKi)                          :: dt, alpha, OMalpha

   ! If only one data point or time is at or below first time, use first sample
   if (UF%DataSize == 1 .or. Time < UF%Time(1)) then

      op%VelH = UF%VelH(1)
      op%AngleH = UF%AngleH(1)
      op%AngleV = UF%AngleV(1)
      op%VelV = UF%VelV(1)
      op%ShrH = UF%ShrH(1)
      op%ShrV = UF%ShrV(1)
      op%LinShrV = UF%LinShrV(1)
      op%VelGust = UF%VelGust(1)

      ! If time is after end time, use last data point
   else if (Time >= UF%Time(UF%DataSize)) then

      op%VelH = UF%VelH(UF%DataSize)
      op%AngleH = UF%AngleH(UF%DataSize)
      op%AngleV = UF%AngleV(UF%DataSize)
      op%VelV = UF%VelV(UF%DataSize)
      op%ShrH = UF%ShrH(UF%DataSize)
      op%ShrV = UF%ShrV(UF%DataSize)
      op%LinShrV = UF%LinShrV(UF%DataSize)
      op%VelGust = UF%VelGust(UF%DataSize)

   else

      ! Find first index where current time is less than Time(i)
      do i = 2, UF%DataSize
         if (Time < UF%Time(i)) exit
      end do

      ! Calculate interval delta time
      dt = UF%Time(i) - UF%Time(i - 1)

      ! Calculate interpolation coefficient [0,1]
      alpha = real((Time - UF%Time(i - 1))/dt, ReKi)
      OMalpha = 1.0_ReKi - alpha

      ! Blend states before and after time based on alpha
      op%VelH = UF%VelH(i - 1)*OMalpha + UF%VelH(i)*alpha
      op%AngleH = UF%AngleH(i - 1)*OMalpha + UF%AngleH(i)*alpha
      op%AngleV = UF%AngleV(i - 1)*OMalpha + UF%AngleV(i)*alpha
      op%VelV = UF%VelV(i - 1)*OMalpha + UF%VelV(i)*alpha
      op%ShrH = UF%ShrH(i - 1)*OMalpha + UF%ShrH(i)*alpha
      op%ShrV = UF%ShrV(i - 1)*OMalpha + UF%ShrV(i)*alpha
      op%LinShrV = UF%LinShrV(i - 1)*OMalpha + UF%LinShrV(i)*alpha
      op%VelGust = UF%VelGust(i - 1)*OMalpha + UF%VelGust(i)*alpha

   end if

   op%CosAngleH = cos(op%AngleH)
   op%SinAngleH = sin(op%AngleH)
   op%CosAngleV = cos(op%AngleV)
   op%SinAngleV = sin(op%AngleV)

end function

pure function UniformField_InterpCubic(UF, Time) result(op)

   type(UniformFieldType), intent(in)  :: UF
   real(DbKi), intent(in)              :: Time
   type(UniformField_Interp)           :: op

   integer(IntKi)                      :: i
   real(ReKi)                          :: C1, C2, C3, C4, h, t

   ! Initialize data index
   i = 0

   ! If time is outside of array
   if (UF%DataSize == 1 .or. Time < UF%Time(1)) then
      ! One data point or time is at or below first time, use first sample
      i = 1
   else if (Time >= UF%Time(UF%DataSize)) then
      ! Time is after end time, use last data point
      i = UF%DataSize
   end if

   ! If time was inside data array
   if (i == 0) then

      ! Find first index where current time is less than Time(i)
      do i = 2, UF%DataSize
         if (Time < UF%Time(i)) exit
      end do

      h = UF%Time(i) - UF%Time(i - 1)
      t = real((Time - UF%Time(i - 1))/h, ReKi)

      C1 = 2.0_ReKi*t*t*t - 3.0_ReKi*t*t + 1.0_ReKi
      C2 = (t*t*t - 2.0_ReKi*t*t + t)*h
      C3 = -2.0_ReKi*t*t*t + 3.0_ReKi*t*t
      C4 = (t*t*t - t*t)*h

      op%VelH = C1*UF%VelH(i - 1) + C2*UF%VelHDot(i - 1) + C3*UF%VelH(i) + C4*UF%VelHDot(i)
      op%AngleH = C1*UF%AngleH(i - 1) + C2*UF%AngleHDot(i - 1) + C3*UF%AngleH(i) + C4*UF%AngleHDot(i)
      op%AngleV = C1*UF%AngleV(i - 1) + C2*UF%AngleVDot(i - 1) + C3*UF%AngleV(i) + C4*UF%AngleVDot(i)
      op%VelV = C1*UF%VelV(i - 1) + C2*UF%VelVDot(i - 1) + C3*UF%VelV(i) + C4*UF%VelVDot(i)
      op%ShrH = C1*UF%ShrH(i - 1) + C2*UF%ShrHDot(i - 1) + C3*UF%ShrH(i) + C4*UF%ShrHDot(i)
      op%ShrV = C1*UF%ShrV(i - 1) + C2*UF%ShrVDot(i - 1) + C3*UF%ShrV(i) + C4*UF%ShrVDot(i)
      op%LinShrV = C1*UF%LinShrV(i - 1) + C2*UF%LinShrVDot(i - 1) + C3*UF%LinShrV(i) + C4*UF%LinShrVDot(i)
      op%VelGust = C1*UF%VelGust(i - 1) + C2*UF%VelGustDot(i - 1) + C3*UF%VelGust(i) + C4*UF%VelGustDot(i)

      C1 = (6.0_ReKi*t*t - 6.0_ReKi*t)/h
      C2 = (3.0_ReKi*t*t - 4.0_ReKi*t + 1.0_ReKi)
      C3 = -C1
      C4 = (3.0_ReKi*t*t - 2.0_ReKi*t)

      op%VelHDot = C1*UF%VelH(i - 1) + C2*UF%VelHDot(i - 1) + C3*UF%VelH(i) + C4*UF%VelHDot(i)
      op%AngleHDot = C1*UF%AngleH(i - 1) + C2*UF%AngleHDot(i - 1) + C3*UF%AngleH(i) + C4*UF%AngleHDot(i)
      op%AngleVDot = C1*UF%AngleV(i - 1) + C2*UF%AngleVDot(i - 1) + C3*UF%AngleV(i) + C4*UF%AngleVDot(i)
      op%VelVDot = C1*UF%VelV(i - 1) + C2*UF%VelVDot(i - 1) + C3*UF%VelV(i) + C4*UF%VelVDot(i)
      op%ShrHDot = C1*UF%ShrH(i - 1) + C2*UF%ShrHDot(i - 1) + C3*UF%ShrH(i) + C4*UF%ShrHDot(i)
      op%ShrVDot = C1*UF%ShrV(i - 1) + C2*UF%ShrVDot(i - 1) + C3*UF%ShrV(i) + C4*UF%ShrVDot(i)
      op%LinShrVDot = C1*UF%LinShrV(i - 1) + C2*UF%LinShrVDot(i - 1) + C3*UF%LinShrV(i) + C4*UF%LinShrVDot(i)
      op%VelGustDot = C1*UF%VelGust(i - 1) + C2*UF%VelGustDot(i - 1) + C3*UF%VelGust(i) + C4*UF%VelGustDot(i)

   else

      ! Set values based on first/last index
      op%VelH = UF%VelH(i)
      op%AngleH = UF%AngleH(i)
      op%AngleV = UF%AngleV(i)
      op%VelV = UF%VelV(i)
      op%ShrH = UF%ShrH(i)
      op%ShrV = UF%ShrV(i)
      op%LinShrV = UF%LinShrV(i)
      op%VelGust = UF%VelGust(i)

      op%VelHDot = 0.0_ReKi
      op%AngleHDot = 0.0_ReKi
      op%AngleVDot = 0.0_ReKi
      op%VelVDot = 0.0_ReKi
      op%ShrHDot = 0.0_ReKi
      op%ShrVDot = 0.0_ReKi
      op%LinShrVDot = 0.0_ReKi
      op%VelGustDot = 0.0_ReKi

   end if

   op%CosAngleH = cos(op%AngleH)
   op%SinAngleH = sin(op%AngleH)
   op%CosAngleV = cos(op%AngleV)
   op%SinAngleV = sin(op%AngleV)

end function

subroutine IfW_UniformField_CalcAccel(UF, ErrStat, ErrMsg)
   type(UniformFieldType), intent(inout)  :: UF
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg

   character(*), parameter             :: RoutineName = "Uniform_CalcAccel"
   integer(IntKi)                      :: TmpErrStat
   character(ErrMsgLen)                :: TmpErrMsg
   real(ReKi), allocatable             :: b(:), u(:), dy2(:)

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Storage for spline fit arrays
   !----------------------------------------------------------------------------

   call AllocAry(B, UF%DataSize, "storage for B", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   call AllocAry(U, UF%DataSize, "storage for U", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   call AllocAry(dy2, UF%DataSize, "storage for dy2", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Storage for derivative arrays
   !----------------------------------------------------------------------------

   call AllocAry(UF%VelHDot, UF%DataSize, 'Uniform wind horizontal wind speed derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%AngleHDot, UF%DataSize, 'Uniform wind direction derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%AngleVDot, UF%DataSize, 'Uniform wind upflow angle derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%VelVDot, UF%DataSize, 'Uniform vertical wind speed derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%ShrHDot, UF%DataSize, 'Uniform horizontal linear shear derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%ShrVDot, UF%DataSize, 'Uniform vertical power-law shear exponent derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%LinShrVDot, UF%DataSize, 'Uniform vertical linear shear derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   call AllocAry(UF%VelGustDot, UF%DataSize, 'Uniform gust velocity derivative', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat > AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Calculate derivatives
   !----------------------------------------------------------------------------

   call CalcCubicSplineDeriv(UF%Time, UF%VelH, UF%VelHDot)
   call CalcCubicSplineDeriv(UF%Time, UF%AngleH, UF%AngleHDot)
   call CalcCubicSplineDeriv(UF%Time, UF%AngleV, UF%AngleVDot)
   call CalcCubicSplineDeriv(UF%Time, UF%VelV, UF%VelVDot)
   call CalcCubicSplineDeriv(UF%Time, UF%ShrH, UF%ShrHDot)
   call CalcCubicSplineDeriv(UF%Time, UF%ShrV, UF%ShrVDot)
   call CalcCubicSplineDeriv(UF%Time, UF%LinShrV, UF%LinShrVDot)
   call CalcCubicSplineDeriv(UF%Time, UF%VelGust, UF%VelGustDot)

contains

   !> CalcCubicSplineDeriv fits a cubic spline through the y array with points
   !! at x locations. It then calculates the corresponding derivative of y
   !! with respect to x at the same x values and returns it in the dy array.
   subroutine CalcCubicSplineDeriv(x, y, dy)
      real(ReKi), intent(in)        :: x(:)
      real(ReKi), intent(in)        :: y(:)
      real(ReKi), intent(out)       :: dy(:)

      integer(IntKi)                :: i, n
      real(ReKi)                    :: p, sig, un

      ! Get size of arrays
      n = size(x)

      ! If 1 or 2 points, set derivatives to zero and return
      if (n < 3) then
         do i = 1, n
            dy(i) = 0.0_ReKi
         end do
         return
      end if

      ! Natural lower and upper boundaries (second derivative = 0)
      ! u(1) = 0.0_ReKi
      ! dy2(1) = 0.0_ReKi
      ! dy2(n) = 0.0_ReKi

      ! First derivative is zero at lower boundary condition
      dy2(1) = -0.5_ReKi
      u(1) = 3.0_ReKi*(y(2) - y(1))/(x(2) - x(1))**2

      ! Calculate slopes
      do i = 1, n - 1
         b(i) = (y(i + 1) - y(i))/(x(i + 1) - x(i))
      end do

      ! Decomposition
      do i = 2, n - 1
         sig = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1))
         p = sig*dy2(i - 1) + 2.0_ReKi
         dy2(i) = (sig - 1.0_ReKi)/p
         u(i) = (6.*((y(i + 1) - y(i))/(x(i + 1) - x(i)) - (y(i) - y(i - 1))/(x(i) - x(i - 1)))/ &
                 (x(i + 1) - x(i - 1)) - sig*u(i - 1))/p
      end do

      ! First derviative is zero at upper boundary condition
      un = -3.0_ReKi*(y(n) - y(n - 1))/(x(n) - x(n - 1))**2
      dy2(n) = (un - 0.5_ReKi*u(n - 1))/(0.5_ReKi*dy2(n - 1) + 1.0_ReKi)

      ! Back substitution and derivative calculation
      do i = n - 1, 1, -1
         dy2(i) = dy2(i)*dy2(i + 1) + u(i)
         dy(i) = real(b(i) - (x(i + 1) - x(i))*(dy2(i)/3.0_ReKi + dy2(i + 1)/6.0_ReKi), SiKi)
      end do
      dy(n) = 0.0_ReKi

   end subroutine

end subroutine

!> Routine to compute the Jacobians of the output (Y) function with respect to the inputs (u). The partial
!! derivative dY/du is returned. This submodule does not follow the modularization framework.
subroutine IfW_UniformWind_GetOP(UF, t, InterpCubic, OP_out)
   type(UniformFieldType), intent(IN)  :: UF             !< Parameters
   real(DbKi), intent(IN)              :: t              !< Current simulation time in seconds
   logical, intent(in)                 :: InterpCubic    !< flag for using cubic interpolation
   real(ReKi), intent(OUT)             :: OP_out(2)      !< operating point (HWindSpeed and PLexp

   type(UniformField_Interp)           :: op         ! interpolated values of InterpParams

   ! Linearly interpolate parameters in time at operating point (or use nearest-neighbor to extrapolate)
   if (InterpCubic) then
      op = UniformField_InterpCubic(UF, t)
   else
      op = UniformField_InterpLinear(UF, t)
   end if

   OP_out(1) = op%VelH
   OP_out(2) = op%ShrV

end subroutine

subroutine Grid3DField_GetCell(G3D, Time, Position, CalcAccel, AllowExtrap, &
                               VelCell, AccCell, Xi, Is3D, ErrStat, ErrMsg)

   type(Grid3DFieldType), intent(in)   :: G3D               !< 3D Grid-Field data
   real(DbKi), intent(in)              :: Time              !< time (s)
   real(ReKi), intent(in)              :: Position(3)       !< position X,Y,Z to get value
   logical, intent(in)                 :: CalcAccel         !< flag to populat AccCell
   logical, intent(in)                 :: AllowExtrap       !< is this point allowed to exceed bounds of wind grid
   real(ReKi), intent(out)             :: VelCell(8, 3)     !< Velocity components at corners of grid cell
   real(ReKi), intent(out)             :: AccCell(8, 3)     !< Acceleration components at corners of grid cell
   real(ReKi), intent(out)             :: Xi(3)             !< isoparametric coord of position in cell (y,z,t) [-1, +1]
   logical, intent(out)                :: Is3D              !< flag indicating if interpolation is 3D or 2D
   integer(IntKi), intent(out)         :: ErrStat           !< error status
   character(*), intent(out)           :: ErrMsg            !< error message

   character(*), parameter             :: RoutineName = "Grid3DField_GetCell"
   integer(IntKi), parameter           :: ExtrapNone = 0
   integer(IntKi), parameter           :: ExtrapYmin = 1
   integer(IntKi), parameter           :: ExtrapYmax = 2
   integer(IntKi), parameter           :: ExtrapZmin = 4
   integer(IntKi), parameter           :: ExtrapZmax = 8
   integer(IntKi)                      :: AllExtrap
   integer(IntKi)                      :: IY_Lo, IY_Hi
   integer(IntKi)                      :: IZ_Lo, IZ_Hi
   integer(IntKi)                      :: IT_Lo, IT_Hi
   logical                             :: InGrid

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Initialize to no extrapolation (modified in bounds routines)
   AllExtrap = ExtrapNone

   !----------------------------------------------------------------------------
   ! Find grid bounds in Time and Z
   !----------------------------------------------------------------------------

   ! Get grid time bounds
   call GetBoundsT(Position(1), Xi(3))
   if (ErrStat >= AbortErrLev) return

   ! Get grid Z bounds
   call GetBoundsZ(Position(3), Xi(2))
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Extract interpolation cells from grids based on poisiont
   !----------------------------------------------------------------------------

   ! If position is inside grid
   if (InGrid) then

      ! Set flag to use 3D interpolation
      Is3D = .true.

      ! Get grid Y bounds
      call GetBoundsY(Position(2), Xi(1))
      if (ErrStat >= AbortErrLev) return

      ! Interpolate within grid (or top, left, right if extrapolation enabled)
      call GetCellInGrid(VelCell, G3D%Vel, G3D%VelAvg)

      ! If acceleration requested, get cell values
      if (CalcAccel) then
         call GetCellInGrid(AccCell, G3D%Acc, G3D%AccAvg)
      end if

   else if (G3D%NTGrids > 0) then

      ! Interpolation is 2D
      Is3D = .false.

      ! Tower grids present and position is below main grid
      call GetCellInTower(VelCell, G3D%Vel, G3D%VelAvg, G3D%VelTower)

      ! If acceleration requested, get cell values
      if (CalcAccel) then
         call GetCellInTower(AccCell, G3D%Acc, G3D%AccAvg, G3D%AccTower)
      end if

   else

      ! Set flag to use 3D interpolation
      Is3D = .true.

      ! Get grid Y bounds
      call GetBoundsY(Position(2), Xi(1))
      if (ErrStat >= AbortErrLev) return

      ! Tower interpolation without tower grids
      call GetCellBelowGrid(VelCell, G3D%Vel)

      ! If acceleration requested, get cell values
      if (CalcAccel) then
         call GetCellBelowGrid(AccCell, G3D%Acc)
      end if

   end if

contains

   subroutine GetCellInGrid(cell, gridVal, gridAvg)

      real(ReKi), intent(out)             :: cell(8, 3)
      real(SiKi), intent(in)              :: gridVal(:, :, :, :)
      real(SiKi), intent(in), allocatable :: gridAvg(:, :, :)

      ! Select based on extrapolation flags
      select case (AllExtrap)

      case (ExtrapNone)                   ! No extrapolation

         cell(1, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Lo)
         cell(2, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Lo)
         cell(3, :) = gridVal(:, IY_Lo, IZ_Hi, IT_Lo)
         cell(4, :) = gridVal(:, IY_Hi, IZ_Hi, IT_Lo)
         cell(5, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Hi)
         cell(6, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Hi)
         cell(7, :) = gridVal(:, IY_Lo, IZ_Hi, IT_Hi)
         cell(8, :) = gridVal(:, IY_Hi, IZ_Hi, IT_Hi)

      case (ior(ExtrapZmax, ExtrapYmax))   ! Extrapolate top right corner

         cell(1, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Lo)
         cell(2, :) = gridAvg(:, IZ_Lo, IT_Lo)
         cell(3, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(4, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(5, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Hi)
         cell(6, :) = gridAvg(:, IZ_Lo, IT_Hi)
         cell(7, :) = gridAvg(:, IZ_Hi, IT_Hi)
         cell(8, :) = gridAvg(:, IZ_Hi, IT_Hi)

      case (ior(ExtrapZmax, ExtrapYmin))! Extrapolate top left corner

         cell(1, :) = gridAvg(:, IZ_Lo, IT_Lo)
         cell(2, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Lo)
         cell(3, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(4, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(5, :) = gridAvg(:, IZ_Lo, IT_Hi)
         cell(6, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Hi)
         cell(7, :) = gridAvg(:, IZ_Hi, IT_Hi)
         cell(8, :) = gridAvg(:, IZ_Hi, IT_Hi)

      case (ExtrapZmax)   ! Extrapolate above grid only

         cell(1, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Lo)
         cell(2, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Lo)
         cell(3, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(4, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(5, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Hi)
         cell(6, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Hi)
         cell(7, :) = gridAvg(:, IZ_Hi, IT_Hi)
         cell(8, :) = gridAvg(:, IZ_Hi, IT_Hi)

      case (ExtrapYmax)   ! Extrapolate to the right of grid only

         cell(1, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Lo)
         cell(2, :) = gridAvg(:, IZ_Lo, IT_Lo)
         cell(3, :) = gridVal(:, IY_Lo, IZ_Hi, IT_Lo)
         cell(4, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(5, :) = gridVal(:, IY_Lo, IZ_Lo, IT_Hi)
         cell(6, :) = gridAvg(:, IZ_Lo, IT_Hi)
         cell(7, :) = gridVal(:, IY_Lo, IZ_Hi, IT_Hi)
         cell(8, :) = gridAvg(:, IZ_Hi, IT_Hi)

      case (ExtrapYmin)   ! Extrapolate to the left of grid only

         cell(1, :) = gridAvg(:, IZ_Lo, IT_Lo)
         cell(2, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Lo)
         cell(3, :) = gridAvg(:, IZ_Hi, IT_Lo)
         cell(4, :) = gridVal(:, IY_Hi, IZ_Hi, IT_Lo)
         cell(5, :) = gridAvg(:, IZ_Lo, IT_Hi)
         cell(6, :) = gridVal(:, IY_Hi, IZ_Lo, IT_Hi)
         cell(7, :) = gridAvg(:, IZ_Hi, IT_Hi)
         cell(8, :) = gridVal(:, IY_Hi, IZ_Hi, IT_Hi)

      case (ExtrapZmin)   ! Extrapolate below grid

         cell(1, :) = 0.0_ReKi                        ! Ground
         cell(2, :) = 0.0_ReKi                        ! Ground
         cell(3, :) = gridVal(:, IY_Lo, 1, IT_Lo)
         cell(4, :) = gridVal(:, IY_Hi, 1, IT_Lo)
         cell(5, :) = 0.0_ReKi                        ! Ground
         cell(6, :) = 0.0_ReKi                        ! Ground
         cell(7, :) = gridVal(:, IY_Lo, 1, IT_Hi)
         cell(8, :) = gridVal(:, IY_Hi, 1, IT_Hi)

      case (ior(ExtrapZmin, ExtrapYmin))   ! Extrapolate lower left of grid

         cell(1, :) = 0.0_ReKi                        ! Ground
         cell(2, :) = 0.0_ReKi                        ! Ground
         cell(3, :) = gridAvg(:, 1, IT_Lo)            ! Average
         cell(4, :) = gridVal(:, 1, 1, IT_Lo)
         cell(5, :) = 0.0_ReKi                        ! Ground
         cell(6, :) = 0.0_ReKi                        ! Ground
         cell(7, :) = gridAvg(:, 1, IT_Hi)            ! Average
         cell(8, :) = gridVal(:, 1, 1, IT_Hi)

      case (ior(ExtrapZmin, ExtrapYmax))   ! Extrapolate lower right of grid

         cell(1, :) = 0.0_ReKi                        ! Ground
         cell(2, :) = 0.0_ReKi                        ! Ground
         cell(3, :) = gridVal(:, G3D%NYGrids, 1, IT_Lo)
         cell(4, :) = gridAvg(:, 1, IT_Lo)            ! Average
         cell(5, :) = 0.0_ReKi                        ! Ground
         cell(6, :) = 0.0_ReKi                        ! Ground
         cell(7, :) = gridVal(:, G3D%NYGrids, 1, IT_Hi)
         cell(8, :) = gridAvg(:, 1, IT_Hi)            ! Average

      end select

   end subroutine

   !> GetCellBelowGrid interpolates between bottom of grid and ground. This
   !! is only called if G3D%InterpTower == .true.
   subroutine GetCellBelowGrid(cell, gridVal)

      real(ReKi), intent(out) :: cell(8, 3)
      real(SiKi), intent(in)  :: gridVal(:, :, :, :)

      cell(1, :) = 0.0_ReKi                           ! Ground
      cell(2, :) = 0.0_ReKi                           ! Ground
      cell(3, :) = gridVal(:, IY_Lo, IZ_Hi, IT_Lo)
      cell(4, :) = gridVal(:, IY_Hi, IZ_Hi, IT_Lo)
      cell(5, :) = 0.0_ReKi                           ! Ground
      cell(6, :) = 0.0_ReKi                           ! Ground
      cell(7, :) = gridVal(:, IY_Lo, IZ_Hi, IT_Hi)
      cell(8, :) = gridVal(:, IY_Hi, IZ_Hi, IT_Hi)

   end subroutine

   subroutine GetCellInTower(cell, gridVal, gridAvg, towerVal)

      real(ReKi), intent(out)             :: cell(8, 3)
      real(SiKi), intent(in)              :: gridVal(:, :, :, :)
      real(SiKi), intent(in), allocatable :: gridAvg(:, :, :)
      real(SiKi), intent(in), allocatable :: towerVal(:, :, :)

      real(ReKi), dimension(2)            :: P, P1, P2, P3, V0, V1, V2
      real(ReKi)                          :: d00, d01, d11, d20, d21
      real(ReKi)                          :: V(3, 3, 2), W(3)
      real(ReKi)                          :: alpha, omalpha, denom
      integer(IntKi)                      :: ic

      !-------------------------------------------------------------------------
      ! If extrapolation is not allowed or Y is nearly zero, only interpolate
      ! along the tower
      !-------------------------------------------------------------------------

      if (.not. AllowExtrap) then
         if (IZ_HI <= G3D%NTGrids) then      ! In tower grid
            cell(1, :) = towerVal(:, IZ_LO, IT_LO)
            cell(2, :) = towerVal(:, IZ_HI, IT_LO)
            cell(3, :) = towerVal(:, IZ_LO, IT_HI)
            cell(4, :) = towerVal(:, IZ_HI, IT_HI)
         else                                ! Between tower grid and ground
            cell(1, :) = towerVal(:, IZ_LO, IT_LO)
            cell(2, :) = 0.0_ReKi
            cell(3, :) = towerVal(:, IZ_LO, IT_HI)
            cell(4, :) = 0.0_ReKi
         end if
         return
      end if

      !-------------------------------------------------------------------------
      ! If Y is beyond grid width from the tower (clamped),
      ! interp between ground and bottom of grid average
      !-------------------------------------------------------------------------

      if (abs(Position(2)) >= 2.0_ReKi*G3D%YHWid) then
         Xi(2) = 2.0_ReKi*Position(3)/G3D%GridBase - 1.0_ReKi
         cell(1, :) = 0.0_ReKi
         cell(2, :) = gridAvg(:, 1, IT_LO)
         cell(3, :) = 0.0_ReKi
         cell(4, :) = gridAvg(:, 1, IT_HI)
         return
      end if

      !-------------------------------------------------------------------------
      ! Otherwise, position is below grid and within +- 2*GridWidth from tower
      ! This section uses Barycentric interpolation of a triangle to get
      ! the wind components at the desired Position. The components on the 
      ! bottom of the grid and on the tower are interpolated to get the first
      ! two points. The third point is on the ground at the GridWidth away
      ! from the tower.
      !-------------------------------------------------------------------------

      ! Get grid Y bounds
      call GetBoundsY(Position(2), Xi(1))
      if (ErrStat >= AbortErrLev) return

      ! Get interpolation point
      P = [abs(Position(2)), Position(3)]

      ! Point 1 (grid bottom point)
      P1 = [abs(Position(2)), G3D%GridBase]
      select case (AllExtrap)
      case (ExtrapNone)
         ! Interpolate between grid points
         alpha = (Xi(1) + 1.0_ReKi)/2.0_ReKi
         omalpha = 1.0_ReKi - alpha
         V(:, 1, 1) = omalpha*gridVal(:, IY_Lo, 1, IT_Lo) + &
                      alpha*gridVal(:, IY_Hi, 1, IT_Lo)
         V(:, 1, 2) = omalpha*gridVal(:, IY_Lo, 1, IT_Hi) + &
                      alpha*gridVal(:, IY_Hi, 1, IT_Hi)
      case (ExtrapYmin, ExtrapYmax)
         ! Interpolate between edge of grid and grid average
         alpha = abs(Position(2))/G3D%YHWid - 1.0_ReKi
         omalpha = 1.0_ReKi - alpha
         V(:, 1, 1) = omalpha*gridVal(:, IY_Lo, 1, IT_Lo) + &
                      alpha*gridAvg(:, 1, IT_Lo)
         V(:, 1, 2) = omalpha*gridVal(:, IY_Lo, 1, IT_Hi) + &
                      alpha*gridAvg(:, 1, IT_Hi)
      end select

      ! Point 2 (tower point)
      P2 = [0.0_ReKi, Position(3)]
      alpha = (Xi(2) + 1.0_ReKi)/2.0_ReKi
      omalpha = 1.0_ReKi - alpha
      if (IZ_HI <= G3D%NTGrids) then   ! Lower point above ground
         V(:, 2, 1) = omalpha*towerVal(:, IZ_Lo, IT_Lo) + &
                      alpha*towerVal(:, IZ_Hi, IT_Lo)
         V(:, 2, 2) = omalpha*towerVal(:, IZ_Lo, IT_Hi) + &
                      alpha*towerVal(:, IZ_Hi, IT_Hi)
      else                             ! Lower point on ground
         V(:, 2, 1) = omalpha*towerVal(:, IZ_Lo, IT_Lo)
         V(:, 2, 2) = omalpha*towerVal(:, IZ_Lo, IT_Hi)
      end if

      ! Point 3 (ground @ grid width away from tower)
      P3 = [2.0_ReKi*G3D%YHWid, 0.0_Reki]
      ! V(:, 3, :) = 0.0_ReKi ! Not used

      ! Calculate Barycentric weights for triangle
      V0 = P1 - P3
      V1 = P2 - P3
      V2 = P - P3
      d00 = dot_product(v0, v0)
      d01 = dot_product(v0, v1)
      d11 = dot_product(v1, v1)
      d20 = dot_product(v2, v0)
      d21 = dot_product(v2, v1)
      denom = d00*d11 - d01*d01
      W(1) = (d11*d20 - d01*d21)/denom
      W(2) = (d00*d21 - d01*d20)/denom
      ! W(3) = 1.0_ReKi - W(1) - W(2) ! Not used

      ! Interpolate wind components based on weights
      do ic = 1, 3
         cell(1, ic) = V(ic, 1, 1) * W(1) + V(ic, 2, 1) * W(2)
         cell(3, ic) = V(ic, 1, 2) * W(1) + V(ic, 2, 2) * W(2)
      end do
      cell(2, :) = cell(1, :)
      cell(4, :) = cell(3, :)

   end subroutine

   !> GetBoundsY populates IY_Lo, IY_Hi, and the interpolant [-1,1]. It also
   !! adds ExtrapYmin or ExtrapYmax to AllExtrap if applicable.
   subroutine GetBoundsY(PosY, DY)

      real(ReKi), intent(in)     :: PosY
      real(ReKi), intent(out)    :: DY

      real(ReKi)                 :: Y_Grid

      ! Calculate position on Y grid
      Y_Grid = (PosY + G3D%YHWid)*G3D%InvDY + 1

      ! Calculate bounding grid indices
      IY_LO = floor(Y_Grid, IntKi)
      IY_HI = IY_LO + 1

      ! Position location within interval [0,1]
      DY = Y_Grid - aint(Y_Grid)

      if (IY_LO >= 1 .and. IY_HI <= G3D%NYGrids) then
         DY = 2.0_ReKi*DY - 1.0_ReKi
      else if (IY_LO == 0 .and. DY >= 1.0_ReKi - GridTol) then
         IY_LO = 1
         IY_HI = 2
         DY = -1.0_ReKi
      else if (IY_LO == G3D%NYGrids .and. DY <= GridTol) then
         IY_LO = G3D%NYGrids - 1
         IY_HI = G3D%NYGrids
         DY = 1.0_ReKi
      else if (AllowExtrap) then
         if (IY_LO <= 0) then
            ! Clamp value at grid width below the low side of grid
            DY = 2.0_ReKi*max(PosY/G3D%YHWid + 2.0_ReKi, 0.0_Reki) - 1.0_ReKi
            IY_LO = 1
            IY_HI = 1
            AllExtrap = ior(AllExtrap, ExtrapYmin)
         else if (IY_LO >= G3D%NYGrids) then
            ! Clamp value at grid width above the high side of grid
            DY = 2.0_ReKi*min(PosY/G3D%YHWid - 1.0_ReKi, 1.0_Reki) - 1.0_ReKi
            IY_LO = G3D%NYGrids
            IY_HI = G3D%NYGrids
            AllExtrap = ior(AllExtrap, ExtrapYmax)
         end if
      else
         ! Position outside
         call SetErrStat(ErrID_Fatal, ' GF wind array boundaries violated: Grid too small in Y direction. Y='// &
                         TRIM(Num2LStr(PosY))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*G3D%YHWid))// &
                         ', '//TRIM(Num2LStr(G3D%YHWid))//']', &
                         ErrStat, ErrMsg, RoutineName)
      end if

   end subroutine

   !> GetBoundsZ populates IZ_Lo, IZ_Hi, and the interpolant [-1,1]. It also
   !! adds ExtrapZmin or ExtrapZmax to AllExtrap if applicable.
   subroutine GetBoundsZ(PosZ, DZ)

      real(ReKi), intent(in)     :: PosZ
      real(ReKi), intent(out)    :: DZ

      real(ReKi)                 :: Z_GRID

      ! Calculate position on Z grid
      Z_GRID = (PosZ - G3D%GridBase)*G3D%InvDZ + 1

      ! Calculate bounding grid indices
      IZ_LO = floor(Z_GRID, IntKi)
      IZ_HI = IZ_LO + 1

      ! Position location within interval [-1,1]
      DZ = Z_GRID - aint(Z_GRID)

      ! If indices are within grid, set in grid to true
      if (IZ_LO >= 1 .and. IZ_HI <= G3D%NZGrids) then
         InGrid = .true.
         DZ = 2.0_ReKi*DZ - 1.0_ReKi
         return
      end if

      ! If below grid
      if (IZ_LO < 1) then
         if (IZ_LO == 0 .and. DZ >= 1.0_ReKi - GridTol) then
            InGrid = .true.
            IZ_LO = 1
            IZ_HI = 2
            DZ = -1.0_ReKi
         else if (G3D%InterpTower) then
            ! Interp from bottom of grid to ground (zero velocity)
            InGrid = .false.
            IZ_LO = 0
            IZ_HI = 1
            DZ = 2.0_ReKi*(PosZ/G3D%GridBase) - 1.0_ReKi
         else if (G3D%NTGrids > 0) then
            ! Interpolate with tower grid
            InGrid = .false.
            ! Tower grid is reversed (lowest index is top of tower)
            IZ_LO = int(-(Z_GRID - 1)) + 1
            if (IZ_LO >= G3D%NTGrids) then
               ! Between end of tower grid and ground (zero velocity)
               IZ_LO = G3D%NTGrids
               DZ = 1.0_ReKi - 2.0_ReKi*(PosZ/(G3D%GridBase - real(IZ_LO - 1, ReKi)/G3D%InvDZ))
            else
               ! Within tower grid
               DZ = 2.0_ReKi*(real(2 - IZ_LO, ReKi) - Z_GRID) - 1.0_ReKi
            end if
            IZ_HI = IZ_LO + 1
         else if (AllowExtrap) then
            InGrid = .true.
            IZ_LO = 1
            IZ_HI = 1
            DZ = 2.0_ReKi*max(PosZ/G3D%GridBase, 0.0_Reki) - 1.0_ReKi
            AllExtrap = ior(AllExtrap, ExtrapZmin)
         else
            ! Position below grid
            call SetErrStat(ErrID_Fatal, ' G3D wind array boundaries violated. '// &
                            'Grid too small in Z direction '// &
                            '(height (Z='//TRIM(Num2LStr(Position(3)))// &
                            ' m) is below the grid and no tower points are defined).', &
                            ErrStat, ErrMsg, RoutineName)
         end if
         return
      end if

      ! If above grid
      if (IZ_HI > G3D%NZGrids) then
         if (IZ_HI == G3D%NZGrids + 1 .and. DZ <= GridTol) then
            InGrid = .true.
            IZ_LO = G3D%NZGrids - 1
            IZ_HI = G3D%NZGrids
            DZ = 1.0_ReKi
         else if (AllowExtrap) then
            InGrid = .true.
            IZ_LO = G3D%NZGrids
            IZ_HI = G3D%NZGrids
            ! Calculate interpolation, limit to value at grid width above top of grid
            DZ = 2.0_ReKi*min((Position(3) - (G3D%GridBase + 2*G3D%ZHWid))/G3D%ZHWid, 1.0_ReKi) - 1.0_ReKi
            AllExtrap = ior(AllExtrap, ExtrapZmax)
         else
            ! Position above grid
            call SetErrStat(ErrID_Fatal, ' G3D wind array boundaries violated. '// &
                            'Grid too small in Z direction '// &
                            '(Z='//TRIM(Num2LStr(Position(3)))//' m is above grid.)', &
                            ErrStat, ErrMsg, RoutineName)
         end if
         return
      end if

   end subroutine

   !> GetBoundsT populates IT_Lo, IT_Hi, and the interpolant [-1,1].
   subroutine GetBoundsT(PosX, DT)

      real(ReKi), intent(in)     :: PosX
      real(ReKi), intent(out)    :: DT

      real(ReKi)                 :: TimeShifted
      real(ReKi)                 :: T_GRID

      ! Perform the time shift. At time=0, a point half the grid width downstream
      ! (p%YHWid) will index into the zero time slice. If we did not do this,
      ! any point downstream of the tower at the beginning of the run would
      ! index outside of the array. This all assumes the grid width is at least as
      ! large as the rotor. If it isn't, then the interpolation will not work.

      ! In distance, X: InputInfo%PosX - p%InitXPosition - TIME*p%MeanWS
      TimeShifted = real(Time, ReKi) + (G3D%InitXPosition - PosX)*G3D%InvMWS

      ! If field is periodic
      if (G3D%Periodic) then
         TimeShifted = MODULO(TimeShifted, G3D%TotalTime)
         ! If TimeShifted is a very small negative number,
         ! modulo returns the incorrect value due to internal rounding errors.
         ! See bug report #471
         if (TimeShifted == G3D%TotalTime) TimeShifted = 0.0_ReKi
      end if

      ! Get position on T grid
      T_GRID = TimeShifted*G3D%Rate + 1

      ! Calculate bounding grid indices
      IT_LO = floor(T_GRID, IntKi)
      IT_HI = ceiling(T_GRID, IntKi)

      ! Position location within interval [0,1]
      DT = T_GRID - aint(T_GRID)

      ! Adjust indices and interpolant
      if (IT_LO >= 1 .and. IT_HI <= G3D%NSteps) then
         ! Point is within grid
         DT = 2.0_ReKi*DT - 1.0_ReKi
      else if (IT_LO == G3D%NSteps) then
         if (G3D%Periodic) then
            ! Time wraps back to beginning
            IT_HI = 1
            DT = 2.0_ReKi*DT - 1.0_ReKi
         else if (DT <= GridTol) then
            ! Within tolerance of last time
            IT_HI = IT_LO
            DT = -1.0_Reki
         else
            ! Extrapolate
            IT_LO = G3D%NSteps - 1
            IT_HI = G3D%NSteps
            DT = DT + 1.0_ReKi
         end if
      else
         ! Time exceeds array bounds
         call SetErrStat(ErrID_Fatal, ' Error: GF wind array was exhausted at '// &
                         TRIM(Num2LStr(TIME))//' seconds (trying to access data at '// &
                         TRIM(Num2LStr(TimeShifted))//' seconds).', &
                         ErrStat, ErrMsg, RoutineName)
      end if

   end subroutine

end subroutine

pure function Grid3DField_GetVelLinear(VelCell, Xi, Is3D) result(Velocity)

   real(ReKi), intent(in)    :: VelCell(8, 3)  !< velocities at corners of grid cell
   real(ReKi), intent(in)    :: Xi(3)          !< isoparametric coordinates in cell (dy, dz, dt) [-1, +1]
   logical, intent(in)       :: Is3D           !< flag for 3D or 2D grid
   real(ReKi)                :: Velocity(3)    !< The U, V, W velocities

   real(ReKi)                :: N(8)           ! Shape function values
   integer(IntKi)            :: IC

   if (Is3D) then

      ! Get 3D interpolation weights
      N(1) = (1.0_ReKi - Xi(1))*(1.0_ReKi - Xi(2))
      N(2) = (1.0_ReKi + Xi(1))*(1.0_ReKi - Xi(2))
      N(3) = (1.0_ReKi - Xi(1))*(1.0_ReKi + Xi(2))
      N(4) = (1.0_ReKi + Xi(1))*(1.0_ReKi + Xi(2))
      N(5) = (1.0_ReKi - Xi(1))*(1.0_ReKi - Xi(2))
      N(6) = (1.0_ReKi + Xi(1))*(1.0_ReKi - Xi(2))
      N(7) = (1.0_ReKi - Xi(1))*(1.0_ReKi + Xi(2))
      N(8) = (1.0_ReKi + Xi(1))*(1.0_ReKi + Xi(2))
      N(1:4) = N(1:4)*(1.0_ReKi - Xi(3))/8.0_ReKi
      N(5:8) = N(5:8)*(1.0_ReKi + Xi(3))/8.0_ReKi

      ! Calculate velocity
      do ic = 1, 3
         Velocity(ic) = dot_product(VelCell(:, ic), N)
      end do

   else

      ! Get 2D interpolation weights
      N(1) = (1.0_ReKi - Xi(2))*(1.0_ReKi - Xi(3))/4.0_ReKi
      N(2) = (1.0_ReKi + Xi(2))*(1.0_ReKi - Xi(3))/4.0_ReKi
      N(3) = (1.0_ReKi - Xi(2))*(1.0_ReKi + Xi(3))/4.0_ReKi
      N(4) = (1.0_ReKi + Xi(2))*(1.0_ReKi + Xi(3))/4.0_ReKi

      ! Calculate velocity
      do ic = 1, 3
         Velocity(ic) = dot_product(VelCell(1:4, ic), N(1:4))
      end do

   end if

end function

subroutine Grid3DField_GetVelAccCubic(DTime, VelCell, AccCell, Xi, Is3D, Velocity, Accel)

   real(ReKi), intent(in)              :: DTime          !< cell delta time
   real(ReKi), intent(in)              :: VelCell(8, 3)  !<
   real(ReKi), intent(in)              :: AccCell(8, 3)  !<
   real(ReKi), intent(in)              :: Xi(3)          !< cell distance (dy, dz, dt) [-1, +1]
   logical, intent(in)                 :: Is3D           !<
   real(ReKi), intent(out), optional   :: Velocity(3)    !<
   real(ReKi), intent(out), optional   :: Accel(3)       !<

   character(*), parameter    :: RoutineName = "Grid3DField_GetVelAccCubic"
   integer(IntKi)             :: IC
   real(ReKi)                 :: N(4)
   real(ReKi)                 :: P(3, 2), PP(3, 2)
   real(ReKi)                 :: t, C1, C2, C3, C4

   ! If 3D interpolation
   if (Is3D) then

      ! Get interpolation weights
      N(1) = (1.0_ReKi - Xi(1))*(1.0_ReKi - Xi(2))/4.0_ReKi
      N(2) = (1.0_ReKi + Xi(1))*(1.0_ReKi - Xi(2))/4.0_ReKi
      N(3) = (1.0_ReKi - Xi(1))*(1.0_ReKi + Xi(2))/4.0_ReKi
      N(4) = (1.0_ReKi + Xi(1))*(1.0_ReKi + Xi(2))/4.0_ReKi

      ! Calculate velocity and acceleration at lo and hi time
      do IC = 1, 3
         P(IC, 1) = dot_product(VelCell(1:4, IC), N)  ! lo time
         P(IC, 2) = dot_product(VelCell(5:8, IC), N)  ! hi time
         PP(IC, 1) = dot_product(AccCell(1:4, IC), N) ! lo time
         PP(IC, 2) = dot_product(AccCell(5:8, IC), N) ! hi time
      end do

   else  ! 2D (Tower)

      ! Get interpolation weights
      N(1) = (1.0_ReKi - Xi(2))/2.0_ReKi
      N(2) = (1.0_ReKi + Xi(2))/2.0_ReKi

      ! Calculate velocity and acceleration at lo and hi time
      do IC = 1, 3
         P(IC, 1) = dot_product(VelCell(1:2, IC), N(1:2))   ! lo time
         P(IC, 2) = dot_product(VelCell(3:4, IC), N(1:2))   ! hi time
         PP(IC, 1) = dot_product(AccCell(1:2, IC), N(1:2))  ! lo time
         PP(IC, 2) = dot_product(AccCell(3:4, IC), N(1:2))  ! hi time
      end do

   end if

   ! Calculate interval percent
   t = (Xi(3) + 1)/2.0_ReKi

   ! If velocity requested
   if (present(Velocity)) then
      C1 = 2.0_ReKi*t*t*t - 3.0_ReKi*t*t + 1.0_ReKi
      C2 = (t*t*t - 2.0_ReKi*t*t + t)*DTime
      C3 = -2.0_ReKi*t*t*t + 3.0_ReKi*t*t
      C4 = (t*t*t - t*t)*DTime
      Velocity = C1*P(:, 1) + C2*PP(:, 1) + C3*P(:, 2) + C4*PP(:, 2)
   end if

   ! If acceleration requested
   if (present(Accel)) then
      C1 = (6.0_ReKi*t*t - 6.0_ReKi*t)/DTime
      C2 = 3.0_ReKi*t*t - 4.0_ReKi*t + 1.0_ReKi
      C3 = -C1
      C4 = 3.0_ReKi*t*t - 2.0_ReKi*t
      Accel = C1*P(:, 1) + C2*PP(:, 1) + C3*P(:, 2) + C4*PP(:, 2)
   end if

end subroutine

subroutine IfW_Grid3DField_CalcAccel(G3D, ErrStat, ErrMsg)
   type(Grid3DFieldType), intent(inout)   :: G3D
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = "Grid3DField_CalcAccel"
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg
   integer(IntKi)                         :: ic, iy, iz
   real(ReKi), allocatable                :: u(:), dy2(:)

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Allocate storage for acceleration grid
   call AllocAry(G3D%Acc, size(G3D%Vel, dim=1), size(G3D%Vel, dim=2), &
                 size(G3D%Vel, dim=3), size(G3D%Vel, dim=4), &
                 'grid-field velocity data', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! If number of time grids is 1 or 2, set all accelerations to zero, return
   if (G3D%NTGrids < 3) then
      G3D%Acc = 0.0_SiKi
      return
   end if

   ! Allocate storage for U used in cubic spline derivative calc
   call AllocAry(U, G3D%NSteps, "storage for U", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for V used in cubic spline derivative calc
   call AllocAry(dy2, G3D%NSteps, "storage for V", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Loop through grid points and calculate derivative using spline
   do iz = 1, G3D%NZGrids
      do iy = 1, G3D%NYGrids
         do ic = 1, G3D%NComp
            call CalcCubicSplineDeriv(G3D%NSteps, G3D%DTime, G3D%Vel(ic, iy, iz, :), G3D%Acc(ic, iy, iz, :))
         end do
      end do
   end do

   ! If grid field does not include tower grids, return
   if (G3D%NTGrids == 0) return

   ! Allocate storage for tower acceleration
   call AllocAry(G3D%AccTower, size(G3D%VelTower, dim=1), &
                 size(G3D%VelTower, dim=2), size(G3D%VelTower, dim=3), &
                 'tower wind acceleration data.', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! If number of time grids is 1 or 2, set all accelerations to zero
   if (G3D%NTGrids < 3) then
      G3D%Acc = 0.0_SiKi
   else  ! Otherwise, calculate acceleration at each grid point
      do iz = 1, G3D%NTGrids
         do ic = 1, G3D%NComp
            call CalcCubicSplineDeriv(G3D%NSteps, G3D%DTime, G3D%VelTower(ic, iz, :), G3D%AccTower(ic, iz, :))
         end do
      end do
   end if

contains

   !> CalcCubicSplineDeriv fits a cubic spline through the y array with points
   !! spaced a constant 'h' apart. It then calculates the corresponding
   !! derivative of y with respect to x at the same x values and returns it
   !! in the dy array.
   subroutine CalcCubicSplineDeriv(n, h, y, dy)
      integer(IntKi), intent(in)    :: n     ! number of points
      real(ReKi), intent(in)        :: h     ! delta time
      real(SiKi), intent(in)        :: y(:)  ! value at each time
      real(SiKi), intent(out)       :: dy(:) ! value derivative at each time

      integer(IntKi)                :: i
      real(ReKi)                    :: p, un

      ! If periodic function, set beginning and end to have same slope
      if (G3D%Periodic) then
         dy(1) = real((y(2) - y(n - 1))/(2.0_ReKi*h), SiKi)
         dy(n) = dy(1)
      else
         dy(1) = 0.0_ReKi
         dy(n) = 0.0_ReKi
      end if

      ! Apply first derivative at lower boundary condition
      dy2(1) = -0.5_ReKi
      u(1) = 3.0_ReKi*((y(2) - y(1))/h - dy(1))/h

      ! Decomposition
      do i = 2, n - 1
         p = 0.5_ReKi*dy2(i - 1) + 2.0_ReKi
         dy2(i) = -0.5_ReKi/p
         u(i) = (6.*((y(i + 1) - 2.0_ReKi*y(i) + y(i - 1))/h)/(2.0_ReKi*h) - 0.5_ReKi*u(i - 1))/p
      end do

      ! Apply first derviative at upper boundary condition
      un = 3.0_ReKi*(dy(n) - (y(n) - y(n - 1))/h)/h
      dy2(n) = (un - 0.5_ReKi*u(n - 1))/(0.5_ReKi*dy2(n - 1) + 1.0_ReKi)

      ! Back substitution and derivative calculation
      do i = n - 1, 1, -1
         dy2(i) = dy2(i)*dy2(i + 1) + u(i)
         dy(i) = real((y(i + 1) - y(i))/h - h*(dy2(i)/3.0_ReKi + dy2(i + 1)/6.0_ReKi), SiKi)
      end do

   end subroutine

end subroutine

!> This subroutine generates the mean wind vector timeseries for each height above the ground.  This
!! is essentially compressing the Y dimension of the wind box leaving a Z-T plane of vectors.  The
!! resulting dimensions will be (NZGrids, NYGrids, NFFComp, NFFSteps)
subroutine IfW_Grid3DField_CalcVelAvgProfile(G3D, CalcAccel, ErrStat, ErrMsg)

   type(Grid3DFieldType), intent(inout)   :: G3D         !< Parameters
   logical, intent(inout)                 :: CalcAccel   !< Flag to calculate acceleration
   integer(IntKi), intent(out)            :: ErrStat     !< Error status of the operation
   character(*), intent(out)              :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter    :: RoutineName = 'IfW_Grid3DField_CalcVelAvgProfile'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Allocate velocity array
   if (.not. allocated(G3D%VelAvg)) then
      call AllocAry(G3D%VelAvg, G3D%NComp, G3D%NZGrids, G3D%NSteps, &
                    'Full-field average wind velocity timeseries data array.', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
      G3D%VelAvg = 0.0_SiKi
   end if

   ! Calculate average velocity for each component across grid (Y)
   G3D%VelAvg = sum(G3D%Vel, dim=2)/G3D%NYGrids

   ! If acceleration calculation not requested, return
   if (.not. CalcAccel) return

   ! Allocate acceleration array
   if (.not. allocated(G3D%AccAvg)) then
      call AllocAry(G3D%AccAvg, G3D%NComp, G3D%NZGrids, G3D%NSteps, &
                    'Full-field average wind acceleration timeseries data array.', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
      G3D%AccAvg = 0.0_SiKi
   end if

   ! Calculate average acceleration for each component across grid (Y)
   G3D%AccAvg = sum(G3D%Acc, dim=2)/G3D%NYGrids

end subroutine

subroutine Grid4DField_GetVel(G4D, Time, Position, Velocity, ErrStat, ErrMsg)

   type(Grid4DFieldType), intent(in)   :: G4D            !< 4D grid-field data
   real(DbKi), intent(in)              :: Time           !< time to get value
   real(ReKi), intent(in)              :: Position(3)    !< position X,Y,Z to get value
   real(ReKi), intent(out)             :: Velocity(3)    !< The U, V, W velocities
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg

   character(*), parameter             :: RoutineName = "Grid4DField_GetVel"

   integer(IntKi)                      :: Indx_Lo(4)  ! index associated with lower bound of dimension 1-4 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   integer(IntKi)                      :: Indx_Hi(4)  ! index associated with upper bound of dimension 1-4 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   real(ReKi)                          :: xi(4)       ! isoparametric coordinates
   real(ReKi)                          :: N(16, 1)    ! Shape function
   real(ReKi)                          :: P(3, 16)    ! Point values
   real(ReKi)                          :: tmp
   integer(IntKi)                      :: i

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Find the bounding indices for XYZ position
   !----------------------------------------------------------------------------

   do i = 1, 3
      tmp = (Position(i) - G4D%pZero(i))/G4D%delta(i)
      Indx_Lo(i) = INT(tmp) + 1                          ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
      xi(i) = 2.0_ReKi*(tmp - aint(tmp)) - 1.0_ReKi   ! convert to value between -1 and 1
   end do

   !----------------------------------------------------------------------------
   ! Find the bounding indices for time
   !----------------------------------------------------------------------------

   i = 4
   tmp = real((Time - G4D%TimeStart)/G4D%delta(i), ReKi)
   Indx_Lo(i) = INT(tmp) + 1     ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
   xi(i) = 2.0_ReKi*(tmp - aint(tmp)) - 1.0_ReKi  ! convert to value between -1 and 1
   if ((Indx_Lo(i) == G4D%n(i))) then
      if (abs(xi(i) + 1.0_SiKi) < 0.001_SiKi) then    ! Allow for the special case where Time = TgridStart + deltat*( n_high_low - 1 )
         Indx_Lo(i) = Indx_Lo(i) - 1
         xi(i) = 1.0_SiKi
      end if
   end if

   !----------------------------------------------------------------------------
   ! Return error if outside bounds
   !----------------------------------------------------------------------------

   do i = 1, 4
      if (Indx_Lo(i) <= 0) then
         Indx_Lo(i) = 1
         call SetErrStat(ErrID_Fatal, 'Outside the grid bounds.', ErrStat, ErrMsg, RoutineName)
         return
      elseif (Indx_Lo(i) >= G4D%n(i)) then
         Indx_Lo(i) = max(G4D%n(i) - 1, 1)           ! make sure it's a valid index
         call SetErrStat(ErrID_Fatal, 'Outside the grid bounds.', ErrStat, ErrMsg, RoutineName)
         return
      end if
      Indx_Hi(i) = min(Indx_Lo(i) + 1, G4D%n(i))     ! make sure it's a valid index
   end do

   !----------------------------------------------------------------------------
   ! Clamp isopc to [-1, 1] so we don't extrapolate (effectively nearest neighbor)
   !----------------------------------------------------------------------------

   xi = min(+1.0_ReKi, max(-1.0_ReKi, xi))

   !----------------------------------------------------------------------------
   ! compute weighting factors
   !----------------------------------------------------------------------------

   N(1, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi - xi(4))
   N(2, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi - xi(4))
   N(3, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi - xi(4))
   N(4, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi - xi(4))
   N(5, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi - xi(4))
   N(6, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi - xi(4))
   N(7, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi - xi(4))
   N(8, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi - xi(4))
   N(9, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi + xi(4))
   N(10, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi + xi(4))
   N(11, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi + xi(4))
   N(12, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi - xi(3))*(1.0_ReKi + xi(4))
   N(13, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi + xi(4))
   N(14, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi - xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi + xi(4))
   N(15, 1) = (1.0_ReKi - xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi + xi(4))
   N(16, 1) = (1.0_ReKi + xi(1))*(1.0_ReKi + xi(2))*(1.0_ReKi + xi(3))*(1.0_ReKi + xi(4))
   N = N/16.0_ReKi

   !----------------------------------------------------------------------------
   ! Get point values
   !----------------------------------------------------------------------------

   P(:, 1) = G4D%Vel(:, Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), Indx_Lo(4))
   P(:, 2) = G4D%Vel(:, Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), Indx_Lo(4))
   P(:, 3) = G4D%Vel(:, Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), Indx_Lo(4))
   P(:, 4) = G4D%Vel(:, Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), Indx_Lo(4))
   P(:, 5) = G4D%Vel(:, Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), Indx_Lo(4))
   P(:, 6) = G4D%Vel(:, Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), Indx_Lo(4))
   P(:, 7) = G4D%Vel(:, Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), Indx_Lo(4))
   P(:, 8) = G4D%Vel(:, Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), Indx_Lo(4))
   P(:, 9) = G4D%Vel(:, Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), Indx_Hi(4))
   P(:, 10) = G4D%Vel(:, Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), Indx_Hi(4))
   P(:, 11) = G4D%Vel(:, Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), Indx_Hi(4))
   P(:, 12) = G4D%Vel(:, Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), Indx_Hi(4))
   P(:, 13) = G4D%Vel(:, Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), Indx_Hi(4))
   P(:, 14) = G4D%Vel(:, Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), Indx_Hi(4))
   P(:, 15) = G4D%Vel(:, Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), Indx_Hi(4))
   P(:, 16) = G4D%Vel(:, Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), Indx_Hi(4))

   !----------------------------------------------------------------------------
   ! Interpolate
   !----------------------------------------------------------------------------

   Velocity = pack(matmul(P, N), .true.)

end subroutine

subroutine UserField_GetVel(UF, Time, Position, Velocity, ErrStat, ErrMsg)

   type(UserFieldType), intent(in)     :: UF             !< user-field data
   real(DbKi), intent(in)              :: Time           !< time to get value
   real(ReKi), intent(in)              :: Position(3)    !< position X,Y,Z to get value
   real(ReKi), intent(out)             :: Velocity(3)    !< The U, V, W velocities
   integer(IntKi), intent(out)         :: ErrStat
   character(*), intent(out)           :: ErrMsg

   character(*), parameter             :: RoutineName = "UserField_GetVel"

   ErrStat = ErrID_None
   ErrMsg = ""

   Velocity = 0.0_ReKi
   call SetErrStat(ErrID_Fatal, "UserField_GetVel not implemented", ErrStat, ErrMsg, RoutineName)

end subroutine

function CalculateMeanVelocity(G3D, z, y) result(u)

   type(Grid3DFieldType), intent(IN)   :: G3D   !< Parameters
   real(ReKi), intent(IN)              :: Z     ! height
   real(ReKi), intent(IN)              :: y     ! lateral location
   real(ReKi)                          :: u     ! mean wind speed at position (y,z)

   if  (Z <= 0.0_ReKi) then
      U = 0.0_ReKi
      return
   end if
   
   select case (G3D%WindProfileType)

   case (WindProfileType_PL)

      U = G3D%MeanWS*(Z/G3D%RefHeight)**G3D%PLExp      ! [IEC 61400-1 6.3.1.2 (10)]

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

   if (G3D%VLinShr .ne. 0.0_ReKi) then  ! Add vertical linear shear, if has
      U = U + G3D%MeanWS*G3D%VLinShr*(Z - G3D%RefHeight)/G3D%RefLength
   end if

   if (G3D%HLinShr .ne. 0.0_ReKi) then  ! Add horizontal linear shear, if has
      U = U + G3D%MeanWS*G3D%HLinShr*y/G3D%RefLength
   end if

end function CalculateMeanVelocity

subroutine Uniform_to_Grid3D(UF, InterpCubic, G3D, ErrStat, ErrMsg)

   type(UniformFieldType), intent(IN)  :: UF                !< UniformWind Parameters
   logical, intent(in)                 :: InterpCubic       !< Flag to use cubic interpolation
   type(Grid3DFieldType), intent(OUT)  :: G3D               !< FF Parameters
   integer(IntKi), intent(OUT)         :: ErrStat           !< error status
   character(*), intent(OUT)           :: ErrMsg            !< error message

   character(*), parameter             :: RoutineName = 'Uniform_to_FFWind'
   integer(ReKi), parameter            :: dz = 5.0
   integer(ReKi), parameter            :: dy = 5.0
   real(DbKi)                          :: Time
   real(ReKi)                          :: PositionXYZ(3)
   type(UniformField_Interp)           :: op
   integer(IntKi)                      :: n, i, it, iy, iz
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   G3D%WindFileFormat = -1      ! Binary file format description number
   G3D%NComp = 3                ! Number of wind components
   G3D%Periodic = .false.
   G3D%InterpTower = .true.
   G3D%RefHeight = UF%RefHeight
   G3D%NTGrids = 0
   G3D%InvDY = 1.0_ReKi/dy       ! reciprocal of delta y (1/meters)
   G3D%InvDZ = 1.0_ReKi/dz       ! reciprocal of delta z (1/meters)

   ! Number of points in the lateral (y) direction of the grids
   n = nint(UF%RefLength*1.1_ReKi*0.5_ReKi/dy)
   G3D%NYGrids = n*2 + 1

   ! Number of points in the vertical (z) direction of the grids
   n = nint(UF%RefLength*1.1_ReKi*0.5_ReKi/dz)
   G3D%NZGrids = nint(G3D%RefHeight/dy) + n + 1

   ! Half the grid width (meters)
   G3D%YHWid = 0.5_ReKi*dy*(G3D%NYGrids - 1)

   ! Half the grid height (meters)
   G3D%ZHWid = 0.5_ReKi*dz*(G3D%NZGrids - 1)

   ! Height of the bottom of the grid (meters)
   G3D%GridBase = G3D%RefHeight + n*dz - G3D%ZHWid*2.0_ReKi

   ! Initial x position of grid (distance in FF is offset) meters)
   G3D%InitXPosition = 0.0_ReKi

   ! time will be the smallest delta t in this Uniform wind file
   if (UF%DataSize < 2) then
      G3D%DTime = 600.0_ReKi     ! doesn't matter what the time step is
      G3D%NSteps = 2             ! "Number of time steps in the FF array
   else
      G3D%DTime = minval(UF%Time(2:) - UF%Time(:size(UF%Time) - 1))   ! Delta time (seconds)
      if (G3D%DTime < 0.0001) then
         call SetErrStat(ErrID_Fatal, "Smallest time step in uniform wind file is less that 0.0001 seconds. "// &
                         "Increase the time step to convert to a FF file.", ErrStat, ErrMsg, RoutineName)
         return
      end if
      G3D%NSteps = NINT(UF%Time(UF%DataSize)/G3D%DTime) + 1
   end if

   G3D%Rate = 1.0_ReKi/G3D%DTime                ! Data rate (1/DTime)
   G3D%AddMeanAfterInterp = .false.             ! Add the mean wind speed after interpolating at a given height?
   G3D%WindProfileType = WindProfileType_PL     ! Wind profile type (0=constant;1=logarithmic;2=power law)
   G3D%PLExp = GetAverageVal(UF%ShrV)           ! Power law exponent (used for PL wind profile type only)
   G3D%Z0 = 0.0_ReKi                            ! Surface roughness length (used for LOG wind profile type only)
   G3D%TotalTime = (G3D%NSteps - 1)*G3D%DTime   ! The total time of the simulation (seconds)

   ! Allocate velocity array
   call AllocAry(G3D%Vel, G3D%NComp, G3D%NYGrids, G3D%NZGrids, G3D%NSteps, 'G3D%Vel', ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Initialize position
   PositionXYZ = 0.0_ReKi

   ! Loop through time steps
   do it = 1, G3D%NSteps

      ! Calculate time
      Time = (it - 1)*G3D%DTime

      ! Get operating point
      if (InterpCubic) then
         op = UniformField_InterpCubic(UF, Time)
      else
         op = UniformField_InterpLinear(UF, Time)
      end if

      ! Loop through y grid
      do iy = 1, G3D%NYGrids

         ! Calculate Y position
         PositionXYZ(2) = (iy - 1)*dy - G3D%YHWid

         ! Loop through z grid
         do iz = 1, G3D%NZGrids

            ! Calculate Z position
            PositionXYZ(3) = (iz - 1)*dz + G3D%GridBase

            ! If Z is zero or less
            if (PositionXYZ(3) <= 0.0_Reki) then
               ! Set wind velocity to zero
               G3D%Vel(:, iy, iz, it) = 0.0_SiKi
            else
               ! Calculate velocity at operating point and position, store in grid
               G3D%Vel(:, iy, iz, it) = real(UniformField_GetVel(UF, op, PositionXYZ), SiKi)
            end if
         end do ! iz
      end do ! iy
   end do ! it

   ! compute some averages for this simulation
   G3D%MeanWS = GetAverageVal(UF%VelH)  ! Mean wind speed (advection speed)
   G3D%InvMWS = 1.0_ReKi/G3D%MeanWS

contains

   function GetAverageVal(Ary) result(Avg)
      real(ReKi), intent(in)  :: Ary(:)
      real(ReKi)              :: Avg

      ! If array has one element, average is first value
      if (UF%DataSize < 2) then
         Avg = Ary(1)
         return
      end if

      Avg = UF%Time(1)*Ary(1) ! in case tData(1) /= 0
      do i = 2, UF%DataSize
         Avg = Avg + (UF%Time(i) - UF%Time(i - 1))*(Ary(i) + Ary(i - 1))
      end do
      Avg = Avg/(UF%Time(UF%DataSize) - UF%Time(1))/2.0_ReKi

   end function GetAverageVal

end subroutine Uniform_to_Grid3D

subroutine Grid3D_to_Uniform(G3D, UF, ErrStat, ErrMsg, SmoothingRadius)

   type(Grid3DFieldType), intent(IN)   :: G3D                  !< FF Parameters
   type(UniformFieldType), intent(OUT) :: UF                   !< UniformWind Parameters
   integer(IntKi), intent(OUT)         :: ErrStat              !< error status
   character(*), intent(OUT)           :: ErrMsg               !< error message
   real(ReKi), optional, intent(IN)    :: SmoothingRadius      !< length of time used for smoothing data, seconds (if omitted, no smoothing will occur)

   character(*), parameter             :: RoutineName = 'FFWind_to_Uniform'
   integer(IntKi)                      :: i
   integer(IntKi)                      :: iy_ref, iz_ref, iz_p1
   integer(IntKi)                      :: iy, iz, ic
   real(ReKi)                          :: meanVel(3)
   real(ReKi)                          :: meanWindDir
   real(ReKi)                          :: u_p1, z_p1
   real(ReKi), parameter               :: HubPositionX = 0.0_ReKi
   real(ReKi)                          :: radius ! length of time to use for smoothing uniform wind data, seconds
   integer(IntKi)                      :: ErrStat2
   character(ErrMsgLen)                :: ErrMsg2

   real(SiKi), allocatable             :: Vel(:, :, :, :)
   real(ReKi), allocatable             :: tmp(:)
   real(R8Ki)                          :: transformMat(3, 3)

   ErrStat = ErrID_None
   ErrMsg = ""

   if (G3D%RefLength > epsilon(0.0_ReKi)) then
      UF%RefLength = G3D%VLinShr
   else
      UF%RefLength = (G3D%nYGrids - 1)/G3D%InvDY        ! width of the FF wind field
   end if

   UF%DataSize = G3D%NSteps

   call AllocAry(UF%Time, UF%DataSize, 'time', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%VelH, UF%DataSize, 'horizontal wind speed', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%AngleH, UF%DataSize, 'direction', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%AngleV, UF%DataSize, 'upflow', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%VelV, UF%DataSize, 'vertical wind speed', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%ShrH, UF%DataSize, 'horizontal linear shear', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%ShrV, UF%DataSize, 'vertical power-law shear exponent', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%LinShrV, UF%DataSize, 'vertical linear shear', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(UF%VelGust, UF%DataSize, 'gust velocity', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call AllocAry(Vel, G3D%NZGrids, G3D%NYGrids, G3D%NComp, G3D%NSteps, 'FFData', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(tmp, G3D%NSteps, 'tmp', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   if (ErrStat >= AbortErrLev .or. UF%DataSize < 1) then
      if (allocated(Vel)) deallocate (Vel)
      if (allocated(tmp)) deallocate (tmp)
      return
   end if

   ! we'll assume these are 0, for simplicity
   UF%ShrH = G3D%HLinShr
   UF%LinShrV = G3D%VLinShr
   UF%VelGust = 0.0_ReKi

   ! fill time array with time at hub position
   do i = 1, UF%DataSize
      UF%Time(i) = (G3D%InitXPosition - HubPositionX)*G3D%InvMWS + (i - 1)*G3D%DTime
   end do

   ! calculate mean velocity at grid point nearest lateral center of grid at reference height:
   iy_ref = nint(G3D%nYGrids/2.0_ReKi)
   iz_ref = nint((G3D%RefHeight - G3D%GridBase)*G3D%InvDZ) + 1
   UF%RefHeight = G3D%GridBase + (iz_ref - 1)/G3D%InvDZ ! make sure RefHt is on the grid

   ! Calculate mean value for each component through
   meanVel = sum(G3D%Vel(:, iy_ref, iz_ref, :), dim=2)/UF%DataSize

   ! calculate the average upflow angle
   UF%AngleV = atan2(meanVel(3), TwoNorm(meanVel(1:2)))
   meanWindDir = atan2(meanVel(2), meanVel(1))

   ! rotate the FF wind to remove the mean upflow and direction
   transformMat(1, 1) = cos(meanWindDir)*cos(UF%AngleV(1))
   transformMat(2, 1) = -sin(meanWindDir)
   transformMat(3, 1) = -cos(meanWindDir)*sin(UF%AngleV(1))

   transformMat(1, 2) = sin(meanWindDir)*cos(UF%AngleV(1))
   transformMat(2, 2) = cos(meanWindDir)
   transformMat(3, 2) = -sin(meanWindDir)*sin(UF%AngleV(1))

   transformMat(1, 3) = sin(UF%AngleV(1))
   transformMat(2, 3) = 0.0_R8Ki
   transformMat(3, 3) = cos(UF%AngleV(1))

   do ic = 1, size(Vel, 4)
      do iy = 1, size(Vel, 2)
         do iz = 1, size(Vel, 1)
            Vel(iz, iy, :, ic) = real(matmul(transformMat, G3D%Vel(:, iy, iz, i)), SiKi)
         end do
      end do
   end do

   ! make sure we have the correct mean, or the direction will also be off here
   if (G3D%AddMeanAfterInterp) then
      Vel(iz_ref, iy_ref, 1, :) = Vel(iz_ref, iy_ref, 1, :) + &
                                  real(CalculateMeanVelocity(G3D, UF%RefHeight, 0.0_ReKi), SiKi)
   end if

   meanVel = 0.0_ReKi
   do i = 1, UF%DataSize
      meanVel = meanVel + Vel(iz_ref, iy_ref, :, i)
   end do
   meanVel = meanVel/UF%DataSize

   ! Fill velocity arrays for uniform wind
   do i = 1, UF%DataSize
      UF%VelH(i) = TwoNorm(Vel(iz_ref, iy_ref, 1:2, i))
   end do
   UF%VelV = Vel(iz_ref, iy_ref, 3, :)

   ! Fill wind direction array
   do i = 1, UF%DataSize
      UF%AngleH(i) = -(meanWindDir + atan2(Vel(iz_ref, iy_ref, 2, i), Vel(iz_ref, iy_ref, 1, i)))
   end do

   ! Now, time average values, if desired:
   if (present(SmoothingRadius)) then
      radius = SmoothingRadius
   else
      radius = 0.0_ReKi
   end if

   tmp = UF%VelH; call kernelSmoothing(UF%Time, tmp, kernelType_TRIWEIGHT, radius, UF%VelH)
   tmp = UF%VelV; call kernelSmoothing(UF%Time, tmp, kernelType_TRIWEIGHT, radius, UF%VelV)
   tmp = UF%AngleH; call kernelSmoothing(UF%Time, tmp, kernelType_TRIWEIGHT, radius, UF%AngleH)

   ! Calculate averaged power law coefficient:
   if (G3D%WindProfileType == WindProfileType_PL) then
      UF%ShrV = G3D%PLExp
   else
      iz_p1 = G3D%nZGrids    ! pick a point to compute the power law exponent (least squares would be better than a single point)
      z_p1 = G3D%GridBase + (iz_p1 - 1)/G3D%InvDZ

      if (G3D%AddMeanAfterInterp) then
         u_p1 = CalculateMeanVelocity(G3D, z_p1, 0.0_ReKi)
      else
         u_p1 = 0.0_ReKi
         do i = 1, UF%DataSize
            u_p1 = u_p1 + Vel(iz_p1, iy_ref, 1, i)
         end do
         u_p1 = u_p1/UF%DataSize
      end if

      if (EqualRealNos(meanVel(1), u_p1) .or. EqualRealNos(u_p1, 0.0_ReKi) .or. EqualRealNos(meanVel(1), 0.0_ReKi)) then
         UF%ShrV = 0.0_ReKi
      else
         UF%ShrV = log(u_p1/meanVel(1))/log(z_p1/UF%RefHeight)
      end if
   end if

end subroutine Grid3D_to_Uniform

end module
