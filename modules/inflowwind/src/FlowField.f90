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

module FlowField

use NWTC_Library
use FlowField_Types

implicit none

public FlowField_GetVelAcc
public UniformField_CalcAccel, Grid3DField_CalcAccel

integer(IntKi), parameter  :: WindProfileType_None = -1     !< don't add wind profile; already included in input
integer(IntKi), parameter  :: WindProfileType_Constant = 0  !< constant wind
integer(IntKi), parameter  :: WindProfileType_Log = 1       !< logarithmic
integer(IntKi), parameter  :: WindProfileType_PL = 2        !< power law

real(ReKi), parameter      :: GridTol = 1.0E-3              ! Tolerance for determining if position is within grid

contains

!> FlowField_GetVelAcc gets the velocities (and accelerations) at the given point positions.
!! Accelerations are only calculated if the AccelUVW array is allocated.
subroutine FlowField_GetVelAcc(FF, IStart, Time, PositionXYZ, VelocityUVW, AccelUVW, ErrStat, ErrMsg)

   type(FlowFieldType), intent(in)           :: FF                !< FlowField data structure
   integer(IntKi), intent(in)                :: IStart            !< Start index for returning velocities for external field
   real(DbKi), intent(in)                    :: Time              !< Time to evaluate velocities/accelerations
   real(ReKi), intent(in)                    :: PositionXYZ(:, :) !< Array of positions to evaluate velocites/accelerations
   real(ReKi), intent(inout)                 :: VelocityUVW(:, :) !< Array of velocity outputs
   real(ReKi), allocatable, intent(inout)    :: AccelUVW(:, :)    !< Array of acceleration outputs
   integer(IntKi), intent(out)               :: ErrStat           !< Error status
   character(*), intent(out)                 :: ErrMsg            !< Error message

   character(*), parameter                   :: RoutineName = "FlowField_GetVelAcc"
   integer(IntKi)                            :: i
   integer(IntKi)                            :: NumPoints
   logical                                   :: CalcAccel
   real(ReKi), allocatable                   :: Position(:, :)
   integer(IntKi)                            :: TmpErrStat
   character(ErrMsgLen)                      :: TmpErrMsg
   ! Uniform Field
   type(UniformField_Interp)                 :: UFopVel, UFopAcc
   ! Grid3D Field
   real(ReKi)                                :: Xi(3)
   real(ReKi)                                :: VelCell(8, 3), AccCell(8, 3)
   logical                                   :: Is3D

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get number of points to evaluate
   NumPoints = size(PositionXYZ, dim=2)

   ! Determine if acceleration should be calculated and returned
   CalcAccel = allocated(AccelUVW)
   if (CalcAccel .and. .not. FF%AccFieldValid) then
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

   ! Copy positions or transform based on wind box rotation
   if (FF%RotateWindBox) then
      do i = 1, NumPoints
         Position(:, i) = GetPrimePosition(PositionXYZ(:, i))
      end do
   else
      Position = PositionXYZ
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

   ! Switch based on flow type
   select case (FF%FieldType)
   case (Uniform_FieldType)

      !-------------------------------------------------------------------------
      ! Uniform Flow Field
      !-------------------------------------------------------------------------

      ! If cubic velocity interp requested, calculate operating point using
      ! cubic interpolation. Otherwise, use linear interpolation
      if (FF%VelInterpCubic) then
         UFopVel = UniformField_IterpCubic(FF%Uniform, Time)
      else
         UFopVel = UniformField_InterpLinear(FF%Uniform, Time)
      end if

      ! If acceleration output is requested
      if (CalcAccel) then

         ! If cubic interpolation was used for the velocity, use same for accel
         ! otherwise, calculate operating point via cubic interpolation
         if (FF%VelInterpCubic) then
            UFopAcc = UFopVel
         else
            UFopAcc = UniformField_IterpCubic(FF%Uniform, Time)
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

      if (CalcAccel) then

         ! Loop through points
         do i = 1, NumPoints

            ! If height < zero, set velocity/acceleration to zero, continue
            if (Position(3, i) <= 0.0_ReKi) then
               VelocityUVW(:, i) = 0.0_ReKi
               AccelUVW(:, i) = 0.0_ReKi
               cycle
            end if

            ! Calculate grid cells for interpolation
            call Grid3DField_GetCell(FF%Grid3D, Time, Position(:, i), CalcAccel, &
                                     VelCell, AccCell, Xi, Is3D, TmpErrStat, TmpErrMsg)
            if (TmpErrStat >= AbortErrLev) then
               call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               return
            end if

            if (FF%VelInterpCubic) then
               ! Cubic velocity and cubic acceleration
               call Grid3DField_GetVelAccCubic(FF%Grid3D%DTime, VelCell, AccCell, Xi, Is3D, &
                                               Velocity=VelocityUVW(:, i), Accel=AccelUVW(:, i))
            else
               ! Linear velocity and cubic acceleration
               VelocityUVW(:, i) = Grid3DField_GetVelLinear(VelCell, Xi, Is3D)
               call Grid3DField_GetVelAccCubic(FF%Grid3D%DTime, VelCell, AccCell, Xi, Is3D, &
                                               Accel=AccelUVW(:, i))
            end if

         end do
      else

         ! Loop through points
         do i = 1, NumPoints

            ! If height < zero, set velocity to zero, continue
            if (Position(3, i) <= 0.0_ReKi) then
               VelocityUVW(:, i) = 0.0_ReKi
               cycle
            end if

            ! Calculate grid cells for interpolation
            call Grid3DField_GetCell(FF%Grid3D, Time, Position(:, i), CalcAccel, &
                                     VelCell, AccCell, Xi, Is3D, TmpErrStat, TmpErrMsg)
            if (TmpErrStat >= AbortErrLev) then
               call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! Cubic or linear velocity interpolation
            if (FF%VelInterpCubic) then
               call Grid3DField_GetVelAccCubic(FF%Grid3D%DTime, VelCell, AccCell, Xi, Is3D, &
                                               Velocity=VelocityUVW(:, i))
            else
               VelocityUVW(:, i) = Grid3DField_GetVelLinear(VelCell, Xi, Is3D)
            end if

         end do
      end if

      ! Add mean wind speed after interpolation if flag is set
      if (FF%Grid3D%AddMeanAfterInterp) then
         do i = 1, NumPoints
            VelocityUVW(1, i) = VelocityUVW(1, i) + GetMeanVelocity(FF%Grid3D, Position(3, i))
         end do
      end if

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
      if (.not. CalcAccel) then
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

   function GetMeanVelocity(GF, PosZ) result(U)
      type(Grid3DFieldType), intent(in)  :: GF
      real(ReKi), intent(in)           :: PosZ
      real(ReKi)                       :: U
      select case (GF%WindProfileType)
      case (WindProfileType_None)
         U = 0.0_ReKi
      case (WindProfileType_PL)
         U = GF%MeanWS*(PosZ/GF%RefHeight)**GF%PLExp ! [IEC 61400-1 6.3.1.2 (10)]
      case (WindProfileType_Log)
         if (.not. EqualRealNos(GF%RefHeight, GF%Z0) .and. PosZ > 0.0_ReKi) then
            U = GF%MeanWS*log(PosZ/GF%Z0)/log(GF%RefHeight/GF%Z0)
         else
            U = 0.0_ReKi
         end if
      case (WindProfileType_Constant)
         U = GF%MeanWS
      case default
         U = 0.0_ReKi
      end select
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

pure function UniformField_IterpCubic(UF, Time) result(op)

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

subroutine UniformField_CalcAccel(UF, ErrStat, ErrMsg)
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

subroutine Grid3DField_GetCell(G3D, Time, Position, CalcAccel, VelCell, AccCell, Xi, Is3D, ErrStat, ErrMsg)

   type(Grid3DFieldType), intent(in)   :: G3D            !< 3D Grid-Field data
   real(DbKi), intent(in)              :: Time           !< time (s)
   real(ReKi), intent(in)              :: Position(3)    !< position X,Y,Z to get value
   logical, intent(in)                 :: CalcAccel      !< flag to populat AccCell
   real(ReKi), intent(out)             :: VelCell(8, 3)  !<
   real(ReKi), intent(out)             :: AccCell(8, 3)  !<
   real(ReKi), intent(out)             :: Xi(3)          !< isoparametric coord of position in cell (y,z,t) [-1, +1]
   logical, intent(out)                :: Is3D           !<
   integer(IntKi), intent(out)         :: ErrStat        !< error status
   character(*), intent(out)           :: ErrMsg         !< error message

   character(*), parameter    :: RoutineName = "Grid3DField_GetSmoothInterp"
   integer(IntKi)             :: IY_Lo, IY_Hi
   integer(IntKi)             :: IZ_Lo, IZ_Hi
   integer(IntKi)             :: IT_Lo, IT_Hi
   logical                    :: OnGrid
   real(ReKi)                 :: TimeShifted
   integer(IntKi)             :: boundStat

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Find grid bounds in Time and Z
   !----------------------------------------------------------------------------

   ! Get grid time bounds
   boundStat = Grid3DField_GetBoundsT(G3D, Time, Position(1), Xi(3), IT_Lo, IT_Hi, TimeShifted)
   if (boundStat > 0) then
      ErrMsg = ' Error: GF wind array was exhausted at '//TRIM(Num2LStr(TIME))// &
               ' seconds (trying to access data at '//TRIM(Num2LStr(TimeShifted))//' seconds).'
      ErrStat = ErrID_Fatal
      return
   end if

   ! Get grid Z bounds
   boundStat = Grid3DField_GetBoundsZ(G3D, Position(3), Xi(2), IZ_Lo, IZ_Hi, OnGrid)
   if (boundStat < 0) then
      ErrMsg = ' GF wind array boundaries violated. Grid too small in Z direction '// &
               '(height (Z='//TRIM(Num2LStr(Position(3)))//' m) is below the grid and no tower points are defined).'
      ErrStat = ErrID_Fatal
      return
   else if (boundStat > 0) then
      ErrMsg = ' GF wind array boundaries violated. Grid too small in Z direction (Z='// &
               TRIM(Num2LStr(Position(3)))//' m is above the grid).'
      ErrStat = ErrID_Fatal
      return
   end if

   !----------------------------------------------------------------------------
   ! Extract cells from grids
   !----------------------------------------------------------------------------

   if (OnGrid .or. G3D%InterpTower) then

      ! Get grid Y bounds
      boundStat = Grid3DField_GetBoundsY(G3D, Position(2), Xi(1), IY_Lo, IY_Hi)
      if (boundStat /= 0) then
         ErrMsg = ' GF wind array boundaries violated: Grid too small in Y direction. Y='// &
                  TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*G3D%YHWid))// &
                  ', '//TRIM(Num2LStr(G3D%YHWid))//']'
         ErrStat = ErrID_Fatal
         return
      end if

      Is3D = .true.

      if (OnGrid) then

         ! Get velocities from the grid
         VelCell(1, :) = G3D%Vel(:, IY_Lo, IZ_Lo, IT_Lo)
         VelCell(2, :) = G3D%Vel(:, IY_Hi, IZ_Lo, IT_Lo)
         VelCell(3, :) = G3D%Vel(:, IY_Lo, IZ_Hi, IT_Lo)
         VelCell(4, :) = G3D%Vel(:, IY_Hi, IZ_Hi, IT_Lo)
         VelCell(5, :) = G3D%Vel(:, IY_Lo, IZ_Lo, IT_Hi)
         VelCell(6, :) = G3D%Vel(:, IY_Hi, IZ_Lo, IT_Hi)
         VelCell(7, :) = G3D%Vel(:, IY_Lo, IZ_Hi, IT_Hi)
         VelCell(8, :) = G3D%Vel(:, IY_Hi, IZ_Hi, IT_Hi)

         ! Get accelerations from the grid
         if (CalcAccel) then
            AccCell(1, :) = G3D%Acc(:, IY_Lo, IZ_Lo, IT_Lo)
            AccCell(2, :) = G3D%Acc(:, IY_Hi, IZ_Lo, IT_Lo)
            AccCell(3, :) = G3D%Acc(:, IY_Lo, IZ_Hi, IT_Lo)
            AccCell(4, :) = G3D%Acc(:, IY_Hi, IZ_Hi, IT_Lo)
            AccCell(5, :) = G3D%Acc(:, IY_Lo, IZ_Lo, IT_Hi)
            AccCell(6, :) = G3D%Acc(:, IY_Hi, IZ_Lo, IT_Hi)
            AccCell(7, :) = G3D%Acc(:, IY_Lo, IZ_Hi, IT_Hi)
            AccCell(8, :) = G3D%Acc(:, IY_Hi, IZ_Hi, IT_Hi)
         end if

      else if (G3D%InterpTower) then

         ! Get velocities from the grid
         VelCell(1, :) = 0.0_ReKi ! GF%Vel(:, IY_Lo, IZ_Lo, IT_Lo)
         VelCell(2, :) = 0.0_ReKi ! GF%Vel(:, IY_Hi, IZ_Lo, IT_Lo)
         VelCell(3, :) = G3D%Vel(:, IY_Lo, IZ_Hi, IT_Lo)
         VelCell(4, :) = G3D%Vel(:, IY_Hi, IZ_Hi, IT_Lo)
         VelCell(5, :) = 0.0_ReKi ! GF%Vel(:, IY_Lo, IZ_Lo, IT_Hi)
         VelCell(6, :) = 0.0_ReKi ! GF%Vel(:, IY_Hi, IZ_Lo, IT_Hi)
         VelCell(7, :) = G3D%Vel(:, IY_Lo, IZ_Hi, IT_Hi)
         VelCell(8, :) = G3D%Vel(:, IY_Hi, IZ_Hi, IT_Hi)

         ! Get accelerations from the grid
         if (CalcAccel) then
            AccCell(1, :) = 0.0_ReKi ! GF%Acc(:, IY_Lo, IZ_Lo, IT_Lo)
            AccCell(2, :) = 0.0_ReKi ! GF%Acc(:, IY_Hi, IZ_Lo, IT_Lo)
            AccCell(3, :) = G3D%Acc(:, IY_Lo, IZ_Hi, IT_Lo)
            AccCell(4, :) = G3D%Acc(:, IY_Hi, IZ_Hi, IT_Lo)
            AccCell(5, :) = 0.0_ReKi ! GF%Acc(:, IY_Lo, IZ_Lo, IT_Hi)
            AccCell(6, :) = 0.0_ReKi ! GF%Acc(:, IY_Hi, IZ_Lo, IT_Hi)
            AccCell(7, :) = G3D%Acc(:, IY_Lo, IZ_Hi, IT_Hi)
            AccCell(8, :) = G3D%Acc(:, IY_Hi, IZ_Hi, IT_Hi)
         end if

      end if

   else

      Is3D = .false.

      ! In tower grid
      if (IZ_HI <= G3D%NTGrids) then
         VelCell(1, :) = G3D%VelTower(:, IZ_LO, IT_LO)
         VelCell(2, :) = G3D%VelTower(:, IZ_HI, IT_LO)
         VelCell(3, :) = G3D%VelTower(:, IZ_LO, IT_HI)
         VelCell(4, :) = G3D%VelTower(:, IZ_HI, IT_HI)

         if (CalcAccel) then
            AccCell(1, :) = G3D%AccTower(:, IZ_LO, IT_LO)
            AccCell(2, :) = G3D%AccTower(:, IZ_HI, IT_LO)
            AccCell(3, :) = G3D%AccTower(:, IZ_LO, IT_HI)
            AccCell(4, :) = G3D%AccTower(:, IZ_HI, IT_HI)
         end if

      else  ! Between tower grid and ground

         VelCell(1, :) = G3D%VelTower(:, IZ_LO, IT_LO)
         VelCell(2, :) = 0.0_ReKi
         VelCell(3, :) = G3D%VelTower(:, IZ_LO, IT_HI)
         VelCell(4, :) = 0.0_ReKi

         if (CalcAccel) then
            AccCell(1, :) = G3D%AccTower(:, IZ_LO, IT_LO)
            AccCell(2, :) = 0.0_ReKi
            AccCell(3, :) = G3D%AccTower(:, IZ_LO, IT_HI)
            AccCell(4, :) = 0.0_ReKi
         end if

      end if

   end if

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
   real(ReKi)                 :: h, t, C1, C2, C3, C4

   ! If 3D interpolation
   if (Is3D) then

      ! Get interpolation weights
      N(1) = (1.0_ReKi - Xi(1))*(1.0_ReKi - Xi(2))/4.0_ReKi
      N(2) = (1.0_ReKi + Xi(1))*(1.0_ReKi - Xi(2))/4.0_ReKi
      N(3) = (1.0_ReKi - Xi(1))*(1.0_ReKi + Xi(2))/4.0_ReKi
      N(4) = (1.0_ReKi + Xi(1))*(1.0_ReKi + Xi(2))/4.0_ReKi

      ! Calculate velocity at lo and hi time
      if (present(Velocity)) then
         do IC = 1, 3
            P(IC, 1) = dot_product(VelCell(1:4, IC), N)  ! low time
            P(IC, 2) = dot_product(VelCell(5:8, IC), N)  ! hi time
         end do
      end if

      ! Calculate acceleration at lo and hi time
      if (present(Accel)) then
         do IC = 1, 3
            PP(IC, 1) = dot_product(AccCell(1:4, IC), N) ! low time
            PP(IC, 2) = dot_product(AccCell(5:8, IC), N) ! hi time
         end do
      end if

   else  ! 2D (Tower)

      ! Get interpolation weights
      N(1) = (1.0_ReKi - Xi(2))/2.0_ReKi
      N(2) = (1.0_ReKi + Xi(2))/2.0_ReKi

      ! Calculate velocity and acceleration at lo and hi time
      if (present(Velocity)) then
         do IC = 1, 3
            P(IC, 1) = dot_product(VelCell(1:2, IC), N(1:2))   ! low time
            P(IC, 2) = dot_product(VelCell(3:4, IC), N(1:2))   ! hi time
         end do
      end if

      ! Calculate acceleration at lo and hi time
      if (present(Accel)) then
         do IC = 1, 3
            PP(IC, 1) = dot_product(AccCell(1:2, IC), N(1:2))  ! low time
            PP(IC, 2) = dot_product(AccCell(3:4, IC), N(1:2))  ! hi time
         end do
      end if

   end if

   !----------------------------------------------------------------------------
   ! Smooth velocity and acceleration using cubic hermite spline
   !----------------------------------------------------------------------------

   h = DTime
   t = (Xi(3) + 1)/2.0_ReKi

   if (present(Velocity)) then
      C1 = 2.0_ReKi*t*t*t - 3.0_ReKi*t*t + 1.0_ReKi
      C2 = (t*t*t - 2.0_ReKi*t*t + t)*h
      C3 = -2.0_ReKi*t*t*t + 3.0_ReKi*t*t
      C4 = (t*t*t - t*t)*h
      Velocity = C1*P(:, 1) + C2*PP(:, 1) + C3*P(:, 2) + C4*PP(:, 2)
   end if

   if (present(Accel)) then
      C1 = (6.0_ReKi*t*t - 6.0_ReKi*t)/h
      C2 = 3.0_ReKi*t*t - 4.0_ReKi*t + 1.0_ReKi
      C3 = -C1
      C4 = 3.0_ReKi*t*t - 2.0_ReKi*t
      Accel = C1*P(:, 1) + C2*PP(:, 1) + C3*P(:, 2) + C4*PP(:, 2)
   end if

end subroutine

function Grid3DField_GetBoundsY(G3D, PosY, DY, IY_LO, IY_HI) result(stat)

   type(Grid3DFieldType), intent(in)   :: G3D
   real(ReKi), intent(in)              :: PosY
   real(ReKi), intent(out)             :: DY
   integer(IntKi), intent(out)         :: IY_LO, IY_HI

   real(ReKi)                          :: Y_Grid
   integer(IntKi)                      :: stat

   ! Calculate position on Y grid
   Y_Grid = (PosY + G3D%YHWid)*G3D%InvDY + 1

   ! Calculate bounding grid indices
   IY_LO = floor(Y_Grid, IntKi)
   IY_HI = ceiling(Y_Grid, IntKi)

   ! Position location within interval [0,1]
   DY = Y_Grid - aint(Y_Grid)

   ! Initialize stat to zero to indicate position is in bounds
   stat = 0

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
   else
      ! Position outside
      stat = -1
   end if

end function

function Grid3DField_GetBoundsZ(G3D, PosZ, DZ, IZ_LO, IZ_HI, OnGrid) result(stat)

   type(Grid3DFieldType), intent(in)   :: G3D
   real(ReKi), intent(in)              :: PosZ
   real(ReKi), intent(out)             :: DZ
   integer(IntKi), intent(out)         :: IZ_LO, IZ_HI
   logical, intent(out)                :: OnGrid

   integer(IntKi)                      :: stat
   real(ReKi)                          :: Z_GRID

   ! Calculate position on Z grid
   Z_GRID = (PosZ - G3D%GridBase)*G3D%InvDZ + 1

   ! Calculate bounding grid indices
   IZ_LO = floor(Z_GRID, IntKi)
   IZ_HI = ceiling(Z_GRID, IntKi)

   ! Position location within interval [-1,1]
   DZ = Z_GRID - aint(Z_GRID)

   ! Initialize stat to zero to indicate position is in bounds
   stat = 0

   ! If indices are within grid, set on grid to true
   if (IZ_LO >= 1 .and. IZ_HI <= G3D%NZGrids) then

      OnGrid = .true.
      DZ = 2.0_ReKi*DZ - 1.0_ReKi

   else if (IZ_LO < 1) then

      if (IZ_LO == 0 .and. DZ >= 1.0_ReKi - GridTol) then
         OnGrid = .true.
         IZ_LO = 1
         IZ_HI = 2
         DZ = -1.0_ReKi
      else if (G3D%InterpTower) then
         ! Interp from bottom of grid to ground (zero velocity)
         OnGrid = .false.
         IZ_LO = 0
         IZ_HI = 1
         DZ = 2.0_ReKi*(PosZ/G3D%GridBase) - 1.0_ReKi
      else if (G3D%NTGrids > 0) then
         ! Interpolate with tower grid
         OnGrid = .false.
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
      else
         ! Position below grid
         stat = -1
      end if

   else if (IZ_HI > G3D%NZGrids) then ! Above Grid

      if (IZ_HI == G3D%NZGrids + 1 .and. DZ <= GridTol) then
         OnGrid = .true.
         IZ_LO = G3D%NZGrids - 1
         IZ_HI = G3D%NZGrids
         DZ = 1.0_ReKi
      else
         ! Position above grid
         stat = 1
      end if

   end if

end function

function Grid3DField_GetBoundsT(G3D, Time, PosX, DT, IT_LO, IT_HI, TimeShifted) result(stat)

   type(Grid3DFieldType), intent(in)   :: G3D
   real(DbKi), intent(in)              :: Time
   real(ReKi), intent(in)              :: PosX
   real(ReKi), intent(out)             :: DT
   integer(IntKi), intent(out)         :: IT_LO, IT_HI
   real(ReKi), intent(out)             :: TimeShifted

   real(ReKi)                          :: T_GRID
   integer(IntKi)                      :: stat

   ! Perform the time shift. At time=0, a point half the grid width downstream
   ! (p%YHWid) will index into the zero time slice. If we did not do this,
   ! any point downstream of the tower at the beginning of the run would
   ! index outside of the array. This all assumes the grid width is at least as
   ! large as the rotor. If it isn't, then the interpolation will not work.

   ! in distance, X: InputInfo%PosX - p%InitXPosition - TIME*p%MeanWS
   TimeShifted = real(Time, ReKi) + (G3D%InitXPosition - PosX)*G3D%InvMWS

   ! If field is periodic and time is after total time, remove total time
   if (G3D%Periodic .and. TimeShifted > G3D%TotalTime) then
      TimeShifted = TimeShifted - G3D%TotalTime
   end if

   ! Get position on T grid
   T_GRID = TimeShifted*G3D%Rate + 1

   ! Calculate bounding grid indices
   IT_LO = floor(T_GRID, IntKi)
   IT_HI = ceiling(T_GRID, IntKi)

   ! Position location within interval [0,1]
   DT = T_GRID - aint(T_GRID)

   ! Initialize stat to indicate position is within grid
   stat = 0

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
      stat = 1
   end if

end function

subroutine Grid3DField_CalcAccel(G3D, ErrStat, ErrMsg)
   type(Grid3DFieldType), intent(inout)   :: G3D
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = "Grid3DField_CalcAccel"
   integer(IntKi)                         :: TmpErrStat
   character(ErrMsgLen)                   :: TmpErrMsg
   integer(IntKi)                         :: ic, iy, iz
   real(ReKi), allocatable                :: b(:), u(:), dy2(:)

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Allocate storage for acceleration grid
   call AllocAry(G3D%Acc, size(G3D%Vel, dim=1), size(G3D%Vel, dim=2), &
                 size(G3D%Vel, dim=3), size(G3D%Vel, dim=4), &
                 'grid-field velocity data', TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for B used in cubic spline derivative calc
   call AllocAry(B, G3D%NSteps, "storage for B", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for U used in cubic spline derivative calc
   call AllocAry(U, G3D%NSteps, "storage for U", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Allocate storage for V used in cubic spline derivative calc
   call AllocAry(dy2, G3D%NSteps, "storage for V", TmpErrStat, TmpErrMsg)
   call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Calculate acceleration at each grid point
   do iz = 1, G3D%NZGrids
      do iy = 1, G3D%NYGrids
         do ic = 1, G3D%NComp
            call CalcCubicSplineDeriv(G3D%DTime, G3D%Vel(ic, iy, iz, :), G3D%Acc(ic, iy, iz, :))
         end do
      end do
   end do

   ! If grid field includes tower grids
   if (G3D%NTGrids > 0) then

      ! Allocate storage for tower acceleration
      call AllocAry(G3D%AccTower, size(G3D%VelTower, dim=1), &
                    size(G3D%VelTower, dim=2), size(G3D%VelTower, dim=3), &
                    'tower wind acceleration data.', TmpErrStat, TmpErrMsg)
      call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Loop through tower grid and calculate acceleration
      do iz = 1, G3D%NTGrids
         do ic = 1, G3D%NComp
            call CalcCubicSplineDeriv(G3D%DTime, G3D%VelTower(ic, iz, :), G3D%AccTower(ic, iz, :))
         end do
      end do
   end if

contains

   !> CalcCubicSplineDeriv fits a cubic spline through the y array with points
   !! spaced a constant 'h' apart. It then calculates the corresponding
   !! derivative of y with respect to x at the same x values and returns it
   !! in the dy array.
   subroutine CalcCubicSplineDeriv(h, y, dy)
      real(ReKi), intent(in)        :: h
      real(SiKi), intent(in)        :: y(:)
      real(SiKi), intent(out)       :: dy(:)

      integer(IntKi)                :: i, n
      real(ReKi)                    :: p, un

      ! Get size of arrays
      n = size(y)

      ! If 1 or 2 points, set derivatives to zero and return
      if (n < 3) then
         do i = 1, n
            dy(i) = 0.0_ReKi
         end do
         return
      end if

      ! First derivative is zero at lower boundary condition
      dy2(1) = -0.5_ReKi
      u(1) = 3.0_ReKi*(y(2) - y(1))/h**2

      ! Calculate slopes
      do i = 1, n - 1
         b(i) = (y(i + 1) - y(i))/h
      end do

      ! Decomposition
      do i = 2, n - 1
         p = 0.5_ReKi*dy2(i - 1) + 2.0_ReKi
         dy2(i) = -0.5_ReKi/p
         u(i) = (6.*((y(i + 1) - y(i))/h - (y(i) - y(i - 1))/h)/(2.0_ReKi*h) - 0.5_ReKi*u(i - 1))/p
      end do

      ! First derviative is zero at upper boundary condition
      un = -3.0_ReKi*(y(n) - y(n - 1))/h**2
      dy2(n) = (un - 0.5_ReKi*u(n - 1))/(0.5_ReKi*dy2(n - 1) + 1.0_ReKi)

      ! Back substitution and derivative calculation
      do i = n - 1, 1, -1
         dy2(i) = dy2(i)*dy2(i + 1) + u(i)
         dy(i) = real(b(i) - h*(dy2(i)/3.0_ReKi + dy2(i + 1)/6.0_ReKi), SiKi)
      end do
      dy(n) = 0.0_ReKi

   end subroutine

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

end module
