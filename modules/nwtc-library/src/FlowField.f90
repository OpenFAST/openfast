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

integer(IntKi), parameter  :: WindProfileType_None = -1     !< don't add wind profile; already included in input
integer(IntKi), parameter  :: WindProfileType_Constant = 0  !< constant wind
integer(IntKi), parameter  :: WindProfileType_Log = 1       !< logarithmic
integer(IntKi), parameter  :: WindProfileType_PL = 2        !< power law

integer(IntKi), parameter  :: ScaleMethod_None = 0          !< no scaling
integer(IntKi), parameter  :: ScaleMethod_Direct = 1        !< direct scaling factors
integer(IntKi), parameter  :: ScaleMethod_StdDev = 2        !< requested standard deviation

real(ReKi), parameter      :: GridTol = 1.0E-3              ! Tolerance for determining if position is within grid

contains

!> FlowField_GetVelAcc gets the velocities (and accelerations) at the given point positions.
!! Accelerations are only calculated if the AccelUVW array is allocated.
subroutine FlowField_GetVelAcc(FF, IStart, Time, PositionXYZ, VelocityUVW, AccelUVW, ErrStat, ErrMsg)
   type(FlowFieldType), intent(in)                 :: FF              !< FlowField data structure
   integer(IntKi), intent(in)                      :: IStart          !< Start index for returning velocities for external field
   real(DbKi), intent(in)                          :: Time            !< Time to evaluate velocities/accelerations
   real(ReKi), dimension(:, :), intent(in)         :: PositionXYZ     !< Array of positions to evaluate velocites/accelerations
   real(ReKi), dimension(:, :), intent(inout)      :: VelocityUVW     !< Array of velocity outputs
   real(ReKi), dimension(:, :), allocatable, intent(inout) :: AccelUVW  !< Array of acceleration outputs
   integer(IntKi), intent(out)                     :: ErrStat         !< Error status
   character(*), intent(out)                       :: ErrMsg          !< Error message

   character(*), parameter                         :: RoutineName = "FlowField_GetVelAcc"
   integer(IntKi)                                  :: i
   integer(IntKi)                                  :: NumPoints
   real(ReKi), dimension(3)                        :: PositionPrime
   type(UniformField_Interp)                       :: UFop
   type(UniformField_Interp)                       :: UFopdt
   real(ReKi)                                      :: DY, DZ, DT
   real(ReKi), dimension(3, 8)                     :: P
   logical                                         :: Interp3D
   integer(IntKi)                                  :: TmpErrStat
   character(ErrMsgLen)                            :: TmpErrMsg

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get number of points to evaluate
   NumPoints = size(PositionXYZ, dim=2)

   !TODO: Check that position, velocity, and acceleration are all the same shape


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

      ! Get uniform flow operating point values at current time
      UFop = UniformField_GetOP(FF%Uniform, Time)

      ! If getting accel, get change in operating point values wrt time
      if (allocated(AccelUVW)) then
         UFopdt%VelH = (FF%Uniform%VelH(Ufop%IData) - FF%Uniform%VelH(Ufop%IData - 1))/Ufop%dt
         UFopdt%AngleH = (FF%Uniform%AngleH(Ufop%IData) - FF%Uniform%AngleH(Ufop%IData - 1))/Ufop%dt
         UFopdt%AngleV = (FF%Uniform%AngleV(Ufop%IData) - FF%Uniform%AngleV(Ufop%IData - 1))/Ufop%dt
         UFopdt%VelV = (FF%Uniform%VelV(Ufop%IData) - FF%Uniform%VelV(Ufop%IData - 1))/Ufop%dt
         UFopdt%ShrH = (FF%Uniform%ShrH(Ufop%IData) - FF%Uniform%ShrH(Ufop%IData - 1))/Ufop%dt
         UFopdt%ShrV = (FF%Uniform%ShrV(Ufop%IData) - FF%Uniform%ShrV(Ufop%IData - 1))/Ufop%dt
         UFopdt%LinShrV = (FF%Uniform%LinShrV(Ufop%IData) - FF%Uniform%LinShrV(Ufop%IData - 1))/Ufop%dt
         UFopdt%VelGust = (FF%Uniform%VelGust(Ufop%IData) - FF%Uniform%VelGust(Ufop%IData - 1))/Ufop%dt
      end if

      if (.not. FF%RotateWindBox) then
         if (.not. allocated(AccelUVW)) then
            ! Velocity
            do i = 1, NumPoints
               VelocityUVW(:, i) = UniformField_GetVel(FF%Uniform, UFop, PositionXYZ(:, i))
            end do
         else
            ! Velocity and acceleration
            do i = 1, NumPoints
               VelocityUVW(:, i) = UniformField_GetVel(FF%Uniform, UFop, PositionXYZ(:, i))
               AccelUVW(:, i) = UniformField_GetAcc(FF%Uniform, UFop, UFopdt, PositionXYZ(:, i))
            end do
         end if

      else ! Rotated wind box

         if (.not. allocated(AccelUVW)) then
            ! Rotation, velocity
            do i = 1, NumPoints
               PositionPrime = GetPrimePosition(PositionXYZ(:, i))
               VelocityUVW(:, i) = UniformField_GetVel(FF%Uniform, UFop, PositionPrime)
            end do
         else
            ! Rotation, velocity, and acceleration
            do i = 1, NumPoints
               PositionPrime = GetPrimePosition(PositionXYZ(:, i))
               VelocityUVW(:, i) = UniformField_GetVel(FF%Uniform, UFop, PositionPrime)
               AccelUVW(:, i) = UniformField_GetAcc(FF%Uniform, UFop, UFopdt, PositionPrime)
            end do
         end if

      end if

   case (Grid_FieldType)

      !-------------------------------------------------------------------------
      ! Grid Flow Field
      !-------------------------------------------------------------------------

      if (.not. FF%RotateWindBox) then

         if (.not. allocated(AccelUVW)) then
            ! Velocity
            do i = 1, NumPoints
               if (PositionXYZ(3, i) <= 0.0_ReKi) then
                  VelocityUVW(:, i) = 0.0_ReKi
                  cycle
               end if
               call GridField_GetInterp(FF%Grid, Time, PositionXYZ(:, i), DY, DZ, DT, P, Interp3D, TmpErrStat, TmpErrMsg)
               if (TmpErrStat >= AbortErrLev) then
                  call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
                  return
               end if
               VelocityUVW(:, i) = GridField_GetVel(FF%Grid, PositionXYZ(:, i), DY, DZ, DT, P, Interp3D)
            end do
         else
            ! Velocity and acceleration
            do i = 1, NumPoints
               if (PositionXYZ(3, i) <= 0.0_ReKi) then
                  VelocityUVW(:, i) = 0.0_ReKi
                  AccelUVW(:, i) = 0.0_Reki
                  cycle
               end if
               call GridField_GetInterp(FF%Grid, Time, PositionXYZ(:, i), DY, DZ, DT, P, Interp3D, TmpErrStat, TmpErrMsg)
               if (TmpErrStat >= AbortErrLev) then
                  call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
                  return
               end if
               VelocityUVW(:, i) = GridField_GetVel(FF%Grid, PositionXYZ(:, i), DY, DZ, DT, P, Interp3D)
               AccelUVW(:, i) = GridField_GetAcc(FF%Grid, PositionXYZ(:, i), DY, DZ, DT, P, Interp3D)
            end do
         end if

      else ! Rotated wind box

         if (.not. allocated(AccelUVW)) then
            ! Velocity (rotated)
            do i = 1, NumPoints
               PositionPrime = GetPrimePosition(PositionXYZ(:, i))
               if (PositionPrime(3) <= ReKi) then
                  VelocityUVW(:, i) = 0.0_ReKi
                  cycle
               end if
               call GridField_GetInterp(FF%Grid, Time, PositionPrime, DY, DZ, DT, P, Interp3D, TmpErrStat, TmpErrMsg)
               if (TmpErrStat >= AbortErrLev) then
                  call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
                  return
               end if
               VelocityUVW(:, i) = matmul(FF%RotFromWind, GridField_GetVel(FF%Grid, PositionPrime, DY, DZ, DT, P, Interp3D))
            end do
         else
            ! Velocity and acceleration (rotated)
            do i = 1, NumPoints
               ! Get rotated coordinates
               PositionPrime = GetPrimePosition(PositionXYZ(:, i))
               if (PositionPrime(3) <= ReKi) then
                  VelocityUVW(:, i) = 0.0_ReKi
                  AccelUVW(:, i) = 0.0_ReKi
                  cycle
               end if
               call GridField_GetInterp(FF%Grid, Time, PositionPrime, DY, DZ, DT, P, Interp3D, TmpErrStat, TmpErrMsg)
               if (TmpErrStat >= AbortErrLev) then
                  call SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
                  return
               end if
               VelocityUVW(:, i) = matmul(FF%RotFromWind, GridField_GetVel(FF%Grid, PositionPrime, DY, DZ, DT, P, Interp3D))
               AccelUVW(:, i) = matmul(FF%RotFromWind, GridField_GetAcc(FF%Grid, PositionPrime, DY, DZ, DT, P, Interp3D))
            end do
         end if
      end if

   case (User_FieldType)

      !-------------------------------------------------------------------------
      ! User Flow Field
      !-------------------------------------------------------------------------

      call SetErrStat(ErrID_Fatal, "User Field not implemented", ErrStat, ErrMsg, RoutineName)
      return

   case (External_FieldType)

      !-------------------------------------------------------------------------
      ! External Flow Field
      !-------------------------------------------------------------------------

      ! If external field is not allocated, return error
      if (.not. allocated(FF%external%Vel)) then
         call SetErrStat(ErrID_Fatal, "External Field not allocated", ErrStat, ErrMsg, RoutineName)
         return
      end if

      VelocityUVW = FF%external%Vel(:, IStart:IStart + NumPoints - 1)

   case default
      call SetErrStat(ErrID_Fatal, "Invalid FieldType "//trim(num2lstr(FF%FieldType)), ErrStat, ErrMsg, RoutineName)
      return
   end select

contains

   pure function GetPrimePosition(Position) result(PrimePosition)
      real(ReKi), dimension(3), intent(in)      :: Position
      real(ReKi), dimension(3)                  :: PrimePosition
      PrimePosition = matmul(FF%RotToWind, (Position - FF%RefPosition)) + FF%RefPosition
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

   if (Position(3) > 0.0_ReKi) then

      ! Calculate horizontal velocity if position is above ground
      V1 = op%VelH*((Position(3)/UF%RefHeight)**op%ShrV &                                 ! power-law wind shear
                    + (op%ShrH*(Position(2)*op%CosAngleH + Position(1)*op%SinAngleH) & ! horizontal linear shear
                       + op%LinShrV*(Position(3) - UF%RefHeight))/UF%RefLength) &         ! vertical linear shear
           + op%VelGust                                                                      ! gust speed

      ! Apply upflow angle
      V1_rotate = op%CosAngleV*V1 - op%SinAngleV*op%VelV
      VZ_rotate = op%SinAngleV*V1 + op%CosAngleV*op%VelV

      ! Apply wind direction
      Velocity(1) = V1_rotate*op%CosAngleH
      Velocity(2) = -V1_rotate*op%SinAngleH
      Velocity(3) = VZ_rotate

   else

      ! Height is below ground, set velocity to zero.
      ! Removed error check to be consistent with full field method
      Velocity = 0.0_ReKi

   end if

end function

function UniformField_GetAcc(UF, op, opdt, Position) result(Accel)
   type(UniformFieldType), intent(in)     :: UF
   type(UniformField_Interp), intent(in)  :: op
   type(UniformField_Interp), intent(in)  :: opdt
   real(ReKi), dimension(3), intent(in)   :: Position
   real(ReKi), dimension(3)               :: Accel

   character(*), parameter                :: RoutineName = "UniformField_GetAcc"
   real(ReKi)                             :: C1, C2, C3, C4, C5

   ! If only one data point or time is outside data range, acceleration is zero
   if (op%IData == 1 .or. op%IData > UF%DataSize .or. Position(3) <= 0.0_Reki) then
      Accel = 0.0_ReKi
   else

      ! If change in time is negative or near zero, acceleration is zero
      if (op%dt < 0.0_ReKi .or. EqualRealNos(op%dt, 0.0_ReKi)) then
         Accel = 0.0_ReKi
         return
      end if

      C1 = (Position(3)/UF%RefHeight)**op%ShrV + &
           (op%LinShrV*(Position(3) - UF%RefHeight) + &
            op%ShrH*(Position(1)*op%SinAngleH + &
                     Position(2)*op%CosAngleH))/UF%RefLength

      C2 = op%CosAngleV*(op%VelGust + op%VelH*(C1))

      C3 = opdt%AngleV*op%SinAngleV*(op%VelGust + op%VelH*(C1))

      C4 = opdt%LinShrV*(Position(3) - UF%RefHeight) + &
           opdt%ShrH*(Position(1)*op%SinAngleH + Position(2)*op%CosAngleH) + &
           op%ShrH*(Position(1)*opdt%AngleH*op%CosAngleH - &
                    Position(2)*opdt%AngleH*op%SinAngleH)

      C5 = opdt%VelGust + opdt%VelH*C1 + &
           op%VelH*(opdt%ShrV*(Position(3)/UF%RefHeight)**op%ShrV* &
                    log(Position(3)/UF%RefHeight) + C4/UF%RefLength)

      Accel(1) = -opdt%AngleH*op%SinAngleH*(C2 - op%SinAngleV*op%VelV) + &
                 op%CosAngleH*(-opdt%AngleV*op%CosAngleV*op%VelV - C3 - &
                               opdt%VelV*op%SinAngleV + op%CosAngleV*C5)

      Accel(2) = opdt%AngleH*op%CosAngleH*(-C2 + op%SinAngleV*op%VelV) + &
                 op%SinAngleH*(opdt%AngleV*op%CosAngleV*op%VelV + C3 + &
                               opdt%VelV*op%SinAngleV - op%CosAngleV*C5)

      Accel(3) = opdt%AngleV*C2 - opdt%AngleV*op%SinAngleV*op%VelV + &
                 opdt%VelV*op%CosAngleV + op%SinAngleV*C5
   end if

end function

pure function UniformField_GetOP(UF, Time) result(op)

   type(UniformFieldType), intent(in)  :: UF
   real(DbKi), intent(in)              :: Time
   type(UniformField_Interp)           :: op
   integer(IntKi)                      :: i
   real(ReKi)                          :: alpha, OMalpha

   ! Save the time in operating point
   op%Time = real(Time, ReKi)

   ! If only one data point or time is at or below first time, use first sample
   if (UF%DataSize == 1 .or. Time < UF%Time(1)) then

      op%IData = 1
      op%dt = 0.0_ReKi
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

      op%IData = UF%DataSize + 1
      op%dt = 0.0_ReKi
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

      op%IData = i

      ! Calculate interval delta time
      op%dt = UF%Time(i) - UF%Time(i - 1)

      ! Calculate interpolation coefficient [0,1]
      alpha = real((Time - UF%Time(i - 1))/op%dt, ReKi)
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

pure function UniformField_AD_DiskVel(op) result(DiskVel)
   type(UniformField_Interp), intent(in)   :: op
   real(ReKi), dimension(3)               :: DiskVel

   real(ReKi)                             :: V1_rotate, Vz_rotate

   V1_rotate = op%CosAngleV*op%VelH - op%SinAngleV*op%VelV
   Vz_rotate = op%SinAngleV*op%VelH + op%CosAngleV*op%VelV

   DiskVel(1) = V1_rotate*op%CosAngleH
   DiskVel(2) = -V1_rotate*op%SinAngleH
   DiskVel(3) = Vz_rotate

end function

function GridField_GetVel(GF, Position, DY, DZ, DT, P, Interp3D) result(Velocity)

   type(GridFieldType), intent(in)        :: GF             !< Grid-Field data
   real(ReKi), intent(in)                 :: Position(3)    !< position X,Y,Z to get value
   real(ReKi), intent(in)                 :: DY, DZ, DT     !< Weights [-1,1]
   real(ReKi), dimension(3, 8), intent(in) :: P              !< Interp points
   logical, intent(in)                    :: Interp3D       !< flag for 3D vs 2D interp
   real(ReKi)                             :: Velocity(3)    !< The U, V, W velocities

   character(*), parameter                :: RoutineName = "GridField_GetVel"

   real(ReKi)                             :: N(8, 1)     ! Shape function values

   if (Interp3D) then

      ! Get 3D interpolation weights
      N(1, 1) = (1.0_ReKi - DY)*(1.0_ReKi - DZ)
      N(2, 1) = (1.0_ReKi + DY)*(1.0_ReKi - DZ)
      N(3, 1) = (1.0_ReKi - DY)*(1.0_ReKi + DZ)
      N(4, 1) = (1.0_ReKi + DY)*(1.0_ReKi + DZ)
      N(5, 1) = (1.0_ReKi - DY)*(1.0_ReKi - DZ)
      N(6, 1) = (1.0_ReKi + DY)*(1.0_ReKi - DZ)
      N(7, 1) = (1.0_ReKi - DY)*(1.0_ReKi + DZ)
      N(8, 1) = (1.0_ReKi + DY)*(1.0_ReKi + DZ)
      N(1:4, 1) = N(1:4, 1)*(1.0_ReKi - DT)/8.0_ReKi
      N(5:8, 1) = N(5:8, 1)*(1.0_ReKi + DT)/8.0_ReKi

      ! Calculate velocity via matrix multiplication
      Velocity = pack(matmul(P, N), .true.)

   else

      ! Get 2D interpolation weights
      N(1, 1) = (1.0_ReKi - DZ)*(1.0_ReKi - DT)/4.0_ReKi
      N(2, 1) = (1.0_ReKi + DZ)*(1.0_ReKi - DT)/4.0_ReKi
      N(3, 1) = (1.0_ReKi - DZ)*(1.0_ReKi + DT)/4.0_ReKi
      N(4, 1) = (1.0_ReKi + DZ)*(1.0_ReKi + DT)/4.0_ReKi

      ! Calculate velocity via matrix multiplication
      Velocity = pack(matmul(P(:, 1:4), N(1:4, 1)), .true.)

   end if

   ! Add mean wind speed after interpolation if flag is set
   if (GF%AddMeanAfterInterp) then
      Velocity(1) = Velocity(1) + GetMeanVelocity(Position(3))
   end if

contains

   function GetMeanVelocity(PosZ) result(U)
      real(ReKi), intent(in)   :: PosZ
      real(ReKi)              :: U
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

end function

function GridField_GetAcc(GF, Position, DY, DZ, DT, P, Interp3D) result(Accel)

   type(GridFieldType), intent(in)           :: GF
   real(ReKi), dimension(3), intent(in)      :: Position
   real(ReKi), intent(in)                    :: DY, DZ, DT     !< Weights [-1,1]
   real(ReKi), dimension(3, 8), intent(in)   :: P              !< Interp points
   logical, intent(in)                       :: Interp3D       !< flag for 3D vs 2D interp
   real(ReKi), dimension(3)                  :: Accel

   character(*), parameter                   :: RoutineName = "GridField_GetAcc"
   real(ReKi)                                :: N(8, 1)

   if (Interp3D) then

      ! Get 3D interpolation weights
      N(1, 1) = (1.0_ReKi - DY)*(1.0_ReKi - DZ)
      N(2, 1) = (1.0_ReKi + DY)*(1.0_ReKi - DZ)
      N(3, 1) = (1.0_ReKi - DY)*(1.0_ReKi + DZ)
      N(4, 1) = (1.0_ReKi + DY)*(1.0_ReKi + DZ)
      N(5, 1) = (1.0_ReKi - DY)*(1.0_ReKi - DZ)
      N(6, 1) = (1.0_ReKi + DY)*(1.0_ReKi - DZ)
      N(7, 1) = (1.0_ReKi - DY)*(1.0_ReKi + DZ)
      N(8, 1) = (1.0_ReKi + DY)*(1.0_ReKi + DZ)
      N(1:4, 1) = N(1:4, 1)/(-4.0_ReKi*GF%DTime)
      N(5:8, 1) = N(5:8, 1)/(4.0_ReKi*GF%DTime)

      ! Calculate Accel via matrix multiplication
      Accel = pack(matmul(P, N), .true.)

   else

      ! Get 2D interpolation weights
      N(1, 1) = (1.0_ReKi - DZ)/(-2.0_ReKi*GF%DTime)
      N(2, 1) = (1.0_ReKi + DZ)/(-2.0_ReKi*GF%DTime)
      N(3, 1) = (1.0_ReKi - DZ)/(2.0_ReKi*GF%DTime)
      N(4, 1) = (1.0_ReKi + DZ)/(2.0_ReKi*GF%DTime)

      ! Calculate Accel via matrix multiplication
      Accel = pack(matmul(P(:, 1:4), N(1:4, 1)), .true.)

   end if

end function

subroutine GridField_GetInterp(GF, Time, Position, DY, DZ, DT, P, Interp3D, ErrStat, ErrMsg)

   type(GridFieldType), intent(in)     :: GF                !< Grid-Field data
   real(DbKi), intent(in)              :: Time              !< time (s)
   real(ReKi), intent(in)              :: Position(3)       !< position X,Y,Z to get value
   real(ReKi), intent(out)             :: DY, DZ, DT        !<
   real(ReKi), intent(out)             :: P(3, 8)
   logical, intent(out)                :: Interp3D
   integer(IntKi), intent(out)         :: ErrStat           !< error status
   character(*), intent(out)           :: ErrMsg            !< error message

   character(*), parameter             :: RoutineName = "GridField_GetInterp"

   integer(IntKi)                      :: IT_LO, IT_HI
   integer(IntKi)                      :: IY_LO, IY_HI
   integer(IntKi)                      :: IZ_LO, IZ_HI

   logical                             :: OnGrid
   real(ReKi)                          :: TimeShifted
   integer(IntKi)                      :: boundStat

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Find grid bounds in Time and Z
   !----------------------------------------------------------------------------

   ! Get grid time bounds
   boundStat = GridField_GetBoundsT(GF, Time, Position(1), DT, IT_LO, IT_HI, TimeShifted)
   if (boundStat > 0) then
      ErrMsg = ' Error: GF wind array was exhausted at '//TRIM(Num2LStr(TIME))// &
               ' seconds (trying to access data at '//TRIM(Num2LStr(TimeShifted))//' seconds).'
      ErrStat = ErrID_Fatal
      return
   end if

   ! Get grid Z bounds
   boundStat = GridField_GetBoundsZ(GF, Position(3), DZ, IZ_LO, IZ_HI, OnGrid)
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
   ! Interpolate
   !----------------------------------------------------------------------------

   if (OnGrid .or. GF%InterpTower) then

      ! Set flag for 3D interpolation to true
      Interp3D = .true.

      ! Get grid Y bounds
      boundStat = GridField_GetBoundsY(GF, Position(2), DY, IY_LO, IY_HI)
      if (boundStat /= 0) then
         ErrMsg = ' GF wind array boundaries violated: Grid too small in Y direction. Y='// &
                  TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*GF%YHWid))// &
                  ', '//TRIM(Num2LStr(GF%YHWid))//']'
         ErrStat = ErrID_Fatal
         return
      end if

      if (OnGrid) then

         ! Get points from the grid
         P(:, 1) = GF%Vel(:, IY_LO, IZ_LO, IT_LO)
         P(:, 2) = GF%Vel(:, IY_HI, IZ_LO, IT_LO)
         P(:, 3) = GF%Vel(:, IY_LO, IZ_HI, IT_LO)
         P(:, 4) = GF%Vel(:, IY_HI, IZ_HI, IT_LO)
         P(:, 5) = GF%Vel(:, IY_LO, IZ_LO, IT_HI)
         P(:, 6) = GF%Vel(:, IY_HI, IZ_LO, IT_HI)
         P(:, 7) = GF%Vel(:, IY_LO, IZ_HI, IT_HI)
         P(:, 8) = GF%Vel(:, IY_HI, IZ_HI, IT_HI)

      else if (GF%InterpTower) then

         ! Get points from grid bottom and ground
         P(:, 1) = 0.0_ReKi !GF%Vel(:, IY_LO, IZ_LO, IT_LO)
         P(:, 2) = 0.0_ReKi !GF%Vel(:, IY_HI, IZ_LO, IT_LO)
         P(:, 3) = GF%Vel(:, IY_LO, IZ_HI, IT_LO)
         P(:, 4) = GF%Vel(:, IY_HI, IZ_HI, IT_LO)
         P(:, 5) = 0.0_ReKi !GF%Vel(:, IY_HI, IZ_LO, IT_HI)
         P(:, 6) = 0.0_ReKi !GF%Vel(:, IY_LO, IZ_LO, IT_HI)
         P(:, 7) = GF%Vel(:, IY_LO, IZ_HI, IT_HI)
         P(:, 8) = GF%Vel(:, IY_HI, IZ_HI, IT_HI)

      end if

   else

      ! Set flag for 3D interpolation to false for 2D
      Interp3D = .false.

      if (IZ_HI <= GF%NTGrids) then
         ! In tower grid
         P(:, 1) = GF%VelTower(:, IZ_LO, IT_LO)
         P(:, 2) = GF%VelTower(:, IZ_HI, IT_LO)
         P(:, 3) = GF%VelTower(:, IZ_LO, IT_HI)
         P(:, 4) = GF%VelTower(:, IZ_HI, IT_HI)
      else
         ! Between tower grid and ground
         P(:, 1) = GF%VelTower(:, IZ_LO, IT_LO)
         P(:, 2) = 0.0_ReKi
         P(:, 3) = GF%VelTower(:, IZ_LO, IT_HI)
         P(:, 4) = 0.0_ReKi
      end if

   end if

end subroutine

function GridField_GetBoundsY(GF, PosY, DY, IY_LO, IY_HI) result(stat)

   type(GridFieldType), intent(in)     :: GF
   real(ReKi), intent(in)              :: PosY
   real(ReKi), intent(out)             :: DY
   integer(IntKi), intent(out)         :: IY_LO, IY_HI

   real(ReKi)                          :: Y_Grid
   integer(IntKi)                      :: stat

   ! Calculate position on Y grid
   Y_Grid = (PosY + GF%YHWid)*GF%InvDY + 1

   ! Calculate bounding grid indices
   IY_LO = floor(Y_Grid, IntKi)
   IY_HI = ceiling(Y_Grid, IntKi)

   ! Position location within interval [0,1]
   DY = Y_Grid - aint(Y_Grid)

   ! Initialize stat to zero to indicate position is in bounds
   stat = 0

   if (IY_LO >= 1 .and. IY_HI <= GF%NYGrids) then
      DY = 2.0_ReKi*DY - 1.0_ReKi
   else if (IY_LO == 0 .and. DY >= 1.0_ReKi - GridTol) then
      IY_LO = 1
      IY_HI = 2
      DY = -1.0_ReKi
   else if (IY_LO == GF%NYGrids .and. DY <= GridTol) then
      IY_LO = GF%NYGrids - 1
      IY_HI = GF%NYGrids
      DY = 1.0_ReKi
   else
      ! Position outside
      stat = -1
   end if

end function

function GridField_GetBoundsZ(GF, PosZ, DZ, IZ_LO, IZ_HI, OnGrid) result(stat)

   type(GridFieldType), intent(in)     :: GF
   real(ReKi), intent(in)              :: PosZ
   real(ReKi), intent(out)             :: DZ
   integer(IntKi), intent(out)         :: IZ_LO, IZ_HI
   logical, intent(out)                :: OnGrid

   integer(IntKi)                      :: stat
   real(ReKi)                          :: Z_GRID
   real(ReKi)                          :: Z_twr(3)

   ! Calculate position on Z grid
   Z_GRID = (PosZ - GF%GridBase)*GF%InvDZ + 1

   ! Calculate bounding grid indices
   IZ_LO = floor(Z_GRID, IntKi)
   IZ_HI = ceiling(Z_GRID, IntKi)

   ! Position location within interval [-1,1]
   DZ = Z_GRID - aint(Z_GRID)

   ! Initialize stat to zero to indicate position is in bounds
   stat = 0

   ! If indices are within grid, set on grid to true
   if (IZ_LO >= 1 .and. IZ_HI <= GF%NZGrids) then

      OnGrid = .true.
      DZ = 2.0_ReKi*DZ - 1.0_ReKi

   else if (IZ_LO < 1) then

      if (IZ_LO == 0 .and. DZ >= 1.0_ReKi - GridTol) then
         OnGrid = .true.
         IZ_LO = 1
         IZ_HI = 2
         DZ = -1.0_ReKi
      else if (GF%InterpTower) then
         ! Interp from bottom of grid to ground (zero velocity)
         OnGrid = .false.
         IZ_LO = 0
         IZ_HI = 1
         DZ = 2.0_ReKi*(PosZ/GF%GridBase) - 1.0_ReKi
      else if (GF%NTGrids > 0) then
         ! Interpolate with tower grid
         OnGrid = .false.
         ! Tower grid is reversed (lowest index is top of tower)
         IZ_LO = int(-(Z_GRID - 1)) + 1
         if (IZ_LO >= GF%NTGrids) then
            ! Between end of tower grid and ground (zero velocity)
            IZ_LO = GF%NTGrids
            DZ = 1.0_ReKi - 2.0_ReKi*(PosZ/(GF%GridBase - real(IZ_LO - 1, ReKi)/GF%InvDZ))
         else
            ! Within tower grid
            DZ = 2.0_ReKi*(real(2 - IZ_LO, ReKi) - Z_GRID) - 1.0_ReKi
         end if
         IZ_HI = IZ_LO + 1
      else
         ! Position below grid
         stat = -1
      end if

   else if (IZ_HI > GF%NZGrids) then ! Above Grid

      if (IZ_HI == GF%NZGrids + 1 .and. DZ <= GridTol) then
         OnGrid = .true.
         IZ_LO = GF%NZGrids - 1
         IZ_HI = GF%NZGrids
         DZ = 1.0_ReKi
      else
         ! Position above grid
         stat = 1
      end if

   end if

end function

function GridField_GetBoundsT(GF, Time, PosX, DT, IT_LO, IT_HI, TimeShifted) result(stat)

   type(GridFieldType), intent(in)     :: GF
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
   TimeShifted = real(Time, ReKi) + (GF%InitXPosition - PosX)*GF%InvMWS

   ! If field is periodic and time is after total time, remove total time
   if (GF%Periodic .and. TimeShifted > GF%TotalTime) then
      TimeShifted = TimeShifted - GF%TotalTime
   end if

   ! Get position on T grid
   T_GRID = TimeShifted*GF%Rate + 1

   ! Calculate bounding grid indices
   IT_LO = floor(T_GRID, IntKi)
   IT_HI = ceiling(T_GRID, IntKi)

   ! Position location within interval [0,1]
   DT = T_GRID - aint(T_GRID)

   ! Initialize stat to indicate position is within grid
   stat = 0

   ! Adjust indices and interpolant
   if (IT_LO >= 1 .and. IT_HI <= GF%NSteps) then
      ! Point is within grid
      DT = 2.0_ReKi*DT - 1.0_ReKi
   else if (IT_LO == GF%NSteps) then
      if (GF%Periodic) then
         ! Time wraps back to beginning
         IT_HI = 1
         DT = 2.0_ReKi*DT - 1.0_ReKi
      else if (DT <= GridTol) then
         ! Within tolerance of last time
         IT_HI = IT_LO
         DT = -1.0_Reki
      else
         ! Extrapolate
         IT_LO = GF%NSteps - 1
         IT_HI = GF%NSteps
         DT = DT + 1.0_ReKi
      end if
   else
      ! Time exceeds array bounds
      stat = 1
   end if

end function

end module
