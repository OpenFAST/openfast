module test_grid3d_field

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use testdrive, only: new_unittest, unittest_type, error_type, check
use ifw_test_tools
use IfW_FlowField
use IfW_FlowField_Types
use NWTC_Library

implicit none
private
public :: test_grid3d_field_suite

integer(IntKi), parameter :: NY = 4, NZ = 4, NT = 10
real(ReKi), parameter     :: DTIME = 0.1_ReKi

contains

!> Collect all exported unit tests
subroutine test_grid3d_field_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_grid3d_cubic_vel_only", test_grid3d_cubic_vel_only), &
               new_unittest("test_grid3d_calcaccel_no_tower", test_grid3d_calcaccel_no_tower), &
               new_unittest("test_grid3d_calcaccel_few_tower_points", test_grid3d_calcaccel_few_tower_points) &
               ]
end subroutine

!> Reproduces a bug in IfW_FlowField_GetVelAcc: when VelInterpCubic is enabled and the
!! caller does not request acceleration output (AccelUVW left unallocated), the local
!! AccCell array used by the cubic Hermite velocity formula is never populated, so the
!! returned velocity is computed from uninitialized memory. This test builds a minimal
!! Grid3D flow field with a velocity that is a known linear ramp in time (spatially
!! uniform), so the exact correct answer at any query time is known, and checks that
!! querying velocity without requesting acceleration gives the same (finite, correct)
!! answer as querying with acceleration requested.
subroutine test_grid3d_cubic_vel_only(error)
   type(error_type), allocatable, intent(out) :: error

   type(FlowFieldType)       :: FF
   real(ReKi), allocatable   :: Position(:, :), VelNoAcc(:, :), VelWithAcc(:, :)
   real(ReKi), allocatable   :: AccelUVW(:, :)
   integer(IntKi)            :: it, iy, iz
   integer(IntKi)            :: TmpErrStat
   character(ErrMsgLen)      :: TmpErrMsg
   real(DbKi)                :: QueryTime
   real(ReKi)                :: Expected

   ! Build minimal Grid3D field: NY x NZ spatial points, NT time steps, no tower grid.
   FF%FieldType = Grid3D_FieldType
   FF%VelInterpCubic = .true.
   FF%RotateWindBox = .false.

   FF%Grid3D%NComp = 3
   FF%Grid3D%NYGrids = NY
   FF%Grid3D%NZGrids = NZ
   FF%Grid3D%NTGrids = 0
   FF%Grid3D%NSteps = NT
   FF%Grid3D%DTime = DTIME
   FF%Grid3D%Rate = 1.0_ReKi/DTIME
   FF%Grid3D%YHWid = 5.0_ReKi
   FF%Grid3D%ZHWid = 5.0_ReKi
   FF%Grid3D%GridBase = 5.0_ReKi
   FF%Grid3D%InvDY = real(NY - 1, ReKi)/(2.0_ReKi*FF%Grid3D%YHWid)
   FF%Grid3D%InvDZ = real(NZ - 1, ReKi)/(2.0_ReKi*FF%Grid3D%ZHWid)
   FF%Grid3D%MeanWS = 8.0_ReKi
   FF%Grid3D%InvMWS = 1.0_ReKi/FF%Grid3D%MeanWS
   FF%Grid3D%InitXPosition = 0.0_ReKi
   FF%Grid3D%TotalTime = real(NT - 1, ReKi)*DTIME
   FF%Grid3D%Periodic = .false.
   FF%Grid3D%InterpTower = .false.

   allocate (FF%Grid3D%Vel(3, NY, NZ, NT))
   FF%Grid3D%Vel = 0.0_SiKi
   ! U component is a spatially-uniform linear ramp in time: U(t) = t (seconds -> m/s)
   do it = 1, NT
      do iz = 1, NZ
         do iy = 1, NY
            FF%Grid3D%Vel(1, iy, iz, it) = real((it - 1), SiKi)*real(DTIME, SiKi)
         end do
      end do
   end do

   ! Set the true time-derivative (dU/dt = 1, dV/dt = dW/dt = 0) directly rather than via
   ! IfW_Grid3DField_CalcAccel so this test isolates the IfW_FlowField_GetVelAcc behavior.
   ! (test_grid3d_calcaccel_no_tower below covers IfW_Grid3DField_CalcAccel itself.)
   allocate (FF%Grid3D%Acc(3, NY, NZ, NT))
   FF%Grid3D%Acc = 0.0_SiKi
   FF%Grid3D%Acc(1, :, :, :) = 1.0_SiKi
   FF%AccFieldValid = .true.

   ! Query at a position centered in the grid, at a time between samples
   allocate (Position(3, 1))
   Position(:, 1) = [0.0_ReKi, 0.0_ReKi, FF%Grid3D%GridBase + FF%Grid3D%ZHWid]
   QueryTime = 0.25_DbKi
   Expected = real(QueryTime, ReKi)

   ! Case 1: acceleration NOT requested (AccelUVW left unallocated) - buggy path
   allocate (VelNoAcc(3, 1))
   call IfW_FlowField_GetVelAcc(FF, 1, QueryTime, Position, VelNoAcc, AccelUVW, TmpErrStat, TmpErrMsg)
   call check(error, TmpErrStat, ErrID_None, message='GetVelAcc (no accel) error: '//trim(TmpErrMsg)); if (allocated(error)) return

   call check(error, ieee_is_finite(VelNoAcc(1, 1)), message='Velocity (no accel requested) is not finite (NaN/Inf)'); if (allocated(error)) return
   call check(error, VelNoAcc(1, 1), Expected, thr=1.0e-3_ReKi); if (allocated(error)) return

   ! Case 2: acceleration requested (AccelUVW allocated) - reference path
   allocate (VelWithAcc(3, 1))
   allocate (AccelUVW(3, 1))
   call IfW_FlowField_GetVelAcc(FF, 1, QueryTime, Position, VelWithAcc, AccelUVW, TmpErrStat, TmpErrMsg)
   call check(error, TmpErrStat, ErrID_None, message='GetVelAcc (with accel) error: '//trim(TmpErrMsg)); if (allocated(error)) return

   call check(error, ieee_is_finite(VelWithAcc(1, 1)), message='Velocity (accel requested) is not finite (NaN/Inf)'); if (allocated(error)) return
   call check(error, VelWithAcc(1, 1), Expected, thr=1.0e-3_ReKi); if (allocated(error)) return

   ! The two calls must agree: requesting acceleration must not change the velocity result
   call check(error, VelNoAcc(1, 1), VelWithAcc(1, 1), thr=1.0e-6_ReKi); if (allocated(error)) return

end subroutine

!> Reproduces a bug in IfW_Grid3DField_CalcAccel: it checks G3D%NTGrids (the tower-grid
!! point count) instead of G3D%NSteps (the number of time samples) to decide whether to
!! compute real cubic-spline time derivatives. For the common case of no tower file
!! (NTGrids=0), this unconditionally zeroes G3D%Acc regardless of how many time steps
!! are actually available, silently defeating cubic-in-time interpolation. This test
!! uses a spatially-uniform linear velocity ramp in time (dU/dt = 1 exactly), so the
!! correct computed derivative is known, and checks that it is recovered when there is
!! no tower grid but plenty of time steps.
subroutine test_grid3d_calcaccel_no_tower(error)
   type(error_type), allocatable, intent(out) :: error

   type(Grid3DFieldType)  :: G3D
   integer(IntKi)          :: it, iy, iz
   integer(IntKi)          :: TmpErrStat
   character(ErrMsgLen)    :: TmpErrMsg

   G3D%NComp = 3
   G3D%NYGrids = NY
   G3D%NZGrids = NZ
   G3D%NTGrids = 0
   G3D%NSteps = NT
   G3D%DTime = DTIME
   G3D%Periodic = .false.

   allocate (G3D%Vel(3, NY, NZ, NT))
   G3D%Vel = 0.0_SiKi
   do it = 1, NT
      do iz = 1, NZ
         do iy = 1, NY
            G3D%Vel(1, iy, iz, it) = real((it - 1), SiKi)*real(DTIME, SiKi)
         end do
      end do
   end do

   call IfW_Grid3DField_CalcAccel(G3D, TmpErrStat, TmpErrMsg)
   call check(error, TmpErrStat, ErrID_None, message='CalcAccel error: '//trim(TmpErrMsg)); if (allocated(error)) return

   ! Interior time points should recover close to the true derivative (dU/dt = 1); a
   ! natural cubic spline has some boundary-driven error, but zero here would indicate
   ! the NTGrids-vs-NSteps bug has regressed (computation skipped entirely).
   call check(error, real(G3D%Acc(1, 2, 2, 5), ReKi), 1.0_ReKi, thr=0.01_ReKi); if (allocated(error)) return

end subroutine

!> Reproduces the mirror-image bug in IfW_Grid3DField_CalcAccel's tower branch: it checked
!! G3D%NTGrids (tower height count) against the same <3 threshold meant for time-step counts,
!! zeroing out tower acceleration whenever there were only 1 or 2 tower grid points, even though
!! each tower height's time derivative is computed independently and only needs G3D%NSteps >= 3
!! (already guaranteed by the earlier NSteps<3 early return). This test uses NTGrids=2 with a
!! spatially-uniform linear velocity ramp in time (dU/dt = 1 exactly) and checks that the tower
!! acceleration is actually computed rather than silently zeroed.
subroutine test_grid3d_calcaccel_few_tower_points(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter :: NTWR = 2
   type(Grid3DFieldType)  :: G3D
   integer(IntKi)          :: it, iy, iz
   integer(IntKi)          :: TmpErrStat
   character(ErrMsgLen)    :: TmpErrMsg

   G3D%NComp = 3
   G3D%NYGrids = NY
   G3D%NZGrids = NZ
   G3D%NTGrids = NTWR
   G3D%NSteps = NT
   G3D%DTime = DTIME
   G3D%Periodic = .false.

   allocate (G3D%Vel(3, NY, NZ, NT))
   G3D%Vel = 0.0_SiKi
   do it = 1, NT
      do iz = 1, NZ
         do iy = 1, NY
            G3D%Vel(1, iy, iz, it) = real((it - 1), SiKi)*real(DTIME, SiKi)
         end do
      end do
   end do

   allocate (G3D%VelTower(3, NTWR, NT))
   G3D%VelTower = 0.0_SiKi
   do it = 1, NT
      do iz = 1, NTWR
         G3D%VelTower(1, iz, it) = real((it - 1), SiKi)*real(DTIME, SiKi)
      end do
   end do

   call IfW_Grid3DField_CalcAccel(G3D, TmpErrStat, TmpErrMsg)
   call check(error, TmpErrStat, ErrID_None, message='CalcAccel error: '//trim(TmpErrMsg)); if (allocated(error)) return

   call check(error, allocated(G3D%AccTower), message='AccTower was not allocated'); if (allocated(error)) return

   ! Interior time point at either tower height should recover close to the true derivative
   ! (dU/dt = 1); zero here would indicate the NTGrids<3 tower guard bug has regressed.
   call check(error, real(G3D%AccTower(1, 1, 5), ReKi), 1.0_ReKi, thr=0.01_ReKi); if (allocated(error)) return
   call check(error, real(G3D%AccTower(1, 2, 5), ReKi), 1.0_ReKi, thr=0.01_ReKi); if (allocated(error)) return

end subroutine

end module
