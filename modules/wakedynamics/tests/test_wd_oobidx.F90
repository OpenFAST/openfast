module test_wd_oobidx

use testdrive, only: new_unittest, unittest_type, error_type, check
use WakeDynamics
use WakeDynamics_Types
use NWTC_Library

implicit none
private
public :: test_wd_oobidx_suite

contains

!> Collect all exported unit tests
subroutine test_wd_oobidx_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_all_planes_simultaneously_out_of_bounds", test_all_planes_simultaneously_out_of_bounds) &
               ]
end subroutine

!> Regression test for a stack array-bounds overflow in WD_UpdateStates: when every
!! existing wake plane (up to p%MaxNumPlanes of them) is simultaneously found to be
!! outside the low-resolution domain bounds, the local automatic array `oobIdx` must
!! have room to record all of them. Before the fix, `oobIdx` was declared
!! `0:MaxNumPlanes-1` (size MaxNumPlanes) -- one element too small to hold
!! `MaxNumPlanes` out-of-bounds indices -- causing an out-of-bounds stack write.
!! (With gfortran's `-fcheck=all` Debug flag, this manifests as a hard runtime
!! crash rather than a graceful test failure.)
subroutine test_all_planes_simultaneously_out_of_bounds(error)
   type(error_type), allocatable, intent(out) :: error

   type(WD_InitInputType)       :: InitInp
   type(WD_InputType)           :: u
   type(WD_ParameterType)       :: p
   type(WD_ContinuousStateType) :: x
   type(WD_DiscreteStateType)   :: xd
   type(WD_ConstraintStateType) :: z
   type(WD_OtherStateType)      :: OtherState
   type(WD_OutputType)          :: y
   type(WD_MiscVarType)         :: m
   type(WD_InitOutputType)      :: InitOut
   integer(IntKi)               :: errStat
   character(ErrMsgLen)         :: errMsg
   real(DbKi), parameter        :: DT_low = 1.0_DbKi
   integer(IntKi), parameter    :: MaxNumPlanes = 4
   integer(IntKi)                :: n

   errStat = ErrID_None
   errMsg = ''

   ! --- Minimal, valid WakeDynamics (Polar model) initialization input
   InitInp%TurbNum = 1
   InitInp%OutFileRoot = 'test_wd_oobidx'
   InitInp%MaxNumPlanes = MaxNumPlanes
   ! Generous low-res domain bounds so no plane is dropped during natural growth
   InitInp%LowResBounds(:, 1) = -1.0e6_ReKi
   InitInp%LowResBounds(:, 2) = 1.0e6_ReKi

   InitInp%InputFileData%Mod_Wake = Mod_Wake_Polar
   InitInp%InputFileData%RotorDiamRef = 20.0_ReKi
   InitInp%InputFileData%dr = 5.0_ReKi
   InitInp%InputFileData%NumRadii = 3
   InitInp%InputFileData%NumDFull = 15.0_ReKi
   InitInp%InputFileData%NumDBuff = 5.0_ReKi
   InitInp%InputFileData%f_c = 0.17_ReKi
   InitInp%InputFileData%C_NearWake = 1.8_ReKi
   InitInp%InputFileData%k_vAmb = 0.05_ReKi
   InitInp%InputFileData%C_vAmb_FMin = 1.0_ReKi
   InitInp%InputFileData%C_vAmb_DMin = 0.0_ReKi
   InitInp%InputFileData%C_vAmb_DMax = 1.0_ReKi
   InitInp%InputFileData%C_vAmb_Exp = 0.01_ReKi
   InitInp%InputFileData%k_vShr = 0.016_ReKi
   InitInp%InputFileData%C_vShr_FMin = 0.2_ReKi
   InitInp%InputFileData%C_vShr_DMin = 3.0_ReKi
   InitInp%InputFileData%C_vShr_DMax = 25.0_ReKi
   InitInp%InputFileData%C_vShr_Exp = 0.1_ReKi
   InitInp%InputFileData%Mod_WakeDiam = 1

   call WD_Init(InitInp, u, p, x, xd, z, OtherState, y, m, DT_low, InitOut, errStat, errMsg)
   call check(error, errStat < AbortErrLev, "WD_Init failed: "//trim(errMsg))
   if (allocated(error)) return

   ! --- Steady, uniform axial inflow so plane advection is simple and predictable
   u%xhat_disk = (/1.0_ReKi, 0.0_ReKi, 0.0_ReKi/)
   u%YawErr = 0.0_ReKi
   u%psi_skew = 0.0_ReKi
   u%chi_skew = 0.0_ReKi
   u%p_hub = (/0.0_ReKi, 0.0_ReKi, 90.0_ReKi/)
   u%Vx_wind_disk = 8.0_ReKi
   u%TI_amb = 0.1_ReKi
   u%D_rotor = 20.0_ReKi
   u%Vx_rel_disk = 8.0_ReKi
   u%Ct_azavg = 0.0_ReKi
   u%Cq_azavg = 0.0_ReKi
   u%V_plane(1, :) = 8.0_ReKi
   u%V_plane(2, :) = 0.0_ReKi
   u%V_plane(3, :) = 0.0_ReKi

   ! --- Grow the number of tracked wake planes up to p%MaxNumPlanes
   do n = 0, MaxNumPlanes - 3
      call WD_UpdateStates(real(n, DbKi)*DT_low, n, u, p, x, xd, z, OtherState, m, errStat, errMsg)
      call check(error, errStat < AbortErrLev, "WD_UpdateStates failed during plane growth: "//trim(errMsg))
      if (allocated(error)) return
   end do
   call check(error, NINT(xd%NumPlanes) == p%MaxNumPlanes, "Expected wake planes to grow to MaxNumPlanes")
   if (allocated(error)) return

   ! --- Tighten the low-resolution domain bounds so every existing plane
   !     (including the disk plane, which stays near Z=90) is now out of bounds.
   p%LowResBounds(:, 1) = 500.0_ReKi
   p%LowResBounds(:, 2) = 600.0_ReKi

   ! With the fix, this call completes without an out-of-bounds write to the
   ! oobIdx local array (previously sized one element too small to hold
   ! MaxNumPlanes out-of-bounds indices).
   call WD_UpdateStates(real(MaxNumPlanes - 2, DbKi)*DT_low, MaxNumPlanes - 2, u, p, x, xd, z, OtherState, m, errStat, errMsg)
   call check(error, errStat < AbortErrLev, &
              "WD_UpdateStates failed when all planes were simultaneously out of bounds: "//trim(errMsg))
   if (allocated(error)) return

end subroutine

end module
