module test_wd_maxplanes_warning

use testdrive, only: new_unittest, unittest_type, error_type, check
use WakeDynamics
use WakeDynamics_Types
use NWTC_Library

implicit none
private
public :: test_wd_maxplanes_warning_suite

contains

!> Collect all exported unit tests
subroutine test_wd_maxplanes_warning_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_warning_issued_only_once", test_warning_issued_only_once) &
               ]
end subroutine

!> Regression test: once the number of tracked wake planes exceeds p%MaxNumPlanes,
!! WD_UpdateStates clamps the count and issues an ErrID_Warn. With a small
!! MaxNumPlanes and no removal (merging/buffer-exit) happening, this condition
!! recurs on every subsequent call. OtherState%MaxPlanesWarned should ensure the
!! warning is only issued on the first occurrence, not on every following call.
subroutine test_warning_issued_only_once(error)
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
   integer(IntKi), parameter    :: MaxNumPlanes = 3
   integer(IntKi)               :: n
   integer(IntKi)               :: warnCount

   errStat = ErrID_None
   errMsg = ''

   ! --- Minimal, valid WakeDynamics (Polar model) initialization input
   InitInp%TurbNum = 1
   InitInp%OutFileRoot = 'test_wd_maxplanes_warning'
   InitInp%MaxNumPlanes = MaxNumPlanes
   ! Generous low-res domain bounds: planes are never dropped via OOB merging,
   ! so the only thing keeping NumPlanes at the cap is the clamp itself.
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

   ! --- Advance well past MaxNumPlanes so the exceeded-cap clamp is hit repeatedly.
   warnCount = 0
   do n = 0, 4*MaxNumPlanes
      errStat = ErrID_None
      errMsg = ''
      call WD_UpdateStates(real(n, DbKi)*DT_low, n, u, p, x, xd, z, OtherState, m, errStat, errMsg)
      call check(error, errStat < AbortErrLev, "WD_UpdateStates failed: "//trim(errMsg))
      if (allocated(error)) return
      if (index(errMsg, 'exceeded the allowed number') > 0) warnCount = warnCount + 1
      call check(error, NINT(xd%NumPlanes) <= p%MaxNumPlanes, "NumPlanes exceeded MaxNumPlanes despite clamp")
      if (allocated(error)) return
   end do

   call check(error, warnCount == 1, "Expected the MaxNumPlanes-exceeded warning exactly once")
   if (allocated(error)) return

end subroutine

end module
