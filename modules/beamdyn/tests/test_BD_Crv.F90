module test_BD_Crv

use test_tools
use BeamDyn_Subs

implicit none

private
public :: test_BD_Crv_suite

contains

!> Collect all exported unit tests
subroutine test_BD_Crv_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_BD_CheckRotMat", test_BD_CheckRotMat), &
               new_unittest("test_BD_ComputeIniNodalCrv", test_BD_ComputeIniNodalCrv), &
               new_unittest("test_BD_CrvCompose", test_BD_CrvCompose), &
               new_unittest("test_BD_CrvExtractCrv", test_BD_CrvExtractCrv), &
               new_unittest("test_BD_CrvMatrixH", test_BD_CrvMatrixH), &
               new_unittest("test_BD_CrvMatrixR", test_BD_CrvMatrixR) &
               ]
end subroutine

subroutine test_BD_CheckRotMat(error)
   type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - known valid rotation matrix: pi about x-axis
   ! - known invalid rotation matrix: halve the angle of the diagonal elements

   real(BDKi)           :: n(3)
   real(BDKi)           :: angle
   real(BDKi)           :: testR(3, 3)
   integer(IntKi)       :: ErrStat
   character(ErrMsgLen) :: ErrMsg
   character(1024)      :: testname

   ! set the rotation axis and angle for all tests
   n = [1, 0, 0] ! x axis
   angle = Pi

   ! --------------------------------------------------------------------------
   testname = "known valid rotation matrix: pi about x-axis:"
   testR = calcRotationMatrix(angle, n)
   call BD_CheckRotMat(testR, ErrStat, ErrMsg)
   call check(error, 0, ErrStat, testname); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "known invalid rotation matrix: halve the angle of the diagonal elements:"
   ! this should produce a fatal error (ErrStat = 4)
   testR(:, 2) = [testR(1, 2), cos(real(Pi / 2, BDKi)), testR(3, 2)]
   testR(:, 3) = [testR(1, 2), testR(2, 2), cos(real(Pi / 2, BDKi))]
   call BD_CheckRotMat(testR, ErrStat, ErrMsg)
   call check(error, 4, ErrStat, testname); if (allocated(error)) return

end subroutine

subroutine test_BD_ComputeIniNodalCrv(error)
   type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - simple rotation with known parameters: Pi on xaxis
   ! - 0 rotation
   ! - small rotation with baseline WM parameters calculated

   real(BDKi), dimension(3, 3) :: r
   real(BDKi), dimension(3)   :: test_wmparams, baseline_wmparams
   real(BDKi)                 :: angle, param, n(3)
   character(1024)            :: testname
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg

   ! --------------------------------------------------------------------------
   testname = "Tangent aligned with z-axis and 0 degree twist:"
   n = [real(0.0, BDKi), real(0.0, BDKi), real(1.0, BDKi)] ! tangent axis
   angle = 0

   ! Baseline Wiener-Milenkovic parameters
   baseline_wmparams = [0.0, 0.0, 0.0]

   call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "Tangent at 45 degree w.r.t. y-axis and 0 degree twist:"
   n = [1.0_BDKi / sqrt(2.0_BDKi), 0.0_BDKi, 1.0_BDKi / sqrt(2.0_BDKi)] ! tangent axis
   angle = 0.0_BDKi

   ! Baseline Wiener-Milenkovic parameters
   param = 4 * tan((Pi_D / 4) / 4)
   baseline_wmparams = [real(0.0, BDKi), param, real(0.0, BDKi)]

   call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "Tangent at -45 degree w.r.t. x-axis and 0 degree twist:"
   n = [0.0_BDKi, 1.0_BDKi / sqrt(2.0_BDKi), 1.0_BDKi / sqrt(2.0_BDKi)] ! tangent axis
   angle = 0.0_BDKi

   ! Baseline Wiener-Milenkovic parameters
   param = 4.*tan((-Pi_D / 4.) / 4.)
   baseline_wmparams = [param, real(0.0, BDKi), real(0.0, BDKi)]

   call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "Tangent along z-axis with 45 degree twist:"
   n = [real(0.0, BDKi), real(0.0, BDKi), 1.0_BDKi] ! tangent axis
   angle = 45.0_BDKi

   ! Baseline Wiener-Milenkovic parameters
   param = 4.*tan((Pi_D / 4.) / 4.)
   baseline_wmparams = [real(0.0, BDKi), real(0.0, BDKi), param]

   call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

end subroutine

subroutine test_BD_CrvCompose(error)
   type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - both rotation angles 0, no transpose of input rotations (flag = 0)
   ! - delta2 > 0, no transpose of input rotations (flag = 0)
   ! - delta2 < 0, no transpose of input rotations (flag = 0)
   ! - flag = 1
   ! - flag = 2
   ! - flag = 3

   ! input rotation axis and angle
   real(BDKi), dimension(3)   :: n1, n2
   real(BDKi)                 :: angle1, angle2

   ! result rotations
   real(BDKi), dimension(3, 3)   :: testrotation, baselinerotation, r1, r2
   real(BDKi), dimension(3)      :: composedparams

   ! other test settings
   integer                    :: flag
   character(1024)            :: testname

   ! set the rotation axes for all tests
   n1 = [1, 0, 0] ! x axis
   n2 = [0, 0, 1] ! z axis

   ! --------------------------------------------------------------------------
   testname = "both rotation angles 0, no transpose of input rotations (flag = 0):"
   angle1 = 0 ! 0 degrees
   angle2 = 0 ! 0 degrees
   flag = 0

   ! both rotations should return an identity matrix
   r1 = calcRotationMatrix(angle1, n1)
   r2 = calcRotationMatrix(angle2, n2)
   baselinerotation = matmul(r1, r2)

   call BD_CrvCompose(composedparams, 4 * tan(angle1 / 4) * n1, 4 * tan(angle2 / 4) * n2, flag)
   call BD_CrvMatrixR(composedparams, testrotation)

   call check_array(error, baselinerotation, testrotation, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "delta2 > 0, no transpose of input rotations (flag = 0):"
   angle1 = PiBy2_D ! 90 degrees
   angle2 = PiBy2_D ! 90 degrees
   flag = 0

   r1 = calcRotationMatrix(angle1, n1)
   r2 = calcRotationMatrix(angle2, n2)
   baselinerotation = matmul(r1, r2)

   call BD_CrvCompose(composedparams, 4 * tan(angle1 / 4) * n1, 4 * tan(angle2 / 4) * n2, flag)
   call BD_CrvMatrixR(composedparams, testrotation)

   call check_array(error, baselinerotation, testrotation, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "delta2 < 0, no transpose of input rotations (flag = 0):"
   angle1 = PiBy2_D ! 90 degrees
   angle2 = 1.5 * Pi  ! 270 degrees
   flag = 0

   r1 = calcRotationMatrix(angle1, n1)
   r2 = calcRotationMatrix(angle2, n2)
   baselinerotation = matmul(r1, r2)

   call BD_CrvCompose(composedparams, 4 * tan(angle1 / 4) * n1, 4 * tan(angle2 / 4) * n2, flag)
   call BD_CrvMatrixR(composedparams, testrotation)

   call check_array(error, baselinerotation, testrotation, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "delta2 > 0, transpose of first rotation (flag = 1):"
   angle1 = PiBy2_D ! 90 degrees
   angle2 = PiBy2_D ! 90 degrees
   flag = 1

   r1 = calcRotationMatrix(-angle1, n1)
   r2 = calcRotationMatrix(angle2, n2)
   baselinerotation = matmul(r1, r2)

   call BD_CrvCompose(composedparams, 4 * tan(angle1 / 4) * n1, 4 * tan(angle2 / 4) * n2, flag)
   call BD_CrvMatrixR(composedparams, testrotation)

   call check_array(error, baselinerotation, testrotation, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "delta2 > 0, transpose of second rotation (flag = 2):"
   angle1 = PiBy2_D ! 90 degrees
   angle2 = PiBy2_D ! 90 degrees
   flag = 2

   r1 = calcRotationMatrix(angle1, n1)
   r2 = calcRotationMatrix(-angle2, n2)
   baselinerotation = matmul(r1, r2)

   call BD_CrvCompose(composedparams, 4 * tan(angle1 / 4) * n1, 4 * tan(angle2 / 4) * n2, flag)
   call BD_CrvMatrixR(composedparams, testrotation)

   call check_array(error, baselinerotation, testrotation, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "delta2 > 0, transpose of both rotations (flag = 3):"
   angle1 = PiBy2_D ! 90 degrees
   angle2 = PiBy2_D ! 90 degrees
   flag = 3

   r1 = calcRotationMatrix(-angle1, n1)
   r2 = calcRotationMatrix(-angle2, n2)
   baselinerotation = matmul(r1, r2)

   call BD_CrvCompose(composedparams, 4 * tan(angle1 / 4) * n1, 4 * tan(angle2 / 4) * n2, flag)
   call BD_CrvMatrixR(composedparams, testrotation)

   call check_array(error, baselinerotation, testrotation, testname, tolerance); if (allocated(error)) return

end subroutine

subroutine test_BD_CrvExtractCrv(error)
   type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - simple rotation with known parameters: Pi on xaxis
   ! - 0 rotation
   ! - small rotation with baseline WM parameters calculated

   real(BDKi), dimension(3, 3) :: r
   real(BDKi), dimension(3)   :: test_wmparams, baseline_wmparams
   real(BDKi)                 :: angle, n(3)
   character(1024)            :: testname
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg

   ! set the rotation axis for all tests
   n = [1, 0, 0] ! x axis

   ! --------------------------------------------------------------------------
   testname = "simple rotation with known parameters: Pi on xaxis:"
   angle = Pi_D

   ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
   baseline_wmparams = [4.0, 0.0, 0.0]

   r = RonXAxis(angle)
   call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "0 rotation:"
   angle = 0

   ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
   baseline_wmparams = [0.0, 0.0, 0.0]

   r = RonXAxis(angle)
   call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "small rotation with baseline WM parameters calculated:"
   angle = 0.1 * Pi_D

   ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
   baseline_wmparams = 4 * tan(angle / 4) * n

   r = RonXAxis(angle)
   call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

   call check(error, ErrID_None, ErrStat, testname); if (allocated(error)) return
   call check_array(error, baseline_wmparams, test_wmparams, testname, tolerance); if (allocated(error)) return

end subroutine

subroutine test_BD_CrvMatrixH(error)
   type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - simple rotation with known parameters: Pi on xaxis
   ! - 0 rotation
   ! - small rotation with baseline WM parameters calculated

   ! TODO
   ! invalid wm parameters (if thats a thing)
   ! does the implemented WM formulation have any boundaries?

   real(BDKi), dimension(3, 3) :: testH, baselineH
   real(BDKi), dimension(3)   :: wmparams
   real(BDKi)                 :: angle, n(3)
   character(1024)            :: testname
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg

   ! set the rotation axis for all tests
   n = [1., 0., 0.] ! x axis

   ! --------------------------------------------------------------------------
   testname = "simple rotation with known parameters: Pi on xaxis:"
   angle = Pi_D

   ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
   wmparams = [4.0, 0.0, 0.0]

   baselineH = H(wmparams)

   call BD_CrvMatrixH(wmparams, testH)

   call check_array(error, baselineH, testH, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "0 rotation:"
   angle = 0

   ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
   wmparams = [0.0, 0.0, 0.0]

   baselineH = H(wmparams)

   call BD_CrvMatrixH(wmparams, testH)

   call check_array(error, baselineH, testH, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "small rotation with baseline WM parameters calculated:"
   angle = 0.1 * Pi_D

   ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
   wmparams = 4.*tan(angle / 4.) * n

   baselineH = H(wmparams)

   call BD_CrvMatrixH(wmparams, testH)

   call check_array(error, baselineH, testH, testname, tolerance); if (allocated(error)) return

contains
   function H(c)
      real(BDKi)          :: c0, c(3)
      real(BDKi)          :: H(3, 3)

      c0 = 2.0 - dot_product(c, c) / 8.0

      H(1, :) = [c0 + c(1) * c(1) / 4., c(1) * c(2) / 4.-c(3), c(1) * c(3) / 4.+c(2)]
      H(2, :) = [c(1) * c(2) / 4.+c(3), c0 + c(2) * c(2) / 4., c(2) * c(3) / 4.-c(1)]
      H(3, :) = [c(1) * c(3) / 4.-c(2), c(2) * c(3) / 4.+c(1), c0 + c(3) * c(3) / 4.]
      H = 2.*H / (4.-c0)**2

   end function
end subroutine

subroutine test_BD_CrvMatrixR(error)
   type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - simple rotation with known parameters: Pi on xaxis
   ! - 0 rotation
   ! - small rotation with baseline WM parameters calculated

   ! TODO
   ! invalid wm parameters (if thats a thing)
   ! does the implemented WM formulation have any boundaries?

   real(BDKi), dimension(3, 3) :: testR, baselineR
   real(BDKi), dimension(3)   :: wmparams
   real(BDKi)                 :: angle, n(3)
   character(1024)            :: testname
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg

   ! set the rotation axis for all tests
   n = [1, 0, 0] ! x axis

   ! --------------------------------------------------------------------------
   testname = "simple rotation with known parameters: Pi on xaxis:"
   angle = Pi_D

   baselineR = RonXAxis(angle)

   ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
   wmparams = [4.0, 0.0, 0.0]
   call BD_CrvMatrixR(wmparams, testR)

   call check_array(error, baselineR, testR, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "0 rotation:"
   angle = 0

   baselineR = RonXAxis(angle)

   ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
   wmparams = [0.0, 0.0, 0.0]
   call BD_CrvMatrixR(wmparams, testR)

   call check_array(error, baselineR, testR, testname, tolerance); if (allocated(error)) return

   ! --------------------------------------------------------------------------
   testname = "small rotation with baseline WM parameters calculated:"
   angle = 0.1 * Pi_D

   baselineR = RonXAxis(angle)

   ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
   wmparams = 4 * tan(angle / 4) * n
   call BD_CrvMatrixR(wmparams, testR)

   call check_array(error, baselineR, testR, testname, tolerance); if (allocated(error)) return

end subroutine

! this is actually an integration test not a unit test...
subroutine test_BD_ExtractRelativeRotation(error)
   type(error_type), allocatable, intent(out) :: error
   real(BDKi), dimension(3)   :: rr
   character(1024)            :: testname
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg

   type(BD_ParameterType)     :: parametertype
   type(BD_OtherStateType)    :: otherstate

   testname = "static simple beam under gravity:"
   otherstate = simpleOtherState()
   parametertype = simpleParameterType(1, 16, 16, 0, 0)
   call ExtractRelativeRotation(identity(), parametertype, otherstate, rr, ErrStat, ErrMsg)
   call check_array(error, rr, [0.0_BDKi, 0.0_BDKi, 0.0_BDKi], testname, tolerance); if (allocated(error)) return
end subroutine

end module
