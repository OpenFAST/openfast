module test_hawc_wind

use testdrive, only: new_unittest, unittest_type, error_type, check
use ifw_test_tools
use InflowWind_Subs
use InflowWind_Types

implicit none
private
public :: test_hawc_wind_suite

contains

!> Collect all exported unit tests
subroutine test_hawc_wind_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [new_unittest("test_hawc_wind_input", test_hawc_wind_input)]
end subroutine

subroutine test_hawc_wind_input(error)
   type(error_type), allocatable, intent(out) :: error

   type(FileInfoType)              :: InFileInfo
   type(InflowWind_InputFile)      :: InputFileData
   character(1024)                 :: PriPath
   integer(IntKi)                  :: TmpErrStat
   character(ErrMsgLen)            :: TmpErrMsg

   character(32)                   :: expected_fnu
   character(32)                   :: expected_fnv
   character(32)                   :: expected_fnw

   PriPath = ""
   expected_fnu = "wasp/Output/basic_5u.bin"
   expected_fnv = "wasp/Output/basic_5v.bin"
   expected_fnw = "wasp/Output/basic_5w.bin"

   InFileInfo = getInputFileData()
   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

   call check(error, TmpErrStat, ErrID_None, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return

   call check(error, InputFileData%HAWC_FileName_u, trim(expected_fnu)); if (allocated(error)) return
   call check(error, InputFileData%HAWC_FileName_v, trim(expected_fnv)); if (allocated(error)) return
   call check(error, InputFileData%HAWC_FileName_w, trim(expected_fnw)); if (allocated(error)) return
   call check(error, InputFileData%HAWC_nx, 64); if (allocated(error)) return
   call check(error, InputFileData%HAWC_ny, 32); if (allocated(error)) return
   call check(error, InputFileData%HAWC_nz, 32); if (allocated(error)) return
   call check(error, InputFileData%HAWC_dx, 16.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%HAWC_dy, 3.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%HAWC_dz, 3.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%RefHt, 90.0_ReKi); if (allocated(error)) return

   call check(error, InputFileData%FF%ScaleMethod, 1); if (allocated(error)) return
   call check(error, InputFileData%FF%SF(1), 1.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%SF(2), 1.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%SF(3), 1.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%SigmaF(1), 12.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%SigmaF(2), 8.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%SigmaF(3), 2.0_ReKi); if (allocated(error)) return

   call check(error, InputFileData%FF%URef, 5.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%WindProfileType, 2); if (allocated(error)) return
   call check(error, InputFileData%FF%PLExp, 0.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%Z0, 0.03_ReKi); if (allocated(error)) return
   call check(error, InputFileData%FF%XOffset, 0.0_ReKi); if (allocated(error)) return

end subroutine

end module
