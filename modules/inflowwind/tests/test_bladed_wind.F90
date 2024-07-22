module test_bladed_wind

use testdrive, only: new_unittest, unittest_type, error_type, check
use ifw_test_tools
use InflowWind_Subs
use InflowWind_Types

implicit none

private
public :: test_bladed_wind_suite

contains

!> Collect all exported unit tests
subroutine test_bladed_wind_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [new_unittest("test_bladed_wind_input", test_bladed_wind_input)]
end subroutine

subroutine test_bladed_wind_input(error)
   type(error_type), allocatable, intent(out) :: error

   type(FileInfoType)              :: InFileInfo
   type(InflowWind_InputFile)      :: InputFileData
   character(1024)                 :: PriPath
   integer(IntKi)                  :: TmpErrStat
   character(ErrMsgLen)            :: TmpErrMsg

   character(16)                   :: expected

   expected = "unused"
   PriPath = ""

   InFileInfo = getInputFileData()
   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

   call check(error, TmpErrStat, 0, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return
   call check(error, InputFileData%BladedFF_FileName, trim(expected)); if (allocated(error)) return
   call check(error, InputFileData%BladedFF_TowerFile, .false.); if (allocated(error)) return

end subroutine

end module
