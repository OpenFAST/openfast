module test_outputs

use testdrive, only: new_unittest, unittest_type, error_type, check
use ifw_test_tools
use InflowWind_Subs
use InflowWind_Types

implicit none
private
public :: test_outputs_suite

contains

!> Collect all exported unit tests
subroutine test_outputs_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_outputs_parsing", test_outputs_parsing), &
               new_unittest("test_outputs_parsing_alternate", test_outputs_parsing_alternate) &
               ]
end subroutine

subroutine test_outputs_parsing(error)
   type(error_type), allocatable, intent(out) :: error
   type(FileInfoType)              :: InFileInfo
   type(InflowWind_InputFile)      :: InputFileData
   character(1024)                 :: PriPath
   integer(IntKi)                  :: TmpErrStat
   character(ErrMsgLen)            :: TmpErrMsg

   PriPath = ""

   InFileInfo = getInputFileData()
   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

   call check(error, TmpErrStat, ErrID_None, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return
   call check(error, InputFileData%SumPrint, .false.); if (allocated(error)) return
   call check(error, InputFileData%OutList(1), "Wind1VelX"); if (allocated(error)) return
   call check(error, InputFileData%OutList(2), "Wind1VelY"); if (allocated(error)) return
   call check(error, InputFileData%OutList(3), "Wind1VelZ"); if (allocated(error)) return
end subroutine

subroutine test_outputs_parsing_alternate(error)
   type(error_type), allocatable, intent(out) :: error
   type(FileInfoType)              :: InFileInfo
   type(InflowWind_InputFile)      :: InputFileData
   character(1024)                 :: PriPath
   integer(IntKi)                  :: TmpErrStat
   character(ErrMsgLen)            :: TmpErrMsg

   PriPath = ""

   InFileInfo = getInputFileData()
   InFileInfo%Lines(65:67) = [ &
                             'True          SumPrint     - Print summary data to <RootName>.IfW.sum (flag)                                                                                   ', &
                             '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)', &
                             '"Wind1VelX,Wind1VelY"      - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.                           ' &
                             ]

   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

   call check(error, TmpErrStat, ErrID_None, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return
   call check(error, InputFileData%SumPrint, .true.); if (allocated(error)) return
   call check(error, InputFileData%OutList(1), "Wind1VelX"); if (allocated(error)) return
   call check(error, InputFileData%OutList(2), "Wind1VelY"); if (allocated(error)) return
end subroutine

end module
