module test_steady_wind

use testdrive, only: new_unittest, unittest_type, error_type, check
use ifw_test_tools
use InflowWind_Subs
use InflowWind_Types

implicit none
private
public :: test_steady_wind_suite

contains

!> Collect all exported unit tests
subroutine test_steady_wind_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_steady_wind_input_single_height", test_steady_wind_input_single_height), &
               new_unittest("test_steady_wind_input_mult_heights", test_steady_wind_input_mult_heights) &
               ]
end subroutine

subroutine test_steady_wind_input_single_height(error)
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
   call check(error, InputFileData%WindType, 1); if (allocated(error)) return
   call check(error, InputFileData%NWindVel, 1); if (allocated(error)) return
   call check(error, InputFileData%WindVziList(1), 90.0_ReKi); if (allocated(error)) return

end subroutine

subroutine test_steady_wind_input_mult_heights(error)
   type(error_type), allocatable, intent(out) :: error
   type(FileInfoType)              :: InFileInfo
   type(InflowWind_InputFile)      :: InputFileData
   character(1024)                 :: PriPath
   integer(IntKi)                  :: TmpErrStat
   character(ErrMsgLen)            :: TmpErrMsg

   PriPath = ""

   InFileInfo = getInputFileData()
   InFileInfo%Lines(9:12) = [ &
                            '          2   NWindVel       - Number of points to output the wind velocity    (0 to 9)                                                                                            ', &
                            '        0,0   WindVxiList    - List of coordinates in the inertial X direction (m)                                                                                                 ', &
                            '        0,0   WindVyiList    - List of coordinates in the inertial Y direction (m)                                                                                                 ', &
                            '     80,100   WindVziList    - List of coordinates in the inertial Z direction (m)                                                                                                 ' &
                            ]

   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

   call check(error, TmpErrStat, ErrID_None, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return
   call check(error, InputFileData%WindType, 1); if (allocated(error)) return
   call check(error, InputFileData%NWindVel, 2); if (allocated(error)) return
   call check(error, InputFileData%WindVziList(1), 80.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%WindVziList(2), 100.0_ReKi); if (allocated(error)) return

end subroutine

end module
