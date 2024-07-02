module test_uniform_wind

use testdrive, only: new_unittest, unittest_type, error_type, check
use ifw_test_tools
use InflowWind
use InflowWind_Subs
use InflowWind_Types

implicit none
private
public :: test_uniform_wind_suite

contains

!> Collect all exported unit tests
subroutine test_uniform_wind_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_uniform_wind_input", test_uniform_wind_input), &
               new_unittest("test_uniform_wind_direct_data", test_uniform_wind_direct_data) &
               ]
end subroutine

subroutine test_uniform_wind_input(error)
   type(error_type), allocatable, intent(out) :: error
   type(FileInfoType)              :: InFileInfo
   type(InflowWind_InputFile)      :: InputFileData
   character(1024)                 :: PriPath
   integer(IntKi)                  :: TmpErrStat
   character(ErrMsgLen)            :: TmpErrMsg

   character(16)                   :: expected

   expected = "Wind/08ms.wnd"
   PriPath = ""

   InFileInfo = getInputFileData()
   call InflowWind_ParseInputFileInfo(InputFileData, InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

   call check(error, TmpErrStat, ErrID_None, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return
   call check(error, InputFileData%Uniform_FileName, trim(expected)); if (allocated(error)) return
   call check(error, InputFileData%Uniform_RefHt, 90.0_ReKi); if (allocated(error)) return
   call check(error, InputFileData%Uniform_RefLength, 125.88_ReKi); if (allocated(error)) return

end subroutine

subroutine test_uniform_wind_direct_data(error)
   type(error_type), allocatable, intent(out) :: error

   ! Types for setting up module
   type(InflowWind_InitInputType)                     :: InitInp           !< Input data for initialization
   type(InflowWind_InputType)                         :: InputGuess        !< An initial guess for the input; the input mesh must be defined
   type(InflowWind_ParameterType)                     :: p                 !< Parameters
   type(InflowWind_ContinuousStateType)               :: ContStates        !< Initial continuous states
   type(InflowWind_DiscreteStateType)                 :: DiscStates        !< Initial discrete states
   type(InflowWind_ConstraintStateType)               :: ConstrStateGuess  !< Initial guess of the constraint states
   type(InflowWind_OtherStateType)                    :: OtherStates       !< Initial other/optimization states
   type(InflowWind_OutputType)                        :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   type(InflowWind_MiscVarType)                       :: m                 !< Misc variables for optimization (not copied in glue code)
   real(DbKi)                                         :: TimeInterval      !< Coupling time interval in seconds: InflowWind does not change this.
   type(InflowWind_InitOutputType)                    :: InitOutData

   ! Variables for testing
   integer                         :: ErrStat
   character(ErrMsgLen)            :: ErrMsg
   type(FileInfoType)              :: InFileInfo
   type(FileInfoType)              :: WindType2Info
   character(1024), dimension(6)   :: data = [ &
                                      '! Wind file for sheared 18 m/s wind with 30 degree direction.    ', &
                                      '! Time Wind Wind  Vert. Horiz. Vert. LinV Gust                   ', &
                                      '!      Speed Dir Speed Shear Shear Shear Speed                   ', &
                                      ' 0.0   12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          ', &
                                      ' 0.1   12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          ', &
                                      ' 999.9 12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          ' &
                                      ]

   ! Error handling
   integer(IntKi)                                        :: TmpErrStat
   character(ErrMsgLen)                                  :: TmpErrMsg         !< temporary error message

   InFileInfo = getInputFileDataWindType2()
   call InitFileInfo(data, WindType2Info, ErrStat, ErrMsg)

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
   ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.

   ! Variable definitions
   InitInp%InputFileName = ""
   InitInp%NumWindPoints = 5
   InitInp%FilePassingMethod = 1_IntKi
   InitInp%RootName = ""
   InitInp%PassedFileInfo = InFileInfo
   InitInp%WindType2UseInputFile = .false.
   InitInp%WindType2Info = WindType2Info

   call InflowWind_Init(InitInp, InputGuess, p, ContStates, DiscStates, &
                        ConstrStateGuess, OtherStates, y, m, TimeInterval, &
                        InitOutData, TmpErrStat, TmpErrMsg)

   ! Results
   call check(error, TmpErrStat, ErrID_None, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: '); if (allocated(error)) return
   call check(error, p%FlowField%Uniform%Time(1), 0.0_ReKi); if (allocated(error)) return
   call check(error, p%FlowField%Uniform%Time(2), 0.1_ReKi); if (allocated(error)) return
   call check(error, p%FlowField%Uniform%Time(3), 999.9_ReKi); if (allocated(error)) return

   call check(error, p%FlowField%Uniform%VelH(1), 12.0_ReKi); if (allocated(error)) return
   call check(error, p%FlowField%Uniform%VelH(2), 12.0_ReKi); if (allocated(error)) return
   call check(error, p%FlowField%Uniform%VelH(3), 12.0_ReKi); if (allocated(error)) return

end subroutine

end module
