module test_uniform_wind

    use pFUnit_mod
    use ifw_test_tools
    use InflowWind
    use InflowWind_Subs
    use InflowWind_Types

    implicit none

contains

    @test
    subroutine test_uniform_wind_input()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData
        CHARACTER(1024)                 :: PriPath 
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        CHARACTER(16)                   :: expected

        expected = "Wind/08ms.wnd"
        PriPath = ""

        InFileInfo = getInputFileData()
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')
        @assertEqual(trim(expected), InputFileData%Uniform_FileName)
        @assertEqual(90, InputFileData%Uniform_RefHt)
        @assertEqual(125.88, InputFileData%Uniform_RefLength)

    end subroutine

    @test
    subroutine test_uniform_wind_direct_data()

            ! Types for setting up module
        TYPE(InflowWind_InitInputType)                     :: InitInp           !< Input data for initialization
        TYPE(InflowWind_InputType)                         :: InputGuess        !< An initial guess for the input; the input mesh must be defined
        TYPE(InflowWind_ParameterType)                     :: p                 !< Parameters
        TYPE(InflowWind_ContinuousStateType)               :: ContStates        !< Initial continuous states
        TYPE(InflowWind_DiscreteStateType)                 :: DiscStates        !< Initial discrete states
        TYPE(InflowWind_ConstraintStateType)               :: ConstrStateGuess  !< Initial guess of the constraint states
        TYPE(InflowWind_OtherStateType)                    :: OtherStates       !< Initial other/optimization states
        TYPE(InflowWind_OutputType)                        :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
        TYPE(InflowWind_MiscVarType)                       :: m                 !< Misc variables for optimization (not copied in glue code)
        REAL(DbKi)                                         :: TimeInterval      !< Coupling time interval in seconds: InflowWind does not change this.
        TYPE(InflowWind_InitOutputType)                    :: InitOutData 

            ! Variables for testing
        INTEGER                         :: ErrStat
        CHARACTER(ErrMsgLen)            :: ErrMsg
        TYPE(FileInfoType)              :: InFileInfo
        TYPE(FileInfoType)              :: WindType2Data
        CHARACTER(1024), DIMENSION(6)   :: data = (/ &
            '! Wind file for sheared 18 m/s wind with 30 degree direction.    ', &
            '! Time Wind Wind  Vert. Horiz. Vert. LinV Gust                   ', &
            '!      Speed Dir Speed Shear Shear Shear Speed                   ', &
            ' 0.0   12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          ', &
            ' 0.1   12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          ', &
            ' 999.9 12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          ' &
        /)

            ! Error handling
        INTEGER(IntKi)                                        :: TmpErrStat
        CHARACTER(ErrMsgLen)                                  :: TmpErrMsg         !< temporary error message

        InFileInfo = getInputFileDataWindType2()
        CALL InitFileInfo(data, WindType2Data, ErrStat, ErrMsg)

      ! For diagnostic purposes, the following can be used to display the contents
      ! of the InFileInfo data structure.
      ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.

            ! Variable definitions
        InitInp%InputFileName = ""
        InitInp%NumWindPoints = 5
        InitInp%UseInputFile = .FALSE.
        InitInp%RootName = ""
        InitInp%PassedFileData = InFileInfo
        InitInp%WindType2UseInputFile = .FALSE.
        InitInp%WindType2Data = WindType2Data

        CALL InflowWind_Init( InitInp, InputGuess, p, ContStates, DiscStates, &
                        ConstrStateGuess, OtherStates, y, m, TimeInterval, &
                        InitOutData, TmpErrStat, TmpErrMsg)

            ! Results
        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')
        @assertEqual(0.0, p%FlowField%Uniform%Time(1))
        @assertEqual(0.1, p%FlowField%Uniform%Time(2))
        @assertEqual(999.9, p%FlowField%Uniform%Time(3))

        @assertEqual(12.0, p%FlowField%Uniform%VelH(1))
        @assertEqual(12.0, p%FlowField%Uniform%VelH(2))
        @assertEqual(12.0, p%FlowField%Uniform%VelH(3))

    end subroutine

end module
