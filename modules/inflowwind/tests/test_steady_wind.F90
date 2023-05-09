module test_steady_wind

    use pFUnit_mod
    use ifw_test_tools
    use InflowWind_Subs
    use InflowWind_Types

    implicit none

contains

    @test
    subroutine test_steady_wind_input_single_height()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData
        CHARACTER(1024)                 :: PriPath 
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        PriPath = ""

        InFileInfo = getInputFileData()
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')
        @assertEqual(1, InputFileData%WindType)
        @assertEqual(1, InputFileData%NWindVel)
        @assertEqual(90, InputFileData%WindVziList(1))

    end subroutine


    @test
    subroutine test_steady_wind_input_mult_heights()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData 
        CHARACTER(1024)                 :: PriPath
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        PriPath = ""

        InFileInfo = getInputFileData()
        InFileInfo%Lines(9:12) = (/ &
            '          2   NWindVel       - Number of points to output the wind velocity    (0 to 9)                                                                                            ', &
            '        0,0   WindVxiList    - List of coordinates in the inertial X direction (m)                                                                                                 ', &
            '        0,0   WindVyiList    - List of coordinates in the inertial Y direction (m)                                                                                                 ', &
            '     80,100   WindVziList    - List of coordinates in the inertial Z direction (m)                                                                                                 ' &
        /)

        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')
        @assertEqual(1, InputFileData%WindType)
        @assertEqual(2, InputFileData%NWindVel)
        @assertEqual(80, InputFileData%WindVziList(1))
        @assertEqual(100, InputFileData%WindVziList(2))

    end subroutine

end module
