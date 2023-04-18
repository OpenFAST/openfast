module test_outputs

    use pFUnit_mod
    use ifw_test_tools
    use InflowWind_Subs
    use InflowWind_Types

    implicit none

contains

    @test
    subroutine test_outputs_parsing()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData
        CHARACTER(1024)                 :: PriPath 
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        PriPath = ""

        InFileInfo = getInputFileData()
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')
        @assertEqual(.FALSE.,     InputFileData%SumPrint)
        @assertEqual("Wind1VelX", InputFileData%OutList(1))
        @assertEqual("Wind1VelY", InputFileData%OutList(2))
        @assertEqual("Wind1VelZ", InputFileData%OutList(3))

    end subroutine


    @test
    subroutine test_outputs_parsing_alternate()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData 
        CHARACTER(1024)                 :: PriPath
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        PriPath = ""

        InFileInfo = getInputFileData()
        InFileInfo%Lines(65:67) = (/ &
            'True          SumPrint     - Print summary data to <RootName>.IfW.sum (flag)                                                                                   ', &
            '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)', &
            '"Wind1VelX,Wind1VelY"      - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.                           ' &
        /)

        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')
        @assertEqual(.TRUE.,      InputFileData%SumPrint)
        @assertEqual("Wind1VelX", InputFileData%OutList(1))
        @assertEqual("Wind1VelY", InputFileData%OutList(2))

    end subroutine

end module
