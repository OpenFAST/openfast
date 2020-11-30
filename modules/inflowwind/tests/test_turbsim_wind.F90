module test_turbsim_wind

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

        CHARACTER(16)                   :: expected

        expected = "Wind/08ms.wnd"
        PriPath = ""

        InFileInfo = getInputFileData()
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, TmpErrStat, TmpErrMsg)

        @assertEqual(TmpErrStat, 0)
        @assertEqual(InputFileData%TSFF_FileName, trim(expected))

    end subroutine

end module
