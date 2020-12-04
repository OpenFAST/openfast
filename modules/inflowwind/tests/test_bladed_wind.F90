module test_bladed_wind

    use pFUnit_mod
    use ifw_test_tools
    use InflowWind_Subs
    use InflowWind_Types

    implicit none

contains

    @test
    subroutine test_bladed_wind_input()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData
        CHARACTER(1024)                 :: PriPath 
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        CHARACTER(16)                   :: expected

        expected = "unused.wnd"
        PriPath = ""

        InFileInfo = getInputFileData()
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, TmpErrStat, TmpErrMsg)

        @assertEqual(TmpErrStat, 0)
        @assertEqual(InputFileData%BladedFF_FileName, trim(expected))
        @assertEqual(InputFileData%BladedFF_TowerFile, .FALSE.)

    end subroutine

end module
