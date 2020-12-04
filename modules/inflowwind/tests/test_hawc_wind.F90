module test_hawc_wind

    use pFUnit_mod
    use ifw_test_tools
    use InflowWind_Subs
    use InflowWind_Types

    implicit none

contains

    @test
    subroutine test_hawc_wind_input()

        TYPE(FileInfoType)              :: InFileInfo
        TYPE(InflowWind_InputFile)      :: InputFileData
        CHARACTER(1024)                 :: PriPath 
        INTEGER(IntKi)                  :: TmpErrStat
        CHARACTER(ErrMsgLen)            :: TmpErrMsg

        CHARACTER(32)                   :: expected_fnu
        CHARACTER(32)                   :: expected_fnv
        CHARACTER(32)                   :: expected_fnw

        PriPath = ""
        expected_fnu = "wasp\Output\basic_5u.bin"
        expected_fnv = "wasp\Output\basic_5v.bin"
        expected_fnw = "wasp\Output\basic_5w.bin"

        InFileInfo = getInputFileData()
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, TmpErrStat, TmpErrMsg)

        @assertEqual(TmpErrStat, 0)

        @assertEqual(InputFileData%HAWC_FileName_u, trim(expected_fnu))
        @assertEqual(InputFileData%HAWC_FileName_v, trim(expected_fnv))
        @assertEqual(InputFileData%HAWC_FileName_w, trim(expected_fnw))
        @assertEqual(InputFileData%HAWC_nx, 64)
        @assertEqual(InputFileData%HAWC_ny, 32)
        @assertEqual(InputFileData%HAWC_nz, 32)
        @assertEqual(InputFileData%HAWC_dx, 16)
        @assertEqual(InputFileData%HAWC_dy, 3)
        @assertEqual(InputFileData%HAWC_dz, 3)
        @assertEqual(InputFileData%HAWC_RefHt, 90)

        @assertEqual(InputFileData%HAWC_ScaleMethod, 1)
        @assertEqual(InputFileData%HAWC_SFx, 1)
        @assertEqual(InputFileData%HAWC_SFy, 1)
        @assertEqual(InputFileData%HAWC_SFz, 1)
        @assertEqual(InputFileData%HAWC_SigmaFx, 12)
        @assertEqual(InputFileData%HAWC_SigmaFy, 8)
        @assertEqual(InputFileData%HAWC_SigmaFz, 2)

        @assertEqual(InputFileData%HAWC_URef, 5)
        @assertEqual(InputFileData%HAWC_ProfileType, 2)
        @assertEqual(InputFileData%HAWC_PLExp, 0)
        @assertEqual(InputFileData%HAWC_Z0, 0.03)
        @assertEqual(InputFileData%HAWC_InitPosition(1), 0)
        @assertEqual(InputFileData%HAWC_InitPosition(2), 0)
        @assertEqual(InputFileData%HAWC_InitPosition(3), 0)

    end subroutine

end module
