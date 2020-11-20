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
        @assertEqual(InputFileData%FF%RefHt, 90)

        @assertEqual(InputFileData%FF%ScaleMethod, 1)
        @assertEqual(InputFileData%FF%SF(1), 1)
        @assertEqual(InputFileData%FF%SF(2), 1)
        @assertEqual(InputFileData%FF%SF(3), 1)
        @assertEqual(InputFileData%FF%SigmaF(1), 12)
        @assertEqual(InputFileData%FF%SigmaF(2), 8)
        @assertEqual(InputFileData%FF%SigmaF(3), 2)

        @assertEqual(InputFileData%FF%URef, 5)
        @assertEqual(InputFileData%FF%WindProfileType, 2)
        @assertEqual(InputFileData%FF%PLExp, 0)
        @assertEqual(InputFileData%FF%Z0, 0.03)
        @assertEqual(InputFileData%FF%XOffset, 0)

    end subroutine

end module
