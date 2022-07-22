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
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, PriPath, "inputFile.inp", "test.ech", .false., -1, TmpErrStat, TmpErrMsg)

        @assertEqual(0, TmpErrStat, message='Error message: '//trim(TmpErrMsg)//NewLine//'ErrStat: ')

        @assertEqual(trim(expected_fnu), InputFileData%HAWC_FileName_u)
        @assertEqual(trim(expected_fnv), InputFileData%HAWC_FileName_v)
        @assertEqual(trim(expected_fnw), InputFileData%HAWC_FileName_w)
        @assertEqual(64, InputFileData%HAWC_nx)
        @assertEqual(32, InputFileData%HAWC_ny)
        @assertEqual(32, InputFileData%HAWC_nz)
        @assertEqual(16, InputFileData%HAWC_dx)
        @assertEqual(3,  InputFileData%HAWC_dy)
        @assertEqual(3,  InputFileData%HAWC_dz)
        @assertEqual(90, InputFileData%FF%RefHt)

        @assertEqual(1,  InputFileData%FF%ScaleMethod)
        @assertEqual(1,  InputFileData%FF%SF(1))
        @assertEqual(1,  InputFileData%FF%SF(2))
        @assertEqual(1,  InputFileData%FF%SF(3))
        @assertEqual(12, InputFileData%FF%SigmaF(1))
        @assertEqual(8,  InputFileData%FF%SigmaF(2))
        @assertEqual(2,  InputFileData%FF%SigmaF(3))

        @assertEqual(5,  InputFileData%FF%URef)
        @assertEqual(2,  InputFileData%FF%WindProfileType)
        @assertEqual(0,  InputFileData%FF%PLExp)
        @assertEqual(0.03, InputFileData%FF%Z0)
        @assertEqual(0,  InputFileData%FF%XOffset)

    end subroutine

end module
