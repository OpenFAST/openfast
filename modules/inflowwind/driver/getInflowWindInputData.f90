    use NWTC_IO
    use InflowWind_Subs
    use InflowWind_Types
	USE iso_c_binding

    implicit none
	
	subroutine getInflowWindInputData(data_in, InputFileData, TmpErrStat, TmpErrMsg)
	
    	CHARACTER(KIND=C_CHAR),DIMENSION(65536), INTENT(IN   )  :: data_in
        TYPE(InflowWind_InputFile)             , INTENT(  OUT)  :: InputFileData
        INTEGER(IntKi)                         , INTENT(  OUT)  :: TmpErrStat
        CHARACTER(1024)                        , INTENT(  OUT)  :: TmpErrMsg

        ! Local variables
		CHARACTER(1024), DIMENSION(55)                          :: data
	    CHARACTER(1024)                                         :: tmp_str
		CHARACTER                                               :: tmp_let
		TYPE(FileInfoType)                                      :: InFileInfo
		INTEGER                                                 :: count
		INTEGER                                                 :: idx
		
		TmpErrStat                 = 0
        TmpErrMsg                  = "Completed Successfully!"

        ! Convert the input data into fortran and parse it using '\n'
		count = 0
		idx = 1
		DO
		    IF (COUNT=54) EXIT
			! Convert to Fortran type
			tmp_let = CHARACTER(data_in(idx))
		    IF ()
		END DO

        ! Call the functions
		CALL InitFileInfo(data, InFileInfo, TmpErrStat, TmpErrMsg) ! (in, out, out, out)
		IF (TmpErrStat /= 0) THEN
		    TmpErrMsg = "Failed in InitFileInfo"
		ENDIF
		
        CALL InflowWind_ParseInputFileInfo(InputFileData , InFileInfo, "", TmpErrStat, TmpErrMsg) ! (inout, in, in, out, out)
		IF (TmpErrStat /= 0) THEN
		    TmpErrMsg = "Failed in InflowWind_ParseInputFileInfo"
		ENDIF
 
    end subroutine
	