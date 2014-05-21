   SUBROUTINE ReadBladeFile(BldFile,BladeInputFileData,UnEc,ErrStat,ErrMsg)

   ! Passed variables:
   TYPE(BladeInputData), INTENT(  OUT):: BladeInputFileData
   CHARACTER(*),         INTENT(IN   ):: BldFile
   INTEGER(IntKi),       INTENT(IN   ):: UnEc
   INTEGER(IntKi),       INTENT(  OUT):: ErrStat                             ! Error status
   CHARACTER(*),         INTENT(  OUT):: ErrMsg                              ! Error message
   
   ! Local variables:
   INTEGER(IntKi)             :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)             :: ErrStat2                                        ! Temporary Error status
   CHARACTER(LEN(ErrMsg))     :: ErrMsg2                                         ! Temporary Err msg   
   REAL(ReKi)                 :: temp_xm2
   REAL(ReKi)                 :: temp_xm3
   REAL(ReKi)                 :: temp_xc2
   REAL(ReKi)                 :: temp_xc3
   REAL(ReKi)                 :: temp_xs2
   REAL(ReKi)                 :: temp_xs3
   INTEGER(IntKi)             :: i
   INTEGER(IntKi)             :: j

   CALL GetNewUnit(UnIn,ErrStat,ErrMsg)

   CALL OpenFInpFile (UnIn,BldFile,ErrStat2,ErrMsg2)

   !  -------------- HEADER -------------------------------------------------------
   ! Skip the header.
   CALL ReadCom(UnIn,BldFile,'unused blade file header line 1',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'unused blade file header line 2',ErrStat2,ErrMsg2,UnEc)

   !  -------------- BLADE PARAMETER-----------------------------------------------
   CALL ReadCom(UnIn,BldFile,'blade parameters',ErrStat2,ErrMsg2,UnEc)

   CALL ReadVar(UnIn,BldFile,BladeInputFileData%station_total,'station_total','Number of blade input stations',ErrStat2,ErrMsg2,UnEc)

   CALL AllocAry(BladeInputFileData%stiff0,6,6,BladeInputFileData%station_total,'Cross-sectional 6 by 6 stiffness matrix',ErrStat2,ErrMsg2)
   BladeInputFileData%stiff0(:,:,:) = 0.0D0 
   CALL AllocAry(BladeInputFileData%mass0,6,6,BladeInputFileData%station_total,'Cross-sectional 6 by 6 mass matrix',ErrStat2,ErrMsg2)
   BladeInputFileData%mass0(:,:,:) = 0.0D0 
   CALL AllocAry(BladeInputFileData%station_eta,BladeInputFileData%station_total,'Station eta array',ErrStat2,ErrMsg2)
   BladeInputFileData%station_eta(:) = 0.0D0 

   CALL ReadVar(UnIn,BldFile,BladeInputFileData%format_index,'format_index','Index of input format',ErrStat2,ErrMsg2,UnEc)

!  -------------- DISTRIBUTED PROPERTIES--------------------------------------------
  CALL ReadCom(UnIn,BldFile,'Distributed properties',ErrStat2,ErrMsg2,UnEc)

   IF(BladeInputFileData%format_index .EQ. 1) THEN
       DO i=1,BladeInputFileData%station_total
           READ(UnIn,*) BladeInputFileData%station_eta(i)
           DO j=1,6
               CALL ReadAry(UnIn,BldFile,BladeInputFileData%stiff0(j,:,i),6,'siffness_matrix',&
                       'Blade C/S stiffness matrix',ErrStat2,ErrMsg2,UnEc)
           ENDDO
           DO j=1,6
               CALL ReadAry(UnIn,BldFile,BladeInputFileData%mass0(j,:,i),6,'siffness_matrix',&
                       'Blade C/S stiffness matrix',ErrStat2,ErrMsg2,UnEc)
           ENDDO
       ENDDO
   ELSEIF(BladeInputFileData%format_index .EQ. 2) THEN
       DO i=1,BladeInputFileData%station_total
           READ(UnIn,*) BladeInputFileData%station_eta(i)
           READ(UnIn,*) BladeInputFileData%stiff0(1,1,i)
           READ(UnIn,*) BladeInputFileData%stiff0(5,5,i),BladeInputFileData%stiff0(6,6,i),&
                        BladeInputFileData%stiff0(5,6,i)
           READ(UnIn,*) BladeInputFileData%stiff0(4,4,i)
           READ(UnIn,*) BladeInputFileData%stiff0(2,2,i),BladeInputFileData%stiff0(3,3,i),&
                        BladeInputFileData%stiff0(2,3,i)
           BladeInputFileData%stiff0(6,5,i) = BladeInputFileData%stiff0(5,6,i)
           BladeInputFileData%stiff0(3,2,i) = BladeInputFileData%stiff0(2,3,i)

           READ(UnIn,*) BladeInputFileData%mass0(1,1,i)
           READ(UnIn,*) BladeInputFileData%mass0(4,4,i),BladeInputFileData%mass0(5,5,i),&
                        BladeInputFileData%mass0(6,6,i)
           temp_xm2 = 0.0D0
           temp_xm3 = 0.0D0
           READ(UnIn,*) temp_xm2,temp_xm3
           BladeInputFileData%mass0(2,2,i) = BladeInputFileData%mass0(1,1,i)
           BladeInputFileData%mass0(3,3,i) = BladeInputFileData%mass0(1,1,i)
           BladeInputFileData%mass0(1,5,i) = BladeInputFileData%mass0(1,1,i)*temp_xm3
           BladeInputFileData%mass0(1,6,i) =-BladeInputFileData%mass0(1,1,i)*temp_xm2
           BladeInputFileData%mass0(2,4,i) =-BladeInputFileData%mass0(1,1,i)*temp_xm3
           BladeInputFileData%mass0(3,4,i) = BladeInputFileData%mass0(1,1,i)*temp_xm2
           BladeInputFileData%mass0(5,1,i) = BladeInputFileData%mass0(1,5,i)
           BladeInputFileData%mass0(6,1,i) = BladeInputFileData%mass0(1,6,i)
           BladeInputFileData%mass0(4,2,i) = BladeInputFileData%mass0(2,4,i)
           BladeInputFileData%mass0(4,3,i) = BladeInputFileData%mass0(3,4,i)
       ENDDO
   ENDIF

   END SUBROUTINE ReadBladeFile
